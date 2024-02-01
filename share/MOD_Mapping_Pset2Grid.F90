#include <define.h>

MODULE MOD_Mapping_Pset2Grid

!----------------------------------------------------------------------------
! DESCRIPTION:
!
!    Mapping data types and subroutines from vector data defined on pixelsets
!    to gridded data.
!
!    Notice that:
!    1. A mapping can be built with method mapping%build.
!    2. Overloaded method "map" can map 1D, 2D or 3D vector data to gridded data 
!       by using area weighted scheme. 
!    3. Method "map_split" can split data in a vector according to pixelset type
!       and map data to 3D gridded data. 
!       The dimensions are from [vector] to [type,lon,lat].
! 
! Created by Shupeng Zhang, May 2023
!----------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_DataType
   IMPLICIT NONE

   ! ------
   type :: mapping_pset2grid_type

      type(grid_type) :: grid
      integer :: npset

      type(grid_list_type), allocatable :: glist (:)

      type(pointer_int32_2d), allocatable :: address(:)
      type(pointer_real8_1d), allocatable :: olparea(:)

   CONTAINS

      procedure, PUBLIC :: build => mapping_pset2grid_build

      procedure, PRIVATE :: map_2d => map_p2g_2d
      procedure, PRIVATE :: map_3d => map_p2g_3d
      procedure, PRIVATE :: map_4d => map_p2g_4d
      generic, PUBLIC :: map => map_2d, map_3d, map_4d

      procedure, PUBLIC  :: map_split => map_p2g_split_to_3d

      final :: mapping_pset2grid_free_mem

   END type mapping_pset2grid_type

!-----------------------
CONTAINS

   !------------------------------------------
   SUBROUTINE mapping_pset2grid_build (this, pixelset, fgrid, pctpset)

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_Block
   USE MOD_Pixel
   USE MOD_Grid
   USE MOD_Pixelset
   USE MOD_DataType
   USE MOD_Mesh
   USE MOD_Utils
   USE MOD_SPMD_Task
   IMPLICIT NONE

   class (mapping_pset2grid_type) :: this

   type(pixelset_type), intent(in) :: pixelset
   type(grid_type),     intent(in) :: fgrid
   real(r8), optional,  intent(in) :: pctpset (:)

   ! Local variables
   type(pointer_real8_1d), allocatable :: afrac(:)
   type(grid_list_type),   allocatable :: gfrom(:)
   type(pointer_int32_1d), allocatable :: list_lat(:)
   integer,  allocatable :: ng_lat(:)
   integer,  allocatable :: ys(:), yn(:), xw(:), xe(:)
   integer,  allocatable :: xlist(:), ylist(:)
   integer,  allocatable :: ipt(:)
   logical,  allocatable :: msk(:)

   integer  :: ie, iset
   integer  :: ng, ig, ng_all, iloc
   integer  :: npxl, ipxl, ilat, ilon
   integer  :: iworker, iproc, idest, isrc, nrecv
   integer  :: rmesg(2), smesg(2)
   integer  :: iy, ix, xblk, yblk, xloc, yloc
   real(r8) :: lat_s, lat_n, lon_w, lon_e, area
   logical  :: is_new

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_master) THEN
         write(*,"('Making mapping from pixel set to grid: ', I7, A, I7, A)") &
            fgrid%nlat, ' grids in latitude', fgrid%nlon, ' grids in longitude'
      ENDIF

      IF (allocated(this%grid%xblk)) deallocate(this%grid%xblk)
      IF (allocated(this%grid%yblk)) deallocate(this%grid%yblk)
      IF (allocated(this%grid%xloc)) deallocate(this%grid%xloc)
      IF (allocated(this%grid%yloc)) deallocate(this%grid%yloc)
      allocate (this%grid%xblk (size(fgrid%xblk)))
      allocate (this%grid%yblk (size(fgrid%yblk)))
      allocate (this%grid%xloc (size(fgrid%xloc)))
      allocate (this%grid%yloc (size(fgrid%yloc)))

      this%grid%xblk = fgrid%xblk
      this%grid%yblk = fgrid%yblk
      this%grid%xloc = fgrid%xloc
      this%grid%yloc = fgrid%yloc

      this%npset = pixelset%nset

      IF (p_is_worker) THEN

         allocate (afrac (pixelset%nset))
         allocate (gfrom (pixelset%nset))

         allocate (ys (pixel%nlat))
         allocate (yn (pixel%nlat))
         allocate (xw (pixel%nlon))
         allocate (xe (pixel%nlon))

         DO ilat = 1, pixel%nlat
            ys(ilat) = find_nearest_south (pixel%lat_s(ilat), fgrid%nlat, fgrid%lat_s)
            yn(ilat) = find_nearest_north (pixel%lat_n(ilat), fgrid%nlat, fgrid%lat_n)
         ENDDO

         DO ilon = 1, pixel%nlon
            xw(ilon) = find_nearest_west (pixel%lon_w(ilon), fgrid%nlon, fgrid%lon_w)
            xe(ilon) = find_nearest_east (pixel%lon_e(ilon), fgrid%nlon, fgrid%lon_e)
         ENDDO

         allocate (list_lat (fgrid%nlat))
         DO iy = 1, fgrid%nlat
            allocate (list_lat(iy)%val (100))
         ENDDO

         allocate (ng_lat (fgrid%nlat))
         ng_lat(:) = 0

         DO iset = 1, pixelset%nset

            ie = pixelset%ielm(iset)
            npxl = pixelset%ipxend(iset) - pixelset%ipxstt(iset) + 1

            allocate (afrac(iset)%val (npxl))
            allocate (gfrom(iset)%ilat(npxl))
            allocate (gfrom(iset)%ilon(npxl))

            gfrom(iset)%ng = 0
            DO ipxl = pixelset%ipxstt(iset), pixelset%ipxend(iset)

               ilat = mesh(ie)%ilat(ipxl)
               ilon = mesh(ie)%ilon(ipxl)

               DO iy = ys(ilat), yn(ilat), fgrid%yinc

                  lat_s = max(fgrid%lat_s(iy), pixel%lat_s(ilat))
                  lat_n = min(fgrid%lat_n(iy), pixel%lat_n(ilat))

                  IF ((lat_n-lat_s) < 1.0e-6_r8) THEN
                     CYCLE
                  ENDIF

                  ix = xw(ilon)
                  DO WHILE (.true.)

                     IF (ix == xw(ilon)) THEN
                        lon_w = pixel%lon_w(ilon)
                     ELSE
                        lon_w = fgrid%lon_w(ix)
                     ENDIF

                     IF (ix == xe(ilon)) THEN
                        lon_e = pixel%lon_e(ilon)
                     ELSE
                        lon_e = fgrid%lon_e(ix)
                     ENDIF

                     IF (lon_e > lon_w) THEN
                        IF ((lon_e-lon_w) < 1.0e-6_r8) THEN
                           IF (ix == xe(ilon))  EXIT
                           ix = mod(ix,fgrid%nlon) + 1
                           CYCLE
                        ENDIF
                     ELSE
                        IF ((lon_e+360.0_r8-lon_w) < 1.0e-6_r8) THEN
                           IF (ix == xe(ilon))  EXIT
                           ix = mod(ix,fgrid%nlon) + 1
                           CYCLE
                        ENDIF
                     ENDIF

                     area = areaquad (lat_s, lat_n, lon_w, lon_e)

                     CALL insert_into_sorted_list2 ( ix, iy, &
                        gfrom(iset)%ng, gfrom(iset)%ilon, gfrom(iset)%ilat, &
                        iloc, is_new)

                     IF (is_new) THEN
                        IF (iloc < gfrom(iset)%ng) THEN
                           afrac(iset)%val(iloc+1:gfrom(iset)%ng) &
                              = afrac(iset)%val(iloc:gfrom(iset)%ng-1)
                        ENDIF

                        afrac(iset)%val(iloc) = area
                     ELSE
                        afrac(iset)%val(iloc) = afrac(iset)%val(iloc) + area
                     ENDIF

                     IF (gfrom(iset)%ng == size(gfrom(iset)%ilat)) THEN
                        CALL expand_list (gfrom(iset)%ilat, 0.2_r8)
                        CALL expand_list (gfrom(iset)%ilon, 0.2_r8)
                        CALL expand_list (afrac(iset)%val,  0.2_r8)
                     ENDIF

                     CALL insert_into_sorted_list1 ( &
                        ix, ng_lat(iy), list_lat(iy)%val, iloc)

                     IF (ng_lat(iy) == size(list_lat(iy)%val)) THEN
                        CALL expand_list (list_lat(iy)%val, 0.2_r8)
                     ENDIF

                     IF (ix == xe(ilon))  EXIT
                     ix = mod(ix,fgrid%nlon) + 1
                  ENDDO
               ENDDO

            ENDDO
         ENDDO

         deallocate (ys)
         deallocate (yn)
         deallocate (xw)
         deallocate (xe)

         ng_all = sum(ng_lat)
         allocate (xlist(ng_all))
         allocate (ylist(ng_all))

         ig = 0
         DO iy = 1, fgrid%nlat
            IF (ng_lat(iy) > 0) THEN
               DO ix = 1, ng_lat(iy)
                  ig = ig + 1
                  xlist(ig) = list_lat(iy)%val(ix)
                  ylist(ig) = iy
               ENDDO
            ENDIF
         ENDDO

         deallocate (ng_lat)
         DO iy = 1, fgrid%nlat
            deallocate (list_lat(iy)%val)
         ENDDO
         deallocate (list_lat)

#ifdef USEMPI
         allocate (ipt (ng_all))
         allocate (msk (ng_all))
         DO ig = 1, ng_all
            xblk = fgrid%xblk(xlist(ig))
            yblk = fgrid%yblk(ylist(ig))
            ipt(ig) = gblock%pio(xblk,yblk)
         ENDDO
#endif

         IF (allocated(this%glist)) deallocate(this%glist)
         allocate (this%glist (0:p_np_io-1))
         DO iproc = 0, p_np_io-1
#ifdef USEMPI
            msk = (ipt == p_address_io(iproc))
            ng  = count(msk)
#else
            ng  = ng_all
#endif

            allocate (this%glist(iproc)%ilat (ng))
            allocate (this%glist(iproc)%ilon (ng))

            this%glist(iproc)%ng = 0
         ENDDO

         DO ig = 1, ng_all
#ifdef USEMPI
            iproc = p_itis_io(ipt(ig))
#else
            iproc = 0
#endif

            this%glist(iproc)%ng = this%glist(iproc)%ng + 1

            ng = this%glist(iproc)%ng
            this%glist(iproc)%ilon(ng) = xlist(ig)
            this%glist(iproc)%ilat(ng) = ylist(ig)
         ENDDO

#ifdef USEMPI
         deallocate (ipt)
         deallocate (msk)
#endif

         IF (allocated(this%address)) deallocate(this%address)
         IF (allocated(this%olparea)) deallocate(this%olparea)
         allocate (this%address (pixelset%nset))
         allocate (this%olparea (pixelset%nset))

         DO iset = 1, pixelset%nset
            ng = gfrom(iset)%ng
            allocate (this%address(iset)%val (2,ng))
            allocate (this%olparea(iset)%val (ng))

            this%olparea(iset)%val = afrac(iset)%val(1:ng)

            IF (present(pctpset)) THEN
               this%olparea(iset)%val = this%olparea(iset)%val * pctpset(iset)
            ENDIF

            DO ig = 1, gfrom(iset)%ng
               ilon = gfrom(iset)%ilon(ig)
               ilat = gfrom(iset)%ilat(ig)
               xblk = fgrid%xblk(ilon)
               yblk = fgrid%yblk(ilat)

#ifdef USEMPI
               iproc = p_itis_io(gblock%pio(xblk,yblk))
#else
               iproc = 0
#endif

               this%address(iset)%val(1,ig) = iproc
               this%address(iset)%val(2,ig) = find_in_sorted_list2 ( &
                  ilon, ilat, this%glist(iproc)%ng, this%glist(iproc)%ilon, this%glist(iproc)%ilat)
            ENDDO
         ENDDO

         deallocate (xlist)
         deallocate (ylist)

         DO iset = 1, pixelset%nset
            deallocate (afrac(iset)%val )
            deallocate (gfrom(iset)%ilon)
            deallocate (gfrom(iset)%ilat)
         ENDDO

         deallocate (afrac)
         deallocate (gfrom)

#ifdef USEMPI
         DO iproc = 0, p_np_io-1
            idest = p_address_io(iproc)
            smesg = (/p_iam_glb, this%glist(iproc)%ng/)

            CALL mpi_send (smesg, 2, MPI_INTEGER, &
               idest, mpi_tag_mesg, p_comm_glb, p_err)
 
            IF (this%glist(iproc)%ng > 0) THEN
               CALL mpi_send (this%glist(iproc)%ilon, this%glist(iproc)%ng, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err)
               CALL mpi_send (this%glist(iproc)%ilat, this%glist(iproc)%ng, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err)
            ENDIF
         ENDDO
#endif

      ENDIF

#ifdef USEMPI
      IF (p_is_io) THEN

         IF (allocated(this%glist)) deallocate(this%glist)
         allocate (this%glist (0:p_np_worker-1))

         DO iworker = 0, p_np_worker-1

            CALL mpi_recv (rmesg, 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = rmesg(1)
            nrecv = rmesg(2)
            iproc = p_itis_worker(isrc)

            this%glist(iproc)%ng = nrecv

            IF (nrecv > 0) THEN
               allocate (this%glist(iproc)%ilon (nrecv))
               allocate (this%glist(iproc)%ilat (nrecv))

               CALL mpi_recv (this%glist(iproc)%ilon, nrecv, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL mpi_recv (this%glist(iproc)%ilat, nrecv, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
            ENDIF
         ENDDO
      ENDIF
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE mapping_pset2grid_build


   !-----------------------------------------------------
   SUBROUTINE map_p2g_2d (this, pdata, gdata, spv, msk)

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_SPMD_Task
   IMPLICIT NONE

   class (mapping_pset2grid_type) :: this

   real(r8), intent(in) :: pdata(:)
   type(block_data_real8_2d), intent(inout) :: gdata

   real(r8), intent(in), optional :: spv
   logical,  intent(in), optional :: msk(:)

   ! Local variables
   integer :: iproc, idest, isrc
   integer :: ig, ilon, ilat, xblk, yblk, xloc, yloc, iloc, iset

   real(r8), allocatable :: gbuff(:)
   type(pointer_real8_1d), allocatable :: pbuff(:)

      IF (p_is_worker) THEN

         allocate (pbuff (0:p_np_io-1))

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               allocate (pbuff(iproc)%val (this%glist(iproc)%ng))

               IF (present(spv)) THEN
                  pbuff(iproc)%val(:) = spv
               ELSE
                  pbuff(iproc)%val(:) = 0.0
               ENDIF
            ENDIF
         ENDDO

         DO iset = 1, this%npset

            IF (present(spv)) THEN
               IF (pdata(iset) == spv) CYCLE
            ENDIF

            IF (present(msk)) THEN
               IF (.not. msk(iset)) CYCLE
            ENDIF

            DO ig = 1, size(this%olparea(iset)%val)
               iproc = this%address(iset)%val(1,ig)
               iloc  = this%address(iset)%val(2,ig)

               IF (present(spv)) THEN
                  IF (pbuff(iproc)%val(iloc) /= spv) THEN
                     pbuff(iproc)%val(iloc) = pbuff(iproc)%val(iloc) &
                        + pdata(iset) * this%olparea(iset)%val(ig)
                  ELSE
                     pbuff(iproc)%val(iloc) = &
                        pdata(iset) * this%olparea(iset)%val(ig)
                  ENDIF
               ELSE
                  pbuff(iproc)%val(iloc) = pbuff(iproc)%val(iloc) &
                     + pdata(iset) * this%olparea(iset)%val(ig)
               ENDIF
            ENDDO
         ENDDO

#ifdef USEMPI
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               idest = p_address_io(iproc)
               CALL mpi_send (pbuff(iproc)%val, this%glist(iproc)%ng, MPI_DOUBLE, &
                  idest, mpi_tag_data, p_comm_glb, p_err)
            ENDIF
         ENDDO
#endif

      ENDIF

      IF (p_is_io) THEN

         IF (present(spv)) THEN
            CALL flush_block_data (gdata, spv)
         ELSE
            CALL flush_block_data (gdata, 0.0_r8)
         ENDIF

         DO iproc = 0, p_np_worker-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (gbuff (this%glist(iproc)%ng))

#ifdef USEMPI
               isrc = p_address_worker(iproc)
               CALL mpi_recv (gbuff, this%glist(iproc)%ng, MPI_DOUBLE, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
#else
               gbuff = pbuff(0)%val
#endif

               DO ig = 1, this%glist(iproc)%ng
                  IF (present(spv)) THEN
                     IF (gbuff(ig) /= spv) THEN
                        ilon = this%glist(iproc)%ilon(ig)
                        ilat = this%glist(iproc)%ilat(ig)
                        xblk = this%grid%xblk (ilon)
                        yblk = this%grid%yblk (ilat)
                        xloc = this%grid%xloc (ilon)
                        yloc = this%grid%yloc (ilat)

                        IF (gdata%blk(xblk,yblk)%val(xloc,yloc) /= spv) THEN
                           gdata%blk(xblk,yblk)%val(xloc,yloc) = &
                              gdata%blk(xblk,yblk)%val(xloc,yloc) + gbuff(ig)
                        ELSE
                           gdata%blk(xblk,yblk)%val(xloc,yloc) = gbuff(ig)
                        ENDIF
                     ENDIF
                  ELSE
                     ilon = this%glist(iproc)%ilon(ig)
                     ilat = this%glist(iproc)%ilat(ig)
                     xblk = this%grid%xblk (ilon)
                     yblk = this%grid%yblk (ilat)
                     xloc = this%grid%xloc (ilon)
                     yloc = this%grid%yloc (ilat)

                     gdata%blk(xblk,yblk)%val(xloc,yloc) = &
                        gdata%blk(xblk,yblk)%val(xloc,yloc) + gbuff(ig)
                  ENDIF
               ENDDO

               deallocate (gbuff)
            ENDIF
         ENDDO

      ENDIF

      IF (p_is_worker) THEN
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               deallocate (pbuff(iproc)%val)
            ENDIF
         ENDDO
         deallocate (pbuff)
      ENDIF


   END SUBROUTINE map_p2g_2d

   !-----------------------------------------------------
   SUBROUTINE map_p2g_3d (this, pdata, gdata, spv, msk)

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_SPMD_Task
   IMPLICIT NONE

   class (mapping_pset2grid_type) :: this

   real(r8), intent(in) :: pdata(:,:)
   type(block_data_real8_3d), intent(inout) :: gdata

   real(r8), intent(in), optional :: spv
   logical,  intent(in), optional :: msk(:)

   ! Local variables
   integer :: iproc, idest, isrc
   integer :: ig, ilon, ilat, iloc, iset
   integer :: xblk, yblk, xloc, yloc
   integer :: lb1, ub1, i1

   real(r8), allocatable :: gbuff(:,:)
   type(pointer_real8_2d), allocatable :: pbuff(:)


      IF (p_is_worker) THEN

         allocate (pbuff (0:p_np_io-1))

         lb1 = lbound(pdata,1)
         ub1 = ubound(pdata,1)

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               allocate (pbuff(iproc)%val (lb1:ub1, this%glist(iproc)%ng))

               IF (present(spv)) THEN
                  pbuff(iproc)%val(:,:) = spv
               ELSE
                  pbuff(iproc)%val(:,:) = 0.0
               ENDIF
            ENDIF
         ENDDO

         DO iset = 1, this%npset
            IF (present(msk)) THEN
               IF (.not. msk(iset)) CYCLE
            ENDIF

            DO ig = 1, size(this%olparea(iset)%val)
               iproc = this%address(iset)%val(1,ig)
               iloc  = this%address(iset)%val(2,ig)

               DO i1 = lb1, ub1
                  IF (present(spv)) THEN
                     IF (pdata(i1,iset) /= spv) THEN
                        IF (pbuff(iproc)%val(i1,iloc) /= spv) THEN
                           pbuff(iproc)%val(i1,iloc) = pbuff(iproc)%val(i1,iloc) &
                              + pdata(i1,iset) * this%olparea(iset)%val(ig)
                        ELSE
                           pbuff(iproc)%val(i1,iloc) = &
                              pdata(i1,iset) * this%olparea(iset)%val(ig)
                        ENDIF
                     ENDIF
                  ELSE
                     pbuff(iproc)%val(i1,iloc) = pbuff(iproc)%val(i1,iloc) &
                        + pdata(i1,iset) * this%olparea(iset)%val(ig)
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

#ifdef USEMPI
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               idest = p_address_io(iproc)
               CALL mpi_send (pbuff(iproc)%val, &
                  (ub1-lb1+1) * this%glist(iproc)%ng, MPI_DOUBLE, &
                  idest, mpi_tag_data, p_comm_glb, p_err)

            ENDIF
         ENDDO
#endif

      ENDIF

      IF (p_is_io) THEN

         lb1 = gdata%lb1
         ub1 = gdata%ub1

         IF (present(spv)) THEN
            CALL flush_block_data (gdata, spv)
         ELSE
            CALL flush_block_data (gdata, 0.0_r8)
         ENDIF

         DO iproc = 0, p_np_worker-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (gbuff (lb1:ub1, this%glist(iproc)%ng))

#ifdef USEMPI
               isrc = p_address_worker(iproc)
               CALL mpi_recv (gbuff, &
                  (ub1-lb1+1) * this%glist(iproc)%ng, MPI_DOUBLE, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
#else
               gbuff = pbuff(0)%val
#endif

               DO ig = 1, this%glist(iproc)%ng
                  ilon = this%glist(iproc)%ilon(ig)
                  ilat = this%glist(iproc)%ilat(ig)
                  xblk = this%grid%xblk (ilon)
                  yblk = this%grid%yblk (ilat)
                  xloc = this%grid%xloc (ilon)
                  yloc = this%grid%yloc (ilat)

                  DO i1 = lb1, ub1
                     IF (present(spv)) THEN
                        IF (gbuff(i1,ig) /= spv) THEN
                           IF (gdata%blk(xblk,yblk)%val(i1,xloc,yloc) /= spv) THEN
                              gdata%blk(xblk,yblk)%val(i1,xloc,yloc) = &
                                 gdata%blk(xblk,yblk)%val(i1,xloc,yloc) + gbuff(i1,ig)
                           ELSE
                              gdata%blk(xblk,yblk)%val(i1,xloc,yloc) = gbuff(i1,ig)
                           ENDIF
                        ENDIF
                     ELSE
                        gdata%blk(xblk,yblk)%val(i1,xloc,yloc) = &
                           gdata%blk(xblk,yblk)%val(i1,xloc,yloc) + gbuff(i1,ig)
                     ENDIF
                  ENDDO
               ENDDO

               deallocate (gbuff)
            ENDIF

         ENDDO

      ENDIF

      IF (p_is_worker) THEN
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               deallocate (pbuff(iproc)%val)
            ENDIF
         ENDDO
         deallocate (pbuff)
      ENDIF

   END SUBROUTINE map_p2g_3d

   !-----------------------------------------------------
   SUBROUTINE map_p2g_4d (this, pdata, gdata, spv, msk)

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_SPMD_Task
   IMPLICIT NONE

   class (mapping_pset2grid_type) :: this

   real(r8), intent(in) :: pdata(:,:,:)
   type(block_data_real8_4d), intent(inout) :: gdata

   real(r8), intent(in), optional :: spv
   logical,  intent(in), optional :: msk(:)

   ! Local variables
   integer :: iproc, idest, isrc
   integer :: ig, ilon, ilat, iloc, iset
   integer :: xblk, yblk, xloc, yloc
   integer :: lb1, ub1, i1, ndim1, lb2, ub2, i2, ndim2

   real(r8), allocatable :: gbuff(:,:,:)
   type(pointer_real8_3d), allocatable :: pbuff(:)

      IF (p_is_worker) THEN

         allocate (pbuff (0:p_np_io-1))

         lb1 = lbound(pdata,1)
         ub1 = ubound(pdata,1)
         ndim1 = ub1 - lb1 + 1

         lb2 = lbound(pdata,2)
         ub2 = ubound(pdata,2)
         ndim2 = ub2 - lb2 + 1

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               allocate (pbuff(iproc)%val (lb1:ub1, lb2:ub2, this%glist(iproc)%ng))

               IF (present(spv)) THEN
                  pbuff(iproc)%val(:,:,:) = spv
               ELSE
                  pbuff(iproc)%val(:,:,:) = 0.0
               ENDIF
            ENDIF
         ENDDO

         DO iset = 1, this%npset
            IF (present(msk)) THEN
               IF (.not. msk(iset)) CYCLE
            ENDIF

            DO ig = 1, size(this%olparea(iset)%val)
               iproc = this%address(iset)%val(1,ig)
               iloc  = this%address(iset)%val(2,ig)

               DO i1 = lb1, ub1
                  DO i2 = lb2, ub2
                     IF (present(spv)) THEN
                        IF (pdata(i1,i2,iset) /= spv) THEN
                           IF (pbuff(iproc)%val(i1,i2,iloc) /= spv) THEN
                              pbuff(iproc)%val(i1,i2,iloc) = pbuff(iproc)%val(i1,i2,iloc) &
                                 + pdata(i1,i2,iset) * this%olparea(iset)%val(ig)
                           ELSE
                              pbuff(iproc)%val(i1,i2,iloc) = &
                                 pdata(i1,i2,iset) * this%olparea(iset)%val(ig)
                           ENDIF
                        ENDIF
                     ELSE
                        pbuff(iproc)%val(i1,i2,iloc) = pbuff(iproc)%val(i1,i2,iloc) &
                           + pdata(i1,i2,iset) * this%olparea(iset)%val(ig)
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

#ifdef USEMPI
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               idest = p_address_io(iproc)
               CALL mpi_send (pbuff(iproc)%val, ndim1 * ndim2 * this%glist(iproc)%ng, MPI_DOUBLE, &
                  idest, mpi_tag_data, p_comm_glb, p_err)
            ENDIF
         ENDDO
#endif

      ENDIF

      IF (p_is_io) THEN

         lb1 = gdata%lb1
         ub1 = gdata%ub1
         ndim1 = ub1 - lb1 + 1

         lb2 = gdata%lb2
         ub2 = gdata%ub2
         ndim2 = ub2 - lb2 + 1

         IF (present(spv)) THEN
            CALL flush_block_data (gdata, spv)
         ELSE
            CALL flush_block_data (gdata, 0.0_r8)
         ENDIF

         DO iproc = 0, p_np_worker-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (gbuff (lb1:ub1, lb2:ub2, this%glist(iproc)%ng))

#ifdef USEMPI
               isrc = p_address_worker(iproc)
               CALL mpi_recv (gbuff, ndim1 * ndim2 * this%glist(iproc)%ng, MPI_DOUBLE, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
#else
               gbuff = pbuff(0)%val
#endif

               DO ig = 1, this%glist(iproc)%ng
                  ilon = this%glist(iproc)%ilon(ig)
                  ilat = this%glist(iproc)%ilat(ig)
                  xblk = this%grid%xblk (ilon)
                  yblk = this%grid%yblk (ilat)
                  xloc = this%grid%xloc (ilon)
                  yloc = this%grid%yloc (ilat)

                  DO i1 = lb1, ub1
                     DO i2 = lb2, ub2
                        IF (present(spv)) THEN
                           IF (gbuff(i1,i2,ig) /= spv) THEN
                              IF (gdata%blk(xblk,yblk)%val(i1,i2,xloc,yloc) /= spv) THEN
                                 gdata%blk(xblk,yblk)%val(i1,i2,xloc,yloc) = &
                                    gdata%blk(xblk,yblk)%val(i1,i2,xloc,yloc) + gbuff(i1,i2,ig)
                              ELSE
                                 gdata%blk(xblk,yblk)%val(i1,i2,xloc,yloc) = gbuff(i1,i2,ig)
                              ENDIF
                           ENDIF
                        ELSE
                           gdata%blk(xblk,yblk)%val(i1,i2,xloc,yloc) = &
                              gdata%blk(xblk,yblk)%val(i1,i2,xloc,yloc) + gbuff(i1,i2,ig)
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO

               deallocate (gbuff)
            ENDIF
         ENDDO
      ENDIF

      IF (p_is_worker) THEN
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               deallocate (pbuff(iproc)%val)
            ENDIF
         ENDDO
         deallocate (pbuff)
      ENDIF

   END SUBROUTINE map_p2g_4d

   !-----------------------------------------------------
   SUBROUTINE map_p2g_split_to_3d (this, pdata, settyp, typidx, gdata, spv)

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_SPMD_Task
   IMPLICIT NONE

   class (mapping_pset2grid_type) :: this

   real(r8), intent(in) :: pdata (:)
   integer , intent(in) :: settyp(:)
   integer , intent(in) :: typidx(:)
   type(block_data_real8_3d), intent(inout) :: gdata

   real(r8), intent(in) :: spv

   ! Local variables
   integer :: iproc, idest, isrc
   integer :: ig, ilon, ilat, iloc, iset, ityp, ntyps
   integer :: xblk, yblk, xloc, yloc

   real(r8), allocatable :: gbuff(:)
   type(pointer_real8_1d), allocatable :: pbuff (:)

      IF (p_is_worker) THEN
         allocate (pbuff (0:p_np_io-1))
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               allocate (pbuff(iproc)%val  (this%glist(iproc)%ng))
            ENDIF
         ENDDO
      ENDIF

      IF (p_is_io) THEN
         CALL flush_block_data (gdata, spv)
      ENDIF

      ntyps = size(typidx)

      DO ityp = 1, ntyps

         IF (p_is_worker) THEN

            DO iproc = 0, p_np_io-1
               IF (this%glist(iproc)%ng > 0) THEN
                  pbuff(iproc)%val(:) = spv
               ENDIF
            ENDDO

            DO iset = 1, this%npset
               IF ((settyp(iset) == typidx(ityp)) .and. (pdata(iset) /= spv)) THEN
                  DO ig = 1, size(this%olparea(iset)%val)
                     iproc = this%address(iset)%val(1,ig)
                     iloc  = this%address(iset)%val(2,ig)

                     IF (pbuff(iproc)%val(iloc) /= spv) THEN
                        pbuff(iproc)%val(iloc) = pbuff(iproc)%val(iloc) &
                           + pdata(iset) * this%olparea(iset)%val(ig)
                     ELSE
                        pbuff(iproc)%val(iloc) = &
                           pdata(iset) * this%olparea(iset)%val(ig)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO

#ifdef USEMPI
            DO iproc = 0, p_np_io-1
               IF (this%glist(iproc)%ng > 0) THEN
                  idest = p_address_io(iproc)
                  CALL mpi_send (pbuff(iproc)%val, this%glist(iproc)%ng, MPI_DOUBLE, &
                     idest, mpi_tag_data, p_comm_glb, p_err)
               ENDIF
            ENDDO
#endif

         ENDIF

         IF (p_is_io) THEN

            DO iproc = 0, p_np_worker-1
               IF (this%glist(iproc)%ng > 0) THEN

                  allocate (gbuff (this%glist(iproc)%ng))

#ifdef USEMPI
                  isrc = p_address_worker(iproc)
                  CALL mpi_recv (gbuff, this%glist(iproc)%ng, MPI_DOUBLE, &
                     isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
#else
                  gbuff = pbuff(0)%val
#endif

                  DO ig = 1, this%glist(iproc)%ng
                     IF (gbuff(ig) /= spv) THEN
                        ilon = this%glist(iproc)%ilon(ig)
                        ilat = this%glist(iproc)%ilat(ig)
                        xblk = this%grid%xblk (ilon)
                        yblk = this%grid%yblk (ilat)
                        xloc = this%grid%xloc (ilon)
                        yloc = this%grid%yloc (ilat)

                        IF (gdata%blk(xblk,yblk)%val(ityp,xloc,yloc) /= spv) THEN
                           gdata%blk(xblk,yblk)%val(ityp,xloc,yloc) = &
                              gdata%blk(xblk,yblk)%val(ityp,xloc,yloc) + gbuff(ig)
                        ELSE
                           gdata%blk(xblk,yblk)%val(ityp,xloc,yloc) = gbuff(ig)
                        ENDIF
                     ENDIF
                  ENDDO

                  deallocate (gbuff)
               ENDIF

            ENDDO

         ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif
      ENDDO

      IF (p_is_worker) THEN
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               deallocate (pbuff(iproc)%val)
            ENDIF
         ENDDO
         deallocate (pbuff)
      ENDIF

   END SUBROUTINE map_p2g_split_to_3d

   !-----------------------------------------------------
   SUBROUTINE mapping_pset2grid_free_mem (this)

   USE MOD_SPMD_Task
   IMPLICIT NONE

   type (mapping_pset2grid_type) :: this

   ! Local variables
   integer :: iproc, iset

      IF (allocated (this%grid%xblk))   deallocate (this%grid%xblk)
      IF (allocated (this%grid%yblk))   deallocate (this%grid%yblk)

      IF (allocated (this%grid%xloc))   deallocate (this%grid%xloc)
      IF (allocated (this%grid%yloc))   deallocate (this%grid%yloc)

      IF (p_is_io) THEN
         IF (allocated(this%glist)) THEN
            DO iproc = 0, p_np_worker-1
               IF (allocated(this%glist(iproc)%ilat)) deallocate (this%glist(iproc)%ilat)
               IF (allocated(this%glist(iproc)%ilon)) deallocate (this%glist(iproc)%ilon)
            ENDDO

            deallocate (this%glist)
         ENDIF
      ENDIF

      IF (p_is_worker) THEN
         IF (allocated(this%glist)) THEN
            DO iproc = 0, p_np_io-1
               IF (allocated(this%glist(iproc)%ilat)) deallocate (this%glist(iproc)%ilat)
               IF (allocated(this%glist(iproc)%ilon)) deallocate (this%glist(iproc)%ilon)
            ENDDO

            deallocate (this%glist)
         ENDIF

         IF (allocated(this%address)) THEN
            DO iset = 1, this%npset
               IF (allocated(this%address(iset)%val)) THEN
                  deallocate (this%address(iset)%val)
               ENDIF
            ENDDO

            deallocate (this%address)
         ENDIF

         IF (allocated(this%olparea)) THEN
            DO iset = 1, this%npset
               IF (allocated(this%olparea(iset)%val)) THEN
                  deallocate (this%olparea(iset)%val)
               ENDIF
            ENDDO

            deallocate (this%olparea)
         ENDIF
      ENDIF

   END SUBROUTINE mapping_pset2grid_free_mem

END MODULE MOD_Mapping_Pset2Grid
