#include <define.h>

MODULE MOD_Mapping_Grid2Pset

   !-----------------------------------------------------------------------
   ! DESCRIPTION:
   !
   !    Mapping data types and subroutines from gridded data to vector data
   !    defined on pixelset.
   !
   !    Notice that:
   !    1. A mapping can be built with method mapping%build.
   !    2. Area weighted mapping is carried out.     
   !    3. For 2D gridded data, dimensions are from [lon,lat] to [vector].
   !    4. For 3D gridded data, dimensions are from [d,lon,lat] to [d,vector].
   ! 
   ! Created by Shupeng Zhang, May 2023
   !-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_DataType
   IMPLICIT NONE

   !------
   TYPE :: mapping_grid2pset_type

      TYPE(grid_type) :: grid
      INTEGER :: npset

      TYPE(grid_list_type), allocatable :: glist (:)

      TYPE(pointer_int32_2d), allocatable :: address(:)
      TYPE(pointer_real8_1d), allocatable :: gweight(:)

   CONTAINS

      procedure, PUBLIC :: build => mapping_grid2pset_build

      procedure, PRIVATE :: map_aweighted_2d => map_g2p_aweighted_2d
      procedure, PRIVATE :: map_aweighted_3d => map_g2p_aweighted_3d
      generic, PUBLIC :: map_aweighted => map_aweighted_2d, map_aweighted_3d

      procedure, PUBLIC :: map_max_frenquency_2d => map_g2p_max_frequency_2d 

      final :: mapping_grid2pset_free_mem

   END TYPE mapping_grid2pset_type

!-------------------------------------------------------------------
CONTAINS

   !------------------------------------------
   SUBROUTINE mapping_grid2pset_build (this, fgrid, pixelset, gfilter, missing_value, pfilter)

      USE MOD_Precision
      USE MOD_Namelist
      USE MOD_Block
      USE MOD_Pixel
      USE MOD_Grid
      USE MOD_DataType
      USE MOD_Mesh
      USE MOD_Pixelset
      USE MOD_Utils
      USE MOD_SPMD_Task
      IMPLICIT NONE

      class (mapping_grid2pset_type) :: this

      TYPE(grid_type),     intent(in) :: fgrid
      TYPE(pixelset_type), intent(in) :: pixelset

      TYPE(block_data_real8_2d), intent(in), optional :: gfilter
      REAL(r8), intent(in),    optional :: missing_value
      LOGICAL,  intent(inout), optional :: pfilter(:)


      ! Local variables
      TYPE(pointer_real8_1d), allocatable :: afrac(:)
      TYPE(grid_list_type),   allocatable :: gfrom(:)
      TYPE(pointer_int32_1d), allocatable :: list_lat(:)
      INTEGER,  allocatable :: ng_lat(:)
      INTEGER,  allocatable :: ys(:), yn(:), xw(:), xe(:)
      INTEGER,  allocatable :: xlist(:), ylist(:)
      INTEGER,  allocatable :: ipt(:)
      LOGICAL,  allocatable :: msk(:)

      INTEGER  :: ie, iset
      INTEGER  :: ng, ig, ng_all, iloc, ng0
      INTEGER  :: npxl, ipxl, ilat, ilon
      INTEGER  :: iworker, iproc, iio, idest, isrc, nrecv
      INTEGER  :: rmesg(2), smesg(2)
      INTEGER  :: iy, ix, xblk, yblk, xloc, yloc
      REAL(r8) :: lat_s, lat_n, lon_w, lon_e, area
      LOGICAL  :: is_new

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_master) THEN
         write(*,100)
         100 format (/, 'Making mapping from grid to pixel set')
         write(*,*) fgrid%nlat, 'grids in latitude'
         write(*,*) fgrid%nlon, 'grids in longitude'
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
                     cycle
                  ENDIF

                  ix = xw(ilon)
                  DO while (.true.)

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
                           IF (ix == xe(ilon))  exit
                           ix = mod(ix,fgrid%nlon) + 1
                           cycle
                        ENDIF
                     ELSE
                        IF ((lon_e+360.0_r8-lon_w) < 1.0e-6_r8) THEN
                           IF (ix == xe(ilon))  exit
                           ix = mod(ix,fgrid%nlon) + 1
                           cycle
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

                     IF (ix == xe(ilon))  exit
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
            DO ix = 1, ng_lat(iy)
               ig = ig + 1
               xlist(ig) = list_lat(iy)%val(ix)
               ylist(ig) = iy
            ENDDO
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

         deallocate (xlist)
         deallocate (ylist)

#ifdef USEMPI
         deallocate (ipt)
         deallocate (msk)
#endif

      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN
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
      ENDIF

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

      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (present(missing_value)) THEN

         IF (p_is_io) THEN
            DO iproc = 0, p_np_worker-1
               IF (this%glist(iproc)%ng > 0) THEN

                  allocate (msk (this%glist(iproc)%ng))

                  DO ig = 1, this%glist(iproc)%ng
                     ilon = this%glist(iproc)%ilon(ig)
                     ilat = this%glist(iproc)%ilat(ig)
                     xblk = this%grid%xblk (ilon)
                     yblk = this%grid%yblk (ilat)
                     xloc = this%grid%xloc (ilon)
                     yloc = this%grid%yloc (ilat)

                     msk(ig) = gfilter%blk(xblk,yblk)%val(xloc,yloc) /= missing_value
                  ENDDO

                  IF (any(.not. msk)) THEN

                     this%glist(iproc)%ng = count(msk)

                     IF (this%glist(iproc)%ng > 0) THEN
                        allocate (xlist(this%glist(iproc)%ng))
                        allocate (ylist(this%glist(iproc)%ng))
                        xlist = pack(this%glist(iproc)%ilon, mask=msk)
                        ylist = pack(this%glist(iproc)%ilat, mask=msk)
                     ENDIF

                     deallocate (this%glist(iproc)%ilon)
                     deallocate (this%glist(iproc)%ilat)

                     IF (this%glist(iproc)%ng > 0) THEN
                        allocate (this%glist(iproc)%ilon(this%glist(iproc)%ng))
                        allocate (this%glist(iproc)%ilat(this%glist(iproc)%ng))
                        this%glist(iproc)%ilon = xlist
                        this%glist(iproc)%ilat = ylist
                     ENDIF

                     IF (allocated(xlist)) deallocate(xlist)
                     IF (allocated(ylist)) deallocate(ylist)
                  ENDIF

                  deallocate(msk)
               ENDIF
            ENDDO
         ENDIF

#ifdef USEMPI
         IF (p_is_io) THEN
            DO iworker = 0, p_np_worker-1

               idest = p_address_worker(iworker)
               smesg = (/p_iam_glb, this%glist(iworker)%ng/)
               CALL mpi_send (smesg, 2, MPI_INTEGER, &
                  idest, mpi_tag_mesg, p_comm_glb, p_err)

               IF (this%glist(iworker)%ng > 0) THEN
                  CALL mpi_send (this%glist(iworker)%ilon, this%glist(iworker)%ng, MPI_INTEGER, &
                     idest, mpi_tag_data, p_comm_glb, p_err)
                  CALL mpi_send (this%glist(iworker)%ilat, this%glist(iworker)%ng, MPI_INTEGER, &
                     idest, mpi_tag_data, p_comm_glb, p_err)
               ENDIF
            ENDDO
         ENDIF

         IF (p_is_worker) THEN
            DO iio = 0, p_np_io-1

               CALL mpi_recv (rmesg, 2, MPI_INTEGER, &
                  MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

               isrc  = rmesg(1)
               nrecv = rmesg(2)
               iproc = p_itis_io(isrc)

               this%glist(iproc)%ng = nrecv

               IF (allocated(this%glist(iproc)%ilon)) deallocate(this%glist(iproc)%ilon)
               IF (allocated(this%glist(iproc)%ilat)) deallocate(this%glist(iproc)%ilat)

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

         CALL mpi_barrier (p_comm_glb, p_err)
#endif

         IF (p_is_worker) THEN
            DO iset = 1, pixelset%nset

               allocate (msk(gfrom(iset)%ng))

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
                  msk(ig) = find_in_sorted_list2 (ilon, ilat, this%glist(iproc)%ng, &
                     this%glist(iproc)%ilon, this%glist(iproc)%ilat) > 0

               ENDDO

               IF (present(pfilter)) THEN
                  pfilter(iset) = any(msk)
               ENDIF

               ng0 = gfrom(iset)%ng
               gfrom(iset)%ng = count(msk)
               IF (any(msk) .and. any(.not. msk)) THEN
                  ng = gfrom(iset)%ng
                  gfrom(iset)%ilon(1:ng) = pack(gfrom(iset)%ilon(1:ng0), mask = msk)
                  gfrom(iset)%ilat(1:ng) = pack(gfrom(iset)%ilat(1:ng0), mask = msk)
                  afrac(iset)%val (1:ng) = pack(afrac(iset)%val (1:ng0), mask = msk)
               ENDIF

               deallocate (msk)
            ENDDO
         ENDIF
      ENDIF

      IF (p_is_worker) THEN

         IF (allocated(this%address)) deallocate(this%address)
         IF (allocated(this%gweight)) deallocate(this%gweight)
         allocate (this%address (pixelset%nset))
         allocate (this%gweight (pixelset%nset))

         DO iset = 1, pixelset%nset

            ng = gfrom(iset)%ng
            IF (ng > 0) THEN
               allocate (this%address(iset)%val (2,ng))
               allocate (this%gweight(iset)%val (ng))

               IF (sum(afrac(iset)%val(1:ng)) < 1.0e-12) THEN
                  this%gweight(iset)%val = 1.0_r8 / ng
               ELSE
                  this%gweight(iset)%val &
                     = afrac(iset)%val(1:ng) / sum(afrac(iset)%val(1:ng))
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
            ENDIF
         ENDDO

         DO iset = 1, pixelset%nset
            deallocate (afrac(iset)%val )
            deallocate (gfrom(iset)%ilon)
            deallocate (gfrom(iset)%ilat)
         ENDDO

         deallocate (afrac)
         deallocate (gfrom)

      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE mapping_grid2pset_build

   !-----------------------------------------------------
   SUBROUTINE map_g2p_aweighted_2d (this, gdata, pdata)

      USE MOD_Precision
      USE MOD_Grid
      USE MOD_Pixelset
      USE MOD_DataType
      USE MOD_SPMD_Task
      USE MOD_Vars_Global, only : spval
      IMPLICIT NONE

      class (mapping_grid2pset_type) :: this

      TYPE(block_data_real8_2d), intent(in) :: gdata
      REAL(r8), intent(out) :: pdata(:)

      ! Local variables
      INTEGER :: iproc, idest, isrc
      INTEGER :: ig, ilon, ilat, xblk, yblk, xloc, yloc, iloc, iset

      REAL(r8), allocatable :: gbuff(:)
      TYPE(pointer_real8_1d), allocatable :: pbuff(:)

      IF (p_is_io) THEN

         DO iproc = 0, p_np_worker-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (gbuff (this%glist(iproc)%ng))

               DO ig = 1, this%glist(iproc)%ng
                  ilon = this%glist(iproc)%ilon(ig)
                  ilat = this%glist(iproc)%ilat(ig)
                  xblk = this%grid%xblk (ilon)
                  yblk = this%grid%yblk (ilat)
                  xloc = this%grid%xloc (ilon)
                  yloc = this%grid%yloc (ilat)

                  gbuff(ig) = gdata%blk(xblk,yblk)%val(xloc,yloc)

               ENDDO

#ifdef USEMPI
               idest = p_address_worker(iproc)
               CALL mpi_send (gbuff, this%glist(iproc)%ng, MPI_DOUBLE, &
                  idest, mpi_tag_data, p_comm_glb, p_err)

               deallocate (gbuff)
#endif
            ENDIF
         ENDDO

      ENDIF

      IF (p_is_worker) THEN

         allocate (pbuff (0:p_np_io-1))

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (pbuff(iproc)%val (this%glist(iproc)%ng))

#ifdef USEMPI
               isrc = p_address_io(iproc)
               CALL mpi_recv (pbuff(iproc)%val, this%glist(iproc)%ng, MPI_DOUBLE, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
#else
               pbuff(0)%val = gbuff
               deallocate (gbuff)
#endif
            ENDIF
         ENDDO

         DO iset = 1, this%npset
            IF (allocated(this%gweight(iset)%val)) THEN
               pdata(iset) = 0._r8
               DO ig = 1, size(this%gweight(iset)%val)
                  iproc = this%address(iset)%val(1,ig)
                  iloc  = this%address(iset)%val(2,ig)

                  pdata(iset) = pdata(iset) &
                     + pbuff(iproc)%val(iloc) * this%gweight(iset)%val(ig)
               ENDDO
            ELSE
               pdata(iset) = spval
            ENDIF
         ENDDO

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               deallocate (pbuff(iproc)%val)
            ENDIF
         ENDDO

         deallocate (pbuff)

      ENDIF

   END SUBROUTINE map_g2p_aweighted_2d

   !-----------------------------------------------------
   SUBROUTINE map_g2p_aweighted_3d (this, gdata, ndim1, pdata)

      USE MOD_Precision
      USE MOD_Grid
      USE MOD_Pixelset
      USE MOD_DataType
      USE MOD_SPMD_Task
      USE MOD_Vars_Global, only : spval
      IMPLICIT NONE

      class (mapping_grid2pset_type) :: this

      TYPE(block_data_real8_3d), intent(in) :: gdata
      INTEGER, intent(in) :: ndim1
      REAL(r8), intent(out) :: pdata(:,:)

      ! Local variables
      INTEGER :: iproc, idest, isrc
      INTEGER :: ig, ilon, ilat, xblk, yblk, xloc, yloc, iloc, iset

      REAL(r8), allocatable :: gbuff(:,:)
      TYPE(pointer_real8_2d), allocatable :: pbuff(:)


      IF (p_is_io) THEN

         DO iproc = 0, p_np_worker-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (gbuff (ndim1, this%glist(iproc)%ng))

               DO ig = 1, this%glist(iproc)%ng
                  ilon = this%glist(iproc)%ilon(ig)
                  ilat = this%glist(iproc)%ilat(ig)
                  xblk = this%grid%xblk (ilon)
                  yblk = this%grid%yblk (ilat)
                  xloc = this%grid%xloc (ilon)
                  yloc = this%grid%yloc (ilat)

                  gbuff(:,ig) = gdata%blk(xblk,yblk)%val(:,xloc,yloc)
               ENDDO

#ifdef USEMPI
               idest = p_address_worker(iproc)
               CALL mpi_send (gbuff, ndim1 * this%glist(iproc)%ng, MPI_DOUBLE, &
                  idest, mpi_tag_data, p_comm_glb, p_err)

               deallocate (gbuff)
#endif
            ENDIF
         ENDDO

      ENDIF

      IF (p_is_worker) THEN

         allocate (pbuff (0:p_np_io-1))

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (pbuff(iproc)%val (ndim1, this%glist(iproc)%ng))

#ifdef USEMPI
               isrc = p_address_io(iproc)
               CALL mpi_recv (pbuff(iproc)%val, ndim1 * this%glist(iproc)%ng, MPI_DOUBLE, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
#else
               pbuff(0)%val = gbuff
               deallocate (gbuff)
#endif
            ENDIF
         ENDDO

         DO iset = 1, this%npset
            IF (allocated(this%gweight(iset)%val)) THEN
               pdata(:,iset) = 0._r8
               DO ig = 1, size(this%gweight(iset)%val)
                  iproc = this%address(iset)%val(1,ig)
                  iloc  = this%address(iset)%val(2,ig)

                  pdata(:,iset) = pdata(:,iset) &
                     + pbuff(iproc)%val(:,iloc) * this%gweight(iset)%val(ig)
               ENDDO
            ELSE
               pdata(:,iset) = spval
            ENDIF
         ENDDO

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               deallocate (pbuff(iproc)%val)
            ENDIF
         ENDDO

         deallocate (pbuff)

      ENDIF

   END SUBROUTINE map_g2p_aweighted_3d

   !-----------------------------------------------------
   SUBROUTINE map_g2p_max_frequency_2d (this, gdata, pdata)

      USE MOD_Precision
      USE MOD_Grid
      USE MOD_Pixelset
      USE MOD_DataType
      USE MOD_SPMD_Task
      USE MOD_Vars_Global, only : spval
      IMPLICIT NONE

      class (mapping_grid2pset_type) :: this

      TYPE(block_data_int32_2d), intent(in) :: gdata
      integer, intent(out) :: pdata(:)

      ! Local variables
      INTEGER :: iproc, idest, isrc
      INTEGER :: ig, ilon, ilat, xblk, yblk, xloc, yloc, iloc, iset

      integer, allocatable :: gbuff(:)
      TYPE(pointer_int32_1d), allocatable :: pbuff(:)

      IF (p_is_io) THEN

         DO iproc = 0, p_np_worker-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (gbuff (this%glist(iproc)%ng))

               DO ig = 1, this%glist(iproc)%ng
                  ilon = this%glist(iproc)%ilon(ig)
                  ilat = this%glist(iproc)%ilat(ig)
                  xblk = this%grid%xblk (ilon)
                  yblk = this%grid%yblk (ilat)
                  xloc = this%grid%xloc (ilon)
                  yloc = this%grid%yloc (ilat)

                  gbuff(ig) = gdata%blk(xblk,yblk)%val(xloc,yloc)

               ENDDO

#ifdef USEMPI
               idest = p_address_worker(iproc)
               CALL mpi_send (gbuff, this%glist(iproc)%ng, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err)

               deallocate (gbuff)
#endif
            ENDIF
         ENDDO

      ENDIF

      IF (p_is_worker) THEN

         allocate (pbuff (0:p_np_io-1))

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (pbuff(iproc)%val (this%glist(iproc)%ng))

#ifdef USEMPI
               isrc = p_address_io(iproc)
               CALL mpi_recv (pbuff(iproc)%val, this%glist(iproc)%ng, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
#else
               pbuff(0)%val = gbuff
               deallocate (gbuff)
#endif
            ENDIF
         ENDDO

         DO iset = 1, this%npset
            IF (allocated(this%gweight(iset)%val)) THEN
               ig = maxloc(this%gweight(iset)%val, dim=1)
               iproc = this%address(iset)%val(1,ig)
               iloc  = this%address(iset)%val(2,ig)
               pdata(iset) = pbuff(iproc)%val(iloc)
            ELSE
               pdata(iset) = -9999
            ENDIF
         ENDDO

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               deallocate (pbuff(iproc)%val)
            ENDIF
         ENDDO

         deallocate (pbuff)

      ENDIF

   END SUBROUTINE map_g2p_max_frequency_2d
   !-----------------------------------------------------
   SUBROUTINE mapping_grid2pset_free_mem (this)

      USE MOD_SPMD_Task
      IMPLICIT NONE

      TYPE(mapping_grid2pset_type) :: this

      ! Local variables
      INTEGER :: iproc, iset

      IF (allocated (this%grid%xblk))   deallocate (this%grid%xblk)
      IF (allocated (this%grid%yblk))   deallocate (this%grid%yblk)

      IF (allocated (this%grid%xloc))   deallocate (this%grid%xloc)
      IF (allocated (this%grid%yloc))   deallocate (this%grid%yloc)

      IF (allocated(this%glist)) THEN
         DO iproc = lbound(this%glist,1), ubound(this%glist,1)
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

      IF (allocated(this%gweight)) THEN
         DO iset = 1, this%npset
            IF (allocated(this%gweight(iset)%val)) THEN
               deallocate (this%gweight(iset)%val)
            ENDIF
         ENDDO

         deallocate (this%gweight)
      ENDIF

   END SUBROUTINE mapping_grid2pset_free_mem

END MODULE MOD_Mapping_Grid2Pset
