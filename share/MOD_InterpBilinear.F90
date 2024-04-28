#include <define.h>

MODULE MOD_InterpBilinear

!-----------------------------------------------------------------------
! DESCRIPTION:
!
!    Bilinear Interpolation module.
!    
! Created by Shupeng Zhang, April 2024
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_DataType
   IMPLICIT NONE

   !------
   type :: interp_bilinear_type

      type(grid_type) :: grid
      integer :: npset
      type(grid_list_type),   allocatable :: glist (:)
      integer,  allocatable :: address(:,:,:)
      real(r8), allocatable :: weight (:,:)

   CONTAINS

      procedure, PUBLIC :: build  => interp_bilinear_build
      procedure, PUBLIC :: interp => interp_bilinear_2d

      final :: interp_bilinear_free_mem

   END type interp_bilinear_type

!-------------------------------------------------------------------
CONTAINS

   !------------------------------------------
   SUBROUTINE interp_bilinear_build (this, fgrid, pixelset, gfilter, missing_value, pfilter)

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
   USE MOD_Vars_Global, only: pi
   IMPLICIT NONE

   class (interp_bilinear_type) :: this

   type(grid_type),     intent(in) :: fgrid
   type(pixelset_type), intent(in) :: pixelset

   type(block_data_real8_2d), intent(in), optional :: gfilter
   real(r8), intent(in),    optional :: missing_value
   logical,  intent(inout), optional :: pfilter(:)


   ! Local variables
   integer,  allocatable :: ys(:), yn(:), xw(:), xe(:)
   integer,  allocatable :: xlist(:), ylist(:), ipt(:)
   
   real(r8), allocatable :: rlon_pset(:), rlat_pset(:)
   real(r8), allocatable :: nwgt(:), swgt(:), wwgt(:), ewgt(:)

   logical,  allocatable :: msk(:)

   integer  :: iset, ilat, ilon, iwest, ieast
   integer  :: nglist, iloc, ng, ig
   integer  :: iworker, iproc, iio, idest, isrc, nrecv
   integer  :: rmesg(2), smesg(2)
   integer  :: iy, ix, xblk, yblk, xloc, yloc

   real(r8) :: lon, lonw, lone, latn, lats
   real(r8) :: distn, dists, distw, diste, diffw, diffe, sumwt

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_master) THEN
         write(*,*) 
         write(*,"(A, I7, A, I7, A)") &
            'Building bilinear interpolation from grid to pixel set: ', &
            fgrid%nlat, ' grids in latitude', fgrid%nlon, ' grids in longitude'
         write(*,*)
      ENDIF

      this%grid%nlat = fgrid%nlat
      this%grid%nlon = fgrid%nlon

      allocate (this%grid%xblk (this%grid%nlon));  this%grid%xblk = fgrid%xblk
      allocate (this%grid%yblk (this%grid%nlat));  this%grid%yblk = fgrid%yblk
      allocate (this%grid%xloc (this%grid%nlon));  this%grid%xloc = fgrid%xloc
      allocate (this%grid%yloc (this%grid%nlat));  this%grid%yloc = fgrid%yloc

      IF (p_is_worker) THEN
      
         allocate (this%grid%lat_s(this%grid%nlat));  this%grid%lat_s = fgrid%lat_s
         allocate (this%grid%lat_n(this%grid%nlat));  this%grid%lat_n = fgrid%lat_n
         allocate (this%grid%lon_w(this%grid%nlon));  this%grid%lon_w = fgrid%lon_w
         allocate (this%grid%lon_e(this%grid%nlon));  this%grid%lon_e = fgrid%lon_e
         allocate (this%grid%rlon (this%grid%nlon));  CALL this%grid%set_rlon ()
         allocate (this%grid%rlat (this%grid%nlat));  CALL this%grid%set_rlat ()
      
         this%npset = pixelset%nset

         allocate (yn (this%npset))
         allocate (ys (this%npset))
         allocate (xw (this%npset))
         allocate (xe (this%npset))
         allocate (rlon_pset (this%npset))
         allocate (rlat_pset (this%npset))

         CALL pixelset%get_lonlat_radian (rlon_pset, rlat_pset)
         
         allocate (xlist(4*this%npset))
         allocate (ylist(4*this%npset))

         allocate (nwgt (this%npset))
         allocate (swgt (this%npset))
         allocate (wwgt (this%npset))
         allocate (ewgt (this%npset))

         nglist = 0
         
         DO iset = 1, this%npset

            IF (this%grid%rlat(1) > this%grid%rlat(this%grid%nlat)) THEN
               ! from north to south
               ilat = 1
               DO WHILE ((rlat_pset(iset) < this%grid%rlat(ilat)) .and. (ilat < this%grid%nlat))
                  ilat = ilat + 1
               ENDDO

               IF (rlat_pset(iset) >= this%grid%rlat(ilat)) THEN
                  yn(iset) = max(ilat-1,1)
                  ys(iset) = ilat
               ELSE
                  yn(iset) = this%grid%nlat
                  ys(iset) = this%grid%nlat
               ENDIF
            ELSE
               ! from south to north
               ilat = this%grid%nlat
               DO WHILE ((rlat_pset(iset) < this%grid%rlat(ilat)) .and. (ilat > 1))
                  ilat = ilat - 1
               ENDDO

               IF (rlat_pset(iset) >= this%grid%rlat(ilat)) THEN
                  yn(iset) = min(ilat+1,this%grid%nlat)
                  ys(iset) = ilat
               ELSE
                  yn(iset) = 1
                  ys(iset) = 1
               ENDIF
            ENDIF 

            IF (yn(iset) /= ys(iset)) THEN
               latn = this%grid%rlat(yn(iset))
               lats = this%grid%rlat(ys(iset))
               distn = arclen(rlat_pset(iset), rlon_pset(iset), latn, rlon_pset(iset))
               dists = arclen(rlat_pset(iset), rlon_pset(iset), lats, rlon_pset(iset))
               nwgt(iset) = dists/(dists+distn)
               swgt(iset) = distn/(dists+distn)
            ELSE
               nwgt(iset) = 1.0
               swgt(iset) = 0.0
            ENDIF
                     

            lon = rlon_pset(iset)*180.0/pi
            CALL normalize_longitude (lon)

            DO iwest = 1, this%grid%nlon
               lonw = this%grid%rlon(iwest) *180.0/pi
               CALL normalize_longitude (lonw)
               
               ieast = mod(iwest,this%grid%nlon) + 1
               lone  = this%grid%rlon(ieast)*180.0/pi
               CALL normalize_longitude (lone)

               IF (lon_between_floor(lon, lonw, lone)) EXIT
            ENDDO 

            xw(iset) = iwest
            xe(iset) = ieast

            ! for the case grid does not cover [-180,180)
            IF ((iwest == this%grid%nlon) .and. (this%grid%nlon > 1)) THEN
               IF (lon_between_floor( &
                  this%grid%lon_e(this%grid%nlon), lonw, this%grid%lon_w(1))) THEN

                  diffw = lon - lonw;  IF (diffw < 0) diffw = diffw + 360.0
                  diffe = lone - lon;  IF (diffe < 0) diffe = diffe + 360.0

                  IF (diffw > diffe) THEN
                     xw(iset) = ieast
                     xe(iset) = ieast
                  ELSE
                     xw(iset) = iwest
                     xe(iset) = iwest
                  ENDIF

               ENDIF
            ENDIF
            
            IF (xw(iset) /= xe(iset)) THEN
               lonw = this%grid%rlon(xw(iset))
               lone = this%grid%rlon(xe(iset))
               distw = arclen(rlat_pset(iset), rlon_pset(iset), rlat_pset(iset), lonw)
               diste = arclen(rlat_pset(iset), rlon_pset(iset), rlat_pset(iset), lone)
               wwgt(iset) = diste/(distw+diste)
               ewgt(iset) = distw/(distw+diste)
            ELSE
               wwgt(iset) = 1.0
               ewgt(iset) = 0.0
            ENDIF
            
            CALL insert_into_sorted_list2 ( xw(iset), yn(iset), nglist, xlist, ylist, iloc)
            CALL insert_into_sorted_list2 ( xe(iset), yn(iset), nglist, xlist, ylist, iloc)
            CALL insert_into_sorted_list2 ( xw(iset), ys(iset), nglist, xlist, ylist, iloc)
            CALL insert_into_sorted_list2 ( xe(iset), ys(iset), nglist, xlist, ylist, iloc)

         ENDDO

#ifdef USEMPI
         allocate (ipt (nglist))
         allocate (msk (nglist))
         DO ig = 1, nglist
            xblk = this%grid%xblk(xlist(ig))
            yblk = this%grid%yblk(ylist(ig))
            ipt(ig) = gblock%pio(xblk,yblk)
         ENDDO
#endif

         allocate (this%glist (0:p_np_io-1))
         DO iproc = 0, p_np_io-1
#ifdef USEMPI
            msk = (ipt == p_address_io(iproc))
            ng  = count(msk)
#else
            ng  = nglist
#endif

            IF (ng > 0) THEN
               allocate (this%glist(iproc)%ilat (ng))
               allocate (this%glist(iproc)%ilon (ng))
            ENDIF

            this%glist(iproc)%ng = 0
         ENDDO

         DO ig = 1, nglist
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
      ENDIF

      IF (p_is_worker) THEN

         allocate (this%address (2,4,this%npset))
         allocate (this%weight  (4,this%npset))

         DO iset = 1, pixelset%nset

            ! northwest grid
#ifdef USEMPI
            ix = xw(iset);  xblk = this%grid%xblk(ix)
            iy = yn(iset);  yblk = this%grid%yblk(iy)
            iproc = p_itis_io(gblock%pio(xblk,yblk))
#else
            iproc = 0
#endif
            this%address(1,1,iset) = iproc
            IF (this%glist(iproc)%ng > 0) THEN
               this%address(2,1,iset) = find_in_sorted_list2 ( ix, iy, &
                  this%glist(iproc)%ng, this%glist(iproc)%ilon, this%glist(iproc)%ilat)
            ELSE
               this%address(2,1,iset) = -1
            ENDIF

            IF (this%address(2,1,iset) > 0) THEN
               this%weight(1,iset) = nwgt(iset) * wwgt(iset)
            ELSE
               this%weight(1,iset) = 0
            ENDIF

            ! northeast grid
#ifdef USEMPI
            ix = xe(iset);  xblk = this%grid%xblk(ix)
            iy = yn(iset);  yblk = this%grid%yblk(iy)
            iproc = p_itis_io(gblock%pio(xblk,yblk))
#else
            iproc = 0
#endif
            this%address(1,2,iset) = iproc
            IF (this%glist(iproc)%ng > 0) THEN
               this%address(2,2,iset) = find_in_sorted_list2 ( ix, iy, &
                  this%glist(iproc)%ng, this%glist(iproc)%ilon, this%glist(iproc)%ilat)
            ELSE
               this%address(2,2,iset) = -1
            ENDIF

            IF (this%address(2,2,iset) > 0) THEN
               this%weight(2,iset) = nwgt(iset) * ewgt(iset)
            ELSE
               this%weight(2,iset) = 0
            ENDIF

            ! southwest
#ifdef USEMPI
            ix = xw(iset);  xblk = this%grid%xblk(ix)
            iy = ys(iset);  yblk = this%grid%yblk(iy)
            iproc = p_itis_io(gblock%pio(xblk,yblk))
#else
            iproc = 0
#endif
            this%address(1,3,iset) = iproc
            IF (this%glist(iproc)%ng > 0) THEN
               this%address(2,3,iset) = find_in_sorted_list2 ( ix, iy, &
                  this%glist(iproc)%ng, this%glist(iproc)%ilon, this%glist(iproc)%ilat)
            ELSE
               this%address(2,3,iset) = -1
            ENDIF

            IF (this%address(2,3,iset) > 0) THEN
               this%weight(3,iset) = swgt(iset) * wwgt(iset)
            ELSE
               this%weight(3,iset) = 0
            ENDIF

            ! southeast
#ifdef USEMPI
            ix = xe(iset);  xblk = this%grid%xblk(ix)
            iy = ys(iset);  yblk = this%grid%yblk(iy)
            iproc = p_itis_io(gblock%pio(xblk,yblk))
#else
            iproc = 0
#endif
            this%address(1,4,iset) = iproc
            IF (this%glist(iproc)%ng > 0) THEN
               this%address(2,4,iset) = find_in_sorted_list2 ( ix, iy, &
                  this%glist(iproc)%ng, this%glist(iproc)%ilon, this%glist(iproc)%ilat)
            ELSE
               this%address(2,4,iset) = -1
            ENDIF

            IF (this%address(2,4,iset) > 0) THEN
               this%weight(4,iset) = swgt(iset) * ewgt(iset)
            ELSE
               this%weight(4,iset) = 0
            ENDIF

            sumwt = sum(this%weight(:,iset))
            IF (sumwt > 0) THEN
               this%weight(:,iset) = this%weight(:,iset) / sumwt
            ENDIF

         ENDDO

         IF (present(pfilter)) THEN
            pfilter = sum(this%weight, dim=1) > 0
         ENDIF

      ENDIF

      IF (allocated(this%grid%lat_s)) deallocate(this%grid%lat_s)
      IF (allocated(this%grid%lat_n)) deallocate(this%grid%lat_n)
      IF (allocated(this%grid%lon_w)) deallocate(this%grid%lon_w)
      IF (allocated(this%grid%lon_e)) deallocate(this%grid%lon_e)
      IF (allocated(this%grid%rlon )) deallocate(this%grid%rlon )
      IF (allocated(this%grid%rlat )) deallocate(this%grid%rlat )
      
      IF (allocated(yn)) deallocate(yn)
      IF (allocated(ys)) deallocate(ys)
      IF (allocated(xw)) deallocate(xw)
      IF (allocated(xe)) deallocate(xe)

      IF (allocated(rlon_pset)) deallocate(rlon_pset)
      IF (allocated(rlat_pset)) deallocate(rlat_pset)
      
      IF (allocated(nwgt)) deallocate(nwgt)
      IF (allocated(swgt)) deallocate(swgt)
      IF (allocated(wwgt)) deallocate(wwgt)
      IF (allocated(ewgt)) deallocate(ewgt)


#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE interp_bilinear_build

   !-----------------------------------------------------
   SUBROUTINE interp_bilinear_2d (this, gdata, pdata)

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_Pixelset
   USE MOD_DataType
   USE MOD_SPMD_Task
   USE MOD_Vars_Global, only : spval
   IMPLICIT NONE

   class (interp_bilinear_type) :: this

   type(block_data_real8_2d), intent(in) :: gdata
   real(r8), intent(out) :: pdata(:)

   ! Local variables
   integer :: iproc, idest, isrc
   integer :: ig, ilon, ilat, xblk, yblk, xloc, yloc, iloc, iset

   real(r8), allocatable :: gbuff(:)
   type(pointer_real8_1d), allocatable :: pbuff(:)

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
            pdata(iset) = spval
            DO ig = 1, 4
               iproc = this%address(1,ig,iset)
               iloc  = this%address(2,ig,iset)
               IF (iloc > 0) THEN
                  IF (pdata(iset) == spval) THEN
                     pdata(iset) = pbuff(iproc)%val(iloc) * this%weight(ig,iset)
                  ELSE
                     pdata(iset) = pdata(iset) + pbuff(iproc)%val(iloc) * this%weight(ig,iset)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               deallocate (pbuff(iproc)%val)
            ENDIF
         ENDDO

         deallocate (pbuff)

      ENDIF

   END SUBROUTINE interp_bilinear_2d

   !-----------------------------------------------------
   SUBROUTINE interp_bilinear_free_mem (this)

   USE MOD_SPMD_Task
   IMPLICIT NONE

   type(interp_bilinear_type) :: this

   ! Local variables
   integer :: iproc

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

      IF (allocated(this%address)) deallocate (this%address)
      IF (allocated(this%weight )) deallocate (this%weight )

   END SUBROUTINE interp_bilinear_free_mem

END MODULE MOD_InterpBilinear
