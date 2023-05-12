#include <define.h>  

#ifdef SrfdataDiag
MODULE mod_srfdata_diag

   USE mod_grid
   USE mod_mapping_pset2grid
   
   IMPLICIT NONE
      
   ! PUBLIC variables and subroutines
   type(grid_type) :: gdiag

   TYPE(mapping_pset2grid_type) :: m_patch2diag
#ifdef PFT_CLASSIFICATION
   TYPE(mapping_pset2grid_type) :: m_pft2diag
#endif
#ifdef URBAN_MODEL
   TYPE(mapping_pset2grid_type) :: m_urb2diag
#endif

   PUBLIC :: srfdata_diag_init
   PUBLIC :: srfdata_map_and_write

   ! PRIVATE
   TYPE(grid_concat_type), PRIVATE :: srf_concat
   INTEGER, PRIVATE :: srf_data_id

CONTAINS

   ! ------ SUBROUTINE ------
   SUBROUTINE srfdata_diag_init (dir_landdata)

      USE spmd_task
      USE mod_landpatch
#ifdef PFT_CLASSIFICATION
      USE mod_landpft
#endif
#ifdef URBAN_MODEL
      USE mod_landurban
#endif

      IMPLICIT NONE

      CHARACTER(len=256), intent(in) :: dir_landdata

      ! Local Variables
      CHARACTER(len=256) :: landdir, landname
      INTEGER :: ityp
      INTEGER :: typindex(N_land_classification+1)

      landdir = trim(dir_landdata) // '/diag/'
      IF (p_is_master) THEN
         CALL system('mkdir -p ' // trim(adjustl(landdir)))
      ENDIF

      call srf_concat%set (gdiag)

#ifndef CROP
      CALL m_patch2diag%build (landpatch, gdiag)
#else
      CALL m_patch2diag%build (landpatch, gdiag, pctcrop)
#endif

#ifdef PFT_CLASSIFICATION
      CALL m_pft2diag%build (landpft, gdiag)
#endif

#ifdef URBAN_MODEL
      CALL m_urb2diag%build (landurban, gdiag)
#endif

      srf_data_id = 666

      typindex = (/(ityp, ityp = 0, N_land_classification)/)
      landname = trim(dir_landdata)//'/diag/patchfrac_elm.nc'
      CALL srfdata_map_and_write (elm_patch%subfrc, landpatch%settyp, typindex, m_patch2diag, &
         -1.0e36_r8, landname, 'patchfrac_elm', compress = 0, write_mode = 'one')

#ifdef CATCHMENT
      typindex = (/(ityp, ityp = 0, N_land_classification)/)
      landname = trim(dir_landdata)//'/diag/patchfrac_hru.nc'
      CALL srfdata_map_and_write (hru_patch%subfrc, landpatch%settyp, typindex, m_patch2diag, &
         -1.0e36_r8, landname, 'patchfrac_hru', compress = 0, write_mode = 'one')
#endif

   END SUBROUTINE srfdata_diag_init

   ! ------ SUBROUTINE ------
   subroutine srfdata_map_and_write ( &
         vsrfdata, settyp, typindex, m_srf, spv, filename, dataname, &
         compress, write_mode)

      use spmd_task
      use mod_namelist
      use mod_block
      use mod_grid
      USE mod_data_type
      USE ncio_serial
      implicit none

      REAL(r8), intent(in) :: vsrfdata (:)
      INTEGER , intent(in) :: settyp   (:)
      INTEGER , intent(in) :: typindex (:)

      TYPE(mapping_pset2grid_type), intent(in) :: m_srf

      REAL(r8), intent(in) :: spv
      
      character (len=*), intent(in) :: filename
      character (len=*), intent(in) :: dataname
      integer, intent(in) :: compress

      character (len=*), intent(in), optional :: write_mode

      ! Local variables
      type(block_data_real8_3d) :: wdata, sumwt 
      REAL(r8), allocatable :: vecone (:)

      CHARACTER(len=10) :: wmode
      integer :: iblkme, ib, jb, iblk, jblk, idata, ixseg, iyseg
      integer :: ntyps, xcnt, ycnt, xbdsp, ybdsp, xgdsp, ygdsp 
      integer :: rmesg(3), smesg(3), isrc
      character(len=256) :: fileblock
      real(r8), allocatable :: rbuf(:,:,:), sbuf(:,:,:), vdata(:,:,:)
      LOGICAL :: fexists

      IF (present(write_mode)) THEN
         wmode = trim(write_mode)
      ELSE
         wmode = 'one'
      ENDIF

      ntyps = size(typindex)

      IF (p_is_io) THEN
         call allocate_block_data (gdiag, sumwt, ntyps)
         call allocate_block_data (gdiag, wdata, ntyps)
      ENDIF

      IF (p_is_worker) THEN
         IF (size(vsrfdata) > 0) THEN
            allocate (vecone (size(vsrfdata)))
            vecone(:) = 1.0
         ENDIF
      ENDIF
   
      CALL m_srf%map_split (vecone  , settyp, typindex, sumwt, spv)
      CALL m_srf%map_split (vsrfdata, settyp, typindex, wdata, spv)

      IF (p_is_io) THEN
         DO iblkme = 1, gblock%nblkme 
            ib = gblock%xblkme(iblkme)
            jb = gblock%yblkme(iblkme)

            where ((sumwt%blk(ib,jb)%val > 0.) .and. (wdata%blk(ib,jb)%val /= spv))
               wdata%blk(ib,jb)%val = wdata%blk(ib,jb)%val / sumwt%blk(ib,jb)%val
            elsewhere
               wdata%blk(ib,jb)%val = spv
            end where
         ENDDO
      ENDIF

      if (trim(wmode) == 'one') then

         if (p_is_master) then
            
            allocate (vdata (ntyps, srf_concat%ginfo%nlon, srf_concat%ginfo%nlat))
            vdata(:,:,:) = spv
            
#ifdef USEMPI
            do idata = 1, srf_concat%ndatablk
            
               call mpi_recv (rmesg, 3, MPI_INTEGER, MPI_ANY_SOURCE, &
                  srf_data_id, p_comm_glb, p_stat, p_err)

               isrc  = rmesg(1)
               ixseg = rmesg(2)
               iyseg = rmesg(3)
               
               xgdsp = srf_concat%xsegs(ixseg)%gdsp
               ygdsp = srf_concat%ysegs(iyseg)%gdsp
               xcnt  = srf_concat%xsegs(ixseg)%cnt
               ycnt  = srf_concat%ysegs(iyseg)%cnt

               allocate (rbuf (ntyps,xcnt,ycnt))

               call mpi_recv (rbuf, ntyps * xcnt * ycnt, MPI_DOUBLE, &
                  isrc, srf_data_id, p_comm_glb, p_stat, p_err)

               vdata (:,xgdsp+1:xgdsp+xcnt,ygdsp+1:ygdsp+ycnt) = rbuf

               deallocate (rbuf)
            end do
#else
            do iyseg = 1, srf_concat%nyseg
               do ixseg = 1, srf_concat%nxseg
                  iblk  = srf_concat%xsegs(ixseg)%blk
                  jblk  = srf_concat%ysegs(iyseg)%blk
                  xbdsp = srf_concat%xsegs(ixseg)%bdsp
                  ybdsp = srf_concat%ysegs(iyseg)%bdsp
                  xgdsp = srf_concat%xsegs(ixseg)%gdsp
                  ygdsp = srf_concat%ysegs(iyseg)%gdsp
                  xcnt  = srf_concat%xsegs(ixseg)%cnt
                  ycnt  = srf_concat%ysegs(iyseg)%cnt

                  vdata (:,xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt) = &
                     wdata%blk(iblk,jblk)%val(:,xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)
               end do
            ENDDO
#endif

            write(*,*) 'Please check gridded data < ', trim(dataname), ' > in ', trim(filename)

            inquire (file=trim(filename), exist=fexists)
            IF (.not. fexists) THEN
               CALL ncio_create_file (filename)

               call ncio_define_dimension (filename, 'TypeIndex', ntyps) 
               call ncio_define_dimension (filename, 'lon' , srf_concat%ginfo%nlon)
               call ncio_define_dimension (filename, 'lat' , srf_concat%ginfo%nlat)

               call ncio_write_serial (filename, 'lat', srf_concat%ginfo%lat_c, 'lat')
               CALL ncio_put_attr (filename, 'lat', 'long_name', 'latitude')
               CALL ncio_put_attr (filename, 'lat', 'units', 'degrees_north')

               call ncio_write_serial (filename, 'lon', srf_concat%ginfo%lon_c, 'lon')
               CALL ncio_put_attr (filename, 'lon', 'long_name', 'longitude')
               CALL ncio_put_attr (filename, 'lon', 'units', 'degrees_east')

               call ncio_write_serial (filename, 'TypeIndex', typindex, 'TypeIndex')
            ENDIF

            call ncio_write_serial (filename, dataname, vdata, 'TypeIndex', 'lon', 'lat', compress)
               
            CALL ncio_put_attr (filename, dataname, 'missing_value', spv)

            deallocate (vdata)

         ENDIF

#ifdef USEMPI
         if (p_is_io) then

            do iyseg = 1, srf_concat%nyseg
               do ixseg = 1, srf_concat%nxseg

                  iblk = srf_concat%xsegs(ixseg)%blk
                  jblk = srf_concat%ysegs(iyseg)%blk

                  if (gblock%pio(iblk,jblk) == p_iam_glb) then

                     xbdsp = srf_concat%xsegs(ixseg)%bdsp
                     ybdsp = srf_concat%ysegs(iyseg)%bdsp
                     xcnt  = srf_concat%xsegs(ixseg)%cnt
                     ycnt  = srf_concat%ysegs(iyseg)%cnt

                     allocate (sbuf (ntyps,xcnt,ycnt))
                     sbuf = wdata%blk(iblk,jblk)%val(:,xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)

                     smesg = (/p_iam_glb, ixseg, iyseg/)
                     call mpi_send (smesg, 3, MPI_INTEGER, &
                        p_root, srf_data_id, p_comm_glb, p_err) 
                     call mpi_send (sbuf, ntyps*xcnt*ycnt, MPI_DOUBLE, &
                        p_root, srf_data_id, p_comm_glb, p_err)

                     deallocate (sbuf)
                  end if
               end do
            end do
         end if
#endif

         srf_data_id = srf_data_id + 1

      elseif (trim(wmode) == 'block') then

         if (p_is_io) then

            DO iblkme = 1, gblock%nblkme 
               iblk = gblock%xblkme(iblkme)
               jblk = gblock%yblkme(iblkme)

               if ((gdiag%xcnt(iblk) == 0) .or. (gdiag%ycnt(jblk) == 0)) cycle

               call get_filename_block (filename, iblk, jblk, fileblock)

               inquire (file=trim(filename), exist=fexists)
               IF (.not. fexists) THEN
                  CALL ncio_create_file (fileblock)
                  call ncio_define_dimension (fileblock, 'TypeIndex', ntyps)
                  CALL srf_write_grid_info   (fileblock, gdiag, iblk, jblk)
               ENDIF

               call ncio_write_serial (fileblock, dataname, &
                  wdata%blk(iblk,jblk)%val, 'TypeIndex', 'lon', 'lat', compress)
               
               CALL ncio_put_attr (fileblock, dataname, 'missing_value', spv)

            end do

         end if
      end if

      IF (allocated(vecone)) deallocate(vecone)

   end subroutine srfdata_map_and_write

!    ! ------ SUBROUTINE ------
!    subroutine srfdata_map_and_write_int ( &
!          vsrfdata, settyp, typindex, m_srf, spv, filename, dataname, &
!          compress, write_mode)

!       use spmd_task
!       use mod_namelist
!       use mod_block
!       use mod_grid
!       USE mod_data_type
!       USE ncio_serial
!       implicit none

!       INTEGER , intent(in) :: vsrfdata (:)
!       INTEGER , intent(in) :: settyp   (:)
!       INTEGER , intent(in) :: typindex (:)

!       TYPE(mapping_pset2grid_type), intent(in) :: m_srf

!       INTEGER , intent(in) :: spv
      
!       character (len=*), intent(in) :: filename
!       character (len=*), intent(in) :: dataname
!       integer, intent(in) :: compress

!       character (len=*), intent(in), optional :: write_mode

!       ! Local variables
!       type(block_data_real8_3d) :: wdata, sumwt 
!       REAL(r8), allocatable :: vecone (:)

!       CHARACTER(len=10) :: wmode
!       integer :: iblkme, ib, jb, iblk, jblk, idata, ixseg, iyseg
!       integer :: ntyps, xcnt, ycnt, xbdsp, ybdsp, xgdsp, ygdsp 
!       integer :: rmesg(3), smesg(3), isrc
!       character(len=256) :: fileblock
!       real(r8), allocatable :: rbuf(:,:,:), sbuf(:,:,:), vdata(:,:,:)
!       LOGICAL :: fexists

!       IF (present(write_mode)) THEN
!          wmode = trim(write_mode)
!       ELSE
!          wmode = 'one'
!       ENDIF

!       ntyps = size(typindex)

!       IF (p_is_io) THEN
!          call allocate_block_data (gdiag, sumwt, ntyps)
!          call allocate_block_data (gdiag, wdata, ntyps)
!       ENDIF

!       IF (p_is_worker) THEN
!          IF (size(vsrfdata) > 0) THEN
!             allocate (vecone (size(vsrfdata)))
!             vecone(:) = 1.0
!          ENDIF
!       ENDIF
   
!       CALL m_srf%map_split (vecone  , settyp, typindex, sumwt, spv)
!       CALL m_srf%map_split (vsrfdata, settyp, typindex, wdata, spv)

!       IF (p_is_io) THEN
!          DO iblkme = 1, gblock%nblkme 
!             ib = gblock%xblkme(iblkme)
!             jb = gblock%yblkme(iblkme)

!             where ((sumwt%blk(ib,jb)%val > 0.) .and. (wdata%blk(ib,jb)%val /= spv))
!                wdata%blk(ib,jb)%val = wdata%blk(ib,jb)%val / sumwt%blk(ib,jb)%val
!             elsewhere
!                wdata%blk(ib,jb)%val = spv
!             end where
!          ENDDO
!       ENDIF

!       if (trim(wmode) == 'one') then

!          if (p_is_master) then
            
!             allocate (vdata (ntyps, srf_concat%ginfo%nlon, srf_concat%ginfo%nlat))
!             vdata(:,:,:) = spv
            
! #ifdef USEMPI
!             do idata = 1, srf_concat%ndatablk
            
!                call mpi_recv (rmesg, 3, MPI_INTEGER, MPI_ANY_SOURCE, &
!                   srf_data_id, p_comm_glb, p_stat, p_err)

!                isrc  = rmesg(1)
!                ixseg = rmesg(2)
!                iyseg = rmesg(3)
               
!                xgdsp = srf_concat%xsegs(ixseg)%gdsp
!                ygdsp = srf_concat%ysegs(iyseg)%gdsp
!                xcnt  = srf_concat%xsegs(ixseg)%cnt
!                ycnt  = srf_concat%ysegs(iyseg)%cnt

!                allocate (rbuf (ntyps,xcnt,ycnt))

!                call mpi_recv (rbuf, ntyps * xcnt * ycnt, MPI_DOUBLE, &
!                   isrc, srf_data_id, p_comm_glb, p_stat, p_err)

!                vdata (:,xgdsp+1:xgdsp+xcnt,ygdsp+1:ygdsp+ycnt) = rbuf

!                deallocate (rbuf)
!             end do
! #else
!             do iyseg = 1, srf_concat%nyseg
!                do ixseg = 1, srf_concat%nxseg
!                   iblk  = srf_concat%xsegs(ixseg)%blk
!                   jblk  = srf_concat%ysegs(iyseg)%blk
!                   xbdsp = srf_concat%xsegs(ixseg)%bdsp
!                   ybdsp = srf_concat%ysegs(iyseg)%bdsp
!                   xgdsp = srf_concat%xsegs(ixseg)%gdsp
!                   ygdsp = srf_concat%ysegs(iyseg)%gdsp
!                   xcnt  = srf_concat%xsegs(ixseg)%cnt
!                   ycnt  = srf_concat%ysegs(iyseg)%cnt

!                   vdata (:,xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt) = &
!                      wdata%blk(iblk,jblk)%val(:,xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)
!                end do
!             ENDDO
! #endif

!             write(*,*) 'Please check gridded data < ', trim(dataname), ' > in ', trim(filename)

!             inquire (file=trim(filename), exist=fexists)
!             IF (.not. fexists) THEN
!                CALL ncio_create_file (filename)

!                call ncio_define_dimension (filename, 'TypeIndex', ntyps) 
!                call ncio_define_dimension (filename, 'lon' , srf_concat%ginfo%nlon)
!                call ncio_define_dimension (filename, 'lat' , srf_concat%ginfo%nlat)

!                call ncio_write_serial (filename, 'lat', srf_concat%ginfo%lat_c, 'lat')
!                CALL ncio_put_attr (filename, 'lat', 'long_name', 'latitude')
!                CALL ncio_put_attr (filename, 'lat', 'units', 'degrees_north')

!                call ncio_write_serial (filename, 'lon', srf_concat%ginfo%lon_c, 'lon')
!                CALL ncio_put_attr (filename, 'lon', 'long_name', 'longitude')
!                CALL ncio_put_attr (filename, 'lon', 'units', 'degrees_east')

!                call ncio_write_serial (filename, 'TypeIndex', typindex, 'TypeIndex')
!             ENDIF

!             call ncio_write_serial (filename, dataname, vdata, 'TypeIndex', 'lon', 'lat', compress)
               
!             CALL ncio_put_attr (filename, dataname, 'missing_value', spv)

!             deallocate (vdata)

!          ENDIF

! #ifdef USEMPI
!          if (p_is_io) then

!             do iyseg = 1, srf_concat%nyseg
!                do ixseg = 1, srf_concat%nxseg

!                   iblk = srf_concat%xsegs(ixseg)%blk
!                   jblk = srf_concat%ysegs(iyseg)%blk

!                   if (gblock%pio(iblk,jblk) == p_iam_glb) then

!                      xbdsp = srf_concat%xsegs(ixseg)%bdsp
!                      ybdsp = srf_concat%ysegs(iyseg)%bdsp
!                      xcnt  = srf_concat%xsegs(ixseg)%cnt
!                      ycnt  = srf_concat%ysegs(iyseg)%cnt

!                      allocate (sbuf (ntyps,xcnt,ycnt))
!                      sbuf = wdata%blk(iblk,jblk)%val(:,xbdsp+1:xbdsp+xcnt,ybdsp+1:ybdsp+ycnt)

!                      smesg = (/p_iam_glb, ixseg, iyseg/)
!                      call mpi_send (smesg, 3, MPI_INTEGER, &
!                         p_root, srf_data_id, p_comm_glb, p_err) 
!                      call mpi_send (sbuf, ntyps*xcnt*ycnt, MPI_DOUBLE, &
!                         p_root, srf_data_id, p_comm_glb, p_err)

!                      deallocate (sbuf)
!                   end if
!                end do
!             end do
!          end if
! #endif

!          srf_data_id = srf_data_id + 1

!       elseif (trim(wmode) == 'block') then

!          if (p_is_io) then

!             DO iblkme = 1, gblock%nblkme 
!                iblk = gblock%xblkme(iblkme)
!                jblk = gblock%yblkme(iblkme)

!                if ((gdiag%xcnt(iblk) == 0) .or. (gdiag%ycnt(jblk) == 0)) cycle

!                call get_filename_block (filename, iblk, jblk, fileblock)

!                inquire (file=trim(filename), exist=fexists)
!                IF (.not. fexists) THEN
!                   CALL ncio_create_file (fileblock)
!                   call ncio_define_dimension (fileblock, 'TypeIndex', ntyps)
!                   CALL srf_write_grid_info   (fileblock, gdiag, iblk, jblk)
!                ENDIF

!                call ncio_write_serial (fileblock, dataname, &
!                   wdata%blk(iblk,jblk)%val, 'TypeIndex', 'lon', 'lat', compress)
               
!                CALL ncio_put_attr (fileblock, dataname, 'missing_value', spv)

!             end do

!          end if
!       end if

!       IF (allocated(vecone)) deallocate(vecone)

!    end subroutine srfdata_map_and_write_int

   !------------------
   subroutine srf_write_grid_info (fileblock, grid, iblk, jblk)

      use mod_block
      use mod_grid
      USE ncio_serial
      implicit none

      character(len=*), intent(in) :: fileblock
      type (grid_type), intent(in) :: grid
      integer, intent(in) :: iblk, jblk

      ! Local variable
      integer :: yl, yu, xl, xu, nx
      real(r8), allocatable :: lat_s(:), lat_n(:), lon_w(:), lon_e(:)

      allocate (lon_w (grid%xcnt(iblk)))
      allocate (lon_e (grid%xcnt(iblk)))
      allocate (lat_s (grid%ycnt(jblk)))
      allocate (lat_n (grid%ycnt(jblk)))

      yl = grid%ydsp(jblk) + 1
      yu = grid%ydsp(jblk) + grid%ycnt(jblk)

      lat_s = grid%lat_s(yl:yu)
      lat_n = grid%lat_n(yl:yu)

      if (grid%xdsp(iblk) + grid%xcnt(iblk) > grid%nlon) then
         xl = grid%xdsp(iblk) + 1
         xu = grid%nlon
         nx = grid%nlon - grid%xdsp(iblk)
         lon_w(1:nx) = grid%lon_w(xl:xu)
         lon_e(1:nx) = grid%lon_e(xl:xu)

         xl = 1
         xu = grid%xcnt(iblk) - nx 
         lon_w(nx+1:grid%xcnt(iblk)) = grid%lon_w(xl:xu)
         lon_e(nx+1:grid%xcnt(iblk)) = grid%lon_e(xl:xu)
      else
         xl = grid%xdsp(iblk) + 1
         xu = grid%xdsp(iblk) + grid%xcnt(iblk)
         lon_w = grid%lon_w(xl:xu)
         lon_e = grid%lon_e(xl:xu)
      end if

      CALL ncio_define_dimension (fileblock, 'lat', grid%ycnt(jblk))
      CALL ncio_define_dimension (fileblock, 'lon', grid%xcnt(iblk))
      call ncio_write_serial (fileblock, 'lat_s', lat_s, 'lat')
      call ncio_write_serial (fileblock, 'lat_n', lat_n, 'lat')
      call ncio_write_serial (fileblock, 'lon_w', lon_w, 'lon')
      call ncio_write_serial (fileblock, 'lon_e', lon_e, 'lon')

   end subroutine srf_write_grid_info

END MODULE mod_srfdata_diag
#endif
