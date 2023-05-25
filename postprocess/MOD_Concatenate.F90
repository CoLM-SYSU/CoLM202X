#include <define.h>

module mod_concatenate

   use MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Block
   USE MOD_NetCDFSerial
   use netcdf
   USE MOD_Vars_Global, only : spval
   implicit none

   type :: segment_type
      integer :: blk
      integer :: cnt
      integer :: bdsp
      integer :: gdsp
   end type segment_type

   type :: grid_info_type
      integer :: nlat, nlon
      real(r8), allocatable :: lat_s(:)
      real(r8), allocatable :: lat_n(:)
      real(r8), allocatable :: lon_w(:)
      real(r8), allocatable :: lon_e(:)
      real(r8), allocatable :: lon_c(:) !grid center
      real(r8), allocatable :: lat_c(:) !grid center
   end type grid_info_type

   TYPE :: block_info_type
      integer :: ndatablk
      integer :: nxseg, nyseg
      type(segment_type), allocatable :: xsegs(:), ysegs(:)
      type(grid_info_type) :: ginfo
   CONTAINS
      final :: block_info_free_mem
   END TYPE block_info_type

   TYPE(block_info_type) :: hist_block_info

   ! --- subroutines ---
   public :: hist_concatenate_var_2d
   public :: hist_concatenate_var_3d
   public :: hist_concatenate_var_4d

   public :: hist_concatenate_time
   public :: hist_concatenate_grid_info

contains

   ! ------
   subroutine set_hist_block_info (ghist, block_info)

      use MOD_Grid
      use MOD_Block
      USE MOD_Utils
      implicit none

      type(grid_type), intent(in) :: ghist
      TYPE(block_info_type), intent(out) :: block_info

      ! Local variables
      integer :: ilat_l, ilat_u, ilat, ilatloc, jblk, iyseg
      integer :: ilon_w, ilon_e, ilon, ilonloc, iblk, ixseg

      ilat_l = findloc(ghist%yblk /= 0, .true., dim=1)
      ilat_u = findloc(ghist%yblk /= 0, .true., dim=1, back=.true.)

      block_info%ginfo%nlat = ilat_u - ilat_l + 1
      allocate (block_info%ginfo%lat_s (block_info%ginfo%nlat))
      allocate (block_info%ginfo%lat_n (block_info%ginfo%nlat))
      allocate (block_info%ginfo%lat_c (block_info%ginfo%nlat))

      block_info%nyseg = 0
      jblk  = 0
      ilatloc = 0
      do ilat = ilat_l, ilat_u
         if (ghist%yblk(ilat) /= jblk) then
            block_info%nyseg = block_info%nyseg + 1
            jblk  = ghist%yblk(ilat)
         end if

         ilatloc = ilatloc + 1
         block_info%ginfo%lat_s(ilatloc) = ghist%lat_s(ilat)
         block_info%ginfo%lat_n(ilatloc) = ghist%lat_n(ilat)
         block_info%ginfo%lat_c(ilatloc) = (ghist%lat_s(ilat)+ghist%lat_n(ilat)) * 0.5
      end do

      allocate (block_info%ysegs (block_info%nyseg))

      iyseg = 0
      jblk  = 0
      do ilat = ilat_l, ilat_u
         if (ghist%yblk(ilat) /= jblk) then
            iyseg = iyseg + 1
            jblk  = ghist%yblk(ilat)
            block_info%ysegs(iyseg)%blk  = jblk
            block_info%ysegs(iyseg)%bdsp = ghist%yloc(ilat) - 1
            block_info%ysegs(iyseg)%gdsp = ilat - ilat_l
            block_info%ysegs(iyseg)%cnt  = 1
         else
            block_info%ysegs(iyseg)%cnt  = block_info%ysegs(iyseg)%cnt + 1
         end if
      end do

      if (all(ghist%xblk > 0)) then
         ilon_w = 1
         ilon_e = ghist%nlon
      else
         ilon_w = findloc(ghist%xblk /= 0, .true., dim=1)
         do while (.true.)
            ilon = ilon_w - 1
            if (ilon == 0) ilon = ghist%nlon

            if (ghist%xblk(ilon) /= 0) then
               ilon_w = ilon
            else
               exit
            end if
         end do

         ilon_e = ilon_w
         do while (.true.)
            ilon = mod(ilon_e,ghist%nlon) + 1

            if (ghist%xblk(ilon) /= 0) then
               ilon_e = ilon
            else
               exit
            end if
         end do
      end if

      block_info%ginfo%nlon = ilon_e - ilon_w + 1
      if (block_info%ginfo%nlon <= 0) THEN
         block_info%ginfo%nlon = block_info%ginfo%nlon + ghist%nlon
      ENDIF

      allocate (block_info%ginfo%lon_w (block_info%ginfo%nlon))
      allocate (block_info%ginfo%lon_e (block_info%ginfo%nlon))
      allocate (block_info%ginfo%lon_c (block_info%ginfo%nlon))

      block_info%nxseg = 0
      ilon = ilon_w - 1
      iblk = 0
      ilonloc = 0
      do while (.true.)
         ilon = mod(ilon,ghist%nlon) + 1
         if (ghist%xblk(ilon) /= iblk) then
            block_info%nxseg = block_info%nxseg + 1
            iblk = ghist%xblk(ilon)
         end if

         ilonloc = ilonloc + 1
         block_info%ginfo%lon_w(ilonloc) = ghist%lon_w(ilon)
         block_info%ginfo%lon_e(ilonloc) = ghist%lon_e(ilon)

         block_info%ginfo%lon_c(ilonloc) = (ghist%lon_w(ilon) + ghist%lon_e(ilon)) * 0.5
         IF (ghist%lon_w(ilon) > ghist%lon_e(ilon)) THEN
            block_info%ginfo%lon_c(ilonloc) = block_info%ginfo%lon_c(ilonloc) + 180.0
            CALL normalize_longitude (block_info%ginfo%lon_c(ilonloc))
         ENDIF

         if (ilon == ilon_e) exit
      end do

      allocate (block_info%xsegs (block_info%nxseg))

      ixseg = 0
      iblk = 0
      ilon = ilon_w - 1
      ilonloc = 0
      do while (.true.)
         ilon = mod(ilon,ghist%nlon) + 1
         ilonloc = ilonloc + 1
         if (ghist%xblk(ilon) /= iblk) then
            ixseg = ixseg + 1
            iblk = ghist%xblk(ilon)
            block_info%xsegs(ixseg)%blk  = iblk
            block_info%xsegs(ixseg)%bdsp = ghist%xloc(ilon) - 1
            block_info%xsegs(ixseg)%gdsp = ilonloc - 1
            block_info%xsegs(ixseg)%cnt = 1
         else
            block_info%xsegs(ixseg)%cnt = block_info%xsegs(ixseg)%cnt + 1
         end if

         if (ilon == ilon_e) exit
      end do

   end subroutine set_hist_block_info

   ! -----------
   SUBROUTINE block_info_free_mem (this)

      IMPLICIT NONE

      TYPE(block_info_type) :: this

      IF (allocated(this%xsegs)) deallocate(this%xsegs)
      IF (allocated(this%ysegs)) deallocate(this%ysegs)

      IF (allocated(this%ginfo%lat_s)) deallocate(this%ginfo%lat_s)
      IF (allocated(this%ginfo%lat_n)) deallocate(this%ginfo%lat_n)
      IF (allocated(this%ginfo%lat_c)) deallocate(this%ginfo%lat_c)
      IF (allocated(this%ginfo%lon_w)) deallocate(this%ginfo%lon_w)
      IF (allocated(this%ginfo%lon_e)) deallocate(this%ginfo%lon_e)
      IF (allocated(this%ginfo%lon_c)) deallocate(this%ginfo%lon_c)

   END SUBROUTINE block_info_free_mem

   !------------------------------
   subroutine hist_concatenate_time (filename, dataname, timelen)

      implicit none

      character (len=*), intent(in) :: filename
      character (len=*), intent(in) :: dataname
      INTEGER, intent(out) :: timelen

      ! Local variables
      CHARACTER(len=256) :: fileblock
      integer :: ixseg, iyseg
      INTEGER, allocatable :: minutes (:)

      iyseg = 1
      DO WHILE (hist_block_info%ysegs(iyseg)%cnt <= 0)
         iyseg = iyseg + 1
      ENDDO

      ixseg = 1
      DO WHILE (hist_block_info%xsegs(ixseg)%cnt <= 0)
         ixseg = ixseg + 1
      ENDDO

      call get_filename_block ( &
         filename, hist_block_info%xsegs(ixseg)%blk, hist_block_info%ysegs(iyseg)%blk, fileblock)

      CALL ncio_read_serial (fileblock, dataname, minutes)

      timelen = size(minutes)
      CALL ncio_define_dimension (filename, 'time', timelen)
      call ncio_write_serial (filename, dataname, minutes, 'time')

      CALL ncio_put_attr (filename, dataname, 'long_name', 'time')
      CALL ncio_put_attr (filename, dataname, 'units', 'minutes since 1900-1-1 0:0:0')

   end subroutine hist_concatenate_time

   !------------------
   subroutine hist_concatenate_grid_info (filename)
      
      use MOD_NetCDFSerial
      implicit none

      character(len=*), intent(in) :: filename

      call ncio_define_dimension(filename, 'lat' , hist_block_info%ginfo%nlat)
      call ncio_define_dimension(filename, 'lon' , hist_block_info%ginfo%nlon)

      call ncio_write_serial (filename, 'lat_s', hist_block_info%ginfo%lat_s, 'lat')
      call ncio_write_serial (filename, 'lat_n', hist_block_info%ginfo%lat_n, 'lat')
      call ncio_write_serial (filename, 'lon_w', hist_block_info%ginfo%lon_w, 'lon')
      call ncio_write_serial (filename, 'lon_e', hist_block_info%ginfo%lon_e, 'lon')
      call ncio_write_serial (filename, 'lat',   hist_block_info%ginfo%lat_c, 'lat')
      call ncio_write_serial (filename, 'lon',   hist_block_info%ginfo%lon_c, 'lon')

      CALL ncio_put_attr (filename, 'lat', 'long_name', 'latitude')
      CALL ncio_put_attr (filename, 'lat', 'units', 'degrees_north')
      CALL ncio_put_attr (filename, 'lat', 'axis', 'X')

      CALL ncio_put_attr (filename, 'lon', 'long_name', 'longitude')
      CALL ncio_put_attr (filename, 'lon', 'units', 'degrees_east')
      CALL ncio_put_attr (filename, 'lon', 'axis', 'Y')

   end subroutine hist_concatenate_grid_info

#ifndef USEMPI
   ! -----
   subroutine hist_concatenate_var_2d (filename, varname, timelen, compress, longname, units)

      implicit none

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: varname
      integer, intent(in) :: timelen
      integer, intent(in) :: compress
      character (len=*), intent(in), optional :: longname
      character (len=*), intent(in), optional :: units

      ! Local variables
      CHARACTER(len=256) :: fileblock
      integer :: ixseg, iyseg, nlon, nlat, xgdsp, ygdsp, xcnt, ycnt
      real(r8), allocatable :: vdata(:,:,:), rcache(:,:,:)
      LOGICAL :: fexists

      write(*,*) 'Concatenate <', trim(varname), '> to <<', trim(filename), '>>'

      nlat = hist_block_info%ginfo%nlat
      nlon = hist_block_info%ginfo%nlon

      allocate (vdata (nlon, nlat, timelen))
      vdata(:,:,:) = spval

      do iyseg = 1, hist_block_info%nyseg
         do ixseg = 1, hist_block_info%nxseg

            call get_filename_block (filename, &
               hist_block_info%xsegs(ixseg)%blk, hist_block_info%ysegs(iyseg)%blk, fileblock)

            inquire(file=fileblock, exist=fexists)
            IF (fexists) THEN

               CALL ncio_read_serial (fileblock, varname, rcache)

               xgdsp = hist_block_info%xsegs(ixseg)%gdsp
               ygdsp = hist_block_info%ysegs(iyseg)%gdsp
               xcnt  = hist_block_info%xsegs(ixseg)%cnt
               ycnt  = hist_block_info%ysegs(iyseg)%cnt

               vdata (xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt,:) = rcache
            end if

         end do
      end do

      call ncio_write_serial (filename, varname, vdata, 'lon', 'lat', 'time', compress)

      IF (present(longname)) THEN
         CALL ncio_put_attr (filename, varname, 'long_name', longname)
      ENDIF
      IF (present(units)) THEN
         CALL ncio_put_attr (filename, varname, 'units', units)
      ENDIF

      CALL ncio_put_attr (filename, varname, 'missing_value', spval)

      IF (allocated(rcache)) deallocate(rcache)
      IF (allocated(vdata )) deallocate(vdata )

   end subroutine hist_concatenate_var_2d

   ! -----
   subroutine hist_concatenate_var_3d (filename, varname, timelen, dim1name, ndim1, compress, &
         longname, units)

      implicit none

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: varname
      INTEGER, intent(in) :: timelen
      character(len=*), intent(in) :: dim1name
      INTEGER, intent(in) :: ndim1
      integer, intent(in) :: compress
      character (len=*), intent(in), optional :: longname
      character (len=*), intent(in), optional :: units

      ! Local variables
      CHARACTER(len=256) :: fileblock
      integer :: ixseg, iyseg, nlon, nlat, xgdsp, ygdsp, xcnt, ycnt
      real(r8), allocatable :: vdata(:,:,:,:), rcache(:,:,:,:)
      LOGICAL :: fexists

      nlat = hist_block_info%ginfo%nlat
      nlon = hist_block_info%ginfo%nlon

      allocate (vdata (ndim1, nlon, nlat, timelen))
      vdata(:,:,:,:) = spval

      write(*,*) 'Concatenate <', trim(varname), '> to <<', trim(filename), '>>'

      do iyseg = 1, hist_block_info%nyseg
         do ixseg = 1, hist_block_info%nxseg

            call get_filename_block (filename, &
               hist_block_info%xsegs(ixseg)%blk, hist_block_info%ysegs(iyseg)%blk, fileblock)

            inquire(file=fileblock, exist=fexists)
            IF (fexists) THEN

               CALL ncio_read_serial (fileblock, varname, rcache)

               xgdsp = hist_block_info%xsegs(ixseg)%gdsp
               ygdsp = hist_block_info%ysegs(iyseg)%gdsp
               xcnt  = hist_block_info%xsegs(ixseg)%cnt
               ycnt  = hist_block_info%ysegs(iyseg)%cnt

               vdata (:,xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt,:) = rcache
            end if

         end do
      end do

      call ncio_write_serial (filename, varname, vdata, dim1name, 'lon', 'lat', 'time', compress)

      IF (present(longname)) THEN
         CALL ncio_put_attr (filename, varname, 'long_name', longname)
      ENDIF
      IF (present(units)) THEN
         CALL ncio_put_attr (filename, varname, 'units', units)
      ENDIF

      CALL ncio_put_attr (filename, varname, 'missing_value', spval)

      IF (allocated(rcache)) deallocate(rcache)
      IF (allocated(vdata )) deallocate(vdata )

   end subroutine hist_concatenate_var_3d

   ! -----
   subroutine hist_concatenate_var_4d ( &
         filename, varname, timelen, dim1name, dim2name, ndim1, ndim2, compress, &
         longname, units)

      implicit none

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: varname
      INTEGER, intent(in) :: timelen
      character(len=*), intent(in) :: dim1name, dim2name
      INTEGER, intent(in) :: ndim1, ndim2
      integer, intent(in) :: compress
      character (len=*), intent(in), optional :: longname
      character (len=*), intent(in), optional :: units

      ! Local variables
      CHARACTER(len=256) :: fileblock
      integer :: ixseg, iyseg, nlon, nlat, xgdsp, ygdsp, xcnt, ycnt
      real(r8), allocatable :: vdata(:,:,:,:,:), rcache(:,:,:,:,:)
      LOGICAL :: fexists

      write(*,*) 'Concatenate <', trim(varname), '> to <<', trim(filename), '>>'

      nlat = hist_block_info%ginfo%nlat
      nlon = hist_block_info%ginfo%nlon

      allocate (vdata (ndim1, ndim2, nlon, nlat, timelen))
      vdata(:,:,:,:,:) = spval

      do iyseg = 1, hist_block_info%nyseg
         do ixseg = 1, hist_block_info%nxseg

            call get_filename_block (filename, &
               hist_block_info%xsegs(ixseg)%blk, hist_block_info%ysegs(iyseg)%blk, fileblock)

            inquire(file=fileblock, exist=fexists)
            IF (fexists) THEN

               CALL ncio_read_serial (fileblock, varname, rcache)

               xgdsp = hist_block_info%xsegs(ixseg)%gdsp
               ygdsp = hist_block_info%ysegs(iyseg)%gdsp
               xcnt  = hist_block_info%xsegs(ixseg)%cnt
               ycnt  = hist_block_info%ysegs(iyseg)%cnt

               vdata (:,:, xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt,:) = rcache
            end if

         end do
      end do

      call ncio_write_serial (filename, varname, vdata, dim1name, dim2name, 'lon', 'lat', 'time', &
            compress)

      IF (present(longname)) THEN
         CALL ncio_put_attr (filename, varname, 'long_name', longname)
      ENDIF
      IF (present(units)) THEN
         CALL ncio_put_attr (filename, varname, 'units', units)
      ENDIF

      CALL ncio_put_attr (filename, varname, 'missing_value', spval)

      IF (allocated(rcache)) deallocate(rcache)
      IF (allocated(vdata )) deallocate(vdata )

   end subroutine hist_concatenate_var_4d
#endif

#ifdef USEMPI
   ! -----
   subroutine hist_concatenate_var_2d (filename, varname, timelen, compress, longname, units)

      implicit none

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: varname
      integer, intent(in) :: timelen
      integer, intent(in) :: compress
      character (len=*), intent(in), optional :: longname
      character (len=*), intent(in), optional :: units

      ! Local variables
      CHARACTER(len=256) :: fileblock
      integer :: ixseg, iyseg, iblock, jblock, nblock
      integer :: nlon, nlat, xgdsp, ygdsp, xcnt, ycnt
      real(r8), allocatable :: vdata(:,:,:), rcache(:,:,:)
      INTEGER :: rmesg(2), smesg(2), idest, isrc
      INTEGER :: nrecv

      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_master) THEN

         write(*,*) 'Concatenate <', trim(varname), '> to <<', trim(filename), '>>'

         nlat = hist_block_info%ginfo%nlat
         nlon = hist_block_info%ginfo%nlon

         allocate (vdata (nlon, nlat, timelen))
         vdata(:,:,:) = spval

         nblock = hist_block_info%nxseg * hist_block_info%nyseg

         iblock = 1
         nrecv  = 0
         DO WHILE (nrecv < nblock)

            CALL mpi_recv (rmesg, 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            IF (rmesg(2) > 0) THEN
               jblock = rmesg(2)
               ixseg = mod(jblock, hist_block_info%nxseg)
               IF (ixseg == 0) ixseg = hist_block_info%nxseg
               iyseg = (jblock-1) / hist_block_info%nxseg + 1

               xgdsp = hist_block_info%xsegs(ixseg)%gdsp
               ygdsp = hist_block_info%ysegs(iyseg)%gdsp
               xcnt  = hist_block_info%xsegs(ixseg)%cnt
               ycnt  = hist_block_info%ysegs(iyseg)%cnt

               IF (allocated(rcache)) deallocate(rcache)
               allocate (rcache (xcnt, ycnt, timelen))

               isrc = rmesg(1)
               CALL mpi_recv (rcache, size(rcache), MPI_DOUBLE, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               vdata (xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt,:) = rcache

               nrecv = nrecv + 1
            ENDIF

            idest  = rmesg(1)
            IF (iblock <= nblock) THEN
               ixseg = mod(iblock, hist_block_info%nxseg)
               IF (ixseg == 0) ixseg = hist_block_info%nxseg
               iyseg = (iblock-1) / hist_block_info%nxseg + 1

               call get_filename_block (filename, &
                  hist_block_info%xsegs(ixseg)%blk, hist_block_info%ysegs(iyseg)%blk, fileblock)

               CALL mpi_send (iblock,      1,   MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)
               CALL mpi_send (fileblock, 256, MPI_CHARACTER, idest, mpi_tag_mesg, p_comm_glb, p_err)

               iblock = iblock + 1
            ELSE
               jblock = 0
               CALL mpi_send (jblock, 1, MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)
            ENDIF

         ENDDO

         call ncio_write_serial (filename, varname, vdata, 'lon', 'lat', 'time', compress)

         IF (present(longname)) THEN
            CALL ncio_put_attr (filename, varname, 'long_name', longname)
         ENDIF
         IF (present(units)) THEN
            CALL ncio_put_attr (filename, varname, 'units', units)
         ENDIF

         CALL ncio_put_attr (filename, varname, 'missing_value', spval)

      ENDIF

      IF (.not. p_is_master) THEN
         smesg(:) = (/p_iam_glb, 0/)
         CALL mpi_send (smesg, 2, MPI_INTEGER, p_root, mpi_tag_mesg, p_comm_glb, p_err)

         DO WHILE (.true.)

            CALL mpi_recv (iblock, 1, MPI_INTEGER, p_root, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            IF (iblock /= 0) THEN
               CALL mpi_recv (fileblock, 256, MPI_CHARACTER, &
                  p_root, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

               CALL ncio_read_serial (fileblock, varname, rcache)

               smesg(:) = (/p_iam_glb, iblock/)
               CALL mpi_send (smesg, 2, MPI_INTEGER, &
                  p_root, mpi_tag_mesg, p_comm_glb, p_err)

               CALL mpi_send (rcache, size(rcache), MPI_DOUBLE, &
                  p_root, mpi_tag_data, p_comm_glb, p_err)
            ELSE
               exit
            ENDIF
         ENDDO
      ENDIF

      IF (allocated(rcache)) deallocate(rcache)
      IF (allocated(vdata )) deallocate(vdata )

      CALL mpi_barrier (p_comm_glb, p_err)

   end subroutine hist_concatenate_var_2d

   ! -----
   subroutine hist_concatenate_var_3d (filename, varname, timelen, dim1name, ndim1, compress, &
         longname, units)

      implicit none

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: varname
      INTEGER, intent(in) :: timelen
      character(len=*), intent(in) :: dim1name
      INTEGER, intent(in) :: ndim1
      integer, intent(in) :: compress
      character(len=*), intent(in), optional :: longname
      character(len=*), intent(in), optional :: units

      ! Local variables
      CHARACTER(len=256) :: fileblock
      integer :: nlon, nlat, xgdsp, ygdsp, xcnt, ycnt
      integer :: ixseg, iyseg, iblock, jblock, nblock
      real(r8), allocatable :: vdata(:,:,:,:), rcache(:,:,:,:)
      INTEGER :: rmesg(2), smesg(2), idest, isrc
      INTEGER :: nrecv

      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_master) THEN

         write(*,*) 'Concatenate <', trim(varname), '> to <<', trim(filename), '>>'

         nlat = hist_block_info%ginfo%nlat
         nlon = hist_block_info%ginfo%nlon

         allocate (vdata (ndim1, nlon, nlat, timelen))
         vdata(:,:,:,:) = spval

         nblock = hist_block_info%nxseg * hist_block_info%nyseg

         iblock = 1
         nrecv  = 0
         DO WHILE (nrecv < nblock)

            CALL mpi_recv (rmesg, 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            IF (rmesg(2) > 0) THEN
               jblock = rmesg(2)
               ixseg = mod(jblock, hist_block_info%nxseg)
               IF (ixseg == 0) ixseg = hist_block_info%nxseg
               iyseg = (jblock-1) / hist_block_info%nxseg + 1

               xgdsp = hist_block_info%xsegs(ixseg)%gdsp
               ygdsp = hist_block_info%ysegs(iyseg)%gdsp
               xcnt  = hist_block_info%xsegs(ixseg)%cnt
               ycnt  = hist_block_info%ysegs(iyseg)%cnt

               IF (allocated(rcache)) deallocate(rcache)
               allocate (rcache (ndim1, xcnt, ycnt, timelen))

               isrc = rmesg(1)
               CALL mpi_recv (rcache, size(rcache), MPI_DOUBLE, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               vdata (:, xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt,:) = rcache

               nrecv = nrecv + 1
            ENDIF

            idest  = rmesg(1)
            IF (iblock <= nblock) THEN
               ixseg = mod(iblock, hist_block_info%nxseg)
               IF (ixseg == 0) ixseg = hist_block_info%nxseg
               iyseg = (iblock-1) / hist_block_info%nxseg + 1
               call get_filename_block (filename, &
                  hist_block_info%xsegs(ixseg)%blk, hist_block_info%ysegs(iyseg)%blk, fileblock)

               CALL mpi_send (iblock,      1,  MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)
               CALL mpi_send (fileblock, 256, MPI_CHARACTER,idest, mpi_tag_mesg, p_comm_glb, p_err)

               iblock = iblock + 1
            ELSE
               jblock = 0
               CALL mpi_send (jblock, 1, MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)
            ENDIF

         ENDDO

         call ncio_write_serial (filename, varname, vdata, dim1name, &
            'lon', 'lat', 'time', compress)

         IF (present(longname)) THEN
            CALL ncio_put_attr (filename, varname, 'long_name', longname)
         ENDIF
         IF (present(units)) THEN
            CALL ncio_put_attr (filename, varname, 'units', units)
         ENDIF

         CALL ncio_put_attr (filename, varname, 'missing_value', spval)

      ENDIF

      IF (.not. p_is_master) THEN
         smesg(:) = (/p_iam_glb, 0/)
         CALL mpi_send (smesg, 2, MPI_INTEGER, &
            p_root, mpi_tag_mesg, p_comm_glb, p_err)

         DO WHILE (.true.)

            CALL mpi_recv (iblock, 1, MPI_INTEGER, &
               p_root, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            IF (iblock /= 0) THEN
               CALL mpi_recv (fileblock, 256, MPI_CHARACTER, &
                  p_root, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

               CALL ncio_read_serial (fileblock, varname, rcache)

               smesg(:) = (/p_iam_glb, iblock/)
               CALL mpi_send (smesg, 2, MPI_INTEGER, &
                  p_root, mpi_tag_mesg, p_comm_glb, p_err)

               CALL mpi_send (rcache, size(rcache), MPI_DOUBLE, &
                  p_root, mpi_tag_data, p_comm_glb, p_err)
            ELSE
               exit
            ENDIF
         ENDDO
      ENDIF

      IF (allocated(rcache)) deallocate(rcache)
      IF (allocated(vdata )) deallocate(vdata )

      CALL mpi_barrier (p_comm_glb, p_err)

   end subroutine hist_concatenate_var_3d

   ! -----
   subroutine hist_concatenate_var_4d (filename, varname, timelen, &
         dim1name, dim2name, ndim1, ndim2, compress, longname, units)

      implicit none

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: varname
      INTEGER, intent(in) :: timelen
      character(len=*), intent(in) :: dim1name, dim2name
      INTEGER, intent(in) :: ndim1, ndim2
      integer, intent(in) :: compress
      character (len=*), intent(in), optional :: longname
      character (len=*), intent(in), optional :: units

      ! Local variables
      CHARACTER(len=256) :: fileblock
      integer :: ixseg, iyseg, iblock, jblock, nblock
      integer :: nlon, nlat, xgdsp, ygdsp, xcnt, ycnt
      real(r8), allocatable :: vdata(:,:,:,:,:), rcache(:,:,:,:,:)
      INTEGER :: rmesg(2), smesg(2), idest, isrc
      INTEGER :: nrecv

      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_master) THEN

         write(*,*) 'Concatenate <', trim(varname), '> to <<', trim(filename), '>>'

         nlat = hist_block_info%ginfo%nlat
         nlon = hist_block_info%ginfo%nlon

         allocate (vdata (ndim1, ndim2, nlon, nlat, timelen))
         vdata(:,:,:,:,:) = spval

         nblock = hist_block_info%nxseg * hist_block_info%nyseg

         iblock = 1
         nrecv  = 0
         DO WHILE (nrecv < nblock)

            CALL mpi_recv (rmesg, 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            IF (rmesg(2) > 0) THEN
               jblock = rmesg(2)
               ixseg = mod(jblock, hist_block_info%nxseg)
               IF (ixseg == 0) ixseg = hist_block_info%nxseg
               iyseg = (jblock-1) / hist_block_info%nxseg + 1

               xgdsp = hist_block_info%xsegs(ixseg)%gdsp
               ygdsp = hist_block_info%ysegs(iyseg)%gdsp
               xcnt  = hist_block_info%xsegs(ixseg)%cnt
               ycnt  = hist_block_info%ysegs(iyseg)%cnt

               IF (allocated(rcache)) deallocate(rcache)
               allocate (rcache (ndim1, ndim2, xcnt, ycnt, timelen))

               isrc = rmesg(1)
               CALL mpi_recv (rcache, size(rcache), MPI_DOUBLE, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               vdata (:, :, xgdsp+1:xgdsp+xcnt, ygdsp+1:ygdsp+ycnt,:) = rcache

               nrecv = nrecv + 1
            ENDIF

            idest  = rmesg(1)
            IF (iblock <= nblock) THEN
               ixseg = mod(iblock, hist_block_info%nxseg)
               IF (ixseg == 0) ixseg = hist_block_info%nxseg
               iyseg = (iblock-1) / hist_block_info%nxseg + 1
               call get_filename_block (filename, &
                  hist_block_info%xsegs(ixseg)%blk, hist_block_info%ysegs(iyseg)%blk, fileblock)

               CALL mpi_send (iblock,      1,   MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)
               CALL mpi_send (fileblock, 256, MPI_CHARACTER, idest, mpi_tag_mesg, p_comm_glb, p_err)

               iblock = iblock + 1
            ELSE
               jblock = 0
               CALL mpi_send (jblock, 1, MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)
            ENDIF

         ENDDO

         call ncio_write_serial (filename, varname, vdata, dim1name, dim2name, &
            'lon', 'lat', 'time', compress)

         IF (present(longname)) THEN
            CALL ncio_put_attr (filename, varname, 'long_name', longname)
         ENDIF
         IF (present(units)) THEN
            CALL ncio_put_attr (filename, varname, 'units', units)
         ENDIF

         CALL ncio_put_attr (filename, varname, 'missing_value', spval)

      ENDIF

      IF (.not. p_is_master) THEN
         smesg(:) = (/p_iam_glb, 0/)
         CALL mpi_send (smesg, 2, MPI_INTEGER, &
            p_root, mpi_tag_mesg, p_comm_glb, p_err)

         DO WHILE (.true.)

            CALL mpi_recv (iblock, 1, MPI_INTEGER, &
               p_root, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            IF (iblock /= 0) THEN
               CALL mpi_recv (fileblock, 256, MPI_CHARACTER, &
                  p_root, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

               CALL ncio_read_serial (fileblock, varname, rcache)

               smesg(:) = (/p_iam_glb, iblock/)
               CALL mpi_send (smesg, 2, MPI_INTEGER, &
                  p_root, mpi_tag_mesg, p_comm_glb, p_err)

               CALL mpi_send (rcache, size(rcache), MPI_DOUBLE, &
                  p_root, mpi_tag_data, p_comm_glb, p_err)
            ELSE
               exit
            ENDIF
         ENDDO
      ENDIF

      IF (allocated(rcache)) deallocate(rcache)
      IF (allocated(vdata )) deallocate(vdata )

      CALL mpi_barrier (p_comm_glb, p_err)

   end subroutine hist_concatenate_var_4d
#endif

END module mod_concatenate
