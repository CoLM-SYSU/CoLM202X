#include <define.h>

module mod_hist_nc

   use precision
   USE spmd_task
   USE mod_block
   USE ncio_serial
   use netcdf
   implicit none

   type :: segment_type
      integer :: blk
      integer :: cnt
      integer :: bdsp
      integer :: gdsp
   end type segment_type

   integer :: nxseg
   type(segment_type), allocatable :: xsegment(:)

   integer :: nyseg
   type(segment_type), allocatable :: ysegment(:)

   integer :: nlat, nlon
   real(r8), allocatable :: lat_s(:)
   real(r8), allocatable :: lat_n(:)
   real(r8), allocatable :: lon_w(:)
   real(r8), allocatable :: lon_e(:)

   ! public subroutines

   public set_segment_info
   
   public :: hist_concatenate_var_2d
   public :: hist_concatenate_var_3d
   public :: hist_concatenate_var_4d

   public :: hist_concatenate_time
   public :: hist_concatenate_grid_info

contains

   ! -----
   subroutine set_segment_info (ghist)

      use mod_grid
      implicit none

      type(grid_type), intent(in) :: ghist

      ! Local variables
      integer :: ilat_l, ilat_u, ilat, ilatloc, jblk, iyseg
      integer :: ilon_w, ilon_e, ilon, ilonloc, iblk, ixseg

      ilat_l = findloc(ghist%yblk /= 0, .true., dim=1)
      ilat_u = findloc(ghist%yblk /= 0, .true., dim=1, back=.true.)

      nlat = ilat_u - ilat_l + 1
      allocate (lat_s(nlat))
      allocate (lat_n(nlat))
      
      nyseg = 0
      jblk  = 0
      ilatloc = 0
      do ilat = ilat_l, ilat_u
         if (ghist%yblk(ilat) /= jblk) then
            nyseg = nyseg + 1
            jblk  = ghist%yblk(ilat)
         end if

         ilatloc = ilatloc + 1
         lat_s(ilatloc) = ghist%lat_s(ilat)
         lat_n(ilatloc) = ghist%lat_n(ilat)
      end do

      allocate (ysegment (nyseg))
      
      iyseg = 0
      jblk  = 0
      do ilat = ilat_l, ilat_u
         if (ghist%yblk(ilat) /= jblk) then
            iyseg = iyseg + 1
            jblk  = ghist%yblk(ilat)
            ysegment(iyseg)%blk  = jblk
            ysegment(iyseg)%bdsp = ghist%yloc(ilat) - 1
            ysegment(iyseg)%gdsp = ilat - ilat_l
            ysegment(iyseg)%cnt  = 1
         else
            ysegment(iyseg)%cnt  = ysegment(iyseg)%cnt + 1
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

      nlon = ilon_e - ilon_w + 1
      if (nlon <= 0) nlon = nlon + ghist%nlon

      allocate (lon_w(nlon))
      allocate (lon_e(nlon))

      nxseg = 0
      ilon = ilon_w - 1
      iblk = 0
      ilonloc = 0
      do while (.true.) 
         ilon = mod(ilon,ghist%nlon) + 1
         if (ghist%xblk(ilon) /= iblk) then
            nxseg = nxseg + 1
            iblk = ghist%xblk(ilon)
         end if

         ilonloc = ilonloc + 1
         lon_w(ilonloc) = ghist%lon_w(ilon)
         lon_e(ilonloc) = ghist%lon_e(ilon)

         if (ilon == ilon_e) exit
      end do

      allocate (xsegment (nxseg))

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
            xsegment(ixseg)%blk  = iblk
            xsegment(ixseg)%bdsp = ghist%xloc(ilon) - 1
            xsegment(ixseg)%gdsp = ilonloc - 1
            xsegment(ixseg)%cnt = 1
         else
            xsegment(ixseg)%cnt = xsegment(ixseg)%cnt + 1
         end if

         if (ilon == ilon_e) exit
      end do

   end subroutine set_segment_info

   !------------------------------
   subroutine hist_concatenate_time (filename, dataname, timelen)

      implicit none

      character (len=*), intent(in) :: filename
      character (len=*), intent(in) :: dataname
      INTEGER, intent(out) :: timelen

      ! Local variables
      CHARACTER(len=256) :: fileblock
      integer :: ixseg, iyseg
      INTEGER, allocatable :: time_component(:,:)

      iyseg = 1
      DO WHILE (ysegment(iyseg)%cnt <= 0)
         iyseg = iyseg + 1
      ENDDO

      ixseg = 1
      DO WHILE (xsegment(ixseg)%cnt <= 0)
         ixseg = ixseg + 1
      ENDDO
               
      call get_filename_block ( &
         filename, xsegment(ixseg)%blk, ysegment(iyseg)%blk, fileblock)

      CALL ncio_read_serial (fileblock, dataname, time_component)
      
      timelen = size(time_component,2)
      CALL ncio_define_dimension (filename, 'time',      timelen)
      CALL ncio_define_dimension (filename, 'component', 3)

      call ncio_write_serial (filename, dataname, time_component, 'component', 'time')

   end subroutine hist_concatenate_time

   !------------------
   subroutine hist_concatenate_grid_info (filename)
      
      use ncio_serial
      implicit none

      character(len=*), intent(in) :: filename

      CALL ncio_define_dimension (filename, 'lat', nlat)
      CALL ncio_define_dimension (filename, 'lon', nlon)

      call ncio_write_serial (filename, 'lat_s', lat_s, 'lat')
      call ncio_write_serial (filename, 'lat_n', lat_n, 'lat')
      call ncio_write_serial (filename, 'lon_w', lon_w, 'lon')
      call ncio_write_serial (filename, 'lon_e', lon_e, 'lon')

   end subroutine hist_concatenate_grid_info

#ifndef USEMPI
   ! -----
   subroutine hist_concatenate_var_2d (filename, varname, timelen, compress)

      implicit none

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: varname
      integer, intent(in) :: timelen 
      integer, intent(in), optional :: compress

      ! Local variables
      CHARACTER(len=256) :: fileblock
      integer :: ixseg, iyseg
      real(r8), allocatable :: vdata(:,:,:), rcache(:,:,:)

      write(*,*) 'Concatenate ', trim(varname), ' to ', trim(filename)

      allocate (vdata (nlon, nlat, timelen))

      do iyseg = 1, nyseg
         do ixseg = 1, nxseg

            if ((xsegment(ixseg)%cnt > 0) .and. (ysegment(iyseg)%cnt > 0)) then

               call get_filename_block (filename, &
                  xsegment(ixseg)%blk, ysegment(iyseg)%blk, fileblock)
               CALL ncio_read_serial (fileblock, varname, rcache)

               vdata (xsegment(ixseg)%gdsp+1:xsegment(ixseg)%gdsp+xsegment(ixseg)%cnt, &
                  ysegment(iyseg)%gdsp+1:ysegment(iyseg)%gdsp+ysegment(iyseg)%cnt, :)  &
                  = rcache
            end if

         end do
      end do

      if (present(compress)) then
         call ncio_write_serial (filename, varname, vdata, 'lon', 'lat', 'time', compress)
      else
         call ncio_write_serial (filename, varname, vdata, 'lon', 'lat', 'time')
      end if

      IF (allocated(rcache)) deallocate(rcache)
      IF (allocated(vdata )) deallocate(vdata )
   
   end subroutine hist_concatenate_var_2d

   ! -----
   subroutine hist_concatenate_var_3d (filename, varname, timelen, dim1name, ndim1, compress)

      implicit none

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: varname
      INTEGER, intent(in) :: timelen
      character(len=*), intent(in) :: dim1name 
      INTEGER, intent(in) :: ndim1
      integer, intent(in), optional :: compress

      ! Local variables
      CHARACTER(len=256) :: fileblock
      integer :: ixseg, iyseg
      real(r8), allocatable :: vdata(:,:,:,:), rcache(:,:,:,:)

      allocate (vdata (ndim1, nlon, nlat, timelen))
         
      write(*,*) 'Concatenate ', trim(varname), ' to ', trim(filename)

      do iyseg = 1, nyseg
         do ixseg = 1, nxseg

            if ((xsegment(ixseg)%cnt > 0) .and. (ysegment(iyseg)%cnt > 0)) then
               call get_filename_block (filename, &
                  xsegment(ixseg)%blk, ysegment(iyseg)%blk, fileblock)
               CALL ncio_read_serial (fileblock, varname, rcache)
               
               vdata (:, &
                  xsegment(ixseg)%gdsp+1:xsegment(ixseg)%gdsp+xsegment(ixseg)%cnt, &
                  ysegment(iyseg)%gdsp+1:ysegment(iyseg)%gdsp+ysegment(iyseg)%cnt, :) &
                  = rcache
            end if

         end do
      end do
      
      if (present(compress)) then
         call ncio_write_serial (filename, varname, vdata, &
            dim1name, 'lon', 'lat', 'time', compress)
      else
         call ncio_write_serial (filename, varname, vdata, &
            dim1name, 'lon', 'lat', 'time')
      end if
      
      IF (allocated(rcache)) deallocate(rcache)
      IF (allocated(vdata )) deallocate(vdata )
   
   end subroutine hist_concatenate_var_3d

   ! -----
   subroutine hist_concatenate_var_4d ( &
         filename, varname, timelen, dim1name, dim2name, ndim1, ndim2, compress)

      implicit none

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: varname
      INTEGER, intent(in) :: timelen 
      character(len=*), intent(in) :: dim1name, dim2name 
      INTEGER, intent(in) :: ndim1, ndim2
      integer, intent(in), optional :: compress

      ! Local variables
      CHARACTER(len=256) :: fileblock
      integer :: ixseg, iyseg
      real(r8), allocatable :: vdata(:,:,:,:,:), rcache(:,:,:,:,:)

      write(*,*) 'Concatenate ', trim(varname), ' to ', trim(filename)

      allocate (vdata (ndim1, ndim2, nlon, nlat, timelen))

      do iyseg = 1, nyseg
         do ixseg = 1, nxseg

            if ((xsegment(ixseg)%cnt > 0) .and. (ysegment(iyseg)%cnt > 0)) then
               call get_filename_block (filename, &
                  xsegment(ixseg)%blk, ysegment(iyseg)%blk, fileblock)
               CALL ncio_read_serial (fileblock, varname, rcache)
               
               vdata (:, :, &
                  xsegment(ixseg)%gdsp+1:xsegment(ixseg)%gdsp+xsegment(ixseg)%cnt, &
                  ysegment(iyseg)%gdsp+1:ysegment(iyseg)%gdsp+ysegment(iyseg)%cnt, :) &
                  = rcache
            end if

         end do
      end do
      
      if (present(compress)) then
         call ncio_write_serial (filename, varname, vdata, &
            dim1name, dim2name, 'lon', 'lat', 'time', compress)
      else
         call ncio_write_serial (filename, varname, vdata, &
            dim1name, dim2name, 'lon', 'lat', 'time')
      end if
   
      IF (allocated(rcache)) deallocate(rcache)
      IF (allocated(vdata )) deallocate(vdata )

   end subroutine hist_concatenate_var_4d
#endif

#ifdef USEMPI
   ! -----
   subroutine hist_concatenate_var_2d (filename, varname, timelen, compress)

      implicit none

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: varname
      integer, intent(in) :: timelen 
      integer, intent(in), optional :: compress

      ! Local variables
      CHARACTER(len=256) :: fileblock
      integer :: ixseg, iyseg, iblock, jblock, nblock
      real(r8), allocatable :: vdata(:,:,:), rcache(:,:,:)
      INTEGER :: rmesg(2), smesg(2), idest, isrc
      INTEGER :: nrecv

      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_master) THEN

         write(*,102) trim(varname), filename
         102 format('Concatenate ', A20, ' to ', A100)

         allocate (vdata (nlon, nlat, timelen))

         iblock = 1
         nblock = nxseg * nyseg
         nrecv  = 0
         DO WHILE (nrecv < nblock)

            CALL mpi_recv (rmesg, 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            IF (rmesg(2) > 0) THEN
               jblock = rmesg(2)
               ixseg = mod(jblock, nxseg)
               IF (ixseg == 0) ixseg = nxseg
               iyseg = (jblock-1) / nxseg + 1
                  
               IF (allocated(rcache)) deallocate(rcache)
               allocate (rcache (xsegment(ixseg)%cnt, ysegment(iyseg)%cnt, timelen))

               isrc = rmesg(1)
               CALL mpi_recv (rcache, size(rcache), MPI_DOUBLE, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               vdata (xsegment(ixseg)%gdsp+1:xsegment(ixseg)%gdsp+xsegment(ixseg)%cnt, &
                  ysegment(iyseg)%gdsp+1:ysegment(iyseg)%gdsp+ysegment(iyseg)%cnt, :)  &
                  = rcache
               
               nrecv = nrecv + 1
            ENDIF

            idest  = rmesg(1)
            IF (iblock <= nblock) THEN
               ixseg = mod(iblock, nxseg)
               IF (ixseg == 0) ixseg = nxseg
               iyseg = (iblock-1) / nxseg + 1
               call get_filename_block (filename, &
                  xsegment(ixseg)%blk, ysegment(iyseg)%blk, fileblock)

               CALL mpi_send (iblock, 1, MPI_INTEGER, &
                  idest, mpi_tag_mesg, p_comm_glb, p_err) 
               CALL mpi_send (fileblock, 256, MPI_CHARACTER, &
                  idest, mpi_tag_mesg, p_comm_glb, p_err) 

               iblock = iblock + 1
            ELSE
               jblock = 0
               CALL mpi_send (jblock, 1, MPI_INTEGER, &
                  idest, mpi_tag_mesg, p_comm_glb, p_err) 
            ENDIF

         ENDDO
      
         if (present(compress)) then
            call ncio_write_serial (filename, varname, vdata, 'lon', 'lat', 'time', compress)
         else
            call ncio_write_serial (filename, varname, vdata, 'lon', 'lat', 'time')
         end if

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
   
   end subroutine hist_concatenate_var_2d

   ! -----
   subroutine hist_concatenate_var_3d (filename, varname, timelen, dim1name, ndim1, compress)

      implicit none

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: varname
      INTEGER, intent(in) :: timelen
      character(len=*), intent(in) :: dim1name 
      INTEGER, intent(in) :: ndim1
      integer, intent(in), optional :: compress

      ! Local variables
      CHARACTER(len=256) :: fileblock
      integer :: ixseg, iyseg, iblock, jblock, nblock
      real(r8), allocatable :: vdata(:,:,:,:), rcache(:,:,:,:)
      INTEGER :: rmesg(2), smesg(2), idest, isrc
      INTEGER :: nrecv

      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_master) THEN

         write(*,103) trim(varname), filename
         103 format('Concatenate ', A20, ' to ', A100)

         allocate (vdata (ndim1, nlon, nlat, timelen))

         iblock = 1
         nblock = nxseg * nyseg
         nrecv  = 0
         DO WHILE (nrecv < nblock)

            CALL mpi_recv (rmesg, 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            IF (rmesg(2) > 0) THEN
               jblock = rmesg(2)
               ixseg = mod(jblock, nxseg)
               IF (ixseg == 0) ixseg = nxseg
               iyseg = (jblock-1) / nxseg + 1
                  
               IF (allocated(rcache)) deallocate(rcache)
               allocate (rcache (ndim1, xsegment(ixseg)%cnt, ysegment(iyseg)%cnt, timelen))

               isrc = rmesg(1)
               CALL mpi_recv (rcache, size(rcache), MPI_DOUBLE, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               vdata (:, xsegment(ixseg)%gdsp+1:xsegment(ixseg)%gdsp+xsegment(ixseg)%cnt, &
                  ysegment(iyseg)%gdsp+1:ysegment(iyseg)%gdsp+ysegment(iyseg)%cnt, :)  &
                  = rcache
               
               nrecv = nrecv + 1
            ENDIF

            idest  = rmesg(1)
            IF (iblock <= nblock) THEN
               ixseg = mod(iblock, nxseg)
               IF (ixseg == 0) ixseg = nxseg
               iyseg = (iblock-1) / nxseg + 1
               call get_filename_block (filename, &
                  xsegment(ixseg)%blk, ysegment(iyseg)%blk, fileblock)

               CALL mpi_send (iblock, 1, MPI_INTEGER, &
                  idest, mpi_tag_mesg, p_comm_glb, p_err) 
               CALL mpi_send (fileblock, 256, MPI_CHARACTER, &
                  idest, mpi_tag_mesg, p_comm_glb, p_err) 

               iblock = iblock + 1
            ELSE
               jblock = 0
               CALL mpi_send (jblock, 1, MPI_INTEGER, &
                  idest, mpi_tag_mesg, p_comm_glb, p_err) 
            ENDIF

         ENDDO
      
         if (present(compress)) then
            call ncio_write_serial (filename, varname, vdata, &
               dim1name, 'lon', 'lat', 'time', compress)
         else
            call ncio_write_serial (filename, varname, vdata, &
               dim1name, 'lon', 'lat', 'time')
         end if

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
         dim1name, dim2name, ndim1, ndim2, compress)

      implicit none

      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: varname
      INTEGER, intent(in) :: timelen
      character(len=*), intent(in) :: dim1name, dim2name
      INTEGER, intent(in) :: ndim1, ndim2
      integer, intent(in), optional :: compress

      ! Local variables
      CHARACTER(len=256) :: fileblock
      integer :: ixseg, iyseg, iblock, jblock, nblock
      real(r8), allocatable :: vdata(:,:,:,:,:), rcache(:,:,:,:,:)
      INTEGER :: rmesg(2), smesg(2), idest, isrc
      INTEGER :: nrecv

      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_master) THEN

         write(*,104) trim(varname), filename
         104 format('Concatenate ', A20, ' to ', A100)

         allocate (vdata (ndim1, ndim2, nlon, nlat, timelen))

         iblock = 1
         nblock = nxseg * nyseg
         nrecv  = 0
         DO WHILE (nrecv < nblock)

            CALL mpi_recv (rmesg, 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            IF (rmesg(2) > 0) THEN
               jblock = rmesg(2)
               ixseg = mod(jblock, nxseg)
               IF (ixseg == 0) ixseg = nxseg
               iyseg = (jblock-1) / nxseg + 1
                  
               IF (allocated(rcache)) deallocate(rcache)
               allocate (rcache (ndim1, ndim2, xsegment(ixseg)%cnt, ysegment(iyseg)%cnt, timelen))

               isrc = rmesg(1)
               CALL mpi_recv (rcache, size(rcache), MPI_DOUBLE, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               vdata (:, :, &
                  xsegment(ixseg)%gdsp+1:xsegment(ixseg)%gdsp+xsegment(ixseg)%cnt, &
                  ysegment(iyseg)%gdsp+1:ysegment(iyseg)%gdsp+ysegment(iyseg)%cnt, :)  &
                  = rcache
               
               nrecv = nrecv + 1
            ENDIF

            idest  = rmesg(1)
            IF (iblock <= nblock) THEN
               ixseg = mod(iblock, nxseg)
               IF (ixseg == 0) ixseg = nxseg
               iyseg = (iblock-1) / nxseg + 1
               call get_filename_block (filename, &
                  xsegment(ixseg)%blk, ysegment(iyseg)%blk, fileblock)

               CALL mpi_send (iblock, 1, MPI_INTEGER, &
                  idest, mpi_tag_mesg, p_comm_glb, p_err) 
               CALL mpi_send (fileblock, 256, MPI_CHARACTER, &
                  idest, mpi_tag_mesg, p_comm_glb, p_err) 

               iblock = iblock + 1
            ELSE
               jblock = 0
               CALL mpi_send (jblock, 1, MPI_INTEGER, &
                  idest, mpi_tag_mesg, p_comm_glb, p_err) 
            ENDIF

         ENDDO
      
         if (present(compress)) then
            call ncio_write_serial (filename, varname, vdata, &
               dim1name, dim2name, 'lon', 'lat', 'time', compress)
         else
            call ncio_write_serial (filename, varname, vdata, &
               dim1name, dim2name, 'lon', 'lat', 'time')
         end if

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

end module mod_hist_nc
