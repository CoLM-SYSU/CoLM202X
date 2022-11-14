#include <define.h>

module mod_modis_data

   implicit none

   integer, parameter :: nxbox = 1200
   integer, parameter :: nybox = 1200
   integer, parameter :: nxglb = 86400
   integer, parameter :: nyglb = 43200
   
   INTEGER, parameter :: N_PFT_modis = 16
   
   interface modis_read_data
      MODULE procedure modis_read_data_int32
      MODULE procedure modis_read_data_real8
   END interface modis_read_data

   PUBLIC :: modis_read_data_pft
   PUBLIC :: modis_read_data_time
   PUBLIC :: modis_read_data_pft_time

contains

   ! ----
   SUBROUTINE nccheck (status)
      USE netcdf

      INTEGER, INTENT(IN) :: status

      IF (status /= NF90_NOERR) THEN
         print *, trim(nf90_strerror(status))
         stop 2
      ENDIF
   END SUBROUTINE nccheck

   ! -----
   subroutine this_block_and_move_to_next (dir_modis, grid, &
         isouth, inorth, iwest, ieast, ibox, jbox, ibox0, &
         i0, i1, j0, j1, il0, il1, jl0, jl1, &
         file_modis)

      use mod_grid
      implicit none

      character (len=*), intent(in) :: dir_modis
      type (grid_type),  intent(in) :: grid
      
      integer, intent(in)    :: isouth, inorth, iwest, ieast
      INTEGER, intent(inout) :: ibox, jbox, ibox0
      integer, intent(out)   :: i0, i1, j0, j1, il0, il1, jl0, jl1

      character (len=*), intent(out) :: file_modis

      ! Local variables
      integer :: xdsp, ydsp
      character(len=4) :: str


      xdsp = (ibox-1) * nxbox
      ydsp = (jbox-1) * nybox

      j0 = max(inorth-ydsp, 1)
      j1 = min(isouth-ydsp, nybox)
      jl0 = j0 + ydsp - inorth + 1
      jl1 = j1 + ydsp - inorth + 1

      if (ieast >= iwest) then
         i0 = max(iwest-xdsp, 1)
         i1 = min(ieast-xdsp, nxbox)
         il0 = i0 + xdsp - iwest + 1
         il1 = i1 + xdsp - iwest + 1
      else
         if (iwest <= xdsp+nxbox) then
            i0 = max(iwest-xdsp, 1)
            i1 = nxbox
            il0 = i0 + xdsp - iwest + 1
            il1 = i1 + xdsp - iwest + 1
         else
            i0 = 1
            i1 = min(ieast-xdsp, nxbox)
            il0 = i0 + xdsp + nxglb - iwest + 1
            il1 = i1 + xdsp + nxglb - iwest + 1
         end if
      end if

      file_modis = trim(dir_modis) // '/RG'
      write(str, '(I4)') (19-jbox)*5
      file_modis = trim(file_modis) // '_' // trim(adjustl(str))
      write(str, '(I4)') (ibox-37)*5
      file_modis = trim(file_modis) // '_' // trim(adjustl(str))
      write(str, '(I4)') (18-jbox)*5
      file_modis = trim(file_modis) // '_' // trim(adjustl(str))
      write(str, '(I4)') (ibox-36)*5
      file_modis = trim(file_modis) // '_' // trim(adjustl(str))
      file_modis = trim(file_modis) // '.MOD2005.nc'

      if ((ieast >= xdsp + 1) .and. (ieast <= xdsp + nxbox)) then 
         if (isouth <= ydsp + nybox) then 
            jbox = -1
         else
            ibox = ibox0
            jbox = jbox + 1
         end if
      else
         ibox = mod(ibox, nxglb/nxbox) + 1
      end if

   end subroutine this_block_and_move_to_next 

   ! -----
   subroutine modis_read_data_int32 (dir_modis, dataname, grid, rdata)

      use spmd_task
      use mod_block
      use mod_grid
      use mod_data_type
      use netcdf
      implicit none

      character (len=*), intent(in) :: dir_modis
      character (len=*), intent(in) :: dataname
      type (grid_type),  intent(in) :: grid
      type (block_data_int32_2d), intent(inout) :: rdata

      ! Local variables
      integer :: iblkme, iblk, jblk, isouth, inorth, iwest, ieast, ibox, jbox, ibox0
      integer :: i0, i1, j0, j1, il0, il1, jl0, jl1
      character(len=256) :: file_modis
      INTEGER :: ncid, varid
      integer, allocatable :: dcache(:,:)
      logical :: fexists


      if (p_is_io) then

         DO iblkme = 1, gblock%nblkme 
            iblk = gblock%xblkme(iblkme)
            jblk = gblock%yblkme(iblkme)
            IF (grid%xcnt(iblk) == 0) cycle
            IF (grid%ycnt(jblk) == 0) cycle
                     
            rdata%blk(iblk,jblk)%val(:,:) = 0

            inorth = grid%ydsp(jblk) + 1
            isouth = grid%ydsp(jblk) + grid%ycnt(jblk)

            iwest = grid%xdsp(iblk) + 1
            ieast = grid%xdsp(iblk) + grid%xcnt(iblk)
            if (ieast > nxglb) ieast = ieast - nxglb

            ibox = grid%xdsp(iblk)/nxbox + 1
            jbox = grid%ydsp(jblk)/nybox + 1
            ibox0 = ibox

            do while (.true.)

               CALL this_block_and_move_to_next (dir_modis, grid, &
                  isouth, inorth, iwest, ieast, ibox, jbox, ibox0, &
                  i0, i1, j0, j1, il0, il1, jl0, jl1, &
                  file_modis)

               inquire(file=file_modis, exist=fexists)
               if (fexists) then
                  allocate (dcache (i1-i0+1,j1-j0+1))

                  CALL nccheck( nf90_open(trim(file_modis), NF90_NOWRITE, ncid) )
                  CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )
                  CALL nccheck( nf90_get_var(ncid, varid, dcache, (/i0,j0/), (/i1-i0+1,j1-j0+1/)) )
                  CALL nccheck( nf90_close(ncid) )

                  rdata%blk(iblk,jblk)%val(il0:il1,jl0:jl1) = dcache

                  deallocate (dcache)
               end if

               IF (jbox == -1) exit

            end do

         end do

      end if

   end subroutine modis_read_data_int32

   ! -----
   subroutine modis_read_data_real8 (dir_modis, dataname, grid, rdata)

      USE precision
      use spmd_task
      use mod_block
      use mod_grid
      use mod_data_type
      use netcdf
      implicit none

      character (len=*), intent(in) :: dir_modis
      character (len=*), intent(in) :: dataname
      type (grid_type),  intent(in) :: grid
      type (block_data_real8_2d), intent(inout) :: rdata

      ! Local variables
      integer :: iblkme, iblk, jblk, isouth, inorth, iwest, ieast, ibox, jbox, ibox0
      integer :: i0, i1, j0, j1, il0, il1, jl0, jl1
      character(len=256) :: file_modis
      INTEGER :: ncid, varid
      REAL(r8), allocatable :: dcache(:,:)
      logical :: fexists


      if (p_is_io) then

         DO iblkme = 1, gblock%nblkme 
            iblk = gblock%xblkme(iblkme)
            jblk = gblock%yblkme(iblkme)
            IF (grid%xcnt(iblk) == 0) cycle
            IF (grid%ycnt(jblk) == 0) cycle
                     
            rdata%blk(iblk,jblk)%val(:,:) = 0

            inorth = grid%ydsp(jblk) + 1
            isouth = grid%ydsp(jblk) + grid%ycnt(jblk)

            iwest = grid%xdsp(iblk) + 1
            ieast = grid%xdsp(iblk) + grid%xcnt(iblk)
            if (ieast > nxglb) ieast = ieast - nxglb

            ibox = grid%xdsp(iblk)/nxbox + 1
            jbox = grid%ydsp(jblk)/nybox + 1
            ibox0 = ibox

            do while (.true.)

               CALL this_block_and_move_to_next (dir_modis, grid, &
                  isouth, inorth, iwest, ieast, ibox, jbox, ibox0, &
                  i0, i1, j0, j1, il0, il1, jl0, jl1, &
                  file_modis)

               inquire(file=file_modis, exist=fexists)
               if (fexists) then
                  allocate (dcache (i1-i0+1,j1-j0+1))

                  CALL nccheck( nf90_open(trim(file_modis), NF90_NOWRITE, ncid) )
                  CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )
                  CALL nccheck( nf90_get_var(ncid, varid, dcache, (/i0,j0/), (/i1-i0+1,j1-j0+1/)) )
                  CALL nccheck( nf90_close(ncid) )

                  rdata%blk(iblk,jblk)%val(il0:il1,jl0:jl1) = dcache

                  deallocate(dcache)
               end if

               IF (jbox == -1) exit

            end do

         end do

      end if

   end subroutine modis_read_data_real8
   
   ! -----
   subroutine modis_read_data_pft (dir_modis, dataname, grid, rdata)

      USE precision
      use spmd_task
      use mod_block
      use mod_grid
      use mod_data_type
      use netcdf
      implicit none

      character (len=*), intent(in) :: dir_modis
      character (len=*), intent(in) :: dataname
      type (grid_type),  intent(in) :: grid
      
      type (block_data_real8_3d), intent(inout) :: rdata

      ! Local variables
      integer :: iblkme, iblk, jblk, isouth, inorth, iwest, ieast, ibox, jbox, ibox0
      integer :: i0, i1, j0, j1, il0, il1, jl0, jl1
      character(len=256) :: file_modis
      INTEGER :: ncid, varid
      REAL(r8), allocatable :: dcache(:,:,:)
      logical :: fexists
      INTEGER :: ipft

      if (p_is_io) then

         DO iblkme = 1, gblock%nblkme 
            iblk = gblock%xblkme(iblkme)
            jblk = gblock%yblkme(iblkme)
            IF (grid%xcnt(iblk) == 0) cycle
            IF (grid%ycnt(jblk) == 0) cycle
         
            rdata%blk(iblk,jblk)%val(:,:,:) = 0

            inorth = grid%ydsp(jblk) + 1
            isouth = grid%ydsp(jblk) + grid%ycnt(jblk)

            iwest = grid%xdsp(iblk) + 1
            ieast = grid%xdsp(iblk) + grid%xcnt(iblk)
            if (ieast > nxglb) ieast = ieast - nxglb

            ibox = grid%xdsp(iblk)/nxbox + 1
            jbox = grid%ydsp(jblk)/nybox + 1
            ibox0 = ibox

            do while (.true.)

               CALL this_block_and_move_to_next (dir_modis, grid, &
                  isouth, inorth, iwest, ieast, ibox, jbox, ibox0, &
                  i0, i1, j0, j1, il0, il1, jl0, jl1, &
                  file_modis)

               inquire(file=file_modis, exist=fexists)
               if (fexists) then
                  allocate (dcache (i1-i0+1,j1-j0+1,0:N_PFT_modis-1))

                  CALL nccheck( nf90_open(trim(file_modis), NF90_NOWRITE, ncid) )
                  CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )
                  CALL nccheck( nf90_get_var(ncid, varid, dcache, &
                     (/i0,j0,1/), (/i1-i0+1,j1-j0+1,N_PFT_modis/)) )
                  CALL nccheck( nf90_close(ncid) )

                  DO ipft = 0, N_PFT_modis-1
                     rdata%blk(iblk,jblk)%val(ipft,il0:il1,jl0:jl1) = dcache(:,:,ipft)
                  ENDDO

                  deallocate (dcache)
               end if

               IF (jbox == -1) exit

            end do

         end do

      end if

   end subroutine modis_read_data_pft

   ! -----
   subroutine modis_read_data_time (dir_modis, dataname, grid, time, rdata)

      USE precision
      use spmd_task
      use mod_block
      use mod_grid
      use mod_data_type
      use netcdf
      implicit none

      character (len=*), intent(in) :: dir_modis
      character (len=*), intent(in) :: dataname
      type (grid_type),  intent(in) :: grid
      integer, intent(in) :: time
      
      type (block_data_real8_2d), intent(inout) :: rdata

      ! Local variables
      integer :: iblkme, iblk, jblk, isouth, inorth, iwest, ieast, ibox, jbox, ibox0
      integer :: i0, i1, j0, j1, il0, il1, jl0, jl1
      character(len=256) :: file_modis
      INTEGER :: ncid, varid
      REAL(r8), allocatable :: dcache(:,:)
      logical :: fexists

      if (p_is_io) then

                     
         DO iblkme = 1, gblock%nblkme 
            iblk = gblock%xblkme(iblkme)
            jblk = gblock%yblkme(iblkme)
            IF (grid%xcnt(iblk) == 0) cycle
            IF (grid%ycnt(jblk) == 0) cycle
         
            rdata%blk(iblk,jblk)%val(:,:) = 0

            inorth = grid%ydsp(jblk) + 1
            isouth = grid%ydsp(jblk) + grid%ycnt(jblk)

            iwest = grid%xdsp(iblk) + 1
            ieast = grid%xdsp(iblk) + grid%xcnt(iblk)
            if (ieast > nxglb) ieast = ieast - nxglb

            ibox = grid%xdsp(iblk)/nxbox + 1
            jbox = grid%ydsp(jblk)/nybox + 1
            ibox0 = ibox

            do while (.true.)

               CALL this_block_and_move_to_next (dir_modis, grid, &
                  isouth, inorth, iwest, ieast, ibox, jbox, ibox0, &
                  i0, i1, j0, j1, il0, il1, jl0, jl1, &
                  file_modis)

               inquire(file=file_modis, exist=fexists)
               if (fexists) then
                  allocate (dcache (i1-i0+1,j1-j0+1))

                  CALL nccheck( nf90_open(trim(file_modis), NF90_NOWRITE, ncid) )
                  CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )
                  CALL nccheck( nf90_get_var(ncid, varid, dcache, &
                     (/i0,j0,time/), (/i1-i0+1,j1-j0+1,1/)) )
                  CALL nccheck( nf90_close(ncid) )

                  rdata%blk(iblk,jblk)%val(il0:il1,jl0:jl1) = dcache

                  deallocate (dcache)
               end if

               IF (jbox == -1) exit

            end do

         end do

      end if

   end subroutine modis_read_data_time

   ! -----
   subroutine modis_read_data_pft_time (dir_modis, dataname, grid, time, rdata)

      USE precision
      use spmd_task
      use mod_block
      use mod_grid
      use mod_data_type
      use netcdf
      implicit none

      character (len=*), intent(in) :: dir_modis
      character (len=*), intent(in) :: dataname
      type (grid_type),  intent(in) :: grid
      INTEGER, intent(in) :: time
      
      type (block_data_real8_3d), intent(inout) :: rdata

      ! Local variables
      integer :: iblkme, iblk, jblk, isouth, inorth, iwest, ieast, ibox, jbox, ibox0
      integer :: i0, i1, j0, j1, il0, il1, jl0, jl1
      character(len=256) :: file_modis
      INTEGER :: ncid, varid
      REAL(r8), allocatable :: dcache(:,:,:)
      logical :: fexists
      INTEGER :: ipft

      if (p_is_io) then

         DO iblkme = 1, gblock%nblkme 
            iblk = gblock%xblkme(iblkme)
            jblk = gblock%yblkme(iblkme)
            IF (grid%xcnt(iblk) == 0) cycle
            IF (grid%ycnt(jblk) == 0) cycle
                     
            rdata%blk(iblk,jblk)%val(:,:,:) = 0

            inorth = grid%ydsp(jblk) + 1
            isouth = grid%ydsp(jblk) + grid%ycnt(jblk)

            iwest = grid%xdsp(iblk) + 1
            ieast = grid%xdsp(iblk) + grid%xcnt(iblk)
            if (ieast > nxglb) ieast = ieast - nxglb

            ibox = grid%xdsp(iblk)/nxbox + 1
            jbox = grid%ydsp(jblk)/nybox + 1
            ibox0 = ibox

            do while (.true.)

               CALL this_block_and_move_to_next (dir_modis, grid, &
                  isouth, inorth, iwest, ieast, ibox, jbox, ibox0, &
                  i0, i1, j0, j1, il0, il1, jl0, jl1, &
                  file_modis)

               inquire(file=file_modis, exist=fexists)
               if (fexists) then
                  allocate (dcache (i1-i0+1,j1-j0+1,0:N_PFT_modis-1))

                  CALL nccheck( nf90_open(trim(file_modis), NF90_NOWRITE, ncid) )
                  CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )
                  CALL nccheck( nf90_get_var(ncid, varid, dcache, &
                     (/i0,j0,1,time/), (/i1-i0+1,j1-j0+1,N_PFT_modis,1/)) )
                  CALL nccheck( nf90_close(ncid) )

                  DO ipft = 0, N_PFT_modis-1
                     rdata%blk(iblk,jblk)%val(ipft,il0:il1,jl0:jl1) = dcache(:,:,ipft)
                  ENDDO

                  deallocate (dcache)
               end if

               IF (jbox == -1) exit

            end do

         end do

      end if

   end subroutine modis_read_data_pft_time

end module mod_modis_data
