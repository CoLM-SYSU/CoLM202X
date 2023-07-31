#include <define.h>

module MOD_5x5DataReadin

   !-----------------------------------------------------------------------
   ! DESCRIPTION:
   !
   !    Reading data in netCDF files by 5 degree blocks.   
   !
   !    The file name gives the boundaries of the block. 
   !    For example, file "RG_65_75_60_80.URB2010.nc" stores data in region 
   !    from 65N to 60N and 75E to 80E. 
   !
   !    Notice that:
   !    1. Subroutines loop over all 5 degree blocks in simulation region.
   !    2. Latitude in files is from north to south. 
   !    3. "read_5x5_data_pft" reads data with dimension "pft" and permute 
   !       dimension (lon,lat,pft) in files to (pft,lon,lat) in variables.
   !    4. "read_5x5_data_time" reads data with dimension "time" 
   !       at given time.
   !    5. "read_5x5_data_pft_time" reads data with dimension "pft" and "time" 
   !       at given time and permute dimension (lon,lat,pft) in files 
   !       to (pft,lon,lat) in variables.
   ! 
   ! Created by Shupeng Zhang, May 2023
   !-----------------------------------------------------------------------

   USE MOD_NetCDFSerial, only : nccheck
   implicit none

   INTEGER, parameter :: N_PFT_modis = 16

   interface read_5x5_data
      MODULE procedure read_5x5_data_int32
      MODULE procedure read_5x5_data_real8
   END interface read_5x5_data

   PUBLIC :: read_5x5_data_pft
   PUBLIC :: read_5x5_data_time
   PUBLIC :: read_5x5_data_pft_time

contains

   ! -----
   subroutine this_block_and_move_to_next ( &
         dir_5x5, sfx, nxbox, nybox, nxglb, isouth, inorth, iwest, ieast, &
         ibox, jbox, ibox0, i0, i1, j0, j1, il0, il1, jl0, jl1, &
         file_5x5)

      use MOD_Grid
      implicit none

      character (len=*), intent(in) :: dir_5x5
      character (len=*), intent(in) :: sfx

      integer, intent(in)    :: nxbox, nybox, nxglb
      integer, intent(in)    :: isouth, inorth, iwest, ieast
      INTEGER, intent(inout) :: ibox, jbox, ibox0
      integer, intent(out)   :: i0,  i1,  j0,  j1
      integer, intent(out)   :: il0, il1, jl0, jl1

      character (len=*), intent(out) :: file_5x5

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

      file_5x5 = trim(dir_5x5) // '/RG'
      write(str, '(I4)') (19-jbox)*5
      file_5x5 = trim(file_5x5) // '_' // trim(adjustl(str))
      write(str, '(I4)') (ibox-37)*5
      file_5x5 = trim(file_5x5) // '_' // trim(adjustl(str))
      write(str, '(I4)') (18-jbox)*5
      file_5x5 = trim(file_5x5) // '_' // trim(adjustl(str))
      write(str, '(I4)') (ibox-36)*5
      file_5x5 = trim(file_5x5) // '_' // trim(adjustl(str))
      file_5x5 = trim(file_5x5) // '.' // trim(sfx) // '.nc'

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
   subroutine read_5x5_data_int32 (dir_5x5, sfx, grid, dataname, rdata)

      use MOD_SPMD_Task
      use MOD_Block
      use MOD_Grid
      use MOD_DataType
      use netcdf
      implicit none

      character (len=*), intent(in) :: dir_5x5
      character (len=*), intent(in) :: sfx
      type (grid_type),  intent(in) :: grid

      character (len=*), intent(in) :: dataname
      type (block_data_int32_2d), intent(inout) :: rdata

      ! Local variables
      integer :: nxbox, nybox, nxglb, nyglb
      integer :: iblkme, iblk, jblk, isouth, inorth, iwest, ieast, ibox, jbox, ibox0
      integer :: i0, i1, j0, j1, il0, il1, jl0, jl1
      character(len=256) :: file_5x5
      INTEGER :: ncid, varid
      integer, allocatable :: dcache(:,:)
      logical :: fexists

      nxglb = grid%nlon
      nyglb = grid%nlat

      nxbox = nxglb / 360 * 5
      nybox = nyglb / 180 * 5

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

               CALL this_block_and_move_to_next ( &
                  dir_5x5, sfx, nxbox, nybox, nxglb, isouth, inorth, iwest, ieast, &
                  ibox, jbox, ibox0, i0, i1, j0, j1, il0, il1, jl0, jl1, &
                  file_5x5)

               inquire(file=file_5x5, exist=fexists)
               if (fexists) then
                  allocate (dcache (i1-i0+1,j1-j0+1))

                  CALL nccheck( nf90_open(trim(file_5x5), NF90_NOWRITE, ncid) )
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

   end subroutine read_5x5_data_int32

   ! -----
   subroutine read_5x5_data_real8 (dir_5x5, sfx, grid, dataname, rdata)

      USE MOD_Precision
      use MOD_SPMD_Task
      use MOD_Block
      use MOD_Grid
      use MOD_DataType
      use netcdf
      implicit none

      character (len=*), intent(in) :: dir_5x5
      character (len=*), intent(in) :: sfx
      type (grid_type),  intent(in) :: grid

      character (len=*), intent(in) :: dataname
      type (block_data_real8_2d), intent(inout) :: rdata

      ! Local variables
      integer :: nxbox, nybox, nxglb, nyglb
      integer :: iblkme, iblk, jblk, isouth, inorth, iwest, ieast, ibox, jbox, ibox0
      integer :: i0, i1, j0, j1, il0, il1, jl0, jl1
      character(len=256) :: file_5x5
      INTEGER :: ncid, varid
      REAL(r8), allocatable :: dcache(:,:)
      logical :: fexists

      nxglb = grid%nlon
      nyglb = grid%nlat

      nxbox = nxglb / 360 * 5
      nybox = nyglb / 180 * 5

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

               CALL this_block_and_move_to_next ( &
                  dir_5x5, sfx, nxbox, nybox, nxglb, isouth, inorth, iwest, ieast, &
                  ibox, jbox, ibox0, i0, i1, j0, j1, il0, il1, jl0, jl1, &
                  file_5x5)

               inquire(file=file_5x5, exist=fexists)
               if (fexists) then
                  allocate (dcache (i1-i0+1,j1-j0+1))

                  CALL nccheck( nf90_open(trim(file_5x5), NF90_NOWRITE, ncid) )
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

   end subroutine read_5x5_data_real8

   ! -----
   subroutine read_5x5_data_pft (dir_5x5, sfx, grid, dataname, rdata)

      USE MOD_Precision
      use MOD_SPMD_Task
      use MOD_Block
      use MOD_Grid
      use MOD_DataType
      use netcdf
      implicit none

      character (len=*), intent(in) :: dir_5x5
      character (len=*), intent(in) :: sfx
      type (grid_type),  intent(in) :: grid

      character (len=*), intent(in) :: dataname
      type (block_data_real8_3d), intent(inout) :: rdata

      ! Local variables
      integer :: nxbox, nybox, nxglb, nyglb
      integer :: iblkme, iblk, jblk, isouth, inorth, iwest, ieast, ibox, jbox, ibox0
      integer :: i0, i1, j0, j1, il0, il1, jl0, jl1
      character(len=256) :: file_5x5
      INTEGER :: ncid, varid
      REAL(r8), allocatable :: dcache(:,:,:)
      logical :: fexists
      INTEGER :: ipft

      nxglb = grid%nlon
      nyglb = grid%nlat

      nxbox = nxglb / 360 * 5
      nybox = nyglb / 180 * 5

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

               CALL this_block_and_move_to_next ( &
                  dir_5x5, sfx, nxbox, nybox, nxglb, isouth, inorth, iwest, ieast, &
                  ibox, jbox, ibox0, i0, i1, j0, j1, il0, il1, jl0, jl1, &
                  file_5x5)

               inquire(file=file_5x5, exist=fexists)
               if (fexists) then
                  allocate (dcache (i1-i0+1,j1-j0+1,0:N_PFT_modis-1))

                  CALL nccheck( nf90_open(trim(file_5x5), NF90_NOWRITE, ncid) )
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

   end subroutine read_5x5_data_pft

   ! -----
   subroutine read_5x5_data_time (dir_5x5, sfx, grid, dataname, time, rdata)

      USE MOD_Precision
      use MOD_SPMD_Task
      use MOD_Block
      use MOD_Grid
      use MOD_DataType
      use netcdf
      implicit none

      character (len=*), intent(in) :: dir_5x5
      character (len=*), intent(in) :: sfx
      type (grid_type),  intent(in) :: grid

      character (len=*), intent(in) :: dataname
      integer, intent(in) :: time
      type (block_data_real8_2d), intent(inout) :: rdata

      ! Local variables
      integer :: nxbox, nybox, nxglb, nyglb
      integer :: iblkme, iblk, jblk, isouth, inorth, iwest, ieast, ibox, jbox, ibox0
      integer :: i0, i1, j0, j1, il0, il1, jl0, jl1
      character(len=256) :: file_5x5
      INTEGER :: ncid, varid
      REAL(r8), allocatable :: dcache(:,:)
      logical :: fexists

      nxglb = grid%nlon
      nyglb = grid%nlat

      nxbox = nxglb / 360 * 5
      nybox = nyglb / 180 * 5

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

               CALL this_block_and_move_to_next ( &
                  dir_5x5, sfx, nxbox, nybox, nxglb, isouth, inorth, iwest, ieast, &
                  ibox, jbox, ibox0, i0, i1, j0, j1, il0, il1, jl0, jl1, &
                  file_5x5)

               inquire(file=file_5x5, exist=fexists)
               if (fexists) then
                  allocate (dcache (i1-i0+1,j1-j0+1))

                  CALL nccheck( nf90_open(trim(file_5x5), NF90_NOWRITE, ncid) )
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

   end subroutine read_5x5_data_time

   ! -----
   subroutine read_5x5_data_pft_time (dir_5x5, sfx, grid, dataname, time, rdata)

      USE MOD_Precision
      use MOD_SPMD_Task
      use MOD_Block
      use MOD_Grid
      use MOD_DataType
      use netcdf
      implicit none

      character (len=*), intent(in) :: dir_5x5
      character (len=*), intent(in) :: sfx
      type (grid_type),  intent(in) :: grid

      character (len=*), intent(in) :: dataname
      INTEGER, intent(in) :: time
      type (block_data_real8_3d), intent(inout) :: rdata

      ! Local variables
      integer :: nxbox, nybox, nxglb, nyglb
      integer :: iblkme, iblk, jblk, isouth, inorth, iwest, ieast, ibox, jbox, ibox0
      integer :: i0, i1, j0, j1, il0, il1, jl0, jl1
      character(len=256) :: file_5x5
      INTEGER :: ncid, varid
      REAL(r8), allocatable :: dcache(:,:,:)
      logical :: fexists
      INTEGER :: ipft

      nxglb = grid%nlon
      nyglb = grid%nlat

      nxbox = nxglb / 360 * 5
      nybox = nyglb / 180 * 5

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

               CALL this_block_and_move_to_next ( &
                  dir_5x5, sfx, nxbox, nybox, nxglb, isouth, inorth, iwest, ieast, &
                  ibox, jbox, ibox0, i0, i1, j0, j1, il0, il1, jl0, jl1, &
                  file_5x5)

               inquire(file=file_5x5, exist=fexists)
               if (fexists) then
                  allocate (dcache (i1-i0+1,j1-j0+1,0:N_PFT_modis-1))

                  CALL nccheck( nf90_open(trim(file_5x5), NF90_NOWRITE, ncid) )
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

   end subroutine read_5x5_data_pft_time

end module MOD_5x5DataReadin
