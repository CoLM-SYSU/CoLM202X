#include <define.h>

MODULE MOD_CatchmentDataReadin

!--------------------------------------------------------------------------------------
! !DESCRIPTION:
!
!    Reading preprocessed MERIT Hydro data and generated catchment data in netcdf files.
!
!    1. If "in_one_file" is false, then the data is orgnized by 5 degree blocks.
!    The file name gives the southwest corner of the block.
!    For example, file "n60e075.nc" stores data in region from 65N to 60N and 75E to 80E,
!    Subroutines loop over all 5 degree blocks in simulation region.
!
!    2. Data is saved in variables with types of "block_data_xxxxx_xd".
!
!    3. Latitude in files is from north to south.
!
!  Created by Shupeng Zhang, May 2023
!--------------------------------------------------------------------------------------

   IMPLICIT NONE

   integer, parameter :: nxhbox = 6000
   integer, parameter :: nyhbox = 6000
   integer, parameter :: nxhglb = 432000
   integer, parameter :: nyhglb = 216000

   INTERFACE catchment_data_read
      MODULE procedure catchment_data_read_int32
      MODULE procedure catchment_data_read_real8
   END INTERFACE catchment_data_read

CONTAINS

   ! -----
   SUBROUTINE catchment_data_read_int32 (file_meshdata_in, dataname, grid, rdata_int32, spv)

   USE MOD_Grid
   USE MOD_DataType
   IMPLICIT NONE

   character (len=*), intent(in) :: file_meshdata_in
   character (len=*), intent(in) :: dataname
   type (grid_type),  intent(in) :: grid

   type (block_data_int32_2d), intent(inout) :: rdata_int32
   integer,  intent(in), optional :: spv

      IF (present(spv)) THEN
         CALL catchment_data_read_general (file_meshdata_in, dataname, grid, &
            rdata_int32 = rdata_int32, spv_i4 = spv)
      ELSE
         CALL catchment_data_read_general (file_meshdata_in, dataname, grid, &
            rdata_int32 = rdata_int32)
      ENDIF

   END SUBROUTINE catchment_data_read_int32

   ! -----
   SUBROUTINE catchment_data_read_real8 (file_meshdata_in, dataname, grid, rdata_real8, spv)

   USE MOD_Grid
   USE MOD_DataType
   IMPLICIT NONE

   character (len=*), intent(in) :: file_meshdata_in
   character (len=*), intent(in) :: dataname
   type (grid_type),  intent(in) :: grid

   type (block_data_real8_2d), intent(inout) :: rdata_real8
   real(r8), intent(in), optional :: spv

      IF (present(spv)) THEN
         CALL catchment_data_read_general (file_meshdata_in, dataname, grid, &
            rdata_real8 = rdata_real8, spv_r8 = spv)
      ELSE
         CALL catchment_data_read_general (file_meshdata_in, dataname, grid, &
            rdata_real8 = rdata_real8)
      ENDIF

   END SUBROUTINE catchment_data_read_real8

   ! -----
   SUBROUTINE catchment_data_read_general (file_meshdata_in, dataname, grid, &
      rdata_int32, spv_i4, rdata_real8, spv_r8)

   USE MOD_SPMD_Task
   USE MOD_Block
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_Utils
   USE MOD_NetCDFSerial
   IMPLICIT NONE

   character (len=*), intent(in) :: file_meshdata_in
   character (len=*), intent(in) :: dataname
   type (grid_type),  intent(in) :: grid

   type (block_data_int32_2d), intent(inout), optional :: rdata_int32
   integer,  intent(in), optional :: spv_i4

   type (block_data_real8_2d), intent(inout), optional :: rdata_real8
   real(r8), intent(in), optional :: spv_r8


   ! Local Variables
   logical :: in_one_file
   integer :: nlat, nlon, ilon
   integer :: iblkme, iblk, jblk, isouth, inorth, iwest, ieast, ibox, jbox
   integer :: xdsp, ydsp, i0, i1, j0, j1, il0, il1, jl0, jl1
   integer :: i0min, i1max, if0, if1, jf0, jf1, i0next, i1next
   character(len=256) :: file_mesh, path_mesh
   character(len=3)   :: pre1
   character(len=4)   :: pre2
   integer,  allocatable :: dcache_i4(:,:)
   real(r8), allocatable :: dcache_r8(:,:)
   real(r8), allocatable :: latitude(:), longitude(:)
   logical :: fexists

      IF (p_is_io) THEN

         IF (p_iam_io == p_root) THEN
            IF (grid%yinc == 1) THEN
               write(*,*) 'Warning: latitude in catchment data should be from north to south.'
            ENDIF
         ENDIF

         IF (p_iam_io == p_root) THEN
            in_one_file = ncio_var_exist (file_meshdata_in, dataname)
         ENDIF
#ifdef USEMPI
         CALL mpi_bcast (in_one_file, 1, mpi_logical, p_root, p_comm_io, p_err)
#endif

         IF (in_one_file) THEN

            file_mesh = file_meshdata_in

            CALL ncio_read_serial (file_mesh, 'lat', latitude)
            CALL ncio_read_serial (file_mesh, 'lon', longitude)

            nlat = size(latitude )
            nlon = size(longitude)

            isouth = find_nearest_south (latitude(nlat), grid%nlat, grid%lat_s)
            inorth = find_nearest_north (latitude(1),    grid%nlat, grid%lat_n)

            DO ilon = 1, nlon
               CALL normalize_longitude (longitude(ilon))
            ENDDO

            iwest = find_nearest_west (longitude(1),    grid%nlon, grid%lon_w)
            ieast = find_nearest_east (longitude(nlon), grid%nlon, grid%lon_e)

            DO iblkme = 1, gblock%nblkme
               iblk = gblock%xblkme(iblkme)
               jblk = gblock%yblkme(iblkme)

               IF (present(rdata_int32)) THEN
                  IF (present(spv_i4)) THEN
                     rdata_int32%blk(iblk,jblk)%val(:,:) = spv_i4
                  ELSE
                     rdata_int32%blk(iblk,jblk)%val(:,:) = 0
                  ENDIF
               ENDIF

               IF (present(rdata_real8)) THEN
                  IF (present(spv_r8)) THEN
                     rdata_real8%blk(iblk,jblk)%val(:,:) = spv_r8
                  ELSE
                     rdata_real8%blk(iblk,jblk)%val(:,:) = -1.e36_r8
                  ENDIF
               ENDIF

               IF ((inorth > grid%ydsp(jblk)+nyhbox) .or. (isouth < grid%ydsp(jblk)+1)) THEN
                  CYCLE
               ENDIF

               j0 = max(inorth, grid%ydsp(jblk)+1)
               j1 = min(isouth, grid%ydsp(jblk)+grid%ycnt(jblk))

               jl0 = j0 - grid%ydsp(jblk)
               jl1 = j1 - grid%ydsp(jblk)
               jf0 = j0 - inorth + 1
               jf1 = j1 - inorth + 1

               i0min = grid%xdsp(iblk) + 1
               i1max = grid%xdsp(iblk) + grid%xcnt(iblk)
               IF (i1max > grid%nlon) i1max = i1max - grid%nlon

               DO WHILE ((i0min /= i1max) .and. (.not. (lon_between_floor(grid%lon_w(i0min), &
                     grid%lon_w(iwest), grid%lon_e(ieast)))))
                  i0min = i0min + 1; IF (i0min > grid%nlon) i0min = 1
               ENDDO

               IF (lon_between_floor(grid%lon_w(i0min), grid%lon_w(iwest), grid%lon_e(ieast))) THEN
                  i0 = i0min
                  i1 = i0
                  i1next = i1 + 1; IF (i1next > grid%nlon) i1next = 1
                  DO WHILE ((i1 /= i1max) .and. &
                        lon_between_floor(grid%lon_w(i1next), grid%lon_w(iwest), grid%lon_e(ieast)))
                     i1 = i1next
                     i1next = i1 + 1; IF (i1next > grid%nlon) i1next = 1
                  ENDDO

                  if0 = i0 - iwest + 1;       IF (if0 <= 0) if0 = if0 + grid%nlon
                  if1 = i1 - iwest + 1;       IF (if1 <= 0) if1 = if1 + grid%nlon
                  il0 = i0 - grid%xdsp(iblk); IF (il0 <= 0) il0 = il0 + grid%nlon
                  il1 = i1 - grid%xdsp(iblk); IF (il1 <= 0) il1 = il1 + grid%nlon

                  IF (present(rdata_int32)) THEN
                     CALL ncio_read_part_serial (file_mesh, dataname, (/jf0,if0/), (/jf1,if1/), dcache_i4)
                     dcache_i4 = transpose(dcache_i4)

                     rdata_int32%blk(iblk,jblk)%val(il0:il1,jl0:jl1) = dcache_i4
                  ENDIF

                  IF (present(rdata_real8)) THEN
                     CALL ncio_read_part_serial (file_mesh, dataname, (/jf0,if0/), (/jf1,if1/), dcache_r8)
                     dcache_r8 = transpose(dcache_r8)

                     rdata_real8%blk(iblk,jblk)%val(il0:il1,jl0:jl1) = dcache_r8
                  ENDIF

               ENDIF

               IF (lon_between_ceil(grid%lon_e(i1max), grid%lon_w(iwest), grid%lon_e(ieast))) THEN
                  i1 = i1max
                  i0 = i1
                  i0next = i0 - 1; IF (i0next == 0) i0next = grid%nlon
                  DO WHILE ((i0 /= i0min) .and. &
                        lon_between_ceil(grid%lon_e(i0next), grid%lon_w(iwest), grid%lon_e(ieast)))
                     i0 = i0next
                     i0next = i0 - 1; IF (i0next == 0) i0next = grid%nlon
                  ENDDO

                  IF (i0 /= i0min) THEN
                     if0 = i0 - iwest + 1;       IF (if0 <= 0) if0 = if0 + grid%nlon
                     if1 = i1 - iwest + 1;       IF (if1 <= 0) if1 = if1 + grid%nlon
                     il0 = i0 - grid%xdsp(iblk); IF (il0 <= 0) il0 = il0 + grid%nlon
                     il1 = i1 - grid%xdsp(iblk); IF (il1 <= 0) il1 = il1 + grid%nlon

                     IF (present(rdata_int32)) THEN
                        CALL ncio_read_part_serial (file_mesh, dataname, (/jf0,if0/), (/jf1,if1/), dcache_i4)
                        dcache_i4 = transpose(dcache_i4)

                        rdata_int32%blk(iblk,jblk)%val(il0:il1,jl0:jl1) = dcache_i4
                     ENDIF

                     IF (present(rdata_real8)) THEN
                        CALL ncio_read_part_serial (file_mesh, dataname, (/jf0,if0/), (/jf1,if1/), dcache_r8)
                        dcache_r8 = transpose(dcache_r8)

                        rdata_real8%blk(iblk,jblk)%val(il0:il1,jl0:jl1) = dcache_r8
                     ENDIF
                  ENDIF
               ENDIF

            ENDDO

         ELSE

            ! remove suffix ".nc"
            path_mesh = file_meshdata_in(1:len_trim(file_meshdata_in)-3)

            DO iblkme = 1, gblock%nblkme
               iblk = gblock%xblkme(iblkme)
               jblk = gblock%yblkme(iblkme)

               IF (present(rdata_int32)) THEN
                  IF (present(spv_i4)) THEN
                     rdata_int32%blk(iblk,jblk)%val(:,:) = spv_i4
                  ELSE
                     rdata_int32%blk(iblk,jblk)%val(:,:) = 0
                  ENDIF
               ENDIF

               IF (present(rdata_real8)) THEN
                  IF (present(spv_r8)) THEN
                     rdata_real8%blk(iblk,jblk)%val(:,:) = spv_r8
                  ELSE
                     rdata_real8%blk(iblk,jblk)%val(:,:) = -1.e36_r8
                  ENDIF
               ENDIF

               inorth = grid%ydsp(jblk) + 1
               isouth = grid%ydsp(jblk) + grid%ycnt(jblk)

               iwest = grid%xdsp(iblk) + 1
               ieast = grid%xdsp(iblk) + grid%xcnt(iblk)
               IF (ieast > nxhglb) ieast = ieast - nxhglb

               ibox = grid%xdsp(iblk)/nxhbox + 1
               jbox = grid%ydsp(jblk)/nyhbox + 1

               DO WHILE (.true.)

                  xdsp = (ibox-1) * nxhbox
                  ydsp = (jbox-1) * nyhbox

                  j0 = max(inorth-ydsp, 1)
                  j1 = min(isouth-ydsp, nyhbox)
                  jl0 = j0 + ydsp - inorth + 1
                  jl1 = j1 + ydsp - inorth + 1

                  IF (ieast >= iwest) THEN
                     i0 = max(iwest-xdsp, 1)
                     i1 = min(ieast-xdsp, nxhbox)
                     il0 = i0 + xdsp - iwest + 1
                     il1 = i1 + xdsp - iwest + 1
                  ELSE
                     IF (iwest <= xdsp+nxhbox) THEN
                        i0 = max(iwest-xdsp, 1)
                        i1 = nxhbox
                        il0 = i0 + xdsp - iwest + 1
                        il1 = i1 + xdsp - iwest + 1
                     ELSE
                        i0 = 1
                        i1 = min(ieast-xdsp, nxhbox)
                        il0 = i0 + xdsp + nxhglb - iwest + 1
                        il1 = i1 + xdsp + nxhglb - iwest + 1
                     ENDIF
                  ENDIF

                  IF (jbox <= 18) THEN
                     write (pre1,'(A1,I2.2)') 'n', (18-jbox)*5
                  ELSE
                     write (pre1,'(A1,I2.2)') 's', (jbox-18)*5
                  ENDIF

                  IF (ibox <= 36) THEN
                     write (pre2,'(A1,I3.3)') 'w', (37-ibox)*5
                  ELSE
                     write (pre2,'(A1,I3.3)') 'e', (ibox-37)*5
                  ENDIF

                  file_mesh = trim(path_mesh) // '/' // trim(pre1) // trim(pre2) // '.nc'

                  inquire(file=file_mesh, exist=fexists)
                  IF (fexists) THEN
                     IF (present(rdata_int32)) THEN
                        CALL ncio_read_part_serial (file_mesh, dataname, (/j0,i0/), (/j1,i1/), dcache_i4)
                        dcache_i4 = transpose(dcache_i4)

                        rdata_int32%blk(iblk,jblk)%val(il0:il1,jl0:jl1) = dcache_i4
                     ENDIF

                     IF (present(rdata_real8)) THEN
                        CALL ncio_read_part_serial (file_mesh, dataname, (/j0,i0/), (/j1,i1/), dcache_r8)
                        dcache_r8 = transpose(dcache_r8)

                        rdata_real8%blk(iblk,jblk)%val(il0:il1,jl0:jl1) = dcache_r8
                     ENDIF
                  ENDIF

                  IF ((ieast >= xdsp + 1) .and. (ieast <= xdsp + nxhbox)) THEN
                     IF (isouth <= ydsp + nyhbox) THEN
                        EXIT
                     ELSE
                        ibox = grid%xdsp(iblk)/nxhbox + 1
                        jbox = jbox + 1
                     ENDIF
                  ELSE
                     ibox = mod(ibox, nxhglb/nxhbox) + 1
                  ENDIF

               ENDDO
            ENDDO
         ENDIF
      ENDIF

   END SUBROUTINE catchment_data_read_general

END MODULE MOD_CatchmentDataReadin
