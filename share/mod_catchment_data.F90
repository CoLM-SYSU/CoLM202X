#include <define.h>

MODULE mod_catchment_data

   IMPLICIT NONE

   INTEGER, parameter :: nxhbox = 6000
   INTEGER, parameter :: nyhbox = 6000
   INTEGER, parameter :: nxhglb = 432000
   INTEGER, parameter :: nyhglb = 216000

CONTAINS

   ! -----
   SUBROUTINE catchment_data_read (path_hydro, dataname, grid, rdata, in_one_file, spv)

      USE spmd_task
      USE mod_block
      USE mod_grid
      USE mod_data_type
      USE mod_utils
      USE ncio_serial
      IMPLICIT NONE

      CHARACTER (len=*), intent(in) :: path_hydro
      CHARACTER (len=*), intent(in) :: dataname
      TYPE (grid_type),  intent(in) :: grid
      TYPE (block_data_int32_2d), intent(inout) :: rdata
      LOGICAL, intent(in) :: in_one_file
      INTEGER, intent(in), optional :: spv

      ! Local Variables
      INTEGER :: nlat, nlon, ilon, ilat
      INTEGER :: iblkme, iblk, jblk, isouth, inorth, iwest, ieast, ibox, jbox
      INTEGER :: xdsp, ydsp, i0, i1, j0, j1, il0, il1, jl0, jl1
      INTEGER :: i0min, i1max, if0, if1, jf0, jf1, i0next, i1next
      CHARACTER(len=256) :: file_hydro
      CHARACTER(len=3)   :: pre1
      CHARACTER(len=4)   :: pre2
      INTEGER :: ncid, varid
      INTEGER,  allocatable :: dcache(:,:)
      REAL(r8), allocatable :: latitude(:), longitude(:)
      LOGICAL :: fexists

      IF (in_one_file) THEN

         file_hydro = path_hydro

         CALL ncio_read_bcast_serial (file_hydro, 'latitude',  latitude)
         CALL ncio_read_bcast_serial (file_hydro, 'longitude', longitude)

         IF (p_is_io) THEN

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

               IF (present(spv)) THEN
                  rdata%blk(iblk,jblk)%val(:,:) = spv
               ELSE
                  rdata%blk(iblk,jblk)%val(:,:) = 0
               ENDIF

               IF ((inorth > grid%ydsp(jblk)+nyhbox) .or. (isouth < grid%ydsp(jblk)+1)) THEN
                  cycle
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

                  CALL ncio_read_part_serial (file_hydro, dataname, (/jf0,if0/), (/jf1,if1/), dcache)
                  dcache = transpose(dcache)
                  rdata%blk(iblk,jblk)%val(il0:il1,jl0:jl1) = dcache
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

                     CALL ncio_read_part_serial (file_hydro, dataname, &
                        (/jf0,if0/), (/jf1,if1/), dcache)
                     dcache = transpose(dcache)
                     rdata%blk(iblk,jblk)%val(il0:il1,jl0:jl1) = dcache
                  ENDIF
               ENDIF

            ENDDO
         ENDIF
                     
      ELSE

         IF (p_is_io) THEN

            DO iblkme = 1, gblock%nblkme 
               iblk = gblock%xblkme(iblkme)
               jblk = gblock%yblkme(iblkme)

               IF (present(spv)) THEN
                  rdata%blk(iblk,jblk)%val(:,:) = spv
               ELSE
                  rdata%blk(iblk,jblk)%val(:,:) = 0
               ENDIF

               inorth = grid%ydsp(jblk) + 1
               isouth = grid%ydsp(jblk) + grid%ycnt(jblk)

               iwest = grid%xdsp(iblk) + 1
               ieast = grid%xdsp(iblk) + grid%xcnt(iblk)
               IF (ieast > nxhglb) ieast = ieast - nxhglb

               ibox = grid%xdsp(iblk)/nxhbox + 1
               jbox = grid%ydsp(jblk)/nyhbox + 1

               DO while (.true.)

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

                  file_hydro = trim(path_hydro) // '/' // trim(pre1) // trim(pre2) // '.nc'

                  inquire(file=file_hydro, exist=fexists)
                  IF (fexists) THEN
                     CALL ncio_read_part_serial (file_hydro, dataname, (/j0,i0/), (/j1,i1/), dcache)
                     dcache = transpose(dcache)
                     rdata%blk(iblk,jblk)%val(il0:il1,jl0:jl1) = dcache
                  ENDIF

                  IF ((ieast >= xdsp + 1) .and. (ieast <= xdsp + nxhbox)) THEN 
                     IF (isouth <= ydsp + nyhbox) THEN 
                        exit
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

   END SUBROUTINE catchment_data_read

END MODULE mod_catchment_data
