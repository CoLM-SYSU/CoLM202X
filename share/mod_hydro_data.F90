#include <define.h>

MODULE mod_hydro_data

   IMPLICIT NONE

   INTEGER, parameter :: nxhbox = 6000
   INTEGER, parameter :: nyhbox = 6000
   INTEGER, parameter :: nxhglb = 432000
   INTEGER, parameter :: nyhglb = 216000

CONTAINS

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
   SUBROUTINE hydro_data_read (dir_hydro, dataname, grid, rdata, spv)

      USE spmd_task
      USE mod_block
      USE mod_grid
      USE mod_data_type
      USE netcdf 
      IMPLICIT NONE

      CHARACTER (len=*), intent(in) :: dir_hydro
      CHARACTER (len=*), intent(in) :: dataname
      TYPE (grid_type),  intent(in) :: grid
      TYPE (block_data_int32_2d), intent(inout) :: rdata
      INTEGER, intent(in), optional :: spv

      ! Local variables
      INTEGER :: iblk, jblk, isouth, inorth, iwest, ieast, ibox, jbox
      INTEGER :: xdsp, ydsp, i0, i1, j0, j1, il0, il1, jl0, jl1
      CHARACTER(len=256) :: file_hydro
      CHARACTER(len=3)   :: pre1
      CHARACTER(len=4)   :: pre2
      INTEGER :: ncid, varid
      INTEGER, allocatable :: dcache(:,:)
      LOGICAL :: fexists


      IF (p_is_io) THEN

         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                     
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

                     file_hydro = trim(dir_hydro) // '/' // trim(pre1) // trim(pre2) // '.nc'

                     inquire(file=file_hydro, exist=fexists)
                     IF (fexists) THEN
                        allocate (dcache (j1-j0+1,i1-i0+1))

                        CALL nccheck( nf90_open(trim(file_hydro), NF90_NOWRITE, ncid) )
                        CALL nccheck( nf90_inq_varid(ncid, trim(dataname), varid) )
                        CALL nccheck( nf90_get_var(ncid, varid, dcache, (/j0,i0/), (/j1-j0+1,i1-i0+1/)) )
                        CALL nccheck( nf90_close(ncid) )

                        dcache = transpose(dcache)
                        rdata%blk(iblk,jblk)%val(il0:il1,jl0:jl1) = dcache

                        deallocate (dcache)
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

               ENDIF
            ENDDO
         ENDDO

      ENDIF

   END SUBROUTINE hydro_data_read

END MODULE mod_hydro_data
