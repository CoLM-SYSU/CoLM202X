#include <define.h>

SUBROUTINE LuLccDRIVER (casename,dir_landdata,dir_restart,&
                        idate,greenwich)

!=======================================================================
! PURPOSE:
!   the main subroutine for Land use and land cover change simulation
!
! Created by Hua Yuan, 04/08/2022
!=======================================================================

   USE precision
   USE spmd_task
   USE MOD_LuLccTimeInvars
   USE MOD_LuLccTimeVars
   USE MOD_TimeVariables
   ! USE MOD_LuLccTMatrix

   IMPLICIT NONE

   CHARACTER(LEN=256), intent(in) :: casename      !casename name
   CHARACTER(LEN=256), intent(in) :: dir_landdata  !surface data directory
   CHARACTER(LEN=256), intent(in) :: dir_restart   !case restart data directory

   LOGICAL, intent(in)    :: greenwich   !true: greenwich time, false: local time
   INTEGER, intent(inout) :: idate(3)    !year, julian day, seconds of the starting time

   ! allocate LuLcc memory
   CALL allocate_LuLccTimeInvars
   CALL allocate_LuLccTimeVars

   ! SAVE variables
   CALL SAVE_LuLccTimeInvars
   CALL SAVE_LuLccTimeVars

   ! cold start for LuLcc
   IF (p_is_master) THEN
      print *, ">>> LULCC: initializing..."
   ENDIF

   CALL LuLccInitialize (casename,dir_landdata,dir_restart,&
                         idate,greenwich)

   ! simple method for variable recovery
   IF (p_is_master) THEN
      print *, ">>> LULCC: simple method for variable recovery..."
   ENDIF
   CALL REST_LuLccTimeVars

   ! conserved method for variable revocery
   !print *, ">>> LULCC: Mass&Energy conserve for variable recovery..."
   !CALL READ_LuLccTMatrix()
   !CALL LuLccEnergyConserve()
   !CALL LuLccWaterConserve()

   ! deallocate LuLcc memory
   CALL deallocate_LuLccTimeInvars()
   CALL deallocate_LuLccTimeVars()

   ! write out state variables
   CALL WRITE_TimeVariables (idate, casename, dir_restart)
END SUBROUTINE LuLccDRIVER
