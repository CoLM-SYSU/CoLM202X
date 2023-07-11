#include <define.h>

#ifdef LULCC

MODULE MOD_Lulcc_Driver

!-----------------------------------------------------------------------
   USE MOD_Precision
   IMPLICIT NONE
   SAVE

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: LulccDriver


!-----------------------------------------------------------------------

   CONTAINS

!-----------------------------------------------------------------------



 SUBROUTINE LulccDriver (casename,dir_landdata,dir_restart,&
                         idate,greenwich)

!=======================================================================
! PURPOSE:
!   the main subroutine for Land use and land cover change simulation
!
! Created by Hua Yuan, 04/08/2022
!=======================================================================

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Lulcc_Vars_TimeInvariants
   USE MOD_Lulcc_Vars_TimeVariables
   USE MOD_Lulcc_Initialize
   USE MOD_Vars_TimeVariables
   USE MOD_Lulcc_TMatrix

   IMPLICIT NONE

   CHARACTER(LEN=256), intent(in) :: casename      !casename name
   CHARACTER(LEN=256), intent(in) :: dir_landdata  !surface data directory
   CHARACTER(LEN=256), intent(in) :: dir_restart   !case restart data directory

   LOGICAL, intent(in)    :: greenwich   !true: greenwich time, false: local time
   INTEGER, intent(inout) :: idate(3)    !year, julian day, seconds of the starting time

   ! allocate Lulcc memory
   CALL allocate_LulccTimeInvariants
   CALL allocate_LulccTimeVariables

   ! SAVE variables
   CALL SAVE_LulccTimeInvariants
   CALL SAVE_LulccTimeVariables

   ! cold start for Lulcc
   IF (p_is_master) THEN
      print *, ">>> LULCC: initializing..."
   ENDIF

   CALL LulccInitialize (casename,dir_landdata,dir_restart,&
                         idate,greenwich)

   ! simple method for variable recovery
   IF (p_is_master) THEN
      print *, ">>> LULCC: simple method for variable recovery..."
   ENDIF
   !CALL REST_LulccTimeVariables

   ! conserved method for variable revocery
   IF (p_is_master) THEN
      print *, ">>> LULCC: Mass&Energy conserve for variable recovery..."
   ENDIF
   CALL allocate_LulccTMatrix()
   CALL READ_LulccTMatrix(idate(1))
   CALL LulccEnergyConserve()
   !CALL LulccWaterConserve()

   ! deallocate Lulcc memory
   CALL deallocate_LulccTimeInvariants()
   CALL deallocate_LulccTimeVariables()
   CALL deallocate_LulccTMatrix()

   ! write out state variables
   CALL WRITE_TimeVariables (idate, idate(1), casename, dir_restart)

 END SUBROUTINE LulccDriver

END MODULE MOD_Lulcc_Driver

#endif
