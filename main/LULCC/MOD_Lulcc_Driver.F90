#include <define.h>

#ifdef LULCC
MODULE MOD_Lulcc_Driver

   USE MOD_Precision
   IMPLICIT NONE
   SAVE

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: LulccDriver

 CONTAINS


   SUBROUTINE LulccDriver (casename,dir_landdata,dir_restart,&
                           idate,greenwich)

!-----------------------------------------------------------------------
!
! !DESCRIPTION:
!  the main subroutine for Land use and land cover change simulation
!
!  Created by Hua Yuan, 04/08/2022
!
! !REVISIONS:
!  07/2023, Wenzong Dong: porting to MPI version.
!  08/2023, Wanyi Lin: add interface for Mass&Energy conserved scheme.
!
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Lulcc_Vars_TimeInvariants
   USE MOD_Lulcc_Vars_TimeVariables
   USE MOD_Lulcc_Initialize
   USE MOD_Vars_TimeVariables
   USE MOD_Lulcc_TransferTrace
   USE MOD_Lulcc_MassEnergyConserve
   USE MOD_Namelist

   IMPLICIT NONE

   character(len=256), intent(in) :: casename      !casename name
   character(len=256), intent(in) :: dir_landdata  !surface data directory
   character(len=256), intent(in) :: dir_restart   !case restart data directory

   logical, intent(in)    :: greenwich   !true: greenwich time, false: local time
   integer, intent(inout) :: idate(3)    !year, julian day, seconds of the starting time

      ! allocate Lulcc memory
      CALL allocate_LulccTimeInvariants
      CALL allocate_LulccTimeVariables

      ! SAVE variables
      CALL SAVE_LulccTimeInvariants
      CALL SAVE_LulccTimeVariables

      ! =============================================================
      ! cold start for Lulcc
      ! =============================================================

      IF (p_is_master) THEN
         print *, ">>> LULCC: initializing..."
      ENDIF

      CALL LulccInitialize (casename,dir_landdata,dir_restart,&
                            idate,greenwich)


      ! =============================================================
      ! 1. Same Type Assignment (SAT) scheme for variable recovery
      ! =============================================================

      IF (DEF_LULCC_SCHEME == 1) THEN
         IF (p_is_master) THEN
            print *, ">>> LULCC: Same Type Assignment (SAT) scheme for variable recovery..."
         ENDIF
         CALL REST_LulccTimeVariables
      ENDIF


      ! =============================================================
      ! 2. Mass and Energy conservation (MEC) scheme for variable revocery
      ! =============================================================

      IF (DEF_LULCC_SCHEME == 2) THEN
         IF (p_is_master) THEN
            print *, ">>> LULCC: Mass&Energy conserve (MEC) for variable recovery..."
         ENDIF
         CALL allocate_LulccTransferTrace()
         CALL REST_LulccTimeVariables
         CALL MAKE_LulccTransferTrace(idate(1))
         CALL LulccMassEnergyConserve()
      ENDIF


      ! deallocate Lulcc memory
      CALL deallocate_LulccTimeInvariants()
      CALL deallocate_LulccTimeVariables()
      IF (DEF_LULCC_SCHEME == 2) THEN
         CALL deallocate_LulccTransferTrace()
      ENDIF

   END SUBROUTINE LulccDriver

END MODULE MOD_Lulcc_Driver
#endif
! ---------- EOP ------------
