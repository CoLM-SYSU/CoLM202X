#include <define.h>

MODULE MOD_Lulcc_TMatrix
! -------------------------------
! Created by Hua Yuan, 04/2022
! -------------------------------

  USE MOD_Precision
  USE MOD_Vars_Global
  IMPLICIT NONE
  SAVE
! -----------------------------------------------------------------
  !TODO: need coding below...

! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: allocate_LulccTMatrix
  PUBLIC :: deallocate_LulccTMatrix
  PUBLIC :: READ_LulccTMatrix

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_LulccTMatrix
  ! --------------------------------------------------------------------
  ! Allocates memory for Lulcc time invariant variables
  ! --------------------------------------------------------------------

     USE MOD_Precision
     USE MOD_Vars_Global

     IMPLICIT NONE
!TODO: need coding below...

  END SUBROUTINE allocate_LulccTMatrix


  SUBROUTINE READ_LulccTMatrix

     USE MOD_Precision
     USE MOD_Vars_Global
     USE MOD_Vars_TimeInvariants
     USE MOD_Vars_PFTimeInvariants
     USE MOD_Vars_PCTimeInvariants
     USE MOD_Urban_Vars_TimeInvariants
     IMPLICIT NONE
!TODO: need coding below...

  END SUBROUTINE READ_LulccTMatrix

  SUBROUTINE deallocate_LulccTMatrix
! --------------------------------------------------
! Deallocates memory for Lulcc time invariant variables
! --------------------------------------------------

!TODO: need coding below...


  END SUBROUTINE deallocate_LulccTMatrix

END MODULE MOD_LulccTMatrix
! ---------- EOP ------------
