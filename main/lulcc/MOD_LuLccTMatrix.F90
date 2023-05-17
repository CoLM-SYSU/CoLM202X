#include <define.h>

MODULE MOD_LuLccTMatrix
! -------------------------------
! Created by Hua Yuan, 04/2022
! -------------------------------

  USE precision
  USE GlobalVars
  IMPLICIT NONE
  SAVE
! -----------------------------------------------------------------
  !TODO: need coding below...

! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: allocate_LuLccTMatrix
  PUBLIC :: deallocate_LuLccTMatrix
  PUBLIC :: READ_LuLccTMatrix

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_LuLccTMatrix
  ! --------------------------------------------------------------------
  ! Allocates memory for LuLcc time invariant variables
  ! --------------------------------------------------------------------

     USE precision
     USE GlobalVars
     IMPLICIT NONE
!TODO: need coding below...

  END SUBROUTINE allocate_LuLccTMatrix


  SUBROUTINE READ_LuLccTMatrix

     USE precision
     USE GlobalVars
     USE MOD_TimeInvariants
     USE MOD_PFTimeInvars
     USE MOD_PCTimeInvars
     USE MOD_UrbanTimeInvars
     IMPLICIT NONE
!TODO: need coding below...

  END SUBROUTINE READ_LuLccTMatrix

  SUBROUTINE deallocate_LuLccTMatrix
! --------------------------------------------------
! Deallocates memory for LuLcc time invariant variables
! --------------------------------------------------

!TODO: need coding below...


  END SUBROUTINE deallocate_LuLccTMatrix

END MODULE MOD_LuLccTMatrix
! ---------- EOP ------------
