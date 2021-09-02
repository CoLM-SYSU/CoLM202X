#include <define.h> 

MODULE MOD_PCTimeInvars 
! -------------------------------
! Created by Hua Yuan, 08/2019
! -------------------------------

  USE precision
  USE GlobalVars
  IMPLICIT NONE
  SAVE

  ! for PC_CLASSIFICATION
  INTEGER , allocatable :: patch2pc(:)    !projection from patch to PC
  INTEGER , allocatable :: pc2patch(:)    !projection from PC to patch
  REAL(r8), allocatable :: pcfrac(:,:)    !PC fractional cover
  REAL(r8), allocatable :: htop_c(:,:)    !canopy top height [m]
  REAL(r8), allocatable :: hbot_c(:,:)    !canopy bottom height [m]   

! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: allocate_PCTimeInvars
  PUBLIC :: deallocate_PCTimeInvars

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_PCTimeInvars
  ! --------------------------------------------------------------------
  ! Allocates memory for CLM Plant Community (PC) [numpc] variables
  ! --------------------------------------------------------------------

     USE precision
     IMPLICIT NONE

     allocate (patch2pc        (numpatch))
     allocate (pc2patch           (numpc))
     allocate (pcfrac   (0:N_PFT-1,numpc))
     allocate (htop_c   (0:N_PFT-1,numpc)) 
     allocate (hbot_c   (0:N_PFT-1,numpc)) 

  END SUBROUTINE allocate_PCTimeInvars


  SUBROUTINE deallocate_PCTimeInvars
! --------------------------------------------------
! Deallocates memory for CLM Plant Community (PC) variables
! --------------------------------------------------

     deallocate (patch2pc )
     deallocate (pc2patch )
     deallocate (pcfrac   )  
     deallocate (htop_c   ) 
     deallocate (hbot_c   ) 

  END SUBROUTINE deallocate_PCTimeInvars

END MODULE MOD_PCTimeInvars
! ---------- EOP ------------
