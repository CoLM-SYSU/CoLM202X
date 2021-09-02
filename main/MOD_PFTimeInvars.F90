#include <define.h> 

MODULE MOD_PFTimeInvars 
! -------------------------------
! Created by Hua Yuan, 08/2019
! -------------------------------

  USE precision
  USE GlobalVars
  IMPLICIT NONE
  SAVE

  ! for PFT_CLASSIFICATION
  INTEGER , allocatable :: pftclass    (:)    !PFT type
  INTEGER , allocatable :: patch_pft_s (:)    !start PFT index of a patch
  INTEGER , allocatable :: patch_pft_e (:)    !end PFT index of a patch
  INTEGER , allocatable :: pft2patch   (:)    !patch index of a PFT
  REAL(r8), allocatable :: pftfrac     (:)    !PFT fractional cover
  REAL(r8), allocatable :: htop_p      (:)    !canopy top height [m]
  REAL(r8), allocatable :: hbot_p      (:)    !canopy bottom height [m]

! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: allocate_PFTimeInvars
  PUBLIC :: deallocate_PFTimeInvars

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_PFTimeInvars
  ! --------------------------------------------------------------------
  ! Allocates memory for CLM PFT 1d [numpft] variables
  ! --------------------------------------------------------------------

     USE precision
     IMPLICIT NONE

     allocate (pftclass      (numpft))   
     allocate (pftfrac       (numpft))   
     allocate (patch_pft_s (numpatch))
     allocate (patch_pft_e (numpatch))
     allocate (pft2patch     (numpft)) 
     allocate (htop_p        (numpft)) 
     allocate (hbot_p        (numpft)) 

  END SUBROUTINE allocate_PFTimeInvars


  SUBROUTINE deallocate_PFTimeInvars
! --------------------------------------------------
! Deallocates memory for CLM PFT 1d [numpft] variables
! --------------------------------------------------

     deallocate (pftclass    )
     deallocate (pftfrac     )
     deallocate (patch_pft_s )
     deallocate (patch_pft_e )
     deallocate (pft2patch   )
     deallocate (htop_p      ) 
     deallocate (hbot_p      ) 

  END SUBROUTINE deallocate_PFTimeInvars

END MODULE MOD_PFTimeInvars
! ---------- EOP ------------
