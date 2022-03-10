#include <define.h>

MODULE MOD_1D_cama_Fluxes
! -------------------------------
! Created by Yongjiu Dai, 03/2014
! -------------------------------

  USE precision
  IMPLICIT NONE
  SAVE

  !REAL(r8), allocatable :: rnof   (:) !total runoff (mm h2o/s)
  TYPE history_var_cama_type
  LOGICAL :: rivout       = .false.
  LOGICAL :: rivsto       = .false.
  LOGICAL :: rivdph       = .false.
  LOGICAL :: rivvel       = .false.
  LOGICAL :: fldout       = .false.
  LOGICAL :: fldsto       = .false.
  LOGICAL :: flddph       = .false.
  LOGICAL :: fldfrc       = .false.
  LOGICAL :: fldare       = .false.
  LOGICAL :: sfcelv       = .false.
  LOGICAL :: totout       = .false.
  LOGICAL :: outflw       = .false.
  LOGICAL :: totsto       = .false.
  LOGICAL :: storge       = .false.
  LOGICAL :: pthout       = .false.
  LOGICAL :: maxflw       = .false.
  LOGICAL :: maxdph       = .false.
  LOGICAL :: maxsto       = .false.
  LOGICAL :: gwsto        = .false.
  LOGICAL :: gdwsto       = .false.
  LOGICAL :: gwout        = .false.
  LOGICAL :: gdwrtn       = .false.
  LOGICAL :: runoff       = .false.
  LOGICAL :: runoffsub    = .false.
  LOGICAL :: rofsfc       = .false.
  LOGICAL :: rofsub       = .false.
  LOGICAL :: damsto       = .false.
  LOGICAL :: daminf       = .false.
END TYPE history_var_cama_type
TYPE (history_var_cama_type) :: DEF_hist_cama_vars

! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: allocate_1D_cama_Fluxes
  PUBLIC :: deallocate_1D_cama_Fluxes

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_1D_cama_Fluxes
  ! --------------------------------------------------------------------
  ! Allocates memory for CLM 1d [numpatch] variables
  ! --------------------------------------------------------------------
     USE precision
     USE GlobalVars
!#ifdef PFT_CLASSIFICATION
!     USE MOD_1D_PFTFluxes
!#endif
!#ifdef PC_CLASSIFICATION
!     USE MOD_1D_PCFluxes
!#endif
     USE spmd_task
     USE mod_landpatch
     IMPLICIT NONE

      if (p_is_worker) then

         if (numpatch > 0) then

     !       allocate ( rnof   (numpatch) )  ! total runoff (mm h2o/s)

         end if
      end if

!#ifdef PFT_CLASSIFICATION
 !     CALL allocate_1D_PFTFluxes
!#endif

!#ifdef PC_CLASSIFICATION
!      CALL allocate_1D_PCFluxes
!#endif

   END SUBROUTINE allocate_1D_cama_Fluxes

   SUBROUTINE deallocate_1D_cama_Fluxes ()
  ! --------------------------------------------------------------------
  ! deallocates memory for CLM 1d [numpatch] variables
  ! --------------------------------------------------------------------
!#ifdef PFT_CLASSIFICATION
!     USE MOD_1D_PFTFluxes
!#endif
!#ifdef PC_CLASSIFICATION
!     USE MOD_1D_PCFluxes
!#endif
     USE spmd_task
     USE mod_landpatch

     if (p_is_worker) then
        
        if (numpatch > 0) then

    
  !         deallocate ( rnof    )  ! total runoff (mm h2o/s)
 

        end if
     end if
   END SUBROUTINE deallocate_1D_cama_Fluxes

END MODULE MOD_1D_cama_Fluxes
! ---------- EOP ------------

