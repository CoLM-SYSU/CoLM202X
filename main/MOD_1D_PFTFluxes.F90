#include <define.h>

MODULE MOD_1D_PFTFluxes
! -------------------------------
! Created by Hua Yuan, 08/2019
! -------------------------------

  USE precision
  IMPLICIT NONE
  SAVE

! -----------------------------------------------------------------
! Fluxes
! -----------------------------------------------------------------
  REAL(r8), allocatable :: taux_p   (:) !wind stress: E-W [kg/m/s2]
  REAL(r8), allocatable :: tauy_p   (:) !wind stress: N-S [kg/m/s2]
  REAL(r8), allocatable :: fsenl_p  (:) !sensible heat from leaves [W/m2]
  REAL(r8), allocatable :: fevpl_p  (:) !evaporation+transpiration from leaves [mm/s]
  REAL(r8), allocatable :: etr_p    (:) !transpiration rate [mm/s]
  REAL(r8), allocatable :: fseng_p  (:) !sensible heat flux from ground [W/m2]
  REAL(r8), allocatable :: fevpg_p  (:) !evaporation heat flux from ground [mm/s]
  REAL(r8), allocatable :: parsun_p (:) !solar absorbed by sunlit vegetation [W/m2]
  REAL(r8), allocatable :: parsha_p (:) !solar absorbed by shaded vegetation [W/m2]
  REAL(r8), allocatable :: sabvsun_p(:) !solar absorbed by sunlit vegetation [W/m2]
  REAL(r8), allocatable :: sabvsha_p(:) !solar absorbed by shaded vegetation [W/m2]
  REAL(r8), allocatable :: qintr_p  (:) !interception (mm h2o/s)
  REAL(r8), allocatable :: assim_p  (:) !canopy assimilation rate (mol m-2 s-1)
  REAL(r8), allocatable :: respc_p  (:) !canopy respiration (mol m-2 s-1)

! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: allocate_1D_PFTFluxes
  PUBLIC :: deallocate_1D_PFTFluxes

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_1D_PFTFluxes
  ! --------------------------------------------------------------------
  ! Allocates memory for CLM PFT 1d [numpft] variables
  ! --------------------------------------------------------------------

     USE precision
     USE GlobalVars
     IMPLICIT NONE

     allocate (taux_p    (numpft))   !wind stress: E-W [kg/m/s2]
     allocate (tauy_p    (numpft))   !wind stress: N-S [kg/m/s2]
     allocate (fsenl_p   (numpft))   !sensible heat from leaves [W/m2]
     allocate (fevpl_p   (numpft))   !evaporation+transpiration from leaves [mm/s]
     allocate (etr_p     (numpft))   !transpiration rate [mm/s]
     allocate (fseng_p   (numpft))   !sensible heat flux from ground [W/m2]
     allocate (fevpg_p   (numpft))   !evaporation heat flux from ground [mm/s]
     allocate (parsun_p  (numpft))   !solar absorbed by sunlit vegetation [W/m2]
     allocate (parsha_p  (numpft))   !solar absorbed by shaded vegetation [W/m2]
     allocate (sabvsun_p (numpft))   !solar absorbed by sunlit vegetation [W/m2]
     allocate (sabvsha_p (numpft))   !solar absorbed by shaded vegetation [W/m2]
     allocate (qintr_p   (numpft))   !interception (mm h2o/s)
     allocate (assim_p   (numpft))   !canopy assimilation rate (mol m-2 s-1)
     allocate (respc_p   (numpft))   !canopy respiration (mol m-2 s-1)

  END SUBROUTINE allocate_1D_PFTFluxes

  SUBROUTINE deallocate_1D_PFTFluxes
  ! --------------------------------------------------------------------
  ! deallocates memory for CLM PFT 1d [numpft] variables
  ! --------------------------------------------------------------------
     deallocate (taux_p    )
     deallocate (tauy_p    )
     deallocate (fsenl_p   )
     deallocate (fevpl_p   )
     deallocate (etr_p     )
     deallocate (fseng_p   )
     deallocate (fevpg_p   )
     deallocate (parsun_p  )
     deallocate (parsha_p  )
     deallocate (sabvsun_p )
     deallocate (sabvsha_p )
     deallocate (qintr_p   )
     deallocate (assim_p   )
     deallocate (respc_p   )
  END SUBROUTINE deallocate_1D_PFTFluxes

END MODULE MOD_1D_PFTFluxes
! ---------- EOP ------------
