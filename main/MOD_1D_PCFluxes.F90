#include <define.h>

MODULE MOD_1D_PCFluxes
! -------------------------------
! Created by Hua Yuan, 08/2019
! -------------------------------

  USE precision
  IMPLICIT NONE
  SAVE

! -----------------------------------------------------------------
! Fluxes
! -----------------------------------------------------------------
  REAL(r8), allocatable :: fsenl_c   (:,:) !sensible heat from leaves [W/m2]
  REAL(r8), allocatable :: fevpl_c   (:,:) !evaporation+transpiration from leaves [mm/s]
  REAL(r8), allocatable :: etr_c     (:,:) !transpiration rate [mm/s]
  REAL(r8), allocatable :: fseng_c   (:,:) !sensible heat flux from ground [W/m2]
  REAL(r8), allocatable :: fevpg_c   (:,:) !evaporation heat flux from ground [mm/s]
  REAL(r8), allocatable :: parsun_c  (:,:) !solar absorbed by sunlit vegetation [W/m2]
  REAL(r8), allocatable :: parsha_c  (:,:) !solar absorbed by shaded vegetation [W/m2]
  REAL(r8), allocatable :: sabvsun_c (:,:) !solar absorbed by sunlit vegetation [W/m2]
  REAL(r8), allocatable :: sabvsha_c (:,:) !solar absorbed by shaded vegetation [W/m2]
  REAL(r8), allocatable :: qintr_c   (:,:) !interception (mm h2o/s)
  REAL(r8), allocatable :: assim_c   (:,:) !canopy assimilation rate (mol m-2 s-1)
  REAL(r8), allocatable :: respc_c   (:,:) !canopy respiration (mol m-2 s-1)

! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: allocate_1D_PCFluxes
  PUBLIC :: deallocate_1D_PCFluxes

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_1D_PCFluxes
  ! --------------------------------------------------------------------
  ! Allocates memory for CLM 1d [numpc] variables
  ! --------------------------------------------------------------------

     USE precision
     USE GlobalVars
     IMPLICIT NONE

     allocate (fsenl_c   (0:N_PFT-1,numpc)) !sensible heat from leaves [W/m2]
     allocate (fevpl_c   (0:N_PFT-1,numpc)) !evaporation+transpiration from leaves [mm/s]
     allocate (etr_c     (0:N_PFT-1,numpc)) !transpiration rate [mm/s]
     allocate (fseng_c   (0:N_PFT-1,numpc)) !sensible heat flux from ground [W/m2]
     allocate (fevpg_c   (0:N_PFT-1,numpc)) !evaporation heat flux from ground [mm/s]
     allocate (parsun_c  (0:N_PFT-1,numpc)) !solar absorbed by sunlit vegetation [W/m2]
     allocate (parsha_c  (0:N_PFT-1,numpc)) !solar absorbed by shaded vegetation [W/m2]
     allocate (sabvsun_c (0:N_PFT-1,numpc)) !solar absorbed by sunlit vegetation [W/m2]
     allocate (sabvsha_c (0:N_PFT-1,numpc)) !solar absorbed by shaded vegetation [W/m2]
     allocate (qintr_c   (0:N_PFT-1,numpc)) !interception (mm h2o/s)
     allocate (assim_c   (0:N_PFT-1,numpc)) !canopy assimilation rate (mol m-2 s-1)
     allocate (respc_c   (0:N_PFT-1,numpc)) !canopy respiration (mol m-2 s-1)

  END SUBROUTINE allocate_1D_PCFluxes

  SUBROUTINE deallocate_1D_PCFluxes
  ! --------------------------------------------------------------------
  ! deallocates memory for CLM 1d [numpc] variables
  ! --------------------------------------------------------------------
     deallocate (fsenl_c   )
     deallocate (fevpl_c   )
     deallocate (etr_c     )
     deallocate (fseng_c   )
     deallocate (fevpg_c   )
     deallocate (parsun_c  )
     deallocate (parsha_c  )
     deallocate (sabvsun_c )
     deallocate (sabvsha_c )
     deallocate (qintr_c   )
     deallocate (assim_c   )
     deallocate (respc_c   )
  END SUBROUTINE deallocate_1D_PCFluxes

END MODULE MOD_1D_PCFluxes
! ---------- EOP ------------
