#include <define.h>

MODULE MOD_1D_Fluxes
! -------------------------------
! Created by Yongjiu Dai, 03/2014
! -------------------------------

  USE precision
  IMPLICIT NONE
  SAVE

! -----------------------------------------------------------------
! Fluxes
! -----------------------------------------------------------------
  REAL(r8), allocatable :: taux   (:) !wind stress: E-W [kg/m/s2]
  REAL(r8), allocatable :: tauy   (:) !wind stress: N-S [kg/m/s2]
  REAL(r8), allocatable :: fsena  (:) !sensible heat from canopy height to atmosphere [W/m2]
  REAL(r8), allocatable :: lfevpa (:) !latent heat flux from canopy height to atmosphere [W/m2]
  REAL(r8), allocatable :: fevpa  (:) !evapotranspiration from canopy to atmosphere [mm/s]
  REAL(r8), allocatable :: fsenl  (:) !sensible heat from leaves [W/m2]
  REAL(r8), allocatable :: fevpl  (:) !evaporation+transpiration from leaves [mm/s]
  REAL(r8), allocatable :: etr    (:) !transpiration rate [mm/s]
  REAL(r8), allocatable :: fseng  (:) !sensible heat flux from ground [W/m2]
  REAL(r8), allocatable :: fevpg  (:) !evaporation heat flux from ground [mm/s]
  REAL(r8), allocatable :: fgrnd  (:) !ground heat flux [W/m2]
  REAL(r8), allocatable :: sabvsun(:) !solar absorbed by sunlit vegetation [W/m2]
  REAL(r8), allocatable :: sabvsha(:) !solar absorbed by shaded vegetation [W/m2]
  REAL(r8), allocatable :: sabg   (:) !solar absorbed by ground  [W/m2]
  REAL(r8), allocatable :: sr     (:) !total reflected solar radiation (W/m2)
  REAL(r8), allocatable :: solvd  (:) !incident direct beam vis solar radiation (W/m2)
  REAL(r8), allocatable :: solvi  (:) !incident diffuse beam vis solar radiation (W/m2)
  REAL(r8), allocatable :: solnd  (:) !incident direct beam nir solar radiation (W/m2)
  REAL(r8), allocatable :: solni  (:) !incident diffuse beam nir solar radiation (W/m2)
  REAL(r8), allocatable :: srvd   (:) !reflected direct beam vis solar radiation (W/m2)
  REAL(r8), allocatable :: srvi   (:) !reflected diffuse beam vis solar radiation (W/m2)
  REAL(r8), allocatable :: srnd   (:) !reflected direct beam nir solar radiation (W/m2)
  REAL(r8), allocatable :: srni   (:) !reflected diffuse beam nir solar radiation (W/m2)
  REAL(r8), allocatable :: solvdln(:) !incident direct beam vis solar radiation at local noon (W/m2)
  REAL(r8), allocatable :: solviln(:) !incident diffuse beam vis solar radiation at local noon (W/m2)
  REAL(r8), allocatable :: solndln(:) !incident direct beam nir solar radiation at local noon (W/m2)
  REAL(r8), allocatable :: solniln(:) !incident diffuse beam nir solar radiation at local noon (W/m2)
  REAL(r8), allocatable :: srvdln (:) !reflected direct beam vis solar radiation at local noon (W/m2)
  REAL(r8), allocatable :: srviln (:) !reflected diffuse beam vis solar radiation at local noon (W/m2)
  REAL(r8), allocatable :: srndln (:) !reflected direct beam nir solar radiation at local noon (W/m2)
  REAL(r8), allocatable :: srniln (:) !reflected diffuse beam nir solar radiation at local noon (W/m2)
  REAL(r8), allocatable :: olrg   (:) !outgoing long-wave radiation from ground+canopy [W/m2]
  REAL(r8), allocatable :: rnet   (:) !net radiation by surface [W/m2]
  REAL(r8), allocatable :: xerr   (:) !the error of water banace [mm/s]
  REAL(r8), allocatable :: zerr   (:) !the error of energy balance [W/m2]

  REAL(r8), allocatable :: rsur   (:) !surface runoff (mm h2o/s)
  REAL(r8), allocatable :: rnof   (:) !total runoff (mm h2o/s)
  REAL(r8), allocatable :: qintr  (:) !interception (mm h2o/s)
  REAL(r8), allocatable :: qinfl  (:) !inflitration (mm h2o/s)
  REAL(r8), allocatable :: qdrip  (:) !throughfall (mm h2o/s)
  REAL(r8), allocatable :: assim  (:) !canopy assimilation rate (mol m-2 s-1)
  REAL(r8), allocatable :: respc  (:) !canopy respiration (mol m-2 s-1)

  REAL(r8), allocatable :: qcharge(:) !groundwater recharge [mm/s]

! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: allocate_1D_Fluxes
  PUBLIC :: deallocate_1D_Fluxes

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_1D_Fluxes
  ! --------------------------------------------------------------------
  ! Allocates memory for CLM 1d [numpatch] variables
  ! --------------------------------------------------------------------
     USE precision
     USE GlobalVars
     USE MOD_1D_PFTFluxes
     USE MOD_1D_PCFluxes
     IMPLICIT NONE

     allocate (taux    (numpatch))  !wind stress: E-W [kg/m/s2]
     allocate (tauy    (numpatch))  !wind stress: N-S [kg/m/s2]
     allocate (fsena   (numpatch))  !sensible heat from canopy height to atmosphere [W/m2]
     allocate (lfevpa  (numpatch))  !latent heat flux from canopy height to atmosphere [W/m2]
     allocate (fevpa   (numpatch))  !evapotranspiration from canopy to atmosphere [mm/s]
     allocate (fsenl   (numpatch))  !sensible heat from leaves [W/m2]
     allocate (fevpl   (numpatch))  !evaporation+transpiration from leaves [mm/s]
     allocate (etr     (numpatch))  !transpiration rate [mm/s]
     allocate (fseng   (numpatch))  !sensible heat flux from ground [W/m2]
     allocate (fevpg   (numpatch))  !evaporation heat flux from ground [mm/s]
     allocate (fgrnd   (numpatch))  !ground heat flux [W/m2]
     allocate (sabvsun (numpatch))  !solar absorbed by sunlit vegetation [W/m2]
     allocate (sabvsha (numpatch))  !solar absorbed by shaded vegetation [W/m2]
     allocate (sabg    (numpatch))  !solar absorbed by ground  [W/m2]
     allocate (sr      (numpatch))  !incident direct beam vis solar radiation (W/m2)
     allocate (solvd   (numpatch))  !incident direct beam vis solar radiation (W/m2)
     allocate (solvi   (numpatch))  !incident diffuse beam vis solar radiation (W/m2)
     allocate (solnd   (numpatch))  !incident direct beam nir solar radiation (W/m2)
     allocate (solni   (numpatch))  !incident diffuse beam nir solar radiation (W/m2)
     allocate (srvd    (numpatch))  !reflected direct beam vis solar radiation (W/m2)
     allocate (srvi    (numpatch))  !reflected diffuse beam vis solar radiation (W/m2)
     allocate (srnd    (numpatch))  !reflected direct beam nir solar radiation (W/m2)
     allocate (srni    (numpatch))  !reflected diffuse beam nir solar radiation (W/m2)
     allocate (solvdln (numpatch))  !incident direct beam vis solar radiation at local noon(W/m2)
     allocate (solviln (numpatch))  !incident diffuse beam vis solar radiation at local noon(W/m2)
     allocate (solndln (numpatch))  !incident direct beam nir solar radiation at local noon(W/m2)
     allocate (solniln (numpatch))  !incident diffuse beam nir solar radiation at local noon(W/m2)
     allocate (srvdln  (numpatch))  !reflected direct beam vis solar radiation at local noon(W/m2)
     allocate (srviln  (numpatch))  !reflected diffuse beam vis solar radiation at local noon(W/m2)
     allocate (srndln  (numpatch))  !reflected direct beam nir solar radiation at local noon(W/m2)
     allocate (srniln  (numpatch))  !reflected diffuse beam nir solar radiation at local noon(W/m2)
     allocate (olrg    (numpatch))  !outgoing long-wave radiation from ground+canopy [W/m2]
     allocate (rnet    (numpatch))  !net radiation by surface [W/m2]
     allocate (xerr    (numpatch))  !the error of water banace [mm/s]
     allocate (zerr    (numpatch))  !the error of energy balance [W/m2]

     allocate (rsur    (numpatch))  !surface runoff (mm h2o/s)
     allocate (rnof    (numpatch))  !total runoff (mm h2o/s)
     allocate (qintr   (numpatch))  !interception (mm h2o/s)
     allocate (qinfl   (numpatch))  !inflitration (mm h2o/s)
     allocate (qdrip   (numpatch))  !throughfall (mm h2o/s)
     allocate (assim   (numpatch))  !canopy assimilation rate (mol m-2 s-1)
     allocate (respc   (numpatch))  !canopy respiration (mol m-2 s-1)

     allocate (qcharge (numpatch))  !groundwater recharge [mm/s]

#ifdef PFT_CLASSIFICATION
     CALL allocate_1D_PFTFluxes
#endif

#ifdef PC_CLASSIFICATION
     CALL allocate_1D_PCFluxes
#endif

  END SUBROUTINE allocate_1D_Fluxes

  SUBROUTINE deallocate_1D_Fluxes
  ! --------------------------------------------------------------------
  ! deallocates memory for CLM 1d [numpatch] variables
  ! --------------------------------------------------------------------
     USE MOD_1D_PFTFluxes
     USE MOD_1D_PCFluxes

     deallocate (taux    )  !wind stress: E-W [kg/m/s2]
     deallocate (tauy    )  !wind stress: N-S [kg/m/s2]
     deallocate (fsena   )  !sensible heat from canopy height to atmosphere [W/m2]
     deallocate (lfevpa  )  !latent heat flux from canopy height to atmosphere [W/m2]
     deallocate (fevpa   )  !evapotranspiration from canopy to atmosphere [mm/s]
     deallocate (fsenl   )  !sensible heat from leaves [W/m2]
     deallocate (fevpl   )  !evaporation+transpiration from leaves [mm/s]
     deallocate (etr     )  !transpiration rate [mm/s]
     deallocate (fseng   )  !sensible heat flux from ground [W/m2]
     deallocate (fevpg   )  !evaporation heat flux from ground [mm/s]
     deallocate (fgrnd   )  !ground heat flux [W/m2]
     deallocate (sabvsun )  !solar absorbed by sunlit vegetation [W/m2]
     deallocate (sabvsha )  !solar absorbed by shaded vegetation [W/m2]
     deallocate (sabg    )  !solar absorbed by ground  [W/m2]
     deallocate (sr      )  !incident direct beam vis solar radiation (W/m2)
     deallocate (solvd   )  !incident direct beam vis solar radiation (W/m2)
     deallocate (solvi   )  !incident diffuse beam vis solar radiation (W/m2)
     deallocate (solnd   )  !incident direct beam nir solar radiation (W/m2)
     deallocate (solni   )  !incident diffuse beam nir solar radiation (W/m2)
     deallocate (srvd    )  !reflected direct beam vis solar radiation (W/m2)
     deallocate (srvi    )  !reflected diffuse beam vis solar radiation (W/m2)
     deallocate (srnd    )  !reflected direct beam nir solar radiation (W/m2)
     deallocate (srni    )  !reflected diffuse beam nir solar radiation (W/m2)
     deallocate (solvdln )  !incident direct beam vis solar radiation at local noon(W/m2)
     deallocate (solviln )  !incident diffuse beam vis solar radiation at local noon(W/m2)
     deallocate (solndln )  !incident direct beam nir solar radiation at local noon(W/m2)
     deallocate (solniln )  !incident diffuse beam nir solar radiation at local noon(W/m2)
     deallocate (srvdln  )  !reflected direct beam vis solar radiation at local noon(W/m2)
     deallocate (srviln  )  !reflected diffuse beam vis solar radiation at local noon(W/m2)
     deallocate (srndln  )  !reflected direct beam nir solar radiation at local noon(W/m2)
     deallocate (srniln  )  !reflected diffuse beam nir solar radiation at local noon(W/m2)
     deallocate (olrg    )  !outgoing long-wave radiation from ground+canopy [W/m2]
     deallocate (rnet    )  !net radiation by surface [W/m2]
     deallocate (xerr    )  !the error of water banace [mm/s]
     deallocate (zerr    )  !the error of energy balance [W/m2]

     deallocate (rsur    )  !surface runoff (mm h2o/s)
     deallocate (rnof    )  !total runoff (mm h2o/s)
     deallocate (qintr   )  !interception (mm h2o/s)
     deallocate (qinfl   )  !inflitration (mm h2o/s)
     deallocate (qdrip   )  !throughfall (mm h2o/s)
     deallocate (assim   )  !canopy assimilation rate (mol m-2 s-1)
     deallocate (respc   )  !canopy respiration (mol m-2 s-1)

     deallocate (qcharge )  !groundwater recharge [mm/s]
     
#ifdef PFT_CLASSIFICATION
     CALL deallocate_1D_PFTFluxes
#endif

#ifdef PC_CLASSIFICATION
     CALL deallocate_1D_PCFluxes
#endif

  END SUBROUTINE deallocate_1D_Fluxes

END MODULE MOD_1D_Fluxes
! ---------- EOP ------------
