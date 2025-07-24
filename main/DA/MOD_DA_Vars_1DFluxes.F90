#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_Vars_1DFluxes
!-----------------------------------------------------------------------------
! DESCRIPTION:
!    Process fluxes variables for diagnostic for data assimilation
!
! AUTHOR:
!   Lu Li, 07/2025: Initial version
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Namelist, only: DEF_DA_ENS
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   USE MOD_Vars_1DPFTFluxes
#endif
#ifdef BGC
   USE MOD_BGC_Vars_1DFluxes
#endif
#ifdef CatchLateralFlow
   USE MOD_Catch_Vars_1DFluxes
#endif
#ifdef URBAN_MODEL
   USE MOD_Urban_Vars_1DFluxes
#endif
   IMPLICIT NONE
   SAVE

   ! public functions
   PUBLIC :: allocate_1D_Fluxes_ens
   PUBLIC :: deallocate_1D_Fluxes_ens

   ! define variables
   real(r8), allocatable :: fsena_ens  (:,:) !sensible heat from canopy height to atmosphere [W/m2]
   real(r8), allocatable :: lfevpa_ens (:,:) !latent heat flux from canopy height to atmosphere [W/m2]
   real(r8), allocatable :: fevpa_ens  (:,:) !evapotranspiration from canopy to atmosphere [mm/s]
   real(r8), allocatable :: rsur_ens   (:,:) !surface runoff (mm h2o/s)

!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE allocate_1D_Fluxes_ens()

!-----------------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Vars_Global
      USE MOD_SPMD_Task
      USE MOD_LandPatch
      IMPLICIT NONE

!-----------------------------------------------------------------------------
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            allocate ( fsena_ens  (DEF_DA_ENS, numpatch) )  ; fsena_ens  (:,:) = spval 
            allocate ( lfevpa_ens (DEF_DA_ENS, numpatch) )  ; lfevpa_ens (:,:) = spval 
            allocate ( fevpa_ens  (DEF_DA_ENS, numpatch) )  ; fevpa_ens  (:,:) = spval 
            allocate ( rsur_ens   (DEF_DA_ENS, numpatch) )  ; rsur_ens   (:,:) = spval 
         ENDIF
      ENDIF

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      CALL allocate_1D_PFTFluxes
#endif

#ifdef BGC
      CALL allocate_1D_BGCFluxes
#endif

#ifdef CatchLateralFlow
      CALL allocate_1D_CatchFluxes
#endif

#ifdef URBAN_MODEL
      CALL allocate_1D_UrbanFluxes
#endif

   END SUBROUTINE allocate_1D_Fluxes_ens

!-----------------------------------------------------------------------------

   SUBROUTINE deallocate_1D_Fluxes_ens()

!-----------------------------------------------------------------------------
      USE MOD_SPMD_Task
      USE MOD_LandPatch

!-----------------------------------------------------------------------------
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            deallocate ( fsena_ens  )
            deallocate ( lfevpa_ens )
            deallocate ( fevpa_ens  )
            deallocate ( rsur_ens   )
         ENDIF
      ENDIF

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      CALL deallocate_1D_PFTFluxes
#endif

#ifdef BGC
      CALL deallocate_1D_BGCFluxes
#endif

#ifdef CatchLateralFlow
      CALL deallocate_1D_CatchFluxes
#endif

#ifdef URBAN_MODEL
      CALL deallocate_1D_UrbanFluxes
#endif

   END SUBROUTINE deallocate_1D_Fluxes_ens

!-----------------------------------------------------------------------------
END MODULE MOD_DA_Vars_1DFluxes
#endif