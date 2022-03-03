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

! bgc variables
  REAL(r8), allocatable :: leafc_xfer_to_leafc_p                (:)
  REAL(r8), allocatable :: frootc_xfer_to_frootc_p              (:)
  REAL(r8), allocatable :: livestemc_xfer_to_livestemc_p        (:)
  REAL(r8), allocatable :: deadstemc_xfer_to_deadstemc_p        (:)
  REAL(r8), allocatable :: livecrootc_xfer_to_livecrootc_p      (:)
  REAL(r8), allocatable :: deadcrootc_xfer_to_deadcrootc_p      (:)
  REAL(r8), allocatable :: grainc_xfer_to_grainc_p              (:)

  REAL(r8), allocatable :: leafc_storage_to_xfer_p              (:)
  REAL(r8), allocatable :: frootc_storage_to_xfer_p             (:)
  REAL(r8), allocatable :: livestemc_storage_to_xfer_p          (:)
  REAL(r8), allocatable :: deadstemc_storage_to_xfer_p          (:)
  REAL(r8), allocatable :: livecrootc_storage_to_xfer_p         (:)
  REAL(r8), allocatable :: deadcrootc_storage_to_xfer_p         (:)
  REAL(r8), allocatable :: grainc_storage_to_xfer_p             (:)
  REAL(r8), allocatable :: gresp_storage_to_xfer_p              (:)

  REAL(r8), allocatable :: leafc_to_litter_p                    (:)
  REAL(r8), allocatable :: frootc_to_litter_p                   (:)
  REAL(r8), allocatable :: grainc_to_food_p                     (:)
  REAL(r8), allocatable :: grainc_to_seed_p                     (:)
  REAL(r8), allocatable :: crop_seedc_to_leaf_p                 (:)
  REAL(r8), allocatable :: livestemc_to_litter_p                (:)
  REAL(r8), allocatable :: livestemc_to_deadstemc_p             (:)
  REAL(r8), allocatable :: livecrootc_to_deadcrootc_p           (:)

  REAL(r8), allocatable :: m_leafc_to_litter_p                  (:)
  REAL(r8), allocatable :: m_frootc_to_litter_p                 (:)
  REAL(r8), allocatable :: m_livestemc_to_litter_p              (:)
  REAL(r8), allocatable :: m_deadstemc_to_litter_p              (:)
  REAL(r8), allocatable :: m_livecrootc_to_litter_p             (:)
  REAL(r8), allocatable :: m_deadcrootc_to_litter_p             (:)

  REAL(r8), allocatable :: m_leafc_storage_to_litter_p          (:)
  REAL(r8), allocatable :: m_frootc_storage_to_litter_p         (:)
  REAL(r8), allocatable :: m_livestemc_storage_to_litter_p      (:)
  REAL(r8), allocatable :: m_deadstemc_storage_to_litter_p      (:)
  REAL(r8), allocatable :: m_livecrootc_storage_to_litter_p     (:)
  REAL(r8), allocatable :: m_deadcrootc_storage_to_litter_p     (:)
  REAL(r8), allocatable :: m_gresp_storage_to_litter_p          (:)

  REAL(r8), allocatable :: m_leafc_xfer_to_litter_p             (:)
  REAL(r8), allocatable :: m_frootc_xfer_to_litter_p            (:)
  REAL(r8), allocatable :: m_livestemc_xfer_to_litter_p         (:)
  REAL(r8), allocatable :: m_deadstemc_xfer_to_litter_p         (:)
  REAL(r8), allocatable :: m_livecrootc_xfer_to_litter_p        (:)
  REAL(r8), allocatable :: m_deadcrootc_xfer_to_litter_p        (:)
  REAL(r8), allocatable :: m_gresp_xfer_to_litter_p             (:)

  REAL(r8), allocatable :: m_leafc_to_fire_p                    (:)
  REAL(r8), allocatable :: m_frootc_to_fire_p                   (:)
  REAL(r8), allocatable :: m_livestemc_to_fire_p                (:)
  REAL(r8), allocatable :: m_deadstemc_to_fire_p                (:)
  REAL(r8), allocatable :: m_livecrootc_to_fire_p               (:)
  REAL(r8), allocatable :: m_deadcrootc_to_fire_p               (:)

  REAL(r8), allocatable :: m_leafc_storage_to_fire_p            (:)
  REAL(r8), allocatable :: m_frootc_storage_to_fire_p           (:)
  REAL(r8), allocatable :: m_livestemc_storage_to_fire_p        (:)
  REAL(r8), allocatable :: m_deadstemc_storage_to_fire_p        (:)
  REAL(r8), allocatable :: m_livecrootc_storage_to_fire_p       (:)
  REAL(r8), allocatable :: m_deadcrootc_storage_to_fire_p       (:)
  REAL(r8), allocatable :: m_gresp_storage_to_fire_p            (:)

  REAL(r8), allocatable :: m_leafc_xfer_to_fire_p               (:)
  REAL(r8), allocatable :: m_frootc_xfer_to_fire_p              (:)
  REAL(r8), allocatable :: m_livestemc_xfer_to_fire_p           (:)
  REAL(r8), allocatable :: m_deadstemc_xfer_to_fire_p           (:)
  REAL(r8), allocatable :: m_livecrootc_xfer_to_fire_p          (:)
  REAL(r8), allocatable :: m_deadcrootc_xfer_to_fire_p          (:)
  REAL(r8), allocatable :: m_gresp_xfer_to_fire_p               (:)

  REAL(r8), allocatable :: m_livestemc_to_deadstemc_fire_p      (:)
  REAL(r8), allocatable :: m_livecrootc_to_deadcrootc_fire_p    (:)

  REAL(r8), allocatable :: m_leafc_to_litter_fire_p             (:)
  REAL(r8), allocatable :: m_frootc_to_litter_fire_p            (:)
  REAL(r8), allocatable :: m_livestemc_to_litter_fire_p         (:)
  REAL(r8), allocatable :: m_deadstemc_to_litter_fire_p         (:)
  REAL(r8), allocatable :: m_livecrootc_to_litter_fire_p        (:)
  REAL(r8), allocatable :: m_deadcrootc_to_litter_fire_p        (:)

  REAL(r8), allocatable :: m_leafc_storage_to_litter_fire_p     (:)
  REAL(r8), allocatable :: m_frootc_storage_to_litter_fire_p    (:)
  REAL(r8), allocatable :: m_livestemc_storage_to_litter_fire_p (:)
  REAL(r8), allocatable :: m_deadstemc_storage_to_litter_fire_p (:)
  REAL(r8), allocatable :: m_livecrootc_storage_to_litter_fire_p(:)
  REAL(r8), allocatable :: m_deadcrootc_storage_to_litter_fire_p(:)
  REAL(r8), allocatable :: m_gresp_storage_to_litter_fire_p     (:)

  REAL(r8), allocatable :: m_leafc_xfer_to_litter_fire_p        (:)
  REAL(r8), allocatable :: m_frootc_xfer_to_litter_fire_p       (:)
  REAL(r8), allocatable :: m_livestemc_xfer_to_litter_fire_p    (:)
  REAL(r8), allocatable :: m_deadstemc_xfer_to_litter_fire_p    (:)
  REAL(r8), allocatable :: m_livecrootc_xfer_to_litter_fire_p   (:)
  REAL(r8), allocatable :: m_deadcrootc_xfer_to_litter_fire_p   (:)
  REAL(r8), allocatable :: m_gresp_xfer_to_litter_fire_p        (:)

  REAL(r8), allocatable :: cpool_to_xsmrpool_p          (:)
  REAL(r8), allocatable :: cpool_to_gresp_storage_p     (:)
  REAL(r8), allocatable :: cpool_to_leafc_p             (:)
  REAL(r8), allocatable :: cpool_to_leafc_storage_p     (:)
  REAL(r8), allocatable :: cpool_to_frootc_p            (:)
  REAL(r8), allocatable :: cpool_to_frootc_storage_p    (:)
  REAL(r8), allocatable :: cpool_to_livestemc_p         (:)
  REAL(r8), allocatable :: cpool_to_livestemc_storage_p (:)
  REAL(r8), allocatable :: cpool_to_deadstemc_p         (:)
  REAL(r8), allocatable :: cpool_to_deadstemc_storage_p (:)
  REAL(r8), allocatable :: cpool_to_livecrootc_p        (:)
  REAL(r8), allocatable :: cpool_to_livecrootc_storage_p(:)
  REAL(r8), allocatable :: cpool_to_deadcrootc_p        (:)
  REAL(r8), allocatable :: cpool_to_deadcrootc_storage_p(:) 
  REAL(r8), allocatable :: cpool_to_grainc_p            (:)
  REAL(r8), allocatable :: cpool_to_grainc_storage_p    (:)

  REAL(r8), allocatable :: leaf_xsmr_p                  (:)
  REAL(r8), allocatable :: froot_xsmr_p                 (:)
  REAL(r8), allocatable :: livestem_xsmr_p              (:)
  REAL(r8), allocatable :: livecroot_xsmr_p             (:)
  REAL(r8), allocatable :: grain_xsmr_p                 (:)

  REAL(r8), allocatable :: cpool_leaf_gr_p                      (:)
  REAL(r8), allocatable :: cpool_froot_gr_p                     (:)
  REAL(r8), allocatable :: cpool_livestem_gr_p                  (:)
  REAL(r8), allocatable :: cpool_deadstem_gr_p                  (:)
  REAL(r8), allocatable :: cpool_livecroot_gr_p                 (:)
  REAL(r8), allocatable :: cpool_deadcroot_gr_p                 (:)
  REAL(r8), allocatable :: cpool_grain_gr_p                     (:)

  REAL(r8), allocatable :: cpool_leaf_storage_gr_p              (:)
  REAL(r8), allocatable :: cpool_froot_storage_gr_p             (:)
  REAL(r8), allocatable :: cpool_livestem_storage_gr_p          (:)
  REAL(r8), allocatable :: cpool_deadstem_storage_gr_p          (:)
  REAL(r8), allocatable :: cpool_livecroot_storage_gr_p         (:)
  REAL(r8), allocatable :: cpool_deadcroot_storage_gr_p         (:)
  REAL(r8), allocatable :: cpool_grain_storage_gr_p             (:)

  REAL(r8), allocatable :: transfer_leaf_gr_p                   (:)
  REAL(r8), allocatable :: transfer_froot_gr_p                  (:)
  REAL(r8), allocatable :: transfer_livestem_gr_p               (:)
  REAL(r8), allocatable :: transfer_deadstem_gr_p               (:)
  REAL(r8), allocatable :: transfer_livecroot_gr_p              (:)
  REAL(r8), allocatable :: transfer_deadcroot_gr_p              (:)
  REAL(r8), allocatable :: transfer_grain_gr_p                  (:)

  REAL(r8), allocatable :: xsmrpool_to_atm_p                    (:)

  REAL(r8), allocatable :: cropprod1c_loss_p                    (:)

  REAL(r8), allocatable :: plant_ndemand_p              (:)

  REAL(r8), allocatable :: leafn_xfer_to_leafn_p                (:)
  REAL(r8), allocatable :: frootn_xfer_to_frootn_p              (:)
  REAL(r8), allocatable :: livestemn_xfer_to_livestemn_p        (:)
  REAL(r8), allocatable :: deadstemn_xfer_to_deadstemn_p        (:)
  REAL(r8), allocatable :: livecrootn_xfer_to_livecrootn_p      (:)
  REAL(r8), allocatable :: deadcrootn_xfer_to_deadcrootn_p      (:)
  REAL(r8), allocatable :: grainn_xfer_to_grainn_p              (:)

  REAL(r8), allocatable :: leafn_storage_to_xfer_p              (:)
  REAL(r8), allocatable :: frootn_storage_to_xfer_p             (:)
  REAL(r8), allocatable :: livestemn_storage_to_xfer_p          (:)
  REAL(r8), allocatable :: deadstemn_storage_to_xfer_p          (:)
  REAL(r8), allocatable :: livecrootn_storage_to_xfer_p         (:)
  REAL(r8), allocatable :: deadcrootn_storage_to_xfer_p         (:)
  REAL(r8), allocatable :: grainn_storage_to_xfer_p             (:)

  REAL(r8), allocatable :: leafn_to_litter_p                    (:)
  REAL(r8), allocatable :: frootn_to_litter_p                   (:)
  REAL(r8), allocatable :: grainn_to_food_p                     (:)
  REAL(r8), allocatable :: grainn_to_seed_p                     (:)
  REAL(r8), allocatable :: crop_seedn_to_leaf_p                 (:)
  REAL(r8), allocatable :: livestemn_to_litter_p                (:)
  REAL(r8), allocatable :: livestemn_to_deadstemn_p             (:)
  REAL(r8), allocatable :: livecrootn_to_deadcrootn_p           (:)

  REAL(r8), allocatable :: leafn_to_retransn_p                  (:)
  REAL(r8), allocatable :: frootn_to_retransn_p                 (:)
  REAL(r8), allocatable :: livestemn_to_retransn_p              (:)
  REAL(r8), allocatable :: livecrootn_to_retransn_p             (:)
  REAL(r8), allocatable :: retransn_to_npool_p                  (:)
  REAL(r8), allocatable :: free_retransn_to_npool_p             (:)

  REAL(r8), allocatable :: m_leafn_to_litter_p                  (:)
  REAL(r8), allocatable :: m_frootn_to_litter_p                 (:)
  REAL(r8), allocatable :: m_livestemn_to_litter_p              (:)
  REAL(r8), allocatable :: m_deadstemn_to_litter_p              (:)
  REAL(r8), allocatable :: m_livecrootn_to_litter_p             (:)
  REAL(r8), allocatable :: m_deadcrootn_to_litter_p             (:)
  REAL(r8), allocatable :: m_retransn_to_litter_p               (:)

  REAL(r8), allocatable :: m_leafn_storage_to_litter_p          (:)
  REAL(r8), allocatable :: m_frootn_storage_to_litter_p         (:)
  REAL(r8), allocatable :: m_livestemn_storage_to_litter_p      (:)
  REAL(r8), allocatable :: m_deadstemn_storage_to_litter_p      (:)
  REAL(r8), allocatable :: m_livecrootn_storage_to_litter_p     (:)
  REAL(r8), allocatable :: m_deadcrootn_storage_to_litter_p     (:)

  REAL(r8), allocatable :: m_leafn_xfer_to_litter_p             (:)
  REAL(r8), allocatable :: m_frootn_xfer_to_litter_p            (:)
  REAL(r8), allocatable :: m_livestemn_xfer_to_litter_p         (:)
  REAL(r8), allocatable :: m_deadstemn_xfer_to_litter_p         (:)
  REAL(r8), allocatable :: m_livecrootn_xfer_to_litter_p        (:)
  REAL(r8), allocatable :: m_deadcrootn_xfer_to_litter_p        (:)

  REAL(r8), allocatable :: m_leafn_to_fire_p                    (:)
  REAL(r8), allocatable :: m_frootn_to_fire_p                   (:)
  REAL(r8), allocatable :: m_livestemn_to_fire_p                (:)
  REAL(r8), allocatable :: m_deadstemn_to_fire_p                (:)
  REAL(r8), allocatable :: m_livecrootn_to_fire_p               (:)
  REAL(r8), allocatable :: m_deadcrootn_to_fire_p               (:)

  REAL(r8), allocatable :: m_leafn_storage_to_fire_p            (:)
  REAL(r8), allocatable :: m_frootn_storage_to_fire_p           (:)
  REAL(r8), allocatable :: m_livestemn_storage_to_fire_p        (:)
  REAL(r8), allocatable :: m_deadstemn_storage_to_fire_p        (:)
  REAL(r8), allocatable :: m_livecrootn_storage_to_fire_p       (:)
  REAL(r8), allocatable :: m_deadcrootn_storage_to_fire_p       (:)

  REAL(r8), allocatable :: m_leafn_xfer_to_fire_p               (:)
  REAL(r8), allocatable :: m_frootn_xfer_to_fire_p              (:)
  REAL(r8), allocatable :: m_livestemn_xfer_to_fire_p           (:)
  REAL(r8), allocatable :: m_deadstemn_xfer_to_fire_p           (:)
  REAL(r8), allocatable :: m_livecrootn_xfer_to_fire_p          (:)
  REAL(r8), allocatable :: m_deadcrootn_xfer_to_fire_p          (:)

  REAL(r8), allocatable :: m_livestemn_to_deadstemn_fire_p      (:)
  REAL(r8), allocatable :: m_livecrootn_to_deadcrootn_fire_p    (:)

  REAL(r8), allocatable :: m_retransn_to_fire_p                 (:)

  REAL(r8), allocatable :: m_leafn_to_litter_fire_p             (:)
  REAL(r8), allocatable :: m_frootn_to_litter_fire_p            (:)
  REAL(r8), allocatable :: m_livestemn_to_litter_fire_p         (:)
  REAL(r8), allocatable :: m_deadstemn_to_litter_fire_p         (:)
  REAL(r8), allocatable :: m_livecrootn_to_litter_fire_p        (:)
  REAL(r8), allocatable :: m_deadcrootn_to_litter_fire_p        (:)

  REAL(r8), allocatable :: m_leafn_storage_to_litter_fire_p     (:)
  REAL(r8), allocatable :: m_frootn_storage_to_litter_fire_p    (:)
  REAL(r8), allocatable :: m_livestemn_storage_to_litter_fire_p (:)
  REAL(r8), allocatable :: m_deadstemn_storage_to_litter_fire_p (:)
  REAL(r8), allocatable :: m_livecrootn_storage_to_litter_fire_p(:)
  REAL(r8), allocatable :: m_deadcrootn_storage_to_litter_fire_p(:)

  REAL(r8), allocatable :: m_leafn_xfer_to_litter_fire_p        (:)
  REAL(r8), allocatable :: m_frootn_xfer_to_litter_fire_p       (:)
  REAL(r8), allocatable :: m_livestemn_xfer_to_litter_fire_p    (:)
  REAL(r8), allocatable :: m_deadstemn_xfer_to_litter_fire_p    (:)
  REAL(r8), allocatable :: m_livecrootn_xfer_to_litter_fire_p   (:)
  REAL(r8), allocatable :: m_deadcrootn_xfer_to_litter_fire_p   (:)

  REAL(r8), allocatable :: m_retransn_to_litter_fire_p          (:)

  REAL(r8), allocatable :: npool_to_leafn_p                     (:)
  REAL(r8), allocatable :: npool_to_leafn_storage_p             (:)
  REAL(r8), allocatable :: npool_to_frootn_p                    (:)
  REAL(r8), allocatable :: npool_to_frootn_storage_p            (:)
  REAL(r8), allocatable :: npool_to_livestemn_p                 (:)
  REAL(r8), allocatable :: npool_to_livestemn_storage_p         (:)
  REAL(r8), allocatable :: npool_to_deadstemn_p                 (:)
  REAL(r8), allocatable :: npool_to_deadstemn_storage_p         (:)
  REAL(r8), allocatable :: npool_to_livecrootn_p                (:)
  REAL(r8), allocatable :: npool_to_livecrootn_storage_p        (:)
  REAL(r8), allocatable :: npool_to_deadcrootn_p                (:)
  REAL(r8), allocatable :: npool_to_deadcrootn_storage_p        (:)
  REAL(r8), allocatable :: npool_to_grainn_p                    (:)
  REAL(r8), allocatable :: npool_to_grainn_storage_p            (:)

  REAL(r8), allocatable :: respcsun_p     (:) !sunlit leaf respiration
  REAL(r8), allocatable :: respcsha_p     (:) !shaded leaf respiration
  REAL(r8), allocatable :: leaf_mr_p      (:) !leaf maintenance respiration
  REAL(r8), allocatable :: froot_mr_p     (:) !fine root maintenance respiration
  REAL(r8), allocatable :: livestem_mr_p  (:) !live stem maintenance respiration
  REAL(r8), allocatable :: livecroot_mr_p (:) !live coarse root maintenance respiration
  REAL(r8), allocatable :: grain_mr_p     (:) !grain maintenance respiration

  REAL(r8), allocatable :: soil_change_p  (:)

  REAL(r8), allocatable :: psn_to_cpool_p               (:)
  REAL(r8), allocatable :: gpp_p                        (:)
  REAL(r8), allocatable :: availc_p                     (:)
  REAL(r8), allocatable :: avail_retransn_p             (:)
  REAL(r8), allocatable :: xsmrpool_recover_p           (:)
  REAL(r8), allocatable :: excess_cflux_p               (:)
  REAL(r8), allocatable :: sminn_to_npool_p             (:)

  REAL(r8), allocatable :: plant_calloc_p               (:)
  REAL(r8), allocatable :: plant_nalloc_p               (:)
  REAL(r8), allocatable :: leaf_curmr_p                 (:)
  REAL(r8), allocatable :: froot_curmr_p                (:)
  REAL(r8), allocatable :: livestem_curmr_p             (:)
  REAL(r8), allocatable :: livecroot_curmr_p            (:)
  REAL(r8), allocatable :: grain_curmr_p                (:)

  REAL(r8), allocatable :: fire_closs_p                 (:)
  REAL(r8), allocatable :: fire_nloss_p                 (:)
  REAL(r8), allocatable :: wood_harvestc_p              (:)
  REAL(r8), allocatable :: wood_harvestn_p              (:)
  REAL(r8), allocatable :: grainc_to_cropprodc_p        (:)
  REAL(r8), allocatable :: grainn_to_cropprodn_p        (:)
  REAL(r8), allocatable :: hrv_xsmrpool_to_atm_p        (:)
  REAL(r8), allocatable :: fert_p                       (:)
  REAL(r8), allocatable :: soyfixn_p                    (:)
!--------


 !   soyfixn                           , &


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

! bgc variables
     allocate (leafc_xfer_to_leafc_p                (numpft))
     allocate (frootc_xfer_to_frootc_p              (numpft))
     allocate (livestemc_xfer_to_livestemc_p        (numpft))
     allocate (deadstemc_xfer_to_deadstemc_p        (numpft))
     allocate (livecrootc_xfer_to_livecrootc_p      (numpft))
     allocate (deadcrootc_xfer_to_deadcrootc_p      (numpft))
     allocate (grainc_xfer_to_grainc_p              (numpft))

     allocate (leafc_storage_to_xfer_p              (numpft))
     allocate (frootc_storage_to_xfer_p             (numpft))
     allocate (livestemc_storage_to_xfer_p          (numpft))
     allocate (deadstemc_storage_to_xfer_p          (numpft))
     allocate (livecrootc_storage_to_xfer_p         (numpft))
     allocate (deadcrootc_storage_to_xfer_p         (numpft))
     allocate (grainc_storage_to_xfer_p             (numpft))
     allocate (gresp_storage_to_xfer_p              (numpft))

     allocate (leafc_to_litter_p                    (numpft))
     allocate (frootc_to_litter_p                   (numpft))
     allocate (grainc_to_food_p                     (numpft))
     allocate (grainc_to_seed_p                     (numpft))
     allocate (crop_seedc_to_leaf_p                 (numpft))
     allocate (livestemc_to_litter_p                (numpft))
     allocate (livestemc_to_deadstemc_p             (numpft))
     allocate (livecrootc_to_deadcrootc_p           (numpft))

     allocate (m_leafc_to_litter_p                  (numpft))
     allocate (m_frootc_to_litter_p                 (numpft))
     allocate (m_livestemc_to_litter_p              (numpft))
     allocate (m_deadstemc_to_litter_p              (numpft))
     allocate (m_livecrootc_to_litter_p             (numpft))
     allocate (m_deadcrootc_to_litter_p             (numpft))

     allocate (m_leafc_storage_to_litter_p          (numpft))
     allocate (m_frootc_storage_to_litter_p         (numpft))
     allocate (m_livestemc_storage_to_litter_p      (numpft))
     allocate (m_deadstemc_storage_to_litter_p      (numpft))
     allocate (m_livecrootc_storage_to_litter_p     (numpft))
     allocate (m_deadcrootc_storage_to_litter_p     (numpft))
     allocate (m_gresp_storage_to_litter_p          (numpft))

     allocate (m_leafc_xfer_to_litter_p             (numpft))
     allocate (m_frootc_xfer_to_litter_p            (numpft))
     allocate (m_livestemc_xfer_to_litter_p         (numpft))
     allocate (m_deadstemc_xfer_to_litter_p         (numpft))
     allocate (m_livecrootc_xfer_to_litter_p        (numpft))
     allocate (m_deadcrootc_xfer_to_litter_p        (numpft))
     allocate (m_gresp_xfer_to_litter_p             (numpft))

     allocate (m_leafc_to_fire_p                    (numpft))
     allocate (m_frootc_to_fire_p                   (numpft))
     allocate (m_livestemc_to_fire_p                (numpft))
     allocate (m_deadstemc_to_fire_p                (numpft))
     allocate (m_livecrootc_to_fire_p               (numpft))
     allocate (m_deadcrootc_to_fire_p               (numpft))

     allocate (m_leafc_storage_to_fire_p            (numpft))
     allocate (m_frootc_storage_to_fire_p           (numpft))
     allocate (m_livestemc_storage_to_fire_p        (numpft))
     allocate (m_deadstemc_storage_to_fire_p        (numpft))
     allocate (m_livecrootc_storage_to_fire_p       (numpft))
     allocate (m_deadcrootc_storage_to_fire_p       (numpft))
     allocate (m_gresp_storage_to_fire_p            (numpft))

     allocate (m_leafc_xfer_to_fire_p               (numpft))
     allocate (m_frootc_xfer_to_fire_p              (numpft))
     allocate (m_livestemc_xfer_to_fire_p           (numpft))
     allocate (m_deadstemc_xfer_to_fire_p           (numpft))
     allocate (m_livecrootc_xfer_to_fire_p          (numpft))
     allocate (m_deadcrootc_xfer_to_fire_p          (numpft))
     allocate (m_gresp_xfer_to_fire_p               (numpft))

     allocate (m_livestemc_to_deadstemc_fire_p      (numpft))
     allocate (m_livecrootc_to_deadcrootc_fire_p    (numpft))

     allocate (m_leafc_to_litter_fire_p             (numpft))
     allocate (m_frootc_to_litter_fire_p            (numpft))
     allocate (m_livestemc_to_litter_fire_p         (numpft))
     allocate (m_deadstemc_to_litter_fire_p         (numpft))
     allocate (m_livecrootc_to_litter_fire_p        (numpft))
     allocate (m_deadcrootc_to_litter_fire_p        (numpft))

     allocate (m_leafc_storage_to_litter_fire_p     (numpft))
     allocate (m_frootc_storage_to_litter_fire_p    (numpft))
     allocate (m_livestemc_storage_to_litter_fire_p (numpft))
     allocate (m_deadstemc_storage_to_litter_fire_p (numpft))
     allocate (m_livecrootc_storage_to_litter_fire_p(numpft))
     allocate (m_deadcrootc_storage_to_litter_fire_p(numpft))
     allocate (m_gresp_storage_to_litter_fire_p     (numpft))

     allocate (m_leafc_xfer_to_litter_fire_p        (numpft))
     allocate (m_frootc_xfer_to_litter_fire_p       (numpft))
     allocate (m_livestemc_xfer_to_litter_fire_p    (numpft))
     allocate (m_deadstemc_xfer_to_litter_fire_p    (numpft))
     allocate (m_livecrootc_xfer_to_litter_fire_p   (numpft))
     allocate (m_deadcrootc_xfer_to_litter_fire_p   (numpft))
     allocate (m_gresp_xfer_to_litter_fire_p        (numpft))

     allocate (cpool_to_xsmrpool_p          (numpft))
     allocate (cpool_to_gresp_storage_p     (numpft))
     allocate (cpool_to_leafc_p             (numpft))
     allocate (cpool_to_leafc_storage_p     (numpft))
     allocate (cpool_to_frootc_p            (numpft))
     allocate (cpool_to_frootc_storage_p    (numpft))
     allocate (cpool_to_livestemc_p         (numpft))
     allocate (cpool_to_livestemc_storage_p (numpft))
     allocate (cpool_to_deadstemc_p         (numpft))
     allocate (cpool_to_deadstemc_storage_p (numpft))
     allocate (cpool_to_livecrootc_p        (numpft))
     allocate (cpool_to_livecrootc_storage_p(numpft))
     allocate (cpool_to_deadcrootc_p        (numpft))
     allocate (cpool_to_deadcrootc_storage_p(numpft)) 
     allocate (cpool_to_grainc_p            (numpft))
     allocate (cpool_to_grainc_storage_p    (numpft))

     allocate (leaf_xsmr_p                  (numpft))
     allocate (froot_xsmr_p                 (numpft))
     allocate (livestem_xsmr_p              (numpft))
     allocate (livecroot_xsmr_p             (numpft))
     allocate (grain_xsmr_p                 (numpft))

     allocate (cpool_leaf_gr_p                      (numpft))
     allocate (cpool_froot_gr_p                     (numpft))
     allocate (cpool_livestem_gr_p                  (numpft))
     allocate (cpool_deadstem_gr_p                  (numpft))
     allocate (cpool_livecroot_gr_p                 (numpft))
     allocate (cpool_deadcroot_gr_p                 (numpft))
     allocate (cpool_grain_gr_p                     (numpft))

     allocate (cpool_leaf_storage_gr_p              (numpft))
     allocate (cpool_froot_storage_gr_p             (numpft))
     allocate (cpool_livestem_storage_gr_p          (numpft))
     allocate (cpool_deadstem_storage_gr_p          (numpft))
     allocate (cpool_livecroot_storage_gr_p         (numpft))
     allocate (cpool_deadcroot_storage_gr_p         (numpft))
     allocate (cpool_grain_storage_gr_p             (numpft))

     allocate (transfer_leaf_gr_p                   (numpft))
     allocate (transfer_froot_gr_p                  (numpft))
     allocate (transfer_livestem_gr_p               (numpft))
     allocate (transfer_deadstem_gr_p               (numpft))
     allocate (transfer_livecroot_gr_p              (numpft))
     allocate (transfer_deadcroot_gr_p              (numpft))
     allocate (transfer_grain_gr_p                  (numpft))

     allocate (xsmrpool_to_atm_p                    (numpft))

     allocate (cropprod1c_loss_p                    (numpft))

     allocate (plant_ndemand_p              (numpft))

     allocate (leafn_xfer_to_leafn_p                (numpft))
     allocate (frootn_xfer_to_frootn_p              (numpft))
     allocate (livestemn_xfer_to_livestemn_p        (numpft))
     allocate (deadstemn_xfer_to_deadstemn_p        (numpft))
     allocate (livecrootn_xfer_to_livecrootn_p      (numpft))
     allocate (deadcrootn_xfer_to_deadcrootn_p      (numpft))
     allocate (grainn_xfer_to_grainn_p              (numpft))

     allocate (leafn_storage_to_xfer_p              (numpft))
     allocate (frootn_storage_to_xfer_p             (numpft))
     allocate (livestemn_storage_to_xfer_p          (numpft))
     allocate (deadstemn_storage_to_xfer_p          (numpft))
     allocate (livecrootn_storage_to_xfer_p         (numpft))
     allocate (deadcrootn_storage_to_xfer_p         (numpft))
     allocate (grainn_storage_to_xfer_p             (numpft))

     allocate (leafn_to_litter_p                    (numpft))
     allocate (frootn_to_litter_p                   (numpft))
     allocate (grainn_to_food_p                     (numpft))
     allocate (grainn_to_seed_p                     (numpft))
     allocate (crop_seedn_to_leaf_p                 (numpft))
     allocate (livestemn_to_litter_p                (numpft))
     allocate (livestemn_to_deadstemn_p             (numpft))
     allocate (livecrootn_to_deadcrootn_p           (numpft))

     allocate (leafn_to_retransn_p                  (numpft))
     allocate (frootn_to_retransn_p                 (numpft))
     allocate (livestemn_to_retransn_p              (numpft))
     allocate (livecrootn_to_retransn_p             (numpft))
     allocate (retransn_to_npool_p                  (numpft))
     allocate (free_retransn_to_npool_p             (numpft))

     allocate (m_leafn_to_litter_p                  (numpft))
     allocate (m_frootn_to_litter_p                 (numpft))
     allocate (m_livestemn_to_litter_p              (numpft))
     allocate (m_deadstemn_to_litter_p              (numpft))
     allocate (m_livecrootn_to_litter_p             (numpft))
     allocate (m_deadcrootn_to_litter_p             (numpft))
     allocate (m_retransn_to_litter_p               (numpft))

     allocate (m_leafn_storage_to_litter_p          (numpft))
     allocate (m_frootn_storage_to_litter_p         (numpft))
     allocate (m_livestemn_storage_to_litter_p      (numpft))
     allocate (m_deadstemn_storage_to_litter_p      (numpft))
     allocate (m_livecrootn_storage_to_litter_p     (numpft))
     allocate (m_deadcrootn_storage_to_litter_p     (numpft))

     allocate (m_leafn_xfer_to_litter_p             (numpft))
     allocate (m_frootn_xfer_to_litter_p            (numpft))
     allocate (m_livestemn_xfer_to_litter_p         (numpft))
     allocate (m_deadstemn_xfer_to_litter_p         (numpft))
     allocate (m_livecrootn_xfer_to_litter_p        (numpft))
     allocate (m_deadcrootn_xfer_to_litter_p        (numpft))

     allocate (m_leafn_to_fire_p                    (numpft))
     allocate (m_frootn_to_fire_p                   (numpft))
     allocate (m_livestemn_to_fire_p                (numpft))
     allocate (m_deadstemn_to_fire_p                (numpft))
     allocate (m_livecrootn_to_fire_p               (numpft))
     allocate (m_deadcrootn_to_fire_p               (numpft))

     allocate (m_leafn_storage_to_fire_p            (numpft))
     allocate (m_frootn_storage_to_fire_p           (numpft))
     allocate (m_livestemn_storage_to_fire_p        (numpft))
     allocate (m_deadstemn_storage_to_fire_p        (numpft))
     allocate (m_livecrootn_storage_to_fire_p       (numpft))
     allocate (m_deadcrootn_storage_to_fire_p       (numpft))

     allocate (m_leafn_xfer_to_fire_p               (numpft))
     allocate (m_frootn_xfer_to_fire_p              (numpft))
     allocate (m_livestemn_xfer_to_fire_p           (numpft))
     allocate (m_deadstemn_xfer_to_fire_p           (numpft))
     allocate (m_livecrootn_xfer_to_fire_p          (numpft))
     allocate (m_deadcrootn_xfer_to_fire_p          (numpft))

     allocate (m_livestemn_to_deadstemn_fire_p      (numpft))
     allocate (m_livecrootn_to_deadcrootn_fire_p    (numpft))

     allocate (m_retransn_to_fire_p                 (numpft))

     allocate (m_leafn_to_litter_fire_p             (numpft))
     allocate (m_frootn_to_litter_fire_p            (numpft))
     allocate (m_livestemn_to_litter_fire_p         (numpft))
     allocate (m_deadstemn_to_litter_fire_p         (numpft))
     allocate (m_livecrootn_to_litter_fire_p        (numpft))
     allocate (m_deadcrootn_to_litter_fire_p        (numpft))

     allocate (m_leafn_storage_to_litter_fire_p     (numpft))
     allocate (m_frootn_storage_to_litter_fire_p    (numpft))
     allocate (m_livestemn_storage_to_litter_fire_p (numpft))
     allocate (m_deadstemn_storage_to_litter_fire_p (numpft))
     allocate (m_livecrootn_storage_to_litter_fire_p(numpft))
     allocate (m_deadcrootn_storage_to_litter_fire_p(numpft))

     allocate (m_leafn_xfer_to_litter_fire_p        (numpft))
     allocate (m_frootn_xfer_to_litter_fire_p       (numpft))
     allocate (m_livestemn_xfer_to_litter_fire_p    (numpft))
     allocate (m_deadstemn_xfer_to_litter_fire_p    (numpft))
     allocate (m_livecrootn_xfer_to_litter_fire_p   (numpft))
     allocate (m_deadcrootn_xfer_to_litter_fire_p   (numpft))

     allocate (m_retransn_to_litter_fire_p          (numpft))

     allocate (npool_to_leafn_p                     (numpft))
     allocate (npool_to_leafn_storage_p             (numpft))
     allocate (npool_to_frootn_p                    (numpft))
     allocate (npool_to_frootn_storage_p            (numpft))
     allocate (npool_to_livestemn_p                 (numpft))
     allocate (npool_to_livestemn_storage_p         (numpft))
     allocate (npool_to_deadstemn_p                 (numpft))
     allocate (npool_to_deadstemn_storage_p         (numpft))
     allocate (npool_to_livecrootn_p                (numpft))
     allocate (npool_to_livecrootn_storage_p        (numpft))
     allocate (npool_to_deadcrootn_p                (numpft))
     allocate (npool_to_deadcrootn_storage_p        (numpft))
     allocate (npool_to_grainn_p                    (numpft))
     allocate (npool_to_grainn_storage_p            (numpft))

     allocate (respcsun_p     (numpft)) !sunlit leaf respiration
     allocate (respcsha_p     (numpft)) !shaded leaf respiration
     allocate (leaf_mr_p      (numpft)) !leaf maintenance respiration
     allocate (froot_mr_p     (numpft)) !fine root maintenance respiration
     allocate (livestem_mr_p  (numpft)) !live stem maintenance respiration
     allocate (livecroot_mr_p (numpft)) !live coarse root maintenance respiration
     allocate (grain_mr_p     (numpft)) !grain maintenance respiration

     allocate (soil_change_p  (numpft))

     allocate (psn_to_cpool_p               (numpft))
     allocate (gpp_p                        (numpft))
     allocate (availc_p                     (numpft))
     allocate (avail_retransn_p             (numpft))
     allocate (xsmrpool_recover_p           (numpft))
     allocate (excess_cflux_p               (numpft))
     allocate (sminn_to_npool_p             (numpft))

     allocate (plant_calloc_p               (numpft))
     allocate (plant_nalloc_p               (numpft))
     allocate (leaf_curmr_p                 (numpft))
     allocate (froot_curmr_p                (numpft))
     allocate (livestem_curmr_p             (numpft))
     allocate (livecroot_curmr_p            (numpft))
     allocate (grain_curmr_p                (numpft))

     allocate (fire_closs_p                 (numpft))
     allocate (fire_nloss_p                 (numpft))
     allocate (wood_harvestc_p              (numpft))
     allocate (wood_harvestn_p              (numpft))
     allocate (grainc_to_cropprodc_p        (numpft))
     allocate (grainn_to_cropprodn_p        (numpft))
     allocate (hrv_xsmrpool_to_atm_p        (numpft))
     allocate (fert_p                       (numpft))
     allocate (soyfixn_p                    (numpft))
!--------

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

! bgc variables
     deallocate (leafc_xfer_to_leafc_p                )
     deallocate (frootc_xfer_to_frootc_p              )
     deallocate (livestemc_xfer_to_livestemc_p        )
     deallocate (deadstemc_xfer_to_deadstemc_p        )
     deallocate (livecrootc_xfer_to_livecrootc_p      )
     deallocate (deadcrootc_xfer_to_deadcrootc_p      )
     deallocate (grainc_xfer_to_grainc_p              )

     deallocate (leafc_storage_to_xfer_p              )
     deallocate (frootc_storage_to_xfer_p             )
     deallocate (livestemc_storage_to_xfer_p          )
     deallocate (deadstemc_storage_to_xfer_p          )
     deallocate (livecrootc_storage_to_xfer_p         )
     deallocate (deadcrootc_storage_to_xfer_p         )
     deallocate (grainc_storage_to_xfer_p             )
     deallocate (gresp_storage_to_xfer_p              )

     deallocate (leafc_to_litter_p                    )
     deallocate (frootc_to_litter_p                   )
     deallocate (grainc_to_food_p                     )
     deallocate (grainc_to_seed_p                     )
     deallocate (crop_seedc_to_leaf_p                 )
     deallocate (livestemc_to_litter_p                )
     deallocate (livestemc_to_deadstemc_p             )
     deallocate (livecrootc_to_deadcrootc_p           )

     deallocate (m_leafc_to_litter_p                  )
     deallocate (m_frootc_to_litter_p                 )
     deallocate (m_livestemc_to_litter_p              )
     deallocate (m_deadstemc_to_litter_p              )
     deallocate (m_livecrootc_to_litter_p             )
     deallocate (m_deadcrootc_to_litter_p             )

     deallocate (m_leafc_storage_to_litter_p          )
     deallocate (m_frootc_storage_to_litter_p         )
     deallocate (m_livestemc_storage_to_litter_p      )
     deallocate (m_deadstemc_storage_to_litter_p      )
     deallocate (m_livecrootc_storage_to_litter_p     )
     deallocate (m_deadcrootc_storage_to_litter_p     )
     deallocate (m_gresp_storage_to_litter_p          )

     deallocate (m_leafc_xfer_to_litter_p             )
     deallocate (m_frootc_xfer_to_litter_p            )
     deallocate (m_livestemc_xfer_to_litter_p         )
     deallocate (m_deadstemc_xfer_to_litter_p         )
     deallocate (m_livecrootc_xfer_to_litter_p        )
     deallocate (m_deadcrootc_xfer_to_litter_p        )
     deallocate (m_gresp_xfer_to_litter_p             )

     deallocate (m_leafc_to_fire_p                    )
     deallocate (m_frootc_to_fire_p                   )
     deallocate (m_livestemc_to_fire_p                )
     deallocate (m_deadstemc_to_fire_p                )
     deallocate (m_livecrootc_to_fire_p               )
     deallocate (m_deadcrootc_to_fire_p               )

     deallocate (m_leafc_storage_to_fire_p            )
     deallocate (m_frootc_storage_to_fire_p           )
     deallocate (m_livestemc_storage_to_fire_p        )
     deallocate (m_deadstemc_storage_to_fire_p        )
     deallocate (m_livecrootc_storage_to_fire_p       )
     deallocate (m_deadcrootc_storage_to_fire_p       )
     deallocate (m_gresp_storage_to_fire_p            )

     deallocate (m_leafc_xfer_to_fire_p               )
     deallocate (m_frootc_xfer_to_fire_p              )
     deallocate (m_livestemc_xfer_to_fire_p           )
     deallocate (m_deadstemc_xfer_to_fire_p           )
     deallocate (m_livecrootc_xfer_to_fire_p          )
     deallocate (m_deadcrootc_xfer_to_fire_p          )
     deallocate (m_gresp_xfer_to_fire_p               )

     deallocate (m_livestemc_to_deadstemc_fire_p      )
     deallocate (m_livecrootc_to_deadcrootc_fire_p    )

     deallocate (m_leafc_to_litter_fire_p             )
     deallocate (m_frootc_to_litter_fire_p            )
     deallocate (m_livestemc_to_litter_fire_p         )
     deallocate (m_deadstemc_to_litter_fire_p         )
     deallocate (m_livecrootc_to_litter_fire_p        )
     deallocate (m_deadcrootc_to_litter_fire_p        )

     deallocate (m_leafc_storage_to_litter_fire_p     )
     deallocate (m_frootc_storage_to_litter_fire_p    )
     deallocate (m_livestemc_storage_to_litter_fire_p )
     deallocate (m_deadstemc_storage_to_litter_fire_p )
     deallocate (m_livecrootc_storage_to_litter_fire_p)
     deallocate (m_deadcrootc_storage_to_litter_fire_p)
     deallocate (m_gresp_storage_to_litter_fire_p     )

     deallocate (m_leafc_xfer_to_litter_fire_p        )
     deallocate (m_frootc_xfer_to_litter_fire_p       )
     deallocate (m_livestemc_xfer_to_litter_fire_p    )
     deallocate (m_deadstemc_xfer_to_litter_fire_p    )
     deallocate (m_livecrootc_xfer_to_litter_fire_p   )
     deallocate (m_deadcrootc_xfer_to_litter_fire_p   )
     deallocate (m_gresp_xfer_to_litter_fire_p        )

     deallocate (cpool_to_xsmrpool_p          )
     deallocate (cpool_to_gresp_storage_p     )
     deallocate (cpool_to_leafc_p             )
     deallocate (cpool_to_leafc_storage_p     )
     deallocate (cpool_to_frootc_p            )
     deallocate (cpool_to_frootc_storage_p    )
     deallocate (cpool_to_livestemc_p         )
     deallocate (cpool_to_livestemc_storage_p )
     deallocate (cpool_to_deadstemc_p         )
     deallocate (cpool_to_deadstemc_storage_p )
     deallocate (cpool_to_livecrootc_p        )
     deallocate (cpool_to_livecrootc_storage_p)
     deallocate (cpool_to_deadcrootc_p        )
     deallocate (cpool_to_deadcrootc_storage_p)
     deallocate (cpool_to_grainc_p            )
     deallocate (cpool_to_grainc_storage_p    )

     deallocate (leaf_xsmr_p                  )
     deallocate (froot_xsmr_p                 )
     deallocate (livestem_xsmr_p              )
     deallocate (livecroot_xsmr_p             )
     deallocate (grain_xsmr_p                 )

     deallocate (cpool_leaf_gr_p                      )
     deallocate (cpool_froot_gr_p                     )
     deallocate (cpool_livestem_gr_p                  )
     deallocate (cpool_deadstem_gr_p                  )
     deallocate (cpool_livecroot_gr_p                 )
     deallocate (cpool_deadcroot_gr_p                 )
     deallocate (cpool_grain_gr_p                     )

     deallocate (cpool_leaf_storage_gr_p              )
     deallocate (cpool_froot_storage_gr_p             )
     deallocate (cpool_livestem_storage_gr_p          )
     deallocate (cpool_deadstem_storage_gr_p          )
     deallocate (cpool_livecroot_storage_gr_p         )
     deallocate (cpool_deadcroot_storage_gr_p         )
     deallocate (cpool_grain_storage_gr_p             )

     deallocate (transfer_leaf_gr_p                   )
     deallocate (transfer_froot_gr_p                  )
     deallocate (transfer_livestem_gr_p               )
     deallocate (transfer_deadstem_gr_p               )
     deallocate (transfer_livecroot_gr_p              )
     deallocate (transfer_deadcroot_gr_p              )
     deallocate (transfer_grain_gr_p                  )

     deallocate (xsmrpool_to_atm_p                    )

     deallocate (cropprod1c_loss_p                    )

     deallocate (plant_ndemand_p              )

     deallocate (leafn_xfer_to_leafn_p                )
     deallocate (frootn_xfer_to_frootn_p              )
     deallocate (livestemn_xfer_to_livestemn_p        )
     deallocate (deadstemn_xfer_to_deadstemn_p        )
     deallocate (livecrootn_xfer_to_livecrootn_p      )
     deallocate (deadcrootn_xfer_to_deadcrootn_p      )
     deallocate (grainn_xfer_to_grainn_p              )

     deallocate (leafn_storage_to_xfer_p              )
     deallocate (frootn_storage_to_xfer_p             )
     deallocate (livestemn_storage_to_xfer_p          )
     deallocate (deadstemn_storage_to_xfer_p          )
     deallocate (livecrootn_storage_to_xfer_p         )
     deallocate (deadcrootn_storage_to_xfer_p         )
     deallocate (grainn_storage_to_xfer_p             )

     deallocate (leafn_to_litter_p                    )
     deallocate (frootn_to_litter_p                   )
     deallocate (grainn_to_food_p                     )
     deallocate (grainn_to_seed_p                     )
     deallocate (crop_seedn_to_leaf_p                 )
     deallocate (livestemn_to_litter_p                )
     deallocate (livestemn_to_deadstemn_p             )
     deallocate (livecrootn_to_deadcrootn_p           )

     deallocate (leafn_to_retransn_p                  )
     deallocate (frootn_to_retransn_p                 )
     deallocate (livestemn_to_retransn_p              )
     deallocate (livecrootn_to_retransn_p             )
     deallocate (retransn_to_npool_p                  )
     deallocate (free_retransn_to_npool_p             )

     deallocate (m_leafn_to_litter_p                  )
     deallocate (m_frootn_to_litter_p                 )
     deallocate (m_livestemn_to_litter_p              )
     deallocate (m_deadstemn_to_litter_p              )
     deallocate (m_livecrootn_to_litter_p             )
     deallocate (m_deadcrootn_to_litter_p             )
     deallocate (m_retransn_to_litter_p               )

     deallocate (m_leafn_storage_to_litter_p          )
     deallocate (m_frootn_storage_to_litter_p         )
     deallocate (m_livestemn_storage_to_litter_p      )
     deallocate (m_deadstemn_storage_to_litter_p      )
     deallocate (m_livecrootn_storage_to_litter_p     )
     deallocate (m_deadcrootn_storage_to_litter_p     )

     deallocate (m_leafn_xfer_to_litter_p             )
     deallocate (m_frootn_xfer_to_litter_p            )
     deallocate (m_livestemn_xfer_to_litter_p         )
     deallocate (m_deadstemn_xfer_to_litter_p         )
     deallocate (m_livecrootn_xfer_to_litter_p        )
     deallocate (m_deadcrootn_xfer_to_litter_p        )

     deallocate (m_leafn_to_fire_p                    )
     deallocate (m_frootn_to_fire_p                   )
     deallocate (m_livestemn_to_fire_p                )
     deallocate (m_deadstemn_to_fire_p                )
     deallocate (m_livecrootn_to_fire_p               )
     deallocate (m_deadcrootn_to_fire_p               )

     deallocate (m_leafn_storage_to_fire_p            )
     deallocate (m_frootn_storage_to_fire_p           )
     deallocate (m_livestemn_storage_to_fire_p        )
     deallocate (m_deadstemn_storage_to_fire_p        )
     deallocate (m_livecrootn_storage_to_fire_p       )
     deallocate (m_deadcrootn_storage_to_fire_p       )

     deallocate (m_leafn_xfer_to_fire_p               )
     deallocate (m_frootn_xfer_to_fire_p              )
     deallocate (m_livestemn_xfer_to_fire_p           )
     deallocate (m_deadstemn_xfer_to_fire_p           )
     deallocate (m_livecrootn_xfer_to_fire_p          )
     deallocate (m_deadcrootn_xfer_to_fire_p          )

     deallocate (m_livestemn_to_deadstemn_fire_p      )
     deallocate (m_livecrootn_to_deadcrootn_fire_p    )

     deallocate (m_retransn_to_fire_p                 )

     deallocate (m_leafn_to_litter_fire_p             )
     deallocate (m_frootn_to_litter_fire_p            )
     deallocate (m_livestemn_to_litter_fire_p         )
     deallocate (m_deadstemn_to_litter_fire_p         )
     deallocate (m_livecrootn_to_litter_fire_p        )
     deallocate (m_deadcrootn_to_litter_fire_p        )

     deallocate (m_leafn_storage_to_litter_fire_p     )
     deallocate (m_frootn_storage_to_litter_fire_p    )
     deallocate (m_livestemn_storage_to_litter_fire_p )
     deallocate (m_deadstemn_storage_to_litter_fire_p )
     deallocate (m_livecrootn_storage_to_litter_fire_p)
     deallocate (m_deadcrootn_storage_to_litter_fire_p)

     deallocate (m_leafn_xfer_to_litter_fire_p        )
     deallocate (m_frootn_xfer_to_litter_fire_p       )
     deallocate (m_livestemn_xfer_to_litter_fire_p    )
     deallocate (m_deadstemn_xfer_to_litter_fire_p    )
     deallocate (m_livecrootn_xfer_to_litter_fire_p   )
     deallocate (m_deadcrootn_xfer_to_litter_fire_p   )

     deallocate (m_retransn_to_litter_fire_p          )

     deallocate (npool_to_leafn_p                     )
     deallocate (npool_to_leafn_storage_p             )
     deallocate (npool_to_frootn_p                    )
     deallocate (npool_to_frootn_storage_p            )
     deallocate (npool_to_livestemn_p                 )
     deallocate (npool_to_livestemn_storage_p         )
     deallocate (npool_to_deadstemn_p                 )
     deallocate (npool_to_deadstemn_storage_p         )
     deallocate (npool_to_livecrootn_p                )
     deallocate (npool_to_livecrootn_storage_p        )
     deallocate (npool_to_deadcrootn_p                )
     deallocate (npool_to_deadcrootn_storage_p        )
     deallocate (npool_to_grainn_p                    )
     deallocate (npool_to_grainn_storage_p            )

     deallocate (respcsun_p     )
     deallocate (respcsha_p     )
     deallocate (leaf_mr_p      )
     deallocate (froot_mr_p     )
     deallocate (livestem_mr_p  )
     deallocate (livecroot_mr_p )
     deallocate (grain_mr_p     )

     deallocate (soil_change_p  )

     deallocate (psn_to_cpool_p               )
     deallocate (gpp_p                        )
     deallocate (availc_p                     )
     deallocate (avail_retransn_p             )
     deallocate (xsmrpool_recover_p           )
     deallocate (excess_cflux_p               )
     deallocate (sminn_to_npool_p             )

     deallocate (plant_calloc_p               )
     deallocate (plant_nalloc_p               )
     deallocate (leaf_curmr_p                 )
     deallocate (froot_curmr_p                )
     deallocate (livestem_curmr_p             )
     deallocate (livecroot_curmr_p            )
     deallocate (grain_curmr_p                )

     deallocate (fire_closs_p                 )
     deallocate (fire_nloss_p                 )
     deallocate (wood_harvestc_p              )
     deallocate (wood_harvestn_p              )
     deallocate (grainc_to_cropprodc_p        )
     deallocate (grainn_to_cropprodn_p        )
     deallocate (hrv_xsmrpool_to_atm_p        )
     deallocate (fert_p                       )
     deallocate (soyfixn_p                    )

  END SUBROUTINE deallocate_1D_PFTFluxes

END MODULE MOD_1D_PFTFluxes
! ---------- EOP ------------
