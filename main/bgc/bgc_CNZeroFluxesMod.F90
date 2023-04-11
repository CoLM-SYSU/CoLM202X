#include <define.h>
#ifdef BGC

module bgc_CNZeroFluxesMod

use precision

use MOD_1D_BGCPFTFluxes, only:            &
   m_leafc_to_litter_p                  , &
   m_frootc_to_litter_p                 , &
   m_leafc_storage_to_litter_p          , &
   m_frootc_storage_to_litter_p         , &
   m_livestemc_storage_to_litter_p      , &
   m_deadstemc_storage_to_litter_p      , &
   m_livecrootc_storage_to_litter_p     , &
   m_deadcrootc_storage_to_litter_p     , &
   m_leafc_xfer_to_litter_p             , &
   m_frootc_xfer_to_litter_p            , &  
   m_livestemc_xfer_to_litter_p         , &
   m_deadstemc_xfer_to_litter_p         , &
   m_livecrootc_xfer_to_litter_p        , &
   m_deadcrootc_xfer_to_litter_p        , &
   m_livestemc_to_litter_p              , &
   m_deadstemc_to_litter_p              , &
   m_livecrootc_to_litter_p             , &
   m_deadcrootc_to_litter_p             , &
   m_gresp_storage_to_litter_p          , &
   m_gresp_xfer_to_litter_p             , &
   m_leafc_to_fire_p                    , &
   m_leafc_storage_to_fire_p            , &
   m_leafc_xfer_to_fire_p               , &
   m_livestemc_to_fire_p                , &
   m_livestemc_storage_to_fire_p        , &
   m_livestemc_xfer_to_fire_p           , &
   m_deadstemc_to_fire_p                , &
   m_deadstemc_storage_to_fire_p        , &
   m_deadstemc_xfer_to_fire_p           , &
   m_frootc_to_fire_p                   , &
   m_frootc_storage_to_fire_p           , &
   m_frootc_xfer_to_fire_p              , &
   m_livecrootc_to_fire_p               , &
   m_livecrootc_storage_to_fire_p       , &
   m_livecrootc_xfer_to_fire_p          , &
   m_deadcrootc_to_fire_p               , &
   m_deadcrootc_storage_to_fire_p       , &
   m_deadcrootc_xfer_to_fire_p          , &
   m_gresp_storage_to_fire_p            , &
   m_gresp_xfer_to_fire_p               , &
   
   m_leafc_to_litter_fire_p             , &
   m_leafc_storage_to_litter_fire_p     , &
   m_leafc_xfer_to_litter_fire_p        , &
   m_livestemc_to_litter_fire_p         , &
   m_livestemc_storage_to_litter_fire_p , &
   m_livestemc_xfer_to_litter_fire_p    , &
   m_livestemc_to_deadstemc_fire_p      , &
   m_deadstemc_to_litter_fire_p         , &
   m_deadstemc_storage_to_litter_fire_p , &
   m_deadstemc_xfer_to_litter_fire_p    , &
   m_frootc_to_litter_fire_p            , &
   m_frootc_storage_to_litter_fire_p    , &
   m_frootc_xfer_to_litter_fire_p       , &
   m_livecrootc_to_litter_fire_p        , &
   m_livecrootc_storage_to_litter_fire_p, &
   m_livecrootc_xfer_to_litter_fire_p   , &
   m_livecrootc_to_deadcrootc_fire_p    , &
   m_deadcrootc_to_litter_fire_p        , &
   m_deadcrootc_storage_to_litter_fire_p, &
   m_deadcrootc_xfer_to_litter_fire_p   , &
   m_gresp_storage_to_litter_fire_p     , &
   m_gresp_xfer_to_litter_fire_p        , &

   leafc_xfer_to_leafc_p                , &
   frootc_xfer_to_frootc_p              , &
   livestemc_xfer_to_livestemc_p        , &
   deadstemc_xfer_to_deadstemc_p        , &
   livecrootc_xfer_to_livecrootc_p      , &
   deadcrootc_xfer_to_deadcrootc_p      , &
   leafc_to_litter_p                    , &
   frootc_to_litter_p                   , &

   leaf_mr_p                            , &
   froot_mr_p                           , &
   livestem_mr_p                        , &
   livecroot_mr_p                       , &
   grain_mr_p                           , &
   leaf_curmr_p                         , &
   froot_curmr_p                        , &
   livestem_curmr_p                     , &
   livecroot_curmr_p                    , &
   grain_curmr_p                        , &
   leaf_xsmr_p                          , &
   froot_xsmr_p                         , &
   livestem_xsmr_p                      , &
   livecroot_xsmr_p                     , &
   grain_xsmr_p                         , &
   psn_to_cpool_p                       , &
   cpool_to_xsmrpool_p                  , &
   cpool_to_leafc_p                     , &
   cpool_to_leafc_storage_p             , &
   cpool_to_frootc_p                    , &
   cpool_to_frootc_storage_p            , &
   cpool_to_livestemc_p                 , &
   cpool_to_livestemc_storage_p         , &
   cpool_to_deadstemc_p                 , &
   cpool_to_deadstemc_storage_p         , &
   cpool_to_livecrootc_p                , &
   cpool_to_livecrootc_storage_p        , &
   cpool_to_deadcrootc_p                , &
   cpool_to_deadcrootc_storage_p        , &
   cpool_to_gresp_storage_p             , &
   cpool_leaf_gr_p                      , &
   cpool_leaf_storage_gr_p              , &
   transfer_leaf_gr_p                   , &
   cpool_froot_gr_p                     , &
   cpool_froot_storage_gr_p             , &
   transfer_froot_gr_p                  , &
   cpool_livestem_gr_p                  , &
   cpool_livestem_storage_gr_p          , &
   transfer_livestem_gr_p               , &
   cpool_deadstem_gr_p                  , &
   cpool_deadstem_storage_gr_p          , &
   transfer_deadstem_gr_p               , &
   cpool_livecroot_gr_p                 , &
   cpool_livecroot_storage_gr_p         , &
   transfer_livecroot_gr_p              , &
   cpool_deadcroot_gr_p                 , &
   cpool_deadcroot_storage_gr_p         , &
   transfer_deadcroot_gr_p              , &
   leafc_storage_to_xfer_p              , &
   frootc_storage_to_xfer_p             , &
   livestemc_storage_to_xfer_p          , &
   deadstemc_storage_to_xfer_p          , &
   livecrootc_storage_to_xfer_p         , &
   deadcrootc_storage_to_xfer_p         , &
   gresp_storage_to_xfer_p              , &
   livestemc_to_deadstemc_p             , &
   livecrootc_to_deadcrootc_p           , &
   crop_seedc_to_leaf_p                 , &

   hrv_xsmrpool_to_atm_p                , &

   xsmrpool_to_atm_p                    , &
   livestemc_to_litter_p                , &
   grainc_to_food_p                     , &
   grainc_to_seed_p                     , &
   grainc_xfer_to_grainc_p              , &
   cpool_to_grainc_p                    , &
   cpool_to_grainc_storage_p            , &
   cpool_grain_gr_p                     , &
   cpool_grain_storage_gr_p             , &
   transfer_grain_gr_p                  , &
   grainc_storage_to_xfer_p             , &

   gpp_p, &
   
   m_leafn_to_litter_p                  , &
   m_frootn_to_litter_p                 , &
   m_leafn_storage_to_litter_p          , &
   m_frootn_storage_to_litter_p         , &
   m_livestemn_storage_to_litter_p      , &
   m_deadstemn_storage_to_litter_p      , &
   m_livecrootn_storage_to_litter_p     , &
   m_deadcrootn_storage_to_litter_p     , &
   m_leafn_xfer_to_litter_p             , &
   m_frootn_xfer_to_litter_p            , &
   m_livestemn_xfer_to_litter_p         , &
   m_deadstemn_xfer_to_litter_p         , &
   m_livecrootn_xfer_to_litter_p        , &
   m_deadcrootn_xfer_to_litter_p        , &
   m_livestemn_to_litter_p              , &
   m_deadstemn_to_litter_p              , &
   m_livecrootn_to_litter_p             , &
   m_deadcrootn_to_litter_p             , &
   m_retransn_to_litter_p               , &

   m_leafn_to_fire_p                    , &
   m_leafn_storage_to_fire_p            , &
   m_leafn_xfer_to_fire_p               , &
   m_livestemn_to_fire_p                , &
   m_livestemn_storage_to_fire_p        , &
   m_livestemn_xfer_to_fire_p           , &
   m_deadstemn_to_fire_p                , &
   m_deadstemn_storage_to_fire_p        , &
   m_deadstemn_xfer_to_fire_p           , &
   m_frootn_to_fire_p                   , &
   m_frootn_storage_to_fire_p           , &
   m_frootn_xfer_to_fire_p              , &
   m_livecrootn_to_fire_p               , &
   m_livecrootn_storage_to_fire_p       , &
   m_livecrootn_xfer_to_fire_p          , &
   m_deadcrootn_to_fire_p               , &
   m_deadcrootn_storage_to_fire_p       , &
   m_deadcrootn_xfer_to_fire_p          , &
   m_retransn_to_fire_p                 , &
   
   
   m_leafn_to_litter_fire_p             , &
   m_leafn_storage_to_litter_fire_p     , &
   m_leafn_xfer_to_litter_fire_p        , &
   m_livestemn_to_litter_fire_p         , &
   m_livestemn_storage_to_litter_fire_p , &
   m_livestemn_xfer_to_litter_fire_p    , &
   m_livestemn_to_deadstemn_fire_p      , &
   m_deadstemn_to_litter_fire_p         , &
   m_deadstemn_storage_to_litter_fire_p , &
   m_deadstemn_xfer_to_litter_fire_p    , &
   m_frootn_to_litter_fire_p            , &
   m_frootn_storage_to_litter_fire_p    , &
   m_frootn_xfer_to_litter_fire_p       , &
   m_livecrootn_to_litter_fire_p        , &
   m_livecrootn_storage_to_litter_fire_p, &
   m_livecrootn_xfer_to_litter_fire_p   , &
   m_livecrootn_to_deadcrootn_fire_p    , &
   m_deadcrootn_to_litter_fire_p        , &
   m_deadcrootn_storage_to_litter_fire_p, &
   m_deadcrootn_xfer_to_litter_fire_p   , &
   m_retransn_to_litter_fire_p          , &

   leafn_xfer_to_leafn_p                , &
   frootn_xfer_to_frootn_p              , &
   livestemn_xfer_to_livestemn_p        , &
   deadstemn_xfer_to_deadstemn_p        , &
   livecrootn_xfer_to_livecrootn_p      , &
   deadcrootn_xfer_to_deadcrootn_p      , &
   leafn_to_litter_p                    , &
   leafn_to_retransn_p                  , &
   frootn_to_litter_p                   , &
   retransn_to_npool_p                  , &
   free_retransn_to_npool_p             , &
   sminn_to_npool_p                     , &
   npool_to_leafn_p                     , &
   npool_to_leafn_storage_p             , &
   npool_to_frootn_p                    , &
   npool_to_frootn_storage_p            , &
   npool_to_livestemn_p                 , &
   npool_to_livestemn_storage_p         , &
   npool_to_deadstemn_p                 , &
   npool_to_deadstemn_storage_p         , &
   npool_to_livecrootn_p                , &
   npool_to_livecrootn_storage_p        , &
   npool_to_deadcrootn_p                , &
   npool_to_deadcrootn_storage_p        , &
   leafn_storage_to_xfer_p              , &
   frootn_storage_to_xfer_p             , &
   livestemn_storage_to_xfer_p          , &
   deadstemn_storage_to_xfer_p          , &
   livecrootn_storage_to_xfer_p         , &
   deadcrootn_storage_to_xfer_p         , &
   livestemn_to_deadstemn_p             , &
   livestemn_to_retransn_p              , &
   livecrootn_to_deadcrootn_p           , &
   livecrootn_to_retransn_p             , &
   
   crop_seedn_to_leaf_p                 , &

   livestemn_to_litter_p                , &
   grainn_to_food_p                     , &
   grainn_to_seed_p                     , &
   grainn_xfer_to_grainn_p              , &
   npool_to_grainn_p                    , &
   npool_to_grainn_storage_p            , &
   grainn_storage_to_xfer_p             , &
!   soyfixn                           , &
   frootn_to_retransn_p                 , &      

   fire_closs_p                         , &
   fire_nloss_p                         , &
   wood_harvestc_p                      , &
   wood_harvestn_p                      , &
   grainc_to_cropprodc_p                , &
   grainn_to_cropprodn_p                , &
   soyfixn_p

use MOD_1D_BGCFluxes, only:             &

      ! phenologgy
   phenology_to_met_c                 , &
   phenology_to_cel_c                 , &
   phenology_to_lig_c                 , &
   
      ! gap mortality
   gap_mortality_to_met_c             , &
   gap_mortality_to_cel_c             , &
   gap_mortality_to_lig_c             , &
   gap_mortality_to_cwdc              , &
   
      ! fire
   fire_mortality_to_cwdc             , &
   fire_mortality_to_met_c            , &
   fire_mortality_to_cel_c            , &
   fire_mortality_to_lig_c            , &

   m_decomp_cpools_to_fire_vr         , &

      ! phenologgy
   phenology_to_met_n                 , &
   phenology_to_cel_n                 , &
   phenology_to_lig_n                 , &

      ! gap mortality
   gap_mortality_to_met_n             , &
   gap_mortality_to_cel_n             , &
   gap_mortality_to_lig_n             , &
   gap_mortality_to_cwdn              , &

      ! fire
   fire_mortality_to_cwdn             , &
   fire_mortality_to_met_n            , &
   fire_mortality_to_cel_n            , &
   fire_mortality_to_lig_n            , &

!      ! harvest
!      this%harvest_n_to_litr_met_n(j,i)          = 0._r8
!      this%harvest_n_to_litr_cel_n(j,i)          = 0._r8
!      this%harvest_n_to_litr_lig_n(j,i)          = 0._r8
!      this%harvest_n_to_cwdn(j,i)                = 0._r8
   m_decomp_npools_to_fire_vr , &

   fire_closs   , &
   fire_nloss   , &
   wood_harvestc, &
   wood_harvestn, &
   grainc_to_cropprodc, &
   grainn_to_cropprodn, &

!   npp_Nactive, &
!   npp_burnedoff, &
!   npp_Nnonmyc, &
!   npp_Nam, &
!   npp_Necm, &
!   npp_Nactive_no3, &
!   npp_Nactive_nh4, &
!   npp_Nnonmyc_no3, &
!   npp_Nnonmyc_nh4, &
!   npp_Nam_no3, &
!   npp_Nam_nh4, &
!   npp_Necm_no3, &
!   npp_Necm_nh4, &
!   npp_Nfix, &
!   npp_Nretrans, &
!   npp_Nuptake, &
!   npp_growth
!   leafc_change, &
!   soilc_change, &

   decomp_hr           , &
   decomp_hr_vr        , &
!   decomp_ctransfer    , &
   decomp_ctransfer_vr , &

!   decomp_cpools_leached            , &
   decomp_cpools_transport_tendency , &
   decomp_cpools_sourcesink         , &

!   hr_vr , &

!   hr            , &
   somc_fire      , &
   som_c_leached , &
!   somhr         , &
!   lithr         

#ifndef NITRIF
    sminn_to_denit_excess_vr       , &
    sminn_leached_vr       , &
    sminn_to_plant_fun_vr       , &
#else
    f_nit_vr       , &
    f_denit_vr       , &
    smin_no3_leached_vr       , &
    smin_no3_runoff_vr       , &
    n2_n2o_ratio_denit_vr       , &
    pot_f_nit_vr       , &
    pot_f_denit_vr       , &
    actual_immob_no3_vr       , &
    actual_immob_nh4_vr       , &
    smin_no3_to_plant_vr       , &
    smin_nh4_to_plant_vr       , &
    f_n2o_denit_vr       , &
    f_n2o_nit_vr       , &

#endif
    potential_immob_vr       , &
    actual_immob_vr       , &
    sminn_to_plant       , &
    sminn_to_plant_vr       , &
    supplement_to_sminn_vr       , &
    gross_nmin_vr       , &
    net_nmin_vr       , &
    sminn_to_plant_fun_no3_vr       , &
    sminn_to_plant_fun_nh4_vr       , &

    nfix_to_sminn       , &
    ffix_to_sminn       , &
    fert_to_sminn       , &
    soyfixn_to_sminn       , &
!    potential_immob       , &
!    actual_immob       , &
    sminn_to_plant       , &
    supplement_to_sminn       , &
    gross_nmin       , &
    net_nmin       , &
    denit       , &
!    sminn_to_plant_fun       , &
!#ifndef NITRIF
!    f_nit       , &
!    pot_f_nit       , &
!    f_denit       , &
!    pot_f_denit       , &
!    f_n2o_denit       , &
    f_n2o_nit       , &
    smin_no3_leached       , &
    smin_no3_runoff       , &
!#else
!    sminn_to_denit_excess       , &
    sminn_leached       , &
!#endif
!    ninputs       , &
!    noutputs       , &
    som_n_leached       , &

!    decomp_npools_leached       , &
   
!    if(use_soil_matrixcn)then
!       call this%matrix_Ninput%SetValueV_scaler(num_column,filter_column(1:num_column),value_column)
!    end if

    decomp_npools_transport_tendency       , &

!    decomp_ntransfer       , &
!    decomp_sminn_flux       , &
#ifndef NITRIF
!    sminn_to_denit_decomp       , &
#endif

    decomp_ntransfer_vr       , &
    decomp_sminn_flux_vr       , &
#ifndef NITRIF
    sminn_to_denit_decomp_vr       , &
#endif

    decomp_npools_sourcesink       
   
use MOD_BGCTimeVars, only: &
    decomp_k                  
    
implicit none

public CNZeroFluxes

contains

subroutine CNZeroFluxes (i,ps,pe,nl_soil,ndecomp_pools,ndecomp_transitions)
   

  integer, intent(in) :: i
  integer, intent(in) :: ps
  integer, intent(in) :: pe
  integer, intent(in) :: nl_soil
  integer, intent(in) :: ndecomp_pools
  integer, intent(in) :: ndecomp_transitions

  integer j,k, m

  do m = ps , pe
! CNVegCarbonFluxes set zero
     m_leafc_to_litter_p(m)  = 0._r8
     m_frootc_to_litter_p(m) = 0._r8
     m_leafc_storage_to_litter_p(m) = 0._r8
     m_frootc_storage_to_litter_p(m) = 0._r8
     m_livestemc_storage_to_litter_p(m) = 0._r8
     m_deadstemc_storage_to_litter_p(m) = 0._r8
     m_livecrootc_storage_to_litter_p(m) = 0._r8
     m_deadcrootc_storage_to_litter_p(m) = 0._r8
     m_leafc_xfer_to_litter_p(m) = 0._r8
     m_frootc_xfer_to_litter_p(m) = 0._r8  
     m_livestemc_xfer_to_litter_p(m) = 0._r8
     m_deadstemc_xfer_to_litter_p(m) = 0._r8
     m_livecrootc_xfer_to_litter_p(m) = 0._r8
     m_deadcrootc_xfer_to_litter_p(m) = 0._r8
     m_livestemc_to_litter_p(m) = 0._r8
     m_deadstemc_to_litter_p(m) = 0._r8
     m_livecrootc_to_litter_p(m) = 0._r8
     m_deadcrootc_to_litter_p(m) = 0._r8
     m_gresp_storage_to_litter_p(m) = 0._r8
     m_gresp_xfer_to_litter_p(m) = 0._r8
     m_leafc_to_fire_p(m) = 0._r8
     m_leafc_storage_to_fire_p(m) = 0._r8
     m_leafc_xfer_to_fire_p(m) = 0._r8
     m_livestemc_to_fire_p(m) = 0._r8
     m_livestemc_storage_to_fire_p(m) = 0._r8
     m_livestemc_xfer_to_fire_p(m) = 0._r8
     m_deadstemc_to_fire_p(m) = 0._r8
     m_deadstemc_storage_to_fire_p(m) = 0._r8
     m_deadstemc_xfer_to_fire_p(m) = 0._r8
     m_frootc_to_fire_p(m) = 0._r8
     m_frootc_storage_to_fire_p(m) = 0._r8
     m_frootc_xfer_to_fire_p(m) = 0._r8
     m_livecrootc_to_fire_p(m) = 0._r8
     m_livecrootc_storage_to_fire_p(m) = 0._r8
     m_livecrootc_xfer_to_fire_p(m) = 0._r8
     m_deadcrootc_to_fire_p(m) = 0._r8
     m_deadcrootc_storage_to_fire_p(m) = 0._r8
     m_deadcrootc_xfer_to_fire_p(m) = 0._r8
     m_gresp_storage_to_fire_p(m) = 0._r8
     m_gresp_xfer_to_fire_p(m) = 0._r8
   
     m_leafc_to_litter_fire_p(m)                    = 0._r8
     m_leafc_storage_to_litter_fire_p(m)            = 0._r8
     m_leafc_xfer_to_litter_fire_p(m)               = 0._r8
     m_livestemc_to_litter_fire_p(m)                = 0._r8
     m_livestemc_storage_to_litter_fire_p(m)        = 0._r8
     m_livestemc_xfer_to_litter_fire_p(m)           = 0._r8
     m_livestemc_to_deadstemc_fire_p(m)             = 0._r8
     m_deadstemc_to_litter_fire_p(m)                = 0._r8
     m_deadstemc_storage_to_litter_fire_p(m)        = 0._r8
     m_deadstemc_xfer_to_litter_fire_p(m)           = 0._r8
     m_frootc_to_litter_fire_p(m)                   = 0._r8
     m_frootc_storage_to_litter_fire_p(m)           = 0._r8
     m_frootc_xfer_to_litter_fire_p(m)              = 0._r8
     m_livecrootc_to_litter_fire_p(m)               = 0._r8
     m_livecrootc_storage_to_litter_fire_p(m)       = 0._r8
     m_livecrootc_xfer_to_litter_fire_p(m)          = 0._r8
     m_livecrootc_to_deadcrootc_fire_p(m)           = 0._r8
     m_deadcrootc_to_litter_fire_p(m)               = 0._r8
     m_deadcrootc_storage_to_litter_fire_p(m)       = 0._r8
     m_deadcrootc_xfer_to_litter_fire_p(m)          = 0._r8
     m_gresp_storage_to_litter_fire_p(m)            = 0._r8
     m_gresp_xfer_to_litter_fire_p(m)               = 0._r8
     leafc_xfer_to_leafc_p(m)                 = 0._r8
     frootc_xfer_to_frootc_p(m)               = 0._r8
     livestemc_xfer_to_livestemc_p(m)         = 0._r8
     deadstemc_xfer_to_deadstemc_p(m)         = 0._r8
     livecrootc_xfer_to_livecrootc_p(m)       = 0._r8
     deadcrootc_xfer_to_deadcrootc_p(m)       = 0._r8
     leafc_to_litter_p(m)                     = 0._r8
     frootc_to_litter_p(m)                    = 0._r8
     leaf_mr_p(m)                             = 0._r8
     froot_mr_p(m)                            = 0._r8
     livestem_mr_p(m)                         = 0._r8
     livecroot_mr_p(m)                        = 0._r8
     grain_mr_p(m)                            = 0._r8
     leaf_curmr_p(m)                          = 0._r8
     froot_curmr_p(m)                         = 0._r8
     livestem_curmr_p(m)                      = 0._r8
     livecroot_curmr_p(m)                     = 0._r8
     grain_curmr_p(m)                         = 0._r8
     leaf_xsmr_p(m)                           = 0._r8
     froot_xsmr_p(m)                          = 0._r8
     livestem_xsmr_p(m)                       = 0._r8
     livecroot_xsmr_p(m)                      = 0._r8
     grain_xsmr_p(m)                          = 0._r8
     psn_to_cpool_p(m)                        = 0._r8
     cpool_to_xsmrpool_p(m)                   = 0._r8
     cpool_to_leafc_p(m)                      = 0._r8
     cpool_to_leafc_storage_p(m)              = 0._r8
     cpool_to_frootc_p(m)                     = 0._r8
     cpool_to_frootc_storage_p(m)             = 0._r8
     cpool_to_livestemc_p(m)                  = 0._r8
     cpool_to_livestemc_storage_p(m)          = 0._r8
     cpool_to_deadstemc_p(m)                  = 0._r8
     cpool_to_deadstemc_storage_p(m)          = 0._r8
     cpool_to_livecrootc_p(m)                 = 0._r8
     cpool_to_livecrootc_storage_p(m)         = 0._r8
     cpool_to_deadcrootc_p(m)                 = 0._r8
     cpool_to_deadcrootc_storage_p(m)         = 0._r8
     cpool_to_gresp_storage_p(m)              = 0._r8
     cpool_leaf_gr_p(m)                       = 0._r8
     cpool_leaf_storage_gr_p(m)               = 0._r8
     transfer_leaf_gr_p(m)                    = 0._r8
     cpool_froot_gr_p(m)                      = 0._r8
     cpool_froot_storage_gr_p(m)              = 0._r8
     transfer_froot_gr_p(m)                   = 0._r8
     cpool_livestem_gr_p(m)                   = 0._r8
     cpool_livestem_storage_gr_p(m)           = 0._r8
     transfer_livestem_gr_p(m)                = 0._r8
     cpool_deadstem_gr_p(m)                   = 0._r8
     cpool_deadstem_storage_gr_p(m)           = 0._r8
     transfer_deadstem_gr_p(m)                = 0._r8
     cpool_livecroot_gr_p(m)                  = 0._r8
     cpool_livecroot_storage_gr_p(m)          = 0._r8
     transfer_livecroot_gr_p(m)               = 0._r8
     cpool_deadcroot_gr_p(m)                  = 0._r8
     cpool_deadcroot_storage_gr_p(m)          = 0._r8
     transfer_deadcroot_gr_p(m)               = 0._r8
     leafc_storage_to_xfer_p(m)               = 0._r8
     frootc_storage_to_xfer_p(m)              = 0._r8
     livestemc_storage_to_xfer_p(m)           = 0._r8
     deadstemc_storage_to_xfer_p(m)           = 0._r8
     livecrootc_storage_to_xfer_p(m)          = 0._r8
     deadcrootc_storage_to_xfer_p(m)          = 0._r8
     gresp_storage_to_xfer_p(m)               = 0._r8
     livestemc_to_deadstemc_p(m)              = 0._r8
     livecrootc_to_deadcrootc_p(m)            = 0._r8
     crop_seedc_to_leaf_p(m)                  = 0._r8

     hrv_xsmrpool_to_atm_p(m)                 = 0._r8 
  
     xsmrpool_to_atm_p(m)                     = 0._r8
     livestemc_to_litter_p(m)                 = 0._r8
     grainc_to_food_p(m)                      = 0._r8
     grainc_to_seed_p(m)                      = 0._r8
     grainc_xfer_to_grainc_p(m)               = 0._r8
     cpool_to_grainc_p(m)                     = 0._r8
     cpool_to_grainc_storage_p(m)             = 0._r8
     cpool_grain_gr_p(m)                      = 0._r8
     cpool_grain_storage_gr_p(m)              = 0._r8
     transfer_grain_gr_p(m)                   = 0._r8
     grainc_storage_to_xfer_p(m)              = 0._r8
  end do

  do j=1,nl_soil
     phenology_to_met_c(j,i)               = 0._r8
     phenology_to_cel_c(j,i)               = 0._r8
     phenology_to_lig_c(j,i)               = 0._r8
   
     gap_mortality_to_met_c(j,i)           = 0._r8
     gap_mortality_to_cel_c(j,i)           = 0._r8
     gap_mortality_to_lig_c(j,i)           = 0._r8
     gap_mortality_to_cwdc(j,i)            = 0._r8
   
     fire_mortality_to_cwdc(j,i)           = 0._r8
     fire_mortality_to_met_c(j,i)          = 0._r8
     fire_mortality_to_cel_c(j,i)          = 0._r8
     fire_mortality_to_lig_c(j,i)          = 0._r8

     do k=1,ndecomp_pools
        m_decomp_cpools_to_fire_vr(j,k,i)= 0._r8
     end do
  end do

  fire_closs(i)          = 0._r8
  fire_nloss(i)                          = 0._r8
  wood_harvestc(i)       = 0._r8
  wood_harvestn(i)                       = 0._r8
  grainc_to_cropprodc(i) = 0._r8
  grainn_to_cropprodn(i) = 0._r8

  do m = ps, pe
     gpp_p(m)                                 = 0._r8
     wood_harvestc_p(m)                       = 0._r8
     wood_harvestn_p(m)                       = 0._r8
     grainc_to_cropprodc_p(m)                 = 0._r8
     grainn_to_cropprodn_p(m)                 = 0._r8
     soyfixn_p(m)                             = 0._r8 
!4386        this%npp_Nactive_patch(i)     = value_patch
!4387        this%npp_burnedoff_patch(i)     = value_patch
!4388        this%npp_Nnonmyc_patch(i)     = value_patch
!4389        this%npp_Nam_patch(i)         = value_patch
!4390        this%npp_Necm_patch(i)        = value_patch
!4391        this%npp_Nactive_no3_patch(i) = value_patch
!4392        this%npp_Nactive_nh4_patch(i) = value_patch
!4393        this%npp_Nnonmyc_no3_patch(i) = value_patch
!4394        this%npp_Nnonmyc_nh4_patch(i) = value_patch
!4395        this%npp_Nam_no3_patch(i)     = value_patch
!4396        this%npp_Nam_nh4_patch(i)     = value_patch
!4397        this%npp_Necm_no3_patch(i)    = value_patch
!4398        this%npp_Necm_nh4_patch(i)    = value_patch
!4399        this%npp_Nfix_patch(i)        = value_patch
!4400        this%npp_Nretrans_patch(i)    = value_patch
!4401        this%npp_Nuptake_patch(i)     = value_patch
!   npp_growth(i)                          = 0._r8
!4403        this%leafc_change_patch(i)    = value_patch
!4404        this%soilc_change_patch(i)    = value_patch

    fire_closs_p(m)                          = 0._r8
    fire_nloss_p(m)                          = 0._r8
        ! Zero p2c column fluxes

! CNVegNitrogenFluxes set zero

     m_leafn_to_litter_p(m)                   = 0._r8
     m_frootn_to_litter_p(m)                  = 0._r8
     m_leafn_storage_to_litter_p(m)           = 0._r8
     m_frootn_storage_to_litter_p(m)          = 0._r8
     m_livestemn_storage_to_litter_p(m)       = 0._r8
     m_deadstemn_storage_to_litter_p(m)       = 0._r8
     m_livecrootn_storage_to_litter_p(m)      = 0._r8
     m_deadcrootn_storage_to_litter_p(m)      = 0._r8
     m_leafn_xfer_to_litter_p(m)              = 0._r8
     m_frootn_xfer_to_litter_p(m)             = 0._r8
     m_livestemn_xfer_to_litter_p(m)          = 0._r8
     m_deadstemn_xfer_to_litter_p(m)          = 0._r8
     m_livecrootn_xfer_to_litter_p(m)         = 0._r8
     m_deadcrootn_xfer_to_litter_p(m)         = 0._r8
     m_livestemn_to_litter_p(m)               = 0._r8
     m_deadstemn_to_litter_p(m)               = 0._r8
     m_livecrootn_to_litter_p(m)              = 0._r8
     m_deadcrootn_to_litter_p(m)              = 0._r8
     m_retransn_to_litter_p(m)                = 0._r8
  
     m_leafn_to_fire_p(m)                     = 0._r8
     m_leafn_storage_to_fire_p(m)             = 0._r8
     m_leafn_xfer_to_fire_p(m)                = 0._r8
     m_livestemn_to_fire_p(m)                 = 0._r8
     m_livestemn_storage_to_fire_p(m)         = 0._r8
     m_livestemn_xfer_to_fire_p(m)            = 0._r8
     m_deadstemn_to_fire_p(m)                 = 0._r8
     m_deadstemn_storage_to_fire_p(m)         = 0._r8
     m_deadstemn_xfer_to_fire_p(m)            = 0._r8
     m_frootn_to_fire_p(m)                    = 0._r8
     m_frootn_storage_to_fire_p(m)            = 0._r8
     m_frootn_xfer_to_fire_p(m)               = 0._r8
     m_livecrootn_to_fire_p(m)                = 0._r8
     m_livecrootn_storage_to_fire_p(m)        = 0._r8
     m_livecrootn_xfer_to_fire_p(m)           = 0._r8
     m_deadcrootn_to_fire_p(m)                = 0._r8
     m_deadcrootn_storage_to_fire_p(m)        = 0._r8
     m_deadcrootn_xfer_to_fire_p(m)           = 0._r8
     m_retransn_to_fire_p(m)                  = 0._r8
     
     
     m_leafn_to_litter_fire_p(m)              = 0._r8
     m_leafn_storage_to_litter_fire_p(m)      = 0._r8
     m_leafn_xfer_to_litter_fire_p(m)         = 0._r8
     m_livestemn_to_litter_fire_p(m)          = 0._r8
     m_livestemn_storage_to_litter_fire_p(m)  = 0._r8
     m_livestemn_xfer_to_litter_fire_p(m)     = 0._r8
     m_livestemn_to_deadstemn_fire_p(m)       = 0._r8
     m_deadstemn_to_litter_fire_p(m)          = 0._r8
     m_deadstemn_storage_to_litter_fire_p(m)  = 0._r8
     m_deadstemn_xfer_to_litter_fire_p(m)     = 0._r8
     m_frootn_to_litter_fire_p(m)             = 0._r8
     m_frootn_storage_to_litter_fire_p(m)     = 0._r8
     m_frootn_xfer_to_litter_fire_p(m)        = 0._r8
     m_livecrootn_to_litter_fire_p(m)         = 0._r8
     m_livecrootn_storage_to_litter_fire_p(m) = 0._r8
     m_livecrootn_xfer_to_litter_fire_p(m)    = 0._r8
     m_livecrootn_to_deadcrootn_fire_p(m)     = 0._r8
     m_deadcrootn_to_litter_fire_p(m)         = 0._r8
     m_deadcrootn_storage_to_litter_fire_p(m) = 0._r8
     m_deadcrootn_xfer_to_litter_fire_p(m)    = 0._r8
     m_retransn_to_litter_fire_p(m)           = 0._r8

     leafn_xfer_to_leafn_p(m)                 = 0._r8
     frootn_xfer_to_frootn_p(m)               = 0._r8
     livestemn_xfer_to_livestemn_p(m)         = 0._r8
     deadstemn_xfer_to_deadstemn_p(m)         = 0._r8
     livecrootn_xfer_to_livecrootn_p(m)       = 0._r8
     deadcrootn_xfer_to_deadcrootn_p(m)       = 0._r8
     leafn_to_litter_p(m)                     = 0._r8
     leafn_to_retransn_p(m)                   = 0._r8
     frootn_to_litter_p(m)                    = 0._r8
     retransn_to_npool_p(m)                   = 0._r8
     free_retransn_to_npool_p(m)              = 0._r8
     sminn_to_npool_p(m)                      = 0._r8
     npool_to_leafn_p(m)                      = 0._r8
     npool_to_leafn_storage_p(m)              = 0._r8
     npool_to_frootn_p(m)                     = 0._r8
     npool_to_frootn_storage_p(m)             = 0._r8
     npool_to_livestemn_p(m)                  = 0._r8
     npool_to_livestemn_storage_p(m)          = 0._r8
     npool_to_deadstemn_p(m)                  = 0._r8
     npool_to_deadstemn_storage_p(m)          = 0._r8
     npool_to_livecrootn_p(m)                 = 0._r8
     npool_to_livecrootn_storage_p(m)         = 0._r8
     npool_to_deadcrootn_p(m)                 = 0._r8
     npool_to_deadcrootn_storage_p(m)         = 0._r8
     leafn_storage_to_xfer_p(m)               = 0._r8
     frootn_storage_to_xfer_p(m)              = 0._r8
     livestemn_storage_to_xfer_p(m)           = 0._r8
     deadstemn_storage_to_xfer_p(m)           = 0._r8
     livecrootn_storage_to_xfer_p(m)          = 0._r8
     deadcrootn_storage_to_xfer_p(m)          = 0._r8
     livestemn_to_deadstemn_p(m)              = 0._r8
     livestemn_to_retransn_p(m)               = 0._r8
     livecrootn_to_deadcrootn_p(m)            = 0._r8
     livecrootn_to_retransn_p(m)              = 0._r8
   
     crop_seedn_to_leaf_p(m)                  = 0._r8

     livestemn_to_litter_p(m)              = 0._r8
     grainn_to_food_p(m)                   = 0._r8
     grainn_to_seed_p(m)                   = 0._r8
     grainn_xfer_to_grainn_p(m)            = 0._r8
     npool_to_grainn_p(m)                  = 0._r8
     npool_to_grainn_storage_p(m)          = 0._r8
     grainn_storage_to_xfer_p(m)           = 0._r8
   !soyfixn(i)                          = 0._r8
     frootn_to_retransn_p(m)               = 0._r8
  end do

  do j=1,nl_soil

      ! phenology: litterfall and crop fluxes associated wit
     phenology_to_met_n(j,i)        = 0._r8
     phenology_to_cel_n(j,i)        = 0._r8
     phenology_to_lig_n(j,i)        = 0._r8

      ! gap mortality
     gap_mortality_to_met_n(j,i)    = 0._r8
     gap_mortality_to_cel_n(j,i)    = 0._r8
     gap_mortality_to_lig_n(j,i)    = 0._r8
     gap_mortality_to_cwdn(j,i)          = 0._r8

      ! fire
     fire_mortality_to_cwdn(j,i)         = 0._r8
     fire_mortality_to_met_n(j,i)         = 0._r8
     fire_mortality_to_cel_n(j,i)         = 0._r8
     fire_mortality_to_lig_n(j,i)         = 0._r8

!      ! harvest
!      this%harvest_n_to_litr_met_n(j,i)          = 0._r8
!      this%harvest_n_to_litr_cel_n(j,i)          = 0._r8
!      this%harvest_n_to_litr_lig_n(j,i)          = 0._r8
!      this%harvest_n_to_cwdn(j,i)                = 0._r8
  end do
         
  do k = 1, ndecomp_pools
     do j = 1, nl_soil
        m_decomp_npools_to_fire_vr(j,k,i) = 0._r8
     end do
  end do

  do k=1,ndecomp_transitions
     do j=1,nl_soil
!         decomp_cascade_hr(k,i)             = 0._r8
        decomp_hr_vr(j,k,i)        = 0._r8
!         decomp_cascade_ctransfer(k,i)      = 0._r8
        decomp_ctransfer_vr(j,k,i) = 0._r8
     end do
  end do

  do k = 1, ndecomp_pools
!      decomp_cpools_leached(k,i) = 0._r8
     do j = 1, nl_soil
        decomp_cpools_transport_tendency(j,k,i) = 0._r8
        decomp_cpools_sourcesink(j,k,i)         = 0._r8
        decomp_k(j,k,i)                         = 0._r8
     end do
  end do

!   do j = 1, nl_soil
!      hr_vr(j,i) = 0._r8
!   end do

!   hr(i)            = 0._r8
   somc_fire(i)     = 0._r8
   som_c_leached(i) = 0._r8
!   somhr(i)         = 0._r8
!   lithr(i)         = 0._r8
!   soilc_change(i)  = 0._r8
   decomp_hr(i)     = 0._r8


   do j = 1, nl_soil
#ifndef NITRIF
      sminn_to_denit_excess_vr(j,i)      = 0._r8
      sminn_leached_vr(j,i)              = 0._r8
      sminn_to_plant_fun_vr(j,i)         = 0._r8
#else
      f_nit_vr(j,i)                      = 0._r8
      f_denit_vr(j,i)                    = 0._r8
      smin_no3_leached_vr(j,i)           = 0._r8
      smin_no3_runoff_vr(j,i)            = 0._r8
      n2_n2o_ratio_denit_vr(j,i)         = 0._r8
      pot_f_nit_vr(j,i)                  = 0._r8
      pot_f_denit_vr(j,i)                = 0._r8
      actual_immob_no3_vr(j,i)           = 0._r8
      actual_immob_nh4_vr(j,i)           = 0._r8
      smin_no3_to_plant_vr(j,i)          = 0._r8
      smin_nh4_to_plant_vr(j,i)          = 0._r8
      f_n2o_denit_vr(j,i)                = 0._r8
      f_n2o_nit_vr(j,i)                  = 0._r8

#endif
      potential_immob_vr(j,i)               = 0._r8
      actual_immob_vr(j,i)                  = 0._r8
      sminn_to_plant(i)                = 0._r8
      sminn_to_plant_vr(j,i)                = 0._r8
      supplement_to_sminn_vr(j,i)           = 0._r8
      gross_nmin_vr(j,i)                    = 0._r8
      net_nmin_vr(j,i)                      = 0._r8
      sminn_to_plant_fun_no3_vr(j,i)        = 0._r8
      sminn_to_plant_fun_nh4_vr(j,i)        = 0._r8
   end do

   nfix_to_sminn(i)             = 0._r8
   ffix_to_sminn(i)             = 0._r8
   fert_to_sminn(i)             = 0._r8
    soyfixn_to_sminn(i)          = 0._r8
!    potential_immob(i)           = 0._r8
!    actual_immob(i)              = 0._r8
!    sminn_to_plant(i)            = 0._r8
   supplement_to_sminn(i)       = 0._r8
   gross_nmin(i)                = 0._r8
   net_nmin(i)                  = 0._r8
    denit(i)                     = 0._r8
!    sminn_to_plant_fun(i)        = 0._r8
!#ifdef NITRIF
!    f_nit(i)                  = 0._r8
!    pot_f_nit(i)              = 0._r8
!    f_denit(i)                = 0._r8
!    pot_f_denit(i)            = 0._r8
!    f_n2o_denit(i)            = 0._r8
    f_n2o_nit(i)              = 0._r8
    smin_no3_leached(i)       = 0._r8
    smin_no3_runoff(i)        = 0._r8
!#else
!    sminn_to_denit_excess(i)  = 0._r8
    sminn_leached(i)          = 0._r8
!#endif
!    ninputs(i)                   = 0._r8
!    noutputs(i)                  = 0._r8
    som_n_leached(i)             = 0._r8

!    do k = 1, ndecomp_pools
!       decomp_npools_leached(i,k) = 0._r8
!    end do
   
!    if(use_soil_matrixcn)then
!       call this%matrix_Ninput%SetValueV_scaler(num_column,filter_column(1:num_column),value_column)
!    end if

   do k = 1, ndecomp_pools
      do j = 1, nl_soil
         decomp_npools_transport_tendency(j,k,i) = 0._r8
      end do
   end do

!    do l = 1, ndecomp_cascade_transitions
!       decomp_ntransfer(i,l) = 0._r8
!       decomp_sminn_flux(i,l) = 0._r8
#ifndef NITRIF
!       sminn_to_denit_decomp_cascade(i,l) = 0._r8
#endif
!    end do

   do k = 1, ndecomp_transitions
      do j = 1, nl_soil
         decomp_ntransfer_vr(j,k,i) = 0._r8
         decomp_sminn_flux_vr(j,k,i) = 0._r8
#ifndef NITRIF
         sminn_to_denit_decomp_vr(j,k,i) = 0._r8
#endif
      end do
   end do

   do k = 1, ndecomp_pools
      do j = 1, nl_soil
         decomp_npools_sourcesink(j,k,i) = 0._r8
      end do
   end do


end subroutine CNZeroFluxes

end module bgc_CNZeroFluxesMod
#endif
