#include <define.h>
#ifdef BGC
#include <define.h>

  SUBROUTINE bgc_driver &
          (i,idate,deltim,dlat,dlon)

    !-----------------------------------------------------------------------------------------------------------
    ! !DESCRIPTION:
    ! The trunk subroutine of the CoLM biogeochemistry module. The bgc_driver link different bgc processes, and
    ! sequentially run each process step by step. bgc_driver was called by CoLMDRIVER includes vegetation
    ! and soil CN cycle processes.
    ! 
    !
    ! !ORIGINAL:
    ! The Community Land Model version 5.0 (CLM5.0)
    !
    ! !REFERENCE:
    ! Lawrence, D.M., Fisher, R.A., Koven, C.D., Oleson, K.W., Swenson, S.C., Bonan, G., Collier, N., 
    ! Ghimire, B., van Kampenhout, L., Kennedy, D. and Kluzek, E., 2019. 
    ! The Community Land Model version 5: Description of new features, benchmarking,
    ! and impact of forcing uncertainty. Journal of Advances in Modeling Earth Systems, 11(12), 4245-4287.
    !
    ! !REVISION:
    ! Xingjie Lu, 2022, modify original CLM5 to be compatible with CoLM code structure.


    use MOD_Precision
    use MOD_Namelist, only : DEF_USE_SASU, DEF_USE_NITRIF, DEF_USE_CNSOYFIXN, DEF_USE_FIRE
    use MOD_Const_Physical, only : tfrz, denh2o, denice
    use MOD_Vars_PFTimeInvariants, only: pftfrac
    use MOD_LandPFT, only: patch_pft_s, patch_pft_e
    use MOD_BGC_Vars_1DFluxes, only: plant_ndemand, ndep_to_sminn
    use MOD_BGC_Vars_1DPFTFluxes, only: plant_ndemand_p, cpool_to_leafc_p, crop_seedc_to_leaf_p
    use MOD_BGC_Veg_CNMResp, only: CNMResp
    use MOD_BGC_Soil_BiogeochemDecompCascadeBGC, only: decomp_rate_constants_bgc
    use MOD_BGC_Soil_BiogeochemPotential, only: SoilBiogeochemPotential
    use MOD_BGC_Soil_BiogeochemVerticalProfile, only: SoilBiogeochemVerticalProfile
    use MOD_BGC_Veg_NutrientCompetition, only: calc_plant_nutrient_demand_CLM45_default,&
                                            calc_plant_nutrient_competition_CLM45_default
    use MOD_BGC_Soil_BiogeochemNitrifDenitrif, only: SoilBiogeochemNitrifDenitrif

    use MOD_BGC_Soil_BiogeochemCompetition, only: SoilBiogeochemCompetition
    use MOD_BGC_Soil_BiogeochemDecomp, only: SoilBiogeochemDecomp
    use MOD_BGC_Veg_CNPhenology, only: CNPhenology
    use MOD_BGC_Veg_CNGResp, only: CNGResp
    use MOD_BGC_CNCStateUpdate1, only: CStateUpdate1
    use MOD_BGC_CNNStateUpdate1, only: NStateUpdate1
    use MOD_BGC_Soil_BiogeochemNStateUpdate1, only: SoilBiogeochemNStateUpdate1
    use MOD_BGC_Soil_BiogeochemLittVertTransp, only: SoilBiogeochemLittVertTransp
    use MOD_BGC_Veg_CNGapMortality, only: CNGapMortality
    use MOD_BGC_CNCStateUpdate2, only: CStateUpdate2
    use MOD_BGC_CNNStateUpdate2, only: NStateUpdate2
    use MOD_BGC_CNCStateUpdate3, only: CStateUpdate3
    use MOD_BGC_Soil_BiogeochemNLeaching, only: SoilBiogeochemNLeaching
    use MOD_BGC_CNNStateUpdate3, only: NstateUpdate3
    use MOD_BGC_CNSummary, only: CNDriverSummarizeStates, CNDriverSummarizeFluxes
    use MOD_BGC_CNAnnualUpdate, only: CNAnnualUpdate
    use MOD_BGC_CNZeroFluxes, only: CNZeroFluxes
    use MOD_BGC_Veg_CNVegStructUpdate, only: CNVegStructUpdate
    use MOD_BGC_CNBalanceCheck, only: BeginCNBalance, CBalanceCheck, NBalanceCheck
    use MOD_BGC_CNSASU, only: CNSASU
    use MOD_BGC_Veg_CNNDynamics, only: CNNFixation
#ifdef CROP
    use MOD_BGC_Veg_CNNDynamics, only: CNNFert, CNSoyfix
#endif
    use MOD_TimeManager
    use MOD_Vars_Global, only: nl_soil, nl_soil_full, ndecomp_pools, ndecomp_pools_vr, ndecomp_transitions, npcropmin, &
                        z_soi,dz_soi,zi_soi,nbedrock,zmin_bedrock

    use MOD_BGC_Vars_TimeVariables, only: sminn_vr, col_begnb, skip_balance_check, decomp_cpools_vr
    use MOD_BGC_Veg_CNFireBase, only: CNFireFluxes
    use MOD_BGC_Veg_CNFireLi2016, only: CNFireArea

    implicit none

    integer ,intent(in) :: i        ! patch index
    real(r8),intent(in) :: deltim   ! time step in seconds
    integer ,intent(in) :: idate(3) ! current date (year, day of the year, seconds of the day)
    real(r8),intent(in) :: dlat     ! latitude (degrees)
    real(r8),intent(in) :: dlon     ! longitude (degrees)

    integer :: ps, pe
    integer j
    ps = patch_pft_s(i)      
    pe = patch_pft_e(i)
    call BeginCNBalance(i)
    call CNZeroFluxes(i, ps, pe, nl_soil, ndecomp_pools, ndecomp_transitions)
    call CNNFixation(i,idate)
    call CNMResp(i, ps, pe, nl_soil, npcropmin)
    call decomp_rate_constants_bgc(i,nl_soil,z_soi)
    call SoilBiogeochemPotential(i,nl_soil,ndecomp_pools,ndecomp_transitions)
    call SoilBiogeochemVerticalProfile(i,ps,pe,nl_soil,nl_soil_full,nbedrock,zmin_bedrock,z_soi,dz_soi)
    if(DEF_USE_NITRIF)then
       call SoilBiogeochemNitrifDenitrif(i,nl_soil,dz_soi)
    end if
    call calc_plant_nutrient_demand_CLM45_default(i,ps,pe,deltim,npcropmin)
  
    plant_ndemand(i) = sum( plant_ndemand_p(ps:pe)*pftfrac(ps:pe) )

    call SoilBiogeochemCompetition(i,deltim,nl_soil,dz_soi)
    call calc_plant_nutrient_competition_CLM45_default(i,ps,pe,npcropmin)
#ifdef CROP
    if(DEF_USE_CNSOYFIXN)then
       call CNSoyfix (i, ps, pe, nl_soil)
    end if
#endif
  
    call SoilBiogeochemDecomp(i,nl_soil,ndecomp_pools,ndecomp_transitions, dz_soi)
    call CNPhenology(i,ps,pe,nl_soil,idate(1:3),dz_soi,deltim,dlat,npcropmin,phase=1)
    call CNPhenology(i,ps,pe,nl_soil,idate(1:3),dz_soi,deltim,dlat,npcropmin,phase=2)
#ifdef CROP
    call CNNFert(i, ps, pe)
#endif
    call CNGResp(i, ps, pe, npcropmin)

    ! update vegetation pools from phenology, allocation and nitrogen uptake
    ! update soil N pools from decomposition and nitrogen competition
    call CStateUpdate1(i, ps, pe, deltim, nl_soil, ndecomp_transitions, npcropmin)
    call NStateUpdate1(i, ps, pe, deltim, nl_soil, ndecomp_transitions, npcropmin,dz_soi)
    call SoilBiogeochemNStateUpdate1(i,deltim,nl_soil,ndecomp_transitions,dz_soi)

    ! update soil pools from soil vertical mixing
    call SoilBiogeochemLittVertTransp(i,deltim,nl_soil,nl_soil_full,ndecomp_pools,nbedrock,z_soi,zi_soi,dz_soi)

    ! update vegetation pools from gap mortality
    call CNGapMortality(i, ps, pe, nl_soil,npcropmin)
    call CStateUpdate2(i, ps, pe, deltim, nl_soil)
    call NStateUpdate2(i, ps, pe, deltim, nl_soil, dz_soi)


    if(DEF_USE_FIRE)then
       ! update vegetation and fire pools from fire
       call CNFireArea(i,ps,pe,dlat,nl_soil,idate,dz_soi)
       call CNFireFluxes(i,ps,pe,dlat,nl_soil,ndecomp_pools)
    end if   
    call CStateUpdate3(i,ps,pe,deltim, nl_soil, ndecomp_pools)
    call CNAnnualUpdate(i,ps,pe,deltim,idate(1:3))

! update soil mineral nitrogen pools leaching
    call SoilBiogeochemNLeaching(i,deltim,nl_soil,zi_soi,dz_soi)
    call NstateUpdate3(i, ps, pe, deltim, nl_soil, ndecomp_pools,dz_soi)

    if(DEF_USE_SASU)then
       call CNSASU(i,ps,pe,deltim,idate(1:3),nl_soil,ndecomp_transitions,ndecomp_pools,ndecomp_pools_vr)! only for spin up
    end if

    call CNDriverSummarizeStates(i,ps,pe,nl_soil,dz_soi,ndecomp_pools)
    call CNDriverSummarizeFluxes(i,ps,pe,nl_soil,dz_soi,ndecomp_transitions,ndecomp_pools,deltim)

    if( .not. skip_balance_check(i) )then
      call CBalanceCheck(i,ps,pe,deltim,dlat,dlon)
      call NBalanceCheck(i,deltim,dlat,dlon)


    else
      skip_balance_check(i) = .false.
    end if

    call CNVegStructUpdate(i,ps,pe,deltim,npcropmin)

  END SUBROUTINE bgc_driver

#endif
