#include <define.h>

SUBROUTINE bgc_driver &
          (i,ps,pe,idate,deltim,dlat)!, woody, &

  use precision
  use PhysicalConstants, only : tfrz, denh2o, denice
  use MOD_PFTimeInvars, only: pftfrac
  use MOD_1D_Fluxes, only: plant_ndemand
  use MOD_1D_PFTFluxes, only: plant_ndemand_p
  use bgc_veg_CNMRespMod, only: CNMResp
  use bgc_soil_SoilBiogeochemDecompCascadeBGCMod, only: decomp_rate_constants_bgc
  use bgc_soil_SoilBiogeochemPotentialMod, only: SoilBiogeochemPotential
  use bgc_soil_SoilBiogeochemVerticalProfileMod, only: SoilBiogeochemVerticalProfile
  use bgc_veg_NutrientCompetitionMod, only: calc_plant_nutrient_demand_CLM45_default,&
                                            calc_plant_nutrient_competition_CLM45_default
  use bgc_soil_SoilBiogeochemCompetitionMod, only: SoilBiogeochemCompetition
  use bgc_soil_SoilBiogeochemDecompMod, only: SoilBiogeochemDecomp
  use bgc_veg_CNPhenologyMod, only: CNPhenology
  use bgc_veg_CNGRespMod, only: CNGResp
  use bgc_CNCStateUpdate1Mod, only: CStateUpdate0, CStateUpdate1
  use bgc_CNNStateUpdate1Mod, only: NStateUpdate1
  use bgc_SoilBiogeochemNStateUpdate1Mod, only: SoilBiogeochemNStateUpdate1
  use bgc_soil_SoilBiogeochemLittVertTranspMod, only: SoilBiogeochemLittVertTransp
  use bgc_veg_CNGapMortalityMod, only: CNGapMortality
  use bgc_CNCStateUpdate2Mod, only: CStateUpdate2
  use bgc_CNNStateUpdate2Mod, only: NStateUpdate2
!  use bgc_veg_CNFireLi2016Mod, only: CNFireArea
!  use bgc_veg_CNFireBaseMod, only: CNFireFluxes
  use bgc_CNCStateUpdate3Mod, only: CStateUpdate3
  use bgc_soil_SoilBiogeochemNLeachingMod, only: SoilBiogeochemNLeaching
  use bgc_CNNStateUpdate3Mod, only: NstateUpdate3
  use bgc_CNSummaryMod, only: CNDriverSummarizeStates, CNDriverSummarizeFluxes
  use bgc_CNAnnualUpdateMod, only: CNAnnualUpdate
  use bgc_CNZeroFluxesMod, only: CNZeroFluxes
  use timemanager
  use GlobalVars, only: nl_soil, nl_soil_full, ndecomp_pools, ndecomp_transitions, npcropmin, &
                      z_soi,dz_soi,zi_soi,nbedrock,zmin_bedrock

  implicit none

  integer ,intent(in) :: i
  integer ,intent(in) :: ps
  integer ,intent(in) :: pe
  real(r8),intent(in) :: deltim
  integer ,intent(in) :: idate(3)
  real(r8),intent(in) :: dlat

! update vegetation pools from phenology, allocation and nitrogen uptake
! update soil pools from decomposition and nitrogen competition
  call CNZeroFluxes(i, ps, pe, nl_soil, ndecomp_pools, ndecomp_transitions)
  call CNMResp(i, ps, pe, nl_soil, npcropmin)
  call decomp_rate_constants_bgc(i,nl_soil,z_soi)
  call SoilBiogeochemPotential(i,nl_soil,ndecomp_pools,ndecomp_transitions)
  call SoilBiogeochemVerticalProfile(i,ps,pe,nl_soil,nl_soil_full,nbedrock,zmin_bedrock,z_soi,dz_soi)
  call calc_plant_nutrient_demand_CLM45_default(i,ps,pe,deltim,npcropmin)
  
  plant_ndemand(i) = sum( plant_ndemand_p(ps:pe)*pftfrac(ps:pe) )

  call SoilBiogeochemCompetition(i,deltim,nl_soil,dz_soi)
  call calc_plant_nutrient_competition_CLM45_default(i,ps,pe,npcropmin)
  call SoilBiogeochemDecomp(i,nl_soil,ndecomp_pools,ndecomp_transitions, dz_soi)
  call CNPhenology(i,ps,pe,nl_soil,idate(1:3),deltim,dlat,npcropmin,phase=1)
  call CNPhenology(i,ps,pe,nl_soil,idate(1:3),deltim,dlat,npcropmin,phase=2)
  call CNGResp(i, ps, pe, npcropmin)
  call CStateUpdate0()
  call CStateUpdate1(i, ps, pe, deltim, nl_soil, ndecomp_transitions, npcropmin)
  call NStateUpdate1(i, ps, pe, deltim, nl_soil, ndecomp_transitions, npcropmin)
  call SoilBiogeochemNStateUpdate1(i,deltim,nl_soil,ndecomp_transitions)

! update soil pools from soil vertical mixing
  call SoilBiogeochemLittVertTransp(i,deltim,nl_soil,nl_soil_full,ndecomp_pools,nbedrock,z_soi,zi_soi,dz_soi)

! update vegetation pools from gap mortality
  call CNGapMortality(i, ps, pe, nl_soil,npcropmin)
  call CStateUpdate2(i, ps, pe, deltim, nl_soil)
  call NStateUpdate2(i, ps, pe, deltim, nl_soil)


#ifdef FIRE
! update vegetation and fire pools from fire
  call CNFireArea(i,ps,pe,dlat,nl_soil,idate,dz_soi)
  call CNFireFluxes(i,ps,pe,dlat,nl_soil,ndecomp_pools)
#endif
  call CStateUpdate3(i,ps,pe,deltim, nl_soil, ndecomp_pools)
  call CNAnnualUpdate(i,ps,pe,deltim,idate(1:3))

! update soil mineral nitrogen pools leaching
  call SoilBiogeochemNLeaching(i,deltim,nl_soil,zi_soi,dz_soi)
  call NstateUpdate3(i, ps, pe, deltim, nl_soil, ndecomp_pools)
  call CNDriverSummarizeStates()
  call CNDriverSummarizeFluxes(i,ps,pe,nl_soil,dz_soi)

END SUBROUTINE bgc_driver


