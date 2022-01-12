#include <define.h>

SUBROUTINE bgc_driver &
          (i,idate,deltim,dlat,dlon)!, woody, &

  use precision
  use PhysicalConstants, only : tfrz, denh2o, denice
  use MOD_PFTimeInvars, only: pftfrac, patch_pft_s, patch_pft_e
  use MOD_1D_Fluxes, only: plant_ndemand, ndep_to_sminn
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
  use bgc_CNBalanceCheckMod, only: BeginCNBalance, CBalanceCheck, NBalanceCheck
  use timemanager
  use GlobalVars, only: nl_soil, nl_soil_full, ndecomp_pools, ndecomp_transitions, npcropmin, &
                      z_soi,dz_soi,zi_soi,nbedrock,zmin_bedrock

  use MOD_TimeVariables, only: sminn_vr, col_begnb
  implicit none

  integer ,intent(in) :: i
  real(r8),intent(in) :: deltim
  integer ,intent(in) :: idate(3)
  real(r8),intent(in) :: dlat
  real(r8),intent(in) :: dlon

  integer :: ps, pe
  call BeginCNBalance(i)
!  if(i .eq. 79738)print*,'ndep in bgcdriver',ndep_to_sminn(i) 
!  if(i .eq. 79738)print*,'col_begnb after BeginCNBalance',col_begnb(i)
  ps = patch_pft_s(i)      
  pe = patch_pft_e(i)
! update vegetation pools from phenology, allocation and nitrogen uptake
! update soil pools from decomposition and nitrogen competition
  call CNZeroFluxes(i, ps, pe, nl_soil, ndecomp_pools, ndecomp_transitions)
!  if(i .eq. 79738)print*,'ndep after 0flux',ndep_to_sminn(i) 
  call CNMResp(i, ps, pe, nl_soil, npcropmin)
  call decomp_rate_constants_bgc(i,nl_soil,z_soi)
  call SoilBiogeochemPotential(i,nl_soil,ndecomp_pools,ndecomp_transitions)
!  if(i .eq. 79738)print*,'ndep after soil pontential',ndep_to_sminn(i) 
  call SoilBiogeochemVerticalProfile(i,ps,pe,nl_soil,nl_soil_full,nbedrock,zmin_bedrock,z_soi,dz_soi)
  call calc_plant_nutrient_demand_CLM45_default(i,ps,pe,deltim,npcropmin)
  
  plant_ndemand(i) = sum( plant_ndemand_p(ps:pe)*pftfrac(ps:pe) )

  call SoilBiogeochemCompetition(i,deltim,nl_soil,dz_soi)
!  if(i .eq. 79738)print*,'ndep after soil competition',ndep_to_sminn(i) 
  call calc_plant_nutrient_competition_CLM45_default(i,ps,pe,npcropmin)
  call SoilBiogeochemDecomp(i,nl_soil,ndecomp_pools,ndecomp_transitions, dz_soi)
  call CNPhenology(i,ps,pe,nl_soil,idate(1:3),deltim,dlat,npcropmin,phase=1)
  call CNPhenology(i,ps,pe,nl_soil,idate(1:3),deltim,dlat,npcropmin,phase=2)
  call CNGResp(i, ps, pe, npcropmin)
  call CStateUpdate0()
  call CStateUpdate1(i, ps, pe, deltim, nl_soil, ndecomp_transitions, npcropmin)
  call NStateUpdate1(i, ps, pe, deltim, nl_soil, ndecomp_transitions, npcropmin,dz_soi)
!  if(i .eq. 79738)print*,'ndep before SoilBiogeochemNStateUpdate1',ndep_to_sminn(i)*deltim, &
!                  sum(sminn_vr(1:nl_soil,i)*dz_soi(1:nl_soil))
  call SoilBiogeochemNStateUpdate1(i,deltim,nl_soil,ndecomp_transitions,dz_soi)
!  if(i .eq. 147958)print*,'ndep after SoilBiogeochemNStateUpdate1',sum(sminn_vr(1:nl_soil,i)*dz_soi(1:nl_soil))

! update soil pools from soil vertical mixing
  call SoilBiogeochemLittVertTransp(i,deltim,nl_soil,nl_soil_full,ndecomp_pools,nbedrock,z_soi,zi_soi,dz_soi)

! update vegetation pools from gap mortality
  call CNGapMortality(i, ps, pe, nl_soil,npcropmin)
  call CStateUpdate2(i, ps, pe, deltim, nl_soil)
  call NStateUpdate2(i, ps, pe, deltim, nl_soil, dz_soi)

!  if(i .eq. 79738)print*,'col_begnb after NStateUpdate2',col_begnb(i)

#ifdef FIRE
! update vegetation and fire pools from fire
  call CNFireArea(i,ps,pe,dlat,nl_soil,idate,dz_soi)
  call CNFireFluxes(i,ps,pe,dlat,nl_soil,ndecomp_pools)
#endif
  call CStateUpdate3(i,ps,pe,deltim, nl_soil, ndecomp_pools)
  call CNAnnualUpdate(i,ps,pe,deltim,idate(1:3))

! update soil mineral nitrogen pools leaching
  call SoilBiogeochemNLeaching(i,deltim,nl_soil,zi_soi,dz_soi)
!if(i .eq. 147958)print*,'sminn after NStateUpdate2',sum(sminn_vr(1:nl_soil,i)*dz_soi(1:nl_soil))
  call NstateUpdate3(i, ps, pe, deltim, nl_soil, ndecomp_pools,dz_soi)
!if(i .eq. 147958)print*,'sminn after NStateUpdate3',sum(sminn_vr(1:nl_soil,i)*dz_soi(1:nl_soil))
  call CNDriverSummarizeStates(i,ps,pe,nl_soil,dz_soi,ndecomp_pools)
  call CNDriverSummarizeFluxes(i,ps,pe,nl_soil,dz_soi,ndecomp_transitions,ndecomp_pools)
  call CBalanceCheck(i,deltim,dlat,dlon)
!  if(i .eq. 79738)print*,'col_begnb before NBalanceCheck',col_begnb(i)
!  if(i .eq. 79738)print*,'ndep before Nbalance check',ndep_to_sminn(i)*deltim, &
!                  sum(sminn_vr(1:nl_soil,i)*dz_soi(1:nl_soil))
  call NBalanceCheck(i,deltim,dlat,dlon)

END SUBROUTINE bgc_driver


