#include <define.h>
#ifdef BGC
module MOD_BGC_Veg_NutrientCompetition

  !----------------------------------------------------------------------------------------------------
  ! !DESCRIPTION
  ! This module simulates the plant growth with regard to the available soil mineral nitrogen.
  ! Allocation of NPP and N uptake to different vegetation CN pools uses allocation scheme from CLM4.5.
  ! Call sequence is: calc_plant_nutrient_demand_CLM45_default => calc_plant_nutrient_competition_CLM45_default
  ! 
  ! !ORIGINAL:
  ! The Community Land Model version 5.0 (CLM5.0)

  ! !REVISION:
  ! Xingjie Lu, 2022, modify original CLM5 to be compatible with CoLM code structure. 
  ! Fang Li, 2022, add GPAM C allocation scheme for crop.

  ! 
    use MOD_Precision
    use MOD_Vars_PFTConst, only: &
        woody, leafcn, frootcn, livewdcn, deadwdcn, graincn, &
        froot_leaf, croot_stem, stem_leaf, flivewd, grperc, grpnow, fcur2, &
! crop variables
        astemf, arooti, arootf, fleafi, bfact, declfact, allconss, allconsl, fleafcn, fstemcn, ffrootcn, &
        lfemerg, grnfill
    
    use MOD_Vars_PFTimeInvars, only: pftclass, pftfrac

    use MOD_BGC_Vars_PFTimeVars, only: &
        xsmrpool_p, retransn_p, &
        tempsum_potential_gpp_p, tempmax_retransn_p, annmax_retransn_p, annsum_potential_gpp_p, &
! crop variables
#ifdef CROP
        croplive_p, hui_p,  peaklai_p, &
        aroot_p, astem_p, arepr_p, aleaf_p, astemi_p, aleafi_p, vf_p, &
#endif
        c_allometry_p, n_allometry_p, downreg_p, grain_flag_p, annsum_npp_p, &
        leafc_p, livestemc_p, frootc_p
    use GlobalVars, only: nwwheat, nirrig_wwheat 

    use MOD_BGC_Vars_TimeVars, only: fpg 
    use GlobalVars, only: ntmp_soybean, ntrp_soybean, nirrig_tmp_soybean, nirrig_trp_soybean

    use MOD_Vars_1DPFTFluxes, only: assim_p

    use MOD_BGC_Vars_1DPFTFluxes, only: &
        leaf_xsmr_p, froot_xsmr_p, livestem_xsmr_p, livecroot_xsmr_p, grain_xsmr_p, cpool_to_xsmrpool_p, &
        leaf_mr_p, froot_mr_p, livestem_mr_p, livecroot_mr_p, grain_mr_p, &
        plant_ndemand_p, retransn_to_npool_p, cpool_to_leafc_p, cpool_to_leafc_storage_p, &
        cpool_to_frootc_p, cpool_to_frootc_storage_p, cpool_to_livestemc_p, cpool_to_livestemc_storage_p, &
        cpool_to_deadstemc_p, cpool_to_deadstemc_storage_p, cpool_to_livecrootc_p, cpool_to_livecrootc_storage_p, &
        cpool_to_deadcrootc_p, cpool_to_deadcrootc_storage_p, cpool_to_grainc_p, cpool_to_grainc_storage_p, &
        cpool_to_gresp_storage_p, npool_to_leafn_p, npool_to_leafn_storage_p, &
        npool_to_frootn_p, npool_to_frootn_storage_p, npool_to_livestemn_p, npool_to_livestemn_storage_p, &
        npool_to_deadstemn_p, npool_to_deadstemn_storage_p, npool_to_livecrootn_p, npool_to_livecrootn_storage_p, &
        npool_to_deadcrootn_p, npool_to_deadcrootn_storage_p, npool_to_grainn_p, npool_to_grainn_storage_p, &
        leafn_to_retransn_p, livestemn_to_retransn_p, frootn_to_retransn_p, &
        plant_calloc_p, plant_nalloc_p, leaf_curmr_p, froot_curmr_p, livestem_curmr_p, livecroot_curmr_p, grain_curmr_p, &
        psn_to_cpool_p, gpp_p, availc_p, avail_retransn_p, xsmrpool_recover_p, sminn_to_npool_p, excess_cflux_p

implicit none

public calc_plant_nutrient_competition_CLM45_default
public calc_plant_nutrient_demand_CLM45_default

contains

 subroutine calc_plant_nutrient_competition_CLM45_default(i,ps,pe,npcropmin)

    !----------------------------------------------------------------------------
    ! !DESCRIPTION
    ! Calulate the nitrogen limitation on the plant growth based on the available 
    ! nitrogen and nitrogen demand from "calc_plant_nutrient_demand_CLM45_default".
    !
    ! !Original:
    ! The Community Land Model version 5.0 (CLM5.0)

    ! !REVISION:
    ! Xingjie Lu, 2022, modify original CLM5 to be compatible with CoLM code structure. 
    ! Fang Li, 2022, add GPAM C allocation scheme for crop.

    integer ,intent(in) :: i
    integer ,intent(in) :: ps
    integer ,intent(in) :: pe
    integer ,intent(in) :: npcropmin

    ! !LOCAL VARIABLES:
    real(r8):: f1,f2,f3,f4,g1,g2  ! allocation parameters
    real(r8):: cnl,cnfr,cnlw,cndw ! C:N ratios for leaf, fine root, and wood
    real(r8):: fcur               ! fraction of current psn displayed as growth
    real(r8):: gresp_storage      ! temporary variable for growth resp to storage
    real(r8):: nlc                ! temporary variable for total new leaf carbon allocation
    real(r8):: f5                 ! grain allocation parameter
    real(r8):: cng                ! C:N ratio for grain (= cnlw for now; slevis)
    integer :: ivt, m

      do m = ps, pe
         ivt = pftclass(m)
         ! set some local allocation variables
         f1 = froot_leaf(ivt)
         f2 = croot_stem(ivt)

         ! modified wood allocation to be 2.2 at npp=800 gC/m2/yr, 0.2 at npp=0,
         ! constrained so that it does not go lower than 0.2 (under negative annsum_npp)
         ! There was an error in this formula in previous version, where the coefficient
         ! was 0.004 instead of 0.0025.
         ! This variable allocation is only for trees. Shrubs have a constant
         ! allocation as specified in the pft-physiology file.  The value is also used
         ! as a trigger here: -1.0 means to use the dynamic allocation (trees).
         if (stem_leaf(ivt) == -1._r8) then
            f3 = (2.7/(1.0+exp(-0.004*(annsum_npp_p(m) - 300.0)))) - 0.4
         else
            f3 = stem_leaf(ivt)
         end if

         f4   = flivewd(ivt)
         g1   = grperc(ivt)
         g2   = grpnow(ivt)
         cnl  = leafcn(ivt)
         cnfr = frootcn(ivt)
         cnlw = livewdcn(ivt)
         cndw = deadwdcn(ivt)
         fcur = fcur2(ivt)

#ifdef CROP
         if (ivt >= npcropmin) then ! skip 2 generic crops
            if (croplive_p(m)) then
               f1 = aroot_p(m) / aleaf_p(m)
               f3 = astem_p(m) / aleaf_p(m)
               f5 = arepr_p(m) / aleaf_p(m)
               g1 = grperc(ivt)
            else
               f1 = 0._r8
               f3 = 0._r8
               f5 = 0._r8
               g1 = grperc(ivt)
            end if
         end if
#endif

         sminn_to_npool_p(m) = plant_ndemand_p(m) * fpg(i)        

         plant_nalloc_p(m) = sminn_to_npool_p(m) + retransn_to_npool_p(m)
         plant_calloc_p(m) = plant_nalloc_p(m) * (c_allometry_p(m)/n_allometry_p(m))


         excess_cflux_p(m) = availc_p(m) - plant_calloc_p(m)
	    ! reduce gpp fluxes due to N limitation
         if (gpp_p(m) > 0.0_r8) then
	    downreg_p(m) = excess_cflux_p(m)/gpp_p(m)
            psn_to_cpool_p(m) = psn_to_cpool_p(m) * (1._r8 - downreg_p(m))

	 end if

         ! calculate the amount of new leaf C dictated by these allocation
         ! decisions, and calculate the daily fluxes of C and N to current
         ! growth and storage pools

         ! fcur is the proportion of this day's growth that is displayed now,
         ! the remainder going into storage for display next year through the
         ! transfer pools

         nlc = plant_calloc_p(m) / c_allometry_p(m)
         cpool_to_leafc_p(m)          = nlc * fcur
         cpool_to_leafc_storage_p(m)  = nlc * (1._r8 - fcur)
         cpool_to_frootc_p(m)         = nlc * f1 * fcur
         cpool_to_frootc_storage_p(m) = nlc * f1 * (1._r8 - fcur)
         if (woody(ivt) == 1._r8) then
            cpool_to_livestemc_p(m)          = nlc * f3 * f4 * fcur
            cpool_to_livestemc_storage_p(m)  = nlc * f3 * f4 * (1._r8 - fcur)
            cpool_to_deadstemc_p(m)          = nlc * f3 * (1._r8 - f4) * fcur
            cpool_to_deadstemc_storage_p(m)  = nlc * f3 * (1._r8 - f4) * (1._r8 - fcur)
            cpool_to_livecrootc_p(m)         = nlc * f2 * f3 * f4 * fcur
            cpool_to_livecrootc_storage_p(m) = nlc * f2 * f3 * f4 * (1._r8 - fcur)
            cpool_to_deadcrootc_p(m)         = nlc * f2 * f3 * (1._r8 - f4) * fcur
            cpool_to_deadcrootc_storage_p(m) = nlc * f2 * f3 * (1._r8 - f4) * (1._r8 - fcur)
         end if
#ifdef CROP
         if (ivt >= npcropmin) then ! skip 2 generic crops
            cpool_to_livestemc_p(m)          = nlc * f3 * f4 * fcur
            cpool_to_livestemc_storage_p(m)  = nlc * f3 * f4 * (1._r8 - fcur)
            cpool_to_deadstemc_p(m)          = nlc * f3 * (1._r8 - f4) * fcur
            cpool_to_deadstemc_storage_p(m)  = nlc * f3 * (1._r8 - f4) * (1._r8 - fcur)
            cpool_to_livecrootc_p(m)         = nlc * f2 * f3 * f4 * fcur
            cpool_to_livecrootc_storage_p(m) = nlc * f2 * f3 * f4 * (1._r8 - fcur)
            cpool_to_deadcrootc_p(m)         = nlc * f2 * f3 * (1._r8 - f4) * fcur
            cpool_to_deadcrootc_storage_p(m) = nlc * f2 * f3 * (1._r8 - f4) * (1._r8 - fcur)
            cpool_to_grainc_p(m)             = nlc * f5 * fcur
            cpool_to_grainc_storage_p(m)     = nlc * f5 * (1._r8 -fcur)
         end if
#endif

         ! corresponding N fluxes
         npool_to_leafn_p(m)          = (nlc / cnl) * fcur
         npool_to_leafn_storage_p(m)  = (nlc / cnl) * (1._r8 - fcur)
         npool_to_frootn_p(m)         = (nlc * f1 / cnfr) * fcur
         npool_to_frootn_storage_p(m) = (nlc * f1 / cnfr) * (1._r8 - fcur)
         if (woody(ivt) == 1._r8) then
            npool_to_livestemn_p(m)          = (nlc * f3 * f4 / cnlw) * fcur
            npool_to_livestemn_storage_p(m)  = (nlc * f3 * f4 / cnlw) * (1._r8 - fcur)
            npool_to_deadstemn_p(m)          = (nlc * f3 * (1._r8 - f4) / cndw) * fcur
            npool_to_deadstemn_storage_p(m)  = (nlc * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
            npool_to_livecrootn_p(m)         = (nlc * f2 * f3 * f4 / cnlw) * fcur
            npool_to_livecrootn_storage_p(m) = (nlc * f2 * f3 * f4 / cnlw) * (1._r8 - fcur)
            npool_to_deadcrootn_p(m)         = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * fcur
            npool_to_deadcrootn_storage_p(m) = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
         end if
#ifdef CROP
         if (ivt >= npcropmin) then ! skip 2 generic crops
            cng = graincn(ivt)
            npool_to_livestemn_p(m)          = (nlc * f3 * f4 / cnlw) * fcur
            npool_to_livestemn_storage_p(m)  = (nlc * f3 * f4 / cnlw) * (1._r8 - fcur)
            npool_to_deadstemn_p(m)          = (nlc * f3 * (1._r8 - f4) / cndw) * fcur
            npool_to_deadstemn_storage_p(m)  = (nlc * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
            npool_to_livecrootn_p(m)         = (nlc * f2 * f3 * f4 / cnlw) * fcur
            npool_to_livecrootn_storage_p(m) = (nlc * f2 * f3 * f4 / cnlw) * (1._r8 - fcur)
            npool_to_deadcrootn_p(m)         = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * fcur
            npool_to_deadcrootn_storage_p(m) = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
            npool_to_grainn_p(m)             = (nlc * f5 / cng) * fcur
            npool_to_grainn_storage_p(m)     = (nlc * f5 / cng) * (1._r8 -fcur)
         end if		
#endif

         ! Calculate the amount of carbon that needs to go into growth
         ! respiration storage to satisfy all of the storage growth demands.
         ! Allows for the fraction of growth respiration that is released at the
         ! time of fixation, versus the remaining fraction that is stored for
         ! release at the time of display. Note that all the growth respiration
         ! fluxes that get released on a given timestep are calculated in growth_resp(),
         ! but that the storage of C for growth resp during display of transferred
         ! growth is assigned here.

         gresp_storage = cpool_to_leafc_storage_p(m) + cpool_to_frootc_storage_p(m)
         if (woody(ivt) == 1._r8) then
            gresp_storage = gresp_storage + cpool_to_livestemc_storage_p(m)
            gresp_storage = gresp_storage + cpool_to_deadstemc_storage_p(m)

            gresp_storage = gresp_storage + cpool_to_livecrootc_storage_p(m)
            gresp_storage = gresp_storage + cpool_to_deadcrootc_storage_p(m)
         end if
         if (ivt >= npcropmin) then ! skip 2 generic crops
            gresp_storage = gresp_storage + cpool_to_livestemc_storage_p(m)
            gresp_storage = gresp_storage + cpool_to_grainc_storage_p(m)
         end if
         cpool_to_gresp_storage_p(m) = gresp_storage * g1 * (1._r8 - g2)

      end do ! end patch loop

 end subroutine calc_plant_nutrient_competition_CLM45_default

 subroutine calc_plant_nutrient_demand_CLM45_default(i,ps,pe,deltim,npcropmin)

    !----------------------------------------------------------------------------
    ! !DESCRIPTION
    ! Calculate allocation fraction and plant nitrogen demand.
    !
    ! !Original:
    ! The Community Land Model version 5.0 (CLM5.0)

    ! !REVISION:
    ! Xingjie Lu, 2022, modify original CLM5 to be compatible with CoLM code structure. 
    ! Fang Li, 2022, add GPAM C allocation scheme for crop.


    integer ,intent(in) :: i
    integer ,intent(in) :: ps
    integer ,intent(in) :: pe
    real(r8),intent(in) :: deltim
    integer ,intent(in) :: npcropmin

    ! !LOCAL VARIABLES:
    integer :: j            ! indices
    real(r8):: mr                 ! maintenance respiration (gC/m2/s)
    real(r8):: f1,f2,f3,f4,g1,g2  ! allocation parameters
    real(r8):: cnl,cnfr,cnlw,cndw ! C:N ratios for leaf, fine root, and wood
    real(r8):: curmr, curmr_ratio ! xsmrpool temporary variables
    real(r8):: f5                 ! grain allocation parameter
    real(r8):: cng                ! C:N ratio for grain (= cnlw for now; slevis)
    real(r8):: fleaf              ! fraction allocated to leaf
    real(r8):: t1                 ! temporary variable
    real(r8):: dayscrecover       ! number of days to recover negative cpool
    integer :: ivt, m
   dayscrecover = 30._r8
   
   do m = ps, pe
      ivt = pftclass(m)
      psn_to_cpool_p(m) = assim_p(m) * 12.011_r8

     gpp_p(m) = psn_to_cpool_p(m)

   ! get the time step total maintenance respiration
   ! These fluxes should already be in gC/m2/s

      mr = leaf_mr_p(m) + froot_mr_p(m)
      if (woody(ivt) == 1.0_r8) then
         mr = mr + livestem_mr_p(m) + livecroot_mr_p(m)
      else if (ivt >= npcropmin) then
#ifdef CROP
         if (croplive_p(m)) mr = mr + livestem_mr_p(m) + grain_mr_p(m)
#endif
      end if

   ! carbon flux available for allocation
      availc_p(m) = gpp_p(m) - mr

   ! new code added for isotope calculations, 7/1/05, PET
   ! If mr > gpp, then some mr comes from gpp, the rest comes from
   ! cpool (xsmr)
      if (mr > 0._r8 .and. availc_p(m) < 0._r8) then
         curmr = gpp_p(m)
         curmr_ratio = curmr / mr
      else
         curmr_ratio = 1._r8
      end if
      leaf_curmr_p(m)      = leaf_mr_p(m) * curmr_ratio
      leaf_xsmr_p(m)       = leaf_mr_p(m) - leaf_curmr_p(m)
      froot_curmr_p(m)     = froot_mr_p(m) * curmr_ratio
      froot_xsmr_p(m)      = froot_mr_p(m) - froot_curmr_p(m)
      livestem_curmr_p(m)  = livestem_mr_p(m) * curmr_ratio
      livestem_xsmr_p(m)   = livestem_mr_p(m) - livestem_curmr_p(m)
      livecroot_curmr_p(m) = livecroot_mr_p(m) * curmr_ratio
      livecroot_xsmr_p(m)  = livecroot_mr_p(m) - livecroot_curmr_p(m)
      grain_curmr_p(m)     = grain_mr_p(m) * curmr_ratio
      grain_xsmr_p(m)      = grain_mr_p(m) - grain_curmr_p(m)

   ! no allocation when available c is negative
      availc_p(m) = max(availc_p(m),0.0_r8)

   ! test for an xsmrpool deficit
      if (xsmrpool_p(m) < 0.0_r8) then
      ! Running a deficit in the xsmrpool, so the first priority is to let
      ! some availc from this timestep accumulate in xsmrpool.
      ! Determine rate of recovery for xsmrpool deficit

         xsmrpool_recover_p(m) = -xsmrpool_p(m)/(dayscrecover*86400._r8)
         if (xsmrpool_recover_p(m) < availc_p(m)) then
         ! available carbon reduced by amount for xsmrpool recovery
            availc_p(m) = availc_p(m) - xsmrpool_recover_p(m)
         else
         ! all of the available carbon goes to xsmrpool recovery
            xsmrpool_recover_p(m) = availc_p(m)
            availc_p(m) = 0.0_r8
         end if
         cpool_to_xsmrpool_p(m) = xsmrpool_recover_p(m)
      end if
   
      f1 = froot_leaf(ivt)
      f2 = croot_stem(ivt)

   ! modified wood allocation to be 2.2 at npp=800 gC/m2/yr, 0.2 at npp=0,
   ! constrained so that it does not go lower than 0.2 (under negative annsum_npp)
   ! This variable allocation is only for trees. Shrubs have a constant
   ! allocation as specified in the pft-physiologfy file.  The value is also used
   ! as a trigger here: -1.0 means to use the dynamic allocation (trees).

      if (stem_leaf(ivt) == -1._r8) then
         f3 = (2.7/(1.0+exp(-0.004*(annsum_npp_p(m) - 300.0)))) - 0.4
      else
         f3 = stem_leaf(ivt)
      end if

      f4   = flivewd(ivt)
      g1   = grperc(ivt)
      g2   = grpnow(ivt)
      cnl  = leafcn(ivt)
      cnfr = frootcn(ivt)
      cnlw = livewdcn(ivt)
      cndw = deadwdcn(ivt)

   ! calculate f1 to f5 for prog crops following AgroIBIS subr phenocrop

      f5 = 0._r8 ! continued intializations from above
#ifdef CROP
     if (ivt >= npcropmin) then ! skip 2 generic crops

        if (croplive_p(m)) then
         ! same phases appear in subroutine CropPhenology

         ! Phase 1 completed:
         ! ==================
         ! Next phase: leaf emergence to start of leaf decline
           
            if (hui_p(m) >= lfemerg(ivt) .and. hui_p(m) < grnfill(ivt)) then
            ! allocation rules for crops based on maturity and linear decrease
            ! of amount allocated to roots over course of the growing season

               if (peaklai_p(m) == 1) then ! lai at maximum allowed
                  arepr_p(m) = 0._r8
                  aleaf_p(m) = 1.e-5_r8
                  astem_p(m) = 0._r8
                  aroot_p(m) = 1._r8 - arepr_p(m) - aleaf_p(m) - astem_p(m)
               else
                  arepr_p(m) = 0._r8
                  aroot_p(m) = arooti(ivt) - (arooti(ivt) - arootf(ivt)) * hui_p(m)
                  fleaf = fleafi(ivt) * (exp(-bfact(ivt)) -         &
                       exp(-bfact(ivt)*hui_p(m)/grnfill(ivt))) / &
                       (exp(-bfact(ivt))-1) ! fraction alloc to leaf (from J Norman alloc curve)
                  aleaf_p(m) = max(1.e-5_r8, (1._r8 - aroot_p(m)) * fleaf)
                  astem_p(m) = 1._r8 - arepr_p(m) - aleaf_p(m) - aroot_p(m)
               end if
            ! AgroIBIS included here an immediate adjustment to aleaf & astem if the 
            ! predicted lai from the above allocation coefficients exceeded laimx.
            ! We have decided to live with lais slightly higher than laimx by
            ! enforcing the cap in the following tstep through the peaklai logic above.

               astemi_p(m) = astem_p(m) ! save for use by equations after shift to reproductive
               grain_flag_p(m) = 0._r8  ! phenology stage begins setting to 0 while in phase 2

            ! Phase 2 completed:
            ! ==================
            ! shift allocation either when enough hui are accumulated or maximum number
            ! of days has elapsed since planting

            else if (hui_p(m) >= grnfill(ivt)) then

               aroot_p(m) = arooti(ivt) - (arooti(ivt) - arootf(ivt)) * min(1._r8, hui_p(m))
                  astem_p(m) = max(astemf(ivt), astem_p(m) * max(0._r8, (1._r8-hui_p(m))/  &
                          (1._r8-grnfill(ivt)))**allconss(ivt))
                  aleaf_p(m) = 1.e-5_r8

            !Beth's retranslocation of leafn, stemn, rootn to organ
            !Filter excess plant N to retransn pool for organ N
            !Only do one time then hold grain_flag till onset next season

               if (astem_p(m) == astemf(ivt) .or. &
                    (ivt /= ntmp_soybean .and. ivt /= nirrig_tmp_soybean .and.&
                     ivt /= ntrp_soybean .and. ivt /= nirrig_trp_soybean)) then
                  if (grain_flag_p(m) == 0._r8)then
                        t1 = 1._r8 / deltim
                        leafn_to_retransn_p(m) = t1 * ((leafc_p(m) / leafcn(ivt)) - (leafc_p(m) / &
                             fleafcn(ivt)))
                        livestemn_to_retransn_p(m) = t1 * ((livestemc_p(m) / livewdcn(ivt)) - (livestemc_p(m) / &
                             fstemcn(ivt)))
                        frootn_to_retransn_p(m) = 0._r8
                        if (ffrootcn(ivt) > 0._r8) then
                           frootn_to_retransn_p(m) = t1 * ((frootc_p(m) / frootcn(ivt)) - (frootc_p(m) / &
                               ffrootcn(ivt)))
                        end if
                     grain_flag_p(m) = 1._r8
                  end if
               end if

               arepr_p(m) = 1._r8 - aroot_p(m) - astem_p(m) - aleaf_p(m)
!F. Li for vernalization effect 2 
             if(ivt == nwwheat .or. ivt == nirrig_wwheat) then
                  arepr_p(m) = arepr_p(m)*vf_p(m)
                  aroot_p(m) = 1._r8 - aleaf_p(m) - astem_p(m) - arepr_p(m)
              end if
           else                   ! pre emergence
              aleaf_p(m) = 1.e-5_r8 ! allocation coefficients should be irrelevant
              astem_p(m) = 0._r8    ! because crops have no live carbon pools;
              aroot_p(m) = 0._r8    ! this applies to this "else" and to the "else"
              arepr_p(m) = 0._r8    ! a few lines down
           end if

           f1 = aroot_p(m) / aleaf_p(m)
           f3 = astem_p(m) / aleaf_p(m)
           f5 = arepr_p(m) / aleaf_p(m)
           g1 = grperc(ivt)

        else   ! .not croplive
           f1 = 0._r8
           f3 = 0._r8
           f5 = 0._r8
           g1 = grperc(ivt)
        end if
     end if
#endif

   ! based on available C, use constant allometric relationships to
   ! determine N requirements
   
   !RF. I removed the growth respiration from this, because it is used to calculate 
   !plantCN for N uptake AND c_allometry for allocation. If we add gresp to the 
   !allometry calculation then we allocate too much carbon since gresp is not allocated here. 
     if (woody(ivt) == 1.0_r8) then
        c_allometry_p(m) = (1._r8+g1)*(1._r8+f1+f3*(1._r8+f2))
        n_allometry_p(m) = 1._r8/cnl + f1/cnfr + (f3*f4*(1._r8+f2))/cnlw + &
              (f3*(1._r8-f4)*(1._r8+f2))/cndw
     else if (ivt >= npcropmin) then ! skip generic crops
#ifdef CROP
        cng = graincn(ivt)
        c_allometry_p(m) = (1._r8+g1)*(1._r8+f1+f5+f3*(1._r8+f2))
        n_allometry_p(m) = 1._r8/cnl + f1/cnfr + f5/cng + (f3*f4*(1._r8+f2))/cnlw + &
              (f3*(1._r8-f4)*(1._r8+f2))/cndw
#endif
     else
        c_allometry_p(m) = 1._r8+g1+f1+f1*g1
        n_allometry_p(m) = 1._r8/cnl + f1/cnfr
     end if           
         
     plant_ndemand_p(m) = availc_p(m)*(n_allometry_p(m)/c_allometry_p(m))
 
   ! retranslocated N deployment depends on seasonal cycle of potential GPP
   ! (requires one year run to accumulate demand)

     tempsum_potential_gpp_p(m) = tempsum_potential_gpp_p(m) + gpp_p(m)

   ! Adding the following line to carry max retransn info to CN Annual Update
     tempmax_retransn_p(m) = max(tempmax_retransn_p(m),retransn_p(m))

   ! Beth's code: crops pull from retransn pool only during grain fill;
   !              retransn pool has N from leaves, stems, and roots for
   !              retranslocation
         
     if (ivt >= npcropmin .and. grain_flag_p(m) == 1._r8) then
        avail_retransn_p(m) = plant_ndemand_p(m)
     else if (ivt < npcropmin .and. annsum_potential_gpp_p(m) > 0._r8) then
        avail_retransn_p(m) = (annmax_retransn_p(m)/2._r8)*(gpp_p(m)/annsum_potential_gpp_p(m))/deltim
     else
        avail_retransn_p(m) = 0.0_r8
     end if
   
      ! make sure available retrans N doesn't exceed storage
     avail_retransn_p(m) = min(avail_retransn_p(m), retransn_p(m)/deltim)

      ! modify plant N demand according to the availability of
      ! retranslocated N
      ! take from retransn pool at most the flux required to meet
      ! plant ndemand

     if (plant_ndemand_p(m) > avail_retransn_p(m)) then
        retransn_to_npool_p(m) = avail_retransn_p(m)
     else
        retransn_to_npool_p(m) = plant_ndemand_p(m)
     end if

     plant_ndemand_p(m) = plant_ndemand_p(m) - retransn_to_npool_p(m)
         
  end do ! end loop pft patch.

 end subroutine calc_plant_nutrient_demand_CLM45_default

end module MOD_BGC_Veg_NutrientCompetition
#endif
