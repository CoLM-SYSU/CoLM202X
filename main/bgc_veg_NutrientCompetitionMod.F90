module bgc_veg_NutrientCompetitionMod

    use precision
    use PFT_Const, only: &
        woody, leafcn, frootcn, livewdcn, deadwdcn, graincn, &
        froot_leaf, croot_stem, stem_leaf, flivewd, grperc, grpnow, fcur2
    
    use MOD_PFTimeInvars, only: pftclass, pftfrac

    use MOD_PFTimeVars, only: &
        xsmrpool_p, retransn_p, &
        tempsum_potential_gpp_p, tempmax_retransn_p, annmax_retransn_p, annsum_potential_gpp_p, &
        c_allometry_p, n_allometry_p, downreg_p, grain_flag_p, annsum_npp_p

    use MOD_TimeVariables, only: fpg 

    use MOD_1D_PFTFluxes, only: &
        assim_p, leaf_xsmr_p, froot_xsmr_p, livestem_xsmr_p, livecroot_xsmr_p, grain_xsmr_p, cpool_to_xsmrpool_p, &
        leaf_mr_p, froot_mr_p, livestem_mr_p, livecroot_mr_p, grain_mr_p, &
        plant_ndemand_p, retransn_to_npool_p, cpool_to_leafc_p, cpool_to_leafc_storage_p, &
        cpool_to_frootc_p, cpool_to_frootc_storage_p, cpool_to_livestemc_p, cpool_to_livestemc_storage_p, &
        cpool_to_deadstemc_p, cpool_to_deadstemc_storage_p, cpool_to_livecrootc_p, cpool_to_livecrootc_storage_p, &
        cpool_to_deadcrootc_p, cpool_to_deadcrootc_storage_p, cpool_to_grainc_p, cpool_to_grainc_storage_p, &
        cpool_to_gresp_storage_p, npool_to_leafn_p, npool_to_leafn_storage_p, &
        npool_to_frootn_p, npool_to_frootn_storage_p, npool_to_livestemn_p, npool_to_livestemn_storage_p, &
        npool_to_deadstemn_p, npool_to_deadstemn_storage_p, npool_to_livecrootn_p, npool_to_livecrootn_storage_p, &
        npool_to_deadcrootn_p, npool_to_deadcrootn_storage_p, npool_to_grainn_p, npool_to_grainn_storage_p, &
        plant_calloc_p, plant_nalloc_p, leaf_curmr_p, froot_curmr_p, livestem_curmr_p, livecroot_curmr_p, grain_curmr_p, &
        psn_to_cpool_p, gpp_p, availc_p, avail_retransn_p, xsmrpool_recover_p, sminn_to_npool_p, excess_cflux_p

implicit none

public calc_plant_nutrient_competition_CLM45_default
public calc_plant_nutrient_demand_CLM45_default

contains

 subroutine calc_plant_nutrient_competition_CLM45_default(i,ps,pe,npcropmin)

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
!    real(r8):: fsmn(bounds%begp:bounds%endp)  ! A emperate variable for adjusting FUN uptakes 
    integer :: ivt, m

!if (i  .eq. 123226)print*,'npool_to_wood before alloc',&
!   sum((npool_to_livestemn_p(ps:pe) + npool_to_livestemn_storage_p(ps:pe) &
!   +npool_to_deadstemn_p(ps:pe) + npool_to_deadstemn_storage_p(ps:pe) &
!   +npool_to_livecrootn_p(ps:pe) + npool_to_livecrootn_storage_p(ps:pe) &
!   +npool_to_deadcrootn_p(ps:pe) + npool_to_deadcrootn_storage_p(ps:pe))*pftfrac(ps:pe)*1800._r8)
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

         if (ivt >= npcropmin) then ! skip 2 generic crops
!            if (croplive(i)) then
!               f1 = aroot(i) / aleaf(i)
!               f3 = astem(i) / aleaf(i)
!               f5 = arepr(i) / aleaf(i)
!               g1 = 0.25_r8
!            else
               f1 = 0._r8
               f3 = 0._r8
               f5 = 0._r8
               g1 = 0.25_r8
!            end if
         end if

!         if(use_fun)then ! if we are using FUN, we get the N available from there.
!            sminn_to_npool(p) = sminn_to_plant_fun(p)
!         else ! no FUN. :( we get N available from the FPG calculation in soilbiogeochemistry competition. 
            sminn_to_npool_p(m) = plant_ndemand_p(m) * fpg(i)        
!         end if
!         print*,'nalloc in nutrient compeition',i,sminn_to_npool(i),plant_ndemand(i),fpg(i)

         plant_nalloc_p(m) = sminn_to_npool_p(m) + retransn_to_npool_p(m)
!if(i .eq. 79738)print*,'nalloc',i,m,plant_nalloc_p(m),sminn_to_npool_p(m),retransn_to_npool_p(m),plant_ndemand_p(m),fpg(i)
         plant_calloc_p(m) = plant_nalloc_p(m) * (c_allometry_p(m)/n_allometry_p(m))
         !if(i .eq. 79738)print*,'calloc',i,m,plant_calloc_p(m),plant_nalloc_p(m),c_allometry_p(m)/n_allometry_p(m)
         !print*,'plant_calloc in nutrient competition',i,plant_calloc(i),plant_nalloc(i),&
!         sminn_to_npool(i),retransn_to_npool(i)
!         if (use_matrixcn)then 
!            matrix_Ninput(p) =  sminn_to_npool(p)! + retransn_to_npool(p) 
!            matrix_Cinput(p) = plant_nalloc(p) * (c_allometry(p)/n_allometry(p))
!         end if


!         if(i .eq.  79738)print*,'availc',i,m,availc_p(m)
!         if(.not.use_fun)then  !ORIGINAL CLM(CN) downregulation code. 
            excess_cflux_p(m) = availc_p(m) - plant_calloc_p(m)
	    ! reduce gpp fluxes due to N limitation
	    if (gpp_p(m) > 0.0_r8) then
	       downreg_p(m) = excess_cflux_p(m)/gpp_p(m)
               !print*,'1-downreg',1._r8-downreg(i)
!	       psnsun_to_cpool(m) = psnsun_to_cpool(m)*(1._r8 - downreg(m))
!	       psnsha_to_cpool(m) = psnsha_to_cpool(m)*(1._r8 - downreg(m))
               psn_to_cpool_p(m) = psn_to_cpool_p(m) * (1._r8 - downreg_p(m))

!	       if ( use_c13 ) then
!	          c13_cnveg_carbonflux_inst%psnsun_to_cpool_patch(p)   = &
!	               c13_cnveg_carbonflux_inst%psnsun_to_cpool_patch(p)  *(1._r8 - downreg(p))
!	          c13_cnveg_carbonflux_inst%psnshade_to_cpool_patch(p) = &
!	               c13_cnveg_carbonflux_inst%psnshade_to_cpool_patch(p)*(1._r8 - downreg(p))
!	       end if
!	       if ( use_c14 ) then
!	          c14_cnveg_carbonflux_inst%psnsun_to_cpool_patch(p)   = &
!	               c14_cnveg_carbonflux_inst%psnsun_to_cpool_patch(p)  *(1._r8 - downreg(p))
!	          c14_cnveg_carbonflux_inst%psnshade_to_cpool_patch(p) = &
!	               c14_cnveg_carbonflux_inst%psnshade_to_cpool_patch(p)*(1._r8 - downreg(p))
!	       end if
               if(m .eq. 22784)print*,'downreg_p',m,downreg_p(m),gpp_p(m),excess_cflux_p(m),availc_p(m),plant_calloc_p(m)
	    end if
!	         
!	 end if !use_fun

!         if(i .eq. 79738)print*,'gpp',i,m,gpp_p(m), downreg_p(m),psn_to_cpool_p(m)

         ! calculate the amount of new leaf C dictated by these allocation
         ! decisions, and calculate the daily fluxes of C and N to current
         ! growth and storage pools

         ! fcur is the proportion of this day's growth that is displayed now,
         ! the remainder going into storage for display next year through the
         ! transfer pools

         nlc = plant_calloc_p(m) / c_allometry_p(m)
         !print*,'nlc',nlc
!         print*,'cpool_to_leafc in nutrient competition',i,cpool_to_leafc(i),nlc,fcur,plant_calloc(i),c_allometry(i),n_allometry(i)
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

         ! corresponding N fluxes
         npool_to_leafn_p(m)          = (nlc / cnl) * fcur
         npool_to_leafn_storage_p(m)  = (nlc / cnl) * (1._r8 - fcur)
         npool_to_frootn_p(m)         = (nlc * f1 / cnfr) * fcur
         npool_to_frootn_storage_p(m) = (nlc * f1 / cnfr) * (1._r8 - fcur)
!         if(i .eq. 123226)print*,'woody(ivt)',ivt, woody(ivt)
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
!         if (use_matrixcn) then
!            matrix_alloc(p,ileaf)             = (1.0_r8) / c_allometry(p) * fcur
!            matrix_alloc(p,ileaf_st)          = (1.0_r8) / c_allometry(p) * (1._r8 - fcur)
!            matrix_alloc(p,ifroot)            = (1.0_r8) / c_allometry(p) * f1 * fcur
!            matrix_alloc(p,ifroot_st)         = (1.0_r8) / c_allometry(p) * f1 * (1._r8 - fcur)

!            matrix_nalloc(p,ileaf)            = ((1.0_r8/cnl)          / n_allometry(p)) * fcur
!            matrix_nalloc(p,ileaf_st)         = ((1.0_r8/cnl)          / n_allometry(p))* (1._r8 - fcur)
!            matrix_nalloc(p,ifroot)           = ((f1/cnfr)             / n_allometry(p)) * fcur
!            matrix_nalloc(p,ifroot_st)        = ((f1/cnfr)             / n_allometry(p)) * (1._r8 - fcur)
!            if (woody(ivt(p)) == 1._r8) then
!               matrix_alloc(p,ilivestem)      = (1.0_r8) / c_allometry(p) * f3 * f4 * fcur 
!               matrix_alloc(p,ilivestem_st)   = (1.0_r8) / c_allometry(p) * f3 * f4 * (1._r8 - fcur)
!               matrix_alloc(p,ideadstem)      = (1.0_r8) / c_allometry(p) * f3 * (1._r8 - f4) * fcur
!               matrix_alloc(p,ideadstem_st)   = (1.0_r8) / c_allometry(p) * f3 * (1._r8 - f4) * (1._r8 - fcur)
!               matrix_alloc(p,ilivecroot)     = (1.0_r8) / c_allometry(p) * f2 * f3 * f4 * fcur
!               matrix_alloc(p,ilivecroot_st)  = (1.0_r8) / c_allometry(p) * f2 * f3 * f4 * (1._r8 - fcur)
!               matrix_alloc(p,ideadcroot)     = (1.0_r8) / c_allometry(p) * f2 * f3 * (1._r8 - f4) * fcur
!               matrix_alloc(p,ideadcroot_st)  = (1.0_r8) / c_allometry(p) * f2 * f3 * (1._r8 - f4) * (1._r8 - fcur)
 
!               matrix_nalloc(p,ilivestem)     = (f3*f4/cnlw)                     / n_allometry(p) * fcur 
!               matrix_nalloc(p,ilivestem_st)  = (f3*f4/cnlw)                     / n_allometry(p) * (1._r8 - fcur)
!               matrix_nalloc(p,ideadstem)     = (f3 * (1._r8 - f4)/cndw)         / n_allometry(p) * fcur
!               matrix_nalloc(p,ideadstem_st)  = (f3 * (1._r8 - f4)/cndw)         / n_allometry(p) * (1._r8 - fcur)
!               matrix_nalloc(p,ilivecroot)    = (f2 * f3 * f4/cnlw)              / n_allometry(p) *  fcur
!               matrix_nalloc(p,ilivecroot_st) = (f2 * f3 * f4 /cnlw)             / n_allometry(p) * (1._r8 - fcur)
!               matrix_nalloc(p,ideadcroot)    = (f2 * f3 * (1._r8 - f4)/cndw)    / n_allometry(p) * fcur
!               matrix_nalloc(p,ideadcroot_st) = (f2 * f3 * (1._r8 - f4)/cndw)    / n_allometry(p) *(1._r8 - fcur)
!            end if
!            if (ivt(p) >= npcropmin) then ! skip 2 generic crops
!               matrix_alloc(p,ilivestem)      = (1.0_r8) / c_allometry(p) * f3 * f4 * fcur 
!               matrix_alloc(p,ilivestem_st)   = (1.0_r8) / c_allometry(p) * f3 * f4 * (1._r8 - fcur)
!               matrix_alloc(p,ideadstem)      = (1.0_r8) / c_allometry(p) * f3 * (1._r8 - f4) * fcur
!               matrix_alloc(p,ideadstem_st)   = (1.0_r8) / c_allometry(p) * f3 * (1._r8 - f4) * (1._r8 - fcur)
!               matrix_alloc(p,ilivecroot)     = (1.0_r8) / c_allometry(p) * f2 * f3 * f4 * fcur
!               matrix_alloc(p,ilivecroot_st)  = (1.0_r8) / c_allometry(p) * f2 * f3 * f4 * (1._r8 - fcur)
!               matrix_alloc(p,ideadcroot)     = (1.0_r8) / c_allometry(p) * f2 * f3 * (1._r8 - f4) * fcur
!               matrix_alloc(p,ideadcroot_st)  = (1.0_r8) / c_allometry(p) * f2 * f3 * (1._r8 - f4) * (1._r8 - fcur)
!               matrix_alloc(p,igrain)         = (1.0_r8) / c_allometry(p) * f5 * fcur
!               matrix_alloc(p,igrain_st)      = (1.0_r8) / c_allometry(p) * f5 * (1._r8 - fcur)
! 
!               matrix_nalloc(p,ilivestem)     = (f3*f4/cnlw)                     / n_allometry(p) * fcur 
!               matrix_nalloc(p,ilivestem_st)  = (f3*f4/cnlw)                     / n_allometry(p) * (1._r8 - fcur)
!               matrix_nalloc(p,ideadstem)     = (f3 * (1._r8 - f4)/cndw)         / n_allometry(p) * fcur
!               matrix_nalloc(p,ideadstem_st)  = (f3 * (1._r8 - f4)/cndw)         / n_allometry(p) * (1._r8 - fcur)
!               matrix_nalloc(p,ilivecroot)    = (f2 * f3 * f4/cnlw)              / n_allometry(p) *  fcur
!               matrix_nalloc(p,ilivecroot_st) = (f2 * f3 * f4 /cnlw)             / n_allometry(p) * (1._r8 - fcur)
!               matrix_nalloc(p,ideadcroot)    = (f2 * f3 * (1._r8 - f4)/cndw)    / n_allometry(p) * fcur
!               matrix_nalloc(p,ideadcroot_st) = (f2 * f3 * (1._r8 - f4)/cndw)    / n_allometry(p) * (1._r8 - fcur)
!               matrix_nalloc(p,igrain)        = (f5 / cng)   / n_allometry(p)    * fcur
!               matrix_nalloc(p,igrain_st)     = (f5 / cng)   / n_allometry(p)    *(1._r8 - fcur)
!            end if
!         end if !end use_matrixcn

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

!         if(use_matrixcn)then
!            matrix_Cinput(p) = plant_calloc(p)
!            if(retransn(p) .ne. 0)then
!               matrix_nphtransfer(p,iretransn_to_ileaf)           = matrix_nphtransfer(p,iretransn_to_ileaf) &
!                                                                  + matrix_nalloc(p,ileaf    )     * retransn_to_npool(p) / retransn(p)
!               matrix_nphtransfer(p,iretransn_to_ileafst)         = matrix_nphtransfer(p,iretransn_to_ileafst) &
!                                                                  + matrix_nalloc(p,ileaf_st )     * retransn_to_npool(p) / retransn(p)
!               matrix_nphtransfer(p,iretransn_to_ifroot)          = matrix_nphtransfer(p,iretransn_to_ifroot) &
!                                                                  + matrix_nalloc(p,ifroot   )     * retransn_to_npool(p) / retransn(p)
!               matrix_nphtransfer(p,iretransn_to_ifrootst)        = matrix_nphtransfer(p,iretransn_to_ifrootst) &
!                                                                  + matrix_nalloc(p,ifroot_st)     * retransn_to_npool(p) / retransn(p)
!               matrix_nphtransfer(p,iretransn_to_ilivestem)       = matrix_nphtransfer(p,iretransn_to_ilivestem) &
!                                                                  + matrix_nalloc(p,ilivestem    ) * retransn_to_npool(p) / retransn(p)
!               matrix_nphtransfer(p,iretransn_to_ilivestemst)     = matrix_nphtransfer(p,iretransn_to_ilivestemst) &
!                                                                  + matrix_nalloc(p,ilivestem_st ) * retransn_to_npool(p) / retransn(p)
!               matrix_nphtransfer(p,iretransn_to_ideadstem)       = matrix_nphtransfer(p,iretransn_to_ideadstem) &
!                                                                  + matrix_nalloc(p,ideadstem    ) * retransn_to_npool(p) / retransn(p)
!               matrix_nphtransfer(p,iretransn_to_ideadstemst)     = matrix_nphtransfer(p,iretransn_to_ideadstemst) &
!                                                                  + matrix_nalloc(p,ideadstem_st ) * retransn_to_npool(p) / retransn(p)
!               matrix_nphtransfer(p,iretransn_to_ilivecroot)      = matrix_nphtransfer(p,iretransn_to_ilivecroot) &
!                                                                  + matrix_nalloc(p,ilivecroot   ) * retransn_to_npool(p) / retransn(p)
!               matrix_nphtransfer(p,iretransn_to_ilivecrootst)    = matrix_nphtransfer(p,iretransn_to_ilivecrootst) &
!                                                                  + matrix_nalloc(p,ilivecroot_st) * retransn_to_npool(p) / retransn(p)
!               matrix_nphtransfer(p,iretransn_to_ideadcroot)      = matrix_nphtransfer(p,iretransn_to_ideadcrootst) &
!                                                                  + matrix_nalloc(p,ideadcroot   ) * retransn_to_npool(p) / retransn(p)
!               matrix_nphtransfer(p,iretransn_to_ideadcrootst)    = matrix_nphtransfer(p,iretransn_to_ideadcrootst) &
!                                                                  + matrix_nalloc(p,ideadcroot_st) * retransn_to_npool(p) / retransn(p)
!               if(ivt(p) >= npcropmin)then 
!                  matrix_nphtransfer(p,iretransn_to_igrain)       = matrix_nphtransfer(p,iretransn_to_igrain) &
!                                                                  + matrix_nalloc(p,igrain       ) * retransn_to_npool(p) / retransn(p)
!                  matrix_nphtransfer(p,iretransn_to_igrainst)     = matrix_nphtransfer(p,iretransn_to_igrainst) &
!                                                                  + matrix_nalloc(p,igrain_st    ) * retransn_to_npool(p) / retransn(p)
!               end if
!            end if
!         end if !end use_matrixcn  
      end do ! end patch loop
!if (i  .eq. 123226)print*,'npool_to_wood',&
!   sum((npool_to_livestemn_p(ps:pe) + npool_to_livestemn_storage_p(ps:pe) &
!   +npool_to_deadstemn_p(ps:pe) + npool_to_deadstemn_storage_p(ps:pe) &
!   +npool_to_livecrootn_p(ps:pe) + npool_to_livecrootn_storage_p(ps:pe) &
!   +npool_to_deadcrootn_p(ps:pe) + npool_to_deadcrootn_storage_p(ps:pe))*pftfrac(ps:pe)*1800._r8)
!
!    end associate 

 end subroutine calc_plant_nutrient_competition_CLM45_default

 subroutine calc_plant_nutrient_demand_CLM45_default(i,ps,pe,deltim,npcropmin)

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
!      psnsun_to_cpool(m) = assimsun_p(m) * laisun_p(m) * 12.011_r8
!      psnsha_to_cpool(m) = assimsha_p(m) * laisha_p(m) * 12.011_r8
      psn_to_cpool_p(m) = assim_p(m) * 12.011_r8
      if(assim_p(m) .lt.  0)print*,'m',assim_p(m)

  !       if ( use_c13 ) then
  !          c13_cnveg_carbonflux_inst%psnsun_to_cpool_patch(p)   = c13_psnsun(p) * laisun(p) * 12.011e-6_r8
  !          c13_cnveg_carbonflux_inst%psnshade_to_cpool_patch(p) = c13_psnsha(p) * laisha(p) * 12.011e-6_r8
  !       endif

  !       if ( use_c14 ) then
  !          c14_cnveg_carbonflux_inst%psnsun_to_cpool_patch(p)   = c14_psnsun(p) * laisun(p) * 12.011e-6_r8
  !          c14_cnveg_carbonflux_inst%psnshade_to_cpool_patch(p) = c14_psnsha(p) * laisha(p) * 12.011e-6_r8
  !       endif

!      gpp(m) = psnsun_to_cpool(m) + psnsha_to_cpool(m)
     gpp_p(m) = psn_to_cpool_p(m)

   ! get the time step total maintenance respiration
   ! These fluxes should already be in gC/m2/s

      mr = leaf_mr_p(m) + froot_mr_p(m)
      if (woody(ivt) == 1.0_r8) then
         mr = mr + livestem_mr_p(m) + livecroot_mr_p(m)
      else if (ivt >= npcropmin) then
!         if (croplive(i)) mr = mr + livestem_mr(i) + grain_mr(i)
      end if

   ! carbon flux available for allocation
      availc_p(m) = gpp_p(m) - mr
   !print*,'gpp',gpp(i),mr,availc(i),leaf_mr(i),froot_mr(i),livestem_mr(i),livecroot_mr(i),laisun(i),laisha(i)

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

!   print*,'availc',i,availc(i),xsmrpool(i),gpp(i),mr
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
!   print*,'avail after xsmr recover',availc(i),xsmrpool(i)

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

      if (ivt >= npcropmin) then ! skip 2 generic crops

!      if (croplive(i)) then
         ! same phases appear in subroutine CropPhenology

         ! Phase 1 completed:
         ! ==================
         ! if hui is less than the number of gdd needed for filling of grain
         ! leaf emergence also has to have taken place for lai changes to occur
         ! and carbon assimilation
         ! Next phase: leaf emergence to start of leaf decline

!         if (leafout(i) >= huileaf(i) .and. hui(i) < huigrain(i)) then

            ! allocation rules for crops based on maturity and linear decrease
            ! of amount allocated to roots over course of the growing season

!            if (peaklai(i) == 1) then ! lai at maximum allowed
!               arepr(i) = 0._r8
!               aleaf(i) = 1.e-5_r8
!               astem(i) = 0._r8
!               aroot(i) = 1._r8 - arepr(i) - aleaf(i) - astem(i)
!            else
!               arepr(i) = 0._r8
!               aroot(i) = max(0._r8, min(1._r8, arooti(ivt) -   &
!                    (arooti(ivt) - arootf(ivt)) *  &
!                    min(1._r8, hui(i)/gddmaturity(i))))
!               fleaf = fleafi(ivt) * (exp(-bfact(ivt)) -         &
!                    exp(-bfact(ivt)*hui(i)/huigrain(i))) / &
!                    (exp(-bfact(ivt))-1) ! fraction alloc to leaf (from J Norman alloc curve)
!               aleaf(i) = max(1.e-5_r8, (1._r8 - aroot(i)) * fleaf)
!               astem(i) = 1._r8 - arepr(i) - aleaf(i) - aroot(i)
!            end if
!
            ! AgroIBIS included here an immediate adjustment to aleaf & astem if the 
            ! predicted lai from the above allocation coefficients exceeded laimx.
            ! We have decided to live with lais slightly higher than laimx by
            ! enforcing the cap in the following tstep through the peaklai logic above.

!            astemi(i) = astem(i) ! save for use by equations after shift
!            aleafi(i) = aleaf(i) ! to reproductive phenology stage begins
!            grain_flag(i) = 0._r8 ! setting to 0 while in phase 2

            ! Phase 2 completed:
            ! ==================
            ! shift allocation either when enough gdd are accumulated or maximum number
            ! of days has elapsed since planting

!         else if (hui(i) >= huigrain(i)) then

!            aroot(i) = max(0._r8, min(1._r8, arooti(ivt) - &
!                 (arooti(ivt) - arootf(ivt)) * min(1._r8, hui(i)/gddmaturity(i))))
!            if (astemi(i) > astemf(ivt)) then
!               astem(i) = max(0._r8, max(astemf(ivt), astem(i) * &
!                    (1._r8 - min((hui(i)-                 &
!                    huigrain(i))/((gddmaturity(i)*declfact(ivt))- &
!                    huigrain(i)),1._r8)**allconss(ivt) )))
!            end if
!            if (aleafi(i) > aleaff(ivt)) then
!               aleaf(i) = max(1.e-5_r8, max(aleaff(ivt), aleaf(i) * &
!                    (1._r8 - min((hui(i)-                    &
!                    huigrain(i))/((gddmaturity(i)*declfact(ivt))- &
!                    huigrain(i)),1._r8)**allconsl(ivt) )))
!            end if

            !Beth's retranslocation of leafn, stemn, rootn to organ
            !Filter excess plant N to retransn pool for organ N
            !Only do one time then hold grain_flag till onset next season

            ! slevis: Will astem ever = astemf exactly?
            ! Beth's response: ...looks like astem can equal astemf under the right circumstances. 
            !It might be worth a rewrite to capture what I was trying to do, but the retranslocation for 
            !corn and wheat begins at the beginning of the grain fill stage, but for soybean I was holding it 
            !until after the leaf and stem decline were complete. Looking at how astem is calculated, once the 
            !stem decline is near complete, astem should (usually) be set to astemf. The reason for holding off 
            !on soybean is that the retranslocation scheme begins at the beginning of the grain phase, when the 
            !leaf and stem are still growing, but declining. Since carbon is still getting allocated and now 
            !there is more nitrogen available, the nitrogen can be diverted from grain. For corn and wheat 
            !the impact was probably enough to boost productivity, but for soybean the nitrogen was better off 
            !fulfilling the grain fill. It seems that if the peak lai is reached for soybean though that this 
            !would be bypassed altogether, not the intended outcome. I checked several of my output files and 
            !they all seemed to be going through the retranslocation loop for soybean - good news.

!            if (astem(i) == astemf(ivt) .or. &
!                 (ivt /= ntmp_soybean .and. ivt /= nirrig_tmp_soybean .and.&
!                  ivt /= ntrp_soybean .and. ivt /= nirrig_trp_soybean)) then
!               if (grain_flag(i) == 0._r8)then
!                  if(.not.use_fun) then
!                     t1 = 1 / dt
!                     leafn_to_retransn(i) = t1 * ((leafc(i) / leafcn(ivt)) - (leafc(i) / &
!                          fleafcn(ivt)))
!                     livestemn_to_retransn(i) = t1 * ((livestemc(i) / livewdcn(ivt)) - (livestemc(i) / &
!                          fstemcn(ivt)))
!                     frootn_to_retransn(i) = 0._r8
!                     if (ffrootcn(ivt) > 0._r8) then
!                        frootn_to_retransn(i) = t1 * ((frootc(i) / frootcn(ivt)) - (frootc(i) / &
!                            ffrootcn(ivt)))
!                     end if
!                  else !leafn retrans flux is handled in phenology
!                     frootn_to_retransn(i) = 0._r8
!                     livestemn_to_retransn(i)=0.0_r8 
!                  end if !fun
!                  grain_flag(i) = 1._r8
!                  if(use_matrixcn)then
!                     if(leafn(i) .ne. 0._r8)then
!                        matrix_nphtransfer(p,ileaf_to_iretransn)  = (leafn_to_retransn(i)) / leafn(i)
!                     end if
!                     if(frootn(i) .ne. 0._r8)then
!                        matrix_nphtransfer(p,ifroot_to_iretransn) = (frootn_to_retransn(i))/ frootn(i)
!                     end if
!                     if(livestemn(i) .ne. 0._r8)then
!                        matrix_nphtransfer(p,ilivestem_to_iretransn) = (livestemn_to_retransn(i))/ livestemn(i)
!                     end if
!                  end if
!               end if
!            end if

!            arepr(i) = 1._r8 - aroot(i) - astem(i) - aleaf(i)

!        else                   ! pre emergence
!           aleaf(i) = 1.e-5_r8 ! allocation coefficients should be irrelevant
!           astem(i) = 0._r8    ! because crops have no live carbon pools;
!           aroot(i) = 0._r8    ! this applies to this "else" and to the "else"
!           arepr(i) = 0._r8    ! a few lines down
!        end if

!        f1 = aroot(i) / aleaf(i)
!        f3 = astem(i) / aleaf(i)
!        f5 = arepr(i) / aleaf(i)
!        g1 = 0.25_r8

!     else   ! .not croplive
            f1 = 0._r8
            f3 = 0._r8
            f5 = 0._r8
            g1 = 0.25_r8
!     end if
      end if

   ! based on available C, use constant allometric relationships to
   ! determine N requirements
   
   !RF. I removed the growth respiration from this, because it is used to calculate 
   !plantCN for N uptake AND c_allometry for allocation. If we add gresp to the 
   !allometry calculation then we allocate too much carbon since gresp is not allocated here. 
!   if(.not.use_fun)then
      if (woody(ivt) == 1.0_r8) then
         c_allometry_p(m) = (1._r8+g1)*(1._r8+f1+f3*(1._r8+f2))
         n_allometry_p(m) = 1._r8/cnl + f1/cnfr + (f3*f4*(1._r8+f2))/cnlw + &
               (f3*(1._r8-f4)*(1._r8+f2))/cndw
      else if (ivt >= npcropmin) then ! skip generic crops
         cng = graincn(ivt)
         c_allometry_p(m) = (1._r8+g1)*(1._r8+f1+f5+f3*(1._r8+f2))
         n_allometry_p(m) = 1._r8/cnl + f1/cnfr + f5/cng + (f3*f4*(1._r8+f2))/cnlw + &
              (f3*(1._r8-f4)*(1._r8+f2))/cndw
      else
         c_allometry_p(m) = 1._r8+g1+f1+f1*g1
         n_allometry_p(m) = 1._r8/cnl + f1/cnfr
      end if           
!   else !no FUN. 
!      if (woody(ivt(p)) == 1.0_r8) then
!         c_allometry(p) = (1._r8)*(1._r8+f1+f3*(1._r8+f2))
!         n_allometry(p) = 1._r8/cnl + f1/cnfr + (f3*f4*(1._r8+f2))/cnlw + &
!              (f3*(1._r8-f4)*(1._r8+f2))/cndw
!      else if (ivt(p) >= npcropmin) then ! skip generic crops
!         cng = graincn(ivt(p))
!         c_allometry(p) = (1._r8)*(1._r8+f1+f5+f3*(1._r8+f2))
!         n_allometry(p) = 1._r8/cnl + f1/cnfr + f5/cng + (f3*f4*(1._r8+f2))/cnlw + &
!              (f3*(1._r8-f4)*(1._r8+f2))/cndw
!      else
!         c_allometry(p) = 1._r8+f1
!         n_allometry(p) = 1._r8/cnl + f1/cnfr
!      end if
!   end if !use_fun
         
      plant_ndemand_p(m) = availc_p(m)*(n_allometry_p(m)/c_allometry_p(m))
!   print*,'plant_ndemand',i,plant_ndemand(i), availc(i),(n_allometry(i)/c_allometry(i))
 
   ! retranslocated N deployment depends on seasonal cycle of potential GPP
   ! (requires one year run to accumulate demand)

      tempsum_potential_gpp_p(m) = tempsum_potential_gpp_p(m) + gpp_p(m)

   ! Adding the following line to carry max retransn info to CN Annual Update
      tempmax_retransn_p(m) = max(tempmax_retransn_p(m),retransn_p(m))

   ! Beth's code: crops pull from retransn pool only during grain fill;
   !              retransn pool has N from leaves, stems, and roots for
   !              retranslocation
         
!   if(.not.use_fun)then

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

!      if ( .not. use_fun ) then
         plant_ndemand_p(m) = plant_ndemand_p(m) - retransn_to_npool_p(m)
!      else
!         if (season_decid(ivt) == 1._r8.or.stress_decid(ivt)==1._r8) then
!            plant_ndemand(m) = plant_ndemand(m) - retransn_to_npool(m)
!         end if
!      end if
         
!   end if !use_fun
    end do ! end loop pft patch.

 end subroutine calc_plant_nutrient_demand_CLM45_default

end module bgc_veg_NutrientCompetitionMod
