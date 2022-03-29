#include <define.h>
#ifdef BGC
module bgc_veg_CNVegStructUpdateMod

use precision
use GlobalVars, only: nc3crop, nc3irrig, nbrdlf_evr_shrub, nbrdlf_dcd_brl_shrub, &
                              npcropmin, ntmp_corn, nirrig_tmp_corn, ntrp_corn, nirrig_trp_corn, &
                              nsugarcane, nirrig_sugarcane, nmiscanthus, nirrig_miscanthus, &
                              nswitchgrass, nirrig_switchgrass, noveg

use MOD_PFTimeVars, only: tlai_p, tsai_p, leafc_p, deadstemc_p, peaklai_p, harvdate_p
use MOD_PFTimeInvars, only: pftclass
use MOD_TimeVariables, only: farea_burned
use PFT_Const, only : dsladlai, slatop, laimx, woody
  !-----------------------------------------------------------------------
  ! Module for vegetation structure updates (LAI, SAI, htop, hbot)
  !
  ! !USES:
  ! !PUBLIC MEMBER FUNCTIONS:
public :: CNVegStructUpdate
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CNVegStructUpdate(i,ps,pe,deltim,npcropmin)
    ! !DESCRIPTION:
    ! On the radiation time step, use C state variables and epc to diagnose
    ! vegetation structure (LAI, SAI, height)
    !
    ! !USES:

    integer,intent(in)  :: i
    integer,intent(in)  :: ps
    integer,intent(in)  :: pe
    real(r8),intent(in) :: deltim
    integer,intent(in)  :: npcropmin
    
    !
    ! !REVISION HISTORY:
    ! 10/28/03: Created by Peter Thornton
    ! 2/29/08, David Lawrence: revised snow burial fraction for short vegetation
    !
    ! !LOCAL VARIABLES:
    integer  :: p,c,g      ! indices
    integer  :: fp         ! lake filter indices
    real(r8) :: stocking                            ! #stems / ha (stocking density)
    real(r8) :: ol         ! thickness of canopy layer covered by snow (m)
    real(r8) :: fb         ! fraction of canopy layer covered by snow
    real(r8) :: tlai_old   ! for use in Zeng tsai formula
    real(r8) :: tsai_old   ! for use in Zeng tsai formula
    real(r8) :: tsai_min   ! PATCH derived minimum tsai
    real(r8) :: tsai_alpha ! monthly decay rate of tsai
    real(r8) :: frac_sno_adjusted ! frac_sno adjusted per frac_sno_threshold

    real(r8), parameter :: dtsmonth = 2592000._r8 ! number of seconds in a 30 day month (60x60x24x30)
    real(r8), parameter :: frac_sno_threshold = 0.999_r8  ! frac_sno values greater than this are treated as 1
    integer m, ivt
    !-----------------------------------------------------------------------
    ! tsai formula from Zeng et. al. 2002, Journal of Climate, p1835
    !
    ! tsai(m) = max( tsai_alpha(ivt(m))*tsai_old + max(tlai_old-tlai(m),0_r8), tsai_min(ivt(m)) )
    ! notes:
    ! * RHS tsai & tlai are from previous timestep
    ! * should create tsai_alpha(ivt(m)) & tsai_min(ivt(m)) in pftconMod.F90 - slevis
    ! * all non-crop patches use same values:
    !   crop    tsai_alpha,tsai_min = 0.0,0.1
    !   noncrop tsai_alpha,tsai_min = 0.5,1.0  (includes bare soil and urban)
    !-------------------------------------------------------------------------------
    
      ! patch loop

     do m = ps, pe
        ivt = pftclass(m)
        if (ivt /= noveg) then

            tlai_old = tlai_p(m) ! n-1 value
            tsai_old = tsai_p(m) ! n-1 value

#ifdef LAIfdbk
            ! update the leaf area index based on leafC and SLA
            ! Eq 3 from Thornton and Zimmerman, 2007, J Clim, 20, 3902-3923. 
!            if (dsladlai(ivt) > 0._r8) then
!               tlai_p(m) = (slatop(ivt)*(exp(leafc_p(m)*dsladlai(ivt)) - 1._r8))/dsladlai(ivt)
!            else
               tlai_p(m) = slatop(ivt) * leafc_p(m)
!            end if
            tlai_p(m) = max(0._r8, tlai_p(m))
#endif

            ! update the stem area index and height based on LAI, stem mass, and veg type.
            ! With the exception of htop for woody vegetation, this follows the DGVM logic.

            ! tsai formula from Zeng et. al. 2002, Journal of Climate, p1835 (see notes)
            ! Assumes doalb time step .eq. CLM time step, SAI min and monthly decay factor
            ! alpha are set by PFT, and alpha is scaled to CLM time step by multiplying by
            ! deltim and dividing by dtsmonth (seconds in average 30 day month)
            ! tsai_min scaled by 0.5 to match MODIS satellite derived values
            if (ivt == nc3crop .or. ivt == nc3irrig) then ! generic crops

               tsai_alpha = 1.0_r8-1.0_r8*deltim/dtsmonth
               tsai_min = 0.1_r8
            else
               tsai_alpha = 1.0_r8-0.5_r8*deltim/dtsmonth
               tsai_min = 1.0_r8
            end if
            tsai_min = tsai_min * 0.5_r8
            tsai_p(m) = max(tsai_alpha*tsai_old+max(tlai_old-tlai_p(m),0._r8),tsai_min)

            ! calculate vegetation physiological parameters used in biomass heat storage
            !
!            if (use_biomass_heat_storage) then
               ! Assumes fbw (fraction of biomass that is water) is the same for leaves and stems
!               leaf_biomass(m) = max(0.0025_r8,leafc(m)) &
!                    * c_to_b * 1.e-3_r8 / (1._r8 - fbw(ivt))

!            else
!               leaf_biomass(m) = 0_r8
!            end if

            if (woody(ivt) == 1._r8) then

               ! trees and shrubs for now have a very simple allometry, with hard-wired
               ! stem taper (height:radius) and nstem from PFT parameter file
!               if (use_cndv) then

!                  if (fpcgrid(m) > 0._r8 .and. nind(m) > 0._r8) then

!                     stocking = nind(m)/fpcgrid(m) !#ind/m2 nat veg area -> #ind/m2 patch area
!                     htop(m) = allom2(ivt) * ( (24._r8 * deadstemc(m) / &
!                          (SHR_CONST_PI * stocking * dwood(ivt) * taper(ivt)))**(1._r8/3._r8) )**allom3(ivt) ! lpj's htop w/ cn's stemdiam

!                  else
!                     htop(m) = 0._r8
!                  end if

!               else
                  !correct height calculation if doing accelerated spinup
!                  htop_p(m) = ((3._r8 * deadstemc_p(m) * taper(ivt) * taper(ivt))/ &
!                         (SHR_CONST_PI * nstem(ivt) * dwood(ivt)))**(1._r8/3._r8)

!               endif

!               if (use_biomass_heat_storage) then
!                  ! Assumes fbw (fraction of biomass that is water) is the same for leaves and stems
!                  stem_biomass(m) = (spinup_factor_deadwood*deadstemc(m) + livestemc(m)) &
!                       * c_to_b * 1.e-3_r8 / (1._r8 - fbw(ivt))
!               else
!                  stem_biomass(m) = 0_r8
!               end if

               !
               ! Peter Thornton, 5/3/2004
               ! Adding test to keep htop from getting too close to forcing height for windspeed
               ! Also added for grass, below, although it is not likely to ever be an issue.
!               htop_p(m) = min(htop_p(m),(forc_hgt_u(i)/(displar(ivt)+z0mr(ivt)))-3._r8)

               ! Peter Thornton, 8/11/2004
               ! Adding constraint to keep htop from going to 0.0.
               ! This becomes an issue when fire mortality is pushing deadstemc
               ! to 0.0.
!               htop_p(m) = max(htop_p(m), 0.01_r8)

!               hbot_p(m) = max(0._r8, min(3._r8, htop_p(m)-1._r8))

            else if (ivt >= npcropmin) then ! prognostic crops

               if (tlai_p(m) >= laimx(ivt)) peaklai_p(m) = 1 ! used in CNAllocation

               if (ivt == ntmp_corn .or. ivt == nirrig_tmp_corn .or. &
                   ivt == ntrp_corn .or. ivt == nirrig_trp_corn .or. &
                   ivt == nsugarcane .or. ivt == nirrig_sugarcane .or. &
                   ivt == nmiscanthus .or. ivt == nirrig_miscanthus .or. &
                   ivt == nswitchgrass .or. ivt == nirrig_switchgrass) then
                  tsai_p(m) = 0.1_r8 * tlai_p(m)
               else
                  tsai_p(m) = 0.2_r8 * tlai_p(m)
               end if

               ! "stubble" after harvest
               if (harvdate_p(m) < 999 .and. tlai_p(m) == 0._r8) then
                  tsai_p(m) = 0.25_r8*(1._r8-farea_burned(i)*0.90_r8)    !changed by F. Li and S. Levis
!                  htmx_p(m) = 0._r8
                  peaklai_p(m) = 0
               end if
               !if (harvdate(m) < 999 .and. tlai(m) > 0._r8) write(iulog,*) 'CNVegStructUpdate: tlai>0 after harvest!' ! remove after initial debugging?

               ! canopy top and bottom heights
!               htop_p(m) = ztopmx(ivt) * (min(tlai_p(m)/(laimx(ivt)-1._r8),1._r8))**2
!               htmx_p(m) = max(htmx_p(m), htop_p(m))
!               htop_p(m) = max(0.05_r8, max(htmx_p(m),htop_p(m)))
!               hbot_p(m) = 0.02_r8

            else ! generic crops and ...

               ! grasses

               ! height for grasses depends only on LAI
!               htop_p(m) = max(0.25_r8, tlai_p(m) * 0.25_r8)

!               htop_p(m) = min(htop_p(m),(forc_hgt_u(i)/(displar(ivt)+z0mr(ivt)))-3._r8)

               ! Peter Thornton, 8/11/2004
               ! Adding constraint to keep htop from going to 0.0.
!               htop_p(m) = max(htop_p(m), 0.01_r8)

!               hbot_p(m) = max(0.0_r8, min(0.05_r8, htop_p(m)-0.20_r8))
            end if

         else

!            tlai(m) = 0._r8
!            tsai(m) = 0._r8
!            htop(m) = 0._r8
!            hbot(m) = 0._r8

         end if

         ! adjust lai and sai for burying by snow. 
         ! snow burial fraction for short vegetation (e.g. grasses, crops) changes with vegetation height 
         ! accounts for a 20% bending factor, as used in Lombardozzi et al. (2018) GRL 45(18), 9889-9897

         ! NOTE: The following snow burial code is duplicated in SatellitePhenologyMod.
         ! Changes in one place should be accompanied by similar changes in the other.

!         if (ivt > noveg .and. ivt <= nbrdlf_dcd_brl_shrub ) then
!            ol = min( max(snowdp(i)-hbot_p(m), 0._r8), htop_p(m)-hbot_p(m))
!            fb = 1._r8 - ol / max(1.e-06_r8, htop_p(m)-hbot_p(m))
!         else
!            fb = 1._r8 - (max(min(snowdp(i),max(0.05,htop_p(m)*0.8_r8)),0._r8)/(max(0.05,htop_p(m)*0.8_r8)))
!            !depth of snow required for complete burial of grasses
!         endif
!
!         if (fsno(i) <= frac_sno_threshold) then
!            frac_sno_adjusted = fsno(i)
!         else
!            ! avoid tiny but non-zero elai and esai that can cause radiation and/or photosynthesis code to blow up
!            frac_sno_adjusted = 1._r8
!         end if

!         elai(m) = max(tlai(m)*(1.0_r8 - frac_sno_adjusted) + tlai(m)*fb*frac_sno_adjusted, 0.0_r8)
!         esai(m) = max(tsai(m)*(1.0_r8 - frac_sno_adjusted) + tsai(m)*fb*frac_sno_adjusted, 0.0_r8)

         ! Fraction of vegetation free of snow
!         if ((elai(m) + esai(m)) > 0._r8) then
!            frac_veg_nosno_alb(m) = 1
!         else
!            frac_veg_nosno_alb(m) = 0
!         end if

      end do

 end subroutine CNVegStructUpdate

end module bgc_veg_CNVegStructUpdateMod
#endif
