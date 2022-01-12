module bgc_CNCStateUpdate1Mod

use precision
use MOD_PFTimeInvars, only: pftclass
use PFT_Const, only: woody
use MOD_TimeInvariants, only: &
  ! bgc constants
           donor_pool, receiver_pool, i_met_lit, i_cel_lit, i_lig_lit, i_cwd, i_soil1, i_soil2, i_soil3

use MOD_1D_Fluxes, only: &
! decomposition pools flux varables (in)
           decomp_cpools_sourcesink, decomp_ctransfer_vr, decomp_hr_vr      , &
           phenology_to_met_c      , phenology_to_cel_c , phenology_to_lig_c

use MOD_PFTimeVars, only: &
! vegetation carbon state variables (inout)
           leafc_p            , leafc_storage_p     , leafc_xfer_p     , &
           frootc_p           , frootc_storage_p    , frootc_xfer_p    , &
           livestemc_p        , livestemc_storage_p , livestemc_xfer_p , &
           deadstemc_p        , deadstemc_storage_p , deadstemc_xfer_p , &
           livecrootc_p       , livecrootc_storage_p, livecrootc_xfer_p, &
           deadcrootc_p       , deadcrootc_storage_p, deadcrootc_xfer_p, &
           grainc_p           , grainc_storage_p    , grainc_xfer_p    , &
           cropseedc_deficit_p, xsmrpool_p          , gresp_storage_p  , gresp_xfer_p, &    

! crop variables (in)
           harvdate_p

use MOD_1D_PFTFluxes, only: &
! vegetation carbon flux variables (in)
! xfer to display
           leafc_xfer_to_leafc_p          , frootc_xfer_to_frootc_p        , &
           livestemc_xfer_to_livestemc_p  , deadstemc_xfer_to_deadstemc_p  , &
           livecrootc_xfer_to_livecrootc_p, deadcrootc_xfer_to_deadcrootc_p, &
           grainc_xfer_to_grainc_p        , &

! storage to xfer (in)
           leafc_storage_to_xfer_p     , frootc_storage_to_xfer_p    , &
           livestemc_storage_to_xfer_p , deadstemc_storage_to_xfer_p , &
           livecrootc_storage_to_xfer_p, deadcrootc_storage_to_xfer_p, &
           grainc_storage_to_xfer_p    , gresp_storage_to_xfer_p     , &

! display to litter & live to dead (in)
           leafc_to_litter_p       , frootc_to_litter_p        , &
           grainc_to_food_p        , grainc_to_seed_p          , &
           crop_seedc_to_leaf_p    , livestemc_to_litter_p     , &
           livestemc_to_deadstemc_p, livecrootc_to_deadcrootc_p, &

! cpool to display/storage (in)
           cpool_to_xsmrpool_p  , cpool_to_gresp_storage_p     , &
           cpool_to_leafc_p     , cpool_to_leafc_storage_p     , &
           cpool_to_frootc_p    , cpool_to_frootc_storage_p    , &
           cpool_to_livestemc_p , cpool_to_livestemc_storage_p , &
           cpool_to_deadstemc_p , cpool_to_deadstemc_storage_p , &
           cpool_to_livecrootc_p, cpool_to_livecrootc_storage_p, &
           cpool_to_deadcrootc_p, cpool_to_deadcrootc_storage_p, &
           cpool_to_grainc_p    , cpool_to_grainc_storage_p    , &

! maintenance respiration fluxes (in)
           leaf_xsmr_p,   froot_xsmr_p,   livestem_xsmr_p,   livecroot_xsmr_p,   grain_xsmr_p,   &

! growth respiration fluxes (in/inout)
           transfer_leaf_gr_p     , transfer_froot_gr_p    , &
           transfer_livestem_gr_p , transfer_deadstem_gr_p , &
           transfer_livecroot_gr_p, transfer_deadcroot_gr_p, &
           transfer_grain_gr_p    , xsmrpool_to_atm_p      

    

implicit none

public :: CStateUpdate0
public :: CStateUpdate1

contains

subroutine CStateUpdate0()

end subroutine CStateUpdate0

subroutine CStateUpdate1 &
! model configurations
          (i, ps, pe, deltim, nl_soil, ndecomp_transitions, npcropmin)

integer ,intent(in) :: i
integer ,intent(in) :: ps
integer ,intent(in) :: pe
real(r8),intent(in) :: deltim
integer ,intent(in) :: nl_soil
integer ,intent(in) :: ndecomp_transitions
integer ,intent(in) :: npcropmin

integer j,k
integer ivt, m

   !print*,'leafc,leafc_storage,leafc_xfer',leafc(i),leafc_storage(i),leafc_xfer(i)
   !print*,'frootc,frootc_storage,frootc_xfer',frootc(i),frootc_storage(i),frootc_xfer(i)
   !print*,'livestemc,livestemc_storage,livestemc_xfer',livestemc(i),livestemc_storage(i),livestemc_xfer(i)
   !print*,'deadstemc,deadstemc_storage,deadstemc_xfer',deadstemc(i),deadstemc_storage(i),deadstemc_xfer(i)
   !print*,'livecrootc,livecrootc_storage,livecrootc_xfer',livecrootc(i),livecrootc_storage(i),livecrootc_xfer(i)
   !print*,'deadcrootc,deadcrootc_storage,deadcrootc_xfer',deadcrootc(i),deadcrootc_storage(i),deadcrootc_xfer(i)

   do j=1,nl_soil
      decomp_cpools_sourcesink(j,i_met_lit,i) = phenology_to_met_c(j,i) *deltim
      decomp_cpools_sourcesink(j,i_cel_lit,i) = phenology_to_cel_c(j,i) *deltim
      decomp_cpools_sourcesink(j,i_lig_lit,i) = phenology_to_lig_c(j,i) *deltim
      decomp_cpools_sourcesink(j,i_cwd    ,i) = 0._r8
      !print*,'phenology source to cel',j,phenology_to_cel_c(j,i) *deltim
   end do

   do k = 1, ndecomp_transitions
      do j = 1,nl_soil
         decomp_cpools_sourcesink(j,donor_pool(k),i) = &
              decomp_cpools_sourcesink(j,donor_pool(k),i) &
              - (decomp_hr_vr(j,k,i) + decomp_ctransfer_vr(j,k,i)) * deltim
      end do
   end do 

   do k = 1,ndecomp_transitions
      if ( receiver_pool(k) /= 0 ) then  ! skip terminal transitions
         do j = 1,nl_soil
            decomp_cpools_sourcesink(j,receiver_pool(k),i) = &
                 decomp_cpools_sourcesink(j,receiver_pool(k),i) &
                 + decomp_ctransfer_vr(j,k,i) * deltim
         end do
      end if
   end do

   do m = ps , pe
      ivt = pftclass(m)
   !print*,'leafc_xfer_to_leafc',leafc_xfer_to_leafc  (m) * deltim
      leafc_p      (m) = leafc_p      (m) + leafc_xfer_to_leafc_p  (m) * deltim
      leafc_xfer_p (m) = leafc_xfer_p (m) - leafc_xfer_to_leafc_p  (m) * deltim
      frootc_p     (m) = frootc_p     (m) + frootc_xfer_to_frootc_p(m) * deltim
      frootc_xfer_p(m) = frootc_xfer_p(m) - frootc_xfer_to_frootc_p(m) * deltim
      if(woody(ivt) == 1)then
         livestemc_p      (m) = livestemc_p      (m) + livestemc_xfer_to_livestemc_p  (m) * deltim
         livestemc_xfer_p (m) = livestemc_xfer_p (m) - livestemc_xfer_to_livestemc_p  (m) * deltim
         deadstemc_p      (m) = deadstemc_p      (m) + deadstemc_xfer_to_deadstemc_p  (m) * deltim
         deadstemc_xfer_p (m) = deadstemc_xfer_p (m) - deadstemc_xfer_to_deadstemc_p  (m) * deltim
         livecrootc_p     (m) = livecrootc_p     (m) + livecrootc_xfer_to_livecrootc_p(m) * deltim
         livecrootc_xfer_p(m) = livecrootc_xfer_p(m) - livecrootc_xfer_to_livecrootc_p(m) * deltim
         deadcrootc_p     (m) = deadcrootc_p     (m) + deadcrootc_xfer_to_deadcrootc_p(m) * deltim
         deadcrootc_xfer_p(m) = deadcrootc_xfer_p(m) - deadcrootc_xfer_to_deadcrootc_p(m) * deltim
      end if
      if (ivt >= npcropmin) then ! skip 2 generic crops
        ! lines here for consistency; the transfer terms are zero
         livestemc_p     (m)  = livestemc_p     (m) + livestemc_xfer_to_livestemc_p(m) * deltim
         livestemc_xfer_p(m)  = livestemc_xfer_p(m) - livestemc_xfer_to_livestemc_p(m) * deltim
         grainc_p        (m)  = grainc_p        (m) + grainc_xfer_to_grainc_p      (m) * deltim
         grainc_xfer_p   (m)  = grainc_xfer_p   (m) - grainc_xfer_to_grainc_p      (m) * deltim
      end if

   !print*,'leafc_xfer_to_leafc',leafc_to_litter  (m) * deltim
 ! phenology: litterfall fluxes
      leafc_p (m) = leafc_p (m) - leafc_to_litter_p (m) * deltim
      frootc_p(m) = frootc_p(m) - frootc_to_litter_p(m) * deltim

   ! livewood turnover fluxes
      if (woody(ivt) == 1) then
         livestemc_p (m) = livestemc_p (m) - livestemc_to_deadstemc_p  (m) * deltim
         deadstemc_p (m) = deadstemc_p (m) + livestemc_to_deadstemc_p  (m) * deltim
         livecrootc_p(m) = livecrootc_p(m) - livecrootc_to_deadcrootc_p(m) * deltim
         deadcrootc_p(m) = deadcrootc_p(m) + livecrootc_to_deadcrootc_p(m) * deltim
      end if
      if (ivt >= npcropmin) then ! skip 2 generic crops
         livestemc_p        (m) = livestemc_p        (m) - livestemc_to_litter_p(m) * deltim
         grainc_p           (m) = grainc_p           (m) - (grainc_to_food_p(m) + grainc_to_seed_p(m)) * deltim
         cropseedc_deficit_p(m) = cropseedc_deficit_p(m) - crop_seedc_to_leaf_p(m) * deltim + grainc_to_seed_p(m) * deltim
      end if

   ! maintenance respiration fluxes from xsmrpool
      xsmrpool_p(m) = xsmrpool_p(m) + cpool_to_xsmrpool_p(m) * deltim
      xsmrpool_p(m) = xsmrpool_p(m) - leaf_xsmr_p        (m) * deltim
      xsmrpool_p(m) = xsmrpool_p(m) - froot_xsmr_p       (m) * deltim
      if (woody(ivt) == 1) then
         xsmrpool_p(m) = xsmrpool_p(m) - livestem_xsmr_p (m) * deltim
         xsmrpool_p(m) = xsmrpool_p(m) - livecroot_xsmr_p(m) * deltim
      end if
      !print*,'cpool_to_leafc in CNCStateUpdate1',i,cpool_to_leafc(m) * deltim
   !print*,'cpool_to_leafc_storage in CNCStateUpdate1',i,cpool_to_leafc_storage(m) * deltim
      leafc_p         (m) = leafc_p         (m) + cpool_to_leafc_p         (m) * deltim
      leafc_storage_p (m) = leafc_storage_p (m) + cpool_to_leafc_storage_p (m) * deltim
      frootc_p        (m) = frootc_p        (m) + cpool_to_frootc_p        (m) * deltim
      frootc_storage_p(m) = frootc_storage_p(m) + cpool_to_frootc_storage_p(m) * deltim
      if (woody(ivt) == 1) then
         livestemc_p         (m) = livestemc_p         (m) + cpool_to_livestemc_p         (m) * deltim
         livestemc_storage_p (m) = livestemc_storage_p (m) + cpool_to_livestemc_storage_p (m) * deltim
         deadstemc_p         (m) = deadstemc_p         (m) + cpool_to_deadstemc_p         (m) * deltim
         deadstemc_storage_p (m) = deadstemc_storage_p (m) + cpool_to_deadstemc_storage_p (m) * deltim
         livecrootc_p        (m) = livecrootc_p        (m) + cpool_to_livecrootc_p        (m) * deltim
         livecrootc_storage_p(m) = livecrootc_storage_p(m) + cpool_to_livecrootc_storage_p(m) * deltim
         deadcrootc_p        (m) = deadcrootc_p        (m) + cpool_to_deadcrootc_p        (m) * deltim
         deadcrootc_storage_p(m) = deadcrootc_storage_p(m) + cpool_to_deadcrootc_storage_p(m) * deltim
      end if
      if (ivt >= npcropmin) then ! skip 2 generic crops
         livestemc_p        (m) = livestemc_p        (m) + cpool_to_livestemc_p        (m) * deltim
         livestemc_storage_p(m) = livestemc_storage_p(m) + cpool_to_livestemc_storage_p(m) * deltim
         grainc_p           (m) = grainc_p           (m) + cpool_to_grainc_p           (m) * deltim
         grainc_storage_p   (m) = grainc_storage_p   (m) + cpool_to_grainc_storage_p   (m) * deltim
      end if
      ! growth respiration for transfer growth
      gresp_xfer_p(m) = gresp_xfer_p(m) - transfer_leaf_gr_p (m) * deltim
      gresp_xfer_p(m) = gresp_xfer_p(m) - transfer_froot_gr_p(m) * deltim
      if (woody(ivt) == 1) then
         gresp_xfer_p(m) = gresp_xfer_p(m) - transfer_livestem_gr_p (m) * deltim
         gresp_xfer_p(m) = gresp_xfer_p(m) - transfer_deadstem_gr_p (m) * deltim
         gresp_xfer_p(m) = gresp_xfer_p(m) - transfer_livecroot_gr_p(m) * deltim
         gresp_xfer_p(m) = gresp_xfer_p(m) - transfer_deadcroot_gr_p(m) * deltim
      end if
      if (ivt >= npcropmin) then ! skip 2 generic crops
         gresp_xfer_p(m) = gresp_xfer_p(m) - transfer_livestem_gr_p(m) * deltim
         gresp_xfer_p(m) = gresp_xfer_p(m) - transfer_grain_gr_p   (m) * deltim
      end if
   ! growth respiration stored for release during transfer growth
      gresp_storage_p(m) = gresp_storage_p(m) + cpool_to_gresp_storage_p(m) * deltim

   !print*,'leaf storage to xfer',leafc_storage_to_xfer (m) * deltim
   ! move storage pools into transfer pools
      leafc_storage_p (m) = leafc_storage_p (m) - leafc_storage_to_xfer_p (m) * deltim
      leafc_xfer_p    (m) = leafc_xfer_p    (m) + leafc_storage_to_xfer_p (m) * deltim
      frootc_storage_p(m) = frootc_storage_p(m) - frootc_storage_to_xfer_p(m) * deltim
      frootc_xfer_p   (m) = frootc_xfer_p   (m) + frootc_storage_to_xfer_p(m) * deltim
      if (woody(ivt) == 1) then
         gresp_storage_p     (m) = gresp_storage_p    (m) - gresp_storage_to_xfer_p(m) * deltim
         gresp_xfer_p        (m) = gresp_xfer_p       (m) + gresp_storage_to_xfer_p(m) * deltim
     
         livestemc_storage_p (m) = livestemc_storage_p (m) - livestemc_storage_to_xfer_p (m) * deltim
         livestemc_xfer_p    (m) = livestemc_xfer_p    (m) + livestemc_storage_to_xfer_p (m) * deltim
         deadstemc_storage_p (m) = deadstemc_storage_p (m) - deadstemc_storage_to_xfer_p (m) * deltim
         deadstemc_xfer_p    (m) = deadstemc_xfer_p    (m) + deadstemc_storage_to_xfer_p (m) * deltim
         livecrootc_storage_p(m) = livecrootc_storage_p(m) - livecrootc_storage_to_xfer_p(m) * deltim
         livecrootc_xfer_p   (m) = livecrootc_xfer_p   (m) + livecrootc_storage_to_xfer_p(m) * deltim
         deadcrootc_storage_p(m) = deadcrootc_storage_p(m) - deadcrootc_storage_to_xfer_p(m) * deltim
         deadcrootc_xfer_p   (m) = deadcrootc_xfer_p   (m) + deadcrootc_storage_to_xfer_p(m) * deltim
      end if
      if (ivt >= npcropmin) then ! skip 2 generic crops
      ! lines here for consistency; the transfer terms are zero
         livestemc_storage_p (m) = livestemc_storage_p(m) - livestemc_storage_to_xfer_p(m) * deltim
         livestemc_xfer_p    (m) = livestemc_xfer_p   (m) + livestemc_storage_to_xfer_p(m) * deltim
         grainc_storage_p    (m) = grainc_storage_p   (m) - grainc_storage_to_xfer_p   (m) * deltim
         grainc_xfer_p       (m) = grainc_xfer_p      (m) + grainc_storage_to_xfer_p   (m) * deltim
      end if
      if (ivt >= npcropmin) then ! skip 2 generic crops
         xsmrpool_p(m) = xsmrpool_p(m) - livestem_xsmr_p(m)*deltim
         xsmrpool_p(m) = xsmrpool_p(m) - grain_xsmr_p   (m)*deltim
         if (harvdate_p(m) < 999) then ! beginning at harvest, send to atm
         ! TODO (mv, 11-02-2014) the following lines are why the cf_veg is
         ! an intent(inout)
         ! fluxes should not be updated in this module - not sure where
         ! this belongs
         ! DML (06-20-2017) While debugging crop isotope code, found that cpool_patch and frootc_patch 
         ! could occasionally be very small but nonzero numbers after crop harvest, which persists 
         ! through to next planting and for reasons that could not 100%
         ! isolate, caused C12/C13 ratios to occasionally go out of
         ! bounds. Zeroing out these small pools and putting them into the flux to the
         ! atmosphere solved many of the crop isotope problems

            xsmrpool_to_atm_p(m) = xsmrpool_to_atm_p(m) + xsmrpool_p(m)/deltim
            xsmrpool_p       (m) = 0._r8
!         xsmrpool_to_atm(m) = xsmrpool_to_atm(m) + cpool(m)/deltim
            xsmrpool_to_atm_p(m) = xsmrpool_to_atm_p(m) + frootc_p(m)/deltim
            frootc_p         (m) = 0._r8
         end if
      end if
      if(m .eq. 22784)print*,'leafc',i,m,ivt,leafc_p(m)+leafc_storage_p(m)+leafc_xfer_p(m)
      if(m .eq. 22784)print*,'frootc',frootc_p(m)+frootc_storage_p(m)+frootc_xfer_p(m)
      if(m .eq. 22784)print*,'livestemc',livestemc_p(m)+livestemc_storage_p(m)+livestemc_xfer_p(m)
      if(m .eq. 22784)print*,'deadstemc',deadstemc_p(m)+deadstemc_storage_p(m)+deadstemc_xfer_p(m)
      if(m .eq. 22784)print*,'livecrootc',livecrootc_p(m)+livecrootc_storage_p(m)+livecrootc_xfer_p(m)
      if(m .eq. 22784)print*,'deadcrootc',deadcrootc_p(m)+deadcrootc_storage_p(m)+deadcrootc_xfer_p(m)
   end do ! end pft loop
    
end subroutine CStateUpdate1

end module bgc_CNCStateUpdate1Mod
