#include <define.h>

module bgc_CNCStateUpdate1Mod

use precision
use MOD_PFTimeInvars, only: pftclass, pftfrac
use PFT_Const, only: woody
use MOD_TimeInvariants, only: &
  ! bgc constants
           donor_pool, receiver_pool, i_met_lit, i_cel_lit, i_lig_lit, i_cwd, i_soil1, i_soil2, i_soil3

use MOD_TimeVariables, only: &
           I_met_c_vr_acc, I_cel_c_vr_acc, I_lig_c_vr_acc, &
           AKX_met_to_soil1_c_vr_acc  , AKX_cel_to_soil1_c_vr_acc  , AKX_lig_to_soil2_c_vr_acc  , AKX_soil1_to_soil2_c_vr_acc, &
           AKX_cwd_to_cel_c_vr_acc    , AKX_cwd_to_lig_c_vr_acc    , AKX_soil1_to_soil3_c_vr_acc, AKX_soil2_to_soil1_c_vr_acc, &
           AKX_soil2_to_soil3_c_vr_acc, AKX_soil3_to_soil1_c_vr_acc, &
           AKX_met_exit_c_vr_acc      , AKX_cel_exit_c_vr_acc      , AKX_lig_exit_c_vr_acc      , AKX_cwd_exit_c_vr_acc      , &
           AKX_soil1_exit_c_vr_acc    , AKX_soil2_exit_c_vr_acc    , AKX_soil3_exit_c_vr_acc 
   

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
           cpool_p, &

! crop variables (in)
           harvdate_p         , cropprod1c_p        , &

! SASU variables           
           I_leafc_p_acc     , I_leafc_st_p_acc     , I_frootc_p_acc    , I_frootc_st_p_acc    , &
           I_livestemc_p_acc , I_livestemc_st_p_acc , I_deadstemc_p_acc , I_deadstemc_st_p_acc , &
           I_livecrootc_p_acc, I_livecrootc_st_p_acc, I_deadcrootc_p_acc, I_deadcrootc_st_p_acc, &
           I_grainc_p_acc    , I_grainc_st_p_acc    , &

           AKX_leafc_xf_to_leafc_p_acc           , AKX_frootc_xf_to_frootc_p_acc           , AKX_livestemc_xf_to_livestemc_p_acc     , &
           AKX_deadstemc_xf_to_deadstemc_p_acc   , AKX_livecrootc_xf_to_livecrootc_p_acc   , AKX_deadcrootc_xf_to_deadcrootc_p_acc   , &
           AKX_grainc_xf_to_grainc_p_acc         , AKX_livestemc_to_deadstemc_p_acc        , AKX_livecrootc_to_deadcrootc_p_acc      , &
           
           AKX_leafc_st_to_leafc_xf_p_acc        , AKX_frootc_st_to_frootc_xf_p_acc        , AKX_livestemc_st_to_livestemc_xf_p_acc  , &
           AKX_deadstemc_st_to_deadstemc_xf_p_acc, AKX_livecrootc_st_to_livecrootc_xf_p_acc, AKX_deadcrootc_st_to_deadcrootc_xf_p_acc, &
           AKX_livestemc_st_to_livestemc_xf_p_acc, AKX_grainc_st_to_grainc_xf_p_acc        , &

           AKX_leafc_exit_p_acc                  , AKX_frootc_exit_p_acc                   , AKX_livestemc_exit_p_acc                , &
           AKX_deadstemc_exit_p_acc              , AKX_livecrootc_exit_p_acc               , AKX_deadcrootc_exit_p_acc               , &
           AKX_grainc_exit_p_acc                 , &

           AKX_leafc_st_exit_p_acc               , AKX_frootc_st_exit_p_acc                , AKX_livestemc_st_exit_p_acc             , &
           AKX_deadstemc_st_exit_p_acc           , AKX_livecrootc_st_exit_p_acc            , AKX_deadcrootc_st_exit_p_acc            , &
           AKX_grainc_st_exit_p_acc              , &

           AKX_leafc_xf_exit_p_acc               , AKX_frootc_xf_exit_p_acc                , AKX_livestemc_xf_exit_p_acc             , &
           AKX_deadstemc_xf_exit_p_acc           , AKX_livecrootc_xf_exit_p_acc            , AKX_deadcrootc_xf_exit_p_acc            , &
           AKX_grainc_xf_exit_p_acc
           

use MOD_1D_PFTFluxes, only: &
! vegetation carbon flux variables (in)
! Vegetation physiology
           psn_to_cpool_p, &

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

! crop
           cropprod1c_loss_p, &

! cpool to display/storage (in)
           cpool_to_xsmrpool_p  , cpool_to_gresp_storage_p     , &
           cpool_to_leafc_p     , cpool_to_leafc_storage_p     , &
           cpool_to_frootc_p    , cpool_to_frootc_storage_p    , &
           cpool_to_livestemc_p , cpool_to_livestemc_storage_p , &
           cpool_to_deadstemc_p , cpool_to_deadstemc_storage_p , &
           cpool_to_livecrootc_p, cpool_to_livecrootc_storage_p, &
           cpool_to_deadcrootc_p, cpool_to_deadcrootc_storage_p, &
           cpool_to_grainc_p    , cpool_to_grainc_storage_p    , &

! cpool to growth repsiration
           cpool_leaf_gr_p             , cpool_froot_gr_p             , &
           cpool_livestem_gr_p         , cpool_deadstem_gr_p          , &
           cpool_livecroot_gr_p        , cpool_deadcroot_gr_p         , &
           cpool_leaf_storage_gr_p     , cpool_froot_storage_gr_p     , &
           cpool_livestem_storage_gr_p , cpool_deadstem_storage_gr_p  , &
           cpool_livecroot_storage_gr_p, cpool_deadcroot_storage_gr_p , &

           cpool_grain_gr_p            , cpool_grain_storage_gr_p     , &

! maintenance respiration fluxes (in)
           leaf_xsmr_p,  froot_xsmr_p,  livestem_xsmr_p,  livecroot_xsmr_p,  grain_xsmr_p , &
           leaf_curmr_p, froot_curmr_p, livestem_curmr_p, livecroot_curmr_p, grain_curmr_p, &

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

   do m = ps, pe
      cpool_p(m) = cpool_p(m) + psn_to_cpool_p(m) * deltim
   end do
   do j=1,nl_soil
      decomp_cpools_sourcesink(j,i_met_lit,i) = phenology_to_met_c(j,i) *deltim
      decomp_cpools_sourcesink(j,i_cel_lit,i) = phenology_to_cel_c(j,i) *deltim
      decomp_cpools_sourcesink(j,i_lig_lit,i) = phenology_to_lig_c(j,i) *deltim
      decomp_cpools_sourcesink(j,i_cwd    ,i) = 0._r8
      !print*,'phenology source to cel',j,phenology_to_cel_c(j,i) *deltim
   end do

#ifdef SASU
   do j=1,nl_soil
      I_met_c_vr_acc(j,i) = I_met_c_vr_acc(j,i) + phenology_to_met_c(j,i) *deltim
      I_cel_c_vr_acc(j,i) = I_cel_c_vr_acc(j,i) + phenology_to_cel_c(j,i) *deltim
      I_lig_c_vr_acc(j,i) = I_lig_c_vr_acc(j,i) + phenology_to_lig_c(j,i) *deltim
   end do
#endif

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

#ifdef SASU
   do j = 1, nl_soil
      AKX_met_to_soil1_c_vr_acc  (j,i) = AKX_met_to_soil1_c_vr_acc  (j,i) + decomp_ctransfer_vr(j, 1,i) * deltim
      AKX_cel_to_soil1_c_vr_acc  (j,i) = AKX_cel_to_soil1_c_vr_acc  (j,i) + decomp_ctransfer_vr(j, 2,i) * deltim
      AKX_lig_to_soil2_c_vr_acc  (j,i) = AKX_lig_to_soil2_c_vr_acc  (j,i) + decomp_ctransfer_vr(j, 3,i) * deltim
      AKX_soil1_to_soil2_c_vr_acc(j,i) = AKX_soil1_to_soil2_c_vr_acc(j,i) + decomp_ctransfer_vr(j, 4,i) * deltim
      AKX_cwd_to_cel_c_vr_acc    (j,i) = AKX_cwd_to_cel_c_vr_acc    (j,i) + decomp_ctransfer_vr(j, 5,i) * deltim
      AKX_cwd_to_lig_c_vr_acc    (j,i) = AKX_cwd_to_lig_c_vr_acc    (j,i) + decomp_ctransfer_vr(j, 6,i) * deltim
      AKX_soil1_to_soil3_c_vr_acc(j,i) = AKX_soil1_to_soil3_c_vr_acc(j,i) + decomp_ctransfer_vr(j, 7,i) * deltim
      AKX_soil2_to_soil1_c_vr_acc(j,i) = AKX_soil2_to_soil1_c_vr_acc(j,i) + decomp_ctransfer_vr(j, 8,i) * deltim
      AKX_soil2_to_soil3_c_vr_acc(j,i) = AKX_soil2_to_soil3_c_vr_acc(j,i) + decomp_ctransfer_vr(j, 9,i) * deltim
      AKX_soil3_to_soil1_c_vr_acc(j,i) = AKX_soil3_to_soil1_c_vr_acc(j,i) + decomp_ctransfer_vr(j,10,i) * deltim

      AKX_met_exit_c_vr_acc      (j,i) = AKX_met_exit_c_vr_acc      (j,i) + (decomp_hr_vr(j, 1,i) + decomp_ctransfer_vr(j, 1,i)) * deltim
      AKX_cel_exit_c_vr_acc      (j,i) = AKX_cel_exit_c_vr_acc      (j,i) + (decomp_hr_vr(j, 2,i) + decomp_ctransfer_vr(j, 2,i)) * deltim
      AKX_lig_exit_c_vr_acc      (j,i) = AKX_lig_exit_c_vr_acc      (j,i) + (decomp_hr_vr(j, 3,i) + decomp_ctransfer_vr(j, 3,i)) * deltim
      AKX_soil1_exit_c_vr_acc    (j,i) = AKX_soil1_exit_c_vr_acc    (j,i) + (decomp_hr_vr(j, 4,i) + decomp_ctransfer_vr(j, 4,i)) * deltim
      AKX_cwd_exit_c_vr_acc      (j,i) = AKX_cwd_exit_c_vr_acc      (j,i) + (decomp_hr_vr(j, 5,i) + decomp_ctransfer_vr(j, 5,i)) * deltim
      AKX_cwd_exit_c_vr_acc      (j,i) = AKX_cwd_exit_c_vr_acc      (j,i) + (decomp_hr_vr(j, 6,i) + decomp_ctransfer_vr(j, 6,i)) * deltim
      AKX_soil1_exit_c_vr_acc    (j,i) = AKX_soil1_exit_c_vr_acc    (j,i) + (decomp_hr_vr(j, 7,i) + decomp_ctransfer_vr(j, 7,i)) * deltim
      AKX_soil2_exit_c_vr_acc    (j,i) = AKX_soil2_exit_c_vr_acc    (j,i) + (decomp_hr_vr(j, 8,i) + decomp_ctransfer_vr(j, 8,i)) * deltim
      AKX_soil2_exit_c_vr_acc    (j,i) = AKX_soil2_exit_c_vr_acc    (j,i) + (decomp_hr_vr(j, 9,i) + decomp_ctransfer_vr(j, 9,i)) * deltim
      AKX_soil3_exit_c_vr_acc    (j,i) = AKX_soil3_exit_c_vr_acc    (j,i) + (decomp_hr_vr(j,10,i) + decomp_ctransfer_vr(j,10,i)) * deltim
   end do
#endif

   do m = ps , pe
      ivt = pftclass(m)
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

#ifdef SASU
      AKX_leafc_xf_to_leafc_p_acc  (m) = AKX_leafc_xf_to_leafc_p_acc  (m) + leafc_xfer_to_leafc_p  (m) * deltim
      AKX_frootc_xf_to_frootc_p_acc(m) = AKX_frootc_xf_to_frootc_p_acc(m) + frootc_xfer_to_frootc_p(m) * deltim
      AKX_leafc_xf_exit_p_acc      (m) = AKX_leafc_xf_exit_p_acc      (m) + leafc_xfer_to_leafc_p  (m) * deltim
      AKX_frootc_xf_exit_p_acc     (m) = AKX_frootc_xf_exit_p_acc     (m) + frootc_xfer_to_frootc_p(m) * deltim
      if(woody(ivt) == 1)then
         AKX_livestemc_xf_to_livestemc_p_acc  (m) = AKX_livestemc_xf_to_livestemc_p_acc  (m) + livestemc_xfer_to_livestemc_p  (m) * deltim
         AKX_livestemc_xf_exit_p_acc          (m) = AKX_livestemc_xf_exit_p_acc          (m) + livestemc_xfer_to_livestemc_p  (m) * deltim
         AKX_deadstemc_xf_to_deadstemc_p_acc  (m) = AKX_deadstemc_xf_to_deadstemc_p_acc  (m) + deadstemc_xfer_to_deadstemc_p  (m) * deltim
         AKX_deadstemc_xf_exit_p_acc          (m) = AKX_deadstemc_xf_exit_p_acc          (m) + deadstemc_xfer_to_deadstemc_p  (m) * deltim
         AKX_livecrootc_xf_to_livecrootc_p_acc(m) = AKX_livecrootc_xf_to_livecrootc_p_acc(m) + livecrootc_xfer_to_livecrootc_p(m) * deltim
         AKX_livecrootc_xf_exit_p_acc         (m) = AKX_livecrootc_xf_exit_p_acc         (m) + livecrootc_xfer_to_livecrootc_p(m) * deltim
         AKX_deadcrootc_xf_to_deadcrootc_p_acc(m) = AKX_deadcrootc_xf_to_deadcrootc_p_acc(m) + deadcrootc_xfer_to_deadcrootc_p(m) * deltim
         AKX_deadcrootc_xf_exit_p_acc         (m) = AKX_deadcrootc_xf_exit_p_acc         (m) + deadcrootc_xfer_to_deadcrootc_p(m) * deltim
      end if
      if(ivt >= npcropmin) then
         AKX_livestemc_xf_to_livestemc_p_acc(m) = AKX_livestemc_xf_to_livestemc_p_acc(m) + livestemc_xfer_to_livestemc_p(m) * deltim
         AKX_livestemc_xf_exit_p_acc        (m) = AKX_livestemc_xf_exit_p_acc        (m) + livestemc_xfer_to_livestemc_p(m) * deltim
         AKX_grainc_xf_to_grainc_p_acc      (m) = AKX_grainc_xf_to_grainc_p_acc      (m) + grainc_xfer_to_grainc_p      (m) * deltim
         AKX_grainc_xf_exit_p_acc           (m) = AKX_grainc_xf_exit_p_acc           (m) + grainc_xfer_to_grainc_p      (m) * deltim
      end if
#endif

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

#ifdef SASU
      AKX_leafc_exit_p_acc (m) = AKX_leafc_exit_p_acc (m) + leafc_to_litter_p (m) * deltim
      AKX_frootc_exit_p_acc(m) = AKX_frootc_exit_p_acc(m) + frootc_to_litter_p(m) * deltim
      if(woody(ivt) == 1) then
         AKX_livestemc_to_deadstemc_p_acc  (m) = AKX_livestemc_to_deadstemc_p_acc  (m) + livestemc_to_deadstemc_p  (m) * deltim
         AKX_livestemc_exit_p_acc          (m) = AKX_livestemc_exit_p_acc          (m) + livestemc_to_deadstemc_p  (m) * deltim
         AKX_livecrootc_to_deadcrootc_p_acc(m) = AKX_livecrootc_to_deadcrootc_p_acc(m) + livecrootc_to_deadcrootc_p(m) * deltim
         AKX_livecrootc_exit_p_acc         (m) = AKX_livecrootc_exit_p_acc         (m) + livecrootc_to_deadcrootc_p(m) * deltim
      end if
      if(ivt >= npcropmin) then
         AKX_livestemc_exit_p_acc          (m) = AKX_livestemc_exit_p_acc          (m) + livestemc_to_litter_p                  (m)  * deltim
         AKX_grainc_exit_p_acc             (m) = AKX_grainc_exit_p_acc             (m) + (grainc_to_food_p(m) + grainc_to_seed_p(m)) * deltim
      end if
#endif
   ! maintenance respiration fluxes from xsmrpool
      cpool_p   (m) = cpool_p   (m) - cpool_to_xsmrpool_p(m) * deltim
      cpool_p   (m) = cpool_p   (m) - leaf_curmr_p       (m) * deltim
      cpool_p   (m) = cpool_p   (m) - froot_curmr_p      (m) * deltim
      if (woody(ivt) == 1) then
         cpool_p(m) = cpool_p   (m) - livestem_curmr_p   (m) * deltim
         cpool_p(m) = cpool_p   (m) - livecroot_curmr_p  (m) * deltim
      end if
      if (ivt >= npcropmin) then
         cpool_p(m) = cpool_p   (m) - livestem_curmr_p   (m) * deltim
         cpool_p(m) = cpool_p   (m) - grain_curmr_p      (m) * deltim
      end if
#ifdef FUN
      cpool_p   (m) = cpool_p   (m) - soilc_change_p     (m) * deltim
#endif
      xsmrpool_p(m) = xsmrpool_p(m) + cpool_to_xsmrpool_p(m) * deltim
      xsmrpool_p(m) = xsmrpool_p(m) - leaf_xsmr_p        (m) * deltim
      xsmrpool_p(m) = xsmrpool_p(m) - froot_xsmr_p       (m) * deltim
      if (woody(ivt) == 1) then
         xsmrpool_p(m) = xsmrpool_p(m) - livestem_xsmr_p (m) * deltim
         xsmrpool_p(m) = xsmrpool_p(m) - livecroot_xsmr_p(m) * deltim
      end if
      cpool_p         (m) = cpool_p         (m) - cpool_to_leafc_p         (m) * deltim
      leafc_p         (m) = leafc_p         (m) + cpool_to_leafc_p         (m) * deltim
      cpool_p         (m) = cpool_p         (m) - cpool_to_leafc_storage_p (m) * deltim
      leafc_storage_p (m) = leafc_storage_p (m) + cpool_to_leafc_storage_p (m) * deltim
      cpool_p         (m) = cpool_p         (m) - cpool_to_frootc_p        (m) * deltim
      frootc_p        (m) = frootc_p        (m) + cpool_to_frootc_p        (m) * deltim
      cpool_p         (m) = cpool_p         (m) - cpool_to_frootc_storage_p(m) * deltim
      frootc_storage_p(m) = frootc_storage_p(m) + cpool_to_frootc_storage_p(m) * deltim
      if (woody(ivt) == 1) then
         cpool_p             (m) = cpool_p             (m) - cpool_to_livestemc_p         (m) * deltim
         livestemc_p         (m) = livestemc_p         (m) + cpool_to_livestemc_p         (m) * deltim
         cpool_p             (m) = cpool_p             (m) - cpool_to_livestemc_storage_p (m) * deltim
         livestemc_storage_p (m) = livestemc_storage_p (m) + cpool_to_livestemc_storage_p (m) * deltim
         cpool_p             (m) = cpool_p             (m) - cpool_to_deadstemc_p         (m) * deltim
         deadstemc_p         (m) = deadstemc_p         (m) + cpool_to_deadstemc_p         (m) * deltim
         cpool_p             (m) = cpool_p             (m) - cpool_to_deadstemc_storage_p (m) * deltim
         deadstemc_storage_p (m) = deadstemc_storage_p (m) + cpool_to_deadstemc_storage_p (m) * deltim
         cpool_p             (m) = cpool_p             (m) - cpool_to_livecrootc_p        (m) * deltim
         livecrootc_p        (m) = livecrootc_p        (m) + cpool_to_livecrootc_p        (m) * deltim
         cpool_p             (m) = cpool_p             (m) - cpool_to_livecrootc_storage_p(m) * deltim
         livecrootc_storage_p(m) = livecrootc_storage_p(m) + cpool_to_livecrootc_storage_p(m) * deltim
         cpool_p             (m) = cpool_p             (m) - cpool_to_deadcrootc_p        (m) * deltim
         deadcrootc_p        (m) = deadcrootc_p        (m) + cpool_to_deadcrootc_p        (m) * deltim
         cpool_p             (m) = cpool_p             (m) - cpool_to_deadcrootc_storage_p(m) * deltim
         deadcrootc_storage_p(m) = deadcrootc_storage_p(m) + cpool_to_deadcrootc_storage_p(m) * deltim
      end if
      if (ivt >= npcropmin) then ! skip 2 generic crops
         cpool_p            (m) = cpool_p            (m) - cpool_to_livestemc_p        (m) * deltim
         livestemc_p        (m) = livestemc_p        (m) + cpool_to_livestemc_p        (m) * deltim
         cpool_p            (m) = cpool_p            (m) - cpool_to_livestemc_storage_p(m) * deltim
         livestemc_storage_p(m) = livestemc_storage_p(m) + cpool_to_livestemc_storage_p(m) * deltim
         cpool_p            (m) = cpool_p            (m) - cpool_to_grainc_p           (m) * deltim
         grainc_p           (m) = grainc_p           (m) + cpool_to_grainc_p           (m) * deltim
         cpool_p            (m) = cpool_p            (m) - cpool_to_grainc_storage_p   (m) * deltim
         grainc_storage_p   (m) = grainc_storage_p   (m) + cpool_to_grainc_storage_p   (m) * deltim
      end if
#ifdef SASU
      I_leafc_p_acc(m) = I_leafc_p_acc(m) + cpool_to_leafc_p         (m) * deltim
      I_leafc_st_p_acc(m) = I_leafc_st_p_acc(m) + cpool_to_leafc_storage_p     (m) * deltim
      I_frootc_p_acc(m) = I_frootc_p_acc(m) + cpool_to_frootc_p         (m) * deltim
      I_frootc_st_p_acc(m) = I_frootc_st_p_acc(m) + cpool_to_frootc_storage_p     (m) * deltim
      if(woody(ivt) == 1) then
         I_livestemc_p_acc    (m) = I_livestemc_p_acc    (m) + cpool_to_livestemc_p         (m) * deltim
         I_livestemc_st_p_acc (m) = I_livestemc_st_p_acc (m) + cpool_to_livestemc_storage_p (m) * deltim
         I_deadstemc_p_acc    (m) = I_deadstemc_p_acc    (m) + cpool_to_deadstemc_p         (m) * deltim
         I_deadstemc_st_p_acc (m) = I_deadstemc_st_p_acc (m) + cpool_to_deadstemc_storage_p (m) * deltim
         I_livecrootc_p_acc   (m) = I_livecrootc_p_acc   (m) + cpool_to_livecrootc_p        (m) * deltim
         I_livecrootc_st_p_acc(m) = I_livecrootc_st_p_acc(m) + cpool_to_livecrootc_storage_p(m) * deltim
         I_deadcrootc_p_acc   (m) = I_deadcrootc_p_acc   (m) + cpool_to_deadcrootc_p        (m) * deltim
         I_deadcrootc_st_p_acc(m) = I_deadcrootc_st_p_acc(m) + cpool_to_deadcrootc_storage_p(m) * deltim
      end if
      if(ivt >= npcropmin) then
         I_livestemc_p_acc    (m) = I_livestemc_p_acc    (m) + cpool_to_livestemc_p         (m) * deltim
         I_livestemc_st_p_acc (m) = I_livestemc_st_p_acc (m) + cpool_to_livestemc_storage_p (m) * deltim
         I_grainc_p_acc       (m) = I_grainc_p_acc       (m) + cpool_to_grainc_p            (m) * deltim
         I_grainc_st_p_acc    (m) = I_grainc_st_p_acc    (m) + cpool_to_grainc_storage_p    (m) * deltim
      end if
#endif
      ! growth respiration for transfer growth
      cpool_p     (m) = cpool_p     (m) - cpool_leaf_gr_p     (m) * deltim
      cpool_p     (m) = cpool_p     (m) - cpool_froot_gr_p    (m) * deltim
      if(woody(ivt) == 1) then
         cpool_p  (m) = cpool_p     (m) - cpool_livestem_gr_p (m) * deltim
         cpool_p  (m) = cpool_p     (m) - cpool_deadstem_gr_p (m) * deltim
         cpool_p  (m) = cpool_p     (m) - cpool_livecroot_gr_p(m) * deltim
         cpool_p  (m) = cpool_p     (m) - cpool_deadcroot_gr_p(m) * deltim
      end if
      if(ivt >= npcropmin)then
         cpool_p  (m) = cpool_p     (m) - cpool_livestem_gr_p (m) * deltim
         cpool_p  (m) = cpool_p     (m) - cpool_grain_gr_p    (m) * deltim
      end if

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
   ! growth respiration at time of storage
      cpool_p        (m) = cpool_p        (m) - cpool_leaf_storage_gr_p (m) * deltim
      cpool_p        (m) = cpool_p        (m) - cpool_froot_storage_gr_p(m) * deltim
      if(woody(ivt) == 1)then
          cpool_p    (m) = cpool_p        (m) - cpool_livestem_storage_gr_p (m) * deltim
          cpool_p    (m) = cpool_p        (m) - cpool_deadstem_storage_gr_p (m) * deltim
          cpool_p    (m) = cpool_p        (m) - cpool_livecroot_storage_gr_p(m) * deltim
          cpool_p    (m) = cpool_p        (m) - cpool_deadcroot_storage_gr_p(m) * deltim
      end if
      if(ivt >= npcropmin)then
          cpool_p    (m) = cpool_p        (m) - cpool_livestem_storage_gr_p (m) * deltim
          cpool_p    (m) = cpool_p        (m) - cpool_grain_storage_gr_p    (m) * deltim
      end if
          
   ! growth respiration stored for release during transfer growth
      cpool_p        (m) = cpool_p        (m) - cpool_to_gresp_storage_p(m) * deltim
      gresp_storage_p(m) = gresp_storage_p(m) + cpool_to_gresp_storage_p(m) * deltim

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
#ifdef SASU
      AKX_leafc_st_to_leafc_xf_p_acc  (m) = AKX_leafc_st_to_leafc_xf_p_acc  (m) + leafc_storage_to_xfer_p (m) * deltim
      AKX_leafc_st_exit_p_acc         (m) = AKX_leafc_st_exit_p_acc         (m) + leafc_storage_to_xfer_p (m) * deltim
      AKX_frootc_st_to_frootc_xf_p_acc(m) = AKX_frootc_st_to_frootc_xf_p_acc(m) + frootc_storage_to_xfer_p(m) * deltim
      AKX_frootc_st_exit_p_acc        (m) = AKX_frootc_st_exit_p_acc        (m) + frootc_storage_to_xfer_p(m) * deltim
      if(woody(ivt) == 1) then
         AKX_livestemc_st_to_livestemc_xf_p_acc  (m) = AKX_livestemc_st_to_livestemc_xf_p_acc  (m) + livestemc_storage_to_xfer_p (m) * deltim
         AKX_livestemc_st_exit_p_acc             (m) = AKX_livestemc_st_exit_p_acc             (m) + livestemc_storage_to_xfer_p (m) * deltim
         AKX_deadstemc_st_to_deadstemc_xf_p_acc  (m) = AKX_deadstemc_st_to_deadstemc_xf_p_acc  (m) + deadstemc_storage_to_xfer_p (m) * deltim
         AKX_deadstemc_st_exit_p_acc             (m) = AKX_deadstemc_st_exit_p_acc             (m) + deadstemc_storage_to_xfer_p (m) * deltim
         AKX_livecrootc_st_to_livecrootc_xf_p_acc(m) = AKX_livecrootc_st_to_livecrootc_xf_p_acc(m) + livecrootc_storage_to_xfer_p(m) * deltim
         AKX_livecrootc_st_exit_p_acc            (m) = AKX_livecrootc_st_exit_p_acc            (m) + livecrootc_storage_to_xfer_p(m) * deltim
         AKX_deadcrootc_st_to_deadcrootc_xf_p_acc(m) = AKX_deadcrootc_st_to_deadcrootc_xf_p_acc(m) + deadcrootc_storage_to_xfer_p(m) * deltim
         AKX_deadcrootc_st_exit_p_acc            (m) = AKX_deadcrootc_st_exit_p_acc            (m) + deadcrootc_storage_to_xfer_p(m) * deltim
      end if
      if( ivt >= npcropmin) then
         AKX_livestemc_st_to_livestemc_xf_p_acc  (m) = AKX_livestemc_st_to_livestemc_xf_p_acc  (m) + livestemc_storage_to_xfer_p (m) * deltim
         AKX_livestemc_st_exit_p_acc             (m) = AKX_livestemc_st_exit_p_acc             (m) + livestemc_storage_to_xfer_p (m) * deltim
         AKX_grainc_st_to_grainc_xf_p_acc        (m) = AKX_grainc_st_to_grainc_xf_p_acc        (m) + grainc_storage_to_xfer_p    (m) * deltim
         AKX_grainc_st_exit_p_acc                (m) = AKX_grainc_st_exit_p_acc                (m) + grainc_storage_to_xfer_p    (m) * deltim
      end if
#endif
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
            xsmrpool_to_atm_p(m) = xsmrpool_to_atm_p(m) + cpool_p (m)/deltim
            cpool_p          (m) = 0._r8
            xsmrpool_to_atm_p(m) = xsmrpool_to_atm_p(m) + frootc_p(m)/deltim
            frootc_p         (m) = 0._r8
         end if
      end if
#ifdef CROP
      cropprod1c_loss_p(m) = cropprod1c_p(m) * 7.2e-8_r8
      cropprod1c_p     (m) = cropprod1c_p(m) + grainc_to_food_p(m) * deltim - cropprod1c_loss_p(m) * deltim
#endif
   end do ! end pft loop

end subroutine CStateUpdate1

end module bgc_CNCStateUpdate1Mod
