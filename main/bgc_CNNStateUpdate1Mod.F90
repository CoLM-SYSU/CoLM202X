#include <define.h>
module bgc_CNNStateUpdate1Mod

use precision
use MOD_PFTimeInvars, only: pftclass
use PFT_Const, only: woody
use MOD_TimeInvariants, only: &
  ! bgc constants
           donor_pool, receiver_pool, i_met_lit, i_cel_lit, i_lig_lit, i_cwd, i_soil1, i_soil2, i_soil3

use MOD_TimeVariables, only: &
           I_met_n_vr_acc, I_cel_n_vr_acc, I_lig_n_vr_acc

use MOD_1D_Fluxes, only: &
! decomposition pools flux varables (in)
           decomp_npools_sourcesink, &
           phenology_to_met_n      , phenology_to_cel_n,  phenology_to_lig_n

use MOD_PFTimeVars, only: &
! vegetation carbon state variables (inout)
           leafn_p            , leafn_storage_p     , leafn_xfer_p     , &
           frootn_p           , frootn_storage_p    , frootn_xfer_p    , &
           livestemn_p        , livestemn_storage_p , livestemn_xfer_p , &
           deadstemn_p        , deadstemn_storage_p , deadstemn_xfer_p , &
           livecrootn_p       , livecrootn_storage_p, livecrootn_xfer_p, &
           deadcrootn_p       , deadcrootn_storage_p, deadcrootn_xfer_p, &
           grainn_p           , grainn_storage_p    , grainn_xfer_p    , &
           cropseedn_deficit_p, retransn_p          , npool_p          , &
     
! SASU variables
           I_leafn_p_acc     , I_leafn_st_p_acc     , I_frootn_p_acc    , I_frootn_st_p_acc    , &
           I_livestemn_p_acc , I_livestemn_st_p_acc , I_deadstemn_p_acc , I_deadstemn_st_p_acc , &
           I_livecrootn_p_acc, I_livecrootn_st_p_acc, I_deadcrootn_p_acc, I_deadcrootn_st_p_acc, &
           I_grainn_p_acc    , I_grainn_st_p_acc    , &

           AKX_leafn_xf_to_leafn_p_acc           , AKX_frootn_xf_to_frootn_p_acc           , AKX_livestemn_xf_to_livestemn_p_acc     , &
           AKX_deadstemn_xf_to_deadstemn_p_acc   , AKX_livecrootn_xf_to_livecrootn_p_acc   , AKX_deadcrootn_xf_to_deadcrootn_p_acc   , &
           AKX_grainn_xf_to_grainn_p_acc         , AKX_livestemn_to_deadstemn_p_acc        , AKX_livecrootn_to_deadcrootn_p_acc      , &

           AKX_leafn_st_to_leafn_xf_p_acc        , AKX_frootn_st_to_frootn_xf_p_acc        , AKX_livestemn_st_to_livestemn_xf_p_acc  , &
           AKX_deadstemn_st_to_deadstemn_xf_p_acc, AKX_livecrootn_st_to_livecrootn_xf_p_acc, AKX_deadcrootn_st_to_deadcrootn_xf_p_acc, &
           AKX_livestemn_st_to_livestemn_xf_p_acc, AKX_grainn_st_to_grainn_xf_p_acc        , &

           AKX_leafn_to_retransn_p_acc           , AKX_frootn_to_retransn_p_acc            , AKX_livestemn_to_retransn_p_acc         , &
           AKX_livecrootn_to_retransn_p_acc      , &
 
           AKX_retransn_to_leafn_p_acc           , AKX_retransn_to_frootn_p_acc            , AKX_retransn_to_livestemn_p_acc         , &
           AKX_retransn_to_deadstemn_p_acc       , AKX_retransn_to_livecrootn_p_acc        , AKX_retransn_to_deadcrootn_p_acc        , &
           AKX_retransn_to_grainn_p_acc          , &

           AKX_retransn_to_leafn_st_p_acc        , AKX_retransn_to_frootn_st_p_acc         , AKX_retransn_to_livestemn_st_p_acc      , &
           AKX_retransn_to_deadstemn_st_p_acc    , AKX_retransn_to_livecrootn_st_p_acc     , AKX_retransn_to_deadcrootn_st_p_acc     , &
           AKX_retransn_to_grainn_st_p_acc       , &

           AKX_leafn_exit_p_acc                  , AKX_frootn_exit_p_acc                   , AKX_livestemn_exit_p_acc                , &
           AKX_deadstemn_exit_p_acc              , AKX_livecrootn_exit_p_acc               , AKX_deadcrootn_exit_p_acc               , &
           AKX_grainn_exit_p_acc                 , AKX_retransn_exit_p_acc                 , &

           AKX_leafn_st_exit_p_acc               , AKX_frootn_st_exit_p_acc                , AKX_livestemn_st_exit_p_acc             , &
           AKX_deadstemn_st_exit_p_acc           , AKX_livecrootn_st_exit_p_acc            , AKX_deadcrootn_st_exit_p_acc            , &
           AKX_grainn_st_exit_p_acc              , & 

           AKX_leafn_xf_exit_p_acc               , AKX_frootn_xf_exit_p_acc                , AKX_livestemn_xf_exit_p_acc             , &
           AKX_deadstemn_xf_exit_p_acc           , AKX_livecrootn_xf_exit_p_acc            , AKX_deadcrootn_xf_exit_p_acc            , &
           AKX_grainn_xf_exit_p_acc               

use MOD_1D_PFTFluxes, only: &
! vegetation carbon flux variables (in)
! xfer to display
           leafn_xfer_to_leafn_p          , frootn_xfer_to_frootn_p        , &
           livestemn_xfer_to_livestemn_p  , deadstemn_xfer_to_deadstemn_p  , &
           livecrootn_xfer_to_livecrootn_p, deadcrootn_xfer_to_deadcrootn_p, &
           grainn_xfer_to_grainn_p        , &

! storage to xfer (in)
           leafn_storage_to_xfer_p     , frootn_storage_to_xfer_p    , &
           livestemn_storage_to_xfer_p , deadstemn_storage_to_xfer_p , &
           livecrootn_storage_to_xfer_p, deadcrootn_storage_to_xfer_p, &
           grainn_storage_to_xfer_p    , &

! display to litter & live to dead (in)
           leafn_to_litter_p       , frootn_to_litter_p        , &
           grainn_to_food_p        , grainn_to_seed_p          , &
           crop_seedn_to_leaf_p    , livestemn_to_litter_p     , &
           livestemn_to_deadstemn_p, livecrootn_to_deadcrootn_p, &

! display to retransn / retransn to npool (in)
           leafn_to_retransn_p     , frootn_to_retransn_p      , &
           livestemn_to_retransn_p , livecrootn_to_retransn_p  , &
           retransn_to_npool_p     , free_retransn_to_npool_p  , &

! npool to display/storage (in)
           npool_to_leafn_p     , npool_to_leafn_storage_p     , &
           npool_to_frootn_p    , npool_to_frootn_storage_p    , &
           npool_to_livestemn_p , npool_to_livestemn_storage_p , &
           npool_to_deadstemn_p , npool_to_deadstemn_storage_p , &
           npool_to_livecrootn_p, npool_to_livecrootn_storage_p, &
           npool_to_deadcrootn_p, npool_to_deadcrootn_storage_p, &
           npool_to_grainn_p    , npool_to_grainn_storage_p    , plant_nalloc_p 

use MOD_PFTimeInvars, only: pftfrac
implicit none

public NStateUpdate1

contains

subroutine NStateUpdate1 &
! model configurations
          (i, ps, pe, deltim, nl_soil, ndecomp_transitions, npcropmin,dz_soi)

integer ,intent(in) :: i
integer ,intent(in) :: ps
integer ,intent(in) :: pe
real(r8),intent(in) :: deltim
integer ,intent(in) :: nl_soil
integer ,intent(in) :: ndecomp_transitions
integer ,intent(in) :: npcropmin
real(r8),intent(in) :: dz_soi(1:nl_soil)

integer j,k
integer ivt, m
real(r8) f_retr_in_nall
!if(i .eq. 123226)print*,'nsource_sink before NStateUpdate1',&
!  sum(decomp_npools_sourcesink(1:nl_soil,i_met_lit,i)*dz_soi(1:nl_soil)) &
! +sum(decomp_npools_sourcesink(1:nl_soil,i_cel_lit,i)*dz_soi(1:nl_soil)) &
! +sum(decomp_npools_sourcesink(1:nl_soil,i_cwd,i)*dz_soi(1:nl_soil)) &
! +sum(decomp_npools_sourcesink(1:nl_soil,i_soil1,i)*dz_soi(1:nl_soil)) &
! +sum(decomp_npools_sourcesink(1:nl_soil,i_soil2,i)*dz_soi(1:nl_soil)) &
! +sum(decomp_npools_sourcesink(1:nl_soil,i_soil3,i)*dz_soi(1:nl_soil))

!if(i .eq. 123226)print*,'vegn before NStateUpdate1',&
!             sum((leafn_p(ps:pe)             + frootn_p(ps:pe)             + livestemn_p(ps:pe) &
!                + deadstemn_p(ps:pe)         + livecrootn_p(ps:pe)         + deadcrootn_p(ps:pe) &
!                + leafn_storage_p(ps:pe)     + frootn_storage_p(ps:pe)     + livestemn_storage_p(ps:pe) &
!                + deadstemn_storage_p(ps:pe) + livecrootn_storage_p(ps:pe) + deadcrootn_storage_p(ps:pe) &
!                + leafn_xfer_p(ps:pe)        + frootn_xfer_p(ps:pe)        + livestemn_xfer_p(ps:pe) &
!                + deadstemn_xfer_p(ps:pe)    + livecrootn_xfer_p(ps:pe)    + deadcrootn_xfer_p(ps:pe) &
!                + npool_p(ps:pe)             + retransn_p(ps:pe) &
!                + grainn_p(ps:pe)            + grainn_storage_p(ps:pe)     + grainn_xfer_p(ps:pe) &
!                + cropseedn_deficit_p(ps:pe))*pftfrac(ps:pe))

      ! soilbiogeochemistry fluxes TODO - this should be moved elsewhere
      ! plant to litter fluxes -  phenology and dynamic landcover fluxes
      do j = 1, nl_soil
         decomp_npools_sourcesink(j,i_met_lit,i) = phenology_to_met_n(j,i) * deltim

         decomp_npools_sourcesink(j,i_cel_lit,i) = phenology_to_cel_n(j,i) * deltim

         decomp_npools_sourcesink(j,i_lig_lit,i) = phenology_to_lig_n(j,i) * deltim

            ! NOTE(wjs, 2017-01-02) This used to be set to a non-zero value, but the
            ! terms have been moved to CStateUpdateDynPatch. I think this is zeroed every
            ! time step, but to be safe, I'm explicitly setting it to zero here.
         decomp_npools_sourcesink(j,i_cwd,i) = 0._r8

      end do
!      if(i .eq. 123226)print*,'phen_to_lit',deltim*sum(phenology_to_met_n(1:nl_soil,i)*dz_soi(1:nl_soil)) &
!                    + deltim*sum(phenology_to_cel_n(1:nl_soil,i)*dz_soi(1:nl_soil)) &
!                    + deltim*sum(phenology_to_lig_n(1:nl_soil,i)*dz_soi(1:nl_soil)) 

#ifdef SASU
   do j=1,nl_soil
      I_met_n_vr_acc(j,i) = I_met_n_vr_acc(j,i) + phenology_to_met_n(j,i) * deltim
      I_cel_n_vr_acc(j,i) = I_cel_n_vr_acc(j,i) + phenology_to_cel_n(j,i) * deltim
      I_lig_n_vr_acc(j,i) = I_lig_n_vr_acc(j,i) + phenology_to_lig_n(j,i) * deltim
   end do
#endif

      do m = ps , pe
         ivt = pftclass(m)
      ! phenology: transfer growth fluxes
         leafn_p(m)       = leafn_p(m)       + leafn_xfer_to_leafn_p(m)*deltim
         leafn_xfer_p(m)  = leafn_xfer_p(m)  - leafn_xfer_to_leafn_p(m)*deltim
         frootn_p(m)      = frootn_p(m)      + frootn_xfer_to_frootn_p(m)*deltim
         frootn_xfer_p(m) = frootn_xfer_p(m) - frootn_xfer_to_frootn_p(m)*deltim

         if (woody(ivt) == 1) then
            livestemn_p(m)       = livestemn_p(m)       + livestemn_xfer_to_livestemn_p(m)*deltim
            livestemn_xfer_p(m)  = livestemn_xfer_p(m)  - livestemn_xfer_to_livestemn_p(m)*deltim
            deadstemn_p(m)       = deadstemn_p(m)       + deadstemn_xfer_to_deadstemn_p(m)*deltim
            deadstemn_xfer_p(m)  = deadstemn_xfer_p(m)  - deadstemn_xfer_to_deadstemn_p(m)*deltim
            livecrootn_p(m)      = livecrootn_p(m)      + livecrootn_xfer_to_livecrootn_p(m)*deltim
            livecrootn_xfer_p(m) = livecrootn_xfer_p(m) - livecrootn_xfer_to_livecrootn_p(m)*deltim
            deadcrootn_p(m)      = deadcrootn_p(m)      + deadcrootn_xfer_to_deadcrootn_p(m)*deltim
            deadcrootn_xfer_p(m) = deadcrootn_xfer_p(m) - deadcrootn_xfer_to_deadcrootn_p(m)*deltim
         end if

         if (ivt >= npcropmin) then ! skip 2 generic crops
            ! lines here for consistency; the transfer terms are zero
            livestemn_p(m)       = livestemn_p(m)      + livestemn_xfer_to_livestemn_p(m)*deltim
            livestemn_xfer_p(m)  = livestemn_xfer_p(m) - livestemn_xfer_to_livestemn_p(m)*deltim
            grainn_p(m)          = grainn_p(m)         + grainn_xfer_to_grainn_p(m)*deltim
            grainn_xfer_p(m)     = grainn_xfer_p(m)    - grainn_xfer_to_grainn_p(m)*deltim
         end if
#ifdef SASU
         AKX_leafn_xf_to_leafn_p_acc  (m) = AKX_leafn_xf_to_leafn_p_acc  (m) + leafn_xfer_to_leafn_p  (m) * deltim
         AKX_frootn_xf_to_frootn_p_acc(m) = AKX_frootn_xf_to_frootn_p_acc(m) + frootn_xfer_to_frootn_p(m) * deltim
         AKX_leafn_xf_exit_p_acc      (m) = AKX_leafn_xf_exit_p_acc      (m) + leafn_xfer_to_leafn_p  (m) * deltim
         AKX_frootn_xf_exit_p_acc     (m) = AKX_frootn_xf_exit_p_acc     (m) + frootn_xfer_to_frootn_p(m) * deltim
         if(woody(ivt) == 1)then
            AKX_livestemn_xf_to_livestemn_p_acc  (m) = AKX_livestemn_xf_to_livestemn_p_acc  (m) + livestemn_xfer_to_livestemn_p  (m) * deltim
            AKX_livestemn_xf_exit_p_acc          (m) = AKX_livestemn_xf_exit_p_acc          (m) + livestemn_xfer_to_livestemn_p  (m) * deltim
            AKX_deadstemn_xf_to_deadstemn_p_acc  (m) = AKX_deadstemn_xf_to_deadstemn_p_acc  (m) + deadstemn_xfer_to_deadstemn_p  (m) * deltim
            AKX_deadstemn_xf_exit_p_acc          (m) = AKX_deadstemn_xf_exit_p_acc          (m) + deadstemn_xfer_to_deadstemn_p  (m) * deltim
            AKX_livecrootn_xf_to_livecrootn_p_acc(m) = AKX_livecrootn_xf_to_livecrootn_p_acc(m) + livecrootn_xfer_to_livecrootn_p(m) * deltim
            AKX_livecrootn_xf_exit_p_acc         (m) = AKX_livecrootn_xf_exit_p_acc         (m) + livecrootn_xfer_to_livecrootn_p(m) * deltim
            AKX_deadcrootn_xf_to_deadcrootn_p_acc(m) = AKX_deadcrootn_xf_to_deadcrootn_p_acc(m) + deadcrootn_xfer_to_deadcrootn_p(m) * deltim
            AKX_deadcrootn_xf_exit_p_acc         (m) = AKX_deadcrootn_xf_exit_p_acc         (m) + deadcrootn_xfer_to_deadcrootn_p(m) * deltim
         end if
         if(ivt >= npcropmin) then
            AKX_livestemn_xf_to_livestemn_p_acc(m) = AKX_livestemn_xf_to_livestemn_p_acc(m) + livestemn_xfer_to_livestemn_p(m) * deltim
            AKX_livestemn_xf_exit_p_acc        (m) = AKX_livestemn_xf_exit_p_acc        (m) + livestemn_xfer_to_livestemn_p(m) * deltim
            AKX_grainn_xf_to_grainn_p_acc      (m) = AKX_grainn_xf_to_grainn_p_acc      (m) + grainn_xfer_to_grainn_p      (m) * deltim
            AKX_grainn_xf_exit_p_acc           (m) = AKX_grainn_xf_exit_p_acc           (m) + grainn_xfer_to_grainn_p      (m) * deltim
         end if
#endif

!      end do
!    
!if(i .eq. 123226)print*,'vegn before phen litfall',&
!             sum((leafn_p(ps:pe)             + frootn_p(ps:pe)             + livestemn_p(ps:pe) &
!                + deadstemn_p(ps:pe)         + livecrootn_p(ps:pe)         + deadcrootn_p(ps:pe) &
!                + leafn_storage_p(ps:pe)     + frootn_storage_p(ps:pe)     + livestemn_storage_p(ps:pe) &
!                + deadstemn_storage_p(ps:pe) + livecrootn_storage_p(ps:pe) + deadcrootn_storage_p(ps:pe) &
!                + leafn_xfer_p(ps:pe)        + frootn_xfer_p(ps:pe)        + livestemn_xfer_p(ps:pe) &
!                + deadstemn_xfer_p(ps:pe)    + livecrootn_xfer_p(ps:pe)    + deadcrootn_xfer_p(ps:pe) &
!                + npool_p(ps:pe)             + retransn_p(ps:pe) &
!                + grainn_p(ps:pe)            + grainn_storage_p(ps:pe)     + grainn_xfer_p(ps:pe) &
!                + cropseedn_deficit_p(ps:pe))*pftfrac(ps:pe))
!
         ! phenology: litterfall and retranslocation fluxes
         leafn_p(m)    = leafn_p(m)    - leafn_to_litter_p(m)*deltim
         frootn_p(m)   = frootn_p(m)   - frootn_to_litter_p(m)*deltim
         leafn_p(m)    = leafn_p(m)    - leafn_to_retransn_p(m)*deltim
         retransn_p(m) = retransn_p(m) + leafn_to_retransn_p(m)*deltim

         ! live wood turnover and retranslocation fluxes
         if (woody(ivt) == 1) then
            livestemn_p(m)    = livestemn_p(m)  - livestemn_to_deadstemn_p(m)*deltim
            deadstemn_p(m)    = deadstemn_p(m)  + livestemn_to_deadstemn_p(m)*deltim
            livecrootn_p(m)   = livecrootn_p(m) - livecrootn_to_deadcrootn_p(m)*deltim
            deadcrootn_p(m)   = deadcrootn_p(m) + livecrootn_to_deadcrootn_p(m)*deltim

            livestemn_p(m)    = livestemn_p(m)  - livestemn_to_retransn_p(m)*deltim
            retransn_p(m)     = retransn_p(m)   + livestemn_to_retransn_p(m)*deltim
            livecrootn_p(m)   = livecrootn_p(m) - livecrootn_to_retransn_p(m)*deltim
            retransn_p(m)     = retransn_p(m)   + livecrootn_to_retransn_p(m)*deltim
         end if 
         if (ivt >= npcropmin) then ! Beth adds retrans from froot
            frootn_p(m)       = frootn_p(m)     - frootn_to_retransn_p(m)*deltim
            retransn_p(m)     = retransn_p(m)   + frootn_to_retransn_p(m)*deltim
            livestemn_p(m)    = livestemn_p(m)  - livestemn_to_litter_p(m)*deltim
            livestemn_p(m)    = livestemn_p(m)  - livestemn_to_retransn_p(m)*deltim
            retransn_p(m)     = retransn_p(m)   + livestemn_to_retransn_p(m)*deltim
            grainn_p(m)       = grainn_p(m) &
                   - (grainn_to_food_p(m) + grainn_to_seed_p(m))*deltim
            cropseedn_deficit_p(m) = cropseedn_deficit_p(m) &
                 - crop_seedn_to_leaf_p(m) * deltim &
                 + grainn_to_seed_p(m) * deltim
         end if
!      end do
#ifdef SASU
         AKX_leafn_exit_p_acc       (m) = AKX_leafn_exit_p_acc       (m) + leafn_to_litter_p  (m) * deltim
         AKX_frootn_exit_p_acc      (m) = AKX_frootn_exit_p_acc      (m) + frootn_to_litter_p (m) * deltim
         AKX_leafn_to_retransn_p_acc(m) = AKX_leafn_to_retransn_p_acc(m) + leafn_to_retransn_p(m) * deltim
         AKX_leafn_exit_p_acc       (m) = AKX_leafn_exit_p_acc       (m) + leafn_to_retransn_p(m) * deltim
         if(woody(ivt) == 1) then
            AKX_livestemn_to_deadstemn_p_acc  (m) = AKX_livestemn_to_deadstemn_p_acc  (m) + livestemn_to_deadstemn_p  (m) * deltim
            AKX_livestemn_exit_p_acc          (m) = AKX_livestemn_exit_p_acc          (m) + livestemn_to_deadstemn_p  (m) * deltim
            AKX_livecrootn_to_deadcrootn_p_acc(m) = AKX_livecrootn_to_deadcrootn_p_acc(m) + livecrootn_to_deadcrootn_p(m) * deltim
            AKX_livecrootn_exit_p_acc         (m) = AKX_livecrootn_exit_p_acc         (m) + livecrootn_to_deadcrootn_p(m) * deltim
 
            AKX_livestemn_to_retransn_p_acc   (m) = AKX_livestemn_to_retransn_p_acc   (m) + livestemn_to_retransn_p   (m) * deltim
            AKX_livestemn_exit_p_acc          (m) = AKX_livestemn_exit_p_acc          (m) + livestemn_to_retransn_p   (m) * deltim
            AKX_livecrootn_to_retransn_p_acc  (m) = AKX_livecrootn_to_retransn_p_acc  (m) + livecrootn_to_retransn_p  (m) * deltim
            AKX_livecrootn_exit_p_acc         (m) = AKX_livecrootn_exit_p_acc         (m) + livecrootn_to_retransn_p  (m) * deltim
         end if
         if(ivt >= npcropmin) then
            AKX_frootn_to_retransn_p_acc      (m) = AKX_frootn_to_retransn_p_acc      (m) + frootn_to_retransn_p      (m) * deltim
            AKX_frootn_exit_p_acc             (m) = AKX_frootn_exit_p_acc             (m) + frootn_to_retransn_p      (m) * deltim
            AKX_livestemn_exit_p_acc          (m) = AKX_livestemn_exit_p_acc          (m) + livestemn_to_litter_p     (m) * deltim
            AKX_livestemn_to_retransn_p_acc   (m) = AKX_livestemn_to_retransn_p_acc   (m) + livestemn_to_retransn_p   (m) * deltim
            AKX_livestemn_exit_p_acc          (m) = AKX_livestemn_exit_p_acc          (m) + livestemn_to_retransn_p   (m) * deltim
            AKX_grainn_exit_p_acc             (m) = AKX_grainn_exit_p_acc             (m) + (grainn_to_food_p(m) + grainn_to_seed_p(m)) * deltim
         end if
#endif
  
!if(i .eq. 123226)print*,'vegn after litter falll',ivt,npcropmin,&
!             sum((leafn_p(ps:pe)             + frootn_p(ps:pe)             + livestemn_p(ps:pe) &
!                + deadstemn_p(ps:pe)         + livecrootn_p(ps:pe)         + deadcrootn_p(ps:pe) &
!                + leafn_storage_p(ps:pe)     + frootn_storage_p(ps:pe)     + livestemn_storage_p(ps:pe) &
!                + deadstemn_storage_p(ps:pe) + livecrootn_storage_p(ps:pe) + deadcrootn_storage_p(ps:pe) &
!                + leafn_xfer_p(ps:pe)        + frootn_xfer_p(ps:pe)        + livestemn_xfer_p(ps:pe) &
!                + deadstemn_xfer_p(ps:pe)    + livecrootn_xfer_p(ps:pe)    + deadcrootn_xfer_p(ps:pe) &
!                + npool_p(ps:pe)             + retransn_p(ps:pe) &
!                + grainn_p(ps:pe)            + grainn_storage_p(ps:pe)     + grainn_xfer_p(ps:pe) &
!                + cropseedn_deficit_p(ps:pe)&
!                 )*pftfrac(ps:pe))
!if(i .eq. 123226)print*,'vegn before allocation',sum(npool_p(ps:pe)*pftfrac(ps:pe)),&
!             sum((leafn_p(ps:pe)              + leafn_storage_p(ps:pe))*pftfrac(ps:pe)),&
!             sum((frootn_p(ps:pe)             + frootn_storage_p(ps:pe))*pftfrac(ps:pe)),&
!             sum((livestemn_p(ps:pe)          + livestemn_storage_p(ps:pe))*pftfrac(ps:pe)),&
!             sum((deadstemn_p(ps:pe)          + deadstemn_storage_p(ps:pe))*pftfrac(ps:pe)),&
!             sum((livecrootn_p(ps:pe)         + livecrootn_storage_p(ps:pe))*pftfrac(ps:pe)),&
!             sum((deadcrootn_p(ps:pe)         + deadcrootn_storage_p(ps:pe))*pftfrac(ps:pe))
!if(i .eq. 123226)print*,'vegn before allocation,total',&
!             sum((leafn_p(ps:pe)              + leafn_storage_p(ps:pe))*pftfrac(ps:pe))+&
!             sum((frootn_p(ps:pe)             + frootn_storage_p(ps:pe))*pftfrac(ps:pe))+&
!             sum((livestemn_p(ps:pe)          + livestemn_storage_p(ps:pe))*pftfrac(ps:pe))+&
!             sum((deadstemn_p(ps:pe)          + deadstemn_storage_p(ps:pe))*pftfrac(ps:pe))+&
!             sum((livecrootn_p(ps:pe)         + livecrootn_storage_p(ps:pe))*pftfrac(ps:pe))+&
!             sum((deadcrootn_p(ps:pe)         + deadcrootn_storage_p(ps:pe))*pftfrac(ps:pe))
!if(i .eq. 123226)print*,'pftfrac',pftfrac(ps:pe)         
!if(i .eq. 123226)print*,'retransn to npool',sum((retransn_to_npool_p(ps:pe)+free_retransn_to_npool_p(ps:pe))*pftfrac(ps:pe)*deltim)
!      do m = ps , pe
!         ivt=pftclass(m)
!         if(i .eq. 123226)print*,'woody',pftclass(m),woody(pftclass(m))
         ! allocation fluxes
         retransn_p(m)        = retransn_p(m)       - retransn_to_npool_p(m)*deltim
         retransn_p(m)        = retransn_p(m)       - free_retransn_to_npool_p(m)*deltim !how is retransn a state? 
         leafn_p(m)           = leafn_p(m)          + npool_to_leafn_p(m)*deltim
         leafn_storage_p(m)   = leafn_storage_p(m)  + npool_to_leafn_storage_p(m)*deltim
         frootn_p(m)          = frootn_p(m)         + npool_to_frootn_p(m)*deltim
         frootn_storage_p(m)  = frootn_storage_p(m) + npool_to_frootn_storage_p(m)*deltim
         
!         if( i .eq. 123226)print*,'livestemn before npool',m,pftclass(m),livestemn_p(m)+livestemn_storage_p(m),&
!                                   (npool_to_livestemn_p(m)+npool_to_livestemn_storage_p(m))*deltim
         if (woody(ivt) == 1) then
            livestemn_p(m)          = livestemn_p(m)          + npool_to_livestemn_p(m)*deltim
            livestemn_storage_p(m)  = livestemn_storage_p(m)  + npool_to_livestemn_storage_p(m)*deltim
            deadstemn_p(m)          = deadstemn_p(m)          + npool_to_deadstemn_p(m)*deltim
            deadstemn_storage_p(m)  = deadstemn_storage_p(m)  + npool_to_deadstemn_storage_p(m)*deltim
            livecrootn_p(m)         = livecrootn_p(m)         + npool_to_livecrootn_p(m)*deltim
            livecrootn_storage_p(m) = livecrootn_storage_p(m) + npool_to_livecrootn_storage_p(m)*deltim
            deadcrootn_p(m)         = deadcrootn_p(m)         + npool_to_deadcrootn_p(m)*deltim
            deadcrootn_storage_p(m) = deadcrootn_storage_p(m) + npool_to_deadcrootn_storage_p(m)*deltim
         end if
!         if( i .eq. 123226)print*,'livestemn after npool',livestemn_p(m)+livestemn_storage_p(m)

         if (ivt >= npcropmin) then ! skip 2 generic crops
!            if(i .eq. 123226)print*,'ivt>npcropmin',ivt,npcropmin
            livestemn_p(m)          = livestemn_p(m)          + npool_to_livestemn_p(m)*deltim
            livestemn_storage_p(m)  = livestemn_storage_p(m)  + npool_to_livestemn_storage_p(m)*deltim
            grainn_p(m)             = grainn_p(m)             + npool_to_grainn_p(m)*deltim
            grainn_storage_p(m)     = grainn_storage_p(m)     + npool_to_grainn_storage_p(m)*deltim
         end if
!      end do
#ifdef SASU
         if(plant_nalloc_p(m) .ne. 0)then
            f_retr_in_nall = retransn_to_npool_p(m) / plant_nalloc_p(m)
            AKX_retransn_exit_p_acc        (m) = AKX_retransn_exit_p_acc        (m) &
                                               + (retransn_to_npool_p       (m) + free_retransn_to_npool_p (m)) * deltim
            I_leafn_p_acc                  (m) = I_leafn_p_acc              (m) + npool_to_leafn_p         (m) * (1._r8 - f_retr_in_nall) * deltim
            AKX_retransn_to_leafn_p_acc    (m) = AKX_retransn_to_leafn_p_acc    (m) + npool_to_leafn_p         (m) * f_retr_in_nall           * deltim
            I_leafn_st_p_acc               (m) = I_leafn_st_p_acc           (m) + npool_to_leafn_storage_p (m) * (1._r8 - f_retr_in_nall) * deltim
            AKX_retransn_to_leafn_st_p_acc (m) = AKX_retransn_to_leafn_st_p_acc (m) + npool_to_leafn_storage_p (m) * f_retr_in_nall           * deltim
            I_frootn_p_acc                 (m) = I_frootn_p_acc             (m) + npool_to_frootn_p        (m) * (1._r8 - f_retr_in_nall) * deltim
            AKX_retransn_to_frootn_p_acc   (m) = AKX_retransn_to_frootn_p_acc   (m) + npool_to_frootn_p        (m) * f_retr_in_nall           * deltim
            I_frootn_st_p_acc              (m) = I_frootn_st_p_acc          (m) + npool_to_frootn_storage_p(m) * (1._r8 - f_retr_in_nall) * deltim
            AKX_retransn_to_frootn_st_p_acc(m) = AKX_retransn_to_frootn_st_p_acc(m) + npool_to_frootn_storage_p(m) * f_retr_in_nall           * deltim
            if(woody(ivt) == 1)then
               I_livestemn_p_acc                  (m) = I_livestemn_p_acc              (m) &
                                                      + npool_to_livestemn_p           (m) * (1._r8 - f_retr_in_nall) * deltim
               AKX_retransn_to_livestemn_p_acc    (m) = AKX_retransn_to_livestemn_p_acc    (m) &
                                                      + npool_to_livestemn_p           (m) * f_retr_in_nall           * deltim
               I_livestemn_st_p_acc               (m) = I_livestemn_st_p_acc           (m) &
                                                      + npool_to_livestemn_storage_p   (m) * (1._r8 - f_retr_in_nall) * deltim
               AKX_retransn_to_livestemn_st_p_acc (m) = AKX_retransn_to_livestemn_st_p_acc (m) &
                                                      + npool_to_livestemn_storage_p   (m) * f_retr_in_nall           * deltim
               I_deadstemn_p_acc                  (m) = I_deadstemn_p_acc              (m) &
                                                      + npool_to_deadstemn_p           (m) * (1._r8 - f_retr_in_nall) * deltim
               AKX_retransn_to_deadstemn_p_acc    (m) = AKX_retransn_to_deadstemn_p_acc    (m) &
                                                      + npool_to_deadstemn_p           (m) * f_retr_in_nall           * deltim
               I_deadstemn_st_p_acc               (m) = I_deadstemn_st_p_acc           (m) &
                                                      + npool_to_deadstemn_storage_p   (m) * (1._r8 - f_retr_in_nall) * deltim
               AKX_retransn_to_deadstemn_st_p_acc (m) = AKX_retransn_to_deadstemn_st_p_acc (m) &
                                                      + npool_to_deadstemn_storage_p   (m) * f_retr_in_nall           * deltim
               I_livecrootn_p_acc                 (m) = I_livecrootn_p_acc             (m) &
                                                      + npool_to_livecrootn_p          (m) * (1._r8 - f_retr_in_nall) * deltim
               AKX_retransn_to_livecrootn_p_acc   (m) = AKX_retransn_to_livecrootn_p_acc   (m) &
                                                      + npool_to_livecrootn_p          (m) * f_retr_in_nall           * deltim
               I_livecrootn_st_p_acc              (m) = I_livecrootn_st_p_acc          (m) &
                                                      + npool_to_livecrootn_storage_p  (m) * (1._r8 - f_retr_in_nall) * deltim
               AKX_retransn_to_livecrootn_st_p_acc(m) = AKX_retransn_to_livecrootn_st_p_acc(m) &
                                                      + npool_to_livecrootn_storage_p  (m) * f_retr_in_nall           * deltim
               I_deadcrootn_p_acc                 (m) = I_deadcrootn_p_acc             (m) &
                                                      + npool_to_deadcrootn_p          (m) * (1._r8 - f_retr_in_nall) * deltim
               AKX_retransn_to_deadcrootn_p_acc   (m) = AKX_retransn_to_deadcrootn_p_acc   (m) &
                                                      + npool_to_deadcrootn_p          (m) * f_retr_in_nall           * deltim
               I_deadcrootn_st_p_acc              (m) = I_deadcrootn_st_p_acc          (m) &
                                                      + npool_to_deadcrootn_storage_p  (m) * (1._r8 - f_retr_in_nall) * deltim
               AKX_retransn_to_deadcrootn_st_p_acc(m) = AKX_retransn_to_deadcrootn_st_p_acc(m) &
                                                      + npool_to_deadcrootn_storage_p  (m) * f_retr_in_nall           * deltim
            end if
            if (ivt >= npcropmin) then ! skip 2 generic crops
               I_livestemn_p_acc                  (m) = I_livestemn_p_acc          (m) &
                                                      + npool_to_livestemn_p       (m) * (1._r8 - f_retr_in_nall) * deltim
               AKX_retransn_to_livestemn_p_acc    (m) = AKX_retransn_to_livestemn_p_acc   (m) &
                                                      + npool_to_livestemn_p       (m) * f_retr_in_nall           * deltim
               I_livestemn_st_p_acc               (m) = I_livestemn_st_p_acc       (m) &
                                                  + npool_to_livestemn_storage_p   (m) * (1._r8 - f_retr_in_nall) * deltim
               AKX_retransn_to_livestemn_st_p_acc (m) = AKX_retransn_to_livestemn_st_p_acc(m) &
                                                  + npool_to_livestemn_storage_p   (m) * f_retr_in_nall           * deltim
               I_grainn_p_acc                     (m) = I_grainn_p_acc             (m) &
                                                  + npool_to_grainn_p              (m) * (1._r8 - f_retr_in_nall) * deltim
               AKX_retransn_to_grainn_p_acc       (m) = AKX_retransn_to_grainn_p_acc      (m) &
                                                  + npool_to_grainn_p              (m) * f_retr_in_nall           * deltim
               I_grainn_st_p_acc                  (m) = I_grainn_st_p_acc          (m) &
                                                  + npool_to_grainn_storage_p      (m) * (1._r8 - f_retr_in_nall) * deltim
               AKX_retransn_to_grainn_st_p_acc    (m) = AKX_retransn_to_grainn_st_p_acc   (m) &
                                                  + npool_to_grainn_storage_p      (m) * f_retr_in_nall           * deltim
            end if
         end if
#endif
!
! if(i .eq.  123226)print*,'veg nitrogen gain from npool',&
!      sum(pftfrac(ps:pe)*deltim*(npool_to_leafn_p(ps:pe)+npool_to_leafn_storage_p(ps:pe))),&
!      sum(pftfrac(ps:pe)*(npool_to_frootn_p(ps:pe)+npool_to_frootn_storage_p(ps:pe)))*deltim,&
!      sum(pftfrac(ps:pe)*(npool_to_livestemn_p(ps:pe)+npool_to_livestemn_storage_p(ps:pe)))*deltim,&
!      sum(pftfrac(ps:pe)*(npool_to_deadstemn_p(ps:pe)+npool_to_deadstemn_storage_p(ps:pe)))*deltim,&
!      sum(pftfrac(ps:pe)*(npool_to_livecrootn_p(ps:pe)+npool_to_livecrootn_storage_p(ps:pe)))*deltim,&
!      sum(pftfrac(ps:pe)*(npool_to_deadcrootn_p(ps:pe)+npool_to_deadcrootn_storage_p(ps:pe)))*deltim
!if(i .eq. 123226)print*,'vegn after allocation,total',&
!             sum((leafn_p(ps:pe)              + leafn_storage_p(ps:pe))*pftfrac(ps:pe))+&
!             sum((frootn_p(ps:pe)             + frootn_storage_p(ps:pe))*pftfrac(ps:pe))+&
!             sum((livestemn_p(ps:pe)          + livestemn_storage_p(ps:pe))*pftfrac(ps:pe))+&
!             sum((deadstemn_p(ps:pe)          + deadstemn_storage_p(ps:pe))*pftfrac(ps:pe))+&
!             sum((livecrootn_p(ps:pe)         + livecrootn_storage_p(ps:pe))*pftfrac(ps:pe))+&
!             sum((deadcrootn_p(ps:pe)         + deadcrootn_storage_p(ps:pe))*pftfrac(ps:pe))
!       +retransn_to_npool_p(ps:pe)+free_retransn_to_npool_p(ps:pe)
!if(i .eq. 123226)print*,'vegn after allocation',sum(npool_p(ps:pe)*pftfrac(ps:pe)),&
!             sum((leafn_p(ps:pe)             + leafn_storage_p(ps:pe))*pftfrac(ps:pe)),&
!             sum((frootn_p(ps:pe)             + frootn_storage_p(ps:pe))*pftfrac(ps:pe)),&
!             sum((livestemn_p(ps:pe)          + livestemn_storage_p(ps:pe))*pftfrac(ps:pe)),&
!             sum((deadstemn_p(ps:pe)          + deadstemn_storage_p(ps:pe))*pftfrac(ps:pe)),&
!             sum((livecrootn_p(ps:pe)         + livecrootn_storage_p(ps:pe))*pftfrac(ps:pe)),&
!             sum((deadcrootn_p(ps:pe)         + deadcrootn_storage_p(ps:pe))*pftfrac(ps:pe))
!                + leafn_xfer_p(ps:pe)        + frootn_xfer_p(ps:pe)        + livestemn_xfer_p(ps:pe) &
!                + deadstemn_xfer_p(ps:pe)    + livecrootn_xfer_p(ps:pe)    + deadcrootn_xfer_p(ps:pe) &
!                + npool_p(ps:pe)             + retransn_p(ps:pe) &
!                + grainn_p(ps:pe)            + grainn_storage_p(ps:pe)     + grainn_xfer_p(ps:pe) &
!                + cropseedn_deficit_p(ps:pe)
!                 )*pftfrac(ps:pe))
!      do m = ps, pe
!         ivt=pftclass(m)
         ! move storage pools into transfer pools
         leafn_storage_p(m)  = leafn_storage_p(m)  - leafn_storage_to_xfer_p(m)*deltim
         leafn_xfer_p(m)     = leafn_xfer_p(m)     + leafn_storage_to_xfer_p(m)*deltim
         frootn_storage_p(m) = frootn_storage_p(m) - frootn_storage_to_xfer_p(m)*deltim
         frootn_xfer_p(m)    = frootn_xfer_p(m)    + frootn_storage_to_xfer_p(m)*deltim

         if (woody(ivt) == 1) then
            livestemn_storage_p(m)  = livestemn_storage_p(m)  - livestemn_storage_to_xfer_p(m)*deltim
            livestemn_xfer_p(m)     = livestemn_xfer_p(m)     + livestemn_storage_to_xfer_p(m)*deltim
            deadstemn_storage_p(m)  = deadstemn_storage_p(m)  - deadstemn_storage_to_xfer_p(m)*deltim
            deadstemn_xfer_p(m)     = deadstemn_xfer_p(m)     + deadstemn_storage_to_xfer_p(m)*deltim
            livecrootn_storage_p(m) = livecrootn_storage_p(m) - livecrootn_storage_to_xfer_p(m)*deltim
            livecrootn_xfer_p(m)    = livecrootn_xfer_p(m)    + livecrootn_storage_to_xfer_p(m)*deltim
            deadcrootn_storage_p(m) = deadcrootn_storage_p(m) - deadcrootn_storage_to_xfer_p(m)*deltim
            deadcrootn_xfer_p(m)    = deadcrootn_xfer_p(m)    + deadcrootn_storage_to_xfer_p(m)*deltim
         end if

         if (ivt >= npcropmin) then ! skip 2 generic crops
         ! lines here for consistency; the transfer terms are zero
            livestemn_storage_p(m)  = livestemn_storage_p(m) - livestemn_storage_to_xfer_p(m)*deltim
            livestemn_xfer_p(m)     = livestemn_xfer_p(m)    + livestemn_storage_to_xfer_p(m)*deltim
            grainn_storage_p(m)     = grainn_storage_p(m)    - grainn_storage_to_xfer_p(m)*deltim
            grainn_xfer_p(m)        = grainn_xfer_p(m)       + grainn_storage_to_xfer_p(m)*deltim
         end if

#ifdef SASU
         AKX_leafn_st_to_leafn_xf_p_acc             (m) = AKX_leafn_st_to_leafn_xf_p_acc          (m) + leafn_storage_to_xfer_p     (m) * deltim
         AKX_leafn_st_exit_p_acc                    (m) = AKX_leafn_st_exit_p_acc                 (m) + leafn_storage_to_xfer_p     (m) * deltim
         AKX_frootn_st_to_frootn_xf_p_acc           (m) = AKX_frootn_st_to_frootn_xf_p_acc        (m) + frootn_storage_to_xfer_p    (m) * deltim
         AKX_frootn_st_exit_p_acc                   (m) = AKX_frootn_st_exit_p_acc                (m) + frootn_storage_to_xfer_p    (m) * deltim
         if(woody(ivt) == 1) then
            AKX_livestemn_st_to_livestemn_xf_p_acc  (m) = AKX_livestemn_st_to_livestemn_xf_p_acc  (m) + livestemn_storage_to_xfer_p (m) * deltim
            AKX_livestemn_st_exit_p_acc             (m) = AKX_livestemn_st_exit_p_acc             (m) + livestemn_storage_to_xfer_p (m) * deltim
            AKX_deadstemn_st_to_deadstemn_xf_p_acc  (m) = AKX_deadstemn_st_to_deadstemn_xf_p_acc  (m) + deadstemn_storage_to_xfer_p (m) * deltim
            AKX_deadstemn_st_exit_p_acc             (m) = AKX_deadstemn_st_exit_p_acc             (m) + deadstemn_storage_to_xfer_p (m) * deltim
            AKX_livecrootn_st_to_livecrootn_xf_p_acc(m) = AKX_livecrootn_st_to_livecrootn_xf_p_acc(m) + livecrootn_storage_to_xfer_p(m) * deltim
            AKX_livecrootn_st_exit_p_acc            (m) = AKX_livecrootn_st_exit_p_acc            (m) + livecrootn_storage_to_xfer_p(m) * deltim
            AKX_deadcrootn_st_to_deadcrootn_xf_p_acc(m) = AKX_deadcrootn_st_to_deadcrootn_xf_p_acc(m) + deadcrootn_storage_to_xfer_p(m) * deltim
            AKX_deadcrootn_st_exit_p_acc            (m) = AKX_deadcrootn_st_exit_p_acc            (m) + deadcrootn_storage_to_xfer_p(m) * deltim
         end if
         if( ivt >= npcropmin) then
            AKX_livestemn_st_to_livestemn_xf_p_acc  (m) = AKX_livestemn_st_to_livestemn_xf_p_acc  (m) + livestemn_storage_to_xfer_p (m) * deltim
            AKX_livestemn_st_exit_p_acc             (m) = AKX_livestemn_st_exit_p_acc             (m) + livestemn_storage_to_xfer_p (m) * deltim
            AKX_grainn_st_to_grainn_xf_p_acc        (m) = AKX_grainn_st_to_grainn_xf_p_acc        (m) + grainn_storage_to_xfer_p    (m) * deltim
            AKX_grainn_st_exit_p_acc                (m) = AKX_grainn_st_exit_p_acc                (m) + grainn_storage_to_xfer_p    (m) * deltim
         end if
#endif
      end do ! end pft loop
!    if(i .eq.  123226)print*,'veg nitrogen loss from phen',sum(pftfrac(ps:pe)*deltim*&
!       (leafn_to_litter_p(ps:pe)+frootn_to_litter_p(ps:pe)+livestemn_to_litter_p(ps:pe)))
!
!    if(i .eq.  123226)print*,'veg nitrogen gain from alloc',sum(pftfrac(ps:pe)*deltim*&
!       (npool_to_leafn_p(ps:pe)+npool_to_leafn_storage_p(ps:pe)&
!       +npool_to_frootn_p(ps:pe)+npool_to_frootn_storage_p(ps:pe)&
!       +npool_to_livestemn_p(ps:pe)+npool_to_livestemn_storage_p(ps:pe)&
!       +npool_to_deadstemn_p(ps:pe)+npool_to_deadstemn_storage_p(ps:pe)&
!       +npool_to_livecrootn_p(ps:pe)+npool_to_livecrootn_storage_p(ps:pe)&
!       +npool_to_deadcrootn_p(ps:pe)+npool_to_deadcrootn_storage_p(ps:pe)&
!       +npool_to_grainn_p(ps:pe)+npool_to_grainn_storage_p(ps:pe)))

!if(i .eq. 123226)print*,'vegn after NStateUpdate1',&
!             sum((leafn_p(ps:pe)             + frootn_p(ps:pe)             + livestemn_p(ps:pe) &
!                + deadstemn_p(ps:pe)         + livecrootn_p(ps:pe)         + deadcrootn_p(ps:pe) &
!                + leafn_storage_p(ps:pe)     + frootn_storage_p(ps:pe)     + livestemn_storage_p(ps:pe) &
!                + deadstemn_storage_p(ps:pe) + livecrootn_storage_p(ps:pe) + deadcrootn_storage_p(ps:pe) &
!                + leafn_xfer_p(ps:pe)        + frootn_xfer_p(ps:pe)        + livestemn_xfer_p(ps:pe) &
!                + deadstemn_xfer_p(ps:pe)    + livecrootn_xfer_p(ps:pe)    + deadcrootn_xfer_p(ps:pe) &
!                + npool_p(ps:pe)             + retransn_p(ps:pe) &
!                + grainn_p(ps:pe)            + grainn_storage_p(ps:pe)     + grainn_xfer_p(ps:pe) &
!                + cropseedn_deficit_p(ps:pe))*pftfrac(ps:pe))

!if(i .eq. 123226)print*,'nsource_sink after NStateUpdate1',&
!  sum(decomp_npools_sourcesink(1:nl_soil,i_met_lit,i)*dz_soi(1:nl_soil)) &
! +sum(decomp_npools_sourcesink(1:nl_soil,i_cel_lit,i)*dz_soi(1:nl_soil)) &
! +sum(decomp_npools_sourcesink(1:nl_soil,i_cwd,i)*dz_soi(1:nl_soil)) &
! +sum(decomp_npools_sourcesink(1:nl_soil,i_soil1,i)*dz_soi(1:nl_soil)) &
! +sum(decomp_npools_sourcesink(1:nl_soil,i_soil2,i)*dz_soi(1:nl_soil)) &
! +sum(decomp_npools_sourcesink(1:nl_soil,i_soil3,i)*dz_soi(1:nl_soil))

end subroutine NStateUpdate1

end module bgc_CNNStateUpdate1Mod
