module bgc_CNNStateUpdate1Mod

use precision
use MOD_PFTimeInvars, only: pftclass
use PFT_Const, only: woody
use MOD_TimeInvariants, only: &
  ! bgc constants
           donor_pool, receiver_pool, i_met_lit, i_cel_lit, i_lig_lit, i_cwd, i_soil1, i_soil2, i_soil3

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
           cropseedn_deficit_p, retransn_p          

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
           npool_to_grainn_p    , npool_to_grainn_storage_p    

implicit none

public NStateUpdate1

contains

subroutine NStateUpdate1 &
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

         
         ! allocation fluxes
         retransn_p(m)        = retransn_p(m)       - retransn_to_npool_p(m)*deltim
         retransn_p(m)        = retransn_p(m)       - free_retransn_to_npool_p(m)*deltim !how is retransn a state? 
         leafn_p(m)           = leafn_p(m)          + npool_to_leafn_p(m)*deltim
         leafn_storage_p(m)   = leafn_storage_p(m)  + npool_to_leafn_storage_p(m)*deltim
         frootn_p(m)          = frootn_p(m)         + npool_to_frootn_p(m)*deltim
         frootn_storage_p(m)  = frootn_storage_p(m) + npool_to_frootn_storage_p(m)*deltim

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

         if (ivt >= npcropmin) then ! skip 2 generic crops
            livestemn_p(m)          = livestemn_p(m)          + npool_to_livestemn_p(m)*deltim
            livestemn_storage_p(m)  = livestemn_storage_p(m)  + npool_to_livestemn_storage_p(m)*deltim
            grainn_p(m)             = grainn_p(m)             + npool_to_grainn_p(m)*deltim
            grainn_storage_p(m)     = grainn_storage_p(m)     + npool_to_grainn_storage_p(m)*deltim
         end if

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
      end do ! end pft loop

end subroutine NStateUpdate1

end module bgc_CNNStateUpdate1Mod
