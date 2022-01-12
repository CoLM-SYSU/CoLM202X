#include <define.h>

module bgc_CNNStateUpdate3Mod
use precision
use MOD_TimeInvariants, only: &
           i_met_lit,i_cel_lit,i_lig_lit ,i_cwd, i_soil1, i_soil2, i_soil3
use MOD_TimeVariables, only: &
    ! decomposition pools & fluxes variables (inout)
           decomp_npools_vr, sminn_vr, smin_no3_vr, smin_nh4_vr

use MOD_1D_Fluxes, only: &
           m_decomp_npools_to_fire_vr, &
           fire_mortality_to_met_n, fire_mortality_to_cel_n, &
           fire_mortality_to_lig_n, fire_mortality_to_cwdn , &

    ! mineral nitrogen pools & fluxes variables (inout)
           sminn_leached_vr, smin_no3_leached_vr, smin_no3_runoff_vr

use MOD_PFTimeVars, only: &
    ! vegetation nitrogen state variables (inout)
           leafn_p            , leafn_storage_p     , leafn_xfer_p     , &
           frootn_p           , frootn_storage_p    , frootn_xfer_p    , &
           livestemn_p        , livestemn_storage_p , livestemn_xfer_p , &
           deadstemn_p        , deadstemn_storage_p , deadstemn_xfer_p , &
           livecrootn_p       , livecrootn_storage_p, livecrootn_xfer_p, &
           deadcrootn_p       , deadcrootn_storage_p, deadcrootn_xfer_p, &
           retransn_p

use MOD_1D_PFTFluxes, only: &
    ! vegetation nitrogen flux variables
           m_leafn_to_fire_p        , m_leafn_storage_to_fire_p     , m_leafn_xfer_to_fire_p     , &
           m_frootn_to_fire_p       , m_frootn_storage_to_fire_p    , m_frootn_xfer_to_fire_p    , &
           m_livestemn_to_fire_p    , m_livestemn_storage_to_fire_p , m_livestemn_xfer_to_fire_p , &
           m_deadstemn_to_fire_p    , m_deadstemn_storage_to_fire_p , m_deadstemn_xfer_to_fire_p , &
           m_livecrootn_to_fire_p   , m_livecrootn_storage_to_fire_p, m_livecrootn_xfer_to_fire_p, &
           m_deadcrootn_to_fire_p   , m_deadcrootn_storage_to_fire_p, m_deadcrootn_xfer_to_fire_p, &
           m_livestemn_to_deadstemn_fire_p , m_livecrootn_to_deadcrootn_fire_p , &
           m_retransn_to_fire_p, &

           m_leafn_to_litter_fire_p        , m_leafn_storage_to_litter_fire_p     , m_leafn_xfer_to_litter_fire_p     , &
           m_frootn_to_litter_fire_p       , m_frootn_storage_to_litter_fire_p    , m_frootn_xfer_to_litter_fire_p    , &
           m_livestemn_to_litter_fire_p    , m_livestemn_storage_to_litter_fire_p , m_livestemn_xfer_to_litter_fire_p , &
           m_deadstemn_to_litter_fire_p    , m_deadstemn_storage_to_litter_fire_p , m_deadstemn_xfer_to_litter_fire_p , &
           m_livecrootn_to_litter_fire_p   , m_livecrootn_storage_to_litter_fire_p, m_livecrootn_xfer_to_litter_fire_p, &
           m_deadcrootn_to_litter_fire_p   , m_deadcrootn_storage_to_litter_fire_p, m_deadcrootn_xfer_to_litter_fire_p, &
           m_retransn_to_litter_fire_p
           

implicit none

public NStateUpdate3

contains

subroutine NStateUpdate3(i, ps, pe, deltim, nl_soil, ndecomp_pools, dz_soi)

integer ,intent(in) :: i
integer ,intent(in) :: ps
integer ,intent(in) :: pe
real(r8),intent(in) :: deltim
integer ,intent(in) :: nl_soil
integer ,intent(in) :: ndecomp_pools
real(r8),intent(in) :: dz_soi(1:nl_soil)

integer j,l

!if(i .eq. 123226)print*,'soiln before NStateUpdate3',&
!      sum(decomp_npools_vr(1:nl_soil,i_met_lit,i)*dz_soi(1:nl_soil)) &
!    + sum(decomp_npools_vr(1:nl_soil,i_cel_lit,i)*dz_soi(1:nl_soil)) &
!    + sum(decomp_npools_vr(1:nl_soil,i_lig_lit,i)*dz_soi(1:nl_soil)) &
!    + sum(decomp_npools_vr(1:nl_soil,i_cwd,i)*dz_soi(1:nl_soil)) &
!    + sum(decomp_npools_vr(1:nl_soil,i_soil1,i)*dz_soi(1:nl_soil)) &
!    + sum(decomp_npools_vr(1:nl_soil,i_soil2,i)*dz_soi(1:nl_soil)) &
!    + sum(decomp_npools_vr(1:nl_soil,i_soil3,i)*dz_soi(1:nl_soil)) 
      do j = 1, nl_soil

#ifndef NITRIF
         ! mineral N loss due to leaching
         sminn_vr(j,i) = sminn_vr(j,i) - sminn_leached_vr(j,i) * deltim
#else
         ! mineral N loss due to leaching and runoff
         smin_no3_vr(j,i) = max( smin_no3_vr(j,i) - &
                    ( smin_no3_leached_vr(j,i) + smin_no3_runoff_vr(j,i) ) * deltim, 0._r8)

         sminn_vr(j,i) = smin_no3_vr(j,i) + smin_nh4_vr(j,i)
#endif
         
            ! column level nitrogen fluxes from fire
            ! patch-level wood to column-level CWD (uncombusted wood)
!            if (.not. use_soil_matrixcn)then
#ifdef FIRE
         decomp_npools_vr(j,i_cwd,i) = decomp_npools_vr(j,i_cwd,i) + &
                 fire_mortality_to_cwdn(j,i) * deltim

            ! patch-level wood to column-level litter (uncombusted wood)
         decomp_npools_vr(j,i_met_lit,i) = decomp_npools_vr(j,i_met_lit,i) + &
                 fire_mortality_to_met_n(j,i)* deltim
         decomp_npools_vr(j,i_cel_lit,i) = decomp_npools_vr(j,i_cel_lit,i) + &
                 fire_mortality_to_cel_n(j,i)* deltim
         decomp_npools_vr(j,i_lig_lit,i) = decomp_npools_vr(j,i_lig_lit,i) + &
                 fire_mortality_to_lig_n(j,i)* deltim
!            else
!               matrix_Ninput%V(c,j+(i_cwd-1)*nlevdecomp) = matrix_Ninput%V(c,j+(i_cwd-1)*nlevdecomp) + &
!                 fire_mortality_n_to_cwdn(j) * deltim

            ! patch-level wood to column-level litter (uncombusted wood)
!               matrix_Ninput%V(c,j+(i_met_lit-1)*nlevdecomp) = matrix_Ninput%V(c,j+(i_met_lit-1)*nlevdecomp) + &
!                 fire_mortality_to_met_n(j)* deltim
!               matrix_Ninput%V(c,j+(i_cel_lit-1)*nlevdecomp) = matrix_Ninput%V(c,j+(i_cel_lit-1)*nlevdecomp) + &
!                 fire_mortality_to_cel_n(j)* deltim
!               matrix_Ninput%V(c,j+(i_lig_lit-1)*nlevdecomp) = matrix_Ninput%V(c,j+(i_lig_lit-1)*nlevdecomp) + &
!                 fire_mortality_to_lig_n(j)* deltim
!            end if ! not use_soil_matrix
#endif
      end do
!      if(i .eq. 147958)print*,'sminn after leached',sum(sminn_vr(1:nl_soil,i)*dz_soi(1:nl_soil))
!      if(i .eq. 147958)print*,'leached sminn',sum(sminn_leached_vr(1:nl_soil,i)*deltim*dz_soi(1:nl_soil))
#ifdef FIRE
      ! litter and CWD losses to fire
!      if(.not. use_soil_matrixcn)then
      do l = 1, ndecomp_pools
         do j = 1, nl_soil
            decomp_npools_vr(j,l,i) = decomp_npools_vr(j,l,i) - &
                    m_decomp_npools_to_fire_vr(j,l,i) * deltim
         end do
      end do
!      end if ! not use_soil_matrixcn

      do m = ps , pe
!        if(.not. use_matrixcn)then 
         !from fire displayed pools
         leafn_p(m) =  leafn_p(m) -                           &
                 m_leafn_to_fire_p(m) * deltim
         frootn_p(m) =  frootn_p(m) -                         &
                 m_frootn_to_fire_p(m) * deltim
         livestemn_p(m) =  livestemn_p(m) -                   &
                 m_livestemn_to_fire_p(m) * deltim
         deadstemn_p(m) =  deadstemn_p(m) -                   &
                 m_deadstemn_to_fire_p(m) * deltim
         livecrootn_p(m) =  livecrootn_p(m) -                 &
                 m_livecrootn_to_fire_p(m) * deltim
         deadcrootn_p(m) =  deadcrootn_p(m) -                 &
                 m_deadcrootn_to_fire_p(m) * deltim
   
         leafn_p(m) =  leafn_p(m) -                             &
                 m_leafn_to_litter_fire_p(m) * deltim
         frootn_p(m) =  frootn_p(m) -                           &
                 m_frootn_to_litter_fire_p(m) * deltim
         livestemn_p(m) =  livestemn_p(m) -                     &
                 m_livestemn_to_litter_fire_p(m) * deltim  -  &
                 m_livestemn_to_deadstemn_fire_p(m) * deltim
         deadstemn_p(m) =  deadstemn_p(m) -                     &
                  m_deadstemn_to_litter_fire_p(m) * deltim +  &
                  m_livestemn_to_deadstemn_fire_p(m) * deltim
         livecrootn_p(m) =  livecrootn_p(m) -                   &
                  m_livecrootn_to_litter_fire_p(m) * deltim - &
                  m_livecrootn_to_deadcrootn_fire_p(m) * deltim
         deadcrootn_p(m) =  deadcrootn_p(m) -                   &
                  m_deadcrootn_to_litter_fire_p(m) * deltim + &
                  m_livecrootn_to_deadcrootn_fire_p(m) * deltim 
   
         ! storage pools
         leafn_storage_p(m) =  leafn_storage_p(m) -           &
                 m_leafn_storage_to_fire_p(m) * deltim
         frootn_storage_p(m) =  frootn_storage_p(m) -         &
                 m_frootn_storage_to_fire_p(m) * deltim
         livestemn_storage_p(m) =  livestemn_storage_p(m) -   &
                 m_livestemn_storage_to_fire_p(m) * deltim
         deadstemn_storage_p(m) =  deadstemn_storage_p(m) -   &
                 m_deadstemn_storage_to_fire_p(m) * deltim
         livecrootn_storage_p(m) =  livecrootn_storage_p(m) - &
                 m_livecrootn_storage_to_fire_p(m) * deltim
         deadcrootn_storage_p(m) =  deadcrootn_storage_p(m) - &
                 m_deadcrootn_storage_to_fire_p(m) * deltim
   
         leafn_storage_p(m) =  leafn_storage_p(m) -           &
                 m_leafn_storage_to_litter_fire_p(m) * deltim
         frootn_storage_p(m) =  frootn_storage_p(m) -         &
                 m_frootn_storage_to_litter_fire_p(m) * deltim
         livestemn_storage_p(m) =  livestemn_storage_p(m) -   &
                 m_livestemn_storage_to_litter_fire_p(m) * deltim
         deadstemn_storage_p(m) =  deadstemn_storage_p(m) -   &
                 m_deadstemn_storage_to_litter_fire_p(m) * deltim
         livecrootn_storage_p(m) =  livecrootn_storage_p(m) - &
                 m_livecrootn_storage_to_litter_fire_p(m) * deltim
         deadcrootn_storage_p(m) =  deadcrootn_storage_p(m) - &
                 m_deadcrootn_storage_to_litter_fire_p(m) * deltim
   

         ! transfer pools
         leafn_xfer_p(m) =  leafn_xfer_p(m) -                 &
                 m_leafn_xfer_to_fire_p(m) * deltim
         frootn_xfer_p(m) =  frootn_xfer_p(m) -               &
                 m_frootn_xfer_to_fire_p(m) * deltim
         livestemn_xfer_p(m) =  livestemn_xfer_p(m) -         &
                 m_livestemn_xfer_to_fire_p(m) * deltim
         deadstemn_xfer_p(m) =  deadstemn_xfer_p(m) -         &
                 m_deadstemn_xfer_to_fire_p(m) * deltim
         livecrootn_xfer_p(m) =  livecrootn_xfer_p(m) -       &
                 m_livecrootn_xfer_to_fire_p(m) * deltim
         deadcrootn_xfer_p(m) =  deadcrootn_xfer_p(m) -       &
                 m_deadcrootn_xfer_to_fire_p(m) * deltim

         leafn_xfer_p(m) =  leafn_xfer_p(m) -                 &
                 m_leafn_xfer_to_litter_fire_p(m) * deltim
         frootn_xfer_p(m) =  frootn_xfer_p(m) -               &
                 m_frootn_xfer_to_litter_fire_p(m) * deltim
         livestemn_xfer_p(m) =  livestemn_xfer_p(m) -         &
                 m_livestemn_xfer_to_litter_fire_p(m) * deltim
         deadstemn_xfer_p(m) =  deadstemn_xfer_p(m) -         &
                 m_deadstemn_xfer_to_litter_fire_p(m) * deltim
         livecrootn_xfer_p(m) =  livecrootn_xfer_p(m) -       &
                 m_livecrootn_xfer_to_litter_fire_p(m) * deltim
         deadcrootn_xfer_p(m) =  deadcrootn_xfer_p(m) -       &
                 m_deadcrootn_xfer_to_litter_fire_p(m) * deltim

      ! retranslocated N pool
         retransn_p(m) =  retransn_p(m) -                     &
                 m_retransn_to_fire_p(m) * deltim
         retransn_p(m) =  retransn_p(m) -                     &
                 m_retransn_to_litter_fire_p(m) * deltim
!         end if !.not. use_matrixcn
      end do
!if(i .eq. 123226)print*,'soiln after NStateUpdate3',&
!      sum(decomp_npools_vr(1:nl_soil,i_met_lit,i)*dz_soi(1:nl_soil)) &
!    + sum(decomp_npools_vr(1:nl_soil,i_cel_lit,i)*dz_soi(1:nl_soil)) &
!    + sum(decomp_npools_vr(1:nl_soil,i_lig_lit,i)*dz_soi(1:nl_soil)) &
!    + sum(decomp_npools_vr(1:nl_soil,i_cwd,i)*dz_soi(1:nl_soil)) &
!    + sum(decomp_npools_vr(1:nl_soil,i_soil1,i)*dz_soi(1:nl_soil)) &
!    + sum(decomp_npools_vr(1:nl_soil,i_soil2,i)*dz_soi(1:nl_soil)) &
!    + sum(decomp_npools_vr(1:nl_soil,i_soil3,i)*dz_soi(1:nl_soil)) 
#endif

end subroutine NStateUpdate3

end module bgc_CNNStateUpdate3Mod
