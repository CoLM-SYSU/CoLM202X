#include <define.h>
#ifdef BGC
module bgc_CNCStateUpdate3Mod
use precision
use MOD_BGCTimeInvars, only: &
           i_met_lit,i_cel_lit,i_lig_lit ,i_cwd 
use MOD_BGCTimeVars, only: &
    ! decomposition pools & fluxes variables (inout)
           decomp_cpools_vr 

use MOD_1D_BGCFluxes, only: &
           m_decomp_cpools_to_fire_vr, &
           fire_mortality_to_met_c, fire_mortality_to_cel_c, &
           fire_mortality_to_lig_c, fire_mortality_to_cwdc

use MOD_BGCPFTimeVars, only: &
    ! vegetation carbon state variables (inout)
           leafc_p            , leafc_storage_p     , leafc_xfer_p     , &
           frootc_p           , frootc_storage_p    , frootc_xfer_p    , &
           livestemc_p        , livestemc_storage_p , livestemc_xfer_p , &
           deadstemc_p        , deadstemc_storage_p , deadstemc_xfer_p , &
           livecrootc_p       , livecrootc_storage_p, livecrootc_xfer_p, &
           deadcrootc_p       , deadcrootc_storage_p, deadcrootc_xfer_p, &
           gresp_storage_p    , gresp_xfer_p      

use MOD_1D_BGCPFTFluxes, only: &
    ! vegetation carbon flux variables
           m_leafc_to_fire_p        , m_leafc_storage_to_fire_p     , m_leafc_xfer_to_fire_p     , &
           m_frootc_to_fire_p       , m_frootc_storage_to_fire_p    , m_frootc_xfer_to_fire_p    , &
           m_livestemc_to_fire_p    , m_livestemc_storage_to_fire_p , m_livestemc_xfer_to_fire_p , &
           m_deadstemc_to_fire_p    , m_deadstemc_storage_to_fire_p , m_deadstemc_xfer_to_fire_p , &
           m_livecrootc_to_fire_p   , m_livecrootc_storage_to_fire_p, m_livecrootc_xfer_to_fire_p, &
           m_deadcrootc_to_fire_p   , m_deadcrootc_storage_to_fire_p, m_deadcrootc_xfer_to_fire_p, &
           m_livestemc_to_deadstemc_fire_p , m_livecrootc_to_deadcrootc_fire_p , &
           m_gresp_storage_to_fire_p       , m_gresp_xfer_to_fire_p            , &

           m_leafc_to_litter_fire_p        , m_leafc_storage_to_litter_fire_p     , m_leafc_xfer_to_litter_fire_p     , &
           m_frootc_to_litter_fire_p       , m_frootc_storage_to_litter_fire_p    , m_frootc_xfer_to_litter_fire_p    , &
           m_livestemc_to_litter_fire_p    , m_livestemc_storage_to_litter_fire_p , m_livestemc_xfer_to_litter_fire_p , &
           m_deadstemc_to_litter_fire_p    , m_deadstemc_storage_to_litter_fire_p , m_deadstemc_xfer_to_litter_fire_p , &
           m_livecrootc_to_litter_fire_p   , m_livecrootc_storage_to_litter_fire_p, m_livecrootc_xfer_to_litter_fire_p, &
           m_deadcrootc_to_litter_fire_p   , m_deadcrootc_storage_to_litter_fire_p, m_deadcrootc_xfer_to_litter_fire_p, &
           m_gresp_storage_to_litter_fire_p, m_gresp_xfer_to_litter_fire_p        
           
implicit none

public CStateUpdate3

contains

subroutine CStateUpdate3(i, ps, pe, deltim, nl_soil, ndecomp_pools)

integer ,intent(in) :: i
integer ,intent(in) :: ps
integer ,intent(in) :: pe
real(r8),intent(in) :: deltim
integer ,intent(in) :: nl_soil
integer ,intent(in) :: ndecomp_pools

integer j,l,m

!       if(use_matrixcn)then
       do j = 1, nl_soil
      ! patch-level wood to column-level CWD (uncombusted wood)
!          if (.not. use_soil_matrixcn) then
          decomp_cpools_vr(j,i_cwd,i) = decomp_cpools_vr(j,i_cwd,i) + &
                  fire_mortality_to_cwdc(j,i) * deltim

      ! patch-level wood to column-level litter (uncombusted wood)
          decomp_cpools_vr(j,i_met_lit,i) = decomp_cpools_vr(j,i_met_lit,i) + &
                  fire_mortality_to_met_c(j,i)* deltim
          decomp_cpools_vr(j,i_cel_lit,i) = decomp_cpools_vr(j,i_cel_lit,i) + &
                  fire_mortality_to_cel_c(j,i)* deltim
          decomp_cpools_vr(j,i_lig_lit,i) = decomp_cpools_vr(j,i_lig_lit,i) + &
                  fire_mortality_to_lig_c(j,i)* deltim
!             else
             ! patch-level wood to column-level CWD (uncombusted wood)
!                  cf_soil%matrix_Cinput%V(c,j+(i_cwd-1)*nlevdecomp) = cf_soil%matrix_Cinput%V(c,j+(i_cwd-1)*nlevdecomp) + &
!                    fire_mortality_c_to_cwdc(j) * deltim

             ! patch-level wood to column-level litter (uncombusted wood)
!                  cf_soil%matrix_Cinput%V(c,j+(i_met_lit-1)*nlevdecomp) = cf_soil%matrix_Cinput%V(c,j+(i_met_lit-1)*nlevdecomp) + &
!                    m_c_to_met_fire(j)* deltim
!                  cf_soil%matrix_Cinput%V(c,j+(i_cel_lit-1)*nlevdecomp) = cf_soil%matrix_Cinput%V(c,j+(i_cel_lit-1)*nlevdecomp) + &
!                    m_c_to_cel_fire(j)* deltim
!                  cf_soil%matrix_Cinput%V(c,j+(i_lig_lit-1)*nlevdecomp) = cf_soil%matrix_Cinput%V(c,j+(i_lig_lit-1)*nlevdecomp) + &
!                    m_c_to_lig_fire(j)* deltim
!         end if
       end do

       ! litter and CWD losses to fire
       do l = 1, ndecomp_pools
          do j = 1, nl_soil
             decomp_cpools_vr(j,l,i) = decomp_cpools_vr(j,l,i) - &
                   m_decomp_cpools_to_fire_vr(j,l,i) * deltim
          end do
       end do

       ! patch-level carbon fluxes from fire
       do m = ps , pe
          gresp_storage_p(m) = gresp_storage_p(m) -           &
               m_gresp_storage_to_fire_p(m) * deltim
          gresp_storage_p(m) = gresp_storage_p(m) -           &
               m_gresp_storage_to_litter_fire_p(m) * deltim
          gresp_xfer_p(m) = gresp_xfer_p(m) -                 &
               m_gresp_xfer_to_fire_p(m) * deltim
          gresp_xfer_p(m) = gresp_xfer_p(m) -                 &
               m_gresp_xfer_to_litter_fire_p(m) * deltim
   !       if(.not. use_matrixcn)then
           ! displayed pools
          leafc_p(m) = leafc_p(m) -                           &
               m_leafc_to_fire_p(m) * deltim
          leafc_p(m) = leafc_p(m) -                           &
               m_leafc_to_litter_fire_p(m) * deltim
          frootc_p(m) = frootc_p(m) -                         &
               m_frootc_to_fire_p(m) * deltim
          frootc_p(m) = frootc_p(m) -                         &
               m_frootc_to_litter_fire_p(m) * deltim
          livestemc_p(m) = livestemc_p(m) -                   &
               m_livestemc_to_fire_p(m) * deltim
          livestemc_p(m) = livestemc_p(m) -                   &
               m_livestemc_to_litter_fire_p(m) * deltim  -                   &
               m_livestemc_to_deadstemc_fire_p(m) * deltim
          deadstemc_p(m) = deadstemc_p(m) -                   &
               m_deadstemc_to_fire_p(m) * deltim
          deadstemc_p(m) = deadstemc_p(m) -                   &
               m_deadstemc_to_litter_fire_p(m) * deltim  +                   &
               m_livestemc_to_deadstemc_fire_p(m) * deltim
          livecrootc_p(m) = livecrootc_p(m) -                 &
               m_livecrootc_to_fire_p(m) * deltim
          livecrootc_p(m) = livecrootc_p(m) -                 &
               m_livecrootc_to_litter_fire_p(m) * deltim   -                 &
               m_livecrootc_to_deadcrootc_fire_p(m) * deltim
          deadcrootc_p(m) = deadcrootc_p(m) -                 &
               m_deadcrootc_to_fire_p(m) * deltim
          deadcrootc_p(m) = deadcrootc_p(m) -                 &
               m_deadcrootc_to_litter_fire_p(m)* deltim    +                 &
               m_livecrootc_to_deadcrootc_fire_p(m) * deltim

     ! storage pools
          leafc_storage_p(m) = leafc_storage_p(m) -           &
               m_leafc_storage_to_fire_p(m) * deltim
          leafc_storage_p(m) = leafc_storage_p(m) -           &
               m_leafc_storage_to_litter_fire_p(m) * deltim
          frootc_storage_p(m) = frootc_storage_p(m) -         &
               m_frootc_storage_to_fire_p(m) * deltim
          frootc_storage_p(m) = frootc_storage_p(m) -         &
               m_frootc_storage_to_litter_fire_p(m) * deltim
          livestemc_storage_p(m) = livestemc_storage_p(m) -   &
               m_livestemc_storage_to_fire_p(m) * deltim
          livestemc_storage_p(m) = livestemc_storage_p(m) -   &
               m_livestemc_storage_to_litter_fire_p(m) * deltim
          deadstemc_storage_p(m) = deadstemc_storage_p(m) -   &
               m_deadstemc_storage_to_fire_p(m) * deltim
          deadstemc_storage_p(m) = deadstemc_storage_p(m) -   &
               m_deadstemc_storage_to_litter_fire_p(m) * deltim
          livecrootc_storage_p(m) = livecrootc_storage_p(m) - &
               m_livecrootc_storage_to_fire_p(m) * deltim
          livecrootc_storage_p(m) = livecrootc_storage_p(m) - &
               m_livecrootc_storage_to_litter_fire_p(m)* deltim
          deadcrootc_storage_p(m) = deadcrootc_storage_p(m) - &
               m_deadcrootc_storage_to_fire_p(m) * deltim
          deadcrootc_storage_p(m) = deadcrootc_storage_p(m) - &
               m_deadcrootc_storage_to_litter_fire_p(m)* deltim

    ! transfer pools
          leafc_xfer_p(m) = leafc_xfer_p(m) -                 &
               m_leafc_xfer_to_fire_p(m) * deltim
          leafc_xfer_p(m) = leafc_xfer_p(m) -                 &
               m_leafc_xfer_to_litter_fire_p(m) * deltim
          frootc_xfer_p(m) = frootc_xfer_p(m) -               &
               m_frootc_xfer_to_fire_p(m) * deltim
          frootc_xfer_p(m) = frootc_xfer_p(m) -               &
               m_frootc_xfer_to_litter_fire_p(m) * deltim
          livestemc_xfer_p(m) = livestemc_xfer_p(m) -         &
               m_livestemc_xfer_to_fire_p(m) * deltim
          livestemc_xfer_p(m) = livestemc_xfer_p(m) -         &
               m_livestemc_xfer_to_litter_fire_p(m) * deltim
          deadstemc_xfer_p(m) = deadstemc_xfer_p(m) -         &
               m_deadstemc_xfer_to_fire_p(m) * deltim
          deadstemc_xfer_p(m) = deadstemc_xfer_p(m) -         &
               m_deadstemc_xfer_to_litter_fire_p(m) * deltim
          livecrootc_xfer_p(m) = livecrootc_xfer_p(m) -       &
               m_livecrootc_xfer_to_fire_p(m) * deltim
          livecrootc_xfer_p(m) = livecrootc_xfer_p(m) -       &
               m_livecrootc_xfer_to_litter_fire_p(m)* deltim
          deadcrootc_xfer_p(m) = deadcrootc_xfer_p(m) -       &
               m_deadcrootc_xfer_to_fire_p(m) * deltim
          deadcrootc_xfer_p(m) = deadcrootc_xfer_p(m) -       &
               m_deadcrootc_xfer_to_litter_fire_p(m)* deltim
!    end if !not use_matrixcn
       end do

end subroutine CStateUpdate3

end module bgc_CNCStateUpdate3Mod
#endif
