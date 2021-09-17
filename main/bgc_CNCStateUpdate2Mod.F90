module bgc_CNCStateUpdate2Mod

use precision
use MOD_TimeInvariants, only: &
           i_met_lit,i_cel_lit,i_lig_lit ,i_cwd 
use MOD_TimeVariables, only: &
    ! decomposition pools & fluxes variables (inout)
           decomp_cpools_vr

use MOD_1D_Fluxes, only: &
           gap_mortality_to_met_c, gap_mortality_to_cel_c , &
           gap_mortality_to_lig_c, gap_mortality_to_cwdc  
  
use MOD_PFTimeVars, only: &
    ! vegetation carbon state variables (inout)
           leafc_p            , leafc_storage_p     , leafc_xfer_p     , &
           frootc_p           , frootc_storage_p    , frootc_xfer_p    , &
           livestemc_p        , livestemc_storage_p , livestemc_xfer_p , &
           deadstemc_p        , deadstemc_storage_p , deadstemc_xfer_p , &
           livecrootc_p       , livecrootc_storage_p, livecrootc_xfer_p, &
           deadcrootc_p       , deadcrootc_storage_p, deadcrootc_xfer_p, &
           gresp_storage_p    , gresp_xfer_p       

use MOD_1D_PFTFluxes, only: &
    ! vegetation carbon flux variables
           m_leafc_to_litter_p        , m_leafc_storage_to_litter_p     , m_leafc_xfer_to_litter_p     , &
           m_frootc_to_litter_p       , m_frootc_storage_to_litter_p    , m_frootc_xfer_to_litter_p    , &
           m_livestemc_to_litter_p    , m_livestemc_storage_to_litter_p , m_livestemc_xfer_to_litter_p , &
           m_deadstemc_to_litter_p    , m_deadstemc_storage_to_litter_p , m_deadstemc_xfer_to_litter_p , &
           m_livecrootc_to_litter_p   , m_livecrootc_storage_to_litter_p, m_livecrootc_xfer_to_litter_p, &
           m_deadcrootc_to_litter_p   , m_deadcrootc_storage_to_litter_p, m_deadcrootc_xfer_to_litter_p, &
           m_gresp_storage_to_litter_p, m_gresp_xfer_to_litter_p

implicit none

public CStateUpdate2

contains

subroutine CStateUpdate2 (i, ps, pe, deltim, nl_soil)

integer ,intent(in) :: i
integer ,intent(in) :: ps
integer ,intent(in) :: pe
real(r8),intent(in) :: deltim
integer ,intent(in) :: nl_soil

integer j
integer m

      ! column level carbon fluxes from gap-phase mortality
       do j = 1,nl_soil
             ! column gap mortality fluxes
!             if (.not. use_soil_matrixcn)then
          decomp_cpools_vr(j,i_met_lit,i) = &
                decomp_cpools_vr(j,i_met_lit,i) + gap_mortality_to_met_c(j,i) * deltim
          decomp_cpools_vr(j,i_cel_lit,i) = &
                decomp_cpools_vr(j,i_cel_lit,i) + gap_mortality_to_cel_c(j,i) * deltim
          decomp_cpools_vr(j,i_lig_lit,i) = &
                decomp_cpools_vr(j,i_lig_lit,i) + gap_mortality_to_lig_c(j,i) * deltim
          decomp_cpools_vr(j,i_cwd,i)     = &
                decomp_cpools_vr(j,i_cwd,i)     + gap_mortality_to_cwdc(j,i) * deltim
!             else
!                cf_soil%matrix_Cinput%V(j+(i_met_lit-1)*nlevdecomp) = &
!                  cf_soil%matrix_Cinput%V(j+(i_met_lit-1)*nlevdecomp) + gap_mortality_c_to_met_c(j) * deltim
!                cf_soil%matrix_Cinput%V(j+(i_cel_lit-1)*nlevdecomp) = &
!                  cf_soil%matrix_Cinput%V(j+(i_cel_lit-1)*nlevdecomp) + gap_mortality_c_to_cel_c(j) * deltim
!                cf_soil%matrix_Cinput%V(j+(i_lig_lit-1)*nlevdecomp) = &
!                  cf_soil%matrix_Cinput%V(j+(i_lig_lit-1)*nlevdecomp) + gap_mortality_c_to_lig_c(j) * deltim
!                cf_soil%matrix_Cinput%V(j+(i_cwd-1)*nlevdecomp) =     &
!                  cf_soil%matrix_Cinput%V(j+(i_cwd-1)*nlevdecomp)     + gap_mortality_c_to_cwdc(j) * deltim
!             end if !soil_matrix
          !print*,'gap_mortality_to_met_c',j,gap_mortality_to_met_c(j,i) * deltim
          !print*,'gap_mortality_to_cel_c',j,gap_mortality_to_cel_c(j,i) * deltim
       end do

       ! patch loop

       do m = ps, pe
          gresp_xfer_p(m) = gresp_xfer_p(m)                 &
                  - m_gresp_xfer_to_litter_p(m) * deltim
          gresp_storage_p(m) = gresp_storage_p(m)           &
                  - m_gresp_storage_to_litter_p(m) * deltim
 !         if(.not.  use_matrixcn)then
          ! patch-level carbon fluxes from gap-phase mortality
          ! displayed pools
          leafc_p(m) = leafc_p(m)                           &
                  - m_leafc_to_litter_p(m) * deltim
          frootc_p(m) = frootc_p(m)                         &
                  - m_frootc_to_litter_p(m) * deltim
          livestemc_p(m) = livestemc_p(m)                   &
                  - m_livestemc_to_litter_p(m) * deltim
          deadstemc_p(m) = deadstemc_p(m)                   &
                  - m_deadstemc_to_litter_p(m) * deltim
          livecrootc_p(m) = livecrootc_p(m)                 &
                  - m_livecrootc_to_litter_p(m) * deltim
          deadcrootc_p(m) = deadcrootc_p(m)                 &
                  - m_deadcrootc_to_litter_p(m) * deltim

          ! storage pools
          leafc_storage_p(m) = leafc_storage_p(m)           &
                  - m_leafc_storage_to_litter_p(m) * deltim
          frootc_storage_p(m) = frootc_storage_p(m)         &
                  - m_frootc_storage_to_litter_p(m) * deltim
          livestemc_storage_p(m) = livestemc_storage_p(m)   &
                  - m_livestemc_storage_to_litter_p(m) * deltim
          deadstemc_storage_p(m) = deadstemc_storage_p(m)   &
                  - m_deadstemc_storage_to_litter_p(m) * deltim
          livecrootc_storage_p(m) = livecrootc_storage_p(m) &
                  - m_livecrootc_storage_to_litter_p(m) * deltim
          deadcrootc_storage_p(m) = deadcrootc_storage_p(m) &
                  - m_deadcrootc_storage_to_litter_p(m) * deltim

          ! transfer pools
          leafc_xfer_p(m) = leafc_xfer_p(m)                 &
                  - m_leafc_xfer_to_litter_p(m) * deltim
          frootc_xfer_p(m) = frootc_xfer_p(m)               &
                  - m_frootc_xfer_to_litter_p(m) * deltim
          livestemc_xfer_p(m) = livestemc_xfer_p(m)         &
                  - m_livestemc_xfer_to_litter_p(m) * deltim
          deadstemc_xfer_p(m) = deadstemc_xfer_p(m)         &
                  - m_deadstemc_xfer_to_litter_p(m) * deltim
          livecrootc_xfer_p(m) = livecrootc_xfer_p(m)       &
                  - m_livecrootc_xfer_to_litter_p(m) * deltim
          deadcrootc_xfer_p(m) = deadcrootc_xfer_p(m)       &
                  - m_deadcrootc_xfer_to_litter_p(m) * deltim
       end do


end subroutine CStateUpdate2

end module bgc_CNCStateUpdate2Mod
