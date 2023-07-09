#include <define.h>
#ifdef BGC
module MOD_BGC_Soil_BiogeochemDecomp

  !-----------------------------------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module caluclates the CN transfer fluxes between different soil and litter pools,
  ! which includes CN transfer fluxes (decomp_ctransfer or decomp_ntransfer), heterotrophic respiration (decomp_hr),
  ! net mineralisation and gross mineralisation. Denitrification flux will be also calculated when nitrification model 
  ! is activated. 
  !
  ! !ORIGINAL:
  ! The Community Land Model version 5.0 (CLM5.0)
  !
  ! !REVISION:
  ! Xingjie Lu, 2021, revised original CLM5 code to be compatible with CoLM code structure.

  use MOD_Precision
  use MOD_Namelist, only : DEF_USE_NITRIF
  use MOD_BGC_Vars_TimeInvariants, only: &
      floating_cn_ratio, initial_cn_ratio, dnp, rf_decomp, receiver_pool, donor_pool, i_atm 

  use MOD_BGC_Vars_TimeVariables, only: &
      ! decomposition carbon & nitrogen pools
      decomp_cpools_vr, decomp_npools_vr, &

      ! other variables
      cn_decomp_pools, fpi_vr

  use MOD_BGC_Vars_1DFluxes, only: &
      ! decomposition fluxes variables
      decomp_sminn_flux_vr, decomp_hr_vr, decomp_ctransfer_vr, decomp_ntransfer_vr, &
      pmnf_decomp, p_decomp_cpool_loss, sminn_to_denit_decomp_vr, &
      net_nmin_vr, gross_nmin_vr, net_nmin, gross_nmin


  implicit none

  public SoilBiogeochemDecomp

contains

  subroutine SoilBiogeochemDecomp(i,nl_soil,ndecomp_pools,ndecomp_transitions, dz_soi)

    integer ,intent(in) :: i                   ! patch index
    integer ,intent(in) :: nl_soil             ! number of total soil layers
    integer ,intent(in) :: ndecomp_pools       ! number of total soil & litter pools in the decompositions
    integer ,intent(in) :: ndecomp_transitions ! number of total transfers between soil and litter pools in the decomposition
    real(r8),intent(in) :: dz_soi(1:nl_soil)   ! thicknesses of each soil layer

    integer j,k,l
      ! calculate c:n ratios of applicable pools
    do l = 1, ndecomp_pools
       if ( floating_cn_ratio(l) ) then
          do j = 1,nl_soil
             if ( decomp_npools_vr(j,l,i) > 0._r8 ) then
                cn_decomp_pools(j,l,i) = decomp_cpools_vr(j,l,i) / decomp_npools_vr(j,l,i)
             end if
          end do
       else
          do j = 1,nl_soil
             cn_decomp_pools(j,l,i) = initial_cn_ratio(l)
          end do
       end if
    end do

    ! column loop to calculate actual immobilization and decomp rates, following
    ! resolution of plant/heterotroph  competition for mineral N

    ! upon return from SoilBiogeochemCompetition, the fraction of potential immobilization
    ! has been set (soilbiogeochem_state_inst%fpi_vr_col). now finish the decomp calculations.
    ! Only the immobilization steps are limited by fpi_vr (pmnf > 0)
    ! Also calculate denitrification losses as a simple proportion
    ! of mineralization flux.

    do k = 1, ndecomp_transitions
       do j = 1,nl_soil
          if (decomp_cpools_vr(j,donor_pool(k),i) > 0._r8) then
             if ( pmnf_decomp(j,k,i) > 0._r8 ) then
                p_decomp_cpool_loss(j,k,i) = p_decomp_cpool_loss(j,k,i) * fpi_vr(j,i)
                pmnf_decomp(j,k,i) = pmnf_decomp(j,k,i) * fpi_vr(j,i)
                if(.not. DEF_USE_NITRIF)then
                   sminn_to_denit_decomp_vr(j,k,i) = 0._r8
                end if
             else
                if(.not. DEF_USE_NITRIF)then
                   sminn_to_denit_decomp_vr(j,k,i) = -dnp * pmnf_decomp(j,k,i)
                end if
             end if
             decomp_hr_vr(j,k,i) = rf_decomp(j,k,i) * p_decomp_cpool_loss(j,k,i)
             decomp_ctransfer_vr(j,k,i) = (1._r8 - rf_decomp(j,k,i)) * p_decomp_cpool_loss(j,k,i)
             if (decomp_npools_vr(j,donor_pool(k),i) > 0._r8 .and. receiver_pool(k) /= i_atm) then
                decomp_ntransfer_vr(j,k,i) = p_decomp_cpool_loss(j,k,i) / cn_decomp_pools(j,donor_pool(k),i)
             else
                decomp_ntransfer_vr(j,k,i) = 0._r8
             endif
             if ( receiver_pool(k) /= 0 ) then
                decomp_sminn_flux_vr(j,k,i) = pmnf_decomp(j,k,i)
             else  ! keep sign convention negative for terminal pools
                decomp_sminn_flux_vr(j,k,i) = - pmnf_decomp(j,k,i)
             endif
             net_nmin_vr(j,i) = net_nmin_vr(j,i) - pmnf_decomp(j,k,i)
          else
             decomp_ntransfer_vr(j,k,i) = 0._r8
             if(.not. DEF_USE_NITRIF)then
                sminn_to_denit_decomp_vr(j,k,i) = 0._r8
             end if
             decomp_sminn_flux_vr(j,k,i) = 0._r8
          end if

       end do
    end do

    do j = 1,nl_soil
       net_nmin(i) = net_nmin(i) + net_nmin_vr(j,i) * dz_soi(j)
       gross_nmin(i) = gross_nmin(i) + gross_nmin_vr(j,i) * dz_soi(j)
    end do

  end subroutine SoilBiogeochemDecomp

end module MOD_BGC_Soil_BiogeochemDecomp
#endif
