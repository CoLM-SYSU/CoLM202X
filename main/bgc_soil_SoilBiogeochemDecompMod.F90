module bgc_soil_SoilBiogeochemDecompMod

use precision
use MOD_TimeInvariants, only: &
    floating_cn_ratio, initial_cn_ratio, dnp, rf_decomp, receiver_pool, donor_pool, i_atm 

use MOD_TimeVariables, only: &
    ! decomposition carbon & nitrogen pools
    decomp_cpools_vr, decomp_npools_vr, &

    ! other variables
    cn_decomp_pools, fpi_vr

use MOD_1D_Fluxes, only: &
    ! decomposition fluxes variables
    decomp_sminn_flux_vr, decomp_hr_vr, decomp_ctransfer_vr, decomp_ntransfer_vr, &
    pmnf_decomp, p_decomp_cpool_loss, sminn_to_denit_decomp_vr, &
    net_nmin_vr, gross_nmin_vr, net_nmin, gross_nmin


implicit none

public SoilBiogeochemDecomp

contains

subroutine SoilBiogeochemDecomp(i,nl_soil,ndecomp_pools,ndecomp_transitions, dz_soi)

  integer ,intent(in) :: i
  integer ,intent(in) :: nl_soil
  integer ,intent(in) :: ndecomp_pools
  integer ,intent(in) :: ndecomp_transitions
  real(r8),intent(in) :: dz_soi(1:nl_soil)

  integer j,k,l
!   if ( .not. use_fates) then
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
!                  if (use_soil_matrixcn)then ! correct only when one transfer from each litter pool
!                     Ksoil%DM(j+nl_soil*(donor_pool(k)-1)) &
!                     = Ksoil%DM(j+nl_soil*(donor_pool(k)-1)) * fpi_vr(j)
!                  end if
#ifndef NITRIF
                  sminn_to_denit_decomp_vr(j,k,i) = 0._r8
#endif
               else
#ifndef NITRIF
                  sminn_to_denit_decomp_vr(j,k,i) = -dnp * pmnf_decomp(j,k,i)
#endif
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
#ifndef NITRIF
               sminn_to_denit_decomp_vr(j,k,i) = 0._r8
#endif
               decomp_sminn_flux_vr(j,k,i) = 0._r8
            end if

         end do
      end do
!else
!   do k = 1, ndecomp_transitions
!      do j = 1,nl_soil
!            !
!            decomp_hr_vr(j,k) = rf_decomp(j,k) * p_decomp_cpool_loss(j,k)
!            !
!            decomp_ctransfer_vr(j,k) = (1._r8 - rf_decomp(j,k)) * p_decomp_cpool_loss(j,k)
!            !
!      end do
!   end do
!end if

do j = 1,nl_soil
!   if(.not.use_fates)then
   net_nmin(i) = net_nmin(i) + net_nmin_vr(j,i) * dz_soi(j)
   gross_nmin(i) = gross_nmin(i) + gross_nmin_vr(j,i) * dz_soi(j)
   ! else
   !   net_nmin( = 0.0_r8
   !   gross_nmin( = 0.0_r8
!   endif
end do

end subroutine SoilBiogeochemDecomp

end module bgc_soil_SoilBiogeochemDecompMod
