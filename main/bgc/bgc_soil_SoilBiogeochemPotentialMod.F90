#include <define.h>
#ifdef BGC
module bgc_soil_SoilBiogeochemPotentialMod

use precision
use MOD_BGCTimeInvars, only: &
    floating_cn_ratio, initial_cn_ratio, rf_decomp, receiver_pool, donor_pool, i_atm, pathfrac_decomp

use MOD_BGCTimeVars, only: &
    ! decomposition carbon & nitrogen pools
    decomp_cpools_vr, decomp_npools_vr, decomp_k, &

    ! other variables
    cn_decomp_pools

use MOD_1D_BGCFluxes, only: &
    ! decomposition fluxes variables
    pmnf_decomp, p_decomp_cpool_loss, gross_nmin_vr, &

    ! mineral N fluxes
    potential_immob_vr, phr_vr


implicit none

public SoilBiogeochemPotential

contains

   subroutine SoilBiogeochemPotential(i,nl_soil,ndecomp_pools,ndecomp_transitions)

     integer ,intent(in) :: i
     integer ,intent(in) :: nl_soil
     integer ,intent(in) :: ndecomp_pools
     integer ,intent(in) :: ndecomp_transitions

     integer j,k,l
     real(r8) immob(1:nl_soil)
     real(r8) ratio

      p_decomp_cpool_loss(:, :, i) = 0._r8
      pmnf_decomp(:, :, i) = 0._r8

      ! column loop to calculate potential decomp rates and total immobilization demand

      !! calculate c:n ratios of applicable pools
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

      ! calculate the non-nitrogen-limited fluxes
      ! these fluxes include the  "/ dt" term to put them on a
      ! per second basis, since the rate constants have been
      ! calculated on a per timestep basis.

      do k = 1, ndecomp_transitions
         do j = 1,nl_soil

            if (decomp_cpools_vr(j,donor_pool(k),i) > 0._r8 .and. &
                 decomp_k(j,donor_pool(k),i) > 0._r8 ) then
               p_decomp_cpool_loss(j,k,i) = decomp_cpools_vr(j,donor_pool(k),i) &
                    * decomp_k(j,donor_pool(k),i)  * pathfrac_decomp(j,k,i)
               if ( .not. floating_cn_ratio(receiver_pool(k)) ) then  !! not transition of cwd to litter

                  if (receiver_pool(k) /= i_atm ) then  ! not 100% respiration
                     ratio = 0._r8

                     if (decomp_npools_vr(j,donor_pool(k),i) > 0._r8) then
                        ratio = cn_decomp_pools(j,receiver_pool(k),i)/cn_decomp_pools(j,donor_pool(k),i)
                     endif

                     pmnf_decomp(j,k,i) = (p_decomp_cpool_loss(j,k,i) * (1.0_r8 - rf_decomp(j,k,i) - ratio) &
                          / cn_decomp_pools(j,receiver_pool(k),i) )

                  else   ! 100% respiration
                     pmnf_decomp(j,k,i) = - p_decomp_cpool_loss(j,k,i) / cn_decomp_pools(j,donor_pool(k),i)
                  endif

               else   ! CWD -> litter
                  pmnf_decomp(j,k,i) = 0._r8
               end if
            end if

         end do
      end do

      ! Sum up all the potential immobilization fluxes (positive pmnf flux)
      ! and all the mineralization fluxes (negative pmnf flux)
      do j = 1,nl_soil
         immob(j) = 0._r8
      end do
      do k = 1, ndecomp_transitions
         do j = 1,nl_soil
            if (pmnf_decomp(j,k,i) > 0._r8) then
               immob(j) = immob(j) + pmnf_decomp(j,k,i)
            else
               gross_nmin_vr(j,i) = gross_nmin_vr(j,i) - pmnf_decomp(j,k,i)
            end if
         end do
      end do

      do j = 1,nl_soil
         potential_immob_vr(j,i) = immob(j)
      end do

      ! Add up potential hr for methane calculations
      do j = 1,nl_soil
         phr_vr(j,i) = 0._r8
      end do
      do k = 1, ndecomp_transitions
         do j = 1,nl_soil
            phr_vr(j,i) = phr_vr(j,i) + rf_decomp(j,k,i) * p_decomp_cpool_loss(j,k,i)
         end do
      end do

   end subroutine SoilBiogeochemPotential

end module bgc_soil_SoilBiogeochemPotentialMod
#endif
