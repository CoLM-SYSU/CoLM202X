#include <define.h>

module bgc_SoilBiogeochemNStateUpdate1Mod
use precision
use MOD_TimeInvariants, only: &
    receiver_pool, donor_pool, nitrif_n2o_loss_frac

use MOD_TimeVariables, only: &
    ! Mineral nitrogen pools (inout)
           sminn_vr                , smin_nh4_vr              , smin_no3_vr              , &
           ndep_prof               , nfixation_prof

use MOD_1D_Fluxes, only: &
    ! Decomposition fluxes variables (inout)
           decomp_npools_sourcesink, decomp_ntransfer_vr      , decomp_sminn_flux_vr     , sminn_to_denit_decomp_vr, &
           gross_nmin_vr           , actual_immob_nh4_vr      , actual_immob_no3_vr      , &
           sminn_to_plant_vr       , smin_nh4_to_plant_vr     , smin_no3_to_plant_vr     , supplement_to_sminn_vr, &
           sminn_to_plant_fun_vr   , sminn_to_plant_fun_nh4_vr, sminn_to_plant_fun_no3_vr, &
           sminn_to_denit_excess_vr, f_nit_vr                 , f_denit_vr               , &
           ndep_to_sminn           , ffix_to_sminn            , nfix_to_sminn            

implicit none

public SoilBiogeochemNStateUpdate1

contains

subroutine SoilBiogeochemNStateUpdate1(i,deltim,nl_soil,ndecomp_transitions)

integer ,intent(in) :: i
real(r8),intent(in) :: deltim
integer ,intent(in) :: nl_soil
integer ,intent(in) :: ndecomp_transitions

integer j,k

      do j = 1, nl_soil
#ifdef FUN
               ! N deposition and fixation (put all into NH4 pool)
               smin_nh4_vr(j,i) = smin_nh4_vr(j,i) + ndep_to_sminn(i)*deltim * ndep_prof(j,i)
               smin_nh4_vr(j,i) = smin_nh4_vr(j,i) + ffix_to_sminn(i)*deltim * nfixation_prof(j,i)
#else
#ifndef NITRIF
               ! N deposition and fixation
               sminn_vr(j,i) = sminn_vr(j,i) + ndep_to_sminn(i)*deltim * ndep_prof(j,i)
               sminn_vr(j,i) = sminn_vr(j,i) + nfix_to_sminn(i)*deltim * nfixation_prof(j,i)

#else
               ! N deposition and fixation (put all into NH4 pool)
               smin_nh4_vr(j,i) = smin_nh4_vr(j,i) + ndep_to_sminn(i)*deltim * ndep_prof(j,i)
               smin_nh4_vr(j,i) = smin_nh4_vr(j,i) + nfix_to_sminn(i)*deltim * nfixation_prof(j,i)
                       
#endif
#endif
!         print*,'after ndep and ffix',i,j,sminn_vr(j,i),ndep_to_sminn(i),deltim , ndep_prof(j,i),&
!                                                        nfix_to_sminn(i)*deltim * nfixation_prof(j,i)          
      end do

!      ! repeating N dep and fixation for crops
!      if ( use_crop )then
!         do j = 1, nl_soil
!
!            ! column loop
!            do fc = 1,num_soilc
!               c = filter_soilc(fc)
!               if (.not. use_nitrif_denitrif) then
!
!                  ! N deposition and fixation
!                  sminn_vr(j) = sminn_vr(j) &
!                       + fert_to_sminn(*deltim * ndep_prof(c,j)
!                  sminn_vr(j) = sminn_vr(j) &
!                       + soyfixn_to_sminn(*deltim * nfixation_prof(c,j)
!
!               else
!
!                  ! N deposition and fixation (put all into NH4 pool)
!                  smin_nh4_vr(j) = smin_nh4_vr(j) &
!                       + fert_to_sminn(*deltim * ndep_prof(c,j)
!                  smin_nh4_vr(j) = smin_nh4_vr(j) &
!                       + soyfixn_to_sminn(*deltim * nfixation_prof(c,j)
!
!               end if
!            end do
!         end do
!      end if

      ! decomposition fluxes
!   if (.not. use_soil_matrixcn) then
      do k = 1, ndecomp_transitions
         do j = 1, nl_soil

               decomp_npools_sourcesink(j,donor_pool(k),i) = &
                    decomp_npools_sourcesink(j,donor_pool(k),i) - &
                    decomp_ntransfer_vr(j,k,i) * deltim
         end do
      end do


      do k = 1, ndecomp_transitions
         if ( receiver_pool(k) /= 0 ) then  ! skip terminal transitions
            do j = 1, nl_soil

                  decomp_npools_sourcesink(j,receiver_pool(k),i) = &
                       decomp_npools_sourcesink(j,receiver_pool(k),i) + &
                       (decomp_ntransfer_vr(j,k,i) + &
                        decomp_sminn_flux_vr(j,k,i)) * deltim
            end do
         else  ! terminal transitions
            do j = 1, nl_soil
                  decomp_npools_sourcesink(j,donor_pool(k),i) = &
                       decomp_npools_sourcesink(j,donor_pool(k),i) - &
                       decomp_sminn_flux_vr(j,k,i) * deltim
            end do
         end if
      end do
!  end if  ! 

#ifndef NITRIF

         !--------------------------------------------------------
         !-------------    NITRIF_DENITRIF OFF -------------------
         !--------------------------------------------------------

         ! immobilization/mineralization in litter-to-SOM and SOM-to-SOM fluxes and denitrification fluxes
         do k = 1, ndecomp_transitions
            if ( receiver_pool(k) /= 0 ) then  ! skip terminal transitions
               do j = 1, nl_soil
                     sminn_vr(j,i)  = sminn_vr(j,i) - &
                          (sminn_to_denit_decomp_vr(j,k,i) + &
                          decomp_sminn_flux_vr(j,k,i))* deltim
               end do
            else
               do j = 1, nl_soil
                     sminn_vr(j,i)  = sminn_vr(j,i) - &
                          sminn_to_denit_decomp_vr(j,k,i)* deltim

                     sminn_vr(j,i)  = sminn_vr(j,i) + &
                          decomp_sminn_flux_vr(j,k,i)* deltim

               end do
            endif
         end do

         do j = 1, nl_soil
               ! "bulk denitrification"
               sminn_vr(j,i) = sminn_vr(j,i) - sminn_to_denit_excess_vr(j,i) * deltim

               ! total plant uptake from mineral N
#if !defined(FUN)
                  sminn_vr(j,i) = sminn_vr(j,i) - sminn_to_plant_vr(j,i)*deltim
#else
                  sminn_vr(j,i) = sminn_vr(j,i) - sminn_to_plant_fun_vr(j,i)*deltim
#endif
               ! flux that prevents N limitation (when Carbon_only is set)
               sminn_vr(j,i) = sminn_vr(j,i) + supplement_to_sminn_vr(j,i)*deltim
         end do

#else

         !--------------------------------------------------------
         !-------------    NITRIF_DENITRIF ON --------------------
         !--------------------------------------------------------

         do j = 1, nl_soil

               ! mineralization fluxes (divert a fraction of this stream to nitrification flux, add the rest to NH4 pool)
               smin_nh4_vr(j,i) = smin_nh4_vr(j,i) + gross_nmin_vr(j,i)*deltim

               ! immobilization fluxes
               smin_nh4_vr(j,i) = smin_nh4_vr(j,i) - actual_immob_nh4_vr(j,i)*deltim

               smin_no3_vr(j,i) = smin_no3_vr(j,i) - actual_immob_no3_vr(j,i)*deltim

               ! plant uptake fluxes
#if !defined(FUN)
                  smin_nh4_vr(j,i) = smin_nh4_vr(j,i) - smin_nh4_to_plant_vr(j,i)*deltim

                  smin_no3_vr(j,i) = smin_no3_vr(j,i) - smin_no3_to_plant_vr(j,i)*deltim
#else
                  smin_nh4_vr(j,i) = smin_nh4_vr(j,i) -  sminn_to_plant_fun_nh4_vr(j,i)*deltim

                  smin_no3_vr(j,i) = smin_no3_vr(j,i) -  sminn_to_plant_fun_no3_vr(j,i)*deltim
#endif
             

               ! Account for nitrification fluxes
               smin_nh4_vr(j,i) = smin_nh4_vr(j,i) - f_nit_vr(j,i) * deltim

               smin_no3_vr(j,i) = smin_no3_vr(j,i) + f_nit_vr(j,i) * deltim &
                    * (1._r8 - nitrif_n2o_loss_frac)

               ! Account for denitrification fluxes
               smin_no3_vr(j,i) = smin_no3_vr(j,i) - f_denit_vr(j,i) * deltim

               ! flux that prevents N limitation (when Carbon_only is set; put all into NH4)
               smin_nh4_vr(j,i) = smin_nh4_vr(j,i) + supplement_to_sminn_vr(j,i)*deltim

               ! update diagnostic total
               sminn_vr(j,i) = smin_nh4_vr(j,i) + smin_no3_vr(j,i)
               
         end do
              
#endif

end subroutine SoilBiogeochemNStateUpdate1

end module bgc_SoilBiogeochemNStateUpdate1Mod
