#include <define.h>
#ifdef BGC
module bgc_soil_SoilBiogeochemCompetitionMod

use precision
use MOD_1D_BGCFluxes, only: &
    pot_f_nit_vr, potential_immob_vr, sminn_to_plant_vr, sminn_to_denit_excess_vr, plant_ndemand, &
    actual_immob_vr, sminn_to_plant, pot_f_nit_vr, actual_immob_nh4_vr, f_nit_vr, &
    smin_nh4_to_plant_vr, pot_f_denit_vr, actual_immob_no3_vr, f_denit_vr, smin_no3_to_plant_vr, &
    n2_n2o_ratio_denit_vr, f_n2o_nit_vr, f_n2o_denit_vr
use MOD_BGCTimeVars, only: &
    sminn_vr, smin_no3_vr, smin_nh4_vr, nfixation_prof, fpi_vr, fpi, fpg
use MOD_BGCTimeInVars,only: &
    bdnr, compet_plant_no3, compet_plant_nh4, compet_decomp_no3, compet_decomp_nh4, compet_denit, compet_nit, &
    nitrif_n2o_loss_frac

implicit none

public SoilBiogeochemCompetition

contains

subroutine SoilBiogeochemCompetition(i,deltim,nl_soil,dz_soi)

  integer ,intent(in) :: i
  real(r8),intent(in) :: deltim
  integer ,intent(in) :: nl_soil
  real(r8),intent(in) :: dz_soi(1:nl_soil)


    integer  :: p,l,pi,j                                            ! indices
!    logical :: local_use_fun                                          ! local version of use_fun
    real(r8) :: fpi_no3_vr(1:nl_soil)      ! fraction of potential immobilization supplied by no3(no units)
    real(r8) :: fpi_nh4_vr(1:nl_soil)      ! fraction of potential immobilization supplied by nh4 (no units)
    real(r8) :: sum_nh4_demand(1:nl_soil)
    real(r8) :: sum_nh4_demand_scaled(1:nl_soil)
    real(r8) :: sum_no3_demand(1:nl_soil)
    real(r8) :: sum_no3_demand_scaled(1:nl_soil)
    real(r8) :: sum_ndemand_vr( 1:nl_soil) !total column N demand (gN/m3/s) at a given level
    real(r8) :: nuptake_prof( 1:nl_soil)
    real(r8) :: sminn_tot
    integer  :: nlimit(1:nl_soil)          !flag for N limitation
    integer  :: nlimit_no3(1:nl_soil)      !flag for NO3 limitation
    integer  :: nlimit_nh4(1:nl_soil)      !flag for NH4 limitation
    real(r8) :: residual_sminn_vr( 1:nl_soil)
    real(r8) :: residual_sminn
    real(r8) :: residual_smin_nh4_vr( 1:nl_soil)
    real(r8) :: residual_smin_no3_vr( 1:nl_soil)
    real(r8) :: residual_smin_nh4
    real(r8) :: residual_smin_no3
    real(r8) :: residual_plant_ndemand
    real(r8) :: sminn_to_plant_new
    real(r8) :: actual_immob = 0
    real(r8) :: potential_immob = 0
    !-----------------------------------------------------------------------

    sminn_to_plant_new  =  0._r8

!      local_use_fun = use_fun

#ifndef NITRIF

         ! init sminn_tot
         sminn_tot = 0.

         do j = 1, nl_soil
            sminn_tot = sminn_tot + sminn_vr(j,i) * dz_soi(j)
         end do

         do j = 1, nl_soil
            if (sminn_tot  >  0.) then
               nuptake_prof(j) = sminn_vr(j,i) / sminn_tot
            else
               nuptake_prof(j) = nfixation_prof(j,i)
            endif
         end do

         do j = 1, nl_soil
            sum_ndemand_vr(j) = plant_ndemand(i) * nuptake_prof(j) + potential_immob_vr(j,i)
         end do

         do j = 1, nl_soil
!               l = col%landunit(c)
!            if(i .eq. 79738)print*,'here1',j,sum_ndemand_vr(j)*deltim, sminn_vr(j,i)
            if (sum_ndemand_vr(j)*deltim < sminn_vr(j,i)) then

               ! N availability is not limiting immobilization or plant
               ! uptake, and both can proceed at their potential rates
               nlimit(j) = 0
               fpi_vr(j,i) = 1.0_r8
               actual_immob_vr(j,i) = potential_immob_vr(j,i)
               sminn_to_plant_vr(j,i) = plant_ndemand(i) * nuptake_prof(j)
!               else if ( cnallocate_carbon_only()) then !.or. &
!                  ! this code block controls the addition of N to sminn pool
!                  ! to eliminate any N limitation, when Carbon_Only is set.  This lets the
!                  ! model behave essentially as a carbon-only model, but with the
!                  ! benefit of keeping track of the N additions needed to
!                  ! eliminate N limitations, so there is still a diagnostic quantity
!                  ! that describes the degree of N limitation at steady-state.
!
!                  nlimit(c,j) = 1
!                  fpi_vr(c,j) = 1.0_r8
!                  actual_immob_vr(c,j) = potential_immob_vr(c,j)
!                  sminn_to_plant_vr(c,j) =  plant_ndemand(c) * nuptake_prof(c,j)
!                  supplement_to_sminn_vr(c,j) = sum_ndemand_vr(c,j) - (sminn_vr(c,j)/deltim)
            else
                  ! N availability can not satisfy the sum of immobilization and
                  ! plant growth demands, so these two demands compete for available
                  ! soil mineral N resource.

               nlimit(j) = 1
               if (sum_ndemand_vr(j) > 0.0_r8) then
                  actual_immob_vr(j,i) = (sminn_vr(j,i)/deltim)*(potential_immob_vr(j,i) / sum_ndemand_vr(j))
               else
                  actual_immob_vr(j,i) = 0.0_r8
               end if

               if (potential_immob_vr(j,i) > 0.0_r8) then
                  fpi_vr(j,i) = actual_immob_vr(j,i) / potential_immob_vr(j,i)
               else
                  fpi_vr(j,i) = 0.0_r8
               end if

               sminn_to_plant_vr(j,i) = (sminn_vr(j,i)/deltim) - actual_immob_vr(j,i)
            end if
!            if(i .eq. 79738)print*,'fpi_vr',i,j,fpi_vr(j,i),actual_immob_vr(j,i),potential_immob_vr(j,i)
!            if(i .eq. 79738)print*,'sminn_to_plant_vr',i,j,sum_ndemand_vr(j)*deltim,sminn_vr(j,i),&
!                                   sminn_to_plant_vr(j,i),plant_ndemand(i),nuptake_prof(j)
         end do

#ifdef FUN
            call t_startf( 'CNFUN' )
            call CNFUN(bounds,num_soilc,filter_soilc,num_soilp,filter_soilp,waterstatebulk_inst, &
                      waterfluxbulk_inst,temperature_inst,soilstate_inst,cnveg_state_inst,cnveg_carbonstate_inst,&
                      cnveg_carbonflux_inst,cnveg_nitrogenstate_inst,cnveg_nitrogenflux_inst                ,&
                      soilbiogeochem_nitrogenflux_inst,soilbiogeochem_carbonflux_inst,canopystate_inst,      &
                      soilbiogeochem_nitrogenstate_inst)
            call p2c(bounds, nl_soil, &
                      cnveg_nitrogenflux_inst%sminn_to_plant_fun_vr_patch(bounds%begp:bounds%endp,1:nl_soil),&
                      soilbiogeochem_nitrogenflux_inst%sminn_to_plant_fun_vr_col(bounds%begc:bounds%endc,1:nl_soil), &
                      'unity')
            call t_stopf( 'CNFUN' )
#endif

         ! sum up N fluxes to plant
         do j = 1, nl_soil
            sminn_to_plant(i) = sminn_to_plant(i) + sminn_to_plant_vr(j,i) * dz_soi(j)
#ifdef FUN
                  if (sminn_to_plant_fun_vr(c,j).gt.sminn_to_plant_vr(c,j)) then
                      sminn_to_plant_fun_vr(c,j)  = sminn_to_plant_vr(c,j)
                  end if
#endif
         end do

         ! give plants a second pass to see if there is any mineral N left over with which to satisfy residual N demand.
         residual_sminn = 0._r8

         ! sum up total N left over after initial plant and immobilization fluxes
         residual_plant_ndemand = plant_ndemand(i) - sminn_to_plant(i)

         do j = 1, nl_soil
            if (residual_plant_ndemand  >  0._r8 ) then
               if (nlimit(j) .eq. 0) then
                  residual_sminn_vr(j) = max(sminn_vr(j,i) - (actual_immob_vr(j,i) + sminn_to_plant_vr(j,i) ) * deltim, 0._r8)
                  residual_sminn       = residual_sminn + residual_sminn_vr(j) * dz_soi(j)
               else
                  residual_sminn_vr(j) = 0._r8
               endif
            endif
         end do

         ! distribute residual N to plants
         do j = 1, nl_soil
            if ( residual_plant_ndemand  >  0._r8 .and. residual_sminn  >  0._r8 .and. nlimit(j) .eq. 0) then
               sminn_to_plant_vr(j,i) = sminn_to_plant_vr(j,i) + residual_sminn_vr(j) * &
                    min(( residual_plant_ndemand *  deltim ) / residual_sminn, 1._r8) / deltim
            endif
         end do

         ! re-sum up N fluxes to plant
         sminn_to_plant(i) = 0._r8
         do j = 1, nl_soil
            sminn_to_plant(i) = sminn_to_plant(i) + sminn_to_plant_vr(j,i) * dz_soi(j)
#ifndef FUN
            sum_ndemand_vr(j) = potential_immob_vr(j,i) + sminn_to_plant_vr(j,i)
#else
            sminn_to_plant_new(i)  = sminn_to_plant_new(i)   + sminn_to_plant_fun_vr(c,j) * dz_soi(j)
            sum_ndemand_vr(c,j)    = potential_immob_vr(c,j) + sminn_to_plant_fun_vr(c,j)
#endif
         end do

         ! under conditions of excess N, some proportion is assumed to
         ! be lost to denitrification, in addition to the constant
         ! proportion lost in the decomposition pathways
         do j = 1, nl_soil
#ifndef FUN
            if ((sminn_to_plant_vr(j,i) + actual_immob_vr(j,i))*deltim < sminn_vr(j,i)) then
               sminn_to_denit_excess_vr(j,i) = max(bdnr*deltim/86400._r8*((sminn_vr(j,i)/deltim) - sum_ndemand_vr(j)),0._r8)
            else
               sminn_to_denit_excess_vr(j,i) = 0._r8
            endif
#else
                  if ((sminn_to_plant_fun_vr(c,j)  + actual_immob_vr(c,j))*deltim < sminn_vr(c,j))  then
                     sminn_to_denit_excess_vr(c,j) = max(bdnr*deltim/86400._r8*((sminn_vr(c,j)/deltim) - sum_ndemand_vr(c,j)),0._r8)
                  else
                     sminn_to_denit_excess_vr(c,j) = 0._r8
                  endif
#endif
         end do

         ! sum up N fluxes to immobilization
         do j = 1, nl_soil
            actual_immob    = actual_immob    + actual_immob_vr(j,i) * dz_soi(j)
            potential_immob = potential_immob + potential_immob_vr(j,i) * dz_soi(j)
         end do

         ! calculate the fraction of potential growth that can be
         ! acheived with the N available to plants      
         if (plant_ndemand(i) > 0.0_r8) then
#ifndef FUN
            fpg(i) = sminn_to_plant(i) / plant_ndemand(i)
#else
            fpg(c) = sminn_to_plant_new(c) / plant_ndemand(c)
#endif
         else
            fpg(i) = 1.0_r8
         end if

         ! calculate the fraction of immobilization realized (for diagnostic purposes)
         if (potential_immob > 0.0_r8) then
!            if(i .eq. 79738)print*,'potential_immob',i,potential_immob
            fpi(i) = actual_immob / potential_immob
         else
            fpi(i) = 1.0_r8
         end if

#else 
! ifdef NITRIF

         ! column loops to resolve plant/heterotroph/nitrifier/denitrifier competition for mineral N
         !read constants from external netcdf file
!         compet_plant_no3  = params_inst%compet_plant_no3
!         compet_plant_nh4  = params_inst%compet_plant_nh4
!         compet_decomp_no3 = params_inst%compet_decomp_no3
!         compet_decomp_nh4 = params_inst%compet_decomp_nh4
!         compet_denit      = params_inst%compet_denit
!         compet_nit        = params_inst%compet_nit

         ! init total mineral N pools
         sminn_tot = 0.

         ! sum up total mineral N pools
         do j = 1, nl_soil
            sminn_tot = sminn_tot + (smin_no3_vr(j,i) + smin_nh4_vr(j,i)) * dz_soi(j)
         end do

         ! define N uptake profile for initial vertical distribution of plant N uptake, assuming plant seeks N from where it is most abundant
         do j = 1, nl_soil
            if (sminn_tot  >  0.) then
               nuptake_prof(j) = sminn_vr(j,i) / sminn_tot
            else
               nuptake_prof(j) = nfixation_prof(j,i)
            endif
         end do

         ! main column/vertical loop
         do j = 1, nl_soil  
!               l = col%landunit(c)

               !  first compete for nh4
               sum_nh4_demand(j) = plant_ndemand(i) * nuptake_prof(j) + potential_immob_vr(j,i) + pot_f_nit_vr(j,i)
               sum_nh4_demand_scaled(j) = plant_ndemand(i)* nuptake_prof(j) * compet_plant_nh4 + &
                    potential_immob_vr(j,i)*compet_decomp_nh4 + pot_f_nit_vr(j,i)*compet_nit

               if (sum_nh4_demand(j)*deltim < smin_nh4_vr(j,i)) then

                  ! NH4 availability is not limiting immobilization or plant
                  ! uptake, and all can proceed at their potential rates
                  nlimit_nh4(j) = 0
                  fpi_nh4_vr(j) = 1.0_r8
                  actual_immob_nh4_vr(j,i) = potential_immob_vr(j,i)
                  !RF added new term. 

                  f_nit_vr(j,i) = pot_f_nit_vr(j,i)
                  
#ifndef FUN
                  smin_nh4_to_plant_vr(j,i) = plant_ndemand(i) * nuptake_prof(j)
#else
                  smin_nh4_to_plant_vr(c,j) = smin_nh4_vr(c,j)/deltim - actual_immob_nh4_vr(c,j) - f_nit_vr(c,j)
#endif

               else

                  ! NH4 availability can not satisfy the sum of immobilization, nitrification, and
                  ! plant growth demands, so these three demands compete for available
                  ! soil mineral NH4 resource.
                  nlimit_nh4(j) = 1
                  if (sum_nh4_demand(j) > 0.0_r8) then
                  ! RF microbes compete based on the hypothesised plant demand. 
                     actual_immob_nh4_vr(j,i) = min((smin_nh4_vr(j,i)/deltim)*(potential_immob_vr(j,i)* &
                          compet_decomp_nh4 / sum_nh4_demand_scaled(j)), potential_immob_vr(j,i))

                     f_nit_vr(j,i) =  min((smin_nh4_vr(j,i)/deltim)*(pot_f_nit_vr(j,i)*compet_nit / &
                          sum_nh4_demand_scaled(j)), pot_f_nit_vr(j,i))
                                                 
#ifndef FUN
                         smin_nh4_to_plant_vr(j,i) = min((smin_nh4_vr(j,i)/deltim)*(plant_ndemand(i)* &
                          nuptake_prof(j)*compet_plant_nh4 / sum_nh4_demand_scaled(j)), plant_ndemand(i)*nuptake_prof(j))
                          
#else
                        ! RF added new term. send rest of N to plant - which decides whether it should pay or not? 
                        smin_nh4_to_plant_vr(c,j) = smin_nh4_vr(c,j)/deltim - actual_immob_nh4_vr(c,j) - f_nit_vr(c,j)
#endif
                    
                  else
                     actual_immob_nh4_vr(j,i) = 0.0_r8
                     smin_nh4_to_plant_vr(j,i) = 0.0_r8
                     f_nit_vr(j,i) = 0.0_r8
                  end if

                  if (potential_immob_vr(j,i) > 0.0_r8) then
                     fpi_nh4_vr(j) = actual_immob_nh4_vr(j,i) / potential_immob_vr(j,i)
                  else
                     fpi_nh4_vr(j) = 0.0_r8
                  end if

               end if
          
              
              
#ifndef FUN
                   sum_no3_demand(j) = (plant_ndemand(i)*nuptake_prof(j)-smin_nh4_to_plant_vr(j,i)) + &
                  (potential_immob_vr(j,i)-actual_immob_nh4_vr(j,i)) + pot_f_denit_vr(j,i)
                   sum_no3_demand_scaled(j) = (plant_ndemand(i)*nuptake_prof(j) &
                                                 -smin_nh4_to_plant_vr(j,i))*compet_plant_no3 + &
                  (potential_immob_vr(j,i)-actual_immob_nh4_vr(j,i))*compet_decomp_no3 + pot_f_denit_vr(j,i)*compet_denit
#else
                  sum_no3_demand(c,j) = plant_ndemand(c)*nuptake_prof(c,j) + &
                  (potential_immob_vr(c,j)-actual_immob_nh4_vr(c,j)) + pot_f_denit_vr(c,j)
                   sum_no3_demand_scaled(c,j) = (plant_ndemand(c)*nuptake_prof(c,j))*compet_plant_no3 + &
                  (potential_immob_vr(c,j)-actual_immob_nh4_vr(c,j))*compet_decomp_no3 + pot_f_denit_vr(c,j)*compet_denit
#endif
                  
          

               if (sum_no3_demand(j)*deltim < smin_no3_vr(j,i)) then

                  ! NO3 availability is not limiting immobilization or plant
                  ! uptake, and all can proceed at their potential rates
                  nlimit_no3(j) = 0
                  fpi_no3_vr(j) = 1.0_r8 -  fpi_nh4_vr(j)
                  actual_immob_no3_vr(j,i) = (potential_immob_vr(j,i)-actual_immob_nh4_vr(j,i))

                  f_denit_vr(j,i) = pot_f_denit_vr(j,i)

#ifndef FUN
                     smin_no3_to_plant_vr(j,i) = (plant_ndemand(i)*nuptake_prof(j)-smin_nh4_to_plant_vr(j,i))
#else
                     ! This restricts the N uptake of a single layer to the value determined from the total demands and the 
                     ! hypothetical uptake profile above. Which is a strange thing to do, since that is independent of FUN
                     ! do we need this at all? 
                     smin_no3_to_plant_vr(c,j) = plant_ndemand(c)*nuptake_prof(c,j)
                     ! RF added new term. send rest of N to plant - which decides whether it should pay or not? 
!                     if ( local_use_fun ) then
                     smin_no3_to_plant_vr(c,j) = smin_no3_vr(c,j)/deltim - actual_immob_no3_vr(c,j) - f_denit_vr(c,j)
!                     end if
#endif
                
               else 

                  ! NO3 availability can not satisfy the sum of immobilization, denitrification, and
                  ! plant growth demands, so these three demands compete for available
                  ! soil mineral NO3 resource.
                  nlimit_no3(j) = 1
                                  
                  if (sum_no3_demand(j) > 0.0_r8) then
#ifndef FUN
                        actual_immob_no3_vr(j,i) = min((smin_no3_vr(j,i)/deltim)*((potential_immob_vr(j,i)- &
                        actual_immob_nh4_vr(j,i))*compet_decomp_no3 / sum_no3_demand_scaled(j)), &
                                  potential_immob_vr(j,i)-actual_immob_nh4_vr(j,i))
        
                        smin_no3_to_plant_vr(j,i) = min((smin_no3_vr(j,i)/deltim)*((plant_ndemand(i)* &
                                  nuptake_prof(j)-smin_nh4_to_plant_vr(j,i))*compet_plant_no3 / sum_no3_demand_scaled(j)), &
                                  plant_ndemand(i)*nuptake_prof(j)-smin_nh4_to_plant_vr(j,i))
        
                        f_denit_vr(j,i) = min((smin_no3_vr(j,i)/deltim)*(pot_f_denit_vr(j,i)*compet_denit / &
                                  sum_no3_demand_scaled(j)), pot_f_denit_vr(j,i))
#else
                        actual_immob_no3_vr(c,j) = min((smin_no3_vr(c,j)/deltim)*((potential_immob_vr(c,j)- &
                        actual_immob_nh4_vr(c,j))*compet_decomp_no3 / sum_no3_demand_scaled(c,j)), &
                                  potential_immob_vr(c,j)-actual_immob_nh4_vr(c,j))

                        f_denit_vr(c,j) = min((smin_no3_vr(c,j)/deltim)*(pot_f_denit_vr(c,j)*compet_denit / &
                        sum_no3_demand_scaled(c,j)), pot_f_denit_vr(c,j))
        
                        smin_no3_to_plant_vr(c,j) = (smin_no3_vr(c,j)/deltim)*((plant_ndemand(c)* &
                                  nuptake_prof(c,j)-smin_nh4_to_plant_vr(c,j))*compet_plant_no3 / sum_no3_demand_scaled(c,j))
                                  
                        ! RF added new term. send rest of N to plant - which decides whether it should pay or not? 
                        smin_no3_to_plant_vr(c,j) = (smin_no3_vr(c,j) / deltim) - actual_immob_no3_vr(c,j) - f_denit_vr(c,j)
                        
  
#endif

                  else ! no no3 demand. no uptake fluxes.
                     actual_immob_no3_vr(j,i) = 0.0_r8
                     smin_no3_to_plant_vr(j,i) = 0.0_r8
                     f_denit_vr(j,i) = 0.0_r8

                  end if !any no3 demand?
                  
                  
                  

                  if (potential_immob_vr(j,i) > 0.0_r8) then
                     fpi_no3_vr(j) = actual_immob_no3_vr(j,i) / potential_immob_vr(j,i)
                  else
                     fpi_no3_vr(j) = 0.0_r8
                  end if

               end if

               
                    

               ! n2o emissions: n2o from nitr is const fraction, n2o from denitr is calculated in nitrif_denitrif
               f_n2o_nit_vr(j,i) = f_nit_vr(j,i) * nitrif_n2o_loss_frac
               f_n2o_denit_vr(j,i) = f_denit_vr(j,i) / (1._r8 + n2_n2o_ratio_denit_vr(j,i))


               ! this code block controls the addition of N to sminn pool
               ! to eliminate any N limitation, when Carbon_Only is set.  This lets the
               ! model behave essentially as a carbon-only model, but with the
               ! benefit of keeping track of the N additions needed to
               ! eliminate N limitations, so there is still a diagnostic quantity
               ! that describes the degree of N limitation at steady-state.

               ! sum up no3 and nh4 fluxes
               fpi_vr(j,i) = fpi_no3_vr(j) + fpi_nh4_vr(j)
               sminn_to_plant_vr(j,i) = smin_no3_to_plant_vr(j,i) + smin_nh4_to_plant_vr(j,i)
               actual_immob_vr(j,i) = actual_immob_no3_vr(j,i) + actual_immob_nh4_vr(j,i)
         end do

#ifdef FUN
            call CNFUN(bounds,num_soilc,filter_soilc,num_soilp,filter_soilp,waterstatebulk_inst,&
                      waterfluxbulk_inst,temperature_inst,soilstate_inst,cnveg_state_inst,cnveg_carbonstate_inst,&
                      cnveg_carbonflux_inst,cnveg_nitrogenstate_inst,cnveg_nitrogenflux_inst                ,&
                      soilbiogeochem_nitrogenflux_inst,soilbiogeochem_carbonflux_inst,canopystate_inst,      &
                      soilbiogeochem_nitrogenstate_inst)
                      
            ! sminn_to_plant_fun is output of actual N uptake from FUN
            call p2c(bounds,nl_soil, &
                       cnveg_nitrogenflux_inst%sminn_to_plant_fun_no3_vr_patch(bounds%begp:bounds%endp,1:nl_soil),&
                       soilbiogeochem_nitrogenflux_inst%sminn_to_plant_fun_no3_vr_col(bounds%begc:bounds%endc,1:nl_soil),&
                       'unity')

            call p2c(bounds,nl_soil, &
                       cnveg_nitrogenflux_inst%sminn_to_plant_fun_nh4_vr_patch(bounds%begp:bounds%endp,1:nl_soil),&
                       soilbiogeochem_nitrogenflux_inst%sminn_to_plant_fun_nh4_vr_col(bounds%begc:bounds%endc,1:nl_soil),&
                       'unity')
#endif


#ifndef FUN
            ! sum up N fluxes to plant after initial competition
            sminn_to_plant(i) = 0._r8
            do j = 1, nl_soil  
               sminn_to_plant(i) = sminn_to_plant(i) + sminn_to_plant_vr(j,i) * dz_soi(j)
            end do
#else
            do fc=1,num_soilc
               c = filter_soilc(fc)
               ! sum up N fluxes to plant after initial competition
               sminn_to_plant(c) = 0._r8 !this isn't use in fun. 
               do j = 1, nl_soil
                  if ((sminn_to_plant_fun_no3_vr(c,j)-smin_no3_to_plant_vr(c,j)).gt.0.0000000000001_r8) then
                      write(iulog,*) 'problem with limitations on no3 uptake', &
                                 sminn_to_plant_fun_no3_vr(c,j),smin_no3_to_plant_vr(c,j)
                      call endrun("too much NO3 uptake predicted by FUN")
                  end if
!KO                  if ((sminn_to_plant_fun_nh4_vr(c,j)-smin_nh4_to_plant_vr(c,j)).gt.0.0000000000001_r8) then
!KO
                  if ((sminn_to_plant_fun_nh4_vr(c,j)-smin_nh4_to_plant_vr(c,j)).gt.0.0000001_r8) then
!KO
                      write(iulog,*) 'problem with limitations on nh4 uptake', &
                                  sminn_to_plant_fun_nh4_vr(c,j),smin_nh4_to_plant_vr(c,j)
                      call endrun("too much NH4 uptake predicted by FUN")
                  end if
               end do
            end do

#endif

#ifndef FUN
            ! give plants a second pass to see if there is any mineral N left over with which to satisfy residual N demand.
            ! first take frm nh4 pool; then take from no3 pool
            residual_plant_ndemand = plant_ndemand(i) - sminn_to_plant(i)
            residual_smin_nh4 = 0._r8
            do j = 1, nl_soil  
               if (residual_plant_ndemand  >  0._r8 ) then
                  if (nlimit_nh4(j) .eq. 0) then
                     residual_smin_nh4_vr(j) = max(smin_nh4_vr(j,i) - (actual_immob_nh4_vr(j,i) + &
                                                 smin_nh4_to_plant_vr(j,i) + f_nit_vr(j,i) ) * deltim, 0._r8)

                     residual_smin_nh4 = residual_smin_nh4 + residual_smin_nh4_vr(j) * dz_soi(j)
                  else
                     residual_smin_nh4_vr(j)  = 0._r8
                  endif
   
                  if ( residual_smin_nh4 > 0._r8 .and. nlimit_nh4(j) .eq. 0 ) then
                     smin_nh4_to_plant_vr(j,i) = smin_nh4_to_plant_vr(j,i) + residual_smin_nh4_vr(j) * &
                          min(( residual_plant_ndemand *  deltim ) / residual_smin_nh4, 1._r8) / deltim
                  endif
               end if
            end do

            ! re-sum up N fluxes to plant after second pass for nh4
            sminn_to_plant(i) = 0._r8
            do j = 1, nl_soil
               sminn_to_plant_vr(j,i) = smin_nh4_to_plant_vr(j,i) + smin_no3_to_plant_vr(j,i)
               sminn_to_plant(i) = sminn_to_plant(i) + (sminn_to_plant_vr(j,i)) * dz_soi(j)
            end do

            !
            ! and now do second pass for no3
            residual_plant_ndemand = plant_ndemand(i) - sminn_to_plant(i)
            residual_smin_no3 = 0._r8

            do j = 1, nl_soil
               if (residual_plant_ndemand > 0._r8 ) then
                  if (nlimit_no3(j) .eq. 0) then
                     residual_smin_no3_vr(j) = max(smin_no3_vr(j,i) - (actual_immob_no3_vr(j,i) + &
                                                smin_no3_to_plant_vr(j,i) + f_denit_vr(j,i) ) * deltim, 0._r8)
                     residual_smin_no3 = residual_smin_no3 + residual_smin_no3_vr(j) * dz_soi(j)
                  else
                     residual_smin_no3_vr(j)  = 0._r8
                  endif

                  if ( residual_smin_no3 > 0._r8 .and. nlimit_no3(j) .eq. 0) then
                     smin_no3_to_plant_vr(j,i) = smin_no3_to_plant_vr(j,i) + residual_smin_no3_vr(j) * &
                          min(( residual_plant_ndemand *  deltim ) / residual_smin_no3, 1._r8) / deltim
                  endif
               endif
            end do

            ! re-sum up N fluxes to plant after second passes of both no3 and nh4
            sminn_to_plant(i) = 0._r8
            do j = 1, nl_soil
               sminn_to_plant_vr(j,i) = smin_nh4_to_plant_vr(j,i) + smin_no3_to_plant_vr(j,i)
               sminn_to_plant(i) = sminn_to_plant(i) + (sminn_to_plant_vr(j,i)) * dz_soi(j)
!               if(i .eq. 79738)print*,'sminn_to_plant',i,j,sminn_to_plant(i),(sminn_to_plant_vr(j,i)) * dz_soi(j)
            end do
   
#else
         !calculate maximum N available to plants. 
            do fc=1,num_soilc
               c = filter_soilc(fc)
               sminn_to_plant(c) = 0._r8
            end do
            do j = 1, nl_soil
               do fc=1,num_soilc
                  c = filter_soilc(fc)
                  sminn_to_plant_vr(c,j) = smin_nh4_to_plant_vr(c,j) + smin_no3_to_plant_vr(c,j)
                  sminn_to_plant(c) = sminn_to_plant(c) + (sminn_to_plant_vr(c,j)) * dz_soi(j)
               end do
            end do
   

             ! add up fun fluxes from SMINN to plant. 
             do j = 1, nl_soil
                do fc=1,num_soilc
                   c = filter_soilc(fc)
                   sminn_to_plant_new(c)  = sminn_to_plant_new(c) + &
                             (sminn_to_plant_fun_no3_vr(c,j) + sminn_to_plant_fun_nh4_vr(c,j)) * dz_soi(j)
                      
                end do
             end do
                             
                              
#endif
         ! sum up N fluxes to immobilization
         actual_immob = 0._r8
         potential_immob = 0._r8
         do j = 1, nl_soil  
            actual_immob = actual_immob + actual_immob_vr(j,i) * dz_soi(j)
            potential_immob = potential_immob + potential_immob_vr(j,i) * dz_soi(j)
         end do
        
        
     
       
            ! calculate the fraction of potential growth that can be
            ! acheived with the N available to plants
            ! calculate the fraction of immobilization realized (for diagnostic purposes)
#ifndef FUN
            
               if (plant_ndemand(i) > 0.0_r8) then
                  fpg(i) = sminn_to_plant(i) / plant_ndemand(i)
               else
                  fpg(i) = 1._r8
               end if
#endif

            if (potential_immob > 0.0_r8) then
               fpi(i) = actual_immob / potential_immob
            else
               fpi(i) = 1._r8
            end if
#endif


  end subroutine SoilBiogeochemCompetition

end module bgc_soil_SoilBiogeochemCompetitionMod
#endif
