#include <define.h>
#ifdef BGC
module bgc_soil_SoilBiogeochemDecompCascadeBGCMod

use precision
use MOD_TimeInvariants, only: &
    Q10, smpmax, smpmin, tau_l1, tau_l2_l3, tau_s1, tau_s2, tau_s3, tau_cwd, froz_q10, &
    i_met_lit,i_cel_lit,i_lig_lit ,i_cwd,i_soil1,i_soil2,i_soil3
use MOD_TimeVariables, only: &
    smp, t_soisno, t_scalar, w_scalar, o_scalar, depth_scalar, decomp_k
use GlobalVars, only: PI

implicit none

public decomp_rate_constants_bgc

contains

   subroutine decomp_rate_constants_bgc(i,nl_soil,z_soi)

integer ,intent(in) :: i
integer ,intent(in) :: nl_soil
real(r8),intent(in) :: z_soi(1:nl_soil)

real(r8) normalization_factor ! factor by which to offset the decomposition rates frm century to a q10 formulation
real(r8),parameter :: decomp_depth_efolding = 10._r8
real(r8) k_l1, k_l2_l3, k_s1, k_s2, k_s3, k_frag
real(r8) psi
integer j
real(r8) catanf
real(r8) catanf_30
real(r8) t1

catanf(t1) = 11.75_r8 +(29.7_r8 / PI) * atan( PI * 0.031_r8  * ( t1 - 15.4_r8 ))


!      tau_l1 = 1./18.5
!      tau_l2_l3 = 1./4.9
!      tau_s1 = 1./7.3
!      tau_s2 = 1./0.2
!      tau_s3 = 1./.0045

      ! century leaves wood decomposition rates open, within range of 0 - 0.5 yr^-1
!      tau_cwd  = 1./0.3

      ! set "Q10" parameter
!      Q10 = CNParamsShareInst%Q10

      ! set "froz_q10" parameter
!      froz_q10  = CNParamsShareInst%froz_q10 ! = 1.5

      ! Set "decomp_depth_efolding" parameter
!      decomp_depth_efolding = CNParamsShareInst%decomp_depth_efolding

      ! translate to per-second time constant
      k_l1 = 1._r8    / (86400._r8 * 365._r8 * tau_l1)
      k_l2_l3 = 1._r8 / (86400._r8 * 365._r8 * tau_l2_l3)
      k_s1 = 1._r8    / (86400._r8 * 365._r8 * tau_s1)
      k_s2 = 1._r8    / (86400._r8 * 365._r8 * tau_s2)
      k_s3 = 1._r8    / (86400._r8 * 365._r8 * tau_s3)
      k_frag = 1._r8  / (86400._r8 * 365._r8 * tau_cwd)

      ! calc ref rate
      catanf_30 = catanf(30._r8)

!      if ( spinup_state >= 1 ) then
!          do fc = 1,num_soilc
!             c = filter_soilc(fc)
!             !
!             if ( abs(spinup_factor(i_litr1) - 1._r8) .gt. .000001_r8) then
!                spinup_geogterm_l1(c) = spinup_factor(i_litr1) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
!             else
!                spinup_geogterm_l1(c) = 1._r8
!             endif
!             !
!             if ( abs(spinup_factor(i_litr2) - 1._r8) .gt. .000001_r8) then
!                spinup_geogterm_l23(c) = spinup_factor(i_litr2) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
!             else
!                spinup_geogterm_l23(c) = 1._r8
!             endif
!             !
!             if ( .not. use_fates ) then
!                if ( abs(spinup_factor(i_cwd) - 1._r8) .gt. .000001_r8) then
!                   spinup_geogterm_cwd(c) = spinup_factor(i_cwd) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
!                else
!                   spinup_geogterm_cwd(c) = 1._r8
!                endif
!             endif
!             !
!             if ( abs(spinup_factor(i_soil1) - 1._r8) .gt. .000001_r8) then
!                spinup_geogterm_s1(c) = spinup_factor(i_soil1) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
!             else
!                spinup_geogterm_s1(c) = 1._r8
!             endif
!             !
!             if ( abs(spinup_factor(i_soil2) - 1._r8) .gt. .000001_r8) then
!                spinup_geogterm_s2(c) = spinup_factor(i_soil2) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
!             else
!                spinup_geogterm_s2(c) = 1._r8
!             endif
!             !
!             if ( abs(spinup_factor(i_soil3) - 1._r8) .gt. .000001_r8) then
!                spinup_geogterm_s3(c) = spinup_factor(i_soil3) * get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
!             else
!                spinup_geogterm_s3(c) = 1._r8
!             endif
!             !
!          end do
!       else
!          do fc = 1,num_soilc
!             c = filter_soilc(fc)
!             spinup_geogterm_l1(c) = 1._r8
!             spinup_geogterm_l23(c) = 1._r8
!             spinup_geogterm_cwd(c) = 1._r8
!             spinup_geogterm_s1(c) = 1._r8
!             spinup_geogterm_s2(c) = 1._r8
!             spinup_geogterm_s3(c) = 1._r8
!          end do
!       endif

      ! calculate rate constant scalar for soil temperature
      ! assuming that the base rate constants are assigned for non-moisture
      ! limiting conditions at 25 C.
      ! Peter Thornton: 3/13/09
      ! Replaced the Lloyd and Taylor function with a Q10 formula, with Q10 = 1.5
      ! as part of the modifications made to improve the seasonal cycle of
      ! atmospheric CO2 concentration in global simulations. This does not impact
      ! the base rates at 25 C, which are calibrated from microcosm studies.

!      do j = 1,nl_soil
!         if (j==1) t_scalar(:) = 0._r8
!         if (t_soisno(j) >= 273.15_r8) then
!            t_scalar(1)=t_scalar(1) + &
!                           (Q10**((t_soisno(j)-(273.15_r8+25._r8))/10._r8))*fr(j)
!         else
!            t_scalar(1)=t_scalar(1) + &
!                           (Q10**(-25._r8/10._r8))*(froz_q10**((t_soisno(j)-273.15_r8)/10._r8))*fr(j)
!         endif
!      end do

       do j = 1, nl_soil
          if (t_soisno(j,i) >= 273.15_r8) then
             t_scalar(j,i)= (Q10**((t_soisno(j,i)-(273.15_r8+25._r8))/10._r8))
          else
             t_scalar(j,i)= (Q10**(-25._r8/10._r8))*(froz_q10**((t_soisno(j,i)-273.15_r8)/10._r8))
          endif
       end do


      ! calculate the rate constant scalar for soil water content.
      ! Uses the log relationship with water potential given in
      ! Andren, O., and K. Paustian, 1987. Barley straw decomposition in the field:
      ! a comparison of models. Ecology, 68(5):1190-1200.
      ! and supported by data in
      ! Orchard, V.A., and F.J. Cook, 1983. Relationship between soil respiration
      ! and soil moisture. Soil Biol. Biochem., 15(4):447-453.

      do j = 1,nl_soil
         psi = min(smp(j,i),smpmax)
         ! decomp only if soilpsi is higher than minpsi
         if (psi > smpmin) then
            w_scalar(j,i) = (log(smpmin/psi)/log(smpmin/smpmax))
         else
            w_scalar(j,i) = 0._r8
         end if
!            if (use_lch4) then
!               if (anoxia_wtsat .and. t_soisno(c,j) > 273.15_r8) then ! wet area will have w_scalar of 1 if unfrozen
!                  w_scalar(c,j) = w_scalar(c,j)*(1._r8 - finundated(c)) + finundated(c)
!               end if
!            end if
      end do

!      if (use_lch4) then
!         ! Calculate ANOXIA
!         ! Check for anoxia w/o LCH4 now done in controlMod.
!
!         if (anoxia) then
!            do j = 1, nl_soil
!               do fc = 1,num_soilc
!                  c = filter_soilc(fc)
!
!                  if (.not. anoxia_wtsat) then
!                     o_scalar(c,j) = max(o2stress_unsat(c,j), mino2lim)
!                  else
!                     o_scalar(c,j) = max(o2stress_unsat(c,j), mino2lim) * (1._r8 - finundated(c)) + &
!                          max(o2stress_sat(c,j), mino2lim) * finundated(c)
!                  end if
!               end do
!            end do
!         else
!            o_scalar(bounds%begc:bounds%endc,1: nl_soil) = 1._r8
!         end if
!      else
         o_scalar(1:nl_soil,i) = 1._r8
!      end if

      ! scale all decomposition rates by a constant to compensate for offset between original CENTURY temp func and Q10
      normalization_factor = (catanf(15._r8)/catanf_30) / (Q10**((15._r8-25._r8)/10._r8))
      do j = 1, nl_soil
         t_scalar(j,i) = t_scalar(j,i) * normalization_factor
      end do

      do j = 1, nl_soil
         depth_scalar(j,i) = exp(-z_soi(j)/decomp_depth_efolding)
      end do

      do j = 1, nl_soil
         decomp_k(j,i_met_lit,i) = k_l1    * t_scalar(j,i) * w_scalar(j,i) * depth_scalar(j,i) * o_scalar(j,i) !&
!                                 * spinup_geogterm_l1(c)
         decomp_k(j,i_cel_lit,i) = k_l2_l3 * t_scalar(j,i) * w_scalar(j,i) * depth_scalar(j,i) * o_scalar(j,i) !&
!                                 * spinup_geogterm_l23(c)
         decomp_k(j,i_lig_lit,i) = k_l2_l3 * t_scalar(j,i) * w_scalar(j,i) * depth_scalar(j,i) * o_scalar(j,i) !&
!                                 * spinup_geogterm_l23(c)
         decomp_k(j,i_soil1  ,i) = k_s1    * t_scalar(j,i) * w_scalar(j,i) * depth_scalar(j,i) * o_scalar(j,i) !&
!                                 * spinup_geogterm_s1(c)
         decomp_k(j,i_soil2  ,i) = k_s2    * t_scalar(j,i) * w_scalar(j,i) * depth_scalar(j,i) * o_scalar(j,i) !&
!                                 * spinup_geogterm_s2(c)
         decomp_k(j,i_soil3  ,i) = k_s3    * t_scalar(j,i) * w_scalar(j,i) * depth_scalar(j,i) * o_scalar(j,i) !&
!                                 * spinup_geogterm_s3(c)
!            if(use_soil_matrixcn)then
!               Ksoil%DM(c,j+ nl_soil*(i_litr1-1)) = k_l1    * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) * spinup_geogterm_l1(     c) * dt
!               Ksoil%DM(c,j+ nl_soil*(i_litr2-1)) = k_l2_l3 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) * spinup_geogterm_l23     (c) * dt
!               Ksoil%DM(c,j+ nl_soil*(i_litr3-1)) = k_l2_l3 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) * spinup_geogterm_l23     (c) * dt
!               Ksoil%DM(c,j+ nl_soil*(i_soil1-1)) = k_s1    * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) * spinup_geogterm_s1(     c) * dt
!               Ksoil%DM(c,j+ nl_soil*(i_soil2-1)) = k_s2    * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) * spinup_geogterm_s2(     c) * dt
!               Ksoil%DM(c,j+ nl_soil*(i_soil3-1)) = k_s3    * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) * spinup_geogterm_s3(     c) * dt
!            end if !use_soil_matrixcn
      end do

      do j = 1,nl_soil
         decomp_k(j,i_cwd,i)   = k_frag  * t_scalar(j,i) * w_scalar(j,i) * depth_scalar(j,i) * &
                 o_scalar(j,i) ! * spinup_geogterm_cwd(i)
!            if(use_soil_matrixcn)then
!               Ksoil%DM(c,j+ nl_soil*(i_cwd-1))   = k_frag  * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * &
!                      o_scalar(c,j) * spinup_geogterm_cwd(c) * dt
!            end if !use_soil_matrixcn
      end do

   end subroutine decomp_rate_constants_bgc

end module bgc_soil_SoilBiogeochemDecompCascadeBGCMod
#endif
