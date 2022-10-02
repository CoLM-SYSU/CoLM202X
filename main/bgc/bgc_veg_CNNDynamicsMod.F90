#include <define.h>
#ifdef BGC
#ifdef CROP
module bgc_veg_CNNDynamicsMod

use precision

use MOD_PFTimeInvars, only: pftclass
use MOD_TimeInvariants, only: porsl, psi0, bsw
use MOD_TimeVariables, only: h2osoi

use MOD_1D_BGCPFTFluxes, only: fert_p  ! intent(in)
 
use MOD_1D_BGCFluxes, only: fert_to_sminn

use MOD_BGCTimeVars, only: sminn, fpg
use MOD_BGCPFTimeVars, only: croplive_p, gddplant_p, gddmaturity_p

use MOD_1D_BGCPFTFluxes, only: plant_ndemand_p, soyfixn_p

use GlobalVars, only: z_soi, dz_soi
implicit none

public CNNFert

contains

 subroutine CNNFert(i,ps,pe)

  integer ,intent(in) :: i
  integer ,intent(in) :: ps
  integer ,intent(in) :: pe

  fert_to_sminn(i) = sum(fert_p(ps:pe))

 end subroutine CNNFert

 subroutine CNSoyfix (i, ps, pe, nl_soil)

   integer, intent(in) :: i
   integer, intent(in) :: ps
   integer, intent(in) :: pe
   integer, intent(in) :: nl_soil
   real(r8):: fxw,fxn,fxg,fxr             ! soil water factor, nitrogen factor, growth stage factor
   real(r8):: soy_ndemand                 ! difference between nitrogen supply and demand
   real(r8):: GDDfrac
   real(r8):: sminnthreshold1, sminnthreshold2
   real(r8):: GDDfracthreshold1, GDDfracthreshold2
   real(r8):: GDDfracthreshold3, GDDfracthreshold4
   integer m, ivt, j
   real(r8) :: rwat, swat, rz, watdry, wf, tsw, stsw

   sminnthreshold1 = 30._r8
   sminnthreshold2 = 10._r8
   GDDfracthreshold1 = 0.15_r8
   GDDfracthreshold2 = 0.30_r8
   GDDfracthreshold3 = 0.55_r8
   GDDfracthreshold4 = 0.75_r8

  rwat = 0._r8
  swat = 0._r8
  rz   = 0._r8
 
   do j = 1, nl_soil
      if (z_soi(j)+0.5_r8*dz_soi(j) <= 0.05_r8) then
         watdry = porsl(j,i) * (316230._r8/psi0(j,i)) ** (-1._r8/bsw(j,i))
         rwat = rwat + (h2osoi(j,i)-watdry) * dz_soi(j)
         swat = swat + (porsl (j,i)-watdry) * dz_soi(j)
         rz   = rz   + dz_soi(j)
      end if
   end do

   tsw  = rwat/rz
   stsw = swat/rz
   wf = tsw/stsw
            
   do m = ps, pe
      ivt = pftclass(m)
      if(croplive_p(m) .and. ivt .eq. 23 .or. ivt .eq. 24 .or. ivt .eq. 77 .or. ivt .eq. 78)then 
        ! & rain fed temperate soybean & irrigated temperate soybean 
        ! & rain fed tropical soybean & irrigated tropical soybean

        ! difference between supply and demand

         if(fpg(i) .lt. 1._r8) then
            soy_ndemand = plant_ndemand_p(m) - plant_ndemand_p(m) * fpg(i)
            
           ! fixation depends on nitrogen, soil water, and growth stage
           ! soil water factor

            fxw = wf / 0.85_r8

           ! soil nitrogen factor (Beth says: CHECK UNITS)

            if (sminn(i) .gt. sminnthreshold1) then
               fxn = 0._r8
            else if (sminn(i) > sminnthreshold2 .and. sminn(i) <= sminnthreshold1) then
               fxn = 1.5_r8 - .005_r8 * (sminn(i) * 10._r8)
            else if (sminn(i) <= sminnthreshold2) then
               fxn = 1._r8
            end if
           
            ! growth stage factor
            ! slevis: to replace GDDfrac, assume...
            ! Beth's crit_offset_gdd_def is similar to my gddmaturity
            ! Beth's ac_gdd (base 5C) similar to my hui=gddplant (base 10
            ! for soy)
            ! Ranges below are not firm. Are they lit. based or tuning based?

            GDDfrac = gddplant_p(m) / gddmaturity_p(m)

            if (GDDfrac <= GDDfracthreshold1) then
               fxg = 0._r8
            else if (GDDfrac > GDDfracthreshold1 .and. GDDfrac <= GDDfracthreshold2) then
               fxg = 6.67_r8 * GDDfrac - 1._r8
            else if (GDDfrac > GDDfracthreshold2 .and. GDDfrac <= GDDfracthreshold3) then
               fxg = 1._r8
            else if (GDDfrac > GDDfracthreshold3 .and. GDDfrac <= GDDfracthreshold4) then
               fxg = 3.75_r8 - 5._r8 * GDDfrac
            else  ! GDDfrac > GDDfracthreshold4
               fxg = 0._r8
            end if

            ! calculate the nitrogen fixed by the soybean

            fxr = min(1._r8, fxw, fxn) * fxg
            fxr = max(0._r8, fxr)
            soyfixn_p(m) = fxr * soy_ndemand
            soyfixn_p(m) = min(soyfixn_p(m), soy_ndemand)

         else ! if nitrogen demand met, no fixation

            soyfixn_p(m) = 0._r8

         end if

      else ! if not live soybean, no fixation

         soyfixn_p(m) = 0._r8

      end if
   end do

 end subroutine CNSoyfix

end module bgc_veg_CNNDynamicsMod
#endif
#endif
