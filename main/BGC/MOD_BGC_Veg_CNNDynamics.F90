#include <define.h>
#ifdef BGC
module MOD_BGC_Veg_CNNDynamics

  !---------------------------------------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! This module simulates the plant biological fixation (Cleveland et al., 1999), 
  ! crop fertilisatoin(Lawrence et al., 2016, and soy nitrogen fixation (Neitsch et al., 2005).
  !
  ! !REFERENCE:
  ! Cleveland, C.C., Townsend, A.R., Schimel, D.S., Fisher, H., Howarth, R.W., Hedin, L.O., Perakis, S.S., Latty, E.F.,
  ! Von Fischer, J.C., Elseroad, A., and Wasson, M.F. 1999. Global patterns of terrestrial biological nitrogen (N2) fixation
  ! in natural ecosystems. Global Biogeochem. Cycles 13:623-645.
  ! Lawrence, D.M., Hurtt, G.C., Arneth, A., Brovkin, V., Calvin, K.V., Jones, A.D., Jones, C.D., Lawrence, P.J., de
  ! Noblet-Ducoudré, N., Pongratz, J., Seneviratne, S.I., and Shevliakova, E. 2016. The Land Use Model Intercomparison
  ! Project (LUMIP) contribution to CMIP6: rationale and experimental design. Geosci. Model Dev. 9:2973-2998.
  ! DOI:10.5194/gmd-9-2973-2016.
  ! Neitsch, S.L., Arnold, J.G., Kiniry, J.R., and Williams J.R. 2005. Soil and Water Assessment Tool, 
  ! Theoretical Documentation: Version 2005. Temple, TX. USDA Agricultural Research Service and 
  ! Texas A&M Blackland Research Center.
  !
  ! !ORIGINAL:
  ! The Community Land Model version 5.0 (CLM5)
  !
  ! !REVISION:
  ! Xingjie Lu, 2021, revised the CLM5 code to be compatible with CoLM code sturcture.

  use MOD_Precision

  use MOD_Vars_PFTimeInvariants, only: pftclass, pftfrac
#ifdef CROP
  use MOD_BGC_Vars_PFTimeVariables, only: fertnitro_p, fert_p
#endif
  use MOD_Vars_TimeInvariants, only: porsl, psi0, bsw
  use MOD_Vars_TimeVariables, only: h2osoi

  use MOD_BGC_Vars_1DFluxes, only: fert_to_sminn, soyfixn_to_sminn, nfix_to_sminn

  use MOD_BGC_Vars_TimeVariables, only: sminn, fpg, lag_npp
#ifdef CROP
  use MOD_BGC_Vars_PFTimeVariables, only: croplive_p, hui_p
#endif

  use MOD_BGC_Vars_1DPFTFluxes, only: plant_ndemand_p, soyfixn_p

  use MOD_Vars_Global, only: z_soi, dz_soi, spval
  use MOD_TimeManager
  implicit none

  public CNNFixation
#ifdef CROP
  public CNNFert
  public CNSoyfix
#endif

contains

  subroutine CNNFixation(i,idate)

    integer ,intent(in) :: i          ! patch index
    integer ,intent(in) :: idate(3)   ! current date (year, day of the year, seconds of the day)
    real(r8) t, dayspyr

    if(isleapyear(idate(1)))then
       dayspyr = 366
    else
       dayspyr = 365
    end if

    if (lag_npp(i) /= spval) then
       ! need to put npp in units of gC/m^2/year here first
       t = (1.8_r8 * (1._r8 - exp(-0.003_r8 * lag_npp(i)*(86400._r8 * dayspyr))))/(86400._r8 * dayspyr)
       nfix_to_sminn(i) = max(0._r8,t)
    else
       nfix_to_sminn(i) = 0._r8
    endif

  end subroutine CNNFixation

#ifdef CROP
  subroutine CNNFert(i,ps,pe)

    integer ,intent(in) :: i  ! patch index
    integer ,intent(in) :: ps ! start pft index
    integer ,intent(in) :: pe ! end pft index

    fert_to_sminn(i) = sum(fert_p(ps:pe))

  end subroutine CNNFert

  subroutine CNSoyfix (i, ps, pe, nl_soil)
    ! GPAM Soybean biological N fixation
    integer, intent(in) :: i       ! patch index
    integer, intent(in) :: ps      ! start pft index
    integer, intent(in) :: pe      ! end pft index
    integer, intent(in) :: nl_soil ! number of total soil layers
    real(r8):: fxw,fxn,fxg,fxr             ! soil water factor, nitrogen factor, growth stage factor
    real(r8):: soy_ndemand                 ! difference between nitrogen supply and demand
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
          watdry = porsl(j,i) * (316230._r8/(-psi0(j,i))) ** (-1._r8/bsw(j,i))
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
       if(croplive_p(m) .and. ivt ==  23 .or. ivt == 24 .or. ivt == 77 .or. ivt == 78)then 

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

             if (hui_p(m) <= GDDfracthreshold1) then
                fxg = 0._r8
             else if (hui_p(m) > GDDfracthreshold1 .and. hui_p(m) <= GDDfracthreshold2) then
                fxg = 6.67_r8 * hui_p(m) - 1._r8
             else if (hui_p(m) > GDDfracthreshold2 .and. hui_p(m) <= GDDfracthreshold3) then
                fxg = 1._r8
             else if (hui_p(m) > GDDfracthreshold3 .and. hui_p(m) <= GDDfracthreshold4) then
                fxg = 3.75_r8 - 5._r8 * hui_p(m)
             else 
                fxg = 0._r8
             end if

             ! calculate the nitrogen fixed by the soybean

             fxr = max(0._r8, min(1._r8, fxw, fxn) * fxg)
             soyfixn_p(m) = min(fxr * soy_ndemand, soy_ndemand)

          else ! if nitrogen demand met, no fixation

             soyfixn_p(m) = 0._r8

          end if

       else ! if not live soybean, no fixation

          soyfixn_p(m) = 0._r8

       end if
    end do

    soyfixn_to_sminn(i) = sum(soyfixn_p(ps:pe)*pftfrac(ps:pe))

  end subroutine CNSoyfix
#endif

end module MOD_BGC_Veg_CNNDynamics
#endif
