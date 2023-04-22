#include <define.h>

 SUBROUTINE rain_snow_temp (itypwat,&
                            forc_t,forc_q,forc_psrf,forc_prc,forc_prl,tcrit,&
                            prc_rain,prc_snow,prl_rain,prl_snow,t_precip,bifall)

!=======================================================================
! define the rate of rainfall and snowfall and precipitation water temp
! Original author : Yongjiu Dai, 09/1999; 08/31/2002, 04/2014, 01/2023
!=======================================================================
!
  use precision
  use PhysicalConstants, only : tfrz
 
  IMPLICIT NONE
 
! ------------------------ Dummy Argument ------------------------------
  integer, INTENT(in) :: itypwat     ! land water type (3=glaciers)


  real(r8), INTENT(in) :: forc_t     ! temperature at agcm reference height [kelvin]
  real(r8), INTENT(in) :: forc_q     ! specific humidity at agcm reference height [kg/kg]
  real(r8), INTENT(in) :: forc_psrf  ! atmosphere pressure at the surface [pa]
  real(r8), INTENT(in) :: forc_prc   ! convective precipitation [mm/s]
  real(r8), INTENT(in) :: forc_prl   ! large scale precipitation [mm/s]

  real(r8), INTENT(in) :: tcrit      ! critical temp. to determine rain or snow

  real(r8), INTENT(out) :: prc_rain  ! convective rainfall [kg/(m2 s)]
  real(r8), INTENT(out) :: prc_snow  ! convective snowfall [kg/(m2 s)]
  real(r8), INTENT(out) :: prl_rain  ! large scale rainfall [kg/(m2 s)]
  real(r8), INTENT(out) :: prl_snow  ! large scale snowfall [kg/(m2 s)]
  real(r8), INTENT(out) :: t_precip  ! snowfall/rainfall temperature [kelvin]
  real(r8), INTENT(out) :: bifall    ! bulk density of newly fallen dry snow [kg/m3]
  real(r8) :: flfall  ! fraction of liquid water within falling precip.

  real(r8) :: all_snow_t   ! temperature at which all precip falls entirely as snow (K)
  real(r8) :: frac_rain_slope ! slope of the frac_rain vs. temperature relationship
  real(r8) :: all_snow_t_c ! Temperature at which precip falls entirely as rain (deg C)
  real(r8) :: all_rain_t_c ! Temperature at which precip falls entirely as snow (deg C)

  logical :: glaciers    ! true: glacier column
!-----------------------------------------------------------------------

! wet-bulb temperature
      call wetbulb(forc_t,forc_psrf,forc_q,t_precip)

#IF(DEFINED option_precip_phase_discrimination_I)
! Wang, Y.H., Broxton, P., Fang, Y., Behrangi, A., Barlage, M., Zeng, X., & Niu, G.Y. (2019).
! A Wet-Bulb Temperature Based Rain-Snow Partitioning Scheme Improves Snowpack Prediction
! Over the Drier Western United States. Geophysical Research Letters, 46, 13,825-13,835.
!
! Behrangi et al. (2018) On distinguishing snowfall from rainfall
! using near-surface atmospheric information: Comparative analysis,
! uncertainties and hydrologic importance. Q J R Meteorol Soc. 144 (Suppl. 1):89-102

      if(t_precip > tfrz + 3.0)then
         flfall = 1.0      ! fraction of liquid water within falling precip
         bifall = 169.15   ! (not used)
      else if (t_precip >= tfrz -2.0)then
         flfall = max(0.0, 1.0 - 1.0/(1.0+5.00e-5*exp(2.0*(t_precip+4.))))   !Figure 5c of Behrangi et al. (2018)
!*       flfall = max(0.0, 1.0 - 1.0/(1.0+6.99e-5*exp(2.0*(t_precip+3.97)))) !Equation 1 of Wang et al. (2019)
         bifall = min(169.15, 50. + 1.7*(max(0.0, t_precip-tfrz+15.))**1.5)
      else
         flfall = 0.0
         bifall = 50.0
      endif

#ELIF(DEFINED option_precip_phase_discrimination_II)
! CLM5.0 
      glaciers = .false.
      if (itypwat == 3) glaciers = .true.

      if(glaciers) then
         all_snow_t_c = -2.0
         all_rain_t_c =  0.0
      else
         all_snow_t_c = 0.0
         all_rain_t_c = 2.0
      endif

      all_snow_t = all_snow_t_c + tfrz
      frac_rain_slope = 1._r8 / (all_rain_t_c - all_snow_t_c)

      ! Re-partition precipitation into rain/snow for a single column.
      ! Rain and snow variables should be set initially, and are updated here

      flfall = min(1.0_r8, max(0.0_r8,(forc_t - all_snow_t)*frac_rain_slope))
      bifall = min(169.0_r8, 50. + 1.7*(max(0.0, forc_t-tfrz+15.))**1.5)

#ELSE
! the upper limit of air temperature is set for snowfall, this cut-off
! was selected based on Fig. 1, Plate 3-1, of Snow Hydrology (1956).
! the percentage of liquid water by mass, which is arbitrarily set to
! vary linearly with air temp, from 0% at 273.16 to 40% max at 275.16.

      if(forc_t>tfrz+2.0)then
         flfall = 1.0     ! fraction of liquid water within falling precip.
         bifall = 169.15  ! (not used)
      else
         flfall = max(0.0, -54.632+0.2*forc_t)
         bifall = 50. + 1.7*(max(0.0, forc_t-tfrz+15.))**1.5
      end if

#ENDIF

      prc_rain = forc_prc*flfall        ! convective rainfall (mm/s)
      prl_rain = forc_prl*flfall        ! large scale rainfall (mm/s)
      prc_snow = forc_prc*(1.-flfall)   ! convective snowfall (mm/s)            
      prl_snow = forc_prl*(1.-flfall)   ! large scale snowfall (mm/s)              

! -------------------------------------------------------------
! temperature of rainfall or snowfall
! -------------------------------------------------------------

      if (forc_t > 275.65) then
         if (t_precip < tfrz) t_precip = tfrz
      else
         t_precip = min(tfrz,t_precip)
         if(flfall > 0.0)then
           t_precip = tfrz - sqrt((1.0/flfall)-1.0)/100.0
         endif
      endif

 END SUBROUTINE rain_snow_temp
