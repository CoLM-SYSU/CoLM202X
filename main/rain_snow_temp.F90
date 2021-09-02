
 SUBROUTINE rain_snow_temp (forc_t,forc_q,forc_psrf,forc_prc,forc_prl,tcrit,&
                            prc_rain,prc_snow,prl_rain,prl_snow,t_precip,bifall)

!=======================================================================
! define the rate of rainfall and snowfall and precipitation water temp
! Original author : Yongjiu Dai, 09/1999; 08/31/2002, 04/2014
!=======================================================================
!
  use precision
  use PhysicalConstants, only : tfrz
 
  IMPLICIT NONE
 
! ------------------------ Dummy Argument ------------------------------

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

!-----------------------------------------------------------------------
! the upper limit of air temperature is set for snowfall, this cut-off 
! was selected based on Fig. 1, Plate 3-1, of Snow Hydrology (1956).
! the percentage of liquid water by mass, which is arbitrarily set to 
! vary linearly with air temp, from 0% at 273.16 to 40% max at 275.16.

      if(forc_t>tfrz+tcrit)then
        prc_rain = forc_prc       ! convective rainfall (mm/s)
        prl_rain = forc_prl       ! large scale rainfall (mm/s)
        prc_snow = 0.             ! convective snowfall (mm/s)
        prl_snow = 0.             ! large scale snowfall (mm/s)
        flfall = 1.               ! fraction of liquid water within falling precip.
        bifall = 169.             ! (not used)
      else
        if(forc_t<=tfrz)then
          flfall = 0.
        else if(forc_t<=tfrz+2.)then
          flfall = -54.632+0.2*forc_t
        else
          flfall = 0.4
        endif
 
! use Alta relationship, Anderson(1976); LaChapelle(1961), 
! U.S.Department of Agriculture Forest Service, Project F, 
! Progress Rep. 1, Alta Avalanche Study Center:Snow Layer Densification.

        if(forc_t>tfrz+2.)then
          bifall = 169.
        else if(forc_t>tfrz-15.)then
          bifall = 50.+1.7*(forc_t-tfrz+15.)**1.5
        else
          bifall = 50.
        endif
 
        prc_rain = forc_prc*flfall
        prl_rain = forc_prl*flfall
        prc_snow = forc_prc*(1.-flfall)                 
        prl_snow = forc_prl*(1.-flfall)                 

      endif


! -------------------------------------------------------------
! temperature of rainfall or snowfall
! -------------------------------------------------------------

      call wetbulb(forc_t,forc_psrf,forc_q,t_precip)

      if (forc_t > 275.65) then
         if (t_precip < tfrz) t_precip = tfrz
      else
         t_precip = min(tfrz,t_precip)
         if(flfall > 0.0)then
           t_precip = tfrz - sqrt((1.0/flfall)-1.0)/100.0
         endif
      endif


 END SUBROUTINE rain_snow_temp
