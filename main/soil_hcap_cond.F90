#include <define.h>
!====================================================================
SUBROUTINE soil_hcap_cond(vf_gravels_s,vf_om_s,vf_sand_s,vf_pores_s,&
                                    wf_gravels_s,wf_sand_s,k_solids,&
                                            csol,kdry,ksat_u,ksat_f,&
#ifdef THERMAL_CONDUCTIVITY_SCHEME_4
                                                   BA_alpha,BA_beta,&
#endif
                               temperature,vf_water,vf_ice,hcap,thk)
!====================================================================
! Reference: Dai et al.,2019:
! Evaluation of Soil Thermal Conductivity Schemes for Use in Land Surface Modeling
!
! Yongjiu Dai, 02/2018, 06/2018
! -------------------------------------------------------------------
use precision
USE PhysicalConstants,only:tfrz

IMPLICIT NONE
      real(r8), intent(in) :: vf_gravels_s ! volumetric fraction of gravels within the soil solids
      real(r8), intent(in) :: vf_om_s      ! volumetric fraction of organic matter within the soil solids
      real(r8), intent(in) :: vf_sand_s    ! volumetric fraction of sand within soil soilds
      real(r8), intent(in) :: vf_pores_s   ! volumetric pore space of the soil

      real(r8), intent(in) :: wf_gravels_s ! gravimetric fraction of gravels
      real(r8), intent(in) :: wf_sand_s    ! gravimetric fraction of sand within soil soilds
      real(r8), intent(in) :: k_solids     ! thermal conductivity of soil solids

      real(r8), intent(in) :: temperature  !
      real(r8), intent(in) :: vf_water     !
      real(r8), intent(in) :: vf_ice       !

      real(r8), intent(in) :: csol         ! heat capacity of dry soil [J/(m3 K)]
      real(r8), intent(in) :: kdry         ! thermal conductivity for dry soil [W/m/K]
      real(r8), intent(in) :: ksat_u       ! thermal conductivity of unfrozen saturated soil [W/m/K]
      real(r8), intent(in) :: ksat_f       ! thermal conductivity of frozen saturated soil [W/m/K]
#ifdef THERMAL_CONDUCTIVITY_SCHEME_4
      real(r8), intent(in) :: BA_alpha     ! alpha in Balland and Arp(2005) thermal conductivity scheme
      real(r8), intent(in) :: BA_beta      ! beta in Balland and Arp(2005) thermal conductivity scheme
#endif

      real(r8), intent(out) :: hcap        ! J/(m3 K)
      real(r8), intent(out) :: thk         ! W/(m K)

      real(r8) c_water, c_ice
      real(r8) a, kappa, alpha, beta
      real(r8) aa,aaa,nwm,nw_nwm,x,ga,gc
      integer i

      real(r8) sr ! wetness or degree of saturation = (vf_water+vf_ice)/vf_pores_s
      real(r8) ke ! Kersten number or normalized thermal conductivity
      real(r8) k_air,k_water,k_ice

! =========================================================================================
! The heat capacity and thermal conductivity [J(m3 K)]
! =========================================================================================
!*    c_water = 4.18e6 ! J/(m3 K)
!*    c_ice = 1.88e6   ! J/(m3 K)
      c_water = 4.188e6   ! J/(m3 K) = 4188[J/(kg K)]*1000(kg/m3)
      c_ice = 1.94153e6   ! J/(m3 K) = 2117.27[J/(kg K)]*917(kg/m3)


      hcap = csol + vf_water*c_water + vf_ice*c_ice 

! -----------------------------------------------------------------------------------------
! Setting
! -----------------------------------------------------------------------------------------
      k_air = 0.024    ! (W/m/K)
      k_water = 0.57   ! (W/m/K)
      k_ice = 2.29     ! (W/m/K)

      a =  vf_gravels_s + vf_sand_s

      sr = (vf_water+vf_ice)/vf_pores_s
!      sr = max(1.0e-6, sr)
      sr = min(1.0, sr)
if(sr >= 1.0e-10) then
#if(defined THERMAL_CONDUCTIVITY_SCHEME_1)
! -----------------------------------------------------------------------------------------
! [1] Oleson et al., 2013: Technical Description of version 4.5 of the Community Land Model
!     (CLM). NCAR/TN-503+STR (Section 6.3: Soil and Snow Thermal Properties)
! -----------------------------------------------------------------------------------------
      if(temperature > tfrz)then ! Unfrozen soil
         ke = log10(sr) + 1.0
      else                         ! Fozen or partially frozen soils
         ke = sr
      end if
#endif

#if(defined THERMAL_CONDUCTIVITY_SCHEME_2)
! -----------------------------------------------------------------------------------------
! [2] Johansen O (1975): Thermal conductivity of soils. PhD Thesis. Trondheim, Norway:
!     University of Trondheim. US army Crops of Engineerings,
!     CRREL English Translation 637.
! -----------------------------------------------------------------------------------------
      if(temperature > tfrz)then ! Unfrozen soils
         if(a > 0.4)then ! coarse-grained
            ke = 0.7*log10(max(sr,0.05)) + 1.0
         else            ! Fine-grained
            ke = log10(max(sr,0.1)) + 1.0
         endif
      else                         ! Fozen or partially frozen soils
         ke = sr
      endif
#endif

#if(defined THERMAL_CONDUCTIVITY_SCHEME_3)
! -----------------------------------------------------------------------------------------
! [3] Cote, J., and J.-M. Konrad (2005), A generalized thermal conductivity model for soils
!     and construction materials. Canadian Geotechnical Journal, 42(2): 443-458.
! -----------------------------------------------------------------------------------------
      if(temperature > tfrz)then ! Unfrozen soils
!        kappa =                       Unfrozen
!        /gravels and coarse sand     /4.60/
!        /medium and fine sands       /3.55/
!        /silty and clayey soils      /1.90/
!        /organic fibrous soils (peat)/0.60/
         if(a > 0.40)then
            kappa = 4.60
         else if(a > 0.25)then
            kappa = 3.55
         else if(a > 0.01)then
            kappa = 1.90
         else
            kappa = 0.60
         endif

      else                         ! Fozen or partially frozen soils
!        kappa =                      Frozen
!       /gravels and coarse sand     /1.70/
!       /medium and fine sands       /0.95/
!       /silty and clayey soils      /0.85/
!       /organic fibrous soils (peat)/0.25/
         if(a > 0.40)then
            kappa = 1.70
         else if(a > 0.25)then
            kappa = 0.95
         else if(a > 0.01)then
            kappa = 0.85
         else
            kappa = 0.25
         endif
      endif
      ke = kappa*sr/(1.0+(kappa-1.0)*sr)
#endif

#if(defined THERMAL_CONDUCTIVITY_SCHEME_4)
! -----------------------------------------------------------------------------------------
! [4] Balland V. and P. A. Arp, 2005: Modeling soil thermal conductivities over a wide
! range of conditions. J. Environ. Eng. Sci. 4: 549-558.
! be careful in specifying all k affecting fractions as VOLUME FRACTION,
! whether these fractions are part of the bulk volume, the pore space, or the solid space.
! -----------------------------------------------------------------------------------------
      if(temperature > tfrz)then ! Unfrozen soil
!         alpha = 0.24 ! adjustable parameter
!         beta = 18.1  ! adjustable parameter

         ke = sr**(0.5*(1.0+vf_om_s-BA_alpha*vf_sand_s-vf_gravels_s)) &
               * ((1.0/(1.0+exp(-BA_beta*sr)))**3-((1.0-sr)/2.0)**3)**(1.0-vf_om_s)
      else                         ! Fozen or partially frozen soils
         ke = sr**(1.0+vf_om_s)
      endif
#endif

#if(defined THERMAL_CONDUCTIVITY_SCHEME_5)
! -----------------------------------------------------------------------------------------
! [5] Lu et al., 2007: An improved model for predicting soil thermal conductivity from
!     water content at room temperature. Soil Sci. Soc. Am. J. 71:8-14
! -----------------------------------------------------------------------------------------
      if(a > 0.4)then ! Coarse-textured soils = soils with sand fractions >40 (%)
         alpha = 0.728
         beta = 1.165
      else ! Fine-textured soils = soils with sand fractions <40 (%)
         alpha = 0.37
         beta = 1.29
      endif

      if(temperature > tfrz)then ! Unfrozen soils
         ke = exp(alpha*(1.0-sr**(alpha-beta)))
      else                         ! Fozen or partially frozen soils
         ke = sr
      endif
#endif
else
      ke = 0.0
endif
#if(defined THERMAL_CONDUCTIVITY_SCHEME_1 || defined THERMAL_CONDUCTIVITY_SCHEME_2 || defined THERMAL_CONDUCTIVITY_SCHEME_3 || defined THERMAL_CONDUCTIVITY_SCHEME_4 || defined THERMAL_CONDUCTIVITY_SCHEME_5)
      ke = max(ke, 0.0)
      ke = min(ke, 1.0)
      if(temperature > tfrz)then ! Unfrozen soil
         thk = (ksat_u-kdry)*ke + kdry
      else                         ! Frozen or partially frozen soils
         thk = (ksat_f-kdry)*ke + kdry
      endif
#endif

#if(defined THERMAL_CONDUCTIVITY_SCHEME_6)
! -----------------------------------------------------------------------------------------
! [6] Series-Parallel Models (Tarnawski and Leong, 2012)
! -----------------------------------------------------------------------------------------
      a = wf_gravels_s+wf_sand_s

! a fitting parameter of the soil solid uniform passage
      aa = 0.0237 - 0.0175*a**3

! a fitting parameter of a minuscule portion of soil water (nw) plus a minuscule portion of soil air (na)
      nwm = 0.088 - 0.037*a**3

! the degree of saturation of the minuscle pore space
      x = 0.6 - 0.3*a**3
      if(sr < 1.0e-6)then
         nw_nwm = 0.0
      else
         nw_nwm = exp(1.0-sr**(-x))
      endif

      if(temperature > tfrz)then ! Unfrozen soil
         thk = k_solids*aa + (1.0-vf_pores_s-aa+nwm)**2 &
                / ((1.0-vf_pores_s-aa)/k_solids+nwm/(k_water*nw_nwm+k_air*(1.0-nw_nwm))) &
                + k_water*(vf_pores_s*sr-nwm*nw_nwm) &
                + k_air*(vf_pores_s*(1.0-sr)-nwm*(1.0-nw_nwm))
      else
         thk = k_solids*aa + (1.0-vf_pores_s-aa+nwm)**2 &
                / ((1.0-vf_pores_s-aa)/k_solids+nwm/(k_ice*nw_nwm+k_air*(1.0-nw_nwm))) &
                + k_ice*(vf_pores_s*sr-nwm*nw_nwm) &
                + k_air*(vf_pores_s*(1.0-sr)-nwm*(1.0-nw_nwm))
      endif
#endif

#if(defined THERMAL_CONDUCTIVITY_SCHEME_7)
! -----------------------------------------------------------------------------------------
! [7] Thermal properties of soils, in Physics of Plant Environment, 
!     ed. by W.R. van Wijk (North-Holland, Amsterdam, 1963), pp. 210-235
! -----------------------------------------------------------------------------------------
      if(sr*vf_pores_s <= 0.09)then
         ga = 0.013+0.944*sr*vf_pores_s
      else
         ga = 0.333 - (1.-sr)*vf_pores_s/vf_pores_s*(0.333-0.035)
      endif
         gc = 1.0-2.0*ga

      if(temperature > tfrz)then ! Unfrozen soil
         aa = (2.0/(1.0+(k_air/k_water-1.0)*ga) &    ! the shape factor
            +  1.0/(1.0+(k_air/k_water-1.0)*gc))/3.0
         aaa = (2.0/(1.0+(k_solids/k_water-1.0)*0.125) &    ! the shape factor
            +  1.0/(1.0+(k_solids/k_water-1.0)*(1.0-2.0*0.125)))/3.0

         thk = (sr*vf_pores_s*k_water + (1.-sr)*vf_pores_s*aa*k_air + (1.-vf_pores_s)*aaa*k_solids) &
                / (sr*vf_pores_s + (1.-sr)*vf_pores_s*aa + (1.-vf_pores_s)*aaa)
      else
         aa = (2.0/(1.0+(k_air/k_ice-1.0)*ga) &    ! the shape factor
            +  1.0/(1.0+(k_air/k_ice-1.0)*gc))/3.0
         aaa = (2.0/(1.0+(k_solids/k_ice-1.0)*0.125) &    ! the shape factor
            +  1.0/(1.0+(k_solids/k_ice-1.0)*(1.0-2.0*0.125)))/3.0

         thk = (sr*vf_pores_s*k_ice + (1.-sr)*vf_pores_s*aa*k_air + (1.-vf_pores_s)*aaa*k_solids) &
                / (sr*vf_pores_s + (1.-sr)*vf_pores_s*aa + (1.-vf_pores_s)*aaa)
      endif
#endif

END SUBROUTINE soil_hcap_cond
