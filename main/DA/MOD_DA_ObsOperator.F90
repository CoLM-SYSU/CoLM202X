#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DA_ObsOperator
!-----------------------------------------------------------------------
! DESCRIPTION:
!    Forward modeling of brightness temperature observations
!
! AUTHOR:
!   Lu Li, 12/2024: Initial version
!   Zhilong Fan, Lu Li, 03/2024: Debug and clean codes
!-----------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Const_Physical
   USE MOD_Vars_1DForcing
   USE MOD_DA_Const
   USE MOD_Vars_Global, only: nl_soil, nl_lake, N_land_classification
   USE MOD_Namelist, only: DEF_DA_CMEM
   IMPLICIT NONE
   SAVE

! public functions
   PUBLIC   :: forward

! local variables
   real(r8) :: fghz                       ! frequency of satellite (GHz)
   real(r8) :: theta                      ! incidence angle of satellite (rad)
   real(r8) :: f                          ! frequency (Hz)
   real(r8) :: omega                      ! radian frequency
   real(r8) :: lam                        ! wavelength (m)
   real(r8) :: k                          ! wave number (rad/m)
   real(r8) :: kcm                        ! wave number (rad/cm)
   real(r8) :: kr                         ! size parameter used in calcuate single-particle albedo
   real(r8) :: hr(N_land_classification)  ! roughness parameter define from CMEM, only support IGBP

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE forward( &
      patchtype, patchclass, dz_sno, &
      forc_topo, htop, &
      tref, t_soisno, tleaf, &
      wliq_soisno, wice_soisno, h2osoi, &
      snowdp, lai, sai, &
      vf_clay, vf_sand, BD_all, porsl, &
      sat_theta, sat_fghz, &
      tb_toa_h, tb_toa_v)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Forward modeling of brightness temperature observations
!
! REFERENCES:
!   [1] Wigneron et al., 2007, "L-band Microwave Emission of the Biosphere (L-MEB) Model:
!       Description and calibration against experimental
!       data sets over crop fields" Remote Sensing of Environment. Vol. 107, pp. 639-655k
!
!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Const_Physical
      USE MOD_Vars_Global, only: nl_soil, nl_lake, maxsnl, spval, dz_soi
      USE MOD_DA_Const
      IMPLICIT NONE

!------------------------ Dummy Argument ------------------------------
      integer, intent(in)   :: patchtype                        ! land cover type
      integer, intent(in)   :: patchclass                       ! land cover class
      real(r8), intent(in)  :: dz_sno(maxsnl + 1:0)             ! layer thickness (m)
      real(r8), intent(in)  :: forc_topo                        ! topography [m]
      real(r8), intent(in)  :: htop                             ! upper height of vegetation [m]
      real(r8), intent(in)  :: tref                             ! 2 m height air temperature [kelvin]
      real(r8), intent(in)  :: tleaf                            ! leaf temperature [K]
      real(r8), intent(in)  :: t_soisno(maxsnl + 1:nl_soil)     ! soil temperature [K]
      real(r8), intent(in)  :: wliq_soisno(maxsnl + 1:nl_soil)  ! liquid water in layers [kg/m2]
      real(r8), intent(in)  :: wice_soisno(maxsnl + 1:nl_soil)  ! ice lens in layers [kg/m2]
      real(r8), intent(in)  :: h2osoi(nl_soil)                  ! volumetric soil water in layers [m3/m3]
      real(r8), intent(in)  :: snowdp                           ! snow depth [meter]
      real(r8), intent(in)  :: lai                              ! leaf area index
      real(r8), intent(in)  :: sai                              ! stem area index
      real(r8), intent(in)  :: vf_clay(nl_soil)                 ! volumetric fraction of clay
      real(r8), intent(in)  :: vf_sand(nl_soil)                 ! volumetric fraction of sand
      real(r8), intent(in)  :: BD_all(nl_soil)                  ! bulk density of soil (GRAVELS + ORGANIC MATTER + Mineral Soils,kg/m3)
      real(r8), intent(in)  :: porsl(nl_soil)                   ! fraction of soil that is voids [-]
      real(r8), intent(in)  :: sat_theta                        ! incidence angle of satellite (rad)
      real(r8), intent(in)  :: sat_fghz                         ! frequency of satellite (GHz)
      real(r8), intent(out) :: tb_toa_h                         ! brightness temperature of top-of-atmosphere for H- polarization
      real(r8), intent(out) :: tb_toa_v                         ! brightness temperature of top-of-atmosphere for V- polarization

!----------------------- Local Variables -------------------------------
      logical  :: is_low_veg                  ! flag for low vegetation
      real(r8) :: dz_soisno(maxsnl+1:nl_soil) ! liquid water in layers [kg/m2]
      integer  :: lb                          ! lower bound of arrays
      real(r8) :: tau_atm                     ! atmospheric optical depth
      real(r8) :: r_r(2)                      ! rough surface reflectivity for H and V polarizations
      real(r8) :: r_sn(2)                     ! reflectivity between the snow and ground for H and V polarizations
      real(r8) :: r_snow(2)                   ! reflectivity of the snow for H and V polarizations
      real(r8) :: tb_soil(2)                  ! brightness temperature of soil for H and V polarizations
      real(r8) :: tb_tos(2)                   ! brightness temperature of snow-covered ground for H and V polarizations
      real(r8) :: tb_tov(2)                   ! brightness temperature of vegetation (consider snow) for H and V polarizations
      real(r8) :: tb_tov_noad(2)              ! brightness temperature of vegetation (no downwelling radiation) for H and V polarizations
      real(r8) :: tb_au(2)                    ! upwelling radiation (brightness temperature) of atmosphere
      real(r8) :: tb_ad(2)                    ! downwelling radiation (brightness temperature) of atmosphere
      real(r8) :: rho_snow                    ! snow density (g/cm3)
      real(r8) :: liq_snow                    ! snow liquid water content (cm3/cm3)
      real(r8) :: gamma_p(2)                  ! vegetation opacity for H- and V- polarization
      real(r8) :: tb_veg(2)                   ! brightness temperature of vegetation for H- and V- polarization
      real(r8) :: tb_2(2)                     ! the downwelling vegetation emission reflected by the soil and attenuated by the canopy layer
      real(r8) :: tb_3(2)                     ! upwelling soil emission attenuated by the canopy
      real(r8) :: tb_4(2)                     ! the downwelling cosmic ray reflected by the soil and attenuated by the canopy layer
      real(r8) :: tb_toa(2)                   ! brightness temperature of top-of-atmosphere for H- and V- polarization
      integer  :: i

!-----------------------------------------------------------------------

      ! calculate depth of soil and snow layers
      dz_soisno(maxsnl+1:0) = dz_sno(maxsnl+1:0)
      dz_soisno(1:nl_soil) = dz_soi(1:nl_soil)

      ! calculate lower bound of snow
      lb = 0
      DO i = 0, maxsnl
         IF (wliq_soisno(i) < 0.0) THEN
            lb = i+1
            EXIT
         ENDIF
      ENDDO

      ! check the patch type
      IF (patchtype >= 3) THEN ! ocean, lake, ice
         tb_toa = spval
      ELSE
         ! calculate parameters used in operator varied with satellite
         CALL calc_parameters(sat_theta, sat_fghz)

!#############################################################################
! atmosphere module
!#############################################################################
         CALL atm(forc_topo, tref, tau_atm, tb_au, tb_ad)

!#############################################################################
! soil module
!#############################################################################
         CALL soil( &
            patchtype, patchclass, dz_soisno(1:nl_soil), &
            t_soisno(1:nl_soil), h2osoi, wliq_soisno(1:nl_soil), &
            vf_sand, vf_clay, BD_all, porsl, &
            r_r, tb_soil)

!#############################################################################
! vegetation and snow module
!    We categorized four different cases for the calculations: 
!    1) no vegetation and no snow
!    2) no vegetation with snow
!    3) vegetation without snow
!    4) vegetation with snow
!#############################################################################
         ! roughly judge low or high vegetation (only for IGBP)
         is_low_veg = .true.
         IF (patchclass >= 1 .and. patchclass <= 5) THEN
            is_low_veg = .false.
         END IF

         ! ensure snow density to <= 1 g/cm3
         IF (snowdp > 0.01) THEN
            rho_snow = (wliq_soisno(lb) + wice_soisno(lb))/(dz_soisno(lb)*1e3)
            liq_snow = wliq_soisno(lb)*rho_snow/(dz_soisno(lb)*1e3)

            IF (liq_snow > 1.0 .or. rho_snow > 1.0) then
               rho_snow = 1.0
               liq_snow = wliq_soisno(lb)*rho_snow/(dz_soisno(lb)*1e3)
            END IF
         END IF

         ! main procedures
         ! --------------------------------------------------------------------
         ! 1) no veg and no snow
         !    two components: 
         !       (1) brightness temperature of soil
         !       (2) the downwelling radiation reflected by the soil
         ! --------------------------------------------------------------------
         IF ((lai + sai < 1e-6) .and. (snowdp < 0.01)) THEN    
            tb_tov = tb_soil + tb_ad*r_r
            tb_tov_noad = tb_soil

         ! --------------------------------------------------------------------
         ! 2) no veg and has snow
         !    two components: 
         !       (1) brightness temperature of snow
         !       (2) the downwelling radiation reflected by the snow
         ! --------------------------------------------------------------------
         ELSE IF ((lai + sai < 1e-6) .and. (snowdp > 0.01)) THEN  
            ! calculate brightness temperature of snow-covered ground
            CALL snow(t_soisno(1), t_soisno(1), snowdp, rho_snow, liq_snow, r_r, r_snow, tb_tos)

            tb_tov = tb_tos + tb_ad*r_snow
            tb_tov_noad = tb_tos

         ! --------------------------------------------------------------------
         ! 3) has veg and no snow
         !    four components: 
         !       (1) the direct upwelling vegetation emission,
         !       (2) the downwelling vegetation emission reflected by the soil and attenuated by the canopy layer
         !       (3) upwelling soil emission attenuated by the canopy
         !       (4) the downwelling reflected by the soil and attenuated by the canopy layer
         ! --------------------------------------------------------------------         
         ELSE IF ((lai + sai > 1e-6) .and. (snowdp < 0.01)) THEN 
            ! calculate brightness temperature of vegetation
            CALL veg(patchclass, lai, htop, 0.0, t_soisno(1), tb_veg, gamma_p)

            DO i = 1, 2
               tb_2(i) = tb_veg(i)*gamma_p(i)*r_r(i)
               tb_3(i) = tb_soil(i)*gamma_p(i)
               tb_4(i) = tb_ad(i)*r_r(i)*(gamma_p(i)**2)
               tb_tov(i) = tb_veg(i) + tb_2(i) + tb_3(i) + tb_4(i)
               tb_tov_noad(i) = tb_veg(i) + tb_2(i) + tb_3(i)
            END DO

         ! --------------------------------------------------------------------
         ! 4) has veg and has snow
         !    We need to determine the positional relationship between vegetation and snow.
         !    
         !    If vegetation is higher than snow, 
         !    we first calculate brightness temperature of snow (soil boundary), then calculate
         !    four components to derive brightness temperature of top of vegetation: 
         !       (1) the direct upwelling vegetation emission,
         !       (2) the downwelling vegetation emission reflected by the snow and attenuated by the canopy layer
         !       (3) upwelling snow emission attenuated by the canopy
         !       (4) the downwelling reflected by the snow and attenuated by the canopy layer
         !
         !    If vegetation is lower than snow
         !    we first calculate brightness temperature of top of vegetation (soil boundary), then calculate
         !    four components to derive brightness temperature of top of snow:
         !       (1) the direct upwelling vegetation emission,
         !       (2) the downwelling vegetation emission reflected by the soil and attenuated by the canopy layer
         !       (3) upwelling soil emission attenuated by the canopy
         !       (4) the downwelling reflected by the soil and attenuated by the canopy layer      
         ! --------------------------------------------------------------------     
         ELSE IF ((lai + sai > 1e-6) .and. (snowdp > 0.01)) THEN
            IF (htop < snowdp) THEN
               ! calculate brightness temperature of low vegetation
               CALL veg(patchclass, lai, htop, snowdp, t_soisno(1), tb_veg, gamma_p)

               ! calculate brightness temperature of top of vegetation
               DO i = 1, 2
                  tb_2(i) = tb_veg(i)*gamma_p(i)*r_r(i)
                  tb_3(i) = tb_soil(i)*gamma_p(i)
                  tb_4(i) = tb_ad(i)*r_r(i)*(gamma_p(i)**2)
                  tb_tov(i) = tb_veg(i) + tb_2(i) + tb_3(i) + tb_4(i)
                  tb_tov_noad(i) = tb_veg(i) + tb_2(i) + tb_3(i)
               END DO

               ! calculate reflectivity between the snow and low veg (adopted from CMEM)
               r_sn(:) = 1.0 - tb_tov_noad(:)/t_soisno(1)

               ! calculate brightness temperature of snow-covered ground
               CALL snow(t_soisno(1), t_soisno(1), snowdp, rho_snow, liq_snow, r_sn, r_snow, tb_tos)

               ! calculate brightness temperature of top of snow
               tb_tov = tb_tos + tb_ad*r_snow
               tb_tov_noad = tb_tos

            ELSE
               ! calculate brightness temperature of snow-covered ground
               CALL snow(t_soisno(1), t_soisno(1), snowdp, rho_snow, liq_snow, r_r, r_snow, tb_tos)

               ! calculate brightness temperature of high vegetation
               CALL veg(patchclass, lai, htop, snowdp, t_soisno(1), tb_veg, gamma_p)

               ! calculate brightness temperature of top of vegetation
               DO i = 1, 2
                  tb_2(i) = tb_veg(i)*gamma_p(i)*r_snow(i)
                  tb_3(i) = tb_tos(i)*gamma_p(i)
                  tb_4(i) = tb_ad(i)*r_snow(i)*(gamma_p(i)**2)
                  tb_tov(i) = tb_veg(i) + tb_2(i) + tb_3(i) + tb_4(i)
                  tb_tov_noad(i) = tb_veg(i) + tb_2(i) + tb_3(i)
               END DO
            END IF
         END IF

!#############################################################################
! Caculate brightness temperature of top-of-atmosphere
!#############################################################################
         tb_toa = tb_tov*exp(-tau_atm) + tb_au

      END IF

      tb_toa_h = tb_toa(1)
      tb_toa_v = tb_toa(2)

   END SUBROUTINE forward


!-----------------------------------------------------------------------

   SUBROUTINE calc_parameters (sat_theta, sat_fghz)

!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_DA_Const
      IMPLICIT NONE

!------------------------ Dummy Argument -------------------------------
      real(r8), intent(in) :: sat_theta, sat_fghz

!-----------------------------------------------------------------------

      theta = sat_theta                   ! incidence angle of satellite (rad)
      fghz = sat_fghz                     ! frequency of satellite (GHz)   
      f = fghz*1e9                        ! frequency (Hz)
      omega = 2.0d0*pi*f                  ! radian frequency (rad/s)
      lam = C/f                           ! wavelength (m)
      k = 2*pi/lam                        ! wave number (rad/m)
      kcm = k*100.0                       ! wave number (rad/cm)
      kr = k*(0.5*1e-3)                   ! size parameter used in calcuate single-particle albedo
      hr(:) = (2.0d0*kcm*rgh_surf)**2.0d0 ! roughness parameter define from CMEM, only support IGBP

   END SUBROUTINE calc_parameters

!-----------------------------------------------------------------------

   SUBROUTINE atm(z, tref, tau_atm, tb_au, tb_ad)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the atmospheric opacity and up and downwelling brightness temperature
!
! REFERENCES:
!   [1] Pellarin, T., et al. (2003), Two-year global simulation of L-band brightness
!       temperature over land, IEEE Trans. Geosci. Remote Sens., 41, 2135–2139.
!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Const_Physical
      IMPLICIT NONE

!------------------------ Dummy Argument ------------------------------
      real(r8), intent(in)  :: z            ! altitude (m)
      real(r8), intent(in)  :: tref         ! 2m air temperature (K)
      real(r8), intent(out) :: tau_atm      ! atmospheric optical depth
      real(r8), intent(out) :: tb_au(2)     ! upwelling radiation (brightness temperature) of atmosphere
      real(r8), intent(out) :: tb_ad(2)     ! downwelling radiation (brightness temperature) of atmosphere

!----------------------- Local Variables -------------------------------
      real(r8) :: t_sky = 2.7    ! cosmic ray radiation (K)
      real(r8) :: t_eq           ! equivalent layer temperature
      real(r8) :: gossat

!-----------------------------------------------------------------------

      ! calculate optical depth of atmosphere [1] eq(A1)
      tau_atm = exp(-3.9262 - 0.2211*z/1000 - 0.00369*tref)/cos(theta)
      gossat = exp(-tau_atm)

      ! calculate equivalent layer temperature
      t_eq = exp(4.9274 + 0.002195*tref)

      ! upwelling radiation (brightness temperature) of atmosphere
      tb_au(:) = t_eq*(1.-gossat)

      ! downwelling radiation (brightness temperature) of atmosphere [1] eq(A2)
      tb_ad(:) = t_eq*(1.-gossat) + t_sky*gossat

   END SUBROUTINE atm

!-----------------------------------------------------------------------

   SUBROUTINE soil(patchtype, patchclass, dz_soisno, t_soisno, h2osoi, wliq_soisno, vf_sand, vf_clay, BD_all, porsl, r_r, tb_soil)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate brightness temperature of soil surface
! REFERENCES:
!   [1] Wigneron et al., 2007, "L-band Microwave Emission of the Biosphere (L-MEB) Model:
!       Description and calibration against experimental
!       data sets over crop fields" Remote Sensing of Environment. Vol. 107, pp. 639-655k
!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Const_Physical
      IMPLICIT NONE

!------------------------ Dummy Argument ------------------------------
      integer, intent(in)   :: patchtype                ! land cover type
      integer, intent(in)   :: patchclass               ! land cover class
      real(r8), intent(in)  :: dz_soisno(1:nl_soil)     ! soil layer thickness (m)
      real(r8), intent(in)  :: t_soisno(1:nl_soil)      ! soil temperature (K)
      real(r8), intent(in)  :: h2osoi(1:nl_soil)        ! soil water content (m3/m3)
      real(r8), intent(in)  :: wliq_soisno(1:nl_soil)   ! liquid water in layers [kg/m2]
      real(r8), intent(in)  :: vf_sand(nl_soil)         ! sand volume fraction
      real(r8), intent(in)  :: vf_clay(nl_soil)         ! clay volume fraction
      real(r8), intent(in)  :: BD_all(nl_soil)          ! bulk density of soil (kg/m3)
      real(r8), intent(in)  :: porsl(nl_soil)           ! soil porosity
      real(r8), intent(out) :: r_r(2)                   ! rough surface reflectivity for H and V polarizations
      real(r8), intent(out) :: tb_soil(2)               ! brightness temperature of soil

!----------------------- Local Variables -------------------------------
      real(r8)    :: liq_h2osoi(1:nl_soil)      ! liquid volumetric water content
      logical     :: is_desert                  ! flag for desert
      real(r8)    :: t_eff(2)                   ! effective temperature for H and V polarizations, [K]
      complex(r8) :: eps_soil                   ! dielectric constant of soil for H and V polarizations
      real(r8)    :: r_s(2)                     ! smooth surface reflectivity for H and V polarizations

!-----------------------------------------------------------------------

      ! liquid volumetric water content
      liq_h2osoi = wliq_soisno/(dz_soisno*denh2o)

      ! whether this patch is desert
      is_desert = .FALSE.
      IF (liq_h2osoi(1) < 0.02 .and. vf_sand(1) > 90) THEN
         is_desert = .TRUE.
      END IF

      ! caculate effective temperature [1](11)(12)
      CALL eff_soil_temp(h2osoi, t_soisno, t_eff)

      ! caculate soil dielectric constant
      CALL diel_soil(is_desert, t_soisno, h2osoi, liq_h2osoi, vf_sand, vf_clay, BD_all, porsl, eps_soil)

      ! caculate smooth surface reflectivity [1](3)
      CALL smooth_reflectivity(eps_soil, r_s)

      ! caculate rough surface reflectivity [1](7)
      CALL rough_reflectivity(is_desert, patchclass, r_s, r_r)

      ! calculate brightness temperature [1](2)
      IF (is_desert) THEN
         CALL desert(t_eff, r_r, eps_soil, tb_soil)
      ELSE
         tb_soil = t_eff*(1 - r_r)
      END IF

   END SUBROUTINE soil

!-----------------------------------------------------------------------

   SUBROUTINE desert(t_soil, r_r, eps, tb_desert)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate desert emissivity using Grody and Weng, 2008
!
! REFERENCES:
!  [1] Grody, N. C., & Weng, F. (2008). Microwave emission and scattering from deserts:
!      Theory compared with satellite measurements.
!      IEEE Transactions on Geoscience and Remote Sensing, 46, 361–375.
!-----------------------------------------------------------------------

      USE MOD_Precision
      USE MOD_Const_Physical
      IMPLICIT NONE

!------------------------ Dummy Argument ------------------------------
      real(r8), intent(in)    :: t_soil(2)     ! desert surface temperature (K) (h-pol and v-pol)
      real(r8), intent(in)    :: r_r(2)        ! reflectivity of rough surface  (h-pol and v-pol)
      complex(r8), intent(in) :: eps           ! diel. const of desert
      real(r8), intent(out)   :: tb_desert(2)  ! brightness temperature of soil for H- and V- polarization

!----------------------- Local Variables -------------------------------
      real(r8) :: f0 = 0.7                      ! the fractional volume of spherical particles
                                                ! (f = (4/3)*pi*r^3*n0),
                                                ! r : the particle radius = 0.5 (mm)
                                                ! n0 : the number of particles per unit volume.
      real(r8) :: w     ! single-particle albedo
      real(r8) :: g     ! asymmetry parameter
      real(r8) :: a     ! similarity parameter
      real(r8) :: em(2) ! desert soil emissivity
      real(r8) :: y_r   ! real part of y-parameters
      real(r8) :: y_i   ! imaginary part of y-parameters

!-----------------------------------------------------------------------
      ! calculate y-parameters (eq.A15)
      y_r = (real(eps) - 1)/(real(eps) + 2)                              ! [1](A15)
      y_i = 3*aimag(eps)/(real(eps) + 2)**2                              ! [1](A15)

      ! calculate single-particle albedo (eq.A16)
      w = (1 - f0)**4*kr**3*y_r**2/ &
         ((1 - f0)**4*kr**3*y_r**2 + 1.5*(1 + 2*f0)**2*y_i)              ! [1](A16)

      ! calculate asymmetry parameter (p.374)
      g = 0.23*kr**2                                                     ! [1]p.374

      ! calculate similarity parameter (eq.3b)
      a = sqrt((1 - w)/(1 - w*g))                                        ! [1](3b)

      ! calculate desert soil emissivity (eq.A13)
      em = (1 - r_r)*(2*a/((1 + a) - (1 - a)*r_r))                       ! [1](13)

      ! calculate brightness temperature of desert
      tb_desert = t_soil*em

   END SUBROUTINE desert

!-----------------------------------------------------------------------

   SUBROUTINE veg(patchclass, lai, htop, snowdp, tleaf, tb_veg, gamma_p)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the brightness temperature on the top of vegetation
!
! REFERENCES:
!   [1] Wigneron et al., 2007, "L-band Microwave Emission of the Biosphere (L-MEB) Model:
!       Description and calibration against experimental
!       data sets over crop fields" Remote Sensing of Environment. Vol. 107, pp. 639-655k
!
!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Const_Physical
      IMPLICIT NONE

!------------------------ Dummy Argument ------------------------------
      integer, intent(in)   :: patchclass     ! land cover class
      real(r8), intent(in)  :: lai            ! leaf area index
      real(r8), intent(in)  :: tleaf          ! leaf temperature (K)
      real(r8), intent(out) :: tb_veg(2)      ! brightness temperature of vegetation for H- and V- polarization
      real(r8), intent(out) :: gamma_p(2)     ! vegetation opacity for H- and V- polarization
      real(r8), intent(in)  :: htop, snowdp

!----------------------- Local Variables -------------------------------
      real(r8) :: tau_nadir    ! vegetation opacity at nadir
      real(r8) :: tau_veg(2)   ! vegetation opacity for H- and V- polarization
      integer  :: i

!-----------------------------------------------------------------------

      ! caculate vegetation opacity (optical depth) at nadir b*VWC
      IF (htop < snowdp) THEN
         tau_nadir = b1(patchclass)*lai + b2(patchclass) ! low veg              ! [1](22)
      ELSE
         tau_nadir = b3(patchclass) ! high veg
      END IF

      ! calculate vegetation optical depth at H- and V- polarizations
      tau_veg(1) = tau_nadir*(cos(theta)**2 + tth(patchclass)*sin(theta)**2)    ! [1](23)
      tau_veg(2) = tau_nadir*(cos(theta)**2 + ttv(patchclass)*sin(theta)**2)    ! [1](24)

      ! calculate brightness temperature of vegetation
      DO i = 1, 2
         gamma_p(i) = exp(-tau_veg(i)/cos(theta))                               !  [1](15)
         IF (DEF_DA_CMEM) THEN
            tb_veg(i) = (1.-w_CMEM(patchclass))*(1.-gamma_p(i))*tleaf
         ELSE
            tb_veg(i) = (1.-w_Wigneron(patchclass))*(1.-gamma_p(i))*tleaf
         END IF
      END DO

   END SUBROUTINE veg

!-----------------------------------------------------------------------

   SUBROUTINE snow(t_snow, t, snowdp, rho_snow, liq_snow, r_sn, r_snow, tb_tos)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the brightness temperature of snow-covered ground
!
! REFERENCES:
!   [1] Christian Mätzler (1987) Applications of the interaction of
!       microwaves with the natural snow cover, Remote Sensing Reviews, 2:2, 259-387, DOI:
!       10.1080/02757258709532086
!
!   [2] Anderson, E. A., 1976: A point energy and mass balance model of a snow cover.
!       NOAA Tech. Rep. NWS 19, 150 pp. U.S. Dept. of Commer., Washington, D.C.(eq.5.1)
!
!   [3] Hallikainen, M. T., F. Ulaby, and T. Deventer. 1987. Extinction behavior of dry snow in the
!       18- to 90-GHz range. IEEE Trans. Geosci. Remote Sens., GE-25, 737–745.
!
!   [4] Microwave remote sensing : active and passive
!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Const_Physical
      IMPLICIT NONE

!------------------------ Dummy Argument ------------------------------
      real(r8), intent(in)  :: t_snow       ! average snow temperature (K)
      real(r8), intent(in)  :: t            ! temperature at bottom of snow (K), i.e., soil or leaf
      real(r8), intent(in)  :: snowdp       ! snow depth (m)
      real(r8), intent(in)  :: rho_snow     ! snow density (g/cm3)
      real(r8), intent(in)  :: liq_snow     ! snow liquid water content (cm3/cm3)
      real(r8), intent(in)  :: r_sn(2)      ! reflectivity between the snow and ground at (1, H-POL. 2, V.)
      real(r8), intent(out) :: r_snow(2)    ! reflectivity between the snow and air for H- and V- polarization
      real(r8), intent(out) :: tb_tos(2)    ! brightness temperature of snow-cover ground for H- and V- polarization

!----------------------- Local Variables -------------------------------
      real(r8) :: sal_snow = 0.0                   ! snow salinity (pmm)
      real(r8) :: eps_i_r                          ! real part of dielectric constant of ice
      real(r8) :: eps_i_i                          ! imaginary part of dielectric constant of ice
      real(r8) :: eps_i_is                         ! imaginary part of dielectric constant of impure ice -5(C)
      real(r8) :: eps_i_ip                         ! imaginary part of dielectric constant of pure ice -5(C)
      real(r8) :: eps_ds_r                         ! real part of dielectric constant of dry snow
      real(r8) :: eps_ds_i                         ! imaginary part of dielectric constant of dry snow
      real(r8) :: eps_ws_i                         ! imaginary part of dielectric constant of wet snow
      real(r8) :: eps_ws_r                         ! real part of dielectric constant of wet snow
      real(r8) :: eps_w_s = 88.                    ! dielectric constant of static water
      real(r8) :: eps_a_inf, eps_b_inf, eps_c_inf  ! infinite frequency dielectric constant of three parts
      real(r8) :: eps_a_s, eps_b_s, eps_c_s        ! static dielectric constant of three parts
      complex(r8) :: eps_i                         ! dielectric constant of ice
      complex(r8) :: eps_a
      complex(r8) :: eps_b
      complex(r8) :: eps_c                         ! dielectric constant of three parts
      complex(r8) :: eps                           ! dielectric constant of wet snow
      real(r8) :: rho_ds                           ! density of dry snow  (g/cm3)
      real(r8) :: rho_i = 0.916                    ! density of ice (g/cm3)
      real(r8) :: aa = 0.005
      real(r8) :: bb = 0.4975
      real(r8) :: cc = 0.4975                      ! fitting parameters
      real(r8) :: fa, fb, fc                       ! relaxation frequency of wet snow
      real(r8) :: d                                ! snow grain size (mm)
      real(r8) :: alpha, beta, pp, qq              ! parameter used to calculate propogation angle in snow
      real(r8) :: theta_s                          ! propogation angle in snow
      complex(r8) :: z_s                           ! wave impedance in snow
      real(r8) :: r_sa                             ! reflectivity between the snow and air for H- and V- polarization
      real(r8) :: tb_2                             ! the net apparent temperature contributions due to emission by layers 2 (snow)
      real(r8) :: tb_3                             ! the net apparent temperature contributions due to emission by layers 3 (soil)
      real(r8) :: l2_apu                               ! extinction coefficient of snow (Beer's Law)
      real(r8) :: q = 0.96                         ! parameter
      real(r8) :: ka_ws                            ! absorption coefficient of wet snow
      real(r8) :: ka_ds                            ! absorption coefficient of dry snow
      real(r8) :: ke                               ! extinction coefficient of wet snow
      real(r8) :: ke_ds                            ! extinction coefficient of dry snow
      real(r8) :: ks                               ! scattering coefficient of snow
      real(r8) :: b_ds, b_ws
      real(r8) :: wk_h                             ! [m], equal to cmem cmem_snow_set_var.F90:476
      integer  :: i

!-----------------------------------------------------------------------
      IF (snowdp > 0.01) THEN  ! > 1cm
         ! calculate dielectric constant of ice
         CALL diel_ice(t_snow - tfrz, eps_i)
         eps_i_r = real(eps_i)
         eps_i_i = aimag(eps_i)

         ! consider the effect of salinity on the dielectric constant of ice
         eps_i_is = 0.0026/fghz + 0.00023*(fghz**0.87) ! impure ice -5(C)
         eps_i_ip = 6.e-4/fghz + 6.5e-5*(fghz**1.07)   ! pure ice -5(C)
         eps_i_i = eps_i_i + (eps_i_is - eps_i_ip)*sal_snow/13.0d0 ! corrected imaginary part of diel cons of ice

         ! calculate dielectric constant of dry snow (mixed by air and ice) (Polder–van Santen mixing model)
         rho_ds = (rho_snow - liq_snow)/(1.0 - liq_snow)   ! caculate density of dry snow,
         eps_ds_r = 1.0 + 1.58*rho_ds/(1.0 - 0.365*rho_ds)
         eps_ds_i = 3.0*(rho_ds/rho_i)*eps_i_i*(eps_ds_r**2)*(2*eps_ds_r + 1)/ &
            ((eps_i_r + 2*eps_ds_r)*(eps_i_r + 2*eps_ds_r**2))   ! negative imaginary part of diel cons of dry snow

         ! calculate dielectric constant of wet snow (Matzler 1987)
         IF (liq_snow > 0.) THEN  ! wet snow
            ! caculate relaxation frequency of three parts (eq.2.26)
            fa = f0w*(1 + (aa*(eps_w_s - eps_w_inf)/(eps_ds_r + (aa*(eps_w_inf - eps_ds_r)))))
            fb = f0w*(1 + (bb*(eps_w_s - eps_w_inf)/(eps_ds_r + (bb*(eps_w_inf - eps_ds_r)))))
            fc = f0w*(1 + (cc*(eps_w_s - eps_w_inf)/(eps_ds_r + (cc*(eps_w_inf - eps_ds_r)))))

            ! caculate infinite frequency dielectric constant of three parts
            eps_a_inf = (liq_snow*(eps_w_inf - eps_ds_r)/3.)/ &
               (1.+aa*((eps_w_inf/eps_ds_r) - 1.))
            eps_b_inf = (liq_snow*(eps_w_inf - eps_ds_r)/3.)/ &
               (1.+bb*((eps_w_inf/eps_ds_r) - 1.))
            eps_c_inf = (liq_snow*(eps_w_inf - eps_ds_r)/3.)/ &
               (1.+cc*((eps_w_inf/eps_ds_r) - 1.))

            ! caculate static dielectric constant of three parts
            eps_a_s = (liq_snow/3.)*(eps_w_s - eps_ds_r)/ &
               (1.+aa*((eps_w_s/eps_ds_r) - 1.))
            eps_b_s = (liq_snow/3.)*(eps_w_s - eps_ds_r)/ &
               (1.+bb*((eps_w_s/eps_ds_r) - 1.))
            eps_c_s = (liq_snow/3.)*(eps_w_s - eps_ds_r)/ &
               (1.+cc*((eps_w_s/eps_ds_r) - 1.))

            ! Debye equations
            eps_a = eps_a_inf + (eps_a_s - eps_a_inf)/(1 + jj*fghz/fa)
            eps_b = eps_b_inf + (eps_b_s - eps_b_inf)/(1 + jj*fghz/fb)
            eps_c = eps_c_inf + (eps_c_s - eps_c_inf)/(1 + jj*fghz/fc)

            ! calculate dielectric constant of wet snow
            eps = eps_a + eps_b + eps_c + (eps_ds_r - jj*eps_ds_i)
            eps_ws_r = real(eps)
            eps_ws_i = -1.*aimag(eps)
         ELSE
            eps_ws_r = eps_ds_r
            eps_ws_i = eps_ds_i
         END IF

         ! caculate propogation angle in snow (change medium from air to snow)
         alpha = k*abs(aimag(sqrt(eps_ws_r - jj*eps_ws_i)))
         beta = k*real(sqrt(eps_ws_r - jj*eps_ws_i))
         pp = 2.*alpha*beta
         qq = beta**2 - alpha**2 - (k*k)*(sin(theta)**2)
         theta_s = atan(k*sin(theta)/((1./sqrt(2.)) &
            *sqrt(sqrt(pp**2.+qq**2.) + qq)))

         ! caclulate wave impedance in snow
         z_s = z0/sqrt(eps_ws_r - jj*eps_ws_i)

         ! calculate brightness temperature above snow for H- V- polarization
         DO i = 1, 2
            ! Fresnel reflection coefficient between snow and air
            IF (i == 1) THEN
               r_sa = abs((z_s*cos(theta) - z0*cos(theta_s))/ &
                  (z_s*cos(theta) + z0*cos(theta_s)))**2
            ELSE
               r_sa = abs((z0*cos(theta) - z_s*cos(theta_s))/ &
                  (z0*cos(theta) + z_s*cos(theta_s)))**2
            END IF

            ! calculate snow grain size (mm) (Anderson 1976, eq.5.1) TODO
            d = min(1000*(1.6e-4 + 1.1e-13*((rho_snow*1000.0)**4)), 3.0)

            ! extinction coefficient of dry snow
            !//TODO: the paper is not focus on L-band, thus the formula is not suitable
            ke_ds = 0.0018*(fghz**2.8)*(d**2)/4.3429    !  [3](14)

            ! absorption coefficient of dry snow
            b_ds = (eps_ds_i/eps_ds_r)**2
            ka_ds = 2.*omega*sqrt(mu0*eps_0*eps_ds_r)* &
               sqrt(b_ds/(2.*(sqrt(1.+b_ds) + 1.)))

            IF (ke_ds < ka_ds) ke_ds = ka_ds

            ! absorption coefficient of wet snow
            b_ws = (eps_ws_i/eps_ws_r)**2
            ka_ws = 2.*omega*sqrt(mu0*eps_0*eps_ws_r)* &
               sqrt(b_ws/(2.*(sqrt(1.+b_ws) + 1)))

            ! total extinction (assuming scattering is the same for dry and wet snow)
            ke = (ke_ds - ka_ds) + ka_ws

            ! scattering coefficient of dry and wet snow
            ks = ke_ds - ka_ds

            wk_h = snowdp/rho_snow
            IF (((ke - q*ks)*(1./cos(theta_s))*wk_h) > (log(HUGE(1.)) - 1.)) THEN
               l2_apu = sqrt(HUGE(1.))
            ELSE
               l2_apu = min(sqrt(HUGE(1.)), exp((ke - q*ks)*(1./cos(theta_s))*wk_h))
            END IF

            ! brightness temperature through snow-air layer ([4] pp.243, eq.4.161)
            tb_2 = (1.+r_sn(i)/l2_apu)*(1.-r_sa)*t_snow* &
               (ka_ws/(ke - q*ks))*(1.-1./l2_apu)/ &
               (1.-r_sn(i)*r_sa/l2_apu**2.)

            ! brightness temperature through soil-snow layer ([4] pp.243, eq.4.162)
            tb_3 = ((1.-r_sn(i))*(1.-r_sa)*t)/ &
               (l2_apu*(1.-r_sn(i)*r_sa/l2_apu**2.))

            ! brightness temperature of snow-cover ground
            tb_tos(i) = tb_2 + tb_3

            ! emission by layers 2 (snow) and 3 (soil)
            ! reflectivity by layers 2 (snow) and 3 (soil)
            r_snow(i) = 1 - (tb_2/t_snow + tb_3/t)
         END DO
      END IF

   END SUBROUTINE snow

!-----------------------------------------------------------------------

   SUBROUTINE eff_soil_temp(wliq_soisno, t_soisno, t_eff)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the effective temperature of soil
!
! REFERENCES:
!   [1] Choudhury, B. J., Schmugge, T. J., & Mo, T. (1982). A parameterization
!       of effective soil temperature for microwave emission.
!       Journal of Geophysical Research: Oceans, 87(C2), 1301-1304.
!       https://doi.org/10.1029/JC087iC02p01301
!   [2] Wigneron, J.-P., L. Laguerre, and Y. Kerr (2001), A simple parmeterization
!       of the L-band microwave emission from rough agricultural soils, IEEE Trans.
!       Geosci. Remote. Sens., 39, 1697–1707.
!   [3] Wigneron et al., 2007, "L-band Microwave Emission of the Biosphere (L-MEB) Model:
!       Description and calibration against experimental
!       data sets over crop fields" Remote Sensing of Environment. Vol. 107, pp. 639-655k
!-----------------------------------------------------------------------
      USE MOD_Precision
      IMPLICIT NONE

! ------------------------ Dummy Argument ------------------------------
      real(r8), intent(in)  :: wliq_soisno(1:nl_soil)      ! soil moisture (m3/m3)
      real(r8), intent(in)  :: t_soisno(1:nl_soil)         ! soil temperature (K) or glacier temperature
      real(r8), intent(out) :: t_eff(2)                    ! effective temperature for H and V polarizations, [K]

!----------------------- Local Variables -------------------------------
      real(r8) :: wc_surf       ! soil moisture in the top 3cm (default in L-MEB)
      real(r8) :: t_deep        ! soil temperature at 50cm (default in L-MEB)
      real(r8) :: t_surf        ! soil temperature at 0-5cm (default in L-MEB)
      real(r8) :: C             ! parameter depending mainly on frequency
                                ! and soil moisture to describe the impact of
                                ! surface temperature on the effective temperature;
                                ! soil moisture increase, C large, teff close to tsurf
                                ! soil moisture decrease, C small, tdeep impact teff more
      real(r8) :: w0 = 0.41     ! parameter (f(texture))
      real(r8) :: bw = 0.35     ! parameter (f(texture)), the CMEM use 0.35, SMAP use 0.266
!-----------------------------------------------------------------------

      wc_surf = wliq_soisno(1)
      t_surf = t_soisno(1)
      t_deep = t_soisno(nl_soil)

      C = max(0.001, (wc_surf/w0)**bw)                !  [2](7)    [3](12)
      t_eff(:) = t_deep + (t_surf - t_deep)*C         !  [2](8)    [3](11)

   END SUBROUTINE eff_soil_temp

!-----------------------------------------------------------------------

   SUBROUTINE diel_soil(is_desert, t_soisno, h2osoi, liq_h2osoi, vf_sand, vf_clay, BD_all, porsl, eps)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the dielectric constant of a wet soil
!
! REFERENCES:
!   [1]  Matzler, C. (1998). Microwave permittivity of dry sand.
!   IEEE Transactions on Geoscience and Remote Sensing, 36(1), 317-319.
!
!   [2]  Wang and Schmugge, 1980: An empirical model for the
!   complex dielectric permittivity of soils as a function of water
!   content. IEEE Trans. Geosci. Rem. Sens., GE-18, No. 4, 288-295.
!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_Const_Physical
      IMPLICIT NONE

! ------------------------ Dummy Argument ------------------------------
      logical, intent(in)   :: is_desert              ! flag for desert soil
      real(r8), intent(in)  :: t_soisno(1:nl_soil)    ! soil temperature (K) or glacier temperature
      real(r8), intent(in)  :: h2osoi(1:nl_soil)      ! soil moisture (m3/m3)
      real(r8), intent(in)  :: liq_h2osoi(1:nl_soil)  ! liquid water content of soil (m3/m3)
      real(r8), intent(in)  :: vf_sand(nl_soil)       ! sand volume fraction
      real(r8), intent(in)  :: vf_clay(nl_soil)       ! clay volume fraction
      real(r8), intent(in)  :: BD_all(nl_soil)        ! bulk density of soil (kg/m3)
      real(r8), intent(in)  :: porsl(nl_soil)         ! soil porosity
      complex(r8), intent(out) :: eps

!----------------------- Local Variables -------------------------------
      integer  :: j
      real(r8) :: tc                       ! soil temperature (C)
      real(r8) :: wc                       ! soil moisture (m3/m3)
      real(r8) :: wp                       ! wilting point
      real(r8) :: ffrz                     ! fraction of frozen soil
      real(r8) :: wt                       ! transition moisture point (cm3/cm3)
      real(r8) :: gamma                    ! fitting parameter
      real(r8) :: ecl                      ! conductivity loss
      real(r8) :: alpha                    ! conductivity loss parameter
      real(r8) :: sal_sea = 32.5           ! sea water salinity (psu)
      real(r8) :: sal_soil = 0.0           ! soil salinity (psu)
      complex(r8) :: eps_x                 ! dielectric constant of the initially absorbed water
      complex(r8) :: eps_w                 ! dielectric constant of water
      complex(r8) :: eps_a = (1.0, 0.0)    ! dielectric constant of air,  [2]IV
      complex(r8) :: eps_r = (5.5, 0.2)    ! dielectric constant of rock, [2]IV
      complex(r8) :: eps_i = (3.2, 0.1)    ! dielectric constant of ice,  [2]IV
      complex(r8) :: eps_f = (5.0, 0.5)    ! dielectric constant of frozen soil
      real(r8)    :: porsl_loc, bd_all_loc

!-----------------------------------------------------------------------

      tc = t_soisno(1) - tfrz
      wc = h2osoi(1)

      IF (is_desert) THEN
         ! Microwave 1-10GHz permittivity of dry sand (matzler '98, eq.1)
         eps = 2.53 + (2.79 - 2.53)/(1 - jj*(fghz/0.27)) + jj*0.002    ! [1](1)
      ELSE
         ! calculate wilting point at the soil layer
         wp = 0.06774 - 0.00064*vf_sand(1) + 0.00478*vf_clay(1)        ! [2](1)

         ! define bulk density and porosity (from CoLM or CMEM)
         IF (DEF_DA_CMEM) THEN
            bd_all_loc = (vf_sand(1)*1.60d0 + vf_clay(1)*1.10d0 + (100.0d0 - vf_sand(1) - vf_clay(1))*1.20d0)/100.0d0
            porsl_loc = 1.0d0 - bd_all_loc/2.660d0
         ELSE
            bd_all_loc = BD_all(1)
            porsl_loc = porsl(1)
         END IF

         ! calculate fitting parameters
         gamma = -0.57*wp + 0.481                                      ! [2](8)

         ! calculate transition moisture point
         wt = 0.49*wp + 0.165                                          ! [2](9)

         IF (tc < -0.5) THEN ! assume all freeze
            ! calculate dielectric constant of pure ice
            CALL diel_ice(tc, eps_w)
         ELSE
            ! calculate dielectric constant of soil-water
            CALL diel_water(2, wc, tc, vf_sand(1), vf_clay(1), bd_all_loc, sal_soil, eps_w)
         END IF

         ! calculate dielectric constant of wet soil (when all soil freeze, eps_x = eps_i)
         ! (NOTE): use soil moisture (ice+liquid) to calculate dielectric constant
         IF (wc <= wt) THEN
            eps_x = eps_i + (eps_w - eps_i)*(wc/wt)*gamma                                        ! [2](3)
            eps = wc*eps_x + (porsl_loc - wc)*eps_a + (1.-porsl_loc)*eps_r                       ! [2](2)
         ELSE
            eps_x = eps_i + (eps_w - eps_i)*gamma                                                ! [2](5)
            eps = wt*eps_x + (wc - wt)*eps_w + (porsl_loc - wc)*eps_a + (1.-porsl_loc)*eps_r     ! [2](4)
         END IF

         ! add conductivity loss for imaginary part
         alpha = min(100.*wp, 26.)
         ecl = alpha*wc**2                          ! [2](6)
         eps = eps + jj*ecl                         ! [2](6)
      END IF

      ! mix dielectric constant of frozen and non-frozen soil
      IF (DEF_DA_CMEM) THEN
         IF (tc < -5.0d0) THEN
            ffrz = 1.0d0
         ELSE IF (tc < -0.5d0) THEN
            ffrz = 0.5d0
         ELSE
            ffrz = 0.0d0
         END IF
      ELSE
         IF (h2osoi(1) < 1e-10) THEN
            ffrz = 1
         ELSE
            ffrz = 1 - (liq_h2osoi(1))/h2osoi(1)
         END IF
      END IF
      eps = eps*(1.-ffrz) + eps_f*ffrz

   END SUBROUTINE diel_soil

!-----------------------------------------------------------------------

   SUBROUTINE diel_ice(t, eps_i)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate dielectric constant of pure ice
!
! REFERENCES:
!   [1] Matzler, C. (2006). Thermal Microwave Radiation: Applications for Remote Sensing p456-461
!-----------------------------------------------------------------------
      USE MOD_Const_Physical
      USE MOD_Precision
      IMPLICIT NONE

!------------------------ Dummy Argument ------------------------------
      real(r8), intent(in)  :: t             ! temperature (C)
      complex(r8), intent(out) :: eps_i      ! dielectric constant of ice water

!----------------------- Local Variables -------------------------------
      real(r8) :: betam                ! beta parameter by Mishima et al. (1983)
      real(r8) :: dbeta                ! corrected delta beta parameter
      real(r8) :: beta                 ! beta parameter
      real(r8) :: t_inv                ! modified inverse temperature
      real(r8) :: tk                   ! temperature (K)
      real(r8) :: alpha                ! alpha parameter
      real(r8) :: eps_i_r              ! real part of pure ice dielectric constant
      real(r8) :: eps_i_i              ! imaginary part of pure ice dielectric constant

!-----------------------------------------------------------------------

      ! C to K
      tk = t + tfrz

      ! eq.(5.33): calculate beta parameter by Mishima et al. (1983)
      betam = (0.0207/tk)*(exp(335./tk)/((exp(335./tk) - 1.)**2.)) + 1.16e-11*(fghz**2.)      !  [1](5.33)

      ! eq.(5.35): calculate delta beta parameter
      dbeta = exp(-10.02 + 0.0364*t)                              ! [1](5.35)

      ! eq.(5.34): calculate beta parameter
      beta = betam + dbeta                                        ! [1](5.34)

      ! eq.(5.32): calculate alpha parameter
      t_inv = 300./tk - 1                                         ! [1](p.457)
      alpha = (0.00504 + 0.0062*t_inv)*exp(-22.1*t_inv) !(GHz)    ! [1](5.32)

      ! eq.(5.30): calculate real part of pure ice dielectric constant
      eps_i_r = 3.1884 + 9.1e-4*t                                 ! [1](5.30)

      ! eq.(5.31): calculate imaginary part of pure ice dielectric constant
      eps_i_i = alpha/fghz + beta*fghz                            ! [1](5.31)

      ! calculate dielectric constant of pure ice
      eps_i = eps_i_r - jj*eps_i_i                                ! [1](5.31)

   END SUBROUTINE diel_ice

!-----------------------------------------------------------------------

   SUBROUTINE diel_water(type, swc, t, vf_sand, vf_clay, BD_all, sal, eps_w)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate dielectric constant of water in water (saline water)
!
! REFERENCES:
!   Ulaby FT, R. K. Moore, and A. K. Fung, Microwave Remote Sensing:
!   Active and Passive. Vol. III. From theory to applications. Artech House,
!   Norwood, MA., 1986
!
!   Stogryn, A. (1971): Equations for calculating the dielectric constant of
!   saline water, IEEE Transactions on Microwave Theory and Techniques,
!   Vol. MTT-19, 733-736.
!
!   Klein, L. A. and C. T. Swift (1977): An improved model
!   for the dielectric constant of sea water at microwave
!   frequencies, IEEE Transactions on  Antennas and Propagation,
!   Vol. AP-25, No. 1, 104-111.
!
!   Dobson, M.C., Ulaby, F.T., Hallikainen, M.T. and El-Rayes, M.A. (1985)
!   Microwave Dielectric Behavior of Wet Soil- Part II: Dielectric
!   Mixing Models. IEEE Transactions on Geoscience and Remote Sensing, GE-23, 35-46.
!   http://dx.doi.org/10.1109/TGRS.1985.289498
!-----------------------------------------------------------------------
      USE MOD_Const_Physical
      USE MOD_Precision
      IMPLICIT NONE

!------------------------ Dummy Argument ------------------------------
      integer, intent(in)      :: type       ! type of water, 0: pure water, 1: sea water, 2: soil water
      real(r8), intent(inout)  :: swc        ! soil water content (m3/m3)
      real(r8), intent(in)     :: t          ! soil temperature (C)
      real(r8), intent(in)     :: vf_sand    ! sand volume fraction
      real(r8), intent(in)     :: vf_clay    ! clay volume fraction
      real(r8), intent(in)     :: BD_all     ! bulk density(kg/m3)
      real(r8), intent(in)     :: sal        ! water salinity (psu)
      complex(r8), intent(out) :: eps_w      ! dielectric constant of soil water

!----------------------- Local Variables -------------------------------
      real(r8) :: sigma               ! ionic conductivity (S/m)
      real(r8) :: a, b                ! parameters
      real(r8) :: tau_w               ! relaxation time of pure water
      real(r8) :: eps_w0              ! static dielectric constant of pure water
      real(r8) :: bd_all_loc          ! bulk density (g/cm3)

!-----------------------------------------------------------------------

      ! [3] eq.16: tau(T, sal) = tau_w(T) * b(sal, T)
      ! calculate relaxation time of pure water (Stogryn)
      tau_w = 1.768e-11 - 6.068e-13*t + 1.104e-14*t**2 - 8.111e-17*t**3               ! [3](17)
      b = 1.000 + 2.282e-5*sal*t - 7.638e-4*sal - 7.760e-6*sal**2 + 1.105e-8*sal**3   ! [3](18)
      tau_w = tau_w*b                                                                 ! [3](16)

      ! [3] eq.13: eps_w0(sal, T) = eps_w0(T) * a(sal, T)
      ! static dielectric constant of pure water (Klein and Swift)
      eps_w0 = 87.134 - 1.949e-1*t - 1.276e-2*t**2 + 2.491e-4*t**3                    ! [3](14)
      a = 1.000 + 1.613e-5*sal*t - 3.656e-3*sal + 3.210e-5*sal**2 - 4.232e-7*sal**3   ! [3](15)
      eps_w0 = eps_w0*a                                                               ! [3](13)

      IF (type == 0) THEN  ! pure water
         ! [1] eq.19
         eps_w0 = 88.045 - 0.4147*t + 6.295e-4*t**2 + 1.075e-5*t**3
         eps_w = eps_w_inf + (eps_w0 - eps_w_inf)/(1 - jj*omega*tau_w)

      ELSEIF (type == 1) THEN  ! sea water
         ! calculate ionic conductivity [1] eq.27, eq.28
         sigma = sal*(0.182521 - 1.46192e-3*sal + 2.09324e-5*sal**2 - 1.28205e-7*sal**3) &
            *exp(-1.*(25 - t)* &
            (2.033e-2 + 1.266e-4*(25 - t) + 2.464e-6*(25 - t)**2 &
            - sal*(1.849e-5 - 2.551e-7*(25 - t) + 2.551e-8*(25 - t)**2)))

         ! diel constant of sea water [1] eq.21
         eps_w = eps_w_inf + (eps_w0 - eps_w_inf)/(1 - jj*omega*tau_w) + jj*sigma/(omega*eps_0)
      ELSE
         ! calculate soil conductivity
         bd_all_loc = (vf_sand*1.60d0 + vf_clay*1.10d0 + (100.0d0 - vf_sand - vf_clay)*1.20d0)/100.0d0
         sigma = -1.645 + 1.939*bd_all_loc - 0.02256*vf_sand + 0.01594*vf_clay
         IF (sigma < 0.) THEN
            sigma = 0. ! negative for very sandy soils with low bulk density
         END IF

         ! calculate dielectric constant of soil-water by modified Debye expression
         swc = max(0.001, swc)
         eps_w = eps_w_inf + (eps_w0 - eps_w_inf)/(1 - jj*omega*tau_w) &
            + jj*sigma/(omega*eps_0)*(rho_soil - bd_all_loc)/(rho_soil*swc)
      END IF

   END SUBROUTINE diel_water

!-----------------------------------------------------------------------

   SUBROUTINE smooth_reflectivity(eps_surf, r_s)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the smooth surface reflectivity by Fresnel Law
!
! REFERENCES:
!   [1] Njoku and Kong, 1977: Theory for passive microwave remote sensing
!       of near-surface soil moisture. Journal of Geophysical Research,
!       Vol. 82, No. 20, 3108-3118.
!-----------------------------------------------------------------------
      USE MOD_Precision
      IMPLICIT NONE

!------------------------ Dummy Argument ------------------------------
      complex(r8), intent(in) :: eps_surf       ! dielectric constant of the surface
      real(r8), intent(out)   :: r_s(2)         ! reflectivities of flat surfaces for H and V polarizations

!----------------------- Local Variables -------------------------------
      complex(r8) :: g                          ! parameter in Fresnel Law

!-----------------------------------------------------------------------

      g = sqrt(eps_surf - sin(theta)**2)                                     ! [1]p4
      r_s(1) = abs((cos(theta) - g)/(cos(theta) + g))**2.                    ! [1]p4
      r_s(2) = abs((cos(theta)*eps_surf - g)/(cos(theta)*eps_surf + g))**2.  ! [1]p4

   END SUBROUTINE smooth_reflectivity

!-----------------------------------------------------------------------

   SUBROUTINE rough_reflectivity(is_desert, patchclass, r_s, r_r)

!-----------------------------------------------------------------------
! DESCRIPTION:
!   Calculate the rough surface reflectivity
!
! REFERENCES:
!   [1] Kerr and Njokui, 1990: A Semiempirical Model For Interpreting Microwave
!       Emission From Semiarid Land Surfaces as Seen From Space
!       IEEE Trans. Geosci. Rem. Sens., Vol.28, No.3, 384-393.
!
!   [2] Choudhury et al., 1979: Effect of surface roughness on the
!       microwave emission from soils, J.Geo.Res. Vol.84, 5699-5706
!
!   [3] Wang and Choudhury, 1981: Remote sensing of soil moisture
!       content over bare field at 1.4 GHz frequency, J.Geo.Res.
!       Vol.86, 5277-5287
!-----------------------------------------------------------------------
      USE MOD_Precision
      USE MOD_DA_Const
      IMPLICIT NONE

!------------------------ Dummy Argument ------------------------------
      logical, intent(in)   :: is_desert   ! flag for desert soil
      integer, intent(in)   :: patchclass  ! patch class
      real(r8), intent(in)  :: r_s(2)      ! reflectivities of flat surfaces for H and V polarizations
      real(r8), intent(out) :: r_r(2)      ! reflectivities of rough surfaces for H and V polarizations

!----------------------- Local Variables -------------------------------
      real(r8) :: Q                        ! parameter for polarization mixing

!-----------------------------------------------------------------------

      IF (is_desert) THEN
         r_r = r_s
      ELSE
         ! calculate parameter for polarization mixing due to surface roughness
         IF (fghz < 2.) THEN
            Q = 0. ! Q is assumed zero at low frequency
         ELSE
            Q = 0.35*(1.0 - exp(-0.6*rgh_surf**2*fghz))     !    [1](16)
         END IF

         ! calculate rough surface reflectivity (using CMEM parameters `hr`)
         r_r(1) = (Q*r_s(2) + (1.-Q)*r_s(1))*exp(-hr(patchclass)*cos(theta)**nrh(patchclass))     !  [1](6)
         r_r(2) = (Q*r_s(1) + (1.-Q)*r_s(2))*exp(-hr(patchclass)*cos(theta)**nrv(patchclass))     !  [1](6)
      END IF

   END SUBROUTINE rough_reflectivity

!-----------------------------------------------------------------------
END MODULE MOD_DA_ObsOperator
#endif
