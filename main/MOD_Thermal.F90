#include <define.h>

MODULE MOD_Thermal

!-----------------------------------------------------------------------
  USE MOD_Precision
  IMPLICIT NONE
  SAVE

! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: THERMAL


!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------


  SUBROUTINE THERMAL (ipatch      ,patchtype   ,lb          ,deltim     ,&
                      trsmx0      ,zlnd        ,zsno        ,csoilc     ,&
                      dewmx       ,capr        ,cnfac       ,vf_quartz  ,&
                      vf_gravels  ,vf_om       ,vf_sand     ,wf_gravels ,&
                      wf_sand     ,csol        ,porsl       ,psi0       ,&
#ifdef Campbell_SOIL_MODEL
                      bsw         ,                                      &
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
                      theta_r     ,alpha_vgm   ,n_vgm       ,L_vgm      ,&
                      sc_vgm      ,fc_vgm      ,                         &
#endif
                      k_solids    ,dksatu      ,dksatf      ,dkdry      ,&
                      BA_alpha    ,BA_beta                              ,&
                      lai         ,laisun      ,laisha                  ,&
                      sai         ,htop        ,hbot        ,sqrtdi     ,&
                      rootfr      ,rstfacsun_out   ,rstfacsha_out       ,&
                      gssun_out   ,gssha_out   ,&
                      assimsun_out,etrsun_out  ,assimsha_out,etrsha_out ,&
! photosynthesis and plant hydraulic variables
                      effcon      ,vmax25      ,hksati      ,smp     ,hk,&
                      kmax_sun    ,kmax_sha    ,kmax_xyl    ,kmax_root  ,&
                      psi50_sun   ,psi50_sha   ,psi50_xyl   ,psi50_root ,&
                      ck          ,vegwp       ,gs0sun      ,gs0sha     ,&
!Ozone stress variables
                      lai_old     ,o3uptakesun ,o3uptakesha ,forc_ozone, &
!end ozone stress variables
                      slti        ,hlti        ,shti        ,hhti       ,&
                      trda        ,trdm        ,trop        ,gradm      ,&
                      binter      ,extkn       ,forc_hgt_u  ,forc_hgt_t ,&
                      forc_hgt_q  ,forc_us     ,forc_vs     ,forc_t     ,&
                      forc_q      ,forc_rhoair ,forc_psrf   ,forc_pco2m ,&
                      forc_hpbl                                         ,&
                      forc_po2m   ,coszen      ,parsun      ,parsha     ,&
                      sabvsun     ,sabvsha     ,sabg        ,frl        ,&
                      extkb       ,extkd       ,thermk      ,fsno       ,&
                      sigf        ,dz_soisno   ,z_soisno    ,zi_soisno  ,&
                      tleaf       ,t_soisno    ,wice_soisno ,wliq_soisno,&
                      ldew,ldew_rain,ldew_snow ,scv,snowdp  ,imelt      ,&
                      taux        ,tauy        ,fsena       ,fevpa      ,&
                      lfevpa      ,fsenl       ,fevpl       ,etr        ,&
                      fseng       ,fevpg       ,olrg        ,fgrnd      ,&
                      rootr       ,qseva       ,qsdew       ,qsubl      ,&
                      qfros       ,sm          ,tref        ,qref       ,&
                      trad        ,rst         ,assim       ,respc      ,&
                      errore      ,emis        ,z0m         ,zol        ,&
                      rib         ,ustar       ,qstar       ,tstar      ,&
                      fm          ,fh          ,fq          ,pg_rain    ,&
                      pg_snow     ,t_precip    ,qintr_rain  ,qintr_snow ,&
                      snofrz      ,sabg_lyr                              )

!=======================================================================
! this is the main subroutine to execute the calculation
! of thermal processes and surface fluxes
!
! Original author : Yongjiu Dai, 09/15/1999; 08/30/2002
!
! FLOW DIAGRAM FOR THERMAL.F90
!
! THERMAL ===> qsadv
!              groundfluxes
!              eroot                             |dewfraction
!              LeafTemperature   |               |qsadv
!              LeafTemperaturePC |  ---------->  |moninobukini
!                                                |moninobuk
!                                                |MOD_AssimStomataConductance
!
!              groundTem            ---------->   meltf
!
!
! REVISIONS:
! Hua Yuan, 08/2019: added initial codes for PFT and Plant Community (PC)
!                    vegetation classification processes
! Nan Wei,  01/2021: added variables passing of plant hydraulics and precipitation sensible heat
!                    with canopy and ground for PFT and Plant Community (PC)
!=======================================================================

  USE MOD_Precision
  USE MOD_Vars_Global
  USE MOD_Const_PFT
  USE MOD_Const_Physical, only: denh2o,roverg,hvap,hsub,rgas,cpair,&
                                stefnc,denice,tfrz,vonkar,grav,cpliq,cpice
  USE MOD_FrictionVelocity
  USE MOD_Eroot
  USE MOD_GroundFluxes
  USE MOD_LeafTemperature
  USE MOD_LeafTemperaturePC
  USE MOD_GroundTemperature
  USE MOD_Qsadv
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
  USE MOD_LandPFT, only: patch_pft_s, patch_pft_e
  USE MOD_Vars_TimeInvariants, only: patchclass
  USE MOD_Vars_PFTimeInvariants
  USE MOD_Vars_PFTimeVariables
  USE MOD_Vars_1DPFTFluxes
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
  USE MOD_Hydro_SoilFunction, only : soil_psi_from_vliq
#endif
  USE MOD_SPMD_Task
  USE MOD_Namelist, only: DEF_USE_PLANTHYDRAULICS,DEF_USE_PFT,DEF_USE_PC

  IMPLICIT NONE

!---------------------Argument------------------------------------------

  integer, intent(in) :: &
        ipatch,      &! patch index
        lb,          &! lower bound of array
        patchtype     ! land water type (0=soil, 1=urban or built-up, 2=wetland,
                      !                  3=glacier/ice sheet, 4=land water bodies)

  real(r8), intent(in) :: &
        deltim,      &! model time step [second]
        trsmx0,      &! max transpiration for moist soil+100% veg.  [mm/s]
        zlnd,        &! roughness length for soil [m]
        zsno,        &! roughness length for snow [m]
        csoilc,      &! drag coefficient for soil under canopy [-]
        dewmx,       &! maximum dew
        capr,        &! tuning factor to turn first layer T into surface T
        cnfac,       &! Crank Nicholson factor between 0 and 1

        ! soil physical parameters
        vf_quartz (1:nl_soil), &! volumetric fraction of quartz within mineral soil
        vf_gravels(1:nl_soil), &! volumetric fraction of gravels
        vf_om     (1:nl_soil), &! volumetric fraction of organic matter
        vf_sand   (1:nl_soil), &! volumetric fraction of sand
        wf_gravels(1:nl_soil), &! gravimetric fraction of gravels
        wf_sand   (1:nl_soil), &! gravimetric fraction of sand
        csol      (1:nl_soil), &! heat capacity of soil solids [J/(m3 K)]
        porsl     (1:nl_soil), &! soil porosity [-]
        psi0      (1:nl_soil), &! soil water suction, negative potential [mm]
#ifdef Campbell_SOIL_MODEL
        bsw(1:nl_soil),    &! clapp and hornbereger "b" parameter [-]
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
        theta_r  (1:nl_soil),  &
        alpha_vgm(1:nl_soil),  &
        n_vgm    (1:nl_soil),  &
        L_vgm    (1:nl_soil),  &
        sc_vgm   (1:nl_soil),  &
        fc_vgm   (1:nl_soil),  &
#endif
        k_solids  (1:nl_soil), &! thermal conductivity of minerals soil [W/m-K]
        dkdry     (1:nl_soil), &! thermal conductivity of dry soil [W/m-K]
        dksatu    (1:nl_soil), &! thermal conductivity of saturated unfrozen soil [W/m-K]
        dksatf    (1:nl_soil), &! thermal conductivity of saturated frozen soil [W/m-K]
        hksati    (1:nl_soil), &! hydraulic conductivity at saturation [mm h2o/s]
        BA_alpha  (1:nl_soil), &! alpha in Balland and Arp(2005) thermal conductivity scheme
        BA_beta   (1:nl_soil), &! beta in Balland and Arp(2005) thermal conductivity scheme

        ! vegetation parameters
        lai,         &! adjusted leaf area index for seasonal variation [-]
        sai,         &! stem area index  [-]
        htop,        &! canopy crown top height [m]
        hbot,        &! canopy crown bottom height [m]
        sqrtdi,      &! inverse sqrt of leaf dimension [m**-0.5]
        rootfr(1:nl_soil),&! root fraction

        effcon,      &! quantum efficiency of RuBP regeneration (mol CO2/mol quanta)
        vmax25,      &! maximum carboxylation rate at 25 C at canopy top
        kmax_sun,    &
        kmax_sha,    &
        kmax_xyl,    &
        kmax_root,   &
        psi50_sun,   &! water potential at 50% loss of sunlit leaf tissue conductance (mmH2O)
        psi50_sha,   &! water potential at 50% loss of shaded leaf tissue conductance (mmH2O)
        psi50_xyl,   &! water potential at 50% loss of xylem tissue conductance (mmH2O)
        psi50_root,  &! water potential at 50% loss of root tissue conductance (mmH2O)
        ck,          &! shape-fitting parameter for vulnerability curve (-)
        slti,        &! slope of low temperature inhibition function      [s3]
        hlti,        &! 1/2 point of low temperature inhibition function  [s4]
        shti,        &! slope of high temperature inhibition function     [s1]
        hhti,        &! 1/2 point of high temperature inhibition function [s2]
        trda,        &! temperature coefficient in gs-a model             [s5]
        trdm,        &! temperature coefficient in gs-a model             [s6]
        trop,        &! temperature coefficient in gs-a model
        gradm,       &! conductance-photosynthesis slope parameter
        binter,      &! conductance-photosynthesis intercept
        extkn,       &! coefficient of leaf nitrogen allocation

        ! atmospherical variables and observational height
        forc_hgt_u,  &! observational height of wind [m]
        forc_hgt_t,  &! observational height of temperature [m]
        forc_hgt_q,  &! observational height of humidity [m]
        forc_us,     &! wind component in eastward direction [m/s]
        forc_vs,     &! wind component in northward direction [m/s]
        forc_t,      &! temperature at agcm reference height [kelvin]
        forc_q,      &! specific humidity at agcm reference height [kg/kg]
        forc_rhoair, &! density air [kg/m3]
        forc_psrf,   &! atmosphere pressure at the surface [pa]
        forc_pco2m,  &! CO2 concentration in atmos. (pascals)
        forc_po2m,   &! O2 concentration in atmos. (pascals)
        forc_hpbl,   &! atmospheric boundary layer height [m]
        pg_rain,     &! rainfall onto ground including canopy runoff [kg/(m2 s)]
        pg_snow,     &! snowfall onto ground including canopy runoff [kg/(m2 s)]
        t_precip,    &! snowfall/rainfall temperature [kelvin]
        qintr_rain,  &! rainfall interception (mm h2o/s)
        qintr_snow,  &! snowfall interception (mm h2o/s)

        ! radiative fluxes
        coszen,      &! cosine of the solar zenith angle
        parsun,      &! photosynthetic active radiation by sunlit leaves (W m-2)
        parsha,      &! photosynthetic active radiation by shaded leaves (W m-2)
        sabvsun,     &! solar radiation absorbed by vegetation [W/m2]
        sabvsha,     &! solar radiation absorbed by vegetation [W/m2]
        sabg,        &! solar radiation absorbed by ground [W/m2]
        frl,         &! atmospheric infrared (longwave) radiation [W/m2]
        extkb,       &! (k, g(mu)/mu) direct solar extinction coefficient
        extkd,       &! diffuse and scattered diffuse PAR extinction coefficient
        thermk,      &! canopy gap fraction for tir radiation

        ! state variable (1)
        fsno,        &! fraction of ground covered by snow
        sigf,        &! fraction of veg cover, excluding snow-covered veg [-]
        dz_soisno(lb:nl_soil),  &! layer thickiness [m]
        z_soisno (lb:nl_soil),  &! node depth [m]
        zi_soisno(lb-1:nl_soil)  ! interface depth [m]

  real(r8), intent(in) :: &
        sabg_lyr(lb:1)           ! snow layer aborption

        ! state variables (2)
  real(r8), intent(inout) :: &
        vegwp(1:nvegwcs),&! vegetation water potential
        gs0sun,      &!
        gs0sha,      &!
!Ozone stress variables
        lai_old    , & ! lai in last time step
        o3uptakesun, & ! Ozone does, sunlit leaf (mmol O3/m^2)
        o3uptakesha, & ! Ozone does, shaded leaf (mmol O3/m^2)
        forc_ozone , & ! Ozone
!end ozone stress variables
        tleaf,       &! shaded leaf temperature [K]
        t_soisno(lb:nl_soil),   &! soil temperature [K]
        wice_soisno(lb:nl_soil),&! ice lens [kg/m2]
        wliq_soisno(lb:nl_soil),&! liqui water [kg/m2]
        smp(1:nl_soil)         ,&! soil matrix potential [mm]
        hk(1:nl_soil)          ,&! hydraulic conductivity [mm h2o/s]

        ldew,        &! depth of water on foliage [kg/(m2 s)]
        ldew_rain,   &! depth of rain on foliage [kg/(m2 s)]
        ldew_snow,   &! depth of rain on foliage [kg/(m2 s)]
        scv,         &! snow cover, water equivalent [mm, kg/m2]
        snowdp        ! snow depth [m]

  real(r8), intent(out) :: &
        snofrz (lb:0) !snow freezing rate (col,lyr) [kg m-2 s-1]

  integer, intent(out) :: &
       imelt(lb:nl_soil) ! flag for melting or freezing [-]

  real(r8), intent(out) :: &
       laisun,       &! sunlit leaf area index
       laisha,       &! shaded leaf area index
       gssun_out,    &! sunlit stomata conductance
       gssha_out,    &! shaded stomata conductance
       rstfacsun_out,&! factor of soil water stress on sunlit leaf
       rstfacsha_out  ! factor of soil water stress on shaded leaf
  real(r8), intent(out) :: &
       assimsun_out ,&
       etrsun_out   ,&
       assimsha_out ,&
       etrsha_out

        ! Output fluxes
  real(r8), intent(out) :: &
        taux,        &! wind stress: E-W [kg/m/s**2]
        tauy,        &! wind stress: N-S [kg/m/s**2]
        fsena,       &! sensible heat from canopy height to atmosphere [W/m2]
        fevpa,       &! evapotranspiration from canopy height to atmosphere [mm/s]
        lfevpa,      &! latent heat flux from canopy height to atmosphere [W/m2]
        fsenl,       &! ensible heat from leaves [W/m2]
        fevpl,       &! evaporation+transpiration from leaves [mm/s]
        etr,         &! transpiration rate [mm/s]
        fseng,       &! sensible heat flux from ground [W/m2]
        fevpg,       &! evaporation heat flux from ground [mm/s]
        olrg,        &! outgoing long-wave radiation from ground+canopy
        fgrnd,       &! ground heat flux [W/m2]
        rootr(1:nl_soil),&! root resistance of a layer, all layers add to 1

        qseva,       &! ground surface evaporation rate (mm h2o/s)
        qsdew,       &! ground surface dew formation (mm h2o /s) [+]
        qsubl,       &! sublimation rate from snow pack (mm h2o /s) [+]
        qfros,       &! surface dew added to snow pack (mm h2o /s) [+]

        sm,          &! rate of snowmelt [kg/(m2 s)]
        tref,        &! 2 m height air temperature [kelvin]
        qref,        &! 2 m height air specific humidity
        trad,        &! radiative temperature [K]

        rst,         &! stomatal resistance (s m-1)
        assim,       &! assimilation
        respc,       &! respiration

        ! additional variables required by coupling with WRF or RSM model
        emis,        &! averaged bulk surface emissivity
        z0m,         &! effective roughness [m]
        zol,         &! dimensionless height (z/L) used in Monin-Obukhov theory
        rib,         &! bulk Richardson number in surface layer
        ustar,       &! u* in similarity theory [m/s]
        qstar,       &! q* in similarity theory [kg/kg]
        tstar,       &! t* in similarity theory [K]
        fm,          &! integral of profile function for momentum
        fh,          &! integral of profile function for heat
        fq            ! integral of profile function for moisture

!---------------------Local Variables-----------------------------------

  integer i,j

  real(r8) :: &
       cgrnd,        &! deriv. of soil energy flux wrt to soil temp [w/m2/k]
       cgrndl,       &! deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
       cgrnds,       &! deriv of soil latent heat flux wrt soil temp [w/m**2/k]
       degdT,        &! d(eg)/dT
       dqgdT,        &! d(qg)/dT
       dlrad,        &! downward longwave radiation blow the canopy [W/m2]
       eg,           &! water vapor pressure at temperature T [pa]
       egsmax,       &! max. evaporation which soil can provide at one time step
       egidif,       &! the excess of evaporation over "egsmax"
       emg,          &! ground emissivity (0.97 for snow,
                      ! glaciers and water surface; 0.96 for soil and wetland)
       errore,       &! energy balnce error [w/m2]
       etrc,         &! maximum possible transpiration rate [mm/s]
       fac,          &! soil wetness of surface layer
       fact(lb:nl_soil), &! used in computing tridiagonal matrix
       fsun,         &! fraction of sunlit canopy
       hr,           &! relative humidity
       htvp,         &! latent heat of vapor of water (or sublimation) [j/kg]
       olru,         &! olrg excluding dwonwelling reflection [W/m2]
       olrb,         &! olrg assuming blackbody emission [W/m2]
       psit,         &! negative potential of soil
       qg,           &! ground specific humidity [kg/kg]
       qsatg,        &! saturated humidity [kg/kg]
       qsatgdT,      &! d(qsatg)/dT
       qred,         &! soil surface relative humidity
       sabv,         &! solar absorbed by canopy [W/m2]
       thm,          &! intermediate variable (forc_t+0.0098*forc_hgt_t)
       th,           &! potential temperature (kelvin)
       thv,          &! virtual potential temperature (kelvin)
       rstfac,       &! factor of soil water stress
       t_grnd,       &! ground surface temperature [K]
       t_grnd_bef,   &! ground surface temperature [K]
       t_soisno_bef(lb:nl_soil), &! soil/snow temperature before update
       tinc,         &! temperature difference of two time step
       ur,           &! wind speed at reference height [m/s]
       ulrad,        &! upward longwave radiation above the canopy [W/m2]
       wice0(lb:nl_soil),&! ice mass from previous time-step
       wliq0(lb:nl_soil),&! liquid mass from previous time-step
       wx,           &! patitial volume of ice and water of surface layer
       xmf,          &! total latent heat of phase change of ground water
       hprl           ! precipitation sensible heat from canopy

  real(r8) :: z0m_g,z0h_g,zol_g,obu_g,rib_g,ustar_g,qstar_g,tstar_g
  real(r8) :: fm10m,fm_g,fh_g,fq_g,fh2m,fq2m,um,obu
!Ozone stress variables
  real(r8) :: o3coefv_sun, o3coefv_sha, o3coefg_sun, o3coefg_sha
!end ozone stress variables

  integer p, ps, pe, pc

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
  real(r8), allocatable :: rootr_p     (:,:)
  real(r8), allocatable :: etrc_p        (:)
  real(r8), allocatable :: rstfac_p      (:)
  real(r8), allocatable :: rstfacsun_p   (:)
  real(r8), allocatable :: rstfacsha_p   (:)
  real(r8), allocatable :: gssun_p       (:)
  real(r8), allocatable :: gssha_p       (:)
  real(r8), allocatable :: fsun_p        (:)
  real(r8), allocatable :: sabv_p        (:)
  real(r8), allocatable :: cgrnd_p       (:)
  real(r8), allocatable :: cgrnds_p      (:)
  real(r8), allocatable :: cgrndl_p      (:)
  real(r8), allocatable :: dlrad_p       (:)
  real(r8), allocatable :: ulrad_p       (:)
  real(r8), allocatable :: zol_p         (:)
  real(r8), allocatable :: rib_p         (:)
  real(r8), allocatable :: ustar_p       (:)
  real(r8), allocatable :: qstar_p       (:)
  real(r8), allocatable :: tstar_p       (:)
  real(r8), allocatable :: fm_p          (:)
  real(r8), allocatable :: fh_p          (:)
  real(r8), allocatable :: fq_p          (:)
  real(r8), allocatable :: hprl_p        (:)
  real(r8), allocatable :: assimsun_p    (:)
  real(r8), allocatable :: etrsun_p      (:)
  real(r8), allocatable :: assimsha_p    (:)
  real(r8), allocatable :: etrsha_p      (:)
#endif

!=======================================================================
! [1] Initial set and propositional variables
!=======================================================================
      ! emissivity
      emg = 0.96
      IF (scv>0. .or. patchtype==3) emg = 0.97

      ! fluxes
      taux   = 0.;  tauy   = 0.
      fsena  = 0.;  fevpa  = 0.
      lfevpa = 0.;  fsenl  = 0.
      fevpl  = 0.;  etr    = 0.
      fseng  = 0.;  fevpg  = 0.
      dlrad  = frl
      ulrad  = frl*(1.-emg) + emg*stefnc*t_soisno(lb)**4
      cgrnds = 0.;  cgrndl = 0.
      cgrnd  = 0.;  tref   = 0.
      qref   = 0.;  rst    = 2.0e4
      assim  = 0.;  respc  = 0.
      hprl   = 0.

      emis   = 0.;  z0m    = 0.
      zol    = 0.;  rib    = 0.
      ustar  = 0.;  qstar  = 0.
      tstar  = 0.;  rootr  = 0.

      ! temperature and water mass from previous time step
      t_grnd = t_soisno(lb)

      t_soisno_bef(lb:) = t_soisno(lb:)
      t_grnd_bef = t_grnd
      wice0(lb:) = wice_soisno(lb:)
      wliq0(lb:) = wliq_soisno(lb:)

      ! latent heat, assumed that the sublimation occured only as wliq_soisno=0
      htvp = hvap
      IF (wliq_soisno(lb)<=0. .and. wice_soisno(lb)>0.) htvp = hsub

      ! potential temperatur at the reference height
      thm = forc_t + 0.0098*forc_hgt_t              !intermediate variable equivalent to
                                                    !forc_t*(pgcm/forc_psrf)**(rgas/cpair)
      th = forc_t*(100000./forc_psrf)**(rgas/cpair) !potential T
      thv = th*(1.+0.61*forc_q)                     !virtual potential T
      ur = max(0.1,sqrt(forc_us*forc_us+forc_vs*forc_vs)) !limit set to 0.1

!=======================================================================
! [2] specific humidity and its derivative at ground surface
!=======================================================================

      qred = 1.
      CALL qsadv(t_grnd,forc_psrf,eg,degdT,qsatg,qsatgdT)

      IF (patchtype<=1) THEN            !soil ground
         wx   = (wliq_soisno(1)/denh2o + wice_soisno(1)/denice)/dz_soisno(1)
         IF (porsl(1) < 1.e-6) THEN     !bed rock
            fac  = 0.001
         ELSE
            fac  = min(1.,wx/porsl(1))
            fac  = max( fac, 0.001 )
         ENDIF

#ifdef Campbell_SOIL_MODEL
         psit = psi0(1) * fac ** (- bsw(1) )   !psit = max(smpmin, psit)
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
         psit = soil_psi_from_vliq ( fac*(porsl(1)-theta_r(1)) + theta_r(1), &
            porsl(1), theta_r(1), psi0(1), &
            5, (/alpha_vgm(1), n_vgm(1), L_vgm(1), sc_vgm(1), fc_vgm(1)/))
#endif
         psit = max( -1.e8, psit )
         hr   = exp(psit/roverg/t_grnd)
         qred = (1.-fsno)*hr + fsno
      ENDIF

      qg = qred*qsatg
      dqgdT = qred*qsatgdT

      IF (qsatg > forc_q .and. forc_q > qred*qsatg) THEN
        qg = forc_q; dqgdT = 0.
      ENDIF

!=======================================================================
! [3] Compute sensible and latent fluxes and their derivatives with respect
!     to ground temperature using ground temperatures from previous time step.
!=======================================================================

IF (patchtype == 0) THEN

!=======================================================================
!=======================================================================
#if(defined LULC_USGS || defined LULC_IGBP)
      CALL groundfluxes (zlnd,zsno,forc_hgt_u,forc_hgt_t,forc_hgt_q,forc_hpbl, &
                         forc_us,forc_vs,forc_t,forc_q,forc_rhoair,forc_psrf, &
                         ur,thm,th,thv,t_grnd,qg,dqgdT,htvp, &
                         fsno,cgrnd,cgrndl,cgrnds, &
                         taux,tauy,fseng,fevpg,tref,qref, &
                         z0m_g,z0h_g,zol_g,rib_g,ustar_g,qstar_g,tstar_g,fm_g,fh_g,fq_g)

      ! SAVE variables for bareground case
      obu_g = forc_hgt_u / zol_g

!=======================================================================
! [4] Canopy temperature, fluxes from the canopy
!=======================================================================

      sabv = sabvsun + sabvsha

      IF (lai+sai > 1e-6) THEN

         ! soil water stress factor on stomatal resistance
         CALL eroot (nl_soil,trsmx0,porsl,&
#ifdef Campbell_SOIL_MODEL
            bsw,&
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
            theta_r, alpha_vgm, n_vgm, L_vgm, sc_vgm, fc_vgm, &
#endif
            psi0,rootfr,dz_soisno,t_soisno,wliq_soisno,rootr,etrc,rstfac)

         ! fraction of sunlit and shaded leaves of canopy
         fsun = ( 1. - exp(-min(extkb*lai,40.))) / max( min(extkb*lai,40.), 1.e-6 )

         IF (coszen<=0.0 .or. sabv<1.) fsun = 0.5

         laisun = lai*fsun
         laisha = lai*(1-fsun)
         rstfacsun_out = rstfac
         rstfacsha_out = rstfac

         CALL LeafTemperature(ipatch,1,deltim,csoilc,dewmx     ,htvp       ,&
                 lai        ,sai        ,htop      ,hbot       ,sqrtdi     ,&
                 effcon     ,vmax25     ,slti      ,hlti       ,shti       ,&
                 hhti       ,trda       ,trdm      ,trop       ,gradm      ,&
                 binter     ,extkn      ,extkb     ,extkd      ,forc_hgt_u ,&
                 forc_hgt_t ,forc_hgt_q ,forc_us   ,forc_vs    ,thm        ,&
                 th         ,thv        ,forc_q    ,forc_psrf  ,forc_rhoair,&
                 parsun     ,parsha     ,sabv      ,frl        ,fsun       ,&
                 thermk     ,rstfacsun_out         ,rstfacsha_out          ,&
                 gssun_out  ,gssha_out  ,forc_po2m ,forc_pco2m ,z0h_g      ,&
                 obu_g      ,ustar_g    ,zlnd      ,zsno       ,fsno       ,&
                 sigf       ,etrc       ,t_grnd    ,qg         ,dqgdT      ,&
                 emg        ,tleaf      ,ldew      ,ldew_rain  ,ldew_snow  ,&
                 taux       ,tauy       ,&
                 fseng      ,fevpg      ,cgrnd     ,cgrndl     ,cgrnds     ,&
                 tref       ,qref       ,rst       ,assim      ,respc      ,&
                 fsenl      ,fevpl      ,etr       ,dlrad      ,ulrad      ,&
                 z0m        ,zol        ,rib       ,ustar      ,qstar      ,&
                 tstar      ,fm         ,fh        ,fq         ,rootfr     ,&
                 kmax_sun    ,kmax_sha  ,kmax_xyl  ,kmax_root  ,psi50_sun  ,&
                 psi50_sha   ,psi50_xyl ,psi50_root,ck         ,vegwp      ,&
                 gs0sun      ,gs0sha                                       ,&
                 assimsun_out,etrsun_out,assimsha_out          ,etrsha_out ,&
!Ozone stress variables
                 o3coefv_sun ,o3coefv_sha ,o3coefg_sun ,o3coefg_sha, &
                 lai_old     ,o3uptakesun ,o3uptakesha ,forc_ozone , &
!end ozone stress variables
                 forc_hpbl                                                 ,&
                 qintr_rain  ,qintr_snow,t_precip  ,hprl       ,smp        ,&
                 hk(1:)      ,hksati(1:),rootr(1:)                         )
      ENDIF

    ! equate canopy temperature to air over bareland.
    ! required as sigf=0 carried over to next time step
      IF (lai+sai <= 1e-6) THEN
         tleaf         = forc_t
         laisun        = 0.
         laisha        = 0.
         ldew_rain     = 0.
         ldew_snow     = 0.
         ldew          = 0.
         rstfacsun_out = 0.
         rstfacsha_out = 0.
         if(DEF_USE_PLANTHYDRAULICS)THEN
            vegwp = -2.5e4
         ENDIF
      ENDIF
#endif


!=======================================================================
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)

      ps = patch_pft_s(ipatch)
      pe = patch_pft_e(ipatch)

      allocate ( rootr_p (nl_soil, ps:pe) )
      allocate ( etrc_p           (ps:pe) )
      allocate ( rstfac_p         (ps:pe) )
      allocate ( rstfacsun_p      (ps:pe) )
      allocate ( rstfacsha_p      (ps:pe) )
      allocate ( gssun_p          (ps:pe) )
      allocate ( gssha_p          (ps:pe) )
      allocate ( fsun_p           (ps:pe) )
      allocate ( sabv_p           (ps:pe) )
IF (DEF_USE_PFT .or. patchclass(ipatch)==CROPLAND) THEN
      allocate ( cgrnd_p          (ps:pe) )
      allocate ( cgrnds_p         (ps:pe) )
      allocate ( cgrndl_p         (ps:pe) )
      allocate ( dlrad_p          (ps:pe) )
      allocate ( ulrad_p          (ps:pe) )
      allocate ( zol_p            (ps:pe) )
      allocate ( rib_p            (ps:pe) )
      allocate ( ustar_p          (ps:pe) )
      allocate ( qstar_p          (ps:pe) )
      allocate ( tstar_p          (ps:pe) )
      allocate ( fm_p             (ps:pe) )
      allocate ( fh_p             (ps:pe) )
      allocate ( fq_p             (ps:pe) )
ENDIF
      allocate ( hprl_p           (ps:pe) )
      allocate ( assimsun_p       (ps:pe) )
      allocate ( etrsun_p         (ps:pe) )
      allocate ( assimsha_p       (ps:pe) )
      allocate ( etrsha_p         (ps:pe) )

      ! always DO CALL groundfluxes
      CALL groundfluxes (zlnd,zsno,forc_hgt_u,forc_hgt_t,forc_hgt_q,forc_hpbl, &
                         forc_us,forc_vs,forc_t,forc_q,forc_rhoair,forc_psrf, &
                         ur,thm,th,thv,t_grnd,qg,dqgdT,htvp, &
                         fsno,cgrnd,cgrndl,cgrnds, &
                         taux,tauy,fseng,fevpg,tref,qref, &
                         z0m_g,z0h_g,zol_g,rib_g,ustar_g,qstar_g,tstar_g,fm_g,fh_g,fq_g)

      obu_g = forc_hgt_u / zol_g

!=======================================================================
! [4] Canopy temperature, fluxes from the canopy
!=======================================================================

      sabv_p(ps:pe) = sabvsun_p(ps:pe) + sabvsha_p(ps:pe)
      sabv = sabvsun + sabvsha

      DO i = ps, pe
         p = pftclass(i)

         IF (lai_p(i)+sai_p(i) > 1e-6) THEN

            CALL eroot (nl_soil,trsmx0,porsl,&
#ifdef Campbell_SOIL_MODEL
               bsw, &
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
               theta_r, alpha_vgm, n_vgm, L_vgm, sc_vgm, fc_vgm, &
#endif
               psi0,rootfr_p(:,p),&
               dz_soisno,t_soisno,wliq_soisno,rootr_p(:,i),etrc_p(i),rstfac_p(i))

            ! fraction of sunlit and shaded leaves of canopy
            fsun_p(i) = ( 1. - exp(-min(extkb_p(i)*lai_p(i),40.))) &
                      / max( min(extkb_p(i)*lai_p(i),40.), 1.e-6 )

            IF (coszen<=0.0 .or. sabv_p(i)<1.) fsun_p(i) = 0.5

            laisun_p(i)    = lai_p(i)*fsun_p(i)
            laisha_p(i)    = lai_p(i)*(1-fsun_p(i))
            rstfacsun_p(i) = rstfac_p(i)
            rstfacsha_p(i) = rstfac_p(i)
         ENDIF
      ENDDO


IF (DEF_USE_PFT .or. patchclass(ipatch)==CROPLAND) THEN

      DO i = ps, pe
         p = pftclass(i)
         IF (lai_p(i)+sai_p(i) > 1e-6) THEN

            CALL LeafTemperature(ipatch,p,deltim,csoilc,dewmx   ,htvp       ,&
                 lai_p(i)   ,sai_p(i)   ,htop_p(i)  ,hbot_p(i)  ,sqrtdi_p(p),&
                 effcon_p(p),vmax25_p(p),slti_p(p)  ,hlti_p(p)  ,shti_p(p)  ,&
                 hhti_p(p)  ,trda_p(p)  ,trdm_p(p)  ,trop_p(p)  ,gradm_p(p) ,&
                 binter_p(p),extkn_p(p) ,extkb_p(i) ,extkd_p(i) ,forc_hgt_u ,&
                 forc_hgt_t ,forc_hgt_q ,forc_us    ,forc_vs    ,thm        ,&
                 th         ,thv        ,forc_q     ,forc_psrf  ,forc_rhoair,&
                 parsun_p(i),parsha_p(i),sabv_p(i)  ,frl        ,fsun_p(i)  ,&
                 thermk_p(i),rstfacsun_p(i)         ,rstfacsha_p(i)         ,&
                 gssun_p(i) ,gssha_p(i) ,forc_po2m  ,forc_pco2m ,z0h_g      ,&
                 obu_g      ,ustar_g    ,zlnd       ,zsno       ,fsno       ,&
                 sigf_p(i)  ,etrc_p(i)  ,t_grnd     ,qg         ,dqgdT      ,&
                 emg        ,tleaf_p(i) ,ldew_p(i)  ,ldew_rain_p(i),ldew_snow_p(i),&
                 taux_p(i)  ,tauy_p(i)  ,&
                 fseng_p(i) ,fevpg_p(i) ,cgrnd_p(i) ,cgrndl_p(i),cgrnds_p(i),&
                 tref_p(i)  ,qref_p(i)  ,rst_p(i)   ,assim_p(i) ,respc_p(i) ,&
                 fsenl_p(i) ,fevpl_p(i) ,etr_p(i)   ,dlrad_p(i) ,ulrad_p(i) ,&
                 z0m_p(i)   ,zol_p(i)   ,rib_p(i)   ,ustar_p(i) ,qstar_p(i) ,&
                 tstar_p(i) ,fm_p(i)    ,fh_p(i)    ,fq_p(i)    ,rootfr_p(:,p),&
                 kmax_sun_p(p) ,kmax_sha_p(p) ,kmax_xyl_p(p)  ,kmax_root_p(p) ,psi50_sun_p(p),&
                 psi50_sha_p(p),psi50_xyl_p(p),psi50_root_p(p),ck_p(p)        ,vegwp_p(:,i)  ,&
                 gs0sun_p(i)   ,gs0sha_p(i)                                                  ,&
                 assimsun_p(i) ,etrsun_p(i)   ,assimsha_p(i)  ,etrsha_p(i)    ,&
!Ozone stress variables
                 o3coefv_sun_p(i) ,o3coefv_sha_p(i) ,o3coefg_sun_p(i) ,o3coefg_sha_p(i),&
                 lai_old_p(i), o3uptakesun_p(i) ,o3uptakesha_p(i) ,forc_ozone ,&
!end ozone stress variables
                 forc_hpbl                                                  ,&
                 qintr_rain_p(i),qintr_snow_p(i),t_precip,hprl_p(i),smp     ,&
                 hk(1:)      ,hksati(1:),rootr_p(1:,i)                      )

         ELSE

            CALL groundfluxes (zlnd,zsno,forc_hgt_u,forc_hgt_t,forc_hgt_q,forc_hpbl, &
                               forc_us,forc_vs,forc_t,forc_q,forc_rhoair,forc_psrf, &
                               ur,thm,th,thv,t_grnd,qg,dqgdT,htvp, &
                               fsno,cgrnd_p(i),cgrndl_p(i),cgrnds_p(i), &
                               taux_p(i),tauy_p(i),fseng_p(i),fevpg_p(i),tref_p(i),qref_p(i), &
                               z0m_p(i),z0h_g,zol_p(i),rib_p(i),ustar_p(i),&
                               qstar_p(i),tstar_p(i),fm_p(i),fh_p(i),fq_p(i))

            tleaf_p      (i) = forc_t
            laisun_p     (i) = 0.
            laisha_p     (i) = 0.
            ldew_rain_p  (i) = 0.
            ldew_snow_p  (i) = 0.
            ldew_p       (i) = 0.
            rootr_p    (:,i) = 0.
            rstfacsun_p  (i) = 0.
            rstfacsha_p  (i) = 0.
            gssun_p      (i) = 0.
            gssha_p      (i) = 0.
            assimsun_p   (i) = 0.
            etrsun_p     (i) = 0.
            assimsha_p   (i) = 0.
            etrsha_p     (i) = 0.
            rst_p        (i) = 2.0e4
            assim_p      (i) = 0.
            respc_p      (i) = 0.
            fsenl_p      (i) = 0.
            fevpl_p      (i) = 0.
            etr_p        (i) = 0.
            dlrad_p      (i) = frl
            ulrad_p      (i) = frl*(1.-emg) + emg*stefnc*t_grnd**4
            hprl_p       (i) = 0.

            IF(DEF_USE_PLANTHYDRAULICS)THEN
               vegwp_p(:,i) = -2.5e4
            ENDIF
         ENDIF
      ENDDO

ENDIF


IF (DEF_USE_PC .and. patchclass(ipatch)/=CROPLAND) THEN

      ! initialization
      tleaf_p (ps:pe) = forc_t  !???
      rst_p   (ps:pe) = 2.0e4
      assim_p (ps:pe) = 0.
      respc_p (ps:pe) = 0.
      fsenl_p (ps:pe) = 0.
      fevpl_p (ps:pe) = 0.
      etr_p   (ps:pe) = 0.
      hprl_p  (ps:pe) = 0.
      z0m_p   (ps:pe) = (1.-fsno)*zlnd + fsno*zsno

      IF (lai+sai > 1e-6) THEN

         CALL LeafTemperaturePC (ipatch,ps,pe  ,deltim        ,csoilc             ,dewmx              ,&
            htvp            ,pftclass(ps:pe)   ,pftfrac(ps:pe),htop_p(ps:pe)      ,hbot_p(ps:pe)      ,&
            lai_p(ps:pe)    ,sai_p(ps:pe)      ,&
            extkb_p(ps:pe)  ,extkd_p(ps:pe)    ,forc_hgt_u    ,forc_hgt_t         ,forc_hgt_q         ,&
            forc_us         ,forc_vs           ,thm           ,th                 ,thv                ,&
            forc_q          ,forc_psrf         ,forc_rhoair   ,parsun_p(ps:pe)    ,parsha_p(ps:pe)    ,&
            fsun_p(:)       ,sabv_p(:)         ,frl           ,thermk_p(ps:pe)    ,fshade_p(ps:pe)    ,&
            rstfacsun_p(:)  ,rstfacsha_p(:)    ,&
            gssun_p(:)      ,gssha_p(:)        ,forc_po2m     ,forc_pco2m         ,z0h_g       ,obu_g ,&
            ustar_g         ,zlnd              ,zsno          ,fsno               ,sigf_p(ps:pe)      ,&
            etrc_p(:)       ,t_grnd            ,qg            ,dqgdT              ,emg                ,&
            z0m_p(ps:pe)    ,tleaf_p(ps:pe)    ,ldew_p(ps:pe) ,ldew_rain_p(ps:pe) ,ldew_snow_p(ps:pe) ,&
            taux            ,tauy              ,&
            fseng           ,fevpg             ,cgrnd         ,cgrndl             ,cgrnds             ,&
            tref            ,qref              ,rst_p(ps:pe)  ,assim_p(ps:pe)     ,respc_p(ps:pe)     ,&
            fsenl_p(ps:pe)  ,fevpl_p(ps:pe)    ,etr_p(ps:pe)  ,dlrad              ,ulrad              ,&
            z0m             ,zol               ,rib           ,ustar              ,qstar              ,&
            tstar           ,fm                ,fh            ,fq                 ,vegwp_p(:,ps:pe)   ,&
            gs0sun_p(ps:pe) ,gs0sha_p(ps:pe)   ,&
            assimsun_p(:)   ,etrsun_p(:)       ,assimsha_p(:) ,etrsha_p(:)        ,&
!Ozone stress variables
            o3coefv_sun_p(ps:pe) ,o3coefv_sha_p(ps:pe) ,o3coefg_sun_p(ps:pe) ,o3coefg_sha_p(ps:pe) ,&
            lai_old_p(ps:pe)     ,o3uptakesun_p(ps:pe) ,o3uptakesha_p(ps:pe) ,forc_ozone           ,&
!End ozone stress variables
            forc_hpbl                                                                  ,&
            qintr_rain_p(ps:pe) ,qintr_snow_p(ps:pe) ,t_precip  ,hprl_p(:)   ,smp      ,&
            hk(1:)              ,hksati(1:)          ,rootr_p(:,:)                                  )

      ELSE
         laisun_p    (ps:pe) = 0.
         laisha_p    (ps:pe) = 0.
         tleaf_p     (ps:pe) = forc_t
         ldew_rain_p (ps:pe) = 0.
         ldew_snow_p (ps:pe) = 0.
         ldew_p      (ps:pe) = 0.
         rst_p       (ps:pe) = 2.0e4
         assim_p     (ps:pe) = 0.
         respc_p     (ps:pe) = 0.
         fsenl_p     (ps:pe) = 0.
         fevpl_p     (ps:pe) = 0.
         etr_p       (ps:pe) = 0.
         hprl_p      (ps:pe) = 0.
         z0m_p       (ps:pe) = (1.-fsno)*zlnd + fsno*zsno

         IF(DEF_USE_PLANTHYDRAULICS)THEN
            vegwp_p (:,ps:pe) = -2.5e4
         ENDIF
      ENDIF

ENDIF

      laisun        = sum( laisun_p    (ps:pe)*pftfrac(ps:pe) )
      laisha        = sum( laisha_p    (ps:pe)*pftfrac(ps:pe) )
      tleaf         = sum( tleaf_p     (ps:pe)*pftfrac(ps:pe) )
      ldew_rain     = sum( ldew_rain_p (ps:pe)*pftfrac(ps:pe) )
      ldew_snow     = sum( ldew_snow_p (ps:pe)*pftfrac(ps:pe) )
      ldew          = sum( ldew_p      (ps:pe)*pftfrac(ps:pe) )
      ! may have problem with rst, but the same for LC
      rst           = sum( rst_p       (ps:pe)*pftfrac(ps:pe) )
      assim         = sum( assim_p     (ps:pe)*pftfrac(ps:pe) )
      respc         = sum( respc_p     (ps:pe)*pftfrac(ps:pe) )
      fsenl         = sum( fsenl_p     (ps:pe)*pftfrac(ps:pe) )
      fevpl         = sum( fevpl_p     (ps:pe)*pftfrac(ps:pe) )
      etr           = sum( etr_p       (ps:pe)*pftfrac(ps:pe) )
IF (DEF_USE_PFT .or. patchclass(ipatch)==CROPLAND) THEN
      dlrad         = sum( dlrad_p     (ps:pe)*pftfrac(ps:pe) )
      ulrad         = sum( ulrad_p     (ps:pe)*pftfrac(ps:pe) )
      tref          = sum( tref_p      (ps:pe)*pftfrac(ps:pe) )
      qref          = sum( qref_p      (ps:pe)*pftfrac(ps:pe) )
      taux          = sum( taux_p      (ps:pe)*pftfrac(ps:pe) )
      tauy          = sum( tauy_p      (ps:pe)*pftfrac(ps:pe) )
      fseng         = sum( fseng_p     (ps:pe)*pftfrac(ps:pe) )
      fevpg         = sum( fevpg_p     (ps:pe)*pftfrac(ps:pe) )
      cgrnd         = sum( cgrnd_p     (ps:pe)*pftfrac(ps:pe) )
      cgrndl        = sum( cgrndl_p    (ps:pe)*pftfrac(ps:pe) )
      cgrnds        = sum( cgrnds_p    (ps:pe)*pftfrac(ps:pe) )
      z0m           = sum( z0m_p       (ps:pe)*pftfrac(ps:pe) )
      zol           = sum( zol_p       (ps:pe)*pftfrac(ps:pe) )
      rib           = sum( rib_p       (ps:pe)*pftfrac(ps:pe) )
      ustar         = sum( ustar_p     (ps:pe)*pftfrac(ps:pe) )
      qstar         = sum( qstar_p     (ps:pe)*pftfrac(ps:pe) )
      tstar         = sum( tstar_p     (ps:pe)*pftfrac(ps:pe) )
      fm            = sum( fm_p        (ps:pe)*pftfrac(ps:pe) )
      fh            = sum( fh_p        (ps:pe)*pftfrac(ps:pe) )
      fq            = sum( fq_p        (ps:pe)*pftfrac(ps:pe) )
ENDIF
      rstfacsun_out = sum( rstfacsun_p (ps:pe)*pftfrac(ps:pe) )
      rstfacsha_out = sum( rstfacsha_p (ps:pe)*pftfrac(ps:pe) )
      gssun_out     = sum( gssun_p     (ps:pe)*pftfrac(ps:pe) )
      gssha_out     = sum( gssha_p     (ps:pe)*pftfrac(ps:pe) )
      assimsun_out  = sum( assimsun_p  (ps:pe)*pftfrac(ps:pe) )
      etrsun_out    = sum( etrsun_p    (ps:pe)*pftfrac(ps:pe) )
      assimsha_out  = sum( assimsha_p  (ps:pe)*pftfrac(ps:pe) )
      etrsha_out    = sum( etrsha_p    (ps:pe)*pftfrac(ps:pe) )
      hprl          = sum( hprl_p      (ps:pe)*pftfrac(ps:pe) )

      IF(DEF_USE_PLANTHYDRAULICS)THEN
         DO j = 1, nvegwcs
            vegwp(j) = sum( vegwp_p(j,ps:pe)*pftfrac(ps:pe) )
         ENDDO

         IF (abs(etr) > 0.) THEN
            DO j = 1, nl_soil
               rootr(j) = sum(rootr_p(j,ps:pe)*pftfrac(ps:pe))
            ENDDO
         ENDIF
      ELSE
         IF (abs(etr) > 0.) THEN
            DO j = 1, nl_soil
               rootr(j) = sum(rootr_p(j,ps:pe)*etr_p(ps:pe)*pftfrac(ps:pe)) / etr
            ENDDO
         ENDIF
      ENDIF

      deallocate ( rootr_p     )
      deallocate ( etrc_p      )
      deallocate ( rstfac_p    )
      deallocate ( rstfacsun_p )
      deallocate ( rstfacsha_p )
      deallocate ( gssun_p     )
      deallocate ( gssha_p     )
      deallocate ( fsun_p      )
      deallocate ( sabv_p      )
IF (DEF_USE_PFT .or. patchclass(ipatch)==CROPLAND) THEN
      deallocate ( cgrnd_p     )
      deallocate ( cgrnds_p    )
      deallocate ( cgrndl_p    )
      deallocate ( dlrad_p     )
      deallocate ( ulrad_p     )
      deallocate ( zol_p       )
      deallocate ( rib_p       )
      deallocate ( ustar_p     )
      deallocate ( qstar_p     )
      deallocate ( tstar_p     )
      deallocate ( fm_p        )
      deallocate ( fh_p        )
      deallocate ( fq_p        )
ENDIF
      deallocate ( hprl_p      )
      deallocate ( assimsun_p  )
      deallocate ( etrsun_p    )
      deallocate ( assimsha_p  )
      deallocate ( etrsha_p    )

#endif

! For patchtype/=0, not a soil patch
ELSE
      CALL groundfluxes (zlnd,zsno,forc_hgt_u,forc_hgt_t,forc_hgt_q,forc_hpbl, &
                         forc_us,forc_vs,forc_t,forc_q,forc_rhoair,forc_psrf, &
                         ur,thm,th,thv,t_grnd,qg,dqgdT,htvp, &
                         fsno,cgrnd,cgrndl,cgrnds, &
                         taux,tauy,fseng,fevpg,tref,qref, &
                         z0m_g,z0h_g,zol_g,rib_g,ustar_g,qstar_g,tstar_g,fm_g,fh_g,fq_g)

      ! SAVE variables for bareground case
      obu_g = forc_hgt_u / zol_g

!=======================================================================
! [4] Canopy temperature, fluxes from the canopy
!=======================================================================

      sabv = sabvsun + sabvsha

      IF (lai+sai > 1e-6) THEN

         ! soil water stress factor on stomatal resistance
         CALL eroot (nl_soil,trsmx0,porsl,&
#ifdef Campbell_SOIL_MODEL
            bsw,&
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
            theta_r, alpha_vgm, n_vgm, L_vgm, sc_vgm, fc_vgm, &
#endif
            psi0,rootfr,dz_soisno,t_soisno,wliq_soisno,rootr,etrc,rstfac)

         ! fraction of sunlit and shaded leaves of canopy
         fsun = ( 1. - exp(-min(extkb*lai,40.))) / max( min(extkb*lai,40.), 1.e-6 )

         IF (coszen<=0.0 .or. sabv<1.) fsun = 0.5

         laisun = lai*fsun
         laisha = lai*(1-fsun)
         rstfacsun_out = rstfac
         rstfacsha_out = rstfac

         CALL LeafTemperature(ipatch,1,deltim,csoilc,dewmx     ,htvp       ,&
                 lai        ,sai        ,htop      ,hbot       ,sqrtdi     ,&
                 effcon     ,vmax25     ,slti      ,hlti       ,shti       ,&
                 hhti       ,trda       ,trdm      ,trop       ,gradm      ,&
                 binter     ,extkn      ,extkb     ,extkd      ,forc_hgt_u ,&
                 forc_hgt_t ,forc_hgt_q ,forc_us   ,forc_vs    ,thm        ,&
                 th         ,thv        ,forc_q    ,forc_psrf  ,forc_rhoair,&
                 parsun     ,parsha     ,sabv      ,frl        ,fsun       ,&
                 thermk     ,rstfacsun_out         ,rstfacsha_out          ,&
                 gssun_out  ,gssha_out  ,forc_po2m ,forc_pco2m ,z0h_g      ,&
                 obu_g      ,ustar_g    ,zlnd      ,zsno       ,fsno       ,&
                 sigf       ,etrc       ,t_grnd    ,qg         ,dqgdT      ,&
                 emg        ,tleaf      ,ldew      ,ldew_rain  ,ldew_snow  ,&
                 taux       ,tauy       ,&
                 fseng      ,fevpg      ,cgrnd     ,cgrndl     ,cgrnds     ,&
                 tref       ,qref       ,rst       ,assim      ,respc      ,&
                 fsenl      ,fevpl      ,etr       ,dlrad      ,ulrad      ,&
                 z0m        ,zol        ,rib       ,ustar      ,qstar      ,&
                 tstar      ,fm         ,fh        ,fq         ,rootfr     ,&
                 kmax_sun    ,kmax_sha  ,kmax_xyl  ,kmax_root  ,psi50_sun  ,&
                 psi50_sha   ,psi50_xyl ,psi50_root,ck         ,vegwp      ,&
                 gs0sun      ,gs0sha                                       ,&
                 assimsun_out,etrsun_out,assimsha_out          ,etrsha_out ,&
! Ozone stress variables
                 o3coefv_sun ,o3coefv_sha ,o3coefg_sun ,o3coefg_sha ,&
                 lai_old     ,o3uptakesun ,o3uptakesha ,forc_ozone  ,&
! End ozone stress variables
                 forc_hpbl                                                 ,&
                 qintr_rain  ,qintr_snow,t_precip  ,hprl       ,smp        ,&
                 hk(1:)      ,hksati(1:),rootr(1:)                         )
      ENDIF

    ! equate canopy temperature to air over bareland.
    ! required as sigf=0 carried over to next time step
      IF (lai+sai <= 1e-6) THEN
         tleaf         = forc_t
         laisun        = 0.
         laisha        = 0.
         ldew_rain     = 0.
         ldew_snow     = 0.
         ldew          = 0.
         rstfacsun_out = 0.
         rstfacsha_out = 0.
         IF(DEF_USE_PLANTHYDRAULICS)THEN
            vegwp = -2.5e4
         ENDIF
      ENDIF

ENDIF

!=======================================================================
! [5] Gound temperature
!=======================================================================

      CALL groundtem (patchtype,lb,nl_soil,deltim,&
                      capr,cnfac,vf_quartz,vf_gravels,vf_om,vf_sand,wf_gravels,wf_sand,&
                      porsl,psi0,&
#ifdef Campbell_SOIL_MODEL
                      bsw,&
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
                      theta_r, alpha_vgm, n_vgm, L_vgm,&
                      sc_vgm , fc_vgm,&
#endif
                      csol,k_solids,dksatu,dksatf,dkdry,&
                      BA_alpha,BA_beta,&
                      sigf,dz_soisno,z_soisno,zi_soisno,&
                      t_soisno,wice_soisno,wliq_soisno,scv,snowdp,&
                      frl,dlrad,sabg,sabg_lyr,fseng,fevpg,cgrnd,htvp,emg,&
                      imelt,snofrz,sm,xmf,fact,pg_rain,pg_snow,t_precip)

!=======================================================================
! [6] Correct fluxes to present soil temperature
!=======================================================================

      t_grnd = t_soisno(lb)
      tinc   = t_soisno(lb) - t_soisno_bef(lb)
      fseng  = fseng + tinc*cgrnds
      fevpg  = fevpg + tinc*cgrndl

! calculation of evaporative potential; flux in kg m-2 s-1.
! egidif holds the excess energy IF all water is evaporated
! during the timestep.  this energy is later added to the sensible heat flux.

      egsmax = (wice_soisno(lb)+wliq_soisno(lb)) / deltim
      egidif = max( 0., fevpg - egsmax )
      fevpg  = min( fevpg, egsmax )
      fseng  = fseng + htvp*egidif

! total fluxes to atmosphere
      fsena  = fsenl + fseng
      fevpa  = fevpl + fevpg
      lfevpa = hvap*fevpl + htvp*fevpg   ! W/m^2 (accouting for sublimation)

      qseva = 0.
      qsubl = 0.
      qfros = 0.
      qsdew = 0.

      IF (fevpg >= 0.) THEN
! not allow for sublimation in melting (melting ==> evap. ==> sublimation)
         qseva = min(wliq_soisno(lb)/deltim, fevpg)
         qsubl = fevpg - qseva
      ELSE
         IF (t_grnd < tfrz) THEN
            qfros = abs(fevpg)
         ELSE
            qsdew = abs(fevpg)
         ENDIF
      ENDIF

! ground heat flux
      fgrnd = sabg + dlrad*emg &
            - emg*stefnc*t_grnd_bef**4 &
            - emg*stefnc*t_grnd_bef**3*(4.*tinc) &
            - (fseng+fevpg*htvp) + cpliq * pg_rain * (t_precip - t_grnd) &
            + cpice * pg_snow * (t_precip - t_grnd)

! outgoing long-wave radiation from canopy + ground
      olrg = ulrad &
! for conservation we put the increase of ground longwave to outgoing
           + 4.*emg*stefnc*t_grnd_bef**3*tinc

! averaged bulk surface emissivity
      olrb = stefnc*t_grnd_bef**3*(4.*tinc)
      olru = ulrad + emg*olrb
      olrb = ulrad + olrb
      emis = olru / olrb

! radiative temperature
      IF (olrg < 0) THEN
         print *, "olrg abnormal value: ",ipatch, olrg, tinc, ulrad
         write(6,*) ipatch,errore,sabv,sabg,frl,olrg,fsenl,fseng,hvap*fevpl,htvp*fevpg,xmf,fgrnd
      ENDIF

      trad = (olrg/stefnc)**0.25

! additonal variables required by WRF and RSM model
      IF (lai+sai <= 1e-6) THEN
         ustar = ustar_g
         tstar = tstar_g
         qstar = qstar_g
         rib   = rib_g
         zol   = zol_g
         z0m   = z0m_g
         fm    = fm_g
         fh    = fh_g
         fq    = fq_g
      ELSE
         ustar = ustar
         tstar = tstar
         qstar = qstar
         rib   = rib
         zol   = zol
         z0m   = z0m
         fm    = fm
         fh    = fh
         fq    = fq
      ENDIF

!=======================================================================
! [7] energy balance error
!=======================================================================

      ! one way to check energy
      errore = sabv + sabg + frl - olrg - fsena - lfevpa - fgrnd

      ! another way to check energy
      errore = sabv + sabg + frl - olrg - fsena - lfevpa - xmf + hprl + &
               cpliq * pg_rain * (t_precip - t_grnd) + cpice * pg_snow * (t_precip - t_grnd)
      DO j = lb, nl_soil
         errore = errore - (t_soisno(j)-t_soisno_bef(j))/fact(j)
      ENDDO

#if (defined CoLMDEBUG)
      IF (abs(errore) > .5) THEN
      write(6,*) 'THERMAL.F90: energy balance violation'
      write(6,*) ipatch,errore,sabv,sabg,frl,olrg,fsenl,fseng,hvap*fevpl,htvp*fevpg,xmf,hprl
      STOP
      ENDIF
100   format(10(f15.3))
#endif

  END SUBROUTINE THERMAL

END MODULE MOD_Thermal
