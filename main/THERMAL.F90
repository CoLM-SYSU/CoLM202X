#include <define.h>

 SUBROUTINE THERMAL (ipatch      ,patchtype   ,lb          ,deltim     ,&
                     trsmx0      ,zlnd        ,zsno        ,csoilc     ,&
                     dewmx       ,capr        ,cnfac       ,csol       ,&
                     porsl       ,psi0        ,bsw         ,dkdry      ,&
                     dksatu      ,lai         ,laisun      ,laisha     ,&
                     sai         ,htop        ,hbot        ,sqrtdi     ,&
                     rootfr      ,rstfac      ,effcon      ,vmax25     ,&
                     slti        ,hlti        ,shti        ,hhti       ,&
                     trda        ,trdm        ,trop        ,gradm      ,&
                     binter      ,extkn       ,forc_hgt_u  ,forc_hgt_t ,&
                     forc_hgt_q  ,forc_us     ,forc_vs     ,forc_t     ,&
                     forc_q      ,forc_rhoair ,forc_psrf   ,forc_pco2m ,&
                     forc_po2m   ,coszen      ,parsun      ,parsha     ,&
                     sabvsun     ,sabvsha     ,sabg        ,frl        ,&
                     extkb       ,extkd       ,thermk      ,fsno       ,&
                     sigf        ,dz_soisno   ,z_soisno    ,zi_soisno  ,&
                     tleaf       ,t_soisno    ,wice_soisno ,wliq_soisno,&
                     ldew        ,scv         ,snowdp      ,imelt      ,&
                     taux        ,tauy        ,fsena       ,fevpa      ,&
                     lfevpa      ,fsenl       ,fevpl       ,etr        ,&
                     fseng       ,fevpg       ,olrg        ,fgrnd      ,&
                     rootr       ,qseva       ,qsdew       ,qsubl      ,&
                     qfros       ,sm          ,tref        ,qref       ,&
                     trad        ,rst         ,assim       ,respc      ,&
                     errore      ,emis        ,z0m         ,zol        ,&
                     rib         ,ustar       ,qstar       ,tstar      ,&
                     fm          ,fh          ,fq                       )

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
!              eroot                      |dewfraction
!              LeafTemp   |               |qsadv
!              LeafTempPC |  ---------->  |moninobukini
!                                         |moninobuk
!                                         |ASSIM_STOMATA_conductance
!
!              groundTem     ---------->   meltf
!                               
!=======================================================================

  USE precision
  USE GlobalVars
  USE PFT_Const
  USE PhysicalConstants, only: denh2o,roverg,hvap,hsub,rgas,cpair,&
                               stefnc,denice,tfrz,vonkar,grav 
  USE FRICTION_VELOCITY
  USE LEAF_temperature
  USE LEAF_temperature_PC
  USE MOD_PFTimeInvars
  USE MOD_PCTimeInvars
  USE MOD_PFTimeVars
  USE MOD_PCTimeVars
  USE MOD_1D_PFTFluxes
  USE MOD_1D_PCFluxes

  IMPLICIT NONE
 
!---------------------Argument------------------------------------------

  INTEGER, intent(in) :: &
        ipatch,      &! patch index
        lb,          &! lower bound of array 
        patchtype     ! land water TYPE (0=soil, 1=urban or built-up, 2=wetland,
                      !                  3=glacier/ice sheet, 4=land water bodies)

  REAL(r8), intent(in) :: &
        deltim,      &! model time step [second]
        trsmx0,      &! max transpiration for moist soil+100% veg.  [mm/s]
        zlnd,        &! roughness length for soil [m]
        zsno,        &! roughness length for snow [m]
        csoilc,      &! drag coefficient for soil under canopy [-]
        dewmx,       &! maximum dew
        capr,        &! tuning factor to turn first layer T into surface T
        cnfac,       &! Crank Nicholson factor between 0 and 1

        ! soil physical parameters
        csol(1:nl_soil),   &! heat capacity of soil solids [J/(m3 K)]
        porsl(1:nl_soil),  &! soil porosity [-]
        psi0(1:nl_soil),   &! soil water suction, negative potential [m]
        bsw(1:nl_soil),    &! clapp and hornbereger "b" parameter [-]
        dkdry(1:nl_soil),  &! thermal conductivity of dry soil [W/m-K]
        dksatu(1:nl_soil), &! thermal conductivity of saturated soil [W/m-K]

        ! vegetation parameters
        lai,         &! adjusted leaf area index for seasonal variation [-]
        sai,         &! stem area index  [-]
        htop,        &! canopy crown top height [m]
        hbot,        &! canopy crown bottom height [m]
        sqrtdi,      &! inverse sqrt of leaf dimension [m**-0.5]
        rootfr(1:nl_soil),&! root fraction 
       
        effcon,      &! quantum efficiency of RuBP regeneration (mol CO2/mol quanta)
        vmax25,      &! maximum carboxylation rate at 25 C at canopy top
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

        ! state variables (2)
  REAL(r8), intent(inout) :: &
        tleaf,       &! shaded leaf temperature [K]
        t_soisno(lb:nl_soil),   &! soil temperature [K]
        wice_soisno(lb:nl_soil),&! ice lens [kg/m2]
        wliq_soisno(lb:nl_soil),&! liqui water [kg/m2]
        ldew,        &! depth of water on foliage [kg/(m2 s)] 
        scv,         &! snow cover, water equivalent [mm, kg/m2]
        snowdp        ! snow depth [m]

  INTEGER, intent(out) :: & 
       imelt(lb:nl_soil) ! flag for melting or freezing [-]
  
  REAL(r8), intent(out) :: &
       laisun,       &! sunlit leaf area index
       laisha,       &! shaded leaf area index
       rstfac         ! factor of soil water stress 
 
        ! Output fluxes
  REAL(r8), intent(out) :: &
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

  INTEGER i,j

  REAL(r8) :: &
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
       par,          &! PAR absorbed by canopy [W/m2]
       qg,           &! ground specific humidity [kg/kg]
       qsatg,        &! saturated humidity [kg/kg]
       qsatgdT,      &! d(qsatg)/dT
       qred,         &! soil surface relative humidity
       sabv,         &! solar absorbed by canopy [W/m2]
       thm,          &! intermediate variable (forc_t+0.0098*forc_hgt_t)
       th,           &! potential temperature (kelvin)
       thv,          &! virtual potential temperature (kelvin)
       t_grnd,       &! ground surface temperature [K]
       t_grnd_bef,   &! ground surface temperature [K]
       t_soisno_bef(lb:nl_soil), &! soil/snow temperature before update
       tinc,         &! temperature difference of two time step
       ur,           &! wind speed at reference height [m/s]
       ulrad,        &! upward longwave radiation above the canopy [W/m2]
       wice0(lb:nl_soil),&! ice mass from previous time-step
       wliq0(lb:nl_soil),&! liquid mass from previous time-step
       wx,           &! patitial volume of ice and water of surface layer
       xmf            ! total latent heat of phase change of ground water

  REAL(r8) :: z0m_g,z0h_g,zol_g,obu_g,rib_g,ustar_g,qstar_g,tstar_g
  REAL(r8) :: fm10m,fm_g,fh_g,fq_g,fh2m,fq2m,um,obu
     
  INTEGER p, ps, pe, pc

#ifdef PFT_CLASSIFICATION
  REAL(r8), allocatable :: rootr_p (:,:)
  REAL(r8), allocatable :: etrc_p  (:)
  REAL(r8), allocatable :: rstfac_p(:)
  REAL(r8), allocatable :: laisun_p(:)
  REAL(r8), allocatable :: laisha_p(:)
  REAL(r8), allocatable :: fsun_p  (:)
  REAL(r8), allocatable :: sabv_p  (:)
  REAL(r8), allocatable :: cgrnd_p (:)
  REAL(r8), allocatable :: cgrnds_p(:)
  REAL(r8), allocatable :: cgrndl_p(:)
  REAL(r8), allocatable :: dlrad_p (:)
  REAL(r8), allocatable :: ulrad_p (:)
  REAL(r8), allocatable :: zol_p   (:)
  REAL(r8), allocatable :: rib_p   (:)
  REAL(r8), allocatable :: ustar_p (:)
  REAL(r8), allocatable :: qstar_p (:)
  REAL(r8), allocatable :: tstar_p (:)
  REAL(r8), allocatable :: fm_p    (:)
  REAL(r8), allocatable :: fh_p    (:)
  REAL(r8), allocatable :: fq_p    (:)
#endif

#ifdef PC_CLASSIFICATION
  REAL(r8) :: rootr_c (nl_soil,0:N_PFT-1)
  REAL(r8) :: etrc_c  (0:N_PFT-1)
  REAL(r8) :: rstfac_c(0:N_PFT-1)
  REAL(r8) :: laisun_c(0:N_PFT-1)
  REAL(r8) :: laisha_c(0:N_PFT-1)
  REAL(r8) :: fsun_c  (0:N_PFT-1)
  REAL(r8) :: sabv_c  (0:N_PFT-1)
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

         psit = psi0(1) * fac ** (- bsw(1) )   !psit = max(smpmin, psit)
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
#if(defined USGS_CLASSIFICATION || defined IGBP_CLASSIFICATION)
      CALL groundfluxes (zlnd,zsno,forc_hgt_u,forc_hgt_t,forc_hgt_q, &
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

         ! soil water strees factor on stomatal resistance
         CALL eroot (nl_soil,trsmx0,porsl,bsw,psi0,rootfr,&
                     dz_soisno,t_soisno,wliq_soisno,rootr,etrc,rstfac)

         ! fraction of sunlit and shaded leaves of canopy
         fsun = ( 1. - exp(-min(extkb*lai,40.))) / max( min(extkb*lai,40.), 1.e-6 )

         IF (coszen<=0.0 .or. sabv<1.) fsun = 0.5
         
         laisun = lai*fsun
         laisha = lai*(1-fsun)

         CALL LeafTemp (ipatch,deltim   ,csoilc    ,dewmx      ,htvp       ,&
                 lai        ,sai        ,htop      ,hbot       ,sqrtdi     ,&
                 effcon     ,vmax25     ,slti      ,hlti       ,shti       ,&
                 hhti       ,trda       ,trdm      ,trop       ,gradm      ,&
                 binter     ,extkn      ,extkb     ,extkd      ,forc_hgt_u ,&
                 forc_hgt_t ,forc_hgt_q ,forc_us   ,forc_vs    ,thm        ,&
                 th         ,thv        ,forc_q    ,forc_psrf  ,forc_rhoair,&
                 parsun     ,parsha     ,sabv      ,frl        ,fsun       ,&
                 thermk     ,rstfac     ,forc_po2m ,forc_pco2m ,z0h_g      ,&
                 obu_g      ,ustar_g    ,zlnd      ,zsno       ,fsno       ,&
                 sigf       ,etrc       ,t_grnd    ,qg         ,dqgdT      ,&
                 emg        ,tleaf      ,ldew      ,taux       ,tauy       ,&
                 fseng      ,fevpg      ,cgrnd     ,cgrndl     ,cgrnds     ,&
                 tref       ,qref       ,rst       ,assim      ,respc      ,&
                 fsenl      ,fevpl      ,etr       ,dlrad      ,ulrad      ,&
                 z0m        ,zol        ,rib       ,ustar      ,qstar      ,&
                 tstar      ,fm         ,fh        ,fq                      )
      ENDIF

    ! equate canopy temperature to air over bareland.
    ! required as sigf=0 carried over to next time step
      IF (lai+sai <= 1e-6) THEN
         tleaf  = forc_t
         laisun = 0.
         laisha = 0.
         ldew   = 0.
         rstfac = 0.
      ENDIF
#endif


!=======================================================================
#ifdef PFT_CLASSIFICATION
      
      ps = patch_pft_s(ipatch)      
      pe = patch_pft_e(ipatch)

      allocate ( rootr_p (nl_soil, ps:pe) )
      allocate ( etrc_p  (ps:pe) )
      allocate ( rstfac_p(ps:pe) ) 
      allocate ( laisun_p(ps:pe) ) 
      allocate ( laisha_p(ps:pe) ) 
      allocate ( fsun_p  (ps:pe) )
      allocate ( sabv_p  (ps:pe) )
      allocate ( cgrnd_p (ps:pe) ) 
      allocate ( cgrnds_p(ps:pe) )
      allocate ( cgrndl_p(ps:pe) )
      allocate ( dlrad_p (ps:pe) )
      allocate ( ulrad_p (ps:pe) )
      allocate ( zol_p   (ps:pe) )
      allocate ( rib_p   (ps:pe) )
      allocate ( ustar_p (ps:pe) )
      allocate ( qstar_p (ps:pe) )
      allocate ( tstar_p (ps:pe) )
      allocate ( fm_p    (ps:pe) )
      allocate ( fh_p    (ps:pe) )
      allocate ( fq_p    (ps:pe) )

      ! always DO CALL groundfluxes
      CALL groundfluxes (zlnd,zsno,forc_hgt_u,forc_hgt_t,forc_hgt_q, &
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

            CALL eroot (nl_soil,trsmx0,porsl,bsw,psi0,rootfr_p(:,p),&
               dz_soisno,t_soisno,wliq_soisno,rootr_p(:,i),etrc_p(i),rstfac_p(i))

            ! fraction of sunlit and shaded leaves of canopy
            fsun_p(i) = ( 1. - exp(-min(extkb_p(i)*lai_p(i),40.))) &
                      / max( min(extkb_p(i)*lai_p(i),40.), 1.e-6 )

            IF (coszen<=0.0 .or. sabv_p(i)<1.) fsun_p(i) = 0.5

            laisun_p(i) = lai_p(i)*fsun_p(i)
            laisha_p(i) = lai_p(i)*(1-fsun_p(i))

            CALL LeafTemp (ipatch,deltim,csoilc     ,dewmx      ,htvp       ,&
                 lai_p(i)   ,sai_p(i)   ,htop_p(i)  ,hbot_p(i)  ,sqrtdi_p(p),&
                 effcon_p(p),vmax25_p(p),slti_p(p)  ,hlti_p(p)  ,shti_p(p)  ,&
                 hhti_p(p)  ,trda_p(p)  ,trdm_p(p)  ,trop_p(p)  ,gradm_p(p) ,&
                 binter_p(p),extkn_p(p) ,extkb_p(i) ,extkd_p(i) ,forc_hgt_u ,&
                 forc_hgt_t ,forc_hgt_q ,forc_us    ,forc_vs    ,thm        ,&
                 th         ,thv        ,forc_q     ,forc_psrf  ,forc_rhoair,&
                 parsun_p(i),parsha_p(i),sabv_p(i)  ,frl        ,fsun_p(i)  ,&
                 thermk_p(i),rstfac_p(i),forc_po2m  ,forc_pco2m ,z0h_g      ,&
                 obu_g      ,ustar_g    ,zlnd       ,zsno       ,fsno       ,&
                 sigf_p(i)  ,etrc_p(i)  ,t_grnd     ,qg         ,dqgdT      ,&
                 emg        ,tleaf_p(i) ,ldew_p(i)  ,taux_p(i)  ,tauy_p(i)  ,&
                 fseng_p(i) ,fevpg_p(i) ,cgrnd_p(i) ,cgrndl_p(i),cgrnds_p(i),&
                 tref_p(i)  ,qref_p(i)  ,rst_p(i)   ,assim_p(i) ,respc_p(i) ,&
                 fsenl_p(i) ,fevpl_p(i) ,etr_p(i)   ,dlrad_p(i) ,ulrad_p(i) ,&
                 z0m_p(i)   ,zol_p(i)   ,rib_p(i)   ,ustar_p(i) ,qstar_p(i) ,&
                 tstar_p(i) ,fm_p(i)    ,fh_p(i)    ,fq_p(i)                 )

         ELSE 

            CALL groundfluxes (zlnd,zsno,forc_hgt_u,forc_hgt_t,forc_hgt_q, &
                               forc_us,forc_vs,forc_t,forc_q,forc_rhoair,forc_psrf, &
                               ur,thm,th,thv,t_grnd,qg,dqgdT,htvp, &
                               fsno,cgrnd_p(i),cgrndl_p(i),cgrnds_p(i), &
                               taux_p(i),tauy_p(i),fseng_p(i),fevpg_p(i),tref_p(i),qref_p(i), &
            z0m_p(i),z0h_g,zol_p(i),rib_p(i),ustar_p(i),qstar_p(i),tstar_p(i),fm_p(i),fh_p(i),fq_p(i))

            tleaf_p(i)   = forc_t
            laisun_p(i)  = 0.
            laisha_p(i)  = 0.
            ldew_p(i)    = 0.
            rootr_p(:,i) = 0.
            rstfac_p(i)  = 0.
            rst_p(i)     = 2.0e4
            assim_p(i)   = 0.
            respc_p(i)   = 0.
            fsenl_p(i)   = 0.
            fevpl_p(i)   = 0.
            etr_p(i)     = 0.
            dlrad_p(i)   = frl  
            ulrad_p(i)   = frl*(1.-emg) + emg*stefnc*t_grnd**4

         ENDIF
      ENDDO 
   
      laisun = sum( laisun_p(ps:pe)*pftfrac(ps:pe) )
      laisha = sum( laisha_p(ps:pe)*pftfrac(ps:pe) )
      dlrad  = sum( dlrad_p (ps:pe)*pftfrac(ps:pe) )
      ulrad  = sum( ulrad_p (ps:pe)*pftfrac(ps:pe) )
      tleaf  = sum( tleaf_p (ps:pe)*pftfrac(ps:pe) )
      ldew   = sum( ldew_p  (ps:pe)*pftfrac(ps:pe) )
      tref   = sum( tref_p  (ps:pe)*pftfrac(ps:pe) )
      qref   = sum( qref_p  (ps:pe)*pftfrac(ps:pe) )
      ! may have problem with rst, but the same for LC
      rst    = sum( rst_p   (ps:pe)*pftfrac(ps:pe) )
      assim  = sum( assim_p (ps:pe)*pftfrac(ps:pe) )
      respc  = sum( respc_p (ps:pe)*pftfrac(ps:pe) )
      taux   = sum( taux_p  (ps:pe)*pftfrac(ps:pe) )
      tauy   = sum( tauy_p  (ps:pe)*pftfrac(ps:pe) )
      fseng  = sum( fseng_p (ps:pe)*pftfrac(ps:pe) )
      fevpg  = sum( fevpg_p (ps:pe)*pftfrac(ps:pe) )
      cgrnd  = sum( cgrnd_p (ps:pe)*pftfrac(ps:pe) )
      cgrndl = sum( cgrndl_p(ps:pe)*pftfrac(ps:pe) )
      cgrnds = sum( cgrnds_p(ps:pe)*pftfrac(ps:pe) )
      fsenl  = sum( fsenl_p (ps:pe)*pftfrac(ps:pe) )
      fevpl  = sum( fevpl_p (ps:pe)*pftfrac(ps:pe) )
      etr    = sum( etr_p   (ps:pe)*pftfrac(ps:pe) )
      z0m    = sum( z0m_p   (ps:pe)*pftfrac(ps:pe) )
      zol    = sum( zol_p   (ps:pe)*pftfrac(ps:pe) )
      rib    = sum( rib_p   (ps:pe)*pftfrac(ps:pe) )
      ustar  = sum( ustar_p (ps:pe)*pftfrac(ps:pe) )
      qstar  = sum( qstar_p (ps:pe)*pftfrac(ps:pe) )
      tstar  = sum( tstar_p (ps:pe)*pftfrac(ps:pe) )
      fm     = sum( fm_p    (ps:pe)*pftfrac(ps:pe) )
      fh     = sum( fh_p    (ps:pe)*pftfrac(ps:pe) )
      fq     = sum( fq_p    (ps:pe)*pftfrac(ps:pe) )

      IF (etr > 0.) THEN
         DO j = 1, nl_soil
            rootr(j) = sum(rootr_p(j,ps:pe)*etr_p(ps:pe)*pftfrac(ps:pe)) / etr
         ENDDO 
      ENDIF

      rstfac = sum( rstfac_p(ps:pe)*pftfrac(ps:pe) )

      deallocate ( rootr_p  )
      deallocate ( etrc_p   )
      deallocate ( rstfac_p )
      deallocate ( fsun_p   )
      deallocate ( sabv_p   )
      deallocate ( cgrnd_p  )
      deallocate ( cgrnds_p )
      deallocate ( cgrndl_p )
      deallocate ( dlrad_p  )
      deallocate ( ulrad_p  )
      deallocate ( zol_p    )
      deallocate ( rib_p    )
      deallocate ( ustar_p  )
      deallocate ( qstar_p  )
      deallocate ( tstar_p  )
      deallocate ( fm_p     )
      deallocate ( fh_p     )
      deallocate ( fq_p     )

#endif

!=======================================================================
#ifdef PC_CLASSIFICATION

      pc = patch2pc(ipatch)      
      
      ! always DO CALL groundfluxes
      CALL groundfluxes (zlnd,zsno,forc_hgt_u,forc_hgt_t,forc_hgt_q, &
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

      sabv_c(:) = sabvsun_c(:,pc) + sabvsha_c(:,pc)
      sabv = sabvsun + sabvsha

      DO p = 0, N_PFT-1

         IF (lai_c(p,pc)+sai_c(p,pc) > 1e-6) THEN

            ! soil water strees factor on stomatal resistance
            CALL eroot (nl_soil,trsmx0,porsl,bsw,psi0,rootfr_p(:,p),&
               dz_soisno,t_soisno,wliq_soisno,rootr_c(:,p),etrc_c(p),rstfac_c(p))

            ! fraction of sunlit and shaded leaves of canopy
            fsun_c(p) = ( 1. - exp(-min(extkb_c(p,pc)*lai_c(p,pc),40.))) &
                       / max( min(extkb_c(p,pc)*lai_c(p,pc),40.), 1.e-6 )

! 01/06/2020, yuan: change to 0.5
            IF (coszen<=0.0 .or. sabv_c(p)<1.) fsun_c(p) = 0.5

            laisun_c(p) = lai_c(p,pc)*fsun_c(p)
            laisha_c(p) = lai_c(p,pc)*(1-fsun_c(p))

         ELSE 
            tleaf_c(p,pc) = forc_t
            laisun_c(p)   = 0.
            laisha_c(p)   = 0.
            ldew_c(p,pc)  = 0.
            rootr_c(:,p)  = 0.
            rstfac_c(p)   = 0.
         ENDIF

      ENDDO 
      
      IF (lai+sai > 1e-6) THEN

         !print *, "CALL LeafTempPC 表层温度t_grnd:", t_grnd  !fordebug
         CALL LeafTempPC (ipatch,N_PFT  ,deltim        ,csoilc        ,dewmx         ,&
           htvp          ,pcfrac(:,pc)  ,canlay(:)     ,htop_c(:,pc)  ,hbot_c(:,pc)  ,&
           lai_c(:,pc)   ,sai_c(:,pc)   ,sqrtdi_p(:)   ,effcon_p(:)   ,vmax25_p(:)   ,&
           slti_p(:)     ,hlti_p(:)     ,shti_p(:)     ,hhti_p(:)     ,trda_p(:)     ,&
           trdm_p(:)     ,trop_p(:)     ,gradm_p(:)    ,binter_p(:)   ,extkn_p(:)    ,&
           extkb_c(:,pc) ,extkd_c(:,pc) ,forc_hgt_u    ,forc_hgt_t    ,forc_hgt_q    ,&
           forc_us       ,forc_vs       ,thm           ,th            ,thv           ,&
           forc_q        ,forc_psrf     ,forc_rhoair   ,parsun_c(:,pc),parsha_c(:,pc),&
           fsun_c(:)     ,sabv_c(:)     ,frl           ,thermk_c(:,pc),fshade_c(:,pc),&
           rstfac_c(:)   ,forc_po2m     ,forc_pco2m    ,z0h_g         ,obu_g         ,&
           ustar_g       ,zlnd          ,zsno          ,fsno          ,sigf_c(:,pc)  ,&
           etrc_c(:)     ,t_grnd        ,qg            ,dqgdT         ,emg           ,&
           z0m_c(:,pc)   ,tleaf_c(:,pc) ,ldew_c(:,pc)  ,taux          ,tauy          ,&
           fseng         ,fevpg         ,cgrnd         ,cgrndl        ,cgrnds        ,&
           tref          ,qref          ,rst_c(:,pc)   ,assim_c(:,pc) ,respc_c(:,pc) ,&
           fsenl_c(:,pc) ,fevpl_c(:,pc) ,etr_c(:,pc)   ,dlrad         ,ulrad         ,&
           z0m           ,zol           ,rib           ,ustar         ,qstar         ,&
           tstar         ,fm            ,fh            ,fq                            )
      ELSE 
         laisun_c(:)    = 0.
         laisha_c(:)    = 0.
         tleaf_c (:,pc) = forc_t
         ldew_c  (:,pc) = 0.
         rst_c   (:,pc) = 2.0e4
         assim_c (:,pc) = 0.
         respc_c (:,pc) = 0.
         fsenl_c (:,pc) = 0.
         fevpl_c (:,pc) = 0.
         etr_c   (:,pc) = 0.
      ENDIF

      laisun = sum( laisun_c(:)   *pcfrac(:,pc) )
      laisha = sum( laisha_c(:)   *pcfrac(:,pc) )
      tleaf  = sum( tleaf_c (:,pc)*pcfrac(:,pc) )
      ldew   = sum( ldew_c  (:,pc)*pcfrac(:,pc) )
      rst    = sum( rst_c   (:,pc)*pcfrac(:,pc) )
      assim  = sum( assim_c (:,pc)*pcfrac(:,pc) )
      respc  = sum( respc_c (:,pc)*pcfrac(:,pc) )
      fsenl  = sum( fsenl_c (:,pc)*pcfrac(:,pc) )
      fevpl  = sum( fevpl_c (:,pc)*pcfrac(:,pc) )
      etr    = sum( etr_c   (:,pc)*pcfrac(:,pc) )

      ! loop for each soil layer
      IF (etr > 0.) THEN
         DO j = 1, nl_soil
            rootr(j) = sum(rootr_c(j,:)*etr_c(:,pc)*pcfrac(:,pc)) / etr
         ENDDO 
      ENDIF

      rstfac = sum( rstfac_c(:)*pcfrac(:,pc) )

#endif

! For patchtype/=0, not a soil patch
ELSE 
      CALL groundfluxes (zlnd,zsno,forc_hgt_u,forc_hgt_t,forc_hgt_q, &
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

         ! soil water strees factor on stomatal resistance
         CALL eroot (nl_soil,trsmx0,porsl,bsw,psi0,rootfr,&
                     dz_soisno,t_soisno,wliq_soisno,rootr,etrc,rstfac)

         ! fraction of sunlit and shaded leaves of canopy
         fsun = ( 1. - exp(-min(extkb*lai,40.))) / max( min(extkb*lai,40.), 1.e-6 )

         IF (coszen<=0.0 .or. sabv<1.) fsun = 0.5
         
         laisun = lai*fsun
         laisha = lai*(1-fsun)

         CALL LeafTemp (ipatch,deltim   ,csoilc    ,dewmx      ,htvp       ,&
                 lai        ,sai        ,htop      ,hbot       ,sqrtdi     ,&
                 effcon     ,vmax25     ,slti      ,hlti       ,shti       ,&
                 hhti       ,trda       ,trdm      ,trop       ,gradm      ,&
                 binter     ,extkn      ,extkb     ,extkd      ,forc_hgt_u ,&
                 forc_hgt_t ,forc_hgt_q ,forc_us   ,forc_vs    ,thm        ,&
                 th         ,thv        ,forc_q    ,forc_psrf  ,forc_rhoair,&
                 parsun     ,parsha     ,sabv      ,frl        ,fsun       ,&
                 thermk     ,rstfac     ,forc_po2m ,forc_pco2m ,z0h_g      ,&
                 obu_g      ,ustar_g    ,zlnd      ,zsno       ,fsno       ,&
                 sigf       ,etrc       ,t_grnd    ,qg         ,dqgdT      ,&
                 emg        ,tleaf      ,ldew      ,taux       ,tauy       ,&
                 fseng      ,fevpg      ,cgrnd     ,cgrndl     ,cgrnds     ,&
                 tref       ,qref       ,rst       ,assim      ,respc      ,&
                 fsenl      ,fevpl      ,etr       ,dlrad      ,ulrad      ,&
                 z0m        ,zol        ,rib       ,ustar      ,qstar      ,&
                 tstar      ,fm         ,fh        ,fq                      )
 
      ENDIF

    ! equate canopy temperature to air over bareland.
    ! required as sigf=0 carried over to next time step
      IF (lai+sai <= 1e-6) THEN
         tleaf  = forc_t
         laisun = 0.
         laisha = 0.
         ldew   = 0.
         rstfac = 0.
      ENDIF

ENDIF

!=======================================================================
! [5] Gound temperature
!=======================================================================

      CALL groundtem (patchtype,lb,nl_soil,deltim,&
                      capr,cnfac,csol,porsl,dkdry,dksatu,&
                      sigf,dz_soisno,z_soisno,zi_soisno,&
                      t_soisno,wice_soisno,wliq_soisno,scv,snowdp,&
                      frl,dlrad,sabg,fseng,fevpg,cgrnd,htvp,emg,&
                      imelt,sm,xmf,fact,psi0,bsw)

!=======================================================================
! [6] Correct fluxes to present soil temperature
!=======================================================================

! 问题：为什么tinc更新fseng, fevpg后能继续保持守恒
! 在植被湍流计算过程中的等式保证了能量的平衡，最对土壤温度求导
      t_grnd = t_soisno(lb)
      tinc   = t_soisno(lb) - t_soisno_bef(lb)
      fseng  = fseng + tinc*cgrnds 
      fevpg  = fevpg + tinc*cgrndl

! calculation of evaporative potential; flux in kg m-2 s-1.  
! egidif holds the excess energy IF all water is evaporated
! during the timestep.  this energy is later added to the sensible heat flux.

      egsmax = (wice_soisno(lb)+wliq_soisno(lb)) / deltim
      egidif = max( 0., fevpg - egsmax )
      fevpg = min ( fevpg, egsmax )
      fseng = fseng + htvp*egidif

! total fluxes to atmosphere
      fsena = fsenl + fseng
      fevpa = fevpl + fevpg
      lfevpa= hvap*fevpl + htvp*fevpg   ! W/m^2 (accouting for sublimation)
      
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
            ! 08/18/2021, yuan: bug, t_grnd->t_grnd_bef
            !- emg*stefnc*t_grnd**4 &
            - emg*stefnc*t_grnd_bef**4 &
            - emg*stefnc*t_grnd_bef**3*(4.*tinc) &
            - (fseng+fevpg*htvp)

! outgoing long-wave radiation from canopy + ground
      olrg = ulrad &
! for conservation we put the increase of ground longwave to outgoing
           + 4.*emg*stefnc*t_grnd_bef**3*tinc

! averaged bulk surface emissivity 
! 03/10/2020, yuan: 
      !olrb = stefnc*t_soisno_bef(lb)**3*(4.*tinc)
      olrb = stefnc*t_grnd_bef**3*(4.*tinc)
      olru = ulrad + emg*olrb
      olrb = ulrad + olrb
      emis = olru / olrb

! radiative temperature
      IF (olrg < 0) THEN !fordebug
         print *, ipatch, olrg, tinc, ulrad
         write(6,*) ipatch,errore,sabv,sabg,frl,olrg,fsenl,fseng,hvap*fevpl,htvp*fevpg,xmf
         olrg = 0.
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

      errore = sabv + sabg + frl - olrg - fsena - lfevpa - fgrnd
      !print *, "errore 1:", errore     !one way to check energy
      errore = sabv + sabg + frl - olrg - fsena - lfevpa - xmf
      DO j = lb, nl_soil
         errore = errore - (t_soisno(j)-t_soisno_bef(j))/fact(j)
      ENDDO
      !print *, "errore 2:", errore     !the other way to check energy
      
#if (defined CLMDEBUG)
      IF (abs(errore) > .5) THEN 
      write(6,*) 'THERMAL.F90: energy balance violation'
      write(6,*) ipatch,errore,'=',sabv+sabg+frl-olrg-fsenl-fseng-hvap*fevpl-htvp*fevpg-xmf
      ENDIF
100   format(10(f15.3))
#endif

 END SUBROUTINE THERMAL
