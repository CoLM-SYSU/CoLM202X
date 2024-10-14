#include <define.h>

MODULE  MOD_Lake_FLake

! ================================== code history =====================================
! Description:
!     The interface for the Flake model is adapted from that provided by Mironov et al. 
!    
!     [1] The input and output variables of FLake model are now normalized to be consistent with CoLM-Lake model.
!     [2] The FLake model is layered as follows: snow, ice, water, bottom-sediments.
!     [3] It's a simplification that FLake's snowpack is divided into 5 layers, but its temperature is uniform.
!     [4] The FLake model does not calculate xwliq and xwice, these values are set to 0.
!     [5] Similar to the snow layer, the temperature of the bottom sediment layer is also set to a common temperature.
!  
!           ###################
!           #####   Snow   ####    ---> tmsno, tskin, snwdp
!           ###################
!           $$$$$$$ Ice $$$$$$$    ---> tmice, icedp, icefr
!           -------------------
!           ---- Mix-Layer ----    ---> tmwml,tmmnw, mldp
!           -------------------
!           ~~~~~~~~~~~~~~~~~~~
!           ~~~~~~~~~~~~~~~~~~~
!           ~~~ thermocline ~~~    ---> lktmp
!           ~~~~~~~~~~~~~~~~~~~
!           ~~~~~~~~~~~~~~~~~~~    ---> (tmbot), CTfrac
!           ///////////////////
!           //upper sediments//    ---> t_ssb(tmups), upsdp
!           ===================
!           ===  outer edge ===    ---> T_bs  <=== from the forcing datazxz
!           ===================
!
!+WMEJ TODO: check function performance between SfcFlx_spechum with qsadv
!
! Original author : Dmitrii Mironov, 2005/11/17
!     Current Code Owner: DWD, Dmitrii Mironov
!     Phone:  +49-69-8062 2705
!     Fax:    +49-69-8062 3721
!     E-mail: dmitrii.mironov@dwd.de
!
! Revisions:
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024: Adjusted code style to unify all lake model interfaces
! =====================================================================================

   USE MOD_Precision
   IMPLICIT NONE

!--------------------------------------------------------------------------
CONTAINS
!--------------------------------------------------------------------------

   SUBROUTINE Lake_FLake ( &
               ! "in" arguments
               ! ---------------------------
               nlake     , nsnow     , nsoil     , scwat     ,&
               zopt      , betaopt   , fetchopt  , etaopt    ,&
               dtlak     , xlat      , xlon      , dplak     ,&
               hwref     , htref     , hqref     , tmref     ,&
               usurf     , vsurf     , qmref     , arhos     ,&
               sols      , soll      , solsd     , solld     ,&
               sabg      , lwdns     , psurf     , hpbl      ,&
               crain     , csnow     , lrain     , lsnow     ,&
               dpsed     , T_bs      , zcbcv     , ipatch    ,&
               ! "inout" arguments
               ! ---------------------------
               dzlak     , zlake     , zilak     , lktmp     ,&
               dzssb     , zssb      , zissb     , t_ssb     ,&
               snlay     , scv       , snwag     , snwdp     ,&
               stke1     , tskin     , xwliq     , xwice     ,&
               z0m       , z0h       , z0q       , felak     ,&
               gamma     , etal      , btpri     , frlak     ,&
               tmsno     , tmice     , tmmnw     , tmwml     ,&  
               tmbot     , tmups     , icedp     , mldp      ,&  
               upsdp     , CTfrac    , icefr     , &   
               ! "out" arguments
               ! ---------------------------
               fsena     , fevpa     , lfevpa    , fseng     ,&
               fevpg     , olrg      , fgrnd     , trad      ,&
               qseva     , qsubl     , qsdew     , qfros     ,&
               taux      , tauy      , ustar     , qstar     ,&
               tstar     , emis      , zol       , snwml     ,&
               rib       , ram       , rah       , raw       ,&
               wdm       , t2m       , q2m       , u10m      ,&
               v10m      , fm10m     , fq10m     , fh2m      ,&
               fq2m      , fm        , fq        , fh        ,&
               shfdt     , urban_call)

! ---------------------------------- code history -------------------------------------
! Description:
!     The interface for the Flake model is adapted from that provided by Mironov et al.
!
! Called: (* means optional)
!    -> flake_radtrans     : Scheme for bulk density of newly fallen dry snow
!    *-> FLake_SurfFlux    : Calculate the surface fluxes (FLake Shecme)
!       |-> LakStaParms    : Calculate the stability parameters for the lake
!    *-> CoLML_SurfFlux    : Calculate the surface fluxes (CoLM Shecme)
!    -> flaketem           : Advance FLake variables
!    -> snowdp2lev4lake    : Update the snow layer level
!    -> snowage            : Update the snow age
!
!
! Original author : Dmitrii Mironov, 2005/11/17
!     Current Code Owner: DWD, Dmitrii Mironov
!     Phone:  +49-69-8062 2705
!     Fax:    +49-69-8062 3721
!     E-mail: dmitrii.mironov@dwd.de
!
! Revisions:
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024: Adjusted code style to unify all lake model interfaces
!
! -------------------------------------------------------------------------------------
   USE MOD_Lake_Const, only: stefnc, omega, c_lwrad_emis, h_Ice_min_flk, h_Snow_min_flk,&
                              tpl_rho_I, tpl_rho_w_r, d2r, tfrz
   USE MOD_Lake_CoLML, only: CoLML_SurfFlux
   USE MOD_Lake_Utils, only:  snowage, snowdp2lev4lake
   USE MOD_Lake_Subs, only: LakStaParms
   USE MOD_Namelist, only: DEF_USE_COLML_FLUX_SCHEME
!==============================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in)      :: &
      ipatch                  ,&! patch index
      nlake                   ,&! Maximum number of lake layer
      nsnow                   ,&! Maximum number of snow layer  !positive
      nsoil                   ,&! Maximum number of soil layer
      scwat                   ,&! surface category of water characteristics:
      zopt                    ,&! option for roughness length, 1: constant, 2: Subin et al. (2012)
      betaopt                 ,&! option for betaprime calculation, 1: constant, 2: equation
      fetchopt                ,&! option for fetch length, 1: constant, 2: equation
      etaopt                    ! option for Extinction coefficient calculation, 1: constant, 2: equation

   real(r8), intent(in)     :: &
      dtlak                   ,&! time step [s]
      xlat                    ,&! latitude in degrees
      xlon                    ,&! longitude in degrees
      dplak                   ,&! lake depth [m]
      hwref                   ,&! observational height of wind [m]
      htref                   ,&! observational height of temperature [m]
      hqref                   ,&! observational height of humidity [m]
      tmref                   ,&! temperature at the reference height [K] 
      qmref                   ,&! specific humidity at the reference height [kg/kg]
      usurf                   ,&! wind component in eastward direction [m/s]
      vsurf                   ,&! wind component in northward direction [m/s]
      arhos                   ,&! surface air density [kg m-3]
      sols                    ,&! atm vis direct beam solar rad onto srf [W m-2]
      soll                    ,&! atm nir direct beam solar rad onto srf [W m-2]
      solsd                   ,&! atm vis diffuse solar rad onto srf [W m-2]
      solld                   ,&! atm nir diffuse solar rad onto srf [W m-2]
      sabg                    ,&! solar absorbed by ground  [W/m2]
      lwdns                   ,&! atmospheric infrared (longwave) radiation [W/m2]
      psurf                   ,&! surface pressure [Pa]
      hpbl                    ,&! downward longwave radiation at surface [wm-2]
      zcbcv                   ,&! convective boundary height [m]
      csnow                   ,&! convective snowfall [kg/(m2 s)]
      lsnow                   ,&! large scale snowfall [kg/(m2 s)]
      crain                   ,&! convective rainfall [kg/(m2 s)]
      lrain                   ,&! large scale rainfall [kg/(m2 s)]
      dpsed                   ,&! bottom sediments depth (m)
      T_bs                      ! Temperature at the outer edge of the thermally active layer of the bottom sediments [K]

!  ------------------------- inout variables ---------------------------
   integer, intent(inout)   :: &
      snlay                     ! number of snow layers

   real(r8), intent(inout)  :: &
      dzlak(nlake)            ,&! lake layer thickness [m]
      zlake(nlake)            ,&! lake layer node depth [m]
      zilak(nlake+1)          ,&! lake layer interface level [m]
      lktmp(nlake)            ,&! lake temperature (K)
      icefr(nlake)              ! lake ice fraction [-]

   real(r8), intent(inout)  :: &
      dzssb (-nsnow+1:nsoil)  ,&! snow + soil + bedrock layer thickness [m]
      zssb  (-nsnow+1:nsoil)  ,&! snow + soil + bedrock node layer depth [m]
      zissb (-nsnow:nsoil+1)  ,&! snow + soil + bedrock layer interface level [m]
      t_ssb (-nsnow+1:nsoil)  ,&! snow + soil + bedrock layer temperature [K]
      xwliq (-nsnow+1:nsoil)  ,&! snow + soil + bedrock layer liquid water (kg/m2)
      xwice (-nsnow+1:nsoil)    ! snow + soil + bedrock layer ice lens (kg/m2)
        
   real(r8), intent(inout)  :: &
      tskin                   ,&! ground surface temperature [k]
      scv                     ,&! snow mass (kg/m2)
      snwag                   ,&! non dimensional snow age [-]
      stke1                   ,&! top level eddy conductivity [W/m/K]
      z0m                     ,&! roughness length over ground, momentum [m]
      z0h                     ,&! roughness length over ground, sensible heat [m]
      z0q                     ,&! roughness length over ground, latent heat [m]
      felak                   ,&! lake fetch (>0), otherwise USE emperical eq by dplak [m]
      gamma                   ,&! Mixing enhancement factorfor eddy diffusion coefficient
      etal                    ,&! extinction coefficient [1/m]
      btpri                   ,&! beta prime in Monin-Obukhov theory
      frlak                   ,&! lake fraction [-]
      tmsno                   ,&! Temperature at the air-snow interface [K] 
      tmice                   ,&! Temperature at the snow-ice or air-ice interface [K]
      tmmnw                   ,&! Mean temperature of the water column [K]
      tmwml                   ,&! Mixed-layer temperature [K]
      tmbot                   ,&! Temperature at the water-bottom sediment interface [K]
      tmups                   ,&! Temperature at the bottom of the upper layer of the sediments [K]
      CTfrac                  ,&! Shape factor (thermocline)
      snwdp                   ,&! snow depth [m]
      icedp                   ,&! ice depth [m]
      mldp                    ,&! mixed layer depth [m]
      upsdp                     ! bottom of the upper layer of the sediments [m]

!  ------------------------- output variables ---------------------------

   real(r8), intent(out)    :: &
      fsena                   ,&! sensible heat from canopy height to atmosphere [W/m2]
      fevpa                   ,&! evapotranspiration from canopy height to atmosphere [mm/s]
      lfevpa                  ,&! latent heat flux from canopy height to atmosphere [W/2]
      fseng                   ,&! sensible heat flux from ground [W/m2]
      fevpg                   ,&! evaporation heat flux from ground [mm/s]
      olrg                    ,&! outgoing long-wave radiation from ground+canopy
      fgrnd                   ,&! ground heat flux [W/m2]
      trad                    ,&! radiative temperature [K]
      qseva                   ,&! ground surface evaporation rate [mm h2o/s]
      qsdew                   ,&! ground surface dew formation [mm h2o/s] (+)
      qsubl                   ,&! sublimation rate from snow pack [mm h2o/s] (+)
      qfros                   ,&! surface dew added to snow pack [mm h2o/s] (+)   
      taux                    ,&! wind stress: E-W [kg/m/s**2]
      tauy                    ,&! wind stress: N-S [kg/m/s**2]
      ustar                   ,&! u* in similarity theory [m/s]
      qstar                   ,&! q* in similarity theory [kg/kg]
      tstar                   ,&! t* in similarity theory [K]
      emis                    ,&! averaged bulk surface emissivity
      zol                     ,&! dimensionless height (z/L) used in Monin-Obukhov theory
      snwml                   ,&! snowmelt rate [mm/s]
      rib                     ,&! bulk Richardson number in surface layer
      ram                     ,&! aerodynamical resistance [s/m]
      rah                     ,&! thermal resistance [s/m]
      raw                     ,&! moisture resistance [s/m]
      wdm                     ,&! lowest level wind speed including the stablity effect [m/s]
      t2m                     ,&! 2 m height air temperature [K]
      q2m                     ,&! 2 m height air specific humidity
      u10m                    ,&! 10 m height wind speed in eastward direction [m/s]
      v10m                    ,&! 10 m height wind speed in northward direction [m/s]
      fm10m                   ,&! integral of profile function for momentum
      fq10m                   ,&! integral of profile function for moisture
      fh2m                    ,&! integral of profile function for heat
      fq2m                    ,&! integral of profile function for moisture
      fm                      ,&! integral of profile function for momentum
      fq                      ,&! integral of profile function for moisture
      fh                      ,&! integral of profile function for heat
      shfdt                     ! derivative of srf sen+lat  heat flux wrt srf temp [W/m2/K]

   logical, optional, intent(in) :: urban_call   ! whether it is a urban CALL

!  ------------------------- local variables ---------------------------
   integer                  :: &
      i                       ,&! loop index
      lb                        ! lower bound for snow layer

   real(r8)                 :: &
      rlat                    ,&! latitude in radians
      par_Coriolis            ,&! The Coriolis parameter [s^{-1}]
      scvold                  ,&! old snow mass (kg/m2)
      iceml                   ,&! ice mass [mm/s]
      phi                     ,&! --
      fgrnd1                  ,&! ground heat flux into the first snow/lake layer [W/m2]
      Q_w_flk                 ,&! Heat flux through the ice-water or air-water interface [W m^{-2}]
      Q_snow_flk              ,&! Heat flux through the air-snow interface [W m^{-2}]
      Q_ice_flk               ,&! Heat flux through the snow-ice or air-ice interface [W m^{-2}]
      u_star_w_flk            ,&! Friction velocity [m s^{-1}]
      mmtum                   ,&! Momentum flux [N m^{-2}]
      Q_sensible              ,&! Sensible heat flux [W m^{-2}]
      Q_latent                ,&! Latent heat flux [W m^{-2}]
      Q_watvap                ,&! Flux of water vapour [kg m^{-2} s^{-1}]
      I_snow_flk              ,&! Radiation flux through the air-snow interface [W m^{-2}]
      I_ice_flk               ,&! Radiation flux through the snow-ice or air-ice interface [W m^{-2}]
      I_w_flk                 ,&! Radiation flux through the ice-water or air-water interface [W m^{-2}]
      I_h_flk                 ,&! Radiation flux through the mixed-layer-thermocline interface [W m^{-2}]
      I_bot_flk               ,&! Radiation flux through the water-bottom sediment interface [W m^{-2}]
      I_intm_0_h_flk          ,&! Mean radiation flux over the mixed layer [W m^{-1}]
      I_intm_h_D_flk          ,&! Mean radiation flux over the thermocline [W m^{-1}]
      pg_snow                 ,&! snowfall onto ground including canopy runoff [kg/(m2 s)]
      zeta                    ,&! dimensionless height (z/L) used in Monin-Obukhov theory
      ztmp                    ,&! temporary variable
      htvp                    ,&! heat transfer coefficient for vapor [W/m2/K]
      ewatdm                  ,&! Equivalent to the water depth after melting.
      watwei                  ,&! water weight [kg]
      rhosnow                 ,&! snow density [kg/m3]
      u2m                     ,&! wind speed at 2m height [m/s]
      U_a                     ,&! wind speed at reference height [m/s]
      wice_lake(nlake)        ,&! ice lens [kg/m2]
      wliq_lake(nlake)          ! liquid water [kg/m2]
!================================================================================

!------------------------------------------------------------------------------
!  Set the initial values
!------------------------------------------------------------------------------
      scvold = scv
      U_a = max(0.1,sqrt(usurf*usurf+vsurf*vsurf))     ! limit set to 0.1
      par_Coriolis= 2*omega*sin(xlat)                  ! Coriolis parameter
      rlat = xlat * d2r
      pg_snow = csnow + lsnow
        
!------------------------------------------------------------------------------
!  Compute solar radiation fluxes (positive downward) at the snow-ice, ice-water, air-water,
!------------------------------------------------------------------------------
      CALL flake_radtrans ( &
            ! "in" arguments
            ! ---------------------------
            dplak         , icedp         , snwdp    ,&
            mldp          , sabg          , etaopt   ,&
            ! "inout" arguments
            ! ---------------------------
            etal          ,&
            ! "out" arguments
            ! ---------------------------
            I_snow_flk    , I_ice_flk     , I_w_flk  ,&
            I_h_flk       , I_bot_flk     ,&
            I_intm_0_h_flk, I_intm_h_D_flk)

!------------------------------------------------------------------------------
!  Compute lake surface fluxes
!------------------------------------------------------------------------------
      lb = snlay + 1
      IF (DEF_USE_COLML_FLUX_SCHEME) THEN

         CALL CoLML_SurfFlux ( &
               ! "in" arguments
               ! -------------------
               scwat       , snlay       , zopt         , betaopt     ,&
               fetchopt    , rlat        , dplak        , dtlak       ,& 
               stke1       , hwref       , htref        , hqref       ,&
               usurf       , vsurf       , tmref        , qmref       ,&
               arhos       , psurf       , sols         , soll        ,&
               solsd       , solld       , lwdns        , sabg        ,&
               hpbl        , dzlak(1)    , zlake(nlake) , dzssb(lb)   ,&
               lktmp(1)    , t_ssb(lb)   , xwliq(lb)    , xwice(lb)   ,&
               icefr(1)    , zcbcv       , ipatch       ,&
               ! "inout" arguments
               ! -------------------
               tskin       , z0m         , z0h          , z0q         ,&
               felak       , btpri       ,&
               ! "out" arguments
               ! -------------------
               fseng       , fevpg       , fsena        , fevpa       ,&
               lfevpa      , olrg        , fgrnd        , fgrnd1      ,&
               zol         , rib         , trad         , htvp        ,&
               emis        , wdm         , ram          , rah         ,&
               raw         , shfdt       , taux         , tauy        ,&
               t2m         , q2m         , u10m         , v10m        ,&
               fh2m        , fq2m        , fm10m        , fq10m       ,&
               fm          , fh          , fq           , ustar       ,&
               qstar       , tstar       , rhosnow      , u2m         )
         u_star_w_flk = ustar**2
         Q_w_flk = lwdns - olrg - fsena - lfevpa  

      ELSE

         CALL FLake_SurfFlux ( &
               ! "in" arguments
               ! ---------------------------
               scwat       , zopt        , fetchopt     , dplak       ,&
               tmref       , qmref       , arhos        , U_a         ,&
               icedp       , hwref       , hqref        , tskin       ,&
               lwdns       , sabg        , psurf        ,&
               ! "inout" arguments
               ! -------------------
               z0m         , z0h         , z0q          , felak       ,&
               ! "out" arguments
               ! ---------------------------
               mmtum       , fsena       , lfevpa       , fevpa       ,&
               olrg        , fgrnd       )

         CALL LakStaParms ( &
               ! "in" arguments
               ! -------------------
               hwref       , htref       , hqref        , usurf       ,&
               vsurf       , tmref       , qmref        , tskin       ,&
               arhos       , psurf       , zcbcv        ,& 
               ! "out" arguments
               ! -------------------
               zol         , rib         , wdm          , shfdt       ,&
               taux        , tauy        , t2m          , q2m         ,&
               u10m        , v10m        , fh2m         , fq2m        ,&
               fm10m       , fq10m       , fm           , fh          ,&
               fq          , ustar       , qstar        , tstar       ,&
               ram         , rah         , raw          , emis        )
         fseng = fsena
         fevpg = fevpa
         u_star_w_flk = SQRT(-mmtum/tpl_rho_w_r)  ! Friction velocity [m s^{-1}]
         Q_w_flk = lwdns - olrg - fsena - lfevpa
         trad = (olrg/stefnc)**0.25

      ENDIF

!------------------------------------------------------------------------------
!  Compute heat fluxes Q_snow_flk, Q_ice_flk, Q_w_flk
!------------------------------------------------------------------------------
      IF(icedp.GE.h_Ice_min_flk) THEN            ! Ice exists
         IF(snwdp.GE.h_Snow_min_flk) THEN       ! There is snow above the ice
               Q_snow_flk = Q_w_flk
               Q_ice_flk  = 0.
               Q_w_flk    = 0.
         ELSE                                   ! No snow above the ice
               Q_snow_flk = 0.
               Q_ice_flk  = Q_w_flk
               Q_w_flk    = 0.
         ENDIF
      ELSE                                       ! No ice-snow cover
         Q_snow_flk = 0.
         Q_ice_flk  = 0.
      ENDIF

!------------------------------------------------------------------------------
!  Advance FLake variables
!------------------------------------------------------------------------------
      CALL flaketem ( &
            ! "in" arguments
            ! ---------------------------
            dplak         , dpsed           , par_Coriolis   ,& 
            etal          , dtlak           , I_w_flk        ,& 
            I_h_flk       , I_intm_0_h_flk  , I_bot_flk      ,&
            Q_snow_flk    , Q_ice_flk       , I_snow_flk     ,&
            I_ice_flk     , pg_snow         , T_bs           ,&
            u_star_w_flk  , I_intm_h_D_flk  ,&
            ! "inout" arguments
            ! -------------------
            Q_w_flk       ,&
            tskin         , tmsno           , tmice          ,&
            tmwml         , tmmnw           , tmbot          ,&
            tmups         , snwdp           , icedp          ,&
            mldp          , upsdp           , CTfrac         ,&
            ! "out" arguments
            ! -------------------
            snwml         , iceml           )

!------------------------------------------------------------------------------
!  Calculate the lake temperature for each layer (Mironov 2005, eq.56)
!------------------------------------------------------------------------------
      DO i = 1, nlake
         IF (zilak(i)<= mldp) THEN
            lktmp(i) = tmwml
         ELSE
            zeta = (zlake(i) - mldp)/(dplak + 0.01 - mldp)
            phi  = (40./3. * CTfrac - 20./3.) * zeta + &
                  (18. - 30. * CTfrac) * zeta**2 +    &
                  (20. * CTfrac - 12.) * zeta**3 +    &
                  (5./3. - 10./3. * CTfrac) * zeta**4
            lktmp(i) = tmwml - (tmwml - tmbot)*phi
         ENDIF
      ENDDO

!------------------------------------------------------------------------------
!  update the snow layer
!------------------------------------------------------------------------------
      CALL snowdp2lev4lake( &
            ! "in" arguments
            ! ---------------------------
            nsnow     , snwdp    ,&
            ! "inout" arguments
            ! ---------------------------
            dzssb(:0) , zssb(:0) , zissb(:0), snlay )

!------------------------------------------------------------------------------
!  Setting snow related variables
!------------------------------------------------------------------------------
      ! calculate sublimation, frosting, dewing, following the CoLM-Lake model
      qseva = 0.
      qsubl = 0.
      qsdew = 0.
      qfros = 0.
      IF (fevpg >= 0.0) THEN
         IF(lb < 0)THEN
            qseva = min(xwliq(lb)/dtlak, fevpg)
            qsubl = fevpg - qseva
         ELSE
            qseva = min((1.-icefr(1))*1000.*dzlak(1)/dtlak, fevpg)
            qsubl = fevpg - qseva
         ENDIF
      ELSE
         IF (tskin < tfrz) THEN
            qfros = abs(fevpg)
         ELSE
            qsdew = abs(fevpg)
         ENDIF
      ENDIF
      IF (snwdp > 0.) scv = scv + (pg_snow-snwml-qsubl+qfros)*dtlak
      scv = max( scv, 0. )
      
      ! no snow IF lake unfrozen
      IF (tskin > tfrz) scv = 0.
      CALL snowage (dtlak,tskin,scv,scvold,snwag)

!------------------------------------------------------------------------------
!  Setting snow + soil temperature
!------------------------------------------------------------------------------
      xwliq = 0.
      xwice = 0.
      t_ssb(snlay+1:0) = tmsno
      DO i = 1, nsoil
         IF (upsdp > zssb(i)) THEN
            t_ssb(i) = tmups
         ELSE
            t_ssb(i) = T_bs
         ENDIF
      ENDDO

!------------------------------------------------------------------------------
!  Setting ice related variables
!------------------------------------------------------------------------------
      !+WMEJ: Here it is assumed that the thickness conversion ratio between lake water and ice is 1:1,
      !       which is incorrect. The water balance of the lake has not been taken into account.
      !+WMEJ TODO: Consider the water balance of the lake
      IF (icedp < h_Ice_min_flk) THEN 
         DO i = 1, nlake
            icefr(i) = 0.
            wice_lake(i) = 0.
            wliq_lake(i) = tpl_rho_w_r*dzlak(i)
         ENDDO        
      ELSE
         DO i = 1, nlake
            IF (icedp > zilak(i+1)) THEN
               icefr(i) = 1.
               wice_lake(i) = tpl_rho_I*dzlak(i)
               wliq_lake(i) = tpl_rho_w_r*dzlak(i)
            ELSE 
               wice_lake(i) = tpl_rho_I*(icedp - zilak(i))
               wliq_lake(i) = tpl_rho_w_r*dzlak(i)
               icefr(i) = wice_lake(i)/(wice_lake(i) + wliq_lake(i))
            ENDIF
         ENDDO
      ENDIF

   END SUBROUTINE Lake_FLake



   SUBROUTINE flake_radtrans ( &
               ! "in" arguments
               ! ---------------------------
               depth_w       , h_ice_p_flk   , h_snow_p_flk  , h_ML_p_flk    ,&
               I_atm_flk     , etaopt        ,&
               ! "inout" arguments
               ! ---------------------------
               etal          ,&
               ! "out" arguments
               ! ---------------------------
               I_snow_flk    , I_ice_flk     , I_w_flk       , I_h_flk       ,&
               I_bot_flk     , I_intm_0_h_flk, I_intm_h_D_flk )

! ---------------------------------- code history -------------------------------------
! Description:
!
!  Computes the radiation fluxes 
!  at the snow-ice, ice-water, air-water, 
!  mixed layer-thermocline and water column-bottom sediment interfaces,
!  the mean radiation flux over the mixed layer,
!  and the mean radiation flux over the thermocline.
!
! ----------------------------------------------------------
! +WMEJ Added for a more visible presentation process
! ----------------------------------------------------------
!
!    | solor-radiation         ^                         
!    v                         | reflected radiation  ----> I_snow_flk    
!   ~~~~~~~~~~      snow       ~~~~~~~~~~
!   ------------------------------------------------------> I_ice_flk
!   ==========       ice       ==========
!   ------------------------------------------------------> I_w_flk
!   ++++++++++   mixed-layer   ++++++++++                  -------------> I_intm_0_h_flk
!   ------------------------------------------------------> I_h_flk   
!   //////////   thermocline   //////////                  -------------> I_intm_h_D_flk
!   ------------------------------------------------------> I_bot_flk
!   ********** bottom-sediment **********
!
! Original author: 
!     Dmitrii Mironov, 2005/11/17
!
! Revisions:
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024: Added single band radiation scheme
!
! -------------------------------------------------------------------------------------
   USE MOD_Lake_Const, ONLY: h_Snow_min_flk, h_Ice_min_flk, h_ML_min_flk, opticpar_medium,&
                              opticpar_water=>opticpar_water_ref, opticpar_ice=>opticpar_ice_opaque,&
                              opticpar_snow=>opticpar_snow_opaque, stefnc, c_lwrad_emis
!==============================================================================
!  ------------------------- inout variables ---------------------------
   integer, intent(in)      :: &
      etaopt                    ! option for Extinction coefficient calculation, 1: constant, 2: equation

   real(r8), intent(in)     :: &
      depth_w                 ,&! The lake depth [m]
      h_ice_p_flk             ,&! Ice depth [m]
      h_snow_p_flk            ,&! Snow depth [m]
      h_ML_p_flk              ,&! Mixed-layer depth [m]
      I_atm_flk                 ! Atmospheric radiation flux [W/m2]

!  ------------------------- inout variables ---------------------------
   real(r8), intent(inout)  :: &
      etal                      ! extinction coefficient [1/m]

!  ------------------------- output variables ---------------------------
   real(r8), intent(out)    :: &
      I_snow_flk              ,&! Radiation flux through the air-snow interface [W m^{-2}]
      I_ice_flk               ,&! Radiation flux through the snow-ice or air-ice interface [W m^{-2}]
      I_w_flk                 ,&! Radiation flux through the ice-water or air-water interface [W m^{-2}]
      I_h_flk                 ,&! Radiation flux through the mixed-layer-thermocline interface [W m^{-2}]
      I_bot_flk               ,&! Radiation flux through the water-bottom sediment interface [W m^{-2}]
      I_intm_0_h_flk          ,&! Mean radiation flux over the mixed layer [W m^{-1}]
      I_intm_h_D_flk            ! Mean radiation flux over the thermocline [W m^{-1}]

!  ------------------------- local variables ---------------------------
   integer  :: i                ! DO loop index
   !+WEMJ Explicitly add intermediate variables to solve the original variable(I_bot_flk) reuse problem
   real(r8) :: radext           ! Intermediate variables in radiation calculations
!================================================================================

   !==============================================================================
   !  Start calculations
   !------------------------------------------------------------------------------
      !+WMEJ Radiation reflection to be done at initialization
      !+WMEJ <I_snow_flk>, <I_ice_flk>, <I_w_flk>
      IF(h_ice_p_flk.GE.h_Ice_min_flk) THEN            ! Ice exists
         IF(h_snow_p_flk.GE.h_Snow_min_flk) THEN      ! There is snow above the ice
            I_snow_flk = I_atm_flk !-WMEJ *(1.-albedo_snow)  !+WMEJ radiation flux through the air-snow interface
            radext = 0.                           
            DO i=1, opticpar_snow%nband_optic
               radext = radext +                    & 
               opticpar_snow%frac_optic(i)*EXP(-opticpar_snow%extincoef_optic(i)*h_snow_p_flk) 
            ENDDO   !+WMEJ It is not recommended to USE irrelevant variables to replace intermediate variables.
            I_ice_flk  = I_snow_flk*radext        
         ELSE                                           ! No snow above the ice 
            I_snow_flk = I_atm_flk  
            I_ice_flk  = I_atm_flk !-WMEJ *(1.-albedo_ice)
         ENDIF 
         radext = 0.
         DO i=1, opticpar_ice%nband_optic
            radext = radext +                      & 
            opticpar_ice%frac_optic(i)*EXP(-opticpar_ice%extincoef_optic(i)*h_ice_p_flk) 
         ENDDO 
         I_w_flk      = I_ice_flk*radext
      ELSE                                             ! No ice-snow cover
         I_snow_flk   = I_atm_flk  
         I_ice_flk    = I_atm_flk
         I_w_flk      = I_atm_flk !-WMEJ *(1.-albedo_water)
      ENDIF 

      !+WMEJ <I_h_flk>
      IF(h_ML_p_flk.GE.h_ML_min_flk) THEN           ! Radiation flux at the bottom of the mixed layer
         IF (etaopt == 2) THEN
            radext = 0.
            DO i=1, opticpar_water%nband_optic
               radext = radext +            & 
               opticpar_water%frac_optic(i)*EXP(-opticpar_water%extincoef_optic(i)*h_ML_p_flk) 
            ENDDO
            I_h_flk = I_w_flk*radext
         ELSE
            I_h_flk = I_w_flk*EXP(-etal*h_ML_p_flk)
         ENDIF
      ELSE                                          ! Mixed-layer depth is less THEN a minimum value
         I_h_flk = I_w_flk
      ENDIF

      !+WMEJ <I_bot_flk>
      IF (etaopt == 2) THEN
         radext = 0.                         ! Radiation flux at the lake bottom
         DO i=1, opticpar_water%nband_optic
            radext = radext +              & 
            opticpar_water%frac_optic(i)*EXP(-opticpar_water%extincoef_optic(i)*depth_w) 
         ENDDO 
         I_bot_flk = I_w_flk*radext
      ELSE
         I_bot_flk = I_w_flk*EXP(-etal*depth_w)
      ENDIF

      !+WMEJ <I_intm_0_h_flk>
      IF(h_ML_p_flk.GE.h_ML_min_flk) THEN   ! Integral-mean radiation flux over the mixed layer
         IF (etaopt == 2) THEN
            radext = 0.
            DO i=1, opticpar_water%nband_optic
               radext = radext +                                &
               opticpar_water%frac_optic(i)/opticpar_water%extincoef_optic(i)*  &
               (1. - EXP(-opticpar_water%extincoef_optic(i)*h_ML_p_flk))
            ENDDO 
            I_intm_0_h_flk = I_w_flk*radext/h_ML_p_flk
         ELSE
            I_intm_0_h_flk = I_w_flk*EXP(-etal*h_ML_p_flk)/h_ML_p_flk
         ENDIF
      ELSE
         I_intm_0_h_flk = I_h_flk
      ENDIF

      !+WMEJ <I_intm_h_D_flk>
      IF(h_ML_p_flk.LE.depth_w-h_ML_min_flk) THEN   ! Integral-mean radiation flux over the thermocline
         IF (etaopt == 1) THEN
            radext = 0. 
            DO i=1, opticpar_water%nband_optic
               radext = radext +                                &
               opticpar_water%frac_optic(i)/opticpar_water%extincoef_optic(i)*  &
               ( EXP(-opticpar_water%extincoef_optic(i)*h_ML_p_flk)             &
               - EXP(-opticpar_water%extincoef_optic(i)*depth_w) )
            ENDDO 
            I_intm_h_D_flk = I_w_flk*radext/(depth_w-h_ML_p_flk)
         ELSE
            I_intm_h_D_flk = I_w_flk*(EXP(-etal*h_ML_p_flk)-EXP(-etal*depth_w))/(depth_w-h_ML_p_flk)
         ENDIF
      ELSE
         I_intm_h_D_flk = I_h_flk
      ENDIF
      
      !------------------------------------------------------------------------------
      !  END calculations
      !==============================================================================

   END SUBROUTINE flake_radtrans



   SUBROUTINE FLake_SurfFlux ( &
               ! "in" arguments
               ! ---------------------------
               scwat        , z_opt        , fetchopt    , dplak      ,&
               T_a          , q_a          , arhos       , U_a        ,&
               h_ice        , height_u     , height_tq   , T_s        ,&
               lwdns        , sabg         , P_a         ,&
               ! "inout" arguments
               ! -------------------
               z0m          , z0h          , z0q         , felak      ,&
               ! "out" arguments
               ! ---------------------------
               Q_momentum   , Q_sensible   , Q_latent    , Q_watvap   ,&
               Q_olrg       , Q_fgrnd      )

! ---------------------------------- code history -------------------------------------
! Description:
!
!     The SfcFlx routine 
!     where fluxes of momentum and of sensible and latent heat 
!     at the air-water or air-ice (air-snow) interface are computed. 
!
!     Lines embraced with "!_tmp" contain temporary parts of the code.
!     Lines embraced/marked with "!_dev" may be replaced
!     as improved parameterizations are developed and tested.
!     Lines embraced/marked with "!_dm" are DM's comments
!     that may be helpful to a user.
!     Lines embraced/marked with "!_dbg" are used 
!     for debugging purposes only.
!
! Original author: 
!     Dmitrii Mironov, 2005/11/17
!
!  Omarjan Obulkasim (WMEJ), March, 2024 : adjusted the in/out variable to make 
!                                          it can be used in the CoLM model
!==============================================================================
   USE MOD_Lake_Const,only: tpsf_nu_u_a, tpsf_kappa_t_a, tpsf_kappa_q_a, &
                                 tpsf_Rd_o_Rv, grav, tpsf_nu_u_a, c_free_conv, &
                                 c_MO_t_stab, c_MO_u_stab, tpsf_alpha_q, u_wind_min_sf, &
                                 Pr_neutral, c_small_sf, c_accur_sf, u_star_min_sf, &
                                 c_MO_u_exp, vonkar, h_Ice_min_flk, c_MO_q_conv,&
                                 c_MO_q_exp, c_MO_t_exp, tpsf_c_a_p, tpsf_L_evap ,&
                                 tpl_L_f, c_MO_t_conv, c_MO_q_stab, Sc_neutral,&
                                 c_lwrad_emis, stefnc, c_MO_u_conv, SHALLOW, DEEP
!==============================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in)      :: &
      scwat                   ,&! land patch type (1=shallow lake, 2=deep lake)
      fetchopt                ,&! option for fetch length, 1: constant, 2: fetch length is calculated
      z_opt                     ! option for roughness length, 1: constant, 2: Subin et al. (2012)

   real(r8), intent(in)     :: &
      dplak                   ,&! lake depth [m]
      T_a                     ,&! Air temperature [K]
      q_a                     ,&! Air specific humidity [-]
      arhos                   ,&! Air density [kg m^{-3}]
      U_a                     ,&! Wind speed [m s^{-1}]
      P_a                     ,&! Air pressure [Pa]
      T_s                     ,&! Surface temperature (water, ice or snow) [K]
      sabg                    ,&! Solar absorbed by ground [W m^{-2}]
      lwdns                   ,&! Atmospheric infrared (longwave) radiation [W m^{-2}]
      h_ice                   ,&! Ice thickness [m]
      height_u                ,&! Height where wind is measured [m]
      height_tq                 ! Height where temperature and humidity are measured [m]

!  ------------------------- inout variables ---------------------------
   real(r8), intent(inout)  :: &
      z0m                     ,&! roughness length over ground, momentum [m]
      z0h                     ,&! roughness length over ground, sensible heat [m]
      z0q                     ,&! roughness length over ground, latent heat [m]
      felak                     ! fetch length of lake [m]

!  ------------------------- output variables ---------------------------
   real(r8), intent(out)    :: &
      Q_olrg                  ,&! outgoing longwave radiation [W/m2]
      Q_fgrnd                 ,&! ground heat flux [W/m2]
      Q_momentum              ,&! Momentum flux [N m^{-2}]  
      Q_sensible              ,&! Sensible heat flux [W m^{-2}]  
      Q_latent                ,&! Laten heat flux [W m^{-2}]
      Q_watvap                  ! Flux of water vapout [kg m^{-2} s^{-1}]

!  ------------------------- local variables ---------------------------
   real(r8)                 :: &
      wvpres_s                ,&! Saturation water vapour pressure at T=T_s [N m^{-2}]
      rho_a                   ,&! Air density [kg m^{-3}]
      q_s                     ,&! Saturation specific humidity at T=T_s [-]
      Q_mom_mol               ,&! Molecular momentum flux [N m^{-2}]
      Q_sen_mol               ,&! Molecular sensible heat flux [W m^{-2}]  
      Q_lat_mol               ,&! Molecular laten heat flux [W m^{-2}]
      Q_mom_con               ,&! Momentum flux in free convection [N m^{-2}]
      Q_sen_con               ,&! Sensible heat flux in free convection [W m^{-2}]  
      Q_lat_con               ,&! Laten heat flux in free convection [W m^{-2}]
      Q_mom_tur               ,&! Turbulent momentum flux [N m^{-2}]
      Q_sen_tur               ,&! Turbulent sensible heat flux [W m^{-2}]  
      Q_lat_tur               ,&! Turbulent laten heat flux [W m^{-2}]
      R_z                     ,&! Ratio of "height_tq" to "height_u"
      ZoL                     ,&! The z/L ratio, z=height_u
      Ri                      ,&! Gradient Richardson number 
      Ri_cr                   ,&! Critical value of Ri 
      Fun                     ,&! A function of generic variable "x"
      Delta                   ,&! Relative error 
      Fun_prime               ,&! Derivative of "Fun" with respect to "x"
      u_star_st               ,&! Friction velocity with due regard for stratification [m s^{-1}]
      c_z0u_fetch             ,&! Fetch-dependent Charnock parameter
      par_conv_visc           ,&! Viscous convection stability parameter
      u_star_previter         ,&! Friction velocity from previous iteration [m s^{-1}]
      u_star_thresh           ,&! Threshld value of friction velocity [m s^{-1}]
      z0u_sf                  ,&! Roughness length for momentum [m]
      z0t_sf                  ,&! Roughness length for temperature [m]
      z0q_sf                  ,&! Roughness length for specific humidity [m]
      psi_u                   ,&! The MO stability function for wind profile
      psi_q                   ,&! The MO stability function for specific humidity profile
      psi_t                   ,&! The MO stability function for temperature profile
      fetch                   ,&! Fetch [m]
      U_a_thresh                ! Threshld value of the wind speed [m s^{-1}] 

   logical  :: l_conv_visc      ! Switch, TRUE = viscous free convection, the Nu=C Ra^(1/3) law is used
   integer  :: i, n_iter        ! Loop index, number of iterations performed
   integer, PARAMETER :: n_iter_max = 24  ! Maximum number of iterations
!================================================================================
      ! Base on lake depth, assuming that small lakes are likely to be shallower
      ! Estimate crudely based on lake depth
      !-WMEJ USE scwat to control lake classification, 
      !      1 = shallow lake, 2 = deep lake, 2024-03-23
      IF (fetchopt == 1) THEN
         fetch = felak
      ELSE
         IF (dplak < 4.) THEN
            fetch = 100. ! shallow lake
         ELSE
            fetch = 25.*dplak ! deep lake
         ENDIF
      ENDIF

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------
      !_dm All fluxes are positive when directed upwards.
!------------------------------------------------------------------------------
!  Compute saturation specific humidity and the air density at T=T_s
!------------------------------------------------------------------------------
      wvpres_s = SfcFlx_satwvpres(T_s, h_ice)  ! Saturation water vapour pressure at T=T_s
      q_s = SfcFlx_spechum (wvpres_s, P_a)     ! Saturation specific humidity at T=T_s
      !+WMEJ: is it necessary to calculate the air density(rho_a) at T=T_s in this SUBROUTINE?
      !       how about provided by the forcing data?
      ! rho_a = arhos  !+WMEJ provided by the forcing data
      rho_a = SfcFlx_rhoair(T_s, q_s, P_a)     ! Air density at T_s and q_s (surface values)
!------------------------------------------------------------------------------
!  Compute molecular fluxes of momentum and of sensible and latent heat
!------------------------------------------------------------------------------
      !_dm The fluxes are in kinematic units
      Q_mom_mol = -tpsf_nu_u_a*U_a/height_u 
      Q_sen_mol = -tpsf_kappa_t_a*(T_a-T_s)/height_tq    
      Q_lat_mol = -tpsf_kappa_q_a*(q_a-q_s)/height_tq  
      
!------------------------------------------------------------------------------
!  Compute fluxes in free convection
!------------------------------------------------------------------------------
      par_conv_visc = (T_s-T_a)/T_s*SQRT(tpsf_kappa_t_a) + (q_s-q_a)*tpsf_alpha_q*SQRT(tpsf_kappa_q_a)
      IF(par_conv_visc.GT.0.) THEN   ! Viscous convection takes place
         l_conv_visc = .TRUE.
         par_conv_visc = (par_conv_visc*grav/tpsf_nu_u_a)**(1./3.)
         Q_sen_con = c_free_conv*SQRT(tpsf_kappa_t_a)*par_conv_visc  
         Q_sen_con = Q_sen_con*(T_s-T_a)
         Q_lat_con = c_free_conv*SQRT(tpsf_kappa_q_a)*par_conv_visc
         Q_lat_con = Q_lat_con*(q_s-q_a)
      ELSE                                  ! No viscous convection, set fluxes to zero
         l_conv_visc = .FALSE.
         Q_sen_con = 0. 
         Q_lat_con = 0.
      ENDIF
      Q_mom_con = 0.                 ! Momentum flux in free (viscous or CBL-scale) convection is zero  

!------------------------------------------------------------------------------
!  Compute turbulent fluxes
!------------------------------------------------------------------------------
      R_z   = height_tq/height_u                        ! Ratio of "height_tq" to "height_u"
      Ri_cr = c_MO_t_stab/c_MO_u_stab**2*R_z  ! Critical Ri
      Ri    = grav*((T_a-T_s)/T_s+tpsf_alpha_q*(q_a-q_s))/MAX(U_a,u_wind_min_sf)**2
      Ri    = Ri*height_u/Pr_neutral                    ! Gradient Richardson number

      Turb_Fluxes: IF(U_a .LT. u_wind_min_sf .OR. Ri .GT. Ri_cr-c_small_sf) THEN  ! Low wind or Ri>Ri_cr 
         u_star_st = 0.                       ! Set turbulent fluxes to zero 
         Q_mom_tur = 0.                       
         Q_sen_tur = 0.   
         Q_lat_tur = 0.  

      ELSE Turb_Fluxes                            ! Compute turbulent fluxes using MO similarity

         ! Compute z/L, where z=height_u
         IF(Ri.GE.0.) THEN   ! Stable stratification
            ZoL = SQRT(1.-4.*(c_MO_u_stab-R_z*c_MO_t_stab)*Ri)
            ZoL = ZoL - 1. + 2.*c_MO_u_stab*Ri
            ZoL = ZoL/2./c_MO_u_stab/c_MO_u_stab/(Ri_cr-Ri)
         ELSE                       ! Convection
            n_iter = 0
            Delta = 1.                ! Set initial error to a large value (as compared to the accuracy)
            u_star_previter = Ri*MAX(1., SQRT(R_z*c_MO_t_conv/c_MO_u_conv)) ! Initial guess for ZoL
            DO WHILE (Delta.GT.c_accur_sf.AND.n_iter.LT.n_iter_max) 
               Fun = u_star_previter**2*(c_MO_u_conv*u_star_previter-1.)  &
                  + Ri**2*(1.-R_z*c_MO_t_conv*u_star_previter)
               Fun_prime = 3.*c_MO_u_conv*u_star_previter**2              &
                        - 2.*u_star_previter - R_z*c_MO_t_conv*Ri**2
               ZoL = u_star_previter - Fun/Fun_prime
               Delta = ABS(ZoL-u_star_previter)/MAX(c_accur_sf, ABS(ZoL+u_star_previter))
               u_star_previter = ZoL
               n_iter = n_iter + 1
            ENDDO 
         ENDIF
            
         !  Compute fetch-dependent Charnock parameter, USE "u_star_min_sf"
         IF (z_opt == 1) THEN
            CALL SfcFlx_roughness (fetch, U_a, u_star_min_sf, h_ice, c_z0u_fetch, u_star_thresh, z0u_sf, z0t_sf, z0q_sf)
            z0u_sf = z0m
            z0t_sf = z0h
            z0q_sf = z0q
         ELSE
            CALL SfcFlx_roughness (fetch, U_a, u_star_st, h_ice, c_z0u_fetch, u_star_thresh, z0u_sf, z0t_sf, z0q_sf)
         ENDIF 

         !  Threshold value of wind speed 
         u_star_st = u_star_thresh
         IF (z_opt == 1) THEN
            CALL SfcFlx_roughness (fetch, U_a, u_star_st, h_ice, c_z0u_fetch, u_star_thresh, z0u_sf, z0t_sf, z0q_sf)
            z0u_sf = z0m
            z0t_sf = z0h
            z0q_sf = z0q
         ELSE
            CALL SfcFlx_roughness (fetch, U_a, u_star_st, h_ice, c_z0u_fetch, u_star_thresh, z0u_sf, z0t_sf, z0q_sf)
         ENDIF 
         
         IF(ZoL.GT.0.) THEN   ! MO function in stable stratification 
            psi_u = c_MO_u_stab*ZoL*(1.-MIN(z0u_sf/height_u, 1.))
         ELSE                        ! MO function in convection
            psi_t = (1.-c_MO_u_conv*ZoL)**c_MO_u_exp
            psi_q = (1.-c_MO_u_conv*ZoL*MIN(z0u_sf/height_u, 1.))**c_MO_u_exp
            psi_u = 2.*(ATAN(psi_t)-ATAN(psi_q))                  &
                  + 2.*LOG((1.+psi_q)/(1.+psi_t))   &
                  + LOG((1.+psi_q*psi_q)/(1.+psi_t*psi_t))   
         ENDIF 
         U_a_thresh = u_star_thresh/vonkar*(LOG(height_u/z0u_sf)+psi_u)

         !  Compute friction velocity 
         n_iter = 0
         Delta = 1.                ! Set initial error to a large value (as compared to the accuracy)
         u_star_previter = u_star_thresh  ! Initial guess for friction velocity  
         IF(U_a.LE.U_a_thresh) THEN  ! Smooth surface
            DO WHILE (Delta.GT.c_accur_sf.AND.n_iter.LT.n_iter_max) 
               IF (z_opt == 1) THEN
                  CALL SfcFlx_roughness (fetch, U_a, MIN(u_star_thresh, u_star_previter), h_ice,   &
                                       c_z0u_fetch, u_star_thresh, z0u_sf, z0t_sf, z0q_sf)
                  z0u_sf = z0m
                  z0t_sf = z0h
                  z0q_sf = z0q
               ELSE
                  CALL SfcFlx_roughness (fetch, U_a, MIN(u_star_thresh, u_star_previter), h_ice,   &
                                       c_z0u_fetch, u_star_thresh, z0u_sf, z0t_sf, z0q_sf)
               ENDIF 
               IF(ZoL.GE.0.) THEN  ! Stable stratification
                  psi_u = c_MO_u_stab*ZoL*(1.-MIN(z0u_sf/height_u, 1.))
                  Fun = LOG(height_u/z0u_sf) + psi_u
                  Fun_prime = (Fun + 1. + c_MO_u_stab*ZoL*MIN(z0u_sf/height_u, 1.))/vonkar
                  Fun = Fun*u_star_previter/vonkar - U_a
               ELSE                       ! Convection 
                  psi_t = (1.-c_MO_u_conv*ZoL)**c_MO_u_exp
                  psi_q = (1.-c_MO_u_conv*ZoL*MIN(z0u_sf/height_u, 1.))**c_MO_u_exp
                  psi_u = 2.*(ATAN(psi_t)-ATAN(psi_q))                  &
                           + 2.*LOG((1.+psi_q)/(1.+psi_t))   &
                           + LOG((1.+psi_q*psi_q)/(1.+psi_t*psi_t))   
                  Fun = LOG(height_u/z0u_sf) + psi_u
                  Fun_prime = (Fun + 1./psi_q)/vonkar
                  Fun = Fun*u_star_previter/vonkar - U_a
               ENDIF
               u_star_st = u_star_previter - Fun/Fun_prime
               Delta = ABS((u_star_st-u_star_previter)/(u_star_st+u_star_previter))
               u_star_previter = u_star_st
               n_iter = n_iter + 1
            ENDDO 
         ELSE                        ! Rough surface
            DO WHILE (Delta.GT.c_accur_sf.AND.n_iter.LT.n_iter_max) 
               IF (z_opt == 1) THEN
                  CALL SfcFlx_roughness (fetch, U_a, MAX(u_star_thresh, u_star_previter), h_ice,   &
                                       c_z0u_fetch, u_star_thresh, z0u_sf, z0t_sf, z0q_sf)
                  z0u_sf = z0m
                  z0t_sf = z0h
                  z0q_sf = z0q
               ELSE
                  CALL SfcFlx_roughness (fetch, U_a, MAX(u_star_thresh, u_star_previter), h_ice,   &
                                       c_z0u_fetch, u_star_thresh, z0u_sf, z0t_sf, z0q_sf)
               ENDIF 
               IF(ZoL.GE.0.) THEN  ! Stable stratification
                  psi_u = c_MO_u_stab*ZoL*(1.-MIN(z0u_sf/height_u, 1.))
                  Fun = LOG(height_u/z0u_sf) + psi_u
                  Fun_prime = (Fun - 2. - 2.*c_MO_u_stab*ZoL*MIN(z0u_sf/height_u, 1.))/vonkar
                  Fun = Fun*u_star_previter/vonkar - U_a
               ELSE                       ! Convection 
                  psi_t = (1.-c_MO_u_conv*ZoL)**c_MO_u_exp
                  psi_q = (1.-c_MO_u_conv*ZoL*MIN(z0u_sf/height_u, 1.))**c_MO_u_exp
                  psi_u = 2.*(ATAN(psi_t)-ATAN(psi_q))                  &
                           + 2.*LOG((1.+psi_q)/(1.+psi_t))   &
                           + LOG((1.+psi_q*psi_q)/(1.+psi_t*psi_t))   
                  Fun = LOG(height_u/z0u_sf) + psi_u
                  Fun_prime = (Fun - 2./psi_q)/vonkar
                  Fun = Fun*u_star_previter/vonkar - U_a
               ENDIF
               IF(h_ice.GE.h_Ice_min_flk) THEN   ! No iteration is required for rough flow over ice
                  u_star_st = vonkar*U_a/MAX(c_small_sf, LOG(height_u/z0u_sf)+psi_u)
                  u_star_previter = u_star_st
               ELSE                              ! Iterate in case of open water
                  u_star_st = u_star_previter - Fun/Fun_prime
               ENDIF
               Delta = ABS((u_star_st-u_star_previter)/(u_star_st+u_star_previter))
               u_star_previter = u_star_st
               n_iter = n_iter + 1
            ENDDO 
         ENDIF

         !  Momentum flux
         Q_mom_tur = -u_star_st*u_star_st

         !  Temperature and specific humidity fluxes
         IF (z_opt == 1) THEN
            CALL SfcFlx_roughness (fetch, U_a, u_star_st, h_ice, c_z0u_fetch, u_star_thresh, z0u_sf, z0t_sf, z0q_sf)
            z0u_sf = z0m
            z0t_sf = z0h
            z0q_sf = z0q
         ELSE
            CALL SfcFlx_roughness (fetch, U_a, u_star_st, h_ice, c_z0u_fetch, u_star_thresh, z0u_sf, z0t_sf, z0q_sf)
         ENDIF 
         IF(ZoL.GE.0.) THEN   ! Stable stratification 
            psi_t = c_MO_t_stab*R_z*ZoL*(1.-MIN(z0t_sf/height_tq, 1.))
            psi_q = c_MO_q_stab*R_z*ZoL*(1.-MIN(z0q_sf/height_tq, 1.))
         ELSE                        ! Convection 
            psi_u = (1.-c_MO_t_conv*R_z*ZoL)**c_MO_t_exp
            psi_t = (1.-c_MO_t_conv*R_z*ZoL*MIN(z0t_sf/height_tq, 1.))**c_MO_t_exp
            psi_t = 2.*LOG((1.+psi_t)/(1.+psi_u))
            psi_u = (1.-c_MO_q_conv*R_z*ZoL)**c_MO_q_exp
            psi_q = (1.-c_MO_q_conv*R_z*ZoL*MIN(z0q_sf/height_tq, 1.))**c_MO_q_exp
            psi_q = 2.*LOG((1.+psi_q)/(1.+psi_u))
         ENDIF 
         Q_sen_tur = -(T_a-T_s)*u_star_st*vonkar/Pr_neutral  &
                  / MAX(c_small_sf, LOG(height_tq/z0t_sf)+psi_t)
         Q_lat_tur = -(q_a-q_s)*u_star_st*vonkar/Sc_neutral  &
                  / MAX(c_small_sf, LOG(height_tq/z0q_sf)+psi_q)

      ENDIF Turb_Fluxes

!------------------------------------------------------------------------------
!  Decide between turbulent, molecular, and convective fluxes
!------------------------------------------------------------------------------
      Q_momentum = MIN(Q_mom_tur, Q_mom_mol, Q_mom_con)  ! Momentum flux is negative          
      IF(l_conv_visc) THEN    ! Convection, take fluxes that are maximal in magnitude 
         IF(ABS(Q_sen_tur).GE.ABS(Q_sen_con)) THEN
            Q_sensible = Q_sen_tur
         ELSE
            Q_sensible = Q_sen_con
         ENDIF
         IF(ABS(Q_sensible).LT.ABS(Q_sen_mol)) THEN
            Q_sensible = Q_sen_mol
         ENDIF
         IF(ABS(Q_lat_tur).GE.ABS(Q_lat_con)) THEN
            Q_latent = Q_lat_tur
         ELSE
            Q_latent = Q_lat_con
         ENDIF
         IF(ABS(Q_latent).LT.ABS(Q_lat_mol)) THEN
            Q_latent = Q_lat_mol
         ENDIF
      ELSE                    ! Stable or neutral stratification, chose fluxes that are maximal in magnitude 
         IF(ABS(Q_sen_tur).GE.ABS(Q_sen_mol)) THEN 
            Q_sensible = Q_sen_tur
         ELSE 
            Q_sensible = Q_sen_mol    
         ENDIF
         IF(ABS(Q_lat_tur).GE.ABS(Q_lat_mol)) THEN 
            Q_latent = Q_lat_tur
         ELSE 
            Q_latent = Q_lat_mol  
         ENDIF
      ENDIF

!------------------------------------------------------------------------------
!  Set output (notice that fluxes are no longer in kinematic units)
!------------------------------------------------------------------------------
      Q_momentum = Q_momentum*rho_a 
      Q_sensible = Q_sensible*rho_a*tpsf_c_a_p
      Q_watvap   = Q_latent*rho_a
      Q_latent = tpsf_L_evap
      IF(h_ice.GE.h_Ice_min_flk) Q_latent = Q_latent + tpl_L_f   ! Add latent heat of fusion over ice
      Q_latent = Q_watvap*Q_latent
      !+WMEJ The net longwave radiation flux at the surface
      Q_olrg     = c_lwrad_emis*stefnc*T_s**4                        

      !+WMEJ The net heat flux at the surface
      Q_fgrnd  = sabg + lwdns - Q_olrg - Q_sensible - Q_latent      
      z0m  = z0u_sf
      z0h  = z0t_sf
      z0q  = z0q_sf
      felak = fetch
      !-WMEJ it is not necessary to calculate
      ! Set "*_sf" variables to make fluxes accessible to driving routines that USE "SfcFlx"
      ! u_star_a_sf     = u_star_st 
      ! Q_mom_a_sf      = Q_momentum  
      ! Q_sens_a_sf     = Q_sensible 
      ! Q_lat_a_sf      = Q_latent
      ! Q_watvap_a_sf   = Q_watvap
!------------------------------------------------------------------------------
!  END calculations
!==============================================================================

   END SUBROUTINE FLake_SurfFlux
       


   SUBROUTINE flaketem ( &
               ! "in" arguments
               ! ---------------------------
               depth_w        ,   depth_bs      ,   par_Coriolis  ,   extincoef_water,&
               del_time       ,   I_w_flk       ,   I_h_flk       ,   I_intm_0_h_flk ,&
               Q_w_flk        ,   Q_snow_flk    ,   Q_ice_flk     ,   I_snow_flk     ,&
               I_ice_flk      ,   dMsnowdt_flk  ,   T_bs          ,   u_star_w_flk   ,&
               I_bot_flk      ,   I_intm_h_D_flk,&
               ! "inout" arguments
               ! -------------------
               T_sfc_p        ,   T_snow_p_flk  ,   T_ice_p_flk   ,   T_wML_p_flk   ,&
               T_mnw_p_flk    ,   T_bot_p_flk   ,   T_B1_p_flk    ,   h_snow_p_flk  ,&
               h_ice_p_flk    ,   h_ML_p_flk    ,   H_B1_p_flk    ,   C_T_p_flk     ,&
               ! "out" arguments
               ! -------------------
               h_snow_dt      ,   h_ice_dt      )

! ---------------------------------- code history -------------------------------------
! Description:
!
!  The main driving routine of the lake model FLake 
!  where computations are performed.
!  Advances the surface temperature
!  and other FLake variables one time step.
!  At the moment, the Euler explicit scheme is used.
!
!  Lines embraced with "!_tmp" contain temporary parts of the code.
!  Lines embraced/marked with "!_dev" may be replaced
!  as improved parameterizations are developed and tested.
!  Lines embraced/marked with "!_dm" are DM's comments
!  that may be helpful to a user.
!  Lines embraced/marked with "!_dbg" are used 
!  for debugging purposes only.
!
! Original author: 
!     Dmitrii Mironov, 2005/11/17
!
!==============================================================================
   USE MOD_Lake_Const, only: tpl_kappa_w       , h_Ice_min_flk  , h_ML_min_flk ,&
                              Phi_T_pr0_1      , Phi_T_pr0_2    ,&
                              tpl_T_f          , c_small_flk    , tpl_rho_I    ,&
                              tpl_kappa_I      , Phi_I_pr0_lin  , tpl_rho_S_min,&
                              Phi_I_pr1_lin    , Phi_I_ast_MR   , H_Ice_max    ,&
                              Phi_S_pr0_lin    , h_Snow_min_flk , tpl_rho_S_max,&
                              tpl_Gamma_rho_S  , tpl_rho_w_r    , C_I_lin      ,&
                              C_I_MR           , tpl_c_I        , tpl_c_S      ,&
                              C_T_min          , lflk_botsed_use, tpl_L_f      ,&
                              rflk_depth_bs_ref, C_S_lin        , H_B1_min_flk ,&
                              Phi_B1_pr0       , tpl_c_w        , tpl_T_r      ,&
                              C_T_max          , u_star_min_flk , c_relax_C    ,&
                              C_TT_1           , C_TT_2         , c_cbl_1      ,&
                              c_cbl_2          , c_sbl_ZM_n     , c_sbl_ZM_i   ,&
                              c_sbl_ZM_s       , h_ML_max_flk   , c_relax_h    ,&
                              C_B1             , C_B2           , c_z0q_rough_1,&
                              c_z0q_rough_2    , c_z0q_rough_3
!================================================================================
!  ------------------------- input variables ---------------------------
   real(r8), intent(in)     :: &
      depth_w                 ,&! The lake depth [m]
      depth_bs                ,&! Depth of the thermally active layer of bottom sediments [m]
      T_bs                    ,&! Temperature at the outer edge of the thermally active layer of bottom sediments [K]
      par_Coriolis            ,&! The Coriolis parameter [s^{-1}]
      extincoef_water         ,&! "Typical" extinction coefficient of the lake water [m^{-1}], used to compute the equilibrium CBL depth
      del_time                ,&! The model time step [s]
      I_w_flk                 ,&! Radiation flux through the ice-water or air-water interface [W m^{-2}]
      I_h_flk                 ,&! Radiation flux through the mixed-layer-thermocline interface [W m^{-2}]
      I_intm_0_h_flk          ,&! Mean radiation flux over the mixed layer [W m^{-1}]
      I_intm_h_D_flk          ,&! Mean radiation flux over the thermocline [W m^{-1}]
      Q_snow_flk              ,&! Radiation flux through the air-snow interface [W m^{-2}]
      Q_ice_flk               ,&! Radiation flux through the snow-ice or air-ice interface [W m^{-2}]
      I_snow_flk              ,&! Radiation flux through the air-snow interface [W m^{-2}]
      I_ice_flk               ,&! Ice radiation [W m^{-2}]
      I_bot_flk               ,&! Radiation flux through the water-bottom sediment interface [W m^{-2}]
      u_star_w_flk            ,&! Friction velocity over water [m s^{-1}]
      dMsnowdt_flk              ! The rate of snow accumulation [kg m^{-2} s^{-1}]

!  ------------------------- inout variables ---------------------------
   real(r8), intent(inout)  :: &
      Q_w_flk                 ,&! Heat flux through the ice-water or air-water interface [W m^{-2}]
      T_sfc_p                 ,&! Previous surface temperature [K] (equal to the updated value of either T_ice, T_snow or T_wML)
      T_snow_p_flk            ,&! Previous Temperature at the air-snow interface [K] 
      T_ice_p_flk             ,&! Previous Temperature at the snow-ice or air-ice interface [K] 
      T_wML_p_flk             ,&! Previous Mixed-layer temperature [K] 
      T_mnw_p_flk             ,&! Previous Mean temperature of the water column [K] 
      T_bot_p_flk             ,&! Previous Temperature at the water-bottom sediment interface [K] 
      T_B1_p_flk              ,&! Previous Temperature at the bottom of the upper layer of the sediments [K] 
      h_snow_p_flk            ,&! Previous Snow thickness [m]
      h_ice_p_flk             ,&! Previous Ice thickness [m]
      h_ML_p_flk              ,&! Previous Thickness of the mixed-layer [m] 
      H_B1_p_flk              ,&! Previous Thickness of the upper layer of bottom sediments [m] 
      C_T_p_flk                 ! Shape factor (thermocline)

!  ------------------------- output variables ---------------------------
   real(r8), intent(out)    :: &
      h_snow_dt               ,&! Time derivative of h_snow [m s^{-1}]
      h_ice_dt                  ! Time derivative of h_ice [m s^{-1}]

!  ------------------------- local variables ---------------------------
   real(r8)                 :: &
      Q_star_flk              ,&! A generalized heat flux scale  [W m^{-2}]
      Q_bot_flk               ,&! Heat flux through the water-bottom sediment interface
      T_sfc_n                 ,&! Updated surface temperature [K] (equal to the updated value of either T_ice, T_snow or T_wML)
      T_snow_n_flk            ,&! Temperature at the air-snow interface [K] 
      T_ice_n_flk             ,&! Temperature at the snow-ice or air-ice interface [K] 
      T_wML_n_flk             ,&! Updated mixed-layer temperature [K]
      T_mnw_n_flk             ,&! Updated Mean temperature of the water column [K] 
      T_bot_n_flk             ,&! Updated Temperature at the water-bottom sediment interface [K] 
      T_B1_n_flk              ,&! Updated Temperature at the bottom of the upper layer of the sediments [K] 
      h_snow_n_flk            ,&! Updated Snow thickness [m]
      h_ice_n_flk             ,&! Updated Ice thickness [m]
      h_ML_n_flk              ,&! Updated Thickness of the mixed-layer [m] 
      H_B1_n_flk              ,&! Updated Thickness of the upper layer of bottom sediments [m] 
      C_T_n_flk                 ! Updated Shape factor (thermocline)

   real(r8)                 :: &
      d_T_mnw_dt              ,&! Time derivative of T_mnw [K s^{-1}] 
      d_T_ice_dt              ,&! Time derivative of T_ice [K s^{-1}] 
      d_T_bot_dt              ,&! Time derivative of T_bot [K s^{-1}] 
      d_T_B1_dt               ,&! Time derivative of T_B1 [K s^{-1}] 
      d_h_snow_dt             ,&! Time derivative of h_snow [m s^{-1}]
      d_h_ice_dt              ,&! Time derivative of h_ice [m s^{-1}]
      d_h_ML_dt               ,&! Time derivative of h_ML [m s^{-1}]
      d_H_B1_dt               ,&! Time derivative of H_B1 [m s^{-1}]
      d_C_T_dt                ,&! Time derivative of C_T [s^{-1}]
      C_TT_flk                ,&! Updated Shape factor (thermocline)
      C_Q_flk                 ,&! Updated Shape factor (thermocline)
      N_T_mean                ,&! The mean buoyancy frequency in the thermocline [s^{-1}] 
      ZM_h_scale              ,&! The ZM96 equilibrium SBL depth scale [m] 
      w_star_sfc_flk          ,&! The surface friction velocity [m s^{-1}]
      conv_equil_h_scale      ,&! The equilibrium CBL depth scale [m]
      h_ice_threshold         ,&! IF h_ice<h_ice_threshold, USE quasi-equilibrium ice model 
      flk_str_1               ,&! Help storage variable
      flk_str_2               ,&! Help storage variable
      R_H_icesnow             ,&! Dimensionless ratio, used to store intermediate results
      R_rho_c_icesnow         ,&! Dimensionless ratio, used to store intermediate results
      R_TI_icesnow            ,&! Dimensionless ratio, used to store intermediate results
      R_Tstar_icesnow         ,&! Dimensionless ratio, used to store intermediate results
      Phi_T_pr0_flk           ,&! d\Phi(0)/d\zeta (thermocline)
      Phi_I_pr0_flk           ,&! d\Phi(0)/d\zeta (ice)
      Phi_I_pr1_flk           ,&! d\Phi(1)/d\zeta (ice)
      C_I_flk                   ! Updated Shape factor (ice)

   logical                  :: &
      l_ice_create            ,&! Switch, .TRUE. = ice does not exist but should be created
      l_snow_exists           ,&! Switch, .TRUE. = there is snow above the ice
      l_ice_meltabove           ! Switch, .TRUE. = snow/ice melting from above takes place

   integer  :: i                ! DO loop index
!================================================================================


!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------
      !_dm 
      ! Security. Set time-rate-of-change of prognostic variables to zero.
      ! Set prognostic variables to their values at the previous time step.
      ! (This is to avoid spurious changes of prognostic variables 
      ! when FLake is used within a 3D model, e.g. to avoid spurious generation of ice 
      ! at the neighbouring lake points as noticed by Burkhardt Rockel.)
      !_dm 

      d_T_mnw_dt   = 0. 
      d_T_ice_dt   = 0. 
      d_T_bot_dt   = 0. 
      d_T_B1_dt    = 0. 
      d_h_snow_dt  = 0. 
      d_h_ice_dt   = 0. 
      d_h_ML_dt    = 0. 
      d_H_B1_dt    = 0. 
      d_C_T_dt     = 0. 
      T_snow_n_flk = T_snow_p_flk   
      T_ice_n_flk  = T_ice_p_flk    
      T_wML_n_flk  = T_wML_p_flk   
      T_mnw_n_flk  = T_mnw_p_flk     
      T_bot_n_flk  = T_bot_p_flk  
      T_B1_n_flk   = T_B1_p_flk      
      h_snow_n_flk = h_snow_p_flk 
      h_ice_n_flk  = h_ice_p_flk   
      h_ML_n_flk   = h_ML_p_flk    
      H_B1_n_flk   = H_B1_p_flk   
      C_T_n_flk    = C_T_p_flk    

!------------------------------------------------------------------------------
!  Compute fluxes, using variables from the previous time step.
!------------------------------------------------------------------------------

   !_dm
   ! At this point, the heat and radiation fluxes, namely,
   ! Q_snow_flk, Q_ice_flk, Q_w_flk, 
   ! I_atm_flk, I_snow_flk, I_ice_flk, I_w_flk, I_h_flk, I_bot_flk,     
   ! the mean radiation flux over the mixed layer, I_intm_0_h_flk, 
   ! and the mean radiation flux over the thermocline, I_intm_h_D_flk, 
   ! should be known.
   ! They are computed within "flake_interface" (or within the driving model)
   ! and are available to "flake_driver"
   ! through the above variables declared in the MODULE  "flake".
   ! in case a lake is ice-covered, Q_w_flk is re-computed below.
   !_dm

      ! Heat flux through the ice-water interface
      IF(h_ice_p_flk.GE.h_Ice_min_flk) THEN    ! Ice exists 
         IF(h_ML_p_flk.LE.h_ML_min_flk) THEN    ! Mixed-layer depth is zero, compute flux 
            Q_w_flk = -tpl_kappa_w*(T_bot_p_flk-T_wML_p_flk)/depth_w  ! Flux with linear T(z) 
            Phi_T_pr0_flk = Phi_T_pr0_1*C_T_p_flk-Phi_T_pr0_2         ! d\Phi(0)/d\zeta (thermocline)
            Q_w_flk = Q_w_flk*MAX(Phi_T_pr0_flk, 1.)           ! Account for an increased d\Phi(0)/d\zeta 
         ELSE                    
            Q_w_flk = 0.                  ! Mixed-layer depth is greater than zero, set flux to zero
         ENDIF   
      ENDIF   

      ! A generalized heat flux scale 
      Q_star_flk = Q_w_flk + I_w_flk + I_h_flk - 2.*I_intm_0_h_flk

      ! Heat flux through the water-bottom sediment interface
      IF(lflk_botsed_use) THEN
         Q_bot_flk = -tpl_kappa_w*(T_B1_p_flk-T_bot_p_flk)/MAX(H_B1_p_flk, H_B1_min_flk)*Phi_B1_pr0
      ELSE  
         Q_bot_flk = 0.   ! The bottom-sediment scheme is not used
      ENDIF

!------------------------------------------------------------------------------
!  Check IF ice exists or should be created.
!  IF so, compute the thickness and the temperature of ice and snow.
!------------------------------------------------------------------------------

      !_dm
      ! Notice that a quasi-equilibrium ice-snow model is used 
      ! to avoid numerical instability when the ice is thin.
      ! This is always the case when new ice is created.
      !_dm

      !_dev
      ! The dependence of snow density and of snow heat conductivity 
      ! on the snow thickness is accounted for parametrically.
      ! That is, the time derivatives of \rho_S and \kappa_S are neglected.
      ! The exception is the equation for the snow thickness 
      ! in case of snow accumulation and no melting, 
      ! where d\rho_S/dt is incorporated.
      ! Furthermore, some (presumably small) correction terms incorporating 
      ! the snow density and the snow heat conductivity are dropped out.
      ! Those terms may be included as better formulations 
      ! for \rho_S and \kappa_S are available.
      !_dev

      ! Default values
      l_ice_create    = .FALSE.  
      l_ice_meltabove = .FALSE.  

      Ice_exist: IF(h_ice_p_flk.LT.h_Ice_min_flk) THEN   ! Ice does not exist 

         l_ice_create = T_wML_p_flk.LE.(tpl_T_f+c_small_flk).AND.Q_w_flk.LT.0.
         IF(l_ice_create) THEN                            ! Ice does not exist but should be created
            d_h_ice_dt = -Q_w_flk/tpl_rho_I/tpl_L_f                                  
            h_ice_n_flk = h_ice_p_flk + d_h_ice_dt*del_time                          ! Advance h_ice 
            T_ice_n_flk = tpl_T_f + h_ice_n_flk*Q_w_flk/tpl_kappa_I/Phi_I_pr0_lin    ! Ice temperature
            d_h_snow_dt = dMsnowdt_flk/tpl_rho_S_min 
            h_snow_n_flk = h_snow_p_flk + d_h_snow_dt*del_time                       ! Advance h_snow
            Phi_I_pr1_flk = Phi_I_pr1_lin                                    & 
                        + Phi_I_ast_MR*MIN(1., h_ice_n_flk/H_Ice_max)       ! d\Phi_I(1)/d\zeta_I (ice)
            R_H_icesnow = Phi_I_pr1_flk/Phi_S_pr0_lin*tpl_kappa_I/flake_snowheatconduct(h_snow_n_flk) &
                        * h_snow_n_flk/MAX(h_ice_n_flk, h_Ice_min_flk)
            T_snow_n_flk = T_ice_n_flk + R_H_icesnow*(T_ice_n_flk-tpl_T_f)           ! Snow temperature
         ENDIF

      ELSE Ice_exist                                     ! Ice exists
         l_snow_exists = h_snow_p_flk.GE.h_Snow_min_flk   ! Check IF there is snow above the ice

         Melting: IF(T_snow_p_flk.GE.(tpl_T_f-c_small_flk)) THEN  ! T_sfc = T_f, check for melting from above
                                                                     ! T_snow = T_ice IF snow is absent 
            IF(l_snow_exists) THEN   ! There is snow above the ice
               flk_str_1 = Q_snow_flk + I_snow_flk - I_ice_flk        ! Atmospheric forcing
               IF(flk_str_1.GE.0.) THEN  ! Melting of snow and ice from above
                  l_ice_meltabove = .TRUE.
                  d_h_snow_dt = (-flk_str_1/tpl_L_f+dMsnowdt_flk)/flake_snowdensity(h_snow_p_flk)
                  d_h_ice_dt  = -(I_ice_flk - I_w_flk - Q_w_flk)/tpl_L_f/tpl_rho_I 
               ENDIF 
            ELSE                     ! No snow above the ice
               flk_str_1 = Q_ice_flk + I_ice_flk - I_w_flk - Q_w_flk  ! Atmospheric forcing + heating from the water
               IF(flk_str_1.GE.0.) THEN  ! Melting of ice from above, snow accumulation may occur
                  l_ice_meltabove = .TRUE.
                  d_h_ice_dt  = -flk_str_1/tpl_L_f/tpl_rho_I 
                  d_h_snow_dt = dMsnowdt_flk/tpl_rho_S_min
               ENDIF 
            ENDIF 
            IF(l_ice_meltabove) THEN  ! Melting from above takes place
               h_ice_n_flk  = h_ice_p_flk  + d_h_ice_dt *del_time  ! Advance h_ice
               h_snow_n_flk = h_snow_p_flk + d_h_snow_dt*del_time  ! Advance h_snow
               T_ice_n_flk  = tpl_T_f                              ! Set T_ice to the freezing point
               T_snow_n_flk = tpl_T_f                              ! Set T_snow to the freezing point
            ENDIF

         ENDIF Melting

         No_Melting: IF(.NOT.l_ice_meltabove) THEN                 ! No melting from above

            d_h_snow_dt = flake_snowdensity(h_snow_p_flk)  
            IF(d_h_snow_dt.LT.tpl_rho_S_max) THEN    ! Account for d\rho_S/dt
               flk_str_1 = h_snow_p_flk*tpl_Gamma_rho_S/tpl_rho_w_r
               flk_str_1 = flk_str_1/(1.-flk_str_1)
            ELSE                                     ! Snow density is equal to its maximum value, d\rho_S/dt=0
               flk_str_1 = 0.
            ENDIF
            d_h_snow_dt = dMsnowdt_flk/d_h_snow_dt/(1.+flk_str_1)       ! Snow accumulation
            h_snow_n_flk = h_snow_p_flk + d_h_snow_dt*del_time                         ! Advance h_snow
            
            Phi_I_pr0_flk = h_ice_p_flk/H_Ice_max                              ! h_ice relative to its maximum value
            C_I_flk = C_I_lin - C_I_MR*(1.+Phi_I_ast_MR)*Phi_I_pr0_flk  ! Shape factor (ice)
            Phi_I_pr1_flk = Phi_I_pr1_lin + Phi_I_ast_MR*Phi_I_pr0_flk         ! d\Phi_I(1)/d\zeta_I (ice)
            Phi_I_pr0_flk = Phi_I_pr0_lin - Phi_I_pr0_flk                      ! d\Phi_I(0)/d\zeta_I (ice)

            h_ice_threshold = MAX(1., 2.*C_I_flk*tpl_c_I*(tpl_T_f-T_ice_p_flk)/tpl_L_f)
            h_ice_threshold = Phi_I_pr0_flk/C_I_flk*tpl_kappa_I/tpl_rho_I/tpl_c_I*h_ice_threshold
            h_ice_threshold = SQRT(h_ice_threshold*del_time)                   ! Threshold value of h_ice
            h_ice_threshold = MIN(0.9       *H_Ice_max, MAX(h_ice_threshold, h_Ice_min_flk))
                                                                           ! h_ice(threshold) < 0.9*H_Ice_max

            IF(h_ice_p_flk.LT.h_ice_threshold) THEN  ! USE a quasi-equilibrium ice model

               IF(l_snow_exists) THEN   ! USE fluxes at the air-snow interface
                  flk_str_1 = Q_snow_flk + I_snow_flk - I_w_flk
               ELSE                     ! USE fluxes at the air-ice interface
                  flk_str_1 = Q_ice_flk + I_ice_flk - I_w_flk
               ENDIF
               d_h_ice_dt = -(flk_str_1-Q_w_flk)/tpl_L_f/tpl_rho_I
               h_ice_n_flk = h_ice_p_flk + d_h_ice_dt *del_time                         ! Advance h_ice
               T_ice_n_flk = tpl_T_f + h_ice_n_flk*flk_str_1/tpl_kappa_I/Phi_I_pr0_flk  ! Ice temperature

            ELSE                                     ! USE a complete ice model

               d_h_ice_dt = tpl_kappa_I*(tpl_T_f-T_ice_p_flk)/h_ice_p_flk*Phi_I_pr0_flk
               d_h_ice_dt = (Q_w_flk+d_h_ice_dt)/tpl_L_f/tpl_rho_I
               h_ice_n_flk = h_ice_p_flk  + d_h_ice_dt*del_time                         ! Advance h_ice

               R_TI_icesnow = tpl_c_I*(tpl_T_f-T_ice_p_flk)/tpl_L_f         ! Dimensionless parameter
               R_Tstar_icesnow = 1. - C_I_flk                        ! Dimensionless parameter
               IF(l_snow_exists) THEN  ! There is snow above the ice
                  R_H_icesnow = Phi_I_pr1_flk/Phi_S_pr0_lin*tpl_kappa_I/flake_snowheatconduct(h_snow_p_flk) &
                              * h_snow_p_flk/h_ice_p_flk
                  R_rho_c_icesnow = flake_snowdensity(h_snow_p_flk)*tpl_c_S/tpl_rho_I/tpl_c_I 
               ! These terms should be included as an improved understanding of the snow scheme is gained, 
               ! of the effect of snow density in particular. 
               !_nu        R_Tstar_icesnow = R_Tstar_icesnow                                                           &
               !_nu                        + (1.+C_S_lin*h_snow_p_flk/h_ice_p_flk)*R_H_icesnow*R_rho_c_icesnow

                  R_Tstar_icesnow = R_Tstar_icesnow*R_TI_icesnow             ! Dimensionless parameter
               !_nu        R_Tstar_icesnow = R_Tstar_icesnow                                                         &
               !_nu                        + (1.-R_rho_c_icesnow)*tpl_c_I*T_ice_p_flk/tpl_L_f
                  flk_str_2 = Q_snow_flk+I_snow_flk-I_w_flk                  ! Atmospheric fluxes
                  flk_str_1  = C_I_flk*h_ice_p_flk + (1.+C_S_lin*R_H_icesnow)*R_rho_c_icesnow*h_snow_p_flk
                  d_T_ice_dt = -(1.-2.*C_S_lin)*R_H_icesnow*(tpl_T_f-T_ice_p_flk)             & 
                           * tpl_c_S*dMsnowdt_flk                          ! Effect of snow accumulation
               ELSE                    ! No snow above the ice
                  R_Tstar_icesnow = R_Tstar_icesnow*R_TI_icesnow             ! Dimensionless parameter
                  flk_str_2 = Q_ice_flk+I_ice_flk-I_w_flk                    ! Atmospheric fluxes
                  flk_str_1  = C_I_flk*h_ice_p_flk
                  d_T_ice_dt = 0.
               ENDIF 
               d_T_ice_dt = d_T_ice_dt + tpl_kappa_I*(tpl_T_f-T_ice_p_flk)/h_ice_p_flk*Phi_I_pr0_flk       &
                           * (1.-R_Tstar_icesnow)                     ! Add flux due to heat conduction
               d_T_ice_dt = d_T_ice_dt - R_Tstar_icesnow*Q_w_flk            ! Add flux from water to ice
               d_T_ice_dt = d_T_ice_dt + flk_str_2                          ! Add atmospheric fluxes
               d_T_ice_dt = d_T_ice_dt/tpl_rho_I/tpl_c_I                    ! Total forcing
               d_T_ice_dt = d_T_ice_dt/flk_str_1                            ! dT_ice/dt 
               T_ice_n_flk = T_ice_p_flk + d_T_ice_dt*del_time                          ! Advance T_ice
            ENDIF

            Phi_I_pr1_flk = MIN(1., h_ice_n_flk/H_Ice_max)          ! h_ice relative to its maximum value
            Phi_I_pr1_flk = Phi_I_pr1_lin + Phi_I_ast_MR*Phi_I_pr1_flk     ! d\Phi_I(1)/d\zeta_I (ice)
            R_H_icesnow = Phi_I_pr1_flk/Phi_S_pr0_lin*tpl_kappa_I/flake_snowheatconduct(h_snow_n_flk) &
                        *h_snow_n_flk/MAX(h_ice_n_flk, h_Ice_min_flk)
            T_snow_n_flk = T_ice_n_flk + R_H_icesnow*(T_ice_n_flk-tpl_T_f)             ! Snow temperature

         ENDIF No_Melting

      ENDIF Ice_exist   

      ! Security, limit h_ice by its maximum value
      h_ice_n_flk = MIN(h_ice_n_flk, H_Ice_max)      

      ! Security, limit the ice and snow temperatures by the freezing point 
      T_snow_n_flk = MIN(T_snow_n_flk, tpl_T_f)  
      T_ice_n_flk =  MIN(T_ice_n_flk,  tpl_T_f)    

      !_tmp
      ! Security, avoid too low values (these constraints are used for debugging purposes)
      T_snow_n_flk = MAX(T_snow_n_flk, 73.15       )  
      T_ice_n_flk =  MAX(T_ice_n_flk,  73.15       )    
      !_tmp

      ! Remove too thin ice and/or snow
      IF(h_ice_n_flk.LT.h_Ice_min_flk)  THEN        ! Check ice
         h_ice_n_flk = 0.       ! Ice is too thin, remove it, and
         T_ice_n_flk = tpl_T_f         ! set T_ice to the freezing point.
         h_snow_n_flk = 0.      ! Remove snow when there is no ice, and
         T_snow_n_flk = tpl_T_f        ! set T_snow to the freezing point.
         l_ice_create = .FALSE.        ! "Exotic" case, ice has been created but proved to be too thin
      ELSE IF(h_snow_n_flk.LT.h_Snow_min_flk) THEN  ! Ice exists, check snow
         h_snow_n_flk = 0.      ! Snow is too thin, remove it, 
         T_snow_n_flk = T_ice_n_flk    ! and set the snow temperature equal to the ice temperature.
      ENDIF
      h_snow_dt = d_h_snow_dt * 1000.!-----------
      h_ice_dt = d_h_ice_dt * 1000.  !-----------
!------------------------------------------------------------------------------
!  Compute the mean temperature of the water column.
!------------------------------------------------------------------------------

      IF(l_ice_create) Q_w_flk = 0.     ! Ice has just been created, set Q_w to zero
      d_T_mnw_dt = (Q_w_flk - Q_bot_flk + I_w_flk - I_bot_flk)/tpl_rho_w_r/tpl_c_w/depth_w
      T_mnw_n_flk = T_mnw_p_flk + d_T_mnw_dt*del_time   ! Advance T_mnw
      T_mnw_n_flk = MAX(T_mnw_n_flk, tpl_T_f)           ! Limit T_mnw by the freezing point 
      !------------------------------------------------------------------------------
      !  Compute the mixed-layer depth, the mixed-layer temperature, 
      !  the bottom temperature and the shape factor
      !  with respect to the temperature profile in the thermocline. 
      !  Different formulations are used, depending on the regime of mixing. 
      !------------------------------------------------------------------------------

      HTC_Water: IF(h_ice_n_flk.GE.h_Ice_min_flk) THEN    ! Ice exists

         T_mnw_n_flk = MIN(T_mnw_n_flk, tpl_T_r) ! Limit the mean temperature under the ice by T_r 
         T_wML_n_flk = tpl_T_f                   ! The mixed-layer temperature is equal to the freezing point 

         IF(l_ice_create) THEN                  ! Ice has just been created 
            IF(h_ML_p_flk.GE.depth_w-h_ML_min_flk) THEN    ! h_ML=D when ice is created 
               h_ML_n_flk = 0.                 ! Set h_ML to zero 
               C_T_n_flk = C_T_min                    ! Set C_T to its minimum value 
            ELSE                                          ! h_ML<D when ice is created 
               h_ML_n_flk = h_ML_p_flk                ! h_ML remains unchanged 
               C_T_n_flk = C_T_p_flk                  ! C_T (thermocline) remains unchanged 
            ENDIF 
            T_bot_n_flk = T_wML_n_flk - (T_wML_n_flk-T_mnw_n_flk)/C_T_n_flk/(1.-h_ML_n_flk/depth_w)
                                                   ! Update the bottom temperature 

         ELSE IF(T_bot_p_flk.LT.tpl_T_r) THEN   ! Ice exists and T_bot < T_r, molecular heat transfer 
            h_ML_n_flk = h_ML_p_flk                  ! h_ML remains unchanged 
            C_T_n_flk = C_T_p_flk                    ! C_T (thermocline) remains unchanged 
            T_bot_n_flk = T_wML_n_flk - (T_wML_n_flk-T_mnw_n_flk)/C_T_n_flk/(1.-h_ML_n_flk/depth_w)
                                                   ! Update the bottom temperature 

         ELSE                                   ! Ice exists and T_bot = T_r, convection due to bottom heating 
            T_bot_n_flk = tpl_T_r                      ! T_bot is equal to the temperature of maximum density 
            IF(h_ML_p_flk.GE.c_small_flk) THEN   ! h_ML > 0 
               C_T_n_flk = C_T_p_flk                     ! C_T (thermocline) remains unchanged 
               h_ML_n_flk = depth_w*(1.-(T_wML_n_flk-T_mnw_n_flk)/(T_wML_n_flk-T_bot_n_flk)/C_T_n_flk)
               h_ML_n_flk = MAX(h_ML_n_flk, 0.)   ! Update the mixed-layer depth  
            ELSE                                 ! h_ML = 0 
               h_ML_n_flk = h_ML_p_flk                   ! h_ML remains unchanged 
               C_T_n_flk = (T_wML_n_flk-T_mnw_n_flk)/(T_wML_n_flk-T_bot_n_flk) 
               C_T_n_flk = MIN(C_T_max, MAX(C_T_n_flk, C_T_min)) ! Update the shape factor (thermocline)  
            ENDIF 
         ENDIF 

         T_bot_n_flk = MIN(T_bot_n_flk, tpl_T_r)    ! Security, limit the bottom temperature by T_r 

      ELSE HTC_Water                                      ! Open water

         ! Generalised buoyancy flux scale and convective velocity scale
         flk_str_1 = flake_buoypar(T_wML_p_flk)*Q_star_flk/tpl_rho_w_r/tpl_c_w                    
         IF(flk_str_1.LT.0.) THEN       
            w_star_sfc_flk = (-flk_str_1*h_ML_p_flk)**(1./3.)  ! Convection     
         ELSE 
            w_star_sfc_flk = 0.                                       ! Neutral or stable stratification
         ENDIF 

         !_dm
         ! The equilibrium depth of the CBL due to surface cooling with the volumetric heating
         ! is not computed as a solution to the transcendental equation.
         ! Instead, an algebraic formula is used
         ! that interpolates between the two asymptotic limits.
         !_dm
         conv_equil_h_scale = -Q_w_flk/MAX(I_w_flk, c_small_flk)
         IF(conv_equil_h_scale.GT.0. .AND. conv_equil_h_scale.LT.1.  &
            .AND. T_wML_p_flk.GT.tpl_T_r) THEN   ! The equilibrium CBL depth scale is only used above T_r
            conv_equil_h_scale = SQRT(6.*conv_equil_h_scale)                 &
                           + 2.*conv_equil_h_scale/(1.-conv_equil_h_scale)
            conv_equil_h_scale = MIN(depth_w, conv_equil_h_scale/extincoef_water    )
         ELSE
            conv_equil_h_scale = 0.       ! Set the equilibrium CBL depth to zero
         ENDIF

         ! Mean buoyancy frequency in the thermocline

         N_T_mean = flake_buoypar(0.5 * (T_wML_p_flk+T_bot_p_flk)) * (T_wML_p_flk-T_bot_p_flk)
         IF(h_ML_p_flk.LE.depth_w-h_ML_min_flk) THEN
            N_T_mean = SQRT(N_T_mean/(depth_w-h_ML_p_flk))  ! Compute N                   
         ELSE 
            N_T_mean = 0.                            ! h_ML=D, set N to zero
         ENDIF 

         ! The rate of change of C_T
         d_C_T_dt = MAX(w_star_sfc_flk, u_star_w_flk, u_star_min_flk)**2
         d_C_T_dt = N_T_mean*(depth_w-h_ML_p_flk)**2       &
                  / c_relax_C/d_C_T_dt                               ! Relaxation time scale for C_T
         d_C_T_dt = (C_T_max-C_T_min)/MAX(d_C_T_dt, c_small_flk)     ! Rate-of-change of C_T 

         ! Compute the shape factor and the mixed-layer depth, 
         ! using different formulations for convection and wind mixing

         C_TT_flk = C_TT_1*C_T_p_flk-C_TT_2         ! C_TT, using C_T at the previous time step
         C_Q_flk = 2.*C_TT_flk/C_T_p_flk     ! C_Q using C_T at the previous time step

         Mixing_regime: IF(flk_str_1.LT.0.) THEN  ! Convective mixing 

            C_T_n_flk = C_T_p_flk + d_C_T_dt*del_time                        ! Update C_T, assuming dh_ML/dt>0
            C_T_n_flk = MIN(C_T_max, MAX(C_T_n_flk, C_T_min))                ! Limit C_T 
            d_C_T_dt = (C_T_n_flk-C_T_p_flk)/del_time                        ! Re-compute dC_T/dt

            IF(h_ML_p_flk.LE.depth_w-h_ML_min_flk) THEN       ! Compute dh_ML/dt
               IF(h_ML_p_flk.LE.h_ML_min_flk) THEN    ! USE a reduced entrainment equation (spin-up)
                  d_h_ML_dt = c_cbl_1/c_cbl_2*MAX(w_star_sfc_flk, c_small_flk)

! #ifdef LKDEBUG
!             PRINT *, ' FLake: reduced entrainment eq. D_time*d_h_ML_dt  = ', d_h_ML_dt*del_time
!             PRINT *, '         w_*       = ', w_star_sfc_flk
!             PRINT *, '         \beta*Q_* = ', flk_str_1
! #endif

               ELSE                                   ! USE a complete entrainment equation 
                  R_H_icesnow     = depth_w/h_ML_p_flk
                  R_rho_c_icesnow = R_H_icesnow-1.
                  R_TI_icesnow    = C_T_p_flk/C_TT_flk
                  R_Tstar_icesnow = (R_TI_icesnow/2.-1.)*R_rho_c_icesnow + 1.
                  d_h_ML_dt = -Q_star_flk*(R_Tstar_icesnow*(1.+c_cbl_1)-1.) - Q_bot_flk
                  d_h_ML_dt = d_h_ML_dt/tpl_rho_w_r/tpl_c_w                        ! Q_* and Q_b flux terms
                  flk_str_2 = (depth_w-h_ML_p_flk)*(T_wML_p_flk-T_bot_p_flk)*C_TT_2/C_TT_flk*d_C_T_dt 
                  d_h_ML_dt = d_h_ML_dt + flk_str_2                                 ! Add dC_T/dt term
                  flk_str_2 = I_bot_flk + (R_TI_icesnow-1.)*I_h_flk - R_TI_icesnow*I_intm_h_D_flk
                  flk_str_2 = flk_str_2 + (R_TI_icesnow-2.)*R_rho_c_icesnow*(I_h_flk-I_intm_0_h_flk)
                  flk_str_2 = flk_str_2/tpl_rho_w_r/tpl_c_w
                  d_h_ML_dt = d_h_ML_dt + flk_str_2                                 ! Add radiation terms
                  flk_str_2 = -c_cbl_2*R_Tstar_icesnow*Q_star_flk/tpl_rho_w_r/tpl_c_w/MAX(w_star_sfc_flk, c_small_flk)
                  flk_str_2 = flk_str_2 + C_T_p_flk*(T_wML_p_flk-T_bot_p_flk)
                  d_h_ML_dt = d_h_ML_dt/flk_str_2                                   ! dh_ML/dt = r.h.s.
               ENDIF 
         !_dm
         ! Notice that dh_ML/dt may appear to be negative  
         ! (e.g. due to buoyancy loss to bottom sediments and/or
         ! the effect of volumetric radiation heating),
         ! although a negative generalized buoyancy flux scale indicates 
         ! that the equilibrium CBL depth has not yet been reached
         ! and convective deepening of the mixed layer should take place.
         ! Physically, this situation reflects an approximate character of the lake model.
         ! Using the self-similar temperature profile in the thermocline, 
         ! there is always communication between the mixed layer, the thermocline 
         ! and the lake bottom. As a result, the rate of change of the CBL depth
         ! is always dependent on the bottom heat flux and the radiation heating of the thermocline.
         ! in reality, convective mixed-layer deepening may be completely decoupled
         ! from the processes underneath. in order to account for this fact,
         ! the rate of CBL deepening is set to a small value
         ! IF dh_ML/dt proves to be negative.
         ! This is "double insurance" however, 
         ! as a negative dh_ML/dt is encountered very rarely.
         !_dm

               d_h_ML_dt = MAX(d_h_ML_dt, c_small_flk)    
               h_ML_n_flk = h_ML_p_flk + d_h_ML_dt*del_time                       ! Update h_ML 
               h_ML_n_flk = MAX(h_ML_min_flk, MIN(h_ML_n_flk, depth_w))           ! Security, limit h_ML
            ELSE                                              ! Mixing down to the lake bottom
               h_ML_n_flk = depth_w
            ENDIF

         ELSE Mixing_regime                              ! Wind mixing

            d_h_ML_dt = MAX(u_star_w_flk, u_star_min_flk)                        ! The surface friction velocity
            ZM_h_scale = (ABS(par_Coriolis)/c_sbl_ZM_n + N_T_mean/c_sbl_ZM_i)*d_h_ML_dt**2
            ZM_h_scale = ZM_h_scale + flk_str_1/c_sbl_ZM_s
            ZM_h_scale = MAX(ZM_h_scale, c_small_flk)
            ZM_h_scale = d_h_ML_dt**3/ZM_h_scale 
            ZM_h_scale = MAX(h_ML_min_flk, MIN(ZM_h_scale, h_ML_max_flk))        ! The ZM96 SBL depth scale 
            ZM_h_scale = MAX(ZM_h_scale, conv_equil_h_scale)                     ! Equilibrium mixed-layer depth 

         !_dm 
         ! in order to avoid numerical discretization problems,
         ! an analytical solution to the evolution equation 
         ! for the wind-mixed layer depth is used.
         ! That is, an exponential relaxation formula is applied
         ! over the time interval equal to the model time step.
         !_dm 

            d_h_ML_dt = c_relax_h*d_h_ML_dt/ZM_h_scale*del_time
            h_ML_n_flk = ZM_h_scale - (ZM_h_scale-h_ML_p_flk)*EXP(-d_h_ML_dt)    ! Update h_ML 
            h_ML_n_flk = MAX(h_ML_min_flk, MIN(h_ML_n_flk, depth_w))             ! Limit h_ML 
            d_h_ML_dt = (h_ML_n_flk-h_ML_p_flk)/del_time                         ! Re-compute dh_ML/dt

            IF(h_ML_n_flk.LE.h_ML_p_flk)           &
            d_C_T_dt = -d_C_T_dt                 ! Mixed-layer retreat or stationary state, dC_T/dt<0
            C_T_n_flk = C_T_p_flk + d_C_T_dt*del_time                            ! Update C_T
            C_T_n_flk = MIN(C_T_max, MAX(C_T_n_flk, C_T_min))                    ! Limit C_T 
            d_C_T_dt = (C_T_n_flk-C_T_p_flk)/del_time                            ! Re-compute dC_T/dt

         ENDIF Mixing_regime

         ! Compute the time-rate-of-change of the the bottom temperature, 
         ! depending on the sign of dh_ML/dt 
         ! Update the bottom temperature and the mixed-layer temperature

         IF(h_ML_n_flk.LE.depth_w-h_ML_min_flk) THEN       ! Mixing did not reach the bottom 

            IF(h_ML_n_flk.GT.h_ML_p_flk) THEN   ! Mixed-layer deepening 
               R_H_icesnow     = h_ML_p_flk/depth_w
               R_rho_c_icesnow = 1.-R_H_icesnow 
               R_TI_icesnow    = 0.5       *C_T_p_flk*R_rho_c_icesnow+C_TT_flk*(2.*R_H_icesnow-1.)
               R_Tstar_icesnow = (0.5       +C_TT_flk-C_Q_flk)/R_TI_icesnow
               R_TI_icesnow    = (1.-C_T_p_flk*R_rho_c_icesnow)/R_TI_icesnow
               
               d_T_bot_dt = (Q_w_flk-Q_bot_flk+I_w_flk-I_bot_flk)/tpl_rho_w_r/tpl_c_w
               d_T_bot_dt = d_T_bot_dt - C_T_p_flk*(T_wML_p_flk-T_bot_p_flk)*d_h_ML_dt
               d_T_bot_dt = d_T_bot_dt*R_Tstar_icesnow/depth_w                   ! Q+I fluxes and dh_ML/dt term

               flk_str_2 = I_intm_h_D_flk - (1.-C_Q_flk)*I_h_flk - C_Q_flk*I_bot_flk
               flk_str_2 = flk_str_2*R_TI_icesnow/(depth_w-h_ML_p_flk)/tpl_rho_w_r/tpl_c_w
               d_T_bot_dt = d_T_bot_dt + flk_str_2                               ! Add radiation-flux term

               flk_str_2 = (1.-C_TT_2*R_TI_icesnow)/C_T_p_flk
               flk_str_2 = flk_str_2*(T_wML_p_flk-T_bot_p_flk)*d_C_T_dt
               d_T_bot_dt = d_T_bot_dt + flk_str_2                               ! Add dC_T/dt term
                  
            ELSE                                ! Mixed-layer retreat or stationary state
                  d_T_bot_dt = 0.                                            ! dT_bot/dt=0
            ENDIF

            T_bot_n_flk = T_bot_p_flk + d_T_bot_dt*del_time                      ! Update T_bot  
            T_bot_n_flk = MAX(T_bot_n_flk, tpl_T_f)           ! Security, limit T_bot by the freezing point
            flk_str_2 = (T_bot_n_flk-tpl_T_r)*flake_buoypar(T_mnw_n_flk)
            IF(flk_str_2.LT.0.) T_bot_n_flk = tpl_T_r  ! Security, avoid T_r crossover 
            T_wML_n_flk = C_T_n_flk*(1.-h_ML_n_flk/depth_w)
            T_wML_n_flk = (T_mnw_n_flk-T_bot_n_flk*T_wML_n_flk)/(1.-T_wML_n_flk)
            T_wML_n_flk = MAX(T_wML_n_flk, tpl_T_f)           ! Security, limit T_wML by the freezing point

         ELSE                                              ! Mixing down to the lake bottom 

            h_ML_n_flk = depth_w
            T_wML_n_flk = T_mnw_n_flk
            T_bot_n_flk = T_mnw_n_flk
            C_T_n_flk = C_T_min

         ENDIF

      ENDIF HTC_Water

!------------------------------------------------------------------------------
!  Compute the depth of the upper layer of bottom sediments
!  and the temperature at that depth.
!------------------------------------------------------------------------------

      Use_sediment: IF(lflk_botsed_use) THEN   ! The bottom-sediment scheme is used
      
         IF(H_B1_p_flk.GE.depth_bs-H_B1_min_flk) THEN   ! No T(z) maximum (no thermal wave) 
            H_B1_p_flk = 0.                       ! Set H_B1_p to zero
            T_B1_p_flk = T_bot_p_flk                     ! Set T_B1_p to the bottom temperature
         ENDIF 

         flk_str_1 = 2.*Phi_B1_pr0/(1.-C_B1)*tpl_kappa_w/tpl_rho_w_r/tpl_c_w*del_time
         h_ice_threshold = SQRT(flk_str_1)                              ! Threshold value of H_B1
         h_ice_threshold = MIN(0.9       *depth_bs, h_ice_threshold)    ! Limit H_B1
         flk_str_2 = C_B2/(1.-C_B2)*(T_bs-T_B1_p_flk)/(depth_bs-H_B1_p_flk)

         IF(H_B1_p_flk.LT.h_ice_threshold) THEN  ! USE a truncated equation for H_B1(t)
            H_B1_n_flk = SQRT(H_B1_p_flk**2+flk_str_1)  ! Advance H_B1
            d_H_B1_dt = (H_B1_n_flk-H_B1_p_flk)/del_time          ! Re-compute dH_B1/dt 
         ELSE                                    ! USE a full equation for H_B1(t)
            flk_str_1 = (Q_bot_flk+I_bot_flk)/H_B1_p_flk/tpl_rho_w_r/tpl_c_w
            flk_str_1 = flk_str_1 - (1.-C_B1)*(T_bot_n_flk-T_bot_p_flk)/del_time
            d_H_B1_dt = (1.-C_B1)*(T_bot_p_flk-T_B1_p_flk)/H_B1_p_flk + C_B1*flk_str_2
            d_H_B1_dt = flk_str_1/d_H_B1_dt
            H_B1_n_flk = H_B1_p_flk + d_H_B1_dt*del_time          ! Advance H_B1
         ENDIF 
         d_T_B1_dt = flk_str_2*d_H_B1_dt
         T_B1_n_flk = T_B1_p_flk + d_T_B1_dt*del_time            ! Advance T_B1

         !_nu  
         ! USE a very simplistic procedure, where only the upper layer profile is used, 
         ! H_B1 is always set to depth_bs, and T_B1 is always set to T_bs.
         ! THEN, the time derivatives are zero, and the sign of the bottom heat flux depends on 
         ! whether T_bot is smaller or greater than T_bs.
         ! This is, of course, an oversimplified scheme.
         !_nu  d_H_B1_dt = 0.
         !_nu  d_T_B1_dt = 0.
         !_nu  H_B1_n_flk = H_B1_p_flk + d_H_B1_dt*del_time   ! Advance H_B1
         !_nu  T_B1_n_flk = T_B1_p_flk + d_T_B1_dt*del_time   ! Advance T_B1

         l_snow_exists = H_B1_n_flk.GE.depth_bs-H_B1_min_flk                    & ! H_B1 reached depth_bs, or
                     .OR. H_B1_n_flk.LT.H_B1_min_flk                             & ! H_B1 decreased to zero, or
                     .OR.(T_bot_n_flk-T_B1_n_flk)*(T_bs-T_B1_n_flk).LE.0.   ! there is no T(z) maximum
         IF(l_snow_exists) THEN      
               H_B1_n_flk = depth_bs                     ! Set H_B1 to the depth of the thermally active layer
               T_B1_n_flk = T_bs                         ! Set T_B1 to the climatological temperature 
         ENDIF

      ELSE Use_sediment                        ! The bottom-sediment scheme is not used

         H_B1_n_flk = rflk_depth_bs_ref              ! H_B1 is set to a reference value 
         T_B1_n_flk = tpl_T_r                        ! T_B1 is set to the temperature of maximum density

      ENDIF Use_sediment

!------------------------------------------------------------------------------
!  Impose additional constraints.
!------------------------------------------------------------------------------

      ! in case of unstable stratification, force mixing down to the bottom
      flk_str_2 = (T_wML_n_flk-T_bot_n_flk)*flake_buoypar(T_mnw_n_flk)
      IF(flk_str_2.LT.0.) THEN 
         h_ML_n_flk = depth_w
         T_wML_n_flk = T_mnw_n_flk
         T_bot_n_flk = T_mnw_n_flk
         C_T_n_flk = C_T_min
      ENDIF

      !------------------------------------------------------------------------------
      !  Update the surface temperature.
      !------------------------------------------------------------------------------

      IF(h_snow_n_flk.GE.h_Snow_min_flk) THEN   
         T_sfc_n = T_snow_n_flk                   ! Snow exists, USE the snow temperature
      ELSE IF(h_ice_n_flk.GE.h_Ice_min_flk) THEN
         T_sfc_n = T_ice_n_flk                    ! Ice exists but there is no snow, USE the ice temperature
      ELSE 
         T_sfc_n = T_wML_n_flk                    ! No ice-snow cover, USE the mixed-layer temperature
      ENDIF

      !------------------------------------------------------------------------------
      ! Update the lake state variables.
      !------------------------------------------------------------------------------
      T_sfc_p      = T_sfc_n                                           
      T_snow_p_flk = T_snow_n_flk                                     
      T_ice_p_flk  = T_ice_n_flk                                      
      T_wML_p_flk  = T_wML_n_flk                                     
      T_mnw_p_flk  = T_mnw_n_flk                                    
      T_bot_p_flk  = T_bot_n_flk                                     
      T_B1_p_flk   = T_B1_n_flk                                    
      h_snow_p_flk = h_snow_n_flk                                  
      h_ice_p_flk  = h_ice_n_flk                                 
      h_ML_p_flk   = h_ML_n_flk                                   
      H_B1_p_flk   = H_B1_n_flk                              
      C_T_p_flk    = C_T_n_flk                                       

      !------------------------------------------------------------------------------
      !  END calculations
      !==============================================================================

   END SUBROUTINE flaketem



   SUBROUTINE SfcFlx_roughness ( &                           
                              ! "in" arguments
                              ! ---------------------------
                              U_a, u_star, fetch, h_ice, &
                              ! "out" arguments
                              ! ---------------------------
                              c_z0u_fetch, u_star_thresh, z0u, z0t, z0q )

! ---------------------------------- code history -------------------------------------
! Description:
!
!  Computes the water-surface or the ice-surface roughness lengths
!  with respect to wind velocity, potential temperature and specific humidity.
!
!  The water-surface roughness lengths with respect to wind velocity is computed
!  from the Charnock formula when the surface is aerodynamically rough.
!  A simple empirical formulation is used to account for the dependence 
!  of the Charnock parameter on the wind fetch. 
!  When the flow is aerodynamically smooth, the roughness length with respect to 
!  wind velocity is proportional to the depth of the viscous sub-layer.
!  The water-surface roughness lengths for scalars are computed using the power-law 
!  formulations in terms of the roughness Reynolds number (Zilitinkevich et al. 2001).
!  The ice-surface aerodynamic roughness is taken to be constant.
!  The ice-surface roughness lengths for scalars 
!  are computed through the power-law formulations 
!  in terms of the roughness Reynolds number (Andreas 2002).
!
! Original author: 
!     Dmitrii Mironov, 2005/11/17
!==============================================================================
   USE MOD_Lake_Const, only: vonkar,  grav, c_accur_sf   , Re_z0s_ice_t , tpsf_nu_u_a  ,&
                              c_z0u_smooth , z0u_ice_rough, c_z0u_ftch_f , c_z0u_ftch_ex,&
                              h_Ice_min_flk, u_wind_min_sf, c_z0t_ice_b0t, c_z0q_ice_b0t,&
                              c_z0t_ice_b0s, c_z0q_ice_b0s, c_z0t_ice_b1t, c_z0q_ice_b1t,&
                              c_z0t_ice_b1r, c_z0q_ice_b1r, c_z0t_ice_b2r, c_z0q_ice_b2r,&
                              c_z0t_ice_b0r, c_z0q_ice_b0r, c_z0u_rough  , c_z0u_rough_L,&
                              c_z0t_rough_1, c_z0t_rough_2, c_z0t_rough_3, c_z0q_rough_1,&
                              c_z0q_rough_2, c_z0q_rough_3, Pr_neutral   , Sc_neutral
!================================================================================
!  ------------------------- input variables --------------------------- 
   real(r8), intent(in)     :: &
      fetch                   ,&! Typical wind fetch [m]
      U_a                     ,&! Wind speed [m s^{-1}]
      u_star                  ,&! Friction velocity in the surface air layer [m s^{-1}]
      h_ice                     ! Ice thickness [m]

!  ------------------------- output variables ---------------------------
   real(r8), intent(out)    :: &
      c_z0u_fetch             ,&! Fetch-dependent Charnock parameter
      u_star_thresh           ,&! Threshold value of friction velocity [m s^{-1}]
      z0u                     ,&! Roughness length with respect to wind velocity [m]
      z0t                     ,&! Roughness length with respect to potential temperature [m]
      z0q                       ! Roughness length with respect to specific humidity [m]
        
!  ------------------------- local variables ---------------------------
   real(r8)                 :: &
      Re_s                    ,&! Surface Reynolds number 
      Re_s_thresh               ! Threshold value of Re_s
!================================================================================

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------
      Water_or_Ice: IF(h_ice.LT.h_Ice_min_flk) THEN  ! Water surface  
            
         ! The Charnock parameter as dependent on dimensionless fetch
         c_z0u_fetch = MAX(U_a, u_wind_min_sf)**2/grav/fetch  ! Inverse dimensionless fetch
         c_z0u_fetch = c_z0u_rough + c_z0u_ftch_f*c_z0u_fetch**c_z0u_ftch_ex
         c_z0u_fetch = MIN(c_z0u_fetch, c_z0u_rough_L)                      ! Limit Charnock parameter

         ! Threshold value of friction velocity
         u_star_thresh = (c_z0u_smooth/c_z0u_fetch*grav*tpsf_nu_u_a)**(1./3.)

         ! Surface Reynolds number and its threshold value
         Re_s = u_star**3/tpsf_nu_u_a/grav
         Re_s_thresh = c_z0u_smooth/c_z0u_fetch

         ! Aerodynamic roughness
         IF(Re_s.LE.Re_s_thresh) THEN                 
            z0u = c_z0u_smooth*tpsf_nu_u_a/u_star     ! Smooth flow
         ELSE
            z0u = c_z0u_fetch*u_star*u_star/grav  ! Rough flow
         ENDIF 

         ! Roughness for scalars  
         z0q = c_z0u_fetch*MAX(Re_s, Re_s_thresh)
         z0t = c_z0t_rough_1*z0q**c_z0t_rough_3 - c_z0t_rough_2
         z0q = c_z0q_rough_1*z0q**c_z0q_rough_3 - c_z0q_rough_2
         z0t = z0u*EXP(-vonkar/Pr_neutral*z0t)
         z0q = z0u*EXP(-vonkar/Sc_neutral*z0q) 
      
      ELSE Water_or_Ice                              ! Ice surface
        
         ! The Charnock parameter is not used over ice, formally set "c_z0u_fetch" to its minimum value
         c_z0u_fetch = c_z0u_rough

         ! Threshold value of friction velocity
         u_star_thresh = c_z0u_smooth*tpsf_nu_u_a/z0u_ice_rough

         ! Aerodynamic roughness
         z0u = MAX(z0u_ice_rough, c_z0u_smooth*tpsf_nu_u_a/u_star)

         ! Roughness Reynolds number 
         Re_s = MAX(u_star*z0u/tpsf_nu_u_a, c_accur_sf)

         ! Roughness for scalars  
         IF(Re_s.LE.Re_z0s_ice_t) THEN 
            z0t = c_z0t_ice_b0t + c_z0t_ice_b1t*LOG(Re_s)
            z0t = MIN(z0t, c_z0t_ice_b0s)
            z0q = c_z0q_ice_b0t + c_z0q_ice_b1t*LOG(Re_s)
            z0q = MIN(z0q, c_z0q_ice_b0s)
         ELSE 
            z0t = c_z0t_ice_b0r + c_z0t_ice_b1r*LOG(Re_s) + c_z0t_ice_b2r*LOG(Re_s)**2
            z0q = c_z0q_ice_b0r + c_z0q_ice_b1r*LOG(Re_s) + c_z0q_ice_b2r*LOG(Re_s)**2
         ENDIF
         z0t = z0u*EXP(z0t)
         z0q = z0u*EXP(z0q)
      ENDIF Water_or_Ice

!------------------------------------------------------------------------------
!  END calculations
!==============================================================================
   END SUBROUTINE SfcFlx_roughness



   real (r8) FUNCTION SfcFlx_satwvpres (T, h_ice)
! ---------------------------------- code history -------------------------------------
! Description:
!     Computes saturation water vapour pressure 
!     over the water surface or over the ice surface
!     as function of temperature. 
!
! Original author: 
!     Dmitrii Mironov, 2005/11/17
!==============================================================================
   USE MOD_Lake_Const, only: h_Ice_min_flk
!================================================================================
!  ------------------------- input variables ---------------------------
   real(r8), intent(in)     :: &
      T                       ,&! Temperature [K]s
      h_ice                     ! Ice thickness [m]
!  ------------------------- local variables ---------------------------
   real(r8), parameter :: b1_vap   = 610.78       ! Coefficient [N m^{-2} = kg m^{-1} s^{-2}]
   real(r8), parameter :: b3_vap   = 273.16       ! Triple point [K]
   real(r8), parameter :: b2w_vap  = 17.2693882   ! Coefficient (water)
   real(r8), parameter :: b2i_vap  = 21.8745584   ! Coefficient (ice) 
   real(r8), parameter :: b4w_vap  = 35.86        ! Coefficient (temperature) [K]
   real(r8), parameter :: b4i_vap  = 7.66         ! Coefficient (temperature) [K]
!================================================================================

   !  Start calculations
      ! Saturation water vapour pressure [N m^{-2} = kg m^{-1} s^{-2}]
      IF(h_ice.LT.h_Ice_min_flk) THEN  ! Water surface
         SfcFlx_satwvpres = b1_vap*EXP(b2w_vap*(T-b3_vap)/(T-b4w_vap))
      ELSE                             ! Ice surface
         SfcFlx_satwvpres = b1_vap*EXP(b2i_vap*(T-b3_vap)/(T-b4i_vap))
      ENDIF 
   !  END calculations
   END FUNCTION SfcFlx_satwvpres



   real (r8) FUNCTION SfcFlx_spechum (wvpres, P)
! ---------------------------------- code history -------------------------------------
! Description:
!     Computes specific humidity as function 
!     of water vapour pressure and air pressure. 
!
! Original author: 
!     Dmitrii Mironov, 2005/11/17
!
!==============================================================================
   USE MOD_Lake_Const, only: tpsf_Rd_o_Rv
!================================================================================
!  ------------------------- input variables ---------------------------
   real(r8), intent(in)     :: &
      wvpres                  ,&! Water vapour pressure [N m^{-2} = kg m^{-1} s^{-2}]
      P                         ! Air pressure [N m^{-2} = kg m^{-1} s^{-2}]
!================================================================================

      ! Specific humidity 
      SfcFlx_spechum = tpsf_Rd_o_Rv*wvpres/(P-(1.-tpsf_Rd_o_Rv)*wvpres)

   END FUNCTION SfcFlx_spechum



   real (r8) FUNCTION SfcFlx_rhoair (T, q, P)
! ---------------------------------- code history -------------------------------------
! Description:
!     Computes the air density as function 
!     of temperature, specific humidity and pressure.
!
! Original author: 
!     Dmitrii Mironov, 2005/11/17
!
!==============================================================================
   USE MOD_Lake_Const, only: tpsf_R_dryair, tpsf_Rd_o_Rv
!================================================================================
!  ------------------------- input variables ---------------------------
   real(r8), intent(in)     :: &
      T                       ,& ! Temperature [K]
      q                       ,& ! Specific humidity 
      P                          ! Pressure [N m^{-2} = kg m^{-1} s^{-2}]
!================================================================================
      
      ! Air density [kg m^{-3}] 
        SfcFlx_rhoair = P/tpsf_R_dryair/T/(1. +(1. /tpsf_Rd_o_Rv-1.)*q)

   END FUNCTION SfcFlx_rhoair



   real (r8) FUNCTION flake_snowheatconduct (h_snow)
! ---------------------------------- code history -------------------------------------
! Description:
!     Computes the snow heat conductivity,
!     using an empirical approximation from Heise et al. (2003).
!
! Original author: 
!     Dmitrii Mironov, 2005/11/17
!
!==============================================================================
   USE MOD_Lake_Const,only: tpl_rho_w_r, tpl_kappa_S_min,&
                                 tpl_kappa_S_max, tpl_Gamma_kappa_S
!================================================================================
!  ------------------------- input variables ---------------------------
   real(r8) , intent(in) :: h_snow  ! Snow thickness [m]
!================================================================================

      ! Snow heat conductivity [J m^{-1} s^{-1} K^{-1} = kg m s^{-3} K^{-1}]
      flake_snowheatconduct = flake_snowdensity( h_snow )   ! Compute snow density
      flake_snowheatconduct = MIN( tpl_kappa_S_max, tpl_kappa_S_min                      &
                              + h_snow*tpl_Gamma_kappa_S*flake_snowheatconduct/tpl_rho_w_r )

   END FUNCTION flake_snowheatconduct



   real (r8) FUNCTION flake_snowdensity (h_snow)
! ---------------------------------- code history -------------------------------------
! Description:
!     Computes the snow density,
!     using an empirical approximation from Heise et al. (2003).
!    
! Original author: 
!     Dmitrii Mironov, 2005/11/17
!
!==============================================================================
   USE MOD_Lake_Const,only: tpl_rho_w_r, tpl_rho_S_min,&
                                 tpl_rho_S_max, tpl_Gamma_rho_S, c_small_flk
!================================================================================
!  ------------------------- input variables ---------------------------
   real(r8) , intent(in) :: h_snow  ! Snow thickness [m]
!================================================================================
       
      ! Snow density [kg m^{-3}]
      !  Security. Ensure that the expression in () does not become negative at a very large h_snow.
      flake_snowdensity = MAX( c_small_flk, (1. - h_snow*tpl_Gamma_rho_S/tpl_rho_w_r) )
      flake_snowdensity = MIN( tpl_rho_S_max, tpl_rho_S_min/flake_snowdensity )

   END FUNCTION flake_snowdensity 



   real (r8) FUNCTION flake_buoypar (T_water)
! ---------------------------------- code history -------------------------------------
! Description:
!     Computes the buoyancy parameter,
!     using a quadratic equation of state for the fresh-water.
!
! Original author: 
!     Dmitrii Mironov, 2005/11/17
!
!==============================================================================
   USE MOD_Lake_Const,only: grav, tpl_T_r,tpl_a_T
!================================================================================
!  ------------------------- input variables ---------------------------
   real(r8) , intent(in) :: T_water! Water temperature [K]
!================================================================================

      ! Buoyancy parameter [m s^{-2} K^{-1}]
      flake_buoypar = grav*tpl_a_T*(T_water-tpl_T_r)

   END FUNCTION flake_buoypar


END MODULE  MOD_Lake_FLake