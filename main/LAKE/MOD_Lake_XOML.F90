#include <define.h>

MODULE MOD_Lake_XOML

   USE MOD_Precision
   IMPLICIT NONE

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: XOML_Driver
   PUBLIC :: XOML_SurfFlux

   ! PRIVATE MEMBER FUNCTIONS:
   PRIVATE :: XOML_solve
   PRIVATE :: moninobuk_xoml
   PRIVATE :: moninobukini_xoml
   PRIVATE :: roughness_lake
   PRIVATE :: radpene
   PRIVATE :: mtmpene
   PRIVATE :: solve_t
   PRIVATE :: update_rho
   PRIVATE :: update_n2_s2
   PRIVATE :: solve_uv
   PRIVATE :: solve_tke
   PRIVATE :: update_diag
   PRIVATE :: update_eddy_pars
   PRIVATE :: update_hobl
   PRIVATE :: cal_elnth
   PRIVATE :: cal_prandt
   PRIVATE :: cal_ri
   PRIVATE :: buoyeval
   PRIVATE :: trid
   PRIVATE :: nhf
   PRIVATE :: phf

!+WMEJ TODO: Add description
!+WMEJ TODO: Add history
!+WMEJ TODO: Add reference 
!+WMEJ TODO: check mtmpene SUBROUTINE and radpene SUBROUTINE, hobl may be equal to lake depth when simulating lake
!+WMEJ TODO: check <tcoef>
!+WMEJ TODO: check the signes of <lhf, shf>
!+WMEJ TODO: Check the parameterization scheme for calculating hobl. No comments were found for this SUBROUTINE <update_hobl>.
!+WMEJ TODO: Check the SUBROUTINE <update_eddy_pars> for the parameter <prtopt> and <riopt>. No comments were found for this SUBROUTINE.
!+WMEJ TODO: Dont know why USE eta in the SUBROUTINE <solve_t>.
!+WMEJ TODO: Check the SUBROUTINE <buoyeval>

!--------------------------------------------------------------------------
CONTAINS
!--------------------------------------------------------------------------

   SUBROUTINE XOML_Driver( &
               ! "in" arguments    
               ! -------------------
               scwat    , nlayer   , nsoil    , nsnow    ,&
               zopt     , betaopt  , fetchopt , etaopt   ,&
               slpeopt  , mtpeopt  , rhoopt   , prtopt   ,&
               elopt    , tskBC    , SSTNudg  ,&
               xlat     , xlon     , dtlak    , dplak    ,&
               hwref    , htref    , hqref    , tmref    ,&
               usurf    , vsurf    , qmref    , arhos    ,&
               sols     , soll     , solsd    , solld    ,&
               sabg     , lwdns    , psurf    , zcbcv    ,&
               crain    , csnow    , lrain    , lsnow    ,&
               obsdsst  , predsst  , hpbl     , prtin    ,&
               ! "inout" arguments
               ! -------------------
               dzlak    , zilak    , zlake    , lktmp    ,&
               uwatv    , vwatv    , lksal    , lkrho   ,&
               dzssb    , zissb    , zssb     , t_ssb    ,&
               icefr    , xwliq    , xwice    , snlay    ,&
               scv      , snwag    , snwdp    , rhosnw   ,&
               icedp    , bicedp   , wicedp   , tskin    ,&
               tke      , eps      , hobl     , stke1    ,&
               Km       , Kh       , Ke       , tmice    ,&
               z0m      , z0h      , z0q      , frlak    ,&
               felak    , gamma    , etal     , btpri    ,&
               ! "out" arguments
               ! -------------------
               fsena    , fevpa    , lfevpa   , fevpg    ,&
               fseng    , olrg     , fgrnd    , trad     ,&
               taux     , tauy     , ustar    , qstar    ,&
               tstar    , emis     , zol      , snwml    ,&
               rib      , ram      , rah      , raw      ,&
               wdm      , t2m      , q2m      , u10m     ,&
               v10m     , fm10m    , fq10m    , fh2m     ,&
               fq2m     , fm       , fq       , fh       ,&
               shfdt    , urban_call)

! ---------------------------------- code history -------------------------------------
! Description:
!     Interface of XOML. Used to initialize and calculate model variables.
!     XOML is a multi-purpose model that can be used for both lakes and oceans.
!     It is based on the Ocean Mixed Layer Model (OMLM, Noh et al., (1996, 1999, 2011)).
!     reference: 
!         Noh, Y. (1996), Dynamics of Diurnal Thermocline Formation in the Oceanic Mixed Layer. https://doi.org/10.1175/1520-0485(1996)026<2183:DODTFI>2.0.CO;2
!         Noh, Y. et al. (1999), Simulations of temperature and turbulence structure of the oceanic boundary layer with the improved near‐surface process. https://doi.org/10.1029/1999JC900068
!         Noh, Y. et al. (2002), Simulation of More Realistic Upper-Ocean Processes from an OGCM with a New Ocean Mixed Layer Model.  https://doi.org/10.1175/1520-0485(2002)032<1284:SOMRUO>2.0.CO;2
!         Noh, Y. et al. (2005), Effect of the Prandtl number in the parameterization of vertical mixing in an OGCM of the tropical Pacific. https://doi.org/10.1029/2005GL024540
!         Noh, Y. et al. (2011a), Influence of Langmuir Circulation on the Deepening of the Wind-Mixed Layer. https://doi.org/10.1175/2010JPO4494.1
!         Noh, Y. et al. (2011b), Prediction of the diurnal warming of sea surface temperature using an atmosphere‐ocean mixed layer coupled model. https://doi.org/10.1029/2011JC006970
!         Ling, T. et al. (2015), A multilevel ocean mixed layer model resolving the diurnal cycle: Development and validation. https://doi.org/10.1002/2015MS000476
!         Sun, L. et al. (2020), Improving a Multilevel Turbulence Closure Model for a Shallow Lake in Comparison With Other 1‐D Models. https://doi.org/10.1029/2019MS001971
!
! Called: (* means optional)
!    *-> XOML_SurfFlux     : Calculate the surface fluxes (cssp Shecme)
!    *-> CoLML_SurfFlux      : Calculate the surface fluxes (CoLM Shecme)
!    -> XOML_solve         : Solve the lake model
!
! Original author: 
!     Tiejun Ling
!     
! Revisions:
!     Xin-Zhong Liang
!     Lei Sun, 2019: Modified for lake
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024: Reorganized the code to adapt to the new LakeDriver
!
! -------------------------------------------------------------------------------------
   USE MOD_Lake_Const, only: PI
   USE MOD_Lake_CoLML, only: CoLML_SurfFlux
   USE MOD_Namelist, only: DEF_USE_COLML_FLUX_SCHEME
!================================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in)      :: &
      scwat                   ,&! surface category of water characteristics:
      nlayer                  ,&! Maximum number of layer (lake or ocean)
      nsoil                   ,&! Maximum number of soil layer
      nsnow                   ,&! Maximum number of snow layer
      zopt                    ,&! option for roughness length, 1: constant, 2: Subin et al. (2012)
      betaopt                 ,&! option for betaprime calculation, 1: constant, 2: equation
      fetchopt                ,&! option for fetch length, 1: constant, 2: equation
      etaopt                  ,&! option for Extinction coefficient calculation, 1: constant, 2: equation
      slpeopt                 ,&! method of solar radiation penetration, penetration is defined in the levels
      mtpeopt                 ,&! method of momentum penetration, penetration is defined in the levels
      rhoopt                  ,&! option for water density calculation
      prtopt                  ,&! option for turbulent Prandtl number calculation
      elopt                     ! option for turbulence energy length scale calculation

   logical, intent(in)      :: &
      tskBC                   ,&! option for lake surface temperature boundary condition, T=> USE SST boundary condition
      SSTNudg                   ! surface temperature nudging option, T=> USE SST Nudging

   real(r8), intent(in)     :: &
      xlat                    ,&! latitude in degrees
      xlon                    ,&! longitude in degrees
      dtlak                   ,&! time step [s]
      dplak                   ,&! lake depth [m]
      prtin                   ,&! User-defined turbulent Prandtl number
      hwref                   ,&! observational height of wind [m]
      htref                   ,&! observational height of temperature [m]
      hqref                   ,&! observational height of humidity [m]
      tmref                   ,&! temperature at the reference height [K] 
      usurf                   ,&! wind component in eastward direction [m/s]
      vsurf                   ,&! wind component in northward direction [m/s]
      qmref                   ,&! specific humidity at the reference height [kg/kg]
      arhos                   ,&! surface air density [kg m-3]
      sols                    ,&! atm vis direct beam solar rad onto srf [W m-2]
      soll                    ,&! atm nir direct beam solar rad onto srf [W m-2]
      solsd                   ,&! atm vis diffuse solar rad onto srf [W m-2]
      solld                   ,&! atm nir diffuse solar rad onto srf [W m-2]
      sabg                    ,&! solar absorbed by ground  [W/m2]
      lwdns                   ,&! atmospheric infrared (longwave) radiation [W/m2]
      psurf                   ,&! surface pressure [Pa]
      hpbl                    ,&! atmospheric boundary layer height, only USE in CoLML_SurFlux [m]
      zcbcv                   ,&! convective boundary height [m]
      obsdsst                 ,&! observed sea surface temperature [K]
      predsst                 ,&! predicted sea surface temperature [K]
      crain                   ,&! convective rainfall [kg/(m2 s)]
      csnow                   ,&! convective snowfall [kg/(m2 s)]
      lrain                   ,&! large scale rainfall [kg/(m2 s)]
      lsnow                     ! large scale snowfall [kg/(m2 s)]

!  ------------------------- inout variables ---------------------------
   integer, intent(inout)   :: &
      snlay                     ! number of snow layers

   real(r8), intent(inout)  :: &
      dzlak(nlayer)           ,&! lake layer thickness (m)
      zlake(nlayer)           ,&! lake layer depth (m)
      zilak(nlayer+1)         ,&! lake layer interface level [m]
      lktmp(nlayer)           ,&! lake temperature [K]
      uwatv(nlayer)           ,&! Water velocity in x-direction [m/s]
      vwatv(nlayer)           ,&! Water velocity in y-direction [m/s]
      lksal(nlayer)           ,&! lksality [‰]
      lkrho(nlayer)           ,&! Water density [kg/m^3]
      Km(nlayer)              ,&! eddy viscosity [m2s-1]
      Kh(nlayer)              ,&! diffusivity for scalar [m2/s]
      Ke(nlayer)              ,&! diffusivity for tke [m2/s]
      tke(nlayer+1)           ,&! Turbulent kinetic energy (TKE) [J/kg]
      eps(nlayer+1)           ,&! TKE dissipation rate [W/kg]
      icefr(nlayer)             ! lake ice fraction [-]
   
   real(r8), intent(inout)  :: &
      dzssb(-nsnow+1:nsoil)   ,&! soil + snow layer thickness [m]
      zssb(-nsnow+1:nsoil)    ,&! soil + snow layer depth [m]
      zissb(-nsnow:nsoil+1)   ,&! soil + snow layer interface level [m]
      t_ssb(-nsnow+1:nsoil)   ,&! soil + snow layer temperature [K]
      xwliq(-nsnow+1:nsoil)   ,&! liquid water (kg/m2)
      xwice(-nsnow+1:nsoil)     ! ice lens (kg/m2)
   
   real(r8), intent(inout)  :: &
      scv                     ,&! snow mass (kg/m2)
      hobl                    ,&! oceanic boundary layer depth [m]
      icedp                   ,&! lake ice thickness (m)
      bicedp                  ,&! black ice depth (m)
      wicedp                  ,&! white ice depth (m)
      snwag                   ,&! non dimensional snow age [-]
      snwdp                   ,&! snow depth (m)
      rhosnw                  ,&! snow density (kg/m3)
      stke1                   ,&! top level eddy conductivity [W/m/K]
      tskin                   ,&! ground surface temperature [k]
      tmice                   ,&! ice temperature [K]
      z0m                     ,&! roughness length over ground, momentum [m]
      z0h                     ,&! roughness length over ground, sensible heat [m]
      z0q                     ,&! roughness length over ground, latent heat [m]
      felak                   ,&! lake fetch (>0), otherwise USE emperical eq by dplak [m]
      gamma                   ,&! Mixing enhancement factorfor eddy diffusion coefficient
      etal                    ,&! extinction coefficient [1/m]
      btpri                   ,&! beta prime in Monin-Obukhov theory
      frlak                     ! lake fraction [-]

!  ------------------------- output variables ---------------------------
   real(r8), intent(out)    :: &
      fsena                   ,&! sensible heat from canopy height to atmosphere [W/m2]
      fevpa                   ,&! evapotranspiration from canopy height to atmosphere [mm/s]
      lfevpa                  ,&! latent heat flux from canopy height to atmosphere [W/2]
      fevpg                   ,&! evaporation heat flux from ground [mm/s]
      fseng                   ,&! sensible heat flux from ground [W/m2]
      olrg                    ,&! outgoing long-wave radiation from ground+canopy
      fgrnd                   ,&! total net heat flux from the atmosphere absorbed by lake [W/m2]
      trad                    ,&! radiative temperature [K]
      taux                    ,&! wind stress: E-W [kg/m/s**2]
      tauy                    ,&! wind stress: N-S [kg/m/s**2]
      ustar                   ,&! uwatv* in similarity theory [m/s]
      qstar                   ,&! q* in similarity theory [kg/kg]
      tstar                   ,&! lktmp* in similarity theory [K]
      emis                    ,&! averaged bulk surface emissivity
      zol                     ,&! dimensionless height (z_soisno/L) used in Monin-Obukhov theory
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
      shfdt                     ! derivative of srf sen+rlat  heat flux wrt srf temp [W/m2/K]
      
   !+WMEJ urban_call is optional argument for SUBROUTINE laketem, it is not used in this version
   logical, optional, intent(in) :: urban_call   ! whether it is a urban CALL
      
!  ------------------------- local variables ---------------------------
   real(r8)                 :: &
      rhosnow                 ,&! snow density [kg/m3]
      rlaxtau                 ,&! relaxation time scale for lake temperature [s]
      fgrnd1                  ,&! net heat flux from the atmosphere absorbed by top layer [W/m2]
      rlat                    ,&! latitude in radians
      b0                      ,&! Surface Bouyancy flux
      htvp                    ,&! latent heat flux parameter
      u2m                     ,&! 2 m height wind speed in eastward direction [m/s]
      buo(nlayer)             ,&! deltarho/lkrho mean buoyancy [ms-2]
      n2(nlayer)              ,&! brunt-vaisala frequency (/s**2)
      s2(nlayer)              ,&! shear number
      Rt(nlayer)              ,&! turbulent Richard number
      p_tke(nlayer)           ,&! TKE production by mean velocity shear [m2s-3]
      b_tke(nlayer)           ,&! TKE production by bouyancy flux [m2s-3]
      d_tke(nlayer)             ! TKE diffusion term [m2s-3]

   integer :: lb, j, i, snl
!================================================================================
      !- degree to radian
      rlat = xlat * PI / 180.0

      !------------------------------------------------------------------------------
      !  Compute lake surface fluxes
      !------------------------------------------------------------------------------
      IF (DEF_USE_COLML_FLUX_SCHEME) THEN
         lb = snlay + 1

         CALL CoLML_SurfFlux ( &
               ! "in" arguments
               ! -------------------
               patchtype = scwat      , snl = snlay         , zopt = zopt             , betaopt = betaopt     ,&
               fetchopt = fetchopt    , rlat = rlat         , lakedepth = dplak       , deltim = dtlak        ,& 
               savedtke1 = stke1      , forc_hgt_u = hwref  , forc_hgt_t = htref      , forc_hgt_q = hqref    ,&
               forc_us = usurf        , forc_vs = vsurf     , forc_t = tmref          , forc_q = qmref        ,&
               forc_rhoair = arhos    , forc_psrf = psurf   , forc_sols = sols        , forc_soll = soll      ,&
               forc_solsd = solsd     , forc_solld = solld  , forc_frl = lwdns        , sabg = sabg           ,&
               hpbl = hpbl            , dzlaktop = dzlak(1) , zlakebot = zlake(nlayer), dzssbtop = dzssb(lb)  ,&
               lktmptop = lktmp(1)    , tssbtop = t_ssb(lb) , xwliqtop = xwliq(lb)    , xwicetop = xwice(lb)  ,&
               icefrtop = icefr(1)    , zcb = zcbcv         ,&
               ! "inout" arguments
               ! -------------------
               t_grnd = tskin         , z0m  = z0m          , z0h =z0h                , z0q = z0q             ,&
               felak = felak          , btpri = btpri       ,&
               ! "out" arguments
               ! -------------------
               fseng = fseng          , fevpg = fevpg       , fsena = fsena           , fevpa = fevpa         ,&
               lfevpa = lfevpa        , olrg = olrg         , fgrnd = fgrnd           , fgrnd1 = fgrnd1       ,&
               zol = zol              , rib = rib           , trad = trad             , htvp = htvp           ,&
               emis = emis            , wdm = wdm           , ram = ram               , rah = rah             ,&
               raw = raw              , shfdt = shfdt       , taux = taux             , tauy = tauy           ,&
               tref = t2m             , qref = q2m          , u10m = u10m             , v10m = v10m           ,&
               fh2m = fh2m            , fq2m = fq2m         , fm10m = fm10m           , fq10m = fq10m         ,&
               fm = fm                , fh = fh             , fq = fq                 , ustar = ustar         ,&
               qstar = qstar          , tstar = tstar       , rhosnow = rhosnow       , u2m = u2m             )
       
      ELSE
         CALL XOML_SurfFlux( &
               ! "in" arguments
               ! -------------------
               scwat = scwat          , betaopt = betaopt   , fetchopt = fetchopt     , zopt = zopt           ,&
               rlat = rlat            , hwref = hwref       , htref = htref           , hqref = hqref         ,&
               tmref = tmref          , usurf = usurf       , vsurf = vsurf           , qmref = qmref         ,&
               psurf = psurf          , arhos = arhos       , sols = sols             , soll = soll           ,&
               solsd = solsd          , solld = solld       , sabg = sabg             , lwdns = lwdns         ,&
               dplak = dplak          , zcb = zcbcv         ,&
               ! "inout" arguments
               ! -------------------
               t_grnd = tskin         , z0m = z0m           , z0h = z0h               , z0q = z0q             ,&
               felak = felak          , btpri = btpri       ,&
               ! "out" arguments
               ! -------------------
               fseng = fseng          , fsena = fsena       , fevpa = fevpa           , fevpg = fevpg         ,&
               lfevpa = lfevpa        , olrg = olrg         , fgrnd = fgrnd           , fgrnd1 = fgrnd1       ,&
               zol = zol              , rib = rib           , trad = trad             , htvp = htvp           ,&
               emis = emis            , wdm = wdm           , ram = ram               , rah = rah             ,&
               raw = raw              , shfdt = shfdt       , taux = taux             , tauy = tauy           ,&
               tref = t2m             , qref = q2m          , u10m = u10m             , v10m = v10m           ,&
               fh2m = fh2m            , fq2m = fq2m         , fm10m = fm10m           , fq10m = fq10m         ,&
               fm = fm                , fh = fh             , fq = fq                 , ustar = ustar         ,&
               qstar = qstar          , tstar = tstar       , u2m = u2m               )

      ENDIF

         rlaxtau = 8640.   !-WMEJ why 8640?
         CALL XOML_solve( &
               ! "in" arguments
               ! -------------------
               etaopt = etaopt        , slpeopt = slpeopt   , mtpeopt = mtpeopt       , rhoopt = rhoopt       ,&
               prtopt = prtopt        , elopt = elopt       , nlayer = nlayer         , tskBC = tskBC         ,&
               SSTNudg = SSTNudg      , scwat = scwat       , betaopt = betaopt       , rlat = rlat           ,&
               dtlak = dtlak          , btpri = btpri       , etal = etal             , fgrnd1 = fgrnd1       ,&
               fgrnd = fgrnd          , sabg = sabg         , taux = taux             , tauy = tauy           ,&
               tskin = tskin          , rlaxtau = rlaxtau   , obsdsst = obsdsst       , predsst = predsst     ,&
               dplak = dplak          , z0m = z0m           , prtin = prtin           ,&
               ! "inout" arguments
               ! -------------------
               dzlak = dzlak          , zlake = zlake       , zilak = zilak           , lktmp = lktmp         ,&
               uwatv = uwatv          , vwatv = vwatv       , lksal = lksal           , lkrho = lkrho       ,&
               tke = tke(1:nlayer)    , eps = eps(1:nlayer) , Km = Km                 , Kh = Kh               ,&
               Ke = Ke                , hobl = hobl         ,&
               ! "out" arguments
               ! -------------------
               b0 = b0                , buo = buo           , n2 = n2                 , s2 = s2               ,&
               Rt = Rt                , p_tke = p_tke       , b_tke = b_tke           , d_tke = d_tke          )

      scv = 0.0
      snwml = 0.0
      snwdp = 0.0
      rhosnw = 0.0
      icedp = 0.0
      bicedp = 0.0
      wicedp = 0.0
      snwag = 0.0

   END SUBROUTINE XOML_Driver



   SUBROUTINE XOML_SurfFlux( &
               ! "in" arguments
               ! -------------------
               scwat    , betaopt  , fetchopt , zopt     ,&
               rlat     , hwref    , htref    , hqref    ,&
               tmref    , usurf    , vsurf    , qmref    ,&
               psurf    , arhos    , sols     , soll     ,&
               solsd    , solld    , sabg     , lwdns    ,&
               dplak    , zcb      ,&
               ! "inout" arguments
               ! -------------------
               t_grnd   , z0m      , z0h      , z0q      ,&
               felak    , btpri    ,&
               ! "out" arguments
               ! -------------------
               fseng    , fsena    , fevpa    , fevpg    ,&
               lfevpa   , olrg     , fgrnd    , fgrnd1   ,&
               zol      , rib      , trad     , htvp     ,&
               emis     , wdm      , ram      , rah      ,&
               raw      , shfdt    , taux     , tauy     ,&
               tref     , qref     , u10m     , v10m     ,&
               fh2m     , fq2m     , fm10m    , fq10m    ,&
               fm       , fh       , fq       , ustar    ,&
               qstar    , tstar    , u2m      )

! ---------------------------------- code history -------------------------------------
! Description:
!     Main control to air-sea exchanges of heat, mass, and momentum.
!     Compute sensible and latent fluxes and their derivatives with respect
!     to ground temperature using ground temperatures from previous time step.
!
!
! Original author:
!     Yongjiu Dai, 09/15/1999; 08/30/2002
!
! Revisions:
!     Xin-Zhong Liang, July 1, 2003
!     Lei Sun, 2019: Modified for lake
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024: Reorganized the code to adapt to the new LakeDriver
!
!--------------------------------------------------------------------------------------
   USE module_sf_lake_subs, only : qsadv
   USE module_sf_lake_table, only: cpair, rair, vonkar, grav, z0sice, one3rd, &
                                    beta_wc, emisw, emisi, hvap, hsub, stefnc, &
                                    SEAICE, OCEAN 
!================================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in)      :: &
      scwat                   ,&! surface category of water characteristics: 0: ocean, 1: shallow lake, 2: deep lake
      betaopt                 ,&! option for betaprime calculation, 1: constant, 2: equation
      fetchopt                ,&! option for fetch length, 1: constant, 2: equation
      zopt                      ! option for roughness length, 1: constant, 2: Subin et al. (2012)

   real(r8), intent(in)     :: &
      rlat                    ,&! latitude in degrees
      hwref                   ,&! observational height of wind [m]
      htref                   ,&! observational height of temperature [m]
      hqref                   ,&! observational height of humidity [m]
      tmref                   ,&! temperature at the reference height [K]
      usurf                   ,&! wind component in eastward direction [m/s]
      vsurf                   ,&! wind component in northward direction [m/s]
      qmref                   ,&! specific humidity at the reference height [kg/kg]
      psurf                   ,&! surface pressure [Pa]
      arhos                   ,&! surface air density [kg m-3]
      sols                    ,&! atm vis direct beam solar rad onto srf [W m-2]
      soll                    ,&! atm nir direct beam solar rad onto srf [W m-2]
      solsd                   ,&! atm vis diffuse solar rad onto srf [W m-2]
      solld                   ,&! atm nir diffuse solar rad onto srf [W m-2]
      sabg                    ,&! solar absorbed by ground [W/m2]
      lwdns                   ,&! atmospheric infrared (longwave) radiation [W/m2]
      dplak                   ,&! lake depth [m]
      zcb                       ! convective boundary height [m]

!  ------------------------- inout variables ---------------------------
   real(r8), intent(inout)  :: &
      t_grnd                  ,&! surface temperature (kelvin)
      z0m                     ,&! roughness length over ground, momentum [m]
      z0h                     ,&! roughness length over ground, sensible heat [m]
      z0q                     ,&! roughness length over ground, latent heat [m]
      felak                   ,&! fetch length of lake [m]
      btpri                     ! fraction of solar radiation in the NIR, only used when MODWMEJ

!  ------------------------- output variables ---------------------------
   real(r8), intent(out)    :: &
      emis                    ,&! averaged bulk surface emissivity
      wdm                     ,&! lowest level wind speed including the stablity effect [m/s]
      ram                     ,&! aerodynamical resistance [s/m]
      rah                     ,&! thermal resistance [s/m]
      raw                     ,&! moisture resistance [s/m]
      shfdt                   ,&! derivative of srf sen+lat  heat flux wrt srf temp [W/m2/K]
      taux                    ,&! wind stress: E-W [kg/m/s**2]
      tauy                    ,&! wind stress: N-S [kg/m/s**2]
      tref                    ,&! 2 m height air temperature [K]
      qref                    ,&! 2 m height air specific humidity
      u10m                    ,&! 10 m height wind speed in eastward direction [m/s]
      v10m                    ,&! 10 m height wind speed in northward direction [m/s]
      u2m                     ,&! 2 m height wind speed in eastward direction [m/s]
      zol                     ,&! dimensionless height (z/L) used in Monin-Obukhov theory
      rib                     ,&! bulk Richardson number in surface layer
      fm                      ,&! integral of profile function for momentum
      fq                      ,&! integral of profile function for moisture
      fh                      ,&! integral of profile function for heat
      fh2m                    ,&! integral of profile function for heat
      fq2m                    ,&! integral of profile function for moisture
      fm10m                   ,&! integral of profile function for momentum
      fq10m                   ,&! integral of profile function for moisture
      ustar                   ,&! friction velocity [m/s]
      tstar                   ,&! temperature scaling parameter
      qstar                   ,&! moisture scaling parameter
      fseng                   ,&! sensible heat flux from ground [W/m2]
      fsena                   ,&! sensible heat from the reference height to atmosphere [W/m2]
      fevpa                   ,&! evaporation from the reference height to atmosphere [mm/s]
      fevpg                   ,&! evaporation heat flux from ground [mm/s]
      htvp                    ,&! latent heat of vapor of water (or sublimation) [J/kg]
      lfevpa                  ,&! laten heat from reference height to atmosphere [W/m2]
      olrg                    ,&! longwave up flux at surface [W/m2]
      fgrnd                   ,&! ground heat flux [W/m2]
      fgrnd1                  ,&! net heat flux from the atmosphere absorbed by top layer [W/m2]
      trad                      ! radiative temperature [K]

!  ------------------------- local variables ---------------------------
   integer                  :: &
      nmozsgn                 ,&! number of times moz changes sign
      niters                  ,&! maximum number of iterations for surface temperature
      iter                    ,&! iteration index
      ii                      ,&! loop index
      i                         ! loop index

   real(r8)                 :: &
      displax                 ,&! zero-displacement height [m]
      dth                     ,&! diff of virtual temp. between ref. height and surface
      dqh                     ,&! diff of humidity between ref. height and surface
      dthv                    ,&! diff of vir. poten. temp. between ref. height and surface
      obu                     ,&! monin-obukhov length [m]
      obuold                  ,&! monin-obukhov length from previous iteration
      degdT                   ,&! d(eg)/dT
      eg                      ,&! water vapor pressure at temperature T [pa]
      qsatg                   ,&! ground saturated specific humidity [kg/kg]
      qsatgdT                 ,&! d(qsatg)/dT
      fetch                   ,&! lake fetch (m)
      raih                    ,&! temporary variable [kg/m2/s]
      raiw                    ,&! temporary variable [kg/m2/s]
      shfdtl                  ,&! derivative of srf sensible heat flux wrt srf temp [W/m2/K]
      shfdts                  ,&! derivative of srf latent   heat flux wrt srf temp [W/m2/K]
      temp1                   ,&! relation for potential temperature profile
      temp2                   ,&! relation for specific humidity profile
      temp12m                 ,&! relation for temperature at 2m
      temp22m                 ,&! relation for specific humidity at 2m
      betaprime               ,&! beta prime in Monin-Obukhov theory
      betavis                 ,&! middle variable for betaprime calculation
      thm                     ,&! intermediate variable (tmref+0.0098*ht)
      th                      ,&! potential temperature [K]
      thv                     ,&! virtual potential temperature [K]
      thvstar                 ,&! virtual potential temperature scaling parameter
      um                      ,&! wind speed including the stablity effect [m/s]
      ur                      ,&! wind speed at reference height [m/s]
      visa                    ,&! kinematic viscosity of dry air [m2/s]
      wc                      ,&! convective velocity [m/s]
      xt                      ,&!
      xq                      ,&!
      ws                      ,&! 
      ks                      ,&! 
      zldis                   ,&! reference height minus zero-displacement height [m]
      z0mg                    ,&! roughness length over ground, momentum [m]
      z0hg                    ,&! roughness length over ground, sensible heat [m]
      z0qg                      ! roughness length over ground, latent heat [m]
!================================================================================

      ! Initialization variables
      nmozsgn = 0
      obuold = 0.

      !-WMEJ qsatg = sathdi(t_grnd, psurf, qsatgdT)
      CALL qsadv(t_grnd, psurf, eg, degdT, qsatg, qsatgdT)    !+WEMJ Unified the calculation method

      ! potential temperatur at the reference height
      thm = tmref + 0.0098 * htref                            ! intermediate variable equivalent to
                                                               ! tmref*(pgcm/psurf)**(rair/cpair)
      th = tmref * (100000. / psurf)**(rair / cpair)           ! potential T
      thv = th * (1. + 0.61 * qmref)                          ! virtual potential T
      ur = max(0.1, sqrt(usurf * usurf + vsurf * vsurf))      ! limit set to 0.1

      dth = thm - t_grnd
      dqh = qmref - qsatg
      dthv = dth * (1. + 0.61 * qmref) - 0.61 * th * dqh
      zldis = hwref - 0.

      ! Base on lake depth, assuming that small lakes are likely to be shallower
      ! Estimate crudely based on lake depth
      !-WMEJ USE scwat to control lake classification, 
      !      4 = shallow lake, 5 = deep lake, 2024-03-23
      IF (fetchopt == 1) THEN
         fetch = felak
      else
         IF (dplak < 4.) THEN
               fetch = 100. ! shallow lake
         else
               fetch = 25.*dplak ! deep lake
         ENDIF
      ENDIF 

      !-WMEJ calculate surface roughness lengths
      IF (scwat == SEAICE) THEN ! sea ice
         htvp = hsub
         emis = emisi

         z0mg  = z0sice
         z0qg  = z0mg
         z0hg  = z0mg
      ELSE IF (scwat == OCEAN) THEN ! ocean
         htvp = hvap
         emis = emisw
         visa  = 1.326e-5*(1. + (6.542e-3 + (8.301e-6 - 4.84e-9*tmref)*tmref)*tmref)
               ! - Andreas (1989) CRREL Rep. 89-11
         ! loop to obtain initial and good ustar and z0
         ustar = 0.06
         wc = 0.5
         IF (dthv .ge. 0.) THEN
               um = ur
         else
               um =  sqrt(ur * ur + wc * wc)
         ENDIF
         IF (zopt == 1) THEN
            z0mg = z0m
            z0hg = z0h
            z0qg = z0q
            ustar = vonkar*um/log(zldis/z0mg)
         else
            DO i = 1,5
               z0mg = 0.013 * ustar * ustar / grav + 0.11 * visa / ustar
               ustar = vonkar * um / log(zldis / z0mg)
            ENDDO
            CALL roughness_lake (0, t_grnd, psurf, ur, ustar, fetch, dplak, &
                                 z0mg,z0hg,z0qg)
         ENDIF 
      ENDIF

      CALL moninobukini_xoml(ur,th,thm,thv,dth,dqh,dthv,zldis,z0mg,um,obu)

      ! Evaluated stability-dependent variables using moz from prior iteration
      niters  = 10
      displax = 0.

      !----------------------------------------------------------------------
      DO iter = 1, niters              ! begin stability iteration
      !----------------------------------------------------------------------
         IF (scwat == OCEAN) THEN         ! ocean
            z0mg = 0.013 * ustar * ustar / grav + 0.11 * visa / ustar
            xq   = 2.67*(ustar*z0mg/visa)**0.25 - 2.57
            xt   = xq
            z0qg = z0mg/exp(xq)
            z0hg = z0mg/exp(xt)

         elseif (scwat == SEAICE) THEN   ! sea ice 
            htvp = hsub
            emis = emisi

            z0mg  = z0sice
            z0qg  = z0mg
            z0hg  = z0mg
         else
            IF (zopt .eq. 0) THEN
               z0mg = z0m
               z0hg = z0h
               z0qg = z0q
            ELSE IF (zopt .eq. 1) THEN
               CALL roughness_lake(0, t_grnd, psurf, ur, ustar, fetch, dplak, &
                                    z0mg,z0hg,z0qg)
            ENDIF
         ENDIF
            ustar = vonkar * um / log(zldis / z0mg)

            CALL moninobuk_xoml(hwref, htref, hqref, displax, z0mg, z0hg, z0qg, obu, um,&
                                ustar, temp1, temp2, temp12m, temp22m, fq10m, fm10m, fm, fh, fq)

            tstar = temp1 * dth
            qstar = temp2 * dqh

            thvstar = tstar + 0.61 * th * qstar
            zol = zldis * vonkar * grav * thvstar / (ustar**2 * thv)

            IF (zol >= 0.) THEN     !stable
               zol = min(2., max(zol,1.e-6))
            else                    !unstable
               zol = max(-100., min(zol,-1.e-6))
            ENDIF

            obu = zldis / zol
        
            IF (zol >= 0.) THEN
               um = ur
            else
               wc = beta_wc * (-grav * ustar * thvstar * zcb / thv)**one3rd
               um = sqrt(ur * ur + wc * wc)
            ENDIF

            IF (obuold * obu < 0.) nmozsgn = nmozsgn+1
            IF (nmozsgn >= 4) exit

            obuold = obu
         !----------------------------------------------------------------------
         ENDDO                       ! end stability iteration
         !----------------------------------------------------------------------
         
         ! Get derivative of fluxes with repect to ground temperature
         ram    = 1./(ustar*ustar/um)
         rah    = 1./(temp1*ustar)
         raw    = 1./(temp2*ustar)

         raih   = arhos * cpair / rah
         raiw   = arhos / raw
         shfdts = raih
         shfdtl = raiw * qsatgdT
         shfdt  = shfdts + htvp * shfdtl

         ! Get surface fluxes of momentum, sensible, latent and outgoing longwave radiation
         !     using ground temperatures from previous time step
         taux   = -arhos * usurf / ram
         tauy   = -arhos * vsurf / ram

         fseng  = -raih * dth
         fevpg  = -raiw * dqh
         !      print *, "dqh=",dqh
         fsena  = fseng
         fevpa  = fevpg
         lfevpa = htvp * fevpa

         ! Get net longwave from ground to atmosphere
         olrg   = stefnc * emis * t_grnd**4 + (1.-emis) * lwdns

         ! Get net ground heat flux
         fgrnd  = sabg + lwdns - olrg - fsena - lfevpa

         ! Get 2 m height air temperature and specific humidity
         tref   = thm + temp1*dth * (1./temp12m - 1./temp1)
         qref   = qmref + temp2*dqh * (1./temp22m - 1./temp2)
         rib    = zol*vonkar**3*ustar**2/(temp1*um**2)

         ! 10 m wind
         u10m   = usurf/ur * ustar/vonkar * fm10m
         v10m   = vsurf/ur * ustar/vonkar * fm10m
         u2m    = max(0.1_r8, ustar / vonkar * log(2._r8 / z0mg))

         !+WMEJ Recalculation of supplementary variables
         IF (betaopt == 1) THEN
            betaprime = btpri
         else
            betaprime = (soll + solld) / max(1.e-5, sols + soll + solsd + solld)
            betavis = 0. ! The fraction of the visible (e.g. vis not nir from atm) sunlight
                        ! absorbed in ~1 m of water (the surface layer za_lake).
                        ! This is roughly the fraction over 700 nm but may depend on the details
                        ! of atmospheric radiative transfer.
                        ! As long as NIR = 700 nm and up, this can be zero.
            betaprime = betaprime + (1.0 - betaprime) * betavis
         ENDIF

         fgrnd1 = betaprime * sabg + lwdns - olrg - fsena - lfevpa
         fgrnd  = sabg + lwdns - olrg - fsena - lfevpa
         ws = 1.2e-03_r8 * u2m
         ks = 6.6_r8 * sqrt(abs(sin(rlat)))*(u2m**(-1.84_r8))
         z0m    = z0mg
         z0h    = z0hg
         z0q    = z0qg
         felak  = fetch
         btpri  = betaprime
         emis   = emis
         wdm    = um
         fh     = vonkar / temp1
         fq     = vonkar / temp2
         fh2m   = vonkar / temp12m
         fq2m   = vonkar / temp22m
         fm10m  = fm10m
         ram    = 1. / (ustar * ustar / um)
         rah    = 1. / (temp1 * ustar)
         raw    = 1. / (temp2 * ustar)
         raih   = arhos * cpair / rah
         raiw   = arhos / raw
         shfdts = raih
         shfdtl = raiw * qsatgdT
         shfdt  = shfdts + htvp * shfdtl
         taux   = -arhos * usurf / ram
         tauy   = -arhos * vsurf / ram
         tref   = thm + temp1*dth * (1./temp12m - 1./temp1)
         qref   = qmref + temp2*dqh * (1./temp22m - 1./temp2)
         u10m   = usurf/ur * ustar/vonkar * fm10m
         v10m   = vsurf/ur * ustar/vonkar * fm10m
         u2m    = max(0.1_r8, ustar / vonkar * log(2._r8 / z0mg))
         zol    = zol
         rib    = zol*vonkar**3*ustar**2/(temp1*um**2)
        
   END SUBROUTINE XOML_SurfFlux



   SUBROUTINE XOML_solve( &
               ! "in" arguments
               ! -------------------
               etaopt    , slpeopt , mtpeopt , rhoopt ,&
               prtopt    , elopt   , nlayer  , tskBC  ,&
               SSTNudg   , scwat   , betaopt , rlat   ,&
               dtlak     , btpri   , etal    , fgrnd1 ,&
               fgrnd     , sabg    , taux    , tauy   ,&
               tskin     , rlaxtau , obsdsst , predsst,&
               dplak     , z0m     , prtin   ,&
               ! "inout" arguments
               ! -------------------
               dzlak     , zlake   , zilak   , lktmp  ,&
               uwatv     , vwatv   , lksal   , lkrho ,&
               tke       , eps     , Km      , Kh     ,&
               Ke        , hobl   ,&
               ! "out" arguments
               ! -------------------
               b0        , buo     , n2      , s2     ,&
               Rt        , p_tke   , b_tke   , d_tke   )

! ---------------------------------- code history -------------------------------------
! Description:
!     Interface of Lake_CoLML. Used to initialize some of Lake_CoLML's variables 
!     and perform calculations of lake snow cover, lake flux, lake temperature, etc.
!
!
! Original author: 
!     Tiejun Ling
!
! Revisions:
!     Xin-Zhong Liang
!     Lei Sun, 2019: Modified for lake
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024: Reorganized the code to adapt to the new LakeDriver
!
! -------------------------------------------------------------------------------------
!================================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in)      :: &
      scwat                   ,&! surface category of water characteristics: 0: ocean, 1: shallow lake, 2: deep lake
      etaopt                  ,&! option for Extinction coefficient calculation
      slpeopt                 ,&! method of solar radiation penetration, penetration is defined in the levels
      mtpeopt                 ,&! method of momentum penetration, penetration is defined in the levels
      rhoopt                  ,&! option for water density calculation
      prtopt                  ,&! option for turbulent Prandtl number calculation
      elopt                   ,&! option for turbulence energy length scale calculation
      betaopt                 ,&! option for betaprime calculation, 1: constant, 2: equation
      nlayer                    ! Maximum number of layer (lake or ocean)

   logical, intent(in)      :: &
      tskBC                   ,&! option for lake surface temperature boundary condition, T=> USE SST boundary condition
      SSTNudg                   ! surface temperature nudging option, T=> USE SST Nudging

   real(r8), intent(in)     :: &
      dplak                   ,&! lake depth [m]
      tskin                   ,&! ground surface temperature [k]
      prtin                   ,&! User-defined turbulent Prandtl number
      rlat                    ,&! latitude in radians
      dtlak                   ,&! time step [s]
      rlaxtau                 ,&! relaxation time scale for lake temperature [s]
      etal                    ,&! extinction coefficient [1/m]
      btpri                   ,&! fraction of solar radiation in the NIR, only used when MODWMEJ
      obsdsst                 ,&! observed daily SST (previous day) [K]
      predsst                 ,&! predicted daily SST (previous day) [K]
      fgrnd1                  ,&! net heat flux from the atmosphere absorbed by top layer [W/m2]
      fgrnd                   ,&! total net heat flux from the atmosphere absorbed by lake [W/m2]
      sabg                    ,&! solar radiation absorbed by ground [W/m2]
      z0m                     ,&! roughness length over ground, momentum [m]
      taux                    ,&! wind stress: E-W [kg/m/s**2]
      tauy                      ! wind stress: N-S [kg/m/s**2]

   ! ----- input/output variables -----
   real(r8), intent(inout)  :: &
      dzlak(nlayer)           ,&! lake layer thickness (m)
      zlake(nlayer)           ,&! lake layer depth (m)
      zilak(nlayer+1)         ,&! lake layer interface level [m]
      lktmp(nlayer)           ,&! lake temperature (K)
      uwatv(nlayer)           ,&! Water velocity in x-direction [m/s]
      vwatv(nlayer)           ,&! Water velocity in y-direction [m/s]
      lksal(nlayer)           ,&! lksality [‰]
      lkrho(nlayer)           ,&! Water density [kg/m^3]
      tke(nlayer)             ,&! Turbulent kinetic energy (TKE) [J/kg]
      eps(nlayer)             ,&! TKE dissipation rate [W/kg]
      Km(nlayer)              ,&! eddy viscosity [m2s-1]
      Kh(nlayer)              ,&! diffusivity for scalar  [m2/s]
      Ke(nlayer)                ! diffusivity for tke  [m2/s]
            
   real(r8), intent(inout)  :: &
      hobl                      ! oceanic boundary layer depth [m]

!  ------------------------- output variables ---------------------------
   real(r8), intent(out)    :: &
      b0                      ,&! Surface Bouyancy flux
      buo(nlayer)             ,&! deltarho/lkrho mean buoyancy [ms-2]
      n2(nlayer)              ,&! brunt-vaisala frequency (/s**2)
      s2(nlayer)              ,&! shear number
      Rt(nlayer)              ,&! turbulent Richard number
      p_tke(nlayer)           ,&! TKE production by mean velocity shear [m2s-3]
      b_tke(nlayer)           ,&! TKE production by bouyancy flux [m2s-3]
      d_tke(nlayer)             ! TKE diffusion term [m2s-3]
        
!  ------------------------- local variables ---------------------------
   integer  :: iobl             ! mixed layer depth index
   real(r8) :: Khmin            ! minimum eddy viscosity [m2s-1]
   real(r8) :: gdbuo(nlayer)    ! gradient buoyancy [s-2] dB/dzlak
   real(r8) :: qcl(nlayer)      ! rms velocity of turbulence [=sqrt(2*tke) [ms-1]]
   real(r8) :: elnth(nlayer)    ! turbulence energy length scale [m]
   real(r8) :: prfmt(nlayer)    ! momentum profile
   real(r8) :: prfsl(nlayer)    ! solar penetration profile
!==============================================================================

      qcl = sqrt(2. * tke)

      !- radiation penetraion
      CALL radpene( &
            ! "in" arguments
            ! -------------------
            nlayer    , slpeopt , dplak   , zilak    ,&
            etaopt    , etal    , betaopt , btpri    ,&
            ! "out" arguments
            ! -------------------
            prfsl )

      ! momentum penetraion
      CALL mtmpene( &
            ! "in" arguments
            mtpeopt   , nlayer  , hobl     , zilak   ,&
            ! "out" arguments
            prfmt     )

      !- update lake temperature
      CALL solve_t( &
            ! "in" arguments
            nlayer    , SSTNudg , tskBC   , dtlak   ,&
            fgrnd1    , sabg    , btpri   , tskin   ,&
            obsdsst   , predsst , dzlak   , Kh      ,&
            prfsl     , lkrho  , rlaxtau ,&
            ! "inout" arguments
            lktmp   )

      !- update water density
      CALL update_rho ( &
            ! "in" arguments
            rhoopt    , nlayer   , lktmp    , lksal    ,&
            ! "out" arguments
            lkrho    )

      !- update buoyancy
      buo = buoyeval(nlayer, lktmp)

      !- update n2, s2
      CALL update_n2_s2( &
            ! "in" arguments
            nlayer    , lkrho   , uwatv   , vwatv   ,&
            dzlak     ,&
            ! "out" arguments
            n2        , s2       )

      ! update mean water velocity in x-direction [m/s]
      CALL solve_uv( &
            ! "in" arguments
            nlayer    , rlat     , taux    , dtlak   ,&
            dzlak     , prfmt    , Km      , lkrho  ,&
            ! "inout" arguments
            uwatv     )

      ! update mean water velocity in y-direction [m/s]
      CALL solve_uv( &
            ! "in" arguments
            nlayer    , rlat     , tauy    , dtlak   ,&
            dzlak     , prfmt    , Km      , lkrho  ,&
            ! "inout" arguments
            vwatv     )

      ! update turbulence kinetic energy (TKE) [J/kg]
      CALL solve_tke( &
            ! "in" arguments
            nlayer    , dtlak    , dzlak   , Kh      ,&
            Km        , Ke       , eps     , buo     ,&
            qcl       , uwatv    , vwatv   , n2      ,&
            ! "inout" arguments
            tke       )

      ! Update the diagnostic variables
      CALL update_diag( &
            ! "in" arguments
            nlayer    , dzlak    , zilak   , zlake   ,&
            uwatv     , vwatv    , lktmp   , lksal   ,&
            tke       , buo      , lkrho  , Km      ,&
            Kh        , Ke       , fgrnd   ,&
            ! "out" arguments
            b0        , n2       , s2      , p_tke   ,&
            b_tke     , d_tke    , qcl     )  

      ! Update the eddy viscosity and diffusivity 
      CALL update_eddy_pars( &
            ! "in" arguments
            ! -------------------
            prtopt   , elopt    , nlayer   ,&! riopt    ,&
            dplak    , dzlak    , zilak    , hobl     ,&
            z0m      , buo      , qcl      , lkrho   ,&
            n2       , s2       , sabg     , prtin    ,&
            ! "out" arguments
            gdbuo     , Km       , Kh      , Ke      ,&
            eps       , Rt       , elnth   , Khmin   )

      ! Update the mixed layer depth
      CALL update_hobl( &
            ! "in" arguments
            nlayer   , scwat    , dplak    , zilak   ,&
            tke      , lktmp    , gdbuo    , Khmin   ,&
            Kh       ,&
            ! "out" arguments
            iobl     , hobl     )

   END SUBROUTINE XOML_solve



   SUBROUTINE radpene( &
               ! "in" arguments
               ! -------------------
               nlayer    , slpeopt , dplak   , zilak    ,&
               etaopt    , etal    , betaopt , btpri    ,&
               ! "out" arguments
               ! -------------------
               prfsl )

! ---------------------------------- code history -------------------------------------
! Description:
!     Calculate the solar radiation penetration profile in the lake or ocean.
!
!
! Original author: 
!     Tiejun Ling, 2016
!
! Revisions:
!     Xin-Zhong Liang
!     Lei Sun, 2019: Modified for lake
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024: Reorganized the code to adapt to the new LakeDriver
!
! -------------------------------------------------------------------------------------
   USE module_sf_lake_table, only: SL_NOPE, SL_JRV, SL_BKS, SL_UNK, SL_SUB
!================================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in)      :: &
      slpeopt                 ,&! method of solar radiation penetration, penetration is defined in the levels
      etaopt                  ,&! option for Extinction coefficient calculation, 1: constant, 2: equation
      betaopt                 ,&! option for betaprime calculation, 1: constant, 2: equation
      nlayer                    ! Maximum number of layer (lake or ocean)

   real(r8), intent(in)     :: &
      dplak                   ,&! lake depth [m]
      btpri                   ,&! beta prime in Monin-Obukhov theory
      etal                    ,&! extinction coefficient [1/m]
      zilak(nlayer+1)           ! lake layer depth (m)

!  ------------------------- output variables ---------------------------
   real(r8), intent(out)    :: &
      prfsl(nlayer)             ! solar penetration profile

!  ------------------------- local variables ---------------------------
   integer,  parameter :: njrv=5+1
   real(r8), dimension(njrv) ::                               & !
                fjrvred=(/0.58, 0.62, 0.67, 0.77, 0.78, 0.78/),& !
                fjrvpar=(/0.42, 0.38, 0.33, 0.23, 0.22, 0.22/),& !
                kjrvred=(/0.35, 0.60, 1.00, 1.50, 1.40, 1.00/),& ! [m]
                kjrvpar=(/23.0, 20.0, 17.0, 14.0, 7.90, 7.90/)   ! [m]
   real(r8) :: fbkspar, fbksred, kbksred, kbkspar, eta, betaprime
   integer  :: k, jt
   real(r8) :: jrvtype= 4                                       ! Jerlov water type
!================================================================================
    
      !*******************************************************************************
      SELECT CASE(slpeopt)    ! method of solar radiation penetration, penetration is defined in the levels (interfaces)
      !*******************************************************************************
         !---------------
         CASE(SL_NOPE) !- no penetration
         !---------------
            prfsl = 0.0
            prfsl(1) = 1.0

         !---------------
         CASE(SL_JRV) !- classical Jerlov (jrv)
         !---------------
            jt = jrvtype
            DO k=1, nlayer
               prfsl(k) = fjrvred(jt) * exp(-zilak(k) / kjrvred(jt)) + &
                           fjrvpar(jt) * exp(-zilak(k) / kjrvpar(jt))
            ENDDO

         !---------------
         CASE(SL_BKS) !- Birol Kara scheme (bks)
         !---------------
            prfsl=0.
            kbksred = 0.5
            kbkspar = kbkspar
            fbkspar = max(0.27, 0.695-5.7*kbkspar)
            fbksred = 1.0 - fbkspar
            DO k=1, nlayer
               prfsl(k) = fbksred * exp(-zilak(k) / kbksred) + &
                           fbkspar * exp(-zilak(k) * kbkspar)
            ENDDO

         !---------------
         CASE(SL_UNK) !- UNKOWN
         !---------------
            DO k=1,nlayer
               prfsl(k)= 0.237 * exp(-zilak(k)*0.132) + 0.36 * exp(-zilak(k)*0.44)
            ENDDO
            prfsl(1)=prfsl(1)+0.403

         !---------------
         CASE(SL_SUB) !- Subin et al. (2012)
         !---------------
            DO k=1,nlayer
               IF (etaopt == 1) THEN
                  eta = etal
               else
                  eta = 1.1925 * max(dplak, 1.)**(-0.424)
               ENDIF
               prfsl(k)= exp( -eta * zilak(k) ) 
            ENDDO

         !---------------
         CASE DEFAULT
         !---------------
            write(*,*) 'Fatal: Wrong solar radiation penetration method selection'
            stop

      END SELECT

   END SUBROUTINE radpene



   SUBROUTINE mtmpene( &
               ! "in" arguments
               ! -------------------
               mtpeopt   , nlayer  , hobl    , zilak    ,&
               ! "out" arguments
               ! -------------------
               prfmt   )

! ---------------------------------- code history -------------------------------------
! Description:
!     Calculate the momentum penetration profile in the lake or ocean.
!
!
! Original author: 
!     Tiejun Ling, 2016
!
! Revisions:
!     Xin-Zhong Liang
!     Lei Sun, 2019: Modified for lake
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024: Reorganized the code to adapt to the new LakeDriver
!
! -------------------------------------------------------------------------------------
   USE module_sf_lake_table, only: MT_NOPE, MT_RHZ
!================================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in) :: mtpeopt                    ! method of momentum penetration, penetration is defined in the levels
   integer, intent(in) :: nlayer                     ! Maximum number of layer (lake or ocean)
   real(r8), intent(in) :: hobl                      ! height of ocean bounday layer
   real(r8), intent(in) :: zilak(nlayer+1)           ! lake layer interface depth (m)

!  ------------------------- output variables ---------------------------
   real(r8), intent(out) :: prfmt(nlayer)            ! momentum penetration profile

!  ------------------------- local variables ---------------------------
   integer :: k
!================================================================================

      ! *****************************************************************************
      SELECT CASE(mtpeopt)    ! sea surface momentum penetration option
      ! *****************************************************************************
         !---------------
         CASE(MT_NOPE) !- no momentum penetration
         !---------------
            prfmt(:) = 0
            prfmt(1) = 1.0  ! add 20130616

         !---------------
         CASE(MT_RHZ)  !- Rong-hua Zhang's scheme
         !---------------
            DO k=1, nlayer
#ifdef _XLAKE
               prfmt(k) = 1 - zilak(k) / hobl     !(hoblmax)  ! (hobl) importrant ?! ltj maybe get fix value.
#else
               prfmt(k) = 1 - zilak(k) / 50     !(hoblmax)  ! (hobl) importrant ?! ltj maybe get fix value.
#endif
               prfmt(k) = max(0.0, prfmt(k))
            ENDDO

         !---------------
         CASE DEFAULT
         !---------------
            write(*,*) 'Fatal: Wrong sea surface momentum penetration method selection'
            stop
      
      ! *****************************************************************************
      END SELECT
      ! *****************************************************************************

   END SUBROUTINE mtmpene



   SUBROUTINE solve_t( &
               ! "in" arguments
               ! -------------------
               nlayer  , SSTNudg , tskBC   , dtlak   ,&
               fgrnd1  , sabg    , btpri   , tskin   ,&
               obsdsst , predsst , dzlak   , Kh      ,&
               prfsl   , lkrho  , rlaxtau ,&
               ! "inout" arguments
               ! -------------------
               lktmp   )
                    
! ---------------------------------- code history -------------------------------------
! Description:
!     Solve the temperature profile in the lake or ocean.
!
! Original author: 
!     Tiejun Ling, 2016
!
! Revisions:
!     Xin-Zhong Liang
!     Lei Sun, 2019: Modified for lake
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024: Reorganized the code to adapt to the new LakeDriver
!
! -------------------------------------------------------------------------------------
      USE module_sf_lake_table, only: cpliq
!================================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in)      :: &
      nlayer                    ! Maximum number of layer (lake or ocean)

   logical, intent(in)      :: &
      SSTNudg                 ,&! SST Nudging
      tskBC                     ! boundary condition for lake temperature

   real(r8), intent(in)     :: &
      dtlak                   ,&! time step [s]
      fgrnd1                  ,&! net heat flux from the atmosphere absorbed by top layer [W/m2]
      sabg                    ,&! solar radiation absorbed by the ground [W/m2]
      btpri                   ,&! beta prime in Monin-Obukhov theory
      rlaxtau                 ,&! relaxation time scale for SST Nudging [s]
      obsdsst                 ,&! observed daily SST (previous day) [K]
      predsst                 ,&! predicted daily SST (previous day) [K]
      tskin                     ! ground temperature [K]
        
   real(r8), intent(in)     :: &
      dzlak(nlayer)           ,&! layer thickness (m)
      Kh(nlayer)              ,&! diffusivity for scalar [m2s-1]
      prfsl(nlayer)           ,&! solar radiation penetration profile
      lkrho(nlayer)             ! density [kg/m3]
      ! eta                      ! extinction coefficient [1/m]
        
!  ------------------------- inout variables ---------------------------
   real(r8), intent(inout)  :: &
      lktmp(nlayer)             ! temperature [K]

!  ------------------------- local variables ---------------------------
   real(r8), parameter :: beta = 0.4                       ! coefficient for the implicit scheme
   real(r8), dimension(nlayer) :: A                        ! Coefficients above the main diagonal
   real(r8), dimension(nlayer) :: B                        ! Coefficients below the main diagonal
   real(r8), dimension(nlayer) :: D                        ! Coefficients on the main diagonal
   real(r8), dimension(nlayer) :: X                        ! right-hand side 
   real(r8) :: dtz, pdz, ndz, rlx                          ! temporary variables
   integer  :: k                                           ! loop index
   real(r8) :: cp                                          ! volumetric heat capacity of water
   integer  :: l1                                          ! low boundary option (air-water interface)
   integer  :: lm                                          ! upper boundary option (water-sediment interface)
   real(r8) :: w1                                          ! low boundary values (air-water interface)
   real(r8) :: wm                                          ! upper boundary values (water-sediment interface)
   real(r8) :: eta                                         !+WMEJ unknown variable
!================================================================================

      ! Since l1=3 and lm=2 in Noh 2011, that
      ! means layer 1 will be used and nlayer not used in Trid.
      !- USE in the tridiagonal matrix solver 
      l1 = 3; w1 = tskin; lm = 2; wm = 0.        ! Assume that the bottom geothermal flux is 0
      eta = 1.0                                  !+WMEJ unknown variable, why is it used here?

      IF(SSTNudg) THEN !-WMEJ: SST Nudging
         rlx =  (obsdsst - predsst) / rlaxtau
      else
         rlx = 0.0
      ENDIF

      !-WMEJ: When correcting for SST Nudging, this cannot be 0, so this check is incorrect.
      ! IF (rlx .ne. 0) THEN
         ! print *, 'Fatal: rleax has been used'
         ! stop
      ! ENDIF

      !- for the top layer
      IF (tskBC) THEN  ! USE the air-water interface temperature as the boundary condition
         
         IF (l1.eq.3)  THEN
            k = 1
            cp  = lkrho(k) * cpliq     
            dtz = 2 * dtlak / dzlak(k)
            pdz = dzlak(k)
            ndz = dzlak(k) + dzlak(k+1)

            B(k) = 0.0
            D(k) = 1.0 + (1-beta) * eta * dtz *  &
                  ( nhf(nlayer, k, dzlak, Kh ) / ndz )
            A(k) = (1-beta) * eta * dtz * nhf(nlayer, k, dzlak, Kh ) / ndz
            X(k) = lktmp(k)                                                                        &!- current temperature
               + beta * dtz * nhf(nlayer, k, dzlak, Kh )/ndz * (lktmp(k+1)-lktmp(k))             &!- Heat flux from the layer below 
               + dtlak / dzlak(k) * Kh(1) * 2 * (tskin-lktmp(1)) / dzlak(1)                      &!- Heat flux from the air-water interface
               + dtlak / dzlak(k) * (1 - btpri) * sabg / cp * (prfsl(k) - prfsl(k+1))            &!- Solar radiation heating term
               + dtlak * rlx                                                                      !- Relaxation term (correction for SST Nudging)
         else
            print *,'Fatal: wrong boundary option for solve_t'
            stop
         ENDIF

      else  ! USE the net heat flux at the air-water interface as the boundary condition

         k = 1
         cp  = lkrho(k) * cpliq
         dtz = 2 * dtlak / dzlak(k)
         pdz = dzlak(k)
         ndz = dzlak(k) + dzlak(k+1)

         B(k) = 0.0
         D(k) = 1.0 + (1-beta) * eta * dtz *  &
                  ( nhf(nlayer, k, dzlak, Kh )/ndz )
         A(k) = (1-beta) * eta * dtz * nhf(nlayer, k, dzlak, Kh ) / ndz
         X(k) = lktmp(k)                                                                              &!- current temperature
                  + beta * dtz * nhf(nlayer, k, dzlak, Kh )/ndz * (lktmp(k+1) - lktmp(k))               &!- Heat flux from the layer below
                  + dtlak / dzlak(k) * fgrnd1 / cp                                                      &!- Heat flux from the air-water interface
                  + dtlak / dzlak(k) * (1 - btpri) * sabg / cp * (prfsl(k) - prfsl(k+1))                &!- Solar radiation heating term
                  + dtlak * rlx                                                                          !- Relaxation term (correction for SST Nudging)
      
      ENDIF

      !- for the other layers
      DO k = 2, nlayer-1

         cp  = lkrho(k) * cpliq
         dtz = 2 * dtlak / dzlak(k)
         pdz = dzlak(k) + dzlak(k-1)
         ndz = dzlak(k) + dzlak(k+1)

         B(k) = (1 - beta)  *  dtz  *  phf(nlayer, k, dzlak, Kh ) / pdz
         D(k) = 1.0 + (1 - beta)  *  dtz  *   &
                  ( nhf(nlayer, k, dzlak, Kh ) / ndz + phf(nlayer, k, dzlak, Kh ) / pdz )
         A(k) = (1 - beta)  *  dtz  *  nhf(nlayer, k, dzlak, Kh ) / ndz
         X(k) = lktmp(k)                                                                             &!- current temperature
                  + beta * dtz * ( nhf(nlayer, k, dzlak, Kh ) / ndz * (lktmp(k+1) - lktmp(k))          &!- Heat flux from the layer below
                  - phf(nlayer, k, dzlak, Kh ) / pdz * (lktmp(k) - lktmp(k-1)) )                       &!- Heat flux from the layer above
                  + dtlak / dzlak(k) * (1 - btpri) * sabg / cp * (prfsl(k) - prfsl(k+1))                !- Solar radiation heating term
      
      ENDDO

      CALL trid(nlayer, l1, w1, lm, wm, B, D, A, X, lktmp)

      return
   END SUBROUTINE solve_t



   SUBROUTINE update_rho( &
               ! "in" arguments
               ! -------------------
               rhoopt    , nlayer   , lktmp    , lksal    ,&
               ! "out" arguments
               ! -------------------
               lkrho    )

! ---------------------------------- code history -------------------------------------
! Description:
!     update the density of the lake water
!
!
! Original author: 
!     Tiejun Ling, 2016
!
! Revisions:
!     Xin-Zhong Liang
!     Lei Sun, 2019: Modified for lake
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024: Reorganized the code to adapt to the new LakeDriver
!
! -------------------------------------------------------------------------------------
   USE module_sf_lake_table, only: denh2o, tfrz, RHO_FIX, RHO_HOST, RHO_UNESCO
!================================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in)      :: &
      rhoopt                  ,&! method of density calculation
      nlayer                    ! Maximum number of layer (lake or ocean)

   real(r8), intent(in)     :: &
      lktmp(nlayer)           ,&! temperature [K]
      lksal(nlayer)             ! lksality [‰]

!  ------------------------- output variables ---------------------------
   real(r8), intent(out)    :: &
      lkrho(nlayer)             ! density [kg/m3]

!  ------------------------- local variables ---------------------------
   integer :: k
!================================================================================

      !- calculate the density of the lake water with the lksality effect
      DO k=1, nlayer
         ! *****************************************************************************
         SELECT CASE(rhoopt)  ! method of density calculation
         ! *****************************************************************************
            !---------------
            CASE(RHO_FIX) ! fixed value
            !---------------
               lkrho(k) = denh2o
            
            !---------------
            CASE(RHO_HOST) ! Hosteller
            !---------------
               lkrho(k) = denh2o* (1-1.9549e-05*(abs(lktmp(k)-277))**1.68)
            
            !---------------
            CASE(RHO_UNESCO) ! UNESCO formula
            !---------------
               lkrho(k) = 1000*(1+8.0e-5+5.88e-5*(lktmp(k)-tfrz)-8.11e-6*(lktmp(k)-tfrz)**2+4.77e-8*(lktmp(k)-tfrz)**3)

            !---------------
            CASE DEFAULT
            !---------------
               write(*,*) 'Fatal: wrong density calculation method selection in update_rho'
               stop
         END SELECT
      ENDDO

      return
   END SUBROUTINE update_rho



   SUBROUTINE update_n2_s2( &
               ! "in" arguments
               ! -------------------
               nlayer   , lkrho   , uwatv    , vwatv    ,&
               dzlak    ,&
               ! "out" arguments
               ! -------------------
               n2       , s2       )

! ---------------------------------- code history -------------------------------------
! Description:
!     Calculate the Brunt-Vaisala frequency and shear number
!
!
! Original author: 
!     Tiejun Ling, 2016
!
! Revisions:
!     Xin-Zhong Liang
!     Lei Sun, 2019: Modified for lake
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024: Reorganized the code to adapt to the new LakeDriver
!
! -------------------------------------------------------------------------------------
   USE module_sf_lake_table, only: grav, denh2o
!================================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in)      :: &
      nlayer                    ! Maximum number of layer (lake or ocean)

   real(r8), intent(in)     :: &
      lkrho(nlayer)           ,&! Density [kg/m3]
      uwatv(nlayer)           ,&! Water velocity in x-direction [m/s]
      vwatv(nlayer)           ,&! Water velocity in y-direction [m/s]
      dzlak(nlayer)             ! Layer thickness [m]

!  ------------------------- output variables ---------------------------
   real(r8), intent(out)    :: &
      n2(nlayer)              ,&! Brunt-Vaisala frequency [/s^2]
      s2(nlayer)                ! Shear number

!  ------------------------- local variables ---------------------------
   integer :: k
   real(r8) :: pdz, ndz
   real(r8), parameter :: beta = 0.4                       ! coefficient for the implicit scheme
!================================================================================

      DO k = 2, nlayer-1
         pdz = (dzlak(k) + dzlak(k-1)) / 2.
         ndz = (dzlak(k) + dzlak(k+1)) / 2.
         n2(k) = grav / lkrho(k) * (nhf(nlayer, k, dzlak, lkrho) - phf(nlayer, k, dzlak, lkrho)) / dzlak(k)
         s2(k) = ((nhf(nlayer, k, dzlak, uwatv) - phf(nlayer, k, dzlak, uwatv)) / dzlak(k))**2 + &
                  ((nhf(nlayer, k, dzlak, vwatv) - phf(nlayer, k, dzlak, vwatv)) / dzlak(k))**2
      ENDDO

      n2(1) = n2(2)
      n2(nlayer) = n2(nlayer-1)
      s2(1) = s2(2)
      s2(nlayer) = s2(nlayer-1)

      return
   END SUBROUTINE update_n2_s2




   SUBROUTINE solve_uv( &
               ! "in" arguments
               ! -------------------
               nlayer   , rlat    , tau     , dtlak   ,&
               dzlak    , prfmt   , Km      , lkrho  ,&
               ! "inout" arguments
               ! -------------------
               watv     )
                  
! ---------------------------------- code history -------------------------------------
! Description:
!     Solve the water velocity in x- and y-directions
!     To simplify the code, the original solve_u and solve_v routines have been merged.
!     Now, the solution for U and V is handled by a single routine.
!
! Original author: 
!     Tiejun Ling, 2016
!
! Revisions:
!     Xin-Zhong Liang
!     Lei Sun, 2019: Modified for lake
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024: Reorganized the code to adapt to the new LakeDriver
!
! -------------------------------------------------------------------------------------
   USE module_sf_lake_table, only: omega
!================================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in)      :: &
      nlayer                    ! Maximum number of layer (lake or ocean)

   real(r8), intent(in)     :: &
      rlat                    ,&! Latitude [rad]
      tau                     ,&! Wind stress [N/m2]
      dtlak                   ,&! time step [s]
      dzlak(nlayer)           ,&! layer thickness [m]
      prfmt(nlayer)           ,&! Momentum penetration profile
      Km(nlayer)              ,&! eddy viscosity [m2s-1]
      lkrho(nlayer)             ! Density [kg/m3]

!  ------------------------- inout variables ---------------------------
   real(r8), intent(inout)  :: &
      watv(nlayer)              ! Water velocity [m/s]

!  ------------------------- local variables ---------------------------
   integer  :: k
   real(r8), dimension(nlayer) :: A       ! Coefficients above the main diagonal
   real(r8), dimension(nlayer) :: B       ! Coefficients below the main diagonal
   real(r8), dimension(nlayer) :: D       ! Coefficients on the main diagonal
   real(r8), dimension(nlayer) :: X       ! right-hand side 
   real(r8) :: dtz, pdz, ndz              ! temporary variables
   integer  :: l1                         ! low boundary option (air-water interface)
   integer  :: lm                         ! upper boundary option (water-sediment interface)
   real(r8) :: w1                         ! low boundary values (air-water interface)
   real(r8) :: wm                         ! upper boundary values (water-sediment interface)
   real(r8) :: fcor                       ! Coriolis parameter [s-1]
   real(r8), parameter :: beta = 0.4      ! coefficient for the implicit scheme
!================================================================================

      !Since l1=2 and lm=2 in Noh 2011, that
      !means both layer 1 and nlayer not were used in Trid.
      !- USE in the tridiagonal matrix solver 
      ! l1 = 2; lm = 2; w1 = - tau / lkrho(1) * (prfmt(1) - prfmt(2)) * dzlak(1) / Km (1); wm = 0_r8     ! Noh 2011
      l1 = 2; lm = 1; w1 = - tau / lkrho(1) * (prfmt(1) - prfmt(2)) * dzlak(1) / Km (1); wm = 0_r8       ! Lei Sun or Tiejun Ling modified

      fcor = 2 * omega * sin(rlat) 
      
      A=0.0
      B=0.0
      D=0.0
      X=0.0
      
      DO k = 2, nlayer-1

         dtz = 2_r8 * dtlak / dzlak(k)
         pdz = dzlak(k) + dzlak(k-1)
         ndz = dzlak(k) + dzlak(k+1)

         B(k) = (1_r8 - beta) * dtz * phf(nlayer, k, dzlak, Km ) / pdz
         D(k) = 1_r8 + (1_r8 - beta) * dtz *  &
                  ( nhf(nlayer, k, dzlak, Km ) / ndz + phf(nlayer, k, dzlak, Km ) / pdz )
         A(k) = (1_r8 - beta) * dtz * nhf(nlayer, k, dzlak, Km ) / ndz
         X(k) = watv(k)                                                                         &!- current velocity
                  + beta * dtz * ( nhf(nlayer, k, dzlak, Km ) / ndz * (watv(k+1) - watv(k))       &!- production by mean velocity shear
                                 - phf(nlayer, k, dzlak, Km ) / pdz * (watv(k) - watv(k-1)))      &
                  + dtlak / dzlak(k) * tau/lkrho(k) * (prfmt(k) - prfmt(k+1))                   !&!- Wind stress term
                  ! + dtlak * fcor * watv(k)
      ENDDO
                
      CALL trid(nlayer, l1, w1, lm, wm, B, D, A, X, watv)

      return
   END SUBROUTINE solve_uv



   SUBROUTINE solve_tke( &
               ! "in" arguments
               ! -------------------
               nlayer   , dtlak   , dzlak   , Kh      ,&
               Km       , Ke      , eps     , buo     ,&
               qcl      , uwatv   , vwatv   , n2      ,&
               ! "inout" arguments
               ! -------------------
               tke      )
                  
! ---------------------------------- code history -------------------------------------
! Description:
!     Solve the turbulence kinetic energy (TKE) equation
!
! Original author: 
!     Tiejun Ling, 2016
!
! Revisions:
!     Xin-Zhong Liang
!     Lei Sun, 2019: Modified for lake
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024: Reorganized the code to adapt to the new LakeDriver
!
! -------------------------------------------------------------------------------------
!================================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in)      :: &
      nlayer                    ! Maximum number of layer (lake or ocean)

   real(r8), intent(in)     :: &
      dtlak                   ,&! time step [s]
      dzlak(nlayer)           ,&! layer thickness [m]
      Km(nlayer)              ,&! eddy viscosity [m2s-1]
      Kh(nlayer)              ,&! diffusivity for scalar [m2s-1]
      Ke(nlayer)              ,&! eddy viscosity for TKE [m2s-1]
      eps(nlayer)             ,&! turbulent dissipation rate [m2s-3]
      buo(nlayer)             ,&! deltarho/lkrho mean buoyancy [ms-2]
      qcl(nlayer)             ,&! rms velocity of turbulence [=sqrt(2*tke) [ms-1]]
      uwatv(nlayer)           ,&! Water velocity in x-direction [m/s]
      vwatv(nlayer)           ,&! Water velocity in y-direction [m/s]
      n2(nlayer)                ! Brunt-Vaisala frequency [/s^2]
            
!  ------------------------- inout variables ---------------------------
   real(r8), intent(inout)  :: &
      tke(nlayer)               ! Turbulence kinetic energy [J/kg]

!  ------------------------- local variables ---------------------------
   integer  :: k                       
   real(r8), dimension(nlayer) :: A          ! Coefficients above the main diagonal
   real(r8), dimension(nlayer) :: B          ! Coefficients below the main diagonal
   real(r8), dimension(nlayer) :: D          ! Coefficients on the main diagonal
   real(r8), dimension(nlayer) :: X          ! right-hand side 
   real(r8) :: dtz, pdz, ndz                 ! temporary variables
   integer  :: l1                            ! low boundary option (air-water interface)
   integer  :: lm                            ! upper boundary option (water-sediment interface)
   real(r8) :: w1                            ! low boundary values (air-water interface)
   real(r8) :: wm                            ! upper boundary values (water-sediment interface)
   real(r8), parameter :: beta = 0.4         ! coefficient for the implicit scheme
!================================================================================

      !Since l1=2 and lm=2 in Noh 2011, that
      !means both layer 1 and nlayer not were used in Trid.
      l1 = 3; lm = 2; w1 = 0_r8; wm = 0_r8       ! Lei Sun or Tiejun Ling modified
      B = 0_r8; D = 0_r8; A = 0_r8; X = 0_r8

      IF (l1 .eq. 3) THEN

         k = 1
         dtz = 2 * dtlak / dzlak(k)
         ndz = dzlak(k) + dzlak(k+1)

         B(k) = 0.0
         D(k) = 1.0 + (1 - beta) *  dtz *  &
                  ( nhf(nlayer, k, dzlak, Ke) / ndz )
         A(k) = (1 - beta) *  dtz * nhf(nlayer, k, dzlak, Ke) / ndz
         X(k) = tke(k)                                                                         &!- current TKE
                  + beta * dtz* nhf(nlayer, k, dzlak, Ke) / ndz * (tke(k+1) - tke(k))            &!- TKE diffusion term
                  + dtlak / dzlak(k) * w1                                                         !- Top boundary condition

      ENDIF

      DO k = 2, nlayer-1

         dtz = 2 * dtlak / dzlak(k)
         pdz = dzlak(k) + dzlak(k-1)
         ndz = dzlak(k) + dzlak(k+1)

         B(k) = (1 - beta) * dtz * phf(nlayer, k, dzlak, Ke) / pdz
         D(k) = 1.0 + (1 - beta) * dtz *  &
                  ( nhf(nlayer, k, dzlak, Ke) / ndz + phf(nlayer, k, dzlak, Ke) / pdz )   &
                  + eps(k) * dtlak / (0.5 * qcl(k) * qcl(k))
         A(k) = (1 - beta) * dtz * nhf(nlayer, k, dzlak, Ke) / ndz
         X(k) = tke(k)                                                                                              &!- current TKE
                  ! + dtlak*Kh (k)*(nhf(nlayer,k,dzlak,buo)-phf(nlayer,k,dzlak,buo))/dzlak(k)  &
                  - dtlak * Kh (k) * n2(k)                                                                            &!- TKE production by bouyancy flux
                  + beta * dtz * ( nhf(nlayer, k, dzlak, Ke) / ndz * (tke(k+1) - tke(k))                              &!- TKE diffusion term
                                 - phf(nlayer, k, dzlak, Ke) / pdz * (tke(k) - tke(k-1)) )                            &
                  + dtlak * Km (k) * ((nhf(nlayer, k, dzlak, uwatv) - phf(nlayer, k, dzlak, uwatv)) / dzlak(k))**2    &!- TKE production by mean velocity shear
                  + dtlak * Km (k) * ((nhf(nlayer, k, dzlak, vwatv) - phf(nlayer, k, dzlak, vwatv)) / dzlak(k))**2

      ENDDO

      CALL trid(nlayer, l1, w1, lm, wm, B, D, A, X, tke)

      DO k=1, nlayer
         IF(tke(k) <= 1.0e-8) THEN
               tke(k) = 1.0e-8
         ENDIF
      ENDDO

      return
   END SUBROUTINE solve_tke



   SUBROUTINE update_diag( &
               ! "in" arguments
               ! -------------------
               nlayer   , dzlak   , zilak   , zlake   ,&
               uwatv    , vwatv   , lktmp   , lksal   ,&
               tke      , buo     , lkrho  , Km      ,&
               Kh       , Ke      , fgrnd   ,&
               ! "out" arguments
               ! -------------------
               b0       , n2      , s2      , p_tke   ,&
               b_tke    , d_tke   , qcl     )                                       

! ---------------------------------- code history -------------------------------------
! Description:
!     Update the diagnostic variables: surface buoyancy flux, Brunt-Vaisala frequency, shear number,
!     TKE production by mean velocity shear, TKE production by buoyancy flux, and TKE diffusion term
!
!
! Original author: 
!     Tiejun Ling, 2016
!
! Revisions:
!     Xin-Zhong Liang
!     Lei Sun, 2019: Modified for lake
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024: Reorganized the code to adapt to the new LakeDriver
!
! -------------------------------------------------------------------------------------
   USE module_sf_lake_table, only: grav, cpliq
!================================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in)      :: &
      nlayer                    ! Maximum number of layer (lake or ocean)

   real(r8), intent(in)     :: &
      fgrnd                   ,&! total net heat flux from the atmosphere absorbed by lake [W/m2]
      dzlak(nlayer)           ,&! layer thickness [m]
      zilak(nlayer)           ,&! layer interface depth [m]
      zlake(nlayer)           ,&! layer node depth [m]
      uwatv(nlayer)           ,&! Water velocity in x-direction [m/s]
      vwatv(nlayer)           ,&! Water velocity in y-direction [m/s]
      lktmp(nlayer)           ,&! temperature [K]
      lksal(nlayer)           ,&! lksality [‰]
      lkrho(nlayer)           ,&! Density [kg/m3]
      tke(nlayer)             ,&! Turbulence kinetic energy [J/kg]
      buo(nlayer)             ,&! deltarho/lkrho mean buoyancy [ms-2]
      Km(nlayer)              ,&! eddy viscosity [m2s-1]
      Kh(nlayer)              ,&! diffusivity for scalar [m2s-1]
      Ke(nlayer)                ! diffusivity for tke [m2s-1]
        
!  ------------------------- output variables ---------------------------
   real(r8), intent(out)    :: &
      b0                      ,&! Surface Bouyancy flux
      n2(nlayer)              ,&! Brunt-Vaisala frequency [/s^2]
      s2(nlayer)              ,&! Shear number
      p_tke(nlayer)           ,&! TKE production by mean velocity shear [m2s-3]
      b_tke(nlayer)           ,&! TKE production by bouyancy flux [m2s-3]
      d_tke(nlayer)           ,&! TKE diffusion term [m2s-3]
      qcl(nlayer)               ! rms velocity of turbulence [=sqrt(2*tke) [ms-1]]

!  ------------------------- local variables ---------------------------
   integer :: k
   real(r8) :: pdz, ndz
   real(r8) :: tcoef = -2.7e-4                             ! thermal expansion(1/rho0*drho/dT) [K-1] in s0=34.2 and t0=29.
!================================================================================

      ! tcoef=(-68.0+18.2091*(t-tdk)-0.30866*(t-tdk)**2+5.3445e-3*(t-tdk)**3 &                             !- Expression in buoyeval_old function
      !         -6.0721e-5*(t-tdk)**4+3.1441e-7*(t-tdk)**5)*1.e-6*-1

      ! b0 = grav * tcoef * (sabg+atmfrcs%nlw+atmfrcs%lhf+atmfrcs%shf)/(lkrho(1)*cpliq)                      !- check the sign of lhf and shf
      b0 = grav * tcoef * (fgrnd) / (lkrho(1) * cpliq)                                                        

      DO k = 2, nlayer-1
         pdz = (dzlak(k) + dzlak(k-1)) / 2.
         ndz = (dzlak(k) + dzlak(k+1)) / 2.

         p_tke(k) = Km (k) * ((nhf(nlayer, k, dzlak, uwatv) - phf(nlayer, k, dzlak, uwatv)) / dzlak(k))**2 + &
                  Km (k) * ((nhf(nlayer, k, dzlak, vwatv) - phf(nlayer, k, dzlak, vwatv)) / dzlak(k))**2
         d_tke(k) = (nhf(nlayer, k, dzlak, Ke) / ndz * (tke(k+1) - tke(k))        &
                  - phf(nlayer, k, dzlak, Ke) / pdz * (tke(k) - tke(k-1))) / dzlak(k)
         n2(k) = grav / lkrho(k) * (nhf(nlayer, k, dzlak, lkrho) - phf(nlayer, k, dzlak, lkrho)) / dzlak(k)
         b_tke(k) = -Kh (k) * n2(k)
         ! b_tke(k) = Kh (k) * (nhf(nlayer, k, dzlak, buo) - phf(nlayer, k, dzlak, buo)) / dzlak(k)
         s2(k) = ((nhf(nlayer, k, dzlak, uwatv) - phf(nlayer, k, dzlak, uwatv)) / dzlak(k))**2    &
               + ((nhf(nlayer, k, dzlak, vwatv) - phf(nlayer, k, dzlak, vwatv)) / dzlak(k))**2

      ENDDO

      ! top layer
      p_tke(1) = p_tke(2)
      b_tke(1) = b_tke(2)
      d_tke(1) = d_tke(2)
      n2(1) = n2(2)
      s2(1) = s2(2)

      ! bottom layer
      p_tke(nlayer) = p_tke(nlayer-1)
      b_tke(nlayer) = p_tke(nlayer-1)
      d_tke(nlayer) = d_tke(nlayer-1)
      n2(nlayer) = n2(nlayer-1)
      s2(nlayer) = s2(nlayer-1)

      ! update qcl
      DO k=1, nlayer
         qcl(k) = sqrt(2 * tke(k))
      ENDDO

   END SUBROUTINE update_diag



   SUBROUTINE update_eddy_pars( &
               ! "in" arguments
               ! -------------------
               prtopt   , elopt    , nlayer   ,&! riopt    ,&
               dplak    , dzlak    , zilak    , hobl     ,&
               z0m      , buo      , qcl      , rho      ,&
               n2       , s2       , sabg     , prtin    ,&
               ! "out" arguments
               ! -------------------
               gdbuo    , Km       , Kh       , Ke       ,&
               eps      , Rt       , elnth    , Khmin    )

! ---------------------------------- code history -------------------------------------
! Description:
!     Update the the vertical eddy viscosity Km, eddy diffusivity of temperature and lksality Kh
!     and the eddy diffusivity of turbulent kinetic energy Ke (Kk)
!
!     In this SUBROUTINE we follow the naming convention of Noh et al. (2011).
!     
!     the orginal Noh et al. (2011) parameterization:
!         Km = Sm * q * l       [Sm is equal to Ck in Sun et al. 2019]
!         Kh = Sh * q * l       [Sh is equal to Ck/PrtMY in Sun et al. 2019, PrtMY=0.8 ref to Noh and Kim 1999 and Mellor and Yamada 1982]
!         Ke = Se * q * l       [Se is equal to Cq in Sun et al. 2019]
!         q = sqrt(2 * tke) is the root mean square velocity scale
!         l = 1/(1/lg+1/h)  is the turbulence energy length scale
!         Sm = 0.39 (=Sm0), Sm/Sh = Prt, Sm/Se = Sigma [=1.95] are empirical constants [no stratification], Prt equals to PrtMY
!         PrtMY = 0.8 following Mellor and Yamada 1982. 
!        Sugested by Noh (1996), PrtMY and Sigma are independent of Rt.
!        It is suggested that the proportions of Sh and Se to Sm are constant.
!         
!         lg = k(z+z0), is geometric scale
!         k=0.4 is the von Karman constant, z0 is the roughness length, h is the mixing layer height
!
!         Note: The stratification (stability) of water bodies has an important inhibitory effect on the scale of turbulence.
!             Therefore, in this case we can parameterize the effects of stratification as:
!                 Sm / Sm0 = (1 + alpha*Rt)^(-1/2), alpha is a tunable parameter (=120 and 1000 in Noh et al. 2001, =50 in Noh and Kim 2011)
!                 same for Sh.
!         The parmesantrization described by Sun et al., 2019, is different here.
!         we have retained both schemes.       
!
! Original author: 
!     Tiejun Ling, 2016
!
! Revisions:
!     Xin-Zhong Liang
!     Lei Sun, 2019: Modified for lake
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024: Reorganized the code to adapt to the new LakeDriver
!
! -------------------------------------------------------------------------------------
   USE module_sf_lake_table, only: grav, cpliq, vonkar
!================================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in)      :: &
      prtopt                  ,&! Turbulent Prandtl number option
      ! riopt                   ,&! Richardson number option
      elopt                   ,&! turbulence length scale option
      nlayer                    ! Maximum number of layer (lake or ocean)

   real(r8), intent(in)     :: &
      dplak                   ,&! depth of the lake [m]
      z0m                     ,&! roughness length for momentum [m]
      sabg                    ,&! solar absorption [W/m2]
      prtin                   ,&! User-defined turbulent Prandtl number
      hobl                      ! ocean mixing layer height [m]

   real(r8), intent(in)     :: &
      dzlak(nlayer)           ,&! layer thickness [m]
      zilak(nlayer+1)         ,&! layer interface depth [m]
      buo(nlayer)             ,&! deltarho/lkrho mean buoyancy [ms-2]
      qcl(nlayer)             ,&! rms velocity of turbulence [=sqrt(2*tke) [ms-1]]
      rho(nlayer)             ,&! density [kg/m3]
      n2(nlayer)              ,&! Brunt-Vaisala frequency [/s^2]
      s2(nlayer)                ! Shear number

!  ------------------------- output variables ---------------------------
   real(r8), intent(out)    :: &
      Khmin                   ,&! minimum value for Kh
      gdbuo(nlayer)           ,&! gradient buoyancy [s-2] dB/dz
      Km(nlayer)              ,&! eddy viscosity [m2s-1]
      Kh(nlayer)              ,&! diffusivity for scalar [m2s-1]
      Ke(nlayer)              ,&! diffusivity for tke [m2s-1]
      eps(nlayer)             ,&! turbulent dissipation rate [m2s-3]
      Rt(nlayer)              ,&! turbulent Richard number
      elnth(nlayer)             ! turbulence length scale [m]
        
!  ------------------------- local variables ---------------------------
   integer :: k
   real(r8), parameter :: Sm0 = 0.39       ! default value for Sm (Ck)
   real(r8), parameter :: PrtMY = 0.8      ! Prandtl number for heat, following Mellor and Yamada 1982
   real(r8), parameter :: Ceps0 = 0.06     ! default value for Ceps
   real(r8), parameter :: Kmmax = 0.5      ! maximum value for Km
   real(r8), parameter :: Kmmin = 1.0e-4   ! minimum value for Km
   real(r8) :: alpha = 50.0    !120        ! this is a tunable parameter. !-be careful with this value
   real(r8) :: tcoef = -2.7e-4             ! thermal expansion(1/rho0*drho/dT) [K-1] in s0=34.2 and t0=29.
   real(r8) :: Sm, Sh, Se                  ! empirical constants, see subroutines description
   real(r8) :: cp                          ! specific heat capacity [J/kg/K]
   real(r8) :: Ri(nlayer)                  ! Richardson number
   real(r8) :: prandt                      ! turbulent Prandtl number
   real(r8) :: Ceps                        ! turbulent dissipation rate coefficient
!================================================================================
      Khmin = Kmmin   !initialize Khmin
      ! calculaiton of diffusivity, dissiation & length scale
      DO k = 2, nlayer-1
         gdbuo(k) = (buo(k+1) - buo(k-1)) / (dzlak(k+1) * 0.5 + dzlak(k) + 0.5 * dzlak(k-1))
      ENDDO

      ! calculate specific heat capacity
      cp=rho(1)*cpliq
      !==============================================================================
      ! tcoef=(-68.0+18.2091*(t-tdk)-0.30866*(t-tdk)**2+5.3445e-3*(t-tdk)**3 &                             !- Expression in buoyeval_old function
      !         -6.0721e-5*(t-tdk)**4+3.1441e-7*(t-tdk)**5)*1.e-6*-1
      gdbuo(1) = grav * tcoef * sabg / cp / Kh(1)
      gdbuo(nlayer) = 0.0

      ! calculate the eddy viscosity and diffusivity
      DO k = 1, nlayer
         ! calculate the turbulence length scale [m]
         elnth(k) = cal_elnth(elopt, hobl, zilak(k), z0m, dplak)

         ! calculate the Richardson number
         Ri(k) = cal_ri(1, n2(k), s2(k), elnth(k), qcl(k))       !- Richardson scheme
         Rt(k) = cal_ri(2, n2(k), s2(k), elnth(k), qcl(k))       !- Noh et al. 2011
         Rt(k) = max(Rt(k), 0.0)                                 !- Noh et al. 2011

         ! calculate the Prandtl number
         prandt = cal_prandt(prtopt, prtin, Ri(k))
         Khmin = min(Kmmin / prandt, Khmin)
         
         !-WMEJ: I don't know why the following code is commented out!
         !==============================================================================================================
         ! IF (gdbuo(k).gt. 1.0e-6) THEN
         !     ! --- unstable
         !     elnth(k) = elnth(k)                            
         !     Sm = Sm0                                                !- empirical constants for vertical eddy viscosity
         !     Sh = Sm / prandt                                        !- empirical constants for vertical eddy diffusivity of temperature and lksality
         !     Se = Sm / 1.95                                          !- empirical constants for vertical eddy diffusivity of turbulent kinetic energy
         !     Ceps = Ceps0                                            !- empirical constants for turbulent dissipation rate
               
         !     Km(k) = Sm * qcl(k) * elnth(k)                          !- calculate the vertival eddy viscosity
         !     Kh(k) = Sh * qcl(k) * elnth(k)                          !- calculate the eddy diffusivity of temperature and lksality
         !     Ke(k) = Se * qcl(k) * elnth(k)                          !- calculate the eddy diffusivity of turbulent kinetic energy
         !     eps(k) = Ceps * (qcl(k)**3) / elnth(k)                  !- calculate the turbulent dissipation rate
         ! else
         !     ! --- stable
         !     elnth(k) = elnth(k) / sqrt(1.0 + alpha * Rt(k))   
         !     Sm = Sm0 / sqrt(1.0 + alpha * Rt(k))                    !- empirical constants for vertical eddy viscosity
         !     Sh = Sm / prandt                                        !- empirical constants for vertical eddy diffusivity of temperature and lksality
         !     Se = Sm / 1.95                                          !- empirical constants for vertical eddy diffusivity of turbulent kinetic energy
         !     Ceps = Ceps0 / sqrt(1.0 + alpha * Rt(k))                !- empirical constants for turbulent dissipation rate
               
         !     Km(k) = Sm * qcl(k) * elnth(k)                          !- calculate the vertival eddy viscosity
         !     Km(k) = max(Kmmin, Km(k))
         !     Kh(k) = Sh * qcl(k) * elnth(k) !/ sqrt(1.0+bl*Rt(k))    !- calculate the eddy diffusivity of temperature and lksality
         !     Ke(k) = Se * qcl(k) * elnth(k)                          !- calculate the eddy diffusivity of turbulent kinetic energy
         !     eps(k) = Ceps * (qcl(k)**3) / elnth(k)                  !- calculate the turbulent dissipation rate
         ! ENDIF
         !==============================================================================================================

         ! calculate empirical constants 
         elnth(k) = elnth(k) / sqrt(1.0 + alpha * Rt(k))   
         Sm = Sm0 / sqrt(1.0 + alpha * Rt(k))                          !- empirical constants for vertical eddy viscosity
         Sh = Sm / prandt                                              !- empirical constants for vertical eddy diffusivity of temperature and lksality
         Se = Sm / 1.95                                                !- empirical constants for vertical eddy diffusivity of turbulent kinetic energy 
         Ceps = Ceps0 / sqrt(1.0 + alpha * Rt(k))                      !- empirical constants for turbulent dissipation rate
         print *, 'update_eddy_pars: prandt = ', prandt
         !- calculate the eddy viscosity
         Km(k) = Sm * qcl(k) * elnth(k)                                !- calculate the vertival eddy viscosity
         Km(k) = max(Kmmin, Km(k))
         Kh(k) = Sh * qcl(k) * elnth(k) !/ sqrt(1.0+bl*Rt(k))          !- calculate the eddy diffusivity of temperature and lksality
         Ke(k) = Se * qcl(k) * elnth(k)                                !- calculate the eddy diffusivity of turbulent kinetic energy
         eps(k) = Ceps * (qcl(k)**3) / elnth(k)                        !- calculate the turbulent dissipation rate

      ENDDO

      Km=Km+1.5e-6
      Kh=Kh+1.5e-7

      return
   END SUBROUTINE update_eddy_pars



   SUBROUTINE update_hobl( &
               ! "in" arguments
               ! -------------------
               nlayer   , scwat    , dplak    , zilak    ,&
               tke      , lktmp    , gdbuo    , Khmin    ,&
               Kh       ,&
               ! "out" arguments
               ! -------------------
               iobl     , hobl     )

! ---------------------------------- code history -------------------------------------
! Description:
!     Update the ocean boundary height or mixing layer height
!
! TODO: NEED CHECK THIS SUBROUTINE
!
!
! Original author: 
!     Tiejun Ling, 2016
!
! Revisions:
!     Xin-Zhong Liang
!     Lei Sun, 2019: Modified for lake
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024: Reorganized the code to adapt to the new LakeDriver
!
!--------------------------------------------------------------------------------------
   USE module_sf_lake_table, only: OCEAN
!==============================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in) :: nlayer                            ! Maximum number of layer (lake or ocean)
   integer, intent(in) :: scwat                             ! water body type
   real(r8), intent(in) :: dplak                            ! maximum ocean mixing layer height [m] [= dplak in lake and = 100.0 in ocean]
   real(r8), intent(in) :: tke(nlayer)                      ! Turbulence kinetic energy [J/kg]
   real(r8), intent(in) :: lktmp(nlayer)                    ! lake temperature [K]
   real(r8), intent(in) :: gdbuo(nlayer)                    ! gradient buoyancy [s-2] dB/dz
   real(r8), intent(in) :: zilak(nlayer+1)                  ! layer interface depth [m]
   real(r8), intent(in) :: Kh(nlayer)                       ! diffusivity for scalar [m2s-1]
   real(r8), intent(in) :: Khmin                            ! minimum value for Kh

!  ------------------------- output variables ---------------------------
   integer, intent(out) :: iobl                             ! mixing layer height index
   real(r8), intent(out) :: hobl                            ! ocean mixing layer height [m]        
   
!  ------------------------- local variables ---------------------------
   integer :: mh1, mh2, mh, k
   real(r8) :: tkemin, gdbuomax, hoblmax
   character(len=12) :: fstr
   real(r8), dimension(nlayer) :: ltke
!================================================================================

      iobl = nlayer - 1
      hobl = zilak(nlayer-1)
      tkemin = tke(1) * 1.0e-4
      mh1 = 1

      IF (scwat == OCEAN) THEN
         hoblmax = 100.0
      else
         hoblmax = dplak
      ENDIF

      DO k = 1, nlayer-1
         IF (tke(k) >= tkemin) THEN
               mh1 = mh1 +1
         else
               exit
         ENDIF
      ENDDO

      mh2 = 1
      gdbuomax = gdbuo(1)
      DO k = 2, nlayer-1
#ifdef _NDBC
         IF (gdbuo(k) >= gdbuomax*1.0e-2) THEN
#else
         IF(gdbuo(k) >= gdbuomax) THEN
#endif
            !gdbuomax = abs(gdbuo(k))
            mh2 = k
            exit
         ENDIF
      ENDDO

      IF (mh2 .ge. 2 .and. mh1 .ge. nlayer-2) THEN
         mh = mh2
      else
         mh = max(mh1, mh2)
      ENDIF

#ifdef _IMET
      DO k = 2, nlayer-1
         IF (abs(lktmp(k) - lktmp(1)) .gt. 0.4) THEN
            mh = k
            exit
         ENDIF
      ENDDO
#endif

#ifdef _NOH
      DO k = 2, nlayer-1
         IF (Kh(k) .lt. Khmin) THEN
            mh = k
            exit
         ENDIF
      ENDDO
#endif

      iobl = max(1, min(mh, nlayer-1))
      hobl = min(hoblmax, zilak(iobl+1))

      return 

   END SUBROUTINE update_hobl



   SUBROUTINE roughness_lake (snl, t_grnd, forc_psrf, ur, ustar, fetch, depth, &
                              z0mg,z0hg,z0qg)

! ---------------------------------- code history -------------------------------------
! Description:
!     Update the roughness length over ground for momentum, sensible heat, and latent heat
!
! -------------------------------------------------------------------------------------
!================================================================================
!  ------------------------- input variables ---------------------------
   integer,  intent(in) :: snl            ! number of snow layers
   real(r8), intent(in) :: t_grnd         ! ground temperature
   real(r8), intent(in) :: forc_psrf      ! atmosphere pressure at the surface [pa]
   real(r8), intent(in) :: ur             ! wind effect
   real(r8), intent(in) :: fetch          ! lake fetch [m]
   real(r8), intent(in) :: depth          ! lake depth [m]
   ! real(r8), intent(in) :: cur          ! Charnock parameter (-)
   real(r8), intent(in) :: ustar          ! u* in similarity theory [m/s]

!  ------------------------- output variables ---------------------------
   real(r8), intent(out) :: z0mg          ! roughness length over ground, momentum [m]
   real(r8), intent(out) :: z0hg          ! roughness length over ground, sensible heat [m]
   real(r8), intent(out) :: z0qg          ! roughness length over ground, latent heat [m]
   real(r8), parameter :: cus = 0.1       ! empirical constant for roughness under smooth flow
   real(r8), parameter :: kva0 = 1.51e-5  ! kinematic viscosity of air (m^2/s) at 20C and 1.013e5 Pa
   real(r8), parameter :: prn = 0.713     ! Prandtl # for air at neutral stability
   real(r8), parameter :: sch = 0.66      ! Schmidt # for water in air at neutral stability
   real(r8), parameter :: vonkar = 0.4    ! von Karman constant [-]
   real(r8), parameter :: tfrz   = 273.16 ! freezing temperature [K]
   real(r8), parameter :: grav   = 9.80616! gravity constant [m/s2]
   real(r8), parameter :: cur0 = 0.01     ! min. Charnock parameter
   real(r8), parameter :: curm = 0.1      ! maximum Charnock parameter
   real(r8), parameter :: fcrit = 22.     ! critical dimensionless fetch for Charnock parameter (Vickers & Mahrt 1997)
   real(r8)  :: cur                       ! Charnock parameter (-)
   real(r8) kva                           ! kinematic viscosity of air at ground temperature and forcing pressure
   real(r8) sqre0                         ! root of roughness Reynolds number
!================================================================================
      cur = cur0 + curm * exp( max( -(fetch*grav/ur/ur)**(1./3.)/fcrit, & ! Fetch-limited
                     -(depth*grav)**0.5/ur ) )   ! depth-limited

      IF ((t_grnd > tfrz .and. snl == 0)) THEN
         kva = kva0 * (t_grnd/293.15)**1.5 * 1.013e5/forc_psrf ! kinematic viscosity of air
         z0mg = max(cus*kva/max(ustar,1.e-4),cur*ustar*ustar/grav) ! momentum roughness length
         z0mg = max(z0mg, 1.0e-5) ! This limit is redundant with current values.
         sqre0 = (max(z0mg*ustar/kva,0.1))**0.5   ! square root of roughness Reynolds number
         z0hg = z0mg * exp( -vonkar/prn*( 4.*sqre0 - 3.2) ) ! SH roughness length
         z0qg = z0mg * exp( -vonkar/sch*( 4.*sqre0 - 4.2) ) ! LH roughness length
         z0qg = max(z0qg, 1.0e-5)  ! Minimum allowed roughness length for unfrozen lakes
         z0hg = max(z0hg, 1.0e-5)  ! set low so it is only to avoid floating point exceptions
      ELSE IF (snl == 0) THEN ! frozen lake with ice, and no snow cover
         z0mg = 0.001              ! z0mg won't have changed
         z0hg = z0mg/exp(0.13 * (ustar*z0mg/1.5e-5)**0.45)
         z0qg = z0hg
      else                          ! USE roughness over snow
         z0mg = 0.0024             ! z0mg won't have changed
         z0hg = z0mg/exp(0.13 * (ustar*z0mg/1.5e-5)**0.45)
         z0qg = z0hg
      ENDIF

   END SUBROUTINE roughness_lake





   function cal_elnth(elopt, hobl, zi, z0m, dplak) result(elnth)
! ---------------------------------- code history -------------------------------------
! Description:
!     Update the turbulence length scale
!     According to Sun et al. (2019) and Elliott and Venayagamoorthy (2011), six options were implemented
!     They are:
!         1. Noh et al. (2011)
!         2. Axell and Liungman（2001）
!
! Original author:
!     Tiejun Ling, 2016
!
! Revisions: 
!     Lei Sun, 2019: Modified for lake
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024 : Separate the function from the update_eddy_pars SUBROUTINE
!
! -------------------------------------------------------------------------------------
   USE module_sf_lake_table, only: vonkar, EL_Noh, EL_AL
!==============================================================================
!  ------------------------- input variables ---------------------------
        integer, intent(in) :: elopt         ! Turbulence length scale parameterization option
        real(r8), intent(in) :: hobl         ! Ocean mixing layer height [m]
        real(r8), intent(in) :: zi           ! layer interface depth [m]
        real(r8), intent(in) :: z0m          ! roughness length [m]
        real(r8), intent(in) :: dplak        ! depth of the lake [m]
        real(r8) :: elnth                    ! Function return value: Turbulence length scale
        real(r8) :: lg                       ! geometric length scale [m]
!================================================================================
   ! *****************************************************************************
   SELECT CASE (elopt)   !- Select the turbulence length scale parameterization
   ! *****************************************************************************
      !---------------
      CASE (EL_Noh) !- Noh et al. (2011)
      !---------------
         !- Noh et al. 2011 suggested that the turbulence length scale is given by:
         lg = vonkar * (zi + z0m) 
         elnth = lg * hobl / (lg + hobl)

      !---------------
      CASE (EL_AL) !- Axell and Liungman（2001）
      !---------------
            !- Sun et al. 2019 suggested that the turbulence length scale is given by:
            lg = sqrt(1 / (1 / (vonkar**2 * (zi + z0m)**2) + 1 / (vonkar**2 * (dplak + z0m - zi)**2)))
            elnth = vonkar * (zi + z0m) / (1.0 + vonkar * (zi + z0m) / hobl)

      ! *****************************************************************************
      END SELECT
      ! *****************************************************************************
   end function cal_elnth



   function cal_prandt(prtopt, prtin,Ri) result(prandt)
! ---------------------------------- code history -------------------------------------
! Description:
!     Update the turbulent Prandtl number
!     According to Sun et al. (2019) and Elliott and Venayagamoorthy (2011), six options were implemented
!     They are:
!         1. User-defined turbulent Prandtl number
!         2. Mellor and Yamada (1982)
!         3. Kim and Mahrt (1982)
!         4. Peters et al. (1988)
!         5. Venayagamoorthy and Stretch (2010)
!         6. Munk and Anderson (1948)
!
! Original author: 
!     Lei Sun, 2019
!
! Revisions:
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024 : Separate the function from the update_eddy_pars SUBROUTINE
!
!--------------------------------------------------------------------------------------
   USE module_sf_lake_table, only: PRT_USER, PRT_MY, PRT_KM, PRT_PGT, PRT_VS, PRT_MA
!================================================================================
   integer, intent(in) :: prtopt   ! Turbulent Prandtl number option
   real(r8), intent(in) :: prtin   ! User-defined turbulent Prandtl number
   real(r8), intent(in) :: Ri      ! Richardson number
   real(r8) :: prandt              ! Function return value: Turbulent Prandtl number
!  ------------------------- local variables ---------------------------
   integer :: k
   real(r8) :: PrtMY    = 0.8      ! Prandtl number in Neutral condition [0.7~1.0], following Mellor and Yamada 1982
   real(r8) :: PrtKM    = 1        ! Prandtl number in Neutral condition, following Kim and Mahrt 1982
   real(r8) :: PrtVS    = 0.7      ! Prandtl number in Neutral condition, following Venayagamoorthy and Stretch 2010
   real(r8) :: Rfinf    = 1./4.    ! Asymptotic value of the flux Richardson number for large Ri, following Venayagamoorthy and Stretch 2010
   real(r8) :: tauinf   = 1./3.    ! Mixing efficiency, tauinf = Rfinf / (1 - Rfinf), following Venayagamoorthy and Stretch 2010
   real(r8) :: betarho  = 10./3.   ! Empirical constants, following Elliott and Venayagamoorthy (2011)
   real(r8) :: alpharho = -3./2.   ! Empirical constants, following Elliott and Venayagamoorthy (2011)
   real(r8) :: beta     = 10.      ! Empirical constants, following Elliott and Venayagamoorthy (2011)
   real(r8) :: alpha    = -1./2.   ! Empirical constants, following Elliott and Venayagamoorthy (2011)
!================================================================================

      ! *****************************************************************************
      SELECT CASE (prtopt)   !- Select the turbulent Prandtl number parameterization
      ! *****************************************************************************
         !---------------
         CASE (PRT_USER) !- User-defined turbulent Prandtl number
         !---------------
            !- No action is taken, as the user-defined Prandtl number is directly used.
            prandt = prtin

         !---------------
         CASE (PRT_MY) !- Mellor and Yamada (1982)
         !---------------
            prandt = PrtMY

         !---------------
         CASE (PRT_KM) !- Kim and Mahrt (1982)
         !---------------
            IF (Ri > 0) THEN
               prandt = PrtKM * (1. + 15. * Ri * sqrt(1. + 5. * Ri)) / (1. + 10. * Ri / sqrt(1. + 5. * Ri)) !- eq. (10, 12) in Kim and Mahrt (1982)
               ! Two approximate options are mentioned in the paper
               ! prandt = 1. + 3.8 * Ri     !- eq. 8 in Kim and Mahrt (1982)
               ! prandt = 1. + 2.1 * Ri     !- after eq. 12 in Kim and Mahrt (1982)
            else
               prandt = PrtKM
            ENDIF

         !---------------
         CASE (PRT_PGT) !- Peters et al. (1988)                                                                    !- eq. (17, 18) in Elliott and Venayagamoorthy (2011)
         !---------------
            ! This model of the turbulent Prandtl number behaves differently from others,
            ! approaching a maximum value of approximately 20 as the Richardson number increases to infinity, 
            ! rather than continuing to grow. It also yields a value of 0 when the Richardson number is 0, 
            ! indicating infinite turbulent mixing in unstratified flows, and is only applicable for positive Richardson numbers.
            IF (Ri > 0) THEN
               IF (Ri .le. 0.25) THEN
                  prandt = 56. / 3. * Ri**1.4
               else
                  prandt = (5 / (sqrt(1 + 5 * Ri)**3) + 0.2) / (5 / (sqrt(1 + 5 * Ri)**5) + 0.01)
               ENDIF
            else
               prandt = PrtVS  !- WMEJ add this line for the case of Ri < 0
            ENDIF

         !---------------
         CASE (PRT_VS) !- Venayagamoorthy and Stretch (2010)
         !---------------
            ! It is questionable whether the second item has PrtVS [-> Ri / (Rfinf * PrtVS)],
            ! as the descriptions in Venayagamoorthy and Stretch (2010) and Elliott and Venayagamoorthy (2011) are different.
            IF (Ri > 0) THEN
               prandt = PrtVS * exp(-1.*Ri / (PrtVS * tauinf)) + Ri / (Rfinf * PrtVS)  !- eq. (3.6) in Venayagamoorthy and Stretch (2010)
            else
               prandt = PrtVS
            ENDIF

         !---------------
         CASE (PRT_MA) !- Munk and Anderson (1948)
         !---------------
            ! Munk and Anderson did not explicitly give the neutral Prandtl number for this scheme.
            ! we USE the same value as Venayagamoorthy and Stretch (2010) for the neutral Prandtl number same as Elliott and Venayagamoorthy (2011)
            IF (Ri > 0) THEN
               prandt = PrtVS * (1 + beta * Ri)**alpha / (1 + betarho * Ri)**alpharho  !- eq. (10) in Elliott and Venayagamoorthy (2011)
            else
               prandt = PrtVS
            ENDIF

         !---------------
         CASE DEFAULT    !- Mellor and Yamada (1982) (default)
         !---------------
            print *, 'Fatal: The default turbulent Prandtl number is used. [ = 0.8, Mellor and Yamada 1982]'
            prandt = PrtMY

      ! *****************************************************************************
      END SELECT
      ! *****************************************************************************

   end function cal_prandt



   function cal_ri(riopt, n2, s2, elnth, qcl) result(Ri)
! ---------------------------------- code history -------------------------------------
! Description:
!     Calculate the turbulent Richardson number.
!     there are two options for the Richardson number parameterization:
!         1. Lewis Fry Richardson : Ri = n^2 / s^2
!         2. Noh et al. (2011)    : Ri =  (nl / q)^2
!
! Original author:
!     Lei Sun, 2019
!
! Revisions:
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024 : Separate the function from the update_eddy_pars SUBROUTINE
! -------------------------------------------------------------------------------------
   USE module_sf_lake_table, only: Richardson, NohRi
!================================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in) :: riopt   ! Richardson number parameterization option
   real(r8), intent(in) :: n2     ! Brunt-Vaisala frequency [/s^2]
   real(r8), intent(in) :: s2     ! Shear number
   real(r8), intent(in) :: elnth  ! turbulence length scale [m]
   real(r8), intent(in) :: qcl    ! rms velocity of turbulence [=sqrt(2*tke) [ms-1]]

   ! ----- output value -----
   real(r8) :: Ri                 ! Function result: Turbulent Richardson number
!================================================================================
      ! *****************************************************************************
      SELECT CASE (riopt)   !- Select the turbulent Richardson number parameterization
      ! *****************************************************************************
         !---------------
         CASE (Richardson) !- Lewis Fry Richardson
         !---------------
            Ri = n2 / s2

         !---------------
         CASE (NohRi) !- Noh et al. (2011)
         !---------------
            Ri = n2 * elnth**2 / qcl**2

         !---------------
         CASE DEFAULT    !- Lewis Fry Richardson (default)
         !---------------
            print *, 'Warning: The default turbulent Richardson number is used. [Ri = n2 / s2]'
            Ri = n2 / s2
      ! *****************************************************************************
      END SELECT
      ! *****************************************************************************

   end function cal_ri



   function buoyeval(nlayer, lktmp) result(buo)
! ---------------------------------- code history -------------------------------------
! Description : 
!     Calculate the buoyancy evaluation
!
! Original author:
!    Tiejun Ling, 2016
!
! Revisions:
!    Lei Sun, 2019: Modified for lake
!    Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024: Reorganized the code to adapt to the new LakeDriver
! -------------------------------------------------------------------------------------
   USE module_sf_lake_table, only: grav
!================================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in)  :: nlayer               ! Maximum number of layer (lake or ocean)
   real(r8), intent(in) :: lktmp(nlayer)        ! lake temperature [K]
!  ------------------------- output variables ---------------------------
   real(r8)             :: buo(nlayer)          ! deltarho/lkrho mean buoyancy [ms-2]
   ! ----- local variables -----
   real(r8)             :: tcl_bkg(nlayer) 
   integer :: k
!================================================================================

      DO k = 1, nlayer
         IF (k == 1) THEN
               tcl_bkg(k) = (lktmp(k) + lktmp(k+1)) / 2.0_r8
         elseif (k == nlayer) THEN
               tcl_bkg(k) = (lktmp(k) + lktmp(k-1)) / 2.0_r8
         else
               tcl_bkg(k) = (lktmp(k-1) + lktmp(k) + lktmp(k+1)) / 3.0_r8
         ENDIF
      ENDDO

      buo = grav * (lktmp / tcl_bkg - 1.0_r8)

   end function buoyeval



   SUBROUTINE trid(m, l1, w1, lm, wm, b, d, a, x, p)
! ---------------------------------- code history -------------------------------------
! Description : tridiagnoal maxtrix solve
!
! Original author:
!    Tiejun Ling, 2016
! -------------------------------------------------------------------------------------
!==============================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in)      :: &
      m                       ,&! dimension of maxtrix
      l1                      ,&! low boundary option
      lm                        ! up boundary option

   real(r8), intent(in)     :: &
      w1                      ,&! low boundary values
      wm                        ! up boundary values

!  ------------------------- inout variables ---------------------------
   real(r8), intent(inout)  :: &
      b(m)                    ,&! low coefs in tri-diagnol maxtrix but the sign is reverse
      d(m)                    ,&! diagnol coefs
      a(m)                    ,&! up coefs but the sign is reverse
      x(m)                      ! RHS vector
!  ------------------------- output variables ---------------------------
   real(r8), intent(out)    :: &
      p(m)                      ! solved vector

!  ------------------------- local variables ---------------------------
   real(r8), dimension(m)  ::  e, f
   integer  :: i, mm
   real(r8) :: den
!================================================================================
      SELECT CASE(l1)
         CASE(1)     ! p(1) = w1
            e(1) = 0.0
            f(1) = w1
         CASE(2)     ! p(1) = p(2)- w1  xxx -k(p(2)-p(1))/dzlak=wx
            e(1) = 1.0
            !f(1) = -w1    ! x switch the sign
            !f(1) = w1 ! ltj comment 20120710
            f(1)= -w1
         CASE(3)
            e(1) = a(1) / d(1)
            f(1) = x(1) / d(1)
      END SELECT

      mm = m-1
      DO i=2, mm
         den  = d(i) - b(i) * e(i-1)
         e(i) = a(i) / den
         f(i) = (x(i) + b(i) * f(i-1)) / den
      ENDDO

      SELECT CASE(lm)
         CASE(1)                          ! fixed value
            p(m) = wm
         CASE(2)
            p(m) = (f(mm) + wm) / (1.0 - e(mm))  ! p(m) = p(mm-1) + wm
         CASE(3)                          ! xxx add the flux in top level in mixing layer
            den  = d(m) - b(m) * e(m-1)
            e(m) = 0.
            f(m) = (x(m) + b(m) * f(m-1)) / den
            p(m) = f(m)
      END SELECT

      DO i=mm, 1, -1
         p(i) = e(i)*p(i+1) + f(i)
      ENDDO

      return
   END SUBROUTINE trid



   real(8) function nhf(mz, k, dz, var)
!-------------------------------------------------------------------------------
! var_{k+1/2}, nhf: next half step
! Description:
!     Calculate the next half step value for a given variable and depth layer.
! Inputs:
!     mz     - total number of layers
!     k      - current layer index
!     dz     - layer thicknesses
!     var    - variable array (e.g., temperature, density)
! Output:
!     nhf    - interpolated value at the next half-step
!
! Original author:
!     Tiejun Ling, 2016
!-------------------------------------------------------------------------------
      ! ----- input variables -----
      integer, intent(in) :: mz, k                   ! mz: total number of layers, k: current index
      real(8), dimension(mz), intent(in) :: dz    ! layer thicknesses
      real(8), dimension(mz), intent(in) :: var      ! variable values (e.g., temperature, density)

      ! ----- function computation -----
      nhf = (var(k+1) * dz(k) + var(k) * dz(k+1)) / (dz(k) + dz(k+1))

   end function nhf


   real(8) function phf(mz, k, dz, var)
!-------------------------------------------------------------------------------
! var_{k-1/2}, phf: previous half step
! Description:
!     Calculate the previous half step value for a given variable and depth layer.
! Inputs:
!     mz     - total number of layers
!     k      - current layer index
!     dz  - layer thicknesses
!     var    - variable array (e.g., temperature, density)
! Output:
!     phf    - interpolated value at the previous half-step
!
! Original author:
!     Tiejun Ling, 2016
!-------------------------------------------------------------------------------
      ! ----- input variables -----
      integer, intent(in) :: mz, k                   ! mz: total number of layers, k: current index
      real(8), dimension(mz), intent(in) :: dz    ! layer thicknesses
      real(8), dimension(mz), intent(in) :: var      ! variable values (e.g., temperature, density)

      ! ----- function computation -----
      phf = (var(k-1) * dz(k) + var(k) * dz(k-1)) / (dz(k) + dz(k-1))
   end function phf



   SUBROUTINE moninobuk_xoml (hw,ht,hq,displa,z0m,z0h,z0q,obu,um,&
                     ustar,temp1,temp2,temp12m,temp22m,fq10m,fm10m,fm,fh,fq)
!============================================================================
!
! Calculation of friction velocity, relation for potential temperatur
! and humidity profiles of surface boundary layer.
! the scheme is based on the work of Zeng et al. (1998):
! Intercomparison of bulk aerodynamic algorithms for the computation
! of sea surface fluxes using TOGA CORE and TAO data. J. Climate, Vol. 11: 2628-2644
!
! Original author: Yongjiu Dai, September 15, 1999
!
!============================================================================
   USE module_sf_lake_table, ONLY: vonkar,one3rd
!  ------------------------- input variables ---------------------------
   real(r8), intent(in) :: hw        ! observational height of wind [m]
   real(r8), intent(in) :: ht        ! observational height of temperature [m]
   real(r8), intent(in) :: hq        ! observational height of humidity [m]
   real(r8), intent(in) :: displa    ! displacement height [m]
   real(r8), intent(in) :: z0m       ! roughness length, momentum [m]
   real(r8), intent(in) :: z0h       ! roughness length, sensible heat [m]
   real(r8), intent(in) :: z0q       ! roughness length, latent heat [m]
   real(r8), intent(in) :: obu       ! monin-obukhov length [m]
   real(r8), intent(in) :: um        ! wind speed including the stablity effect [m/s]

!  ------------------------- output variables ---------------------------
   real(r8), intent(out) :: ustar    ! friction velocity [m/s]
   real(r8), intent(out) :: temp1    ! relation for potential temperature profile
   real(r8), intent(out) :: temp2    ! relation for specific humidity profile
   real(r8), intent(out) :: temp12m  ! relation for temperature at 2m
   real(r8), intent(out) :: temp22m  ! relation for specific humidity at 2m
   real(r8), intent(out) :: fm10m    ! integral of profile function for momentum at 10m
   real(r8), intent(out) :: fq10m    ! integral of profile function for moisture at 10m  +WMEJ add for CWRF
   real(r8), intent(out) :: fm       ! integral of profile function for momentum
   real(r8), intent(out) :: fh       ! integral of profile function for heat
   real(r8), intent(out) :: fq       ! integral of profile function for moisture

!  ------------------------- local variables ---------------------------
   real(r8) zldis  ! reference height minus zero-displacement height [m]
   ! real(r8) psi    ! stability function for unstable case
   real(r8) zetam  ! transition point of flux-gradient relation (wind profile)
   real(r8) zetat  ! transition point of flux-gradient relation (temp. profile)
   real(r8) zeta   ! dimensionless height used in Monin-Obukhov theory
   real(r8) f2     ! integral of profile function at 2m
!================================================================================
   ! adjustment factors for unstable (moz < 0) or stable (moz > 0) conditions.

      ! wind profile
      zldis = hw-displa
      zeta  = zldis/obu
      zetam = 1.574
      IF (zeta < -zetam) THEN
         fm = log(-zetam*obu/z0m) - psi(1,-zetam) &
            + psi(1,z0m/obu) + 1.14*((-zeta)**one3rd-(zetam)**one3rd)
      elseif (zeta < 0.) THEN
         fm = log(zldis/z0m) - psi(1,zeta) + psi(1,z0m/obu)
      elseif (zeta <= 1.) THEN
         fm = log(zldis/z0m) + 5.*zeta - 5.*z0m/obu
      else
         fm = log(obu/z0m) + 5. - 5.*z0m/obu + (5.*log(zeta)+zeta-1.)
      ENDIF
      ustar = vonkar*um / fm

      ! wind at 10m
      zldis = 10.+z0m
      zeta  = zldis/obu
      zetam = 1.574
      IF (zeta < -zetam) THEN
         fm10m = log(-zetam*obu/z0m) - psi(1,-zetam) &
               + psi(1,z0m/obu) + 1.14*((-zeta)**one3rd-(zetam)**one3rd)
      elseif (zeta < 0.) THEN
         fm10m = log(zldis/z0m) - psi(1,zeta) + psi(1,z0m/obu)
      elseif (zeta <= 1.) THEN
         fm10m = log(zldis/z0m) + 5.*zeta - 5.*z0m/obu
      else
         fm10m = log(obu/z0m) + 5. - 5.*z0m/obu + (5.*log(zeta)+zeta-1.)
      ENDIF

      ! temperature profile
      zldis = ht-displa
      zeta  = zldis/obu
      zetat = 0.465
      IF (zeta < -zetat) THEN
         fh = log(-zetat*obu/z0h)-psi(2,-zetat) &
            + psi(2,z0h/obu) + 0.8*((zetat)**(-one3rd)-(-zeta)**(-one3rd))
      elseif (zeta < 0.) THEN
         fh = log(zldis/z0h) - psi(2,zeta) + psi(2,z0h/obu)
      elseif (zeta <= 1.) THEN
         fh = log(zldis/z0h) + 5.*zeta - 5.*z0h/obu
      else
         fh = log(obu/z0h) + 5. - 5.*z0h/obu + (5.*log(zeta)+zeta-1.)
      ENDIF
      temp1 = vonkar / fh

      ! for 2 meter screen temperature
      zldis = 2.+z0h  ! ht-displa
      zeta  = zldis/obu
      zetat = 0.465
      IF (zeta < -zetat) THEN
         f2 = log(-zetat*obu/z0h)-psi(2,-zetat) &
            + psi(2,z0h/obu) + 0.8*((zetat)**(-one3rd)-(-zeta)**(-one3rd))
      elseif (zeta < 0.) THEN
         f2 = log(zldis/z0h) - psi(2,zeta) + psi(2,z0h/obu)
      elseif (zeta <= 1.) THEN
         f2 = log(zldis/z0h) + 5.*zeta - 5.*z0h/obu
      else
         f2 = log(obu/z0h) + 5. - 5.*z0h/obu + (5.*log(zeta)+zeta-1.)
      ENDIF
      temp12m = vonkar / f2

      ! humidity profile
      zldis = hq-displa
      zeta  = zldis/obu
      zetat = 0.465
      IF (zeta < -zetat) THEN
         fq = log(-zetat*obu/z0q) - psi(2,-zetat) &
            + psi(2,z0q/obu) + 0.8*((zetat)**(-one3rd)-(-zeta)**(-one3rd))
      elseif (zeta < 0.) THEN
         fq = log(zldis/z0q) - psi(2,zeta) + psi(2,z0q/obu)
      elseif (zeta <= 1.) THEN
         fq = log(zldis/z0q) + 5.*zeta - 5.*z0q/obu
      else
         fq = log(obu/z0q) + 5. - 5.*z0q/obu + (5.*log(zeta)+zeta-1.)
      ENDIF
      temp2 = vonkar / fq

      ! for 2 meter screen humidity
      zldis = 2.+z0h
      zeta  = zldis/obu
      zetat = 0.465
      IF (zeta < -zetat) THEN
         f2 = log(-zetat*obu/z0q) - psi(2,-zetat) &
            + psi(2,z0q/obu) + 0.8*((zetat)**(-one3rd)-(-zeta)**(-one3rd))
      elseif (zeta < 0.) THEN
         f2 = log(zldis/z0q) - psi(2,zeta) + psi(2,z0q/obu)
      elseif (zeta <= 1.) THEN
         f2 = log(zldis/z0q) + 5.*zeta - 5.*z0q/obu
      else
         f2 = log(obu/z0q) + 5. - 5.*z0q/obu + (5.*log(zeta)+zeta-1.)
      ENDIF
      temp22m = vonkar / f2

      ! for 10 meter humidity
      zldis = 10.+z0h
      zeta  = zldis/obu
      zetat = 0.465
      IF (zeta < -zetat) THEN
         fq10m = log(-zetat*obu/z0q) - psi(2,-zetat) &
               + psi(2,z0q/obu) + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333))
      elseif (zeta < 0.) THEN
         fq10m = log(zldis/z0q) - psi(2,zeta) + psi(2,z0q/obu)
      elseif (zeta <= 1.) THEN
         fq10m = log(zldis/z0q) + 5.*zeta - 5.*z0q/obu
      else
         fq10m = log(obu/z0q) + 5. - 5.*z0q/obu + (5.*log(zeta)+zeta-1.)
      ENDIF

   END SUBROUTINE moninobuk_xoml



   SUBROUTINE moninobukini_xoml (ur,th,thm,thv,dth,dqh,dthv,zldis,z0m,um,obu)
!============================================================================
!
! Initialzation of Monin-Obukhov length,
! the scheme is based on the work of Zeng et al. (1998):
! Intercomparison of bulk aerodynamic algorithms for the computation
! of sea surface fluxes using TOGA CORE and TAO data. J. Climate, Vol. 11: 2628-2644
!
! Original author: Yongjiu Dai, September 15, 1999
!
!============================================================================
   USE module_sf_lake_table, ONLY: grav,vonkar
!  ------------------------- input variables ---------------------------
   real(r8), intent(in) :: ur    ! wind speed at reference height [m/s]
   real(r8), intent(in) :: thm   ! intermediate variable (tm+0.0098*ht)
   real(r8), intent(in) :: th    ! potential temperature [K]
   real(r8), intent(in) :: thv   ! virtual potential temperature [K]
   real(r8), intent(in) :: dth   ! diff of virtual temp. between ref. height and surface
   real(r8), intent(in) :: dthv  ! diff of vir. poten. temp. between ref. height and surface
   real(r8), intent(in) :: dqh   ! diff of humidity between ref. height and surface
   real(r8), intent(in) :: zldis ! reference height minus zero-displacement height [m]
   real(r8), intent(in) :: z0m   ! roughness length, momentum [m]

!  ------------------------- output variables ---------------------------
   real(r8), intent(out) :: um   ! wind speed including the stablity effect [m/s]
   real(r8), intent(out) :: obu  ! monin-obukhov length [m]

!  ------------------------- local variables ---------------------------
   real(r8) wc     ! convective velocity [m/s]
   real(r8) rib    ! bulk Richardson number
   real(r8) zeta   ! dimensionless height used in Monin-Obukhov theory
   real(r8) ustar  ! friction velocity [m/s]
!================================================================================
      ! Initial values of u* and convective velocity

      ustar = 0.06
      wc = 0.5
      IF (dthv >= 0.) THEN
         um = max(ur,0.1)
      else
         um = sqrt(ur*ur+wc*wc)
      ENDIF

      rib = grav*zldis*dthv/(thv*um*um)

      IF (rib >= 0.) THEN      ! neutral or stable
         zeta = rib*log(zldis/z0m)/(1.-5.*min(rib,0.19))
         zeta = min(+2.00,max(+1.e-6,zeta))
      else                   ! unstable
         zeta = rib*log(zldis/z0m)
         zeta = max(-100.,min(-1.e-6,zeta))
      ENDIF
      obu = zldis/zeta

   END SUBROUTINE moninobukini_xoml



   FUNCTION psi(k,zeta)
!============================================================================
   ! stability function for rib < 0
   integer k
   real(r8) psi   ! stability function for unstable case
   real(r8) zeta  ! dimensionless height used in Monin-Obukhov theory
   real(r8) chik  !
!================================================================================
      chik = (1.-16.*zeta)**0.25
      IF (k == 1) THEN
            psi = 2.*log((1.+chik)*0.5)+log((1.+chik*chik)*0.5)-2.*atan(chik)+2.*atan(1.)
      else
            psi = 2.*log((1.+chik*chik)*0.5)
      ENDIF

   END FUNCTION psi


END MODULE MOD_Lake_XOML