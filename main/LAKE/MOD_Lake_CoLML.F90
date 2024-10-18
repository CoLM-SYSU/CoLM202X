#include <define.h>

MODULE  MOD_Lake_CoLML
!-----------------------------------------------------------------------
! DESCRIPTION:
!     Simulating energy balance processes of land water body
!
! REFERENCE:
!     Dai et al, 2018, The lake scheme of the common land model and its performance evaluation.
!     Chinese Science Bulletin, 63(28-29), 3002–3021, https://doi.org/10.1360/N972018-00609
!
! Original author: 
!     Yongjiu Dai 04/2014/
!
! Revisions:
!     Nan Wei,  01/2018: interaction btw prec and lake surface including phase change of prec and water body
!     Nan Wei,  06/2018: update heat conductivity of water body and soil below and snow hydrology
!     Hua Yuan, 01/2023: added snow layer absorption in melting calculation
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024: revised the code structure
!-----------------------------------------------------------------------

   USE MOD_Precision
   IMPLICIT NONE
   SAVE

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: Lake_CoLML
   PUBLIC :: CoLML_SurfFlux
   PUBLIC :: roughness_lake

   ! PRIVATE MEMBER FUNCTIONS:
   PRIVATE :: newsnow_lake
   PRIVATE :: CoLMLTem
   PRIVATE :: snowwater_lake
   PRIVATE :: hConductivity_lake
   PRIVATE :: tridia


!+WMEJ TODO: Adding layered extinction coefficients in the future
!+WMEJ TODO: CoLM may have problems with radiation transmission inside the lake. 
!            The possible reason is that the maximum incident depth of solar radiation is set.
!+WMEJ TODO: USE table parameters to define lake type (Shallow, Deep)
!+WMEJ TODO: Delete the unnecessary RETURN value of the SUBROUTINE 

!--------------------------------------------------------------------------
CONTAINS
!--------------------------------------------------------------------------

   SUBROUTINE Lake_CoLML ( &
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
            porsl     , dksatu    , dkdry     , k_solids  ,&
            vf_quartz , vf_gravels, vf_om     , vf_sand   ,&
            wf_gravels, wf_sand   , csol      , dksatf    ,&
            BA_alpha  , BA_beta   , zcbcv     , tprec     ,&
            ipatch    , bifall    ,&
            ! "inout" arguments
            ! ---------------------------
            dzlak     , zlake     , zilak     , lktmp     ,&
            dzssb     , zssb      , zissb     , t_ssb     ,&
            snlay     , scv       , snwag     , snwdp     ,&
            stke1     , tskin     , xwliq     , xwice     ,&
            z0m       , z0h       , z0q       , felak     ,&
            gamma     , etal      , btpri     , frlak     ,&
            icefr     , tmice     , snout     , lkrho    ,& 
! SNICAR model variables
            forc_aer  , sabg_lyr  , snofrz    ,&
            mss_bcpho , mss_bcphi , mss_ocpho , mss_ocphi ,&
            mss_dst1  , mss_dst2  , mss_dst3  , mss_dst4  ,&
! END SNICAR model variables
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
!     Interface of Lake_CoLML. Used to initialize some of Lake_CoLML's variables 
!     and perform calculations of lake snow cover, lake flux, lake temperature, etc.
!
! Called:
!    -> newsnow_lake             : Add new snow nodes and interaction btw prec and lake surface 
!                                      including phase change of prec and water body
!    -> CoLML_SurfFlux           : Lake surface flux calculations
!    -> CoLMLTem                  : Lake temperature and snow on frozen lake
!    -> snowwater_lake           : Calculation of Lake Hydrology. 
!
! Original author: 
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024
!
! -------------------------------------------------------------------------------------
   USE MOD_Lake_Const, only : tfrz, d2r
   USE MOD_Namelist, only: DEF_USE_CBL_HEIGHT, DEF_USE_SNICAR
!================================================================================
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
      zcbcv                   ,&! convective boundary height [m]
      hpbl                    ,&! atmospheric boundary layer height [m]
      crain                   ,&! convective rainfall [kg/(m2 s)]
      csnow                   ,&! convective snowfall [kg/(m2 s)]
      tprec                   ,&! snowfall/rainfall temperature [K]
      bifall                  ,&! bulk density of newly fallen dry snow [kg/m3]
      lrain                   ,&! large scale rainfall [kg/(m2 s)]
      lsnow                     ! large scale snowfall [kg/(m2 s)]
      
   real(r8), intent(in)     :: &
      vf_quartz (nsoil)       ,&! volumetric fraction of quartz within mineral soil               
      vf_gravels(nsoil)       ,&! volumetric fraction of gravels
      vf_om     (nsoil)       ,&! volumetric fraction of organic matter
      vf_sand   (nsoil)       ,&! volumetric fraction of sand
      wf_gravels(nsoil)       ,&! gravimetric fraction of gravels
      wf_sand   (nsoil)       ,&! gravimetric fraction of sand
      porsl     (nsoil)       ,&! volumetric soil water at saturation (porosity)
      csol      (nsoil)       ,&! heat capacity of soil solids [J/(m3 K)]
      k_solids  (nsoil)       ,&! thermal conductivity of minerals soil [W/m-K]
      dksatf    (nsoil)       ,&! thermal conductivity of saturated frozen soil [W/m-K]
      dkdry     (nsoil)       ,&! thermal conductivity for dry soil  [J/(K s m)]
      dksatu    (nsoil)       ,&! thermal conductivity of saturated unfrozen soil [W/m-K]
      BA_alpha  (nsoil)       ,&! alpha in Balland and Arp(2005) thermal conductivity scheme
      BA_beta   (nsoil)         ! beta in Balland and Arp(2005) thermal conductivity scheme

!  ------------------------- inout variables ---------------------------
   integer, intent(inout)   :: &
      snlay                     ! number of snow layers

   real(r8), intent(inout)  :: &
      dzlak(nlake)            ,&! lake layer thickness [m]
      zlake(nlake)            ,&! lake layer depth [m]
      zilak(nlake+1)          ,&! lake layer interface level [m]
      lktmp(nlake)            ,&! lake temperature [K]
      lkrho(nlake)            ,&! water density [kg/m3]
      icefr(nlake)              ! lake ice fraction [-]
   
   real(r8), intent(inout)  :: &
      dzssb(-nsnow+1:nsoil)   ,&! soil + snow layer thickness [m]
      zssb (-nsnow+1:nsoil)   ,&! soil + snow layer depth [m]
      zissb(-nsnow:nsoil+1)   ,&! soil + snow layer interface level [m]
      t_ssb(-nsnow+1:nsoil)   ,&! soil + snow layer temperature [K]
      xwliq(-nsnow+1:nsoil)   ,&! soil + snow layer liquid water [kg/m2]
      xwice(-nsnow+1:nsoil)     ! soil + snow layer ice lens [kg/m2]
   
   real(r8), intent(inout)  :: &
      scv                     ,&! snow mass [kg/m2]
      snwag                   ,&! non dimensional snow age [-]
      snwdp                   ,&! snow depth [m]
      snout                   ,&! rate of water out of snow bottom [mm/s]
      stke1                   ,&! top level eddy conductivity [W/m/K]
      tskin                   ,&! ground surface temperature [k]
      tmice                   ,&! ice temperature [K]
      z0m                     ,&! roughness length over ground, momentum [m]
      z0h                     ,&! roughness length over ground, sensible heat [m]
      z0q                     ,&! roughness length over ground, latent heat [m]
      felak                   ,&! lake fetch (>0), otherwise USE emperical eq by dplak [m]
      gamma                   ,&! Mixing enhancement factorfor eddy diffusion coefficient [-]
      etal                    ,&! extinction coefficient [1/m]
      btpri                   ,&! beta prime in Monin-Obukhov theory [-]
      frlak                     ! lake fraction [-]
        
   ! ----- SNICAR variables -----
   real(r8), intent(in)     :: &
      forc_aer(14)              ! aerosol deposition from atmosphere model (grd,aer) [kg m-1 s-1]

   real(r8), intent(in)     :: &
      sabg_lyr(-nsnow+1:1)      ! solar radiation absorbed by ground [W/m2]
   
   real(r8), intent(inout)  :: &
      mss_bcpho(-nsnow+1:0)   ,&! mass of hydrophobic BC in snow  (col,lyr) [kg]
      mss_bcphi(-nsnow+1:0)   ,&! mass of hydrophillic BC in snow (col,lyr) [kg]
      mss_ocpho(-nsnow+1:0)   ,&! mass of hydrophobic OC in snow  (col,lyr) [kg]
      mss_ocphi(-nsnow+1:0)   ,&! mass of hydrophillic OC in snow (col,lyr) [kg]
      mss_dst1 (-nsnow+1:0)   ,&! mass of dust species 1 in snow  (col,lyr) [kg]
      mss_dst2 (-nsnow+1:0)   ,&! mass of dust species 2 in snow  (col,lyr) [kg]
      mss_dst3 (-nsnow+1:0)   ,&! mass of dust species 3 in snow  (col,lyr) [kg]
      mss_dst4 (-nsnow+1:0)     ! mass of dust species 4 in snow  (col,lyr) [kg]
   
   real(r8), intent(out)    :: &
      snofrz(-nsnow+1:0)        ! snow freezing flag (col,lyr) [0/1]

!  ------------------------- output variables ---------------------------
   real(r8), intent(out)    :: &
      fsena                   ,&! sensible heat from canopy height to atmosphere [W/m2]
      fevpa                   ,&! evapotranspiration from canopy height to atmosphere [mm/s]
      lfevpa                  ,&! latent heat flux from canopy height to atmosphere [W/2]
      fseng                   ,&! sensible heat flux from ground [W/m2]
      fevpg                   ,&! evaporation heat flux from ground [mm/s]
      olrg                    ,&! outgoing long-wave radiation from ground+canopy [ W/m2]
      fgrnd                   ,&! ground heat flux [W/m2]
      trad                    ,&! radiative temperature [K]
      qseva                   ,&! ground surface evaporation rate (mm h2o/s)
      qsdew                   ,&! ground surface dew formation (mm h2o /s) [+]
      qsubl                   ,&! sublimation rate from snow pack (mm h2o /s) [+]
      qfros                   ,&! surface dew added to snow pack (mm h2o /s) [+]
      taux                    ,&! wind stress: E-W [kg/m/s**2]
      tauy                    ,&! wind stress: N-S [kg/m/s**2]
      ustar                   ,&! u* in similarity theory [m/s]
      qstar                   ,&! q* in similarity theory [kg/kg]
      tstar                   ,&! t* in similarity theory [K]
      emis                    ,&! averaged bulk surface emissivity [-]
      zol                     ,&! dimensionless height (z/L) used in Monin-Obukhov theory [-]
      snwml                   ,&! snowmelt rate [mm/s]
      rib                     ,&! bulk Richardson number in surface layer [-]
      ram                     ,&! aerodynamical resistance [s/m]
      rah                     ,&! thermal resistance [s/m]
      raw                     ,&! moisture resistance [s/m]
      wdm                     ,&! lowest level wind speed including the stablity effect [m/s]
      t2m                     ,&! 2 m height air temperature [K]
      q2m                     ,&! 2 m height air specific humidity [kg/kg]
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
      j                       ,&! loop counter
      i                       ,&! loop counter
      nsl                     ,&! number of soil layers
      lb                      ,&! number of snow layers
      imelt(-nsnow+1:nsoil)     ! flag for: melting=1, freezing=2, Nothing happended=0

   real(r8)                 :: &
      rlat                    ,&! latitude in radians
      htvp                    ,&! latent heat of vapor of water (or sublimation) [j/kg]
      fgrnd1                  ,&! ground heat flux [W/m2]
      rhosnow                 ,&! snow density [kg/m3]
      scvold                  ,&! snow cover for previous time step [mm]
      w_old                   ,&! liquid water mass of the column at the previous time step (mm)
      frain                   ,&! rainfall onto ground including canopy runoff [kg/(m2 s)]
      fsnow                   ,&! snowfall onto ground including canopy runoff [kg/(m2 s)]
      u2m                     ,&! 2 m height wind speed [m/s]
      a                       ,&! fraction of snowfall that is new snow
      aa                      ,&! fraction of snowfall that is new snow
      fiold(-nsnow+1:nsoil)     ! fraction of ice relative to the total water
!================================================================================

      ! --- initialize variables ---
      rlat = xlat * d2r

      scvold = scv          !snow mass at previous time step
      fiold(:) = 0.0
      IF (snlay < 0) THEN
         fiold(snlay+1:0) = xwice(snlay+1:0) / (xwliq(snlay+1:0) + xwice(snlay+1:0))
      ENDIF

      w_old = sum(xwliq(1:)) + sum(xwice(1:))

      frain = crain + lrain
      fsnow   = csnow + lsnow
        
      CALL newsnow_lake ( &
               ! "in" arguments
               ! ---------------
               -nsnow        , nlake         , dtlak        , dzlak        ,&
               frain         , fsnow         , tprec        , bifall       ,&
               ! "inout" arguments
               ! ------------------
               lktmp         , zissb(:0)     , zssb(:0)     , dzssb(:0)    ,&
               t_ssb(:0)     , xwliq(:0)     , xwice(:0)    , fiold(:0)    ,&
               snlay         , snwag         , scv          , snwdp        ,&
               icefr         )  

      ! define snow layer on ice lake
      snlay = 0
      DO j = -nsnow+1, 0
         IF(xwliq(j) + xwice(j) > 0.) snlay = snlay-1
      ENDDO 
      lb = snlay+1

      !- Lake surface flux calculations
      CALL CoLML_SurfFlux ( &
               ! "in" arguments
               ! -------------------
               scwat         , snlay         , zopt          , betaopt     ,&
               fetchopt      , rlat          , dplak         , dtlak       ,& 
               stke1         , hwref         , htref         , hqref       ,&
               usurf         , vsurf         , tmref         , qmref       ,&
               arhos         , psurf         , sols          , soll        ,&
               solsd         , solld         , lwdns         , sabg        ,&
               hpbl          , dzlak(1)      , zlake(nlake)  , dzssb(lb)   ,&
               lktmp(1)      , t_ssb(lb)     , xwliq(lb)     , xwice(lb)   ,&
               icefr(1)      , zcbcv         , ipatch        ,&
               ! "inout" arguments
               ! -------------------
               tskin         , z0m           , z0h           , z0q         ,&
               felak         , btpri         ,&
               ! "out" arguments
               ! -------------------
               fseng         , fevpg         , fsena         , fevpa       ,&
               lfevpa        , olrg          , fgrnd         , fgrnd1      ,&
               zol           , rib           , trad          , htvp        ,&
               emis          , wdm           , ram           , rah         ,&
               raw           , shfdt         , taux          , tauy        ,&
               t2m           , q2m           , u10m          , v10m        ,&
               fh2m          , fq2m          , fm10m         , fq10m       ,&
               fm            , fh            , fq            , ustar       ,&
               qstar         , tstar         , rhosnow       , u2m          )

      !- Lake temperature
      CALL CoLMLTem ( &
               ! "in" arguments
               ! -------------------
               scwat         , -nsnow        , nsoil         , nlake       ,&
               rlat          , dtlak         , dplak         ,&
               hwref         , htref         , hqref         , hpbl        ,&
               usurf         , vsurf         , qmref         , tmref       ,&
               arhos         , psurf         , lwdns         , sabg        ,&
               sols          , soll          , solsd         , solld       ,&
               vf_quartz     , vf_gravels    , vf_om         , vf_sand     ,&
               wf_gravels    , wf_sand       , porsl         , csol        ,&
               k_solids      , dksatu        , dksatf        , dkdry       ,&
               BA_alpha      , BA_beta       ,& 
               zopt          , betaopt       , fetchopt      , etaopt      ,&      
               ! "inout" arguments
               ! -------------------
               tskin         , scv           , icefr         , snwdp       ,&
               dzlak         , zlake         , zilak         , lktmp       ,&
               dzssb         , zssb          , zissb         , t_ssb       ,&
               z0m           , z0h           , z0q           , stke1       ,&
               gamma         , btpri         , felak         , etal        ,&
               xwliq         , xwice         , imelt         , snlay       ,&
               fseng         , fevpg         , fsena         , olrg        ,&
               lfevpa        , fgrnd         , fgrnd1        , htvp        ,&
               ustar         , rhosnow       , lkrho         ,&
               ! SNICAR model variables
               snofrz        , sabg_lyr      ,&
               ! END SNICAR model variables
               ! "out" arguments
               ! -------------------
               qseva         , qsubl         , qsdew         , qfros       ,&
               snwml         )

      !- Lake Hydrology
      CALL snowwater_lake ( &
               ! "in"  arguments
               ! ---------------------------
               -nsnow        , nsoil         , nlake         , dtlak       ,&
               qseva         , qsubl         , qsdew         , qfros       ,&
               porsl         , frain         , fsnow         , dzlak       ,&
               imelt(-nsnow+1:0), fiold(-nsnow+1:0),&
               ! "inout"  arguments
               ! ---------------------------
               zssb          , dzssb         , zissb         , t_ssb       ,&
               xwice         , xwliq         , lktmp         , icefr       ,&
               fseng         , fgrnd         , snlay         , scv         ,&
               snout         ,&
               ! SNICAR model variables
               forc_aer      ,&
               mss_bcpho     , mss_bcphi     , mss_ocpho     , mss_ocphi   ,&
               mss_dst1      , mss_dst2      , mss_dst3      , mss_dst4    ,&
               ! END SNICAR model variables
               snwdp         , snwml         , usurf         , vsurf       )

      !-WMEJ: The variable saturated flow scheme does not need to be considered at present,
      !       so the code related to this scheme has been deleted.

      ! Set zero to the empty node
      IF (snlay > -nsnow) THEN
         xwice(-nsnow+1:snlay) = 0.
         xwliq(-nsnow+1:snlay) = 0.
         t_ssb(-nsnow+1:snlay) = 0.
         zssb (-nsnow+1:snlay) = 0.
         dzssb(-nsnow+1:snlay) = 0.
      ENDIF
      tmice = tfrz
        
   END SUBROUTINE Lake_CoLML

    

   SUBROUTINE newsnow_lake ( &
            ! "in" arguments
            ! ---------------
            maxsnl    , nl_lake      , deltim      , dz_lake     ,&
            pg_rain   , pg_snow      , t_precip    , bifall      ,&
            ! "inout" arguments
            ! ------------------
            t_lake    , zi_soisno    , z_soisno    ,&
            dz_soisno , t_soisno     , wliq_soisno , wice_soisno ,&
            fiold     , snl          , sag         , scv         ,&
            snowdp    , lake_icefrac )
                            
!-----------------------------------------------------------------------
! DESCRIPTION:
!     Add new snow nodes and interaction btw prec and lake surface 
!     including phase change of prec and water body
!
! Original author : 
!     Yongjiu Dai, 04/2014
!
! Revisions:
!     Nan Wei,  01/2018: update interaction btw prec and lake surface
!-----------------------------------------------------------------------
   USE MOD_Lake_Const, only : tfrz, denh2o, cpliq, cpice, hfus
   USE MOD_Namelist, only: DEF_USE_SNICAR
!================================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in)     :: maxsnl                   ! maximum number of snow layers
   integer, intent(in)     :: nl_lake                  ! number of soil layers
   real(r8), intent(in)    :: deltim                   ! seconds in a time step [second]
   real(r8), intent(in)    :: t_precip                 ! snowfall/rainfall temperature [kelvin]
   real(r8), intent(in)    :: bifall                   ! bulk density of newly fallen dry snow [kg/m3]
   real(r8), intent(in)    :: dz_lake(1:nl_lake)       ! lake layer thickness (m)
   
!  ------------------------- inout variables ---------------------------
   real(r8), intent(inout) :: pg_rain                  ! liquid water onto ground [kg/(m2 s)]
   real(r8), intent(inout) :: pg_snow                  ! ice onto ground [kg/(m2 s)]
   real(r8), intent(inout) :: zi_soisno(maxsnl:0)      ! interface level below a "z" level (m)
   real(r8), intent(inout) :: z_soisno(maxsnl+1:0)     ! snow layer depth (m)
   real(r8), intent(inout) :: dz_soisno(maxsnl+1:0)    ! snow layer thickness (m)
   real(r8), intent(inout) :: t_soisno(maxsnl+1:0)     ! snow layer temperature [K]
   real(r8), intent(inout) :: wliq_soisno(maxsnl+1:0)  ! snow layer liquid water (kg/m2)
   real(r8), intent(inout) :: wice_soisno(maxsnl+1:0)  ! snow layer ice lens (kg/m2)
   real(r8), intent(inout) :: fiold(maxsnl+1:0)        ! fraction of ice relative to the total water
   integer, intent(inout)  :: snl                      ! number of snow layers
   real(r8), intent(inout) :: sag                      ! non dimensional snow age [-]
   real(r8), intent(inout) :: scv                      ! snow mass (kg/m2)
   real(r8), intent(inout) :: snowdp                   ! snow depth (m)
   real(r8), intent(inout) :: lake_icefrac(1:nl_lake)  ! mass fraction of lake layer that is frozen
   real(r8), intent(inout) :: t_lake(1:nl_lake)        ! lake layer temperature (m)
   
!  ------------------------- local variables ---------------------------
   integer lb
   integer newnode                                     ! signification when new snow node is set, (1=yes, 0=non)
   real(r8) dz_snowf                                   ! layer thickness rate change due to precipitation [m/s]
   real(r8) a, b, c, d, e, f, g, h
   real(r8) wice_lake(1:nl_lake), wliq_lake(1:nl_lake), tw
!================================================================================

      newnode = 0
      dz_snowf = pg_snow/bifall
      snowdp = snowdp + dz_snowf*deltim
      scv = scv + pg_snow*deltim                 ! snow water equivalent (mm)

      zi_soisno(0) = 0.

      IF (snl==0 .and. snowdp < 0.01) THEN       ! no snow layer, energy exchange between prec and lake surface
         a = cpliq*pg_rain*deltim*(t_precip-tfrz)                          !cool down rainfall to tfrz
         b = pg_rain*deltim*hfus                                           !all rainfall frozen
         c = cpice*denh2o*dz_lake(1)*lake_icefrac(1)*(tfrz-t_lake(1))      !warm up lake surface ice to tfrz
         d = denh2o*dz_lake(1)*lake_icefrac(1)*hfus                        !all lake surface ice melt
         e = cpice*pg_snow*deltim*(tfrz-t_precip)                          !warm up snowfall to tfrz
         f = pg_snow*deltim*hfus                                           !all snowfall melt
         g = cpliq*denh2o*dz_lake(1)*(1-lake_icefrac(1))*(t_lake(1)-tfrz)  !cool down lake surface water to tfrz
         h = denh2o*dz_lake(1)*(1-lake_icefrac(1))*hfus                    !all lake surface water frozen
         sag = 0.0

         IF (lake_icefrac(1) > 0.999) THEN
            ! all rainfall frozen, release heat to warm up frozen lake surface
            IF (a+b<=c) THEN
               tw=min(tfrz,t_precip)
               t_lake(1)=(a+b+cpice*(pg_rain+pg_snow)*deltim*tw+cpice*denh2o*dz_lake(1)*t_lake(1)*lake_icefrac(1))/&
                           (cpice*denh2o*dz_lake(1)*lake_icefrac(1)+cpice*(pg_rain+pg_snow)*deltim)
               scv = scv+pg_rain*deltim
               snowdp = snowdp + pg_rain*deltim/bifall
               pg_snow = pg_snow+pg_rain
               pg_rain = 0.0
            ! prec tem at tfrz, partial rainfall frozen ->release heat -> warm up lake surface to tfrz (no latent heat)
            ELSE IF (a<=c) THEN
               t_lake(1)=tfrz
               scv = scv + (c-a)/hfus
               snowdp = snowdp + (c-a)/(hfus*bifall)
               pg_snow = pg_snow + min(pg_rain,(c-a)/(hfus*deltim))
               pg_rain = max(0.0,pg_rain - (c-a)/(hfus*deltim))
            ! lake surface tem at tfrz, partial lake surface melt -> absorb heat -> cool down rainfall to tfrz (no latent heat)
            ELSE IF (a<=c+d) THEN
               t_lake(1)=tfrz
               wice_lake(1) = denh2o*dz_lake(1) - (a-c)/hfus
               wliq_lake(1) = (a-c)/hfus
               lake_icefrac(1) = wice_lake(1)/(wice_lake(1) + wliq_lake(1))
            ! all lake surface melt, absorb heat to cool down rainfall
            ELSE  !(a>c+d)
               t_lake(1)=(cpliq*pg_rain*deltim*t_precip+cpliq*denh2o*dz_lake(1)*tfrz-c-d)/&
                           (cpliq*denh2o*dz_lake(1)+cpliq*pg_rain*deltim)
               lake_icefrac(1) = 0.0
            ENDIF

            IF (snowdp>=0.01) THEN  !frozen rain may make new snow layer
               snl = -1
               newnode = 1
               dz_soisno(0)  = snowdp             ! meter
               z_soisno (0)  = -0.5*dz_soisno(0)
               zi_soisno(-1) = -dz_soisno(0)
               sag = 0.                           ! snow age

               t_soisno (0) = t_lake(1)           ! K
               wice_soisno(0) = scv               ! kg/m2
               wliq_soisno(0) = 0.                ! kg/m2
               fiold(0) = 1.
            ENDIF

         ELSE IF (lake_icefrac(1) >= 0.001) THEN
            IF (pg_rain > 0.0 .and. pg_snow > 0.0) THEN
               t_lake(1)=tfrz
            ELSE IF (pg_rain > 0.0) THEN
               IF (a>=d) THEN
                  t_lake(1)=(cpliq*pg_rain*deltim*t_precip+cpliq*denh2o*dz_lake(1)*tfrz-d)/&
                              (cpliq*denh2o*dz_lake(1)+cpliq*pg_rain*deltim)
                  lake_icefrac(1) = 0.0
               ELSE
                  t_lake(1)=tfrz
                  wice_lake(1) = denh2o*dz_lake(1)*lake_icefrac(1) - a/hfus
                  wliq_lake(1) = denh2o*dz_lake(1)*(1-lake_icefrac(1)) + a/hfus
                  lake_icefrac(1) = wice_lake(1)/(wice_lake(1) + wliq_lake(1))
               ENDIF
            ELSE IF (pg_snow > 0.0) THEN
               IF (e>=h) THEN
                  t_lake(1)=(h+cpice*denh2o*dz_lake(1)*tfrz+cpice*pg_snow*deltim*t_precip)/&
                              (cpice*pg_snow*deltim+cpice*denh2o*dz_lake(1))
                  lake_icefrac(1) = 1.0
               ELSE
                  t_lake(1)=tfrz
                  wice_lake(1) = denh2o*dz_lake(1)*lake_icefrac(1) + e/hfus
                  wliq_lake(1) = denh2o*dz_lake(1)*(1-lake_icefrac(1)) - e/hfus
                  lake_icefrac(1) = wice_lake(1)/(wice_lake(1) + wliq_lake(1))
               ENDIF
            ENDIF

         ELSE
            ! all snowfall melt, absorb heat to cool down lake surface water
            IF (e+f<=g) THEN
               tw=max(tfrz,t_precip)
               t_lake(1)=(cpliq*denh2o*dz_lake(1)*t_lake(1)*(1-lake_icefrac(1))+cpliq*(pg_rain+pg_snow)*deltim*tw-e-f)/&
                           (cpliq*(pg_rain+pg_snow)*deltim+cpliq*denh2o*dz_lake(1)*(1-lake_icefrac(1)))
               scv = scv - pg_snow*deltim
               snowdp = snowdp - dz_snowf*deltim
               pg_rain = pg_rain + pg_snow
               pg_snow = 0.0
            ! prec tem at tfrz, partial snowfall melt ->absorb heat -> cool down lake surface to tfrz (no latent heat)
            ELSE IF (e<=g) THEN
               t_lake(1) = tfrz
               scv = scv - (g-e)/hfus
               snowdp = snowdp - (g-e)/(hfus*bifall)
               pg_rain = pg_rain + min(pg_snow, (g-e)/(hfus*deltim))
               pg_snow = max(0.0, pg_snow - (g-e)/(hfus*deltim))
            ! lake surface tem at tfrz, partial lake surface frozen -> release heat -> warm up snowfall to tfrz (no latent heat)
            ELSE IF (e<=g+h) THEN
               t_lake(1) = tfrz
               wice_lake(1) = (e-g)/hfus
               wliq_lake(1) = denh2o*dz_lake(1) - (e-g)/hfus
               lake_icefrac(1) = wice_lake(1)/(wice_lake(1) + wliq_lake(1))
            ! all lake surface frozen, release heat to warm up snowfall
            ELSE       !(e>g+h)
               t_lake(1) = (g+h+cpice*denh2o*dz_lake(1)*tfrz+cpice*pg_snow*deltim*t_precip)/&
                           (cpice*pg_snow*deltim+cpice*denh2o*dz_lake(1))
               lake_icefrac(1) = 1.0
            ENDIF
         ENDIF

      ELSE IF (snl==0 .and. snowdp >= 0.01) THEN
         ! only ice part of snowfall is added here, the liquid part will be added later
         snl = -1
         newnode = 1
         dz_soisno(0)  = snowdp             ! meter
         z_soisno (0)  = -0.5*dz_soisno(0)
         zi_soisno(-1) = -dz_soisno(0)
         sag = 0.                           ! snow age

         t_soisno (0) = min(tfrz, t_precip) ! K
         wice_soisno(0) = scv               ! kg/m2
         wliq_soisno(0) = 0.                ! kg/m2
         fiold(0) = 1.

      ELSE                                   ! ( snl<0 .and. newnode ==0 )
         lb = snl + 1
         t_soisno(lb) = ( (wice_soisno(lb)*cpice+wliq_soisno(lb)*cpliq)*t_soisno(lb) &
                     +   (pg_rain*cpliq + pg_snow*cpice)*deltim*t_precip ) &
                     / ( wice_soisno(lb)*cpice + wliq_soisno(lb)*cpliq &
                     +   pg_rain*deltim*cpliq + pg_snow*deltim*cpice )

         t_soisno(lb) = min(tfrz, t_soisno(lb))
         wice_soisno(lb) = wice_soisno(lb)+deltim*pg_snow
         dz_soisno(lb) = dz_soisno(lb)+dz_snowf*deltim
         z_soisno(lb) = zi_soisno(lb) - 0.5*dz_soisno(lb)
         zi_soisno(lb-1) = zi_soisno(lb) - dz_soisno(lb)

      ENDIF

   END SUBROUTINE newsnow_lake



   SUBROUTINE CoLMLTem ( &
            ! "in" arguments
            ! -------------------
            patchtype    , maxsnl      , nl_soil      , nl_lake   ,&
            rlat         , deltim      , lakedepth    ,&
            forc_hgt_u   , forc_hgt_t  , forc_hgt_q   , hpbl      ,&
            forc_us      , forc_vs     , forc_q       , forc_t    ,&
            forc_rhoair  , forc_psrf   , forc_frl     , sabg      ,&
            forc_sols    , forc_soll   , forc_solsd   , forc_solld,&
            vf_quartz    , vf_gravels  , vf_om        , vf_sand   ,&
            wf_gravels   , wf_sand     , porsl        , csol      ,&
            k_solids     , dksatu      , dksatf       , dkdry     ,&
            BA_alpha     , BA_beta     ,&
            zopt         , betaopt     , fetchopt     , etaopt    ,&
            ! "inout" arguments
            ! -------------------
            t_grnd       , scv         , lake_icefrac , snowdp    ,&
            dz_lake      , z_lake      , zi_lake      , t_lake    ,&
            dz_soisno    , z_soisno    , zi_soisno    , t_soisno  ,&
            z0m          , z0h         , z0q          , savedtke1 ,&
            gamma        , btpri       , felak        , etal      ,&
            wliq_soisno  , wice_soisno , imelt_soisno , snl       ,&
            fseng        , fevpg       , fsena        , olrg      ,&
            lfevpa       , fgrnd       , fgrnd1       , htvp      ,&
            ustar        , rhosnow     , lkrho       ,&
! SNICAR model variables
            snofrz       , sabg_lyr    ,& 
! END SNICAR model variables
            ! "out" arguments
            ! -------------------
            qseva       , qsubl        , qsdew        , qfros     ,&
            sm          , urban_call  )

! ------------------------ code history ---------------------------
! purpose: lake temperature and snow on frozen lake
! initial  Yongjiu Dai, 2000
!          Zack Subin, 2009
!          Yongjiu Dai, /12/2012/, /04/2014/, 06/2018
!          Nan Wei, /06/2018/
!
! ------------------------ notes ----------------------------------
! Lakes have variable depth, possible snow layers above, freezing & thawing of lake water,
! and soil layers with active temperature and gas diffusion below.
!
! Calculates temperatures in the 25-30 layer column of (possible) snow,
! lake water, soil, and bedrock beneath lake.
! Snow and soil temperatures are determined as in SoilTemperature, except
! for appropriate boundary conditions at the top of the snow (the flux is fixed
! to be the ground heat flux), the bottom of the snow (adjacent to top lake layer),
! and the top of the soil (adjacent to the bottom lake layer).
! Also, the soil is kept fully saturated.
! The whole column is solved simultaneously as one tridiagonal matrix.
!
! calculate lake temperatures from one-dimensional thermal
! stratification model based on eddy diffusion concepts to
! represent vertical mixing of heat
!
! d ts    d            d ts     1 ds
! ---- = -- [(km + ke) ----] + -- --
!  dt    dz             dz     cw dz
! where: ts = temperature (kelvin)
!         t = time (s)
!         z = depth (m)
!        km = molecular diffusion coefficient (m**2/s)
!        ke = eddy diffusion coefficient (m**2/s)
!        cw = heat capacity (j/m**3/kelvin)
!         s = heat source term (w/m**2)
!
! USE crank-nicholson method to set up tridiagonal system of equations to
! solve for ts at time n+1, where the temperature equation for layer i is
! r_i = a_i [ts_i-1] n+1 + b_i [ts_i] n+1 + c_i [ts_i+1] n+1
! the solution conserves energy as
! cw*([ts(  1)] n+1 - [ts(  1)] n)*dz(  1)/dt + ... +
! cw*([ts(nl_lake)] n+1 - [ts(nl_lake)] n)*dz(nl_lake)/dt = fin
! where
! [ts] n   = old temperature (kelvin)
! [ts] n+1 = new temperature (kelvin)
! fin      = heat flux into lake (w/m**2)
!          = beta*sabg_lyr(1)+forc_frl-olrg-fsena-lfevpa-hm + phi(1) + ... + phi(nl_lake)
!
! REVISIONS:
!     Yongjiu Dai and Hua Yuan, 01/2023: added SNICAR for layer solar absorption, ground heat
!                                        flux, temperature and freezing mass calculations
!     Shaofeng Liu, 05/2023: add option to CALL moninobuk_leddy, the LargeEddy
!                            surface turbulence scheme (LZD2022);
!                            make a proper update of um.
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024 : Moved surface flux calculations to CoLML_SurfFlux 
!                                                            so that other lake models can CALL it
!
! -----------------------------------------------------------------
   USE MOD_Lake_Const, only : tfrz, hvap, hfus, hsub, tkwat, tkice, tkair, stefnc,&
                                    vonkar, grav, cpliq, cpice, cpair, denh2o, denice, rgas, SHALLOW, DEEP
   USE MOD_SoilThermalParameters, only : soil_hcap_cond
   USE MOD_Namelist, only: DEF_USE_SNICAR
!================================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in)     :: patchtype               ! it is lake type, not land type (1: Shallow, 2: Deep)
   integer, intent(in)     :: maxsnl                  ! maximum number of snow layers
   integer, intent(in)     :: nl_soil                 ! number of soil layers
   integer, intent(in)     :: nl_lake                 ! number of lake layers
   real(r8), intent(in)    :: rlat                    ! latitude (radians)
   real(r8), intent(in)    :: deltim                  ! seconds in a time step (s)
   real(r8), intent(in)    :: lakedepth               ! column lake depth (m)
   real(r8), intent(in)    :: forc_hgt_u              ! observational height of wind [m]
   real(r8), intent(in)    :: forc_hgt_t              ! observational height of temperature [m]
   real(r8), intent(in)    :: forc_hgt_q              ! observational height of humidity [m]
   real(r8), intent(in)    :: hpbl                    ! atmospheric boundary layer height [m]
   real(r8), intent(in)    :: forc_us                 ! wind component in eastward direction [m/s]
   real(r8), intent(in)    :: forc_vs                 ! wind component in northward direction [m/s]
   real(r8), intent(in)    :: forc_t                  ! temperature at agcm reference height [kelvin]
   real(r8), intent(in)    :: forc_q                  ! specific humidity at agcm reference height [kg/kg]
   real(r8), intent(in)    :: forc_rhoair             ! density air [kg/m3]
   real(r8), intent(in)    :: forc_psrf               ! atmosphere pressure at the surface [pa]
   real(r8), intent(in)    :: forc_frl                ! atmospheric infrared (longwave) radiation [W/m2]
   real(r8), intent(in)    :: sabg                    ! solar radiation absorbed by ground [W/m2]
   real(r8), intent(in)    :: forc_sols               ! atm vis direct beam solar rad onto srf [W/m2]
   real(r8), intent(in)    :: forc_soll               ! atm nir direct beam solar rad onto srf [W/m2]
   real(r8), intent(in)    :: forc_solsd              ! atm vis diffuse solar rad onto srf [W/m2]
   real(r8), intent(in)    :: forc_solld              ! atm nir diffuse solar rad onto srf [W/m2]
   real(r8), intent(in)    :: vf_quartz (1:nl_soil)   ! volumetric fraction of quartz within mineral soil
   real(r8), intent(in)    :: vf_gravels(1:nl_soil)   ! volumetric fraction of gravels
   real(r8), intent(in)    :: vf_om     (1:nl_soil)   ! volumetric fraction of organic matter
   real(r8), intent(in)    :: vf_sand   (1:nl_soil)   ! volumetric fraction of sand
   real(r8), intent(in)    :: wf_gravels(1:nl_soil)   ! gravimetric fraction of gravels
   real(r8), intent(in)    :: wf_sand   (1:nl_soil)   ! gravimetric fraction of sand
   real(r8), intent(in)    :: porsl     (1:nl_soil)   ! soil porosity [-]
   real(r8), intent(in)    :: csol      (1:nl_soil)   ! heat capacity of soil solids [J/(m3 K)]
   real(r8), intent(in)    :: k_solids  (1:nl_soil)   ! thermal conductivity of mineralssoil [W/m-K]
   real(r8), intent(in)    :: dksatu    (1:nl_soil)   ! thermal conductivity of saturated unfrozen soil [W/m-K]
   real(r8), intent(in)    :: dksatf    (1:nl_soil)   ! thermal conductivity of saturated frozen soil [W/m-K]
   real(r8), intent(in)    :: dkdry     (1:nl_soil)   ! thermal conductivity of dry soil [W/m-K]
   real(r8), intent(in)    :: BA_alpha  (1:nl_soil)   ! alpha in Balland and Arp(2005) thermal conductivity scheme
   real(r8), intent(in)    :: BA_beta   (1:nl_soil)   ! beta in Balland and Arp(2005) thermal conductivity scheme
!-WMEJ add new variables <gamma, btpri, zopt, z0m, z0h, z0q> 
   integer , intent(in)    :: zopt                    ! option for roughness length, 1: constant, 2: Subin et al. (2012)
   integer , intent(in)    :: betaopt                 ! option for betaprime calculation, 1: constant, 2: equation
   integer , intent(in)    :: fetchopt                ! option for fetch length, 1: constant, 2: equation
   integer , intent(in)    :: etaopt                  ! option for Extinction coefficient calculation, 1: constant, 2: equation
   real(r8), intent(inout) :: gamma                   ! Mixing enhancement factorfor eddy diffusion coefficient
   real(r8), intent(inout) :: btpri                   ! Effective beta
   real(r8), intent(inout) :: felak                   ! fetch length of lake [m] 
   real(r8), intent(inout) :: etal                    ! light extinction coefficient [/m]
   real(r8), intent(inout) :: z0m                     ! roughness length over ground, momentum [m]
   real(r8), intent(inout) :: z0h                     ! roughness length over ground, sensible heat [m]
   real(r8), intent(inout) :: z0q                     ! roughness length over ground, latent heat [m]
   real(r8), intent(inout) :: lkrho(nl_lake)          ! density of water [kg/m3]
!-WMEJ END new variables

!  ------------------------- inout variables ---------------------------
   real(r8), intent(inout) :: t_grnd                  ! surface temperature [K]
   real(r8), intent(inout) :: scv                     ! snow water equivalent [mm]
   real(r8), intent(inout) :: lake_icefrac(nl_lake)   ! lake mass fraction of lake layer that is frozen
   real(r8), intent(inout) :: snowdp                  ! snow depth [mm]
   real(r8), intent(inout) :: dz_lake(nl_lake)        ! lake layer thickness (m)
   real(r8), intent(inout) :: z_lake(nl_lake)         ! lake layer depth (m)
   real(r8), intent(inout) :: zi_lake(nl_lake+1)      ! interface level below a "z" level (m)
   real(r8), intent(inout) :: t_lake(nl_lake)         ! lake temperature [K]
   real(r8), intent(inout) :: dz_soisno(maxsnl+1:nl_soil) ! soil/snow layer thickness (m)
   real(r8), intent(inout) :: z_soisno(maxsnl+1:nl_soil)  ! soil/snow layer depth (m)
   real(r8), intent(inout) :: zi_soisno(maxsnl:nl_soil+1) ! interface level below a "z" level (m)
   real(r8), intent(inout) :: t_soisno(maxsnl+1:nl_soil)  ! soil/snow temperature [K]
   real(r8), intent(inout) :: savedtke1               ! top level eddy conductivity (W/m K)
   real(r8), intent(inout) :: wliq_soisno (maxsnl+1:nl_soil) ! soil/snow liquid water (kg/m2)
   real(r8), intent(inout) :: wice_soisno (maxsnl+1:nl_soil) ! soil/snow ice lens (kg/m2)
   integer,  intent(inout) :: imelt_soisno(maxsnl+1:nl_soil) ! soil/snow flag for melting (=1), freezing (=2), Not=0 (new)
   integer,  intent(inout) :: snl                     ! number of snow layers
   real(r8), intent(inout) :: fseng                   ! sensible heat flux from ground [W/m2]
   real(r8), intent(inout) :: fevpg                   ! evaporation heat flux from ground [mm/s]
   real(r8), intent(inout) :: fsena                   ! sensible heat from canopy height to atmosphere [W/m2]
   real(r8), intent(inout) :: lfevpa                  ! latent heat flux from canopy height to atmosphere [W/m2]
   real(r8), intent(inout) :: olrg                    ! outgoing long-wave radiation from ground+canopy
   real(r8), intent(inout) :: fgrnd                   ! ground heat flux [W/m2]
   real(r8), intent(inout) :: ustar                   ! u* in similarity theory [m/s]
   real(r8), intent(inout) :: htvp                    ! latent heat of vapor of water (or sublimation) [j/kg]
   real(r8), intent(inout) :: rhosnow                 ! partitial density of water (ice + liquid)
   real(r8), intent(inout) :: fgrnd1                  ! ground heat flux into the first snow/lake layer [W/m2]
! SNICAR model variables
   real(r8), intent(out)   :: snofrz   (maxsnl+1:0)   ! snow freezing rate (col,lyr) [kg m-2 s-1]
   real(r8), intent(in)    :: sabg_lyr (maxsnl+1:1)   ! solar radiation absorbed by ground [W/m2]
! END SNICAR model variables

!  ------------------------- output variables ---------------------------
   real(r8), intent(out)   :: qseva                   ! ground surface evaporation rate (mm h2o/s)
   real(r8), intent(out)   :: qsubl                   ! sublimation rate from snow pack (mm H2O /s) [+]
   real(r8), intent(out)   :: qsdew                   ! surface dew added to snow pack (mm H2O /s) [+]
   real(r8), intent(out)   :: qfros                   ! ground surface frosting formation (mm H2O /s) [+]
   real(r8), intent(out)   :: sm                      ! rate of snowmelt [mm/s, kg/(m2 s)]

   logical, optional, intent(in) :: urban_call        ! whether it is a urban CALL

!  ------------------------- local variables ---------------------------
   integer  :: idlak                                  ! index of lake, 1 = deep lake, 2 = shallow lake
   real(r8) :: errore                                 ! lake temperature energy conservation error (w/m**2)
   real(r8) :: hm                                     ! energy residual [W/m2]
   real(r8) :: wliq_lake(nl_lake)                     ! lake liquid water (kg/m2)
   real(r8) :: wice_lake(nl_lake)                     ! lake ice lens (kg/m2)
   real(r8) :: vf_water(1:nl_soil)                    ! volumetric fraction liquid water within underlying soil
   real(r8) :: vf_ice(1:nl_soil)                      ! volumetric fraction ice len within underlying soil
   real(r8) :: betaprime                              ! Effective beta
   real(r8) :: cfus                                   ! effective heat of fusion per unit volume
   real(r8) :: tkice_eff                              ! effective conductivity since layer depth is constant
   real(r8) :: cice_eff                               ! effective heat capacity of ice (using density of
                                                      ! water because layer depth is not adjusted when freezing
   real(r8) :: cwat                                   ! specific heat capacity of water (j/m**3/kelvin)
   real(r8) :: u2m                                    ! 2 m height wind speed [m/s]
   real(r8) :: rhow(nl_lake)                          ! density of water (kg/m**3)
   real(r8) :: fin                                    ! heat flux into lake - flux out of lake (w/m**2)
   real(r8) :: phi(nl_lake)                           ! solar radiation absorbed by layer (w/m**2)
   real(r8) :: phi_soil                               ! solar radiation into top soil layer (W/m^2)
   real(r8) :: phidum                                 ! temporary value of phi
   real(r8) :: tktopsoil                              ! thermal conductivity of the top soil layer [W/(m K)]
   integer  :: imelt_lake  (1:nl_lake)                ! lake flag for melting or freezing snow and soil layer [-]
   real(r8) :: cv_lake     (1:nl_lake)                ! heat capacity [J/(m2 K)]
   real(r8) :: t_lake_bef  (1:nl_lake)                ! beginning lake temp for energy conservation check [K]
   real(r8) :: tk_lake     (1:nl_lake)                ! thermal conductivity at layer node [W/(m K)]
   real(r8) :: hcap        (1:nl_soil)                ! J/(m3 K)
   real(r8) :: cv_soisno   (maxsnl+1:nl_soil)         ! heat capacity of soil/snow [J/(m2 K)]
   real(r8) :: tk_soisno   (maxsnl+1:nl_soil)         ! thermal conductivity of soil/snow [W/(m K)] (at interface below, except for j=0)
   real(r8) :: thk         (maxsnl+1:nl_soil)         ! W/(m K)
   real(r8) :: t_soisno_bef(maxsnl+1:nl_soil)         ! beginning soil/snow temp for E cons. check [K]
   real(r8) :: wice_soisno_bef(maxsnl+1:0)            ! ice lens [kg/m2]
   real(r8) :: cvx    (maxsnl+1:nl_lake+nl_soil)      ! heat capacity for whole column [J/(m2 K)]
   real(r8) :: tkix   (maxsnl+1:nl_lake+nl_soil)      ! thermal conductivity at layer interfaces for whole column [W/(m K)]
   real(r8) :: phix   (maxsnl+1:nl_lake+nl_soil)      ! solar source term for whole column [W/m**2]
   real(r8) :: zx     (maxsnl+1:nl_lake+nl_soil)      ! interface depth (+ below surface) for whole column [m]
   real(r8) :: tx     (maxsnl+1:nl_lake+nl_soil)      ! temperature of whole column [K]
   real(r8) :: tx_bef (maxsnl+1:nl_lake+nl_soil)      ! beginning lake/snow/soil temp for energy conservation check [K]
   real(r8) :: factx  (maxsnl+1:nl_lake+nl_soil)      ! coefficient used in computing tridiagonal matrix
   real(r8) :: fnx    (maxsnl+1:nl_lake+nl_soil)      ! heat diffusion through the layer interface below [W/m2]
   real(r8) :: a      (maxsnl+1:nl_lake+nl_soil)      ! "a" vector for tridiagonal matrix
   real(r8) :: b      (maxsnl+1:nl_lake+nl_soil)      ! "b" vector for tridiagonal matrix
   real(r8) :: c      (maxsnl+1:nl_lake+nl_soil)      ! "c" vector for tridiagonal matrix
   real(r8) :: r      (maxsnl+1:nl_lake+nl_soil)      ! "r" vector for tridiagonal solution
   real(r8) :: fn1    (maxsnl+1:nl_lake+nl_soil)      ! heat diffusion through the layer interface below [W/m2]
   real(r8) :: brr    (maxsnl+1:nl_lake+nl_soil)      !
   integer  :: imelt_x(maxsnl+1:nl_lake+nl_soil)      ! flag for melting (=1), freezing (=2), Not=0 (new)
   real(r8) :: dzm                                    ! used in computing tridiagonal matrix [m]
   real(r8) :: dzp                                    ! used in computing tridiagonal matrix [m]
   real(r8) :: zin                                    ! depth at top of layer (m)
   real(r8) :: zout                                   ! depth at bottom of layer (m)
   real(r8) :: rsfin                                  ! relative flux of solar radiation into layer
   real(r8) :: rsfout                                 ! relative flux of solar radiation out of layer
   real(r8) :: eta                                    ! light extinction coefficient (/m): depends on lake type
   real(r8) :: za(2)                                  ! base of surface absorption layer (m): depends on lake type
   ! real(r8) :: sm                                     ! rate of snowmelt [mm/s, kg/(m2 s)]
   real(r8) :: hs                                     ! net ground heat flux into the surface
   real(r8) :: dhsdT                                  ! temperature derivative of "hs"
   real(r8) :: heatavail                              ! available energy for melting or freezing (J/m^2)
   real(r8) :: heatrem                                ! energy residual or loss after melting or freezing
   real(r8) :: melt                                   ! actual melting (+) or freezing (-) [kg/m2]
   real(r8) :: xmf                                    ! total per-column latent heat abs. from phase change  (J/m^2)
   real(r8) :: ocvts                                  ! (cwat*(t_lake[n  ])*dz_lake
   real(r8) :: ncvts                                  ! (cwat*(t_lake[n+1])*dz_lake
   real(r8) :: esum1                                  ! temp for checking energy (J/m^2)
   real(r8) :: esum2                                  ! ""
   real(r8) :: zsum                                   ! temp for putting ice at the top during convection (m)
   real(r8) :: errsoi                                 ! soil/lake energy conservation error (W/m^2)
   real(r8) :: iceav                                  ! used in calc aver ice for convectively mixed layers
   real(r8) :: qav                                    ! used in calc aver heat content for conv. mixed layers
   real(r8) :: tav                                    ! used in aver temp for convectively mixed layers
   real(r8) :: tav_froz                               ! used in aver temp for convectively mixed layers (C)
   real(r8) :: tav_unfr                               ! "
   real(r8) :: nav                                    ! used in aver temp for convectively mixed layers
   real(r8) :: fevpg_lim                              ! temporary evap_soi limited by top snow layer content [mm/s]
   real(r8) :: scv_temp                               ! temporary h2osno [kg/m^2]
   real(r8) :: tmp                                    !
   real(r8) :: h_fin                                  !
   real(r8) :: h_finDT                                !
   integer  :: iter                                   ! iteration index
   integer  :: convernum                              ! number of time when del_T_grnd < 0.01
   integer  :: nl_sls                                 ! abs(snl)+nl_lake+nl_soil
   integer  :: lb                                     ! lower bound of arrays
   integer  :: jprime                                 ! j - nl_lake
   integer  :: i,j                                    ! DO loop or array index

!  ------------------------- constants variables ---------------------------
   integer, parameter  :: itmax  = 40                 ! maximum number of iteration
   integer, parameter  :: itmin  = 6                  ! minimum number of iteration
   real(r8), parameter :: delmax = 3.0                ! maximum change in lake temperature [K]
   real(r8), parameter :: dtmin  = 0.01               ! max limit for temperature convergence [K]
   real(r8), parameter :: dlemin = 0.1                ! max limit for energy flux convergence [w/m2]
   real(r8), parameter :: depthcrit = 25.             ! (m) Depth beneath which to enhance mixing
   real(r8), parameter :: fangmult = 5.               ! Multiplier for unfrozen diffusivity
   real(r8), parameter :: minmultdepth = 20.          ! (m) Minimum depth for imposing fangmult
   real(r8), parameter :: cnfac  = 0.5                ! Crank Nicholson factor between 0 and 1
!================================================================================

! ======================================================================
!*[1] constants and model parameters
! ======================================================================
      ! constants for lake temperature model
      za = (/0.5, 0.6/)     !- old scheme, with problems
      ! za = (/0.1, 0.1/)
      cwat = cpliq*denh2o     ! water heat capacity per unit volume
      cice_eff = cpice*denh2o ! USE water density because layer depth is not adjusted for freezing
      cfus = hfus*denh2o      ! latent heat per unit volume
      tkice_eff = tkice * denice/denh2o ! effective conductivity since layer depth is constant

      ! define snow layer on ice lake
      snl = 0
      DO j=maxsnl+1,0
         IF(wliq_soisno(j)+wice_soisno(j)>0.) snl=snl-1
      ENDDO
      lb = snl+1

      ! Base on lake depth, assuming that small lakes are likely to be shallower
      ! Estimate crudely based on lake depth
      IF (z_lake(nl_lake) < 4.) THEN
         idlak = 1
      ELSE
         idlak = 2
      ENDIF

      !+WMEJ betaprime is calculated in CoLML_SurfFlux
      betaprime = btpri

! ======================================================================
!*[2] Begin stability iteration and fluxes calculation
! ======================================================================
      ! -WMEJ move surface flux calculations to CoLML_SurfFlux

      ! January 12, 2023 by Yongjiu Dai
      IF (DEF_USE_SNICAR .and. .not. present(urban_call)) THEN
         hs = sabg_lyr(lb) + forc_frl - olrg - fseng - htvp*fevpg
         dhsdT = 0.0
      ENDIF

! ======================================================================
!*[3] Begin temperature calculation
! ======================================================================
      !------------------------------------------------------------
      ! Set up vector r and vectors a, b, c that define tridiagonal matrix
      ! snow and lake and soil layer temperature
      !------------------------------------------------------------

      !------------------------------------------------------------
      ! Lake density
      !------------------------------------------------------------

      DO j = 1, nl_lake
         rhow(j) = (1.-lake_icefrac(j))*denh2o*(1.0-1.9549e-05*(abs(t_lake(j)-277.))**1.68) &
                     + lake_icefrac(j)*denice
         ! allow for ice fraction; assume constant ice density.
         ! this is not the correct average-weighting but that's OK because the density will only
         ! be used for convection for lakes with ice, and the ice fraction will dominate the
         ! density differences between layers.
         ! using this average will make sure that surface ice is treated properly during
         ! convective mixing.
      ENDDO

      !------------------------------------------------------------
      ! Diffusivity and implied thermal "conductivity" = diffusivity * cwat
      !------------------------------------------------------------
      DO j = 1, nl_lake
         cv_lake(j) = dz_lake(j) * (cwat*(1.-lake_icefrac(j)) + cice_eff*lake_icefrac(j))
      ENDDO

      CALL hConductivity_lake(nl_lake,snl,t_grnd,&
                              z_lake,t_lake,lake_icefrac,rhow,&
                              rlat,ustar,z0m,lakedepth,depthcrit,gamma,tk_lake,savedtke1)

      !------------------------------------------------------------
      ! Set the thermal properties of the snow above frozen lake and underlying soil
      ! and check initial energy content.
      !------------------------------------------------------------
      lb = snl+1
      DO i = 1, nl_soil
         vf_water(i) = wliq_soisno(i)/(dz_soisno(i)*denh2o)
         vf_ice(i) = wice_soisno(i)/(dz_soisno(i)*denice)
         CALL soil_hcap_cond(vf_gravels(i),vf_om(i),vf_sand(i),porsl(i),&
                              wf_gravels(i),wf_sand(i),k_solids(i),&
                              csol(i),dkdry(i),dksatu(i),dksatf(i),&
                              BA_alpha(i),BA_beta(i),&
                              t_soisno(i),vf_water(i),vf_ice(i),hcap(i),thk(i))
         cv_soisno(i) = hcap(i)*dz_soisno(i)
      ENDDO

      ! Snow heat capacity and conductivity
      IF(lb <=0 )THEN
         DO j = lb, 0
               cv_soisno(j) = cpliq*wliq_soisno(j) + cpice*wice_soisno(j)
               rhosnow = (wice_soisno(j)+wliq_soisno(j))/dz_soisno(j)
               thk(j) = tkair + (7.75e-5*rhosnow + 1.105e-6*rhosnow*rhosnow)*(tkice-tkair)
         ENDDO
      ENDIF

      ! Thermal conductivity at the layer interface
      DO i = lb, nl_soil-1

      ! the following consideration is try to avoid the snow conductivity
      ! to be dominant in the thermal conductivity of the interface.
      ! Because when the distance of bottom snow node to the interfacee
      ! is larger than that of interface to top soil node,
      ! the snow thermal conductivity will be dominant, and the result is that
      ! lees heat tranfer between snow and soil

      ! modified by Nan Wei, 08/25/2014
         IF (i /= 0) THEN
               tk_soisno(i) = thk(i)*thk(i+1)*(z_soisno(i+1)-z_soisno(i)) &
                  /(thk(i)*(z_soisno(i+1)-zi_soisno(i))+thk(i+1)*(zi_soisno(i)-z_soisno(i)))
         ELSE
               tk_soisno(i) = thk(i)
         ENDIF
      ENDDO
      tk_soisno(nl_soil) = 0.
      tktopsoil = thk(1)

      ! Sum cv_lake*t_lake for energy check
      ! Include latent heat term, and USE tfrz as reference temperature
      ! to prevent abrupt change in heat content due to changing heat capacity with phase change.

      ! This will need to be over all soil / lake / snow layers. Lake is below.
      ocvts = 0.
      DO j = 1, nl_lake
         ocvts = ocvts + cv_lake(j)*(t_lake(j)-tfrz) + cfus*dz_lake(j)*(1.-lake_icefrac(j))
      ENDDO

      ! Now DO for soil / snow layers
      DO j = lb, nl_soil
         ocvts = ocvts + cv_soisno(j)*(t_soisno(j)-tfrz) + hfus*wliq_soisno(j)
         IF (j == 1 .and. scv > 0. .and. j == lb) THEN
               ocvts = ocvts - scv*hfus
         ENDIF
      ENDDO

      ! Set up solar source terms (phix)
      ! Modified January 12, 2023 by Yongjiu Dai
      IF (.not. DEF_USE_SNICAR .or. present(urban_call)) THEN
         IF ((t_grnd > tfrz .and. t_lake(1) > tfrz .and. snl == 0)) THEN      !no snow cover, unfrozen layer lakes
            DO j = 1, nl_lake
               ! extinction coefficient from surface data (1/m), IF no eta from surface data,
               ! set eta, the extinction coefficient, according to L Hakanson, Aquatic Sciences, 1995
               ! (regression of secchi depth with lake depth for small glacial basin lakes), and the
               ! Poole & Atkins expression for extinction coeffient of 1.7 / secchi Depth (m).

               IF (etaopt == 1) THEN
                  eta = etal
               ELSE
                  eta = 1.1925*max(lakedepth,1.)**(-0.424)
               ENDIF
               zin  = z_lake(j) - 0.5*dz_lake(j)
               zout = z_lake(j) + 0.5*dz_lake(j)
               rsfin  = exp( -eta*max(  zin-za(idlak),0. ) )  ! the radiation within surface layer (z<za)
               rsfout = exp( -eta*max( zout-za(idlak),0. ) )  ! is considered fixed at (1-beta)*sabg
                                                              ! i.e, max(z-za, 0)
               ! Let rsfout for bottom layer go into soil.
               ! This looks like it should be robust even for pathological cases,
               ! like lakes thinner than za(idlak).

               phi(j) = (rsfin-rsfout) * sabg * (1.-betaprime)
               IF (j == nl_lake) phi_soil = rsfout * sabg * (1.-betaprime)
            ENDDO
         ELSE IF (snl == 0) THEN     !no snow-covered layers, but partially frozen
            phi(1) = sabg * (1.-betaprime)
            phi(2:nl_lake) = 0.
            phi_soil = 0.
         ELSE   ! snow covered, this should be improved upon; Mironov 2002 suggests that SW can penetrate thin ice and may
            ! cause spring convection.
            phi(:) = 0.
            phi_soil = 0.
         ENDIF
      ELSE
         DO j = 1, nl_lake
            ! extinction coefficient from surface data (1/m), IF no eta from surface data,
            ! set eta, the extinction coefficient, according to L Hakanson, Aquatic Sciences, 1995
            ! (regression of secchi depth with lake depth for small glacial basin lakes), and the
            ! Poole & Atkins expression for extinction coeffient of 1.7 / secchi Depth (m).
            IF (etaopt == 1) THEN
               eta = etal
            ELSE
               eta = 1.1925*max(lakedepth,1.)**(-0.424)
            ENDIF
            zin  = z_lake(j) - 0.5*dz_lake(j)
            zout = z_lake(j) + 0.5*dz_lake(j)
            rsfin  = exp( -eta*max(  zin-za(idlak),0. ) )  ! the radiation within surface layer (z<za)
            rsfout = exp( -eta*max( zout-za(idlak),0. ) )  ! is considered fixed at (1-beta)*sabg
                                                           ! i.e, max(z-za, 0)
            ! Let rsfout for bottom layer go into soil.
            ! This looks like it should be robust even for pathological cases,
            ! like lakes thinner than za(idlak).

            phi(j) = (rsfin-rsfout) * sabg_lyr(1) * (1.-betaprime)
            IF (j == nl_lake) phi_soil = rsfout * sabg_lyr(1) * (1.-betaprime)
         ENDDO
      ENDIF

      phix(:) = 0.
      phix(1:nl_lake) = phi(1:nl_lake)         !lake layer
      phix(nl_lake+1) = phi_soil               !top soil layer

      ! Set up interface depths(zx), and temperatures (tx).
      DO j = lb, nl_lake+nl_soil
         jprime = j - nl_lake
         IF (j <= 0) THEN                      !snow layer
               zx(j) = z_soisno(j)
               tx(j) = t_soisno(j)
         ELSE IF (j <= nl_lake) THEN           !lake layer
               zx(j) = z_lake(j)
               tx(j) = t_lake(j)
         ELSE                                  !soil layer
               zx(j) = z_lake(nl_lake) + dz_lake(nl_lake)/2. + z_soisno(jprime)
               tx(j) = t_soisno(jprime)
         ENDIF
      ENDDO

      tx_bef = tx

      ! Heat capacity and resistance of snow without snow layers (<1cm) is ignored during diffusion,
      ! but its capacity to absorb latent heat may be used during phase change.

      DO j = lb, nl_lake+nl_soil
         jprime = j - nl_lake

         ! heat capacity [J/(m2 K)]
         IF (j <= 0) THEN                      !snow layer
            cvx(j) = cv_soisno(j)
         ELSE IF (j <= nl_lake) THEN           !lake layer
            cvx(j) = cv_lake(j)
         ELSE                                  !soil layer
            cvx(j) = cv_soisno(jprime)
         ENDIF

         ! Determine interface thermal conductivities at layer interfaces [W/(m K)]
         IF (j < 0) THEN                       !non-bottom snow layer
            tkix(j) = tk_soisno(j)
         ELSE IF (j == 0) THEN                 !bottom snow layer
            dzp = zx(j+1) - zx(j)
            tkix(j) = tk_lake(1)*tk_soisno(j)*dzp &
                    /(tk_soisno(j)*z_lake(1) + tk_lake(1)*(-zx(j)))
                           ! tk_soisno(0) is the conductivity at the middle of that layer
         ELSE IF (j < nl_lake) THEN            !non-bottom lake layer
            tkix(j) = (tk_lake(j)*tk_lake(j+1) * (dz_lake(j+1)+dz_lake(j))) &
                    / (tk_lake(j)*dz_lake(j+1) + tk_lake(j+1)*dz_lake(j))
         ELSE IF (j == nl_lake) THEN           !bottom lake layer
            dzp = zx(j+1) - zx(j)
            tkix(j) = (tktopsoil*tk_lake(j)*dzp &
                    / (tktopsoil*dz_lake(j)/2. + tk_lake(j)*z_soisno(1)))
         ELSE !soil layer
            tkix(j) = tk_soisno(jprime)
         ENDIF

      ENDDO

      ! Determine heat diffusion through the layer interface and factor used in computing
      ! tridiagonal matrix and set up vector r and vectors a, b, c that define tridiagonal
      ! matrix and solve system
      DO j = lb, nl_lake+nl_soil
         factx(j) = deltim/cvx(j)
         IF (j < nl_lake+nl_soil) THEN         !top or interior layer
            fnx(j) = tkix(j)*(tx(j+1)-tx(j))/(zx(j+1)-zx(j))
         ELSE                                  !bottom soil layer
            fnx(j) = 0. !not used
         ENDIF
      ENDDO

      IF (DEF_USE_SNICAR .and. .not. present(urban_call)) THEN
         IF (lb <= 0) THEN                        ! snow covered
            DO j = lb, 1
               IF (j == lb) THEN                  ! top snow layer
                  dzp  = zx(j+1)-zx(j)
                  a(j) = 0.0
                  b(j) = 1. + (1.-cnfac)*factx(j)*tkix(j)/dzp
                  c(j) = -(1.-cnfac)*factx(j)* tkix(j)/dzp
                  r(j) = tx_bef(j)+ factx(j)*(hs - dhsdT*tx_bef(j) + cnfac*fnx(j))
               ELSE IF (j <= 0) THEN              ! non-top snow layers
                  dzm  = (zx(j)-zx(j-1))
                  dzp  = (zx(j+1)-zx(j))
                  a(j) =   - (1.-cnfac)*factx(j)* tkix(j-1)/dzm
                  b(j) = 1.+ (1.-cnfac)*factx(j)*(tkix(j)/dzp + tkix(j-1)/dzm)
                  c(j) =   - (1.-cnfac)*factx(j)* tkix(j)/dzp
                  r(j) = tx_bef(j) + cnfac*factx(j)*(fnx(j) - fnx(j-1)) + factx(j)*sabg_lyr(j)
               ELSE                               ! snow covered top lake layer
                  dzm  = (zx(j)-zx(j-1))
                  dzp  = (zx(j+1)-zx(j))
                  a(j) =   - (1.-cnfac)*factx(j)* tkix(j-1)/dzm
                  b(j) = 1.+ (1.-cnfac)*factx(j)*(tkix(j)/dzp + tkix(j-1)/dzm)
                  c(j) =   - (1.-cnfac)*factx(j)* tkix(j)/dzp
                  r(j) = tx_bef(j) + cnfac*factx(j)*(fnx(j) - fnx(j-1)) + factx(j)*(phix(j) + betaprime*sabg_lyr(j))
               ENDIF
            ENDDO
         ELSE
            j = 1                                 ! no snow covered top lake layer
            dzp  = zx(j+1)-zx(j)
            a(j) = 0.0
            b(j) = 1. + (1.-cnfac)*factx(j)*tkix(j)/dzp
            c(j) = -(1.-cnfac)*factx(j)* tkix(j)/dzp
            r(j) = tx_bef(j)+ factx(j)*(cnfac*fnx(j)+ phix(j)+fgrnd1)
         ENDIF

         DO j = 2, nl_lake+nl_soil
               IF (j < nl_lake+nl_soil) THEN         ! middle lake and soil layers
               dzm  = (zx(j)-zx(j-1))
               dzp  = (zx(j+1)-zx(j))
               a(j) =   - (1.-cnfac)*factx(j)* tkix(j-1)/dzm
               b(j) = 1.+ (1.-cnfac)*factx(j)*(tkix(j)/dzp + tkix(j-1)/dzm)
               c(j) =   - (1.-cnfac)*factx(j)* tkix(j)/dzp
               r(j) = tx_bef(j) + cnfac*factx(j)*(fnx(j) - fnx(j-1)) + factx(j)*phix(j)
               ELSE                                  ! bottom soil layer
               dzm  = (zx(j)-zx(j-1))
               a(j) =   - (1.-cnfac)*factx(j)*tkix(j-1)/dzm
               b(j) = 1.+ (1.-cnfac)*factx(j)*tkix(j-1)/dzm
               c(j) = 0.
               r(j) = tx_bef(j) - cnfac*factx(j)*fnx(j-1)
               ENDIF
         ENDDO
      ! January 12, 2023
      ELSE
         DO j = lb, nl_lake+nl_soil
               IF (j == lb) THEN                     ! top layer
                  dzp  = zx(j+1)-zx(j)
                  a(j) = 0.0
                  b(j) = 1. + (1.-cnfac)*factx(j)*tkix(j)/dzp
                  c(j) = -(1.-cnfac)*factx(j)* tkix(j)/dzp
                  r(j) = tx_bef(j)+ factx(j)*(cnfac*fnx(j)+ phix(j)+fgrnd1)
               ELSE IF (j < nl_lake+nl_soil) THEN    ! middle layer
                  dzm  = (zx(j)-zx(j-1))
                  dzp  = (zx(j+1)-zx(j))
                  a(j) =   - (1.-cnfac)*factx(j)* tkix(j-1)/dzm
                  b(j) = 1.+ (1.-cnfac)*factx(j)*(tkix(j)/dzp + tkix(j-1)/dzm)
                  c(j) =   - (1.-cnfac)*factx(j)* tkix(j)/dzp
                  r(j) = tx_bef(j) + cnfac*factx(j)*(fnx(j) - fnx(j-1)) + factx(j)*phix(j)
               ELSE                                  ! bottom soil layer
                  dzm  = (zx(j)-zx(j-1))
                  a(j) =   - (1.-cnfac)*factx(j)*tkix(j-1)/dzm
                  b(j) = 1.+ (1.-cnfac)*factx(j)*tkix(j-1)/dzm
                  c(j) = 0.
                  r(j) = tx_bef(j) - cnfac*factx(j)*fnx(j-1)
               ENDIF
         ENDDO
      ENDIF

      !------------------------------------------------------------
      ! Solve for tdsolution
      !------------------------------------------------------------

      nl_sls = abs(snl) + nl_lake + nl_soil

      CALL tridia (nl_sls, a(lb:), b(lb:), c(lb:), r(lb:), tx(lb:))

      DO j = lb, nl_lake + nl_soil
         jprime = j - nl_lake
         IF (j < 1) THEN               ! snow layer
               t_soisno(j) = tx(j)
         ELSE IF (j <= nl_lake) THEN   ! lake layer
               t_lake(j) = tx(j)
         ELSE                          ! soil layer
               t_soisno(jprime) = tx(j)
         ENDIF
      ENDDO

      ! calculate sublimation, frosting, dewing
      qseva = 0.
      qsubl = 0.
      qsdew = 0.
      qfros = 0.
      IF (fevpg >= 0.0) THEN
         IF(lb < 0)THEN
            qseva = min(wliq_soisno(lb)/deltim, fevpg)
            qsubl = fevpg - qseva
         ELSE
            qseva = min((1.-lake_icefrac(1))*1000.*dz_lake(1)/deltim, fevpg)
            qsubl = fevpg - qseva
         ENDIF
      ELSE
         IF (t_grnd < tfrz) THEN
            qfros = abs(fevpg)
         ELSE
            qsdew = abs(fevpg)
         ENDIF
      ENDIF


#if(defined CoLMDEBUG)
      ! sum energy content and total energy into lake for energy check. 
      !     any errors will be from the tridiagonal solution.
      esum1 = 0.0
      esum2 = 0.0
      DO j = lb, nl_lake + nl_soil
         esum1 = esum1 + (tx(j)-tx_bef(j))*cvx(j)
         esum2 = esum2 + (tx(j)-tfrz)*cvx(j)
      ENDDO
      ! fgrnd includes all the solar radiation absorbed in the lake,
      errsoi = esum1/deltim - fgrnd
      IF(abs(errsoi) > 0.1) THEN
         write(6,*)'energy conservation error in LAND WATER COLUMN during tridiagonal solution,', &
                  'error (W/m^2):', errsoi, fgrnd
      ENDIF
#endif

!------------------------------------------------------------
!*[4] Phase change
!------------------------------------------------------------
      sm = 0.0
      xmf = 0.0
      imelt_soisno(:) = 0
      imelt_lake(:) = 0

      IF (DEF_USE_SNICAR .and. .not. present(urban_call)) THEN
         wice_soisno_bef(lb:0) = wice_soisno(lb:0)
      ENDIF

      ! Check for case of snow without snow layers and top lake layer temp above freezing.
      IF (snl == 0 .and. scv > 0. .and. t_lake(1) > tfrz) THEN
         heatavail = (t_lake(1) - tfrz) * cv_lake(1)
         melt = min(scv, heatavail/hfus)
         heatrem = max(heatavail - melt*hfus, 0.) !catch small negative value to keep t at tfrz
         t_lake(1) = tfrz + heatrem/(cv_lake(1))

         snowdp = max(0., snowdp*(1. - melt/scv))
         scv = scv - melt

         IF (scv < 1.e-12) scv = 0.        ! prevent tiny residuals
         IF (snowdp < 1.e-12) snowdp = 0.  ! prevent tiny residuals
         sm = sm + melt/deltim
         xmf = xmf + melt*hfus
      ENDIF

      ! Lake phase change
      DO j = 1,nl_lake
         IF (t_lake(j) > tfrz .and. lake_icefrac(j) > 0.) THEN ! melting
            imelt_lake(j) = 1
            heatavail = (t_lake(j) - tfrz) * cv_lake(j)
            melt = min(lake_icefrac(j)*denh2o*dz_lake(j), heatavail/hfus)
                  !denh2o is used because layer thickness is not adjusted for freezing
            heatrem = max(heatavail - melt*hfus, 0.)  !catch small negative value to keep t at tfrz
         ELSE IF (t_lake(j) < tfrz .and. lake_icefrac(j) < 1.) THEN !freezing
            imelt_lake(j) = 2
            heatavail = (t_lake(j) - tfrz) * cv_lake(j)
            melt = max(-(1.-lake_icefrac(j))*denh2o*dz_lake(j), heatavail/hfus)
                  !denh2o is used because layer thickness is not adjusted for freezing
            heatrem = min(heatavail - melt*hfus, 0.)  !catch small positive value to keep t at tfrz
         ENDIF
         ! Update temperature and ice fraction.
         IF (imelt_lake(j) > 0) THEN
            lake_icefrac(j) = lake_icefrac(j) - melt/(denh2o*dz_lake(j))
            IF (lake_icefrac(j) > 1.-1.e-12) lake_icefrac(j) = 1.  ! prevent tiny residuals
            IF (lake_icefrac(j) < 1.e-12)    lake_icefrac(j) = 0.  ! prevent tiny residuals
            cv_lake(j) = cv_lake(j) + melt*(cpliq-cpice)           ! update heat capacity
            t_lake(j) = tfrz + heatrem/cv_lake(j)
            xmf = xmf + melt*hfus
         ENDIF
      ENDDO

      ! snow & soil phase change. currently, does not DO freezing point depression.
      DO j = snl+1,nl_soil
         IF (t_soisno(j) > tfrz .and. wice_soisno(j) > 0.) THEN ! melting
            imelt_soisno(j) = 1
            heatavail = (t_soisno(j) - tfrz) * cv_soisno(j)
            melt = min(wice_soisno(j), heatavail/hfus)
            heatrem = max(heatavail - melt*hfus, 0.) !catch small negative value to keep t at tfrz
            IF (j <= 0) sm = sm + melt/deltim
         ELSE IF (t_soisno(j) < tfrz .and. wliq_soisno(j) > 0.) THEN !freezing
            imelt_soisno(j) = 2
            heatavail = (t_soisno(j) - tfrz) * cv_soisno(j)
            melt = max(-wliq_soisno(j), heatavail/hfus)
            heatrem = min(heatavail - melt*hfus, 0.) !catch small positive value to keep t at tfrz
         ENDIF

         ! Update temperature and soil components.
         IF (imelt_soisno(j) > 0) THEN
            wice_soisno(j) = wice_soisno(j) - melt
            wliq_soisno(j) = wliq_soisno(j) + melt
            IF (wice_soisno(j) < 1.e-12) wice_soisno(j) = 0. ! prevent tiny residuals
            IF (wliq_soisno(j) < 1.e-12) wliq_soisno(j) = 0. ! prevent tiny residuals
            cv_soisno(j) = cv_soisno(j) + melt*(cpliq-cpice) ! update heat capacity
            t_soisno(j) = tfrz + heatrem/cv_soisno(j)
            xmf = xmf + melt*hfus
         ENDIF
      ENDDO
      !------------------------------------------------------------

      IF (DEF_USE_SNICAR .and. .not. present(urban_call)) THEN
         !for SNICAR: layer freezing mass flux (positive):
         DO j = lb, 0
            IF (imelt_soisno(j)==2 .and. j<1) THEN
               snofrz(j) = max(0._r8,(wice_soisno(j)-wice_soisno_bef(j)))/deltim
            ENDIF
         ENDDO
      ENDIF

#if(defined CoLMDEBUG)
      ! second energy check and water check. now check energy balance before and after phase
      ! change, considering the possibility of changed heat capacity during phase change, by
      ! using initial heat capacity in the first step, final heat capacity in the second step,
      ! and differences from tfrz only to avoid enthalpy correction for (cpliq-cpice)*melt*tfrz.
      ! also check soil water sum.
      DO j = 1, nl_lake
         esum2 = esum2 - (t_lake(j)-tfrz)*cv_lake(j)
      ENDDO

      DO j = lb, nl_soil
         esum2 = esum2 - (t_soisno(j)-tfrz)*cv_soisno(j)
      ENDDO

      esum2 = esum2 - xmf
      errsoi = esum2/deltim

      IF(abs(errsoi) > 0.1) THEN
         write(6,*) 'energy conservation error in LAND WATER COLUMN during phase change, error (W/m^2):', errsoi
      ENDIF
#endif

      !------------------------------------------------------------
      !*[5] Convective mixing: make sure fracice*dz is conserved, heat content c*dz*T is conserved, and
      ! all ice ends up at the top. Done over all lakes even IF frozen.
      ! Either an unstable density profile or ice in a layer below an incompletely frozen layer will trigger.
      !------------------------------------------------------------
      ! recalculate density
      DO j = 1, nl_lake
         rhow(j) = (1.-lake_icefrac(j))*1000.*(1.0-1.9549e-05*(abs(t_lake(j)-277.))**1.68) &
                     + lake_icefrac(j)*denice
      ENDDO

      DO j = 1, nl_lake-1
         qav = 0.
         nav = 0.
         iceav = 0.

         IF (rhow(j)>rhow(j+1) .or. (lake_icefrac(j)<1.0 .and. lake_icefrac(j+1)>0.)) THEN
            DO i = 1, j+1
               qav = qav + dz_lake(i)*(t_lake(i)-tfrz) * &
                        ((1. - lake_icefrac(i))*cwat + lake_icefrac(i)*cice_eff)
               iceav = iceav + lake_icefrac(i)*dz_lake(i)
               nav = nav + dz_lake(i)
            ENDDO

            qav = qav/nav
            iceav = iceav/nav
            !IF the average temperature is above freezing, put the extra energy into the water.
            !IF it is below freezing, take it away from the ice.
            IF (qav > 0.) THEN
               tav_froz = 0. !Celsius
               tav_unfr = qav / ((1. - iceav)*cwat)
            ELSE IF (qav < 0.) THEN
               tav_froz = qav / (iceav*cice_eff)
               tav_unfr = 0. !Celsius
            ELSE
               tav_froz = 0.
               tav_unfr = 0.
            ENDIF
         ENDIF

         IF (nav > 0.) THEN
            DO i = 1, j+1

               !put all the ice at the top.
               !IF the average temperature is above freezing, put the extra energy into the water.
               !IF it is below freezing, take it away from the ice.
               !for the layer with both ice & water, be careful to USE the average temperature
               !that preserves the correct total heat content given what the heat capacity of that
               !layer will actually be.

               IF (i == 1) zsum = 0.
               IF ((zsum+dz_lake(i))/nav <= iceav) THEN
                  lake_icefrac(i) = 1.
                  t_lake(i) = tav_froz + tfrz
               ELSE IF (zsum/nav < iceav) THEN
                  lake_icefrac(i) = (iceav*nav - zsum) / dz_lake(i)
                  ! Find average value that preserves correct heat content.
                  t_lake(i) = ( lake_icefrac(i)*tav_froz*cice_eff &
                              + (1. - lake_icefrac(i))*tav_unfr*cwat ) &
                              / ( lake_icefrac(i)*cice_eff + (1-lake_icefrac(i))*cwat ) + tfrz
               ELSE
                  lake_icefrac(i) = 0.
                  t_lake(i) = tav_unfr + tfrz
               ENDIF
               zsum = zsum + dz_lake(i)

               rhow(i) = (1.-lake_icefrac(i))*1000.*(1.-1.9549e-05*(abs(t_lake(i)-277.))**1.68) &
                           + lake_icefrac(i)*denice
            ENDDO
         ENDIF
      ENDDO

      lkrho = rhow
!------------------------------------------------------------
!*[6] Re-evaluate thermal properties and sum energy content.
!------------------------------------------------------------
      ! for lake
      DO j = 1, nl_lake
         cv_lake(j) = dz_lake(j) * (cwat*(1.-lake_icefrac(j)) + cice_eff*lake_icefrac(j))
      ENDDO

      ! DO as above to sum energy content
      ncvts = 0.
      DO j = 1, nl_lake
         ncvts = ncvts + cv_lake(j)*(t_lake(j)-tfrz) + cfus*dz_lake(j)*(1.-lake_icefrac(j))
      ENDDO

      DO j = lb, nl_soil
         ncvts = ncvts + cv_soisno(j)*(t_soisno(j)-tfrz) + hfus*wliq_soisno(j)
         IF (j == 1 .and. scv > 0. .and. j == lb) THEN
               ncvts = ncvts - scv*hfus
         ENDIF
      ENDDO

      ! check energy conservation.
      errsoi = (ncvts-ocvts)/deltim - fgrnd
      IF (abs(errsoi) < 0.10) THEN
         fseng = fseng - errsoi
         fsena = fseng
         fgrnd = fgrnd + errsoi
         errsoi = 0.
      ELSE
         print*, "energy conservation error in LAND WATER COLUMN during convective mixing", errsoi,fgrnd,ncvts,ocvts
      ENDIF


   END SUBROUTINE CoLMLTem



   SUBROUTINE CoLML_SurfFlux ( &
            ! "in" arguments
            ! -------------------
            patchtype   , snl          , zopt         , betaopt   ,&
            fetchopt    , rlat         , lakedepth    , deltim    ,& 
            savedtke1   , forc_hgt_u   , forc_hgt_t   , forc_hgt_q,&
            forc_us     , forc_vs      , forc_t       , forc_q    ,&
            forc_rhoair , forc_psrf    , forc_sols    , forc_soll ,&
            forc_solsd  , forc_solld   , forc_frl     , sabg      ,&
            hpbl        , dzlaktop     , zlakebot     , dzssbtop  ,&
            lktmptop    , tssbtop      , xwliqtop     , xwicetop  ,&
            icefrtop    , zcb          , ipatch       ,&
            ! "inout" arguments
            ! -------------------
            t_grnd      , z0m          , z0h          , z0q       ,&
            felak       , btpri        ,&
            ! "out" arguments
            ! -------------------
            fseng       , fevpg        , fsena        , fevpa     ,&
            lfevpa      , olrg         , fgrnd        , fgrnd1    ,&
            zol         , rib          , trad         , htvp      ,&
            emis        , wdm          , ram          , rah       ,&
            raw         , shfdt        , taux         , tauy      ,&
            tref        , qref         , u10m         , v10m      ,&
            fh2m        , fq2m         , fm10m        , fq10m     ,&
            fm          , fh           , fq           , ustar     ,&
            qstar       , tstar        , rhosnow      , u2m       ,&
            urban_call)

! ------------------------ code history ---------------------------
! Description:
!     Lake surface flux calculations
!
! Called:
!    -> roughness_lake    : Calculating lake surface roughness
!    -> qsadv:            : computes saturation mixing ratio and change in saturation mixing ratio with respect to temperature.
!    -> moninobuk:        : calculation of friction velocity, relation for potential temperatur and humidity profiles of surface boundary layer.
!    -> moninobukini:     : initialzation of Monin-Obukhov length. 
!
! Original author:
!     Yongjiu Dai, 2000
!
! Revisions:
!     Yongjiu Dai, 2000, 12/2012, 04/2014, 06/2018
!     Zack Subin, 2009
!     Nan Wei, 06/2018
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024: Separated this SUBROUTINE in order for it to be called by other lake models
!------------------------------------------------------------------
   USE MOD_Lake_Const, only : tfrz, hvap, hfus, hsub, grav, vonkar, tkice, tkair, stefnc, cpair,&
                              rgas, emisw, cur0, curm, fcrit, SHALLOW, DEEP
   USE MOD_Qsadv, only : qsadv
   USE MOD_FrictionVelocity, only: moninobukini, moninobuk
   USE MOD_TurbulenceLEddy, only: moninobuk_leddy
   USE MOD_Namelist, only: DEF_USE_CBL_HEIGHT, DEF_USE_SNICAR
!================================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in)      :: &
      ipatch                  ,&! patch index
      patchtype               ,&! land patch type (1=deep lake, 2=shallow lake)
      snl                     ,&! number of snow layers
      betaopt                 ,&! option for betaprime calculation, 1: constant, 2: equation
      fetchopt                ,&! option for fetch length, 1: constant, 2: equation
      zopt                      ! option for roughness length, 1: constant, 2: Subin et al. (2012)

   real(r8), intent(in)     :: &
      rlat                    ,&! latitude in radians
      lakedepth               ,&! column lake depth (m)
      deltim                  ,&! time step length
      forc_hgt_u              ,&! observational height of wind [m]
      forc_hgt_t              ,&! observational height of temperature [m]
      forc_hgt_q              ,&! observational height of humidity [m]
      forc_us                 ,&! wind component in eastward direction [m/s]
      forc_vs                 ,&! wind component in northward direction [m/s]
      forc_t                  ,&! temperature at agcm reference height [kelvin]
      forc_q                  ,&! specific humidity at agcm reference height [kg/kg]
      forc_rhoair             ,&! density air [kg/m3]
      forc_psrf               ,&! atmosphere pressure at the surface [pa]
      forc_sols               ,&! atm vis direct beam solar rad onto srf [W/m2]
      forc_soll               ,&! atm nir direct beam solar rad onto srf [W/m2]
      forc_solsd              ,&! atm vis diffuse solar rad onto srf [W/m2]
      forc_solld              ,&! atm nir diffuse solar rad onto srf [W/m2]
      forc_frl                ,&! atmospheric infrared (longwave) radiation [W/m2]
      hpbl                    ,&! atmospheric boundary layer height [m]
      sabg                    ,&! solar radiation absorbed by ground [W/m2]
      zcb                     ,&! convective boundary height [m]
      savedtke1               ,&! saved value of TKE
      dzlaktop                ,&! thickness of top lake layer [m]
      zlakebot                ,&! depth of top lake layer [m]
      lktmptop                ,&! temperature of top lake layer [K]
      icefrtop                ,&! lake mass fraction of top lake layer that is frozen
      dzssbtop                ,&! thickness of top lake layer [m]
      xwicetop                ,&! ice mass in top lake layer [kg/m2]
      xwliqtop                ,&! liquid water mass in top lake layer [kg/m2]
      tssbtop                   ! temperature of top lake layer [K]

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
      tref                    ,&! 2 m height air temperature [kelvin]
      qref                    ,&! 2 m height specific humidity [kg/kg]
      u10m                    ,&! wind speed at 10m [m/s]
      v10m                    ,&! wind speed at 10m [m/s]
      zol                     ,&! dimensionless height (z/L) used in Monin-Obukhov theory
      rib                     ,&! bulk Richardson number in surface layer
      fm                      ,&! integral of profile function for momentum
      fh                      ,&! integral of profile function for heat
      fq                      ,&! integral of profile function for moisture
      fh2m                    ,&! relation for temperature at 2m
      fq2m                    ,&! relation for specific humidity at 2m
      fm10m                   ,&! integral of profile function for momentum at 10m
      fq10m                   ,&! integral of profile function for moisture at 10m
      ustar                   ,&! u* in similarity theory [m/s]
      qstar                   ,&! humidity scaling parameter [kg/kg]
      tstar                   ,&! temperature scaling parameter [kelvin]
      fseng                   ,&! sensible heat flux into the ground [W/m2]
      fevpg                   ,&! latent heat flux into the ground [W/m2]
      fsena                   ,&! sensible heat from canopy height to atmosphere [W/m2]
      fevpa                   ,&! evapotranspiration from canopy height to atmosphere [mm/s]
      htvp                    ,&! heat flux into the lake [W/m2]
      lfevpa                  ,&! latent heat flux from canopy height to atmosphere [W/m2]
      olrg                    ,&! outgoing longwave radiation [W/m2]
      fgrnd                   ,&! ground heat flux [W/m2]
      fgrnd1                  ,&! ground heat flux [W/m2]
      rhosnow                 ,&! partitial density of water (ice + liquid)
      trad                      ! radiative temperature [K]

   logical, optional, intent(in) :: urban_call   ! whether it is a urban CALL

!  ------------------------- local variables ---------------------------
   real(r8)                 :: &
      fetch                   ,&! lake fetch (m)
      z0mg                    ,&! roughness length over ground, momentum [m]
      z0hg                    ,&! roughness length over ground, sensible heat [m]
      z0qg                    ,&! roughness length over ground, humidity [m]
      betaprime               ,&! fraction of solar radiation in the NIR
      betavis                 ,&! middle variable for betaprime calculation
      eg                      ,&! water vapor pressure at temperature T [pa]
      qsatg                   ,&! saturated humidity [kg/kg]
      qsatgdT                 ,&! d(qsatg)/dT
      degdT                   ,&! d(eg)/dT
      beta1                   ,&! coefficient of conective velocity [-]
      zii                     ,&! convective boundary height [m]
      u2m                     ,&! wind speed at 2m [m/s]
      um                      ,&! wind speed including the stablity effect [m/s]
      ws                      ,&! surface friction velocity (m/s)
      ks                      ,&! coefficient passed to SLakeTemperature for calculation of decay of eddy diffusivity with depth
      th                      ,&! potential temperature (kelvin)
      thv                     ,&! virtual potential temperature (kelvin)
      ur                      ,&! wind speed at reference height [m/s]
      wc                      ,&! convective velocity [m/s]
      obuold                  ,&! monin-obukhov length of previous iteration
      dthv                    ,&! diff of vir. poten. temp. between ref. height and surface
      zldis                   ,&! reference height "minus" zero displacement heght [m]
      visa                    ,&! kinematic viscosity of dry air [m2/s]
      cur                     ,&! Charnock parameter (-)
      obu                     ,&! monin-obukhov length (m)
      del_T_grnd              ,&! t_grnd diff
      t_grnd_bef              ,&! initial ground temperature
      tksur                   ,&! thermal conductivity of snow/soil (w/m/kelvin)
      tsur                    ,&! top layer temperature
      dzsur                   ,&! top layer thickness
      displax                 ,&! zero- displacement height [m]
      emg                     ,&! emissivity of the surface
      stftg3                  ,&! emg*sb*t_grnd*t_grnd*t_grnd
      ax                      ,&! used in iteration loop for calculating t_grnd (numerator of NR solution)
      bx                      ,&! used in iteration loop for calculating t_grnd (denomin. of NR solution)
      thvstar                 ,&! virtual potential temperature scaling parameter
      zeta                    ,&! dimensionless height used in Monin-Obukhov theory
      wc2                     ,&! wc*wc
      raih                    ,&! temporary variable [kg/m2/s]
      raiw                    ,&! temporary variable [kg/m2/s]
      shfdtl                  ,&! derivative of srf sensible heat flux wrt srf temp [W/m2/K]
      shfdts                  ,&! derivative of srf latent   heat flux wrt srf temp [W/m2/K]
      thm                     ,&! intermediate variable (forc_t+0.0098*forc_hgt_t)
      dth                     ,&! diff of virtual temp. between ref. height and surface
      dqh                     ,&! diff of humidity between ref. height and surface
      tdmax                     ! temperature of maximum water density

   integer                  :: &
      lb                      ,&! lower bound of arrays
      i,j                     ,&! DO loop or array index
      iter                    ,&! iteration index
      convernum               ,&! number of time when del_T_grnd < 0.01
      nmozsgn                   ! number of times moz changes sign

!  ------------------------- constants variables ---------------------------
   integer,  parameter :: itmax  = 40         ! maximum number of iteration
   integer,  parameter :: itmin  = 6          ! minimum number of iteration
   real(r8), parameter :: delmax = 3.0        ! maximum change in lake temperature [K]
   real(r8), parameter :: dtmin  = 0.01       ! max limit for temperature convergence [K]
   real(r8), parameter :: dlemin = 0.1        ! max limit for energy flux convergence [w/m2]
!================================================================================
      emg = emisw    ! surface emissivity for water and ice [0.97]

      ! latent heat
      IF (t_grnd > tfrz )THEN
         htvp = hvap
      ELSE
         htvp = hsub
      ENDIF

      ! Base on lake depth, assuming that small lakes are likely to be shallower
      ! Estimate crudely based on lake depth
      !-WMEJ TODO: USE scwat to control lake classification, 
      !            1 = shallow lake, 2 = deep lake, 2024-03-23
      IF (fetchopt == 1) THEN
         fetch = felak
      ELSE
         IF (lakedepth < 4.) THEN
               fetch = 100. ! shallow lake
         ELSE
               fetch = 25.*lakedepth ! deep lake
         ENDIF
      ENDIF 

! ======================================================================
! pre-processing for the calcilation of the surface fluxes
! ======================================================================
      !+WMEJ add a new option for the calculation of betaprime, 1: constant, 2: equation
      !      There may be other calculation methods in the future,
      !      Can also be used to make adjustments for specific lakes.
      IF (betaopt == 1) THEN
         betaprime = btpri
      ELSE
         IF (.not. DEF_USE_SNICAR .or. present(urban_call)) THEN
            IF (snl == 0) THEN
               ! calculate the nir fraction of absorbed solar.
               betaprime = (forc_soll+forc_solld)/max(1.e-5,forc_sols+forc_soll+forc_solsd+forc_solld)
               betavis = 0. ! The fraction of the visible (e.g. vis not nir from atm) sunlight
                           ! absorbed in ~1 m of water (the surface layer za_lake).
                           ! This is roughly the fraction over 700 nm but may depend on the details
                           ! of atmospheric radiative transfer.
                           ! As long as NIR = 700 nm and up, this can be zero.
               betaprime = betaprime + (1.0-betaprime)*betavis
            ELSE
               ! or frozen but no snow layers or
               ! currently ignor the transmission of solar in snow and ice layers
               ! to be updated in the future version
               betaprime = 1.0
            ENDIF
         ELSE
            ! calculate the nir fraction of absorbed solar.
            betaprime = (forc_soll+forc_solld)/max(1.e-5,forc_sols+forc_soll+forc_solsd+forc_solld)
            betavis = 0. ! The fraction of the visible (e.g. vis not nir from atm) sunlight
                        ! absorbed in ~1 m of water (the surface layer za_lake).
                        ! This is roughly the fraction over 700 nm but may depend on the details
                        ! of atmospheric radiative transfer.
                        ! As long as NIR = 700 nm and up, this can be zero.
            betaprime = betaprime + (1.0-betaprime)*betavis
         ENDIF
      ENDIF

      CALL qsadv(t_grnd,forc_psrf,eg,degdT,qsatg,qsatgdT)
      ! potential temperatur at the reference height
      beta1 = 1.       ! -  (in computing W_*)
      zii   = zcb      ! m  (pbl height)  ! -WMEJ 1000.
      thm   = forc_t + 0.0098*forc_hgt_t  ! intermediate variable equivalent to
                                          ! forc_t*(pgcm/forc_psrf)**(rgas/cpair)
      th    = forc_t*(100000./forc_psrf)**(rgas/cpair) ! potential T
      thv   = th*(1.+0.61*forc_q)         ! virtual potential T
      ur    = max(0.1,sqrt(forc_us*forc_us+forc_vs*forc_vs))   ! limit set to 0.1

      ! Initialization variables
      nmozsgn = 0
      obuold  = 0.
      dth     = thm-t_grnd
      dqh     = forc_q-qsatg
      dthv    = dth*(1.+0.61*forc_q)+0.61*th*dqh
      zldis   = forc_hgt_u-0.

      ! Roughness lengths, allow all roughness lengths to be prognostic
      ustar = 0.06
      wc    = 0.5

      ! Kinematic viscosity of dry air (m2/s)- Andreas (1989) CRREL Rep. 89-11
      visa = 1.326e-5*(1.+6.542e-3*(forc_t-tfrz) &
           + 8.301e-6*(forc_t-tfrz)**2 - 4.84e-9*(forc_t-tfrz)**3)

      cur = cur0 + curm * exp( max( -(fetch*grav/ur/ur)**(1./3.)/fcrit, & ! Fetch-limited
                                       -(zlakebot*grav)**0.5/ur ) )   ! depth-limited

      IF(dthv.ge.0.) THEN
         um = max(ur,0.1)
      ELSE
         um = sqrt(ur*ur+wc*wc)
      ENDIF

      IF (zopt == 1) THEN
         z0mg = z0m
         z0hg = z0h
         z0qg = z0q
         ustar = vonkar*um/log(zldis/z0mg)
      ELSE
         DO i = 1,5
            z0mg = 0.013*ustar*ustar/grav+0.11*visa/ustar
            ustar = vonkar*um/log(zldis/z0mg)
         ENDDO
         !+WMEJ USE tmtop to instead t_lake(1), USE icefrac_top to instead lake_icefrac(1)
         CALL roughness_lake (snl,t_grnd,lktmptop,forc_psrf,&
                              cur,ustar,z0mg,z0hg,z0qg)
      ENDIF 

      CALL moninobukini(ur,th,thm,thv,dth,dqh,dthv,zldis,z0mg,um,obu)

      !+WMEJ USE dzlktop to instead dz_lake(1), USE dzssbtop to instead z_soisno(lb)-zi_soisno(lb-1)
      IF (snl == 0) THEN
         dzsur = dzlaktop/2.
      ELSE
         dzsur = dzssbtop/2. !z_soisno(lb)-zi_soisno(lb-1)
      ENDIF

      iter = 1
      del_T_grnd = 1.0    ! t_grnd diff
      convernum = 0       ! number of time when del_T_grnd <= 0.01

! ======================================================================
! Begin stability iteration and temperature and fluxes calculation
! ======================================================================
      ! =====================================
      ITERATION : DO WHILE (iter <= itmax)
      ! =====================================

         t_grnd_bef = t_grnd

         IF (t_grnd_bef > tfrz .and. lktmptop > tfrz .and. snl == 0) THEN
            tksur = savedtke1       !water molecular conductivity
            tsur = lktmptop
            htvp = hvap
         ELSE IF (snl == 0) THEN !frozen but no snow layers
            tksur = tkice       ! This is an approximation because the whole layer may not be frozen, and it is not
                              ! accounting for the physical (but not nominal) expansion of the frozen layer.
            tsur = lktmptop
            htvp = hsub
         ELSE
         ! need to calculate thermal conductivity of the top snow layer
            !+WMEJ USE xwicetop to instead wice_soisno(lb), USE xwliqtop to instead wliq_soisno(lb)
            rhosnow = (xwicetop+xwliqtop)/dzssbtop
            tksur = tkair + (7.75e-5*rhosnow + 1.105e-6*rhosnow*rhosnow)*(tkice-tkair)
            tsur = tssbtop
            htvp = hsub
         ENDIF

         ! Evaluated stability-dependent variables using moz from prior iteration
         displax = 0.
         ! fq10m needs to be calculated in moninobuk (IF needed)
         IF (DEF_USE_CBL_HEIGHT) THEN    
            CALL moninobuk_leddy(forc_hgt_u,forc_hgt_t,forc_hgt_q,displax,z0mg,z0hg,z0qg,obu,um, hpbl, &
                        ustar,fh2m,fq2m,fm10m,fm,fh,fq)
         ELSE
            CALL moninobuk(forc_hgt_u,forc_hgt_t,forc_hgt_q,displax,z0mg,z0hg,z0qg,obu,um,&
                        ustar,fh2m,fq2m,fm10m,fm,fh,fq)
         ENDIF

         ! Get derivative of fluxes with repect to ground temperature
         ram    = 1./(ustar*ustar/um)
         rah    = 1./(vonkar/fh*ustar)
         raw    = 1./(vonkar/fq*ustar)
         stftg3 = emg*stefnc*t_grnd_bef*t_grnd_bef*t_grnd_bef

         ! +WMEJ TODO: Soil layers are not considered in Simstrat, 
         !             therefore the geothermal heat flux(G, tksur*tsur/dzsur) should be replaced by the constant in Simstrat.
         !             Geothermal fluxes have minimal impact on lake temperatures, making this practice acceptable
         ax  = betaprime*sabg + emg*forc_frl + 3.*stftg3*t_grnd_bef &
               + forc_rhoair*cpair/rah*thm &
               - htvp*forc_rhoair/raw*(qsatg-qsatgdT*t_grnd_bef - forc_q) &
               + tksur*tsur/dzsur

         bx  = 4.*stftg3 + forc_rhoair*cpair/rah &
               + htvp*forc_rhoair/raw*qsatgdT + tksur/dzsur

         t_grnd = ax/bx
         
         !-----------------------------------------------------------------
         ! h_fin = betaprime*sabg + emg*forc_frl + 3.*stftg3*t_grnd_bef & !
         !     + forc_rhoair*cpair/rah*thm &                              !
         !     - htvp*forc_rhoair/raw*(qsatg-qsatgdT*t_grnd_bef - forc_q) !
         ! h_finDT = 4.*stftg3 + forc_rhoair*cpair/rah &                  !
         !     + htvp*forc_rhoair/raw*qsatgdT                             !
         ! del_T_grnd = t_grnd - t_grnd_bef                               !
         !----------------------------------------------------------------!

         ! surface fluxes of momentum, sensible and latent
         ! using ground temperatures from previous time step

         fseng = forc_rhoair*cpair*(t_grnd-thm)/rah
         fevpg = forc_rhoair*(qsatg+qsatgdT*(t_grnd-t_grnd_bef)-forc_q)/raw

         CALL qsadv(t_grnd,forc_psrf,eg,degdT,qsatg,qsatgdT)
         dth = thm-t_grnd
         dqh = forc_q-qsatg
         tstar = vonkar/fh*dth
         qstar = vonkar/fq*dqh
         thvstar = tstar*(1.+0.61*forc_q)+0.61*th*qstar
         zeta = zldis*vonkar*grav*thvstar/(ustar**2*thv)

         IF(zeta >= 0.) THEN     !stable
            zeta = min(2.,max(zeta,1.e-6))
         ELSE                    !unstable
            zeta = max(-100.,min(zeta,-1.e-6))
         ENDIF

         obu = zldis/zeta

         IF(zeta >= 0.)THEN
            um = max(ur,0.1)
         ELSE
            IF (DEF_USE_CBL_HEIGHT) THEN !//TODO: Shaofeng, 2023.05.18
               zii = max(5.*forc_hgt_u,hpbl)
            ENDIF !//TODO: Shaofeng, 2023.05.18
            wc = (-grav*ustar*thvstar*zii/thv)**(1./3.)
            wc2 = beta1*beta1*(wc*wc)
            um = sqrt(ur*ur+wc2)
         ENDIF

         IF (zopt == 1) THEN
            z0mg = z0m
            z0hg = z0h
            z0qg = z0q
         ELSE
            CALL roughness_lake (snl,t_grnd,lktmptop,forc_psrf,&
                                 cur,ustar,z0mg,z0hg,z0qg)
         ENDIF 

         iter = iter + 1
         del_T_grnd = abs(t_grnd - t_grnd_bef)

         IF(iter .gt. itmin) THEN
            IF(del_T_grnd <= dtmin) THEN
            convernum = convernum + 1
            ENDIF
            IF(convernum >= 4) EXIT
         ENDIF

      ! ===============================================
      ENDDO ITERATION   ! END of stability iteration
      ! ===============================================
      !*----------------------------------------------------------------------
      !*Zack Subin, 3/27/09
      !*Since they are now a function of whatever t_grnd was before cooling
      !*to freezing temperature, THEN this value should be used in the derivative correction term.
      !*Allow convection IF ground temp is colder than lake but warmer than 4C, or warmer than
      !*lake which is warmer than freezing but less than 4C.
      tdmax = tfrz + 4.0
      IF ( (snl < 0 .or. lktmptop <= tfrz) .and. t_grnd > tfrz) THEN
         t_grnd_bef = t_grnd
         t_grnd = tfrz
         fseng = forc_rhoair*cpair*(t_grnd-thm)/rah
         fevpg = forc_rhoair*(qsatg+qsatgdT*(t_grnd-t_grnd_bef)-forc_q)/raw
      
      ELSE IF ( (lktmptop > t_grnd .and. t_grnd > tdmax) .or. &
               (lktmptop < t_grnd .and. lktmptop > tfrz .and. t_grnd < tdmax) ) THEN
               ! Convective mixing will occur at surface
         t_grnd_bef = t_grnd
         t_grnd = lktmptop
         fseng = forc_rhoair*cpair*(t_grnd-thm)/rah
         fevpg = forc_rhoair*(qsatg+qsatgdT*(t_grnd-t_grnd_bef)-forc_q)/raw
      ENDIF
      !*----------------------------------------------------------------------

      ! net longwave from ground to atmosphere
      stftg3 = emg*stefnc*t_grnd_bef*t_grnd_bef*t_grnd_bef
      olrg = (1.-emg)*forc_frl + emg*stefnc*t_grnd_bef**4 + 4.*stftg3*(t_grnd - t_grnd_bef)
      IF (t_grnd > tfrz )THEN
         htvp = hvap
      ELSE
         htvp = hsub
      ENDIF

      u2m = max(0.1_r8,ustar/vonkar*log(2._r8/z0mg))

      !The actual heat flux from the ground interface into the lake, not including the light that penetrates the surface.
      fgrnd1 = betaprime*sabg + forc_frl - olrg - fseng - htvp*fevpg

      ! radiative temperature
      trad = (olrg/stefnc)**0.25

      ! solar absorption below the surface.
      fgrnd = sabg + forc_frl - olrg - fseng - htvp*fevpg
      fsena = fseng
      fevpa = fevpg
      lfevpa = htvp*fevpg

      !+WMEJ Recalculation of supplementary variables
      ws = 1.2e-03_r8 * u2m
      ks = 6.6_r8*sqrt(abs(sin(rlat)))*(u2m**(-1.84_r8))
      z0m    = z0mg
      z0h    = z0hg
      z0q    = z0qg
      felak  = fetch
      btpri  = betaprime
      emis   = emg
      wdm    = um
      ram    = 1./(ustar*ustar/um)
      rah    = 1./(vonkar/fh*ustar)
      raw    = 1./(vonkar/fq*ustar)
      raih   = forc_rhoair*cpair/rah
      raiw   = forc_rhoair/raw
      shfdts = raih
      shfdtl = raiw*qsatgdT
      shfdt  = shfdts + htvp*shfdtl
      taux   = -forc_rhoair*forc_us/ram
      tauy   = -forc_rhoair*forc_vs/ram
      tref   = thm + vonkar/fh*dth * (fh2m/vonkar - fh/vonkar)
      qref   = forc_q + vonkar/fq*dqh * (fq2m/vonkar - fq/vonkar)
      u10m   = forc_us/ur * ustar/vonkar * fm10m
      v10m   = forc_vs/ur * ustar/vonkar * fm10m
      u2m    = max(0.1_r8,ustar/vonkar*log(2._r8/z0mg))
      zol    = zeta
      rib    = min(5.,zol*ustar**2/(vonkar*vonkar/fh*um**2))

   END SUBROUTINE CoLML_SurfFlux



   SUBROUTINE snowwater_lake ( &
            ! "in" arguments
            ! ---------------------------
            maxsnl      , nl_soil     , nl_lake   , deltim       ,&
            qseva       , qsubl       , qsdew     , qfros        ,&
            porsl       , pg_rain     , pg_snow   , dz_lake      ,&
            imelt       , fiold      ,&
            ! "inout" arguments
            ! ---------------------------
            z_soisno    , dz_soisno   , zi_soisno , t_soisno     ,&
            wice_soisno , wliq_soisno , t_lake    , lake_icefrac ,&
            fseng       , fgrnd       , snl       , scv          ,&
            qout_snowb  ,&
            ! SNICAR model variables
            forc_aer    ,&
            mss_bcpho   , mss_bcphi   , mss_ocpho , mss_ocphi    ,&
            mss_dst1    , mss_dst2    , mss_dst3  , mss_dst4     ,&
            ! END SNICAR model variables
            snowdp      , sm          , forc_us   , forc_vs      ,&
            urban_call)

!-----------------------------------------------------------------------------------------------
! Calculation of Lake Hydrology. Lake water mass is kept constant. The soil is simply maintained at
! volumetric saturation IF ice melting frees up pore space.
!
! Called:
!    -> snowwater:                  change of snow mass and snow water onto soil
!    -> snowcompaction:             compaction of snow layers
!    -> snowlayerscombine:          combine snow layers that are thinner than minimum     
!    -> snowlayersdivide:           subdivide snow layers that are thicker than maximum
!
! Initial: Yongjiu Dai, December, 2012
!                          April, 2014
! REVISIONS:
!     Nan Wei, 06/2018: update snow hydrology above lake
!     Yongjiu Dai, 01/2023: added for SNICAR model effects for snowwater,
!                           combinesnowlayers, dividesnowlayers processes by calling snowwater_snicar(),
!                           SnowLayersCombine_snicar, SnowLayersDivide_snicar()
!-----------------------------------------------------------------------------------------------
   USE MOD_Lake_Const, only : denh2o, denice, hfus, tfrz, cpliq, cpice, ssi, wimp 
   USE MOD_SoilThermalParameters, only: soil_hcap_cond
   USE MOD_SoilSnowHydrology, only : snowwater, SnowWater_snicar
   USE MOD_SnowLayersCombineDivide, only: SnowLayersCombine_snicar, SnowLayersDivide_snicar, snowcompaction, snowlayerscombine, snowlayersdivide
   USE MOD_Namelist, only: DEF_USE_SNICAR
!================================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in)     :: maxsnl                        ! maximum number of snow layers
   integer, intent(in)     :: nl_soil                       ! number of soil layers
   integer, intent(in)     :: nl_lake                       ! number of soil layers
   real(r8), intent(in)    :: deltim                        ! seconds in a time step (sec)
   real(r8), intent(in)    :: forc_us                       ! wind speed at reference height [m/s]
   real(r8), intent(in)    :: forc_vs                       ! wind speed at reference height [m/s]
   real(r8), intent(in)    :: pg_rain                       ! rainfall incident on ground [mm/s]
   real(r8), intent(in)    :: pg_snow                       ! snowfall incident on ground [mm/s]
   real(r8), intent(in)    :: porsl(1:nl_soil)              ! volumetric soil water at saturation (porosity)
   real(r8), intent(in)    :: dz_lake(1:nl_lake)            ! layer thickness for lake (m)
   integer,  intent(in)    :: imelt(maxsnl+1:0)             ! signifies IF node in melting (imelt = 1)
   real(r8), intent(in)    :: fiold(maxsnl+1:0)             ! fraction of ice relative to the total water content at the previous time step
   real(r8), intent(in)    :: qseva                         ! ground surface evaporation rate (mm h2o/s)
   real(r8), intent(in)    :: qsubl                         ! sublimation rate from snow pack (mm H2O /s) [+]
   real(r8), intent(in)    :: qsdew                         ! surface dew added to snow pack (mm H2O /s) [+]
   real(r8), intent(in)    :: qfros                         ! ground surface frosting formation (mm H2O /s) [+]

!  ------------------------- inout arguments ---------------------------
   real(r8), intent(inout) :: z_soisno(maxsnl+1:nl_soil)    ! layer depth  (m)
   real(r8), intent(inout) :: dz_soisno(maxsnl+1:nl_soil)   ! layer thickness depth (m)
   real(r8), intent(inout) :: zi_soisno(maxsnl:nl_soil)     ! interface depth (m)
   real(r8), intent(inout) :: t_soisno(maxsnl+1:nl_soil)    ! snow temperature (Kelvin)
   real(r8), intent(inout) :: wice_soisno(maxsnl+1:nl_soil) ! ice lens (kg/m2)
   real(r8), intent(inout) :: wliq_soisno(maxsnl+1:nl_soil) ! liquid water (kg/m2)
   real(r8), intent(inout) :: t_lake(1:nl_lake)             ! lake temperature (Kelvin)
   real(r8), intent(inout) :: lake_icefrac(1:nl_lake)       ! mass fraction of lake layer that is frozen
   real(r8), intent(inout) :: qout_snowb                    ! rate of water out of snow bottom (mm/s)
   real(r8), intent(inout) :: fseng                         ! total sensible heat flux (W/m**2) [+ to atm]
   real(r8), intent(inout) :: fgrnd                         ! heat flux into snow / lake (W/m**2) [+ = into soil]
   integer , intent(inout) :: snl                           ! number of snow layers
   real(r8), intent(inout) :: scv                           ! snow water (mm H2O)
   real(r8), intent(inout) :: snowdp                        ! snow height (m)
   real(r8), intent(inout) :: sm                            ! rate of snow melt (mm H2O /s)

   logical, optional, intent(in) :: urban_call              ! whether it is a urban CALL

! SNICAR model variables
   ! Aerosol Fluxes (Jan. 07, 2023 by Yongjiu Dai)
   real(r8), intent(in)    :: forc_aer (14)                 ! aerosol deposition from atmosphere model (grd,aer) [kg m-1 s-1]
   real(r8), intent(inout) :: mss_bcpho(maxsnl+1:0)         ! mass of hydrophobic BC in snow  (col,lyr) [kg]
   real(r8), intent(inout) :: mss_bcphi(maxsnl+1:0)         ! mass of hydrophillic BC in snow (col,lyr) [kg]
   real(r8), intent(inout) :: mss_ocpho(maxsnl+1:0)         ! mass of hydrophobic OC in snow  (col,lyr) [kg]
   real(r8), intent(inout) :: mss_ocphi(maxsnl+1:0)         ! mass of hydrophillic OC in snow (col,lyr) [kg]
   real(r8), intent(inout) :: mss_dst1 (maxsnl+1:0)         ! mass of dust species 1 in snow  (col,lyr) [kg]
   real(r8), intent(inout) :: mss_dst2 (maxsnl+1:0)         ! mass of dust species 2 in snow  (col,lyr) [kg]
   real(r8), intent(inout) :: mss_dst3 (maxsnl+1:0)         ! mass of dust species 3 in snow  (col,lyr) [kg]
   real(r8), intent(inout) :: mss_dst4 (maxsnl+1:0)         ! mass of dust species 4 in snow  (col,lyr) [kg]
   ! Aerosol Fluxes (Jan. 07, 2023)
! END SNICAR model variables

!  ------------------------- local variables ---------------------------
   integer  :: j                                            ! indices
   integer  :: lb                                           ! lower bound of array
   real(r8) :: xmf                                          ! snow melt heat flux (W/m**2)
   real(r8) :: sumsnowice                                   ! sum of snow ice IF snow layers found above unfrozen lake [kg/m&2]
   real(r8) :: sumsnowliq                                   ! sum of snow liquid IF snow layers found above unfrozen lake [kg/m&2]
   logical  :: unfrozen                                     ! true IF top lake layer is unfrozen with snow layers above
   real(r8) :: heatsum                                      ! used in case above [J/m^2]
   real(r8) :: heatrem                                      ! used in case above [J/m^2]
   real(r8) :: a, b, c, d   
   real(r8) :: wice_lake(1:nl_lake)                         ! ice lens (kg/m2)
   real(r8) :: wliq_lake(1:nl_lake)                         ! liquid water (kg/m2)
   real(r8) :: t_ave, frac_
!================================================================================

      ! for runoff calculation (assumed no mass change in the land water bodies)
      lb = snl + 1
      qout_snowb = 0.0

! ----------------------------------------------------------
!*[1] snow layer on frozen lake
! ----------------------------------------------------------
      IF (snl < 0) THEN
         lb = snl + 1

         IF (DEF_USE_SNICAR .and. .not. present(urban_call)) THEN
            CALL snowwater_SNICAR (lb,deltim,ssi,wimp,&
                  pg_rain,qseva,qsdew,qsubl,qfros,&
                  dz_soisno(lb:0),wice_soisno(lb:0),wliq_soisno(lb:0),qout_snowb,    &
                  forc_aer,&
                  mss_bcpho(lb:0), mss_bcphi(lb:0), mss_ocpho(lb:0), mss_ocphi(lb:0),&
                  mss_dst1(lb:0),  mss_dst2(lb:0),  mss_dst3(lb:0),  mss_dst4(lb:0) )
         ELSE
            CALL snowwater (lb,deltim,ssi,wimp,&
                  pg_rain,qseva,qsdew,qsubl,qfros,&
                  dz_soisno(lb:0),wice_soisno(lb:0),wliq_soisno(lb:0),qout_snowb)
         ENDIF

         ! Natural compaction and metamorphosis.
         lb = snl + 1
         CALL snowcompaction (lb,deltim, &
               imelt(lb:0),fiold(lb:0),t_soisno(lb:0),&
               wliq_soisno(lb:0),wice_soisno(lb:0),forc_us,forc_vs,dz_soisno(lb:0))

         ! Combine thin snow elements
         lb = maxsnl + 1
         IF (DEF_USE_SNICAR .and. .not. present(urban_call)) THEN
            CALL snowlayerscombine_SNICAR (lb, snl,&
                  z_soisno(lb:1),dz_soisno(lb:1),zi_soisno(lb-1:0),&
                  wliq_soisno(lb:1),wice_soisno(lb:1), t_soisno(lb:1),scv,snowdp, &
                  mss_bcpho(lb:0), mss_bcphi(lb:0), mss_ocpho(lb:0), mss_ocphi(lb:0),&
                  mss_dst1(lb:0),  mss_dst2(lb:0),  mss_dst3(lb:0),  mss_dst4(lb:0))
         ELSE
            CALL snowlayerscombine (lb, snl,&
                  z_soisno(lb:1),dz_soisno(lb:1),zi_soisno(lb-1:0),&
                  wliq_soisno(lb:1),wice_soisno(lb:1),&
                  t_soisno(lb:1),scv,snowdp)
         ENDIF

         ! Divide thick snow elements
         IF (snl < 0) THEN
            IF (DEF_USE_SNICAR .and. .not. present(urban_call)) THEN
               CALL snowlayersdivide_SNICAR (lb,snl,z_soisno(lb:0),dz_soisno(lb:0),zi_soisno(lb-1:0),&
                     wliq_soisno(lb:0),wice_soisno(lb:0),t_soisno(lb:0)     ,&
                     mss_bcpho(lb:0), mss_bcphi(lb:0), mss_ocpho(lb:0), mss_ocphi(lb:0),&
                     mss_dst1(lb:0),  mss_dst2(lb:0),  mss_dst3(lb:0),  mss_dst4(lb:0) )
            ELSE
               CALL snowlayersdivide (lb,snl,z_soisno(lb:0),dz_soisno(lb:0),zi_soisno(lb-1:0),&
                     wliq_soisno(lb:0),wice_soisno(lb:0),t_soisno(lb:0))
            ENDIF
         ENDIF

! ----------------------------------------------------------
!*[2] check for single completely unfrozen snow layer over lake.
!     Modeling this ponding is unnecessary and can cause instability after the timestep
!     when melt is completed, as the temperature after melt can be excessive
!     because the fluxes were calculated with a fixed ground temperature of freezing, but the
!     phase change was unable to restore the temperature to freezing.  (Zack Subnin 05/2010)
! ----------------------------------------------------------
         IF (snl == -1 .and. wice_soisno(0) == 0.) THEN
            ! Remove layer
            ! Take extra heat of layer and release to sensible heat in order to maintain energy conservation.
            heatrem = cpliq*wliq_soisno(0)*(t_soisno(0) - tfrz)
            fseng = fseng + heatrem/deltim
            fgrnd = fgrnd - heatrem/deltim

            snl = 0
            scv = 0.
            snowdp = 0.
         ENDIF

      ENDIF

! ----------------------------------------------------------
!*[3] check for snow layers above lake with unfrozen top layer. Mechanically,
!     the snow will fall into the lake and melt or turn to ice. IF the top layer has
!     sufficient heat to melt the snow without freezing, THEN that will be done.
!     Otherwise, the top layer will undergo freezing, but only IF the top layer will
!     not freeze completely. Otherwise, let the snow layers persist and melt by diffusion.
! ----------------------------------------------------------
      IF (t_lake(1) > tfrz .and. snl < 0 .and. lake_icefrac(1) < 0.001) THEN ! for unfrozen lake
         unfrozen = .true.
      ELSE
         unfrozen = .false.
      ENDIF

      sumsnowice = 0.
      sumsnowliq = 0.
      heatsum = 0.0
      DO j = snl+1,0
         IF (unfrozen) THEN
            sumsnowice = sumsnowice + wice_soisno(j)
            sumsnowliq = sumsnowliq + wliq_soisno(j)
            heatsum = heatsum + wice_soisno(j)*cpice*(tfrz-t_soisno(j)) &
                              + wliq_soisno(j)*cpliq*(tfrz-t_soisno(j))
         ENDIF
      ENDDO

      IF (unfrozen) THEN
         ! changed by weinan as the SUBROUTINE newsnow_lake
         ! Remove snow and subtract the latent heat from the top layer.

         t_ave = tfrz - heatsum/(sumsnowice*cpice + sumsnowliq*cpliq)

         a = heatsum
         b = sumsnowice*hfus
         c = (t_lake(1) - tfrz)*cpliq*denh2o*dz_lake(1)
         d = denh2o*dz_lake(1)*hfus

         ! all snow melt
         IF (c>=a+b)THEN
            t_lake(1) = (cpliq*(denh2o*dz_lake(1)*t_lake(1) + (sumsnowice+sumsnowliq)*tfrz) - a - b) / &
                        (cpliq*(denh2o*dz_lake(1) + sumsnowice+ sumsnowice))
            sm = sm + scv/deltim
            scv = 0.
            snowdp = 0.
            snl = 0
         ! lake partially freezing to melt all snow
         ELSE IF(c+d >= a+b)THEN
            t_lake(1) = tfrz
            sm = sm + scv/deltim
            scv = 0.
            snowdp = 0.
            snl = 0
            lake_icefrac(1) = (a+b-c)/d

      !  snow do not melt while all lake freezing
      !    ELSE IF(c+d < a) THEN
      !     t_lake(1) = (c+d + cpice*(sumsnowice*t_ave+denh2o*dz_lake(1)*tfrz) + cpliq*sumsnowliq*t_ave)/&
      !                (cpice*(sumsnowice+denh2o*dz_lake(1))+cpliq*sumsnowliq)
      !     lake_icefrac(1) = 1.0
         ENDIF
      ENDIF

! ----------------------------------------------------------
!*[4] Soil water and ending water balance
! ----------------------------------------------------------
      ! Here this consists only of making sure that soil is saturated even as it melts and
      ! pore space opens up. Conversely, IF excess ice is melting and the liquid water exceeds the
      ! saturation value, THEN remove water.
      DO j = 1, nl_soil
         a = wliq_soisno(j)/(dz_soisno(j)*denh2o) + wice_soisno(j)/(dz_soisno(j)*denice)

         IF (a < porsl(j)) THEN
            wliq_soisno(j) = max( 0., (porsl(j)*dz_soisno(j) - wice_soisno(j)/denice)*denh2o )
            wice_soisno(j) = max( 0., (porsl(j)*dz_soisno(j) - wliq_soisno(j)/denh2o)*denice )
         ELSE
            wliq_soisno(j) = max(0., wliq_soisno(j) - (a - porsl(j))*denh2o*dz_soisno(j) )
            wice_soisno(j) = max( 0., (porsl(j)*dz_soisno(j) - wliq_soisno(j)/denh2o)*denice )
         ENDIF

         IF (wliq_soisno(j) > porsl(j)*denh2o*dz_soisno(j)) THEN
               wliq_soisno(j) = porsl(j)*denh2o*dz_soisno(j)
               wice_soisno(j) = 0.0
         ENDIF
      ENDDO

   END SUBROUTINE snowwater_lake



   SUBROUTINE roughness_lake (snl,t_grnd,lktmptop,forc_psrf,&
                              cur,ustar,z0mg,z0hg,z0qg)

!-----------------------------------------------------------------------
! DESCRIPTION:
! Calculate lake surface roughness
!
! Original:
! The Community Land Model version 4.5 (CLM4.5)
!
! Revisions:
! Yongjiu Dai, Nan Wei, 01/2018
!-----------------------------------------------------------------------
   USE MOD_Lake_Const, only : tfrz, vonkar, grav, cus, kva0, prn, sch
!================================================================================
!  -------------------------- "in" arguments ---------------------------
   integer,  intent(in) :: snl            ! number of snow layers
   real(r8), intent(in) :: t_grnd         ! ground temperature
   real(r8), intent(in) :: lktmptop       ! surface lake layer temperature [K]
   real(r8), intent(in) :: forc_psrf      ! atmosphere pressure at the surface [pa]
   real(r8), intent(in) :: cur            ! Charnock parameter (-)
   real(r8), intent(in) :: ustar          ! u* in similarity theory [m/s]

!  ------------------------- "out" arguments ---------------------------
   real(r8), intent(out) :: z0mg          ! roughness length over ground, momentum [m]
   real(r8), intent(out) :: z0hg          ! roughness length over ground, sensible heat [m]
   real(r8), intent(out) :: z0qg          ! roughness length over ground, latent heat [m]

! ------------------------- local variables ---------------------------
   real(r8) :: kva                        ! kinematic viscosity of air at ground temperature and forcing pressure
   real(r8) :: sqre0                      ! root of roughness Reynolds number
!================================================================================
      IF (t_grnd > tfrz .and. lktmptop > tfrz .and. snl == 0) THEN
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
      ELSE                          ! USE roughness over snow
         z0mg = 0.0024             ! z0mg won't have changed
         z0hg = z0mg/exp(0.13 * (ustar*z0mg/1.5e-5)**0.45)
         z0qg = z0hg
      ENDIF

   END SUBROUTINE roughness_lake



   SUBROUTINE hConductivity_lake(nl_lake,snl,t_grnd,&
                                z_lake,t_lake,lake_icefrac,rhow,&
                                rlat,ustar,z0mg,lakedepth, depthcrit, gamma, tk_lake, savedtke1)

! ------------------------ code history ---------------------------
! Description:
!     Calculate the thermal conductivity of the lake layers.
!     Diffusivity and implied thermal "conductivity" = diffusivity * cwat
!     old : tau = md * (Ke + Ked + Km) * Cwat * Rhow
!     fixbug: tau = [md * (Ke + Ked) + Km] * Cwat * Rhow
!     [Ke ] Represents the diffusivity of large vortices driven by wind
!     [Ked] Represents the enhanced diffusivity resulting from some unexpressible mixing process
!     [Km ] Represents the molecular diffusivity of liquid water
!     [md ] Represents an enhancement factor dependent on lake depth
!     [Cwat] Represents the specific heat capacity of water
!
! Original author:
!     Yongjiu Dai, 2000
!
! Revisions:
!     Nan Wei, 2018: added new mixing scheme for lake thermal diffusivity
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024: 
!                    Fixed the issue where the molecular diffusivity of liquid water was incorrectly enhanced.
!                    new: tau = [md * (Ke + Ked) + Km] * Cwat * Rhow
! -------------------------------------------------------------------------
   USE MOD_Lake_Const, only : tfrz, tkwat, tkice, tkair, p0, vonkar, &
                              grav, cpliq, cpice, cpair, denh2o, denice
!================================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in) :: nl_lake                ! number of soil layers
   integer, intent(in) :: snl                    ! number of snow layers
   real(r8), intent(in) :: t_grnd                ! ground surface temperature [k]
   real(r8), intent(in) :: z_lake(nl_lake)       ! lake node depth (middle point of layer) (m)
   real(r8), intent(in) :: t_lake(nl_lake)       ! lake temperature (kelvin)
   real(r8), intent(in) :: lake_icefrac(nl_lake) ! lake mass fraction of lake layer that is frozen
   real(r8), intent(in) :: rhow(nl_lake)         ! density of water (kg/m**3)
   real(r8), intent(in) :: rlat                  ! latitude (radians)
   real(r8), intent(in) :: ustar                 ! u* in similarity theory [m/s]
   real(r8), intent(in) :: z0mg                  ! roughness length over ground, momentum [m]
   real(r8), intent(in) :: lakedepth             ! column lake depth (m)
   real(r8), intent(in) :: depthcrit             ! (m) Depth beneath which to enhance mixing
   real(r8), intent(in) :: gamma                 ! Mixing enhancement factor

!  ------------------------- output variables ---------------------------
   real(r8), intent(out) :: tk_lake(nl_lake)     ! thermal conductivity at layer node [W/(m K)]
   real(r8), intent(out) :: savedtke1            ! top level eddy conductivity (W/mK)

!  ------------------------- local variables ---------------------------
   real(r8) :: kme(nl_lake)                      ! molecular + eddy diffusion coefficient (m**2/s)
   real(r8) :: cwat                              ! specific heat capacity of water (j/m**3/kelvin)
   real(r8) :: den                               ! used in calculating ri
   real(r8) :: drhodz                            ! d [rhow] /dz (kg/m**4)
   real(r8) :: fangkm                            ! (m^2/s) extra diffusivity based on Fang & Stefan 1996
   real(r8) :: ke                                ! eddy diffusion coefficient (m**2/s)
   real(r8) :: km                                ! molecular diffusion coefficient (m**2/s)
   real(r8) :: ks                                ! coefficient for calculation of decay of eddy diffusivity with depth
   real(r8) :: n2                                ! brunt-vaisala frequency (/s**2)
   real(r8) :: num                               ! used in calculating ri
   real(r8) :: ri                                ! richardson number
   real(r8) :: tkice_eff                         ! effective conductivity since layer depth is constant
   real(r8) :: tmp                               !
   real(r8) :: u2m                               ! 2 m wind speed (m/s)
   real(r8) :: ws                                ! surface friction velocity (m/s)
   integer  :: j                                 ! indices
   real(r8) :: mixfact = 5.                      ! Mixing enhancement factor.
!================================================================================

      mixfact = 5.0 * gamma ! Mixing enhancement factor. gamma is a parameter in the namelist

      cwat = cpliq*denh2o
      tkice_eff = tkice * denice/denh2o ! effective conductivity since layer depth is constant
      km = tkwat/cwat                   ! a constant (molecular diffusivity)
      u2m = max(0.1,ustar/vonkar*log(2./z0mg))
      ws = 1.2e-03 * u2m
      ks = 6.6 * sqrt( abs(sin(rlat)) ) * (u2m**(-1.84))

      DO j = 1, nl_lake-1
         drhodz = (rhow(j+1)-rhow(j)) / (z_lake(j+1)-z_lake(j))
         n2 = max(7.5e-5, grav / rhow(j) * drhodz)
         num = 40. * n2 * (vonkar*z_lake(j))**2
         tmp = -2.*ks*z_lake(j)        ! to avoid underflow computing
         IF(tmp < -40.) tmp = -40.     !
         den = max( (ws**2) * exp(tmp), 1.e-10 )
         ri = ( -1. + sqrt( max(1.+num/den, 0.) ) ) / 20.

         IF ((t_grnd > tfrz .and. t_lake(1) > tfrz .and. snl == 0) ) THEN
            tmp = -ks*z_lake(j)        ! to avoid underflow computing
            IF(tmp < -40.) tmp = -40.  !
            ke = vonkar*ws*z_lake(j)/p0 * exp(tmp) / (1.+37.*ri*ri)
            kme(j) = km + ke

            fangkm = 1.039e-8_r8 * max(n2,7.5e-5)**(-0.43)  ! Fang & Stefan 1996, citing Ellis et al 1991
            kme(j) = kme(j) + fangkm

            IF (lakedepth >= depthcrit) THEN
               kme(j) = kme(j) * mixfact    ! Mixing enhancement factor for lake deep than 25m.
            ENDIF
            tk_lake(j) = kme(j)*cwat
         ELSE
            kme(j) = km
            fangkm = 1.039e-8 * max(n2,7.5e-5)**(-0.43)
            kme(j) = kme(j) + fangkm
            IF (lakedepth >= depthcrit) THEN
               kme(j) = kme(j) * mixfact
            ENDIF
            tk_lake(j) = kme(j)*cwat*tkice_eff / ((1.-lake_icefrac(j))*tkice_eff &
                        + kme(j)*cwat*lake_icefrac(j))
         ENDIF
      ENDDO

      kme(nl_lake) = kme(nl_lake-1)
      savedtke1 = kme(1)*cwat

      IF ((t_grnd > tfrz .and. t_lake(1) > tfrz .and. snl == 0) ) THEN
         tk_lake(nl_lake) = tk_lake(nl_lake-1)
      ELSE
         tk_lake(nl_lake) = kme(nl_lake)*cwat*tkice_eff / ( (1.-lake_icefrac(nl_lake))*tkice_eff &
                          + kme(nl_lake)*cwat*lake_icefrac(nl_lake) )
      ENDIF

   END SUBROUTINE hConductivity_lake


    
   SUBROUTINE tridia (n, a, b, c, r, u)
!================================================================================
!  ------------------------- input variables ---------------------------
   integer,  intent(in) :: n       ! length of diagonal element vector
   real(r8), intent(in) :: a(1:n)  ! subdiagonal elements
   real(r8), intent(in) :: b(1:n)  ! diagonal elements
   real(r8), intent(in) :: c(1:n)  ! superdiagonal elements
   real(r8), intent(in) :: r(1:n)  ! right hand side

!  ------------------------- output variables ---------------------------
   real(r8), intent(out) :: u(1:n) ! solution vector

!  ------------------------- local variables ---------------------------
   integer  :: j                    ! indices  
   real(r8) :: gam(1:n)
   real(r8) :: bet
!================================================================================

      bet = b(1)
      u(1) = r(1) / bet
      DO j = 2, n
         gam(j) = c(j-1) / bet
         bet = b(j) - a(j) * gam(j)
         u(j) = (r(j) - a(j)*u(j-1)) / bet
      ENDDO
      DO j = n-1, 1, -1
         u(j) = u(j) - gam(j+1) * u(j+1)
      ENDDO

   END SUBROUTINE tridia

END MODULE  MOD_Lake_CoLML
