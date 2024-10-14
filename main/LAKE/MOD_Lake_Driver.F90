#include <define.h>

MODULE  MOD_Lake_Driver

! ================================== code history =====================================
! Description:
!     This is the new version of the CWRF lake simulation driver, which includes a total of 7 lake models.
!     The CALL of the lake model is controlled by <lkopt>:
!         (1) CoLML, Yongjiu Dai et al.       
!         (2) FLake, Dmitrii Mironov et al.   
!         (3) Simstrat, Goudsmit et al.       
!         (4) XOML, Tiejun Ling et al.        
!
!     -*- list of key output variables from Lake Driver -------------------------------------------------
!     ---------------------------------------------------------------------------------------------------      
!              Lake model:    COLML     FLake     Simstrat      XOML      note
!                   tmsno:        x        io           io         x
!                   tmice:        d        io           io         x
!                   lktmp:       io         o           io        io
!                   tmmnw:        x        io            x         x
!                   tmwml:        x        io            x         x                  
!                   tmbot:        x        io            x         x
!                   tmups:        x        io            x         x
!                   t_ssb:       io        io          sio         x      
!                   xwliq:       io         x            x         x
!                   xwice:       io         x            x         x                                                 
!                   icefr:       io         o            o         x                
!                   stke1:       io         i            i         i              
!                   icedp:       io        io           io         x                           
!                  bicedp:        x         x           io         x
!                  wicedp:        x         x           io         x
!                   snwdp:       io        io           io         x
!                   snwml:        o         o            o         x
!                   uwatv:        x         x           io        io     
!                   vwatv:        x         x           io        io
!                   lksal:        x         x           io        io
!                     tke:        x         x           io        io
!                    etke:       io         x           io         x        
!                     eps:        x         x           io        io
!                     num:        x         x           io         x
!                     nuh:        x         x           io         x
!                      Km:        x         x            x        io
!                      Kh:        x         x            x        io
!                      Ke:        x         x            x        io
!                   lkrho:        x         x           io        io
!                  rhosnw:        x         x           io         x
!                    mldp:        x        io            x         x
!                   upsdp:        x        io            x         x
!                  CTfrac:        x        io            x         x
!                   -----:   o=output    d=diagnostic output
!                            i=input     x=not set nor used (=spv_water)
!                            s=only snow layers otherwise [snow + sediment] layers
!                            r=required when restarting
!     ---------------------------------------------------------------------------------------------------      
!     -*- list of key output variables from Lake Driver -------------------------------------------------
!
!-lxz  TODO: make the whole lakectl to consider ifdef FrcICE
!-lxz  TODO: combine with seaice using xice to link with icefr(1)
!-WMEJ TODO: check all variable comments.
!-WMEJ TODO: organize variables in the order of input, output, and inout.
!-WMEJ TODO: add the necessary variables.
!-WMEJ TODO: review the variable names and comments.
! * WMEJICE: Options for calculating the ice fraction in the lake model. [OFF by default]
!            IF defined, the mass fraction will be used to calculate the ice fraction.
!
! Original author: 
!     Yongjiu Dai, 2000
!
! Revisions:
!     Min Xu, 02/2012
!     Xin-Zhong Liang,
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024 
!
! =====================================================================================

   USE MOD_Precision
   USE MOD_SPMD_Task
   IMPLICIT NONE

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: Lake_Driver

   ! PRIVATE MEMBER FUNCTIONS:
   PRIVATE :: LakeCtl

!--------------------------------------------------------------------------
CONTAINS
!--------------------------------------------------------------------------
   SUBROUTINE Lake_Driver( &
            ! "in" arguments    
            ! -------------------
            nlake     , nsnow     , nsoil     , nlice     ,& 
            dtlak     , rlat      , rlon      ,&
            hwref     , htref     , hqref     , usurf     ,&
            vsurf     , tmref     , qmref     , arhos     ,&
            psurf     , lwdns     , sabg      , zcbcv     ,& 
            crain     , lrain     , csnow     , lsnow     ,&
            hpbl      , vf_quartz , vf_gravels, vf_om     ,&
            vf_sand   , wf_gravels, wf_sand   , porsl     ,&
            csol      , k_solids  , dksatf    , dkdry     ,&
            dksatu    , BA_alpha  , BA_beta   , tprec     ,&
            sols      , soll      , solsd     , solld     ,&
            ipatch    , bifall    ,&
            ! "inout" arguments
            ! -------------------
            tskin     , lktmp     , t_ssb     ,&
            zlake     , zilak     , dzlak     , dplak     ,&
            zssb      , zissb     , dzssb     ,&
            stke1     , snwcv     , snwag     ,&
            snwdp     , xwliq     , xwice     , z0m       ,&
            z0h       , z0q       , felak     , gamma     ,&
            etal      , btpri     , icefr     , snlay     ,&
            tmsno     , tmice     , tmmnw     , tmwml     ,&
            tmbot     , tmups     , icedp     , mldp      ,&
            upsdp     , CTfrac    , frlak     , ziarea    ,&
            uwatv     , vwatv     , lksal     , tke       ,&
            eps       , etke      , num       , nuh       ,&
            bicedp    , wicedp    , lkrho     , rhosnw    ,&
            snout     ,&
! SNICAR model variables
            forc_aer  , sabg_lyr  , snofrz    ,&
            mss_bcpho , mss_bcphi , mss_ocpho , mss_ocphi ,&
            mss_dst1  , mss_dst2  , mss_dst3  , mss_dst4  ,&
! END SNICAR model variables
            ! "out" arguments
            ! -------------------
            fsena     , fevpa     , lfevpa    , fseng     ,&
            fevpg     , lwups     , fgrnd     , trad      ,&
            qseva     , qsubl     , qsdew     , qfros     ,&
            taux      , tauy      , ustar     , qstar     ,&
            tstar     , lkems     , snwml     , zol       ,&
            t2m       , q2m       , fm        , fq        ,&
            rib       , fh        , urban_call)
                  
! ---------------------------------- code history -------------------------------------
! Description:
!     This SUBROUTINE is responsible for managing the energy exchange between
!     lakes and the atmosphere, and for preparing certain variables for the lake models.
!     All variable processing related to the lakes should be completed here.
!
! Called: (* means optional)
!    -> LakeCtl             : Lake Model Controller
!
! Original author: 
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024
!
! -------------------------------------------------------------------------------------
   USE MOD_Lake_Const, only: rflk_depth_bs_ref, d2r, ShallowDeepPartition, SHALLOW, DEEP, PI
   USE MOD_Namelist, only: DEF_LAKE_MODEL_SCHEME, DEF_LAKE_ROUGHNESSS_SCHEME, DEF_LAKE_BETAPRIME_SCHEME,&
                           DEF_LAKE_FETCH_SCHEME, DEF_LAKE_ETA_SCHEME
   USE MOD_Vars_Global, only: spval
!==============================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in)      :: &
      ipatch                  ,&! patch index
      nlake                   ,&! Maximum number of lake layer
      nsnow                   ,&! Maximum number of snow layer (positive value)
      nsoil                   ,&! Maximum number of soil layer
      nlice                     ! Maximum number of ice layer

   real(r8), intent(in)     :: &
      dtlak                   ,&! time step [s]
      rlat                    ,&! latitude [radians]
      rlon                    ,&! longitude [radians]
      hwref                   ,&! observational height of wind [m]
      htref                   ,&! observational height of temperature [m]
      hqref                   ,&! observational height of humidity [m]
      usurf                   ,&! wind component in eastward direction [m/s]
      vsurf                   ,&! wind component in northward direction [m/s]
      tmref                   ,&! temperature at the reference height [K]
      qmref                   ,&! specific humidity at the reference height [kg/kg]
      arhos                   ,&! air density [kg/m3]
      psurf                   ,&! surface pressure [Pa]
      sabg                    ,&! solar absorbed by ground, the reflection has been removed [W/m2]
      sols                    ,&! atm vis direct beam solar rad onto srf [W m-2]
      soll                    ,&! atm nir direct beam solar rad onto srf [W m-2]
      solsd                   ,&! atm vis diffuse solar rad onto srf [W m-2]
      solld                   ,&! atm nir diffuse solar rad onto srf [W m-2]
      lwdns                   ,&! downward longwave radiation at surface [w/m2]
      hpbl                    ,&! atmospheric boundary layer height, only USE in CoLML_SurFlux [m]
      zcbcv                   ,&! convective boundary height [m]
      tprec                   ,&! snowfall/rainfall temperature [kelvin]
      bifall                  ,&! bulk density of newly fallen dry snow [kg/m3]
      crain                   ,&! convective rain [mm/s, kg/(m2 s)]
      lrain                   ,&! large scale rain [mm/s, kg/(m2 s)]
      csnow                   ,&! convective snow [mm/s, kg/(m2 s)]
      lsnow                     ! large scale snow [mm/s, kg/(m2 s)]

   real(r8), intent(in)     :: &
      porsl     (nsoil)       ,&! soil porosity [-]
      dksatu    (nsoil)       ,&! thermal conductivity of saturated unfrozen soil [W/m-K]
      k_solids  (nsoil)       ,&! thermal conductivity of mineralssoil [W/m-K]
      dkdry     (nsoil)       ,&! thermal conductivity of dry soil [W/m-K]
      vf_quartz (nsoil)       ,&! volumetric fraction of quartz within mineral soil [-]
      vf_gravels(nsoil)       ,&! volumetric fraction of gravels [-]
      vf_om     (nsoil)       ,&! volumetric fraction of organic matter [-]
      vf_sand   (nsoil)       ,&! volumetric fraction of sand [-]
      wf_gravels(nsoil)       ,&! gravimetric fraction of gravels [-]
      wf_sand   (nsoil)       ,&! gravimetric fraction of sand [-]
      csol      (nsoil)       ,&! heat capacity of soil solids [J/(m3 K)]
      dksatf    (nsoil)       ,&! thermal conductivity of saturated frozen soil [W/m-K]
      BA_alpha  (nsoil)       ,&! alpha in Balland and Arp(2005) thermal conductivity scheme
      BA_beta   (nsoil)         ! beta in Balland and Arp(2005) thermal conductivity scheme

!  ------------------------- inout variables ---------------------------
   integer, intent(inout)   :: &
      snlay                     ! number of snow layers

   real(r8), intent(inout)  :: &
      dzlak(nlake)            ,&! lake layer thickness [m]
      zlake(nlake)            ,&! lake layer node depth [m]
      zilak(nlake+1)          ,&! lake layer interface level [m]
      ziarea(nlake+1)         ,&! lake layer interface area [m2]
      lktmp(nlake)            ,&! lake temperature (K)
      uwatv(nlake)            ,&! Water velocity in x-direction [m/s]
      vwatv(nlake)            ,&! Water velocity in y-direction [m/s]
      lksal(nlake)            ,&! Salinity [‰]
      lkrho(nlake)            ,&! Water density [kg/m^3]
      icefr(nlake)            ,&! lake ice fraction [-]
      tke(nlake+1)            ,&! Turbulent kinetic energy (TKE) [J/kg]
      eps(nlake+1)            ,&! TKE dissipation rate [W/kg]
      num(nlake+1)            ,&! Turbulent viscosity (momentum)
      nuh(nlake+1)              ! Turbulent diffusivity (heat)


   real(r8), intent(inout)  :: &
      dzssb (-nsnow+1:nsoil)  ,&! snow + soil + bedrock layer thickness [m]
      zssb  (-nsnow+1:nsoil)  ,&! snow + soil + bedrock node layer depth [m]
      zissb (-nsnow:nsoil+1)  ,&! snow + soil + bedrock layer interface level [m]
      t_ssb (-nsnow+1:nsoil)  ,&! snow + soil + bedrock layer temperature [K]
      xwliq (-nsnow+1:nsoil)  ,&! snow + soil + bedrock layer liquid water (kg/m2)
      xwice (-nsnow+1:nsoil)    ! snow + soil + bedrock layer ice lens (kg/m2)
        
   real(r8), intent(inout)  :: &
      dplak                   ,&! lake depth [m]
      snwcv                   ,&! snow mass (kg/m2)
      snwag                   ,&! non dimensional snow age [-]
      snwdp                   ,&! snow depth [m]
      snout                   ,&! rate of water out of snow bottom (mm/s)
      tmsno                   ,&! snow temperature [K]
      tmice                   ,&! Temperature at the snow-ice or air-ice interface [K]
      tmmnw                   ,&! Mean temperature of the water column [K]
      tmwml                   ,&! Mixed-layer temperature [K]
      tmbot                   ,&! Temperature at the water-bottom sediment interface [K]
      tmups                   ,&! Temperature at the bottom of the upper layer of the sediments [K]
      mldp                    ,&! mixed layer depth [m]
      upsdp                   ,&! bottom of the upper layer of the sediments [m]
      icedp                   ,&! ice depth [m]
      bicedp                  ,&! black snow depth (m)
      wicedp                  ,&! white snow depth (m)
      rhosnw                  ,&! snow density (kg/m3)
      stke1                   ,&! top level eddy conductivity [W/m/K]
      etke                    ,&! Seiche energy [J]
      CTfrac                  ,&! Shape factor (thermocline)
      tskin                   ,&! ground surface temperature [k]
      z0m                     ,&! roughness length over ground, momentum [m]
      z0h                     ,&! roughness length over ground, sensible heat [m]
      z0q                     ,&! roughness length over ground, latent heat [m]
      felak                   ,&! lake fetch (>0), otherwise USE emperical eq by dplak [m]
      gamma                   ,&! Mixing enhancement factorfor eddy diffusion coefficient, 1 means no enhancement [-]
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
      lwups                   ,&! outgoing long-wave radiation from ground+canopy
      fgrnd                   ,&! ground heat flux [W/m2]
      qseva                   ,&! ground surface evaporation rate (mm h2o/s)
      qsdew                   ,&! ground surface dew formation (mm h2o /s) [+]
      qsubl                   ,&! sublimation rate from snow pack (mm h2o /s) [+]
      qfros                   ,&! surface dew added to snow pack (mm h2o /s) [+]
      taux                    ,&! wind stress: E-W [kg/m/s**2]
      tauy                    ,&! wind stress: N-S [kg/m/s**2]
      ustar                   ,&! u* in similarity theory [m/s]
      qstar                   ,&! q* in similarity theory [kg/kg]
      tstar                   ,&! t* in similarity theory [K]
      lkems                   ,&! averaged bulk surface emissivity
      zol                     ,&! dimensionless height (z/L) used in Monin-Obukhov theory
      snwml                   ,&! snowmelt rate [mm/s]
      rib                     ,&! bulk Richardson number in surface layer
      t2m                     ,&! 2 m height air temperature [K]
      q2m                     ,&! 2 m height air specific humidity
      fm                      ,&! integral of profile function for momentum
      fq                      ,&! integral of profile function for moisture
      fh                      ,&! integral of profile function for heat
      trad                      ! radiative temperature [K]
        
   logical, optional, intent(in) :: urban_call   ! whether it is a urban CALL

!  ------------------------- local variables ---------------------------
   integer                  :: &
      lkopt                   ,&! lake models: 1) CoLML, 2) FLake, 3) Simstrat, 4) XOML
      zopt                    ,&! lake surface roughness : 1) external constant, 2) wind-dependent(Subin 2012)
      betaopt                 ,&! fraction of solar radiation in the NIR : 1) external constant, 2) equation(CoLM)
      fetchopt                ,&! lake fetch : 1) external constant, 2) lake depth-dependent(Subin)
      etaopt                  ,&! extinction coefficient : 1) external constant, 2) lake depth-dependent(Subin)
      scwat                   ,&! surface category of water characteristics: (1) shallow lake, (2) deep lake
      j                       ,&! loop counter
      i                         ! loop counter

   real(r8)                 :: &
      xlat                    ,&! latitude [degrees]
      xlon                    ,&! longitude [degrees]
      dpsed                   ,&! sediment depth [m]
      T_bs                    ,&! bottom sediment temperature condition [K]
      wdm                     ,&! lowest level wind speed including the stablity effect [m/s]
      ram                     ,&! aerodynamical resistance [s/m]
      rah                     ,&! thermal resistance [s/m]
      raw                     ,&! moisture resistance [s/m]
      u10m                    ,&! 10 m height wind speed in eastward direction [m/s]
      v10m                    ,&! 10 m height wind speed in northward direction [m/s]
      fm10m                   ,&! integral of profile function for momentum
      fq10m                   ,&! integral of profile function for moisture
      fh2m                    ,&! integral of profile function for heat
      fq2m                    ,&! integral of profile function for moisture
      rnird                   ,&! reflected direct beam nir solar radiation (W/m**2)
      rniri                   ,&! reflected diffuse nir solar radiation (W/m**2]
      shfdt                     ! derivative of srf sen+lat  heat flux wrt srf temp [W/m2/K]
!================================================================================
      !- Solar radiation related (not used)
      rnird = 0. ! reflected direct beam nir solar radiation (W/m**2) 
      rniri = 0. ! reflected diffuse nir solar radiation (W/m**2)

      T_bs  = t_ssb(nsoil)         ! +WMEJ Assume the temperature is consistent with the soil layers of other models. 
      ! dpsed = zissb(nsoil+1)       ! +WMEJ Assume the thickness is consistent with the soil layers of other models.
      dpsed = rflk_depth_bs_ref  ! +WMEJ The default recommended thickness in the model is 10m., only used in CoLML
      IF (dplak <= ShallowDeepPartition) THEN
         scwat = SHALLOW
      ELSE
         scwat = DEEP
      ENDIF

      !- Convert latitude and longitude to degrees
      xlat = rlat/d2r
      xlon = rlon/d2r

      !- Lake model controller
      !   Avoid using the USE MODULE  method to import, 
      !   as it may need to be adjusted for different lakes in the future,
      !   and is not conducive to subsequent testing.
      lkopt = DEF_LAKE_MODEL_SCHEME
      zopt = DEF_LAKE_ROUGHNESSS_SCHEME
      betaopt = DEF_LAKE_BETAPRIME_SCHEME
      fetchopt = DEF_LAKE_FETCH_SCHEME
      etaopt = DEF_LAKE_ETA_SCHEME

      CALL LakeCtl ( &
            ! "in" arguments    
            ! -------------------
            nlake     , nsnow     , nsoil     , nlice     ,& 
            scwat     , dtlak     , xlat      , xlon      , dplak     ,& 
            lkopt     , zopt      , betaopt   , fetchopt  , etaopt    ,& 
            hwref     , htref     , hqref     , usurf     , vsurf     ,& 
            tmref     , qmref     , arhos     , psurf     , lwdns     ,&
            sols      , soll      , solsd     , solld     , sabg      ,&
            rnird     , rniri     , zcbcv     , hpbl      , tprec     ,& 
            crain     , lrain     , csnow     , lsnow     , dpsed     ,&           
            vf_quartz , vf_gravels, vf_om     , vf_sand   , wf_gravels,&
            wf_sand   , porsl     , csol      , k_solids  , dksatf    ,&
            dkdry     , dksatu    , BA_alpha  , BA_beta   , T_bs      ,&
            ipatch    , bifall    ,&
            ! "inout" arguments
            ! -------------------
            zlake     , zilak     , dzlak     , lktmp     , tskin     ,& 
            dzssb     , zssb      , zissb     , t_ssb     , stke1     ,& 
            snwcv     , snwag     , snwdp     , xwliq     , xwice     ,& 
            icefr     , snlay     , icedp     , etal      , btpri     ,&
            z0m       , z0h       , z0q       , felak     , gamma     ,&
            tmsno     , tmice     , tmmnw     , tmwml     , tmbot     ,&
            tmups     , mldp      , upsdp     , CTfrac    , frlak     ,&
            ziarea    , uwatv     , vwatv     , lksal     , tke       ,&
            eps       , etke      , num       , nuh       , bicedp    ,&
            wicedp    , lkrho     , rhosnw    , snout     ,&
! SNICAR model variables
            forc_aer  , sabg_lyr  , snofrz    ,&
            mss_bcpho , mss_bcphi , mss_ocpho , mss_ocphi ,&
            mss_dst1  , mss_dst2  , mss_dst3  , mss_dst4  ,&
! END SNICAR model variables
            ! "out" arguments
            ! -------------------
            fsena     , fevpa     , lfevpa    , fseng     , fevpg     ,&
            qseva     , qsubl     , qsdew     , qfros     ,&
            lwups     , fgrnd     , trad      , taux      , tauy      ,&
            ustar     , qstar     , tstar     , lkems     , snwml     ,&
            zol       , rib       , ram       , rah       , raw       ,&
            wdm       , t2m       , q2m       , u10m      , v10m      ,&
            fm10m     , fq10m     , fh2m      , fq2m      , fm        ,&
            fq        , fh        , shfdt     )

   END SUBROUTINE Lake_Driver



   SUBROUTINE LakeCtl ( &
               ! "in" arguments    
               ! -------------------
               nlake     , nsnow     , nsoil     , nlice     ,& 
               scwat     , dtlak     , xlat      , xlon      , dplak     ,& 
               lkopt     , zopt      , betaopt   , fetchopt  , etaopt    ,& 
               hwref     , htref     , hqref     , usurf     , vsurf     ,& 
               tmref     , qmref     , arhos     , psurf     , lwdns     ,&
               sols      , soll      , solsd     , solld     , sabg      ,&
               rnird     , rniri     , zcbcv     , hpbl      , tprec     ,& 
               crain     , lrain     , csnow     , lsnow     , dpsed     ,&
               vf_quartz , vf_gravels, vf_om     , vf_sand   , wf_gravels,&
               wf_sand   , porsl     , csol      , k_solids  , dksatf    ,&
               dkdry     , dksatu    , BA_alpha  , BA_beta   , T_bs      ,&
               ipatch    , bifall    ,&
               ! "inout" arguments
               ! -------------------
               zlake     , zilak     , dzlak     , lktmp     , tskin     ,& 
               dzssb     , zssb      , zissb     , t_ssb     , stke1     ,& 
               snwcv     , snwag     , snwdp     , xwliq     , xwice     ,& 
               icefr     , snlay     , icedp     , etal      , btpri     ,&
               z0m       , z0h       , z0q       , felak     , gamma     ,&
               tmsno     , tmice     , tmmnw     , tmwml     , tmbot     ,&
               tmups     , mldp      , upsdp     , CTfrac    , frlak     ,&
               ziarea    , uwatv     , vwatv     , lksal     , tke       ,&
               eps       , etke      , num       , nuh       , bicedp    ,&
               wicedp    , lkrho     , rhosnw    , snout     ,&
! SNICAR model variables
               forc_aer  , sabg_lyr  , snofrz    ,&
               mss_bcpho , mss_bcphi , mss_ocpho , mss_ocphi ,&
               mss_dst1  , mss_dst2  , mss_dst3  , mss_dst4  ,&
! END SNICAR model variables
               ! "out" arguments
               ! -------------------
               fsena     , fevpa     , lfevpa    , fseng     , fevpg     ,&
               qseva     , qsubl     , qsdew     , qfros     ,&
               lwups     , fgrnd     , trad      , taux      , tauy      ,&
               ustar     , qstar     , tstar     , lkems     , snwml     ,&
               zol       , rib       , ram       , rah       , raw       ,&
               wdm       , t2m       , q2m       , u10m      , v10m      ,&
               fm10m     , fq10m     , fh2m      , fq2m      , fm        ,&
               fq        , fh        , shfdt      ,urban_call)

! ---------------------------------- code history -------------------------------------
! Description:
!     Lake Model Controller. This should be used as an isolator between the lake models
!
! *************************************************************************************  
! TO ENSURE CODE READABILITY, DO NOT PERFORM ANY CALCULATIONS OR MACRO DEFINITIONS HERE.
! *************************************************************************************  
!
! Called: (* means optional)
!    *-> Lake_CoLML        : CoLM-Lake 
!    *-> Lake_FLake        : FLake 
!    *-> Lake_SIMSTRAT     : Simstrat
!    *-> Lake_XOML         : XOML
!
! Original author:
!     Min Xu, 02/2012
!
! Revisions:
!     Xin-Zhong Liang, 
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024: Reconstructed the lake model driver
! -------------------------------------------------------------------------------------
   USE MOD_Lake_Const, only: COLML, FLAKE, SIMSTRAT, XOML
   USE MOD_Lake_CoLML, only: Lake_CoLML
   USE MOD_Lake_FLake, only: Lake_FLake
   USE MOD_Lake_Simstrat, only: Lake_Simstrat
!================================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in)      :: &
      ipatch                  ,&! patch index
      lkopt                   ,&! lake models: 1) CoLML, 2) FLake, 3) Simstrat, 4) XOML
      nlake                   ,&! Maximum number of lake layer
      nsnow                   ,&! Maximum number of snow layer
      nsoil                   ,&! Maximum number of soil layer
      nlice                   ,&! Maximum number of ice layer, not used in this version
      scwat                   ,&! water body type 
      zopt                    ,&! lake surface roughness : 1) external constant, 2) wind-dependent(Subin 2012)
      betaopt                 ,&! fraction of solar radiation in the NIR : 1) external constant, 2) equation(CoLM)
      fetchopt                ,&! lake fetch : 1) external constant, 2) lake depth-dependent(Subin)
      etaopt                    ! lake option : 1) external constant, 2) lake depth-dependent(Subin)
        

   real(r8), intent(in)     :: &
      dtlak                   ,&! time step [s]
      xlat                    ,&! latitude [degree]
      xlon                    ,&! longitude [degree]
      dplak                   ,&! lake depth [m]
      dpsed                   ,&! sediment depth [m]
      T_bs                    ,&! bottom sediment temperature condition [K]
      hwref                   ,&! observational height of wind [m]
      htref                   ,&! observational height of temperature [m]
      hqref                   ,&! observational height of humidity [m]
      usurf                   ,&! wind component in eastward direction [m/s]
      vsurf                   ,&! wind component in northward direction [m/s]
      tmref                   ,&! temperature at the reference height [K]
      qmref                   ,&! specific humidity at the reference height [kg/kg]
      arhos                   ,&! air density [kg/m3]
      psurf                   ,&! surface pressure [Pa]
      sols                    ,&! atm vis direct beam solar rad onto srf [W m-2]
      soll                    ,&! atm nir direct beam solar rad onto srf [W m-2]
      solsd                   ,&! atm vis diffuse solar rad onto srf [W m-2]
      solld                   ,&! atm nir diffuse solar rad onto srf [W m-2]
      rnird                   ,&! reflected direct beam nir solar radiation (W/m**2)
      rniri                   ,&! reflected diffuse nir solar radiation (W/m**2)
      sabg                    ,&! solar absorbed by ground  [W/m2]
      lwdns                   ,&! downward longwave radiation at surface [w/m2]
      hpbl                    ,&! downward longwave radiation at surface [wm-2]
      zcbcv                   ,&! convective boundary height [m]
      tprec                   ,&! snowfall/rainfall temperature [kelvin]
      bifall                  ,&! bulk density of newly fallen dry snow [kg/m3]
      crain                   ,&! convective rain [mm/s, kg/(m2 s)]
      lrain                   ,&! large scale rain [mm/s, kg/(m2 s)]
      csnow                   ,&! convective snow [mm/s, kg/(m2 s)]
      lsnow                     ! large scale snow [mm/s, kg/(m2 s)]

   real(r8), intent(in)     :: &
      porsl     (nsoil)       ,&! soil porosity [-]
      dksatu    (nsoil)       ,&! thermal conductivity of saturated unfrozen soil [W/m-K]
      k_solids  (nsoil)       ,&! thermal conductivity of mineralssoil [W/m-K]
      dkdry     (nsoil)       ,&! thermal conductivity of dry soil [W/m-K]
      vf_quartz (nsoil)       ,&! volumetric fraction of quartz within mineral soil [-]
      vf_gravels(nsoil)       ,&! volumetric fraction of gravels [-]
      vf_om     (nsoil)       ,&! volumetric fraction of organic matter [-]
      vf_sand   (nsoil)       ,&! volumetric fraction of sand [-]
      wf_gravels(nsoil)       ,&! gravimetric fraction of gravels [-]
      wf_sand   (nsoil)       ,&! gravimetric fraction of sand [-]
      csol      (nsoil)       ,&! heat capacity of soil solids [J/(m3 K)]
      dksatf    (nsoil)       ,&! thermal conductivity of saturated frozen soil [W/m-K]
      BA_alpha  (nsoil)       ,&! alpha in Balland and Arp(2005) thermal conductivity scheme
      BA_beta   (nsoil)         ! beta in Balland and Arp(2005) thermal conductivity scheme

!  ------------------------- inout variables ---------------------------
   integer, intent(inout)   :: &
      snlay                     ! number of snow layers

   real(r8), intent(inout)  :: &
      dzlak(nlake)            ,&! lake layer thickness [m]
      zlake(nlake)            ,&! lake layer node depth [m]
      zilak(nlake+1)          ,&! lake layer interface level [m]
      ziarea(nlake+1)         ,&! lake layer interface area [m2]
      lktmp(nlake)            ,&! lake temperature (K)
      uwatv(nlake)            ,&! Water velocity in x-direction [m/s]
      vwatv(nlake)            ,&! Water velocity in y-direction [m/s]
      lksal(nlake)            ,&! Salinity [‰]
      lkrho(nlake)            ,&! Water density [kg/m^3]
      icefr(nlake)            ,&! lake ice fraction [-]
      tke(nlake+1)            ,&! Turbulent kinetic energy (TKE) [J/kg]
      eps(nlake+1)            ,&! TKE dissipation rate [W/kg]
      num(nlake+1)            ,&! Turbulent viscosity (momentum)
      nuh(nlake+1)              ! Turbulent diffusivity (heat)

   real(r8), intent(inout)  :: &
      dzssb (-nsnow+1:nsoil)  ,&! snow + soil + bedrock layer thickness [m]
      zssb  (-nsnow+1:nsoil)  ,&! snow + soil + bedrock node layer depth [m]
      zissb (-nsnow:nsoil+1)  ,&! snow + soil + bedrock layer interface level [m]
      t_ssb (-nsnow+1:nsoil)  ,&! snow + soil + bedrock layer temperature [K]
      xwliq (-nsnow+1:nsoil)  ,&! snow + soil + bedrock layer liquid water (kg/m2)
      xwice (-nsnow+1:nsoil)    ! snow + soil + bedrock layer ice lens (kg/m2)

   real(r8), intent(inout)  :: &
      snwcv                   ,&! snow mass (kg/m2)
      snwag                   ,&! non dimensional snow age [-]
      snwdp                   ,&! snow depth [m]
      snout                   ,&! rate of water out of snow bottom (mm/s)
      tmsno                   ,&! snow temperature [K]
      tmice                   ,&! Temperature at the snow-ice or air-ice interface [K]
      tmmnw                   ,&! Mean temperature of the water column [K]
      tmwml                   ,&! Mixed-layer temperature [K]
      tmbot                   ,&! Temperature at the water-bottom sediment interface [K]
      tmups                   ,&! Temperature at the bottom of the upper layer of the sediments [K]
      mldp                    ,&! mixed layer depth [m]
      upsdp                   ,&! bottom of the upper layer of the sediments [m]
      icedp                   ,&! ice depth [m]
      bicedp                  ,&! black snow depth (m)
      wicedp                  ,&! white snow depth (m)
      rhosnw                  ,&! snow density (kg/m3)
      stke1                   ,&! top level eddy conductivity [W/m/K]
      etke                    ,&! Seiche energy [J]
      CTfrac                  ,&! Shape factor (thermocline)
      tskin                   ,&! ground surface temperature [k]
      z0m                     ,&! roughness length over ground, momentum [m]
      z0h                     ,&! roughness length over ground, sensible heat [m]
      z0q                     ,&! roughness length over ground, latent heat [m]
      felak                   ,&! lake fetch (>0), otherwise USE emperical eq by dplak [m]
      gamma                   ,&! Mixing enhancement factorfor eddy diffusion coefficient
      etal                    ,&! extinction coefficient [1/m]
      btpri                   ,&! beta prime in Monin-Obukhov theory
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
      lwups                   ,&! outgoing long-wave radiation from ground+canopy
      fgrnd                   ,&! ground heat flux [W/m2]
      qseva                   ,&! ground surface evaporation rate (mm h2o/s)
      qsdew                   ,&! ground surface dew formation (mm h2o /s) [+]
      qsubl                   ,&! sublimation rate from snow pack (mm h2o /s) [+]
      qfros                   ,&! surface dew added to snow pack (mm h2o /s) [+]
      taux                    ,&! wind stress: E-W [kg/m/s**2]
      tauy                    ,&! wind stress: N-S [kg/m/s**2]
      ustar                   ,&! u* in similarity theory [m/s]
      qstar                   ,&! q* in similarity theory [kg/kg]
      tstar                   ,&! t* in similarity theory [K]
      lkems                   ,&! averaged bulk surface emissivity
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
      trad                    ,&! radiative temperature [K]
      shfdt                     ! derivative of srf sen+lat  heat flux wrt srf temp [W/m2/K]

   logical, optional, intent(in) :: urban_call   ! whether it is a urban CALL
!================================================================================

      ! *************************************************************************************
      SELECT CASE (lkopt)   !- Select the lake model  !
      ! *************************************************************************************
         !-----------------------
         CASE (COLML) !- CoLM-Lake
         !-----------------------
            ! print *, 'Lake Model: CoLM-Lake'
            CALL Lake_CoLML ( &
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
                  snlay     , snwcv     , snwag     , snwdp     ,&
                  stke1     , tskin     , xwliq     , xwice     ,&
                  z0m       , z0h       , z0q       , felak     ,&
                  gamma     , etal      , btpri     , frlak     ,&
                  icefr     , tmice     , snout     , lkrho     ,& 
! SNICAR model variables
                  forc_aer  , sabg_lyr  , snofrz    ,&
                  mss_bcpho , mss_bcphi , mss_ocpho , mss_ocphi ,&
                  mss_dst1  , mss_dst2  , mss_dst3  , mss_dst4  ,&
! END SNICAR model variables
                  ! "out" arguments
                  ! ---------------------------
                  fsena     , fevpa     , lfevpa    , fseng     ,&
                  fevpg     , lwups     , fgrnd     , trad      ,&
                  qseva     , qsubl     , qsdew     , qfros     ,&
                  taux      , tauy      , ustar     , qstar     ,&
                  tstar     , lkems     , zol       , snwml     ,&
                  rib       , ram       , rah       , raw       ,&
                  wdm       , t2m       , q2m       , u10m      ,&  
                  v10m      , fm10m     , fq10m     , fh2m      ,&
                  fq2m      , fm        , fq        , fh        ,&
                  shfdt     )

         !-----------------------
         CASE (FLAKE) !- FLake
         !-----------------------
            ! print *, 'Lake Model: FLake'
            CALL Lake_FLake ( &
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
                  snlay     , snwcv     , snwag     , snwdp     ,&
                  stke1     , tskin     , xwliq     , xwice     ,&
                  z0m       , z0h       , z0q       , felak     ,&
                  gamma     , etal      , btpri     , frlak     ,&
                  tmsno     , tmice     , tmmnw     , tmwml     ,&  
                  tmbot     , tmups     , icedp     , mldp      ,&  
                  upsdp     , CTfrac    , icefr     , &              
                  ! "out" arguments
                  ! ---------------------------
                  fsena     , fevpa     , lfevpa    , fseng     ,&
                  fevpg     , lwups     , fgrnd     , trad      ,&
                  qseva     , qsubl     , qsdew     , qfros     ,&
                  taux      , tauy      , ustar     , qstar     ,&
                  tstar     , lkems     , zol       , snwml     ,&
                  rib       , ram       , rah       , raw       ,&
                  wdm       , t2m       , q2m       , u10m      ,&
                  v10m      , fm10m     , fq10m     , fh2m      ,&
                  fq2m      , fm        , fq        , fh        ,&
                  shfdt     )

         !-----------------------
         CASE (SIMSTRAT) !- Simstrat
         !-----------------------
            ! print *, 'Lake Model: Simstrat'
            CALL Lake_Simstrat( &
                  ! "in" arguments    
                  ! -------------------
                  nlake     , nsnow     , scwat     , nsoil     ,&
                  zopt      , betaopt   , fetchopt  , etaopt    ,&
                  dtlak     , xlat      , xlon      , dplak     ,&
                  hwref     , htref     , hqref     , tmref     ,&
                  usurf     , vsurf     , qmref     , arhos     ,&
                  sols      , soll      , solsd     , solld     ,&
                  sabg      , lwdns     , psurf     , hpbl      ,&
                  crain     , csnow     , lrain     , lsnow     ,&
                  zcbcv     , ipatch    ,&
                  ! "inout" arguments
                  ! ---------------------------
                  dzlak     , zlake     , zilak     , lktmp     ,&
                  dzssb     , zssb      , zissb     , t_ssb     ,&
                  snlay     , snwcv     , snwag     , snwdp     ,&
                  stke1     , tskin     , xwliq     , xwice     ,&
                  z0m       , z0h       , z0q       , felak     ,&
                  gamma     , etal      , btpri     , frlak     ,&
                  icefr     , icedp     , tmice     , ziarea    ,&
                  uwatv     , vwatv     , lksal     , tke       ,&
                  eps       , etke      , num       , nuh       ,&
                  bicedp    , wicedp    , lkrho     , rhosnw    ,&
                  ! "out" arguments
                  ! ---------------------------
                  fsena     , fevpa     , lfevpa    , fseng     ,&
                  fevpg     , lwups     , fgrnd     , trad      ,&
                  qseva     , qsubl     , qsdew     , qfros     ,&
                  taux      , tauy      , ustar     , qstar     ,&
                  tstar     , lkems     , zol       , snwml     ,&
                  rib       , ram       , rah       , raw       ,&
                  wdm       , t2m       , q2m       , u10m      ,&
                  v10m      , fm10m     , fq10m     , fh2m      ,&
                  fq2m      , fm        , fq        , fh        ,&
                  shfdt     )

         !-----------------------
         CASE (XOML) !- XOML
         !-----------------------
            print *, '[Error] XOML is not implemented yet'
            CALL CoLM_stop ()
        END SELECT 

    END SUBROUTINE LakeCtl


END MODULE  MOD_Lake_Driver

