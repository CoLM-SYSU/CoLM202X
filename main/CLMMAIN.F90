#include <define.h>

SUBROUTINE CLMMAIN ( &
   
         ! 模型运行信息
           ipatch,       idate,        coszen,       deltim,        &
           patchlonr,    patchlatr,    patchclass,   patchtype,     &
           doalb,        dolai,        dosst,        oro,           &


         ! soil information and lake depth
           soil_s_v_alb, soil_d_v_alb, soil_s_n_alb, soil_d_n_alb,  &
           porsl,        psi0,         bsw,          hksati,        &
           csol,         dksatu,       dkdry,        rootfr,        &
           lakedepth,    dz_lake,                                   &

         ! vegetation information
           htop,         hbot,         sqrtdi,                      &
           effcon,       vmax25,       slti,         hlti,          &
           shti,         hhti,         trda,         trdm,          &
           trop,         gradm,        binter,       extkn,         &
           chil,         rho,          tau,                         &

         ! atmospheric forcing
           forc_pco2m,   forc_po2m,    forc_us,      forc_vs,       &
           forc_t,       forc_q,       forc_prc,     forc_prl,      &
           forc_rain,    forc_snow,    forc_psrf,    forc_pbot,     &
           forc_sols,    forc_soll,    forc_solsd,   forc_solld,    &
           forc_frl,     forc_hgt_u,   forc_hgt_t,   forc_hgt_q,    &
           forc_rhoair,                                             &

         ! land surface variables required for restart
           z_soisno,     dz_soisno,    t_soisno,     wliq_soisno,   &
           wice_soisno,  t_grnd,       tleaf,        ldew,          &
           sag,          scv,          snowdp,       fveg,          &
           fsno,         sigf,         green,        lai,           &
           sai,          alb,          ssun,         ssha,          &
           thermk,       extkb,        extkd,                       &

           zwt,          wa,                                        &
           t_lake,       lake_icefrac,                              &

         ! additional diagnostic variables for output
           laisun,       laisha,                                    &
           rstfac,       h2osoi,       wat,                         &

         ! FLUXES
           taux,         tauy,         fsena,        fevpa,         &
           lfevpa,       fsenl,        fevpl,        etr,           &
           fseng,        fevpg,        olrg,         fgrnd,         &
           trad,         tref,         qref,         rsur,          &
           rnof,         qintr,        qinfl,        qdrip,         &
           rst,          assim,        respc,        sabvsun,       &
           sabvsha,      sabg,         sr,           solvd,         &
           solvi,        solnd,        solni,        srvd,          &
           srvi,         srnd,         srni,         solvdln,       &
           solviln,      solndln,      solniln,      srvdln,        &
           srviln,       srndln,       srniln,       qcharge,       &
           xerr,         zerr,                                      &

         ! TUNABLE modle constants
           zlnd,         zsno,         csoilc,       dewmx,         &
           wtfact,       capr,         cnfac,        ssi,           &
           wimp,         pondmx,       smpmax,       smpmin,        &
           trsmx0,       tcrit,                                     & 

         ! additional variables required by coupling with WRF model 
           emis,         z0m,          zol,          rib,           &
           ustar,        qstar,        tstar,        fm,            &
           fh,           fq                                         )

!=======================================================================
!
! Main subroutine, advance time information
!
! Original author : Yongjiu Dai, 09/30/1999; 08/30/2002, 12/2012, 02/2014
!
!    FLOW DIAGRAM FOR CLMMAIN
!
!    CLMMAIN ===> netsolar                 |> all surface
!                 rain_snow_temp           !> all surface
!
!                 LEAF_interception        |]  
!                 newsnow                  |] patchtype = 0 (soil ground)
!                 THERMAL                  |]           = 1 (urban & built-up)
!                 WATER                    |]           = 2 (wetland)
!                 snowcompaction           |]           = 3 (land ice)
!                 snowlayerscombine        |]           = 4 (lake)
!                 snowlayersdivide         |]
!                 snowage                  |]  
!
!                 newsnow_lake             |]
!                 laketem                  |] lake scheme
!                 snowwater_lake           |]
!
!                 SOCEAN                   |> ocean and sea ice
!                
!                 orb_coszen               |> all surface
!                 EcoModel (LAI_empirical) |> land 
!                 snowfraction             |> land
!                 albland                  |> land 
!                 albocean                 |> ocean & sea ice
!
!=======================================================================

  USE precision
  USE PhysicalConstants, only: tfrz, denh2o, denice
  USE MOD_TimeInvariants, only: gridlond
  USE MOD_TimeVariables, only: tlai, tsai
  USE MOD_PFTimeInvars
  USE MOD_PFTimeVars
  USE MOD_PCTimeInvars
  USE MOD_PCTimeVars
  USE SOIL_SNOW_hydrology
  USE SNOW_Layers_CombineDivide
  USE GLACIER
  USE LAKE
  USE SIMPLE_OCEAN
  USE ALBEDO
  USE timemanager
 
  IMPLICIT NONE
 
! ------------------------ Dummy Argument ------------------------------
  REAL(r8),intent(in) :: deltim  !seconds in a time step [second]
  LOGICAL, intent(in) :: doalb   !true if time for surface albedo calculation
  LOGICAL, intent(in) :: dolai   !true if time for leaf area index calculation
  LOGICAL, intent(in) :: dosst   !true to update sst/ice/snow before calculation

  INTEGER, intent(in) :: & 
        ipatch        ! patch index

  REAL(r8), intent(in) :: &
        patchlonr   ,&! logitude in radians
        patchlatr     ! latitude in radians

  INTEGER, intent(in) :: & 
        patchclass  ,&! land cover type of USGS classification or others
        patchtype     ! land water type (0=soil, 1=urban and built-up, 
                      ! 2=wetland, 3=land ice, 4=land water bodies, 99 = ocean)
! Parameters
! ----------------------
  REAL(r8), intent(in) :: &
        lakedepth        ,&! lake depth (m)
        dz_lake(nl_lake) ,&! lake layer thickness (m)

        ! soil physical parameters and lake info
        soil_s_v_alb     ,&! albedo of visible of the saturated soil
        soil_d_v_alb     ,&! albedo of visible of the dry soil
        soil_s_n_alb     ,&! albedo of near infrared of the saturated soil
        soil_d_n_alb     ,&! albedo of near infrared of the dry soil
        porsl(nl_soil)   ,&! fraction of soil that is voids [-]
        psi0(nl_soil)    ,&! minimum soil suction [mm]
        bsw(nl_soil)     ,&! clapp and hornbereger "b" parameter [-]
        hksati(nl_soil)  ,&! hydraulic conductivity at saturation [mm h2o/s]
        csol(nl_soil)    ,&! heat capacity of soil solids [J/(m3 K)]
        dksatu(nl_soil)  ,&! thermal conductivity of saturated soil [W/m-K]
        dkdry(nl_soil)   ,&! thermal conductivity for dry soil  [J/(K s m)]
        rootfr(nl_soil)  ,&! fraction of roots in each soil layer

        ! vegetation static, dynamic, derived parameters
        htop        ,&! canopy top height [m]
        hbot        ,&! canopy bottom height [m]
        sqrtdi      ,&! inverse sqrt of leaf dimension [m**-0.5]
        effcon      ,&! quantum efficiency of RuBP regeneration (mol CO2/mol quanta)
        vmax25      ,&! maximum carboxylation rate at 25 C at canopy top
        slti        ,&! slope of low temperature inhibition function      [s3] 
        hlti        ,&! 1/2 point of low temperature inhibition function  [s4]
        shti        ,&! slope of high temperature inhibition function     [s1]
        hhti        ,&! 1/2 point of high temperature inhibition function [s2]
        trda        ,&! temperature coefficient in gs-a model             [s5]
        trdm        ,&! temperature coefficient in gs-a model             [s6]
        trop        ,&! temperature coefficient in gs-a model          
        gradm       ,&! conductance-photosynthesis slope parameter
        binter      ,&! conductance-photosynthesis intercep
        extkn       ,&! coefficient of leaf nitrogen allocation
        chil        ,&! leaf angle distribution factor
        rho(2,2)    ,&! leaf reflectance (iw=iband, il=life and dead)
        tau(2,2)    ,&! leaf transmittance (iw=iband, il=life and dead)

        ! tunable parameters
        zlnd        ,&! roughness length for soil [m]
        zsno        ,&! roughness length for snow [m]
        csoilc      ,&! drag coefficient for soil under canopy [-]
        dewmx       ,&! maximum dew
        wtfact      ,&! fraction of model area with high water table
        capr        ,&! tuning factor to turn first layer T into surface T
        cnfac       ,&! Crank Nicholson factor between 0 and 1
        ssi         ,&! irreducible water saturation of snow
        wimp        ,&! water impremeable if porosity less than wimp
        pondmx      ,&! ponding depth (mm)
        smpmax      ,&! wilting point potential in mm
        smpmin      ,&! restriction for min of soil poten.  (mm)
        trsmx0      ,&! max transpiration for moist soil+100% veg.  [mm/s]
        tcrit         ! critical temp. to determine rain or snow

! Forcing
! ----------------------
  REAL(r8), intent(in) :: &
        forc_pco2m  ,&! partial pressure of CO2 at observational height [pa]
        forc_po2m   ,&! partial pressure of O2 at observational height [pa]
        forc_us     ,&! wind speed in eastward direction [m/s]
        forc_vs     ,&! wind speed in northward direction [m/s]
        forc_t      ,&! temperature at agcm reference height [kelvin]
        forc_q      ,&! specific humidity at agcm reference height [kg/kg]
        forc_prc    ,&! convective precipitation [mm/s]
        forc_prl    ,&! large scale precipitation [mm/s]
        forc_psrf   ,&! atmosphere pressure at the surface [pa]
        forc_pbot   ,&! atmosphere pressure at the bottom of the atmos. model level [pa]
        forc_sols   ,&! atm vis direct beam solar rad onto srf [W/m2]
        forc_soll   ,&! atm nir direct beam solar rad onto srf [W/m2]
        forc_solsd  ,&! atm vis diffuse solar rad onto srf [W/m2]
        forc_solld  ,&! atm nir diffuse solar rad onto srf [W/m2]
        forc_frl    ,&! atmospheric infrared (longwave) radiation [W/m2]
        forc_hgt_u  ,&! observational height of wind [m]
        forc_hgt_t  ,&! observational height of temperature [m]
        forc_hgt_q  ,&! observational height of humidity [m]
        forc_rhoair   ! density air [kg/m3]

! Variables required for restart run
! ----------------------------------------------------------------------
  INTEGER, intent(in) :: &
        idate(3)      ! next time-step /year/julian day/second in a day/

  REAL(r8), intent(inout) :: oro  ! ocean(0)/seaice(2)/ flag
  REAL(r8), intent(inout) :: &
        z_soisno(maxsnl+1:nl_soil)    ,&! layer depth (m)
        dz_soisno(maxsnl+1:nl_soil)   ,&! layer thickness (m)
        t_soisno(maxsnl+1:nl_soil)    ,&! soil + snow layer temperature [K]
        wliq_soisno(maxsnl+1:nl_soil) ,&! liquid water (kg/m2)
        wice_soisno(maxsnl+1:nl_soil) ,&! ice lens (kg/m2)

        t_lake(nl_lake)       ,&! lake temperature (kelvin)
        lake_icefrac(nl_lake) ,&! lake mass fraction of lake layer that is frozen

        t_grnd      ,&! ground surface temperature [k]
        tleaf       ,&! leaf temperature [K]
        ldew        ,&! depth of water on foliage [kg/m2/s]
        sag         ,&! non dimensional snow age [-]
        scv         ,&! snow mass (kg/m2)
        snowdp      ,&! snow depth (m)
        zwt         ,&! the depth to water table [m]
        wa          ,&! water storage in aquifer [mm]

        fveg        ,&! fraction of vegetation cover
        fsno        ,&! fractional snow cover
        sigf        ,&! fraction of veg cover, excluding snow-covered veg [-]
        green       ,&! greenness
        lai         ,&! leaf area index
        sai         ,&! stem area index
 
        coszen      ,&! cosine of solar zenith angle
        alb(2,2)    ,&! averaged albedo [-]
        ssun(2,2)   ,&! sunlit canopy absorption for solar radiation
        ssha(2,2)   ,&! shaded canopy absorption for solar radiation
        thermk      ,&! canopy gap fraction for tir radiation
        extkb       ,&! (k, g(mu)/mu) direct solar extinction coefficient
        extkd         ! diffuse and scattered diffuse PAR extinction coefficient
        
     
! additional diagnostic variables for output
  REAL(r8), intent(out) :: &
        laisun      ,&! sunlit leaf area index
        laisha      ,&! shaded leaf area index
        rstfac      ,&! factor of soil water stress 
        wat         ,&! total water storage
        h2osoi(nl_soil)! volumetric soil water in layers [m3/m3]

! Fluxes
! ----------------------------------------------------------------------
  REAL(r8), intent(out) :: &
        taux        ,&! wind stress: E-W [kg/m/s**2]
        tauy        ,&! wind stress: N-S [kg/m/s**2]
        fsena       ,&! sensible heat from canopy height to atmosphere [W/m2]
        fevpa       ,&! evapotranspiration from canopy height to atmosphere [mm/s]
        lfevpa      ,&! latent heat flux from canopy height to atmosphere [W/2]
        fsenl       ,&! ensible heat from leaves [W/m2]
        fevpl       ,&! evaporation+transpiration from leaves [mm/s]
        etr         ,&! transpiration rate [mm/s]
        fseng       ,&! sensible heat flux from ground [W/m2]
        fevpg       ,&! evaporation heat flux from ground [mm/s]
        olrg        ,&! outgoing long-wave radiation from ground+canopy
        fgrnd       ,&! ground heat flux [W/m2]
        xerr        ,&! water balance error at current time-step [mm/s]
        zerr        ,&! energy balnce errore at current time-step [W/m2]

        tref        ,&! 2 m height air temperature [K]
        qref        ,&! 2 m height air specific humidity
        trad        ,&! radiative temperature [K]
        rsur        ,&! surface runoff (mm h2o/s)
        rnof        ,&! total runoff (mm h2o/s)
        qintr       ,&! interception (mm h2o/s)
        qinfl       ,&! inflitration (mm h2o/s)
        qdrip       ,&! throughfall (mm h2o/s)
        qcharge     ,&! groundwater recharge [mm/s]
       
        rst         ,&! canopy stomatal resistance 
        assim       ,&! canopy assimilation
        respc       ,&! canopy respiration

        sabvsun     ,&! solar absorbed by sunlit vegetation [W/m2]
        sabvsha     ,&! solar absorbed by shaded vegetation [W/m2]
        sabg        ,&! solar absorbed by ground  [W/m2]
        sr          ,&! total reflected solar radiation (W/m2)
        solvd       ,&! incident direct beam vis solar radiation (W/m2)
        solvi       ,&! incident diffuse beam vis solar radiation (W/m2)
        solnd       ,&! incident direct beam nir solar radiation (W/m2)
        solni       ,&! incident diffuse beam nir solar radiation (W/m2)
        srvd        ,&! reflected direct beam vis solar radiation (W/m2)
        srvi        ,&! reflected diffuse beam vis solar radiation (W/m2)
        srnd        ,&! reflected direct beam nir solar radiation (W/m2)
        srni        ,&! reflected diffuse beam nir solar radiation (W/m2)
        solvdln     ,&! incident direct beam vis solar radiation at local noon(W/m2)
        solviln     ,&! incident diffuse beam vis solar radiation at local noon(W/m2)
        solndln     ,&! incident direct beam nir solar radiation at local noon(W/m2)
        solniln     ,&! incident diffuse beam nir solar radiation at local noon(W/m2)
        srvdln      ,&! reflected direct beam vis solar radiation at local noon(W/m2)
        srviln      ,&! reflected diffuse beam vis solar radiation at local noon(W/m2)
        srndln      ,&! reflected direct beam nir solar radiation at local noon(W/m2)
        srniln      ,&! reflected diffuse beam nir solar radiation at local noon(W/m2)

        forc_rain   ,&! rain [mm/s]
        forc_snow   ,&! snow [mm/s]

        emis        ,&! averaged bulk surface emissivity
        z0m         ,&! effective roughness [m]
        zol         ,&! dimensionless height (z/L) used in Monin-Obukhov theory
        rib         ,&! bulk Richardson number in surface layer
        ustar       ,&! u* in similarity theory [m/s]
        qstar       ,&! q* in similarity theory [kg/kg]
        tstar       ,&! t* in similarity theory [K]
        fm          ,&! integral of profile function for momentum
        fh          ,&! integral of profile function for heat
        fq            ! integral of profile function for moisture

! ----------------------- Local  Variables -----------------------------
   REAL(r8) :: &
        calday      ,&! Julian cal day (1.xx to 365.xx)
        endwb       ,&! water mass at the end of time step
        errore      ,&! energy balnce errore (Wm-2)
        errorw      ,&! water balnce errore (mm)
        fiold(maxsnl+1:nl_soil), &! fraction of ice relative to the total water
        w_old       ,&! liquid water mass of the column at the previous time step (mm)
        orb_coszen  ,&! cosine of the solar zenith angle
        sabvg       ,&! solar absorbed by ground + vegetation [W/m2]
        parsun      ,&! PAR by sunlit leaves [W/m2]
        parsha      ,&! PAR by shaded leaves [W/m2]
        qseva       ,&! ground surface evaporation rate (mm h2o/s)
        qsdew       ,&! ground surface dew formation (mm h2o /s) [+]
        qsubl       ,&! sublimation rate from snow pack (mm h2o /s) [+]
        qfros       ,&! surface dew added to snow pack (mm h2o /s) [+]
        rootr(1:nl_soil) ,&! root resistance of a layer, all layers add to 1.0
        scvold      ,&! snow cover for previous time step [mm]
        sm          ,&! rate of snowmelt [kg/(m2 s)]
        ssw         ,&! water volumetric content of soil surface layer [m3/m3]
        tssub(7)    ,&! surface/sub-surface temperatures [K]
        tssea       ,&! sea surface temperature [K]
        totwb       ,&! water mass at the begining of time step
        wt          ,&! fraction of vegetation buried (covered) by snow [-]
        zi_soisno(maxsnl:nl_soil) ! interface level below a "z" level (m)

   REAL(r8) :: &
        prc_rain    ,&! convective rainfall [kg/(m2 s)]
        prc_snow    ,&! convective snowfall [kg/(m2 s)]
        prl_rain    ,&! large scale rainfall [kg/(m2 s)]
        prl_snow    ,&! large scale snowfall [kg/(m2 s)]
        t_precip    ,&! snowfall/rainfall temperature [kelvin]
        bifall      ,&! bulk density of newly fallen dry snow [kg/m3]
        pg_rain     ,&! rainfall onto ground including canopy runoff [kg/(m2 s)]
        pg_snow       ! snowfall onto ground including canopy runoff [kg/(m2 s)]

  INTEGER snl       ,&! number of snow layers
        imelt(maxsnl+1:nl_soil), &! flag for: melting=1, freezing=2, Nothing happended=0
        lb          ,&! lower bound of arrays
        j             ! do looping index

     REAL(r8) :: a, aa 
     INTEGER ps, pe, pc

!======================================================================
!  [1] Solar absorbed by vegetation and ground
!      and precipitation information (rain/snow fall and precip temperature
!======================================================================

      CALL netsolar (ipatch,idate,deltim,patchlonr,patchtype,&
                     forc_sols,forc_soll,forc_solsd,forc_solld,&
                     alb,ssun,ssha,lai,sai,rho,tau,&
                     parsun,parsha,sabvsun,sabvsha,sabg,sabvg,sr,&
                     solvd,solvi,solnd,solni,srvd,srvi,srnd,srni,&
                     solvdln,solviln,solndln,solniln,srvdln,srviln,srndln,srniln)

      CALL rain_snow_temp (forc_t,forc_q,forc_psrf,forc_prc,forc_prl,tcrit,&
                           prc_rain,prc_snow,prl_rain,prl_snow,t_precip,bifall)

      forc_rain = prc_rain + prl_rain
      forc_snow = prc_snow + prl_snow

!======================================================================

                         !         / SOIL GROUND          (patchtype = 0)
IF (patchtype <= 2) THEN ! <=== is - URBAN and BUILT-UP   (patchtype = 1)
                         !         \ WETLAND              (patchtype = 2)

! NOTE: PFT and PC are only for soil patches, i.e., patchtype=0.
!======================================================================
                         !initial set
      scvold = scv       !snow mass at previous time step

      snl = 0
      DO j=maxsnl+1,0
         IF(wliq_soisno(j)+wice_soisno(j)>0.) snl=snl-1
      ENDDO

      zi_soisno(0)=0.
      IF (snl < 0) THEN
      DO j = -1, snl, -1
         zi_soisno(j)=zi_soisno(j+1)-dz_soisno(j+1)
      ENDDO
      ENDIF
      DO j = 1,nl_soil
         zi_soisno(j)=zi_soisno(j-1)+dz_soisno(j)
      ENDDO

      totwb = ldew + scv + sum(wice_soisno(1:)+wliq_soisno(1:)) + wa
      fiold(:) = 0.0
      IF (snl <0 ) THEN
         fiold(snl+1:0)=wice_soisno(snl+1:0)/(wliq_soisno(snl+1:0)+wice_soisno(snl+1:0))
      ENDIF

!----------------------------------------------------------------------
! [2] Canopy interception and precipitation onto ground surface
!----------------------------------------------------------------------

IF (patchtype == 0) THEN

#if(defined USGS_CLASSIFICATION || defined IGBP_CLASSIFICATION)
      CALL LEAF_interception (deltim,dewmx,chil,sigf,lai,sai,tleaf,&
                              prc_rain,prc_snow,prl_rain,prl_snow,&
                              ldew,pg_rain,pg_snow,qintr)
#endif

#ifdef PFT_CLASSIFICATION
      CALL LEAF_interception_pftwrap (ipatch,deltim,dewmx,&
                              prc_rain,prc_snow,prl_rain,prl_snow,&
                              ldew,pg_rain,pg_snow,qintr)
#endif

#ifdef PC_CLASSIFICATION
      CALL LEAF_interception_pcwrap (ipatch,deltim,dewmx,&
                              prc_rain,prc_snow,prl_rain,prl_snow,&
                              ldew,pg_rain,pg_snow,qintr)
#endif

ELSE
      CALL LEAF_interception (deltim,dewmx,chil,sigf,lai,sai,tleaf,&
                              prc_rain,prc_snow,prl_rain,prl_snow,&
                              ldew,pg_rain,pg_snow,qintr)
ENDIF
      
      qdrip = pg_rain + pg_snow

!----------------------------------------------------------------------
! [3] Initilize new snow nodes for snowfall / sleet
!----------------------------------------------------------------------

      CALL newsnow (patchtype,maxsnl,deltim,t_grnd,pg_rain,pg_snow,bifall,&
                    t_precip,zi_soisno(:0),z_soisno(:0),dz_soisno(:0),t_soisno(:0),&
                    wliq_soisno(:0),wice_soisno(:0),fiold(:0),snl,sag,scv,snowdp,fsno)

!----------------------------------------------------------------------
! [4] Energy and Water balance 
!----------------------------------------------------------------------
      lb  = snl + 1           !lower bound of array 

      CALL THERMAL (ipatch   ,patchtype         ,lb                ,deltim            ,&
           trsmx0            ,zlnd              ,zsno              ,csoilc            ,&
           dewmx             ,capr              ,cnfac             ,csol              ,&
           porsl             ,psi0              ,bsw               ,dkdry             ,&
           dksatu            ,lai               ,laisun            ,laisha            ,&
           sai               ,htop              ,hbot              ,sqrtdi            ,&
           rootfr            ,rstfac            ,effcon            ,vmax25            ,&
           slti              ,hlti              ,shti              ,hhti              ,&
           trda              ,trdm              ,trop              ,gradm             ,&
           binter            ,extkn             ,forc_hgt_u        ,forc_hgt_t        ,&
           forc_hgt_q        ,forc_us           ,forc_vs           ,forc_t            ,&
           forc_q            ,forc_rhoair       ,forc_psrf         ,forc_pco2m        ,&
           forc_po2m         ,coszen            ,parsun            ,parsha            ,&
           sabvsun           ,sabvsha           ,sabg              ,forc_frl          ,&
           extkb             ,extkd             ,thermk            ,fsno              ,&
           sigf              ,dz_soisno(lb:)    ,z_soisno(lb:)     ,zi_soisno(lb-1:)  ,&
           tleaf             ,t_soisno(lb:)     ,wice_soisno(lb:)  ,wliq_soisno(lb:)  ,&
           ldew              ,scv               ,snowdp            ,imelt(lb:)        ,&
           taux              ,tauy              ,fsena             ,fevpa             ,&
           lfevpa            ,fsenl             ,fevpl             ,etr               ,&
           fseng             ,fevpg             ,olrg              ,fgrnd             ,&
           rootr             ,qseva             ,qsdew             ,qsubl             ,&
           qfros             ,sm                ,tref              ,qref              ,&
           trad              ,rst               ,assim             ,respc             ,&
           errore            ,emis              ,z0m               ,zol               ,&
           rib               ,ustar             ,qstar             ,tstar             ,&
           fm                ,fh                ,fq                                    )

      CALL WATER (ipatch     ,patchtype         ,lb                ,nl_soil           ,&
           deltim            ,z_soisno(lb:)     ,dz_soisno(lb:)    ,zi_soisno(lb-1:)  ,&
           bsw               ,porsl             ,psi0              ,hksati            ,&
           rootr             ,t_soisno(lb:)     ,wliq_soisno(lb:)  ,wice_soisno(lb:)  ,&
           pg_rain           ,sm                ,etr               ,qseva             ,&
           qsdew             ,qsubl             ,qfros             ,rsur              ,&
           rnof              ,qinfl             ,wtfact            ,pondmx            ,&
           ssi               ,wimp              ,smpmin            ,zwt               ,&
           wa                ,qcharge                                                  )

      IF (snl < 0) THEN
         ! Compaction rate for snow 
         ! Natural compaction and metamorphosis. The compaction rate
         ! is recalculated for every new timestep
         lb  = snl + 1   !lower bound of array 
         CALL snowcompaction (lb,deltim,&
                         imelt(lb:0),fiold(lb:0),t_soisno(lb:0),&
                         wliq_soisno(lb:0),wice_soisno(lb:0),dz_soisno(lb:0))

         ! Combine thin snow elements
         lb = maxsnl + 1
         CALL snowlayerscombine (lb,snl,&
                         z_soisno(lb:1),dz_soisno(lb:1),zi_soisno(lb-1:1),&
                         wliq_soisno(lb:1),wice_soisno(lb:1),t_soisno(lb:1),scv,snowdp)

         ! Divide thick snow elements
         IF(snl<0) &
         CALL snowlayersdivide (lb,snl,&
                         z_soisno(lb:0),dz_soisno(lb:0),zi_soisno(lb-1:0),&
                         wliq_soisno(lb:0),wice_soisno(lb:0),t_soisno(lb:0))
      ENDIF
      
      ! Set zero to the empty node
      IF (snl > maxsnl) THEN
         wice_soisno(maxsnl+1:snl) = 0.
         wliq_soisno(maxsnl+1:snl) = 0.
         t_soisno   (maxsnl+1:snl) = 0.
         z_soisno   (maxsnl+1:snl) = 0.
         dz_soisno  (maxsnl+1:snl) = 0.
      ENDIF

      lb = snl + 1
      !TODO: double check
      t_grnd = t_soisno(lb)
      !t_grnd = fsno*t_soisno(lb) + (1.0-fsno)*t_soisno(1)

      ! ----------------------------------------
      ! energy balance
      ! ----------------------------------------
      zerr=errore
#if(defined CLMDEBUG)
      IF (abs(errore) > .5) THEN
         write(6,*) 'Warning: energy balance violation ',errore,patchclass
      ENDIF
#endif

      ! ----------------------------------------
      ! water balance
      ! ----------------------------------------
      endwb=sum(wice_soisno(1:)+wliq_soisno(1:))+ldew+scv + wa
      errorw=(endwb-totwb)-(forc_prc+forc_prl-fevpa-rnof)*deltim
      IF(patchtype==2) errorw=0.    !wetland
      xerr=errorw/deltim

#if(defined CLMDEBUG)
      IF (abs(errorw) > 1.e-3) THEN
         write(6,*) 'Warning: water balance violation', errorw,patchclass
         !stop
      ENDIF
#endif

!======================================================================

ELSE IF(patchtype == 3)THEN   ! <=== is LAND ICE (glacier/ice sheet) (patchtype = 3) 

!======================================================================
                            !initial set
      scvold = scv          !snow mass at previous time step

      snl = 0
      DO j=maxsnl+1,0
         IF(wliq_soisno(j)+wice_soisno(j)>0.) snl=snl-1
      ENDDO

      zi_soisno(0)=0.
      IF (snl < 0) THEN
      DO j = -1, snl, -1
         zi_soisno(j)=zi_soisno(j+1)-dz_soisno(j+1)
      ENDDO
      ENDIF
      DO j = 1,nl_soil
         zi_soisno(j)=zi_soisno(j-1)+dz_soisno(j)
      ENDDO

      totwb = scv + sum(wice_soisno(1:)+wliq_soisno(1:))
      fiold(:) = 0.0
      IF (snl <0 ) THEN
         fiold(snl+1:0)=wice_soisno(snl+1:0)/(wliq_soisno(snl+1:0)+wice_soisno(snl+1:0))
      ENDIF

      pg_rain = prc_rain + prl_rain
      pg_snow = prc_snow + prl_snow

      !----------------------------------------------------------------
      ! Initilize new snow nodes for snowfall / sleet
      !----------------------------------------------------------------

      CALL newsnow (patchtype,maxsnl,deltim,t_grnd,pg_rain,pg_snow,bifall,&
                    t_precip,zi_soisno(:0),z_soisno(:0),dz_soisno(:0),t_soisno(:0),&
                    wliq_soisno(:0),wice_soisno(:0),fiold(:0),snl,sag,scv,snowdp,fsno)

      !----------------------------------------------------------------
      ! Energy and Water balance 
      !----------------------------------------------------------------
      lb  = snl + 1            !lower bound of array 

      CALL GLACIER_TEMP (lb    ,nl_soil     ,deltim     ,&
                   zlnd        ,zsno        ,capr       ,cnfac       ,&
                   forc_hgt_u  ,forc_hgt_t  ,forc_hgt_q ,forc_us     ,&
                   forc_vs     ,forc_t      ,forc_q     ,forc_rhoair ,&
                   forc_psrf   ,coszen      ,sabg       ,forc_frl    ,&
                   fsno,dz_soisno(lb:),z_soisno(lb:),zi_soisno(lb-1:),&
                   t_soisno(lb:),wice_soisno(lb:),wliq_soisno(lb:)   ,&
                   scv         ,snowdp      ,imelt(lb:) ,taux        ,&
                   tauy        ,fsena       ,fevpa      ,lfevpa      ,&
                   fseng       ,fevpg       ,olrg       ,fgrnd       ,&
                   qseva       ,qsdew       ,qsubl      ,qfros       ,&
                   sm          ,tref        ,qref       ,trad        ,&
                   errore      ,emis        ,z0m        ,zol         ,&
                   rib         ,ustar       ,qstar      ,tstar       ,&
                   fm          ,fh          ,fq                       )


      CALL GLACIER_WATER (nl_soil,maxsnl,deltim                      ,&
                   z_soisno    ,dz_soisno   ,zi_soisno  ,t_soisno    ,&
                   wliq_soisno ,wice_soisno ,pg_rain    ,pg_snow     ,&
                   sm          ,scv         ,snowdp     ,imelt       ,&
                   fiold       ,snl         ,qseva      ,qsdew       ,&
                   qsubl       ,qfros       ,rsur       ,rnof        ,&
                   ssi         ,wimp                                  )


      lb = snl + 1
      t_grnd = t_soisno(lb)

      ! ----------------------------------------
      ! energy and water balance check
      ! ----------------------------------------
      zerr=errore

      endwb=scv+sum(wice_soisno(1:)+wliq_soisno(1:))
      errorw=(endwb-totwb)-(pg_rain+pg_snow-fevpa-rnof)*deltim
      xerr=errorw/deltim

!======================================================================

ELSE IF(patchtype == 4) THEN   ! <=== is LAND WATER BODIES (lake, reservior and river) (patchtype = 4) 

!======================================================================

      snl = 0
      DO j = maxsnl+1, 0
         IF (wliq_soisno(j)+wice_soisno(j) > 0.) THEN
            snl=snl-1
         ENDIF
      ENDDO

      zi_soisno(0) = 0.
      IF (snl < 0) THEN
         DO j = -1, snl, -1
            zi_soisno(j)=zi_soisno(j+1)-dz_soisno(j+1)
         ENDDO
      ENDIF

      DO j = 1,nl_soil
         zi_soisno(j)=zi_soisno(j-1)+dz_soisno(j)
      ENDDO

      scvold = scv          !snow mass at previous time step
      fiold(:) = 0.0
      IF (snl < 0) THEN
         fiold(snl+1:0)=wice_soisno(snl+1:0)/(wliq_soisno(snl+1:0)+wice_soisno(snl+1:0))
      ENDIF

      w_old = sum(wliq_soisno(snl+1:))

      pg_rain = prc_rain + prl_rain
      pg_snow = prc_snow + prl_snow

      CALL newsnow_lake ( &
           ! "in" arguments
           ! ---------------
           maxsnl       ,nl_lake      ,deltim          ,dz_lake         ,&
           pg_rain      ,pg_snow      ,t_precip        ,bifall          ,&

           ! "inout" arguments
           ! ------------------
           t_lake       ,zi_soisno(:0),z_soisno(:0)    ,&
           dz_soisno(:0),t_soisno(:0) ,wliq_soisno(:0) ,wice_soisno(:0) ,&
           fiold(:0)    ,snl          ,sag             ,scv             ,&
           snowdp       ,lake_icefrac )

 
      CALL laketem ( &
           ! "in" laketem arguments
           ! ---------------------------
           patchtype    ,maxsnl       ,nl_soil         ,nl_lake         ,&
           patchlatr    ,deltim       ,forc_hgt_u      ,forc_hgt_t      ,&
           forc_hgt_q   ,forc_us      ,forc_vs         ,forc_t          ,&
           forc_q       ,forc_rhoair  ,forc_psrf       ,forc_sols       ,&
           forc_soll    ,forc_solsd   ,forc_solld      ,sabg            ,&
           forc_frl     ,dz_soisno    ,z_soisno        ,zi_soisno       ,&
           dz_lake      ,lakedepth    ,csol            ,porsl           ,&
           dkdry        ,dksatu       ,&

           ! "inout" laketem arguments
           ! ---------------------------
           t_grnd       ,scv          ,snowdp          ,t_soisno        ,&
           wliq_soisno  ,wice_soisno  ,imelt           ,t_lake          ,&
           lake_icefrac ,&

           ! "out" laketem arguments
           ! ---------------------------
           taux         ,tauy         ,fsena                            ,&
           fevpa        ,lfevpa       ,fseng           ,fevpg           ,&
           qseva        ,qsubl        ,qsdew           ,qfros           ,&
           olrg         ,fgrnd        ,tref            ,qref            ,&
           trad         ,emis         ,z0m             ,zol             ,&
           rib          ,ustar        ,qstar           ,tstar           ,&
           fm           ,fh           ,fq              ,sm )

      CALL snowwater_lake ( &
           ! "in" snowater_lake arguments
           ! ---------------------------
           maxsnl       ,nl_soil      ,nl_lake         ,deltim          ,&
           ssi          ,wimp         ,porsl           ,pg_rain         ,&
           pg_snow      ,dz_lake      ,imelt(:0)       ,fiold(:0)       ,&
           qseva        ,qsubl        ,qsdew           ,qfros           ,&

           ! "inout" snowater_lake arguments
           ! ---------------------------
           z_soisno     ,dz_soisno    ,zi_soisno       ,t_soisno        ,&
           wice_soisno  ,wliq_soisno  ,t_lake          ,lake_icefrac    ,&
           fseng        ,fgrnd        ,snl             ,scv             ,&
           snowdp       ,sm )

      ! We assume the land water bodies have zero extra liquid water capacity 
      ! (i.e.,constant capacity), all excess liquid water are put into the runoff,
      ! this unreasonable assumption should be updated in the future version
      a = (sum(wliq_soisno(snl+1:))-w_old)/deltim
      aa = qseva-(qsubl-qsdew)
      rsur = max(0., pg_rain - aa - a)
      rnof = rsur

      ! Set zero to the empty node
      IF (snl > maxsnl) THEN
         wice_soisno(maxsnl+1:snl) = 0.
         wliq_soisno(maxsnl+1:snl) = 0.
         t_soisno   (maxsnl+1:snl) = 0.
         z_soisno   (maxsnl+1:snl) = 0.
         dz_soisno  (maxsnl+1:snl) = 0.
      ENDIF


!======================================================================

ELSE                     ! <=== is OCEAN (patchtype >= 99) 

!======================================================================
! simple ocean-sea ice model

    tssea = t_grnd
    tssub (1:7) = t_soisno (1:7) 
    CALL SOCEAN (dosst,deltim,oro,forc_hgt_u,forc_hgt_t,forc_hgt_q,&
                 forc_us,forc_vs,forc_t,forc_t,forc_rhoair,forc_psrf,&
                 sabg,forc_frl,tssea,tssub(1:7),scv,&
                 taux,tauy,fsena,fevpa,lfevpa,fseng,fevpg,tref,qref,&
                 z0m,zol,rib,ustar,qstar,tstar,fm,fh,fq,emis,olrg)
                 
               ! null data for sea component
                 z_soisno   (:) = 0.0
                 dz_soisno  (:) = 0.0
                 t_soisno   (:) = 0.0
                 t_soisno (1:7) = tssub(1:7)
                 wliq_soisno(:) = 0.0
                 wice_soisno(:) = 0.0
                 t_grnd  = tssea
                 snowdp  = scv/1000.*20. 

                 trad    = tssea
                 fgrnd   = 0.0
                 rsur    = 0.0
                 rnof    = 0.0

!======================================================================

ENDIF


!======================================================================
! Preparation for the next time step
! 1) time-varying parameters for vegatation
! 2) fraction of snow cover 
! 3) solar zenith angle and
! 4) albedos 
!======================================================================

    ! cosine of solar zenith angle 
    calday = calendarday(idate, gridlond(1))
    coszen = orb_coszen(calday,patchlonr,patchlatr)

    IF (patchtype <= 5) THEN   !LAND
#if(defined DYN_PHENOLOGY)
       ! need to update lai and sai, fveg, green, they are done once in a day only
       IF (dolai) THEN
          CALL LAI_empirical(patchclass,nl_soil,rootfr,t_soisno(1:),lai,sai,fveg,green)
       ENDIF
#endif

! ONLY for soil patches
!NOTE: lai from remove sensing has already considered snow coverage
IF (patchtype == 0) THEN

#if(defined USGS_CLASSIFICATION || defined IGBP_CLASSIFICATION)
       CALL snowfraction (tlai(ipatch),tsai(ipatch),z0m,zlnd,scv,snowdp,wt,sigf,fsno)
       lai = tlai(ipatch)
       sai = tsai(ipatch) * sigf
#endif

#ifdef PFT_CLASSIFICATION
       ps = patch_pft_s(ipatch)      
       pe = patch_pft_e(ipatch)
       CALL snowfraction_pftwrap (ipatch,zlnd,scv,snowdp,wt,sigf,fsno)
       lai_p(ps:pe) = tlai_p(ps:pe)
       sai_p(ps:pe) = tsai_p(ps:pe) * sigf_p(ps:pe)
       lai = tlai(ipatch)
       sai = sum(sai_p(ps:pe)*pftfrac(ps:pe))
#endif

#ifdef PC_CLASSIFICATION
       pc = patch2pc(ipatch)
       CALL snowfraction_pcwrap (ipatch,zlnd,scv,snowdp,wt,sigf,fsno)
       lai_c(:,pc) = tlai_c(:,pc)
       sai_c(:,pc) = tsai_c(:,pc) * sigf_c(:,pc)
       lai = tlai(ipatch)
       sai = sum(sai_c(:,pc)*pcfrac(:,pc))
#endif

ELSE
       CALL snowfraction (tlai(ipatch),tsai(ipatch),z0m,zlnd,scv,snowdp,wt,sigf,fsno)
       lai = tlai(ipatch)
       sai = tsai(ipatch) * sigf
ENDIF

       ! update the snow age 
       IF (snl == 0) sag=0.
       CALL snowage (deltim,t_grnd,scv,scvold,sag)

       ! water volumetric content of soil surface layer [m3/m3]
       ssw = min(1.,1.e-3*wliq_soisno(1)/dz_soisno(1))
       IF (patchtype >= 3) ssw = 1.0

       ! albedos 
       ! we supposed CALL it every time-step, because 
       ! other vegeation related parameters are needed to create
       IF (doalb) THEN
            CALL albland (ipatch, patchtype,&
                 soil_s_v_alb,soil_d_v_alb,soil_s_n_alb,soil_d_n_alb,&
                 chil,rho,tau,fveg,green,lai,sai,coszen,wt,fsno,scv,sag,ssw,t_grnd,&
                 alb,ssun,ssha,thermk,extkb,extkd)
       ENDIF
    ELSE                   !OCEAN
       sag = 0.0
       IF(doalb)THEN
            CALL albocean (oro,scv,coszen,alb)
       ENDIF
    ENDIF

    ! zero-filling set for glacier/ice-sheet/land water bodies/ocean components
    IF (patchtype > 2) THEN
       lai = 0.0
       sai = 0.0
       laisun = 0.0
       laisha = 0.0
       green = 0.0
       fveg = 0.0
       sigf = 0.0

       ssun(:,:) = 0.0
       ssha(:,:) = 0.0
       thermk = 0.0
       extkb = 0.0
       extkd = 0.0

       tleaf = forc_t
       ldew  = 0.0
       fsenl = 0.0
       fevpl = 0.0
       etr   = 0.0
       assim = 0.0
       respc = 0.0

       zerr=0.
       xerr=0.

       qinfl = 0.
       qdrip = forc_rain + forc_snow
       qintr = 0.
       h2osoi = 0. 
       rstfac = 0.
       zwt = 0.
       wa = 4800.
       qcharge = 0.
    ENDIF
      
    h2osoi = wliq_soisno(1:)/(dz_soisno(1:)*denh2o) + wice_soisno(1:)/(dz_soisno(1:)*denice)
    wat = sum(wice_soisno(1:)+wliq_soisno(1:))+ldew+scv + wa
!----------------------------------------------------------------------

END SUBROUTINE CLMMAIN
