#include <define.h>

MODULE LEAF_temperature_PC

!-----------------------------------------------------------------------
 USE precision
 IMPLICIT NONE
 SAVE
 
! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: LeafTempPC

! PRIVATE MEMBER FUNCTIONS:
  PRIVATE :: dewfraction
  PRIVATE :: cal_z0_displa


!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE  LeafTempPC (ipatch,npft,deltim,csoilc,dewmx ,htvp    ,&
              fcover  ,canlev  ,htop    ,hbot    ,lai     ,sai     ,&
              sqrtdi  ,effcon  ,vmax25  ,slti    ,hlti    ,shti    ,&
              hhti    ,trda    ,trdm    ,trop    ,gradm   ,binter  ,&
              extkn   ,extkb   ,extkd   ,hu      ,ht      ,hq      ,&
              us      ,vs      ,thm     ,th      ,thv     ,qm      ,&
              psrf    ,rhoair  ,parsun  ,parsha  ,fsun    ,sabv    ,&
              frl     ,thermk  ,fshade  ,rstfac  ,po2m    ,pco2m   ,&
              z0h_g   ,obug    ,ustarg  ,zlnd    ,zsno    ,fsno    ,&                               
              sigf    ,etrc    ,tg      ,qg      ,dqgdT   ,emg     ,&
              z0mpc   ,tl      ,ldew    ,taux    ,tauy    ,fseng   ,&
              fevpg   ,cgrnd   ,cgrndl  ,cgrnds  ,tref    ,qref    ,&
              rst     ,assim   ,respc   ,fsenl   ,fevpl   ,etr     ,&
              dlrad   ,ulrad   ,z0m     ,zol     ,rib     ,ustar   ,&
              qstar   ,tstar   ,fm      ,fh      ,fq                ) 
 
!=======================================================================
! Original author : Yongjiu Dai, August 15, 2001
!
! Foliage energy conservation is given by foliage energy budget equation
!                      Rnet - Hf - LEf = 0
! The equation is solved by Newton-Raphson iteration, in which this iteration
! includes the calculation of the photosynthesis and stomatal resistance, and the
! integration of turbulent flux profiles. The sensible and latent heat
! transfer between foliage and atmosphere and ground is linked by the equations:
!                      Ha = Hf + Hg and Ea = Ef + Eg
!
! ________________
! REVISION HISTORY:
! 07/09/2014, Hua Yuan: imbalanced energy due to T/q adjustment is
!                       allocated to sensible heat flux.
!
!=======================================================================

  USE precision
  USE PhysicalConstants, only: vonkar, grav, hvap, cpair, stefnc
  USE FRICTION_VELOCITY
  USE ASSIM_STOMATA_conductance
  IMPLICIT NONE
 
!-----------------------Arguments---------------------------------------

  INTEGER,  intent(in) :: ipatch
  INTEGER , intent(in) :: &
        npft,       &! potential PFT number in a column
        canlev(npft) ! potential canopy layer in a column

  REAL(r8), intent(in) :: &
        deltim,     &! seconds in a time step [second]
        csoilc,     &! drag coefficient for soil under canopy [-]
        dewmx,      &! maximum dew
        htvp         ! latent heat of evaporation (/sublimation) [J/kg]

! vegetation parameters
  REAL(r8), dimension(npft), intent(in) :: &
        fcover,     &! PFT fractiona coverage [-]
        htop,       &! PFT crown top height [m]
        hbot,       &! PFT crown bottom height [m]
        lai,        &! adjusted leaf area index for seasonal variation [-]
        sai,        &! stem area index  [-]
        sqrtdi,     &! inverse sqrt of leaf dimension [m**-0.5]

        effcon,     &! quantum efficiency of RuBP regeneration (mol CO2 / mol quanta)
        vmax25,     &! maximum carboxylation rate at 25 C at canopy top
                     ! the range : 30.e-6 <-> 100.e-6 (mol co2 m-2 s-1)
        shti,       &! slope of high temperature inhibition function     (s1)
        hhti,       &! 1/2 point of high temperature inhibition function (s2)
        slti,       &! slope of low temperature inhibition function      (s3)
        hlti,       &! 1/2 point of low temperature inhibition function  (s4)
        trda,       &! temperature coefficient in gs-a model             (s5)
        trdm,       &! temperature coefficient in gs-a model             (s6)
        trop,       &! temperature coefficient in gs-a model         (273+25)
        gradm,      &! conductance-photosynthesis slope parameter
        binter,     &! conductance-photosynthesis intercept
        extkn        ! coefficient of leaf nitrogen allocation

! input variables
  REAL(r8), intent(in) :: &
        hu,         &! observational height of wind [m]
        ht,         &! observational height of temperature [m]
        hq,         &! observational height of humidity [m]
        us,         &! wind component in eastward direction [m/s]
        vs,         &! wind component in northward direction [m/s]
        thm,        &! intermediate variable (tm+0.0098*ht)
        th,         &! potential temperature (kelvin)
        thv,        &! virtual potential temperature (kelvin)
        qm,         &! specific humidity at reference height [kg/kg]
        psrf,       &! pressure at reference height [pa]
        rhoair,     &! density air [kg/m**3]

        parsun(npft),&! par absorbed per unit sunlit lai [w/m**2]
        parsha(npft),&! par absorbed per unit shaded lai [w/m**2]
        fsun(npft),  &! sunlit fraction of canopy
        sabv(npft),  &! solar radiation absorbed by vegetation [W/m2]
        frl,         &! atmospheric infrared (longwave) radiation [W/m2]

        extkb(npft), &! (k, g(mu)/mu) direct solar extinction coefficient
        extkd(npft), &! diffuse and scattered diffuse PAR extinction coefficient
        thermk(npft),&! canopy gap fraction for tir radiation
        fshade(npft),&! shadow for each PFT
        rstfac(npft),&! factor of soil water stress to plant physiologocal processes

        po2m,       &! atmospheric partial pressure  o2 (pa)
        pco2m,      &! atmospheric partial pressure co2 (pa)

        z0h_g,      &! bare soil roughness length, sensible heat [m]
        obug,       &! bare soil obu
        ustarg,     &! bare soil ustar
        zlnd,       &! roughness length for soil [m]
        zsno,       &! roughness length for snow [m]
        fsno,       &! fraction of snow cover on ground
        
        sigf(npft), &! fraction of veg cover, excluding snow-covered veg [-]
        etrc(npft), &! maximum possible transpiration rate (mm/s)
        tg,         &! ground surface temperature [K]
        qg,         &! specific humidity at ground surface [kg/kg]
        dqgdT,      &! temperature derivative of "qg"
        emg          ! vegetation emissivity

  REAL(r8), dimension(npft), intent(inout) :: &
        tl,         &! leaf temperature [K]
        ldew         ! depth of water on foliage [mm]
  
  REAL(r8), intent(inout) :: &
        dlrad,      &! downward longwave radiation blow the canopy [W/m2]
        ulrad,      &! upward longwave radiation above the canopy [W/m2]
        taux,       &! wind stress: E-W [kg/m/s**2]
        tauy,       &! wind stress: N-S [kg/m/s**2]
        fseng,      &! sensible heat flux from ground [W/m2]
        fevpg,      &! evaporation heat flux from ground [mm/s]
        tref,       &! 2 m height air temperature (kelvin)
        qref         ! 2 m height air specific humidity

  REAL(r8), dimension(npft), intent(out) :: &
        z0mpc,      &! z0m for individual PFT
        rst,        &! stomatal resistance
        assim,      &! rate of assimilation
        respc,      &! rate of respiration
        fsenl,      &! sensible heat from leaves [W/m2]
        fevpl,      &! evaporation+transpiration from leaves [mm/s]
        etr          ! transpiration rate [mm/s]

  REAL(r8), intent(inout) :: &
        z0m,        &! effective roughness [m]
        zol,        &! dimensionless height (z/L) used in Monin-Obukhov theory
        rib,        &! bulk Richardson number in surface layer
        ustar,      &! friction velocity [m/s]
        tstar,      &! temperature scaling parameter
        qstar,      &! moisture scaling parameter
        fm,         &! integral of profile function for momentum
        fh,         &! integral of profile function for heat
        fq           ! integral of profile function for moisture

  REAL(r8), intent(inout) :: &
        cgrnd,      &! deriv. of soil energy flux wrt to soil temp [w/m2/k]
        cgrndl,     &! deriv, of soil latent heat flux wrt soil temp [w/m2/k]
        cgrnds       ! deriv of soil sensible heat flux wrt soil temp [w/m**2/k]

!-----------------------Local Variables---------------------------------
! assign iteration parameters
   INTEGER, parameter :: itmax  = 40   !maximum number of iteration
   INTEGER, parameter :: itmin  = 6    !minimum number of iteration
   REAL(r8),parameter :: delmax = 3.0  !maximum change in leaf temperature [K]
   REAL(r8),parameter :: dtmin  = 0.01 !max limit for temperature convergence [K]
   REAL(r8),parameter :: dlemin = 0.1  !max limit for energy flux convergence [w/m2]

   REAL(r8) dtl(0:itmax+1,npft)     !difference of tl between two iterative step

   REAL(r8) :: &
        zldis,      &! reference height "minus" zero displacement heght [m]
        zii,        &! convective boundary layer height [m]
        z0mv,       &! roughness length, momentum [m]
        z0hv,       &! roughness length, sensible heat [m]
        z0qv,       &! roughness length, latent heat [m]
        zeta,       &! dimensionless height used in Monin-Obukhov theory
        beta,       &! coefficient of conective velocity [-]
        wc,         &! convective velocity [m/s]
        wc2,        &! wc**2
        dth,        &! diff of virtual temp. between ref. height and surface 
        dthv,       &! diff of vir. poten. temp. between ref. height and surface
        dqh,        &! diff of humidity between ref. height and surface
        obu,        &! monin-obukhov length (m)
        um,         &! wind speed including the stablity effect [m/s]
        ur,         &! wind speed at reference height [m/s]
        uaf,        &! velocity of air within foliage [m/s]
        fh2m,       &! relation for temperature at 2m
        fq2m,       &! relation for specific humidity at 2m
        fm10m,      &! integral of profile function for momentum at 10m
        thvstar,    &! virtual potential temperature scaling parameter
        eah,        &! canopy air vapor pressure (pa)
        pco2g,      &! co2 pressure (pa) at ground surface (pa)
        pco2a,      &! canopy air co2 pressure (pa)

        fdry(npft), &! fraction of foliage that is green and dry [-]
        fwet(npft), &! fraction of foliage covered by water [-]
        cf,         &! heat transfer coefficient from leaves [-]
        rb(npft),   &! leaf boundary layer resistance [s/m]
        rbone,      &! canopy bulk boundary layer resistance 
        rbsun,      &! bulk boundary layer resistance of sunlit fraction of canopy
        rbsha,      &! bulk boundary layer resistance of shaded fraction of canopy
        ram,        &! aerodynamical resistance [s/m]
        rah,        &! thermal resistance [s/m]
        raw,        &! moisture resistance [s/m]
        clai,       &! canopy heat capacity [Jm-2K-1]
        cfh(npft),  &! heat conductance for leaf [m/s]
        cfw(npft),  &! latent heat conductance for leaf [m/s]
        wtl0(npft), &! normalized heat conductance for air and leaf [-]
        wtlq0(npft),&! normalized latent heat cond. for air and leaf [-]

        ei(npft),   &! vapor pressure on leaf surface [pa]
        deidT(npft),&! derivative of "ei" on "tl" [pa/K]
        qsatl(npft),&! leaf specific humidity [kg/kg]
        qsatldT(npft),&! derivative of "qsatl" on "tlef"

        del(npft),  &! absolute change in leaf temp in current iteration [K]
        del2(npft), &! change in leaf temperature in previous iteration [K]
        dele(npft), &! change in heat fluxes from leaf [K]
        dele2(npft),&! change in heat fluxes from leaf [K]
        det,        &! maximum leaf temp. change in two consecutive iter [K]
        dee,        &! maximum leaf temp. change in two consecutive iter [K]
 
        obuold,     &! monin-obukhov length from previous iteration
        tlbef(npft),&! leaf temperature from previous iteration [K]
        err,        &! balance error

        fsha(npft),    &! shaded fraction of canopy
        laisun(npft),  &! sunlit leaf area index, one-sided
        laisha(npft),  &! shaded leaf area index, one-sided
        rssun(npft),   &! sunlit leaf stomatal resistance [s/m]
        rssha(npft),   &! shaded leaf stomatal resistance [s/m]
        assimsun(npft),&! sunlit leaf assimilation rate [umol co2 /m**2/ s] [+]
        assimsha(npft),&! shaded leaf assimilation rate [umol co2 /m**2/ s] [+]
        respcsun(npft),&! sunlit leaf respiration rate [umol co2 /m**2/ s] [+]
        respcsha(npft),&! shaded leaf respiration rate [umol co2 /m**2/ s] [+]
        
        rsoil,      &! soil respiration
        gah2o,      &! conductance between canopy and atmosphere
        gdh2o,      &! conductance between canopy and ground
        tprcor,     &! tf*psur*100./1.013e5

        fht,        &! integral of profile function for heat at the top layer
        fqt,        &! integral of profile function for moisture at the top layer
        phih         ! phi(h), similarity function for sensible heat


   INTEGER it, nmozsgn 

   REAL(r8) delta(npft), fac(npft)
   REAL(r8) evplwet(npft), evplwet_dtl(npft), etr_dtl(npft), elwmax, elwdif
   REAL(r8) irab(npft), dirab_dtl(npft), fsenl_dtl(npft), fevpl_dtl(npft)  
   REAL(r8) w, csoilcn, z0mg, z0hg, z0qg, cintsun(3, npft), cintsha(3, npft)
   REAL(r8), dimension(npft) :: fevpl_bef, fevpl_noadj, dtl_noadj, erre


   ! .................................................................
   ! defination for 3d run
   ! .................................................................

   INTEGER , parameter :: nlay = 3
   REAL(r8), parameter :: pi   = 3.14159265358979323846_r8  !pi

   REAL(r8), parameter :: &
        c1   = 0.320,  &! parameter to calculate drag coefficients of Massman's method
        c2   = 0.264,  &! parameter to calculate drag coefficients of Massman's method
        c3   = 15.1,   &! parameter to calculate drag coefficients of Massman's method
        iw   = 0.5,    &! parameter to calculate alpha of Goudriaa's method  
        Cd   = 0.2,    &! leaf drag coefficient
        cd1  = 7.5,    &! a free parameter for d/h calculation, Raupach 1992, 1994
        psih = 0.193    ! psih = ln(cw) - 1 + cw^-1, cw = 2, Raupach 1994

   REAL(r8) :: sqrtdragc! sqrt(drag coefficient)
   REAL(r8) :: lm       ! mix length within canopy
   REAL(r8) :: fai      ! canopy frontal area index

   REAL(r8), dimension(0:nlay) :: &
        z0m_lays,      &! roughness length for momentum for the layer and below
        z0h_lays,      &! roughness length for SH for the layer and below
        z0q_lays,      &! roughness length for LH for the layer and below
        displa_lays,   &! displacement height for the layer and below
        fcover_lays     ! vegetation fractional cover for this layer and above
   
   REAL(r8), dimension(npft) :: &
        lsai            ! lai + sai

   REAL(r8), dimension(nlay) :: &
        htop_lay,      &! canopy crown top for each layer
        hbot_lay,      &! canopy crown bottom for each layer
        fcover_lay,    &! vegetation fractional coverage for each layer
        lsai_lay,      &! (lai+sai) for each layer
        a_lay,         &! exponential extinction factor for u/k decline within canopy 
        a_lay_i63,     &! exponential extinction factor for u/k decline within canopy (Inoue 1963)
        a_lay_k71,     &! exponential extinction factor for u/k decline within canopy (Kondo 1971)
        a_lay_g77,     &! exponential extinction factor for u/k decline within canopy (Groudrian 1977)
        a_lay_m97,     &! exponential extinction factor for u/k decline within canopy (Massman 1997)
        utop_lay,      &! wind speed at layer top [m/s]
        ubot_lay,      &! wind speed at layer bottom [m/s]
        ueff_lay,      &! effective wind speed within canopy layer [m/s]
        ueff_lay_,     &! effective wind speed within canopy layer [m/s]
        ueff_lay_norm, &! normalized effective wind speed within canopy layer [m/s]
        ktop_lay,      &! eddy coefficient at layer top 
        kbot_lay,      &! eddy coefficient at layer bottom
        z0m_lay,       &! roughness length for the vegetation covered area
        displa_lay,    &! displacement height for the vegetaion covered area
        taf,           &! air temperature within canopy space [K]
        qaf,           &! humidity of canopy air [kg/kg]
        rd,            &! aerodynamic resistance between layers [s/m]
        rd_,           &! aerodynamic resistance between layers [s/m]
        cah,           &! heat conductance for air [m/s]
        cgh,           &! heat conductance for ground [m/s]
        caw,           &! latent heat conductance for air [m/s]
        cgw,           &! latent heat conductance for ground [m/s]
        wtshi,         &! sensible heat resistance for air, grd and leaf [-]
        wtsqi,         &! latent heat resistance for air, grd and leaf [-]
        wta0,          &! normalized heat conductance for air [-]
        wtg0,          &! normalized heat conductance for ground [-]
        wtaq0,         &! normalized latent heat conductance for air [-]
        wtgq0,         &! normalized heat conductance for ground [-]
        wtll,          &! sum of normalized heat conductance for air and leaf 
        wtlql           ! sum of normalized heat conductance for air and leaf 
 
   REAL(r8) :: ktop, utop, fmtop, bee, tmpw1, tmpw2, fact, facq 
   
   INTEGER i, clev
   INTEGER toplay, botlay, upplay, numlay
   INTEGER d_opt, rb_opt, rd_opt

   REAL(r8) :: displa 

   ! variables for longwave transfer calculation
   ! .................................................................
   REAL(r8) :: tdn(0:4,0:4)     !downward transfer coefficient matrix for LW
   REAL(r8) :: tup(0:4,0:4)     !upward transfer coefficient matrix for LW
   REAL(r8) :: thermk_lay(nlay) !transmittance of longwave radiation for each layer
   REAL(r8) :: fshade_lay(nlay) !shadow of each layer
   REAL(r8) :: L(nlay)          !longwave radiation emitted by canopy layer
   REAL(r8) :: Ltd(nlay)        !trasmitted downward longwave radiation from canopy layer
   REAL(r8) :: Ltu(nlay)        !trasmitted upward longwave radiation from canopy layer
   REAL(r8) :: Lin(0:4)         !incomming longwave radiation for each layer
   REAL(r8) :: Ld(0:4)          !total downward longwave radiation for each layer
   REAL(r8) :: Lu(0:4)          !total upward longwave radiation for each layer
   REAL(r8) :: Lg               !emitted longwave radiation from ground
   REAL(r8) :: Lv(npft)         !absorbed longwave raidation for each pft
   REAL(r8) :: dLv(npft)        !LW change due to temperature change
   REAL(r8) :: dLvpar(nlay)     !temporal variable for calcualting dLv

!-----------------------End Variable List-------------------------------

! initialization of errors and  iteration parameters
       it       = 1    !counter for leaf temperature iteration
       del(:)   = 0.0  !change in leaf temperature from previous iteration
       dele(:)  = 0.0  !latent head flux from leaf for previous iteration

       dtl(:,:) = 0.
       fevpl_bef(:) = 0.

       d_opt  = 2
       rd_opt = 3
       rb_opt = 3

! initial values for z0hg, z0qg
       z0mg = (1.-fsno)*zlnd + fsno*zsno
       z0hg = z0mg
       z0qg = z0mg
 
!-----------------------------------------------------------------------
! scaling-up coefficients from leaf to canopy
!-----------------------------------------------------------------------

! note: need to sperate to sunlit/shaded pars
!-----------------------------------------------------------------------

! partion visible canopy absorption to sunlit and shaded fractions
! to get average absorbed par for sunlit and shaded leaves
       fsha(:)   = 1. - fsun(:)
       laisun(:) = lai(:)*fsun(:)
       laisha(:) = lai(:)*fsha(:)

       cintsun(1,:) = (1.-exp(-(0.110+extkb)*lai))/(0.110+extkb)
       cintsun(2,:) = (1.-exp(-(extkb+extkd)*lai))/(extkb+extkd)
       cintsun(3,:) = (1.-exp(-extkb*lai))/extkb

       cintsha(1,:) = (1.-exp(-0.110*lai))/0.110 - cintsun(1,:)
       cintsha(2,:) = (1.-exp(-extkd*lai))/extkd - cintsun(2,:)
       cintsha(3,:) = lai(:) - cintsun(3,:)

!-----------------------------------------------------------------------
! get fraction of wet and dry canopy surface (fwet & fdry)
! initial saturated vapor pressure and humidity and their derivation
!-----------------------------------------------------------------------

       !clai = 4.2 * 1000. * 0.2
       clai = 0.0
       lsai(:) = lai(:) + sai(:)

       DO i = 1, npft
          IF (fcover(i) > 0 .and. lsai(i)>0) THEN
             CALL dewfraction (sigf(i),lai(i),sai(i),dewmx,ldew(i),fwet(i),fdry(i))
             CALL qsadv(tl(i),psrf,ei(i),deiDT(i),qsatl(i),qsatlDT(i))
          ENDIF
       ENDDO

!-----------------------------------------------------------------------
! initial for fluxes profile
!-----------------------------------------------------------------------

       nmozsgn = 0     !number of times moz changes sign
       obuold  = 0.    !monin-obukhov length from previous iteration
       zii     = 1000. !m (pbl height)
       beta    = 1.    !- (in computing W_*)

!-----------------------------------------------------------------------
! calculate layer average propeties: height (htop_lay, hbot_lay), lsai_lay, ...
! !!NOTE: adjustment may needed for htop_lay/hbot_lay
!-----------------------------------------------------------------------
       htop_lay(:)   = 0
       hbot_lay(:)   = 0
       lsai_lay(:)   = 0
       fcover_lay(:) = 0

       DO i = 1, npft
          IF (fcover(i)>0 .and. lsai(i)>0) THEN
             clev = canlev(i)
             htop_lay(clev) = htop_lay(clev) + htop(i) * fcover(i)
             hbot_lay(clev) = hbot_lay(clev) + hbot(i) * fcover(i)
             lsai_lay(clev) = lsai_lay(clev) + lsai(i) * fcover(i)
             fcover_lay(clev) = fcover_lay(clev) + fcover(i)
          ENDIF
       ENDDO

       DO i = 1, nlay
          IF (fcover_lay(i) > 0) THEN
             htop_lay(i) = htop_lay(i) / fcover_lay(i)
             hbot_lay(i) = hbot_lay(i) / fcover_lay(i)
             lsai_lay(i) = lsai_lay(i) / fcover_lay(i)
          ENDIF
       ENDDO

       ! 12/12/2017: calculate fcover_lays
! 03/16/2020, yuan: TODO, determine to set fc=0 or fcover above for
! gaps between layers, 0 maybe more consistent
       fcover_lays(0) = sum(fcover_lay(:))
       fcover_lays(1) = sum(fcover_lay(1:3))
       fcover_lays(2) = sum(fcover_lay(2:3))
       fcover_lays(3) = sum(fcover_lay(3:3))
       fcover_lays(:) = 0.

!-----------------------------------------------------------------------
! scaling factor bee
!-----------------------------------------------------------------------
! 09/26/2017
! NOTE: bee value, the default is 1
       !bee = 1. / sum(fcover_lay(1:3))
       !bee = 2./ (1+sum(fcover_lay(1:3)))
       bee = 1.

!-----------------------------------------------------------------------
! calculate z0m and displa for PFTs
!-----------------------------------------------------------------------
       DO i = 1, npft
          IF (lsai(i) > 1.e-6) THEN
             CALL cal_z0_displa(lsai(i), htop(i), 1., z0mpc(i), displa)
          ELSE 
             z0mpc(i) = z0mg
          ENDIF
       ENDDO 
       
!-----------------------------------------------------------------------
! calculate z0m and displa for layers
!-----------------------------------------------------------------------

       displa_lay (:) = 0.
       displa_lays(:) = 0.
       z0m_lay    (:) = 0. 
       z0m_lays   (:) = 0. 

       DO i = 1, nlay
          IF (fcover_lay(i)>0 .and. lsai_lay(i)>0) THEN
             
             CALL cal_z0_displa(lsai_lay(i), htop_lay(i), 1., z0m_lay(i), displa_lay(i))
             
             CALL cal_z0_displa(lsai_lay(i), htop_lay(i), fcover_lay(i), z0m_lays(i), displa_lays(i))

          ENDIF
       ENDDO

       ! ground
       z0m_lays(0)    = z0mg
       displa_lays(0) = 0.

       ! 10/05/2017: robust check
       WHERE (z0m_lays(:) < z0mg) z0m_lays(:) = z0mg
       WHERE (z0m_lay(:)  < z0mg) z0m_lay(:)  = z0mg

       ! maximum assumption
       z0m_lays(1) = maxval(z0m_lays(0:1))
       z0m_lays(2) = maxval(z0m_lays(0:2))
       z0m_lays(3) = maxval(z0m_lays(0:3))

       displa_lays(1) = maxval(displa_lays(0:1))
       displa_lays(2) = maxval(displa_lays(0:2))
       displa_lays(3) = maxval(displa_lays(0:3))
       
       ! roughness length and displacement height for sensible 
       ! and latent heat transfer
       z0h_lays(:) = z0m_lays(:)
       z0q_lays(:) = z0m_lays(:)
      
!-----------------------------------------------------------------------
! calculate layer a_lay
!-----------------------------------------------------------------------
       ! initialization
       a_lay(:)     = 0.
       a_lay_i63(:) = 0.
       a_lay_k71(:) = 0.
       a_lay_g77(:) = 0.
       a_lay_m97(:) = 0.

       DO i = 1, nlay
          IF (fcover_lay(i)>0 .and. lsai_lay(i)>0) THEN

             ! mixing length and sqrt(drag coefficient)
             lm = vonkar*(htop_lay(i) - displa_lay(i))

             ! Raupach, 1992
             fai   = 1. - exp(-0.5*lsai_lay(i))
             sqrtdragc = min( (0.003+0.3*fai)**0.5, 0.3 )

             ! Inoue, 1963
             a_lay_i63(i) = htop_lay(i) * &
                           (Cd*lsai_lay(i)/(2.*htop_lay(i)*lm**2))**(1./3.)

             ! Kondo, 1971
             a_lay_k71(i) = htop_lay(i)/(htop_lay(i)-displa_lay(i))/ &
                            (vonkar/sqrtdragc)

             ! Goudriaan, 1977
             a_lay_g77(i) = (Cd*lsai_lay(i)*htop_lay(i)/lm)**0.5

             ! Massman, 1997
             a_lay_m97(i) = Cd*lsai_lay(i) / (2.*sqrtdragc**2)
             
             a_lay(i) = a_lay_k71(i)

             displa_lay(i) = max(htop_lay(i)/2., displa_lay(i))

          ENDIF
       ENDDO

!-----------------------------------------------------------------------
! claculate layer info
! how may layers, top layer and bottom layer number
!-----------------------------------------------------------------------

       toplay = 0
       botlay = 0
       numlay = 0

       DO i = nlay, 1, -1
          IF (fcover_lay(i)>0 .and. lsai_lay(i)>0) THEN
             
             ! to count the layer number
             numlay = numlay + 1
             IF (toplay .eq. 0) THEN
                ! set the top layer to current layer
                toplay = i
             ENDIF

             ! set this layer to be the bottom layer 
             botlay = i

             displa_lay(i) = max(displa_lay(i), hbot_lay(i))
          ENDIF
       ENDDO

!-----------------------------------------------------------------------
! calculate transmittance of longwave radiation for each layer
! diffuse case
!-----------------------------------------------------------------------
       
       thermk_lay(:) = 0.
       fshade_lay(:) = 0.

       DO i = 1, npft
          IF (fshade(i) > 0) THEN
             clev = canlev(i)
             thermk_lay(clev) = thermk_lay(clev) + fshade(i) * thermk(i)
             fshade_lay(clev) = fshade_lay(clev) + fshade(i)
          ENDIF
       ENDDO
       
       DO i = 1, nlay
          IF (fshade_lay(i) > 0) THEN
             thermk_lay(i) = thermk_lay(i) / fshade_lay(i)
          ELSE
             thermk_lay(i) = 1.
          ENDIF
       ENDDO

!-----------------------------------------------------------------------
! calculate the transfer matrix for long-wave radiation transfer
! direct case
! NOTE: don't need to calculate at each step
!-----------------------------------------------------------------------

       tdn(:,:) = 0.
       tup(:,:) = 0.

       tdn(1,0) = 1.
       tdn(2,0) = 1 - fshade_lay(1)
       tdn(3,0) = 1 - fshade_lay(1) - fshade_lay(2) + fshade_lay(1)*fshade_lay(2)
       tdn(4,0) = 1 - fshade_lay(1) - fshade_lay(2) - fshade_lay(3) &
                + fshade_lay(1)*fshade_lay(2) &
                + fshade_lay(1)*fshade_lay(3) &
                + fshade_lay(2)*fshade_lay(3) &
                - fshade_lay(1)*fshade_lay(2)*fshade_lay(3)

       tdn(2,1) = fshade_lay(1)
       tdn(3,1) = (1 - fshade_lay(2))*fshade_lay(1)
       tdn(4,1) = (1 - fshade_lay(2) - fshade_lay(3) + fshade_lay(2)*fshade_lay(3))*fshade_lay(1)

       tdn(3,2) = fshade_lay(2)
       tdn(4,2) = (1 - fshade_lay(3))*fshade_lay(2)
       tdn(4,3) = fshade_lay(3)
         
       tup(0,1) = fshade_lay(1)
       tup(0,2) = (1 - fshade_lay(1))*fshade_lay(2)
       tup(1,2) = fshade_lay(2)

       tup(0,3) = (1 - fshade_lay(1) - fshade_lay(2) + fshade_lay(1)*fshade_lay(2))*fshade_lay(3)
       tup(1,3) = (1 - fshade_lay(2))*fshade_lay(3)
       tup(2,3) = fshade_lay(3)

       tup(0,4) = tdn(4,0)
       tup(1,4) = 1 - fshade_lay(2) - fshade_lay(3) + fshade_lay(2)*fshade_lay(3)
       tup(2,4) = 1 - fshade_lay(3)
       tup(3,4) = 1.

!-----------------------------------------------------------------------
! calculate parameters for delta(Lv) for LW radiation transfer
!-----------------------------------------------------------------------
       dLvpar(1) = 1.
       dLvpar(2) = ( (1-fshade_lay(1)) + thermk_lay(1)*fshade_lay(1) )**2
       dLvpar(3) = ( tdn(3,0) + thermk_lay(2)*fshade_lay(2)*(1-fshade_lay(1)+thermk_lay(1)*fshade_lay(1)) &
                 + (1-fshade_lay(2))*thermk_lay(1)*fshade_lay(1) )**2

!-----------------------------------------------------------------------
! first guess for taf and qaf for each layer
! a large differece from previous schemes
!-----------------------------------------------------------------------
       taf(:) = 0.
       qaf(:) = 0.

       ! 05/02/2016: set taf/qaf according to layer number
       IF (numlay .eq. 1) THEN 
          taf(toplay) = 0.5 * (tg + thm)
          qaf(toplay) = 0.5 * (qm + qg )
       ENDIF

       IF (numlay .eq. 2) THEN 
          taf(botlay) = (2.*tg + thm)/3.
          qaf(botlay) = (2.*qg + qm )/3.
          taf(toplay) = (tg + 2.*thm)/3.
          qaf(toplay) = (qg + 2.*qm )/3.
       ENDIF

       IF (numlay .eq. 3) THEN 
          taf(1) = (3.*tg + thm)/4.
          qaf(1) = (3.*qg + qm )/4.
          taf(2) = (tg + thm )/2.
          qaf(2) = (qg + qm  )/2.
          taf(3) = (tg + 3.*thm)/4.
          qaf(3) = (qg + 3.*qm )/4.
       ENDIF

!-----------------------------------------------------------------------
! some environment variables
! how to calculate rsoil and what is its usage?
!-----------------------------------------------------------------------
       pco2a = pco2m
       tprcor = 44.6*273.16*psrf/1.013e5
       rsoil = 0.   !respiration (mol m-2 s-1)
!      rsoil = 1.22e-6*exp(308.56*(1./56.02-1./(tg-227.13)))
!      rsoil = rstfac * 0.23 * 15. * 2.**((tg-273.16-10.)/10.) * 1.e-6
!      rsoil = 5.22 * 1.e-6
       rsoil = 0.22 * 1.e-6

! initialization and input values for Monin-Obukhov
       ! have been set before
       z0mv = z0m_lays(3); z0hv = z0m_lays(3); z0qv = z0m_lays(3)
       ur = max(0.1, sqrt(us*us+vs*vs))    !limit set to 0.1
       dth = thm - taf(toplay)
       dqh = qm - qaf(toplay)
       dthv = dth*(1.+0.61*qm) + 0.61*th*dqh
       zldis = hu - displa_lays(3)

       IF(zldis <= 0.0) THEN
          write(6,*) 'the obs height of u less than the zero displacement heght'
          CALL abort
       ENDIF

       CALL moninobukini(ur,th,thm,thv,dth,dqh,dthv,zldis,z0mv,um,obu)

! ======================================================================
!      BEGIN stability iteration 
! ======================================================================

       DO WHILE (it .le. itmax) 

          !print *,"iteration index:", it
          tlbef = tl

          del2  = del
          dele2 = dele
 
!-----------------------------------------------------------------------
! Aerodynamical resistances
!-----------------------------------------------------------------------
! Evaluate stability-dependent variables using moz from prior iteration

          CALL moninobukm(hu,ht,hq,displa_lays(toplay),z0mv,z0hv,z0qv,obu,um, &
                          displa_lay(toplay),z0m_lay(toplay),ustar,fh2m,fq2m, &
                          htop_lay(toplay),fmtop,fm,fh,fq,fht,fqt,phih)
 
! Aerodynamic resistance
          ! 09/16/2017:
          ! note that for ram, it is the resistance from Href to z0mv+displa
          ! however, for rah and raw is only from Href to canopy effective 
          ! exchange height.
          ! so rah/raw is not comparable with that of 1D case
          ram = 1./(ustar*ustar/um)
          
          ! 05/02/2016: calculate resistance from the top layer (effective exchange 
          ! height) to reference height
          rah = 1./(vonkar/(fh-fht)*ustar) 
          raw = 1./(vonkar/(fq-fqt)*ustar) 

! update roughness length for sensible/latent heat
          z0hg = z0mg/exp(0.13 * (ustar*z0mg/1.5e-5)**0.45)
          z0qg = z0hg

          z0h_lays(0) = z0hg
          z0q_lays(0) = z0qg

          z0h_lays(1) = maxval(z0h_lays(0:1))
          z0h_lays(2) = maxval(z0h_lays(0:2))
          z0h_lays(3) = maxval(z0h_lays(0:3))

          z0q_lays(:) = z0h_lays(:)
          z0hv = z0h_lays(3)
          z0qv = z0q_lays(3)

! ......................................................................
! new method to calculate rd and ueffect
! the kernel part of 3d model
! ......................................................................

          ! initialization
          rd(:)  = 0.
          rd_(:) = 0.
          upplay = 0
          
          ! calculate canopy top wind speed (utop) and exchange coefficient (ktop)
          ! need to update each time as obu changed after each iteration
          utop = ustar/vonkar * fmtop
          ktop = vonkar * (htop_lay(toplay)-displa_lays(toplay)) * ustar / phih
          !print *, "ktop:", ktop, htop_lay(toplay),displa_lays(toplay),ustar,phih ! fordebug

          ! start layer loop
          DO i = toplay, 1, -1

             IF (fcover_lay(i)>0 .and. lsai_lay(i)>0) THEN

                IF (i .eq. toplay) THEN
                   utop_lay(i) = utop
                   ktop_lay(i) = ktop
                ELSE
                   ! calculate utop of this layer
                   utop_lay(i) = uprofile(ubot_lay(upplay), fcover_lays(upplay), bee, 0., &
                      z0mg, hbot_lay(upplay), htop_lay(i), htop_lay(i))

                   ! calculate ktop of this layer
                   ktop_lay(i) = kprofile(kbot_lay(upplay), fcover_lays(upplay), bee, 0., &
                      displa_lays(toplay)/htop_lay(toplay), &
                      hbot_lay(upplay), htop_lay(i), obug, ustarg, htop_lay(i)) 
  
                   ! areodynamic resistance between this layer top and above layer bottom
! 03/15/2020, yuan: TODO, vertical gaps between layers, fc = fcover_lays(upplay) or just 0? 
                   !rd(upplay) = rd(upplay) + kintegral(kbot_lay(upplay), fcover_lays(upplay), bee, 0., &
                   !   z0mg, displa_lays(toplay)/htop_lay(toplay), &
                   !   hbot_lay(upplay), htop_lay(i), obug, ustarg, hbot_lay(upplay), htop_lay(i) )

                   rd(upplay) = rd(upplay) + frd(kbot_lay(upplay), hbot_lay(upplay), htop_lay(i), &
                      hbot_lay(upplay), htop_lay(i), displa_lays(toplay)/htop_lay(toplay), &
                      z0h_g, obug, ustarg, z0mg, 0., bee, fcover_lays(upplay))

                ENDIF

                ! for robust check
                hbot_lay(i) = max(hbot_lay(i), displa_lays(i-1)+z0m_lays(i-1))
                
                ! wind speed at layer bottom
                ubot_lay(i) = uprofile(utop_lay(i), fcover_lay(i), bee, a_lay(i), &
                   z0mg, htop_lay(i), hbot_lay(i), hbot_lay(i))

                ! effective wind speed for rb calculation
                !ueff_lay(i) = uintegral(utop_lay(i), fcover_lay(i), bee, a_lay(i), &
                !   z0mg, htop_lay(i), hbot_lay(i))

                IF (it == 1) THEN

                   !REAL(r8) FUNCTION ueffect(utop, htop, hbot, &
                   !         z0mg, alpha, bee, fc)
                   ueff_lay_norm(i) = ueffect(1., htop_lay(i), hbot_lay(i), &
                      z0mg, a_lay(i), bee, fcover_lay(i))
                ENDIF
                ueff_lay(i) = utop_lay(i)*ueff_lay_norm(i)

                ! normalized eddy coefficient (K) at layer bottom
                kbot_lay(i) = kprofile(ktop_lay(i), fcover_lay(i), bee, a_lay(i), &
                   displa_lays(toplay)/htop_lay(toplay), &
                   htop_lay(i), hbot_lay(i), obug, ustarg, hbot_lay(i)) 

                ! areodynamic resistance from effective fluxes exchange height of 
                ! of this layer to the top of this layer
                IF (upplay > 0) THEN
                   !rd(upplay) = rd(upplay) + kintegral(ktop_lay(i), fcover_lay(i), bee, &
                   !   a_lay(i), z0mg, displa_lays(toplay)/htop_lay(toplay), &
                   !   htop_lay(i), hbot_lay(i), obug, ustarg, htop_lay(i), displa_lay(i)+z0m_lay(i))

                   !REAL(r8) FUNCTION frd(ktop, htop, hbot, &
                   !      ztop, zbot, displah, z0h, obu, ustar, &
                   !      z0mg, alpha, bee, fc)
                   rd(upplay) = rd(upplay) + frd(ktop_lay(i), htop_lay(i), hbot_lay(i), &
                      htop_lay(i), displa_lay(i)+z0m_lay(i), displa_lays(toplay)/htop_lay(toplay), &
                      z0h_g, obug, ustarg, z0mg, a_lay(i), bee, fcover_lay(i))
                ENDIF

                ! areodynamic resistance from effective 'sink' height of 
                ! this layer to the bottom of this layer
                ! can not be lower than hbot_lay, otherwise, may introduce
                ! some contradiction.
                !rd(i) = rd(i) + kintegral(ktop_lay(i), fcover_lay(i), bee, a_lay(i), &
                !   z0mg, displa_lays(toplay)/htop_lay(toplay), &
                   ! 11/20/2017: integrate to z0qg, ONLY right for hbot=0
                   ! otherwise change to hbot_lay(i)
                   !htop_lay(i), hbot_lay(i), obug, ustarg, displa_lay(i)+z0m_lay(i), z0qg)
                   !htop_lay(i), hbot_lay(i), obug, ustarg, displa_lay(i)+z0m_lay(i), hbot_lay(i) )
! 01/06/2020, yuan: adjust to 0 hbot
                !   htop_lay(i), hbot_lay(i), obug, ustarg, displa_lay(i)+z0m_lay(i), max(z0qg,hbot_lay(i)) )

                rd(i) = rd(i) + frd(ktop_lay(i), htop_lay(i), hbot_lay(i), &
                   displa_lay(i)+z0m_lay(i), max(z0qg,hbot_lay(i)), &
                   displa_lays(toplay)/htop_lay(toplay), z0h_g, obug, ustarg, &
                   z0mg, a_lay(i), bee, fcover_lay(i))

                upplay = i
                !print *, "rd(i)", rd(i), i  ! fordebug
             
             ENDIF
          ENDDO

! ......................................................................
! areodynamic resistance between ground and the upper layer bottom
! ......................................................................

          ! uncomment the below when the upper codes change to hbot_lay
          !rd(botlay) = rd(botlay) + kintegral(kbot_lay(botlay), fcover_lays(botlay), bee, 0., &
          !   z0mg, displa_lays(toplay)/htop_lay(toplay), &
          !   hbot_lay(botlay), z0qg, obug, ustarg, hbot_lay(botlay), z0qg )

          rd(botlay) = rd(botlay) + frd(kbot_lay(botlay), hbot_lay(botlay), z0qg, &
             hbot_lay(botlay), z0qg, displa_lays(toplay)/htop_lay(toplay), &
             z0h_g, obug, ustarg, z0mg, 0., bee, fcover_lays(botlay))
 
! ......................................................................
! Bulk boundary layer resistance of leaves
! ......................................................................
          rb(:) = 0.
        
          DO i = 1, npft
             IF (fcover(i)>0 .and. lsai(i)>0) THEN
                clev = canlev(i)
                cf = 0.01*sqrtdi(i)*sqrt(ueff_lay(clev))
                rb(i) = 1./cf
             ENDIF
          ENDDO

          ! 10/01/2017, back to 1D case, for test ONLY
          IF (rb_opt == 1) THEN
             uaf   = ustar
             cf    = 0.01*sqrtdi(2)/sqrt(uaf)
             rb(:) = 1/(cf*uaf)
          ENDIF

!         rd = 1./(csoilc*uaf)                 ! BATS legacy
!         w = exp(-0.5*(lai+sai))              ! Dickinson's modification :
!         csoilc = ( 1.-w + w*um/uaf)/rah      ! "rah" here is the resistance over
!         rd = 1./(csoilc*uaf)                 ! bare ground fraction

          ! 10/01/2017, back to 1D case, for test ONLY
          IF (rd_opt == 1 ) THEN
! modified by Xubin Zeng's suggestion at 08-07-2002
             uaf   = ustar
             w = exp(-(lai(2)+sai(2)))
             csoilcn = (vonkar/(0.13*(z0mg*uaf/1.5e-5)**0.45))*w + csoilc*(1.-w)
             rd(:) = 1./(csoilcn*uaf)
          ENDIF

!-----------------------------------------------------------------------
! stomatal resistances 
!-----------------------------------------------------------------------

          DO i = 1, npft
             IF(fcover(i)>0 .and. lai(i)>0.001) THEN

                rbsun = rb(i) / laisun(i)
                rbsha = rb(i) / laisha(i)
          
                clev = canlev(i)
                eah = qaf(clev) * psrf / ( 0.622 + 0.378 * qaf(clev) )    !pa
               
! note: calculate resistance for sunlit/shaded leaves
!-----------------------------------------------------------------------
                CALL stomata (vmax25(i)   ,effcon(i) ,slti(i)   ,hlti(i)   ,&
                   shti(i)    ,hhti(i)    ,trda(i)   ,trdm(i)   ,trop(i)   ,&
                   gradm(i)   ,binter(i)  ,thm       ,psrf      ,po2m      ,&
                   pco2m      ,pco2a      ,eah       ,ei(i)     ,tl(i)     ,&
                   parsun(i)  ,rbsun      ,raw       ,rstfac(i) ,cintsun(:,i),&
                   assimsun(i),respcsun(i),rssun(i)     )

                CALL stomata (vmax25(i)   ,effcon(i) ,slti(i)   ,hlti(i)   ,&
                   shti(i)    ,hhti(i)    ,trda(i)   ,trdm(i)   ,trop(i)   ,&
                   gradm(i)   ,binter(i)  ,thm       ,psrf      ,po2m      ,&
                   pco2m      ,pco2a      ,eah       ,ei(i)     ,tl(i)     ,&
                   parsha(i)  ,rbsha      ,raw       ,rstfac(i) ,cintsha(:,i),&
                   assimsha(i),respcsha(i),rssha(i)     )
             ELSE
                rssun(i) = 2.e4; assimsun(i) = 0.; respcsun(i) = 0.
                rssha(i) = 2.e4; assimsha(i) = 0.; respcsha(i) = 0.
             ENDIF
          ENDDO

! above stomatal resistances are for the canopy, the stomatal rsistances 
! and the "rb" in the following calculations are the average for single leaf. thus,
          rssun = rssun * laisun
          rssha = rssha * laisha

!-----------------------------------------------------------------------
! dimensional and non-dimensional sensible and latent heat conductances
! for canopy and soil flux calculations.
!-----------------------------------------------------------------------
          
          cfh(:) = 0.
          cfw(:) = 0.
 
          DO i = 1, npft
             IF (fcover(i)>0 .and. lsai(i)>0) THEN
 
                clev = canlev(i)
                delta(i) = 0.0
                IF(qsatl(i)-qaf(clev) .gt. 0.) delta(i) = 1.0
 
                cfh(i) = lsai(i) / rb(i)
 
! note: combine sunlit and shaded leaves
!-----------------------------------------------------------------------
                cfw(i) = (1.-delta(i)*(1.-fwet(i)))*lsai(i)/rb(i) + &
                         (1.-fwet(i))*delta(i)* &
                         ( laisun(i)/(rb(i)+rssun(i)) + laisha(i)/(rb(i)+rssha(i)) )
             ENDIF
          ENDDO
          
          ! initialization
          cah(:) = 0.
          caw(:) = 0.
          cgh(:) = 0.
          cgw(:) = 0.

          DO i = 1, nlay
             IF (fcover_lay(i)>0 .and. lsai_lay(i)>0) THEN
                IF (i == toplay) THEN
                   cah(i) = 1. / rah
                   caw(i) = 1. / raw
                ELSE
                   cah(i) = 1. / rd(i+1)
                   caw(i) = 1. / rd(i+1) 
                ENDIF

                cgh(i) = 1. / rd(i)
                cgw(i) = 1. / rd(i)
             ENDIF
          ENDDO
          
          ! claculate wtshi, wtsqi
          wtshi(:) = cah(:) + cgh(:)
          wtsqi(:) = caw(:) + cgw(:)
          
          DO i = 1, npft
             IF (fcover(i)>0 .and. lsai(i)>0) THEN
                clev = canlev(i)
                wtshi(clev) = wtshi(clev) + fcover(i)*cfh(i)
                wtsqi(clev) = wtsqi(clev) + fcover(i)*cfw(i)
             ENDIF
          ENDDO

          DO i = 1, nlay
             IF (fcover_lay(i)>0 .and. lsai_lay(i)>0) THEN
                wtshi(i) = 1./wtshi(i)
                wtsqi(i) = 1./wtsqi(i)
             ENDIF
          ENDDO

          wta0(:) = cah(:) * wtshi(:)
          wtg0(:) = cgh(:) * wtshi(:)
 
          wtaq0(:) = caw(:) * wtsqi(:)
          wtgq0(:) = cgw(:) * wtsqi(:)
          
          ! calculate wtl0, wtll, wtlq0, wtlql
          wtll(:)  = 0.
          wtlql(:) = 0.

          DO i = 1, npft
             IF (fcover(i)>0 .and. lsai(i)>0) THEN
                clev = canlev(i)

                wtl0(i)  = cfh(i) * wtshi(clev) * fcover(i)
                wtll(clev) = wtll(clev) + wtl0(i)*tl(i)
                
                wtlq0(i) = cfw(i) * wtsqi(clev) * fcover(i)
                wtlql(clev) = wtlql(clev) + wtlq0(i)*qsatl(i)
             ENDIF
          ENDDO
 
          ! to solve taf(:) and qaf(:)
          IF (numlay .eq. 1) THEN 

             taf(toplay) = wta0(toplay)*thm + wtg0(toplay)*tg + wtll(toplay)
             qaf(toplay) = wtaq0(toplay)*qm + wtgq0(toplay)*qg + wtlql(toplay)
             fact = 1.
             facq = 1.

          ENDIF
          
          IF (numlay .eq. 2) THEN 

             tmpw1 = wtg0(botlay)*tg + wtll(botlay)
             fact  = 1. - wtg0(toplay)*wta0(botlay)
             taf(toplay) = ( wta0(toplay)*thm + wtg0(toplay)*tmpw1 + wtll(toplay) ) /  fact

             tmpw1 = wtgq0(botlay)*qg + wtlql(botlay)
             facq  = 1. - wtgq0(toplay)*wtaq0(botlay)
             qaf(toplay) = ( wtaq0(toplay)*qm + wtgq0(toplay)*tmpw1 + wtlql(toplay) ) / facq
             
             taf(botlay) = wta0(botlay)*taf(toplay) + wtg0(botlay)*tg + wtll(botlay)
             qaf(botlay) = wtaq0(botlay)*qaf(toplay) + wtgq0(botlay)*qg + wtlql(botlay)

          ENDIF
         
          IF (numlay .eq. 3) THEN 

             tmpw1 = wta0(3)*thm + wtll(3)
             tmpw2 = wtg0(1)*tg + wtll(1)
             fact  = 1. - wta0(2)*wtg0(3) - wtg0(2)*wta0(1) 
             taf(2) = ( wta0(2)*tmpw1 + wtg0(2)*tmpw2 + wtll(2) ) / fact

             tmpw1 = wtaq0(3)*qm + wtlql(3)
             tmpw2 = wtgq0(1)*qg + wtlql(1)
             facq  = 1. - wtaq0(2)*wtgq0(3) - wtgq0(2)*wtaq0(1)
             qaf(2) = ( wtaq0(2)*tmpw1 + wtgq0(2)*tmpw2 + wtlql(2) ) / facq

             taf(1) = wta0(1)*taf(2) + wtg0(1)*tg + wtll(1)
             qaf(1) = wtaq0(1)*qaf(2) + wtgq0(1)*qg + wtlql(1)
             
             taf(3) = wta0(3)*thm + wtg0(3)*taf(2) + wtll(3)
             qaf(3) = wtaq0(3)*qm + wtgq0(3)*qaf(2) + wtlql(3)

          ENDIF
          
!-----------------------------------------------------------------------
! IR radiation, sensible and latent heat fluxes and their derivatives
!-----------------------------------------------------------------------
! the partial derivatives of areodynamical resistance are ignored 
! which cannot be determined analtically

! calculate L for each canopy layer
          L(:) = 0.
          DO i = 1, npft
             IF (fcover(i)>0 .and. lsai(i)>0) THEN
                clev = canlev(i)
                ! according to absorption = emissivity, fcover -> fshade
                L(clev) = L(clev) + fshade(i) * (1-thermk(i)) * stefnc * tl(i)**4
             ENDIF
          ENDDO

! calculate Ltd
          Ltd(:) = 0.
          Ltd(3) = thermk_lay(3) * tdn(4,3) * frl
          Ltd(2) = thermk_lay(2) * ( tdn(4,2)*frl + tdn(3,2)*(Ltd(3) + L(3)) )
          Ltd(1) = thermk_lay(1) * ( tdn(4,1)*frl + tdn(3,1)*(Ltd(3) + L(3)) + &
                              tdn(2,1)*(Ltd(2) + L(2)) )

! calculate Ld = Ltd + L
          Ld(0) = 0.
          Ld(4) = frl
          Ld(1:3) = Ltd + L

! calculate Lin = Ld * tdn          
          Lin(:) = matmul(Ld(:), tdn(:,:))

! calcilate Lg = (1-emg)*dlrad + emg*stefnc*tg**4
! dlrad = Lin(0)
          Lg = (1 - emg)*Lin(0) + emg*stefnc*tg**4

! calculate Ltu
          Ltu(1) = thermk_lay(1) * tup(0,1) * Lg
          Ltu(2) = thermk_lay(2) * ( tup(0,2)*Lg + tup(1,2)*(Ltu(1) + L(1)) )
          Ltu(3) = thermk_lay(3) * ( tup(0,3)*Lg + tup(1,3)*(Ltu(1) + L(1)) + &
                                     tup(2,3)*(Ltu(2) + L(2)) )

! calculate Lu = Ltu + L
          Lu(0) = Lg
          Lu(4) = 0.
          Lu(1:3) = Ltu + L

! calculate Lin = Lin + Lu*tup
          Lin(:) = Lin(:) + matmul(Lu(:), tup(:,:))

! calculate Lv
          Lv(:) = 0.
          DO i = 1, npft
             IF (fshade(i) > 0) THEN
                clev  = canlev(i)
                Lv(i) = fshade(i)/fshade_lay(clev) * (1-thermk(i)) * Lin(clev) / fcover(i) &
                      - 2. * fshade(i) * (1-thermk(i)) * stefnc * tl(i)**4 / fcover(i)
             ENDIF
          ENDDO

! calculate delata(Lv)
          dLv(:) = 0.
          DO i = 1, npft
             IF (fshade(i) > 0) THEN
                clev   = canlev(i)
                dLv(i) = (4.*dLvpar(clev)*(1-emg)*fshade(i)*(1-thermk(i)) - 8.) &
                       * fshade(i) * (1-thermk(i)) * stefnc *  tl(i)**3 / fcover(i)
             ENDIF
          ENDDO

!-----------------------------------------------------------------------

          irab(:)      = Lv(:)
          dirab_dtl(:) = dLv(:)

          DO i = 1, npft
            
             IF (fcover(i)>0 .and. lsai(i)>0) THEN

                clev = canlev(i)
                fac(i) = 1. - thermk(i)

! sensible heat fluxes and their derivatives
                fsenl(i) = rhoair * cpair * cfh(i) * (tl(i) - taf(clev))
                
                ! 09/24/2017: why fact/facq here? bugs? YES
                ! 09/25/2017: re-written, check it clearfully 
                IF (numlay < 3 .or. clev == 2) THEN
                   fsenl_dtl(i) = rhoair * cpair * cfh(i) * (1 - wtl0(i)/fact)
                ELSE
                   IF (clev == 1) THEN
                      fsenl_dtl(i) = rhoair * cpair * cfh(i) * &
                         (1 - (1-wta0(2)*wtg0(3))*wtl0(i)/fact)
                   ENDIF
                   IF (clev == 3) THEN
                      fsenl_dtl(i) = rhoair * cpair * cfh(i) * &
                         (1 - (1-wtg0(2)*wta0(1))*wtl0(i)/fact)
                   ENDIF
                ENDIF

! latent heat fluxes and their derivatives
                etr(i) = rhoair * (1.-fwet(i)) * delta(i) &
                       * ( laisun(i)/(rb(i)+rssun(i)) + laisha(i)/(rb(i)+rssha(i)) ) &
                       * ( qsatl(i) - qaf(clev) )
                
                ! 09/25/2017: re-written 
                IF (numlay < 3 .or. clev == 2) THEN
                   etr_dtl(i) = rhoair * (1.-fwet(i)) * delta(i) &
                      * ( laisun(i)/(rb(i)+rssun(i)) + laisha(i)/(rb(i)+rssha(i)) ) &
                      * (1 - wtlq0(i)/facq)*qsatlDT(i) 
                ELSE
                   IF (clev == 1) THEN
                      etr_dtl(i) = rhoair * (1.-fwet(i)) * delta(i) &
                         * ( laisun(i)/(rb(i)+rssun(i)) + laisha(i)/(rb(i)+rssha(i)) ) &
                         * (1 - (1-wtaq0(2)*wtgq0(3))*wtlq0(i)/facq)*qsatlDT(i) 
                   ENDIF
                   IF (clev == 3) THEN
                      etr_dtl(i) = rhoair * (1.-fwet(i)) * delta(i) &
                         * ( laisun(i)/(rb(i)+rssun(i)) + laisha(i)/(rb(i)+rssha(i)) ) &
                         * (1 - (1-wtgq0(2)*wtaq0(1))*wtlq0(i)/facq)*qsatlDT(i) 
                   ENDIF
                ENDIF

                IF(etr(i).ge.etrc(i))THEN
                   etr(i) = etrc(i)
                   etr_dtl(i) = 0.
                ENDIF
 
                evplwet(i) = rhoair * (1.-delta(i)*(1.-fwet(i))) * lsai(i)/rb(i) &
                           * ( qsatl(i) - qaf(clev) )

                ! 09/25/2017: re-written 
                IF (numlay < 3 .or. clev == 2) THEN
                   evplwet_dtl(i) = rhoair * (1.-delta(i)*(1.-fwet(i))) * lsai(i)/rb(i) &
                      * (1 - wtlq0(i)/facq)*qsatlDT(i) 
                ELSE
                   IF (clev == 1) THEN
                      evplwet_dtl(i) = rhoair * (1.-delta(i)*(1.-fwet(i))) * lsai(i)/rb(i) &
                         * (1 - (1-wtaq0(2)*wtgq0(3))*wtlq0(i)/facq)*qsatlDT(i) 
                   ENDIF
                   IF (clev == 3) THEN
                      evplwet_dtl(i) = rhoair * (1.-delta(i)*(1.-fwet(i))) * lsai(i)/rb(i) &
                         * (1 - (1-wtgq0(2)*wtaq0(1))*wtlq0(i)/facq)*qsatlDT(i) 
                   ENDIF
                ENDIF
 
                ! 03/02/2018: convert evplwet from fc to whole area
                ! because ldew right now is for the whole area
                ! 09/05/2019: back to fc area
                IF(evplwet(i).ge.ldew(i)/deltim)THEN
                   evplwet(i) = ldew(i)/deltim
! 01/07/2020, yuan: bug? 不会产生计算的错误
                   evplwet_dtl(i) = 0.
                ENDIF
  
                fevpl(i) = etr(i) + evplwet(i)
                fevpl_dtl(i) = etr_dtl(i) + evplwet_dtl(i)
 
                erre(i) = 0.
                fevpl_noadj(i) = fevpl(i)
                IF ( fevpl(i)*fevpl_bef(i) < 0. ) THEN 
                   erre(i)  = -0.9*fevpl(i)
                   fevpl(i) =  0.1*fevpl(i)
                ENDIF
 
!-----------------------------------------------------------------------
! difference of temperatures by quasi-newton-raphson method for the non-linear system equations
!-----------------------------------------------------------------------

                dtl(it,i) = (sabv(i) + irab(i) - fsenl(i) - hvap*fevpl(i)) &
                   / (lsai(i)*clai/deltim - dirab_dtl(i) + fsenl_dtl(i) + hvap*fevpl_dtl(i))
                dtl_noadj(i) = dtl(it,i)
 
                ! check magnitude of change in leaf temperature limit to maximum allowed value
 
                IF(it .le. itmax) THEN
 
                  ! put brakes on large temperature excursions
                   IF(abs(dtl(it,i)).gt.delmax)THEN
                      dtl(it,i) = delmax*dtl(it,i)/abs(dtl(it,i))
                   ENDIF
 
                   ! 09/26/2017
                   ! could be a bug IF dtl*dtl==0? 
                   ! should le -> lt
                   IF((it.ge.2) .and. (dtl(it-1,i)*dtl(it,i).le.0.))THEN
                      dtl(it,i) = 0.5*(dtl(it-1,i) + dtl(it,i))
                   ENDIF
 
                ENDIF
 
                tl(i) = tlbef(i) + dtl(it,i)

!-----------------------------------------------------------------------
! square roots differences of temperatures and fluxes for use as the condition of convergences
!-----------------------------------------------------------------------

                del(i)  = sqrt( dtl(it,i)*dtl(it,i) )
                dele(i) = dtl(it,i) * dtl(it,i) * &
                   ( dirab_dtl(i)**2 + fsenl_dtl(i)**2 + hvap*fevpl_dtl(i)**2 ) 
                dele(i) = sqrt(dele(i))
 
!-----------------------------------------------------------------------
!  saturated vapor pressures and canopy air temperature, canopy air humidity
!-----------------------------------------------------------------------
! Recalculate leaf saturated vapor pressure (ei_)for updated leaf temperature
! and adjust specific humidity (qsatl_) proportionately
                CALL qsadv(tl(i),psrf,ei(i),deiDT(i),qsatl(i),qsatlDT(i))
            
             ENDIF
          ENDDO !end pft loop

! update vegetation/ground surface temperature, canopy air temperature, 
! canopy air humidity
          
          ! calculate wtll, wtlql
          wtll(:)  = 0.
          wtlql(:) = 0.

          DO i = 1, npft
             IF (fcover(i)>0 .and. lsai(i)>0) THEN
                clev = canlev(i)
                wtll(clev)  =  wtll(clev) +  wtl0(i)*tl(i)
                wtlql(clev) = wtlql(clev) + wtlq0(i)*qsatl(i)
             ENDIF
          ENDDO
 
          IF (numlay .eq. 1) THEN 

             taf(toplay) =  wta0(toplay)*thm + wtg0(toplay)*tg + wtll(toplay)
             qaf(toplay) = wtaq0(toplay)*qm + wtgq0(toplay)*qg + wtlql(toplay)

          ENDIF
          
          IF (numlay .eq. 2) THEN 

             tmpw1 = wtg0(botlay)*tg + wtll(botlay)
             fact  = 1. - wtg0(toplay)*wta0(botlay)
             taf(toplay) = (wta0(toplay)*thm + wtg0(toplay)*tmpw1 + wtll(toplay)) / fact

             tmpw1 = wtgq0(botlay)*qg + wtlql(botlay)
             facq  = 1. - wtgq0(toplay)*wtaq0(botlay)
             qaf(toplay) = (wtaq0(toplay)*qm + wtgq0(toplay)*tmpw1 + wtlql(toplay)) / facq
             
             taf(botlay) =  wta0(botlay)*taf(toplay) +  wtg0(botlay)*tg +  wtll(botlay)
             qaf(botlay) = wtaq0(botlay)*qaf(toplay) + wtgq0(botlay)*qg + wtlql(botlay)

          ENDIF
         
          IF (numlay .eq. 3) THEN 

             tmpw1 = wta0(3)*thm + wtll(3)
             tmpw2 = wtg0(1)*tg + wtll(1)
             fact  = 1. - wta0(2)*wtg0(3) - wtg0(2)*wta0(1) 
             taf(2) = (wta0(2)*tmpw1 + wtg0(2)*tmpw2 + wtll(2)) / fact

             tmpw1 = wtaq0(3)*qm + wtlql(3)
             tmpw2 = wtgq0(1)*qg + wtlql(1)
             facq  = 1. - wtaq0(2)*wtgq0(3) - wtgq0(2)*wtaq0(1)
             qaf(2) = (wtaq0(2)*tmpw1 + wtgq0(2)*tmpw2 + wtlql(2)) / facq

             taf(1) =  wta0(1)*taf(2) +  wtg0(1)*tg + wtll(1)
             qaf(1) = wtaq0(1)*qaf(2) + wtgq0(1)*qg + wtlql(1)
             
             taf(3) =  wta0(3)*thm + wtg0(3)*taf(2) + wtll(3)
             qaf(3) = wtaq0(3)*qm + wtgq0(3)*qaf(2) + wtlql(3)

          ENDIF

          ! fordebug
          !print *, "taf", taf
 
! update co2 partial pressure within canopy air
          ! 05/02/2016: may have some problem with gdh2o, however,
          ! this variable seems never used here. Different height
          ! level vegetation should have different gdh2o, i.e., 
          ! different rd(layer) values.
          gah2o = 1.0/raw * tprcor/thm                     !mol m-2 s-1
          gdh2o = 1.0/rd(botlay) * tprcor/thm              !mol m-2 s-1
          
          pco2a = pco2m - 1.37*psrf/max(0.446,gah2o) * &
                  sum(fcover*(assimsun + assimsha - respcsun - respcsha - rsoil))

!-----------------------------------------------------------------------
! Update monin-obukhov length and wind speed including the stability effect
!-----------------------------------------------------------------------

          ! 这里使用的是最高层的taf和qaf
          ! 如何进行限制?是不是梯度太大的问题?运行单点模型测试
          dth = thm - taf(toplay)
          dqh =  qm - qaf(toplay)
 
          tstar = vonkar/(fh-fht)*dth
          qstar = vonkar/(fq-fqt)*dqh
 
          thvstar = tstar + 0.61*th*qstar
          zeta = zldis*vonkar*grav*thvstar / (ustar**2*thv)
          IF(zeta .ge. 0.)THEN                             !stable
             zeta = min(2.,max(zeta,1.e-6))
          ELSE                                             !unstable
             zeta = max(-100.,min(zeta,-1.e-6))
          ENDIF
          obu = zldis/zeta
 
          IF(zeta .ge. 0.)THEN
             um = max(ur,.1)
          ELSE
             wc = (-grav*ustar*thvstar*zii/thv)**(1./3.)
             wc2 = beta*beta*(wc*wc)
             um = sqrt(ur*ur+wc2)
          ENDIF
 
          IF(obuold*obu .lt. 0.) nmozsgn = nmozsgn+1
          IF(nmozsgn .ge. 4) obu = zldis/(-0.01)
          obuold = obu

!-----------------------------------------------------------------------
! Test for convergence
!-----------------------------------------------------------------------

          it = it+1
 
          IF(it .gt. itmin) THEN
             fevpl_bef = fevpl
             det = maxval(max(del,del2))
             ! 10/03/2017: possible diff
             dee = maxval(max(dele,dele2))
             IF(det .lt. dtmin .and. dee .lt. dlemin) EXIT 
          ENDIF
 
       ENDDO 

       !IF (it > itmax) print *, "*** NOTE: it = 41! ***"

! ======================================================================
!     END stability iteration 
! ======================================================================

       z0m  = z0mv
       zol  = zeta
       rib  = min(5.,zol*ustar**2/(vonkar**2/fh*um**2))

! canopy fluxes and total assimilation amd respiration

       DO i = 1, npft
          IF (fcover(i)>0 .and. lsai(i)>0) THEN
             
             IF(lai(i) .gt. 0.001) THEN
                rst(i) = 1./(laisun(i)/rssun(i) + laisha(i)/rssha(i))
             ELSE
                rssun(i) = 2.0e4 ; rssha(i) = 2.0e4
                assimsun(i) = 0. ; assimsha(i) = 0.
                respcsun(i) = 0. ; respcsha(i) = 0.
                rst(i) = 2.0e4
             ENDIF
             assim(i) = assimsun(i) + assimsha(i)
             respc(i) = respcsun(i) + respcsha(i) + rsoil

! canopy fluxes and total assimilation amd respiration
             fsenl(i) = fsenl(i) + fsenl_dtl(i)*dtl(it-1,i) &
               ! add the imbalanced energy below due to T adjustment to sensibel heat
               + (dtl_noadj(i)-dtl(it-1,i)) * (lsai(i)*clai/deltim - dirab_dtl(i) &
               + fsenl_dtl(i) + hvap*fevpl_dtl(i)) &
               ! add the imbalanced energy below due to q adjustment to sensibel heat
               + hvap*erre(i)
 
             etr(i)     = etr(i)     +     etr_dtl(i)*dtl(it-1,i)
             evplwet(i) = evplwet(i) + evplwet_dtl(i)*dtl(it-1,i)
             fevpl(i)   = fevpl_noadj(i)
             fevpl(i)   = fevpl(i)   +   fevpl_dtl(i)*dtl(it-1,i)
 
             elwmax  = ldew(i)/deltim
             
             ! 03/02/2018: convert fc to whole area
             ! because ldew now is for the whole area
             ! may need to change to canopy covered area
             ! 09/14/2019: change back to canopy area
             elwdif  = max(0., evplwet(i)-elwmax)
             evplwet(i) = min(evplwet(i), elwmax)
             
             fevpl(i) = fevpl(i) - elwdif
             fsenl(i) = fsenl(i) + hvap*elwdif
             
!-----------------------------------------------------------------------
! Update dew accumulation (kg/m2)
!-----------------------------------------------------------------------

             ldew(i) = max(0., ldew(i)-evplwet(i)*deltim)

!-----------------------------------------------------------------------
! balance check
! (the computational error was created by the assumed 'dtl' in line 406-408) 
!-----------------------------------------------------------------------

             err = sabv(i) + irab(i) + dirab_dtl(i)*dtl(it-1,i) &
                 - fsenl(i) - hvap*fevpl(i)

#if(defined CLMDEBUG)
             IF(abs(err) .gt. .2) &
                write(6,*) 'energy imbalance in LeafTempPC.F90', &
                           i,it-1,err,sabv(i),irab(i),fsenl(i),hvap*fevpl(i)
#endif
 
          ENDIF
       ENDDO

!-----------------------------------------------------------------------
! downward (upward) longwave radiation below (above) the canopy
!-----------------------------------------------------------------------
       dlrad = Lin(0) &
             + sum( 4.* fshade * (1-thermk) * stefnc * tlbef**3 * dtl(it-1,:) )
       ulrad = Lin(4) - sum( fcover * dLv * dtl(it-1,:) ) &
             - emg * sum( 4.* fshade * (1-thermk) * stefnc * tlbef**3 * dtl(it-1,:) )

!-----------------------------------------------------------------------
! wind stresses 
!-----------------------------------------------------------------------

       taux = - rhoair*us/ram
       tauy = - rhoair*vs/ram

!-----------------------------------------------------------------------
! fluxes from ground to canopy space
!-----------------------------------------------------------------------

! for check purpose ONLY
! taf = wta0*thm + wtg0*tg + wtl0*tl 
! taf(1) = wta0(1)*taf(2) + wtg0(1)*tg + wtll(1)
! qaf(1) = wtaq0(1)*qaf(2) + wtgq0(1)*qg + wtlql(1)
! taf(botlay) = wta0(botlay)*taf(toplay) + wtg0(botlay)*tg + wtll(botlay)
! qaf(botlay) = wtaq0(botlay)*qaf(toplay) + wtgq0(botlay)*qg + wtlql(botlay)
! taf(toplay) = wta0(toplay)*thm +  wtg0(toplay)*tg + wtll(toplay)
! qaf(toplay) = wtaq0(toplay)*qm + wtgq0(toplay)*qg + wtlql(toplay)
       
       fseng = cpair*rhoair*cgh(botlay)*(tg-taf(botlay))
       fevpg = rhoair*cgw(botlay)*(qg-qaf(botlay))

!-----------------------------------------------------------------------
! Derivative of soil energy flux with respect to soil temperature (cgrnd)
!-----------------------------------------------------------------------

       cgrnds = cpair*rhoair*cgh(botlay)*(1.-wtg0(botlay))
       cgrndl = rhoair*cgw(botlay)*(1.-wtgq0(botlay))*dqgdT
       cgrnd  = cgrnds + cgrndl*htvp

!-----------------------------------------------------------------------
! 2 m height air temperature
!-----------------------------------------------------------------------

       tref = thm + vonkar/(fh-fht)*dth * (fh2m/vonkar - fh/vonkar) 
       qref =  qm + vonkar/(fq-fqt)*dqh * (fq2m/vonkar - fq/vonkar)

  END SUBROUTINE LeafTempPC
!----------------------------------------------------------------------         



  SUBROUTINE dewfraction (sigf,lai,sai,dewmx,ldew,fwet,fdry)
       
!=======================================================================
! Original author: Yongjiu Dai, September 15, 1999
!
! determine fraction of foliage covered by water and
! fraction of foliage that is dry and transpiring
!
!=======================================================================

 USE precision
 IMPLICIT NONE

  REAL(r8), intent(in) :: sigf   !fraction of veg cover, excluding snow-covered veg [-]
  REAL(r8), intent(in) :: lai    !leaf area index  [-]
  REAL(r8), intent(in) :: sai    !stem area index  [-]
  REAL(r8), intent(in) :: dewmx  !maximum allowed dew [0.1 mm]
  REAL(r8), intent(in) :: ldew   !depth of water on foliage [kg/m2/s]

  REAL(r8), intent(out) :: fwet  !fraction of foliage covered by water [-]
  REAL(r8), intent(out) :: fdry  !fraction of foliage that is green and dry [-]

  REAL(r8) lsai                  !lai + sai
  REAL(r8) dewmxi                !inverse of maximum allowed dew [1/mm]
  REAL(r8) vegt                  !sigf*lsai
!
!-----------------------------------------------------------------------
! Fwet is the fraction of all vegetation surfaces which are wet 
! including stem area which contribute to evaporation
      lsai = lai + sai
      dewmxi = 1.0/dewmx
      ! why * sigf? may have bugs
      ! 06/17/2018:
      ! for ONLY one PFT, there may be no problem
      ! but for multiple PFTs, bugs exist!!!
      ! convert the whole area ldew to sigf ldew
      vegt   =  lsai

      fwet = 0
      IF(ldew > 0.) THEN
         fwet = ((dewmxi/vegt)*ldew)**.666666666666

! Check for maximum limit of fwet
         fwet = min(fwet,1.0)

      ENDIF

! fdry is the fraction of lai which is dry because only leaves can 
! transpire. Adjusted for stem area which does not transpire
      fdry = (1.-fwet)*lai/lsai


  END SUBROUTINE dewfraction
!----------------------------------------------------------------------         

  REAL(r8) FUNCTION uprofile(utop, fc, bee, alpha, z0mg, htop, hbot, z)
     
     USE precision
     USE FRICTION_VELOCITY
     IMPLICIT NONE

     REAL(r8), intent(in) :: utop
     REAL(r8), intent(in) :: fc
     REAL(r8), intent(in) :: bee 
     REAL(r8), intent(in) :: alpha
     REAL(r8), intent(in) :: z0mg
     REAL(r8), intent(in) :: htop
     REAL(r8), intent(in) :: hbot
     REAL(r8), intent(in) :: z
     
     REAL(r8) :: ulog,uexp

     ! when canopy LAI->0, z0->zs, fac->1, u->umoninobuk
     ! canopy LAI->large, fac->0 or=0, u->log profile
     ulog = utop*log(z/z0mg)/log(htop/z0mg)
     uexp = utop*exp(-alpha*(1-(z-hbot)/(htop-hbot)))
     
     uprofile = bee*fc*min(uexp,ulog) + (1-bee*fc)*ulog
     
     RETURN
  END FUNCTION uprofile

  REAL(r8) FUNCTION kprofile(ktop, fc, bee, alpha, &
                    displah, htop, hbot, obu, ustar, z)
     
     USE precision
     USE FRICTION_VELOCITY
     IMPLICIT NONE
     
     REAL(r8), parameter :: com1 = 0.4
     REAL(r8), parameter :: com2 = 0.08
     
     REAL(r8), intent(in) :: ktop
     REAL(r8), intent(in) :: fc
     REAL(r8), intent(in) :: bee
     REAL(r8), intent(in) :: alpha
     REAL(r8), intent(in) :: displah
     REAL(r8), intent(in) :: htop
     REAL(r8), intent(in) :: hbot
     REAL(r8), intent(in) :: obu
     REAL(r8), intent(in) :: ustar
     REAL(r8), intent(in) :: z
     
     REAL(r8) :: fac
     REAL(r8) :: kcob, klin, kexp

     klin = ktop*z/htop
     
     fac  = 1. / (1.+exp(-(displah-com1)/com2))
     kcob = 1. / (fac/klin + (1.-fac)/kmoninobuk(0.,obu,ustar,z))
     
     kexp = ktop*exp(-alpha*(1-(z-hbot)/(htop-hbot)))

     kprofile = 1./( bee*fc/min(kexp,kcob) + (1-bee*fc)/kcob )

     RETURN
  END FUNCTION kprofile

  REAL(r8) FUNCTION uintegral(utop, fc, bee, alpha, z0mg, htop, hbot)

     USE precision
     IMPLICIT NONE
     
     REAL(r8), intent(in) :: utop
     REAL(r8), intent(in) :: fc
     REAL(r8), intent(in) :: bee
     REAL(r8), intent(in) :: alpha
     REAL(r8), intent(in) :: z0mg
     REAL(r8), intent(in) :: htop
     REAL(r8), intent(in) :: hbot
     
     INTEGER  :: i, n
     REAL(r8) :: dz, z, u

     ! 09/26/2017: change fixed n -> fixed dz
     dz = 0.01 !fordebug
     n  = int( (htop-hbot) / dz ) + 1

     uintegral = 0.

     DO i = 1, n
        IF (i < n) THEN
           z = htop - (i-0.5)*dz
        ELSE
           dz = htop - hbot - (n-1)*dz
           z  = hbot + 0.5*dz
        ENDIF

        u = uprofile(utop, fc, bee, alpha, z0mg, htop, hbot, z)

        u = max(0._r8, u)
        !uintegral = uintegral + sqrt(u)*dz / (htop-hbot)
        ! 03/04/2020, yuan: TODO-hard to solve
        ! u开根号后不能解析求解积分，可近似直接对u积分
        ! 如此，最后就不用平方
        uintegral = uintegral + u*dz / (htop-hbot)
     ENDDO
      
     !uintegral = uintegral * uintegral

     RETURN
  END FUNCTION uintegral


  REAL(r8) FUNCTION ueffect(utop, htop, hbot, &
                            z0mg, alpha, bee, fc)
     USE precision
     IMPLICIT NONE
     
     REAL(r8), intent(in) :: utop
     REAL(r8), intent(in) :: htop
     REAL(r8), intent(in) :: hbot
     REAL(r8), intent(in) :: z0mg
     REAL(r8), intent(in) :: alpha
     REAL(r8), intent(in) :: bee
     REAL(r8), intent(in) :: fc
     
     REAL(r8) :: roots(2), uint
     INTEGER  :: rootn

     rootn = 0
     uint  = 0.

     CALL ufindroots(htop,hbot,(htop+hbot)/2., &
        utop, htop, hbot, z0mg, alpha, roots, rootn)

     IF (rootn == 0) THEN !no root
        uint = uint + fuint(utop, htop, hbot, &
           htop, hbot, z0mg, alpha, bee, fc)
     ENDIF
     
     IF (rootn == 1) THEN
        uint = uint + fuint(utop, htop, roots(1), &
           htop, hbot, z0mg, alpha, bee, fc)
        uint = uint + fuint(utop, roots(1), hbot, &
           htop, hbot, z0mg, alpha, bee, fc)
     ENDIF

     IF (rootn == 2) THEN
        uint = uint + fuint(utop, htop,     roots(1), &
           htop, hbot, z0mg, alpha, bee, fc)
        uint = uint + fuint(utop, roots(1), roots(2), &
           htop, hbot, z0mg, alpha, bee, fc)
        uint = uint + fuint(utop, roots(2), hbot,     &
           htop, hbot, z0mg, alpha, bee, fc)
     ENDIF 

     ueffect = uint / (htop-hbot)

     RETURN
  END FUNCTION ueffect


  REAL(r8) FUNCTION fuint(utop, ztop, zbot, &
        htop, hbot, z0mg, alpha, bee, fc)

     USE precision
     IMPLICIT NONE

     REAL(r8), intent(in) :: utop, ztop, zbot
     REAL(r8), intent(in) :: htop, hbot
     REAL(r8), intent(in) :: z0mg, alpha
     REAL(r8), intent(in) :: bee, fc

     ! local variables
     REAL(r8) :: fuexpint, fulogint
     
     fulogint = utop/log(htop/z0mg) *&
        (ztop*log(ztop/z0mg) - zbot*log(zbot/z0mg) + zbot - ztop)

     IF (udif((ztop+zbot)/2.,utop,htop,hbot,z0mg,alpha) <= 0) THEN
        ! uexp is smaller
        fuexpint = utop*(htop-hbot)/alpha*( &
           exp(-alpha*(htop-ztop)/(htop-hbot))-&
           exp(-alpha*(htop-zbot)/(htop-hbot)) )

        fuint = bee*fc*fuexpint + (1.-bee*fc)*fulogint
     ELSE
        ! ulog is smaller
        fuint = fulogint
     ENDIF
     
     RETURN
  END FUNCTION fuint


  RECURSIVE SUBROUTINE ufindroots(ztop,zbot,zmid, &
     utop, htop, hbot, z0mg, alpha, roots, rootn)

     USE precision
     IMPLICIT NONE 

     REAL(r8), intent(in) :: ztop, zbot, zmid
     REAL(r8), intent(in) :: utop, htop, hbot
     REAL(r8), intent(in) :: z0mg, alpha
     
     REAL(r8), intent(inout) :: roots(2)
     INTEGER,  intent(inout) :: rootn

     ! local variables
     REAL(r8) :: udif_ub, udif_lb
     
     udif_ub = udif(ztop,utop,htop,hbot,z0mg,alpha)
     udif_lb = udif(zmid,utop,htop,hbot,z0mg,alpha)

     IF (udif_ub*udif_lb == 0) THEN
        IF (udif_lb == 0) THEN !root found
           rootn = rootn + 1
           IF (rootn > 2) THEN 
              print *, "U root number > 2, abort!"
              CALL abort
           ENDIF
           roots(rootn) = zmid
        ENDIF
     ELSE IF (udif_ub*udif_lb < 0) THEN
        IF (ztop-zmid < 0.01) THEN
           rootn = rootn + 1 !root found
           IF (rootn > 2) THEN 
              print *, "U root number > 2, abort!"
              CALL abort
           ENDIF
           roots(rootn) = (ztop+zmid)/2.
        ELSE
           CALL ufindroots(ztop,zmid,(ztop+zmid)/2., &
              utop, htop, hbot, z0mg, alpha, roots, rootn)
        ENDIF
     ENDIF
     
     udif_ub = udif(zmid,utop,htop,hbot,z0mg,alpha)
     udif_lb = udif(zbot,utop,htop,hbot,z0mg,alpha)

     IF (udif_ub*udif_lb == 0) THEN
        IF (udif_ub == 0) THEN !root found
           rootn = rootn + 1
           IF (rootn > 2) THEN 
              print *, "U root number > 2, abort!"
              CALL abort
           ENDIF
           roots(rootn) = zmid
        ENDIF
     ELSE IF (udif_ub*udif_lb < 0) THEN
        IF (zmid-zbot < 0.01) THEN
           rootn = rootn + 1 !root found
           IF (rootn > 2) THEN 
              print *, "U root number > 2, abort!"
              CALL abort
           ENDIF
           roots(rootn) = (zmid+zbot)/2.
        ELSE
           CALL ufindroots(zmid,zbot,(zmid+zbot)/2., &
              utop, htop, hbot, z0mg, alpha, roots, rootn)
        ENDIF
     ENDIF

  END SUBROUTINE ufindroots


  REAL(r8) FUNCTION udif(z, utop, htop, hbot, z0mg, alpha)
     
     USE precision
     IMPLICIT NONE

     REAL(r8), intent(in) :: z, utop, htop, hbot
     REAL(r8), intent(in) :: z0mg, alpha

     REAL(r8) :: uexp, ulog

     uexp = utop*exp(-alpha*(1-(z-hbot)/(htop-hbot)))
     ulog = utop*log(z/z0mg)/log(htop/z0mg)

     udif = uexp - ulog

     RETURN
  END FUNCTION udif

  
  REAL(r8) FUNCTION kintegral(ktop, fc, bee, alpha, z0mg, &
                    displah, htop, hbot, obu, ustar, ztop, zbot)
     USE precision
     IMPLICIT NONE            
     
     REAL(r8), intent(in) :: ktop
     REAL(r8), intent(in) :: fc
     REAL(r8), intent(in) :: bee
     REAL(r8), intent(in) :: alpha
     REAL(r8), intent(in) :: z0mg
     REAL(r8), intent(in) :: displah
     REAL(r8), intent(in) :: htop
     REAL(r8), intent(in) :: hbot
     REAL(r8), intent(in) :: obu
     REAL(r8), intent(in) :: ustar
     REAL(r8), intent(in) :: ztop
     REAL(r8), intent(in) :: zbot
     
     INTEGER  :: i, n
     REAL(r8) :: dz, z, k
     
     kintegral = 0.

     IF (ztop <= zbot) THEN
        RETURN
     ENDIF
        
     ! 09/26/2017: change fixed n -> fixed dz
     ! 10/05/2017: need to improve
     dz = 0.01 !fordebug
     n  = int( (ztop-zbot) / dz ) + 1

     DO i = 1, n
        IF (i < n) THEN
           z  = ztop - (i-0.5)*dz
        ELSE
           dz = ztop - zbot - (n-1)*dz
           z  = zbot + 0.5*dz
        ENDIF

        k = kprofile(ktop, fc, bee, alpha, &
           displah, htop, hbot, obu, ustar, z) 

        kintegral = kintegral + 1./k * dz

     ENDDO
      
     RETURN
  END FUNCTION kintegral
  
  REAL(r8) FUNCTION frd(ktop, htop, hbot, &
        ztop, zbot, displah, z0h, obu, ustar, &
        z0mg, alpha, bee, fc)
     
     USE precision
     IMPLICIT NONE
     
     REAL(r8), intent(in) :: ktop, htop, hbot
     REAL(r8), intent(in) :: ztop, zbot
     REAL(r8), intent(in) :: displah, z0h, obu, ustar
     REAL(r8), intent(in) :: z0mg, alpha, bee, fc
     
     ! local parameters
     REAL(r8), parameter :: com1 = 0.4
     REAL(r8), parameter :: com2 = 0.08
     
     REAL(r8) :: roots(2), fac, kint
     INTEGER  :: rootn

     rootn = 0
     kint  = 0.

     ! calculate fac
     fac = 1. / (1.+exp(-(displah-com1)/com2))

     CALL kfindroots(ztop,zbot,(ztop+zbot)/2., &
        ktop, htop, hbot, obu, ustar, fac, alpha, roots, rootn)

     !print *, roots, rootn
     IF (rootn == 0) THEN !no root
        kint = kint + fkint(ktop, ztop, zbot, htop, hbot, &
           z0h, obu, ustar, fac, alpha, bee, fc)
     ENDIF
     
     IF (rootn == 1) THEN
        kint = kint + fkint(ktop, ztop, roots(1), htop, hbot, &
           z0h, obu, ustar, fac, alpha, bee, fc)
        kint = kint + fkint(ktop, roots(1), zbot, htop, hbot, &
           z0h, obu, ustar, fac, alpha, bee, fc)
     ENDIF

     IF (rootn == 2) THEN
        kint = kint + fkint(ktop, ztop, roots(1), htop, hbot, &
           z0h, obu, ustar, fac, alpha, bee, fc)
        kint = kint + fkint(ktop, roots(1), roots(2), htop, hbot, &
           z0h, obu, ustar, fac, alpha, bee, fc)
        kint = kint + fkint(ktop, roots(2), zbot, htop, hbot, &
           z0h, obu, ustar, fac, alpha, bee, fc)
     ENDIF 

     frd = kint

     RETURN
  END FUNCTION frd


  REAL(r8) FUNCTION fkint(ktop, ztop, zbot, htop, hbot, &
        z0h, obu, ustar, fac, alpha, bee, fc)

     USE precision
     USE FRICTION_VELOCITY
     IMPLICIT NONE

     REAL(r8), intent(in) :: ktop, ztop, zbot
     REAL(r8), intent(in) :: htop, hbot
     REAL(r8), intent(in) :: z0h, obu, ustar, fac, alpha
     REAL(r8), intent(in) :: bee, fc

     ! local variables
     REAL(r8) :: fkexpint, fkcobint
     
     !klin = ktop*z/htop
     !kcob = 1./(fac/klin + (1.-fac)/kmoninobuk(0.,obu,ustar,z))
     fkcobint = fac*htop/ktop*(log(ztop)-log(zbot)) +&
        (1.-fac)*kintmoninobuk(0.,z0h,obu,ustar,ztop,zbot)

     IF (kdif((ztop+zbot)/2.,ktop,htop,hbot,obu,ustar,fac,alpha) <= 0) THEN
        ! kexp is smaller
        IF (alpha > 0) THEN
           fkexpint = -(htop-hbot)/alpha/ktop*( &
              exp(alpha*(htop-ztop)/(htop-hbot))-&
              exp(alpha*(htop-zbot)/(htop-hbot)) )
        ELSE
           fkexpint = (ztop-zbot)/ktop
        ENDIF 

        fkint = bee*fc*fkexpint + (1.-bee*fc)*fkcobint
     ELSE
        ! kcob is smaller
        fkint = fkcobint
     ENDIF
     
     RETURN
  END FUNCTION fkint


  RECURSIVE SUBROUTINE kfindroots(ztop,zbot,zmid, &
     ktop, htop, hbot, obu, ustar, fac, alpha, roots, rootn)

     USE precision
     IMPLICIT NONE 

     REAL(r8), intent(in) :: ztop, zbot, zmid
     REAL(r8), intent(in) :: ktop, htop, hbot
     REAL(r8), intent(in) :: obu, ustar, fac, alpha
     
     REAL(r8), intent(inout) :: roots(2)
     INTEGER,  intent(inout) :: rootn

     ! local variables
     REAL(r8) :: kdif_ub, kdif_lb
     
     !print *, "*** CALL recursive SUBROUTINE kfindroots!!"
     kdif_ub = kdif(ztop,ktop,htop,hbot,obu,ustar,fac,alpha)
     kdif_lb = kdif(zmid,ktop,htop,hbot,obu,ustar,fac,alpha)

     IF (kdif_ub*kdif_lb == 0) THEN
        IF (kdif_lb == 0) THEN !root found
           rootn = rootn + 1
           IF (rootn > 2) THEN 
              print *, "K root number > 2, abort!"
              CALL abort
           ENDIF
           roots(rootn) = zmid
        ENDIF
     ELSE IF (kdif_ub*kdif_lb < 0) THEN
        IF (ztop-zmid < 0.01) THEN
           rootn = rootn + 1 !root found
           IF (rootn > 2) THEN 
              print *, "K root number > 2, abort!"
              CALL abort
           ENDIF
           roots(rootn) = (ztop+zmid)/2.
        ELSE
           CALL kfindroots(ztop,zmid,(ztop+zmid)/2., &
              ktop, htop, hbot, obu, ustar, fac, alpha, roots, rootn)
        ENDIF
     ENDIF
     
     kdif_ub = kdif(zmid,ktop,htop,hbot,obu,ustar,fac,alpha)
     kdif_lb = kdif(zbot,ktop,htop,hbot,obu,ustar,fac,alpha)

     IF (kdif_ub*kdif_lb == 0) THEN
        IF (kdif_ub == 0) THEN !root found
           rootn = rootn + 1
           IF (rootn > 2) THEN 
              print *, "K root number > 2, abort!"
              CALL abort
           ENDIF
           roots(rootn) = zmid
        ENDIF
     ELSE IF (kdif_ub*kdif_lb < 0) THEN
        IF (zmid-zbot < 0.01) THEN
           rootn = rootn + 1 !root found
           IF (rootn > 2) THEN 
              print *, "K root number > 2, abort!"
              CALL abort
           ENDIF
           roots(rootn) = (zmid+zbot)/2.
        ELSE
           CALL kfindroots(zmid,zbot,(zmid+zbot)/2., &
              ktop, htop, hbot, obu, ustar, fac, alpha, roots, rootn)
        ENDIF
     ENDIF

  END SUBROUTINE kfindroots


  REAL(r8) FUNCTION kdif(z, ktop, htop, hbot, &
        obu, ustar, fac, alpha)
     
     USE precision
     USE FRICTION_VELOCITY
     IMPLICIT NONE

     REAL(r8), intent(in) :: z, ktop, htop, hbot
     REAL(r8), intent(in) :: obu, ustar, fac, alpha 

     REAL(r8) :: kexp, klin, kcob

     kexp = ktop*exp(-alpha*(1-(z-hbot)/(htop-hbot)))
     
     klin = ktop*z/htop
     kcob = 1./(fac/klin + (1.-fac)/kmoninobuk(0.,obu,ustar,z))

     kdif = kexp - kcob

     RETURN
  END FUNCTION kdif


  SUBROUTINE cal_z0_displa (lai, h, fc, z0, displa)

     USE PhysicalConstants, only: vonkar
     IMPLICIT NONE

     REAL(r8), intent(in)  :: lai
     REAL(r8), intent(in)  :: h
     REAL(r8), intent(in)  :: fc
     REAL(r8), intent(out) :: z0
     REAL(r8), intent(out) :: displa

     REAL(r8), parameter :: Cd   = 0.2   !leaf drag coefficient
     REAL(r8), parameter :: cd1  = 7.5   !a free parameter for d/h calculation, Raupach 1992, 1994
     REAL(r8), parameter :: psih = 0.193 !psih = ln(cw) - 1 + cw^-1, cw = 2, Raupach 1994

     ! local variables
     REAL(r8) :: fai, sqrtdragc, temp1, delta , lai0

     ! when assume z0=0.01, displa=0
     ! to calculate lai0, delta displa
     !----------------------------------------------------
     sqrtdragc = -vonkar/(log(0.01/h) - psih)
     sqrtdragc = max(sqrtdragc, 0.0031**0.5)
     IF (sqrtdragc .le. 0.3) THEN
        fai = (sqrtdragc**2-0.003) / 0.3
        fai = min(fai, fc*(1-exp(-20.)))
     ELSE
        fai = 0.29
        print *, "z0m, displa error!"
     ENDIF

     ! calculate delta displa when z0 = 0.01
     lai0  = -log(1.-fai/fc)/0.5
     temp1 = (2.*cd1*fai)**0.5
     delta = -h * ( fc*1.1*log(1. + (Cd*lai0*fc)**0.25) + &
        (1.-fc)*(1.-(1.-exp(-temp1))/temp1) )
     
     ! calculate z0m, displa
     !----------------------------------------------------
     ! NOTE: potential bug below, ONLY apply for spheric
     ! crowns. For other cases, fc*(...) ==> a*fc*(...)
     fai   = fc*(1. - exp(-0.5*lai))
     sqrtdragc = min( (0.003+0.3*fai)**0.5, 0.3 )
     temp1 = (2.*cd1*fai)**0.5
     
     IF (lai > lai0) THEN
        displa = delta + h*( &
           (  fc)*1.1*log(1. + (Cd*lai*fc)**0.25) + &
           (1-fc)*(1.-(1.-exp(-temp1))/temp1) )
     ELSE 
        displa = h*( &
           (  fc)*1.1*log(1. + (Cd*lai*fc)**0.25) + &
           (1-fc)*(1.-(1.-exp(-temp1))/temp1) )
     ENDIF

     displa = max(displa, 0.)
     z0 = (h-displa) * exp(-vonkar/sqrtdragc + psih)

     IF (z0 < 0.01) THEN
        z0 = 0.01
        displa = 0.
     ENDIF

  END SUBROUTINE cal_z0_displa

END MODULE LEAF_temperature_PC
