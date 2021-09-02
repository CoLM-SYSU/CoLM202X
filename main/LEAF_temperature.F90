#include <define.h>

MODULE LEAF_temperature

!-----------------------------------------------------------------------
 USE precision
 IMPLICIT NONE
 SAVE

! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: LeafTemp

! PRIVATE MEMBER FUNCTIONS:
  PRIVATE :: dewfraction
  PRIVATE :: cal_z0_displa


!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE  LeafTemp(ipatch,deltim,csoilc,dewmx,htvp    ,lai     ,&
              sai     ,htop    ,hbot    ,sqrtdi  ,effcon  ,vmax25  ,&
              slti    ,hlti    ,shti    ,hhti    ,trda    ,trdm    ,&
              trop    ,gradm   ,binter  ,extkn   ,extkb   ,extkd   ,&
              hu      ,ht      ,hq      ,us      ,vs      ,thm     ,&
              th      ,thv     ,qm      ,psrf    ,rhoair  ,parsun  ,&
              parsha  ,sabv    ,frl     ,fsun    ,thermk  ,rstfac  ,&
              po2m    ,pco2m   ,z0h_g   ,obug    ,ustarg  ,zlnd    ,&
              zsno    ,fsno    ,sigf    ,etrc    ,tg      ,qg      ,&
              dqgdT   ,emg     ,tl      ,ldew    ,taux    ,tauy    ,&
              fseng   ,fevpg   ,cgrnd   ,cgrndl  ,cgrnds  ,tref    ,&
              qref    ,rst     ,assim   ,respc   ,fsenl   ,fevpl   ,&
              etr     ,dlrad   ,ulrad   ,z0m     ,zol     ,rib     ,&
              ustar   ,qstar   ,tstar   ,fm      ,fh      ,fq       ) 
 
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
  USE MOD_TimeInvariants, only: patchclass
  USE LC_Const, only: z0mr, displar
  IMPLICIT NONE
 
!-----------------------Arguments---------------------------------------

  INTEGER,  intent(in) :: ipatch
  REAL(r8), intent(in) :: &
        deltim,     &! seconds in a time step [second]
        csoilc,     &! drag coefficient for soil under canopy [-]
        dewmx,      &! maximum dew
        htvp         ! latent heat of evaporation (/sublimation) [J/kg]

! vegetation parameters
  REAL(r8), intent(in) :: &
        sai,        &! stem area index  [-]
        sqrtdi,     &! inverse sqrt of leaf dimension [m**-0.5]
        htop,       &! PFT crown top height [m]
        hbot,       &! PFT crown bot height [m]

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

        lai,        &! adjusted leaf area index for seasonal variation [-]
        parsun,     &! par absorbed per unit lai [w/m**2]
        parsha,     &! par absorbed per unit lai [w/m**2]
        sabv,       &! solar radiation absorbed by vegetation [W/m2]
        frl,        &! atmospheric infrared (longwave) radiation [W/m2]
        fsun,       &! sunlit fraction of canopy

        extkb,      &! (k, g(mu)/mu) direct solar extinction coefficient
        extkd,      &! diffuse and scattered diffuse PAR extinction coefficient
        thermk,     &! canopy gap fraction for tir radiation
        rstfac,     &! factor of soil water stress to plant physiologocal processes

        po2m,       &! atmospheric partial pressure  o2 (pa)
        pco2m,      &! atmospheric partial pressure co2 (pa)

        z0h_g,      &! bare soil roughness length, sensible heat [m]
        obug,       &! bare soil obu
        ustarg,     &! bare soil ustar
        zlnd,       &! roughness length for soil [m]
        zsno,       &! roughness length for snow [m]
        fsno,       &! fraction of snow cover on ground

        sigf,       &! fraction of veg cover, excluding snow-covered veg [-]
        etrc,       &! maximum possible transpiration rate (mm/s)
        tg,         &! ground surface temperature [K]
        qg,         &! specific humidity at ground surface [kg/kg]
        dqgdT,      &! temperature derivative of "qg"
        emg          ! vegetation emissivity

  REAL(r8), intent(inout) :: &
        tl,         &! leaf temperature [K]
        ldew,       &! depth of water on foliage [mm]
        taux,       &! wind stress: E-W [kg/m/s**2]
        tauy,       &! wind stress: N-S [kg/m/s**2]
        fseng,      &! sensible heat flux from ground [W/m2]
        fevpg,      &! evaporation heat flux from ground [mm/s]
        cgrnd,      &! deriv. of soil energy flux wrt to soil temp [w/m2/k]
        cgrndl,     &! deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
        cgrnds,     &! deriv of soil latent heat flux wrt soil temp [w/m**2/k]
        tref,       &! 2 m height air temperature (kelvin)
        qref         ! 2 m height air specific humidity

  REAL(r8), intent(out) :: &
        rst,        &! stomatal resistance
        assim,      &! rate of assimilation
        respc,      &! rate of respiration
        fsenl,      &! sensible heat from leaves [W/m2]
        fevpl,      &! evaporation+transpiration from leaves [mm/s]
        etr,        &! transpiration rate [mm/s]
        dlrad,      &! downward longwave radiation blow the canopy [W/m2]
        ulrad,      &! upward longwave radiation above the canopy [W/m2]

        z0m,        &! effective roughness [m]
        zol,        &! dimensionless height (z/L) used in Monin-Obukhov theory
        rib,        &! bulk Richardson number in surface layer
        ustar,      &! friction velocity [m/s]
        tstar,      &! temperature scaling parameter
        qstar,      &! moisture scaling parameter
        fm,         &! integral of profile function for momentum
        fh,         &! integral of profile function for heat
        fq           ! integral of profile function for moisture

!-----------------------Local Variables---------------------------------
! assign iteration parameters
   INTEGER, parameter :: itmax  = 40   !maximum number of iteration
   INTEGER, parameter :: itmin  = 6    !minimum number of iteration
   REAL(r8),parameter :: delmax = 3.0  !maximum change in leaf temperature [K]
   REAL(r8),parameter :: dtmin  = 0.01 !max limit for temperature convergence [K]
   REAL(r8),parameter :: dlemin = 0.1  !max limit for energy flux convergence [w/m2]

   REAL(r8) dtl(0:itmax+1)     !difference of tl between two iterative step

   REAL(r8) :: &
        displa,     &! displacement height [m]
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
        taf,        &! air temperature within canopy space [K]
        qaf,        &! humidity of canopy air [kg/kg]
        eah,        &! canopy air vapor pressure (pa)
        pco2g,      &! co2 pressure (pa) at ground surface (pa)
        pco2a,      &! canopy air co2 pressure (pa)

        fdry,       &! fraction of foliage that is green and dry [-]
        fwet,       &! fraction of foliage covered by water [-]
        cf,         &! heat transfer coefficient from leaves [-]
        rb,         &! leaf boundary layer resistance [s/m]
        rbone,      &! canopy bulk boundary layer resistance 
        rbsun,      &! canopy bulk boundary layer resistance 
        rbsha,      &! canopy bulk boundary layer resistance 
        rd,         &! aerodynamical resistance between ground and canopy air
        ram,        &! aerodynamical resistance [s/m]
        rah,        &! thermal resistance [s/m]
        raw,        &! moisture resistance [s/m]
        clai,       &! canopy heat capacity [Jm-2K-1]
        cah,        &! heat conductance for air [m/s]
        cgh,        &! heat conductance for ground [m/s]
        cfh,        &! heat conductance for leaf [m/s]
        caw,        &! latent heat conductance for air [m/s]
        cgw,        &! latent heat conductance for ground [m/s]
        cfw,        &! latent heat conductance for leaf [m/s]
        wtshi,      &! sensible heat resistance for air, grd and leaf [-]
        wtsqi,      &! latent heat resistance for air, grd and leaf [-]
        wta0,       &! normalized heat conductance for air [-]
        wtg0,       &! normalized heat conductance for ground [-]
        wtl0,       &! normalized heat conductance for air and leaf [-]
        wtaq0,      &! normalized latent heat conductance for air [-]
        wtgq0,      &! normalized heat conductance for ground [-]
        wtlq0,      &! normalized latent heat cond. for air and leaf [-]

        ei,         &! vapor pressure on leaf surface [pa]
        deidT,      &! derivative of "ei" on "tl" [pa/K]
        qsatl,      &! leaf specific humidity [kg/kg]
        qsatldT,    &! derivative of "qsatl" on "tlef"

        del,        &! absolute change in leaf temp in current iteration [K]
        dee,        &! maximum leaf temp. change in two consecutive iter [K]
        del2,       &! change in leaf temperature in previous iteration [K]
        dele,       &! change in heat fluxes from leaf [K]
        dele2,      &! change in heat fluxes from leaf [K]
        det,        &! maximum leaf temp. change in two consecutive iter [K]
 
        obuold,     &! monin-obukhov length from previous iteration
        tlbef,      &! leaf temperature from previous iteration [K]
        ecidif,     &! excess energies [W/m2]
        err,        &! balance error

        rssun,      &! sunlit leaf stomatal resistance [s/m]
        rssha,      &! shaded leaf stomatal resistance [s/m]
        fsha,       &! shaded fraction of canopy
        laisun,     &! sunlit leaf area index, one-sided
        laisha,     &! shaded leaf area index, one-sided
        assimsun,   &! sunlit leaf assimilation rate [umol co2 /m**2/ s] [+]
        assimsha,   &! shaded leaf assimilation rate [umol co2 /m**2/ s] [+]
        respcsun,   &! sunlit leaf respiration rate [umol co2 /m**2/ s] [+]
        respcsha,   &! shaded leaf respiration rate [umol co2 /m**2/ s] [+]
        rsoil,      &! soil respiration
        gah2o,      &! conductance between canopy and atmosphere
        gdh2o,      &! conductance between canopy and ground
        tprcor       ! tf*psur*100./1.013e5

   INTEGER it, nmozsgn 

   REAL(r8) delta, fac
   REAL(r8) evplwet, evplwet_dtl, etr_dtl, elwmax, elwdif
   REAL(r8) irab, dirab_dtl, fsenl_dtl, fevpl_dtl  
   REAL(r8) w, csoilcn, z0mg, cint(3), cintsun(3), cintsha(3)
   REAL(r8) fevpl_bef, fevpl_noadj, dtl_noadj, errt, erre

   REAL(r8) lt, egvf
   
   REAL(r8) :: sqrtdragc !sqrt(drag coefficient)
   REAL(r8) :: fai       !canopy frontal area index
   REAL(r8) :: a_k71     !exponential extinction factor for u/k decline within canopy (Kondo 1971)
   REAL(r8) :: fqt, fht, fmtop
   REAL(r8) :: utop, ueff, ktop
   REAL(r8) :: phih, z0qg, z0hg 
   REAL(r8) :: hsink, displasink

   INTEGER,  parameter :: zd_opt = 3   
   INTEGER,  parameter :: rb_opt = 3
   INTEGER,  parameter :: rd_opt = 3
 
!-----------------------End Variable List-------------------------------

! initialization of errors and  iteration parameters
       it     = 1    !counter for leaf temperature iteration
       del    = 0.0  !change in leaf temperature from previous iteration
       dele   = 0.0  !latent head flux from leaf for previous iteration

       dtl(0) = 0.
       fevpl_bef = 0.

       fht  = 0.     !integral of profile function for heat
       fqt  = 0.     !integral of profile function for moisture
          
!-----------------------------------------------------------------------
! scaling-up coefficients from leaf to canopy
!-----------------------------------------------------------------------

       fsha   = 1. - fsun
       laisun = lai*fsun
       laisha = lai*fsha

! scaling-up coefficients from leaf to canopy
       cintsun(1) = (1.-exp(-(0.110+extkb)*lai))/(0.110+extkb)
       cintsun(2) = (1.-exp(-(extkb+extkd)*lai))/(extkb+extkd)
       cintsun(3) = (1.-exp(-extkb*lai))/extkb

       cintsha(1) = (1.-exp(-0.110*lai))/0.110 - cintsun(1)
       cintsha(2) = (1.-exp(-extkd*lai))/extkd - cintsun(2)
       cintsha(3) = lai - cintsun(3)

!-----------------------------------------------------------------------
! get fraction of wet and dry canopy surface (fwet & fdry)
! initial saturated vapor pressure and humidity and their derivation
!-----------------------------------------------------------------------

       !clai = 4.2 * 1000. * 0.2
       clai = 0.0

       ! loop
       CALL dewfraction (sigf,lai,sai,dewmx,ldew,fwet,fdry)

       ! loop
       CALL qsadv(tl,psrf,ei,deiDT,qsatl,qsatlDT)

!-----------------------------------------------------------------------
! initial for fluxes profile
!-----------------------------------------------------------------------

       nmozsgn = 0    !number of times moz changes sign
       obuold = 0.    !monin-obukhov length from previous iteration
       zii = 1000.    !m (pbl height)
       beta = 1.      !- (in computing W_*)
       z0mg = (1.-fsno)*zlnd + fsno*zsno
       z0hg = z0mg
       z0qg = z0mg
      
! 12/27/2019, yuan: bug found
       ! initialize z0m for 1D case
       z0m    = htop * z0mr(patchclass(ipatch))
       displa = htop * displar(patchclass(ipatch))

       z0mv = z0m; z0hv = z0m; z0qv = z0m

       ! Modify aerodynamic parameters for sparse/dense canopy (X. Zeng)
       lt     = min(lai+sai, 2.)
       egvf   = (1._r8 - exp(-lt)) / (1._r8 - exp(-2.))
       displa = egvf * displa
       z0mv   = exp(egvf * log(z0mv) + (1._r8 - egvf) * log(z0mg))
       
       z0hv   = z0mv
       z0qv   = z0mv

       !print *, ipatch
       ! 10/17/2017: 3D z0m and displa
       IF (zd_opt == 3) THEN 

          CALL cal_z0_displa(lai+sai, htop, 1., z0mv, displa)
             
          ! NOTE: adjusted for samll displa
          displasink = max(htop/2., displa)
          hsink = z0mv + displasink
          
          z0hv   = z0mv
          z0qv   = z0mv
       ENDIF 
       
       fai    = 1. - exp(-0.5*(lai+sai))
       sqrtdragc = min( (0.003+0.3*fai)**0.5, 0.3 )
       
       a_k71 = htop/(htop-displa)/(vonkar/sqrtdragc)

       taf = 0.5 * (tg + thm)
       qaf = 0.5 * (qm + qg)

       pco2a = pco2m
       tprcor = 44.6*273.16*psrf/1.013e5
       rsoil = 0.  !respiration (mol m-2 s-1)
!      rsoil = 1.22e-6*exp(308.56*(1./56.02-1./(tg-227.13)))
!      rsoil = rstfac * 0.23 * 15. * 2.**((tg-273.16-10.)/10.) * 1.e-6
!      rsoil = 5.22 * 1.e-6
       rsoil = 0.22 * 1.e-6

       ur = max(0.1, sqrt(us*us+vs*vs))    !limit set to 0.1
       dth = thm - taf
       dqh = qm - qaf
       dthv = dth*(1.+0.61*qm) + 0.61*th*dqh
       zldis = hu - displa

       IF(zldis <= 0.0) THEN
          write(6,*) 'the obs height of u less than the zero displacement heght'
          CALL abort
       ENDIF

       CALL moninobukini(ur,th,thm,thv,dth,dqh,dthv,zldis,z0mv,um,obu)

! ======================================================================
!      BEGIN stability iteration 
! ======================================================================

       DO WHILE (it .le. itmax) 

          tlbef = tl

          del2  = del
          dele2 = dele

!-----------------------------------------------------------------------
! Aerodynamical resistances
!-----------------------------------------------------------------------
! Evaluate stability-dependent variables using moz from prior iteration
         IF (rd_opt == 3) THEN
           CALL moninobukm(hu,ht,hq,displa,z0mv,z0hv,z0qv,obu,um, &
               displasink,z0mv,ustar,fh2m,fq2m, &
               htop,fmtop,fm,fh,fq,fht,fqt,phih)
            ! Aerodynamic resistance
            ram = 1./(ustar*ustar/um)
            rah = 1./(vonkar/(fh-fht)*ustar) 
            raw = 1./(vonkar/(fq-fqt)*ustar) 
         ELSE 
            CALL moninobuk(hu,ht,hq,displa,z0mv,z0hv,z0qv,obu,um,&
               ustar,fh2m,fq2m,fm10m,fm,fh,fq)
            ! Aerodynamic resistance
            ram = 1./(ustar*ustar/um)
            rah = 1./(vonkar/fh*ustar) 
            raw = 1./(vonkar/fq*ustar) 
         ENDIF

         z0hg = z0mg/exp(0.13 * (ustar*z0mg/1.5e-5)**0.45)
         z0qg = z0hg

! Bulk boundary layer resistance of leaves
         uaf = ustar
         cf = 0.01*sqrtdi/sqrt(uaf)
         rb = 1/(cf*uaf)

        ! 11/17/2017: 3D rb calculation
! 03/13/2020, yuan: TODO, add analytical solution
         IF (rb_opt == 3) THEN 
            utop = ustar/vonkar * fmtop
            !REAL(r8) FUNCTION uintegral(utop, fc, bee, alpha, z0mg, htop, hbot)
            !ueff = uintegral(utop, 1._r8, 1._r8, a_k71, z0mg, htop, z0mg)
            !REAL(r8) FUNCTION ueffect(utop, htop, hbot, z0mg, alpha, bee, fc)
            ueff = ueffect(utop, htop, z0mg, z0mg, a_k71, 1._r8, 1._r8)
            cf = 0.01*sqrtdi*sqrt(ueff)
            rb = 1./cf
         ENDIF
         
       ! rd = 1./(csoilc*uaf)                 ! BATS legacy
       ! w = exp(-0.5*(lai+sai))              ! Dickinson's modification :
       ! csoilc = ( 1.-w + w*um/uaf)/rah      ! "rah" here is the resistance over
       ! rd = 1./(csoilc*uaf)                 ! bare ground fraction

! modified by Xubin Zeng's suggestion at 08-07-2002
         w = exp(-(lai+sai))
         csoilcn = (vonkar/(0.13*(z0mg*uaf/1.5e-5)**0.45))*w + csoilc*(1.-w)
         rd = 1./(csoilcn*uaf)

         ! 11/17/2017: 3D rd calculation
! 03/13/2020, yuan: TODO, add analytical solution
         IF (rd_opt == 3) THEN 
            ktop = vonkar * (htop-displa) * ustar / phih
            !REAL(r8) FUNCTION kintegral(ktop, fc, bee, alpha, z0mg, &
            !  displah, htop, hbot, obu, ustar, ztop, zbot)
            !rd = kintegral(ktop, 1._r8, 1._r8, a_k71, z0mg, &
            !   displa/htop, htop, z0qg, obug, ustarg, hsink, z0qg )
 
            !REAL(r8) FUNCTION frd(ktop, htop, hbot, ztop, zbot, displah, &
            !  z0h, obu, ustar, z0mg, alpha, bee, fc)
            rd = frd(ktop, htop, z0qg, hsink, z0qg, displa/htop, &
               z0qg, obug, ustar, z0mg, a_k71, 1._r8, 1._r8)
            !print *, "ktop:", ktop, htop,z0mv,ustar,phih !fordebug
         ENDIF

!-----------------------------------------------------------------------
! stomatal resistances 
!-----------------------------------------------------------------------

         IF(lai .gt. 0.001) THEN
            
            rbsun = rb / laisun
            rbsha = rb / laisha
 
            eah = qaf * psrf / ( 0.622 + 0.378 * qaf )    !pa

! Sunlit leaves
            CALL stomata  (vmax25   ,effcon ,slti   ,hlti    ,&
                 shti     ,hhti     ,trda   ,trdm   ,trop    ,&
                 gradm    ,binter   ,thm    ,psrf   ,po2m    ,&
                 pco2m    ,pco2a    ,eah    ,ei     ,tl      ,&
                 parsun   ,rbsun    ,raw    ,rstfac ,cintsun ,&
                 assimsun ,respcsun ,rssun  )

! Shaded leaves
            CALL stomata  (vmax25   ,effcon ,slti   ,hlti    ,&
                 shti     ,hhti     ,trda   ,trdm   ,trop    ,&
                 gradm    ,binter   ,thm    ,psrf   ,po2m    ,&
                 pco2m    ,pco2a    ,eah    ,ei     ,tl      ,&
                 parsha   ,rbsha    ,raw    ,rstfac ,cintsha ,&
                 assimsha ,respcsha ,rssha  )

         ELSE
            rssun = 2.e4; assimsun = 0.; respcsun = 0.
            rssha = 2.e4; assimsha = 0.; respcsha = 0.
         ENDIF

! above stomatal resistances are for the canopy, the stomatal rsistances 
! and the "rb" in the following calculations are the average for single leaf. thus,
         rssun = rssun * laisun
         rssha = rssha * laisha

!-----------------------------------------------------------------------
! dimensional and non-dimensional sensible and latent heat conductances
! for canopy and soil flux calculations.
!-----------------------------------------------------------------------

         delta = 0.0
         IF(qsatl-qaf .gt. 0.) delta = 1.0

         cah = 1. / rah
         cgh = 1. / rd
         cfh = (lai + sai) / rb

         caw = 1. / raw
         cgw = 1. / rd
         cfw = (1.-delta*(1.-fwet))*(lai+sai)/rb + (1.-fwet)*delta* &
            ( laisun/(rb+rssun) + laisha/(rb+rssha) ) 

         wtshi = 1. / ( cah + cgh + cfh )
         wtsqi = 1. / ( caw + cgw + cfw )

         wta0 = cah * wtshi
         wtg0 = cgh * wtshi
         wtl0 = cfh * wtshi

         wtaq0 = caw * wtsqi
         wtgq0 = cgw * wtsqi
         wtlq0 = cfw * wtsqi
 
!-----------------------------------------------------------------------
! IR radiation, sensible and latent heat fluxes and their derivatives
!-----------------------------------------------------------------------
! the partial derivatives of areodynamical resistance are ignored 
! which cannot be determined analtically
         fac = 1. - thermk

! longwave absorption and their derivatives 
         irab = (frl - 2. * stefnc * tl**4 + emg*stefnc*tg**4 ) * fac &
            + (1-emg)*thermk*fac*frl + (1-emg)*(1-thermk)*fac*stefnc*tl**4
         dirab_dtl = - 8. * stefnc * tl**3                      * fac &
            + 4.*(1-emg)*(1-thermk)*fac*stefnc*tl**3

! sensible heat fluxes and their derivatives
         fsenl = rhoair * cpair * cfh * ( (wta0 + wtg0)*tl - wta0*thm - wtg0*tg )
         fsenl_dtl = rhoair * cpair * cfh * (wta0 + wtg0)

! latent heat fluxes and their derivatives
         etr = rhoair * (1.-fwet) * delta &
             * ( laisun/(rb+rssun) + laisha/(rb+rssha) ) &
             * ( (wtaq0 + wtgq0)*qsatl - wtaq0*qm - wtgq0*qg )
         etr_dtl = rhoair * (1.-fwet) * delta &
             * ( laisun/(rb+rssun) + laisha/(rb+rssha) ) &
             * (wtaq0 + wtgq0)*qsatlDT 
      
         IF(etr.ge.etrc)THEN
            etr = etrc
            etr_dtl = 0.
         ENDIF
 
         evplwet = rhoair * (1.-delta*(1.-fwet)) * (lai+sai) / rb &
                 * ( (wtaq0 + wtgq0)*qsatl - wtaq0*qm - wtgq0*qg )
         evplwet_dtl = rhoair * (1.-delta*(1.-fwet)) * (lai+sai) / rb &
                     * (wtaq0 + wtgq0)*qsatlDT 
         IF(evplwet.ge.ldew/deltim)THEN
            evplwet = ldew/deltim
            evplwet_dtl = 0.
         ENDIF
  
         fevpl = etr + evplwet
         fevpl_dtl = etr_dtl + evplwet_dtl
 
         erre = 0.
         fevpl_noadj = fevpl
         IF ( fevpl*fevpl_bef < 0. ) THEN 
            erre  = -0.9*fevpl
            fevpl =  0.1*fevpl
         ENDIF

!-----------------------------------------------------------------------
! difference of temperatures by quasi-newton-raphson method for the non-linear system equations
!-----------------------------------------------------------------------

         dtl(it) = (sabv + irab - fsenl - hvap*fevpl) &
             / ((lai+sai)*clai/deltim - dirab_dtl + fsenl_dtl + hvap*fevpl_dtl)
         dtl_noadj = dtl(it)
 
         ! check magnitude of change in leaf temperature limit to maximum allowed value
 
         IF(it .le. itmax) THEN
 
         ! put brakes on large temperature excursions
           IF(abs(dtl(it)).gt.delmax)THEN
               dtl(it) = delmax*dtl(it)/abs(dtl(it))
           ENDIF
 
           IF((it.ge.2) .and. (dtl(it-1)*dtl(it).le.0.))THEN
               dtl(it) = 0.5*(dtl(it-1) + dtl(it))
           ENDIF
 
         ENDIF
 
         tl = tlbef + dtl(it)

!-----------------------------------------------------------------------
! square roots differences of temperatures and fluxes for use as the condition of convergences
!-----------------------------------------------------------------------

         del  = sqrt( dtl(it)*dtl(it) )
         dele = dtl(it) * dtl(it) * ( dirab_dtl**2 + fsenl_dtl**2 + hvap*fevpl_dtl**2 ) 
         dele = sqrt(dele)

!-----------------------------------------------------------------------
!  saturated vapor pressures and canopy air temperature, canopy air humidity
!-----------------------------------------------------------------------
! Recalculate leaf saturated vapor pressure (ei_)for updated leaf temperature
! and adjust specific humidity (qsatl_) proportionately
         CALL qsadv(tl,psrf,ei,deiDT,qsatl,qsatlDT)

! update vegetation/ground surface temperature, canopy air temperature, 
! canopy air humidity
         taf = wta0*thm + wtg0*tg + wtl0*tl 
 
         qaf = wtaq0*qm + wtgq0*qg + wtlq0*qsatl
 
! update co2 partial pressure within canopy air
         gah2o = 1.0/raw * tprcor/thm                     !mol m-2 s-1
         gdh2o = 1.0/rd  * tprcor/thm                     !mol m-2 s-1
         pco2a = pco2m - 1.37*psrf/max(0.446,gah2o) * &
            (assimsun + assimsha  - respcsun -respcsha - rsoil)

!-----------------------------------------------------------------------
! Update monin-obukhov length and wind speed including the stability effect
!-----------------------------------------------------------------------

         dth = thm - taf       
         dqh = qm - qaf
 
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
            det = max(del,del2)
            dee = max(dele,dele2)
            IF(det .lt. dtmin .and. dee .lt. dlemin) exit 
         ENDIF
 
       ENDDO 
        
       !IF (it > itmax) print *, "*** NOTE: it = 41! ***"

! ======================================================================
!      END stability iteration 
! ======================================================================

       z0m = z0mv
       zol = zeta
       rib = min(5.,zol*ustar**2/(vonkar**2/fh*um**2))

! canopy fluxes and total assimilation amd respiration

       IF(lai .gt. 0.001) THEN
          rst = 1./(laisun/rssun + laisha/rssha)
       ELSE
          rssun = 2.0e4 ; rssha = 2.0e4
          assimsun = 0. ; assimsha = 0.
          respcsun = 0. ; respcsha = 0.
          rst = 2.0e4
       ENDIF
       assim = assimsun + assimsha
       respc = respcsun + respcsha + rsoil

! canopy fluxes and total assimilation amd respiration
       fsenl = fsenl + fsenl_dtl*dtl(it-1) &
               ! add the imbalanced energy below due to T adjustment to sensibel heat
               + (dtl_noadj-dtl(it-1)) * ((lai+sai)*clai/deltim - dirab_dtl + fsenl_dtl + hvap*fevpl_dtl) &
               ! add the imbalanced energy below due to q adjustment to sensibel heat
               + hvap*erre
 
       etr     = etr     +     etr_dtl*dtl(it-1)
       evplwet = evplwet + evplwet_dtl*dtl(it-1)
       fevpl   = fevpl_noadj
       fevpl   = fevpl   +   fevpl_dtl*dtl(it-1)
 
       elwmax  = ldew/deltim
       elwdif  = max(0., evplwet-elwmax)
       evplwet = min(evplwet, elwmax)
 
       fevpl = fevpl - elwdif
       fsenl = fsenl + hvap*elwdif 
 
       taux = - rhoair*us/ram
       tauy = - rhoair*vs/ram

!-----------------------------------------------------------------------
! fluxes from ground to canopy space
!-----------------------------------------------------------------------

       fseng = cpair*rhoair*cgh*(tg-taf)
       fevpg = rhoair*cgw*(qg-qaf)

!-----------------------------------------------------------------------
! downward (upward) longwave radiation below (above) the canopy
!-----------------------------------------------------------------------

       ! 10/16/2017: add soil reflectance
       ! ONLY for vegetation covered area
       dlrad = thermk * frl &
             + stefnc * fac * tlbef**3 * (tlbef + 4.*dtl(it-1)) 
       ulrad = stefnc * ( fac * tlbef**3 * (tlbef + 4.*dtl(it-1)) &
             + thermk*emg*tg**4 ) &
             + (1-emg)*thermk*thermk*frl &
             + (1-emg)*thermk*fac*stefnc*tlbef**4 &
             + 4.*(1-emg)*thermk*fac*stefnc*tlbef**3*dtl(it-1)

!-----------------------------------------------------------------------
! Derivative of soil energy flux with respect to soil temperature (cgrnd)
!-----------------------------------------------------------------------

       cgrnds = cpair*rhoair*cgh*(1.-wtg0)
       cgrndl = rhoair*cgw*(1.-wtgq0)*dqgdT
       cgrnd  = cgrnds + cgrndl*htvp

!-----------------------------------------------------------------------
! balance check
! (the computational error was created by the assumed 'dtl' in line 406-408) 
!-----------------------------------------------------------------------

       err = sabv + irab + dirab_dtl*dtl(it-1) - fsenl - hvap*fevpl

#if(defined CLMDEBUG)
       IF(abs(err) .gt. .2) &
       write(6,*) 'energy imbalance in leaftemone.F90',it-1,err,sabv,irab,fsenl,hvap*fevpl
#endif

!-----------------------------------------------------------------------
! Update dew accumulation (kg/m2)
!-----------------------------------------------------------------------

       ldew = max(0., ldew-evplwet*deltim)

!-----------------------------------------------------------------------
! 2 m height air temperature
!-----------------------------------------------------------------------

       tref = thm + vonkar/(fh-fht)*dth * (fh2m/vonkar - fh/vonkar) 
       qref =  qm + vonkar/(fq-fqt)*dqh * (fq2m/vonkar - fq/vonkar)

  END SUBROUTINE LeafTemp
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
     !dz = 0.001 !fordebug
     dz = 0.05 !fordebug
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
     !dz = 0.001 !fordebug
     dz = 0.05 !fordebug
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

END MODULE LEAF_temperature
