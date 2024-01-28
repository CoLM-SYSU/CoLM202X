#include <define.h>

MODULE MOD_Urban_Flux

!-----------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Vars_Global
   USE MOD_Qsadv, only: qsadv
   IMPLICIT NONE
   SAVE

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: UrbanOnlyFlux
   PUBLIC :: UrbanVegFlux
   PUBLIC :: dewfraction

! Exponential extinction factor (alpha) options:
!   1. Masson, 2000; Oleson et al., 2008
!   2. Swaid, 1993; Kusaka, 2001; Lee and Park, 2008
!   3. Macdonald, 2000
   integer, parameter :: alpha_opt = 3

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE UrbanOnlyFlux ( &
         ! Model running information
         ipatch      ,deltim      ,lbr         ,lbi         ,&
         ! Forcing
         hu          ,ht          ,hq          ,us          ,&
         vs          ,thm         ,th          ,thv         ,&
         qm          ,psrf        ,rhoair      ,Fhac        ,&
         Fwst        ,Fach        ,vehc        ,meta        ,&
         ! Urban parameters
         hroof       ,hwr         ,nurb        ,fcover      ,&
         ! Status of surface
         z0h_g       ,obug        ,ustarg      ,zlnd        ,&
         zsno        ,fsno_roof   ,fsno_gimp   ,fsno_gper   ,&
         wliq_roofsno,wliq_gimpsno,wice_roofsno,wice_gimpsno,&
         htvp_roof   ,htvp_gimp   ,htvp_gper   ,troof       ,&
         twsun       ,twsha       ,tgimp       ,tgper       ,&
         qroof       ,qgimp       ,qgper       ,dqroofdT    ,&
         dqgimpdT    ,dqgperdT    ,rsr                      ,&
         ! Output
         taux        ,tauy        ,fsenroof    ,fsenwsun    ,&
         fsenwsha    ,fsengimp    ,fsengper    ,fevproof    ,&
         fevpgimp    ,fevpgper    ,croofs      ,cwalls      ,&
         cgrnds      ,croofl      ,cgimpl      ,cgperl      ,&
         croof       ,cgimp       ,cgper       ,tref        ,&
         qref        ,z0m         ,zol         ,rib         ,&
         ustar       ,qstar       ,tstar       ,fm          ,&
         fh          ,fq          ,tafu                      )

!=======================================================================
   USE MOD_Precision
   USE MOD_Const_Physical, only: cpair,vonkar,grav
   USE MOD_FrictionVelocity
   USE MOD_CanopyLayerProfile
   IMPLICIT NONE

!----------------------- Dummy argument --------------------------------
   integer, intent(in) :: &
        ipatch,   &! patch index [-]
        lbr,      &! lower bound of array
        lbi        ! lower bound of array

   real(r8), intent(in) :: &
        deltim     ! seconds in a time step [second]

   ! atmospherical variables and observational height
   real(r8), intent(in) :: &
        hu,       &! observational height of wind [m]
        ht,       &! observational height of temperature [m]
        hq,       &! observational height of humidity [m]
        us,       &! wind component in eastward direction [m/s]
        vs,       &! wind component in northward direction [m/s]
        thm,      &! intermediate variable (tm+0.0098*ht) [K]
        th,       &! potential temperature (kelvin)
        thv,      &! virtual potential temperature (kelvin)
        qm,       &! specific humidity at agcm reference height [kg/kg]
        psrf,     &! atmosphere pressure at the surface [pa] [not used]
        rhoair     ! density air [kg/m3]

   real(r8), intent(in) :: &
        vehc,     &! flux from vehicle
        meta,     &! flux from metabolic
        Fhac,     &! flux from heat or cool AC
        Fwst,     &! waste heat from cool or heat
        Fach       ! flux from air exchange

   integer, intent(in) :: &
        nurb       ! number of aboveground urban components [-]

   real(r8), intent(in) :: &
        hroof,    &! average building height [m]
        hwr,      &! average building height to their distance [-]
        fcover(0:4)! coverage of aboveground urban components [-]

   real(r8), intent(in) :: &
        rsr,      &! bare soil resistance for evaporation
        z0h_g,    &! roughness length for bare ground, sensible heat [m]
        obug,     &! monin-obukhov length for bare ground (m)
        ustarg,   &! friction velocity for bare ground [m/s]
        zlnd,     &! roughness length for soil [m]
        zsno,     &! roughness length for snow [m]
        fsno_roof,&! fraction of ground covered by snow [-]
        fsno_gimp,&! fraction of ground covered by snow [-]
        fsno_gper,&! fraction of ground covered by snow [-]
        wliq_roofsno,&! liqui water [kg/m2]
        wliq_gimpsno,&! liqui water [kg/m2]
        wice_roofsno,&! ice lens [kg/m2]
        wice_gimpsno,&! ice lens [kg/m2]
        htvp_roof,&! latent heat of vapor of water (or sublimation) [j/kg]
        htvp_gimp,&! latent heat of vapor of water (or sublimation) [j/kg]
        htvp_gper,&! latent heat of vapor of water (or sublimation) [j/kg]

        troof,    &! temperature of roof [K]
        twsun,    &! temperature of sunlit wall [K]
        twsha,    &! temperature of shaded wall [K]
        tgimp,    &! temperature of impervious road [K]
        tgper,    &! pervious ground temperature [K]

        qroof,    &! roof specific humidity [kg/kg]
        qgimp,    &! imperivous road specific humidity [kg/kg]
        qgper,    &! pervious ground specific humidity [kg/kg]
        dqroofdT, &! d(qroof)/dT
        dqgimpdT, &! d(qgimp)/dT
        dqgperdT   ! d(qgper)/dT

   ! Output
   real(r8), intent(out) :: &
        taux,     &! wind stress: E-W [kg/m/s**2]
        tauy,     &! wind stress: N-S [kg/m/s**2]
        fsenroof, &! sensible heat flux from roof [W/m2]
        fsenwsun, &! sensible heat flux from snulit wall [W/m2]
        fsenwsha, &! sensible heat flux from shaded wall [W/m2]
        fsengimp, &! sensible heat flux from impervious road [W/m2]
        fsengper, &! sensible heat flux from pervious ground [W/m2]
        fevproof, &! evaperation heat flux from roof [W/m2]
        fevpgimp, &! evaperation heat flux from impervious road [W/m2]
        fevpgper, &! evaporation heat flux from pervious ground [mm/s]

        croofs,   &! deriv of roof sensible heat flux wrt soil temp [w/m**2/k]
        cwalls,   &! deriv of wall sensible heat flux wrt soil temp [w/m**2/k]
        cgrnds,   &! deriv of soil sensible heat flux wrt soil temp [w/m**2/k]
        croofl,   &! deriv of roof latent heat flux wrt soil temp [w/m**2/k]
        cgimpl,   &! deriv of gimp latent heat flux wrt soil temp [w/m**2/k]
        cgperl,   &! deriv of soil latent heat flux wrt soil temp [w/m**2/k]
        croof,    &! deriv of roof total heat flux wrt soil temp [w/m**2/k]
        cgimp,    &! deriv of gimp total heat flux wrt soil temp [w/m**2/k]
        cgper,    &! deriv of soil total heat flux wrt soil temp [w/m**2/k]

        tref,     &! 2 m height air temperature [kelvin]
        qref,     &! 2 m height air humidity [kg/kg]

        z0m,      &! effective roughness [m]
        zol,      &! dimensionless height (z/L) used in Monin-Obukhov theory
        rib,      &! bulk Richardson number in surface layer
        ustar,    &! friction velocity [m/s]
        tstar,    &! temperature scaling parameter
        qstar,    &! moisture scaling parameter
        fm,       &! integral of profile function for momentum
        fh,       &! integral of profile function for heat
        fq,       &! integral of profile function for moisture
        tafu       ! effective urban air temperature (2nd layer, walls)

!------------------------ LOCAL VARIABLES ------------------------------
   integer ::   &
        niters,   &! maximum number of iterations for surface temperature
        iter,     &! iteration index
        nmozsgn    ! number of times moz changes sign

   real(r8) ::  &
        beta,     &! coefficient of conective velocity [-]
        dth,      &! diff of virtual temp. between ref. height and surface
        dqh,      &! diff of humidity between ref. height and surface
        dthv,     &! diff of vir. poten. temp. between ref. height and surface
        obu,      &! monin-obukhov length (m)
        obuold,   &! monin-obukhov length from previous iteration
        ram,      &! aerodynamical resistance [s/m]
        rah,      &! thermal resistance [s/m]
        raw,      &! moisture resistance [s/m]
        fh2m,     &! relation for temperature at 2m
        fq2m,     &! relation for specific humidity at 2m
        fm10m,    &! integral of profile function for momentum at 10m
        thvstar,  &! virtual potential temperature scaling parameter
        um,       &! wind speed including the stablity effect [m/s]
        ur,       &! wind speed at reference height [m/s]
        wc,       &! convective velocity [m/s]
        wc2,      &! wc**2
        zeta,     &! dimensionless height used in Monin-Obukhov theory
        zii,      &! convective boundary height [m]
        zldis,    &! reference height "minus" zero displacement heght [m]
        z0mg,     &! roughness length over ground, momentum [m]
        z0hg,     &! roughness length over ground, sensible heat [m]
        z0qg       ! roughness length over ground, latent heat [m]

   real(r8) evplwet, evplwet_dtl, elwmax, elwdif

!----------------------- defination for 3d run ------------------------ !

   integer, parameter :: nlay = 3  ! potential layer number

   integer :: &
        clev,     &! current layer index
        numlay     ! available layer number

   real(r8) :: &
        huu,      &! observational height of wind [m]
        htu,      &! observational height of temperature [m]
        hqu,      &! observational height of humidity [m]
        ktop,     &! K value at a specific height
        utop,     &! u value at a specific height
        fht,      &! integral of profile function for heat at the top layer
        fqt,      &! integral of profile function for moisture at the top layer
        fmtop,    &! fm value at a specific height
        phih,     &! phi(h), similarity function for sensible heat
        displa,   &! displacement height for urban
        displau,  &! displacement height for urban building
        z0mu,     &! roughless length for urban building only
        z0h,      &! roughless length for sensible heat
        z0q,      &! roughless length for latent heat
        tg,       &! ground temperature
        qg         ! ground specific humidity

   real(r8) :: &
        fg,       &! ground fractional cover
        fgimp,    &! weight of impervious ground
        fgper,    &! weight of pervious ground
        hlr,      &! average building height to their length of edge [-]
        sqrtdragc,&! sqrt(drag coefficient)
        lm,       &! mix length within canopy
        fai,      &! frontal area index
        fwet,     &! fractional wet area
        delta,    &! 0 or 1
        alpha      ! exponential extinction factor for u/k decline within urban

   real(r8), dimension(0:nurb) :: &
        tu,       &! termperature array
        fc,       &! fractional cover array
        canlev,   &! urban canopy layer lookup table
        rb,       &! leaf boundary layer resistance [s/m]
        cfh,      &! heat conductance for leaf [m/s]
        cfw,      &! latent heat conductance for leaf [m/s]
        wtl0,     &! normalized heat conductance for air and leaf [-]
        wtlq0,    &! normalized latent heat cond. for air and leaf [-]

        ei,       &! vapor pressure on leaf surface [pa]
        deidT,    &! derivative of "ei" on "tl" [pa/K]
        qsatl,    &! leaf specific humidity [kg/kg]
        qsatldT    ! derivative of "qsatl" on "tlef"

   real(r8), dimension(nlay) :: &
        fah,      &! weight for thermal resistance to upper layer
        faw,      &! weight for moisture resistance to upper layer
        fgh,      &! weight for thermal resistance to lower layer
        fgw,      &! weight for moisture resistance to lower layer
        ueff_lay, &! effective wind speed within canopy layer [m/s]
        ueff_lay_,&! effective wind speed within canopy layer [m/s]
        taf,      &! air temperature within canopy space [K]
        qaf,      &! humidity of canopy air [kg/kg]
        rd,       &! aerodynamic resistance between layers [s/m]
        rd_,      &! aerodynamic resistance between layers [s/m]
        cah,      &! heat conductance for air [m/s]
        cgh,      &! heat conductance for ground [m/s]
        caw,      &! latent heat conductance for air [m/s]
        cgw,      &! latent heat conductance for ground [m/s]
        wtshi,    &! sensible heat resistance for air, grd and leaf [-]
        wtsqi,    &! latent heat resistance for air, grd and leaf [-]
        wta0,     &! normalized heat conductance for air [-]
        wtg0,     &! normalized heat conductance for ground [-]
        wtaq0,    &! normalized latent heat conductance for air [-]
        wtgq0,    &! normalized heat conductance for ground [-]
        wtll,     &! sum of normalized heat conductance for air and leaf
        wtlql      ! sum of normalized heat conductance for air and leaf

   real(r8) ::  &
        ra2m,     &! aerodynamic resistance between 2m and bottom layer [s/m]
        rd2m       ! aerodynamic resistance between bottom layer and ground [s/m]

   ! temporal
   integer i
   real(r8) h_vec, l_vec, tmpw3, cgw_per, cgw_imp
   real(r8) bee, tmpw1, tmpw2, fact, facq
   real(r8) fwet_roof, fwet_roof_, fwet_gimp, fwet_gimp_
   real(r8) fwetfac

!-----------------------End Variable List-------------------------------

! initialization
      tu(0) = troof; tu(1) = twsun; tu(2) = twsha

      fc(:)  = fcover(0:nurb)
      fg     = 1 - fcover(0)
      fgimp  = fcover(3)/fg
      fgper  = fcover(4)/fg
      hlr    = hwr*(1-sqrt(fcover(0)))/sqrt(fcover(0))
      canlev = (/3, 2, 2/)
      numlay = 2

!-----------------------------------------------------------------------
! initial roughness length for z0mg, z0hg, z0qg
! Roughness of the city ground only (excluding buildings and vegetation)

      !NOTE: change to original
      !z0mg = (1.-fsno)*zlnd + fsno*zsno
      IF (fsno_gper > 0) THEN
         z0mg = zsno
      ELSE
         z0mg = zlnd
      ENDIF
      z0hg = z0mg
      z0qg = z0mg

!-----------------------------------------------------------------------
! initial saturated vapor pressure and humidity and their derivation
!    0: roof, 1: sunlit wall, 2: shaded wall
!-----------------------------------------------------------------------

      qsatl(0) = qroof
      qsatldT(0) = dqroofdT
      DO i = 1, nurb
         CALL qsadv(tu(i),psrf,ei(i),deiDT(i),qsatl(i),qsatldT(i))
      ENDDO

!-----------------------------------------------------------------------
! set weight
!-----------------------------------------------------------------------

      ! set weighting factor
      fah(1) = 1.; fah(2) = 1.; fah(3) = 1.
      faw(1) = 1.; faw(2) = 1.; faw(3) = 1.
      fgh(1) = 1.; fgh(2) = fg; fgh(3) = 1.
      fgw(1) = 1.; fgw(2) = fg; fgw(3) = 1.

      ! weighted tg
      tg = tgimp*fgimp + tgper*fgper

      ! wet fraction for roof and impervious ground
      !-------------------------------------------
      ! roof
      IF (lbr < 1) THEN
         fwet_roof_ = fsno_roof !for snow layer exist
      ELSE
         ! surface wet fraction. assuming max ponding = 1 kg/m2
         fwet_roof_ = (max(0., wliq_roofsno+wice_roofsno))**(2/3.)
         fwet_roof_ = min(1., fwet_roof_)
      ENDIF

      ! impervious ground
      IF (lbi < 1) THEN
         fwet_gimp_ = fsno_gimp !for snow layer exist
      ELSE
         ! surface wet fraction. assuming max ponding = 1 kg/m2
         fwet_gimp_ = (max(0., wliq_gimpsno+wice_gimpsno))**(2/3.)
         fwet_gimp_ = min(1., fwet_gimp_)
      ENDIF

      ! dew case
      IF (qm > qroof) THEN
         fwet_roof = 1.
      ELSE
         fwet_roof = fwet_roof_
      ENDIF

      ! ! dew case
      IF (qm > qgimp) THEN
         fwet_gimp = 1.
      ELSE
         fwet_gimp = fwet_gimp_
      ENDIF

      ! weighted qg
      ! NOTE: IF fwet_gimp=1, same as previous
      fwetfac = fgimp*fwet_gimp + fgper
      qg = (qgimp*fgimp*fwet_gimp + qgper*fgper) / fwetfac

      fgw(2) = fg*fwetfac

!-----------------------------------------------------------------------
! initial for fluxes profile
!-----------------------------------------------------------------------

      nmozsgn = 0     !number of times moz changes sign
      obuold  = 0.    !monin-obukhov length from previous iteration
      zii     = 1000. !m (pbl height)
      beta    = 1.    !- (in computing W_*)

!-----------------------------------------------------------------------
! scaling factor bee
!-----------------------------------------------------------------------
!NOTE: bee value, the default is 1
      bee = 1.

!-----------------------------------------------------------------------
! calculate z0m and displa
!-----------------------------------------------------------------------

      ! Macdonald et al., 1998, Eq. (23), A=4.43
      displau = hroof * (1 + 4.43**(-fcover(0))*(fcover(0) - 1))
      fai  = 4/PI*hlr*fcover(0)
      z0mu = (hroof - displau) * &
           exp( -(0.5*1.2/vonkar/vonkar*(1-displau/hroof)*fai)**(-0.5) )

      ! to compare z0 of urban and only the surface
      ! maximum assumption
      IF (z0mu < z0mg) z0mu = z0mg

      ! roughness length and displacement height for sensible
      ! and latent heat transfer
      z0m = z0mu

      displa  = displau
      displau = max(hroof/2., displau)

!-----------------------------------------------------------------------
! calculate layer decay coefficient
!-----------------------------------------------------------------------

      !NOTE: the below is for vegetation, may not suitable for urban
      ! Raupach, 1992
      !sqrtdragc = min( (0.003+0.3*fai)**0.5, 0.3 )

      ! Kondo, 1971
      !alpha = hroof/(hroof-displa)/(vonkar/sqrtdragc)

      ! Masson, 2000; Oleson et al., 2008
      IF (alpha_opt == 1) alpha = 0.5*hwr

      ! Swaid, 1993; Kusaka, 2001; Lee and Park, 2008
      IF (alpha_opt == 2) alpha = 0.772*hwr

      ! Macdonald, 2000
      IF (alpha_opt == 3) alpha = 9.6*fai

!-----------------------------------------------------------------------
! first guess for taf and qaf for each layer
! a large differece from previous schemes
!-----------------------------------------------------------------------

      IF (numlay .eq. 2) THEN
         taf(3) = (tg + 2.*thm)/3.
         qaf(3) = (qg + 2.*qm )/3.
         taf(2) = (2.*tg + thm)/3.
         qaf(2) = (2.*qg + qm )/3.
      ENDIF

! initialization and input values for Monin-Obukhov
      ! have been set before
      z0h  = z0m; z0q = z0m
      ur   = max(0.1, sqrt(us*us+vs*vs))    !limit set to 0.1
      dth  = thm - taf(2)
      dqh  =  qm - qaf(2)
      dthv = dth*(1.+0.61*qm) + 0.61*th*dqh

      ! to ensure the obs height >= hroof+10.
      huu = max(hroof+10., hu)
      htu = max(hroof+10., ht)
      hqu = max(hroof+10., hq)

      zldis = huu - displa

      IF (zldis <= 0.0) THEN
         write(6,*) 'the obs height of u less than the zero displacement heght'
         CALL abort
      ENDIF

      CALL moninobukini(ur,th,thm,thv,dth,dqh,dthv,zldis,z0m,um,obu)

      niters=6

! ======================================================================
!     BEGIN stability iteration
! ======================================================================

      ITERATION : DO iter = 1, niters         !begin stability iteration

!-----------------------------------------------------------------------
! Aerodynamical resistances
!-----------------------------------------------------------------------
! Evaluate stability-dependent variables using moz from prior iteration

         !NOTE: displat=hroof, z0mt=0, are set for roof
         ! fmtop is calculated at the same height of fht, fqt
         CALL moninobukm(huu,htu,hqu,displa,z0m,z0h,z0q,obu,um, &
            hroof,0.,ustar,fh2m,fq2m,hroof,fmtop,fm,fh,fq,fht,fqt,phih)

! Aerodynamic resistance
         ! 09/16/2017:
         ! note that for ram, it is the resistance from Href to z0mv+displa
         ! however, for rah and raw is only from Href to canopy effective
         ! exchange height.
         ! for Urban: from Href to roof height
         ! so rah/raw is not comparable with that of 1D case
         ram = 1./(ustar*ustar/um)

         ! 05/02/2016: calculate resistance from the top layer (effective exchange
         ! height) to reference height
         ! for Urban: from roof height to reference height
         rah = 1./(vonkar/(fh-fht)*ustar)
         raw = 1./(vonkar/(fq-fqt)*ustar)

         ! update roughness length for sensible/latent heat
         z0hg = z0mg/exp(0.13 * (ustar*z0mg/1.5e-5)**0.45)
         z0qg = z0hg

         z0h  = max(z0hg, z0h)
         z0q  = max(z0qg, z0q)

!-----------------------------------------------------------------------
! new method to calculate rd and ueffect
! the kernel part of 3d model
!-----------------------------------------------------------------------

         ! initialization
         rd(:)  = 0.
         rd_(:) = 0.
         ueff_lay(:)  = 0.
         ueff_lay_(:) = 0.

         ! calculate canopy top wind speed (utop) and exchange coefficient (ktop)
         ! need to update each time as obu changed after each iteration
         utop = ustar/vonkar * fmtop
         ktop = vonkar * (hroof-displa) * ustar / phih

         ueff_lay(3) = utop

         !real(r8) FUNCTION kintegral(ktop, fc, bee, alpha, z0mg, &
         !      displah, htop, hbot, obu, ustar, ztop, zbot)
         !rd(3)  = kintegral(ktop, 1., bee, alpha, z0mg, displa/hroof, &
         !   hroof, 0., obug, ustarg, hroof, displa+z0m)

         !real(r8) FUNCTION frd(ktop, htop, hbot, &
         !      ztop, zbot, displah, z0h, obu, ustar, &
         !      z0mg, alpha, bee, fc)
         rd(3) = frd(ktop, hroof, 0., hroof, displau+z0mu, displa/hroof, z0h_g, &
            obug, ustarg, z0mg, alpha, bee, 1.)

         !real(r8) FUNCTION uintegralz(utop, fc, bee, alpha, z0mg, htop, hbot, ztop, zbot)
         !ueff_lay(2)  = uintegralz(utop, 1., bee, alpha, z0mg, hroof, 0., hroof, z0mg)

         !real(r8) FUNCTION ueffectz(utop, htop, hbot, ztop, zbot, z0mg, alpha, bee, fc)
         ueff_lay(2) = ueffectz(utop, hroof, 0., hroof, z0mg, z0mg, alpha, bee, 1.)

         !rd(2)  = kintegral(ktop, 1., bee, alpha, z0mg, displa/hroof, &
         !   hroof, 0., obug, ustarg, displau+z0mu, z0qg)
         rd(2) = frd(ktop, hroof, 0., displau+z0mu, z0qg, displa/hroof, z0h_g, &
            obug, ustarg, z0mg, alpha, bee, 1.)

         !print *, "------------------------"
         !print *, "rd :", rd
         !print *, "rd_:", rd_

         ! calculate ra2m, rd2m
         ra2m = frd(ktop, hroof, 0., displau+z0mu, 2., displa/hroof, z0h_g, &
            obug, ustarg, z0mg, alpha, bee, 1.)

         rd2m = frd(ktop, hroof, 0., 2., z0qg, displa/hroof, z0h_g, &
            obug, ustarg, z0mg, alpha, bee, 1.)

         ! Masson, 2000: Account for different canyon orientations
         ! 2/PI is a factor derived from 0-360deg integration
         IF (alpha_opt == 1) THEN
            ueff_lay(2) = 2/PI*ueff_lay(2)
            rd(:)       = PI/2*rd(:)
         ENDIF

!-----------------------------------------------------------------------
! Bulk boundary layer resistance of leaves
!-----------------------------------------------------------------------

         rb(:) = 0.

         DO i = 0, nurb
            clev  = canlev(i)
            rb(i) = rhoair * cpair / ( 11.8 + 4.2*ueff_lay(clev) )
         ENDDO

!-----------------------------------------------------------------------
! dimensional and non-dimensional sensible and latent heat conductances
! for canopy and soil flux calculations.
!-----------------------------------------------------------------------

         !NOTE: 0: roof, 1: sunlit wall, 2: shaded wall,
         !      3: impervious road, 4: pervious road, 5: vegetation
         cfh(:) = 0.
         cfw(:) = 0.

         DO i = 0, nurb
            cfh(i) = 1 / rb(i)

            IF (i == 0) THEN !roof
               ! account for fwet
               cfw(i) = fwet_roof / rb(i)
            ELSE
               cfw(i) = 1 / rb(i)
            ENDIF
         ENDDO

         ! For simplicity, there is no water exchange on the wall
         cfw(1:2) = 0.

         ! initialization
         cah(:) = 0.
         caw(:) = 0.
         cgh(:) = 0.
         cgw(:) = 0.

         ! conductance for each layer
         DO i = 3, 2, -1
            IF (i == 3) THEN
               cah(i) = 1. / rah
               caw(i) = 1. / raw
            ELSE
               cah(i) = 1. / rd(i+1)
               caw(i) = 1. / rd(i+1)
            ENDIF

            cgh(i) = 1. / rd(i)
            cgw(i) = 1. / rd(i)
         ENDDO

         ! claculate wtshi, wtsqi
         wtshi(:) = cah(:)*fah(:) + cgh(:)*fgh(:)
         wtsqi(:) = caw(:)*faw(:) + cgw(:)*fgw(:)

         DO i = 0, nurb
            clev = canlev(i)
            wtshi(clev) = wtshi(clev) + fc(i)*cfh(i)
            wtsqi(clev) = wtsqi(clev) + fc(i)*cfw(i)
         ENDDO

         DO i = 3, 2, -1
            wtshi(i) = 1./wtshi(i)
            wtsqi(i) = 1./wtsqi(i)
         ENDDO

         wta0(:)  = cah(:) * wtshi(:) * fah(:)
         wtg0(:)  = cgh(:) * wtshi(:) * fgh(:)

         wtaq0(:) = caw(:) * wtsqi(:) * faw(:)
         wtgq0(:) = cgw(:) * wtsqi(:) * fgw(:)

         ! calculate wtl0, wtll, wtlq0, wtlql
         wtll(:)  = 0.
         wtlql(:) = 0.

         DO i = 0, nurb
            clev = canlev(i)

            wtl0(i)    = cfh(i) * wtshi(clev) * fc(i)
            wtll(clev) = wtll(clev) + wtl0(i)*tu(i)

            wtlq0(i)    = cfw(i) * wtsqi(clev) * fc(i)
            wtlql(clev) = wtlql(clev) + wtlq0(i)*qsatl(i)
         ENDDO

         IF (numlay .eq. 2) THEN

            ! - Equations:
            ! taf(3) = (1/rah*thm + 1/rd(3)*taf(2) + 1/rb(0)*troof*fc(0) + AHE/(rho*cp))/(1/rah + 1/rd(3) + 1/rb(0)*fc(0))
            ! taf(2) = (1/rd(3)*taf(3) + 1/rd(2)*tg*fg + 1/rb(1)*twsun*fc(1) + 1/rb(2)*twsha*fc(2) + AHE/(rho*cp))/ &
            !          (1/rd(3) + 1/rd(2)*fg + 1/rb(1)*fc(1) + 1/rb(2)*fc(2))
            ! Also written as:
            ! taf(3) = (cah(3)*thm + cah(2)*taf(2) + cfh(0)*troof*fc(0))/(cah(3) + cah(2) + cfh(0)*fc(0))
            ! taf(2) = (cah(2)*taf(3) + cgh(2)*tg*fg + cfh(1)*twsun*fc(1) + cfh(2)*twsha*fc(2) + AHE/(rho*cp))/ &
            !          (cah(2) + cgh(2)*fg + cfh(1)*fc(1) + cfh(2)*fc(2))

            ! - Equations:
            ! qaf(3) = (1/raw*qm + 1/rd(3)*qaf(2) + 1/rb(0)*qroof*fc(0))/(1/raw + 1/rd(3) + 1/rb(0)*fc(0))
            ! qaf(2) = (1/rd(3)*qaf(3) + 1/(rd(2)+rsr)*qper*fgper*fg + fwetimp/rd(2)*qimp*fgimp*fg + AHE/rho)/ &
            !          (1/rd(3) + 1/(rd(2)+rsr)*fgper*fg + fwetimp/rd(2)*fgimp*fg)
            ! Also written as:
            ! qaf(3) = (caw(3)*qm + caw(2)*qaf(2) + cfw(0)*qroof*fc(0))/(caw(3) + caw(2) + cfw(0)*fc(0))
            ! qaf(2) = (caw(2)*qaf(3) + cgwper*qper*fgper*fg + cgwimp*qimp*fgimp*fg + AHE/rho)/ &
            !          (caw(2) + cgwper*fgper*fg + cgwimp*fgimp*fg)

            ! 06/20/2021, yuan: account for Anthropogenic heat
            ! 92% heat release as SH, Pigeon et al., 2007

            h_vec  = vehc
            tmpw1  = cah(2)*((cah(3)*thm + cfh(0)*tu(0)*fc(0) + 1/(4*hlr+1)*(Fhac+Fwst)/(rhoair*cpair))/&
                     (cah(3) + cah(2) + cfh(0)*fc(0)))
            tmpw2  = (4*hlr/(4*hlr+1)*(Fhac+Fwst)+Fach)/(rhoair*cpair) + (h_vec+meta)/(rhoair*cpair)
            tmpw3  = cgh(2)*fg*tg + cfh(1)*tu(1)*fc(1) + cfh(2)*tu(2)*fc(2)
            fact   = 1. - (cah(2)*cah(2)/(cah(3) + cah(2) + cfh(0)*fc(0))/&
                     (cah(2) + cgh(2)*fg + cfh(1)*fc(1) + cfh(2)*fc(2)))
            taf(2) = (tmpw1 + tmpw2 + tmpw3) / &
                     (cah(2) + cgh(2)*fg + cfh(1)*fc(1) + cfh(2)*fc(2)) / &
                     fact

            IF (qgper < qaf(2)) THEN
              ! dew case. no soil resistance
              cgw_per= cgw(2)
            ELSE
              cgw_per= 1/(1/cgw(2)+rsr)
            ENDIF

            cgw_imp= fwet_gimp*cgw(2)

            ! account for soil resistance, qgper and qgimp are calculated separately
            l_vec  = 0
            tmpw1  = caw(2)*((caw(3)*qm + cfw(0)*qsatl(0)*fc(0))/&
                     (caw(3) + caw(2) + cfw(0)*fc(0)))
            tmpw2  = l_vec/(rhoair)
            tmpw3  = cgw_per*qgper*fgper*fg + cgw_imp*qgimp*fgimp*fg
            facq   = 1. - (caw(2)*caw(2)/&
                     (caw(3) + caw(2) + cfw(0)*fc(0))/&
                     (caw(2) + cgw_per*fgper*fg + cgw_imp*fgimp*fg))
            qaf(2) = (tmpw1 + tmpw2 + tmpw3)/&
                     (caw(2) + cgw_per*fgper*fg + cgw_imp*fgimp*fg)/&
                     facq

            tmpw1  = 1/(4*hlr+1)*(Fhac+Fwst)/(rhoair*cpair)
            taf(3) = (cah(3)*thm + cah(2)*taf(2) + cfh(0)*tu(0)*fc(0) + tmpw1)/&
                     (cah(3) + cah(2) + cfh(0)*fc(0))
            qaf(3) = (caw(3)*qm  + caw(2)*qaf(2) + cfw(0)*qsatl(0)*fc(0))/&
                     (caw(3) + caw(2) + cfw(0)*fc(0))
         ENDIF

         !------------------------------------------------
         ! update fwet for roof and impervious ground
         ! to check whether dew happens
         IF (qaf(3) > qroof) THEN
            fwet_roof = 1. !dew case
         ELSE
            fwet_roof = fwet_roof_
         ENDIF

         ! to check whether dew happens
         IF (qaf(2) > qgimp) THEN
            fwet_gimp = 1. !dew case
         ELSE
            fwet_gimp = fwet_gimp_
         ENDIF

         ! weighted qg
         ! NOTE: IF fwet_gimp=1, same as previous
         fwetfac = fgimp*fwet_gimp + fgper
         qg = (qgimp*fgimp*fwet_gimp + qgper*fgper) / fwetfac

         fgw(2) = fg*fwetfac

!-----------------------------------------------------------------------
! Update monin-obukhov length and wind speed including the stability effect
!-----------------------------------------------------------------------

         ! USE the top layer taf and qaf
         !TODO: need more check
         dth = thm - taf(2)
         dqh =  qm - qaf(2)

         tstar = vonkar/(fh)*dth
         qstar = vonkar/(fq)*dqh

         thvstar = tstar*(1.+0.61*qm)+0.61*th*qstar
         zeta = zldis*vonkar*grav*thvstar / (ustar**2*thv)
         IF (zeta .ge. 0.) THEN                             !stable
            zeta = min(2.,max(zeta,1.e-6))
         ELSE                                             !unstable
            zeta = max(-100.,min(zeta,-1.e-6))
         ENDIF
         obu = zldis/zeta

         IF (zeta .ge. 0.) THEN
            um  = max(ur,.1)
         ELSE
            wc  = (-grav*ustar*thvstar*zii/thv)**(1./3.)
            wc2 = beta*beta*(wc*wc)
            um  = sqrt(ur*ur+wc2)
         ENDIF

         IF (obuold*obu .lt. 0.) nmozsgn = nmozsgn+1
         IF (nmozsgn >= 4) EXIT

         obuold = obu

      ENDDO ITERATION                         !end stability iteration

! ======================================================================
!     END stability iteration
! ======================================================================

      zol = zeta
      rib = min(5.,zol*ustar**2/(vonkar**2/fh*um**2))

      ! sensible heat fluxes
      fsenroof = rhoair*cpair*cfh(0)*(troof-taf(3))
      fsenwsun = rhoair*cpair*cfh(1)*(twsun-taf(2))
      fsenwsha = rhoair*cpair*cfh(2)*(twsha-taf(2))

      ! latent heat fluxes
      fevproof = rhoair*cfw(0)*(qsatl(0)-qaf(3))
      fevproof = fevproof*fwet_roof

      ! fact   = 1. - wta0(2)*wtg0(3)
      ! facq   = 1. - wtaq0(2)*wtgq0(3)
      ! deduce: croofs = rhoair*cpair*cfh(0)*(1.-wtg0(3)*wta0(2)*wtl0(0)/fact-wtl0(0))
      croofs = rhoair*cpair*cfh(0)*(1.-wtl0(0)/fact)
      cwalls = rhoair*cpair*cfh(1)*(1.-wtl0(1)/fact)
      ! deduce: croofl = rhoair*cfw(0)*(1.-wtgq0(3)*wtaq0(2)*wtlq0(0)/facq-wtlq0(0))*qsatldT(0)
      ! croofl = rhoair*cfw(0)*(1.-wtlq0(0)/facq)*qsatldT(0)
      croofl = rhoair*cfw(0)*(1.-cfw(0)*fc(0)/(caw(3)+cgw(3)+cfw(0)*fc(0))-cgw(3) &
               /(caw(3)+cgw(3)+cfw(0)*fc(0)) &
               /(cgw(3)+cgw_per*fgper*fg+cgw_imp*fgimp*fg)* &
               cfw(0)*fc(0)*cgw(3)/(caw(3)+cgw(3)+cfw(0)*fc(0))/facq)*qsatldT(0)
      croofl = croofl*fwet_roof

      croof = croofs + croofl*htvp_roof

#if(defined CoLMDEBUG)
#endif

      tafu = taf(2)

!-----------------------------------------------------------------------
! wind stresses
!-----------------------------------------------------------------------

      taux = - rhoair*us/ram
      tauy = - rhoair*vs/ram

!-----------------------------------------------------------------------
! fluxes from urban ground to canopy space
!-----------------------------------------------------------------------

      fsengper = cpair*rhoair*cgh(2)*(tgper-taf(2))
      fsengimp = cpair*rhoair*cgh(2)*(tgimp-taf(2))

      fevpgper = rhoair*cgw_per*(qgper-qaf(2))
      fevpgimp = rhoair*cgw_imp*(qgimp-qaf(2))
      fevpgimp = fevpgimp*fwet_gimp

!-----------------------------------------------------------------------
! Derivative of soil energy flux with respect to soil temperature (cgrnd)
!-----------------------------------------------------------------------

      cgrnds = cpair*rhoair*cgh(2)*(1.-wtg0(2)/fact)
      ! cgperl = rhoair*cgw(2)*(1.-wtgq0(2)/facq)*dqgperdT
      ! cgimpl = rhoair*cgw(2)*(1.-wtgq0(2)/facq)*dqgimpdT

      cgperl = rhoair*cgw_per*(dqgperdT &
               - (dqgperdT*cgw_per*fgper*fg) &
               /(caw(2) + cgw_per*fgper*fg + cgw_imp*fgimp*fg) &
               /facq)
      cgimpl = rhoair*cgw_imp*(dqgimpdT &
               - (dqgimpdT*cgw_imp*fgimp*fg) &
               /(caw(2) + cgw_per*fgper*fg + cgw_imp*fgimp*fg) &
               /facq)
      cgimpl = cgimpl*fwet_gimp

      cgimp  = cgrnds + cgimpl*htvp_gimp
      cgper  = cgrnds + cgperl*htvp_gper

!-----------------------------------------------------------------------
! 2 m height air temperature above apparent sink height
!-----------------------------------------------------------------------

      !tref = thm + vonkar/(fh-fht)*dth * (fh2m/vonkar - fh/vonkar)
      !qref =  qm + vonkar/(fq-fqt)*dqh * (fq2m/vonkar - fq/vonkar)

   END SUBROUTINE UrbanOnlyFlux


   SUBROUTINE  UrbanVegFlux ( &
         ! Model running information
         ipatch      ,deltim      ,lbr         ,lbi         ,&
         ! Forcing
         hu          ,ht          ,hq          ,us          ,&
         vs          ,thm         ,th          ,thv         ,&
         qm          ,psrf        ,rhoair      ,frl         ,&
         po2m        ,pco2m       ,par         ,sabv        ,&
         rstfac      ,Fhac        ,Fwst        ,Fach        ,&
         vehc        ,meta                                  ,&
         ! Urban and vegetation parameters
         hroof       ,hwr         ,nurb        ,fcover      ,&
         ewall       ,egimp       ,egper       ,ev          ,&
         htop        ,hbot        ,lai         ,sai         ,&
         sqrtdi      ,effcon      ,vmax25      ,slti        ,&
         hlti        ,shti        ,hhti        ,trda        ,&
         trdm        ,trop        ,g1          ,g0          ,&
         gradm       ,binter      ,extkn       ,extkd       ,&
         dewmx       ,etrc        ,&
         ! Status of surface
         z0h_g       ,obug        ,ustarg      ,zlnd        ,&
         zsno        ,fsno_roof   ,fsno_gimp   ,fsno_gper   ,&
         wliq_roofsno,wliq_gimpsno,wice_roofsno,wice_gimpsno,&
         htvp_roof   ,htvp_gimp   ,htvp_gper   ,troof       ,&
         twsun       ,twsha       ,tgimp       ,tgper       ,&
         qroof       ,qgimp       ,qgper       ,dqroofdT    ,&
         dqgimpdT    ,dqgperdT    ,sigf        ,tl          ,&
         ldew        ,rsr                                   ,&
         ! Longwave information
         Ainv        ,B           ,B1          ,dBdT        ,&
         SkyVF       ,VegVF                                 ,&
         ! Output
         taux        ,tauy        ,fsenroof    ,fsenwsun    ,&
         fsenwsha    ,fsengimp    ,fsengper    ,fevproof    ,&
         fevpgimp    ,fevpgper    ,croofs      ,cwalls      ,&
         cgrnds      ,croofl      ,cgimpl      ,cgperl      ,&
         croof       ,cgimp       ,cgper       ,fsenl       ,&
         fevpl       ,etr         ,rst         ,assim       ,&
         respc       ,lwsun       ,lwsha       ,lgimp       ,&
         lgper       ,lveg        ,lout        ,tref        ,&
         qref        ,z0m         ,zol         ,rib         ,&
         ustar       ,qstar       ,tstar       ,fm          ,&
         fh          ,fq          ,tafu                      )

!=======================================================================

   USE MOD_Precision
   USE MOD_Const_Physical, only: vonkar,grav,hvap,cpair,stefnc
   USE MOD_FrictionVelocity
   USE MOD_CanopyLayerProfile
   USE MOD_AssimStomataConductance
   IMPLICIT NONE

!-----------------------Arguments---------------------------------------
   integer,  intent(in) :: &
        ipatch,   &! patch index [-]
        lbr,      &! lower bound of array
        lbi        ! lower bound of array

   real(r8), intent(in) :: &
        deltim     ! seconds in a time step [second]

   ! Forcing
   real(r8), intent(in) :: &
        hu,       &! observational height of wind [m]
        ht,       &! observational height of temperature [m]
        hq,       &! observational height of humidity [m]
        us,       &! wind component in eastward direction [m/s]
        vs,       &! wind component in northward direction [m/s]
        thm,      &! intermediate variable (tm+0.0098*ht)
        th,       &! potential temperature (kelvin)
        thv,      &! virtual potential temperature (kelvin)
        qm,       &! specific humidity at reference height [kg/kg]
        psrf,     &! pressure at reference height [pa]
        rhoair,   &! density air [kg/m**3]

        frl,      &! atmospheric infrared (longwave) radiation [W/m2]
        par,      &! par absorbed per unit sunlit lai [w/m**2]
        sabv,     &! solar radiation absorbed by vegetation [W/m2]
        rstfac,   &! factor of soil water stress to plant physiologocal processes

        po2m,     &! atmospheric partial pressure  o2 (pa)
        pco2m,    &! atmospheric partial pressure co2 (pa)

        vehc,     &! flux from vehicle
        meta,     &! flux from metabolic
        Fhac,     &! flux from heat or cool AC
        Fwst,     &! waste heat from cool or heat
        Fach       ! flux from air exchange

   ! Urban and vegetation parameters
   integer,  intent(in) :: &
        nurb       ! number of aboveground urban components [-]

   real(r8), intent(in) :: &
        hroof,    &! average building height [m]
        hwr,      &! average building height to their distance [-]
        fcover(0:5)! coverage of aboveground urban components [-]

   real(r8), intent(in) :: &
        ewall,    &! emissivity of walls
        egimp,    &! emissivity of impervious road
        egper,    &! emissivity of pervious road
        ev         ! emissivity of vegetation

   real(r8), intent(in) :: &
        htop,     &! PFT crown top height [m]
        hbot,     &! PFT crown bottom height [m]
        lai,      &! adjusted leaf area index for seasonal variation [-]
        sai,      &! stem area index  [-]
        sqrtdi,   &! inverse sqrt of leaf dimension [m**-0.5]

        effcon,   &! quantum efficiency of RuBP regeneration (mol CO2 / mol quanta)
        vmax25,   &! maximum carboxylation rate at 25 C at canopy top
                   ! the range : 30.e-6 <-> 100.e-6 (mol co2 m-2 s-1)
        shti,     &! slope of high temperature inhibition function     (s1)
        hhti,     &! 1/2 point of high temperature inhibition function (s2)
        slti,     &! slope of low temperature inhibition function      (s3)
        hlti,     &! 1/2 point of low temperature inhibition function  (s4)
        trda,     &! temperature coefficient in gs-a model             (s5)
        trdm,     &! temperature coefficient in gs-a model             (s6)
        trop,     &! temperature coefficient in gs-a model         (273+25)
        g1,       &! conductance-photosynthesis slope parameter for medlyn model
        g0,       &! conductance-photosynthesis intercept for medlyn model
        gradm,    &! conductance-photosynthesis slope parameter
        binter,   &! conductance-photosynthesis intercept

        extkn,    &! coefficient of leaf nitrogen allocation
        extkd,    &! diffuse and scattered diffuse PAR extinction coefficient
        dewmx,    &! maximum dew
        etrc       ! maximum possible transpiration rate (mm/s)

   ! Status of surface
   real(r8), intent(in) :: &
        rsr,      &! bare soil resistance for evaporation
        z0h_g,    &! roughness length for bare ground, sensible heat [m]
        obug,     &! monin-obukhov length for bare ground (m)
        ustarg,   &! friction velocity for bare ground [m/s]
        zlnd,     &! roughness length for soil [m]
        zsno,     &! roughness length for snow [m]
        fsno_roof,&! fraction of ground covered by snow
        fsno_gimp,&! fraction of ground covered by snow
        fsno_gper,&! fraction of ground covered by snow
        wliq_roofsno,&! liqui water [kg/m2]
        wliq_gimpsno,&! liqui water [kg/m2]
        wice_roofsno,&! ice lens [kg/m2]
        wice_gimpsno,&! ice lens [kg/m2]
        htvp_roof,&! latent heat of vapor of water (or sublimation) [j/kg]
        htvp_gimp,&! latent heat of vapor of water (or sublimation) [j/kg]
        htvp_gper,&! latent heat of vapor of water (or sublimation) [j/kg]

        troof,    &! temperature of roof [K]
        twsun,    &! temperature of sunlit wall [K]
        twsha,    &! temperature of shaded wall [K]
        tgimp,    &! temperature of impervious road [K]
        tgper,    &! pervious ground temperature [K]

        qroof,    &! roof specific humidity [kg/kg]
        qgimp,    &! imperivous road specific humidity [kg/kg]
        qgper,    &! pervious ground specific humidity [kg/kg]
        dqroofdT, &! d(qroof)/dT
        dqgimpdT, &! d(qgimp)/dT
        dqgperdT, &! d(qgper)/dT
        sigf       !

   real(r8), intent(inout) :: &
        tl,       &! leaf temperature [K]
        ldew       ! depth of water on foliage [mm]

   real(r8), intent(in)    :: Ainv(5,5)  !Inverse of Radiation transfer matrix
   real(r8), intent(in)    :: SkyVF (5)  !View factor to sky
   real(r8), intent(in)    :: VegVF (5)  !View factor to veg
   real(r8), intent(inout) :: B     (5)  !Vectors of incident radition on each surface
   real(r8), intent(inout) :: B1    (5)  !Vectors of incident radition on each surface
   real(r8), intent(inout) :: dBdT  (5)  !Vectors of incident radition on each surface

   real(r8), intent(out) :: &
        taux,     &! wind stress: E-W [kg/m/s**2]
        tauy,     &! wind stress: N-S [kg/m/s**2]
        fsenroof, &! sensible heat flux from roof [W/m2]
        fsenwsun, &! sensible heat flux from sunlit wall [W/m2]
        fsenwsha, &! sensible heat flux from shaded wall [W/m2]
        fsengimp, &! sensible heat flux from impervious road [W/m2]
        fsengper, &! sensible heat flux from pervious ground [W/m2]
        fevproof, &! evaporation heat flux from roof [mm/s]
        fevpgimp, &! evaporation heat flux from impervious road [mm/s]
        fevpgper, &! evaporation heat flux from pervious ground [mm/s]

        croofs,   &! deriv of roof sensible heat flux wrt soil temp [w/m**2/k]
        cwalls,   &! deriv of wall sensible heat flux wrt soil temp [w/m**2/k]
        cgrnds,   &! deriv of ground latent heat flux wrt soil temp [w/m**2/k]
        croofl,   &! deriv of roof latent heat flux wrt soil temp [w/m**2/k]
        cgimpl,   &! deriv of impervious latent heat flux wrt soil temp [w/m**2/k]
        cgperl,   &! deriv of soil atent heat flux wrt soil temp [w/m**2/k]
        croof,    &! deriv of roof total flux wrt soil temp [w/m**2/k]
        cgimp,    &! deriv of impervious total heat flux wrt soil temp [w/m**2/k]
        cgper,    &! deriv of soil total heat flux wrt soil temp [w/m**2/k]

        tref,     &! 2 m height air temperature [kelvin]
        qref       ! 2 m height air humidity

   real(r8), intent(out) :: &
        fsenl,    &! sensible heat from leaves [W/m2]
        fevpl,    &! evaporation+transpiration from leaves [mm/s]
        etr,      &! transpiration rate [mm/s]
        rst,      &! stomatal resistance
        assim,    &! rate of assimilation
        respc      ! rate of respiration

   real(r8), intent(inout) :: &
        lwsun,    &! net longwave radiation of sunlit wall
        lwsha,    &! net longwave radiation of shaded wall
        lgimp,    &! net longwave radiation of impervious road
        lgper,    &! net longwave radiation of pervious road
        lveg,     &! net longwave radiation of vegetation
        lout       ! out-going longwave radiation

   real(r8), intent(inout) :: &
        z0m,      &! effective roughness [m]
        zol,      &! dimensionless height (z/L) used in Monin-Obukhov theory
        rib,      &! bulk Richardson number in surface layer
        ustar,    &! friction velocity [m/s]
        tstar,    &! temperature scaling parameter
        qstar,    &! moisture scaling parameter
        fm,       &! integral of profile function for momentum
        fh,       &! integral of profile function for heat
        fq,       &! integral of profile function for moisture
        tafu       ! effective urban air temperature (2nd layer, walls)

!-----------------------Local Variables---------------------------------
! assign iteration parameters
   integer, parameter :: itmax  = 40   !maximum number of iteration
   integer, parameter :: itmin  = 6    !minimum number of iteration
   real(r8),parameter :: delmax = 3.0  !maximum change in leaf temperature [K]
   real(r8),parameter :: dtmin  = 0.01 !max limit for temperature convergence [K]
   real(r8),parameter :: dlemin = 0.1  !max limit for energy flux convergence [w/m2]

   real(r8) dtl(0:itmax+1)             !difference of tl between two iterative step

   real(r8) ::  &
        zldis,    &! reference height "minus" zero displacement heght [m]
        zii,      &! convective boundary layer height [m]
        z0mv,     &! roughness length of vegetation only, momentum [m]
        z0mu,     &! roughness length of building only, momentum [m]
        z0h,      &! roughness length, sensible heat [m]
        z0q,      &! roughness length, latent heat [m]
        zeta,     &! dimensionless height used in Monin-Obukhov theory
        beta,     &! coefficient of conective velocity [-]
        wc,       &! convective velocity [m/s]
        wc2,      &! wc**2
        dth,      &! diff of virtual temp. between ref. height and surface
        dthv,     &! diff of vir. poten. temp. between ref. height and surface
        dqh,      &! diff of humidity between ref. height and surface
        obu,      &! monin-obukhov length (m)
        um,       &! wind speed including the stablity effect [m/s]
        ur,       &! wind speed at reference height [m/s]
        uaf,      &! velocity of air within foliage [m/s]
        fh2m,     &! relation for temperature at 2m
        fq2m,     &! relation for specific humidity at 2m
        fm10m,    &! integral of profile function for momentum at 10m
        thvstar,  &! virtual potential temperature scaling parameter
        eah,      &! canopy air vapor pressure (pa)
        pco2g,    &! co2 pressure (pa) at ground surface (pa)
        pco2a,    &! canopy air co2 pressure (pa)

        ram,      &! aerodynamical resistance [s/m]
        rah,      &! thermal resistance [s/m]
        raw,      &! moisture resistance [s/m]
        clai,     &! canopy heat capacity [Jm-2K-1]
        del,      &! absolute change in leaf temp in current iteration [K]
        del2,     &! change in leaf temperature in previous iteration [K]
        dele,     &! change in heat fluxes from leaf [K]
        dele2,    &! change in heat fluxes from leaf [K]
        det,      &! maximum leaf temp. change in two consecutive iter [K]
        dee,      &! maximum leaf temp. change in two consecutive iter [K]

        obuold,   &! monin-obukhov length from previous iteration
        tlbef,    &! leaf temperature from previous iteration [K]
        err,      &! balance error

        rs,       &! sunlit leaf stomatal resistance [s/m]
        rsoil,    &! soil respiration
        gah2o,    &! conductance between canopy and atmosphere
        gdh2o,    &! conductance between canopy and ground
        tprcor     ! tf*psur*100./1.013e5

   integer it, nmozsgn

   real(r8) evplwet, evplwet_dtl, etr_dtl, elwmax, elwdif
   real(r8) irab, dirab_dtl, fsenl_dtl, fevpl_dtl
   real(r8) z0mg, z0hg, z0qg, cint(3)
   real(r8) fevpl_bef, fevpl_noadj, dtl_noadj, erre

!----------------------- defination for 3d run ------------------------ !
   integer, parameter :: nlay = 3
   integer, parameter :: uvec(5) = (/0,0,0,0,1/) !unit vector

   integer :: &
        clev,     &! current layer index
        botlay,   &! botom layer index
        numlay     ! available layer number

   real(r8) :: &
        huu,      &! observational height of wind [m]
        htu,      &! observational height of temperature [m]
        hqu,      &! observational height of humidity [m]
        ktop,     &! K value at a specific height
        utop,     &! u value at a specific height
        fht,      &! integral of profile function for heat at the top layer
        fqt,      &! integral of profile function for moisture at the top layer
        fmtop,    &! fm value at a specific height
        phih,     &! phi(h), similarity function for sensible heat
        displa,   &! displacement height for urban
        displau,  &! displacement height for urban building
        displav,  &! displacement height for urban vegetation
        displav_lay,&!displacement height for urban vegetation layer
        z0mv_lay, &! roughless length for vegetation
        ueff_veg, &! effective wind speed within canopy layer [m/s]
        tg,       &! ground temperature
        qg         ! ground specific humidity

   real(r8) :: &
        fg,       &! ground fractional cover
        fgimp,    &! weight of impervious ground
        fgper,    &! weight of pervious ground
        hlr,      &! average building height to their length of edge [-]
        sqrtdragc,&! sqrt(drag coefficient)
        lm,       &! mix length within canopy
        fai,      &! frontal area index for urban
        faiv,     &! frontal area index for trees
        lsai,     &! lai+sai
        fwet,     &! fractional wet area
        delta,    &! 0 or 1
        alpha,    &! exponential extinction factor for u/k decline within urban
        alphav     ! exponential extinction factor for u/k decline within trees

   real(r8) :: &
        lwsun_bef,&! change of lw for the last time
        lwsha_bef,&! change of lw for the last time
        lgimp_bef,&! change of lw for the last time
        lgper_bef,&! change of lw for the last time
        lveg_bef   ! change of lw for the last time

   real(r8), dimension(0:nurb) :: &
        tu,       &! termperature array
        fc,       &! fractional cover array
        canlev,   &! urban canopy layer lookup table
        rb,       &! leaf boundary layer resistance [s/m]
        cfh,      &! heat conductance for leaf [m/s]
        cfw,      &! latent heat conductance for leaf [m/s]
        wtl0,     &! normalized heat conductance for air and leaf [-]
        wtlq0,    &! normalized latent heat cond. for air and leaf [-]

        ei,       &! vapor pressure on leaf surface [pa]
        deidT,    &! derivative of "ei" on "tl" [pa/K]
        qsatl,    &! leaf specific humidity [kg/kg]
        qsatldT    ! derivative of "qsatl" on "tlef"

   real(r8), dimension(nlay) :: &
        fah,      &! weight for thermal resistance to upper layer
        faw,      &! weight for moisture resistance to upper layer
        fgh,      &! weight for thermal resistance to lower layer
        fgw,      &! weight for moisture resistance to lower layer
        ueff_lay, &! effective wind speed within canopy layer [m/s]
        ueff_lay_,&! effective wind speed within canopy layer [m/s]
        taf,      &! air temperature within canopy space [K]
        qaf,      &! humidity of canopy air [kg/kg]
        rd,       &! aerodynamic resistance between layers [s/m]
        rd_,      &! aerodynamic resistance between layers [s/m]
        cah,      &! heat conductance for air [m/s]
        cgh,      &! heat conductance for ground [m/s]
        caw,      &! latent heat conductance for air [m/s]
        cgw,      &! latent heat conductance for ground [m/s]
        wtshi,    &! sensible heat resistance for air, grd and leaf [-]
        wtsqi,    &! latent heat resistance for air, grd and leaf [-]
        wta0,     &! normalized heat conductance for air [-]
        wtg0,     &! normalized heat conductance for ground [-]
        wtaq0,    &! normalized latent heat conductance for air [-]
        wtgq0,    &! normalized heat conductance for ground [-]
        wtll,     &! sum of normalized heat conductance for air and leaf
        wtlql      ! sum of normalized heat conductance for air and leaf

   real(r8) :: &
        ra2m,     &! aerodynamic resistance between 2m and bottom layer [s/m]
        rd2m       ! aerodynamic resistance between bottom layer and ground [s/m]

   ! temporal
   integer i
   real(r8) bee, cf, tmpw1, tmpw2, tmpw3, tmpw4, fact, facq, taftmp
   real(r8) B_5, B1_5, dBdT_5, X(5), dX(5)
   real(r8) fwet_roof, fwet_roof_, fwet_gimp, fwet_gimp_
   real(r8) fwetfac, lambda
   real(r8) cgw_imp, cgw_per
   real(r8) h_vec, l_vec

   ! for interface
   real(r8) o3coefv,o3coefg,assim_RuBP, assim_Rubisco, ci, vpd, gammas

!-----------------------End Variable List-------------------------------

! initialization of errors and  iteration parameters
      it    = 1    !counter for leaf temperature iteration
      del   = 0.0  !change in leaf temperature from previous iteration
      dele  = 0.0  !latent head flux from leaf for previous iteration

      dtl   = 0.
      fevpl_bef = 0.

! initial values for z0hg, z0qg

      !TODO: change to original
      !z0mg = (1.-fsno)*zlnd + fsno*zsno
      IF (fsno_gper > 0) THEN
         z0mg = zsno
      ELSE
         z0mg = zlnd
      ENDIF
      z0hg = z0mg
      z0qg = z0mg

!-----------------------------------------------------------------------
! scaling-up coefficients from leaf to canopy
!-----------------------------------------------------------------------

      cint(1) = (1.-exp(-0.110*lai))/0.110
      cint(2) = (1.-exp(-extkd*lai))/extkd
      cint(3) = lai

!-----------------------------------------------------------------------
! initial saturated vapor pressure and humidity and their derivation
!-----------------------------------------------------------------------

      !clai = 4.2 * 1000. * 0.2
      clai = 0.0
      lsai = lai + sai

      ! index 0:roof, 1:sunlit wall, 2:shaded wall, 3: vegetation
      tu(0) = troof; tu(1) = twsun; tu(2) = twsha; tu(3) = tl

      fg     = 1 - fcover(0)
      fc(:)  = fcover(0:nurb)
      fc(3)  = fcover(5)
      fgimp  = fcover(3)/fg
      fgper  = fcover(4)/fg
      hlr    = hwr*(1-sqrt(fcover(0)))/sqrt(fcover(0))
      canlev = (/3, 2, 2, 1/)

      B_5    = B(5)
      B1_5   = B1(5)
      dBdT_5 = dBdT(5)

      CALL dewfraction (sigf,lai,sai,dewmx,ldew,fwet)

      qsatl(0) = qroof
      qsatldT(0) = dqroofDT
      DO i = 1, nurb
         CALL qsadv(tu(i),psrf,ei(i),deiDT(i),qsatl(i),qsatldT(i))
      ENDDO

      ! Save the longwave for the last time
      lwsun_bef = lwsun
      lwsha_bef = lwsha
      lgimp_bef = lgimp
      lgper_bef = lgper
      lveg_bef  = lveg

!-----------------------------------------------------------------------
! Calculate the weighted qg, tg
!-----------------------------------------------------------------------

      ! set weghting factor
      fah(1) = 1.; fah(2) = 1.; fah(3) = 1.
      faw(1) = 1.; faw(2) = 1.; faw(3) = 1.
      fgh(1) = 1.; fgh(2) = 1.; fgh(3) = 1.
      fgw(1) = 1.; fgw(2) = 1.; fgw(3) = 1.

      ! weighted tg and qg
      tg = tgimp*fgimp + tgper*fgper

      ! wet fraction for roof and impervious ground
      !-------------------------------------------
      ! roof
      IF (lbr < 1) THEN
         fwet_roof_ = fsno_roof !for snow layer exist
      ELSE
         ! surface wet fraction. assuming max ponding = 1 kg/m2
         fwet_roof_ = (max(0., wliq_roofsno+wice_roofsno))**(2/3.)
         fwet_roof_ = min(1., fwet_roof_)
      ENDIF

      ! impervious ground
      IF (lbi < 1) THEN
         fwet_gimp_ = fsno_gimp !for snow layer exist
      ELSE
         ! surface wet fraction. assuming max ponding = 1 kg/m2
         fwet_gimp_ = (max(0., wliq_gimpsno+wice_gimpsno))**(2/3.)
         fwet_gimp_ = min(1., fwet_gimp_)
      ENDIF

      ! dew case
      IF (qm > qroof) THEN
         fwet_roof = 1.
      ELSE
         fwet_roof = fwet_roof_
      ENDIF

      ! dew case
      IF (qm > qgimp) THEN
         fwet_gimp = 1.
      ELSE
         fwet_gimp = fwet_gimp_
      ENDIF

      ! weighted qg
      ! NOTE: IF fwet_gimp=1, same as previous
      fwetfac = fgimp*fwet_gimp + fgper
      qg = (qgimp*fgimp*fwet_gimp + qgper*fgper) / fwetfac

      fgw(2) = fg*fwetfac

!-----------------------------------------------------------------------
! initial for fluxes profile
!-----------------------------------------------------------------------

      nmozsgn = 0     !number of times moz changes sign
      obuold  = 0.    !monin-obukhov length from previous iteration
      zii     = 1000. !m (pbl height)
      beta    = 1.    !- (in computing W_*)

!-----------------------------------------------------------------------
! scaling factor bee
!-----------------------------------------------------------------------
!NOTE: bee value, the default is 1
      bee = 1.

!-----------------------------------------------------------------------
! calculate z0m and displa for layers
!-----------------------------------------------------------------------

      ! Calculate z0 and displa for vegetation only and the whole area
      CALL cal_z0_displa(lsai, htop, 1., z0mv, displav)
      CALL cal_z0_displa(lsai, htop, fc(3), z0mv_lay, displav_lay)

      ! For building only below
      ! Macdonald et al., 1998, Eq. (23), A=4.43
      lambda  = fcover(0)
      displau = hroof * (1 + 4.43**(-lambda)*(lambda - 1))
      fai     = 4/PI*hlr*fcover(0)
      z0mu    = (hroof - displau) * &
                exp( -(0.5*1.2/vonkar/vonkar*(1-displau/hroof)*fai)**(-0.5) )

      ! account for vegetation
      faiv    = fc(3)*(1. - exp(-0.5*lsai))
      lambda  = fcover(0) + faiv*htop/hroof
      displa  = hroof * (1 + 4.43**(-lambda)*(lambda - 1))
      displa  = min(0.95*hroof, displa)
      z0m     = (hroof - displa) * &
                exp( -(0.5*1.2/vonkar/vonkar*(1-displa/hroof)*(fai+faiv*htop/hroof))**(-0.5) )

      ! to compare z0 of urban and only the surface
      ! maximum assumption
      ! 11/26/2021, yuan: remove the below
      !IF (z0mu < z0mv_lay) z0mu = z0mv_lay
      !IF (displau < displav_lay) displau = displav_lay
      IF (z0m < z0mg) z0m = z0mg
      IF (displa >= hroof-z0mg) displa = hroof-z0mg

      ! minimum building displa limit
      displau = max(hroof/2., displau)

      ! Layer setting
      ! NOTE: right now only for 2 layers
      !IF (z0mv+displav > z0mu+displau) THEN
         numlay = 2; botlay = 2; canlev(3) = 2
         fgh(2) = fg; fgw(2) = fg;
      !ELSE
      !   numlay = 3; botlay = 1
      !   fgh(1) = fg; fgw(1) = fg;
      !ENDIF

!-----------------------------------------------------------------------
! calculate layer decay coefficient
!-----------------------------------------------------------------------

      ! Raupach, 1992
      sqrtdragc = min( (0.003+0.3*faiv)**0.5, 0.3 )

      ! Kondo, 1971
      alphav = htop/(htop-displav_lay)/(vonkar/sqrtdragc)
      alphav = alphav*htop/hroof

      ! Masson, 2000; Oleson et al., 2008 plus tree (+)
      IF (alpha_opt == 1) alpha = 0.5*hwr + alphav

      ! Swaid, 1993; Kusaka, 2001; Lee and Park, 2008. plus tree (+)
      IF (alpha_opt == 2) alpha = 0.772*hwr + alphav

      ! Macdonald, 2000 plus tree (+)
      IF (alpha_opt == 3) alpha = 9.6*fai + alphav

!-----------------------------------------------------------------------
! first guess for taf and qaf for each layer
! a large differece from previous schemes
!-----------------------------------------------------------------------
      taf(:) = 0.
      qaf(:) = 0.

      IF (numlay .eq. 2) THEN
         taf(3) = (tg + 2.*thm)/3.
         qaf(3) = (qg + 2.*qm )/3.
         taf(2) = (2.*tg + thm)/3.
         qaf(2) = (2.*qg + qm )/3.
      ENDIF

      IF (numlay .eq. 3) THEN
         taf(3) = (tg + 3.*thm)/4.
         qaf(3) = (qg + 3.*qm )/4.
         taf(2) = (tg + thm   )/2.
         qaf(2) = (qg + qm    )/2.
         taf(1) = (3.*tg + thm)/4.
         qaf(1) = (3.*qg + qm )/4.
      ENDIF

!-----------------------------------------------------------------------
! some environment variables
! how to calculate rsoil and what is its usage?
!-----------------------------------------------------------------------
      pco2a = pco2m
      tprcor = 44.6*273.16*psrf/1.013e5
      rsoil = 0.   !respiration (mol m-2 s-1)
      !rsoil = 1.22e-6*exp(308.56*(1./56.02-1./(tg-227.13)))
      !rsoil = rstfac * 0.23 * 15. * 2.**((tg-273.16-10.)/10.) * 1.e-6
      !rsoil = 5.22 * 1.e-6
      rsoil = 0.22 * 1.e-6

! initialization and input values for Monin-Obukhov
      ! have been set before
      z0h = z0m; z0q = z0m
      ur  = max(0.1, sqrt(us*us+vs*vs))    !limit set to 0.1
      dth = thm - taf(2)
      dqh =  qm - qaf(2)
      dthv = dth*(1.+0.61*qm) + 0.61*th*dqh

      ! To ensure the obs height >= hroof+10.
      huu = max(hroof+10., hu)
      htu = max(hroof+10., ht)
      hqu = max(hroof+10., hq)

      zldis = huu - displa

      IF (zldis <= 0.0) THEN
         write(6,*) 'the obs height of u less than the zero displacement heght'
         CALL abort
      ENDIF

      CALL moninobukini(ur,th,thm,thv,dth,dqh,dthv,zldis,z0m,um,obu)

! ======================================================================
!     BEGIN stability iteration
! ======================================================================

      DO WHILE (it .le. itmax)

         tlbef = tl

         del2  = del
         dele2 = dele

!-----------------------------------------------------------------------
! Aerodynamical resistances
!-----------------------------------------------------------------------
! Evaluate stability-dependent variables using moz from prior iteration

         CALL moninobukm(huu,htu,hqu,displa,z0m,z0h,z0q,obu,um, &
            hroof,0.,ustar,fh2m,fq2m,hroof,fmtop,fm,fh,fq,fht,fqt,phih)

! Aerodynamic resistance
         ! 09/16/2017:
         ! note that for ram, it is the resistance from Href to z0m+displa
         ! however, for rah and raw is only from Href to canopy effective
         ! exchange height.
         ! so rah/raw is not comparable with that of 1D case
         ram = 1./(ustar*ustar/um)

         ! 05/02/2016: calculate resistance from the top layer (effective exchange
         ! height) to reference height
         ! for urban, from roof height to reference height
         rah = 1./(vonkar/(fh-fht)*ustar)
         raw = 1./(vonkar/(fq-fqt)*ustar)

! update roughness length for sensible/latent heat
         z0hg = z0mg/exp(0.13 * (ustar*z0mg/1.5e-5)**0.45)
         z0qg = z0hg

         z0h = max(z0hg, z0h)
         z0q = max(z0qg, z0q)

!-----------------------------------------------------------------------
! new method to calculate rd and ueffect
! the kernel part of 3d model
!-----------------------------------------------------------------------

         ! initialization
         rd(:)  = 0.
         rd_(:) = 0.
         ueff_lay(:)  = 0.
         ueff_lay_(:) = 0.

         ! calculate canopy top wind speed (utop) and exchange coefficient (ktop)
         ! need to update each time as obu changed after each iteration
         ! print*, ustar, fmtop
         utop = ustar/vonkar * fmtop
         ktop = vonkar * (hroof-displa) * ustar / phih

         ueff_lay(3)  = utop
         ueff_lay_(3) = utop

         ! real(r8) FUNCTION kintegral(ktop, fc, bee, alpha, z0mg, &
         !      displah, htop, hbot, obu, ustar, ztop, zbot)
         !rd_(3)  = kintegral(ktop, 1., bee, alpha, z0mg, displa/hroof, &
         !   hroof, 0., obug, ustarg, hroof, displau+z0mu)

         ! real(r8) FUNCTION frd(ktop, htop, hbot, &
         !      ztop, zbot, displah, z0h, obu, ustar, &
         !      z0mg, alpha, bee, fc)
         rd(3) = frd(ktop, hroof, 0., hroof, displau+z0mu, displa/hroof, z0h_g, &
            obug, ustarg, z0mg, alpha, bee, 1.)

         ! real(r8) FUNCTION uintegralz(utop, fc, bee, alpha, z0mg, htop, hbot, ztop, zbot)
         !ueff_lay(2)  = uintegralz(utop, 1., bee, alpha, z0mg, hroof, 0., hroof, z0mg)

         ! real(r8) FUNCTION ueffectz(utop, htop, hbot, &
         !      ztop, zbot, z0mg, alpha, bee, fc)
         ueff_lay(2) = ueffectz(utop, hroof, 0., hroof, z0mg, z0mg, alpha, bee, 1.)

         IF (numlay == 3) THEN
            ! real(r8) FUNCTION kintegral(ktop, fc, bee, alpha, z0mg, &
            !      displah, htop, hbot, obu, ustar, ztop, zbot)
            !rd(2)  = kintegral(ktop, 1., bee, alpha, z0mg, displa/hroof, &
            !   hroof, 0., obug, ustarg, displau+z0mu, displav+z0mv)
            rd(2) = frd(ktop, hroof, 0., displau+z0mu, displav+z0mv, displa/hroof, z0h_g, &
               obug, ustarg, z0mg, alpha, bee, 1.)

            !rd(1)  = kintegral(ktop, 1., bee, alpha, z0mg, displa/hroof, &
            !   hroof, 0., obug, ustarg, displav+z0mv, z0qg)
            rd(1) = frd(ktop, hroof, 0., displav+z0mv, z0qg, displa/hroof, z0h_g, &
               obug, ustarg, z0mg, alpha, bee, 1.)

            ! calculate ra2m, rd2m
            ra2m = frd(ktop, hroof, 0., displav+z0mv, 2., displa/hroof, z0h_g, &
               obug, ustarg, z0mg, alpha, bee, 1.)

            rd2m = frd(ktop, hroof, 0., 2., z0qg, displa/hroof, z0h_g, &
               obug, ustarg, z0mg, alpha, bee, 1.)
         ELSE
            !rd_(2)  = kintegral(ktop, 1., bee, alpha, z0mg, displa/hroof, &
            !   hroof, 0., obug, ustarg, displau+z0mu, z0qg)
            rd(2) = frd(ktop, hroof, 0., displau+z0mu, z0qg, displa/hroof, z0h_g, &
               obug, ustarg, z0mg, alpha, bee, 1.)

            ! calculate ra2m, rd2m
            ra2m = frd(ktop, hroof, 0., displau+z0mu, 2., displa/hroof, z0h_g, &
               obug, ustarg, z0mg, alpha, bee, 1.)

            rd2m = frd(ktop, hroof, 0., 2., z0qg, displa/hroof, z0h_g, &
               obug, ustarg, z0mg, alpha, bee, 1.)
         ENDIF

         !ueff_lay(2)  = uintegralz(utop, 1., bee, alpha, z0mg, hroof, 0., hroof, z0mg)
         !print *, "htop/hbot:", htop, hbot  !fordebug
         !ueff_veg  = uintegralz(utop, 1., bee, alpha, z0mg, hroof, 0., htop, hbot)

         !ueff_lay_(2) = ueffectz(utop, hroof, 0., hroof, z0mg, z0mg, alpha, bee, 1.)
         ueff_veg = ueffectz(utop, hroof, 0., htop, hbot, z0mg, alpha, bee, 1.)

         ! Masson, 2000: Account for different canyon orientations
         ! 2/PI is a factor derived from 0-360deg integration
         IF (alpha_opt == 1) THEN
            ueff_lay(2) = 2/PI*ueff_lay(2)
            ueff_veg    = 2/PI*ueff_veg
            rd(:)       = PI/2*rd(:)
         ENDIF

         ! ueff_lay(3) = ueff_lay(2)

         !print *, "ueff_lay :", ueff_lay
         !print *, "ueff_lay_:", ueff_lay_
         !print *, "------------------------"
         !print *, "rd :", rd
         !print *, "rd_:", rd_

!-----------------------------------------------------------------------
! Bulk boundary layer resistance of leaves
!-----------------------------------------------------------------------
         rb(:) = 0.

         DO i = 0, nurb

            IF (i == 3) THEN
               cf = 0.01*sqrtdi*sqrt(ueff_veg)
               rb(i) = 1./cf
               CYCLE
            ENDIF

            clev = canlev(i)
            rb(i) = rhoair * cpair / ( 11.8 + 4.2*ueff_lay(clev) )

            ! Cole & Sturrock (1977) Building and Environment, 12, 207–214.
            ! rb(i) = rhoair * cpair / ( 5.8 + 4.1*ueff_lay(clev) )
            !IF (ueff_lay(clev) > 5.) THEN
            !   rb(i) = rhoair * cpair / (7.51*ueff_lay(clev)**0.78)
            !ELSE
            !   rb(i) = rhoair * cpair / (5.8 + 4.1*ueff_lay(clev))
            !ENDIF
            !rb(i) = rhoair * cpair / (cpair*vonkar*vonkar*ueff_lay(clev)/(log(0.1*hroof/)*(2.3+log(0.1*hroof/))))
         ENDDO

!-----------------------------------------------------------------------
! stomatal resistances
!-----------------------------------------------------------------------

         IF (lai > 0.) THEN

            ! only for vegetation
            ! rb(3) = rb(3)

            clev = canlev(3)
            eah = qaf(clev) * psrf / ( 0.622 + 0.378 * qaf(clev) )    !pa

!-----------------------------------------------------------------------
! note: calculate resistance for leaves
!-----------------------------------------------------------------------
            CALL stomata (vmax25,effcon ,slti   ,hlti   ,&
               shti    ,hhti    ,trda   ,trdm   ,trop   ,&
               g1      ,g0      ,gradm  ,binter ,thm    ,&
               psrf    ,po2m    ,pco2m  ,pco2a  ,eah    ,&
               ei(3)   ,tu(3)   ,par    ,&
               o3coefv ,o3coefg ,&
               rb(3)/lai,raw    ,rstfac ,cint(:),&
               assim   ,respc   ,rs     &
               )
         ELSE
            rs = 2.e4; assim = 0.; respc = 0.
         ENDIF

! above stomatal resistances are for the canopy, the stomatal rsistances
! and the "rb" in the following calculations are the average for single leaf. thus,
         rs = rs * lai

!-----------------------------------------------------------------------
! dimensional and non-dimensional sensible and latent heat conductances
! for canopy and soil flux calculations.
!-----------------------------------------------------------------------

         cfh(:) = 0.
         cfw(:) = 0.

         DO i = 0, nurb

            IF (i == 3) THEN

               clev = canlev(i)
               delta = 0.0
               IF (qsatl(i)-qaf(clev) .gt. 0.) delta = 1.0

               ! calculate sensible heat conductance
               cfh(i) = lsai / rb(i)

               ! for building walls, cfw=0., no water transfer
               ! for canopy, keep the same but for one leaf
               ! calculate latent heat conductance
               cfw(i) = (1.-delta*(1.-fwet))*lsai/rb(i) + &
                  (1.-fwet)*delta* ( lai/(rb(i)+rs) )
            ELSE
               cfh(i) = 1 / rb(i)

               IF (i == 0) THEN !roof
                  ! account for fwet
                  cfw(i) = fwet_roof / rb(i)
               ELSE
                  cfw(i) = 1 / rb(i)
               ENDIF
            ENDIF
         ENDDO

         ! For simplicity, there is no water exchange on the wall
         cfw(1:2) = 0.

         ! initialization
         cah(:) = 0.
         caw(:) = 0.
         cgh(:) = 0.
         cgw(:) = 0.

         ! conductance for each layer
         DO i = 3, botlay, -1
            IF (i == 3) THEN
               cah(i) = 1. / rah
               caw(i) = 1. / raw
            ! ELSE IF (i == 2) THEN
            !    cah(i) = 1e6
            !    caw(i) = 1e6
            ELSE
               cah(i) = 1. / rd(i+1)
               caw(i) = 1. / rd(i+1)
            ENDIF

            ! IF (i == 3) THEN
            !    cgh(i) = 1e6
            !    cgw(i) = 1e6
            ! ELSE
               cgh(i) = 1. / rd(i)
               cgw(i) = 1. / rd(i)
            ! ENDIF
         ENDDO

         ! claculate wtshi, wtsqi
         wtshi(:) = cah(:)*fah(:) + cgh(:)*fgh(:)
         wtsqi(:) = caw(:)*faw(:) + cgw(:)*fgw(:)

         DO i = 0, nurb
            clev = canlev(i)
            wtshi(clev) = wtshi(clev) + fc(i)*cfh(i)
            wtsqi(clev) = wtsqi(clev) + fc(i)*cfw(i)
         ENDDO

         DO i = 3, 3-numlay+1, -1
            wtshi(i) = 1./wtshi(i)
            wtsqi(i) = 1./wtsqi(i)
         ENDDO

         wta0(:) = cah(:) * wtshi(:) * fah(:)
         wtg0(:) = cgh(:) * wtshi(:) * fgh(:)

         wtaq0(:) = caw(:) * wtsqi(:) * faw(:)
         wtgq0(:) = cgw(:) * wtsqi(:) * fgw(:)

         ! calculate wtl0, wtll, wtlq0, wtlql
         wtll(:)  = 0.
         wtlql(:) = 0.

         DO i = 0, nurb
            clev = canlev(i)

            wtl0(i)  = cfh(i) * wtshi(clev) * fc(i)
            wtll(clev) = wtll(clev) + wtl0(i)*tu(i)

            wtlq0(i) = cfw(i) * wtsqi(clev) * fc(i)
            wtlql(clev) = wtlql(clev) + wtlq0(i)*qsatl(i)
         ENDDO

         ! to solve taf(:) and qaf(:)

         IF (numlay .eq. 2) THEN

            ! - Equations:
            ! taf(3) = (1/rah*thm + 1/rd(3)*taf(2) + 1/rb(0)*troof*fc(0) + AHE/(rho*cp))/(1/rah + 1/rd(3) + 1/rb(0)*fc(0))
            ! taf(2) = (1/rd(3)*taf(3) + 1/rd(2)*tg*fg + 1/rb(1)*twsun*fc(1) + 1/rb(2)*twsha*fc(2) + lsai/rb(3)*tl*fc(3) + AHE/(rho*cp))/ &
            !          (1/rd(3) + 1/rd(2)*fg + 1/rb(1)*fc(1) + 1/rb(2)*fc(2) + lsai/rb(3)*fc(3))
            ! Also written as:
            ! taf(3) = (cah(3)*thm + cah(2)*taf(2) + cfh(0)*troof*fc(0))/(cah(3) + cah(2) + cfh(0)*fc(0))
            ! taf(2) = (cah(2)*taf(3) + cgh(2)*tg*fg + cfh(1)*twsun*fc(1) + cfh(2)*twsha*fc(2) + cfh(3)*tl*fc(3) + AHE/(rho*cp))/ &
            !          (cah(2) + cgh(2)*fg + cfh(1)*fc(1) + cfh(2)*fc(2) + cfh(3)*fc(3))
            ! - Equations:
            ! qaf(3) = (1/raw*qm + 1/rd(3)*qaf(2) + 1/rb(0)*qroof*fc(0))/(1/raw + 1/rd(3) + 1/rb(0)*fc(0))
            ! qaf(2) = (1/rd(3)*qaf(3) + 1/(rd(2)+rsr)*qper*fgper*fg + fwetimp/rd(2)*qimp*fgimp*fg + lsai/(rb(3)+rs)*ql*fc(3) + AHE/rho)/ &
            !          (1/rd(3) + 1/(rd(2)+rsr)*fgper*fg + fwetimp/rd(2)*fgimp*fg + lsai/(rb(3)+rs)*fc(3))
            ! Also written as:
            ! qaf(3) = (caw(3)*qm + caw(2)*qaf(2) + cfw(0)*qroof*fc(0))/(caw(3) + caw(2) + cfw(0)*fc(0))
            ! qaf(2) = (caw(2)*qaf(3) + cgwper*qper*fgper*fg + cgwimp*qimp*fgimp*fg + cfw(3)*ql*fc(3) + AHE/rho)/ &
            !          (caw(2) + cgwper*fgper*fg + cgwimp*fgimp*fg + cfw(3)*fc(3))

            ! 06/20/2021, yuan: account for Anthropogenic heat
            ! 92% heat release as SH, Pigeon et al., 2007

            h_vec  = vehc!
            tmpw1  = cah(2)*((cah(3)*thm + cfh(0)*tu(0)*fc(0) + 1/(4*hlr+1)*(Fhac+Fwst)/(rhoair*cpair))/&
                     (cah(3) + cah(2) + cfh(0)*fc(0)))
            tmpw2  = (4*hlr/(4*hlr+1)*(Fhac+Fwst)+Fach)/(rhoair*cpair) + (h_vec+meta)/(rhoair*cpair)
            tmpw3  = cgh(2)*fg*tg + cfh(1)*tu(1)*fc(1) + cfh(2)*tu(2)*fc(2) + cfh(3)*tu(3)*fc(3)
            fact   = 1. - (cah(2)*cah(2)/(cah(3) + cah(2) + cfh(0)*fc(0))/&
                     (cah(2) + cgh(2)*fg + cfh(1)*fc(1) + cfh(2)*fc(2) + cfh(3)*fc(3)))
            taf(2) = (tmpw1 + tmpw2 + tmpw3) / &
                     (cah(2) + cgh(2)*fg + cfh(1)*fc(1) + cfh(2)*fc(2) + cfh(3)*fc(3)) / &
                     fact

            IF (qgper < qaf(2)) THEN
              ! dew case. no soil resistance
              cgw_per= cgw(2)
            ELSE
              cgw_per= 1/(1/cgw(2)+rsr)
            ENDIF

            cgw_imp= fwet_gimp*cgw(2)

            ! account for soil resistance, qgper and qgimp are calculated separately
            l_vec  = 0
            tmpw1  = caw(2)*((caw(3)*qm + cfw(0)*qsatl(0)*fc(0))/&
                     (caw(3) + caw(2) + cfw(0)*fc(0)))
            tmpw2  = l_vec/(rhoair)
            tmpw3  = cgw_per*qgper*fgper*fg + cgw_imp*qgimp*fgimp*fg + cfw(3)*qsatl(3)*fc(3)
            facq   = 1. - (caw(2)*caw(2)/&
                     (caw(3) + caw(2) + cfw(0)*fc(0))/&
                     (caw(2) + cgw_per*fgper*fg + cgw_imp*fgimp*fg + cfw(3)*fc(3)))
            qaf(2) = (tmpw1 + tmpw2 + tmpw3)/&
                     (caw(2) + cgw_per*fgper*fg + cgw_imp*fgimp*fg + cfw(3)*fc(3))/&
                     facq

            tmpw1  = 1/(4*hlr+1)*(Fhac+Fwst)/(rhoair*cpair)
            taf(3) = (cah(3)*thm + cah(2)*taf(2) + cfh(0)*tu(0)*fc(0) + tmpw1)/&
                     (cah(3) + cah(2) + cfh(0)*fc(0))
            qaf(3) = (caw(3)*qm  + caw(2)*qaf(2) + cfw(0)*qsatl(0)*fc(0))/&
                     (caw(3) + caw(2) + cfw(0)*fc(0))
         ENDIF

         IF (numlay .eq. 3) THEN

            ! - Equations:
            ! taf(3) = (thm/rah+1/rd(3)*taf(2)+AHE2/rho/cpair+1/rb(0)*troof*fc(0))/&
            !          (1/rah+1/rd(3)+1/rb(0)*fc(0))
            ! taf(2) = (1/rd(3)*taf(3)+1/rd(2)*taf1+1/rb(1)*twsun*fc(1)+1/rb(2)*twsha*fc(2)+AHE1/rho/cpair)/&
            !          (1/rd(3)+1/rd(2)+1/rb(1)*fc(1)+1/rb(2)*fc(2))
            ! taf(1) = (1/rd(2)*taf(2)+1/rd(1)*tg*fg+1/rb(3)*tl*fc(3)+Hveh/rhoair/cpair)/&
            !          (1/rd(2)+1/rd(1)*fg+1/rb(3)*fc(3))
            ! - Equations:
            ! qaf(3) = (1/raw*qm+1/rd(3)*qaf(2)+1/rb(0)*qroof*fc(0))/&
            !          (1/raw+1/rd(3)+1/rb(0)*fc(0))
            ! qaf(2) = (1/rd(3)*qaf(3)+1/rd(2)*qaf(1))/&
            !          (1/rd(3) + 1/rd(2))
            ! qaf(1) = (1/rd(2)*qaf(2)+1/(rd(1)+rsr)*qgper*fgper*fg+1/rd(1)*qimp*fgimp*fg+1/(rb(3)+rs)*ql*fc(3)+h_veh/rho))/&
            !          (1/rd(2)+1/(rd(1)+rsr)*fgper*fg+1/rd(1)*fgimp*fg+1/(rb(3)+rs)*fc(3))

            tmpw1  = cah(1)*(cgh(1)*tg*fg + cfh(3)*tu(3)*fc(3) + (vehc+meta)/rhoair/cpair)/&
                     (cah(1) + cgh(1)*fg + cfh(3)*fc(3))
            tmpw2  = cah(2)*(cah(3)*thm + cfh(0)*tu(0)*fc(0) + 1/(4*hlr+1)*(Fhac+Fwst)/(rhoair*cpair))/&
                     (cah(3) + cah(2) + cfh(0)*fc(0))
            tmpw3  = cah(1)*cah(1)/&
                     (cah(1) + cgh(1)*fg + cfh(3)*fc(3))/&
                     (cah(1) + cah(2) + cfh(1)*fc(1) + cfh(2)*fc(2))
            tmpw4  = cah(2)*cah(2)/&
                     (cah(3) + cah(2) + cfh(0)*fc(0))/&
                     (cah(1) + cah(2) + cfh(1)*fc(1) + cfh(2)*fc(2))
            fact   = 1. - tmpw3 - tmpw4

            taf(2) = (tmpw1 + tmpw2 + cfh(1)*tu(1)*fc(1) + cfh(2)*tu(2)*fc(2) + (4*hlr/(4*hlr+1)*(Fhac+Fwst)+Fach)/(rhoair*cpair))/&
                     (cah(1) + cah(2) + cfh(1)*fc(1) + cfh(2)*fc(2))/&
                     fact

            taf(1) = (cah(1)*taf(2) + cgh(1)*tg*fg + cfh(3)*tu(3)*fc(3) + (vehc+meta)/rhoair/cpair)/&
                     (cah(1) + cgh(1)*fg + cfh(3)*fc(3))
            tmpw1  = 1/(4*hlr+1)*(Fhac+Fwst)/(rhoair*cpair)
            taf(3) = (cah(3)*thm + cah(2)*taf(2) + cfh(0)*tu(0)*fc(0) + tmpw1)/&
                     (cah(3) + cah(2) + cfh(0)*fc(0))

            IF (qgper < qaf(1)) THEN
              ! dew case. no soil resistance
              cgw_per= cgw(1)
            ELSE
              cgw_per= 1/(1/cgw(1)+rsr)
            ENDIF

            cgw_imp= fwet_gimp*cgw(1)

            l_vec   = 0
            tmpw1  = caw(1)*(cgw_per*qgper*fgper*fg + cgw_imp*qgimp*fgimp*fg + cfw(3)*qsatl(3)*fc(3) + l_vec/(rhoair))/&
                     (caw(1) + cgw_per*fgper*fg + cgw_imp*fgimp*fg + cfw(3)*fc(3))
            tmpw2  = caw(2)*(caw(3)*qm + cfw(0)*qsatl(0)*fc(0))/&
                     (caw(3) + caw(2) + cfw(0)*fc(0))
            tmpw3  = caw(1)*caw(1)/&
                     (caw(1) + cgw_per*fgper*fg + cgw_imp*fgimp*fg + cfw(3)*fc(3))/&
                     (caw(2) + caw(1))
            tmpw4  = caw(2)*caw(2)/&
                     (caw(3) + caw(2) + cfw(0)*fc(0))/&
                     (caw(2) + caw(1))
            facq   = 1. - tmpw3 - tmpw4

            qaf(2) = (tmpw1 + tmpw2)/&
                     (caw(2) + caw(1))/&
                     facq

            tmpw1  = l_vec/(rhoair)
            qaf(1) = (caw(1)*qaf(2) + qgper*cgw_per*fgper*fg + qgimp*cgw_imp*fgimp*fg + cfw(3)*qsatl(3)*fc(3) + tmpw1)/&
                     (caw(1) + cgw_per*fgper*fg + cgw_imp*fgimp*fg + cfw(3)*fc(3))
            qaf(3) = (caw(3)*qm + caw(2)*qaf(2) + cfw(0)*qsatl(0)*fc(0))/&
                     (caw(3) + caw(2) + cfw(0)*fc(0))
         ENDIF

!-----------------------------------------------------------------------
! IR radiation, sensible and latent heat fluxes and their derivatives
!-----------------------------------------------------------------------
! the partial derivatives of areodynamical resistance are ignored
! which cannot be determined analtically

         !NOTE: ONLY for vegetation
         i = 3

! sensible heat fluxes and their derivatives
         fsenl = rhoair * cpair * cfh(i) * (tl - taf(botlay))

         ! 09/24/2017: why fact/facq here? bugs? YES
         ! 09/25/2017: re-written, check it clearfully
         ! 11/25/2021: re-written, double check
         IF (botlay == 2) THEN
            fsenl_dtl = rhoair * cpair * cfh(i) * (1.-wtl0(i)/fact)
         ELSE
            fsenl_dtl = rhoair * cpair * cfh(i) * (1.-wta0(1)*wtg0(2)*wtl0(i)/fact-wtl0(i))
         ENDIF


! latent heat fluxes and their derivatives
         etr = rhoair * (1.-fwet) * delta * lai/(rb(i)+rs) &
             * (qsatl(i) - qaf(botlay))

         IF (botlay == 2) THEN
            etr_dtl = rhoair * (1.-fwet) * delta * lai/(rb(i)+rs) &
                    * (1.-wtlq0(i)/facq)*qsatldT(i)
         ELSE
            etr_dtl = rhoair * (1.-fwet) * delta * lai/(rb(i)+rs) &
                    * (1.-wtaq0(1)*wtgq0(2)*wtlq0(i)/facq-wtlq0(i))*qsatldT(i)
         ENDIF

         IF (etr.ge.etrc) THEN
            etr = etrc
            etr_dtl = 0.
         ENDIF

         evplwet = rhoair * (1.-delta*(1.-fwet)) * lsai/rb(i) &
                 * (qsatl(i) - qaf(botlay))

         IF (botlay == 2) THEN
            evplwet_dtl = rhoair * (1.-delta*(1.-fwet)) * lsai/rb(i) &
                        * (1.-wtlq0(i)/facq)*qsatldT(i)
         ELSE
            evplwet_dtl = rhoair * (1.-delta*(1.-fwet)) * lsai/rb(i) &
                        * (1.-wtaq0(1)*wtgq0(2)*wtlq0(i)/facq-wtlq0(i))*qsatldT(i)
         ENDIF

         IF (evplwet.ge.ldew/deltim) THEN
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

         ! calculate irab, dirab_dtl
         B(5)    = B_5*tl**4
         B1(5)   = B1_5*tl**4
         dBdT(5) = dBdT_5*tl**3
         X  = matmul(Ainv, B)
         ! first 5 items of dBdT is 0, dBdT*(0,0,0,0,0,1)
         dX = matmul(Ainv, dBdT*uvec)

         ! calculate longwave for vegetation
         irab = ( (sum(X(1:4)*VegVF(1:4)) + frl*VegVF(5))*ev - B1(5))/fcover(5)*fg
         dirab_dtl = ( sum(dX(1:4)*VegVF(1:4))*ev - dBdT(5) )/fcover(5)*fg

         ! solve for leaf temperature
         dtl(it) = (sabv + irab - fsenl - hvap*fevpl) &
            / (lsai*clai/deltim - dirab_dtl + fsenl_dtl + hvap*fevpl_dtl)
         dtl_noadj = dtl(it)

         ! check magnitude of change in leaf temperature limit to maximum allowed value

         IF (it .le. itmax) THEN

            ! put brakes on large temperature excursions
            IF (abs(dtl(it)).gt.delmax) THEN
               dtl(it) = delmax*dtl(it)/abs(dtl(it))
            ENDIF

            IF ((it.ge.2) .and. (dtl(it-1)*dtl(it).le.0.)) THEN
               dtl(it) = 0.5*(dtl(it-1) + dtl(it))
            ENDIF

         ENDIF

         tl = tlbef + dtl(it)
         tu(3) = tl

!-----------------------------------------------------------------------
! square roots differences of temperatures and fluxes for use as the condition of convergences
!-----------------------------------------------------------------------

         del  = sqrt( dtl(it)*dtl(it) )
         dele = dtl(it) * dtl(it) * &
            ( dirab_dtl**2 + fsenl_dtl**2 + hvap*fevpl_dtl**2 )
         dele = sqrt(dele)

!-----------------------------------------------------------------------
!  saturated vapor pressures and canopy air temperature, canopy air humidity
!-----------------------------------------------------------------------
! Recalculate leaf saturated vapor pressure (ei_)for updated leaf temperature
! and adjust specific humidity (qsatl_) proportionately
         CALL qsadv(tu(i),psrf,ei(i),deiDT(i),qsatl(i),qsatldT(i))

! update vegetation/ground surface temperature, canopy air temperature,
! canopy air humidity

         ! calculate wtll, wtlql
         wtll(:)  = 0.
         wtlql(:) = 0.

         DO i = 0, nurb
            clev = canlev(i)
            wtll(clev)  =  wtll(clev) + wtl0(i)*tu(i)
            wtlql(clev) = wtlql(clev) + wtlq0(i)*qsatl(i)
         ENDDO

         IF (numlay .eq. 2) THEN

            ! - Equations:
            ! taf(3) = (1/rah*thm + 1/rd(3)*taf(2) + 1/rb(0)*troof*fc(0) + AHE/(rho*cp))/(1/rah + 1/rd(3) + 1/rb(0)*fc(0))
            ! taf(2) = (1/rd(3)*taf(3) + 1/rd(2)*tg*fg + 1/rb(1)*twsun*fc(1) + 1/rb(2)*twsha*fc(2) + lsai/rb(3)*tl*fc(3) + AHE/(rho*cp))/ &
            !          (1/rd(3) + 1/rd(2)*fg + 1/rb(1)*fc(1) + 1/rb(2)*fc(2) + lsai/rb(3)*fc(3))
            ! Also written as:
            ! taf(3) = (cah(3)*thm + cah(2)*taf(2) + cfh(0)*troof*fc(0))/(cah(3) + cah(2) + cfh(0)*fc(0))
            ! taf(2) = (cah(2)*taf(3) + cgh(2)*tg*fg + cfh(1)*twsun*fc(1) + cfh(2)*twsha*fc(2) + cfh(3)*tl*fc(3) + AHE/(rho*cp))/ &
            !          (cah(2) + cgh(2)*fg + cfh(1)*fc(1) + cfh(2)*fc(2) + cfh(3)*fc(3))
            ! - Equations:
            ! qaf(3) = (1/raw*qm + 1/rd(3)*qaf(2) + 1/rb(0)*qroof*fc(0))/(1/raw + 1/rd(3) + 1/rb(0)*fc(0))
            ! qaf(2) = (1/rd(3)*qaf(3) + 1/(rd(2)+rsr)*qper*fgper*fg + fwetimp/rd(2)*qimp*fgimp*fg + lsai/(rb(3)+rs)*ql*fc(3) + AHE/rho)/ &
            !          (1/rd(3) + 1/(rd(2)+rsr)*fgper*fg + fwetimp/rd(2)*fgimp*fg + lsai/(rb(3)+rs)*fc(3))
            ! Also written as:
            ! qaf(3) = (caw(3)*qm + caw(2)*qaf(2) + cfw(0)*qroof*fc(0))/(caw(3) + caw(2) + cfw(0)*fc(0))
            ! qaf(2) = (caw(2)*qaf(3) + cgwper*qper*fgper*fg + cgwimp*qimp*fgimp*fg + cfw(3)*ql*fc(3) + AHE/rho)/ &
            !          (caw(2) + cgwper*fgper*fg + cgwimp*fgimp*fg + cfw(3)*fc(3))

            ! 06/20/2021, yuan: account for AH
            ! 92% heat release as SH, Pigeon et al., 2007

            h_vec  = vehc
            tmpw1  = cah(2)*((cah(3)*thm + cfh(0)*tu(0)*fc(0) + 1/(4*hlr+1)*(Fhac+Fwst)/(rhoair*cpair))/&
                     (cah(3) + cah(2) + cfh(0)*fc(0)))
            tmpw2  = (4*hlr/(4*hlr+1)*(Fhac+Fwst)+Fach)/(rhoair*cpair) + (h_vec+meta)/(rhoair*cpair)
            tmpw3  = cgh(2)*fg*tg + cfh(1)*tu(1)*fc(1) + cfh(2)*tu(2)*fc(2) + cfh(3)*tu(3)*fc(3)
            fact   = 1. - (cah(2)*cah(2)/(cah(3) + cah(2) + cfh(0)*fc(0))/&
                     (cah(2) + cgh(2)*fg + cfh(1)*fc(1) + cfh(2)*fc(2) + cfh(3)*fc(3)))
            taf(2) = (tmpw1 + tmpw2 + tmpw3) / &
                     (cah(2) + cgh(2)*fg + cfh(1)*fc(1) + cfh(2)*fc(2) + cfh(3)*fc(3)) / &
                     fact

            IF (qgper < qaf(2)) THEN
              ! dew case. no soil resistance
              cgw_per= cgw(2)
            ELSE
              cgw_per= 1/(1/cgw(2)+rsr)
            ENDIF

            cgw_imp= fwet_gimp*cgw(2)

            ! account for soil resistance, qgper and qgimp are calculated separately
            l_vec  = 0
            tmpw1  = caw(2)*((caw(3)*qm + cfw(0)*qsatl(0)*fc(0))/&
                     (caw(3) + caw(2) + cfw(0)*fc(0)))
            tmpw2  = l_vec/(rhoair)
            tmpw3  = cgw_per*qgper*fgper*fg + cgw_imp*qgimp*fgimp*fg + cfw(3)*qsatl(3)*fc(3)
            facq   = 1. - (caw(2)*caw(2)/&
                     (caw(3) + caw(2) + cfw(0)*fc(0))/&
                     (caw(2) + cgw_per*fgper*fg + cgw_imp*fgimp*fg + cfw(3)*fc(3)))
            qaf(2) = (tmpw1 + tmpw2 + tmpw3)/&
                     (caw(2) + cgw_per*fgper*fg + cgw_imp*fgimp*fg + cfw(3)*fc(3))/&
                     facq

            tmpw1  = 1/(4*hlr+1)*(Fhac+Fwst)/(rhoair*cpair)
            taf(3) = (cah(3)*thm + cah(2)*taf(2) + cfh(0)*tu(0)*fc(0) + tmpw1)/&
                     (cah(3) + cah(2) + cfh(0)*fc(0))
            qaf(3) = (caw(3)*qm  + caw(2)*qaf(2) + cfw(0)*qsatl(0)*fc(0))/&
                     (caw(3) + caw(2) + cfw(0)*fc(0))
         ENDIF

         IF (numlay .eq. 3) THEN

            ! - Equations:
            ! taf(3) = (thm/rah+1/rd(3)*taf(2)+AHE2/rho/cpair+1/rb(0)*troof*fc(0))/&
            !          (1/rah+1/rd(3)+1/rb(0)*fc(0))
            ! taf(2) = (1/rd(3)*taf(3)+1/rd(2)*taf1+1/rb(1)*twsun*fc(1)+1/rb(2)*twsha*fc(2)+AHE1/rho/cpair)/&
            !          (1/rd(3)+1/rd(2)+1/rb(1)*fc(1)+1/rb(2)*fc(2))
            ! taf(1) = (1/rd(2)*taf(2)+1/rd(1)*tg*fg+1/rb(3)*tl*fc(3)+Hveh/rhoair/cpair)/&
            !          (1/rd(2)+1/rd(1)*fg+1/rb(3)*fc(3))
            ! - Equations:
            ! qaf(3) = (1/raw*qm+1/rd(3)*qaf(2)+1/rb(0)*qroof*fc(0))/&
            !          (1/raw+1/rd(3)+1/rb(0)*fc(0))
            ! qaf(2) = (1/rd(3)*qaf(3)+1/rd(2)*qaf(1))/&
            !          (1/rd(3) + 1/rd(2))
            ! qaf(1) = (1/rd(2)*qaf(2)+1/(rd(1)+rsr)*qgper*fgper*fg+1/rd(1)*qimp*fgimp*fg+1/(rb(3)+rs)*ql*fc(3)+h_veh/rho))/&
            !          (1/rd(2)+1/(rd(1)+rsr)*fgper*fg+1/rd(1)*fgimp*fg+1/(rb(3)+rs)*fc(3))

            tmpw1  = cah(1)*(cgh(1)*tg*fg + cfh(3)*tu(3)*fc(3) + (vehc+meta)/rhoair/cpair)/&
                     (cah(1) + cgh(1)*fg + cfh(3)*fc(3))
            tmpw2  = cah(2)*(cah(3)*thm + cfh(0)*tu(0)*fc(0) + 1/(4*hlr+1)*(Fhac+Fwst)/(rhoair*cpair))/&
                     (cah(3) + cah(2) + cfh(0)*fc(0))
            tmpw3  = cah(1)*cah(1)/&
                     (cah(1) + cgh(1)*fg + cfh(3)*fc(3))/&
                     (cah(1) + cah(2) + cfh(1)*fc(1) + cfh(2)*fc(2))
            tmpw4  = cah(2)*cah(2)/&
                     (cah(3) + cah(2) + cfh(0)*fc(0))/&
                     (cah(1) + cah(2) + cfh(1)*fc(1) + cfh(2)*fc(2))
            fact   = 1. - tmpw3 - tmpw4

            taf(2) = (tmpw1 + tmpw2 + cfh(1)*tu(1)*fc(1) + cfh(2)*tu(2)*fc(2) + (4*hlr/(4*hlr+1)*(Fhac+Fwst)+Fach)/(rhoair*cpair))/&
                     (cah(1) + cah(2) + cfh(1)*fc(1) + cfh(2)*fc(2))/&
                     fact

            taf(1) = (cah(1)*taf(2) + cgh(1)*tg*fg + cfh(3)*tu(3)*fc(3) + (vehc+meta)/rhoair/cpair)/&
                     (cah(1) + cgh(1)*fg + cfh(3)*fc(3))
            tmpw1  = 1/(4*hlr+1)*(Fhac+Fwst)/(rhoair*cpair)
            taf(3) = (cah(3)*thm + cah(2)*taf(2) + cfh(0)*tu(0)*fc(0) + tmpw1)/&
                     (cah(3) + cah(2) + cfh(0)*fc(0))

            IF (qgper < qaf(1)) THEN
              ! dew case. no soil resistance
              cgw_per= cgw(1)
            ELSE
              cgw_per= 1/(1/cgw(1)+rsr)
            ENDIF

            cgw_imp= fwet_gimp*cgw(1)

            l_vec   = 0!vehc*0.08
            tmpw1  = caw(1)*(cgw_per*qgper*fgper*fg + cgw_imp*qgimp*fgimp*fg + cfw(3)*qsatl(3)*fc(3) + l_vec/(rhoair))/&
                     (caw(1) + cgw_per*fgper*fg + cgw_imp*fgimp*fg + cfw(3)*fc(3))
            tmpw2  = caw(2)*(caw(3)*qm + cfw(0)*qsatl(0)*fc(0))/&
                     (caw(3) + caw(2) + cfw(0)*fc(0))
            tmpw3  = caw(1)*caw(1)/&
                     (caw(1) + cgw_per*fgper*fg + cgw_imp*fgimp*fg + cfw(3)*fc(3))/&
                     (caw(2) + caw(1))
            tmpw4  = caw(2)*caw(2)/&
                     (caw(3) + caw(2) + cfw(0)*fc(0))/&
                     (caw(2) + caw(1))
            facq   = 1. - tmpw3 - tmpw4

            qaf(2) = (tmpw1 + tmpw2)/&
                     (caw(2) + caw(1))/&
                     facq

            tmpw1  = l_vec/(rhoair)
            qaf(1) = (caw(1)*qaf(2) + qgper*cgw_per*fgper*fg + qgimp*cgw_imp*fgimp*fg + cfw(3)*qsatl(3)*fc(3) + tmpw1)/&
                     (caw(1) + cgw_per*fgper*fg + cgw_imp*fgimp*fg + cfw(3)*fc(3))
            qaf(3) = (caw(3)*qm + caw(2)*qaf(2) + cfw(0)*qsatl(0)*fc(0))/&
                     (caw(3) + caw(2) + cfw(0)*fc(0))

         ENDIF

         !------------------------------------------------
         ! account for fwet for roof and impervious ground
         IF (qaf(3) > qroof) THEN
            fwet_roof = 1. !dew case
         ELSE
            fwet_roof = fwet_roof_
         ENDIF

         IF (qaf(botlay) > qgimp) THEN
            fwet_gimp = 1. !dew case
         ELSE
            fwet_gimp = fwet_gimp_
         ENDIF

         ! weighted qg
         ! NOTE: IF fwet_gimp=1, same as previous
         fwetfac = fgimp*fwet_gimp + fgper
         qg = (qgimp*fgimp*fwet_gimp + qgper*fgper) / fwetfac

         fgw(2) = fg*fwetfac

! update co2 partial pressure within canopy air
         ! 05/02/2016: may have some problem with gdh2o, however,
         ! this variable seems never used here. Different height
         ! level vegetation should have different gdh2o, i.e.,
         ! different rd(layer) values.
         gah2o = 1.0/raw * tprcor/thm                     !mol m-2 s-1
         gdh2o = 1.0/rd(botlay) * tprcor/thm              !mol m-2 s-1

         pco2a = pco2m - 1.37*psrf/max(0.446,gah2o) * &
            (assim - respc - rsoil)

!-----------------------------------------------------------------------
! Update monin-obukhov length and wind speed including the stability effect
!-----------------------------------------------------------------------

         ! USE the top layer taf and qaf
         !TODO: need more check
         dth = thm - taf(2)
         dqh =  qm - qaf(2)

         tstar = vonkar/(fh)*dth
         qstar = vonkar/(fq)*dqh

         thvstar = tstar*(1.+0.61*qm)+0.61*th*qstar
         zeta = zldis*vonkar*grav*thvstar / (ustar**2*thv)
         IF (zeta .ge. 0.) THEN                             !stable
            zeta = min(2.,max(zeta,1.e-6))
         ELSE                                             !unstable
            zeta = max(-100.,min(zeta,-1.e-6))
         ENDIF
         obu = zldis/zeta

         IF (zeta .ge. 0.) THEN
            um = max(ur,.1)
         ELSE
            wc = (-grav*ustar*thvstar*zii/thv)**(1./3.)
            wc2 = beta*beta*(wc*wc)
            um = sqrt(ur*ur+wc2)
         ENDIF

         IF (obuold*obu .lt. 0.) nmozsgn = nmozsgn+1
         IF (nmozsgn .ge. 4) obu = zldis/(-0.01)
         obuold = obu

!-----------------------------------------------------------------------
! Test for convergence
!-----------------------------------------------------------------------

         it = it+1

         IF (it .gt. itmin) THEN
            fevpl_bef = fevpl
            det = max(del,del2)
            dee = max(dele,dele2)
            IF (det .lt. dtmin .and. dee .lt. dlemin) EXIT
         ENDIF

      ENDDO

! ======================================================================
!     END stability iteration
! ======================================================================

      zol = zeta
      rib = min(5.,zol*ustar**2/(vonkar**2/fh*um**2))

! canopy fluxes and total assimilation amd respiration

      IF (lai .gt. 0.001) THEN
         rst = rs/lai
      ELSE
         rs    = 2.0e4
         assim = 0.
         respc = 0.
         rst   = 2.0e4
      ENDIF
      respc = respc + rsoil

! canopy fluxes and total assimilation amd respiration

      fsenl = fsenl + fsenl_dtl*dtl(it-1) &
         ! add the imbalanced energy below due to T adjustment to sensibel heat
         + (dtl_noadj-dtl(it-1)) * (lsai*clai/deltim - dirab_dtl &
         + fsenl_dtl + hvap*fevpl_dtl) &
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

!-----------------------------------------------------------------------
! Update dew accumulation (kg/m2)
!-----------------------------------------------------------------------

      ldew = max(0., ldew-evplwet*deltim)

!-----------------------------------------------------------------------
! balance check
!-----------------------------------------------------------------------

      err = sabv + irab + dirab_dtl*dtl(it-1) &
          - fsenl - hvap*fevpl

#if(defined CLMDEBUG)
      IF (abs(err) .gt. .2) THEN
         write(6,*) 'energy imbalance in UrbanVegFlux.F90', &
         i,it-1,err,sabv,irab,fsenl,hvap*fevpl
         CALL CoLM_stop()
      ENDIF
#endif


      ! calculate longwave absorption
      lwsun = ( ewall*X(1) - B1(1) ) / (1-ewall)
      lwsha = ( ewall*X(2) - B1(2) ) / (1-ewall)
      lgimp = ( egimp*X(3) - B1(3) ) / (1-egimp)
      lgper = ( egper*X(4) - B1(4) ) / (1-egper)
      lveg  = ( (sum(X(1:4)*VegVF(1:4)) + frl*VegVF(5))*ev - B1(5) )
      lout  = sum( X * SkyVF )

      ! longwave absorption due to leaf temperature change
      lwsun = lwsun + ( ewall*dX(1) ) / (1-ewall) * dtl(it-1)
      lwsha = lwsha + ( ewall*dX(2) ) / (1-ewall) * dtl(it-1)
      lgimp = lgimp + ( egimp*dX(3) ) / (1-egimp) * dtl(it-1)
      lgper = lgper + ( egper*dX(4) ) / (1-egper) * dtl(it-1)
      lveg  = lveg  + ( sum(dX(1:4)*VegVF(1:4))*ev - dBdT(5) ) * dtl(it-1)
      lout  = lout  + sum( dX * SkyVF * dtl(it-1) )

      ! Energy balance check
      err = lwsun + lwsha + lgimp + lgper + lveg + lout

      IF (abs(err-frl) > 1e-6) THEN
         print *, "Longwave - Energy Balance Check error!", err-frl
      ENDIF

      ! convert to per unit area
      IF (fcover(1) > 0.) lwsun = lwsun / fcover(1) * fg !/ (4*fwsun*HL*fb/fg)
      IF (fcover(2) > 0.) lwsha = lwsha / fcover(2) * fg !/ (4*fwsha*HL*fb/fg)
      IF (fcover(3) > 0.) lgimp = lgimp / fcover(3) * fg !/ fgimp
      IF (fcover(4) > 0.) lgper = lgper / fcover(4) * fg !/ fgper
      IF (fcover(5) > 0.) lveg  = lveg  / fcover(5) * fg !/ fv/fg

      ! add previous longwave
      lwsun = lwsun + lwsun_bef
      lwsha = lwsha + lwsha_bef
      lgimp = lgimp + lgimp_bef
      lgper = lgper + lgper_bef
      lveg  = lveg  + lveg_bef

      tafu = taf(2)

!-----------------------------------------------------------------------
! wind stresses
!-----------------------------------------------------------------------

      taux = - rhoair*us/ram
      tauy = - rhoair*vs/ram

!-----------------------------------------------------------------------
! fluxes from roof, walls to canopy space
!-----------------------------------------------------------------------

      ! sensible heat fluxes
      fsenroof = rhoair*cpair*cfh(0)*(troof-taf(3))
      fsenwsun = rhoair*cpair*cfh(1)*(twsun-taf(2))
      fsenwsha = rhoair*cpair*cfh(2)*(twsha-taf(2))

      ! latent heat fluxes
      fevproof = rhoair*cfw(0)*(qsatl(0)-qaf(3))
      fevproof = fevproof*fwet_roof

      croofs = rhoair*cpair*cfh(0)*(1.-wtg0(3)*wta0(2)*wtl0(0)/fact-wtl0(0))
      cwalls = rhoair*cpair*cfh(1)*(1.-wtl0(1)/fact)
      ! croofl = rhoair*cfw(0)*(1.-wtgq0(3)*wtaq0(2)*wtlq0(0)/facq-wtlq0(0))*qsatldT(0)
      croofl = rhoair*cfw(0)*(1.-cfw(0)*fc(0)/(caw(3)+cgw(3)+cfw(0)*fc(0))-cgw(3) &
               /(caw(3)+cgw(3)+cfw(0)*fc(0)) &
               /(cgw(3)+cgw_per*fgper*fg+cgw_imp*fgimp*fg+cfw(3)*fc(3))* &
               cfw(0)*fc(0)*cgw(3)/(caw(3)+cgw(3)+cfw(0)*fc(0))/facq)*qsatldT(0)
      croofl = croofl*fwet_roof

      croof = croofs + croofl*htvp_roof

!-----------------------------------------------------------------------
! fluxes from urban ground to canopy space
!-----------------------------------------------------------------------

      fsengimp = cpair*rhoair*cgh(botlay)*(tgimp-taf(botlay))
      fsengper = cpair*rhoair*cgh(botlay)*(tgper-taf(botlay))

      fevpgimp = rhoair*cgw_imp*(qgimp-qaf(botlay))
      fevpgper = rhoair*cgw_per*(qgper-qaf(botlay))

      fevpgimp = fevpgimp*fwet_gimp

!-----------------------------------------------------------------------
! Derivative of soil energy flux with respect to soil temperature
!-----------------------------------------------------------------------

      IF (botlay == 2) THEN
         cgrnds = cpair*rhoair*cgh(2)*(1.-wtg0(2)/fact)
         ! cgperl = rhoair*cgw(2)*(1.-wtgq0(2)/facq)*dqgperdT
         ! cgimpl = rhoair*cgw(2)*(1.-wtgq0(2)/facq)*dqgimpdT
         cgperl = rhoair*cgw_per*(dqgperdT &
                  - (dqgperdT*cgw_per*fgper*fg) &
                  /(caw(2) + cgw_per*fgper*fg + cgw_imp*fgimp*fg + cfw(3)*fc(3)) &
                  /facq)
         cgimpl = rhoair*cgw_imp*(dqgimpdT &
                  - (dqgimpdT*cgw_imp*fgimp*fg) &
                  /(caw(2) + cgw_per*fgper*fg + cgw_imp*fgimp*fg + cfw(3)*fc(3)) &
                  /facq)
         cgimpl = cgimpl*fwet_gimp
      ELSE !botlay == 1
         cgrnds = cpair*rhoair*cgh(1)*(1.-wta0(1)*wtg0(2)*wtg0(1)/fact-wtg0(1))
         cgperl = rhoair*cgw_per*(1.-wtaq0(1)*wtgq0(2)*wtgq0(1)/facq-wtgq0(1))*dqgperdT
         cgimpl = rhoair*cgw_imp*(1.-wtaq0(1)*wtgq0(2)*wtgq0(1)/facq-wtgq0(1))*dqgimpdT
         cgimpl = cgimpl*fwet_gimp
      ENDIF

      cgimp = cgrnds + cgimpl*htvp_gimp
      cgper = cgrnds + cgperl*htvp_gper

!-----------------------------------------------------------------------
! 2 m height air temperature above apparent sink height
!-----------------------------------------------------------------------

      !tref = thm + vonkar/(fh)*dth * (fh2m/vonkar - fh/vonkar)
      !qref =  qm + vonkar/(fq)*dqh * (fq2m/vonkar - fq/vonkar)

  END SUBROUTINE UrbanVegFlux
!----------------------------------------------------------------------


   SUBROUTINE dewfraction (sigf,lai,sai,dewmx,ldew,fwet)

!=======================================================================
! Original author: Yongjiu Dai, September 15, 1999
!
! determine fraction of foliage covered by water and
! fraction of foliage that is dry and transpiring
!
!=======================================================================

   USE MOD_Precision
   IMPLICIT NONE

   real(r8), intent(in) :: sigf   ! fraction of veg cover, excluding snow-covered veg [-]
   real(r8), intent(in) :: lai    ! leaf area index  [-]
   real(r8), intent(in) :: sai    ! stem area index  [-]
   real(r8), intent(in) :: dewmx  ! maximum allowed dew [0.1 mm]
   real(r8), intent(in) :: ldew   ! depth of water on foliage [kg/m2/s]

   real(r8), intent(out) :: fwet  ! fraction of foliage covered by water [-]

   real(r8) lsai                  ! lai + sai
   real(r8) dewmxi                ! inverse of maximum allowed dew [1/mm]
   real(r8) vegt                  ! sigf*lsai
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
      IF (ldew > 0.) THEN
         fwet = ((dewmxi/vegt)*ldew)**.666666666666

! Check for maximum limit of fwet
         fwet = min(fwet,1.0)

      ENDIF

   END SUBROUTINE dewfraction

END MODULE MOD_Urban_Flux
