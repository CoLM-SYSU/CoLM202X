#include <define.h>

 MODULE GLACIER

!-----------------------------------------------------------------------
! Energy and Mass Balance Model of LAND ICE (GLACIER / ICE SHEET)
! 
! Original author: Yongjiu Dai, /05/2014/
!-----------------------------------------------------------------------
 use precision
 IMPLICIT NONE
 SAVE

! PUBLIC MEMBER FUNCTIONS:
  public :: GLACIER_TEMP
  public :: GLACIER_WATER


! PRIVATE MEMBER FUNCTIONS:
  private :: groundfluxes_glacier
  private :: groundtem_glacier


!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------


 subroutine GLACIER_TEMP (lb    ,nl_ice     ,deltim      ,&
                    zlnd        ,zsno       ,capr        ,cnfac       ,&
                    forc_hgt_u ,forc_hgt_t  ,forc_hgt_q  ,&
                    forc_us     ,forc_vs    ,forc_t      ,forc_q      ,&
                    forc_rhoair ,forc_psrf  ,coszen      ,sabg        ,&
                    forc_frl    ,fsno       ,dz_icesno   ,z_icesno    ,&
                    zi_icesno   ,t_icesno   ,wice_icesno ,wliq_icesno ,&
                    scv         ,snowdp     ,imelt       ,taux        ,&
                    tauy        ,fsena      ,fevpa       ,lfevpa      ,&
                    fseng       ,fevpg      ,olrg        ,fgrnd       ,&
                    qseva       ,qsdew      ,qsubl       ,qfros       ,&
                    sm          ,tref       ,qref        ,trad        ,&
                    errore      ,emis       ,z0m         ,zol         ,&
                    rib         ,ustar      ,qstar       ,tstar       ,&
                    fm          ,fh         ,fq)

!=======================================================================
! this is the main subroutine to execute the calculation 
! of thermal processes and surface fluxes of the land ice (glacier and ice sheet)
!
! Original author : Yongjiu Dai and Nan Wei, /05/2014/
! 
! FLOW DIAGRAM FOR GLACIER_TEMP.F90
! 
! GLACIER_TEMP ===> qsadv
!                   groundfluxes | --------->  |moninobukini
!                                |             |moninobuk
!
!                   groundTem    | --------->  |meltf
!                               
!=======================================================================

  use precision
  use PhysicalConstants, only : hvap,hsub,rgas,cpair,stefnc,tfrz
  use FRICTION_VELOCITY

  IMPLICIT NONE
 
!---------------------Argument------------------------------------------

  integer, INTENT(in) :: &
        lb,          &! lower bound of array 
        nl_ice        ! upper bound of array

  real(r8), INTENT(in) :: &
        deltim,      &! model time step [second]
        zlnd,        &! roughness length for ice surface [m]
        zsno,        &! roughness length for snow [m]
        capr,        &! tuning factor to turn first layer T into surface T
        cnfac,       &! Crank Nicholson factor between 0 and 1

        ! Atmospherical variables and observational height
        forc_hgt_u,  &! observational height of wind [m]
        forc_hgt_t,  &! observational height of temperature [m]
        forc_hgt_q,  &! observational height of humidity [m]
        forc_us,     &! wind component in eastward direction [m/s]
        forc_vs,     &! wind component in northward direction [m/s]
        forc_t,      &! temperature at agcm reference height [kelvin]
        forc_q,      &! specific humidity at agcm reference height [kg/kg]
        forc_rhoair, &! density air [kg/m3]
        forc_psrf,   &! atmosphere pressure at the surface [pa]

        ! Radiative fluxes
        coszen,      &! cosine of the solar zenith angle
        sabg,        &! solar radiation absorbed by ground [W/m2]
        forc_frl,    &! atmospheric infrared (longwave) radiation [W/m2]

        ! State variable (1)
        fsno,        &! fraction of ground covered by snow
        dz_icesno(lb:nl_ice),  &! layer thickiness [m]
        z_icesno (lb:nl_ice),  &! node depth [m]
        zi_icesno(lb-1:nl_ice)  ! interface depth [m]

        ! State variables (2)
  real(r8), INTENT(inout) :: &
        t_icesno(lb:nl_ice),   &! snow/ice temperature [K]
        wice_icesno(lb:nl_ice),&! ice lens [kg/m2]
        wliq_icesno(lb:nl_ice),&! liqui water [kg/m2]
        scv,                   &! snow cover, water equivalent [mm, kg/m2]
        snowdp                  ! snow depth [m]

  integer, INTENT(out) :: & 
       imelt(lb:nl_ice)  ! flag for melting or freezing [-]
 
        ! Output fluxes
  real(r8), INTENT(out) :: &
        taux,        &! wind stress: E-W [kg/m/s**2]
        tauy,        &! wind stress: N-S [kg/m/s**2]
        fsena,       &! sensible heat to atmosphere [W/m2]
        lfevpa,      &! latent heat flux to atmosphere [W/m2]
        fseng,       &! sensible heat flux from ground [W/m2]
        fevpg,       &! evaporation heat flux from ground [mm/s]
        olrg,        &! outgoing long-wave radiation to atmosphere
        fgrnd,       &! ground heat flux [W/m2]

        fevpa,       &! evapotranspiration to atmosphere (mm h2o/s)
        qseva,       &! ground surface evaporation rate (mm h2o/s)
        qsdew,       &! ground surface dew formation (mm h2o /s) [+]
        qsubl,       &! sublimation rate from snow pack (mm h2o /s) [+]
        qfros,       &! surface dew added to snow pack (mm h2o /s) [+]

        sm,          &! rate of snowmelt [kg/(m2 s)]
        tref,        &! 2 m height air temperature [kelvin]
        qref,        &! 2 m height air specific humidity
        trad,        &! radiative temperature [K]

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
       cgrnd,        &! deriv. of ice energy flux wrt to ice temp [w/m2/k]
       cgrndl,       &! deriv, of ice sensible heat flux wrt ice temp [w/m2/k]
       cgrnds,       &! deriv of ice latent heat flux wrt ice temp [w/m**2/k]
       degdT,        &! d(eg)/dT
       dqgdT,        &! d(qg)/dT
       eg,           &! water vapor pressure at temperature T [pa]
       egsmax,       &! max. evaporation which ice can provide at one time step
       egidif,       &! the excess of evaporation over "egsmax"
       emg,          &! ground emissivity (0.96) 
       errore,       &! energy balnce error [w/m2]
       fact(lb:nl_ice), &! used in computing tridiagonal matrix
       htvp,         &! latent heat of vapor of water (or sublimation) [j/kg]
       qg,           &! ground specific humidity [kg/kg]
       qsatg,        &! saturated humidity [kg/kg]
       qsatgdT,      &! d(qsatg)/dT
       qred,         &! ice surface relative humidity
       thm,          &! intermediate variable (forc_t+0.0098*forc_hgt_t)
       th,           &! potential temperature (kelvin)
       thv,          &! virtual potential temperature (kelvin)
       t_grnd,       &! ground surface temperature [K]
       t_icesno_bef(lb:nl_ice), &! ice/snow temperature before update
       tinc,         &! temperature difference of two time step
       ur,           &! wind speed at reference height [m/s]
       xmf            ! total latent heat of phase change of ground water

!=======================================================================
! [1] Initial set and propositional variables
!=======================================================================

      ! temperature and water mass from previous time step
      t_grnd = t_icesno(lb)
      t_icesno_bef(lb:) = t_icesno(lb:)

      ! emissivity
      emg = 0.97

      ! latent heat, assumed that the sublimation occured only as wliq_icesno=0
      htvp = hvap
      if(wliq_icesno(lb)<=0. .and. wice_icesno(lb)>0.) htvp = hsub

      ! potential temperatur at the reference height
      thm = forc_t + 0.0098*forc_hgt_t  ! intermediate variable equivalent to
                                        ! forc_t*(pgcm/forc_psrf)**(rgas/cpair)
      th = forc_t*(100000./forc_psrf)**(rgas/cpair) ! potential T
      thv = th*(1.+0.61*forc_q)         ! virtual potential T
      ur = max(0.1,sqrt(forc_us*forc_us+forc_vs*forc_vs))   ! limit set to 0.1

!=======================================================================
! [2] specific humidity and its derivative at ground surface
!=======================================================================

      qred = 1.
      call qsadv(t_grnd,forc_psrf,eg,degdT,qsatg,qsatgdT)

      qg = qred*qsatg  
      dqgdT = qred*qsatgdT

!=======================================================================
! [3] Compute sensible and latent fluxes and their derivatives with respect 
!     to ground temperature using ground temperatures from previous time step.
!=======================================================================

      call groundfluxes_glacier (zlnd,zsno,forc_hgt_u,forc_hgt_t,forc_hgt_q,&
                        forc_us,forc_vs,forc_t,forc_q,forc_rhoair,forc_psrf, &
                        ur,thm,th,thv,t_grnd,qg,dqgdT,htvp,&
                        fsno,cgrnd,cgrndl,cgrnds,&
                        taux,tauy,fsena,fevpa,fseng,fevpg,tref,qref,&
                        z0m,zol,rib,ustar,qstar,tstar,fm,fh,fq)

!=======================================================================
! [4] Gound temperature
!=======================================================================

      call groundtem_glacier (lb,nl_ice,deltim,&
                     capr,cnfac,dz_icesno,z_icesno,zi_icesno,&
                     t_icesno,wice_icesno,wliq_icesno,scv,snowdp,&
                     forc_frl,sabg,fseng,fevpg,cgrnd,htvp,emg,&
                     imelt,sm,xmf,fact)

!=======================================================================
! [5] Correct fluxes to present ice temperature
!=======================================================================

      t_grnd = t_icesno(lb)
      tinc = t_icesno(lb) - t_icesno_bef(lb)
      fseng = fseng + tinc*cgrnds 
      fevpg = fevpg + tinc*cgrndl

! calculation of evaporative potential; flux in kg m-2 s-1.  
! egidif holds the excess energy if all water is evaporated
! during the timestep. this energy is later added to the sensible heat flux.

      egsmax = (wice_icesno(lb)+wliq_icesno(lb)) / deltim

      egidif = max( 0., fevpg - egsmax )
      fevpg = min ( fevpg, egsmax )
      fseng = fseng + htvp*egidif

! total fluxes to atmosphere
      fsena = fseng
      fevpa = fevpg
      lfevpa= htvp*fevpg   ! W/m^2 (accouting for sublimation)
      
      qseva = 0.
      qsubl = 0.
      qfros = 0.
      qsdew = 0.

      if(fevpg >= 0)then 
         qseva = min(wliq_icesno(lb)/deltim, fevpg)
         qsubl = fevpg - qseva
      else 
        if(t_grnd < tfrz)then
           qfros = abs(fevpg)
        else
           qsdew = abs(fevpg)
        endif
      end if

! ground heat flux
      fgrnd = sabg + emg*forc_frl &
            - emg*stefnc*t_icesno_bef(lb)**3*(t_icesno_bef(lb) + 4.*tinc) &
            - (fseng+fevpg*htvp)

! outgoing long-wave radiation from ground
      olrg = (1.-emg)*forc_frl + emg*stefnc * t_icesno_bef(lb)**4 &
! for conservation we put the increase of ground longwave to outgoing
           + 4.*emg*stefnc*t_icesno_bef(lb)**3*tinc

! averaged bulk surface emissivity 
      emis = emg

! radiative temperature
      trad = (olrg/stefnc)**0.25

!=======================================================================
! [6] energy balance error
!=======================================================================

      errore = sabg + forc_frl - olrg - fsena - lfevpa - xmf
      do j = lb, nl_ice
         errore = errore - (t_icesno(j)-t_icesno_bef(j))/fact(j)
      enddo

!#if (defined CLMDEBUG)
     if(abs(errore)>.2)then 
     write(6,*) 'GLACIER_TEMP.F90 : energy  balance violation'
     write(6,100) errore,sabg,forc_frl,olrg,fsena,lfevpa,xmf
     endif
100  format(10(f7.3))
!#endif

 end subroutine GLACIER_TEMP



 subroutine groundfluxes_glacier (zlnd,zsno,hu,ht,hq,&
                                  us,vs,tm,qm,rhoair,psrf,&
                                  ur,thm,th,thv,t_grnd,qg,dqgdT,htvp,&
                                  fsno,cgrnd,cgrndl,cgrnds,& 
                                  taux,tauy,fsena,fevpa,fseng,fevpg,tref,qref,&
                                  z0m,zol,rib,ustar,qstar,tstar,fm,fh,fq)

!=======================================================================
! this is the main subroutine to execute the calculation of thermal processes
! and surface fluxes of land ice (glacier and ice sheet)
!
! Original author : Yongjiu Dai and Nan Wei, /05/2014/
!=======================================================================

  use precision
  use PhysicalConstants, only : cpair,vonkar,grav
  use FRICTION_VELOCITY
  implicit none
 
!----------------------- Dummy argument --------------------------------
  real(r8), INTENT(in) :: &
        zlnd,     &! roughness length for ice [m]
        zsno,     &! roughness length for snow [m]

        ! atmospherical variables and observational height
        hu,       &! observational height of wind [m]
        ht,       &! observational height of temperature [m]
        hq,       &! observational height of humidity [m]
        us,       &! wind component in eastward direction [m/s]
        vs,       &! wind component in northward direction [m/s]
        tm,       &! temperature at agcm reference height [kelvin] [not used]
        qm,       &! specific humidity at agcm reference height [kg/kg]
        rhoair,   &! density air [kg/m3]
        psrf,     &! atmosphere pressure at the surface [pa] [not used]

        fsno,     &! fraction of ground covered by snow

        ur,       &! wind speed at reference height [m/s]
        thm,      &! intermediate variable (tm+0.0098*ht)
        th,       &! potential temperature (kelvin)
        thv,      &! virtual potential temperature (kelvin)

        t_grnd,   &! ground surface temperature [K]
        qg,       &! ground specific humidity [kg/kg]
        dqgdT,    &! d(qg)/dT
        htvp       ! latent heat of vapor of water (or sublimation) [j/kg]

  real(r8), INTENT(out) :: &
        taux,     &! wind stress: E-W [kg/m/s**2]
        tauy,     &! wind stress: N-S [kg/m/s**2]
        fsena,    &! sensible heat to atmosphere [W/m2]
        fevpa,    &! evapotranspiration to atmosphere [mm/s]
        fseng,    &! sensible heat flux from ground [W/m2]
        fevpg,    &! evaporation heat flux from ground [mm/s]
        cgrnd,    &! deriv. of ice energy flux wrt to ice temp [w/m2/k]
        cgrndl,   &! deriv, of ice sensible heat flux wrt ice temp [w/m2/k]
        cgrnds,   &! deriv of ice latent heat flux wrt ice temp [w/m**2/k]
        tref,     &! 2 m height air temperature [kelvin]
        qref,     &! 2 m height air humidity

        z0m,      &! effective roughness [m]
        zol,      &! dimensionless height (z/L) used in Monin-Obukhov theory
        rib,      &! bulk Richardson number in surface layer
        ustar,    &! friction velocity [m/s]
        tstar,    &! temperature scaling parameter
        qstar,    &! moisture scaling parameter
        fm,       &! integral of profile function for momentum
        fh,       &! integral of profile function for heat
        fq         ! integral of profile function for moisture

!------------------------ LOCAL VARIABLES ------------------------------
  integer niters, &! maximum number of iterations for surface temperature
       iter,      &! iteration index
       nmozsgn     ! number of times moz changes sign

  real(r8) :: &
       beta,      &! coefficient of conective velocity [-]
       displax,   &! zero-displacement height [m]
       dth,       &! diff of virtual temp. between ref. height and surface
       dqh,       &! diff of humidity between ref. height and surface
       dthv,      &! diff of vir. poten. temp. between ref. height and surface
       obu,       &! monin-obukhov length (m)
       obuold,    &! monin-obukhov length from previous iteration
       ram,       &! aerodynamical resistance [s/m]
       rah,       &! thermal resistance [s/m]
       raw,       &! moisture resistance [s/m]
       raih,      &! temporary variable [kg/m2/s]
       raiw,      &! temporary variable [kg/m2/s]
       fh2m,      &! relation for temperature at 2m
       fq2m,      &! relation for specific humidity at 2m
       fm10m,     &! integral of profile function for momentum at 10m
       thvstar,   &! virtual potential temperature scaling parameter
       um,        &! wind speed including the stablity effect [m/s]
       wc,        &! convective velocity [m/s]
       wc2,       &! wc**2
       zeta,      &! dimensionless height used in Monin-Obukhov theory
       zii,       &! convective boundary height [m]
       zldis,     &! reference height "minus" zero displacement heght [m]
       z0mg,      &! roughness length over ground, momentum [m]
       z0hg,      &! roughness length over ground, sensible heat [m]
       z0qg        ! roughness length over ground, latent heat [m]

!----------------------- Dummy argument --------------------------------
! initial roughness length
      if(fsno > 0.)then
       ! z0mg = zsno
         z0mg = 0.002 ! Table 1 of Brock et al., (2006)
         z0hg = z0mg
         z0qg = z0mg
      else
       ! z0mg = zlnd
         z0mg = 0.001 ! Table 1 of Brock et al., (2006)
         z0hg = z0mg
         z0qg = z0mg
      endif

! potential temperatur at the reference height
      beta = 1.      ! -  (in computing W_*)
      zii = 1000.    ! m  (pbl height)
      z0m = z0mg

!-----------------------------------------------------------------------
!     Compute sensible and latent fluxes and their derivatives with respect 
!     to ground temperature using ground temperatures from previous time step.
!-----------------------------------------------------------------------
! Initialization variables
      nmozsgn = 0
      obuold = 0.

      dth   = thm-t_grnd
      dqh   = qm-qg
      dthv  = dth*(1.+0.61*qm)+0.61*th*dqh
      zldis = hu-0.

      call moninobukini(ur,th,thm,thv,dth,dqh,dthv,zldis,z0mg,um,obu)
 
! Evaluated stability-dependent variables using moz from prior iteration
      niters=6

      !----------------------------------------------------------------
      ITERATION : do iter = 1, niters         ! begin stability iteration
      !----------------------------------------------------------------
         displax = 0.
         call moninobuk(hu,ht,hq,displax,z0mg,z0hg,z0qg,obu,um,&
                        ustar,fh2m,fq2m,fm10m,fm,fh,fq)

         tstar = vonkar/fh*dth
         qstar = vonkar/fq*dqh

         z0hg = z0mg/exp(0.13 * (ustar*z0mg/1.5e-5)**0.45)
         z0qg = z0hg

         thvstar=tstar+0.61*th*qstar
         zeta=zldis*vonkar*grav*thvstar/(ustar**2*thv)
         if(zeta >= 0.) then     !stable
           zeta = min(2.,max(zeta,1.e-6))
         else                    !unstable
           zeta = max(-100.,min(zeta,-1.e-6))
         endif
         obu = zldis/zeta

         if(zeta >= 0.)then
           um = max(ur,0.1)
         else
           wc = (-grav*ustar*thvstar*zii/thv)**(1./3.)
          wc2 = beta*beta*(wc*wc)
           um = sqrt(ur*ur+wc2)
         endif

         if (obuold*obu < 0.) nmozsgn = nmozsgn+1
         if(nmozsgn >= 4) EXIT

         obuold = obu

      !----------------------------------------------------------------
      enddo ITERATION                         ! end stability iteration
      !----------------------------------------------------------------

! Get derivative of fluxes with repect to ground temperature
      ram    = 1./(ustar*ustar/um)
      rah    = 1./(vonkar/fh*ustar) 
      raw    = 1./(vonkar/fq*ustar) 

      raih   = rhoair*cpair/rah
      raiw   = rhoair/raw          
      cgrnds = raih
      cgrndl = raiw*dqgdT
      cgrnd  = cgrnds + htvp*cgrndl

      zol = zeta
      rib = min(5.,zol*ustar**2/(vonkar**2/fh*um**2))

! surface fluxes of momentum, sensible and latent 
! using ground temperatures from previous time step
      taux   = -rhoair*us/ram        
      tauy   = -rhoair*vs/ram
      fseng  = -raih*dth
      fevpg  = -raiw*dqh 

      fsena  = fseng
      fevpa  = fevpg

! 2 m height air temperature
      tref   = (thm + vonkar/fh*dth * (fh2m/vonkar - fh/vonkar))
      qref   = ( qm + vonkar/fq*dqh * (fq2m/vonkar - fq/vonkar))

 end subroutine groundfluxes_glacier



 subroutine groundtem_glacier (lb,nl_ice,deltim,&
                      capr,cnfac,dz_icesno,z_icesno,zi_icesno,&
                      t_icesno,wice_icesno,wliq_icesno,scv,snowdp,&
                      forc_frl,sabg,fseng,fevpg,cgrnd,htvp,emg,&
                      imelt,sm,xmf,fact)

!=======================================================================
! SNOW and LAND ICE temperatures
! o The volumetric heat capacity is calculated as a linear combination
!   in terms of the volumetric fraction of the constituent phases.
! o The thermal conductivity of snow/ice is computed from
!   the formulation used in SNTHERM (Jordan 1991) and Yen (1981), respectively.
! o Boundary conditions:
!   F = Rnet - Hg - LEg (top), F= 0 (base of the land ice column).
! o Ice/snow temperature is predicted from heat conduction
!   in 10 ice layers and up to 5 snow layers.
!   The thermal conductivities at the interfaces between two neighbor layers
!   (j, j+1) are derived from an assumption that the flux across the interface
!   is equal to that from the node j to the interface and the flux from the
!   interface to the node j+1. The equation is solved using the Crank-Nicholson
!   method and resulted in a tridiagonal system equation.
!
! Phase change (see meltf.F90)
! 
! Original author : Yongjiu Dai, /05/2014/
!=======================================================================

  use precision
  use PhysicalConstants, only : stefnc,cpice,cpliq,denh2o,denice,tfrz,tkwat,tkice,tkair

  IMPLICIT NONE

  integer, INTENT(in) :: lb         !lower bound of array
  integer, INTENT(in) :: nl_ice     !upper bound of array
  real(r8), INTENT(in) :: deltim    !seconds in a time step [second]
  real(r8), INTENT(in) :: capr      !tuning factor to turn first layer T into surface T
  real(r8), INTENT(in) :: cnfac     !Crank Nicholson factor between 0 and 1

  real(r8), INTENT(in) :: dz_icesno(lb:nl_ice)   !layer thickiness [m]
  real(r8), INTENT(in) :: z_icesno (lb:nl_ice)   !node depth [m]
  real(r8), INTENT(in) :: zi_icesno(lb-1:nl_ice) !interface depth [m]

  real(r8), INTENT(in) :: sabg      !solar radiation absorbed by ground [W/m2]
  real(r8), INTENT(in) :: forc_frl  !atmospheric infrared (longwave) radiation [W/m2]
  real(r8), INTENT(in) :: fseng     !sensible heat flux from ground [W/m2]
  real(r8), INTENT(in) :: fevpg     !evaporation heat flux from ground [mm/s]
  real(r8), INTENT(in) :: cgrnd     !deriv. of ice energy flux wrt to ice temp [W/m2/k]
  real(r8), INTENT(in) :: htvp      !latent heat of vapor of water (or sublimation) [J/kg]
  real(r8), INTENT(in) :: emg       !ground emissivity (0.97 for snow,

  real(r8), INTENT(inout) :: t_icesno (lb:nl_ice) !snow and ice temperature [K]
  real(r8), INTENT(inout) :: wice_icesno(lb:nl_ice) !ice lens [kg/m2]
  real(r8), INTENT(inout) :: wliq_icesno(lb:nl_ice) !liqui water [kg/m2]
  real(r8), INTENT(inout) :: scv    !snow cover, water equivalent [mm, kg/m2]
  real(r8), INTENT(inout) :: snowdp !snow depth [m]

  real(r8), INTENT(out) :: sm       !rate of snowmelt [kg/(m2 s)]
  real(r8), INTENT(out) :: xmf      !total latent heat of phase change of ground water
  real(r8), INTENT(out) :: fact(lb:nl_ice) !used in computing tridiagonal matrix
  integer, INTENT(out)  :: imelt(lb:nl_ice)    !flag for melting or freezing [-]

!------------------------ local variables ------------------------------
  real(r8) rhosnow         ! partitial density of water (ice + liquid)
  real(r8) cv(lb:nl_ice)   ! heat capacity [J/(m2 K)]
  real(r8) thk(lb:nl_ice)  ! thermal conductivity of layer
  real(r8) tk(lb:nl_ice)   ! thermal conductivity [W/(m K)]

  real(r8) at(lb:nl_ice)   !"a" vector for tridiagonal matrix
  real(r8) bt(lb:nl_ice)   !"b" vector for tridiagonal matrix
  real(r8) ct(lb:nl_ice)   !"c" vector for tridiagonal matrix
  real(r8) rt(lb:nl_ice)   !"r" vector for tridiagonal solution

  real(r8) fn  (lb:nl_ice) ! heat diffusion through the layer interface [W/m2]
  real(r8) fn1 (lb:nl_ice) ! heat diffusion through the layer interface [W/m2]
  real(r8) dzm             ! used in computing tridiagonal matrix
  real(r8) dzp             ! used in computing tridiagonal matrix

  real(r8) t_icesno_bef(lb:nl_ice) ! snow/ice temperature before update
  real(r8) hs              ! net energy flux into the surface (w/m2)
  real(r8) dhsdt           ! d(hs)/dT
  real(r8) brr(lb:nl_ice)  ! temporay set

  integer i,j

!=======================================================================
! SNOW and LAND ICE heat capacity 
      cv(1:) = wice_icesno(1:)*cpice + wliq_icesno(1:)*cpliq
      if(lb==1 .and. scv>0.) cv(1) = cv(1) + cpice*scv

      if(lb<=0)then
        cv(:0) = cpliq*wliq_icesno(:0) + cpice*wice_icesno(:0)
      endif

! SNOW and LAND ICE thermal conductivity [W/(m K)]
      do j = lb, nl_ice
         thk(j) = tkwat
         if(t_icesno(j)<=tfrz) thk(j) = 9.828*exp(-0.0057*t_icesno(j))
      enddo

      if(lb < 1)then
        do j = lb, 0
          rhosnow = (wice_icesno(j)+wliq_icesno(j))/dz_icesno(j)

        ! presently option [1] is the default option
        ! [1] Jordan (1991) pp. 18
          thk(j) = tkair+(7.75e-5*rhosnow+1.105e-6*rhosnow*rhosnow)*(tkice-tkair)

        ! [2] Sturm et al (1997)
        ! thk(j) = 0.0138 + 1.01e-3*rhosnow + 3.233e-6*rhosnow**2
        ! [3] Ostin and Andersson presented in Sturm et al., (1997)
        ! thk(j) = -0.871e-2 + 0.439e-3*rhosnow + 1.05e-6*rhosnow**2
        ! [4] Jansson(1901) presented in Sturm et al. (1997)
        ! thk(j) = 0.0293 + 0.7953e-3*rhosnow + 1.512e-12*rhosnow**2
        ! [5] Douville et al., (1995)
        ! thk(j) = 2.2*(rhosnow/denice)**1.88
        ! [6] van Dusen (1992) presented in Sturm et al. (1997)
        ! thk(j) = 0.021 + 0.42e-3*rhosnow + 0.22e-6*rhosnow**2
        enddo
      endif

! Thermal conductivity at the layer interface
      do j = lb, nl_ice-1

! the following consideration is try to avoid the snow conductivity
! to be dominant in the thermal conductivity of the interface.
! Because when the distance of bottom snow node to the interfacee
! is larger than that of interface to top ice node,
! the snow thermal conductivity will be dominant, and the result is that
! lees heat tranfer between snow and ice 
         if((j==0) .and. (z_icesno(j+1)-zi_icesno(j)<zi_icesno(j)-z_icesno(j)))then
            tk(j) = 2.*thk(j)*thk(j+1)/(thk(j)+thk(j+1))
            tk(j) = max(0.5*thk(j+1),tk(j))
         else
            tk(j) = thk(j)*thk(j+1)*(z_icesno(j+1)-z_icesno(j)) &
                  /(thk(j)*(z_icesno(j+1)-zi_icesno(j))+thk(j+1)*(zi_icesno(j)-z_icesno(j)))
         endif
      enddo
      tk(nl_ice) = 0.



! net ground heat flux into the surface and its temperature derivative
      hs = sabg + emg*forc_frl - emg*stefnc*t_icesno(lb)**4 - (fseng+fevpg*htvp) 

      dhsdT = - cgrnd - 4.*emg * stefnc * t_icesno(lb)**3
      t_icesno_bef(lb:) = t_icesno(lb:)

      j       = lb
      fact(j) = deltim / cv(j) * dz_icesno(j) &
              / (0.5*(z_icesno(j)-zi_icesno(j-1)+capr*(z_icesno(j+1)-zi_icesno(j-1))))

      do j = lb + 1, nl_ice
         fact(j) = deltim/cv(j)
      enddo

      do j = lb, nl_ice - 1
        fn(j) = tk(j)*(t_icesno(j+1)-t_icesno(j))/(z_icesno(j+1)-z_icesno(j))
      enddo
      fn(nl_ice) = 0.

! set up vector r and vectors a, b, c that define tridiagonal matrix
      j     = lb
      dzp   = z_icesno(j+1)-z_icesno(j)
      at(j) = 0.
      bt(j) = 1+(1.-cnfac)*fact(j)*tk(j)/dzp-fact(j)*dhsdT
      ct(j) =  -(1.-cnfac)*fact(j)*tk(j)/dzp
      rt(j) = t_icesno(j) + fact(j)*( hs - dhsdT*t_icesno(j) + cnfac*fn(j) )


      do j = lb + 1, nl_ice - 1
         dzm   = (z_icesno(j)-z_icesno(j-1))
         dzp   = (z_icesno(j+1)-z_icesno(j))
         at(j) =   - (1.-cnfac)*fact(j)* tk(j-1)/dzm
         bt(j) = 1.+ (1.-cnfac)*fact(j)*(tk(j)/dzp + tk(j-1)/dzm)
         ct(j) =   - (1.-cnfac)*fact(j)* tk(j)/dzp
         rt(j) = t_icesno(j) + cnfac*fact(j)*( fn(j) - fn(j-1) )
      end do

      j     =  nl_ice
      dzm   = (z_icesno(j)-z_icesno(j-1))
      at(j) =   - (1.-cnfac)*fact(j)*tk(j-1)/dzm
      bt(j) = 1.+ (1.-cnfac)*fact(j)*tk(j-1)/dzm
      ct(j) = 0.
      rt(j) = t_icesno(j) - cnfac*fact(j)*fn(j-1)

! solve for t_icesno
      i = size(at)
      call tridia (i ,at ,bt ,ct ,rt ,t_icesno) 

!=======================================================================
! melting or freezing 
!=======================================================================

      do j = lb, nl_ice - 1
         fn1(j) = tk(j)*(t_icesno(j+1)-t_icesno(j))/(z_icesno(j+1)-z_icesno(j))
      enddo
      fn1(nl_ice) = 0.

      j = lb
      brr(j) = cnfac*fn(j) + (1.-cnfac)*fn1(j)

      do j = lb + 1, nl_ice
         brr(j) = cnfac*(fn(j)-fn(j-1)) + (1.-cnfac)*(fn1(j)-fn1(j-1))
      enddo

      call meltf (lb,nl_ice,deltim, &
                  fact(lb:),brr(lb:),hs,dhsdT, &
                  t_icesno_bef(lb:),t_icesno(lb:),wliq_icesno(lb:),wice_icesno(lb:),imelt(lb:), &
                  scv,snowdp,sm,xmf)

!-----------------------------------------------------------------------

 end subroutine groundtem_glacier



 subroutine GLACIER_WATER ( nl_ice,maxsnl,deltim,&
                    z_icesno    ,dz_icesno   ,zi_icesno ,t_icesno,&
                    wliq_icesno ,wice_icesno ,pg_rain   ,pg_snow ,&
                    sm          ,scv         ,snowdp    ,imelt   ,&
                    fiold       ,snl         ,qseva     ,qsdew   ,&
                    qsubl       ,qfros       ,rsur      ,rnof    ,&
                    ssi         ,wimp    )
	           
!=======================================================================
  use precision
  use PhysicalConstants, only : denice, denh2o, tfrz
  use SNOW_Layers_CombineDivide
  use SOIL_SNOW_hydrology

  IMPLICIT NONE
    
!-----------------------Argument---------- ------------------------------
  integer, INTENT(in) :: nl_ice  ! upper bound of array
  integer, INTENT(in) :: maxsnl  ! maximum number of snow layers

  real(r8), INTENT(in) :: &
       deltim    , &! time step (s)
       ssi       , &! irreducible water saturation of snow
       wimp      , &! water impremeable if porosity less than wimp
       pg_rain   , &! rainfall (mm h2o/s)
       pg_snow   , &! snowfall (mm h2o/s)
       sm        , &! snow melt (mm h2o/s)
       qseva     , &! ground surface evaporation rate (mm h2o/s)
       qsdew     , &! ground surface dew formation (mm h2o /s) [+]
       qsubl     , &! sublimation rate from snow pack (mm h2o /s) [+]
       qfros     , &! surface dew added to snow pack (mm h2o /s) [+]
       fiold(maxsnl+1:nl_ice)  ! fraction of ice relative to the total water

  integer, INTENT(in) :: imelt(maxsnl+1:nl_ice)  ! flag for: melting=1, freezing=2, nothing happended=0
  integer, INTENT(inout) :: snl ! lower bound of array

  real(r8), INTENT(inout) :: &
       z_icesno   (maxsnl+1:nl_ice) , &! layer depth (m)
       dz_icesno  (maxsnl+1:nl_ice) , &! layer thickness (m)
       zi_icesno  (maxsnl  :nl_ice) , &! interface level below a "z" level (m)
       t_icesno   (maxsnl+1:nl_ice) , &! snow/ice skin temperature (K)
       wice_icesno(maxsnl+1:nl_ice) , &! ice lens (kg/m2)
       wliq_icesno(maxsnl+1:nl_ice) , &! liquid water (kg/m2)
       scv       , &! snow mass (kg/m2)
       snowdp       ! snow depth (m)

  real(r8), INTENT(out) :: &
       rsur      , &! surface runoff (mm h2o/s)
       rnof         ! total runoff (mm h2o/s)
!                    
!-----------------------Local Variables------------------------------
!                   
  integer lb, j

  real(r8) :: gwat   ! net water input from top (mm/s)
  real(r8) :: rsubst ! subsurface runoff (mm h2o/s)

!=======================================================================
! [1] update the liquid water within snow layer and the water onto the ice surface
!
! Snow melting is treated in a realistic fashion, with meltwater
! percolating downward through snow layers as long as the snow is unsaturated.
! Once the underlying snow is saturated, any additional meltwater runs off.
! When glacier ice melts, however, the meltwater is assumed to remain in place until it refreezes.
! In warm parts of the ice sheet, the meltwater does not refreeze, but stays in place indefinitely.
!=======================================================================

      lb = snl + 1 
      if (lb>=1)then
         gwat = pg_rain + sm - qseva
      else
         call snowwater (lb,deltim,ssi,wimp,&
                         pg_rain,qseva,qsdew,qsubl,qfros,&
                         dz_icesno(lb:0),wice_icesno(lb:0),wliq_icesno(lb:0),gwat)
      endif

!=======================================================================
! [2] surface runoff and infiltration
!=======================================================================

      rsur = max(0.0,gwat)
      rsubst = 0.
      rnof = rsur

      if(snl<0)then
         ! Compaction rate for snow
         ! Natural compaction and metamorphosis. The compaction rate
         ! is recalculated for every new timestep
         lb  = snl + 1   ! lower bound of array
         call snowcompaction (lb,deltim,&
                         imelt(lb:0),fiold(lb:0),t_icesno(lb:0),&
                         wliq_icesno(lb:0),wice_icesno(lb:0),dz_icesno(lb:0))

         ! Combine thin snow elements
         lb = maxsnl + 1
         call snowlayerscombine (lb,snl,&
                         z_icesno(lb:1),dz_icesno(lb:1),zi_icesno(lb-1:1),&
                         wliq_icesno(lb:1),wice_icesno(lb:1),t_icesno(lb:1),scv,snowdp)

         ! Divide thick snow elements
         if(snl<0) &
         call snowlayersdivide (lb,snl,&
                         z_icesno(lb:0),dz_icesno(lb:0),zi_icesno(lb-1:0),&
                         wliq_icesno(lb:0),wice_icesno(lb:0),t_icesno(lb:0))
      endif

      if (snl > maxsnl) then
         wice_icesno(maxsnl+1:snl) = 0.
         wliq_icesno(maxsnl+1:snl) = 0.
         t_icesno   (maxsnl+1:snl) = 0.
         z_icesno   (maxsnl+1:snl) = 0.
         dz_icesno  (maxsnl+1:snl) = 0.
      endif


 end subroutine GLACIER_WATER



 END MODULE GLACIER
