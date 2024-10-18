#include <define.h>

MODULE  MOD_Lake_Subs

   USE MOD_Precision
   IMPLICIT NONE
   
   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: LakStaParms
   PUBLIC :: LakStaParmsSubs

!-----------------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------------

   SUBROUTINE LakStaParms ( &
               ! "in" arguments
               ! -------------------
               hwref    , htref    , hqref    , usurf    ,&
               vsurf    , tmref    , qmref    , tskin    ,&
               arhos    , psurf    , zcb      ,& 
               ! "out" arguments
               ! -------------------
               zol      , rib      , wdm      , shfdt    ,&
               taux     , tauy     , t2m      , q2m      ,&
               u10m     , v10m     , fh2m     , fq2m     ,&
               fm10m    , fq10m    , fm       , fh       ,&
               fq       , ustar    , qstar    , tstar    ,&
               ram      , rah      , raw      , emis      )

! ------------------------ code history -------------------------------------
! Oldname: lakesp
! Description:
!     Calculating lake surface stability-related parameters
!
! Original author:  
!     Xin-Zhong Liang, 2021-11-01
!
! Revisions:
!     Omarjan Obulkasim, 04/2024: Updated variable names and renamed subrountines
!----------------------------------------------------------------------------
   USE MOD_Qsadv, only : qsadv
   USE MOD_FrictionVelocity, only : moninobukini, moninobuk
   USE MOD_Lake_Const, only :tfrz, hvap, hfus, hsub, grav,&
                                 vonkar,cpair,beta_wc,emisw,rgas
!================================================================================
!  ------------------------- input variables ---------------------------
   real(r8), intent(in)     :: &
      hwref                   ,&! observational height of wind [m]
      htref                   ,&! observational height of temperature [m]
      hqref                   ,&! observational height of humidity [m]
      usurf                   ,&! wind component in eastward direction [m/s]
      vsurf                   ,&! wind component in northward direction [m/s]
      tmref                   ,&! temperature at the reference height [K]
      qmref                   ,&! specific humidity at the reference height [kg/kg]
      arhos                   ,&! air density [kg/m3]
      psurf                   ,&! surface pressure [Pa]
      zcb                     ,&! convective boundary height [m]
      tskin                     ! surface temperature [K]

!  ------------------------- inout variables ---------------------------

   real(r8), intent(out)    :: &
      taux                    ,&! x-wind stress [kg/m/s2]
      tauy                    ,&! y-wind stress [kg/m/s2]
      t2m                     ,&! 2 m height air temperature [K]
      q2m                     ,&! 2 m height air specific humidity
      shfdt                   ,&! derivative of srf sen+lat  heat flux wrt srf temp [W/m2/K]
      zol                     ,&! dimensionless height (z/L) used in Monin-Obukhov theory
      rib                     ,&! bulk Richardson number in surface layer
      ustar                   ,&! u* in similarity theory [m/s]
      qstar                   ,&! q* in similarity theory [kg/kg]
      tstar                   ,&! t* in similarity theory [K]
      ram                     ,&! aerodynamic resistance [s/m]
      rah                     ,&! thermal resistance [s/m]
      raw                     ,&! moisture resistance [s/m]
      u10m                    ,&! u-velocity at 10m [m/s]
      v10m                    ,&! v-velocity at 10m [m/s]
      wdm                     ,&! lowest level wind speed including the stablity effect [m/s]
      emis                    ,&! averaged bulk surface emissivity
      fm10m                   ,&! integral of profile function for momentum at 10m
      fq10m                   ,&! integral of profile function for moisture at 10m
      fh2m                    ,&! integral of profile function for heat at 2m
      fq2m                    ,&! integral of profile function for moisture at 2m
      fm                      ,&! integral of profile function for momentum
      fh                      ,&! integral of profile function for heat
      fq                        ! integral of profile function for moisture

!  ------------------------- local variables ---------------------------
   integer                  :: &
      iter                    ,&! iteration index
      nmozsgn                   ! number of times moz changes sign

   real(r8)                 :: &
      eg                      ,&! water vapor pressure at temperature T [pa]
      qsatg                   ,&! saturated humidity [kg/kg]
      qsatgdT                 ,&! d(qsatg)/dT
      degdT                   ,&! d(eg)/dT
      thm                     ,&! potential temperaure at htref with reference to surface pressure [K]
      th                      ,&! atmospheric potential temperaure [K]
      displax                 ,&! zero- displacement height [m]
      dqh                     ,&! diff of humidity between ref. height and surface
      dth                     ,&! diff of virtual temp. between ref. height and surface
      dthv                    ,&! diff of vir. poten. temp. between ref. height and surface
      htvp                    ,&! latent heat of vapor of water (or sublimation) [J/kg]
      obu                     ,&! monin-obukhov length [m]
      obuold                  ,&! monin-obukhov length of previous iteration
      raih                    ,&! temporary variable [kg/m2/s]
      raiw                    ,&! temporary variable [kg/m2/s]
      shfdtl                  ,&! derivative of srf sensible heat flux wrt srf temp [W/m2/K]
      shfdts                  ,&! derivative of srf latent   heat flux wrt srf temp [W/m2/K]
      temp1                   ,&! relation for potential temperature profile
      temp2                   ,&! relation for specific humidity profile
      temp12m                 ,&! relation for temperature at 2m
      temp22m                 ,&! relation for specific humidity at 2m
      thv                     ,&! virtual potential temperature [K]
      thvstar                 ,&! virtual potential temperature scaling parameter
      um                      ,&! wind speed including the stablity effect [m/s]
      ur                      ,&! wind speed at reference height [m/s]
      tmc                     ,&! temperature at the reference height [C]
      visa                    ,&! kinematic viscosity of dry air [m2/s]
      wc                      ,&! convective velocity [m/s]
      xt                      ,&!
      xq                      ,&!
      emg                     ,&! ground emissivity
      zeta                    ,&! dimensionless height used in Monin-Obukhov theory
      zldis                   ,&! reference height minus zero-displacement height [m]
      z0mg                    ,&! roughness length over ground, momentum [m]
      z0hg                    ,&! roughness length over ground, sensible heat [m]
      z0qg                      ! roughness length over ground, latent heat [m]
   ! real(r8), parameter :: zcb = 1000.                      
   integer,  parameter :: niters  = 3    ! maximum number of iteration
   integer ii                            ! DO loop or array index
!================================================================================

      !----------------------------------------------------------------------------
      ! +WMEJ fix a bug : should USE tskin instead of tmref
      IF (tskin > tfrz) THEN
         htvp = hvap
      ELSE
         htvp = hsub
      ENDIF
      emg = emisw
      CALL qsadv(tskin,psurf,eg,degdT,qsatg,qsatgdT)
      thm   = tmref + 0.0098*hwref  ! intermediate variable equivalent to
                           ! tmref*(pgcm/forc_psrf)**(rgas/cpair)
      th    = tmref*(100000./psurf)**(rgas/cpair) ! potential T
      thv   = th*(1.+0.61*qmref)         ! virtual potential T
      ur    = max(0.1,sqrt(usurf*usurf+vsurf*vsurf))   ! limit set to 0.1

      ! Initialization variables
      nmozsgn = 0
      obuold  = 0.
      dth     = thm-tskin
      dqh     = qmref-qsatg
      dthv    = dth*(1.+0.61*qmref)+0.61*th*dqh
      zldis   = hwref-0.

      ! aerodynamic roughness
      ! Kinematic viscosity of dry air (m2/s)- Andreas (1989) CRREL Rep. 89-11
      visa = 1.326e-5*(1.+6.542e-3*(tmref-tfrz) &
      + 8.301e-6*(tmref-tfrz)**2 - 4.84e-9*(tmref-tfrz)**3)

      ! loop to obtain initial and good ustar and zo
      ustar = 0.06
      wc    = 0.5
      ! zcb   = 1000.    ! m  (pbl height)
      IF (tskin.ge.tfrz) THEN        ! unfrozen lake
         IF (dthv.ge.0.) THEN
            um = ur
         ELSE
            um = sqrt(ur*ur+wc*wc)
         ENDIF
         DO ii = 1,5
            z0mg = 0.013*ustar*ustar/grav + 0.11*visa/ustar
            ustar = vonkar*um/log(zldis/z0mg)
         ENDDO
      ELSE                        ! frozen lake
         z0mg = 0.04
      ENDIF
      z0qg = z0mg
      z0hg = z0mg
      CALL moninobukini(ur,th,thm,thv,dth,dqh,dthv,zldis,z0mg,um,obu)
      !----------------------------------------------------------------------------
      DO iter = 1, niters         ! begin stability iteration
         IF (tskin.ge.tfrz) THEN     ! unfrozen lake
            z0mg = 0.013*ustar*ustar/grav + 0.11*visa/ustar
            xq = 2.67*(ustar*z0mg/visa)**0.25 - 2.57
            xt = xq
            z0qg = z0mg/exp(xq)
            z0hg = z0mg/exp(xt)
         ENDIF
         ! Evaluated stability-dependent variables using moz from prior iteration
         displax = 0.
         ! fq10m needs to be calculated in moninobuk (IF needed)
         CALL moninobuk(hwref,htref,hqref,displax,z0mg,z0hg,z0qg,obu,um,&
               ustar,fh2m,fq2m,fm10m,fm,fh,fq)
         obuold = obu
         CALL qsadv(tskin,psurf,eg,degdT,qsatg,qsatgdT)
         dth = thm-tskin
         dqh = qmref-qsatg
         tstar = vonkar/fh*dth
         qstar = vonkar/fq*dqh

         thvstar = tstar*(1.+0.61*qmref)+0.61*th*qstar
         zeta  = zldis*vonkar*grav*thvstar/(ustar**2*thv)
         IF (zeta >= 0.) THEN               ! stable
            zeta = min(2.,max(zeta,0.01))
         ELSE                               ! unstable
            zeta = max(-100.,min(zeta,-0.01))
         ENDIF
         obu   = zldis/zeta
         IF (zeta >= 0.) THEN
            um = ur
         ELSE
            wc = beta_wc*(-grav*ustar*thvstar*zcb/thv)**(1./3.)
            um = sqrt(ur*ur+wc*wc)
         ENDIF
         IF (obuold*obu < 0.) nmozsgn = nmozsgn+1
         IF (nmozsgn >= 4) exit
      ENDDO   ! END iter


      !*----------------------------------------------------------------------
      !*----------------------------------------------------------------------
      emis   = emg
      wdm    = um
      ram    = 1./(ustar*ustar/um)
      rah    = 1./(vonkar/fh*ustar)
      raw    = 1./(vonkar/fq*ustar)
      raih   = arhos*cpair/rah
      raiw   = arhos/raw
      shfdts = raih                 ! --
      shfdtl = raiw*qsatgdT         ! --
      shfdt  = shfdts + htvp*shfdtl
      taux   = -arhos*usurf/ram
      tauy   = -arhos*vsurf/ram
      t2m    = thm + vonkar/fh*dth * (fh2m/vonkar - fh/vonkar)
      q2m    = qmref + vonkar/fq*dqh * (fq2m/vonkar - fq/vonkar)
      u10m   = usurf/ur * ustar/vonkar * fm10m
      v10m   = vsurf/ur * ustar/vonkar * fm10m
      zol    = zeta
      rib    = min(5.,zol*ustar**2/(vonkar*vonkar/fh*um**2))
      !*----------------------------------------------------------------------
      !*----------------------------------------------------------------------

   END SUBROUTINE LakStaParms 



   SUBROUTINE LakStaParmsSubs( &
               ! "in" arguments
               ! -------------------
               hwref    , htref    , hqref    , usurf    ,&
               vsurf    , tmref    , qmref    , tskin    ,&
               arhos    , psurf    , zcb      ,& 
               ! "out" arguments
               ! -------------------
               ustar    , qstar    , tstar    , wdm      ,&
               fm10m    , fq10m    , fh2m     , fq2m     ,&
               fm       , fh       , fq       )

! ------------------------ code history -------------------------------------
! Oldname: lakess
! Purpose:          lake surface stability (a subset of lakesp)
!
! Original author : 
!     Xin-Zhong Liang, 2021-11-02
!
! Revisions:
!     Omarjan Obulkasim, 04/2024: Updated variable names and renamed subrountines
! 
!----------------------------------------------------------------------------
   USE MOD_Lake_Const, only :tfrz, hvap, hfus, hsub, grav,&
                                    vonkar,cpair,beta_wc,emisw,rgas
   USE MOD_Qsadv, only : qsadv
   USE MOD_FrictionVelocity, only : moninobukini, moninobuk
!================================================================================
!  ------------------------- input variables ---------------------------
   real(r8), intent(in)     :: &
      hwref                   ,&! observational height of wind [m]
      htref                   ,&! observational height of temperature [m]
      hqref                   ,&! observational height of humidity [m]
      usurf                   ,&! wind component in eastward direction [m/s]
      vsurf                   ,&! wind component in northward direction [m/s]
      tmref                   ,&! temperature at the reference height [K]
      qmref                   ,&! specific humidity at the reference height [kg/kg]
      arhos                   ,&! air density [kg/m3]
      psurf                   ,&! surface pressure [Pa]
      zcb                     ,&! convective boundary height [m]
      tskin                     ! surface temperature [K]

!  ------------------------- output variables ---------------------------
   real(r8), intent(out)    :: &
      ustar                   ,&! u* in similarity theory [m/s]
      qstar                   ,&! q* in similarity theory [kg/kg]
      tstar                   ,&! t* in similarity theory [K]
      wdm                     ,&! lowest level wind speed including the stablity effect [m/s]
      fm10m                   ,&! integral of profile function for momentum at 10m
      fq10m                   ,&! integral of profile function for moisture at 10m
      fh2m                    ,&! integral of profile function for heat at 2m
      fq2m                    ,&! integral of profile function for moisture at 2m
      fm                      ,&! integral of profile function for momentum
      fh                      ,&! integral of profile function for heat
      fq                        ! integral of profile function for moisture
       
!  ------------------------- local variables ---------------------------
   integer                  :: &
      iter                    ,&! iteration index
      nmozsgn                   ! number of times moz changes sign

   real(r8)                 :: &
      eg                      ,&! water vapor pressure at temperature T [pa]
      qsatg                   ,&! saturated humidity [kg/kg]
      qsatgdT                 ,&! d(qsatg)/dT
      degdT                   ,&! d(eg)/dT
      thm                     ,&! potential temperaure at htref with reference to surface pressure [K]
      th                      ,&! atmospheric potential temperaure [K]
      displax                 ,&! zero- displacement height [m]
      dqh                     ,&! diff of humidity between ref. height and surface
      dth                     ,&! diff of virtual temp. between ref. height and surface
      dthv                    ,&! diff of vir. poten. temp. between ref. height and surface
      htvp                    ,&! latent heat of vapor of water (or sublimation) [J/kg]
      obu                     ,&! monin-obukhov length [m]
      obuold                  ,&! monin-obukhov length of previous iteration
      thv                     ,&! virtual potential temperature [K]
      thvstar                 ,&! virtual potential temperature scaling parameter
      um                      ,&! wind speed including the stablity effect [m/s]
      ur                      ,&! wind speed at reference height [m/s]
      tmc                     ,&! temperature at the reference height [C]
      visa                    ,&! kinematic viscosity of dry air [m2/s]
      wc                      ,&! convective velocity [m/s]
      xt                      ,&!
      xq                      ,&!
      zeta                    ,&! dimensionless height used in Monin-Obukhov theory
      zldis                   ,&! reference height minus zero-displacement height [m]
      z0mg                    ,&! roughness length over ground, momentum [m]
      z0hg                    ,&! roughness length over ground, sensible heat [m]
      z0qg                      ! roughness length over ground, latent heat [m]
   integer ii                   ! DO loop or array index
   integer,  parameter :: niters  = 3     ! maximum number of iteration
!================================================================================

      !----------------------------------------------------------------------------
      ! +WMEJ fix a bug : should USE tskin instead of tmref
      IF (tskin > tfrz) THEN
         htvp = hvap
      ELSE
         htvp = hsub
      ENDIF
      ! IF (tmref > tfrz) THEN
      ! htvp = hvap
      ! ELSE
      ! htvp = hsub
      ! ENDIF

      CALL qsadv(tskin,psurf,eg,degdT,qsatg,qsatgdT) !+WMEJ 
      ! qsatg = sathdi(tskin,psurf,qsatgdT)
      thm   = tmref + 0.0098*hwref  ! intermediate variable equivalent to
                  ! tmref*(pgcm/forc_psrf)**(rgas/cpair)
      th    = tmref*(100000./psurf)**(rgas/cpair) ! potential T
      thv   = th*(1.+0.61*qmref)         ! virtual potential T
      ur  = max(0.1,sqrt(usurf*usurf+vsurf*vsurf))

      ! Initialization variables
      nmozsgn = 0
      obuold = 0.
      dth   = thm-tskin
      dqh   = qmref-qsatg
      dthv  = dth*(1.+0.61*qmref)+0.61*th*dqh
      zldis = hwref-0.
      ! aerodynamic roughness
      tmc   = tmref-tfrz
      visa  = 1.326e-5*(1. + (6.542e-3 + (8.301e-6 - 4.84e-9*tmc)*tmc)*tmc)
      ! - Andreas (1989) CRREL Rep. 89-11
      ! loop to obtain initial and good ustar and zo
      ustar = 0.06
      wc    = 0.5
      IF (tskin.ge.tfrz) THEN        ! unfrozen lake
         IF (dthv.ge.0.) THEN
            um = ur
         ELSE
            um = sqrt(ur*ur+wc*wc)
         ENDIF
         DO ii = 1,5
            z0mg = 0.013*ustar*ustar/grav + 0.11*visa/ustar
            ustar = vonkar*um/log(zldis/z0mg)
         ENDDO
      ELSE                        ! frozen lake
         z0mg = 0.04
      ENDIF
      z0qg = z0mg
      z0hg = z0mg
      CALL moninobukini(ur,th,thm,thv,dth,dqh,dthv,zldis,z0mg,um,obu)
      ! CALL moninobukini(ur,thv,dthv,zldis,z0mg,um,obu)
      !----------------------------------------------------------------------------
      DO iter = 1, niters         ! begin stability iteration
         IF (tskin.ge.tfrz) THEN     ! unfrozen lake
            z0mg = 0.013*ustar*ustar/grav + 0.11*visa/ustar
            xq = 2.67*(ustar*z0mg/visa)**0.25 - 2.57
            xt = xq
            z0qg = z0mg/exp(xq)
            z0hg = z0mg/exp(xt)
         ENDIF
         ! Evaluated stability-dependent variables using moz from prior iteration
         displax = 0.
         ! CALL moninobuk(hw,ht,hq,displax,z0mg,z0hg,z0qg,obu,um,&
         ! ustar,temp1,temp2,temp12m,temp22m,&
         ! fh2m,fq2m,f10m,fq10m,fm,fh,fq)
         CALL moninobuk(hwref,htref,hqref,displax,z0mg,z0hg,z0qg,obu,um,&
               ustar,fh2m,fq2m,fm10m,fm,fh,fq)

         obuold = obu
         ! qsatg = sathdi(tskin,psurf,qsatgdT)
         CALL qsadv(tskin,psurf,eg,degdT,qsatg,qsatgdT)
         dth = thm-tskin
         dqh = qmref-qsatg
         tstar = vonkar/fh*dth
         qstar = vonkar/fq*dqh
         thvstar = tstar*(1.+0.61*qmref) + 0.61*th*qstar
         zeta  = zldis*vonkar*grav*thvstar/(ustar**2*thv)
         IF (zeta >= 0.) THEN               ! stable
            zeta = min(2.,max(zeta,0.01))
         ELSE                               ! unstable
            zeta = max(-100.,min(zeta,-0.01))
         ENDIF
         obu   = zldis/zeta
         IF (zeta >= 0.) THEN
            um = ur
         ELSE
            wc = beta_wc*(-grav*ustar*thvstar*zcb/thv)**(1./3.)
            um = sqrt(ur*ur+wc*wc)
         ENDIF
         IF (obuold*obu < 0.) nmozsgn = nmozsgn+1
         IF (nmozsgn >= 4) exit
      ENDDO   ! END iter
      wdm = um

   END SUBROUTINE LakStaParmsSubs


END MODULE  MOD_Lake_Subs