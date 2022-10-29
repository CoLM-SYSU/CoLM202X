#include <define.h>
SUBROUTINE LEAF_interception_CoLM (deltim,dewmx,forc_us,forc_vs,chil,sigf,lai,sai,tair,tleaf,&
                               prc_rain,prc_snow,prl_rain,prl_snow,&
                              ldew,ldew_rain,ldew_snow,z0m,hu,pg_rain,pg_snow,qintr,qintr_rain,qintr_snow)

!=======================================================================
!
! calculation of  interception and drainage of precipitation
! the treatment are based on Sellers et al. (1996)
!
! modified by Yongjiu Dai, 08/31/2002, /04/2014
!
!----------------------------------------------------------------------

  USE precision
  USE PhysicalConstants, only: tfrz
  IMPLICIT NONE

!-----------------------Arguments---------------------------------------

  REAL(r8), intent(in) :: deltim    !seconds in a time step [second]
  REAL(r8), intent(in) :: dewmx     !maximum dew [mm]
  REAL(r8), intent(in) :: forc_us   !wind speed
  REAL(r8), intent(in) :: forc_vs   !wind speed
  REAL(r8), intent(in) :: chil      !leaf angle distribution factor
  REAL(r8), intent(in) :: prc_rain  !convective ranfall [mm/s]
  REAL(r8), intent(in) :: prc_snow  !convective snowfall [mm/s]
  REAL(r8), intent(in) :: prl_rain  !large-scale rainfall [mm/s]
  REAL(r8), intent(in) :: prl_snow  !large-scale snowfall [mm/s]
  REAL(r8), intent(in) :: sigf      !fraction of veg cover, excluding snow-covered veg [-]
  REAL(r8), intent(in) :: lai       !leaf area index [-]
  REAL(r8), intent(in) :: sai       !stem area index [-]
  REAL(r8), intent(in) :: tair     !air temperature [K]
  REAL(r8), intent(in) :: tleaf     !sunlit canopy leaf temperature [K]

  REAL(r8), intent(inout) :: ldew   !depth of water on foliage [mm]
  REAL(r8), intent(inout) :: ldew_rain   !depth of water on foliage [mm]
  REAL(r8), intent(inout) :: ldew_snow   !depth of water on foliage [mm]
  REAL(r8), intent(in) :: z0m            !roughness length
  REAL(r8), intent(in) :: hu             !forcing height of U

  REAL(r8), intent(out) :: pg_rain  !rainfall onto ground including canopy runoff [kg/(m2 s)]
  REAL(r8), intent(out) :: pg_snow  !snowfall onto ground including canopy runoff [kg/(m2 s)]
  REAL(r8), intent(out) :: qintr    !interception [kg/(m2 s)]
  real(r8), INTENT(out) :: qintr_rain ! rainfall interception (mm h2o/s)
  real(r8), INTENT(out) :: qintr_snow ! snowfall interception (mm h2o/s)

!-----------------------Local Variabes---------------------------------

  REAL(r8) :: satcap    !maximum allowed water on canopy [mm]
  REAL(r8) :: satcap_rain    !maximum allowed rain on canopy [mm]
  REAL(r8) :: satcap_snow    !maximum allowed snow on canopy [mm]
  REAL(r8) :: lsai      !sum of leaf area index and stem area index [-]
  REAL(r8) :: chiv      !leaf angle distribution factor
  REAL(r8) :: ppc       !convective precipitation in time-step [mm]
  REAL(r8) :: ppl       !large-scale precipitation in time-step [mm]
  REAL(r8) :: p0        !precipitation in time-step [mm]
  REAL(r8) :: fpi       !coefficient of interception
  REAL(r8) :: fpi_rain      !coefficient of interception of rain
  REAL(r8) :: fpi_snow      !coefficient of interception of snow
  REAL(r8) :: alpha_rain
  REAL(r8) :: alpha_snow
  REAL(r8) :: pinf      !interception of precipitation in time step [mm]
  REAL(r8) :: tti_rain  !direct rain throughfall in time step [mm]
  REAL(r8) :: tti_snow  !direct snow throughfall in time step [mm]
  REAL(r8) :: tex_rain  !canopy rain drainage in time step [mm]
  REAL(r8) :: tex_snow  !canopy snow drainage in time step [mm]
  REAL(r8) :: vegt      !sigf*lsai
  REAL(r8) :: xs        !proportion of the grid area where the intercepted rainfall 
                        !plus the preexisting canopy water storage
  REAL(r8) :: unl_snow_temp,U10,unl_snow_wind,unl_snow

  REAL(r8) :: ap, cp, bp, aa, bb, exrain, arg, alpha, w
  REAL(r8) :: thru_rain, thru_snow
  REAL(r8) :: xsc_rain, xsc_snow
  REAL pcoefs (2,2)

!-----------------------End Variable List-------------------------------
 IF (lai+sai > 1e-6) THEN
    pcoefs(1,1) = 20.
    pcoefs(1,2) = 0.206e-8
    pcoefs(2,1) = 0.0001 
    pcoefs(2,2) = 0.9999 
    bp = 20. 

    lsai   = lai + sai
    vegt   = lsai
    satcap = dewmx*vegt

    p0  = (prc_rain + prc_snow + prl_rain + prl_snow)*deltim
    ppc = (prc_rain+prc_snow)*deltim
    ppl = (prl_rain+prl_snow)*deltim

    w = ldew+p0


    ! 06/08/2019, yuan: why excessed rain calculated here
   IF (tleaf > tfrz) THEN
      xsc_rain = max(0., ldew-satcap)
      xsc_snow = 0.
   ELSE
      xsc_rain = 0.
      xsc_snow = max(0., ldew-satcap)
   ENDIF
   ! 06/08/2019, yuan: ??
   ldew = ldew - (xsc_rain + xsc_snow)


    ap = pcoefs(2,1)
    cp = pcoefs(2,2)

    IF (p0 > 1.e-8) THEN
      ap = ppc/p0 * pcoefs(1,1) + ppl/p0 * pcoefs(2,1)
      cp = ppc/p0 * pcoefs(1,2) + ppl/p0 * pcoefs(2,2)

      !----------------------------------------------------------------------
      !      proportional saturated area (xs) and leaf drainage(tex)
      !-----------------------------------------------------------------------

      chiv = chil
      IF ( abs(chiv) .le. 0.01 ) chiv = 0.01
      aa = 0.5 - 0.633 * chiv - 0.33 * chiv * chiv
      bb = 0.877 * ( 1. - 2. * aa )
      exrain = aa + bb

      ! coefficient of interception
      ! set fraction of potential interception to max 0.25 (Lawrence et al. 2007)
      alpha = 0.25
      fpi = alpha * ( 1.-exp(-exrain*lsai) )
      tti_rain = (prc_rain+prl_rain)*deltim * ( 1.-fpi )
      tti_snow = (prc_snow+prl_snow)*deltim * ( 1.-fpi )

      xs = 1.
      IF (p0*fpi>1.e-9) THEN
         arg = (satcap-ldew)/(p0*fpi*ap) - cp/ap
         IF (arg>1.e-9) THEN
            xs = -1./bp * log( arg )
            xs = min( xs, 1. )
            xs = max( xs, 0. )
         ENDIF
      ENDIF


       ! assume no fall down of the intercepted snowfall in a time step
       ! drainage

      tex_rain = (prc_rain+prl_rain)*deltim * fpi * (ap/bp*(1.-exp(-bp*xs))+cp*xs) - (satcap-ldew) * xs


!       tex_rain = (prc_rain+prl_rain)*deltim * fpi * (ap/bp*(1.-exp(-bp*xs))+cp*xs) &
!                - (satcap-ldew) * xs
      tex_rain = max( tex_rain, 0. )
      tex_snow = 0.

#if(defined CLMDEBUG)
       IF (tex_rain+tex_snow+tti_rain+tti_snow-p0 > 1.e-10) THEN
         write(6,*) 'tex_ + tti_ > p0 in interception code : '
       ENDIF
#endif

    ELSE
       ! all intercepted by canopy leves for very small precipitation
       tti_rain = 0.
       tti_snow = 0.
       tex_rain = 0.
       tex_snow = 0.
    ENDIF

!----------------------------------------------------------------------
!   total throughfall (thru) and store augmentation
!----------------------------------------------------------------------

    thru_rain = tti_rain + tex_rain
    thru_snow = tti_snow + tex_snow
    pinf = p0 - (thru_rain + thru_snow)
    ldew = ldew + pinf

    pg_rain = (xsc_rain + thru_rain) / deltim
    pg_snow = (xsc_snow + thru_snow) / deltim
    qintr   = pinf / deltim


    qintr_rain = prc_rain + prl_rain - thru_rain / deltim
    qintr_snow = prc_snow + prl_snow - thru_snow / deltim


#if(defined CLMDEBUG)
    w = w - ldew - (pg_rain+pg_snow)*deltim
    IF (abs(w) > 1.e-6) THEN
       write(6,*) 'something wrong in interception code : '
       write(6,*) w, ldew, (pg_rain+pg_snow)*deltim, satcap
       CALL abort
    ENDIF
#endif

 ELSE

    ldew = 0.
    pg_rain = prc_rain + prl_rain
    pg_snow = prc_snow + prl_snow
    qintr   = 0.
    qintr_rain = 0.
    qintr_snow = 0.

 ENDIF
 END SUBROUTINE LEAF_interception_CoLM


 SUBROUTINE LEAF_interception_CLM4 (deltim,dewmx,forc_us,forc_vs,chil,sigf,lai,sai,tair,tleaf,&
   prc_rain,prc_snow,prl_rain,prl_snow,&
  ldew,ldew_rain,ldew_snow,z0m,hu,pg_rain,pg_snow,qintr,qintr_rain,qintr_snow)

!=======================================================================
!
! calculation of  interception and drainage of precipitation
! the treatment are based on Sellers et al. (1996)
!
! modified by Yongjiu Dai, 08/31/2002, /04/2014
!
!----------------------------------------------------------------------

USE precision
USE PhysicalConstants, only: tfrz
IMPLICIT NONE

!-----------------------Arguments---------------------------------------

REAL(r8), intent(in) :: deltim    !seconds in a time step [second]
REAL(r8), intent(in) :: dewmx     !maximum dew [mm]
REAL(r8), intent(in) :: forc_us   !wind speed
REAL(r8), intent(in) :: forc_vs   !wind speed
REAL(r8), intent(in) :: chil      !leaf angle distribution factor
REAL(r8), intent(in) :: prc_rain  !convective ranfall [mm/s]
REAL(r8), intent(in) :: prc_snow  !convective snowfall [mm/s]
REAL(r8), intent(in) :: prl_rain  !large-scale rainfall [mm/s]
REAL(r8), intent(in) :: prl_snow  !large-scale snowfall [mm/s]
REAL(r8), intent(in) :: sigf      !fraction of veg cover, excluding snow-covered veg [-]
REAL(r8), intent(in) :: lai       !leaf area index [-]
REAL(r8), intent(in) :: sai       !stem area index [-]
REAL(r8), intent(in) :: tair     !air temperature [K]
REAL(r8), intent(in) :: tleaf     !sunlit canopy leaf temperature [K]

REAL(r8), intent(inout) :: ldew   !depth of water on foliage [mm]
REAL(r8), intent(inout) :: ldew_rain   !depth of water on foliage [mm]
REAL(r8), intent(inout) :: ldew_snow   !depth of water on foliage [mm]
REAL(r8), intent(in) :: z0m            !roughness length
REAL(r8), intent(in) :: hu             !forcing height of U

REAL(r8), intent(out) :: pg_rain  !rainfall onto ground including canopy runoff [kg/(m2 s)]
REAL(r8), intent(out) :: pg_snow  !snowfall onto ground including canopy runoff [kg/(m2 s)]
REAL(r8), intent(out) :: qintr    !interception [kg/(m2 s)]
real(r8), INTENT(out) :: qintr_rain ! rainfall interception (mm h2o/s)
real(r8), INTENT(out) :: qintr_snow ! snowfall interception (mm h2o/s)

!-----------------------Local Variabes---------------------------------

REAL(r8) :: satcap    !maximum allowed water on canopy [mm]
REAL(r8) :: satcap_rain    !maximum allowed rain on canopy [mm]
REAL(r8) :: satcap_snow    !maximum allowed snow on canopy [mm]
REAL(r8) :: lsai      !sum of leaf area index and stem area index [-]
REAL(r8) :: chiv      !leaf angle distribution factor
REAL(r8) :: ppc       !convective precipitation in time-step [mm]
REAL(r8) :: ppl       !large-scale precipitation in time-step [mm]
REAL(r8) :: p0        !precipitation in time-step [mm]
REAL(r8) :: fpi       !coefficient of interception
REAL(r8) :: fpi_rain      !coefficient of interception of rain
REAL(r8) :: fpi_snow      !coefficient of interception of snow
REAL(r8) :: alpha_rain
REAL(r8) :: alpha_snow
REAL(r8) :: pinf      !interception of precipitation in time step [mm]
REAL(r8) :: tti_rain  !direct rain throughfall in time step [mm]
REAL(r8) :: tti_snow  !direct snow throughfall in time step [mm]
REAL(r8) :: tex_rain  !canopy rain drainage in time step [mm]
REAL(r8) :: tex_snow  !canopy snow drainage in time step [mm]
REAL(r8) :: vegt      !sigf*lsai
REAL(r8) :: xs        !proportion of the grid area where the intercepted rainfall 
!plus the preexisting canopy water storage
REAL(r8) :: unl_snow_temp,U10,unl_snow_wind,unl_snow

REAL(r8) :: ap, cp, bp, aa, bb, exrain, arg, alpha, w
REAL(r8) :: thru_rain, thru_snow
REAL(r8) :: xsc_rain, xsc_snow
REAL pcoefs (2,2)

!-----------------------End Variable List-------------------------------
IF (lai+sai > 1e-6) THEN
pcoefs(1,1) = 20.
pcoefs(1,2) = 0.206e-8
pcoefs(2,1) = 0.0001 
pcoefs(2,2) = 0.9999 
bp = 20. 

lsai   = lai + sai
vegt   = lsai
satcap = dewmx*vegt

p0  = (prc_rain + prc_snow + prl_rain + prl_snow)*deltim
ppc = (prc_rain+prc_snow)*deltim
ppl = (prl_rain+prl_snow)*deltim

w = ldew+p0


! 06/08/2019, yuan: why excessed rain calculated here
IF (tleaf > tfrz) THEN
xsc_rain = max(0., ldew-satcap)
xsc_snow = 0.
ELSE
xsc_rain = 0.
xsc_snow = max(0., ldew-satcap)
ENDIF
! 06/08/2019, yuan: ??
ldew = ldew - (xsc_rain + xsc_snow)


ap = pcoefs(2,1)
cp = pcoefs(2,2)

IF (p0 > 1.e-8) THEN
ap = ppc/p0 * pcoefs(1,1) + ppl/p0 * pcoefs(2,1)
cp = ppc/p0 * pcoefs(1,2) + ppl/p0 * pcoefs(2,2)

!----------------------------------------------------------------------
!      proportional saturated area (xs) and leaf drainage(tex)
!-----------------------------------------------------------------------

chiv = chil
IF ( abs(chiv) .le. 0.01 ) chiv = 0.01
aa = 0.5 - 0.633 * chiv - 0.33 * chiv * chiv
bb = 0.877 * ( 1. - 2. * aa )
exrain = aa + bb
exrain =0.5
! coefficient of interception
! set fraction of potential interception to max 0.25 (Lawrence et al. 2007)
alpha = 0.25
fpi = alpha * ( 1.-exp(-exrain*lsai) )
tti_rain = (prc_rain+prl_rain)*deltim * ( 1.-fpi )
tti_snow = (prc_snow+prl_snow)*deltim * ( 1.-fpi )

xs = 1.
IF (p0*fpi>1.e-9) THEN
arg = (satcap-ldew)/(p0*fpi*ap) - cp/ap
IF (arg>1.e-9) THEN
xs = -1./bp * log( arg )
xs = min( xs, 1. )
xs = max( xs, 0. )
ENDIF
ENDIF


! assume no fall down of the intercepted snowfall in a time step
! drainage

tex_rain = (prc_rain+prl_rain)*deltim * fpi * (ap/bp*(1.-exp(-bp*xs))+cp*xs) - (satcap-ldew) * xs


!       tex_rain = (prc_rain+prl_rain)*deltim * fpi * (ap/bp*(1.-exp(-bp*xs))+cp*xs) &
!                - (satcap-ldew) * xs
tex_rain = max( tex_rain, 0. )
tex_snow = 0.

#if(defined CLMDEBUG)
IF (tex_rain+tex_snow+tti_rain+tti_snow-p0 > 1.e-10) THEN
write(6,*) 'tex_ + tti_ > p0 in interception code : '
ENDIF
#endif

ELSE
! all intercepted by canopy leves for very small precipitation
tti_rain = 0.
tti_snow = 0.
tex_rain = 0.
tex_snow = 0.
ENDIF

!----------------------------------------------------------------------
!   total throughfall (thru) and store augmentation
!----------------------------------------------------------------------

thru_rain = tti_rain + tex_rain
thru_snow = tti_snow + tex_snow
pinf = p0 - (thru_rain + thru_snow)
ldew = ldew + pinf

pg_rain = (xsc_rain + thru_rain) / deltim
pg_snow = (xsc_snow + thru_snow) / deltim
qintr   = pinf / deltim


qintr_rain = prc_rain + prl_rain - thru_rain / deltim
qintr_snow = prc_snow + prl_snow - thru_snow / deltim


#if(defined CLMDEBUG)
w = w - ldew - (pg_rain+pg_snow)*deltim
IF (abs(w) > 1.e-6) THEN
write(6,*) 'something wrong in interception code : '
write(6,*) w, ldew, (pg_rain+pg_snow)*deltim, satcap
CALL abort
ENDIF
#endif

ELSE

ldew = 0.
pg_rain = prc_rain + prl_rain
pg_snow = prc_snow + prl_snow
qintr   = 0.
qintr_rain = 0.
qintr_snow = 0.

ENDIF
END SUBROUTINE LEAF_interception_CLM4


 SUBROUTINE LEAF_interception_CLM5 (deltim,dewmx,forc_us,forc_vs,chil,sigf,lai,sai,tair,tleaf,&
   prc_rain,prc_snow,prl_rain,prl_snow,&
  ldew,ldew_rain,ldew_snow,z0m,hu,pg_rain,pg_snow,qintr,qintr_rain,qintr_snow)

!=======================================================================
!
! calculation of  interception and drainage of precipitation
! the treatment are based on Sellers et al. (1996)
!
! modified by Yongjiu Dai, 08/31/2002, /04/2014
!
!----------------------------------------------------------------------

USE precision
USE PhysicalConstants, only: tfrz
IMPLICIT NONE

!-----------------------Arguments---------------------------------------

REAL(r8), intent(in) :: deltim    !seconds in a time step [second]
REAL(r8), intent(in) :: dewmx     !maximum dew [mm]
REAL(r8), intent(in) :: forc_us   !wind speed
REAL(r8), intent(in) :: forc_vs   !wind speed
REAL(r8), intent(in) :: chil      !leaf angle distribution factor
REAL(r8), intent(in) :: prc_rain  !convective ranfall [mm/s]
REAL(r8), intent(in) :: prc_snow  !convective snowfall [mm/s]
REAL(r8), intent(in) :: prl_rain  !large-scale rainfall [mm/s]
REAL(r8), intent(in) :: prl_snow  !large-scale snowfall [mm/s]
REAL(r8), intent(in) :: sigf      !fraction of veg cover, excluding snow-covered veg [-]
REAL(r8), intent(in) :: lai       !leaf area index [-]
REAL(r8), intent(in) :: sai       !stem area index [-]
REAL(r8), intent(in) :: tair      !air temperature [K]
REAL(r8), intent(in) :: tleaf     !sunlit canopy leaf temperature [K]

REAL(r8), intent(inout) :: ldew   !depth of water on foliage [mm]
REAL(r8), intent(inout) :: ldew_rain   !depth of water on foliage [mm]
REAL(r8), intent(inout) :: ldew_snow   !depth of water on foliage [mm]
REAL(r8), intent(in) :: z0m            !roughness length
REAL(r8), intent(in) :: hu             !forcing height of U

REAL(r8), intent(out) :: pg_rain  !rainfall onto ground including canopy runoff [kg/(m2 s)]
REAL(r8), intent(out) :: pg_snow  !snowfall onto ground including canopy runoff [kg/(m2 s)]
REAL(r8), intent(out) :: qintr    !interception [kg/(m2 s)]
real(r8), INTENT(out) :: qintr_rain ! rainfall interception (mm h2o/s)
real(r8), INTENT(out) :: qintr_snow ! snowfall interception (mm h2o/s)

!-----------------------Local Variabes---------------------------------

REAL(r8) :: satcap    !maximum allowed water on canopy [mm]
REAL(r8) :: satcap_rain    !maximum allowed rain on canopy [mm]
REAL(r8) :: satcap_snow    !maximum allowed snow on canopy [mm]
REAL(r8) :: lsai      !sum of leaf area index and stem area index [-]
REAL(r8) :: chiv      !leaf angle distribution factor
REAL(r8) :: ppc       !convective precipitation in time-step [mm]
REAL(r8) :: ppl       !large-scale precipitation in time-step [mm]
REAL(r8) :: p0        !precipitation in time-step [mm]
REAL(r8) :: fpi       !coefficient of interception
REAL(r8) :: fpi_rain      !coefficient of interception of rain
REAL(r8) :: fpi_snow      !coefficient of interception of snow
REAL(r8) :: alpha_rain
REAL(r8) :: alpha_snow
REAL(r8) :: pinf      !interception of precipitation in time step [mm]
REAL(r8) :: tti_rain  !direct rain throughfall in time step [mm]
REAL(r8) :: tti_snow  !direct snow throughfall in time step [mm]
REAL(r8) :: tex_rain  !canopy rain drainage in time step [mm]
REAL(r8) :: tex_snow  !canopy snow drainage in time step [mm]
REAL(r8) :: vegt      !sigf*lsai
REAL(r8) :: xs        !proportion of the grid area where the intercepted rainfall 
!plus the preexisting canopy water storage
REAL(r8) :: unl_snow_temp,U10,unl_snow_wind,unl_snow

REAL(r8) :: ap, cp, bp, aa, bb, exrain, arg, alpha, w
REAL(r8) :: thru_rain, thru_snow
REAL(r8) :: xsc_rain, xsc_snow
REAL pcoefs (2,2)

!-----------------------End Variable List-------------------------------
IF (lai+sai > 1e-6) THEN
lsai   = lai + sai
vegt   = lsai
p0  = (prc_rain + prc_snow + prl_rain + prl_snow)*deltim
ppc = (prc_rain+prc_snow)*deltim
ppl = (prl_rain+prl_snow)*deltim
w = ldew+p0
satcap_rain = dewmx*vegt
satcap_snow = satcap_rain*60.0

xsc_rain      = max(0., ldew_rain-satcap_rain)
xsc_snow      = max(0., ldew_snow-satcap_snow)
!xsc_snow      = min(ldew_snow, xsc_snow)
!ldew_rain     = ldew_rain - xsc_rain
!ldew_snow     = ldew_snow - xsc_snow
!ldew          = ldew_rain+ldew_snow
!xsc_rain=0.0

!unload due to wind and temperature
U10= sqrt(forc_us*forc_us+forc_vs*forc_vs)*log(10.0/z0m)/log(hu/z0m)
unl_snow_temp = ldew_snow*(tleaf-270.0)/(1.87*1.e5)
unl_snow_temp =max(unl_snow_temp,0.0)
unl_snow_wind = U10*ldew_snow/(1.56*1.e5)
unl_snow      = min(unl_snow_temp+unl_snow_wind,ldew_snow)

xsc_snow=min(xsc_snow+xsc_snow,ldew_snow)

IF(p0 > 1.e-8) THEN
alpha_rain = 1.0
alpha_snow = 1.0
fpi_rain   = alpha_rain * tanh(lsai)  
fpi_snow   = alpha_snow * ( 1.-exp(-0.5*lsai) )
tti_rain   = (prc_rain+prl_rain)*deltim * ( 1.-fpi_rain )
tti_snow   = (prc_snow+prl_snow)*deltim * ( 1.-fpi_snow )
tex_rain = (prc_rain+prl_rain)*deltim * fpi_rain -satcap_rain         !*(prc_rain+prl_rain)/p0 !(satcap-ldew) * xs
tex_snow = (prc_snow+prl_snow)*deltim * fpi_snow -satcap_snow         ! (ap/bp*(1.-exp(-bp*xs))+cp*xs) - (satcap-ldew) * xs
tex_rain = max( tex_rain, 0. )
tex_snow = max( tex_snow, 0. )

#if(defined CLMDEBUG)
IF (tex_rain+tex_snow+tti_rain+tti_snow-p0 > 1.e-10) THEN
write(6,*) 'tex_ + tti_ > p0 in interception code : '
ENDIF
#endif
ELSE
! all intercepted by canopy leves for very small precipitation
tti_rain = 0.
tti_snow = 0.
tex_rain = 0.
tex_snow = 0.
ENDIF



!----------------------------------------------------------------------
!   total throughfall (thru) and store augmentation
!----------------------------------------------------------------------

thru_rain = tti_rain + tex_rain
thru_snow = tti_snow + tex_snow
pinf = p0 - (thru_rain + thru_snow)
ldew_rain=ldew_rain+ (prc_rain + prl_rain)*deltim - thru_rain -xsc_rain
ldew_snow=ldew_snow+ (prc_snow + prl_snow)*deltim  - thru_snow-xsc_snow
ldew_snow=max(0.0,ldew_snow)
ldew_rain=max(0.0,ldew_rain)
ldew = ldew_rain+ldew_snow !+ pinf

pg_rain = (xsc_rain + thru_rain) / deltim
pg_snow = (xsc_snow + thru_snow) / deltim
qintr   = pinf / deltim
qintr_rain = prc_rain + prl_rain - thru_rain / deltim
qintr_snow = prc_snow + prl_snow - thru_snow / deltim

#if(defined CLMDEBUG)
w = w - ldew - (pg_rain+pg_snow)*deltim
IF (abs(w) > 1.e-6) THEN
write(6,*) 'something wrong in interception code : '
write(6,*) 'w, ldew, (pg_rain+pg_snow)*deltim, satcap_rain,satcap_snow:',w, ldew, (pg_rain+pg_snow)*deltim, satcap_rain,satcap_snow
CALL abort
ENDIF
#endif

ELSE

ldew = 0.
ldew_rain = 0.
ldew_snow = 0.
pg_rain = prc_rain + prl_rain
pg_snow = prc_snow + prl_snow
qintr   = 0.
qintr_rain = 0.
qintr_snow = 0.

ENDIF

END SUBROUTINE LEAF_interception_CLM5

SUBROUTINE LEAF_interception_noahmp (deltim,dewmx,forc_us,forc_vs,chil,sigf,lai,sai,tair,tleaf, &
   prc_rain,prc_snow,prl_rain,prl_snow,&         
  ldew,ldew_rain,ldew_snow,z0m,hu,pg_rain,pg_snow,qintr,qintr_rain,qintr_snow)
!=======================================================================
!
! calculation of  interception and drainage of precipitation
! the treatment are based on Sellers et al. (1996)
!
! modified by Yongjiu Dai, 08/31/2002, /04/2014
!
!----------------------------------------------------------------------

USE precision
USE PhysicalConstants, only: tfrz
IMPLICIT NONE

!-----------------------Arguments---------------------------------------

REAL(r8), intent(in) :: deltim    !seconds in a time step [second]
REAL(r8), intent(in) :: dewmx     !maximum dew [mm]
REAL(r8), intent(in) :: forc_us   !wind speed
REAL(r8), intent(in) :: forc_vs   !wind speed
REAL(r8), intent(in) :: chil      !leaf angle distribution factor
REAL(r8), intent(in) :: prc_rain  !convective ranfall [mm/s]
REAL(r8), intent(in) :: prc_snow  !convective snowfall [mm/s]
REAL(r8), intent(in) :: prl_rain  !large-scale rainfall [mm/s]
REAL(r8), intent(in) :: prl_snow  !large-scale snowfall [mm/s]
REAL(r8), intent(in) :: sigf      !fraction of veg cover, excluding snow-covered veg [-]
REAL(r8), intent(in) :: lai       !leaf area index [-]
REAL(r8), intent(in) :: sai       !stem area index [-]
REAL(r8), intent(in) :: tair     !air temperature [K]
REAL(r8), intent(inout) :: tleaf   !sunlit canopy leaf temperature [K]

REAL(r8), intent(inout) :: ldew   !depth of water on foliage [mm]
REAL(r8), intent(inout) :: ldew_rain   !depth of liquid on foliage [mm]
REAL(r8), intent(inout) :: ldew_snow   !depth of liquid on foliage [mm]
REAL(r8), intent(in) :: z0m            !roughness length
REAL(r8), intent(in) :: hu             !forcing height of U

REAL(r8), intent(out) :: pg_rain  !rainfall onto ground including canopy runoff [kg/(m2 s)]
REAL(r8), intent(out) :: pg_snow  !snowfall onto ground including canopy runoff [kg/(m2 s)]
REAL(r8), intent(out) :: qintr    !interception [kg/(m2 s)]
real(r8), INTENT(out) :: qintr_rain ! rainfall interception (mm h2o/s)
real(r8), INTENT(out) :: qintr_snow ! snowfall interception (mm h2o/s)

!-----------------------Local Variabes---------------------------------

REAL(r8) :: satcap_rain    !maximum allowed liquid water on canopy [mm]
REAL(r8) :: satcap_snow    !maximum allowed solid water on canopy [mm]
REAL(r8) :: fveg           !vegetation fraction
REAL(r8) :: FT             !The temperature factor for snow unloading 
REAL(r8) :: FV             !The wind factor for snow unloading 
REAL(r8) :: ICEDRIP        ! snow unloading


REAL(r8) :: lsai      !sum of leaf area index and stem area index [-]
REAL(r8) :: chiv      !leaf angle distribution factor
REAL(r8) :: ppc       !convective precipitation in time-step [mm]
REAL(r8) :: ppl       !large-scale precipitation in time-step [mm]
REAL(r8) :: p0        !precipitation in time-step [mm]
REAL(r8) :: fpi       !coefficient of interception
REAL(r8) :: pinf      !interception of precipitation in time step [mm]
REAL(r8) :: tti_rain  !direct rain throughfall in time step [mm]
REAL(r8) :: tti_snow  !direct snow throughfall in time step [mm]
REAL(r8) :: tex_rain  !canopy rain drainage in time step [mm]
REAL(r8) :: tex_snow  !canopy snow drainage in time step [mm]
REAL(r8) :: vegt      !sigf*lsai
REAL(r8) :: xs        !proportion of the grid area where the intercepted rainfall 
!plus the preexisting canopy water storage
REAL(r8) :: unl_snow_temp,U10,unl_snow_wind,unl_snow

REAL(r8) :: ap, cp, bp, aa, bb, exrain, arg, alpha, w
REAL(r8) :: thru_rain, thru_snow
REAL(r8) :: xsc_rain, xsc_snow
REAL pcoefs (2,2)
REAL(r8) :: CICE,DENICE,HFUS,ldew_smelt,ldew_frzc,CWAT,DENH2O,FP,int_rain,int_snow
!-----------------------End Variable List-------------------------------
CICE   = 2.094E06  !specific heat capacity of ice (j/m3/k)
DENICE = 917.      !density of ice (kg/m3)
HFUS   = 0.3336E06 !latent heat of fusion (j/kg)
CWAT   = 4.188E06  !specific heat capacity of water (j/m3/k)
DENH2O = 1000.     !density of water (kg/m3)
IF (lai+sai > 1e-6) THEN

   lsai   = lai + sai
   vegt   = lsai
   satcap_rain = dewmx*vegt
   satcap_snow = satcap_rain*60.0
   fveg=max(0.05,1.0-exp(-0.52*lsai))

   !snow unloading
   if (ldew_snow>1.e-8) then 
      FT = MAX(0.0,(tair - 270.15) / 1.87E5)
      FV = SQRT(forc_us*forc_us + forc_vs*forc_vs) / 1.56E5
      ICEDRIP = MAX(0.,ldew_snow) * (FV+FT)    !MB: removed /DT
   else
      ICEDRIP = 0.
   endif

   ! phase change and excess !
   IF (tleaf > tfrz) THEN
      ldew_smelt = MIN(ldew_snow,(tleaf-tfrz)*CICE*ldew_snow/DENICE/(HFUS))
      ldew_smelt = max(ldew_smelt,0.0)
      ldew_snow  = ldew_snow-ldew_smelt
      ldew_rain  = ldew_rain+ldew_smelt
      xsc_rain   = max(0., ldew_rain-satcap_rain)
      xsc_snow   = ICEDRIP+max(0., ldew_snow-satcap_snow)
      ! tleaf      = fveg*tfrz+ (1.0-fwet)*tleaf
   ELSE
      ldew_frzc  = MIN(ldew_rain,(tfrz-tleaf)*CWAT*ldew_rain/DENH2O/(HFUS))
      ldew_frzc  = max(ldew_smelt,0.0)
      ldew_snow  = ldew_snow+ldew_frzc
      ldew_rain  = ldew_rain-ldew_frzc
      xsc_rain   = max(0., ldew_rain-satcap_rain)
      xsc_snow   = ICEDRIP+max(0., ldew_snow-satcap_snow)
      !tleaf      = fveg*tfrz+ (1.0-fwet)*tleaf
   ENDIF

      
   IF (p0 > 1.e-8) THEN
      p0  = (prc_rain + prc_snow + prl_rain + prl_snow)*deltim
      ppc = (prc_rain+prc_snow)*deltim
      ppl = (prl_rain+prl_snow)*deltim

      w = ldew+p0

      tti_rain = (prc_rain+prl_rain)*deltim * ( 1.-fveg )
      tti_snow = (prc_snow+prl_snow)*deltim * ( 1.-fveg )

      FP=p0/(10.*ppc+ppl)
      int_rain=min(fveg*FP,(satcap_rain-ldew_rain)/((prc_rain+prl_rain)*deltim)*(1.0-exp(-(prc_rain+prl_rain)*deltim/satcap_rain)))
      int_snow=min(fveg*FP,(satcap_snow-ldew_snow)/((prc_snow+prl_snow)*deltim)*(1.0-exp(-(prc_snow+prl_snow)*deltim/satcap_snow)))
      int_rain=max(0.,int_rain)
      int_snow=max(0.,int_snow)

      tex_rain = (prc_rain+prl_rain)*deltim * ( 1. - int_rain )
      tex_snow = (prc_snow+prl_snow)*deltim * ( 1. - int_snow )
#if(defined CLMDEBUG)
   IF (tex_rain+tex_snow+tti_rain+tti_snow-p0 > 1.e-10) THEN
      write(6,*) 'tex_ + tti_ > p0 in interception code : '
   ENDIF
#endif

   else
      ! all intercepted by canopy leves for very small precipitation
      tti_rain = 0.
      tti_snow = 0.
      tex_rain = 0.
      tex_snow = 0.
   endif

   !BDFALL = 67.92+51.25*EXP(MIN(2.5,(SFCTMP-TFRZ))/2.59)

!----------------------------------------------------------------------
!   total throughfall (thru) and store augmentation
!----------------------------------------------------------------------

   thru_rain = tti_rain + tex_rain
   thru_snow = tti_snow + tex_snow
   pinf = p0 - (thru_rain + thru_snow)
   ldew = ldew + pinf
   
   pg_rain = (xsc_rain + thru_rain) / deltim
   pg_snow = (xsc_snow + thru_snow) / deltim
   qintr   = pinf / deltim
   
   
   qintr_rain = prc_rain + prl_rain - thru_rain / deltim
   qintr_snow = prc_snow + prl_snow - thru_snow / deltim
   
   
#if(defined CLMDEBUG)
   w = w - ldew - (pg_rain+pg_snow)*deltim
   IF (abs(w) > 1.e-6) THEN
      write(6,*) 'something wrong in interception code : '
      write(6,*) w, ldew, (pg_rain+pg_snow)*deltim  !, satcap
      CALL abort
   ENDIF
#endif
   
ELSE
   
   ldew = 0.
   pg_rain = prc_rain + prl_rain
   pg_snow = prc_snow + prl_snow
   qintr   = 0.
   qintr_rain = 0.
   qintr_snow = 0.
ENDIF

END SUBROUTINE LEAF_interception_noahmp


SUBROUTINE LEAF_interception_matsiro (deltim,dewmx,forc_us,forc_vs,chil,sigf,lai,sai,tair,tleaf, &
   prc_rain,prc_snow,prl_rain,prl_snow,&         
  ldew,ldew_rain,ldew_snow,z0m,hu,pg_rain,pg_snow,qintr,qintr_rain,qintr_snow)
!=======================================================================
!
! calculation of  interception and drainage of precipitation
! the treatment are based on Sellers et al. (1996)
!
! modified by Yongjiu Dai, 08/31/2002, /04/2014
!
!----------------------------------------------------------------------

USE precision
USE PhysicalConstants, only: tfrz
IMPLICIT NONE

!-----------------------Arguments---------------------------------------

REAL(r8), intent(in) :: deltim    !seconds in a time step [second]
REAL(r8), intent(in) :: dewmx     !maximum dew [mm]
REAL(r8), intent(in) :: forc_us   !wind speed
REAL(r8), intent(in) :: forc_vs   !wind speed
REAL(r8), intent(in) :: chil      !leaf angle distribution factor
REAL(r8), intent(in) :: prc_rain  !convective ranfall [mm/s]
REAL(r8), intent(in) :: prc_snow  !convective snowfall [mm/s]
REAL(r8), intent(in) :: prl_rain  !large-scale rainfall [mm/s]
REAL(r8), intent(in) :: prl_snow  !large-scale snowfall [mm/s]
REAL(r8), intent(in) :: sigf      !fraction of veg cover, excluding snow-covered veg [-]
REAL(r8), intent(in) :: lai       !leaf area index [-]
REAL(r8), intent(in) :: sai       !stem area index [-]
REAL(r8), intent(in) :: tair     !air temperature [K]
REAL(r8), intent(inout) :: tleaf   !sunlit canopy leaf temperature [K]

REAL(r8), intent(inout) :: ldew   !depth of water on foliage [mm]
REAL(r8), intent(inout) :: ldew_rain   !depth of liquid on foliage [mm]
REAL(r8), intent(inout) :: ldew_snow   !depth of liquid on foliage [mm]
REAL(r8), intent(in) :: z0m            !roughness length
REAL(r8), intent(in) :: hu             !forcing height of U


REAL(r8), intent(out) :: pg_rain  !rainfall onto ground including canopy runoff [kg/(m2 s)]
REAL(r8), intent(out) :: pg_snow  !snowfall onto ground including canopy runoff [kg/(m2 s)]
REAL(r8), intent(out) :: qintr    !interception [kg/(m2 s)]
real(r8), INTENT(out) :: qintr_rain ! rainfall interception (mm h2o/s)
real(r8), INTENT(out) :: qintr_snow ! snowfall interception (mm h2o/s)

!-----------------------Local Variabes---------------------------------

REAL(r8) :: satcap_rain    !maximum allowed liquid water on canopy [mm]
REAL(r8) :: satcap_snow    !maximum allowed solid water on canopy [mm]
REAL(r8) :: fveg           !vegetation fraction
REAL(r8) :: FT             !The temperature factor for snow unloading 
REAL(r8) :: FV             !The wind factor for snow unloading 
REAL(r8) :: ICEDRIP        ! snow unloading


REAL(r8) :: lsai      !sum of leaf area index and stem area index [-]
REAL(r8) :: chiv      !leaf angle distribution factor
REAL(r8) :: ppc       !convective precipitation in time-step [mm]
REAL(r8) :: ppl       !large-scale precipitation in time-step [mm]
REAL(r8) :: p0        !precipitation in time-step [mm]
REAL(r8) :: fpi       !coefficient of interception
REAL(r8) :: pinf      !interception of precipitation in time step [mm]
REAL(r8) :: tti_rain  !direct rain throughfall in time step [mm]
REAL(r8) :: tti_snow  !direct snow throughfall in time step [mm]
REAL(r8) :: tex_rain  !canopy rain drainage in time step [mm]
REAL(r8) :: tex_snow  !canopy snow drainage in time step [mm]
REAL(r8) :: vegt      !sigf*lsai
REAL(r8) :: xs        !proportion of the grid area where the intercepted rainfall 
!plus the preexisting canopy water storage
REAL(r8) :: unl_snow_temp,U10,unl_snow_wind,unl_snow

REAL(r8) :: ap, cp, bp, aa, bb, exrain, arg, alpha, w
REAL(r8) :: thru_rain, thru_snow
REAL(r8) :: xsc_rain, xsc_snow
REAL pcoefs (2,2)
REAL(r8) :: CICE,DENICE,HFUS,ldew_smelt,ldew_frzc,CWAT,DENH2O,fpi_rain,fpi_snow
!-----------------------End Variable List-------------------------------
CICE   = 2.094E06  !specific heat capacity of ice (j/m3/k)
DENICE = 917.      !density of ice (kg/m3)
HFUS   = 0.3336E06 !latent heat of fusion (j/kg)
CWAT   = 4.188E06  !specific heat capacity of water (j/m3/k)
DENH2O = 1000.     !density of water (kg/m3)
IF (lai+sai > 1e-6) THEN

   lsai   = lai + sai
   vegt   = lsai
   satcap_rain = 0.2*vegt
   satcap_snow = 0.2*vegt


   !fveg=max(0.05,1.0-exp(-0.52*lsai))

   ! phase change and excess !
   IF (tleaf > tfrz) THEN
      ldew_smelt = MIN(ldew_snow,(tleaf-tfrz)*CICE*ldew_snow/DENICE/(HFUS))
      ldew_smelt = max(ldew_smelt,0.0)
      ldew_snow  = ldew_snow-ldew_smelt
      ldew_rain  = ldew_rain+ldew_smelt
      xsc_rain   = max(0., ldew_rain-satcap_rain)
      xsc_snow   = max(0., ldew_snow-satcap_snow)
   ELSE
      ldew_frzc  = MIN(ldew_rain,(tfrz-tleaf)*CWAT*ldew_rain/DENH2O/(HFUS))
      ldew_frzc  = max(ldew_smelt,0.0)
      ldew_snow  = ldew_snow+ldew_frzc
      ldew_rain  = ldew_rain-ldew_frzc
      xsc_rain   = max(0., ldew_rain-satcap_rain)
      xsc_snow   = max(0., ldew_snow-satcap_snow)
   ENDIF

      
   IF (p0 > 1.e-8) THEN
      p0  = (prc_rain + prc_snow + prl_rain + prl_snow)*deltim
      ppc = (prc_rain+prc_snow)*deltim
      ppl = (prl_rain+prl_snow)*deltim

      w = ldew+p0
      fpi_rain   = max(min(lsai, 1.0),0.0)
      fpi_snow   = max(min(lsai, 1.0),0.0)

      tti_rain = fpi_rain * (prc_rain/0.1+prl_rain)*deltim + fpi_rain * (prl_rain)*deltim
      tti_snow = fpi_snow * (prc_snow/0.1+prl_snow)*deltim + fpi_snow * (prl_rain)*deltim
      tti_rain = min(tti_rain,(prc_rain+prl_rain)*deltim)
      tti_snow = min(tti_snow,(prc_snow+prl_snow)*deltim)

      tex_rain=max(ldew_rain+(prc_rain+prl_rain)*deltim-tti_rain-satcap_rain,0.0) + (1.14d-11)*exp(3.7d3*(min(ldew_rain+(prc_rain+prl_rain)*deltim-tti_rain,satcap_rain)/deltim))*deltim
      tex_snow=max(ldew_snow+(prc_snow+prl_snow)*deltim-tti_snow-satcap_snow,0.0) + (1.14d-11)*exp(3.7d3*(min(ldew_snow+(prc_snow+prl_snow)*deltim-tti_snow,satcap_snow)/deltim))*deltim

      
#if(defined CLMDEBUG)
   IF (tex_rain+tex_snow+tti_rain+tti_snow-p0 > 1.e-10) THEN
      write(6,*) 'tex_ + tti_ > p0 in interception code : '
   ENDIF
#endif

   else
      ! all intercepted by canopy leves for very small precipitation
      tti_rain = 0.
      tti_snow = 0.
      tex_rain = 0.
      tex_snow = 0.
   endif

   !BDFALL = 67.92+51.25*EXP(MIN(2.5,(SFCTMP-TFRZ))/2.59)

!----------------------------------------------------------------------
!   total throughfall (thru) and store augmentation
!----------------------------------------------------------------------
   
   thru_rain = tti_rain + tex_rain
   thru_snow = tti_snow + tex_snow
   pinf = p0 - (thru_rain + thru_snow)
   ldew = ldew + pinf
   ldew_rain= ldew_rain+(prc_rain+prl_rain)*deltim- thru_rain 
   ldew_snow= ldew_snow+(prc_snow+prl_snow)*deltim- thru_snow 


   pg_rain = (xsc_rain + thru_rain) / deltim
   pg_snow = (xsc_snow + thru_snow) / deltim
   qintr   = pinf / deltim

   qintr_rain = prc_rain + prl_rain - thru_rain / deltim
   qintr_snow = prc_snow + prl_snow - thru_snow / deltim
#if(defined CLMDEBUG)
   w = w - ldew - (pg_rain+pg_snow)*deltim
   IF (abs(w) > 1.e-6) THEN
      write(6,*) 'something wrong in interception code : '
      write(6,*) w, ldew, (pg_rain+pg_snow)*deltim !, satcap
      CALL abort
   ENDIF
#endif
   
ELSE
   
   ldew = 0.
   ldew_rain = 0.
   ldew_snow = 0.
   pg_rain = prc_rain + prl_rain
   pg_snow = prc_snow + prl_snow
   qintr   = 0.
   qintr_rain = 0.
   qintr_snow = 0.
ENDIF

END SUBROUTINE LEAF_interception_matsiro


SUBROUTINE LEAF_interception_vic (deltim,dewmx,forc_us,forc_vs,chil,sigf,lai,sai,tair,tleaf, &
   prc_rain,prc_snow,prl_rain,prl_snow,&         
  ldew,ldew_rain,ldew_snow,z0m,hu,pg_rain,pg_snow,qintr,qintr_rain,qintr_snow)
!=======================================================================
!
! calculation of  interception and drainage of precipitation
! the treatment are based on Sellers et al. (1996)
!
! modified by Yongjiu Dai, 08/31/2002, /04/2014
!
!----------------------------------------------------------------------

USE precision
USE PhysicalConstants, only: tfrz
IMPLICIT NONE

!-----------------------Arguments---------------------------------------

REAL(r8), intent(in) :: deltim    !seconds in a time step [second]
REAL(r8), intent(in) :: dewmx     !maximum dew [mm]
REAL(r8), intent(in) :: forc_us   !wind speed
REAL(r8), intent(in) :: forc_vs   !wind speed
REAL(r8), intent(in) :: chil      !leaf angle distribution factor
REAL(r8), intent(in) :: prc_rain  !convective ranfall [mm/s]
REAL(r8), intent(in) :: prc_snow  !convective snowfall [mm/s]
REAL(r8), intent(in) :: prl_rain  !large-scale rainfall [mm/s]
REAL(r8), intent(in) :: prl_snow  !large-scale snowfall [mm/s]
REAL(r8), intent(in) :: sigf      !fraction of veg cover, excluding snow-covered veg [-]
REAL(r8), intent(in) :: lai       !leaf area index [-]
REAL(r8), intent(in) :: sai       !stem area index [-]
REAL(r8), intent(in) :: tair     !air temperature [K]
REAL(r8), intent(inout) :: tleaf   !sunlit canopy leaf temperature [K]

REAL(r8), intent(inout) :: ldew   !depth of water on foliage [mm]
REAL(r8), intent(inout) :: ldew_rain   !depth of liquid on foliage [mm]
REAL(r8), intent(inout) :: ldew_snow   !depth of liquid on foliage [mm]
REAL(r8), intent(in) :: z0m            !roughness length
REAL(r8), intent(in) :: hu             !forcing height of U


REAL(r8), intent(out) :: pg_rain  !rainfall onto ground including canopy runoff [kg/(m2 s)]
REAL(r8), intent(out) :: pg_snow  !snowfall onto ground including canopy runoff [kg/(m2 s)]
REAL(r8), intent(out) :: qintr    !interception [kg/(m2 s)]
real(r8), INTENT(out) :: qintr_rain ! rainfall interception (mm h2o/s)
real(r8), INTENT(out) :: qintr_snow ! snowfall interception (mm h2o/s)

!-----------------------Local Variabes---------------------------------

REAL(r8) :: satcap_rain    !maximum allowed liquid water on canopy [mm]
REAL(r8) :: satcap_snow    !maximum allowed solid water on canopy [mm]
REAL(r8) :: fveg           !vegetation fraction
REAL(r8) :: FT             !The temperature factor for snow unloading 
REAL(r8) :: FV             !The wind factor for snow unloading 
REAL(r8) :: ICEDRIP        ! snow unloading


REAL(r8) :: lsai      !sum of leaf area index and stem area index [-]
REAL(r8) :: chiv      !leaf angle distribution factor
REAL(r8) :: ppc       !convective precipitation in time-step [mm]
REAL(r8) :: ppl       !large-scale precipitation in time-step [mm]
REAL(r8) :: p0        !precipitation in time-step [mm]
REAL(r8) :: fpi       !coefficient of interception
REAL(r8) :: pinf      !interception of precipitation in time step [mm]
REAL(r8) :: tti_rain  !direct rain throughfall in time step [mm]
REAL(r8) :: tti_snow  !direct snow throughfall in time step [mm]
REAL(r8) :: tex_rain  !canopy rain drainage in time step [mm]
REAL(r8) :: tex_snow  !canopy snow drainage in time step [mm]
REAL(r8) :: vegt      !sigf*lsai
REAL(r8) :: xs        !proportion of the grid area where the intercepted rainfall 
!plus the preexisting canopy water storage
REAL(r8) :: unl_snow_temp,U10,unl_snow_wind,unl_snow

REAL(r8) :: ap, cp, bp, aa, bb, exrain, arg, alpha, w
REAL(r8) :: thru_rain, thru_snow
REAL(r8) :: xsc_rain, xsc_snow
REAL pcoefs (2,2)
REAL(r8) :: CICE,DENICE,HFUS,CWAT,DENH2O
real(r8) :: Imax1,Lr,ldew_max_snow,Snow,Rain,DeltaSnowInt,Wind,BlownSnow,SnowThroughFall
real(r8) :: MaxInt,MaxWaterInt,RainThroughFall,Overload,IntRainFract,IntSnowFract,ldew_smelt
real(r8) :: ldew_frzc,drip
!-----------------------End Variable List-------------------------------
CICE   = 2.094E06  !specific heat capacity of ice (j/m3/k)
DENICE = 917.      !density of ice (kg/m3)
HFUS   = 0.3336E06 !latent heat of fusion (j/kg)
CWAT   = 4.188E06  !specific heat capacity of water (j/m3/k)
DENH2O = 1000.     !density of water (kg/m3)
IF (lai+sai > 1e-6) THEN

   lsai   = lai + sai
   vegt   = lsai
   !satcap_rain = dewmx*vegt
   !satcap_snow = satcap_rain*60.0
   fveg=max(0.05,1.0-exp(-0.52*lsai))
   !the maximum bearing  capacity of the tree regardless of air temp (Imax1)
   Imax1=4.0*lsai*0.0005 *1000.0 ! in mm

   if (tair>-272.15) then 
      Lr=4.0
   else if (tair<=-272.15 .and. tair>=-270.15) then 
      Lr=1.5*(tair-273.15)+5.5
   else
      Lr=1.0
   endif

   ldew_max_snow=0.0005 *Lr *lsai * 1000.0  ! in mm !!!
   Snow=(prc_snow+prl_snow)*deltim
   if (ldew_max_snow>0.0) then 
      DeltaSnowInt=(1.0-ldew_snow/ldew_max_snow)*Snow
      if ((DeltaSnowInt+ldew_snow)>ldew_max_snow) THEN
         DeltaSnowInt=ldew_max_snow-ldew_snow
      endif
      if (DeltaSnowInt<0.0) then
         DeltaSnowInt=0.0
      endif
   else
      DeltaSnowInt=0.0
   endif 

   !* Reduce the amount of intercepted snow if windy and cold.
   !Ringyo Shikenjo Tokyo, #54, 1952.
   !Bulletin of the Govt. Forest Exp. Station,
   !Govt. Forest Exp. Station, Meguro, Tokyo, Japan.
   !FORSTX 634.9072 R475r #54.
   !Page 146, Figure 10.

   !Reduce the amount of intercepted snow if snowing, windy, and
   !cold (< -3 to -5 C).
   !Schmidt and Troendle 1992 western snow conference paper. */
   Wind= SQRT(forc_us*forc_us + forc_vs*forc_vs)
   if (tleaf<-3.0 .and. DeltaSnowInt>0.0 .and. Wind> 1.0) then
      BlownSnow=(0.2*Wind -0.2)*DeltaSnowInt
      if (BlownSnow>=DeltaSnowInt) then 
         BlownSnow = DeltaSnowInt
      endif
      DeltaSnowInt=DeltaSnowInt-BlownSnow
   endif
   
   
   ! now update snowfall and total accumulated intercepted snow amounts */
   if ((DeltaSnowInt+ldew_snow)>Imax1) THEN 
      DeltaSnowInt = 0.0
   endif

   !/* pixel depth    */
   SnowThroughFall = (Snow - DeltaSnowInt) * fveg + Snow * (1 - fveg);

   !/* Snow in canopy too thin for EB calculations; let it fall through */
   !param.SNOW_MIN_SWQ_EB_THRES = 0.0010; m*1000.
   if (Snow == 0.0  .and. ldew_snow < 0.0010*1000.0) then
      SnowThroughFall  = ldew_snow+SnowThroughFall
      DeltaSnowInt =DeltaSnowInt-ldew_snow;
   endif 

   !/* physical depth */
   ldew_snow = ldew_snow+ DeltaSnowInt;
    if (ldew_snow < 1.e-8) then
    ldew_snow = 0.0;
    endif


   !  /* Calculate amount of rain intercepted on branches and stored in  intercepted snow. */

   ! /* physical depth */
    MaxInt=exp(-4.0)*lsai !need check the unit!!  maximum interception capacity!!1
    MaxWaterInt =0.035 * (ldew_snow) + MaxInt

   Rain=(prc_rain+prl_rain)*deltim
   if (ldew_rain+Rain <=MaxWaterInt) then 
      !/* physical depth */
      ldew_rain=ldew_rain+Rain
      ! /* pixel depth */
      RainThroughFall = Rain * (1 -fveg)
   ELSE
       ! /* pixel depth */
      RainThroughFall = (ldew_rain +Rain - MaxWaterInt) * fveg +(Rain * (1.0 - fveg));
      !/* physical depth */
      ldew_rain = MaxWaterInt
   endif

   !// Liquid water in canopy too thin for EB calculations; let it fall through
   if (Rain <= 1.e-8 .and.ldew_rain < 0.0010/1000.) then
       RainThroughFall =RainThroughFall+ldew_rain;
       ldew_rain = 0.0;
   endif


   !/* at this point we have calculated the amount of snowfall intercepted and
   !/* the amount of rainfall intercepted.  These values have been
   !/* appropriately subtracted from SnowFall and RainFall to determine
   !/* SnowThroughfall and RainThroughfall.  However, we can end up with the
   !/* condition that the total intercepted rain plus intercepted snow is
   !/* greater than the maximum bearing capacity of the tree regardless of air
   !/* temp (Imax1).  The following routine will adjust ldew_rain and ldew_snow
   !/* by triggering mass release due to overloading.  Of course since ldew_rain
   !/* and ldew_snow are mixed, we need to slough them of as fixed fractions  */

   if (ldew_rain + ldew_snow > Imax1) then
      ! /*then trigger structural unloading*/
      Overload = (ldew_snow + ldew_rain) - Imax1
      IntRainFract = ldew_rain / (ldew_rain + ldew_snow)
      IntSnowFract = ldew_snow / (ldew_rain + ldew_snow)
      ldew_rain = ldew_rain - Overload * IntRainFract
      ldew_snow = ldew_snow - Overload * IntSnowFract
      RainThroughFall = RainThroughFall + (Overload * IntRainFract) * fveg
      SnowThroughFall = SnowThroughFall + (Overload * IntSnowFract) * fveg
   endif



   !here is diff from original vic 
   ! phase change and excess !
   IF (tleaf > tfrz) THEN
      ldew_smelt = MIN(ldew_snow,(tleaf-tfrz)*CICE*ldew_snow/DENICE/(HFUS))
      ldew_smelt = max(ldew_smelt,0.0)
      ldew_snow  = ldew_snow-ldew_smelt
      ldew_rain  = ldew_rain+ldew_smelt
   ELSE
      ldew_frzc  = MIN(ldew_rain,(tfrz-tleaf)*CWAT*ldew_rain/DENH2O/(HFUS))
      ldew_frzc  = max(ldew_smelt,0.0)
      ldew_snow  = ldew_snow+ldew_frzc
      ldew_rain  = ldew_rain-ldew_frzc
   ENDIF
   
   !/* Update maximum water interception storage */
   MaxInt=exp(-4.0)*lsai !need check the unit!!  maximum interception capacity!!1
   MaxWaterInt =0.035 * (ldew_snow) + MaxInt

   drip=max(0.0,ldew_rain-MaxWaterInt)
   ldew_rain=ldew_rain-drip

   ldew=ldew_rain+ldew_snow

   pg_rain=drip+RainThroughFall
   pg_snow=SnowThroughFall
   qintr_snow=-0.0
   qintr_rain=-0.0
ELSE 
   ldew = 0.
   pg_rain = prc_rain + prl_rain
   pg_snow = prc_snow + prl_snow
   qintr   = 0.
   qintr_rain = 0.
   qintr_snow = 0.
ENDIF

END SUBROUTINE LEAF_interception_vic



#ifdef PFT_CLASSIFICATION
 SUBROUTINE LEAF_interception_pftwrap (ipatch,deltim,dewmx,forc_us,forc_vs,forc_t,&
                               prc_rain,prc_snow,prl_rain,prl_snow,&
                              ldew,ldew_rain,ldew_snow,z0m,hu,pg_rain,pg_snow,qintr,qintr_rain,qintr_snow)


!=======================================================================

     USE precision
     USE PhysicalConstants, only: tfrz
     USE mod_landpft
     USE MOD_PFTimeInvars
     USE MOD_PFTimeVars
     USE MOD_1D_PFTFluxes
     USE PFT_Const

     IMPLICIT NONE

!-----------------------Arguments---------------------------------------

     INTEGER,  intent(in) :: ipatch    !patch index
     REAL(r8), intent(in) :: deltim    !seconds in a time step [second]
     REAL(r8), intent(in) :: dewmx     !maximum dew [mm]
     REAL(r8), intent(in) :: forc_us   !wind speed
     REAL(r8), intent(in) :: forc_vs   !wind speed
     REAL(r8), intent(in) :: forc_t    !air temperature
     REAL(r8), intent(in) :: z0m       ! roughness length
     REAL(r8), intent(in) :: hu        ! forcing height of U
     REAL(r8), intent(in) :: ldew_rain ! depth of water on foliage [mm]
     REAL(r8), intent(in) :: ldew_snow ! depth of water on foliage [mm]
     REAL(r8), intent(in) :: prc_rain  !convective ranfall [mm/s]
     REAL(r8), intent(in) :: prc_snow  !convective snowfall [mm/s]
     REAL(r8), intent(in) :: prl_rain  !large-scale rainfall [mm/s]
     REAL(r8), intent(in) :: prl_snow  !large-scale snowfall [mm/s]

     REAL(r8), intent(inout) :: ldew   !depth of water on foliage [mm]
     REAL(r8), intent(out) :: pg_rain  !rainfall onto ground including canopy runoff [kg/(m2 s)]
     REAL(r8), intent(out) :: pg_snow  !snowfall onto ground including canopy runoff [kg/(m2 s)]
     REAL(r8), intent(out) :: qintr    !interception [kg/(m2 s)]
     real(r8), INTENT(out) :: qintr_rain ! rainfall interception (mm h2o/s)
     real(r8), INTENT(out) :: qintr_snow ! snowfall interception (mm h2o/s)
     ! local variables
     INTEGER i, p, ps, pe
     REAL(r8) pg_rain_tmp, pg_snow_tmp

     pg_rain_tmp = 0.
     pg_snow_tmp = 0.

     ps = patch_pft_s(ipatch)      
     pe = patch_pft_e(ipatch)

     DO i = ps, pe
        p = pftclass(i)
        CALL LEAF_interception_CoLM (deltim,dewmx,forc_us,forc_vs,chil_p(p),sigf_p(i),lai_p(i),sai_p(i),forc_t,tleaf_p(i),&
           prc_rain,prc_snow,prl_rain,prl_snow,&
           ldew_p(i),ldew_p(i),ldew_p(i),z0m_p(i),hu,pg_rain,pg_snow,qintr_p(i),qintr_rain_p(i),qintr_snow_p(i))
        pg_rain_tmp = pg_rain_tmp + pg_rain*pftfrac(i)
        pg_snow_tmp = pg_snow_tmp + pg_snow*pftfrac(i)
     ENDDO 

     pg_rain = pg_rain_tmp
     pg_snow = pg_snow_tmp
     ldew  = sum( ldew_p(ps:pe) * pftfrac(ps:pe))
     qintr = sum(qintr_p(ps:pe) * pftfrac(ps:pe))
     qintr_rain = sum(qintr_rain_p(ps:pe) * pftfrac(ps:pe))
     qintr_snow = sum(qintr_snow_p(ps:pe) * pftfrac(ps:pe))
 END SUBROUTINE LEAF_interception_pftwrap
#endif
 
#ifdef PC_CLASSIFICATION
 SUBROUTINE LEAF_interception_pcwrap (ipatch,deltim,dewmx,forc_us,forc_vs,forc_t,chil,&
                               prc_rain,prc_snow,prl_rain,prl_snow,&
                              ldew,ldew_rain, ldew_snow,hu,pg_rain,pg_snow,qintr,qintr_rain,qintr_snow)


!=======================================================================

     USE precision
     USE GlobalVars
     USE PhysicalConstants, only: tfrz
     USE MOD_PCTimeInvars
     USE MOD_PCTimeVars
     USE MOD_1D_PCFluxes
     USE PFT_Const
     USE mod_landpc

     IMPLICIT NONE

!-----------------------Arguments---------------------------------------

     INTEGER,  intent(in) :: ipatch    !patch index
     REAL(r8), intent(in) :: deltim    !seconds in a time step [second]
     REAL(r8), intent(in) :: dewmx     !maximum dew [mm]
     REAL(r8), intent(in) :: forc_us   !wind speed
     REAL(r8), intent(in) :: forc_vs   !wind speed
     REAL(r8), intent(in) :: forc_t    !air temperature
     REAL(r8), intent(in) :: chil
     REAL(r8), intent(in) :: prc_rain  !convective ranfall [mm/s]
     REAL(r8), intent(in) :: prc_snow  !convective snowfall [mm/s]
     REAL(r8), intent(in) :: prl_rain  !large-scale rainfall [mm/s]
     REAL(r8), intent(in) :: prl_snow  !large-scale snowfall [mm/s]
     REAL(r8), intent(in) :: hu

     REAL(r8), intent(inout) :: ldew        !depth of water on foliage [mm]
     REAL(r8), intent(inout) :: ldew_rain   !depth of water on foliage [mm]
     REAL(r8), intent(inout) :: ldew_snow   !depth of water on foliage [mm]
     REAL(r8), intent(out) :: pg_rain  !rainfall onto ground including canopy runoff [kg/(m2 s)]
     REAL(r8), intent(out) :: pg_snow  !snowfall onto ground including canopy runoff [kg/(m2 s)]
     REAL(r8), intent(out) :: qintr    !interception [kg/(m2 s)]
     real(r8), INTENT(out) :: qintr_rain ! rainfall interception (mm h2o/s)
     real(r8), INTENT(out) :: qintr_snow ! snowfall interception (mm h2o/s)
     ! local variables
     INTEGER p, pc
     REAL(r8) pg_rain_tmp, pg_snow_tmp

     pg_rain_tmp = 0.
     pg_snow_tmp = 0.

     pc = patch2pc(ipatch)      

     DO p = 0, N_PFT-1
        CALL LEAF_interception_CoLM (deltim,dewmx,forc_us,forc_vs,&
           chil,sigf_c(p,pc),lai_c(p,pc),sai_c(p,pc),forc_t,tleaf_c(p,pc),&
           prc_rain,prc_snow,prl_rain,prl_snow,&
          ldew_c(p,pc),ldew_rain_c(p,pc),ldew_snow_c(p,pc),z0m_c(p,pc),hu,pg_rain,pg_snow,qintr_c(p,pc),qintr_rain_c(p,pc),qintr_snow_c(p,pc))

           
        pg_rain_tmp = pg_rain_tmp + pg_rain*pcfrac(p,pc)
        pg_snow_tmp = pg_snow_tmp + pg_snow*pcfrac(p,pc)
     ENDDO 

     pg_rain = pg_rain_tmp
     pg_snow = pg_snow_tmp
     ldew  = sum( ldew_c(:,pc) * pcfrac(:,pc))
     qintr = sum(qintr_c(:,pc) * pcfrac(:,pc))
     qintr = sum(qintr_rain_c(:,pc) * pcfrac(:,pc))
     qintr = sum(qintr_snow_c(:,pc) * pcfrac(:,pc))

 END SUBROUTINE LEAF_interception_pcwrap
#endif




