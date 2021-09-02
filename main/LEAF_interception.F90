#include <define.h>

 SUBROUTINE LEAF_interception (deltim,dewmx,chil,sigf,lai,sai,tleaf,&
                               prc_rain,prc_snow,prl_rain,prl_snow,&
                               ldew,pg_rain,pg_snow,qintr)

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
  REAL(r8), intent(in) :: chil      !leaf angle distribution factor
  REAL(r8), intent(in) :: prc_rain  !convective ranfall [mm/s]
  REAL(r8), intent(in) :: prc_snow  !convective snowfall [mm/s]
  REAL(r8), intent(in) :: prl_rain  !large-scale rainfall [mm/s]
  REAL(r8), intent(in) :: prl_snow  !large-scale snowfall [mm/s]
  REAL(r8), intent(in) :: sigf      !fraction of veg cover, excluding snow-covered veg [-]
  REAL(r8), intent(in) :: lai       !leaf area index [-]
  REAL(r8), intent(in) :: sai       !stem area index [-]
  REAL(r8), intent(in) :: tleaf     !sunlit canopy leaf temperature [K]

  REAL(r8), intent(inout) :: ldew   !depth of water on foliage [mm]
  REAL(r8), intent(out) :: pg_rain  !rainfall onto ground including canopy runoff [kg/(m2 s)]
  REAL(r8), intent(out) :: pg_snow  !snowfall onto ground including canopy runoff [kg/(m2 s)]
  REAL(r8), intent(out) :: qintr    !interception [kg/(m2 s)]

!-----------------------Local Variabes---------------------------------

  REAL(r8) :: satcap    !maximum allowed water on canopy [mm]
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

  REAL(r8) :: ap, cp, bp, aa, bb, exrain, arg, alpha, w
  REAL(r8) :: thru_rain, thru_snow
  REAL(r8) :: xsc_rain, xsc_snow
  REAL pcoefs (2,2)

!-----------------------End Variable List-------------------------------

 IF (lai+sai >= 1e-6) THEN

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
       tex_rain = (prc_rain+prl_rain)*deltim * fpi * (ap/bp*(1.-exp(-bp*xs))+cp*xs) &
                - (satcap-ldew) * xs
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

 ENDIF

 END SUBROUTINE LEAF_interception

#ifdef PFT_CLASSIFICATION
 SUBROUTINE LEAF_interception_pftwrap (ipatch,deltim,dewmx,&
                               prc_rain,prc_snow,prl_rain,prl_snow,&
                               ldew,pg_rain,pg_snow,qintr)

!=======================================================================

     USE precision
     USE PhysicalConstants, only: tfrz
     USE MOD_PFTimeInvars
     USE MOD_PFTimeVars
     USE MOD_1D_PFTFluxes
     USE PFT_Const

     IMPLICIT NONE

!-----------------------Arguments---------------------------------------

     INTEGER,  intent(in) :: ipatch    !patch index
     REAL(r8), intent(in) :: deltim    !seconds in a time step [second]
     REAL(r8), intent(in) :: dewmx     !maximum dew [mm]
     REAL(r8), intent(in) :: prc_rain  !convective ranfall [mm/s]
     REAL(r8), intent(in) :: prc_snow  !convective snowfall [mm/s]
     REAL(r8), intent(in) :: prl_rain  !large-scale rainfall [mm/s]
     REAL(r8), intent(in) :: prl_snow  !large-scale snowfall [mm/s]

     REAL(r8), intent(inout) :: ldew   !depth of water on foliage [mm]
     REAL(r8), intent(out) :: pg_rain  !rainfall onto ground including canopy runoff [kg/(m2 s)]
     REAL(r8), intent(out) :: pg_snow  !snowfall onto ground including canopy runoff [kg/(m2 s)]
     REAL(r8), intent(out) :: qintr    !interception [kg/(m2 s)]

     ! local variables
     INTEGER i, p, ps, pe
     REAL(r8) pg_rain_tmp, pg_snow_tmp

     pg_rain_tmp = 0.
     pg_snow_tmp = 0.

     ps = patch_pft_s(ipatch)      
     pe = patch_pft_e(ipatch)

     DO i = ps, pe
        p = pftclass(i)
        CALL LEAF_interception (deltim,dewmx,chil_p(p),sigf_p(i),lai_p(i),sai_p(i),tleaf_p(i),&
           prc_rain,prc_snow,prl_rain,prl_snow,&
           ldew_p(i),pg_rain,pg_snow,qintr_p(i))
        pg_rain_tmp = pg_rain_tmp + pg_rain*pftfrac(i)
        pg_snow_tmp = pg_snow_tmp + pg_snow*pftfrac(i)
     ENDDO 

     pg_rain = pg_rain_tmp
     pg_snow = pg_snow_tmp
     ldew  = sum( ldew_p(ps:pe) * pftfrac(ps:pe))
     qintr = sum(qintr_p(ps:pe) * pftfrac(ps:pe))

 END SUBROUTINE LEAF_interception_pftwrap
#endif
 
#ifdef PC_CLASSIFICATION
 SUBROUTINE LEAF_interception_pcwrap (ipatch,deltim,dewmx,&
                               prc_rain,prc_snow,prl_rain,prl_snow,&
                               ldew,pg_rain,pg_snow,qintr)

!=======================================================================

     USE precision
     USE GlobalVars
     USE PhysicalConstants, only: tfrz
     USE MOD_PCTimeInvars
     USE MOD_PCTimeVars
     USE MOD_1D_PCFluxes
     USE PFT_Const

     IMPLICIT NONE

!-----------------------Arguments---------------------------------------

     INTEGER,  intent(in) :: ipatch    !patch index
     REAL(r8), intent(in) :: deltim    !seconds in a time step [second]
     REAL(r8), intent(in) :: dewmx     !maximum dew [mm]
     REAL(r8), intent(in) :: prc_rain  !convective ranfall [mm/s]
     REAL(r8), intent(in) :: prc_snow  !convective snowfall [mm/s]
     REAL(r8), intent(in) :: prl_rain  !large-scale rainfall [mm/s]
     REAL(r8), intent(in) :: prl_snow  !large-scale snowfall [mm/s]

     REAL(r8), intent(inout) :: ldew   !depth of water on foliage [mm]
     REAL(r8), intent(out) :: pg_rain  !rainfall onto ground including canopy runoff [kg/(m2 s)]
     REAL(r8), intent(out) :: pg_snow  !snowfall onto ground including canopy runoff [kg/(m2 s)]
     REAL(r8), intent(out) :: qintr    !interception [kg/(m2 s)]

     ! local variables
     INTEGER p, pc
     REAL(r8) pg_rain_tmp, pg_snow_tmp

     pg_rain_tmp = 0.
     pg_snow_tmp = 0.

     pc = patch2pc(ipatch)      

     DO p = 0, N_PFT-1
        CALL LEAF_interception (deltim,dewmx,&
           chil_p(p),sigf_c(p,pc),lai_c(p,pc),sai_c(p,pc),tleaf_c(p,pc),&
           prc_rain,prc_snow,prl_rain,prl_snow,&
           ldew_c(p,pc),pg_rain,pg_snow,qintr_c(p,pc))
        pg_rain_tmp = pg_rain_tmp + pg_rain*pcfrac(p,pc)
        pg_snow_tmp = pg_snow_tmp + pg_snow*pcfrac(p,pc)
     ENDDO 

     pg_rain = pg_rain_tmp
     pg_snow = pg_snow_tmp
     ldew  = sum( ldew_c(:,pc) * pcfrac(:,pc))
     qintr = sum(qintr_c(:,pc) * pcfrac(:,pc))

 END SUBROUTINE LEAF_interception_pcwrap
#endif
