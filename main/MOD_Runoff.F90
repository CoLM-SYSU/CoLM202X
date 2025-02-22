#include <define.h>

MODULE MOD_Runoff

!-----------------------------------------------------------------------
   USE MOD_Precision
   IMPLICIT NONE
   SAVE

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: SurfaceRunoff_SIMTOP
   PUBLIC :: SubsurfaceRunoff_SIMTOP
   PUBLIC :: Runoff_XinAnJiang


!-----------------------------------------------------------------------

CONTAINS

   SUBROUTINE SurfaceRunoff_SIMTOP (nl_soil,wimp,porsl,psi0,hksati,&
                                    fsatmax,fsatdcf,&
                                    z_soisno,dz_soisno,zi_soisno,&
                                    eff_porosity,icefrac,zwt,gwat,&
                                    rsur,rsur_se,rsur_ie)

!=======================================================================
!  the original code was provide by Robert E. Dickinson based on
!  following clues: a water table level determination level added
!  including highland and lowland levels and fractional area of wetland
!  (water table above the surface.  Runoff is parametrized from the
!  lowlands in terms of precip incident on wet areas and a base flow,
!  where these are estimated using ideas from TOPMODEL.
!
!  Author : Yongjiu Dai, 07/29/2002, Guoyue Niu, 06/2012
!=======================================================================

   IMPLICIT NONE

!-------------------------- Dummy Arguments ----------------------------

   integer, intent(in) :: nl_soil   ! number of soil layers
   real(r8), intent(in) :: &
        ! wtfact,                 &! (updated to gridded 'fsatmax' data) fraction of model area with high water table
        wimp,                     &! water impermeable if porosity less than wimp
        porsl(1:nl_soil),         &! saturated volumetric soil water content(porosity)
        psi0(1:nl_soil),          &! saturated soil suction (mm) (NEGATIVE)
        hksati(1:nl_soil),        &! hydraulic conductivity at saturation (mm h2o/s)
        fsatmax,                  &! maximum fraction of saturation area [-]
        fsatdcf,                  &! decay factor in calculation of fraction of saturation area [1/m]
        z_soisno(1:nl_soil),      &! layer depth (m)
        dz_soisno(1:nl_soil),     &! layer thickness (m)
        zi_soisno(0:nl_soil),     &! interface level below a "z" level (m)
        eff_porosity(1:nl_soil),  &! effective porosity = porosity - vol_ice
        icefrac(1:nl_soil),       &! ice fraction (-)
        gwat,                     &! net water input from top
        zwt                        ! the depth from ground (soil) surface to water table [m]

   real(r8), intent(out) :: rsur   ! surface runoff (mm h2o/s)
   real(r8), intent(out), optional :: rsur_se! saturation excess surface runoff (mm h2o/s)
   real(r8), intent(out), optional :: rsur_ie! infiltration excess surface runoff (mm h2o/s)

!-------------------------- Local Variables ----------------------------

   real(r8) qinmax       ! maximum infiltration capability
   real(r8) fsat         ! fractional area with water table at surface

   ! updated to gridded 'fsatdcf' (by Shupeng Zhang)
   ! real(r8), parameter :: fff = 0.5   ! runoff decay factor (m-1)

!-----------------------------------------------------------------------

!  fraction of saturated area (updated to gridded 'fsatmax' and 'fsatdcf')
      !fsat = wtfact*min(1.0,exp(-0.5*fff*zwt))
      fsat = fsatmax*min(1.0,exp(-2.0*fsatdcf*zwt))

! Maximum infiltration capacity
      qinmax = minval(10.**(-6.0*icefrac(1:min(3,nl_soil)))*hksati(1:min(3,nl_soil)))
      IF(eff_porosity(1)<wimp) qinmax = 0.

! Surface runoff
      rsur = fsat*max(0.0,gwat) + (1.-fsat)*max(0.,gwat-qinmax)

      IF (present(rsur_se)) THEN
         rsur_se = fsat*max(0.0,gwat)
      ENDIF

      IF (present(rsur_ie)) THEN
         rsur_ie = (1.-fsat)*max(0.,gwat-qinmax)
      ENDIF

   END SUBROUTINE SurfaceRunoff_SIMTOP

! -------------------------------------------------------------------------
   SUBROUTINE SubsurfaceRunoff_SIMTOP (nl_soil, icefrac, dz_soisno, zi_soisno, zwt, rsubst)

   IMPLICIT NONE

!-------------------------- Dummy Arguments ----------------------------
   integer,  intent(in) :: nl_soil                 !
   real(r8), intent(in) :: icefrac(1:nl_soil)      ! ice fraction (-)

   real(r8), intent(in) :: dz_soisno  (1:nl_soil)  ! layer depth (m)
   real(r8), intent(in) :: zi_soisno  (0:nl_soil)  ! interface level below a "z" level (m)

   real(r8), intent(in)  :: zwt    ! the depth from ground (soil) surface to water table [m]
   real(r8), intent(out) :: rsubst ! subsurface runoff (positive = out of soil column) (mm H2O /s)

!-------------------------- Local Variables ----------------------------

   integer  :: j                ! indices
   integer  :: jwt              ! index of the soil layer right above the water table (-)
   real(r8) :: dzmm(1:nl_soil)  ! layer thickness (mm)

   real(r8) :: dzsum
   real(r8) :: icefracsum
   real(r8) :: fracice_rsub
   real(r8) :: imped
!-----------------------------------------------------------------------

      DO j = 1,nl_soil
         dzmm(j) = dz_soisno(j)*1000.
      ENDDO

      jwt = nl_soil
      ! allow jwt to equal zero when zwt is in top layer
      DO j = 1, nl_soil
         IF(zwt <= zi_soisno(j)) THEN
            jwt = j-1
            EXIT
         ENDIF
      ENDDO

      !-- Topographic runoff  --
      dzsum = 0.
      icefracsum = 0.
      DO j = max(jwt,1), nl_soil
         dzsum = dzsum + dzmm(j)
         icefracsum = icefracsum + icefrac(j) * dzmm(j)
      ENDDO
      ! add ice impedance factor to baseflow
      fracice_rsub = max(0.,exp(-3.*(1.-(icefracsum/dzsum)))-exp(-3.))/(1.0-exp(-3.))
      imped = max(0.,1.-fracice_rsub)
      rsubst = imped * 5.5e-3 * exp(-2.5*zwt)

   END SUBROUTINE SubsurfaceRunoff_SIMTOP

! -------------------------------------------------------------------------
   SUBROUTINE Runoff_XinAnJiang ( &
         nl_soil, dz_soisno, eff_porosity, vol_liq, topostd, gwat, deltim, &
         rsur, rsubst)

   USE MOD_Precision
   IMPLICIT NONE

   integer,  intent(in) :: nl_soil ! number of soil layers

   real(r8), intent(in) :: &
        dz_soisno   (1:nl_soil),  &! layer thickness (m)
        eff_porosity(1:nl_soil),  &! effective porosity = porosity - vol_ice
        vol_liq     (1:nl_soil),  &! partial volume of liquid water in layer
        topostd,                  &! standard deviation of elevation (m)
        gwat,                     &! net water input from top
        deltim                     ! time step (s)

   real(r8), intent(out) :: rsur   ! surface runoff (mm h2o/s)
   real(r8), intent(out) :: rsubst ! subsurface runoff (mm h2o/s)

   ! Local Variables
   real(r8) :: btopo, watin, w_int, wsat_int, wtmp, infil
   real(r8), parameter :: sigmin = 100.
   real(r8), parameter :: sigmax = 1000.

      watin = gwat * deltim / 1000.

      IF (watin <= 0.) THEN

         rsur   = 0.
         rsubst = 0.

      ELSE

         btopo = (topostd - sigmin) / (topostd + sigmax)
         btopo = min(max(btopo, 0.01), 0.5)

         w_int    = sum(vol_liq     (1:6) * dz_soisno(1:6))
         wsat_int = sum(eff_porosity(1:6) * dz_soisno(1:6))

         wtmp  = (1-w_int/wsat_int)**(1/(btopo+1)) - watin/((btopo+1)*wsat_int)
         infil = wsat_int - w_int - wsat_int * (max(0., wtmp))**(btopo+1)

         infil = min(infil, watin)

         rsur   = (watin - infil) * 1000. / deltim
         rsubst = 0.

      ENDIF

   END SUBROUTINE Runoff_XinAnJiang


   ! -------------------------------------------------------------------------
   SUBROUTINE Runoff_SimpleVIC ( &
      nl_soil, dz_soisno, eff_porosity, vol_liq, BVIC, gwat, deltim, &
      rsur, rsubst)

   USE MOD_Precision
   IMPLICIT NONE

   integer,  intent(in) :: nl_soil ! number of soil layers

   real(r8), intent(in) :: &
      dz_soisno   (1:nl_soil),  &  ! layer thickness (m)
      eff_porosity(1:nl_soil),  &  ! effective porosity = porosity - vol_ice
      vol_liq     (1:nl_soil),  &  ! partial volume of liquid water in layer
      BVIC,                     &  ! VIC infiltration parameter
      gwat,                     &  ! net water input from top
      deltim                       ! time step (s)

   real(r8), intent(out) :: rsur   ! surface runoff (mm h2o/s)
   real(r8), intent(out) :: rsubst ! subsurface runoff (mm h2o/s)

   ! Local Variables
   real(r8) :: btopo, watin, w_int, wsat_int, wtmp, infil
   real(r8) :: InfilExpFac, WaterDepthMax, WaterDepthInit, RunoffSurface, InfilVarTmp
   real(r8) :: SoilSaturateFrac

      watin = gwat * deltim / 1000. ! convert mm/s to m

      IF (watin <= 0.) THEN
         rsur   = 0.
         rsubst = 0.
      ELSE
         w_int    = sum(vol_liq     (1:6) * dz_soisno(1:6))
         wsat_int = sum(eff_porosity(1:6) * dz_soisno(1:6))
         ! fractional saturated area from soil moisture
         InfilExpFac            = BVIC / ( 1.0 + BVIC )
         SoilSaturateFrac = 1.0 - (max(0.0, (1.0-(w_int/wsat_int))))**InfilExpFac
         SoilSaturateFrac = max(0.0, SoilSaturateFrac)
         SoilSaturateFrac = min(1.0, SoilSaturateFrac)

         ! Infiltration for the previous time-step soil moisture based on SoilSaturateFrac
         WaterDepthMax  = (1.0 + BVIC) * wsat_int
         WaterDepthInit = WaterDepthMax * (1.0 - (1.0 - SoilSaturateFrac)**(1.0/BVIC))

         ! Solve for surface runoff
         if ( WaterDepthMax <= 0.0 ) then
            RunoffSurface = watin
         ELSEIF   ( (WaterDepthInit + watin) > WaterDepthMax ) then
            !RunoffSurface = (WaterDepthInit + w_int) - WaterDepthMax
            RunoffSurface = watin - wsat_int + w_int
         ELSE
            InfilVarTmp  = 1.0 - ((WaterDepthInit +watin ) / WaterDepthMax)
            RunoffSurface =watin - wsat_int + w_int + wsat_int * (InfilVarTmp**(1.0+BVIC))
         ENDIF

         IF ( RunoffSurface < 0.0 ) RunoffSurface = 0.0
         IF ( RunoffSurface > watin) RunoffSurface = watin

         infil = watin - RunoffSurface
         rsur= RunoffSurface * 1000. / deltim
         rsubst = 0.
      ENDIF

   END SUBROUTINE Runoff_SimpleVIC

END MODULE MOD_Runoff
! ---------- EOP ------------
