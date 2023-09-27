#include <define.h>

MODULE MOD_SoilSnowHydrology

!-----------------------------------------------------------------------
  use MOD_Precision
  use MOD_Namelist, only: DEF_USE_PLANTHYDRAULICS, DEF_USE_SNICAR, &
                          DEF_URBAN_RUN, DEF_USE_IRRIGATION
#if(defined CaMa_Flood)
   USE YOS_CMF_INPUT,      ONLY: LWINFILT
#endif
#ifdef CROP
   use MOD_LandPFT, only: patch_pft_s, patch_pft_e
   use MOD_Irrigation, only: CalIrrigationApplicationFluxes
#endif
  IMPLICIT NONE
  SAVE

! PUBLIC MEMBER FUNCTIONS:
  public :: WATER
  PUBLIC :: WATER_VSF
  public :: snowwater
  public :: soilwater
  public :: snowwater_snicar


! PRIVATE MEMBER FUNCTIONS:
  private :: surfacerunoff
  private :: subsurfacerunoff


!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------



  subroutine WATER (ipatch,patchtype  ,lb          ,nl_soil ,deltim,&
             z_soisno    ,dz_soisno   ,zi_soisno                   ,&
             bsw         ,porsl       ,psi0        ,hksati  ,rootr ,&
             t_soisno    ,wliq_soisno ,wice_soisno ,smp     ,hk    ,pg_rain ,sm    ,&
             etr         ,qseva       ,qsdew       ,qsubl   ,qfros ,&
             rsur        ,rnof        ,qinfl       ,wtfact  ,pondmx,&
             ssi         ,wimp        ,smpmin      ,zwt     ,wa    ,&
             qcharge     ,errw_rsub     &
#if(defined CaMa_Flood)
            ,flddepth,fldfrc,qinfl_fld  &
#endif
! SNICAR model variables
             ,forc_aer   ,&
             mss_bcpho   ,mss_bcphi   ,mss_ocpho   ,mss_ocphi   ,&
             mss_dst1    ,mss_dst2    ,mss_dst3    ,mss_dst4     &
! END SNICAR model variables
            )
!=======================================================================
! this is the main subroutine to execute the calculation of
! hydrological processes
!
! Original author : Yongjiu Dai, /09/1999/, /08/2002/, /04/2014/
!
! FLOW DIAGRAM FOR WATER.F90
!
! WATER ===> snowwater
!            surfacerunoff
!            soilwater
!            subsurfacerunoff
!
!=======================================================================

  use MOD_Precision
  use MOD_Const_Physical, only : denice, denh2o, tfrz

  implicit none

!-----------------------Argument---------- ------------------------------
  integer, INTENT(in) :: &
        ipatch           ,& ! patch index
        patchtype           ! land patch type (0=soil, 1=urban or built-up, 2=wetland,
                            ! 3=land ice, 4=land water bodies, 99=ocean

  integer, INTENT(in) :: &
        lb               , &! lower bound of array
        nl_soil             ! upper bound of array

  real(r8), INTENT(in) :: &
        deltim           , &! time step (s)
        wtfact           , &! fraction of model area with high water table
        pondmx           , &! ponding depth (mm)
        ssi              , &! irreducible water saturation of snow
        wimp             , &! water impremeable if porosity less than wimp
        smpmin           , &! restriction for min of soil poten. (mm)
        z_soisno (lb:nl_soil)   , &! layer depth (m)
        dz_soisno(lb:nl_soil)   , &! layer thickness (m)
        zi_soisno(lb-1:nl_soil) , &! interface level below a "z" level (m)

        bsw(1:nl_soil)   , &! Clapp-Hornberger "B"
        porsl(1:nl_soil) , &! saturated volumetric soil water content(porosity)
        psi0(1:nl_soil)  , &! saturated soil suction (mm) (NEGATIVE)
        hksati(1:nl_soil), &! hydraulic conductivity at saturation (mm h2o/s)
        rootr(1:nl_soil) , &! root resistance of a layer, all layers add to 1.0

        t_soisno(lb:nl_soil), &! soil/snow skin temperature (K)
        pg_rain          , &! rainfall after removal of interception (mm h2o/s)
        sm               , &! snow melt (mm h2o/s)
        etr              , &! actual transpiration (mm h2o/s)
        qseva            , &! ground surface evaporation rate (mm h2o/s)
        qsdew            , &! ground surface dew formation (mm h2o /s) [+]
        qsubl            , &! sublimation rate from snow pack (mm h2o /s) [+]
        qfros               ! surface dew added to snow pack (mm h2o /s) [+]
#if(defined CaMa_Flood)
         real(r8), INTENT(inout) :: flddepth  ! inundation water depth [mm]
         real(r8), INTENT(in)    :: fldfrc    ! inundation water depth   [0-1]
         real(r8), INTENT(out)   :: qinfl_fld ! grid averaged inundation water input from top (mm/s)
#endif
  real(r8), INTENT(inout) :: &
        wice_soisno(lb:nl_soil) , &! ice lens (kg/m2)
        wliq_soisno(lb:nl_soil) , &! liquid water (kg/m2)
        smp(1:nl_soil)   , &! soil matrix potential [mm]
        hk (1:nl_soil)   , &! hydraulic conductivity [mm h2o/m]
        zwt              , &! the depth from ground (soil) surface to water table [m]
        wa                  ! water storage in aquifer [mm]

  real(r8), INTENT(out) :: &
        rsur             , &! surface runoff (mm h2o/s)
        rnof             , &! total runoff (mm h2o/s)
        qinfl            , &! infiltration rate (mm h2o/s)
        qcharge          , &! groundwater recharge (positive to aquifer) [mm/s]
        errw_rsub

! SNICAR model variables
! Aerosol Fluxes (Jan. 07, 2023)
  real(r8), intent(in) :: forc_aer ( 14 )  ! aerosol deposition from atmosphere model (grd,aer) [kg m-1 s-1]

  real(r8), INTENT(inout) :: &
        mss_bcpho (lb:0), &! mass of hydrophobic BC in snow  (col,lyr) [kg]
        mss_bcphi (lb:0), &! mass of hydrophillic BC in snow (col,lyr) [kg]
        mss_ocpho (lb:0), &! mass of hydrophobic OC in snow  (col,lyr) [kg]
        mss_ocphi (lb:0), &! mass of hydrophillic OC in snow (col,lyr) [kg]
        mss_dst1  (lb:0), &! mass of dust species 1 in snow  (col,lyr) [kg]
        mss_dst2  (lb:0), &! mass of dust species 2 in snow  (col,lyr) [kg]
        mss_dst3  (lb:0), &! mass of dust species 3 in snow  (col,lyr) [kg]
        mss_dst4  (lb:0)   ! mass of dust species 4 in snow  (col,lyr) [kg]
! Aerosol Fluxes (Jan. 07, 2023)
! END SNICAR model variables

!-----------------------Local Variables------------------------------
!
  integer j                 ! loop counter

  real(r8) :: &
  eff_porosity(1:nl_soil), &! effective porosity = porosity - vol_ice
       dwat(1:nl_soil)   , &! change in soil water
       gwat              , &! net water input from top (mm/s)
       rsubst            , &! subsurface runoff (mm h2o/s)
       vol_liq(1:nl_soil), &! partitial volume of liquid water in layer
       vol_ice(1:nl_soil), &! partitial volume of ice lens in layer
       icefrac(1:nl_soil), &! ice fraction (-)
       zmm (1:nl_soil)   , &! layer depth (mm)
       dzmm(1:nl_soil)   , &! layer thickness (mm)
       zimm(0:nl_soil)      ! interface level below a "z" level (mm)

  real(r8) :: err_solver, w_sum
#if(defined CaMa_Flood)
  real(r8) ::gfld ,rsur_fld, qinfl_fld_subgrid ! inundation water input from top (mm/s)
#endif

#ifdef CROP
   integer  :: ps, pe
   integer  :: irrig_flag  ! 1 if sprinker, 2 if others
   real(r8) :: qflx_irrig_drip
   real(r8) :: qflx_irrig_sprinkler
   real(r8) :: qflx_irrig_flood
   real(r8) :: qflx_irrig_paddy
#endif

!=======================================================================
! [1] update the liquid water within snow layer and the water onto soil
!=======================================================================

      if (lb>=1)then
         gwat = pg_rain + sm - qseva
      else
         IF (.not. DEF_USE_SNICAR .or. (patchtype==1 .and. DEF_URBAN_RUN)) THEN
            call snowwater (lb,deltim,ssi,wimp,&
                         pg_rain,qseva,qsdew,qsubl,qfros,&
                         dz_soisno(lb:0),wice_soisno(lb:0),wliq_soisno(lb:0),gwat)
         ELSE
            call snowwater_snicar (lb,deltim,ssi,wimp,&
                         pg_rain,qseva,qsdew,qsubl,qfros,&
                         dz_soisno(lb:0),wice_soisno(lb:0),wliq_soisno(lb:0),gwat,&
                         forc_aer,&
                         mss_bcpho(lb:0), mss_bcphi(lb:0), mss_ocpho(lb:0), mss_ocphi(lb:0),&
                         mss_dst1(lb:0), mss_dst2(lb:0), mss_dst3(lb:0), mss_dst4(lb:0) )
         ENDIF
      endif
#ifdef CROP
      if(DEF_USE_IRRIGATION)then
         if(patchtype==0)then
            ps = patch_pft_s(ipatch)
            pe = patch_pft_e(ipatch)
            call CalIrrigationApplicationFluxes(ipatch,ps,pe,deltim,qflx_irrig_drip,qflx_irrig_sprinkler,qflx_irrig_flood,qflx_irrig_paddy,irrig_flag=2)
            gwat = gwat + qflx_irrig_drip + qflx_irrig_flood + qflx_irrig_paddy
         end if
      end if
#endif
!=======================================================================
! [2] surface runoff and infiltration
!=======================================================================

  if(patchtype<=1)then   ! soil ground only

      ! For water balance check, the sum of water in soil column before the calcultion
      w_sum = sum(wliq_soisno(1:)) + sum(wice_soisno(1:)) + wa

      ! porosity of soil, partitial volume of ice and liquid
      do j = 1, nl_soil
         vol_ice(j) = min(porsl(j), wice_soisno(j)/(dz_soisno(j)*denice))
         eff_porosity(j) =  max(0.01, porsl(j)-vol_ice(j))
         vol_liq(j) = min(eff_porosity(j), wliq_soisno(j)/(dz_soisno(j)*denh2o))
         if(porsl(j) < 1.e-6)then
            icefrac(j) = 0.
         else
            icefrac(j) = min(1.,vol_ice(j)/porsl(j))
         endif
      enddo

      ! surface runoff including water table and surface staturated area

      rsur = 0.
      if (gwat > 0.) then
      call surfacerunoff (nl_soil,wtfact,wimp,porsl,psi0,hksati,&
                          z_soisno(1:),dz_soisno(1:),zi_soisno(0:),&
                          eff_porosity,icefrac,zwt,gwat,rsur)
      else
           rsur = 0.
      endif

      ! infiltration into surface soil layer
      qinfl = gwat - rsur
#if(defined CaMa_Flood)
   IF (LWINFILT) then
         !  re-infiltration [mm/s] calculation.
         ! if surface runoff is ocurred (rsur != 0.), flood depth <1.e-6  and flood frction <0.05,
         ! the re-infiltration will not be calculated.
      IF ((flddepth .GT. 1.e-6).and.(fldfrc .GT. 0.05) .and. (patchtype == 0) ) then
!         write(6,*) 'flddepth=',flddepth,'fldfrc=',fldfrc
         gfld=flddepth/deltim ! [mm/s]
         ! surface runoff from inundation, this should not be added to the surface runoff from soil
         ! otherwise, the surface runoff will be double counted.
         ! only the re-infiltration is added to water balance calculation.
         CALL surfacerunoff (nl_soil,1.0,wimp,porsl,psi0,hksati,&
                    z_soisno(1:),dz_soisno(1:),zi_soisno(0:),&
                    eff_porosity,icefrac,zwt,gfld,rsur_fld)
         ! infiltration into surface soil layer
         qinfl_fld_subgrid = gfld - rsur_fld !assume the re-infiltration is occured in whole patch area.
!         write(6,*) 'gfld=',gfld,'   qinfl_fld_subgrid=',qinfl_fld_subgrid
      ELSE
         qinfl_fld_subgrid=0.0d0
         gfld=0.0d0
         rsur_fld=0.0d0

      ENDIF
         qinfl_fld=qinfl_fld_subgrid*fldfrc ! [mm/s] re-infiltration in grid.
         qinfl=qinfl_fld+qinfl ! [mm/s] total infiltration in grid.
         flddepth=flddepth-deltim*qinfl_fld_subgrid ! renew flood depth [mm], the flood depth is reduced by re-infiltration but only in inundation area.
   ENDIF
#endif

!=======================================================================
! [3] determine the change of soil water
!=======================================================================

      ! convert length units from m to mm
      zmm(1:) = z_soisno(1:)*1000.
      dzmm(1:) = dz_soisno(1:)*1000.
      zimm(0:) = zi_soisno(0:)*1000.

      call soilwater(patchtype,nl_soil,deltim,wimp,smpmin,&
                     qinfl,etr,z_soisno(1:),dz_soisno(1:),zi_soisno(0:),&
                     t_soisno(1:),vol_liq,vol_ice,smp,hk,icefrac,eff_porosity,&
                     porsl,hksati,bsw,psi0,rootr,&
                     zwt,dwat,qcharge)

      ! update the mass of liquid water
      do j= 1, nl_soil
         wliq_soisno(j) = wliq_soisno(j)+dwat(j)*dzmm(j)
      enddo

!=======================================================================
! [4] subsurface runoff and the corrections
!=======================================================================

      call subsurfacerunoff (nl_soil,deltim,pondmx,&
                             eff_porosity,icefrac,dz_soisno(1:),zi_soisno(0:),&
                             wice_soisno(1:),wliq_soisno(1:),&
                             porsl,psi0,bsw,zwt,wa,&
                             qcharge,rsubst,errw_rsub)

      ! total runoff (mm/s)
      rnof = rsubst + rsur


      ! Renew the ice and liquid mass due to condensation
      if(lb >= 1)then
         ! make consistent with how evap_grnd removed in infiltration
         wliq_soisno(1) = max(0., wliq_soisno(1) + qsdew * deltim)
         wice_soisno(1) = max(0., wice_soisno(1) + (qfros-qsubl) * deltim)
      end if

      err_solver = (sum(wliq_soisno(1:))+sum(wice_soisno(1:))+wa) - w_sum &
                 - (gwat-etr-rnof-errw_rsub)*deltim

      if(lb >= 1)then
         err_solver = err_solver-(qsdew+qfros-qsubl)*deltim
      endif
#if(defined CaMa_Flood)
      IF (LWINFILT) THEN
         err_solver = err_solver-(gfld-rsur_fld)*fldfrc*deltim
      ENDIF
#endif


#if(defined CoLMDEBUG)
     if(abs(err_solver) > 1.e-3)then
        write(6,*) 'Warning: water balance violation after all soilwater calculation', err_solver
     endif
#endif


!=======================================================================
! [6] assumed hydrological scheme for the wetland and glacier
!=======================================================================

  else
      if(patchtype==2)then        ! WETLAND
         ! 09/20/2019, by Chaoqun Li: a potential bug below
         ! surface runoff could > total runoff
         ! original CoLM: rusr=0., qinfl=gwat, rsubst=0., rnof=0.
         ! i.e., all water to be infiltration
         qinfl = 0.
         rsur = max(0.,gwat)
         rsubst = 0.
         rnof = 0.
         do j = 1, nl_soil
            if(t_soisno(j)>tfrz)then
               wice_soisno(j) = 0.0
               wliq_soisno(j) = porsl(j)*dz_soisno(j)*1000.
            endif
         enddo
      endif
      if(patchtype==3)then        ! LAND ICE
         rsur = max(0.0,gwat)
         qinfl = 0.
         rsubst = 0.
         rnof = rsur
         wice_soisno(1:nl_soil) = dz_soisno(1:nl_soil)*1000.
         wliq_soisno(1:nl_soil) = 0.0
      endif

      wa = 4800.
      zwt = 0.
      qcharge = 0.
      errw_rsub = 0.

  endif

  end subroutine WATER

!-----------------------------------------------------------------------
  subroutine WATER_VSF (ipatch,  patchtype,lb      ,nl_soil ,deltim ,&
             z_soisno    ,dz_soisno   ,zi_soisno                    ,&
#ifdef Campbell_SOIL_MODEL
             bsw         ,                                           &
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
             theta_r     ,alpha_vgm   ,n_vgm       ,L_vgm   ,        &
             sc_vgm      ,fc_vgm      ,                              &
#endif
             porsl       ,psi0        ,hksati      ,rootr   ,        &
             t_soisno    ,wliq_soisno ,wice_soisno ,smp     ,hk     ,&
             pg_rain     ,sm          ,                              &
             etr         ,qseva       ,qsdew       ,qsubl   ,qfros  ,&
             rsur        ,rnof        ,qinfl       ,wtfact  ,ssi    ,&
             pondmx      ,                                           &
             wimp        ,zwt         ,wdsrf       ,wa      ,wetwat ,&
             qcharge     ,errw_rsub                                  &
#if(defined CaMa_Flood)
             ,flddepth,fldfrc,qinfl_fld&
#endif
! SNICAR model variables
             ,forc_aer   ,&
             mss_bcpho   ,mss_bcphi   ,mss_ocpho   ,mss_ocphi ,&
             mss_dst1    ,mss_dst2    ,mss_dst3    ,mss_dst4   &
! END SNICAR model variables
             )

!===================================================================================
! this is the main subroutine to execute the calculation of soil water processes
!
! Original author : Yongjiu Dai, /09/1999/, /08/2002/, /04/2014/
!
! Modified by Shupeng Zhang /07/2023/ to use Variably Saturated Flow algorithm
!    Reference :
!    Dai, Y., Zhang, S., Yuan, H., & Wei, N. (2019).
!    Modeling Variably Saturated Flow in Stratified Soils
!      With Explicit Tracking of Wetting Front and Water Table Locations.
!    Water Resources Research. doi:10.1029/2019wr025368
!
! FLOW DIAGRAM FOR WATER_VSF.F90
!
! WATER ===> snowwater
!            surfacerunoff     [caculated by lateral flow when defined Lateral_Flow]
!            subsurfacerunoff  [caculated by lateral flow when defined Lateral_Flow]
!            soilwater         [Variably Saturated Flow algorithm]
!
!===================================================================================

  use MOD_Precision
  USE MOD_Hydro_SoilWater
  USE MOD_Vars_TimeInvariants, only : wetwatmax
  use MOD_Const_Physical, only : denice, denh2o, tfrz
#ifdef DataAssimilation
  USE MOD_DA_GRACE, only : fslp_patch
#endif

  implicit none

!-----------------------Argument---------- ------------------------------
  integer, INTENT(in) :: &
        ipatch           ,& ! patch index
        patchtype           ! land patch type (0=soil, 1=urban or built-up, 2=wetland,
                            ! 3=land ice, 4=land water bodies, 99=ocean

  integer, INTENT(in) :: &
        lb               , &! lower bound of array
        nl_soil            ! upper bound of array

  real(r8), INTENT(in) :: &
        deltim           , &! time step (s)
        wtfact           , &! fraction of model area with high water table
        ssi              , &! irreducible water saturation of snow
        pondmx           , &! ponding depth (mm)
        wimp             , &! water impremeable if porosity less than wimp
        z_soisno (lb:nl_soil)   , &! layer depth (m)
        dz_soisno(lb:nl_soil)   , &! layer thickness (m)
        zi_soisno(lb-1:nl_soil) , &! interface level below a "z" level (m)
#ifdef Campbell_SOIL_MODEL
        bsw      (1:nl_soil), &! clapp and hornbereger "b" parameter [-]
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
        theta_r  (1:nl_soil), & ! residual moisture content [-]
        alpha_vgm(1:nl_soil), & ! a parameter corresponding approximately to the inverse of the air-entry value
        n_vgm    (1:nl_soil), & ! a shape parameter [dimensionless]
        L_vgm    (1:nl_soil), & ! pore-connectivity parameter [dimensionless]
        sc_vgm   (1:nl_soil), & ! saturation at the air entry value in the classical vanGenuchten model [-]
        fc_vgm   (1:nl_soil), & ! a scaling factor by using air entry value in the Mualem model [-]
#endif
        porsl(1:nl_soil) , &! saturated volumetric soil water content(porosity)
        psi0(1:nl_soil)  , &! saturated soil suction (mm) (NEGATIVE)
        hksati(1:nl_soil), &! hydraulic conductivity at saturation (mm h2o/s)
        rootr(1:nl_soil) , &! root resistance of a layer, all layers add to 1.0

        t_soisno(lb:nl_soil), &! soil/snow skin temperature (K)
        pg_rain          , &! rainfall after removal of interception (mm h2o/s)
        sm               , &! snow melt (mm h2o/s)
        etr              , &! actual transpiration (mm h2o/s)
        qseva            , &! ground surface evaporation rate (mm h2o/s)
        qsdew            , &! ground surface dew formation (mm h2o /s) [+]
        qsubl            , &! sublimation rate from snow pack (mm h2o /s) [+]
        qfros               ! surface dew added to snow pack (mm h2o /s) [+]
#if(defined CaMa_Flood)
  real(r8), INTENT(inout) :: flddepth  ! inundation water input from top (mm/s)
  real(r8), INTENT(in)    :: fldfrc    ! inundation water input from top (mm/s)
  real(r8), INTENT(out)   :: qinfl_fld ! inundation water input from top (mm/s)
#endif

  real(r8), INTENT(inout) :: &
        wice_soisno(lb:nl_soil) , &! ice lens (kg/m2)
        wliq_soisno(lb:nl_soil) , &! liquid water (kg/m2)
        smp(1:nl_soil)   , &! soil matrix potential [mm]
        hk (1:nl_soil)   , &! hydraulic conductivity [mm h2o/m]
        zwt              , &! the depth from ground (soil) surface to water table [m]
        wdsrf            , &! depth of surface water [mm]
        wa               , &! water storage in aquifer [mm]
        wetwat              ! water storage in wetland [mm]

  real(r8), INTENT(out) :: &
        rsur             , &! surface runoff (mm h2o/s)
        rnof             , &! total runoff (mm h2o/s)
        qinfl            , &! infiltration rate (mm h2o/s)
        qcharge             ! groundwater recharge (positive to aquifer) [mm/s]

  real(r8), INTENT(out) :: errw_rsub ! the possible subsurface runoff dificit after PHS is included

! SNICAR model variables
! Aerosol Fluxes (Jan. 07, 2023)
  real(r8), intent(in) :: forc_aer ( 14 )  ! aerosol deposition from atmosphere model (grd,aer) [kg m-1 s-1]

  real(r8), INTENT(inout) :: &
        mss_bcpho (lb:0), &! mass of hydrophobic BC in snow  (col,lyr) [kg]
        mss_bcphi (lb:0), &! mass of hydrophillic BC in snow (col,lyr) [kg]
        mss_ocpho (lb:0), &! mass of hydrophobic OC in snow  (col,lyr) [kg]
        mss_ocphi (lb:0), &! mass of hydrophillic OC in snow (col,lyr) [kg]
        mss_dst1  (lb:0), &! mass of dust species 1 in snow  (col,lyr) [kg]
        mss_dst2  (lb:0), &! mass of dust species 2 in snow  (col,lyr) [kg]
        mss_dst3  (lb:0), &! mass of dust species 3 in snow  (col,lyr) [kg]
        mss_dst4  (lb:0)   ! mass of dust species 4 in snow  (col,lyr) [kg]
! Aerosol Fluxes (Jan. 07, 2023)
! END SNICAR model variables

!-----------------------Local Variables------------------------------
!
  integer j                 ! loop counter

  real(r8) :: &
  eff_porosity(1:nl_soil), &! effective porosity = porosity - vol_ice
       gwat              , &! net water input from top (mm/s)
       drainmax          , &! drainage max (mm h2o/s)
       rsubst            , &! subsurface runoff (mm h2o/s)
       vol_liq(1:nl_soil), &! partial volume of liquid water in layer
       vol_ice(1:nl_soil), &! partial volume of ice lens in layer
       icefrac(1:nl_soil)   ! ice fraction (-)

  real(r8) :: err_solver, w_sum, wresi(1:nl_soil)
  REAL(r8) :: qgtop

  REAL(r8) :: zwtmm
  REAL(r8) :: sp_zc(1:nl_soil), sp_zi(0:nl_soil), sp_dz(1:nl_soil) ! in mm
  LOGICAL  :: is_permeable(1:nl_soil)
  REAL(r8) :: dzsum, dz
  REAL(r8) :: icefracsum, fracice_rsub, imped

#ifdef CROP
   integer  :: ps, pe
   integer  :: irrig_flag  ! 1 if sprinker, 2 if others
   real(r8) :: qflx_irrig_drip
   real(r8) :: qflx_irrig_sprinkler
   real(r8) :: qflx_irrig_flood
   real(r8) :: qflx_irrig_paddy
#endif

#ifdef Campbell_SOIL_MODEL
  real(r8) :: theta_r(1:nl_soil)
#endif

#ifdef Campbell_SOIL_MODEL
  INTEGER, parameter :: nprms = 1
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
  INTEGER, parameter :: nprms = 5
#endif
  REAL(r8) :: prms(nprms, 1:nl_soil)
#if(defined CaMa_Flood)
  real(r8) :: gfld,qinfl_all,rsur_fld, qinfl_fld_subgrid! inundation water input from top (mm/s)
#endif

#ifdef Campbell_SOIL_MODEL
      theta_r(1:nl_soil) = 0._r8
#endif

!=======================================================================
! [1] update the liquid water within snow layer and the water onto soil
!=======================================================================

      if (lb>=1)then
         gwat = pg_rain + sm - qseva
      else

         IF (.not. DEF_USE_SNICAR .or. (patchtype==1 .and. DEF_URBAN_RUN)) THEN
            call snowwater (lb,deltim,ssi,wimp,&
                         pg_rain,qseva,qsdew,qsubl,qfros,&
                         dz_soisno(lb:0),wice_soisno(lb:0),wliq_soisno(lb:0),gwat)
         ELSE
            call snowwater_snicar (lb,deltim,ssi,wimp,&
                         pg_rain,qseva,qsdew,qsubl,qfros,&
                         dz_soisno(lb:0),wice_soisno(lb:0),wliq_soisno(lb:0),gwat,&
                         forc_aer,&
                         mss_bcpho(lb:0), mss_bcphi(lb:0), mss_ocpho(lb:0), mss_ocphi(lb:0),&
                         mss_dst1(lb:0), mss_dst2(lb:0), mss_dst3(lb:0), mss_dst4(lb:0) )
         ENDIF
      endif

#ifdef CROP
      if(DEF_USE_IRRIGATION)then
         if(patchtype==0)then
            ps = patch_pft_s(ipatch)
            pe = patch_pft_e(ipatch)
            call CalIrrigationApplicationFluxes(ipatch,ps,pe,deltim,qflx_irrig_drip,qflx_irrig_sprinkler,qflx_irrig_flood,qflx_irrig_paddy,irrig_flag=2)
            gwat = gwat + qflx_irrig_drip + qflx_irrig_flood + qflx_irrig_paddy
         end if
      end if
#endif
!=======================================================================
! [2] surface runoff and infiltration
!=======================================================================

  if(patchtype<=1)then   ! soil ground only

      ! For water balance check, the sum of water in soil column before the calcultion
      w_sum = sum(wliq_soisno(1:nl_soil)) + sum(wice_soisno(1:nl_soil)) + wa + wdsrf

      ! Due to the increase in volume after freezing, the total volume of water and
      ! ice may exceed the porosity of the soil. This excess water is temporarily
      ! stored in "wresi". After calculating the movement of soil water, "wresi"
      ! is added back to "wliq_soisno".
      wresi(1:nl_soil) = 0.
      ! porosity of soil, partitial volume of ice and liquid
      do j = 1, nl_soil
         vol_ice(j) = min(porsl(j), wice_soisno(j)/(dz_soisno(j)*denice))
         if(porsl(j) < 1.e-6)then
            icefrac(j) = 0.
         else
            icefrac(j) = min(1.,vol_ice(j)/porsl(j))
         endif

         eff_porosity(j) = max(wimp, porsl(j)-vol_ice(j))
         is_permeable(j) = eff_porosity(j) > max(wimp, theta_r(j))
         IF (is_permeable(j)) THEN
            vol_liq(j) = wliq_soisno(j)/(dz_soisno(j)*denh2o)
            vol_liq(j) = min(eff_porosity(j), max(0., vol_liq(j)))
            wresi(j) = wliq_soisno(j) - dz_soisno(j) * denh2o * vol_liq(j)
         ELSE
            vol_liq(j) = 0.
         ENDIF
      enddo

      ! surface runoff including water table and surface staturated area

#ifndef LATERAL_FLOW
      if (gwat > 0.) then
         call surfacerunoff (nl_soil,wtfact,wimp,porsl,psi0,hksati,&
                             z_soisno(1:),dz_soisno(1:),zi_soisno(0:),&
                             eff_porosity,icefrac,zwt,gwat,rsur)
      else
         rsur = 0.
      endif

      ! infiltration into surface soil layer
      qgtop = gwat - rsur
#else
      ! for lateral flow, "rsur" is calculated in HYDRO/MOD_Hydro_SurfaceFlow.F90
      ! and is removed from surface water there.
      qgtop = gwat
#endif

#if(defined CaMa_Flood)
      IF (LWINFILT) then
         ! \ re-infiltration [mm/s] calculation.
         ! if surface runoff is ocurred (rsur != 0.), flood depth <1.e-6  and flood frction <0.05,
         ! the re-infiltration will not be calculated.
         IF ((flddepth .GT. 1.e-6).and.(fldfrc .GT. 0.05) .and. (patchtype == 0) ) then
            !           write(6,*) 'flddepth=',flddepth,'fldfrc=',fldfrc
            gfld=flddepth/deltim ! [mm/s]
            ! surface runoff from inundation, this should not be added to the surface runoff from soil
            ! otherwise, the surface runoff will be double counted.
            ! only the re-infiltration is added to water balance calculation.
            CALL surfacerunoff (nl_soil,1.0,wimp,porsl,psi0,hksati,&
               z_soisno(1:),dz_soisno(1:),zi_soisno(0:),&
               eff_porosity,icefrac,zwt,gfld,rsur_fld)
            ! infiltration into surface soil layer
            qinfl_fld_subgrid = gfld - rsur_fld !assume the re-infiltration is occured in whole patch area.
            !            write(6,*) 'gfld=',gfld,'   qinfl_fld_subgrid=',qinfl_fld_subgrid
         ELSE
            qinfl_fld_subgrid=0.0d0
            gfld=0.0d0
            rsur_fld=0.0d0

         ENDIF
         qinfl_fld=qinfl_fld_subgrid*fldfrc ! [mm/s] re-infiltration in grid.
         qgtop=qinfl_fld+qgtop ! [mm/s] total infiltration in grid.
         flddepth=flddepth-deltim*qinfl_fld_subgrid ! renew flood depth [mm], the flood depth is reduced by re-infiltration but only in inundation area.
      ENDIF
#endif
!=======================================================================
! [3] determine the change of soil water
!=======================================================================

      ! convert length units from m to mm
      zwtmm = zwt * 1000.0
      sp_zc(1:nl_soil) = z_soisno (1:nl_soil) * 1000.0   ! from meter to mm
      sp_zi(0:nl_soil) = zi_soisno(0:nl_soil) * 1000.0   ! from meter to mm

      ! check consistancy between water table location and liquid water content
      DO j = 1, nl_soil
         IF ((vol_liq(j) < eff_porosity(j)-1.e-8) .and. (zwtmm <= sp_zi(j-1))) THEN
            zwtmm = sp_zi(j)
         ENDIF
      ENDDO

#ifndef LATERAL_FLOW
      !-- Topographic runoff  ----------------------------------------------------------
      imped = 1.0
      IF (zwtmm < sp_zi(nl_soil)) THEN
         dzsum = 0.
         icefracsum = 0.
         do j = nl_soil, 1, -1
            sp_dz(j) = sp_zi(j) - sp_zi(j-1)
            IF (zwtmm < sp_zi(j)) THEN
               dz = min(sp_dz(j), sp_zi(j)-zwtmm)
               dzsum = dzsum + dz
               icefracsum = icefracsum + icefrac(j) * dz
            ENDIF
         end do
         ! add ice impedance factor to baseflow
         IF (dzsum > 0) THEN
            fracice_rsub = max(0.,exp(-3.*(1.-(icefracsum/dzsum)))-exp(-3.))/(1.0-exp(-3.))
            imped = max(0.,1.-fracice_rsub)
         ENDIF
      ENDIF

      rsubst = imped * 5.5e-3 * exp(-2.5*zwt)  ! drainage (positive = out of soil column)
#ifdef DataAssimilation
      rsubst = rsubst * fslp_patch(ipatch)
#endif

#else
      ! for lateral flow:
      ! "rsub" is calculated and removed from soil water in HYDRO/MOD_Hydro_SubsurfaceFlow.F90
      rsubst = 0
#endif


#ifdef Campbell_SOIL_MODEL
      prms(1,:) = bsw(1:nl_soil)
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
      prms(1,1:nl_soil) = alpha_vgm(1:nl_soil)
      prms(2,1:nl_soil) = n_vgm    (1:nl_soil)
      prms(3,1:nl_soil) = L_vgm    (1:nl_soil)
      prms(4,1:nl_soil) = sc_vgm   (1:nl_soil)
      prms(5,1:nl_soil) = fc_vgm   (1:nl_soil)
#endif

      ! update "vol_liq" in the level containing water table
      ! "vol_liq" in this level refers to volume content in unsaturated part
      IF (zwtmm < sp_zi(nl_soil)) THEN
         DO j = nl_soil, 1, -1
            IF ((zwtmm >= sp_zi(j-1)) .and. (zwtmm < sp_zi(j))) THEN

               IF ((zwtmm > sp_zi(j-1)) .and. (is_permeable(j))) THEN
                  vol_liq(j) = (wliq_soisno(j)*1000.0/denh2o - eff_porosity(j)*(sp_zi(j)-zwtmm))  &
                     / (zwtmm - sp_zi(j-1))
                  IF (vol_liq(j) < 0.) THEN
                     zwtmm = sp_zi(j)
                     vol_liq(j) = wliq_soisno(j)*1000.0/denh2o / (sp_zi(j) - sp_zi(j-1))
                  ENDIF

                  vol_liq(j) = max(0., min(eff_porosity(j), vol_liq(j)))
                  wresi(j) = wliq_soisno(j)*1000.0/denh2o - eff_porosity(j)*(sp_zi(j)-zwtmm) &
                     - vol_liq(j) * (zwtmm - sp_zi(j-1))
               ENDIF

               exit
            ENDIF
         ENDDO
      ENDIF

      wdsrf = max(0., wdsrf)

      IF ((.not. is_permeable(1)) .and. (qgtop < 0.)) THEN
         IF (wdsrf > 0) THEN
            wdsrf = wdsrf + qgtop * deltim
            IF (wdsrf < 0) THEN
               wliq_soisno(1) = max(0., wliq_soisno(1) + wdsrf)
               wdsrf = 0
            ENDIF
         ELSE
            wliq_soisno(1) = max(0., wliq_soisno(1) + qgtop * deltim)
         ENDIF

         qgtop = 0.

      ENDIF

      CALL soil_water_vertical_movement ( &
         nl_soil, deltim, sp_zc(1:nl_soil), sp_zi(0:nl_soil), is_permeable(1:nl_soil),    &
         eff_porosity(1:nl_soil), theta_r(1:nl_soil), psi0(1:nl_soil), hksati(1:nl_soil), &
         nprms, prms(:,1:nl_soil), porsl(nl_soil),     &
         qgtop, etr, rootr(1:nl_soil), rsubst, qinfl, &
         wdsrf, zwtmm, wa, vol_liq(1:nl_soil), smp(1:nl_soil), hk(1:nl_soil), 1.e-3)

      ! update the mass of liquid water
      DO j = nl_soil, 1, -1
         IF (is_permeable(j)) THEN
            IF (zwtmm < sp_zi(j)) THEN
               IF (zwtmm >= sp_zi(j-1)) THEN
                  wliq_soisno(j)  = denh2o * ((eff_porosity(j)*(sp_zi(j)-zwtmm))  &
                     + vol_liq(j) * (zwtmm - sp_zi(j-1)))/1000.0
               ELSE
                  wliq_soisno(j)  = denh2o * (eff_porosity(j)*(sp_zi(j)-sp_zi(j-1)))/1000.0
               ENDIF
            ELSE
               wliq_soisno(j) = denh2o * (vol_liq(j)*(sp_zi(j)-sp_zi(j-1)))/1000.0
            ENDIF

            wliq_soisno(j) = wliq_soisno(j) + wresi(j)
         ENDIF
      ENDDO

      zwt = zwtmm/1000.0

      ! Renew the ice and liquid mass due to condensation
      if(lb >= 1)then
         ! make consistent with how evap_grnd removed in infiltration
         wliq_soisno(1) = max(0., wliq_soisno(1) + qsdew * deltim)
         wice_soisno(1) = max(0., wice_soisno(1) + (qfros-qsubl) * deltim)
      end if

#ifndef LATERAL_FLOW
     IF (wdsrf > pondmx) THEN
        rsur = rsur + (wdsrf - pondmx) / deltim
        wdsrf = pondmx
     ENDIF

     ! total runoff (mm/s)
     rnof = rsubst + rsur
#endif

#ifndef LATERAL_FLOW
      err_solver = (sum(wliq_soisno(1:))+sum(wice_soisno(1:))+wa+wdsrf) - w_sum &
         - (gwat-etr-rsur-rsubst)*deltim
#else
      err_solver = (sum(wliq_soisno(1:))+sum(wice_soisno(1:))+wa+wdsrf) - w_sum &
         - (gwat-etr)*deltim
#endif
      if(lb >= 1)then
         err_solver = err_solver - (qsdew+qfros-qsubl)*deltim
      endif
#if(defined CaMa_Flood)
      IF (LWINFILT) THEN
         err_solver = err_solver-(gfld-rsur_fld)*fldfrc*deltim
      ENDIF
#endif
#if(defined CoLMDEBUG)
      if(abs(err_solver) > 1.e-3)then
         write(6,'(A,E20.5)') 'Warning (WATER_VSF): water balance violation', err_solver,ipatch
      endif
      IF (any(wliq_soisno < -1.e-3)) THEN
         write(6,'(A,10E20.5)') 'Warning (WATER_VSF): negative soil water', wliq_soisno(1:nl_soil)
      ENDIF
#endif

!=======================================================================
! [6] assumed hydrological scheme for the wetland
!=======================================================================

  else
      if(patchtype==2)then        ! WETLAND
         qinfl = 0.
         zwt = 0.
         qcharge = 0.
         
         IF (lb >= 1) THEN
            wetwat = wdsrf + wa + wetwat + (gwat - etr + qsdew + qfros - qsubl) * deltim
         ELSE
            wetwat = wdsrf + wa + wetwat + (gwat - etr) * deltim
         ENDIF

         IF (wetwat > wetwatmax) THEN
            wdsrf  = wetwat - wetwatmax
            wetwat = wetwatmax
            wa     = 0.
         ELSEIF (wetwat < 0) THEN
            wa     = wetwat
            wdsrf  = 0.
            wetwat = 0.
         ELSE
            wdsrf = 0.
            wa    = 0.
         ENDIF
      
#ifndef LATERAL_FLOW
         IF (wdsrf > pondmx) THEN
            rsur = rsur + (wdsrf - pondmx) / deltim
            wdsrf = pondmx
         ELSE
            rsur = 0.
         ENDIF
         rnof = rsur
#endif
      endif

  endif

  errw_rsub = 0.

  end subroutine WATER_VSF


  subroutine snowwater (lb,deltim,ssi,wimp, &
                        pg_rain,qseva,qsdew,qsubl,qfros, &
                        dz_soisno,wice_soisno,wliq_soisno,qout_snowb)

!-----------------------------------------------------------------------
! Original author : Yongjiu Dai, /09/1999; /04/2014
!
! Water flow wihtin snow is computed by an explicit and non-physical based scheme,
! which permits a part of liquid water over the holding capacity (a tentative value
! is used, i.e., equal to 0.033*porosity) to percolate into the underlying layer,
! except the case of that the porosity of one of the two neighboring layers is
! less than 0.05, the zero flow is assumed. The water flow out of the bottom
! snow pack will participate as the input of the soil water and runoff.
!
!-----------------------------------------------------------------------

  use MOD_Precision
  use MOD_Const_Physical, only : denice, denh2o  ! physical constant
  implicit none

!----------------------- dummy argument --------------------------------
  integer, INTENT(in) :: &
        lb          ! lower bound of array

  real(r8), INTENT(in) :: &
        deltim,    &! seconds in a time step (s)
        ssi,       &! irreducible water saturation of snow
        wimp,      &! water impremeable if porosity less than wimp
        dz_soisno(lb:0),  &! layer thickness (m)

        pg_rain,   &! rainfall after removal of interception (mm h2o/s)
        qseva,     &! ground surface evaporation rate (mm h2o/s)
        qsdew,     &! ground surface dew formation (mm h2o /s) [+]
        qsubl,     &! sublimation rate from snow pack (mm h2o /s) [+]
        qfros       ! surface dew added to snow pack (mm h2o /s) [+]

  real(r8), INTENT(inout) :: &
        wice_soisno(lb:0),&! ice lens (kg/m2)
        wliq_soisno(lb:0)  ! liquid water (kg/m2)

  real(r8), INTENT(out) :: &
        qout_snowb  ! rate of water out of snow bottom (mm/s)

!----------------------- local variables --------------------------------
  integer j         ! k do loop/array indices

  real(r8) :: &
       qin,        &! water flow into the elmement (mm/s)
       qout,       &! water flow out of the elmement (mm/s)
       zwice,      &! the sum of ice mass of snow cover (kg/m2)
       wgdif,      &! ice mass after minus sublimation
    vol_liq(lb:0), &! partitial volume of liquid water in layer
    vol_ice(lb:0), &! partitial volume of ice lens in layer
 eff_porosity(lb:0) ! effective porosity = porosity - vol_ice

!=======================================================================
! renew the mass of ice lens (wice_soisno) and liquid (wliq_soisno) in the surface snow layer,
! resulted by sublimation (frost) / evaporation (condense)

      wgdif = wice_soisno(lb) + (qfros - qsubl)*deltim
      wice_soisno(lb) = wgdif
      if(wgdif < 0.)then
         wice_soisno(lb) = 0.
         wliq_soisno(lb) = wliq_soisno(lb) + wgdif
      endif
      wliq_soisno(lb) = wliq_soisno(lb) + (pg_rain + qsdew - qseva)*deltim
      wliq_soisno(lb) = max(0., wliq_soisno(lb))

! Porosity and partitial volume
      do j = lb, 0
         vol_ice(j) = min(1., wice_soisno(j)/(dz_soisno(j)*denice))
         eff_porosity(j) = max(0.01, 1. - vol_ice(j))
         vol_liq(j) = min(eff_porosity(j), wliq_soisno(j)/(dz_soisno(j)*denh2o))
      enddo

! Capillary force within snow could be two or more orders of magnitude
! less than those of gravity, this term may be ignored.
! Here we could keep the garavity term only. The genernal expression
! for water flow is "K * ss**3", however, no effective paramterization
! for "K". Thus, a very simple treatment (not physical based) is introduced:
! when the liquid water of layer exceeds the layer's holding
! capacity, the excess meltwater adds to the underlying neighbor layer.

      qin = 0.
      do j= lb, 0
         wliq_soisno(j) = wliq_soisno(j) + qin

         if(j <= -1)then
         ! no runoff over snow surface, just ponding on surface
           if(eff_porosity(j)<wimp .OR. eff_porosity(j+1)<wimp)then
             qout = 0.
           else
             qout = max(0.,(vol_liq(j)-ssi*eff_porosity(j))*dz_soisno(j))
             qout = min(qout,(1.-vol_ice(j+1)-vol_liq(j+1))*dz_soisno(j+1))
           endif
         else
           qout = max(0.,(vol_liq(j)-ssi*eff_porosity(j))*dz_soisno(j))
         endif

         qout = qout*1000.
         wliq_soisno(j) = wliq_soisno(j) - qout
         qin = qout

      enddo

      qout_snowb = qout/deltim

  end subroutine snowwater


!-----------------------------------------------------------------------
  subroutine SnowWater_snicar (lb,deltim,ssi,wimp, &
                        pg_rain,qseva,qsdew,qsubl,qfros, &
                        dz_soisno,wice_soisno,wliq_soisno,qout_snowb, &

! Aerosol Fluxes (Jan. 07, 2023)
                        forc_aer, &
                        mss_bcpho, mss_bcphi, mss_ocpho, mss_ocphi, &
                        mss_dst1,  mss_dst2,  mss_dst3,  mss_dst4 )
! Aerosol Fluxes (Jan. 07, 2023)


!-----------------------------------------------------------------------
! Original author : Yongjiu Dai, /09/1999, /04/2014, /01/2023/
!
! Water flow wihtin snow is computed by an explicit and non-physical based scheme,
! which permits a part of liquid water over the holding capacity (a tentative value
! is used, i.e., equal to 0.033*porosity) to percolate into the underlying layer,
! except the case of that the porosity of one of the two neighboring layers is
! less than 0.05, the zero flow is assumed. The water flow out of the bottom
! snow pack will participate as the input of the soil water and runoff.
!
! REVISIONS:
! Yongjiu Dai, 01/2023: added Aerosol fluxes from SNICAR model
!-----------------------------------------------------------------------

  IMPLICIT NONE

  real(r8), parameter :: denice = 917.0_r8  ! density of ice [kg/m3]
  real(r8), parameter :: denh2o = 1000.0_r8 ! density of liquid water [kg/m3]

!----------------------- dummy argument --------------------------------
  integer, INTENT(in) :: &
        lb          ! lower bound of array

  real(r8), INTENT(in) :: &
        deltim,    &! seconds in a time step (s)
        ssi,       &! irreducible water saturation of snow
        wimp,      &! water impremeable if porosity less than wimp
        dz_soisno(lb:0),  &! layer thickness (m)

        pg_rain,   &! rainfall after removal of interception (mm h2o/s)
        qseva,     &! ground surface evaporation rate (mm h2o/s)
        qsdew,     &! ground surface dew formation (mm h2o /s) [+]
        qsubl,     &! sublimation rate from snow pack (mm h2o /s) [+]
        qfros       ! surface dew added to snow pack (mm h2o /s) [+]

  real(r8), INTENT(inout) :: &
        wice_soisno(lb:0),&! ice lens (kg/m2)
        wliq_soisno(lb:0)  ! liquid water (kg/m2)

  real(r8), INTENT(out) :: &
        qout_snowb  ! rate of water out of snow bottom (mm/s)

! Aerosol Fluxes (Jan. 07, 2023)
  real(r8), intent(in) :: forc_aer ( 14 )  ! aerosol deposition from atmosphere model (grd,aer) [kg m-1 s-1]

  real(r8), INTENT(inout) :: &
        mss_bcpho (lb:0), &! mass of hydrophobic BC in snow  (col,lyr) [kg]
        mss_bcphi (lb:0), &! mass of hydrophillic BC in snow (col,lyr) [kg]
        mss_ocpho (lb:0), &! mass of hydrophobic OC in snow  (col,lyr) [kg]
        mss_ocphi (lb:0), &! mass of hydrophillic OC in snow (col,lyr) [kg]
        mss_dst1  (lb:0), &! mass of dust species 1 in snow  (col,lyr) [kg]
        mss_dst2  (lb:0), &! mass of dust species 2 in snow  (col,lyr) [kg]
        mss_dst3  (lb:0), &! mass of dust species 3 in snow  (col,lyr) [kg]
        mss_dst4  (lb:0)   ! mass of dust species 4 in snow  (col,lyr) [kg]
! Aerosol Fluxes (Jan. 07, 2023)

!----------------------- local variables --------------------------------
  integer j         ! do loop/array indices

  real(r8) :: &
       qin,        &! water flow into the elmement (mm/s)
       qout,       &! water flow out of the elmement (mm/s)
       zwice,      &! the sum of ice mass of snow cover (kg/m2)
       wgdif,      &! ice mass after minus sublimation
    vol_liq(lb:0), &! partitial volume of liquid water in layer
    vol_ice(lb:0), &! partitial volume of ice lens in layer
 eff_porosity(lb:0) ! effective porosity = porosity - vol_ice

! Aerosol Fluxes (Jan. 07, 2023)
  !  Aerosol species indices:
  !  1= hydrophillic (bulk model) or within-ice (modal model) black carbon
  !  2= hydrophobic (bulk model) or external (modal model) black carbon
  !  3= hydrophilic organic carbon
  !  4= hydrophobic organic carbon
  !  5= dust species 1
  !  6= dust species 2
  !  7= dust species 3
  !  8= dust species 4
  !
  real(r8), parameter :: scvng_fct_mlt_bcphi = 0.20 ! scavenging factor for hydrophillic BC inclusion in meltwater [frc]
  real(r8), parameter :: scvng_fct_mlt_bcpho = 0.03 ! scavenging factor for hydrophobic BC inclusion in meltwater  [frc]
  real(r8), parameter :: scvng_fct_mlt_ocphi = 0.20 ! scavenging factor for hydrophillic OC inclusion in meltwater [frc]
  real(r8), parameter :: scvng_fct_mlt_ocpho = 0.03 ! scavenging factor for hydrophobic OC inclusion in meltwater  [frc]
  real(r8), parameter :: scvng_fct_mlt_dst1  = 0.02 ! scavenging factor for dust species 1 inclusion in meltwater  [frc]
  real(r8), parameter :: scvng_fct_mlt_dst2  = 0.02 ! scavenging factor for dust species 2 inclusion in meltwater  [frc]
  real(r8), parameter :: scvng_fct_mlt_dst3  = 0.01 ! scavenging factor for dust species 3 inclusion in meltwater  [frc]
  real(r8), parameter :: scvng_fct_mlt_dst4  = 0.01 ! scavenging factor for dust species 4 inclusion in meltwater  [frc]

  ! !LOCAL VARIABLES:
  real(r8) :: qin_bc_phi       ! flux of hydrophilic BC into   layer [kg]
  real(r8) :: qout_bc_phi      ! flux of hydrophilic BC out of layer [kg]
  real(r8) :: qin_bc_pho       ! flux of hydrophobic BC into   layer [kg]
  real(r8) :: qout_bc_pho      ! flux of hydrophobic BC out of layer [kg]
  real(r8) :: qin_oc_phi       ! flux of hydrophilic OC into   layer [kg]
  real(r8) :: qout_oc_phi      ! flux of hydrophilic OC out of layer [kg]
  real(r8) :: qin_oc_pho       ! flux of hydrophobic OC into   layer [kg]
  real(r8) :: qout_oc_pho      ! flux of hydrophobic OC out of layer [kg]
  real(r8) :: qin_dst1         ! flux of dust species 1 into   layer [kg]
  real(r8) :: qout_dst1        ! flux of dust species 1 out of layer [kg]
  real(r8) :: qin_dst2         ! flux of dust species 2 into   layer [kg]
  real(r8) :: qout_dst2        ! flux of dust species 2 out of layer [kg]
  real(r8) :: qin_dst3         ! flux of dust species 3 into   layer [kg]
  real(r8) :: qout_dst3        ! flux of dust species 3 out of layer [kg]
  real(r8) :: qin_dst4         ! flux of dust species 4 into   layer [kg]
  real(r8) :: qout_dst4        ! flux of dust species 4 out of layer [kg]
  real(r8) :: mss_liqice(lb:0) ! mass of liquid+ice in a layer

  real(r8) :: subsnow          ! sublimated snow [kg m-2]
  real(r8) :: frc_sub          ! fraction of layer mass that has sublimated [frc]
  real(r8) :: frc_transfer     ! frc_refrz + frc_sub
  real(r8) :: dm_int           ! mass transfer [kg]

  ! !LOCAL VARIABLES for AerosolFluxes
  real(r8) :: flx_bc_dep       ! total BC deposition        (col) [kg m-2 s-1]
  real(r8) :: flx_bc_dep_phi   ! hydrophillic BC deposition (col) [kg m-1 s-1]
  real(r8) :: flx_bc_dep_pho   ! hydrophobic BC deposition  (col) [kg m-1 s-1]
  real(r8) :: flx_oc_dep       ! total OC deposition        (col) [kg m-2 s-1]
  real(r8) :: flx_oc_dep_phi   ! hydrophillic OC deposition (col) [kg m-1 s-1]
  real(r8) :: flx_oc_dep_pho   ! hydrophobic OC deposition  (col) [kg m-1 s-1]
  real(r8) :: flx_dst_dep      ! total dust deposition      (col) [kg m-2 s-1]

  real(r8) :: flx_dst_dep_wet1 ! wet dust (species 1) deposition (col) [kg m-2 s-1]
  real(r8) :: flx_dst_dep_dry1 ! dry dust (species 1) deposition (col) [kg m-2 s-1]
  real(r8) :: flx_dst_dep_wet2 ! wet dust (species 2) deposition (col) [kg m-2 s-1]
  real(r8) :: flx_dst_dep_dry2 ! dry dust (species 2) deposition (col) [kg m-2 s-1]
  real(r8) :: flx_dst_dep_wet3 ! wet dust (species 3) deposition (col) [kg m-2 s-1]
  real(r8) :: flx_dst_dep_dry3 ! dry dust (species 3) deposition (col) [kg m-2 s-1]
  real(r8) :: flx_dst_dep_wet4 ! wet dust (species 4) deposition (col) [kg m-2 s-1]
  real(r8) :: flx_dst_dep_dry4 ! dry dust (species 4) deposition (col) [kg m-2 s-1]
! Aerosol Fluxes (Jan. 07, 2023)

!=======================================================================
! renew the mass of ice lens (wice_soisno) and liquid (wliq_soisno) in the surface snow layer,
! resulted by sublimation (frost) / evaporation (condense)

      wgdif = wice_soisno(lb) + (qfros - qsubl)*deltim
      wice_soisno(lb) = wgdif
      if(wgdif < 0.)then
         wice_soisno(lb) = 0.
         wliq_soisno(lb) = wliq_soisno(lb) + wgdif
      endif
      wliq_soisno(lb) = wliq_soisno(lb) + (pg_rain + qsdew - qseva)*deltim
      wliq_soisno(lb) = max(0., wliq_soisno(lb))

! Porosity and partitial volume
      do j = lb, 0
         vol_ice(j) = min(1., wice_soisno(j)/(dz_soisno(j)*denice))
         eff_porosity(j) = max(0.01, 1. - vol_ice(j))
         vol_liq(j) = min(eff_porosity(j), wliq_soisno(j)/(dz_soisno(j)*denh2o))
      enddo

! Capillary force within snow could be two or more orders of magnitude
! less than those of gravity, this term may be ignored.
! Here we could keep the garavity term only. The genernal expression
! for water flow is "K * ss**3", however, no effective paramterization
! for "K". Thus, a very simple treatment (not physical based) is introduced:
! when the liquid water of layer exceeds the layer's holding
! capacity, the excess meltwater adds to the underlying neighbor layer.

! Aerosol Fluxes (Jan. 07, 2023)
! Also compute aerosol fluxes through snowpack in this loop:
! 1) compute aerosol mass in each layer
! 2) add aerosol mass flux from above layer to mass of this layer
! 3) qout_xxx is mass flux of aerosol species xxx out bottom of
!    layer in water flow, proportional to (current) concentration
!    of aerosol in layer multiplied by a scavenging ratio.
! 4) update mass of aerosol in top layer, accordingly
! 5) update mass concentration of aerosol accordingly

      qin        = 0._r8

! Aerosol Fluxes (Jan. 07, 2023)
      qin_bc_phi = 0._r8
      qin_bc_pho = 0._r8
      qin_oc_phi = 0._r8
      qin_oc_pho = 0._r8
      qin_dst1   = 0._r8
      qin_dst2   = 0._r8
      qin_dst3   = 0._r8
      qin_dst4   = 0._r8
! Aerosol Fluxes (Jan. 07, 2023)

      do j= lb, 0

         wliq_soisno(j) = wliq_soisno(j) + qin

! Aerosol Fluxes (Jan. 07, 2023)
         mss_bcphi(j) = mss_bcphi(j) + qin_bc_phi
         mss_bcpho(j) = mss_bcpho(j) + qin_bc_pho
         mss_ocphi(j) = mss_ocphi(j) + qin_oc_phi
         mss_ocpho(j) = mss_ocpho(j) + qin_oc_pho

         mss_dst1(j)  = mss_dst1(j) + qin_dst1
         mss_dst2(j)  = mss_dst2(j) + qin_dst2
         mss_dst3(j)  = mss_dst3(j) + qin_dst3
         mss_dst4(j)  = mss_dst4(j) + qin_dst4
! Aerosol Fluxes (Jan. 07, 2023)

         if(j <= -1)then
         ! no runoff over snow surface, just ponding on surface
           if(eff_porosity(j)<wimp .OR. eff_porosity(j+1)<wimp)then
             qout = 0._r8
           else
             qout = max(0._r8,(vol_liq(j)-ssi*eff_porosity(j))*dz_soisno(j))
             qout = min(qout,(1.-vol_ice(j+1)-vol_liq(j+1))*dz_soisno(j+1))
           endif
         else
           qout = max(0._r8,(vol_liq(j)-ssi*eff_porosity(j))*dz_soisno(j))
         endif

         qout = qout*1000._r8
         wliq_soisno(j) = wliq_soisno(j) - qout
         qin = qout

! Aerosol Fluxes (Jan. 07, 2023)
         ! mass of ice+water: in extremely rare circumstances, this can
         ! be zero, even though there is a snow layer defined. In
         ! this case, set the mass to a very small value to
         ! prevent division by zero.

         mss_liqice(j) = wliq_soisno(j)+wice_soisno(j)
         if (mss_liqice(j) < 1E-30_r8) then
            mss_liqice(j) = 1E-30_r8
         endif

         ! BCPHI:
         ! 1. flux with meltwater:
         qout_bc_phi = qout*scvng_fct_mlt_bcphi*(mss_bcphi(j)/mss_liqice(j))
         if (qout_bc_phi > mss_bcphi(j)) then
            qout_bc_phi = mss_bcphi(j)
         endif
         mss_bcphi(j) = mss_bcphi(j) - qout_bc_phi
         qin_bc_phi = qout_bc_phi

         ! BCPHO:
         ! 1. flux with meltwater:
         qout_bc_pho = qout*scvng_fct_mlt_bcpho*(mss_bcpho(j)/mss_liqice(j))
         if (qout_bc_pho > mss_bcpho(j)) then
            qout_bc_pho = mss_bcpho(j)
         endif
         mss_bcpho(j) = mss_bcpho(j) - qout_bc_pho
         qin_bc_pho = qout_bc_pho

         ! OCPHI:
         ! 1. flux with meltwater:
         qout_oc_phi = qout*scvng_fct_mlt_ocphi*(mss_ocphi(j)/mss_liqice(j))
         if (qout_oc_phi > mss_ocphi(j)) then
            qout_oc_phi = mss_ocphi(j)
         endif
         mss_ocphi(j) = mss_ocphi(j) - qout_oc_phi
         qin_oc_phi = qout_oc_phi

         ! OCPHO:
         ! 1. flux with meltwater:
         qout_oc_pho = qout*scvng_fct_mlt_ocpho*(mss_ocpho(j)/mss_liqice(j))
         if (qout_oc_pho > mss_ocpho(j)) then
            qout_oc_pho = mss_ocpho(j)
         endif
         mss_ocpho(j) = mss_ocpho(j) - qout_oc_pho
         qin_oc_pho = qout_oc_pho

         ! DUST 1:
         ! 1. flux with meltwater:
         qout_dst1 = qout*scvng_fct_mlt_dst1*(mss_dst1(j)/mss_liqice(j))
         if (qout_dst1 > mss_dst1(j)) then
            qout_dst1 = mss_dst1(j)
         endif
         mss_dst1(j) = mss_dst1(j) - qout_dst1
         qin_dst1 = qout_dst1

         ! DUST 2:
         ! 1. flux with meltwater:
         qout_dst2 = qout*scvng_fct_mlt_dst2*(mss_dst2(j)/mss_liqice(j))
         if (qout_dst2 > mss_dst2(j)) then
            qout_dst2 = mss_dst2(j)
         endif
         mss_dst2(j) = mss_dst2(j) - qout_dst2
         qin_dst2 = qout_dst2

         ! DUST 3:
         ! 1. flux with meltwater:
         qout_dst3 = qout*scvng_fct_mlt_dst3*(mss_dst3(j)/mss_liqice(j))
         if (qout_dst3 > mss_dst3(j)) then
            qout_dst3 = mss_dst3(j)
         endif
         mss_dst3(j) = mss_dst3(j) - qout_dst3
         qin_dst3 = qout_dst3

         ! DUST 4:
         ! 1. flux with meltwater:
         qout_dst4 = qout*scvng_fct_mlt_dst4*(mss_dst4(j)/mss_liqice(j))
         if (qout_dst4 > mss_dst4(j)) then
            qout_dst4 = mss_dst4(j)
         endif
         mss_dst4(j) = mss_dst4(j) - qout_dst4
         qin_dst4 = qout_dst4
! Aerosol Fluxes (Jan. 07, 2023)

      enddo

      qout_snowb = qout/deltim


! Aerosol Fluxes (Jan. 07, 2023)
! Compute aerosol fluxes through snowpack and aerosol deposition fluxes into top layere
!-----------------------------------------------------------------------
! set aerosol deposition fluxes from forcing array
! The forcing array is either set from an external file
! or from fluxes received from the atmosphere model
#ifdef MODAL_AER
    ! Mapping for modal aerosol scheme where within-hydrometeor and
    ! interstitial aerosol fluxes are differentiated. Here, "phi"
    ! flavors of BC and OC correspond to within-hydrometeor
    ! (cloud-borne) aerosol, and "pho" flavors are interstitial
    ! aerosol. "wet" and "dry" fluxes of BC and OC specified here are
    ! purely diagnostic
    !
    ! NOTE: right now the macro 'MODAL_AER' is not defined anywhere, i.e.,
    ! the below (modal aerosol scheme) is not available and can not be
    ! active either. It depends on the specific input aerosol deposition
    ! data which is suitable for modal scheme. [06/15/2023, Hua Yuan]

    flx_bc_dep_phi   = forc_aer(3)
    flx_bc_dep_pho   = forc_aer(1) + forc_aer(2)
    flx_bc_dep       = forc_aer(1) + forc_aer(2) + forc_aer(3)

    flx_oc_dep_phi   = forc_aer(6)
    flx_oc_dep_pho   = forc_aer(4) + forc_aer(5)
    flx_oc_dep       = forc_aer(4) + forc_aer(5) + forc_aer(6)

    flx_dst_dep_wet1 = forc_aer(7)
    flx_dst_dep_dry1 = forc_aer(8)
    flx_dst_dep_wet2 = forc_aer(9)
    flx_dst_dep_dry2 = forc_aer(10)
    flx_dst_dep_wet3 = forc_aer(11)
    flx_dst_dep_dry3 = forc_aer(12)
    flx_dst_dep_wet4 = forc_aer(13)
    flx_dst_dep_dry4 = forc_aer(14)
    flx_dst_dep      = forc_aer(7)  + forc_aer(8)  + forc_aer(9) + &
                       forc_aer(10) + forc_aer(11) + forc_aer(12) + &
                       forc_aer(13) + forc_aer(14)
#else
    ! Original mapping for bulk aerosol deposition. phi and pho BC/OC
    ! species are distinguished in model, other fluxes (e.g., dry and
    ! wet BC/OC) are purely diagnostic.

    flx_bc_dep_phi   = forc_aer(1) + forc_aer(3)
    flx_bc_dep_pho   = forc_aer(2)
    flx_bc_dep       = forc_aer(1) + forc_aer(2) + forc_aer(3)

    flx_oc_dep_phi   = forc_aer(4) + forc_aer(6)
    flx_oc_dep_pho   = forc_aer(5)
    flx_oc_dep       = forc_aer(4) + forc_aer(5) + forc_aer(6)

    flx_dst_dep_wet1 = forc_aer(7)
    flx_dst_dep_dry1 = forc_aer(8)
    flx_dst_dep_wet2 = forc_aer(9)
    flx_dst_dep_dry2 = forc_aer(10)
    flx_dst_dep_wet3 = forc_aer(11)
    flx_dst_dep_dry3 = forc_aer(12)
    flx_dst_dep_wet4 = forc_aer(13)
    flx_dst_dep_dry4 = forc_aer(14)
    flx_dst_dep      = forc_aer(7)  + forc_aer(8)  + forc_aer(9) + &
                       forc_aer(10) + forc_aer(11) + forc_aer(12) + &
                       forc_aer(13) + forc_aer(14)
#endif

    ! aerosol deposition fluxes into top layer
    ! This is done after the inter-layer fluxes so that some aerosol
    ! is in the top layer after deposition, and is not immediately
    ! washed out before radiative calculations are done

    mss_bcphi(lb) = mss_bcphi(lb) + (flx_bc_dep_phi*deltim)
    mss_bcpho(lb) = mss_bcpho(lb) + (flx_bc_dep_pho*deltim)
    mss_ocphi(lb) = mss_ocphi(lb) + (flx_oc_dep_phi*deltim)
    mss_ocpho(lb) = mss_ocpho(lb) + (flx_oc_dep_pho*deltim)

    mss_dst1(lb) = mss_dst1(lb) + (flx_dst_dep_dry1 + flx_dst_dep_wet1)*deltim
    mss_dst2(lb) = mss_dst2(lb) + (flx_dst_dep_dry2 + flx_dst_dep_wet2)*deltim
    mss_dst3(lb) = mss_dst3(lb) + (flx_dst_dep_dry3 + flx_dst_dep_wet3)*deltim
    mss_dst4(lb) = mss_dst4(lb) + (flx_dst_dep_dry4 + flx_dst_dep_wet4)*deltim

#ifdef MODAL_AER
    !
    ! Transfer BC and OC from the within-ice state to the external
    ! state based on snow sublimation and re-freezing of liquid water.
    ! Re-freezing effect is inactived by default because of
    ! uncertainty in how this process operates.

    do j= lb, 0
       if (j >= lb) then
          if (j == lb) then
             ! snow that has sublimated [kg/m2] (top layer only)
             subsnow = max(0._r8, (qsubl*deltim))

             ! fraction of layer mass that has sublimated:
             if ((wliq_soisno(j) + wice_soisno(j)) > 0._r8) then
                frc_sub = subsnow / (wliq_soisno(j) + wice_soisno(j))
             else
                frc_sub = 0._r8
             endif
          else
             ! prohibit sublimation effect to operate on sub-surface layers:
             frc_sub = 0._r8
          endif

          ! fraction of layer mass transformed (sublimation only)
          frc_transfer = frc_sub

          ! cap the fraction at 1
          if (frc_transfer > 1._r8) then
             frc_transfer = 1._r8
          endif

          ! transfer proportionate mass of BC and OC:
          dm_int       = mss_bcphi(j)*frc_transfer
          mss_bcphi(j) = mss_bcphi(j) - dm_int
          mss_bcpho(j) = mss_bcpho(j) + dm_int

          dm_int       = mss_ocphi(j)*frc_transfer
          mss_ocphi(j) = mss_ocphi(j) - dm_int
          mss_ocpho(j) = mss_ocpho(j) + dm_int

       end if
    end do
#endif
! Aerosol Fluxes (Jan. 7, 2023)

  end subroutine SnowWater_snicar

  subroutine surfacerunoff (nl_soil,wtfact,wimp,porsl,psi0,hksati,&
                            z_soisno,dz_soisno,zi_soisno,&
                            eff_porosity,icefrac,zwt,gwat,rsur)

!=======================================================================
! the original code was provide by Robert E. Dickinson based on following clues:
! a water table level determination level added including highland and
! lowland levels and fractional area of wetland (water table above the surface.
! Runoff is parametrized from the lowlands in terms of precip incident on
! wet areas and a base flow, where these are estimated using ideas from TOPMODEL.
!
! Author : Yongjiu Dai, 07/29/2002, Guoyue Niu, 06/2012
!=======================================================================

  use MOD_Precision
  implicit none

!-----------------------Arguments---------------------------------------

  integer, INTENT(in) :: nl_soil   ! number of soil layers
  real(r8), INTENT(in) :: &
        wtfact,        &! fraction of model area with high water table
        wimp,          &! water impremeable if porosity less than wimp
     porsl(1:nl_soil), &! saturated volumetric soil water content(porosity)
      psi0(1:nl_soil), &! saturated soil suction (mm) (NEGATIVE)
    hksati(1:nl_soil), &! hydraulic conductivity at saturation (mm h2o/s)
  z_soisno(1:nl_soil), &! layer depth (m)
 dz_soisno(1:nl_soil), &! layer thickness (m)
 zi_soisno(0:nl_soil), &! interface level below a "z" level (m)
 eff_porosity(1:nl_soil), &! effective porosity = porosity - vol_ice
   icefrac(1:nl_soil), &! ice fraction (-)
        gwat,          &! net water input from top
         zwt            ! the depth from ground (soil) surface to water table [m]

  real(r8), INTENT(out) :: rsur    ! surface runoff (mm h2o/s)

!-----------------------Local Variables---------------------------------

  real(r8) qinmax       ! maximum infiltration capability
  real(r8) fsat         ! fractional area with water table at surface

  real(r8), parameter :: fff = 0.5   ! runoff decay factor (m-1)

!-----------------------End Variable List-------------------------------

!  fraction of saturated area
      fsat = wtfact*min(1.0,exp(-0.5*fff*zwt))

! Maximum infiltration capacity
      qinmax = minval(10.**(-6.0*icefrac(1:min(3,nl_soil)))*hksati(1:min(3,nl_soil)))
      if(eff_porosity(1)<wimp) qinmax = 0.

! Surface runoff
      rsur = fsat*max(0.0,gwat) + (1.-fsat)*max(0.,gwat-qinmax)


  end subroutine surfacerunoff



  subroutine soilwater(patchtype,nl_soil,deltim,wimp,smpmin,&
                       qinfl,etr,z_soisno,dz_soisno,zi_soisno,&
                       t_soisno,vol_liq,vol_ice,smp,hk,icefrac,eff_porosity,&
                       porsl,hksati,bsw,psi0,rootr,&
                       zwt,dwat,qcharge)

!-----------------------------------------------------------------------
! Original author : Yongjiu Dai, 09/1999, 04/2014, 07/2014
!
! some new parameterization are added, which are based on CLM4.5
!
! Soil moisture is predicted from a 10-layer model (as with soil
! temperature), in which the vertical soil moisture transport is governed
! by infiltration, runoff, gradient diffusion, gravity, and root
! extraction through canopy transpiration. The net water applied to the
! surface layer is the snowmelt plus precipitation plus the throughfall
! of canopy dew minus surface runoff and evaporation.
!
! The vertical water flow in an unsaturated porous media is described by
! Darcy's law, and the hydraulic conductivity and the soil negative
! potential vary with soil water content and soil texture based on the work
! of Clapp and Hornberger (1978) and Cosby et al. (1984). The equation is
! integrated over the layer thickness, in which the time rate of change in
! water mass must equal the net flow across the bounding interface, plus the
! rate of internal source or sink. The terms of water flow across the layer
! interfaces are linearly expanded by using first-order Taylor expansion.
! The equations result in a tridiagonal system equation.
!
! Note: length units here are all millimeter
! (in temperature subroutine uses same soil layer
! structure required but lengths are m)
!
! Richards equation:
!
! d wat     d     d psi
! ----- =  -- [ k(----- - 1) ] + S
!   dt     dz       dz
!
! where: wat = volume of water per volume of soil (mm**3/mm**3)
! psi = soil matrix potential (mm)
! dt  = time step (s)
! z   = depth (mm) (positive downward)
! dz  = thickness (mm)
! qin = inflow at top (mm h2o /s)
! qout= outflow at bottom (mm h2o /s)
! s   = source/sink flux (mm h2o /s)
! k   = hydraulic conductivity (mm h2o /s)
!
!                       d qin                  d qin
! qin[n+1] = qin[n] +  --------  d wat(j-1) + --------- d wat(j)
!                       d wat(j-1)             d wat(j)
!                ==================|=================
!                                  < qin
!
!                 d wat(j)/dt * dz = qin[n+1] - qout[n+1] + S(j)
!
!                                  > qout
!                ==================|=================
!                        d qout               d qout
! qout[n+1] = qout[n] + --------- d wat(j) + --------- d wat(j+1)
!                        d wat(j)             d wat(j+1)
!
!
! Solution: linearize k and psi about d wat and use tridiagonal
! system of equations to solve for d wat,
! where for layer j
!
!
! r_j = a_j [d wat_j-1] + b_j [d wat_j] + c_j [d wat_j+1]
!
!-----------------------------------------------------------------------
    use MOD_Precision
    use MOD_Const_Physical , only : grav,hfus,tfrz,denh2o,denice
    USE MOD_Utils

    IMPLICIT NONE

    INTEGER , intent(in) :: patchtype ! land patch type
    integer , INTENT(in) :: nl_soil   ! number of soil layers
    real(r8), INTENT(in) :: deltim    ! land model time step (sec)
    real(r8), INTENT(in) :: wimp      ! water impremeable if porosity less than wimp
    real(r8), INTENT(in) :: smpmin    ! restriction for min of soil potential (mm)

    real(r8), INTENT(in) :: qinfl     ! infiltration (mm H2O /s)
    real(r8), INTENT(in) :: etr       ! vegetation transpiration (mm H2O/s) (+ = to atm)

    real(r8), INTENT(in) :: z_soisno (1:nl_soil) ! layer depth (m)
    real(r8), INTENT(in) :: dz_soisno(1:nl_soil) ! layer thickness (m)
    real(r8), INTENT(in) :: zi_soisno(0:nl_soil) ! interface level below a "z" level (m)

    real(r8), INTENT(in) :: t_soisno (1:nl_soil) ! soil temperature (Kelvin)
    real(r8), INTENT(in) :: vol_liq  (1:nl_soil) ! liquid volumetric water content
    real(r8), INTENT(in) :: vol_ice  (1:nl_soil) ! ice volumetric water content
    real(r8), INTENT(in) :: icefrac  (1:nl_soil)
    real(r8), INTENT(in) :: eff_porosity(1:nl_soil) ! effective porosity = porosity - vol_ice

    real(r8), INTENT(in) :: porsl  (1:nl_soil) ! volumetric soil water at saturation (porosity)
    real(r8), INTENT(in) :: hksati (1:nl_soil) ! hydraulic conductivity at saturation (mm H2O /s)
    real(r8), INTENT(in) :: bsw    (1:nl_soil) ! Clapp and Hornberger "b"
    real(r8), INTENT(in) :: psi0   (1:nl_soil) ! minimum soil suction (mm) [-]
    real(r8), INTENT(in) :: rootr  (1:nl_soil) ! effective fraction of roots in each soil layer
    real(r8), INTENT(in) :: zwt                ! the depth from ground (soil) surface to water table [m]

    real(r8), intent(out) :: dwat(1:nl_soil)   ! change of soil water [m3/m3]
    real(r8), INTENT(out) :: qcharge           ! aquifer recharge rate (positive to aquifer) (mm/s)
    real(r8), INTENT(inout) :: smp(1:nl_soil)  ! soil matrix potential [mm]
    real(r8), INTENT(inout) :: hk (1:nl_soil)  ! hydraulic conductivity [mm h2o/s]

!
! local arguments
!
    integer  :: j                 ! do loop indices
    real(r8) :: amx(1:nl_soil)    ! "a" left off diagonal of tridiagonal matrix
    real(r8) :: bmx(1:nl_soil)    ! "b" diagonal column for tridiagonal matrix
    real(r8) :: cmx(1:nl_soil)    ! "c" right off diagonal tridiagonal matrix
    real(r8) :: rmx(1:nl_soil)    ! "r" forcing term of tridiagonal matrix
    real(r8) :: zmm(1:nl_soil)    ! layer depth [mm]
    real(r8) :: dzmm(1:nl_soil)   ! layer thickness [mm]
    real(r8) :: zimm(0:nl_soil)   ! layer interface depth [mm]
    real(r8) :: den(1:nl_soil)    ! used in calculating qin, qout
    real(r8) :: alpha(1:nl_soil)  ! used in calculating qin, qout
    real(r8) :: qin(1:nl_soil)    ! flux of water into soil layer [mm h2o/s]
    real(r8) :: qout(1:nl_soil)   ! flux of water out of soil layer [mm h2o/s]
    real(r8) :: dqidw0(1:nl_soil) ! d(qin)/d(vol_liq(j-1))
    real(r8) :: dqidw1(1:nl_soil) ! d(qin)/d(vol_liq(j))
    real(r8) :: dqodw1(1:nl_soil) ! d(qout)/d(vol_liq(j))
    real(r8) :: dqodw2(1:nl_soil) ! d(qout)/d(vol_liq(j+1))
    real(r8) :: dsmpdw(1:nl_soil) ! d(smp)/d(vol_liq)
    real(r8) :: s_node            ! soil wetness
    real(r8) :: s1                ! "s" at interface of layer
    real(r8) :: s2                ! k*s**(2b+2)
    real(r8) :: dhkdw1(1:nl_soil) ! d(hk)/d(vol_liq(j))
    real(r8) :: dhkdw2(1:nl_soil) ! d(hk)/d(vol_liq(j+1))
    real(r8) :: imped(1:nl_soil)  !
    real(r8) :: errorw            ! mass balance error for this time step

    integer  :: jwt               ! index of the soil layer right above the water table (-)

    real(r8), parameter :: e_ice=6.0      !soil ice impedance factor
!-----------------------------------------------------------------------

    !compute jwt index
    ! The layer index of the first unsaturated layer,
    ! i.e., the layer right above the water table

    jwt = nl_soil
    ! allow jwt to equal zero when zwt is in top layer
    do j = 1, nl_soil
       if(zwt <= zi_soisno(j)) then
          jwt = j-1
          exit
       end if
    enddo

    ! Because the depths in this routine are in mm, use local
    ! variable arrays instead of pointers
    do j = 1, nl_soil
       zmm(j) = z_soisno(j)*1000.
       dzmm(j) = dz_soisno(j)*1000.
       zimm(j) = zi_soisno(j)*1000.
    end do

    zimm(0) = 0.0

    ! Compute matric potential and derivative based on liquid water content only
    do j = 1, nl_soil
       if(DEF_USE_PLANTHYDRAULICS .and. (patchtype/=1 .or. .not.DEF_URBAN_RUN))then
          if(t_soisno(j)>=tfrz) then
             if(porsl(j)<1.e-6)then     ! bed rock
                s_node = 0.001
                smp(j) = psi0(j)
                dsmpdw(j) = 0.
             else
                s_node = max(vol_liq(j)/porsl(j),0.01)
                s_node = min(1.0,s_node)
                smp(j) = psi0(j)*s_node**(-bsw(j))
                smp(j) = max(smpmin,smp(j))
                dsmpdw(j) = -bsw(j)*smp(j)/(s_node*porsl(j))
             endif
          else
             ! when ice is present, the matric potential is only related to temperature
             ! by (Fuchs et al., 1978: Soil Sci. Soc. Amer. J. 42(3):379-385)
             ! Unit 1 Joule = 1 (kg m2/s2), J/kg /(m/s2) ==> m ==> 1e3 mm
             smp(j) = 1.e3 * 0.3336e6/9.80616*(t_soisno(j)-tfrz)/t_soisno(j)
             smp(j) = max(smpmin, smp(j))        ! Limit soil suction
             dsmpdw(j) = 0.
          endif
       else
          if(t_soisno(j)>tfrz) then
             if(porsl(j)<1.e-6)then     ! bed rock
                s_node = 0.001
                smp(j) = psi0(j)
                dsmpdw(j) = 0.
             else
                s_node = max(vol_liq(j)/porsl(j),0.01)
                s_node = min(1.0,s_node)
                smp(j) = psi0(j)*s_node**(-bsw(j))
                smp(j) = max(smpmin,smp(j))
                dsmpdw(j) = -bsw(j)*smp(j)/(s_node*porsl(j))
             endif
          else
             ! when ice is present, the matric potential is only related to temperature
             ! by (Fuchs et al., 1978: Soil Sci. Soc. Amer. J. 42(3):379-385)
             ! Unit 1 Joule = 1 (kg m2/s2), J/kg /(m/s2) ==> m ==> 1e3 mm
             smp(j) = 1.e3 * 0.3336e6/9.80616*(t_soisno(j)-tfrz)/t_soisno(j)
             smp(j) = max(smpmin, smp(j))        ! Limit soil suction
             dsmpdw(j) = 0.
          endif
       end if
    end do

    ! Hydraulic conductivity and soil matric potential and their derivatives
    do j = 1, nl_soil

       if(j < nl_soil)then
          den(j) = (zmm(j+1)-zmm(j))
          alpha(j) = (smp(j+1)-smp(j))/den(j) - 1.
       else
          den(j) = 0.        ! not used
          alpha(j) = 0.      ! not used
       endif

       if((eff_porosity(j) < wimp) .OR. (eff_porosity(min(nl_soil,j+1)) < wimp) &
                                   .OR. (vol_liq(j) <= 1.e-3))then
          imped(j) = 0.
          hk(j) = 0.
          dhkdw1(j) = 0.
          dhkdw2(j) = 0.
       else
          ! The average conductivity between two heterogeneous medium layers (j and j + 1),
          ! are computed using different methods
          if(j < nl_soil)then
! Method I: UPSTREAM MEAN
             if(alpha(j) <= 0.)then
                hk(j) = hksati(j) * (vol_liq(j)/porsl(j))**(2.*bsw(j)+3.)
                dhkdw1(j) = hksati(j) * (2.*bsw(j)+3.)*(vol_liq(j)/porsl(j))**(2.*bsw(j)+2.)/porsl(j)
                dhkdw2(j) = 0.
             else
                hk(j) = hksati(j+1) * (vol_liq(j+1)/porsl(j+1))**(2.*bsw(j+1)+3.)
                dhkdw1(j) = 0.
                dhkdw2(j) = hksati(j+1) * (2.*bsw(j+1)+3.)*(vol_liq(j+1)/porsl(j+1))**(2.*bsw(j+1)+2.)/porsl(j+1)
             endif
! Method II:
          !  ! The harmonic averaging of the saturated conductivities
          !  hksat_interface = (zmm(j+1)-zmm(j))/((zimm(j)-zmm(j))/hksati(j)+(zmm(j+1)-zimm(j))/hksati(j+1))
          !  s1 = (vol_liq(j)*(zimm(j)-zmm(j)) + vol_liq(j+1)*(zmm(j+1)-zimm(j))) &
          !     / (porsl(j)*(zimm(j)-zmm(j)) + porsl(j+1)*(zmm(j+1)-zimm(j)))
          !  s1 = min(1.,s1)
          !  s2 = hksat_interface*s1**(2.*bsw(j)+2.)
          !  hk(j) = s1*s2
          !  dhkdw1(j) = (2.*bsw(j)+3.)*s2*(zimm(j)-zmm(j))/(porsl(j)*(zimm(j)-zmm(j))+porsl(j+1)*(zmm(j+1)-zimm(j)))
          !  dhkdw2(j) = (2.*bsw(j)+3.)*s2*(zmm(j+1)-zimm(j))/(porsl(j)*(zimm(j)-zmm(j))+porsl(j+1)*(zmm(j+1)-zimm(j)))

          else
             hk(j) = hksati(j) * (vol_liq(j)/porsl(j))**(2.*bsw(j)+3.)
             dhkdw1(j) = hksati(j) * (2.*bsw(j)+3.)*(vol_liq(j)/porsl(j))**(2.*bsw(j)+2.)/porsl(j)
             dhkdw2(j) = 0.
          endif

          ! replace fracice with impedance factor
          imped(j)=10.**(-e_ice*(0.5*(icefrac(j)+icefrac(min(nl_soil,j+1)))))
          hk(j) = imped(j) * hk(j)
          dhkdw1(j) = imped(j) * dhkdw1(j)
          dhkdw2(j) = imped(j) * dhkdw2(j)
       endif
    end do


    ! Set up r, a, b, and c vectors for tridiagonal solution

    ! Node j=1 (top)

    j = 1
    qin(j) = qinfl

    qout(j) = -hk(j)*alpha(j)
    dqodw1(j) = -(alpha(j)*dhkdw1(j) - hk(j)*dsmpdw(j)/den(j))
    dqodw2(j) = -(alpha(j)*dhkdw2(j) + hk(j)*dsmpdw(j+1)/den(j))

    amx(j) = 0.
    bmx(j) = dzmm(j)/deltim + dqodw1(j)
    cmx(j) = dqodw2(j)
    if(DEF_USE_PLANTHYDRAULICS .and. (patchtype/=1 .or. .not.DEF_URBAN_RUN))then
       rmx(j) =  qin(j) - qout(j) - rootr(j)
    else
       rmx(j) =  qin(j) - qout(j) - etr*rootr(j)
    end if

    ! Nodes j=2 to j=nl_soil-1

    do j = 2, nl_soil - 1
       qin(j) = -hk(j-1)*alpha(j-1)
       dqidw0(j) = -(alpha(j-1)*dhkdw1(j-1) - hk(j-1)*dsmpdw(j-1)/den(j-1))
       dqidw1(j) = -(alpha(j-1)*dhkdw2(j-1) + hk(j-1)*dsmpdw(j)/den(j-1))

       qout(j) = -hk(j)*alpha(j)
       dqodw1(j) = -(alpha(j)*dhkdw1(j) - hk(j)*dsmpdw(j)/den(j))
       dqodw2(j) = -(alpha(j)*dhkdw2(j) + hk(j)*dsmpdw(j+1)/den(j))

       amx(j) = -dqidw0(j)
       bmx(j) =  dzmm(j)/deltim - dqidw1(j) + dqodw1(j)
       cmx(j) =  dqodw2(j)
       if(DEF_USE_PLANTHYDRAULICS .and. (patchtype/=1 .or. .not.DEF_URBAN_RUN))then
          rmx(j) =  qin(j) - qout(j) - rootr(j)
       else
          rmx(j) =  qin(j) - qout(j) - etr*rootr(j)
       end if
    end do

    ! Node j=nl_soil (bottom)

    j = nl_soil
    qin(j) = -hk(j-1)*alpha(j-1)
    dqidw0(j) = -(alpha(j-1)*dhkdw1(j-1) - hk(j-1)*dsmpdw(j-1)/den(j-1))
    dqidw1(j) = -(alpha(j-1)*dhkdw2(j-1) + hk(j-1)*dsmpdw(j)/den(j-1))

!   if(j > jwt) then ! water table is in soil column
!      qout(j) = 0.
!      dqodw1(j) = 0.
!      dqodw2(j) = 0.
!   else
       qout(j) = hk(j)
       dqodw1(j) = dhkdw1(j)
       dqodw2(j) = 0.
!   endif

    amx(j) = -dqidw0(j)
    bmx(j) =  dzmm(j)/deltim - dqidw1(j) + dqodw1(j)
    cmx(j) =  dqodw2(j)
    if(DEF_USE_PLANTHYDRAULICS .and. (patchtype/=1 .or. .not.DEF_URBAN_RUN))then
       rmx(j) =  qin(j) - qout(j) - rootr(j)
    else
       rmx(j) =  qin(j) - qout(j) - etr*rootr(j)
    end if

    ! Solve for dwat

    call tridia (nl_soil, amx, bmx, cmx, rmx, dwat )

#if(defined CoLMDEBUG)
! The mass balance error (mm) for this time step is
    errorw = -deltim*(qin(1)-qout(nl_soil)-dqodw1(nl_soil)*dwat(nl_soil))
    do j = 1, nl_soil
       if(DEF_USE_PLANTHYDRAULICS .and. (patchtype/=1 .or. .not.DEF_URBAN_RUN))then
          errorw = errorw+dwat(j)*dzmm(j)+rootr(j)*deltim
       else
          errorw = errorw+dwat(j)*dzmm(j)+etr*rootr(j)*deltim
       end if
    enddo

    if(abs(errorw) > 1.e-3)then
       write(6,*) 'mass balance error in time step =',errorw
    endif
#endif

    ! Recharge rate qcharge to groundwater (positive to aquifer)
    qcharge = qout(nl_soil) + dqodw1(nl_soil)*dwat(nl_soil)


  end subroutine soilwater


  subroutine subsurfacerunoff(nl_soil,deltim,pondmx,&
                              eff_porosity,icefrac,&
                              dz_soisno,zi_soisno,wice_soisno,wliq_soisno,&
                              porsl,psi0,bsw,zwt,wa,&
                              qcharge,rsubst,errw_rsub)

! -------------------------------------------------------------------------


    use MOD_Precision
    use MOD_Const_Physical, only : tfrz
!
! ARGUMENTS:
    IMPLICIT NONE

    integer, INTENT(in) :: nl_soil       !
    real(r8), INTENT(in) :: deltim       ! land model time step (sec)
    real(r8), INTENT(in) :: pondmx       !

    real(r8), INTENT(in) :: eff_porosity(1:nl_soil) ! effective porosity = porosity - vol_ice
    real(r8), INTENT(in) :: icefrac(1:nl_soil)      ! ice fraction (-)

    real(r8), INTENT(in) :: dz_soisno  (1:nl_soil)  ! layer depth (m)
    real(r8), INTENT(in) :: zi_soisno  (0:nl_soil)  ! interface level below a "z" level (m)
    real(r8), INTENT(inout) :: wice_soisno(1:nl_soil)  ! ice lens (kg/m2)
    real(r8), INTENT(inout) :: wliq_soisno(1:nl_soil)  ! liquid water (kg/m2)

    real(r8), INTENT(in) :: porsl(1:nl_soil)        ! volumetric soil water at saturation (porosity)
    real(r8), INTENT(in) :: psi0(1:nl_soil)         ! minimum soil suction (mm) [-]
    real(r8), INTENT(in) :: bsw(1:nl_soil)          ! Clapp and Hornberger "b"

    real(r8), INTENT(inout) :: zwt       ! the depth from ground (soil) surface to water table [m]
    real(r8), INTENT(inout) :: wa        ! water in the unconfined aquifer (mm)
    real(r8), INTENT(in)    :: qcharge   ! aquifer recharge rate (positive to aquifer) (mm/s)
    real(r8), INTENT(out)   :: rsubst    ! drainage drainage (positive = out of soil column) (mm H2O /s)
    real(r8), INTENT(out)   :: errw_rsub ! the possible subsurface runoff dificit after PHS is included

!
! LOCAL ARGUMENTS
!

    integer  :: j                ! indices
    integer  :: jwt              ! index of the soil layer right above the water table (-)
    real(r8) :: xs               ! water needed to bring soil moisture to watmin (mm)
    real(r8) :: dzmm(1:nl_soil)  ! layer thickness (mm)
    real(r8) :: xsi              ! excess soil water above saturation at layer i (mm)
    real(r8) :: xsia             ! available pore space at layer i (mm)
    real(r8) :: xs1              ! excess soil water above saturation at layer 1 (mm)
    real(r8) :: ws               ! summation of pore space of layers below water table (mm)
    real(r8) :: s_node           ! soil wetness (-)
    real(r8) :: available_wliq_soisno     ! available soil liquid water in a layer
    real(r8) :: qcharge_tot      !
    real(r8) :: qcharge_layer    !
    real(r8) :: drainage         !
    real(r8) :: drainage_tot     !
    real(r8) :: drainage_layer   !
    real(r8) :: s_y              !
    real(r8) :: rous             ! specific yield [-]

    real(r8) :: wt
    real(r8) :: wtsub
    real(r8) :: dzsum
    real(r8) :: icefracsum
    real(r8) :: fracice_rsub
    real(r8) :: imped


    real(r8), parameter :: watmin = 0.01  ! Limit irreduciable wrapping liquid water
                                          ! a tunable constant
    real(r8), parameter :: rsbmx  = 5.0   ! baseflow coefficient [mm/s]
    real(r8), parameter :: timean = 10.5  ! global mean topographic index


! -------------------------------------------------------------------------

!   ! Convert layer thicknesses from m to mm

    do j = 1,nl_soil
       dzmm(j) = dz_soisno(j)*1000.
    end do

!   ! The layer index of the first unsaturated layer,
!   ! i.e., the layer right above the water table

    jwt = nl_soil
    ! allow jwt to equal zero when zwt is in top layer
    do j = 1, nl_soil
       if(zwt <= zi_soisno(j)) then
          jwt = j-1
          exit
       end if
    enddo

!============================== QCHARGE =========================================
! Water table changes due to qcharge
! use analytical expression for aquifer specific yield
    rous = porsl(nl_soil)*(1.-(1.-1.e3*zwt/psi0(nl_soil))**(-1./bsw(nl_soil)))
    rous = max(rous,0.02)

    wa = wa + qcharge*deltim
!
!---------------------------------------
    ! water table is below the soil column
    if(jwt == nl_soil) then
       zwt = max(0.,zwt - (qcharge*deltim)/1000./rous)
    else
    ! water table within soil layers 1-9
    ! try to raise water table to account for qcharge

       qcharge_tot = qcharge * deltim

       if(qcharge_tot > 0.) then ! rising water table
          do j = jwt+1, 1,-1
             ! use analytical expression for specific yield

             s_y = porsl(j) * (1.-(1.-1.e3*zwt/psi0(j))**(-1./bsw(j)))
             s_y=max(s_y,0.02)

             qcharge_layer = min(qcharge_tot,(s_y*(zwt-zi_soisno(j-1))*1.e3))
             qcharge_layer = max(qcharge_layer,0.)

             zwt = max(0.,zwt - qcharge_layer/s_y/1000.)

             qcharge_tot = qcharge_tot - qcharge_layer
             if (qcharge_tot <= 0.) exit
          enddo
       else ! deepening water table (negative qcharge)
          do j = jwt+1, nl_soil
             ! use analytical expression for specific yield
             s_y = porsl(j) * (1.-(1.-1.e3*zwt/psi0(j))**(-1./bsw(j)))
             s_y=max(s_y,0.02)
             qcharge_layer = max(qcharge_tot,-(s_y*(zi_soisno(j) - zwt)*1.e3))
             qcharge_layer = min(qcharge_layer,0.)
             qcharge_tot = qcharge_tot - qcharge_layer

             if (qcharge_tot >= 0.) then
                zwt = max(0.,zwt - qcharge_layer/s_y/1000.)
                exit
             else
                zwt = zi_soisno(j)
             endif
          enddo
          if (qcharge_tot > 0.) zwt = max(0.,zwt - qcharge_tot/1000./rous)
       endif
    endif

!-- Topographic runoff  ----------------------------------------------------------
    dzsum = 0.
    icefracsum = 0.
    do j = max(jwt,1), nl_soil
       dzsum = dzsum + dzmm(j)
       icefracsum = icefracsum + icefrac(j) * dzmm(j)
    end do
    ! add ice impedance factor to baseflow
    fracice_rsub = max(0.,exp(-3.*(1.-(icefracsum/dzsum)))-exp(-3.))/(1.0-exp(-3.))
    imped = max(0.,1.-fracice_rsub)
    drainage = imped * 5.5e-3 * exp(-2.5*zwt)  ! drainage (positive = out of soil column)

!-- Water table is below the soil column  ----------------------------------------
    if(jwt == nl_soil) then
       wa = wa - drainage * deltim
       zwt = max(0.,zwt + (drainage * deltim)/1000./rous)
       wliq_soisno(nl_soil) = wliq_soisno(nl_soil) + max(0.,(wa-5000.))
       wa = min(wa, 5000.)
    else
!-- Water table within soil layers 1-9  ------------------------------------------
!============================== RSUB_TOP =========================================
       !-- Now remove water via drainage
       drainage_tot = - drainage * deltim
       do j = jwt+1, nl_soil
          ! use analytical expression for specific yield
          s_y = porsl(j) * ( 1. - (1.-1.e3*zwt/psi0(j))**(-1./bsw(j)))
          s_y = max(s_y,0.02)

          drainage_layer = max(drainage_tot, -(s_y*(zi_soisno(j)-zwt)*1.e3))
          drainage_layer = min(drainage_layer,0.)
          wliq_soisno(j) = wliq_soisno(j) + drainage_layer

          drainage_tot = drainage_tot - drainage_layer

          if(drainage_tot >= 0.)then
             zwt = max(0.,zwt - drainage_layer/s_y/1000.)
             exit
          else
             zwt = zi_soisno(j)
          endif
       enddo

!-- Remove residual drainage  ------------------------------------------------
       zwt = max(0.,zwt - drainage_tot/1000./rous)
       wa = wa + drainage_tot

!-- Recompute jwt  ---------------------------------------------------------------
       ! allow jwt to equal zero when zwt is in top layer
       jwt = nl_soil
       do j = 1, nl_soil
          if(zwt <= zi_soisno(j)) then
             jwt = j-1
             exit
          end if
       enddo

    end if   ! end of jwt if construct

    zwt = max(0.0,zwt)
    zwt = min(80.,zwt)

    rsubst = drainage


    ! Correction [1]
    ! NON-physically based corection on wliq_soisno
    ! excessive water above saturation added to the above unsaturated layer like a bucket
    ! if column over saturated, excess water goes to runoff

    do j = nl_soil,2,-1
       xsi = max(wliq_soisno(j)-eff_porosity(j)*dzmm(j),0.)
       wliq_soisno(j) = min(eff_porosity(j)*dzmm(j), wliq_soisno(j))
       wliq_soisno(j-1) = wliq_soisno(j-1) + xsi
    end do

    ! 12/2022, note by yuan: a potential bug below which needs check,
    ! if wice_soisno(1) > pondmx + porsl*dzmm, so xs1>0, in that case,
    ! wliq_soisno(1) will be nagtive, and xs1 is positive.
    xs1 = wliq_soisno(1) - (pondmx+porsl(1)*dzmm(1)-wice_soisno(1))
    if(xs1 > 0.)then
       wliq_soisno(1) = pondmx+porsl(1)*dzmm(1)-wice_soisno(1)
    else
       xs1 = 0.
    endif

    rsubst = rsubst + xs1 / deltim


    ! Correction [2]
    ! NON-physically based corection on wliq_soisno
    ! Limit wliq_soisno to be greater than or equal to watmin.
    ! Get water needed to bring wliq_soisno equal watmin from lower layer.
    ! If insufficient water in soil layers, get from aquifer water

    xs = 0.
    do j = 1, nl_soil
       if (wliq_soisno(j) < 0.) then
          xs = xs + wliq_soisno(j)
          wliq_soisno(j) = 0.
       endif
    enddo

    ! Sub-surface runoff and drainage
    errw_rsub = min(0., rsubst + xs/deltim)
    rsubst = max(0., rsubst + xs/deltim)


!   do j = 1, nl_soil-1
!      if (wice_soisno(j)*wice_soisno(j+1) < 1.e-6)then
!         if (wliq_soisno(j) < watmin) then
!            xs = watmin - wliq_soisno(j)
!            ! deepen water table if water is passed from below zwt layer
!            if(j == jwt) then
!               zwt = max(0.,zwt + xs/eff_porosity(j)/1000.)
!            endif
!         else
!            xs = 0.
!         end if
!         wliq_soisno(j  ) = wliq_soisno(j  ) + xs
!         wliq_soisno(j+1) = wliq_soisno(j+1) - xs
!      endif
!   end do

!   ! Get water for bottom layer from layers above if possible
!   if(wliq_soisno(nl_soil) < watmin)then
!      xs = watmin-wliq_soisno(nl_soil)
!      do j = nl_soil-1, 1, -1
!         available_wliq_soisno = max(wliq_soisno(j)-watmin-xs,0.)
!         if(available_wliq_soisno >= xs)then
!            wliq_soisno(nl_soil) = wliq_soisno(nl_soil) + xs
!            wliq_soisno(j      ) = wliq_soisno(j      ) - xs
!            xs = 0.
!            exit
!         else
!            wliq_soisno(nl_soil) = wliq_soisno(nl_soil) + available_wliq_soisno
!            wliq_soisno(j      ) = wliq_soisno(j      ) - available_wliq_soisno
!            xs = xs - available_wliq_soisno
!         end if
!      end do
!   else
!      xs = 0.
!   end if

!   ! Needed in case there is no water to be found
!   wliq_soisno(nl_soil) = wliq_soisno(nl_soil) + xs

!   ! Sub-surface runoff and drainage
!   rsubst = rsubst - xs/deltim

  end subroutine subsurfacerunoff


END MODULE MOD_SoilSnowHydrology
! --------- EOP ----------
