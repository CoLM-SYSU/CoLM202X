#include <define.h>

MODULE SOIL_SNOW_hydrology

!-----------------------------------------------------------------------
  use precision
  IMPLICIT NONE
  SAVE

! PUBLIC MEMBER FUNCTIONS:
  public :: WATER
  public :: snowwater
  public :: soilwater


! PRIVATE MEMBER FUNCTIONS:
  private :: surfacerunoff
  private :: subsurfacerunoff


!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------



  subroutine WATER (ipatch,patchtype  ,lb          ,nl_soil ,deltim,&
             z_soisno    ,dz_soisno   ,zi_soisno                   ,&
             bsw         ,porsl       ,psi0        ,hksati  ,rootr ,&
             t_soisno    ,wliq_soisno ,wice_soisno ,pg_rain ,sm    ,&
             etr         ,qseva       ,qsdew       ,qsubl   ,qfros ,&
             rsur        ,rnof        ,qinfl       ,wtfact  ,pondmx,&
             ssi         ,wimp        ,smpmin      ,zwt     ,wa    ,&
             qcharge     )  

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

  use precision
  use PhysicalConstants, only: denice, denh2o, tfrz

  implicit none
    
!-----------------------Argument---------- ------------------------------
  integer, INTENT(in) :: &
        ipatch           ,& ! patch index
        patchtype           ! land water type (0=soil, 1=urban or built-up, 2=wetland, 
                            ! 3=land ice, 4=land water bodies, 99=ocean
  
  integer, INTENT(in) :: &
        lb               , &! lower bound of array
        nl_soil            ! upper bound of array

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

  real(r8), INTENT(inout) :: &
        wice_soisno(lb:nl_soil) , &! ice lens (kg/m2)
        wliq_soisno(lb:nl_soil) , &! liquid water (kg/m2)
        zwt              , &! the depth from ground (soil) surface to water table [m]
        wa                  ! water storage in aquifer [mm]

  real(r8), INTENT(out) :: &
        rsur             , &! surface runoff (mm h2o/s)
        rnof             , &! total runoff (mm h2o/s)
        qinfl            , &! infiltration rate (mm h2o/s)
        qcharge             ! groundwater recharge (positive to aquifer) [mm/s]
  
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

!=======================================================================
! [1] update the liquid water within snow layer and the water onto soil
!=======================================================================

      if (lb>=1)then
         gwat = pg_rain + sm - qseva
      else
         call snowwater (lb,deltim,ssi,wimp,&
                         pg_rain,qseva,qsdew,qsubl,qfros,&
                         dz_soisno(lb:0),wice_soisno(lb:0),wliq_soisno(lb:0),gwat)
      endif

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
      call surfacerunoff (nl_soil,wtfact,wimp,bsw,porsl,psi0,hksati,&
                          z_soisno(1:),dz_soisno(1:),zi_soisno(0:),&
                          eff_porosity,icefrac,zwt,gwat,rsur)
      else
           rsur = 0.
      endif

      ! infiltration into surface soil layer 
      qinfl = gwat - rsur 

!=======================================================================
! [3] determine the change of soil water
!=======================================================================

      ! convert length units from m to mm
      zmm(1:) = z_soisno(1:)*1000.
      dzmm(1:) = dz_soisno(1:)*1000.
      zimm(0:) = zi_soisno(0:)*1000.

      call soilwater(nl_soil,deltim,wimp,smpmin,&
                     qinfl,etr,z_soisno(1:),dz_soisno(1:),zi_soisno(0:),&
                     t_soisno(1:),vol_liq,vol_ice,icefrac,eff_porosity,&
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
                             qcharge,rsubst)

      ! total runoff (mm/s)
      rnof = rsubst + rsur                
 
      ! Renew the ice and liquid mass due to condensation
      if(lb >= 1)then
         ! make consistent with how evap_grnd removed in infiltration
         wliq_soisno(1) = max(0., wliq_soisno(1) + qsdew * deltim)
         wice_soisno(1) = max(0., wice_soisno(1) + (qfros-qsubl) * deltim)
      end if


      if(lb >= 1)then
         err_solver = (sum(wliq_soisno(1:))+sum(wice_soisno(1:))+wa) - w_sum &
                    - (gwat+qsdew+qfros-qsubl-etr-rnof)*deltim
      else
         err_solver = (sum(wliq_soisno(1:))+sum(wice_soisno(1:))+wa) - w_sum &
                    - (gwat-etr-rnof)*deltim
      endif

#if(defined CLMDEBUG)
     if(abs(err_solver) > 1.e-3)then
        write(6,*) 'Warning: water balance violation after all soilwater calculation', err_solver
     endif
#endif


!=======================================================================
! [6] assumed hydrological scheme for the wetland and glacier
!=======================================================================

  else                          
      if(patchtype==2)then        ! WETLAND
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

  endif

!-----------------------------------------------------------------------

  end subroutine WATER



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

  use precision
  use PhysicalConstants, only: denice, denh2o  ! physical constant
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



  subroutine surfacerunoff (nl_soil,wtfact,wimp,bsw,porsl,psi0,hksati,&
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

  use precision
  implicit none

!-----------------------Arguments---------------------------------------

  integer, INTENT(in) :: nl_soil   ! number of soil layers
  real(r8), INTENT(in) :: &
        wtfact,        &! fraction of model area with high water table
        wimp,          &! water impremeable if porosity less than wimp
       bsw(1:nl_soil), &! Clapp-Hornberger "B"
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
      qinmax = minval(10.**(-6.0*icefrac(1:3))*hksati(1:3))
      if(eff_porosity(1)<wimp) qinmax = 0.

! Surface runoff
      rsur = fsat*max(0.0,gwat) + (1.-fsat)*max(0.,gwat-qinmax)


  end subroutine surfacerunoff



  subroutine soilwater(nl_soil,deltim,wimp,smpmin,&
                       qinfl,etr,z_soisno,dz_soisno,zi_soisno,&
                       t_soisno,vol_liq,vol_ice,icefrac,eff_porosity,&
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
    use precision
    use PhysicalConstants , only: grav,hfus,tfrz,denh2o,denice

    IMPLICIT NONE

    integer, INTENT(in) :: nl_soil  ! number of soil layers
    real(r8), INTENT(in) :: deltim  ! land model time step (sec)
    real(r8), INTENT(in) :: wimp    ! water impremeable if porosity less than wimp
    real(r8), INTENT(in) :: smpmin  ! restriction for min of soil potential (mm)

    real(r8), INTENT(in) :: qinfl   ! infiltration (mm H2O /s)
    real(r8), INTENT(in) :: etr     ! vegetation transpiration (mm H2O/s) (+ = to atm)

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
    real(r8) :: smp(1:nl_soil)    ! soil matrix potential [mm]
    real(r8) :: hk(1:nl_soil)     ! hydraulic conductivity [mm h2o/s]
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
    rmx(j) = qin(j) - qout(j) - etr * rootr(j)

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
       rmx(j) =  qin(j) - qout(j) - etr*rootr(j)
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
    rmx(j) =  qin(j) - qout(j) - etr*rootr(j)

    ! Solve for dwat

    call tridia (nl_soil, amx, bmx, cmx, rmx, dwat )

#if(defined CLMDEBUG)
! The mass balance error (mm) for this time step is
    errorw = -deltim*(qin(1)-qout(nl_soil)-dqodw1(nl_soil)*dwat(nl_soil))
    do j = 1, nl_soil
       errorw = errorw+dwat(j)*dzmm(j)+etr*rootr(j)*deltim
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
                              qcharge,rsubst)

! -------------------------------------------------------------------------


    use precision
    use PhysicalConstants, only: tfrz
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


END MODULE SOIL_SNOW_hydrology
! --------- EOP ----------
