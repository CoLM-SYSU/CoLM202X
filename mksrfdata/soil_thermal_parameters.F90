!====================================================================
SUBROUTINE soil_thermal_parameters(gravel,sand,clay,SOC,BD,porsl,soildepth,&
                                   csol,ksatu,kdry)
!====================================================================
! Reference: Dai et al.,2014: Implementation of a New Global Soil Dataset in the Common Land Model.
!
! Created by Yongjiu Dai, 12/2013
! ----------------------------------------------------
use precision

IMPLICIT NONE
      real(r8), intent(in) :: gravel ! gravel percentage (% of volume)
      real(r8), intent(in) :: sand   ! sand percentage   (% of weight)
      real(r8), intent(in) :: clay   ! clay percentage   (% of weight)
      real(r8), intent(in) :: SOC    ! soil organic carbon (% of weight)
      real(r8), intent(in) :: BD     ! bulk density of dry soil material [g/cm^3]
      real(r8), intent(in) :: porsl  ! fraction of soil that is voids [-]
      real(r8), intent(in) :: soildepth ! (cm)

      real(r8), intent(out) :: csol  ! heat capacity of soil solids [J/(m3 K)]
      real(r8), intent(out) :: ksatu ! thermal conductivity of saturated unforzen soil [W/m-K]
      real(r8), intent(out) :: kdry  ! thermal conductivity for dry soil  [W/(m-K)]

      real(r8) csol_om      ! heat capacity of peat soil *10^6 (J//m3/K) (Farouki, 1986)
      real(r8) csol_gravel  ! heat capacity of gravel *10^6 (J//m3/K) (Farouki, 1986)
      real(r8) k_om         ! thermal conductivity of organic soil (Farouki, 1986) [W/m/K]
      real(r8) k_om_dry     ! thermal conductivity of dry organic soil (Farouki, 1981)
      real(r8) k_quartz     ! thermal conductivity of quartz (W/m/K)
      real(r8) k_minerals   ! thermal conductivity of non-quartz minerals (W/m/K)
      real(r8) k_water      ! thermal conductivity of liquid water (W/m/K)
      real(r8) k_ice        ! thermal conductivity of ice (W/m/K)
      real(r8) k_air        ! thermal conductivity of air (W/m/K)
      real(r8) zsapric      ! depth (m) that organic matter takes on characteristics of sapric peat
      real(r8) kmg          ! thermal conductivity of soil minerals [W/m-K]

      real(r8) om_watsat    ! porosity of organic soil
      real(r8) watsat       !
      real(r8) wf_om_s      ! weight fraction of soil organic matter iwithin the soil solids 
      real(r8) vf_om_s      ! volumetric fraction of soil organic matter iwithin the soil solids
      real(r8) vf_quartz_s  ! volumetric fraction of quartz within soil solids
      real(r8) vf_gravel    ! volumetric fraction of gravel 
      real(r8) v_pores      ! volumetric pore space of the soil

      real(r8) BD_minerals  ! (g/cm3)
      real(r8) BD_particle  ! particle density of the soil (g/cm3)
      real(r8) k_soild      ! thermal conductivity of the solid fraction of the soil (W/m/K)
      real(r8) csol_soil    ! volumetric heat capacity of soil (J/m3/K)

      real(r8) lambda_q     ! thermal conductivities of quartz (W /(m K))
      real(r8) lambda_o     ! thermal conductivities of other minerals (W /(m K))
      real(r8) lambda_w     ! thermal conductivity of water at 20C (W /(m K))
      real(r8) lambda_s     ! effective thermal conductivity of soil solids (W /(m K))

      real(r8) silt

!-----------------------------------------------------------------------
! The volumetric heat capacity (J/m3/K)
!-----------------------------------------------------------------------
      ! The weight fraction of the soil organic matter
      wf_om_s = 1.72*SOC/100.  
      ! The volumetric fraction of soil organic matter   
      vf_om_s = min(wf_om_s * BD / 1.3, 1.)   ! 1.3 (g/cm3) density of soil organic matter

      ! The volumetric fraction of soil gravel
      vf_gravel = gravel/100.

      csol_om = 2.51e6  
      csol_gravel = 2.01e6  ! (J/m3/K) volumetric  heat capacity of gravel
      csol_soil = (1.-vf_om_s)*( (2.128*sand+2.385*clay)/(sand+clay) )*1.e6 + vf_om_s*csol_om
      csol = vf_gravel*csol_gravel + (1.-vf_gravel)*csol_soil

! ---------------------------------------------------------------------------------------------------
! [1]
! Oleson K.W. et al., 2013: Technical Description of version 4.5 of the Community Land Model (CLM)
! NCAR/TN-503+STR (Section 6.3: Soil and Snow Thermal Properties)
! ---------------------------------------------------------------------------------------------------
      ! The weight fraction of the soil organic matter
      wf_om_s = 1.72*SOC/100.  
      ! The volumetric fraction of soil organic matter   
      vf_om_s = min(wf_om_s * BD / 1.3, 1.)   ! 1.3 (g/cm3) density of soil organic matter

      zsapric = 0.5      ! (m)
      om_watsat = max(0.93 - 0.1*(0.01*soildepth/zsapric), 0.83)
      watsat = (1.-vf_om_s)*porsl + vf_om_s*om_watsat
     
      k_om = 0.25        ! (W/m/K)
      k_om_dry = 0.05    ! (W/m/K)

      kmg   = ((1.-vf_om_s)*(8.80*sand+2.92*clay)/(sand+clay) + vf_om_s*k_om) ** (1.- watsat)  ! W/(m K)
      ksatu = kmg*0.57**watsat
      kdry  = (1.-vf_om_s)*(0.135*BD+0.0647)/(2.7-0.947*BD) + vf_om_s*k_om_dry

! --------------------------------------------------------------------------------------------------!
! [2]
! Balland V. and P. A. Arp, 2005: Modeling soil thermal conductivities over a wide
! range of conditions. J. Environ. Eng. Sci. 4: 549-558.
! ---------------------------------------------------------------------------------------------------
!     v_pores = porsl     ! volumetric pore space of the soil
!
!     ! The weight fraction of the soil organic matter
!     wf_om_s = 1.72*SOC/100.  
!     ! The volumetric fraction of soil organic matter   
!     vf_om_s = min(wf_om_s * BD / 1.3, 1.)   ! 1.3 (g/cm3) density of soil organic matter
!
!     ! The volumetric fraction of quartz (vf_quartz_s) within the soil solids
!     CALL vf_quartz(sand,clay,vf_quartz_s)
!
!     BD_minerals = 2.65  ! (g/cm3)
!     BD_particle = 1. / (wf_om_s/1.3 + (1.-wf_om_s)/BD_minerals )  ! particle density of the soil
!     
!     k_om = 0.25        ! (W/m/K)
!     k_quartz = 8.0     ! (W/m/K)
!     k_minerals = 2.5   ! (W/m/K)
!     k_water = 0.57     ! (W/m/K)
!     k_ice = 2.21       ! (W/m/K)
!     k_air = 0.024      ! (W/m/K)
!
!     k_soild = k_om**vf_om_s * k_quartz**vf_quartz_s * k_minerals**(1.-vf_om_s-vf_quartz_s)
!     ksatu = k_soild**(1.-v_pores) * k_water**v_pores  ! for saturated unforzen soil
!     kdry = ((0.053*k_soild-k_air)*BD + k_air*BD_particle) / (BD_particle - (1.-0.053)*BD)
!
!     for saturated unfrozen soils
!         k_sat = ksatu     
!     for saturated frozen soils
!         k_sat = k_soild**(1.-v_pores) * k_ice**(v_pores-v_water) * k_water**v_water  
!   
!         for unfrozen soil:
!             ke = w**(0.5*(1.+vf_om_s-alpha*v_sand-v_cf_s))
!                  *( (1./(1.+exp(-beta*w)))**3 - ((1.-w)/2.)**3))**(1.-vf_om_s)  
!             w = wetness or saturation (theta/theta_s)
!             alpha = adjustable parameter
!             beta = adjustable parameter
!             v_sand = volumetric fraction of sand within the soil solids
!             v_cf_s = volumetric fraction of coarse fragments within the soil solids
!         for fozen or partially frozen soils:
!             ke = w**(1.+vf_om_s)
!
! ---------------------------------------------------------------------------------------------------
! [3]
! Tarnawski et al (2012, 2009,...) gave a higher recommendation on Lu's (2006) equations. 
! (1) Tarnawski VR and WH Leong, 2012: A series-parallel model for estimating
! the thermal conductivity of unsaturated soils. Int J Thermophys (2012) 33:1191-1218]
! (2) Lu et al., 2007: An improved model for predicting soil thermal conductivity from water content 
! at room temperature. Soil Sci. Soc. Am. J. 71:8-14
! ---------------------------------------------------------------------------------------------------
!     v_pores = porsl     
!     CALL vf_quartz(sand,clay,vf_quartz_s)
!
!     kdry = 0.56*porsl + 0.51
!
!   ! ksatu followed the same procedure as Johansen model 
!   ! lambda_s = thermal conductivities of quartz: lambda_q = 7.7 W/(m K)
!   ! and other minerals: where lambda_o = 2.0 W/(m K) for soils with q > 0.2, and
!                               lambda_o = 3.0 W/(m K) for soils with q <= 0.2
!     lambda_q = 7.7
!     if(vf_quartz_s > 0.2)then
!        lambda_o = 2.0
!     else
!        lambda_o = 3.0
!     endif
!     lambda_w = 0.594  ! W /(m K) at 20C
!     lambda_s = lambda_q**vf_quartz_s*lambda_o**(1.0-vf_quartz_s)
!     ksatu = lambda_s**(1.-v_pores)*lambda_w**v_pores
!
!     ke = exp(alpha1*(1.-w**(alpha1-beta1)))
!          w = wetness or saturation (theta/theta_s)
!          alpha1 = 0.728; beta1 = 1.165 for coarse soils
!          alpha1 = 0.37 ; beta1 = 1.29 for fine soils 
!          coarse-textured soils = soils with sand fractions >40 (%), and 
!          fiine-textured soils = soils with sand fractions <40 (%)
!     k_soil = (k_sat - k_dry)*ke + k_dry

END SUBROUTINE soil_thermal_parameters



! ========================================
SUBROUTINE vf_quartz(sand,clay,vf_quartz_s)
! ========================================
! Table 2 (page 1212) of  Peters-Lidard, et al., 1998, The effect of soil thermal conductivity
! parameterization on surface energy fluxes and temperatures. J Atmos. Sci., Vol.55, 1209-1224.
!
! Created by Yongjiu Dai, 02/2014
! --------------------------------------------------------------------------------------------
use precision

IMPLICIT NONE
      real(r8), intent(in) :: sand
      real(r8), intent(in) :: clay
      real(r8), intent(out) :: vf_quartz_s  ! volumetric fraction of quartz within the soil solids

      real(r8) silt 
      integer, parameter :: PNUM=12  ! number of polygons(texture classes)
      logical c(PNUM)   ! indicate wheather a soil is in an class

      vf_quartz_s = 0.0
      silt = 100.-sand-clay

      if(sand<0. .or. silt<0. .or. clay<0.)then
         print*,'Each of the 3 variables should be >= 0: check the data'
         call abort
      end if
      if(sand+silt+clay<99. .or. sand+silt+clay>101.) then   ! the tolerance of the sum may be changed
         print*,'The sum of the 3 plotted variables should be around 100: check the data.'
         call abort
      end if

      CALL USDA_soil_classes(silt,clay,c)

! Quartz content 
      if(c(1))  vf_quartz_s = 0.25   ! clay
      if(c(2))  vf_quartz_s = 0.1    ! silty clay
      if(c(3))  vf_quartz_s = 0.52   ! sandy clay
      if(c(4))  vf_quartz_s = 0.35   ! clay loam
      if(c(5))  vf_quartz_s = 0.1    ! silty clay loam
      if(c(6))  vf_quartz_s = 0.6    ! sandy clay loam
      if(c(7))  vf_quartz_s = 0.4    ! loam
      if(c(8))  vf_quartz_s = 0.25   ! silty loam
      if(c(9))  vf_quartz_s = 0.6    ! sandy loam
      if(c(10)) vf_quartz_s = 0.1    ! silt
      if(c(11)) vf_quartz_s = 0.82   ! loamy sand
      if(c(12)) vf_quartz_s = 0.92   ! sand

END SUBROUTINE vf_quartz



SUBROUTINE USDA_soil_classes(x,y,c)
! ---------------------------------------------------------------------------------------------------
! USDA major soil textural classes based on the relative percentage of sand, silt and clay in the soil
!
! Initial Author : Wei Shangguan, 02/2014
! ---------------------------------------------------------------------------------------------------
use precision

IMPLICIT NONE
   integer, parameter :: TNUM=26   ! number of points in the triangle
   integer, parameter :: PNUM=12   ! number of polygons(texture class) in the triangle
   integer, parameter :: PONUM(PNUM)=(/5,3,4,6,4,5,5,8,7,4,4,3/)  ! number of points in a polygon (texture class)
   real(r8), intent(in) :: x           ! x(silt) of a soil
   real(r8), intent(in) :: y           ! y(clay) of a soil
   logical, intent(out) :: c(PNUM) ! indicate wheather a soil is in an class

   integer i,j
   real(r8) :: xpos(TNUM)              ! x(silt) coordinates of the  points in the triangle
   real(r8) :: ypos(TNUM)              ! y(clay) coordinates of the  points in the triangle
   integer :: points(PNUM,8)       ! sequence number of the points in a poygon (texture class)
                                   ! 8 is the maximun number of the points
   character(len=15) :: tnames(PNUM)  ! name of a texture class
   integer :: tcodes(PNUM)         ! code of a texture class, may be change accordingly
   real(r8) :: xpol(8)                 ! x(silt) coordinates of the  points in a poygon
   real(r8) :: ypol(8)                 ! y(clay) coordinates of the  points in a poygon

   xpos = (/ 0.0,  40.0,   0.0,  20.0,  15.0,  40.0,  60.0,   0.0,  27.5,  27.5,  50.0,  52.5,&
            72.5,   0.0,   0.0,  40.0,  50.0,  80.0,  87.5,  15.0,  30.0,  50.0,  80.0,   0.0,&
             0.0, 100.0/)
   ypos = (/55.0,  60.0,  35.0,  35.0,  40.0,  40.0,  40.0,  20.0,  20.0,  27.5,  27.5,  27.5,&
            27.5,  15.0,  10.0,   7.5,   7.5,  12.5,  12.5,   0.0,   0.0,   0.0,   0.0, 100.0,&
             0.0,   0.0/)

   points(1,1:PONUM(1))   = (/24,  1,  5,  6,  2/)
   points(2,1:PONUM(2))   = (/2, 6, 7/)
   points(3,1:PONUM(3))   = (/1, 3, 4, 5/)
   points(4,1:PONUM(4))   = (/5,  4, 10, 11, 12,  6/)
   points(5,1:PONUM(5))   = (/6, 12, 13,  7/)
   points(6,1:PONUM(6))   = (/3,  8,  9, 10,  4/)
   points(7,1:PONUM(7))   = (/10,  9, 16, 17, 11/)
   points(8,1:PONUM(8))   = (/11, 17, 22, 23, 18, 19, 13, 12/)
   points(9,1:PONUM(9))   = (/8, 14, 21, 22, 17, 16,  9/)
   points(10,1:PONUM(10)) = (/18, 23, 26, 19/)
   points(11,1:PONUM(11)) = (/14, 15, 20, 21/)
   points(12,1:PONUM(12)) = (/15, 25, 20/)

   tnames( 1) = 'clay           '
   tnames( 2) = 'silty clay     '
   tnames( 3) = 'sandy clay     '
   tnames( 4) = 'clay loam      '
   tnames( 5) = 'silty clay loam'
   tnames( 6) = 'sandy clay loam'
   tnames( 7) = 'loam           '
   tnames( 8) = 'silty loam     '
   tnames( 9) = 'sandy loam     '
   tnames(10) = 'silt           '
   tnames(11) = 'loamy sand     '
   tnames(12) = 'sand           '

   tcodes=(/1,2,3,4,5,6,7,8,9,10,11,12/)
!  -------------------------------------
   do i = 1, PNUM
      xpol(:) = 0
      do j = 1, PONUM(i)
         xpol(j) = xpos(points(i,j))
         ypol(j) = ypos(points(i,j))
      end do

      call pointinpolygon(x,y,xpol(1:PONUM(i)),ypol(1:PONUM(i)),PONUM(i),c(i))
   end do

END SUBROUTINE USDA_soil_classes



SUBROUTINE pointinpolygon(xp,yp,xpol,ypol,ponum,c)
! --------------------------------------------------------
! For each query point q, InPoly returns one of four char's:
!    i : q is strictly interior to P
!    o : q is strictly exterior to P
!    v : q is a vertex of P
!    e : q lies on the relative interior of an edge of P
!
! Initial Author :  Wei Shangguan, 02/2014
! --------------------------------------------------------
use precision

IMPLICIT NONE

   integer, intent(in) :: ponum ! number of points in a polygon
   real(r8), intent(in) :: xp, yp   ! x, y of a point
   real(r8), intent(in) :: xpol(ponum), ypol(ponum)
   logical, intent(out) :: c    ! indicate wheather a soil is in an class

   integer i, i1   ! point index; i1 = i-1 mod n 
   real(r8) x      ! x intersection of e with ray
   integer Rcross  ! number of right edge/ray crossings 
   integer Lcross  ! number of left edge/ray crossings
   character c2

   Rcross = 0
   Lcross = 0
   c2 = ''

! For each edge e=(i-1,i), see if crosses ray. 
   do i = 1, ponum
! First see if q=(0,0) is a vertex. 
      if(( xpol(i) - xp )==0 .AND. ( ypol(i) - yp )==0 )then
           c2 = 'v'
           exit
      end if
      i1 = mod(( i-2 + ponum ), ponum) + 1

! if e "straddles" the x-axis... 
      if( (( ypol(i) - yp ) > 0 ) .NEQV. (( ypol(i1) - yp ) > 0 ) )then
          ! e straddles ray, so compute intersection with ray. 
          x = ( (xpol(i)-xp)*(ypol(i1)-yp) - (xpol(i1)-xp )*(ypol(i)-yp) ) &
              / (ypol(i1)-ypol(i))
          ! crosses ray if strictly positive intersection. 
          if(x > 0)then
             Rcross=Rcross+1
          end if
      end if

! if e straddles the x-axis when reversed... 
      if( (( ypol(i) - yp ) < 0 ) .NEQV. (( ypol(i1) - yp ) < 0 ) )then
    ! e straddles ray, so compute intersection with ray.
          x = ( (xpol(i)-xp)*(ypol(i1)-yp) - (xpol(i1)-xp)*(ypol(i)-yp) ) &
              / (ypol(i1)-ypol(i))
    ! crosses ray if strictly positive intersection. 
          if(x < 0)then 
             Lcross=Lcross+1
          end if
      end if

   end do  
  
    ! q on the edge if left and right cross are not the same parity /
    if(c2=='v')then
       c = .true.
    else if( mod(Rcross,2) .NE. mod(Lcross, 2) )then
       c = .true.
       c2 = 'e'  
    ! q inside iff an odd number of crossings.
    else if( mod(Rcross,2) == 1 )then
       c = .true.
       c2 = 'i'
    else 
       c = .false.
       c2 = 'o'
    end if 

END SUBROUTINE pointinpolygon
