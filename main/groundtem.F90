
 subroutine groundtem (patchtype,lb,nl_soil,deltim, &
                       capr,cnfac,csol,porsl,dkdry,dksatu, &
                       sigf,dz_soisno,z_soisno,zi_soisno,&
                       t_soisno,wice_soisno,wliq_soisno,scv,snowdp, &
                       frl,dlrad,sabg,fseng,fevpg,cgrnd,htvp,emg, &
                       imelt,sm,xmf,fact,psi0,bsw)
                       !TODO: not used, psi0, bsw

!=======================================================================
! Snow and soil temperatures
! o The volumetric heat capacity is calculated as a linear combination
!   in terms of the volumetric fraction of the constituent phases.
! o The thermal conductivity of soil is computed from
!   the algorithm of Johansen (as reported by Farouki 1981), and of snow is from
!   the formulation used in SNTHERM (Jordan 1991).
! o Boundary conditions:
!   F = Rnet - Hg - LEg (top),  F= 0 (base of the soil column).
! o Soil / snow temperature is predicted from heat conduction
!   in 10 soil layers and up to 5 snow layers.
!   The thermal conductivities at the interfaces between two neighbor layers
!   (j, j+1) are derived from an assumption that the flux across the interface
!   is equal to that from the node j to the interface and the flux from the
!   interface to the node j+1. The equation is solved using the Crank-Nicholson
!   method and resulted in a tridiagonal system equation.
!
! Phase change (see meltf.F90)
! 
! Original author : Yongjiu Dai, 09/15/1999; 08/30/2002
!=======================================================================

  use precision
  use PhysicalConstants, only: stefnc
  use SOIL_thermal_parameters

  implicit none

  integer, INTENT(in) :: lb           !lower bound of array
  integer, INTENT(in) :: nl_soil      !upper bound of array
  integer, INTENT(in) :: patchtype    !land water type (0=soil,1=urban or built-up,2=wetland,
                                      !3=land ice, 4=deep lake, 5=shallow lake)
  real(r8), INTENT(in) :: deltim      !seconds in a time step [second]
  real(r8), INTENT(in) :: capr        !tuning factor to turn first layer T into surface T
  real(r8), INTENT(in) :: cnfac       !Crank Nicholson factor between 0 and 1

  real(r8), INTENT(in) :: csol (1:nl_soil) !heat capacity of soil solids [J/(m3 K)]
  real(r8), INTENT(in) :: porsl(1:nl_soil) !soil porosity [-]
  real(r8), INTENT(in) :: psi0 (1:nl_soil) !soil water suction, negative potential [m]
  real(r8), INTENT(in) :: bsw  (1:nl_soil) !clapp and hornbereger "b" parameter [-]

  real(r8), INTENT(in) :: dkdry(1:nl_soil) !thermal conductivity of dry soil [W/m-K]
  real(r8), INTENT(in) :: dksatu(1:nl_soil)!thermal conductivity of saturated soil [W/m-K]

  real(r8), INTENT(in) :: sigf     !fraction of veg cover, excluding snow-covered veg [-]
  real(r8), INTENT(in) :: dz_soisno(lb:nl_soil)   !layer thickiness [m]
  real(r8), INTENT(in) :: z_soisno (lb:nl_soil)   !node depth [m]
  real(r8), INTENT(in) :: zi_soisno(lb-1:nl_soil) !interface depth [m]

  real(r8), INTENT(in) :: sabg     !solar radiation absorbed by ground [W/m2]
  real(r8), INTENT(in) :: frl      !atmospheric infrared (longwave) radiation [W/m2]
  real(r8), INTENT(in) :: dlrad    !downward longwave radiation blow the canopy [W/m2]
  real(r8), INTENT(in) :: fseng    !sensible heat flux from ground [W/m2]
  real(r8), INTENT(in) :: fevpg    !evaporation heat flux from ground [mm/s]
  real(r8), INTENT(in) :: cgrnd    !deriv. of soil energy flux wrt to soil temp [w/m2/k]
  real(r8), INTENT(in) :: htvp     !latent heat of vapor of water (or sublimation) [j/kg]
  real(r8), INTENT(in) :: emg      !ground emissivity (0.97 for snow,

  real(r8), INTENT(inout) :: t_soisno (lb:nl_soil)   !soil temperature [K]
  real(r8), INTENT(inout) :: wice_soisno(lb:nl_soil) !ice lens [kg/m2]
  real(r8), INTENT(inout) :: wliq_soisno(lb:nl_soil) !liqui water [kg/m2]
  real(r8), INTENT(inout) :: scv      !snow cover, water equivalent [mm, kg/m2]
  real(r8), INTENT(inout) :: snowdp   !snow depth [m]

  real(r8), INTENT(out) :: sm         !rate of snowmelt [kg/(m2 s)]
  real(r8), INTENT(out) :: xmf        !total latent heat of phase change of ground water
  real(r8), INTENT(out) :: fact(lb:nl_soil)  !used in computing tridiagonal matrix
  integer,  INTENT(out) :: imelt(lb:nl_soil) !flag for melting or freezing [-]

!------------------------ local variables ------------------------------
  real(r8) cv(lb:nl_soil)     !heat capacity [J/(m2 K)]
  real(r8) tk(lb:nl_soil)     !thermal conductivity [W/(m K)]

  real(r8) at(lb:nl_soil)     !"a" vector for tridiagonal matrix
  real(r8) bt(lb:nl_soil)     !"b" vector for tridiagonal matrix
  real(r8) ct(lb:nl_soil)     !"c" vector for tridiagonal matrix
  real(r8) rt(lb:nl_soil)     !"r" vector for tridiagonal solution

  real(r8) fn  (lb:nl_soil)   !heat diffusion through the layer interface [W/m2]
  real(r8) fn1 (lb:nl_soil)   !heat diffusion through the layer interface [W/m2]
  real(r8) dzm                !used in computing tridiagonal matrix
  real(r8) dzp                !used in computing tridiagonal matrix

  real(r8) t_soisno_bef(lb:nl_soil) !soil/snow temperature before update
  real(r8) hs                 !net energy flux into the surface (w/m2)
  real(r8) dhsdt              !d(hs)/dT
  real(r8) brr(lb:nl_soil)    !temporay set

  integer i,j

!=======================================================================
! heat capacity 
      call hCapacity (patchtype,lb,nl_soil,csol,porsl,wice_soisno,wliq_soisno,scv,dz_soisno,cv)

! thermal conductivity
      if(zi_soisno(0) < 0.)then
         print*,'[groundtem],zi_soisno(0)',zi_soisno(0)
         stop
      endif
      
      call hConductivity (patchtype,lb,nl_soil,&
                          dkdry,dksatu,porsl,dz_soisno,z_soisno,zi_soisno,&
                          t_soisno,wice_soisno,wliq_soisno,tk)

! net ground heat flux into the surface and its temperature derivative
      hs = sabg + dlrad*emg &
! 08/19/2021, yuan: remove sigf, LAI->100% cover
         !+ (1.-sigf)*emg*frl - emg*stefnc*t_soisno(lb)**4 &
         - emg*stefnc*t_soisno(lb)**4 &
         - (fseng+fevpg*htvp) 

      dhsdT = - cgrnd - 4.*emg * stefnc * t_soisno(lb)**3
      t_soisno_bef(lb:) = t_soisno(lb:)

      j       = lb
      fact(j) = deltim / cv(j) &
              * dz_soisno(j) / (0.5*(z_soisno(j)-zi_soisno(j-1)+capr*(z_soisno(j+1)-zi_soisno(j-1))))

      do j = lb + 1, nl_soil
         fact(j) = deltim/cv(j)
      enddo

      do j = lb, nl_soil - 1
        fn(j) = tk(j)*(t_soisno(j+1)-t_soisno(j))/(z_soisno(j+1)-z_soisno(j))
      enddo
      fn(nl_soil) = 0.

! set up vector r and vectors a, b, c that define tridiagonal matrix
      j     = lb
      dzp   = z_soisno(j+1)-z_soisno(j)
      at(j) = 0.
      bt(j) = 1+(1.-cnfac)*fact(j)*tk(j)/dzp-fact(j)*dhsdT
      ct(j) =  -(1.-cnfac)*fact(j)*tk(j)/dzp
      rt(j) = t_soisno(j) + fact(j)*( hs - dhsdT*t_soisno(j) + cnfac*fn(j) )


      do j = lb + 1, nl_soil - 1
         dzm   = (z_soisno(j)-z_soisno(j-1))
         dzp   = (z_soisno(j+1)-z_soisno(j))
         at(j) =   - (1.-cnfac)*fact(j)* tk(j-1)/dzm
         bt(j) = 1.+ (1.-cnfac)*fact(j)*(tk(j)/dzp + tk(j-1)/dzm)
         ct(j) =   - (1.-cnfac)*fact(j)* tk(j)/dzp
         rt(j) = t_soisno(j) + cnfac*fact(j)*( fn(j) - fn(j-1) )
      end do

      j     =  nl_soil
      dzm   = (z_soisno(j)-z_soisno(j-1))
      at(j) =   - (1.-cnfac)*fact(j)*tk(j-1)/dzm
      bt(j) = 1.+ (1.-cnfac)*fact(j)*tk(j-1)/dzm
      ct(j) = 0.
      rt(j) = t_soisno(j) - cnfac*fact(j)*fn(j-1)

! solve for t_soisno
      i = size(at)
      call tridia (i ,at ,bt ,ct ,rt ,t_soisno) 

!=======================================================================
! melting or freezing 
!=======================================================================

      do j = lb, nl_soil - 1
         fn1(j) = tk(j)*(t_soisno(j+1)-t_soisno(j))/(z_soisno(j+1)-z_soisno(j))
      enddo
      fn1(nl_soil) = 0.

      j = lb
      brr(j) = cnfac*fn(j) + (1.-cnfac)*fn1(j)

      do j = lb + 1, nl_soil
         brr(j) = cnfac*(fn(j)-fn(j-1)) + (1.-cnfac)*(fn1(j)-fn1(j-1))
      enddo

      call meltf (lb,nl_soil,deltim, &
                  fact(lb:),brr(lb:),hs,dhsdT, &
                  t_soisno_bef(lb:),t_soisno(lb:),wliq_soisno(lb:),wice_soisno(lb:),imelt(lb:), &
                  scv,snowdp,sm,xmf)

!-----------------------------------------------------------------------

 end subroutine groundtem
