
 SUBROUTINE UrbanPerviousTem (patchtype,lb,deltim, &
                              capr,cnfac,csol,porsl,dkdry,dksatu,&
                              dz_gpersno,z_gpersno,zi_gpersno,&
                              t_gpersno,wice_gpersno,wliq_gpersno,scv_gper,snowdp_gper,&
                              lgper,clgper,sabgper,fsengper,fevpgper,cgper,htvp,&
                              imelt,sm,xmf,fact)

!=======================================================================
! Snow and road temperatures
! o The volumetric heat capacity is calculated as a linear combination
!   in terms of the volumetric fraction of the constituent phases.
! o The thermal conductivity of road soil is computed from
!   the algorithm of Johansen (as reported by Farouki 1981), impervious and perivious from
!   LOOK-UP table and of snow is from the formulation used in SNTHERM (Jordan 1991).
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
! Original author : Yongjiu Dai, 09/15/1999; 08/30/2002; 05/2020
!=======================================================================

  USE precision
  USE GlobalVars
  USE PhysicalConstants
  USE SOIL_thermal_parameters

  IMPLICIT NONE

  INTEGER, intent(in) :: lb        !lower bound of array
  INTEGER, intent(in) :: patchtype !land water TYPE (0=soil,1=urban or built-up,2=wetland,
                                   !3=land ice, 4=deep lake, 5=shallow lake)
  REAL(r8), intent(in) :: deltim   !seconds in a time step [second]
  REAL(r8), intent(in) :: capr     !tuning factor to turn first layer T into surface T
  REAL(r8), intent(in) :: cnfac    !Crank Nicholson factor between 0 and 1

  REAL(r8), intent(in) :: csol  (1:nl_soil) !heat capacity of soil solids [J/(m3 K)]
  REAL(r8), intent(in) :: porsl (1:nl_soil) !soil porosity [-]

  REAL(r8), intent(in) :: dkdry (1:nl_soil) !thermal conductivity of dry soil [W/m-K]
  REAL(r8), intent(in) :: dksatu(1:nl_soil) !thermal conductivity of saturated soil [W/m-K]

  REAL(r8), intent(in) :: dz_gpersno(lb:nl_soil)   !layer thickiness [m]
  REAL(r8), intent(in) :: z_gpersno (lb:nl_soil)   !node depth [m]
  REAL(r8), intent(in) :: zi_gpersno(lb-1:nl_soil) !interface depth [m]

  REAL(r8), intent(in) :: sabgper  !solar radiation absorbed by ground [W/m2]
  REAL(r8), intent(in) :: lgper    !atmospheric infrared (longwave) radiation [W/m2]
  REAL(r8), intent(in) :: clgper   !deriv. of longwave wrt to soil temp [w/m2/k]
  REAL(r8), intent(in) :: fsengper !sensible heat flux from ground [W/m2]
  REAL(r8), intent(in) :: fevpgper !evaporation heat flux from ground [mm/s]
  REAL(r8), intent(in) :: cgper    !deriv. of soil energy flux wrt to soil temp [w/m2/k]
  REAL(r8), intent(in) :: htvp     !latent heat of vapor of water (or sublimation) [j/kg]

  REAL(r8), intent(inout) :: t_gpersno (lb:nl_soil)   !soil temperature [K]
  REAL(r8), intent(inout) :: wice_gpersno(lb:nl_soil) !ice lens [kg/m2]
  REAL(r8), intent(inout) :: wliq_gpersno(lb:nl_soil) !liqui water [kg/m2]
  REAL(r8), intent(inout) :: scv_gper                 !snow cover, water equivalent [mm, kg/m2]
  REAL(r8), intent(inout) :: snowdp_gper              !snow depth [m]

  REAL(r8), intent(out) :: sm                         !rate of snowmelt [kg/(m2 s)]
  REAL(r8), intent(out) :: xmf                        !total latent heat of phase change of ground water
  REAL(r8), intent(out) :: fact(lb:nl_soil)           !used in computing tridiagonal matrix
  INTEGER,  intent(out) :: imelt(lb:nl_soil)          !flag for melting or freezing [-]

!------------------------ local variables ------------------------------
  REAL(r8) cv(lb:nl_soil)     !heat capacity [J/(m2 K)]
  REAL(r8) tk(lb:nl_soil)     !thermal conductivity [W/(m K)]

  REAL(r8) at(lb:nl_soil)     !"a" vector for tridiagonal matrix
  REAL(r8) bt(lb:nl_soil)     !"b" vector for tridiagonal matrix
  REAL(r8) ct(lb:nl_soil)     !"c" vector for tridiagonal matrix
  REAL(r8) rt(lb:nl_soil)     !"r" vector for tridiagonal solution

  REAL(r8) fn (lb:nl_soil)    !heat diffusion through the layer interface [W/m2]
  REAL(r8) fn1(lb:nl_soil)    !heat diffusion through the layer interface [W/m2]
  REAL(r8) dzm                !used in computing tridiagonal matrix
  REAL(r8) dzp                !used in computing tridiagonal matrix

  REAL(r8) t_gpersno_bef(lb:nl_soil) !soil/snow temperature before update
  REAL(r8) hs                 !net energy flux into the surface (w/m2)
  REAL(r8) dhsdt              !d(hs)/dT
  REAL(r8) brr(lb:nl_soil)    !temporay set

  INTEGER i,j

!=======================================================================
! heat capacity
      CALL hCapacity (patchtype,lb,nl_soil,csol,porsl,wice_gpersno,wliq_gpersno,scv_gper,dz_gpersno,cv)

! thermal conductivity
      CALL hConductivity (patchtype,lb,nl_soil,&
                          dkdry,dksatu,porsl,dz_gpersno,z_gpersno,zi_gpersno,&
                          t_gpersno,wice_gpersno,wliq_gpersno,tk)

! net ground heat flux into the surface and its temperature derivative
      hs = sabgper + lgper - (fsengper+fevpgper*htvp)
      dhsdT = - cgper + clgper

      t_gpersno_bef(lb:) = t_gpersno(lb:)

      j       = lb
      fact(j) = deltim / cv(j) &
              * dz_gpersno(j) / (0.5*(z_gpersno(j)-zi_gpersno(j-1)+capr*(z_gpersno(j+1)-zi_gpersno(j-1))))

      DO j = lb + 1, nl_soil
         fact(j) = deltim/cv(j)
      ENDDO

      DO j = lb, nl_soil - 1
        fn(j) = tk(j)*(t_gpersno(j+1)-t_gpersno(j))/(z_gpersno(j+1)-z_gpersno(j))
      ENDDO
      fn(nl_soil) = 0.

! set up vector r and vectors a, b, c that define tridiagonal matrix
      j     = lb
      dzp   = z_gpersno(j+1)-z_gpersno(j)
      at(j) = 0.
      bt(j) = 1+(1.-cnfac)*fact(j)*tk(j)/dzp-fact(j)*dhsdT
      ct(j) =  -(1.-cnfac)*fact(j)*tk(j)/dzp
      rt(j) = t_gpersno(j) + fact(j)*( hs - dhsdT*t_gpersno(j) + cnfac*fn(j) )


      DO j = lb + 1, nl_soil - 1
         dzm   = (z_gpersno(j)-z_gpersno(j-1))
         dzp   = (z_gpersno(j+1)-z_gpersno(j))
         at(j) =   - (1.-cnfac)*fact(j)* tk(j-1)/dzm
         bt(j) = 1.+ (1.-cnfac)*fact(j)*(tk(j)/dzp + tk(j-1)/dzm)
         ct(j) =   - (1.-cnfac)*fact(j)* tk(j)/dzp
         rt(j) = t_gpersno(j) + cnfac*fact(j)*( fn(j) - fn(j-1) )
      ENDDO

      j     =  nl_soil
      dzm   = (z_gpersno(j)-z_gpersno(j-1))
      at(j) =   - (1.-cnfac)*fact(j)*tk(j-1)/dzm
      bt(j) = 1.+ (1.-cnfac)*fact(j)*tk(j-1)/dzm
      ct(j) = 0.
      rt(j) = t_gpersno(j) - cnfac*fact(j)*fn(j-1)

! solve for t_gpersno
      i = size(at)
      CALL tridia (i ,at ,bt ,ct ,rt ,t_gpersno)

!=======================================================================
! melting or freezing
!=======================================================================

      DO j = lb, nl_soil - 1
         fn1(j) = tk(j)*(t_gpersno(j+1)-t_gpersno(j))/(z_gpersno(j+1)-z_gpersno(j))
      ENDDO
      fn1(nl_soil) = 0.

      j = lb
      brr(j) = cnfac*fn(j) + (1.-cnfac)*fn1(j)

      DO j = lb + 1, nl_soil
         brr(j) = cnfac*(fn(j)-fn(j-1)) + (1.-cnfac)*(fn1(j)-fn1(j-1))
      ENDDO

      CALL meltf (lb,nl_soil,deltim, &
                  fact(lb:),brr(lb:),hs,dhsdT, &
                  t_gpersno_bef(lb:),t_gpersno(lb:),wliq_gpersno(lb:),wice_gpersno(lb:),imelt(lb:), &
                  scv_gper,snowdp_gper,sm,xmf)

 END SUBROUTINE UrbanPerviousTem
! ---------- EOP ------------
