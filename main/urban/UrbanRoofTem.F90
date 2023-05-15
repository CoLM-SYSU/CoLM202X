
 SUBROUTINE UrbanRoofTem (lb,deltim,capr,cnfac,&
                          cv_roof,tk_roof,dz_roofsno,z_roofsno,zi_roofsno,&
                          t_roofsno,wice_roofsno,wliq_roofsno,scv_roof,snowdp_roof,&
                          troof_inner,lroof,clroof,sabroof,fsenroof,fevproof,croof,htvp,&
                          imelt_roof,sm_roof,xmf_roof,fact,tkdz_roof)

!=======================================================================
! Snow and roof temperatures
! o The volumetric heat capacity is calculated as a linear combination
!   in terms of the volumetric fraction of the constituent phases.
! o The thermal conductivity of roof is given by LOOK-UP table, and of snow is from
!   the formulation used in SNTHERM (Jordan 1991).
! o Boundary conditions:
!   F = Rnet - Hg - LEg (top),
!   For urban sunwall, shadewall, and roof columns, there is a non-zero heat flux across
!   the bottom "building inner surface" layer and the equations are derived assuming
!   a prescribed or adjusted internal building temperature.
!   T = T_roof_inner (at the roof inner surface).
! o Roof / snow temperature is predicted from heat conduction
!   in N roof layers and up to 5 snow layers.
!   The thermal conductivities at the interfaces between two neighbor layers
!   (j, j+1) are derived from an assumption that the flux across the interface
!   is equal to that from the node j to the interface and the flux from the
!   interface to the node j+1. The equation is solved using the Crank-Nicholson
!   method and resulted in a tridiagonal system equation.
!
! Phase change (see meltf.F90)
!
! Original author : Yongjiu Dai, 05/2020
!=======================================================================

  USE precision
  USE GlobalVars
  USE PhysicalConstants

  IMPLICIT NONE

  INTEGER , intent(in) :: lb       !lower bound of array
  REAL(r8), intent(in) :: deltim   !seconds in a time step [second]
  REAL(r8), intent(in) :: capr     !tuning factor to turn first layer T into surface T
  REAL(r8), intent(in) :: cnfac    !Crank Nicholson factor between 0 and 1

  REAL(r8), intent(in) :: cv_roof(1:nl_roof)       !heat capacity of urban roof [J/m3/K]
  REAL(r8), intent(in) :: tk_roof(1:nl_roof)       !thermal conductivity of urban roof [W/m/K]

  REAL(r8), intent(in) :: dz_roofsno(lb:nl_roof)   !layer thickiness [m]
  REAL(r8), intent(in) :: z_roofsno (lb:nl_roof)   !node depth [m]
  REAL(r8), intent(in) :: zi_roofsno(lb-1:nl_roof) !interface depth [m]

  REAL(r8), intent(in) :: troof_inner !temperature at the roof inner surface [K]
  REAL(r8), intent(in) :: lroof    !atmospheric infrared (longwave) radiation [W/m2]
  REAL(r8), intent(in) :: clroof   !atmospheric infrared (longwave) radiation [W/m2]
  REAL(r8), intent(in) :: sabroof  !solar radiation absorbed by roof [W/m2]
  REAL(r8), intent(in) :: fsenroof !sensible heat flux from roof [W/m2]
  REAL(r8), intent(in) :: fevproof !evaporation heat flux from roof [mm/s]
  REAL(r8), intent(in) :: croof    !deriv. of roof energy flux wrt to roof temp [w/m2/k]
  REAL(r8), intent(in) :: htvp     !latent heat of vapor of water (or sublimation) [j/kg]

  REAL(r8), intent(inout) :: t_roofsno   (lb:nl_roof) !roof layers' temperature [K]
  REAL(r8), intent(inout) :: wice_roofsno(lb:nl_roof) !ice lens [kg/m2]
  REAL(r8), intent(inout) :: wliq_roofsno(lb:nl_roof) !liqui water [kg/m2]
  REAL(r8), intent(inout) :: scv_roof                 !snow cover, water equivalent [mm, kg/m2]
  REAL(r8), intent(inout) :: snowdp_roof              !snow depth [m]

  REAL(r8), intent(out) :: sm_roof                    !rate of snowmelt [kg/(m2 s)]
  REAL(r8), intent(out) :: xmf_roof                   !total latent heat of phase change of roof residual water
  REAL(r8), intent(out) :: fact(lb:nl_roof)           !used in computing tridiagonal matrix
  REAL(r8), intent(out) :: tkdz_roof                  !heat diffusion with inner room space
  INTEGER , intent(out) :: imelt_roof(lb:nl_roof)     !flag for melting or freezing [-]

!------------------------ local variables ------------------------------
  REAL(r8) cv (lb:nl_roof)  !heat capacity [J/(m2 K)]
  REAL(r8) thk(lb:nl_roof)  !thermal conductivity of layer
  REAL(r8) tk (lb:nl_roof)  !thermal conductivity [W/(m K)]

  REAL(r8) at (lb:nl_roof)  !"a" vector for tridiagonal matrix
  REAL(r8) bt (lb:nl_roof)  !"b" vector for tridiagonal matrix
  REAL(r8) ct (lb:nl_roof)  !"c" vector for tridiagonal matrix
  REAL(r8) rt (lb:nl_roof)  !"r" vector for tridiagonal solution

  REAL(r8) fn (lb:nl_roof)  !heat diffusion through the layer interface [W/m2]
  REAL(r8) fn1(lb:nl_roof)  !heat diffusion through the layer interface [W/m2]
  REAL(r8) dzm              !used in computing tridiagonal matrix
  REAL(r8) dzp              !used in computing tridiagonal matrix

  REAL(r8) t_roofsno_bef(lb:nl_roof) !roof/snow temperature before update
  REAL(r8) hs               !net energy flux into the surface (w/m2)
  REAL(r8) dhsdt            !d(hs)/dT
  REAL(r8) brr(lb:nl_roof)  !temporay set
  REAL(r8) bw               !snow density [kg/m3]

  INTEGER i,j

!=======================================================================

      wice_roofsno(2:) = 0.0 !ice lens [kg/m2]
      wliq_roofsno(2:) = 0.0 !liquid water [kg/m2]

! heat capacity
      IF (lb <= 0) THEN
         DO j = lb, 0
            cv(j) = max(1.0e-6_r8,(cpliq*wliq_roofsno(j) + cpice*wice_roofsno(j)))
         ENDDO
      ENDIF

      cv(1:) = cv_roof(1:)*dz_roofsno(1:)
      IF (lb == 1 .and. scv_roof > 0.0) THEN
         cv(1) = cv(1) + cpice*scv_roof
      ELSE
         !ponding water
         cv(1) = cv(1) + cpliq*wliq_roofsno(1) + cpice*wice_roofsno(1)
      ENDIF

! thermal conductivity
      ! Thermal conductivity of snow, which from Yen (1980)
      IF (lb <= 0) THEN
         DO j = lb, 0
         bw = (wice_roofsno(j)+wliq_roofsno(j))/(dz_roofsno(j))
         thk(j) = tkair + (7.75e-5_r8 *bw + 1.105e-6_r8*bw*bw)*(tkice-tkair) ! Yen, 1980
         !thk(j) = 0.024 - 1.23e-4_r8*bw + 2.5e-6_r8*bw*bw ! Calonne et al., 2011
         ENDDO
      ENDIF

    ! thermal conductivity at the layer interface
      thk(1:) = tk_roof(1:)
      IF (lb <= 0) THEN
         DO j = lb, 0
         tk(j) = thk(j)*thk(j+1)*(z_roofsno(j+1)-z_roofsno(j)) &
               /(thk(j)*(z_roofsno(j+1)-zi_roofsno(j))+thk(j+1)*(zi_roofsno(j)-z_roofsno(j)))
         ENDDO
      ENDIF

      DO j = 1, nl_roof-1
         tk(j) = thk(j)*thk(j+1)*(z_roofsno(j+1)-z_roofsno(j)) &
               /(thk(j)*(z_roofsno(j+1)-zi_roofsno(j))+thk(j+1)*(zi_roofsno(j)-z_roofsno(j)))
      ENDDO
      tk(nl_roof) = thk(nl_roof)

!???
!     ! update thermal conductivity of the ponding water
!     zh2osfc=1.0e-3*(0.5*h2osfc(c)) !convert to [m] from [mm]
!     tk(1)= tkwat*thk(1)*(z(1)+zh2osfc) &
!             /(tkwat*z(1)+thk(1)*zh2osfc)

! net ground heat flux into the roof surface and its temperature derivative
      hs = sabroof + lroof - (fsenroof+fevproof*htvp)
      dhsdT = - croof + clroof

      t_roofsno_bef(lb:) = t_roofsno(lb:)

      j       = lb
      fact(j) = deltim / cv(j) * dz_roofsno(j) &
              / (0.5*(z_roofsno(j)-zi_roofsno(j-1)+capr*(z_roofsno(j+1)-zi_roofsno(j-1))))

      DO j = lb + 1, nl_roof
         fact(j) = deltim/cv(j)
      ENDDO

      DO j = lb, nl_roof - 1
        fn(j) = tk(j)*(t_roofsno(j+1)-t_roofsno(j))/(z_roofsno(j+1)-z_roofsno(j))
      ENDDO

      j     = nl_roof
      fn(j) = tk(j)*(troof_inner - cnfac*t_roofsno(j))/(zi_roofsno(j)-z_roofsno(j))
      tkdz_roof = tk(j)/(zi_roofsno(j)-z_roofsno(j))

! set up vector r and vectors a, b, c that define tridiagonal matrix
      j     = lb
      dzp   = z_roofsno(j+1)-z_roofsno(j)
      at(j) = 0.
      bt(j) = 1+(1.-cnfac)*fact(j)*tk(j)/dzp-fact(j)*dhsdT
      ct(j) =  -(1.-cnfac)*fact(j)*tk(j)/dzp
      rt(j) = t_roofsno(j) + fact(j)*( hs - dhsdT*t_roofsno(j) + cnfac*fn(j) )

      DO j = lb + 1, nl_roof - 1
         dzm   = (z_roofsno(j)-z_roofsno(j-1))
         dzp   = (z_roofsno(j+1)-z_roofsno(j))
         at(j) =   - (1.-cnfac)*fact(j)* tk(j-1)/dzm
         bt(j) = 1.+ (1.-cnfac)*fact(j)*(tk(j)/dzp + tk(j-1)/dzm)
         ct(j) =   - (1.-cnfac)*fact(j)* tk(j)/dzp
         rt(j) = t_roofsno(j) + cnfac*fact(j)*( fn(j) - fn(j-1) )
      ENDDO

      j     =  nl_roof
      dzm   = (z_roofsno(j)-z_roofsno(j-1))
      dzp   = (zi_roofsno(j)-z_roofsno(j))
      at(j) =   - (1.-cnfac)*fact(j)*tk(j-1)/dzm
      bt(j) = 1.+ (1.-cnfac)*fact(j)*(tk(j-1)/dzm+tk(j)/dzp)
      ct(j) = 0.
      rt(j) = t_roofsno(j) + fact(j)*(fn(j) - cnfac*fn(j-1))

! solve for t_roofsno
      i = size(at)
      CALL tridia (i ,at ,bt ,ct ,rt ,t_roofsno)

!=======================================================================
! melting or freezing
!=======================================================================

      DO j = lb, nl_roof - 1
         fn1(j) = tk(j)*(t_roofsno(j+1)-t_roofsno(j))/(z_roofsno(j+1)-z_roofsno(j))
      ENDDO

      j = nl_roof
      fn1(j) = tk(j)*(troof_inner - cnfac*t_roofsno(j))/(zi_roofsno(j)-z_roofsno(j))

      j = lb
      brr(j) = cnfac*fn(j) + (1.-cnfac)*fn1(j)

      DO j = lb + 1, nl_roof
         brr(j) = cnfac*(fn(j)-fn(j-1)) + (1.-cnfac)*(fn1(j)-fn1(j-1))
      ENDDO

      CALL meltf_urban (lb,1,deltim, &
                  fact(lb:1),brr(lb:1),hs,dhsdT, &
                  t_roofsno_bef(lb:1),t_roofsno(lb:1),wliq_roofsno(lb:1),wice_roofsno(lb:1),imelt_roof(lb:1), &
                  scv_roof,snowdp_roof,sm_roof,xmf_roof)

 END SUBROUTINE UrbanRoofTem
! ---------- EOP ------------
