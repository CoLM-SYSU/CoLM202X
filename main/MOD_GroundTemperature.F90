#include <define.h>

MODULE MOD_GroundTemperature

!-----------------------------------------------------------------------
   USE MOD_Precision
   IMPLICIT NONE
   SAVE

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: GroundTemperature


!-----------------------------------------------------------------------

   CONTAINS

!-----------------------------------------------------------------------


   SUBROUTINE GroundTemperature (patchtype,lb,nl_soil,deltim,&
                         capr,cnfac,vf_quartz,vf_gravels,vf_om,vf_sand,wf_gravels,wf_sand,&
                         porsl,psi0,&
#ifdef Campbell_SOIL_MODEL
                         bsw,&
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
                         theta_r, alpha_vgm, n_vgm, L_vgm,&
                         sc_vgm , fc_vgm,&
#endif
                         csol,k_solids,dksatu,dksatf,dkdry,&
                         BA_alpha,BA_beta,&
                         sigf,dz_soisno,z_soisno,zi_soisno,&
                         t_soisno,wice_soisno,wliq_soisno,scv,snowdp,&
                         frl,dlrad,sabg,sabg_lyr,fseng,fevpg,cgrnd,htvp,emg,&
                         imelt,snofrz,sm,xmf,fact,pg_rain,pg_snow,t_precip)

!=======================================================================
! Snow and soil temperatures
! o The volumetric heat capacity is calculated as a linear combination
!   in terms of the volumetric fraction of the constituent phases.
! o The thermal conductivity of soil is computed from
!   the algorithm of Johansen (as reported by Farouki 1981), and of snow is from
!   the formulation used in SNTHERM (Jordan 1991).
! o Boundary conditions:
!   F = Rnet - Hg - LEg + Hpr(top),  F= 0 (base of the soil column).
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
! Original author : Yongjiu Dai, 09/15/1999; 08/30/2002; 05/2018
!
! REVISIONS:
! Nan Wei,  07/2017: interaction btw prec and land surface
! Nan Wei,  01/2019: USE the new version of soil thermal parameters to calculate soil temperature
! Hua Yuan, 01/2023: modified ground heat flux, temperature and meltf
!                    calculation for SNICAR model
!=======================================================================

   USE MOD_Precision
   USE MOD_Const_Physical, only: stefnc,denh2o,denice,tfrz,cpice,cpliq,tkwat,tkice,tkair
   USE MOD_Namelist, only: DEF_USE_SNICAR
   USE MOD_PhaseChange
   USE MOD_SoilThermalParameters
   USE MOD_Utils

   IMPLICIT NONE

   integer, intent(in) :: lb           !lower bound of array
   integer, intent(in) :: nl_soil      !upper bound of array
   integer, intent(in) :: patchtype    !land patch type (0=soil,1=urban or built-up,2=wetland,
                                       !3=land ice, 4=deep lake, 5=shallow lake)
   real(r8), intent(in) :: deltim      !seconds in a time step [second]
   real(r8), intent(in) :: capr        !tuning factor to turn first layer T into surface T
   real(r8), intent(in) :: cnfac       !Crank Nicholson factor between 0 and 1

   real(r8), intent(in) :: vf_quartz (1:nl_soil) ! volumetric fraction of quartz within mineral soil
   real(r8), intent(in) :: vf_gravels(1:nl_soil) ! volumetric fraction of gravels
   real(r8), intent(in) :: vf_om     (1:nl_soil) ! volumetric fraction of organic matter
   real(r8), intent(in) :: vf_sand   (1:nl_soil) ! volumetric fraction of sand
   real(r8), intent(in) :: wf_gravels(1:nl_soil) ! gravimetric fraction of gravels
   real(r8), intent(in) :: wf_sand   (1:nl_soil) ! gravimetric fraction of sand

   real(r8), intent(in) :: porsl(1:nl_soil) ! soil porosity [-]
   real(r8), intent(in) :: psi0 (1:nl_soil) ! soil water suction, negative potential [mm]
#ifdef Campbell_SOIL_MODEL
   real(r8), intent(in) :: bsw(1:nl_soil)   ! clapp and hornbereger "b" parameter [-]
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
   real(r8), intent(in) :: theta_r  (1:nl_soil), &
                           alpha_vgm(1:nl_soil), &
                           n_vgm    (1:nl_soil), &
                           L_vgm    (1:nl_soil), &
                           sc_vgm   (1:nl_soil), &
                           fc_vgm   (1:nl_soil)
#endif
   real(r8), intent(in) :: csol(1:nl_soil)   ! heat capacity of soil solids [J/(m3 K)]
   real(r8), intent(in) :: k_solids(1:nl_soil) ! thermal conductivity of minerals soil [W/m-K]
   real(r8), intent(in) :: dksatu(1:nl_soil) ! thermal conductivity of saturated unfrozen soil [W/m-K]
   real(r8), intent(in) :: dksatf(1:nl_soil) ! thermal conductivity of saturated frozen soil [W/m-K]
   real(r8), intent(in) :: dkdry(1:nl_soil)  ! thermal conductivity of dry soil [W/m-K]
   real(r8), intent(in) :: BA_alpha(1:nl_soil) ! alpha in Balland and Arp(2005) thermal conductivity scheme
   real(r8), intent(in) :: BA_beta(1:nl_soil)  ! beta in Balland and Arp(2005) thermal conductivity scheme

   real(r8), intent(in) :: sigf     !fraction of veg cover, excluding snow-covered veg [-]
   real(r8), intent(in) :: dz_soisno(lb:nl_soil)   !layer thickiness [m]
   real(r8), intent(in) :: z_soisno (lb:nl_soil)   !node depth [m]
   real(r8), intent(in) :: zi_soisno(lb-1:nl_soil) !interface depth [m]

   real(r8), intent(in) :: sabg_lyr(lb:1)          !snow layer absorption [W/m-2]

   real(r8), intent(in) :: sabg     !solar radiation absorbed by ground [W/m2]
   real(r8), intent(in) :: frl      !atmospheric infrared (longwave) radiation [W/m2]
   real(r8), intent(in) :: dlrad    !downward longwave radiation blow the canopy [W/m2]
   real(r8), intent(in) :: fseng    !sensible heat flux from ground [W/m2]
   real(r8), intent(in) :: fevpg    !evaporation heat flux from ground [mm/s]
   real(r8), intent(in) :: cgrnd    !deriv. of soil energy flux wrt to soil temp [w/m2/k]
   real(r8), intent(in) :: htvp     !latent heat of vapor of water (or sublimation) [j/kg]
   real(r8), intent(in) :: emg      !ground emissivity (0.97 for snow,
   real(r8), intent(in) :: pg_rain  ! rainfall onto ground including canopy runoff [kg/(m2 s)]
   real(r8), intent(in) :: pg_snow  ! snowfall onto ground including canopy runoff [kg/(m2 s)]
   real(r8), intent(in) :: t_precip ! snowfall/rainfall temperature [kelvin]

   real(r8), intent(inout) :: t_soisno (lb:nl_soil)   !soil temperature [K]
   real(r8), intent(inout) :: wice_soisno(lb:nl_soil) !ice lens [kg/m2]
   real(r8), intent(inout) :: wliq_soisno(lb:nl_soil) !liqui water [kg/m2]
   real(r8), intent(inout) :: scv      !snow cover, water equivalent [mm, kg/m2]
   real(r8), intent(inout) :: snowdp   !snow depth [m]

   real(r8), intent(out) :: sm       !rate of snowmelt [kg/(m2 s)]
   real(r8), intent(out) :: xmf      !total latent heat of phase change of ground water
   real(r8), intent(out) :: fact(lb:nl_soil) !used in computing tridiagonal matrix
   integer,  intent(out) :: imelt(lb:nl_soil)!flag for melting or freezing [-]

   real(r8), intent(out) :: snofrz(lb:0)      !snow freezing rate (lyr) [kg m-2 s-1]

!------------------------ local variables ------------------------------
   real(r8) cv(lb:nl_soil)     ! heat capacity [J/(m2 K)]
   real(r8) tk(lb:nl_soil)     ! thermal conductivity [W/(m K)]
   real(r8) hcap(1:nl_soil)    ! J/(m3 K)
   real(r8) thk(lb:nl_soil)     ! W/(m K)

   real(r8) at(lb:nl_soil)     !"a" vector for tridiagonal matrix
   real(r8) bt(lb:nl_soil)     !"b" vector for tridiagonal matrix
   real(r8) ct(lb:nl_soil)     !"c" vector for tridiagonal matrix
   real(r8) rt(lb:nl_soil)     !"r" vector for tridiagonal solution

   real(r8) fn  (lb:nl_soil)   !heat diffusion through the layer interface [W/m2]
   real(r8) fn1 (lb:nl_soil)   !heat diffusion through the layer interface [W/m2]
   real(r8) dzm                !used in computing tridiagonal matrix
   real(r8) dzp                !used in computing tridiagonal matrix

   real(r8) t_soisno_bef(lb:nl_soil) !soil/snow temperature before update
   real(r8) wice_soisno_bef(lb:0)    !ice lens [kg/m2]
   real(r8) hs                 !net energy flux into the surface (w/m2)
   real(r8) dhsdt              !d(hs)/dT
   real(r8) brr(lb:nl_soil)    !temporay set
   real(r8) vf_water(1:nl_soil) ! volumetric fraction liquid water within soil
   real(r8) vf_ice(1:nl_soil)   ! volumetric fraction ice len within soil
   real(r8) rhosnow  ! partitial density of water (ice + liquid)
   integer i,j

!=======================================================================
! soil ground and wetland heat capacity
   DO i = 1, nl_soil
      vf_water(i) = wliq_soisno(i)/(dz_soisno(i)*denh2o)
      vf_ice(i) = wice_soisno(i)/(dz_soisno(i)*denice)
      CALL soil_hcap_cond(vf_gravels(i),vf_om(i),vf_sand(i),porsl(i),&
                          wf_gravels(i),wf_sand(i),k_solids(i),&
                          csol(i),dkdry(i),dksatu(i),dksatf(i),&
                          BA_alpha(i),BA_beta(i),&
                          t_soisno(i),vf_water(i),vf_ice(i),hcap(i),thk(i))
      cv(i) = hcap(i)*dz_soisno(i)
   ENDDO
   IF(lb==1 .and. scv>0.) cv(1) = cv(1) + cpice*scv

! Snow heat capacity
   IF(lb <= 0)THEN
      cv(:0) = cpliq*wliq_soisno(:0) + cpice*wice_soisno(:0)
   ENDIF

! Snow thermal conductivity
   IF(lb <= 0)THEN
      DO i = lb, 0
      rhosnow = (wice_soisno(i)+wliq_soisno(i))/dz_soisno(i)

      ! presently option [1] is the default option
      ! [1] Jordan (1991) pp. 18
      thk(i) = tkair+(7.75e-5*rhosnow+1.105e-6*rhosnow*rhosnow)*(tkice-tkair)

      ! [2] Sturm et al (1997)
      ! thk(i) = 0.0138 + 1.01e-3*rhosnow + 3.233e-6*rhosnow**2
      ! [3] Ostin and Andersson presented in Sturm et al., (1997)
      ! thk(i) = -0.871e-2 + 0.439e-3*rhosnow + 1.05e-6*rhosnow**2
      ! [4] Jansson(1901) presented in Sturm et al. (1997)
      ! thk(i) = 0.0293 + 0.7953e-3*rhosnow + 1.512e-12*rhosnow**2
      ! [5] Douville et al., (1995)
      ! thk(i) = 2.2*(rhosnow/denice)**1.88
      ! [6] van Dusen (1992) presented in Sturm et al. (1997)
      ! thk(i) = 0.021 + 0.42e-3*rhosnow + 0.22e-6*rhosnow**2

      ENDDO
   ENDIF

! Thermal conductivity at the layer interface
   DO i = lb, nl_soil-1

! the following consideration is try to avoid the snow conductivity
! to be dominant in the thermal conductivity of the interface.
! Because when the distance of bottom snow node to the interfacee
! is larger than that of interface to top soil node,
! the snow thermal conductivity will be dominant, and the result is that
! lees heat tranfer between snow and soil
      IF((i==0) .and. (z_soisno(i+1)-zi_soisno(i)<zi_soisno(i)-z_soisno(i)))THEN
         tk(i) = 2.*thk(i)*thk(i+1)/(thk(i)+thk(i+1))
         tk(i) = max(0.5*thk(i+1),tk(i))
      ELSE
         tk(i) = thk(i)*thk(i+1)*(z_soisno(i+1)-z_soisno(i)) &
               /(thk(i)*(z_soisno(i+1)-zi_soisno(i))+thk(i+1)*(zi_soisno(i)-z_soisno(i)))
      ENDIF
   ENDDO
   tk(nl_soil) = 0.

! net ground heat flux into the surface and its temperature derivative
   IF (DEF_USE_SNICAR) THEN
      hs = sabg_lyr(lb) + dlrad*emg &
         ! 08/19/2021, yuan: NOTE! removed sigf, LAI->100% cover
         - emg*stefnc*t_soisno(lb)**4 &
         - (fseng+fevpg*htvp) + cpliq * pg_rain * (t_precip - t_soisno(lb)) &
         + cpice * pg_snow * (t_precip - t_soisno(lb))
   ELSE
      hs = sabg + dlrad*emg &
         ! 08/19/2021, yuan: NOTE! removed sigf, LAI->100% cover
         - emg*stefnc*t_soisno(lb)**4 &
         - (fseng+fevpg*htvp) + cpliq * pg_rain * (t_precip - t_soisno(lb)) &
         + cpice * pg_snow * (t_precip - t_soisno(lb))
   ENDIF

   dhsdT = - cgrnd - 4.*emg * stefnc * t_soisno(lb)**3 - cpliq * pg_rain - cpice * pg_snow
   t_soisno_bef(lb:) = t_soisno(lb:)

   j       = lb
   fact(j) = deltim / cv(j) &
           * dz_soisno(j) / (0.5*(z_soisno(j)-zi_soisno(j-1)+capr*(z_soisno(j+1)-zi_soisno(j-1))))

   DO j = lb + 1, nl_soil
      fact(j) = deltim/cv(j)
   ENDDO

   DO j = lb, nl_soil - 1
      fn(j) = tk(j)*(t_soisno(j+1)-t_soisno(j))/(z_soisno(j+1)-z_soisno(j))
   ENDDO
   fn(nl_soil) = 0.

! set up vector r and vectors a, b, c that define tridiagonal matrix
   j     = lb
   dzp   = z_soisno(j+1)-z_soisno(j)
   at(j) = 0.
   bt(j) = 1+(1.-cnfac)*fact(j)*tk(j)/dzp-fact(j)*dhsdT
   ct(j) =  -(1.-cnfac)*fact(j)*tk(j)/dzp
   rt(j) = t_soisno(j) + fact(j)*( hs - dhsdT*t_soisno(j) + cnfac*fn(j) )

! January 12, 2023
   IF (lb <= 0) THEN
       DO j = lb + 1, 1
          dzm   = (z_soisno(j)-z_soisno(j-1))
          dzp   = (z_soisno(j+1)-z_soisno(j))
          at(j) =   - (1.-cnfac)*fact(j)* tk(j-1)/dzm
          bt(j) = 1.+ (1.-cnfac)*fact(j)*(tk(j)/dzp + tk(j-1)/dzm)
          ct(j) =   - (1.-cnfac)*fact(j)* tk(j)/dzp
          rt(j) = t_soisno(j) + fact(j)*sabg_lyr(j) + cnfac*fact(j)*( fn(j) - fn(j-1) )
       ENDDO
   ENDIF

   DO j = 2, nl_soil - 1
! January 12, 2023
      dzm   = (z_soisno(j)-z_soisno(j-1))
      dzp   = (z_soisno(j+1)-z_soisno(j))
      at(j) =   - (1.-cnfac)*fact(j)* tk(j-1)/dzm
      bt(j) = 1.+ (1.-cnfac)*fact(j)*(tk(j)/dzp + tk(j-1)/dzm)
      ct(j) =   - (1.-cnfac)*fact(j)* tk(j)/dzp
      rt(j) = t_soisno(j) + cnfac*fact(j)*( fn(j) - fn(j-1) )
   ENDDO

   j     =  nl_soil
   dzm   = (z_soisno(j)-z_soisno(j-1))
   at(j) =   - (1.-cnfac)*fact(j)*tk(j-1)/dzm
   bt(j) = 1.+ (1.-cnfac)*fact(j)*tk(j-1)/dzm
   ct(j) = 0.
   rt(j) = t_soisno(j) - cnfac*fact(j)*fn(j-1)

! solve for t_soisno
   i = size(at)
   CALL tridia (i ,at ,bt ,ct ,rt ,t_soisno)
!=======================================================================
! melting or freezing
!=======================================================================

   DO j = lb, nl_soil - 1
      fn1(j) = tk(j)*(t_soisno(j+1)-t_soisno(j))/(z_soisno(j+1)-z_soisno(j))
   ENDDO
   fn1(nl_soil) = 0.

   j = lb
   brr(j) = cnfac*fn(j) + (1.-cnfac)*fn1(j)

   DO j = lb + 1, nl_soil
      brr(j) = cnfac*(fn(j)-fn(j-1)) + (1.-cnfac)*(fn1(j)-fn1(j-1))
   ENDDO


   IF (DEF_USE_SNICAR) THEN

      wice_soisno_bef(lb:0) = wice_soisno(lb:0)

      CALL meltf_snicar (patchtype,lb,nl_soil,deltim, &
               fact(lb:),brr(lb:),hs,dhsdT,sabg_lyr, &
               t_soisno_bef(lb:),t_soisno(lb:),wliq_soisno(lb:),wice_soisno(lb:),imelt(lb:), &
               scv,snowdp,sm,xmf,porsl,psi0,&
#ifdef Campbell_SOIL_MODEL
               bsw,&
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
               theta_r,alpha_vgm,n_vgm,L_vgm,&
               sc_vgm,fc_vgm,&
#endif
               dz_soisno(1:nl_soil))

      ! layer freezing mass flux (positive):
      DO j = lb, 0
         IF (imelt(j)==2 .and. j<1) THEN
             snofrz(j) = max(0._r8,(wice_soisno(j)-wice_soisno_bef(j)))/deltim
         ENDIF
      ENDDO

   ELSE
      CALL meltf (patchtype,lb,nl_soil,deltim, &
               fact(lb:),brr(lb:),hs,dhsdT, &
               t_soisno_bef(lb:),t_soisno(lb:),wliq_soisno(lb:),wice_soisno(lb:),imelt(lb:), &
               scv,snowdp,sm,xmf,porsl,psi0,&
#ifdef Campbell_SOIL_MODEL
               bsw,&
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
               theta_r,alpha_vgm,n_vgm,L_vgm,&
               sc_vgm,fc_vgm,&
#endif
               dz_soisno(1:nl_soil))
   ENDIF

!-----------------------------------------------------------------------

   END SUBROUTINE GroundTemperature

END MODULE MOD_GroundTemperature
