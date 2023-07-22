#include <define.h>

MODULE MOD_SoilSurfaceResistance
   ! -----------------------------------------------------------------------
   ! !DESCRIPTION:
   ! Calculate the soil surface resistance by using three parameterization schemes
   !
   ! ORIGINAL:
   ! Zhuo Liu, June, 2023
   !
   ! -----------------------------------------------------------------------
   ! !USE

   USE MOD_Precision
   IMPLICIT NONE
   SAVE

   PUBLIC :: SoilSurfaceResistance

   !TODO@Zhuo: need complement below
   ! soil-gas diffusivity schemes:
   ! 1: BBC (Buckingham-Burdine-Campbell Model), Moldrup et al., 1999.
   ! 2: P_WLR (Penman Water Linear Reduction Model), Moldrup et al., 2000
   ! 3: MI_WLR (Millington Water Linear Reduction Model), Moldrup et al., 2000
   ! 4: MA_WLR (Marshal Water Linear Reduction Model), Moldrup et al., 2000
   ! 5: M_Q, Millington and Quirk, 1961
   ! 6: 3POE (Three-Porosity-Encased), Moldrup et al., 2005
   integer, parameter :: soil_gas_diffusivity_scheme = 1


CONTAINS
!-----------------------------------------------------------------------

 SUBROUTINE SoilSurfaceResistance (nl_soil,forc_rhoair,hksati,porsl,bsw,psi0,&
                   dz_soisno,t_soisno,wliq_soisno,wice_soisno,fsno,qg,rss)

  !=======================================================================
  ! !DESCRIPTION:
  !
  ! REFERENCES:
  !
  !
  !=======================================================================

   USE MOD_Precision
   USE MOD_Const_Physical, only: denice, denh2o  ! physical constant
   USE MOD_Namelist, only: DEF_RSS_SCHEME
   IMPLICIT NONE


!-----------------------Argument-----------------------------------------

   integer, intent(in) :: &
        nl_soil                       ! upper bound of array

   real(r8), intent(in) :: &
        porsl           (1:nl_soil), &! soil porosity [-]
        bsw             (1:nl_soil), &! Clapp-Hornberger "B"
        psi0            (1:nl_soil), &! saturated soil suction [mm] (NEGATIVE)
        dz_soisno       (1:nl_soil), &! layer thickness [m]
        t_soisno        (1:nl_soil), &! soil/snow skin temperature [K]
        wliq_soisno     (1:nl_soil), &! liquid water [kg/m2]
        wice_soisno     (1:nl_soil), &! ice lens [kg/m2]
        hksati          (1:nl_soil), &! hydraulic conductivity at saturation [mm h2o/s] 
        qg,                          &! ground specific humidity [kg/kg]      
        fsno,                        &! snow cover 
        forc_rhoair                   ! density air [kg/m**3]
   real(r8), intent(out) :: &
        rss                           ! soil surface resistance [s/m]

!-----------------------Local Variables------------------------------

   REAL(r8) :: &
        wx,               & ! patitial volume of ice and water of surface layer   
        vol_liq,          & ! water content by volume [m3/m3]
        smp_node,         & ! matrix potential [m]
        eff_porosity,     & ! effective porosity = porosity - vol_ice
        aird,             & ! “air-dry” soil moisture value 
        d0,               & ! water vapor diffusivity in open air [m2/s]
        eps,              & ! air filled pore space
        dg,               & ! gaseous diffusivity [m2/s]
        dsl,              & ! soil dry surface layer thickness [m]
        dw,               & ! aqueous diffusivity [m2/s]
        hk,               & ! hydraulic conductivity [m h2o/s]
        wfc,              & ! field capacity of the first layer soil
        rg_1,             & ! inverse of vapor diffusion resistance [m/s]
        rw_1,             & ! inverse of volatilization resistance [m/s]
        rss_1,            & ! inverse of soil surface resistance [m/s]
        tao,              & ! tortuosity of the vapor flow paths through the soil matrix
        eps100,           & ! air-filled porosity (cm3 soil-air cm−3 soil) at −1000 mm of water matric potential
        fac,              & ! temporal variable for calculating wx/porsl 
        fac_fc,           & ! temporal variable for calculating wx/wfc
        B                   ! bunsen solubility coefficient
       
        !TODO@Zhuo Liu:  need descriptions for vars.
        ! if it is a temporal varialbe, add descriptions like 'temporal variable for calculating ***'
!-----------------------End Variables list---------------------------


   !calculate the top soil volumetric water content (m3/m3), soil matrix potential and soil hydraulic conductivity
   vol_liq      = max(wliq_soisno(1),1.0e-6_r8)/(denh2o*dz_soisno(1))
   smp_node     = (psi0(1)/1000.)*(vol_liq/porsl(1))**(-bsw(1))
   hk           = (hksati(1)/1000.)*(vol_liq/porsl(1))**(2.*bsw(1)+3.)

   !eff_porosity not calculated til SoilHydrolog
   eff_porosity = max(0.01_r8,porsl(1)-min(porsl(1), wice_soisno(1)/(dz_soisno(1)*denice)))
   
   !calculate diffusivity (dg, dw) and air free pore space
   aird = porsl(1)*(psi0(1)/-1.e7_r8)**(1./bsw(1))
   d0   = 2.12e-5*(t_soisno(1)/273.15)**1.75
   eps  = porsl(1) - aird
   
   select case (soil_gas_diffusivity_scheme)

   ! 1: BBC
   case (1)
      tao = eps*eps*(eps/porsl(1))**(3._r8/max(3._r8,bsw(1)))

   ! 2: P_WLR
   case (2)
      tao = 0.66*eps*(eps/porsl(1))

   ! 3: MI_WLR
   case (3)
      tao = eps**(4._r8/3._r8)*(eps/porsl(1))

   ! 4: MA_WLR
   case (4)
      tao = eps**(3./2.)*(eps/porsl(1))

   ! 5: M_Q
   case (5)
      tao = eps**(4._r8/3._r8)*(eps/porsl(1))**(2.0_r8)

   ! 6: 3POE
   case (6)
      eps100 = porsl(1) - porsl(1)*(psi0(1)/-1000.)**(1./bsw(1))
      tao    = porsl(1)*porsl(1)*(eps/porsl(1))**(2.+log(eps100**0.25_r8)/log(eps100/porsl(1)))

   endselect

   dg = d0*tao
   dw = -hk*bsw(1)*smp_node/vol_liq


   select case (DEF_RSS_SCHEME)

   ! calculate rss by SL14
   case (1)
      dsl = dz_soisno(1)*max(1.e-6_r8,(0.8*eff_porosity - vol_liq)) &
                        /max(1.e-6_r8,(0.8*porsl(1)- aird)) ! SL14

      dsl = max(dsl,0._r8)
      dsl = min(dsl,0.2_r8)
      
      rss = dsl/dg

   ! calculate rss by SZ09   
   case (2)
      dsl = dz_soisno(1)*(exp((1._r8 - vol_liq/eff_porosity)**5) - 1._r8)/ (exp(1._r8) - 1._r8)
      dsl = min(dsl,0.2_r8)
      dsl = max(dsl,0._r8)

      rss = dsl/dg

   ! calculate rss by TR13
   case (3)
      B           = denh2o/(qg*forc_rhoair)                       ! Eq.(12)
      rg_1        = 2.0_r8*dg*eps/dz_soisno(1)
      rw_1        = 2.0_r8*dw*B*vol_liq/dz_soisno(1)
      rss_1       = rg_1 +rw_1
      rss         = 1.0/rss_1
     
   ! LP92 beta scheme
   case (4)
      wx  = (max(wliq_soisno(1),1.e-6)/denh2o+wice_soisno(1)/denice)/dz_soisno(1)
      fac = min(1._r8, wx/porsl(1))
      fac = max(fac , 0.001_r8)
      wfc = porsl(1)*(-3399._r8/psi0(1))**(-1./bsw(1))
      !! Lee and Pielke 1992 beta
      IF (wx < wfc ) THEN  !when water content of ths top layer is less than that at F.C.
         fac_fc = min(1._r8, wx/wfc)
         fac_fc = max(fac_fc,0.001_r8)
         ! modify soil beta by snow cover. soilbeta for snow surface is one
         rss  = (1._r8-fsno) &
                              *0.25_r8*(1._r8 - cos(fac_fc*3.1415926))**2._r8
      ELSE   ! when water content of ths top layer is more than that at F.C.
         rss  = 1._r8
      ENDIF
      ! write(*,*) wfc,porsl(1),psi0(1),bsw(1)
      ! raw = 50m/s NOTE: raw = 50m/s
      ! rss = 50._r8 * (1._r8/beta - 1._r8)

   ! Sellers 1992
   case (5)
      wx  = (max(wliq_soisno(1),1.e-6)/denh2o+wice_soisno(1)/denice)/dz_soisno(1)
      fac = min(1._r8, wx/porsl(1))
      fac = max(fac , 0.001_r8)
      !! Lee and Pielke 1992 beta
      rss = (1-fsno)*exp(8.206-4.255*fac)

   endselect

   rss = min(1.e6_r8,rss)
   write(*,*) porsl(1)

 END Subroutine SoilSurfaceResistance

END MODULE MOD_SoilSurfaceResistance
