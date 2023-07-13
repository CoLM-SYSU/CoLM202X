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
   !USE MOD_Vars_Global
   !USE MOD_Const_Physical, only : denice, denh2o  ! physical constant
   !USE MOD_Vars_TimeInvariants
   IMPLICIT NONE
   SAVE

   PUBLIC :: SoilSurfaceResistance


CONTAINS
!-----------------------------------------------------------------------

 SUBROUTINE SoilSurfaceResistance (nl_soil,forc_rhoair,hksati,porsl,bsw,psi0,&
                   dz_soisno,t_soisno,wliq_soisno,wice_soisno,fsno,wfc,qg,rss)
                   
  !=======================================================================
  ! !DESCRIPTION:
  !
  ! REFERENCES:
  ! 
  !
  !=======================================================================

   USE MOD_Precision
   USE MOD_Const_Physical, only : denice, denh2o  ! physical constant
   IMPLICIT NONE

    
!-----------------------Argument-----------------------------------------

   integer, intent(in) :: &
        nl_soil                       ! upper bound of array

   real(r8), intent(in) :: &
        porsl           (1:nl_soil), &! soil porosity [-]
        bsw             (1:nl_soil), &! Clapp-Hornberger "B"
        psi0            (1:nl_soil), &! saturated soil suction (mm) (NEGATIVE)
        dz_soisno       (1:nl_soil), &! layer thickness (m)
        t_soisno        (1:nl_soil), &! soil/snow skin temperature (K)
        wliq_soisno     (1:nl_soil), &! liquid water (kg/m2)
        wice_soisno     (1:nl_soil), &! ice lens [kg/m2]
        hksati          (1:nl_soil), &! hydraulic conductivity at saturation (mm h2o/s) 
        wfc             (1:nl_soil), &! fc
        qg,                          &! ground specific humidity [kg/kg]      
        fsno,                        &! snow cover  
        forc_rhoair                   ! density air [kg/m**3]
   real(r8), intent(out) :: &
        rss                           ! soil surface resistance(m/s)


                   
!-----------------------Local Variables------------------------------
                   


   REAL(r8) :: &
        wx,               & ! soil wetness   
        vol_liq,          & ! vol_liq
        smp_node,         & ! matrix potential
        eff_porosity,     & ! effective porosity = porosity - vol_ice
        aird,             & ! air free pore space
        d0,               & ! molecular diffusivity of water vapor in air (m2/s)
        eps,              & ! air filled pore space
        dg,               & ! d0*(tortuosity of the vapor flow paths through the soil matrix)
        dsl,              & ! soil dry surface layer thickness
        dw,               & ! aqueous diffusivity (m2/s)
        hk,               & ! hydraulic conductivity
        d1,               & !
        beta,             & !
        tao,              & !
        eps100,           & !
        fac,              & !
        fac_fc,           & ! 
        B                   ! liquid water density / water vapor density
       
!-----------------------End Variables list---------------------------
 

   !calculate the top soil volumetric water content (m3/m3), soil matrix potential and soil hydraulic conductivity
  
   vol_liq      = max(wliq_soisno(1),1.0e-6_r8)/(denh2o*dz_soisno(1))
   smp_node     = (psi0(1)/1000.)*(vol_liq/porsl(1))**(-bsw(1))
   hk           = (hksati(1)/1000.)*(vol_liq/porsl(1))**(2.*bsw(1)+3.)

   !eff_porosity not calculated til SoilHydrolog

   eff_porosity = max(0.01_r8,porsl(1)-min(porsl(1), wice_soisno(1)/(dz_soisno(1)*denice)))


   !calculate diffusivity (dg, dw) and air free pore space
   aird        = porsl(1)*(psi0(1)/-1.e7_r8)**(1./bsw(1))
   d0          = 2.12e-5*(t_soisno(1)/273.15)**1.75 ![Bitelli et al., JH, 08]
   eps         = porsl(1) - aird

#ifdef BBC
   tao         = eps*eps*(eps/porsl(1))**(3._r8/max(3._r8,bsw(1)))
#endif

#ifdef P_WLR
   tao         = 0.66*eps*(eps/porsl(1))
#endif

#ifdef MI_WLR
   tao         = eps**(4._r8/3._r8)*(eps/porsl(1))
#endif

#ifdef MA_WLR
   tao         = eps**(3./2.)*(eps/porsl(1))
#endif

#ifdef M_Q
   tao         = eps**(4._r8/3._r8)*(eps/porsl(1))**(2.0_r8)
#endif

#ifdef POE 
   eps100      = porsl(1) - porsl(1)*(psi0(1)/-1000.)**(1./bsw(1))
   tao         = porsl(1)*porsl(1)*(eps/porsl(1))**(2.+log(eps100**0.25_r8)/log(eps100/porsl(1)))
#endif

   dg          = d0*tao
   dw          = -hk*bsw(1)*smp_node/vol_liq


  ! calculate dsl by SL14
#ifdef RSS_SL14
   dsl         = dz_soisno(1)*max(1.e-6_r8,(0.8*eff_porosity - vol_liq)) &
   /max(1.e-6_r8,(0.8*porsl(1)- aird)) ! SL14

   dsl         = max(dsl,0._r8)
   dsl         = min(dsl,0.2_r8)
   rss         = dsl/dg + 20._r8
#endif
  
#ifdef RSS_SZ09
   ! calculate dsl by SZ09
   dsl         = dz_soisno(1)*(exp((1._r8 - vol_liq/eff_porosity)**5) - 1._r8)/ (exp(1) - 1._r8)   ! SZ09
   dsl         = min(dsl,0.2_r8)
   dsl         = max(dsl,0._r8)

   rss         = dsl/dg

#endif

#ifdef RSS_TR13
   ! calculate dsl by TR13 
   B           = denh2o/(qg*forc_rhoair)                       ! Eq.(12)
   d1          = (B*vol_liq*dw +eps*dg)/(B*vol_liq+eps)        ! Eq.(9)
   dsl         = dz_soisno(1)*(1._r8/(2._r8*(B*vol_liq+eps)))  ! Eq.(10) 
   dsl         = min(dsl,0.2_r8)
   dsl         = max(dsl,0._r8)
   rss         = dsl/d1

#endif

#ifdef beta
   
  
   wx      = (max(wliq_soisno(1),1.e-6)/denh2o+wice_soisno(1)/denice)/dz_soisno(1)
   fac     = min(1._r8, wx/porsl(1))
   fac     = max(fac , 0.001_r8)
   !! Lee and Pielke 1992 beta
   IF (wx < wfc(1) ) THEN  !when water content of ths top layer is less than that at F.C.
   fac_fc  = min(1._r8, wx/wfc(1))
   fac_fc  = max(fac_fc,0.001_r8)
          ! modify soil beta by snow cover. soilbeta for snow surface is one
   beta    = (1._r8-fsno) &
                        *0.25_r8*(1._r8 - cos(fac_fc*3.1415926))**2._r8 
   ELSE   ! when water content of ths top layer is more than that at F.C.
   beta    = 1._r8
   ENDIF
   ! raw = 50m/s NOTE: raw = 50m/s  
   rss     = 50._r8 * (1._r8/soilbeta - 1._r8)

#endif

 
#ifdef Sellers
   wx      = (max(wliq_soisno(1),1.e-6)/denh2o+wice_soisno(1)/denice)/dz_soisno(1)
   fac     = min(1._r8, wx/porsl(1))
   fac     = max(fac , 0.001_r8)
   !! Lee and Pielke 1992 beta
   rss     = (1-fsno)*exp(8.206-4.255*fac)
#endif

   rss     = min(1.e6_r8,rss)

 END Subroutine SoilSurfaceResistance

END MODULE MOD_SoilSurfaceResistance
