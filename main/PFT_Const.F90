#include <define.h>

MODULE PFT_Const
! -------------------------------
! Created by Hua Yuan, 08/2019
! -------------------------------

   USE precision
   USE GlobalVars

   IMPLICIT NONE
   SAVE

! Plant Functional Type classification 
!---------------------------
! 0  not vegetated
! 1  needleleaf evergreen temperate tree 
! 2  needleleaf evergreen boreal tree
! 3  needleleaf deciduous boreal tree
! 4  broadleaf evergreen tropical tree
! 5  broadleaf evergreen temperate tree
! 6  broadleaf deciduous tropical tree
! 7  broadleaf deciduous temperate tree
! 8  broadleaf deciduous boreal tree
! 9  broadleaf evergreen shrub
!10  broadleaf deciduous temperate shrub
!11  broadleaf deciduous boreal shrub
!12  c3 arctic grass
!13  c3 non-arctic grass
!14  c4 grass
!15  c3 crop
   
   ! canopy layer number
   INTEGER , parameter :: canlay(0:15) &
      = (/0, 2, 2, 2, 2, 2, 2, 2, &
          2, 1, 1, 1, 1, 1, 1, 1 /)

   ! canopy top height
   REAL(r8), parameter :: htop0_p(0:15) &
      =(/ 0.5,  17.0,  17.0,  14.0,  35.0,  35.0,  18.0,  20.0,&
         20.0,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5/)

   ! canopy bottom height
   REAL(r8), parameter :: hbot0_p(0:15) &
! 01/06/2020, yuan: adjust htop: grass/shrub -> 0, tree->1
      !=(/0.01,   8.5,   8.5,   7.0,   1.0,   1.0,  10.0,  11.5,&
      !   11.5,   0.1,   0.1,   0.1,  0.01,  0.01,  0.01,  0.01/)
      =(/0.00,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,&
          1.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0/)

   ! defulat vegetation fractional cover
   REAL(r8), parameter :: fveg0_p(0:15) &
      = 1.0 !(/.../)

   ! default stem area index
   REAL(r8), parameter :: sai0_p(0:15) &
      =(/0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,&
         2.0, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.2/)

   ! ratio to calculate roughness length z0m
   REAL(r8), parameter :: z0mr_p(0:15) = 0.1 
   
   ! ratio to calculate displacement height d 
   REAL(r8), parameter :: displar_p(0:15) = 0.667

   ! inverse&sqrt leaf specific dimension size 4 cm
   REAL(r8), parameter :: sqrtdi_p(0:15) = 5.0
      
   ! leaf angle distribution parameter
   REAL(r8), parameter :: chil_p(0:15) &
      = (/-0.300,  0.010,  0.010,  0.010,  0.100,  0.100,  0.010,  0.250,&
           0.250,  0.010,  0.250,  0.250, -0.300, -0.300, -0.300, -0.300/)

   ! reflectance of green leaf in virsible band 
   REAL(r8), parameter :: rhol_vis_p(0:15) &
      = (/0.110,  0.070,  0.070,  0.070,  0.100,  0.100,  0.100,  0.100,&
          0.100,  0.070,  0.100,  0.100,  0.110,  0.110,  0.110,  0.110/)

   ! reflectance of dead leaf in virsible band 
   REAL(r8), parameter :: rhos_vis_p(0:15) &
      = (/0.310,  0.160,  0.160,  0.160,  0.160,  0.160,  0.160,  0.160,&
          0.160,  0.160,  0.160,  0.160,  0.310,  0.310,  0.310,  0.310/)

   ! reflectance of green leaf in near infrared band 
   REAL(r8), parameter :: rhol_nir_p(0:15) &
      = (/0.350,  0.350,  0.350,  0.350,  0.450,  0.450,  0.450,  0.450,&
          0.450,  0.350,  0.450,  0.450,  0.350,  0.350,  0.350,  0.350/)

   ! reflectance of dead leaf in near infrared band 
   REAL(r8), parameter :: rhos_nir_p(0:15) &
      = (/0.530,  0.390,  0.390,  0.390,  0.390,  0.390,  0.390,  0.390,&
          0.390,  0.390,  0.390,  0.390,  0.530,  0.530,  0.530,  0.530/)

   ! transmittance of green leaf in visible band 
   REAL(r8), parameter :: taul_vis_p(0:15) &
      = (/0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,&
          0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050/)

   ! transmittance of dead leaf in visible band 
   REAL(r8), parameter :: taus_vis_p(0:15) &
      = (/0.120,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,&
          0.001,  0.001,  0.001,  0.001,  0.120,  0.120,  0.120,  0.120/)

   ! transmittance of green leaf in near infrared band 
   REAL(r8), parameter :: taul_nir_p(0:15) &
      = (/0.340,  0.100,  0.100,  0.100,  0.250,  0.250,  0.250,  0.250,&
          0.250,  0.100,  0.250,  0.250,  0.340,  0.340,  0.340,  0.340/)

   ! transmittance of dead leaf in near infrared band 
   REAL(r8), parameter :: taus_nir_p(0:15) &
      = (/0.250,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,&
          0.001,  0.001,  0.001,  0.001,  0.250,  0.250,  0.250,  0.250/)

   ! maximum carboxylation rate at 25 C at canopy top
   ! /06/03/2014/ based on Bonan et al., 2010 (Table 2)
   REAL(r8), parameter :: vmax25_p(0:15) &
      = (/ 52.0, 61.0, 54.0, 57.0, 72.0, 72.0, 52.0, 52.0,&
           52.0, 72.0, 52.0, 52.0, 52.0, 52.0, 52.0, 57.0/) * 1.e-6

   ! quantum efficiency
   REAL(r8), parameter :: effcon_p(0:15) &
      = (/0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08,&
          0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.05, 0.08/)

   ! conductance-photosynthesis slope parameter
   REAL(r8), parameter :: gradm_p(0:15) &
      = (/9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0,&
          9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 4.0, 9.0/)

   ! conductance-photosynthesis intercept
   REAL(r8), parameter :: binter_p(0:15) &
      = (/0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,&
          0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.04, 0.01/)

   ! respiration fraction
   REAL(r8), parameter :: respcp_p(0:15) &
      = (/0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015,&
          0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.025, 0.015/)

   ! slope of high temperature inhibition FUNCTION (s1)
   REAL(r8), parameter :: shti_p(0:15) = 0.3

   ! slope of low temperature inhibition FUNCTION (s3) 
   REAL(r8), parameter :: slti_p(0:15) = 0.2

   ! temperature coefficient in gs-a model (s5)
   REAL(r8), parameter :: trda_p(0:15) = 1.3
   
   ! temperature coefficient in gs-a model (s6)
   REAL(r8), parameter :: trdm_p(0:15) = 328.0

   ! temperature coefficient in gs-a model (273.16+25)
   REAL(r8), parameter :: trop_p(0:15) = 298.0

   ! 1/2 point of high temperature inhibition FUNCTION (s2)
   REAL(r8), parameter :: hhti_p(0:15) &
      =(/308.0, 303.0, 303.0, 303.0, 313.0, 313.0, 311.0, 311.0,&
         311.0, 313.0, 313.0, 303.0, 303.0, 308.0, 313.0, 308.0/)

   ! 1/2 point of low temperature inhibition FUNCTION (s4)
   REAL(r8), parameter :: hlti_p(0:15) &
      =(/281.0, 278.0, 278.0, 278.0, 288.0, 288.0, 283.0, 283.0,&
         283.0, 283.0, 283.0, 278.0, 278.0, 281.0, 288.0, 281.0/)

   ! coefficient of leaf nitrogen allocation
   REAL(r8), parameter :: extkn_p(0:15) = 0.5

   REAL(r8) :: &
      rho_p(2,2,0:15), &!leaf reflectance
      tau_p(2,2,0:15)   !leaf transmittance

   ! depth at 50% roots
   REAL(r8), parameter, dimension(0:15) :: d50_p &
      =(/27.0,  21.0,  12.0,  12.0,  15.0,  23.0,  16.0,  23.0,&
         12.0,  23.5,  23.5,  23.5,   9.0,   7.0,  16.0,  22.0/)

   ! coefficient of root profile
   REAL(r8), parameter, dimension(0:15) :: beta_p &
      =(/-2.051, -1.835, -1.880, -1.880, -1.632, -1.757, -1.681, -1.757,&
         -1.880, -1.623, -1.623, -1.623, -2.621, -1.176, -1.452, -1.796/)

   !设定PFT的根分布参数 
   REAL(r8), PRIVATE, parameter :: roota(0:15) &
      =(/  0.0,   7.0,   7.0,   7.0,   7.0,   7.0,   6.0,   6.0,&
           6.0,   7.0,   7.0,   7.0,  11.0,  11.0,  11.0,   6.0/)
   
   REAL(r8), PRIVATE, parameter :: rootb(0:15) &
      =(/  0.0,   2.0,   2.0,   2.0,   1.0,   1.0,   2.0,   2.0,&
           2.0,   1.5,   1.5,   1.5,   2.0,   2.0,   2.0,   3.0/)

   ! scheme 1: Zeng 2001, 2: Schenk and Jackson, 2002
   INTEGER, PRIVATE :: ROOTFR_SCHEME = 1 
   
   !fraction of roots in each soil layer 
   REAL(r8), dimension(nl_soil,16) :: &
      rootfr_p(1:nl_soil, 0:15) 

   INTEGER, PRIVATE :: i, nsl
   
   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: Init_PFT_Const

CONTAINS

   SUBROUTINE Init_PFT_Const

      IMPLICIT NONE

      rho_p(1,1,:) = rhol_vis_p(:)
      rho_p(2,1,:) = rhol_nir_p(:)
      rho_p(1,2,:) = rhos_vis_p(:)
      rho_p(2,2,:) = rhos_nir_p(:)
      tau_p(1,1,:) = taul_vis_p(:)
      tau_p(2,1,:) = taul_nir_p(:)
      tau_p(1,2,:) = taus_vis_p(:)
      tau_p(2,2,:) = taus_nir_p(:)

IF (ROOTFR_SCHEME == 1) THEN
      DO i = 0, N_PFT-1
         rootfr_p(1,i)=1./(1.+(z_soih(1)*100./d50_p(i))**beta_p(i)) 
         rootfr_p(nl_soil,i)=1.-1./(1.+(z_soih(nl_soil-1)*100./d50_p(i))**beta_p(i)) 

         DO nsl=2,nl_soil-1
            rootfr_p(nsl,i)=1./(1.+(z_soih(nsl)*100./d50_p(i))**beta_p(i)) &
               -1./(1.+(z_soih(nsl-1)*100./d50_p(i))**beta_p(i))
         ENDDO
      ENDDO 
ELSE
      ! PFT rootfr_p (Zeng, 2001)
      DO i = 0, N_PFT-1
         rootfr_p(1,i) = 1. - 0.5*( &
              exp(-roota(i) * z_soih(1)) &
            + exp(-rootb(i) * z_soih(1)) )
            
         rootfr_p(nl_soil,i) = 0.5*( &
              exp(-roota(i) * z_soih(nl_soil)) &
            + exp(-rootb(i) * z_soih(nl_soil)) )
         
         DO nsl = 2, nl_soil-1
            rootfr_p(nsl,i) = 0.5*( &
                 exp(-roota(i) * z_soih(nsl-1)) &
               + exp(-rootb(i) * z_soih(nsl-1)) &
               - exp(-roota(i) * z_soih(nsl)) &
               - exp(-rootb(i) * z_soih(nsl)) )
         ENDDO
      ENDDO
ENDIF
   
   END SUBROUTINE Init_PFT_Const

END MODULE PFT_Const
