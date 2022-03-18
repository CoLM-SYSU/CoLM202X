#include <define.h>

MODULE PFT_Const
! -------------------------------
! Created by Hua Yuan, 08/2019
! -------------------------------

   USE precision
   USE GlobalVars
   USE timemanager, only: get_calday

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
!16  c3_irrigated
!17  temperate_corn
!18  irrigated_temperate_corn
!19  spring_wheat
!20  irrigated_spring_wheat
!21  winter_wheat
!22  irrigated_winter_wheat
!23  temperate_soybean
!24  irrigated_temperate_soybean
!25  barley
!26  irrigated_barley
!27  winter_barley
!28  irrigated_winter_barley
!29  rye
!30  irrigated_rye
!31  winter_rye
!32  irrigated_winter_rye
!33  cassava
!34  irrigated_cassava
!35  citrus
!36  irrigated_citrus
!37  cocoa
!38  irrigated_cocoa
!39  coffee
!40  irrigated_coffee
!41  cotton
!42  irrigated_cotton
!43  datepalm
!44  irrigated_datepalm
!45  foddergrass
!46  irrigated_foddergrass
!47  grapes
!48  irrigated_grapes
!49  groundnuts
!50  irrigated_groundnuts
!51  millet
!52  irrigated_millet
!53  oilpalm
!54  irrigated_oilpalm
!55  potatoes
!56  irrigated_potatoes
!57  pulses
!58  irrigated_pulses
!59  rapeseed
!60  irrigated_rapeseed
!61  rice
!62  irrigated_rice
!63  sorghum
!64  irrigated_sorghum
!65  sugarbeet
!66  irrigated_sugarbeet
!67  sugarcane
!68  irrigated_sugarcane
!69  sunflower
!70  irrigated_sunflower
!71  miscanthus
!72  irrigated_miscanthus
!73  switchgrass
!74  irrigated_switchgrass
!75  tropical_corn
!76  irrigated_tropical_corn
!77  tropical_soybean
!78  irrigated_tropical_soybean
   
   ! canopy layer number
   INTEGER , parameter :: canlay(0:78) &
      = (/0, 2, 2, 2, 2, 2, 2, 2, &
          2, 1, 1, 1, 1, 1, 1, 1, &
          1, 1, 1, 1, 1, 1, 1, 1, &
          1, 1, 1, 1, 1, 1, 1, 1, &
          1, 1, 1, 1, 1, 1, 1, 1, &
          1, 1, 1, 1, 1, 1, 1, 1, &
          1, 1, 1, 1, 1, 1, 1, 1, &
          1, 1, 1, 1, 1, 1, 1, 1, &
          1, 1, 1, 1, 1, 1, 1, 1, &
          1, 1, 1, 1, 1, 1, 1/)

   ! canopy top height
   REAL(r8), parameter :: htop0_p(0:78) &
      =(/ 0.5,  17.0,  17.0,  14.0,  35.0,  35.0,  18.0,  20.0,&
         20.0,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,&
          0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,&
          0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,&
          0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,&
          0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,&
          0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,&
          0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,&
          0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,&
          0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5/)

   ! canopy bottom height
   REAL(r8), parameter :: hbot0_p(0:78) &
! 01/06/2020, yuan: adjust htop: grass/shrub -> 0, tree->1
      !=(/0.01,   8.5,   8.5,   7.0,   1.0,   1.0,  10.0,  11.5,&
      !   11.5,   0.1,   0.1,   0.1,  0.01,  0.01,  0.01,  0.01/)
      =(/0.00,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,&
          1.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,&
          0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,&
          0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,&
          0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,&
          0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,&
          0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,&
          0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,&
          0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,&
          0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0/)

   ! defulat vegetation fractional cover
   REAL(r8), parameter :: fveg0_p(0:78) &
      = 1.0 !(/.../)

   ! default stem area index
   REAL(r8), parameter :: sai0_p(0:78) &
      =(/0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,&
         2.0, 0.5, 0.5, 0.5, 0.2, 0.2, 0.2, 0.2,&
         0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,&
         0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,&
         0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,&
         0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,&
         0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,&
         0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,&
         0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,&
         0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2/)

   ! ratio to calculate roughness length z0m
   REAL(r8), parameter :: z0mr_p(0:78) = 0.1 
   
   ! ratio to calculate displacement height d 
   REAL(r8), parameter :: displar_p(0:78) = 0.667

   ! inverse&sqrt leaf specific dimension size 4 cm
   REAL(r8), parameter :: sqrtdi_p(0:78) = 5.0
      
   ! leaf angle distribution parameter
   REAL(r8), parameter :: chil_p(0:N_PFT+N_CFT-1) &
      = (/-0.300,  0.010,  0.010,  0.010,  0.100,  0.100,  0.010,  0.250,&
           0.250,  0.010,  0.250,  0.250, -0.300, -0.300, -0.300, -0.300,&
          -0.300, -0.300, -0.300, -0.300, -0.300, -0.300, -0.300, -0.300,&
          -0.300, -0.300, -0.300, -0.300, -0.300, -0.300, -0.300, -0.300,&
          -0.300, -0.300, -0.300, -0.300, -0.300, -0.300, -0.300, -0.300,&
          -0.300, -0.300, -0.300, -0.300, -0.300, -0.300, -0.300, -0.300,&
          -0.300, -0.300, -0.300, -0.300, -0.300, -0.300, -0.300, -0.300,&
          -0.300, -0.300, -0.300, -0.300, -0.300, -0.300, -0.300, -0.300,&
          -0.300, -0.300, -0.300, -0.300, -0.300, -0.300, -0.300, -0.300,&
          -0.300, -0.300, -0.300, -0.300, -0.300, -0.300, -0.300/)

   ! reflectance of green leaf in virsible band 
   REAL(r8), parameter :: rhol_vis_p(0:78) &
      = (/0.110,  0.070,  0.070,  0.070,  0.100,  0.100,  0.100,  0.100,&
          0.100,  0.070,  0.100,  0.100,  0.110,  0.110,  0.110,  0.110,&
          0.110,  0.110,  0.110,  0.110,  0.110,  0.110,  0.110,  0.110,&
          0.110,  0.110,  0.110,  0.110,  0.110,  0.110,  0.110,  0.110,&
          0.110,  0.110,  0.110,  0.110,  0.110,  0.110,  0.110,  0.110,&
          0.110,  0.110,  0.110,  0.110,  0.110,  0.110,  0.110,  0.110,&
          0.110,  0.110,  0.110,  0.110,  0.110,  0.110,  0.110,  0.110,&
          0.110,  0.110,  0.110,  0.110,  0.110,  0.110,  0.110,  0.110,&
          0.110,  0.110,  0.110,  0.110,  0.110,  0.110,  0.110,  0.110,&
          0.110,  0.110,  0.110,  0.110,  0.110,  0.110,  0.110/)

   ! reflectance of dead leaf in virsible band 
   REAL(r8), parameter :: rhos_vis_p(0:78) &
      = (/0.310,  0.160,  0.160,  0.160,  0.160,  0.160,  0.160,  0.160,&
          0.160,  0.160,  0.160,  0.160,  0.310,  0.310,  0.310,  0.310,&
          0.310,  0.310,  0.310,  0.310,  0.310,  0.310,  0.310,  0.310,&
          0.310,  0.310,  0.310,  0.310,  0.310,  0.310,  0.310,  0.310,&
          0.310,  0.310,  0.310,  0.310,  0.310,  0.310,  0.310,  0.310,&
          0.310,  0.310,  0.310,  0.310,  0.310,  0.310,  0.310,  0.310,&
          0.310,  0.310,  0.310,  0.310,  0.310,  0.310,  0.310,  0.310,&
          0.310,  0.310,  0.310,  0.310,  0.310,  0.310,  0.310,  0.310,&
          0.310,  0.310,  0.310,  0.310,  0.310,  0.310,  0.310,  0.310,&
          0.310,  0.310,  0.310,  0.310,  0.310,  0.310,  0.310/)

   ! reflectance of green leaf in near infrared band 
   REAL(r8), parameter :: rhol_nir_p(0:78) &
      = (/0.350,  0.350,  0.350,  0.350,  0.450,  0.450,  0.450,  0.450,&
          0.450,  0.350,  0.450,  0.450,  0.350,  0.350,  0.350,  0.350,&
          0.350,  0.350,  0.350,  0.350,  0.350,  0.350,  0.350,  0.350,&
          0.350,  0.350,  0.350,  0.350,  0.350,  0.350,  0.350,  0.350,&
          0.350,  0.350,  0.350,  0.350,  0.350,  0.350,  0.350,  0.350,&
          0.350,  0.350,  0.350,  0.350,  0.350,  0.350,  0.350,  0.350,&
          0.350,  0.350,  0.350,  0.350,  0.350,  0.350,  0.350,  0.350,&
          0.350,  0.350,  0.350,  0.350,  0.350,  0.350,  0.350,  0.350,&
          0.350,  0.350,  0.350,  0.350,  0.350,  0.350,  0.350,  0.350,&
          0.350,  0.350,  0.350,  0.350,  0.350,  0.350,  0.350/)

   ! reflectance of dead leaf in near infrared band 
   REAL(r8), parameter :: rhos_nir_p(0:78) &
      = (/0.530,  0.390,  0.390,  0.390,  0.390,  0.390,  0.390,  0.390,&
          0.390,  0.390,  0.390,  0.390,  0.530,  0.530,  0.530,  0.530,&
          0.530,  0.530,  0.530,  0.530,  0.530,  0.530,  0.530,  0.530,&
          0.530,  0.530,  0.530,  0.530,  0.530,  0.530,  0.530,  0.530,&
          0.530,  0.530,  0.530,  0.530,  0.530,  0.530,  0.530,  0.530,&
          0.530,  0.530,  0.530,  0.530,  0.530,  0.530,  0.530,  0.530,&
          0.530,  0.530,  0.530,  0.530,  0.530,  0.530,  0.530,  0.530,&
          0.530,  0.530,  0.530,  0.530,  0.530,  0.530,  0.530,  0.530,&
          0.530,  0.530,  0.530,  0.530,  0.530,  0.530,  0.530,  0.530,&
          0.530,  0.530,  0.530,  0.530,  0.530,  0.530,  0.530/)

   ! transmittance of green leaf in visible band 
   REAL(r8), parameter :: taul_vis_p(0:78) &
      = (/0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,&
          0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,&
          0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,&
          0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,&
          0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,&
          0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,&
          0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,&
          0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,&
          0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,&
          0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050/)

   ! transmittance of dead leaf in visible band 
   REAL(r8), parameter :: taus_vis_p(0:78) &
      = (/0.120,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,&
          0.001,  0.001,  0.001,  0.001,  0.120,  0.120,  0.120,  0.120,&
          0.120,  0.120,  0.120,  0.120,  0.120,  0.120,  0.120,  0.120,&
          0.120,  0.120,  0.120,  0.120,  0.120,  0.120,  0.120,  0.120,&
          0.120,  0.120,  0.120,  0.120,  0.120,  0.120,  0.120,  0.120,&
          0.120,  0.120,  0.120,  0.120,  0.120,  0.120,  0.120,  0.120,&
          0.120,  0.120,  0.120,  0.120,  0.120,  0.120,  0.120,  0.120,&
          0.120,  0.120,  0.120,  0.120,  0.120,  0.120,  0.120,  0.120,&
          0.120,  0.120,  0.120,  0.120,  0.120,  0.120,  0.120,  0.120,&
          0.120,  0.120,  0.120,  0.120,  0.120,  0.120,  0.120/)

   ! transmittance of green leaf in near infrared band 
   REAL(r8), parameter :: taul_nir_p(0:78) &
      = (/0.340,  0.100,  0.100,  0.100,  0.250,  0.250,  0.250,  0.250,&
          0.250,  0.100,  0.250,  0.250,  0.340,  0.340,  0.340,  0.340,&
          0.340,  0.340,  0.340,  0.340,  0.340,  0.340,  0.340,  0.340,&
          0.340,  0.340,  0.340,  0.340,  0.340,  0.340,  0.340,  0.340,&
          0.340,  0.340,  0.340,  0.340,  0.340,  0.340,  0.340,  0.340,&
          0.340,  0.340,  0.340,  0.340,  0.340,  0.340,  0.340,  0.340,&
          0.340,  0.340,  0.340,  0.340,  0.340,  0.340,  0.340,  0.340,&
          0.340,  0.340,  0.340,  0.340,  0.340,  0.340,  0.340,  0.340,&
          0.340,  0.340,  0.340,  0.340,  0.340,  0.340,  0.340,  0.340,&
          0.340,  0.340,  0.340,  0.340,  0.340,  0.340,  0.340/)

   ! transmittance of dead leaf in near infrared band 
   REAL(r8), parameter :: taus_nir_p(0:78) &
      = (/0.250,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,&
          0.001,  0.001,  0.001,  0.001,  0.250,  0.250,  0.250,  0.250,&
          0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,&
          0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,&
          0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,&
          0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,&
          0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,&
          0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,&
          0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,&
          0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250/)

   ! maximum carboxylation rate at 25 C at canopy top
   ! /06/03/2014/ based on Bonan et al., 2010 (Table 2)
   REAL(r8), parameter :: vmax25_p(0:78) &
!      = (/ 52.0, 61.0, 54.0, 57.0, 72.0, 72.0, 52.0, 52.0,&
      = (/ 52.0, 61.0, 54.0, 57.0, 30.0, 72.0, 52.0, 52.0,&
           52.0, 72.0, 52.0, 52.0, 52.0, 52.0, 52.0, 57.0,&
           57.0, 57.0, 57.0, 57.0, 57.0, 57.0, 57.0, 57.0,&
           57.0, 57.0, 57.0, 57.0, 57.0, 57.0, 57.0, 57.0,&
           57.0, 57.0, 57.0, 57.0, 57.0, 57.0, 57.0, 57.0,&
           57.0, 57.0, 57.0, 57.0, 57.0, 57.0, 57.0, 57.0,&
           57.0, 57.0, 57.0, 57.0, 57.0, 57.0, 57.0, 57.0,&
           57.0, 57.0, 57.0, 57.0, 57.0, 57.0, 57.0, 57.0,&
           57.0, 57.0, 57.0, 57.0, 57.0, 57.0, 57.0, 57.0,&
           57.0, 57.0, 57.0, 57.0, 57.0, 57.0, 57.0/) * 1.e-6

   ! quantum efficiency
   REAL(r8), parameter :: effcon_p(0:78) &
      = (/0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08,&
          0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.05, 0.08,&
          0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08,&
          0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08,&
          0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08,&
          0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08,&
          0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08,&
          0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08,&
          0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08,&
          0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08/)

   ! conductance-photosynthesis slope parameter
   REAL(r8), parameter :: gradm_p(0:78) &
      = (/9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0,&
          9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 4.0, 9.0,&
          9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0,&
          9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0,&
          9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0,&
          9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0,&
          9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0,&
          9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0,&
          9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0,&
          9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0/)

   ! conductance-photosynthesis intercept
   REAL(r8), parameter :: binter_p(0:78) &
      = (/0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,&
          0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.04, 0.01,&
          0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,&
          0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,&
          0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,&
          0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,&
          0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,&
          0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,&
          0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,&
          0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01/)

   ! respiration fraction
   REAL(r8), parameter :: respcp_p(0:78) &
      = (/0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015,&
          0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.025, 0.015,&
          0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015,&
          0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015,&
          0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015,&
          0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015,&
          0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015,&
          0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015,&
          0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015,&
          0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015/)

   ! slope of high temperature inhibition FUNCTION (s1)
   REAL(r8), parameter :: shti_p(0:78) = 0.3

   ! slope of low temperature inhibition FUNCTION (s3) 
   REAL(r8), parameter :: slti_p(0:78) = 0.2

   ! temperature coefficient in gs-a model (s5)
   REAL(r8), parameter :: trda_p(0:78) = 1.3
   
   ! temperature coefficient in gs-a model (s6)
   REAL(r8), parameter :: trdm_p(0:78) = 328.0

   ! temperature coefficient in gs-a model (273.16+25)
   REAL(r8), parameter :: trop_p(0:78) = 298.0

   ! 1/2 point of high temperature inhibition FUNCTION (s2)
   REAL(r8), parameter :: hhti_p(0:78) &
      =(/308.0, 303.0, 303.0, 303.0, 313.0, 313.0, 311.0, 311.0,&
         311.0, 313.0, 313.0, 303.0, 303.0, 308.0, 313.0, 308.0,&
         308.0, 308.0, 308.0, 308.0, 308.0, 308.0, 308.0, 308.0,&
         308.0, 308.0, 308.0, 308.0, 308.0, 308.0, 308.0, 308.0,&
         308.0, 308.0, 308.0, 308.0, 308.0, 308.0, 308.0, 308.0,&
         308.0, 308.0, 308.0, 308.0, 308.0, 308.0, 308.0, 308.0,&
         308.0, 308.0, 308.0, 308.0, 308.0, 308.0, 308.0, 308.0,&
         308.0, 308.0, 308.0, 308.0, 308.0, 308.0, 308.0, 308.0,&
         308.0, 308.0, 308.0, 308.0, 308.0, 308.0, 308.0, 308.0,&
         308.0, 308.0, 308.0, 308.0, 308.0, 308.0, 308.0/)

   ! 1/2 point of low temperature inhibition FUNCTION (s4)
   REAL(r8), parameter :: hlti_p(0:78) &
      =(/281.0, 278.0, 278.0, 278.0, 288.0, 288.0, 283.0, 283.0,&
         283.0, 283.0, 283.0, 278.0, 278.0, 281.0, 288.0, 281.0,&
         281.0, 281.0, 281.0, 281.0, 281.0, 281.0, 281.0, 281.0,&
         281.0, 281.0, 281.0, 281.0, 281.0, 281.0, 281.0, 281.0,&
         281.0, 281.0, 281.0, 281.0, 281.0, 281.0, 281.0, 281.0,&
         281.0, 281.0, 281.0, 281.0, 281.0, 281.0, 281.0, 281.0,&
         281.0, 281.0, 281.0, 281.0, 281.0, 281.0, 281.0, 281.0,&
         281.0, 281.0, 281.0, 281.0, 281.0, 281.0, 281.0, 281.0,&
         281.0, 281.0, 281.0, 281.0, 281.0, 281.0, 281.0, 281.0,&
         281.0, 281.0, 281.0, 281.0, 281.0, 281.0, 281.0/)

   ! coefficient of leaf nitrogen allocation
   REAL(r8), parameter :: extkn_p(0:78) = 0.5

   REAL(r8) :: &
#ifndef CROP
      rho_p(2,2,0:N_PFT-1), &!leaf reflectance
      tau_p(2,2,0:N_PFT-1)   !leaf transmittance
#else
      rho_p(2,2,0:N_PFT+N_CFT-1), &!leaf reflectance
      tau_p(2,2,0:N_PFT+N_CFT-1)   !leaf transmittance
#endif

   ! depth at 50% roots
   REAL(r8), parameter, dimension(0:78) :: d50_p &
      =(/27.0,  21.0,  12.0,  12.0,  15.0,  23.0,  16.0,  23.0,&
         12.0,  23.5,  23.5,  23.5,   9.0,   7.0,  16.0,  22.0,&
         22.0,  22.0,  22.0,  22.0,  22.0,  22.0,  22.0,  22.0,&
         22.0,  22.0,  22.0,  22.0,  22.0,  22.0,  22.0,  22.0,&
         22.0,  22.0,  22.0,  22.0,  22.0,  22.0,  22.0,  22.0,&
         22.0,  22.0,  22.0,  22.0,  22.0,  22.0,  22.0,  22.0,&
         22.0,  22.0,  22.0,  22.0,  22.0,  22.0,  22.0,  22.0,&
         22.0,  22.0,  22.0,  22.0,  22.0,  22.0,  22.0,  22.0,&
         22.0,  22.0,  22.0,  22.0,  22.0,  22.0,  22.0,  22.0,&
         22.0,  22.0,  22.0,  22.0,  22.0,  22.0,  22.0/)

   ! coefficient of root profile
   REAL(r8), parameter, dimension(0:78) :: beta_p &
      =(/-2.051, -1.835, -1.880, -1.880, -1.632, -1.757, -1.681, -1.757,&
         -1.880, -1.623, -1.623, -1.623, -2.621, -1.176, -1.452, -1.796,&
         -1.796, -1.796, -1.796, -1.796, -1.796, -1.796, -1.796, -1.796,&
         -1.796, -1.796, -1.796, -1.796, -1.796, -1.796, -1.796, -1.796,&
         -1.796, -1.796, -1.796, -1.796, -1.796, -1.796, -1.796, -1.796,&
         -1.796, -1.796, -1.796, -1.796, -1.796, -1.796, -1.796, -1.796,&
         -1.796, -1.796, -1.796, -1.796, -1.796, -1.796, -1.796, -1.796,&
         -1.796, -1.796, -1.796, -1.796, -1.796, -1.796, -1.796, -1.796,&
         -1.796, -1.796, -1.796, -1.796, -1.796, -1.796, -1.796, -1.796,&
         -1.796, -1.796, -1.796, -1.796, -1.796, -1.796, -1.796/)

   ! woody (1) or grass (0)
   INTEGER , parameter, dimension(0:78) :: woody &
      =(/0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, &
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/)

   !设定PFT的根分布参数 
   REAL(r8), PRIVATE, parameter :: roota(0:78) &
      =(/  0.0,   7.0,   7.0,   7.0,   7.0,   7.0,   6.0,   6.0,&
           6.0,   7.0,   7.0,   7.0,  11.0,  11.0,  11.0,   6.0,&
           6.0,   6.0,   6.0,   6.0,   6.0,   6.0,   6.0,   6.0,&
           6.0,   6.0,   6.0,   6.0,   6.0,   6.0,   6.0,   6.0,&
           6.0,   6.0,   6.0,   6.0,   6.0,   6.0,   6.0,   6.0,&
           6.0,   6.0,   6.0,   6.0,   6.0,   6.0,   6.0,   6.0,&
           6.0,   6.0,   6.0,   6.0,   6.0,   6.0,   6.0,   6.0,&
           6.0,   6.0,   6.0,   6.0,   6.0,   6.0,   6.0,   6.0,&
           6.0,   6.0,   6.0,   6.0,   6.0,   6.0,   6.0,   6.0,&
           6.0,   6.0,   6.0,   6.0,   6.0,   6.0,   6.0/)
   
   REAL(r8), PRIVATE, parameter :: rootb(0:78) &
      =(/  0.0,   2.0,   2.0,   2.0,   1.0,   1.0,   2.0,   2.0,&
           2.0,   1.5,   1.5,   1.5,   2.0,   2.0,   2.0,   3.0,&
           3.0,   3.0,   3.0,   3.0,   3.0,   3.0,   3.0,   3.0,&
           3.0,   3.0,   3.0,   3.0,   3.0,   3.0,   3.0,   3.0,&
           3.0,   3.0,   3.0,   3.0,   3.0,   3.0,   3.0,   3.0,&
           3.0,   3.0,   3.0,   3.0,   3.0,   3.0,   3.0,   3.0,&
           3.0,   3.0,   3.0,   3.0,   3.0,   3.0,   3.0,   3.0,&
           3.0,   3.0,   3.0,   3.0,   3.0,   3.0,   3.0,   3.0,&
           3.0,   3.0,   3.0,   3.0,   3.0,   3.0,   3.0,   3.0,&
           3.0,   3.0,   3.0,   3.0,   3.0,   3.0,   3.0/)


!   bgc PFT constants

   REAL(r8), parameter, dimension(0:78) :: grperc = 0.11_r8


   REAL(r8), parameter, dimension(0:78) :: grpnow = 1._r8


   REAL(r8), parameter, dimension(0:78) :: lf_flab = 0.25_r8


   REAL(r8), parameter, dimension(0:78) :: lf_fcel = 0.5_r8


   REAL(r8), parameter, dimension(0:78) :: lf_flig = 0.25_r8


   REAL(r8), parameter, dimension(0:78) :: fr_flab = 0.25_r8


   REAL(r8), parameter, dimension(0:78) :: fr_fcel = 0.5_r8


   REAL(r8), parameter, dimension(0:78) :: fr_flig = 0.25_r8


   LOGICAL , parameter, dimension(0:78) :: isshrub & ! True => is a shrub
      =(/.False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .True.,  .True.,  .True.,  .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False./)

   LOGICAL , parameter, dimension(0:78) :: isgrass & ! True => is a grass
      =(/.False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .True.,  .True.,  .True.,  .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False./)

   LOGICAL , parameter, dimension(0:78) :: isbetr  & ! True => is tropical broadleaf evergreen tree
      =(/.False., .False., .False., .False., .True.,  .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False./)

   LOGICAL , parameter, dimension(0:78) :: isbdtr  & ! True => is a broadleaf deciduous tree
      =(/.False., .False., .False., .False., .False., .False., .True.,  .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False./)

   LOGICAL , parameter, dimension(0:78) :: isevg   & ! True => is a evergreen tree
      =(/.False., .True.,  .True.,  .False., .True.,  .True.,  .False., .False., &
         .False., .True.,  .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False./)

   LOGICAL , parameter, dimension(0:78) :: issed   & ! True => is a seasonal deciduous tree
      =(/.False., .False., .False., .True.,  .False., .False., .False., .True.,  &
         .True.,  .False., .False., .True.,  .True.,  .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False./)

   LOGICAL , parameter, dimension(0:78) :: isstd   & ! True => is a stress deciduous tree
      =(/.False., .False., .False., .False., .False., .False., .True.,  .False., &
         .False., .False., .True.,  .False., .False., .True.,  .True.,  .True. , & 
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False./)

   LOGICAL , parameter, dimension(0:78) :: isbare  & ! True => is a bare land
      =(/.True.,  .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False./)

   LOGICAL , parameter, dimension(0:78) :: iscrop  & ! True => is a crop land
      =(/.False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .True. , &
          .True.,  .True.,  .True.,  .True.,  .True.,  .True.,  .True.,  .True., &
          .True.,  .True.,  .True.,  .True.,  .True.,  .True.,  .True.,  .True., &
          .True.,  .True.,  .True.,  .True.,  .True.,  .True.,  .True.,  .True., &
          .True.,  .True.,  .True.,  .True.,  .True.,  .True.,  .True.,  .True., &
          .True.,  .True.,  .True.,  .True.,  .True.,  .True.,  .True.,  .True., &
          .True.,  .True.,  .True.,  .True.,  .True.,  .True.,  .True.,  .True., &
          .True.,  .True.,  .True.,  .True.,  .True.,  .True.,  .True.,  .True., &
          .True.,  .True.,  .True.,  .True.,  .True.,  .True.,  .True./)

   LOGICAL , parameter, dimension(0:78) :: isnatveg &! True => is a natural vegetation
      =(/.False., .True.,  .True.,  .True.,  .True.,  .True.,  .True.,  .True., &
         .True.,  .True.,  .True.,  .True.,  .True.,  .True.,  .True.,  .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False., .False., &
         .False., .False., .False., .False., .False., .False., .False./)

   REAL(r8), parameter, dimension(0:78) :: fsr_pft &
      =(/   0.,   0.26,   0.26,   0.26,   0.25,   0.25,   0.25,   0.25, &
          0.25,   0.28,   0.28,   0.28,   0.33,   0.33,   0.33,   0.33, &
          0.33,   0.33,   0.33,   0.33,   0.33,   0.33,   0.33,   0.33, &
          0.33,   0.33,   0.33,   0.33,   0.33,   0.33,   0.33,   0.33, &
          0.33,   0.33,   0.33,   0.33,   0.33,   0.33,   0.33,   0.33, &
          0.33,   0.33,   0.33,   0.33,   0.33,   0.33,   0.33,   0.33, &
          0.33,   0.33,   0.33,   0.33,   0.33,   0.33,   0.33,   0.33, &
          0.33,   0.33,   0.33,   0.33,   0.33,   0.33,   0.33,   0.33, &
          0.33,   0.33,   0.33,   0.33,   0.33,   0.33,   0.33,   0.33, &
          0.33,   0.33,   0.33,   0.33,   0.33,   0.33,   0.33/)

   REAL(r8), parameter, dimension(0:78) :: fd_pft &
      =(/   0.,     24.,     24.,     24.,     24.,     24.,     24.,     24., &
           24.,     24.,     24.,     24.,     24.,     24.,     24.,     24., &
           24.,      0.,      0.,      0.,      0.,      0.,      0.,      0., &
            0.,      0.,      0.,      0.,      0.,      0.,      0.,      0., &
            0.,      0.,      0.,      0.,      0.,      0.,      0.,      0., &
            0.,      0.,      0.,      0.,      0.,      0.,      0.,      0., &
            0.,      0.,      0.,      0.,      0.,      0.,      0.,      0., &
            0.,      0.,      0.,      0.,      0.,      0.,      0.,      0., &
            0.,      0.,      0.,      0.,      0.,      0.,      0.,      0., &
            0.,      0.,      0.,      0.,      0.,      0.,      0./)

   REAL(r8), parameter, dimension(0:78) :: leafcn &
      =(/              1.,              58.,              58., 25.8131130614352, &
          29.603315571344,  29.603315571344, 23.4521575984991, 23.4521575984991, &
         23.4521575984991, 36.4166059723234, 23.2558139534884, 23.2558139534884, &
         28.0269058295964, 28.0269058295964, 35.3606789250354, 28.0269058295964, &
                       25.,               25.,               25.,               20., &
                       20.,               20.,               20.,               20., &
                       20.,               20.,               20.,               20., &
                       20.,               20.,               20.,               20., &
                       20.,               20.,               20.,               20., &
                       20.,               20.,               20.,               20., &
                       20.,               20.,               20.,               20., &
                       20.,               20.,               20.,               20., &
                       20.,               20.,               20.,               20., &
                       20.,               20.,               20.,               20., &
                       20.,               20.,               20.,               20., &
                       20.,               20.,               20.,               20., &
                       20.,               20.,               20.,               25., &
                       25.,               20.,               20.,               20., &
                       20.,               20.,               20.,               25., &
                       25.,               20.,               20./)

   REAL(r8), parameter, dimension(0:78) :: frootcn &
      =(/   1.,     42.,     42.,     42.,     42.,     42.,     42.,     42.,&
           42.,     42.,     42.,     42.,     42.,     42.,     42.,     42.,&
           42.,     42.,     42.,     42.,     42.,     42.,     42.,     42.,&
           42.,     42.,     42.,     42.,     42.,     42.,     42.,     42.,&
           42.,     42.,     42.,     42.,     42.,     42.,     42.,     42.,&
           42.,     42.,     42.,     42.,     42.,     42.,     42.,     42.,&
           42.,     42.,     42.,     42.,     42.,     42.,     42.,     42.,&
           42.,     42.,     42.,     42.,     42.,     42.,     42.,     42.,&
           42.,     42.,     42.,     42.,     42.,     42.,     42.,     42.,&
           42.,     42.,     42.,     42.,     42.,     42.,     42./)

   REAL(r8), parameter, dimension(0:78) :: livewdcn &
      =(/   1.,     50.,     50.,     50.,     50.,     50.,     50.,     50.,&
           50.,     50.,     50.,     50.,      0.,      0.,      0.,      0.,&
            0.,     50.,     50.,     50.,     50.,     50.,     50.,     50.,&
           50.,     50.,     50.,     50.,     50.,     50.,     50.,     50.,&
           50.,     50.,     50.,     50.,     50.,     50.,     50.,     50.,&
           50.,     50.,     50.,     50.,     50.,     50.,     50.,     50.,&
           50.,     50.,     50.,     50.,     50.,     50.,     50.,     50.,&
           50.,     50.,     50.,     50.,     50.,     50.,     50.,     50.,&
           50.,     50.,     50.,     50.,     50.,     50.,     50.,     50.,&
           50.,     50.,     50.,     50.,     50.,     50.,     50./)

   REAL(r8), parameter, dimension(0:78) :: deadwdcn &
      =(/   1.,    500.,    500.,    500.,    500.,    500.,    500.,    500.,&
          500.,    500.,    500.,    500.,      0.,      0.,      0.,      0.,& 
            0.,    500.,    500.,    500.,    500.,    500.,    500.,    500.,&
          500.,    500.,    500.,    500.,    500.,    500.,    500.,    500.,&
          500.,    500.,    500.,    500.,    500.,    500.,    500.,    500.,&
          500.,    500.,    500.,    500.,    500.,    500.,    500.,    500.,&
          500.,    500.,    500.,    500.,    500.,    500.,    500.,    500.,&
          500.,    500.,    500.,    500.,    500.,    500.,    500.,    500.,&
          500.,    500.,    500.,    500.,    500.,    500.,    500.,    500.,&
          500.,    500.,    500.,    500.,    500.,    500.,    500./)

   REAL(r8), parameter, dimension(0:78) :: graincn &
      =(/-999.,   -999.,   -999.,   -999.,   -999.,   -999.,   -999.,   -999.,&
         -999.,   -999.,   -999.,   -999.,   -999.,   -999.,   -999.,   -999.,&
         -999.,     50.,     50.,     50.,     50.,     50.,     50.,     50.,&
           50.,     50.,     50.,     50.,     50.,     50.,     50.,     50.,&
           50.,     50.,     50.,     50.,     50.,     50.,     50.,     50.,&
           50.,     50.,     50.,     50.,     50.,     50.,     50.,     50.,&
           50.,     50.,     50.,     50.,     50.,     50.,     50.,     50.,&
           50.,     50.,     50.,     50.,     50.,     50.,     50.,     50.,&
           50.,     50.,     50.,     50.,     50.,     50.,     50.,     50.,&
           50.,     50.,     50.,     50.,     50.,     50.,     50./)

   REAL(r8), parameter, dimension(0:78) :: lflitcn &
      =(/   1.,     70.,     80.,     50.,     60.,     60.,     50.,     50.,&
           50.,     60.,     50.,     50.,     50.,     50.,     50.,     50.,&
           50.,     25.,     25.,     25.,     25.,     25.,     25.,     25.,&
           25.,     25.,     25.,     25.,     25.,     25.,     25.,     25.,&
           25.,     25.,     25.,     25.,     25.,     25.,     25.,     25.,&
           25.,     25.,     25.,     25.,     25.,     25.,     25.,     25.,&
           25.,     25.,     25.,     25.,     25.,     25.,     25.,     25.,&
           25.,     25.,     25.,     25.,     25.,     25.,     25.,     25.,&
           25.,     25.,     25.,     25.,     25.,     25.,     25.,     25.,&
           25.,     25.,     25.,     25.,     25.,     25.,     25./)

   REAL(r8), parameter, dimension(0:78) :: leaf_long &
      =(/            0., 3.30916666666667, 3.30916666666667, 0.506666666666667,&
                 1.4025,           1.4025, 0.48333333333333, 0.483333333333333,&
      0.483333333333333, 1.32333333333333,             0.39,              0.39,&
      0.320833333333333, 0.32083333333333,             0.14, 0.320833333333333,&
                     1.,               1.,               1.,                1.,&
                     1.,               1.,               1.,                1.,&
                     1.,               1.,               1.,                1.,&
                     1.,               1.,               1.,                1.,&
                     1.,               1.,               1.,                1.,&
                     1.,               1.,               1.,                1.,&
                     1.,               1.,               1.,                1.,&
                     1.,               1.,               1.,                1.,&
                     1.,               1.,               1.,                1.,&
                     1.,               1.,               1.,                1.,&
                     1.,               1.,               1.,                1.,&
                     1.,               1.,               1.,                1.,&
                     1.,               1.,               1.,                1.,&
                     1.,               1.,               1.,                1.,&
                     1.,               1.,               1.,                1.,&
                     1.,               1.,               1./)

   REAL(r8), parameter, dimension(0:78) :: cc_leaf  &
      =(/   0.,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8/)

   REAL(r8), parameter, dimension(0:78) :: cc_lstem &
      =(/   0.,     0.3,     0.3,     0.3,    0.27,    0.27,    0.27,    0.27,&
          0.27,    0.35,    0.35,    0.35,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8/)

   REAL(r8), parameter, dimension(0:78) :: cc_dstem &
      =(/   0.,     0.3,     0.3,     0.3,    0.27,    0.27,    0.27,    0.27,&
          0.27,    0.35,    0.35,    0.35,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8/)

   REAL(r8), parameter, dimension(0:78) :: cc_other &
      =(/   0.,     0.5,     0.5,     0.5,    0.45,    0.45,    0.45,    0.45,&
          0.45,    0.55,    0.55,    0.55,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8/)

   REAL(r8), parameter, dimension(0:78) :: fm_leaf  &
      =(/   0.,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8/)

   REAL(r8), parameter, dimension(0:78) :: fm_lstem &
      =(/   0.,     0.5,     0.5,     0.5,    0.45,    0.45,    0.35,    0.35,&
          0.45,    0.55,    0.55,    0.55,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8/)

   REAL(r8), parameter, dimension(0:78) :: fm_lroot &
      =(/   0.,    0.15,    0.15,    0.15,    0.13,    0.13,     0.1,     0.1,&
          0.13,    0.17,    0.17,    0.17,     0.2,     0.2,     0.2,     0.2,&
           0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,&
           0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,&
           0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,&
           0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,&
           0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,&
           0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,&
           0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,&
           0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2/)

   REAL(r8), parameter, dimension(0:78) :: fm_root  &
      =(/   0.,    0.15,    0.15,    0.15,    0.13,    0.13,     0.1,     0.1,&
          0.13,    0.17,    0.17,    0.17,     0.2,     0.2,     0.2,     0.2,&
           0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,&
           0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,&
           0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,&
           0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,&
           0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,&
           0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,&
           0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,&
           0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2/)

   REAL(r8), parameter, dimension(0:78) :: fm_droot &
      =(/   0.,    0.15,    0.15,    0.15,    0.13,    0.13,     0.1,     0.1,&
          0.13,    0.17,    0.17,    0.17,     0.2,     0.2,     0.2,     0.2,&
           0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,&
           0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,&
           0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,&
           0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,&
           0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,&
           0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,&
           0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2,&
           0.2,     0.2,     0.2,     0.2,     0.2,     0.2,     0.2/)

   REAL(r8), parameter, dimension(0:78) :: fm_other &
      =(/   0.,     0.5,     0.5,     0.5,    0.45,    0.45,    0.35,    0.35,&
          0.45,    0.55,    0.55,    0.55,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8,&
           0.8,     0.8,     0.8,     0.8,     0.8,     0.8,     0.8/)

   REAL(r8), parameter, dimension(0:78) :: froot_leaf         &
      =(/   0.,     1.5,     1.5,     1.5,     1.5,     1.5,     1.5,     1.5,&
           1.5,     1.5,     1.5,     1.5,     1.5,     1.5,     1.5,     1.5,&
            1.,      2.,      2.,      2.,      2.,      2.,      2.,      2.,&
            2.,      2.,      2.,      2.,      2.,      2.,      2.,      2.,&
            2.,      2.,      2.,      2.,      2.,      2.,      2.,      2.,&
            2.,      2.,      2.,      2.,      2.,      2.,      2.,      2.,&
            2.,      2.,      2.,      2.,      2.,      2.,      2.,      2.,&
            2.,      2.,      2.,      2.,      2.,      2.,      2.,      2.,&
            2.,      2.,      2.,      2.,      2.,      2.,      2.,      2.,&
            2.,      2.,      2.,      2.,      2.,      2.,      2./)

   REAL(r8), parameter, dimension(0:78) :: croot_stem         &
      =(/  0.3,     0.3,     0.3,     0.3,     0.3,     0.3,     0.3,     0.3,&
           0.3,     0.3,     0.3,     0.3,      0.,      0.,      0.,      0.,&
            0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.,&
            0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.,&
            0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.,&
            0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.,&
            0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.,&
            0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.,&
            0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.,&
            0.,      0.,      0.,      0.,      0.,      0.,      0./)

   REAL(r8), parameter, dimension(0:78) :: stem_leaf          &
      =(/   0.,     2.3,     2.3,      1.,     2.3,     1.5,      1.,     2.3,&
           2.3,     1.4,    0.24,    0.24,      0.,      0.,      0.,      0.,&
            0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.,&
            0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.,&
            0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.,&
            0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.,&
            0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.,&
            0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.,&
            0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.,&
            0.,      0.,      0.,      0.,      0.,      0.,      0./)

   REAL(r8), parameter, dimension(0:78) :: flivewd            &
      =(/   0.,     0.1,     0.1,     0.1,     0.1,     0.1,     0.1,     0.1,&
           0.1,     0.5,     0.5,     0.1,      0.,      0.,      0.,      0.,&
            0.,      1.,      1.,      1.,      1.,      1.,      1.,      1.,&
            1.,      1.,      1.,      1.,      1.,      1.,      1.,      1.,&
            1.,      1.,      1.,      1.,      1.,      1.,      1.,      1.,&
            1.,      1.,      1.,      1.,      1.,      1.,      1.,      1.,&
            1.,      1.,      1.,      1.,      1.,      1.,      1.,      1.,&
            1.,      1.,      1.,      1.,      1.,      1.,      1.,      1.,&
            1.,      1.,      1.,      1.,      1.,      1.,      1.,      1.,&
            1.,      1.,      1.,      1.,      1.,      1.,      1./)

   REAL(r8), parameter, dimension(0:78) :: fcur2              &
      =(/   0.,      1.,      1.,      0.,      1.,      1.,      0.,      0.,&
            0.,      1.,      0.,      0.,      0.,      0.,      0.,      0.,&
            0.,      1.,      1.,      1.,      1.,      1.,      1.,      1.,&
            1.,      1.,      1.,      1.,      1.,      1.,      1.,      1.,&
            1.,      1.,      1.,      1.,      1.,      1.,      1.,      1.,&
            1.,      1.,      1.,      1.,      1.,      1.,      1.,      1.,&
            1.,      1.,      1.,      1.,      1.,      1.,      1.,      1.,&
            1.,      1.,      1.,      1.,      1.,      1.,      1.,      1.,&
            1.,      1.,      1.,      1.,      1.,      1.,      1.,      1.,&
            1.,      1.,      1.,      1.,      1.,      1.,      1./)

   REAL(r8), parameter, dimension(0:78) :: dsladlai             &
      =(/   0., 0.00125,   0.001,   0.003, 0.00122,  0.0015,  0.0027,  0.0027,&
        0.0027,      0.,      0.,      0.,      0.,      0.,      0.,      0.,&
            0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.,&
            0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.,&
            0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.,&
            0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.,&
            0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.,&
            0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.,&
            0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.,&
            0.,      0.,      0.,      0.,      0.,      0.,      0./)

   REAL(r8), parameter, dimension(0:78) :: slatop             &
      =(/   0.,    0.01,    0.01, 0.02018,   0.019,   0.019,  0.0308,  0.0308,&
        0.0308, 0.01798, 0.03072, 0.03072, 0.04024, 0.04024, 0.03846, 0.04024,&
         0.035,    0.05,    0.05,   0.035,   0.035,   0.035,   0.035,   0.035,&
         0.035,   0.035,   0.035,   0.035,   0.035,   0.035,   0.035,   0.035,&
         0.035,   0.035,   0.035,   0.035,   0.035,   0.035,   0.035,   0.035,&
         0.035,   0.035,   0.035,   0.035,   0.035,   0.035,   0.035,   0.035,&
         0.035,   0.035,   0.035,   0.035,   0.035,   0.035,   0.035,   0.035,&
         0.035,   0.035,   0.035,   0.035,   0.035,   0.035,   0.035,   0.035,&
         0.035,   0.035,   0.035,    0.05,    0.05,   0.035,   0.035,   0.035,&
         0.035,   0.035,   0.035,    0.05,    0.05,   0.035,   0.035/)
!--- crop variables

   REAL(r8), parameter, dimension(0:78) :: minplanttemp       &
     =(/-999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9,&
        -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9,&
        -999.9, 279.15, 279.15, 272.15, 272.15, 278.15, 278.15, 279.15,&
        279.15, 272.15, 272.15, 278.15, 278.15, 272.15, 272.15, 278.15,&
        278.15, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9,&
        -999.9, 283.15, 283.15, -999.9, -999.9, -999.9, -999.9, -999.9,&
        -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9,&
        -999.9, -999.9, -999.9, -999.9, -999.9, 283.15, 283.15, -999.9,&
        -999.9, -999.9, -999.9, 283.15, 283.15, -999.9, -999.9, -999.9,&
        -999.9, -999.9, -999.9, 283.15, 283.15, 283.15, 283.15 /)

   INTEGER, parameter, dimension(0:78,1:2) :: minplantymd   &
     =(/(/-999, -999, -999, -999, -999, -999, -999, -999, & ! north hemisphere
          -999, -999, -999, -999, -999, -999, -999, -999, &
          -999,  401,  401,  401,  401,  901,  901,  501, &
           501,  401,  401,  901,  901,  401,  401,  901, &
           901, -999, -999, -999, -999, -999, -999, -999, &
          -999,  401,  401, -999, -999, -999, -999, -999, &
          -999, -999, -999, -999, -999, -999, -999, -999, &
          -999, -999, -999, -999, -999,  101,  101, -999, &
          -999, -999, -999,  101,  101, -999, -999, -999, &
          -999, -999, -999,  320,  320,  415,  415/)    , &
        (/-999, -999, -999, -999, -999, -999, -999, -999, & ! south hemisphere
          -999, -999, -999, -999, -999, -999, -999, -999, &
          -999, 1001, 1001, 1001, 1001,  301,  301, 1101, &
          1101, 1001, 1001,  301,  301, 1001, 1001,  301, &
           301, -999, -999, -999, -999, -999, -999, -999, &
          -999,  901,  901, -999, -999, -999, -999, -999, &
          -999, -999, -999, -999, -999, -999, -999, -999, &
          -999, -999, -999, -999, -999, 1015, 1015, -999, &
          -999, -999, -999,  801,  801, -999, -999, -999, &
          -999, -999, -999,  920,  920, 1015, 1015/)/)

   INTEGER, dimension(0:78,1:2) :: minplantjday

   INTEGER, parameter, dimension(0:78,1:2) :: maxplantymd   &
     =(/(/-999, -999, -999, -999, -999, -999, -999, -999, & ! north hemisphere
          -999, -999, -999, -999, -999, -999, -999, -999, &
          -999,  615,  615,  615,  615, 1130, 1130,  615, &
           615,  615,  615, 1130, 1130,  615,  615, 1130, &
          1130, -999, -999, -999, -999, -999, -999, -999, &
          -999,  531,  531, -999, -999, -999, -999, -999, &
          -999, -999, -999, -999, -999, -999, -999, -999, &
          -999, -999, -999, -999, -999,  228,  228, -999, &
          -999, -999, -999,  331,  331, -999, -999, -999, &
          -999, -999, -999,  415,  415,  701,  701/)    , &
        (/-999, -999, -999, -999, -999, -999, -999, -999, & ! south hemispehre
          -999, -999, -999, -999, -999, -999, -999, -999, &
          -999, 1215, 1215, 1215, 1215,  530,  530, 1215, &
          1215, 1215, 1215,  530,  530, 1215, 1215,  530, &
           530, -999, -999, -999, -999, -999, -999, -999, &
          -999, 1130, 1130, -999, -999, -999, -999, -999, &
          -999, -999, -999, -999, -999, -999, -999, -999, &
          -999, -999, -999, -999, -999, 1231, 1231, -999, &
          -999, -999, -999, 1031, 1031, -999, -999, -999, &
          -999, -999, -999, 1015, 1015, 1231, 1231/)/)

   INTEGER, dimension(0:78,1:2) :: maxplantjday

    REAL(r8),parameter, dimension(0:78) :: hybgdd &
      = (/-999, -999, -999, -999, -999, -999, -999, -999, &
          -999, -999, -999, -999, -999, -999, -999, -999, &
          -999, 1700, 1700, 1700, 1700, 1700, 1700, 1900, &
          1900, 1700, 1700, 1700, 1700, 1700, 1700, 1700, &
          1700, -999, -999, -999, -999, -999, -999, -999, &
          -999, 1700, 1700, -999, -999, -999, -999, -999, &
          -999, -999, -999, -999, -999, -999, -999, -999, &
          -999, -999, -999, -999, -999, 2100, 2100, -999, &
          -999, -999, -999, 4300, 4300, -999, -999, -999, &
          -999, -999, -999, 1800, 1800, 2100, 2100 /)

    REAL(r8),parameter, dimension(0:78) :: manunitro  &   ! Max fertilizer to be applied in total (kg N/m2)
      = (/  0.,     0.,     0.,     0.,     0.,     0.,     0.,     0., &
            0.,     0.,     0.,     0.,     0.,     0.,     0.,     0., &
            0., 0.0150, 0.0150, 0.0080, 0.0080, 0.0080, 0.0080, 0.0025, &
        0.0025, 0.0080, 0.0080, 0.0080, 0.0080, 0.0080, 0.0080, 0.0080, &
        0.0080,     0.,     0.,     0.,     0.,     0.,     0.,     0., &
            0.,   0.02,   0.02,     0.,     0.,     0.,     0.,     0., &
            0.,     0.,     0.,     0.,     0.,     0.,     0.,     0., &
            0.,     0.,     0.,     0.,     0.,   0.02,   0.02,     0., &
            0.,     0.,     0.,   0.04,   0.04,     0.,     0.,     0., &
            0.,     0.,     0.,   0.03,   0.03,   0.05,   0.05/)

    REAL(r8),parameter, dimension(0:78) :: planttemp  &  ! planting temperature used in CNPhenology (K)
      = (/-999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, 283.15, 283.15, 280.15, 280.15, -999.9, -999.9, 286.15, &
          286.15, 280.15, 280.15, -999.9, -999.9, 280.15, 280.15, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, 294.15, 294.15, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, 294.15, 294.15, -999.9, &
          -999.9, -999.9, -999.9, 294.15, 294.15, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, 294.15, 294.15, 294.15, 294.15/)

    REAL(r8),parameter, dimension(0:78) :: lfemerg   & ! parameter used in CNPhenology
      = (/-999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9,   0.03,   0.03,   0.05,   0.05,   0.03,   0.03,   0.03, &
            0.03,   0.05,   0.05,   0.05,   0.05,   0.05,   0.05,   0.05, &
            0.05, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9,   0.03,   0.03, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9,   0.01,   0.01, -999.9, &
          -999.9, -999.9, -999.9,   0.03,   0.03, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9,   0.03,   0.03,   0.03,   0.03/)

    INTEGER, parameter, dimension(0:78) :: mxmat   & ! parameter used in CNPhenology
      = (/-999, -999, -999, -999, -999, -999, -999, -999, &
          -999, -999, -999, -999, -999, -999, -999, -999, &
          -999,  165,  165,  150,  150,  330,  330,  150, &
           150,  150,  150,  265,  265,  150,  150,  265, &
           265, -999, -999, -999, -999, -999, -999, -999, &
          -999,  160,  160, -999, -999, -999, -999, -999, &
          -999, -999, -999, -999, -999, -999, -999, -999, &
          -999, -999, -999, -999, -999,  150,  150, -999, &
          -999, -999, -999,  300,  300, -999, -999, -999, &
          -999, -999, -999,  160,  160,  150,  150/)

    REAL(r8),parameter, dimension(0:78) :: grnfill  & ! parameter used in CNPhenology
      = (/-999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9,   0.65,   0.65,    0.6,    0.6,    0.4,    0.4,    0.5, &
             0.5,    0.6,    0.6,    0.4,    0.4,    0.6,    0.6,    0.4, &
             0.4, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9,    0.5,    0.5, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9,    0.4,    0.4, -999.9, &
          -999.9, -999.9, -999.9,   0.65,   0.65, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9,    0.5,    0.5,    0.5,    0.5/)

    REAL(r8),parameter, dimension(0:78) :: mxtmp   & ! parmaeter used in accFlds
      = (/-999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9,    30.,    30.,    26.,    26.,    26.,    26.,    30., &
             30.,    26.,    26.,    26.,    26.,    26.,    26.,    26., &
             26., -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9,    30.,    30., -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9,    30.,    30., -999.9, &
          -999.9, -999.9, -999.9,    30.,    30., -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9,    30.,    30.,    30.,    30. /)

    REAL(r8),parameter, dimension(0:78) :: baset   & ! parameter used in accFlds
      = (/0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., &
          0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., &
          0.,  8.,  8.,  0.,  0.,  0.,  0., 10., &
          10., 0.,  0.,  0.,  0.,  0.,  0.,  0., &
          0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., &
          0., 10., 10.,  0.,  0.,  0.,  0.,  0., &
          0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., &
          0.,  0.,  0.,  0.,  0., 10., 10.,  0., &
          0.,  0.,  0., 10., 10.,  0.,  0.,  0., &
          0.,  0.,  0., 10., 10., 10., 10./)     

    REAL(r8),parameter, dimension(0:78) :: gddmin  & ! parameter used in CNPhenology
      = (/-999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9,    50.,    50.,    50.,    50.,    50.,    50.,    50., &
             50.,    50.,    50.,    50.,    50.,    50.,    50.,    50., &
             50.,    50.,    50.,    50.,    50.,    50.,    50.,    50., &
             50.,    50.,    50.,    50.,    50.,    50.,    50.,    50., &
             50.,    50.,    50.,    50.,    50.,    50.,    50.,    50., &
             50.,    50.,    50.,    50.,    50.,    50.,    50.,    50., &
             50.,    50.,    50.,    50.,    50.,    50.,    50.,    50., &
             50.,    50.,    50.,    50.,    50.,    50.,    50./)

    REAL(r8),parameter, dimension(0:78) :: aleaff = 0._r8 ! ! parameter used in CNAllocation

    REAL(r8),parameter, dimension(0:78) :: astemf  & ! parameter used in CNAllocation
      = (/-999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9,   0.05,   0.05,   0.05,   0.05,    0.3, &
             0.3,   0.05,   0.05,   0.05,   0.05,   0.05,   0.05,   0.05, &
            0.05, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9,    0.3,    0.3, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9,   0.05,   0.05, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9,    0.3,    0.3/)

    REAL(r8),parameter, dimension(0:78) :: arooti  & ! parameter used in CNAllocation
      = (/-999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9,    0.1,    0.1,   0.05,   0.05,    0.2,    0.2,    0.2, &
             0.2,    0.3,    0.3,    0.3,    0.3,    0.3,    0.3,    0.3, &
             0.3, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9,    0.2,    0.2, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9,    0.1,    0.1, -999.9, &
          -999.9, -999.9, -999.9,    0.1,    0.1, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9,    0.1,    0.1,    0.2,    0.2/)

    REAL(r8),parameter, dimension(0:78) :: arootf  & ! parameter used in CNAllocation
      = (/-999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9,   0.05,   0.05, -999.9, -999.9, -999.9, -999.9,    0.2, &
             0.2, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9,    0.2,    0.2, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9,   0.05,   0.05, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9,   0.05,   0.05,    0.2,    0.2/)

    REAL(r8),parameter, dimension(0:78) ::fleafi   & ! parameter used in CNAllocation
      = (/-999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9,    0.6,    0.6,    0.9,    0.9,    0.8,    0.8,   0.85, &
            0.85,   0.75,   0.75,  0.425,  0.425,   0.75,   0.75,  0.425, &
           0.425, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9,   0.85,   0.85, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9,   0.75,   0.75, -999.9, &
          -999.9, -999.9, -999.9,    0.6,    0.6, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9,    0.6,    0.6,   0.85,   0.85/)

    REAL(r8),parameter, dimension(0:78) :: bfact   & ! parameter used in CNAllocation
      = (/-999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9,    0.1,    0.1,    0.1,    0.1,    0.1,    0.1,    0.1, &
             0.1,    0.1,    0.1,    0.1,    0.1,    0.1,    0.1,    0.1, &
             0.1,    0.1,    0.1,    0.1,    0.1,    0.1,    0.1,    0.1, &
             0.1,    0.1,    0.1,    0.1,    0.1,    0.1,    0.1,    0.1, &
             0.1,    0.1,    0.1,    0.1,    0.1,    0.1,    0.1,    0.1, &
             0.1,    0.1,    0.1,    0.1,    0.1,    0.1,    0.1,    0.1, &
             0.1,    0.1,    0.1,    0.1,    0.1,    0.1,    0.1,    0.1, &
             0.1,    0.1,    0.1,    0.1,    0.1,    0.1,    0.1/)

    REAL(r8),parameter, dimension(0:78) :: declfact & ! parameter used in CNAllocation
      = (/-999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9,   1.05,   1.05,   1.05,   1.05,   1.05,   1.05,   1.05, &
            1.05,   1.05,   1.05,   1.05,   1.05,   1.05,   1.05,   1.05, &
            1.05,   1.05,   1.05,   1.05,   1.05,   1.05,   1.05,   1.05, &
            1.05,   1.05,   1.05,   1.05,   1.05,   1.05,   1.05,   1.05, &
            1.05,   1.05,   1.05,   1.05,   1.05,   1.05,   1.05,   1.05, &
            1.05,   1.05,   1.05,   1.05,   1.05,   1.05,   1.05,   1.05, &
            1.05,   1.05,   1.05,   1.05,   1.05,   1.05,   1.05,   1.05, &
            1.05,   1.05,   1.05,   1.05,   1.05,   1.05,   1.05/)

    REAL(r8),parameter, dimension(0:78) :: allconss & ! parameter used in CNAllocation
      = (/-999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9,     2.,     2.,     1.,     1.,     1.,     1.,     5., &
              5.,     1.,     1.,     1.,     1.,     1.,     1.,     1., &
              1., -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9,     5.,     5., -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9,     1.,     1., -999.9, &
          -999.9, -999.9, -999.9,     2.,     2., -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9,     2.,     2.,     5.,     5./)

    REAL(r8),parameter, dimension(0:78) :: allconsl & ! parameter used in CNAllocation
      = (/-999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9,     5.,     5.,     3.,     3.,     3.,     3.,     2., &
              2.,     3.,     3.,     3.,     3.,     3.,     3.,     3., &
              3., -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9,     2.,     2., -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9,     3.,     3., -999.9, &
          -999.9, -999.9, -999.9,     5.,     5., -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9,     5.,     5.,     2.,     2./)


    REAL(r8),parameter, dimension(0:78) :: fleafcn & ! C:N during grain fill; leaf
      = (/999., 999., 999., 999., 999., 999., 999., 999., &
          999., 999., 999., 999., 999., 999., 999., 999., &
          999.,  65.,  65.,  65.,  65.,  65.,  65.,  65., &
           65.,  65.,  65.,  65.,  65.,  65.,  65.,  65., &
           65.,  65.,  65.,  65.,  65.,  65.,  65.,  65., &
           65.,  65.,  65.,  65.,  65.,  65.,  65.,  65., &
           65.,  65.,  65.,  65.,  65.,  65.,  65.,  65., &
           65.,  65.,  65.,  65.,  65.,  65.,  65.,  65., &
           65.,  65.,  65.,  65.,  65.,  65.,  65.,  65., &
           65.,  65.,  65.,  65.,  65.,  65.,  65./)

    REAL(r8),parameter, dimension(0:78) :: fstemcn & ! C:N during grain fill; stem
      = (/999., 999., 999., 999., 999., 999., 999., 999., &
          999., 999., 999., 999., 999., 999., 999., 999., &
          999., 120., 120., 100., 100., 100., 100., 130., &
          130., 100., 100., 100., 100., 100., 100., 100., &
          100., 999., 999., 999., 999., 999., 999., 999., &
          999., 130., 130., 999., 999., 999., 999., 999., &
          999., 999., 999., 999., 999., 999., 999., 999., &
          999., 999., 999., 999., 999., 100., 100., 999., &
          999., 999., 999., 120., 120., 999., 999., 999., &
          999., 999., 999., 120., 120., 130., 130./)

    REAL(r8),parameter, dimension(0:78) :: ffrootcn & ! C:N during grain fill; fine root
      = (/999., 999., 999., 999., 999., 999., 999., 999., &
          999., 999., 999., 999., 999., 999., 999., 999., &
          999.,   0.,   0.,  40.,  40.,  40.,  40.,   0., &
            0.,  40.,  40.,  40.,  40.,  40.,  40.,  40., &
           40., 999., 999., 999., 999., 999., 999., 999., &
          999.,   0.,   0., 999., 999., 999., 999., 999., &
          999., 999., 999., 999., 999., 999., 999., 999., &
          999., 999., 999., 999., 999.,  40.,  40., 999., &
          999., 999., 999.,   0.,   0., 999., 999., 999., &
          999., 999., 999.,   0.,   0.,   0.,   0./)

    REAL(r8),parameter, dimension(0:78) :: laimx    & ! maximum leaf area index
      = (/-999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9,     5.,     5.,     7.,     7.,     7.,     7.,     6., &
              6.,     7.,     7.,     7.,     7.,     7.,     7.,     7., &
              7., -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9,     6.,     6., -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9, -999.9, -999.9,     7.,     7., -999.9, &
          -999.9, -999.9, -999.9,     5.,     5., -999.9, -999.9, -999.9, &
          -999.9, -999.9, -999.9,     5.,     5.,     6.,      6./)

    INTEGER, parameter, dimension(0:78) :: mergetoclmpft & ! merge crop functional types
      = (/0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, &
         19, 20, 21, 22, 23, 24, 19, 20, 21, 22, 19, 20, 21, 22, 61, 62, 19, 20, 61, &
         62, 61, 62, 41, 42, 41, 42, 19, 20, 19, 20, 61, 62, 75, 76, 61, 62, 19, 20, &
         19, 20, 19, 20, 61, 62, 75, 76, 19, 20, 67, 68, 19, 20, 75, 76, 75, 76, 75, &
         76, 77, 78/)
!   end bgc variables

   ! scheme 1: Zeng 2001, 2: Schenk and Jackson, 2002
   INTEGER, PRIVATE :: ROOTFR_SCHEME = 1 
   
   !fraction of roots in each soil layer 
#ifdef CROP
   REAL(r8), dimension(nl_soil,N_PFT+N_CFT) :: &
      rootfr_p(1:nl_soil, 0:N_PFT+N_CFT-1)
#else
   REAL(r8), dimension(nl_soil,16) :: &
      rootfr_p(1:nl_soil, 0:15) 
#endif

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
#ifdef CROP
      DO i = 0, N_PFT+N_CFT-1
#else
      DO i = 0, N_PFT-1
#endif
         rootfr_p(1,i)=1./(1.+(z_soih(1)*100./d50_p(i))**beta_p(i)) 
         rootfr_p(nl_soil,i)=1.-1./(1.+(z_soih(nl_soil-1)*100./d50_p(i))**beta_p(i)) 

         DO nsl=2,nl_soil-1
            rootfr_p(nsl,i)=1./(1.+(z_soih(nsl)*100./d50_p(i))**beta_p(i)) &
               -1./(1.+(z_soih(nsl-1)*100./d50_p(i))**beta_p(i))
         ENDDO
      ENDDO 
ELSE
      ! PFT rootfr_p (Zeng, 2001)
#ifdef CROP
      DO i = 0, N_PFT+N_CFT-1
#else
      DO i = 0, N_PFT-1
#endif
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

#ifdef CROP
   DO i = 0, N_PFT+N_CFT-1
      minplantjday(i,1) = get_calday(minplantymd(i,1),.false.)
      minplantjday(i,2) = get_calday(minplantymd(i,2),.false.)
      maxplantjday(i,1) = get_calday(maxplantymd(i,1),.false.)
      maxplantjday(i,2) = get_calday(maxplantymd(i,2),.false.)
   END DO
#endif
   
   END SUBROUTINE Init_PFT_Const

END MODULE PFT_Const
