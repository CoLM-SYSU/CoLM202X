#include <define.h>

MODULE LC_Const
! -------------------------------
! Created by Hua Yuan, 08/2019
! -------------------------------

   USE precision
   USE GlobalVars

   IMPLICIT NONE
   SAVE

! GLCC USGS Land Use/Land Cover System Legend 
!---------------------------
! 0  Ocean
! 1  Urban and Built-Up Land
! 2  Dryland Cropland and Pasture
! 3  Irrigated Cropland and Pasture
! 4  Mixed Dryland/Irrigated Cropland and Pasture
! 5  Cropland/Grassland Mosaic
! 6  Cropland/Woodland Mosaic
! 7  Grassland
! 8  Shrubland
! 9  Mixed Shrubland/Grassland
!10  Savanna
!11  Deciduous Broadleaf Forest 
!12  Deciduous Needleleaf Forest 
!13  Evergreen Broadleaf Forest
!14  Evergreen Needleleaf Forest
!15  Mixed Forest
!16  Inland Water
!17  Herbaceous Wetland
!18  Wooded Wetland
!19  Barren or Sparsely Vegetated
!20  Herbaceous Tundra
!21  Wooded Tundra
!22  Mixed Tundra
!23  Bare Ground Tundra
!24  Snow or Ice

   ! land water types 
   ! 0: soil, 1: urban, 2: wetland, 3: ice, 4: lake
   INTEGER , parameter, dimension(24) :: patchtypes_usgs &
      = (/1, 0, 0, 0, 0, 0, 0, 0,&
          0, 0, 0, 0, 0, 0, 0, 4,&
          2, 2, 0, 0, 0, 0, 0, 3 /)

   ! default canopy top height
   !NOTE: woody wetland 35m?
   ! shrub land 0.5m? grass like land 1m? all set to 0.5
   REAL(r8), parameter, dimension(24) :: htop0_usgs &
      !=(/ 1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   0.5,&
      !    0.5,   1.0,  20.0,  17.0,  35.0,  17.0,  20.0,   1.0,&
      !    1.0,  35.0,   0.5,   1.0,   1.0,   1.0,   1.0,   1.0/)
      =(/ 1.0,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,&
          0.5,   0.5,  20.0,  17.0,  35.0,  17.0,  20.0,   0.5,&
          0.5,  17.0,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5/)

   ! default canopy bottom height
   REAL(r8), parameter, dimension(24) :: hbot0_usgs &
! 01/06/2020, yuan: adjust htop: grass/shrub -> 0, tree->1
      !=(/0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,   0.1,&
      !    0.1,   0.1,  11.5,   8.5,   1.0,   8.5,  10.0,   0.1,&
      !    0.1,   1.0,   0.1,  0.01,  0.01,  0.01,  0.01,   0.01/)
      =(/ 0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,&
          0.0,   0.0,   1.0,   1.0,   1.0,   1.0,   1.0,   0.0,&
          0.0,   1.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0/)

   ! defulat vegetation fractional cover
   REAL(r8), parameter, dimension(24) :: fveg0_usgs &
      = 1.0 !(/.../)

   ! default stem area index
   ! NOW read from nc file
   REAL(r8), parameter, dimension(24) :: sai0_usgs &
      !=(/0.2, 0.2, 0.3, 0.3, 0.5, 0.5, 1.0, 0.5,&
      !   1.0, 0.5, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0,&
      !   2.0, 2.0, 0.0, 0.1, 0.1, 0.1, 0.0, 0.0/)
      =(/0.2, 0.2, 0.3, 0.3, 0.5, 0.5, 1.0, 0.5,&
         1.0, 0.5, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0,&
         0.2, 2.0, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0/)

   ! ratio to calculate roughness length z0m
   REAL(r8), parameter, dimension(24) :: z0mr_usgs = 0.1 
   
   ! ratio to calculate displacement height d 
   REAL(r8), parameter, dimension(24) :: displar_usgs = 0.667

   ! inverse sqrt of leaf dimension [m**-0.5, m=4 cm]
   REAL(r8), parameter, dimension(24) :: sqrtdi_usgs = 5.0
      
   ! leaf angle distribution parameter
   REAL(r8), parameter, dimension(24) :: chil_usgs &
      = (/-0.300, -0.300, -0.300, -0.300, -0.300, -0.300, -0.300,  0.010,&
           0.010, -0.300,  0.250,  0.010,  0.100,  0.010,  0.125, -0.300,&
          -0.300,  0.100,  0.010, -0.300, -0.300, -0.300, -0.300, -0.300/)

   ! reflectance of green leaf in virsible band 
   REAL(r8), parameter, dimension(24) :: rhol_vis_usgs &
      = (/0.105,  0.105,  0.105,  0.105,  0.105,  0.105,  0.105,  0.100,&
          0.100,  0.105,  0.100,  0.070,  0.100,  0.070,  0.070,  0.105,&
          0.105,  0.100,  0.100,  0.105,  0.105,  0.105,  0.105,  0.105/)

   ! reflectance of dead leaf in virsible band 
   REAL(r8), parameter, dimension(24) :: rhos_vis_usgs &
      = (/0.360,  0.360,  0.360,  0.360,  0.360,  0.360,  0.360,  0.160,&
          0.160,  0.360,  0.160,  0.160,  0.160,  0.160,  0.160,  0.360,&
          0.360,  0.160,  0.160,  0.360,  0.360,  0.360,  0.360,  0.360/)

   ! reflectance of green leaf in near infrared band 
   REAL(r8), parameter, dimension(24) :: rhol_nir_usgs &
      = (/0.580,  0.580,  0.580,  0.580,  0.580,  0.580,  0.580,  0.450,&
          0.450,  0.580,  0.450,  0.350,  0.450,  0.350,  0.400,  0.580,&
          0.580,  0.450,  0.450,  0.580,  0.580,  0.580,  0.580,  0.580/)

   ! reflectance of dead leaf in near infrared band 
   REAL(r8), parameter, dimension(24) :: rhos_nir_usgs &
      = (/0.580,  0.580,  0.580,  0.580,  0.580,  0.580,  0.580,  0.390,&
          0.390,  0.580,  0.390,  0.390,  0.390,  0.390,  0.390,  0.580,&
          0.580,  0.390,  0.390,  0.580,  0.580,  0.580,  0.580,  0.580/)

   ! transmittance of green leaf in visible band 
   REAL(r8), parameter, dimension(24) :: taul_vis_usgs &
      = (/0.070,  0.070,  0.070,  0.070,  0.070,  0.070,  0.070,  0.070,&
          0.070,  0.070,  0.050,  0.050,  0.050,  0.050,  0.050,  0.070,&
          0.070,  0.050,  0.070,  0.070,  0.070,  0.070,  0.070,  0.070/)

   ! transmittance of dead leaf in visible band 
   REAL(r8), parameter, dimension(24) :: taus_vis_usgs &
      = (/0.220,  0.220,  0.220,  0.220,  0.220,  0.220,  0.220,  0.001,&
          0.001,  0.220,  0.001,  0.001,  0.001,  0.001,  0.001,  0.220,&
          0.220,  0.001,  0.001,  0.220,  0.220,  0.220,  0.220,  0.220/)

   ! transmittance of green leaf in near infrared band 
   REAL(r8), parameter, dimension(24) :: taul_nir_usgs &
      = (/0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,&
          0.250,  0.250,  0.250,  0.100,  0.250,  0.100,  0.150,  0.250,&
          0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250/)

   ! transmittance of dead leaf in near infrared band 
   REAL(r8), parameter, dimension(24) :: taus_nir_usgs &
      = (/0.380,  0.380,  0.380,  0.380,  0.380,  0.380,  0.380,  0.001,&
          0.001,  0.380,  0.001,  0.001,  0.001,  0.001,  0.001,  0.380,&
          0.380,  0.001,  0.001,  0.380,  0.380,  0.380,  0.380,  0.380/)

   ! maximum carboxylation rate at 25 C at canopy top
   ! /06/03/2014/ based on Bonan et al., 2010 (Table 2)
   REAL(r8), parameter, dimension(24) :: vmax25_usgs &
      = (/100.0, 57.0, 57.0, 57.0, 52.0, 52.0, 52.0, 52.0,&
           52.0, 52.0, 52.0, 57.0, 72.0, 54.0, 52.0, 57.0,&
           52.0, 52.0, 52.0, 52.0, 52.0, 52.0, 52.0, 52.0/)

   ! quantum efficiency
   !TODO: no C4, 0.05 may have problem
   REAL(r8), parameter, dimension(24) :: effcon_usgs &
      = (/0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08,&
          0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08,&
          0.08, 0.08, 0.08, 0.05, 0.05, 0.05, 0.05, 0.05/)

   ! conductance-photosynthesis slope parameter
   !TODO: no C4, 4.0 may have problem
   REAL(r8), parameter, dimension(24) :: gradm_usgs &
      = (/9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0,&
          9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0,&
          9.0, 9.0, 9.0, 4.0, 4.0, 4.0, 4.0, 4.0/)

   ! conductance-photosynthesis intercept
   REAL(r8), parameter, dimension(24) :: binter_usgs &
      = (/0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,&
          0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,&
          0.01, 0.01, 0.01, 0.04, 0.04, 0.04, 0.04, 0.04/)

   ! respiration fraction
   REAL(r8), parameter, dimension(24) :: respcp_usgs &
      = (/0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015,&
          0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015,&
          0.015, 0.015, 0.015, 0.025, 0.025, 0.025, 0.025, 0.025/)

   ! slope of high temperature inhibition FUNCTION (s1)
   REAL(r8), parameter, dimension(24) :: shti_usgs = 0.3

   ! slope of low temperature inhibition FUNCTION (s3) 
   REAL(r8), parameter, dimension(24) :: slti_usgs = 0.2

   ! temperature coefficient in gs-a model (s5)
   REAL(r8), parameter, dimension(24) :: trda_usgs = 1.3
   
   ! temperature coefficient in gs-a model (s6)
   REAL(r8), parameter, dimension(24) :: trdm_usgs = 328.0

   ! temperature coefficient in gs-a model (273.16+25)
   REAL(r8), parameter, dimension(24) :: trop_usgs = 298.0

   ! 1/2 point of high temperature inhibition FUNCTION (s2)
   REAL(r8), parameter, dimension(24) :: hhti_usgs &
      =(/308.0, 308.0, 308.0, 308.0, 308.0, 308.0, 308.0, 313.0,&
         313.0, 308.0, 311.0, 303.0, 313.0, 303.0, 307.0, 308.0,&
         308.0, 313.0, 313.0, 313.0, 313.0, 313.0, 313.0, 308.0/)

   ! 1/2 point of low temperature inhibition FUNCTION (s4)
   REAL(r8), parameter, dimension(24) :: hlti_usgs &
      =(/281.0, 281.0, 281.0, 281.0, 281.0, 281.0, 281.0, 283.0,&
         283.0, 281.0, 283.0, 278.0, 288.0, 278.0, 281.0, 281.0,&
         281.0, 288.0, 283.0, 288.0, 288.0, 288.0, 288.0, 281.0/)

   ! coefficient of leaf nitrogen allocation
   REAL(r8), parameter, dimension(24) :: extkn_usgs = 0.5

   ! depth at 50% roots
   REAL(r8), parameter, dimension(24) :: d50_usgs &
      =(/23.0,  21.0,  23.0,  22.0,  15.7,  19.0,   9.3,  47.0,&
         28.2,  21.7,  16.0,  16.0,  15.0,  15.0,  15.5,   1.0,&
          9.3,  15.5,  27.0,   9.0,   9.0,   9.0,   9.0,   1.0/)

   ! coefficient of root profile
   REAL(r8), parameter, dimension(24) :: beta_usgs &
      =(/-1.757, -1.835, -1.757, -1.796, -1.577, -1.738, -1.359, -3.245,&
         -2.302, -1.654, -1.681, -1.681, -1.632, -1.632, -1.656, -1.000,&
         -1.359, -1.656, -2.051, -2.621, -2.621, -2.621, -2.621, -1.000/)

   ! Table 2. Zeng, 2001
   ! urban ==> cropland
   ! water/glacier ==> grass
   REAL(r8), parameter, dimension(24) :: roota_usgs &
      =(/ 5.558,  5.558,  5.558,  5.558,  8.149,  5.558, 10.740,  7.022,&
          8.881,  7.920,  5.990,  7.066,  7.344,  7.706,  4.453, 10.740,&
         10.740,  4.453,  8.992,  8.992,  8.992,  8.992,  4.372, 10.740/)

   REAL(r8), parameter, dimension(24) :: rootb_usgs &
      =(/ 2.614,  2.614,  2.614,  2.614,  2.611,  2.614,  2.608,  1.415,&
          2.012,  1.964,  1.955,  1.953,  1.303,  2.175,  1.631,  2.608,&
          2.608,  1.631,  8.992,  8.992,  8.992,  8.992,  0.978,  2.608/)


! MODIS IGBP Land Use/Land Cover System Legend 
!---------------------------
! 0  Ocean
! 1  Evergreen Needleleaf Forests  
! 2  Evergreen Broadleaf Forests 
! 3  Deciduous Needleleaf Forests 
! 4  Deciduous Broadleaf Forests
! 5  Mixed Forests 
! 6  Closed Shrublands 
! 7  Open Shrublands 
! 8  Woody Savannas 
! 9  Savannas 
!10  Grasslands 
!11  Permanent Wetlands 
!12  Croplands 
!13  Urban and Built-up Lands 
!14  Cropland/Natural Vegetation Mosaics 
!15  Permanent Snow and Ice 
!16  Barren 
!17  Water Bodies 

   ! land water types 
   ! 0: soil, 1: urban, 2: wetland, 3: ice, 4: lake
   INTEGER , parameter, dimension(24) :: patchtypes_igbp &
      = (/0, 0, 0, 0, 0, 0, 0, 0,&
          0, 0, 2, 0, 1, 0, 3, 0,&
          4, 0, 0, 0, 0, 0, 0, 0 /)

   ! canopy top height
   REAL(r8), parameter, dimension(24) :: htop0_igbp &
      !=(/17.0,  35.0,  17.0,  20.0,  20.0,   0.5,   0.5,   1.0,&
      !    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,&
      !    1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0/)
      =(/17.0,  35.0,  17.0,  20.0,  20.0,   0.5,   0.5,   0.5,&
          0.5,   0.5,   0.5,   0.5,   1.0,   0.5,   0.5,   0.5,&
          0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5/)

   ! canopy bottom height
   REAL(r8), parameter, dimension(24) :: hbot0_igbp &
! 01/06/2020, yuan: adjust htop: grass/shrub -> 0, tree->1
      !=(/ 8.5,   1.0,   8.5,  11.5,  10.0,   0.1,   0.1,   0.1,&
      !   0.1,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,&
      !  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01/)
      =(/ 1.0,   1.0,   1.0,   1.0,   1.0,   0.0,   0.0,   0.0,&
          0.1,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,&
          0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0/)
   
   ! defulat vegetation fractional cover
   REAL(r8), parameter, dimension(24) :: fveg0_igbp &
      = 1.0 !(/.../)

   ! default stem area index
   REAL(r8), parameter, dimension(24) :: sai0_igbp &
      =(/2.0, 2.0, 2.0, 2.0, 2.0, 0.5, 0.5, 0.5,&
         0.5, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0,&
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)

   ! ratio to calculate roughness length z0m
   REAL(r8), parameter, dimension(24) :: z0mr_igbp = 0.1 
   
   ! ratio to calculate displacement height d 
   REAL(r8), parameter, dimension(24) :: displar_igbp = 0.667

   ! inverse&sqrt leaf specific dimension size 4 cm
   REAL(r8), parameter, dimension(24) :: sqrtdi_igbp = 5.0
      
   ! leaf angle distribution parameter
   REAL(r8), parameter, dimension(24) :: chil_igbp &
      = (/ 0.010,  0.100,  0.010,  0.250,  0.125,  0.010,  0.010,  0.010,&
           0.010, -0.300,  0.100, -0.300,  0.010, -0.300,  0.010,  0.010,&
           0.010,  0.010,  0.010,  0.010,  0.010,  0.010,  0.010,  0.010/)

   ! reflectance of green leaf in virsible band 
   REAL(r8), parameter, dimension(24) :: rhol_vis_igbp &
      = (/0.070,  0.100,  0.070,  0.100,  0.070,  0.105,  0.105,  0.105,&
          0.105,  0.105,  0.105,  0.105,  0.105,  0.105,  0.105,  0.105,&
          0.105,  0.105,  0.105,  0.105,  0.105,  0.105,  0.105,  0.105/)

   ! reflectance of dead leaf in virsible band 
   REAL(r8), parameter, dimension(24) :: rhos_vis_igbp &
      = (/0.160,  0.160,  0.160,  0.160,  0.160,  0.160,  0.160,  0.160,&
          0.160,  0.360,  0.160,  0.360,  0.160,  0.360,  0.160,  0.160,&
          0.160,  0.160,  0.160,  0.160,  0.160,  0.160,  0.160,  0.160/)

   ! reflectance of green leaf in near infrared band 
   REAL(r8), parameter, dimension(24) :: rhol_nir_igbp &
      = (/0.350,  0.450,  0.350,  0.450,  0.400,  0.450,  0.450,  0.580,&
          0.580,  0.580,  0.450,  0.580,  0.450,  0.580,  0.450,  0.450,&
          0.580,  0.580,  0.580,  0.580,  0.580,  0.580,  0.580,  0.580/)

   ! reflectance of dead leaf in near infrared band 
   REAL(r8), parameter, dimension(24) :: rhos_nir_igbp &
      = (/0.390,  0.390,  0.390,  0.390,  0.390,  0.390,  0.390,  0.390,&
          0.390,  0.580,  0.390,  0.580,  0.390,  0.580,  0.390,  0.390,&
          0.580,  0.580,  0.580,  0.580,  0.580,  0.580,  0.580,  0.580/)

   ! transmittance of green leaf in visible band 
   REAL(r8), parameter, dimension(24) :: taul_vis_igbp &
      = (/0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,&
          0.050,  0.070,  0.050,  0.070,  0.050,  0.070,  0.050,  0.050,&
          0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050,  0.050/)

   ! transmittance of dead leaf in visible band 
   REAL(r8), parameter, dimension(24) :: taus_vis_igbp &
      = (/0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,&
          0.001,  0.220,  0.001,  0.220,  0.001,  0.220,  0.001,  0.001,&
          0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001/)

   ! transmittance of green leaf in near infrared band 
   REAL(r8), parameter, dimension(24) :: taul_nir_igbp &
      = (/0.100,  0.250,  0.100,  0.250,  0.150,  0.250,  0.250,  0.250,&
          0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,&
          0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250/)

   ! transmittance of dead leaf in near infrared band 
   REAL(r8), parameter, dimension(24) :: taus_nir_igbp &
      = (/0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,&
          0.001,  0.380,  0.001,  0.380,  0.001,  0.380,  0.001,  0.001,&
          0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001,  0.001/)

   ! maximum carboxylation rate at 25 C at canopy top
   ! /06/03/2014/ based on Bonan et al., 2010 (Table 2)
   REAL(r8), parameter, dimension(24) :: vmax25_igbp &
      = (/ 54.0, 72.0, 57.0, 52.0, 52.0, 52.0, 52.0, 52.0,&
           52.0, 52.0, 52.0, 57.0,100.0, 57.0, 52.0, 52.0,&
           52.0, 52.0, 52.0, 52.0, 52.0, 52.0, 52.0, 52.0/)

   ! quantum efficiency
   !TODO: no C4 
   REAL(r8), parameter, dimension(24) :: effcon_igbp &
      = (/0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08,&
          0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08,&
          0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08/)

   ! conductance-photosynthesis slope parameter
   REAL(r8), parameter, dimension(24) :: gradm_igbp &
      = (/9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0,&
          9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0,&
          9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0, 9.0/)

   ! conductance-photosynthesis intercept
   REAL(r8), parameter, dimension(24) :: binter_igbp &
      = (/0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,&
          0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,&
          0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01/)

   ! respiration fraction
   REAL(r8), parameter, dimension(24) :: respcp_igbp &
      = (/0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015,&
          0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015,&
          0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015, 0.015/)

   ! slope of high temperature inhibition FUNCTION (s1)
   REAL(r8), parameter, dimension(24) :: shti_igbp = 0.3

   ! slope of low temperature inhibition FUNCTION (s3) 
   REAL(r8), parameter, dimension(24) :: slti_igbp = 0.2

   ! temperature coefficient in gs-a model (s5)
   REAL(r8), parameter, dimension(24) :: trda_igbp = 1.3
   
   ! temperature coefficient in gs-a model (s6)
   REAL(r8), parameter, dimension(24) :: trdm_igbp = 328.0

   ! temperature coefficient in gs-a model (273.16+25)
   REAL(r8), parameter, dimension(24) :: trop_igbp = 298.0

   ! 1/2 point of high temperature inhibition FUNCTION (s2)
   REAL(r8), parameter, dimension(24) :: hhti_igbp &
      =(/303.0, 313.0, 303.0, 311.0, 307.0, 308.0, 313.0, 313.0,&
         313.0, 308.0, 313.0, 308.0, 308.0, 308.0, 303.0, 313.0,&
         308.0, 313.0, 313.0, 313.0, 313.0, 313.0, 313.0, 313.0/)

   ! 1/2 point of low temperature inhibition FUNCTION (s4)
   REAL(r8), parameter, dimension(24) :: hlti_igbp &
      =(/278.0, 288.0, 278.0, 283.0, 281.0, 281.0, 288.0, 288.0,&
         288.0, 281.0, 283.0, 281.0, 281.0, 281.0, 278.0, 288.0,&
         281.0, 288.0, 288.0, 288.0, 288.0, 288.0, 288.0, 288.0/)

   ! coefficient of leaf nitrogen allocation
   REAL(r8), parameter, dimension(24) :: extkn_igbp = 0.5

   ! depth at 50% roots
   REAL(r8), parameter, dimension(24) :: d50_igbp &
      =(/15.0,  15.0,  16.0,  16.0,  15.5,  19.0,  28.0,  18.5,&
         28.0,   9.0,   9.0,  22.0,  23.0,  22.0,   1.0,   9.0,&
          1.0,   9.0,   9.0,   9.0,   9.0,   9.0,   9.0,   9.0/)

   ! coefficient of root profile
   REAL(r8), parameter, dimension(24) :: beta_igbp &
      =(/-1.623, -1.623, -1.681, -1.681, -1.652, -1.336, -1.909, -1.582,&
         -1.798, -1.359, -1.359, -1.796, -1.757, -1.796, -1.000, -2.261,&
         -1.000, -2.621, -2.621, -2.621, -2.621, -2.621, -2.621, -2.621/)
   
   ! Table 2. Zeng, 2001
   ! water/glacier ==> grass
   ! urban ==> cropland
   REAL(r8), parameter, dimension(24) :: roota_igbp &
      =(/ 6.706,  7.344,  7.066,  5.990,  4.453,  6.326,  7.718,  7.604,&
          8.235, 10.740, 10.740,  5.558,  5.558,  5.558, 10.740,  4.372,&
         10.740, 10.740, 10.740, 10.740, 10.740, 10.740, 10.740, 10.740/)

   REAL(r8), parameter, dimension(24) :: rootb_igbp &
      =(/ 2.175,  1.303,  1.953,  1.955,  1.631,  1.567,  1.262,  2.300,&
          1.627,  2.608,  2.608,  2.614,  2.614,  2.614,  2.608,  0.978,&
          2.608,  2.608,  2.608,  2.608,  2.608,  2.608,  2.608,  2.608/)

   REAL(r8), dimension(24) :: &
      patchtypes, &! land water types
      htop0,      &! canopy top height
      hbot0,      &! canopy bottom height
      fveg0,      &! canopy vegetation fractional cover
      sai0,       &! canopy stem area index
      chil,       &! leaf angle distribution factor
      z0mr,       &! ratio to calculate roughness length z0m
      displar,    &! ratio to calculate displacement height d 
      sqrtdi,     &! inverse sqrt of leaf dimension [m**-0.5]

      vmax25,     &! maximum carboxylation rate at 25 C at canopy top
      effcon,     &! quantum efficiency
      gradm,      &! conductance-photosynthesis slope parameter
      binter,     &! conductance-photosynthesis intercept
      respcp,     &! respiration fraction
      shti,       &! slope of high temperature inhibition function (s1)
      slti,       &! slope of low temperature inhibition function (s3)
      trda,       &! temperature coefficient in gs-a model (s5)
      trdm,       &! temperature coefficient in gs-a model (s6)
      trop,       &! temperature coefficient in gs-a model (273.16+25)
      hhti,       &! 1/2 point of high temperature inhibition function (s2)
      hlti,       &! 1/2 point of low temperature inhibition function (s4)
      extkn,      &! coefficient of leaf nitrogen allocation

      d50,        &! depth at 50% roots
      beta         ! coefficient of root profile
      
   REAL(r8), PRIVATE, dimension(24) :: &
      roota,      &! root fraction para
      rootb        ! root fraction para
   
   REAL(r8) ::    &
      rho(2,2,24),&! leaf reflectance
      tau(2,2,24)  ! leaf transmittance

   ! scheme 1: Schenk and Jackson, 2002, 2: Zeng 2001
   INTEGER, PRIVATE :: ROOTFR_SCHEME = 1 

   ! fraction of roots in each soil layer 
   REAL(r8), dimension(nl_soil,24) :: rootfr(1:nl_soil, 1:24)   

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: Init_LC_Const

CONTAINS

   SUBROUTINE Init_LC_Const

      IMPLICIT NONE

      INTEGER :: i, nsl

#ifdef USGS_CLASSIFICATION
      patchtypes (:) = patchtypes_usgs (:)
      htop0      (:) = htop0_usgs      (:)
      hbot0      (:) = hbot0_usgs      (:)
      fveg0      (:) = fveg0_usgs      (:)
      sai0       (:) = sai0_usgs       (:)
      z0mr       (:) = z0mr_usgs       (:)
      displar    (:) = displar_usgs    (:)
      sqrtdi     (:) = sqrtdi_usgs     (:)
      chil       (:) = chil_usgs       (:)
      vmax25     (:) = vmax25_usgs     (:) * 1.e-6
      effcon     (:) = effcon_usgs     (:)
      gradm      (:) = gradm_usgs      (:)
      binter     (:) = binter_usgs     (:)
      respcp     (:) = respcp_usgs     (:)
      shti       (:) = shti_usgs       (:)
      slti       (:) = slti_usgs       (:)
      trda       (:) = trda_usgs       (:)
      trdm       (:) = trdm_usgs       (:)
      trop       (:) = trop_usgs       (:)
      hhti       (:) = hhti_usgs       (:)
      hlti       (:) = hlti_usgs       (:)
      extkn      (:) = extkn_usgs      (:)
      d50        (:) = d50_usgs        (:)
      beta       (:) = beta_usgs       (:)
      roota      (:) = roota_usgs      (:)
      rootb      (:) = rootb_usgs      (:)
      rho    (1,1,:) = rhol_vis_usgs   (:)
      rho    (2,1,:) = rhol_nir_usgs   (:)
      rho    (1,2,:) = rhos_vis_usgs   (:)
      rho    (2,2,:) = rhos_nir_usgs   (:)
      tau    (1,1,:) = taul_vis_usgs   (:)
      tau    (2,1,:) = taul_nir_usgs   (:)
      tau    (1,2,:) = taus_vis_usgs   (:)
      tau    (2,2,:) = taus_nir_usgs   (:)
#else 
      patchtypes (:) = patchtypes_igbp (:)
      htop0      (:) = htop0_igbp      (:)
      hbot0      (:) = hbot0_igbp      (:)
      fveg0      (:) = fveg0_igbp      (:)
      sai0       (:) = sai0_igbp       (:)
      z0mr       (:) = z0mr_igbp       (:)
      displar    (:) = displar_igbp    (:)
      sqrtdi     (:) = sqrtdi_igbp     (:)
      chil       (:) = chil_igbp       (:)
      vmax25     (:) = vmax25_igbp     (:) * 1.e-6
      effcon     (:) = effcon_igbp     (:)
      gradm      (:) = gradm_igbp      (:)
      binter     (:) = binter_igbp     (:)
      respcp     (:) = respcp_igbp     (:)
      shti       (:) = shti_igbp       (:)
      slti       (:) = slti_igbp       (:)
      trda       (:) = trda_igbp       (:)
      trdm       (:) = trdm_igbp       (:)
      trop       (:) = trop_igbp       (:)
      hhti       (:) = hhti_igbp       (:)
      hlti       (:) = hlti_igbp       (:)
      extkn      (:) = extkn_igbp      (:)
      d50        (:) = d50_igbp        (:)
      beta       (:) = beta_igbp       (:)
      roota      (:) = roota_igbp      (:)
      rootb      (:) = rootb_igbp      (:)
      rho    (1,1,:) = rhol_vis_igbp   (:)
      rho    (2,1,:) = rhol_nir_igbp   (:)
      rho    (1,2,:) = rhos_vis_igbp   (:)
      rho    (2,2,:) = rhos_nir_igbp   (:)
      tau    (1,1,:) = taul_vis_igbp   (:)
      tau    (2,1,:) = taul_nir_igbp   (:)
      tau    (1,2,:) = taus_vis_igbp   (:)
      tau    (2,2,:) = taus_nir_igbp   (:)
#endif

      ! ----------------------------------------------------------
      ! The definition of global root distribution is based on
      ! Schenk and Jackson, 2002: The Global Biogeography of Roots.
      ! Ecological Monagraph 72(3): 311-328.
      ! ----------------------------------------------------------
      IF (ROOTFR_SCHEME == 1) THEN
         DO i = 1, N_land_classification
            rootfr(1,i)=1./(1.+(z_soih(1)*100./d50(i))**beta(i)) 
            rootfr(nl_soil,i)=1.-1./(1.+(z_soih(nl_soil-1)*100./d50(i))**beta(i)) 

            DO nsl=2,nl_soil-1
               rootfr(nsl,i)=1./(1.+(z_soih(nsl)*100./d50(i))**beta(i)) &
                  -1./(1.+(z_soih(nsl-1)*100./d50(i))**beta(i))
            ENDDO
         ENDDO 
      ELSE 
         DO i = 1, N_land_classification
            rootfr(1,i) = 1. - 0.5*( &
                 exp(-roota(i) * z_soih(1)) &
               + exp(-rootb(i) * z_soih(1)) )
               
            rootfr(nl_soil,i) = 0.5*( &
                 exp(-roota(i) * z_soih(nl_soil)) &
               + exp(-rootb(i) * z_soih(nl_soil)) )
            
            DO nsl = 2, nl_soil-1
               rootfr(nsl,i) = 0.5*( &
                    exp(-roota(i) * z_soih(nsl-1)) &
                  + exp(-rootb(i) * z_soih(nsl-1)) &
                  - exp(-roota(i) * z_soih(nsl)) &
                  - exp(-rootb(i) * z_soih(nsl)) )
            ENDDO
         ENDDO
      ENDIF

   END SUBROUTINE Init_LC_Const

END MODULE LC_Const
! ---------- EOP ------------
