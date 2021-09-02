#include <define.h>

SUBROUTINE IniTimeConst (nl_soil , zsoih , ivt   ,patchtype, &
                         z0m     , displa, sqrtdi, effcon, vmax25, slti  , &
                         hlti    , shti  , hhti  , trda  , trdm  , trop  , &
                         gradm   , binter, extkn , chil  , ref   , tran  , &
                         rootfr  ) 

!===========================================================================
! Initialize time invariant model variables
! Original author: Yongjiu Dai, 09/15/1999; 08/30/2002, 02/2014
!===========================================================================
  use precision
  implicit none

!--------------------------- Input 
  integer, INTENT(in) :: nl_soil  !number of model soil layers
  integer, INTENT(in) :: ivt       !index for vegetation type [-]
  real(r8), INTENT(in) :: zsoih(0:nl_soil) !interface level below a zsoi level [m]

!--------------------------- Output
  integer, INTENT(out) :: patchtype  !3=glacier/ice sheet, 4=deep lake, 5=shallow lake)
                                     !land water type (0=soil, 1=urban and built-up, 2=wetland,
  real(r8), INTENT(out) :: &!--- Vegetation static parameters
        z0m              , &!aerodynamic roughness length [m]
        displa           , &!displacement height [m]
        sqrtdi           , &!inverse sqrt of leaf dimension [m**-0.5]
        effcon           , &!quantum efficiency of RuBP regeneration (mol CO2/mol quanta)
        vmax25           , &!maximum carboxylation rate at 25 C at canopy top (mol CO2/m2s)
        shti             , &!slope of high temperature inhibition function     (s1)
        hhti             , &!1/2 point of high temperature inhibition function (s2)
        slti             , &!slope of low temperature inhibition function      (s3)
        hlti             , &!1/2 point of low temperature inhibition function  (s4)
        trda             , &!temperature coefficient in gs-a model             (s5)
        trdm             , &!temperature coefficient in gs-a model             (s6)
        trop             , &!temperature coefficient in gs-a model         (273+25)
        gradm            , &!conductance-photosynthesis slope parameter
        binter           , &!conductance-photosynthesis intercep
        extkn            , &!coefficient of leaf nitrogen allocation
        chil             , &!leaf angle distribution factor
        ref(2,2)         , &!leaf reflectance (iw=iband, il=life and dead)
        tran(2,2)        , &!leaf transmittance (iw=iband, il=life and dead)
        rootfr(1:nl_soil)   !fraction of roots in each soil layer 

!--------------------------- Local variables
  integer i, nsl             !indices

#if (defined USGS_CLASSIFICATION)
!---------------------------
! GLCC USGS Land Use/Land Cover System Legend 
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

  integer, parameter :: N_land_classification = 24 ! GLCC USGS number of land cover category
  real(r8), dimension(N_land_classification) :: &!
  z0m_usgs        , &!roughness
  displa_usgs     , &!zero-plane-distance
  sqrtdi_usgs     , &!inverse sqrt of leaf dimension [m**-0.5]
              
  chil_usgs       , &!leaf angle distribution factor
  ref_s_usgs      , &!leaf reflectance
  ref_sd_usgs     , &!leaf reflectance
  ref_l_usgs      , &!leaf reflectance
  ref_ld_usgs     , &!leaf reflectance
  tran_s_usgs     , &!leaf transmittance
  tran_sd_usgs    , &!leaf transmittance
  tran_l_usgs     , &!leaf transmittance
  tran_ld_usgs    , &!leaf transmittance
 
  vmax0_usgs      , &!maximum carboxylation rate at 25 C at canopy top
  effcon_usgs     , &!quantum efficiency
  gradm_usgs      , &!conductance-photosynthesis slope parameter
  binter_usgs     , &!conductance-photosynthesis intercept
  respcp_usgs     , &!respiration fraction
  shti_usgs       , &!slope of high temperature inhibition function (s1)
  slti_usgs       , &!slope of low temperature inhibition function (s3)
  trda_usgs       , &!temperature coefficient in gs-a model (s5)
  trdm_usgs       , &!temperature coefficient in gs-a model (s6)
  trop_usgs       , &!temperature coefficient in gs-a model (273.16+25)
  hhti_usgs       , &!1/2 point of high temperature inhibition function (s2)
  hlti_usgs       , &!1/2 point of low temperature inhibition function (s4)
  extkn_usgs      , &!coefficient of leaf nitrogen allocation
 
  d50_usgs        , &!depth at 50% roots
  d95_usgs        , &!depth at 95% roots
  beta_usgs          !coefficient of root profile

!----------------------------------------------------------------------
  z0m_usgs    =(/0.100,  0.100,  0.100,  0.100,  0.100,  0.100,  0.100,  0.050,&
                 0.050,  0.100,  2.000,  1.700,  3.500,  1.700,  2.000,  0.100,&
                 0.100,  3.500,  0.050,  0.100,  0.100,  0.100,  0.100,  0.100/)
  displa_usgs =(/0.667,  0.667,  0.667,  0.667,  0.667,  0.667,  0.667,  0.333,&
                 0.333,  0.667, 13.333, 11.333, 23.333, 11.333, 13.333,  0.667,&
                 0.667, 23.333,  0.333,  0.667,  0.667,  0.667,  0.667,  0.667/)
  sqrtdi_usgs(:)=5.0
  chil_usgs  =(/-0.300, -0.300, -0.300, -0.300, -0.300, -0.300, -0.300,  0.010,&
                 0.010, -0.300,  0.250,  0.010,  0.100,  0.010,  0.125, -0.300,&
                -0.300,  0.100,  0.010, -0.300, -0.300, -0.300, -0.300, -0.300/)
  ref_s_usgs  =(/0.105,  0.105,  0.105,  0.105,  0.105,  0.105,  0.105,  0.100,&
                 0.100,  0.105,  0.100,  0.070,  0.100,  0.070,  0.070,  0.105,&
                 0.105,  0.100,  0.100,  0.105,  0.105,  0.105,  0.105,  0.105/)
  ref_sd_usgs =(/0.360,  0.360,  0.360,  0.360,  0.360,  0.360,  0.360,  0.160,&
                 0.160,  0.360,  0.160,  0.160,  0.160,  0.160,  0.160,  0.360,&
                 0.360,  0.160,  0.160,  0.360,  0.360,  0.360,  0.360,  0.360/)
  ref_l_usgs  =(/0.580,  0.580,  0.580,  0.580,  0.580,  0.580,  0.580,  0.450,&
                 0.450,  0.580,  0.450,  0.350,  0.450,  0.350,  0.400,  0.580,&
                 0.580,  0.450,  0.450,  0.580,  0.580,  0.580,  0.580,  0.580/)
  ref_ld_usgs =(/0.580,  0.580,  0.580,  0.580,  0.580,  0.580,  0.580,  0.390,&
                 0.390,  0.580,  0.390,  0.390,  0.390,  0.390,  0.390,  0.580,&
                 0.580,  0.390,  0.390,  0.580,  0.580,  0.580,  0.580,  0.580/)
  tran_s_usgs =(/0.070,  0.070,  0.070,  0.070,  0.070,  0.070,  0.070,  0.070,&
                 0.070,  0.070,  0.050,  0.050,  0.050,  0.050,  0.050,  0.070,&
                 0.070,  0.050,  0.070,  0.070,  0.070,  0.070,  0.070,  0.070/)
  tran_sd_usgs=(/0.220,  0.220,  0.220,  0.220,  0.220,  0.220,  0.220,  0.001,&
                 0.001,  0.220,  0.001,  0.001,  0.001,  0.001,  0.001,  0.220,&
                 0.220,  0.001,  0.001,  0.220,  0.220,  0.220,  0.220,  0.220/)
  tran_l_usgs =(/0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,&
                 0.250,  0.250,  0.250,  0.100,  0.250,  0.100,  0.150,  0.250,&
                 0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250/)
  tran_ld_usgs=(/0.380,  0.380,  0.380,  0.380,  0.380,  0.380,  0.380,  0.001,&
                 0.001,  0.380,  0.001,  0.001,  0.001,  0.001,  0.001,  0.380,&
                 0.380,  0.001,  0.001,  0.380,  0.380,  0.380,  0.380,  0.380/)

! /06/03/2014/
! vmax0_usgs( 1: 7)=100.0;   vmax0_usgs( 8: 9)=60.0  
! vmax0_usgs(10:13)=100.0;   vmax0_usgs   (14)=60.0 
! vmax0_usgs   (15)=80.0;    vmax0_usgs(16:18)=100.0 
! vmax0_usgs   (19)=60.0;    vmax0_usgs(20:24)=30.0  

! /06/03/2014/ based on Bonan et al., 2010 (Table 2)
  vmax0_usgs=(/100.0, 57.0, 57.0, 57.0, 52.0, 52.0, 52.0, 52.0,&
                52.0, 52.0, 52.0, 57.0, 72.0, 54.0, 52.0, 57.0,&
                52.0, 52.0, 52.0, 52.0, 52.0, 52.0, 52.0, 52.0/)

! /06/03/2014/ based on Kattge et al., 2009 from Bonan et al., 2010 (Table 2)
! vmax0_usgs=(/100.0, 100.0, 100.0, 100.0, 78.0, 54.0, 78.0, 54.0,&
!               54.0,  54.0,  58.0,  39.0, 41.0, 62.0, 54.0, 100.0,&
!               78.0,  54.0,  54.0,  78.0, 78.0, 78.0, 78.0, 78.0/) 

  effcon_usgs(1:19)=0.08;  effcon_usgs(20:24)=0.05
  gradm_usgs (1:19)=9.0;   gradm_usgs (20:24)=4.0
  binter_usgs(1:19)=0.01;  binter_usgs(20:24)=0.04
  respcp_usgs(1:19)=0.015; respcp_usgs(20:24)=0.025

  shti_usgs(:)=0.3
  slti_usgs(:)=0.2
  trda_usgs(:)=1.3
  trdm_usgs(:)=328.0
  trop_usgs(:)=298.0

  hhti_usgs=(/308.0, 308.0, 308.0, 308.0, 308.0, 308.0, 308.0, 313.0,&
              313.0, 308.0, 311.0, 303.0, 313.0, 303.0, 307.0, 308.0,&
              308.0, 313.0, 313.0, 313.0, 313.0, 313.0, 313.0, 308.0/)
  hlti_usgs=(/281.0, 281.0, 281.0, 281.0, 281.0, 281.0, 281.0, 283.0,&
              283.0, 281.0, 283.0, 278.0, 288.0, 278.0, 281.0, 281.0,&
              281.0, 288.0, 283.0, 288.0, 288.0, 288.0, 288.0, 281.0/)
  extkn_usgs(:)=0.5

  d50_usgs  =(/23.0,  21.0,  23.0,  22.0,  15.7,  19.0,   9.3,  47.0,&
	       28.2,  21.7,  16.0,  16.0,  15.0,  15.0,  15.5,   1.0,&
	        9.3,  15.5,  27.0,   9.0,   9.0,   9.0,   9.0,   1.0/)
  d95_usgs =(/121.0, 104.0, 121.0, 112.5,  80.8, 103.8,  49.0, 302.0,&
	      175.5,  99.3,  95.0,  95.0,  91.0,  91.0,  93.0,   1.0,&
               49.0,  93.0, 112.0,  29.0,  29.0,  29.0,  29.0,   1.0/)
  beta_usgs=(/-1.757, -1.835, -1.757, -1.796, -1.577, -1.738,&
              -1.359, -3.245, -2.302, -1.654, -1.681, -1.681,&
              -1.632, -1.632, -1.656, -1.000, -1.359, -1.656,&
              -2.051, -2.621, -2.621, -2.621, -2.621, -1.000 /)
#endif

#if (defined IGBP_CLASSIFICATION)

!!! CODING soon !

#endif

#if (defined PFT_CLASSIFICATION)

!!! CODING soon !

#endif

!-----------------------------------------------------------------------
! land water type for land classification
!-----------------------------------------------------------------------
#if (defined USGS_CLASSIFICATION)
         i=ivt
                         patchtype=0  ! soil
      if(i==1)           patchtype=1  ! urban and built-up
      if(i==17.or.i==18) patchtype=2  ! wetland
      if(i==24)          patchtype=3  ! land ice
      if(i==16)          patchtype=4  ! land water bodies
      if(i==0)           patchtype=99 ! ocean
#endif
#if (defined IGBP_CLASSIFICATION)
         i=ivt
                         patchtype=0  ! soil
      if(i==13)          patchtype=1  ! urban and built-up
      if(i==11)          patchtype=2  ! wetland
      if(i==15)          patchtype=3  ! land ice
      if(i==17)          patchtype=4  ! land water bodies
      if(i==0)           patchtype=99 ! ocean
#endif

#if (defined PFT_CLASSIFICATION)

!!! CODING soon !

#endif

!-----------------------------------------------------------------------
! vegetation static parameters
! the values for glacier and lake are assigned arbitrarily (not used)
!-----------------------------------------------------------------------
#if (defined USGS_CLASSIFICATION)
           i = ivt
    if( i > 0 )then             ! land grids
         z0m = z0m_usgs   (i)  
      displa = displa_usgs(i)
      sqrtdi = sqrtdi_usgs(i)  

      vmax25 = vmax0_usgs(i)*1.e-6       
      effcon = effcon_usgs(i)       
      slti   =   slti_usgs(i)       
      hlti   =   hlti_usgs(i)       
      shti   =   shti_usgs(i)       
      hhti   =   hhti_usgs(i)       
      trda   =   trda_usgs(i)       
      trdm   =   trdm_usgs(i)       
      trop   =   trop_usgs(i)       
      gradm  =  gradm_usgs(i)       
      binter = binter_usgs(i)       
      extkn  =  extkn_usgs(i)       
 
      chil = chil_usgs(i)             
      ref(1,1) = ref_s_usgs(i)        
      ref(2,1) = ref_l_usgs(i)        
      ref(1,2) = ref_sd_usgs(i)       
      ref(2,2) = ref_ld_usgs(i)       
      tran(1,1) = tran_s_usgs(i)      
      tran(2,1) = tran_l_usgs(i)      
      tran(1,2) = tran_sd_usgs(i)     
      tran(2,2) = tran_ld_usgs(i)     
 
      ! ----------------------------------------------------------
      ! The definition of global root distribution is based on
      ! Schenk and Jackson, 2002: The Global Biogeography of Roots.
      ! Ecological Monagraph 72(3): 311-328.
      ! ----------------------------------------------------------
      if(patchtype>=3)then  !glacier/ice sheet or land water bodies or ocean
         rootfr(:)=0.
      else
         rootfr(1)=1./(1.+(zsoih(  1)*100./d50_usgs(i))**beta_usgs(i)) 
         rootfr(nl_soil)=1.-1./(1.+(zsoih(nl_soil-1)*100./d50_usgs(i))**beta_usgs(i)) 

         do nsl=2,nl_soil-1
            rootfr(nsl)=1./(1.+(zsoih(nsl  )*100./d50_usgs(i))**beta_usgs(i)) &
                       -1./(1.+(zsoih(nsl-1)*100./d50_usgs(i))**beta_usgs(i))
         enddo
      endif

    else                        ! ocean grids
      z0m       = -1.0e36
      displa    = -1.0e36
      sqrtdi    = -1.0e36
      vmax25    = -1.0e36
      effcon    = -1.0e36
      slti      = -1.0e36
      hlti      = -1.0e36
      shti      = -1.0e36
      hhti      = -1.0e36
      trda      = -1.0e36
      trdm      = -1.0e36
      trop      = -1.0e36
      gradm     = -1.0e36
      binter    = -1.0e36
      extkn     = -1.0e36
      chil      = -1.0e36
      ref(:,:)  = -1.0e36
      tran(:,:) = -1.0e36
      rootfr(:) = -1.0e36
    endif
#endif

#if (defined IGBP_CLASSIFICATION)

!!! COMING soon !

#endif

#if (defined PFT_CLASSIFICATION)

!!! COMING soon !

#endif



END SUBROUTINE IniTimeConst
