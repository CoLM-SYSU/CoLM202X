! 1. Spatial structure: 
!    Select one of the following options.
#define GRIDBASED
#undef CATCHMENT 
#undef UNSTRUCTURED
#undef SinglePoint

! 2. Land TYPE classification : 
!    Select one of the following options.
#define USGS_CLASSIFICATION       
#undef IGBP_CLASSIFICATION       
#undef PFT_CLASSIFICATION       
#undef PC_CLASSIFICATION       

! 3. If defined, debug information is output.
#define CLMDEBUG                  

! 4. If defined, MPI parallelization is enabled.
#define  USEMPI
!    Conflict: not used when defined SingPoint.
#if (defined SinglePoint)
#undef USEMPI
#endif

! 5. If defined, depth to bedrock data is included.
#undef USE_DEPTH_TO_BEDROCK

! 6. Hydrological process options.
! 6.1 Two soil hydraulic models can be used.   
#define  Campbell_SOIL_MODEL
#undef   vanGenuchten_Mualem_SOIL_MODEL

! 7. If defined, soil temperature, wetness and snow depth
!     are initialized by input data.
#undef SOILINI                   

! 8. Soil reflectance can be predefined values or load from files.
#undef SOIL_REFL_GUESSED         
#define SOIL_REFL_READ            

! 9. Soil parameter options:
! 9.1 If defined, soil parameters are upscaled from rawdata (1 km grids) to model pixels through
!      FIT algorithm (Montzka et al., 2017), otherwise use median of rawdata within a model pixel.
#undef SOILPAR_UPS_FIT
! 9.2 Options for soil thermal conductivity schemes.
!      THERMAL_CONDUCTIVITY_SCHEME_1: Farouki (1981)
!      THERMAL_CONDUCTIVITY_SCHEME_2: Johansen(1975)
!      THERMAL_CONDUCTIVITY_SCHEME_3: Cote and Konrad (2005)
!      THERMAL_CONDUCTIVITY_SCHEME_4: Balland and Arp (2005) (default)
!      THERMAL_CONDUCTIVITY_SCHEME_5: Lu et al. (2007)
!      THERMAL_CONDUCTIVITY_SCHEME_6: Tarnawski and Leong (2012)
!      THERMAL_CONDUCTIVITY_SCHEME_7: De Vries (1963)
#define THERMAL_CONDUCTIVITY_SCHEME_4

! 10. If defined, plant hydraulic scheme is used
#define PLANT_HYDRAULIC_STRESS

#undef  DYN_PHENOLOGY             

! 11. If defined, CaMa-Flood model will be used.
#define CaMa_Flood

! 12. If defined, BGC model is used.
#undef BGC
!    Conflicts :  only used when PFT_CLASSIFICATION is defined.
#ifndef PFT_CLASSIFICATION
#undef BGC
#endif
! 12.1 If defined, LAI is prognostically calculated from leaf carbon and specific leaf area
#undef LAIfdbk
!    Conflicts : only used when BGC is defined
#ifndef BGC
#undef LAIfdbk
#endif
! 12.2 If defined, CROP model is used
#undef CROP
!    Conflicts : only used when BGC is defined
#ifndef BGC
#undef CROP
#endif
! 12.3 If defined, Semi-Analytic-Spin-Up (SASU) is used
#undef SASU
!    Conflicts : only used when BGC is defined
#ifndef BGC
#undef SASU
#endif
