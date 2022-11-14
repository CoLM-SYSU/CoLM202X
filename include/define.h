! 1. Spatial structure: 
!    Select one of the following options.
#define GRIDBASED
#undef CATCHMENT 
#undef UNSTRUCTURED
#undef SinglePoint

! 2. Land TYPE classification : 
!    Select one of the following options.
#undef USGS_CLASSIFICATION       
#undef IGBP_CLASSIFICATION       
#define PFT_CLASSIFICATION       
#undef PC_CLASSIFICATION       

! 3. If defined, ocean area is excluded.
#define LANDONLY

! 4. If defined, there is only one patch in each cell,
!    instead of mosaics of different land types.
#undef USE_DOMINANT_PATCHTYPE

! 5. If defined, Surface DATA with dimensions [patch,lon,lat]
!    can be output by subroutines. 
#define MAP_PATCH_TO_GRID
!    Conflict: Only used when defined GRIDBASED.
#ifndef GRIDBASED
#undef MAP_PATCH_TO_GRID
#endif

! 6. If defined, variable histories are output in vectors.
#define HISTORY_IN_VECTOR
!    Conflict: Only used when defined CATCHMENT or UNSTRUCTURED.
#if (!defined CATCHMENT && !defined UNSTRUCTURED)
#undef HISTORY_IN_VECTOR
#endif

! 7. If defined, debug information is output.
#define CLMDEBUG                  

! 8. If defined, MPI parallelization is enabled.
#define  USEMPI
!    Conflict: not used when defined SingPoint.
#if (defined SinglePoint)
#undef USEMPI
#endif

! 9. Hydrological process options.
! 9.1 Two soil hydraulic models can be used.   
#define  Campbell_SOIL_MODEL
#undef vanGenuchten_Mualem_SOIL_MODEL
! 9.2 If defined, variably saturated flow solver is used.
#define VARIABLY_SATURATED_FLOW
! 9.3 If defined, bedrock is defined as lower boudaries of soil water.
#undef USE_DEPTH_TO_BEDROCK
!    Conflicts : 
#ifndef VARIABLY_SATURATED_FLOW
#undef  USE_DEPTH_TO_BEDROCK
#define Campbell_SOIL_MODEL
#undef  vanGenuchten_Mualem_SOIL_MODEL
#endif

! 10. If defined, soil temperature, wetness and snow depth
!     are initialized by input data.
#undef SOILINI                   

! 11. Soil reflectance can be predefined values or load from files.
#undef SOIL_REFL_GUESSED         
#define SOIL_REFL_READ            

! 12. Soil parameter options:
! 12.1 If defined, soil parameters are upscaled from rawdata (1 km grids) to model pixels through
!      FIT algorithm (Montzka et al., 2017), otherwise use median of rawdata within a model pixel.
#define SOILPAR_UPS_FIT
! 12.2 Options for soil thermal conductivity schemes.
!      THERMAL_CONDUCTIVITY_SCHEME_1: Farouki (1981)
!      THERMAL_CONDUCTIVITY_SCHEME_2: Johansen(1975)
!      THERMAL_CONDUCTIVITY_SCHEME_3: Cote and Konrad (2005)
!      THERMAL_CONDUCTIVITY_SCHEME_4: Balland and Arp (2005) (default)
!      THERMAL_CONDUCTIVITY_SCHEME_5: Lu et al. (2007)
!      THERMAL_CONDUCTIVITY_SCHEME_6: Tarnawski and Leong (2012)
!      THERMAL_CONDUCTIVITY_SCHEME_7: De Vries (1963)
#define THERMAL_CONDUCTIVITY_SCHEME_4

! 13. If defined, plant hydraulic scheme is used
#define PLANT_HYDRAULIC_STRESS


! 14. If defined, CaMa-Flood model will be used.
#undef CaMa_Flood

! 15. If defined, BGC model is used.
#define BGC
!    Conflicts :  only used when PFT_CLASSIFICATION is defined.
#ifndef PFT_CLASSIFICATION
#undef BGC
#endif
! 15.1 If defined, LAI is prognostically calculated from leaf carbon and specific leaf area
#define LAIfdbk
!    Conflicts : only used when BGC is defined
#ifndef BGC
#undef LAIfdbk
#endif
! 15.2 If defined, CROP model is used
#define CROP
!    Conflicts : only used when BGC is defined
#ifndef BGC
#undef CROP
#endif
! 15.3 If defined, Semi-Analytic-Spin-Up (SASU) is used
#undef SASU
!    Conflicts : only used when BGC is defined
#ifndef BGC
#undef SASU
#endif
! 15.4 If defined, Fertlization on crop is used
#define FERT
!    Conflicts : only used when CROP is defined
#ifndef CROP
#undef FERT
#endif
! 15.5 If defined, Nitrification-Denitrification is used
#define NITRIF
!    Conflicts : only used when BGC is defined
#ifndef BGC
#undef NITRIF
#endif

! 16 If defined, OzoneStress on plant physiology is used
#define OzoneStress

! 17 if defined, dynamic phenology is used, not used now.
#undef DYN_PHENOLOGY

