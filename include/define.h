! 1. Spatial structure:
!    Select one of the following options.
#define GRIDBASED
#undef CATCHMENT
#undef UNSTRUCTURED
#undef SinglePoint

! 2. Land TYPE classification :
!    Select one of the following options.
#undef USGS_CLASSIFICATION       
#define IGBP_CLASSIFICATION       
#undef PFT_CLASSIFICATION       
#undef PC_CLASSIFICATION       

! 2.1 Urban model setting (put it temporarily here):
#define URBAN_MODEL
#undef URBAN_LCZ

! 3. If defined, debug information is output.
#undef CoLMDEBUG
! 3.1 If defined, surface data in vector is mapped to gridded data for checking.
#undef SrfdataDiag

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
! 6.2 If defined, lateral flow is modeled.
#define  LATERAL_FLOW
!    Conflicts :
#ifndef CATCHMENT
#undef LATERAL_FLOW
#endif

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
!      THERMAL_CONDUCTIVITY_SCHEME_8: Yan Hengnian, He Hailong et al. (2019)
#define THERMAL_CONDUCTIVITY_SCHEME_4

! 11. If defined, CaMa-Flood model will be used.
#undef CaMa_Flood

! 12. If defined, BGC model is used.
#define BGC
!    Conflicts :  only used when PFT_CLASSIFICATION is defined.
#ifndef PFT_CLASSIFICATION
#undef BGC
#endif
! 12.1 If defined, CROP model is used
#define CROP
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
! 12.4 If defined, Fertlization on crop is used
#define FERT
!    Conflicts : only used when CROP is defined
#ifndef CROP
#undef FERT
#endif
! 12.5 If defined, Nitrification-Denitrification is used
#define NITRIF
!    Conflicts : only used when BGC is defined
#ifndef BGC
#undef NITRIF
#endif

! 13 If defined, Fire is on
#undef Fire
!    Conflicts : only used when BGC is defined
#ifndef BGC
#undef Fire
#endif

! 14 If defined, OzoneStress on plant physiology is used
#undef OzoneStress
!    Conflicts : only used when PFT_CLASSIFICATION is used
#ifndef PFT_CLASSIFICATION
#undef OzoneStress
#endif
! 14.1 If defined, Ozone Data is used instead of constant ozone concentration
#undef OzoneData
!    Conflicts : only used when OzoneStress is defined
#ifndef OzoneStress
#undef OzoneData
#endif
! 15 If defined, SNICAR is on
#undef   SNICAR
! 16 If defined, ... need some one to finish here
#undef   Forcing_Downscaling
#define option_precipitation_adjust_II
#define option_longwave_adjust_II
#define option_precip_phase_discrimination_II
! 17. If defined, diagnostics in wue model will be output
#undef WUEdiag
! 18. If defined, supercooled soil water is implemented, Niu & Yang (2006)
#undef supercool_water
!#ifdef BGC
!#define supercool_water
!#endif
