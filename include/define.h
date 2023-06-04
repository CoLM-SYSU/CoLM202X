! 1. Spatial structure:
!    Select one of the following options.
#define GRIDBASED
#undef CATCHMENT
#undef UNSTRUCTURED
#undef SinglePoint

! 2. Land TYPE classification :
!    Select one of the following options.
#undef LULC_USGS
#define LULC_IGBP
#undef LULC_IGBP_PFT
#undef LULC_IGBP_PC

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

! 5. Hydrological process options.
! 5.1 Two soil hydraulic models can be used.
#define  Campbell_SOIL_MODEL
#undef   vanGenuchten_Mualem_SOIL_MODEL
! 5.2 If defined, lateral flow is modeled.
#define  LATERAL_FLOW
!    Conflicts :
#ifndef CATCHMENT
#undef LATERAL_FLOW
#endif

! 6. Soil reflectance can be predefined values or load from files.
#undef SOIL_REFL_GUESSED
#define SOIL_REFL_READ

! 7. Soil parameter options:
! 7.1 If defined, soil parameters are upscaled from rawdata (1 km grids) to model pixels through
!      FIT algorithm (Montzka et al., 2017), otherwise use median of rawdata within a model pixel.
#undef SOILPAR_UPS_FIT
! 7.2 Options for soil thermal conductivity schemes.
!      THERMAL_CONDUCTIVITY_SCHEME_1: Farouki (1981)
!      THERMAL_CONDUCTIVITY_SCHEME_2: Johansen(1975)
!      THERMAL_CONDUCTIVITY_SCHEME_3: Cote and Konrad (2005)
!      THERMAL_CONDUCTIVITY_SCHEME_4: Balland and Arp (2005) (default)
!      THERMAL_CONDUCTIVITY_SCHEME_5: Lu et al. (2007)
!      THERMAL_CONDUCTIVITY_SCHEME_6: Tarnawski and Leong (2012)
!      THERMAL_CONDUCTIVITY_SCHEME_7: De Vries (1963)
!      THERMAL_CONDUCTIVITY_SCHEME_8: Yan Hengnian, He Hailong et al. (2019)
#define THERMAL_CONDUCTIVITY_SCHEME_4

! 8. If defined, CaMa-Flood model will be used.
#undef CaMa_Flood

! 9. If defined, BGC model is used.
#define BGC
!    Conflicts :  only used when LULC_IGBP_PFT is defined.
#ifndef LULC_IGBP_PFT
#undef BGC
#endif
! 9.1 If defined, CROP model is used
#define CROP
!    Conflicts : only used when BGC is defined
#ifndef BGC
#undef CROP
#endif
! 9.3 If defined, Semi-Analytic-Spin-Up (SASU) is used
#undef SASU
!    Conflicts : only used when BGC is defined
#ifndef BGC
#undef SASU
#endif
! 9.4 If defined, Fertlization on crop is used
#define FERT
!    Conflicts : only used when CROP is defined
#ifndef CROP
#undef FERT
#endif
! 9.5 If defined, Nitrification-Denitrification is used
#define NITRIF
!    Conflicts : only used when BGC is defined
#ifndef BGC
#undef NITRIF
#endif

! 10 If defined, Fire is on
#undef Fire
!    Conflicts : only used when BGC is defined
#ifndef BGC
#undef Fire
#endif

! 11 If defined, SNICAR is on
#undef   SNICAR
! 12. If defined, diagnostics in wue model will be output
#undef WUEdiag
! 13. If defined, supercooled soil water is implemented, Niu & Yang (2006)
#undef supercool_water
!#ifdef BGC
!#define supercool_water
!#endif

! 14. If defined, open Land use and land cover change mode.
#undef LULCC
