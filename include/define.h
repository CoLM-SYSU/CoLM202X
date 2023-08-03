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
#undef URBAN_MODEL
#undef URBAN_LCZ

! 3. If defined, debug information is output.
#define CoLMDEBUG
! 3.1 If defined, range of variables is checked.
#define RangeCheck
! 3.1 If defined, surface data in vector is mapped to gridded data for checking.
#undef SrfdataDiag

! 4. If defined, MPI parallelization is enabled.
#define USEMPI
!    Conflict: not used when defined SingPoint.
#if (defined SinglePoint)
#undef USEMPI
#endif

! 5. Hydrological process options.
! 5.1 Two soil hydraulic models can be used.
#undef   Campbell_SOIL_MODEL
#define  vanGenuchten_Mualem_SOIL_MODEL
! 5.2 If defined, lateral flow is modeled.
#define LATERAL_FLOW
!    Conflicts :
#ifndef CATCHMENT
#undef LATERAL_FLOW
#endif

! 6. Soil reflectance can be predefined values or load from files.
! Soil reflectance is now a namelist DEF_SOIL_REFL_SCHEME
! The below will be removed later
!#undef SOIL_REFL_GUESSED
!#define SOIL_REFL_READ

! 7. If defined, CaMa-Flood model will be used.
#undef CaMa_Flood

! 8. If defined, BGC model is used.
#define BGC
!    Conflicts :  only used when LULC_IGBP_PFT is defined.
#ifndef LULC_IGBP_PFT
#undef BGC
#endif
! 8.1 If defined, CROP model is used
#define CROP
!    Conflicts : only used when BGC is defined
#ifndef BGC
#undef CROP
#endif
! 8.2 If defined, Semi-Analytic-Spin-Up (SASU) is used
!#undef SASU
!    Conflicts : only used when BGC is defined
!#ifndef BGC
!#undef SASU
!#endif
!SASU switch has been moved to namelist: DEF_USE_SASU
! 8.3 If defined, Fertlization on crop is used
!@#define FERT
!    Conflicts : only used when CROP is defined
!#ifndef CROP
!#undef FERT
!#endif
!FERT has been moved to namelist: DEF_USE_FERT
! 8.4 If defined, Nitrification-Denitrification is used
!#define NITRIF
!    Conflicts : only used when BGC is defined
!#ifndef BGC
!#undef NITRIF
!#endif
!NITRIF switch has been moved to namelist: DEF_USE_NITRIF

!! 9 If defined, Fire is on
!#undef Fire
!    Conflicts : only used when BGC is defined
!#ifndef BGC
!#undef Fire
!#endif
!FIRE switch has been moved to namelist: DEF_USE_FIRE

! 10 If defined, SNICAR is on
! NOTE: SNICAR is now a namelist DEF_USE_SNICAR
! This macro will be removed later
!#undef SNICAR

! 11. If defined, diagnostics in wue model will be output
!#undef WUEdiag

! 12. If defined, open Land use and land cover change mode.
#undef LULCC
