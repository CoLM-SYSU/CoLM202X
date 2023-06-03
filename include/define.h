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

! 8. Soil reflectance can be predefined values or load from files.
#undef SOIL_REFL_GUESSED
#define SOIL_REFL_READ


! 10. If defined, plant hydraulic scheme is used
#define PLANT_HYDRAULIC_STRESS


! 11. If defined, CaMa-Flood model will be used.
#undef CaMa_Flood

! 12. If defined, BGC model is used.
#define BGC
!    Conflicts :  only used when LULC_IGBP_PFT is defined.
#ifndef LULC_IGBP_PFT
#undef BGC
#endif
! 12.1 If defined, LAI is prognostically calculated from leaf carbon and specific leaf area
#define LAIfdbk
!    Conflicts : only used when BGC is defined
#ifndef BGC
#undef LAIfdbk
#endif
! 12.2 If defined, CROP model is used
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
!    Conflicts : only used when LULC_IGBP_PFT is used
#ifndef LULC_IGBP_PFT
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
! 17. If defined, diagnostics in wue model will be output
#undef WUEdiag

! 19. If defined, open Land use and land cover change mode.
#undef LULCC
