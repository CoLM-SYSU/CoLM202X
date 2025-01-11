#include <define.h>

#ifdef SinglePoint
MODULE MOD_SingleSrfdata
!-----------------------------------------------------------------------------------------
! DESCRIPTION:
!
!    This module includes subroutines to read or write surface data for "SinglePoint".
!
! Created by Shupeng Zhang, May 2023
!-----------------------------------------------------------------------------------------

   USE MOD_Precision, only: r8
   USE MOD_Vars_Global
   USE MOD_Const_LC
   USE MOD_Namelist
   USE MOD_SPMD_Task
   IMPLICIT NONE
   SAVE

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   integer,  allocatable :: SITE_pfttyp  (:)
   real(r8), allocatable :: SITE_pctpfts (:)
#endif

#ifdef CROP
   integer,  allocatable :: SITE_croptyp (:)
   real(r8), allocatable :: SITE_pctcrop (:)
#endif

   real(r8) :: SITE_htop
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   real(r8), allocatable :: SITE_htop_pfts (:)
#endif

   real(r8), allocatable :: SITE_LAI_monthly (:,:)
   real(r8), allocatable :: SITE_SAI_monthly (:,:)
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   real(r8), allocatable :: SITE_LAI_pfts_monthly (:,:,:)
   real(r8), allocatable :: SITE_SAI_pfts_monthly (:,:,:)
#endif

   integer,  allocatable :: SITE_LAI_year (:)
   real(r8), allocatable :: SITE_LAI_8day (:,:)

   real(r8) :: SITE_lakedepth = 1.

   real(r8) :: SITE_soil_s_v_alb
   real(r8) :: SITE_soil_d_v_alb
   real(r8) :: SITE_soil_s_n_alb
   real(r8) :: SITE_soil_d_n_alb

   real(r8), allocatable :: SITE_soil_vf_quartz_mineral (:)
   real(r8), allocatable :: SITE_soil_vf_gravels        (:)
   real(r8), allocatable :: SITE_soil_vf_sand           (:)
   real(r8), allocatable :: SITE_soil_vf_om             (:)
   real(r8), allocatable :: SITE_soil_wf_gravels        (:)
   real(r8), allocatable :: SITE_soil_wf_sand           (:)
   real(r8), allocatable :: SITE_soil_OM_density        (:)
   real(r8), allocatable :: SITE_soil_BD_all            (:)
   real(r8), allocatable :: SITE_soil_theta_s           (:)
   real(r8), allocatable :: SITE_soil_k_s               (:)
   real(r8), allocatable :: SITE_soil_csol              (:)
   real(r8), allocatable :: SITE_soil_tksatu            (:)
   real(r8), allocatable :: SITE_soil_tksatf            (:)
   real(r8), allocatable :: SITE_soil_tkdry             (:)
   real(r8), allocatable :: SITE_soil_k_solids          (:)
   real(r8), allocatable :: SITE_soil_psi_s             (:)
   real(r8), allocatable :: SITE_soil_lambda            (:)
   real(r8), allocatable :: SITE_soil_theta_r           (:)
#ifdef vanGenuchten_Mualem_SOIL_MODEL
   real(r8), allocatable :: SITE_soil_alpha_vgm         (:)
   real(r8), allocatable :: SITE_soil_L_vgm             (:)
   real(r8), allocatable :: SITE_soil_n_vgm             (:)
#endif
   real(r8), allocatable :: SITE_soil_BA_alpha          (:)
   real(r8), allocatable :: SITE_soil_BA_beta           (:)

   integer :: SITE_soil_texture

   real(r8) :: SITE_dbedrock = 0.

   real(r8) :: SITE_topography = 0.
   real(r8) :: SITE_topostd    = 0.

   ! topography factors used for downscaling
   real(r8) :: SITE_svf = 0.
   real(r8) :: SITE_cur = 0.
   real(r8), allocatable :: SITE_slp_type(:)
   real(r8), allocatable :: SITE_asp_type(:)
   real(r8), allocatable :: SITE_area_type(:)
   real(r8), allocatable :: SITE_sf_lut(:,:)


   integer , allocatable :: SITE_urbtyp   (:)

   real(r8), allocatable :: SITE_lucyid   (:)

   real(r8), allocatable :: SITE_fveg_urb (:)
   real(r8), allocatable :: SITE_htop_urb (:)
   real(r8), allocatable :: SITE_flake_urb(:)
   real(r8), allocatable :: SITE_froof    (:)
   real(r8), allocatable :: SITE_hroof    (:)
   real(r8), allocatable :: SITE_fgimp    (:)
   real(r8), allocatable :: SITE_fgper    (:)
   real(r8), allocatable :: SITE_hlr      (:)
   real(r8), allocatable :: SITE_lambdaw  (:)
   real(r8), allocatable :: SITE_popden   (:)

   real(r8), allocatable :: SITE_em_roof  (:)
   real(r8), allocatable :: SITE_em_wall  (:)
   real(r8), allocatable :: SITE_em_gimp  (:)
   real(r8), allocatable :: SITE_em_gper  (:)
   real(r8), allocatable :: SITE_t_roommax(:)
   real(r8), allocatable :: SITE_t_roommin(:)
   real(r8), allocatable :: SITE_thickroof(:)
   real(r8), allocatable :: SITE_thickwall(:)

   real(r8), allocatable :: SITE_cv_roof  (:)
   real(r8), allocatable :: SITE_cv_wall  (:)
   real(r8), allocatable :: SITE_cv_gimp  (:)
   real(r8), allocatable :: SITE_tk_roof  (:)
   real(r8), allocatable :: SITE_tk_wall  (:)
   real(r8), allocatable :: SITE_tk_gimp  (:)

   real(r8), allocatable :: SITE_alb_roof (:,:)
   real(r8), allocatable :: SITE_alb_wall (:,:)
   real(r8), allocatable :: SITE_alb_gimp (:,:)
   real(r8), allocatable :: SITE_alb_gper (:,:)

CONTAINS

   ! -----
   SUBROUTINE read_surface_data_single (fsrfdata, mksrfdata)

   USE MOD_TimeManager
   USE MOD_NetCDFSerial
   USE MOD_Namelist
   USE MOD_Utils
   USE MOD_Vars_Global, only : PI
   IMPLICIT NONE

   character(len=*), intent(in) :: fsrfdata
   logical, intent(in) :: mksrfdata

   ! Local Variables
   real(r8) :: lat_in, lon_in
   integer  :: iyear, itime

      IF (ncio_var_exist(fsrfdata, 'latitude')) THEN
         CALL ncio_read_serial (fsrfdata, 'latitude',  lat_in)
         IF (lat_in /= SITE_lat_location) THEN
            write(*,*) 'Warning: Latitude mismatch: ', &
               lat_in, ' in data file and ', SITE_lat_location, 'in namelist.'
         ENDIF
         SITE_lat_location = lat_in
      ENDIF

      IF (ncio_var_exist(fsrfdata, 'longitude')) THEN
         CALL ncio_read_serial (fsrfdata, 'longitude', lon_in)
         IF (lon_in /= SITE_lon_location) THEN
            write(*,*) 'Warning: Longitude mismatch: ', &
               lon_in, ' in data file and ', SITE_lon_location, 'in namelist.'
         ENDIF
         SITE_lon_location = lon_in
      ENDIF

      CALL normalize_longitude (SITE_lon_location)

      IF (USE_SITE_landtype) THEN
         IF (trim(fsrfdata) /= 'null') THEN
#ifdef LULC_USGS
            CALL ncio_read_serial (fsrfdata, 'USGS_classification', SITE_landtype)
#else
            CALL ncio_read_serial (fsrfdata, 'IGBP_classification', SITE_landtype)
#endif
         ENDIF

         IF (SITE_landtype < 0) THEN
            write(*,*) 'Error! Please set namelist SITE_landtype first!'
            CALL CoLM_stop()
         ENDIF
      ENDIF

      DEF_domain%edges = floor(SITE_lat_location)
      DEF_domain%edgen = DEF_domain%edges + 1.0
      DEF_domain%edgew = floor(SITE_lon_location)
      DEF_domain%edgee = DEF_domain%edgew + 1.0

      IF (.not. isgreenwich) THEN
         LocalLongitude = SITE_lon_location
      ENDIF

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      IF ((.not. mksrfdata) .or. USE_SITE_pctpfts) THEN
         CALL ncio_read_serial (fsrfdata, 'pfttyp', SITE_pfttyp )
         ! otherwise, retrieve from database by MOD_LandPFT.F90
      ENDIF
#endif
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      IF ((.not. mksrfdata) .or. USE_SITE_pctpfts) THEN
         CALL ncio_read_serial (fsrfdata, 'pctpfts', SITE_pctpfts)
         ! otherwise, retrieve from database by Aggregation_PercentagesPFT.F90
      ENDIF
#endif

#ifdef CROP
      IF ((.not. mksrfdata) .or. USE_SITE_pctcrop) THEN
         CALL ncio_read_serial (fsrfdata, 'croptyp', SITE_croptyp)
         CALL ncio_read_serial (fsrfdata, 'pctcrop', SITE_pctcrop)
         ! otherwise, retrieve from database by MOD_LandPatch.F90
      ENDIF
#endif

      IF ((.not. mksrfdata) .or. USE_SITE_htop) THEN
         ! otherwise, retrieve from database by Aggregation_ForestHeight.F90
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
         CALL ncio_read_serial (fsrfdata, 'canopy_height_pfts', SITE_htop_pfts)
#else
         CALL ncio_read_serial (fsrfdata, 'canopy_height', SITE_htop)
#endif
      ENDIF

      IF ((.not. mksrfdata) .or. USE_SITE_LAI) THEN
         ! otherwise, retrieve from database by Aggregation_LAI.F90
         CALL ncio_read_serial (fsrfdata, 'LAI_year', SITE_LAI_year)
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
         IF (DEF_LAI_MONTHLY) THEN
            CALL ncio_read_serial (fsrfdata, 'LAI_pfts_monthly', SITE_LAI_pfts_monthly)
            CALL ncio_read_serial (fsrfdata, 'SAI_pfts_monthly', SITE_SAI_pfts_monthly)
         ENDIF
#else
         IF (DEF_LAI_MONTHLY) THEN
            CALL ncio_read_serial (fsrfdata, 'LAI_monthly', SITE_LAI_monthly)
            CALL ncio_read_serial (fsrfdata, 'SAI_monthly', SITE_SAI_monthly)
         ELSE
            CALL ncio_read_serial (fsrfdata, 'LAI_8day', SITE_LAI_8day)
         ENDIF
#endif
      ENDIF

      IF ((.not. mksrfdata) .or. USE_SITE_lakedepth) THEN
         ! otherwise, retrieve from database by Aggregation_LakeDepth.F90
         CALL ncio_read_serial (fsrfdata, 'lakedepth', SITE_lakedepth)
      ENDIF

      IF ((.not. mksrfdata) .or. USE_SITE_soilreflectance) THEN
         ! otherwise, retrieve from database by Aggregation_SoilBrightness.F90
         CALL ncio_read_serial (fsrfdata, 'soil_s_v_alb', SITE_soil_s_v_alb)
         CALL ncio_read_serial (fsrfdata, 'soil_d_v_alb', SITE_soil_d_v_alb)
         CALL ncio_read_serial (fsrfdata, 'soil_s_n_alb', SITE_soil_s_n_alb)
         CALL ncio_read_serial (fsrfdata, 'soil_d_n_alb', SITE_soil_d_n_alb)
      ENDIF

      IF ((.not. mksrfdata) .or. USE_SITE_soilparameters) THEN
         ! otherwise, retrieve from database by Aggregation_SoilParameters.F90
         CALL ncio_read_serial (fsrfdata, 'soil_vf_quartz_mineral', SITE_soil_vf_quartz_mineral)
         CALL ncio_read_serial (fsrfdata, 'soil_vf_gravels       ', SITE_soil_vf_gravels       )
         CALL ncio_read_serial (fsrfdata, 'soil_vf_sand          ', SITE_soil_vf_sand          )
         CALL ncio_read_serial (fsrfdata, 'soil_vf_om            ', SITE_soil_vf_om            )
         CALL ncio_read_serial (fsrfdata, 'soil_wf_gravels       ', SITE_soil_wf_gravels       )
         CALL ncio_read_serial (fsrfdata, 'soil_wf_sand          ', SITE_soil_wf_sand          )
         CALL ncio_read_serial (fsrfdata, 'soil_OM_density       ', SITE_soil_OM_density       )
         CALL ncio_read_serial (fsrfdata, 'soil_BD_all           ', SITE_soil_BD_all           )
         CALL ncio_read_serial (fsrfdata, 'soil_theta_s          ', SITE_soil_theta_s          )
         CALL ncio_read_serial (fsrfdata, 'soil_k_s              ', SITE_soil_k_s              )
         CALL ncio_read_serial (fsrfdata, 'soil_csol             ', SITE_soil_csol             )
         CALL ncio_read_serial (fsrfdata, 'soil_tksatu           ', SITE_soil_tksatu           )
         CALL ncio_read_serial (fsrfdata, 'soil_tksatf           ', SITE_soil_tksatf           )
         CALL ncio_read_serial (fsrfdata, 'soil_tkdry            ', SITE_soil_tkdry            )
         CALL ncio_read_serial (fsrfdata, 'soil_k_solids         ', SITE_soil_k_solids         )
         CALL ncio_read_serial (fsrfdata, 'soil_psi_s            ', SITE_soil_psi_s            )
         CALL ncio_read_serial (fsrfdata, 'soil_lambda           ', SITE_soil_lambda           )
#ifdef vanGenuchten_Mualem_SOIL_MODEL
         CALL ncio_read_serial (fsrfdata, 'soil_theta_r          ', SITE_soil_theta_r          )
         CALL ncio_read_serial (fsrfdata, 'soil_alpha_vgm        ', SITE_soil_alpha_vgm        )
         CALL ncio_read_serial (fsrfdata, 'soil_L_vgm            ', SITE_soil_L_vgm            )
         CALL ncio_read_serial (fsrfdata, 'soil_n_vgm            ', SITE_soil_n_vgm            )
#else
         !SITE_soil_theta_r(:) = 0.
#endif
         CALL ncio_read_serial (fsrfdata, 'soil_BA_alpha         ', SITE_soil_BA_alpha         )
         CALL ncio_read_serial (fsrfdata, 'soil_BA_beta          ', SITE_soil_BA_beta          )
         CALL ncio_read_serial (fsrfdata, 'soil_texture          ', SITE_soil_texture          )
      ENDIF

      IF (DEF_USE_BEDROCK) THEN
         IF ((.not. mksrfdata) .or. USE_SITE_dbedrock) THEN
            ! otherwise, retrieve from database by Aggregation_DBedrock.F90
            CALL ncio_read_serial (fsrfdata, 'depth_to_bedrock', SITE_dbedrock)
         ENDIF
      ENDIF

      IF ((.not. mksrfdata) .or. USE_SITE_topography) THEN
         ! otherwise, retrieve from database by Aggregation_Topography.F90
         CALL ncio_read_serial (fsrfdata, 'elevation', SITE_topography)
         CALL ncio_read_serial (fsrfdata, 'elvstd   ', SITE_topostd   )

         IF (DEF_USE_Forcing_Downscaling) THEN
            CALL ncio_read_serial (fsrfdata, 'SITE_svf' , SITE_svf            )
            CALL ncio_read_serial (fsrfdata, 'SITE_cur' , SITE_cur            )
            CALL ncio_read_serial (fsrfdata, 'SITE_slp_type' , SITE_slp_type  )
            CALL ncio_read_serial (fsrfdata, 'SITE_asp_type' , SITE_asp_type  )
            CALL ncio_read_serial (fsrfdata, 'SITE_area_type', SITE_area_type )
            CALL ncio_read_serial (fsrfdata, 'SITE_sf_lut'   , SITE_sf_lut    )
         ENDIF
      ENDIF

   END SUBROUTINE read_surface_data_single

   ! -----
   SUBROUTINE read_urban_surface_data_single (fsrfdata, mksrfdata, mkrun)

   USE MOD_TimeManager
   USE MOD_NetCDFSerial
   USE MOD_Namelist
   USE MOD_Utils
   USE MOD_Vars_Global, only: PI, URBAN
   IMPLICIT NONE

   character(len=*), intent(in) :: fsrfdata
   logical, intent(in) :: mksrfdata
   logical, intent(in), optional :: mkrun

   ! Local Variables
   real(r8) :: lat_in, lon_in

      IF (ncio_var_exist(fsrfdata, 'latitude')) THEN
         CALL ncio_read_serial (fsrfdata, 'latitude',  lat_in)
         IF (lat_in /= SITE_lat_location) THEN
            write(*,*) 'Warning: Latitude mismatch: ', &
               lat_in, ' in data file and ', SITE_lat_location, 'in namelist.'
         ENDIF
         SITE_lat_location = lat_in
      ENDIF

      IF (ncio_var_exist(fsrfdata, 'longitude')) THEN
         CALL ncio_read_serial (fsrfdata, 'longitude', lon_in)
         IF (lon_in /= SITE_lon_location) THEN
            write(*,*) 'Warning: Longitude mismatch: ', &
               lon_in, ' in data file and ', SITE_lon_location, 'in namelist.'
         ENDIF
         SITE_lon_location = lon_in
      ENDIF

      CALL normalize_longitude (SITE_lon_location)

      IF (trim(fsrfdata) /= 'null') THEN
         SITE_landtype = URBAN
      ELSEIF (SITE_landtype /= URBAN) THEN
         write(*,*) 'Error! Please set namelist SITE_landtype first!'
         CALL CoLM_stop()
      ENDIF

      DEF_domain%edges = floor(SITE_lat_location)
      DEF_domain%edgen = DEF_domain%edges + 1.0
      DEF_domain%edgew = floor(SITE_lon_location)
      DEF_domain%edgee = DEF_domain%edgew + 1.0

      IF (.not. isgreenwich) THEN
         LocalLongitude = SITE_lon_location
      ENDIF

      IF (.not. present(mkrun)) THEN
         IF ((.not. mksrfdata) .or. USE_SITE_urban_paras) THEN

            CALL ncio_read_serial (fsrfdata, 'tree_area_fraction'         , SITE_fveg_urb  )
            CALL ncio_read_serial (fsrfdata, 'tree_mean_height'           , SITE_htop_urb  )
            CALL ncio_read_serial (fsrfdata, 'water_area_fraction'        , SITE_flake_urb )
            CALL ncio_read_serial (fsrfdata, 'roof_area_fraction'         , SITE_froof     )
            CALL ncio_read_serial (fsrfdata, 'building_mean_height'       , SITE_hroof     )
            CALL ncio_read_serial (fsrfdata, 'impervious_area_fraction'   , SITE_fgimp     )
            CALL ncio_read_serial (fsrfdata, 'wall_to_plan_area_ratio'    , SITE_lambdaw   )
            CALL ncio_read_serial (fsrfdata, 'resident_population_density', SITE_popden    )
IF (DEF_USE_CANYON_HWR) THEN
            CALL ncio_read_serial (fsrfdata, 'canyon_height_width_ratio'  , SITE_hlr       )
ELSE
            CALL ncio_read_serial (fsrfdata, 'wall_to_plan_area_ratio'    , SITE_lambdaw   )
            SITE_hlr       = SITE_lambdaw/4/SITE_froof
ENDIF
            SITE_fgper     = 1 - (SITE_fgimp-SITE_froof)/(1-SITE_froof-SITE_flake_urb)
            SITE_froof     = SITE_froof/(1-SITE_flake_urb)
            SITE_fveg_urb  = SITE_fveg_urb  * 100
            SITE_flake_urb = SITE_flake_urb * 100
         ENDIF
      ELSE
         CALL ncio_read_serial (fsrfdata, 'LAI_year'      , SITE_LAI_year   )
         CALL ncio_read_serial (fsrfdata, 'TREE_LAI'      , SITE_LAI_monthly)
         CALL ncio_read_serial (fsrfdata, 'TREE_SAI'      , SITE_SAI_monthly)

         CALL ncio_read_serial (fsrfdata, 'URBAN_TYPE'    , SITE_urbtyp     )
         CALL ncio_read_serial (fsrfdata, 'LUCY_id'       , SITE_lucyid     )
         CALL ncio_read_serial (fsrfdata, 'PCT_Tree'      , SITE_fveg_urb   )
         CALL ncio_read_serial (fsrfdata, 'URBAN_TREE_TOP', SITE_htop_urb   )
         CALL ncio_read_serial (fsrfdata, 'PCT_Water'     , SITE_flake_urb  )
         CALL ncio_read_serial (fsrfdata, 'WT_ROOF'       , SITE_froof      )
         CALL ncio_read_serial (fsrfdata, 'HT_ROOF'       , SITE_hroof      )
         CALL ncio_read_serial (fsrfdata, 'WTROAD_PERV'   , SITE_fgper      )
         CALL ncio_read_serial (fsrfdata, 'BUILDING_HLR'  , SITE_hlr        )
         CALL ncio_read_serial (fsrfdata, 'POP_DEN'       , SITE_popden     )

         CALL ncio_read_serial (fsrfdata, 'EM_ROOF'       , SITE_em_roof    )
         CALL ncio_read_serial (fsrfdata, 'EM_WALL'       , SITE_em_wall    )
         CALL ncio_read_serial (fsrfdata, 'EM_IMPROAD'    , SITE_em_gimp    )
         CALL ncio_read_serial (fsrfdata, 'EM_PERROAD'    , SITE_em_gper    )
         CALL ncio_read_serial (fsrfdata, 'T_BUILDING_MAX', SITE_t_roommax  )
         CALL ncio_read_serial (fsrfdata, 'T_BUILDING_MIN', SITE_t_roommin  )
         CALL ncio_read_serial (fsrfdata, 'THICK_ROOF'    , SITE_thickroof  )
         CALL ncio_read_serial (fsrfdata, 'THICK_WALL'    , SITE_thickwall  )

         CALL ncio_read_serial (fsrfdata, 'ALB_ROOF'      , SITE_alb_roof   )
         CALL ncio_read_serial (fsrfdata, 'ALB_WALL'      , SITE_alb_wall   )
         CALL ncio_read_serial (fsrfdata, 'ALB_IMPROAD'   , SITE_alb_gimp   )
         CALL ncio_read_serial (fsrfdata, 'ALB_PERROAD'   , SITE_alb_gper   )

         CALL ncio_read_serial (fsrfdata, 'CV_ROOF'       , SITE_cv_roof    )
         CALL ncio_read_serial (fsrfdata, 'CV_WALL'       , SITE_cv_wall    )
         CALL ncio_read_serial (fsrfdata, 'CV_IMPROAD'    , SITE_cv_gimp    )
         CALL ncio_read_serial (fsrfdata, 'TK_ROOF'       , SITE_tk_roof    )
         CALL ncio_read_serial (fsrfdata, 'TK_WALL'       , SITE_tk_wall    )
         CALL ncio_read_serial (fsrfdata, 'TK_IMPROAD'    , SITE_tk_gimp    )
      ENDIF

      IF ((.not. mksrfdata) .or. USE_SITE_lakedepth) THEN
         ! otherwise, retrieve from database by Aggregation_LakeDepth.F90
         CALL ncio_read_serial (fsrfdata, 'lakedepth', SITE_lakedepth)
      ENDIF

      IF ((.not. mksrfdata) .or. USE_SITE_soilreflectance) THEN
         ! otherwise, retrieve from database by Aggregation_SoilBrightness.F90
         CALL ncio_read_serial (fsrfdata, 'soil_s_v_alb', SITE_soil_s_v_alb)
         CALL ncio_read_serial (fsrfdata, 'soil_d_v_alb', SITE_soil_d_v_alb)
         CALL ncio_read_serial (fsrfdata, 'soil_s_n_alb', SITE_soil_s_n_alb)
         CALL ncio_read_serial (fsrfdata, 'soil_d_n_alb', SITE_soil_d_n_alb)
      ENDIF

      IF ((.not. mksrfdata) .or. USE_SITE_soilparameters) THEN
         ! otherwise, retrieve from database by Aggregation_SoilParameters.F90
         CALL ncio_read_serial (fsrfdata, 'soil_vf_quartz_mineral', SITE_soil_vf_quartz_mineral)
         CALL ncio_read_serial (fsrfdata, 'soil_vf_gravels       ', SITE_soil_vf_gravels       )
         CALL ncio_read_serial (fsrfdata, 'soil_vf_sand          ', SITE_soil_vf_sand          )
         CALL ncio_read_serial (fsrfdata, 'soil_vf_om            ', SITE_soil_vf_om            )
         CALL ncio_read_serial (fsrfdata, 'soil_wf_gravels       ', SITE_soil_wf_gravels       )
         CALL ncio_read_serial (fsrfdata, 'soil_wf_sand          ', SITE_soil_wf_sand          )
         CALL ncio_read_serial (fsrfdata, 'soil_OM_density       ', SITE_soil_OM_density       )
         CALL ncio_read_serial (fsrfdata, 'soil_BD_all           ', SITE_soil_BD_all           )
         CALL ncio_read_serial (fsrfdata, 'soil_theta_s          ', SITE_soil_theta_s          )
         CALL ncio_read_serial (fsrfdata, 'soil_k_s              ', SITE_soil_k_s              )
         CALL ncio_read_serial (fsrfdata, 'soil_csol             ', SITE_soil_csol             )
         CALL ncio_read_serial (fsrfdata, 'soil_tksatu           ', SITE_soil_tksatu           )
         CALL ncio_read_serial (fsrfdata, 'soil_tksatf           ', SITE_soil_tksatf           )
         CALL ncio_read_serial (fsrfdata, 'soil_tkdry            ', SITE_soil_tkdry            )
         CALL ncio_read_serial (fsrfdata, 'soil_k_solids         ', SITE_soil_k_solids         )
         CALL ncio_read_serial (fsrfdata, 'soil_psi_s            ', SITE_soil_psi_s            )
         CALL ncio_read_serial (fsrfdata, 'soil_lambda           ', SITE_soil_lambda           )
#ifdef vanGenuchten_Mualem_SOIL_MODEL
         CALL ncio_read_serial (fsrfdata, 'soil_theta_r          ', SITE_soil_theta_r          )
         CALL ncio_read_serial (fsrfdata, 'soil_alpha_vgm        ', SITE_soil_alpha_vgm        )
         CALL ncio_read_serial (fsrfdata, 'soil_L_vgm            ', SITE_soil_L_vgm            )
         CALL ncio_read_serial (fsrfdata, 'soil_n_vgm            ', SITE_soil_n_vgm            )
#else
         !SITE_soil_theta_r(:) = 0.
#endif
         CALL ncio_read_serial (fsrfdata, 'soil_BA_alpha         ', SITE_soil_BA_alpha         )
         CALL ncio_read_serial (fsrfdata, 'soil_BA_beta          ', SITE_soil_BA_beta          )
         CALL ncio_read_serial (fsrfdata, 'soil_texture          ', SITE_soil_texture          )
      ENDIF

      IF (DEF_USE_BEDROCK) THEN
         IF ((.not. mksrfdata) .or. USE_SITE_dbedrock) THEN
            ! otherwise, retrieve from database by Aggregation_DBedrock.F90
            CALL ncio_read_serial (fsrfdata, 'depth_to_bedrock', SITE_dbedrock)
         ENDIF
      ENDIF

      IF ((.not. mksrfdata) .or. USE_SITE_topography) THEN
         ! otherwise, retrieve from database by Aggregation_Topography.F90
         CALL ncio_read_serial (fsrfdata, 'elevation', SITE_topography)
         CALL ncio_read_serial (fsrfdata, 'elvstd   ', SITE_topostd   )

         IF (DEF_USE_Forcing_Downscaling) THEN
            CALL ncio_read_serial (fsrfdata, 'SITE_svf', SITE_svf             )
            CALL ncio_read_serial (fsrfdata, 'SITE_cur', SITE_cur             )
            CALL ncio_read_serial (fsrfdata, 'SITE_slp_type' , SITE_slp_type  )
            CALL ncio_read_serial (fsrfdata, 'SITE_asp_type' , SITE_asp_type  )
            CALL ncio_read_serial (fsrfdata, 'SITE_area_type', SITE_area_type )
            CALL ncio_read_serial (fsrfdata, 'SITE_sf_lut'   , SITE_sf_lut    )
         ENDIF
      ENDIF

   END SUBROUTINE read_urban_surface_data_single

   ! -----
   SUBROUTINE write_surface_data_single (numpatch, numpft)

   USE MOD_NetCDFSerial
   USE MOD_Namelist
   USE MOD_Const_LC
   IMPLICIT NONE

   integer, intent(in) :: numpatch
   integer, intent(in), optional :: numpft

   ! Local Variables
   character(len=256) :: fsrfdata
   integer :: ipft, iyear, itime
   character(len=8)   :: source

      fsrfdata = trim(DEF_dir_landdata) // '/srfdata.nc'

      CALL ncio_create_file (fsrfdata)

      CALL ncio_define_dimension (fsrfdata, 'soil',  nl_soil )
      CALL ncio_define_dimension (fsrfdata, 'azi', num_azimuth)
      CALL ncio_define_dimension (fsrfdata, 'zen', num_zenith)
      CALL ncio_define_dimension (fsrfdata, 'slope_type', num_slope_type)
      CALL ncio_define_dimension (fsrfdata, 'patch', numpatch)
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      CALL ncio_define_dimension (fsrfdata, 'pft', numpft)
#endif

      CALL ncio_define_dimension (fsrfdata, 'LAI_year', size(SITE_LAI_year))
      IF (DEF_LAI_MONTHLY) THEN
         CALL ncio_define_dimension (fsrfdata, 'month', 12)
      ELSE
         CALL ncio_define_dimension (fsrfdata, 'J8day', 46)
      ENDIF

      CALL ncio_write_serial (fsrfdata, 'latitude',  SITE_lat_location)
      CALL ncio_write_serial (fsrfdata, 'longitude', SITE_lon_location)

#ifdef LULC_USGS
      CALL ncio_write_serial (fsrfdata, 'USGS_classification', SITE_landtype)
#else
      CALL ncio_write_serial (fsrfdata, 'IGBP_classification', SITE_landtype)
#endif

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      CALL ncio_write_serial (fsrfdata, 'pfttyp',  SITE_pfttyp,  'pft')
      CALL ncio_put_attr     (fsrfdata, 'pfttyp', 'source', datasource(USE_SITE_pctpfts))
#endif
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      CALL ncio_write_serial (fsrfdata, 'pctpfts', SITE_pctpfts, 'pft')
      CALL ncio_put_attr     (fsrfdata, 'pctpfts', 'source', datasource(USE_SITE_pctpfts))
#endif
#if (defined CROP)
      IF (SITE_landtype == CROPLAND) THEN
         CALL ncio_write_serial (fsrfdata, 'croptyp', SITE_croptyp, 'patch')
         CALL ncio_write_serial (fsrfdata, 'pctcrop', SITE_pctcrop, 'patch')
         CALL ncio_put_attr     (fsrfdata, 'croptyp', 'source', datasource(USE_SITE_pctcrop))
         CALL ncio_put_attr     (fsrfdata, 'pctcrop', 'source', datasource(USE_SITE_pctcrop))
      ENDIF
#endif

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      CALL ncio_write_serial (fsrfdata, 'canopy_height_pfts', SITE_htop_pfts, 'pft')
      CALL ncio_put_attr     (fsrfdata, 'canopy_height_pfts', 'source', datasource(USE_SITE_htop))
#else
      CALL ncio_write_serial (fsrfdata, 'canopy_height', SITE_htop)
      CALL ncio_put_attr     (fsrfdata, 'canopy_height', 'source', datasource(USE_SITE_htop))
#endif

      source = datasource(USE_SITE_LAI)
      CALL ncio_write_serial (fsrfdata, 'LAI_year', SITE_LAI_year, 'LAI_year')
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      IF (DEF_LAI_MONTHLY) THEN
         CALL ncio_write_serial (fsrfdata, 'LAI_pfts_monthly', SITE_LAI_pfts_monthly, 'pft', 'month', 'LAI_year')
         CALL ncio_write_serial (fsrfdata, 'SAI_pfts_monthly', SITE_SAI_pfts_monthly, 'pft', 'month', 'LAI_year')
         CALL ncio_put_attr     (fsrfdata, 'LAI_pfts_monthly', 'source', source)
         CALL ncio_put_attr     (fsrfdata, 'SAI_pfts_monthly', 'source', source)
      ENDIF
#else
      IF (DEF_LAI_MONTHLY) THEN
         CALL ncio_write_serial (fsrfdata, 'LAI_monthly', SITE_LAI_monthly, 'month', 'LAI_year')
         CALL ncio_write_serial (fsrfdata, 'SAI_monthly', SITE_SAI_monthly, 'month', 'LAI_year')
         CALL ncio_put_attr     (fsrfdata, 'LAI_monthly', 'source', source)
         CALL ncio_put_attr     (fsrfdata, 'SAI_monthly', 'source', source)
      ELSE
         CALL ncio_write_serial (fsrfdata, 'LAI_8day', SITE_LAI_8day, 'J8day', 'LAI_year')
         CALL ncio_put_attr     (fsrfdata, 'LAI_8day', 'source', source)
      ENDIF
#endif

      CALL ncio_write_serial (fsrfdata, 'lakedepth', SITE_lakedepth)
      CALL ncio_put_attr     (fsrfdata, 'lakedepth', 'source', datasource(USE_SITE_lakedepth))

      source = datasource(USE_SITE_soilreflectance)
      CALL ncio_write_serial (fsrfdata, 'soil_s_v_alb', SITE_soil_s_v_alb)
      CALL ncio_put_attr     (fsrfdata, 'soil_s_v_alb', 'source', source)
      CALL ncio_write_serial (fsrfdata, 'soil_d_v_alb', SITE_soil_d_v_alb)
      CALL ncio_put_attr     (fsrfdata, 'soil_d_v_alb', 'source', source)
      CALL ncio_write_serial (fsrfdata, 'soil_s_n_alb', SITE_soil_s_n_alb)
      CALL ncio_put_attr     (fsrfdata, 'soil_s_n_alb', 'source', source)
      CALL ncio_write_serial (fsrfdata, 'soil_d_n_alb', SITE_soil_d_n_alb)
      CALL ncio_put_attr     (fsrfdata, 'soil_d_n_alb', 'source', source)

      source = datasource(USE_SITE_soilparameters)
      CALL ncio_write_serial (fsrfdata, 'soil_vf_quartz_mineral', SITE_soil_vf_quartz_mineral, 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_vf_gravels       ', SITE_soil_vf_gravels       , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_vf_sand          ', SITE_soil_vf_sand          , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_vf_om            ', SITE_soil_vf_om            , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_wf_gravels       ', SITE_soil_wf_gravels       , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_wf_sand          ', SITE_soil_wf_sand          , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_OM_density       ', SITE_soil_OM_density       , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_BD_all           ', SITE_soil_BD_all           , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_theta_s          ', SITE_soil_theta_s          , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_k_s              ', SITE_soil_k_s              , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_csol             ', SITE_soil_csol             , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_tksatu           ', SITE_soil_tksatu           , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_tksatf           ', SITE_soil_tksatf           , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_tkdry            ', SITE_soil_tkdry            , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_k_solids         ', SITE_soil_k_solids         , 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_vf_quartz_mineral', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_vf_gravels       ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_vf_sand          ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_vf_om            ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_wf_gravels       ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_wf_sand          ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_OM_density       ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_BD_all           ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_theta_s          ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_k_s              ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_csol             ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_tksatu           ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_tksatf           ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_tkdry            ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_k_solids         ', 'source', source)
      CALL ncio_write_serial (fsrfdata, 'soil_psi_s ', SITE_soil_psi_s , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_lambda', SITE_soil_lambda, 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_psi_s ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_lambda', 'source', source)
#ifdef vanGenuchten_Mualem_SOIL_MODEL
      CALL ncio_write_serial (fsrfdata, 'soil_theta_r  ', SITE_soil_theta_r  , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_alpha_vgm', SITE_soil_alpha_vgm, 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_L_vgm    ', SITE_soil_L_vgm    , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_n_vgm    ', SITE_soil_n_vgm    , 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_theta_r  ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_alpha_vgm', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_L_vgm    ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_n_vgm    ', 'source', source)
#endif
      CALL ncio_write_serial (fsrfdata, 'soil_BA_alpha', SITE_soil_BA_alpha, 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_BA_beta ', SITE_soil_BA_beta , 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_BA_alpha', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_BA_beta ', 'source', source)

      CALL ncio_write_serial (fsrfdata, 'soil_texture ', SITE_soil_texture)
      CALL ncio_put_attr     (fsrfdata, 'soil_texture ', 'source', source)

      IF(DEF_USE_BEDROCK)THEN
         CALL ncio_write_serial (fsrfdata, 'depth_to_bedrock', SITE_dbedrock)
         CALL ncio_put_attr     (fsrfdata, 'depth_to_bedrock', 'source', datasource(USE_SITE_dbedrock))
      ENDIF

      CALL ncio_write_serial (fsrfdata, 'elevation', SITE_topography)
      CALL ncio_put_attr     (fsrfdata, 'elevation', 'source', datasource(USE_SITE_topography))

      CALL ncio_write_serial (fsrfdata, 'elvstd', SITE_topostd)
      CALL ncio_put_attr     (fsrfdata, 'elvstd', 'source', datasource(USE_SITE_topography))

      ! used for downscaling
      IF (DEF_USE_Forcing_Downscaling) THEN
         CALL ncio_write_serial (fsrfdata, 'SITE_svf', SITE_svf)
         CALL ncio_write_serial (fsrfdata, 'SITE_cur', SITE_cur)
         CALL ncio_write_serial (fsrfdata, 'SITE_sf_lut'   , SITE_sf_lut, 'azi', 'zen')
         CALL ncio_write_serial (fsrfdata, 'SITE_slp_type' , SITE_slp_type , 'type')
         CALL ncio_write_serial (fsrfdata, 'SITE_asp_type' , SITE_asp_type , 'type')
         CALL ncio_write_serial (fsrfdata, 'SITE_area_type', SITE_area_type, 'type')
      ENDIF

   END SUBROUTINE write_surface_data_single

   ! -----
   SUBROUTINE write_urban_surface_data_single (numurban)

   USE MOD_NetCDFSerial
   USE MOD_Namelist
   USE MOD_Const_LC
   IMPLICIT NONE

   integer, intent(in) :: numurban

   ! Local Variables
   character(len=256) :: fsrfdata
   integer :: iyear, itime
   character(len=8)   :: source

      fsrfdata = trim(DEF_dir_landdata) // '/srfdata.nc'

      CALL ncio_create_file (fsrfdata)

      CALL ncio_define_dimension (fsrfdata, 'soil',  nl_soil )

      CALL ncio_define_dimension (fsrfdata, 'azi', num_azimuth)
      CALL ncio_define_dimension (fsrfdata, 'zen', num_zenith)
      CALL ncio_define_dimension (fsrfdata, 'slope_type', num_slope_type)

      CALL ncio_define_dimension (fsrfdata, 'patch', numurban)

      CALL ncio_define_dimension (fsrfdata, 'LAI_year', size(SITE_LAI_year))
      CALL ncio_define_dimension (fsrfdata, 'month', 12)

      CALL ncio_define_dimension (fsrfdata, 'ulev'    , 10)
      CALL ncio_define_dimension (fsrfdata, 'numsolar', 2 )
      CALL ncio_define_dimension (fsrfdata, 'numrad'  , 2 )

      CALL ncio_write_serial (fsrfdata, 'latitude',  SITE_lat_location)
      CALL ncio_write_serial (fsrfdata, 'longitude', SITE_lon_location)

      source = datasource(USE_SITE_urban_LAI)
      CALL ncio_write_serial (fsrfdata, 'LAI_year', SITE_LAI_year, 'LAI_year')
      CALL ncio_write_serial (fsrfdata, 'TREE_LAI', SITE_LAI_monthly, 'month', 'LAI_year')
      CALL ncio_write_serial (fsrfdata, 'TREE_SAI', SITE_SAI_monthly, 'month', 'LAI_year')
      CALL ncio_put_attr     (fsrfdata, 'TREE_LAI', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'TREE_SAI', 'source', source)

      CALL ncio_write_serial (fsrfdata, 'lakedepth', SITE_lakedepth)
      CALL ncio_put_attr     (fsrfdata, 'lakedepth', 'source', datasource(USE_SITE_lakedepth))

      CALL ncio_write_serial (fsrfdata, 'URBAN_TYPE'    , SITE_urbtyp     , 'patch')
      CALL ncio_write_serial (fsrfdata, 'LUCY_id'       , SITE_lucyid     , 'patch')
      !CALL ncio_put_attr    (fsrfdata, 'LUCY_id'       , 'source', source)
      source = datasource(USE_SITE_urban_paras)
      CALL ncio_write_serial (fsrfdata, 'PCT_Tree'      , SITE_fveg_urb   , 'patch')
      CALL ncio_write_serial (fsrfdata, 'URBAN_TREE_TOP', SITE_htop_urb   , 'patch')
      CALL ncio_write_serial (fsrfdata, 'PCT_Water'     , SITE_flake_urb  , 'patch')
      CALL ncio_write_serial (fsrfdata, 'WT_ROOF'       , SITE_froof      , 'patch')
      CALL ncio_write_serial (fsrfdata, 'HT_ROOF'       , SITE_hroof      , 'patch')
      CALL ncio_write_serial (fsrfdata, 'WTROAD_PERV'   , SITE_fgper      , 'patch')
      CALL ncio_write_serial (fsrfdata, 'BUILDING_HLR'  , SITE_hlr        , 'patch')
      CALL ncio_write_serial (fsrfdata, 'POP_DEN'       , SITE_popden     , 'patch')

      CALL ncio_put_attr     (fsrfdata, 'PCT_Tree'      , 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'URBAN_TREE_TOP', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'PCT_Water'     , 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'WT_ROOF'       , 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'HT_ROOF'       , 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'WTROAD_PERV'   , 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'BUILDING_HLR'  , 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'POP_DEN'       , 'source', source)

      source = datasource(USE_SITE_thermal_paras)
      CALL ncio_write_serial (fsrfdata, 'EM_ROOF'       , SITE_em_roof   , 'patch')
      CALL ncio_write_serial (fsrfdata, 'EM_WALL'       , SITE_em_wall   , 'patch')
      CALL ncio_write_serial (fsrfdata, 'EM_IMPROAD'    , SITE_em_gimp   , 'patch')
      CALL ncio_write_serial (fsrfdata, 'EM_PERROAD'    , SITE_em_gper   , 'patch')
      CALL ncio_write_serial (fsrfdata, 'T_BUILDING_MAX', SITE_t_roommax , 'patch')
      CALL ncio_write_serial (fsrfdata, 'T_BUILDING_MIN', SITE_t_roommin , 'patch')
      CALL ncio_write_serial (fsrfdata, 'THICK_ROOF'    , SITE_thickroof , 'patch')
      CALL ncio_write_serial (fsrfdata, 'THICK_WALL'    , SITE_thickwall , 'patch')

      CALL ncio_put_attr     (fsrfdata, 'EM_ROOF'       , 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'EM_WALL'       , 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'EM_IMPROAD'    , 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'EM_PERROAD'    , 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'T_BUILDING_MAX', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'T_BUILDING_MIN', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'THICK_ROOF'    , 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'THICK_WALL'    , 'source', source)

      CALL ncio_write_serial (fsrfdata, 'ALB_ROOF'      , SITE_alb_roof   , 'numrad', 'numsolar')
      CALL ncio_write_serial (fsrfdata, 'ALB_WALL'      , SITE_alb_wall   , 'numrad', 'numsolar')
      CALL ncio_write_serial (fsrfdata, 'ALB_IMPROAD'   , SITE_alb_gimp   , 'numrad', 'numsolar')
      CALL ncio_write_serial (fsrfdata, 'ALB_PERROAD'   , SITE_alb_gper   , 'numrad', 'numsolar')

      CALL ncio_put_attr     (fsrfdata, 'ALB_ROOF'      , 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'ALB_WALL'      , 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'ALB_IMPROAD'   , 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'ALB_PERROAD'   , 'source', source)

      CALL ncio_write_serial (fsrfdata, 'CV_ROOF'       , SITE_cv_roof    , 'ulev')
      CALL ncio_write_serial (fsrfdata, 'CV_WALL'       , SITE_cv_wall    , 'ulev')
      CALL ncio_write_serial (fsrfdata, 'CV_IMPROAD'    , SITE_cv_gimp    , 'ulev')
      CALL ncio_write_serial (fsrfdata, 'TK_ROOF'       , SITE_tk_roof    , 'ulev')
      CALL ncio_write_serial (fsrfdata, 'TK_WALL'       , SITE_tk_wall    , 'ulev')
      CALL ncio_write_serial (fsrfdata, 'TK_IMPROAD'    , SITE_tk_gimp    , 'ulev')
      CALL ncio_put_attr     (fsrfdata, 'CV_ROOF'       , 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'CV_WALL'       , 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'CV_IMPROAD'    , 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'TK_ROOF'       , 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'TK_WALL'       , 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'TK_IMPROAD'    , 'source', source)


      source = datasource(USE_SITE_soilreflectance)
      CALL ncio_write_serial (fsrfdata, 'soil_s_v_alb', SITE_soil_s_v_alb)
      CALL ncio_put_attr     (fsrfdata, 'soil_s_v_alb', 'source', source)
      CALL ncio_write_serial (fsrfdata, 'soil_d_v_alb', SITE_soil_d_v_alb)
      CALL ncio_put_attr     (fsrfdata, 'soil_d_v_alb', 'source', source)
      CALL ncio_write_serial (fsrfdata, 'soil_s_n_alb', SITE_soil_s_n_alb)
      CALL ncio_put_attr     (fsrfdata, 'soil_s_n_alb', 'source', source)
      CALL ncio_write_serial (fsrfdata, 'soil_d_n_alb', SITE_soil_d_n_alb)
      CALL ncio_put_attr     (fsrfdata, 'soil_d_n_alb', 'source', source)

      source = datasource(USE_SITE_soilparameters)
      CALL ncio_write_serial (fsrfdata, 'soil_vf_quartz_mineral', SITE_soil_vf_quartz_mineral, 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_vf_gravels       ', SITE_soil_vf_gravels       , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_vf_sand          ', SITE_soil_vf_sand          , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_vf_om            ', SITE_soil_vf_om            , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_wf_gravels       ', SITE_soil_wf_gravels       , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_wf_sand          ', SITE_soil_wf_sand          , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_OM_density       ', SITE_soil_OM_density       , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_BD_all           ', SITE_soil_BD_all           , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_theta_s          ', SITE_soil_theta_s          , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_k_s              ', SITE_soil_k_s              , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_csol             ', SITE_soil_csol             , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_tksatu           ', SITE_soil_tksatu           , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_tksatf           ', SITE_soil_tksatf           , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_tkdry            ', SITE_soil_tkdry            , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_k_solids         ', SITE_soil_k_solids         , 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_vf_quartz_mineral', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_vf_gravels       ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_vf_sand          ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_vf_om            ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_wf_gravels       ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_wf_sand          ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_OM_density       ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_BD_all           ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_theta_s          ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_k_s              ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_csol             ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_tksatu           ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_tksatf           ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_tkdry            ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_k_solids         ', 'source', source)
      CALL ncio_write_serial (fsrfdata, 'soil_psi_s ', SITE_soil_psi_s , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_lambda', SITE_soil_lambda, 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_psi_s ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_lambda', 'source', source)
#ifdef vanGenuchten_Mualem_SOIL_MODEL
      CALL ncio_write_serial (fsrfdata, 'soil_theta_r  ', SITE_soil_theta_r  , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_alpha_vgm', SITE_soil_alpha_vgm, 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_L_vgm    ', SITE_soil_L_vgm    , 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_n_vgm    ', SITE_soil_n_vgm    , 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_theta_r  ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_alpha_vgm', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_L_vgm    ', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_n_vgm    ', 'source', source)
#endif
      CALL ncio_write_serial (fsrfdata, 'soil_BA_alpha', SITE_soil_BA_alpha, 'soil')
      CALL ncio_write_serial (fsrfdata, 'soil_BA_beta ', SITE_soil_BA_beta , 'soil')
      CALL ncio_put_attr     (fsrfdata, 'soil_BA_alpha', 'source', source)
      CALL ncio_put_attr     (fsrfdata, 'soil_BA_beta ', 'source', source)

      CALL ncio_write_serial (fsrfdata, 'soil_texture ', SITE_soil_texture)
      CALL ncio_put_attr     (fsrfdata, 'soil_texture ', 'source', source)

      IF(DEF_USE_BEDROCK)THEN
         CALL ncio_write_serial (fsrfdata, 'depth_to_bedrock', SITE_dbedrock)
         CALL ncio_put_attr     (fsrfdata, 'depth_to_bedrock', 'source', datasource(USE_SITE_dbedrock))
      ENDIF

      CALL ncio_write_serial (fsrfdata, 'elevation', SITE_topography)
      CALL ncio_put_attr     (fsrfdata, 'elevation', 'source', datasource(USE_SITE_topography))

      CALL ncio_write_serial (fsrfdata, 'elvstd', SITE_topostd)
      CALL ncio_put_attr     (fsrfdata, 'elvstd', 'source', datasource(USE_SITE_topography))

      IF (DEF_USE_Forcing_Downscaling) THEN
         ! used for downscaling
         CALL ncio_write_serial (fsrfdata, 'SITE_svf', SITE_svf)
         CALL ncio_write_serial (fsrfdata, 'SITE_cur', SITE_cur)
         CALL ncio_write_serial (fsrfdata, 'SITE_sf_lut'   , SITE_sf_lut, 'azi', 'zen')
         CALL ncio_write_serial (fsrfdata, 'SITE_slp_type' , SITE_slp_type , 'type')
         CALL ncio_write_serial (fsrfdata, 'SITE_asp_type' , SITE_asp_type , 'type')
         CALL ncio_write_serial (fsrfdata, 'SITE_area_type', SITE_area_type, 'type')
      ENDIF

      ! used for downscaling
      IF (DEF_USE_Forcing_Downscaling) THEN
         CALL ncio_write_serial (fsrfdata, 'SITE_svf', SITE_svf)
         CALL ncio_write_serial (fsrfdata, 'SITE_cur', SITE_cur)
         CALL ncio_write_serial (fsrfdata, 'SITE_sf_lut', SITE_sf_lut, 'azi', 'zen')
         CALL ncio_write_serial (fsrfdata, 'SITE_slp_type', SITE_slp_type, 'type')
         CALL ncio_write_serial (fsrfdata, 'SITE_asp_type', SITE_asp_type, 'type')
         CALL ncio_write_serial (fsrfdata, 'SITE_area_type', SITE_area_type, 'type')
      ENDIF

   END SUBROUTINE write_urban_surface_data_single

   ! ---------
   character(len=8) FUNCTION datasource (is_site)

   IMPLICIT NONE
   logical, intent(in) :: is_site

      IF (is_site) THEN
         datasource = 'SITE'
      ELSE
         datasource = 'DATABASE'
      ENDIF

   END FUNCTION datasource

   ! ------
   SUBROUTINE single_srfdata_final ()

   IMPLICIT NONE

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      IF (allocated(SITE_pfttyp )) deallocate(SITE_pfttyp )
      IF (allocated(SITE_pctpfts)) deallocate(SITE_pctpfts)
#endif

#ifdef CROP
      IF (allocated(SITE_croptyp)) deallocate(SITE_croptyp)
      IF (allocated(SITE_pctcrop)) deallocate(SITE_pctcrop)
#endif

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      IF (allocated(SITE_htop_pfts)) deallocate(SITE_htop_pfts)
#endif

      IF (allocated(SITE_LAI_monthly)) deallocate(SITE_LAI_monthly)
      IF (allocated(SITE_SAI_monthly)) deallocate(SITE_SAI_monthly)
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      IF (allocated(SITE_LAI_pfts_monthly)) deallocate(SITE_LAI_pfts_monthly)
      IF (allocated(SITE_SAI_pfts_monthly)) deallocate(SITE_SAI_pfts_monthly)
#endif

      IF (allocated(SITE_LAI_year)) deallocate(SITE_LAI_year)
      IF (allocated(SITE_LAI_8day)) deallocate(SITE_LAI_8day)

      IF (allocated(SITE_soil_vf_quartz_mineral)) deallocate(SITE_soil_vf_quartz_mineral)
      IF (allocated(SITE_soil_vf_gravels       )) deallocate(SITE_soil_vf_gravels       )
      IF (allocated(SITE_soil_vf_sand          )) deallocate(SITE_soil_vf_sand          )
      IF (allocated(SITE_soil_vf_om            )) deallocate(SITE_soil_vf_om            )
      IF (allocated(SITE_soil_wf_gravels       )) deallocate(SITE_soil_wf_gravels       )
      IF (allocated(SITE_soil_wf_sand          )) deallocate(SITE_soil_wf_sand          )
      IF (allocated(SITE_soil_OM_density       )) deallocate(SITE_soil_OM_density       )
      IF (allocated(SITE_soil_BD_all           )) deallocate(SITE_soil_BD_all           )
      IF (allocated(SITE_soil_theta_s          )) deallocate(SITE_soil_theta_s          )
      IF (allocated(SITE_soil_k_s              )) deallocate(SITE_soil_k_s              )
      IF (allocated(SITE_soil_csol             )) deallocate(SITE_soil_csol             )
      IF (allocated(SITE_soil_tksatu           )) deallocate(SITE_soil_tksatu           )
      IF (allocated(SITE_soil_tksatf           )) deallocate(SITE_soil_tksatf           )
      IF (allocated(SITE_soil_tkdry            )) deallocate(SITE_soil_tkdry            )
      IF (allocated(SITE_soil_k_solids         )) deallocate(SITE_soil_k_solids         )
      IF (allocated(SITE_soil_psi_s            )) deallocate(SITE_soil_psi_s            )
      IF (allocated(SITE_soil_lambda           )) deallocate(SITE_soil_lambda           )
      IF (allocated(SITE_soil_theta_r          )) deallocate(SITE_soil_theta_r          )
#ifdef vanGenuchten_Mualem_SOIL_MODEL
      IF (allocated(SITE_soil_alpha_vgm        )) deallocate(SITE_soil_alpha_vgm        )
      IF (allocated(SITE_soil_L_vgm            )) deallocate(SITE_soil_L_vgm            )
      IF (allocated(SITE_soil_n_vgm            )) deallocate(SITE_soil_n_vgm            )
#endif
      IF (allocated(SITE_soil_BA_alpha         )) deallocate(SITE_soil_BA_alpha         )
      IF (allocated(SITE_soil_BA_beta          )) deallocate(SITE_soil_BA_beta          )

      IF (allocated(SITE_sf_lut                )) deallocate(SITE_sf_lut                )
      IF (allocated(SITE_slp_type              )) deallocate(SITE_slp_type              )
      IF (allocated(SITE_asp_type              )) deallocate(SITE_asp_type              )
      IF (allocated(SITE_area_type             )) deallocate(SITE_area_type             )

   END SUBROUTINE single_srfdata_final

END MODULE MOD_SingleSrfdata
#endif
