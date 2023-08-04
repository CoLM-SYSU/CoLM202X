#include <define.h>

MODULE MOD_Namelist

   !-----------------------------------------------------------------------
   ! DESCRIPTION:
   !
   !    Variables in namelist files and subrroutines to read namelist files.
   !
   ! Initial Authors: Shupeng Zhang, Zhongwang Wei, Xingjie Lu, Nan Wei,
   !                  Hua Yuan, Wenzong Dong et al., May 2023
   !-----------------------------------------------------------------------

   USE MOD_Precision, only: r8
   IMPLICIT NONE
   SAVE

   CHARACTER(len=256) :: DEF_CASE_NAME = 'CASENAME'

   ! ----- domain TYPE -----
   TYPE nl_domain_type
      REAL(r8) :: edges = -90.0
      REAL(r8) :: edgen = 90.0
      REAL(r8) :: edgew = -180.0
      REAL(r8) :: edgee = 180.0
   END TYPE nl_domain_type

   TYPE (nl_domain_type) :: DEF_domain

   INTEGER :: DEF_nx_blocks = 72
   INTEGER :: DEF_ny_blocks = 36
   INTEGER :: DEF_PIO_groupsize = 12

   ! ----- For Single Point -----
#ifdef SinglePoint

   CHARACTER(len=256) :: SITE_fsrfdata  = 'null'

   LOGICAL  :: USE_SITE_pctpfts         = .true.
   LOGICAL  :: USE_SITE_pctcrop         = .true.
   LOGICAL  :: USE_SITE_htop            = .true.
   LOGICAL  :: USE_SITE_LAI             = .true.
   LOGICAL  :: USE_SITE_lakedepth       = .true.
   LOGICAL  :: USE_SITE_soilreflectance = .true.
   LOGICAL  :: USE_SITE_soilparameters  = .true.
   LOGICAL  :: USE_SITE_dbedrock        = .true.
   LOGICAL  :: USE_SITE_topography      = .true.
   logical  :: USE_SITE_HistWriteBack   = .true.
#endif

   ! ----- simulation time type -----
   TYPE nl_simulation_time_type
      LOGICAL  :: greenwich   = .TRUE.
      INTEGER  :: start_year  = 2000
      INTEGER  :: start_month = 1
      INTEGER  :: start_day   = 1
      INTEGER  :: start_sec   = 0
      INTEGER  :: end_year    = 2003
      INTEGER  :: end_month   = 1
      INTEGER  :: end_day     = 1
      INTEGER  :: end_sec     = 0
      INTEGER  :: spinup_year = 2000
      INTEGER  :: spinup_month= 1
      INTEGER  :: spinup_day  = 1
      INTEGER  :: spinup_sec  = 0
      INTEGER  :: spinup_repeat = 1
      REAL(r8) :: timestep    = 3600.
   END TYPE nl_simulation_time_type

   TYPE (nl_simulation_time_type) :: DEF_simulation_time

   ! ----- directories -----
   CHARACTER(len=256) :: DEF_dir_rawdata  = 'path/to/rawdata/'
   CHARACTER(len=256) :: DEF_dir_runtime  = 'path/to/runtime/'
   CHARACTER(len=256) :: DEF_dir_output   = 'path/to/output/data'
   CHARACTER(len=256) :: DEF_dir_forcing  = 'path/to/forcing/data'

   CHARACTER(len=256) :: DEF_dir_landdata = 'path/to/landdata'
   CHARACTER(len=256) :: DEF_dir_restart  = 'path/to/restart'
   CHARACTER(len=256) :: DEF_dir_history  = 'path/to/history'

   CHARACTER(len=256) :: DEF_file_mesh    = 'path/to/mesh/file'
   REAL(r8) :: DEF_GRIDBASED_lon_res = 0.5
   REAL(r8) :: DEF_GRIDBASED_lat_res = 0.5

   CHARACTER(len=256) :: DEF_CatchmentMesh_data = 'path/to/catchment/data'

   CHARACTER(len=256) :: DEF_file_mesh_filter = 'path/to/mesh/filter'

   ! ----- Use surface data from existing dataset -----
   CHARACTER(len=256) :: DEF_dir_existing_srfdata = 'path/to/landdata'
   ! case 1: from a larger region
   LOGICAL :: USE_srfdata_from_larger_region   = .false.
   ! case 2: from gridded data with dimensions [patch,lon,lat] or [pft,lon,lat]
   !         only available for USGS/IGBP/PFT CLASSIFICATION
   LOGICAL :: USE_srfdata_from_3D_gridded_data = .false.

   ! ----- Subgrid scheme -----
   logical :: DEF_USE_PFT  = .false.
   logical :: DEF_USE_PC   = .false.
   logical :: DEF_SOLO_PFT = .false.
   logical :: DEF_FAST_PC  = .false.
   CHARACTER(len=256) :: DEF_SUBGRID_SCHEME = 'LCT'

   ! ----- Leaf Area Index -----
   !add by zhongwang wei @ sysu 2021/12/23
   !To allow read satellite observed LAI
   ! 06/2023, note by hua yuan: change DEF_LAI_CLIM to DEF_LAI_MONTHLY
   logical :: DEF_LAI_MONTHLY = .true.
   INTEGER :: DEF_Interception_scheme = 1  !1:CoLM；2:CLM4.5; 3:CLM5; 4:Noah-MP; 5:MATSIRO; 6:VIC

   ! ------LAI change and Land cover year setting ----------
   ! 06/2023, add by wenzong dong and hua yuan: use for updating LAI with simulation year
   LOGICAL :: DEF_LAI_CHANGE_YEARLY = .true.
   ! 05/2023, add by Xingjie Lu: use for updating LAI with leaf carbon
   LOGICAL :: DEF_USE_LAIFEEDBACK = .true.

   ! 06/2023, add by hua yuan and wenzong dong
   ! ------ Land use and land cover (LULC) related -------

   ! USE a static year land cover type
   INTEGER :: DEF_LC_YEAR      = 2005

   ! Options for LULCC year-to-year transfer schemes
   ! 1: Same Type Assignment scheme (STA), state variables assignment for the same TYPE (LC, PFT or PC)
   ! 2: Energy and Mass Conservation scheme (EMC), DO energy and mass conservation calculation
   INTEGER :: DEF_LULCC_SCHEME = 1

   ! ------ Urban model related -------
   ! Options for urban type scheme
   ! 1: NCAR Urban Classification, 3 urban type with Tall Building, High Density and Medium Density
   ! 2: LCZ Classification, 10 urban type with LCZ 1-10
   INTEGER :: DEF_URBAN_type_scheme = 1
   LOGICAL :: DEF_URBAN_ONLY   = .false.
   logical :: DEF_URBAN_RUN    = .false.
   LOGICAL :: DEF_URBAN_BEM    = .true.
   LOGICAL :: DEF_URBAN_TREE   = .true.
   LOGICAL :: DEF_URBAN_WATER  = .true.
   LOGICAL :: DEF_URBAN_LUCY   = .true.

   ! ------ SOIL parameters and supercool water setting -------
   LOGICAL :: DEF_USE_SOILPAR_UPS_FIT = .true.     ! soil hydraulic parameters are upscaled from rawdata (1km resolution)
                                                   ! to model patches through FIT algorithm (Montzka et al., 2017).
   INTEGER :: DEF_THERMAL_CONDUCTIVITY_SCHEME = 4  ! Options for soil thermal conductivity schemes
                                                   ! 1: Farouki (1981)
                                                   ! 2: Johansen(1975)
                                                   ! 3: Cote and Konrad (2005)
                                                   ! 4: Balland and Arp (2005)
                                                   ! 5: Lu et al. (2007)
                                                   ! 6: Tarnawski and Leong (2012)
                                                   ! 7: De Vries (1963)
                                                   ! 8: Yan Hengnian, He Hailong et al.(2019)
   LOGICAL :: DEF_USE_SUPERCOOL_WATER = .true.     ! supercooled soil water scheme, Niu & Yang (2006)


   ! Options for soil reflectance setting schemes
   ! 1: Guessed soil color type according to land cover classes
   ! 2: Read a global soil color map from CLM
   INTEGER :: DEF_SOIL_REFL_SCHEME = 2

   ! ----- Model settings -----
   LOGICAL :: DEF_LANDONLY                    = .true.
   LOGICAL :: DEF_USE_DOMINANT_PATCHTYPE      = .false.
   LOGICAL :: DEF_USE_VARIABLY_SATURATED_FLOW = .true.
   LOGICAL :: DEF_USE_BEDROCK                 = .false.
   LOGICAL :: DEF_USE_OZONESTRESS             = .false.
   LOGICAL :: DEF_USE_OZONEDATA               = .false.

   ! .true. for running SNICAR model
   logical :: DEF_USE_SNICAR                  = .true.

   ! .true. read aerosol deposition data from file or .false. set in the code
   logical :: DEF_Aerosol_Readin              = .true.

   ! .true. Read aerosol deposition climatology data or .false. yearly changed
   logical :: DEF_Aerosol_Clim                = .false.

   CHARACTER(len=5)   :: DEF_precip_phase_discrimination_scheme = 'II'
   CHARACTER(len=256) :: DEF_SSP='585' ! Co2 path for CMIP6 future scenario.

   ! ----- Initialization -----
   LOGICAL            :: DEF_USE_SOIL_INIT  = .false.
   CHARACTER(len=256) :: DEF_file_soil_init = 'null'

   CHARACTER(len=256) :: DEF_file_snowoptics = 'null'
   CHARACTER(len=256) :: DEF_file_snowaging  = 'null'

   ! ----- history -----
   LOGICAL  :: DEF_HISTORY_IN_VECTOR = .false.

   LOGICAL  :: DEF_hist_grid_as_forcing   = .false.
   REAL(r8) :: DEF_hist_lon_res = 0.5
   REAL(r8) :: DEF_hist_lat_res = 0.5

   CHARACTER(len=256) :: DEF_WRST_FREQ    = 'none'  ! write restart file frequency: TIMESTEP/HOURLY/DAILY/MONTHLY/YEARLY
   CHARACTER(len=256) :: DEF_HIST_FREQ    = 'none'  ! write history file frequency: TIMESTEP/HOURLY/DAILY/MONTHLY/YEARLY
   CHARACTER(len=256) :: DEF_HIST_groupby = 'MONTH' ! history file in one file: DAY/MONTH/YEAR
   CHARACTER(len=256) :: DEF_HIST_mode    = 'one'
   INTEGER :: DEF_REST_COMPRESS_LEVEL = 1
   INTEGER :: DEF_HIST_COMPRESS_LEVEL = 1

   CHARACTER(len=256) :: DEF_hist_vars_namelist = 'null'
   LOGICAL :: DEF_hist_vars_turnon_all = .true.

   ! ----- forcing -----
   CHARACTER(len=256) :: DEF_forcing_namelist = 'null'

   LOGICAL          :: DEF_USE_Forcing_Downscaling = .false.
   CHARACTER(len=5) :: DEF_DS_precipitation_adjust_scheme = 'II'
   CHARACTER(len=5) :: DEF_DS_longwave_adjust_scheme      = 'II'

   TYPE nl_forcing_type

      CHARACTER(len=256) :: dataset            = 'CRUNCEP'
      LOGICAL            :: solarin_all_band   = .true.
      REAL(r8)           :: HEIGHT_V           = 100.0
      REAL(r8)           :: HEIGHT_T           = 50.
      REAL(r8)           :: HEIGHT_Q           = 50.

      LOGICAL            :: regional           = .false.
      REAL(r8)           :: regbnd(4)          = (/-90.0, 90.0, -180.0, 180.0/)
      LOGICAL            :: has_missing_value  = .false.

      INTEGER            :: NVAR               = 8              ! variable number of forcing data
      INTEGER            :: startyr            = 2000           ! start year of forcing data        <MARK #1>
      INTEGER            :: startmo            = 1              ! start month of forcing data
      INTEGER            :: endyr              = 2003           ! end year of forcing data
      INTEGER            :: endmo              = 12             ! end month of forcing data
      INTEGER            :: dtime(8)           = (/21600,21600,21600,21600,0,21600,21600,21600/)
      INTEGER            :: offset(8)          = (/10800,10800,10800,10800,0,10800,0,10800/)
      INTEGER            :: nlands             = 1              ! land grid number in 1d

      LOGICAL            :: leapyear           = .false.        ! leapyear calendar
      LOGICAL            :: data2d             = .true.         ! data in 2 dimension (lon, lat)
      LOGICAL            :: hightdim           = .false.        ! have "z" dimension
      LOGICAL            :: dim2d              = .true.         ! lat/lon value in 2 dimension (lon, lat)

      CHARACTER(len=256) :: latname            = 'LATIXY'       ! dimension name of latitude
      CHARACTER(len=256) :: lonname            = 'LONGXY'       ! dimension name of longitude

      CHARACTER(len=256) :: groupby            = 'month'        ! file grouped by year/month

      CHARACTER(len=256) :: fprefix(8)          = (/ &
         'TPHWL6Hrly/clmforc.cruncep.V4.c2011.0.5d.TPQWL.', &
         'TPHWL6Hrly/clmforc.cruncep.V4.c2011.0.5d.TPQWL.', &
         'TPHWL6Hrly/clmforc.cruncep.V4.c2011.0.5d.TPQWL.', &
         'Precip6Hrly/clmforc.cruncep.V4.c2011.0.5d.Prec.', &
         'NULL                                           ', &
         'TPHWL6Hrly/clmforc.cruncep.V4.c2011.0.5d.TPQWL.', &
         'Solar6Hrly/clmforc.cruncep.V4.c2011.0.5d.Solr. ', &
         'TPHWL6Hrly/clmforc.cruncep.V4.c2011.0.5d.TPQWL.' /)
      CHARACTER(len=256) :: vname(8)           = (/ &
         'TBOT    ','QBOT    ','PSRF    ','PRECTmms', &
         'NULL    ','WIND    ','FSDS    ','FLDS    ' /)
      CHARACTER(len=256) :: tintalgo(8)        = (/ &
         'linear ','linear ','linear ','nearest', &
         'NULL   ','linear ','coszen ','linear ' /)

      CHARACTER(len=256) :: CBL_fprefix        = 'TPHWL6Hrly/clmforc.cruncep.V4.c2011.0.5d.TPQWL.'
      CHARACTER(len=256) :: CBL_vname          = 'blh'
      CHARACTER(len=256) :: CBL_tintalgo       = 'linear'
      INTEGER            :: CBL_dtime          = 21600
      INTEGER            :: CBL_offset         = 10800
   END TYPE nl_forcing_type

   TYPE (nl_forcing_type) :: DEF_forcing

   !CBL height
   LOGICAL            :: DEF_USE_CBL_HEIGHT = .false.
   !Plant Hydraulics
   LOGICAL            :: DEF_USE_PLANTHYDRAULICS = .true.
   !Semi-Analytic-Spin-Up
   LOGICAL            :: DEF_USE_SASU = .false.
   !Punctuated nitrogen addition Spin up
   LOGICAL            :: DEF_USE_PN   = .false.
   !Fertilisation on crop
   LOGICAL            :: DEF_USE_FERT = .true.
   !Nitrification and denitrification switch
   LOGICAL            :: DEF_USE_NITRIF = .true.
   !Soy nitrogen fixation
   LOGICAL            :: DEF_USE_CNSOYFIXN = .true.
   !Fire module
   LOGICAL            :: DEF_USE_FIRE = .false.


   ! ----- history variables -----
   TYPE history_var_type

      LOGICAL :: xy_us        = .true.
      LOGICAL :: xy_vs        = .true.
      LOGICAL :: xy_t         = .true.
      LOGICAL :: xy_q         = .true.
      LOGICAL :: xy_prc       = .true.
      LOGICAL :: xy_prl       = .true.
      LOGICAL :: xy_pbot      = .true.
      LOGICAL :: xy_frl       = .true.
      LOGICAL :: xy_solarin   = .true.
      LOGICAL :: xy_rain      = .true.
      LOGICAL :: xy_snow      = .true.
      LOGICAL :: xy_ozone     = .true.

      LOGICAL :: xy_hpbl      = .true.

      LOGICAL :: taux         = .true.
      LOGICAL :: tauy         = .true.
      LOGICAL :: fsena        = .true.
      LOGICAL :: lfevpa       = .true.
      LOGICAL :: fevpa        = .true.
      LOGICAL :: fsenl        = .true.
      LOGICAL :: fevpl        = .true.
      LOGICAL :: etr          = .true.
      LOGICAL :: fseng        = .true.
      LOGICAL :: fevpg        = .true.
      LOGICAL :: fgrnd        = .true.
      LOGICAL :: sabvsun      = .true.
      LOGICAL :: sabvsha      = .true.
      LOGICAL :: sabg         = .true.
      LOGICAL :: olrg         = .true.
      LOGICAL :: rnet         = .true.
      LOGICAL :: xerr         = .true.
      LOGICAL :: zerr         = .true.
      LOGICAL :: rsur         = .true.
      LOGICAL :: rsub         = .true.
      LOGICAL :: rnof         = .true.
      LOGICAL :: qintr        = .true.
      LOGICAL :: qinfl        = .true.
      LOGICAL :: qdrip        = .true.
      LOGICAL :: wat          = .true.
      LOGICAL :: assim        = .true.
      LOGICAL :: respc        = .true.
      LOGICAL :: qcharge      = .true.
      LOGICAL :: t_grnd       = .true.
      LOGICAL :: tleaf        = .true.
      LOGICAL :: ldew         = .true.
      LOGICAL :: scv          = .true.
      LOGICAL :: snowdp       = .true.
      LOGICAL :: fsno         = .true.
      LOGICAL :: sigf         = .true.
      LOGICAL :: green        = .true.
      LOGICAL :: lai          = .true.
      LOGICAL :: laisun       = .true.
      LOGICAL :: laisha       = .true.
      LOGICAL :: sai          = .true.
      LOGICAL :: alb          = .true.
      LOGICAL :: emis         = .true.
      LOGICAL :: z0m          = .true.
      LOGICAL :: trad         = .true.
      LOGICAL :: tref         = .true.
      LOGICAL :: qref         = .true.
#ifdef URBAN_MODEL
      LOGICAL :: fsen_roof    = .true.
      LOGICAL :: fsen_wsun    = .true.
      LOGICAL :: fsen_wsha    = .true.
      LOGICAL :: fsen_gimp    = .true.
      LOGICAL :: fsen_gper    = .true.
      LOGICAL :: fsen_urbl    = .true.
      LOGICAL :: lfevp_roof   = .true.
      LOGICAL :: lfevp_gimp   = .true.
      LOGICAL :: lfevp_gper   = .true.
      LOGICAL :: lfevp_urbl   = .true.
      LOGICAL :: fhac         = .true.
      LOGICAL :: fwst         = .true.
      LOGICAL :: fach         = .true.
      LOGICAL :: fhah         = .true.
      LOGICAL :: meta         = .true.
      LOGICAL :: vehc         = .true.
      LOGICAL :: t_room       = .true.
      LOGICAL :: tafu         = .true.
      LOGICAL :: t_roof       = .true.
      LOGICAL :: t_wall       = .true.
#endif
      LOGICAL :: assimsun      = .true. !1
      LOGICAL :: assimsha      = .true. !1
      LOGICAL :: etrsun        = .true. !1
      LOGICAL :: etrsha        = .true. !1
#ifdef BGC
      LOGICAL :: leafc              = .true.
      LOGICAL :: leafc_storage      = .true.
      LOGICAL :: leafc_xfer         = .true.
      LOGICAL :: frootc             = .true.
      LOGICAL :: frootc_storage     = .true.
      LOGICAL :: frootc_xfer        = .true.
      LOGICAL :: livestemc          = .true.
      LOGICAL :: livestemc_storage  = .true.
      LOGICAL :: livestemc_xfer     = .true.
      LOGICAL :: deadstemc          = .true.
      LOGICAL :: deadstemc_storage  = .true.
      LOGICAL :: deadstemc_xfer     = .true.
      LOGICAL :: livecrootc         = .true.
      LOGICAL :: livecrootc_storage = .true.
      LOGICAL :: livecrootc_xfer    = .true.
      LOGICAL :: deadcrootc         = .true.
      LOGICAL :: deadcrootc_storage = .true.
      LOGICAL :: deadcrootc_xfer    = .true.
      LOGICAL :: grainc             = .true.
      LOGICAL :: grainc_storage     = .true.
      LOGICAL :: grainc_xfer        = .true.
      LOGICAL :: leafn              = .true.
      LOGICAL :: leafn_storage      = .true.
      LOGICAL :: leafn_xfer         = .true.
      LOGICAL :: frootn             = .true.
      LOGICAL :: frootn_storage     = .true.
      LOGICAL :: frootn_xfer        = .true.
      LOGICAL :: livestemn          = .true.
      LOGICAL :: livestemn_storage  = .true.
      LOGICAL :: livestemn_xfer     = .true.
      LOGICAL :: deadstemn          = .true.
      LOGICAL :: deadstemn_storage  = .true.
      LOGICAL :: deadstemn_xfer     = .true.
      LOGICAL :: livecrootn         = .true.
      LOGICAL :: livecrootn_storage = .true.
      LOGICAL :: livecrootn_xfer    = .true.
      LOGICAL :: deadcrootn         = .true.
      LOGICAL :: deadcrootn_storage = .true.
      LOGICAL :: deadcrootn_xfer    = .true.
      LOGICAL :: grainn             = .true.
      LOGICAL :: grainn_storage     = .true.
      LOGICAL :: grainn_xfer        = .true.
      LOGICAL :: retrasn            = .true.
      LOGICAL :: gpp                = .true.
      LOGICAL :: downreg            = .true.
      LOGICAL :: ar                 = .true.
      LOGICAL :: cwdprod            = .true.
      LOGICAL :: cwddecomp          = .true.
      LOGICAL :: hr                 = .true.
      LOGICAL :: fpg                = .true.
      LOGICAL :: fpi                = .true.
      LOGICAL :: gpp_enftemp        = .false. !1
      LOGICAL :: gpp_enfboreal      = .false. !2
      LOGICAL :: gpp_dnfboreal      = .false. !3
      LOGICAL :: gpp_ebftrop        = .false. !4
      LOGICAL :: gpp_ebftemp        = .false. !5
      LOGICAL :: gpp_dbftrop        = .false. !6
      LOGICAL :: gpp_dbftemp        = .false. !7
      LOGICAL :: gpp_dbfboreal      = .false. !8
      LOGICAL :: gpp_ebstemp        = .false. !9
      LOGICAL :: gpp_dbstemp        = .false. !10
      LOGICAL :: gpp_dbsboreal      = .false. !11
      LOGICAL :: gpp_c3arcgrass     = .false. !12
      LOGICAL :: gpp_c3grass        = .false. !13
      LOGICAL :: gpp_c4grass        = .false. !14
      LOGICAL :: leafc_enftemp      = .false. !1
      LOGICAL :: leafc_enfboreal    = .false. !2
      LOGICAL :: leafc_dnfboreal    = .false. !3
      LOGICAL :: leafc_ebftrop      = .false. !4
      LOGICAL :: leafc_ebftemp      = .false. !5
      LOGICAL :: leafc_dbftrop      = .false. !6
      LOGICAL :: leafc_dbftemp      = .false. !7
      LOGICAL :: leafc_dbfboreal    = .false. !8
      LOGICAL :: leafc_ebstemp      = .false. !9
      LOGICAL :: leafc_dbstemp      = .false. !10
      LOGICAL :: leafc_dbsboreal    = .false. !11
      LOGICAL :: leafc_c3arcgrass   = .false. !12
      LOGICAL :: leafc_c3grass      = .false. !13
      LOGICAL :: leafc_c4grass      = .false. !14
#ifdef CROP
      LOGICAL :: cphase             = .true.
      LOGICAL :: gddmaturity        = .true.
      LOGICAL :: gddplant           = .true.
      LOGICAL :: vf                 = .true.
      LOGICAL :: hui                = .true.
      LOGICAL :: cropprod1c         = .true.
      LOGICAL :: cropprod1c_loss    = .true.
      LOGICAL :: cropseedc_deficit  = .true.
      LOGICAL :: grainc_to_cropprodc= .true.
      LOGICAL :: plantdate_rainfed_temp_corn= .true.
      LOGICAL :: plantdate_irrigated_temp_corn= .true.
      LOGICAL :: plantdate_rainfed_spwheat= .true.
      LOGICAL :: plantdate_irrigated_spwheat= .true.
      LOGICAL :: plantdate_rainfed_wtwheat= .true.
      LOGICAL :: plantdate_irrigated_wtwheat= .true.
      LOGICAL :: plantdate_rainfed_temp_soybean= .true.
      LOGICAL :: plantdate_irrigated_temp_soybean= .true.
      LOGICAL :: plantdate_rainfed_cotton= .true.
      LOGICAL :: plantdate_irrigated_cotton= .true.
      LOGICAL :: plantdate_rainfed_rice= .true.
      LOGICAL :: plantdate_irrigated_rice= .true.
      LOGICAL :: plantdate_rainfed_sugarcane= .true.
      LOGICAL :: plantdate_irrigated_sugarcane= .true.
      LOGICAL :: plantdate_rainfed_trop_corn= .true.
      LOGICAL :: plantdate_irrigated_trop_corn= .true.
      LOGICAL :: plantdate_rainfed_trop_soybean= .true.
      LOGICAL :: plantdate_irrigated_trop_soybean= .true.
      LOGICAL :: plantdate_unmanagedcrop= .true.
      LOGICAL :: cropprodc_rainfed_temp_corn= .true.
      LOGICAL :: cropprodc_irrigated_temp_corn= .true.
      LOGICAL :: cropprodc_rainfed_spwheat= .true.
      LOGICAL :: cropprodc_irrigated_spwheat= .true.
      LOGICAL :: cropprodc_rainfed_wtwheat= .true.
      LOGICAL :: cropprodc_irrigated_wtwheat= .true.
      LOGICAL :: cropprodc_rainfed_temp_soybean= .true.
      LOGICAL :: cropprodc_irrigated_temp_soybean= .true.
      LOGICAL :: cropprodc_rainfed_cotton= .true.
      LOGICAL :: cropprodc_irrigated_cotton= .true.
      LOGICAL :: cropprodc_rainfed_rice= .true.
      LOGICAL :: cropprodc_irrigated_rice= .true.
      LOGICAL :: cropprodc_rainfed_sugarcane= .true.
      LOGICAL :: cropprodc_irrigated_sugarcane= .true.
      LOGICAL :: cropprodc_rainfed_trop_corn= .true.
      LOGICAL :: cropprodc_irrigated_trop_corn= .true.
      LOGICAL :: cropprodc_rainfed_trop_soybean= .true.
      LOGICAL :: cropprodc_irrigated_trop_soybean= .true.
      LOGICAL :: cropprodc_unmanagedcrop= .true.

      LOGICAL :: grainc_to_seed     = .true.
      LOGICAL :: fert_to_sminn      = .true.

      LOGICAL :: huiswheat          = .true.
      LOGICAL :: pdcorn             = .true.
      LOGICAL :: pdswheat           = .true.
      LOGICAL :: pdwwheat           = .true.
      LOGICAL :: pdsoybean          = .true.
      LOGICAL :: pdcotton           = .true.
      LOGICAL :: pdrice1            = .true.
      LOGICAL :: pdrice2            = .true.
      LOGICAL :: pdsugarcane        = .true.
      LOGICAL :: fertnitro_corn     = .true.
      LOGICAL :: fertnitro_swheat   = .true.
      LOGICAL :: fertnitro_wwheat   = .true.
      LOGICAL :: fertnitro_soybean  = .true.
      LOGICAL :: fertnitro_cotton   = .true.
      LOGICAL :: fertnitro_rice1    = .true.
      LOGICAL :: fertnitro_rice2    = .true.
      LOGICAL :: fertnitro_sugarcane= .true.
#endif
      LOGICAL :: ndep_to_sminn      = .true.
      LOGICAL :: CONC_O2_UNSAT      = .true.
      LOGICAL :: O2_DECOMP_DEPTH_UNSAT = .true.
      LOGICAL :: abm                = .true.
      LOGICAL :: gdp                = .true.
      LOGICAL :: peatf              = .true.
      LOGICAL :: hdm                = .true.
      LOGICAL :: lnfm               = .true.
#endif

      LOGICAL :: t_soisno     = .true.
      LOGICAL :: wliq_soisno  = .true.
      LOGICAL :: wice_soisno  = .true.

      LOGICAL :: h2osoi       = .true.
      LOGICAL :: rstfacsun    = .true.
      LOGICAL :: rstfacsha    = .true.
      LOGICAL :: gssun        = .true.
      LOGICAL :: gssha        = .true.
      LOGICAL :: rootr        = .true.
      LOGICAL :: vegwp        = .true.
      LOGICAL :: BD_all       = .true.
      LOGICAL :: wfc          = .true.
      LOGICAL :: OM_density   = .true.
      LOGICAL :: wdsrf        = .true.
      LOGICAL :: zwt          = .true.
      LOGICAL :: wa           = .true.

      LOGICAL :: t_lake       = .true.
      LOGICAL :: lake_icefrac = .true.

#ifdef BGC
      LOGICAL :: litr1c_vr    = .true.
      LOGICAL :: litr2c_vr    = .true.
      LOGICAL :: litr3c_vr    = .true.
      LOGICAL :: soil1c_vr    = .true.
      LOGICAL :: soil2c_vr    = .true.
      LOGICAL :: soil3c_vr    = .true.
      LOGICAL :: cwdc_vr      = .true.
      LOGICAL :: litr1n_vr    = .true.
      LOGICAL :: litr2n_vr    = .true.
      LOGICAL :: litr3n_vr    = .true.
      LOGICAL :: soil1n_vr    = .true.
      LOGICAL :: soil2n_vr    = .true.
      LOGICAL :: soil3n_vr    = .true.
      LOGICAL :: cwdn_vr      = .true.
      LOGICAL :: sminn_vr     = .true.
#endif

      LOGICAL :: ustar        = .true.
      LOGICAL :: ustar2       = .true.
      LOGICAL :: tstar        = .true.
      LOGICAL :: qstar        = .true.
      LOGICAL :: zol          = .true.
      LOGICAL :: rib          = .true.
      LOGICAL :: fm           = .true.
      LOGICAL :: fh           = .true.
      LOGICAL :: fq           = .true.
      LOGICAL :: us10m        = .true.
      LOGICAL :: vs10m        = .true.
      LOGICAL :: fm10m        = .true.
      LOGICAL :: sr           = .true.
      LOGICAL :: solvd        = .true.
      LOGICAL :: solvi        = .true.
      LOGICAL :: solnd        = .true.
      LOGICAL :: solni        = .true.
      LOGICAL :: srvd         = .true.
      LOGICAL :: srvi         = .true.
      LOGICAL :: srnd         = .true.
      LOGICAL :: srni         = .true.

      LOGICAL :: solvdln      = .true.
      LOGICAL :: solviln      = .true.
      LOGICAL :: solndln      = .true.
      LOGICAL :: solniln      = .true.
      LOGICAL :: srvdln       = .true.
      LOGICAL :: srviln       = .true.
      LOGICAL :: srndln       = .true.
      LOGICAL :: srniln       = .true.

      LOGICAL :: rsubs_bsn    = .true.
      LOGICAL :: rsubs_hru    = .true.
      LOGICAL :: riv_height   = .true.
      LOGICAL :: riv_veloct   = .true.
      LOGICAL :: wdsrf_hru    = .true.
      LOGICAL :: veloc_hru    = .true.

   END TYPE history_var_type

   TYPE (history_var_type) :: DEF_hist_vars

CONTAINS

   SUBROUTINE read_namelist (nlfile)

      USE MOD_SPMD_Task
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: nlfile

      ! Local variables
      LOGICAL :: fexists
      INTEGER :: ivar
      INTEGER :: ierr

      namelist /nl_colm/          &
         DEF_CASE_NAME,           &
         DEF_domain,              &
#ifdef SinglePoint
         SITE_fsrfdata,            &
         USE_SITE_pctpfts,         &
         USE_SITE_pctcrop,         &
         USE_SITE_htop,            &
         USE_SITE_LAI,             &
         USE_SITE_lakedepth,       &
         USE_SITE_soilreflectance, &
         USE_SITE_soilparameters,  &
         USE_SITE_dbedrock,        &
         USE_SITE_topography,      &
         USE_SITE_HistWriteBack,   &
#endif
         DEF_nx_blocks,                   &
         DEF_ny_blocks,                   &
         DEF_PIO_groupsize,               &
         DEF_simulation_time,             &
         DEF_dir_rawdata,                 &
         DEF_dir_runtime,                 &
         DEF_dir_output,                  &
         DEF_file_mesh,                   &
         DEF_GRIDBASED_lon_res,           &
         DEF_GRIDBASED_lat_res,           &
#ifdef CATCHMENT
         DEF_CatchmentMesh_data,          &
#endif
         DEF_file_mesh_filter,            &

         DEF_USE_PFT,                     &
         DEF_USE_PC,                      &
         DEF_FAST_PC,                     &
         DEF_SUBGRID_SCHEME,              &

         DEF_LAI_MONTHLY,                 &   !add by zhongwang wei @ sysu 2021/12/23
         DEF_Interception_scheme,         &   !add by zhongwang wei @ sysu 2022/05/23
         DEF_SSP,                         &   !add by zhongwang wei @ sysu 2023/02/07

         DEF_LAI_CHANGE_YEARLY,           &
         DEF_USE_LAIFEEDBACK,             &   !add by Xingjie Lu, use for updating LAI with leaf carbon

         DEF_LC_YEAR,                     &
         DEF_LULCC_SCHEME,                &

         DEF_URBAN_type_scheme,           &
         DEF_URBAN_ONLY,                  &
         DEF_URBAN_RUN,                   &   !add by hua yuan, open urban model or not
         DEF_URBAN_BEM,                   &   !add by hua yuan, open urban BEM model or not
         DEF_URBAN_TREE,                  &   !add by hua yuan, modeling urban tree or not
         DEF_URBAN_WATER,                 &   !add by hua yuan, modeling urban water or not
         DEF_URBAN_LUCY,                  &

         DEF_USE_SOILPAR_UPS_FIT,         &
         DEF_THERMAL_CONDUCTIVITY_SCHEME, &
         DEF_USE_SUPERCOOL_WATER,         &
         DEF_SOIL_REFL_SCHEME,            &

         DEF_dir_existing_srfdata,        &
         USE_srfdata_from_larger_region,  &
         USE_srfdata_from_3D_gridded_data,&

         DEF_USE_CBL_HEIGHT,              &   !add by zhongwang wei @ sysu 2022/12/31
         DEF_USE_PLANTHYDRAULICS,         &   !add by xingjie lu @ sysu 2023/05/28
         DEF_USE_SASU,                    &   !add by Xingjie Lu @ sysu 2023/06/27
         DEF_USE_PN,                      &   !add by Xingjie Lu @ sysu 2023/06/27
         DEF_USE_FERT,                    &   !add by Xingjie Lu @ sysu 2023/06/27
         DEF_USE_NITRIF,                  &   !add by Xingjie Lu @ sysu 2023/06/27
         DEF_USE_CNSOYFIXN,               &   !add by Xingjie Lu @ sysu 2023/06/27
         DEF_USE_FIRE,                    &   !add by Xingjie Lu @ sysu 2023/06/27

         DEF_LANDONLY,                    &
         DEF_USE_DOMINANT_PATCHTYPE,      &
         DEF_USE_VARIABLY_SATURATED_FLOW, &
         DEF_USE_BEDROCK,                 &
         DEF_USE_OZONESTRESS,             &
         DEF_USE_OZONEDATA,               &
         DEF_USE_SNICAR,                  &
         DEF_Aerosol_Readin,              &
         DEF_Aerosol_Clim,                &

         DEF_precip_phase_discrimination_scheme, &

         DEF_USE_SOIL_INIT,               &
         DEF_file_soil_init,              &

         DEF_file_snowoptics,             &
         DEF_file_snowaging ,             &

         DEF_forcing_namelist,            &

         DEF_USE_Forcing_Downscaling,        &
         DEF_DS_precipitation_adjust_scheme, &
         DEF_DS_longwave_adjust_scheme,      &

         DEF_HISTORY_IN_VECTOR,           &
         DEF_hist_lon_res,                &
         DEF_hist_lat_res,                &
         DEF_hist_grid_as_forcing,        &
         DEF_WRST_FREQ,                   &
         DEF_HIST_FREQ,                   &
         DEF_HIST_groupby,                &
         DEF_HIST_mode,                   &
         DEF_REST_COMPRESS_LEVEL,         &
         DEF_HIST_COMPRESS_LEVEL,         &
         DEF_hist_vars_namelist,          &
         DEF_hist_vars_turnon_all

      namelist /nl_colm_forcing/ DEF_dir_forcing, DEF_forcing
      namelist /nl_colm_history/ DEF_hist_vars

      ! ----- open the namelist file -----
      IF (p_is_master) THEN

         open(10, status='OLD', file=nlfile, form="FORMATTED")
         read(10, nml=nl_colm, iostat=ierr)
         IF (ierr /= 0) THEN
            write(*,*) ' ***** ERROR: Problem reading namelist.'
            write(*,*) trim(nlfile), ierr
#ifdef USEMPI
            CALL mpi_abort (p_comm_glb, p_err)
#endif
         ENDIF
         close(10)

         open(10, status='OLD', file=trim(DEF_forcing_namelist), form="FORMATTED")
         read(10, nml=nl_colm_forcing, iostat=ierr)
         IF (ierr /= 0) THEN
            write(*,*) ' ***** ERROR: Problem reading forcing namelist.'
            write(*,*) trim(DEF_forcing_namelist), ierr
         ENDIF
         close(10)
#ifdef SinglePoint
         DEF_forcing%has_missing_value = .false.
#endif

         DEF_dir_landdata = trim(DEF_dir_output) // '/' // trim(adjustl(DEF_CASE_NAME)) // '/landdata'
         DEF_dir_restart  = trim(DEF_dir_output) // '/' // trim(adjustl(DEF_CASE_NAME)) // '/restart'
         DEF_dir_history  = trim(DEF_dir_output) // '/' // trim(adjustl(DEF_CASE_NAME)) // '/history'

         CALL system('mkdir -p ' // trim(adjustl(DEF_dir_output  )))
         CALL system('mkdir -p ' // trim(adjustl(DEF_dir_landdata)))
         CALL system('mkdir -p ' // trim(adjustl(DEF_dir_restart )))
         CALL system('mkdir -p ' // trim(adjustl(DEF_dir_history )))

#ifdef SinglePoint
         DEF_nx_blocks = 360
         DEF_ny_blocks = 180
         DEF_HIST_mode = 'one'
#endif


! ===============================================================
! ----- Macros&Namelist conflicts and dependency management -----
! ===============================================================


! ----- SOIL model related ------ Macros&Namelist conflicts and dependency management
#if (defined vanGenuchten_Mualem_SOIL_MODEL)
         write(*,*) '                  *****                  '
         write(*,*) 'Note: DEF_USE_VARIABLY_SATURATED_FLOW is automaticlly set to .true.  '
         write(*,*) 'when using vanGenuchten_Mualem_SOIL_MODEL. '
         DEF_USE_VARIABLY_SATURATED_FLOW = .true.
#endif


! ----- subgrid type related ------ Macros&Namelist conflicts and dependency management

#ifdef LULC_IGBP_PFT
         DEF_USE_PFT = .true.
         DEF_USE_PC  = .false.
         DEF_FAST_PC = .false.
#endif

#ifdef LULC_IGBP_PC
         DEF_USE_PC   = .true.
         DEF_USE_PFT  = .false.
         DEF_SOLO_PFT = .false.
#endif

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
         IF (.not.DEF_LAI_MONTHLY) THEN
            write(*,*) '                  *****                  '
            write(*,*) 'Warning: 8-day LAI data is not supported for '
            write(*,*) 'LULC_IGBP_PFT and LULC_IGBP_PC.'
            write(*,*) 'Changed to monthly data, set DEF_LAI_MONTHLY = .true.'
            DEF_LAI_MONTHLY = .true.
         ENDIF
#endif


! ----- BGC and CROP model related ------ Macros&Namelist conflicts and dependency management

#ifndef BGC
         IF(DEF_USE_LAIFEEDBACK)then
            DEF_USE_LAIFEEDBACK = .false.
            write(*,*) '                  *****                  '
            write(*,*) 'Warning: LAI feedback is not supported for BGC off.'
            write(*,*) 'DEF_USE_LAIFEEDBACK is set to false automatically when BGC is turned off.'
         ENDIF

         IF(DEF_USE_SASU)then
            DEF_USE_SASU = .false.
            write(*,*) '                  *****                  '
            write(*,*) 'Warning: Semi-Analytic Spin-up is on when BGC is off.'
            write(*,*) 'DEF_USE_SASU is set to false automatically when BGC is turned off.'
         ENDIF

         IF(DEF_USE_PN)then
            DEF_USE_PN = .false.
            write(*,*) '                  *****                  '
            write(*,*) 'Warning: Punctuated nitrogen addition spin up is on when BGC is off.'
            write(*,*) 'DEF_USE_PN is set to false automatically when BGC is turned off.'
         ENDIF

         IF(DEF_USE_NITRIF)then
            DEF_USE_NITRIF = .false.
            write(*,*) '                  *****                  '
            write(*,*) 'Warning: Nitrification-Denitrification is on when BGC is off.'
            write(*,*) 'DEF_USE_NITRIF is set to false automatically when BGC is turned off.'
         ENDIF

         IF(DEF_USE_FIRE)then
            DEF_USE_FIRE = .false.
            write(*,*) '                  *****                  '
            write(*,*) 'Warning: Fire model is on when BGC is off.'
            write(*,*) 'DEF_USE_FIRE is set to false automatically when BGC is turned off.'
         ENDIF
#endif

#ifndef CROP
         IF(DEF_USE_FERT)then
            DEF_USE_FERT = .false.
            write(*,*) '                  *****                  '
            write(*,*) 'Warning: Fertilization is on when CROP is off.'
            write(*,*) 'DEF_USE_FERT is set to false automatically when CROP is turned off.'
         ENDIF

         IF(DEF_USE_CNSOYFIXN)then
            DEF_USE_CNSOYFIXN = .false.
            write(*,*) '                  *****                  '
            write(*,*) 'Warning: Soy nitrogen fixation is on when CROP is off.'
            write(*,*) 'DEF_USE_CNSOYFIXN is set to false automatically when CROP is turned off.'
         ENDIF
#endif

         IF(.not. DEF_USE_OZONESTRESS)then
            IF(DEF_USE_OZONEDATA)then
               DEF_USE_OZONEDATA = .false.
               write(*,*) '                  *****                  '
               write(*,*) 'Warning: DEF_USE_OZONEDATA is not supported for OZONESTRESS off.'
               write(*,*) 'DEF_USE_OZONEDATA is set to false automatically.'
            ENDIF
         ENDIF


! ----- SNICAR model ------ Macros&Namelist conflicts and dependency management

         DEF_file_snowoptics = trim(DEF_dir_runtime)//'/snicar/snicar_optics_5bnd_mam_c211006.nc'
         DEF_file_snowaging  = trim(DEF_dir_runtime)//'/snicar/snicar_drdt_bst_fit_60_c070416.nc'

         IF (.not. DEF_USE_SNICAR) THEN
            IF (DEF_Aerosol_Readin) THEN
               DEF_Aerosol_Readin = .false.
               write(*,*) '                  *****                  '
               write(*,*) 'Warning: DEF_Aerosol_Readin is not needed for DEF_USE_SNICAR off. '
               write(*,*) 'DEF_Aerosol_Readin is set to false automatically.'
            ENDIF
         ENDIF


! ----- Urban model ----- Macros&Namelist conflicts and dependency management

#ifdef URBAN_MODEL
         DEF_URBAN_RUN = .true.

         IF (DEF_USE_SNICAR) THEN
            write(*,*) '                  *****                  '
            write(*,*) 'Note: SNICAR is not applied for URBAN model, but for other land covers. '
         ENDIF
#else
         IF (DEF_URBAN_RUN) then
            write(*,*) '                  *****                  '
            write(*,*) 'Note: The Urban model is not opened. IF you want to run Urban model '
            write(*,*) 'please #define URBAN_MODEL in define.h. otherwise DEF_URBAN_RUN will '
            write(*,*) 'be set to false automatically.'
            DEF_URBAN_RUN = .false.
         ENDIF
#endif


! ----- LULCC ----- Macros&Namelist conflicts and dependency management

#ifdef LULCC

#if (defined LULC_USGS || defined BGC)
         write(*,*) '                  *****                  '
         write(*,*) 'Fatal ERROR: LULCC is not supported for LULC_USGS/BGC at present. STOP! '
         STOP
#endif
         IF (.not.DEF_LAI_MONTHLY) THEN
            write(*,*) '                  *****                  '
            write(*,*) 'Note: When LULCC is opened, DEF_LAI_MONTHLY '
            write(*,*) 'will be set to true automatically.'
            DEF_LAI_MONTHLY = .true.
         ENDIF

         IF (.not.DEF_LAI_CHANGE_YEARLY) THEN
            write(*,*) '                  *****                  '
            write(*,*) 'Note: When LULCC is opened, DEF_LAI_CHANGE_YEARLY '
            write(*,*) 'will be set to true automatically.'
            DEF_LAI_CHANGE_YEARLY = .true.
         ENDIF

#if (defined LULC_IGBP_PC || defined URBAN)
         write(*,*) '                  *****                  '
         write(*,*) 'Fatal ERROR: LULCC is not supported for LULC_IGBP_PC/URBAN at present. STOP! '
         write(*,*) 'It is coming soon. '
         STOP
#endif

#endif


! ----- [Complement IF needed] ----- Macros&Namelist conflicts and dependency management



! -----END Macros&Namelist conflicts and dependency management -----
! ===============================================================


      ENDIF


#ifdef USEMPI
      CALL mpi_bcast (DEF_CASE_NAME,    256, mpi_character, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_domain%edges,   1, mpi_real8,     p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_domain%edgen,   1, mpi_real8,     p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_domain%edgew,   1, mpi_real8,     p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_domain%edgee,   1, mpi_real8,     p_root, p_comm_glb, p_err)

      CALL mpi_bcast (DEF_nx_blocks,     1, mpi_integer, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_ny_blocks,     1, mpi_integer, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_PIO_groupsize, 1, mpi_integer, p_root, p_comm_glb, p_err)

      CALL mpi_bcast (DEF_simulation_time%greenwich,     1, mpi_logical, p_root, p_comm_glb, p_err)

      CALL mpi_bcast (DEF_simulation_time%start_year,    1, mpi_integer, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_simulation_time%start_month,   1, mpi_integer, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_simulation_time%start_day,     1, mpi_integer, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_simulation_time%start_sec,     1, mpi_integer, p_root, p_comm_glb, p_err)

      CALL mpi_bcast (DEF_simulation_time%end_year,      1, mpi_integer, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_simulation_time%end_month,     1, mpi_integer, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_simulation_time%end_day,       1, mpi_integer, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_simulation_time%end_sec,       1, mpi_integer, p_root, p_comm_glb, p_err)

      CALL mpi_bcast (DEF_simulation_time%spinup_year,   1, mpi_integer, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_simulation_time%spinup_month,  1, mpi_integer, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_simulation_time%spinup_day,    1, mpi_integer, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_simulation_time%spinup_sec,    1, mpi_integer, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_simulation_time%spinup_repeat, 1, mpi_integer, p_root, p_comm_glb, p_err)

      CALL mpi_bcast (DEF_simulation_time%timestep,     1, mpi_real8,   p_root, p_comm_glb, p_err)

      CALL mpi_bcast (DEF_dir_rawdata,  256, mpi_character, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_dir_runtime,  256, mpi_character, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_dir_output,   256, mpi_character, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_dir_forcing,  256, mpi_character, p_root, p_comm_glb, p_err)

      CALL mpi_bcast (DEF_dir_landdata, 256, mpi_character, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_dir_restart,  256, mpi_character, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_dir_history,  256, mpi_character, p_root, p_comm_glb, p_err)

#if (defined GRIDBASED || defined UNSTRUCTURED)
      CALL mpi_bcast (DEF_file_mesh,    256, mpi_character, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_GRIDBASED_lon_res, 1, mpi_real8, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_GRIDBASED_lat_res, 1, mpi_real8, p_root, p_comm_glb, p_err)
#endif

#ifdef CATCHMENT
      CALL mpi_bcast (DEF_CatchmentMesh_data, 256, mpi_character, p_root, p_comm_glb, p_err)
#endif

      CALL mpi_bcast (DEF_file_mesh_filter, 256, mpi_character, p_root, p_comm_glb, p_err)

      CALL mpi_bcast (DEF_dir_existing_srfdata, 256, mpi_character, p_root, p_comm_glb, p_err)
      call mpi_bcast (USE_srfdata_from_larger_region,   1, mpi_logical, p_root, p_comm_glb, p_err)
      call mpi_bcast (USE_srfdata_from_3D_gridded_data, 1, mpi_logical, p_root, p_comm_glb, p_err)

      ! 07/2023, added by yuan: subgrid setting related
      CALL mpi_bcast (DEF_USE_PFT,          1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_USE_PC,           1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_FAST_PC,          1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_SUBGRID_SCHEME, 256, mpi_character, p_root, p_comm_glb, p_err)

      CALL mpi_bcast (DEF_LAI_CHANGE_YEARLY,  1, mpi_logical, p_root, p_comm_glb, p_err)

      ! 05/2023, added by Xingjie lu
      CALL mpi_bcast (DEF_USE_LAIFEEDBACK,    1, mpi_logical, p_root, p_comm_glb, p_err)

      ! LULC related
      CALL mpi_bcast (DEF_LC_YEAR,         1, mpi_integer, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_LULCC_SCHEME,    1, mpi_integer, p_root, p_comm_glb, p_err)

      CALL mpi_bcast (DEF_URBAN_type_scheme, 1, mpi_integer, p_root, p_comm_glb, p_err)
      ! 05/2023, added by yuan
      CALL mpi_bcast (DEF_URBAN_ONLY,      1, mpi_logical, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_URBAN_RUN,       1, mpi_logical, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_URBAN_BEM,       1, mpi_logical, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_URBAN_TREE,      1, mpi_logical, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_URBAN_WATER,     1, mpi_logical, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_URBAN_LUCY,      1, mpi_logical, p_root, p_comm_glb, p_err)

      ! 06/2023, added by weinan
      CALL mpi_bcast (DEF_USE_SOILPAR_UPS_FIT,          1, mpi_logical, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_THERMAL_CONDUCTIVITY_SCHEME,  1, mpi_integer, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_USE_SUPERCOOL_WATER,          1, mpi_logical, p_root, p_comm_glb, p_err)

      ! 06/2023, added by hua yuan
      CALL mpi_bcast (DEF_SOIL_REFL_SCHEME,             1, mpi_integer, p_root, p_comm_glb, p_err)

      call mpi_bcast (DEF_LAI_MONTHLY,         1, mpi_logical, p_root, p_comm_glb, p_err)
      call mpi_bcast (DEF_Interception_scheme, 1, mpi_integer, p_root, p_comm_glb, p_err)
      call mpi_bcast (DEF_SSP, 256, mpi_character, p_root, p_comm_glb, p_err)

      call mpi_bcast (DEF_USE_CBL_HEIGHT, 1, mpi_logical, p_root, p_comm_glb, p_err)
      call mpi_bcast (DEF_USE_PLANTHYDRAULICS, 1, mpi_logical, p_root, p_comm_glb, p_err)
      call mpi_bcast (DEF_USE_SASU, 1, mpi_logical, p_root, p_comm_glb, p_err)
      call mpi_bcast (DEF_USE_PN, 1, mpi_logical, p_root, p_comm_glb, p_err)
      call mpi_bcast (DEF_USE_FERT, 1, mpi_logical, p_root, p_comm_glb, p_err)
      call mpi_bcast (DEF_USE_NITRIF, 1, mpi_logical, p_root, p_comm_glb, p_err)
      call mpi_bcast (DEF_USE_CNSOYFIXN, 1, mpi_logical, p_root, p_comm_glb, p_err)
      call mpi_bcast (DEF_USE_FIRE, 1, mpi_logical, p_root, p_comm_glb, p_err)

      call mpi_bcast (DEF_LANDONLY,                   1, mpi_logical, p_root, p_comm_glb, p_err)
      call mpi_bcast (DEF_USE_DOMINANT_PATCHTYPE,     1, mpi_logical, p_root, p_comm_glb, p_err)
      call mpi_bcast (DEF_USE_VARIABLY_SATURATED_FLOW,1, mpi_logical, p_root, p_comm_glb, p_err)
      call mpi_bcast (DEF_USE_BEDROCK                ,1, mpi_logical, p_root, p_comm_glb, p_err)
      call mpi_bcast (DEF_USE_OZONESTRESS            ,1, mpi_logical, p_root, p_comm_glb, p_err)
      call mpi_bcast (DEF_USE_OZONEDATA              ,1, mpi_logical, p_root, p_comm_glb, p_err)

      CALL mpi_bcast (DEF_precip_phase_discrimination_scheme, 5, mpi_character, p_root, p_comm_glb, p_err)

      call mpi_bcast (DEF_USE_SOIL_INIT,    1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_file_soil_init, 256, mpi_character, p_root, p_comm_glb, p_err)

      call mpi_bcast (DEF_USE_SNICAR,        1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_file_snowoptics, 256, mpi_character, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_file_snowaging , 256, mpi_character, p_root, p_comm_glb, p_err)

      call mpi_bcast (DEF_Aerosol_Readin,    1, mpi_logical,   p_root, p_comm_glb, p_err)
      call mpi_bcast (DEF_Aerosol_Clim,      1, mpi_logical,   p_root, p_comm_glb, p_err)

      CALL mpi_bcast (DEF_HISTORY_IN_VECTOR, 1, mpi_logical,  p_root, p_comm_glb, p_err)

      CALL mpi_bcast (DEF_hist_lon_res,  1, mpi_real8, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_lat_res,  1, mpi_real8, p_root, p_comm_glb, p_err)

      CALL mpi_bcast (DEF_hist_grid_as_forcing, 1, mpi_logical, p_root, p_comm_glb, p_err)

      CALL mpi_bcast (DEF_WRST_FREQ,         256, mpi_character, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_HIST_FREQ,         256, mpi_character, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_HIST_groupby,      256, mpi_character, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_HIST_mode,         256, mpi_character, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_REST_COMPRESS_LEVEL, 1, mpi_integer,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_HIST_COMPRESS_LEVEL, 1, mpi_integer,   p_root, p_comm_glb, p_err)

      CALL mpi_bcast (DEF_USE_Forcing_Downscaling,        1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_DS_precipitation_adjust_scheme, 5, mpi_character, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_DS_longwave_adjust_scheme,      5, mpi_character, p_root, p_comm_glb, p_err)

      CALL mpi_bcast (DEF_forcing%dataset,          256, mpi_character, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_forcing%solarin_all_band,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_forcing%HEIGHT_V,           1, mpi_real8,     p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_forcing%HEIGHT_T,           1, mpi_real8,     p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_forcing%HEIGHT_Q,           1, mpi_real8,     p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_forcing%regional,           1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_forcing%regbnd,             4, mpi_real8,     p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_forcing%has_missing_value,  1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_forcing%NVAR,               1, mpi_integer,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_forcing%startyr,            1, mpi_integer,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_forcing%startmo,            1, mpi_integer,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_forcing%endyr,              1, mpi_integer,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_forcing%endmo,              1, mpi_integer,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_forcing%dtime,              8, mpi_integer,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_forcing%offset,             8, mpi_integer,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_forcing%nlands,             1, mpi_integer,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_forcing%leapyear,           1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_forcing%data2d,             1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_forcing%hightdim,           1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_forcing%dim2d,              1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_forcing%latname,          256, mpi_character, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_forcing%lonname,          256, mpi_character, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_forcing%groupby,          256, mpi_character, p_root, p_comm_glb, p_err)
      DO ivar = 1, 8
         CALL mpi_bcast (DEF_forcing%fprefix(ivar),  256, mpi_character, p_root, p_comm_glb, p_err)
         CALL mpi_bcast (DEF_forcing%vname(ivar),    256, mpi_character, p_root, p_comm_glb, p_err)
         CALL mpi_bcast (DEF_forcing%tintalgo(ivar), 256, mpi_character, p_root, p_comm_glb, p_err)
      ENDDO
      call mpi_bcast (DEF_forcing%CBL_fprefix,      256, mpi_character, p_root, p_comm_glb, p_err)
      call mpi_bcast (DEF_forcing%CBL_vname,        256, mpi_character, p_root, p_comm_glb, p_err)
      call mpi_bcast (DEF_forcing%CBL_tintalgo,     256, mpi_character, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_forcing%CBL_dtime,          1, mpi_integer,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_forcing%CBL_offset,         1, mpi_integer,   p_root, p_comm_glb, p_err)

      CALL mpi_bcast (DEF_file_snowoptics,  256, mpi_character, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_file_snowaging,   256, mpi_character, p_root, p_comm_glb, p_err)
#endif

      CALL sync_hist_vars (set_defaults = .true.)

      IF (p_is_master) THEN

         inquire (file=trim(DEF_hist_vars_namelist), exist=fexists)
         IF (.not. fexists) THEN
            write(*,*) 'History namelist file: ', trim(DEF_hist_vars_namelist), ' does not exist.'
         ELSE
            open(10, status='OLD', file=trim(DEF_hist_vars_namelist), form="FORMATTED")
            read(10, nml=nl_colm_history, iostat=ierr)
            close(10)
         ENDIF

      ENDIF

      CALL sync_hist_vars (set_defaults = .false.)

   END SUBROUTINE read_namelist

   ! ---------------
   SUBROUTINE sync_hist_vars (set_defaults)

      IMPLICIT NONE

      LOGICAL, intent(in) :: set_defaults

      CALL sync_hist_vars_one (DEF_hist_vars%xy_us       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%xy_vs       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%xy_t        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%xy_q        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%xy_prc      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%xy_prl      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%xy_pbot     ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%xy_frl      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%xy_solarin  ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%xy_rain     ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%xy_snow     ,  set_defaults)

      CALL sync_hist_vars_one (DEF_hist_vars%xy_hpbl     ,  set_defaults)

      CALL sync_hist_vars_one (DEF_hist_vars%taux        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%tauy        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%fsena       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%lfevpa      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%fevpa       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%fsenl       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%fevpl       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%etr         ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%fseng       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%fevpg       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%fgrnd       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%sabvsun     ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%sabvsha     ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%sabg        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%olrg        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%rnet        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%xerr        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%zerr        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%rsur        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%rsub        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%rnof        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%qintr       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%qinfl       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%qdrip       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%wat         ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%assim       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%respc       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%qcharge     ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%t_grnd      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%tleaf       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%ldew        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%scv         ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%snowdp      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%fsno        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%sigf        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%green       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%lai         ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%laisun      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%laisha      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%sai         ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%alb         ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%emis        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%z0m         ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%trad        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%tref        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%qref        ,  set_defaults)
#ifdef URBAN_MODEL
      CALL sync_hist_vars_one (DEF_hist_vars%fsen_roof   ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%fsen_wsun   ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%fsen_wsha   ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%fsen_gimp   ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%fsen_gper   ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%fsen_urbl   ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%lfevp_roof  ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%lfevp_gimp  ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%lfevp_gper  ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%lfevp_urbl  ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%fhac        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%fwst        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%fach        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%fhah        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%meta        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%vehc        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%t_room      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%tafu        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%t_roof      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%t_wall      ,  set_defaults)
#endif
#ifdef BGC
      CALL sync_hist_vars_one (DEF_hist_vars%leafc              ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%leafc_storage      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%leafc_xfer         ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%frootc             ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%frootc_storage     ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%frootc_xfer        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%livestemc          ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%livestemc_storage  ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%livestemc_xfer     ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%deadstemc          ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%deadstemc_storage  ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%deadstemc_xfer     ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%livecrootc         ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%livecrootc_storage ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%livecrootc_xfer    ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%deadcrootc         ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%deadcrootc_storage ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%deadcrootc_xfer    ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%grainc             ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%grainc_storage     ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%grainc_xfer        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%leafn              ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%leafn_storage      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%leafn_xfer         ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%frootn             ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%frootn_storage     ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%frootn_xfer        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%livestemn          ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%livestemn_storage  ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%livestemn_xfer     ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%deadstemn          ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%deadstemn_storage  ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%deadstemn_xfer     ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%livecrootn         ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%livecrootn_storage ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%livecrootn_xfer    ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%deadcrootn         ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%deadcrootn_storage ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%deadcrootn_xfer    ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%grainn             ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%grainn_storage     ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%grainn_xfer        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%retrasn            ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%gpp                ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%downreg            ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%ar                 ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%cwdprod            ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%cwddecomp          ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%hr                 ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%fpg                ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%fpi                ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%gpp_enftemp        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%gpp_enfboreal      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%gpp_dnfboreal      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%gpp_ebftrop        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%gpp_ebftemp        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%gpp_dbftrop        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%gpp_dbftemp        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%gpp_dbfboreal      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%gpp_ebstemp        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%gpp_dbstemp        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%gpp_dbsboreal      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%gpp_c3arcgrass     ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%gpp_c3grass        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%gpp_c4grass        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%leafc_enftemp      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%leafc_enfboreal    ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%leafc_dnfboreal    ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%leafc_ebftrop      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%leafc_ebftemp      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%leafc_dbftrop      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%leafc_dbftemp      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%leafc_dbfboreal    ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%leafc_ebstemp      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%leafc_dbstemp      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%leafc_dbsboreal    ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%leafc_c3arcgrass   ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%leafc_c3grass      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%leafc_c4grass      ,  set_defaults)
         CALL sync_hist_vars_one (DEF_hist_vars%assimsun        ,  set_defaults)
         CALL sync_hist_vars_one (DEF_hist_vars%assimsha        ,  set_defaults)
         CALL sync_hist_vars_one (DEF_hist_vars%etrsun        ,  set_defaults)
         CALL sync_hist_vars_one (DEF_hist_vars%etrsha        ,  set_defaults)
#ifdef CROP
      CALL sync_hist_vars_one (DEF_hist_vars%cphase                          , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%cropprod1c                      , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%cropprod1c_loss                 , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%cropseedc_deficit               , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%grainc_to_cropprodc             , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%grainc_to_seed                  , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%hui                             , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%vf                              , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%gddmaturity                     , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%gddplant                        , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%plantdate_rainfed_temp_corn     , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%plantdate_irrigated_temp_corn   , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%plantdate_rainfed_spwheat       , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%plantdate_irrigated_spwheat     , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%plantdate_rainfed_wtwheat       , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%plantdate_irrigated_wtwheat     , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%plantdate_rainfed_temp_soybean  , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%plantdate_irrigated_temp_soybean, set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%plantdate_rainfed_cotton        , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%plantdate_irrigated_cotton      , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%plantdate_rainfed_rice          , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%plantdate_irrigated_rice        , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%plantdate_rainfed_sugarcane     , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%plantdate_irrigated_sugarcane   , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%plantdate_rainfed_trop_corn     , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%plantdate_irrigated_trop_corn   , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%plantdate_rainfed_trop_soybean  , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%plantdate_irrigated_trop_soybean, set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%plantdate_unmanagedcrop         , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%cropprodc_rainfed_temp_corn     , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%cropprodc_irrigated_temp_corn   , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%cropprodc_rainfed_spwheat       , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%cropprodc_irrigated_spwheat     , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%cropprodc_rainfed_wtwheat       , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%cropprodc_irrigated_wtwheat     , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%cropprodc_rainfed_temp_soybean  , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%cropprodc_irrigated_temp_soybean, set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%cropprodc_rainfed_cotton        , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%cropprodc_irrigated_cotton      , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%cropprodc_rainfed_rice          , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%cropprodc_irrigated_rice        , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%cropprodc_rainfed_sugarcane     , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%cropprodc_irrigated_sugarcane   , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%cropprodc_rainfed_trop_corn     , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%cropprodc_irrigated_trop_corn   , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%cropprodc_rainfed_trop_soybean  , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%cropprodc_irrigated_trop_soybean, set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%cropprodc_unmanagedcrop         , set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%fert_to_sminn                   , set_defaults)
#endif
      CALL sync_hist_vars_one (DEF_hist_vars%ndep_to_sminn                   , set_defaults)
      if(DEF_USE_FIRE)then
         CALL sync_hist_vars_one (DEF_hist_vars%abm                          , set_defaults)
         CALL sync_hist_vars_one (DEF_hist_vars%gdp                          , set_defaults)
         CALL sync_hist_vars_one (DEF_hist_vars%peatf                        , set_defaults)
         CALL sync_hist_vars_one (DEF_hist_vars%hdm                          , set_defaults)
         CALL sync_hist_vars_one (DEF_hist_vars%lnfm                         , set_defaults)
      endif
#endif

      CALL sync_hist_vars_one (DEF_hist_vars%t_soisno    ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%wliq_soisno ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%wice_soisno ,  set_defaults)

      CALL sync_hist_vars_one (DEF_hist_vars%h2osoi      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%rstfacsun   ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%rstfacsha   ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%gssun   ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%gssha   ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%rootr       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%vegwp       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%BD_all      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%wfc         ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%OM_density  ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%wdsrf       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%zwt         ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%wa          ,  set_defaults)

      CALL sync_hist_vars_one (DEF_hist_vars%t_lake      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%lake_icefrac,  set_defaults)

#ifdef BGC
      CALL sync_hist_vars_one (DEF_hist_vars%litr1c_vr   ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%litr2c_vr   ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%litr3c_vr   ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%soil1c_vr   ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%soil2c_vr   ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%soil3c_vr   ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%cwdc_vr     ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%litr1n_vr   ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%litr2n_vr   ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%litr3n_vr   ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%soil1n_vr   ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%soil2n_vr   ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%soil3n_vr   ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%cwdn_vr     ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%sminn_vr    ,  set_defaults)
#endif

      CALL sync_hist_vars_one (DEF_hist_vars%ustar       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%ustar2      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%tstar       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%qstar       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%zol         ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%rib         ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%fm          ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%fh          ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%fq          ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%us10m       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%vs10m       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%fm10m       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%sr          ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%solvd       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%solvi       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%solnd       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%solni       ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%srvd        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%srvi        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%srnd        ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%srni        ,  set_defaults)

      CALL sync_hist_vars_one (DEF_hist_vars%solvdln     ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%solviln     ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%solndln     ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%solniln     ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%srvdln      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%srviln      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%srndln      ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%srniln      ,  set_defaults)

      CALL sync_hist_vars_one (DEF_hist_vars%rsubs_bsn   ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%rsubs_hru   ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%riv_height  ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%riv_veloct  ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%wdsrf_hru   ,  set_defaults)
      CALL sync_hist_vars_one (DEF_hist_vars%veloc_hru   ,  set_defaults)

   END SUBROUTINE sync_hist_vars

   SUBROUTINE sync_hist_vars_one (onoff, set_defaults)

      USE MOD_SPMD_Task
      IMPLICIT NONE

      LOGICAL, intent(inout) :: onoff
      LOGICAL, intent(in)    :: set_defaults

      IF (p_is_master) THEN
         IF (set_defaults) THEN
            onoff = DEF_hist_vars_turnon_all
         ENDIF
      ENDIF

#ifdef USEMPI
      CALL mpi_bcast (onoff, 1, mpi_logical, p_root, p_comm_glb, p_err)
#endif

   END SUBROUTINE sync_hist_vars_one

END MODULE MOD_Namelist
