#include <define.h>

MODULE mod_namelist

   USE precision, only: r8
   IMPLICIT NONE
   SAVE

   CHARACTER(len=256) :: DEF_CASE_NAME = 'CLMCRU'         

   ! ----- domain TYPE -----
   TYPE nl_domain_type
      REAL(r8) :: edges = -90.0      
      REAL(r8) :: edgen = 90.0      
      REAL(r8) :: edgew = -180.0   
      REAL(r8) :: edgee = 180.0      
   END TYPE nl_domain_type

   TYPE (nl_domain_type) :: DEF_domain

   INTEGER :: DEF_nx_blocks = 1
   INTEGER :: DEF_ny_blocks = 1
   INTEGER :: DEF_PIO_groupsize

#ifdef SinglePoint
   REAL(r8) :: SITE_lon_location = 0.
   REAL(r8) :: SITE_lat_location = 0.
   
   INTEGER  :: SITE_landtype = 1
   CHARACTER(len=256) :: SITE_fsrfdata

#if (defined PFT_CLASSIFICATION || defined PC_CLASSIFICATION)
   REAL(r8) :: SITE_pct_pfts(0:N_PFT-1) = 1./N_PFT
#ifdef CROP
   REAL(r8) :: SITE_pctcrop (N_CFT) = 1./N_CFT
#endif
#endif
   LOGICAL  :: USE_SITE_htop            = .true.
   LOGICAL  :: USE_SITE_LAI             = .true.
   LOGICAL  :: USE_SITE_lakedepth       = .true.
   LOGICAL  :: USE_SITE_soilreflectance = .true.
   LOGICAL  :: USE_SITE_soilparameters  = .true.
#ifdef USE_DEPTH_TO_BEDROCK
   LOGICAL  :: USE_SITE_dbedrock = .true.
#endif
   LOGICAL  :: USE_SITE_Forcing  = .false.
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
      INTEGER  :: spinup_repeat = 2
      REAL(r8) :: timestep    = 3600.
   END TYPE nl_simulation_time_type

   TYPE (nl_simulation_time_type) :: DEF_simulation_time

   ! ----- directories -----
   CHARACTER(len=256) :: DEF_dir_rawdata  = 'path/to/rawdata/'
   CHARACTER(len=256) :: DEF_dir_output   = 'path/to/output/data' 
   CHARACTER(len=256) :: DEF_dir_forcing  = 'path/to/forcing/data' 

   CHARACTER(len=256) :: DEF_dir_landdata = 'path/to/landdata'
   CHARACTER(len=256) :: DEF_dir_restart  = 'path/to/restart'
   CHARACTER(len=256) :: DEF_dir_history  = 'path/to/history'

#ifdef GRIDBASED
   CHARACTER(len=256) :: DEF_file_landgrid = 'path/to/landmask/file'
#endif

#ifdef CATCHMENT
   CHARACTER(len=256) :: DEF_dir_hydrodata = 'path/to/hydrodata'
   INTEGER            :: DEF_max_hband = 25
#endif 

   !add by zhongwang wei @ sysu 2021/12/23 
   !To allow read satellite observed LAI        
   logical :: DEF_LAI_CLIM = .FALSE.      
   INTEGER            :: DEF_Interception_scheme    = 1  !1:CoLM；2:CLM4.5; 3:CLM5; 4:Noah-MP; 5:MATSIRO; 6:VIC
                 

   ! ----- history -----
   REAL(r8) :: DEF_hist_lon_res = 0.5
   REAL(r8) :: DEF_hist_lat_res = 0.5       
   CHARACTER(len=256) :: DEF_WRST_FREQ    = 'none'  ! write restart file frequency: HOURLY/DAILY/MONTHLY/YEARLY
   CHARACTER(len=256) :: DEF_HIST_FREQ    = 'none'  ! write history file frequency: HOURLY/DAILY/MONTHLY/YEARLY
   CHARACTER(len=256) :: DEF_HIST_groupby = 'MONTH' ! history file in one file: DAY/MONTH/YEAR
   CHARACTER(len=256) :: DEF_HIST_mode    ='block'
   INTEGER :: DEF_REST_COMPRESS_LEVEL = 1
   INTEGER :: DEF_HIST_COMPRESS_LEVEL = 1

   ! ----- forcing -----
   TYPE nl_forcing_type

      CHARACTER(len=256) :: dataset            = 'CRUNCEP' 
      LOGICAL            :: solarin_all_band   = .true.  
      REAL(r8)           :: HEIGHT_V           = 100.0    
      REAL(r8)           :: HEIGHT_T           = 50.
      REAL(r8)           :: HEIGHT_Q           = 50.

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

   END TYPE nl_forcing_type

   TYPE (nl_forcing_type) :: DEF_forcing

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
#ifdef CROP
      LOGICAL :: cphase             = .true.
      LOGICAL :: cropprod1c         = .true.
      LOGICAL :: cropprod1c_loss    = .true.
      LOGICAL :: cropseedc_deficit  = .true.
      LOGICAL :: grainc_to_cropprodc= .true.
      LOGICAL :: grainc_to_seed     = .true.
#endif
#endif

      LOGICAL :: t_soisno     = .true. 
      LOGICAL :: wliq_soisno  = .true. 
      LOGICAL :: wice_soisno  = .true. 
                                       
      LOGICAL :: h2osoi       = .true. 
      LOGICAL :: rstfacsun    = .true.
      LOGICAL :: rstfacsha    = .true.
      LOGICAL :: rootr        = .true.
      LOGICAL :: vegwp        = .true.
#ifdef VARIABLY_SATURATED_FLOW
      LOGICAL :: dpond        = .true. 
#ifdef USE_DEPTH_TO_BEDROCK
      LOGICAL :: dwatsub      = .true. 
#endif
#endif
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

   END TYPE history_var_type

   TYPE (history_var_type) :: DEF_hist_vars

CONTAINS

   SUBROUTINE read_namelist (nlfile)

      USE spmd_task
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
         SITE_lon_location,       &
         SITE_lat_location,       &
         SITE_fsrfdata,           &
         SITE_landtype,           &
#if (defined PFT_CLASSIFICATION || defined PC_CLASSIFICATION)
         SITE_pct_pfts,           &
#ifdef CROP
         SITE_pctcrop,            &
#endif
#endif
         USE_SITE_htop,            &
         USE_SITE_LAI,             &
         USE_SITE_lakedepth,       &
         USE_SITE_soilreflectance, &
         USE_SITE_soilparameters,  &
#ifdef USE_DEPTH_TO_BEDROCK
         USE_SITE_dbedrock,       &
#endif
         USE_SITE_Forcing,        &
#endif
         DEF_nx_blocks,                   &
         DEF_ny_blocks,                   &
         DEF_PIO_groupsize,               &
         DEF_simulation_time,             &
         DEF_dir_rawdata,                 &  
         DEF_dir_output,                  &  
         DEF_dir_forcing,                 &  
#ifdef GRIDBASED
         DEF_file_landgrid,               &
#endif
#ifdef CATCHMENT
         DEF_dir_hydrodata,               &
         DEF_max_hband,                   &
#endif
         DEF_LAI_CLIM,                    &   !add by zhongwang wei @ sysu 2021/12/23        
         DEF_Interception_scheme,         &   !add by zhongwang wei @ sysu 2022/05/23    
        
         DEF_hist_lon_res,                &
         DEF_hist_lat_res,                &
         DEF_WRST_FREQ,                   &
         DEF_HIST_FREQ,                   &
         DEF_HIST_groupby,                &
         DEF_HIST_mode,                   &  
         DEF_REST_COMPRESS_LEVEL,         & 
         DEF_HIST_COMPRESS_LEVEL,         & 
         DEF_forcing,                     &
         DEF_hist_vars


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

         DEF_dir_landdata = trim(DEF_dir_output) // '/' // trim(adjustl(DEF_CASE_NAME)) // '/landdata'
         DEF_dir_restart  = trim(DEF_dir_output) // '/' // trim(adjustl(DEF_CASE_NAME)) // '/restart'
         DEF_dir_history  = trim(DEF_dir_output) // '/' // trim(adjustl(DEF_CASE_NAME)) // '/history'

         CALL system('mkdir -p ' // trim(adjustl(DEF_dir_output  )))
         CALL system('mkdir -p ' // trim(adjustl(DEF_dir_landdata)))
         CALL system('mkdir -p ' // trim(adjustl(DEF_dir_restart )))
         CALL system('mkdir -p ' // trim(adjustl(DEF_dir_history )))

#ifdef SinglePoint
         DEF_domain%edges = floor(SITE_lat_location)
         DEF_domain%edgen = DEF_domain%edges + 1.0
         DEF_domain%edgew = floor(SITE_lon_location)
         DEF_domain%edgee = DEF_domain%edgew + 1.0

         DEF_nx_blocks = 360
         DEF_ny_blocks = 180
#endif

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
      
      CALL mpi_bcast (DEF_simulation_time%greenwich,    1, mpi_logical, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_simulation_time%start_year,   1, mpi_integer, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_simulation_time%start_month,  1, mpi_integer, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_simulation_time%start_day,    1, mpi_integer, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_simulation_time%start_sec,    1, mpi_integer, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_simulation_time%end_year,     1, mpi_integer, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_simulation_time%end_month,    1, mpi_integer, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_simulation_time%end_day,      1, mpi_integer, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_simulation_time%end_sec,      1, mpi_integer, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_simulation_time%spinup_year,  1, mpi_integer, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_simulation_time%spinup_month, 1, mpi_integer, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_simulation_time%spinup_day,   1, mpi_integer, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_simulation_time%spinup_sec,   1, mpi_integer, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_simulation_time%timestep,     1, mpi_real8,   p_root, p_comm_glb, p_err)

      CALL mpi_bcast (DEF_dir_rawdata,  256, mpi_character, p_root, p_comm_glb, p_err)  
      CALL mpi_bcast (DEF_dir_output,   256, mpi_character, p_root, p_comm_glb, p_err)  
      CALL mpi_bcast (DEF_dir_forcing,  256, mpi_character, p_root, p_comm_glb, p_err)
      
      CALL mpi_bcast (DEF_dir_landdata, 256, mpi_character, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_dir_restart,  256, mpi_character, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_dir_history,  256, mpi_character, p_root, p_comm_glb, p_err)
      
#ifdef GRIDBASED
      CALL mpi_bcast (DEF_file_landgrid, 256, mpi_character, p_root, p_comm_glb, p_err)
#endif

#ifdef CATCHMENT
      CALL mpi_bcast (DEF_dir_hydrodata, 256, mpi_character, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_max_hband,       1, mpi_integer,   p_root, p_comm_glb, p_err)
#endif

      !zhongwang wei, 20210927: add option to read non-climatological mean LAI 
      call mpi_bcast (DEF_LAI_CLIM,        1, mpi_logical, p_root, p_comm_glb, p_err)
      !zhongwang wei, 20220520: add option to choose different canopy interception schemes
      call mpi_bcast (DEF_Interception_scheme, 1, mpi_integer, p_root, p_comm_glb, p_err)

      CALL mpi_bcast (DEF_hist_lon_res,  1, mpi_real8, p_root, p_comm_glb, p_err) 
      CALL mpi_bcast (DEF_hist_lat_res,  1, mpi_real8, p_root, p_comm_glb, p_err)

      CALL mpi_bcast (DEF_WRST_FREQ,         256, mpi_character, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_HIST_FREQ,         256, mpi_character, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_HIST_groupby,      256, mpi_character, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_HIST_mode,      256, mpi_character, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_REST_COMPRESS_LEVEL, 1, mpi_integer,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_HIST_COMPRESS_LEVEL, 1, mpi_integer,   p_root, p_comm_glb, p_err)

      CALL mpi_bcast (DEF_forcing%dataset,          256, mpi_character, p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_forcing%solarin_all_band,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_forcing%HEIGHT_V,           1, mpi_real8,     p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_forcing%HEIGHT_T,           1, mpi_real8,     p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_forcing%HEIGHT_Q,           1, mpi_real8,     p_root, p_comm_glb, p_err)
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

      CALL mpi_bcast (DEF_hist_vars%xy_us       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%xy_vs       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%xy_t        ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%xy_q        ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%xy_prc      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%xy_prl      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%xy_pbot     ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%xy_frl      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%xy_solarin  ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%xy_rain     ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%xy_snow     ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      
      CALL mpi_bcast (DEF_hist_vars%taux        ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%tauy        ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%fsena       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%lfevpa      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%fevpa       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%fsenl       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%fevpl       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%etr         ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%fseng       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%fevpg       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%fgrnd       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%sabvsun     ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%sabvsha     ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%sabg        ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%olrg        ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%rnet        ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%xerr        ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%zerr        ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%rsur        ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%rnof        ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%qintr       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%qinfl       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%qdrip       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%wat         ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%assim       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%respc       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%qcharge     ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%t_grnd      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%tleaf       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%ldew        ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%scv         ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%snowdp      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%fsno        ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%sigf        ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%green       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%lai         ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%laisun      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%laisha      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%sai         ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%alb         ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%emis        ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%z0m         ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%trad        ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%tref        ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%qref        ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
#ifdef BGC
      CALL mpi_bcast (DEF_hist_vars%leafc              ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%leafc_storage      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%leafc_xfer         ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%frootc             ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%frootc_storage     ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%frootc_xfer        ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%livestemc          ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%livestemc_storage  ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%livestemc_xfer     ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%deadstemc          ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%deadstemc_storage  ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%deadstemc_xfer     ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%livecrootc         ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%livecrootc_storage ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%livecrootc_xfer    ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%deadcrootc         ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%deadcrootc_storage ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%deadcrootc_xfer    ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%grainc             ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%grainc_storage     ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%grainc_xfer        ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%leafn              ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%leafn_storage      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%leafn_xfer         ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%frootn             ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%frootn_storage     ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%frootn_xfer        ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%livestemn          ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%livestemn_storage  ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%livestemn_xfer     ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%deadstemn          ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%deadstemn_storage  ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%deadstemn_xfer     ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%livecrootn         ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%livecrootn_storage ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%livecrootn_xfer    ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%deadcrootn         ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%deadcrootn_storage ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%deadcrootn_xfer    ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%grainn             ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%grainn_storage     ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%grainn_xfer        ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%retrasn            ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%gpp                ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%downreg            ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%ar                 ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
#ifdef CROP
      CALL mpi_bcast (DEF_hist_vars%cphase             ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%cropprod1c         ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%cropprod1c_loss    ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%cropseedc_deficit  ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%grainc_to_cropprodc,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%grainc_to_seed     ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
#endif
#endif
      
      CALL mpi_bcast (DEF_hist_vars%t_soisno    ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%wliq_soisno ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%wice_soisno ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      
      CALL mpi_bcast (DEF_hist_vars%h2osoi      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%rstfacsun   ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%rstfacsha   ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%rootr       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%vegwp       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
#ifdef VARIABLY_SATURATED_FLOW
      CALL mpi_bcast (DEF_hist_vars%dpond       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
#ifdef USE_DEPTH_TO_BEDROCK
      CALL mpi_bcast (DEF_hist_vars%dwatsub     ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
#endif
#endif
      CALL mpi_bcast (DEF_hist_vars%zwt         ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%wa          ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      
      CALL mpi_bcast (DEF_hist_vars%t_lake      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%lake_icefrac,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      
#ifdef BGC
      CALL mpi_bcast (DEF_hist_vars%litr1c_vr   ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%litr2c_vr   ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%litr3c_vr   ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%soil1c_vr   ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%soil2c_vr   ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%soil3c_vr   ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%cwdc_vr     ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%litr1n_vr   ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%litr2n_vr   ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%litr3n_vr   ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%soil1n_vr   ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%soil2n_vr   ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%soil3n_vr   ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%cwdn_vr     ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%sminn_vr    ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
#endif

      CALL mpi_bcast (DEF_hist_vars%ustar       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%tstar       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%qstar       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%zol         ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%rib         ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%fm          ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%fh          ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%fq          ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%us10m       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%vs10m       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%fm10m       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%sr          ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%solvd       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%solvi       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%solnd       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%solni       ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%srvd        ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%srvi        ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%srnd        ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%srni        ,   1, mpi_logical,   p_root, p_comm_glb, p_err)

      CALL mpi_bcast (DEF_hist_vars%solvdln     ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%solviln     ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%solndln     ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%solniln     ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%srvdln      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%srviln      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%srndln      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
      CALL mpi_bcast (DEF_hist_vars%srniln      ,   1, mpi_logical,   p_root, p_comm_glb, p_err)
#endif 

   END SUBROUTINE read_namelist

END MODULE mod_namelist
