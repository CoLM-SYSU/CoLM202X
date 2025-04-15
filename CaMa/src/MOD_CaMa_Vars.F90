#include <define.h>

MODULE MOD_CaMa_Vars
!DESCRIPTION
!===========
   !---This MODULE is the coupler for the colm and CaMa-Flood model.

!ANCILLARY FUNCTIONS AND SUBROUTINES
!-------------------
   !* :SUBROUTINE:"allocate_acc_cama_fluxes"   :  Initilization Accumulation of cama-flood variables
   !* :SUBROUTINE:"deallocate_acc_cama_fluxes" :  deallocate Accumulation of cama-flood variables
   !* :SUBROUTINE:"FLUSH_acc_cama_fluxes"      :  Reset Accumulation of cama-flood variables
   !* :SUBROUTINE:"accumulate_cama_fluxes"     :  Get accumulated cama-flood variables
   !* :SUBROUTINE:"allocate_2D_cama_Fluxes"    :  Get floodplain evaporation
   !* :SUBROUTINE:"hist_out_cama"              :  Average camaflood history variables and write out
   !* :SUBROUTINE:"hist_write_cama_time"       :  write out cama-flood history variables
   !* :SUBROUTINE:"flux_map_and_write_2d_cama" :  map cama output and write 1 or 2D variables, e.g. lon, lat
   !* :SUBROUTINE:"colm2cama_real8"            :  Send variables from worker processors to master processors
   !* :SUBROUTINE:"cama2colm_real8"            :  Send variables from master processors (cama) to worker processors (colm)

!REVISION HISTORY
!----------------
   !---2023.02.23  Zhongwang Wei @ SYSU
   !---2021.12.12  Zhongwang Wei @ SYSU
   !---2020.10.21  Zhongwang Wei @ SYSU

#if(defined CaMa_Flood)
   USE MOD_Precision
   USE MOD_Grid
   USE MOD_DataType

   USE MOD_SpatialMapping
   USE YOS_CMF_INPUT,            only: RMIS, DMIS
   USE MOD_Vars_Global,    only: spval


   real(r8) :: nacc                                        ! number of accumulation
   real(r8), allocatable         :: a_rnof_cama (:)        ! on worker : total runoff [mm/s]
   type(block_data_real8_2d)     :: f_rnof_cama            ! on IO     : total runoff [mm/s]
   real(r8), allocatable         :: runoff_2d (:,:)        ! on Master : total runoff [mm/s]

   real(r8), allocatable         :: flddepth_cama (:)      ! on worker : flddepth [m]
   type(block_data_real8_2d)     :: f_flddepth_cama        ! on IO     : flddepth [m]
   real(r8), allocatable         :: flddepth_tmp(:,:)

   real(r8), allocatable         :: fldfrc_cama (:)        ! on worker : flddepth [m]
   type(block_data_real8_2d)     :: f_fldfrc_cama          ! on IO     : flddepth [m]
   real(r8), allocatable         :: fldfrc_tmp (:,:)       ! on Master : total runoff [mm/s]

   real(r8), allocatable         :: fevpg_fld(:)           ! m/s
   real(r8), allocatable         :: finfg_fld(:)           ! m/s

   real(r8), allocatable         :: a_fevpg_fld (:)        ! on worker : flddepth [m]
   type(block_data_real8_2d)     :: f_fevpg_fld            ! on IO : total runoff [mm/s]
   real(r8), allocatable         :: fevpg_2d (:,:)         ! on Master : total runoff [mm/s]

   real(r8), allocatable         :: a_finfg_fld (:)        ! on worker : flddepth [m]
   type(block_data_real8_2d)     :: f_finfg_fld            ! on IO : total runoff [mm/s]
   real(r8), allocatable         :: finfg_2d (:,:)         ! on Master : total runoff [mm/s]
   type(grid_type) :: gcama

   type (spatial_mapping_type) :: mp2g_cama               ! mapping pset to grid
   type (spatial_mapping_type) :: mg2p_cama               ! mapping grid to pset

   type (grid_concat_type)       :: cama_gather            ! gather grid

   type(block_data_real8_2d)     :: IO_Effdepth            ! inundation to water depth [m]
   type(block_data_real8_2d)     :: IO_Effarea             ! inundation to water area [m2]

   type history_var_cama_type
      logical :: rivout       = .false.
      logical :: rivsto       = .false.
      logical :: rivdph       = .false.
      logical :: rivvel       = .false.
      logical :: fldout       = .false.
      logical :: fldsto       = .false.
      logical :: flddph       = .false.
      logical :: fldfrc       = .false.
      logical :: fldare       = .false.
      logical :: sfcelv       = .false.
      logical :: totout       = .false.
      logical :: outflw       = .false.
      logical :: totsto       = .false.
      logical :: storge       = .false.
      logical :: pthflw       = .false.
      logical :: pthout       = .false.
      logical :: maxflw       = .false.
      logical :: maxdph       = .false.
      logical :: maxsto       = .false.
      logical :: gwsto        = .false.
      logical :: gdwsto       = .false.
      logical :: gwout        = .false.
      logical :: gdwrtn       = .false.
      logical :: runoff       = .false.
      logical :: runoffsub    = .false.
      logical :: rofsfc       = .false.
      logical :: rofsub       = .false.
      logical :: damsto       = .false.
      logical :: daminf       = .false.
      logical :: wevap        = .false.
      logical :: winfilt      = .false.
      logical :: levsto       = .false.
      logical :: levdph       = .false.
      logical :: outflw_ocean = .false.
   END type history_var_cama_type

   type (history_var_cama_type) :: DEF_hist_cama_vars

   ! --- subroutines ---
   PUBLIC :: allocate_acc_cama_fluxes
   PUBLIC :: deallocate_acc_cama_fluxes
   PUBLIC :: flush_acc_cama_fluxes
   PUBLIC :: accumulate_cama_fluxes
   PUBLIC :: allocate_2D_cama_Fluxes
   PUBLIC :: colm2cama_real8
   PUBLIC :: cama2colm_real8
   PUBLIC :: hist_out_cama

CONTAINS

   SUBROUTINE allocate_acc_cama_fluxes
!DESCRIPTION
!===========
   ! This subrountine is used for initilization Accumulation of cama-flood variables

!ANCILLARY FUNCTIONS AND SUBROUTINES
!-------------------

!REVISION HISTORY
!----------------
   ! 2021.12.12  Zhongwang Wei @ SYSU


   USE MOD_SPMD_Task !spmd_task
   USE MOD_LandPatch, only: numpatch
   USE MOD_Vars_Global

   IMPLICIT NONE

      !allocate cama-flood variables on worker
      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            allocate (a_rnof_cama(numpatch))
            allocate (a_fevpg_fld(numpatch))
            allocate (a_finfg_fld(numpatch))
         ENDIF
      ENDIF

   END SUBROUTINE allocate_acc_cama_fluxes

   SUBROUTINE deallocate_acc_cama_fluxes()
!DESCRIPTION
!===========
! This subrountine is used for deallocate Accumulation of cama-flood variables

!ANCILLARY FUNCTIONS AND SUBROUTINES
   !-------------------

!REVISION HISTORY
   !----------------
   ! 2020.10.21  Zhongwang Wei @ SYSU

   USE MOD_SPMD_Task
   USE MOD_LandPatch, only: numpatch

   IMPLICIT NONE

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            deallocate (a_rnof_cama)
            deallocate (a_fevpg_fld)
            deallocate (a_finfg_fld)
         ENDIF
      ENDIF

   END SUBROUTINE deallocate_acc_cama_fluxes

   SUBROUTINE FLUSH_acc_cama_fluxes ()
!DESCRIPTION
!===========
! This subrountine is used for reset Accumulation of cama-flood variables

!ANCILLARY FUNCTIONS AND SUBROUTINES
   !-------------------

!REVISION HISTORY
   !----------------
   ! 2020.10.21  Zhongwang Wei @ SYSU
   USE MOD_SPMD_Task
   USE MOD_LandPatch,      only: numpatch
   USE MOD_Vars_Global,    only: spval

   IMPLICIT NONE

      IF (p_is_worker) THEN

         nacc = 0

         IF (numpatch > 0) THEN
            ! flush the Fluxes for accumulation
            a_rnof_cama (:) = spval
            a_fevpg_fld (:) = spval
            a_finfg_fld (:) = spval
         ENDIF
      ENDIF

   END SUBROUTINE FLUSH_acc_cama_fluxes

   SUBROUTINE accumulate_cama_fluxes
!DESCRIPTION
!===========
   ! This subrountine is used for accumulating  cama-flood variables

!ANCILLARY FUNCTIONS AND SUBROUTINES
   !-------------------
   !* :SUBROUTINE:"acc1d_cama"            :  accumulating 1D cama-flood variables

!REVISION HISTORY
   !----------------
   ! 2020.10.21  Zhongwang Wei @ SYSU

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Vars_1DFluxes, only: rnof
   USE MOD_LandPatch, only: numpatch

   IMPLICIT NONE

      IF (p_is_worker) THEN
         IF (numpatch > 0) THEN
            nacc = nacc + 1
            CALL acc1d_cama (rnof, a_rnof_cama)
            CALL acc1d_cama (fevpg_fld, a_fevpg_fld)
            CALL acc1d_cama (finfg_fld, a_finfg_fld)
         ENDIF
      ENDIF

   END SUBROUTINE accumulate_cama_fluxes

   SUBROUTINE acc1d_cama (var, s)
!DESCRIPTION
!===========
   ! This subrountine is used for accumulating 1D cama-flood variables

!ANCILLARY FUNCTIONS AND SUBROUTINES
   !-------------------

!REVISION HISTORY
   !----------------
   ! 2020.10.21  Zhongwang Wei @ SYSU

   USE MOD_Precision
   USE MOD_Vars_Global, only: spval

   IMPLICIT NONE

      real(r8), intent(in)    :: var(:) ! variable to be accumulated
      real(r8), intent(inout) :: s  (:) ! new added value
!----------------------- Dummy argument --------------------------------
      integer :: i

      DO i = lbound(var,1), ubound(var,1)
         IF (var(i) /= spval) THEN
            IF (s(i) /= spval) THEN
               s(i) = s(i) + var(i)
            ELSE
               s(i) = var(i)
            ENDIF
         ENDIF
      ENDDO

   END SUBROUTINE acc1d_cama

   SUBROUTINE allocate_2D_cama_Fluxes (grid)
!DESCRIPTION
!===========
   ! This subrountine is used for accumulating 2D cama-flood variables

!ANCILLARY FUNCTIONS AND SUBROUTINES
!-------------------
   !* :SUBROUTINE:"allocate_block_data"            :  allocate 2D cama-flood variables to colm block

!REVISION HISTORY
   !----------------
   ! 2020.10.21  Zhongwang Wei @ SYSU

   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_DataType

   IMPLICIT NONE

      type(grid_type), intent(in) :: grid

      IF (p_is_io) THEN
         CALL allocate_block_data (grid, f_rnof_cama)      ! total runoff         [m/s]
         CALL allocate_block_data (grid, f_flddepth_cama)  ! inundation depth     [m/s]
         CALL allocate_block_data (grid, f_fldfrc_cama)    ! inundation fraction  [m/s]
         !TODO: check the following variables
         CALL allocate_block_data (grid, f_fevpg_fld)      ! inundation evaporation [m/s]
         CALL allocate_block_data (grid, f_finfg_fld)      ! inundation re-infiltration [m/s]
      ENDIF

   END SUBROUTINE allocate_2D_cama_Fluxes

   SUBROUTINE hist_out_cama (file_hist, itime_in_file)
!DESCRIPTION
!===========
   ! This subrountine is used for averaging and writing 2D cama-flood variables out

!ANCILLARY FUNCTIONS AND SUBROUTINES
   !-------------------
   !* :SUBROUTINE:"CMF_DIAG_AVERAGE"                      :  averaging the diagnostic variables of cama-flood
   !* :SUBROUTINE:"flux_map_and_write_2d_cama"            :  map camaflood variables to colm block and write out
   !* :SUBROUTINE:"CMF_DIAG_RESET"                        :  reset diagnostic variables of cama-flood
!REVISION HISTORY
   !----------------
   ! 2020.10.21  Zhongwang Wei @ SYSU

   USE MOD_SPMD_Task
   USE CMF_CALC_DIAG_MOD,  only: CMF_DIAG_AVERAGE, CMF_DIAG_RESET
   USE YOS_CMF_PROG,       only: P2RIVSTO,     P2FLDSTO,     P2GDWSTO, &
         P2damsto,P2LEVSTO !!! added
   USE YOS_CMF_DIAG,       only: D2RIVDPH,     D2FLDDPH,     D2FLDFRC,     D2FLDARE,     &
         D2SFCELV,     D2STORGE,                                                         &
         D2OUTFLW_AVG, D2RIVOUT_AVG, D2FLDOUT_AVG, D2PTHOUT_AVG, D1PTHFLW_AVG,           &
         D2RIVVEL_AVG, D2GDWRTN_AVG, D2RUNOFF_AVG, D2ROFSUB_AVG,                         &
         D2OUTFLW_MAX, D2STORGE_MAX, D2RIVDPH_MAX,                                       &
         d2daminf_avg, D2WEVAPEX_AVG,D2WINFILTEX_AVG, D2LEVDPH !!! added
 !      USE MOD_Vars_2DFluxes

   IMPLICIT NONE

   character(LEN=*), intent(in) :: file_hist
   integer, intent(in)          :: itime_in_file


      !*** average variable
      CALL CMF_DIAG_AVERAGE

      !*** write output data
      CALL flux_map_and_write_2d_cama (DEF_hist_cama_vars%rivout, &
      real(D2RIVOUT_AVG), file_hist, 'rivout', itime_in_file,'river discharge','m3/s')

      CALL flux_map_and_write_2d_cama (DEF_hist_cama_vars%rivsto, &
      real(P2RIVSTO), file_hist, 'rivsto', itime_in_file,'river storage','m3')

      CALL flux_map_and_write_2d_cama (DEF_hist_cama_vars%rivdph, &
      real(D2RIVDPH), file_hist, 'rivdph', itime_in_file,'river depth','m')

      CALL flux_map_and_write_2d_cama (DEF_hist_cama_vars%rivvel, &
      real(D2RIVVEL_AVG), file_hist, 'rivvel', itime_in_file,'river velocity','m/s')

      CALL flux_map_and_write_2d_cama(DEF_hist_cama_vars%fldout, &
      real(D2FLDOUT_AVG), file_hist, 'fldout', itime_in_file,'floodplain discharge','m3/s')

      CALL flux_map_and_write_2d_cama(DEF_hist_cama_vars%fldsto, &
      real(P2FLDSTO), file_hist, 'fldsto', itime_in_file,'floodplain storage','m3')

      CALL flux_map_and_write_2d_cama(DEF_hist_cama_vars%flddph, &
      real(D2FLDDPH), file_hist, 'flddph', itime_in_file,'floodplain depth','m')

      CALL flux_map_and_write_2d_cama(DEF_hist_cama_vars%fldfrc, &
      real(D2FLDFRC), file_hist, 'fldfrc', itime_in_file,'flooded fraction','0-1')

      CALL flux_map_and_write_2d_cama(DEF_hist_cama_vars%fldare, &
      real(D2FLDARE), file_hist, 'fldare', itime_in_file,'flooded area','m2')

      CALL flux_map_and_write_2d_cama(DEF_hist_cama_vars%sfcelv, &
      real(D2SFCELV), file_hist, 'sfcelv', itime_in_file,'water surface elevation','m')

      CALL flux_map_and_write_2d_cama(DEF_hist_cama_vars%totout, &
      real(D2OUTFLW_AVG), file_hist, 'totout', itime_in_file,'discharge (river+floodplain)','m3/s')

      CALL flux_map_and_write_2d_cama(DEF_hist_cama_vars%outflw, &
      real(D2OUTFLW_AVG), file_hist, 'outflw', itime_in_file,'discharge (river+floodplain)','m3/s')

      CALL flux_map_and_write_2d_cama(DEF_hist_cama_vars%totsto, &
      real(D2STORGE), file_hist, 'totsto', itime_in_file,'total storage (river+floodplain)','m3')

      CALL flux_map_and_write_2d_cama(DEF_hist_cama_vars%storge, &
      real(D2STORGE), file_hist, 'storge', itime_in_file,'total storage (river+floodplain)','m3')

      CALL flux_map_and_write_2d_cama(DEF_hist_cama_vars%pthflw, &
      real(D1PTHFLW_AVG), file_hist, 'pthflw', itime_in_file,'bifurcation channel discharge ','m3/s')

      CALL flux_map_and_write_2d_cama(DEF_hist_cama_vars%pthout, &
      real(D2PTHOUT_AVG), file_hist, 'pthout', itime_in_file,'net bifurcation discharge','m3/s')

      CALL flux_map_and_write_2d_cama(DEF_hist_cama_vars%gdwsto, &
      real(P2GDWSTO), file_hist, 'gdwsto', itime_in_file,'ground water storage','m3')

      CALL flux_map_and_write_2d_cama(DEF_hist_cama_vars%gwsto, &
      real(P2GDWSTO), file_hist, 'gwsto', itime_in_file,'ground water storage','m3')

      CALL flux_map_and_write_2d_cama(DEF_hist_cama_vars%gwout, &
      real(D2GDWRTN_AVG), file_hist, 'gwout', itime_in_file,'ground water discharge','m3/s')

      CALL flux_map_and_write_2d_cama(DEF_hist_cama_vars%runoff, &
      real(D2RUNOFF_AVG), file_hist, 'runoff', itime_in_file,'Surface runoff','m3/s')

      CALL flux_map_and_write_2d_cama(DEF_hist_cama_vars%runoffsub, &
      real(D2ROFSUB_AVG)  , file_hist, 'runoffsub', itime_in_file,'sub-surface runoff','m3/s')

      CALL flux_map_and_write_2d_cama(DEF_hist_cama_vars%maxsto, &
      real(D2STORGE_MAX), file_hist, 'maxsto', itime_in_file,'daily maximum storage','m3')

      CALL flux_map_and_write_2d_cama(DEF_hist_cama_vars%maxflw, &
      real(D2OUTFLW_MAX), file_hist, 'maxflw', itime_in_file,'daily maximum discharge','m3/s')

      CALL flux_map_and_write_2d_cama(DEF_hist_cama_vars%maxdph, &
      real(D2RIVDPH_MAX), file_hist, 'maxdph', itime_in_file,'daily maximum river depth','m')

      IF (DEF_hist_cama_vars%damsto) THEN
      CALL flux_map_and_write_2d_cama(DEF_hist_cama_vars%damsto, &
      real(p2damsto), file_hist, 'damsto', itime_in_file,'reservoir storage','m3')
      ENDIF

      IF (DEF_hist_cama_vars%daminf) THEN
      CALL flux_map_and_write_2d_cama(DEF_hist_cama_vars%daminf, &
      real(d2daminf_avg), file_hist, 'daminf', itime_in_file,'reservoir inflow','m3/s')
      ENDIF

      IF (DEF_hist_cama_vars%wevap) THEN
      CALL flux_map_and_write_2d_cama(DEF_hist_cama_vars%wevap, &
      real(D2WEVAPEX_AVG), file_hist, 'wevap', itime_in_file,'inundation water evaporation','m/s')
      ENDIF

      IF (DEF_hist_cama_vars%winfilt) THEN
      CALL flux_map_and_write_2d_cama(DEF_hist_cama_vars%winfilt, &
      real(D2WINFILTEX_AVG), file_hist, 'winfilt', itime_in_file,'inundation water infiltration','m/s')
      ENDIF

      IF (DEF_hist_cama_vars%levsto) THEN
      CALL flux_map_and_write_2d_cama(DEF_hist_cama_vars%levsto, &
      real(P2LEVSTO), file_hist, 'levsto', itime_in_file,'protected area storage','m3')
      ENDIF

      IF (DEF_hist_cama_vars%levdph) THEN
      CALL flux_map_and_write_2d_cama(DEF_hist_cama_vars%levdph, &
      real(D2LEVDPH), file_hist, 'levdph', itime_in_file,'protected area depth','m')
      ENDIF

      IF (DEF_hist_cama_vars%outflw_ocean) THEN
      CALL flux_map_and_write_2d_cama_ocean(DEF_hist_cama_vars%outflw_ocean, &
      real(D2OUTFLW_AVG), file_hist, 'outflw_ocean', itime_in_file,'discharge to ocean','m3/s')
      ENDIF
      
      !*** reset variable
      CALL CMF_DIAG_RESET

   END SUBROUTINE hist_out_cama


   SUBROUTINE hist_write_cama_time (filename, dataname, time, itime)
!DESCRIPTION
!===========
   ! This subrountine is used for writing time,longitude and latitude of cama-flood output using netcdf format.

!ANCILLARY FUNCTIONS AND SUBROUTINES
!-------------------
   !* :SUBROUTINE:"ncio_create_file"                      :  create netcdf file, see ncio_serial.F90
   !* :SUBROUTINE:"ncio_define_dimension"                 :  define dimension of netcdf file, see ncio_serial.F90
   !* :SUBROUTINE:"ncio_write_serial"                     :  write serial data into netcdf file (lon, lat), see ncio_serial.F90
   !* :SUBROUTINE:"ncio_write_time"                       :  write time serial into netcdf file (lon, lat), see ncio_serial.F90

!REVISION HISTORY
!----------------
   ! 2023.02.23  Zhongwang Wei @ SYSU

   USE MOD_SPMD_Task
   USE MOD_NetCDFSerial
   USE YOS_CMF_INPUT, only: NX, NY
   USE YOS_CMF_MAP,   only: D1LON, D1LAT

   IMPLICIT NONE

   character (len=*), intent(in) :: filename ! file name
   character (len=*), intent(in) :: dataname ! data name
   integer, intent(in)  :: time(3)           ! time (year, month, day)
   integer, intent(out) :: itime             ! number of time step

   ! Local variables
   logical :: fexists

      IF (p_is_master) THEN
         inquire (file=filename, exist=fexists)
         IF (.not. fexists) THEN
            CALL ncio_create_file (trim(filename))
            CALL ncio_define_dimension(filename, 'time', 0)
            CALL ncio_define_dimension(filename,'lat_cama', NY)
            CALL ncio_define_dimension(filename,'lon_cama', NX)
            CALL ncio_write_serial (filename, 'lat_cama', D1LAT,'lat_cama')
            CALL ncio_write_serial (filename, 'lon_cama', D1LON,'lon_cama')
         ENDIF

         CALL ncio_write_time (filename, dataname, time, itime)
      ENDIF

   END SUBROUTINE hist_write_cama_time

   SUBROUTINE flux_map_and_write_2d_cama_ocean (is_hist, &
         var_in, file_hist, varname, itime_in_file,longname,units)
         !DESCRIPTION
         !===========
            ! This subrountine is used for mapping cama-flood output using netcdf format.

         !ANCILLARY FUNCTIONS AND SUBROUTINES
         !-------------------
            !* :SUBROUTINE:"ncio_put_attr"                      :  write netcdf attribute, see ncio_serial.F90
            !* :SUBROUTINE:"vecP2mapR"                          :  convert 1D vector data -> 2D map data (real*4), CAMA/cmf_utils_mod.F90
            !* :SUBROUTINE:"ncio_write_serial_time"             :  define dimension of netcdf file, see ncio_serial.F90

         !REVISION HISTORY
         !----------------
            ! 2023.02.23  Zhongwang Wei @ SYSU

         USE MOD_Namelist
         USE YOS_CMF_INPUT,  only: NX, NY
         USE YOS_CMF_MAP,    only: NSEQALL
         USE PARKIND1,       only: JPRM
         USE CMF_UTILS_MOD,  only: vecP2mapR
         USE MOD_NetCDFSerial,    only: ncio_write_serial_time, ncio_put_attr
         USE YOS_CMF_MAP,        only: I2NEXTX, I2NEXTY
         IMPLICIT NONE
         logical, intent(in)          :: is_hist
         real(r8), intent(in)         ::  var_in (NSEQALL, 1)
         character(len=*), intent(in) :: file_hist
         character(len=*), intent(in) :: varname
         integer, intent(in)          :: itime_in_file
      
         character (len=*), intent(in),optional :: longname
         character (len=*), intent(in),optional :: units
      
         real(KIND=JPRM)             :: R2OUT(NX,NY)
         integer  :: i,j
         integer  :: compress
      
            IF (.not. is_hist) RETURN
            CALL vecP2mapR(var_in,R2OUT)
            do j = 1,NY
               do i = 1,NX
                  if (I2NEXTX(i,j) .NE. -9) then
                     R2OUT(i,j) = real(real(spval,kind=JPRM),kind=8)
                  endif
               enddo
            enddo
 
            compress = DEF_HIST_CompressLevel
            CALL ncio_write_serial_time (file_hist, varname,  &
               itime_in_file, real(R2OUT,kind=8), 'lon_cama', 'lat_cama', 'time',compress)
            IF (itime_in_file == 1) THEN
               CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
               CALL ncio_put_attr (file_hist, varname, 'units', units)
               CALL ncio_put_attr (file_hist, varname, 'missing_value',real(real(spval,kind=JPRM),kind=8))
            ENDIF
   end subroutine flux_map_and_write_2d_cama_ocean


   SUBROUTINE flux_map_and_write_2d_cama (is_hist, &
         var_in, file_hist, varname, itime_in_file,longname,units)

!DESCRIPTION
!===========
   ! This subrountine is used for mapping cama-flood output using netcdf format.

!ANCILLARY FUNCTIONS AND SUBROUTINES
!-------------------
   !* :SUBROUTINE:"ncio_put_attr"                      :  write netcdf attribute, see ncio_serial.F90
   !* :SUBROUTINE:"vecP2mapR"                          :  convert 1D vector data -> 2D map data (real*4), CAMA/cmf_utils_mod.F90
   !* :SUBROUTINE:"ncio_write_serial_time"             :  define dimension of netcdf file, see ncio_serial.F90

!REVISION HISTORY
!----------------
   ! 2023.02.23  Zhongwang Wei @ SYSU

   USE MOD_Namelist
   USE YOS_CMF_INPUT,  only: NX, NY
   USE YOS_CMF_MAP,    only: NSEQALL
   USE PARKIND1,       only: JPRM
   USE CMF_UTILS_MOD,  only: vecP2mapR
   USE MOD_NetCDFSerial,    only: ncio_write_serial_time, ncio_put_attr

   IMPLICIT NONE
   logical, intent(in)          :: is_hist
   real(r8), intent(in)         ::  var_in (NSEQALL, 1)
   character(len=*), intent(in) :: file_hist
   character(len=*), intent(in) :: varname
   integer, intent(in)          :: itime_in_file
   
   character (len=*), intent(in),optional :: longname
   character (len=*), intent(in),optional :: units

   real(KIND=JPRM)             :: R2OUT(NX,NY)

   integer  :: compress

      IF (.not. is_hist) RETURN

      CALL vecP2mapR(var_in,R2OUT)
      compress = DEF_HIST_CompressLevel
      CALL ncio_write_serial_time (file_hist, varname,  &
         itime_in_file, real(R2OUT,kind=8), 'lon_cama', 'lat_cama', 'time',compress)
      IF (itime_in_file == 1) THEN
         CALL ncio_put_attr (file_hist, varname, 'long_name', longname)
         CALL ncio_put_attr (file_hist, varname, 'units', units)
         CALL ncio_put_attr (file_hist, varname, 'missing_value',real(real(spval,kind=JPRM),kind=8))
      ENDIF

   END SUBROUTINE flux_map_and_write_2d_cama

   SUBROUTINE colm2cama_real8 (WorkerVar, IOVar, MasterVar)
!DESCRIPTION
!===========
   ! This subrountine is used for mapping colm output to cama input.

!ANCILLARY FUNCTIONS AND SUBROUTINES
!-------------------
   !* :SUBROUTINE:"allocate_block_data"                      :  allocate data into block

!REVISION HISTORY
!----------------
   ! 2023.02.23  Zhongwang Wei @ SYSU

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_TimeManager
   USE MOD_SPMD_Task
   USE MOD_Block
   USE MOD_DataType
   USE MOD_LandPatch
   USE MOD_Vars_TimeInvariants, only: patchtype
   USE MOD_Forcing, only: forcmask_pch

   IMPLICIT NONE

   real(r8),                  intent(inout) :: WorkerVar(:)    !varialbe on worker processer
   type(block_data_real8_2d), intent(inout) :: IOVar           !varialbe on IO processer
   real(r8),                  intent(inout) :: MasterVar(:,:)  !varialbe on master processer

   type(block_data_real8_2d) :: sumwt                          !sum of weight
   logical,  allocatable     :: filter(:)                      !filter for patchtype
   !----------------------- Dummy argument --------------------------------
   integer :: xblk, yblk, xloc, yloc
   integer :: iblk, jblk, idata, ixseg, iyseg
   integer :: rmesg(3), smesg(3), isrc
   integer :: xdsp, ydsp, xcnt, ycnt
   real(r8), allocatable :: rbuf(:,:), sbuf(:,:), vdata(:,:)

      IF(p_is_master)THEN
         MasterVar(:,:) = spval
      ENDIF

      IF (p_is_worker) THEN
         WHERE (WorkerVar /= spval)
            WorkerVar = WorkerVar / nacc
         endwhere

         IF (numpatch > 0) THEN
            allocate (filter (numpatch))

            filter(:) = patchtype < 99
            IF (DEF_forcing%has_missing_value) THEN
               filter = filter .and. forcmask_pch
            ENDIF
         ENDIF
      ENDIF

      CALL mp2g_cama%pset2grid (WorkerVar, IOVar, spv = spval, msk = filter)

      IF (p_is_io) CALL allocate_block_data (gcama, sumwt)
      CALL mp2g_cama%get_sumarea (sumwt, filter)


      IF (p_is_io) THEN
         DO yblk = 1, gblock%nyblk
            DO xblk = 1, gblock%nxblk
               IF (gblock%pio(xblk,yblk) == p_iam_glb) THEN
                  DO yloc = 1, gcama%ycnt(yblk)
                     DO xloc = 1, gcama%xcnt(xblk)

                        IF (sumwt%blk(xblk,yblk)%val(xloc,yloc) > 0.00001) THEN
                           IF (IOVar%blk(xblk,yblk)%val(xloc,yloc) /= spval) THEN
                              IOVar%blk(xblk,yblk)%val(xloc,yloc) &
                                 = IOVar%blk(xblk,yblk)%val(xloc,yloc) &
                                 / sumwt%blk(xblk,yblk)%val(xloc,yloc)
                           ENDIF
                        ELSE
                           IOVar%blk(xblk,yblk)%val(xloc,yloc) = spval
                        ENDIF

                     ENDDO
                  ENDDO

               ENDIF
            ENDDO
         ENDDO
      ENDIF

      IF (p_is_master) THEN
         DO idata = 1, cama_gather%ndatablk
            CALL mpi_recv (rmesg, 3, MPI_INTEGER, MPI_ANY_SOURCE, 10011, p_comm_glb, p_stat, p_err)
            isrc  = rmesg(1)
            ixseg = rmesg(2)
            iyseg = rmesg(3)

            xdsp = cama_gather%xsegs(ixseg)%gdsp
            ydsp = cama_gather%ysegs(iyseg)%gdsp
            xcnt = cama_gather%xsegs(ixseg)%cnt
            ycnt = cama_gather%ysegs(iyseg)%cnt

            allocate (rbuf(xcnt,ycnt))
            CALL mpi_recv (rbuf, xcnt * ycnt, MPI_DOUBLE, &
               isrc, 10011, p_comm_glb, p_stat, p_err)
            MasterVar (xdsp+1:xdsp+xcnt,ydsp+1:ydsp+ycnt) = rbuf
            deallocate (rbuf)
         ENDDO

      ELSEIF (p_is_io) THEN
         DO iyseg = 1, cama_gather%nyseg
            DO ixseg = 1, cama_gather%nxseg

               iblk = cama_gather%xsegs(ixseg)%blk
               jblk = cama_gather%ysegs(iyseg)%blk

               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                  xdsp = cama_gather%xsegs(ixseg)%bdsp
                  ydsp = cama_gather%ysegs(iyseg)%bdsp
                  xcnt = cama_gather%xsegs(ixseg)%cnt
                  ycnt = cama_gather%ysegs(iyseg)%cnt

                  allocate (sbuf (xcnt,ycnt))
                  sbuf = IOVar%blk(iblk,jblk)%val(xdsp+1:xdsp+xcnt,ydsp+1:ydsp+ycnt)

                  smesg = (/p_iam_glb, ixseg, iyseg/)
                  CALL mpi_send (smesg, 3, MPI_INTEGER, &
                     p_address_master, 10011, p_comm_glb, p_err)
                  CALL mpi_send (sbuf, xcnt*ycnt, MPI_DOUBLE, &
                     p_address_master, 10011, p_comm_glb, p_err)

                  deallocate (sbuf)

               ENDIF
            ENDDO
         ENDDO
      ENDIF

      IF (allocated(filter)) deallocate(filter)

   END SUBROUTINE colm2cama_real8

   SUBROUTINE cama2colm_real8 (MasterVar, IOVar, WorkerVar)

!DESCRIPTION
!===========
   ! This subrountine is used for mapping cama-flood output to colm input

!ANCILLARY FUNCTIONS AND SUBROUTINES
!-------------------
   !* :SUBROUTINE:"mg2p_cama%map_aweighted"                 :  mapping grid to pset_type

!REVISION HISTORY
   !----------------
   ! 2023.02.23  Zhongwang Wei @ SYSU
   ! 2022.?      Zhongwang Wei and ShuPeng Zhang @ SYSU

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_TimeManager
   USE MOD_SPMD_Task
   USE MOD_Block
   USE MOD_DataType
   USE MOD_LandPatch
   USE MOD_Vars_TimeInvariants, only: patchtype
   USE MOD_Grid

   IMPLICIT NONE

   real(r8),                  intent(in)    :: MasterVar (:,:) ! Variable at master processor
   type(block_data_real8_2d), intent(inout) :: IOVar           ! Variable at io processor
   real(r8),                  intent(inout) :: WorkerVar (:)   ! Variable at worker processor

   integer :: xblk    , yblk    , xloc ,  yloc
   integer :: iblk    , jblk    , idata,  ixseg,  iyseg
   integer :: rmesg(2), smesg(2), isrc ,  iproc
   integer :: xdsp    , ydsp    , xcnt ,  ycnt
   real(r8), allocatable :: rbuf(:,:), sbuf(:,:), vdata(:,:)


      IF (p_is_master) THEN
         DO iyseg = 1, cama_gather%nyseg
            DO ixseg = 1, cama_gather%nxseg
               iblk = cama_gather%xsegs(ixseg)%blk
               jblk = cama_gather%ysegs(iyseg)%blk

               IF (gblock%pio(iblk,jblk) >= 0) THEN
                  xdsp = cama_gather%xsegs(ixseg)%gdsp
                  ydsp = cama_gather%ysegs(iyseg)%gdsp
                  xcnt = cama_gather%xsegs(ixseg)%cnt
                  ycnt = cama_gather%ysegs(iyseg)%cnt

                  allocate (sbuf (xcnt,ycnt))
                  sbuf = MasterVar (xdsp+1:xdsp+xcnt, ydsp+1:ydsp+ycnt)
                  smesg = (/ixseg, iyseg/)
                  CALL mpi_send (smesg, 2, MPI_INTEGER, &
                     gblock%pio(iblk,jblk), 10000, p_comm_glb, p_err)
                  CALL mpi_send (sbuf, xcnt*ycnt, MPI_DOUBLE, &
                     gblock%pio(iblk,jblk), 10000, p_comm_glb, p_err)
                  deallocate (sbuf)
               ENDIF
            ENDDO
         ENDDO

         DO iproc = 0, p_np_io-1
            smesg = (/0, 0/)
            CALL mpi_send(smesg, 2, MPI_INTEGER, p_address_io(iproc), 10000, p_comm_glb, p_err)
         ENDDO
      ELSEIF  (p_is_io) THEN
         DO WHILE (.true.)
            CALL mpi_recv (rmesg, 2, MPI_INTEGER, p_address_master, 10000, p_comm_glb, p_stat, p_err)
            ixseg = rmesg(1)
            iyseg = rmesg(2)

            IF ((ixseg > 0) .and. (iyseg > 0)) THEN
               iblk = cama_gather%xsegs(ixseg)%blk
               jblk = cama_gather%ysegs(iyseg)%blk
               xdsp = cama_gather%xsegs(ixseg)%bdsp
               ydsp = cama_gather%ysegs(iyseg)%bdsp
               xcnt = cama_gather%xsegs(ixseg)%cnt
               ycnt = cama_gather%ysegs(iyseg)%cnt

               allocate (rbuf(xcnt,ycnt))
               CALL mpi_recv (rbuf, xcnt*ycnt, MPI_DOUBLE, &
                  p_address_master, 10000, p_comm_glb, p_stat, p_err)
               IOVar%blk(iblk,jblk)%val(xdsp+1:xdsp+xcnt,ydsp+1:ydsp+ycnt)= rbuf
               deallocate (rbuf)
            ELSE
               EXIT
            ENDIF
         ENDDO
      ENDIF

      CALL mg2p_cama%grid2pset (IOVar, WorkerVar) !mapping grid to pset_type

   END SUBROUTINE cama2colm_real8

#endif

END MODULE MOD_CaMa_Vars
! ----- EOP ---------
