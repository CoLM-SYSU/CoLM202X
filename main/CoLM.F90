#include <define.h>

PROGRAM CoLM
   ! ======================================================================
   ! Reference:
   !     [1] Dai et al., 2003: The Common Land Model (CoLM).
   !         Bull. of Amer. Meter. Soc., 84: 1013-1023
   !     [2] Dai et al., 2004: A two-big-leaf model for canopy temperature,
   !         photosynthesis and stomatal conductance. J. Climate, 17: 2281-2299.
   !     [3] Dai et al., 2014: The Terrestrial Modeling System (TMS).
   !     [4] Dai Yamazaki, 2014: The global river model CaMa-Flood (version 3.6.2)
   !
   !     Created by Yongjiu Dai, Februay 2004
   !     Revised by Yongjiu Dai and Hua Yuan, April 2014
   ! ======================================================================

   use precision
   use spmd_task
   use mod_namelist
   USE MOD_Vars_Global
   USE MOD_Vars_LCConst
   USE MOD_Vars_PFTConst
   use MOD_Vars_PhysicalConst
   use MOD_Vars_TimeInvariants
   use MOD_Vars_TimeVariables
   use MOD_Vars_1DForcing
   use MOD_Vars_2DForcing
   use MOD_Vars_1DFluxes
   use MOD_Vars_2DFluxes
   use MOD_Forcing
   use MOD_Hist
   use timemanager
   use mod_colm_debug

   use mod_block
   use mod_pixel
   USE mod_mesh
   use mod_landelm
#ifdef CATCHMENT
   USE mod_landhru
#endif
   use mod_landpatch
#ifdef URBAN_MODEL
   USE mod_landurban
   USE MOD_Urban_LAIReadin
#endif
#ifdef PFT_CLASSIFICATION
   USE mod_landpft
#endif
#ifdef PC_CLASSIFICATION
   USE mod_landpc
#endif
#if (defined UNSTRUCTURED || defined CATCHMENT)
   USE mod_elm_vector
#endif
#ifdef CATCHMENT
   USE mod_hru_vector
#endif
#if(defined CaMa_Flood)
   use MOD_CaMa_colmCaMa ! whether cama-flood is used
#endif
#ifdef SinglePoint
   USE mod_single_srfdata
#endif

#if (defined LATERAL_FLOW)
   USE mod_lateral_flow
#endif

#ifdef Fire
   USE MOD_LightningData, only: init_lightning_data, update_lightning_data
   USE MOD_FireReadin, only: Fire_readin
#endif

#ifdef OzoneData
   USE MOD_OzoneData, only: init_ozone_data, update_ozone_data
#endif
   use mod_srfdata_restart
   USE MOD_LAIReadin
   USE MOD_NitrifReadin

#ifdef BGC
   USE MOD_NdepReadin
#endif

   ! SNICAR
   USE MOD_SnowSnicar , only:SnowAge_init, SnowOptics_init

   IMPLICIT NONE

   character(LEN=256) :: nlfile
   character(LEN=256) :: casename
   character(len=256) :: dir_landdata
   character(len=256) :: dir_forcing
   character(len=256) :: dir_hist
   character(len=256) :: dir_restart
   character(len=256) :: fsrfdata

   real(r8) :: deltim       ! time step (senconds)
   integer  :: sdate(3)     ! calendar (year, julian day, seconds)
   integer  :: idate(3)     ! calendar (year, julian day, seconds)
   integer  :: edate(3)     ! calendar (year, julian day, seconds)
   integer  :: pdate(3)     ! calendar (year, julian day, seconds)
   logical  :: greenwich    ! greenwich time

   logical :: doalb         ! true => start up the surface albedo calculation
   logical :: dolai         ! true => start up the time-varying vegetation paramter
   logical :: dosst         ! true => update sst/ice/snow

   integer :: Julian_1day_p, Julian_1day
   integer :: Julian_8day_p, Julian_8day
   integer :: s_year, s_month, s_day, s_seconds, s_julian
   integer :: e_year, e_month, e_day, e_seconds, e_julian
   integer :: p_year, p_month, p_day, p_seconds, p_julian
   INTEGER :: month, mday, month_p, mday_p
   INTEGER :: spinup_repeat

   type(timestamp) :: ststamp, itstamp, etstamp, ptstamp

   integer*8 :: start_time, end_time, c_per_sec, time_used
   logical isread

#ifdef USEMPI
   call spmd_init ()
#endif

   if (p_is_master) then
      call system_clock (start_time)
   end if

   call getarg (1, nlfile)

   call read_namelist (nlfile)

   casename     = DEF_CASE_NAME
   dir_landdata = DEF_dir_landdata
   dir_forcing  = DEF_dir_forcing
   dir_hist     = DEF_dir_history
   dir_restart  = DEF_dir_restart

#ifdef SinglePoint
   fsrfdata = trim(dir_landdata) // '/srfdata.nc'
   CALL read_surface_data_single (fsrfdata, mksrfdata=.false.)
#endif

   deltim    = DEF_simulation_time%timestep
   greenwich = DEF_simulation_time%greenwich
   s_year    = DEF_simulation_time%start_year
   s_month   = DEF_simulation_time%start_month
   s_day     = DEF_simulation_time%start_day
   s_seconds = DEF_simulation_time%start_sec
   e_year    = DEF_simulation_time%end_year
   e_month   = DEF_simulation_time%end_month
   e_day     = DEF_simulation_time%end_day
   e_seconds = DEF_simulation_time%end_sec
   p_year    = DEF_simulation_time%spinup_year
   p_month   = DEF_simulation_time%spinup_month
   p_day     = DEF_simulation_time%spinup_day
   p_seconds = DEF_simulation_time%spinup_sec

   spinup_repeat = DEF_simulation_time%spinup_repeat

   call initimetype(greenwich)
   call monthday2julian(s_year,s_month,s_day,s_julian)
   call monthday2julian(e_year,e_month,e_day,e_julian)
   call monthday2julian(p_year,p_month,p_day,p_julian)

   sdate(1) = s_year; sdate(2) = s_julian; sdate(3) = s_seconds
   edate(1) = e_year; edate(2) = e_julian; edate(3) = e_seconds
   pdate(1) = p_year; pdate(2) = p_julian; pdate(3) = p_seconds

   CALL Init_GlovalVars
   CAll Init_LC_Const
   CAll Init_PFT_Const

   call pixel%load_from_file    (dir_landdata)
   call gblock%load_from_file   (dir_landdata)

   call mesh_load_from_file (DEF_LC_YEAR, dir_landdata)

   call pixelset_load_from_file (DEF_LC_YEAR, dir_landdata, 'landelm', landelm, numelm)

#ifdef CATCHMENT
   CALL pixelset_load_from_file (DEF_LC_YEAR, dir_landdata, 'landhru', landhru, numhru)
#endif

   call pixelset_load_from_file (DEF_LC_YEAR, dir_landdata, 'landpatch', landpatch, numpatch)

#ifdef PFT_CLASSIFICATION
   call pixelset_load_from_file (DEF_LC_YEAR, dir_landdata, 'landpft', landpft, numpft)
   CALL map_patch_to_pft
#endif

#ifdef PC_CLASSIFICATION
   call pixelset_load_from_file (DEF_LC_YEAR, dir_landdata, 'landpc', landpc, numpc)
   CALL map_patch_to_pc
#endif

#ifdef URBAN_MODEL
   call pixelset_load_from_file (DEF_LC_YEAR, dir_landdata, 'landurban', landurban, numurban)
   CALL map_patch_to_urban
#endif

#if (defined UNSTRUCTURED || defined CATCHMENT)
   CALL elm_vector_init ()
#ifdef CATCHMENT
   CALL hru_vector_init ()
#endif
#endif

   call adj2end(sdate)
   call adj2end(edate)
   call adj2end(pdate)

   ststamp = sdate
   etstamp = edate
   ptstamp = pdate

   IF (ptstamp <= ststamp) THEN
      spinup_repeat = 0
   ELSE
      spinup_repeat = max(0, spinup_repeat)
   ENDIF

   ! ----------------------------------------------------------------------
   ! Read in the model time invariant constant data
   CALL allocate_TimeInvariants ()
   CALL READ_TimeInvariants (DEF_LC_YEAR, casename, dir_restart)

   ! Read in the model time varying data (model state variables)
   CALL allocate_TimeVariables  ()
   CALL READ_TimeVariables (sdate, casename, dir_restart)

   ! Read in SNICAR optical and aging parameters
   CALL SnowOptics_init( DEF_file_snowoptics ) ! SNICAR optical parameters
   CALL SnowAge_init( DEF_file_snowaging )     ! SNICAR aging   parameters

   !-----------------------
   doalb = .true.
   dolai = .true.
   dosst = .false.

   ! Initialize meteorological forcing data module
   call allocate_1D_Forcing ()
   CALL forcing_init (dir_forcing, deltim, sdate)
   call allocate_2D_Forcing (gforc)

   ! Initialize history data module
   call hist_init (dir_hist, DEF_hist_lon_res, DEF_hist_lat_res)
   call allocate_2D_Fluxes (ghist)
   call allocate_1D_Fluxes ()


#if(defined CaMa_Flood)
   call colm_CaMa_init !zhongwang wei, 20210927: initialize CaMa-Flood
#endif
#ifdef OzoneData
   CALL init_Ozone_data(itstamp,sdate)
#endif
#ifdef Fire
   CALL init_lightning_data (itstamp,sdate)
#endif

#if (defined LATERAL_FLOW)
   CALL lateral_flow_init ()
#endif

   ! ======================================================================
   ! begin time stepping loop
   ! ======================================================================

   idate   = sdate
   itstamp = ststamp

   TIMELOOP : DO while (itstamp < etstamp)

      CALL julian2monthday (idate(1), idate(2), month_p, mday_p)

      if (p_is_master) then
         IF (itstamp < ptstamp) THEN
            write(*, 99) idate(1), month_p, mday_p, idate(3), spinup_repeat
         ELSE
            write(*,100) idate(1), month_p, mday_p, idate(3)
         ENDIF
      end if


      Julian_1day_p = int(calendarday(idate)-1)/1*1 + 1
      Julian_8day_p = int(calendarday(idate)-1)/8*8 + 1

      ! Read in the meteorological forcing
      ! ----------------------------------------------------------------------
      CALL read_forcing (idate, dir_forcing)

#ifdef OzoneData
      CALL update_Ozone_data(itstamp, deltim)
#endif
#ifdef Fire
      CALL update_lightning_data (itstamp, deltim)
#endif

      ! Calendar for NEXT time step
      ! ----------------------------------------------------------------------
      CALL TICKTIME (deltim,idate)
      itstamp = itstamp + int(deltim)

      ! lateral flow
#if (defined LATERAL_FLOW)
      CALL lateral_flow (deltim)
#endif

      ! Call colm driver
      ! ----------------------------------------------------------------------
      IF (p_is_worker) THEN
         CALL CoLMDRIVER (idate,deltim,dolai,doalb,dosst,oroflag)
      ENDIF

      ! Get leaf area index
      ! ----------------------------------------------------------------------
#if(defined DYN_PHENOLOGY)
      ! Update once a day
      dolai = .false.
      Julian_1day = int(calendarday(idate)-1)/1*1 + 1
      if(Julian_1day /= Julian_1day_p)then
         dolai = .true.
      endif
#else
      ! READ in Leaf area index and stem area index
      ! Update every 8 days (time interval of the MODIS LAI data)
      ! ----------------------------------------------------------------------
      !zhongwang wei, 20210927: add option to read non-climatological mean LAI
      IF (DEF_LAI_CLIM) then
         ! yuan, 08/03/2019: read global LAI/SAI data
         CALL julian2monthday (idate(1), idate(2), month, mday)
         IF (month /= month_p) THEN
            IF (DEF_LAICHANGE) THEN
               CALL LAI_readin (idate(1), month, dir_landdata)
#ifdef URBAN_MODEL
               CALL UrbanLAI_readin(idate(1), month, dir_landdata)
#endif
            ELSE
               CALL LAI_readin (DEF_LC_YEAR, month, dir_landdata)
#ifdef URBAN_MODEL
               CALL UrbanLAI_readin(DEF_LC_YEAR, month, dir_landdata)
#endif
            ENDIF
         ENDIF
      ELSE
         Julian_8day = int(calendarday(idate)-1)/8*8 + 1
         if(Julian_8day /= Julian_8day_p)then
            CALL LAI_readin (idate(1), Julian_8day, dir_landdata)
         ENDIF
      ENDIF
#endif

#ifdef BGC
#ifdef NITRIF
      CALL julian2monthday (idate(1), idate(2), month, mday)
      if(mday .eq. 1)then
         CALL NITRIF_readin(month, dir_landdata)
      end if
#endif
      if(idate(2) .eq. 1)then
         isread = .true.
      else
         isread = .false.
      end if
      CALL NDEP_readin(idate(1), dir_landdata, isread, .true.)
#ifdef Fire
      if(idate(2)  .eq. 1 .and. idate(3) .eq. 1800)then
         CALL Fire_readin(idate(1), dir_landdata)
      end if
#endif
#endif

#if(defined CaMa_Flood)
   call colm_CaMa_drv(idate(3)) !zhongwang wei, 20210927: run CaMa-Flood
#endif

      ! Write out the model variables for restart run and the histroy file
      ! ----------------------------------------------------------------------
      call hist_out (idate, deltim, itstamp, ptstamp, dir_hist, casename)

      if (save_to_restart (idate, deltim, itstamp, ptstamp)) then
         call WRITE_TimeVariables (idate, casename, dir_restart)
      endif

#ifdef CoLMDEBUG
      call check_TimeVariables ()
#endif

#ifdef USEMPI
      call mpi_barrier (p_comm_glb, p_err)
#endif

      if (p_is_master) then
         call system_clock (end_time, count_rate = c_per_sec)
         time_used = (end_time - start_time) / c_per_sec
         if (time_used >= 3600) then
            write(*,101) time_used/3600, mod(time_used,3600)/60, mod(time_used,60)
         elseif (time_used >= 60) then
            write(*,102) time_used/60, mod(time_used,60)
         else
            write(*,103) time_used
         end if
      end if

      IF ((spinup_repeat > 1) .and. (ptstamp <= itstamp)) THEN
         spinup_repeat = spinup_repeat - 1
         idate   = sdate
         itstamp = ststamp
         CALL forcing_reset ()
      ENDIF

   END DO TIMELOOP

   call deallocate_TimeInvariants ()
   call deallocate_TimeVariables  ()
   call deallocate_1D_Forcing     ()
   call deallocate_1D_Fluxes      ()

#if (defined LATERAL_FLOW)
   CALL lateral_flow_final ()
#endif

   call hist_final ()

#ifdef SinglePoint
   CALL single_srfdata_final ()
#endif

#ifdef USEMPI
   call mpi_barrier (p_comm_glb, p_err)
#endif

#if(defined CaMa_Flood)
   call colm_cama_exit !zhongwang wei, 20210927: finalize CaMa-Flood
#endif

   if (p_is_master) then
      write(*,'(/,A25)') 'CoLM Execution Completed.'
   end if

   99  format(/, 'TIMELOOP = ', I4.4, '-', I2.2, '-', I2.2, '-', I5.5, ' Spinup (', I3, ' repeat left)')
   100 format(/, 'TIMELOOP = ', I4.4, '-', I2.2, '-', I2.2, '-', I5.5)
   101 format (/, 'Time elapsed : ', I4, ' hours', I3, ' minutes', I3, ' seconds.')
   102 format (/, 'Time elapsed : ', I3, ' minutes', I3, ' seconds.')
   103 format (/, 'Time elapsed : ', I3, ' seconds.')

#ifdef USEMPI
   CALL spmd_exit
#endif

END PROGRAM CoLM
! ----------------------------------------------------------------------
! EOP
