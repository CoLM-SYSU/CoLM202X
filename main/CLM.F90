#include <define.h>

PROGRAM CLM
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
   USE GlobalVars
   USE LC_Const
   USE PFT_Const
   use PhysicalConstants
   use MOD_TimeInvariants
   use MOD_TimeVariables
   use MOD_1D_Forcing
   use MOD_2D_Forcing
   use MOD_1D_Fluxes
   use MOD_2D_Fluxes
   use mod_forcing
   use mod_hist
   use timemanager
   use mod_colm_debug

   use mod_block
   use mod_pixel
   use mod_hydrounit
   use mod_landpatch
   use mod_srfdata_restart
#ifdef PFT_CLASSIFICATION
   USE mod_landpft
#endif
#ifdef PC_CLASSIFICATION
   USE mod_landpc
#endif
#if(defined CaMa_Flood)
   use colm_CaMaMod
#endif
#ifdef vsf_statistics
   USE mod_soil_water, only : count_iters
#endif
#ifdef SinglePoint
   USE mod_single_srfdata
#endif

   IMPLICIT NONE

   ! ----------------local variables ---------------------------------
   character(LEN=256) :: nlfile
   character(LEN=256) :: casename   
   character(len=256) :: dir_landdata
   character(len=256) :: dir_forcing
   character(len=256) :: dir_hist
   character(len=256) :: dir_restart

   real(r8) :: deltim           ! time step (senconds)
   integer  :: idate(3)         ! calendar (year, julian day, seconds)
   integer  :: edate(3)         ! calendar (year, julian day, seconds)
   integer  :: pdate(3)         ! calendar (year, julian day, seconds)
   logical  :: greenwich        ! greenwich time
   
   logical :: doalb            ! true => start up the surface albedo calculation
   logical :: dolai            ! true => start up the time-varying vegetation paramter
   logical :: dosst            ! true => update sst/ice/snow

   integer :: Julian_1day_p, Julian_1day 
   integer :: Julian_8day_p, Julian_8day 
   integer :: s_year, s_month, s_day, s_seconds, s_julian
   integer :: e_year, e_month, e_day, e_seconds, e_julian
   integer :: p_year, p_month, p_day, p_seconds, p_julian
   INTEGER :: month, mday, month_p, mday_p

   type(timestamp) :: itstamp, etstamp, ptstamp
   
   integer*8 :: start_time, end_time, c_per_sec, time_used

#ifdef USEMPI
   call spmd_init ()
#endif

#ifdef vsf_statistics
   count_iters(:) = 0
#endif

   if (p_is_master) then
      call system_clock (start_time)
   end if

   call getarg (1, nlfile)

   call read_namelist (nlfile)

#ifdef SinglePoint
   CALL read_surface_data_single (SITE_fsrfdata)
#endif

   casename     = DEF_CASE_NAME
   dir_landdata = DEF_dir_landdata
   dir_forcing  = DEF_dir_forcing
   dir_hist     = DEF_dir_history
   dir_restart  = DEF_dir_restart

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

   call initimetype(greenwich)
   call monthday2julian(s_year,s_month,s_day,s_julian)
   call monthday2julian(e_year,e_month,e_day,e_julian)
   call monthday2julian(p_year,p_month,p_day,p_julian)

   idate(1) = s_year; idate(2) = s_julian; idate(3) = s_seconds
   edate(1) = e_year; edate(2) = e_julian; edate(3) = e_seconds
   pdate(1) = p_year; pdate(2) = p_julian; pdate(3) = p_seconds
   
   CALL Init_GlovalVars
   CAll Init_LC_Const
   CAll Init_PFT_Const

   call pixel%load_from_file    (dir_landdata)
   call gblock%load_from_file   (dir_landdata)

   call landbasin_load_from_file (dir_landdata)

   call pixelset_load_from_file (dir_landdata, 'hydrounit', hydrounit, numhru)
   
   call pixelset_load_from_file (dir_landdata, 'landpatch', landpatch, numpatch)

#ifdef PFT_CLASSIFICATION
   call pixelset_load_from_file (dir_landdata, 'landpft', landpft, numpft)
   CALL map_patch_to_pft
#endif

#ifdef PC_CLASSIFICATION
   call pixelset_load_from_file (dir_landdata, 'landpc', landpc, numpc)
   CALL map_patch_to_pc
#endif

   call adj2end(edate)
   call adj2end(pdate)

   itstamp = idate
   etstamp = edate
   ptstamp = pdate

   ! ----------------------------------------------------------------------
   ! Read in the model time invariant constant data
   CALL allocate_TimeInvariants ()
   CALL READ_TimeInvariants (casename, dir_restart)

   ! Read in the model time varying data (model state variables)
   CALL allocate_TimeVariables  ()
   CALL READ_TimeVariables (idate, casename, dir_restart)

   !-----------------------

   doalb = .true.
   dolai = .true.
   dosst = .false.

   ! Initialize meteorological forcing data module
   CALL forcing_init (dir_forcing, deltim, idate)
   call allocate_2D_Forcing (gforc)
   call allocate_1D_Forcing ()

   ! Initialize history data module
   print*, dir_hist, DEF_hist_lon_res, DEF_hist_lat_res
   call hist_init (dir_hist, DEF_hist_lon_res, DEF_hist_lat_res)
   call allocate_2D_Fluxes (ghist)
   call allocate_1D_Fluxes ()
#if(defined CaMa_Flood)
    call colm_CaMa_init
#endif

   ! ======================================================================
   ! begin time stepping loop
   ! ======================================================================

   TIMELOOP : DO while (itstamp < etstamp)
      
      CALL julian2monthday (idate(1), idate(2), month_p, mday_p)

      if (p_is_master) then
         write(*,100) idate(1), month_p, mday_p, idate(3)
         100 format(/, 'TIMELOOP = ', I4.4, '-', I2.2, '-', I2.2, '-', I5.5)
      end if

      Julian_1day_p = int(calendarday(idate)-1)/1*1 + 1
      Julian_8day_p = int(calendarday(idate)-1)/8*8 + 1

      ! Read in the meteorological forcing
      ! ----------------------------------------------------------------------
      CALL read_forcing (idate, dir_forcing)


      ! Calendar for NEXT time step
      ! ----------------------------------------------------------------------
      CALL TICKTIME (deltim,idate)
      itstamp = itstamp + int(deltim)

      ! Call clm driver
      ! ----------------------------------------------------------------------
      IF (p_is_worker) THEN
         CALL CLMDRIVER (idate,deltim,dolai,doalb,dosst,oroflag)
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
            CALL LAI_readin (idate(1), month, dir_landdata)
         END IF
      ELSE
         Julian_8day = int(calendarday(idate)-1)/8*8 + 1
         if(Julian_8day /= Julian_8day_p)then
            CALL LAI_readin (idate(1), Julian_8day, dir_landdata)
         ENDIF
      ENDIF
#endif

#ifdef NITRIF
      CALL julian2monthday (idate(1), idate(2), month, mday)
      if(mday .eq. 1)then
         CALL NITRIF_readin(month, dir_landdata)
      end if
#endif

!!!! need to acc runoff here!!!
#if(defined CaMa_Flood)
call colm_CaMa_drv
#endif

   
      ! Write out the model variables for restart run and the histroy file
      ! ----------------------------------------------------------------------
      call hist_out (idate, deltim, itstamp, ptstamp, dir_hist, casename)

      if (save_to_restart (idate, deltim, itstamp, ptstamp)) then
         call WRITE_TimeVariables (idate, casename, dir_restart)
      end if

#ifdef CLMDEBUG
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

   END DO TIMELOOP

   call deallocate_TimeInvariants ()  
   call deallocate_TimeVariables  ()
   call deallocate_1D_Forcing     ()
   call deallocate_1D_Fluxes      ()

   call hist_final ()


#ifdef USEMPI
   call mpi_barrier (p_comm_glb, p_err)
#endif

#if(defined CaMa_Flood) 
   call colm_cama_exit
#endif

   if (p_is_master) then
      write(*,'(/,A25)') 'CLM Execution Completed.'
   end if

#ifdef vsf_statistics
   IF (p_is_worker) THEN
#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, count_iters, size(count_iters), &
         MPI_INTEGER, MPI_SUM, p_comm_worker, p_err)
#endif
      IF (p_iam_worker == 0) THEN
         write(*,*) 'iteration stat ', count_iters(:)
      ENDIF
   ENDIF
#endif
         
   101 format (/, 'Time elapsed : ', I4, ' hours', I3, ' minutes', I3, ' seconds.') 
   102 format (/, 'Time elapsed : ', I3, ' minutes', I3, ' seconds.')
   103 format (/, 'Time elapsed : ', I3, ' seconds.')

#ifdef USEMPI
   CALL spmd_exit
#endif

END PROGRAM CLM
! ----------------------------------------------------------------------
! EOP
