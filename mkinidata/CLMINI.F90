#include <define.h>

PROGRAM CLMINI
   ! ======================================================================
   ! Initialization of Land Characteristic Parameters and Initial State Variables
   !
   ! Reference:
   !     [1] Dai et al., 2003: The Common Land Model (CoLM).
   !         Bull. of Amer. Meter. Soc., 84: 1013-1023
   !     [2] Dai et al., 2004: A two-big-leaf model for canopy temperature,
   !         photosynthesis and stomatal conductance. Journal of Climate
   !     [3] Dai et al., 2014: The Terrestrial Modeling System (TMS).
   !
   !     Created by Yongjiu Dai Februay 2004
   !     Revised by Yongjiu Dai Februay 2014
   ! ======================================================================

   use precision
   use mod_namelist
   use spmd_task
   use mod_block
   use mod_pixel
   use mod_mesh
   USE mod_landelm
#ifdef CATCHMENT
   USE mod_landhru
#endif
   use mod_landpatch
   use mod_srfdata_restart
   USE GlobalVars
   USE LC_Const
   USE PFT_Const
   USE timemanager
#ifdef PFT_CLASSIFICATION
   USE mod_landpft
#endif
#ifdef PC_CLASSIFICATION
   USE mod_landpc
#endif
#ifdef URBAN_MODEL
   USE mod_landurban
#endif
#ifdef SinglePoint
   USE mod_single_srfdata
#endif
#if (defined UNSTRUCTURED || defined CATCHMENT) 
   USE mod_elm_vector
#endif
#ifdef CATCHMENT
   USE mod_hru_vector
#endif
   implicit none

   ! ----------------local variables ---------------------------------
   character(len=256) :: nlfile
   character(LEN=256) :: casename ! case name
   character(LEN=256) :: dir_landdata
   character(LEN=256) :: dir_restart
   CHARACTER(LEN=256) :: fsrfdata
   integer  :: s_year      ! starting date for run in year
   integer  :: s_month     ! starting date for run in month
   integer  :: s_day       ! starting date for run in day
   integer  :: s_julian    ! starting date for run in julian day
   integer  :: s_seconds   ! starting time of day for run in seconds
   integer  :: idate(3)    ! starting date
   logical  :: greenwich   ! true: greenwich time, false: local time

   integer*8 :: start_time, end_time, c_per_sec, time_used

#ifdef USEMPI
   call spmd_init ()
#endif
      
   if (p_is_master) then
      call system_clock (start_time)
   end if

   ! ----------------------------------------------------------------------
   call getarg (1, nlfile)
   call read_namelist (nlfile)

   casename     = DEF_CASE_NAME        
   dir_landdata = DEF_dir_landdata 
   dir_restart  = DEF_dir_restart  
   greenwich    = DEF_simulation_time%greenwich    
   s_year       = DEF_simulation_time%start_year 
   s_month      = DEF_simulation_time%start_month
   s_day        = DEF_simulation_time%start_day
   s_seconds    = DEF_simulation_time%start_sec

#ifdef SinglePoint
   fsrfdata = trim(dir_landdata) // '/srfdata.nc'
   CALL read_surface_data_single (fsrfdata, mksrfdata=.false.)
#endif

   CALL monthday2julian(s_year,s_month,s_day,s_julian)
   idate(1) = s_year; idate(2) = s_julian; idate(3) = s_seconds
   CALL adj2end(idate)

   CALL Init_GlovalVars
   CAll Init_LC_Const
   CAll Init_PFT_Const

   call pixel%load_from_file  (dir_landdata)
   call gblock%load_from_file (dir_landdata)
   call mesh_load_from_file   (s_year, dir_landdata)

   CALL pixelset_load_from_file (DEF_LC_YEAR, dir_landdata, 'landelm', landelm, numelm)

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
   CALL pixelset_load_from_file (DEF_LC_YEAR, dir_landdata, 'landurban', landurban, numurban)
   CALL map_patch_to_urban
#endif
#if (defined UNSTRUCTURED || defined CATCHMENT) 
   CALL elm_vector_init ()
#ifdef CATCHMENT
   CALL hru_vector_init ()
#endif
#endif

   CALL initialize (casename, dir_landdata, dir_restart, idate, greenwich)

#ifdef SinglePoint
   CALL single_srfdata_final ()
#endif

#ifdef USEMPI
   call mpi_barrier (p_comm_glb, p_err)
#endif

   if (p_is_master) then
      call system_clock (end_time, count_rate = c_per_sec)
      time_used = (end_time - start_time) / c_per_sec
      if (time_used >= 3600) then
         write(*,101) time_used/3600, mod(time_used,3600)/60, mod(time_used,60)
         101 format (/,'Overall system time used:', I4, ' hours', I3, ' minutes', I3, ' seconds.')
      elseif (time_used >= 60) then
         write(*,102) time_used/60, mod(time_used,60)
         102 format (/,'Overall system time used:', I3, ' minutes', I3, ' seconds.')
      else 
         write(*,103) time_used
         103 format (/,'Overall system time used:', I3, ' seconds.')
      end if

      write(*,*) 'CLM Initialization Execution Completed'
   end if
   
#ifdef USEMPI
   call spmd_exit
#endif
   
END PROGRAM CLMINI
! ----------------------------------------------------------------------
! EOP
