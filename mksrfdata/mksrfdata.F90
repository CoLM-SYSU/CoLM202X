#include <define.h>

PROGRAM mksrfdata
   ! ======================================================================
   ! Surface grid edges: 
   ! The model domain was defined with the north, east, south, west edges:
   !          edgen: northern edge of grid : > -90 and <= 90 (degrees)
   !          edgee: eastern edge of grid  : > western edge and <= 180
   !          edges: southern edge of grid : >= -90  and <  90
   !          edgew: western edge of grid  : >= -180 and < 180
   !
   ! Region (global) latitude grid goes from:
   !                 NORTHERN edge (POLE) to SOUTHERN edge (POLE)
   ! Region (global) longitude grid starts at:
   !                 WESTERN edge (DATELINE with western edge)
   !                 West of Greenwich defined negative for global grids, 
   !                 the western edge of the longitude grid starts at the dateline 
   !
   ! Land characteristics at the 30 arc-seconds grid resolution (RAW DATA):
   !              1. Global Terrain Dataset (elevation height,...)
   !              2. Global Land Cover Characteristics (land cover TYPE, plant leaf area index, Forest Height, ...)
   !              3. Global Lakes and Wetlands Characteristics (lake and wetlands types, lake coverage and lake depth)
   !              4. Global Glacier Characteristics
   !              5. Global Urban Characteristics (urban extent, ...)
   !              6. Global Soil Characteristics (...)
   !              7. Global Cultural Characteristics (ON-GONG PROJECT)
   !
   ! Land charateristics at the model grid resolution (CREATED):
   !              1. Model grid (longitude, latitude)
   !              2. Fraction (area) of patches of grid (0-1)
   !                 2.1 Fraction of land water bodies (lake, reservoir, river) 
   !                 2.2 Fraction of wetland
   !                 2.3 Fraction of glacier
   !                 2.4 Fraction of urban and built-up
   !                 ......
   !              3. Plant leaf area index
   !              4. Tree height
   !              5. Lake depth
   !              6. Soil thermal and hydraulic parameters
   !
   ! Created by Yongjiu Dai, 02/2014
   !
   ! ________________
   ! REVISION HISTORY:
   !   /07/2014, Siguang Zhu & Xiangxiang Zhang: modifiy the aggregation_xxx
   !               interfaces for reading a user defined domain file.
   !
   ! ======================================================================
   USE precision
   USE spmd_task
   USE mod_namelist
   USE mod_block
   USE mod_pixel
   USE mod_grid
   USE mod_landunit
   USE mod_landcell
   USE mod_landpatch
   USE LC_Const
   USE mod_srfdata_restart
#ifdef PFT_CLASSIFICATION
   USE mod_landpft
#endif
#ifdef PC_CLASSIFICATION
   USE mod_landpc
#endif

   IMPLICIT NONE

   CHARACTER(len=256) :: nlfile

   CHARACTER(LEN=256) :: dir_rawdata
   CHARACTER(LEN=256) :: dir_landdata
   REAL(r8) :: edgen  ! northern edge of grid (degrees)
   REAL(r8) :: edgee  ! eastern edge of grid (degrees)
   REAL(r8) :: edges  ! southern edge of grid (degrees)
   REAL(r8) :: edgew  ! western edge of grid (degrees)

   TYPE (grid_type) :: gunit

   INTEGER :: start_time, end_time, c_per_sec, time_used


#ifdef USEMPI
   CALL spmd_init ()
#endif

   IF (p_is_master) THEN
      CALL system_clock (start_time)
   ENDIF

   CALL getarg(1, nlfile)

   CALL read_namelist (nlfile)
   
   ! define blocks
   CALL gblock%set_by_size (DEF_nx_blocks, DEF_ny_blocks)

   dir_rawdata  = DEF_dir_rawdata
   dir_landdata = DEF_dir_landdata
   edges = DEF_domain%edges
   edgen = DEF_domain%edgen
   edgew = DEF_domain%edgew
   edgee = DEF_domain%edgee

   CAll Init_LC_Const

   ! ...........................................................................
   ! 1. Read in or create the modeling grids coordinates and related information
   ! ...........................................................................

   ! define domain in pixel coordinate
   CALL pixel%set_edges (edges, edgen, edgew, edgee) 
   CALL pixel%assimilate_gblock ()

   ! define grid coordinates of land units
#ifdef GRIDBASED
   CALL gunit%define_from_file (DEF_file_landgrid)
#endif
#ifdef CATCHMENT
   CALL gunit%define_by_name ('merit_90m')
#endif

   ! define grid coordinates of height bands in catchment
#ifdef CATCHMENT
   CALL ghband%define_by_name ('merit_90m')
#endif

   ! define grid coordinates of land types
#ifdef USGS_CLASSIFICATION
   CALL gpatch%define_by_name ('colm_1km')
#endif
#ifdef IGBP_CLASSIFICATION
   CALL gpatch%define_by_name ('colm_500m')
#endif
#ifdef PFT_CLASSIFICATION
   CALL gpatch%define_by_name ('colm_500m')
#if (defined CROP) 
   CALL gcrop%define_by_ndims (720,360)
#endif
#endif
#ifdef PC_CLASSIFICATION
   CALL gpatch%define_by_name ('colm_500m')
#endif

   ! assimilate grids to build pixels
   CALL pixel%assimilate_grid (gunit)
#ifdef CATCHMENT
   CALL pixel%assimilate_grid (ghband)
#endif
   CALL pixel%assimilate_grid (gpatch)
#if (defined CROP) 
   CALL pixel%assimilate_grid (gcrop )
#endif

   ! map pixels to grid coordinates
   CALL pixel%map_to_grid (gunit)
#ifdef CATCHMENT
   CALL pixel%map_to_grid (ghband)
#endif
   CALL pixel%map_to_grid (gpatch)
#if (defined CROP) 
   CALL pixel%map_to_grid (gcrop )
#endif

   ! build land units 
   CALL landunit_build (gunit)
   
   ! build land cells
   CALL landcell_build 

   ! build land patches
   CALL landpatch_build

#ifdef PFT_CLASSIFICATION
   CALL landpft_build
#endif
   
#ifdef PC_CLASSIFICATION
   CALL landpc_build
#endif

   ! ................................................................
   ! 2. SAVE land surface tessellation information 
   ! ................................................................

   CALL gblock%save_to_file    (dir_landdata)
   
   CALL pixel%save_to_file     (dir_landdata)

   CALL landunit_save_to_file  (dir_landdata)

   CALL pixelset_save_to_file  (dir_landdata, 'landcell',  landcell)
   
   CALL pixelset_save_to_file  (dir_landdata, 'landpatch', landpatch)

#ifdef PFT_CLASSIFICATION
   CALL pixelset_save_to_file  (dir_landdata, 'landpft'  , landpft  )
#endif

#ifdef PC_CLASSIFICATION
   CALL pixelset_save_to_file  (dir_landdata, 'landpc'   , landpc   )
#endif

   ! ................................................................
   ! 3. Mapping land characteristic parameters to the model grids
   ! ................................................................

   CALL aggregation_soil_parameters (gpatch, dir_rawdata, dir_landdata)

   CALL aggregation_soil_brightness (gpatch, dir_rawdata, dir_landdata)

   CALL aggregation_lakedepth       (gpatch, dir_rawdata, dir_landdata)
   
#ifdef USE_DEPTH_TO_BEDROCK
   CALL aggregation_dbedrock        (gpatch, dir_rawdata, dir_landdata)
#endif

   CALL aggregation_percentages     (gpatch, dir_rawdata, dir_landdata)
   
   CALL aggregation_LAI             (gpatch, dir_rawdata, dir_landdata)

   CALL aggregation_forest_height   (gpatch, dir_rawdata, dir_landdata)

   ! ................................................................
   ! 4. Free memories. 
   ! ................................................................

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif
   
   IF (p_is_master) THEN
      CALL system_clock (end_time, count_rate = c_per_sec)
      time_used = (end_time - start_time) / c_per_sec
      IF (time_used >= 3600) THEN
         write(*,101) time_used/3600, mod(time_used,3600)/60, mod(time_used,60)
         101 format (/, 'Overall system time used:', I4, ' hours', I3, ' minutes', I3, ' seconds.')
      ELSEIF (time_used >= 60) THEN
         write(*,102) time_used/60, mod(time_used,60)
         102 format (/, 'Overall system time used:', I3, ' minutes', I3, ' seconds.')
      ELSE 
         write(*,103) time_used
         103 format (/, 'Overall system time used:', I3, ' seconds.')
      ENDIF

      write(*,*)  'Successful in surface data making.'
   ENDIF

#ifdef USEMPI
   CALL spmd_exit          
#endif

END PROGRAM mksrfdata
! ----------------------------------------------------------------------
! EOP
