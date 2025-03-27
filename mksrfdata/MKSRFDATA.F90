#include <define.h>

PROGRAM MKSRFDATA

!=======================================================================
!  Surface grid edges:
!  The model domain was defined with the north, east, south, west edges:
!           edgen: northern edge of grid : > -90 and <= 90 (degrees)
!           edgee: eastern edge of grid  : > western edge and <= 180
!           edges: southern edge of grid : >= -90  and <  90
!           edgew: western edge of grid  : >= -180 and < 180
!
!  Region (global) latitude grid goes from:
!                  NORTHERN edge (POLE) to SOUTHERN edge (POLE)
!  Region (global) longitude grid starts at:
!                  WESTERN edge (DATELINE with western edge)
!                  West of Greenwich defined negative for global grids,
!                  the western edge of the longitude grid starts at the dateline
!
!  Land characteristics at the 30 arc-seconds grid resolution (RAW DATA):
!               1. Global Terrain Dataset (elevation height, topography-based
!                  factors)
!               2. Global Land Cover Characteristics (land cover type, plant
!                  leaf area index, Forest Height, ...)
!               3. Global Lakes and Wetlands Characteristics (lake and wetlands
!                  types, lake coverage and lake depth)
!               4. Global Glacier Characteristics
!               5. Global Urban Characteristics (urban extent, ...)
!               6. Global Soil Characteristics (...)
!               7. Global Cultural Characteristics (ON-GONG PROJECT)
!
!  Land characteristics at the model grid resolution (CREATED):
!               1. Model grid (longitude, latitude)
!               2. Fraction (area) of patches of grid (0-1)
!                  2.1 Fraction of land water bodies (lake, reservoir, river)
!                  2.2 Fraction of wetland
!                  2.3 Fraction of glacier
!                  2.4 Fraction of urban and built-up
!                  ......
!               3. Plant leaf area index
!               4. Tree height
!               5. Lake depth
!               6. Soil thermal and hydraulic parameters
!
!  Created by Yongjiu Dai, 02/2014
!
! !REVISIONS:
!  Shupeng Zhang, 01/2022: porting codes to MPI parallel version
!
!=======================================================================

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Block
   USE MOD_Pixel
   USE MOD_Grid
   USE MOD_Mesh
   USE MOD_MeshFilter
   USE MOD_TimeManager
   USE MOD_LandElm
#ifdef CATCHMENT
   USE MOD_LandHRU
#endif
   USE MOD_LandPatch
   USE MOD_SrfdataRestart
   USE MOD_Const_LC
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   USE MOD_LandPFT
#endif
#ifdef URBAN_MODEL
   USE MOD_LandUrban
#endif
#ifdef CROP
   USE MOD_LandCrop
#endif
   USE MOD_RegionClip
#ifdef SrfdataDiag
   USE MOD_SrfdataDiag, only: gdiag, srfdata_diag_init
#endif
#ifdef SinglePoint
   USE MOD_SingleSrfdata
#endif

   USE MOD_RegionClip


   IMPLICIT NONE

   character(len=256) :: nlfile
   character(len=256) :: lndname
   character(len=256) :: dir_rawdata
   character(len=256) :: dir_landdata
   real(r8) :: edgen  ! northern edge of grid (degrees)
   real(r8) :: edgee  ! eastern edge of grid (degrees)
   real(r8) :: edges  ! southern edge of grid (degrees)
   real(r8) :: edgew  ! western edge of grid (degrees)

   type (grid_type) :: grid_500m, grid_htop, grid_soil, grid_lai, grid_topo, grid_topo_factor
   type (grid_type) :: grid_urban_5km, grid_urban_500m

   integer   :: lc_year
   character(len=4) :: cyear
   integer*8 :: start_time, end_time, c_per_sec, time_used


#ifdef USEMPI
   CALL spmd_init ()
#endif

   IF (p_is_master) THEN
      CALL system_clock (start_time)
   ENDIF

   CALL getarg(1, nlfile)

   CALL read_namelist (nlfile)

   CALL initimetype (DEF_simulation_time%greenwich)

#ifdef SinglePoint
#ifndef URBAN_MODEL

   CALL read_surface_data_single (SITE_fsitedata, mksrfdata = .true.)
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   CALL write_surface_data_single (numpatch, numpft)
#else
   CALL write_surface_data_single (numpatch)
#endif

#else

   CALL read_urban_surface_data_single (SITE_fsitedata, mksrfdata=.true.)
   CALL write_urban_surface_data_single(numurban)

#endif

   CALL single_srfdata_final ()
   write(*,*)  'Successful in surface data making.'
   CALL CoLM_stop()
#endif

   IF (USE_srfdata_from_larger_region) THEN

      CALL srfdata_region_clip (DEF_dir_existing_srfdata, DEF_dir_landdata)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
      CALL spmd_exit
#endif
      CALL EXIT()
   ENDIF

   IF (USE_srfdata_from_3D_gridded_data) THEN

      ! TODO
      ! CALL srfdata_retrieve_from_3D_data (DEF_dir_existing_srfdata, DEF_dir_landdata)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
      CALL spmd_exit
#endif
      CALL EXIT()
   ENDIF

   dir_rawdata  = DEF_dir_rawdata
   dir_landdata = DEF_dir_landdata
   edges = DEF_domain%edges
   edgen = DEF_domain%edgen
   edgew = DEF_domain%edgew
   edgee = DEF_domain%edgee
   lc_year = DEF_LC_YEAR

   ! define blocks
   CALL gblock%set ()

   CALL Init_GlobalVars
   CAll Init_LC_Const

! ...........................................................................
! 1. Read in or create the modeling grids coordinates and related information
! ...........................................................................

   ! define domain in pixel coordinate
   CALL pixel%set_edges (edges, edgen, edgew, edgee)
   CALL pixel%assimilate_gblock ()

   ! define grid coordinates of mesh
#ifdef GRIDBASED
   CALL init_gridbased_mesh_grid ()
#endif

#ifdef CATCHMENT
   CALL gridmesh%define_by_name ('merit_90m')
#endif

#ifdef UNSTRUCTURED
   CALL gridmesh%define_from_file (DEF_file_mesh)
#endif

   ! define grid coordinates of mesh filter
   has_mesh_filter = inquire_mesh_filter ()
   IF (has_mesh_filter) THEN
      CALL grid_filter%define_from_file (DEF_file_mesh_filter)
   ENDIF

   ! define grid coordinates of hydro units in catchment
#ifdef CATCHMENT
   CALL grid_hru%define_by_name ('merit_90m')
#endif

   CALL grid_500m%define_by_name ('colm_500m')

   ! define grid coordinates of land types
#ifdef LULC_USGS
   CALL grid_patch%define_by_name ('colm_1km')
#endif
#ifdef LULC_IGBP
   CALL grid_patch%define_by_name ('colm_500m')
#endif
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   CALL grid_patch%define_by_name ('colm_500m')
#endif
#if (defined CROP)
   ! define grid for crop parameters
   CALL grid_crop%define_from_file (trim(DEF_dir_rawdata)//'/global_CFT_surface_data.nc', 'lat', 'lon')
#endif

   ! define grid for forest height
#ifdef LULC_USGS
   CALL grid_htop%define_by_name ('colm_1km')
#else
   CALL grid_htop%define_by_name ('colm_500m')
#endif

   ! define grid for soil parameters raw data
   CALL grid_soil%define_by_name ('colm_500m')

   ! define grid for LAI raw data
   CALL grid_lai%define_by_name ('colm_500m')

   ! define grid for topography
   CALL grid_topo%define_by_name ('colm_500m')

   ! define grid for topography factors
   IF (DEF_USE_Forcing_Downscaling) THEN
      lndname = trim(DEF_DS_HiresTopographyDataDir) // '/slope.nc'
      CALL grid_topo_factor%define_from_file (lndname,"lat","lon")
   ENDIF

   ! add by dong, only test for making urban data
#ifdef URBAN_MODEL
   CALL grid_urban%define_by_name      ('colm_500m')
   CALL grid_urban_500m%define_by_name ('colm_500m')
   CALL grid_urban_5km%define_by_name  ('colm_5km' )
#endif

   ! assimilate grids to build pixels
#ifndef SinglePoint
   CALL pixel%assimilate_grid (gridmesh)
#endif
   IF (has_mesh_filter) THEN
      CALL pixel%assimilate_grid (grid_filter)
   ENDIF
#ifdef CATCHMENT
   CALL pixel%assimilate_grid (grid_hru  )
#endif
   CALL pixel%assimilate_grid (grid_500m )
   CALL pixel%assimilate_grid (grid_patch)
#if (defined CROP)
   CALL pixel%assimilate_grid (grid_crop )
#endif
   CALL pixel%assimilate_grid (grid_htop )
   CALL pixel%assimilate_grid (grid_soil )
   CALL pixel%assimilate_grid (grid_lai  )
   CALL pixel%assimilate_grid (grid_topo )

   IF (DEF_USE_Forcing_Downscaling) THEN
      CALL pixel%assimilate_grid (grid_topo_factor)
   ENDIF

#ifdef URBAN_MODEL
   CALL pixel%assimilate_grid (grid_urban     )
   CALL pixel%assimilate_grid (grid_urban_500m)
   CALL pixel%assimilate_grid (grid_urban_5km )
#endif

   ! map pixels to grid coordinates
#ifndef SinglePoint
   CALL pixel%map_to_grid (gridmesh)
#endif
   IF (has_mesh_filter) THEN
      CALL pixel%map_to_grid (grid_filter)
   ENDIF
#ifdef CATCHMENT
   CALL pixel%map_to_grid (grid_hru  )
#endif
   CALL pixel%map_to_grid (grid_500m )
   CALL pixel%map_to_grid (grid_patch)
#if (defined CROP)
   CALL pixel%map_to_grid (grid_crop )
#endif
   CALL pixel%map_to_grid (grid_htop )
   CALL pixel%map_to_grid (grid_soil )
   CALL pixel%map_to_grid (grid_lai  )
   CALL pixel%map_to_grid (grid_topo )

   IF (DEF_USE_Forcing_Downscaling) THEN
      CALL pixel%map_to_grid (grid_topo_factor)
   ENDIF

#ifdef URBAN_MODEL
   CALL pixel%map_to_grid (grid_urban     )
   CALL pixel%map_to_grid (grid_urban_500m)
   CALL pixel%map_to_grid (grid_urban_5km )
#endif


   ! build land elms
   CALL mesh_build ()
   CALL landelm_build

#if (defined GRIDBASED || defined UNSTRUCTURED)
   IF (DEF_LANDONLY) THEN
      !TODO: distinguish USGS and IGBP land cover
#ifndef LULC_USGS
      write(cyear,'(i4.4)') lc_year
      lndname = trim(DEF_dir_rawdata)//'/landtypes/landtype-igbp-modis-'//trim(cyear)//'.nc'
#else
      lndname = trim(DEF_dir_rawdata)//'/landtypes/landtype-usgs-update.nc'
#endif
      CALL mesh_filter (grid_patch, lndname, 'landtype')
   ENDIF
#endif

   ! Filtering pixels
   IF (has_mesh_filter) THEN
      CALL mesh_filter (grid_filter, DEF_file_mesh_filter, 'mesh_filter')
   ENDIF

#ifdef CATCHMENT
   CALL landhru_build
#endif

   ! build land patches
   CALL landpatch_build(lc_year)

#ifdef URBAN_MODEL
   CALL landurban_build(lc_year)
#endif

#ifdef CROP
   CALL landcrop_build (lc_year)
#endif

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   CALL landpft_build  (lc_year)
#endif

! ................................................................
! 2. SAVE land surface tessellation information
! ................................................................

   CALL gblock%save_to_file    (dir_landdata)

   CALL pixel%save_to_file     (dir_landdata)

   CALL mesh_save_to_file      (dir_landdata, lc_year)

   CALL pixelset_save_to_file  (dir_landdata, 'landelm'  , landelm  , lc_year)

#ifdef CATCHMENT
   CALL pixelset_save_to_file  (dir_landdata, 'landhru'  , landhru  , lc_year)
#endif

   CALL pixelset_save_to_file  (dir_landdata, 'landpatch', landpatch, lc_year)

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   CALL pixelset_save_to_file  (dir_landdata, 'landpft'  , landpft  , lc_year)
#endif

#ifdef URBAN_MODEL
   CALL pixelset_save_to_file  (dir_landdata, 'landurban', landurban, lc_year)
#endif

! ................................................................
! 3. Mapping land characteristic parameters to the model grids
! ................................................................
#ifdef SrfdataDiag
   CALL elm_patch%build (landelm, landpatch, use_frac = .true.)
#ifdef GRIDBASED
   CALL gdiag%define_by_copy (gridmesh)
#else
   CALL gdiag%define_by_ndims(3600,1800)
#endif

   CALL srfdata_diag_init (dir_landdata)
#endif

   !TODO: for lulcc, need to run for each year and SAVE to different subdirs

   CALL Aggregation_PercentagesPFT  (grid_500m, dir_rawdata, dir_landdata, lc_year)

   CALL Aggregation_LakeDepth       (grid_500m, dir_rawdata, dir_landdata, lc_year)

   CALL Aggregation_SoilParameters  (grid_soil, dir_rawdata, dir_landdata, lc_year)

   CALL Aggregation_SoilBrightness  (grid_500m, dir_rawdata, dir_landdata, lc_year)

   IF (DEF_USE_BEDROCK) THEN
      CALL Aggregation_DBedrock     (grid_500m, dir_rawdata, dir_landdata)
   ENDIF

   CALL Aggregation_LAI             (grid_lai,  dir_rawdata, dir_landdata, lc_year)

   CALL Aggregation_ForestHeight    (grid_htop, dir_rawdata, dir_landdata, lc_year)

   CALL Aggregation_Topography      (grid_topo, dir_rawdata, dir_landdata, lc_year)

   IF (DEF_USE_Forcing_Downscaling) THEN
      CALL Aggregation_TopographyFactors (grid_topo_factor, &
         trim(DEF_DS_HiresTopographyDataDir), dir_landdata, lc_year)
   ENDIF

#ifdef URBAN_MODEL
   CALL Aggregation_urban (dir_rawdata, dir_landdata, lc_year, &
                           grid_urban_5km, grid_urban_500m)
#endif

   CALL Aggregation_SoilTexture     (grid_soil, dir_rawdata, dir_landdata, lc_year)

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

END PROGRAM MKSRFDATA
! ----------------------------------------------------------------------
! EOP
