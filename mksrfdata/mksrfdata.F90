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
!              2. Global Land Cover Characteristics (land cover type, plant leaf area index, Forest Height, ...)
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
USE GlobalVars
USE omp_lib

IMPLICIT NONE

      character(LEN=256) :: dir_rawdata
      character(LEN=256) :: dir_model_landdata
      integer :: lon_points  ! number of input data longitudes
      integer :: lat_points  ! number of input data latitudes
      real(r8) :: edgen      ! northern edge of grid (degrees)
      real(r8) :: edgee      ! eastern edge of grid (degrees)
      real(r8) :: edges      ! southern edge of grid (degrees)
      real(r8) :: edgew      ! western edge of grid (degrees)

      real(r8),allocatable :: latn(:)       ! grid cell latitude, northern edge (deg)
      real(r8),allocatable :: lats(:)       ! grid cell latitude, sourthern edge (deg)
      real(r8),allocatable :: lonw(:)       ! grid cell longitude, western edge (deg)
      real(r8),allocatable :: lone(:)       ! grid cell longitude, eastern edge (deg)
      real(r8),allocatable :: sinn(:)       ! grid cell latitude, northern edge(sin)  
      real(r8),allocatable :: sins(:)       ! grid cell latitude, northern edge(sin)
      real(r8),allocatable :: lonw_rad(:)   ! grid cell longitude, western edge (radian)
      real(r8),allocatable :: lone_rad(:)   ! grid cell longitude, eastern edge (radian)
      real(r8) :: sinn_i(nlat)              ! fine grid cell latitude, northern edge(sin)
      real(r8) :: sins_i(nlat)              ! fine grid cell latitude, northern edge(sin)
      real(r8) :: lonw_rad_i(nlon)          ! fine grid cell longitude, western edge (radian)
      real(r8) :: lone_rad_i(nlon)          ! fine grid cell longitude, eastern edge (radian)
      integer,allocatable :: READ_row_UB(:) ! north boundary index for fine gird cell  
      integer,allocatable :: READ_col_UB(:) ! west boundary index for fine gird cell  
      integer,allocatable :: READ_row_LB(:) ! south boundary index for fine gird cell
      integer,allocatable :: READ_col_LB(:) ! east boundary index for fine gird cell

      integer :: nrow_start
      integer :: nrow_end
      integer :: ncol_start
      integer :: ncol_end
      integer :: nx_fine_gridcell
      integer :: ny_fine_gridcell

      real(r8), allocatable :: area_fine_gridcell(:,:) ! rwadata fine cell area (km**2)

      namelist /mksrfexp/ dir_rawdata,dir_model_landdata, &
                          lon_points,lat_points,edgen,edgee,edges,edgew

      read(5,mksrfexp)
      
      allocate(latn(lat_points))
      allocate(lats(lat_points))
      allocate(lonw(lon_points))
      allocate(lone(lon_points))
      allocate(sinn(lat_points))
      allocate(sins(lat_points))
      allocate(lonw_rad(lon_points))
      allocate(lone_rad(lon_points)) 
      allocate(READ_row_UB(lat_points))
      allocate(READ_row_LB(lat_points))
      allocate(READ_col_UB(lon_points))
      allocate(READ_col_LB(lon_points))
      allocate(area_fine_gridcell(nlon,nlat))

! ...........................................................................
! 1. Read in or create the modeling grids coordinates and related information
! ...........................................................................

#if(defined USER_GRID)
      CALL rdgrid (dir_model_landdata,edgen,edgee,edges,edgew,&
                   lon_points,lat_points,latn,lats,lonw,lone)
#else
      CALL crgrid (dir_model_landdata,edgen,edgee,edges,edgew,&
                   lon_points,lat_points,latn,lats,lonw,lone)
#endif

! ...................................................................
! 2. Read in and estimate the land characteristic data at 30 arc second resolution
! ...................................................................

#if(defined RAWdata_update)

! ONLY update for USGS landcover
#ifdef USGS_CLASSIFICATION
      CALL rd_land_types ( dir_rawdata )
#endif

      CALL rd_soil_properties ( dir_rawdata )
#endif

! ................................................................
! 3. Mapping land characteristic parameters to the model grids
! ................................................................

      CALL info_gridcell ( lon_points,lat_points,edgen,edgee,edges,edgew, &
                           nrow_start,nrow_end,ncol_start,ncol_end, &
                           nx_fine_gridcell,ny_fine_gridcell,area_fine_gridcell,latn,lats,lonw,lone,&
                           sinn,sins,lonw_rad,lone_rad,sinn_i,sins_i,lonw_rad_i,lone_rad_i,&
                           READ_row_UB,READ_row_LB,READ_col_UB,READ_col_LB )

! read sub-grid structure data
#ifdef USGS_CLASSIFICATION
      CALL aggregation_landtypes ( dir_rawdata,dir_model_landdata, &
                           lon_points,lat_points, &
                           nrow_start,nrow_end,ncol_start,ncol_end, &
                           nx_fine_gridcell,ny_fine_gridcell,area_fine_gridcell,&
                           sinn,sins,lonw_rad,lone_rad,sinn_i,sins_i,lonw_rad_i,lone_rad_i,&
                           READ_row_UB,READ_row_LB,READ_col_UB,READ_col_LB )

      CALL aggregation_LAI ( dir_rawdata,dir_model_landdata, &
                           lon_points,lat_points, &
                           nrow_start,nrow_end,ncol_start,ncol_end, &
                           nx_fine_gridcell,ny_fine_gridcell,area_fine_gridcell,&
                           sinn,sins,lonw_rad,lone_rad,sinn_i,sins_i,lonw_rad_i,lone_rad_i,&
                           READ_row_UB,READ_row_LB,READ_col_UB,READ_col_LB )

      CALL aggregation_forest_height ( dir_rawdata,dir_model_landdata, &
                           lon_points,lat_points, &
                           nrow_start,nrow_end,ncol_start,ncol_end, &
                           nx_fine_gridcell,ny_fine_gridcell,area_fine_gridcell,&
                           READ_row_UB,READ_row_LB,READ_col_UB,READ_col_LB )
#endif


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! yuan, 07/30/2019: landtype, LAI, forest_height read from NC files in mkinidata
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


      CALL aggregation_soil_parameters ( dir_rawdata,dir_model_landdata, &
                           lon_points,lat_points, &
                           nrow_start,nrow_end,ncol_start,ncol_end, &
                           nx_fine_gridcell,ny_fine_gridcell,area_fine_gridcell,&
                           READ_row_UB,READ_row_LB,READ_col_UB,READ_col_LB )

      CALL aggregation_soil_brightness ( dir_rawdata,dir_model_landdata, &
                           lon_points,lat_points, &
                           nrow_start,nrow_end,ncol_start,ncol_end, &
                           nx_fine_gridcell,ny_fine_gridcell,area_fine_gridcell,&
                           READ_row_UB,READ_row_LB,READ_col_UB,READ_col_LB )

      CALL aggregation_lakedepth( dir_rawdata,dir_model_landdata, &
                           lon_points,lat_points, &
                           nrow_start,nrow_end,ncol_start,ncol_end, &
                           nx_fine_gridcell,ny_fine_gridcell,area_fine_gridcell,&
                           READ_row_UB,READ_row_LB,READ_col_UB,READ_col_LB )

1000  PRINT*, 'Successful in surface data making'

      deallocate(latn)
      deallocate(lats)
      deallocate(lonw)
      deallocate(lone)
      deallocate(sinn)
      deallocate(sins)
      deallocate(lonw_rad)
      deallocate(lone_rad) 
      deallocate(READ_row_UB)
      deallocate(READ_row_LB)
      deallocate(READ_col_UB)
      deallocate(READ_col_LB)
      deallocate(area_fine_gridcell)

END PROGRAM mksrfdata
! ----------------------------------------------------------------------
! EOP
