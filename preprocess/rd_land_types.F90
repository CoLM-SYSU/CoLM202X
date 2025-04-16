#include <define.h>
SUBROUTINE rd_land_types(dir_rawdata)
!-----------------------------------------------------------------------
!  => Read in land cover dataset from original "raw" data files -
!      data with 30 arc seconds resolution
!  => Fill the missing data
!  => Correct and update the land types with the specific datasets
!
!  1. Global Elevation Dataset (GTOPO30)
!     (http://webgis.wr.usgs.gov/globalgis/gtopo30/)
!     The elevation values range from -407 to 8,752 meters.
!     ocean areas have been assigned a value of -9999.
!
!  2. Global Land Cover Characteristics
!  2.1 GLCC USGS Land cover /land use types
!      (http://edc2.usgs.gov/glcc/)
!        Value   Description
!        0       Ocean
!        1       Urban and Built-Up Land
!        2       Dryland Cropland and Pasture
!        3       Irrigated Cropland and Pasture
!        4       Mixed Dryland/Irrigated Cropland and Pasture
!        5       Cropland/Grassland Mosaic
!        6       Cropland/Woodland Mosaic
!        7       Grassland
!        8       Shrubland
!        9       Mixed Shrubland/Grassland
!        10      Savanna
!        11      Deciduous Broadleaf Forest
!        12      Deciduous Needleleaf Forest
!        13      Evergreen Broadleaf Forest
!        14      Evergreen Needleleaf Forest
!        15      Mixed Forest
!        16      Land Water Bodies
!        17      Herbaceous Wetland
!        18      Wooded Wetland
!        19      Barren or Sparsely Vegetated
!        20      Herbaceous Tundra
!        21      Wooded Tundra
!        22      Mixed Tundra
!        23      Bare Ground Tundra
!        24      Snow or Ice
!
!  2.2 MODIS Collection 5 IGBP global land cover
!        Value   Description
!        0       Ocean
!        1       Evergreen Needleleaf Forest
!        2       Evergreen Broadleaf Forest
!        3       Deciduous Needleleaf Forest
!        4       Deciduous Broadleaf Forest
!        5       Mixed Forest
!        6       Closed Shrublands
!        7       Open Shrublands
!        8       Woody Savannas
!        9       Savannas
!        10      Grasslands
!        11      Permanent Wetlands
!        12      Croplands
!        13      Urban and Built-Up
!        14      Cropland/Natural Vegetation Mosaic
!        15      Snow and Ice
!        16      Barren or Sparsely Vegetated
!        17      Land Water Bodies
!
!  2.3 PFT classification (ON GOING PROJECT ......)
!
!  3. Global Glacier/Ice Sheet Characteristics
!     (http://www.glims.org/RGI/; http://glims.colorado.edu/glacierdata/)
!
!  4. Global Lake and Wetlands Types
!     (http://www.wwfus.org/science/data.cfm)
!        1       Lake
!        2       Reservoir
!        3       River
!        4       Freshwater Marsh, Floodplain
!        5       Swamp Forest, Flooded Forest
!        6       Coastal Wetland (incl. Mangrove, Estuary, Delta, Lagoon)
!        7       Pan, Brackish/Saline Wetland
!        8       Bog, Fen, Mire (Peatland)
!        9       Intermittent Wetland/Lake
!        10      50-100% Wetland
!        11      25-50% Wetland
!        12      Wetland Complex (0-25% Wetland)
!
!  5. Global Urban Extent
!      (http://sage.wisc.edu/urbanenvironment.html)
!
!  6. Global Cropping Types (ON GOING PROJECT ......)
!
!  Created by Yongjiu Dai, 12/2013
!-----------------------------------------------------------------------
USE MOD_Precision
USE spmd_TM
#if (defined usempi)
USE spmd_io
#endif

IMPLICIT NONE

! arguments:
      character(len=256), intent(in) :: dir_rawdata

! local variables:
      integer, parameter :: nlat=21600    ! 180*(60*2)
      integer, parameter :: nlon=43200    ! 360*(60*2)

      character(len=256) lndname

#if (defined usempi)
      character(len=1), allocatable :: land_chr1(:,:)
      character(len=2), allocatable :: land_chr2(:,:)
      integer(kind=1), allocatable ::  land_int1(:,:)
      integer(kind=2), allocatable ::  land_int2(:,:)
#else
      character(len=1) land_chr1(nlon)
      character(len=2) land_chr2(nlon)
      integer(kind=1)  land_int1(nlon)
      integer(kind=2)  land_int2(nlon)
#endif

      ! (1) global topography
      ! -----------------
      integer, allocatable :: elevation(:,:)

      ! (2) global land cover characteristics
      ! ---------------------------------
#if (defined LULC_USGS)
      integer, allocatable :: landtypes_usgs(:,:)  ! GLCC USGS land cover types
#endif
#if (defined LULC_IGBP)
      integer, allocatable :: landtypes_igbp(:,:)  ! MODIS IGBP land cover types
#endif

      ! (3) global glacier characteristics
      ! ------------------------------
      real(r8), allocatable :: glacier(:,:)            ! glacier coverage (%)

      ! (4) global lakes and wetlands characteristics
      ! -----------------------------------------
      integer, allocatable :: lakewetland(:,:)     ! land water and wetland types

      ! (5) global urban and build-up land characteristics
      ! -----------------------------------------
      integer, allocatable :: urban(:,:)           ! urban and built-up land

!   ---------------------------------------------------------------
      integer ia
      integer i, j, L
      integer nrow, ncol
      integer iunit
      integer length

      character c
      real(r8) a
      integer ii, iii, iiii, jj, jjj, jjjj
      integer nl, np

      integer, allocatable :: exclude(:)
      integer, allocatable :: ntmp(:,:)

      integer  :: buff_lb
      integer  :: buff_ub
      integer  :: nrow_start
      integer  :: nrow_end

      integer  :: int_min, int_max
      real(r8) :: r8_min, r8_max

#if (defined usempi)
      integer, parameter :: rmax = 100  ! searching radius (grid cells)
      type (map_type) :: buff_fine_lat_map
      integer(kind=MPI_OFFSET_KIND) :: fdisp
#endif


! Initialize MPI tasks
#if (defined usempi)
      nrow_start = fine_lat_map%bdisp(1) + 1
      nrow_end   = fine_lat_map%bdisp(1) + fine_lat_map%bstrd(1)
      buff_lb = max(nrow_start-rmax, 1)
      buff_ub = min(nrow_end  +rmax, nlat)

      allocate (buff_fine_lat_map%bstrd(1))
      allocate (buff_fine_lat_map%bdisp(1))

      buff_fine_lat_map%total = nlat
      buff_fine_lat_map%num   = buff_ub - buff_lb + 1

      buff_fine_lat_map%nblocks  = 1
      buff_fine_lat_map%bstrd(1) = buff_ub - buff_lb + 1
      buff_fine_lat_map%bdisp(1) = buff_lb - 1
#else
      nrow_start = 1
      nrow_end   = nlat
      buff_lb = 1
      buff_ub = nlat
#endif
!   ---------------------------------------------------------------
!
      allocate ( elevation      (nlon,buff_lb:buff_ub) ,&
#if (defined LULC_USGS)
                 landtypes_usgs (nlon,buff_lb:buff_ub) ,&
#endif
#if (defined LULC_IGBP)
                 landtypes_igbp (nlon,buff_lb:buff_ub) ,&
#endif

                 lakewetland    (nlon,buff_lb:buff_ub) ,&
                 glacier        (nlon,buff_lb:buff_ub) ,&
                 urban          (nlon,buff_lb:buff_ub)  )

! ----------------------------------------------------------------------
! ... (1) global digital elevation model (DEM) data
! ................................................
!     (1.1) elevation height (median)

      lndname = trim(dir_rawdata)//'terrain/dem-usgs.bin'
      IF (p_master) print*,trim(lndname)
#if (defined usempi)
      allocate (character(len=2) :: land_chr2 (nlon, buff_lb:buff_ub))
      fdisp = 0
      CALL mpi_rdwr_data (lndname, 'read', fdisp, buff_fine_lat_map, nlon*2, land_chr2)

      DO nrow = buff_lb, buff_ub
         DO ncol = 1, nlon
            elevation(ncol,nrow) = ia(land_chr2(ncol,nrow),2,-9999)
         ENDDO
      ENDDO
      deallocate (land_chr2)
#else
      iunit = 100
      inquire(iolength=length) land_chr2
      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
      DO nrow = 1, nlat
         read(iunit,rec=nrow,err=100) land_chr2
         DO ncol = 1, nlon
            elevation(ncol,nrow) = ia(land_chr2(ncol),2,-9999)
         ENDDO
      ENDDO
      close (iunit)
#endif

      ii=0
      DO nrow = nrow_start, nrow_end
         DO ncol = 1, nlon
             IF( elevation(ncol,nrow) < -9990)THEN
                 ii = ii + 1
             ENDIF
         ENDDO
      ENDDO

      int_min = minval(elevation(:,nrow_start:nrow_end))
      int_max = maxval(elevation(:,nrow_start:nrow_end))
#if (defined usempi)
      IF (p_iam_slave == 0) THEN
          CALL mpi_reduce (MPI_IN_PLACE, ii,   1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, int_min, 1, MPI_INTEGER, MPI_MIN, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, int_max, 1, MPI_INTEGER, MPI_MAX, 0, p_comm_slave, p_err)
      ELSE
          CALL mpi_reduce (ii,   0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (int_min, 0, 1, MPI_INTEGER, MPI_MIN, 0, p_comm_slave, p_err)
          CALL mpi_reduce (int_max, 0, 1, MPI_INTEGER, MPI_MAX, 0, p_comm_slave, p_err)
      ENDIF
#endif
      IF (p_master) print*,int_min,int_max
      IF (p_master) print*,'ii=', ii

! ........................................
! ... (2) gloabl land cover characteristics
! ........................................
!     (2.1) global land cover type (version 2.0) (USGS)

#if (defined LULC_USGS)
     ! GLCC USGS classification
     ! -------------------

      lndname = trim(dir_rawdata)//'landtypes/landtypes-usgs.bin'
      IF (p_master) print*,trim(lndname)

#if (defined usempi)
      allocate (land_chr1 (nlon, buff_lb:buff_ub))
      fdisp = 0
      CALL mpi_rdwr_data (lndname, 'read', fdisp, buff_fine_lat_map, nlon, land_chr1)
      landtypes_usgs = ichar(land_chr1)
      deallocate (land_chr1)
#else
      iunit = 100
      inquire(iolength=length) land_chr1
      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
      DO nrow = 1, nlat
         read(iunit,rec=nrow,err=100) land_chr1
         DO ncol = 1, nlon
            landtypes_usgs(ncol,nrow) = ichar(land_chr1(ncol))
         ENDDO
      ENDDO
      close (iunit)
#endif

      ii = 0
      iii = 0
      iiii = 0
      jjj = 0
      DO nrow = nrow_start, nrow_end
         DO ncol = 1, nlon
            IF(elevation(ncol,nrow) > -9990) THEN
               IF(landtypes_usgs(ncol,nrow) == 16) THEN
                  ii = ii + 1
               ENDIF
               IF(landtypes_usgs(ncol,nrow) == 17 .or. landtypes_usgs(ncol,nrow) == 18)THEN
                  landtypes_usgs(ncol,nrow) = 17  ! combine herbaceous / wooded wetland => permanent wetland
                  iii = iii + 1
               ENDIF
               IF(landtypes_usgs(ncol,nrow) == 24) THEN
                  iiii = iiii + 1
               ENDIF
               IF(landtypes_usgs(ncol,nrow) == 1) THEN
                  jjj = jjj + 1
               ENDIF
            ENDIF
            IF(elevation(ncol,nrow) < -9990) THEN
               landtypes_usgs(ncol,nrow) = 0
            ENDIF
         ENDDO
      ENDDO
      int_min = minval(landtypes_usgs(:,nrow_start:nrow_end))
      int_max = maxval(landtypes_usgs(:,nrow_start:nrow_end))
#if (defined usempi)
      IF (p_iam_slave == 0) THEN
          CALL mpi_reduce (MPI_IN_PLACE, ii,   1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, iii,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, iiii, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, jjj,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, int_min, 1, MPI_INTEGER,   MPI_MIN, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, int_max, 1, MPI_INTEGER,   MPI_MAX, 0, p_comm_slave, p_err)
      ELSE
          CALL mpi_reduce (ii,   0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (iii,  0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (iiii, 0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (jjj,  0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (int_min, 0, 1, MPI_INTEGER,   MPI_MIN, 0, p_comm_slave, p_err)
          CALL mpi_reduce (int_max, 0, 1, MPI_INTEGER,   MPI_MAX, 0, p_comm_slave, p_err)
      ENDIF
#endif

      IF (p_master) print*, int_min, int_max
      IF (p_master) print*,' USGS GLCC land cover '
      IF (p_master) print*,'land water points =', ii, 'wetland points=', iii, 'glacier points=', iiii, 'urban points=', jjj
#endif

#if (defined LULC_IGBP)
     ! MODIS IGBP classification
     ! -------------------

      lndname = trim(dir_rawdata)//'landtypes/landtypes-modis-igbp-2012.bin'
      IF (p_master) print*,trim(lndname)

#if (defined usempi)
      allocate (land_chr1 (nlon, buff_lb:buff_ub))
      fdisp = 0
      CALL mpi_rdwr_data (lndname, 'read', fdisp, buff_fine_lat_map, nlon, land_chr1)
      landtypes_igbp = ichar(land_chr1)
      deallocate (land_chr1)
#else
      iunit = 100
      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
      DO nrow = 1, nlat
         read(iunit,rec=nrow,err=100) land_chr1
         DO ncol = 1, nlon
            landtypes_igbp(ncol,nrow) = ichar(land_chr1(ncol))
         ENDDO
      ENDDO
      close (iunit)
#endif

      ii = 0
      iii = 0
      iiii = 0
      jjj = 0
      DO nrow = nrow_start, nrow_end
         DO ncol = 1, nlon
            IF(elevation(ncol,nrow) > -9990) THEN
               IF(landtypes_igbp(ncol,nrow) == 0) THEN
                 landtypes_igbp(ncol,nrow) = 17
                 ii = ii + 1
               ENDIF
               IF(landtypes_igbp(ncol,nrow) == 11) THEN
                 iii = iii + 1
               ENDIF
               IF(landtypes_igbp(ncol,nrow) == 15) THEN
                 iiii = iiii + 1
               ENDIF
               IF(landtypes_igbp(ncol,nrow) == 13) THEN
                 jjj = jjj + 1
               ENDIF
            ENDIF
            IF(elevation(ncol,nrow) < -9990) THEN
               landtypes_igbp(ncol,nrow) = 0
            ENDIF
         ENDDO
      ENDDO
      int_min = minval(landtypes_igbp(:,nrow_start:nrow_end))
      int_max = maxval(landtypes_igbp(:,nrow_start:nrow_end))
#if (defined usempi)
      IF (p_iam_slave == 0) THEN
          CALL mpi_reduce (MPI_IN_PLACE, ii,   1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, iii,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, iiii, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, jjj,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, int_min, 1, MPI_INTEGER,   MPI_MIN, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, int_max, 1, MPI_INTEGER,   MPI_MAX, 0, p_comm_slave, p_err)
      ELSE
          CALL mpi_reduce (ii,   0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (iii,  0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (iiii, 0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (jjj,  0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (int_min, 0, 1, MPI_INTEGER,   MPI_MIN, 0, p_comm_slave, p_err)
          CALL mpi_reduce (int_max, 0, 1, MPI_INTEGER,   MPI_MAX, 0, p_comm_slave, p_err)
      ENDIF
#endif

      IF (p_master) print*, int_min, int_max
      IF (p_master) print*,' MODIS IGBP land cover '
      IF (p_master) print*,'land water points =', ii, 'wetland points=', iii, 'glacier points=', iiii, 'urban points=', jjj
#endif

! ................................................
! ... (3) global lakes and wetland characteristics
! ................................................
!     (3.1) global lakes and wetland

      lndname = trim(dir_rawdata)//'lake_wetland/glwd.bin'
      IF (p_master) print*,trim(lndname)

#if (defined usempi)
      allocate (land_chr1 (nlon, buff_lb:buff_ub))
      fdisp = 0
      CALL mpi_rdwr_data (lndname, 'read', fdisp, buff_fine_lat_map, nlon, land_chr1)
      lakewetland = ichar(land_chr1)
      deallocate (land_chr1)
#else
      iunit = 100
      inquire(iolength=length) land_chr1
      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
      DO nrow = 1, nlat
         read(iunit,rec=nrow,err=100) land_chr1
         DO ncol = 1, nlon
            lakewetland(ncol,nrow) = ichar(land_chr1(ncol))
         ENDDO
      ENDDO
      close (iunit)
#endif

      ii = 0
      iii = 0
      DO nrow = nrow_start, nrow_end
         DO ncol = 1, nlon
            IF(elevation(ncol,nrow)>-9990) THEN
               ! Replace GLCC_USGS and MODIS_IGBP water bodies with GLWD classification
               IF(lakewetland(ncol,nrow) <= 3) THEN
#if (defined LULC_USGS)
                  landtypes_usgs(ncol,nrow) = 16
#endif
#if (defined LULC_IGBP)
                  landtypes_igbp(ncol,nrow) = 17
#endif
                  ii = ii + 1
               ENDIF
               ! Replace GLCC_USGS and MODIS_IGBP wetland with GLWD classification
! modified by dai, 06/02/2016
               !IF(lakewetland(ncol,nrow) > 3 .and. lakewetland(ncol,nrow) < 13)THEN
! deleted by dai, 08/15/2016
!               IF(lakewetland(ncol,nrow) > 3 .and. lakewetland(ncol,nrow) < 11)THEN
!#if (defined LULC_USGS)
!                  landtypes_usgs(ncol,nrow) = 17
!#endif
!#if (defined LULC_IGBP)
!                  landtypes_igbp(ncol,nrow) = 11
!#endif
!                  iii = iii + 1
!               ENDIF
               IF(nrow > 18000)THEN  ! (90N + 60S) * 120
                  lakewetland(ncol,nrow) = 255
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      int_min = minval(lakewetland(:,nrow_start:nrow_end))
      int_max = maxval(lakewetland(:,nrow_start:nrow_end))
#if (defined usempi)
      IF (p_iam_slave == 0) THEN
          CALL mpi_reduce (MPI_IN_PLACE, ii,   1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, iii,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, int_min, 1, MPI_INTEGER,   MPI_MIN, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, int_max, 1, MPI_INTEGER,   MPI_MAX, 0, p_comm_slave, p_err)
      ELSE
          CALL mpi_reduce (ii,   0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (iii,  0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (int_min, 0, 1, MPI_INTEGER,   MPI_MIN, 0, p_comm_slave, p_err)
          CALL mpi_reduce (int_max, 0, 1, MPI_INTEGER,   MPI_MAX, 0, p_comm_slave, p_err)
      ENDIF
#endif

      IF (p_master) print*,int_min, int_max
      IF (p_master) print*,'GLWD'
      IF (p_master) print*,'land water (1-3) =', ii, 'wetland (4-12) =',iii

! ......................................
! ... (4) global glacier and ice sheet characteristics
! ......................................

      lndname = trim(dir_rawdata)//'glacier/glacier.bin'
      IF (p_master) print*,trim(lndname)

#if (defined usempi)
      allocate (land_int2 (nlon, buff_lb:buff_ub))
      fdisp = 0
      CALL mpi_rdwr_data (lndname, 'read', fdisp, buff_fine_lat_map, nlon, land_int2)
      glacier = land_int2 * 0.1
      deallocate (land_int2)
#else
      iunit = 100
      inquire(iolength=length) land_int2
      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
      DO nrow = 1, nlat
         read(iunit,rec=nrow,err=100) land_int2
         DO ncol = 1, nlon
            glacier(ncol,nrow) = land_int2(ncol) * 0.1
         ENDDO
      ENDDO
      close (iunit)
#endif

      ii = 0
      DO nrow = nrow_start, nrow_end
         DO ncol = 1, nlon
            IF(elevation(ncol,nrow)>-9990) THEN
               ! Replace GLCC and MODIS glacier/ice sheet with glacier specific dataset
               IF(glacier(ncol,nrow) > 0.5 .and. glacier(ncol,nrow) < 101.)THEN
#if (defined LULC_USGS)
                  landtypes_usgs(ncol,nrow) = 24
#endif
#if (defined LULC_IGBP)
                  landtypes_igbp(ncol,nrow) = 15
#endif
                  ii = ii + 1
               ENDIF
               ! Antarctic (ice sheet / barren ONLY)
               IF(nrow > 18000)THEN  ! (90N + 60S) * 120
#if (defined LULC_USGS)
                  IF(landtypes_usgs(ncol,nrow)/=23)THEN
                     landtypes_usgs(ncol,nrow) = 24
                     glacier(ncol,nrow) = 100.
                  ENDIF
#endif
#if (defined LULC_IGBP)
                  IF(landtypes_igbp(ncol,nrow)/=16)THEN
                     landtypes_igbp(ncol,nrow) = 15
                     glacier(ncol,nrow) = 100.
                  ENDIF
#endif
               ENDIF
! modified by dai, 06/02/2016
!               ! Greenland  (ice sheet / barren / built-up ONLY)
!               ! BETWEEN (59-83N, 74-11W) =>| 0<nrow<(90-31)*120, (180-74)*120<ncol<(180-11)*120
!               IF(nrow<3720 .and. (ncol>12720 .and. ncol<20280))THEN
!#if (defined LULC_USGS)
!                  IF(landtypes_usgs(ncol,nrow)/=23 .or. landtypes_usgs(ncol,nrow)/=1)THEN
!                     landtypes_usgs(ncol,nrow) = 24
!                     glacier(ncol,nrow) = 100.
!                  ENDIF
!#endif
!#if (defined LULC_IGBP)
!                  IF(landtypes_igbp(ncol,nrow)/=16 .or. landtypes_igbp(ncol,nrow)/=13)THEN
!                     landtypes_igbp(ncol,nrow) = 15
!                     glacier(ncol,nrow) = 100.
!                  ENDIF
!#endif
!               ENDIF
            ENDIF
         ENDDO
      ENDDO


      r8_min = minval(glacier(:,nrow_start:nrow_end))
      r8_max = maxval(glacier(:,nrow_start:nrow_end))
#if (defined usempi)
      IF (p_iam_slave == 0) THEN
          CALL mpi_reduce (MPI_IN_PLACE, ii,   1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, r8_min, 1, MPI_INTEGER,   MPI_MIN, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, r8_max, 1, MPI_INTEGER,   MPI_MAX, 0, p_comm_slave, p_err)
      ELSE
          CALL mpi_reduce (ii,   0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (r8_min, 0, 1, MPI_INTEGER,   MPI_MIN, 0, p_comm_slave, p_err)
          CALL mpi_reduce (r8_max, 0, 1, MPI_INTEGER,   MPI_MAX, 0, p_comm_slave, p_err)
      ENDIF
#endif

      IF (p_master) print*, r8_min, r8_max
      IF (p_master) print*,'DATA glacier points'
      IF (p_master) print*,'glacier (> 0.5%) =', ii

! ....................................
! ... (5) global urban characteristics
! ....................................

      lndname = trim(dir_rawdata)//'urban/urban-builtup.bin'
      IF (p_master) print*,trim(lndname)

#if (defined usempi)
      allocate (land_chr1 (nlon, buff_lb:buff_ub))
      fdisp = 0
      CALL mpi_rdwr_data (lndname, 'read', fdisp, buff_fine_lat_map, nlon, land_chr1)
      urban = ichar(land_chr1)
      deallocate (land_chr1)
#else
      iunit = 100
      inquire(iolength=length) land_chr1
      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
      DO nrow = 1, nlat
         read(iunit,rec=nrow,err=100) land_chr1
         DO ncol = 1, nlon
            urban(ncol,nrow) = ichar(land_chr1(ncol))
         ENDDO
      ENDDO
      close (iunit)
#endif

      ii = 0
      DO nrow = nrow_start, nrow_end
         DO ncol = 1, nlon
            IF(elevation(ncol,nrow)>-9990) THEN
               ! Replace GLCC and MODIS urban/built-up with urban specific dataset
               IF(urban(ncol,nrow) == 1)THEN
#if (defined LULC_USGS)
                  landtypes_usgs(ncol,nrow) = 1
#endif
#if (defined LULC_IGBP)
                  landtypes_igbp(ncol,nrow) = 13
#endif
                  ii = ii + 1
               ENDIF
               IF(nrow > 18000)THEN  ! (90N + 60S) * 120
                  urban(ncol,nrow) = 0
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      int_min = minval(urban(:,nrow_start:nrow_end))
      int_max = maxval(urban(:,nrow_start:nrow_end))
#if (defined usempi)
      IF (p_iam_slave == 0) THEN
          CALL mpi_reduce (MPI_IN_PLACE, ii,   1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, int_min, 1, MPI_INTEGER,   MPI_MIN, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, int_max, 1, MPI_INTEGER,   MPI_MAX, 0, p_comm_slave, p_err)
      ELSE
          CALL mpi_reduce (ii,   0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (int_min, 0, 1, MPI_INTEGER,   MPI_MIN, 0, p_comm_slave, p_err)
          CALL mpi_reduce (int_max, 0, 1, MPI_INTEGER,   MPI_MAX, 0, p_comm_slave, p_err)
      ENDIF
#endif

      IF (p_master) print*, int_min, int_max
      IF (p_master) print*,'URBAN'
      IF (p_master) print*,'Urban and Built-up Land Points =', ii

#if (defined CoLMDEBUG)
#if (defined LULC_USGS)
      ii = 0
      iii = 0
      iiii = 0
      jjj = 0

      DO nrow = nrow_start, nrow_end
         DO ncol = 1, nlon
            IF(elevation(ncol,nrow) > -9990) THEN
               IF(landtypes_usgs(ncol,nrow) == 16) THEN
                  ii = ii + 1
               ENDIF
               IF(landtypes_usgs(ncol,nrow) == 17 .or. landtypes_usgs(ncol,nrow) == 18)THEN
                  iii = iii + 1
               ENDIF
               IF(landtypes_usgs(ncol,nrow) == 24) THEN
                  iiii = iiii + 1
               ENDIF
               IF(landtypes_usgs(ncol,nrow) == 1) THEN
                  jjj = jjj + 1
               ENDIF
            ENDIF
         ENDDO
      ENDDO
#if (defined usempi)
      IF (p_iam_slave == 0) THEN
          CALL mpi_reduce (MPI_IN_PLACE, ii,   1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, iii,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, iiii, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, jjj,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
      ELSE
          CALL mpi_reduce (ii,  0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (iii, 0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (iiii,0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (jjj, 0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
      ENDIF
#endif
      IF (p_master) print*,' SECOND USGS GLCC land cover '
      IF (p_master) print*,'land water points =', ii, 'wetland points=', iii, 'glacier points=', iiii, 'urban points=', jjj
#endif

#if (defined LULC_IGBP)
      ii = 0
      iii = 0
      iiii = 0
      jjj = 0

      DO nrow = nrow_start, nrow_end
         DO ncol = 1, nlon
            IF(elevation(ncol,nrow) > -9990) THEN
               IF(landtypes_igbp(ncol,nrow) == 17) THEN
                 ii = ii + 1
               ENDIF
               IF(landtypes_igbp(ncol,nrow) == 11) THEN
                 iii = iii + 1
               ENDIF
               IF(landtypes_igbp(ncol,nrow) == 15) THEN
                 iiii = iiii + 1
               ENDIF
               IF(landtypes_igbp(ncol,nrow) == 13) THEN
                 jjj = jjj + 1
               ENDIF
            ENDIF
         ENDDO
      ENDDO
#if (defined usempi)
      IF (p_iam_slave == 0) THEN
          CALL mpi_reduce (MPI_IN_PLACE, ii,   1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, iii,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, iiii, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, jjj,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
      ELSE
          CALL mpi_reduce (ii,  0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (iii, 0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (iiii,0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (jjj, 0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
      ENDIF
#endif
      IF (p_master) print*,' SECOND MODIS IGBP land cover '
      IF (p_master) print*,'land water points =', ii, 'wetland points=', iii, 'glacier points=', iiii, 'urban points=', jjj
#endif
#endif

#if (defined usempi)
#if (defined LULC_USGS)
      CALL update_buff (landtypes_usgs, nlon, nlat, buff_lb, nrow_start, nrow_end, buff_ub)
#endif
#if (defined LULC_IGBP)
      CALL update_buff (landtypes_igbp, nlon, nlat, buff_lb, nrow_start, nrow_end, buff_ub)
#endif
#endif
!---------------------------------------------------------------------------------------------
! Correct the land cover types on which the types (GLCC_USGS/MODIS_IGBP) were classified
! water bodies/wetland/glacier/urban, but NOT classified in GLWD/GLACIER/URBAN specific dataset
!---------------------------------------------------------------------------------------------
      allocate (ntmp(nlon,buff_lb:buff_ub))
!      ! WATER BODIES
!      ! ------------
!#if (defined LULC_USGS)
!      iii = 0
!
!      ! GLCC USGS land cover classification (WATER BODIES)
!      nl = 24
!!     np = 5
!!     allocate (exclude(5))
!!     exclude = (/1,16,17,18,24/) !/urban and built-up(1),water bodies(16),herbaceous wetland(17),wooded wetland(18),snow and ice(24)/
!      np = 4
!      allocate (exclude(4))
!      exclude = (/1,16,17,18/)
!      ntmp = landtypes_usgs
!      DO j = nrow_start, nrow_end
!         DO i = 1, nlon
!            IF(elevation(i,j)>-9990) THEN
!               IF(lakewetland(i,j)>3 .and. landtypes_usgs(i,j)==16) THEN
!                  CALL substitute(i,j,nlon,buff_lb,buff_ub,nl,np,landtypes_usgs,exclude,L)
!                  ntmp(i,j) = L
!                  iii = iii + 1
!               ENDIF
!            ENDIF
!         ENDDO
!      ENDDO
!      landtypes_usgs = ntmp
!      deallocate (exclude)
!#if (defined usempi)
!      CALL update_buff (landtypes_usgs, nlon, nlat, buff_lb, nrow_start, nrow_end, buff_ub)
!      IF (p_iam_slave == 0) THEN
!          CALL mpi_reduce (MPI_IN_PLACE, iii,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
!      ELSE
!          CALL mpi_reduce (iii, 0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
!      ENDIF
!#endif
!      IF (p_master) print*, 'GLCC WATER BODIES','iii=',iii
!#endif
!
!#if (defined LULC_IGBP)
!      iii = 0
!
!      ! MODIS IGBP land cover classification (WATER BODIES)
!      nl = 17
!!     np = 4
!!     allocate (exclude(4))
!!     exclude = (/11,13,15,17/)  !/permanent wetland(11),urban and built-up(13),snow and ice(15),water bodies(17)/
!      np = 3
!      allocate (exclude(3))
!      exclude = (/11,13,17/)
!      ntmp = landtypes_igbp
!      DO j = nrow_start, nrow_end
!         DO i = 1, nlon
!            IF(elevation(i,j)>-9990) THEN
!               IF(lakewetland(i,j)>3 .and. landtypes_igbp(i,j)==17) THEN
!                  CALL substitute(i,j,nlon,buff_lb,buff_ub,nl,np,landtypes_igbp,exclude,L)
!                  ntmp(i,j) = L
!                  iii = iii + 1
!               ENDIF
!            ENDIF
!         ENDDO
!      ENDDO
!      landtypes_igbp = ntmp
!      deallocate (exclude)
!#if (defined usempi)
!      CALL update_buff (landtypes_igbp, nlon, nlat, buff_lb, nrow_start, nrow_end, buff_ub)
!      IF (p_iam_slave == 0) THEN
!          CALL mpi_reduce (MPI_IN_PLACE, iii,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
!      ELSE
!          CALL mpi_reduce (iii, 0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
!      ENDIF
!#endif
!      IF (p_master) print*, 'MODIS IGBP WATER BODIES','iii=',iii
!#endif
!
!#if (defined CoLMDEBUG)
!#if (defined LULC_USGS)
!      iiii = 0
!
!      DO j = nrow_start, nrow_end
!         DO i = 1, nlon
!            IF(elevation(i,j)>-9990) THEN
!               IF(lakewetland(i,j)>3 .and. landtypes_usgs(i,j)==16) THEN
!                  iiii = iiii + 1
!               ENDIF
!            ENDIF
!         ENDDO
!      ENDDO
!#if (defined usempi)
!      IF (p_iam_slave == 0) THEN
!          CALL mpi_reduce (MPI_IN_PLACE, iiii,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
!      ELSE
!          CALL mpi_reduce (iiii, 0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
!      ENDIF
!#endif
!      IF (p_master) print*, 'GLCC WATER BODIES','iiii=',iiii
!#endif
!#if (defined LULC_IGBP)
!      iiii = 0
!
!      DO j = nrow_start, nrow_end
!         DO i = 1, nlon
!            IF(elevation(i,j)>-9990) THEN
!               IF(lakewetland(i,j)>3 .and. landtypes_igbp(i,j)==17) THEN
!                  iiii = iiii + 1
!               ENDIF
!            ENDIF
!         ENDDO
!      ENDDO
!#if (defined usempi)
!      IF (p_iam_slave == 0) THEN
!          CALL mpi_reduce (MPI_IN_PLACE, iiii,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
!      ELSE
!          CALL mpi_reduce (iiii, 0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
!      ENDIF
!#endif
!      IF (p_master) print*, 'MODIS IGBP WATER BODIES','iiii=',iiii
!#endif
!#endif

! deleted by dai, 07/27/2016
! The Global Lake and Wetlands Types (http://www.wwfus.org/science/data.cfm)
! may have some problems in presenting the wetland type in some regions.
!
!      ! WETLAND
!      ! -------
!#if (defined LULC_USGS)
!      iii = 0
!
!      ! GLCC USGS land cover classification (WETLAND)
!      nl = 24
!      np = 5
!      allocate (exclude(5))
!      exclude = (/1,16,17,18,24/)  !/urban and built-up(1),water bodies(16),herbaceous wetland(17),wooded wetland(18),snow and ice(24)/
!      ntmp = landtypes_usgs
!      DO j = nrow_start, nrow_end
!         DO i = 1, nlon
!            IF(elevation(i,j)>-9990) THEN
!               IF((lakewetland(i,j)<=3 .or. lakewetland(i,j)>12) .and. landtypes_usgs(i,j)==17) THEN
!                  CALL substitute(i,j,nlon,buff_lb,buff_ub,nl,np,landtypes_usgs,exclude,L)
!                  ntmp(i,j) = L
!                  iii = iii + 1
!               ENDIF
!            ENDIF
!         ENDDO
!      ENDDO
!      landtypes_usgs = ntmp
!      deallocate (exclude)
!#if (defined usempi)
!      CALL update_buff (landtypes_usgs, nlon, nlat, buff_lb, nrow_start, nrow_end, buff_ub)
!      IF (p_iam_slave == 0) THEN
!          CALL mpi_reduce (MPI_IN_PLACE, iii,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
!      ELSE
!          CALL mpi_reduce (iii, 0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
!      ENDIF
!#endif
!      IF (p_master) print*, 'GLCC WETLAND','iii=',iii
!#endif
!
!#if (defined LULC_IGBP)
!      iii = 0
!
!      ! MODIS land cover classification (WETLAND)
!      nl = 17
!      np = 4
!      allocate (exclude(4))
!      exclude = (/11,13,15,17/)   !/permanent wetland(11),urban and built-up(13),snow and ice(15),water bodies(17)/
!      ntmp = landtypes_igbp
!      DO j = nrow_start, nrow_end
!         DO i = 1, nlon
!            IF(elevation(i,j)>-9990) THEN
!               IF((lakewetland(i,j)<=3 .or. lakewetland(i,j)>12) .and. landtypes_igbp(i,j)==11) THEN
!                  CALL substitute(i,j,nlon,buff_lb,buff_ub,nl,np,landtypes_igbp,exclude,L)
!                  ntmp(i,j) = L
!                  iii = iii + 1
!               ENDIF
!            ENDIF
!         ENDDO
!      ENDDO
!      landtypes_igbp = ntmp
!      deallocate (exclude)
!#if (defined usempi)
!      CALL update_buff (landtypes_igbp, nlon, nlat, buff_lb, nrow_start, nrow_end, buff_ub)
!      IF (p_iam_slave == 0) THEN
!          CALL mpi_reduce (MPI_IN_PLACE, iii,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
!      ELSE
!          CALL mpi_reduce (iii, 0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
!      ENDIF
!#endif
!      IF (p_master) print*, 'MODIS IGBP WETLAND','iii=',iii
!#endif
!
!#if (defined CoLMDEBUG)
!#if (defined LULC_USGS)
!      iiii = 0
!
!      DO j = nrow_start, nrow_end
!         DO i = 1, nlon
!            IF(elevation(i,j)>-9990) THEN
!               IF((lakewetland(i,j)<=3 .or. lakewetland(i,j)>12) .and. landtypes_usgs(i,j)==17) THEN
!                  iiii = iiii + 1
!               ENDIF
!            ENDIF
!         ENDDO
!      ENDDO
!#if (defined usempi)
!      IF (p_iam_slave == 0) THEN
!          CALL mpi_reduce (MPI_IN_PLACE, iiii,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
!      ELSE
!          CALL mpi_reduce (iiii, 0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
!      ENDIF
!#endif
!      IF (p_master) print*, 'GLCC WETLAND', 'iiii=',iiii
!#endif
!#if (defined LULC_IGBP)
!      iiii = 0
!
!      DO j = nrow_start, nrow_end
!         DO i = 1, nlon
!            IF(elevation(i,j)>-9990) THEN
!               IF((lakewetland(i,j)<=3 .or. lakewetland(i,j)>12) .and. landtypes_igbp(i,j)==11) THEN
!                  iiii = iiii + 1
!               ENDIF
!            ENDIF
!         ENDDO
!      ENDDO
!#if (defined usempi)
!      IF (p_iam_slave == 0) THEN
!          CALL mpi_reduce (MPI_IN_PLACE, iiii,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
!      ELSE
!          CALL mpi_reduce (iiii, 0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
!      ENDIF
!#endif
!      IF (p_master) print*, 'MODIS IGBP WETLAND', 'iiii=',iiii
!#endif
!#endif

      ! GLACIER/ICESHEET
      ! ------------------
#if (defined LULC_USGS)
      iii = 0

      ! GLCC USGS land cover classification (GLACIER/ICESHEET)
      nl = 24
!     np = 5
!     allocate (exclude(5))
!     exclude = (/1,16,17,18,24/) !/urban and built-up(1),water bodies(16),herbaceous wetland(17),wooded wetland(18),snow and ice(24)/
      np = 4
      allocate (exclude(4))
      exclude = (/1,16,17,18/)
      ntmp = landtypes_usgs
      DO j = nrow_start, nrow_end
         DO i = 1, nlon
            IF(elevation(i,j)>-9990) THEN
               IF((glacier(i,j)<=0.5 .or. glacier(i,j)> 100.1) .and. landtypes_usgs(i,j)==24) THEN
                  CALL substitute(i,j,nlon,buff_lb,buff_ub,nl,np,landtypes_usgs,exclude,L)
                  ntmp(i,j) = L
                  iii = iii + 1
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      landtypes_usgs = ntmp
      deallocate (exclude)
#if (defined usempi)
      CALL update_buff (landtypes_usgs, nlon, nlat, buff_lb, nrow_start, nrow_end, buff_ub)
      IF (p_iam_slave == 0) THEN
          CALL mpi_reduce (MPI_IN_PLACE, iii,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
      ELSE
          CALL mpi_reduce (iii, 0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
      ENDIF
#endif
      IF (p_master) print*, 'GLCC GLACIER/ICESHEET','iii=',iii
#endif
#if (defined LULC_IGBP)
      iii = 0

      ! MODIS IGBP land cover classification (GLACIER/ICESHEET)
      nl = 17
!     np = 4
!     allocate (exclude(4))
!     exclude = (/11,13,15,17/)  !/permanent wetland(11),urban and built-up(13),snow and ice(15),water bodies(17)/
      np = 3
      allocate (exclude(3))
      exclude = (/11,13,17/)
      ntmp = landtypes_igbp
      DO j = nrow_start, nrow_end
         DO i = 1, nlon
            IF(elevation(i,j)>-9990) THEN
               IF((glacier(i,j)<=0.5 .or. glacier(i,j)>100.1) .and. landtypes_igbp(i,j)==15) THEN
                  CALL substitute(i,j,nlon,buff_lb,buff_ub,nl,np,landtypes_igbp,exclude,L)
                  ntmp(i,j) = L
                  iii = iii + 1
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      landtypes_igbp = ntmp
      deallocate (exclude)
#if (defined usempi)
      CALL update_buff (landtypes_igbp, nlon, nlat, buff_lb, nrow_start, nrow_end, buff_ub)
      IF (p_iam_slave == 0) THEN
          CALL mpi_reduce (MPI_IN_PLACE, iii,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
      ELSE
          CALL mpi_reduce (iii, 0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
      ENDIF
#endif
      IF (p_master) print*, 'MODIS IGBP GLACIER/ICESHEET','iii=',iii
#endif

#if (defined CoLMDEBUG)
#if (defined LULC_USGS)
      iiii = 0

      DO j = nrow_start, nrow_end
         DO i = 1, nlon
            IF(elevation(i,j)>-9990) THEN
               IF((glacier(i,j)<=0.5 .or. glacier(i,j)> 100.1) .and. landtypes_usgs(i,j)==24) THEN
                  iiii = iiii + 1
               ENDIF
            ENDIF
         ENDDO
      ENDDO
#if (defined usempi)
      IF (p_iam_slave == 0) THEN
          CALL mpi_reduce (MPI_IN_PLACE, iiii,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
      ELSE
          CALL mpi_reduce (iiii, 0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
      ENDIF
#endif
      IF (p_master) print*, 'GLCC GLACIER/ICESHEET', 'iiii=',iiii
#endif
#if (defined LULC_IGBP)
      iiii = 0

      DO j = nrow_start, nrow_end
         DO i = 1, nlon
            IF(elevation(i,j)>-9990) THEN
               IF((glacier(i,j)<=0.5 .or. glacier(i,j)>100.1) .and. landtypes_igbp(i,j)==15) THEN
                  iiii = iiii + 1
               ENDIF
            ENDIF
         ENDDO
      ENDDO
#if (defined usempi)
      IF (p_iam_slave == 0) THEN
          CALL mpi_reduce (MPI_IN_PLACE, iiii,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
      ELSE
          CALL mpi_reduce (iiii, 0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
      ENDIF
#endif
      IF (p_master) print*, 'MODIS IGBP GLACIER/ICESHEET', 'iiii=',iiii
#endif
#endif

      ! URBAN and BUILT-UP LAND
      ! ------------------
#if (defined LULC_USGS)
      iii = 0

      ! GLCC USGS land cover classification (URBAN and BUILT-UP LAND)
      nl = 24
!     np = 5
!     allocate (exclude(5))
!     exclude = (/1,16,17,18,24/) !/urban and built-up(1),water bodies(16),herbaceous wetland(17),wooded wetland(18),snow and ice(24)/
      np = 3
      allocate (exclude(3))
      exclude = (/1,16,24/)
      ntmp = landtypes_usgs
      DO j = nrow_start, nrow_end
         DO i = 1, nlon
            IF(elevation(i,j)>-9990) THEN
               IF(urban(i,j)==0 .and. landtypes_usgs(i,j)==1) THEN
                  CALL substitute(i,j,nlon,buff_lb,buff_ub,nl,np,landtypes_usgs,exclude,L)
                  ntmp(i,j) = L
                  iii = iii + 1
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      landtypes_usgs = ntmp
      deallocate (exclude)
#if (defined usempi)
      CALL update_buff (landtypes_usgs, nlon, nlat, buff_lb, nrow_start, nrow_end, buff_ub)
      IF (p_iam_slave == 0) THEN
          CALL mpi_reduce (MPI_IN_PLACE, iii,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
      ELSE
          CALL mpi_reduce (iii, 0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
      ENDIF
#endif
      IF (p_master) print*, 'GLCC URBAN and BUILT-UP LAND','iii=',iii
#endif

#if (defined LULC_IGBP)
      iii = 0

      ! MODIS IGBP land cover classification (URBAN and BUILT-UP LAND)
      nl = 17
!     np = 4
!     allocate (exclude(4))
!     exclude = (/11,13,15,17/)  !/permanent wetland(11),urban and built-up(13),snow and ice(15),water bodies(17)/
      np = 3
      allocate (exclude(3))
      exclude = (/13,15,17/)
      ntmp = landtypes_igbp
      DO j = nrow_start, nrow_end
         DO i = 1, nlon
            IF(elevation(i,j)>-9990) THEN
               IF(urban(i,j)==0 .and. landtypes_igbp(i,j)==13) THEN
                  CALL substitute(i,j,nlon,buff_lb,buff_ub,nl,np,landtypes_igbp,exclude,L)
                  ntmp(i,j) = L
                  iii = iii + 1
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      landtypes_igbp = ntmp
      deallocate (exclude)
#if (defined usempi)
      CALL update_buff (landtypes_igbp, nlon, nlat, buff_lb, nrow_start, nrow_end, buff_ub)
      IF (p_iam_slave == 0) THEN
          CALL mpi_reduce (MPI_IN_PLACE, iii,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
      ELSE
          CALL mpi_reduce (iii, 0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
      ENDIF
#endif
      IF (p_master) print*, 'MODIS IGBP URBAN and BUILT-UP LAND','iii=',iii
#endif

#if (defined CoLMDEBUG)
#if (defined LULC_USGS)
      iiii = 0

      DO j = nrow_start, nrow_end
         DO i = 1, nlon
            IF(elevation(i,j)>-9990) THEN
               IF(urban(i,j)==0 .and. landtypes_usgs(i,j)==1) THEN
                  iiii = iiii + 1
               ENDIF
            ENDIF
         ENDDO
      ENDDO
#if (defined usempi)
      IF (p_iam_slave == 0) THEN
          CALL mpi_reduce (MPI_IN_PLACE, iiii,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
      ELSE
          CALL mpi_reduce (iiii, 0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
      ENDIF
#endif
      IF (p_master) print*, 'GLCC URBAN and BUILT-UP LAND', 'iiii=',iiii
#endif

#if (defined LULC_IGBP)
      iiii = 0
      DO j = nrow_start, nrow_end
         DO i = 1, nlon
            IF(elevation(i,j)>-9990) THEN
               IF(urban(i,j)==0 .and. landtypes_igbp(i,j)==13) THEN
                  iiii = iiii + 1
               ENDIF
            ENDIF
         ENDDO
      ENDDO
#if (defined usempi)
      IF (p_iam_slave == 0) THEN
          CALL mpi_reduce (MPI_IN_PLACE, iiii,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
      ELSE
          CALL mpi_reduce (iiii, 0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
      ENDIF
#endif
      IF (p_master) print*, 'MODIS IGBP URBAN and BUILT-UP LAND', 'iiii=',iiii
#endif
#endif

      !--------------------------------------------------------------------
#if (defined CoLMDEBUG)
#if (defined LULC_USGS)
      ii = 0
      iii = 0
      iiii = 0
      jj = 0
      jjj = 0
      jjjj = 0

      DO nrow = nrow_start, nrow_end
         DO ncol = 1, nlon
            IF(elevation(ncol,nrow) > -9990) THEN
               IF(landtypes_usgs(ncol,nrow) == 16) THEN
                  ii = ii + 1
               ENDIF
               IF(landtypes_usgs(ncol,nrow) == 17 .or. landtypes_usgs(ncol,nrow) == 18)THEN
                  iii = iii + 1
               ENDIF
               IF(landtypes_usgs(ncol,nrow) == 24) THEN
                  iiii = iiii + 1
               ENDIF
               IF(landtypes_usgs(ncol,nrow) == 1) THEN
                  jjj = jjj + 1
               ENDIF
            ENDIF
         ENDDO
      ENDDO
#if (defined usempi)
      IF (p_iam_slave == 0) THEN
          CALL mpi_reduce (MPI_IN_PLACE, ii,   1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, iii,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, iiii, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, jjj,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
      ELSE
          CALL mpi_reduce (ii,   0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (iii,  0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (iiii, 0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (jjj,  0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
      ENDIF
#endif
      IF (p_master) print*,' FINAL USGS GLCC land cover '
      IF (p_master) print*,'land water points =', ii, 'wetland points=', iii, 'glacier points=', iiii, 'urban points=', jjj
#endif

#if (defined LULC_IGBP)
      ii = 0
      iii = 0
      iiii = 0
      jj = 0
      jjj = 0
      jjjj = 0

      DO nrow = nrow_start, nrow_end
         DO ncol = 1, nlon
            IF(elevation(ncol,nrow) > -9990) THEN
               IF(landtypes_igbp(ncol,nrow) == 17) THEN
                 ii = ii + 1
               ENDIF
               IF(landtypes_igbp(ncol,nrow) == 11) THEN
                 iii = iii + 1
               ENDIF
               IF(landtypes_igbp(ncol,nrow) == 15) THEN
                 iiii = iiii + 1
               ENDIF
               IF(landtypes_igbp(ncol,nrow) == 13) THEN
                 jjj = jjj + 1
               ENDIF
            ENDIF
         ENDDO
      ENDDO
#if (defined usempi)
      IF (p_iam_slave == 0) THEN
          CALL mpi_reduce (MPI_IN_PLACE, ii,   1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, iii,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, iiii, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (MPI_IN_PLACE, jjj,  1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
      ELSE
          CALL mpi_reduce (ii,   0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (iii,  0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (iiii, 0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
          CALL mpi_reduce (jjj,  0, 1, MPI_INTEGER, MPI_SUM, 0, p_comm_slave, p_err)
      ENDIF
#endif
      IF (p_master) print*,' FINAL MODIS IGBP land cover '
      IF (p_master) print*,'land water points =', ii, 'wetland points=', iii, 'glacier points=', iiii, 'urban points=', jjj
#endif
#endif

!#if (defined ongoing)
! ..............................................
! ... (6) global cultural characteristics (crops)
! ..............................................
!
!#endif

! Write out the land cover types
#if (defined LULC_USGS)
     ! GLCC USGS classification
     ! -------------------
      lndname = trim(dir_rawdata)//'RAW_DATA_updated/landtypes_usgs_update.bin'
      IF (p_master) print*,trim(lndname)
#if (defined usempi)
      allocate (land_chr1 (nlon, nrow_start:nrow_end))
      land_chr1 = char(landtypes_usgs(:,nrow_start:nrow_end))
      fdisp = 0
      CALL mpi_rdwr_data (lndname, 'write', fdisp, fine_lat_map, nlon, land_chr1)
      deallocate (land_chr1)
#else
      iunit = 100
      inquire(iolength=length) land_chr1
      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='unknown')
      DO nrow = 1, nlat
! modifiedy by yuan, 06/02/2016
         !DO ncol = 1, nlon
         !   land_chr1(ncol) = char(landtypes_usgs(ncol,nrow))
         !ENDDO
         land_chr1(:) = char(landtypes_usgs(:,nrow))
         write(iunit,rec=nrow,err=100) land_chr1
      ENDDO
      close (iunit)
#endif
#endif
#if (defined LULC_IGBP)
     ! MODIS IGBP classification
     ! -------------------
      lndname = trim(dir_rawdata)//'RAW_DATA_updated/landtypes_igbp_update.bin'
      IF (p_master) print*,trim(lndname)
#if (defined usempi)
      allocate (land_chr1 (nlon, nrow_start:nrow_end))
      land_chr1 = char(landtypes_igbp(:,nrow_start:nrow_end))
      fdisp = 0
      CALL mpi_rdwr_data (lndname, 'write', fdisp, fine_lat_map, nlon, land_chr1)
      deallocate (land_chr1)
#else
      iunit = 100
      inquire(iolength=length) land_chr1
      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='unknown')
      DO nrow = 1, nlat
! modifiedy by yuan, 06/02/2016
         !DO ncol = 1, nlon
         !   land_chr1(ncol) = char(landtypes_igbp(ncol,nrow))
         !ENDDO
         land_chr1(:) = char(landtypes_igbp(:,nrow))
         write(iunit,rec=nrow,err=100) land_chr1
      ENDDO
      close (iunit)
#endif
#endif

! Deallocate the allocatable variables

#if (defined usempi)
      deallocate (buff_fine_lat_map%bstrd)
      deallocate (buff_fine_lat_map%bdisp)
#endif

      deallocate ( ntmp )
      deallocate ( elevation      ,&
#if (defined LULC_USGS)
                   landtypes_usgs ,&
#endif
#if (defined LULC_IGBP)
                   landtypes_igbp ,&
#endif
                   lakewetland    ,&
                   glacier        ,&
                   urban           )

      go to 1000
100   print 101,nrow,lndname
101   format(' record =',i8,',  error occured on file: ',a50)
1000  continue

IF (p_master) print*,'------ END rd_land_types ------'

END SUBROUTINE rd_land_types



!-----------------------------------------------------------------------
integer FUNCTION ia(chr,n,ispval)

!  purpose: to convert a n-bytes character (chr) to integer ia.
!        ** the integer data file is saved as a n-byte character
!           data file. this function is used to recover the
!           character data to the integer data.
!
!  n      --- the number of bytes in chr
!  ispval --- default value for the negative integer.

      character*(*) chr
      integer bit_1, bit_2

      bit_1 = '200'O     ! BINARY '10000000'
      bit_2 = '377'O     ! BINARY '11111111'
      ia    = 0

      ii1 = ichar(chr(1:1))
! .. get the sign -- isn=0 positive, isn=1 negative:
      jj  = iand(ii1,bit_1)
      isn = ishft(jj,-7)

! .. for negative number:
!    because the negative integers are represented by the supplementary
!    binary code inside machine.

        IF (isn.eq.1) THEN
          DO m = n+1,4
             nbit = (m-1)*8
             jj = ishft(bit_2,nbit)
             ia = ieor(jj,ia)
          ENDDO
        ENDIF

!   .. get the byte from chr:
         DO m = 1,n
           ii2 = ichar(chr(m:m))
           IF (ii2.lt.0) ii2 = ii2 + 256
           mshft = (n-m)*8
           ia2   = ishft(ii2,mshft)
!   .. the abs(integer):
           ia = ieor(ia,ia2)
         ENDDO

      IF (ia.lt.0) ia = ispval

      RETURN
END FUNCTION ia



!-----------------------------------------------------------------------
SUBROUTINE substitute(i0,j0,ni,njl,nju,nl,np,a,c,L)
!-----------------------------------------------------------------------
!
! Initial Author : Ji Duoying, 02/2014
!-----------------------------------------------------------------------
IMPLICIT NONE

      integer, intent(in)  :: i0
      integer, intent(in)  :: j0
      integer, intent(in)  :: ni
      integer, intent(in)  :: njl, nju
      integer, intent(in)  :: nl
      integer, intent(in)  :: np
      integer, intent(in)  :: a(ni,njl:nju)
      integer, intent(in)  :: c(np)
      integer, intent(out) :: L

      integer i1, j1, i2, j2, r, loc1(1), loc2(1)
      integer num1, num2, num(nl)
      logical msk(nl)

      integer, parameter :: rmax = 100  ! searching radius (grid cells)
! --------------------------------------------------------------------
      num = 0

      L = a(i0,j0)
      DO r = 1, rmax
         DO j1 = j0-r, j0+r
            DO i1 = i0-r, i0+r
               IF(abs(j1-j0)<r .and. abs(i1-i0)<r) CYCLE
               IF(j1<njl .or. j1>nju) CYCLE
               i2 = i1
               j2 = j1
               IF(i1<1)  i2 = i1+ni
               IF(i1>ni) i2 = i1-ni
               IF(a(i2,j2) > 0)THEN
                  num(a(i2,j2)) = num(a(i2,j2)) + 1
               ENDIF
            ENDDO
         ENDDO

         msk = .true.
         msk(c(1:np)) = .false.

         loc1 = maxloc(num,mask=msk)
         num1 = num(loc1(1))

         msk(loc1(1)) = .false.

         loc2 = maxloc(num,mask=msk)
         num2 = num(loc2(1))

         IF(num1 /= num2) THEN
            L = loc1(1)
            EXIT
         ENDIF
      ENDDO

!     IF(r.eq.(rmax+1))THEN
!        print*, 'failed to find suitable type in 100 grid cells' radius'
!     ENDIF

END SUBROUTINE substitute


#if (defined usempi)
SUBROUTINE update_buff (buff, nlon, nlat, buff_lb, nrow_start, nrow_end, buff_ub)

    USE spmd_TM

    IMPLICIT NONE
    integer, intent(in)    :: nlon
    integer, intent(in)    :: nlat

    integer, intent(in)    :: buff_lb
    integer, intent(in)    :: buff_ub
    integer, intent(in)    :: nrow_start
    integer, intent(in)    :: nrow_end

    integer, intent(inout) :: buff(nlon,buff_lb:buff_ub)

    integer,   allocatable :: nlat_proc(:)
    integer,   allocatable :: lat_disp (:)
    character, allocatable :: cbuff  (:,:)
    character, allocatable :: cbuff_g(:,:)

    integer, parameter :: rmax = 100  ! searching radius (grid cells)
    integer :: iproc
    integer, allocatable :: reqs(:)
    integer :: blb, bub


    allocate (cbuff (nlon,buff_lb:buff_ub))
    cbuff = char(buff)

    IF (p_iam_slave == 0) THEN
        allocate (reqs     (0:p_nslaves-1))
        allocate (nlat_proc(0:p_nslaves-1))
        allocate (lat_disp (0:p_nslaves-1))
        CALL mpi_gather (nrow_end-nrow_start+1, 1, MPI_INTEGER, nlat_proc, 1, MPI_INTEGER, &
                         0, p_comm_slave, p_err)

        lat_disp(0) = 0
        DO iproc = 1, p_nslaves-1
            lat_disp(iproc) = sum(nlat_proc(0:iproc-1))
        ENDDO

        allocate (cbuff_g (nlon, nlat))
        CALL mpi_gatherv (cbuff(:,nrow_start:nrow_end), nlon*(nrow_end-nrow_start+1), MPI_CHARACTER, &
                          cbuff_g, nlon*nlat_proc, nlon*lat_disp, MPI_CHARACTER, &
                          0, p_comm_slave, p_err)

        cbuff(:,:) = cbuff_g(:,buff_lb:buff_ub)
        DO iproc = 1, p_nslaves-1
            blb = max(lat_disp(iproc)+1-rmax, 1)
            bub = min(lat_disp(iproc)+nlat_proc(iproc)+rmax, nlat)
            CALL mpi_isend (cbuff_g(:,blb:bub), nlon*(bub-blb+1), MPI_CHARACTER, &
                            iproc, iproc, p_comm_slave, reqs(iproc), p_err)
        ENDDO
        CALL mpi_waitall (p_nslaves-1, reqs(1:p_nslaves-1), MPI_STATUSES_IGNORE, p_err)
    ELSE
        CALL mpi_gather  (nrow_end-nrow_start+1, 1, MPI_INTEGER, 0, 1, MPI_INTEGER, &
                          0, p_comm_slave, p_err)
        CALL mpi_gatherv (cbuff(:,nrow_start:nrow_end), nlon*(nrow_end-nrow_start+1), MPI_CHARACTER, &
                          0, 0, 0, MPI_CHARACTER, 0, p_comm_slave, p_err)
        CALL mpi_recv (cbuff, nlon*(buff_ub-buff_lb+1), MPI_CHARACTER, 0, p_iam_slave, &
                       p_comm_slave, p_stat, p_err)
    ENDIF

    buff = ichar(cbuff)

    deallocate (cbuff)
    IF (p_iam_slave == 0) THEN
        deallocate (reqs)
        deallocate (nlat_proc)
        deallocate (lat_disp )
        deallocate (cbuff_g)
    ENDIF

END SUBROUTINE update_buff
#endif
!-----------------------------------------------------------------------
!EOP
