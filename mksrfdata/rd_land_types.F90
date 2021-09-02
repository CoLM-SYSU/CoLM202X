#include <define.h>
SUBROUTINE rd_land_types(dir_rawdata)
! ----------------------------------------------------------------------
! => Read in land cover dataset from original "raw" data files -
!     data with 30 arc seconds resolution
! => Fill the missing data 
! => Correct and update the land types with the specific datasets
!
! 1. Global Elevation Dataset (GTOPO30)
!    (http://webgis.wr.usgs.gov/globalgis/gtopo30/)
!    The elevation values range from -407 to 8,752 meters.
!    ocean areas have been assigned a value of -9999. 
!
! 2. Global Land Cover Characteristics
! 2.1 GLCC USGS Land cover /land use types
!     (http://edc2.usgs.gov/glcc/)
!	Value	Description
!       0       Ocean
! 	1 	Urban and Built-Up Land
! 	2 	Dryland Cropland and Pasture
! 	3 	Irrigated Cropland and Pasture
! 	4 	Mixed Dryland/Irrigated Cropland and Pasture
! 	5 	Cropland/Grassland Mosaic
! 	6 	Cropland/Woodland Mosaic
! 	7 	Grassland
! 	8 	Shrubland
! 	9 	Mixed Shrubland/Grassland
! 	10 	Savanna
! 	11 	Deciduous Broadleaf Forest
! 	12 	Deciduous Needleleaf Forest
! 	13 	Evergreen Broadleaf Forest
! 	14 	Evergreen Needleleaf Forest
! 	15 	Mixed Forest
! 	16 	Land Water Bodies
! 	17 	Herbaceous Wetland
! 	18 	Wooded Wetland
! 	19 	Barren or Sparsely Vegetated
! 	20 	Herbaceous Tundra
! 	21 	Wooded Tundra
! 	22 	Mixed Tundra
! 	23 	Bare Ground Tundra
! 	24 	Snow or Ice
!
! 2.2 MODIS Collection 5 IGBP global land cover
!	Value	Description
!       0       Ocean
!	1 	Evergreen Needleleaf Forest
!	2 	Evergreen Broadleaf Forest
!	3 	Deciduous Needleleaf Forest
!	4 	Deciduous Broadleaf Forest
!	5 	Mixed Forest
!	6 	Closed Shrublands
!	7 	Open Shrublands
!	8 	Woody Savannas
!	9 	Savannas
!	10 	Grasslands
!	11 	Permanent Wetlands
!	12 	Croplands
!	13 	Urban and Built-Up
!	14 	Cropland/Natural Vegetation Mosaic
!	15 	Snow and Ice
!	16 	Barren or Sparsely Vegetated
!	17 	Land Water Bodies
!
! 2.3 PFT classification (ON GOING PROJECT ......) 
!
! 3. Global Glacier/Ice Sheet Characteristics
!    (http://www.glims.org/RGI/; http://glims.colorado.edu/glacierdata/)
!
! 4. Global Lake and Wetlands Types
!    (http://www.wwfus.org/science/data.cfm)
!	1	Lake
!	2	Reservoir
!	3	River
!	4	Freshwater Marsh, Floodplain
!	5	Swamp Forest, Flooded Forest
!	6	Coastal Wetland (incl. Mangrove, Estuary, Delta, Lagoon)
!	7	Pan, Brackish/Saline Wetland
!	8	Bog, Fen, Mire (Peatland)
!	9	Intermittent Wetland/Lake
!	10	50-100% Wetland
!	11	25-50% Wetland
!	12	Wetland Complex (0-25% Wetland)
!
! 5. Global Urban Extent
!     (http://sage.wisc.edu/urbanenvironment.html)

! 6. Global Cropping Types (ON GOING PROJECT ......)
!
! Created by Yongjiu Dai, 12/2013
! ----------------------------------------------------------------------
use precision

IMPLICIT NONE

! arguments:
      character(len=256), intent(in) :: dir_rawdata 

! local variables:
      integer, parameter :: nlat=21600    ! 180*(60*2)
      integer, parameter :: nlon=43200    ! 360*(60*2)

      character(len=256) lndname

      character(len=1) land_chr1(nlon)
      character(len=2) land_chr2(nlon)
      integer(kind=1)  land_int1(nlon)
      integer(kind=2)  land_int2(nlon)

      ! (1) global topography
      ! -----------------
      integer, allocatable :: elevation(:,:)

      ! (2) global land cover characteristics
      ! ---------------------------------
#if(defined USGS_CLASSIFICATION)
      integer, allocatable :: landtypes_usgs(:,:)  ! GLCC USGS land cover types 
#endif
#if(defined IGBP_CLASSIFICATION)
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

!   ---------------------------------------------------------------
!   
      allocate ( elevation      (nlon,nlat) ,&
#if(defined USGS_CLASSIFICATION)
                 landtypes_usgs (nlon,nlat) ,& 
#endif
#if(defined IGBP_CLASSIFICATION)
                 landtypes_igbp (nlon,nlat) ,& 
#endif

                 lakewetland    (nlon,nlat) ,&
                 glacier        (nlon,nlat) ,& 
                 urban          (nlon,nlat)  ) 

! ----------------------------------------------------------------------
! ... (1) global digital elevation model (DEM) data
! ................................................
!     (1.1) elevation height (median)
      iunit = 100
      inquire(iolength=length) land_chr2

      lndname = trim(dir_rawdata)//'terrain/dem-usgs.bin'
      print*,lndname
      ii=0

      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old') 
      do nrow = 1, nlat 
         read(iunit,rec=nrow,err=100) land_chr2

! added by yuan, 06/02/2016
#ifdef OPENMP
!print *, 'OPENMP enabled, threads num = ', OPENMP, "dem_usgs..."
!$OMP PARALLEL DO NUM_THREADS(OPENMP) SCHEDULE(DYNAMIC,1) &
!$OMP REDUCTION(+:ii)
#endif
         do ncol = 1, nlon
            elevation(ncol,nrow) = ia(land_chr2(ncol),2,-9999)
            if( elevation(ncol,nrow) < -9990)then
               ii = ii + 1
            endif
         enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      enddo 
      close (iunit)

      print*,minval(elevation), maxval(elevation)
      print*,'ii=', ii

! ........................................
! ... (2) gloabl land cover characteristics  
! ........................................
!     (2.1) global land cover type (version 2.0) (USGS) 
      iunit = 100 
      inquire(iolength=length) land_chr1 

#if(defined USGS_CLASSIFICATION)
     ! GLCC USGS classification
     ! -------------------
      ii = 0
      iii = 0
      iiii = 0
      jjj = 0

      lndname = trim(dir_rawdata)//'landtypes/landtypes-usgs.bin' 
      print*,lndname
      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old') 
      do nrow = 1, nlat 
         read(iunit,rec=nrow,err=100) land_chr1 

! added by yuan, 06/02/2016
#ifdef OPENMP
!print *, 'OPENMP enabled, threads num = ', OPENMP, "landtypes-usgs..."
!$OMP PARALLEL DO NUM_THREADS(OPENMP) SCHEDULE(DYNAMIC,1) &
!$OMP REDUCTION(+:ii,iii,iiii,jjj)
#endif
         do ncol = 1, nlon 
            landtypes_usgs(ncol,nrow) = ichar(land_chr1(ncol)) 
            if(elevation(ncol,nrow) > -9990) then
               if(landtypes_usgs(ncol,nrow) == 16) then
                  ii = ii + 1
               endif
               if(landtypes_usgs(ncol,nrow) == 17 .or. landtypes_usgs(ncol,nrow) == 18)then
                  landtypes_usgs(ncol,nrow) = 17  ! combine herbaceous / wooded wetland => permanent wetland
                  iii = iii + 1
               endif
               if(landtypes_usgs(ncol,nrow) == 24) then
                  iiii = iiii + 1
               endif
               if(landtypes_usgs(ncol,nrow) == 1) then
                  jjj = jjj + 1
               endif
            endif
            if(elevation(ncol,nrow) < -9990) then
               landtypes_usgs(ncol,nrow) = 0
            endif
         enddo 
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      enddo 
      close (iunit)

      print*,minval(landtypes_usgs), maxval(landtypes_usgs)
      print*,' USGS GLCC land cover '
      print*,'land water points =', ii, 'wetland points=', iii, 'glacier points=', iiii, 'urban points=', jjj
#endif

#if(defined IGBP_CLASSIFICATION)
     ! MODIS IGBP classification
     ! -------------------
      ii = 0
      iii = 0
      iiii = 0
      jjj = 0

      lndname = trim(dir_rawdata)//'landtypes/landtypes-modis-igbp-2012.bin'
      print*,lndname
      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
      do nrow = 1, nlat
         read(iunit,rec=nrow,err=100) land_chr1

! added by yuan, 06/02/2016
#ifdef OPENMP
!print *, 'OPENMP enabled, threads num = ', OPENMP, "landtypes-modis-igbp..."
!$OMP PARALLEL DO NUM_THREADS(OPENMP) SCHEDULE(DYNAMIC,1) &
!$OMP REDUCTION(+:ii,iii,iiii,jjj)
#endif
         do ncol = 1, nlon
            landtypes_igbp(ncol,nrow) = ichar(land_chr1(ncol))
            if(elevation(ncol,nrow) > -9990) then
               if(landtypes_igbp(ncol,nrow) == 0) then
                 landtypes_igbp(ncol,nrow) = 17
                 ii = ii + 1
               endif
               if(landtypes_igbp(ncol,nrow) == 11) then
                 iii = iii + 1
               endif
               if(landtypes_igbp(ncol,nrow) == 15) then
                 iiii = iiii + 1
               endif
               if(landtypes_igbp(ncol,nrow) == 13) then
                 jjj = jjj + 1
               endif
            endif
            if(elevation(ncol,nrow) < -9990) then
               landtypes_igbp(ncol,nrow) = 0
            endif
         enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      enddo
      close (iunit)

      print*,minval(landtypes_igbp), maxval(landtypes_igbp)
      print*,' MODIS IGBP land cover '
      print*,'land water points =', ii, 'wetland points=', iii, 'glacier points=', iiii, 'urban points=', jjj
#endif 

! ................................................
! ... (3) global lakes and wetland characterristics
! ................................................
!     (3.1) global lakes and wetland 
      iunit = 100
      inquire(iolength=length) land_chr1

      lndname = trim(dir_rawdata)//'lake_wetland/glwd.bin' 
      print*,lndname
      ii = 0
      iii = 0

      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old') 
      do nrow = 1, nlat 
         read(iunit,rec=nrow,err=100) land_chr1

! added by yuan, 06/02/2016
#ifdef OPENMP
!print *, 'OPENMP enabled, threads num = ', OPENMP, "lakes and wetland..."
!$OMP PARALLEL DO NUM_THREADS(OPENMP) SCHEDULE(DYNAMIC,1) &
!$OMP REDUCTION(+:ii,iii)
#endif
         do ncol = 1, nlon
            lakewetland(ncol,nrow) = ichar(land_chr1(ncol))
            if(elevation(ncol,nrow)>-9990) then
               ! Replace GLCC_USGS and MODIS_IGBP water bodies with GLWD classification
               if(lakewetland(ncol,nrow) <= 3) then
#if(defined USGS_CLASSIFICATION)
                  landtypes_usgs(ncol,nrow) = 16
#endif
#if(defined IGBP_CLASSIFICATION)
                  landtypes_igbp(ncol,nrow) = 17
#endif
                  ii = ii + 1
               endif
               ! Replace GLCC_USGS and MODIS_IGBP wetland with GLWD classification
! modified by dai, 06/02/2016
               !if(lakewetland(ncol,nrow) > 3 .and. lakewetland(ncol,nrow) < 13)then
! deleted by dai, 08/15/2016
!               if(lakewetland(ncol,nrow) > 3 .and. lakewetland(ncol,nrow) < 11)then
!#if(defined USGS_CLASSIFICATION)
!                  landtypes_usgs(ncol,nrow) = 17
!#endif
!#if(defined IGBP_CLASSIFICATION)
!                  landtypes_igbp(ncol,nrow) = 11
!#endif
!                  iii = iii + 1
!               endif
               if(nrow > 18000)then  ! (90N + 60S) * 120  
                  lakewetland(ncol,nrow) = 255
               endif
            endif
         enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      enddo 
      close (iunit)
      print*,minval(lakewetland), maxval(lakewetland)
      print*,'GLWD'
      print*,'land water (1-3) =', ii, 'wetland (4-12) =',iii

! ......................................
! ... (4) global glacier and ice sheet characterristics
! ......................................
      iunit = 100
      inquire(iolength=length) land_int2

      lndname = trim(dir_rawdata)//'glacier/glacier.bin' 
      print*,lndname
      ii = 0

      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old') 
      do nrow = 1, nlat 
         read(iunit,rec=nrow,err=100) land_int2

! added by yuan, 06/02/2016
#ifdef OPENMP
!print *, 'OPENMP enabled, threads num = ', OPENMP, "glacier..."
!$OMP PARALLEL DO NUM_THREADS(OPENMP) SCHEDULE(DYNAMIC,1) &
!$OMP REDUCTION(+:ii)
#endif
         do ncol = 1, nlon
            glacier(ncol,nrow) = land_int2(ncol) * 0.1
            if(elevation(ncol,nrow)>-9990) then
               ! Replace GLCC and MODIS glacier/ice sheet with glacier specific dataset
               if(glacier(ncol,nrow) > 0.5 .and. glacier(ncol,nrow) < 101.)then
#if(defined USGS_CLASSIFICATION)
                  landtypes_usgs(ncol,nrow) = 24
#endif
#if(defined IGBP_CLASSIFICATION)
                  landtypes_igbp(ncol,nrow) = 15
#endif
                  ii = ii + 1
               endif
               ! Antarctic (ice sheet / baren ONLY)
               if(nrow > 18000)then  ! (90N + 60S) * 120  
#if(defined USGS_CLASSIFICATION)
                  if(landtypes_usgs(ncol,nrow)/=23)then
                     landtypes_usgs(ncol,nrow) = 24
                     glacier(ncol,nrow) = 100.
                  endif
#endif
#if(defined IGBP_CLASSIFICATION)
                  if(landtypes_igbp(ncol,nrow)/=16)then
                     landtypes_igbp(ncol,nrow) = 15
                     glacier(ncol,nrow) = 100.
                  endif
#endif
               endif
! modified by dai, 06/02/2016
!               ! Greenland  (ice sheet / baren / built-up ONLY)
!               ! BETWEEN (59-83N, 74-11W) =>| 0<nrow<(90-31)*120, (180-74)*120<ncol<(180-11)*120
!               if(nrow<3720 .and. (ncol>12720 .and. ncol<20280))then
!#if(defined USGS_CLASSIFICATION)
!                  if(landtypes_usgs(ncol,nrow)/=23 .or. landtypes_usgs(ncol,nrow)/=1)then
!                     landtypes_usgs(ncol,nrow) = 24  
!                     glacier(ncol,nrow) = 100.
!                  endif
!#endif
!#if(defined IGBP_CLASSIFICATION)
!                  if(landtypes_igbp(ncol,nrow)/=16 .or. landtypes_igbp(ncol,nrow)/=13)then
!                     landtypes_igbp(ncol,nrow) = 15
!                     glacier(ncol,nrow) = 100.
!                  endif
!#endif
!               endif
            endif
         enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      enddo 
      close (iunit)
      print*, minval(glacier), maxval(glacier)
      print*,'DATA glacier points'
      print*,'glacier (> 0.5%) =', ii

! ....................................
! ... (5) global urban characterristics
! ....................................
      iunit = 100
      inquire(iolength=length) land_chr1

      lndname = trim(dir_rawdata)//'urban/urban-builtup.bin'
      print*,lndname
      ii = 0

      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
      do nrow = 1, nlat
         read(iunit,rec=nrow,err=100) land_chr1

! added by yuan, 06/02/2016
#ifdef OPENMP
!print *, 'OPENMP enabled, threads num = ', OPENMP, "urban-builtup..."
!$OMP PARALLEL DO NUM_THREADS(OPENMP) SCHEDULE(DYNAMIC,1) &
!$OMP REDUCTION(+:ii)
#endif
         do ncol = 1, nlon
            urban(ncol,nrow) = ichar(land_chr1(ncol)) 
            if(elevation(ncol,nrow)>-9990) then
               ! Replace GLCC and MODIS urban/built-up with urban specific dataset
               if(urban(ncol,nrow) == 1)then
#if(defined USGS_CLASSIFICATION)
                  landtypes_usgs(ncol,nrow) = 1
#endif
#if(defined IGBP_CLASSIFICATION)
                  landtypes_igbp(ncol,nrow) = 13
#endif
                  ii = ii + 1
               endif
               if(nrow > 18000)then  ! (90N + 60S) * 120  
                  urban(ncol,nrow) = 0
               endif
            endif
         enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      enddo
      close (iunit)
      print*, minval(urban), maxval(urban)
      print*,'URBAN'
      print*,'Urban and Built-up Land Points =', ii

#if(defined CLMDEBUG)
#if(defined USGS_CLASSIFICATION)
      ii = 0
      iii = 0
      iiii = 0
      jj = 0
      jjj = 0
      jjjj = 0

      do nrow = 1, nlat
         do ncol = 1, nlon
            if(elevation(ncol,nrow) > -9990) then
               if(landtypes_usgs(ncol,nrow) == 16) then
                  ii = ii + 1
               endif
               if(landtypes_usgs(ncol,nrow) == 17 .or. landtypes_usgs(ncol,nrow) == 18)then
                  iii = iii + 1
               endif
               if(landtypes_usgs(ncol,nrow) == 24) then
                  iiii = iiii + 1
               endif
               if(landtypes_usgs(ncol,nrow) == 1) then
                  jjj = jjj + 1
               endif
            endif
         enddo
      enddo
      print*,' SECOND USGS GLCC land cover '
      print*,'land water points =', ii, 'wetland points=', iii, 'glacier points=', iiii, 'urban points=', jjj
#endif

#if(defined IGBP_CLASSIFICATION)
      ii = 0
      iii = 0
      iiii = 0
      jj = 0
      jjj = 0
      jjjj = 0

      do nrow = 1, nlat
         do ncol = 1, nlon
            if(elevation(ncol,nrow) > -9990) then
               if(landtypes_igbp(ncol,nrow) == 17) then
                 ii = ii + 1
               endif
               if(landtypes_igbp(ncol,nrow) == 11) then
                 iii = iii + 1
               endif
               if(landtypes_igbp(ncol,nrow) == 15) then
                 iiii = iiii + 1
               endif
               if(landtypes_igbp(ncol,nrow) == 13) then
                 jjj = jjj + 1
               endif
            endif
         enddo
      enddo
      print*,' SECOND MODIS IGBP land cover '
      print*,'land water points =', ii, 'wetland points=', iii, 'glacier points=', iiii, 'urban points=', jjj
#endif
#endif

! deleted by dai, 08/15/2016
!---------------------------------------------------------------------------------------------
! Correct the land cover types on which the types (GLCC_USGS/MODIS_IGBP) were classified 
! water bodies/wetland/glacier/urban, but NOT classified in GLWD/GLACIER/URBAN specific dataset
!---------------------------------------------------------------------------------------------
!      allocate (ntmp(nlon,nlat))
!      ! WATER BODIES
!      ! ------------
!#if(defined USGS_CLASSIFICATION)
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
!
!! added by yuan, 06/02/2016
!#ifdef OPENMP
!print *, 'OPENMP enabled, threads num = ', OPENMP, "water bodies substitute..."
!!$OMP PARALLEL DO NUM_THREADS(OPENMP) SCHEDULE(DYNAMIC,1) &
!!$OMP PRIVATE(L) REDUCTION(+:iii)
!#endif
!      do j = 1, nlat
!         do i = 1, nlon
!            if(elevation(i,j)>-9990) then
!               if(lakewetland(i,j)>3 .and. landtypes_usgs(i,j)==16) then
!                  call substitute(i,j,nlon,nlat,nl,np,landtypes_usgs,exclude,L)
!                  ntmp(i,j) = L
!                  iii = iii + 1
!               endif
!            endif
!         enddo
!      enddo
!#ifdef OPENMP
!!$OMP END PARALLEL DO
!#endif
!
!      landtypes_usgs = ntmp
!      deallocate (exclude)
!      print*, 'GLCC WATER BODIES','iii=',iii
!#endif
!
!#if(defined IGBP_CLASSIFICATION)
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
!
!! added by yuan, 06/02/2016
!#ifdef OPENMP
!print *, 'OPENMP enabled, threads num = ', OPENMP, "water bodies substitute..."
!!$OMP PARALLEL DO NUM_THREADS(OPENMP) SCHEDULE(DYNAMIC,1) &
!!$OMP PRIVATE(L) REDUCTION(+:iii)
!#endif
!      do j = 1, nlat
!         do i = 1, nlon
!            if(elevation(i,j)>-9990) then
!               if(lakewetland(i,j)>3 .and. landtypes_igbp(i,j)==17) then
!                  call substitute(i,j,nlon,nlat,nl,np,landtypes_igbp,exclude,L)
!                  ntmp(i,j) = L
!                  iii = iii + 1
!               endif
!            endif
!         enddo
!      enddo
!#ifdef OPENMP
!!$OMP END PARALLEL DO
!#endif
!
!      landtypes_igbp = ntmp
!      deallocate (exclude)
!      print*, 'MODIS IGBP WATER BODIES','iii=',iii
!#endif
!
!#if(defined CLMDEBUG)
!#if(defined USGS_CLASSIFICATION)
!      iiii = 0
!
!      do j = 1, nlat
!         do i = 1, nlon
!            if(elevation(i,j)>-9990) then
!               if(lakewetland(i,j)>3 .and. landtypes_usgs(i,j)==16) then
!                  iiii = iiii + 1
!               endif
!            endif
!         enddo
!      enddo
!      print*, 'GLCC WATER BODIES','iiii=',iiii
!#endif
!#if(defined IGBP_CLASSIFICATION)
!      iiii = 0
!
!      do j = 1, nlat
!         do i = 1, nlon
!            if(elevation(i,j)>-9990) then
!               if(lakewetland(i,j)>3 .and. landtypes_igbp(i,j)==17) then
!                  iiii = iiii + 1
!               endif
!            endif
!         enddo
!      enddo
!      print*, 'MODIS IGBP WATER BODIES','iiii=',iiii
!#endif
!#endif
  
! deleted by dai, 07/27/2016
! The Global Lake and Wetlands Types (http://www.wwfus.org/science/data.cfm)
! may have some problems in presenting the wetland type in some regions.
! 
!      ! WETLAND
!      ! -------
!#if(defined USGS_CLASSIFICATION)
!      iii = 0
!
!      ! GLCC USGS land cover classification (WETLAND)
!      nl = 24
!      np = 5
!      allocate (exclude(5))
!      exclude = (/1,16,17,18,24/)  !/urban and built-up(1),water bodies(16),herbaceous wetland(17),wooded wetland(18),snow and ice(24)/
!      ntmp = landtypes_usgs
!
!! added by yuan, 06/02/2016
!#ifdef OPENMP
!print *, 'OPENMP enabled, threads num = ', OPENMP, "wetland substitute..."
!!$OMP PARALLEL DO NUM_THREADS(OPENMP) SCHEDULE(DYNAMIC,1) &
!!$OMP PRIVATE(L) REDUCTION(+:iii)
!#endif
!      do j = 1, nlat
!         do i = 1, nlon
!            if(elevation(i,j)>-9990) then
!               if((lakewetland(i,j)<=3 .or. lakewetland(i,j)>12) .and. landtypes_usgs(i,j)==17) then
!                  call substitute(i,j,nlon,nlat,nl,np,landtypes_usgs,exclude,L)
!                  ntmp(i,j) = L
!                  iii = iii + 1
!               endif
!            endif
!         enddo
!      enddo
!#ifdef OPENMP
!!$OMP END PARALLEL DO
!#endif
!
!      landtypes_usgs = ntmp
!      deallocate (exclude)
!      print*, 'GLCC WETLAND','iii=',iii
!#endif
!
!#if(defined IGBP_CLASSIFICATION)
!      iii = 0
!
!      ! MODIS land cover classification (WETLAND)
!      nl = 17
!      np = 4
!      allocate (exclude(4))
!      exclude = (/11,13,15,17/)   !/permanent wetland(11),urban and built-up(13),snow and ice(15),water bodies(17)/
!      ntmp = landtypes_igbp
!
!! added by yuan, 06/02/2016
!#ifdef OPENMP
!print *, 'OPENMP enabled, threads num = ', OPENMP, "wetland substitute..."
!!$OMP PARALLEL DO NUM_THREADS(OPENMP) SCHEDULE(DYNAMIC,1) &
!!$OMP PRIVATE(L) REDUCTION(+:iii)
!#endif
!      do j = 1, nlat
!         do i = 1, nlon
!            if(elevation(i,j)>-9990) then
!               if((lakewetland(i,j)<=3 .or. lakewetland(i,j)>12) .and. landtypes_igbp(i,j)==11) then
!                  call substitute(i,j,nlon,nlat,nl,np,landtypes_igbp,exclude,L)
!                  ntmp(i,j) = L
!                  iii = iii + 1
!               endif
!            endif
!         enddo
!      enddo
!#ifdef OPENMP
!!$OMP END PARALLEL DO
!#endif
!
!      landtypes_igbp = ntmp
!      deallocate (exclude)
!      print*, 'MODIS IGBP WETLAND','iii=',iii
!#endif
!
!#if(defined CLMDEBUG)
!#if(defined USGS_CLASSIFICATION)
!      iiii = 0
!
!      do j = 1, nlat
!         do i = 1, nlon
!            if(elevation(i,j)>-9990) then
!               if((lakewetland(i,j)<=3 .or. lakewetland(i,j)>12) .and. landtypes_usgs(i,j)==17) then
!                  iiii = iiii + 1
!               endif
!            endif
!         enddo
!      enddo
!      print*, 'GLCC WETLAND', 'iiii=',iiii
!#endif
!#if(defined IGBP_CLASSIFICATION)
!      iiii = 0
!
!      do j = 1, nlat
!         do i = 1, nlon
!            if(elevation(i,j)>-9990) then
!               if((lakewetland(i,j)<=3 .or. lakewetland(i,j)>12) .and. landtypes_igbp(i,j)==11) then
!                  iiii = iiii + 1
!               endif
!            endif
!         enddo
!      enddo
!      print*, 'MODIS IGBP WETLAND', 'iiii=',iiii
!#endif
!#endif

      ! GLACIER/ICESHEET
      ! ------------------
#if(defined USGS_CLASSIFICATION)
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

! added by yuan, 06/02/2016
#ifdef OPENMP
print *, 'OPENMP enabled, threads num = ', OPENMP, "glacier substitute..."
!$OMP PARALLEL DO NUM_THREADS(OPENMP) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(L) REDUCTION(+:iii)
#endif
      do j = 1, nlat
         do i = 1, nlon
            if(elevation(i,j)>-9990) then
               if((glacier(i,j)<=0.5 .or. glacier(i,j)> 100.1) .and. landtypes_usgs(i,j)==24) then
                  call substitute(i,j,nlon,nlat,nl,np,landtypes_usgs,exclude,L)
                  ntmp(i,j) = L
                  iii = iii + 1
               endif
            endif
         enddo
      enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      landtypes_usgs = ntmp
      deallocate (exclude)
      print*, 'GLCC GLACIER/ICESHEET','iii=',iii
#endif
#if(defined IGBP_CLASSIFICATION)
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

! added by yuan, 06/02/2016
#ifdef OPENMP
print *, 'OPENMP enabled, threads num = ', OPENMP, "glacier substitute..."
!$OMP PARALLEL DO NUM_THREADS(OPENMP) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(L) REDUCTION(+:iii)
#endif
      do j = 1, nlat
         do i = 1, nlon
            if(elevation(i,j)>-9990) then
               if((glacier(i,j)<=0.5 .or. glacier(i,j)>100.1) .and. landtypes_igbp(i,j)==15) then
                  call substitute(i,j,nlon,nlat,nl,np,landtypes_igbp,exclude,L)
                  ntmp(i,j) = L
                  iii = iii + 1
               endif
            endif
         enddo
      enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      landtypes_igbp = ntmp
      deallocate (exclude)
      print*, 'MODIS IGBP GLACIER/ICESHEET','iii=',iii
#endif

#if(defined CLMDEBUG)
#if(defined USGS_CLASSIFICATION)
      iiii = 0

      do j = 1, nlat
         do i = 1, nlon
            if(elevation(i,j)>-9990) then
               if((glacier(i,j)<=0.5 .or. glacier(i,j)> 100.1) .and. landtypes_usgs(i,j)==24) then
                  iiii = iiii + 1
               endif
            endif
         enddo
      enddo
      print*, 'GLCC GLACIER/ICESHEET', 'iiii=',iiii
#endif
#if(defined IGBP_CLASSIFICATION)
      iiii = 0

      do j = 1, nlat
         do i = 1, nlon
            if(elevation(i,j)>-9990) then
               if((glacier(i,j)<=0.5 .or. glacier(i,j)>100.1) .and. landtypes_igbp(i,j)==15) then
                  iiii = iiii + 1
               endif
            endif
         enddo
      enddo
      print*, 'MODIS IGBP GLACIER/ICESHEET', 'iiii=',iiii
#endif
#endif

      ! URBAN and BUILT-UP LAND
      ! ------------------
#if(defined USGS_CLASSIFICATION)
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

! added by yuan, 06/02/2016
#ifdef OPENMP
print *, 'OPENMP enabled, threads num = ', OPENMP, "urban substitute..."
!$OMP PARALLEL DO NUM_THREADS(OPENMP) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(L) REDUCTION(+:iii)
#endif
      do j = 1, nlat
         do i = 1, nlon
            if(elevation(i,j)>-9990) then
               if(urban(i,j)==0 .and. landtypes_usgs(i,j)==1) then
                  call substitute(i,j,nlon,nlat,nl,np,landtypes_usgs,exclude,L)
                  ntmp(i,j) = L
                  iii = iii + 1
               endif
            endif
         enddo
      enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      landtypes_usgs = ntmp
      deallocate (exclude)
      print*, 'GLCC URBAN and BUILT-UP LAND','iii=',iii
#endif

#if(defined IGBP_CLASSIFICATION)
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

! added by yuan, 06/02/2016
#ifdef OPENMP
print *, 'OPENMP enabled, threads num = ', OPENMP, "urban substitute..."
!$OMP PARALLEL DO NUM_THREADS(OPENMP) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(L) REDUCTION(+:iii)
#endif
      do j = 1, nlat
         do i = 1, nlon
            if(elevation(i,j)>-9990) then
               if(urban(i,j)==0 .and. landtypes_igbp(i,j)==13) then
                  call substitute(i,j,nlon,nlat,nl,np,landtypes_igbp,exclude,L)
                  ntmp(i,j) = L
                  iii = iii + 1
               endif
            endif
         enddo
      enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      landtypes_igbp = ntmp
      deallocate (exclude)
      print*, 'MODIS IGBP URBAN and BUILT-UP LAND','iii=',iii
#endif

#if(defined CLMDEBUG)
#if(defined USGS_CLASSIFICATION)
      iiii = 0

      do j = 1, nlat
         do i = 1, nlon
            if(elevation(i,j)>-9990) then
               if(urban(i,j)==0 .and. landtypes_usgs(i,j)==1) then
                  iiii = iiii + 1
               endif
            endif
         enddo
      enddo
      print*, 'GLCC URBAN and BUILT-UP LAND', 'iiii=',iiii
#endif

#if(defined IGBP_CLASSIFICATION)
iiii = 0
      do j = 1, nlat
         do i = 1, nlon
            if(elevation(i,j)>-9990) then
               if(urban(i,j)==0 .and. landtypes_igbp(i,j)==13) then
                  iiii = iiii + 1
               endif
            endif
         enddo
      enddo
      print*, 'MODIS IGBP URBAN and BUILT-UP LAND', 'iiii=',iiii
#endif
#endif

      !--------------------------------------------------------------------
#if(defined CLMDEBUG)
#if(defined USGS_CLASSIFICATION)
      ii = 0
      iii = 0
      iiii = 0
      jj = 0
      jjj = 0
      jjjj = 0

      do nrow = 1, nlat
         do ncol = 1, nlon
            if(elevation(ncol,nrow) > -9990) then
               if(landtypes_usgs(ncol,nrow) == 16) then
                  ii = ii + 1
               endif
               if(landtypes_usgs(ncol,nrow) == 17 .or. landtypes_usgs(ncol,nrow) == 18)then
                  iii = iii + 1
               endif
               if(landtypes_usgs(ncol,nrow) == 24) then
                  iiii = iiii + 1
               endif
               if(landtypes_usgs(ncol,nrow) == 1) then
                  jjj = jjj + 1
               endif
            endif
         enddo
      enddo
      print*,' FINAL USGS GLCC land cover '
      print*,'land water points =', ii, 'wetland points=', iii, 'glacier points=', iiii, 'urban points=', jjj
#endif

#if(defined IGBP_CLASSIFICATION)
      ii = 0
      iii = 0
      iiii = 0
      jj = 0
      jjj = 0
      jjjj = 0

      do nrow = 1, nlat
         do ncol = 1, nlon
            if(elevation(ncol,nrow) > -9990) then
               if(landtypes_igbp(ncol,nrow) == 17) then
                 ii = ii + 1
               endif
               if(landtypes_igbp(ncol,nrow) == 11) then
                 iii = iii + 1
               endif
               if(landtypes_igbp(ncol,nrow) == 15) then
                 iiii = iiii + 1
               endif
               if(landtypes_igbp(ncol,nrow) == 13) then
                 jjj = jjj + 1
               endif
            endif
         enddo
      enddo
      print*,' FINAL MODIS IGBP land cover '
      print*,'land water points =', ii, 'wetland points=', iii, 'glacier points=', iiii, 'urban points=', jjj
#endif
#endif

!#if(defined ongoing)
! ..............................................
! ... (6) global cultural characteristics (crops)
! ..............................................
!
!#endif 

! Write out the land cover types
      iunit = 100
      inquire(iolength=length) land_chr1
#if(defined USGS_CLASSIFICATION)
     ! GLCC USGS classification
     ! -------------------
      lndname = trim(dir_rawdata)//'RAW_DATA_updated_with_igbp/landtypes_usgs_update.bin'
      print*,lndname
      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='unknown')
      do nrow = 1, nlat
! modifiedy by yuan, 06/02/2016
         !do ncol = 1, nlon
         !   land_chr1(ncol) = char(landtypes_usgs(ncol,nrow))
         !enddo
         land_chr1(:) = char(landtypes_usgs(:,nrow))
         write(iunit,rec=nrow,err=100) land_chr1
      enddo
      close (iunit)
#endif
#if(defined IGBP_CLASSIFICATION)
     ! MODIS IGBP classification
     ! -------------------
      lndname = trim(dir_rawdata)//'RAW_DATA_updated_with_igbp/landtypes_igbp_update.bin'
      print*,lndname
      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='unknown')
      do nrow = 1, nlat
! modifiedy by yuan, 06/02/2016
         !do ncol = 1, nlon
         !   land_chr1(ncol) = char(landtypes_igbp(ncol,nrow))
         !enddo
         land_chr1(:) = char(landtypes_igbp(:,nrow))
         write(iunit,rec=nrow,err=100) land_chr1
      enddo
      close (iunit)
#endif

! Deallocate the allocatable variables
      deallocate ( ntmp )
      deallocate ( elevation      ,&
#if(defined USGS_CLASSIFICATION)
                   landtypes_usgs ,&
#endif
#if(defined IGBP_CLASSIFICATION)
                   landtypes_igbp ,&
#endif
                   lakewetland    ,&
                   glacier        ,&
                   urban           )

      go to 1000
100   print 101,nrow,lndname
101   format(' record =',i8,',  error occured on file: ',a50)
1000  continue

print*,'------ END rd_land_types ------'

END SUBROUTINE rd_land_types



!-----------------------------------------------------------------------
INTEGER FUNCTION ia(chr,n,ispval)

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
                                                                                
        if (isn.eq.1) then
          do m = n+1,4
             nbit = (m-1)*8
             jj = ishft(bit_2,nbit)
             ia = ieor(jj,ia)
          end do
        endif
                                                                                
!   .. get the byte from chr:
         do m = 1,n
           ii2 = ichar(chr(m:m))
           if (ii2.lt.0) ii2 = ii2 + 256
           mshft = (n-m)*8
           ia2   = ishft(ii2,mshft)
!   .. the abs(integer):
           ia = ieor(ia,ia2)
         end do
                                                                                
      if (ia.lt.0) ia = ispval

      return
END FUNCTION ia



!-----------------------------------------------------------------------
SUBROUTINE substitute(i0,j0,ni,nj,nl,np,a,c,L)
!-----------------------------------------------------------------------
!
! Initial Author : Ji Duoying, 02/2014
!-----------------------------------------------------------------------
IMPLICIT NONE

      integer, intent(in)  :: i0
      integer, intent(in)  :: j0
      integer, intent(in)  :: ni
      integer, intent(in)  :: nj
      integer, intent(in)  :: nl
      integer, intent(in)  :: np
      integer, intent(in)  :: a(ni,nj)
      integer, intent(in)  :: c(np)
      integer, intent(out) :: L

      integer i1, j1, i2, j2, r, loc1(1), loc2(1)
      integer num1, num2, num(nl)
      logical msk(nl)

      integer, parameter :: rmax = 100  ! searching radius (grid cells)
! --------------------------------------------------------------------
      num = 0

      L = a(i0,j0)
      do r = 1, rmax
         do j1 = j0-r, j0+r
            do i1 = i0-r, i0+r
               if(abs(j1-j0)<r .and. abs(i1-i0)<r) cycle
               if(j1<1 .or. j1>nj) cycle
               i2 = i1
               j2 = j1
               if(i1<1)  i2 = i1+ni
               if(i1>ni) i2 = i1-ni
               if(a(i2,j2) > 0)then
                  num(a(i2,j2)) = num(a(i2,j2)) + 1
               endif
            end do
         end do

         msk = .true.
         msk(c(1:np)) = .false.

         loc1 = maxloc(num,mask=msk)
         num1 = num(loc1(1))

         msk(loc1(1)) = .false.

         loc2 = maxloc(num,mask=msk)
         num2 = num(loc2(1))

         if(num1 /= num2) then
            L = loc1(1)
            exit
         endif
      end do

!     if(r.eq.(rmax+1))then
!        print*, 'failed to find suitable type in 100 grid cells' radius'
!     endif

END SUBROUTINE substitute
!-----------------------------------------------------------------------
!EOP
