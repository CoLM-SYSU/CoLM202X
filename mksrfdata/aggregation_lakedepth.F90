
#include <define.h>

SUBROUTINE aggregation_lakedepth( dir_rawdata, dir_model_landdata, &
                                  lon_points,lat_points, &
                                  nrow_start,nrow_end,ncol_start,ncol_end, &
                                  nx_fine_gridcell,ny_fine_gridcell,area_fine_gridcell,&
                                  READ_row_UB,READ_row_LB,READ_col_UB,READ_col_LB )

! ----------------------------------------------------------------------
! 1. Global land cover types (updated with the specific dataset)
!
! 2. Global Lake Coverage and Lake Depth
!   (http://nwpi.krc.karelia.run/flake/)      
!    Kourzeneva, E., H. Asensio, E. Martin, and S. Faroux, 2012: Global
!    gridded dataset of lake coverage and lake depth for use in numerical
!    weather prediction and climate modelling. Tellus A, 64, 15640.
!
!    Lake depth data legend
!    Value   Description
! 0       no lake indicated in this pixel
! 1       no any information about this lake and set the default value of 10 m
! 2       no information about depth for this lake and set the default value of 10 m
! 3       have the information about lake depth in this pixel
! 4       this is the river pixel according to our map, set the default value of 3 m
!
! Created by Yongjiu Dai, 02/2014
! ----------------------------------------------------------------------
USE precision
USE GlobalVars

IMPLICIT NONE
! arguments:
      character(LEN=256), intent(in) :: dir_rawdata
      character(LEN=256), intent(in) :: dir_model_landdata

      integer, intent(in) :: lon_points ! number of model longitude grid points
      integer, intent(in) :: lat_points ! model  of model latitude grid points
      integer, intent(in) :: nrow_start
      integer, intent(in) :: nrow_end
      integer, intent(in) :: ncol_start
      integer, intent(in) :: ncol_end
      integer, intent(in) :: nx_fine_gridcell
      integer, intent(in) :: ny_fine_gridcell
      integer, intent(in) :: READ_row_UB(lat_points)  ! north boundary index for fine gird cell
      integer, intent(in) :: READ_col_UB(lon_points)  ! west boundary index for fine gird cell  
      integer, intent(in) :: READ_row_LB(lat_points)  ! south boundary index for fine gird cell
      integer, intent(in) :: READ_col_LB(lon_points)  ! east boundary index for fine gird cell

      real(r8), intent(in) :: area_fine_gridcell(nlon,nlat)  ! rwadata fine cell area (km**2)

! local variables:
! ---------------------------------------------------------------
      character(len=256) lndname
      character(len=1) land_chr1(nlon)
      character(len=2) land_chr2(nlon)
      integer(kind=1)  land_int1(nlon)
      integer(kind=2)  land_int2(nlon30s)

      CHARACTER(len=256) suffix
      integer iunit
      integer length
      integer i, j, L
      integer i1, i2, j1, j2
      integer nrow, ncol, ncol_mod, nrow_mod
      integer LL, np
      integer n_fine_gridcell

      INTEGER nrow30s_start, nrow30s_end
      INTEGER ncol30s_start, ncol30s_end

      integer , allocatable :: landtypes(:,:)
      real(r8), allocatable :: lakedepth(:,:)
      real(r8), allocatable :: a_lakedepth(:)
      real(r8), allocatable :: lakedepth_patches(:,:)
      integer , allocatable :: num_patches(:)

      real(r8), external :: median

! ........................................
! ... (1) gloabl land cover types
! ........................................
      iunit = 100
      inquire(iolength=length) land_chr1 
      allocate ( landtypes (nlon,nlat) ) 

#if(defined USE_POINT_DATA)

#if(defined USGS_CLASSIFICATION)
      landtypes(ncol_start,nrow_start) = USGS_CLASSIFICATION
#endif

#if(defined IGBP_CLASSIFICATION)
      landtypes(ncol_start,nrow_start) = IGBP_CLASSIFICATION
#endif

#else

#if(defined USGS_CLASSIFICATION)
     ! GLCC USGS classification
     ! -------------------
      lndname = trim(dir_rawdata)//'RAW_DATA_updated_with_igbp/landtypes_usgs_update.bin' 
#else
      lndname = trim(dir_rawdata)//'landtypes/landtypes-modis-igbp-2005.bin'
#endif
 
      print*,lndname
      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old') 
      do nrow = nrow_start, nrow_end
         read(iunit,rec=nrow,err=100) land_chr1 
         landtypes(:,nrow) = ichar(land_chr1(:)) 
      enddo 
      close (iunit)

#endif

#ifdef USGS_CLASSIFICATION
      ! for 1km resolution
      nrow30s_start = nrow_start
      nrow30s_end   = nrow_end
      ncol30s_start = ncol_start
      ncol30s_end   = ncol_end
      suffix        = ''
#else
      ! for 500m resolution, to get the 1km
      ! data, need to modify the index
      nrow30s_start = int((nrow_start+1)/2)
      nrow30s_end   = int((nrow_end+1)/2)
      ncol30s_start = int((ncol_start+1)/2)
      ncol30s_end   = int((ncol_end+1)/2)
      suffix        = '.igbp'
#endif


! ................................................
! ... (2) global lake coverage and lake depth
! ................................................
      iunit = 100
      inquire(iolength=length) land_int2
      lndname = trim(dir_rawdata)//'lake_depth/GlobalLakeDepth.bin' 
      print*,lndname
      allocate ( lakedepth (nlon30s,nlat30s) )

      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old') 
      do nrow = nrow30s_start, nrow30s_end
         read(iunit,rec=nrow,err=100) land_int2
         lakedepth(:,nrow) = land_int2(:) * 0.1
      enddo 
      close (iunit)
      print*,minval(lakedepth(:,nrow30s_start:nrow30s_end)), &
             maxval(lakedepth(:,nrow30s_start:nrow30s_end))

!   ---------------------------------------------------------------
!   aggregate the lake depth from the resolution of raw data to modelling resolution
!   ---------------------------------------------------------------

      n_fine_gridcell = nx_fine_gridcell * ny_fine_gridcell
      allocate ( num_patches(0:N_land_classification) )
      allocate ( a_lakedepth(1:n_fine_gridcell) )
      allocate ( lakedepth_patches(1:lon_points,1:lat_points) )

#ifdef OPENMP
print *, 'OPENMP enabled, threads num = ', OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(i,j,i1,i2,j1,j2,nrow,ncol,ncol_mod,nrow_mod,L,LL,num_patches,np) &
!$OMP PRIVATE(a_lakedepth) 
#endif
      do j = 1, lat_points

#if(defined USER_GRID)
         j1 = READ_row_UB(j)    ! read upper boundary of latitude
         j2 = READ_row_LB(j)    ! read lower boundary of latitude
#else
         j1 = nrow_start + (j-1)*ny_fine_gridcell
         j2 = nrow_start - 1 + j*ny_fine_gridcell
#endif

         do i = 1, lon_points

#if(defined USER_GRID)
            i1 = READ_col_UB(i)   ! read upper boundary of longitude
            i2 = READ_col_LB(i)   ! read lower boundary of longitude
#else            
            i1 = ncol_start + (i-1)*nx_fine_gridcell 
            i2 = ncol_start -1 + i*nx_fine_gridcell
#endif
            num_patches(:) = 0
            a_lakedepth (:) = 0.

            do nrow = j1, j2          
               if(i1 > i2) i2 = i2 + nlon   ! for coarse grid crosses the dateline
               do ncol = i1, i2
                  ncol_mod = mod(ncol,nlon)   
                  nrow_mod = nrow
                  if(ncol_mod == 0) ncol_mod = nlon

                  ! mapping land cover type from "raw data" resolution to modelling resolution
                  L = landtypes(ncol_mod,nrow)

#if(defined USGS_CLASSIFICATION)
                  if(L==16)then  ! LAND WATER BODIES (16)
#else
                  if(L==17)then  ! LAND WATER BODIES (17)
#endif
                     num_patches(L) = num_patches(L) + 1
                     LL = num_patches(L)

! yuan, 07/30/2019: 500m (15s) ==> 1km (30s)
#ifndef USGS_CLASSIFICATION
                     ncol_mod = int((ncol_mod+1)/2)
                     nrow_mod = int((nrow+1)/2)
#endif
                     a_lakedepth (LL) = lakedepth(ncol_mod,nrow_mod)
                  endif
               enddo
            enddo

#if(defined USGS_CLASSIFICATION)
            np = num_patches(16)  ! GLCC USGS LAND WATER BODIES (16)
#else
            np = num_patches(17)  ! MODIS IGBP LAND WATER BODIES (17)
#endif
            if(np == 0)then
               lakedepth_patches(i,j) = -1.0e36
            else if(np == 1)then
               lakedepth_patches(i,j) = a_lakedepth(1)
            else
               lakedepth_patches(i,j) = median ( a_lakedepth(1:np), np)
            endif

         enddo
      enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

! Write-out the lake depth of the lake pacth in the gridcell
      lndname = trim(dir_model_landdata)//'model_GlobalLakeDepth'//trim(suffix)//'.bin' 
      print*,lndname
      open(iunit,file=trim(lndname),form='unformatted',status='unknown')
      write(iunit,err=100) lakedepth_patches
      close(iunit)

      deallocate ( landtypes )
      deallocate ( lakedepth )
      deallocate ( num_patches )
      deallocate ( a_lakedepth )
      deallocate ( lakedepth_patches )

      go to 1000
100   print 101,nrow,lndname
101   format(' record =',i8,',  error occured on file: ',a50)
1000  continue

END SUBROUTINE aggregation_lakedepth
