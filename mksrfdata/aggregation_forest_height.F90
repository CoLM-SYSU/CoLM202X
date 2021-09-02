
#include <define.h>

SUBROUTINE aggregation_forest_height ( dir_rawdata,dir_model_landdata, &
                                       lon_points,lat_points, &
                                       nrow_start,nrow_end,ncol_start,ncol_end, &
                                       nx_fine_gridcell,ny_fine_gridcell,area_fine_gridcell,&
                                       READ_row_UB,READ_row_LB,READ_col_UB,READ_col_LB )
! ----------------------------------------------------------------------
! 1. Global land cover types (updated with the specific dataset)
!
! 2. Global Forest Height
!    (http://lidarradar.jpl.nasa.gov/)
!     Simard, M., N. Pinto, J. B. Fisher, and A. Baccini, 2011: Mapping
!     forest canopy height globally with spaceborne lidar.
!     J. Geophys. Res., 116, G04021.
!
! Created by Yongjiu Dai, 02/2014
! ----------------------------------------------------------------------
use precision

IMPLICIT NONE

! arguments:
#if(defined USGS_CLASSIFICATION)
      integer, parameter :: N_land_classification = 24 ! GLCC USGS number of land cover category
#else
      integer, parameter :: N_land_classification = 17 ! MODIS IGBP number of land cover category
#endif
      integer, parameter :: nlat = 21600  ! 180*(60*2)
      integer, parameter :: nlon = 43200  ! 360*(60*2)

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
      integer(kind=2)  land_int2(nlon)

      integer iunit
      integer length
      integer i, j, L
      integer i1, i2, j1, j2
      integer nrow, ncol, ncol_mod
      integer LL, np
      integer n_fine_gridcell

      integer, allocatable :: landtypes(:,:) ! GLCC USGS / MODIS IGBP land cover types 
      integer, allocatable :: num_patches(:)
      real(r8), allocatable :: tree_height(:,:)  ! forest canopy height (m)
      real(r8), allocatable :: a_tree_height_patches(:,:)
      real(r8), allocatable :: tree_height_patches(:,:,:)

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
      print*,lndname
      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old') 
      do nrow = nrow_start, nrow_end
         read(iunit,rec=nrow,err=100) land_chr1 
         landtypes(:,nrow) = ichar(land_chr1(:)) 
      enddo 
      close (iunit)
#endif

#if(defined IGBP_CLASSIFICATION)
     ! MODIS IGBP classification
     ! -------------------
      lndname = trim(dir_rawdata)//'RAW_DATA_updated_with_igbp/landtypes_igbp_update.bin'
      print*,lndname
      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
      do nrow = nrow_start, nrow_end
         read(iunit,rec=nrow,err=100) land_chr1
         landtypes(:,nrow) = ichar(land_chr1(:))
      enddo
      close (iunit)
#endif 

#endif

! ................................................
! ... (2) global forest canopy height (m)
! ................................................
      allocate ( tree_height (nlon,nlat) )

      iunit = 100
      inquire(iolength=length) land_chr1
      lndname = trim(dir_rawdata)//'forest_height/Forest_Height.bin'
      print*,lndname

      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old') 
      do nrow = nrow_start, nrow_end
         read(iunit,rec=nrow,err=100) land_chr1
         tree_height(:,nrow) = ichar(land_chr1(:))
      enddo 
      close (iunit)
      print*, minval(tree_height(:,nrow_start:nrow_end)), maxval(tree_height(:,nrow_start:nrow_end))

!   ---------------------------------------------------------------
!   aggregate the forest canopy height from the resolution of raw data to modelling resolution
!   ---------------------------------------------------------------
      n_fine_gridcell = nx_fine_gridcell * ny_fine_gridcell
      allocate ( num_patches(0:N_land_classification) )
      allocate ( a_tree_height_patches(0:N_land_classification,1:n_fine_gridcell) )
      allocate ( tree_height_patches(0:N_land_classification,1:lon_points,1:lat_points) )

#ifdef OPENMP
print *, 'OPENMP enabled, threads num = ', OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(i,j,i1,i2,j1,j2,nrow,ncol,ncol_mod,L,LL,num_patches,np) &
!$OMP PRIVATE(a_tree_height_patches) 
#endif
      do j = 1, lat_points

#if(defined USER_GRID)
         j1 = READ_row_UB(j)  ! read upper boundary of latitude
         j2 = READ_row_LB(j)  ! read lower boundary of latitude
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

            do nrow = j1, j2            
               if(i1 > i2) i2 = i2 + nlon   ! for coarse grid crosses the dateline    
               do ncol = i1, i2
                  ncol_mod = mod(ncol,nlon)
                  if(ncol_mod == 0) ncol_mod = nlon  
                  L = landtypes(ncol_mod,nrow)
                  ! mapping forest canopy height from "raw data" resolution to modelling resolution

#if(defined USGS_CLASSIFICATION)
                  if(L/=0 .and. L/=1 .and. L/=16 .and. L/=24)then   ! NOT OCEAN(0)/URBAN and BUILT-UP(1)/WATER BODIES(16)/ICE(24)
                     num_patches(L) = num_patches(L) + 1
                     LL = num_patches(L)
                     a_tree_height_patches (L,LL) = tree_height(ncol_mod,nrow)
                 endif
#endif
#if(defined IGBP_CLASSIFICATION)
                  if(L/=0 .and. L/=13 .and. L/=17 .and. L/=15)then  ! NOT OCEAN(0)/URBAN and BUILT-UP(13)/WATER BODIES(17)/ICE(15)
                     num_patches(L) = num_patches(L) + 1
                     LL = num_patches(L)
                     a_tree_height_patches (L,LL) = tree_height(ncol_mod,nrow)
                 endif
#endif
               enddo
            enddo
            
            do L = 0, N_land_classification
#if(defined USGS_CLASSIFICATION)
               if(L/=0 .and. L/=1 .and. L/=16 .and. L/=24)then   ! NOT OCEAN(0)/URBAN and BUILT-UP(1)/WATER BODIES(16)/ICE(24)
                  np = num_patches(L)
                  if(np == 0)then
                     tree_height_patches (L,i,j) = -1.e36
                  else if(np == 1)then
                     tree_height_patches(L,i,j) = a_tree_height_patches(L,1)
                  else
                     tree_height_patches(L,i,j) = median( a_tree_height_patches(L,1:np), np )
                  endif
               else
                  tree_height_patches(L,i,j) = -1.e36
               endif

#endif
#if(defined IGBP_CLASSIFICATION)
               if(L/=0 .and. L/=13 .and. L/=17 .and. L/=15)then  ! NOT OCEAN(0)/URBAN and BUILT-UP(13)/WATER BODIES(17)/ICE(15)
                  np = num_patches(L)
                  if(np == 0)then
                     tree_height_patches (L,i,j) = -1.e36
                  else if(np == 1)then
                     tree_height_patches(L,i,j) = a_tree_height_patches(L,1)
                  else
                     tree_height_patches(L,i,j) = median( a_tree_height_patches(L,1:np), np )
                  endif
               else
                  tree_height_patches(L,i,j) = -1.e36
               endif
#endif
            enddo
         enddo
      enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

! Write-out the forest height (m) 
      lndname = trim(dir_model_landdata)//'model_forest_height.bin'
      print*,lndname
      open(iunit,file=trim(lndname),form='unformatted',status='unknown')
      write(iunit,err=100) tree_height_patches
      close(iunit)

      deallocate ( landtypes )
      deallocate ( tree_height )
      deallocate ( num_patches )
      deallocate ( a_tree_height_patches )
      deallocate ( tree_height_patches )

      go to 1000
100   print 101,nrow,lndname
101   format(' record =',i8,',  error occured on file: ',a50)
1000  continue

END SUBROUTINE aggregation_forest_height
