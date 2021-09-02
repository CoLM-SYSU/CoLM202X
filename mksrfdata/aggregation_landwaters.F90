#include <define.h>
SUBROUTINE aggregation_landwaters( dir_rawdata,dir_model_landdata, &
                                lon_points,lat_points, &
                                nrow_start,nrow_end,ncol_start,ncol_end, &
                                nx_fine_gridcell,ny_fine_gridcell,area_fine_gridcell,&
                                sinn,sins,lonw_rad,lone_rad,sinn_i,sins_i,lonw_rad_i,lone_rad_i,&
                                READ_row_UB,READ_row_LB,READ_col_UB,READ_col_LB) 
! ----------------------------------------------------------------------
! 1. Global land cover types (updated with the specific dataset)
!
! 2. Global Lake and Wetlands Types
!     (http://www.wwfus.org/science/data.cfm)
!       1       Lake
!       2       Reservoir
!       3       River
!       4       Freshwater Marsh, Floodplain
!       5       Swamp Forest, Flooded Forest
!       6       Coastal Wetland (incl. Mangrove, Estuary, Delta, Lagoon)
!       7       Pan, Brackish/Saline Wetland
!       8       Bog, Fen, Mire (Peatland)
!       9       Intermittent Wetland/Lake
!       10      50-100% Wetland
!       11      25-50% Wetland
!       12      Wetland Complex (0-25% Wetland)
!
! Created by Yongjiu Dai, 02/2014
! ________________
! REVISION HISTORY:
!   /07/2014, Siguang Zhu & Xiangxiang Zhang: weight average considering 
!               partial overlap between fine grid and model grid for a user
!               defined domain file.
!
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
      
      real(r8), intent(in) :: sinn(lat_points)        ! grid cell latitude, northern edge(sin)  
      real(r8), intent(in) :: sins(lat_points)        ! grid cell latitude, northern edge(sin)
      real(r8), intent(in) :: lonw_rad(lon_points)    ! grid cell longitude, western edge (radian)
      real(r8), intent(in) :: lone_rad(lon_points)    ! grid cell longitude, eastern edge (radian)
      real(r8), intent(in) :: sinn_i(nlat)            ! fine grid cell latitude, northern edge(sin)
      real(r8), intent(in) :: sins_i(nlat)            ! fine grid cell latitude, northern edge(sin)
      real(r8), intent(in) :: lonw_rad_i(nlon)        ! fine grid cell longitude, western edge (radian)
      real(r8), intent(in) :: lone_rad_i(nlon)        ! fine grid cell longitude, eastern edge (radian)
      integer,  intent(in) :: READ_row_UB(lat_points) ! north boundary index for fine gird cell
      integer,  intent(in) :: READ_col_UB(lon_points) ! west boundary index for fine gird cell  
      integer,  intent(in) :: READ_row_LB(lat_points) ! south boundary index for fine gird cell
      integer,  intent(in) :: READ_col_LB(lon_points) ! east boundary index for fine gird cell
      
      real(r8), parameter ::pi = 4.*atan(1.)

      real(r8), intent(in) :: area_fine_gridcell(nlon,nlat)  ! rwadata fine cell area (km**2)

! local variables:
! ----------------------------------------------------------------------
      character(len=256) lndname
      character(len=1) land_chr1(nlon)
      character(len=2) land_chr2(nlon)
      integer(kind=1)  land_int1(nlon)
      integer(kind=2)  land_int2(nlon)

      integer iunit
      integer length
      integer i, j, L, i1, i2, j1, j2
      integer nrow, ncol, ncol_mod
      integer LL, np, n, nn, nnn
      integer n_fine_gridcell

      integer, allocatable :: landtypes(:,:) ! GLCC USGS / MODIS IGBP land cover types 
      integer, allocatable :: lakewetland(:,:)  ! lake and wetland types
      integer, allocatable :: num_patches(:)
      integer, allocatable :: n_landwaters_patches(:)
      integer, allocatable :: num_landwaters(:)
      real(r8), allocatable :: f_landwaters(:)
      real(r8), allocatable :: area_landwaters_patches(:)
      real(r8), allocatable :: fraction_landwaters_patches(:,:,:)
      real(r8) area_landwaters_grids
      real(r8) area_for_sum  

      real(r8), external :: find_min_area
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
! ... (2) global lakes and wetland
! ................................................
      iunit = 100
      inquire(iolength=length) land_chr1
      lndname = trim(dir_rawdata)//'lake_wetland/glwd.bin'
      print*,lndname
      allocate ( lakewetland (nlon,nlat) )

      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old') 
      do nrow = nrow_start, nrow_end
         read(iunit,rec=nrow,err=100) land_chr1
         lakewetland(:,nrow) = ichar(land_chr1(:))
      enddo 
      close (iunit)
      print*,minval(lakewetland(:,nrow_start:nrow_end)), maxval(lakewetland(:,nrow_start:nrow_end))

!   ---------------------------------------------------------------
!   aggregate the wetland from the resolution of raw data to modelling resolution
!   ---------------------------------------------------------------
      n_fine_gridcell = nx_fine_gridcell * ny_fine_gridcell
      allocate ( num_patches(0:N_land_classification) )
      allocate ( fraction_landwaters_patches(1:3,1:lon_points,1:lat_points) )
      allocate ( n_landwaters_patches(n_fine_gridcell) )
      allocate ( num_landwaters(1:3) )
      allocate ( f_landwaters(1:3) )
      allocate ( area_landwaters_patches(n_fine_gridcell) ) 

#ifdef OPENMP
print *, 'OPENMP enabled, threads num = ', OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(i,j,i1,i2,j1,j2,nrow,ncol,ncol_mod,L,LL,num_patches,np,n,nn,nnn) &
!$OMP PRIVATE(n_landwaters_patches,area_landwaters_patches,area_landwaters_grids) &
!$OMP PRIVATE(num_landwaters,f_landwaters,area_for_sum) 
#endif
      do j = 1, lat_points

#if(defined USER_GRID)
         j1 = READ_row_UB(j)   ! read upper boundary of latitude 
         j2 = READ_row_LB(j)   ! read lower boundary of latitude
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
            n_landwaters_patches(:) = 0
            area_landwaters_patches(:) = 0.
            area_landwaters_grids = 0.

            do nrow = j1, j2
               if(i1 > i2) i2 = i2 + nlon   ! for coarse grid crosses the dateline  
               do ncol = i1, i2
                  ncol_mod = mod(ncol,nlon)
                  if(ncol_mod == 0) ncol_mod = nlon
                  
#if(defined USER_GRID)
                  !-------find out the minimum distance for area weighting--------  
                  area_for_sum = find_min_area(lone_rad(i),lonw_rad(i),lone_rad_i(ncol_mod),&
                                 lonw_rad_i(ncol_mod),sinn(j),sins(j),sinn_i(nrow),sins_i(nrow))
#else
                  area_for_sum = area_fine_gridcell(ncol_mod,nrow)
#endif                  

                  L = landtypes(ncol_mod,nrow)
#if(defined USGS_CLASSIFICATION)
                  if(L==16)then    ! LAND WATER BODIES (16)
                     num_patches(L) = num_patches(L) + 1
                     LL = num_patches(L)
                     n_landwaters_patches (LL) = lakewetland(ncol_mod,nrow)
                     area_landwaters_patches(LL) = area_for_sum
                     area_landwaters_grids = area_landwaters_grids + area_for_sum  
                  endif
#endif
#if(defined IGBP_CLASSIFICATION)
                  if(L==17)then    ! LAND WATER BODIES (17)
                     num_patches(L) = num_patches(L) + 1
                     LL = num_patches(L)
                     n_landwaters_patches (LL) = lakewetland(ncol_mod,nrow)
                     area_landwaters_patches(LL) = area_for_sum
                     area_landwaters_grids = area_landwaters_grids + area_for_sum  
                  endif
#endif
               enddo
            enddo

#if(defined USGS_CLASSIFICATION)
            np = num_patches(16)
#endif
#if(defined IGBP_CLASSIFICATION)
            np = num_patches(17)
#endif
          ! [Lake (1)] -> [reservoir (2)] -> [river (3)]
            fraction_landwaters_patches(:,i,j) = 0.0
            if(np == 0)then
               fraction_landwaters_patches(:,i,j) = 0.0
            else
               num_landwaters(:) = 0
               f_landwaters(:) = 0.
               do n = 1, np
                  nn = n_landwaters_patches(n)
                  if(nn>=1 .and. nn<=3)then
                     num_landwaters(nn) = num_landwaters(nn) + 1
                     f_landwaters(nn) = f_landwaters(nn) + area_landwaters_patches(n)
                  endif
               enddo
               do nnn = 1, 3 
                  fraction_landwaters_patches(nnn,i,j) = f_landwaters(nnn) / area_landwaters_grids
               enddo
            endif

         enddo
      enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

! ---------------------------------------------------
! write out the fraction of wetland patches
! ---------------------------------------------------
      lndname = trim(dir_model_landdata)//'model_landwaters_types.bin'
      print*,lndname
      open(iunit,file=trim(lndname),form='unformatted',status='unknown')
      write(iunit,err=100) fraction_landwaters_patches
      close (iunit)

      deallocate ( landtypes )
      deallocate ( lakewetland )
      deallocate ( num_patches )
      deallocate ( n_landwaters_patches )
      deallocate ( num_landwaters )
      deallocate ( f_landwaters )
      deallocate ( area_landwaters_patches ) 
      deallocate ( fraction_landwaters_patches )

      go to 1000
100   print 101,nrow,lndname
101   format(' record =',i8,',  error occured on file: ',a50)
1000  continue

END SUBROUTINE aggregation_landwaters
