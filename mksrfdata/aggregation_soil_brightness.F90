
#include <define.h>

SUBROUTINE aggregation_soil_brightness ( dir_rawdata,dir_model_landdata, &
                                         lon_points,lat_points, &
                                         nrow_start,nrow_end,ncol_start,ncol_end, &
                                         nx_fine_gridcell,ny_fine_gridcell,area_fine_gridcell,&
                                         READ_row_UB,READ_row_LB,READ_col_UB,READ_col_LB )
! ----------------------------------------------------------------------
! Creates land model surface dataset from original "raw" data files -
!     data with 30 arc seconds resolution
!
! Created by Yongjiu Dai, 03/2014
! ----------------------------------------------------------------------
use precision
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
      character(len=1) soil_chr1(nlon30s)
      character(len=256) c
      CHARACTER(len=256) suffix

      integer iunit
      integer length
      integer i, j, L
      integer i1, i2, j1, j2
      integer nrow, ncol, ncol_mod, nrow_mod
      integer LL, LL0, np, ii
      integer n_fine_gridcell

      INTEGER nrow30s_start, nrow30s_end
      INTEGER ncol30s_start, ncol30s_end

      integer, allocatable :: landtypes (:,:)
      integer, allocatable :: isc(:,:) 
      integer, allocatable :: num_patches(:) 

      real(r8), allocatable :: a_s_v_refl (:,:)
      real(r8), allocatable :: a_d_v_refl (:,:)
      real(r8), allocatable :: a_s_n_refl (:,:)
      real(r8), allocatable :: a_d_n_refl (:,:)

      real(r8), allocatable :: soil_s_v_alb (:,:,:)
      real(r8), allocatable :: soil_d_v_alb (:,:,:)
      real(r8), allocatable :: soil_s_n_alb (:,:,:)
      real(r8), allocatable :: soil_d_n_alb (:,:,:)

      real(r8), external :: median

! ----------------------------------------------------------------------
! The soil color and reflectance is from the work:
! Peter J. Lawrence and Thomas N. Chase, 2007:
! Representing a MODIS consistent land surface in the Community Land Model (CLM 3.0):
! Part 1 generating MODIS consistent land surface parameters
      real(r8) soil_s_v_refl(20) ! Saturated visible soil reflectance
      real(r8) soil_d_v_refl(20) ! Dry visible soil reflectance
      real(r8) soil_s_n_refl(20) ! Saturated near infrared soil reflectance
      real(r8) soil_d_n_refl(20) ! Dry near infrared soil reflectance

      soil_s_v_refl = (/ 0.26, 0.24, 0.22, 0.20, 0.19, 0.18, 0.17, 0.16, 0.15, 0.14, &
                         0.13, 0.12, 0.11, 0.10, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04 /)

      soil_d_v_refl = (/ 0.37, 0.35, 0.33, 0.31, 0.30, 0.29, 0.28, 0.27, 0.26, 0.25, &
                         0.24, 0.23, 0.22, 0.21, 0.20, 0.19, 0.18, 0.17, 0.16, 0.15 /)

      soil_s_n_refl = (/ 0.52, 0.48, 0.44, 0.40, 0.38, 0.36, 0.34, 0.32, 0.30, 0.28, &
                         0.26, 0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10, 0.08 /)

      soil_d_n_refl = (/ 0.63, 0.59, 0.55, 0.51, 0.49, 0.47, 0.45, 0.43, 0.41, 0.39, &
                         0.37, 0.35, 0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.19 /)

! ........................................
! ... [1] gloabl land cover types
! ........................................
      iunit = 100
      inquire(iolength=length) land_chr1
      allocate (landtypes(nlon,nlat))

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
      nrow30s_start = nrow_Start
      nrow30s_end   = nrow_end
      ncol30s_start = ncol_start
      ncol30s_end   = ncol_end
      suffix        = ''
#else
      nrow30s_start = int((nrow_Start+1)/2)
      nrow30s_end   = int((nrow_end+1)/2)
      ncol30s_start = int((ncol_start+1)/2)
      ncol30s_end   = int((ncol_end+1)/2)
      suffix        = '.igbp'
#endif

! ........................................
! ... [2] aggregate the soil parameters from the resolution of raw data to modelling resolution
! ........................................
      n_fine_gridcell = nx_fine_gridcell * ny_fine_gridcell

      allocate ( isc (nlon30s,nlat30s) )
      allocate ( num_patches (0:N_land_classification) )
      allocate ( a_s_v_refl (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_d_v_refl (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_s_n_refl (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_d_n_refl (0:N_land_classification,1:n_fine_gridcell) )

      allocate ( soil_s_v_alb (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_d_v_alb (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_s_n_alb (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_d_n_alb (0:N_land_classification,1:lon_points,1:lat_points) )

! Read in the index of soil brightness (color)
      iunit = 100
      inquire(iolength=length) soil_chr1
      lndname = trim(dir_rawdata)//'soil_brightness/soilcol_clm_30s.bin'
      print*,lndname
      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
      do nrow = nrow30s_start, nrow30s_end
         read(iunit,rec=nrow,err=100) soil_chr1
         isc(:,nrow) = ichar(soil_chr1(:))
      enddo
      close(iunit)

#ifdef OPENMP
print *, 'OPENMP enabled, threads num = ', OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(i,j,i1,i2,j1,j2,nrow,ncol,ncol_mod,nrow_mod,L,LL,LL0,num_patches,ii,np) &
!$OMP PRIVATE(a_s_v_refl,a_d_v_refl,a_s_n_refl,a_d_n_refl)
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
            do nrow = j1, j2            
               if(i1 > i2) i2 = i2 + nlon   ! for coarse grid crosses the dateline    
               do ncol = i1, i2
                  ncol_mod = mod(ncol,nlon)
                  nrow_mod = nrow
                  if(ncol_mod == 0) ncol_mod = nlon
                  L = landtypes(ncol_mod,nrow)

                  IF (L == 0) cycle

                  num_patches(L) = num_patches(L) + 1
                  LL = num_patches(L) 
! yuan, 07/30/2019: total grid soil info stored to position 0 (ocean)
                  num_patches(0) = num_patches(0) + 1
                  LL0 = num_patches(0) 

! yuan, 07/30/2019: 500m (15s) ==> 1km (30s)
#ifndef USGS_CLASSIFICATION
                  ncol_mod = int((ncol_mod+1)/2)
                  nrow_mod = int((nrow+1)/2)
#endif
                  ii = isc(ncol_mod,nrow_mod)

                  a_s_v_refl (L,LL) = soil_s_v_refl( ii )
                  a_d_v_refl (L,LL) = soil_d_v_refl( ii )
                  a_s_n_refl (L,LL) = soil_s_n_refl( ii )
                  a_d_n_refl (L,LL) = soil_d_n_refl( ii )

                  a_s_v_refl (0,LL0) = soil_s_v_refl( ii )
                  a_d_v_refl (0,LL0) = soil_d_v_refl( ii )
                  a_s_n_refl (0,LL0) = soil_s_n_refl( ii )
                  a_d_n_refl (0,LL0) = soil_d_n_refl( ii )
               enddo
            enddo
            
            do L = 0, N_land_classification 
#if(defined USGS_CLASSIFICATION)
               if(L/=16 .and. L/=24)then  ! NOT OCEAN(0)/WATER BODIES(16)/GLACIER and ICESHEET(24)
#else
               if(L/=17 .and. L/=15)then  ! NOT OCEAN(0)/WATER BODIES(17)/GLACIER and ICE SHEET(15)
#endif
                  np = num_patches(L) 
                  if(np == 0)then
! yuan, 12/29/2019: bug, inconsistant with land fraction and land cover data
                     ! set to color index 15 (for ocean fill value)
                     ! need to check with soil hydro/thermal paras
                     !soil_s_v_alb (L,i,j) = -1.0e36
                     !soil_d_v_alb (L,i,j) = -1.0e36
                     !soil_s_n_alb (L,i,j) = -1.0e36
                     !soil_d_n_alb (L,i,j) = -1.0e36
                     soil_s_v_alb (L,i,j) = soil_s_v_refl( 15 )
                     soil_d_v_alb (L,i,j) = soil_s_v_refl( 15 )
                     soil_s_n_alb (L,i,j) = soil_s_v_refl( 15 )
                     soil_d_n_alb (L,i,j) = soil_s_v_refl( 15 )
                  else if(np == 1) then
                     soil_s_v_alb (L,i,j) = a_s_v_refl(L,1)
                     soil_d_v_alb (L,i,j) = a_d_v_refl(L,1)
                     soil_s_n_alb (L,i,j) = a_s_n_refl(L,1)
                     soil_d_n_alb (L,i,j) = a_d_n_refl(L,1)
                  else
                     soil_s_v_alb (L,i,j) = median ( a_s_v_refl(L,1:np), np) 
                     soil_d_v_alb (L,i,j) = median ( a_d_v_refl(L,1:np), np)
                     soil_s_n_alb (L,i,j) = median ( a_s_n_refl(L,1:np), np)
                     soil_d_n_alb (L,i,j) = median ( a_d_n_refl(L,1:np), np)
                  endif

               else
                  soil_s_v_alb (L,i,j) = -1.0e36
                  soil_d_v_alb (L,i,j) = -1.0e36
                  soil_s_n_alb (L,i,j) = -1.0e36
                  soil_d_n_alb (L,i,j) = -1.0e36
               endif
            enddo

         enddo
      enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      print*,'s_v_alb =', minval(soil_s_v_alb, mask = soil_s_v_alb .gt. -1.e30), &
                          maxval(soil_s_v_alb, mask = soil_s_v_alb .gt. -1.e30)
      print*,'d_v_alb =', minval(soil_d_v_alb, mask = soil_d_v_alb .gt. -1.e30), &
                          maxval(soil_d_v_alb, mask = soil_d_v_alb .gt. -1.e30)
      print*,'s_n_alb =', minval(soil_s_n_alb, mask = soil_s_n_alb .gt. -1.e30), &
                          maxval(soil_s_n_alb, mask = soil_s_n_alb .gt. -1.e30)
      print*,'d_n_alb =', minval(soil_d_n_alb, mask = soil_d_n_alb .gt. -1.e30), &
                          maxval(soil_d_n_alb, mask = soil_d_n_alb .gt. -1.e30)

! (1) Write-out the albedo of visible of the saturated soil
      lndname = trim(dir_model_landdata)//'soil_s_v_alb'//trim(suffix)//'.bin'
      print*,lndname
      open(iunit,file=trim(lndname),form='unformatted',status='unknown')
      write(iunit,err=100) soil_s_v_alb
      close(iunit)

! (1) Write-out the albedo of visible of the dry soil
      lndname = trim(dir_model_landdata)//'soil_d_v_alb'//trim(suffix)//'.bin'
      print*,lndname
      open(iunit,file=trim(lndname),form='unformatted',status='unknown')
      write(iunit,err=100) soil_d_v_alb
      close(iunit)

! (3) Write-out the albedo of near infrared of the saturated soil
      lndname = trim(dir_model_landdata)//'soil_s_n_alb'//trim(suffix)//'.bin'
      print*,lndname
      open(iunit,file=trim(lndname),form='unformatted',status='unknown')
      write(iunit,err=100) soil_s_n_alb
      close(iunit)

! (4) Write-out the albedo of near infrared of the dry soil
      lndname = trim(dir_model_landdata)//'soil_d_n_alb'//trim(suffix)//'.bin'
      print*,lndname
      open(iunit,file=trim(lndname),form='unformatted',status='unknown')
      write(iunit,err=100) soil_d_n_alb
      close(iunit)

! Deallocate the allocatable array
! --------------------------------
      deallocate ( landtypes )
      deallocate ( isc )
      deallocate ( num_patches )

      deallocate ( a_s_v_refl )
      deallocate ( a_d_v_refl )
      deallocate ( a_s_n_refl )
      deallocate ( a_d_n_refl )

      deallocate ( soil_s_v_alb )
      deallocate ( soil_d_v_alb )
      deallocate ( soil_s_n_alb )
      deallocate ( soil_d_n_alb )

      go to 1000
100   print 101,nrow,lndname
101   format(' record =',i8,',  error occured on file: ',a50)
1000  continue

END SUBROUTINE aggregation_soil_brightness
!-----------------------------------------------------------------------
!EOP
