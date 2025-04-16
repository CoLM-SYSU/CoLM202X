
#include <define.h>

SUBROUTINE aggregation_landtypes ( dir_rawdata,dir_model_landdata, &
                                   lon_points,lat_points, &
                                   nrow_start,nrow_end,ncol_start,ncol_end, &
                                   nx_fine_gridcell,ny_fine_gridcell,area_fine_gridcell,&
                                   sinn,sins,lonw_rad,lone_rad,sinn_i,sins_i,lonw_rad_i,lone_rad_i,&
                                   READ_row_UB,READ_row_LB,READ_col_UB,READ_col_LB)
!-----------------------------------------------------------------------
!  Creates land model surface dataset from original "raw" data files -
!     data with 30 arc seconds resolution
!
!  Created by Yongjiu Dai, 02/2014
!
! !REVISIONS:
!  /07/2014, Siguang Zhu & Xiangxiang Zhang: weight average considering
!            partial overlap between fine grid and model grid for a user
!            defined domain file.
!
!-----------------------------------------------------------------------
use MOD_Precision
use MOD_SPMD_Task

IMPLICIT NONE

! arguments:
#if (defined LULC_USGS)
      integer, parameter :: N_land_classification = 24 ! GLCC USGS number of land cover category
#endif
#if (defined LULC_IGBP)
      integer, parameter :: N_land_classification = 17 ! MODIS IGBP number of land cover category
#endif
      integer, parameter :: nlat = 21600  ! 180*(60*2)
      integer, parameter :: nlon = 43200  ! 360*(60*2)

      character(len=256), intent(in) :: dir_rawdata
      character(len=256), intent(in) :: dir_model_landdata

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

      real(r8), intent(in) :: area_fine_gridcell(nlon,nrow_start:nrow_end)  ! rwadata fine cell area (km**2)

! local variables:
! ---------------------------------------------------------------
      character(len=256) lndname
#if (defined usempi)
      character, allocatable :: land_chr1(:,:)
      integer(kind=2), allocatable ::  land_int2(:,:)
#else
      character(len=1) land_chr1(nlon)
      character(len=2) land_chr2(nlon)
      integer(kind=1)  land_int1(nlon)
      integer(kind=2)  land_int2(nlon)
#endif

      integer iunit
      integer length
      integer i, j, L, i1, i2, j1, j2
      integer nrow, ncol, ncol_mod
      integer Loca(1), nn

      integer, allocatable  :: landtypes (:,:)
      real(r8), allocatable :: glacier(:,:)
      real(r8), allocatable :: fraction_patches(:,:,:)

      real(r8) area_grids
      real(r8) area_for_sum
      real(r8), allocatable :: area_patches(:)

      real(r8) g_patches
      real(r8) a_glacier_patches
      real(r8) f_glacier_patches
      real(r8) err_f_glacier

      real(r8), external :: find_min_area

#if (defined usempi)
      integer(kind=MPI_OFFSET_KIND) :: fdisp
      integer(kind=MPI_OFFSET_KIND) :: fdlast
#endif

#if (defined USE_POINT_DATA)

   allocate (fraction_patches(0:N_land_classification,1:lon_points,1:lat_points))
   fraction_patches(:,:,:) = 0.0

#if (defined LULC_USGS)
   fraction_patches(LULC_USGS,:,:) = 1.0
#else
   fraction_patches(LULC_IGBP,:,:) = 1.0
#endif

#else
! ........................................
! ... (1) gloabl land cover types
! ........................................

      allocate (landtypes(nlon,nrow_start:nrow_end))

#if (defined LULC_USGS)
     ! GLCC USGS classification
      lndname = trim(dir_rawdata)//'RAW_DATA_updated/landtypes_usgs_update.bin'
#endif
#if (defined LULC_IGBP)
     ! MODIS IGBP classification
      lndname = trim(dir_rawdata)//'RAW_DATA_updated/landtypes_igbp_update.bin'
#endif

      if (p_master) print*,trim(lndname)

#if (defined usempi)
      allocate (land_chr1(nlon,nrow_start:nrow_end))
      fdisp = 0
      call mpi_rdwr_data (lndname, 'read', fdisp, fine_lat_map, nlon, land_chr1)
      landtypes(:,:) = ichar(land_chr1(:,:))
      deallocate (land_chr1)
#else
      inquire(iolength=length) land_chr1
      iunit = 100
      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
      do nrow = nrow_start, nrow_end
         read(iunit,rec=nrow,err=100) land_chr1
         landtypes(:,nrow) = ichar(land_chr1(:))
      enddo
      close (iunit)
#endif


! ................................................
! ... (2) specific global lake/wetland/urban extent/glacier and ice sheet coverage
! ................................................
!   WAITING        WAITING            WAITING
!!  WAITING        WAITING            WAITING
!!! WAITING FOR THE coverage data OF GLOBAL lake/wetland/urban and built-up!!!

      allocate (glacier(nlon,nrow_start:nrow_end))

      lndname = trim(dir_rawdata)//'glacier/glacier.bin'
      if (p_master) print*,trim(lndname)

#if (defined usempi)
      allocate (land_int2(nlon,nrow_start:nrow_end))
      fdisp = 0
      call mpi_rdwr_data (lndname, 'read', fdisp, fine_lat_map, nlon, land_int2)
      glacier(:,:) = land_int2(:,:) * 0.1
      deallocate (land_int2)
#else
      inquire(iolength=length) land_int2
      iunit = 100
      open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
      do nrow = nrow_start, nrow_end
         read(iunit,rec=nrow,err=100) land_int2
         glacier(:,nrow) = land_int2(:) * 0.1
      enddo
      close (iunit)
#endif

      call print_minmax('glacier = ', glacier)

! ........................................................................................
! ... (3) aggregate the land types from the resolution of raw data to modelling resolution
! ........................................................................................

      allocate (area_patches(0:N_land_classification))
      allocate (fraction_patches(0:N_land_classification,1:lon_points,1:lat_points))
      fraction_patches(:,:,:) = 0.

#ifdef OPENMP
print *, 'OPENMP enabled, threads num = ', OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(i,j,i1,i2,j1,j2,nrow,ncol,ncol_mod,L) &
!$OMP PRIVATE(area_patches,area_grids,area_for_sum) &
!$OMP PRIVATE(a_glacier_patches,g_patches,f_glacier_patches,err_f_glacier,Loca,nn)
#endif
      do j = 1, lat_points
!----------------------modified by zsg--------------------------------------
#if (defined USER_GRID)
         j1 = READ_row_UB(j)   ! read upper boundary of latitude
         j2 = READ_row_LB(j)   ! read lower boundary of latitude
#else
         j1 = nrow_start + (j-1)*ny_fine_gridcell
         j2 = nrow_start - 1 + j*ny_fine_gridcell
#endif
         do i = 1, lon_points

#if (defined USER_GRID)
            i1 = READ_col_UB(i)   ! read upper boundary of longitude
            i2 = READ_col_LB(i)   ! read lower boundary of longitude
#else
            i1 = ncol_start + (i-1)*nx_fine_gridcell
            i2 = ncol_start -1 + i*nx_fine_gridcell
#endif
            area_patches(:) = 0.
            area_grids = 0.
            a_glacier_patches = 0.

            do nrow = j1, j2
               if(i1 > i2) i2 = i2 + nlon   ! for coarse grid crosses the dateline
               do ncol = i1, i2
                  ncol_mod = mod(ncol,nlon)
                  if(ncol_mod == 0) ncol_mod = nlon

#if (defined USER_GRID)
                  !-------find out the minimum distance for area weighting--------
                  area_for_sum = find_min_area(lone_rad(i),lonw_rad(i),lone_rad_i(ncol_mod),&
                                 lonw_rad_i(ncol_mod),sinn(j),sins(j),sinn_i(nrow),sins_i(nrow))
#else
                  area_for_sum = area_fine_gridcell(ncol_mod,nrow)
#endif

                  L = landtypes(ncol_mod,nrow)
                  area_patches(L) = area_patches(L) + area_for_sum
                  area_grids = area_grids + area_for_sum

#if (defined LULC_USGS)
                  if(L==24)then  ! GLACIER/ICE SHEET (24)
#endif
#if (defined LULC_IGBP)
                  if(L==15)then  ! GLACIER/ICE SHEET (15)
#endif
                     g_patches = 100.
                     if(glacier(ncol_mod,nrow) >= 0.1)then
                        g_patches = glacier(ncol_mod,nrow)
                     endif
                     a_glacier_patches = a_glacier_patches + g_patches/100. * area_for_sum
                  endif
               enddo
            enddo

            do L = 0, N_land_classification
               fraction_patches(L,i,j) = area_patches(L) / area_grids
            enddo

            f_glacier_patches = a_glacier_patches / area_grids
            if(f_glacier_patches > 0.001)then     ! 0.1/100
#if (defined LULC_USGS)
               err_f_glacier = fraction_patches(24,i,j) - f_glacier_patches
               fraction_patches(24,i,j) = f_glacier_patches  ! GLCC USGS GLACIER/ICE SHEET (24)
#endif
#if (defined LULC_IGBP)
               err_f_glacier = fraction_patches(15,i,j) - f_glacier_patches
               fraction_patches(15,i,j) = f_glacier_patches  ! MODIS IGBP GLACIER/ICE SHEET (15)
#endif
               Loca = maxloc(fraction_patches(:,i,j))    ! maxloc get the Loca: 1 - N_land_classification + 1
               nn = Loca(1) - 1                          ! the definition of dimension of fraction_patches: 0 - N_land_classification
               fraction_patches(nn,i,j) = err_f_glacier + fraction_patches(nn,i,j)
            endif

         enddo
      enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      deallocate (landtypes)
      deallocate (glacier)
      deallocate (area_patches)

#endif

! Write-out the fraction of the land types in the gridcells
      lndname = trim(dir_model_landdata)//'model_landtypes.bin'
      if (p_master) print*,trim(lndname)

#if (defined usempi)
      fdisp = 4
      fdlast = 0
      call mpi_rdwr_data (lndname, 'write', fdisp, crse_lat_map, &
                           (N_land_classification+1)*lon_points, fraction_patches)
      call move_next_for_unformatted (lndname, 'write', fdlast, fdisp)
#else
      iunit = 100
      open(iunit,file=trim(lndname),form='unformatted',status='unknown')
      write(iunit,err=100) fraction_patches
      close(iunit)
#endif

      deallocate (fraction_patches)

      go to 1000
100   print 101,nrow,lndname
101   format(' record =',i8,',  error occured on file: ',a50)
1000  continue

END SUBROUTINE aggregation_landtypes
!-----------------------------------------------------------------------
!EOP
