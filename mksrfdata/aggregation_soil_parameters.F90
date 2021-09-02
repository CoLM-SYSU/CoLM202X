
#include <define.h>

SUBROUTINE aggregation_soil_parameters ( dir_rawdata,dir_model_landdata, &
                                         lon_points,lat_points, &
                                         nrow_start,nrow_end,ncol_start,ncol_end, &
                                         nx_fine_gridcell,ny_fine_gridcell,area_fine_gridcell,&
                                         READ_row_UB,READ_row_LB,READ_col_UB,READ_col_LB )
! ----------------------------------------------------------------------
! Creates land model surface dataset from original "raw" data files -
!     data with 30 arc seconds resolution
!
! Created by Yongjiu Dai, 02/2014
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
      character(len=2) land_chr2(nlon)
      character(len=256) c
      CHARACTER(len=256) suffix

      real(r8) tmp(nlon)
      integer iunit
      integer length
      integer i, j, L
      integer i1, i2, j1, j2
      integer nrow, ncol, ncol_mod, nrow_mod
      integer LL, LL0, np
      integer n_fine_gridcell
      integer nsl, MODEL_SOIL_LAYER

      integer, allocatable :: landtypes (:,:)
      integer, allocatable :: num_patches(:) 

      real(r8), allocatable :: theta_s_l      (:,:)
      real(r8), allocatable :: psi_s_l        (:,:)
      real(r8), allocatable :: lambda_l       (:,:)
      real(r8), allocatable :: k_s_l          (:,:)
      real(r8), allocatable :: csol_l         (:,:)
      real(r8), allocatable :: tksatu_l       (:,:)
      real(r8), allocatable :: tkdry_l        (:,:)
     
      real(r8), allocatable :: a_theta_s_l    (:,:) 
      real(r8), allocatable :: a_psi_s_l      (:,:) 
      real(r8), allocatable :: a_lambda_l     (:,:) 
      real(r8), allocatable :: a_k_s_l        (:,:) 
      real(r8), allocatable :: a_csol_l       (:,:) 
      real(r8), allocatable :: a_tksatu_l     (:,:) 
      real(r8), allocatable :: a_tkdry_l      (:,:)

      real(r8), allocatable :: soil_theta_s_l (:,:,:) 
      real(r8), allocatable :: soil_psi_s_l   (:,:,:) 
      real(r8), allocatable :: soil_lambda_l  (:,:,:) 
      real(r8), allocatable :: soil_k_s_l     (:,:,:) 
      real(r8), allocatable :: soil_csol_l    (:,:,:) 
      real(r8), allocatable :: soil_tksatu_l  (:,:,:) 
      real(r8), allocatable :: soil_tkdry_l   (:,:,:) 

      real(r8), external :: median

! ........................................
! ... [1] gloabl land cover types
! ........................................
      iunit = 100
      inquire(iolength=length) land_chr1
      allocate (landtypes(nlon,nlat))

#if(defined USE_POINT_DATA)

!TODO: need modification for point case
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
      suffix        = ''
#else
      suffix        = '.igbp'
#endif


! ........................................
! ... [2] aggregate the soil parameters from the resolution of raw data to modelling resolution
! ........................................
      n_fine_gridcell = nx_fine_gridcell * ny_fine_gridcell

      allocate ( num_patches (0:N_land_classification) )

      allocate ( theta_s_l (nlon,nlat) )
      allocate ( psi_s_l   (nlon,nlat) )
      allocate ( lambda_l  (nlon,nlat) )
      allocate ( k_s_l     (nlon,nlat) )
      allocate ( csol_l    (nlon,nlat) )
      allocate ( tksatu_l  (nlon,nlat) )
      allocate ( tkdry_l   (nlon,nlat) )

      allocate ( a_theta_s_l (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_psi_s_l   (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_lambda_l  (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_k_s_l     (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_csol_l    (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_tksatu_l  (0:N_land_classification,1:n_fine_gridcell) )
      allocate ( a_tkdry_l   (0:N_land_classification,1:n_fine_gridcell) )

      allocate ( soil_theta_s_l (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_psi_s_l   (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_lambda_l  (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_k_s_l     (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_csol_l    (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_tksatu_l  (0:N_land_classification,1:lon_points,1:lat_points) )
      allocate ( soil_tkdry_l   (0:N_land_classification,1:lon_points,1:lat_points) )

      iunit = 100
      DO nsl = 1, 8
         MODEL_SOIL_LAYER = nsl
         write(c,'(i1)') MODEL_SOIL_LAYER

! (1) Read in the saturated water content [cm3/cm3]
         inquire(iolength=length) tmp
         lndname = trim(dir_rawdata)//'RAW_DATA_updated_with_igbp/theta_s_l'//trim(c)//trim(suffix)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) theta_s_l(:,nrow)
         enddo
         close(iunit)

! (2) Read in the matric potential at saturation [cm]
         inquire(iolength=length) tmp
         lndname = trim(dir_rawdata)//'RAW_DATA_updated_with_igbp/psi_s_l'//trim(c)//trim(suffix)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) psi_s_l(:,nrow)
         enddo
         close(iunit)

! (3) Read in the pore size distribution index [dimensionless]
         inquire(iolength=length) tmp
         lndname = trim(dir_rawdata)//'RAW_DATA_updated_with_igbp/lambda_l'//trim(c)//trim(suffix)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) lambda_l(:,nrow)
         enddo
         close(iunit)

! (4) Read in the saturated hydraulic conductivity [cm/day]
         inquire(iolength=length) tmp
         lndname = trim(dir_rawdata)//'RAW_DATA_updated_with_igbp/k_s_l'//trim(c)//trim(suffix)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) k_s_l(:,nrow)
         enddo
         close(iunit)

! (5) Read in the heat capacity of soil solids [J/(m3 K)]
         inquire(iolength=length) tmp
         lndname = trim(dir_rawdata)//'RAW_DATA_updated_with_igbp/csol_l'//trim(c)//trim(suffix)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) csol_l(:,nrow)
         enddo
         close(iunit)

! (6) Read in the thermal conductivity of saturated soil [W/m-K]
         inquire(iolength=length) tmp
         lndname = trim(dir_rawdata)//'RAW_DATA_updated_with_igbp/tksatu_l'//trim(c)//trim(suffix)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) tksatu_l(:,nrow)
         enddo
         close(iunit)

! (7) Read in the thermal conductivity for dry soil [W/(m-K)]
         inquire(iolength=length) tmp
         lndname = trim(dir_rawdata)//'RAW_DATA_updated_with_igbp/tkdry_l'//trim(c)//trim(suffix)
         print*,lndname
         open(iunit,file=trim(lndname),access='direct',recl=length,form='unformatted',status='old')
         do nrow = nrow_start, nrow_end
            read(iunit,rec=nrow,err=100) tkdry_l(:,nrow)
         enddo
         close(iunit)


#ifdef OPENMP
print *, 'OPENMP enabled, threads num = ', OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) SCHEDULE(DYNAMIC,1) &
!$OMP PRIVATE(i,j,i1,i2,j1,j2,nrow,ncol,ncol_mod,nrow_mod,L,LL,LL0,num_patches,np) &
!$OMP PRIVATE(a_theta_s_l,a_psi_s_l,a_lambda_l,a_k_s_l,a_csol_l,a_tksatu_l,a_tkdry_l)  
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
! yuan, 07/30/2019: total grid soil info stored to position 0 (original: ocean)
                     num_patches(0) = num_patches(0) + 1
                     LL0 = num_patches(0) 

                     a_theta_s_l (L,LL) = theta_s_l(ncol_mod,nrow_mod) 
                     a_psi_s_l   (L,LL) = psi_s_l  (ncol_mod,nrow_mod)
                     a_lambda_l  (L,LL) = lambda_l (ncol_mod,nrow_mod)
                     a_k_s_l     (L,LL) = k_s_l    (ncol_mod,nrow_mod)
                     a_csol_l    (L,LL) = csol_l   (ncol_mod,nrow_mod)
                     a_tksatu_l  (L,LL) = tksatu_l (ncol_mod,nrow_mod)
                     a_tkdry_l   (L,LL) = tkdry_l  (ncol_mod,nrow_mod)

                     a_theta_s_l (0,LL0) = theta_s_l(ncol_mod,nrow_mod) 
                     a_psi_s_l   (0,LL0) = psi_s_l  (ncol_mod,nrow_mod)
                     a_lambda_l  (0,LL0) = lambda_l (ncol_mod,nrow_mod)
                     a_k_s_l     (0,LL0) = k_s_l    (ncol_mod,nrow_mod)
                     a_csol_l    (0,LL0) = csol_l   (ncol_mod,nrow_mod)
                     a_tksatu_l  (0,LL0) = tksatu_l (ncol_mod,nrow_mod)
                     a_tkdry_l   (0,LL0) = tkdry_l  (ncol_mod,nrow_mod)

                  enddo
               enddo

               do L = 0, N_land_classification 

                  if(L.ge.0)then  ! include (0) to calculate the whole grid
                     np = num_patches(L) 
                     if(np == 0)then
                        soil_theta_s_l (L,i,j) = -1.0e36
                        soil_psi_s_l   (L,i,j) = -1.0e36
                        soil_lambda_l  (L,i,j) = -1.0e36
                        soil_k_s_l     (L,i,j) = -1.0e36
                        soil_csol_l    (L,i,j) = -1.0e36
                        soil_tksatu_l  (L,i,j) = -1.0e36
                        soil_tkdry_l   (L,i,j) = -1.0e36
                     else if(np == 1) then
                        soil_theta_s_l (L,i,j) = a_theta_s_l(L,1)
                        soil_psi_s_l   (L,i,j) = a_psi_s_l  (L,1)
                        soil_lambda_l  (L,i,j) = a_lambda_l (L,1)
                        soil_k_s_l     (L,i,j) = a_k_s_l    (L,1)
                        soil_csol_l    (L,i,j) = a_csol_l   (L,1)
                        soil_tksatu_l  (L,i,j) = a_tksatu_l (L,1)
                        soil_tkdry_l   (L,i,j) = a_tkdry_l  (L,1)
                     else
                        soil_theta_s_l (L,i,j) = median ( a_theta_s_l(L,1:np), np) 
                        soil_psi_s_l   (L,i,j) = median ( a_psi_s_l  (L,1:np), np)
                        soil_lambda_l  (L,i,j) = median ( a_lambda_l (L,1:np), np)
                        soil_k_s_l     (L,i,j) = median ( a_k_s_l    (L,1:np), np)
                        soil_csol_l    (L,i,j) = median ( a_csol_l   (L,1:np), np)
                        soil_tksatu_l  (L,i,j) = median ( a_tksatu_l (L,1:np), np)
                        soil_tkdry_l   (L,i,j) = median ( a_tkdry_l  (L,1:np), np)
                     endif

                  else          ! OCEAN
                     soil_theta_s_l (L,i,j) = -1.0e36
                     soil_psi_s_l   (L,i,j) = -1.0e36
                     soil_lambda_l  (L,i,j) = -1.0e36
                     soil_k_s_l     (L,i,j) = -1.0e36
                     soil_csol_l    (L,i,j) = -1.0e36
                     soil_tksatu_l  (L,i,j) = -1.0e36
                     soil_tkdry_l   (L,i,j) = -1.0e36
                  endif
               enddo

            enddo
         enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

! (1) Write-out the saturated water content [cm3/cm3]
         lndname = trim(dir_model_landdata)//'model_theta_s_l'//trim(c)//trim(suffix)//'.bin'
         print*,lndname
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) soil_theta_s_l
         close(iunit)

! (2) Write-out the matric potential at saturation [cm]
         lndname = trim(dir_model_landdata)//'model_psi_s_l'//trim(c)//trim(suffix)//'.bin'
         print*,lndname
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) soil_psi_s_l
         close(iunit)

! (3) Write-out the pore size distribution index [dimensionless]
         lndname = trim(dir_model_landdata)//'model_lambda_l'//trim(c)//trim(suffix)//'.bin'
         print*,lndname
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) soil_lambda_l
         close(iunit)

! (4) Write-out the saturated hydraulic conductivity [cm/day]
         lndname = trim(dir_model_landdata)//'model_k_s_l'//trim(c)//trim(suffix)//'.bin'
         print*,lndname
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) soil_k_s_l
         close(iunit)

! (5) Write-out the heat capacity of soil solids [J/(m3 K)]
         lndname = trim(dir_model_landdata)//'model_csol_l'//trim(c)//trim(suffix)//'.bin'
         print*,lndname
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) soil_csol_l
         close(iunit)

! (6) Write-out the thermal conductivity of saturated soil [W/m-K]
         lndname = trim(dir_model_landdata)//'model_tksatu_l'//trim(c)//trim(suffix)//'.bin'
         print*,lndname
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) soil_tksatu_l
         close(iunit)

! (7) Write-out the thermal conductivity for dry soil [W/(m-K)]
         lndname = trim(dir_model_landdata)//'model_tkdry_l'//trim(c)//trim(suffix)//'.bin'
         print*,lndname
         open(iunit,file=trim(lndname),form='unformatted',status='unknown')
         write(iunit,err=100) soil_tkdry_l
         close(iunit)

      ENDDO

! Deallocate the allocatable array
! --------------------------------
      deallocate ( num_patches )

      deallocate ( theta_s_l   )
      deallocate ( psi_s_l     )
      deallocate ( lambda_l    )
      deallocate ( k_s_l       )
      deallocate ( csol_l      )
      deallocate ( tksatu_l    )
      deallocate ( tkdry_l     )

      deallocate ( a_theta_s_l )
      deallocate ( a_psi_s_l   )
      deallocate ( a_lambda_l  )
      deallocate ( a_k_s_l     )
      deallocate ( a_csol_l    )
      deallocate ( a_tksatu_l  )
      deallocate ( a_tkdry_l   )

      deallocate ( soil_theta_s_l )
      deallocate ( soil_psi_s_l   )
      deallocate ( soil_lambda_l  )
      deallocate ( soil_k_s_l     )
      deallocate ( soil_csol_l    )
      deallocate ( soil_tksatu_l  )
      deallocate ( soil_tkdry_l   )

      deallocate ( landtypes )

      go to 1000
100   print 101,nrow,lndname
101   format(' record =',i8,',  error occured on file: ',a50)
1000  continue

END SUBROUTINE aggregation_soil_parameters
!-----------------------------------------------------------------------
!EOP
