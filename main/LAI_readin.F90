#include <define.h>

SUBROUTINE LAI_readin (lon_points,lat_points,&
                       Julian_8day,numpatch,dir_model_landdata)
! ===========================================================
! Read in the LAI, the LAI dataset was created by Yuan et al. (2011)
! http://globalchange.bnu.edu.cn
!
! Created by Yongjiu Dai, March, 2014
! ===========================================================

      use precision
      use MOD_TimeInvariants
      use MOD_TimeVariables
      use omp_lib

      IMPLICIT NONE

      integer, INTENT(in) :: lon_points
      integer, INTENT(in) :: lat_points
      integer, INTENT(in) :: Julian_8day
      integer, INTENT(in) :: numpatch
      character(LEN=256), INTENT(in) :: dir_model_landdata

      character(LEN=256) :: c
      character(LEN=256) :: lndname
      integer :: iunit
      integer :: i, j, m, npatch

      real(r8), allocatable :: LAI_patches(:,:,:)

#if(defined USGS_CLASSIFICATION)
      integer, parameter :: N_land_classification = 24 ! GLCC USGS number of land cover category
      real(r8), dimension(24), parameter :: &   ! Maximum fractional cover of vegetation [-]
      vegc=(/1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, &
             1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, &
             1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0 /)
      real(r8), dimension(24), parameter :: &   ! Stem area index [-]
      sai0=(/0.20, 0.20, 0.30, 0.30, 0.50, 0.50, 1.00, 0.50, &
             1.00, 0.50, 2.00, 2.00, 2.00, 2.00, 2.00, 0.00, &
             2.00, 2.00, 0.00, 0.10, 0.10, 0.10, 0.00, 0.00 /)
#endif
#if(defined IGBP_CLASSIFICATION)
      integer, parameter :: N_land_classification = 17 ! MODIS IGBP number of land cover category
      real(r8), dimension(17), parameter :: &
      vegc=(/1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,&
             1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0/)
      real(r8), dimension(17), parameter :: &
      sai0=(/1.6, 1.8, 1.6, 1.6, 1.5, 1.5, 0.45, 1.4, 1.6, 3.1,&
             1.6, 0.4, 1.1, 1.3, 0.0, 0.14, 0.0/)
#endif

#if(defined USGS_CLASSIFICATION || defined IGBP_CLASSIFICATION)
! READ in Leaf area index and stem area index
      allocate ( LAI_patches(0:N_land_classification,1:lon_points,1:lat_points) )

      iunit = 100
      write(c,'(i3.3)') Julian_8day
      lndname = trim(dir_model_landdata)//'model_LAI_patches.'//trim(c)//'.bin'
      print*,trim(lndname)

      OPEN(iunit,file=trim(lndname),form='unformatted',status='old')
      READ(iunit,err=100) LAI_patches
      CLOSE(iunit)

#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(i,j,m)
#endif
      do npatch = 1, numpatch
         i = patch2lon(npatch)
         j = patch2lat(npatch)
         m = patchclass(npatch)
         if( m == 0 )then
             fveg(npatch)  = 0.
             tlai(npatch)  = 0.
             tsai(npatch)  = 0.
             green(npatch) = 0.
         else
             fveg(npatch)  = vegc(m)    !fraction of veg. cover
             IF (vegc(m) > 0) THEN 
                tlai(npatch)  = LAI_patches(m,i,j)/vegc(m) !leaf area index
                tsai(npatch)  = sai0(m) !stem are index
                green(npatch) = 1.      !fraction of green leaf
             ELSE 
                tlai(npatch)  = 0.  
                tsai(npatch)  = 0.   
                green(npatch) = 0.    
             ENDIF 
         endif
      end do
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      deallocate ( LAI_patches )

      go to 1000
100   print 101,lndname
101   format(' error occured on file: ',a50)
1000  continue
#endif

END SUBROUTINE LAI_readin
