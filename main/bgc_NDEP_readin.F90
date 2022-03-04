#include <define.h>

SUBROUTINE NDEP_readin (lon_points,lat_points,&
                       numpatch,dir_model_landdata,idate)
! ===========================================================
! Read in the LAI, the LAI dataset was created by Yuan et al. (2011)
! http://globalchange.bnu.edu.cn
!
! Created by Yongjiu Dai, March, 2014
! ===========================================================

      use precision
      use MOD_TimeInvariants
      use MOD_TimeVariables
      use MOD_1D_Fluxes, only: ndep_to_sminn
      use omp_lib

      IMPLICIT NONE

      integer, INTENT(in) :: lon_points
      integer, INTENT(in) :: lat_points
      integer, INTENT(in) :: numpatch
      character(LEN=256), INTENT(in) :: dir_model_landdata
      integer, INTENT(in) :: idate(3)

      character(LEN=256) :: c
      character(LEN=256) :: lndname
      integer :: iunit
      integer :: i, j, m, npatch

      real(r8), allocatable :: NDEP_patches(:,:)

! READ in nitrogen deposition data
      allocate (NDEP_patches(1:lon_points,1:lat_points) )

!      iunit = 100
!      write(c,'(i3.3)') Julian_8day
!      lndname = trim(dir_model_landdata)//'model_LAI_patches.'//trim(c)//'.bin'
!      print*,trim(lndname)

!      OPEN(iunit,file=trim(lndname),form='unformatted',status='old')
!      READ(iunit,err=100) LAI_patches
!      CLOSE(iunit)

! added by yuan, 06/20/2016
#ifdef OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(i,j)
#endif
      do npatch = 1, numpatch
         i = patch2lon(npatch)
         j = patch2lat(npatch)
!         ndep_to_smin (npatch) = NDEP_patches(i,j)       !leaf area index
         ndep_to_sminn (npatch) = 4.912128313723134E-009_r8      !leaf area index
!         if(npatch .eq. 79738)print*,'ndep_to_sminn after NDEP_readin',ndep_to_sminn(npatch)
      end do
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      deallocate (NDEP_patches )

END SUBROUTINE NDEP_readin
