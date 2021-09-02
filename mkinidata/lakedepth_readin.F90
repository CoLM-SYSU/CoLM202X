#include <define.h>

SUBROUTINE lakedepth_readin (lon_points,lat_points,dir_model_landdata)

   use precision
   use MOD_TimeInvariants
   USE GlobalVars

   IMPLICIT NONE

   integer, INTENT(in) :: lon_points ! number of longitude points on model grid
   integer, INTENT(in) :: lat_points ! number of latitude points on model grid
   character(LEN=256), INTENT(in) :: dir_model_landdata

   character(len=256) :: lndname
   integer iunit
   real(r8) :: depthratio            ! ratio of lake depth to standard deep lake depth
   integer :: i, j, npatch
   real(r8), allocatable :: lakedepth_patches_in(:,:)

! -----------------------------------------------------------
! For global simulations with 10 body layers,
! the default (50 m lake) body layer thicknesses are given by :
! The node depths zlak located at the center of each layer
  real(r8), dimension(10) :: dzlak
  real(r8), dimension(10) :: zlak

  dzlak = (/0.1, 1., 2., 3., 4., 5., 7., 7., 10.45, 10.45/)  ! m
  zlak  = (/0.05, 0.6, 2.1, 4.6, 8.1, 12.6, 18.6, 25.6, 34.325, 44.775/)

! For site simulations with 25 layers, the default thicknesses are (m):
! real(r8), dimension(25) :: dzlak
! dzlak = (/0.1,                         & ! 0.1 for layer 1;
!           0.25, 0.25, 0.25, 0.25,      & ! 0.25 for layers 2-5;
!           0.50, 0.50, 0.50, 0.50,      & ! 0.5 for layers 6-9;
!           0.75, 0.75, 0.75, 0.75,      & ! 0.75 for layers 10-13;
!           2.00, 2.00,                  & ! 2 for layers 14-15;
!           2.50, 2.50,                  & ! 2.5 for layers 16-17;
!           3.50, 3.50, 3.50, 3.50,      & ! 2.5 for layers 16-17;
!           5.225, 5.225, 5.225, 5.225/)   ! 5.225 for layers 22-25.
!
! For lakes with depth d /= 50 m and d >= 1 m,
!                       the top layer is kept at 10 cm and
!                       the other 9 layer thicknesses are adjusted to maintain fixed proportions.
!
! For lakes with d < 1 m, all layers have equal thickness.
! -----------------------------------------------------------

    ! Read lakedepth
      allocate ( lakedepth_patches_in(lon_points,lat_points) )

      iunit = 100
#ifdef USGS_CLASSIFICATION
      lndname = trim(dir_model_landdata)//'model_GlobalLakeDepth.bin'
#else
      lndname = trim(dir_model_landdata)//'model_GlobalLakeDepth.igbp.bin'
#endif
      print*,lndname
      OPEN(iunit,file=trim(lndname),form='unformatted',status='old')
      read(iunit,err=100) lakedepth_patches_in
      close(iunit)
      print*,'lakedepth =', minval(lakedepth_patches_in, mask = lakedepth_patches_in .gt. -1.0e30), &
                            maxval(lakedepth_patches_in, mask = lakedepth_patches_in .gt. -1.0e30)

    ! Define lake levels
      do npatch = 1, numpatch
         i = patch2lon(npatch)
         j = patch2lat(npatch)
         lakedepth(npatch) = lakedepth_patches_in(i,j)

! testing 04/07/2014
! if(lakedepth(npatch) > 50.) lakedepth(npatch) = 50.

         if(lakedepth(npatch) > 1. .and. lakedepth(npatch) < 1000.)then
            depthratio = lakedepth(npatch) / sum(dzlak(1:nl_lake))
            dz_lake(1,npatch) = dzlak(1)
            dz_lake(2:nl_lake-1,npatch) = dzlak(2:nl_lake-1)*depthratio
            dz_lake(nl_lake,npatch) = dzlak(nl_lake)*depthratio - (dz_lake(1,npatch) - dzlak(1)*depthratio)
         else if(lakedepth(npatch) > 0. .and. lakedepth(npatch) <= 1.)then
            dz_lake(:,npatch) = lakedepth(npatch) / nl_lake
         else   ! non land water bodies or missing value of the lake depth 
            lakedepth(npatch) = sum(dzlak(1:nl_lake))
            dz_lake(1:nl_lake,npatch) = dzlak(1:nl_lake)
         end if
      enddo

      deallocate ( lakedepth_patches_in )

      go to 1000
100   print 101,lndname
101   format(' error occured on file: ',a50)
1000  continue


END SUBROUTINE lakedepth_readin
