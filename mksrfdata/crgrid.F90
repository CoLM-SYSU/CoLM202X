#include <define.h>

SUBROUTINE crgrid(dir_model_landdata,edgen,edgee,edges,edgew,&
                  lon_points,lat_points,latn,lats,lonw,lone)
! ----------------------------------------------------------------------
! generate land model grid when mode is offline.
! surface grid edges -- grids do not have to be global. 
! to allow this, grids must define the north, east, south, west edges:
!    edgen: northern edge of grid : > -90 and <= 90 (degrees)
!    edgee: eastern edge of grid  : > western edge and <= 180
!    edges: southern edge of grid : >= -90  and <  90
!    edgew: western edge of grid  : >= -180 and < 180
!
! region (global) latitude grid goes from:
!                          NORTHERN edge (POLE) to SOUTHERN edge (POLE)
! region (global) longitude grid starts at:
!                          WESTERN edge 
!                         (DATELINE with western edge)
!
! west of Greenwich defined negative
! for global grids, the western edge of the longitude grid starts 
! at the dateline 
!
! ----------------------------------------------------------------------
use precision

IMPLICIT NONE

! arguments:
      character(len=256), intent(in) :: dir_model_landdata
      integer, intent(in) :: lon_points  ! input number of longitudes
      integer, intent(in) :: lat_points  ! input number of latitudes
      real(r8), intent(in) :: edgen      ! northern edge of grid (degrees)
      real(r8), intent(in) :: edgee      ! eastern edge of grid (degrees)
      real(r8), intent(in) :: edges      ! southern edge of grid (degrees)
      real(r8), intent(in) :: edgew      ! western edge of grid (degrees)

! local variables:
      real(r8) :: dx                     ! land model cell width
      real(r8) :: dy                     ! land model cell length
      integer :: i,j                     ! indices
      logical :: center_at_dateline      ! logical for 1st grid center at dateline
      real(r8), intent(out) :: latn(lat_points)     ! grid cell latitude, northern edge (degrees)
      real(r8), intent(out) :: lonw(lon_points)     ! grid cell longitude, western edge (degrees)
      real(r8), intent(out) :: lats(lat_points)     ! grid cell latitude, northern edge (degrees)
      real(r8), intent(out) :: lone(lon_points)     ! grid cell longitude, western edge (degrees)

      integer :: iunit 
      character(len=256) lndname
      real(r8) :: latixy(lon_points,lat_points)  ! latitude (degrees) of the center of the model grids
      real(r8) :: longxy(lon_points,lat_points)  ! longitude (degrees) of the center of the model grids
      real(r8) :: area(lon_points,lat_points)    ! grid cell area
 
      real(r8) dx2

    ! ------------------------------------------------------------------
    ! determine grid longitudes and latitudes in increments of dx and dy
    ! ------------------------------------------------------------------

      dx = (edgee-edgew)/lon_points
      dy = (edgen-edges)/lat_points

      do j = 1, lat_points
         do i = 1, lon_points
            latixy(i,j) = edgen - (2*j-1)*dy/2.
            longxy(i,j) = edgew + (2*i-1)*dx/2.
         end do
      end do
      center_at_dateline = .false.
    
    ! ------------------------------------------
    ! determine the model grid edges 
    ! ------------------------------------------
    ! latitudes
      latn(1) = edgen
      lats(lat_points) = edges
      do j = 2, lat_points
         latn(j) = (latixy(1,j-1) + latixy(1,j)) / 2.
         lats(j-1) = (latixy(1,j-1) + latixy(1,j)) / 2.
      end do

    ! longitudes
    ! western edge of first grid cell -- since grid starts with western
    ! edge on Dateline, lonw(1,j)=-180.
      lonw(1) = edgew
      lone(lon_points) = edgee
      if(center_at_dateline)then
         if(lon_points > 1) then  ! not single-point
            lonw(2) = (longxy(1,1)+longxy(2,1))/2.
            lone(1) = (longxy(1,1)+longxy(2,1))/2.
         endif
      else
         dx2 = longxy(1,1) - lonw(1)
         if(lon_points > 1) then ! not single-point
            lonw(2) = longxy(1,1) + dx2
            lone(1) = longxy(1,1) + dx2
         endif   
      endif

      do i = 3, lon_points
         dx2 = longxy(i-1,1) - lonw(i-1)
         lonw(i) = longxy(i-1,1) + dx2
         lone(i-1) = longxy(i-1,1) + dx2
      end do

    ! ------------------------------------------
    ! get the model grid cell areas
    ! ------------------------------------------

#if(defined USE_POINT_DATA)
      area(lon_points,lat_points) = 1.
#else
      call cellarea (lat_points,lon_points,latn,lats,lonw,lone,&
                     edgen,edgee,edges,edgew,area)
#endif    
    
    ! ------------------------------------------
    ! write out the coordinate of the center of the model grids and area of grid cells
    ! ------------------------------------------
      iunit = 100
      lndname = trim(dir_model_landdata)//'model_lonlat_gridcell.bin'
      print*,trim(lndname)
      open(iunit,file=trim(lndname),form='unformatted',status='unknown')
      write(iunit) latixy
      write(iunit) longxy
      write(iunit) area
      close(iunit)

#if (defined CLMDEBUG)
      print*,'latixy =', minval(latixy, mask = latixy .gt. -1.0e30), &
                         maxval(latixy, mask = latixy .gt. -1.0e30)
      print*,'lpngxy =', minval(longxy, mask = longxy .gt. -1.0e30), &
                         maxval(longxy, mask = longxy .gt. -1.0e30)
#endif

END SUBROUTINE crgrid

! ----------------------------------------------------------------------
! EOP

