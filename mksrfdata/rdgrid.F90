
#include <define.h>

SUBROUTINE rdgrid(dir_model_landdata,edgen,edgee,edges,edgew,&
                  lon_points,lat_points,latn,lats,lonw,lone)
! ----------------------------------------------------------------------
! read land model grid from a user defined domain file
! region (global) latitude grid goes from:
!                          NORTHERN edge (POLE) to SOUTHERN edge (POLE)
! region (global) longitude grid starts at:
!                          WESTERN edge 
!                         (DATELINE with western edge)
! surface grid edges -- grids do not have to be global. 
! to allow this, grids must define the north, east, south, west edges:
!    edgen: northern edge of grid : > -90 and <= 90 (degrees)
!    edgee: eastern edge of grid  : > western edge and <= 180
!    edges: southern edge of grid : >= -90  and <  90
!    edgew: western edge of grid  : >= -180 and < 180
!
! Created by Yongjiu Dai, /02/2014
! ________________
! REVISION HISTORY:
!   /07/2014, Siguang Zhu & Xiangxiang Zhang: read the grid info from a 
!               user defined domain file.
!
! ----------------------------------------------------------------------
use precision
use netcdf

IMPLICIT NONE

! arguments:
      character(len=256), intent(in) :: dir_model_landdata
      integer, intent(in) :: lon_points ! number of input data longitudes
      integer, intent(in) :: lat_points ! number of input data latitudes
      real(r8), intent(in) :: edgen     ! northern edge of grid (degrees)
      real(r8), intent(in) :: edges     ! southern edge of grid (degrees)
      real(r8), intent(in) :: edgee     ! eastern edge of grid (degrees)
      real(r8), intent(in) :: edgew     ! western edge of grid (degrees)

! local variables:
      integer  :: iunit                 ! logical unit number of grid file
      integer  :: i,j                   ! indices
      real(r8) :: lon(lon_points)       ! read-in longitude array (full grid)
      real(r8) :: lat(lat_points)       ! read-in latitude array (full grid)
      logical  :: center_at_dateline    ! input grid: cell area
      
      real(r8),intent(out) :: latn(lat_points)    ! grid cell latitude, northern edge (deg)
      real(r8),intent(out) :: lats(lat_points)    ! grid cell latitude, southern edge (deg)
      real(r8),intent(out) :: lonw(lon_points)    ! grid cell longitude, western edge (deg)
      real(r8),intent(out) :: lone(lon_points)    ! grid cell longitude, eastern edge (deg)

      character(len=256) lndname
      real(r8) :: latixy(lon_points,lat_points) ! latitude (deg) of the center of the model grid
      real(r8) :: longxy(lon_points,lat_points) ! longitude (deg) of the center of the model grid
      real(r8) :: area(lon_points,lat_points)   ! grid cell area

      real(r8) :: dx2
      integer fid, vid

! ----------------------------------------------------------------------
      
      iunit = 100
#ifdef USER_GRID
      lndname = USER_GRID
      print*,lndname
#endif

      print*, 'Attempting to read land grid data .....'
      
      call sanity( nf90_open(path=lndname, mode=nf90_nowrite, ncid=fid) )
      call sanity( nf90_inq_varid(fid, "latixy", vid) )
      call sanity( nf90_get_var(fid, vid, latixy(:,:),start=(/1,1/), count=(/lon_points,lat_points/)) )
      call sanity( nf90_inq_varid(fid, "longxy", vid) )
      call sanity( nf90_get_var(fid, vid, longxy(:,:),start=(/1,1/), count=(/lon_points,lat_points/)) )
      call sanity( nf90_inq_varid(fid, "latn", vid) )
      call sanity( nf90_get_var(fid, vid, latn) )
      call sanity( nf90_inq_varid(fid, "lats", vid) )
      call sanity( nf90_get_var(fid, vid, lats) )
      call sanity( nf90_inq_varid(fid, "lonw", vid) )
      call sanity( nf90_get_var(fid, vid, lonw) ) 
      call sanity( nf90_inq_varid(fid, "lone", vid) )
      call sanity( nf90_get_var(fid, vid, lone) ) 

    ! --------------------------------------
    ! determine the longitude and latitudes of the center of the model grids
    ! --------------------------------------

    ! determine if grid has pole points - if so, make sure that north pole
    ! is non-land and south pole is land
      if(abs((latixy(1,lat_points)-90.)) < 1.e-6) then
         print*, 'GRID: model has pole_points' 
      else
         print*, 'GRID: model does not have pole_points' 
      endif

      center_at_dateline = .false.

    ! for grids starting at greenwich and centered on Greenwich'
      if(abs(longxy(1,1)+180.).lt.1.e-6)then
         center_at_dateline = .true.
      endif

    ! ------------------------------------------
    ! calculate the model grid cell areas
    ! ------------------------------------------
      call cellarea (lat_points,lon_points,latn,lats,lonw,lone,&
                     edgen,edgee,edges,edgew,area)

    ! ------------------------------------------
    ! write out the model grid coordinate and area
    ! ------------------------------------------
      iunit = 100
      lndname = trim(dir_model_landdata)//'model_lonlat_gridcell.bin'
      print*,lndname
      open(iunit,file=trim(lndname),form='unformatted',status='unknown')
      write(iunit) latixy
      write(iunit) longxy
      write(iunit) area
      close(iunit)

END SUBROUTINE rdgrid

! nc file open check
! ------------------------------------------------------------
SUBROUTINE sanity(ret)
      use netcdf
      implicit none
      integer, intent(in) :: ret

      if (ret .ne. nf90_noerr) then
         write(6, *) trim(nf90_strerror(ret)); stop
      end if

END SUBROUTINE sanity
! ----------------------------------------------------------------------
! EOP
