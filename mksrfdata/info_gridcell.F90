
#include <define.h>

SUBROUTINE info_gridcell ( lon_points,lat_points,edgen,edgee,edges,edgew, &
                           nrow_start,nrow_end,ncol_start,ncol_end, &
                           nx_fine_gridcell,ny_fine_gridcell,area_fine_gridcell,&
                           latn,lats,lonw,lone,&
                           sinn,sins,lonw_rad,lone_rad,sinn_i,sins_i,lonw_rad_i,lone_rad_i,&
                           READ_row_UB,READ_row_LB,READ_col_UB,READ_col_LB )
! ----------------------------------------------------------------------
! Creates land gricell infomation of original "raw" data files -
!     data with 30 x 30 arc seconds resolution
!     and model gridcell
!
! Created by Yongjiu Dai, 02/2014
!  
! ________________
! REVISION HISTORY:
!   /07/2014, Siguang Zhu & Xiangxiang Zhang: calculate the start/end
!               row/column number of fine resolution grid for user defined
!               model grid (low resolution).
!
! ----------------------------------------------------------------------
USE precision
USE GlobalVars

IMPLICIT NONE
! arguments:

      integer,  intent(in) :: lon_points ! number of model longitude grid points
      integer,  intent(in) :: lat_points ! model  of model latitude grid points
      real(r8), intent(in) :: edgen ! northern edge of grid (degrees)
      real(r8), intent(in) :: edgee ! eastern edge of grid (degrees)
      real(r8), intent(in) :: edges ! southern edge of grid (degrees)
      real(r8), intent(in) :: edgew ! western edge of grid (degrees)

      real(r8), parameter :: edgen_i = 90.   ! northern edge of grid (deg)
      real(r8), parameter :: edges_i = -90.  ! southern edge of grid (deg)
      real(r8), parameter :: edgew_i = -180. ! western edge of grid (deg)
      real(r8), parameter :: edgee_i = 180.  ! eastern edge of grid (deg)

      integer, intent(out) :: nrow_start
      integer, intent(out) :: nrow_end
      integer, intent(out) :: ncol_start
      integer, intent(out) :: ncol_end
      integer, intent(out) :: nx_fine_gridcell
      integer, intent(out) :: ny_fine_gridcell

      real(r8), intent(in) :: latn(lat_points)  ! grid cell latitude, northern edge (deg)
      real(r8), intent(in) :: lats(lat_points)  ! grid cell latitude, sourthern edge (deg)
      real(r8), intent(in) :: lonw(lon_points)  ! grid cell longitude, western edge (deg)
      real(r8), intent(in) :: lone(lon_points)  ! grid cell longitude, eastern edge (deg)
      real(r8) :: latn_i(nlat)                  ! fine grid cell latitude, northern edge (deg)
      real(r8) :: lats_i(nlat)                  ! fine grid cell latitude, sourthern edge (deg)
      real(r8) :: lonw_i(nlon)                  ! fine grid cell longitude, western edge (deg)
      real(r8) :: lone_i(nlon)                  ! fine grid cell longitude, eastern edge (deg)
      real(r8), intent(out) :: sinn(lat_points) ! grid cell latitude, northern edge(sin)  
      real(r8), intent(out) :: sins(lat_points) ! grid cell latitude, northern edge(sin)
      real(r8), intent(out) :: lonw_rad(lon_points)   ! grid cell longitude, western edge (radian)
      real(r8), intent(out) :: lone_rad(lon_points)   ! grid cell longitude, eastern edge (radian)
      real(r8), intent(out) :: sinn_i(nlat)           ! fine grid cell latitude, northern edge(sin)
      real(r8), intent(out) :: sins_i(nlat)           ! fine grid cell latitude, northern edge(sin)
      real(r8), intent(out) :: lonw_rad_i(nlon)       ! fine grid cell longitude, western edge (radian)
      real(r8), intent(out) :: lone_rad_i(nlon)       ! fine grid cell longitude, eastern edge (radian)
      integer, intent(out) :: READ_row_UB(lat_points) ! north boundary index for fine gird cell
      integer, intent(out) :: READ_col_UB(lon_points) ! west boundary index for fine gird cell  
      integer, intent(out) :: READ_row_LB(lat_points) ! south boundary index for fine gird cell
      integer, intent(out) :: READ_col_LB(lon_points) ! east boundary index for fine gird cell

      real(r8), intent(out) :: area_fine_gridcell(nlon,nlat)  ! cell area (km**2)

! ---------------------------------------------------------------
      real(r8) dx
      real(r8) dy
      real(r8), allocatable :: lat_i(:)
      real(r8), allocatable :: lon_i(:) 
 
      real(r8) deg2rad                          ! pi/180
      real(r8) pi                               ! 3.14159265358979323846

      integer i, j
      integer tmp_lat(lat_points) 
      integer tmp_lon(lon_points) 

      integer, external :: nearest_boundary

! ---------------------------------------------------------------
! define the gridcell edges and area of the RAW DATA at 30 arc-seconds resolution
! ---------------------------------------------------------------
      allocate ( lat_i(nlat) )
      allocate ( lon_i(nlon) )

      pi = 4.*atan(1.)
      deg2rad = pi/180.
      dx = (edgee_i-edgew_i)/nlon   ! = 1./120.
      dy = (edgen_i-edges_i)/nlat   ! = 1./120.

#ifdef OPENMP
print *, 'OPENMP enabled, threads num = ', OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(j)
#endif
      do j = 1, nlat
         lat_i(j)  = edgen_i - (2*j-1)*dy/2.
         latn_i(j) = edgen_i - (j-1)*dy
         lats_i(j) = edgen_i - j*dy
         sinn_i(j) = sin(latn_i(j)*deg2rad)
         sins_i(j) = sin(lats_i(j)*deg2rad)
      enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif
      
#ifdef OPENMP
print *, 'OPENMP enabled, threads num = ', OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(i)
#endif
      do i = 1, nlon
         lon_i(i)  = edgew_i + (2*i-1)*dx/2.
         lonw_i(i) = edgew_i + (i-1)*dx
         lone_i(i) = edgew_i + i*dx
         lonw_rad_i(i) =  lonw_i(i)*deg2rad
         lone_rad_i(i) =  lone_i(i)*deg2rad
      enddo
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif

      do j = 1,lat_points
         sinn(j) = sin(latn(j)*deg2rad)
         sins(j) = sin(lats(j)*deg2rad)
      enddo

      do i = 1,lon_points
         lonw_rad(i) =  lonw(i)*deg2rad
         lone_rad(i) =  lone(i)*deg2rad
      enddo


! calculate the area of the grid-cell of RAW DATA (fine-gridcell)
      call cellarea(nlat,nlon,latn_i,lats_i,lonw_i,lone_i,&
                    edgen_i,edgee_i,edges_i,edgew_i,&
                    area_fine_gridcell)

! --------------------------------------------------------
! define the starting and ending points and the numbers of 
! the RAW DATA fine gridcell of the model grids
! --------------------------------------------------------

#if(defined USE_POINT_DATA)

! modified by yuan, 07/15/2016
! bug may exist, nint -> int
      !nrow_start = nint((90.-edgen)/dy) + 1
      nrow_start = int((90.-edgen)/dy) + 1
      nrow_end   = nrow_start
      if(nrow_start.lt.1) nrow_start = 1
      if(nrow_end.gt.nlat) nrow_end = nlat

      !ncol_start = nint((180.+edgew)/dx) + 1
      ncol_start = int((180.+edgew)/dx) + 1
      ncol_end   = ncol_start 
      if(ncol_start.lt.1) ncol_start = 1
      if(ncol_end.gt.nlon) ncol_end = nlon

      nx_fine_gridcell = 1
      ny_fine_gridcell = 1

#elif(defined USER_GRID)
        
      do i = 1,lat_points
         READ_row_UB(i) = nearest_boundary(latn(i),lat_i,nlat)
      enddo
      
      do i = 1,lat_points
         READ_row_LB(i) = nearest_boundary(lats(i),lat_i,nlat)
      enddo
    
      do i = 1,lon_points
         READ_col_UB(i) = nearest_boundary(lonw(i),lon_i,nlon)
      enddo
      
      do i = 1,lon_points
         READ_col_LB(i) = nearest_boundary(lone(i),lon_i,nlon)
      enddo
      
      ! find the max interval of fine grids
      do i = 1,lat_points
         tmp_lat(i) = READ_row_LB(i) - READ_row_UB(i)
      end do
      
      do i = 1,lon_points
         tmp_lon(i) = READ_col_LB(i) - READ_col_UB(i)
         if(READ_col_LB(i) < READ_col_UB(i)) then  ! gridcell across the dateline
            tmp_lon(i) = tmp_lon(i) + nlon
         endif
      end do

      nx_fine_gridcell = maxval(tmp_lon) + 1    ! max interval of fine grids in lon gird 
      ny_fine_gridcell = maxval(tmp_lat) + 1    ! max interval of fine grids in lat gird

      nrow_start = READ_row_UB(1)
      nrow_end   = READ_row_LB(lat_points)
      ncol_start = READ_col_UB(1)
      ncol_end   = READ_col_LB(lon_points)
#else
      nrow_start = nint((90.-edgen)/dy) + 1
      nrow_end   = nint((90.-edges)/dy)
      if(nrow_start.lt.1) nrow_start = 1
      if(nrow_end.gt.nlat) nrow_end = nlat

      ncol_start = nint((180.+edgew)/dx) + 1
      ncol_end   = nint((180.+edgee)/dx)
      if(ncol_start.lt.1) ncol_start = 1
      if(ncol_end.gt.nlon) ncol_end = nlon

      nx_fine_gridcell = nint( (edgee-edgew)/lon_points/dx )
      ny_fine_gridcell = nint( (edgen-edges)/lat_points/dy )
#endif

!----------------add end-------------------------------------
      deallocate ( lat_i )
      deallocate ( lon_i )

END SUBROUTINE info_gridcell


integer function nearest_boundary(degree,degree_i,dim_i)

!=======================================================================
! find the boundary index for aggregation
! edit by zsg 20140823
!=======================================================================
   use precision
   implicit none
   
   integer  j 
   integer  dim_i
   integer  a_tmp
   integer  minloc_tmp(1)
   real(r8) degree   
   real(r8) degree_i(dim_i)
   real(r8) diff_degree(dim_i)

   real(r8), parameter :: edgen_i = 90.   ! northern edge of grid (deg)
   real(r8), parameter :: edgew_i = -180. ! western edge of grid (deg)
   
   
   ! find the naerest fine grids boundary's index
   diff_degree(:) = abs(degree - degree_i(:))
   
   if(degree_i(1) == edgen_i .or. degree_i(1) == edgew_i) then  ! for upper boundary 
      ! adjust the index to make sure all related fine grids are 
      ! included if the edge of coarse grid coincide with edge of fine gird 
      diff_degree(1:dim_i) = diff_degree(dim_i:1:-1)
      
      minloc_tmp = minloc(diff_degree)
      a_tmp = minloc_tmp(1)
      nearest_boundary = dim_i - a_tmp + 1
   else  ! for lower boundary of latitude
      minloc_tmp = minloc(diff_degree)
      a_tmp = minloc_tmp(1)
      nearest_boundary = a_tmp 
   endif

end function nearest_boundary 
!-----------------------------------------------------------------------
!EOP
