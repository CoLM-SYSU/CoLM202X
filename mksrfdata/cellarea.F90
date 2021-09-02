  
#include <define.h>

subroutine cellarea(lat_points,lon_points,latn,lats,lonw,lone,&
                                 edgen,edgee,edges,edgew,area)

! ----------------------------------------------------------------------
! Area of grid cells (square kilometers) - regional grid
! (can become global as special case)
! ----------------------------------------------------------------------
use precision

IMPLICIT NONE

! arguments
      integer, intent(in) :: lat_points     ! number of latitude points
      integer, intent(in) :: lon_points     ! number of longitude points

      real(r8), intent(in) :: latn(lat_points)! grid cell latitude, northern edge (deg)
      real(r8), intent(in) :: lats(lat_points)! grid cell latitude, sourthern edge (deg)
      real(r8), intent(in) :: lonw(lon_points)! grid cell longitude, western edge (deg)
      real(r8), intent(in) :: lone(lon_points)! grid cell longitude, eastern edge (deg)

      real(r8), intent(in) :: edgen             ! northern edge of grid (deg)
      real(r8), intent(in) :: edges             ! southern edge of grid (deg)
      real(r8), intent(in) :: edgew             ! western edge of grid (deg)
      real(r8), intent(in) :: edgee             ! eastern edge of grid (deg)
      real(r8), intent(out):: area(lon_points,lat_points) ! cell area (km**2)

! local variables
      integer i,j                               ! indices
      real(r8) deg2rad                          ! pi/180
      real(r8) global                           ! summed area
      real(r8) dx                               ! cell width: E-W
      real(r8) dy                               ! cell width: N-S
      real(r8) pi                               ! 3.14159265358979323846
      real(r8) re                               ! radius of earth (km)
      real(r8) error                            ! true area for error check
! -----------------------------------------------------------------------
      re = 6.37122e6 * 0.001                    ! kilometer
      pi = 4.*atan(1.)
      deg2rad = pi/180.
      global = 0.

#ifdef OPENMP
print *, 'OPENMP enabled, threads num = ', OPENMP
!$OMP PARALLEL DO NUM_THREADS(OPENMP) &
!$OMP PRIVATE(i,j,dx,dy)
#endif
      do j = 1, lat_points
         do i = 1, lon_points
            if(lone(i)*lonw(i)<0 .and. lone(i)<lonw(i))then   ! west edge is more western than data line
               dx = (lone(i)-lonw(i)+360.0)*deg2rad
            else
               dx = (lone(i)-lonw(i))*deg2rad
            endif
            if(latn(j)>lats(j)) then          ! north to south grid
               dy = sin(latn(j)*deg2rad) - sin(lats(j)*deg2rad)
            else                              ! south to north grid
               dy = sin(lats(j)*deg2rad) - sin(latn(j)*deg2rad)
            end if
            area(i,j) = dx*dy*re*re
         end do
      end do
#ifdef OPENMP
!$OMP END PARALLEL DO
#endif
            
      global = sum(area(:,:))

    ! make sure total area from grid cells is same as area of grid
    ! as defined by its edges
      dx = (edgee - edgew) * deg2rad
      dy = sin(edgen*deg2rad) - sin(edges*deg2rad)
      error = dx*dy*re*re
      if(abs(global-error)/error > 0.00001) then
         print*, 'CELLAREA error: correct area is ',error, &
              ' but summed area of grid cells is ', global
!        call abort
      end if
      return
  end subroutine cellarea
  
  
  function find_min_area(edgee_rad_coarse_gridcell,edgew_rad_coarse_gridcell,edgee_rad_fine_gridcell, & 
                         edgew_rad_fine_gridcell,edgen_sin_coarse_gridcell,edges_sin_coarse_gridcell, & 
                         edgen_sin_fine_gridcell,edges_sin_fine_gridcell)
!-------------------------------------------------------------------------------
! FUNCTION to find out the minimum area for weighting. 
!-------------------------------------------------------------------------------

  use precision
  implicit none

  real(r8), intent(in) :: edgee_rad_coarse_gridcell  ! grid cell longitude, western edge (radian)
  real(r8), intent(in) :: edgew_rad_coarse_gridcell  ! grid cell longitude, eastern edge (radian)
  real(r8), intent(in) :: edgee_rad_fine_gridcell    ! fine grid cell longitude, western edge (radian)
  real(r8), intent(in) :: edgew_rad_fine_gridcell    ! fine grid cell longitude, eastern edge (radian)
  real(r8), intent(in) :: edgen_sin_coarse_gridcell  ! grid cell latitude, northern edge(sin)  
  real(r8), intent(in) :: edges_sin_coarse_gridcell  ! grid cell latitude, northern edge(sin)
  real(r8), intent(in) :: edgen_sin_fine_gridcell    ! fine grid cell latitude, northern edge(sin)
  real(r8), intent(in) :: edges_sin_fine_gridcell    ! fine grid cell latitude, northern edge(sin)
  real(r8) :: find_min_area

!  --- Local variables ---
  real(r8) diff_lon1  ! interval of coarse grid cell                                                                     
  real(r8) diff_lon2  ! interval of fine grid cell
  real(r8) diff_lon3  ! interval between eastern edge of fine grid cell and western edge of coarse grid cell
  real(r8) diff_lon4  ! interval between eastern edge of coarse grid cell and western edge of fine grid cell
  real(r8) dx         ! cell width: E-W
  real(r8) dy         ! cell width: N-S
  real(r8) pi         ! 3.14159265358979323846
  real(r8) re         ! radius of earth (km)

! -------------------------------------------------------------------------------

     re = 6.37122e6 * 0.001                    ! kilometer
     pi = 4.*atan(1.)

     diff_lon1 = edgee_rad_coarse_gridcell - edgew_rad_coarse_gridcell
     diff_lon2 = edgee_rad_fine_gridcell - edgew_rad_fine_gridcell
     diff_lon3 = edgee_rad_fine_gridcell - edgew_rad_coarse_gridcell
     diff_lon4 = edgee_rad_coarse_gridcell - edgew_rad_fine_gridcell

     if(edgee_rad_coarse_gridcell < edgew_rad_coarse_gridcell) then   
        diff_lon1 = 2 * pi + diff_lon1
        if(edgee_rad_fine_gridcell * edgew_rad_coarse_gridcell < 0 .and. & 
           edgee_rad_fine_gridcell < edgew_rad_coarse_gridcell) then
           diff_lon3 = diff_lon3 + 2 * pi
        endif

        if(edgee_rad_coarse_gridcell * edgew_rad_fine_gridcell < 0 .and. &
           edgee_rad_coarse_gridcell < edgew_rad_fine_gridcell) then
           diff_lon4 = diff_lon4 + 2 * pi
        endif
     endif

     dx = min(min(min(diff_lon1,diff_lon2),diff_lon3),diff_lon4)
     dy = min(min(min(edgen_sin_coarse_gridcell - edges_sin_coarse_gridcell,edgen_sin_fine_gridcell &
          - edges_sin_fine_gridcell ), edgen_sin_fine_gridcell - edges_sin_coarse_gridcell), &
          edgen_sin_coarse_gridcell - edges_sin_fine_gridcell)
      
     find_min_area = dx * dy * re * re

   end function find_min_area
! -----------------------------------------------------------------------
! EOP
