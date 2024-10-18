#include <define.h>

! ---------------------------------- code history -------------------------------------
! Copyright (C) 2020, Eawag
!     Developed by:  Group of Applied System Analysis
!                    Dept. of Surface Waters - Research and Management
!                    Eawag - Swiss Federal institute of Aquatic Science and Technology
!
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General PUBLIC License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General PUBLIC License for more details.
!
!     You should have received a copy of the GNU General PUBLIC License
!     along with this program.  IF not, see <http://www.gnu.org/licenses/>.
!
!
! Description:
!     (add by WMEJ)
!     Simstrat a physical 1D model for lakes and reservoirs.
!     Compared to the other 1D models, Simstrat requires two additional predictive equations
!     to calculate the turbulent kinetic energy k and dissipation rate ε for turbulent mixing
!     to parameterize the vertical transport of energy inside the lake. 
!     The equations for calculating the energy in the lake take into account 
!     the depth dependence of the horizontal cross-sectional area of the lake, 
!     while the momentum equations take into account the Coriolis force and vertical viscosity. 
!     It should be noted that the Simstrat does not take into account sediment at the lake bottom. 
!     Typically, energy exchange at the bottom of lakes is achieved through bottom flux, which is set as a constant(lake_table).
! 
! Simstrat's lake grid discretization is as follows (add by WMEJ):
!     
!     [  z ] - layer node depth (m)
!     [ zi ] - layer interface depth (m)
!     [  h ] - layer thickness (m)
!     [ Az ] - layer area (m^2) 
!
!     centres                                                         faces
!     T, U, V                                                         K, eps
!                                      |
!                 ____________Az+2_____|________#___________ -------- zi = 4
!                  \                   |        |h3       /
!     z = 3 --------\------------------* z + 1  |        /
!                    \________Az+1_____|________#_______/------------ zi = 3
!                     \                |        |h2    /
!     z = 2 -----------\---------------* z      |     /
!                       \_____Az_______|________#____/--------------- zi = 2
!                        \             |        |h1 /
!     z = 1 --------------\------------* z - 1  |  /
!                          \__Az-1_____|________#_/------------------ zi = 1
!                                 
!
!    !**************************************************************************!
!    !* Simstrat's lake depth is measured from the bottom of the lake          *!
!    !* The counting of the lake layer also starts from the bottom of the lake *!
!    !**************************************************************************!
!
!
! -WMEJ in order to ensure that Simstrat can continue to receive updates from Eawag in the future,
!       I did not modify the structure of Simstrat. 
!       IF you are a beginner, the following information should be helpful to you.
!   
!       1) strat_kinds             : Type definitions for some global interfaces
!       2) strat_utilities         : Generic strat_utilities that are used throughout the code
!       3) strat_grid              : modify the grid, update the grid, and interpolate the grid
!       4) strat_simdata           : Data structure definitions for simulation data
!       5) strat_solver            : CONTAINS implementation of a tridiagonal matrix solver
!       6) strat_discretization    : Interface and implementation for discretization methods
!       7) strat_statevar          : Abstract base class for simulation variables (K, EPS, T, S etc)
!       8) strat_turbulence        : Turbulence MODULE
!       9) strat_absorption        : [USELESS] reads and processes absorption input file  
!      10) strat_advection         : [NOT USED NOW] Based on already read in/outflows, calculates Advection
!      11) strat_lateral           : [NOT USED NOW] Reads and processes inflows/outflows such that
!      12) strat_stability         : CONTAINS methods to update cmue_cn /qe and NN
!      13) strat_keps              : ModelVariable implementation for K and EPS variable
!      14) strat_transport         : Implementation of a statevar for a generic transport variable
!      15) strat_forcing           : calculates the lake surface fluxes
!      16) strat_ice               : Ice MODULE
!      17) strat_temp              : Implementation of a statevar for a temperature variable
!      18) strat_inputfile         : set configuration and initial conditions
!      19) MOD_Lake_Simstrat : Interface for the Simstrat model
!       
!       
!      Source : https://github.com/Eawag-AppliedSystemAnalysis/Simstrat
!      More information should be consulted in the Simstrat documentation.
!
!      *Please focus on the following variables when testing Simstrat:
!           - T, U, V, S, K, eps, E_Seiche, num, nuh, NN, ice_h, snow_h
!
!
! Revisions: 
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024, Modified for new Lake Driver in CWRF 
!
! -------------------------------------------------------------------------------------






!<    +---------------------------------------------------------------+
!     |  Type definitions for some global interfaces
!<    +---------------------------------------------------------------+
MODULE  strat_kinds
   
   IMPLICIT NONE

   integer, parameter, public  :: RK = selected_real_kind(12) ! RK = RK

!-----------------------------------------------------------
PRIVATE
!-----------------------------------------------------------
   ! Common constants
   character(len=5), parameter, public :: version = '2.4.2'

   ! Common Types
   type, abstract, public :: LinSysSolver
   CONTAINS
      procedure(generic_linsyssolver_solve), deferred, nopass :: solve
   END type

   type, abstract, public :: SimModule
   CONTAINS
      procedure(generic_simmodule_update), deferred, nopass :: update
   END type

   type, abstract, extends(SimModule), public :: StateVariable
   CONTAINS
      procedure(generic_statevariable_solve), deferred, nopass :: solve
   END type

!-----------------------------------------------------------
CONTAINS
!-----------------------------------------------------------

   SUBROUTINE generic_linsyssolver_solve(ld, md, ud, rhs, x, ubnd)
      ! Arguments
      real(RK), intent(in) :: ld(:), md(:), ud(:) ! Diagonals (A)
      real(RK), intent(in) :: rhs(:) ! right-hand side (b)
      real(RK), intent(out) :: x(:) ! solution (x)
      integer, intent(in) :: ubnd
      integer :: i
   
      DO i = 1, ubnd
         x(i) = 9999.00_RK
      ENDDO

   END SUBROUTINE

   SUBROUTINE generic_statevariable_solve()
   END SUBROUTINE

   SUBROUTINE generic_simmodule_update()
   END SUBROUTINE

END MODULE  strat_kinds





!<    +---------------------------------------------------------------+
!     | Generic strat_utilities that are used throughout the code
!<    +---------------------------------------------------------------+
MODULE  strat_utilities

   USE strat_kinds, only: RK
   IMPLICIT NONE

   type, public :: string
      character(len=:), allocatable :: str
   END type

   interface toStr
      MODULE  procedure str_int, str_real
   END interface

!-----------------------------------------------------------
CONTAINS
!-----------------------------------------------------------
   !> Interpolation of yi on grid zi (based on the given y on grid z)
   SUBROUTINE Interp(z, y, num_z, zi, yi, num_zi)
        IMPLICIT NONE 
        real(RK), DIMENSION(:), intent(in) :: z, y, zi
        real(RK), DIMENSION(:), intent(out) :: yi
        integer, intent(in) :: num_z, num_zi
        integer :: posk1, posk2, posi, i

        IF (num_z == 1) THEN ! only one value
            yi(1:num_zi) = y(1)
            RETURN
        ENDIF

        ! Assign closest value IF out of given grid
        posk1 = 1
        DO WHILE (zi(posk1) <= z(1))
            yi(posk1) = y(1)
            posk1 = posk1 + 1
        ENDDO
        posk2 = num_zi
        DO WHILE (zi(posk2) >= z(num_z))
            yi(posk2) = y(num_z)
            posk2 = posk2 - 1
        ENDDO

        ! Linear interpolation
        posi = 1
        DO i = posk1, posk2
            DO WHILE (zi(i) > z(posi + 1))
                posi = posi + 1
            ENDDO
            yi(i) = y(posi) + (zi(i) - z(posi))*(y(posi + 1) - y(posi))/(z(posi + 1) - z(posi))
        ENDDO

        RETURN
   END SUBROUTINE Interp


   ! Interpolation of yi on grid zi (based on the given y on grid z)
   pure SUBROUTINE Interp_nan(z, y, num_z, zi, yi, num_zi)
      USE strat_kinds, only: RK
      ! USE, intrinsic :: iso_fortran_env
      USE, intrinsic :: ieee_arithmetic
      IMPLICIT NONE

      real(RK), DIMENSION(:), intent(in) :: z, y, zi
      real(RK), DIMENSION(:), intent(out) :: yi
      integer, intent(in) :: num_z, num_zi

      integer posk1, posk2, posi, i

      ! Assign NaN IF out of given grid
      posk1 = 1
      DO WHILE (zi(posk1) < z(1))
         yi(posk1) = 0.0_RK
         yi(posk1) = ieee_value(yi(posk1), ieee_quiet_nan) ! NaN
         posk1 = posk1 + 1
      ENDDO
      posk2 = num_zi
      DO WHILE (zi(posk2) > z(num_z))
         yi(posk2) = 0.0_RK
         yi(posk2) = ieee_value(yi(posk2), ieee_quiet_nan) ! NaN
         posk2 = posk2 - 1
      ENDDO

      ! Linear interpolation
      posi = 1
      DO i = posk1, posk2
         DO WHILE (zi(i) > z(posi + 1))
            posi = posi + 1
         ENDDO
         yi(i) = y(posi) + (zi(i) - z(posi))*(y(posi + 1) - y(posi))/(z(posi + 1) - z(posi))
      ENDDO

      RETURN
   END SUBROUTINE Interp_nan


   !! Integrate discrete function y[x] using the trapezoidal rule
   SUBROUTINE Integrate(x, y, inty, num)
      IMPLICIT NONE

      integer :: num, i
      real(RK), DIMENSION(:), intent(in) :: x, y
      real(RK), DIMENSION(:), intent(inout) :: inty

      inty(1) = 0
      DO i = 2, num
         inty(i) = inty(i - 1) + 0.5_RK*(x(i) - x(i - 1))*(y(i) + y(i - 1))
      ENDDO
      RETURN
   END SUBROUTINE Integrate


   ! Assign nan to values out of current grid
   pure SUBROUTINE Assign_nan(y, ubnd, ubnd_grid)
      USE strat_kinds, only: RK
      ! USE, intrinsic :: iso_fortran_env
      USE, intrinsic :: ieee_arithmetic
      !   IMPLICIT NONE

      real(RK), DIMENSION(:), intent(out) :: y
      integer, intent(in) :: ubnd, ubnd_grid

      integer :: i

      ! Assign NaN IF out of given grid
      i = ubnd_grid
      DO WHILE (i > ubnd)
         y(i) = 0.0_RK
         y(i) = ieee_value(y(i), ieee_quiet_nan) ! NaN
         i = i - 1
      ENDDO

      RETURN
   END SUBROUTINE Assign_nan


   pure function linspace(x0, xend, n, endpoint) result(x)
      IMPLICIT NONE

      integer, intent(in) :: n
      real(RK), intent(in) :: x0, xend
      logical, optional, intent(in) :: endpoint
      real(RK), DIMENSION(n) :: x

      real(RK) :: dx
      real(RK) :: denom
      integer :: i

      denom = real(n - 1, RK)
      IF (present(endpoint) .and. (.not. endpoint)) THEN
         denom = real(n, RK)
      ENDIF

      dx = (xend - x0)/denom

      x = [(real(i, RK)*dx + x0, i=0, n - 1)]

      RETURN
   END function linspace


   pure SUBROUTINE diff(d, a, N)
      IMPLICIT NONE

      integer, intent(in) :: N
      real(RK), DIMENSION(N - 1), intent(out) :: d
      real(RK), DIMENSION(N), intent(in) :: a

      d = a(2:N) - a(1:N - 1)
      RETURN
   END SUBROUTINE diff


   SUBROUTINE check_file_exists(fname)
      IMPLICIT NONE
      character(len=*), intent(in) :: fname

      logical :: file_exists
      IF (fname == '') THEN
         CALL error('Filename is empty')
      ELSE
         inquire (file=fname, exist=file_exists)
         IF (.not. file_exists) THEN
            CALL error('File '//fname//' does not exist')
         ENDIF
      ENDIF
   END SUBROUTINE check_file_exists


   SUBROUTINE ok(message)
      IMPLICIT NONE
      character(len=*), intent(in) :: message
      write(*, *) '[OK] '//message
   END SUBROUTINE ok


   SUBROUTINE error(message)
      IMPLICIT NONE
      character(len=*), intent(in) :: message
      write(6, *) '[ERROR] '//message
      stop
   END SUBROUTINE error


   SUBROUTINE warn(message)
      IMPLICIT NONE
      character(len=*), intent(in) :: message
      write(*, *) '[WARNING] '//message
   END SUBROUTINE warn


   pure function find_index_ordered(array, target_value) result(idx)
      IMPLICIT NONE
      real(RK), DIMENSION(:), intent(in) :: array
      real(RK), intent(in) :: target_value

      integer :: idx

      DO idx = 1, size(array)
         IF (array(idx) > target_value) exit
      ENDDO
   END function


   pure function linear_interpolate(t_start, t_end, v_start, v_end, t) result(v)
      IMPLICIT NONE
      real(RK), intent(in) :: t_start, t_end, v_start, v_end, t
      real(RK) :: v

      v = v_start + t*(v_end - v_start)/(t_end - t_start)
   END function


   pure function convert2height_above_sed(z, z_zero) result(h)
      IMPLICIT NONE
      real(RK), DIMENSION(:), intent(in) :: z
      real(RK), intent(in) :: z_zero

      real(RK), DIMENSION(size(z)) :: h
      integer :: n

      n = size(z)
      h = -z_zero + z(n:1:-1)
   END function


   ! Reverse an array without allocating a second array
   SUBROUTINE reverse_in_place(in_arr)
      IMPLICIT NONE
      real(RK), intent(inout) :: in_arr(:)
      real(RK) :: temp

      integer :: first, last, i, len

      first = lbound(in_arr, dim=1)
      last = ubound(in_arr, dim=1)
      len = size(in_arr)

      ! Works for even and odd sized arrays
      !(as len/2 is always integer and not rounded, but cutoff)
      DO i = last, first + int(len/2), -1
         temp = in_arr(i)
         in_arr(i) = in_arr(len + 1 - i)
         in_arr(len + 1 - i) = temp
      ENDDO

   END SUBROUTINE


   character(len=20) function str_int(k)
      IMPLICIT NONE
      ! "Convert an integer to string."
      integer, intent(in) :: k
      write (str_int, '(a)') k
      str_int = adjustl(str_int)
   END function str_int


   character(len=20) function str_real(k)
      IMPLICIT NONE
      ! "Convert an integer to string."
      real(RK), intent(in) :: k
      write (str_real, '(a)') k
      str_real = adjustl(str_real)
   END function str_real


   character(len=20) function real_to_str(k, fmt)
      IMPLICIT NONE
      real(RK), intent(in) :: k
      character(len=*), intent(in) :: fmt
      write (real_to_str, fmt) k
      real_to_str = adjustl(real_to_str)
   END function

END MODULE  strat_utilities






!<    +---------------------------------------------------------------+
!     | Grid MODULE
!     |  - CONTAINS class to store, USE and modify the grid
!<    +---------------------------------------------------------------+
MODULE  strat_grid

   USE strat_kinds, only: RK
   USE strat_utilities
   IMPLICIT NONE

!-----------------------------------------------------------
PRIVATE
!-----------------------------------------------------------
   ! Class that holds initial configuration for grid setup
   type, public :: GridConfig
      integer :: max_length_input_data
      integer :: nz_grid
      real(RK) :: max_depth
      real(RK), DIMENSION(:), allocatable :: grid_read ! Grid definition
      real(RK), DIMENSION(:), allocatable :: A_read ! Area definition
      real(RK), DIMENSION(:), allocatable :: z_A_read ! Height of area definitions
      logical :: equidistant_grid
   END type

   ! StaggeredGrid implementation
   type, public :: StaggeredGrid
      real(RK), DIMENSION(:), allocatable :: h        ! Box height
      real(RK), DIMENSION(:), allocatable :: z_face   ! Holds z-values of faces
      real(RK), DIMENSION(:), allocatable :: z_volume ! Holds z-values of volume centers
      real(RK), DIMENSION(:), allocatable :: Az       ! Areas
      real(RK), DIMENSION(:), allocatable :: dAz      ! Difference of areas
      real(RK), DIMENSION(:), allocatable :: meanint  ! ?
      real(RK) :: volume, h_old

      !Area factors
      real(RK), DIMENSION(:), allocatable :: AreaFactor_1
      real(RK), DIMENSION(:), allocatable :: AreaFactor_2
      real(RK), DIMENSION(:), allocatable :: AreaFactor_k1
      real(RK), DIMENSION(:), allocatable :: AreaFactor_k2
      real(RK), DIMENSION(:), allocatable :: AreaFactor_eps

      integer :: nz_grid      ! Number of allocated grid cells
      integer :: nz_occupied  ! number of grid cells in USE as per current lake depth
      integer :: max_length_input_data  ! Hard limit of grid cells for reading files of unknnown length etc

      integer :: ubnd_vol, ubnd_fce, length_vol, length_fce   ! Upper and lenght for volume (vol) and face(fce) grids
      real(RK) :: z_zero
      real(RK) :: lake_level, lake_level_old
      real(RK) :: max_depth

   CONTAINS   

      !Many methods...
      ! Init methods used for setup
      procedure, pass :: init => grid_init                              ! Main initialization method
      procedure, pass :: memory_init => grid_memory_init
      procedure, pass :: init_grid_points => grid_init_grid_points
      procedure, pass :: init_z_axes => grid_init_z_axes
      procedure, pass :: init_morphology => grid_init_morphology
      procedure, pass :: init_areas => grid_init_areas

      ! Update methods, used for internal update of parameters
      procedure, pass :: update_area_factors => grid_update_area_factors  ! Called to recalculate area factors
      procedure, pass :: update_depth => grid_update_depth                ! Called for recalculating depth
      procedure, pass :: update_nz => grid_update_nz                      ! updates lower and upper bounds based on nz_occupied

      ! Interpolation methods
      procedure, pass :: interpolate_to_face => grid_interpolate_to_face    ! Interpolate quantitiy onto face grid
      procedure, pass :: interpolate_to_face_from_second => grid_interpolate_to_face_from_second  ! Interpolate quantity onto face grid, ignoring first value
      procedure, pass :: interpolate_to_vol => grid_interpolate_to_vol      ! Interpolate quantity onto volume grid
      procedure, pass :: interpolate_from_face => grid_interpolate_from_face  ! Interpolate quantity that is stored on face grid onto output grid
      procedure, pass :: interpolate_from_vol => grid_interpolate_from_vol    ! Interpolate quantity that is stored on volume grid onto output grid

      ! Manipulation methods
      procedure, pass :: grow => grid_grow          ! Add a new box
      procedure, pass :: shrink => grid_shrink      ! Shrink grid by one box (= merge uppermost 2 boxes)
      procedure, pass :: modify_top_box => grid_modify_top_box  ! Change size of topmost box to reflect niveau change
   END type

!-----------------------------------------------------------
CONTAINS
!-----------------------------------------------------------
   ! Set up grid at program start
   SUBROUTINE grid_init(self, config)
      IMPLICIT NONE
      class(StaggeredGrid), intent(inout) :: self
      class(GridConfig), intent(inout) :: config

      ! Assign config
      self%nz_grid = config%nz_grid
      self%max_length_input_data = config%max_length_input_data
      self%max_depth = config%max_depth
      
      ! USE read config to determine grid size
      CALL self%init_morphology(config)
      CALL self%init_grid_points(config)
   
      ! Allocate arrays according to size
      CALL self%memory_init()

      ! CALL init functions
      CALL self%init_z_axes()
      CALL self%init_areas(config)

   END SUBROUTINE grid_init


   ! Allocate necessary memory based on sizes
   SUBROUTINE grid_memory_init(self)
      IMPLICIT NONE
      class(StaggeredGrid), intent(inout) :: self

      associate (nz_grid=>self%nz_grid)
         ! Allocate memory for initialization using nz_grid
         ! h is already allocated by grid point init
         allocate (self%z_volume(0:nz_grid)) ! Depth axis with center of boxes
         self%z_volume(0) = 0 ! trick for array access - index not in USE

         allocate (self%z_face(nz_grid + 1)) ! Depth axis with faceer border of boxes
         allocate (self%Az(nz_grid + 1)) ! Az is defined on the faces
         allocate (self%dAz(nz_grid)) ! dAz is the difference between Az and thus defined on the volume

         ! Area factors used in calculations
         allocate (self%AreaFactor_1(nz_grid)) ! defined on faces
         allocate (self%AreaFactor_2(nz_grid)) ! defined on faces
         allocate (self%AreaFactor_k1(nz_grid)) ! defined on volumes
         allocate (self%AreaFactor_k2(nz_grid)) ! defined on volumes
         allocate (self%AreaFactor_eps(nz_grid)) ! defined on volumes
         allocate (self%meanint(nz_grid)) ! Inverse ratio of mean height of two adjacent boxes

      END associate
   END SUBROUTINE grid_memory_init


   ! Calculates h and definite number of nz_grid
   ! Depending on configuration
   SUBROUTINE grid_init_grid_points(self, config)
      IMPLICIT NONE
      integer i
      class(StaggeredGrid), intent(inout) :: self
      class(GridConfig), intent(inout) :: config

      self%nz_grid = int(config%nz_grid)

      IF (.not. config%equidistant_grid) THEN ! variable spacing

         ! Grid size given through number of grid_read
         self%nz_grid = config%nz_grid

         ! IF top value not included
         IF (config%grid_read(1) /= (config%max_depth - self%z_zero)) THEN
            CALL error('Top value '//trim(toStr(config%max_depth - self%z_zero))//' is not included in grid file!')
         ENDIF

         !IF maxdepth of grid larger than morphology
         IF (config%grid_read(self%nz_grid + 1) < -self%z_zero) THEN
            CALL error('Grid invalid: maxdepth of grid larger than morphology!')
         ENDIF

         ! IF bottom value not included
         IF (config%grid_read(self%nz_grid + 1) > -self%z_zero) THEN
            CALL error('Bottom value '//trim(toStr(-self%z_zero))//' is not included in grid file!')
         ENDIF

         ! Check for monotonously decreasing grid values
         DO i=2,self%nz_grid
            IF ((config%grid_read(i) - config%grid_read(i-1)) > 0) THEN
            CALL error('Grid input values are not monotonously decreasing')
            ENDIF
         ENDDO
      ENDIF

      ! Construct h
      allocate (self%h(0:self%nz_grid + 1))
      self%h(0) = 0 ! Note that h(0) has no physical meaning but helps with some calculations
      self%h(self%nz_grid + 1) = 0 ! Note that h(nz_grid + 1) has no physical meaning but helps with some calculations
      IF (config%equidistant_grid) THEN
         ! Equidistant grid
         self%h(1:self%nz_grid) = config%max_depth/self%nz_grid
      ELSE
         ! Set up h according to configuration
         DO i = 2, self%nz_grid + 1
            self%h(2 + self%nz_grid - i) = config%grid_read(i - 1) - config%grid_read(i)
         ENDDO
      ENDIF
   END SUBROUTINE grid_init_grid_points


   ! Initializes z_volume and z_face
   SUBROUTINE grid_init_z_axes(self)
      IMPLICIT NONE
      class(StaggeredGrid), intent(inout) :: self
      integer :: i

      ! Compute position of layer center and top
      self%z_volume(1) = 1.0_RK
      self%z_face(1) = 0.0_RK
      DO i = 1, self%nz_grid
         self%z_volume(i) = self%z_volume(i - 1) + 0.5_RK*(self%h(i - 1) + self%h(i))
         self%z_volume(i-1) = nint(1e6_RK*self%z_volume(i-1))/1e6_RK
      ENDDO
      self%z_volume(self%nz_grid) = nint(1e6_RK*self%z_volume(self%nz_grid))/1e6_RK

      ! DO i = 2, self%nz_grid + 1
      !     self%z_face(i) = self%z_face(i - 1) + self%h(i - 1)
      !     self%z_face(i - 1) = nint(1e6_RK*self%z_face(i - 1))/1e6_RK
      ! ENDDO
      ! self%z_face(self%nz_grid + 1) = nint(1e6_RK*self%z_face(self%nz_grid + 1))/1e6_RK
      
      DO i = 2, self%nz_grid + 1
         self%z_face(i) = self%z_face(i - 1) + self%h(i - 1)
      ENDDO
   END SUBROUTINE grid_init_z_axes


   ! Init morphology variables
   SUBROUTINE grid_init_morphology(self, config)
      IMPLICIT NONE
      class(StaggeredGrid), intent(inout) :: self
      class(GridConfig), intent(inout) :: config
      integer :: num_read

      num_read = size(config%A_read)

      self%z_zero = config%z_A_read(1) ! z_zero is the negative value of the lowest depth in the morphology file

      config%z_A_read = self%z_zero - config%z_A_read ! z-coordinate is positive upwards, zero point is at reservoir bottom

   END SUBROUTINE grid_init_morphology


   ! Init areas...
   SUBROUTINE grid_init_areas(self, config)
      IMPLICIT NONE
      class(StaggeredGrid), intent(inout) :: self
      class(GridConfig), intent(inout) :: config

      integer :: num_read
      associate (nz_grid=>self%nz_grid, dAz=>self%dAz, z_face=>self%z_face, Az=>self%Az)
         num_read = size(config%A_read)

         ! Interpolate area (A) at all depths (z_face)
         CALL Interp(config%z_A_read, config%A_read, num_read, self%z_face, Az, nz_grid + 1)

         ! Compute area derivative (= projected sediment area over layer thickness)
         dAz(1:nz_grid) = (Az(2:nz_grid + 1) - Az(1:nz_grid))/(z_face(2:nz_grid + 1) - z_face(1:nz_grid))
      END associate
   END SUBROUTINE


   ! Updates area factors (good method names are self-explanatory)
   SUBROUTINE grid_update_area_factors(self)
      IMPLICIT NONE
      class(StaggeredGrid), intent(inout) :: self

      integer :: i
      associate (Az=>self%Az, &
                  h=>self%h, &
                  nz=>self%nz_occupied)

         self%AreaFactor_1(1:nz) = -4*Az(1:nz)/(h(1:nz) + h(0:nz - 1))/h(1:nz)/(Az(2:nz + 1) + Az(1:nz))
         self%AreaFactor_2(1:nz) = -4*Az(2:nz + 1)/(h(1:nz) + h(2:nz + 1))/h(1:nz)/(Az(2:nz + 1) + Az(1:nz))
         self%AreaFactor_k1(1:nz - 1) = -(Az(2:nz) + Az(1:nz - 1))/(h(1:nz - 1) + h(2:nz))/h(1:nz - 1)/Az(2:nz)
         self%AreaFactor_k1(nz) = 0
         self%AreaFactor_k2(1:nz - 1) = -(Az(2:nz) + Az(3:nz + 1))/(h(1:nz - 1) + h(2:nz))/h(2:nz)/Az(2:nz)
         self%AreaFactor_k2(nz) = 0
         self%AreaFactor_eps(1:nz - 1) = 0.5_RK*((Az(2:nz) - Az(1:nz - 1))/h(1:nz - 1) + (Az(3:nz + 1) - Az(2:nz))/h(2:nz))/Az(2:nz)
         self%AreaFactor_eps(nz) = 0
         self%meanint(1:nz) = 2.0_RK/(h(1:nz) + h(2:nz + 1))

         self%volume = 0
         DO i = 1, nz
            self%volume = self%volume + 0.5_RK*h(i)*(Az(i) + Az(i + 1))
         ENDDO
      END associate
   END SUBROUTINE grid_update_area_factors


   ! Modify size of topmost box
   SUBROUTINE grid_modify_top_box(self, dh)
      IMPLICIT NONE
      class(StaggeredGrid), intent(inout) :: self
      real(RK) :: dh

      associate (z_face=>self%z_face, &
                  ubnd_fce=>self%ubnd_fce, &
                  z_volume=>self%z_volume, &
                  ubnd_vol=>self%ubnd_vol, &
                  h=>self%h, &
                  Az=>self%Az, &
                  dAz=>self%dAz)

         IF (h(ubnd_vol) > self%h_old) THEN
            Az(ubnd_fce) = Az(ubnd_fce) + dAz(ubnd_vol + 1)*dh
         ELSE
            Az(ubnd_fce) = Az(ubnd_fce) + dAz(ubnd_vol)*dh
         ENDIF

         h(ubnd_vol) = h(ubnd_vol) + dh
         z_volume(ubnd_vol) = z_volume(ubnd_vol) + 0.5_RK*dh
         z_face(ubnd_fce) = z_face(ubnd_fce) + dh

      END associate
   END SUBROUTINE


   !Shrink grid (mainly used by advection)
   SUBROUTINE grid_shrink(self, dh)
      IMPLICIT NONE
      class(StaggeredGrid), intent(inout) :: self
      real(RK) :: dh

      associate (z_face=>self%z_face, &
                  ubnd_fce=>self%ubnd_fce, &
                  z_volume=>self%z_volume, &
                  ubnd_vol=>self%ubnd_vol, &
                  h=>self%h, &
                  nz_occupied=>self%nz_occupied, &
                  Az=>self%Az)

         ! Update grid
         z_face(ubnd_fce - 1) = z_face(ubnd_fce) + dh
         z_volume(ubnd_vol - 1) = 0.5_RK*(z_face(ubnd_fce - 1) + z_face(ubnd_fce - 2))

         Az(ubnd_fce - 1) = Az(ubnd_fce)

         self%h_old = h(ubnd_vol - 1)
         h(ubnd_vol - 1) = (z_face(ubnd_fce - 1) - z_face(ubnd_fce - 2))

         ! Update number of occupied cells
         nz_occupied = nz_occupied - 1
         ! Update boundaries (ubnd_fce and ubnd_vol)
         CALL self%update_nz()

      END associate
   END SUBROUTINE


   ! Grow grid (mainly used by advection)
   SUBROUTINE grid_grow(self, dh)
      IMPLICIT NONE
      class(StaggeredGrid), intent(inout) :: self
      real(RK) :: dh

      associate (z_face=>self%z_face, &
                  ubnd_fce=>self%ubnd_fce, &
                  z_volume=>self%z_volume, &
                  ubnd_vol=>self%ubnd_vol, &
                  h=>self%h, &
                  Az=>self%Az, &
                  dAz=>self%dAz, &
                  nz_occupied=>self%nz_occupied)

         h(ubnd_vol + 1) = (h(ubnd_vol) + dh)/2
         h(ubnd_vol) = (h(ubnd_vol) + dh)/2
         self%h_old = h(ubnd_vol + 1)

         z_face(ubnd_fce + 1) = z_face(ubnd_fce) + dh
         z_face(ubnd_fce) = z_face(ubnd_fce) + dh - h(ubnd_vol)

         z_volume(ubnd_vol + 1) = z_volume(ubnd_vol) + (h(ubnd_vol + 1) + dh)/2
         z_volume(ubnd_vol) = z_volume(ubnd_vol) - (h(ubnd_vol) - dh)/2

         Az(ubnd_fce) = Az(ubnd_fce - 1) + h(ubnd_vol)*dAz(ubnd_vol)
         Az(ubnd_fce + 1) = Az(ubnd_fce) + h(ubnd_vol + 1)*dAz(ubnd_vol + 1)

         ! Update number of occupied cells
         nz_occupied = nz_occupied + 1

         CALL self%update_nz()

      END associate
   END SUBROUTINE


   ! Recalculate depth
   SUBROUTINE grid_update_depth(self, new_depth)
      IMPLICIT NONE
      class(StaggeredGrid), intent(inout) :: self
      real(RK), intent(inout) :: new_depth
      integer :: i
      real(RK) ::  zmax
      real(RK), parameter :: epsilon = 1.0e-4_RK
      character(len=32) :: str1, str2
      associate (nz_grid=>self%nz_grid, &
                  nz_occupied=>self%nz_occupied, &
                  z_face=>self%z_face, &
                  z_volume=>self%z_volume, &
                  h=>self%h, &
                  z_zero=>self%z_zero, &
                  Az=>self%Az)

         DO i = 1, nz_grid + 1

               !+WMEJ: The original judgment required too high model accuracy, 
               !       but it has now been modified to only consider four decimal places,
               !       which is completely acceptable for lake processes.
               write(str1, '(F15.4)') z_face(i)
               write(str2, '(F15.4)') z_zero - new_depth
               ! IF (z_face(i) >= (z_zero - new_depth)) THEN ! IF above initial water level
               IF (str1 >= str2) THEN ! IF above initial water level
                  zmax = z_face(i)
                  ! Set top face to new water level
                  z_face(i) = z_zero - new_depth

                  ! Adjust new volume center of cell i-1 (belonging to upper face i)
                  z_volume(i - 1) = (z_face(i) + z_face(i - 1))/2
                  h(i - 1) = z_face(i) - z_face(i - 1)

                  Az(i) = Az(i-1) + h(i - 1)/(zmax-z_face(i-1))*(Az(i)-Az(i-1))
                  nz_occupied = i - 1
                  self%h_old = h(nz_occupied)
                  ! IF (h(nz_occupied) <= 0.5*h(nz_occupied - 1)) THEN ! IF top box is too small
                  !     z_face(nz_occupied) = z_face(nz_occupied + 1) ! Combine the two upper boxes
                  !     z_volume(nz_occupied - 1) = (z_face(nz_occupied) + z_face(nz_occupied - 1))/2
                  !     h(nz_occupied - 1) = h(nz_occupied) + h(nz_occupied - 1)
                  !     self%h_old = h(nz_occupied - 1)
                  !     nz_occupied = nz_occupied - 1 ! Reduce number of boxes
                  ! ENDIF
                  exit
               ENDIF
         ENDDO

         CALL self%update_nz()

      END associate
   END SUBROUTINE


   ! Interpolates values of y on grid z onto array yi and grid z_volume
   SUBROUTINE grid_interpolate_to_vol(self, z, y, num_z, yi)
      IMPLICIT NONE
      class(StaggeredGrid), intent(in) :: self
      real(RK), DIMENSION(:), intent(in) :: z, y
      real(RK), DIMENSION(:), intent(out) :: yi

      integer, intent(in) :: num_z
      CALL Interp(z, y, num_z, self%z_volume(1:self%nz_grid), yi, self%nz_grid)
   END SUBROUTINE


   SUBROUTINE grid_interpolate_to_face(self, z, y, num_z, yi)
      IMPLICIT NONE
      class(StaggeredGrid), intent(in) :: self
      real(RK), DIMENSION(:), intent(in) :: z, y
      real(RK), DIMENSION(:), intent(out) :: yi

      integer, intent(in) :: num_z
      CALL Interp(z, y, num_z, self%z_face, yi, self%nz_grid)
   END SUBROUTINE


   SUBROUTINE grid_interpolate_to_face_from_second(self, z, y, num_z, yi)
      IMPLICIT NONE
      class(StaggeredGrid), intent(in) :: self
      real(RK), DIMENSION(:), intent(in) :: z, y
      real(RK), DIMENSION(:), intent(out) :: yi

      integer, intent(in) :: num_z
      CALL Interp(z, y, num_z, self%z_face(2:self%nz_grid + 1), yi(2:self%nz_grid + 1), self%nz_grid)
   END SUBROUTINE


   SUBROUTINE grid_interpolate_from_vol(self, y, zi, yi, num_zi, output_depth_reference)
      class(StaggeredGrid), intent(in) :: self
      real(RK), DIMENSION(:), intent(in) :: zi, y
      real(RK), DIMENSION(:), intent(out) :: yi
      integer, intent(in) :: num_zi
      character(len=:), allocatable :: output_depth_reference

      real(RK), DIMENSION(self%ubnd_vol) :: z_volume_mod
      integer :: i

      ! Transform z_volume for interpolation on zout grid depending on grid reference chosen in par-file
      IF (output_depth_reference == 'bottom') THEN
         z_volume_mod(1) = self%z_face(1)
         z_volume_mod(2:self%ubnd_vol - 1) = self%z_volume(2:self%ubnd_vol - 1)
         z_volume_mod(self%ubnd_vol) = self%z_face(self%ubnd_fce)
      ELSE IF (output_depth_reference == 'surface') THEN
         z_volume_mod(1) = self%z_face(1) - self%z_face(self%ubnd_fce)
         DO i = 2, self%ubnd_vol-1
            z_volume_mod(i) = self%z_volume(i) - self%z_face(self%ubnd_fce)
         ENDDO
         z_volume_mod(self%ubnd_vol) = 0
      ENDIF

      CALL Interp_nan(z_volume_mod(1:self%ubnd_vol), y(1:self%ubnd_vol), self%ubnd_vol, zi, yi, num_zi)
   END SUBROUTINE


   SUBROUTINE grid_interpolate_from_face(self, y, zi, yi, num_zi, output_depth_reference)
      class(StaggeredGrid), intent(in) :: self
      real(RK), DIMENSION(:), intent(in) :: zi, y
      real(RK), DIMENSION(:), intent(out) :: yi
      integer, intent(in) :: num_zi
      character(len=:), allocatable :: output_depth_reference

      real(RK), DIMENSION(self%ubnd_fce) :: z_face_mod
      integer :: i

      ! Transform z_face for interpolation on zout grid depending on grid reference chosen in par-file
      IF (output_depth_reference == 'bottom') THEN
         z_face_mod(1:self%ubnd_fce) = self%z_face(1:self%ubnd_fce)
      ELSE IF (output_depth_reference == 'surface') THEN
         DO i = 1, self%ubnd_fce
            z_face_mod(i) = self%z_face(i) - self%z_face(self%ubnd_fce)
         ENDDO
      ENDIF

      CALL Interp_nan(z_face_mod(1:self%ubnd_fce), y(1:self%ubnd_fce), self%ubnd_fce, zi, yi, num_zi)
   END SUBROUTINE


   ! Update all upper bounds and lengths
   SUBROUTINE grid_update_nz(self)
      IMPLICIT NONE
      class(StaggeredGrid), intent(inout) :: self

      self%ubnd_vol = self%nz_occupied
      self%ubnd_fce = self%nz_occupied + 1
      self%length_vol = self%nz_occupied
      self%length_fce = self%nz_occupied + 1
   END SUBROUTINE


   pure function grid_convert2height_above_sed(z, z_zero) result(h)
      IMPLICIT NONE
      real(RK), DIMENSION(:), intent(in) :: z
      real(RK), intent(in) :: z_zero

      real(RK), DIMENSION(size(z)) :: h
      integer :: n

      n = size(z)
      h = -z_zero + z(n:1:-1)
   END function

END MODULE  strat_grid





!     +---------------------------------------------------------------+
!     |  Data structure definitions for simulation data
!     +---------------------------------------------------------------+
MODULE  strat_simdata
   USE strat_kinds, only: RK
   USE strat_grid
   USE MOD_Lake_Const 
   USE strat_utilities

   IMPLICIT NONE

!-----------------------------------------------------------
PRIVATE
!-----------------------------------------------------------

   ! Simulation configuration
   type, public :: SimConfig
      integer :: timestep
      integer :: start_year
      real(RK) :: start_datum
      real(RK) :: end_datum
      integer :: disp_simulation
      logical :: continue_from_snapshot = .false.
   END type

   ! Model configuration (read from file)
   type, public :: ModelConfig
      integer :: max_length_input_data
      logical :: couple_aed2
      logical :: has_advection
      integer :: turbulence_model
      logical :: split_a_seiche
      integer :: stability_func
      integer :: flux_condition
      integer :: forcing_mode
      logical :: use_filtered_wind
      integer :: seiche_normalization
      integer :: wind_drag_model
      integer :: inflow_placement
      integer :: pressure_gradients
      logical :: salinity_transport
      integer :: ice_model
      integer :: snow_model
   END type

   ! Model params (read from file)
   type, public :: ModelParam
      real(RK) :: Lat
      real(RK) :: p_air
      real(RK) :: a_seiche
      real(RK) :: a_seiche_w
      real(RK) :: strat_sumr
      real(RK) :: q_NN
      real(RK) :: f_wind
      real(RK) :: C10_constant
      real(RK) :: CD
      real(RK) :: fgeo
      real(RK) :: p_sw
      real(RK) :: p_lw
      real(RK) :: p_windf
      real(RK) :: beta_sol
      real(RK) :: wat_albedo
      real(RK) :: p_albedo
      real(RK) :: freez_temp
      real(RK) :: snow_temp
   END type

   ! Model state (this is actually the simulation data!!!)
   type, public :: ModelState
      ! Iteration variables
      integer :: current_year ! Current year of simulation, used for zenith angle dependent water albedo
      integer :: current_month ! Current month of simulation, used for zenith angle dependent water albedo
      real(RK) :: current_day ! Current day of simulation, used for zenith angle dependent water albedo
      real(RK) :: datum, dt
      integer(8), DIMENSION(2) :: simulation_time, simulation_time_old
      logical :: first_timestep = .true.

      ! Variables located on z_cent grid
      ! Note that for these variables the value at 0 z.b. U(0) is not used
      real(RK), DIMENSION(:), allocatable :: U, V ! Water velocities [m/s]
      real(RK), DIMENSION(:), allocatable :: T, S ! Temperature [°C], Salinity [‰]
      real(RK), DIMENSION(:), allocatable :: dS ! Source/sink for salinity
      real(RK), DIMENSION(:, :), allocatable :: Q_inp ! Horizontal inflow [m^3/s]
      real(RK), DIMENSION(:), allocatable :: rho ! Water density [kg/m^3]

      ! Variables located on z_upp grid
      real(RK), DIMENSION(:), allocatable :: k, ko ! Turbulent kinetic energy (TKE) [J/kg]
      real(RK), DIMENSION(:), allocatable :: avh
      real(RK), DIMENSION(:), allocatable :: eps ! TKE dissipation rate [W/kg]
      real(RK), DIMENSION(:), allocatable :: num, nuh ! Turbulent viscosity (momentum) and diffusivity (temperature)
      real(RK), DIMENSION(:), allocatable :: P, B ! Shear stress production [W/kg], buoyancy production [W/kg]
      real(RK), DIMENSION(:), allocatable :: NN ! Brunt-Väisälä frequency [s-2]
      real(RK), DIMENSION(:), allocatable :: cmue1, cmue2 ! Model constants
      real(RK), DIMENSION(:), allocatable :: P_Seiche ! Production of TKE [W/kg] and seiche energy [J]
      real(RK) :: E_Seiche
      real(RK) :: gamma ! Proportionality constant for loss of seiche energy

      real(RK), DIMENSION(:), allocatable :: absorb ! Absorption coeff [m-1]
      real(RK) :: u10, v10, uv10, Wf ! Wind speeds, wind factor
      real(RK) :: u_taub, drag, u_taus ! Drag
      real(RK) :: tx, ty ! Shear stress
      real(RK) :: C10 ! Wind drag coefficient
      real(RK) :: SST, heat, heat_snow, heat_ice, heat_snowice! Sea surface temperature and heat flux
      real(RK) :: T_atm ! Air temp at surface
      real(RK), DIMENSION(:), allocatable :: rad ! Solar radiation (in water)
      real(RK), DIMENSION(:), allocatable :: Q_vert ! Vertical exchange between boxes
      real(RK), DIMENSION(9,12) :: albedo_data  ! Experimental monthly albedo data for determination of current water albedo
      real(RK) :: albedo_water   ! Current water albedo
      integer :: lat_number ! Latitude band (used for determination of albedo)

      ! Snow and Ice
      real(RK), allocatable :: snow_h ! Snow layer height [m]
      real(RK), allocatable :: total_ice_h ! Total ice layer height [m]
      real(RK), allocatable :: black_ice_h ! Black ice layer height [m]
      real(RK), allocatable :: white_ice_h ! Snowice layer height [m]
      real(RK), allocatable :: snwml ! Snow layer melt [mm/s]
      real(RK), allocatable :: iceml ! Total ice layer melt [mm/s]
      real(RK) :: snow_dens ! Snow density [kg m-3]
      real(RK) :: ice_temp ! Ice temperature [°C]
      real(RK) :: T_surf ! Surface temperature [°C]
      real(RK) :: precip ! Precipiation in water eqvivalent hight [m]

      !For saving heatflux
      real(RK), allocatable :: ha ! Incoming long wave [W m-2]
      real(RK), allocatable :: hw ! Outgoing long wave [W m-2]
      real(RK), allocatable :: hk ! Sensible flux [W m-2]
      real(RK), allocatable :: hv ! Latent heat [W m-2]
      real(RK), allocatable :: rad0 !  Solar radiation at surface  [W m-2]



      real(RK) :: cde, cm0
      real(RK) ::  fsed
      real(RK), DIMENSION(:), allocatable     :: fgeo_add
      logical, DIMENSION(1:4) :: has_surface_input, has_deep_input
      integer :: nz_input
   CONTAINS
      procedure, pass :: init => model_state_init
      procedure, pass :: update => update_model_state
   END type

   ! Structure that encapsulates a full program state
   type, public :: SimulationData
      ! type(InputConfig), public   :: input_cfg
      ! type(OutputConfig), public  :: output_cfg
      type(SimConfig), public     :: sim_cfg
      type(ModelConfig), public   :: model_cfg
      type(ModelParam), public    :: model_param
      type(ModelState), public    :: model
      type(StaggeredGrid), public :: grid
   CONTAINS
      procedure, pass :: init => simulation_data_init
   END type

!-----------------------------------------------------------
CONTAINS
!-----------------------------------------------------------

   SUBROUTINE simulation_data_init(self, state_size)
      class(SimulationData), intent(inout) :: self
      integer, intent(in) :: state_size
      ! Init model data structures
      CALL self%model%init(state_size)
   END SUBROUTINE


   ! Allocates all arrays of the model state in the correct size
   SUBROUTINE model_state_init(self, state_size)
      class(ModelState), intent(inout) :: self
      integer, intent(in) :: state_size

      ! Values on volume grid
      ! Important: Size is smaller than vars on upper grid.
      !            https://en.wikipedia.org/wiki/Off-by-one_error#Fencepost_error ;-)
      allocate (self%U(state_size))
      allocate (self%V(state_size))
      allocate (self%T(state_size))
      allocate (self%S(state_size))
      allocate (self%dS(state_size))
      allocate (self%Q_inp(1:4, state_size + 1))
      allocate (self%rho(state_size))
      allocate (self%avh(state_size))

      ! Values on z_upp grid
      allocate (self%K(state_size + 1))
      allocate (self%ko(state_size + 1))
      allocate (self%eps(state_size + 1))
      allocate (self%num(state_size + 1))
      allocate (self%nuh(state_size + 1))
      allocate (self%P(state_size + 1))
      allocate (self%B(state_size + 1))
      allocate (self%NN(state_size + 1))
      allocate (self%cmue1(state_size + 1))
      allocate (self%cmue2(state_size + 1))
      allocate (self%P_Seiche(state_size + 1))

      allocate (self%absorb(state_size + 1))
      allocate (self%rad(state_size + 1))
      allocate (self%Q_vert(state_size + 1))

      allocate (self%snow_h)
      allocate (self%total_ice_h)
      allocate (self%black_ice_h)
      allocate (self%white_ice_h)
      allocate (self%snwml)    
      allocate (self%iceml)

      allocate (self%ha)
      allocate (self%hw)
      allocate (self%hk)
      allocate (self%hv)
      allocate (self%rad0)

      ! Init to zero
      self%U = 0.0_RK
      self%V = 0.0_RK
      self%T = 0.0_RK
      self%S = 0.0_RK
      self%dS = 0.0_RK
      self%Q_inp = 0.0_RK
      self%rho = 0.0_RK

      self%K = 0.0_RK
      self%ko = 0.0_RK
      self%eps = 0.0_RK
      self%num = 0.0_RK
      self%nuh = 0.0_RK
      self%P = 0.0_RK
      self%B = 0.0_RK
      self%NN = 0.0_RK
      self%cmue1 = 0.0_RK
      self%cmue2 = 0.0_RK
      self%P_Seiche = 0.0_RK
      self%E_Seiche = 0.0_RK

      self%absorb = 0.0_RK
      self%rad = 0.0_RK
      self%Q_vert = 0.0_RK

      self%snow_h = 0.0_RK
      self%total_ice_h = 0.0_RK
      self%black_ice_h = 0.0_RK
      self%white_ice_h = 0.0_RK
      self%ice_temp = 0.0_RK
      self%snow_dens = rho_s_0
      self%precip = 0.0_RK
      self%snwml  = 0.0_RK 
      self%iceml  = 0.0_RK

      self%ha = 0.0_RK
      self%hw = 0.0_RK
      self%hk = 0.0_RK
      self%hv = 0.0_RK
      self%rad0 = 0.0_RK

      self%simulation_time_old = 0

   END SUBROUTINE


   ! load model state unformatted
   SUBROUTINE update_model_state ( self      ,&
               state_size, snwdp     , bsnwdp   , wsnwdp   ,&
               rhosnw    , uwatv     , vwatv    , lktmp    ,&
               lksal     , tke       , eps      , etke     ,&
               lkrho    , num       , nuh      , tmice    ,& 
               Qinp      )


! ---------------------------------- code history -------------------------------------
! Description:
!     Update the model state variables at the previous moment,
!
! Original author: 
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024
!
! -------------------------------------------------------------------------------------
   USE MOD_Lake_Const, only: tfrz
!================================================================================
!  ------------------------- inout variables ---------------------------
   class(ModelState), intent(inout) :: self

!  ------------------------- input variables ---------------------------
   integer, intent(in)      :: &
      state_size                ! number of lake layers

   real(RK), intent(in)     :: &
      etke                    ,&! Seiche energy [J]
      snwdp                   ,&! Snow layer height [m]
      bsnwdp                  ,&! Black ice layer height [m]
      wsnwdp                  ,&! Snowice layer height [m]
      tmice                   ,&! Ice temperature [°C]
      rhosnw                    ! Snow density [kg m-3]

   real(RK), intent(in)     :: &
      uwatv(state_size)       ,&! Water velocity in x-direction [m/s]
      vwatv(state_size)       ,&! Water velocity in y-direction [m/s]
      lktmp(state_size)       ,&! Temperature [°C]
      lksal(state_size)       ,&! Salinity [‰]
      lkrho(state_size)       ,&! Water density [kg/m^3]
      tke(state_size+1)       ,&! Turbulent kinetic energy (TKE) [J/kg]
      eps(state_size+1)       ,&! TKE dissipation rate [W/kg]
      num(state_size+1)       ,&! Turbulent viscosity (momentum)
      nuh(state_size+1)         ! Turbulent diffusivity (heat)

   real(RK), intent(in)     :: &
      Qinp(1:4, state_size+1)   ! Horizontal inflow [m^3/s]
!================================================================================
      self%U = uwatv
      self%V = vwatv
      self%T = lktmp
      self%S = lksal 
      self%Q_inp = Qinp
      self%rho = lkrho

      self%k = tke
      self%eps = eps
      self%num = num
      self%nuh = nuh
      self%E_Seiche = etke

      self%snow_h = snwdp
      self%total_ice_h = snwdp
      self%black_ice_h = bsnwdp
      self%white_ice_h = wsnwdp
      self%ice_temp = tmice
      self%snow_dens = rhosnw

      !-WMEJ maybe need to be updated, 
      !-     The status of the previous moment will not affect the next moment.
      !-     Not sure IF the following variables are needed.
      ! self%precip =  precip
      ! self%avh = avh
      ! self%ko = tkeo
      ! self%P = shearp
      ! self%B = buoy
      ! self%NN = NN 
      ! self%cmue1 = cmue1
      ! self%cmue2 = cmue2
      ! self%P_Seiche = ptke
      ! self%absorb = absorb
      ! self%Q_vert = Qvert
      ! self%rad = radin  
      ! self%dS = dsali
      ! self%ha = lwdns
      ! self%hw = lwups
      ! self%hk = fsena
      ! self%hv = lfevpa
      ! self%rad0 = sabg

   END SUBROUTINE update_model_state

END MODULE  strat_simdata





!<    +---------------------------------------------------------------+
!     |  Solver MODULE
!     | CONTAINS implementation of a tridiagonal matrix solver
!<    +---------------------------------------------------------------+
MODULE  strat_solver

   USE strat_kinds, only: RK, LinSysSolver

   IMPLICIT NONE

!-----------------------------------------------------------
PRIVATE
!-----------------------------------------------------------

   ! Type that represents a solver based on the Thomas algorithm
   ! https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
   type, extends(LinSysSolver), public :: ThomasAlgSolver
   CONTAINS
      procedure, nopass :: solve => solve_tridiag_thomas
   END type

!-----------------------------------------------------------
CONTAINS
!-----------------------------------------------------------

   ! Implementation of Thomas algorithm
   SUBROUTINE solve_tridiag_thomas(ld, md, ud, rhs, x, ubnd)
      IMPLICIT NONE

      ! Arguments
      real(RK), intent(in) :: ld(:), md(:), ud(:) ! Diagonals (A)
      real(RK), intent(in) :: rhs(:) ! right-hand side (b)
      real(RK), intent(out) :: x(:) ! solution (x)
      integer, intent(in) :: ubnd ! upper bound of x

      ! Local variables
      real(RK), DIMENSION(size(md)) :: ru, qu
      integer :: i

      ru(ubnd) = ud(ubnd)/md(ubnd)
      qu(ubnd) = rhs(ubnd)/md(ubnd)

      DO i = ubnd - 1, 2, -1
         ru(i) = ud(i)/(md(i) - ld(i)*ru(i + 1))
         qu(i) = (rhs(i) - ld(i)*qu(i + 1))/(md(i) - ld(i)*ru(i + 1))
      ENDDO

      qu(1) = (rhs(1) - ld(1)*qu(2))/(md(1) - ld(1)*ru(2))

      x(1) = qu(1)
      DO i = 2, ubnd
         x(i) = qu(i) - ru(i)*x(i - 1)
      ENDDO

      RETURN
   END SUBROUTINE

END MODULE  strat_solver





!<    +---------------------------------------------------------------+
!     |  Interface and implementation for discretization methods
!     |  - Currently, there are two impelmentations:
!     |  - EulerIDiscretizationMFQ: Euler implicit discretization for Mean Flow Quantities
!     |  - EulerIDiscretizationKEPS: Euler Implicit discretization for k and epsilon
!<    +---------------------------------------------------------------+
MODULE  strat_discretization
   USE strat_kinds, only: RK
   USE strat_simdata
   USE strat_grid

   IMPLICIT NONE

!-----------------------------------------------------------
PRIVATE
!-----------------------------------------------------------

   ! Interface
   type, abstract, public :: Discretization
      class(StaggeredGrid), pointer :: grid
   CONTAINS
      procedure, pass(self), public :: init => generic_init
      procedure(generic_create_LES), deferred, pass(self), public ::create_LES
   END type

   ! Implementation for Mean flow Quantities
   type, extends(Discretization), public :: EulerIDiscretizationMFQ
   CONTAINS
      procedure, pass, public :: create_LES => euleri_create_LES_MFQ
   END type

   ! Implementation for K and eps
   type, extends(Discretization), public :: EulerIDiscretizationKEPS
   CONTAINS
      procedure, pass, public :: create_LES => euleri_create_LES_KEPS
   END type

!-----------------------------------------------------------
CONTAINS
!-----------------------------------------------------------

   SUBROUTINE generic_init(self, grid)
      class(Discretization), intent(inout) :: self
      type(StaggeredGrid), target :: grid
      self%grid => grid
   END SUBROUTINE


   ! Abstract method definition for creating a linear equation system (LES)
   SUBROUTINE generic_create_LES(self, var, nu, sources, boundaries, lower_diag, main_diag, upper_diag, rhs, dt)
      class(Discretization), intent(inout) :: self
      real(RK), DIMENSION(:), intent(inout) :: var, sources, boundaries, lower_diag, upper_diag, main_diag, rhs, nu
      real(RK), intent(inout) :: dt
   END SUBROUTINE


   ! Implementation for Mean Flow Quantities
   SUBROUTINE euleri_create_LES_MFQ(self, var, nu, sources, boundaries, lower_diag, main_diag, upper_diag, rhs, dt)
      class(EulerIDiscretizationMFQ), intent(inout) :: self
      real(RK), DIMENSION(:), intent(inout) :: var, sources, boundaries, lower_diag, upper_diag, main_diag, rhs, nu
      real(RK), intent(inout) :: dt
      integer :: n

      n=self%grid%ubnd_vol

      ! Build diagonals
      upper_diag(1) = 0.0_RK
      upper_diag(2:n) = dt*nu(2:n)*self%grid%AreaFactor_1(2:n)
      lower_diag(1:n - 1) = dt*nu(2:n)*self%grid%AreaFactor_2(1:n-1)
      lower_diag(n) = 0.0_RK
      main_diag(1:n) = 1.0_RK - upper_diag(1:n) - lower_diag(1:n) + boundaries(1:n)*dt

      ! Calculate RHS
      ! A*phi^{n+1} = phi^{n}+dt*S^{n}
      rhs(1:n) = var(1:n) + dt*sources(1:n)
   END SUBROUTINE


   ! Implementation for K and EPS (note different usage of nu)
   SUBROUTINE euleri_create_LES_KEPS(self, var, nu, sources, boundaries, lower_diag, main_diag, upper_diag, rhs, dt)
      class(EulerIDiscretizationKEPS), intent(inout) :: self
      real(RK), DIMENSION(:), intent(inout) :: var, sources, boundaries, lower_diag, upper_diag, main_diag, rhs, nu
      real(RK), intent(inout) :: dt
      integer :: n

      n=self%grid%ubnd_fce

      ! Build diagonals
      upper_diag(1) = 0
      upper_diag(n) = 0
      lower_diag(1) = 0
      lower_diag(n) = 0
      upper_diag(2:n - 1) = dt*nu(1:n - 2)*self%grid%AreaFactor_k1(1:n - 2)
      lower_diag(2:n - 1) = dt*nu(2:n - 1)*self%grid%AreaFactor_k2(1:n - 2)
      main_diag(1:n) = 1.0_RK - upper_diag(1:n) - lower_diag(1:n) + boundaries(1:n)*dt

      ! Calculate RHS
      ! A*phi^{n+1} = phi^{n}+dt*S^{n}
      rhs(1:n) = var(1:n) + dt*sources(1:n)

   END SUBROUTINE

END MODULE  strat_discretization





!<    +---------------------------------------------------------------+
!     | Abstract base class for simulation variables (K, EPS, T, S etc)
!     | The calc_terms and pos_solve methods should be overwritten according
!<    +---------------------------------------------------------------+
MODULE  strat_statevar
   USE strat_kinds, only: RK, LinSysSolver
   USE strat_simdata
   USE strat_grid
   USE strat_discretization
   USE strat_solver

   IMPLICIT NONE

!-----------------------------------------------------------
PRIVATE
!-----------------------------------------------------------

   type, abstract, public :: ModelVariable
      ! Things a model variable needs:
      class(LinSysSolver), pointer :: solver      ! Solver
      class(Discretization), pointer :: disc      ! Discretization scheme
      class(StaggeredGrid), pointer :: grid       ! Grid it lives on
      class(ModelConfig), pointer :: cfg          ! Configuration of the model
      real(RK), DIMENSION(:), pointer :: nu       ! nu concerning this variable
      real(RK), DIMENSION(:), pointer :: var      ! pointer to variable data (usually points to some state var in modelstate)
      integer, pointer :: ubnd                    ! Current upper bound of variable
   CONTAINS
      procedure, pass(self), public :: init => generic_var_init
      procedure(generic_var_calc_terms), deferred, pass(self), public :: calc_terms
      procedure, pass(self), public :: update => generic_var_update
      procedure, pass(self), public :: post_solve => generic_var_post_solve
   END type

!-----------------------------------------------------------
CONTAINS
!-----------------------------------------------------------

   SUBROUTINE generic_var_init(self, cfg, grid, solver, disc, nu, var, ubnd)
      class(ModelVariable), intent(inout) :: self
      class(LinSysSolver), target :: solver
      class(StaggeredGrid), target :: grid
      class(Discretization), target :: disc
      class(ModelConfig), target :: cfg

      real(RK), DIMENSION(:), target :: nu, var
      integer, target :: ubnd

      ! Assign pointers
      self%cfg => cfg
      self%grid => grid
      self%solver => solver
      self%disc => disc
      self%nu => nu
      self%var => var
      self%ubnd => ubnd
   END SUBROUTINE


   ! Generic base method for update of a variable
   SUBROUTINE generic_var_update(self, state, param)
      class(ModelVariable), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(ModelParam), intent(inout) :: param
      real(RK), DIMENSION(size(self%var)) :: main_diag, rhs, sources, lower_diag, upper_diag, boundaries

      ! Calculate source terms
      CALL self%calc_terms(state, param, sources, boundaries)

      ! Create linear system of equations
      CALL self%disc%create_LES(self%var, self%nu, sources, boundaries, lower_diag, main_diag, upper_diag, rhs, state%dt)

      ! Solve LES
      CALL self%solver%solve(lower_diag, main_diag, upper_diag, rhs, self%var, self%ubnd)

      ! DO post processing (e.g. set boundary values)
      CALL self%post_solve(state)

   END SUBROUTINE


   ! Generic interface for method that calculates source terms
   SUBROUTINE generic_var_calc_terms(self, state, param, sources, boundaries)
      class(ModelVariable), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(ModelParam), intent(inout) :: param
      real(RK), DIMENSION(:) ::  sources, boundaries
   END SUBROUTINE


   ! Generic interface for method that does post solve processing
   SUBROUTINE generic_var_post_solve(self, state)
      class(ModelVariable), intent(inout) :: self
      class(ModelState), intent(inout) :: state
   END SUBROUTINE

END MODULE  strat_statevar







!<    +---------------------------------------------------------------+
!     | Turbulence MODULE
!     | Updates state variables
!<    +---------------------------------------------------------------+
MODULE  strat_turbulence
   USE strat_kinds, only: RK
   USE MOD_Lake_Const
   USE strat_grid
   USE strat_simdata

   IMPLICIT NONE

!-----------------------------------------------------------
PRIVATE
!-----------------------------------------------------------

   ! Common Types
   type, public :: TurbulenceModule
      class(StaggeredGrid), pointer :: grid
      class(ModelConfig), pointer :: model_cfg
      class(ModelParam), pointer :: model_param

   CONTAINS
      procedure, pass :: init => turbulence_module_init
      procedure, pass :: update => turbulence_module_update

      procedure, pass :: do_production => turbulence_module_do_production
      procedure, pass :: do_seiche => turbulence_module_do_seiche

   END type

!-----------------------------------------------------------
CONTAINS
!-----------------------------------------------------------

   SUBROUTINE turbulence_module_init(self, grid, model_cfg, model_param)
      IMPLICIT NONE
      class(TurbulenceModule) :: self
      class(StaggeredGrid), target :: grid
      class(ModelConfig), target :: model_cfg
      class(ModelParam), target :: model_param

      self%grid => grid
      self%model_cfg => model_cfg
      self%model_param => model_param
   END SUBROUTINE


   SUBROUTINE turbulence_module_update(self, state, param)
      IMPLICIT NONE
      class(TurbulenceModule) :: self
      class(ModelState) :: state
      class(ModelParam) :: param

      CALL self%do_production(state)

      CALL self%do_seiche(state, param)

   END SUBROUTINE


   SUBROUTINE turbulence_module_do_production(self, state)
      IMPLICIT NONE
      class(TurbulenceModule) :: self
      class(ModelState) :: state

      associate (grid=>self%grid, &
                  ubnd_vol=>self%grid%ubnd_vol, &
                  ubnd_fce=>self%grid%ubnd_fce)

         ! Equation 5 (left) of Goudsmit, 2002
         ! P is defined on the inner faces
         state%P = 0
         state%P(2:ubnd_fce - 1) = (state%U(2:ubnd_vol) - state%U(1:ubnd_vol - 1))**2 + (state%V(2:ubnd_vol) - state%V(1:ubnd_vol - 1))**2
         state%P(2:ubnd_fce - 1) = state%P(2:ubnd_fce - 1)*state%num(2:ubnd_fce - 1)*grid%meanint(1:ubnd_vol - 1)**2
         ! Equation 5 (right) of Goudsmit, 2002
         state%B = 0
         state%B(2:ubnd_fce - 1) = -state%nuh(2:ubnd_fce - 1)*state%NN(2:ubnd_fce - 1)
         RETURN
      END associate
   END SUBROUTINE


   SUBROUTINE turbulence_module_do_seiche(self, state, param)
      IMPLICIT NONE
      class(TurbulenceModule) :: self
      class(ModelState) :: state
      class(ModelParam) :: param

      ! Local variables
      real(RK) :: W10, PS, PW, f_norm, minNN, a_seiche_local
      real(RK) :: distrib(self%grid%ubnd_fce)
      integer :: i

      associate (grid=>self%grid, &
                  ubnd_vol=>self%grid%ubnd_vol, &
                  ubnd_fce=>self%grid%ubnd_fce)
         minNN = 0
         distrib = 0

         ! Update distrib on inner faces
         DO i = 2, ubnd_fce - 1
            IF (state%NN(i)<0) state%NN(i) = 0
            distrib(i) = max(state%NN(i)**param%q_NN, minNN)/grid%Az(i)*grid%dAz(i - 1)
         ENDDO

         ! Determine a_seiche
         IF (self%model_cfg%split_a_seiche) THEN   ! IF a_seiche is splitted seasonally

            ! IF maximum stratification (N2) is higher than threshold
            IF (maxval(state%NN(2:ubnd_fce - 1)) >= param%strat_sumr) THEN
               a_seiche_local = param%a_seiche
            ! IF maximum stratification (N2) is lower than threshold
            ELSE IF (maxval(state%NN(2:ubnd_fce - 1)) < param%strat_sumr) THEN
               a_seiche_local = param%a_seiche_w
            ENDIF
         ELSE  ! IF a_seiche is not splitted
            a_seiche_local = param%a_seiche
         ENDIF

         ! Exit function IF a_seiche is 0
         IF (a_seiche_local == 0) THEN
            state%P_seiche = 0.0_RK
            RETURN
         ENDIF

         ! Calculate Seiche normalization factor
         f_norm = 0.0_RK
         IF (self%model_cfg%seiche_normalization == 1) THEN ! max NN
            f_norm = maxval(state%NN(2:ubnd_fce - 1))
            f_norm = (f_norm**param%q_NN)*grid%Az(ubnd_fce)*rho_0
         ELSE IF (self%model_cfg%seiche_normalization == 2) THEN ! integral
            DO i = 2, ubnd_fce - 1
               f_norm = f_norm + distrib(i)*grid%Az(i)*grid%h(i - 1)
            ENDDO
            f_norm = f_norm*rho_0
         ENDIF

         ! todo: direct float comparison...? OK?
         ! why is this code here?
         IF (f_norm == 0.) THEN
            DO i = 2, ubnd_fce - 1
               distrib(i) = 1/grid%h(i - 1)
            ENDDO
            f_norm = grid%Az(ubnd_fce)*rho_0
         ENDIF

         ! Adjust wind params based on configuration
         IF (self%model_cfg%use_filtered_wind) THEN !USE filtered wind (AG 2014)
            PW = a_seiche_local*grid%Az(ubnd_fce)*rho_air*state%C10*state%Wf**3
         ELSE ! USE real wind
            W10 = sqrt(state%u10**2 + state%v10**2)
            PW = a_seiche_local*grid%Az(ubnd_fce)*rho_air*state%C10*W10**3
         ENDIF

         ! Update E_Seiche
         PS = state%E_Seiche**(1.5_RK)*state%gamma
         state%E_Seiche = state%E_Seiche + (PW - PS)*state%dt

         ! Limit so that E_Seiche does not become negative
         IF (state%E_Seiche < 0.) THEN
            PS = (PS*state%dt + state%E_Seiche)/state%dt
            state%E_Seiche = 0.0_RK
         ENDIF

         ! Update P_Seiche
         ! Equation 24 in Goudsmit, 2002
         DO i = 2, ubnd_fce - 1
               state%P_Seiche(i) = 1.0_RK/f_norm*distrib(i)*PS*(1.0_RK - 10*sqrt(param%CD))
         ENDDO
         state%P_Seiche(1) = 0.0_RK
         state%P_Seiche(ubnd_fce) = 0.0_RK

         RETURN
      END associate
   END SUBROUTINE

END MODULE  strat_turbulence





!<    +---------------------------------------------------------------+
!     |  Absorption MODULE
!     |  - reads and processes absorption input file
!<    +---------------------------------------------------------------+
MODULE  strat_absorption
   USE strat_kinds, only: RK
   USE strat_simdata
   USE MOD_Lake_Const
   USE strat_grid
   USE strat_utilities
   IMPLICIT NONE

!-----------------------------------------------------------
PRIVATE
!-----------------------------------------------------------

   type, public :: AbsorptionModule
      class(ModelConfig), pointer :: cfg
      class(StaggeredGrid), pointer :: grid
      class(ModelParam), pointer :: param
      character(len=:), allocatable  :: file

      ! Variables that are used in between iteration. These used to be "save" variables
      real(RK) :: tb_start, tb_end !Start and END time
      real(RK), DIMENSION(:), allocatable :: z_absorb !Read depths
      real(RK), DIMENSION(:), allocatable :: absorb_start, absorb_end !Interpolated start and END values
      integer :: number_of_lines_read = 0
      integer :: eof, nval

   CONTAINS
      procedure, pass :: init => absorption_init
      procedure, pass :: update => absorption_update
   END type

!-----------------------------------------------------------
CONTAINS
!-----------------------------------------------------------

      SUBROUTINE absorption_init(self, model_config, model_param, grid)
         ! Initialize absorption
         IMPLICIT NONE
         class(AbsorptionModule) :: self
         class(StaggeredGrid), target :: grid
         class(ModelConfig), target :: model_config
         class(ModelParam), target :: model_param
      
         self%cfg => model_config
         self%param => model_param
         self%grid => grid
      
      
         ! Allocate arrays (used to be "save" variables)
         allocate (self%z_absorb(grid%max_length_input_data))
         allocate (self%absorb_start(grid%nz_grid))
         allocate (self%absorb_end(grid%nz_grid))
      END SUBROUTINE


      SUBROUTINE absorption_update(self, state)
      ! Update absorption
         IMPLICIT NONE
         class(AbsorptionModule) :: self
         class(ModelState) :: state
   
         ! Local Variables
         real(RK) :: dummy !Read depths
         real(RK) :: absorb_read_start(self%grid%max_length_input_data), absorb_read_end(self%grid%max_length_input_data) !Read start and END values
         integer :: i
   
         ! Associations for easier readability / comparability to old code
         associate (tb_start=>self%tb_start, &
                     tb_end=>self%tb_end, &
                     z_absorb=>self%z_absorb, &
                     absorb_start=>self%absorb_start, &
                     absorb_end=>self%absorb_end, &
                     eof=>self%eof, &
                     nval=>self%nval, &
                     nz=>self%grid%nz_occupied)
   
            IF (state%first_timestep) THEN ! First iteration
               open (30, status='old', file=self%file)
               IF (self%number_of_lines_read > 0) THEN
                  DO i = 1, self%number_of_lines_read
                     read (30, *, END=9) ! skip over already read lines
                  ENDDO
               ELSE
                  eof = 0
   
                  !Read depths: first line are names, second line is number of depths available
                  read (30, *, END=9)
                  CALL count_read(self)
                  read (30, *, END=9) nval
                  CALL count_read(self)
                  read (30, *, END=9) dummy, (z_absorb(i), i=1, nval)
                  CALL count_read(self)
   
                  !Make depths positives
                  DO i = 1, nval
                     z_absorb(i) = abs(z_absorb(i))
                  ENDDO
   
                  !Read first values
                  read (30, *, END=9) tb_start, (absorb_read_start(i), i=1, nval)
                  CALL count_read(self)
                  IF (state%datum < tb_start) CALL warn('First light attenuation date after simulation start time.')
   
                  !Interpolate absorb_read_start on z_absorb onto faces of grid
                  CALL self%grid%interpolate_to_face(z_absorb, absorb_read_start, nval, absorb_start)
   
                  read (30, *, END=7) tb_end, (absorb_read_end(i), i=1, nval)
                  CALL count_read(self)
   
                  ! Write to console that file was successfully read
                  CALL ok('Absorption input file successfully read')
   
                  ! DO the same for absorb_read_end
                  CALL self%grid%interpolate_to_face(z_absorb, absorb_read_end, nval, absorb_end)
               ENDIF
            ENDIF
   
            IF (state%datum <= tb_start .or. eof == 1) THEN !IF datum before first date or END of file reached
               goto 8
            ELSE
               DO WHILE (state%datum > tb_end) !Move to appropriate interval to get correct value
                  tb_start = tb_end
                  absorb_start(1:nz) = absorb_end(1:nz)
                  !Read next value
                  read (30, *, END=7) tb_end, (absorb_read_end(i), i=1, nval)
                  CALL count_read(self)
                  CALL self%grid%interpolate_to_face(z_absorb, absorb_read_end, nval, absorb_end)
               ENDDO
               !Linearly interpolate value at correct datum (for all depths)
               state%absorb(1:nz) = absorb_start(1:nz) + (state%datum - tb_start)/(tb_end - tb_start)*(absorb_end(1:nz) - absorb_start(1:nz))
            ENDIF
            RETURN
   
   7        eof = 1
            IF(state%datum>tb_start) CALL warn('Last light attenuation date before simulation END time.')
   
   8        state%absorb(1:nz) = absorb_start(1:nz)           !Take first value of current interval
            RETURN
   
   9        CALL error('Reading absorption file (no data found).')
   
         END associate
      END SUBROUTINE


   SUBROUTINE count_read(self)
      IMPLICIT NONE
      class(AbsorptionModule) :: self
      self%number_of_lines_read = self%number_of_lines_read + 1
   END SUBROUTINE

END MODULE  strat_absorption





!<    +---------------------------------------------------------------+
!     |  Advection MODULE
!     |  - Based on already read in/outflows, calculates Advection
!     |  - Might grow or shrink grid (methods merge/add_box)
!<    +---------------------------------------------------------------+
MODULE  strat_advection
   USE strat_kinds
   USE strat_simdata
   USE MOD_Lake_Const
   USE strat_grid
   USE strat_utilities
   IMPLICIT NONE
   PRIVATE

   type, public :: AdvectionModule
      class(ModelConfig), pointer :: cfg
      class(StaggeredGrid), pointer :: grid
      class(ModelParam), pointer :: param

   CONTAINS
      procedure, pass :: init => advection_init
      procedure, pass :: update => advection_update
      procedure, pass :: merge_box => advection_merge_box
      procedure, pass :: add_box => advection_add_box
   END type

CONTAINS

   SUBROUTINE advection_init(self, model_config, model_param, grid)
      IMPLICIT NONE
      class(AdvectionModule) :: self
      class(StaggeredGrid), target :: grid
      class(ModelConfig), target :: model_config
      class(ModelParam), target :: model_param

      self%cfg => model_config
      self%param => model_param
      self%grid => grid

   END SUBROUTINE

   ! A lot of code that is hard to test - might be refactored in the future
   SUBROUTINE advection_update(self, state)
      IMPLICIT NONE
      class(AdvectionModule) :: self
      class(ModelState) :: state

      real(RK) :: top_z, top_h
      real(RK) :: top
      real(RK) :: dh, dh_i(1:2), h_div_2, h_mult_2 ! depth differences
      real(RK) :: dU(self%grid%nz_grid), dV(self%grid%nz_grid), dTemp(self%grid%nz_grid), dS(self%grid%nz_grid)
      real(RK) :: dt_i(1:2) ! first and second time step
      real(RK) :: AreaFactor_adv(1:self%grid%nz_grid)
      integer :: i, t_i, outflow_above, outflow_below

      associate (grid=>self%grid, &
               nz_occupied=>self%grid%nz_occupied, &
               dt=>state%dt, &
               h=>self%grid%h, &
               Q_vert=>state%Q_vert, &
               ubnd_vol=>self%grid%ubnd_vol, &
               ubnd_fce=>self%grid%ubnd_fce)

         !Depth difference compared to previous timestep
         top_z = grid%z_face(ubnd_fce)
         top_h = grid%h(ubnd_vol)

         dh = state%Q_vert(ubnd_fce)/grid%Az(ubnd_fce)*state%dt
         h_div_2 = 0.5_RK*h(nz_occupied - 1) ! Take second highest box since the top box might not be at the full height
         h_mult_2 = 2_RK*h(nz_occupied - 1)

         ! Calculate timestep splitting
         !Split timestep depending on situation
         IF (dh == 0.) THEN ! IF volume does not change, take one normal time step
            dt_i(1) = dt
         ELSE IF (top_z == grid%max_depth) THEN ! IF we are already at the maximum lake level
            dt_i(1) = dt
         ELSE IF ((dh + top_z) >= grid%max_depth) THEN ! IF the full timestep would lead to a lake level higher than maximum allowed lake level, split the timestep.
            dt_i(1) = (grid%max_depth - top_z)/dh*dt
         ELSE IF (((dh + top_h) > h_div_2) .and. & ! IF top box>0.5*lower box and <2*lower box, take one time step
                  ((dh + top_h) < h_mult_2)) THEN
            dt_i(1) = dt
         ELSE IF ((dh + top_h) <= h_div_2) THEN ! IF top box<=0.5*lower box, first step until top box=0.5*lower box
            dt_i(1) = abs((top_h - h_div_2)/dh)*dt
         ELSE ! IF top box>=2*lower box, first step until top box = 2*lower box
            dt_i(1) = abs((2*h(nz_occupied - 1) - top_h)/dh)*dt
         ENDIF
         dt_i(2) = dt - dt_i(1) ! Rest of timestep

         ! FB 2016/2019: Revisions
         DO t_i = 1, 2 !First and (IF needed) second timestep
            AreaFactor_adv(1:nz_occupied) = dt_i(t_i)/((grid%Az(1:nz_occupied) + grid%Az(2:nz_occupied+1))/2*grid%h(1:nz_occupied)) ! Area factor for dt(t_i)
            dh_i(t_i) = dh*dt_i(t_i)/dt ! Depth difference for dt(t_i)

            ! Calculate changes
            DO i = 1, ubnd_vol
               ! For the top-most cell, IF Q_vert at the upper face is positive, there is still no outflow (the cell is simply growing, but this is done elsewhere)
               IF ((i == ubnd_vol) .and. Q_vert(i + 1) > 0) THEN
                  top = 0
               ELSE
                  top = 1
               ENDIF

               ! IF Q_vert at the upper face of cell i is positive, THEN there is outflow to the cell above
               IF (Q_vert(i + 1) > 0) THEN
                  outflow_above = 1
               ELSE 
                  outflow_above = 0
               ENDIF

               ! IF Q_vert at the lower face of cell i is negative, THEN there is outflow to the cell below
               IF (Q_vert(i) < 0) THEN
                  outflow_below = 1
               ELSE
                  outflow_below = 0
               ENDIF

               ! Calculate advective flow out of cell (thus negative sign in the front) i to the cells above and below
               dU(i) = -(top*outflow_above*Q_vert(i + 1) - outflow_below*Q_vert(i))*state%U(i)
               dV(i) = -(top*outflow_above*Q_vert(i + 1) - outflow_below*Q_vert(i))*state%V(i)
               dTemp(i) = -(top*outflow_above*Q_vert(i + 1) - outflow_below*Q_vert(i))*state%T(i)
               dS(i) = -(top*outflow_above*Q_vert(i + 1) - outflow_below*Q_vert(i))*state%S(i)

               ! Calculate the advective flow into cell i from below
               IF (i > 1 .and. Q_vert(i ) > 0) THEN
                  dU(i) = dU(i) + Q_vert(i)*state%U(i - 1)
                  dV(i) = dV(i) + Q_vert(i)*state%V(i - 1)
                  dTemp(i) = dTemp(i) + Q_vert(i)*state%T(i - 1)
                  dS(i) = dS(i) + Q_vert(i)*state%S(i - 1)
               ENDIF

               ! Calculate the advective flow into cell i from above (- sign in front because Q_vert is negative IF there is inflow)
               IF (i < ubnd_vol .and. Q_vert(i + 1) < 0) THEN
                  dU(i) = dU(i) - Q_vert(i + 1)*state%U(i + 1)
                  dV(i) = dV(i) - Q_vert(i + 1)*state%V(i + 1)
                  dTemp(i) = dTemp(i) - Q_vert(i + 1)*state%T(i + 1)
                  dS(i) = dS(i) - Q_vert(i + 1)*state%S(i + 1)
               ENDIF
            ENDDO

            ! Add change to state variables
            ! dT = dT(vertical advection) + dT(inflow) + dT(outflow), units: °C*m^3/s
            dTemp(1:ubnd_vol) = dTemp(1:ubnd_vol) + state%Q_inp(3, 1:ubnd_vol) + state%Q_inp(2, 1:ubnd_vol)*state%T(1:ubnd_vol)
            ! dS = dS(vertical advection) + dS(inflow) + dS(outflow), units: ‰*m^3/s
            dS(1:ubnd_vol) = dS(1:ubnd_vol) + state%Q_inp(4, 1:ubnd_vol) + state%Q_inp(2, 1:ubnd_vol)*state%S(1:ubnd_vol)

            ! Add change to the state variable
            state%U(1:ubnd_vol) = state%U(1:ubnd_vol) + AreaFactor_adv(1:ubnd_vol)*dU(1:ubnd_vol)
            state%V(1:ubnd_vol) = state%V(1:ubnd_vol) + AreaFactor_adv(1:ubnd_vol)*dV(1:ubnd_vol)
            state%T(1:ubnd_vol) = state%T(1:ubnd_vol) + AreaFactor_adv(1:ubnd_vol)*dTemp(1:ubnd_vol)
            state%S(1:ubnd_vol) = state%S(1:ubnd_vol) + AreaFactor_adv(1:ubnd_vol)*dS(1:ubnd_vol)

            ! Variation of variables due to change in volume
            state%U(ubnd_vol) = state%U(ubnd_vol)*h(ubnd_vol)/(h(ubnd_vol) + dh_i(t_i))
            state%V(ubnd_vol) = state%V(ubnd_vol)*h(ubnd_vol)/(h(ubnd_vol) + dh_i(t_i))
            state%T(ubnd_vol) = state%T(ubnd_vol)*h(ubnd_vol)/(h(ubnd_vol) + dh_i(t_i))
            state%S(ubnd_vol) = state%S(ubnd_vol)*h(ubnd_vol)/(h(ubnd_vol) + dh_i(t_i))

            ! Adjust boxes (Horrible IF/ELSE construction - replace!)
            IF (t_i == 1) THEN
               IF (dh == 0) THEN ! IF volume does not change, RETURN
                  RETURN
               ELSE IF ((dh + top_z) >= grid%max_depth) THEN ! IF surface level reached
                  CALL grid%modify_top_box(grid%max_depth - top_z)
                  RETURN
               ELSE IF (((dh_i(t_i) + top_h) > h_div_2) .and. &
                        ((dh_i(t_i) + top_h) < (h_mult_2))) THEN ! and top box<2*lower box
                  CALL grid%modify_top_box(dh_i(t_i))
                  RETURN
               ELSE IF (t_i == 1 .and. (dh + top_h) <= h_div_2) THEN ! IF top box<=0.5*lower box, merge 2 boxes
                  CALL self%merge_box(state, dh_i(t_i))
               ELSE IF (t_i == 1 .and. (dh + top_h) >= h_mult_2) THEN ! IF top box>=2*lower box, add one box
                  CALL self%add_box(state, dh_i(t_i))
               ENDIF ! dh==0
            ENDIF

         ENDDO !ENDDO t_i=1,2
      END associate
   END SUBROUTINE

   ! Merges two boxes
   ! - Takes care of calculating the new state variable for this box
   ! - Calls grid methods to modify grid spacing etc
   SUBROUTINE advection_merge_box(self, state, dh)
       IMPLICIT NONE
       class(AdvectionModule) :: self
       class(ModelState) :: state
       real(RK) :: dh
       real(RK) :: w_a, w_b
       associate (ubnd_fce=>self%grid%ubnd_fce, ubnd_vol=>self%grid%ubnd_vol)

         ! New values of the state variables are weighted averages
         !determine weighting an normalization connstant
         w_a = 0.5_RK*self%grid%Az(ubnd_fce)
         w_b = self%grid%Az(ubnd_fce - 1)

         ! shrink grid by one (this also updates ubnd_fce/vol)
         CALL self%grid%shrink(dh)

         ! update quantities in new top box (based on former top box and current value)
         state%U(ubnd_vol) = (w_a*state%U(ubnd_vol + 1) + w_b*state%U(ubnd_vol))/(w_a + w_b)
         state%V(ubnd_vol) = (w_a*state%V(ubnd_vol + 1) + w_b*state%V(ubnd_vol))/(w_a + w_b)
         state%T(ubnd_vol) = (w_a*state%T(ubnd_vol + 1) + w_b*state%T(ubnd_vol))/(w_a + w_b)
         state%S(ubnd_vol) = (w_a*state%S(ubnd_vol + 1) + w_b*state%S(ubnd_vol))/(w_a + w_b)

         state%k(ubnd_fce) = (w_a*state%k(ubnd_fce + 1) + w_b*state%k(ubnd_fce))/(w_a + w_b)
         state%eps(ubnd_fce) = (w_a*state%eps(ubnd_fce + 1) + w_b*state%eps(ubnd_fce))/(w_a + w_b)
         state%Q_vert(ubnd_fce) = (w_a*state%Q_vert(ubnd_fce + 1) + w_b*state%Q_vert(ubnd_fce))/(w_a + w_b)

         ! update area factors
         CALL self%grid%update_area_factors()

       END associate
   END SUBROUTINE

   ! Adds a new box
   ! - Takes care of calculating the new state variable for this box
   ! - Calls grid methods to modify grid spacing etc
   SUBROUTINE advection_add_box(self, state, dh)
      IMPLICIT NONE
      class(AdvectionModule) :: self
      class(ModelState) :: state
      real(RK) :: dh
      associate (ubnd_fce=>self%grid%ubnd_fce, ubnd_vol=>self%grid%ubnd_vol)

         ! extend grid by one (also updates ubnd_vol etc)
         CALL self%grid%grow(dh)

         ! Update quantities in new grid element
         state%U(ubnd_vol) = state%U(ubnd_vol - 1)
         state%V(ubnd_vol) = state%V(ubnd_vol - 1)
         state%T(ubnd_vol) = state%T(ubnd_vol - 1)
         state%S(ubnd_vol) = state%S(ubnd_vol - 1)
         state%Q_vert(ubnd_fce) = state%Q_vert(ubnd_fce - 1) ! Vertical discharge of new box

         state%k(ubnd_fce) = state%k(ubnd_fce - 1)
         state%eps(ubnd_fce) = state%eps(ubnd_fce - 1)

         CALL self%grid%update_area_factors()

      END associate
    END SUBROUTINE
 
 END MODULE strat_advection
 




!<    +---------------------------------------------------------------+
!     |  Lateral MODULE
!     |     Reads and processes inflows/outflows such that
!     |     Advection can be calculated in the next step
!     |     There are at least two different implementations possible:
!     |         - LateralRhoModule : Inflow plunges according to density
!     |         - LateralModule:     Inflow affects layer as configured in file
!<    +---------------------------------------------------------------+
MODULE  strat_lateral
   USE strat_kinds
   USE strat_simdata
   USE MOD_Lake_Const
   USE strat_grid
   USE strat_utilities
   IMPLICIT NONE

!-----------------------------------------------------------
PRIVATE 
!-----------------------------------------------------------

   ! Generic base class for both modules (LateralRho and Normal)
   type, abstract, public :: GenericLateralModule
      class(ModelConfig), pointer :: cfg
      class(StaggeredGrid), pointer :: grid
      class(ModelParam), pointer :: param

      ! Variables that where either marked with "save" before, or that have been
      ! global, but only used in the lateral environment:
      real(RK), DIMENSION(:, :), allocatable   :: z_Inp, Q_start, Qs_start, Q_end, Qs_end, Q_read_start, Q_read_end
      real(RK), DIMENSION(:, :), allocatable   :: Inp_read_start, Inp_read_end, Qs_read_start, Qs_read_end
      real(RK) :: tb_start(1:4), tb_end(1:4) ! Input depths, start time, END time
      integer :: number_of_lines_read(1:4) = 0
      integer :: eof(1:4)
      integer :: nval(1:4), nval_deep(1:4), nval_surface(1:4) ! Number of values

   CONTAINS
      procedure, pass :: init => lateral_generic_init
      procedure(lateral_generic_update), deferred, pass :: update
   END type

   ! Subclasses
   type, extends(GenericLateralModule), public :: LateralRhoModule
   CONTAINS
      procedure, pass, public :: update => lateral_rho_update
   END type

   type, extends(GenericLateralModule), public:: LateralModule
   CONTAINS
      procedure, pass, public :: update => lateral_update
   END type
 
!-----------------------------------------------------------
CONTAINS   
!-----------------------------------------------------------
   SUBROUTINE lateral_generic_update(self, state)
      IMPLICIT NONE
      class(GenericLateralModule) :: self
      class(ModelState) :: state
   END SUBROUTINE

   SUBROUTINE lateral_generic_init(self, model_config, model_param, grid)
      IMPLICIT NONE
      class(GenericLateralModule) :: self
      class(StaggeredGrid), target :: grid
      class(ModelConfig), target :: model_config
      class(ModelParam), target :: model_param

      self%cfg => model_config
      self%param => model_param
      self%grid => grid
   END SUBROUTINE

   ! Implementation for lateral rho
   SUBROUTINE lateral_rho_update(self, state)
      IMPLICIT NONE
      class(LateralRhoModule) :: self
      class(ModelState) :: state

      ! Local Declarations
      real(RK) :: Inp(1:4,1:state%nz_input)
      real(RK) :: dummy
      real(RK) :: Q_in(1:self%grid%ubnd_vol), h_in(1:self%grid%ubnd_vol)
      real(RK) :: T_in, S_in, rho_in, CD_in, g_red, slope, Ri, E, Q_inp_inc

      integer :: i, j, k, i1, i2, l
      integer :: fnum(1:4) ! File number
      logical :: at_start(1:4) = .true.
      character(len=20) :: fname(1:4) ! File name
      logical :: first

      associate (datum=>state%datum, &
               idx=>state%first_timestep, &
               number_of_lines_read=>self%number_of_lines_read, &
               Q_inp=>state%Q_inp, & ! Q_inp is the input at each depth for each time step
               Q_vert=>state%Q_vert, & ! Q_vert is the integrated net water input at each depth (integrated inflow - outflow)
               grid=>self%grid, &
               ubnd_vol=>self%grid%ubnd_vol, &
               ubnd_fce=>self%grid%ubnd_fce)

         fname = ['inflow           ','outflow          ','input temperature','input salinity   ']
         fnum = [41,42,43,44]

         CALL detect_first(self, state, first)
         DO i = 1, 4 ! DO this for inflow, outflow, temperature and salinity
            IF (first) THEN ! First iteration
               IF (self%number_of_lines_read(i) == 0) THEN
                  ! Default values
                  self%Q_start(i,:) = 0.0_RK
                  self%Q_end(i,:) = 0.0_RK
                  self%Qs_start(i,:) = 0.0_RK
                  self%Qs_end(i,:) = 0.0_RK

                  ! END of file is not reached
                  self%eof(i) = 0

                  ! Read input depths
                  read(fnum(i),*,END=9)
                  CALL count_read(self, i)
                  ! Check whether there are any surface inputs, otherwise assume that all columns are for deep inputs (backwards compatibility)
                  IF (state%has_surface_input(i)) THEN
                     ! Read number of deep and surface columns
                     read(fnum(i), *, END=9) self%nval_deep(i), self%nval_surface(i)
                     CALL count_read(self, i)
                     ! Total number of values to read
                     self%nval(i) = self%nval_deep(i) + self%nval_surface(i)
                     ! Read input depths
                     read(fnum(i),*,END=9) dummy, (self%z_Inp(i,j),j=1,self%nval(i))
                     CALL count_read(self, i)
                     ! Convert deep input depths
                     self%z_Inp(i,1:self%nval_deep(i)) = grid%z_zero + self%z_Inp(i,1:self%nval_deep(i))
                     ! Convert surface input depths
                     self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)) = grid%lake_level + self%z_Inp(i, self%nval_deep(i) + 1 :self%nval(i))
                  ELSE
                     ! Read number of deep columns
                     read(fnum(i),*,END=9) self%nval_deep(i)
                     CALL count_read(self, i)
                     ! Total number of values to read
                     self%nval(i) = self%nval_deep(i)
                     ! Read input depths
                     read(fnum(i),*,END=9) dummy, (self%z_Inp(i,j),j=1,self%nval(i))
                     CALL count_read(self, i)
                     ! Convert deep input depths
                     self%z_Inp(i,1:self%nval(i)) = grid%z_zero + self%z_Inp(i,1:self%nval(i))
                  ENDIF

                  !Read first input values
                  read(fnum(i),*,END=9) self%tb_start(i),(self%Inp_read_start(i,j),j=1,self%nval(i))
                  CALL count_read(self, i)

                  ! IF there is deep outflow (i==2)
                  IF (i==2 .and. state%has_deep_input(i)) THEN
                     CALL Integrate(self%z_Inp(i,1:self%nval_deep(i)),self%Inp_read_start(i,1:self%nval_deep(i)),self%Q_read_start(i,:),self%nval_deep(i))
                     CALL grid%interpolate_to_face_from_second(self%z_Inp(i,1:self%nval_deep(i)),self%Q_read_start(i,:),self%nval_deep(i),self%Q_start(i,:))
                  ENDIF
                  ! IF there is any surface inflow
                  IF (state%has_surface_input(i)) THEN
                     CALL Integrate(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Inp_read_start(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_start(i, :), self%nval_surface(i))
                     CALL grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_start(i, :), self%nval_surface(i), self%Qs_start(i, :))
                  ENDIF


                  ! Read next line
                  read(fnum(i),*,END=7) self%tb_end(i),(self%Inp_read_end(i,j),j=1,self%nval(i))
                  CALL count_read(self, i)

                  ! IF there is deep outflow (i==2)
                  IF (i==2 .and. state%has_deep_input(i)) THEN
                     CALL Integrate(self%z_Inp(i,1:self%nval_deep(i)),self%Inp_read_end(i,1:self%nval_deep(i)),self%Q_read_end(i,:),self%nval_deep(i))
                     CALL grid%interpolate_to_face_from_second(self%z_Inp(i,1:self%nval_deep(i)),self%Q_read_end(i,:),self%nval_deep(i),self%Q_end(i,:))
                  ENDIF

                  ! IF there is any surface inflow
                  IF (state%has_surface_input(i)) THEN
                     CALL Integrate(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Inp_read_end(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i))
                     CALL grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i), self%Qs_end(i, :))
                  ENDIF
                  CALL ok('Input file successfully read: '//fname(i))
               ENDIF
            ELSE ! first
               IF (at_start(i)) THEN
                  DO l = 1, self%number_of_lines_read(i)
                     read (fnum(i), *, END=9) ! Skip already read and processed lines
                  ENDDO
               ENDIF
            ENDIF
            at_start(i) = .false.

            ! IF lake level changes and IF there is surface inflow, adjust inflow depth to keep relative inflow depth constant
            IF ((.not. grid%lake_level == grid%lake_level_old) .and. state%has_surface_input(i)) THEN

               ! Readjust surface input depths
               self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)) = self%z_Inp(i, self%nval_deep(i) + 1 :self%nval(i)) - grid%lake_level_old + grid%lake_level

               ! Adjust surface inflow to new lake level
               IF (state%has_surface_input(i)) THEN
                  CALL grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_start(i, :), self%nval_surface(i), self%Qs_start(i, :))
                  CALL grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i), self%Qs_end(i, :))
               ENDIF
            ENDIF ! ENDIF not lake_level...


            IF ((datum<=self%tb_start(i)).or.(self%eof(i)==1)) THEN    ! IF datum before first date or END of file reached
               goto 8
            ELSE
               DO WHILE (.not.((datum>=self%tb_start(i)).and.(datum<=self%tb_end(i)))) ! DO until datum between dates
                  ! Move one step in time
                  self%tb_start(i) = self%tb_end(i)
                  self%Qs_start(i, :) = self%Qs_end(i, :)
                  self%Qs_read_start(i, :) = self%Qs_read_end(i, :)

                  ! For outflow, take Q_start; for inflow, temperature and salinity, take Inp for the plunging algorithm
                  IF (i==2) THEN
                     self%Q_start(i, :) = self%Q_end(i, :)
                  ELSE
                     self%Inp_read_start(i,1:self%nval_deep(i)) = self%Inp_read_end(i,1:self%nval_deep(i))
                  ENDIF

                  ! Read next line
                  read(fnum(i),*,END=7) self%tb_end(i),(self%Inp_read_end(i,j),j=1,self%nval(i))
                  CALL count_read(self, i)

                  ! IF there is deep outflow (i==2)
                  IF (i==2 .and. state%has_deep_input(i)) THEN
                     CALL Integrate(self%z_Inp(i,1:self%nval_deep(i)),self%Inp_read_end(i,1:self%nval_deep(i)),self%Q_read_end(i,:),self%nval_deep(i))
                     CALL grid%interpolate_to_face_from_second(self%z_Inp(i,1:self%nval_deep(i)),self%Q_read_end(i,:),self%nval_deep(i),self%Q_end(i,:))
                  ENDIF

                  ! IF there is any surface inflow
                  IF (state%has_surface_input(i)) THEN
                     CALL Integrate(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Inp_read_end(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i))
                     CALL grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i), self%Qs_end(i, :))
                  ENDIF
               ENDDO

               IF(self%tb_end(i)<=self%tb_start(i)) THEN
                  CALL error('Dates in '//trim(fname(i))//' file must always be increasing.')
               ENDIF

               ! Linearly interpolate value at correct datum
               IF (i/=2) THEN
                  ! For plunging input, Inp will be needed later
                  Inp(i,1:self%nval_deep(i)) = self%Inp_read_start(i,1:self%nval_deep(i)) + (datum-self%tb_start(i))/(self%tb_end(i) - self%tb_start(i))* &
                  (self%Inp_read_end(i,1:self%nval_deep(i)) - self%Inp_read_start(i,1:self%nval_deep(i)))

                  ! Surface input is already added to Q_inp; the plunging algorithm will add the deep input further below
                  DO j=1,ubnd_fce
                     Q_inp(i,j) = (self%Qs_start(i,j)) + (datum - self%tb_start(i))/(self%tb_end(i) - self%tb_start(i))* &
                     (self%Qs_end(i,j) - self%Qs_start(i,j))
                  ENDDO
               ELSE
                  ! For outflow (i==2), both surface and deep inputs are added to Q_inp
                  DO j=1,ubnd_fce
                     Q_inp(i,j) = (self%Q_start(i,j) + self%Qs_start(i,j)) + (datum-self%tb_start(i))/(self%tb_end(i)-self%tb_start(i))* &
                     (self%Q_end(i,j) + self%Qs_end(i,j) - self%Q_start(i,j) - self%Qs_start(i,j))
                  ENDDO
               ENDIF
            ENDIF
            goto 11

7           self%eof(i) = 1
8           IF(i/=2) Inp(i,:) = self%Inp_read_start(i,:) ! Set to closest available value
            Q_inp(i,1:ubnd_vol) = self%Q_start(i,1:ubnd_vol) + self%Qs_start(i,1:ubnd_vol) ! Set to closest available value
            goto 11

9           write(6,*) '[WARNING] ','No data found in ',trim(fname(i)),' file. Check number of depths. Values set to zero.'
            self%eof(i) = 1
            IF(i/=2) Inp(i,1:self%nval_deep(i)) = 0.0_RK
            IF(i/=2) self%Inp_read_start(i,1) = 0.0_RK
            IF(i==2) Q_inp(i,1:ubnd_vol) = 0.0_RK
            IF(i==2) self%Q_start(i,1:ubnd_fce) = 0.0_RK
            self%Qs_start(i,1:ubnd_fce) = 0.0_RK

11          continue
         ENDDO      ! ENDDO i=1,4

         !Set Q_inp to the differences (from the integrals)
         DO i = 1, 4
            DO j = 1, ubnd_vol
               Q_inp(i, j) = Q_inp(i, j + 1) - Q_inp(i, j)
            ENDDO
            Q_inp(i,ubnd_vol + 1) = 0
         ENDDO

         ! Plunging algorithm
         DO j = 1,self%nval_deep(1)  ! nval_deep needs to be the same for i=1,3,4
            IF (Inp(1,j) > 1E-15) THEN
               k = ubnd_vol
               DO WHILE (grid%z_volume(k) > self%z_Inp(1,j)) ! Find the place where the plunging inflow enters the lake
                  k = k - 1
               ENDDO
               Q_in(k) = Inp(1,j) !Inflow flow rate [m3/s]
               T_in = Inp(3,j) !Inflow temperature [°C]
               S_in = Inp(4,j) !Inflow salinity [‰]
               rho_in = rho_0*(0.9998395 + T_in*(6.7914e-5 + T_in*(-9.0894e-6 + T_in*&
                     (1.0171e-7 + T_in*(-1.2846e-9 + T_in*(1.1592e-11 + T_in*(-5.0125e-14)))))) + &
                     (8.181e-4 + T_in*(-3.85e-6 + T_in*(4.96e-8)))*S_in) !Inflow density [kg/m3]
               g_red = g*(rho_in - state%rho(k))/rho_in !Reduced gravity [m/s2]

               slope = pi/72 !Slope of inflow
               !hang = pi/3 !Stream half-angle
               CD_in = self%param%CD*10 !Inflow drag coefficient
               !Ri = CD_in*(1+0.21*CD_in**0.5*sin(hang))/(sin(hang)*tan(slope)) !Richardson number
               !Ri = CD_in/tan(slope)*(1/sin(hang)+0.21*CD_in**0.5) !Richardson number
               Ri = CD_in/tan(slope)*(1.15 + 0.21*CD_in**0.5) !Richardson number (assuming an inflow half-angle of pi/3)
               E = 1.6*CD_in**1.5/Ri !Entrainment coefficient
               h_in(k) = (2*Q_in(k)**2*Ri*tan(slope)**2/abs(g_red))**0.2 !Inflow thickness [m]

               IF (g_red > 0) THEN !Inflow plunges
                  DO WHILE ((rho_in > state%rho(k)).and.(k > 1))
                     h_in(k - 1) = 1.2*E*(grid%z_volume(k) - grid%z_volume(k-1))/sin(slope) + h_in(k)
                     Q_in(k - 1) = Q_in(k)*(h_in(k - 1)/h_in(k))**(5./3.)
                     Q_inp(2,k) = Q_inp(2,k) - (Q_in(k-1) - Q_in(k))
                     T_in = (T_in*Q_in(k) + state%T(k)*(Q_in(k - 1) - Q_in(k)))/Q_in(k - 1)
                     S_in = (S_in*Q_in(k) + state%S(k)*(Q_in(k - 1) - Q_in(k)))/Q_in(k - 1)
                     rho_in = (rho_in*Q_in(k) + state%rho(k)*(Q_in(k - 1) - Q_in(k)))/Q_in(k - 1)
                     k = k - 1
                  ENDDO
                  i2 = k
                  DO i1 = k,ubnd_vol !extend upwards
                     IF(i1 == ubnd_vol) exit
                     IF(grid%z_volume(i1 + 1) > (grid%z_volume(k) + h_in(k))) exit
                  ENDDO
               ELSE IF (g_red < 0) THEN !Inflow rises
                  DO WHILE ((rho_in < state%rho(k)) .and. (k < ubnd_vol))
                     h_in(k + 1) = 1.2*E*(grid%z_volume(k + 1) - grid%z_volume(k))/sin(slope) + h_in(k)
                     Q_in(k + 1) = Q_in(k)*(h_in(k + 1)/h_in(k))**(5./3.)
                     Q_inp(2,k) = Q_inp(2,k) - (Q_in(k + 1) - Q_in(k))
                     T_in = (T_in*Q_in(k) + state%T(k)*(Q_in(k + 1) - Q_in(k)))/Q_in(k + 1)
                     S_in = (S_in*Q_in(k) + state%S(k)*(Q_in(k + 1) - Q_in(k)))/Q_in(k + 1)
                     rho_in = (rho_in*Q_in(k) + state%rho(k)*(Q_in(k + 1) - Q_in(k)))/Q_in(k + 1)
                     k = k + 1
                  ENDDO
                  i1 = k
                  DO i2 = k,1,-1 !extend downwards
                     IF(i2 == 1) exit
                     IF(grid%z_volume(i2 - 1) < (grid%z_volume(k) - h_in(k))) exit
                  ENDDO
               ENDIF

               ! Deep plunging input is added to Q_inp for i=1,3,4 (inflow, temperature, salinity)
               DO i = i2,i1
                  Q_inp_inc = Q_in(k)/(grid%z_face(i1 + 1) - grid%z_face(i2))*grid%h(i)
                  Q_inp(1,i) = Q_inp(1,i) + Q_inp_inc
                  Q_inp(3,i) = Q_inp(3,i) + T_in*Q_inp_inc
                  Q_inp(4,i) = Q_inp(4,i) + S_in*Q_inp_inc
               ENDDO
            ENDIF
         ENDDO

         ! Q_vert is the integrated difference between in- and outflow (starting at the lake bottom)
         ! Q_vert is located on the face grid, m^3/s
         Q_vert(1) = 0
         DO i = 2,ubnd_fce
            Q_vert(i) = Q_vert(i - 1) + Q_inp(1,i - 1) + Q_inp(2,i - 1)
         ENDDO
      END associate
   END SUBROUTINE

   ! "Normal" Implementation
   SUBROUTINE lateral_update(self, state)
      IMPLICIT NONE
      class(LateralModule) :: self
      class(ModelState) :: state

      ! Local Declarations
      real(RK) :: dummy

      integer :: i, j, l
      integer :: fnum(1:4) ! File number
      logical :: at_start(1:4) = .true.
      character(len=20) :: fname(1:4)
      logical :: first

      associate (datum=>state%datum, &
               idx=>state%first_timestep, &
               number_of_lines_read=>self%number_of_lines_read, &
               Q_inp=>state%Q_inp, & ! Q_inp is the input at each depth for each time step
               Q_vert=>state%Q_vert, & ! Q_vert is the integrated net water input
               grid=>self%grid, &
               ubnd_vol=>self%grid%ubnd_vol, &
               ubnd_Fce=>self%grid%ubnd_fce)

         fname = ['inflow           ', 'outflow          ', 'input temperature', 'input salinity   ']
         fnum = [41, 42, 43, 44]

         ! FB 2016: Major revision to include surface inflow
         CALL detect_first(self, state, first)
         DO i = 1, 4 ! DO this for inflow, outflow, temperature and salinity
            IF (first) THEN ! First iteration
               IF (self%number_of_lines_read(i) == 0) THEN
                  ! Default values
                  self%Q_start(i,:) = 0.0_RK
                  self%Q_end(i,:) = 0.0_RK
                  self%Qs_start(i, :) = 0.0_RK
                  self%Qs_end(i, :) = 0.0_RK

                  ! Open file and start to read
                  self%eof(i) = 0
                  read (fnum(i), *, END=9) ! Skip first row: description of columns
                  CALL count_read(self, i)

                  IF (state%has_surface_input(i)) THEN
                     ! Read number of deep and surface inflows
                     read (fnum(i), *, END=9) self%nval_deep(i), self%nval_surface(i)
                     CALL count_read(self, i)
                     ! Total number of values to read
                     self%nval(i) = self%nval_deep(i) + self%nval_surface(i)
                  ELSE
                     read (fnum(i), *, END=9) self%nval_deep(i)
                     CALL count_read(self, i)
                     ! Total number of values to read
                     self%nval(i) = self%nval_deep(i)
                  ENDIF

                  ! Read input depths
                  read (fnum(i), *, END=9) dummy, (self%z_Inp(i, j), j=1, self%nval(i))
                  CALL count_read(self, i)

                  ! Convert input depths
                  self%z_Inp(i, 1:self%nval_deep(i)) = grid%z_zero + self%z_Inp(i, 1:self%nval_deep(i))

                  IF (state%has_surface_input(i)) THEN
                     ! Convert surface input depths
                     self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)) = grid%lake_level + self%z_Inp(i, self%nval_deep(i) + 1 :self%nval(i))
                  ENDIF

                  ! Read first input line
                  read (fnum(i), *, END=9) self%tb_start(i), (self%Inp_read_start(i, j), j=1, self%nval(i))
                  CALL count_read(self, i)

                  IF (state%has_deep_input(i)) THEN
                     ! Cumulative integration of input
                     CALL Integrate(self%z_Inp(i, :), self%Inp_read_start(i, :), self%Q_read_start(i, :), self%nval_deep(i))
                     ! Interpolation on face grid
                     CALL grid%interpolate_to_face_from_second(self%z_Inp(i, :), self%Q_read_start(i, :), self%nval_deep(i), self%Q_start(i, :))
                  ENDIF

                  ! IF there is surface input, integrate and interpolate
                  IF (state%has_surface_input(i)) THEN
                     CALL Integrate(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Inp_read_start(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_start(i, :), self%nval_surface(i))
                     CALL grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_start(i, :), self%nval_surface(i), self%Qs_start(i, :))
                  ENDIF


                  ! Read second line and treatment of deep inflow
                  read (fnum(i), *, END=7) self%tb_end(i), (self%Inp_read_end(i, j), j=1, self%nval(i))
                  CALL count_read(self, i)
                  IF (state%has_deep_input(i)) THEN
                     CALL Integrate(self%z_Inp(i, :), self%Inp_read_end(i, :), self%Q_read_end(i, :), self%nval_deep(i))
                     CALL grid%interpolate_to_face_from_second(self%z_Inp(i, :), self%Q_read_end(i, :), self%nval_deep(i), self%Q_end(i, :))
                  ENDIF
                  ! IF there is surface input, integrate and interpolate
                  IF (state%has_surface_input(i)) THEN
                     CALL Integrate(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Inp_read_end(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i))
                     CALL grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i), self%Qs_end(i, :))
                  ENDIF

                  CALL ok('Input file successfully read: '//fname(i))
               ENDIF
            ELSE ! first
               IF (at_start(i)) THEN
                  DO l = 1, self%number_of_lines_read(i)
                     read (fnum(i), *, END=9) ! Skip already read and processed lines
                  ENDDO
               ENDIF
            ENDIF
            at_start(i) = .false.

            ! IF lake level changes and IF there is surface inflow, adjust inflow depth to keep them at the surface
            IF ((.not. grid%lake_level == grid%lake_level_old) .and. (state%has_surface_input(i))) THEN

               ! Readjust surface input depths
               self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)) = self%z_Inp(i, self%nval_deep(i) + 1 :self%nval(i)) - grid%lake_level_old + grid%lake_level

               ! Adjust surface inflow to new lake level
               CALL grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_start(i, :), self%nval_surface(i), self%Qs_start(i, :))
               CALL grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i), self%Qs_end(i, :))

            ENDIF ! ENDIF not lake_level...

            ! Temporal treatment of inflow
            IF ((datum <= self%tb_start(i)) .or. (self%eof(i) == 1)) THEN ! IF datum before first date or END of file reached
               goto 8
            ELSE
               DO WHILE (.not. ((datum >= self%tb_start(i)) .and. (datum <= self%tb_end(i)))) ! DO until datum between dates
                  self%tb_start(i) = self%tb_end(i) ! Move one step in time
                  self%Q_start(i, :) = self%Q_end(i, :)
                  self%Qs_start(i, :) = self%Qs_end(i, :)
                  self%Qs_read_start(i, :) = self%Qs_read_end(i, :)
                  read (fnum(i), *, END=7) self%tb_end(i), (self%Inp_read_end(i, j), j=1, self%nval(i))
                  CALL count_read(self, i)

                  IF (state%has_deep_input(i)) THEN
                  CALL Integrate(self%z_Inp(i, :), self%Inp_read_end(i, :), self%Q_read_end(i, :), self%nval_deep(i))
                  CALL grid%interpolate_to_face_from_second(self%z_Inp(i, :), self%Q_read_end(i, :), self%nval_deep(i), self%Q_end(i, :))
                  ENDIF

                  IF (state%has_surface_input(i)) THEN
                     CALL Integrate(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Inp_read_end(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i))
                     CALL grid%interpolate_to_face_from_second(self%z_Inp(i, self%nval_deep(i) + 1:self%nval(i)), self%Qs_read_end(i, :), self%nval_surface(i), self%Qs_end(i, :))
                  ENDIF
               ENDDO ! ENDDO WHILE
            ENDIF

            ! Linearly interpolate value at correct datum (Q_inp is on face grid)
            DO j = 1, ubnd_fce
               Q_inp(i,j) = (self%Q_start(i,j) + self%Qs_start(i,j)) + (datum-self%tb_start(i))/(self%tb_end(i)-self%tb_start(i))* &
               (self%Q_end(i,j) + self%Qs_end(i,j) - self%Q_start(i,j) - self%Qs_start(i,j))
            ENDDO
            goto 11

            ! IF END of file reached, set to closest available value
7          self%eof(i) = 1
8          Q_inp(i,1:ubnd_fce) = self%Q_start(i,1:ubnd_fce) + self%Qs_start(i,1:ubnd_fce)
            goto 11

            ! IF no data available
9          write(*,*) '[WARNING] ','No data found in ',trim(fname(i)),' file. Check number of depths. Values set to zero.'
            self%eof(i) = 1
            Q_inp(i, 1:ubnd_fce) = 0.0_RK
            self%Q_start(i, 1:ubnd_fce) = 0.0_RK
            self%Qs_start(i, 1:ubnd_fce) = 0.0_RK
11          continue

         ENDDO ! ENDDO i=1,4
         ! Q_vert is the integrated difference between in- and outflow (starting at the lake bottom)
         ! Q_vert is located on the face grid, m^3/s
         Q_vert(1)=0
         Q_vert(2:ubnd_fce) = Q_inp(1, 2:ubnd_fce) + Q_inp(2, 2:ubnd_fce)

         ! The final Q_inp is located on the volume grid
         DO i = 1, 4
            DO j = 1, ubnd_vol
               Q_inp(i, j) = Q_inp(i, j + 1) - Q_inp(i, j)
            ENDDO
            Q_inp(i,ubnd_vol + 1) = 0
         ENDDO

      END associate
   END SUBROUTINE

   SUBROUTINE count_read(self, i)
      IMPLICIT NONE
      class(GenericLateralModule) :: self
      integer i
      self%number_of_lines_read(i) = self%number_of_lines_read(i) + 1
   END SUBROUTINE

   SUBROUTINE detect_first(self, state, first)
      IMPLICIT NONE
      class(GenericLateralModule) :: self
      class(ModelState) :: state
      logical, intent(out) :: first

      first = .not. allocated(self%z_Inp)
      IF (first) THEN
         ! Allocate arrays for very first iteration
         allocate (self%z_Inp(1:4, 1:state%nz_input)) ! Input depths
         allocate (self%Inp_read_start(1:4, 1:state%nz_input)) ! Raw input read
         allocate (self%Inp_read_end(1:4, 1:state%nz_input)) ! Raw input read
         allocate (self%Q_read_start(1:4, 1:state%nz_input)) ! Integrated input
         allocate (self%Q_read_end(1:4, 1:state%nz_input)) ! Integrated input
         allocate (self%Qs_read_start(1:4, 1:state%nz_input))  ! Integrated surface input
         allocate (self%Qs_read_end(1:4, 1:state%nz_input))  ! Integrated surface input
         allocate (self%Q_start(1:4, 1:self%grid%nz_grid+1)) ! Input interpolated on grid
         allocate (self%Q_end(1:4, 1:self%grid%nz_grid+1)) ! Input interpolated on grid
         allocate (self%Qs_start(1:4, 1:self%grid%nz_grid+1)) ! Surface input interpolated on grid
         allocate (self%Qs_end(1:4, 1:self%grid%nz_grid+1)) ! Surface input interpolated on grid
      ENDIF
   END SUBROUTINE

END MODULE 





!<    +---------------------------------------------------------------+
!     |  Stability MODULE
!     |  - CONTAINS methods to update cmue_cn /qe and NN
!<    +---------------------------------------------------------------+
MODULE  strat_stability
   USE strat_kinds, only: RK
   USE MOD_Lake_Const
   USE strat_grid
   USE strat_simdata

   IMPLICIT NONE

!-----------------------------------------------------------
PRIVATE
!-----------------------------------------------------------

   ! Common Types
   type, public :: StabilityModule
      class(StaggeredGrid), pointer :: grid
      class(ModelConfig), pointer :: model_cfg
      class(ModelParam), pointer :: model_param

   CONTAINS
      procedure, pass :: init => stability_module_init
      procedure, pass :: update => stability_module_update
      procedure, pass :: update_cmue_cn => stability_module_update_cmue_cn
      procedure, pass :: update_cmue_qe => stability_module_update_cmue_qe
      procedure, pass :: update_NN => stability_module_update_NN

   END type

!-----------------------------------------------------------
CONTAINS
!-----------------------------------------------------------

   SUBROUTINE stability_module_init(self, grid, model_cfg, model_param)
      IMPLICIT NONE
      class(StabilityModule) :: self
      class(StaggeredGrid), target :: grid
      class(ModelConfig), target :: model_cfg
      class(ModelParam), target :: model_param

      self%grid => grid
      self%model_cfg => model_cfg
      self%model_param => model_param
   END SUBROUTINE


   ! Update state variables
   SUBROUTINE stability_module_update(self, state)
      IMPLICIT NONE
      class(StabilityModule) :: self
      class(ModelState) :: state
      real(RK), DIMENSION(self%grid%ubnd_fce) :: beta

      ! DO buoyancy update (update NN)
      CALL self%update_NN(state%T, state%S, state%rho, state%NN)

      ! Update cmue depending on selected stabilty function
      IF (self%model_cfg%stability_func == 1) THEN
         CALL self%update_cmue_cn(state%cmue1, state%cmue2)

      ELSE IF (self%model_cfg%stability_func == 2) THEN
         beta(1:self%grid%ubnd_fce) = state%NN(1:self%grid%ubnd_fce)*(state%k(1:self%grid%ubnd_fce)/state%eps(1:self%grid%ubnd_fce))**2

         beta(1) = 0
         beta(self%grid%ubnd_fce) = 0

         CALL self%update_cmue_qe(beta, state%cmue1, state%cmue2, state%cde)

      ENDIF

   END SUBROUTINE


   ! Compute NN from T and salinity
   SUBROUTINE stability_module_update_NN(self, T, S, rho, NN)
      IMPLICIT NONE
      class(StabilityModule) :: self

      ! Global variables
      real(RK), DIMENSION(:), intent(in) :: T, S
      real(RK), DIMENSION(:), intent(inout) :: NN, rho

      ! Local variables
      real(RK) :: buoy(self%grid%length_fce)
      real(RK) :: rho0t(self%grid%length_fce), rho0st(self%grid%length_fce)
      integer :: i

      associate (grd=>self%grid)

         DO i = 1, grd%ubnd_fce - 1
               rho0t(i) = 0.9998395_RK + T(i)*(6.7914e-5_RK + T(i)*(-9.0894e-6_RK + T(i)* &
                                             (1.0171e-7_RK + T(i)*(-1.2846e-9_RK + T(i)*(1.1592e-11_RK + T(i)*(-5.0125e-14_RK))))))
               rho0st(i) = (8.181e-4_RK + T(i)*(-3.85e-6_RK + T(i)*(4.96e-8_RK)))*S(i)
               rho(i) = rho_0*(rho0t(i) + rho0st(i))

               buoy(i) = -g*(rho(i) - rho_0)/rho_0
         ENDDO

         NN(2:grd%ubnd_fce - 1) = grd%meanint(1:grd%ubnd_vol - 1)*(buoy(2:grd%ubnd_fce - 1) - buoy(1:grd%ubnd_fce - 2))
         NN(1) = NN(2)
         NN(grd%ubnd_fce) = NN(grd%ubnd_fce - 1)

      END associate
   END SUBROUTINE


   SUBROUTINE stability_module_update_cmue_cn(self, cmue1, cmue2)
      IMPLICIT NONE

      ! Global variables
      class(StabilityModule) :: self
      real(RK), DIMENSION(:), intent(inout) :: cmue1, cmue2

      ! Standard version of k-eps model
      cmue1 = cmue
      cmue2 = cmue/Prndtl

      ! Set boundaries
      cmue1(1) = cmue1(2)
      cmue2(1) = cmue2(2)
      cmue1(self%grid%ubnd_fce) = cmue1(self%grid%ubnd_fce - 1)
      cmue2(self%grid%ubnd_fce) = cmue2(self%grid%ubnd_fce - 1)
   END SUBROUTINE


   SUBROUTINE stability_module_update_cmue_qe(self, beta, cmue1, cmue2, cde)
      IMPLICIT NONE
      class(StabilityModule) :: self
      real(RK), DIMENSION(:), intent(in) :: beta
      real(RK), DIMENSION(:), intent(inout) ::cmue1, cmue2
      real(RK) :: cde
      real(RK) :: gh, sm, sh
      integer :: i

      DO i = 2, self%grid%ubnd_fce - 1
         gh = -cde**2*0.5_RK*beta(i)
         IF (gh > 0.02) gh = gh - (gh - 0.02_RK)**2/(gh + 0.0233_RK - 2*0.02_RK)
         IF (gh < -0.28) gh = -0.28_RK

         sm = 1.0_RK - 3*c1 - 6*a1/b1 - 3*a2*gh*((b2 - 3*a2)*(1.0_RK - 6*a1/b1) - 3*c1*(b2 + 6*a1))
         sm = a1*sm/((1.0_RK - 3*a2*gh*(6*a1 + b2))*(1.0_RK - 9*a1*a2*gh))
         sh = a2*(1.0_RK - 6*a1/b1)/(1.0_RK - 3*a2*gh*(6*a1 + b2))

         cmue1(i) = sqrt(2.0_RK)*cde*sm
         cmue2(i) = sqrt(2.0_RK)*cde*sh
      ENDDO

      ! Set boundaries
      cmue1(1) = cmue1(2)
      cmue2(1) = cmue2(2)
      cmue1(self%grid%ubnd_fce) = cmue1(self%grid%ubnd_fce - 1)
      cmue2(self%grid%ubnd_fce) = cmue2(self%grid%ubnd_fce - 1)
   END SUBROUTINE

END MODULE  strat_stability





!<    +---------------------------------------------------------------+
!     |  ModelVariable implementation for K and EPS variable
!<    +---------------------------------------------------------------+
MODULE  strat_keps
   USE strat_kinds, only: RK
   USE MOD_Lake_Const
   USE strat_simdata
   USE strat_statevar
   USE strat_grid
   USE strat_solver

   IMPLICIT NONE

!-----------------------------------------------------------
PRIVATE
!-----------------------------------------------------------

   ! Class for K variable
   type, extends(ModelVariable), public :: KModelVar
   CONTAINS
      procedure, pass(self), public :: calc_terms => k_var_calc_terms
      procedure, pass(self), public :: post_solve => k_var_post_solve
   END type

   ! Class for Eps variable
   type, extends(ModelVariable), public :: EpsModelVar
   CONTAINS
      procedure, pass(self), public :: calc_terms => eps_var_calc_terms
      procedure, pass(self), public :: post_solve => eps_var_post_solve
   END type

!-----------------------------------------------------------
CONTAINS
!-----------------------------------------------------------

   ! Calculates terms for linear system for a K variable
   SUBROUTINE k_var_calc_terms(self, state, param, sources, boundaries)
      class(KModelVar), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(ModelParam), intent(inout) :: param
      real(RK), DIMENSION(:) ::  sources, boundaries
      real(RK) :: rhs_0, rhs_ubnd
      real(RK) :: pminus(2:self%grid%ubnd_fce - 1), pplus(2:self%grid%ubnd_fce - 1), Prod, Buoy, Diss
      integer :: i
      associate (grid=>self%grid, &
                  ubnd_fce=>self%grid%ubnd_fce, &
                  ubnd_vol=>self%grid%ubnd_vol)

         !!!!!!!! Precalculations !!!!!!!!
         sources = 0.0_RK
         boundaries = 0.0_RK
         state%ko(1:ubnd_fce) = state%k(1:ubnd_fce) ! ko = TKE at old time step

         ! Diffusivity for k is located on volume grid
         state%avh(2:ubnd_vol - 1) = 0.5_RK/sig_k*(state%num(2:ubnd_fce - 2) + state%num(3:ubnd_fce - 1))

         IF (self%cfg%flux_condition == 1 .and. self%cfg%turbulence_model == 1) THEN
            state%avh(1) = 0.0_RK
            state%avh(ubnd_vol) = 0.0_RK
         ELSE
            state%avh(1) = 2*state%u_taub**4/(state%eps(1) + state%eps(2)) ! = 0 for no shear stress
            state%avh(ubnd_vol) = 2*state%u_taus**4/(state%eps(ubnd_fce) + state%eps(ubnd_fce - 1)) ! = 0 for no shear stress
         ENDIF

         DO i = 2, ubnd_fce - 1
            Prod = state%P(i) + state%P_Seiche(i) ! Add seiche energy
            Buoy = state%B(i)
            Diss = state%eps(i)
            IF (Prod + Buoy > 0) THEN
               pplus(i) = Prod + Buoy
               pminus(i) = Diss
            ELSE
               pplus(i) = Prod
               pminus(i) = Diss - Buoy
            ENDIF
         ENDDO

         !!!!!!!! Define sources !!!!!!!!
         sources(2:ubnd_fce - 1) = pplus(2:ubnd_fce - 1)

         !!!!!! Define boundary conditions !!!!
         boundaries(2:ubnd_fce - 1) = pminus(2:ubnd_fce - 1)/state%k(2:ubnd_fce - 1)

         IF (self%cfg%flux_condition == 1 .and. self%cfg%turbulence_model == 1) THEN
            ! K(0) and K(ubnd_fce) are assigned in the post processing function
         ELSE ! no fluxes, unity A-matrix + condition on RHS
            rhs_0 = state%u_taub**2/sqrt(state%cm0*state%cde)
            rhs_ubnd = state%u_taus**2/sqrt(state%cm0*state%cde)

            ! Trick to have rhs(1) = rhs_0 and rhs(ubnd_fce) = rhs_ubnd
            ! in discretization, rhs is calculated as rhs(1) = var(1) + sources(1)*dt
            ! Given the equation below, the following results:
            ! rhs(1) = var(1) + (-var(1)/dt + rhs_0/dt)*dt = var(1) -var(1) + rhs_0 = rhs_0
            sources(1) = -(self%var(1)/state%dt) + (rhs_0/state%dt)
            sources(ubnd_fce) = -(self%var(ubnd_fce)/state%dt) + (rhs_ubnd/state%dt)
         ENDIF
      END associate
   END SUBROUTINE


   ! Updates variable boundaries after solve has been executed
   SUBROUTINE k_var_post_solve(self, state)
      class(KModelVar), intent(inout) :: self
      class(ModelState), intent(inout) :: state

      integer :: i
      associate (grid=>self%grid, &
                  ubnd_fce=>self%grid%ubnd_fce, &
                  ubnd_vol=>self%grid%ubnd_vol)

         IF (self%cfg%flux_condition == 1 .and. self%cfg%turbulence_model == 1) THEN
            ! Define TKE at boundary (no flux)
            self%var(1) = self%var(2)
            self%var(self%grid%ubnd_fce) = self%var(self%grid%ubnd_fce - 1)
         ENDIF

         ! check lower limit of k
         DO i = 1, ubnd_fce
            IF (self%var(i) < k_min) self%var(i) = k_min ! Lower limit of TKE
         ENDDO

      END associate
   END SUBROUTINE


   ! Calculates terms for linear system for a EPS variable
   SUBROUTINE eps_var_calc_terms(self, state, param, sources, boundaries)
      class(EpsModelVar), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(ModelParam), intent(inout) :: param
      real(RK), DIMENSION(:) ::  sources, boundaries
      real(RK) :: cee3, rhs_0, rhs_ubnd
      real(RK) :: flux(1:self%grid%ubnd_fce)
      real(RK) :: pminus(2:self%grid%ubnd_fce - 1), pplus(2:self%grid%ubnd_fce - 1), Prod, Buoy, Diss

      integer :: i
      associate (grid=>self%grid, &
                  ubnd_fce=>self%grid%ubnd_fce, &
                  ubnd_vol=>self%grid%ubnd_vol)

         !!!!!!!! Precalculations !!!!!!!!
         sources = 0.0_RK
         boundaries = 0.0_RK

         ! Diffusivity for eps is located on volume grid
         state%avh(1:ubnd_vol) = 0.5_RK/sig_e*(state%num(1:ubnd_fce - 1) + state%num(2:ubnd_fce))

         IF (self%cfg%flux_condition == 1 .and. self%cfg%turbulence_model == 1) THEN
            flux(1) = state%avh(1)*(state%cde*((state%ko(2))**1.5_RK))/(kappa*(K_s + 0.5_RK*grid%h(1)))**2
            flux(ubnd_fce) = state%avh(ubnd_vol)*(state%cde*((state%ko(ubnd_fce - 1))**1.5_RK))/(kappa*(z0 + 0.5_RK*grid%h(ubnd_vol)))**2

            DO i = 2, ubnd_fce - 1
               flux(i) = state%num(i)/sig_e*(state%cde*((state%ko(i))**1.5_RK)/(kappa*(z0 + 0.25_RK*(grid%h(i - 1) + grid%h(i))))**2)
            ENDDO

            state%avh(1) = 0.0_RK
            state%avh(ubnd_vol) = 0.0_RK
         ELSE
            state%avh(1) = 2*state%u_taub**4/sig_e/(state%eps(1) + state%eps(2)) ! = 0 for no shear stress
            state%avh(ubnd_vol) = 2*state%u_taus**4/sig_e/(state%eps(ubnd_fce) + state%eps(ubnd_fce - 1)) ! = 0 for no shear stress
         ENDIF

         DO i = 2, ubnd_fce - 1
            IF (state%B(i) > 0) THEN
               cee3 = 1.
            ELSE
               cee3 = ce3
            ENDIF
            Prod = ce1*state%eps(i)/state%ko(i)*(state%P(i) + state%P_Seiche(i)) ! New code plus seiche
            Buoy = cee3*state%eps(i)/state%ko(i)*state%B(i)
            Diss = ce2*state%eps(i)*state%eps(i)/state%ko(i)

            IF (Prod + Buoy > 0) THEN
               pplus(i) = Prod + Buoy
               pminus(i) = Diss
            ELSE
               pplus(i) = Prod
               pminus(i) = Diss - Buoy
            ENDIF
         ENDDO

         !!!!!!!! Define sources !!!!!!!!
         sources(2:ubnd_fce - 1) = pplus(2:ubnd_fce - 1)

         !!!!!! Define boundary conditions !!!!
         boundaries(2:ubnd_fce - 1) = pminus(2:ubnd_fce - 1)/state%eps(2:ubnd_fce - 1)
         IF (self%cfg%flux_condition == 1 .and. self%cfg%turbulence_model == 1) THEN
            sources(2:ubnd_fce - 1) = sources(2:ubnd_fce - 1) + flux(2:ubnd_fce - 1)*grid%AreaFactor_eps(1:ubnd_fce - 2) ! AreaFactor_eps = 1/A * dA/dz (at epsilon posions)
            IF (grid%Az(1) /= 0) THEN ! Flux from bottom only!
               sources(2) = sources(2) + flux(1)*(grid%Az(1) + grid%Az(2))/(grid%Az(2)*(grid%h(1) + grid%h(2)))
            ENDIF
            sources(ubnd_fce-1)= sources(ubnd_fce-1)+flux(ubnd_fce)*(grid%Az(ubnd_fce)+grid%Az(ubnd_fce-1))/(grid%Az(ubnd_fce-1)*(grid%h(ubnd_vol)+grid%h(ubnd_vol-1)))
            ! eps(0) and eps(ubnd_fce) are assigned in the post processing function
         ELSE ! no fluxes, unity A-matrix + condition on RHS
            rhs_0 = state%u_taub**2/sqrt(state%cm0*state%cde)
            rhs_ubnd = state%u_taus**2/sqrt(state%cm0*state%cde)

            ! Trick to have rhs(1) = rhs_0 and rhs(ubnd_fce) = rhs_ubnd
            ! in discretization, rhs is calculated as rhs(1) = var(1) + sources(1)*dt
            ! Given the equation below, the following results:
            ! rhs(1) = var(1) + (-var(1)/dt + rhs_0/dt)*dt = var(1) -var(1) + rhs_0 = rhs_0
            sources(1) = -(self%var(1)/state%dt) + (rhs_0/state%dt)
            sources(ubnd_fce) = -(self%var(ubnd_fce)/state%dt) + (rhs_ubnd/state%dt)
         ENDIF
      END associate
   END SUBROUTINE


   SUBROUTINE eps_var_post_solve(self, state)
      class(EpsModelVar), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      integer :: i
      real(RK) ::epslim
      associate (grid=>self%grid, &
                  ubnd_fce=>self%grid%ubnd_fce, &
                  ubnd_vol=>self%grid%ubnd_vol)

         IF (self%cfg%flux_condition == 1 .and. self%cfg%turbulence_model == 1) THEN
               ! Define eps at boundary (no flux)
               self%var(1) = self%var(2) + (state%cde*((state%ko(2))**1.5_RK)/(kappa*(K_s + grid%h(1)))**2)*grid%h(1)
               self%var(ubnd_fce)=   self%var(ubnd_fce-1)+(state%cde*((state%ko(ubnd_fce-1))**1.5_RK)/(kappa*(z0 +grid%h(ubnd_vol)))**2)*grid%h(ubnd_vol)
         ENDIF

         ! check lower limit of eps / update num
         DO i = 1, ubnd_fce
            IF (state%NN(i) > 0) THEN
               epslim = 0.212_RK*state%k(i)*sqrt(state%NN(i))
            ELSE
               epslim = eps_min
            ENDIF
            IF (state%eps(i) < epslim) THEN
               state%eps(i) = epslim
            ENDIF
            IF (state%eps(i) < 0) THEN
               write(*, *) 'Dissipation negative'
            ENDIF

            state%num(i) = state%cmue1(i)*state%k(i)*state%k(i)/state%eps(i) + 1.5e-6_RK
            state%nuh(i) = state%cmue2(i)*state%k(i)*state%k(i)/state%eps(i) + 1.5e-7_RK

         ENDDO

         state%num(1) = kappa*state%u_taub*K_s + avh_min
         state%num(ubnd_fce) = kappa*state%u_taus*z0 + avh_min
      END associate
   END SUBROUTINE

END MODULE  strat_keps





!<    +---------------------------------------------------------------+
!     | Implementation of a statevar for a generic transport variable
!     | At the moment used for S, but could be used for any other biogeochemical type
!<    +---------------------------------------------------------------+
MODULE  strat_transport
   USE strat_kinds, only: RK
   USE MOD_Lake_Const
   USE strat_simdata
   USE strat_statevar
   USE strat_grid
   USE strat_solver

   IMPLICIT NONE

!-----------------------------------------------------------
PRIVATE
!-----------------------------------------------------------

    type, extends(ModelVariable), public :: TransportModVar
        real(RK), DIMENSION(:), pointer :: dVar
    CONTAINS
        procedure, pass(self), public :: assign_external_source => transport_assign_external_source
        procedure, pass(self), public :: calc_terms => transport_var_calc_terms
    END type

!-----------------------------------------------------------
CONTAINS
!-----------------------------------------------------------

   ! Method to assign external source
   SUBROUTINE transport_assign_external_source(self, dVar)
      class(TransportModVar), intent(inout) :: self
      real(RK), DIMENSION(:), target :: dVar
      self%dVar => dVar
   END SUBROUTINE

   ! Calculate source terms according to external source
   SUBROUTINE transport_var_calc_terms(self, state, param, sources, boundaries)
      class(TransportModVar), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(ModelParam), intent(inout) :: param
      real(RK), DIMENSION(:) ::  sources, boundaries
      associate (grid=>self%grid, &
                  ubnd_fce=>self%grid%ubnd_fce, &
                  ubnd_vol=>self%grid%ubnd_vol)

         !!!!!!!! Define sources !!!!!!!!
         sources = self%dVar ! We only have a generic, external source!

         ! No explicit boundary conditions
         boundaries(1:ubnd_vol) = 0

      END associate
   END SUBROUTINE

END MODULE  strat_transport





!<    +---------------------------------------------------------------+
!     | Implementation of a statevar for U/V variable
!<    +---------------------------------------------------------------+
MODULE  strat_windshear
   USE strat_kinds, only: RK
   USE MOD_Lake_Const
   USE strat_simdata
   USE strat_statevar
   USE strat_grid
   USE strat_solver

   IMPLICIT NONE

!-----------------------------------------------------------
PRIVATE
!-----------------------------------------------------------

   type, extends(ModelVariable), public :: UVModelVar
      real(RK), pointer :: stress_t
   CONTAINS
      procedure, pass(self), public :: calc_terms => uv_var_calc_terms
      procedure, pass(self), public :: assign_shear_stress => uv_var_assign_shear_stress
   END type

!-----------------------------------------------------------
CONTAINS
!-----------------------------------------------------------

   ! Method to assign shear stress variable acting on this state variable
   SUBROUTINE uv_var_assign_shear_stress(self, stress_t)
      class(UVModelVar), intent(inout) :: self
      real(RK), target ::stress_t
      self%stress_t => stress_t
   END SUBROUTINE


   ! Calculate source terms for u or v
   SUBROUTINE uv_var_calc_terms(self, state, param, sources, boundaries)
      class(UVModelVar), intent(inout) :: self
      class(ModelState), intent(inout) :: state
      class(ModelParam), intent(inout) :: param
      real(RK) :: integr, length
      real(RK), DIMENSION(:) ::  sources, boundaries
      real(RK), DIMENSION(size(self%var)) ::  uv_norm
      associate (grid=>self%grid, &
                  ubnd_fce=>self%grid%ubnd_fce, &
                  ubnd_vol=>self%grid%ubnd_vol)

         !!!!!!!! Precalculations !!!!!!!!
         IF (self%cfg%pressure_gradients == 1) THEN
               integr = sum(self%var(1:ubnd_vol))
               length = sqrt(grid%Az(ubnd_fce))
         ELSE
               integr = 0
               length = 0
         ENDIF
         
         uv_norm = sqrt(state%U**2 + state%V**2)
         
         !!!!!!!! Define sources !!!!!!!!
         sources = 0
         IF (self%cfg%pressure_gradients == 1) THEN !Svensson 1978
               sources(2:ubnd_vol - 1) = -pi**2*rho_0*g*integr/ubnd_vol*grid%max_depth/length**2
         elseif (self%cfg%pressure_gradients == 2) THEN !???
               sources(2:ubnd_vol - 1) = -state%drag*self%var(2:ubnd_vol - 1)* &
                                       uv_norm(2:ubnd_vol - 1)* &
                                       grid%dAz(2:ubnd_vol - 1)/grid%Az(2:ubnd_vol - 1)
         ENDIF
         
         ! Set surface condition based on shear stress variable
         sources(ubnd_vol) = self%stress_t/grid%h(ubnd_vol)

         !!!!!!!! Explicit boundary conditions !!!!!
         boundaries(1) = state%drag*uv_norm(1)/grid%h(1)
         boundaries(2:ubnd_vol) = 0

      END associate
   END SUBROUTINE

END MODULE  strat_windshear





!<    +---------------------------------------------------------------+
!     |  Forcing MODULE
!     |  - Adjusted heat fluxes to snow- and ice-covered lake
!     |  - Reads forcing file and updates state
!     |  - Updates Coriolis force terms
!<    +---------------------------------------------------------------+
MODULE  strat_forcing
   USE strat_kinds, only: RK
   USE strat_simdata
   USE MOD_Lake_Const
   USE strat_grid
   USE strat_utilities
   IMPLICIT NONE

!-----------------------------------------------------------
PRIVATE
!-----------------------------------------------------------

   type, public :: ForcingModule
      integer :: nz_grid_max
      class(ModelConfig), pointer :: cfg
      class(StaggeredGrid), pointer :: grid
      class(ModelParam), pointer :: param
      character(len=:), allocatable  :: file  ! Forcing file name
   CONTAINS
      procedure, pass :: init => forcing_init
      procedure, pass :: update => Simstrat_SurfFlux
      procedure, pass :: update_coriolis => forcing_update_coriolis

   END type

!-----------------------------------------------------------
CONTAINS
!-----------------------------------------------------------

   SUBROUTINE forcing_init(self, model_config, model_param, forcing_file, grid)
      IMPLICIT NONE
      class(ForcingModule) :: self
      class(StaggeredGrid), target :: grid
      class(ModelConfig), target :: model_config
      class(ModelParam), target :: model_param
      character(len=:), allocatable :: forcing_file

      self%cfg => model_config
      self%param => model_param
      self%grid => grid
      self%file = forcing_file

      IF(self%cfg%snow_model == 1) CALL warn('The snow MODULE  is turned on. This MODULE  needs precipitation data, note that the last column in the forcing file will be interpreted as precipitation.')
      !IF precipitation column is missing in forcing file, Simstrat will read dates as precipitation
   END SUBROUTINE


    ! Compute appropriate forcing parameters at given datum
    ! Copied from old simstrat, NEEDS refactoring!
    ! AG 2014: revision
    !######################################################################
   SUBROUTINE Simstrat_SurfFlux(self, state,&
               ! "in" arguments
               ! ---------------------------
               scwat     , snlay     , betaopt   , fetchopt  ,&
               zopt      , rlat      , dplak     , dtlak     ,&
               hwref     , htref     , hqref     , usurf     ,&
               vsurf     , tmref     , qmref     , arhos     ,&
               psurf     , sols      , soll      , solsd     ,&
               solld     , lwdns     , hpbl      , sabg      ,&
               stke1     , dzlaktop  , zlakebot  , lktmptop  ,&
               icefrtop  , dzssbtop  , xwicetop  , xwliqtop  ,&
               crain     , csnow     , lrain     , lsnow     ,&
               tssbtop   , zcbcv     , ipatch    , &
               ! "inout" arguments
               ! -------------------
               tskin     , z0m       , z0h       , z0q       ,&
               felak     , btpri     ,&
               ! "out" arguments
               ! ---------------------------
               fsena     , fevpa     , lfevpa    , fseng     ,&
               fevpg     , olrg      , fgrnd     , trad      ,&
               taux      , tauy      , ustar     , qstar     ,&
               tstar     , emis      , zol       , &
               rib       , ram       , rah       , raw       ,&
               wdm       , t2m       , q2m       , u10m      ,&
               v10m      , fm10m     , fq10m     , fh2m      ,&
               fq2m      , fm        , fq        , fh        ,&
               shfdt     , urban_call)
        
!================================================================================
   USE MOD_Lake_CoLML, only : CoLML_SurfFlux
   USE MOD_Lake_Subs, only : LakStaParms
   USE MOD_Lake_Const, only : tfrz, hvap, hsub
   USE MOD_Namelist, only : DEF_USE_COLML_FLUX_SCHEME
   IMPLICIT NONE
   class(ForcingModule) :: self
   class(ModelState) :: state
!==============================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in)      :: &
      ipatch                  ,&! patch index
      scwat                   ,&! land patch type (4=deep lake, 5=shallow lake)
      snlay                   ,&! number of snow layers
      betaopt                 ,&! option for betaprime calculation, 1: constant, 2: equation
      fetchopt                ,&! option for fetch length, 1: constant, 2: equation
      zopt                      ! option for roughness length, 1: constant, 2: Subin et al. (2012)

   real(RK), intent(in)     :: &
      rlat                    ,&! latitude in radians
      dtlak                   ,&! time step length
      dplak                   ,&! column lake depth (m)
      hwref                   ,&! observational height of wind [m]
      htref                   ,&! observational height of temperature [m]
      hqref                   ,&! observational height of humidity [m]
      usurf                   ,&! wind component in eastward direction [m/s]
      vsurf                   ,&! wind component in northward direction [m/s]
      tmref                   ,&! temperature at the reference height [K]
      qmref                   ,&! specific humidity at the reference height [kg/kg]
      arhos                   ,&! air density [kg/m3]
      psurf                   ,&! surface pressure [Pa]
      sols                    ,&! atm vis direct beam solar rad onto srf [W m-2]
      soll                    ,&! atm nir direct beam solar rad onto srf [W m-2]
      solsd                   ,&! atm vis diffuse solar rad onto srf [W m-2]
      solld                   ,&! atm nir diffuse solar rad onto srf [W m-2]
      lwdns                   ,&! downward longwave radiation at surface [w/m2]
      hpbl                    ,&! atmospheric boundary layer height [m]
      zcbcv                   ,&! convective boundary height [m]
      csnow                   ,&! convective snowfall [kg/(m2 s)]
      lsnow                   ,&! large scale snowfall [kg/(m2 s)]
      crain                   ,&! convective rainfall [kg/(m2 s)]
      lrain                   ,&! large scale rainfall [kg/(m2 s)]
      sabg                    ,&! solar radiation absorbed by ground [W/m2]
      stke1                   ,&! top level eddy conductivity [W/m/K]
      dzlaktop                ,&! thickness of top lake layer [m]
      zlakebot                ,&! depth of top lake layer [m]
      lktmptop                ,&! temperature of top lake layer [K]
      icefrtop                ,&! lake mass fraction of top lake layer that is frozen
      dzssbtop                ,&! thickness of top lake layer [m]
      xwicetop                ,&! ice mass in top lake layer [kg/m2]
      xwliqtop                ,&! liquid water mass in top lake layer [kg/m2]
      tssbtop                   ! temperature of top lake layer [K]
        
!  ------------------------- inout variables ---------------------------
   real(RK), intent(inout)  :: &
      tskin                   ,&! surface temperature (kelvin)
      z0m                     ,&! roughness length over ground, momentum [m]
      z0h                     ,&! roughness length over ground, sensible heat [m]
      z0q                     ,&! roughness length over ground, latent heat [m]
      felak                   ,&! fetch length of lake [m]
      btpri                     ! fraction of solar radiation in the NIR, only used when COLMFLUX is defined

!  ------------------------- output variables ---------------------------
   real(RK), intent(out)    :: &
      emis                    ,&! averaged bulk surface emissivity
      wdm                     ,&! lowest level wind speed including the stablity effect [m/s]
      ram                     ,&! aerodynamical resistance [s/m]
      rah                     ,&! thermal resistance [s/m]
      raw                     ,&! moisture resistance [s/m]
      shfdt                   ,&! derivative of srf sen+lat  heat flux wrt srf temp [W/m2/K]
      taux                    ,&! wind stress: E-W [kg/m/s**2]
      tauy                    ,&! wind stress: N-S [kg/m/s**2]
      t2m                     ,&! 2 m height air temperature [kelvin]
      q2m                     ,&! 2 m height specific humidity [kg/kg]
      u10m                    ,&! wind speed at 10m [m/s]
      v10m                    ,&! wind speed at 10m [m/s]
      zol                     ,&! dimensionless height (z/L) used in Monin-Obukhov theory
      rib                     ,&! bulk Richardson number in surface layer
      fm                      ,&! integral of profile function for momentum
      fh                      ,&! integral of profile function for heat
      fq                      ,&! integral of profile function for moisture
      fh2m                    ,&! relation for temperature at 2m
      fq2m                    ,&! relation for specific humidity at 2m
      fm10m                   ,&! integral of profile function for momentum at 10m
      fq10m                   ,&! integral of profile function for moisture at 10m
      ustar                   ,&! u* in similarity theory [m/s]
      qstar                   ,&! humidity scaling parameter [kg/kg]
      tstar                   ,&! temperature scaling parameter [kelvin]
      fseng                   ,&! sensible heat flux into the ground [W/m2]
      fevpg                   ,&! latent heat flux into the ground [W/m2]
      fsena                   ,&! sensible heat from canopy height to atmosphere [W/m2]
      fevpa                   ,&! evapotranspiration from canopy height to atmosphere [mm/s]
      lfevpa                  ,&! latent heat flux from canopy height to atmosphere [W/m2]
      olrg                    ,&! outgoing longwave radiation [W/m2]
      fgrnd                   ,&! ground heat flux [W/m2]
      trad                      ! radiative temperature [K]

   !+WMEJ urban_call is optional argument for SUBROUTINE laketem, it is not used in this version
   logical, optional, intent(in) :: urban_call   ! whether it is a urban CALL

!  ------------------------- local variables ---------------------------
   real(r8)                 :: &
      fgrnd1                  ,&! ground heat flux [W/m2]
      htvp                    ,&! heat transfer coefficient for vapor [W/m2/K]
      Clouds                  ,&! Cloud cover
      u2m                     ,&! 2 m height wind speed [m/s]
      fu                      ,&! ---
      Vap_wat                 ,&!  Water vapor saturation pressure in air at water temperature [milibar]
      Vap_atm                 ,&! water vapor pressure in the atmosphere [milibar]
      heat0                   ,&! net heat flux at the surface [W/m2] form forcing ??
      emissivity              ,&! emissivity of the surface
      Cloud                   ,&! Cloud cover
      H_A                     ,&! longwave radiation [W/m2]
      H_W                     ,&! outgoing longwave radiation [W/m2]
      H_K                     ,&! sensible heat flux [W/m2]
      H_V                     ,&! latent heat flux [W/m2]
      H_tot                   ,&! total heat flux [W/m2]
      T_surf                  ,&! surface temperature [K]
      F_glob                  ,&! solar radiation absorbed by ground [W/m2]
      F_snow                  ,&! solar radiation absorbed by snow [W/m2]
      F_snowice               ,&! solar radiation absorbed by snow-ice [W/m2]
      F_ice                   ,&! solar radiation absorbed by ice [W/m2]            
      rhosnow                   ! snow density [kg/m3]

   integer :: nval_offset
   real(RK) :: tau
   real(RK) :: A_s(8), A_e(8), A_cur(8) ! 8 is the maximum number of rows in the forcing file
   real(RK) :: qa, q0, sh_c
   save A_s, A_e

   integer :: lb
!================================================================================

      IF (.NOT.  DEF_USE_COLML_FLUX_SCHEME) THEN

         IF((state%snow_h + state%white_ice_h + state%black_ice_h) > 0) THEN ! Ice cover
            T_surf = self%param%Freez_Temp
         ELSE
            T_surf = state%T(self%grid%ubnd_vol) ! Free water
         ENDIF

         ! Number of values to read, depending on filtered wind and precipitation
         IF (self%cfg%use_filtered_wind .and. self%cfg%ice_model == 0) THEN
            nval_offset = 1
         ELSE IF (self%cfg%ice_model == 1 .and. self%cfg%use_filtered_wind) THEN
            nval_offset = 2
         ELSE IF (self%cfg%ice_model == 1) THEN
            nval_offset = 1
         ELSE
            nval_offset = 0
         ENDIF

         IF (self%cfg%forcing_mode == 3 .and. self%cfg%ice_model == 4) THEN
            CALL error('Forcing mode 3 requires additional Cloud as input variable, mode 4 requires Hnet input, please modify this [Simstrat_SurfFlux].')
         ENDIF

         !*** Forcing mode 1 (date, U10, V10, T_lake, H_sol)
         IF (self%cfg%forcing_mode == 1) THEN
            IF (self%cfg%ice_model == 1) THEN
               CALL error('Ice MODULE  not compatible with forcing mode 1, USE 2, 3 or 5.')
            ENDIF
            !===========================================================
            !-WMEJ Modified the way to read data
            A_cur(1) = usurf ! U10
            A_cur(2) = vsurf ! U10
            A_cur(3) = tskin - tfrz ! Lake surface temperature, from forcing ??
            A_cur(4) = sabg  ! Solar radiation absorbed by ground
            !-----------------------------------------------------------
            !-WMEJ Original way to read data
            ! CALL self%read (state%datum, A_s, A_e, A_cur, 4 + nval_offset, state%first_timestep) !-WMEJ No longer needed
            ! CALL self%read (state%datum, A_s, A_e, A_cur, 4 + nval_offset, state%first_timestep) !-WMEJ No longer needed
            !===========================================================
            
            state%u10 = A_cur(1)*self%param%f_wind ! MS 2014: added f_wind
            state%v10 = A_cur(2)*self%param%f_wind ! MS 2014: added f_wind
            state%uv10 = sqrt(state%u10**2 + state%v10**2) ! AG 2014
            state%SST = A_cur(3) ! Lake surface temperature
            ! state%rad0 = max(A_cur(4)*(1 - state%albedo_water)*(1 - self%param%beta_sol),0.0_RK) ! MS: added beta_sol and albedo_water
            state%rad0 = max(A_cur(4)*(1 - self%param%beta_sol),0.0_RK) ! MS: added beta_sol and albedo_water  !-WMEJ: The reflection of solar radiation must be completed in the lake driver
            state%heat = 0.0_RK
            state%T_atm = 0.0_RK
            state%precip = 0.0_RK
            IF (self%cfg%use_filtered_wind) THEN
               state%Wf = A_cur(5) ! AG 2014
            ENDIF

         !*** Forcing mode 2 (date, U10, V10, T_atm, H_sol, Vap)
         ELSE IF (self%cfg%forcing_mode >= 2) THEN
            IF (self%cfg%forcing_mode == 2) THEN
               !===========================================================
               !-WMEJ Modified the way to read data
               A_cur(1) = usurf ! U10
               A_cur(2) = vsurf ! U10
               A_cur(3) = tmref  - tfrz ! Lake surface temperature, from forcing ??
               A_cur(4) = sabg  ! Solar radiation absorbed by ground
               A_cur(5) = qmref * psurf / (0.622 + 0.378 * qmref) ! Vapour pressure
               A_cur(6) = crain + csnow + lrain + lsnow ! Precipitation [m/h]      !-WMEJ: use_filtered_wind = .false.
               !-----------------------------------------------------------
               !-WMEJ Original way to read data
               ! CALL self%read (state%datum, A_s, A_e, A_cur, 5 + nval_offset, state%first_timestep)   !-WMEJ No longer needed
               !===========================================================
               
               state%u10 = A_cur(1)*self%param%f_wind ! MS 2014: added f_wind
               state%v10 = A_cur(2)*self%param%f_wind ! MS 2014: added f_wind
               state%T_atm = A_cur(3)
               A_cur(4) = max(A_cur(4),0.0_RK)  ! To avoid negative values because of numerical problems

               IF (state%black_ice_h > 0 .and. state%white_ice_h == 0 .and. state%snow_h == 0) THEN !-- Ice
                  ! F_glob = A_cur(4)*(1 - ice_albedo) * self%param%p_albedo               !-WMEJ The reflection of solar radiation must be completed in the lake driver
                  F_glob = A_cur(4) * self%param%p_albedo
               ELSE IF (state%white_ice_h > 0 .and. state%snow_h == 0) THEN !-- Snowice
                  ! F_glob = A_cur(4)*(1 - snowice_albedo) * self%param%p_albedo           !-WMEJ The reflection of solar radiation must be completed in the lake driver
                  F_glob = A_cur(4) * self%param%p_albedo
               ELSE IF (state%snow_h > 0) THEN !-- Snow
                  ! F_glob = A_cur(4)*(1 - snow_albedo) * self%param%p_albedo              !-WMEJ The reflection of solar radiation must be completed in the lake driver
                  F_glob = A_cur(4)* self%param%p_albedo
               ELSE !-- Water
                  ! F_glob = A_cur(4)*(1 - state%albedo_water) * self%param%p_sw           !-WMEJ The reflection of solar radiation must be completed in the lake driver
                  F_glob = A_cur(4) * self%param%p_sw
               ENDIF

               Vap_atm = A_cur(5)
               Cloud = 0.5
               IF (self%cfg%use_filtered_wind) state%Wf = A_cur(6) ! AG 2014
               IF (self%cfg%snow_model == 1 .and. self%cfg%use_filtered_wind) THEN
                  state%precip = max(A_cur(7),0.0_RK)
               ELSE IF (self%cfg%snow_model == 1) THEN
                  state%precip = max(A_cur(6),0.0_RK)
               ENDIF

            !*** Forcing mode 3 (date, U10, V10, T_atm, H_sol, Vap, Clouds)
            ELSE IF (self%cfg%forcing_mode == 3) THEN
               !===========================================================
               !-WMEJ Modified the way to read data
               A_cur(1) = usurf ! U10
               A_cur(2) = vsurf ! U10
               A_cur(3) = tmref  - tfrz ! Lake surface temperature, from forcing ??
               A_cur(4) = sabg  ! Solar radiation absorbed by ground
               A_cur(5) = qmref * psurf / (0.622 + 0.378 * qmref) ! Vapour pressure
               A_cur(6) = Clouds ! Cloudiness [0-1]
               A_cur(7) = crain + csnow + lrain + lsnow ! Precipitation [m/h]           !-WMEJ: use_filtered_wind = .false.
               !-----------------------------------------------------------
               !-WMEJ Original way to read data
               ! CALL self%read (state%datum, A_s, A_e, A_cur, 6 + nval_offset, state%first_timestep)  !-WMEJ No longer needed
               !===========================================================

               state%u10 = A_cur(1)*self%param%f_wind ! MS 2014: added f_wind
               state%v10 = A_cur(2)*self%param%f_wind ! MS 2014: added f_wind
               state%T_atm = A_cur(3)
               A_cur(4) = max(A_cur(4),0.0_RK)  ! To avoid negative values because of numerical problems

               IF (state%black_ice_h > 0 .and. state%white_ice_h == 0 .and. state%snow_h == 0) THEN ! Ice
                  ! F_glob = A_cur(4)*(1 - ice_albedo) * self%param%p_albedo                !-WMEJ The reflection of solar radiation must be completed in the lake driver
                  F_glob = A_cur(4) * self%param%p_albedo
               ELSE IF (state%white_ice_h > 0 .and. state%snow_h == 0) THEN ! Snowice
                  ! F_glob = A_cur(4)*(1 - snowice_albedo) * self%param%p_albedo            !-WMEJ The reflection of solar radiation must be completed in the lake driver
                  F_glob = A_cur(4) * self%param%p_albedo
               ELSE IF (state%snow_h > 0) THEN ! Snow
                  ! F_glob = A_cur(4)*(1 - snow_albedo) * self%param%p_albedo               !-WMEJ The reflection of solar radiation must be completed in the lake driver
                  F_glob = A_cur(4) * self%param%p_albedo
               ELSE ! Water
                  ! F_glob = A_cur(4)*(1 - state%albedo_water) * self%param%p_sw            !-WMEJ The reflection of solar radiation must be completed in the lake driver
                  F_glob = A_cur(4) * self%param%p_sw
               ENDIF

               Vap_atm = A_cur(5)
               Cloud = A_cur(6)
               IF (Cloud < -1e-6 .or. Cloud > 1.000001) THEN
                  write (*,'(A,F12.6)') 'Cloud : ' , Cloud
                  write (*,'(A,F12.6)') 'Date  : ' , state%datum
                  CALL error('Cloudiness should always be between 0 and 1.')
               ENDIF
               IF (self%cfg%use_filtered_wind) state%Wf = A_cur(7) !AG 2014
               IF (self%cfg%snow_model == 1 .and. self%cfg%use_filtered_wind) THEN
                  state%precip = max(A_cur(8),0.0_RK)
               ELSE IF (self%cfg%snow_model == 1) THEN
                  state%precip = max(A_cur(7),0.0_RK)
               ENDIF

            !*** Forcing mode 4 (date, U10, V10, H_net, H_sol)
            ELSE IF (self%cfg%forcing_mode == 4) THEN
               IF (self%cfg%ice_model == 1) THEN
                  CALL error('Ice MODULE  not compatible with forcing mode 4, USE 2, 3 or 5.')
               ENDIF
               !===========================================================
               !-WMEJ Modified the way to read data
               A_cur(1) = usurf ! U10
               A_cur(2) = vsurf ! U10
               A_cur(3) = heat0 ! heat flux from forcing ??
               A_cur(4) = sabg  ! Solar radiation absorbed by ground
               !-----------------------------------------------------------
               !-WMEJ Original way to read data
               ! CALL self%read (state%datum, A_s, A_e, A_cur, 4 + nval_offset, state%first_timestep)  !-WMEJ No longer needed
               !===========================================================
               
               state%u10 = A_cur(1)*self%param%f_wind ! MS 2014: added f_wind
               state%v10 = A_cur(2)*self%param%f_wind ! MS 2014: added f_wind
               heat0 = A_cur(3) ! MS 2014
               ! F_glob = max(A_cur(4)*(1 - state%albedo_water) * self%param%p_sw,0.0_RK)      !-WMEJ The reflection of solar radiation must be completed in the lake driver
               F_glob = max(A_cur(4) * self%param%p_sw,0.0_RK) 
               state%T_atm = 0.0_RK
               IF (self%cfg%use_filtered_wind) state%Wf = A_cur(5) ! AG 2014

            !*** Forcing 5 (date, U10, V10, T_atm, H_sol, Vap, ILWR)
            ELSE IF (self%cfg%forcing_mode == 5) THEN
               !===========================================================
               !-WMEJ Modified the way to read data
               A_cur(1) = usurf ! U10
               A_cur(2) = vsurf ! U10
               A_cur(3) = tmref  - tfrz ! Lake surface temperature, from forcing ??
               A_cur(4) = sabg  ! Solar radiation absorbed by ground
               A_cur(5) = qmref * psurf/100.0_RK / (0.622 + 0.378 * qmref) ! Vapour pressure
               A_cur(6) = lwdns ! downward longwave radiation at surface [w/m2]
               A_cur(7) = crain + csnow + lrain + lsnow ! Precipitation [m/h]      !-WMEJ: use_filtered_wind = .false.
               !-----------------------------------------------------------
               !-WMEJ Original way to read data
               ! CALL self%read (state%datum, A_s, A_e, A_cur, 6 + nval_offset, state%first_timestep)  !-WMEJ No longer needed
               !===========================================================
               
               state%u10 = A_cur(1)*self%param%f_wind !MS 2014: added f_wind
               state%v10 = A_cur(2)*self%param%f_wind !MS 2014: added f_wind
               state%T_atm = A_cur(3)
               A_cur(4) = max(A_cur(4),0.0_RK)  ! To avoid negative values because of numerical problems

               IF (state%black_ice_h > 0 .and. state%white_ice_h == 0 .and. state%snow_h == 0) THEN ! Ice
                  ! F_glob = A_cur(4)*(1 - ice_albedo) * self%param%p_albedo               !-WMEJ The reflection of solar radiation must be completed in the lake driver
                  F_glob = A_cur(4) * self%param%p_albedo
               ELSE IF (state%white_ice_h > 0 .and. state%snow_h == 0) THEN ! Snowice
                  ! F_glob = A_cur(4)*(1 - snowice_albedo) * self%param%p_albedo           !-WMEJ The reflection of solar radiation must be completed in the lake driver
                  F_glob = A_cur(4) * self%param%p_albedo
               ELSE IF (state%snow_h > 0) THEN ! Snow
                  ! F_glob = A_cur(4)*(1 - snow_albedo) * self%param%p_albedo              !-WMEJ The reflection of solar radiation must be completed in the lake driver
                  F_glob = A_cur(4) * self%param%p_albedo
               ELSE ! Water
                  ! F_glob = A_cur(4)*(1 - state%albedo_water) * self%param%p_sw           !-WMEJ The reflection of solar radiation must be completed in the lake driver
                  F_glob = A_cur(4) * self%param%p_sw
               ENDIF

               Vap_atm = A_cur(5)
               H_A = A_cur(6)
               IF (self%cfg%use_filtered_wind) state%Wf = A_cur(7) ! AG 2014
               IF (self%cfg%snow_model == 1 .and. self%cfg%use_filtered_wind) THEN
                  state%precip = max(A_cur(8),0.0_RK)
               ELSE IF (self%cfg%snow_model == 1) THEN
                  state%precip = max(A_cur(7),0.0_RK)
               ENDIF
            ELSE
               CALL error('Wrong forcing type (must be 1, 2, 3, 4 or 5).')
            ENDIF
            state%uv10 = sqrt(state%u10**2 + state%v10**2) ! AG 2014

            self%param%p_air = psurf/100.0_RK !-WMEJ: p_air = psurf, add in Lake driver
            IF (self%cfg%forcing_mode /= 4) THEN ! Heat fluxes calculations (forcing 2, 3 and 5)
               IF (state%black_ice_h + state%white_ice_h == 0) THEN !Free water
                  ! Factor 0.6072 to account for changing wind height from 10 to 2 m
                  ! Further evaluation of evaporation algorithm may be required.
                  fu = sqrt((2.7_RK*max(0.0_RK,(T_surf-state%T_atm)/(1 - 0.378_RK*Vap_atm/self%param%p_air))**0.333_RK)**2 + (0.6072_RK*3.1_RK*state%uv10)**2)
                  fu = fu*self%param%p_windf ! Provided fitting factor p_windf (~1)

                  ! Water vapor saturation pressure in air at water temperature (Gill 1992) [millibar]
                  Vap_wat = 10**((0.7859_RK + 0.03477_RK*T_surf)/(1 + 0.00412_RK*T_surf))
                  Vap_wat = Vap_wat*(1 + 1e-6_RK*self%param%p_air*(4.5_RK + 0.00006_RK*T_surf**2))     !-WMEJ : p_air = psurf, add in Lake driver

                  ! Long-wave radiation according to Dilley and O'Brien
                  ! see Flerchinger et al. (2009)
                  IF (self%cfg%forcing_mode /= 5) THEN
                     H_A = (1 - r_a)*((1 - 0.84_RK*Cloud)*(59.38_RK + 113.7_RK*((state%T_atm + 273.15_RK)/273.16_RK)**6 &
                        + 96.96_RK*sqrt(465*Vap_atm/(state%T_atm + 273.15_RK)*0.04_RK))/5.67e-8_RK/ &
                        (state%T_atm + 273.15_RK)**4 + 0.84_RK*Cloud)*5.67e-8_RK*(state%T_atm + 273.15_RK)**4
                  ENDIF

                  H_A = H_A*self%param%p_lw ! Provided fitting factor p_radin (~1)

                  H_W = -emiss_water*stefnc*(T_surf + 273.15_RK)**4

                  ! Flux of sensible heat (convection)
                  H_K = -B0*fu*(T_surf - state%T_atm)

                  ! Flux of latent heat (evaporation, condensation)
                  H_V = -fu*(Vap_wat - Vap_atm)

                  ! Global heat flux (positive: air to water, negative: water to air)
                  state%heat = H_A + H_W + H_K + H_V + F_glob * self%param%beta_sol !MS: added term with beta_sol
                  ! Removal of solar short-wave radiation absorbed in first water cell
                  state%rad0 = F_glob * (1 - self%param%beta_sol) !MS: added beta_sol

                  state%heat_snow = 0 ! Heat snow
                  state%heat_snowice = 0 ! Heat snowice
                  state%heat_ice = 0 ! Heat ice

               ELSE IF (state%black_ice_h > 0 .or. state%white_ice_h > 0) THEN ! Ice Cover
                  ! Light penetration in snow, ice and snowice as well as wind blocking added 2018 by Love Raaman

                  IF (state%T_atm >= self%param%Freez_Temp) THEN! Melting occures when air temp above freezing point,
                     ! THEN activate surface heat fluxes following
                     ! Matti Leppäranta (2009), Modelling the Formation and Decay of Lake Ice DOI: 10.1007/978-90-481-2945-4_5 in book: The Impact of Climate Change on European Lakes
                     ! in G. George (Ed.), The Impact of Climate Change on European Lakes (pp. 63–83).
                     ! Dordrecht: Springer Netherlands. https://doi.org/10.1007/978-90-481-2945-4_5
                     ! and with corrections in
                     ! Leppäranta, M. (2014). Freezing of lakes and the evolution of their ice cover.
                     ! New York: Springer. ISBN 978-3-642-29080-0
                     IF (state%snow_h == 0) THEN ! Ice Cover (ice and snowice)
                        emissivity = emiss_ice
                     ELSE ! Snow Cover
                        ! Varies from 0.8 to 0.9 depending on snow density
                        emissivity = 5.0e-4_RK * state%snow_dens + 6.75e-1_RK
                     ENDIF
                     ! obs fitting factors self%param%p_lw and self%param%p_windf not applied to ice covered lake
                     IF (self%cfg%forcing_mode /= 5) THEN
                        H_A = (Ha_a + Ha_b * (Vap_atm**(1.0_RK/2.0_RK))) * (1 + Ha_c * Cloud**2) * stefnc * (state%T_atm + 273.15_RK)**4
                     ENDIF
                     H_W = -emissivity * stefnc * (T_surf + 273.15_RK)**4
                     H_K = rho_air * cp_air * Hk_CH * (state%T_atm - T_surf) * state%uv10
                     sh_c = 0.622/self%param%p_air! Converter from absolut vapour pressure to specific humidity (Lepparanta 2015)
                     qa  = sh_c * Vap_atm! Specific humidity air
                     q0  = sh_c * 6.11! Specific humidity for saturation levels (6.11 mbar) at surface (ice) at 0°C (Lepparanta 2015)
                     H_V = rho_air * (l_h + l_e) * Hv_CE * (qa - q0) * state%uv10 ! Through sublimation (solid to gas)
                  ELSE ! IF no melting only considering penetraiting solar radiation since ice formation is parameterised towards temperature
                     H_A = 0
                     H_W = 0
                     H_K = 0
                     H_V = 0
                  ENDIF

                  H_tot = H_A + H_W + H_K + H_V
                  IF (H_tot < 0) THEN
                     H_tot = 0
                  ENDIF

                  ! Global heat flux (positive: air to water, negative: water to air) !MS: added beta_sol ; LRV added lambda_snow_ice
                  ! Leppäranta, M. (2014), Eq. 6.12
                  ! Removal of solar short-wave radiation absorbed in snow, snowice, ice and first water cell (works also when x_h = 0)
                  state%heat = F_glob * exp(-lambda_snow*state%snow_h - lambda_snowice*state%white_ice_h - lambda_ice*state%black_ice_h) * self%param%beta_sol
                  state%rad0 = F_glob * exp(-lambda_snow*state%snow_h - lambda_snowice*state%white_ice_h - lambda_ice*state%black_ice_h) * (1 - self%param%beta_sol)

                  ! Heat flux into snow, ice or snowice layer.
                  ! Light absorption each layer
                  ! Leppäranta, M. (2014), Eq. 6.12
                  F_snow    = F_glob * (1 - exp(-lambda_snow*state%snow_h))
                  F_snowice = F_glob * (exp(-lambda_snow*state%snow_h) - exp(-lambda_snow*state%snow_h - lambda_snowice*state%white_ice_h))
                  F_ice     = F_glob * (exp(-lambda_snow*state%snow_h - lambda_snowice*state%white_ice_h) - exp(-lambda_snow*state%snow_h - lambda_snowice*state%white_ice_h - lambda_ice*state%black_ice_h))
                  IF (F_snow < 0 .or. F_snowice < 0 .or. F_ice < 0 .or. state%heat < 0 .or. state%rad0 < 0) THEN
                     CALL error('Negative heat flux not allowed for melting.')
                  ENDIF

                  !Light + other heat fluxes into top layer
                  IF (state%snow_h > 0 .and. state%white_ice_h > 0 .and. state%black_ice_h > 0) THEN ! snow, snowice and ice
                     state%heat_snow = H_tot + F_snow
                     state%heat_snowice = F_snowice
                     state%heat_ice = F_ice
                  ELSE IF (state%snow_h > 0 .and. state%white_ice_h > 0 .and. state%black_ice_h == 0) THEN! snow and snowice
                     state%heat_snow = H_tot + F_snow
                     state%heat_snowice = F_snowice
                     state%heat_ice = 0
                  ELSE IF (state%snow_h > 0 .and. state%white_ice_h == 0 .and. state%black_ice_h > 0) THEN! snow and ice
                     state%heat_snow = H_tot + F_snow
                     state%heat_snowice = 0
                     state%heat_ice = F_ice
                  ELSE IF (state%snow_h == 0 .and. state%white_ice_h > 0 .and. state%black_ice_h > 0) THEN! snowice and ice
                     state%heat_snow = 0
                     state%heat_snowice = H_tot + F_snowice
                     state%heat_ice = F_ice
                  ELSE IF (state%snow_h == 0 .and. state%white_ice_h == 0 .and. state%black_ice_h > 0) THEN! ice
                     state%heat_snow = 0
                     state%heat_snowice = 0
                     state%heat_ice = H_tot + F_ice
                  ELSE IF (state%snow_h == 0 .and. state%white_ice_h > 0 .and. state%black_ice_h == 0) THEN! snowice
                     state%heat_snow = 0
                     state%heat_snowice = H_tot + F_snowice
                     state%heat_ice = 0
                  ENDIF

                  ! Suppress wind turbulence with wind lid (heat flux affected by wind still active on snow and ice)
                  state%u10 = 0
                  state%v10 = 0
                  state%uv10 = 0
               ENDIF
               ! save for output, not used in calculations
               state%ha = H_A
               state%hw = H_W
               state%hk = H_K
               state%hv = H_V
            ELSE !Forcing mode 4
               state%heat = heat0 + F_glob*self%param%beta_sol !MS: added term with beta_sol
               state%rad0 = F_glob * (1 - self%param%beta_sol) !FB, 2021: added term with beta_sol
               state%heat_snow = 0 ! Heat snow
               state%heat_snowice = 0 ! Heat snowice
               state%heat_ice = 0 ! Heat ice
            ENDIF
            IF (self%cfg%ice_model == 0) THEN
               IF (state%T(self%grid%ubnd_vol) < 0 .and. state%heat < 0) THEN
                  state%heat = 0
                  state%T(self%grid%ubnd_vol) = 0
               ENDIF
            ENDIF
         ENDIF

         ! Drag coefficient as a function of wind speed (AG 2014)
         IF (self%cfg%wind_drag_model == 1) THEN ! Constant wind drag coefficient
            state%C10 = self%param%C10_constant
         ELSE IF (self%cfg%wind_drag_model == 2) THEN ! Ocean model
            state%C10 = self%param%C10_constant*(-0.000000712_RK*state%uv10**2 + 0.00007387_RK*state%uv10 + 0.0006605_RK)
         ELSE IF (self%cfg%wind_drag_model == 3) THEN ! Lake model (Wüest and Lorke 2003)
            IF (state%uv10 <= 0.1) THEN
               state%C10 = self%param%C10_constant*0.06215_RK
            ELSE IF (state%uv10 <= 3.85_RK) THEN
               state%C10 = self%param%C10_constant*0.0044_RK*state%uv10**(-1.15_RK)
            ELSE ! Polynomial approximation of Charnock's law
               state%C10 = self%param%C10_constant*(-0.000000712_RK*state%uv10**2 + 0.00007387_RK*state%uv10 + 0.0006605_RK)
            ENDIF
         ENDIF

         tau = state%C10*rho_air/rho_0*state%uv10**2
         state%u_taus = sqrt(tau)

         state%tx = state%C10*rho_air/rho_0*state%uv10*state%u10
         state%ty = state%C10*rho_air/rho_0*state%uv10*state%v10

         CALL LakStaParms ( &
               ! "in" arguments
               ! -------------------
               hwref     , htref   , hqref   , usurf   ,&
               vsurf     , tmref   , qmref   , tskin   ,&
               arhos     , psurf   , zcbcv   ,&
               ! "out" arguments
               ! -------------------
               zol       , rib     , wdm     , shfdt   ,&
               taux      , tauy    , t2m     , q2m     ,&
               u10m      , v10m    , fh2m    , fq2m    ,&
               fm10m     , fq10m   , fm      , fh      ,&
               fq        , ustar   , qstar   , tstar   ,&
               ram       , rah     , raw     , emis    )

         IF (tskin > tfrz )THEN
               htvp = hvap
         ELSE
               htvp = hsub
         ENDIF

         olrg   = -H_W 
         fsena  = -H_K 
         fseng  = -H_K
         lfevpa = -H_V 
         fevpa  = lfevpa / htvp
         fevpg  = fevpa
         F_glob = sabg
         fgrnd  = F_glob + H_A + H_W + H_K + H_V
         trad   = (olrg/stefnc)**0.25
      ELSE

         lb = snlay + 1
         CALL CoLML_SurfFlux ( &
               ! "in" arguments
               ! -------------------
               scwat      , snlay      , zopt       , betaopt    ,&
               fetchopt   , rlat       , dplak      , dtlak      ,& 
               stke1      , hwref      , htref      , hqref      ,&
               usurf      , vsurf      , tmref      , qmref      ,&
               arhos      , psurf      , sols       , soll       ,&
               solsd      , solld      , lwdns      , sabg       ,&
               hpbl       , dzlaktop   , zlakebot   , dzssbtop   ,&
               lktmptop   , tssbtop    , xwliqtop   , xwicetop   ,&
               icefrtop   , zcbcv      , ipatch     ,&
               ! "inout" arguments
               ! -------------------
               tskin      , z0m        , z0h        , z0q        ,&
               felak      , btpri      ,&
               ! "out" arguments
               ! -------------------
               fseng      , fevpg      , fsena      , fevpa      ,&
               lfevpa     , olrg       , fgrnd      , fgrnd1     ,&
               zol        , rib        , trad       , htvp       ,&
               emis       , wdm        , ram        , rah        ,&
               raw        , shfdt      , taux       , tauy       ,&
               t2m        , q2m        , u10m       , v10m       ,&
               fh2m       , fq2m       , fm10m      , fq10m      ,&
               fm         , fh         , fq         , ustar      ,&
               qstar      , tstar      , rhosnow    , u2m        )

         H_A = lwdns
         H_W = -olrg
         H_K = -fsena
         H_V = -lfevpa
         F_glob = sabg
         H_tot = H_A + H_W + H_K + H_V
         ! Global heat flux (positive: air to water, negative: water to air)
         state%heat = H_A + H_W + H_K + H_V + F_glob * btpri
         state%heat = fgrnd1

         ! Removal of solar short-wave radiation absorbed in first water cell
         state%rad0 = F_glob * (1 - btpri) !MS: added beta_sol

         ! Heat flux into snow, ice or snowice layer.
         ! Light absorption each layer
         ! Leppäranta, M. (2014), Eq. 6.12
         F_snow    = F_glob * (1 - exp(-lambda_snow*state%snow_h))
         F_snowice = F_glob * (exp(-lambda_snow*state%snow_h) - exp(-lambda_snow*state%snow_h - lambda_snowice*state%white_ice_h))
         F_ice     = F_glob * (exp(-lambda_snow*state%snow_h - lambda_snowice*state%white_ice_h) - exp(-lambda_snow*state%snow_h - lambda_snowice*state%white_ice_h - lambda_ice*state%black_ice_h))

         ! -WMEJ there is a problem with the heat flux, they are not correct
         ! IF (F_snow < 0 .or. F_snowice < 0 .or. F_ice < 0 .or. state%heat < 0 .or. state%rad0 < 0) THEN
         !     CALL error('Negative heat flux not allowed for melting.')
         ! ENDIF
         
         ! Light + other heat fluxes into top layer
         IF (state%snow_h > 0 .and. state%white_ice_h > 0 .and. state%black_ice_h > 0) THEN ! snow, snowice and ice
            state%heat_snow = H_tot + F_snow
            state%heat_snowice = F_snowice
            state%heat_ice = F_ice
         ELSE IF (state%snow_h > 0 .and. state%white_ice_h > 0 .and. state%black_ice_h == 0) THEN! snow and snowice
            state%heat_snow = H_tot + F_snow
            state%heat_snowice = F_snowice
            state%heat_ice = 0
         ELSE IF (state%snow_h > 0 .and. state%white_ice_h == 0 .and. state%black_ice_h > 0) THEN! snow and ice
            state%heat_snow = H_tot + F_snow
            state%heat_snowice = 0
            state%heat_ice = F_ice
         ELSE IF (state%snow_h == 0 .and. state%white_ice_h > 0 .and. state%black_ice_h > 0) THEN! snowice and ice
            state%heat_snow = 0
            state%heat_snowice = H_tot + F_snowice
            state%heat_ice = F_ice
         ELSE IF (state%snow_h == 0 .and. state%white_ice_h == 0 .and. state%black_ice_h > 0) THEN! ice
            state%heat_snow = 0
            state%heat_snowice = 0
            state%heat_ice = H_tot + F_ice
         ELSE IF (state%snow_h == 0 .and. state%white_ice_h > 0 .and. state%black_ice_h == 0) THEN! snowice
            state%heat_snow = 0
            state%heat_snowice = H_tot + F_snowice
            state%heat_ice = 0
         ENDIF

         state%u10 = usurf
         state%v10 = vsurf
         state%uv10 = sqrt(state%u10**2 + state%v10**2)
         state%T_atm = tmref
         state%precip = max(crain + csnow + lrain + lsnow,0.0_RK)

         ! save for output, not used in calculations
         state%ha = H_A
         state%hw = H_W
         state%hk = H_K
         state%hv = H_V

         ! Drag coefficient as a function of wind speed (AG 2014)
         IF (self%cfg%wind_drag_model == 1) THEN ! Constant wind drag coefficient
            state%C10 = self%param%C10_constant
         ELSE IF (self%cfg%wind_drag_model == 2) THEN ! Ocean model
            state%C10 = self%param%C10_constant*(-0.000000712_RK*state%uv10**2 + 0.00007387_RK*state%uv10 + 0.0006605_RK)
         ELSE IF (self%cfg%wind_drag_model == 3) THEN ! Lake model (Wüest and Lorke 2003)
            IF (state%uv10 <= 0.1) THEN
               state%C10 = self%param%C10_constant*0.06215_RK
            ELSE IF (state%uv10 <= 3.85_RK) THEN
               state%C10 = self%param%C10_constant*0.0044_RK*state%uv10**(-1.15_RK)
            ELSE ! Polynomial approximation of Charnock's law
               state%C10 = self%param%C10_constant*(-0.000000712_RK*state%uv10**2 + 0.00007387_RK*state%uv10 + 0.0006605_RK)
            ENDIF
         ENDIF

         !-WMEJ Is the 'taux' in CoLML_SurfFlux the same as 'tx' here? There is considerable controversy here.
         ! state%u_taus = ustar
         ! state%tx = taux  !state%C10*rho_air/rho_0*state%uv10*state%u10
         ! state%ty = tauy  !state%C10*rho_air/rho_0*state%uv10*state%v10
         tau = state%C10*rho_air/rho_0*state%uv10**2
         state%u_taus = sqrt(tau)
         state%tx = state%C10*rho_air/rho_0*state%uv10*state%u10
         state%ty = state%C10*rho_air/rho_0*state%uv10*state%v10

      ENDIF

   END SUBROUTINE


   ! Coriolis forces
   SUBROUTINE forcing_update_coriolis(self, state)
      IMPLICIT NONE
      class(ForcingModule) :: self
      class(ModelState) :: state
      real(RK) :: cori
      real(RK), DIMENSION(size(state%U)) :: u_temp
      associate (grid=>self%grid, dt=>state%dt, param=>self%param)

         ! Calculate u_taub before changing U resp V
         state%u_taub = sqrt(state%drag*(state%U(1)**2 + state%V(1)**2))

         ! Calculate coriolis parameter based on latitude
         cori = 2.0_RK*7.292e-5_RK*sin(param%Lat*pi/180.0_RK)

         ! Update state based on coriolis parameter
         u_temp = state%U

         state%U(1:grid%ubnd_vol) = state%U(1:grid%ubnd_vol)*cos(Cori*dt) + state%V(1:grid%ubnd_vol)*sin(Cori*dt)
         state%V(1:grid%ubnd_vol) = -u_temp(1:grid%ubnd_vol)*sin(Cori*dt) + state%V(1:grid%ubnd_vol)*cos(Cori*dt)

         RETURN
      END associate
   END SUBROUTINE

END MODULE  strat_forcing





!------------!
! Ice MODULE
! Following:
! ------ !
! Saloranta, T. M., & Andersen, T. (2007).
! MyLake—A multi-year lake simulation model code suitable for uncertainty and sensitivity analysis simulations.
! Ecological Modelling, 207(1), 45–60.
! ------ !
! Yen (1981)
! Review of thermal properties of snow, ice and sea ice
! ------ !
! Saloranta, T. M. (2000).
! Modeling the evolution of snow, snow ice and ice in the Baltic Sea. Tellus A, 52(1), 93–108.
!--------------------------------------------------
! Ice MODULE  constructed in 2018 by:
! Dr. Love Råman Vinnå
! Department Surface Waters Research & Management
! Eawag, Seestrase 79, 6047 Kastanienbaum, Switzerland
MODULE  strat_ice
   USE strat_kinds
   USE strat_forcing
   USE MOD_Lake_Const
   USE strat_simdata
   USE strat_grid
   IMPLICIT NONE

!-----------------------------------------------------------
PRIVATE
!-----------------------------------------------------------

   ! Common Types
   type, public :: IceModule
      class(ModelConfig), pointer :: model_cfg
      class(ModelParam), pointer :: model_param
      class(StaggeredGrid), pointer :: grid

   CONTAINS
      ! To initialize and run the ice MODULE  with subroutines
      procedure, pass :: init => ice_module_init
      procedure, pass :: update => ice_module_update

      ! The model itself
      procedure, pass :: do_ice_freezing        => ice_formation
      procedure, pass :: do_ice_melting         => ice_melting
      procedure, pass :: do_snow_melting        => snow_melting
      procedure, pass :: do_snow_build          => snow_build
      procedure, pass :: do_underneath_melting  => underneath_melting

   END type

!-----------------------------------------------------------
CONTAINS
!-----------------------------------------------------------

   ! Initiate ice model
   SUBROUTINE ice_module_init(self, model_cfg, model_param, grid)
      IMPLICIT NONE
      class(IceModule) :: self
      class(ModelConfig), target :: model_cfg
      class(ModelParam), target :: model_param
      class(StaggeredGrid), target :: grid

      self%model_cfg => model_cfg
      self%model_param => model_param
      self%grid => grid

   END SUBROUTINE


   ! SUBROUTINE ice_module_update
   SUBROUTINE ice_module_update(self, state, param)
      IMPLICIT NONE
      class(IceModule)  :: self
      class(ModelState) :: state
      class(ModelParam) :: param

      !-------------------
      ! Below freezing
      !-------------------
      IF (param%Freez_Temp >= state%T(self%grid%ubnd_vol) .and. (state%black_ice_h + state%white_ice_h) == 0) THEN
         ! Ice expanding (air temp & water temp less THEN Freez_Temp) no ice present
         CALL self%do_ice_freezing(state, param)
      ELSE IF ((state%black_ice_h + state%white_ice_h) > 0 .and. param%Freez_Temp >= state%T_atm) THEN
         ! IF ice exist and air temp < freez temp, initiate ice formation
         CALL self%do_ice_freezing(state, param)
      ENDIF

      ! Snow fall addition onto ice
      IF (self%model_cfg%snow_model == 1 .and. param%snow_temp >= state%T_atm .and. (state%black_ice_h + state%white_ice_h) > 0 .and. state%precip > 0) THEN
         CALL self%do_snow_build(state)
      ENDIF

      !------------------
      ! Above freezing
      !-------------------
      IF (param%Freez_Temp < state%T_atm .and. (state%black_ice_h + state%white_ice_h) > 0) THEN
         ! Melt snow
         IF (state%snow_h > 0 .and. param%snow_temp < state%T_atm .and. self%model_cfg%snow_model == 1) THEN
            CALL self%do_snow_melting(state)
         ENDIF
         ! Melt ice from above
         IF (state%white_ice_h + state%black_ice_h > 0) THEN
            CALL self%do_ice_melting(state, param)
         ENDIF
         ! Melt ice from underneath
         IF (state%T(self%grid%ubnd_vol) > param%freez_temp .and. (state%white_ice_h + state%black_ice_h) > 0) THEN
            CALL self%do_underneath_melting(state, param)
         ENDIF
      ENDIF

      ! Update total ice thickness
      state%total_ice_h = state%white_ice_h + state%black_ice_h

   END SUBROUTINE

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! The Ice Model
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! %%%%%%%%%%%%%%%%%%%%
    ! Below freezing point
    ! %%%%%%%%%%%%%%%%%%%%

   ! Ice/snowice formation and growth
   SUBROUTINE ice_formation(self, state, param)
      IMPLICIT NONE
      class(IceModule), intent(inout)  :: self
      class(ModelState), intent(inout) :: state
      class(ModelParam), intent(inout) :: param

      ! Define variables only used in ice_model
      real(RK) :: P
      real(RK) :: P_ice
      real(RK) :: P_snow
      real(RK) :: Freez_energy
      real(RK) :: snow_weight
      real(RK) :: buoyancy_ice
      real(RK) :: snow_height_ice_mass

      IF (state%black_ice_h == 0) THEN
         ! Energy released for first ice forming
         Freez_energy = (param%Freez_Temp - state%T(self%grid%ubnd_vol))  * cp * (self%grid%h(self%grid%ubnd_vol) * rho_0)  ![J/kg/K]*[K]*[kg] = [J]

         ! First ice thickness
         state%black_ice_h = Freez_energy / l_h / (ice_dens*1*1) ![J] / [J/kg] / [kg/m3] / [m2] = [m]

         ! Set surface temperature to freezing point
         state%T(self%grid%ubnd_vol) = param%Freez_Temp ![°C]

         ! Ice temp eq 17 in Saloranta et al. (2007).
         state%ice_temp = 0 ![°C]

      ELSE   ! When ice cover exists
         ! Snow and ice insulating effect eq 17 and 18 Saloranta et al. (2007).
         P_ice = 1 / (10 * state%black_ice_h)

         IF (self%model_cfg%snow_model == 1) THEN
            P_snow = (k_ice * state%snow_h)/(k_snow * state%black_ice_h)
         ELSE
            P_snow = 0
         ENDIF

         P = max(P_snow,P_ice)

         ! Ice temp eq 17 in Saloranta et al. (2007).
         state%ice_temp = (P*param%freez_temp + state%T_atm)/(1 + P)

         ! Snow-Ice formation, IF weight of snow exceeds ice buoyancy
         IF (self%model_cfg%snow_model == 1 .and. state%snow_h > 0) THEN
            snow_weight = state%snow_h * state%snow_dens*1*1 ! kg
            buoyancy_ice = state%black_ice_h*1*1 * (rho_0 - ice_dens) + state%white_ice_h*1*1 * (rho_0 - snowice_dens)!kg

            IF (snow_weight > buoyancy_ice) THEN
               snow_height_ice_mass = snow_weight - buoyancy_ice
               state%snow_h = state%snow_h - snow_height_ice_mass/state%snow_dens
               state%white_ice_h  = state%white_ice_h  + snow_height_ice_mass/state%snow_dens ! Assuming water from below fills up zone needed to be flooded to achieve buoyant stability
            ENDIF
         ENDIF
         ! Ice thickness eq. 16 in Saloranta et al. (2007).
         state%black_ice_h = sqrt(state%black_ice_h**2 + (2*k_ice)/(ice_dens * l_h) * (param%freez_temp - state%ice_temp) * state%dt) 

         ! Set surface temperature to freezing point
         state%T(self%grid%ubnd_vol) = param%freez_temp ![°C]
      ENDIF

      ! IF melt larger than ice height, put remaining energy to water and set ice height to zero
      IF (state%black_ice_h < 0) THEN
         state%heat = state%heat + (l_h * ice_dens * (-1 * state%black_ice_h) / state%dt) ![J/kg] * [kg/m3] * [m] / [s] = [J/sm2] = [W/m2]
         state%black_ice_h = 0
      ENDIF

   END SUBROUTINE


   ! Snow layer build-up
   SUBROUTINE snow_build(self, state)
      IMPLICIT NONE
      class(IceModule) :: self
      class(ModelState) :: state
      real(RK) :: snow_h_new
      real(RK) :: Ws, ChangSnowDens

      ! Calculate new snow height
      snow_h_new = (state%precip / 3600) * state%dt * (rho_0 / rho_s_0) ! Go from m/h to m/s and increase volume from water to snow

      ! Calculate snow density due to compression of snow layer, Yen (1981)
      Ws = (snow_h_new * rho_s_0) / (rho_0 * 1 * 1)  ! [m water equivalent]
      ! Compression Yen (1981) eq. 7
      ChangSnowDens = state%snow_dens * C01 * Ws * exp(-C02 * state%snow_dens) * state%dt
      ! Compress old snow layer
      state%snow_h = state%snow_h * (state%snow_dens /(state%snow_dens + ChangSnowDens))

      ! Adjust density
      state%snow_dens = state%snow_dens + ChangSnowDens

      ! Combine the old and the new snow layer
      ! Create new density from old and new density
      state%snow_dens  = (state%snow_dens * state%snow_h + rho_s_0 * snow_h_new) / (state%snow_h + snow_h_new)
      ! Update snow height
      state%snow_h = state%snow_h + snow_h_new

      ! Maximum allowed snow density
      IF (state%snow_dens > rho_s_max) THEN
         state%snow_dens = rho_s_max
      ENDIF

   END SUBROUTINE

    ! %%%%%%%%%%%%%%%%%%%%
    ! Above freezing point
    ! %%%%%%%%%%%%%%%%%%%%

   ! Snow melting from above through sublimation (l_e ~= 0) or non-sublimation (l_e = 0) (i.e. solid to gas (l_h + l_e))
   SUBROUTINE snow_melting(self, state)
      IMPLICIT NONE
      class(IceModule) :: self
      class(ModelState) :: state
      ! Define variables only used in snow_model
      real(RK) :: Melt_energy1
      real(RK) :: MeltHeight1
      real(RK) :: Melt_ice
      real(RK) :: MeltHeightIce

      ! Melt Snow from atmosphere
      Melt_energy1 = state%heat_snow * state%dt * 1 * 1 ! [W/m2] * [s] * [m2] = [J]
      MeltHeight1 = Melt_energy1 / (l_h + l_e) / (state%snow_dens * 1 * 1) ! [J] / [J/kg] / [kg/m3] / [m2] = [m]
      IF (MeltHeight1 < 0) THEN ! DO not add snow from above, only melt snow
         Melt_energy1 = 0
         MeltHeight1 = 0
      ENDIF
      ! Change the snow height
      state%snow_h = state%snow_h - MeltHeight1! [m]

      ! IF melting energy larger than required, put remaining energy to melting of snowice / ice and set snow height to zero
      IF (state%snow_h < 0) THEN
         Melt_ice = (l_h + l_e) * state%snow_dens * (1 * 1) * (-1 * state%snow_h)! [J/kg]  [kg/m3]  [m2]  [m] = [J]
         state%snow_h = 0
         IF (state%white_ice_h > 0) THEN
            MeltHeightIce = Melt_ice / (l_h + l_e) / (snowice_dens * 1 * 1) ! [J] / [J/kg] / [kg/m3] / [m2] = [m]
            state%white_ice_h = state%white_ice_h - MeltHeightIce! [m]
            IF (state%white_ice_h < 0) THEN
               Melt_ice = (l_h + l_e) * snowice_dens * (1 * 1) * (-1 * state%white_ice_h)! [J/kg]  [kg/m3]  [m2]  [m] = [J]
               state%white_ice_h = 0
               MeltHeightIce = Melt_ice / (l_h + l_e) / (ice_dens * 1 * 1) ! [J] / [J/kg] / [kg/m3] / [m2] = [m]
               state%black_ice_h = state%black_ice_h - MeltHeightIce ! [m]
               ! IF melt larger than ice height, put remaining energy to water and set ice height to zero
               IF (state%black_ice_h < 0) THEN
                  state%heat = state%heat + (l_h * ice_dens * (-1 * state%black_ice_h) / state%dt) ! [J/kg] * [kg/m3] * [m] / [s] = [J/sm2] = [W/m2]
                  state%black_ice_h = 0
               ENDIF
            ENDIF
         ELSE IF (state%white_ice_h == 0 .and. state%black_ice_h > 0 ) THEN
            MeltHeightIce = Melt_ice / (l_h + l_e) / (ice_dens * 1 * 1) ! [J] / [J/kg] / [kg/m3] / [m2] = [m]
            state%black_ice_h = state%black_ice_h - MeltHeightIce ! [m]
            ! IF melt larger than ice height, put remaining energy to water and set ice height to zero
            IF (state%black_ice_h < 0) THEN
               state%heat = state%heat + (l_h  * ice_dens * (-1 * state%black_ice_h) / state%dt) ! [J/kg] * [kg/m3] * [m] / [s] = [J/sm2] = [W/m2]
               state%black_ice_h = 0
            ENDIF
         ENDIF
      ELSE
         Melt_ice = 0
      ENDIF
      state%snwml = - (MeltHeightIce + MeltHeight1 ) / state%dt * 1000._RK ! [mm/s]
   END SUBROUTINE


   ! Ice melting from above through sublimation (l_e ~= 0) or non-sublimation (l_e = 0) (i.e. solid to gas (l_h + l_e))
   SUBROUTINE ice_melting(self, state, param)
      IMPLICIT NONE
      class(IceModule) :: self
      class(ModelState) :: state
      class(ModelParam) :: param
      ! Define variables only used in ice_model
      real(RK) :: Melt_energy2, Melt_energy3
      real(RK) :: MeltHeight2, MeltHeight3
      real(RK) :: Melt_ice, MeltHeightIce

      ! Melt snowice from atmosphere
      Melt_energy2 = state%heat_snowice * state%dt * 1 * 1 ! [W/m2] * [s] * [m2] = [J]
      ! Melt ice from atmosphere
      Melt_energy3 = state%heat_ice * state%dt * 1 * 1 ! [W/m2] * [s] * [m2] = [J]

      ! Snowice
      IF (state%snow_h == 0) THEN ! Free surface
         MeltHeight2 = Melt_energy2 / (l_h + l_e) / (snowice_dens * 1 * 1) ! [J] / [J/kg] / [kg/m3] / [m2] = [m]
      ELSE IF (state%snow_h > 0 .and. state%black_ice_h == 0) THEN ! Layer above and none below
         MeltHeight2 = Melt_energy2 / (l_h) / (snowice_dens * 1 * 1) ! [J] / [J/kg] / [kg/m3] / [m2] = [m]
      ELSE ! Layer above and below, melt snow
         state%heat_snow = Melt_energy2 / state%dt
         CALL self%do_snow_melting(state)
         MeltHeight2 = 0
      ENDIF
      IF (MeltHeight2 < 0) THEN ! DO not add ice, only melt.
         Melt_energy2 = 0;
         MeltHeight2 = 0;
      ENDIF
      ! Ice
      IF (state%snow_h + state%white_ice_h == 0) THEN ! Free surface
         MeltHeight3 = Melt_energy3 / (l_h + l_e) / (ice_dens * 1 * 1) ! [J] / [J/kg] / [kg/m3] / [m2] = [m]
      ELSE IF (state%snow_h + state%white_ice_h > 0) THEN ! Layer above and none below
         MeltHeight3 = Melt_energy3 / (l_h) / (ice_dens * 1 * 1) ! [J] / [J/kg] / [kg/m3] / [m2] = [m]
      ENDIF
      IF (MeltHeight3 < 0) THEN ! DO not add ice from above, only melt
         Melt_energy3 = 0;
         MeltHeight3 = 0;
      ENDIF

      ! New snowice height
      state%white_ice_h = state%white_ice_h - MeltHeight2! [m]
      MeltHeightIce = 0
      ! IF melting energy larger than required, put remaining energy to melting of ice and set snowice height to zero
      IF (state%white_ice_h < 0 .and. state%black_ice_h > 0) THEN
         Melt_ice = (l_h + l_e) * snowice_dens * (1 * 1) * (-1 * state%white_ice_h) ! [J/kg]  [kg/m3] [m2] [m] = [J]
         state%white_ice_h = 0
         MeltHeightIce = Melt_ice / (l_h + l_e) / (ice_dens * 1 * 1) ! [J] / [J/kg] / [kg/m3] / [m2] = [m]
      ELSE IF (state%white_ice_h < 0 .and. state%black_ice_h == 0) THEN ! IF melt larger than snowice hight, put remaining energy to water and set snowice hight to zero
         state%heat = state%heat + (l_h * snowice_dens * (-1 * state%white_ice_h) / state%dt) ! [J/kg] * [kg/m3] * [m] / [s] = [J/sm2] = [W/m2]
         state%white_ice_h = 0
      ENDIF

      ! New ice height
      state%black_ice_h = state%black_ice_h - MeltHeight3 - MeltHeightIce ! [m]
      ! IF melt larger than ice height, put remaining energy to water and set ice height to zero
      IF (state%black_ice_h < 0) THEN
         state%heat = state%heat + (l_h * ice_dens * (-1 * state%black_ice_h) / state%dt) ! [J/kg] * [kg/m3] * [m] / [s] = [J/sm2] = [W/m2]
         state%black_ice_h = 0
      ENDIF

      ! Set ice temp to freez point
      state%ice_temp = 0 ![°C]

      state%iceml = - (MeltHeightIce + MeltHeight3 + MeltHeight2 ) / state%dt * 1000._RK ! [mm/s]
   END SUBROUTINE


   ! Melt ice from below (l_h), while keeping freez temperature in the interface between ice and water
   SUBROUTINE underneath_melting(self, state, param)
      IMPLICIT NONE
      class(IceModule) :: self
      class(ModelState) :: state
      class(ModelParam) :: param
      ! Define variables only used in snow_model
      real(RK) :: Melt_energy2
      real(RK) :: MeltHeight2

      Melt_energy2 = (state%T(self%grid%ubnd_vol) - param%freez_temp)  * cp * (self%grid%h(self%grid%ubnd_vol) * rho_0)  ! [J/kg/K]*[K]*[kg] = [J]
      IF (state%black_ice_h > 0) THEN
         MeltHeight2 = Melt_energy2 / l_h / (ice_dens * 1 * 1) ! [J] / [J/kg] / [kg/m3] / [m2] = [m]
      ELSE IF (state%black_ice_h == 0 .and. state%white_ice_h > 0) THEN
         MeltHeight2 = Melt_energy2 / l_h / (snowice_dens * 1 * 1) ! [J] / [J/kg] / [kg/m3] / [m2] = [m]
      ENDIF
      IF (MeltHeight2 < 0) THEN !DO not add ice from below, only melt
         Melt_energy2 = 0;
         MeltHeight2 = 0;
      ENDIF

      ! Set surface water temperature to freezing point
      state%T(self%grid%ubnd_vol) = param%freez_temp ![°C]

      ! New ice height
      IF (state%black_ice_h > 0) THEN
         state%black_ice_h = state%black_ice_h - MeltHeight2 ! [m]
      ELSE IF (state%black_ice_h == 0 .and. state%white_ice_h > 0) THEN
         state%white_ice_h = state%white_ice_h - MeltHeight2 ! [m]
      ENDIF

      ! IF melt larger than ice height, put remaining energy to water and set ice height to zero
      IF (state%black_ice_h < 0) THEN
         state%heat = state%heat + (l_h * ice_dens * (-1 * state%black_ice_h) / state%dt) ! [J/kg] * [kg/m3] * [m] / [s] = [J/sm2] = [W/m2]
         state%black_ice_h = 0
      ENDIF
      IF (state%white_ice_h < 0) THEN ! Open water
         state%heat = state%heat + l_h * (snowice_dens * (-1 * state%white_ice_h)) / state%dt ! [J/kg] * [kg/m3] * [m] / [s] = [J/sm2] = [W/m2]
         state%white_ice_h = 0
      ENDIF
      IF (state%white_ice_h == 0 .and. state%black_ice_h == 0 .and. state%snow_h > 0) THEN
         state%heat = state%heat + l_h * (state%snow_dens * (state%snow_h)) / state%dt ! [J/kg] * [kg/m3] * [m] / [s] = [J/sm2] = [W/m2]
         state%snow_h = 0
         state%snow_dens = rho_s_0
      ENDIF
   END SUBROUTINE

END MODULE  strat_ice





!<    +---------------------------------------------------------------+
!     | Implementation of a statevar for a temperature variable
!<    +---------------------------------------------------------------+
MODULE  strat_temp
   USE strat_kinds, only: RK
   USE MOD_Lake_Const
   USE strat_simdata
   USE strat_statevar
   USE strat_grid
   USE strat_solver
   IMPLICIT NONE

!-----------------------------------------------------------
   PRIVATE
!-----------------------------------------------------------

   type, extends(ModelVariable), public :: TempModelVar
   CONTAINS
      procedure, pass(self), public :: calc_terms => temp_var_calc_terms
      procedure, pass(self), public :: post_solve => temp_var_post_solve
   END type

!-----------------------------------------------------------
CONTAINS
!-----------------------------------------------------------

! Only calc_terms overwritten
SUBROUTINE temp_var_calc_terms(self, state, param, sources, boundaries)
   class(TempModelVar), intent(inout) :: self
   class(ModelState), intent(inout) :: state
   class(ModelParam), intent(inout) :: param
   real(RK), DIMENSION(:) ::  sources, boundaries
   integer :: i
   associate (grid=>self%grid, &
               ubnd_fce=>self%grid%ubnd_fce, &
               ubnd_vol=>self%grid%ubnd_vol)

      !!!!!!!! Precalculations !!!!!!!!
      ! Radiation reaching top layer
      state%rad(self%grid%ubnd_fce) = state%rad0/rho_0/cp ![°C*m/s]  , rad is on faces

      ! Radiation reaching a layer is equal to radiation in the layer above minus absorption
      DO i = ubnd_fce - 1, 1, -1
         state%rad(i) = state%rad(i + 1)*exp(-grid%h(i)*(state%absorb(ubnd_fce - i)+state%absorb(ubnd_fce + 1 - i))/2) !Attenuated by absorption
      ENDDO

      !!!!!!!! Define sources !!!!!!!!
      ! Add Hsol Term to sources (Eq 1, Goudsmit(2002))
      sources(1:ubnd_vol) = (state%rad(2:ubnd_fce) - state%rad(1:ubnd_fce - 1))/grid%h(1:ubnd_vol)

      ! Set boundary heat flux at surface (Eq 25, Goudsmit(2002))
      sources(ubnd_vol) = sources(ubnd_vol) + state%heat/rho_0/cp/grid%h(ubnd_vol)

      ! No explicit boundary conditions
      boundaries(1:ubnd_vol) = 0

      ! Forcing mode 1 for temp is done in post solve method
      IF (self%cfg%forcing_mode==1) THEN
      !   sources(ubnd_vol) = (state%SST-state%T(ubnd_vol))/state%dt
      !   boundaries(ubnd_vol) = 0
      !    bu(nz) = 1.0_dp
      !    au(nz) = 0.0_dp
      !    cu(nz) = 0.0_dp
      !    du(nz) = SST
      ENDIF

      ! Add geothermal flux to sources (Eq 1, Goudsmit(2002))
      IF (param%fgeo /= 0) sources(1:ubnd_vol) = sources(1:ubnd_vol) + state%fgeo_add(1:ubnd_vol)
   END associate
END SUBROUTINE


   SUBROUTINE temp_var_post_solve(self, state)
      class(TempModelVar), intent(inout) :: self
      class(ModelState), intent(inout) :: state

      IF (self%cfg%forcing_mode==1) THEN
         state%T(self%grid%ubnd_vol) = state%SST
      ENDIF

   END SUBROUTINE

END MODULE  strat_temp





!     +---------------------------------------------------------------+
!     | Inputfile  MODULE
!     |  - Reads configuration and initial conditions
!     |  - Sets up simulation data structure!
!     +---------------------------------------------------------------+
MODULE  strat_inputfile
   USE strat_kinds, only: RK
   USE strat_simdata
   USE strat_grid
   USE MOD_Lake_Const
   USE strat_utilities
   IMPLICIT NONE

!-----------------------------------------------------------
PRIVATE
!-----------------------------------------------------------

   !##################################################
   !# Inputfile
   !##################################################
   type, public :: SimstratSimulationFactory
      PRIVATE
      class(SimulationData), pointer :: simdata

   CONTAINS
      procedure, pass(self), public :: initialize => initialize_model
      procedure, pass(self), public :: setup_modcof
      procedure, pass(self), public :: setup_model
   END type SimstratSimulationFactory

!-----------------------------------------------------------
CONTAINS
!-----------------------------------------------------------

   SUBROUTINE initialize_model( self,&
               ! "inout" arguments
               ! ---------------------------
               simdata        , euler_i_disc   , euler_i_disc_keps, mod_forcing   ,&
               mod_stability  , mod_temperature, mod_u            , mod_v         ,&
               mod_k          , mod_eps        , mod_s            , mod_turbulence,&
               mod_ice        , mod_absorption , mod_advection    ,&
               ! "in" arguments
               ! ---------------------------
               nlake          , dtlak          , dplak            , ziarea        ,&
               zilak          , zlake          , dzlak            )

! ---------------------------------- code history -------------------------------------
! Description:  Initialize Simstrat, 
!               Replace the original initialize_model SUBROUTINE
!
! Original author: Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024
! -------------------------------------------------------------------------------------
   USE strat_simdata, only: SimulationData
   USE strat_forcing
   USE strat_utilities
   USE strat_stability, only: StabilityModule
   USE strat_windshear
   USE strat_statevar
   USE strat_temp
   USE strat_solver
   USE strat_discretization
   USE strat_keps
   USE strat_turbulence
   USE strat_ice
   USE strat_transport
   USE strat_absorption
   USE strat_advection
   USE strat_lateral
   USE strat_grid
   USE, intrinsic :: ieee_arithmetic

   IMPLICIT NONE
!================================================================================
   ! Instantiate all modules
   ! note that some are pointers/targets for polymorphism reasons
   class(SimstratSimulationFactory) :: self
   class(SimulationData), pointer, intent(out) :: simdata
   type(ThomasAlgSolver) :: solver
   type(EulerIDiscretizationMFQ) , intent(out) :: euler_i_disc
   type(EulerIDiscretizationKEPS) , intent(out) :: euler_i_disc_keps
   type(ForcingModule) , intent(out) :: mod_forcing
   type(StabilityModule), intent(out)  :: mod_stability
   type(TempModelVar) , intent(out) :: mod_temperature
   type(UVModelVar) , intent(out) :: mod_u, mod_v
   type(KModelVar) , intent(out) :: mod_k
   type(EpsModelVar) , intent(out) :: mod_eps
   type(TransportModVar) , intent(out) :: mod_s
   type(TurbulenceModule), intent(out)  :: mod_turbulence
   type(IceModule) , intent(out) :: mod_ice
   type(AbsorptionModule), intent(out)  :: mod_absorption
   type(AdvectionModule), intent(out)  :: mod_advection
   type(LateralModule), target :: mod_lateral_normal
   type(LateralRhoModule), target :: mod_lateral_rho
   class(GenericLateralModule), pointer :: mod_lateral
   type(GridConfig) :: grid_config
!  ------------------------- input variables ---------------------------
   integer, intent(in)      :: &
      nlake                     ! Number of layers

   real(RK), intent(in)     :: &
      dplak                   ,&! lake depth [m]
      dtlak                     ! Time step [s]

!  ------------------------- inout variables ---------------------------
   real(r8), intent(inout)  :: &
      zlake(nlake)            ,&! Depth of each node layer [m]
      dzlak(nlake)            ,&! center depth of each layer [m]
      zilak(nlake+1)          ,&! Depth of each layer interface [m]
      ziarea(nlake+1)           ! Morphology of each layer interface [m2]

!  ------------------------- local variables ---------------------------
   real(RK) :: z_zero           ! New depth of the lake [m]
   integer ::  i 
!================================================================================

      allocate (SimulationData :: self%simdata)
      simdata => self%simdata
      allocate (grid_config%z_A_read(nlake+1), grid_config%A_read(nlake+1),grid_config%grid_read(nlake+1))

      ! ----- Set up grid
      ! Variable spacing
      grid_config%nz_grid = nlake
      grid_config%equidistant_grid = .FALSE.
      ! Reverse order of values, because the grid is defined from bottom to top 

      grid_config%z_A_read = zilak                     ! z-coordinate is positive upwards, zero point is at reservoir bottom
      grid_config%A_read   = ziarea
      grid_config%grid_read = -zilak
      grid_config%max_depth = grid_config%z_A_read(1) - grid_config%z_A_read(nlake+1)     !-WMEJ: depth = max - min depth
      CALL reverse_in_place (grid_config%grid_read)
      CALL reverse_in_place (grid_config%A_read)
      
      self%simdata%model_cfg%max_length_input_data  = 1000   !-WMEJ not used

      ! Initialize Grid of simdata
      CALL simdata%grid%init(grid_config)

      ! Initialize Model of simdata
      CALL simdata%model%init(simdata%grid%nz_grid)           !-WMEJ  -----> IF in doubt, please pay attention

      ! Update actual filled z in grid
      z_zero = 0.0
      CALL simdata%grid%update_depth(z_zero)

      ! Set initial lake level
      simdata%grid%lake_level = simdata%grid%z_face(simdata%grid%ubnd_fce)    ! grid%z_face set in strat_grid ,grid%ubnd_fce set in strat_grid
      simdata%grid%lake_level_old = simdata%grid%lake_level           

      ! Set initial data
      CALL simdata%grid%update_area_factors()
      simdata%sim_cfg%timestep = dtlak

      !-WMEJ Model parameter settings
      CALL self%setup_modcof()

      ! Init rest of model
      CALL self%setup_model()

      !-WMEJ: Initialization of all modules should be done here
      ! Initialize Discretization     
      CALL euler_i_disc%init(simdata%grid)
      CALL euler_i_disc_keps%init(simdata%grid)

      ! Initialize forcing MODULE   
      mod_forcing%cfg => simdata%model_cfg
      mod_forcing%param => simdata%model_param
      mod_forcing%grid => simdata%grid
      IF(mod_forcing%cfg%snow_model == 1) CALL warn('The snow MODULE  is turned on. This MODULE  needs precipitation data, note that the last column in the forcing file will be interpreted as precipitation.')

      ! Initialize absorption MODULE        
      CALL mod_absorption%init(simdata%model_cfg, &
                                 simdata%model_param, &
                                 simdata%grid)

      !-------------------------------------------------------------------
      ! +WMEJ in the current version, we have disabled inflow and outflow
      !-------------------------------------------------------------------
      ! IF there is advection (due to inflow)                                                                                                         
      IF (simdata%model_cfg%has_advection) THEN
         ! Initialize advection MODULE
         CALL mod_advection%init(simdata%model_cfg, &
                              simdata%model_param, &
                              simdata%grid)

         ! Initialize lateral MODULE  based on configuration
         IF (simdata%model_cfg%inflow_placement == 1) THEN
            ! Gravity based inflow
            mod_lateral => mod_lateral_rho
         ELSE
            ! User defined inflow depths
            mod_lateral => mod_lateral_normal
         ENDIF
         CALL mod_lateral%init(simdata%model_cfg, &
                              simdata%model_param, &
                              simdata%grid)
      ENDIF
      !-------------------------------------------------------------------

      ! Initialize simulation modules    
      CALL mod_stability%init(simdata%grid, simdata%model_cfg, simdata%model_param)
      CALL mod_turbulence%init(simdata%grid, simdata%model_cfg, simdata%model_param)
      CALL mod_ice%init(simdata%model_cfg, simdata%model_param, simdata%grid)

      ! Set temperature state var to have nu_h as nu and T as model variable    
      CALL mod_temperature%init(simdata%model_cfg, simdata%grid, solver, euler_i_disc, simdata%model%nuh, simdata%model%T, simdata%grid%ubnd_vol)

      ! Set U and V var to have num as nu and U reps V as model variable
      ! also, assign shear stress in model for this variable
      CALL mod_u%init(simdata%model_cfg, simdata%grid, solver, euler_i_disc, simdata%model%num, simdata%model%U, simdata%grid%ubnd_vol)
      CALL mod_u%assign_shear_stress(simdata%model%tx)

      CALL mod_v%init(simdata%model_cfg, simdata%grid, solver, euler_i_disc, simdata%model%num, simdata%model%V, simdata%grid%ubnd_vol)
      CALL mod_v%assign_shear_stress(simdata%model%ty)

      ! Set mod_s (transport MODULE) to have nuh as nu and to manipulate S based on dS
      CALL mod_s%init(simdata%model_cfg, simdata%grid, solver, euler_i_disc, simdata%model%nuh, simdata%model%S, simdata%grid%ubnd_vol)
      CALL mod_s%assign_external_source(simdata%model%dS)

      ! Set up K and eps state vars with keps discretization and avh as nu
      CALL mod_k%init(simdata%model_cfg, simdata%grid, solver, euler_i_disc_keps, simdata%model%avh, simdata%model%K, simdata%grid%ubnd_fce)
      CALL mod_eps%init(simdata%model_cfg, simdata%grid, solver, euler_i_disc_keps, simdata%model%avh, simdata%model%eps, simdata%grid%ubnd_fce)

   END SUBROUTINE initialize_model


    ! Set up model configuration and model parameters
   SUBROUTINE setup_modcof(self)
! ---------------------------------- code history -------------------------------------
! Description:  Set up model configuration and model parameters, 
!               Replace the original way of setting parameters by reading files
!
! Original author: Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024
! -------------------------------------------------------------------------------------
      USE MOD_Lake_Const, only: CoupleAED2, SplitSeicheParameter, UseFilteredWind, Has_advection,&
                                       TurbulenceModel, StabilityFunction, FluxCondition,&
                                       ForcingMode, SeicheNormalization, WindDragModel,&
                                       InflowPlacement, PressureGradients, IceModel, SnowModel,&
                                       a_seiche, a_seiche_w, strat_sumr, q_nn, f_wind, c10, cd,&
                                       hgeo, p_sw, p_lw, p_windf, beta_sol, freez_temp, snow_temp
      USE strat_simdata, only: SimulationData
!================================================================================
      class(SimulationData), pointer:: simdata
      class(SimstratSimulationFactory) :: self
      type(GridConfig) :: grid_config
        
      self%simdata%model_cfg%couple_aed2               = CoupleAED2                        
      self%simdata%model_cfg%split_a_seiche            = SplitSeicheParameter
      self%simdata%model_cfg%use_filtered_wind         = UseFilteredWind
      self%simdata%model_cfg%has_advection             = Has_advection    
      self%simdata%model_cfg%turbulence_model          = TurbulenceModel                   
      self%simdata%model_cfg%stability_func            = StabilityFunction
      self%simdata%model_cfg%flux_condition            = FluxCondition
      self%simdata%model_cfg%forcing_mode              = ForcingMode
      self%simdata%model_cfg%seiche_normalization      = SeicheNormalization
      self%simdata%model_cfg%wind_drag_model           = WindDragModel
      self%simdata%model_cfg%inflow_placement          = InflowPlacement
      self%simdata%model_cfg%pressure_gradients        = PressureGradients
      self%simdata%model_cfg%ice_model                 = IceModel
      self%simdata%model_cfg%snow_model                = SnowModel

      ! self%simdata%model_param%Lat                     = Slat
      self%simdata%model_param%a_seiche                = a_seiche
      self%simdata%model_param%a_seiche_w              = a_seiche_w
      self%simdata%model_param%strat_sumr              = strat_sumr
      self%simdata%model_param%q_NN                    = q_nn
      self%simdata%model_param%f_wind                  = f_wind
      self%simdata%model_param%C10_constant            = c10
      self%simdata%model_param%CD                      = cd
      self%simdata%model_param%fgeo                    = hgeo
      self%simdata%model_param%p_sw                    = p_sw
      self%simdata%model_param%p_lw                    = p_lw
      self%simdata%model_param%p_windf                 = p_windf
      self%simdata%model_param%beta_sol                = beta_sol
      self%simdata%model_param%p_albedo                = p_albedo
      self%simdata%model_param%freez_temp              = freez_temp
      self%simdata%model_param%snow_temp               = snow_temp

   END SUBROUTINE setup_modcof


   ! Setup model configuration and state vars
   SUBROUTINE setup_model(self)
      IMPLICIT NONE
      class(SimstratSimulationFactory) :: self
      ! integer :: i
      associate ( simdata=>self%simdata, &
                  model_cfg=>self%simdata%model_cfg, &
                  model_param=>self%simdata%model_param, &
                  model=>self%simdata%model, &
                  grid=>self%simdata%grid)

         ! Initialize some more values
         IF (model_cfg%stability_func == 1) model%cm0 = 0.5625_RK
         IF (model_cfg%stability_func == 2) model%cm0 = 0.556171_RK
         model%cde = model%cm0**3
         sig_e = (kappa/model%cm0)**2/(ce2 - ce1)

         model%num(1:grid%nz_grid + 1) = 0.0_RK
         model%nuh(1:grid%nz_grid + 1) = 0.0_RK

         model%tx = 0.0_RK
         model%ty = 0.0_RK

         model%drag = (kappa/log(1.0_RK + 30/K_s*grid%h(1)/2))**2

         model%gamma = grid%Az(grid%ubnd_fce)/(grid%volume**1.5_RK)/sqrt(rho_0)*model_param%CD

         ! Geothermal heat flux
         IF (model_param%fgeo /= 0) THEN
            allocate(model%fgeo_add(grid%nz_grid))
            model%fgeo_add(1:grid%nz_grid) = model_param%fgeo/rho_0/cp*grid%dAz(1:grid%nz_grid)/grid%Az(2:grid%nz_grid + 1) ! calculation per kg
            IF (grid%Az(1) /= 0) THEN
               model%fgeo_add(1) = model%fgeo_add(1) + 2*model_param%fgeo/rho_0/cp*grid%Az(1)/((grid%Az(1) + grid%Az(2))*grid%h(1))
            ENDIF
         ENDIF

         ! Set up timing
         model%dt = self%simdata%sim_cfg%timestep
         model%simulation_time_old(2) = -model%dt     ! This is a workaround for correct implementation of simulation time, FB 2020
         model%simulation_time = 0
      END associate
   END SUBROUTINE setup_model

END MODULE  strat_inputfile




MODULE  MOD_Lake_Simstrat

   IMPLICIT NONE

!--------------------------------------------------------------------------
CONTAINS
!--------------------------------------------------------------------------

   SUBROUTINE Lake_Simstrat( &
               ! "in" arguments    
               ! -------------------
               nlake     , nsnow     , scwat     , nsoil     ,&
               zopt      , betaopt   , fetchopt  , etaopt    ,&
               dtlak     , xlat      , xlon      , dplak     ,&
               hwref     , htref     , hqref     , tmref     ,&
               usurf     , vsurf     , qmref     , arhos     ,&
               sols      , soll      , solsd     , solld     ,&
               sabg      , lwdns     , psurf     , hpbl      ,&
               crain     , csnow     , lrain     , lsnow     ,&
               zcbcv     , ipatch    ,&
               ! "inout" arguments
               ! ---------------------------
               dzlak     , zlake     , zilak     , lktmp     ,&
               dzssb     , zssb      , zissb     , t_ssb     ,&
               snlay     , snwcv     , snwag     , snwdp     ,&
               stke1     , tskin     , xwliq     , xwice     ,&
               z0m       , z0h       , z0q       , felak     ,&
               gamma     , etal      , btpri     , frlak     ,&
               icefr     , icedp     , tmice     , ziarea    ,&
               uwatv     , vwatv     , lksal     , tke       ,&
               eps       , etke      , num       , nuh       ,&
               bicedp    , wicedp    , lkrho    , rhosnw    ,&
               ! "out" arguments
               ! ---------------------------
               fsena     , fevpa     , lfevpa    , fseng     ,&
               fevpg     , olrg      , fgrnd     , trad      ,&
               qseva     , qsubl     , qsdew     , qfros     ,&
               taux      , tauy      , ustar     , qstar     ,&
               tstar     , emis      , zol       , snwml     ,&
               rib       , ram       , rah       , raw       ,&
               wdm       , t2m       , q2m       , u10m      ,&
               v10m      , fm10m     , fq10m     , fh2m      ,&
               fq2m      , fm        , fq        , fh        ,&
               shfdt     , urban_call)

! ---------------------------------- code history -------------------------------------
! Description:
!     This MODULE is the main MODULE of Simstrat, responsible for controlling the Simstrat model.
!     Simstrat counts layers starting from the lake bottom upwards,
!     so model data needs to be reversed before running Simstrat. 
!     The initialize MODULE must be used to initialize the Simstrat model, 
!     and the update MODULE is used to update the model's state from the previous timestep.
!     Currently, the forcing update for Simstrat is handled in the mod_forcing MODULE as a temporary solution,
!     and it will be modified later. Subsequently, the lake state will be updated in sequence.
!
! Original author: 
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024
! -------------------------------------------------------------------------------------
   USE strat_kinds, only: RK
   USE strat_inputfile, only: SimstratSimulationFactory
   USE strat_simdata, only: SimulationData
   USE strat_forcing
   USE strat_utilities
   USE strat_stability, only: StabilityModule
   USE strat_windshear
   USE strat_statevar
   USE strat_temp
   USE strat_solver
   USE strat_discretization
   USE strat_keps
   USE strat_turbulence
   USE strat_ice
   USE strat_transport
   USE strat_absorption
   USE strat_advection
   USE strat_lateral
   USE, intrinsic :: ieee_arithmetic
   USE MOD_Lake_Const, only : tfrz
   USE MOD_Lake_Utils, only: LakeIceFrc, snowage, snowdp2lev4lake
!================================================================================
!  ------------------------- input variables ---------------------------
   integer, intent(in)      :: &
      ipatch                  ,&! patch index
      nlake                   ,&! Maximum number of lake layer
      nsnow                   ,&! Maximum number of snow layer  !positive
      scwat                   ,&! surface category of water characteristics:
      nsoil                   ,&! Maximum number of soil layer
      zopt                    ,&! option for roughness length, 1: constant, 2: Subin et al. (2012)
      betaopt                 ,&! option for betaprime calculation, 1: constant, 2: equation
      fetchopt                ,&! option for fetch length, 1: constant, 2: equation
      etaopt                    ! option for Extinction coefficient calculation, 1: constant, 2: equation
       
   real(RK), intent(in)     :: &
      dtlak                   ,&! time step [s]
      xlat                    ,&! latitude in degrees
      xlon                    ,&! longitude in degrees
      dplak                   ,&! lake depth [m]
      hwref                   ,&! observational height of wind [m]
      htref                   ,&! observational height of temperature [m]
      hqref                   ,&! observational height of humidity [m]
      tmref                   ,&! temperature at the reference height [K] 
      usurf                   ,&! wind component in eastward direction [m/s]
      vsurf                   ,&! wind component in northward direction [m/s]
      qmref                   ,&! specific humidity at the reference height [kg/kg]
      arhos                   ,&! surface air density [kg m-3]
      sols                    ,&! atm vis direct beam solar rad onto srf [W m-2]
      soll                    ,&! atm nir direct beam solar rad onto srf [W m-2]
      solsd                   ,&! atm vis diffuse solar rad onto srf [W m-2]
      solld                   ,&! atm nir diffuse solar rad onto srf [W m-2]
      sabg                    ,&! solar absorbed by ground  [W/m2]
      lwdns                   ,&! atmospheric infrared (longwave) radiation [W/m2]
      psurf                   ,&! surface pressure [Pa]
      hpbl                    ,&! atmospheric boundary layer height [m]
      zcbcv                   ,&! convective boundary height [m]
      crain                   ,&! convective rainfall [kg/(m2 s)]
      csnow                   ,&! convective snowfall [kg/(m2 s)]
      lrain                   ,&! large scale rainfall [kg/(m2 s)]
      lsnow                     ! large scale snowfall [kg/(m2 s)]

!  ------------------------- inout variables ---------------------------
   integer, intent(inout)   :: &
      snlay                     ! number of snow layers

   real(RK), intent(inout)  :: &
      dzlak(nlake)            ,&! lake layer thickness (m)
      zlake(nlake)            ,&! lake layer depth (m)
      zilak(nlake+1)          ,&! interface level below a "z_soisno" level (m)
      lktmp(nlake)            ,&! lake temperature (K)
      uwatv(nlake)            ,&! Water velocity in x-direction [m/s]
      vwatv(nlake)            ,&! Water velocity in y-direction [m/s]
      lksal(nlake)            ,&! Salinity [‰]
      lkrho(nlake)            ,&! Water density [kg/m^3]
      icefr(nlake)            ,&! lake ice fraction [-]
      ziarea(nlake+1)         ,&! lake area (m2)
      tke(nlake+1)            ,&! Turbulent kinetic energy (TKE) [J/kg]
      eps(nlake+1)            ,&! TKE dissipation rate [W/kg]
      num(nlake+1)            ,&! Turbulent viscosity (momentum)
      nuh(nlake+1)              ! Turbulent diffusivity (heat)

   real(RK), intent(inout)  :: &
      dzssb(-nsnow+1:nsoil)   ,&! layer thickness (m)
      zssb(-nsnow+1:nsoil)    ,&! layer depth (m)
      zissb(-nsnow:nsoil+1)   ,&! interface level below a "z_soisno" level (m)
      t_ssb(-nsnow+1:nsoil)   ,&! soil + snow layer temperature [K]
      xwliq(-nsnow+1:nsoil)   ,&! liquid water (kg/m2)
      xwice(-nsnow+1:nsoil)     ! ice lens (kg/m2)

   real(RK), intent(inout)  :: &
      snwcv                   ,&! snow mass (kg/m2)
      icedp                   ,&! lake ice thickness (m)
      bicedp                  ,&! black ice depth (m)
      wicedp                  ,&! white ice depth (m)
      snwag                   ,&! non dimensional snow age [-]
      snwdp                   ,&! snow depth (m)
      rhosnw                  ,&! snow density (kg/m3)
      stke1                   ,&! top level eddy conductivity [W/m/K]
      etke                    ,&! Seiche energy [J]
      tskin                   ,&! ground surface temperature [k]
      tmice                   ,&! ice temperature [K]
      z0m                     ,&! roughness length over ground, momentum [m]
      z0h                     ,&! roughness length over ground, sensible heat [m]
      z0q                     ,&! roughness length over ground, latent heat [m]
      felak                   ,&! lake fetch (>0), otherwise USE emperical eq by dplak [m]
      gamma                   ,&! Mixing enhancement factorfor eddy diffusion coefficient
      etal                    ,&! extinction coefficient [1/m]
      btpri                   ,&! beta prime in Monin-Obukhov theory
      frlak                     ! lake fraction [-]

!  ------------------------- output variables ---------------------------
   real(RK), intent(out)  :: &
      fsena                 ,&! sensible heat from canopy height to atmosphere [W/m2]
      fevpa                 ,&! evapotranspiration from canopy height to atmosphere [mm/s]
      lfevpa                ,&! latent heat flux from canopy height to atmosphere [W/2]
      fevpg                 ,&! evaporation heat flux from ground [mm/s]
      fseng                 ,&! sensible heat flux from ground [W/m2]
      olrg                  ,&! outgoing long-wave radiation from ground+canopy
      fgrnd                 ,&! ground heat flux [W/m2]
      trad                  ,&! radiative temperature [K]
      qseva                 ,&! ground surface evaporation rate (mm h2o/s)
      qsdew                 ,&! ground surface dew formation (mm h2o /s) [+]
      qsubl                 ,&! sublimation rate from snow pack (mm h2o /s) [+]
      qfros                 ,&! surface dew added to snow pack (mm h2o /s) [+]
      taux                  ,&! wind stress: E-W [kg/m/s**2]
      tauy                  ,&! wind stress: N-S [kg/m/s**2]
      ustar                 ,&! u* in similarity theory [m/s]
      qstar                 ,&! q* in similarity theory [kg/kg]
      tstar                 ,&! t* in similarity theory [K]
      emis                  ,&! averaged bulk surface emissivity
      zol                   ,&! dimensionless height (z_soisno/L) used in Monin-Obukhov theory
      snwml                 ,&! snowmelt rate [mm/s]
      rib                   ,&! bulk Richardson number in surface layer
      ram                   ,&! aerodynamical resistance [s/m]
      rah                   ,&! thermal resistance [s/m]
      raw                   ,&! moisture resistance [s/m]
      wdm                   ,&! lowest level wind speed including the stablity effect [m/s]
      t2m                   ,&! 2 m height air temperature [K]
      q2m                   ,&! 2 m height air specific humidity
      u10m                  ,&! 10 m height wind speed in eastward direction [m/s]
      v10m                  ,&! 10 m height wind speed in northward direction [m/s]
      fm10m                 ,&! integral of profile function for momentum
      fq10m                 ,&! integral of profile function for moisture
      fh2m                  ,&! integral of profile function for heat
      fq2m                  ,&! integral of profile function for moisture
      fm                    ,&! integral of profile function for momentum
      fq                    ,&! integral of profile function for moisture
      fh                    ,&! integral of profile function for heat
      shfdt                   ! derivative of srf sen+rlat  heat flux wrt srf temp [W/m2/K]
        
   !+WMEJ urban_call is optional argument for SUBROUTINE laketem, it is not used in this version
   logical, optional, intent(in) :: urban_call   ! whether it is a urban CALL

!  ------------------------- local variables ---------------------------
   real(RK)               :: &
      dzlaktop              ,&! thickness of top lake layer [m]
      zlakebot              ,&! depth of top lake layer [m]
      lktmptop              ,&! temperature of top lake layer [K]
      icefrtop              ,&! lake mass fraction of top lake layer that is frozen
      dzssbtop              ,&! thickness of top lake layer [m]
      xwicetop              ,&! ice mass in top lake layer [kg/m2]
      xwliqtop              ,&! liquid water mass in top lake layer [kg/m2]
      tssbtop               ,&! temperature of top lake layer [K]
      rlat                  ,&! latitude in radians
      htvp                  ,&! latent heat flux parameter
      oldsnwdp              ,&! snow depth at previous time step [m]
      scvold                ,&! snow cover for previous time step [mm]
      w_old                 ,&! liquid water mass of the column at the previous time step (mm)
      frain                 ,&! rainfall onto ground including canopy runoff [kg/(m2 s)]
      fsnow                 ,&! snowfall onto ground including canopy runoff [kg/(m2 s)]
      t_precip              ,&! snowfall/rainfall temperature [kelvin]
      a                     ,&! fraction of snowfall that is new snow
      aa                    ,&! fraction of snowfall that is new snow
      eta                   ,&! extinction coefficient [1/m]
      iceml                 ,&! ice melt rate [mm/s]
      pg_snow               ,&! snow density [kg/m3]
      pg_rain               ,&! rain density [kg/m3]
      bifall                ,&! bulk density of newly fallen dry snow [kg/m3]
      fiold(-nsnow+1:nsoil) ,&! fraction of ice relative to the total water
      Qinp(1:4, nlake+1)    ,&! Horizontal inflow [m^3/s]
      wice_lake(nlake)      ,&! ice lens [kg/m2]
      wliq_lake(nlake)        ! liquid water [kg/m2]
       
   integer :: imelt(-nsnow+1:nsoil) ! flag for: melting=1, freezing=2, Nothing happended=0
   integer :: lb, j, i, snl
   
   ! Instantiate all modules
   ! note that some are pointers/targets for polymorphism reasons
   type(SimstratSimulationFactory) :: factory
   class(SimulationData), pointer :: simdata
   type(EulerIDiscretizationMFQ) :: euler_i_disc
   type(EulerIDiscretizationKEPS) :: euler_i_disc_keps
   type(ForcingModule) :: mod_forcing
   type(StabilityModule) :: mod_stability
   type(TempModelVar) :: mod_temperature
   type(UVModelVar) :: mod_u, mod_v
   type(KModelVar) :: mod_k
   type(EpsModelVar) :: mod_eps
   type(TransportModVar) :: mod_s
   type(TurbulenceModule) :: mod_turbulence
   type(IceModule) :: mod_ice
   type(AbsorptionModule) :: mod_absorption
   type(AdvectionModule) :: mod_advection
   class(GenericLateralModule), pointer :: mod_lateral
!================================================================================

      Qinp = 0.0_RK

      lb = snlay + 1
      dzlaktop = dzlak(1)
      zlakebot = zlake(nlake)
      lktmptop = lktmp(1)
      icefrtop = icefr(1)
      dzssbtop = dzssb(lb)
      xwicetop = xwice(lb)
      xwliqtop = xwliq(lb)
      tssbtop  = t_ssb(lb)
      oldsnwdp = snwdp
      scvold = snwcv
      lktmp = lktmp - tfrz
      tmice = tmice - tfrz

      ! Initialize all modules
      CALL reverse_variables ( &
            nlake   , dzlak   , zlake   , zilak   , lktmp   , uwatv   , vwatv   ,&
            lksal   , icefr   , tke     , eps     , num     , nuh     )
      
      CALL factory%initialize (  & 
            ! "inout" arguments
            ! ---------------------------
            simdata        , euler_i_disc   , euler_i_disc_keps, mod_forcing   ,&
            mod_stability  , mod_temperature, mod_u            , mod_v         ,&
            mod_k          , mod_eps        , mod_s            , mod_turbulence,&
            mod_ice        , mod_absorption , mod_advection    ,&
            ! "in" arguments
            ! ---------------------------
            nlake          , dtlak          , dplak            , ziarea        ,&
            zilak          , zlake          , dzlak            )

      CALL simdata%model%update ( &
            nlake          , snwdp          , bicedp           , wicedp        ,&
            rhosnw         , uwatv          , vwatv            , lktmp         ,&
            lksal          , tke            , eps              , etke          ,&
            lkrho          , num            , nuh              , tmice         ,& 
            Qinp           )

      simdata%model_param%Lat = xlat

      IF((simdata%model%snow_h + simdata%model%white_ice_h + simdata%model%black_ice_h) > 0) THEN ! Ice cover
         simdata%model%T_surf = simdata%model_param%Freez_Temp
      ELSE
         simdata%model%T_surf = simdata%model%T(simdata%grid%ubnd_vol) ! Free water
      ENDIF

      !- WMEJ simstrat surface flux
      CALL mod_forcing%update(simdata%model,&
            ! "in" arguments
            ! ---------------------------
            scwat     , snlay     , betaopt   , fetchopt  ,&
            zopt      , rlat      , dplak     , dtlak     ,&
            hwref     , htref     , hqref     , usurf     ,&
            vsurf     , tmref     , qmref     , arhos     ,&
            psurf     , sols      , soll      , solsd     ,&
            solld     , lwdns     , hpbl      , sabg      ,&
            stke1     , dzlaktop  , zlakebot  , lktmptop  ,&
            icefrtop  , dzssbtop  , xwicetop  , xwliqtop  ,&
            crain     , csnow     , lrain     , lsnow     ,&
            tssbtop   , zcbcv     , ipatch    ,&
            ! "inout" arguments
            ! -------------------
            tskin     , z0m       , z0h       , z0q       ,&
            felak     , btpri     ,&
            ! "out" arguments
            ! ---------------------------
            fsena     , fevpa     , lfevpa    , fseng     ,&
            fevpg     , olrg      , fgrnd     , trad      ,&
            taux      , tauy      , ustar     , qstar     ,&
            tstar     , emis      , zol       ,&
            rib       , ram       , rah       , raw       ,&
            wdm       , t2m       , q2m       , u10m      ,&
            v10m      , fm10m     , fq10m     , fh2m      ,&
            fq2m      , fm        , fq        , fh        ,&
            shfdt     )

      !- WMEJ set Extinction coefficient for each layer 
      IF (etaopt == 1) THEN
         eta = etal
      ELSE
         eta = 1.1925*max(dplak,1.)**(-0.424)
      ENDIF      

      ! Set absorption coefficient for each layer
      DO i = 1, nlake+1
         simdata%model%absorb(i) = eta
      ENDDO
      CALL reverse_in_place(simdata%model%absorb)

      ! ************************************
      ! ***** Compute next model state *****                                                                                     
      ! ************************************
      ! Update physics                                                                                                          
      CALL mod_stability%update(simdata%model)

      !---------------------------------------------------------------
      ! +WMEJ in the current version, we have disabled inflow and outflow
      ! +WMEJ This step will be implemented in the future
      ! +WMEJ 
      !---------------------------------------------------------------
      ! IF there is inflow/outflow DO advection part                                                                             
      IF (simdata%model_cfg%has_advection) THEN                        
         ! Treat inflow/outflow
         CALL mod_lateral%update(simdata%model)
         ! Set old lake level (before it is changed by advection MODULE)
         simdata%grid%lake_level_old = simdata%grid%z_face(simdata%grid%ubnd_fce)
         ! Update lake advection using the inflow/outflow data
         CALL mod_advection%update(simdata%model)
         ! Update lake level
         simdata%grid%lake_level = simdata%grid%z_face(simdata%grid%ubnd_fce)
      ENDIF
      !---------------------------------------------------------------

      ! Update Coriolis                                                                                                          
      CALL mod_forcing%update_coriolis(simdata%model)

      ! Update and solve U and V - terms 
      CALL mod_u%update(simdata%model, simdata%model_param)
      CALL mod_v%update(simdata%model, simdata%model_param)

      ! Update and solve T - terms          
      CALL mod_temperature%update(simdata%model, simdata%model_param)

      ! Update and solve transportation terms (here: Salinity S only)                                                                              
      CALL mod_S%update(simdata%model, simdata%model_param) 

      ! update turbulence states                                                                                        
      CALL mod_turbulence%update(simdata%model, simdata%model_param)

      ! Solve k & eps        
      CALL mod_k%update(simdata%model, simdata%model_param)
      CALL mod_eps%update(simdata%model, simdata%model_param)
      
      ! Update ice                                                                                                           
      IF (simdata%model_cfg%ice_model == 1) THEN
         CALL mod_ice%update(simdata%model, simdata%model_param)
      ENDIF

      IF((simdata%model%snow_h + simdata%model%white_ice_h + simdata%model%black_ice_h) > 0) THEN ! Ice cover
         simdata%model%T_surf = simdata%model_param%Freez_Temp
      ELSE
         simdata%model%T_surf = simdata%model%T(simdata%grid%ubnd_vol) ! Free water
      ENDIF

      ! **********************************
      ! ***** Write output variables *****                                      
      ! **********************************
      lktmp = simdata%model%T + tfrz
      uwatv = simdata%model%U  
      vwatv = simdata%model%V  
      lksal = simdata%model%S  
      tke = simdata%model%K  
      eps = simdata%model%eps  
      etke = simdata%model%E_Seiche  
      num = simdata%model%num  
      nuh = simdata%model%nuh  
      lkrho = simdata%model%rho  
      tmice = simdata%model%ice_temp + tfrz
      rhosnw = simdata%model%snow_dens  
      snwdp = simdata%model%snow_h  
      wicedp = simdata%model%white_ice_h  
      bicedp = simdata%model%black_ice_h  
      icedp = simdata%model%total_ice_h  
      snwml = simdata%model%snwml
      iceml = simdata%model%iceml

      CALL reverse_variables (    nlake   , dzlak   , zlake   , zilak   , lktmp   , uwatv   , vwatv   ,&
                                 lksal   , icefr   , tke     , eps     , num     , nuh     )

      tskin = simdata%model%T_surf + tfrz

      
      ! **********************************
      ! **** Cal additional variables ****
      ! **********************************
      CALL LakeIceFrc( &
            ! "in" arguments
            ! -------------------
            nlake     , dzlak     , zilak     , icedp     ,&
            ! "inout" arguments
            ! -------------------
            wice_lake , wliq_lake , icefr     )

      CALL snowdp2lev4lake( &
            ! "in" arguments
            ! ---------------------------
            nsnow = nsnow   , snowdp = snwdp   ,&
            ! "inout" arguments
            ! ---------------------------
            dzsno = dzssb(:0), zsno  = zssb(:0),& 
            zisno = zissb(:0), snlay = snlay    )
      lb = snlay + 1
      xwliq = 0.
      xwice = 0.
      t_ssb(-nsnow+1:0) = lktmp(1) + tfrz
      t_ssb(snlay+1:0) = simdata%model_param%snow_temp + tfrz

      ! +WMEJ update snow pack, following the cssplake
      qseva = 0.
      qsubl = 0.
      qsdew = 0.
      qfros = 0.
      IF (fevpg >= 0.) THEN
      ! sublimation. DO not allow for more sublimation than there is snow
      ! after melt. remaining surface evaporation used for infiltration
         ! qsubl = min( fevpg, snwcv/dtlak-snwml )
         IF(lb < 0)THEN
            qseva = min(xwliq(lb)/dtlak, fevpg)
            qsubl = fevpg - qseva
         ELSE
            qseva = min((1.-icefr(1))*1000.*dzlak(1)/dtlak, fevpg)
            qsubl = fevpg - qseva
         ENDIF
      ELSE
         IF (tskin < tfrz-0.1) THEN
            qfros = abs(fevpg)
         ELSE
            qsdew = abs(fevpg)
         ENDIF
      ENDIF
      snwcv = snwcv + (pg_snow-snwml-qsubl+qfros)*dtlak
      snwcv = max( snwcv, 0. )
      
      ! no snow IF lake unfrozen
      IF (tskin > tfrz) snwcv = 0.
      CALL snowage (dtlak,tskin,snwcv,scvold,snwag)

   END SUBROUTINE Lake_Simstrat



   SUBROUTINE reverse_variables ( &
               nlake   ,&
               dzlak   , zlake   , zilak   , lktmp   ,&
               uwatv   , vwatv   , lksal   , icefr   ,&
               tke     , eps     , num     , nuh     )

! ---------------------------------- code history -------------------------------------
! DESCRIPTION:
!     Reverse arrays to fit other models
!
! Original author : 
!     Shulei Zhang(zsl) & Omarjan Obulkasim(WMEJ), 09/2024 
!
!--------------------------------------------------------------------------------------
   USE strat_kinds, only: RK
   USE strat_utilities, only: reverse_in_place
!================================================================================
   integer, intent(in)      :: &
      nlake                     ! Maximum number of lake layer

   real(RK), intent(inout)  :: &
      dzlak(nlake)            ,&! lake layer thickness [m]
      zlake(nlake)            ,&! lake layer node depth [m]
      lktmp(nlake)            ,&! lake temperature (K)
      uwatv(nlake)            ,&! Water velocity in x-direction [m/s]
      vwatv(nlake)            ,&! Water velocity in y-direction [m/s]
      lksal(nlake)            ,&! Salinity [‰]
      icefr(nlake)              ! lake ice fraction [-]

   real(RK), intent(inout)  :: &
      tke(nlake+1)            ,&! Turbulent kinetic energy (TKE) [J/kg]
      eps(nlake+1)            ,&! TKE dissipation rate [W/kg]
      num(nlake+1)            ,&! Turbulent viscosity (momentum)
      nuh(nlake+1)            ,&! Turbulent diffusivity (heat)
      zilak(nlake+1)           ! lake layer interface level [m]
!================================================================================
      CALL reverse_in_place(lktmp)
      CALL reverse_in_place(zilak)
      CALL reverse_in_place(zlake)
      CALL reverse_in_place(dzlak)
      CALL reverse_in_place(uwatv)
      CALL reverse_in_place(vwatv)
      CALL reverse_in_place(lksal)
      CALL reverse_in_place(tke)
      CALL reverse_in_place(eps)
      CALL reverse_in_place(num)
      CALL reverse_in_place(nuh)

   END SUBROUTINE reverse_variables


END MODULE  MOD_Lake_Simstrat


