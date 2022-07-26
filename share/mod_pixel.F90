#include <define.h>

MODULE mod_pixel

   USE precision
   IMPLICIT NONE

   ! ---- data types ----
   TYPE :: pixel_type

      REAL(r8) :: edges  ! southern edge (degrees)
      REAL(r8) :: edgen  ! northern edge (degrees)
      REAL(r8) :: edgew  ! western  edge (degrees)
      REAL(r8) :: edgee  ! eastern  edge (degrees)
      
      INTEGER :: nlon, nlat
      REAL(r8), allocatable :: lat_s (:)
      REAL(r8), allocatable :: lat_n (:)
      REAL(r8), allocatable :: lon_w (:)
      REAL(r8), allocatable :: lon_e (:)

   CONTAINS 
      procedure, PUBLIC :: set_edges   => pixel_set_edges

      procedure, PRIVATE :: assimilate_latlon => pixel_assimilate_latlon 
      procedure, PUBLIC  :: assimilate_gblock => pixel_assimilate_gblock 
      procedure, PUBLIC  :: assimilate_grid   => pixel_assimilate_grid 

      procedure, PUBLIC :: map_to_grid => pixel_map_to_grid

      procedure, PUBLIC :: save_to_file   => pixel_save_to_file
      procedure, PUBLIC :: load_from_file => pixel_load_from_file
      
      final :: pixel_free_mem

   END TYPE pixel_type
   
   ! ---- Instance ----
   TYPE(pixel_type) :: pixel

CONTAINS
   
   ! --------------------------------
   SUBROUTINE pixel_set_edges (this, & 
         edges_in, edgen_in, edgew_in, edgee_in)
      
      USE precision
      USE spmd_task
      USE mod_utils
      IMPLICIT NONE

      class(pixel_type) :: this

      REAL(r8), intent(in) :: edges_in, edgen_in
      REAL(r8), intent(in) :: edgew_in, edgee_in

      this%nlon = 1
      this%nlat = 1

      this%edges = edges_in
      this%edgen = edgen_in
      this%edgew = edgew_in
      this%edgee = edgee_in
      
      CALL normalize_longitude (this%edgew)
      CALL normalize_longitude (this%edgee)

      allocate (this%lat_s (1))
      allocate (this%lat_n (1))
      allocate (this%lon_w (1))
      allocate (this%lon_e (1))

      this%lat_s(1) = this%edges
      this%lat_n(1) = this%edgen
      this%lon_w(1) = this%edgew
      this%lon_e(1) = this%edgee

      IF (p_is_master) THEN
         write(*,*) 'Region information:'
         write(*,'(A28,F10.4,A1,F10.4,A1,F10.4,A1,F10.4,A1)') ' (south,north,west,east) = (', &
            this%edges, ',', this%edgen, ',', this%edgew, ',', this%edgee, ')'
         write(*,*)
      ENDIF

   END SUBROUTINE pixel_set_edges

   ! --------------------------------
   SUBROUTINE pixel_assimilate_latlon (this, &
         nlat, lat_s, lat_n, nlon, lon_w, lon_e)

      USE precision
      USE mod_utils
      IMPLICIT NONE
      class(pixel_type) :: this

      INTEGER,  intent(in) :: nlat
      REAL(r8), intent(in) :: lat_s(nlat), lat_n(nlat)
      INTEGER,  intent(in) :: nlon
      REAL(r8), intent(in) :: lon_w(nlon), lon_e(nlon)

      ! Local variables
      INTEGER :: ny, yinc
      INTEGER :: iy1, iy2, ys2, yn2
      REAL(r8), allocatable :: ytmp(:)

      INTEGER :: nx
      INTEGER :: ix1, ix2, xw2, xe2
      REAL(r8), allocatable :: xtmp(:)

      IF (lat_s(1) <= lat_s(nlat)) THEN
         yinc = 1
      ELSE
         yinc = -1
      ENDIF

      allocate (ytmp (this%nlat+nlat+2))

      ny = 0
      DO iy1 = 1, this%nlat   
         ys2 = find_nearest_south (this%lat_s(iy1), nlat, lat_s) 
         yn2 = find_nearest_north (this%lat_n(iy1), nlat, lat_n) 

         ny = ny + 1
         ytmp(ny) = this%lat_s(iy1)

         DO iy2 = ys2, yn2, yinc
            IF (lat_s(iy2) > this%lat_s(iy1)) THEN
               ny = ny + 1
               ytmp(ny) = lat_s(iy2)
            ENDIF
         ENDDO
      ENDDO

      ny = ny + 1
      ytmp(ny) = this%lat_n(this%nlat)

      deallocate (this%lat_s)
      deallocate (this%lat_n)

      this%nlat = ny - 1
      allocate (this%lat_s (this%nlat))
      allocate (this%lat_n (this%nlat))

      this%lat_s = ytmp(1:ny-1)
      this%lat_n = ytmp(2:ny)

      deallocate (ytmp)
      

      allocate (xtmp (this%nlon+nlon+2))

      nx = 0
      DO ix1 = 1, this%nlon
         xw2 = find_nearest_west (this%lon_w(ix1), nlon, lon_w) 
         xe2 = find_nearest_east (this%lon_e(ix1), nlon, lon_e) 

         nx = nx + 1
         xtmp(nx) = this%lon_w(ix1)

         IF (xw2 /= xe2) THEN
            ix2 = mod(xw2,nlon) + 1
            DO while (.true.) 
               nx = nx + 1
               xtmp(nx) = lon_w(ix2)

               IF (ix2 == xe2)  exit
               ix2 = mod(ix2,nlon) + 1
            ENDDO
         ENDIF
      ENDDO

      nx = nx + 1
      xtmp(nx) = this%lon_e(this%nlon)

      deallocate (this%lon_w)
      deallocate (this%lon_e)

      this%nlon = nx - 1
      allocate (this%lon_w (this%nlon))
      allocate (this%lon_e (this%nlon))

      this%lon_w = xtmp(1:nx-1)
      this%lon_e = xtmp(2:nx)

      deallocate (xtmp)

   END SUBROUTINE pixel_assimilate_latlon

   ! --------------------------------
   SUBROUTINE pixel_assimilate_gblock (this)

      USE mod_block, only : gblock
      IMPLICIT NONE
      class(pixel_type) :: this

      CALL this%assimilate_latlon ( &
         gblock%nyblk, gblock%lat_s, gblock%lat_n, &
         gblock%nxblk, gblock%lon_w, gblock%lon_e)
      
   END SUBROUTINE pixel_assimilate_gblock

   ! --------------------------------
   SUBROUTINE pixel_assimilate_grid (this, grid)

      USE mod_grid
      IMPLICIT NONE
      class(pixel_type) :: this

      TYPE(grid_type), intent(in) :: grid

      CALL this%assimilate_latlon ( &
         grid%nlat, grid%lat_s, grid%lat_n, &
         grid%nlon, grid%lon_w, grid%lon_e)
      
   END SUBROUTINE pixel_assimilate_grid

   ! --------------------------------
   SUBROUTINE pixel_map_to_grid (this, grd)

      USE mod_grid
      USE mod_utils
      IMPLICIT NONE
      class(pixel_type) :: this

      TYPE(grid_type), intent(inout) :: grd 

      ! Local variables
      INTEGER :: iy1, iy2, ix1, ix2

      IF (allocated (grd%xgrd))  deallocate (grd%xgrd)
      IF (allocated (grd%ygrd))  deallocate (grd%ygrd)

      allocate (grd%ygrd (this%nlat))

      iy1 = 1
      DO while (.true.)
         iy2 = find_nearest_south (this%lat_s(iy1), grd%nlat, grd%lat_s) 
         IF (this%lat_n(iy1) <= grd%lat_n(iy2)) THEN
            DO while (this%lat_n(iy1) <= grd%lat_n(iy2))
               grd%ygrd(iy1) = iy2
               iy1 = iy1 + 1
               IF (iy1 > this%nlat) exit
            ENDDO
         ELSE
            grd%ygrd(iy1) = iy2
            iy1 = iy1 + 1
         ENDIF

         IF (iy1 > this%nlat) exit
      ENDDO

      allocate (grd%xgrd (this%nlon))
      
      ix1 = 1
      DO while (.true.)
         ix2 = find_nearest_west (this%lon_w(ix1), grd%nlon, grd%lon_w) 
         DO while (lon_between_ceil(this%lon_e(ix1), grd%lon_w(ix2), grd%lon_e(ix2)))
            grd%xgrd(ix1) = ix2
            ix1 = ix1 + 1
            IF (ix1 > this%nlon) exit
         ENDDO
         IF (ix1 > this%nlon) exit
      ENDDO

   END SUBROUTINE pixel_map_to_grid 

   ! --------------------------------
   SUBROUTINE pixel_save_to_file (this, dir_landdata)

      USE spmd_task
      USE ncio_serial
      IMPLICIT NONE
      class(pixel_type) :: this

      CHARACTER(len=*), intent(in) :: dir_landdata 

      ! Local variables
      CHARACTER(len=256) :: filename

      IF (p_is_master) THEN
         
         filename = trim(dir_landdata) // '/pixel.nc'

         CALL ncio_create_file (filename)
         
         CALL ncio_write_serial (filename, 'edges', this%edges)
         CALL ncio_write_serial (filename, 'edgen', this%edgen)
         CALL ncio_write_serial (filename, 'edgew', this%edgew)
         CALL ncio_write_serial (filename, 'edgee', this%edgee)

         CALL ncio_define_dimension (filename, 'latitude',  this%nlat)
         CALL ncio_define_dimension (filename, 'longitude', this%nlon)

         CALL ncio_write_serial (filename, 'lat_s', this%lat_s, 'latitude' )
         CALL ncio_write_serial (filename, 'lat_n', this%lat_n, 'latitude' )
         CALL ncio_write_serial (filename, 'lon_w', this%lon_w, 'longitude')
         CALL ncio_write_serial (filename, 'lon_e', this%lon_e, 'longitude')

      ENDIF

   END SUBROUTINE pixel_save_to_file

   ! --------------------------------
   SUBROUTINE pixel_load_from_file (this, dir_landdata)

      USE spmd_task
      USE ncio_serial 
      IMPLICIT NONE

      class(pixel_type) :: this

      CHARACTER(len=*), intent(in) :: dir_landdata
      ! Local variables
      CHARACTER(len=256) :: filename

      filename = trim(dir_landdata) // '/pixel.nc'

      CALL ncio_read_bcast_serial (filename, 'edges', this%edges)
      CALL ncio_read_bcast_serial (filename, 'edgen', this%edgen)
      CALL ncio_read_bcast_serial (filename, 'edgew', this%edgew)
      CALL ncio_read_bcast_serial (filename, 'edgee', this%edgee)

      CALL ncio_read_bcast_serial (filename, 'lat_s', this%lat_s)
      CALL ncio_read_bcast_serial (filename, 'lat_n', this%lat_n)
      CALL ncio_read_bcast_serial (filename, 'lon_w', this%lon_w)
      CALL ncio_read_bcast_serial (filename, 'lon_e', this%lon_e)

      this%nlon = size(this%lon_w)
      this%nlat = size(this%lat_s)

   END SUBROUTINE pixel_load_from_file

   ! --------------------------------
   SUBROUTINE pixel_free_mem (this)
      
      IMPLICIT NONE
      TYPE (pixel_type) :: this

      IF (allocated(this%lat_s))  deallocate(this%lat_s)
      IF (allocated(this%lat_n))  deallocate(this%lat_n)
      IF (allocated(this%lon_w))  deallocate(this%lon_w)
      IF (allocated(this%lon_e))  deallocate(this%lon_e)

   END SUBROUTINE pixel_free_mem 

END MODULE mod_pixel
