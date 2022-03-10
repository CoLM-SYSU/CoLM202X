#include <define.h>

MODULE mod_grid

   USE precision
   IMPLICIT NONE

   ! ---- data types ----
   TYPE :: grid_type

      INTEGER :: nlat
      INTEGER :: nlon

      ! Latitude direction. (yinc = 1) means south to north.  
      INTEGER :: yinc

      ! Coordinates.
      REAL(r8), allocatable :: lat_s (:)
      REAL(r8), allocatable :: lat_n (:)
      REAL(r8), allocatable :: lon_w (:)
      REAL(r8), allocatable :: lon_e (:)
      
      ! Blocks.
      INTEGER, allocatable :: xdsp(:), ydsp(:)
      INTEGER, allocatable :: xcnt(:), ycnt(:)
      INTEGER, allocatable :: xblk(:), yblk(:)
      INTEGER, allocatable :: xloc(:), yloc(:)

      ! Mapping to pixels.
      INTEGER, allocatable :: xgrd(:), ygrd(:)

      ! Grid info.
      REAL(r8), allocatable :: rlon(:)
      REAL(r8), allocatable :: rlat(:)

   CONTAINS
      procedure, PUBLIC :: define_by_name   => grid_define_by_name  
      procedure, PUBLIC :: define_by_ndims  => grid_define_by_ndims  
      procedure, PUBLIC :: define_by_center => grid_define_by_center
      procedure, PUBLIC :: define_from_file => grid_define_from_file

      procedure, PUBLIC :: set_rlon => grid_set_rlon
      procedure, PUBLIC :: set_rlat => grid_set_rlat

      procedure, PRIVATE :: init => grid_init
      procedure, PRIVATE :: normalize  => grid_normalize
      procedure, PRIVATE :: set_blocks => grid_set_blocks
      
      final :: grid_free_mem

   END TYPE grid_type
   
   ! ---- data types ----
   TYPE :: grid_list_type
      INTEGER :: ng
      INTEGER, allocatable :: ilat(:)
      INTEGER, allocatable :: ilon(:)
   END TYPE grid_list_type

CONTAINS

   ! --------------------------------
   SUBROUTINE grid_init (this, nlon, nlat)
      
      IMPLICIT NONE
      class (grid_type) :: this

      INTEGER, intent(in) :: nlon
      INTEGER, intent(in) :: nlat

      this%nlat = nlat
      this%nlon = nlon

      allocate (this%lat_s (nlat))
      allocate (this%lat_n (nlat))
      allocate (this%lon_w (nlon))
      allocate (this%lon_e (nlon))

   END SUBROUTINE grid_init
   
   ! --------------------------------
   SUBROUTINE grid_define_by_name (this, gridname)

      IMPLICIT NONE
      class (grid_type) :: this

      CHARACTER(len=*), intent(in) :: gridname

      ! Local variables
      INTEGER  :: nlat, nlon, ilat, ilon
      REAL(r8) :: del_lat, del_lon
      REAL(r8), allocatable :: lon(:), lat(:)

      IF (trim(gridname) == 'merit_90m') THEN 

         nlat = 180*60*20
         nlon = 360*60*20

         this%nlat = nlat
         this%nlon = nlon

         CALL this%init (this%nlon, this%nlat)

         del_lat = 180.0 / nlat
         DO ilat = 1, this%nlat
            this%lat_s(ilat) = 90.0 - del_lat * ilat - del_lat/2.0 
            this%lat_n(ilat) = 90.0 - del_lat * (ilat-1) - del_lat/2.0
         ENDDO

         del_lon = 360.0 / nlon
         DO ilon = 1, this%nlon
            this%lon_w(ilon) = -180.0 + del_lon * (ilon-1) - del_lon/2.0
            this%lon_e(ilon) = -180.0 + del_lon * ilon - del_lon/2.0
         ENDDO

         CALL this%normalize  ()
         CALL this%set_blocks ()

      ENDIF

      IF (trim(gridname) == 'colm_1km') THEN

         CALL this%define_by_ndims (43200,21600)

      ENDIF

      IF (trim(gridname) == 'colm_500m') THEN

         CALL this%define_by_ndims (86400,43200)

      ENDIF

      IF (trim(gridname) == 'ERA5LAND') THEN

         nlat = 1801
         nlon = 3600

         allocate (lon(nlon))
         allocate (lat(nlat))

         del_lat = 180.0 / (nlat-1)
         DO ilat = 1, nlat
            lat(ilat) = 90.0 - del_lat * (ilat-1)
         ENDDO

         del_lon = 360.0 / nlon
         DO ilon = 1, nlon
            lon(ilon) = (ilon-1) * del_lon
         ENDDO

         call this%define_by_center (lat, lon)

      ENDIF

      IF (trim(gridname) == 'ERA5') THEN
         nlat = 721
         nlon = 1440

         allocate (lon(nlon))
         allocate (lat(nlat))

         del_lat = 180.0 / (nlat-1)
         DO ilat = 1, nlat
            lat(ilat) = 90.0 - del_lat * (ilat-1)
         ENDDO

         del_lon = 360.0 / nlon
         DO ilon = 1, nlon
            lon(ilon) = (ilon-1) * del_lon
         ENDDO

         call this%define_by_center (lat, lon)

      ENDIF

   END SUBROUTINE grid_define_by_name

   !-----------------------------------------------------
   SUBROUTINE grid_define_by_ndims (this, lon_points, lat_points)

      IMPLICIT NONE
      class (grid_type) :: this

      INTEGER, intent(in) :: lon_points
      INTEGER, intent(in) :: lat_points

      ! Local variables
      INTEGER  :: ilat, ilon
      REAL(r8) :: del_lat, del_lon

      this%nlat = lat_points
      this%nlon = lon_points

      CALL this%init (this%nlon, this%nlat)

      del_lat = 180.0 / lat_points
      DO ilat = 1, this%nlat
         this%lat_s(ilat) = 90.0 - del_lat * ilat
         this%lat_n(ilat) = 90.0 - del_lat * (ilat-1)
      ENDDO

      del_lon = 360.0 / lon_points
      DO ilon = 1, this%nlon
         this%lon_w(ilon) = -180.0 + del_lon * (ilon-1)
         this%lon_e(ilon) = -180.0 + del_lon * ilon
      ENDDO

      CALL this%normalize  ()
      CALL this%set_blocks ()

   END SUBROUTINE grid_define_by_ndims

   !---------------------------------------------
   SUBROUTINE grid_define_by_center (this, lat_in, lon_in)
      
      USE precision
      USE mod_utils
      IMPLICIT NONE
      class (grid_type) :: this

      REAL(r8), intent(in) :: lat_in(:), lon_in(:)

      ! Local variables
      INTEGER :: ilat, ilon, ilone, ilonw
      REAL(r8), allocatable :: lon_in_n(:)

      this%nlat = size(lat_in)
      this%nlon = size(lon_in)

      CALL this%init (this%nlon, this%nlat)

      IF (lat_in(1) > lat_in(this%nlat)) THEN
         this%yinc = -1
      ELSE
         this%yinc = 1
      ENDIF

      DO ilat = 1, this%nlat
         IF (this%yinc == 1) THEN
            IF (ilat < this%nlat) THEN
               this%lat_n(ilat) = (lat_in(ilat) + lat_in(ilat+1)) * 0.5
            ELSE
               this%lat_n(ilat) = 90.0
            ENDIF

            IF (ilat > 1) THEN
               this%lat_s(ilat) = (lat_in(ilat-1) + lat_in(ilat)) * 0.5
            ELSE
               this%lat_s(ilat) = -90.0
            ENDIF
         ELSEIF (this%yinc == -1) THEN
            IF (ilat > 1) THEN
               this%lat_n(ilat) = (lat_in(ilat-1) + lat_in(ilat)) * 0.5
            ELSE
               this%lat_n(ilat) = 90.0
            ENDIF

            IF (ilat < this%nlat) THEN
               this%lat_s(ilat) = (lat_in(ilat) + lat_in(ilat+1)) * 0.5
            ELSE
               this%lat_s(ilat) = -90.0
            ENDIF
         ENDIF
      ENDDO

      allocate (lon_in_n (size(lon_in)))
      
      lon_in_n = lon_in
      DO ilon = 1, size(lon_in_n)
         CALL normalize_longitude (lon_in_n(ilon))
      ENDDO

      DO ilon = 1, this%nlon
         ilone = mod(ilon,this%nlon) + 1
         IF (lon_in_n(ilon) > lon_in_n(ilone)) THEN
            this%lon_e(ilon) = (lon_in_n(ilon) + lon_in_n(ilone) + 360.0) * 0.5
         ELSE
            this%lon_e(ilon) = (lon_in_n(ilon) + lon_in_n(ilone)) * 0.5
         ENDIF

         ilonw = ilon - 1
         IF (ilonw == 0) ilonw = this%nlon
         IF (lon_in_n(ilonw) > lon_in_n(ilon)) THEN
            this%lon_w(ilon) = (lon_in_n(ilonw) + lon_in_n(ilon) + 360.0) * 0.5
         ELSE
            this%lon_w(ilon) = (lon_in_n(ilonw) + lon_in_n(ilon)) * 0.5
         ENDIF
      ENDDO

      deallocate (lon_in_n)

      CALL this%normalize  ()
      CALL this%set_blocks ()

   END SUBROUTINE grid_define_by_center

   !-----------------------------------------------------
   SUBROUTINE grid_define_from_file (this, filename)
      
      USE ncio_serial
      IMPLICIT NONE
      class (grid_type) :: this

      CHARACTER(len=*), intent(in) :: filename

      CALL ncio_read_bcast_serial (filename, 'lat_s', this%lat_s)
      CALL ncio_read_bcast_serial (filename, 'lat_n', this%lat_n)
      CALL ncio_read_bcast_serial (filename, 'lon_w', this%lon_w)
      CALL ncio_read_bcast_serial (filename, 'lon_e', this%lon_e)

      this%nlat = size(this%lat_s)
      this%nlon = size(this%lon_w)

      CALL this%normalize  ()
      CALL this%set_blocks ()
      
   END SUBROUTINE grid_define_from_file

   !-----------------------------------------------------
   SUBROUTINE grid_normalize (this)

      USE mod_utils
      IMPLICIT NONE
      class(grid_type) :: this

      ! Local variable
      INTEGER :: ilon, ilat

      DO ilon = 1, this%nlon
         CALL normalize_longitude (this%lon_w(ilon))
         CALL normalize_longitude (this%lon_e(ilon))
      ENDDO

      DO ilat = 1, this%nlat
         this%lat_s(ilat) = max(-90.0, min(90.0, this%lat_s(ilat)))
         this%lat_n(ilat) = max(-90.0, min(90.0, this%lat_n(ilat)))
      ENDDO

      IF (this%lat_s(1) <= this%lat_s(this%nlat)) THEN
         this%yinc = 1
      ELSE
         this%yinc = -1
      ENDIF

   END SUBROUTINE grid_normalize

   !-----------------------------------------------------
   SUBROUTINE grid_set_blocks (this)
      
      USE mod_namelist
      USE mod_block
      USE mod_utils
      IMPLICIT NONE
     
      class (grid_type) :: this
      
      ! Local variables
      INTEGER  :: ilat, ilon, iblk, jblk, ilon_e
      INTEGER  :: ilat_north, ilat_south, ilon_east, ilon_west
      REAL(r8) :: edges, edgen, edgew, edgee 
   
      allocate (this%xcnt (gblock%nxblk))
      allocate (this%xdsp (gblock%nxblk))
      allocate (this%ycnt (gblock%nyblk))
      allocate (this%ydsp (gblock%nyblk))
      
      allocate (this%xblk (this%nlon))
      allocate (this%yblk (this%nlat))
      
      allocate (this%xloc (this%nlon))
      allocate (this%yloc (this%nlat))

      edges = DEF_domain%edges
      edgen = DEF_domain%edgen
      edgew = DEF_domain%edgew
      edgee = DEF_domain%edgee
         
      CALL normalize_longitude (edgew)
      CALL normalize_longitude (edgee)

      IF (this%yinc == 1) THEN
         
         this%ycnt(:) = 0
         this%yblk(:) = 0

         jblk = find_nearest_south (edges, gblock%nyblk, gblock%lat_s)
         ilat = find_nearest_south (edges, this%nlat, this%lat_s)

         this%ydsp(jblk) = ilat - 1

         DO while (ilat <= this%nlat)
            IF (this%lat_s(ilat) < edgen) THEN
               IF (this%lat_s(ilat) < gblock%lat_n(jblk)) THEN

                  this%ycnt(jblk) = this%ycnt(jblk) + 1

                  this%yblk(ilat) = jblk
                  this%yloc(ilat) = this%ycnt(jblk)

                  ilat = ilat + 1
               ELSE
                  jblk = jblk + 1
                  IF (jblk <= gblock%nyblk) THEN
                     this%ydsp(jblk) = ilat - 1
                  ELSE
                     exit
                  ENDIF
               ENDIF
            ELSE
               exit
            ENDIF
         ENDDO

      ELSE

         this%ycnt(:) = 0
         this%yblk(:) = 0

         jblk = find_nearest_north (edgen, gblock%nyblk, gblock%lat_n)
         ilat = find_nearest_north (edgen, this%nlat, this%lat_n)

         this%ydsp(jblk) = ilat - 1

         DO while (ilat <= this%nlat)
            IF (this%lat_n(ilat) > edges) THEN
               IF (this%lat_n(ilat) > gblock%lat_s(jblk)) THEN

                  this%ycnt(jblk) = this%ycnt(jblk) + 1

                  this%yblk(ilat) = jblk
                  this%yloc(ilat) = this%ycnt(jblk)

                  ilat = ilat + 1
               ELSE
                  jblk = jblk - 1
                  IF (jblk >= 1) THEN
                     this%ydsp(jblk) = ilat - 1
                  ELSE
                     exit
                  ENDIF
               ENDIF
            ELSE
               exit
            ENDIF
         ENDDO

      ENDIF


      this%xcnt(:) = 0
      this%xblk(:) = 0

      iblk = find_nearest_west (edgew, gblock%nxblk, gblock%lon_w)
      ilon = find_nearest_west (edgew, this%nlon, this%lon_w)

      this%xdsp(iblk) = ilon - 1
      this%xcnt(iblk) = 1
      this%xblk(ilon) = iblk
      this%xloc(ilon) = 1

      ilon_e = ilon - 1
      IF (ilon_e == 0) ilon_e = this%nlon

      ilon = mod(ilon,this%nlon) + 1
      DO while (.true.)
         IF (lon_between_floor(this%lon_w(ilon), edgew, edgee)) THEN
            IF (lon_between_floor(this%lon_w(ilon), gblock%lon_w(iblk), gblock%lon_e(iblk))) THEN

               this%xcnt(iblk) = this%xcnt(iblk) + 1

               this%xblk(ilon) = iblk
               this%xloc(ilon) = this%xcnt(iblk)

               IF (ilon /= ilon_e) THEN
                  ilon = mod(ilon,this%nlon) + 1
               ELSE
                  exit
               ENDIF
            ELSE
               iblk = mod(iblk,gblock%nxblk) + 1
               IF (this%xcnt(iblk) == 0) THEN
                  this%xdsp(iblk) = ilon - 1
               ELSE
                  ilon_e = this%xdsp(iblk) + this%xcnt(iblk)
                  IF (ilon_e > this%nlon) ilon_e = ilon_e - this%nlon

                  this%xdsp(iblk) = ilon - 1
                  this%xcnt(iblk) = 0
                  DO while (.true.)
                     this%xcnt(iblk) = this%xcnt(iblk) + 1
                     this%xblk(ilon) = iblk
                     this%xloc(ilon) = this%xcnt(iblk)

                     IF (ilon /= ilon_e) THEN
                        ilon = mod(ilon,this%nlon) + 1
                     ELSE
                        exit
                     ENDIF
                  ENDDO
                  
                  exit
               ENDIF
            ENDIF
         ELSE
            exit
         ENDIF
      ENDDO

   END SUBROUTINE grid_set_blocks
   
   !-----------
   SUBROUTINE grid_set_rlon (this)

      USE precision
      USE mod_utils
      USE MathConstants, only : pi
      IMPLICIT NONE

      class (grid_type) :: this

      ! Local variables
      INTEGER  :: ix
      REAL(r8) :: lon

      IF (.not. allocated(this%rlon)) THEN
         allocate (this%rlon(this%nlon))
      ENDIF

      DO ix = 1, this%nlon
         IF (this%lon_w(ix) <= this%lon_e(ix)) THEN
            lon = (this%lon_w(ix) + this%lon_e(ix)) * 0.5
         ELSE
            lon = (this%lon_w(ix) + this%lon_e(ix)) * 0.5 + 180.0
         ENDIF

         CALL normalize_longitude (lon)

         this%rlon(ix) = lon / 180.0_r8 * pi
      ENDDO

   END SUBROUTINE grid_set_rlon
   
   !-----------
   SUBROUTINE grid_set_rlat (this)

      USE precision
      USE mod_utils
      USE MathConstants, only : pi
      IMPLICIT NONE

      class (grid_type) :: this
      
      ! Local variables
      INTEGER :: iy

      IF (.not. allocated(this%rlat)) THEN
         allocate (this%rlat(this%nlat))
      ENDIF

      DO iy = 1, this%nlat
         this%rlat(iy) = &
            (this%lat_s(iy) + this%lat_n(iy)) * 0.5 / 180.0_r8 * pi
      ENDDO

   END SUBROUTINE grid_set_rlat

   !---------
   SUBROUTINE grid_free_mem (this)

      IMPLICIT NONE
      TYPE (grid_type) :: this

      IF (allocated (this%lat_s))  deallocate (this%lat_s)
      IF (allocated (this%lat_n))  deallocate (this%lat_n)
      IF (allocated (this%lon_w))  deallocate (this%lon_w)
      IF (allocated (this%lon_e))  deallocate (this%lon_e)
      
      IF (allocated (this%xdsp))   deallocate (this%xdsp) 
      IF (allocated (this%ydsp))   deallocate (this%ydsp) 

      IF (allocated (this%xcnt))   deallocate (this%xcnt) 
      IF (allocated (this%ycnt))   deallocate (this%ycnt) 

      IF (allocated (this%xblk))   deallocate (this%xblk) 
      IF (allocated (this%yblk))   deallocate (this%yblk) 
      
      IF (allocated (this%xloc))   deallocate (this%xloc) 
      IF (allocated (this%yloc))   deallocate (this%yloc) 

      IF (allocated (this%xgrd))   deallocate (this%xgrd) 
      IF (allocated (this%ygrd))   deallocate (this%ygrd) 
      
      IF (allocated (this%rlon))   deallocate (this%rlon) 
      IF (allocated (this%rlat))   deallocate (this%rlat) 

   END SUBROUTINE grid_free_mem

END MODULE mod_grid
