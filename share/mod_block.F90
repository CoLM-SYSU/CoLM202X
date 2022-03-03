#include <define.h>

MODULE mod_block

   USE precision
   IMPLICIT NONE

   ! ---- data types ----
   TYPE :: block_type

      ! Coordinates.
      INTEGER :: nxblk, nyblk
      REAL(r8), allocatable :: lat_s (:)
      REAL(r8), allocatable :: lat_n (:)
      REAL(r8), allocatable :: lon_w (:)
      REAL(r8), allocatable :: lon_e (:)

      ! IO.
      INTEGER, allocatable :: pio(:,:)

   CONTAINS

      procedure, PUBLIC :: set_by_size => block_set_by_size

      procedure, PUBLIC :: save_to_file   => block_save_to_file
      procedure, PUBLIC :: load_from_file => block_load_from_file
      
      final :: block_free_mem

   END TYPE block_type

   ! ---- Instance ----
   TYPE (block_type) :: gblock
   
   ! ---- PUBLIC SUBROUTINE ----
   PUBLIC :: get_filename_block

CONTAINS

   ! --------------------------------
   SUBROUTINE block_set_by_size (this, nxblk_in, nyblk_in)
      
      USE precision
      USE mod_namelist
      USE mod_utils
      USE spmd_task
      IMPLICIT NONE

      class (block_type) :: this
      INTEGER,  intent(in) :: nxblk_in, nyblk_in

      ! Local variables
      INTEGER  :: iblk, jblk, iproc
      INTEGER  :: iblk_south, iblk_north, iblk_west,  iblk_east
      REAL(r8) :: edges, edgen, edgew, edgee
      INTEGER  :: numblocks, numblocks_x, numblocks_y
      INTEGER, parameter :: iset(23) = &
         (/1,  2,  3,  4,  5,  6,  9, 10, 12,  15,  18, &
          20, 24, 30, 36, 40, 45, 60, 72, 90, 120, 180, 360/)
   
      IF ((findloc(iset,nyblk_in,dim=1) <= 0) .or. &
         (findloc(iset,nxblk_in,dim=1) <= 0) ) THEN

         IF (p_is_master) THEN
            write(*,*) 'Number of blocks should be in the set (', iset, ')'
         ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
         CALL mpi_abort   (p_comm_glb, p_err)
#endif

      ENDIF

      this%nxblk = nxblk_in
      this%nyblk = nyblk_in

      allocate (this%lon_w (this%nxblk))
      allocate (this%lon_e (this%nxblk))
      
      DO iblk = 1, this%nxblk
         this%lon_w(iblk) = -180.0 + 360.0/this%nxblk * (iblk-1) 
         this%lon_e(iblk) = -180.0 + 360.0/this%nxblk * iblk

         CALL normalize_longitude (this%lon_w(iblk))
         CALL normalize_longitude (this%lon_e(iblk))
      ENDDO

      allocate (this%lat_s (this%nyblk))
      allocate (this%lat_n (this%nyblk))

      DO jblk = 1, this%nyblk
         this%lat_s(jblk) = -90.0 + 180.0/this%nyblk * (jblk-1) 
         this%lat_n(jblk) = -90.0 + 180.0/this%nyblk * jblk
      ENDDO

      IF (p_is_master) THEN
         write (*,*) 'Block information:'
         write (*,'(I3,A21,I3,A20)') this%nxblk, ' blocks in longitude,', &
            this%nyblk, ' blocks in latitude.'
         write (*,*)
      ENDIF

      allocate (this%pio (this%nxblk,this%nyblk))
      
      IF (p_is_master) THEN
         
         edges = DEF_domain%edges
         edgen = DEF_domain%edgen
         edgew = DEF_domain%edgew
         edgee = DEF_domain%edgee

         iblk_south = find_nearest_south (edges, this%nyblk, this%lat_s)
         iblk_north = find_nearest_north (edgen, this%nyblk, this%lat_n)

         numblocks_y = iblk_north - iblk_south + 1

         CALL normalize_longitude (edgew)
         CALL normalize_longitude (edgee)

         IF (edgew == edgee) THEN
            iblk_west = 1
            iblk_east = this%nxblk
         ELSE
            iblk_west = find_nearest_west (edgew, this%nxblk, this%lon_w)
            iblk_east = find_nearest_east (edgee, this%nxblk, this%lon_e)

            IF (iblk_west == iblk_east) THEN
               IF (lon_between_floor(edgee,this%lon_w(iblk_west),edgew)) THEN
                  iblk_west = 1
                  iblk_east = this%nxblk
               ENDIF
            ENDIF
         ENDIF

         IF (iblk_east >= iblk_west) THEN
            numblocks_x = iblk_east - iblk_west + 1
         ELSE
            numblocks_x = this%nxblk - iblk_west + 1 + iblk_east
         ENDIF

         numblocks = numblocks_x * numblocks_y

      ENDIF
      
#ifdef USEMPI
      CALL mpi_bcast (numblocks, 1, MPI_INTEGER, p_root, p_comm_glb, p_err)
      CALL divide_processes_into_groups (numblocks, DEF_PIO_groupsize)
#endif 
         
      this%pio(:,:) = p_root

#ifdef USEMPI
      IF (p_is_master) THEN

         iproc = -1
         DO jblk = iblk_south, iblk_north

            iblk = iblk_west
            DO while (.true.)
               iproc = mod(iproc+1, p_np_io)
               this%pio(iblk,jblk) = p_address_io(iproc)

               IF (iblk /= iblk_east) THEN
                  iblk = mod(iblk,this%nxblk) + 1
               ELSE
                  exit
               ENDIF
            ENDDO
         ENDDO

      ENDIF

      CALL mpi_bcast (this%pio, this%nxblk * this%nyblk, MPI_INTEGER, &
         p_root, p_comm_glb, p_err)
#endif 

   END SUBROUTINE block_set_by_size

   ! --------------------------------
   SUBROUTINE block_save_to_file (this, dir_landdata)

      USE ncio_serial
      USE spmd_task
      IMPLICIT NONE

      class (block_type) :: this

      CHARACTER(len=*), intent(in) :: dir_landdata

      ! Local variables
      CHARACTER(len=256) :: filename
      
      IF (p_is_master) THEN
         
         filename = trim(dir_landdata) // '/block.nc'

         CALL ncio_create_file (filename)

         CALL ncio_define_dimension (filename, 'longitude', this%nxblk)
         CALL ncio_define_dimension (filename, 'latitude',  this%nyblk)

         CALL ncio_write_serial (filename, 'lat_s', this%lat_s, 'latitude' )
         CALL ncio_write_serial (filename, 'lat_n', this%lat_n, 'latitude' )
         CALL ncio_write_serial (filename, 'lon_w', this%lon_w, 'longitude')
         CALL ncio_write_serial (filename, 'lon_e', this%lon_e, 'longitude')

      ENDIF

   END SUBROUTINE block_save_to_file

   ! --------------------------------
   SUBROUTINE block_load_from_file (this, dir_landdata)

      USE mod_namelist
      USE spmd_task
      USE ncio_serial
      IMPLICIT NONE

      class (block_type) :: this
      CHARACTER(len=*), intent(in) :: dir_landdata

      ! Local variables
      CHARACTER(len=256) :: filename
      INTEGER, allocatable :: nunits_io(:), nunitblk(:,:)
      INTEGER :: numblocks, iblk, jblk, iproc, jproc
         
      filename = trim(dir_landdata) // '/block.nc'
         
      CALL ncio_read_bcast_serial (filename, 'lat_s', this%lat_s)
      CALL ncio_read_bcast_serial (filename, 'lat_n', this%lat_n)
      CALL ncio_read_bcast_serial (filename, 'lon_w', this%lon_w)
      CALL ncio_read_bcast_serial (filename, 'lon_e', this%lon_e)
         
      this%nyblk = size(this%lat_s)
      this%nxblk = size(this%lon_w)

      allocate (this%pio (this%nxblk,this%nyblk))
      this%pio(:,:) = p_root

#ifdef USEMPI
      IF (p_is_master) THEN
         filename = trim(dir_landdata) // '/landunit.nc'
         CALL ncio_read_serial (filename, 'nunits_blk', nunitblk)
         numblocks = count(nunitblk > 0)
      ENDIF 

      CALL mpi_bcast (numblocks, 1, MPI_INTEGER, p_root, p_comm_glb, p_err)
         
      CALL divide_processes_into_groups (numblocks, DEF_PIO_groupsize)

      IF (p_is_master) THEN

         allocate (nunits_io (0:p_np_io-1))
         nunits_io(:) = 0

         jproc = -1
         DO jblk = 1, this%nyblk
            DO iblk = 1, this%nxblk
               IF (nunitblk(iblk,jblk) > 0) THEN
                  iproc = minloc(nunits_io, dim=1) - 1
                  this%pio(iblk,jblk) = p_address_io(iproc)
                  nunits_io(iproc) = nunits_io(iproc) + nunitblk(iblk,jblk)
               ELSEIF (nunitblk(iblk,jblk) == 0) THEN
                  jproc = mod(jproc+1, p_np_io)
                  this%pio(iblk,jblk) = p_address_io(jproc)
               ENDIF
            ENDDO
         ENDDO

         deallocate (nunitblk )
         deallocate (nunits_io)
      ENDIF

      CALL mpi_bcast (this%pio, this%nxblk * this%nyblk, MPI_INTEGER, &
         p_root, p_comm_glb, p_err)
#endif

   END SUBROUTINE block_load_from_file

   ! --------------------------------
   SUBROUTINE block_free_mem (this)

      IMPLICIT NONE
      TYPE (block_type) :: this

      IF (allocated (this%lat_s))  deallocate (this%lat_s)
      IF (allocated (this%lat_n))  deallocate (this%lat_n)
      IF (allocated (this%lon_w))  deallocate (this%lon_w)
      IF (allocated (this%lon_e))  deallocate (this%lon_e)
      
      IF (allocated (this%pio)  )  deallocate (this%pio  )
      
   END SUBROUTINE block_free_mem
   
   ! --------------------------------
   SUBROUTINE get_filename_block (filename, iblk, jblk, fileblock)
      
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      INTEGER, intent(in) :: iblk, jblk

      CHARACTER(len=*), intent(out) :: fileblock

      ! Local variables
      CHARACTER(len=4) :: cx
      CHARACTER(len=3) :: cy
      INTEGER :: n, i

      IF (gblock%lat_s(jblk) < 0) THEN
         write (cy, '(A1,I2.2)') 's', - floor(gblock%lat_s(jblk))
      ELSE
         write (cy, '(A1,I2.2)') 'n',   floor(gblock%lat_s(jblk))
      ENDIF

      IF (gblock%lon_w(iblk) < 0) THEN
         write (cx, '(A1,I3.3)') 'w', - floor(gblock%lon_w(iblk))
      ELSE
         write (cx, '(A1,I3.3)') 'e',   floor(gblock%lon_w(iblk))
      ENDIF

      i = len_trim (filename) 
      DO while (i > 0)
         IF (filename(i:i) == '.') exit
         i = i - 1
      ENDDO

      IF (i > 0) THEN
         fileblock = filename(1:i-1) // '_' // trim(cx) // '_' // trim(cy) // '.nc'
      ELSE
         fileblock = filename // '_' // trim(cx) // '_' // trim(cy) // '.nc'
      ENDIF

   END SUBROUTINE get_filename_block 

END MODULE mod_block
