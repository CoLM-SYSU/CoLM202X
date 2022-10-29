#include <define.h>

MODULE ncio_vector 

   USE mod_data_type
   IMPLICIT NONE

   ! PUBLIC subroutines

   interface ncio_read_vector
      MODULE procedure ncio_read_vector_logical_1d 
      MODULE procedure ncio_read_vector_int32_1d 
      MODULE procedure ncio_read_vector_real8_1d 
      MODULE procedure ncio_read_vector_real8_2d 
      MODULE procedure ncio_read_vector_real8_3d 
      MODULE procedure ncio_read_vector_real8_4d 
   END interface ncio_read_vector

   PUBLIC :: ncio_create_file_vector 
   PUBLIC :: ncio_define_pixelset_dimension 
   PUBLIC :: ncio_define_dimension_vector 

   interface ncio_write_vector
      MODULE procedure ncio_write_vector_logical_1d
      MODULE procedure ncio_write_vector_int32_1d 
      MODULE procedure ncio_write_vector_int32_3d 
      MODULE procedure ncio_write_vector_real8_1d 
      MODULE procedure ncio_write_vector_real8_2d 
      MODULE procedure ncio_write_vector_real8_3d 
      MODULE procedure ncio_write_vector_real8_4d 
   END interface ncio_write_vector

CONTAINS
   
   !---------------------------------------------------------
   SUBROUTINE ncio_read_vector_int32_1d ( &
         filename, dataname, pixelset, rdata, defval)

      USE ncio_serial
      USE spmd_task
      USE mod_block
      USE mod_pixelset
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      TYPE(pixelset_type), intent(in) :: pixelset

      INTEGER, allocatable, intent(inout) :: rdata (:)
      INTEGER, intent(in), optional :: defval

      ! Local variables
      INTEGER :: iblk, jblk, istt, iend
      CHARACTER(len=256) :: fileblock
      INTEGER, allocatable :: sbuff(:), rbuff(:)

      IF (p_is_io) THEN

         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                  IF (pixelset%nonzero(iblk,jblk)) THEN

                     allocate (sbuff (pixelset%vecgs%vlen(iblk,jblk)))
                     CALL get_filename_block (filename, iblk, jblk, fileblock)
                     
                     IF (ncio_var_exist(fileblock,dataname)) THEN 
                        CALL ncio_read_serial (fileblock, dataname, sbuff)
                     ELSEIF (present(defval)) THEN
                        sbuff(:) = defval
                     ENDIF

#ifdef USEMPI
                     CALL mpi_scatterv ( &
                        sbuff, pixelset%vecgs%vcnt(:,iblk,jblk), &
                        pixelset%vecgs%vdsp(:,iblk,jblk), MPI_INTEGER, &
                        MPI_IN_PLACE, 0, MPI_INTEGER, &
                        p_root, p_comm_group, p_err)
#else
                     IF ((pixelset%nset > 0) .and. (.not. allocated(rdata))) THEN
                        allocate (rdata (pixelset%nset))
                     ENDIF

                     istt = pixelset%vecgs%vstt(iblk,jblk)
                     iend = pixelset%vecgs%vend(iblk,jblk)
                     rdata(istt:iend) = sbuff
#endif

                     deallocate (sbuff)

                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN
         
         IF ((pixelset%nset > 0) .and. (.not. allocated(rdata))) THEN
            allocate (rdata (pixelset%nset))
         ENDIF

         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_address_io(p_my_group)) THEN
                  IF (pixelset%nonzero(iblk,jblk)) THEN

                     IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
                        allocate (rbuff (pixelset%vecgs%vlen(iblk,jblk)))
                     ELSE
                        allocate (rbuff(1))
                     ENDIF

                     CALL mpi_scatterv ( &
                        MPI_INULL_P, MPI_INULL_P, MPI_INULL_P, MPI_INTEGER, & ! insignificant on workers
                        rbuff, pixelset%vecgs%vlen(iblk,jblk), MPI_INTEGER, &
                        p_root, p_comm_group, p_err)

                     IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
                        istt = pixelset%vecgs%vstt(iblk,jblk)
                        iend = pixelset%vecgs%vend(iblk,jblk)
                        rdata(istt:iend) = rbuff
                     ENDIF
                        
                     IF (allocated(rbuff)) deallocate (rbuff)

                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF
#endif

   END SUBROUTINE ncio_read_vector_int32_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_vector_logical_1d (filename, dataname, pixelset, rdata, &
         defval)

      USE ncio_serial
      USE spmd_task
      USE mod_block
      USE mod_pixelset
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      TYPE(pixelset_type), intent(in) :: pixelset

      LOGICAL, allocatable, intent(inout) :: rdata (:)
      LOGICAL, intent(in), optional :: defval

      ! Local variables
      INTEGER :: iblk, jblk, istt, iend, i
      CHARACTER(len=256) :: fileblock
      INTEGER(1), allocatable :: sbuff(:), rbuff(:)

      IF (p_is_io) THEN

         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                  IF (pixelset%nonzero(iblk,jblk)) THEN

                     allocate (sbuff (pixelset%vecgs%vlen(iblk,jblk)))
                     CALL get_filename_block (filename, iblk, jblk, fileblock)

                     IF (ncio_var_exist(fileblock,dataname)) THEN 
                        CALL ncio_read_serial (fileblock, dataname, sbuff)
                     ELSEIF (present(defval)) THEN
                        sbuff(:) = defval
                     ENDIF

#ifdef USEMPI
                     CALL mpi_scatterv ( &
                        sbuff, pixelset%vecgs%vcnt(:,iblk,jblk), &
                        pixelset%vecgs%vdsp(:,iblk,jblk), MPI_INTEGER1, &
                        MPI_IN_PLACE, 0, MPI_INTEGER1, &
                        p_root, p_comm_group, p_err)
#else
                     IF ((pixelset%nset > 0) .and. (.not. allocated(rdata))) THEN
                        allocate (rdata (pixelset%nset))
                     ENDIF

                     istt = pixelset%vecgs%vstt(iblk,jblk)
                     iend = pixelset%vecgs%vend(iblk,jblk)
                     rdata(istt:iend) = (sbuff == 1)
#endif

                     deallocate (sbuff)

                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN
         
         IF ((pixelset%nset > 0) .and. (.not. allocated(rdata))) THEN
            allocate (rdata (pixelset%nset))
         ENDIF

         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_address_io(p_my_group)) THEN
                  IF (pixelset%nonzero(iblk,jblk)) THEN

                     IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
                        allocate (rbuff (pixelset%vecgs%vlen(iblk,jblk)))
                     ELSE
                        allocate (rbuff(1))
                     ENDIF

                     CALL mpi_scatterv ( &
                        MPI_INULL_P, MPI_INULL_P, MPI_INULL_P, MPI_INTEGER1, & ! insignificant on workers
                        rbuff, pixelset%vecgs%vlen(iblk,jblk), MPI_INTEGER1, &
                        p_root, p_comm_group, p_err)

                     IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
                        istt = pixelset%vecgs%vstt(iblk,jblk)
                        iend = pixelset%vecgs%vend(iblk,jblk)
                        rdata(istt:iend) = (rbuff == 1)
                     ENDIF
                     
                     IF (allocated(rbuff)) deallocate (rbuff)

                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF
#endif

   END SUBROUTINE ncio_read_vector_logical_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_vector_real8_1d (filename, dataname, pixelset, rdata, &
         defval)

      USE precision
      USE ncio_serial
      USE spmd_task
      USE mod_block
      USE mod_pixelset
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      TYPE(pixelset_type), intent(in) :: pixelset

      REAL(r8), allocatable, intent(inout) :: rdata (:)
      REAL(r8), intent(in), optional :: defval

      ! Local variables
      INTEGER :: iblk, jblk, istt, iend
      CHARACTER(len=256) :: fileblock
      REAL(r8), allocatable :: sbuff(:), rbuff(:)
         
      IF (p_is_worker) THEN
         IF ((pixelset%nset > 0) .and. (.not. allocated(rdata))) THEN
            allocate (rdata (pixelset%nset))
         ENDIF
      ENDIF

      IF (p_is_io) THEN

         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                  IF (pixelset%nonzero(iblk,jblk)) THEN

                     allocate (sbuff (pixelset%vecgs%vlen(iblk,jblk)))
                     CALL get_filename_block (filename, iblk, jblk, fileblock)
                     
                     IF (ncio_var_exist(fileblock,dataname)) THEN 
                        CALL ncio_read_serial (fileblock, dataname, sbuff)
                     ELSEIF (present(defval)) THEN
                        sbuff(:) = defval
                     ENDIF

#ifdef USEMPI
                     CALL mpi_scatterv ( &
                        sbuff, pixelset%vecgs%vcnt(:,iblk,jblk), &
                        pixelset%vecgs%vdsp(:,iblk,jblk), MPI_REAL8, &
                        MPI_IN_PLACE, 0, MPI_REAL8, &
                        p_root, p_comm_group, p_err)
#else
                     istt = pixelset%vecgs%vstt(iblk,jblk)
                     iend = pixelset%vecgs%vend(iblk,jblk)
                     rdata(istt:iend) = sbuff
#endif

                     deallocate (sbuff)

                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN
         
         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_address_io(p_my_group)) THEN
                  IF (pixelset%nonzero(iblk,jblk)) THEN

                     IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
                        allocate (rbuff (pixelset%vecgs%vlen(iblk,jblk)))
                     ELSE
                        allocate (rbuff(1))
                     ENDIF

                     CALL mpi_scatterv ( &
                        MPI_RNULL_P, MPI_INULL_P, MPI_INULL_P, MPI_REAL8, & ! insignificant on workers
                        rbuff, pixelset%vecgs%vlen(iblk,jblk), MPI_REAL8, &
                        p_root, p_comm_group, p_err)

                     IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
                        istt = pixelset%vecgs%vstt(iblk,jblk)
                        iend = pixelset%vecgs%vend(iblk,jblk)
                        rdata(istt:iend) = rbuff
                     ENDIF
                        
                     IF (allocated(rbuff)) deallocate (rbuff)

                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF
#endif

   END SUBROUTINE ncio_read_vector_real8_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_vector_real8_2d ( &
         filename, dataname, ndim1, pixelset, rdata, defval)

      USE precision
      USE ncio_serial
      USE spmd_task
      USE mod_block
      USE mod_pixelset
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      INTEGER, intent(in) :: ndim1
      TYPE(pixelset_type), intent(in) :: pixelset

      REAL(r8), allocatable, intent(inout) :: rdata (:,:)
      REAL(r8), intent(in), optional :: defval

      ! Local variables
      INTEGER :: iblk, jblk, istt, iend
      CHARACTER(len=256) :: fileblock
      REAL(r8), allocatable :: sbuff(:,:), rbuff(:,:)

      IF (p_is_worker) THEN
         IF ((pixelset%nset > 0) .and. (.not. allocated(rdata))) THEN
            allocate (rdata (ndim1, pixelset%nset))
         ENDIF
      ENDIF

      IF (p_is_io) THEN

         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                  IF (pixelset%nonzero(iblk,jblk)) THEN

                     allocate (sbuff (ndim1, pixelset%vecgs%vlen(iblk,jblk)))
                     CALL get_filename_block (filename, iblk, jblk, fileblock)
                     
                     IF (ncio_var_exist(fileblock,dataname)) THEN 
                        CALL ncio_read_serial (fileblock, dataname, sbuff)
                     ELSEIF (present(defval)) THEN
                        sbuff(:,:) = defval
                     ENDIF

#ifdef USEMPI
                     CALL mpi_scatterv ( &
                        sbuff, ndim1 * pixelset%vecgs%vcnt(:,iblk,jblk), &
                        ndim1 * pixelset%vecgs%vdsp(:,iblk,jblk), MPI_REAL8, &
                        MPI_IN_PLACE, 0, MPI_REAL8, &
                        p_root, p_comm_group, p_err)
#else
                     istt = pixelset%vecgs%vstt(iblk,jblk)
                     iend = pixelset%vecgs%vend(iblk,jblk)
                     rdata(:,istt:iend) = sbuff
#endif

                     deallocate (sbuff)

                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN
         
         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_address_io(p_my_group)) THEN
                  IF (pixelset%nonzero(iblk,jblk)) THEN

                     IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
                        allocate (rbuff (ndim1, pixelset%vecgs%vlen(iblk,jblk)))
                     ELSE
                        allocate (rbuff(1,1))
                     ENDIF

                     CALL mpi_scatterv ( &
                        MPI_RNULL_P, MPI_INULL_P, MPI_INULL_P, MPI_REAL8, & ! insignificant on workers
                        rbuff, ndim1 * pixelset%vecgs%vlen(iblk,jblk), MPI_REAL8, &
                        p_root, p_comm_group, p_err)

                     IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
                        istt = pixelset%vecgs%vstt(iblk,jblk)
                        iend = pixelset%vecgs%vend(iblk,jblk)
                        rdata(:,istt:iend) = rbuff
                     ENDIF

                     IF (allocated(rbuff)) deallocate (rbuff)

                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF
#endif

   END SUBROUTINE ncio_read_vector_real8_2d

   !---------------------------------------------------------
   SUBROUTINE ncio_read_vector_real8_3d ( &
         filename, dataname, ndim1, ndim2, pixelset, rdata, defval)

      USE precision
      USE ncio_serial
      USE spmd_task
      USE mod_block
      USE mod_pixelset
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      INTEGER, intent(in) :: ndim1, ndim2
      TYPE(pixelset_type), intent(in) :: pixelset

      REAL(r8), allocatable, intent(inout) :: rdata (:,:,:)
      REAL(r8), intent(in), optional :: defval

      ! Local variables
      INTEGER :: iblk, jblk, istt, iend
      CHARACTER(len=256) :: fileblock
      REAL(r8), allocatable :: sbuff(:,:,:), rbuff(:,:,:)

      IF (p_is_worker) THEN
         IF ((pixelset%nset > 0) .and. (.not. allocated(rdata))) THEN
            allocate (rdata (ndim1,ndim2, pixelset%nset))
         ENDIF
      ENDIF

      IF (p_is_io) THEN

         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                  IF (pixelset%nonzero(iblk,jblk)) THEN

                     allocate (sbuff (ndim1,ndim2, pixelset%vecgs%vlen(iblk,jblk)))
                     CALL get_filename_block (filename, iblk, jblk, fileblock)
                     
                     IF (ncio_var_exist(fileblock,dataname)) THEN 
                        CALL ncio_read_serial (fileblock, dataname, sbuff)
                     ELSEIF (present(defval)) THEN
                        sbuff(:,:,:) = defval
                     ENDIF

#ifdef USEMPI
                     CALL mpi_scatterv ( &
                        sbuff, ndim1 * ndim2 * pixelset%vecgs%vcnt(:,iblk,jblk), &
                        ndim1 * ndim2 * pixelset%vecgs%vdsp(:,iblk,jblk), MPI_REAL8, &
                        MPI_IN_PLACE, 0, MPI_REAL8, &
                        p_root, p_comm_group, p_err)
#else
                     istt = pixelset%vecgs%vstt(iblk,jblk)
                     iend = pixelset%vecgs%vend(iblk,jblk)
                     rdata(:,:,istt:iend) = sbuff
#endif

                     deallocate (sbuff)

                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN

         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_address_io(p_my_group)) THEN
                  IF (pixelset%nonzero(iblk,jblk)) THEN

                     IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
                        allocate (rbuff (ndim1,ndim2, pixelset%vecgs%vlen(iblk,jblk)))
                     ELSE
                        allocate (rbuff(1,1,1))
                     ENDIF

                     CALL mpi_scatterv ( &
                        MPI_RNULL_P, MPI_INULL_P, MPI_INULL_P, MPI_REAL8, & ! insignificant on workers
                        rbuff, ndim1 * ndim2 * pixelset%vecgs%vlen(iblk,jblk), MPI_REAL8, &
                        p_root, p_comm_group, p_err)

                     IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
                        istt = pixelset%vecgs%vstt(iblk,jblk)
                        iend = pixelset%vecgs%vend(iblk,jblk)
                        rdata(:,:,istt:iend) = rbuff
                     ENDIF
                     
                     IF (allocated(rbuff)) deallocate (rbuff)

                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF
#endif

   END SUBROUTINE ncio_read_vector_real8_3d
   
   !---------------------------------------------------------
   SUBROUTINE ncio_read_vector_real8_4d ( &
         filename, dataname, ndim1, ndim2, ndim3, pixelset, rdata, defval)

      USE precision
      USE ncio_serial
      USE spmd_task
      USE mod_block
      USE mod_pixelset
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      INTEGER, intent(in) :: ndim1, ndim2, ndim3
      TYPE(pixelset_type), intent(in) :: pixelset

      REAL(r8), allocatable, intent(inout) :: rdata (:,:,:,:)
      REAL(r8), intent(in), optional :: defval

      ! Local variables
      INTEGER :: iblk, jblk, istt, iend
      CHARACTER(len=256) :: fileblock
      REAL(r8), allocatable :: sbuff(:,:,:,:), rbuff(:,:,:,:)

      IF (p_is_worker) THEN
         IF ((pixelset%nset > 0) .and. (.not. allocated(rdata))) THEN
            allocate (rdata (ndim1,ndim2,ndim3, pixelset%nset))
         ENDIF
      ENDIF

      IF (p_is_io) THEN

         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                  IF (pixelset%nonzero(iblk,jblk)) THEN

                     allocate (sbuff (ndim1,ndim2,ndim3, pixelset%vecgs%vlen(iblk,jblk)))
                     CALL get_filename_block (filename, iblk, jblk, fileblock)
                     
                     IF (ncio_var_exist(fileblock,dataname)) THEN 
                        CALL ncio_read_serial (fileblock, dataname, sbuff)
                     ELSEIF (present(defval)) THEN
                        sbuff(:,:,:,:) = defval
                     ENDIF

#ifdef USEMPI
                     CALL mpi_scatterv ( &
                        sbuff, ndim1 * ndim2 * ndim3 * pixelset%vecgs%vcnt(:,iblk,jblk), &
                        ndim1 * ndim2 * ndim3 * pixelset%vecgs%vdsp(:,iblk,jblk), MPI_REAL8, &
                        MPI_IN_PLACE, 0, MPI_REAL8, &
                        p_root, p_comm_group, p_err)
#else
                     istt = pixelset%vecgs%vstt(iblk,jblk)
                     iend = pixelset%vecgs%vend(iblk,jblk)
                     rdata(:,:,:,istt:iend) = sbuff
#endif

                     deallocate (sbuff)

                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN
         
         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_address_io(p_my_group)) THEN
                  IF (pixelset%nonzero(iblk,jblk)) THEN

                     IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
                        allocate (rbuff (ndim1,ndim2,ndim3, pixelset%vecgs%vlen(iblk,jblk)))
                     ELSE
                        allocate (rbuff(1,1,1,1))
                     ENDIF

                     CALL mpi_scatterv ( &
                        MPI_RNULL_P, MPI_INULL_P, MPI_INULL_P, MPI_REAL8, & ! insignificant on workers
                        rbuff, ndim1 * ndim2 * ndim3 * pixelset%vecgs%vlen(iblk,jblk), MPI_REAL8, &
                        p_root, p_comm_group, p_err)

                     IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
                        istt = pixelset%vecgs%vstt(iblk,jblk)
                        iend = pixelset%vecgs%vend(iblk,jblk)
                        rdata(:,:,:,istt:iend) = rbuff
                     ENDIF
                     
                     IF (allocated(rbuff)) deallocate (rbuff)

                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF
#endif

   END SUBROUTINE ncio_read_vector_real8_4d
   
   !---------------------------------------------------------
   SUBROUTINE ncio_create_file_vector (filename, pixelset)

      USE ncio_serial
      USE spmd_task
      USE mod_block
      USE mod_pixelset
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      TYPE(pixelset_type), intent(in) :: pixelset

      ! Local variables
      INTEGER :: iblk, jblk
      CHARACTER(len=256) :: fileblock

      IF (p_is_io) THEN

         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                  IF (pixelset%nonzero(iblk,jblk)) THEN

                     CALL get_filename_block (filename, iblk, jblk, fileblock)
                     CALL ncio_create_file (fileblock)

                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF
               
   END SUBROUTINE ncio_create_file_vector

   !---------------------------------------------------------
   SUBROUTINE ncio_define_pixelset_dimension (filename, pixelset)

      USE ncio_serial
      USE spmd_task
      USE mod_block
      USE mod_pixelset
      IMPLICIT NONE

      CHARACTER(len=*),    intent(in) :: filename
      TYPE(pixelset_type), intent(in) :: pixelset

      ! Local variables
      INTEGER :: iblk, jblk, dimlen
      CHARACTER(len=256) :: fileblock

      IF (p_is_io) THEN

         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                  IF (pixelset%nonzero(iblk,jblk)) THEN

                     CALL get_filename_block (filename, iblk, jblk, fileblock)
                     dimlen = pixelset%vecgs%vlen(iblk,jblk)
                     CALL ncio_define_dimension (fileblock, 'vector', dimlen)

                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF
               
   END SUBROUTINE ncio_define_pixelset_dimension

   !---------------------------------------------------------
   SUBROUTINE ncio_define_dimension_vector (filename, dimname, dimlen)

      USE ncio_serial
      USE spmd_task
      USE mod_block
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dimname
      INTEGER, intent(in) :: dimlen

      ! Local variables
      INTEGER :: iblk, jblk
      CHARACTER(len=256) :: fileblock
      LOGICAL :: fexists

      IF (p_is_io) THEN
         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN

                  CALL get_filename_block (filename, iblk, jblk, fileblock)
                  inquire (file=trim(fileblock), exist=fexists)
                  IF (fexists) THEN 
                     CALL ncio_define_dimension (fileblock, trim(dimname), dimlen)
                  ENDIF

               ENDIF
            ENDDO
         ENDDO
      ENDIF
               
   END SUBROUTINE ncio_define_dimension_vector

   !---------------------------------------------------------
   SUBROUTINE ncio_write_vector_int32_1d ( &
         filename, dataname, dimname, pixelset, wdata, compress_level)

      USE ncio_serial
      USE spmd_task
      USE mod_block
      USE mod_pixelset
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      CHARACTER(len=*), intent(in) :: dimname
      TYPE(pixelset_type), intent(in) :: pixelset
      INTEGER, intent(in) :: wdata (:)

      INTEGER, intent(in), optional :: compress_level

      ! Local variables
      INTEGER :: iblk, jblk, istt, iend
      CHARACTER(len=256) :: fileblock
      INTEGER, allocatable :: sbuff(:), rbuff(:)

      IF (p_is_io) THEN

         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                  IF (pixelset%nonzero(iblk,jblk)) THEN

                     allocate (rbuff (pixelset%vecgs%vlen(iblk,jblk)))
#ifdef USEMPI
                     CALL mpi_gatherv (MPI_IN_PLACE, 0, MPI_INTEGER, &
                        rbuff, pixelset%vecgs%vcnt(:,iblk,jblk), &
                        pixelset%vecgs%vdsp(:,iblk,jblk), MPI_INTEGER, &
                        p_root, p_comm_group, p_err)
#else
                     istt = pixelset%vecgs%vstt(iblk,jblk)
                     iend = pixelset%vecgs%vend(iblk,jblk)
                     rbuff = wdata(istt:iend)
#endif

                     CALL get_filename_block (filename, iblk, jblk, fileblock)

                     IF (present(compress_level)) THEN
                        CALL ncio_write_serial (fileblock, dataname, rbuff, dimname, &
                           compress = compress_level)
                     ELSE
                        CALL ncio_write_serial (fileblock, dataname, rbuff, dimname)
                     ENDIF

                     deallocate (rbuff)

                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF
               
#ifdef USEMPI
      IF (p_is_worker) THEN

         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_address_io(p_my_group)) THEN
                  IF (pixelset%nonzero(iblk,jblk)) THEN

                     IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
                        allocate (sbuff (pixelset%vecgs%vlen(iblk,jblk)))
                        istt = pixelset%vecgs%vstt(iblk,jblk)
                        iend = pixelset%vecgs%vend(iblk,jblk)
                        sbuff = wdata(istt:iend)
                     ELSE
                        allocate (sbuff (1))
                     ENDIF

                     CALL mpi_gatherv ( &
                        sbuff, pixelset%vecgs%vlen(iblk,jblk), MPI_INTEGER, &
                        MPI_INULL_P, MPI_INULL_P, MPI_INULL_P, MPI_INTEGER, & ! insignificant on workers
                        p_root, p_comm_group, p_err)

                     IF (allocated(sbuff)) deallocate (sbuff)

                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF
#endif

   END SUBROUTINE ncio_write_vector_int32_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_vector_logical_1d ( &
         filename, dataname, dimname, pixelset, wdata, compress_level)

      USE ncio_serial
      USE spmd_task
      USE mod_block
      USE mod_pixelset
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      CHARACTER(len=*), intent(in) :: dimname
      TYPE(pixelset_type), intent(in) :: pixelset
      LOGICAL, intent(in) :: wdata (:)

      INTEGER, intent(in), optional :: compress_level

      ! Local variables
      INTEGER :: iblk, jblk, istt, iend, i
      CHARACTER(len=256) :: fileblock
      INTEGER(1), allocatable :: sbuff(:), rbuff(:)

      IF (p_is_io) THEN

         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                  IF (pixelset%nonzero(iblk,jblk)) THEN

                     allocate (rbuff (pixelset%vecgs%vlen(iblk,jblk)))
#ifdef USEMPI
                     CALL mpi_gatherv (MPI_IN_PLACE, 0, MPI_INTEGER1, &
                        rbuff, pixelset%vecgs%vcnt(:,iblk,jblk), &
                        pixelset%vecgs%vdsp(:,iblk,jblk), MPI_INTEGER1, &
                        p_root, p_comm_group, p_err)
#else
                     istt = pixelset%vecgs%vstt(iblk,jblk)
                     iend = pixelset%vecgs%vend(iblk,jblk)
                     do i = istt, iend
                        if(wdata(i))then
                           rbuff(i-istt+1) = 1
                        else
                           rbuff(i-istt+1) = 0
                        end if
                     end do
#endif

                     CALL get_filename_block (filename, iblk, jblk, fileblock)

                     IF (present(compress_level)) THEN
                        CALL ncio_write_serial (fileblock, dataname, rbuff, dimname, &
                           compress = compress_level)
                     ELSE
                        CALL ncio_write_serial (fileblock, dataname, rbuff, dimname)
                     ENDIF

                     deallocate (rbuff)

                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF
               
#ifdef USEMPI
      IF (p_is_worker) THEN

         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_address_io(p_my_group)) THEN
                  IF (pixelset%nonzero(iblk,jblk)) THEN

                     IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
                        allocate (sbuff (pixelset%vecgs%vlen(iblk,jblk)))
                        istt = pixelset%vecgs%vstt(iblk,jblk)
                        iend = pixelset%vecgs%vend(iblk,jblk)
                        do i = istt, iend
                           if(wdata(i))then
                             sbuff(i-istt+1) = 1
                           else
                             sbuff(i-istt+1) = 0
                           end if
                        end do
                     ELSE
                        allocate (sbuff (1))
                     ENDIF

                     CALL mpi_gatherv ( &
                        sbuff, pixelset%vecgs%vlen(iblk,jblk), MPI_INTEGER1, &
                        MPI_INULL_P, MPI_INULL_P, MPI_INULL_P, MPI_INTEGER1, & ! insignificant on workers
                        p_root, p_comm_group, p_err)

                     IF (allocated(sbuff)) deallocate (sbuff)

                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF
#endif

   END SUBROUTINE ncio_write_vector_logical_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_vector_int32_3d ( &
         filename, dataname, dim1name, ndim1, dim2name, ndim2, &
         dim3name, pixelset, wdata, compress_level)

      USE ncio_serial
      USE spmd_task
      USE mod_block
      USE mod_pixelset
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      CHARACTER(len=*), intent(in) :: dim1name, dim2name, dim3name
      TYPE(pixelset_type), intent(in) :: pixelset
      INTEGER, intent(in) :: ndim1, ndim2
      INTEGER, intent(in) :: wdata (:,:,:)

      INTEGER, intent(in), optional :: compress_level

      ! Local variables
      INTEGER :: iblk, jblk, iproc, istt, iend
      CHARACTER(len=256) :: fileblock
      INTEGER, allocatable :: sbuff(:,:,:), rbuff(:,:,:)

      IF (p_is_io) THEN

         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                  IF (pixelset%nonzero(iblk,jblk)) THEN
                     
                     allocate (rbuff (ndim1,ndim2,pixelset%vecgs%vlen(iblk,jblk)))
#ifdef USEMPI
                     CALL mpi_gatherv (MPI_IN_PLACE, 0, MPI_INTEGER, &
                        rbuff, ndim1*ndim2*pixelset%vecgs%vcnt(:,iblk,jblk), &
                        ndim1*ndim2*pixelset%vecgs%vdsp(:,iblk,jblk), MPI_INTEGER, &
                        p_root, p_comm_group, p_err)
#else
                     istt = pixelset%vecgs%vstt(iblk,jblk)
                     iend = pixelset%vecgs%vend(iblk,jblk)
                     rbuff = wdata(:,:,istt:iend)
#endif

                     CALL get_filename_block (filename, iblk, jblk, fileblock)

                     IF (present(compress_level)) THEN
                        CALL ncio_write_serial (fileblock, dataname, rbuff, &
                           dim1name, dim2name, dim3name, compress = compress_level)
                     ELSE
                        CALL ncio_write_serial (fileblock, dataname, rbuff, &
                           dim1name, dim2name, dim3name)
                     ENDIF

                     deallocate (rbuff)

                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF
               
#ifdef USEMPI
      IF (p_is_worker) THEN

         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_address_io(p_my_group)) THEN
                  IF (pixelset%nonzero(iblk,jblk)) THEN

                     IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
                        allocate (sbuff (ndim1,ndim2, pixelset%vecgs%vlen(iblk,jblk)))
                        istt = pixelset%vecgs%vstt(iblk,jblk)
                        iend = pixelset%vecgs%vend(iblk,jblk)
                        sbuff = wdata(:,:,istt:iend)
                     ELSE
                        allocate (sbuff (1,1,1))
                     ENDIF

                     CALL mpi_gatherv ( &
                        sbuff, ndim1*ndim2*pixelset%vecgs%vlen(iblk,jblk), MPI_INTEGER, &
                        MPI_INULL_P, MPI_INULL_P, MPI_INULL_P, MPI_INTEGER, & ! insignificant on workers
                        p_root, p_comm_group, p_err)

                     IF (allocated(sbuff)) deallocate (sbuff)

                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF
#endif

   END SUBROUTINE ncio_write_vector_int32_3d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_vector_real8_1d ( &
         filename, dataname, dimname, pixelset, wdata, compress_level)

      USE precision
      USE ncio_serial
      USE spmd_task
      USE mod_block
      USE mod_pixelset
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      CHARACTER(len=*), intent(in) :: dimname
      TYPE(pixelset_type), intent(in) :: pixelset
      REAL(r8), intent(in) :: wdata (:)
      
      INTEGER, intent(in), optional :: compress_level

      ! Local variables
      INTEGER :: iblk, jblk, istt, iend
      CHARACTER(len=256) :: fileblock
      REAL(r8), allocatable :: sbuff(:), rbuff(:)

      IF (p_is_io) THEN

         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                  IF (pixelset%nonzero(iblk,jblk)) THEN

                     allocate (rbuff (pixelset%vecgs%vlen(iblk,jblk)))
#ifdef USEMPI
                     CALL mpi_gatherv ( MPI_IN_PLACE, 0, MPI_REAL8, &
                        rbuff, pixelset%vecgs%vcnt(:,iblk,jblk), &
                        pixelset%vecgs%vdsp(:,iblk,jblk), MPI_REAL8, &
                        p_root, p_comm_group, p_err)
#else
                     istt = pixelset%vecgs%vstt(iblk,jblk)
                     iend = pixelset%vecgs%vend(iblk,jblk)
                     rbuff = wdata(istt:iend)
#endif

                     CALL get_filename_block (filename, iblk, jblk, fileblock)
                     IF (present(compress_level)) THEN
                        CALL ncio_write_serial (fileblock, dataname, rbuff, &
                           dimname, compress = compress_level)
                     ELSE
                        CALL ncio_write_serial (fileblock, dataname, rbuff, dimname)
                     ENDIF

                     deallocate (rbuff)

                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF
               
#ifdef USEMPI
      IF (p_is_worker) THEN

         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_address_io(p_my_group)) THEN
                  IF (pixelset%nonzero(iblk,jblk)) THEN

                     IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
                        allocate (sbuff (pixelset%vecgs%vlen(iblk,jblk)))
                        istt = pixelset%vecgs%vstt(iblk,jblk)
                        iend = pixelset%vecgs%vend(iblk,jblk)
                        sbuff = wdata(istt:iend)
                     ELSE
                        allocate (sbuff (1))
                     ENDIF

                     CALL mpi_gatherv ( &
                        sbuff, pixelset%vecgs%vlen(iblk,jblk), MPI_REAL8, &
                        MPI_RNULL_P, MPI_INULL_P, MPI_INULL_P, MPI_REAL8, & ! insignificant on workers
                        p_root, p_comm_group, p_err)

                     IF (allocated(sbuff)) deallocate (sbuff)

                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF
#endif

   END SUBROUTINE ncio_write_vector_real8_1d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_vector_real8_2d ( &
         filename, dataname, dim1name, ndim1, &
         dim2name, pixelset, wdata, compress_level)

      USE precision
      USE ncio_serial
      USE spmd_task
      USE mod_block
      USE mod_pixelset
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      CHARACTER(len=*), intent(in) :: dim1name, dim2name
      INTEGER,  intent(in) :: ndim1
      TYPE(pixelset_type), intent(in) :: pixelset
      REAL(r8), intent(in) :: wdata (:,:)

      INTEGER,  intent(in), optional :: compress_level

      ! Local variables
      INTEGER :: iblk, jblk, istt, iend
      CHARACTER(len=256) :: fileblock
      REAL(r8), allocatable :: sbuff(:,:), rbuff(:,:)

      IF (p_is_io) THEN

         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                  IF (pixelset%nonzero(iblk,jblk)) THEN

                     allocate (rbuff (ndim1, pixelset%vecgs%vlen(iblk,jblk)))
#ifdef USEMPI
                     CALL mpi_gatherv (MPI_IN_PLACE, 0, MPI_REAL8, &
                        rbuff, ndim1 * pixelset%vecgs%vcnt(:,iblk,jblk), &
                        ndim1 * pixelset%vecgs%vdsp(:,iblk,jblk), MPI_REAL8, &
                        p_root, p_comm_group, p_err)
#else
                     istt = pixelset%vecgs%vstt(iblk,jblk)
                     iend = pixelset%vecgs%vend(iblk,jblk)
                     rbuff = wdata(:,istt:iend)
#endif

                     CALL get_filename_block (filename, iblk, jblk, fileblock)

                     IF (present(compress_level)) THEN
                        CALL ncio_write_serial (fileblock, dataname, rbuff, &
                           dim1name, dim2name, compress = compress_level)
                     ELSE
                        CALL ncio_write_serial (fileblock, dataname, rbuff, &
                           dim1name, dim2name)
                     ENDIF

                     deallocate (rbuff)

                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF
               
#ifdef USEMPI
      IF (p_is_worker) THEN

         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_address_io(p_my_group)) THEN
                  IF (pixelset%nonzero(iblk,jblk)) THEN

                     IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
                        allocate (sbuff (ndim1,pixelset%vecgs%vlen(iblk,jblk)))
                        istt = pixelset%vecgs%vstt(iblk,jblk)
                        iend = pixelset%vecgs%vend(iblk,jblk)
                        sbuff = wdata(:,istt:iend)
                     ELSE
                        allocate (sbuff (1,1))
                     ENDIF

                     CALL mpi_gatherv ( &
                        sbuff, ndim1 * pixelset%vecgs%vlen(iblk,jblk), MPI_REAL8, &
                        MPI_RNULL_P, MPI_INULL_P, MPI_INULL_P, MPI_REAL8, & ! insignificant on workers
                        p_root, p_comm_group, p_err)

                     IF (allocated(sbuff)) deallocate (sbuff)

                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF
#endif

   END SUBROUTINE ncio_write_vector_real8_2d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_vector_real8_3d ( &
         filename, dataname, dim1name, ndim1, dim2name, ndim2, &
         dim3name, pixelset, wdata, compress_level)

      USE precision
      USE ncio_serial
      USE spmd_task
      USE mod_block
      USE mod_pixelset
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      CHARACTER(len=*), intent(in) :: dim1name, dim2name, dim3name
      TYPE(pixelset_type), intent(in) :: pixelset
      INTEGER,  intent(in) :: ndim1, ndim2
      REAL(r8), intent(in) :: wdata (:,:,:)
      
      INTEGER,  intent(in), optional :: compress_level

      ! Local variables
      INTEGER :: iblk, jblk, istt, iend
      CHARACTER(len=256) :: fileblock
      REAL(r8), allocatable :: sbuff(:,:,:), rbuff(:,:,:)

      IF (p_is_io) THEN

         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                  IF (pixelset%nonzero(iblk,jblk)) THEN

                     allocate (rbuff (ndim1, ndim2, pixelset%vecgs%vlen(iblk,jblk)))
#ifdef USEMPI
                     CALL mpi_gatherv (MPI_IN_PLACE, 0, MPI_REAL8, &
                        rbuff, ndim1 * ndim2 * pixelset%vecgs%vcnt(:,iblk,jblk), &
                        ndim1 * ndim2 * pixelset%vecgs%vdsp(:,iblk,jblk), MPI_REAL8, &
                        p_root, p_comm_group, p_err)
#else
                     istt = pixelset%vecgs%vstt(iblk,jblk)
                     iend = pixelset%vecgs%vend(iblk,jblk)
                     rbuff = wdata(:,:,istt:iend)
#endif

                     CALL get_filename_block (filename, iblk, jblk, fileblock)
                     IF (present(compress_level)) THEN
                        CALL ncio_write_serial (fileblock, dataname, rbuff, &
                           dim1name, dim2name, dim3name, compress = compress_level)
                     ELSE
                        CALL ncio_write_serial (fileblock, dataname, rbuff, &
                           dim1name, dim2name, dim3name)
                     ENDIF

                     deallocate (rbuff)

                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF
               
#ifdef USEMPI
      IF (p_is_worker) THEN

         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_address_io(p_my_group)) THEN
                  IF (pixelset%nonzero(iblk,jblk)) THEN

                     IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
                        allocate (sbuff (ndim1,ndim2,pixelset%vecgs%vlen(iblk,jblk)))
                        istt = pixelset%vecgs%vstt(iblk,jblk)
                        iend = pixelset%vecgs%vend(iblk,jblk)
                        sbuff = wdata(:,:,istt:iend)
                     ELSE
                        allocate (sbuff (1,1,1))
                     ENDIF

                     CALL mpi_gatherv ( sbuff, &
                        ndim1 * ndim2 * pixelset%vecgs%vlen(iblk,jblk), MPI_REAL8, &
                        MPI_RNULL_P, MPI_INULL_P, MPI_INULL_P, MPI_REAL8, & ! insignificant on workers
                        p_root, p_comm_group, p_err)

                     IF (allocated(sbuff)) deallocate (sbuff)

                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF
#endif

   END SUBROUTINE ncio_write_vector_real8_3d

   !---------------------------------------------------------
   SUBROUTINE ncio_write_vector_real8_4d ( &
         filename, dataname, dim1name, ndim1, dim2name, ndim2, dim3name, ndim3, &
         dim4name, pixelset, wdata, compress_level)

      USE precision
      USE ncio_serial
      USE spmd_task
      USE mod_block
      USE mod_pixelset
      IMPLICIT NONE

      CHARACTER(len=*), intent(in) :: filename
      CHARACTER(len=*), intent(in) :: dataname
      CHARACTER(len=*), intent(in) :: dim1name, dim2name, dim3name, dim4name
      INTEGER,  intent(in) :: ndim1, ndim2, ndim3
      TYPE(pixelset_type), intent(in) :: pixelset
      REAL(r8), intent(in) :: wdata (:,:,:,:)
      
      INTEGER,  intent(in), optional :: compress_level

      ! Local variables
      INTEGER :: iblk, jblk, istt, iend
      CHARACTER(len=256) :: fileblock
      REAL(r8), allocatable :: sbuff(:,:,:,:), rbuff(:,:,:,:)

      IF (p_is_io) THEN

         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_iam_glb) THEN
                  IF (pixelset%nonzero(iblk,jblk)) THEN

                     allocate (rbuff (ndim1, ndim2, ndim3, pixelset%vecgs%vlen(iblk,jblk)))
#ifdef USEMPI
                     CALL mpi_gatherv (MPI_IN_PLACE, 0, MPI_REAL8, &
                        rbuff, ndim1 * ndim2 * ndim3 * pixelset%vecgs%vcnt(:,iblk,jblk), &
                        ndim1 * ndim2 * ndim3 * pixelset%vecgs%vdsp(:,iblk,jblk), MPI_REAL8, &
                        p_root, p_comm_group, p_err)
#else
                     istt = pixelset%vecgs%vstt(iblk,jblk)
                     iend = pixelset%vecgs%vend(iblk,jblk)
                     rbuff = wdata(:,:,:,istt:iend)
#endif

                     CALL get_filename_block (filename, iblk, jblk, fileblock)
                     IF (present(compress_level)) THEN
                        CALL ncio_write_serial (fileblock, dataname, rbuff, &
                           dim1name, dim2name, dim3name, dim4name, compress = compress_level)
                     ELSE
                        CALL ncio_write_serial (fileblock, dataname, rbuff, &
                           dim1name, dim2name, dim3name, dim4name)
                     ENDIF

                     deallocate (rbuff)

                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF
               
#ifdef USEMPI
      IF (p_is_worker) THEN

         DO jblk = 1, gblock%nyblk
            DO iblk = 1, gblock%nxblk
               IF (gblock%pio(iblk,jblk) == p_address_io(p_my_group)) THEN
                  IF (pixelset%nonzero(iblk,jblk)) THEN

                     IF (pixelset%vecgs%vlen(iblk,jblk) > 0) THEN
                        allocate (sbuff (ndim1,ndim2,ndim3,pixelset%vecgs%vlen(iblk,jblk)))
                        istt = pixelset%vecgs%vstt(iblk,jblk)
                        iend = pixelset%vecgs%vend(iblk,jblk)
                        sbuff = wdata(:,:,:,istt:iend)
                     ELSE
                        allocate (sbuff (1,1,1,1))
                     ENDIF

                     CALL mpi_gatherv ( sbuff, &
                        ndim1 * ndim2 * ndim3 * pixelset%vecgs%vlen(iblk,jblk), MPI_REAL8, &
                        MPI_RNULL_P, MPI_INULL_P, MPI_INULL_P, MPI_REAL8, & ! insignificant on workers
                        p_root, p_comm_group, p_err)

                     IF (allocated(sbuff)) deallocate (sbuff)

                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF
#endif

   END SUBROUTINE ncio_write_vector_real8_4d

END MODULE ncio_vector
