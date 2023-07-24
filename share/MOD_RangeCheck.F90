#include <define.h>

MODULE MOD_RangeCheck

   !-----------------------------------------------------------------------
   ! DESCRIPTION:
   !
   !    Subroutines show the range of values in block data or vector data.
   !
   !    Notice that:
   !    1. "check_block_data"  can only be called by IO     processes.
   !    2. "check_vector_data" can only be called by worker processes.
   !
   ! Created by Shupeng Zhang, May 2023
   !-----------------------------------------------------------------------

#ifdef RangeCheck
   IMPLICIT NONE

   interface check_block_data
      MODULE procedure check_block_data_real8_2d
   END interface check_block_data

   interface check_vector_data
      MODULE procedure check_vector_data_real8_1d
      MODULE procedure check_vector_data_real8_2d
      MODULE procedure check_vector_data_real8_3d
      MODULE procedure check_vector_data_real8_4d
      MODULE procedure check_vector_data_int32_1d
   END interface check_vector_data

CONTAINS

   ! ----------
   SUBROUTINE check_block_data_real8_2d (varname, gdata, spv_in, largevalue)

      USE MOD_Precision
      USE MOD_SPMD_Task
      USE MOD_Block
      USE MOD_DataType
      IMPLICIT NONE

      CHARACTER(len=*), intent(in)   :: varname
      TYPE(block_data_real8_2d), intent(in) :: gdata
      REAL(r8), intent(in), optional :: spv_in
      REAL(r8), intent(in), optional :: largevalue

      ! Local variables
      REAL(r8) :: gmin, gmax, spv
      REAL(r8), allocatable :: gmin_all(:), gmax_all(:)
      LOGICAL,  allocatable :: msk2(:,:)
      INTEGER :: iblkme, ib, jb, ix, iy
      LOGICAL :: has_nan
      character(len=256) :: wfmt, ss, info

      IF (p_is_io) THEN

         IF (present(spv_in)) THEN
            spv = spv_in
         ELSE
            spv = -1.0e36_r8
         ENDIF

         gmin = spv
         gmax = spv

         has_nan = .false.
         DO iblkme = 1, gblock%nblkme
            ib = gblock%xblkme(iblkme)
            jb = gblock%yblkme(iblkme)

            IF (.not. allocated(gdata%blk(ib,jb)%val)) cycle

            allocate(msk2 (size(gdata%blk(ib,jb)%val,1), size(gdata%blk(ib,jb)%val,2)))
            msk2 = gdata%blk(ib,jb)%val /= spv

            IF (any(msk2)) THEN
               IF (gmin == spv) THEN
                  gmin = minval(gdata%blk(ib,jb)%val, mask = msk2)
               ELSE
                  gmin = min(gmin, minval(gdata%blk(ib,jb)%val, mask = msk2))
               ENDIF

               IF (gmax == spv) THEN
                  gmax = maxval(gdata%blk(ib,jb)%val, mask = msk2)
               ELSE
                  gmax = max(gmax, maxval(gdata%blk(ib,jb)%val, mask = msk2))
               ENDIF
            ENDIF

            DO iy = 1, size(gdata%blk(ib,jb)%val,2)
               DO ix = 1, size(gdata%blk(ib,jb)%val,1)
                  has_nan = has_nan .or. isnan(gdata%blk(ib,jb)%val(ix,iy))
               ENDDO
            ENDDO

            deallocate(msk2)

         ENDDO

#ifdef USEMPI
         IF (p_iam_io == p_root) THEN
            allocate (gmin_all (0:p_np_io-1))
            allocate (gmax_all (0:p_np_io-1))
            CALL mpi_gather (gmin, 1, MPI_REAL8, gmin_all, 1, MPI_REAL8, p_root, p_comm_io, p_err)
            CALL mpi_gather (gmax, 1, MPI_REAL8, gmax_all, 1, MPI_REAL8, p_root, p_comm_io, p_err)
         ELSE
            CALL mpi_gather (gmin, 1, MPI_REAL8, MPI_RNULL_P, 1, MPI_REAL8, p_root, p_comm_io, p_err)
            CALL mpi_gather (gmax, 1, MPI_REAL8, MPI_RNULL_P, 1, MPI_REAL8, p_root, p_comm_io, p_err)
         ENDIF

         CALL mpi_allreduce (MPI_IN_PLACE, has_nan, 1, MPI_LOGICAL, MPI_LOR, p_comm_io, p_err)

         IF (p_iam_io == p_root) THEN
            IF (any(gmin_all /= spv)) THEN
               gmin = minval(gmin_all, mask = (gmin_all /= spv))
            ELSE
               gmin = spv
            ENDIF

            IF (any(gmax_all /= spv)) THEN
               gmax = maxval(gmax_all, mask = (gmax_all /= spv))
            ELSE
               gmax = spv
            ENDIF

            deallocate (gmin_all)
            deallocate (gmax_all)
         ENDIF
#endif
         IF (p_iam_io == p_root) THEN

            info = ''

            IF (has_nan) THEN
               info = trim(info) // ' with NAN'
            ENDIF

            IF (present(largevalue)) THEN
               IF (max(abs(gmin),abs(gmax)) > largevalue) THEN
                  write(ss,'(e12.2)') largevalue
                  info = trim(info) // ' with value > ' // trim(ss)
               ENDIF
            ENDIF

            wfmt = "('Check block  data:', A25, ' is in (', e20.10, ',', e20.10, ')', A)"
            write(*,wfmt) varname, gmin, gmax, trim(info)

#if(defined CoLMDEBUG)
            IF (len_trim(info) > 0) THEN
#ifdef USEMPI
               CALL mpi_abort (p_comm_glb, p_err)
#else
               STOP
#endif
            ENDIF
#endif

         ENDIF

      ENDIF

   END SUBROUTINE check_block_data_real8_2d


   ! ----------
   SUBROUTINE check_vector_data_real8_1d (varname, vdata, spv_in, largevalue)

      USE MOD_Precision
      USE MOD_SPMD_Task
      IMPLICIT NONE

      CHARACTER(len=*), intent(in)   :: varname
      REAL(r8), intent(in)           :: vdata(:)
      REAL(r8), intent(in), optional :: spv_in
      REAL(r8), intent(in), optional :: largevalue

      ! Local variables
      REAL(r8) :: vmin, vmax, spv
      REAL(r8), allocatable :: vmin_all(:), vmax_all(:)
      INTEGER  :: i
      LOGICAL  :: has_nan
      character(len=256) :: wfmt, ss, info

      IF (p_is_worker) THEN

         IF (present(spv_in)) THEN
            spv = spv_in
         ELSE
            spv = -1.0e36_r8
         ENDIF

         IF (any(vdata /= spv)) THEN
            vmin = minval(vdata, mask = vdata /= spv)
            vmax = maxval(vdata, mask = vdata /= spv)
         ELSE
            vmin = spv
            vmax = spv
         ENDIF

         has_nan = .false.
         DO i = 1, size(vdata)
            has_nan = has_nan .or. isnan(vdata(i))
         ENDDO

#ifdef USEMPI
         IF (p_iam_worker == p_root) THEN
            allocate (vmin_all (0:p_np_worker-1))
            allocate (vmax_all (0:p_np_worker-1))
            CALL mpi_gather (vmin, 1, MPI_REAL8, vmin_all, 1, MPI_REAL8, p_root, p_comm_worker, p_err)
            CALL mpi_gather (vmax, 1, MPI_REAL8, vmax_all, 1, MPI_REAL8, p_root, p_comm_worker, p_err)
         ELSE
            CALL mpi_gather (vmin, 1, MPI_REAL8, MPI_RNULL_P, 1, MPI_REAL8, p_root, p_comm_worker, p_err)
            CALL mpi_gather (vmax, 1, MPI_REAL8, MPI_RNULL_P, 1, MPI_REAL8, p_root, p_comm_worker, p_err)
         ENDIF

         CALL mpi_allreduce (MPI_IN_PLACE, has_nan, 1, MPI_LOGICAL, MPI_LOR, p_comm_worker, p_err)

         IF (p_iam_worker == p_root) THEN
            IF (any(vmin_all /= spv)) THEN
               vmin = minval(vmin_all, mask = (vmin_all /= spv))
            ELSE
               vmin = spv
            ENDIF

            IF (any(vmax_all /= spv)) THEN
               vmax = maxval(vmax_all, mask = (vmax_all /= spv))
            ELSE
               vmax = spv
            ENDIF

            deallocate (vmin_all)
            deallocate (vmax_all)
         ENDIF
#endif

         IF (p_iam_worker == p_root) THEN

            info = ''

            IF (has_nan) THEN
               info = trim(info) // ' with NAN'
            ENDIF

            IF (present(largevalue)) THEN
               IF (max(abs(vmin),abs(vmax)) > largevalue) THEN
                  write(ss,'(e12.2)') largevalue
                  info = trim(info) // ' with value > ' // trim(ss)
               ENDIF
            ENDIF

            wfmt = "('Check vector data:', A25, ' is in (', e20.10, ',', e20.10, ')', A)"
            write(*,wfmt) varname, vmin, vmax, trim(info)

#if(defined CoLMDEBUG)
            IF (len_trim(info) > 0) THEN
#ifdef USEMPI
               CALL mpi_abort (p_comm_glb, p_err)
#else
               STOP
#endif
            ENDIF
#endif

         ENDIF

      ENDIF

   END SUBROUTINE check_vector_data_real8_1d

   ! ----------
   SUBROUTINE check_vector_data_real8_2d (varname, vdata, spv_in, largevalue)

      USE MOD_Precision
      USE MOD_SPMD_Task
      IMPLICIT NONE

      CHARACTER(len=*), intent(in)   :: varname
      REAL(r8), intent(in)           :: vdata(:,:)
      REAL(r8), intent(in), optional :: spv_in
      REAL(r8), intent(in), optional :: largevalue

      ! Local variables
      REAL(r8) :: vmin, vmax, spv
      REAL(r8), allocatable :: vmin_all(:), vmax_all(:)
      INTEGER  :: i, j
      LOGICAL  :: has_nan
      character(len=256) :: wfmt, ss, info

      IF (p_is_worker) THEN

         IF (present(spv_in)) THEN
            spv = spv_in
         ELSE
            spv = -1.0e36_r8
         ENDIF

         IF (any(vdata /= spv)) THEN
            vmin = minval(vdata, mask = vdata /= spv)
            vmax = maxval(vdata, mask = vdata /= spv)
         ELSE
            vmin = spv
            vmax = spv
         ENDIF

         has_nan = .false.
         DO j = 1, size(vdata,2)
            DO i = 1, size(vdata,1)
               has_nan = has_nan .or. isnan(vdata(i,j))
            ENDDO
         ENDDO

#ifdef USEMPI
         IF (p_iam_worker == p_root) THEN
            allocate (vmin_all (0:p_np_worker-1))
            allocate (vmax_all (0:p_np_worker-1))
            CALL mpi_gather (vmin, 1, MPI_REAL8, vmin_all, 1, MPI_REAL8, p_root, p_comm_worker, p_err)
            CALL mpi_gather (vmax, 1, MPI_REAL8, vmax_all, 1, MPI_REAL8, p_root, p_comm_worker, p_err)
         ELSE
            CALL mpi_gather (vmin, 1, MPI_REAL8, MPI_RNULL_P, 1, MPI_REAL8, p_root, p_comm_worker, p_err)
            CALL mpi_gather (vmax, 1, MPI_REAL8, MPI_RNULL_P, 1, MPI_REAL8, p_root, p_comm_worker, p_err)
         ENDIF

         CALL mpi_allreduce (MPI_IN_PLACE, has_nan, 1, MPI_LOGICAL, MPI_LOR, p_comm_worker, p_err)

         IF (p_iam_worker == p_root) THEN
            IF (any(vmin_all /= spv)) THEN
               vmin = minval(vmin_all, mask = (vmin_all /= spv))
            ELSE
               vmin = spv
            ENDIF

            IF (any(vmax_all /= spv)) THEN
               vmax = maxval(vmax_all, mask = (vmax_all /= spv))
            ELSE
               vmax = spv
            ENDIF

            deallocate (vmin_all)
            deallocate (vmax_all)
         ENDIF
#endif

         IF (p_iam_worker == p_root) THEN

            info = ''

            IF (has_nan) THEN
               info = trim(info) // ' with NAN'
            ENDIF

            IF (present(largevalue)) THEN
               IF (max(abs(vmin),abs(vmax)) > largevalue) THEN
                  write(ss,'(e12.2)') largevalue
                  info = trim(info) // ' with value > ' // trim(ss)
               ENDIF
            ENDIF

            wfmt = "('Check vector data:', A25, ' is in (', e20.10, ',', e20.10, ')', A)"
            write(*,wfmt) varname, vmin, vmax, trim(info)

#if(defined CoLMDEBUG)
            IF (len_trim(info) > 0) THEN
#ifdef USEMPI
               CALL mpi_abort (p_comm_glb, p_err)
#else
               STOP
#endif
            ENDIF
#endif

         ENDIF

      ENDIF

   END SUBROUTINE check_vector_data_real8_2d

   ! ----------
   SUBROUTINE check_vector_data_real8_3d (varname, vdata, spv_in, largevalue)

      USE MOD_Precision
      USE MOD_SPMD_Task
      IMPLICIT NONE

      CHARACTER(len=*), intent(in)   :: varname
      REAL(r8), intent(in)           :: vdata(:,:,:)
      REAL(r8), intent(in), optional :: spv_in
      REAL(r8), intent(in), optional :: largevalue

      ! Local variables
      REAL(r8) :: vmin, vmax, spv
      REAL(r8), allocatable :: vmin_all(:), vmax_all(:)
      INTEGER  :: i, j, k
      LOGICAL  :: has_nan
      character(len=256) :: wfmt, ss, info

      IF (p_is_worker) THEN

         IF (present(spv_in)) THEN
            spv = spv_in
         ELSE
            spv = -1.0e36_r8
         ENDIF

         IF (any(vdata /= spv)) THEN
            vmin = minval(vdata, mask = vdata /= spv)
            vmax = maxval(vdata, mask = vdata /= spv)
         ELSE
            vmin = spv
            vmax = spv
         ENDIF

         has_nan = .false.
         DO k = 1, size(vdata,3)
            DO j = 1, size(vdata,2)
               DO i = 1, size(vdata,1)
                  has_nan = has_nan .or. isnan(vdata(i,j,k))
               ENDDO
            ENDDO
         ENDDO

#ifdef USEMPI
         IF (p_iam_worker == p_root) THEN
            allocate (vmin_all (0:p_np_worker-1))
            allocate (vmax_all (0:p_np_worker-1))
            CALL mpi_gather (vmin, 1, MPI_REAL8, vmin_all, 1, MPI_REAL8, p_root, p_comm_worker, p_err)
            CALL mpi_gather (vmax, 1, MPI_REAL8, vmax_all, 1, MPI_REAL8, p_root, p_comm_worker, p_err)
         ELSE
            CALL mpi_gather (vmin, 1, MPI_REAL8, MPI_RNULL_P, 1, MPI_REAL8, p_root, p_comm_worker, p_err)
            CALL mpi_gather (vmax, 1, MPI_REAL8, MPI_RNULL_P, 1, MPI_REAL8, p_root, p_comm_worker, p_err)
         ENDIF


         CALL mpi_allreduce (MPI_IN_PLACE, has_nan, 1, MPI_LOGICAL, MPI_LOR, p_comm_worker, p_err)

         IF (p_iam_worker == p_root) THEN
            IF (any(vmin_all /= spv)) THEN
               vmin = minval(vmin_all, mask = (vmin_all /= spv))
            ELSE
               vmin = spv
            ENDIF

            IF (any(vmax_all /= spv)) THEN
               vmax = maxval(vmax_all, mask = (vmax_all /= spv))
            ELSE
               vmax = spv
            ENDIF

            deallocate (vmin_all)
            deallocate (vmax_all)
         ENDIF
#endif

         IF (p_iam_worker == p_root) THEN

            info = ''

            IF (has_nan) THEN
               info = trim(info) // ' with NAN'
            ENDIF

            IF (present(largevalue)) THEN
               IF (max(abs(vmin),abs(vmax)) > largevalue) THEN
                  write(ss,'(e12.2)') largevalue
                  info = trim(info) // ' with value > ' // trim(ss)
               ENDIF
            ENDIF

            wfmt = "('Check vector data:', A25, ' is in (', e20.10, ',', e20.10, ')', A)"
            write(*,wfmt) varname, vmin, vmax, trim(info)

#if(defined CoLMDEBUG)
            IF (len_trim(info) > 0) THEN
#ifdef USEMPI
               CALL mpi_abort (p_comm_glb, p_err)
#else
               STOP
#endif
            ENDIF
#endif

         ENDIF

      ENDIF

   END SUBROUTINE check_vector_data_real8_3d

   ! ----------
   SUBROUTINE check_vector_data_real8_4d (varname, vdata, spv_in, largevalue)

      USE MOD_Precision
      USE MOD_SPMD_Task
      IMPLICIT NONE

      CHARACTER(len=*), intent(in)   :: varname
      REAL(r8), intent(in)           :: vdata(:,:,:,:)
      REAL(r8), intent(in), optional :: spv_in
      REAL(r8), intent(in), optional :: largevalue

      ! Local variables
      REAL(r8) :: vmin, vmax, spv
      REAL(r8), allocatable :: vmin_all(:), vmax_all(:)
      INTEGER  :: i, j, k, l
      LOGICAL  :: has_nan
      character(len=256) :: wfmt, ss, info

      IF (p_is_worker) THEN

         IF (present(spv_in)) THEN
            spv = spv_in
         ELSE
            spv = -1.0e36_r8
         ENDIF

         IF (any(vdata /= spv)) THEN
            vmin = minval(vdata, mask = vdata /= spv)
            vmax = maxval(vdata, mask = vdata /= spv)
         ELSE
            vmin = spv
            vmax = spv
         ENDIF

         has_nan = .false.
         DO l = 1, size(vdata,4)
            DO k = 1, size(vdata,3)
               DO j = 1, size(vdata,2)
                  DO i = 1, size(vdata,1)
                     has_nan = has_nan .or. isnan(vdata(i,j,k,l))
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

#ifdef USEMPI
         IF (p_iam_worker == p_root) THEN
            allocate (vmin_all (0:p_np_worker-1))
            allocate (vmax_all (0:p_np_worker-1))
            CALL mpi_gather (vmin, 1, MPI_REAL8, vmin_all, 1, MPI_REAL8, p_root, p_comm_worker, p_err)
            CALL mpi_gather (vmax, 1, MPI_REAL8, vmax_all, 1, MPI_REAL8, p_root, p_comm_worker, p_err)
         ELSE
            CALL mpi_gather (vmin, 1, MPI_REAL8, MPI_RNULL_P, 1, MPI_REAL8, p_root, p_comm_worker, p_err)
            CALL mpi_gather (vmax, 1, MPI_REAL8, MPI_RNULL_P, 1, MPI_REAL8, p_root, p_comm_worker, p_err)
         ENDIF


         CALL mpi_allreduce (MPI_IN_PLACE, has_nan, 1, MPI_LOGICAL, MPI_LOR, p_comm_worker, p_err)

         IF (p_iam_worker == p_root) THEN
            IF (any(vmin_all /= spv)) THEN
               vmin = minval(vmin_all, mask = (vmin_all /= spv))
            ELSE
               vmin = spv
            ENDIF

            IF (any(vmax_all /= spv)) THEN
               vmax = maxval(vmax_all, mask = (vmax_all /= spv))
            ELSE
               vmax = spv
            ENDIF

            deallocate (vmin_all)
            deallocate (vmax_all)
         ENDIF
#endif

         IF (p_iam_worker == p_root) THEN

            info = ''

            IF (has_nan) THEN
               info = trim(info) // ' with NAN'
            ENDIF

            IF (present(largevalue)) THEN
               IF (max(abs(vmin),abs(vmax)) > largevalue) THEN
                  write(ss,'(e12.2)') largevalue
                  info = trim(info) // ' with value > ' // trim(ss)
               ENDIF
            ENDIF

            wfmt = "('Check vector data:', A25, ' is in (', e20.10, ',', e20.10, ')', A)"
            write(*,wfmt) varname, vmin, vmax, info

#if(defined CoLMDEBUG)
            IF (len_trim(info) > 0) THEN
#ifdef USEMPI
               CALL mpi_abort (p_comm_glb, p_err)
#else
               STOP
#endif
            ENDIF
#endif

         ENDIF

      ENDIF

   END SUBROUTINE check_vector_data_real8_4d

   ! ----------
   SUBROUTINE check_vector_data_int32_1d (varname, vdata, spv_in)

      USE MOD_Precision
      USE MOD_SPMD_Task
      IMPLICIT NONE

      CHARACTER(len=*), intent(in)  :: varname
      INTEGER, intent(in)           :: vdata(:)
      INTEGER, intent(in), optional :: spv_in

      ! Local variables
      INTEGER :: vmin, vmax
      INTEGER, allocatable :: vmin_all(:), vmax_all(:)
      character(len=256) :: wfmt

      IF (p_is_worker) THEN

         IF (present(spv_in)) THEN
            IF (any(vdata /= spv_in)) THEN
               vmin = minval(vdata, mask = vdata /= spv_in)
               vmax = maxval(vdata, mask = vdata /= spv_in)
            ELSE
               vmin = spv_in
               vmax = spv_in
            ENDIF
         ELSE
            vmin = minval(vdata)
            vmax = maxval(vdata)
         ENDIF

#ifdef USEMPI
         IF (p_iam_worker == p_root) THEN
            allocate (vmin_all (0:p_np_worker-1))
            allocate (vmax_all (0:p_np_worker-1))
            CALL mpi_gather (vmin, 1, MPI_INTEGER, vmin_all, 1, MPI_INTEGER, p_root, p_comm_worker, p_err)
            CALL mpi_gather (vmax, 1, MPI_INTEGER, vmax_all, 1, MPI_INTEGER, p_root, p_comm_worker, p_err)
         ELSE
            CALL mpi_gather (vmin, 1, MPI_INTEGER, MPI_INULL_P, 1, MPI_INTEGER, p_root, p_comm_worker, p_err)
            CALL mpi_gather (vmax, 1, MPI_INTEGER, MPI_INULL_P, 1, MPI_INTEGER, p_root, p_comm_worker, p_err)
         ENDIF

         IF (p_iam_worker == p_root) THEN
            IF (present(spv_in)) THEN
               IF (any(vmin_all /= spv_in)) THEN
                  vmin = minval(vmin_all, mask = (vmin_all /= spv_in))
               ELSE
                  vmin = spv_in
               ENDIF

               IF (any(vmax_all /= spv_in)) THEN
                  vmax = maxval(vmax_all, mask = (vmax_all /= spv_in))
               ELSE
                  vmax = spv_in
               ENDIF
            ELSE
               vmin = minval(vmin_all)
               vmax = maxval(vmax_all)
            ENDIF

            deallocate (vmin_all)
            deallocate (vmax_all)
         ENDIF
#endif

         IF (p_iam_worker == p_root) THEN
            wfmt = "('Check vector data:', A25, ' is in (', I20, ',', I20, ')')"
            write(*,wfmt) varname, vmin, vmax
         ENDIF

      ENDIF

   END SUBROUTINE check_vector_data_int32_1d

#endif

END MODULE MOD_RangeCheck
