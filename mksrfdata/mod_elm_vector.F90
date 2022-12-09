#include <define.h>

#if (defined UNSTRUCTURED || defined CATCHMENT) 
MODULE mod_elm_vector

   USE precision
   USE mod_data_type
   IMPLICIT NONE
   
   INTEGER :: totalnumelm
   TYPE(pointer_int32_1d), allocatable :: elm_data_address (:)

   INTEGER, allocatable :: eindex_glb (:)
   
CONTAINS
   
   ! --------
   SUBROUTINE elm_vector_init 

      USE spmd_task
      USE mod_utils
      USE mod_pixelset
      USE mod_utils
      USE mod_mesh
      USE mod_landelm
      USE mod_landpatch
      IMPLICIT NONE

      ! Local Variables
      INTEGER   :: mesg(2), iwork, isrc, ndata
      INTEGER, allocatable :: numelm_worker (:)

      INTEGER :: i, idsp
      INTEGER, allocatable :: vec_worker_dsp (:)
      INTEGER, allocatable :: indexelm (:)
      INTEGER, allocatable :: order    (:)
      
      IF (p_is_worker) THEN
         CALL elm_patch%build (landelm, landpatch, use_frac = .true.)
      ENDIF

      IF (p_is_worker) THEN
#ifdef USEMPI
         IF (numelm > 0) THEN
            allocate (indexelm (numelm))
            indexelm = landelm%eindex 
         ENDIF
         
         IF (p_iam_worker == 0) allocate (numelm_worker (0:p_np_worker-1))
         CALL mpi_gather (numelm, 1, MPI_INTEGER, &
            numelm_worker, 1, MPI_INTEGER, p_root, p_comm_worker, p_err)

         IF (p_iam_worker == 0) THEN
            call mpi_send (numelm_worker, p_np_worker, MPI_INTEGER, &
               p_root, mpi_tag_size, p_comm_glb, p_err) 
         ENDIF

         mesg = (/p_iam_glb, numelm/)
         call mpi_send (mesg, 2, MPI_INTEGER, p_root, mpi_tag_mesg, p_comm_glb, p_err) 
         IF (numelm > 0) THEN
            call mpi_send (indexelm, numelm, MPI_INTEGER, p_root, mpi_tag_data, p_comm_glb, p_err) 
         ENDIF
#else
         IF (numelm > 0) THEN
            allocate (eindex_glb (numelm))
            eindex_glb = landelm%eindex 
         ENDIF
#endif
      ENDIF

      IF (p_is_master) THEN
#ifdef USEMPI
         allocate (numelm_worker (0:p_np_worker-1))
         call mpi_recv (numelm_worker, p_np_worker, MPI_INTEGER, p_address_worker(0), &
            mpi_tag_size, p_comm_glb, p_stat, p_err)

         allocate (vec_worker_dsp (0:p_np_worker-1))
         vec_worker_dsp(0) = 0
         DO iwork = 1, p_np_worker-1
            vec_worker_dsp(iwork) = vec_worker_dsp(iwork-1) + numelm_worker(iwork-1)
         ENDDO

         totalnumelm = sum(numelm_worker)
         
         allocate (eindex_glb (totalnumelm))

         allocate (elm_data_address(0:p_np_worker-1))
         DO iwork = 0, p_np_worker-1
            IF (numelm_worker(iwork) > 0) THEN
               allocate (elm_data_address(iwork)%val (numelm_worker(iwork)))
            ENDIF
         ENDDO
         
         DO iwork = 0, p_np_worker-1
            call mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
               mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = mesg(1)
            ndata = mesg(2)
            IF (ndata > 0) THEN
               idsp = vec_worker_dsp(p_itis_worker(isrc))
               call mpi_recv (eindex_glb(idsp+1:idsp+ndata), ndata, MPI_INTEGER, isrc, &
                  mpi_tag_data, p_comm_glb, p_stat, p_err)
            ENDIF
         ENDDO
#else
         totalnumelm = numelm
         allocate (elm_data_address(0:0))
         allocate (elm_data_address(0)%val (totalnumelm))
#endif
      ENDIF 
      
#ifdef USEMPI
      CALL mpi_bcast (totalnumelm, 1, MPI_INTEGER, p_root, p_comm_glb, p_err)
#endif

      IF (p_is_master) THEN
         allocate (order (totalnumelm))
         order = (/(i, i=1,totalnumelm)/)

         CALL quicksort (totalnumelm, eindex_glb, order)

#ifdef USEMPI
         DO i = 1, totalnumelm
            iwork = findloc(order(i) > vec_worker_dsp, .true., dim=1, back=.true.) - 1
            elm_data_address(iwork)%val(order(i)-vec_worker_dsp(iwork)) = i
         ENDDO
#else
         elm_data_address(0)%val (order) = (/(i, i=1,totalnumelm)/)
#endif
      ENDIF

      IF (allocated(numelm_worker))  deallocate(numelm_worker)
      IF (allocated(vec_worker_dsp)) deallocate(vec_worker_dsp)
      IF (allocated(indexelm))       deallocate(indexelm)
      IF (allocated(order))          deallocate(order)

   END SUBROUTINE elm_vector_init 

   ! ----------
   SUBROUTINE elm_vector_final ()

      IMPLICIT NONE

      IF (allocated(elm_data_address)) deallocate (elm_data_address)
      IF (allocated(eindex_glb))       deallocate (eindex_glb)
      
   END SUBROUTINE elm_vector_final

END MODULE mod_elm_vector
#endif
