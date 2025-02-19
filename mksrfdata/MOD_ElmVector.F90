#include <define.h>

#if (defined UNSTRUCTURED || defined CATCHMENT)
MODULE MOD_ElmVector

!-----------------------------------------------------------------------
! !DESCRIPTION:
!
!    Address of Data associated with land element.
!
!    To output a vector, Data is gathered from worker processes directly to
!    master.  "elm_data_address" stores information on how to reorganize data
!    gathered.  The output data in vector is sorted by global element index.
!
!  Created by Shupeng Zhang, May 2023
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_DataType
   IMPLICIT NONE

   integer :: totalnumelm
   type(pointer_int32_1d), allocatable :: elm_data_address (:)

   integer*8, allocatable :: eindex_glb (:)

CONTAINS

   ! --------
   SUBROUTINE elm_vector_init

   USE MOD_SPMD_Task
   USE MOD_Utils
   USE MOD_Pixelset
   USE MOD_Utils
   USE MOD_UserDefFun
   USE MOD_Mesh
   USE MOD_LandElm
   USE MOD_LandPatch
#ifdef CROP
   USE MOD_LandCrop
#endif
   IMPLICIT NONE

   ! Local Variables
   integer   :: mesg(2), iwork, isrc, ndata
   integer, allocatable :: numelm_worker (:)

   integer :: i, idsp
   integer, allocatable :: vec_worker_dsp (:)

   integer*8, allocatable :: indexelm (:)
   integer,   allocatable :: order    (:)

      IF (p_is_worker) THEN
         CALL elm_patch%build (landelm, landpatch, use_frac = .true.)
      ENDIF

      IF (p_is_worker) THEN
#ifdef USEMPI
         IF (numelm > 0) THEN
            allocate (indexelm (numelm))
            indexelm = landelm%eindex
         ENDIF

         IF (p_iam_worker == p_root) allocate (numelm_worker (0:p_np_worker-1))
         CALL mpi_gather (numelm, 1, MPI_INTEGER, &
            numelm_worker, 1, MPI_INTEGER, p_root, p_comm_worker, p_err)

         IF (p_iam_worker == p_root) THEN
            CALL mpi_send (numelm_worker, p_np_worker, MPI_INTEGER, &
               p_address_master, mpi_tag_size, p_comm_glb, p_err)
         ENDIF

         mesg = (/p_iam_glb, numelm/)
         CALL mpi_send (mesg, 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)
         IF (numelm > 0) THEN
            CALL mpi_send (indexelm, numelm, MPI_INTEGER8, p_address_master, mpi_tag_data, p_comm_glb, p_err)
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
         CALL mpi_recv (numelm_worker, p_np_worker, MPI_INTEGER, p_address_worker(p_root), &
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
            CALL mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
               mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = mesg(1)
            ndata = mesg(2)
            IF (ndata > 0) THEN
               idsp = vec_worker_dsp(p_itis_worker(isrc))
               CALL mpi_recv (eindex_glb(idsp+1:idsp+ndata), ndata, MPI_INTEGER8, isrc, &
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
      CALL mpi_bcast (totalnumelm, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
#endif

      IF (p_is_master) THEN
         allocate (order (totalnumelm))
         order = (/(i, i=1,totalnumelm)/)

         CALL quicksort (totalnumelm, eindex_glb, order)

#ifdef USEMPI
         DO i = 1, totalnumelm
            iwork = findloc_ud(order(i) > vec_worker_dsp, back=.true.) - 1
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

END MODULE MOD_ElmVector
#endif
