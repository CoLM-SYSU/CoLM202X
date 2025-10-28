#include <define.h>

MODULE MOD_WorkerPushData
!--------------------------------------------------------------------------------
! DESCRIPTION:
!--------------------------------------------------------------------------------

   USE MOD_DataType
   IMPLICIT NONE

   type :: worker_pushdata_type

      integer :: num_req_uniq

      integer, allocatable :: addr_single (:)
      integer, allocatable :: addr_multi  (:,:)

      ! data is on the same processor
      integer :: nself
      integer, allocatable :: self_from (:)
      integer, allocatable :: self_to   (:)
#ifdef USEMPI
      ! data is on other processors
      integer, allocatable :: n_to_other   (:)
      integer, allocatable :: n_from_other (:)
      type(pointer_int32_1d), allocatable :: to_other (:)
      type(pointer_int32_1d), allocatable :: other_to (:)
#endif
   CONTAINS
      final :: worker_pushdata_free_mem
   END type worker_pushdata_type

   ! -- public subroutines --
   INTERFACE build_worker_pushdata
      MODULE procedure build_worker_pushdata_single
      MODULE procedure build_worker_pushdata_multi
   END INTERFACE build_worker_pushdata

   INTERFACE worker_push_data
      MODULE procedure worker_push_data_single_real8
      MODULE procedure worker_push_data_single_int32
      MODULE procedure worker_push_data_multi_real8
      MODULE procedure worker_push_data_multi_int32
   END INTERFACE worker_push_data

CONTAINS

   ! ----------
   SUBROUTINE build_worker_pushdata_single (num_me, ids_me, num_req, ids_req, pushdata)

   USE MOD_SPMD_Task
   USE MOD_Utils
   IMPLICIT NONE

   integer, intent(in) :: num_me,  ids_me  (:)
   integer, intent(in) :: num_req, ids_req (:)
   type(worker_pushdata_type), intent(inout) :: pushdata

   ! Local Variables
   integer :: n_req_uniq, iloc, i
   integer, allocatable :: ids_req_uniq (:)

      IF (p_is_worker) THEN

         n_req_uniq = 0

         IF (num_req > 0) THEN
            allocate (ids_req_uniq (num_req))
            DO i = 1, size(ids_req)
               CALL insert_into_sorted_list1 (ids_req(i), n_req_uniq, ids_req_uniq, iloc)
            ENDDO

            allocate (pushdata%addr_single (num_req))
            DO i = 1, size(ids_req)
               pushdata%addr_single(i) = &
                  find_in_sorted_list1 (ids_req(i), n_req_uniq, ids_req_uniq(1:n_req_uniq))
            ENDDO
         ENDIF

         pushdata%num_req_uniq = n_req_uniq

         CALL build_worker_pushdata_uniq (num_me, ids_me, n_req_uniq, ids_req_uniq, pushdata)

         IF (allocated (ids_req_uniq)) deallocate(ids_req_uniq)

      ENDIF

   END SUBROUTINE build_worker_pushdata_single

   ! ----------
   SUBROUTINE build_worker_pushdata_multi (num_me, ids_me, num_req, ids_req, pushdata)

   USE MOD_SPMD_Task
   USE MOD_Utils
   IMPLICIT NONE

   integer, intent(in) :: num_me,  ids_me  (:)
   integer, intent(in) :: num_req, ids_req (:,:)
   type(worker_pushdata_type), intent(inout) :: pushdata

   ! Local Variables
   integer :: n_req_uniq, iloc, i, j
   integer, allocatable :: ids_req_uniq (:)

      IF (p_is_worker) THEN

         n_req_uniq = 0

         IF (num_req > 0) THEN
            allocate (ids_req_uniq (size(ids_req)))
            DO i = 1, size(ids_req,1)
               DO j = 1, size(ids_req,2)
                  IF (ids_req(i,j) > 0) THEN
                     CALL insert_into_sorted_list1 (ids_req(i,j), n_req_uniq, ids_req_uniq, iloc)
                  ENDIF
               ENDDO
            ENDDO

            allocate (pushdata%addr_multi (size(ids_req,1),size(ids_req,2)))
            pushdata%addr_multi(:,:) = 0

            DO i = 1, size(ids_req,1)
               DO j = 1, size(ids_req,2)
                  IF (ids_req(i,j) > 0) THEN
                     pushdata%addr_multi(i,j) = &
                        find_in_sorted_list1 (ids_req(i,j), n_req_uniq, ids_req_uniq(1:n_req_uniq))
                  ENDIF
               ENDDO
            ENDDO
         ENDIF

         pushdata%num_req_uniq = n_req_uniq

         CALL build_worker_pushdata_uniq (num_me, ids_me, n_req_uniq, ids_req_uniq, pushdata)

         IF (allocated (ids_req_uniq)) deallocate(ids_req_uniq)

      ENDIF

   END SUBROUTINE build_worker_pushdata_multi

   ! ----------
   SUBROUTINE build_worker_pushdata_uniq (num_me, ids_me, n_req_uniq, ids_req_uniq, pushdata)

   USE MOD_SPMD_Task
   USE MOD_Utils
   IMPLICIT NONE

   integer, intent(in) :: num_me,     ids_me       (:)
   integer, intent(in) :: n_req_uniq, ids_req_uniq (:)
   type(worker_pushdata_type), intent(inout) :: pushdata

   ! Local Variables
   integer, allocatable :: ids_me_sorted(:), order_ids(:), self_from(:)
#ifdef USEMPI
   integer, allocatable :: other_to(:), ids(:), loc_from_me(:), addrfrom(:)
#endif
   integer :: i, iloc, iworker, jworker, nreq_other


      IF (p_is_worker) THEN

         pushdata%nself = 0

         IF (num_me > 0) THEN
            allocate (ids_me_sorted (num_me))
            ids_me_sorted = ids_me
            allocate (order_ids (num_me))
            order_ids = (/(i, i = 1, num_me)/)

            CALL quicksort (num_me, ids_me_sorted, order_ids)
         ENDIF

         IF (n_req_uniq > 0) THEN
            allocate (self_from (n_req_uniq))
            self_from(:) = -1

            DO i = 1, n_req_uniq
               iloc = find_in_sorted_list1 (ids_req_uniq(i), num_me, ids_me_sorted)
               IF (iloc > 0) THEN
                  self_from(i) = order_ids(iloc)
               ENDIF
            ENDDO

            pushdata%nself = count(self_from > 0)
            IF (pushdata%nself > 0) THEN
               allocate (pushdata%self_from (pushdata%nself))
               pushdata%self_from = pack(self_from, self_from > 0)
               allocate (pushdata%self_to (pushdata%nself))
               pushdata%self_to = pack((/(i,i=1,n_req_uniq)/), self_from > 0)
            ENDIF
         ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_worker, p_err)

         allocate (pushdata%n_to_other   (0:p_np_worker-1))
         allocate (pushdata%to_other     (0:p_np_worker-1))

         allocate (pushdata%n_from_other (0:p_np_worker-1))
         allocate (pushdata%other_to     (0:p_np_worker-1))

         DO iworker = 0, p_np_worker-1

            IF (p_iam_worker == iworker) THEN
               IF (n_req_uniq > 0) THEN
                  nreq_other = count(self_from <= 0)
               ELSE
                  nreq_other = 0
               ENDIF

               IF (nreq_other > 0) THEN
                  allocate (ids (nreq_other))
                  ids = pack(ids_req_uniq, self_from <= 0)
                  allocate (other_to (nreq_other))
                  other_to = pack((/(i,i=1,n_req_uniq)/), self_from <= 0)
               ENDIF
            ENDIF

            CALL mpi_bcast (nreq_other, 1, MPI_INTEGER, iworker, p_comm_worker, p_err)

            IF (nreq_other > 0) THEN

               IF (p_iam_worker /= iworker) THEN
                  allocate (ids (nreq_other))
               ENDIF

               CALL mpi_bcast (ids, nreq_other, MPI_INTEGER, iworker, p_comm_worker, p_err)

               allocate (loc_from_me (nreq_other))
               loc_from_me(:) = -1

               IF (num_me > 0) THEN
                  DO i = 1, nreq_other
                     iloc = find_in_sorted_list1 (ids(i), num_me, ids_me_sorted)
                     IF (iloc > 0) THEN
                        loc_from_me(i) = order_ids(iloc)
                     ENDIF
                  ENDDO
               ENDIF

               pushdata%n_to_other(iworker) = count(loc_from_me > 0)
               IF (pushdata%n_to_other(iworker) > 0) THEN
                  allocate (pushdata%to_other(iworker)%val (pushdata%n_to_other(iworker)))
                  pushdata%to_other(iworker)%val = pack(loc_from_me, loc_from_me > 0)
               ENDIF

               allocate (addrfrom (nreq_other))
               addrfrom(:) = -1
               WHERE (loc_from_me > 0)
                  addrfrom = p_iam_worker
               END WHERE

               CALL mpi_allreduce (addrfrom, nreq_other, MPI_INTEGER, MPI_MAX, p_comm_worker, p_err)

               IF (p_iam_worker == iworker) THEN
                  DO jworker = 0, p_np_worker-1
                     pushdata%n_from_other(jworker) = count(addrfrom == jworker)
                     IF (pushdata%n_from_other(jworker) > 0) THEN
                        allocate (pushdata%other_to(jworker)%val (pushdata%n_from_other(jworker)))
                        pushdata%other_to(jworker)%val = pack(other_to, addrfrom == jworker)
                     ENDIF
                  ENDDO
               ENDIF

               IF (allocated(other_to     )) deallocate(other_to     )
               IF (allocated(ids          )) deallocate(ids          )
               IF (allocated(loc_from_me  )) deallocate(loc_from_me  )
               IF (allocated(addrfrom     )) deallocate(addrfrom     )

            ELSE
               pushdata%n_to_other(iworker) = 0
            ENDIF

         ENDDO

         CALL mpi_barrier (p_comm_worker, p_err)
#endif

         IF (allocated(ids_me_sorted)) deallocate(ids_me_sorted)
         IF (allocated(order_ids    )) deallocate(order_ids    )
         IF (allocated(self_from    )) deallocate(self_from    )

      ENDIF

   END SUBROUTINE build_worker_pushdata_uniq


   ! ----------
   SUBROUTINE worker_push_data_single_real8 (pushdata, vec_send, vec_recv, fillvalue)

   USE MOD_Precision
   USE MOD_SPMD_Task
   IMPLICIT NONE

   type(worker_pushdata_type) :: pushdata

   real(r8), intent(in)    :: vec_send (:)
   real(r8), intent(inout) :: vec_recv (:)
   real(r8), intent(in)    :: fillvalue

   ! Local Variables
   real(r8), allocatable   :: vec_recv_uniq (:)

      IF (p_is_worker) THEN

         IF (pushdata%num_req_uniq > 0) THEN
            allocate (vec_recv_uniq (pushdata%num_req_uniq))
            vec_recv_uniq(:) = fillvalue
         ENDIF

         CALL worker_push_data_uniq_real8 (pushdata, vec_send, vec_recv_uniq)

         IF (pushdata%num_req_uniq > 0) THEN
            vec_recv = vec_recv_uniq(pushdata%addr_single)
            deallocate (vec_recv_uniq)
         ENDIF

      ENDIF

   END SUBROUTINE worker_push_data_single_real8

   ! ----------
   SUBROUTINE worker_push_data_multi_real8 (pushdata, vec_send, vec_recv, fillvalue)

   USE MOD_Precision
   USE MOD_SPMD_Task
   IMPLICIT NONE

   type(worker_pushdata_type) :: pushdata

   real(r8), intent(in)    :: vec_send (:)
   real(r8), intent(inout) :: vec_recv (:,:)
   real(r8), intent(in)    :: fillvalue

   ! Local Variables
   real(r8), allocatable   :: vec_recv_uniq (:)
   integer :: i, j

      IF (p_is_worker) THEN

         IF (pushdata%num_req_uniq > 0) THEN
            allocate (vec_recv_uniq (pushdata%num_req_uniq))
            vec_recv_uniq(:) = fillvalue
         ENDIF

         CALL worker_push_data_uniq_real8 (pushdata, vec_send, vec_recv_uniq)

         IF (pushdata%num_req_uniq > 0) THEN

            vec_recv(:,:) = fillvalue
            DO i = lbound(pushdata%addr_multi,1), ubound(pushdata%addr_multi,1)
               DO j = lbound(pushdata%addr_multi,2), ubound(pushdata%addr_multi,2)
                  IF (pushdata%addr_multi(i,j) > 0) THEN
                     vec_recv(i,j) = vec_recv_uniq(pushdata%addr_multi(i,j))
                  ENDIF
               ENDDO
            ENDDO

            deallocate (vec_recv_uniq)
         ENDIF

      ENDIF

   END SUBROUTINE worker_push_data_multi_real8

   ! ----------
   SUBROUTINE worker_push_data_uniq_real8 (pushdata, vec_send, vec_recv_uniq)

   USE MOD_Precision
   USE MOD_SPMD_Task
   IMPLICIT NONE

   type(worker_pushdata_type) :: pushdata

   real(r8), intent(in)    :: vec_send      (:)
   real(r8), intent(inout) :: vec_recv_uniq (:)

   ! Local Variables
   integer :: ndatasend
   integer,  allocatable :: req_send  (:)
   real(r8), allocatable :: sendcache (:)

   integer :: ndatarecv
   integer,  allocatable :: req_recv  (:)
   real(r8), allocatable :: recvcache (:)

   integer :: iworker, iproc, istt, iend

      IF (p_is_worker) THEN

         IF (pushdata%nself > 0) THEN
            vec_recv_uniq(pushdata%self_to) = vec_send(pushdata%self_from)
         ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_worker, p_err)

         ndatasend = sum(pushdata%n_to_other)
         IF (ndatasend > 0) THEN
            allocate (sendcache(ndatasend))
            allocate (req_send (count(pushdata%n_to_other > 0)))

            iproc = 0
            iend  = 0
            DO iworker = 0, p_np_worker-1
               IF (pushdata%n_to_other(iworker) > 0) THEN
                  iproc = iproc + 1
                  istt  = iend + 1
                  iend  = iend + pushdata%n_to_other(iworker)

                  sendcache(istt:iend) = vec_send(pushdata%to_other(iworker)%val)

                  CALL mpi_isend(sendcache(istt:iend), pushdata%n_to_other(iworker), MPI_REAL8, &
                     iworker, 101, p_comm_worker, req_send(iproc), p_err)
               ENDIF
            ENDDO
         ENDIF

         ndatarecv = sum(pushdata%n_from_other)
         IF (ndatarecv > 0) THEN

            allocate (recvcache(ndatarecv))
            allocate (req_recv (count(pushdata%n_from_other > 0)))

            iproc = 0
            iend  = 0
            DO iworker = 0, p_np_worker-1
               IF (pushdata%n_from_other(iworker) > 0) THEN
                  iproc = iproc + 1
                  istt  = iend + 1
                  iend  = iend + pushdata%n_from_other(iworker)

                  CALL mpi_irecv(recvcache(istt:iend), pushdata%n_from_other(iworker), MPI_REAL8, &
                     iworker, 101, p_comm_worker, req_recv(iproc), p_err)
               ENDIF
            ENDDO
         ENDIF

         IF (ndatarecv > 0) THEN

            CALL mpi_waitall(size(req_recv), req_recv, MPI_STATUSES_IGNORE, p_err)

            iend = 0
            DO iworker = 0, p_np_worker-1
               IF (pushdata%n_from_other(iworker) > 0) THEN
                  istt = iend + 1
                  iend = iend + pushdata%n_from_other(iworker)
                  vec_recv_uniq(pushdata%other_to(iworker)%val) = recvcache(istt:iend)
               ENDIF
            ENDDO
         ENDIF

         IF (ndatasend > 0) THEN
            CALL mpi_waitall(size(req_send), req_send, MPI_STATUSES_IGNORE, p_err)
         ENDIF

         IF (allocated(req_send )) deallocate(req_send )
         IF (allocated(sendcache)) deallocate(sendcache)
         IF (allocated(req_recv )) deallocate(req_recv )
         IF (allocated(recvcache)) deallocate(recvcache)

         CALL mpi_barrier (p_comm_worker, p_err)
#endif

      ENDIF

   END SUBROUTINE worker_push_data_uniq_real8

   ! ----------
   SUBROUTINE worker_push_data_single_int32 (pushdata, vec_send, vec_recv, fillvalue)

   USE MOD_Precision
   USE MOD_SPMD_Task
   IMPLICIT NONE

   type(worker_pushdata_type) :: pushdata

   integer, intent(in)    :: vec_send (:)
   integer, intent(inout) :: vec_recv (:)
   integer, intent(in)    :: fillvalue

   ! Local Variables
   integer, allocatable   :: vec_recv_uniq (:)

      IF (p_is_worker) THEN

         IF (pushdata%num_req_uniq > 0) THEN
            allocate (vec_recv_uniq (pushdata%num_req_uniq))
            vec_recv_uniq(:) = fillvalue
         ENDIF

         CALL worker_push_data_uniq_int32 (pushdata, vec_send, vec_recv_uniq)

         IF (pushdata%num_req_uniq > 0) THEN
            vec_recv = vec_recv_uniq(pushdata%addr_single)
            deallocate (vec_recv_uniq)
         ENDIF

      ENDIF

   END SUBROUTINE worker_push_data_single_int32

   ! ----------
   SUBROUTINE worker_push_data_multi_int32 (pushdata, vec_send, vec_recv, fillvalue)

   USE MOD_Precision
   USE MOD_SPMD_Task
   IMPLICIT NONE

   type(worker_pushdata_type) :: pushdata

   integer, intent(in)    :: vec_send (:)
   integer, intent(inout) :: vec_recv (:,:)
   integer, intent(in)    :: fillvalue

   ! Local Variables
   integer, allocatable   :: vec_recv_uniq (:)
   integer :: i, j

      IF (p_is_worker) THEN

         IF (pushdata%num_req_uniq > 0) THEN
            allocate (vec_recv_uniq (pushdata%num_req_uniq))
            vec_recv_uniq(:) = fillvalue
         ENDIF

         CALL worker_push_data_uniq_int32 (pushdata, vec_send, vec_recv_uniq)

         IF (pushdata%num_req_uniq > 0) THEN

            vec_recv(:,:) = fillvalue
            DO i = lbound(pushdata%addr_multi,1), ubound(pushdata%addr_multi,1)
               DO j = lbound(pushdata%addr_multi,2), ubound(pushdata%addr_multi,2)
                  IF (pushdata%addr_multi(i,j) > 0) THEN
                     vec_recv(i,j) = vec_recv_uniq(pushdata%addr_multi(i,j))
                  ENDIF
               ENDDO
            ENDDO

            deallocate (vec_recv_uniq)
         ENDIF

      ENDIF

   END SUBROUTINE worker_push_data_multi_int32

   ! ----------
   SUBROUTINE worker_push_data_uniq_int32 (pushdata, vec_send, vec_recv_uniq)

   USE MOD_Precision
   USE MOD_SPMD_Task
   IMPLICIT NONE

   type(worker_pushdata_type) :: pushdata

   integer, intent(in)    :: vec_send      (:)
   integer, intent(inout) :: vec_recv_uniq (:)

   ! Local Variables
   integer :: ndatasend
   integer, allocatable :: req_send  (:)
   integer, allocatable :: sendcache (:)

   integer :: ndatarecv
   integer, allocatable :: req_recv  (:)
   integer, allocatable :: recvcache (:)

   integer :: iworker, iproc, istt, iend

      IF (p_is_worker) THEN

         IF (pushdata%nself > 0) THEN
            vec_recv_uniq(pushdata%self_to) = vec_send(pushdata%self_from)
         ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_worker, p_err)

         ndatasend = sum(pushdata%n_to_other)
         IF (ndatasend > 0) THEN
            allocate (sendcache(ndatasend))
            allocate (req_send (count(pushdata%n_to_other > 0)))

            iproc = 0
            iend  = 0
            DO iworker = 0, p_np_worker-1
               IF (pushdata%n_to_other(iworker) > 0) THEN
                  iproc = iproc + 1
                  istt  = iend + 1
                  iend  = iend + pushdata%n_to_other(iworker)

                  sendcache(istt:iend) = vec_send(pushdata%to_other(iworker)%val)

                  CALL mpi_isend(sendcache(istt:iend), pushdata%n_to_other(iworker), MPI_INTEGER, &
                     iworker, 101, p_comm_worker, req_send(iproc), p_err)
               ENDIF
            ENDDO
         ENDIF

         ndatarecv = sum(pushdata%n_from_other)
         IF (ndatarecv > 0) THEN

            allocate (recvcache(ndatarecv))
            allocate (req_recv (count(pushdata%n_from_other > 0)))

            iproc = 0
            iend  = 0
            DO iworker = 0, p_np_worker-1
               IF (pushdata%n_from_other(iworker) > 0) THEN
                  iproc = iproc + 1
                  istt  = iend + 1
                  iend  = iend + pushdata%n_from_other(iworker)

                  CALL mpi_irecv(recvcache(istt:iend), pushdata%n_from_other(iworker), MPI_INTEGER, &
                     iworker, 101, p_comm_worker, req_recv(iproc), p_err)
               ENDIF
            ENDDO
         ENDIF

         IF (ndatarecv > 0) THEN

            CALL mpi_waitall(size(req_recv), req_recv, MPI_STATUSES_IGNORE, p_err)

            iend = 0
            DO iworker = 0, p_np_worker-1
               IF (pushdata%n_from_other(iworker) > 0) THEN
                  istt = iend + 1
                  iend = iend + pushdata%n_from_other(iworker)
                  vec_recv_uniq(pushdata%other_to(iworker)%val) = recvcache(istt:iend)
               ENDIF
            ENDDO
         ENDIF

         IF (ndatasend > 0) THEN
            CALL mpi_waitall(size(req_send), req_send, MPI_STATUSES_IGNORE, p_err)
         ENDIF

         IF (allocated(req_send )) deallocate(req_send)
         IF (allocated(sendcache)) deallocate(sendcache)
         IF (allocated(req_recv )) deallocate(req_recv)
         IF (allocated(recvcache)) deallocate(recvcache)

         CALL mpi_barrier (p_comm_worker, p_err)
#endif

      ENDIF

   END SUBROUTINE worker_push_data_uniq_int32

   ! ---------
   SUBROUTINE worker_pushdata_free_mem (this)

   IMPLICIT NONE
   type(worker_pushdata_type) :: this

      IF (allocated(this%addr_single )) deallocate(this%addr_single )
      IF (allocated(this%addr_multi  )) deallocate(this%addr_multi  )
      IF (allocated(this%self_from   )) deallocate(this%self_from   )
      IF (allocated(this%self_to     )) deallocate(this%self_to     )
      IF (allocated(this%n_to_other  )) deallocate(this%n_to_other  )
      IF (allocated(this%n_from_other)) deallocate(this%n_from_other)
      IF (allocated(this%to_other    )) deallocate(this%to_other    )
      IF (allocated(this%other_to    )) deallocate(this%other_to    )

   END SUBROUTINE worker_pushdata_free_mem

END MODULE MOD_WorkerPushData
