#include <define.h>

MODULE MOD_WorkerPushData
!--------------------------------------------------------------------------------
! DESCRIPTION:
!--------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_DataType
   USE MOD_SPMD_Task
   USE MOD_Utils
   IMPLICIT NONE

   ! -- Data Type : push data between workers --
   type :: worker_pushdata_type

      integer :: num_req_uniq

      integer,  allocatable :: addr_single (:)

      integer,  allocatable :: addr_multi  (:,:)
      real(r8), allocatable :: area_multi  (:,:)
      real(r8), allocatable :: sum_area    (:)

      ! data is on the same processor
      integer :: nself
      integer,  allocatable :: self_from (:)
      integer,  allocatable :: self_to   (:)
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


   ! -- Data Type : remap data on workers --
   type :: worker_remapdata_type

      integer :: num_grid
      integer, allocatable :: ilon_me (:)
      integer, allocatable :: ilat_me (:)
      integer, allocatable :: ids_me  (:)

      integer :: npset
      real(r8), allocatable :: sum_area (:)
      integer,  allocatable :: npart    (:)
      type(pointer_int32_1d), allocatable :: part_to  (:) !
      type(pointer_real8_1d), allocatable :: areapart (:) ! intersection area

   CONTAINS
      final :: worker_remapdata_free_mem
   END type worker_remapdata_type

   ! -- public subroutines --
   INTERFACE build_worker_pushdata
      MODULE procedure build_worker_pushdata_id2id_single
      MODULE procedure build_worker_pushdata_id2id_multi
   END INTERFACE build_worker_pushdata

   PUBLIC :: build_worker_remapdata

   INTERFACE worker_push_data
      MODULE procedure worker_push_data_single_real8
      MODULE procedure worker_push_data_single_int32
      MODULE procedure worker_push_data_multi_real8
   END INTERFACE worker_push_data

   INTERFACE worker_remap_data_pset2grid
      MODULE procedure worker_remap_data_pset2grid_real8
   END INTERFACE worker_remap_data_pset2grid

   INTERFACE worker_remap_data_grid2pset
      MODULE procedure worker_remap_data_grid2pset_real8
   END INTERFACE worker_remap_data_grid2pset

CONTAINS

   ! ----------
   SUBROUTINE build_worker_pushdata_id2id_single (num_me, ids_me, num_req, ids_req, pushdata)

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

         CALL build_worker_pushdata_id2id_uniq ( &
            num_me, ids_me, n_req_uniq, ids_req_uniq(1:n_req_uniq), pushdata)

         IF (allocated (ids_req_uniq)) deallocate(ids_req_uniq)

      ENDIF

   END SUBROUTINE build_worker_pushdata_id2id_single

   ! ----------
   SUBROUTINE build_worker_pushdata_id2id_multi ( &
         num_me, ids_me, num_req, ids_req, area_req, pushdata)

   IMPLICIT NONE

   integer,  intent(in) :: num_me,  ids_me  (:)
   integer,  intent(in) :: num_req, ids_req (:,:)
   real(r8), intent(in) :: area_req(:,:)
   type(worker_pushdata_type), intent(inout) :: pushdata

   ! Local Variables
   integer :: ndim1, n_req_uniq, iloc, i, j, iworker
   integer, allocatable :: ids_req_uniq (:)
   logical, allocatable :: id_found     (:)

      IF (p_is_worker) THEN

         n_req_uniq = 0

         IF (num_req > 0) THEN

            ndim1 = size(ids_req,1)

            allocate (ids_req_uniq (ndim1*num_req))
            DO j = 1, num_req
               DO i = 1, ndim1
                  CALL insert_into_sorted_list1 (ids_req(i,j), n_req_uniq, ids_req_uniq, iloc)
               ENDDO
            ENDDO

            allocate (pushdata%addr_multi (ndim1,num_req))

            DO j = 1, num_req
               DO i = 1, ndim1
                  pushdata%addr_multi(i,j) = &
                     find_in_sorted_list1 (ids_req(i,j), n_req_uniq, ids_req_uniq(1:n_req_uniq))
               ENDDO
            ENDDO
         ENDIF

         pushdata%num_req_uniq = n_req_uniq

         CALL build_worker_pushdata_id2id_uniq ( &
            num_me, ids_me, n_req_uniq, ids_req_uniq(1:n_req_uniq), pushdata)

         IF (num_req > 0) THEN
            allocate (pushdata%area_multi (ndim1,num_req))
            allocate (pushdata%sum_area   (num_req))

            pushdata%area_multi = area_req

            WHERE ((area_req <= 0.) .or. (ids_req <= 0))
               pushdata%area_multi = 0.
            END WHERE

            allocate (id_found (n_req_uniq))
            id_found(:) = .false.

            IF (pushdata%nself > 0) id_found(pushdata%self_to) = .true.
#ifdef USEMPI
            DO iworker = 0, p_np_worker-1
               IF (pushdata%n_from_other(iworker) > 0) THEN
                  id_found(pushdata%other_to(iworker)%val) = .true.
               ENDIF
            ENDDO
#endif

            DO j = 1, num_req
               DO i = 1, ndim1
                  IF (.not. id_found(pushdata%addr_multi(i,j))) then
                     pushdata%area_multi(i,j) = 0.
                  ENDIF
               ENDDO
            ENDDO

            pushdata%sum_area = sum(pushdata%area_multi, dim = 1)

            deallocate (id_found)
         ENDIF

         IF (allocated (ids_req_uniq)) deallocate(ids_req_uniq)

      ENDIF

   END SUBROUTINE build_worker_pushdata_id2id_multi

   ! ----------
   SUBROUTINE build_worker_pushdata_id2id_uniq (num_me, ids_me, n_req_uniq, ids_req_uniq, pushdata)

   IMPLICIT NONE

   integer, intent(in) :: num_me,     ids_me       (:)
   integer, intent(in) :: n_req_uniq, ids_req_uniq (:)
   type(worker_pushdata_type), intent(inout) :: pushdata

   ! Local Variables
   integer, allocatable :: ids_me_sorted(:), order_ids(:), self_from(:)
#ifdef USEMPI
   integer, allocatable :: ids(:), loc_from_me(:)
#endif
   integer :: i, iloc, iworker, jworker, n_req_other


      IF (p_is_worker) THEN

         IF (num_me > 0) THEN
            allocate (ids_me_sorted (num_me));  ids_me_sorted = ids_me
            allocate (order_ids     (num_me));  order_ids = (/(i, i = 1, num_me)/)

            CALL quicksort (num_me, ids_me_sorted, order_ids)
         ENDIF

         pushdata%nself = 0

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
               allocate (pushdata%self_to   (pushdata%nself))
               pushdata%self_from = pack(self_from, self_from > 0)
               pushdata%self_to   = pack((/(i,i=1,n_req_uniq)/), self_from > 0)
            ENDIF
         ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_worker, p_err)

         allocate (pushdata%n_to_other   (0:p_np_worker-1))
         allocate (pushdata%to_other     (0:p_np_worker-1))

         allocate (pushdata%n_from_other (0:p_np_worker-1))
         allocate (pushdata%other_to     (0:p_np_worker-1))

         DO iworker = 0, p_np_worker-1

            IF (p_iam_worker == iworker) n_req_other = n_req_uniq
            CALL mpi_bcast (n_req_other, 1, MPI_INTEGER, iworker, p_comm_worker, p_err)

            IF (n_req_other > 0) THEN

               allocate (ids         (n_req_other))
               allocate (loc_from_me (n_req_other))

               IF (p_iam_worker == iworker) ids = ids_req_uniq
               CALL mpi_bcast (ids, n_req_other, MPI_INTEGER, iworker, p_comm_worker, p_err)

               IF (p_iam_worker /= iworker) THEN

                  loc_from_me(:) = -1

                  IF (num_me > 0) THEN
                     DO i = 1, n_req_other
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

                  CALL mpi_send (loc_from_me, n_req_other, MPI_INTEGER, &
                     iworker, mpi_tag_data, p_comm_worker, p_err)

               ELSE

                  pushdata%n_to_other  (iworker) = 0
                  pushdata%n_from_other(iworker) = 0

                  DO jworker = 0, p_np_worker-1
                     IF (jworker /= iworker) THEN

                        CALL mpi_recv (loc_from_me, n_req_other, MPI_INTEGER, &
                           jworker, mpi_tag_data, p_comm_worker, p_stat, p_err)

                        pushdata%n_from_other(jworker) = count(loc_from_me > 0)
                        IF (pushdata%n_from_other(jworker) > 0) THEN
                           allocate (pushdata%other_to(jworker)%val (pushdata%n_from_other(jworker)))
                           pushdata%other_to(jworker)%val = pack((/(i,i=1,n_req_other)/), loc_from_me > 0)
                        ENDIF

                     ENDIF
                  ENDDO
               ENDIF

               deallocate(ids        )
               deallocate(loc_from_me)

            ELSE
               IF (p_iam_worker == iworker) THEN
                  pushdata%n_from_other(:) = 0
               ENDIF
               pushdata%n_to_other(iworker) = 0
            ENDIF

         ENDDO

         CALL mpi_barrier (p_comm_worker, p_err)
#endif

         IF (allocated(ids_me_sorted)) deallocate(ids_me_sorted)
         IF (allocated(order_ids    )) deallocate(order_ids    )
         IF (allocated(self_from    )) deallocate(self_from    )

      ENDIF

   END SUBROUTINE build_worker_pushdata_id2id_uniq

   ! ----------
   SUBROUTINE build_worker_remapdata (pixelset, grid, remapdata)

   USE MOD_Grid
   USE MOD_Pixelset
   USE MOD_SpatialMapping
   IMPLICIT NONE

   type(pixelset_type), intent(in) :: pixelset
   type(grid_type),     intent(in) :: grid

   type(worker_remapdata_type), intent(inout) :: remapdata

   ! Local Variables
   type(spatial_mapping_type) :: mapping
   integer, allocatable :: ilon_me(:), ilat_me(:)
   integer :: ngrid, iproc, ig, iloc, iset, ipart


      CALL mapping%build_arealweighted (grid, pixelset)

      IF (p_is_worker) THEN

         ngrid = 0
         DO iproc = 0, p_np_io-1
            ngrid = ngrid + mapping%glist(iproc)%ng
         ENDDO

         IF (ngrid > 0) THEN
            allocate (ilon_me (ngrid))
            allocate (ilat_me (ngrid))
         ENDIF

         ngrid = 0
         DO iproc = 0, p_np_io-1
            DO ig = 1, mapping%glist(iproc)%ng
               CALL insert_into_sorted_list2 ( &
                  mapping%glist(iproc)%ilon(ig), mapping%glist(iproc)%ilat(ig), &
                  ngrid, ilon_me, ilat_me, iloc)
            ENDDO
         ENDDO

         remapdata%num_grid = ngrid

         IF (ngrid > 0) THEN
            allocate (remapdata%ilon_me (ngrid))
            allocate (remapdata%ilat_me (ngrid))
            allocate (remapdata%ids_me  (ngrid))
            remapdata%ilon_me = ilon_me(1:ngrid)
            remapdata%ilat_me = ilat_me(1:ngrid)
            remapdata%ids_me  = (ilat_me(1:ngrid)-1) * grid%nlon + ilon_me(1:ngrid)
         ENDIF

         remapdata%npset = mapping%npset

         IF (remapdata%npset > 0) THEN

            allocate (remapdata%sum_area (remapdata%npset))
            allocate (remapdata%npart    (remapdata%npset))
            allocate (remapdata%part_to  (remapdata%npset))
            allocate (remapdata%areapart (remapdata%npset))

            remapdata%npart = mapping%npart

            DO iset = 1, remapdata%npset
               IF (remapdata%npart(iset) > 0) THEN
                  allocate (remapdata%part_to(iset)%val  (remapdata%npart(iset)))
                  allocate (remapdata%areapart(iset)%val (remapdata%npart(iset)))
               ENDIF

               DO ipart = 1, remapdata%npart(iset)
                  iproc = mapping%address(iset)%val(1,ipart)
                  iloc  = mapping%address(iset)%val(2,ipart)
                  remapdata%part_to(iset)%val(ipart)  = find_in_sorted_list2 ( &
                     mapping%glist(iproc)%ilon(iloc), mapping%glist(iproc)%ilat(iloc), ngrid, ilon_me, ilat_me)
                  ! from km^2 to m^2
                  remapdata%areapart(iset)%val(ipart) = mapping%areapart(iset)%val(ipart) * 1.e6
               ENDDO

               IF (remapdata%npart(iset) > 0) THEN
                  remapdata%sum_area(iset) = sum(remapdata%areapart(iset)%val)
               ENDIF
            ENDDO
         ENDIF

         IF (allocated(ilon_me)) deallocate(ilon_me)
         IF (allocated(ilat_me)) deallocate(ilat_me)

      ENDIF

   END SUBROUTINE build_worker_remapdata


   ! ----------
   SUBROUTINE worker_push_data_uniq ( pushdata,      &
         vec_send_real8, vec_recv_real8, rfillvalue, &
         vec_send_int32, vec_recv_int32, ifillvalue)

   IMPLICIT NONE

   type(worker_pushdata_type), intent(in) :: pushdata

   real(r8), intent(in)   , optional :: vec_send_real8 (:)
   real(r8), intent(inout), optional :: vec_recv_real8 (:)
   real(r8), intent(in)   , optional :: rfillvalue

   integer,  intent(in)   , optional :: vec_send_int32 (:)
   integer,  intent(inout), optional :: vec_recv_int32 (:)
   integer,  intent(in)   , optional :: ifillvalue

   ! Local Variables
   logical :: has_real8, has_int32

   integer :: ndatasend
   integer,  allocatable :: req_send   (:)
   real(r8), allocatable :: rsendcache (:)
   integer,  allocatable :: isendcache (:)

   integer :: ndatarecv
   integer,  allocatable :: req_recv   (:)
   real(r8), allocatable :: rrecvcache (:)
   integer,  allocatable :: irecvcache (:)

   integer :: iworker, iproc, idsp, istt, iend, i, i_to


      IF (p_is_worker) THEN

         has_real8 = present(vec_recv_real8) .and. present(vec_send_real8) .and. present(rfillvalue)
         has_int32 = present(vec_recv_int32) .and. present(vec_send_int32) .and. present(ifillvalue)

         IF ((has_real8 .and. has_int32) .or. (.not. (has_real8 .or. has_int32))) THEN
            CALL CoLM_Stop ('Worker Push Data: Push one data once.')
         ENDIF

         IF (pushdata%num_req_uniq > 0) THEN
            IF (has_real8) vec_recv_real8 = rfillvalue
            IF (has_int32) vec_recv_int32 = ifillvalue
         ENDIF

         IF (pushdata%nself > 0) THEN
            IF (has_real8) THEN
               vec_recv_real8(pushdata%self_to) = vec_send_real8(pushdata%self_from)
            ENDIF
            IF (has_int32) THEN
               vec_recv_int32(pushdata%self_to) = vec_send_int32(pushdata%self_from)
            ENDIF
         ENDIF

#ifdef USEMPI
         ndatasend = sum(pushdata%n_to_other)
         IF (ndatasend > 0) THEN

            IF (has_real8) allocate (rsendcache(ndatasend))
            IF (has_int32) allocate (isendcache(ndatasend))

            allocate (req_send (count(pushdata%n_to_other > 0)))

            iproc = 0
            idsp  = 0
            DO iworker = 0, p_np_worker-1
               IF (pushdata%n_to_other(iworker) > 0) THEN
                  iproc = iproc + 1
                  istt  = idsp + 1
                  iend  = idsp + pushdata%n_to_other(iworker)

                  IF (has_real8) THEN
                     rsendcache(istt:iend) = vec_send_real8(pushdata%to_other(iworker)%val)
                     CALL mpi_isend(rsendcache(istt:iend), pushdata%n_to_other(iworker), MPI_REAL8, &
                        iworker, 101, p_comm_worker, req_send(iproc), p_err)
                  ENDIF

                  IF (has_int32) THEN
                     isendcache(istt:iend) = vec_send_int32(pushdata%to_other(iworker)%val)
                     CALL mpi_isend(isendcache(istt:iend), pushdata%n_to_other(iworker), MPI_INTEGER, &
                        iworker, 101, p_comm_worker, req_send(iproc), p_err)
                  ENDIF

                  idsp = iend
               ENDIF
            ENDDO
         ENDIF

         ndatarecv = sum(pushdata%n_from_other)
         IF (ndatarecv > 0) THEN

            IF (has_real8) allocate (rrecvcache(ndatarecv))
            IF (has_int32) allocate (irecvcache(ndatarecv))

            allocate (req_recv (count(pushdata%n_from_other > 0)))

            iproc = 0
            idsp  = 0
            DO iworker = 0, p_np_worker-1
               IF (pushdata%n_from_other(iworker) > 0) THEN
                  iproc = iproc + 1
                  istt  = idsp + 1
                  iend  = idsp + pushdata%n_from_other(iworker)

                  IF (has_real8) THEN
                     CALL mpi_irecv(rrecvcache(istt:iend), pushdata%n_from_other(iworker), MPI_REAL8, &
                        iworker, 101, p_comm_worker, req_recv(iproc), p_err)
                  ENDIF

                  IF (has_int32) THEN
                     CALL mpi_irecv(irecvcache(istt:iend), pushdata%n_from_other(iworker), MPI_INTEGER, &
                        iworker, 101, p_comm_worker, req_recv(iproc), p_err)
                  ENDIF

                  idsp = iend
               ENDIF
            ENDDO
         ENDIF

         IF (ndatarecv > 0) THEN

            CALL mpi_waitall(size(req_recv), req_recv, MPI_STATUSES_IGNORE, p_err)

            idsp = 0
            DO iworker = 0, p_np_worker-1

               DO i = 1, pushdata%n_from_other(iworker)

                  IF (has_real8) THEN
                     IF (rrecvcache(idsp+i) /= rfillvalue) THEN
                        i_to = pushdata%other_to(iworker)%val(i)
                        IF (vec_recv_real8(i_to) == rfillvalue) THEN
                           vec_recv_real8(i_to) = rrecvcache(idsp+i)
                        ELSE
                           vec_recv_real8(i_to) = vec_recv_real8(i_to) + rrecvcache(idsp+i)
                        ENDIF
                     ENDIF
                  ENDIF

                  IF (has_int32) THEN
                     IF (irecvcache(idsp+i) /= ifillvalue) THEN
                        i_to = pushdata%other_to(iworker)%val(i)
                        IF (vec_recv_int32(i_to) == ifillvalue) THEN
                           vec_recv_int32(i_to) = irecvcache(idsp+i)
                        ELSE
                           vec_recv_int32(i_to) = vec_recv_int32(i_to) + irecvcache(idsp+i)
                        ENDIF
                     ENDIF
                  ENDIF

               ENDDO

               idsp = idsp + pushdata%n_from_other(iworker)
            ENDDO

         ENDIF

         IF (ndatasend > 0) THEN
            CALL mpi_waitall(size(req_send), req_send, MPI_STATUSES_IGNORE, p_err)
         ENDIF

         IF (allocated(req_send  )) deallocate(req_send  )
         IF (allocated(rsendcache)) deallocate(rsendcache)
         IF (allocated(isendcache)) deallocate(isendcache)
         IF (allocated(req_recv  )) deallocate(req_recv  )
         IF (allocated(rrecvcache)) deallocate(rrecvcache)
         IF (allocated(irecvcache)) deallocate(irecvcache)
#endif

      ENDIF

   END SUBROUTINE worker_push_data_uniq


   ! ----------
   SUBROUTINE worker_push_data_single_real8 (pushdata, vec_send, vec_recv, fillvalue)

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

         CALL worker_push_data_uniq (pushdata, &
            vec_send_real8 = vec_send, vec_recv_real8 = vec_recv_uniq, rfillvalue = fillvalue)

         IF (pushdata%num_req_uniq > 0) THEN
            vec_recv = vec_recv_uniq(pushdata%addr_single)
            deallocate (vec_recv_uniq)
         ENDIF

      ENDIF

   END SUBROUTINE worker_push_data_single_real8

   ! ----------
   SUBROUTINE worker_push_data_multi_real8 (pushdata, vec_send, vec_recv, fillvalue, mode)

   IMPLICIT NONE

   type(worker_pushdata_type) :: pushdata

   real(r8), intent(in)    :: vec_send (:)
   real(r8), intent(inout) :: vec_recv (:)
   real(r8), intent(in)    :: fillvalue

   character(len=*), intent(in) :: mode

   ! Local Variables
   real(r8), allocatable :: vec_recv_uniq (:)
   integer  :: i, j
   real(r8) :: val, sumarea

      IF (p_is_worker) THEN

         IF (pushdata%num_req_uniq > 0) THEN
            allocate (vec_recv_uniq (pushdata%num_req_uniq))
            vec_recv_uniq(:) = fillvalue
         ENDIF

         CALL worker_push_data_uniq (pushdata, &
            vec_send_real8 = vec_send, vec_recv_real8 = vec_recv_uniq, rfillvalue = fillvalue)

         IF (pushdata%num_req_uniq > 0) THEN

            vec_recv(:) = fillvalue

            DO j = 1, size(pushdata%addr_multi,2)

               sumarea = 0.

               DO i = 1, size(pushdata%addr_multi,1)
                  val = vec_recv_uniq(pushdata%addr_multi(i,j))
                  IF (val /= fillvalue) THEN
                     IF (vec_recv(j) == fillvalue) THEN
                        vec_recv(j) = val * pushdata%area_multi(i,j)
                     ELSE
                        vec_recv(j) = vec_recv(j) + val * pushdata%area_multi(i,j)
                     ENDIF
                     sumarea = sumarea + pushdata%area_multi(i,j)
                  ENDIF
               ENDDO

               IF (trim(mode) == 'average') THEN
                  IF (vec_recv(j) /= fillvalue) THEN
                     vec_recv(j) = vec_recv(j) / sumarea
                  ENDIF
               ENDIF
            ENDDO

            deallocate (vec_recv_uniq)
         ENDIF

      ENDIF

   END SUBROUTINE worker_push_data_multi_real8

   ! ----------
   SUBROUTINE worker_push_data_single_int32 (pushdata, vec_send, vec_recv, fillvalue)

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

         CALL worker_push_data_uniq (pushdata, &
            vec_send_int32 = vec_send, vec_recv_int32 = vec_recv_uniq, ifillvalue = fillvalue)

         IF (pushdata%num_req_uniq > 0) THEN
            vec_recv = vec_recv_uniq(pushdata%addr_single)
            deallocate (vec_recv_uniq)
         ENDIF

      ENDIF

   END SUBROUTINE worker_push_data_single_int32

   ! ---------
   SUBROUTINE worker_remap_data_pset2grid_real8 (remapdata, vec_in, vec_out, fillvalue, filter)

   IMPLICIT NONE

   type(worker_remapdata_type), intent(in) :: remapdata

   real(r8), intent(in)    :: vec_in (:)
   real(r8), intent(inout) :: vec_out(:)
   real(r8), intent(in)    :: fillvalue
   logical,  intent(in)    :: filter (:)

   ! Local Variables
   integer  :: iset, ipart, iloc
   real(r8) :: area


      IF (p_is_worker) THEN
         IF (remapdata%num_grid > 0) THEN

            vec_out(:) = fillvalue

            DO iset = 1, remapdata%npset
               IF (filter(iset) .and. (vec_in(iset) /= fillvalue)) THEN
                  DO ipart = 1, remapdata%npart(iset)

                     iloc = remapdata%part_to(iset)%val(ipart)
                     area = remapdata%areapart(iset)%val(ipart)

                     IF (vec_out(iloc) == fillvalue) THEN
                        vec_out(iloc) = vec_in(iset) * area
                     ELSE
                        vec_out(iloc) = vec_out(iloc) + vec_in(iset) * area
                     ENDIF

                  ENDDO
               ENDIF
            ENDDO

         ENDIF
      ENDIF

   END SUBROUTINE worker_remap_data_pset2grid_real8

   ! ---------
   SUBROUTINE worker_remap_data_grid2pset_real8 (remapdata, vec_in, vec_out, fillvalue, mode)

   IMPLICIT NONE

   type(worker_remapdata_type), intent(in) :: remapdata

   real(r8), intent(in)    :: vec_in (:)
   real(r8), intent(inout) :: vec_out(:)
   real(r8), intent(in)    :: fillvalue

   character(len=*), intent(in) :: mode

   ! Local Variables
   integer  :: iset, ipart, iloc
   real(r8) :: area, sumarea

      IF (p_is_worker) THEN
         IF (remapdata%npset > 0) THEN

            vec_out(:) = fillvalue

            DO iset = 1, remapdata%npset

               sumarea = 0.

               DO ipart = 1, remapdata%npart(iset)
                  iloc = remapdata%part_to(iset)%val(ipart)
                  area = remapdata%areapart(iset)%val(ipart)

                  IF (vec_in(iloc) /= fillvalue) THEN
                     IF (vec_out(iset) == fillvalue) THEN
                        vec_out(iset) = vec_in(iloc) * area
                     ELSE
                        vec_out(iset) = vec_out(iset) + vec_in(iloc) * area
                     ENDIF
                     sumarea = sumarea + area
                  ENDIF
               ENDDO

               IF (trim(mode) == 'average') THEN
                  IF (vec_out(iset) /= fillvalue) THEN
                     vec_out(iset) = vec_out(iset) / sumarea
                  ENDIF
               ENDIF
            ENDDO

         ENDIF
      ENDIF

   END SUBROUTINE worker_remap_data_grid2pset_real8

   ! ---------
   SUBROUTINE worker_pushdata_free_mem (this)

   IMPLICIT NONE
   type(worker_pushdata_type) :: this

      IF (allocated(this%addr_single )) deallocate(this%addr_single )
      IF (allocated(this%addr_multi  )) deallocate(this%addr_multi  )
      IF (allocated(this%area_multi  )) deallocate(this%area_multi  )
      IF (allocated(this%sum_area    )) deallocate(this%sum_area    )
      IF (allocated(this%self_from   )) deallocate(this%self_from   )
      IF (allocated(this%self_to     )) deallocate(this%self_to     )
      IF (allocated(this%n_to_other  )) deallocate(this%n_to_other  )
      IF (allocated(this%n_from_other)) deallocate(this%n_from_other)
      IF (allocated(this%to_other    )) deallocate(this%to_other    )
      IF (allocated(this%other_to    )) deallocate(this%other_to    )

   END SUBROUTINE worker_pushdata_free_mem

   ! ---------
   SUBROUTINE worker_remapdata_free_mem (this)

   IMPLICIT NONE
   type(worker_remapdata_type) :: this

      IF (allocated(this%ilon_me )) deallocate(this%ilon_me )
      IF (allocated(this%ilat_me )) deallocate(this%ilat_me )
      IF (allocated(this%ids_me  )) deallocate(this%ids_me  )
      IF (allocated(this%npart   )) deallocate(this%npart   )
      IF (allocated(this%sum_area)) deallocate(this%sum_area)
      IF (allocated(this%part_to )) deallocate(this%part_to )
      IF (allocated(this%areapart)) deallocate(this%areapart)

   END SUBROUTINE worker_remapdata_free_mem

END MODULE MOD_WorkerPushData
