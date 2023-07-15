#include <define.h>

#ifdef LATERAL_FLOW
MODULE MOD_Hydro_RiverNetwork
   !--------------------------------------------------------------------------------
   ! DESCRIPTION:
   ! 
   !    River networks: data and communication subroutines.
   !
   ! Created by Shupeng Zhang, May 2023
   !--------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Vars_Global, only : spval
   IMPLICIT NONE
   
   ! -- river parameters --
   REAL(r8), allocatable :: riverlen  (:)
   REAL(r8), allocatable :: riverelv  (:)
   REAL(r8), allocatable :: riverarea (:)
   REAL(r8), allocatable :: riverwdth (:)
   REAL(r8), allocatable :: riverdpth (:)
   
   INTEGER, allocatable :: riverdown  (:)  

   ! address of downstream river 
   ! > 0 on this process; 0 on other processes; -1 not found
   INTEGER, allocatable :: addrdown (:)

   REAL(r8), allocatable :: riverlen_ds  (:)
   REAL(r8), allocatable :: riverelv_ds  (:)
   REAL(r8), allocatable :: riverfac_ds  (:)

   TYPE :: river_sendrecv_type
      INTEGER :: nproc
      INTEGER, allocatable :: iproc (:)
      INTEGER, allocatable :: wdsp  (:)
      INTEGER, allocatable :: ndata (:)
      INTEGER, allocatable :: ups   (:)
      INTEGER, allocatable :: down  (:)
      INTEGER, allocatable :: iloc  (:)
   CONTAINS
      final :: river_sendrecv_free_mem
   END TYPE river_sendrecv_type

   TYPE(river_sendrecv_type), target :: river_up
   TYPE(river_sendrecv_type), target :: river_dn

   INTEGER, parameter :: SEND_DATA_DOWN_TO_UP = 1
   INTEGER, parameter :: SEND_DATA_UP_TO_DOWN = 2

CONTAINS
   
   ! ----------
   SUBROUTINE river_network_init ()

      USE MOD_SPMD_Task
      USE MOD_Namelist
      USE MOD_NetCDFSerial
      USE MOD_Mesh
      USE MOD_Hydro_SurfaceNetwork
      USE MOD_DataType
      USE MOD_Utils
      IMPLICIT NONE

      ! Local Variables
      CHARACTER(len=256) :: river_file

      INTEGER :: numbasin, ibasin, nbasin
      INTEGER :: iworker, mesg(4), isrc, idest, iproc
      INTEGER :: irecv, ifrom, ito, iup, idn, idata
      INTEGER :: nrecv, ndata, nup, ndn
      INTEGER :: iloc, iloc1, iloc2
   
      INTEGER , allocatable :: bindex (:)
      INTEGER , allocatable :: icache (:)
      REAL(r8), allocatable :: rcache (:)
      
      INTEGER , allocatable :: addrbasin (:,:)
      INTEGER , allocatable :: ndata_w (:)

      TYPE(pointer_int32_2d), allocatable :: exchange_w (:)
      INTEGER, allocatable :: exchange(:,:)
      INTEGER, allocatable :: basin_sorted(:), order(:)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      numbasin = numelm

      river_file = DEF_path_Catchment_data 

      IF (p_is_master) THEN
         CALL ncio_read_serial (river_file, 'river_downstream', riverdown)
         CALL ncio_read_serial (river_file, 'river_length'   ,  riverlen )
         CALL ncio_read_serial (river_file, 'river_elevation',  riverelv )
         CALL ncio_read_serial (river_file, 'river_depth    ',  riverdpth)

         riverlen = riverlen * 1.e3 ! km to m
      ENDIF

#ifdef USEMPI
      IF (p_is_master) THEN

         nbasin = size(riverdown)
         allocate (addrbasin (2,nbasin))

         addrbasin(:,:) = -1

         DO iworker = 1, p_np_worker

            CALL mpi_recv (mesg(1:2), 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
            isrc  = mesg(1)
            nrecv = mesg(2)

            IF (nrecv > 0) THEN
               
               allocate (bindex (nrecv))
               allocate (icache (nrecv))
               allocate (rcache (nrecv))

               CALL mpi_recv (bindex, nrecv, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               DO irecv = 1, nrecv
                  addrbasin(1,bindex(irecv)) = isrc
               ENDDO
               
               idest = isrc

               DO irecv = 1, nrecv
                  icache(irecv) = riverdown(bindex(irecv))
               ENDDO
               CALL mpi_send (icache, nrecv, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err) 

               DO irecv = 1, nrecv
                  rcache(irecv) = riverlen(bindex(irecv))
               ENDDO
               CALL mpi_send (rcache, nrecv, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err) 

               DO irecv = 1, nrecv
                  rcache(irecv) = riverelv(bindex(irecv))
               ENDDO
               CALL mpi_send (rcache, nrecv, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err) 

               DO irecv = 1, nrecv
                  rcache(irecv) = riverdpth(bindex(irecv))
               ENDDO
               CALL mpi_send (rcache, nrecv, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err) 

               deallocate (bindex)
               deallocate (icache)
               deallocate (rcache)

            ENDIF
         ENDDO

      ENDIF
#endif

      IF (p_is_worker) THEN
               
         IF (numbasin > 0) THEN
            allocate (bindex (numbasin))
            DO ibasin = 1, numbasin
               bindex(ibasin) = mesh(ibasin)%indx
            ENDDO
         ENDIF 

#ifdef USEMPI
         mesg(1:2) = (/p_iam_glb, numbasin/)
         CALL mpi_send (mesg(1:2), 2, MPI_INTEGER, p_root, mpi_tag_mesg, p_comm_glb, p_err) 

         IF (numbasin > 0) THEN
            CALL mpi_send (bindex, numbasin, MPI_INTEGER, &
               p_root, mpi_tag_data, p_comm_glb, p_err) 

            allocate (riverdown (numbasin))
            CALL mpi_recv (riverdown, numbasin, MPI_INTEGER, &
               p_root, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (riverlen (numbasin))
            CALL mpi_recv (riverlen, numbasin, MPI_REAL8, &
               p_root, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (riverelv (numbasin))
            CALL mpi_recv (riverelv, numbasin, MPI_REAL8, &
               p_root, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (riverdpth (numbasin))
            CALL mpi_recv (riverdpth, numbasin, MPI_REAL8, &
               p_root, mpi_tag_data, p_comm_glb, p_stat, p_err)
         ENDIF
#else
         IF (numbasin > 0) THEN

            riverdown = riverdown(bindex)
            riverlen  = riverlen (bindex)
            riverelv  = riverelv (bindex)
            riverdpth = riverdpth(bindex)

         ENDIF
#endif

      ENDIF 

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef USEMPI
      IF (p_is_master) THEN

         allocate (ndata_w (0:p_np_worker-1))

         ndata_w(:) = 0
         DO ibasin = 1, nbasin
            IF (riverdown(ibasin) >= 1) THEN
               addrbasin(2,ibasin) = addrbasin(1,riverdown(ibasin))
               IF ((addrbasin(1,ibasin) /= -1) .and. (addrbasin(2,ibasin) /= -1) &
                  .and. (addrbasin(1,ibasin) /= addrbasin(2,ibasin))) THEN
                  ifrom = p_itis_worker(addrbasin(1,ibasin))
                  ito   = p_itis_worker(addrbasin(2,ibasin))
                  ndata_w(ifrom) = ndata_w(ifrom) + 1
                  ndata_w(ito)   = ndata_w(ito)   + 1
               ENDIF
            ENDIF
         ENDDO

         allocate (exchange_w (0:p_np_worker-1))
         DO iworker = 0, p_np_worker-1
            IF (ndata_w(iworker) > 0) THEN
               allocate (exchange_w(iworker)%val (4,ndata_w(iworker)))
            ENDIF
         ENDDO

         ndata_w(:) = 0
         DO ibasin = 1, nbasin
            IF ((addrbasin(1,ibasin) /= -1) .and. (addrbasin(2,ibasin) /= -1) &
               .and. (addrbasin(1,ibasin) /= addrbasin(2,ibasin))) THEN
               ifrom = p_itis_worker(addrbasin(1,ibasin))
               ito   = p_itis_worker(addrbasin(2,ibasin))
               ndata_w(ifrom) = ndata_w(ifrom) + 1
               ndata_w(ito)   = ndata_w(ito)   + 1

               exchange_w(ifrom)%val(:,ndata_w(ifrom)) = &
                  (/addrbasin(1,ibasin), ibasin, addrbasin(2,ibasin), riverdown(ibasin)/)
               exchange_w(ito)%val(:,ndata_w(ito)) = &
                  (/addrbasin(1,ibasin), ibasin, addrbasin(2,ibasin), riverdown(ibasin)/)
            ENDIF
         ENDDO

         DO iworker = 0, p_np_worker-1
            CALL mpi_send (ndata_w(iworker), 1, MPI_INTEGER, &
               p_address_worker(iworker), mpi_tag_size, p_comm_glb, p_err) 
            IF (ndata_w(iworker) > 0) THEN
               CALL mpi_send (exchange_w(iworker)%val, 4*ndata_w(iworker), MPI_INTEGER, &
                  p_address_worker(iworker), mpi_tag_data, p_comm_glb, p_err) 
            ENDIF
         ENDDO
      ENDIF
#endif

      IF (p_is_worker) THEN
#ifdef USEMPI
         CALL mpi_recv (ndata, 1, MPI_INTEGER, p_root, mpi_tag_size, p_comm_glb, p_stat, p_err)
         IF (ndata > 0) THEN
            allocate (exchange(4,ndata))
            CALL mpi_recv (exchange, 4*ndata, MPI_INTEGER, &
               p_root, mpi_tag_data, p_comm_glb, p_stat, p_err)
         ENDIF
#endif

         IF (numbasin > 0) THEN
            
            allocate (basin_sorted (numbasin))
            allocate (order (numbasin))
            basin_sorted = bindex
            order = (/(ibasin, ibasin = 1, numbasin)/)

            CALL quicksort (numbasin, basin_sorted, order)

            allocate (addrdown (numbasin))
            addrdown(:) = -1

            DO ibasin = 1, numbasin
               IF (riverdown(ibasin) > 0) THEN
                  iloc = find_in_sorted_list1 (riverdown(ibasin), numbasin, basin_sorted)
                  IF (iloc > 0) THEN
                     addrdown(ibasin) = order(iloc)
                  ENDIF
               ENDIF
            ENDDO

#ifdef USEMPI
            IF (ndata > 0) THEN
               nup = count(exchange(3,:) == p_iam_glb)
               ndn = count(exchange(1,:) == p_iam_glb)

               IF (nup > 0) allocate (river_up%iproc (nup))
               IF (nup > 0) allocate (river_up%ups   (nup))
               IF (nup > 0) allocate (river_up%down  (nup))
               IF (nup > 0) allocate (river_up%iloc  (nup))
               
               IF (ndn > 0) allocate (river_dn%iproc (ndn))
               IF (ndn > 0) allocate (river_dn%ups   (ndn))
               IF (ndn > 0) allocate (river_dn%down  (ndn))
               IF (ndn > 0) allocate (river_dn%iloc  (ndn))
               
               iup = 0
               idn = 0
               DO idata = 1, ndata
                  IF (exchange(3,idata) == p_iam_glb) THEN
                     CALL insert_into_sorted_list2 (exchange(2,idata), exchange(1,idata), &
                        iup, river_up%ups, river_up%iproc, iloc)
                  ELSEIF (exchange(1,idata) == p_iam_glb) THEN
                     CALL insert_into_sorted_list2 (exchange(2,idata), exchange(3,idata), &
                        idn, river_dn%ups, river_dn%iproc, iloc)
                  ENDIF
               ENDDO

               DO idata = 1, ndata
                  IF (exchange(3,idata) == p_iam_glb) THEN
                     
                     iloc1 = find_in_sorted_list2 (exchange(2,idata), exchange(1,idata), &
                        nup, river_up%ups, river_up%iproc)

                     river_up%down(iloc1) = exchange(4,idata)

                     iloc2 = find_in_sorted_list1 (exchange(4,idata), numbasin, basin_sorted)
                     river_up%iloc(iloc1) = order(iloc2)

                  ELSEIF (exchange(1,idata) == p_iam_glb) THEN
                     
                     iloc1 = find_in_sorted_list2 (exchange(2,idata), exchange(3,idata), &
                        ndn, river_dn%ups, river_dn%iproc)
                     
                     river_dn%down(iloc1) = exchange(4,idata)
                     
                     iloc2 = find_in_sorted_list1 (exchange(2,idata), numbasin, basin_sorted)
                     river_dn%iloc(iloc1) = order(iloc2)
                  ENDIF
               ENDDO

               IF (nup > 0) THEN

                  river_up%nproc = 1
                  DO iup = 2, nup
                     IF (river_up%iproc(iup) /= river_up%iproc(iup-1)) THEN
                        river_up%nproc = river_up%nproc + 1
                     ENDIF
                  ENDDO

                  allocate (river_up%wdsp (river_up%nproc))
                  allocate (river_up%ndata(river_up%nproc))

                  river_up%ndata(:) = 0

                  iproc = 1
                  river_up%wdsp (1) = 0
                  river_up%ndata(1) = 1
                  DO iup = 2, nup
                     IF (river_up%iproc(iup) /= river_up%iproc(iup-1)) THEN
                        iproc = iproc + 1
                        river_up%wdsp (iproc) = iup - 1
                        river_up%ndata(iproc) = 1
                     ELSE
                        river_up%ndata(iproc) = river_up%ndata(iproc) + 1
                     ENDIF
                  ENDDO

               ELSE
                  river_up%nproc = 0
               ENDIF
               
               IF (ndn > 0) THEN

                  river_dn%nproc = 1
                  DO idn = 2, ndn
                     IF (river_dn%iproc(idn) /= river_dn%iproc(idn-1)) THEN
                        river_dn%nproc = river_dn%nproc + 1
                     ENDIF
                  ENDDO

                  allocate (river_dn%wdsp (river_dn%nproc))
                  allocate (river_dn%ndata(river_dn%nproc))

                  river_dn%ndata(:) = 0

                  iproc = 1
                  river_dn%wdsp (1) = 0
                  river_dn%ndata(1) = 1
                  DO idn = 2, ndn
                     IF (river_dn%iproc(idn) /= river_dn%iproc(idn-1)) THEN
                        iproc = iproc + 1
                        river_dn%wdsp (iproc) = idn - 1
                        river_dn%ndata(iproc) = 1
                     ELSE
                        river_dn%ndata(iproc) = river_dn%ndata(iproc) + 1
                     ENDIF
                  ENDDO
               ELSE
                  river_dn%nproc = 0
               ENDIF
            ENDIF

            addrdown(river_dn%iloc) = 0

#endif
         ENDIF
      ENDIF

      IF (p_is_worker) THEN

         IF (numbasin > 0) THEN

            allocate (riverarea (numbasin))
            allocate (riverwdth (numbasin))

            DO ibasin = 1, numbasin
               riverarea(ibasin) = surface_network(ibasin)%area(1)
               riverwdth(ibasin) = riverarea(ibasin) / riverlen(ibasin)

               ! modify height above nearest drainage data to consider river depth
               surface_network(ibasin)%hand(1) = &
                  surface_network(ibasin)%hand(1) + riverdpth(ibasin)
            ENDDO

         ENDIF

      ENDIF

      IF (allocated(bindex      )) deallocate(bindex      )
      IF (allocated(addrbasin   )) deallocate(addrbasin   )
      IF (allocated(ndata_w     )) deallocate(ndata_w     )
      IF (allocated(exchange_w  )) deallocate(exchange_w  )
      IF (allocated(exchange    )) deallocate(exchange    )
      IF (allocated(basin_sorted)) deallocate(basin_sorted)
      IF (allocated(order       )) deallocate(order       )

      IF (p_is_worker) THEN
         IF (numbasin > 0) THEN
            allocate (riverlen_ds  (numbasin))
            allocate (riverelv_ds  (numbasin))
            allocate (riverfac_ds  (numbasin))

            DO ibasin = 1, numbasin
               IF (addrdown(ibasin) > 0) THEN
                  riverlen_ds (ibasin) = riverlen (addrdown(ibasin)) 
                  riverelv_ds (ibasin) = riverelv (addrdown(ibasin)) 
                  riverfac_ds (ibasin) = (riverwdth(ibasin) + riverwdth(addrdown(ibasin))) * 0.5
               ELSE
                  riverlen_ds (ibasin) = spval
                  riverelv_ds (ibasin) = spval
                  riverfac_ds (ibasin) = spval
               ENDIF
            ENDDO
         ENDIF

#ifdef USEMPI
         CALL river_data_exchange (SEND_DATA_DOWN_TO_UP, accum = .false., &
            vec_send1 = riverlen,  vec_recv1 = riverlen_ds, &
            vec_send2 = riverelv,  vec_recv2 = riverelv_ds, &
            vec_send3 = riverwdth, vec_recv3 = riverfac_ds)

         DO ibasin = 1, numbasin
            IF (addrdown(ibasin) == 0) THEN
               riverfac_ds(ibasin) = (riverwdth(ibasin) + riverfac_ds(ibasin)) * 0.5
            ENDIF
         ENDDO
#endif
         
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE river_network_init

   ! ----------
#ifdef USEMPI
   SUBROUTINE river_data_exchange (direction, accum, &
         vec_send1, vec_recv1, vec_send2, vec_recv2, &
         vec_send3, vec_recv3, vec_send4, vec_recv4 )

      USE MOD_Precision
      USE MOD_SPMD_Task
      IMPLICIT NONE

      INTEGER, intent(in) :: direction
      LOGICAL, intent(in) :: accum

      REAL(r8), intent(inout) :: vec_send1(:), vec_recv1(:)
      REAL(r8), intent(inout), optional :: vec_send2(:), vec_recv2(:)
      REAL(r8), intent(inout), optional :: vec_send3(:), vec_recv3(:)
      REAL(r8), intent(inout), optional :: vec_send4(:), vec_recv4(:)

      ! Local Variables

      TYPE(river_sendrecv_type), pointer :: send_pointer
      INTEGER :: nproc_send, ndatasend, idest
      INTEGER,  allocatable :: req_send(:,:)
      REAL(r8), allocatable :: sendcache1(:)
      REAL(r8), allocatable :: sendcache2(:)
      REAL(r8), allocatable :: sendcache3(:)
      REAL(r8), allocatable :: sendcache4(:)

      TYPE(river_sendrecv_type), pointer :: recv_pointer
      INTEGER :: nproc_recv, ndatarecv, isrc
      INTEGER,  allocatable :: req_recv(:,:)
      REAL(r8), allocatable :: recvcache1(:)
      REAL(r8), allocatable :: recvcache2(:)
      REAL(r8), allocatable :: recvcache3(:)
      REAL(r8), allocatable :: recvcache4(:)

      INTEGER :: nvec, iproc, i, istt, iend, ndata

      IF (p_is_worker) THEN

         CALL mpi_barrier (p_comm_worker, p_err)

         IF (direction == SEND_DATA_DOWN_TO_UP) THEN
            send_pointer => river_up
            recv_pointer => river_dn
         elseif (direction == SEND_DATA_UP_TO_DOWN) THEN
            send_pointer => river_dn
            recv_pointer => river_up
         ENDIF

         nproc_send = send_pointer%nproc
         IF (nproc_send > 0) THEN

            ndatasend = sum(send_pointer%ndata)
            
            nvec = 1
            allocate (sendcache1(ndatasend))
            DO i = 1, ndatasend
               sendcache1(i) = vec_send1(send_pointer%iloc(i))
            ENDDO

            IF (present(vec_send2) .and. present(vec_recv2)) THEN
               nvec = nvec + 1
               allocate (sendcache2(ndatasend))
               DO i = 1, ndatasend
                  sendcache2(i) = vec_send2(send_pointer%iloc(i))
               ENDDO

               IF (present(vec_send3) .and. present(vec_recv3)) THEN
                  nvec = nvec + 1
                  allocate (sendcache3(ndatasend))
                  DO i = 1, ndatasend
                     sendcache3(i) = vec_send3(send_pointer%iloc(i))
                  ENDDO

                  IF (present(vec_send4) .and. present(vec_recv4)) THEN
                     nvec = nvec + 1
                     allocate (sendcache4(ndatasend))
                     DO i = 1, ndatasend
                        sendcache4(i) = vec_send4(send_pointer%iloc(i))
                     ENDDO
                  ENDIF
               ENDIF
            ENDIF

            allocate (req_send(nvec,nproc_send))
           
            DO iproc = 1, nproc_send
               ndata = send_pointer%ndata(iproc)
               istt  = send_pointer%wdsp (iproc) + 1
               iend  = send_pointer%wdsp (iproc) + ndata
               idest = send_pointer%iproc(istt)

               CALL mpi_isend(sendcache1(istt:iend), ndata, MPI_REAL8, &
                  idest, 101, p_comm_glb, req_send(1,iproc), p_err)
               IF (present(vec_send2) .and. present(vec_recv2)) THEN
                  CALL mpi_isend(sendcache2(istt:iend), ndata, MPI_REAL8, &
                     idest, 102, p_comm_glb, req_send(2,iproc), p_err)
                  IF (present(vec_send3) .and. present(vec_recv3)) THEN
                     CALL mpi_isend(sendcache3(istt:iend), ndata, MPI_REAL8, &
                        idest, 103, p_comm_glb, req_send(3,iproc), p_err)
                     IF (present(vec_send4) .and. present(vec_recv4)) THEN
                        CALL mpi_isend(sendcache4(istt:iend), ndata, MPI_REAL8, &
                           idest, 104, p_comm_glb, req_send(4,iproc), p_err)
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO

         ENDIF

         nproc_recv = recv_pointer%nproc
         IF (nproc_recv > 0) THEN

            ndatarecv = sum(recv_pointer%ndata)
            
            nvec = 1
            allocate (recvcache1(ndatarecv))
            IF (present(vec_send2) .and. present(vec_recv2)) THEN
               nvec = nvec + 1
               allocate (recvcache2(ndatarecv))
               IF (present(vec_send3) .and. present(vec_recv3)) THEN
                  nvec = nvec + 1
                  allocate (recvcache3(ndatarecv))
                  IF (present(vec_send4) .and. present(vec_recv4)) THEN
                     nvec = nvec + 1
                     allocate (recvcache4(ndatarecv))
                  ENDIF
               ENDIF
            ENDIF
            
            allocate (req_recv(nvec,nproc_recv))
            
            DO iproc = 1, nproc_recv
               ndata = recv_pointer%ndata(iproc)
               istt  = recv_pointer%wdsp(iproc) + 1
               iend  = recv_pointer%wdsp(iproc) + ndata
               isrc  = recv_pointer%iproc(istt)

               CALL mpi_irecv(recvcache1(istt:iend), ndata, MPI_REAL8, &
                  isrc, 101, p_comm_glb, req_recv(1,iproc), p_err)
               IF (present(vec_send2) .and. present(vec_recv2)) THEN
                  CALL mpi_irecv(recvcache2(istt:iend), ndata, MPI_REAL8, &
                     isrc, 102, p_comm_glb, req_recv(2,iproc), p_err)
                  IF (present(vec_send3) .and. present(vec_recv3)) THEN
                     CALL mpi_irecv(recvcache3(istt:iend), ndata, MPI_REAL8, &
                        isrc, 103, p_comm_glb, req_recv(3,iproc), p_err)
                     IF (present(vec_send4) .and. present(vec_recv4)) THEN
                        CALL mpi_irecv(recvcache4(istt:iend), ndata, MPI_REAL8, &
                           isrc, 104, p_comm_glb, req_recv(4,iproc), p_err)
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO

         ENDIF

         IF (nproc_recv > 0) THEN

            CALL mpi_waitall(nvec*nproc_recv, req_recv, MPI_STATUSES_IGNORE, p_err)
            ! write(*,*) 'p error', p_err

            IF (accum) THEN
               DO i = 1, ndatarecv
                  vec_recv1(recv_pointer%iloc(i)) = &
                     vec_recv1(recv_pointer%iloc(i)) + recvcache1(i)
               ENDDO

               IF (present(vec_send2) .and. present(vec_recv2)) THEN
                  DO i = 1, ndatarecv
                     vec_recv2(recv_pointer%iloc(i)) = &
                        vec_recv2(recv_pointer%iloc(i)) + recvcache2(i)
                  ENDDO

                  IF (present(vec_send3) .and. present(vec_recv3)) THEN
                     DO i = 1, ndatarecv
                        vec_recv3(recv_pointer%iloc(i)) = &
                           vec_recv3(recv_pointer%iloc(i)) + recvcache3(i)
                     ENDDO

                     IF (present(vec_send4) .and. present(vec_recv4)) THEN
                        DO i = 1, ndatarecv
                           vec_recv4(recv_pointer%iloc(i)) = &
                              vec_recv4(recv_pointer%iloc(i)) + recvcache4(i)
                        ENDDO
                     ENDIF
                  ENDIF
               ENDIF
            ELSE
               DO i = 1, ndatarecv
                  vec_recv1(recv_pointer%iloc(i)) = recvcache1(i)
               ENDDO

               IF (present(vec_send2) .and. present(vec_recv2)) THEN
                  DO i = 1, ndatarecv
                     vec_recv2(recv_pointer%iloc(i)) = recvcache2(i)
                  ENDDO

                  IF (present(vec_send3) .and. present(vec_recv3)) THEN
                     DO i = 1, ndatarecv
                        vec_recv3(recv_pointer%iloc(i)) = recvcache3(i)
                     ENDDO

                     IF (present(vec_send4) .and. present(vec_recv4)) THEN
                        DO i = 1, ndatarecv
                           vec_recv4(recv_pointer%iloc(i)) = recvcache4(i)
                        ENDDO
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         IF (nproc_send > 0) THEN
            CALL mpi_waitall(nvec*nproc_send, req_send, MPI_STATUSES_IGNORE, p_err)
         ENDIF

         IF (allocated(req_send  )) deallocate(req_send)
         IF (allocated(sendcache1)) deallocate(sendcache1)
         IF (allocated(sendcache2)) deallocate(sendcache2)
         IF (allocated(sendcache3)) deallocate(sendcache3)
         IF (allocated(sendcache4)) deallocate(sendcache4)
         
         IF (allocated(req_recv  )) deallocate(req_recv)
         IF (allocated(recvcache1)) deallocate(recvcache1)
         IF (allocated(recvcache2)) deallocate(recvcache2)
         IF (allocated(recvcache3)) deallocate(recvcache3)
         IF (allocated(recvcache4)) deallocate(recvcache4)

         CALL mpi_barrier (p_comm_worker, p_err)

      ENDIF

   END SUBROUTINE river_data_exchange
#endif

   ! ----------
   SUBROUTINE river_network_final ()

      IMPLICIT NONE

      IF (allocated(riverlen )) deallocate(riverlen )
      IF (allocated(riverelv )) deallocate(riverelv )
      IF (allocated(riverarea)) deallocate(riverarea)
      IF (allocated(riverwdth)) deallocate(riverwdth)
      IF (allocated(riverdpth)) deallocate(riverdpth)
      
      IF (allocated(riverdown )) deallocate(riverdown )
      IF (allocated(addrdown  )) deallocate(addrdown  )

      IF (allocated(riverlen_ds))  deallocate(riverlen_ds)
      IF (allocated(riverelv_ds))  deallocate(riverelv_ds)
      IF (allocated(riverfac_ds))  deallocate(riverfac_ds)

   END SUBROUTINE river_network_final

   ! ---------
   SUBROUTINE river_sendrecv_free_mem (this)
      
      IMPLICIT NONE
      TYPE(river_sendrecv_type) :: this

      IF (allocated(this%iproc)) deallocate(this%iproc)
      IF (allocated(this%wdsp )) deallocate(this%wdsp )
      IF (allocated(this%ndata)) deallocate(this%ndata)
      IF (allocated(this%ups  )) deallocate(this%ups  )
      IF (allocated(this%down )) deallocate(this%down )
      IF (allocated(this%iloc )) deallocate(this%iloc )

   END SUBROUTINE river_sendrecv_free_mem

END MODULE MOD_Hydro_RiverNetwork
#endif
