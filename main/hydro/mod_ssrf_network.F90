#include <define.h>

#ifdef LATERAL_FLOW
MODULE mod_ssrf_network

   USE precision
   USE mod_data_type
   IMPLICIT NONE
   
   ! -- ssrf parameters --
   INTEGER, allocatable :: num_nb (:)

   TYPE(pointer_int32_1d), allocatable :: idxbsn_nb (:)
   TYPE(pointer_int32_2d), allocatable :: addr_nb   (:)

   TYPE(pointer_real8_1d), allocatable :: dist_nb   (:)  ! m
   TYPE(pointer_real8_1d), allocatable :: lenbdr_nb (:)  ! m
   TYPE(pointer_real8_1d), allocatable :: area_nb   (:)  ! m^2
   TYPE(pointer_real8_1d), allocatable :: elva_nb   (:)  ! m
   TYPE(pointer_real8_1d), allocatable :: slope_nb  (:)  ! unitless
      
   REAL(r8), allocatable :: area_b(:)
   REAL(r8), allocatable :: elva_b(:)

   ! -- ssrf variables --
   TYPE(pointer_real8_1d), allocatable :: theta_a_nb (:)
   TYPE(pointer_real8_1d), allocatable :: zwt_nb     (:)
   TYPE(pointer_real8_1d), allocatable :: Ks_nb      (:)

   TYPE ssrf_sendrecv_type
      INTEGER :: ndata
      INTEGER, allocatable :: bindx (:)
      INTEGER, allocatable :: ibsn  (:)
   END TYPE ssrf_sendrecv_type

   TYPE(ssrf_sendrecv_type), allocatable :: recvaddr(:)
   TYPE(ssrf_sendrecv_type), allocatable :: sendaddr(:)

CONTAINS
   
   ! ----------
   SUBROUTINE ssrf_init ()

      USE spmd_task
      USE mod_namelist
      USE ncio_serial
      USE mod_mesh
      USE MOD_LandElm
      USE mod_colm_debug
      USE mod_drainage_network
      USE mod_utils
      IMPLICIT NONE

      ! Local Variables
      CHARACTER(len=256) :: ssrf_file

      INTEGER :: numbasin, ibasin
      INTEGER :: iwork, mesg(2), isrc, idest
      INTEGER :: nrecv, irecv
      INTEGER :: iloc, iloc1, iloc2
      INTEGER :: nnb, nnbinq, inb, ndata
   
      INTEGER :: maxnnb
      INTEGER , allocatable :: nnball   (:)
      INTEGER , allocatable :: idxnball (:,:)
      REAL(r8), allocatable :: lenbdall (:,:)

      INTEGER , allocatable :: addrbasin(:)

      INTEGER , allocatable :: bindex  (:)
      INTEGER , allocatable :: icache1 (:)
      INTEGER , allocatable :: icache2 (:,:)
      REAL(r8), allocatable :: rcache2 (:,:)

      INTEGER, allocatable :: basin_sorted(:), order(:)
      INTEGER, allocatable :: idxinq(:), addrinq(:)

      LOGICAL, allocatable :: mask(:)

      REAL(r8), allocatable :: rlon_b(:), rlat_b(:)
      TYPE(pointer_real8_1d), allocatable :: rlon_nb(:), rlat_nb(:)


#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      numbasin = numelm

      ssrf_file = DEF_path_Catchment_data 

      IF (p_is_master) THEN
         CALL ncio_read_serial (ssrf_file, 'basin_num_neighbour', nnball  )
         CALL ncio_read_serial (ssrf_file, 'basin_idx_neighbour', idxnball)
         CALL ncio_read_serial (ssrf_file, 'basin_len_border'   , lenbdall)

         maxnnb = size(idxnball,1)

         lenbdall = lenbdall * 1.e3 ! km to m
      ENDIF

#ifdef USEMPI

      CALL mpi_bcast (maxnnb, 1, MPI_INTEGER, p_root, p_comm_glb, p_err)

      IF (p_is_master) THEN

         allocate (addrbasin (size(nnball)))
         addrbasin(:) = -1

         DO iwork = 0, p_np_worker-1

            CALL mpi_recv (mesg(1:2), 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
            isrc  = mesg(1)
            nrecv = mesg(2)

            IF (nrecv > 0) THEN
               
               allocate (bindex  (nrecv))
               allocate (icache1 (nrecv))
               allocate (icache2 (maxnnb,nrecv))
               allocate (rcache2 (maxnnb,nrecv))

               CALL mpi_recv (bindex, nrecv, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               addrbasin(bindex) = isrc
               
               idest = isrc

               icache1 = nnball(bindex)
               CALL mpi_send (icache1, nrecv, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err) 

               DO irecv = 1, nrecv
                  icache2(:,irecv) = idxnball(:,bindex(irecv))
               ENDDO
               CALL mpi_send (icache2, maxnnb*nrecv, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err) 

               DO irecv = 1, nrecv
                  rcache2(:,irecv) = lenbdall(:,bindex(irecv))
               ENDDO
               CALL mpi_send (rcache2, maxnnb*nrecv, MPI_REAL8, &
                  idest, mpi_tag_data, p_comm_glb, p_err) 

               deallocate (bindex )
               deallocate (icache1)
               deallocate (icache2)
               deallocate (rcache2)

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

            allocate (nnball (numbasin))
            CALL mpi_recv (nnball, numbasin, MPI_INTEGER, &
               p_root, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (idxnball (maxnnb,numbasin))
            CALL mpi_recv (idxnball, maxnnb*numbasin, MPI_INTEGER, &
               p_root, mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (lenbdall (maxnnb,numbasin))
            CALL mpi_recv (lenbdall, maxnnb*numbasin, MPI_REAL8, &
               p_root, mpi_tag_data, p_comm_glb, p_stat, p_err)
         ENDIF
#else
         allocate (icache1 (numbasin))
         allocate (icache2 (maxnnb,numbasin))
         allocate (rcache2 (maxnnb,numbasin))

         icache1 = nnball  
         icache2 = idxnball
         rcache2 = lenbdall

         DO ibasin = 1, numbasin
            nnball   (ibasin)   = icache1   (bindex(ibasin))
            idxnball (:,ibasin) = icache2 (:,bindex(ibasin))
            lenbdall (:,ibasin) = rcache2 (:,bindex(ibasin))
         ENDDO

         deallocate (icache1, icache2, rcache2)
#endif

         IF (numbasin > 0) THEN
            allocate (num_nb    (numbasin))
            allocate (addr_nb   (numbasin))
            allocate (idxbsn_nb (numbasin))
            allocate (lenbdr_nb (numbasin))

            DO ibasin = 1, numbasin
               nnb = nnball(ibasin)
               num_nb(ibasin) = nnb

               IF (nnb > 0) THEN
                  allocate (idxbsn_nb(ibasin)%val (nnb))
                  idxbsn_nb(ibasin)%val = idxnball(1:nnb,ibasin)

                  allocate (lenbdr_nb(ibasin)%val (nnb))
                  lenbdr_nb(ibasin)%val = lenbdall(1:nnb,ibasin)

                  allocate (addr_nb(ibasin)%val (2,nnb))
                  addr_nb(ibasin)%val(1,:) = 0
               ENDIF
            ENDDO
         ENDIF

      ENDIF 

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef USEMPI
      IF (p_is_master) THEN
         DO iwork = 0, p_np_worker-1

            CALL mpi_recv (mesg(1:2), 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
            isrc  = mesg(1)
            nrecv = mesg(2)

            IF (nrecv > 0) THEN
               allocate (bindex  (nrecv))
               allocate (icache1 (nrecv))
               
               CALL mpi_recv (bindex, nrecv, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               icache1 = addrbasin(bindex)

               idest = isrc
               CALL mpi_send (icache1, nrecv, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err) 

               deallocate(bindex, icache1)
            ENDIF
         ENDDO
      ENDIF 
#endif

      IF (p_is_worker) THEN

         IF (numbasin > 0) THEN
            allocate (basin_sorted (numbasin))
            allocate (order (numbasin))

            basin_sorted = bindex
            order = (/(ibasin, ibasin = 1, numbasin)/)

            CALL quicksort (numbasin, basin_sorted, order)

#ifdef USEMPI
            allocate(idxinq (numbasin)) 
#endif

            nnbinq = 0
            DO ibasin = 1, numbasin
               DO inb = 1, num_nb(ibasin)
                  iloc = find_in_sorted_list1 (idxbsn_nb(ibasin)%val(inb), numbasin, basin_sorted)
                  IF (iloc > 0) THEN
                     addr_nb(ibasin)%val(1,inb) = -1
                     addr_nb(ibasin)%val(2,inb) = order(iloc)
#ifdef USEMPI
                  ELSE
                     CALL insert_into_sorted_list1 (idxbsn_nb(ibasin)%val(inb), nnbinq, idxinq, iloc)

                     IF (nnbinq == size(idxinq)) THEN
                        allocate (icache1(nnbinq))
                        icache1 = idxinq
                        deallocate (idxinq)
                        allocate (idxinq (ceiling(1.2*nnbinq)))
                        idxinq(1:nnbinq) = icache1
                        deallocate (icache1)
                     ENDIF
#endif
                  ENDIF
               ENDDO
            ENDDO
         ELSE
            nnbinq = 0
         ENDIF 

#ifdef USEMPI
         mesg(1:2) = (/p_iam_glb, nnbinq/)
         CALL mpi_send (mesg(1:2), 2, MPI_INTEGER, p_root, mpi_tag_mesg, p_comm_glb, p_err) 

         IF (nnbinq > 0) THEN
            
            CALL mpi_send (idxinq(1:nnbinq), nnbinq, MPI_INTEGER, &
               p_root, mpi_tag_data, p_comm_glb, p_err) 

            allocate (addrinq (nnbinq))
            CALL mpi_recv (addrinq, nnbinq, MPI_INTEGER, &
               p_root, mpi_tag_data, p_comm_glb, p_stat, p_err)

         ENDIF

         IF (nnbinq > 0) allocate(mask (nnbinq))

         allocate (recvaddr (0:p_np_worker-1))
         DO iwork = 0, p_np_worker-1
            IF (nnbinq > 0) THEN
               mask = (addrinq == p_address_worker(iwork))
               ndata = count(mask)
            ELSE
               ndata = 0
            ENDIF 

            recvaddr(iwork)%ndata = ndata
            IF (ndata > 0) THEN
               allocate (recvaddr(iwork)%bindx (ndata))
               recvaddr(iwork)%bindx = pack(idxinq, mask)
            ENDIF
         ENDDO

         IF (nnbinq > 0) deallocate(mask)
      
         DO ibasin = 1, numbasin
            DO inb = 1, num_nb(ibasin)
               IF (addr_nb(ibasin)%val(1,inb) == 0) THEN

                  iloc = find_in_sorted_list1 (idxbsn_nb(ibasin)%val(inb), nnbinq, idxinq(1:nnbinq))

                  iwork = p_itis_worker(addrinq(iloc))
                  iloc1 = find_in_sorted_list1 (idxbsn_nb(ibasin)%val(inb), &
                     recvaddr(iwork)%ndata, recvaddr(iwork)%bindx)
                  
                  addr_nb(ibasin)%val(1,inb) = iwork
                  addr_nb(ibasin)%val(2,inb) = iloc1
               ENDIF
            ENDDO
         ENDDO

         allocate (sendaddr (0:p_np_worker-1))
         DO iwork = 0, p_np_worker-1
            sendaddr(iwork)%ndata = 0
         ENDDO

         DO ibasin = 1, numbasin
            DO inb = 1, num_nb(ibasin)
               IF (addr_nb(ibasin)%val(1,inb) >= 0) THEN
                  iwork = addr_nb(ibasin)%val(1,inb)
                  sendaddr(iwork)%ndata = sendaddr(iwork)%ndata + 1
               ENDIF
            ENDDO
         ENDDO

         DO iwork = 0, p_np_worker-1
            IF (sendaddr(iwork)%ndata > 0) THEN 
               allocate (sendaddr(iwork)%bindx (sendaddr(iwork)%ndata))
               sendaddr(iwork)%ndata = 0
            ENDIF
         ENDDO

         DO ibasin = 1, numbasin
            DO inb = 1, num_nb(ibasin)
               IF (addr_nb(ibasin)%val(1,inb) >= 0) THEN
                  iwork = addr_nb(ibasin)%val(1,inb)
                  CALL insert_into_sorted_list1 (bindex(ibasin), &
                     sendaddr(iwork)%ndata, sendaddr(iwork)%bindx, iloc)
               ENDIF
            ENDDO
         ENDDO
         
         DO iwork = 0, p_np_worker-1
            IF (sendaddr(iwork)%ndata > 0) THEN 
               IF (sendaddr(iwork)%ndata < size(sendaddr(iwork)%bindx)) THEN
                  allocate (icache1 (sendaddr(iwork)%ndata))
                  icache1 = sendaddr(iwork)%bindx(1:sendaddr(iwork)%ndata)

                  deallocate (sendaddr(iwork)%bindx)
                  allocate (sendaddr(iwork)%bindx (sendaddr(iwork)%ndata))
                  sendaddr(iwork)%bindx = icache1

                  deallocate (icache1)
               ENDIF
            ENDIF
         ENDDO
         
         DO iwork = 0, p_np_worker-1
            IF (sendaddr(iwork)%ndata > 0) THEN
               allocate (sendaddr(iwork)%ibsn (sendaddr(iwork)%ndata))

               DO inb = 1, sendaddr(iwork)%ndata
                  iloc = find_in_sorted_list1 (sendaddr(iwork)%bindx(inb), numbasin, basin_sorted)
                  sendaddr(iwork)%ibsn(inb) = order(iloc)
               ENDDO
            ENDIF
         ENDDO
#endif
      ENDIF

      IF (allocated(nnball   )) deallocate(nnball   )
      IF (allocated(idxnball )) deallocate(idxnball )
      IF (allocated(lenbdall )) deallocate(lenbdall )
      IF (allocated(addrbasin)) deallocate(addrbasin)
      IF (allocated(bindex ))   deallocate(bindex )
      IF (allocated(icache1))   deallocate(icache1)
      IF (allocated(icache2))   deallocate(icache2)
      IF (allocated(rcache2))   deallocate(rcache2)

      IF (allocated(basin_sorted)) deallocate(basin_sorted)
      IF (allocated(order  )) deallocate(order  )
      IF (allocated(idxinq )) deallocate(idxinq )
      IF (allocated(addrinq)) deallocate(addrinq)
      IF (allocated(mask   )) deallocate(mask   )

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_worker) THEN

         IF (numbasin > 0) THEN
            allocate (rlon_b(numbasin))
            allocate (rlat_b(numbasin))
            CALL landelm%get_lonlat_radian (rlon_b, rlat_b)
         ENDIF

         CALL allocate_neighbour_data (rlon_nb)
         CALL allocate_neighbour_data (rlat_nb)

         CALL retrieve_neighbour_data (rlon_b, rlon_nb)
         CALL retrieve_neighbour_data (rlat_b, rlat_nb)
         
         CALL allocate_neighbour_data (dist_nb)
         
         DO ibasin = 1, numbasin
            DO inb = 1, num_nb(ibasin)
               dist_nb(ibasin)%val(inb) = 1.0e3 * arclen ( &
                  rlat_b (ibasin),          rlon_b (ibasin), &
                  rlat_nb(ibasin)%val(inb), rlon_nb(ibasin)%val(inb))
            ENDDO
         ENDDO
         
         IF (numbasin > 0) THEN
            allocate (area_b(numbasin))
            allocate (elva_b(numbasin))
            DO ibasin = 1, numbasin
               area_b(ibasin) = sum(drainagenetwork(ibasin)%area)
               elva_b(ibasin) = sum(drainagenetwork(ibasin)%area * drainagenetwork(ibasin)%elva) / area_b(ibasin)
            ENDDO
         ENDIF
         
         CALL allocate_neighbour_data (area_nb)
         CALL retrieve_neighbour_data (area_b, area_nb)
         
         CALL allocate_neighbour_data (elva_nb)
         CALL retrieve_neighbour_data (elva_b, elva_nb)
         
         CALL allocate_neighbour_data (slope_nb)
         DO ibasin = 1, numbasin
            DO inb = 1, num_nb(ibasin)
               slope_nb(ibasin)%val(inb) = &
                  abs(elva_nb(ibasin)%val(inb) - elva_b(ibasin)) / dist_nb(ibasin)%val(inb)
            ENDDO
         ENDDO
         
         CALL allocate_neighbour_data (theta_a_nb)
         CALL allocate_neighbour_data (zwt_nb    )
         CALL allocate_neighbour_data (Ks_nb     )

         IF (allocated(rlon_b)) deallocate(rlon_b)
         IF (allocated(rlat_b)) deallocate(rlat_b)
         
         IF (allocated(rlon_nb)) deallocate(rlon_nb)
         IF (allocated(rlat_nb)) deallocate(rlat_nb)
         
      ENDIF

   END SUBROUTINE ssrf_init

   ! ----------
   SUBROUTINE retrieve_neighbour_data (vec_in, nbdata)

      USE precision
      USE spmd_task
      USE mod_mesh, only : numelm
      IMPLICIT NONE

      REAL(r8), intent(inout) :: vec_in (:)
      TYPE(pointer_real8_1d)  :: nbdata (:)

      ! Local Variables
      LOGICAL, allocatable :: smask(:), rmask(:)
      INTEGER, allocatable :: req_send(:), req_recv(:)
      TYPE(pointer_real8_1d), allocatable :: sbuff(:), rbuff(:)
      INTEGER :: numbasin, iwork, ibasin, inb, iloc

      IF (p_is_worker) THEN
      
         numbasin = numelm

         DO ibasin = 1, numbasin
            DO inb = 1, num_nb(ibasin)
               IF (addr_nb(ibasin)%val(1,inb)== -1) THEN
                  iloc = addr_nb(ibasin)%val(2,inb)
                  nbdata(ibasin)%val(inb) = vec_in(iloc)
               ENDIF
            ENDDO
         ENDDO

#ifdef USEMPI
         CALL mpi_barrier (p_comm_worker, p_err)

         allocate (smask    (0:p_np_worker-1))
         allocate (req_send (0:p_np_worker-1))
         allocate (sbuff    (0:p_np_worker-1))
         smask(:) = .false.

         DO iwork = 0, p_np_worker-1
            IF (sendaddr(iwork)%ndata > 0) THEN
               smask(iwork) = .true.

               allocate (sbuff(iwork)%val (sendaddr(iwork)%ndata))
               sbuff(iwork)%val = vec_in(sendaddr(iwork)%ibsn)

               CALL mpi_isend(sbuff(iwork)%val, sendaddr(iwork)%ndata, MPI_REAL8, &
                  p_address_worker(iwork), 101, p_comm_glb, req_send(iwork), p_err)
            ENDIF
         ENDDO
         
         allocate (rmask    (0:p_np_worker-1))
         allocate (req_recv (0:p_np_worker-1))
         allocate (rbuff    (0:p_np_worker-1))
         rmask(:) = .false.

         DO iwork = 0, p_np_worker-1
            IF (recvaddr(iwork)%ndata > 0) THEN
               rmask(iwork) = .true.

               allocate (rbuff(iwork)%val (recvaddr(iwork)%ndata))

               CALL mpi_irecv(rbuff(iwork)%val, recvaddr(iwork)%ndata, MPI_REAL8, &
                  p_address_worker(iwork), 101, p_comm_glb, req_recv(iwork), p_err)
            ENDIF
         ENDDO

         IF (any(rmask)) THEN

            CALL mpi_waitall(count(rmask), pack(req_recv,rmask), MPI_STATUSES_IGNORE, p_err)

            DO ibasin = 1, numbasin
               DO inb = 1, num_nb(ibasin)
                  IF (addr_nb(ibasin)%val(1,inb) >= 0) THEN
                     iwork = addr_nb(ibasin)%val(1,inb)
                     iloc  = addr_nb(ibasin)%val(2,inb)
                     nbdata(ibasin)%val(inb) = rbuff(iwork)%val(iloc)
                  ENDIF
               ENDDO
            ENDDO
         ENDIF

         IF (any(smask)) THEN
            CALL mpi_waitall(count(smask), pack(req_send,smask), MPI_STATUSES_IGNORE, p_err)
         ENDIF

         IF (allocated(smask)) deallocate(smask)
         IF (allocated(rmask)) deallocate(rmask)
         
         IF (allocated(req_send)) deallocate(req_send)
         IF (allocated(req_recv)) deallocate(req_recv)
         
         IF (allocated(sbuff)) deallocate(sbuff)
         IF (allocated(rbuff)) deallocate(rbuff)

         CALL mpi_barrier (p_comm_worker, p_err)
#endif

      ENDIF

   END SUBROUTINE retrieve_neighbour_data

   ! ---
   SUBROUTINE allocate_neighbour_data (nbdata)
      
      USE mod_mesh, only : numelm
      IMPLICIT NONE

      TYPE(pointer_real8_1d), allocatable :: nbdata(:)
      INTEGER :: ibasin
            
      IF (numelm > 0) THEN
         allocate (nbdata(numelm))
         DO ibasin = 1, numelm
            IF (num_nb(ibasin) > 0) THEN
               allocate (nbdata(ibasin)%val (num_nb(ibasin)))
            ENDIF
         ENDDO
      ENDIF 

   END SUBROUTINE allocate_neighbour_data 

   ! ----------
   SUBROUTINE ssrf_final ()

      IMPLICIT NONE
      INTEGER :: i

      IF (allocated(num_nb)) deallocate(num_nb)
      IF (allocated(idxbsn_nb )) deallocate(idxbsn_nb )
      IF (allocated(addr_nb   )) deallocate(addr_nb   )
      IF (allocated(dist_nb   )) deallocate(dist_nb   )
      IF (allocated(lenbdr_nb )) deallocate(lenbdr_nb )
      IF (allocated(elva_nb   )) deallocate(elva_nb   )
      IF (allocated(slope_nb  )) deallocate(slope_nb  )
      IF (allocated(elva_b    )) deallocate(elva_b    )
      IF (allocated(area_b    )) deallocate(area_b    )
      IF (allocated(area_nb   )) deallocate(area_nb   )
      
      IF (allocated(theta_a_nb)) deallocate(theta_a_nb)
      IF (allocated(zwt_nb    )) deallocate(zwt_nb    )
      IF (allocated(Ks_nb     )) deallocate(Ks_nb     )
      
      DO i = lbound(recvaddr,1), ubound(recvaddr,1)
         IF (allocated(recvaddr(i)%bindx)) deallocate(recvaddr(i)%bindx)
         IF (allocated(recvaddr(i)%ibsn )) deallocate(recvaddr(i)%ibsn )
      ENDDO 

      DO i = lbound(sendaddr,1), ubound(sendaddr,1)
         IF (allocated(sendaddr(i)%bindx)) deallocate(sendaddr(i)%bindx)
         IF (allocated(sendaddr(i)%ibsn )) deallocate(sendaddr(i)%ibsn )
      ENDDO
         
      IF (allocated(recvaddr)) deallocate(recvaddr)
      IF (allocated(sendaddr)) deallocate(sendaddr)

   END SUBROUTINE ssrf_final

END MODULE mod_ssrf_network
#endif
