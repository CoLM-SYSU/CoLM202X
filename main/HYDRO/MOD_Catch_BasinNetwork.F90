#include <define.h>

#ifdef CatchLateralFlow
MODULE MOD_Catch_BasinNetwork
!--------------------------------------------------------------------------------
! DESCRIPTION:
! 
! Created by Shupeng Zhang, Feb 2025
!--------------------------------------------------------------------------------

   USE MOD_Pixelset
   IMPLICIT NONE

   ! -- instances --
   integer :: numbasin
   integer, allocatable :: basinindex(:)

   integer :: numbsnhru
   type(subset_type) :: basin_hru
   
   integer :: numrivmth
   integer, allocatable :: rivermouth(:)

   ! -- communications --
   type :: basin_pushdata_type
      ! data is on the same processor
      integer :: nself
      integer, allocatable :: iself (:)
      ! data is on other processors
      integer :: nproc
      integer, allocatable :: paddr (:)
      integer, allocatable :: ndata (:)
      integer, allocatable :: ipush (:)
   CONTAINS
      final :: basin_pushdata_free_mem
   END type basin_pushdata_type
   
   type(basin_pushdata_type), target :: iam_bsn
   type(basin_pushdata_type), target :: iam_elm

   ! -- public subroutines --
   interface worker_push_data
      MODULE procedure worker_push_data_real8
      MODULE procedure worker_push_data_int32
   END interface worker_push_data

   PUBLIC :: worker_push_subset_data

CONTAINS
   
   ! ----------
   SUBROUTINE build_basin_network ()

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_NetCDFSerial
   USE MOD_Mesh
   USE MOD_LandElm
   USE MOD_Utils
   IMPLICIT NONE

   ! Local Variables
   character(len=256)   :: basin_file
   integer, allocatable :: basindown(:), nhru_all(:), nhru_in_bsn(:)
   
   integer, allocatable :: nups_nst(:), iups_nst(:), nups_all(:), b_up2down(:), orderbsn(:)
   
   integer , allocatable :: nb_rs(:), iwrk_rs(:), nwrk_rs(:), nave_rs(:), nb_wrk(:)
   real(r8), allocatable :: wtbsn(:), wt_rs  (:), wt_wrk (:)

   integer  :: totalnumbasin, ibasin, nbasin, iriv
   integer  :: iworker, iwrkdsp, mesg(2), isrc, nrecv, idata, ndatall
   integer  :: ip, iloc, ielm, i, j, ithis
   real(r8) :: sumwt
   
   integer, allocatable :: eindex(:), bindex(:), addrelm(:), addrbasin(:)
   integer, allocatable :: paddr (:), icache(:)

   integer, allocatable :: basin_sorted(:), element_sorted(:)
   integer, allocatable :: basin_order (:), element_order (:)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      basin_file = DEF_CatchmentMesh_data 

      ! step 1: read in parameters from file.
      IF (p_is_master) THEN
         CALL ncio_read_serial (basin_file, 'basin_downstream', basindown)
         CALL ncio_read_serial (basin_file, 'basin_numhru',     nhru_all )
         totalnumbasin = size(basindown)
      ENDIF

#ifdef USEMPI
      ! 3-1: get address of elements
      IF (p_is_master) THEN
         
         allocate (addrelm (totalnumbasin))

         DO iworker = 0, p_np_worker-1

            CALL mpi_recv (mesg(1:2), 2, MPI_INTEGER, MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
            isrc  = mesg(1)
            nrecv = mesg(2)

            IF (nrecv > 0) THEN
               allocate (eindex (nrecv))

               CALL mpi_recv (eindex, nrecv, MPI_INTEGER, isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               addrelm(eindex) = isrc

               deallocate(eindex)
            ENDIF

         ENDDO

      ELSEIF (p_is_worker) THEN
      
         mesg(1:2) = (/p_iam_glb, numelm/)
         CALL mpi_send (mesg(1:2), 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err) 

         IF (numelm > 0) THEN
            allocate (eindex (numelm))

            eindex = landelm%eindex

            CALL mpi_send (eindex, numelm, MPI_INTEGER, p_address_master, mpi_tag_data, p_comm_glb, p_err) 

            deallocate (eindex)
         ENDIF

      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)

      ! 3-2: divide basins into groups and assign to workers
      IF (p_is_master) THEN
         
         IF (ncio_var_exist(basin_file, 'weightbasin')) THEN
            CALL ncio_read_serial (basin_file, 'weightbasin', wtbsn)
         ELSE
            allocate (wtbsn (totalnumbasin));  wtbsn(:) = 1.
         ENDIF

         ! sort basins from up to down, recorded by "b_up2down"

         allocate (nups_nst (totalnumbasin));  nups_nst(:) = 0
         allocate (iups_nst (totalnumbasin));  iups_nst(:) = 0
         allocate (b_up2down(totalnumbasin))
         
         DO i = 1, totalnumbasin
            j = basindown(i)
            IF (j > 0) THEN
               nups_nst(j) = nups_nst(j) + 1
            ENDIF
         ENDDO

         ithis = 0
         DO i = 1, totalnumbasin
            IF (iups_nst(i) == nups_nst(i)) THEN

               ithis = ithis + 1
               b_up2down(ithis) = i

               j = basindown(i)
               DO WHILE (j > 0) 

                  iups_nst(j) = iups_nst(j) + 1

                  IF (iups_nst(j) == nups_nst(j)) THEN
                     ithis = ithis + 1
                     b_up2down(ithis) = j
                     j = basindown(j)
                  ELSE
                     EXIT
                  ENDIF
               ENDDO
            ENDIF
         ENDDO

         deallocate (nups_nst)
         deallocate (iups_nst)

         allocate (rivermouth (totalnumbasin))
         numrivmth = 0
         DO i = totalnumbasin, 1, -1
            j = basindown(b_up2down(i))
            IF (j <= 0) THEN
               numrivmth = numrivmth + 1
               rivermouth(b_up2down(i)) = numrivmth
            ELSE
               rivermouth(b_up2down(i)) = rivermouth(j)
            ENDIF
         ENDDO

         allocate (nb_rs (numrivmth)); nb_rs(:) = 0
         allocate (wt_rs (numrivmth)); wt_rs(:) = 0
         DO i = 1, totalnumbasin
            nb_rs(rivermouth(i)) = nb_rs(rivermouth(i)) + 1 
            wt_rs(rivermouth(i)) = wt_rs(rivermouth(i)) + wtbsn(i)
         ENDDO
         
         sumwt = sum(wt_rs)

         allocate (iwrk_rs (numrivmth))
         allocate (nwrk_rs (numrivmth))
         allocate (nave_rs (numrivmth))

         iwrkdsp = -1
         DO i = 1, numrivmth
            nwrk_rs(i) = floor(wt_rs(i)/sumwt * p_np_worker)
            IF (nwrk_rs(i) > 1) THEN
               
               nave_rs(i) = nb_rs(i) / nwrk_rs(i)
               IF (mod(nb_rs(i), nwrk_rs(i)) /= 0) THEN
                  nave_rs(i) = nave_rs(i) + 1
               ENDIF

               iwrk_rs(i) = iwrkdsp + 1
               iwrkdsp = iwrkdsp + nwrk_rs(i)
            ENDIF
         ENDDO

         allocate (nups_all (totalnumbasin));  nups_all(:) = 1

         DO i = 1, totalnumbasin
            j = basindown(b_up2down(i))
            IF (j > 0) THEN
               nups_all(j) = nups_all(j) + nups_all(b_up2down(i))
            ENDIF
         ENDDO

         allocate (addrbasin (totalnumbasin));  addrbasin(:) = -1
            
         allocate (wt_wrk (0:p_np_worker-1));  wt_wrk(:) = 0
         allocate (nb_wrk (0:p_np_worker-1));  nb_wrk(:) = 0

         allocate (orderbsn(totalnumbasin))
         orderbsn(b_up2down) = (/(i, i = 1, totalnumbasin)/)

         ithis = totalnumbasin
         DO WHILE (ithis > 0)

            i = b_up2down(ithis)

            IF (addrbasin(i) >= 0) THEN
               ithis = ithis - 1
               CYCLE
            ENDIF

            j = basindown(i)
            IF (j > 0) THEN
               IF (addrbasin(j) >= 0) THEN
                  addrbasin(i) = addrbasin(j)
                  ithis = ithis - 1
                  CYCLE
               ENDIF
            ENDIF

            iriv = rivermouth(i)
            IF (nwrk_rs(iriv) > 1) THEN
               iworker = iwrk_rs(iriv)
               IF (nups_all(i) <= nave_rs(iriv)-nb_wrk(iworker)) THEN

                  addrbasin(i) = p_address_worker(iworker)

                  nb_wrk(iworker) = nb_wrk(iworker) + nups_all(i)
                  IF (nb_wrk(iworker) == nave_rs(iriv)) THEN
                     iwrk_rs(iriv) = iwrk_rs(iriv) + 1
                  ENDIF
               
                  j = basindown(i)
                  IF (j > 0) THEN
                     DO WHILE (j > 0) 
                        nups_all(j) = nups_all(j) - nups_all(i)
                        ithis = orderbsn(j)
                        j = basindown(j)
                     ENDDO
                  ELSE
                     ithis = ithis - 1
                  ENDIF 
               ELSE
                  ithis = ithis - 1
               ENDIF
            ELSE
               iworker = minloc(wt_wrk(iwrkdsp+1:p_np_worker-1), dim=1) + iwrkdsp

               addrbasin(i) = p_address_worker(iworker)

               wt_wrk(iworker) = wt_wrk(iworker) + wt_rs(iriv)
               ithis = ithis - 1
            ENDIF

         ENDDO

         deallocate (wtbsn    )
         deallocate (b_up2down)
         deallocate (nups_all )
         deallocate (orderbsn )
         deallocate (nb_rs    )
         deallocate (wt_rs    )
         deallocate (iwrk_rs  )
         deallocate (nwrk_rs  )
         deallocate (nave_rs  )
         deallocate (wt_wrk   )
         deallocate (nb_wrk   )

      ENDIF


      ! 3-3: send basin index to workers
      IF (p_is_master) THEN
      
         allocate(basinindex (totalnumbasin))
         basinindex = (/(i, i = 1, totalnumbasin)/)
         
         DO iworker = 0, p_np_worker-1
            
            nbasin = count(addrbasin == p_address_worker(iworker))
            CALL mpi_send (nbasin, 1, MPI_INTEGER, p_address_worker(iworker), mpi_tag_mesg, p_comm_glb, p_err) 

            IF (nbasin > 0) THEN
               allocate (bindex      (nbasin))
               allocate (icache      (nbasin))
               allocate (nhru_in_bsn (nbasin))

               bindex = pack(basinindex, mask = (addrbasin == p_address_worker(iworker)))
               CALL mpi_send (bindex, nbasin, MPI_INTEGER, p_address_worker(iworker), &
                  mpi_tag_data, p_comm_glb, p_err) 
               
               icache = addrelm(bindex)
               CALL mpi_send (icache, nbasin, MPI_INTEGER, p_address_worker(iworker), &
                  mpi_tag_data, p_comm_glb, p_err) 
               
               icache = basindown(bindex)
               CALL mpi_send (icache, nbasin, MPI_INTEGER, p_address_worker(iworker), &
                  mpi_tag_data, p_comm_glb, p_err) 
               
               nhru_in_bsn = nhru_all(bindex)
               CALL mpi_send (nhru_in_bsn, nbasin, MPI_INTEGER, p_address_worker(iworker), &
                  mpi_tag_data, p_comm_glb, p_err) 

               deallocate (bindex)
               deallocate (icache)
               deallocate (nhru_in_bsn)
            ENDIF
            
         ENDDO

         deallocate (basinindex)

      ELSEIF (p_is_worker) THEN

         CALL mpi_recv (numbasin, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

         IF (numbasin > 0) THEN

            allocate (basinindex (numbasin))
            CALL mpi_recv (basinindex, numbasin, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (addrelm (numbasin))
            CALL mpi_recv (addrelm, numbasin, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)
               
            allocate (basindown (numbasin))
            CALL mpi_recv (basindown, numbasin, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)

            allocate (nhru_in_bsn (numbasin))
            CALL mpi_recv (nhru_in_bsn, numbasin, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)

         ENDIF

      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)

      ! 3-4: send basin index of elements to workers
      IF (p_is_master) THEN
         
         DO iworker = 0, p_np_worker-1

            CALL mpi_recv (mesg(1:2), 2, MPI_INTEGER, MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
            isrc  = mesg(1)
            nrecv = mesg(2)

            IF (nrecv > 0) THEN
               allocate (eindex (nrecv))
               allocate (icache (nrecv))
               
               CALL mpi_recv (eindex, nrecv, MPI_INTEGER, isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               icache = addrbasin(eindex)
               CALL mpi_send (icache, nrecv, MPI_INTEGER, isrc, mpi_tag_data, p_comm_glb, p_err) 

               deallocate(eindex)
               deallocate(icache)
            ENDIF

         ENDDO

      ELSEIF (p_is_worker) THEN
      
         mesg(1:2) = (/p_iam_glb, numelm/)
         CALL mpi_send (mesg(1:2), 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err) 

         IF (numelm > 0) THEN
            allocate (eindex (numelm))
            eindex = landelm%eindex

            CALL mpi_send (eindex, numelm, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_err) 

            allocate(addrbasin (numelm))
            CALL mpi_recv (addrbasin, numelm, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)

            deallocate (eindex)
         ENDIF

      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)

      ! step 4: building push data type
      IF (p_is_worker) THEN
            
         IF (numbasin > 0) THEN
            allocate (basin_sorted (numbasin))
            allocate (basin_order  (numbasin))
            basin_sorted = basinindex
            basin_order  = (/(ibasin, ibasin = 1, numbasin)/)

            CALL quicksort (numbasin, basin_sorted, basin_order)
            
         ENDIF

         IF (numelm > 0) THEN
            allocate (element_sorted (numelm))
            allocate (element_order  (numelm))
            element_sorted = landelm%eindex
            element_order  = (/(ielm, ielm = 1, numelm)/)

            CALL quicksort (numelm, element_sorted, element_order)
            
         ENDIF


         iam_bsn%nself = 0
         iam_bsn%nproc = 0
         iam_elm%nself = 0
         iam_elm%nproc = 0

         IF (numelm > 0) THEN

            allocate (bindex(numelm))
            allocate (paddr (numelm))

            ndatall = 0
            DO ielm = 1, numelm
               IF (addrbasin(ielm) /= p_iam_glb) THEN
                  CALL insert_into_sorted_list2 (int(landelm%eindex(ielm)), addrbasin(ielm), &
                     ndatall, bindex, paddr, iloc)
               ENDIF
            ENDDO

            IF (ndatall > 0) THEN

               DO idata = 1, ndatall
                  IF (idata == 1) THEN
                     iam_elm%nproc = 1
                  ELSEIF (paddr(idata) /= paddr(idata-1)) THEN
                     iam_elm%nproc = iam_elm%nproc + 1
                  ENDIF
               ENDDO

               allocate (iam_elm%paddr (iam_elm%nproc))
               allocate (iam_elm%ndata (iam_elm%nproc))
               allocate (iam_elm%ipush (ndatall))

               DO idata = 1, ndatall

                  IF (idata == 1) THEN
                     ip = 1
                     iam_elm%paddr(ip) = paddr(idata)
                     iam_elm%ndata(ip) = 1
                  ELSEIF (paddr(idata) /= paddr(idata-1)) THEN
                     ip = ip + 1
                     iam_elm%paddr(ip) = paddr(idata)
                     iam_elm%ndata(ip) = 1
                  ELSE
                     iam_elm%ndata(ip) = iam_elm%ndata(ip) + 1
                  ENDIF

                  iloc = find_in_sorted_list1 (bindex(idata), numelm, element_sorted)
                  iam_elm%ipush(idata) = element_order(iloc)

               ENDDO

            ENDIF
            
            deallocate (bindex)
            deallocate (paddr )

         ENDIF


         IF (numbasin > 0) THEN

            iam_bsn%nself = count(addrelm == p_iam_glb)

            IF (iam_bsn%nself > 0) THEN

               allocate (iam_bsn%iself (iam_bsn%nself))

               iam_elm%nself = iam_bsn%nself
               allocate (iam_elm%iself (iam_elm%nself))

               idata = 0
               DO ibasin = 1, numbasin
                  IF (addrelm(ibasin) == p_iam_glb) THEN
                     idata = idata + 1
                     iloc = find_in_sorted_list1 (basinindex(ibasin), numelm, element_sorted)
                     iam_bsn%iself(idata) = ibasin
                     iam_elm%iself(idata) = element_order(iloc)
                  ENDIF
               ENDDO

            ENDIF


            allocate (eindex(numbasin))
            allocate (paddr (numbasin))

            ndatall = 0
            DO ibasin = 1, numbasin
               IF (addrelm(ibasin) /= p_iam_glb) THEN
                  CALL insert_into_sorted_list2 (basinindex(ibasin), addrelm(ibasin), &
                     ndatall, eindex, paddr, iloc)
               ENDIF
            ENDDO

            IF (ndatall > 0) THEN

               DO idata = 1, ndatall
                  IF (idata == 1) THEN
                     iam_bsn%nproc = 1
                  ELSEIF (paddr(idata) /= paddr(idata-1)) THEN
                     iam_bsn%nproc = iam_bsn%nproc + 1
                  ENDIF
               ENDDO

               allocate (iam_bsn%paddr (iam_bsn%nproc))
               allocate (iam_bsn%ndata (iam_bsn%nproc))
               allocate (iam_bsn%ipush (ndatall))

               DO idata = 1, ndatall

                  IF (idata == 1) THEN
                     ip = 1
                     iam_bsn%paddr(ip) = paddr(idata)
                     iam_bsn%ndata(ip) = 1
                  ELSEIF (paddr(idata) /= paddr(idata-1)) THEN
                     ip = ip + 1
                     iam_bsn%paddr(ip) = paddr(idata)
                     iam_bsn%ndata(ip) = 1
                  ELSE
                     iam_bsn%ndata(ip) = iam_bsn%ndata(ip) + 1
                  ENDIF

                  iloc = find_in_sorted_list1 (eindex(idata), numbasin, basin_sorted)
                  iam_bsn%ipush(idata) = basin_order(iloc)

               ENDDO

            ENDIF

            deallocate (eindex)
            deallocate (paddr )

         ENDIF

      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else

      numbasin = numelm
      allocate(basinindex (numbasin))
      basinindex = landelm%eindex

      iam_bsn%nself = numbasin
      allocate(iam_bsn%iself (numbasin))
      iam_bsn%iself = (/(ibasin, ibasin = 1, numbasin)/)

      iam_bsn%nproc = 0

      iam_elm%nself = numelm
      allocate(iam_elm%iself (numelm))
      iam_elm%iself = (/(ielm, ielm = 1, numelm)/)

      iam_elm%nproc = 0

#endif

      IF (p_is_worker) THEN

         numbsnhru = 0

         IF (numbasin > 0) THEN

            numbsnhru = sum(nhru_in_bsn)

            allocate (basin_hru%substt (numbasin))
            allocate (basin_hru%subend (numbasin))

            DO ibasin = 1, numbasin
               IF (ibasin == 1) THEN
                  basin_hru%substt(1) = 1
               ELSE
                  basin_hru%substt(ibasin) = basin_hru%subend(ibasin-1) + 1
               ENDIF
               basin_hru%subend(ibasin) = basin_hru%substt(ibasin) + nhru_in_bsn(ibasin) - 1
            ENDDO
         ENDIF

      ENDIF

#ifdef CoLMDEBUG
      ! check basin network.
      IF (p_is_worker) THEN
         nbasin = 0
         DO ibasin = 1, numbasin
            IF (basindown(ibasin) > 0) THEN
               iloc = find_in_sorted_list1 (basindown(ibasin), numbasin, basin_sorted)
               IF (iloc <= 0) nbasin = nbasin + 1
            ENDIF
         ENDDO

         write(*,'(A,I6,A,I8,A,I8,A)') 'Check basin network: worker ', p_iam_glb, &
            ' has total ', numbasin, ' basins with ', nbasin, ' downstream on other processors.'
      ENDIF
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (allocated(addrbasin     )) deallocate(addrbasin     )
      IF (allocated(addrelm       )) deallocate(addrelm       )
      IF (allocated(basindown     )) deallocate(basindown     )
      IF (allocated(nhru_all      )) deallocate(nhru_all      )  
      IF (allocated(nhru_in_bsn   )) deallocate(nhru_in_bsn   )  
      IF (allocated(basin_sorted  )) deallocate(basin_sorted  )
      IF (allocated(basin_order   )) deallocate(basin_order   )
      IF (allocated(element_sorted)) deallocate(element_sorted)
      IF (allocated(element_order )) deallocate(element_order )

   END SUBROUTINE build_basin_network

   
   ! ----------
   SUBROUTINE worker_push_data_real8 (send_pointer, recv_pointer, vec_send, vec_recv)

   USE MOD_Precision
   USE MOD_SPMD_Task
   IMPLICIT NONE

   type(basin_pushdata_type) :: send_pointer
   type(basin_pushdata_type) :: recv_pointer

   real(r8), intent(in)    :: vec_send(:)
   real(r8), intent(inout) :: vec_recv(:)

   ! Local Variables
   integer :: ndatasend, idest
   integer,  allocatable :: req_send (:)
   real(r8), allocatable :: sendcache(:)

   integer :: ndatarecv, isrc
   integer,  allocatable :: req_recv (:)
   real(r8), allocatable :: recvcache(:)

   integer :: iproc, i, istt, iend, ndata

      IF (p_is_worker) THEN

         IF (send_pointer%nself > 0) THEN
            vec_recv(recv_pointer%iself) = vec_send(send_pointer%iself)
         ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_worker, p_err)

         IF (send_pointer%nproc > 0) THEN

            ndatasend = sum(send_pointer%ndata)
            
            allocate (sendcache(ndatasend))
            sendcache = vec_send(send_pointer%ipush)

            allocate (req_send(send_pointer%nproc))
           
            iend = 0
            DO iproc = 1, send_pointer%nproc
               ndata = send_pointer%ndata(iproc)
               idest = send_pointer%paddr(iproc)
               istt  = iend + 1
               iend  = iend + ndata

               CALL mpi_isend(sendcache(istt:iend), ndata, MPI_REAL8, &
                  idest, 101, p_comm_glb, req_send(iproc), p_err)
            ENDDO

         ENDIF

         IF (recv_pointer%nproc > 0) THEN

            ndatarecv = sum(recv_pointer%ndata)
            
            allocate (recvcache(ndatarecv))
            allocate (req_recv (recv_pointer%nproc))
            
            iend = 0
            DO iproc = 1, recv_pointer%nproc
               ndata = recv_pointer%ndata(iproc)
               isrc  = recv_pointer%paddr(iproc)
               istt  = iend + 1
               iend  = iend + ndata

               CALL mpi_irecv(recvcache(istt:iend), ndata, MPI_REAL8, &
                  isrc, 101, p_comm_glb, req_recv(iproc), p_err)
            ENDDO

         ENDIF

         IF (recv_pointer%nproc > 0) THEN
            CALL mpi_waitall(recv_pointer%nproc, req_recv, MPI_STATUSES_IGNORE, p_err)
            vec_recv(recv_pointer%ipush) = recvcache
         ENDIF

         IF (send_pointer%nproc > 0) THEN
            CALL mpi_waitall(send_pointer%nproc, req_send, MPI_STATUSES_IGNORE, p_err)
         ENDIF

         IF (allocated(req_send )) deallocate(req_send)
         IF (allocated(sendcache)) deallocate(sendcache)
         IF (allocated(req_recv )) deallocate(req_recv)
         IF (allocated(recvcache)) deallocate(recvcache)

         CALL mpi_barrier (p_comm_worker, p_err)
#endif

      ENDIF

   END SUBROUTINE worker_push_data_real8

   ! ----------
   SUBROUTINE worker_push_data_int32 (send_pointer, recv_pointer, vec_send, vec_recv)

   USE MOD_Precision
   USE MOD_SPMD_Task
   IMPLICIT NONE

   type(basin_pushdata_type) :: send_pointer
   type(basin_pushdata_type) :: recv_pointer

   integer, intent(in)    :: vec_send(:)
   integer, intent(inout) :: vec_recv(:)

   ! Local Variables
   integer :: ndatasend, idest
   integer, allocatable :: req_send (:)
   integer, allocatable :: sendcache(:)

   integer :: ndatarecv, isrc
   integer, allocatable :: req_recv (:)
   integer, allocatable :: recvcache(:)

   integer :: iproc, i, istt, iend, ndata

      IF (p_is_worker) THEN

         IF (send_pointer%nself > 0) THEN
            vec_recv(recv_pointer%iself) = vec_send(send_pointer%iself)
         ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_worker, p_err)

         IF (send_pointer%nproc > 0) THEN

            ndatasend = sum(send_pointer%ndata)
            
            allocate (sendcache(ndatasend))
            sendcache = vec_send(send_pointer%ipush)

            allocate (req_send(send_pointer%nproc))
           
            iend = 0
            DO iproc = 1, send_pointer%nproc
               ndata = send_pointer%ndata(iproc)
               idest = send_pointer%paddr(iproc)
               istt  = iend + 1
               iend  = iend + ndata

               CALL mpi_isend(sendcache(istt:iend), ndata, MPI_INTEGER, &
                  idest, 102, p_comm_glb, req_send(iproc), p_err)
            ENDDO

         ENDIF

         IF (recv_pointer%nproc > 0) THEN

            ndatarecv = sum(recv_pointer%ndata)
            
            allocate (recvcache(ndatarecv))
            allocate (req_recv (recv_pointer%nproc))
            
            iend = 0
            DO iproc = 1, recv_pointer%nproc
               ndata = recv_pointer%ndata(iproc)
               isrc  = recv_pointer%paddr(iproc)
               istt  = iend + 1
               iend  = iend + ndata

               CALL mpi_irecv(recvcache(istt:iend), ndata, MPI_INTEGER, &
                  isrc, 102, p_comm_glb, req_recv(iproc), p_err)
            ENDDO

         ENDIF

         IF (recv_pointer%nproc > 0) THEN
            CALL mpi_waitall(recv_pointer%nproc, req_recv, MPI_STATUSES_IGNORE, p_err)
            vec_recv(recv_pointer%ipush) = recvcache
         ENDIF

         IF (send_pointer%nproc > 0) THEN
            CALL mpi_waitall(send_pointer%nproc, req_send, MPI_STATUSES_IGNORE, p_err)
         ENDIF

         IF (allocated(req_send )) deallocate(req_send)
         IF (allocated(sendcache)) deallocate(sendcache)
         IF (allocated(req_recv )) deallocate(req_recv)
         IF (allocated(recvcache)) deallocate(recvcache)

         CALL mpi_barrier (p_comm_worker, p_err)
#endif

      ENDIF

   END SUBROUTINE worker_push_data_int32

   ! ----------
   SUBROUTINE worker_push_subset_data (send_pointer, recv_pointer, &
         subset_send, subset_recv, vec_send, vec_recv)

   USE MOD_Precision
   USE MOD_SPMD_Task
   IMPLICIT NONE

   type(basin_pushdata_type), intent(in) :: send_pointer
   type(basin_pushdata_type), intent(in) :: recv_pointer

   type(subset_type),intent(in) :: subset_send
   type(subset_type),intent(in) :: subset_recv

   real(r8), intent(in)    :: vec_send(:)
   real(r8), intent(inout) :: vec_recv(:)

   ! Local Variables
   integer :: ndatasend, idest, isend, istt_send, iend_send, nreq_send
   integer,  allocatable :: req_send (:)
   real(r8), allocatable :: sendcache(:)

   integer :: ndatarecv, isrc, irecv, istt_recv, iend_recv, nreq_recv
   integer,  allocatable :: req_recv (:)
   real(r8), allocatable :: recvcache(:)

   integer :: iproc, i, istt, iend, isup 
   logical :: has_data


      IF (p_is_worker) THEN

         DO i = 1, send_pointer%nself
            isend = send_pointer%iself(i)
            istt_send = subset_send%substt(isend)
            iend_send = subset_send%subend(isend)

            IF (istt_send <= iend_send) THEN
               irecv = recv_pointer%iself(i)
               istt_recv = subset_recv%substt(irecv)
               iend_recv = subset_recv%subend(irecv)

               vec_recv(istt_recv:iend_recv) = vec_send(istt_send:iend_send)
            ENDIF
         ENDDO

#ifdef USEMPI
         CALL mpi_barrier (p_comm_worker, p_err)

         nreq_send = 0
            
         IF (send_pointer%nproc > 0) THEN

            ndatasend = 0
            DO i = 1, size(send_pointer%ipush)
               isend = send_pointer%ipush(i)
               ndatasend = ndatasend + subset_send%subend(isend) - subset_send%substt(isend) + 1
            ENDDO
               
            IF (ndatasend > 0) THEN

               allocate (sendcache(ndatasend))
               iend = 0
               DO i = 1, size(send_pointer%ipush)
                  isend = send_pointer%ipush(i)
                  istt_send = subset_send%substt(isend)
                  iend_send = subset_send%subend(isend)
                  IF (istt_send <= iend_send) THEN
                     istt = iend + 1
                     iend = iend + iend_send - istt_send + 1
                     sendcache(istt:iend) = vec_send(istt_send:iend_send)
                  ENDIF
               ENDDO

               allocate (req_send(send_pointer%nproc))
           
               isup = 0
               iend = 0
               DO iproc = 1, send_pointer%nproc
                  has_data = .false.
                  DO i = isup+1, isup+send_pointer%ndata(iproc)
                     isend = send_pointer%ipush(i) 
                     istt_send = subset_send%substt(isend)
                     iend_send = subset_send%subend(isend)
                     IF (istt_send <= iend_send) THEN
                        IF (.not. has_data) THEN
                           has_data = .true.
                           istt = iend + 1
                        ENDIF
                        iend = iend + iend_send - istt_send + 1
                     ENDIF
                  ENDDO

                  IF (has_data) THEN
                     nreq_send = nreq_send + 1
                     idest = send_pointer%paddr(iproc)
                     CALL mpi_isend(sendcache(istt:iend), iend-istt+1, MPI_REAL8, &
                        idest, 101, p_comm_glb, req_send(nreq_send), p_err)
                  ENDIF

                  isup = isup + send_pointer%ndata(iproc)
               ENDDO

            ENDIF

         ENDIF

         nreq_recv = 0

         IF (recv_pointer%nproc > 0) THEN

            ndatarecv = 0
            DO i = 1, size(recv_pointer%ipush)
               irecv = recv_pointer%ipush(i)
               ndatarecv = ndatarecv + subset_recv%subend(irecv) - subset_recv%substt(irecv) + 1
            ENDDO
            
            IF (ndatarecv > 0) THEN

               allocate (recvcache(ndatarecv))
               allocate (req_recv (recv_pointer%nproc))
           
               isup = 0
               iend = 0
               DO iproc = 1, recv_pointer%nproc
                  has_data = .false.
                  DO i = isup+1, isup+recv_pointer%ndata(iproc)
                     irecv = recv_pointer%ipush(i) 
                     istt_recv = subset_recv%substt(irecv)
                     iend_recv = subset_recv%subend(irecv)
                     IF (istt_recv <= iend_recv) THEN
                        IF (.not. has_data) THEN
                           has_data = .true.
                           istt = iend + 1
                        ENDIF
                        iend = iend + iend_recv - istt_recv + 1
                     ENDIF
                  ENDDO

                  IF (has_data) THEN
                     nreq_recv = nreq_recv + 1
                     isrc = recv_pointer%paddr(iproc)
                     CALL mpi_irecv(recvcache(istt:iend), iend-istt+1, MPI_REAL8, &
                        isrc, 101, p_comm_glb, req_recv(nreq_recv), p_err)
                  ENDIF

                  isup = isup + recv_pointer%ndata(iproc)
               ENDDO

            ENDIF

         ENDIF

         IF (nreq_recv > 0) THEN

            CALL mpi_waitall(nreq_recv, req_recv(1:nreq_recv), MPI_STATUSES_IGNORE, p_err)

            iend = 0
            DO i = 1, size(recv_pointer%ipush)
               irecv = recv_pointer%ipush(i)
               istt_recv = subset_recv%substt(irecv)
               iend_recv = subset_recv%subend(irecv)
               IF (istt_recv <= iend_recv) THEN
                  istt = iend + 1
                  iend = iend + iend_recv - istt_recv + 1
                  vec_recv(istt_recv:iend_recv) = recvcache(istt:iend)
               ENDIF
            ENDDO
         ENDIF

         IF (nreq_send > 0) THEN
            CALL mpi_waitall(nreq_send, req_send(1:nreq_send), MPI_STATUSES_IGNORE, p_err)
         ENDIF

         IF (allocated(req_send )) deallocate(req_send )
         IF (allocated(sendcache)) deallocate(sendcache)
         IF (allocated(req_recv )) deallocate(req_recv )
         IF (allocated(recvcache)) deallocate(recvcache)

         CALL mpi_barrier (p_comm_worker, p_err)
#endif

      ENDIF

   END SUBROUTINE worker_push_subset_data

   ! ---------
   SUBROUTINE basin_network_final ()

   IMPLICIT NONE 

      IF (allocated(basinindex)) deallocate(basinindex)
      IF (allocated(rivermouth)) deallocate(rivermouth)

   END SUBROUTINE basin_network_final

   ! ---------
   SUBROUTINE basin_pushdata_free_mem (this)
      
   IMPLICIT NONE
   type(basin_pushdata_type) :: this

      IF (allocated(this%iself)) deallocate(this%iself)
      IF (allocated(this%paddr)) deallocate(this%paddr)
      IF (allocated(this%ndata)) deallocate(this%ndata)
      IF (allocated(this%ipush)) deallocate(this%ipush)

   END SUBROUTINE basin_pushdata_free_mem

END MODULE MOD_Catch_BasinNetwork
#endif
