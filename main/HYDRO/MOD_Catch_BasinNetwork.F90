#include <define.h>

#ifdef CatchLateralFlow
MODULE MOD_Catch_BasinNetwork
!--------------------------------------------------------------------------------
! DESCRIPTION:
!
! Created by Shupeng Zhang, Feb 2025
!--------------------------------------------------------------------------------

   USE MOD_Pixelset
   USE MOD_WorkerPushData
   IMPLICIT NONE

   ! -- instances --
   integer :: numbasin
   integer, allocatable :: basinindex(:)

   integer :: numbsnhru
   type(subset_type) :: basin_hru

   integer :: numrivmth
   integer, allocatable :: rivermouth(:)

   integer :: numlake, numresv
   integer, allocatable :: lake_id  (:)
   integer, allocatable :: lake_type(:)  ! lake type:
                                         !   0: not lake;  1: natural lake;
                                         !   2: reservoir; 3: controlled lake.
   integer, allocatable :: bsn2lake (:)
   integer, allocatable :: bsn2resv (:)

   type(worker_pushdata_type) :: push_elm2bsn
   type(worker_pushdata_type) :: push_bsn2elm

   type(worker_pushdata_type) :: push_elmhru2bsnhru
   type(worker_pushdata_type) :: push_bsnhru2elmhru

CONTAINS

   ! ----------
   SUBROUTINE build_basin_network ()

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Utils
   USE MOD_NetCDFSerial
   USE MOD_Mesh,    only: numelm
   USE MOD_LandElm, only: landelm
   USE MOD_LandHRU, only: elm_hru
   IMPLICIT NONE

   ! Local Variables
   character(len=256)   :: basin_file
   integer, allocatable :: basindown(:)

   integer, allocatable :: nups_nst(:), iups_nst(:), nups_all(:), b_up2down(:), orderbsn(:)

   integer , allocatable :: nb_rs(:), iwrk_rs(:), nwrk_rs(:), nave_rs(:), nb_wrk(:)
   real(r8), allocatable :: wtbsn(:), wt_rs  (:), wt_wrk (:)

   integer  :: totalnumbasin, ibasin, nbasin, iriv
   integer  :: iworker, iwrkdsp, mesg(2), isrc, nrecv
   integer  :: iloc, i, j, ithis
   real(r8) :: sumwt

   integer, allocatable :: bindex(:), addrbasin(:), elmindex(:)

   ! lake and reservoir
   character(len=256)   :: lake_info_file
   integer, allocatable :: all_lake_id (:), all_lake_type (:), order(:), icache(:)
   integer, allocatable :: lake_id_resv(:), lake_type_resv(:)
   integer :: ilake, iresv



#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      basin_file = DEF_CatchmentMesh_data

      ! read in parameters from file.
      IF (p_is_master) THEN
         CALL ncio_read_serial (basin_file, 'basin_downstream', basindown)
         totalnumbasin = size(basindown)
      ENDIF

#ifdef USEMPI
      ! divide basins into groups and assign to workers
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
               iups_nst(i) = -1

               j = basindown(i)
               DO WHILE (j > 0)

                  iups_nst(j) = iups_nst(j) + 1

                  IF (iups_nst(j) == nups_nst(j)) THEN
                     ithis = ithis + 1
                     b_up2down(ithis) = j
                     iups_nst(j) = -1

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

         deallocate (basindown)
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


      ! send basin index to workers
      IF (p_is_master) THEN

         allocate(basinindex (totalnumbasin))
         basinindex = (/(i, i = 1, totalnumbasin)/)

         DO iworker = 0, p_np_worker-1

            nbasin = count(addrbasin == p_address_worker(iworker))
            CALL mpi_send (nbasin, 1, MPI_INTEGER, p_address_worker(iworker), mpi_tag_mesg, p_comm_glb, p_err)

            IF (nbasin > 0) THEN
               allocate (bindex (nbasin))

               bindex = pack(basinindex, mask = (addrbasin == p_address_worker(iworker)))
               CALL mpi_send (bindex, nbasin, MPI_INTEGER, p_address_worker(iworker), &
                  mpi_tag_data, p_comm_glb, p_err)

               deallocate (bindex)
            ENDIF

         ENDDO

         deallocate (basinindex)

      ELSEIF (p_is_worker) THEN

         CALL mpi_recv (numbasin, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

         IF (numbasin > 0) THEN

            allocate (basinindex (numbasin))
            CALL mpi_recv (basinindex, numbasin, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_stat, p_err)

         ENDIF

      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      numbasin = totalnumbasin
      allocate(basinindex (totalnumbasin))
      basinindex = (/(i, i = 1, totalnumbasin)/)
#endif


      IF (p_is_worker) THEN
         IF (numelm > 0) THEN
            allocate (elmindex (numelm))
            elmindex = landelm%eindex
         ENDIF
      ENDIF

      CALL build_worker_pushdata (numelm, elmindex, numbasin, basinindex, push_elm2bsn)
      CALL build_worker_pushdata (numbasin, basinindex, numelm, elmindex, push_bsn2elm)

      CALL build_worker_pushdata_subset ( &
         numelm, numbasin, push_elm2bsn, elm_hru, push_elmhru2bsnhru, basin_hru, numbsnhru)

      CALL build_worker_pushdata_subset ( &
         numbasin, numelm, push_bsn2elm, basin_hru, push_bsnhru2elmhru)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (allocated(addrbasin)) deallocate(addrbasin)
      IF (allocated(elmindex )) deallocate(elmindex )


      IF (p_is_master) THEN

         lake_info_file = DEF_CatchmentMesh_data
         CALL ncio_read_serial (lake_info_file, 'lake_id', all_lake_id)

         lake_info_file = trim(DEF_dir_runtime)//'/HydroLAKES_Reservoir.nc'
         CALL ncio_read_serial (lake_info_file, 'hylak_id',  lake_id_resv  )
         CALL ncio_read_serial (lake_info_file, 'lake_type', lake_type_resv)

         allocate (order (size(lake_id_resv)))
         order = (/(iresv, iresv=1,size(lake_id_resv))/)

         CALL quicksort (size(lake_id_resv), lake_id_resv, order)

         lake_type_resv = lake_type_resv(order)

         allocate (all_lake_type (size(all_lake_id)))
         all_lake_type(:) = 0

         DO ibasin = 1, size(all_lake_id)
            IF (all_lake_id(ibasin) > 0) THEN
               all_lake_type(ibasin) = 1
               iresv = find_in_sorted_list1 (all_lake_id(ibasin), size(lake_id_resv), lake_id_resv)
               IF (iresv > 0) THEN
                  all_lake_type(ibasin) = lake_type_resv(iresv)
               ENDIF
            ENDIF
         ENDDO

         deallocate (lake_id_resv  )
         deallocate (lake_type_resv)
         deallocate (order)

      ENDIF

      IF (p_is_worker) THEN
         IF (numbasin > 0) THEN
            allocate (lake_id   (numbasin));   lake_id  (:) = 0
            allocate (lake_type (numbasin));   lake_type(:) = 0
         ENDIF
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_worker) THEN
         mesg = (/p_iam_glb, numbasin/)
         CALL mpi_send (mesg, 2, MPI_INTEGER, p_address_master, &
            mpi_tag_mesg, p_comm_glb, p_err)

         IF (numbasin > 0) THEN
            CALL mpi_send (basinindex, numbasin, MPI_INTEGER, &
               p_address_master, mpi_tag_data, p_comm_glb, p_err)

            CALL mpi_recv (lake_id, numbasin, MPI_INTEGER, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            CALL mpi_recv (lake_type, numbasin, MPI_INTEGER, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)
         ENDIF
      ENDIF

      IF (p_is_master) THEN
         DO iworker = 0, p_np_worker-1

            CALL mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
               mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc   = mesg(1)
            nbasin = mesg(2)
            IF (nbasin > 0) THEN

               allocate(bindex (nbasin))
               CALL mpi_recv (bindex, nbasin, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               allocate(icache (nbasin))

               icache = all_lake_id (bindex)
               CALL mpi_send (icache, nbasin, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_err)

               icache = all_lake_type(bindex)
               CALL mpi_send (icache, nbasin, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_err)

               deallocate (bindex)
               deallocate (icache)
            ENDIF

         ENDDO
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      lake_id   = all_lake_id  (basinindex)
      lake_type = all_lake_type(basinindex)
#endif

      IF (p_is_worker) THEN
         IF (numbasin > 0) THEN
            numlake = count(lake_type /= 0)
            numresv = count(lake_type >= 2)
         ELSE
            numlake = 0
            numresv = 0
         ENDIF

         IF (numlake > 0) THEN
            allocate (bsn2lake (numbasin))
            ilake = 0
            DO ibasin = 1, numbasin
               IF (lake_type(ibasin) /= 0) THEN
                  ilake = ilake + 1
                  bsn2lake(ibasin) = ilake
               ENDIF
            ENDDO
         ENDIF

         IF (numresv > 0) THEN
            allocate (bsn2resv (numbasin))
            iresv = 0
            DO ibasin = 1, numbasin
               IF (lake_type(ibasin) >= 2) THEN
                  iresv = iresv + 1
                  bsn2resv(ibasin) = iresv
               ENDIF
            ENDDO
         ENDIF
      ENDIF

      IF (allocated(all_lake_id  )) deallocate(all_lake_id  )
      IF (allocated(all_lake_type)) deallocate(all_lake_type)


   END SUBROUTINE build_basin_network

   ! ---------
   SUBROUTINE basin_network_final ()

   IMPLICIT NONE

      IF (allocated(basinindex)) deallocate(basinindex)
      IF (allocated(rivermouth)) deallocate(rivermouth)

      IF (allocated(lake_id   )) deallocate(lake_id   )
      IF (allocated(lake_type )) deallocate(lake_type )

      IF (allocated(bsn2lake  )) deallocate(bsn2lake  )
      IF (allocated(bsn2resv  )) deallocate(bsn2resv  )

   END SUBROUTINE basin_network_final

END MODULE MOD_Catch_BasinNetwork
#endif
