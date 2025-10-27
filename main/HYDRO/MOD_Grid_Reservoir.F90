#include <define.h>

#ifdef GridRiverLakeFlow
MODULE MOD_Grid_Reservoir
!-----------------------------------------------------------------------
! DESCRIPTION:
!
!    Reservoir module in gridded mesh.
!
! Created by Shupeng Zhang, Oct 2025
!-----------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Grid_RiverLakeNetwork, only : numucat, ucat_ucid, lake_type

   integer :: totalnumresv
   integer :: numresv
   integer,  allocatable :: ucat2resv   (:)
#ifdef USEMPI
   integer :: sum_numresv_all
   integer,  allocatable :: numresv_all (:)
   integer,  allocatable :: numresv_dsp (:)
   integer,  allocatable :: loc2all_all (:)
#endif

   ! parameters
   integer,  allocatable :: dam_GRAND_ID  (:)  ! GRAND dam ID

   integer,  allocatable :: dam_build_year(:)  ! year in which the dam/barrier was built

   real(r8), allocatable :: volresv_total (:)  ! total reservoir volume      [m^3]
   real(r8), allocatable :: volresv_emerg (:)  ! emergency reservoir volume  [m^3]
   real(r8), allocatable :: volresv_adjust(:)  ! adjustment reservoir volume [m^3]
   real(r8), allocatable :: volresv_normal(:)  ! normal reservoir volume     [m^3]

   real(r8), allocatable :: qresv_flood   (:)  ! flood reservoir outflow      [m^3/s]
   real(r8), allocatable :: qresv_adjust  (:)  ! adjustment reservoir outflow [m^3/s]
   real(r8), allocatable :: qresv_normal  (:)  ! normal reservoir outflow     [m^3/s]

   ! fluxes
   real(r8), allocatable :: qresv_in      (:)  ! reservoir inflow  [m^3/s]
   real(r8), allocatable :: qresv_out     (:)  ! reservoir outflow [m^3/s]

   ! -- PUBLIC SUBROUTINEs --
   PUBLIC :: reservoir_init
   PUBLIC :: reservoir_operation
   PUBLIC :: reservoir_gather_var
   PUBLIC :: reservoir_final

CONTAINS

   ! -------
   SUBROUTINE reservoir_init ( )

   USE MOD_SPMD_Task
   USE MOD_Namelist, only: DEF_ReservoirPara_file
   USE MOD_NetCDFSerial
   USE MOD_Utils

   IMPLICIT NONE

   ! Local variables
   character(len=256) :: parafile

   integer,  allocatable :: dam_seq(:), order(:), loc2all(:)
   real(r8), allocatable :: rcache (:)
   integer,  allocatable :: icache (:)

   integer :: i, iloc, irsv


      parafile = DEF_ReservoirPara_file

      IF (p_is_master) THEN
         CALL ncio_read_serial (parafile, 'dam_GRAND_ID', dam_GRAND_ID)
      ENDIF

      CALL ncio_read_bcast_serial (parafile, 'dam_seq', dam_seq)

      totalnumresv = size(dam_seq)

      IF (p_is_worker) THEN

         allocate (order (totalnumresv))
         order = (/(i, i = 1, totalnumresv)/)

         CALL quicksort (totalnumresv, dam_seq, order)

         allocate (ucat2resv (numucat))
         allocate (loc2all   (numucat))

         numresv = 0
         DO i = 1, numucat
            iloc = find_in_sorted_list1 (ucat_ucid(i), totalnumresv, dam_seq)
            IF (iloc > 0) THEN
               numresv = numresv + 1
               lake_type(i) = 2
               ucat2resv(i) = numresv
               loc2all  (numresv) = order(iloc)
            ENDIF
         ENDDO

#ifdef USEMPI
         allocate (numresv_all (0:p_np_worker-1))
         CALL mpi_allgather (numresv, 1, MPI_INTEGER, &
            numresv_all, 1, MPI_INTEGER, p_comm_worker, p_err)

         sum_numresv_all = sum(numresv_all)

         IF (sum_numresv_all > 0) THEN

            IF (p_iam_worker == p_root) THEN
               allocate (numresv_dsp (0:p_np_worker-1))
               allocate (loc2all_all (sum_numresv_all))

               numresv_dsp(0) = 0
               DO i = 1, p_np_worker-1
                  numresv_dsp(i) = numresv_dsp(i-1) + numresv_all(i-1)
               ENDDO
            ENDIF

            CALL mpi_gatherv (loc2all, numresv, MPI_INTEGER, &
               loc2all_all, numresv_all, numresv_dsp, MPI_INTEGER, p_root, p_comm_worker, p_err)

         ENDIF
#endif

         IF (numresv > 0) THEN

            allocate (dam_build_year (numresv))

            allocate (volresv_total  (numresv))
            allocate (volresv_emerg  (numresv))
            allocate (volresv_adjust (numresv))
            allocate (volresv_normal (numresv))

            allocate (qresv_flood    (numresv))
            allocate (qresv_adjust   (numresv))
            allocate (qresv_normal   (numresv))

            allocate (qresv_in       (numresv))
            allocate (qresv_out      (numresv))

         ENDIF

      ENDIF

      CALL ncio_read_bcast_serial (parafile, 'dam_year', icache)
      IF (p_is_worker .and. (numresv > 0)) THEN
         dam_build_year = icache(loc2all(1:numresv))
      ENDIF

      CALL ncio_read_bcast_serial (parafile, 'dam_TotalVol_mcm', rcache)
      IF (p_is_worker .and. (numresv > 0)) THEN
         volresv_total = rcache(loc2all(1:numresv))*1.e6
      ENDIF

      CALL ncio_read_bcast_serial (parafile, 'dam_ConVol_mcm', rcache)
      IF (p_is_worker .and. (numresv > 0)) THEN
         volresv_normal = rcache(loc2all(1:numresv))*1.e6
      ENDIF

      CALL ncio_read_bcast_serial (parafile, 'dam_Qn', rcache)
      IF (p_is_worker .and. (numresv > 0)) THEN
         qresv_normal = rcache(loc2all(1:numresv))
      ENDIF

      CALL ncio_read_bcast_serial (parafile, 'dam_Qf', rcache)
      IF (p_is_worker .and. (numresv > 0)) THEN
         qresv_flood = rcache(loc2all(1:numresv))
      ENDIF


      IF (p_is_worker) THEN
         DO irsv = 1, numresv
            volresv_emerg (irsv) = volresv_total(irsv) * 0.94
            volresv_adjust(irsv) = volresv_total(irsv) * 0.77
            volresv_normal(irsv) = min(volresv_total(irsv)*0.7, volresv_normal(irsv))
            qresv_adjust  (irsv) = (qresv_normal(irsv) + qresv_flood(irsv)) * 0.5
         ENDDO
      ENDIF

      IF (allocated(dam_seq)) deallocate(dam_seq)
      IF (allocated(order  )) deallocate(order  )
      IF (allocated(loc2all)) deallocate(loc2all)
      IF (allocated(rcache )) deallocate(rcache )
      IF (allocated(icache )) deallocate(icache )

   END SUBROUTINE reservoir_init


   SUBROUTINE reservoir_operation (method, irsv, qin, vol, qout)

   IMPLICIT NONE
   integer,  intent(in)  :: method
   integer,  intent(in)  :: irsv
   real(r8), intent(in)  :: qin, vol
   real(r8), intent(out) :: qout

   ! local variables
   real(r8) :: q1

      IF (method == 1) THEN
         ! *** Reference ***
         ! [1] Mizuki Funato, Dai Yamazaki, Dung Trung Vu.
         ! Development of an Improved Reservoir Operation Scheme for Global Flood Modeling (CaMa-Flood v4.20).
         ! ESS Open Archive . October 24, 2024.

         IF (vol > volresv_emerg(irsv)) THEN
            qout = max(qin, qresv_flood(irsv))
         ELSEIF (vol > volresv_adjust(irsv)) THEN
            qout = qresv_adjust(irsv) + (qresv_flood(irsv)-qresv_adjust(irsv)) &
               * ((vol-volresv_adjust(irsv))/(volresv_emerg(irsv)-volresv_adjust(irsv)))**0.1
            IF (qin > qresv_flood(irsv)) THEN
               q1 = qresv_normal(irsv) + (qin-qresv_normal(irsv)) &
                  * (vol-volresv_normal(irsv))/(volresv_emerg(irsv)-volresv_normal(irsv))
               qout = max(q1, qout)
            ENDIF
         ELSEIF (vol > volresv_normal(irsv)) THEN
            qout = qresv_normal(irsv) + (qresv_adjust(irsv)-qresv_normal(irsv)) &
               * ((vol-volresv_normal(irsv))/(volresv_adjust(irsv)-volresv_normal(irsv)))**3.
         ELSE
            qout = (vol/volresv_normal(irsv))**0.5 * qresv_normal(irsv)
         ENDIF

      ENDIF

   END SUBROUTINE reservoir_operation


   SUBROUTINE reservoir_gather_var (varin, varout)

   USE MOD_SPMD_Task
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   real(r8), intent(in)  :: varin  (:)
   real(r8), intent(out) :: varout (:)

   ! local variables
   integer :: irsv
   real(r8), allocatable :: varall (:)

      IF (p_is_worker) THEN

         IF (totalnumresv == 0) RETURN

#ifdef USEMPI
         IF (sum_numresv_all == 0) RETURN

         IF (p_iam_worker == p_root) THEN
            allocate (varall (sum_numresv_all))
         ENDIF

         CALL mpi_gatherv (varin, numresv, MPI_REAL8, &
            varall, numresv_all, numresv_dsp, MPI_REAL8, p_root, p_comm_worker, p_err)

         IF (p_iam_worker == p_root) THEN

            varout(:) = spval
            DO irsv = 1, sum_numresv_all
               varout(loc2all_all(irsv)) = varall(irsv)
            ENDDO

            deallocate (varall)
         ENDIF

      ENDIF

      IF (p_iam_worker == p_root) THEN
         CALL mpi_send (varout, totalnumresv, MPI_REAL8, p_address_master, &
            mpi_tag_data, p_comm_glb, p_err)
      ENDIF
      IF (p_is_master) THEN
         CALL mpi_recv (varout, totalnumresv, MPI_REAL8, p_address_worker(p_root), &
            mpi_tag_data, p_comm_glb, p_stat, p_err)
      ENDIF
#else
      varout(:) = spval
      DO irsv = 1, numresv
         varout(loc2all(irsv)) = varin(irsv)
      ENDDO
#endif

   END SUBROUTINE reservoir_gather_var


   SUBROUTINE reservoir_final ()

   IMPLICIT NONE

      IF (allocated(ucat2resv     )) deallocate (ucat2resv     )
#ifdef USEMPI
      IF (allocated(numresv_all   )) deallocate (numresv_all   )
      IF (allocated(numresv_dsp   )) deallocate (numresv_dsp   )
      IF (allocated(loc2all_all   )) deallocate (loc2all_all   )
#endif

      IF (allocated(dam_GRAND_ID  )) deallocate (dam_GRAND_ID  )
      IF (allocated(dam_build_year)) deallocate (dam_build_year)

      IF (allocated(volresv_total )) deallocate (volresv_total )
      IF (allocated(volresv_emerg )) deallocate (volresv_emerg )
      IF (allocated(volresv_adjust)) deallocate (volresv_adjust)
      IF (allocated(volresv_normal)) deallocate (volresv_normal)

      IF (allocated(qresv_flood   )) deallocate (qresv_flood   )
      IF (allocated(qresv_adjust  )) deallocate (qresv_adjust  )
      IF (allocated(qresv_normal  )) deallocate (qresv_normal  )

      IF (allocated(qresv_in      )) deallocate (qresv_in      )
      IF (allocated(qresv_out     )) deallocate (qresv_out     )

   END SUBROUTINE reservoir_final

END MODULE MOD_Grid_Reservoir
#endif
