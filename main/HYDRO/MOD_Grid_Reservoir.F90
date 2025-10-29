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
   USE MOD_DataType

   integer :: totalnumresv
   integer :: numresv
   integer,  allocatable :: ucat2resv   (:)
   type(pointer_int32_1d), allocatable :: resv_data_address (:)


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
   PUBLIC :: reservoir_final

CONTAINS

   ! -------
   SUBROUTINE reservoir_init ( )

   USE MOD_SPMD_Task
   USE MOD_NetCDFSerial
   USE MOD_Utils
   USE MOD_Namelist,              only: DEF_ReservoirPara_file
   USE MOD_Grid_RiverLakeNetwork, only: numucat, ucat_ucid, lake_type

   IMPLICIT NONE

   ! Local variables
   character(len=256) :: parafile

   integer,  allocatable :: dam_seq(:), order(:), loc2all(:)
   real(r8), allocatable :: rcache (:)
   integer,  allocatable :: icache (:)

   integer :: i, iloc, irsv, nresv, iworker


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

      ENDIF

#ifdef USEMPI
      IF (p_is_master) THEN

         allocate (resv_data_address (0:p_np_worker-1))

         DO iworker = 0, p_np_worker-1

            CALL mpi_recv (nresv, 1, MPI_INTEGER, &
               p_address_worker(iworker), mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            IF (nresv > 0) THEN
               allocate (resv_data_address(iworker)%val (nresv))
               CALL mpi_recv (resv_data_address(iworker)%val, nresv, MPI_INTEGER, &
                  p_address_worker(iworker), mpi_tag_data, p_comm_glb, p_stat, p_err)
            ENDIF
         ENDDO

      ELSEIF (p_is_worker) THEN

         CALL mpi_send (numresv, 1, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)

         IF (numresv > 0) THEN
            CALL mpi_send (loc2all(1:numresv), numresv, MPI_INTEGER, p_address_master, &
               mpi_tag_data, p_comm_glb, p_err)
         ENDIF

      ENDIF
#else
      IF (numresv > 0) THEN
         allocate (resv_data_address (0:0))
         allocate (resv_data_address(0)%val (numresv))
         resv_data_address(0)%val = loc2all(1:numresv)
      ENDIF
#endif

      IF (p_is_worker) THEN

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


   SUBROUTINE reservoir_final ()

   IMPLICIT NONE

      IF (allocated(ucat2resv        )) deallocate (ucat2resv        )
      IF (allocated(resv_data_address)) deallocate (resv_data_address)

      IF (allocated(dam_GRAND_ID     )) deallocate (dam_GRAND_ID     )
      IF (allocated(dam_build_year   )) deallocate (dam_build_year   )

      IF (allocated(volresv_total    )) deallocate (volresv_total    )
      IF (allocated(volresv_emerg    )) deallocate (volresv_emerg    )
      IF (allocated(volresv_adjust   )) deallocate (volresv_adjust   )
      IF (allocated(volresv_normal   )) deallocate (volresv_normal   )

      IF (allocated(qresv_flood      )) deallocate (qresv_flood      )
      IF (allocated(qresv_adjust     )) deallocate (qresv_adjust     )
      IF (allocated(qresv_normal     )) deallocate (qresv_normal     )

      IF (allocated(qresv_in         )) deallocate (qresv_in         )
      IF (allocated(qresv_out        )) deallocate (qresv_out        )

   END SUBROUTINE reservoir_final

END MODULE MOD_Grid_Reservoir
#endif
