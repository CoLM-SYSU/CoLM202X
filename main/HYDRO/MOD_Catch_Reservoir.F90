#include <define.h>

#ifdef CatchLateralFlow
MODULE MOD_Catch_Reservoir
!-----------------------------------------------------------------------
! DESCRIPTION:
!
!    Reservoir module in catchment mesh.
!
! Created by Shupeng Zhang, July 2025
!-----------------------------------------------------------------------
   USE MOD_Precision

   integer,  allocatable :: lake_type     (:)  ! lake type: 0: not lake;
                                               !            1: natural lake;
                                               !            2: reservoir;
                                               !            3: natural lake with dam
   real(r8), allocatable :: dam_elv       (:)  ! dam elevation [m]

   integer :: numresv
   integer,  allocatable :: bsn2resv      (:)

   integer :: numresv_uniq
   integer,  allocatable :: resv_hylak_id (:)  ! HydroLAKE ID of reservoir
   integer,  allocatable :: resv_loc2glb  (:)  ! global index of reservoir

   ! parameters
   real(r8), allocatable :: volresv_total (:)  ! total reservoir volume      [m^3]
   real(r8), allocatable :: volresv_emerg (:)  ! emergency reservoir volume  [m^3]
   real(r8), allocatable :: volresv_adjust(:)  ! adjustment reservoir volume [m^3]
   real(r8), allocatable :: volresv_normal(:)  ! normal reservoir volume     [m^3]

   real(r8), allocatable :: qresv_mean    (:)  ! mean natural reservoir outflow [m^3/s]
   real(r8), allocatable :: qresv_flood   (:)  ! flood reservoir outflow        [m^3/s]
   real(r8), allocatable :: qresv_adjust  (:)  ! adjustment reservoir outflow   [m^3/s]
   real(r8), allocatable :: qresv_normal  (:)  ! normal reservoir outflow       [m^3/s]

   ! time variables
   real(r8), allocatable :: volresv       (:)  ! reservoir water volume  [m^3]
   real(r8), allocatable :: qresv_in      (:)  ! reservoir inflow        [m^3/s]
   real(r8), allocatable :: qresv_out     (:)  ! reservoir outflow       [m^3/s]

   ! time average variables for output
   real(r8), allocatable :: volresv_ta    (:)  ! reservoir water volume [m^3]
   real(r8), allocatable :: qresv_in_ta   (:)  ! inflow to reservoir    [m^3/s]
   real(r8), allocatable :: qresv_out_ta  (:)  ! outflow from reservoir [m^3/s]

   ! -- PUBLIC SUBROUTINEs --
   PUBLIC :: reservoir_init
   PUBLIC :: reservoir_operation
   PUBLIC :: reservoir_gather_var
   PUBLIC :: reservoir_final

CONTAINS

   ! -------
   SUBROUTINE reservoir_init ( )

   USE MOD_Vars_Global, only: spval
   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_NetCDFSerial
   USE MOD_Utils
   USE MOD_Catch_BasinNetwork,     only : numbasin, basinindex
   USE MOD_Catch_RiverLakeNetwork, only : lake_id,  wtsrfelv, bedelv, lakeinfo

   IMPLICIT NONE

   ! Local variables
   character(len=256) :: lake_info_file, basin_info_file, resv_info_file
   integer :: funit, maxlakeid, this_lak_id, this_lak_typ, numlake, ierr
   integer :: mesg(2), isrc, idest, iworker

   integer,  allocatable :: all_lak_typ(:), senddata(:), recvdata(:)

   integer,  allocatable :: lake_id_basin(:), ilat_outlet_basin(:), ilon_outlet_basin(:)
   integer,  allocatable :: lake_id_resv (:), ilat_outlet_resv (:), ilon_outlet_resv (:)

   real(r8), allocatable :: all_qmean_basin(:), all_qflood_basin(:)
   real(r8), allocatable :: all_qresv_mean (:), all_qresv_flood (:)

   integer,  allocatable :: resv_bsn_id(:), order(:)
   real(r8), allocatable :: rcache(:)

   integer :: nbasin, ibasin, irsv, nrecv



      IF (p_is_master) THEN

         lake_info_file = trim(DEF_dir_runtime)//'/reservoir/HydroLAKES_Reservoir.txt'

         funit = 101
         open (funit, file = trim(lake_info_file), action = 'READ')

         maxlakeid = 0
         read(funit,*,iostat=ierr) ! skip first header line
         DO WHILE (.true.)
            read(funit,*,iostat=ierr) this_lak_id
            IF (ierr /= 0) exit
            maxlakeid = max(maxlakeid, this_lak_id)
         ENDDO

         allocate (all_lak_typ (maxlakeid))
         all_lak_typ(:) = 0

         rewind (funit)
         read(funit,*,iostat=ierr)
         DO WHILE (.true.)
            read(funit,*,iostat=ierr) this_lak_id, this_lak_typ
            IF (ierr /= 0) exit
            all_lak_typ(this_lak_id) = this_lak_typ
         ENDDO

         close (funit)

      ENDIF

      IF (p_is_worker) THEN
         allocate (lake_type (numbasin));   lake_type(:) = 0
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_worker) THEN

         IF (numbasin > 0) THEN
            numlake = count(lake_id > 0)
         ELSE
            numlake = 0
         ENDIF

         mesg = (/p_iam_glb, numlake/)
         CALL mpi_send (mesg, 2, MPI_INTEGER, p_address_master, &
            mpi_tag_mesg, p_comm_glb, p_err)

         IF (numlake > 0) THEN

            allocate (senddata (numlake))
            senddata = pack(lake_id, lake_id > 0)
            CALL mpi_send (senddata, numlake, MPI_INTEGER, &
               p_address_master, mpi_tag_data, p_comm_glb, p_err)

            allocate (recvdata (numlake))
            CALL mpi_recv (recvdata, numlake, MPI_INTEGER, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)
            CALL unpack_inplace (recvdata, lake_id > 0, lake_type)

            deallocate (senddata)
            deallocate (recvdata)
         ENDIF

      ENDIF

      IF (p_is_master) THEN
         DO iworker = 0, p_np_worker-1

            CALL mpi_recv (mesg, 2, MPI_INTEGER, MPI_ANY_SOURCE, &
               mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc    = mesg(1)
            numlake = mesg(2)
            IF (numlake > 0) THEN

               allocate(recvdata (numlake))
               CALL mpi_recv (recvdata, numlake, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               allocate(senddata (numlake))

               senddata = all_lak_typ(recvdata)
               CALL mpi_send (senddata, numlake, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_err)

               deallocate (recvdata)
               deallocate (senddata)
            ENDIF

         ENDDO
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else

      DO ibasin = 1, numbasin
         IF (lake_id(ibasin) > 0) THEN
            lake_type(ibasin) = all_lak_typ(lake_id(ibasin))
         ENDIF
      ENDDO

#endif


      IF (p_is_worker) THEN

         allocate (dam_elv (numbasin))
         dam_elv(:) = spval

         DO ibasin = 1, numbasin
            IF ((lake_type(ibasin) == 2) .or. (lake_type(ibasin) == 3)) THEN
               dam_elv(ibasin) = wtsrfelv(ibasin)
            ENDIF
         ENDDO
      ENDIF



      IF (DEF_Reservoir_Method == 1) THEN

         IF (p_is_worker) THEN
            IF (numbasin > 0) THEN

               numresv = count(lake_type == 2)

               IF (numresv > 0) THEN
                  allocate (volresv_total (numresv))
                  allocate (volresv_emerg (numresv))
                  allocate (volresv_adjust(numresv))
                  allocate (volresv_normal(numresv))

                  allocate (qresv_mean    (numresv))
                  allocate (qresv_flood   (numresv))
                  allocate (qresv_adjust  (numresv))
                  allocate (qresv_normal  (numresv))

                  allocate (volresv       (numresv))
                  allocate (qresv_in      (numresv))
                  allocate (qresv_out     (numresv))

                  allocate (volresv_ta    (numresv))
                  allocate (qresv_in_ta   (numresv))
                  allocate (qresv_out_ta  (numresv))
               ENDIF

            ELSE
               numresv = 0
            ENDIF
         ENDIF

         IF (p_is_worker) THEN
            IF (numresv > 0) THEN

               allocate (resv_bsn_id (numresv))
               resv_bsn_id = pack(basinindex, lake_type == 2)

               allocate (bsn2resv (numbasin))

               irsv = 0
               DO ibasin = 1, numbasin
                  IF (lake_type(ibasin) == 2) THEN
                     irsv = irsv + 1
                     bsn2resv (ibasin) = irsv
                     volresv_total (irsv) = lakeinfo(ibasin)%volume( dam_elv(ibasin)-bedelv(ibasin) )
                     volresv_emerg (irsv) = volresv_total(irsv) * 0.94
                     volresv_adjust(irsv) = volresv_total(irsv) * 0.77
                     volresv_normal(irsv) = volresv_total(irsv) * 0.7
                  ENDIF
               ENDDO

            ENDIF
         ENDIF

         ! read in parameters from file.
         IF (p_is_master) THEN

            basin_info_file = DEF_CatchmentMesh_data
            CALL ncio_read_serial (basin_info_file, 'lake_id',     lake_id_basin    )
            CALL ncio_read_serial (basin_info_file, 'ilat_outlet', ilat_outlet_basin)
            CALL ncio_read_serial (basin_info_file, 'ilon_outlet', ilon_outlet_basin)

            resv_info_file = trim(DEF_dir_runtime)//'/reservoir/qmean_qflood_of_reservoir.nc'
            CALL ncio_read_serial (resv_info_file, 'hylak_id',    lake_id_resv    )
            CALL ncio_read_serial (resv_info_file, 'ilat_outlet', ilat_outlet_resv)
            CALL ncio_read_serial (resv_info_file, 'ilon_outlet', ilon_outlet_resv)
            CALL ncio_read_serial (resv_info_file, 'qmean',  all_qresv_mean )
            CALL ncio_read_serial (resv_info_file, 'qflood', all_qresv_flood)

            nbasin = size(lake_id_basin)

            allocate(all_qmean_basin  (nbasin));  all_qmean_basin (:) = spval
            allocate(all_qflood_basin (nbasin));  all_qflood_basin(:) = spval

            DO ibasin = 1, nbasin
               IF (all_lak_typ(lake_id_basin(ibasin)) == 2) THEN
                  DO irsv = 1, size(lake_id_resv)
                     IF (lake_id_basin(ibasin) == lake_id_resv(irsv)) THEN
                        IF (ilat_outlet_basin(ibasin) == ilat_outlet_resv(irsv)) THEN
                           IF (ilon_outlet_basin(ibasin) == ilon_outlet_resv(irsv)) THEN
                              all_qmean_basin (ibasin) = all_qresv_mean (irsv)
                              all_qflood_basin(ibasin) = all_qresv_flood(irsv)
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO

         ENDIF

#ifdef USEMPI
         IF (p_is_master) THEN
            DO iworker = 0, p_np_worker-1

               CALL mpi_recv (mesg(1:2), 2, MPI_INTEGER, &
                  MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
               isrc  = mesg(1)
               nrecv = mesg(2)

               IF (nrecv > 0) THEN

                  allocate (resv_bsn_id (nrecv))
                  allocate (rcache (nrecv))

                  CALL mpi_recv (resv_bsn_id, nrecv, MPI_INTEGER, isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

                  idest = isrc

                  rcache = all_qmean_basin(resv_bsn_id)
                  CALL mpi_send (rcache, nrecv, MPI_REAL8, idest, mpi_tag_data, p_comm_glb, p_err)

                  rcache = all_qflood_basin(resv_bsn_id)
                  CALL mpi_send (rcache, nrecv, MPI_REAL8, idest, mpi_tag_data, p_comm_glb, p_err)

                  deallocate (resv_bsn_id)
                  deallocate (rcache)

               ENDIF
            ENDDO
         ENDIF

         IF (p_is_worker) THEN

            mesg(1:2) = (/p_iam_glb, numresv/)
            CALL mpi_send (mesg(1:2), 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)

            IF (numresv > 0) THEN
               CALL mpi_send (resv_bsn_id, numresv, MPI_INTEGER, &
                  p_address_master, mpi_tag_data, p_comm_glb, p_err)

               CALL mpi_recv (qresv_mean, numresv, MPI_REAL8, &
                  p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

               CALL mpi_recv (qresv_flood, numresv, MPI_REAL8, &
                  p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)
            ENDIF
         ENDIF
#else
         IF (numresv > 0) THEN
            qresv_mean  = all_qmean_basin (resv_bsn_id)
            qresv_flood = all_qflood_basin(resv_bsn_id)
         ENDIF
#endif

         IF (p_is_worker) THEN
            IF (numresv > 0) THEN
               qresv_normal = volresv_normal*0.7/(180*86400) + qresv_mean*0.25
               qresv_adjust = (qresv_normal + qresv_flood) * 0.5
            ENDIF
         ENDIF


         IF (p_is_master) THEN

            allocate (order (nbasin)); order = (/(ibasin,ibasin=1,nbasin)/)

            CALL quicksort (nbasin, lake_id_basin, order)

            deallocate (order)

            numresv_uniq = 0
            DO ibasin = 1, nbasin
               IF (lake_id_basin(ibasin) > 0) THEN
                  IF (all_lak_typ(lake_id_basin(ibasin)) == 2) THEN
                     IF (ibasin == 1) THEN
                        numresv_uniq = 1
                     ELSE
                        IF (lake_id_basin(ibasin) /= lake_id_basin(ibasin-1)) THEN
                           numresv_uniq = numresv_uniq + 1
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO

            IF (numresv_uniq > 0) THEN
               allocate (resv_hylak_id (numresv_uniq))
               numresv_uniq = 0
               DO ibasin = 1, nbasin
                  IF (lake_id_basin(ibasin) > 0) THEN
                     IF (all_lak_typ(lake_id_basin(ibasin)) == 2) THEN
                        IF (ibasin == 1) THEN
                           numresv_uniq = 1
                           resv_hylak_id(numresv_uniq) = lake_id_basin(ibasin)
                        ELSE
                           IF (lake_id_basin(ibasin) /= lake_id_basin(ibasin-1)) THEN
                              numresv_uniq = numresv_uniq + 1
                              resv_hylak_id(numresv_uniq) = lake_id_basin(ibasin)
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF

         ENDIF

#ifdef USEMPI
         CALL mpi_bcast (numresv_uniq, 1, MPI_INTEGER, p_address_master, p_comm_glb, p_err)

         IF (numresv_uniq > 0) THEN
            IF (.not. p_is_master)  allocate (resv_hylak_id (numresv_uniq))
            CALL mpi_bcast (resv_hylak_id, numresv_uniq, MPI_INTEGER, p_address_master, p_comm_glb, p_err)
         ENDIF
#endif

         IF (p_is_worker) THEN
            IF (numresv > 0) THEN
               allocate (resv_loc2glb (numresv))
               DO ibasin = 1, numbasin
                  IF (lake_type(ibasin) == 2) THEN
                     resv_loc2glb(bsn2resv(ibasin)) = &
                        find_in_sorted_list1 (lake_id(ibasin), numresv_uniq, resv_hylak_id)
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

      ENDIF

      IF (allocated(all_lak_typ      )) deallocate (all_lak_typ      )
      IF (allocated(lake_id_basin    )) deallocate (lake_id_basin    )
      IF (allocated(ilat_outlet_basin)) deallocate (ilat_outlet_basin)
      IF (allocated(ilon_outlet_basin)) deallocate (ilon_outlet_basin)
      IF (allocated(lake_id_resv     )) deallocate (lake_id_resv     )
      IF (allocated(ilat_outlet_resv )) deallocate (ilat_outlet_resv )
      IF (allocated(ilon_outlet_resv )) deallocate (ilon_outlet_resv )
      IF (allocated(all_qmean_basin  )) deallocate (all_qmean_basin  )
      IF (allocated(all_qflood_basin )) deallocate (all_qflood_basin )
      IF (allocated(all_qresv_mean   )) deallocate (all_qresv_mean   )
      IF (allocated(all_qresv_flood  )) deallocate (all_qresv_flood  )
      IF (allocated(resv_bsn_id      )) deallocate (resv_bsn_id      )

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
   IMPLICIT NONE

   real(r8), intent(in)  :: varin  (:)
   real(r8), intent(out) :: varout (:)

   ! local variables
   integer :: irsv

      IF (numresv_uniq == 0) RETURN

      IF (p_is_worker) THEN
         varout(:) = 0.
         DO irsv = 1, numresv
            varout(resv_loc2glb(irsv)) = varout(resv_loc2glb(irsv)) + varin(irsv)
         ENDDO
#ifdef USEMPI
         CALL mpi_reduce (MPI_IN_PLACE, varout, numresv_uniq, MPI_REAL8, &
            MPI_SUM, p_root, p_comm_worker, p_err)
#endif
      ENDIF

#ifdef USEMPI
      IF (p_iam_worker == p_root) THEN
         CALL mpi_send (varout, numresv_uniq, MPI_REAL8, p_address_master, &
            mpi_tag_data, p_comm_glb, p_err)
      ENDIF
      IF (p_is_master) THEN
         CALL mpi_recv (varout, numresv_uniq, MPI_REAL8, p_address_worker(p_root), &
            mpi_tag_data, p_comm_glb, p_stat, p_err)
      ENDIF
#endif

   END SUBROUTINE reservoir_gather_var


   SUBROUTINE reservoir_final ()

   IMPLICIT NONE

      IF (allocated(lake_type     )) deallocate (lake_type     )
      IF (allocated(dam_elv       )) deallocate (dam_elv       )

      IF (allocated(bsn2resv      )) deallocate (bsn2resv      )
      IF (allocated(resv_hylak_id )) deallocate (resv_hylak_id )

      IF (allocated(volresv_total )) deallocate (volresv_total )
      IF (allocated(volresv_emerg )) deallocate (volresv_emerg )
      IF (allocated(volresv_adjust)) deallocate (volresv_adjust)
      IF (allocated(volresv_normal)) deallocate (volresv_normal)

      IF (allocated(qresv_mean    )) deallocate (qresv_mean    )
      IF (allocated(qresv_flood   )) deallocate (qresv_flood   )
      IF (allocated(qresv_adjust  )) deallocate (qresv_adjust  )
      IF (allocated(qresv_normal  )) deallocate (qresv_normal  )

      IF (allocated(volresv       )) deallocate (volresv       )
      IF (allocated(qresv_in      )) deallocate (qresv_in      )
      IF (allocated(qresv_out     )) deallocate (qresv_out     )

      IF (allocated(volresv_ta    )) deallocate (volresv_ta    )
      IF (allocated(qresv_in_ta   )) deallocate (qresv_in_ta   )
      IF (allocated(qresv_out_ta  )) deallocate (qresv_out_ta  )

   END SUBROUTINE reservoir_final

END MODULE MOD_Catch_Reservoir
#endif
