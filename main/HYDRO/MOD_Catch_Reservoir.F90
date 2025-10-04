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
   USE MOD_Catch_BasinNetwork, only : numbasin, basinindex, numresv, bsn2resv, lake_id, lake_type

   real(r8), allocatable :: dam_elv       (:)  ! dam elevation [m]

   ! parameters
   real(r8), allocatable :: volresv_total (:)  ! total reservoir volume      [m^3]
   real(r8), allocatable :: volresv_emerg (:)  ! emergency reservoir volume  [m^3]
   real(r8), allocatable :: volresv_adjust(:)  ! adjustment reservoir volume [m^3]
   real(r8), allocatable :: volresv_normal(:)  ! normal reservoir volume     [m^3]

   real(r8), allocatable :: qresv_mean    (:)  ! mean natural reservoir outflow [m^3/s]
   real(r8), allocatable :: qresv_flood   (:)  ! flood reservoir outflow        [m^3/s]
   real(r8), allocatable :: qresv_adjust  (:)  ! adjustment reservoir outflow   [m^3/s]
   real(r8), allocatable :: qresv_normal  (:)  ! normal reservoir outflow       [m^3/s]

   integer,  allocatable :: dam_build_year(:)  ! year in which the dam/barrier was built

   ! time variables
   real(r8), allocatable :: volresv       (:)  ! reservoir water volume  [m^3]
   real(r8), allocatable :: qresv_in      (:)  ! reservoir inflow        [m^3/s]
   real(r8), allocatable :: qresv_out     (:)  ! reservoir outflow       [m^3/s]

   ! time average variables for output
   real(r8), allocatable :: volresv_ta    (:)  ! reservoir water volume [m^3]
   real(r8), allocatable :: qresv_in_ta   (:)  ! inflow to reservoir    [m^3/s]
   real(r8), allocatable :: qresv_out_ta  (:)  ! outflow from reservoir [m^3/s]

   integer :: numresv_uniq
   integer,  allocatable :: resv_hylak_id (:)  ! HydroLAKE ID of reservoir
   integer,  allocatable :: resv_loc2glb  (:)  ! global index of reservoir

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
   USE MOD_Catch_RiverLakeNetwork, only : wtsrfelv, bedelv, lakeinfo

   IMPLICIT NONE

   ! Local variables
   character(len=256) :: basin_info_file, resv_info_file
   integer :: maxlakeid, mesg(2), isrc, idest, iworker

   integer,  allocatable :: lake_id_basin (:), ilat_outlet_basin(:), ilon_outlet_basin(:)
   integer,  allocatable :: lake_id_resv  (:), ilat_outlet_resv (:), ilon_outlet_resv (:)
   integer,  allocatable :: lake_type_resv(:), lake_id2typ      (:)
   integer,  allocatable :: resv_bsn_id   (:), idlist           (:)
   integer,  allocatable :: all_year_basin(:), all_year_resv    (:)

   real(r8), allocatable :: all_vol_basin (:), all_qmean_basin  (:), all_qflood_basin (:)
   real(r8), allocatable :: all_vol_resv  (:), all_qmean_resv   (:), all_qflood_resv  (:)
   real(r8), allocatable :: all_dhgt_basin(:), all_dhgt_resv    (:), dam_height       (:)

   real(r8), allocatable :: rcache(:)
   integer,  allocatable :: icache(:)
   integer :: nbasin, ibasin, irsv, nrecv, nrsv, iloc


      IF (p_is_worker) THEN

         IF (numresv > 0) THEN

            allocate (dam_elv       (numresv))

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

            allocate (dam_build_year(numresv))

            allocate (dam_height    (numresv))

         ENDIF
      ENDIF

      ! read in parameters from file.
      IF (p_is_master) THEN

         basin_info_file = DEF_CatchmentMesh_data
         CALL ncio_read_serial (basin_info_file, 'lake_id',     lake_id_basin    )
         CALL ncio_read_serial (basin_info_file, 'ilat_outlet', ilat_outlet_basin)
         CALL ncio_read_serial (basin_info_file, 'ilon_outlet', ilon_outlet_basin)

         maxlakeid = maxval(lake_id_basin)
         IF (maxlakeid > 0) THEN

            resv_info_file = trim(DEF_dir_runtime)//'/HydroLAKES_Reservoir.nc'
            CALL ncio_read_serial (resv_info_file, 'hylak_id',    lake_id_resv    )
            CALL ncio_read_serial (resv_info_file, 'lake_type',   lake_type_resv  )
            CALL ncio_read_serial (resv_info_file, 'ilat_outlet', ilat_outlet_resv)
            CALL ncio_read_serial (resv_info_file, 'ilon_outlet', ilon_outlet_resv)
            CALL ncio_read_serial (resv_info_file, 'volresv', all_vol_resv   )
            CALL ncio_read_serial (resv_info_file, 'qmean',   all_qmean_resv )
            CALL ncio_read_serial (resv_info_file, 'qflood',  all_qflood_resv)

            CALL ncio_read_serial (resv_info_file, 'build_year', all_year_resv)
            CALL ncio_read_serial (resv_info_file, 'dam_height', all_dhgt_resv)

            allocate (lake_id2typ (-1:maxlakeid))
            lake_id2typ(:) = 0

            DO irsv = 1, size(lake_id_resv)
               IF (lake_id_resv(irsv) <= maxlakeid) THEN
                  lake_id2typ(lake_id_resv(irsv)) = lake_type_resv(irsv)
               ENDIF
            ENDDO


            nbasin = size(lake_id_basin)

            allocate(all_vol_basin    (nbasin));  all_vol_basin   (:) = spval
            allocate(all_qmean_basin  (nbasin));  all_qmean_basin (:) = spval
            allocate(all_qflood_basin (nbasin));  all_qflood_basin(:) = spval
            allocate(all_year_basin   (nbasin));  all_year_basin  (:) = -99
            allocate(all_dhgt_basin   (nbasin));  all_dhgt_basin  (:) = spval

            DO ibasin = 1, nbasin
               IF (lake_id2typ(lake_id_basin(ibasin)) >= 2) THEN
                  DO irsv = 1, size(lake_id_resv)
                     IF (lake_id_basin(ibasin) == lake_id_resv(irsv)) THEN
                        IF (ilat_outlet_basin(ibasin) == ilat_outlet_resv(irsv)) THEN
                           IF (ilon_outlet_basin(ibasin) == ilon_outlet_resv(irsv)) THEN
                              all_vol_basin   (ibasin) = all_vol_resv   (irsv)
                              all_qmean_basin (ibasin) = all_qmean_resv (irsv)
                              all_qflood_basin(ibasin) = all_qflood_resv(irsv)
                              all_year_basin  (ibasin) = all_year_resv  (irsv)
                              all_dhgt_basin  (ibasin) = all_dhgt_resv  (irsv)
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO

         ENDIF

      ENDIF

      IF (p_is_worker) THEN
         IF (numresv > 0) THEN
            allocate (resv_bsn_id (numresv))
            resv_bsn_id = pack(basinindex, lake_type >= 2)
         ENDIF
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
               allocate (icache (nrecv))

               CALL mpi_recv (resv_bsn_id, nrecv, MPI_INTEGER, isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               idest = isrc

               rcache = all_vol_basin(resv_bsn_id)
               CALL mpi_send (rcache, nrecv, MPI_REAL8, idest, mpi_tag_data, p_comm_glb, p_err)

               rcache = all_qmean_basin(resv_bsn_id)
               CALL mpi_send (rcache, nrecv, MPI_REAL8, idest, mpi_tag_data, p_comm_glb, p_err)

               rcache = all_qflood_basin(resv_bsn_id)
               CALL mpi_send (rcache, nrecv, MPI_REAL8, idest, mpi_tag_data, p_comm_glb, p_err)

               icache = all_year_basin(resv_bsn_id)
               CALL mpi_send (icache, nrecv, MPI_INTEGER, idest, mpi_tag_data, p_comm_glb, p_err)

               rcache = all_dhgt_basin(resv_bsn_id)
               CALL mpi_send (rcache, nrecv, MPI_REAL8, idest, mpi_tag_data, p_comm_glb, p_err)

               deallocate (resv_bsn_id)
               deallocate (rcache)
               deallocate (icache)

            ENDIF
         ENDDO
      ENDIF

      IF (p_is_worker) THEN

         mesg(1:2) = (/p_iam_glb, numresv/)
         CALL mpi_send (mesg(1:2), 2, MPI_INTEGER, p_address_master, mpi_tag_mesg, p_comm_glb, p_err)

         IF (numresv > 0) THEN
            CALL mpi_send (resv_bsn_id, numresv, MPI_INTEGER, &
               p_address_master, mpi_tag_data, p_comm_glb, p_err)

            CALL mpi_recv (volresv_total, numresv, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            CALL mpi_recv (qresv_mean, numresv, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            CALL mpi_recv (qresv_flood, numresv, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            CALL mpi_recv (dam_build_year, numresv, MPI_INTEGER, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)

            CALL mpi_recv (dam_height, numresv, MPI_REAL8, &
               p_address_master, mpi_tag_data, p_comm_glb, p_stat, p_err)
         ENDIF
      ENDIF
#else
      IF (numresv > 0) THEN
         volresv_total  = all_vol_basin   (resv_bsn_id)
         qresv_mean     = all_qmean_basin (resv_bsn_id)
         qresv_flood    = all_qflood_basin(resv_bsn_id)
         dam_build_year = all_year_basin  (resv_bsn_id)
         dam_height     = all_dhgt_basin  (resv_bsn_id)
      ENDIF
#endif

      IF (p_is_worker) THEN
         DO ibasin = 1, numbasin
            IF (lake_type(ibasin) >= 2) THEN

               irsv = bsn2resv(ibasin)

               dam_height(irsv) = max(dam_height(irsv), wtsrfelv(ibasin)-bedelv(ibasin))
               dam_height(irsv) = max(dam_height(irsv), lakeinfo(ibasin)%surface(volresv_total(irsv)))
               dam_height(irsv) = min(dam_height(irsv), 335.)

               dam_elv       (irsv) = bedelv(ibasin) + dam_height(irsv)
               volresv_total (irsv) = lakeinfo(ibasin)%volume( dam_height(irsv) )

               volresv_emerg (irsv) = volresv_total(irsv) * 0.94
               volresv_adjust(irsv) = volresv_total(irsv) * 0.77
               volresv_normal(irsv) = volresv_total(irsv) * 0.7

               qresv_normal  (irsv) = volresv_normal(irsv)*0.7/(180*86400) + qresv_mean(irsv)*0.25
               qresv_adjust  (irsv) = (qresv_normal(irsv) + qresv_flood(irsv)) * 0.5
            ENDIF
         ENDDO
      ENDIF


      IF (p_is_master) THEN

         nrsv = count(lake_id2typ(lake_id_basin) >= 2)

         IF (nrsv > 0) THEN

            allocate (idlist (nrsv))

            numresv_uniq = 0
            DO ibasin = 1, size(lake_id_basin)
               IF (lake_id2typ(lake_id_basin(ibasin)) >= 2) THEN
                  CALL insert_into_sorted_list1 ( &
                     lake_id_basin(ibasin), numresv_uniq, idlist, iloc)
               ENDIF
            ENDDO

            allocate (resv_hylak_id (numresv_uniq))
            resv_hylak_id = idlist(1:numresv_uniq)

            deallocate (idlist)

         ELSE
            numresv_uniq = 0
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
               IF (lake_type(ibasin) >= 2) THEN
                  resv_loc2glb(bsn2resv(ibasin)) = &
                     find_in_sorted_list1 (lake_id(ibasin), numresv_uniq, resv_hylak_id)
               ENDIF
            ENDDO
         ENDIF
      ENDIF


      IF (allocated(lake_id2typ      )) deallocate (lake_id2typ      )
      IF (allocated(lake_id_basin    )) deallocate (lake_id_basin    )
      IF (allocated(ilat_outlet_basin)) deallocate (ilat_outlet_basin)
      IF (allocated(ilon_outlet_basin)) deallocate (ilon_outlet_basin)
      IF (allocated(lake_id_resv     )) deallocate (lake_id_resv     )
      IF (allocated(lake_type_resv   )) deallocate (lake_type_resv   )
      IF (allocated(ilat_outlet_resv )) deallocate (ilat_outlet_resv )
      IF (allocated(ilon_outlet_resv )) deallocate (ilon_outlet_resv )
      IF (allocated(all_qmean_basin  )) deallocate (all_qmean_basin  )
      IF (allocated(all_qflood_basin )) deallocate (all_qflood_basin )
      IF (allocated(all_qmean_resv   )) deallocate (all_qmean_resv   )
      IF (allocated(all_qflood_resv  )) deallocate (all_qflood_resv  )
      IF (allocated(all_year_basin   )) deallocate (all_year_basin   )
      IF (allocated(all_year_resv    )) deallocate (all_year_resv    )
      IF (allocated(all_dhgt_basin   )) deallocate (all_dhgt_basin   )
      IF (allocated(all_dhgt_resv    )) deallocate (all_dhgt_resv    )
      IF (allocated(dam_height       )) deallocate (dam_height       )
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
   integer :: irsv, iworker
   real(r8), allocatable :: varall (:,:)

      IF (numresv_uniq == 0) RETURN

      IF (p_is_worker) THEN
         varout(:) = spval
         DO irsv = 1, numresv
            IF (varin(irsv) /= spval) THEN
               IF (varout(resv_loc2glb(irsv)) /= spval) THEN
                  varout(resv_loc2glb(irsv)) = varout(resv_loc2glb(irsv)) + varin(irsv)
               ELSE
                  varout(resv_loc2glb(irsv)) = varin(irsv)
               ENDIF
            ENDIF
         ENDDO
#ifdef USEMPI
         IF (p_iam_worker == p_root) THEN
            allocate (varall (numresv_uniq,0:p_np_worker-1))
         ENDIF

         CALL mpi_gather (varout, numresv_uniq, MPI_REAL8, &
            varall, numresv_uniq, MPI_REAL8, p_root, p_comm_worker, p_err)

         IF (p_iam_worker == p_root) THEN

            DO irsv = 1, numresv_uniq
               DO iworker = 0, p_np_worker-1
                  IF (iworker /= p_root) THEN
                     IF (varall(irsv,iworker) /= spval) THEN
                        IF (varout(irsv) /= spval) THEN
                           varout(irsv) = varout(irsv) + varall(irsv,iworker)
                        ELSE
                           varout(irsv) = varall(irsv,iworker)
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO

            deallocate (varall)
         ENDIF
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

      IF (allocated(dam_elv       )) deallocate (dam_elv       )

      IF (allocated(resv_hylak_id )) deallocate (resv_hylak_id )
      IF (allocated(resv_loc2glb  )) deallocate (resv_loc2glb  )

      IF (allocated(volresv_total )) deallocate (volresv_total )
      IF (allocated(volresv_emerg )) deallocate (volresv_emerg )
      IF (allocated(volresv_adjust)) deallocate (volresv_adjust)
      IF (allocated(volresv_normal)) deallocate (volresv_normal)

      IF (allocated(qresv_mean    )) deallocate (qresv_mean    )
      IF (allocated(qresv_flood   )) deallocate (qresv_flood   )
      IF (allocated(qresv_adjust  )) deallocate (qresv_adjust  )
      IF (allocated(qresv_normal  )) deallocate (qresv_normal  )

      IF (allocated(dam_build_year)) deallocate (dam_build_year)

      IF (allocated(volresv       )) deallocate (volresv       )
      IF (allocated(qresv_in      )) deallocate (qresv_in      )
      IF (allocated(qresv_out     )) deallocate (qresv_out     )

      IF (allocated(volresv_ta    )) deallocate (volresv_ta    )
      IF (allocated(qresv_in_ta   )) deallocate (qresv_in_ta   )
      IF (allocated(qresv_out_ta  )) deallocate (qresv_out_ta  )

   END SUBROUTINE reservoir_final

END MODULE MOD_Catch_Reservoir
#endif
