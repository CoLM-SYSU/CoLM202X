#include <define.h>

#ifdef GridRiverLakeFlow
MODULE MOD_Grid_RiverLakeHist
!--------------------------------------------------------------------------------
! DESCRIPTION:
!
!     Write out model results in lateral hydrological processes to history files.
!
! Created by Shupeng Zhang, May 2023
!--------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_DataType

   ! -- ACC Fluxes --
   real(r8), allocatable :: acctime       (:)

   real(r8), allocatable :: a_wdsrf_ucat  (:)
   real(r8), allocatable :: a_veloc_riv   (:)
   real(r8), allocatable :: a_discharge   (:)
   real(r8), allocatable :: a_floodarea   (:)  ! flooded area [m^2]

   ! for reservoirs
   real(r8), allocatable :: acctime_resv  (:)

   real(r8), allocatable :: a_volresv     (:)  ! reservoir water volume [m^3]
   real(r8), allocatable :: a_qresv_in    (:)  ! inflow to reservoir    [m^3/s]
   real(r8), allocatable :: a_qresv_out   (:)  ! outflow from reservoir [m^3/s]

   ! auxiliary data
   type(block_data_real8_2d) :: sumarea_ucat,   sumarea_inpmat
   logical, allocatable      :: filter_ucat(:), filter_inpmat(:)

   ! -- PUBLIC SUBROUTINEs --
   PUBLIC :: hist_grid_riverlake_init
   PUBLIC :: hist_grid_riverlake_out
   PUBLIC :: hist_grid_riverlake_final

!--------------------------------------------------------------------------
CONTAINS

   !---------------------------------------
   SUBROUTINE hist_grid_riverlake_init ()

   USE MOD_WorkerPushData
   USE MOD_HistGridded
   USE MOD_Grid_RiverLakeNetwork
   USE MOD_Grid_Reservoir, only: numresv
   USE MOD_Mesh,           only: numelm
   USE MOD_LandPatch,      only: numpatch, elm_patch
   IMPLICIT NONE

   ! Local Variables
   real(r8), allocatable :: vec_ucat(:), vec_elm(:)
   integer :: ielm, istt, iend

      IF (p_is_worker) THEN

         IF (numucat > 0) THEN
            allocate (acctime      (numucat))
            allocate (a_wdsrf_ucat (numucat))
            allocate (a_veloc_riv  (numucat))
            allocate (a_discharge  (numucat))
            allocate (a_floodarea  (numucat))
         ENDIF

         IF (numresv > 0) THEN
            allocate (acctime_resv (numresv))
            allocate (a_volresv    (numresv))
            allocate (a_qresv_in   (numresv))
            allocate (a_qresv_out  (numresv))
         ENDIF

      ENDIF

      IF (p_is_io) THEN
         CALL allocate_block_data (ghist, sumarea_ucat  )
         CALL allocate_block_data (ghist, sumarea_inpmat)
      ENDIF
      IF (p_is_worker) THEN
         IF (numpatch > 0) allocate (filter_ucat   (numpatch))
         IF (numpatch > 0) allocate (filter_inpmat (numpatch))
      ENDIF

      IF (p_is_worker) THEN
         IF (numucat  > 0) allocate (vec_ucat (numucat ))
         IF (numelm   > 0) allocate (vec_elm  (numelm  ))
         IF (numucat  > 0) vec_ucat = 1.
         IF (numelm   > 0) vec_elm  = 0.

         CALL worker_push_data (push_ucat2elm, vec_ucat, vec_elm, fillvalue = 0.)

         DO ielm = 1, numelm
            istt = elm_patch%substt(ielm)
            iend = elm_patch%subend(ielm)
            filter_ucat(istt:iend) = vec_elm(ielm) > 0.
         ENDDO
      ENDIF

      CALL mp2g_hist%get_sumarea (sumarea_ucat, filter_ucat)

      IF (p_is_worker) THEN
         DO ielm = 1, numelm
            istt = elm_patch%substt(ielm)
            iend = elm_patch%subend(ielm)
            filter_inpmat(istt:iend) = sum(inpmat_area_u2e(:,ielm)) > 0.
         ENDDO
      ENDIF

      CALL mp2g_hist%get_sumarea (sumarea_inpmat, filter_inpmat)

      IF (allocated (vec_ucat)) deallocate (vec_ucat)
      IF (allocated (vec_elm )) deallocate (vec_elm )

   END SUBROUTINE hist_grid_riverlake_init

   !---------------------------------------
   SUBROUTINE hist_grid_riverlake_out (file_hist, idate, itime_in_file)

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_NetCDFSerial
   USE MOD_WorkerPushData
   USE MOD_Grid_RiverLakeNetwork
   USE MOD_Grid_Reservoir
   USE MOD_SpatialMapping
   USE MOD_HistGridded, only: flux_map_and_write_2d
   USE MOD_LandPatch,   only: numpatch, elm_patch
   USE MOD_Mesh,        only: numelm
   IMPLICIT NONE

   character(len=*), intent(in) :: file_hist
   integer, intent(in) :: idate(3)
   integer, intent(in) :: itime_in_file

   ! Local variables
   character(len=256) :: file_hist_resv
   logical :: fexists
   integer :: itime_in_file_resv, ielm, istt, iend, i

   real(r8), allocatable :: acc_vec     (:)
   real(r8), allocatable :: acc_vec_pch (:)
   real(r8), allocatable :: a_floodfrc        (:)    ! flooded area fraction
   real(r8), allocatable :: a_floodfrc_inpmat (:,:)  ! flooded area fraction


      IF (p_is_worker) THEN
         IF (numucat > 0) THEN
            WHERE (acctime > 0)
               a_wdsrf_ucat = a_wdsrf_ucat / acctime
               a_veloc_riv  = a_veloc_riv  / acctime
               a_discharge  = a_discharge  / acctime
               a_floodarea  = a_floodarea  / acctime
            END WHERE
         ENDIF
      ENDIF

      IF (p_is_worker) THEN
         IF (numelm   > 0) allocate (acc_vec     (numelm  ))
         IF (numpatch > 0) allocate (acc_vec_pch (numpatch))
         IF (numelm   > 0) acc_vec     = 0.
         IF (numpatch > 0) acc_vec_pch = 0.
      ENDIF

      IF (DEF_hist_vars%riv_height) THEN
         IF (p_is_worker) THEN
            CALL worker_push_data (push_ucat2elm, a_wdsrf_ucat, acc_vec, fillvalue = spval)
            DO ielm = 1, numelm
               istt = elm_patch%substt(ielm)
               iend = elm_patch%subend(ielm)
               acc_vec_pch(istt:iend) = acc_vec(ielm)
            ENDDO
         ENDIF

         CALL flux_map_and_write_2d ( &
            acc_vec_pch, file_hist, 'f_wdpth_riv', itime_in_file, sumarea_ucat, filter_ucat, &
            'deepest water depth in river and flood plain', 'm')
      ENDIF

      IF (DEF_hist_vars%riv_veloct) THEN
         IF (p_is_worker) THEN
            CALL worker_push_data (push_ucat2elm, a_veloc_riv, acc_vec, fillvalue = spval)
            DO ielm = 1, numelm
               istt = elm_patch%substt(ielm)
               iend = elm_patch%subend(ielm)
               acc_vec_pch(istt:iend) = acc_vec(ielm)
            ENDDO
         ENDIF

         CALL flux_map_and_write_2d ( &
            acc_vec_pch, file_hist, 'f_veloc_riv', itime_in_file, sumarea_ucat, filter_ucat, &
            'water velocity in river', 'm/s')
      ENDIF

      IF (DEF_hist_vars%discharge) THEN
         IF (p_is_worker) THEN
            CALL worker_push_data (push_ucat2elm, a_discharge, acc_vec, fillvalue = spval)
            DO ielm = 1, numelm
               istt = elm_patch%substt(ielm)
               iend = elm_patch%subend(ielm)
               acc_vec_pch(istt:iend) = acc_vec(ielm)
            ENDDO
         ENDIF

         CALL flux_map_and_write_2d ( &
            acc_vec_pch, file_hist, 'f_discharge', itime_in_file, sumarea_ucat, filter_ucat, &
            'discharge in river and flood plain', 'm^3/s')
      ENDIF

      IF (DEF_hist_vars%floodarea) THEN
         IF (p_is_worker) THEN
            IF (numucat > 0) THEN
               allocate (a_floodfrc (numucat))
               WHERE (acctime > 0.)
                  a_floodfrc = a_floodarea  / topo_area
               END WHERE
            ENDIF

            IF (numelm > 0) then
               allocate (a_floodfrc_inpmat (nucpart,numelm))
            ENDIF

            CALL worker_push_data (push_ucat2inpmat, a_floodfrc, a_floodfrc_inpmat, fillvalue = spval)

            DO ielm = 1, numelm
               acc_vec(ielm) = sum(a_floodfrc_inpmat(:,ielm) * inpmat_area_u2e(:,ielm), &
                  mask = inpmat_area_u2e(:,ielm) > 0.)
               istt = elm_patch%substt(ielm)
               iend = elm_patch%subend(ielm)
               acc_vec_pch(istt:iend) = acc_vec(ielm)
            ENDDO

            IF (numucat > 0) deallocate (a_floodfrc       )
            IF (numelm > 0)  deallocate (a_floodfrc_inpmat)
         ENDIF

         CALL flux_map_and_write_2d ( &
            acc_vec_pch, file_hist, 'f_floodarea', itime_in_file, sumarea_inpmat, filter_inpmat, &
            'flooded area', 'm^2')

         IF (p_is_worker) THEN
            DO ielm = 1, numelm
               istt = elm_patch%substt(ielm)
               iend = elm_patch%subend(ielm)
               acc_vec_pch(istt:iend) = acc_vec(ielm) / sum(inpmat_area_u2e(:,ielm))
            ENDDO
         ENDIF

         CALL flux_map_and_write_2d ( &
            acc_vec_pch, file_hist, 'f_floodfrc', itime_in_file, sumarea_inpmat, filter_inpmat, &
            'flooded area fraction', '100%')
      ENDIF

      IF (allocated (acc_vec    )) deallocate (acc_vec    )
      IF (allocated (acc_vec_pch)) deallocate (acc_vec_pch)


      ! ----- reservoir variables -----
      IF (DEF_Reservoir_Method > 0) THEN
         IF (totalnumresv > 0) THEN

            IF (p_is_master) THEN

               i = len_trim (file_hist)
               DO WHILE (file_hist(i:i) /= '_')
                  i = i - 1
               ENDDO
               file_hist_resv = file_hist(1:i) // 'reservoir_' // file_hist(i+1:)

               inquire (file=file_hist_resv, exist=fexists)
               IF (.not. fexists) THEN
                  CALL ncio_define_dimension(file_hist_resv, 'reservoir', totalnumresv)
                  CALL ncio_write_serial (file_hist_resv, 'resv_GRAND_ID' , dam_GRAND_ID, 'reservoir')
                  CALL ncio_put_attr (file_hist_resv, 'resv_GRAND_ID', 'long_name', 'reservoir GRAND ID')
               ENDIF

               CALL ncio_write_time (file_hist_resv, 'time', idate, itime_in_file_resv, DEF_HIST_FREQ)

            ENDIF

            allocate (acc_vec (totalnumresv))
            acc_vec(:) = spval

            IF (DEF_hist_vars%volresv) THEN

               IF (p_is_worker) THEN
                  IF (numresv > 0) THEN
                     WHERE (acctime_resv > 0)
                        a_volresv = a_volresv / acctime_resv
                     ELSEWHERE
                        a_volresv = spval
                     END WHERE
                  ENDIF
               ENDIF

               CALL reservoir_gather_var (a_volresv, acc_vec)

               IF (p_is_master) THEN
                  CALL ncio_write_serial_time (file_hist_resv, 'volresv', &
                     itime_in_file_resv, acc_vec, 'reservoir', 'time', DEF_HIST_CompressLevel)
                  IF (itime_in_file_resv == 1) THEN
                     CALL ncio_put_attr (file_hist_resv, 'volresv', 'long_name', 'reservoir water volume')
                     CALL ncio_put_attr (file_hist_resv, 'volresv', 'units',     'm^3')
                     CALL ncio_put_attr (file_hist_resv, 'volresv', 'missing_value', spval)
                  ENDIF
               ENDIF
            ENDIF

            IF (DEF_hist_vars%qresv_in) THEN

               IF (p_is_worker) THEN
                  IF (numresv > 0) THEN
                     WHERE (acctime_resv > 0)
                        a_qresv_in = a_qresv_in / acctime_resv
                     ELSEWHERE
                        a_qresv_in = spval
                     END WHERE
                  ENDIF
               ENDIF

               CALL reservoir_gather_var (a_qresv_in, acc_vec)

               IF (p_is_master) THEN
                  CALL ncio_write_serial_time (file_hist_resv, 'qresv_in', &
                     itime_in_file_resv, acc_vec, 'reservoir', 'time', DEF_HIST_CompressLevel)
                  IF (itime_in_file_resv == 1) THEN
                     CALL ncio_put_attr (file_hist_resv, 'qresv_in', 'long_name', 'reservoir inflow')
                     CALL ncio_put_attr (file_hist_resv, 'qresv_in', 'units',     'm^3/s')
                     CALL ncio_put_attr (file_hist_resv, 'qresv_in', 'missing_value', spval)
                  ENDIF
               ENDIF
            ENDIF

            IF (DEF_hist_vars%qresv_out) THEN

               IF (p_is_worker) THEN
                  IF (numresv > 0) THEN
                     WHERE (acctime_resv > 0)
                        a_qresv_out = a_qresv_out / acctime_resv
                     ELSEWHERE
                        a_qresv_out = spval
                     END WHERE
                  ENDIF
               ENDIF

               CALL reservoir_gather_var (a_qresv_out, acc_vec)

               IF (p_is_master) THEN
                  CALL ncio_write_serial_time (file_hist_resv, 'qresv_out', &
                     itime_in_file_resv, acc_vec, 'reservoir', 'time', DEF_HIST_CompressLevel)
                  IF (itime_in_file_resv == 1) THEN
                     CALL ncio_put_attr (file_hist_resv, 'qresv_out', 'long_name', 'reservoir outflow')
                     CALL ncio_put_attr (file_hist_resv, 'qresv_out', 'units',     'm^3/s')
                     CALL ncio_put_attr (file_hist_resv, 'qresv_out', 'missing_value', spval)
                  ENDIF
               ENDIF
            ENDIF

            deallocate (acc_vec)

         ENDIF
      ENDIF

      CALL flush_acc_fluxes_riverlake ()

   END SUBROUTINE hist_grid_riverlake_out

   !-----------------------
   SUBROUTINE flush_acc_fluxes_riverlake ()

   USE MOD_SPMD_Task
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   USE MOD_Grid_Reservoir,        only: numresv
   IMPLICIT NONE

      IF (p_is_worker) THEN

         IF (numucat > 0) THEN
            acctime      (:) = 0.
            a_wdsrf_ucat (:) = 0.
            a_veloc_riv  (:) = 0.
            a_discharge  (:) = 0.
            a_floodarea  (:) = 0.
         ENDIF

         IF (numresv > 0) THEN
            acctime_resv (:) = 0.
            a_volresv    (:) = 0.
            a_qresv_in   (:) = 0.
            a_qresv_out  (:) = 0.
         ENDIF

      ENDIF

   END SUBROUTINE flush_acc_fluxes_riverlake

   !---------------------------------------
   SUBROUTINE hist_grid_riverlake_final ()

   IMPLICIT NONE

      IF (allocated(acctime      )) deallocate (acctime      )
      IF (allocated(a_wdsrf_ucat )) deallocate (a_wdsrf_ucat )
      IF (allocated(a_veloc_riv  )) deallocate (a_veloc_riv  )
      IF (allocated(a_discharge  )) deallocate (a_discharge  )
      IF (allocated(a_floodarea  )) deallocate (a_floodarea  )

      IF (allocated(acctime_resv )) deallocate (acctime_resv )
      IF (allocated(a_volresv    )) deallocate (a_volresv    )
      IF (allocated(a_qresv_in   )) deallocate (a_qresv_in   )
      IF (allocated(a_qresv_out  )) deallocate (a_qresv_out  )

      IF (allocated(filter_ucat  )) deallocate (filter_ucat  )
      IF (allocated(filter_inpmat)) deallocate (filter_inpmat)

   END SUBROUTINE hist_grid_riverlake_final

END MODULE MOD_Grid_RiverLakeHist
#endif
