#include <define.h>

#ifdef CatchLateralFlow
MODULE MOD_Catch_Hist
!--------------------------------------------------------------------------------
! DESCRIPTION:
!
!     Write out model results in lateral hydrological processes to history files.
!
! Created by Shupeng Zhang, May 2023
!--------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_NetCDFSerial
   USE MOD_Vars_Global, only: spval

   USE MOD_Mesh,    only: numelm
   USE MOD_LandHRU, only: numhru
   USE MOD_Catch_BasinNetwork, only: numbasin, numbsnhru
   USE MOD_Catch_Reservoir
   USE MOD_Catch_Vars_1DFluxes
   USE MOD_Vector_ReadWrite

   ! -- ACC Fluxes --
   integer :: nac_basin

   real(r8), allocatable :: a_wdsrf_bsnhru (:)
   real(r8), allocatable :: a_veloc_bsnhru (:)

   real(r8), allocatable :: a_wdsrf_hru (:)
   real(r8), allocatable :: a_veloc_hru (:)

   real(r8), allocatable :: a_wdsrf_bsn (:)
   real(r8), allocatable :: a_veloc_riv (:)
   real(r8), allocatable :: a_discharge (:)

   real(r8), allocatable :: a_wdsrf_elm (:)
   real(r8), allocatable :: a_veloc_elm (:)
   real(r8), allocatable :: a_dschg_elm (:)

   real(r8), allocatable :: a_xsubs_elm (:)
   real(r8), allocatable :: a_xsubs_hru (:)

   real(r8), allocatable :: a_volresv   (:)
   real(r8), allocatable :: a_qresv_in  (:)
   real(r8), allocatable :: a_qresv_out (:)

   real(r8), allocatable :: ntacc_elm   (:)

   ! -- PUBLIC SUBROUTINEs --
   PUBLIC :: hist_basin_init
   PUBLIC :: hist_basin_out
   PUBLIC :: hist_basin_final

!--------------------------------------------------------------------------
CONTAINS

   SUBROUTINE hist_basin_init

   IMPLICIT NONE

      IF (p_is_worker) THEN
         IF (numbsnhru > 0) THEN
            allocate (a_wdsrf_bsnhru (numbsnhru))
            allocate (a_veloc_bsnhru (numbsnhru))
         ENDIF

         IF (numbasin > 0) THEN
            allocate (a_wdsrf_bsn (numbasin))
            allocate (a_veloc_riv (numbasin))
            allocate (a_discharge (numbasin))
         ENDIF

         IF (numelm > 0) allocate (a_xsubs_elm (numelm))
         IF (numhru > 0) allocate (a_xsubs_hru (numhru))

         IF (DEF_Reservoir_Method > 0) THEN
            IF (numresv > 0) allocate (a_volresv   (numresv))
            IF (numresv > 0) allocate (a_qresv_in  (numresv))
            IF (numresv > 0) allocate (a_qresv_out (numresv))
         ENDIF
      ENDIF

      CALL FLUSH_acc_fluxes_basin ()

   END SUBROUTINE hist_basin_init

   !--------------------------------------
   SUBROUTINE hist_basin_final ()

   IMPLICIT NONE

      IF (allocated(a_wdsrf_bsnhru)) deallocate(a_wdsrf_bsnhru)
      IF (allocated(a_veloc_bsnhru)) deallocate(a_veloc_bsnhru)

      IF (allocated(a_wdsrf_bsn)) deallocate(a_wdsrf_bsn)
      IF (allocated(a_veloc_riv)) deallocate(a_veloc_riv)
      IF (allocated(a_discharge)) deallocate(a_discharge)

      IF (allocated(a_xsubs_elm)) deallocate(a_xsubs_elm)
      IF (allocated(a_xsubs_hru)) deallocate(a_xsubs_hru)

      IF (allocated(a_volresv  )) deallocate(a_volresv  )
      IF (allocated(a_qresv_in )) deallocate(a_qresv_in )
      IF (allocated(a_qresv_out)) deallocate(a_qresv_out)

   END SUBROUTINE hist_basin_final

   !---------------------------------------
   SUBROUTINE hist_basin_out (file_hist, idate)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_ElmVector
   USE MOD_HRUVector
   USE MOD_LandHRU
   USE MOD_Catch_BasinNetwork
   IMPLICIT NONE

   character(len=*), intent(in) :: file_hist
   integer,  intent(in) :: idate(3)

   ! Local variables
   character(len=256) :: file_hist_basin
   logical :: fexists
   integer :: itime_in_file
   logical,  allocatable :: filter (:)
   real(r8), allocatable :: varhist(:)
   integer :: i

      IF (p_is_master) THEN

         i = len_trim (file_hist)
         DO WHILE (file_hist(i:i) /= '_')
            i = i - 1
         ENDDO
         file_hist_basin = file_hist(1:i) // 'basin_' // file_hist(i+1:)

         inquire (file=file_hist_basin, exist=fexists)
         IF (.not. fexists) THEN
            CALL ncio_create_file (trim(file_hist_basin))
            CALL ncio_define_dimension(file_hist_basin, 'time', 0)
            CALL ncio_define_dimension(file_hist_basin, 'basin',     totalnumelm)
            CALL ncio_define_dimension(file_hist_basin, 'hydrounit', totalnumhru)

            CALL ncio_write_serial (file_hist_basin, 'basin', eindex_glb, 'basin')
            CALL ncio_put_attr (file_hist_basin, 'basin', 'long_name', 'basin index')

            CALL ncio_write_serial (file_hist_basin, 'basin_hru', eindx_hru, 'hydrounit')
            CALL ncio_put_attr (file_hist_basin, 'basin_hru', 'long_name', &
               'basin index of hydrological units')

            CALL ncio_write_serial (file_hist_basin, 'hru_type' , htype_hru, 'hydrounit')
            CALL ncio_put_attr (file_hist_basin, 'hru_type' , 'long_name', &
               'index of hydrological units inside basin')

            IF (DEF_Reservoir_Method > 0) THEN
               IF (numresv_uniq > 0) THEN
                  CALL ncio_define_dimension(file_hist_basin, 'reservoir', numresv_uniq)
                  CALL ncio_write_serial (file_hist_basin, 'resv_hylak_id' , resv_hylak_id, 'reservoir')
                  CALL ncio_put_attr (file_hist_basin, 'resv_hylak_id' , 'long_name', &
                     'HydroLAKE ID of reservoirs')
               ENDIF
            ENDIF
         ENDIF

         CALL ncio_write_time (file_hist_basin, 'time', idate, itime_in_file, DEF_HIST_FREQ)

      ENDIF


      IF (p_is_worker) THEN
         IF (numhru > 0) THEN
            allocate (a_wdsrf_hru (numhru))
            allocate (a_veloc_hru (numhru))
         ENDIF

         IF (numelm > 0) THEN
            allocate (a_wdsrf_elm (numelm))
            allocate (a_veloc_elm (numelm))
            allocate (a_dschg_elm (numelm))
            allocate (ntacc_elm   (numelm))
         ENDIF
      ENDIF

      ! ----- water depth in basin -----
      IF (DEF_hist_vars%riv_height) THEN
         IF ((p_is_worker) .and. allocated(a_wdsrf_bsn)) THEN
            WHERE(a_wdsrf_bsn /= spval)
               a_wdsrf_bsn = a_wdsrf_bsn / nac_basin
            END WHERE
         ENDIF

         CALL worker_push_data (push_bsn2elm, a_wdsrf_bsn, a_wdsrf_elm, spval)

         CALL vector_gather_and_write (&
            a_wdsrf_elm, numelm, totalnumelm, elm_data_address, file_hist_basin, 'v_wdsrf_bsn', 'basin', &
            itime_in_file, 'River Height', 'm')
      ENDIF

      ! ----- water velocity in river -----
      IF (DEF_hist_vars%riv_veloct) THEN
         IF ((p_is_worker) .and. allocated(a_veloc_riv)) THEN
            WHERE(a_veloc_riv /= spval)
               a_veloc_riv = a_veloc_riv / nac_basin
            END WHERE
         ENDIF

         CALL worker_push_data (push_bsn2elm, a_veloc_riv, a_veloc_elm, spval)

         CALL vector_gather_and_write (&
            a_veloc_elm, numelm, totalnumelm, elm_data_address, file_hist_basin, 'v_veloc_riv', 'basin', &
            itime_in_file, 'River Velocity', 'm/s')
      ENDIF

      ! ----- discharge in river -----
      IF (DEF_hist_vars%discharge) THEN
         IF ((p_is_worker) .and. allocated(a_discharge)) THEN
            WHERE(a_discharge /= spval)
               a_discharge = a_discharge / nac_basin
            END WHERE
         ENDIF

         CALL worker_push_data (push_bsn2elm, a_discharge, a_dschg_elm, spval)

         CALL vector_gather_and_write (&
            a_dschg_elm, numelm, totalnumelm, elm_data_address, file_hist_basin, 'v_discharge', 'basin', &
            itime_in_file, 'River Discharge', 'm^3/s')
      ENDIF

      ! ----- number of time steps for each basin -----
      CALL worker_push_data (push_bsn2elm, ntacc_bsn, ntacc_elm, spval)

      CALL vector_gather_and_write (&
         ntacc_elm, numelm, totalnumelm, elm_data_address, file_hist_basin, 'timesteps', 'basin', &
         itime_in_file, 'Number of accumulated timesteps for each basin', '-')

      IF (p_is_worker .and. (numbasin > 0)) ntacc_bsn(:) = 0.

      ! ----- water depth in hydro unit -----
      IF (DEF_hist_vars%wdsrf_hru) THEN
         IF ((p_is_worker) .and. allocated(a_wdsrf_bsnhru)) THEN
            WHERE(a_wdsrf_bsnhru /= spval)
               a_wdsrf_bsnhru = a_wdsrf_bsnhru / nac_basin
            END WHERE
         ENDIF

         CALL worker_push_data (push_bsnhru2elmhru, a_wdsrf_bsnhru, a_wdsrf_hru, spval)

         CALL vector_gather_and_write (&
            a_wdsrf_hru, numhru, totalnumhru, hru_data_address, file_hist_basin, 'v_wdsrf_hru', 'hydrounit', &
            itime_in_file, 'Depth of Surface Water in Hydro unit', 'm')
      ENDIF

      ! ----- water velocity in hydro unit -----
      IF (DEF_hist_vars%veloc_hru) THEN
         IF ((p_is_worker) .and. allocated(a_veloc_bsnhru)) THEN
            WHERE(a_veloc_bsnhru /= spval)
               a_veloc_bsnhru = a_veloc_bsnhru / nac_basin
            END WHERE
         ENDIF

         CALL worker_push_data (push_bsnhru2elmhru, a_veloc_bsnhru, a_veloc_hru, spval)

         CALL vector_gather_and_write (&
            a_veloc_hru, numhru, totalnumhru, hru_data_address, file_hist_basin, 'v_veloc_hru', 'hydrounit', &
            itime_in_file, 'Surface Flow Velocity in Hydro unit', 'm/s')
      ENDIF

      ! ----- subsurface water flow between elements -----
      IF (DEF_hist_vars%xsubs_bsn) THEN
         IF (p_is_worker) THEN
            WHERE(a_xsubs_elm /= spval)
               a_xsubs_elm = a_xsubs_elm / nac_basin
            END WHERE
         ENDIF

         CALL vector_gather_and_write (&
            a_xsubs_elm, numelm, totalnumelm, elm_data_address, file_hist_basin, 'v_xsubs_bsn', 'basin', &
            itime_in_file, 'Subsurface lateral flow between basins', 'm/s')
      ENDIF

      ! ----- subsurface water flow between hydro units -----
      IF (DEF_hist_vars%xsubs_hru) THEN
         IF (p_is_worker) THEN
            WHERE(a_xsubs_hru /= spval)
               a_xsubs_hru = a_xsubs_hru / nac_basin
            END WHERE
         ENDIF

         CALL vector_gather_and_write (&
            a_xsubs_hru, numhru, totalnumhru, hru_data_address, file_hist_basin, 'v_xsubs_hru', 'hydrounit', &
            itime_in_file, 'SubSurface lateral flow between HRUs', 'm/s')
      ENDIF

      ! ----- reservoir variables -----
      IF (DEF_Reservoir_Method > 0) THEN
         IF (numresv_uniq > 0) THEN

            allocate (varhist (numresv_uniq))

            IF (DEF_hist_vars%volresv) THEN

               IF (p_is_worker) THEN
                  IF (numresv > 0) THEN
                     WHERE (a_volresv /= spval)
                        a_volresv = a_volresv / nac_basin
                     END WHERE
                  ENDIF
               ENDIF

               CALL reservoir_gather_var (a_volresv, varhist)

               IF (p_is_master) THEN
                  CALL ncio_write_serial_time (file_hist_basin, 'volresv', &
                     itime_in_file, varhist, 'reservoir', 'time', DEF_HIST_CompressLevel)
                  IF (itime_in_file == 1) THEN
                     CALL ncio_put_attr (file_hist_basin, 'volresv', 'long_name', 'reservoir water volume')
                     CALL ncio_put_attr (file_hist_basin, 'volresv', 'units',     'm^3')
                     CALL ncio_put_attr (file_hist_basin, 'volresv', 'missing_value', spval)
                  ENDIF
               ENDIF
            ENDIF

            IF (DEF_hist_vars%qresv_in) THEN

               IF (p_is_worker) THEN
                  IF (numresv > 0) THEN
                     WHERE (a_qresv_in /= spval)
                        a_qresv_in = a_qresv_in / nac_basin
                     END WHERE
                  ENDIF
               ENDIF

               CALL reservoir_gather_var (a_qresv_in, varhist)

               IF (p_is_master) THEN
                  CALL ncio_write_serial_time (file_hist_basin, 'qresv_in', &
                     itime_in_file, varhist, 'reservoir', 'time', DEF_HIST_CompressLevel)
                  IF (itime_in_file == 1) THEN
                     CALL ncio_put_attr (file_hist_basin, 'qresv_in', 'long_name', 'reservoir inflow')
                     CALL ncio_put_attr (file_hist_basin, 'qresv_in', 'units',     'm^3/s')
                     CALL ncio_put_attr (file_hist_basin, 'qresv_in', 'missing_value', spval)
                  ENDIF
               ENDIF
            ENDIF

            IF (DEF_hist_vars%qresv_out) THEN

               IF (p_is_worker) THEN
                  IF (numresv > 0) THEN
                     WHERE (a_qresv_out /= spval)
                        a_qresv_out = a_qresv_out / nac_basin
                     END WHERE
                  ENDIF
               ENDIF

               CALL reservoir_gather_var (a_qresv_out, varhist)

               IF (p_is_master) THEN
                  CALL ncio_write_serial_time (file_hist_basin, 'qresv_out', &
                     itime_in_file, varhist, 'reservoir', 'time', DEF_HIST_CompressLevel)
                  IF (itime_in_file == 1) THEN
                     CALL ncio_put_attr (file_hist_basin, 'qresv_out', 'long_name', 'reservoir outflow')
                     CALL ncio_put_attr (file_hist_basin, 'qresv_out', 'units',     'm^3/s')
                     CALL ncio_put_attr (file_hist_basin, 'qresv_out', 'missing_value', spval)
                  ENDIF
               ENDIF
            ENDIF

            deallocate (varhist)

         ENDIF
      ENDIF


      CALL FLUSH_acc_fluxes_basin ()

      IF (allocated(a_wdsrf_hru)) deallocate(a_wdsrf_hru)
      IF (allocated(a_veloc_hru)) deallocate(a_veloc_hru)

      IF (allocated(a_wdsrf_elm)) deallocate(a_wdsrf_elm)
      IF (allocated(a_veloc_elm)) deallocate(a_veloc_elm)
      IF (allocated(a_dschg_elm)) deallocate(a_dschg_elm)
      IF (allocated(ntacc_elm  )) deallocate(ntacc_elm  )


   END SUBROUTINE hist_basin_out

   !-----------------------
   SUBROUTINE FLUSH_acc_fluxes_basin ()

   USE MOD_SPMD_Task
   USE MOD_Mesh,    only: numelm
   USE MOD_LandHRU, only: numhru
   USE MOD_Vars_Global,  only: spval
   IMPLICIT NONE

      IF (p_is_worker) THEN

         nac_basin = 0

         IF (numbasin > 0) THEN
            a_wdsrf_bsn(:) = spval
            a_veloc_riv(:) = spval
            a_discharge(:) = spval
         ENDIF

         IF (numelm > 0) a_xsubs_elm (:) = spval

         IF (numbsnhru > 0) THEN
            a_wdsrf_bsnhru(:) = spval
            a_veloc_bsnhru(:) = spval
         ENDIF

         IF (numhru > 0) a_xsubs_hru(:) = spval

         IF (DEF_Reservoir_Method > 0) THEN
            IF (numresv > 0) a_volresv  (:) = spval
            IF (numresv > 0) a_qresv_in (:) = spval
            IF (numresv > 0) a_qresv_out(:) = spval
         ENDIF

      ENDIF

   END SUBROUTINE FLUSH_acc_fluxes_basin

   ! -------
   SUBROUTINE accumulate_fluxes_basin

   IMPLICIT NONE

      IF (p_is_worker) THEN

         nac_basin = nac_basin + 1

         IF (numbasin > 0) THEN
            CALL acc1d_basin (wdsrf_bsn_ta, a_wdsrf_bsn)
            CALL acc1d_basin (veloc_riv_ta, a_veloc_riv)
            CALL acc1d_basin (discharge_ta, a_discharge )
         ENDIF

         IF (numelm > 0) CALL acc1d_basin (xsubs_elm, a_xsubs_elm)

         IF (numbsnhru > 0) THEN
            CALL acc1d_basin (wdsrf_bsnhru_ta, a_wdsrf_bsnhru)
            CALL acc1d_basin (veloc_bsnhru_ta, a_veloc_bsnhru)
         ENDIF

         IF (numhru > 0) CALL acc1d_basin (xsubs_hru, a_xsubs_hru)

         IF (DEF_Reservoir_Method > 0) THEN
            IF (numresv > 0) CALL acc1d_basin (volresv_ta,   a_volresv  )
            IF (numresv > 0) CALL acc1d_basin (qresv_in_ta,  a_qresv_in )
            IF (numresv > 0) CALL acc1d_basin (qresv_out_ta, a_qresv_out)
         ENDIF

      ENDIF

   END SUBROUTINE accumulate_fluxes_basin

   ! -------
   SUBROUTINE acc1d_basin (var, s)

   USE MOD_Precision
   USE MOD_Vars_Global, only: spval

   IMPLICIT NONE

   real(r8), intent(in)    :: var(:)
   real(r8), intent(inout) :: s  (:)
   ! Local variables
   integer :: i

      DO i = lbound(var,1), ubound(var,1)
         IF (var(i) /= spval) THEN
            IF (s(i) /= spval) THEN
               s(i) = s(i) + var(i)
            ELSE
               s(i) = var(i)
            ENDIF
         ENDIF
      ENDDO

   END SUBROUTINE acc1d_basin

END MODULE MOD_Catch_Hist
#endif
