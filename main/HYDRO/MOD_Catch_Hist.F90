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
   USE MOD_SPMD_Task
   USE MOD_NetCDFSerial
   USE MOD_Vars_Global, only: spval

   USE MOD_Mesh,    only: numelm
   USE MOD_LandHRU, only: numhru
   USE MOD_Catch_BasinNetwork, only: numbasin, numbsnhru
   USE MOD_Catch_Vars_1DFluxes
   USE MOD_Catch_IO

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

   END SUBROUTINE hist_basin_final

   !---------------------------------------
   SUBROUTINE hist_basin_out (file_hist, idate)

   USE MOD_Precision
   USE MOD_Namelist
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
   logical,  allocatable ::  filter(:)
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
      IF ((p_is_worker) .and. allocated(a_wdsrf_bsn)) THEN
         WHERE(a_wdsrf_bsn /= spval)
            a_wdsrf_bsn = a_wdsrf_bsn / nac_basin
         END WHERE
      ENDIF

      CALL worker_push_data (iam_bsn, iam_elm, a_wdsrf_bsn, a_wdsrf_elm)

      CALL vector_write_basin (&
         file_hist_basin, a_wdsrf_elm, numelm, totalnumelm, 'v_wdsrf_bsn', 'basin', elm_data_address, &
         DEF_hist_vars%riv_height, itime_in_file, 'River Height', 'm')

      ! ----- water velocity in river -----
      IF ((p_is_worker) .and. allocated(a_veloc_riv)) THEN
         WHERE(a_veloc_riv /= spval)
            a_veloc_riv = a_veloc_riv / nac_basin
         END WHERE
      ENDIF

      CALL worker_push_data (iam_bsn, iam_elm, a_veloc_riv, a_veloc_elm)

      CALL vector_write_basin (&
         file_hist_basin, a_veloc_elm, numelm, totalnumelm, 'v_veloc_riv', 'basin', elm_data_address, &
         DEF_hist_vars%riv_veloct, itime_in_file, 'River Velocity', 'm/s')

      ! ----- discharge in river -----
      IF ((p_is_worker) .and. allocated(a_discharge)) THEN
         WHERE(a_discharge /= spval)
            a_discharge = a_discharge / nac_basin
         END WHERE
      ENDIF

      CALL worker_push_data (iam_bsn, iam_elm, a_discharge, a_dschg_elm)

      CALL vector_write_basin (&
         file_hist_basin, a_dschg_elm, numelm, totalnumelm, 'v_discharge', 'basin', elm_data_address, &
         DEF_hist_vars%discharge, itime_in_file, 'River Discharge', 'm^3/s')

      ! ----- number of time steps for each basin -----
      CALL worker_push_data (iam_bsn, iam_elm, ntacc_bsn, ntacc_elm)

      CALL vector_write_basin (&
         file_hist_basin, ntacc_elm, numelm, totalnumelm, 'timesteps', 'basin', elm_data_address, &
         .true., itime_in_file, 'Number of accumulated timesteps for each basin', '-')

      IF (p_is_worker .and. (numbasin > 0)) ntacc_bsn(:) = 0.

      ! ----- water depth in hydro unit -----
      IF ((p_is_worker) .and. allocated(a_wdsrf_bsnhru)) THEN
         WHERE(a_wdsrf_bsnhru /= spval)
            a_wdsrf_bsnhru = a_wdsrf_bsnhru / nac_basin
         END WHERE
      ENDIF

      CALL worker_push_subset_data (iam_bsn, iam_elm, basin_hru, elm_hru, a_wdsrf_bsnhru, a_wdsrf_hru)

      CALL vector_write_basin (&
         file_hist_basin, a_wdsrf_hru, numhru, totalnumhru, 'v_wdsrf_hru', 'hydrounit', hru_data_address, &
         DEF_hist_vars%wdsrf_hru, itime_in_file, 'Depth of Surface Water in Hydro unit', 'm')

      ! ----- water velocity in hydro unit -----
      IF ((p_is_worker) .and. allocated(a_veloc_bsnhru)) THEN
         WHERE(a_veloc_bsnhru /= spval)
            a_veloc_bsnhru = a_veloc_bsnhru / nac_basin
         END WHERE
      ENDIF

      CALL worker_push_subset_data (iam_bsn, iam_elm, basin_hru, elm_hru, a_veloc_bsnhru, a_veloc_hru)

      CALL vector_write_basin (&
         file_hist_basin, a_veloc_hru, numhru, totalnumhru, 'v_veloc_hru', 'hydrounit', hru_data_address, &
         DEF_hist_vars%veloc_hru, itime_in_file, 'Surface Flow Velocity in Hydro unit', 'm/s')

      ! ----- subsurface water flow between elements -----
      IF (p_is_worker) THEN
         WHERE(a_xsubs_elm /= spval)
            a_xsubs_elm = a_xsubs_elm / nac_basin
         END WHERE
      ENDIF

      CALL vector_write_basin (&
         file_hist_basin, a_xsubs_elm, numelm, totalnumelm, 'v_xsubs_bsn', 'basin', elm_data_address, &
         DEF_hist_vars%xsubs_bsn, itime_in_file, 'Subsurface lateral flow between basins', 'm/s')

      ! ----- subsurface water flow between hydro units -----
      IF (p_is_worker) THEN
         WHERE(a_xsubs_hru /= spval)
            a_xsubs_hru = a_xsubs_hru / nac_basin
         END WHERE
      ENDIF

      CALL vector_write_basin (&
         file_hist_basin, a_xsubs_hru, numhru, totalnumhru, 'v_xsubs_hru', 'hydrounit', hru_data_address, &
         DEF_hist_vars%xsubs_hru, itime_in_file, 'SubSurface lateral flow between HRUs', 'm/s')


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
