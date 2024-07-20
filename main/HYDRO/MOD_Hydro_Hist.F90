#include <define.h>

#ifdef CatchLateralFlow
MODULE MOD_Hydro_Hist
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
   USE MOD_Vars_Global,  only : spval

   USE MOD_Mesh,    only : numelm
   USE MOD_LandHRU, only : numhru
   USE MOD_Hydro_Vars_TimeVariables
   USE MOD_Hydro_Vars_1DFluxes
   USE MOD_Hydro_IO

   ! -- ACC Fluxes --
   integer :: nac_basin

   real(r8), allocatable :: a_wdsrf_hru  (:)
   real(r8), allocatable :: a_veloc_hru  (:)

   real(r8), allocatable :: a_xsubs_bsn  (:)
   real(r8), allocatable :: a_xsubs_hru  (:)

   real(r8), allocatable :: a_height_riv (:)
   real(r8), allocatable :: a_veloct_riv (:)
   real(r8), allocatable :: a_discharge  (:)

   ! -- PUBLIC SUBROUTINEs --
   PUBLIC :: hist_basin_init
   PUBLIC :: hist_basin_out
   PUBLIC :: hist_basin_final

!--------------------------------------------------------------------------
CONTAINS

   SUBROUTINE hist_basin_init

   IMPLICIT NONE

   integer :: numbasin

      numbasin = numelm

      IF (p_is_worker) THEN
         IF (numhru > 0) THEN
            allocate ( a_wdsrf_hru  (numhru))
            allocate ( a_veloc_hru  (numhru))
            allocate ( a_xsubs_hru  (numhru))
         ENDIF

         IF (numbasin > 0) THEN
            allocate ( a_height_riv (numbasin))
            allocate ( a_veloct_riv (numbasin))
            allocate ( a_discharge  (numbasin))
            allocate ( a_xsubs_bsn  (numbasin))
         ENDIF
      ENDIF

      CALL FLUSH_acc_fluxes_basin ()

   END SUBROUTINE hist_basin_init

   !--------------------------------------
   SUBROUTINE hist_basin_final ()

   IMPLICIT NONE

      IF (allocated(a_wdsrf_hru )) deallocate(a_wdsrf_hru )
      IF (allocated(a_veloc_hru )) deallocate(a_veloc_hru )

      IF (allocated(a_xsubs_bsn )) deallocate(a_xsubs_bsn )
      IF (allocated(a_xsubs_hru )) deallocate(a_xsubs_hru )

      IF (allocated(a_height_riv)) deallocate(a_height_riv)
      IF (allocated(a_veloct_riv)) deallocate(a_veloct_riv)
      IF (allocated(a_discharge )) deallocate(a_discharge )

   END SUBROUTINE hist_basin_final

   !---------------------------------------
   SUBROUTINE hist_basin_out (file_hist, idate)

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_ElmVector
   USE MOD_HRUVector
   IMPLICIT NONE

   character(len=*), intent(in) :: file_hist
   integer,  intent(in) :: idate(3)

   ! Local variables
   character(len=256) :: file_hist_basin
   logical :: fexists
   integer :: itime_in_file
   logical,  allocatable ::  filter(:)
   integer :: numbasin, i

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

      numbasin = numelm

      IF (p_is_worker) THEN
         WHERE(a_height_riv /= spval)
            a_height_riv = a_height_riv / nac_basin
         END WHERE
      ENDIF

      CALL vector_write_basin (&
         file_hist_basin, a_height_riv, numbasin, totalnumelm, 'wdsrf_bsn', 'basin', elm_data_address, &
         DEF_hist_vars%riv_height, itime_in_file, 'River Height', 'm')

      IF (p_is_worker) THEN
         WHERE(a_veloct_riv /= spval)
            a_veloct_riv = a_veloct_riv / nac_basin
         END WHERE
      ENDIF

      CALL vector_write_basin (&
         file_hist_basin, a_veloct_riv, numbasin, totalnumelm, 'veloc_riv', 'basin', elm_data_address, &
         DEF_hist_vars%riv_veloct, itime_in_file, 'River Velocity', 'm/s')

      IF (p_is_worker) THEN
         WHERE(a_discharge /= spval)
            a_discharge = a_discharge / nac_basin
         END WHERE
      ENDIF

      CALL vector_write_basin (&
         file_hist_basin, a_discharge, numbasin, totalnumelm, 'discharge', 'basin', elm_data_address, &
         DEF_hist_vars%discharge, itime_in_file, 'River Discharge', 'm^3/s')

      IF (p_is_worker) THEN
         WHERE(a_wdsrf_hru /= spval)
            a_wdsrf_hru = a_wdsrf_hru / nac_basin
         END WHERE
      ENDIF

      CALL vector_write_basin (&
         file_hist_basin, a_wdsrf_hru, numhru, totalnumhru, 'wdsrf_hru', 'hydrounit', hru_data_address, &
         DEF_hist_vars%wdsrf_hru, itime_in_file, 'Depth of Surface Water in Hydro unit', 'm')

      IF (p_is_worker) THEN
         WHERE(a_veloc_hru /= spval)
            a_veloc_hru = a_veloc_hru / nac_basin
         END WHERE
      ENDIF

      CALL vector_write_basin (&
         file_hist_basin, a_veloc_hru, numhru, totalnumhru, 'veloc_hru', 'hydrounit', hru_data_address, &
         DEF_hist_vars%veloc_hru, itime_in_file, 'Surface Flow Velocity in Hydro unit', 'm/s')

      IF (p_is_worker) THEN
         WHERE(a_xsubs_bsn /= spval)
            a_xsubs_bsn = a_xsubs_bsn / nac_basin
         END WHERE
      ENDIF

      CALL vector_write_basin (&
         file_hist_basin, a_xsubs_bsn, numbasin, totalnumelm, 'xsubs_bsn', 'basin', elm_data_address, &
         DEF_hist_vars%xsubs_bsn, itime_in_file, 'Subsurface lateral flow between basins', 'm/s')

      IF (p_is_worker) THEN
         WHERE(a_xsubs_hru /= spval)
            a_xsubs_hru = a_xsubs_hru / nac_basin
         END WHERE
      ENDIF

      CALL vector_write_basin (&
         file_hist_basin, a_xsubs_hru, numhru, totalnumhru, 'xsubs_hru', 'hydrounit', hru_data_address, &
         DEF_hist_vars%xsubs_hru, itime_in_file, 'SubSurface lateral flow between HRUs', 'm/s')

      CALL FLUSH_acc_fluxes_basin ()

   END SUBROUTINE hist_basin_out

   !-----------------------
   SUBROUTINE FLUSH_acc_fluxes_basin ()

   USE MOD_SPMD_Task
   USE MOD_Mesh,    only : numelm
   USE MOD_LandHRU, only : numhru
   USE MOD_Vars_Global,  only : spval
   IMPLICIT NONE

   integer :: numbasin

      IF (p_is_worker) THEN

         numbasin = numelm

         nac_basin = 0

         IF (numbasin > 0) THEN
            a_height_riv(:) = spval
            a_veloct_riv(:) = spval
            a_discharge (:) = spval
            a_xsubs_bsn (:) = spval
         ENDIF

         IF (numhru > 0) THEN
            a_wdsrf_hru(:) = spval
            a_veloc_hru(:) = spval
            a_xsubs_hru(:) = spval
         ENDIF

      ENDIF

   END SUBROUTINE FLUSH_acc_fluxes_basin

   ! -------
   SUBROUTINE accumulate_fluxes_basin

   IMPLICIT NONE

   integer :: numbasin

      IF (p_is_worker) THEN

         nac_basin = nac_basin + 1

         numbasin = numelm

         IF (numbasin > 0) THEN
            CALL acc1d_basin (wdsrf_bsn_ta, a_height_riv)
            CALL acc1d_basin (veloc_riv_ta, a_veloct_riv)
            CALL acc1d_basin (discharge   , a_discharge )
            CALL acc1d_basin (xsubs_bsn   , a_xsubs_bsn )
         ENDIF

         IF (numhru > 0) THEN
            CALL acc1d_basin (wdsrf_hru_ta, a_wdsrf_hru)
            CALL acc1d_basin (veloc_hru_ta, a_veloc_hru)
            CALL acc1d_basin (xsubs_hru   , a_xsubs_hru)
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

END MODULE MOD_Hydro_Hist
#endif
