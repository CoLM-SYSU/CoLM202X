#include <define.h>

#ifdef LATERAL_FLOW
module MOD_Hydro_Hist
   !--------------------------------------------------------------------------------
   ! DESCRIPTION:
   ! 
   !     Write out model results in lateral hydrological processes to history files.
   !
   ! Created by Shupeng Zhang, May 2023
   !--------------------------------------------------------------------------------

   use MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_NetCDFSerial
   USE MOD_Vars_Global,  only : spval

   USE MOD_Mesh,    only : numelm
   USE MOD_LandHRU, only : numhru
   USE MOD_Hydro_Vars_TimeVariables
   USE MOD_Hydro_Vars_1DFluxes
   USE MOD_Hydro_IO

   ! -- ACC Fluxes --
   INTEGER :: nac_basin

   REAL(r8), allocatable :: a_wdsrf_hru  (:)
   REAL(r8), allocatable :: a_veloc_hru  (:)

   REAL(r8), allocatable :: a_rsubs_bsn  (:)
   REAL(r8), allocatable :: a_rsubs_hru  (:)

   REAL(r8), allocatable :: a_height_riv (:)
   REAL(r8), allocatable :: a_veloct_riv (:)

   ! -- PUBLIC SUBROUTINEs --
   public :: hist_basin_init
   public :: hist_basin_out
   public :: hist_basin_final

!--------------------------------------------------------------------------
CONTAINS

   SUBROUTINE hist_basin_init

      IMPLICIT NONE

      INTEGER :: numbasin

      numbasin = numelm

      IF (p_is_worker) THEN
         IF (numhru > 0) THEN
            allocate ( a_wdsrf_hru  (numhru))
            allocate ( a_veloc_hru  (numhru))
            allocate ( a_rsubs_hru  (numhru))
         ENDIF

         IF (numbasin > 0) THEN
            allocate ( a_height_riv (numbasin))
            allocate ( a_veloct_riv (numbasin))
            allocate ( a_rsubs_bsn  (numbasin))
         ENDIF
      ENDIF

      call FLUSH_acc_fluxes_basin ()

   END SUBROUTINE hist_basin_init

   !--------------------------------------
   subroutine hist_basin_final ()

      implicit none

      IF (allocated(a_wdsrf_hru )) deallocate(a_wdsrf_hru )
      IF (allocated(a_veloc_hru )) deallocate(a_veloc_hru )

      IF (allocated(a_rsubs_bsn )) deallocate(a_rsubs_bsn )
      IF (allocated(a_rsubs_hru )) deallocate(a_rsubs_hru )

      IF (allocated(a_height_riv)) deallocate(a_height_riv)
      IF (allocated(a_veloct_riv)) deallocate(a_veloct_riv)

   end subroutine hist_basin_final

   !---------------------------------------
   SUBROUTINE hist_basin_out (file_hist, idate)

      use MOD_Precision
      use MOD_Namelist
      use MOD_SPMD_Task
      USE MOD_ElmVector
      USE MOD_HRUVector
      IMPLICIT NONE

      character(LEN=*), intent(in) :: file_hist
      integer,  INTENT(in) :: idate(3)

      ! Local variables
      CHARACTER(len=256) :: file_hist_basin
      LOGICAL :: fexists
      integer :: itime_in_file
      logical,  allocatable ::  filter(:)
      INTEGER :: numbasin, i

      if (p_is_master) then

         i = len_trim (file_hist)
         DO while (file_hist(i:i) /= '_')
            i = i - 1
         ENDDO
         file_hist_basin = file_hist(1:i) // 'basin_' // file_hist(i+1:)

         inquire (file=file_hist_basin, exist=fexists)
         if (.not. fexists) then
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
         endif

         call ncio_write_time (file_hist_basin, 'time', idate, itime_in_file, DEF_HIST_FREQ)

      ENDIF

      numbasin = numelm

      IF (p_is_worker) THEN
         where(a_height_riv /= spval)
            a_height_riv = a_height_riv / nac_basin
         END where
      ENDIF

      CALL vector_write_basin (&
         file_hist_basin, a_height_riv, numbasin, totalnumelm, 'wdsrf_bsn', 'basin', elm_data_address, &
         DEF_hist_vars%riv_height, itime_in_file, 'River Height', 'm')

      IF (p_is_worker) THEN
         where(a_veloct_riv /= spval)
            a_veloct_riv = a_veloct_riv / nac_basin
         END where
      ENDIF

      CALL vector_write_basin (&
         file_hist_basin, a_veloct_riv, numbasin, totalnumelm, 'veloc_riv', 'basin', elm_data_address, &
         DEF_hist_vars%riv_veloct, itime_in_file, 'River Velocity', 'm/s')

      IF (p_is_worker) THEN
         where(a_wdsrf_hru /= spval)
            a_wdsrf_hru = a_wdsrf_hru / nac_basin
         END where
      ENDIF

      CALL vector_write_basin (&
         file_hist_basin, a_wdsrf_hru, numhru, totalnumhru, 'wdsrf_hru', 'hydrounit', hru_data_address, &
         DEF_hist_vars%wdsrf_hru, itime_in_file, 'Depth of Surface Water in Hydro unit', 'm')

      IF (p_is_worker) THEN
         where(a_veloc_hru /= spval)
            a_veloc_hru = a_veloc_hru / nac_basin
         END where
      ENDIF

      CALL vector_write_basin (&
         file_hist_basin, a_veloc_hru, numhru, totalnumhru, 'veloc_hru', 'hydrounit', hru_data_address, &
         DEF_hist_vars%veloc_hru, itime_in_file, 'Surface Flow Velocity in Hydro unit', 'm/s')

      IF (p_is_worker) THEN
         where(a_rsubs_bsn /= spval)
            a_rsubs_bsn = a_rsubs_bsn / nac_basin
         END where
      ENDIF

      CALL vector_write_basin (&
         file_hist_basin, a_rsubs_bsn, numbasin, totalnumelm, 'rsubs_bsn', 'basin', elm_data_address, &
         DEF_hist_vars%rsubs_bsn, itime_in_file, 'Subsurface lateral flow between basins', 'm/s')

      IF (p_is_worker) THEN
         where(a_rsubs_hru /= spval)
            a_rsubs_hru = a_rsubs_hru / nac_basin
         END where
      ENDIF

      CALL vector_write_basin (&
         file_hist_basin, a_rsubs_hru, numhru, totalnumhru, 'rsubs_hru', 'hydrounit', hru_data_address, &
         DEF_hist_vars%rsubs_hru, itime_in_file, 'SubSurface lateral flow between HRUs', 'm/s')

      call FLUSH_acc_fluxes_basin ()

   END SUBROUTINE hist_basin_out

   !-----------------------
   SUBROUTINE FLUSH_acc_fluxes_basin ()

      use MOD_SPMD_Task
      USE MOD_Mesh,    only : numelm
      use MOD_LandHRU, only : numhru
      use MOD_Vars_Global,  only : spval 
      implicit none

      INTEGER :: numbasin

      if (p_is_worker) then

         numbasin = numelm

         nac_basin = 0

         IF (numbasin > 0) THEN
            a_height_riv(:) = spval
            a_veloct_riv(:) = spval
            a_rsubs_bsn (:) = spval
         ENDIF

         IF (numhru > 0) THEN
            a_wdsrf_hru(:) = spval
            a_veloc_hru(:) = spval
            a_rsubs_hru(:) = spval
         ENDIF

      ENDIF

   END SUBROUTINE FLUSH_acc_fluxes_basin

   ! -------
   SUBROUTINE accumulate_fluxes_basin

      IMPLICIT NONE

      INTEGER :: numbasin

      IF (p_is_worker) THEN

         nac_basin = nac_basin + 1

         numbasin = numelm

         IF (numbasin > 0) THEN
            CALL acc1d_basin (wdsrf_bsn_ta, a_height_riv)
            CALL acc1d_basin (veloc_riv_ta, a_veloct_riv)
            CALL acc1d_basin (rsubs_bsn     , a_rsubs_bsn )
         ENDIF

         IF (numhru > 0) THEN
            CALL acc1d_basin (wdsrf_hru_ta, a_wdsrf_hru)
            CALL acc1d_basin (veloc_hru_ta, a_veloc_hru)
            CALL acc1d_basin (rsubs_hru   , a_rsubs_hru)
         ENDIF
      ENDIF

   END SUBROUTINE accumulate_fluxes_basin

   ! -------
   SUBROUTINE acc1d_basin (var, s)
      
      use MOD_Precision
      use MOD_Vars_Global, only: spval

      IMPLICIT NONE

      real(r8), intent(in)    :: var(:)
      real(r8), intent(inout) :: s  (:)
      ! Local variables
      integer :: i

      do i = lbound(var,1), ubound(var,1)
         if (var(i) /= spval) then
            if (s(i) /= spval) then
               s(i) = s(i) + var(i)
            else
               s(i) = var(i)
            end if
         end if
      end do

   END SUBROUTINE acc1d_basin


end module MOD_Hydro_Hist
#endif
