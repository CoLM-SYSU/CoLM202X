#include <define.h>

#ifdef CatchLateralFlow
MODULE MOD_Hydro_Vars_TimeVariables
!-------------------------------------------------------------------------------------
! DESCRIPTION:
!
!   Time Variables in lateral hydrological processes.
!
! Created by Shupeng Zhang, May 2023
!-------------------------------------------------------------------------------------

   USE MOD_Precision
   IMPLICIT NONE

   ! -- state variables --
   real(r8), allocatable :: wdsrf_bsn (:) ! river or lake water depth [m]
   real(r8), allocatable :: veloc_riv (:) ! river velocity [m/s]
   real(r8), allocatable :: momen_riv (:) ! unit river momentum [m^2/s]
   real(r8), allocatable :: wdsrf_hru (:) ! surface water depth [m]
   real(r8), allocatable :: veloc_hru (:) ! surface water velocity [m/s]
   real(r8), allocatable :: momen_hru (:) ! unit surface water momentum [m^2/s]

   real(r8), allocatable :: wdsrf_bsn_prev (:) ! river or lake water depth at previous time step [m]
   real(r8), allocatable :: wdsrf_hru_prev (:) ! surface water depth at previous time step [m]

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: allocate_HydroTimeVariables
   PUBLIC :: deallocate_HydroTimeVariables

   PUBLIC :: read_HydroTimeVariables
   PUBLIC :: write_HydroTimeVariables

CONTAINS

   SUBROUTINE allocate_HydroTimeVariables

   USE MOD_SPMD_Task
   USE MOD_Mesh,    only : numelm
   USE MOD_LandHRU, only : numhru
   IMPLICIT NONE

   integer :: numbasin

      numbasin = numelm

      IF (p_is_worker) THEN
         IF (numbasin > 0) THEN
            allocate (wdsrf_bsn (numbasin))
            allocate (veloc_riv (numbasin))
            allocate (momen_riv (numbasin))
            allocate (wdsrf_bsn_prev (numbasin))
         ENDIF

         IF (numhru > 0) THEN
            allocate (wdsrf_hru (numhru))
            allocate (veloc_hru (numhru))
            allocate (momen_hru (numhru))
            allocate (wdsrf_hru_prev (numhru))
         ENDIF
      ENDIF

   END SUBROUTINE allocate_HydroTimeVariables

   SUBROUTINE READ_HydroTimeVariables (file_restart)

   USE MOD_Mesh
   USE MOD_LandHRU
   USE MOD_Hydro_IO
   USE MOD_ElmVector
   USE MOD_HRUVector
   IMPLICIT NONE

   integer :: numbasin
   character(len=*), intent(in) :: file_restart

      numbasin = numelm

      CALL vector_read_basin (file_restart, wdsrf_bsn, numbasin, 'wdsrf_bsn', elm_data_address)
      CALL vector_read_basin (file_restart, veloc_riv, numbasin, 'veloc_riv', elm_data_address)
      CALL vector_read_basin (file_restart, wdsrf_bsn_prev, numbasin, 'wdsrf_bsn_prev', elm_data_address)

      CALL vector_read_basin (file_restart, wdsrf_hru, numhru, 'wdsrf_hru', hru_data_address)
      CALL vector_read_basin (file_restart, veloc_hru, numhru, 'veloc_hru', hru_data_address)
      CALL vector_read_basin (file_restart, wdsrf_hru_prev, numhru, 'wdsrf_hru_prev', hru_data_address)

   END SUBROUTINE READ_HydroTimeVariables

   SUBROUTINE WRITE_HydroTimeVariables (file_restart)

   USE MOD_SPMD_Task
   USE MOD_NetCDFSerial
   USE MOD_Mesh
   USE MOD_LandHRU
   USE MOD_Hydro_IO
   USE MOD_ElmVector
   USE MOD_HRUVector
   IMPLICIT NONE

   integer :: numbasin, iwork
   character(len=*), intent(in) :: file_restart

      numbasin = numelm

      IF (p_is_master) THEN
         CALL ncio_create_file (trim(file_restart))
         CALL ncio_define_dimension(file_restart, 'basin',     totalnumelm)
         CALL ncio_define_dimension(file_restart, 'hydrounit', totalnumhru)

         CALL ncio_write_serial (file_restart, 'basin', eindex_glb, 'basin')
         CALL ncio_put_attr (file_restart, 'basin', 'long_name', 'basin index')

         CALL ncio_write_serial (file_restart, 'bsn_hru', eindx_hru, 'hydrounit')
         CALL ncio_put_attr (file_restart, 'bsn_hru', &
            'long_name', 'basin index of hydrological units')

         CALL ncio_write_serial (file_restart, 'hru_type' , htype_hru, 'hydrounit')
         CALL ncio_put_attr (file_restart, 'hru_type' , &
            'long_name', 'index of hydrological units inside basin')
      ENDIF

      CALL vector_write_basin (&
         file_restart, wdsrf_bsn, numbasin, totalnumelm, 'wdsrf_bsn', 'basin', elm_data_address)

      CALL vector_write_basin (&
         file_restart, veloc_riv, numbasin, totalnumelm, 'veloc_riv', 'basin', elm_data_address)

      CALL vector_write_basin (&
         file_restart, wdsrf_hru, numhru, totalnumhru, 'wdsrf_hru', 'hydrounit', hru_data_address)

      CALL vector_write_basin (&
         file_restart, veloc_hru, numhru, totalnumhru, 'veloc_hru', 'hydrounit', hru_data_address)

      CALL vector_write_basin (&
         file_restart, wdsrf_bsn_prev, numbasin, totalnumelm, 'wdsrf_bsn_prev', 'basin', elm_data_address)

      CALL vector_write_basin (&
         file_restart, wdsrf_hru_prev, numhru, totalnumhru, 'wdsrf_hru_prev', 'hydrounit', hru_data_address)

   END SUBROUTINE WRITE_HydroTimeVariables

   SUBROUTINE deallocate_HydroTimeVariables

   IMPLICIT NONE

      IF (allocated(wdsrf_bsn)) deallocate(wdsrf_bsn)
      IF (allocated(veloc_riv)) deallocate(veloc_riv)
      IF (allocated(momen_riv)) deallocate(momen_riv)

      IF (allocated(wdsrf_hru)) deallocate(wdsrf_hru)
      IF (allocated(veloc_hru)) deallocate(veloc_hru)
      IF (allocated(momen_hru)) deallocate(momen_hru)

      IF (allocated(wdsrf_bsn_prev)) deallocate(wdsrf_bsn_prev)
      IF (allocated(wdsrf_hru_prev)) deallocate(wdsrf_hru_prev)

   END SUBROUTINE deallocate_HydroTimeVariables

END MODULE MOD_Hydro_Vars_TimeVariables
#endif
