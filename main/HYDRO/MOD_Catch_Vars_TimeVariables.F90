#include <define.h>

#ifdef CatchLateralFlow
MODULE MOD_Catch_Vars_TimeVariables
!-------------------------------------------------------------------------------------
! DESCRIPTION:
!
!   Time Variables in lateral hydrological processes.
!
! Created by Shupeng Zhang, May 2023
!-------------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Catch_BasinNetwork
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   ! -- state variables (1): necessary for restart --
   real(r8), allocatable :: veloc_elm (:) ! river velocity [m/s]
   real(r8), allocatable :: veloc_hru (:) ! surface water velocity [m/s]
   real(r8), allocatable :: wdsrf_hru (:) ! surface water depth [m]

   real(r8), allocatable :: wdsrf_elm_prev (:) ! river or lake water depth at previous time step [m]
   real(r8), allocatable :: wdsrf_hru_prev (:) ! surface water depth at previous time step [m]

   ! -- state variables (2): only in model --
   real(r8), allocatable :: wdsrf_bsn    (:) ! river or lake water depth [m]
   real(r8), allocatable :: veloc_riv    (:) ! river velocity [m/s]
   real(r8), allocatable :: momen_riv    (:) ! unit river momentum [m^2/s]
   real(r8), allocatable :: wdsrf_bsnhru (:) ! surface water depth [m]
   real(r8), allocatable :: veloc_bsnhru (:) ! surface water velocity [m/s]
   real(r8), allocatable :: momen_bsnhru (:) ! unit surface water momentum [m^2/s]

   real(r8), allocatable :: wdsrf_bsn_prev    (:) ! river or lake water depth at previous time step [m]
   real(r8), allocatable :: wdsrf_bsnhru_prev (:) ! surface water depth at previous time step [m]

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: allocate_CatchTimeVariables
   PUBLIC :: deallocate_CatchTimeVariables

   PUBLIC :: read_CatchTimeVariables
   PUBLIC :: write_CatchTimeVariables

CONTAINS

   SUBROUTINE allocate_CatchTimeVariables

   USE MOD_SPMD_Task
   USE MOD_Mesh,    only: numelm
   USE MOD_LandHRU, only: numhru
   IMPLICIT NONE

      IF (p_is_worker) THEN

         IF (numelm > 0) THEN
            allocate (veloc_elm      (numelm))
            allocate (wdsrf_elm_prev (numelm))
         ENDIF

         IF (numhru > 0) THEN
            allocate (veloc_hru      (numhru))
            allocate (wdsrf_hru      (numhru))
            allocate (wdsrf_hru_prev (numhru))
         ENDIF

         IF (numbasin > 0)   allocate (wdsrf_bsn     (numbasin))
         IF (numbasin > 0)   allocate (veloc_riv     (numbasin))
         IF (numbasin > 0)   allocate (momen_riv     (numbasin))
         IF (numbasin > 0)   allocate (wdsrf_bsn_prev(numbasin))

         IF (numbsnhru > 0)  allocate (wdsrf_bsnhru      (numbsnhru))
         IF (numbsnhru > 0)  allocate (veloc_bsnhru      (numbsnhru))
         IF (numbsnhru > 0)  allocate (momen_bsnhru      (numbsnhru))
         IF (numbsnhru > 0)  allocate (wdsrf_bsnhru_prev (numbsnhru))

      ENDIF

   END SUBROUTINE allocate_CatchTimeVariables

   SUBROUTINE READ_CatchTimeVariables (file_restart)

   USE MOD_Mesh
   USE MOD_LandHRU
   USE MOD_Vector_ReadWrite
   USE MOD_ElmVector
   USE MOD_HRUVector
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart

      CALL vector_read_and_scatter (file_restart, veloc_elm, numelm, 'veloc_riv', elm_data_address)
      CALL worker_push_data (push_elm2bsn, veloc_elm, veloc_riv, spval)

      CALL vector_read_and_scatter (file_restart, wdsrf_elm_prev, numelm, 'wdsrf_bsn_prev', elm_data_address)
      CALL worker_push_data (push_elm2bsn, wdsrf_elm_prev, wdsrf_bsn_prev, spval)

      CALL vector_read_and_scatter (file_restart, veloc_hru, numhru, 'veloc_hru', hru_data_address)
      CALL worker_push_data (push_elmhru2bsnhru, veloc_hru, veloc_bsnhru, spval)

      CALL vector_read_and_scatter (file_restart, wdsrf_hru_prev, numhru, 'wdsrf_hru_prev', hru_data_address)
      CALL worker_push_data (push_elmhru2bsnhru, wdsrf_hru_prev, wdsrf_bsnhru_prev, spval)

   END SUBROUTINE READ_CatchTimeVariables

   SUBROUTINE WRITE_CatchTimeVariables (file_restart)

   USE MOD_SPMD_Task
   USE MOD_NetCDFSerial
   USE MOD_Mesh
   USE MOD_LandHRU
   USE MOD_Vector_ReadWrite
   USE MOD_ElmVector
   USE MOD_HRUVector
   IMPLICIT NONE

   integer :: iwork
   character(len=*), intent(in) :: file_restart

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

      CALL worker_push_data (push_bsn2elm, veloc_riv, veloc_elm, spval)
      CALL vector_gather_and_write (&
         veloc_elm, numelm, totalnumelm, elm_data_address, file_restart, 'veloc_riv', 'basin')

      CALL worker_push_data (push_bsn2elm, wdsrf_bsn_prev, wdsrf_elm_prev, spval)
      CALL vector_gather_and_write (&
         wdsrf_elm_prev, numelm, totalnumelm, elm_data_address, file_restart, 'wdsrf_bsn_prev', 'basin')

      CALL worker_push_data (push_bsnhru2elmhru, veloc_bsnhru, veloc_hru, spval)
      CALL vector_gather_and_write (&
         veloc_hru, numhru, totalnumhru, hru_data_address, file_restart, 'veloc_hru', 'hydrounit')

      CALL worker_push_data (push_bsnhru2elmhru, wdsrf_bsnhru_prev, wdsrf_hru_prev, spval)
      CALL vector_gather_and_write (&
         wdsrf_hru_prev, numhru, totalnumhru, hru_data_address, file_restart, 'wdsrf_hru_prev', 'hydrounit')

   END SUBROUTINE WRITE_CatchTimeVariables

   SUBROUTINE deallocate_CatchTimeVariables

   IMPLICIT NONE

      IF (allocated(veloc_elm)) deallocate(veloc_elm)
      IF (allocated(veloc_hru)) deallocate(veloc_hru)
      IF (allocated(wdsrf_hru)) deallocate(wdsrf_hru)

      IF (allocated(wdsrf_elm_prev)) deallocate(wdsrf_elm_prev)
      IF (allocated(wdsrf_hru_prev)) deallocate(wdsrf_hru_prev)

      IF (allocated (wdsrf_bsn     )) deallocate (wdsrf_bsn     )
      IF (allocated (veloc_riv     )) deallocate (veloc_riv     )
      IF (allocated (momen_riv     )) deallocate (momen_riv     )
      IF (allocated (wdsrf_bsn_prev)) deallocate (wdsrf_bsn_prev)

      IF (allocated (wdsrf_bsnhru     )) deallocate (wdsrf_bsnhru     )
      IF (allocated (veloc_bsnhru     )) deallocate (veloc_bsnhru     )
      IF (allocated (momen_bsnhru     )) deallocate (momen_bsnhru     )
      IF (allocated (wdsrf_bsnhru_prev)) deallocate (wdsrf_bsnhru_prev)

   END SUBROUTINE deallocate_CatchTimeVariables

END MODULE MOD_Catch_Vars_TimeVariables
#endif
