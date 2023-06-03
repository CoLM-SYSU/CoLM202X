#include <define.h>

#ifdef LATERAL_FLOW
MODULE MOD_Hydro_Vars_TimeVariables

   USE MOD_Precision
   IMPLICIT NONE


   ! -- state variables --
   REAL(r8), allocatable :: riverheight (:)
   REAL(r8), allocatable :: riverveloct (:)
   REAL(r8), allocatable :: dpond_hru   (:) ! depth of ponding water [m]
   REAL(r8), allocatable :: veloc_hru   (:) ! [m/s]
   REAL(r8), allocatable :: zwt_hru     (:)

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: allocate_HydroTimeVariables
   PUBLIC :: deallocate_HydroTimeVariables

   PUBLIC :: read_HydroTimeVariables
   PUBLIC :: write_HydroTimeVariables

CONTAINS

  SUBROUTINE allocate_HydroTimeVariables

     USE MOD_SPMD_Task
     USE MOD_Mesh
     USE MOD_LandHRU
     IMPLICIT NONE

     INTEGER :: numbasin

     numbasin = numelm

     IF (p_is_worker) THEN
        IF (numbasin > 0) THEN
           allocate (riverheight (numbasin))
           allocate (riverveloct (numbasin))
        ENDIF

        IF (numhru > 0) THEN
           allocate (dpond_hru (numhru))
           allocate (veloc_hru (numhru))
           allocate (zwt_hru   (numhru))
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

     INTEGER :: numbasin
     character(LEN=*), intent(in) :: file_restart

     numbasin = numelm

     CALL vector_read_basin (file_restart, riverheight, numbasin, 'riverheight', elm_data_address)
     CALL vector_read_basin (file_restart, riverveloct, numbasin, 'riverveloct', elm_data_address)

     CALL vector_read_basin (file_restart, dpond_hru, numhru, 'dpond_hru', hru_data_address)
     CALL vector_read_basin (file_restart, veloc_hru, numhru, 'veloc_hru', hru_data_address)
     CALL vector_read_basin (file_restart, zwt_hru  , numhru, 'zwt_hru'  , hru_data_address)

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

     INTEGER :: numbasin, iwork
     character(LEN=*), intent(in) :: file_restart

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
        file_restart, riverheight, numbasin, totalnumelm, 'riverheight', 'basin', elm_data_address)

     CALL vector_write_basin (&
        file_restart, riverveloct, numbasin, totalnumelm, 'riverveloct', 'basin', elm_data_address)

     CALL vector_write_basin (&
        file_restart, dpond_hru, numhru, totalnumhru, 'dpond_hru', 'hydrounit', hru_data_address)

     CALL vector_write_basin (&
        file_restart, veloc_hru, numhru, totalnumhru, 'veloc_hru', 'hydrounit', hru_data_address)

     CALL vector_write_basin (&
        file_restart, zwt_hru, numhru, totalnumhru, 'zwt_hru', 'hydrounit', hru_data_address)

  END SUBROUTINE WRITE_HydroTimeVariables

  SUBROUTINE deallocate_HydroTimeVariables

     IMPLICIT NONE

     IF (allocated(riverheight)) deallocate(riverheight)
     IF (allocated(riverveloct)) deallocate(riverveloct)

     IF (allocated(dpond_hru)) deallocate(dpond_hru)
     IF (allocated(veloc_hru)) deallocate(veloc_hru)
     IF (allocated(zwt_hru  )) deallocate(zwt_hru  )

  END SUBROUTINE deallocate_HydroTimeVariables

END MODULE MOD_Hydro_Vars_TimeVariables
#endif
