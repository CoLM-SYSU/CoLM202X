#include <define.h>

#ifdef LATERAL_FLOW
MODULE MOD_HydroTimeVars
   
   USE precision
   IMPLICIT NONE


   ! -- state variables --
   REAL(r8), allocatable :: riverheight (:)
   REAL(r8), allocatable :: riverveloct (:)
   REAL(r8), allocatable :: dpond_hru   (:) ! depth of ponding water [m]
   REAL(r8), allocatable :: veloc_hru   (:) ! [m/s]
   REAL(r8), allocatable :: zwt_hru     (:)

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: allocate_HydroTimeVars
   PUBLIC :: deallocate_HydroTimeVars

   PUBLIC :: read_HydroTimeVars
   PUBLIC :: write_HydroTimeVars

CONTAINS 

  SUBROUTINE allocate_HydroTimeVars

     USE spmd_task
     USE mod_mesh
     USE mod_landhru
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

  END SUBROUTINE allocate_HydroTimeVars

  SUBROUTINE READ_HydroTimeVars (file_restart)

     USE mod_mesh
     USE mod_landhru
     USE mod_io_basin
     USE mod_elm_vector
     USE mod_hru_vector
     IMPLICIT NONE
     
     INTEGER :: numbasin
     character(LEN=*), intent(in) :: file_restart

     numbasin = numelm

     CALL vector_read_basin (file_restart, riverheight, numbasin, 'riverheight', elm_data_address)
     CALL vector_read_basin (file_restart, riverveloct, numbasin, 'riverveloct', elm_data_address)
     
     CALL vector_read_basin (file_restart, dpond_hru, numhru, 'dpond_hru', hru_data_address)
     CALL vector_read_basin (file_restart, veloc_hru, numhru, 'veloc_hru', hru_data_address)
     CALL vector_read_basin (file_restart, zwt_hru  , numhru, 'zwt_hru'  , hru_data_address)

  END SUBROUTINE READ_HydroTimeVars 

  SUBROUTINE WRITE_HydroTimeVars (file_restart)

     USE spmd_task
     USE ncio_serial
     USE mod_mesh
     USE mod_landhru
     USE mod_io_basin
     USE mod_elm_vector
     USE mod_hru_vector
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

  END SUBROUTINE WRITE_HydroTimeVars

  SUBROUTINE deallocate_HydroTimeVars

     IMPLICIT NONE

     IF (allocated(riverheight)) deallocate(riverheight)
     IF (allocated(riverveloct)) deallocate(riverveloct)

     IF (allocated(dpond_hru)) deallocate(dpond_hru)
     IF (allocated(veloc_hru)) deallocate(veloc_hru)
     IF (allocated(zwt_hru  )) deallocate(zwt_hru  )

  END SUBROUTINE deallocate_HydroTimeVars

END MODULE MOD_HydroTimeVars
#endif
