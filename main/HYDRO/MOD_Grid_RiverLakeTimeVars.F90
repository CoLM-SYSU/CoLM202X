#include <define.h>

#ifdef GridRiverLakeFlow
MODULE MOD_Grid_RiverLakeTimeVars
!-------------------------------------------------------------------------------------
! DESCRIPTION:
!
!   Time Variables in gridded hydrological processes.
!
! Created by Shupeng Zhang, Oct 2025
!-------------------------------------------------------------------------------------

   USE MOD_Precision
   IMPLICIT NONE

   ! -- state variables --
   real(r8), allocatable :: wdsrf_ucat (:) ! river or lake water depth [m]
   real(r8), allocatable :: veloc_riv  (:) ! river velocity            [m/s]
   real(r8), allocatable :: momen_riv  (:) ! unit river momentum       [m^2/s]
   real(r8), allocatable :: volresv    (:) ! reservoir water volume    [m^3]

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: allocate_GridRiverLakeTimeVars
   PUBLIC :: deallocate_GridRiverLakeTimeVars

   PUBLIC :: read_GridRiverLakeTimeVars
   PUBLIC :: write_GridRiverLakeTimeVars

CONTAINS

   SUBROUTINE allocate_GridRiverLakeTimeVars

   USE MOD_SPMD_Task
   USE MOD_Grid_RiverLakeNetwork, only: numucat
   USE MOD_Grid_Reservoir,        only: numresv
   IMPLICIT NONE

      IF (p_is_worker) THEN

         IF (numucat > 0)  allocate (wdsrf_ucat (numucat))
         IF (numucat > 0)  allocate (veloc_riv  (numucat))
         IF (numucat > 0)  allocate (momen_riv  (numucat))
         IF (numresv > 0)  allocate (volresv    (numresv))

      ENDIF

   END SUBROUTINE allocate_GridRiverLakeTimeVars


   SUBROUTINE READ_GridRiverLakeTimeVars (file_restart)

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_Vector_ReadWrite
   USE MOD_WorkerPushData
   USE MOD_ElmVector,             only: elm_data_address, totalnumelm
   USE MOD_Mesh,                  only: numelm
   USE MOD_Grid_RiverLakeNetwork, only: numucat, push_elm2ucat
   USE MOD_Grid_Reservoir,        only: ucat2resv
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart

   ! -- auxiliary variables --
   real(r8), allocatable :: wdsrf_elm    (:)
   real(r8), allocatable :: veloc_elm    (:)
   real(r8), allocatable :: volresv_elm  (:)
   real(r8), allocatable :: volresv_ucat (:)
   integer :: i

      IF (p_is_worker) THEN
         IF (numelm  > 0)  allocate (wdsrf_elm    (numelm ))
         IF (numelm  > 0)  allocate (veloc_elm    (numelm ))
         IF (numelm  > 0)  allocate (volresv_elm  (numelm ))
         IF (numucat > 0)  allocate (volresv_ucat (numucat))
      ENDIF

      CALL vector_read_and_scatter (file_restart, wdsrf_elm, numelm, 'wdsrf_ucat', elm_data_address)
      CALL worker_push_data (push_elm2ucat, wdsrf_elm, wdsrf_ucat)

      CALL vector_read_and_scatter (file_restart, veloc_elm, numelm, 'veloc_riv', elm_data_address)
      CALL worker_push_data (push_elm2ucat, veloc_elm, veloc_riv)

      IF (DEF_Reservoir_Method > 0) THEN
         CALL vector_read_and_scatter (file_restart, volresv_elm, numelm, 'volresv', elm_data_address)
         CALL worker_push_data (push_elm2ucat, volresv_elm, volresv_ucat)
         IF (p_is_worker) THEN
            DO i = 1, numucat
               IF (ucat2resv(i) > 0) volresv(ucat2resv(i)) = volresv_ucat(i)
            ENDDO
         ENDIF
      ENDIF

      IF (allocated (veloc_elm   )) deallocate (veloc_elm   )
      IF (allocated (wdsrf_elm   )) deallocate (wdsrf_elm   )
      IF (allocated (volresv_elm )) deallocate (volresv_elm )
      IF (allocated (volresv_ucat)) deallocate (volresv_ucat)

   END SUBROUTINE READ_GridRiverLakeTimeVars


   SUBROUTINE WRITE_GridRiverLakeTimeVars (file_restart)

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_NetCDFSerial
   USE MOD_WorkerPushData
   USE MOD_Vector_ReadWrite
   USE MOD_ElmVector,             only: totalnumelm, eindex_glb, elm_data_address
   USE MOD_Mesh,                  only: numelm
   USE MOD_Grid_RiverLakeNetwork, only: numucat, push_ucat2elm
   USE MOD_Grid_Reservoir,        only: ucat2resv
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart

   ! -- auxiliary variables --
   real(r8), allocatable :: wdsrf_elm    (:)
   real(r8), allocatable :: veloc_elm    (:)
   real(r8), allocatable :: volresv_elm  (:)
   real(r8), allocatable :: volresv_ucat (:)
   integer :: i

      IF (p_is_worker) THEN
         IF (numelm  > 0)  allocate (wdsrf_elm    (numelm ))
         IF (numelm  > 0)  allocate (veloc_elm    (numelm ))
         IF (numelm  > 0)  allocate (volresv_elm  (numelm ))
         IF (numucat > 0)  allocate (volresv_ucat (numucat))
      ENDIF

      IF (p_is_master) THEN
         CALL ncio_create_file (trim(file_restart))
         CALL ncio_define_dimension(file_restart, 'ucatch', totalnumelm)

         CALL ncio_write_serial (file_restart, 'ucatch', eindex_glb, 'ucatch')
         CALL ncio_put_attr (file_restart, 'ucatch', 'long_name', 'element index of unit catchment')
      ENDIF

      CALL worker_push_data (push_ucat2elm, wdsrf_ucat, wdsrf_elm)
      CALL vector_gather_and_write (&
         file_restart, wdsrf_elm, numelm, totalnumelm, 'wdsrf_ucat', 'ucatch', elm_data_address)

      CALL worker_push_data (push_ucat2elm, veloc_riv, veloc_elm)
      CALL vector_gather_and_write (&
         file_restart, veloc_elm, numelm, totalnumelm, 'veloc_riv', 'ucatch', elm_data_address)

      IF (DEF_Reservoir_Method > 0) THEN
         IF (p_is_worker) THEN
            DO i = 1, numucat
               IF (ucat2resv(i) > 0) volresv_ucat(i) = volresv(ucat2resv(i))
            ENDDO
         ENDIF
         CALL worker_push_data (push_ucat2elm, volresv_ucat, volresv_elm)
         CALL vector_gather_and_write (&
            file_restart, volresv_elm, numelm, totalnumelm, 'volresv', 'ucatch', elm_data_address)
      ENDIF

      IF (allocated (wdsrf_elm   )) deallocate (wdsrf_elm   )
      IF (allocated (veloc_elm   )) deallocate (veloc_elm   )
      IF (allocated (volresv_elm )) deallocate (volresv_elm )
      IF (allocated (volresv_ucat)) deallocate (volresv_ucat)

   END SUBROUTINE WRITE_GridRiverLakeTimeVars

   SUBROUTINE deallocate_GridRiverLakeTimeVars

   IMPLICIT NONE

      IF (allocated (wdsrf_ucat)) deallocate (wdsrf_ucat)
      IF (allocated (veloc_riv )) deallocate (veloc_riv )
      IF (allocated (momen_riv )) deallocate (momen_riv )
      IF (allocated (volresv   )) deallocate (volresv   )

   END SUBROUTINE deallocate_GridRiverLakeTimeVars

END MODULE MOD_Grid_RiverLakeTimeVars
#endif
