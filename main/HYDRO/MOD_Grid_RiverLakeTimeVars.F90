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
   USE MOD_Vars_Global, only: spval
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
   USE MOD_Grid_RiverLakeNetwork, only: numucat, ucat_data_address
   USE MOD_Grid_Reservoir,        only: ucat2resv
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart

   ! -- auxiliary variables --
   real(r8), allocatable :: volresv_ucat (:)
   integer :: i

      IF (p_is_worker) THEN
         IF (numucat > 0)  allocate (volresv_ucat (numucat))
      ENDIF

      CALL vector_read_and_scatter (file_restart, wdsrf_ucat, numucat, 'wdsrf_ucat', ucat_data_address)
      CALL vector_read_and_scatter (file_restart, veloc_riv,  numucat, 'veloc_riv',  ucat_data_address)

      IF (DEF_Reservoir_Method > 0) THEN
         CALL vector_read_and_scatter (file_restart, volresv_ucat, numucat, 'volresv', ucat_data_address)
         IF (p_is_worker) THEN
            DO i = 1, numucat
               IF (ucat2resv(i) > 0) volresv(ucat2resv(i)) = volresv_ucat(i)
            ENDDO
         ENDIF
      ENDIF

      IF (allocated (volresv_ucat)) deallocate (volresv_ucat)

   END SUBROUTINE READ_GridRiverLakeTimeVars


   SUBROUTINE WRITE_GridRiverLakeTimeVars (file_restart)

   USE MOD_SPMD_Task
   USE MOD_Namelist
   USE MOD_NetCDFSerial
   USE MOD_Vector_ReadWrite
   USE MOD_Grid_RiverLakeNetwork, only: numucat, totalnumucat, ucat_data_address
   USE MOD_Grid_Reservoir,        only: ucat2resv
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart

   ! -- auxiliary variables --
   real(r8), allocatable :: volresv_ucat (:)
   integer :: i

      IF (p_is_worker) THEN
         IF (numucat > 0)  allocate (volresv_ucat (numucat))
      ENDIF

      IF (p_is_master) THEN
         CALL ncio_create_file (trim(file_restart))
         CALL ncio_define_dimension(file_restart, 'ucatch', totalnumucat)
      ENDIF

      CALL vector_gather_and_write (&
         file_restart, wdsrf_ucat, numucat, totalnumucat, 'wdsrf_ucat', 'ucatch', ucat_data_address)

      CALL vector_gather_and_write (&
         file_restart, veloc_riv, numucat, totalnumucat, 'veloc_riv', 'ucatch', ucat_data_address)

      IF (DEF_Reservoir_Method > 0) THEN
         IF (p_is_worker) THEN
            DO i = 1, numucat
               IF (ucat2resv(i) > 0) volresv_ucat(i) = volresv(ucat2resv(i))
            ENDDO
         ENDIF
         CALL vector_gather_and_write (&
            file_restart, volresv_ucat, numucat, totalnumucat, 'volresv', 'ucatch', ucat_data_address)
      ENDIF

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
