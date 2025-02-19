#include <define.h>

MODULE MOD_LandElm

!-----------------------------------------------------------------------
! !DESCRIPTION:
!
!    Build pixelset "landelm".
!
!    In CoLM, the global/regional area is divided into a hierarchical structure:
!    1. If GRIDBASED or UNSTRUCTURED is defined, it is
!       ELEMENT >>> PATCH
!    2. If CATCHMENT is defined, it is
!       ELEMENT >>> HRU >>> PATCH
!    If Plant Function Type classification is used, PATCH is further divided into PFT.
!    If Plant Community classification is used,     PATCH is further divided into PC.
!
!    "landelm" refers to pixelset ELEMENT.
!
!  Created by Shupeng Zhang, May 2023
!-----------------------------------------------------------------------

   USE MOD_Pixelset
   IMPLICIT NONE

   ! ---- Instance ----
   type(pixelset_type) :: landelm

CONTAINS

   ! -------------------------------
   SUBROUTINE landelm_build

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Mesh
   IMPLICIT NONE

   ! Local Variables
   integer :: ielm, nelm_glb

      IF (p_is_master) THEN
         write(*,'(A)') 'Making land elements:'
      ENDIF

      IF (p_is_worker) THEN

         allocate (landelm%eindex (numelm))
         allocate (landelm%ipxstt (numelm))
         allocate (landelm%ipxend (numelm))
         allocate (landelm%settyp (numelm))
         allocate (landelm%ielm   (numelm))

         DO ielm = 1, numelm
            landelm%eindex(ielm) = mesh(ielm)%indx
            landelm%ipxstt(ielm) = 1
            landelm%ipxend(ielm) = mesh(ielm)%npxl
            landelm%settyp(ielm) = 0
            landelm%ielm  (ielm) = ielm
         ENDDO

      ENDIF

      landelm%nset = numelm
      CALL landelm%set_vecgs

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_worker) THEN
         CALL mpi_reduce (numelm, nelm_glb, 1, MPI_INTEGER, MPI_SUM, p_root, p_comm_worker, p_err)
         IF (p_iam_worker == 0) THEN
            write(*,'(A,I12,A)') 'Total: ', nelm_glb, ' elements.'
         ENDIF
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      write(*,'(A,I12,A)') 'Total: ', numelm, ' elements.'
#endif

   END SUBROUTINE landelm_build

END MODULE MOD_LandElm
