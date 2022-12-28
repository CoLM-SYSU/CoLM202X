#include <define.h>

MODULE mod_landelm

   USE precision
   USE mod_pixelset
   USE mod_grid
   IMPLICIT NONE

   ! ---- Instance ----
   TYPE(pixelset_type) :: landelm

CONTAINS

   ! -------------------------------
   SUBROUTINE landelm_build 

      USE precision
      USE spmd_task
      USE mod_mesh
      IMPLICIT NONE

      ! Local Variables
      INTEGER :: ielm, nelm_glb

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

END MODULE mod_landelm
