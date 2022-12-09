#include <define.h>

MODULE mod_mesh_filter

   USE mod_grid
   IMPLICIT NONE

   LOGICAL :: has_mesh_filter
   TYPE(grid_type) :: grid_filter

CONTAINS

   LOGICAL FUNCTION inquire_mesh_filter ()

      USE spmd_task
      USE mod_namelist
      IMPLICIT NONE
      LOGICAL :: fexists
      
      IF (p_is_master) THEN
      
         inquire (file=trim(DEF_file_mesh_filter), exist=fexists)

         IF (.not. fexists) THEN
            write(*,'(/, 2A)') 'Mesh Filter not used: file ', trim(DEF_file_mesh_filter)
         ELSE
            write(*,'(/, 2A)') 'Mesh Filter from file ', trim(DEF_file_mesh_filter)
         ENDIF
      ENDIF
      
#ifdef USEMPI
      call mpi_bcast (fexists, 1, MPI_LOGICAL, p_root, p_comm_glb, p_err)
#endif

      inquire_mesh_filter = fexists

   END FUNCTION inquire_mesh_filter

   ! -------------
   SUBROUTINE mesh_filter ()
   
      USE precision
      USE mod_namelist
      USE spmd_task
      USE ncio_block
      USE mod_data_type
      USE mod_mesh
      USE mod_aggregation_generic
      USE mod_block
      IMPLICIT NONE
   
      ! local variables:
      ! ---------------------------------------------------------------
      TYPE (block_data_int32_2d) :: datafilter
      INTEGER, allocatable :: ifilter(:), xtemp(:), ytemp(:)
      LOGICAL, allocatable :: filter(:)
      INTEGER :: ielm, jelm, npxl, nelm_glb
   
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         write(*,'(/, A)') 'Filtering pixels ...'
      ENDIF
   
      IF (p_is_io) THEN
         CALL allocate_block_data (grid_filter, datafilter)
         CALL ncio_read_block (DEF_file_mesh_filter, 'mesh_filter', grid_filter, datafilter)
   
#ifdef USEMPI
         CALL aggregation_gen_data_daemon (grid_filter, data_i4 = datafilter)
#endif
      ENDIF
   
      IF (p_is_worker) THEN
   
         jelm = 0
         DO ielm = 1, numelm
            CALL aggregation_gen_request_data (grid_filter, &
               mesh(ielm)%ilon, mesh(ielm)%ilat, &
               data_i4 = datafilter, out_i4 = ifilter)
   
            allocate (filter (mesh(ielm)%npxl))
            filter = ifilter > 0
   
            IF (any(filter)) THEN
               jelm = jelm + 1
               IF (.not. all(filter)) THEN
   
                  npxl = count(filter)
   
                  allocate (xtemp(npxl))
                  allocate (ytemp(npxl))
                  xtemp = pack(mesh(ielm)%ilon, filter)
                  ytemp = pack(mesh(ielm)%ilat, filter)
   
                  deallocate(mesh(ielm)%ilon)
                  deallocate(mesh(ielm)%ilat)

                  mesh(ielm)%npxl = npxl

                  allocate(mesh(ielm)%ilon(npxl))
                  allocate(mesh(ielm)%ilat(npxl))
                  mesh(ielm)%ilon = xtemp
                  mesh(ielm)%ilat = ytemp 
   
                  deallocate (xtemp)
                  deallocate (ytemp)
               ENDIF
   
               IF (jelm /= ielm) THEN
                  CALL copy_elm (mesh(ielm), mesh(jelm))
               ENDIF
   
            ENDIF
               
            deallocate (filter)
         ENDDO
   
         numelm = jelm
   
#ifdef USEMPI
         CALL aggregation_gen_worker_done ()
#endif
      ENDIF
      
#ifdef USEMPI
      IF (p_is_worker) THEN
         CALL mpi_reduce (numelm, nelm_glb, 1, MPI_INTEGER, MPI_SUM, p_root, p_comm_worker, p_err)
         IF (p_iam_worker == 0) THEN
            write(*,'(A,I12,A)') 'Total: ', nelm_glb, ' elements after mesh filtering.'
         ENDIF
      ENDIF
#else
      write(*,'(A,I12,A)') 'Total: ', numelm, ' elements after mesh filtering.'
#endif
      
      ! Update nelm_blk
      nelm_blk(:,:) = 0
      IF (p_is_worker) THEN 
         DO ielm = 1, numelm
            nelm_blk(mesh(ielm)%xblk,mesh(ielm)%yblk) = &
               nelm_blk(mesh(ielm)%xblk,mesh(ielm)%yblk) + 1
         ENDDO
      ENDIF 
#ifdef USEMPI
      CALL mpi_allreduce (MPI_IN_PLACE, nelm_blk, gblock%nxblk*gblock%nyblk, &
         MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
   
   END SUBROUTINE mesh_filter

END MODULE mod_mesh_filter
