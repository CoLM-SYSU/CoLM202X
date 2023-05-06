#include <define.h>

#ifdef USE_DEPTH_TO_BEDROCK

SUBROUTINE aggregation_dbedrock ( &
      gland, dir_rawdata, dir_model_landdata)

   USE precision
   USE mod_namelist
   USE spmd_task
   USE mod_grid
   USE mod_landpatch
   USE ncio_vector
   USE ncio_block
   USE mod_colm_debug
   USE mod_aggregation
#ifdef SinglePoint
   USE mod_single_srfdata
#endif

#ifdef SrfdataDiag
   USE mod_srfdata_diag
#endif

   IMPLICIT NONE
   ! arguments:

   TYPE(grid_type),  intent(in) :: gland
   CHARACTER(LEN=*), intent(in) :: dir_rawdata
   CHARACTER(LEN=*), intent(in) :: dir_model_landdata

   ! local variables:
   CHARACTER(len=256) :: landdir, lndname

   TYPE (block_data_real8_2d) :: dbedrock
   REAL(r8), allocatable :: dbedrock_patches(:)
   REAL(r8), allocatable :: dbedrock_one(:), area_one(:)
   INTEGER :: ipatch

   landdir = trim(dir_model_landdata) // '/dbedrock/'

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif
   IF (p_is_master) THEN
      write(*,'(/, A)') 'Aggregate depth to bedrock ...'
      CALL system('mkdir -p ' // trim(adjustl(landdir)))
   ENDIF
#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef SinglePoint
   IF (USE_SITE_dbedrock) THEN
      RETURN
   ENDIF
#endif

   ! ................................................
   ! ... (2) global depth to bedrock
   ! ................................................

   IF (p_is_io) THEN
      
      CALL allocate_block_data (gland, dbedrock)
         
      lndname = trim(dir_rawdata)//'/dbedrock.nc' 
      CALL ncio_read_block (lndname, 'dbedrock', gland, dbedrock)

#ifdef USEMPI
      CALL aggregation_data_daemon (gland, data_r8_2d_in_1 = dbedrock)
#endif
   ENDIF

   !   ---------------------------------------------------------------
   !   aggregate the depth to bedrock from the resolution of raw data to modelling resolution
   !   ---------------------------------------------------------------

   IF (p_is_worker) THEN

      allocate (dbedrock_patches (numpatch))
      
      DO ipatch = 1, numpatch
         CALL aggregation_request_data (landpatch, ipatch, gland, area = area_one, &
            data_r8_2d_in1 = dbedrock, data_r8_2d_out1 = dbedrock_one)
         dbedrock_patches (ipatch) = sum(dbedrock_one * area_one) / sum(area_one)
      ENDDO
   
#ifdef USEMPI
      CALL aggregation_worker_done ()
#endif
   ENDIF

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CLMDEBUG
   CALL check_vector_data ('dbedrock_patches ', dbedrock_patches, -9999.0)
#endif

   ! Write-out the depth of the pacth in the gridcell
#ifndef SinglePoint
   lndname = trim(landdir)//'/dbedrock_patches.nc'
   CALL ncio_create_file_vector (lndname, landpatch)
   CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
   CALL ncio_write_vector (lndname, 'dbedrock_patches', 'patch', landpatch, dbedrock_patches, 1)
#else
   SITE_dbedrock = dbedrock_patches(1)
#endif

   IF (p_is_worker) THEN
      deallocate ( dbedrock_patches )
   ENDIF

END SUBROUTINE aggregation_dbedrock

#endif
