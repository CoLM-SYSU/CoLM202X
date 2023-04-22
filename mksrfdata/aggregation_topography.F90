#include <define.h>

SUBROUTINE aggregation_topography ( &
      gtopo, dir_rawdata, dir_model_landdata)

   ! ----------------------------------------------------------------------
   ! 1. topography
   ! ----------------------------------------------------------------------
   USE precision
   USE mod_namelist
   USE spmd_task
   USE mod_grid
   USE mod_landpatch
   USE ncio_vector
   USE ncio_block
#ifdef CLMDEBUG 
   USE mod_colm_debug
#endif
   USE mod_aggregation_lc
   USE mod_utils

   IMPLICIT NONE
   ! arguments:

   TYPE(grid_type),  intent(in) :: gtopo
   CHARACTER(LEN=*), intent(in) :: dir_rawdata
   CHARACTER(LEN=*), intent(in) :: dir_model_landdata

   ! local variables:
   ! ---------------------------------------------------------------
   CHARACTER(len=256) :: landdir, lndname
   INTEGER :: ipatch

   TYPE (block_data_real8_2d) :: topography
   REAL(r8), allocatable :: topography_patches(:)
   REAL(r8), allocatable :: topography_one(:), area_one(:)

   landdir = trim(dir_model_landdata) // '/topography/'

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif
   IF (p_is_master) THEN
      write(*,'(/, A)') 'Aggregate topography ...'
      CALL system('mkdir -p ' // trim(adjustl(landdir)))
   ENDIF
#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef SinglePoint
   RETURN
#endif

   lndname = trim(dir_rawdata)//'/elevation.nc' 

   IF (p_is_io) THEN
      CALL allocate_block_data (gtopo, topography)
      CALL ncio_read_block (lndname, 'elevation', gtopo, topography)

#ifdef USEMPI
      CALL aggregation_lc_data_daemon (gtopo, topography)
#endif
   ENDIF

   !   ---------------------------------------------------------------
   !   aggregate the elevation from the resolution of raw data to modelling resolution
   !   ---------------------------------------------------------------

   IF (p_is_worker) THEN

      allocate (topography_patches (numpatch))
   
      DO ipatch = 1, numpatch
         CALL aggregation_lc_request_data (ipatch, gtopo, topography, topography_one, area_one)
         topography_patches (ipatch) = sum(topography_one * area_one) / sum(area_one)
      ENDDO
      
#ifdef USEMPI
      CALL aggregation_lc_worker_done ()
#endif
   ENDIF

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CLMDEBUG
   CALL check_vector_data ('topography_patches ', topography_patches)
#endif

   lndname = trim(landdir)//'/topography_patches.nc'
   CALL ncio_create_file_vector (lndname, landpatch)
   CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
   CALL ncio_write_vector (lndname, 'topography_patches', 'patch', landpatch, topography_patches, 1)

   IF (p_is_worker) THEN
      deallocate ( topography_patches )
   ENDIF

END SUBROUTINE aggregation_topography
