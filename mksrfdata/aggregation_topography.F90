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
   USE mod_aggregation
   USE mod_utils

#ifdef SrfdataDiag
   USE mod_srfdata_diag
#endif

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
#ifdef SrfdataDiag
   INTEGER :: ityp
   INTEGER :: typindex(N_land_classification+1)
#endif

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
      CALL aggregation_data_daemon (gtopo, data_r8_2d_in1 = topography)
#endif
   ENDIF

   !   ---------------------------------------------------------------
   !   aggregate the elevation from the resolution of raw data to modelling resolution
   !   ---------------------------------------------------------------

   IF (p_is_worker) THEN

      allocate (topography_patches (numpatch))
   
      DO ipatch = 1, numpatch
         CALL aggregation_request_data (landpatch, ipatch, gtopo, area = area_one, &
            data_r8_2d_in1 = topography, data_r8_2d_out1 = topography_one)
         IF (any(topography_one /= -9999.0)) THEN
            topography_patches (ipatch) = &
               sum(topography_one * area_one, mask = topography_one /= -9999.0) &
               / sum(area_one, mask = topography_one /= -9999.0)
         ELSE
            topography_patches (ipatch) = -1.0e36
         ENDIF
      ENDDO
      
#ifdef USEMPI
      CALL aggregation_worker_done ()
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

#ifdef SrfdataDiag
   typindex = (/(ityp, ityp = 0, N_land_classification)/)
   lndname  = trim(dir_model_landdata) // '/diag/topo.nc'
   CALL srfdata_map_and_write (topography_patches, landpatch%settyp, typindex, m_patch2diag, &
      -1.0e36_r8, lndname, 'topography', compress = 0, write_mode = 'one')
#endif

   IF (p_is_worker) THEN
      deallocate ( topography_patches )
   ENDIF

END SUBROUTINE aggregation_topography
