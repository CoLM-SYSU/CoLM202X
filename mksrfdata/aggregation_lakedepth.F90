#include <define.h>

SUBROUTINE aggregation_lakedepth ( &
      gland, dir_rawdata, dir_model_landdata)

   ! ----------------------------------------------------------------------
   ! 1. Global land cover types (updated with the specific dataset)
   !
   ! 2. Global Lake Coverage and Lake Depth
   !   (http://nwpi.krc.karelia.run/flake/)      
   !    Kourzeneva, E., H. Asensio, E. Martin, and S. Faroux, 2012: Global
   !    gridded dataset of lake coverage and lake depth for USE in numerical
   !    weather prediction and climate modelling. Tellus A, 64, 15640.
   !
   !    Lake depth data legend
   !    Value   Description
   ! 0       no lake indicated in this pixel
   ! 1       no any information about this lake and set the default value of 10 m
   ! 2       no information about depth for this lake and set the default value of 10 m
   ! 3       have the information about lake depth in this pixel
   ! 4       this is the river pixel according to our map, set the default value of 3 m
   !
   ! Created by Yongjiu Dai, 02/2014
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
   ! ---------------------------------------------------------------
   CHARACTER(len=256) :: landdir, lndname
   INTEGER :: L, ipatch

   TYPE (block_data_real8_2d) :: lakedepth
   REAL(r8), allocatable :: lakedepth_patches(:), lakedepth_one(:)
#ifdef SrfdataDiag
   INTEGER :: typindex(1) = (/17/)   
#endif

   landdir = trim(dir_model_landdata) // '/lakedepth/'

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif
   IF (p_is_master) THEN
      write(*,'(/, A)') 'Aggregate lake depth ...'
      CALL system('mkdir -p ' // trim(adjustl(landdir)))
   ENDIF
#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef SinglePoint
   IF (USE_SITE_lakedepth) THEN
      RETURN
   ENDIF
#endif

   ! ................................................
   ! ... (2) global lake coverage and lake depth
   ! ................................................
   lndname = trim(dir_rawdata)//'/lake_depth.nc' 

   IF (p_is_io) THEN
      CALL allocate_block_data (gland, lakedepth)
      CALL ncio_read_block (lndname, 'lake_depth', gland, lakedepth)
      CALL block_data_linear_transform (lakedepth, scl = 0.1)

#ifdef USEMPI
      CALL aggregation_data_daemon (gland, data_r8_2d_in1 = lakedepth)
#endif
   ENDIF

   !   ---------------------------------------------------------------
   !   aggregate the lake depth from the resolution of raw data to modelling resolution
   !   ---------------------------------------------------------------

   IF (p_is_worker) THEN

      allocate (lakedepth_patches (numpatch))
   
      DO ipatch = 1, numpatch
         L = landpatch%settyp(ipatch)
#ifdef USGS_CLASSIFICATION
         IF(L==16)THEN  ! LAND WATER BODIES (16)
#endif
#ifdef IGBP_CLASSIFICATION
         IF(L==17)THEN  ! LAND WATER BODIES (17)
#endif
#ifdef PFT_CLASSIFICATION
         IF(L==17)THEN  ! LAND WATER BODIES (17)
#endif
#ifdef PC_CLASSIFICATION
         IF(L==17)THEN  ! LAND WATER BODIES (17)
#endif
            CALL aggregation_request_data (landpatch, ipatch, gland, &
               data_r8_2d_in1 = lakedepth, data_r8_2d_out1 = lakedepth_one)
            lakedepth_patches (ipatch) = median (lakedepth_one, size(lakedepth_one))
         ELSE
            lakedepth_patches (ipatch) = -1.0e36_r8
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
   CALL check_vector_data ('lakedepth_patches ', lakedepth_patches)
#endif

   ! Write-out the lake depth of the lake pacth in the gridcell
#ifndef SinglePoint
   lndname = trim(landdir)//'/lakedepth_patches.nc'
   CALL ncio_create_file_vector (lndname, landpatch)
   CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
   CALL ncio_write_vector (lndname, 'lakedepth_patches', 'patch', landpatch, lakedepth_patches, 1)

#ifdef SrfdataDiag
   lndname = trim(dir_model_landdata)//'/diag/lakedepth.nc'
   CALL srfdata_map_and_write (lakedepth_patches, landpatch%settyp, typindex, m_patch2diag, &
      -1.0e36_r8, lndname, 'lakedepth', compress = 0, write_mode = 'one')
#endif
#else
   SITE_lakedepth = lakedepth_patches(1)
#endif


   IF (p_is_worker) THEN
      deallocate ( lakedepth_patches )
   ENDIF

END SUBROUTINE aggregation_lakedepth
