#include <define.h>

SUBROUTINE aggregation_soil_parameters ( &
      gland, dir_rawdata, dir_model_landdata)
   ! ----------------------------------------------------------------------
   ! Creates land model surface dataset from original "raw" data files -
   !     data with 30 arc seconds resolution
   !
   ! Created by Yongjiu Dai, 02/2014
   ! ----------------------------------------------------------------------
   USE precision
   USE GlobalVars
   USE mod_namelist
   USE spmd_task
   USE mod_grid
   USE mod_landpatch
   USE ncio_block
   USE ncio_vector
   USE mod_aggregation_lc
   USE mod_colm_debug
   USE mod_utils

   IMPLICIT NONE

   ! arguments:
   TYPE(grid_type),  intent(in) :: gland
   CHARACTER(LEN=*), intent(in) :: dir_rawdata
   CHARACTER(LEN=*), intent(in) :: dir_model_landdata

   ! local variables:
   ! ---------------------------------------------------------------
   CHARACTER(len=256) :: lndname
   CHARACTER(len=256) :: c
   INTEGER :: nsl, ipatch, L

   TYPE (block_data_real8_2d) :: theta_s_grid
   TYPE (block_data_real8_2d) :: psi_s_grid  
#ifdef Campbell_SOIL_MODEL
   TYPE (block_data_real8_2d) :: lambda_grid 
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
   TYPE (block_data_real8_2d) :: theta_r_grid 
   TYPE (block_data_real8_2d) :: alpha_vgm_grid 
   TYPE (block_data_real8_2d) :: L_vgm_grid 
   TYPE (block_data_real8_2d) :: n_vgm_grid 
#endif
   TYPE (block_data_real8_2d) :: k_s_grid    
   TYPE (block_data_real8_2d) :: csol_grid   
   TYPE (block_data_real8_2d) :: tksatu_grid 
   TYPE (block_data_real8_2d) :: tkdry_grid  

   REAL(r8), allocatable :: theta_s_patches (:) 
   REAL(r8), allocatable :: psi_s_patches   (:) 
#ifdef Campbell_SOIL_MODEL
   REAL(r8), allocatable :: lambda_patches  (:) 
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
   REAL(r8), allocatable :: theta_r_patches (:) 
   REAL(r8), allocatable :: alpha_vgm_patches  (:) 
   REAL(r8), allocatable :: L_vgm_patches  (:) 
   REAL(r8), allocatable :: n_vgm_patches  (:) 
#endif
   REAL(r8), allocatable :: k_s_patches     (:) 
   REAL(r8), allocatable :: csol_patches    (:) 
   REAL(r8), allocatable :: tksatu_patches  (:) 
   REAL(r8), allocatable :: tkdry_patches   (:) 
   
   REAL(r8), allocatable :: theta_s_one (:) 
   REAL(r8), allocatable :: psi_s_one   (:) 
#ifdef Campbell_SOIL_MODEL
   REAL(r8), allocatable :: lambda_one  (:) 
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
   REAL(r8), allocatable :: theta_r_one  (:) 
   REAL(r8), allocatable :: alpha_vgm_one  (:) 
   REAL(r8), allocatable :: L_vgm_one  (:) 
   REAL(r8), allocatable :: n_vgm_one  (:) 
#endif
   REAL(r8), allocatable :: k_s_one     (:) 
   REAL(r8), allocatable :: csol_one    (:) 
   REAL(r8), allocatable :: tksatu_one  (:) 
   REAL(r8), allocatable :: tkdry_one   (:) 
   
   ! ........................................
   ! ... [2] aggregate the soil parameters from the resolution of raw data to modelling resolution
   ! ........................................
#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

   IF (p_is_master) THEN
      write(*,'(/, A29)') 'Aggregate Soil Parameters ...'
   ENDIF
      
   IF (p_is_worker) THEN

      allocate ( theta_s_patches(numpatch) )
      allocate ( psi_s_patches  (numpatch) )
#ifdef Campbell_SOIL_MODEL
      allocate ( lambda_patches (numpatch) )
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
      allocate ( theta_r_patches   (numpatch) )
      allocate ( alpha_vgm_patches (numpatch) )
      allocate ( L_vgm_patches     (numpatch) )
      allocate ( n_vgm_patches     (numpatch) )
#endif
      allocate ( k_s_patches    (numpatch) )
      allocate ( csol_patches   (numpatch) )
      allocate ( tksatu_patches (numpatch) )
      allocate ( tkdry_patches  (numpatch) )

   ENDIF

   DO nsl = 1, 8

      write(c,'(i1)') nsl

      ! (1) saturated water content [cm3/cm3]
      IF (p_is_io) THEN
        
         CALL allocate_block_data (gland, theta_s_grid)
         lndname = trim(dir_rawdata)//'/soil/theta_s.nc'
         CALL ncio_read_block (lndname, 'theta_s_l'//trim(c), gland, theta_s_grid)
#ifdef USEMPI
         CALL aggregation_lc_data_daemon (gland, theta_s_grid)
#endif
      ENDIF
      
      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch
            L = landpatch%ltyp(ipatch)

            IF (L /= 0) THEN
               CALL aggregation_lc_request_data (ipatch, gland, theta_s_grid, theta_s_one)
               theta_s_patches (ipatch) = median (theta_s_one, size(theta_s_one), spval)
            ELSE
               theta_s_patches (ipatch) = -1.0e36_r8
            ENDIF

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_lc_worker_done ()
#endif
      ENDIF

#ifdef CLMDEBUG
      CALL check_vector_data ('theta_s lev '//trim(c), theta_s_patches)
#endif

      lndname = trim(dir_model_landdata)//'/theta_s_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_pixelset_dimension (lndname, landpatch)
      CALL ncio_write_vector (lndname, 'theta_s_l'//trim(c)//'_patches', 'vector', landpatch, theta_s_patches, 1)

#ifdef vanGenuchten_Mualem_SOIL_MODEL
      ! (1-1) residual water content [cm3/cm3]
      IF (p_is_io) THEN
        
         CALL allocate_block_data (gland, theta_r_grid)
         lndname = trim(dir_rawdata)//'/soil/VGM_theta_r.nc'
         CALL ncio_read_block (lndname, 'VGM_theta_r_l'//trim(c), gland, theta_r_grid)
#ifdef USEMPI
         CALL aggregation_lc_data_daemon (gland, theta_r_grid)
#endif
      ENDIF
      
      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch
            L = landpatch%ltyp(ipatch)

            IF (L /= 0) THEN
               CALL aggregation_lc_request_data (ipatch, gland, theta_r_grid, theta_r_one)
               theta_r_patches (ipatch) = median (theta_r_one, size(theta_r_one), spval)
            ELSE
               theta_r_patches (ipatch) = -1.0e36_r8
            ENDIF

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_lc_worker_done ()
#endif
      ENDIF

#ifdef CLMDEBUG
      CALL check_vector_data ('theta_r lev '//trim(c), theta_r_patches)
#endif

      lndname = trim(dir_model_landdata)//'/theta_r_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_pixelset_dimension (lndname, landpatch)
      CALL ncio_write_vector (lndname, 'theta_r_l'//trim(c)//'_patches', 'vector', landpatch, theta_r_patches, 1)
#endif

      ! (2) matric potential at saturation [cm]
      IF (p_is_io) THEN
         CALL allocate_block_data (gland, psi_s_grid)
         lndname = trim(dir_rawdata)//'/soil/psi_s.nc'
         CALL ncio_read_block (lndname, 'psi_s_l'//trim(c), gland, psi_s_grid)
#ifdef USEMPI
         CALL aggregation_lc_data_daemon (gland, psi_s_grid)
#endif
      ENDIF
      
      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch
            L = landpatch%ltyp(ipatch)

            IF (L /= 0) THEN
               CALL aggregation_lc_request_data (ipatch, gland, psi_s_grid, psi_s_one)
               psi_s_patches (ipatch) = median (psi_s_one, size(psi_s_one), spval)
            ELSE
               psi_s_patches (ipatch) = -1.0e36_r8
            ENDIF

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_lc_worker_done ()
#endif
      ENDIF

#ifdef CLMDEBUG
      CALL check_vector_data ('psi_s lev '//trim(c), psi_s_patches)
#endif

      lndname = trim(dir_model_landdata)//'/psi_s_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_pixelset_dimension (lndname, landpatch)
      CALL ncio_write_vector (lndname, 'psi_s_l'//trim(c)//'_patches', 'vector', landpatch, psi_s_patches, 1)
      
#ifdef Campbell_SOIL_MODEL
      ! (3-1) pore size distribution index [dimensionless]
      IF (p_is_io) THEN
         CALL allocate_block_data (gland, lambda_grid)
         lndname = trim(dir_rawdata)//'/soil/lambda.nc'
         CALL ncio_read_block (lndname, 'lambda_l'//trim(c), gland, lambda_grid)
#ifdef USEMPI
         CALL aggregation_lc_data_daemon (gland, lambda_grid)
#endif
      ENDIF
      
      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch
            L = landpatch%ltyp(ipatch)

            IF (L /= 0) THEN
               CALL aggregation_lc_request_data (ipatch, gland, lambda_grid, lambda_one)
               lambda_patches (ipatch) = median (lambda_one, size(lambda_one), spval)
            ELSE
               lambda_patches (ipatch) = -1.0e36_r8
            ENDIF

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_lc_worker_done ()
#endif
      ENDIF

#ifdef CLMDEBUG
      CALL check_vector_data ('lambda lev '//trim(c), lambda_patches)
#endif

      lndname = trim(dir_model_landdata)//'/lambda_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_pixelset_dimension (lndname, landpatch)
      CALL ncio_write_vector (lndname, 'lambda_l'//trim(c)//'_patches', 'vector', landpatch, lambda_patches, 1)
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
      ! (3-2) alpha in VGM model
      IF (p_is_io) THEN
         CALL allocate_block_data (gland, alpha_vgm_grid)
         lndname = trim(dir_rawdata)//'/soil/VGM_alpha.nc'
         CALL ncio_read_block (lndname, 'VGM_alpha_l'//trim(c), gland, alpha_vgm_grid)
#ifdef USEMPI
         CALL aggregation_lc_data_daemon (gland, alpha_vgm_grid)
#endif
      ENDIF
      
      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch
            L = landpatch%ltyp(ipatch)

            IF (L /= 0) THEN
               CALL aggregation_lc_request_data (ipatch, gland, alpha_vgm_grid, alpha_vgm_one)
               alpha_vgm_patches (ipatch) = median (alpha_vgm_one, size(alpha_vgm_one), spval)
            ELSE
               alpha_vgm_patches (ipatch) = -1.0e36_r8
            ENDIF

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_lc_worker_done ()
#endif
      ENDIF

#ifdef CLMDEBUG
      CALL check_vector_data ('alpha VGM lev '//trim(c), alpha_vgm_patches)
#endif

      lndname = trim(dir_model_landdata)//'/alpha_vgm_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_pixelset_dimension (lndname, landpatch)
      CALL ncio_write_vector (lndname, 'alpha_vgm_l'//trim(c)//'_patches', 'vector', landpatch, alpha_vgm_patches, 1)
      
      ! (3-3) L in VGM model
      IF (p_is_io) THEN
         CALL allocate_block_data (gland, L_vgm_grid)
         lndname = trim(dir_rawdata)//'/soil/VGM_L.nc'
         CALL ncio_read_block (lndname, 'VGM_L_l'//trim(c), gland, L_vgm_grid)
#ifdef USEMPI
         CALL aggregation_lc_data_daemon (gland, L_vgm_grid)
#endif
      ENDIF
      
      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch
            L = landpatch%ltyp(ipatch)

            IF (L /= 0) THEN
               CALL aggregation_lc_request_data (ipatch, gland, L_vgm_grid, L_vgm_one)
               L_vgm_patches (ipatch) = median (L_vgm_one, size(L_vgm_one), spval)
            ELSE
               L_vgm_patches (ipatch) = -1.0e36_r8
            ENDIF

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_lc_worker_done ()
#endif
      ENDIF

#ifdef CLMDEBUG
      CALL check_vector_data ('L VGM lev '//trim(c), L_vgm_patches)
#endif

      lndname = trim(dir_model_landdata)//'/L_vgm_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_pixelset_dimension (lndname, landpatch)
      CALL ncio_write_vector (lndname, 'L_vgm_l'//trim(c)//'_patches', 'vector', landpatch, L_vgm_patches, 1)
      
      ! (3-4) n in VGM model
      IF (p_is_io) THEN
         CALL allocate_block_data (gland, n_vgm_grid)
         lndname = trim(dir_rawdata)//'/soil/VGM_n.nc'
         CALL ncio_read_block (lndname, 'VGM_n_l'//trim(c), gland, n_vgm_grid)
#ifdef USEMPI
         CALL aggregation_lc_data_daemon (gland, n_vgm_grid)
#endif
      ENDIF
      
      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch
            L = landpatch%ltyp(ipatch)

            IF (L /= 0) THEN
               CALL aggregation_lc_request_data (ipatch, gland, n_vgm_grid, n_vgm_one)
               n_vgm_patches (ipatch) = median (n_vgm_one, size(n_vgm_one), spval)
            ELSE
               n_vgm_patches (ipatch) = -1.0e36_r8
            ENDIF

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_lc_worker_done ()
#endif
      ENDIF

#ifdef CLMDEBUG
      CALL check_vector_data ('n VGM lev '//trim(c), n_vgm_patches)
#endif

      lndname = trim(dir_model_landdata)//'/n_vgm_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_pixelset_dimension (lndname, landpatch)
      CALL ncio_write_vector (lndname, 'n_vgm_l'//trim(c)//'_patches', 'vector', landpatch, n_vgm_patches, 1)
#endif


      ! (4) saturated hydraulic conductivity [cm/day]
      IF (p_is_io) THEN
         CALL allocate_block_data (gland, k_s_grid)
         lndname = trim(dir_rawdata)//'/soil/k_s.nc'
         CALL ncio_read_block (lndname, 'k_s_l'//trim(c), gland, k_s_grid)
#ifdef USEMPI
         CALL aggregation_lc_data_daemon (gland, k_s_grid)
#endif
      ENDIF
      
      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch
            L = landpatch%ltyp(ipatch)

            IF (L /= 0) THEN
               CALL aggregation_lc_request_data (ipatch, gland, k_s_grid, k_s_one)
               k_s_patches (ipatch) = median (k_s_one, size(k_s_one), spval)
            ELSE
               k_s_patches (ipatch) = -1.0e36_r8
            ENDIF

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_lc_worker_done ()
#endif
      ENDIF

#ifdef CLMDEBUG
      CALL check_vector_data ('k_s lev '//trim(c), k_s_patches)
#endif

      lndname = trim(dir_model_landdata)//'/k_s_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_pixelset_dimension (lndname, landpatch)
      CALL ncio_write_vector (lndname, 'k_s_l'//trim(c)//'_patches', 'vector', landpatch, k_s_patches, 1)

      ! (5) heat capacity of soil solids [J/(m3 K)]
      IF (p_is_io) THEN
         CALL allocate_block_data (gland, csol_grid)
         lndname = trim(dir_rawdata)//'/soil/csol.nc'
         CALL ncio_read_block (lndname, 'csol_l'//trim(c), gland, csol_grid)
#ifdef USEMPI
         CALL aggregation_lc_data_daemon (gland, csol_grid)
#endif
      ENDIF
      
      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch
            L = landpatch%ltyp(ipatch)

            IF (L /= 0) THEN
               CALL aggregation_lc_request_data (ipatch, gland, csol_grid, csol_one)
               csol_patches (ipatch) = median (csol_one, size(csol_one), spval)
            ELSE
               csol_patches (ipatch) = -1.0e36_r8
            ENDIF

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_lc_worker_done ()
#endif
      ENDIF

#ifdef CLMDEBUG
      CALL check_vector_data ('csol lev '//trim(c), csol_patches)
#endif

      lndname = trim(dir_model_landdata)//'/csol_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_pixelset_dimension (lndname, landpatch)
      CALL ncio_write_vector (lndname, 'csol_l'//trim(c)//'_patches', 'vector', landpatch, csol_patches, 1)

      ! (6) thermal conductivity of saturated soil [W/m-K]
      IF (p_is_io) THEN
         CALL allocate_block_data (gland, tksatu_grid)
         lndname = trim(dir_rawdata)//'/soil/tksatu.nc'
         CALL ncio_read_block (lndname, 'tksatu_l'//trim(c), gland, tksatu_grid)
#ifdef USEMPI
         CALL aggregation_lc_data_daemon (gland, tksatu_grid)
#endif
      ENDIF
      
      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch
            L = landpatch%ltyp(ipatch)

            IF (L /= 0) THEN
               CALL aggregation_lc_request_data (ipatch, gland, tksatu_grid, tksatu_one)
               tksatu_patches (ipatch) = median (tksatu_one, size(tksatu_one), spval)
            ELSE
               tksatu_patches (ipatch) = -1.0e36_r8
            ENDIF

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_lc_worker_done ()
#endif
      ENDIF

#ifdef CLMDEBUG
      CALL check_vector_data ('tksatu lev '//trim(c), tksatu_patches)
#endif

      lndname = trim(dir_model_landdata)//'/tksatu_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_pixelset_dimension (lndname, landpatch)
      CALL ncio_write_vector (lndname, 'tksatu_l'//trim(c)//'_patches', 'vector', landpatch, tksatu_patches, 1)

      ! (7) thermal conductivity for dry soil [W/(m-K)]
      IF (p_is_io) THEN
         CALL allocate_block_data (gland, tkdry_grid)
         lndname = trim(dir_rawdata)//'/soil/tkdry.nc'
         CALL ncio_read_block (lndname, 'tkdry_l'//trim(c), gland, tkdry_grid)
#ifdef USEMPI
         CALL aggregation_lc_data_daemon (gland, tkdry_grid)
#endif
      ENDIF
      
      IF (p_is_worker) THEN

         DO ipatch = 1, numpatch
            L = landpatch%ltyp(ipatch)

            IF (L /= 0) THEN
               CALL aggregation_lc_request_data (ipatch, gland, tkdry_grid, tkdry_one)
               tkdry_patches (ipatch) = median (tkdry_one, size(tkdry_one), spval)
            ELSE
               tkdry_patches (ipatch) = -1.0e36_r8
            ENDIF

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_lc_worker_done ()
#endif
      ENDIF

#ifdef CLMDEBUG
      CALL check_vector_data ('tkdry lev '//trim(c), tkdry_patches)
#endif

      lndname = trim(dir_model_landdata)//'/tkdry_l'//trim(c)//'_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_pixelset_dimension (lndname, landpatch)
      CALL ncio_write_vector (lndname, 'tkdry_l'//trim(c)//'_patches', 'vector', landpatch, tkdry_patches, 1)

   ENDDO


   ! Deallocate the allocatable array
   ! --------------------------------

   IF (p_is_worker) THEN

      deallocate ( theta_s_patches )
      deallocate ( psi_s_patches   )
#ifdef Campbell_SOIL_MODEL
      deallocate ( lambda_patches  )
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
      deallocate ( theta_r_patches  )
      deallocate ( alpha_vgm_patches)
      deallocate ( L_vgm_patches    )
      deallocate ( n_vgm_patches    )
#endif
      deallocate ( k_s_patches     )
      deallocate ( csol_patches    )
      deallocate ( tksatu_patches  )
      deallocate ( tkdry_patches   )

      IF (allocated(theta_s_one)) deallocate (theta_s_one)
      IF (allocated(psi_s_one  )) deallocate (psi_s_one  )
#ifdef Campbell_SOIL_MODEL
      IF (allocated(lambda_one )) deallocate (lambda_one )
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
      IF (allocated ( theta_r_one  )) deallocate ( theta_r_one  )
      IF (allocated ( alpha_vgm_one)) deallocate ( alpha_vgm_one)
      IF (allocated ( L_vgm_one    )) deallocate ( L_vgm_one    )
      IF (allocated ( n_vgm_one    )) deallocate ( n_vgm_one    )
#endif
      IF (allocated(k_s_one    )) deallocate (k_s_one    )
      IF (allocated(csol_one   )) deallocate (csol_one   )
      IF (allocated(tksatu_one )) deallocate (tksatu_one )
      IF (allocated(tkdry_one  )) deallocate (tkdry_one  )

   ENDIF

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

END SUBROUTINE aggregation_soil_parameters
!-----------------------------------------------------------------------
!EOP
