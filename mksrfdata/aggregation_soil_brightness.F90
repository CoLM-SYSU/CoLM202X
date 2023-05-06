#include <define.h>

SUBROUTINE aggregation_soil_brightness ( &
      gland, dir_rawdata, dir_model_landdata)
   ! ----------------------------------------------------------------------
   ! Creates land model surface dataset from original "raw" data files -
   !     data with 30 arc seconds resolution
   !
   ! Created by Yongjiu Dai, 03/2014
   ! ----------------------------------------------------------------------
   USE precision
   USE mod_namelist
   USE spmd_task
   USE mod_block
   USE mod_grid
   USE mod_landpatch
   USE ncio_block
   USE ncio_vector
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

   TYPE (block_data_int32_2d) :: isc 
   TYPE (block_data_real8_2d) :: a_s_v_refl
   TYPE (block_data_real8_2d) :: a_d_v_refl
   TYPE (block_data_real8_2d) :: a_s_n_refl
   TYPE (block_data_real8_2d) :: a_d_n_refl

   REAL(r8), allocatable :: soil_s_v_alb (:)
   REAL(r8), allocatable :: soil_d_v_alb (:)
   REAL(r8), allocatable :: soil_s_n_alb (:)
   REAL(r8), allocatable :: soil_d_n_alb (:)

   INTEGER :: ii, L
   INTEGER :: ipatch, iblkme, iblk, jblk, ix, iy
   REAL(r8), allocatable :: soil_one(:)

   ! ----------------------------------------------------------------------
   ! The soil color and reflectance is from the work:
   ! Peter J. Lawrence and Thomas N. Chase, 2007:
   ! Representing a MODIS consistent land surface in the Community Land Model (CLM 3.0):
   ! Part 1 generating MODIS consistent land surface parameters
   REAL(r8) soil_s_v_refl(20) ! Saturated visible soil reflectance
   REAL(r8) soil_d_v_refl(20) ! Dry visible soil reflectance
   REAL(r8) soil_s_n_refl(20) ! Saturated near infrared soil reflectance
   REAL(r8) soil_d_n_refl(20) ! Dry near infrared soil reflectance

   soil_s_v_refl = (/ 0.26, 0.24, 0.22, 0.20, 0.19, 0.18, 0.17, 0.16, 0.15, 0.14, &
      0.13, 0.12, 0.11, 0.10, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04 /)

   soil_d_v_refl = (/ 0.37, 0.35, 0.33, 0.31, 0.30, 0.29, 0.28, 0.27, 0.26, 0.25, &
      0.24, 0.23, 0.22, 0.21, 0.20, 0.19, 0.18, 0.17, 0.16, 0.15 /)

   soil_s_n_refl = (/ 0.52, 0.48, 0.44, 0.40, 0.38, 0.36, 0.34, 0.32, 0.30, 0.28, &
      0.26, 0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.10, 0.08 /)

   soil_d_n_refl = (/ 0.63, 0.59, 0.55, 0.51, 0.49, 0.47, 0.45, 0.43, 0.41, 0.39, &
      0.37, 0.35, 0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.19 /)

   landdir = trim(dir_model_landdata) // '/soil/'

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif
   IF (p_is_master) THEN
      write(*,'(/, A29)') 'Aggregate Soil Brightness ...'
      CALL system('mkdir -p ' // trim(adjustl(landdir)))
   ENDIF
#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef SinglePoint
   IF (USE_SITE_soilreflectance) THEN
      RETURN
   ENDIF
#endif

   ! ........................................
   ! ...  aggregate the soil parameters from the resolution of raw data to modelling resolution
   ! ........................................

   lndname = trim(dir_rawdata)//'/soil_brightness.nc'

   IF (p_is_io) THEN

      CALL allocate_block_data (gland, isc)
      CALL allocate_block_data (gland, a_s_v_refl)
      CALL allocate_block_data (gland, a_d_v_refl)
      CALL allocate_block_data (gland, a_s_n_refl)
      CALL allocate_block_data (gland, a_d_n_refl)
   
      ! Read in the index of soil brightness (color)
      CALL ncio_read_block (lndname, 'soil_brightness', gland, isc)

      DO iblkme = 1, gblock%nblkme 
         iblk = gblock%xblkme(iblkme)
         jblk = gblock%yblkme(iblkme)

         DO iy = 1, gland%ycnt(jblk)
            DO ix = 1, gland%xcnt(iblk)

               ii = isc%blk(iblk,jblk)%val(ix,iy)
               IF ((ii >= 1) .and. (ii <= 20)) THEN 
                  a_s_v_refl%blk(iblk,jblk)%val(ix,iy) = soil_s_v_refl( ii )
                  a_d_v_refl%blk(iblk,jblk)%val(ix,iy) = soil_d_v_refl( ii )
                  a_s_n_refl%blk(iblk,jblk)%val(ix,iy) = soil_s_n_refl( ii )
                  a_d_n_refl%blk(iblk,jblk)%val(ix,iy) = soil_d_n_refl( ii )
               ENDIF

            ENDDO
         ENDDO

      ENDDO

   ENDIF

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef USEMPI
   IF (p_is_io) THEN
      CALL aggregation_data_daemon (gland, data_r8_2d_in1 = a_s_v_refl)
   ENDIF
#endif

   IF (p_is_worker) THEN
      
      allocate ( soil_s_v_alb (numpatch) )

      DO ipatch = 1, numpatch
         L = landpatch%settyp(ipatch)
#ifdef USGS_CLASSIFICATION
         IF(L/=16 .and. L/=24)THEN  ! NOT OCEAN(0)/WATER BODIES(16)/GLACIER and ICESHEET(24)
#else
         IF(L/=17 .and. L/=15)THEN  ! NOT OCEAN(0)/WATER BODIES(17)/GLACIER and ICE SHEET(15)
#endif
            CALL aggregation_request_data (landpatch, ipatch, gland, &
               data_r8_2d_in1 = a_s_v_refl, data_r8_2d_out1 = soil_one)
            soil_s_v_alb (ipatch) = median (soil_one, size(soil_one))

         ELSE
            soil_s_v_alb (ipatch) = -1.0e36_r8 
         ENDIF

      ENDDO

#ifdef USEMPI
      CALL aggregation_worker_done ()
#endif
   ENDIF

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef USEMPI
   IF (p_is_io) THEN
      CALL aggregation_data_daemon (gland, data_r8_2d_in1 = a_d_v_refl)
   ENDIF
#endif

   IF (p_is_worker) THEN

      allocate ( soil_d_v_alb (numpatch) )

      DO ipatch = 1, numpatch
         L = landpatch%settyp(ipatch)
#ifdef USGS_CLASSIFICATION
         IF(L/=16 .and. L/=24)THEN  ! NOT OCEAN(0)/WATER BODIES(16)/GLACIER and ICESHEET(24)
#else
         IF(L/=17 .and. L/=15)THEN  ! NOT OCEAN(0)/WATER BODIES(17)/GLACIER and ICE SHEET(15)
#endif
            CALL aggregation_request_data (landpatch, ipatch, gland, &
               data_r8_2d_in1 = a_d_v_refl, data_r8_2d_out1 = soil_one)
            soil_d_v_alb (ipatch) = median (soil_one, size(soil_one))

         ELSE
            soil_d_v_alb (ipatch) = -1.0e36_r8 
         ENDIF

      ENDDO

#ifdef USEMPI
      CALL aggregation_worker_done ()
#endif
   ENDIF

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef USEMPI
   IF (p_is_io) THEN
      CALL aggregation_data_daemon (gland, data_r8_2d_in1 = a_s_n_refl)
   ENDIF
#endif

   IF (p_is_worker) THEN

      allocate ( soil_s_n_alb (numpatch) )

      DO ipatch = 1, numpatch
         L = landpatch%settyp(ipatch)
#ifdef USGS_CLASSIFICATION
         IF(L/=16 .and. L/=24)THEN  ! NOT OCEAN(0)/WATER BODIES(16)/GLACIER and ICESHEET(24)
#else
         IF(L/=17 .and. L/=15)THEN  ! NOT OCEAN(0)/WATER BODIES(17)/GLACIER and ICE SHEET(15)
#endif
            CALL aggregation_request_data (landpatch, ipatch, gland, &
               data_r8_2d_in1 = a_s_n_refl, data_r8_2d_out1 = soil_one)
            soil_s_n_alb (ipatch) = median (soil_one, size(soil_one))

         ELSE
            soil_s_n_alb (ipatch) = -1.0e36_r8 
         ENDIF

      ENDDO

#ifdef USEMPI
      CALL aggregation_worker_done ()
#endif
   ENDIF

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef USEMPI
   IF (p_is_io) THEN
      CALL aggregation_data_daemon (gland, data_r8_2d_in1 = a_d_n_refl)
   ENDIF
#endif

   IF (p_is_worker) THEN
      
      allocate ( soil_d_n_alb (numpatch) )

      DO ipatch = 1, numpatch
         L = landpatch%settyp(ipatch)
#ifdef USGS_CLASSIFICATION
         IF(L/=16 .and. L/=24)THEN  ! NOT OCEAN(0)/WATER BODIES(16)/GLACIER and ICESHEET(24)
#else
         IF(L/=17 .and. L/=15)THEN  ! NOT OCEAN(0)/WATER BODIES(17)/GLACIER and ICE SHEET(15)
#endif
            CALL aggregation_request_data (landpatch, ipatch, gland, &
               data_r8_2d_in1 = a_d_n_refl, data_r8_2d_out1 = soil_one)
            soil_d_n_alb (ipatch) = median (soil_one, size(soil_one))

         ELSE
            soil_d_n_alb (ipatch) = -1.0e36_r8 
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
   CALL check_vector_data ('s_v_alb ', soil_s_v_alb, -1.e36_r8)
   CALL check_vector_data ('d_v_alb ', soil_d_v_alb, -1.e36_r8)
   CALL check_vector_data ('s_n_alb ', soil_s_n_alb, -1.e36_r8)
   CALL check_vector_data ('d_n_alb ', soil_d_n_alb, -1.e36_r8)
#endif

#ifndef SinglePoint
   ! (1) Write-out the albedo of visible of the saturated soil
   lndname = trim(landdir)//'/soil_s_v_alb_patches.nc'
   CALL ncio_create_file_vector (lndname, landpatch)
   CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
   CALL ncio_write_vector (lndname, 'soil_s_v_alb', 'patch', landpatch, soil_s_v_alb, 1)

   ! (2) Write-out the albedo of visible of the dry soil
   lndname = trim(landdir)//'/soil_d_v_alb_patches.nc'
   CALL ncio_create_file_vector (lndname, landpatch)
   CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
   CALL ncio_write_vector (lndname, 'soil_d_v_alb', 'patch', landpatch, soil_d_v_alb, 1)

   ! (3) Write-out the albedo of near infrared of the saturated soil
   lndname = trim(landdir)//'/soil_s_n_alb_patches.nc'
   CALL ncio_create_file_vector (lndname, landpatch)
   CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
   CALL ncio_write_vector (lndname, 'soil_s_n_alb', 'patch', landpatch, soil_s_n_alb, 1)

   ! (4) Write-out the albedo of near infrared of the dry soil
   lndname = trim(landdir)//'/soil_d_n_alb_patches.nc'
   CALL ncio_create_file_vector (lndname, landpatch)
   CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
   CALL ncio_write_vector (lndname, 'soil_d_n_alb', 'patch', landpatch, soil_d_n_alb, 1)
#else
   SITE_soil_s_v_alb = soil_s_v_alb(1) 
   SITE_soil_d_v_alb = soil_d_v_alb(1)
   SITE_soil_s_n_alb = soil_s_n_alb(1)
   SITE_soil_d_n_alb = soil_d_n_alb(1)
#endif

   ! Deallocate the allocatable array
   ! --------------------------------
   IF (p_is_worker) THEN
      deallocate ( soil_s_v_alb )
      deallocate ( soil_d_v_alb )
      deallocate ( soil_s_n_alb )
      deallocate ( soil_d_n_alb )
   ENDIF

END SUBROUTINE aggregation_soil_brightness
!-----------------------------------------------------------------------
!EOP
