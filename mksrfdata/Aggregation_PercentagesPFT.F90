#include <define.h>

SUBROUTINE Aggregation_PercentagesPFT (gland, dir_rawdata, dir_model_landdata)

   USE precision
   USE MOD_Vars_Global
   USE mod_namelist
   USE spmd_task
   USE mod_grid
   USE MOD_LandPatch
   USE ncio_block
   USE ncio_vector
#ifdef CoLMDEBUG
   USE mod_colm_debug
#endif
   USE MOD_AggregationRequestData

   USE MOD_Const_LC
   USE mod_5x5_data
#ifdef PFT_CLASSIFICATION
   USE MOD_LandPFT
#endif
#ifdef PC_CLASSIFICATION
   USE MOD_LandPC
#endif
#ifdef SinglePoint
   USE MOD_SingleSrfdata
#endif

#ifdef SrfdataDiag
   USE MOD_SrfdataDiag
#endif

   IMPLICIT NONE

   ! arguments:

   TYPE(grid_type),  intent(in) :: gland
   CHARACTER(LEN=*), intent(in) :: dir_rawdata
   CHARACTER(LEN=*), intent(in) :: dir_model_landdata

   ! local variables:
   ! ----------------------------------------------------------------------
   CHARACTER(len=256) :: landdir, lndname

   ! for IGBP data
   CHARACTER(len=256) :: dir_5x5, suffix
   ! for PFT
   TYPE (block_data_real8_3d) :: pftPCT
   REAL(r8), allocatable :: pct_one(:), area_one(:)
#ifdef PFT_CLASSIFICATION
   REAL(r8), allocatable :: pct_pft_one(:,:)
   REAL(r8), allocatable :: pct_pfts(:)
#endif
#ifdef PC_CLASSIFICATION
   REAL(r8), allocatable :: pct_pft_one(:,:)
   REAL(r8), allocatable :: pct_pcs(:,:)
#endif
   INTEGER  :: ipatch, ipc, ipft, p
   REAL(r8) :: sumarea
#ifdef SrfdataDiag
#ifdef CROP
   INTEGER :: typcrop(N_CFT), ityp
   INTEGER :: typpft(0:N_PFT+N_CFT-1)
#else
   INTEGER :: typpft(0:N_PFT-1)
#endif
#endif

   landdir = trim(dir_model_landdata) // '/pctpft/'

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif
   IF (p_is_master) THEN
      write(*,'(/, A43)') 'Aggregate plant function type fractions ...'
      CALL system('mkdir -p ' // trim(adjustl(landdir)))
   ENDIF
#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif


#ifdef PFT_CLASSIFICATION

#ifdef SinglePoint
   IF (USE_SITE_pctpfts) THEN
      RETURN
   ENDIF
#endif

   dir_5x5 = trim(dir_rawdata) // '/plant_15s_clim'
   suffix  = 'MOD2005'

   IF (p_is_io) THEN
      CALL allocate_block_data (gland, pftPCT, N_PFT_modis, lb1 = 0)
      CALL read_5x5_data_pft   (dir_5x5, suffix, gland, 'PCT_PFT', pftPCT)
#ifdef USEMPI
      CALL aggregation_data_daemon (gland, data_r8_3d_in1 = pftPCT, n1_r8_3d_in1 = N_PFT_modis)
#endif
   ENDIF

   IF (p_is_worker) THEN

      allocate(pct_pfts (numpft))

      DO ipatch = 1, numpatch
         IF (landpatch%settyp(ipatch) == 1) THEN
            CALL aggregation_request_data (landpatch, ipatch, gland, area = area_one, &
               data_r8_3d_in1 = pftPCT, data_r8_3d_out1 = pct_pft_one, n1_r8_3d_in1 = N_PFT_modis, lb1_r8_3d_in1 = 0)

            pct_one = sum(pct_pft_one, dim=1)
            pct_one = max(pct_one, 1.0e-6)
            sumarea = sum(area_one)

            DO ipft = patch_pft_s(ipatch), patch_pft_e(ipatch)
               p = landpft%settyp(ipft)
               pct_pfts(ipft) = sum(pct_pft_one(p,:) / pct_one * area_one) / sumarea
            ENDDO

            pct_pfts(patch_pft_s(ipatch):patch_pft_e(ipatch)) =    &
               pct_pfts(patch_pft_s(ipatch):patch_pft_e(ipatch))   &
               / sum(pct_pfts(patch_pft_s(ipatch):patch_pft_e(ipatch)))
#ifdef CROP
         ELSEIF (landpatch%settyp(ipatch) == 12) THEN
            pct_pfts(patch_pft_s(ipatch):patch_pft_e(ipatch)) = 1.
#endif
         ENDIF
      ENDDO

#ifdef USEMPI
      CALL aggregation_worker_done ()
#endif
   ENDIF

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif


#ifdef CoLMDEBUG
   CALL check_vector_data ('PCT_PFTs ', pct_pfts)
#endif

#ifndef SinglePoint
   lndname = trim(landdir)//'/pct_pfts.nc'
   CALL ncio_create_file_vector (lndname, landpatch)
   CALL ncio_define_dimension_vector (lndname, landpft, 'pft')
   CALL ncio_write_vector (lndname, 'pct_pfts', 'pft', landpft, pct_pfts, 1)
#ifdef SrfdataDiag
#ifdef CROP
   typpft = (/(ipft, ipft = 0, N_PFT+N_CFT-1)/)
#else
   typpft = (/(ipft, ipft = 0, N_PFT-1)/)
#endif
   lndname = trim(dir_model_landdata)//'/diag/pct_pfts.nc'
   CALL srfdata_map_and_write (pct_pfts, landpft%settyp, typpft, m_pft2diag, &
      -1.0e36_r8, lndname, 'pctpfts', compress = 1, write_mode = 'one')
#endif
#else
   allocate (SITE_pctpfts(numpft))
   SITE_pctpfts = pct_pfts
#endif

   IF (p_is_worker) THEN
      IF (allocated(pct_pfts   )) deallocate(pct_pfts   )
      IF (allocated(pct_one    )) deallocate(pct_one    )
      IF (allocated(area_one   )) deallocate(area_one   )
      IF (allocated(pct_pft_one)) deallocate(pct_pft_one)
   ENDIF

#if (defined CROP)
#ifndef SinglePoint
   lndname = trim(landdir)//'/pct_crops.nc'
   CALL ncio_create_file_vector (lndname, landpatch)
   CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
   CALL ncio_write_vector (lndname, 'pct_crops', 'patch', landpatch, pctcrop, 1)

#ifdef SrfdataDiag
   typcrop = (/(ityp, ityp = 1, N_CFT)/)
   lndname = trim(dir_model_landdata) // '/diag/pct_crops_patch.nc'
   CALL srfdata_map_and_write (pctcrop, cropclass, typcrop, m_patch2diag, &
      -1.0e36_r8, lndname, 'pctcrop', compress = 1, write_mode = 'one')
#endif
#else
   allocate (SITE_croptyp(numpatch))
   allocate (SITE_pctcrop(numpatch))
   SITE_croptyp = cropclass
   SITE_pctcrop = pctcrop
#endif
#endif

#endif

#ifdef PC_CLASSIFICATION

#ifdef SinglePoint
   IF (USE_SITE_pctpfts) THEN
      RETURN
   ENDIF
#endif

   dir_5x5 = trim(dir_rawdata) // '/plant_15s_clim'
   suffix  = 'MOD2005'

   IF (p_is_io) THEN
      CALL allocate_block_data (gland, pftPCT, N_PFT_modis, lb1 = 0)
      CALL read_5x5_data_pft   (dir_5x5, suffix, gland, 'PCT_PFT', pftPCT)
#ifdef USEMPI
      CALL aggregation_data_daemon (gland, data_r8_3d_in1 = pftPCT, n1_r8_3d_in1 = N_PFT_modis)
#endif
   ENDIF

   IF (p_is_worker) THEN
      allocate(pct_pcs (0:N_PFT-1, numpc))

      DO ipatch = 1, numpatch

         IF (patchtypes(landpatch%settyp(ipatch)) == 0) THEN
            CALL aggregation_request_data (landpatch, ipatch, gland, area = area_one, &
               data_r8_3d_in1 = pftPCT, data_r8_3d_out1 = pct_pft_one, n1_r8_3d_in1 = N_PFT_modis, lb1_r8_3d_in1 = 0)

            pct_pft_one = max(pct_pft_one, 0.)

            pct_one = sum(pct_pft_one, dim=1)
            pct_one = max(pct_one, 1.0e-6)

            ipc = patch2pc(ipatch)
            DO ipft = 0, N_PFT-1
               sumarea = sum(area_one)
               pct_pcs(ipft,ipc) = sum(pct_pft_one(ipft,:) / pct_one * area_one) / sumarea
            ENDDO
         ENDIF
      ENDDO

#ifdef USEMPI
      CALL aggregation_worker_done ()
#endif
   ENDIF

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

   ! ---------------------------------------------------
   ! write out the plant leaf area index of grid patches
   ! ---------------------------------------------------
#ifdef CoLMDEBUG
   CALL check_vector_data ('PCT_PCs ', pct_pcs)
#endif

#ifndef SinglePoint
   lndname = trim(landdir)//'/pct_pcs.nc'
   CALL ncio_create_file_vector (lndname, landpatch)
   CALL ncio_define_dimension_vector (lndname, landpc, 'pc')
   CALL ncio_define_dimension_vector (lndname, landpc, 'pft', N_PFT)
   CALL ncio_write_vector (lndname, 'pct_pcs', 'pft', N_PFT, 'pc', landpc, pct_pcs, 1)
#else
   allocate (SITE_pctpfts (N_PFT))
   SITE_pctpfts = pct_pcs(:,1)
#endif

   IF (p_is_worker) THEN
      IF (allocated(pct_pcs    )) deallocate(pct_pcs    )
      IF (allocated(pct_one    )) deallocate(pct_one    )
      IF (allocated(area_one   )) deallocate(area_one   )
      IF (allocated(pct_pft_one)) deallocate(pct_pft_one)
   ENDIF

#endif

END SUBROUTINE Aggregation_PercentagesPFT