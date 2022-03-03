#include <define.h>

SUBROUTINE aggregation_LAI (gland, dir_rawdata, dir_model_landdata)
   ! ----------------------------------------------------------------------
   ! 1. Global land cover types (updated with the specific dataset)
   !
   ! 2. Global Plant Leaf Area Index
   !    (http://globalchange.bnu.edu.cn)
   !    Yuan H., et al., 2011:
   !    Reprocessing the MODIS Leaf Area Index products for land surface 
   !    and climate modelling. Remote Sensing of Environment, 115: 1171-1187.
   !
   ! Created by Yongjiu Dai, 02/2014
   !
   ! ________________
   ! REVISION HISTORY:
   !   /07/2014, Siguang Zhu & Xiangxiang Zhang: weight average considering 
   !               partial overlap between fine grid and model grid for a user
   !               defined domain file.
   !
   ! ----------------------------------------------------------------------
   USE precision
   USE GlobalVars
   USE mod_namelist
   USE spmd_task
   USE mod_grid
   USE mod_landpatch
   USE ncio_block
   USE ncio_vector
   USE mod_colm_debug
   USE mod_aggregation_lc

   USE LC_Const
   USE mod_modis_data
#ifdef PFT_CLASSIFICATION
   USE mod_landpft
   USE mod_aggregation_pft
#endif
#ifdef PC_CLASSIFICATION
   USE mod_landpc
   USE mod_aggregation_pft
#endif

   IMPLICIT NONE

   ! arguments:

   TYPE(grid_type),  intent(in) :: gland
   CHARACTER(LEN=*), intent(in) :: dir_rawdata
   CHARACTER(LEN=*), intent(in) :: dir_model_landdata

   ! local variables:
   ! ----------------------------------------------------------------------
   CHARACTER(len=256) :: lndname

   TYPE (block_data_real8_2d) :: LAI          ! plant leaf area index (m2/m2)
   REAL(r8), allocatable :: LAI_patches(:), lai_one(:), area_one(:)
   INTEGER :: N8, Julian_day, ipatch
   CHARACTER(LEN=256) :: c

   ! for IGBP data
   CHARACTER(len=256) :: dir_modis
   INTEGER :: month
   TYPE (block_data_real8_2d) :: SAI          ! plant stem area index (m2/m2)
   REAL(r8), allocatable :: SAI_patches(:), sai_one(:)

   ! for PFT
   TYPE (block_data_real8_3d) :: pftLAI, pftSAI, pftPCT
   REAL(r8), allocatable :: pct_one (:), pct_pft_one(:,:)
   REAL(r8), allocatable :: LAI_pfts(:), lai_pft_one(:,:)
   REAL(r8), allocatable :: SAI_pfts(:), sai_pft_one(:,:)
   INTEGER :: p, ip

   REAL(r8), allocatable :: LAI_pcs(:,:), SAI_pcs(:,:)
   INTEGER :: ipc, ipft
   REAL(r8) :: sumarea
      
#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

   IF (p_is_master) THEN
      write(*,'(/, A17)') 'Aggregate LAI ...'
   ENDIF

   ! ................................................
   ! ... global plant leaf area index
   ! ................................................

#ifdef USGS_CLASSIFICATION
   IF (p_is_io) THEN
      CALL allocate_block_data (gland, LAI)
   ENDIF
   
   IF (p_is_worker) THEN
      allocate (LAI_patches (numpatch))
   ENDIF


   DO N8 = 1, 46

      ! -----------------------
      ! read in leaf area index
      ! -----------------------
      Julian_day = 1 + (N8-1)*8
      write(c,'(i3.3)') Julian_day

      IF (p_is_io) THEN
         lndname = trim(dir_rawdata)//'/lai/global_30s_10_year_avg/LAI_BNU_'//trim(c)//'.nc'
         CALL ncio_read_block (lndname, 'lai', gland, LAI)
         CALL block_data_linear_transform (LAI, scl = 0.1)

#ifdef USEMPI
         CALL aggregation_lc_data_daemon (gland, LAI)
#endif
      ENDIF


      ! ---------------------------------------------------------------
      ! aggregate the plant leaf area index from the resolution of raw data to modelling resolution
      ! ---------------------------------------------------------------

      IF (p_is_worker) THEN
         DO ipatch = 1, numpatch

            CALL aggregation_lc_request_data (ipatch, gland, LAI, lai_one, area_one)
            LAI_patches(ipatch) = sum(lai_one * area_one) / sum(area_one)

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_lc_worker_done ()
#endif
      ENDIF

#ifdef CLMDEBUG 
      CALL check_vector_data ('LAI_patches ' // trim(c), LAI_patches)
#endif

      ! ---------------------------------------------------
      ! write out the plant leaf area index of grid patches
      ! ---------------------------------------------------

      lndname = trim(dir_model_landdata)//'/LAI_patches'//trim(c)//'.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_pixelset_dimension (lndname, landpatch)
      CALL ncio_write_vector (lndname, 'LAI_patches', 'vector', landpatch, LAI_patches, 1)

   ENDDO

   IF (p_is_worker) THEN
      IF (allocated(LAI_patches)) deallocate(LAI_patches)
      IF (allocated(lai_one    )) deallocate(lai_one    )
      IF (allocated(area_one   )) deallocate(area_one   )
   ENDIF
#endif

#ifdef IGBP_CLASSIFICATION
   IF (p_is_io) THEN
      CALL allocate_block_data (gland, LAI)
   ENDIF
   
   IF (p_is_worker) THEN
      allocate (LAI_patches (numpatch))
   ENDIF

   DO month = 1, 12

      IF (p_is_io) THEN
         dir_modis = trim(DEF_dir_rawdata) // '/srf_5x5' 
         CALL modis_read_data_time (dir_modis, 'MONTHLY_LC_LAI', gland, month, LAI)
#ifdef USEMPI
         CALL aggregation_lc_data_daemon (gland, LAI)
#endif
      ENDIF

      ! ---------------------------------------------------------------
      ! aggregate the plant leaf area index from the resolution of raw data to modelling resolution
      ! ---------------------------------------------------------------

      IF (p_is_worker) THEN
         DO ipatch = 1, numpatch

            CALL aggregation_lc_request_data (ipatch, gland, LAI, lai_one, area_one)
            LAI_patches(ipatch) = sum(lai_one * area_one) / sum(area_one)

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_lc_worker_done ()
#endif
      ENDIF

#ifdef CLMDEBUG 
      write(c,'(i2.2)') month
      CALL check_vector_data ('LAI_patches ' // trim(c), LAI_patches)
#endif

      ! ---------------------------------------------------
      ! write out the plant leaf area index of grid patches
      ! ---------------------------------------------------

      write(c,'(i2.2)') month
      lndname = trim(dir_model_landdata)//'/LAI_patches'//trim(c)//'.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_pixelset_dimension (lndname, landpatch)
      CALL ncio_write_vector (lndname, 'LAI_patches', 'vector', landpatch, LAI_patches, 1)

   ENDDO

   IF (p_is_worker) THEN
      IF (allocated(LAI_patches)) deallocate(LAI_patches)
      IF (allocated(lai_one    )) deallocate(lai_one    )
      IF (allocated(area_one   )) deallocate(area_one   )
   ENDIF
   
   IF (p_is_io) THEN
      CALL allocate_block_data (gland, SAI)
   ENDIF
   IF (p_is_worker) THEN
      allocate (SAI_patches (numpatch))
   ENDIF
   
   DO month = 1, 12

      IF (p_is_io) THEN
         dir_modis = trim(DEF_dir_rawdata) // '/srf_5x5' 
         CALL modis_read_data_time (dir_modis, 'MONTHLY_LC_SAI', gland, month, SAI)
#ifdef USEMPI
         CALL aggregation_lc_data_daemon (gland, SAI)
#endif
      ENDIF

      ! ---------------------------------------------------------------
      ! aggregate the plant leaf area index from the resolution of raw data to modelling resolution
      ! ---------------------------------------------------------------

      IF (p_is_worker) THEN
         DO ipatch = 1, numpatch

            CALL aggregation_lc_request_data (ipatch, gland, SAI, sai_one, area_one)
            SAI_patches(ipatch) = sum(sai_one * area_one) / sum(area_one)

         ENDDO
      
#ifdef USEMPI
         CALL aggregation_lc_worker_done ()
#endif
      ENDIF


#ifdef CLMDEBUG 
      write(c,'(i2.2)') month
      CALL check_vector_data ('SAI_patches ' // trim(c), SAI_patches)
#endif

      ! ---------------------------------------------------
      ! write out the plant leaf area index of grid patches
      ! ---------------------------------------------------

      write(c,'(i2.2)') month
      lndname = trim(dir_model_landdata)//'/SAI_patches'//trim(c)//'.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_pixelset_dimension (lndname, landpatch)
      CALL ncio_write_vector (lndname, 'SAI_patches', 'vector', landpatch, SAI_patches, 1)

   ENDDO

   IF (p_is_worker) THEN
      IF (allocated(SAI_patches)) deallocate(SAI_patches)
      IF (allocated(sai_one    )) deallocate(sai_one    )
      IF (allocated(area_one   )) deallocate(area_one   )
   ENDIF
#endif

#ifdef PFT_CLASSIFICATION
   IF (p_is_io) THEN
      CALL allocate_block_data (gland, pftLAI, N_PFT, lb1 = 0)
      CALL allocate_block_data (gland, pftSAI, N_PFT, lb1 = 0)
      CALL allocate_block_data (gland, pftPCT, N_PFT, lb1 = 0)
   ENDIF
     
   dir_modis = trim(DEF_dir_rawdata) // '/srf_5x5' 
      
   IF (p_is_io) THEN
      CALL modis_read_data_pft (dir_modis, 'PCT_PFT', gland, pftPCT)
   ENDIF

   IF (p_is_worker) THEN
      allocate(LAI_patches (numpatch))
      allocate(LAI_pfts    (numpft  ))
   ENDIF
      
   DO month = 1, 12
      IF (p_is_io) THEN
         CALL modis_read_data_pft_time (dir_modis, 'MONTHLY_LAI', gland, month, pftLAI)
#ifdef USEMPI
         CALL aggregation_pft_data_daemon (gland, pftPCT, data3 = pftLAI)
#endif
      ENDIF

      ! ---------------------------------------------------------------
      ! aggregate the plant leaf area index from the resolution of raw data to modelling resolution
      ! ---------------------------------------------------------------

      IF (p_is_worker) THEN
         DO ipatch = 1, numpatch

            CALL aggregation_pft_request_data (ipatch, gland, pftPCT, pct_pft_one, &
               area = area_one, data3 = pftLAI, dout3 = lai_pft_one)
               
            IF (allocated(lai_one)) deallocate(lai_one)
            allocate(lai_one(size(area_one)))

            IF (allocated(pct_one)) deallocate(pct_one)
            allocate(pct_one(size(area_one)))

            pct_one = sum(pct_pft_one,dim=1)
            pct_one = max(pct_one, 1.0e-6)

            lai_one = sum(lai_pft_one * pct_pft_one, dim=1) / pct_one
            LAI_patches(ipatch) = sum(lai_one * area_one) / sum(area_one)
               
            IF (landpatch%ltyp(ipatch) == 1) THEN
               DO ip = patch_pft_s(ipatch), patch_pft_e(ipatch)
                  p = landpft%ltyp(ip)
                  sumarea = sum(pct_pft_one(p,:) * area_one)
                  IF (sumarea > 0) THEN
                     LAI_pfts(ip) = sum(lai_pft_one(p,:) * pct_pft_one(p,:) * area_one) / sumarea
                  ELSE
                     LAI_pfts(ip) = LAI_patches(ipatch)
                  ENDIF
               ENDDO
            ENDIF
         ENDDO

#ifdef USEMPI
         CALL aggregation_pft_worker_done ()
#endif
      ENDIF

      ! ---------------------------------------------------
      ! write out the plant leaf area index of grid patches
      ! ---------------------------------------------------
      write(c,'(i2.2)') month
#ifdef CLMDEBUG 
      CALL check_vector_data ('LAI_patches ' // trim(c), LAI_patches)
      CALL check_vector_data ('LAI_pfts    ' // trim(c), LAI_pfts   )
#endif

      lndname = trim(dir_model_landdata)//'/LAI_patches'//trim(c)//'.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_pixelset_dimension (lndname, landpatch)
      CALL ncio_write_vector (lndname, 'LAI_patches', 'vector', landpatch, LAI_patches, 1)
      
      lndname = trim(dir_model_landdata)//'/LAI_pfts'//trim(c)//'.nc'
      CALL ncio_create_file_vector (lndname, landpft)
      CALL ncio_define_pixelset_dimension (lndname, landpft)
      CALL ncio_write_vector (lndname, 'LAI_pfts', 'vector', landpft, LAI_pfts, 1)

   ENDDO

   IF (p_is_worker) THEN
      IF (allocated(LAI_patches)) deallocate(LAI_patches)
      IF (allocated(LAI_pfts   )) deallocate(LAI_pfts   )
      IF (allocated(lai_one    )) deallocate(lai_one    )
      IF (allocated(pct_one    )) deallocate(pct_one    )
      IF (allocated(pct_pft_one)) deallocate(pct_pft_one)
      IF (allocated(area_one   )) deallocate(area_one   )
   ENDIF
   
   IF (p_is_worker) THEN
      allocate(SAI_patches (numpatch))
      allocate(SAI_pfts    (numpft  ))
   ENDIF
   
   DO month = 1, 12
      IF (p_is_io) THEN
         CALL modis_read_data_pft_time (dir_modis, 'MONTHLY_SAI', gland, month, pftSAI)
#ifdef USEMPI
         CALL aggregation_pft_data_daemon (gland, pftPCT, data3 = pftSAI)
#endif
      ENDIF

      ! ---------------------------------------------------------------
      ! aggregate the plant leaf area index from the resolution of raw data to modelling resolution
      ! ---------------------------------------------------------------

      IF (p_is_worker) THEN
         DO ipatch = 1, numpatch

            CALL aggregation_pft_request_data (ipatch, gland, pftPCT, pct_pft_one, &
               area = area_one, data3 = pftSAI, dout3 = sai_pft_one)
               
            IF (allocated(sai_one)) deallocate(sai_one)
            allocate(sai_one(size(area_one)))
            
            IF (allocated(pct_one)) deallocate(pct_one)
            allocate(pct_one(size(area_one)))

            pct_one = sum(pct_pft_one,dim=1)
            pct_one = max(pct_one, 1.0e-6)

            sai_one = sum(sai_pft_one * pct_pft_one, dim=1) / pct_one
            SAI_patches(ipatch) = sum(sai_one * area_one) / sum(area_one)
               
            IF (landpatch%ltyp(ipatch) == 1) THEN
               DO ip = patch_pft_s(ipatch), patch_pft_e(ipatch)
                  p = landpft%ltyp(ip)
                  sumarea = sum(pct_pft_one(p,:) * area_one)
                  IF (sumarea > 0) THEN
                     SAI_pfts(ip) = sum(sai_pft_one(p,:) * pct_pft_one(p,:) * area_one) / sumarea
                  ELSE
                     SAI_pfts(ip) = SAI_patches(ipatch)
                  ENDIF
               ENDDO
            ENDIF
         ENDDO

#ifdef USEMPI
         CALL aggregation_pft_worker_done ()
#endif
      ENDIF

      ! ---------------------------------------------------
      ! write out the plant stem area index of grid patches
      ! ---------------------------------------------------
      write(c,'(i2.2)') month
#ifdef CLMDEBUG 
      CALL check_vector_data ('SAI_patches ' // trim(c), SAI_patches)
      CALL check_vector_data ('SAI_pfts    ' // trim(c), SAI_pfts   )
#endif
      
      lndname = trim(dir_model_landdata)//'/SAI_patches'//trim(c)//'.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_pixelset_dimension (lndname, landpatch)
      CALL ncio_write_vector (lndname, 'SAI_patches', 'vector', landpatch, SAI_patches, 1)
      
      lndname = trim(dir_model_landdata)//'/SAI_pfts'//trim(c)//'.nc'
      CALL ncio_create_file_vector (lndname, landpft)
      CALL ncio_define_pixelset_dimension (lndname, landpft)
      CALL ncio_write_vector (lndname, 'SAI_pfts', 'vector', landpft, SAI_pfts, 1)

   ENDDO

   IF (p_is_worker) THEN
      IF (allocated(SAI_patches)) deallocate(SAI_patches)
      IF (allocated(SAI_pfts   )) deallocate(SAI_pfts   )
      IF (allocated(sai_one    )) deallocate(sai_one    )
      IF (allocated(pct_one    )) deallocate(pct_one    )
      IF (allocated(pct_pft_one)) deallocate(pct_pft_one)
      IF (allocated(area_one   )) deallocate(area_one   )
   ENDIF
   
#endif

#ifdef PC_CLASSIFICATION
   IF (p_is_io) THEN
      CALL allocate_block_data (gland, pftLAI, N_PFT, lb1 = 0)
      CALL allocate_block_data (gland, pftSAI, N_PFT, lb1 = 0)
      CALL allocate_block_data (gland, pftPCT, N_PFT, lb1 = 0)
   ENDIF
     
   dir_modis = trim(DEF_dir_rawdata) // '/srf_5x5' 
      
   IF (p_is_io) THEN
      CALL modis_read_data_pft (dir_modis, 'PCT_PFT', gland, pftPCT)
   ENDIF

   IF (p_is_worker) THEN
      allocate(LAI_patches (numpatch))
      allocate(LAI_pcs (0:N_PFT-1, numpc))
   ENDIF
      
   DO month = 1, 12
      IF (p_is_io) THEN
         CALL modis_read_data_pft_time (dir_modis, 'MONTHLY_LAI', gland, month, pftLAI)
#ifdef USEMPI
         CALL aggregation_pft_data_daemon (gland, pftPCT, data3 = pftLAI)
#endif
      ENDIF

      ! ---------------------------------------------------------------
      ! aggregate the plant leaf area index from the resolution of raw data to modelling resolution
      ! ---------------------------------------------------------------

      IF (p_is_worker) THEN
         DO ipatch = 1, numpatch

            CALL aggregation_pft_request_data (ipatch, gland, pftPCT, pct_pft_one, &
               area = area_one, data3 = pftLAI, dout3 = lai_pft_one)
               
            IF (allocated(lai_one)) deallocate(lai_one)
            allocate(lai_one(size(area_one)))

            IF (allocated(pct_one)) deallocate(pct_one)
            allocate(pct_one(size(area_one)))

            pct_one = sum(pct_pft_one,dim=1)
            pct_one = max(pct_one, 1.0e-6)

            lai_one = sum(lai_pft_one * pct_pft_one, dim=1) / pct_one
            LAI_patches(ipatch) = sum(lai_one * area_one) / sum(area_one)
               
            IF (patchtypes(landpatch%ltyp(ipatch)) == 0) THEN
               ipc = patch2pc(ipatch)
               DO ipft = 0, N_PFT-1
                  sumarea = sum(pct_pft_one(ipft,:) * area_one)
                  IF (sumarea > 0) THEN
                     LAI_pcs(ipft,ipc) = sum(lai_pft_one(ipft,:) * pct_pft_one(ipft,:) * area_one) / sumarea
                  ELSE
                     LAI_pcs(ipft,ipc) = LAI_patches(ipatch)
                  ENDIF
               ENDDO 
            ENDIF
         ENDDO

#ifdef USEMPI
         CALL aggregation_pft_worker_done ()
#endif
      ENDIF

      ! ---------------------------------------------------
      ! write out the plant leaf area index of grid patches
      ! ---------------------------------------------------
      write(c,'(i2.2)') month
#ifdef CLMDEBUG 
      CALL check_vector_data ('LAI_patches ' // trim(c), LAI_patches)
      CALL check_vector_data ('LAI_pcs     ' // trim(c), LAI_pcs   )
#endif

      lndname = trim(dir_model_landdata)//'/LAI_patches'//trim(c)//'.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_pixelset_dimension (lndname, landpatch)
      CALL ncio_write_vector (lndname, 'LAI_patches', 'vector', landpatch, LAI_patches, 1)

      lndname = trim(dir_model_landdata)//'/LAI_pcs'//trim(c)//'.nc'
      CALL ncio_create_file_vector (lndname, landpc)
      CALL ncio_define_pixelset_dimension (lndname, landpc)
      CALL ncio_define_dimension_vector (lndname, 'pft', N_PFT)
      CALL ncio_write_vector (lndname, 'LAI_pcs', 'pft', 'vector', landpc, N_PFT, LAI_pcs, 1)

   ENDDO
   
   IF (p_is_worker) THEN
      IF (allocated(LAI_patches)) deallocate(LAI_patches)
      IF (allocated(LAI_pcs    )) deallocate(LAI_pcs    )
      IF (allocated(lai_one    )) deallocate(lai_one    )
      IF (allocated(pct_one    )) deallocate(pct_one    )
      IF (allocated(pct_pft_one)) deallocate(pct_pft_one)
      IF (allocated(area_one   )) deallocate(area_one   )
   ENDIF
   
   IF (p_is_worker) THEN
      allocate(SAI_patches (numpatch))
      allocate(SAI_pcs (0:N_PFT-1, numpc))
   ENDIF
      
   DO month = 1, 12
      IF (p_is_io) THEN
         CALL modis_read_data_pft_time (dir_modis, 'MONTHLY_SAI', gland, month, pftSAI)
#ifdef USEMPI
         CALL aggregation_pft_data_daemon (gland, pftPCT, data3 = pftSAI)
#endif
      ENDIF

      ! ---------------------------------------------------------------
      ! aggregate the plant leaf area index from the resolution of raw data to modelling resolution
      ! ---------------------------------------------------------------

      IF (p_is_worker) THEN
         DO ipatch = 1, numpatch

            CALL aggregation_pft_request_data (ipatch, gland, pftPCT, pct_pft_one, &
               area = area_one, data3 = pftSAI, dout3 = sai_pft_one)
               
            IF (allocated(sai_one)) deallocate(sai_one)
            allocate(sai_one(size(area_one)))

            IF (allocated(pct_one)) deallocate(pct_one)
            allocate(pct_one(size(area_one)))

            pct_one = sum(pct_pft_one,dim=1)
            pct_one = max(pct_one, 1.0e-6)

            sai_one = sum(sai_pft_one * pct_pft_one, dim=1) / pct_one
            SAI_patches(ipatch) = sum(sai_one * area_one) / sum(area_one)
               
            IF (patchtypes(landpatch%ltyp(ipatch)) == 0) THEN
               ipc = patch2pc(ipatch)
               DO ipft = 0, N_PFT-1
                  sumarea = sum(pct_pft_one(ipft,:) * area_one)
                  IF (sumarea > 0) THEN
                     SAI_pcs(ipft,ipc) = sum(sai_pft_one(ipft,:) * pct_pft_one(ipft,:) * area_one) / sumarea
                  ELSE
                     SAI_pcs(ipft,ipc) = SAI_patches(ipatch)
                  ENDIF
               ENDDO 
            ENDIF
         ENDDO

#ifdef USEMPI
         CALL aggregation_pft_worker_done ()
#endif
      ENDIF

      ! ---------------------------------------------------
      ! write out the plant stem area index of grid patches
      ! ---------------------------------------------------
      write(c,'(i2.2)') month
#ifdef CLMDEBUG 
      CALL check_vector_data ('SAI_patches ' // trim(c), SAI_patches)
      CALL check_vector_data ('SAI_pcs     ' // trim(c), SAI_pcs   )
#endif

      lndname = trim(dir_model_landdata)//'/SAI_patches'//trim(c)//'.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_pixelset_dimension (lndname, landpatch)
      CALL ncio_write_vector (lndname, 'SAI_patches', 'vector', landpatch, SAI_patches, 1)
      
      lndname = trim(dir_model_landdata)//'/SAI_pcs'//trim(c)//'.nc'
      CALL ncio_create_file_vector (lndname, landpc)
      CALL ncio_define_pixelset_dimension (lndname, landpc)
      CALL ncio_define_dimension_vector (lndname, 'pft', N_PFT)
      CALL ncio_write_vector (lndname, 'SAI_pcs', 'pft', 'vector', landpc, N_PFT, SAI_pcs, 1)

   ENDDO

   IF (p_is_worker) THEN
      IF (allocated(SAI_patches)) deallocate(SAI_patches)
      IF (allocated(SAI_pcs    )) deallocate(SAI_pcs    )
      IF (allocated(sai_one    )) deallocate(sai_one    )
      IF (allocated(pct_one    )) deallocate(pct_one    )
      IF (allocated(pct_pft_one)) deallocate(pct_pft_one)
      IF (allocated(area_one   )) deallocate(area_one   )
   ENDIF
   
   
#endif

END SUBROUTINE aggregation_LAI
