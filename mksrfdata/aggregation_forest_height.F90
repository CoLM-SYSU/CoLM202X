#include <define.h>

SUBROUTINE aggregation_forest_height ( &
      gland, dir_rawdata, dir_model_landdata)

   ! ----------------------------------------------------------------------
   ! 1. Global land cover types (updated with the specific dataset)
   !
   ! 2. Global Forest Height
   !    (http://lidarradar.jpl.nasa.gov/)
   !     Simard, M., N. Pinto, J. B. Fisher, and A. Baccini, 2011: Mapping
   !     forest canopy height globally with spaceborne lidar.
   !     J. Geophys. Res., 116, G04021.
   !
   ! Created by Yongjiu Dai, 02/2014
   ! ----------------------------------------------------------------------
   use precision
   use mod_namelist
   use spmd_task
   use mod_grid
   use mod_landpatch
   use ncio_vector
   use ncio_block
   use mod_colm_debug
   use mod_aggregation_lc
   USE mod_utils

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

   type(grid_type),  intent(in) :: gland
   character(LEN=*), intent(in) :: dir_rawdata
   character(LEN=*), intent(in) :: dir_model_landdata

   ! local variables:
   ! ---------------------------------------------------------------
   character(len=256) :: lndname
   integer :: L, ipatch, p

   type (block_data_real8_2d) :: tree_height
   real(r8), allocatable :: tree_height_patches(:), tree_height_one(:)

   ! for IGBP data
   character(len=256) :: dir_modis
   type (block_data_real8_2d) :: htop
   type (block_data_real8_3d) :: pftPCT
   real(r8), allocatable :: htop_patches(:), htop_pfts(:), htop_pcs(:,:)
   real(r8), allocatable :: htop_one(:), area_one(:), pct_one(:,:)
   INTEGER  :: ip, ipft
   REAL(r8) :: sumarea

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

   if (p_is_master) then
      write(*,'(/, A24)') 'Aggregate forest height ...'
   end if

#ifdef USGS_CLASSIFICATION
   lndname = trim(dir_rawdata)//'/Forest_Height.nc' 

   if (p_is_io) then
      call allocate_block_data (gland, tree_height)
      call ncio_read_block (lndname, 'forest_height', gland, tree_height)

#ifdef USEMPI
      CALL aggregation_lc_data_daemon (gland, tree_height)
#endif
   end if

   if (p_is_worker) then

      allocate (tree_height_patches (numpatch))
   
      do ipatch = 1, numpatch
         L = landpatch%ltyp(ipatch)
         if(L/=0 .and. L/=1 .and. L/=16 .and. L/=24)then   ! NOT OCEAN(0)/URBAN and BUILT-UP(1)/WATER            BODIES(16)/ICE(24)
            CALL aggregation_lc_request_data (ipatch, gland, tree_height, tree_height_one)
            tree_height_patches (ipatch) = median (tree_height_one, size(tree_height_one))
         ELSE
            tree_height_patches (ipatch) = -1.0e36_r8
         ENDIF
      end do
      
#ifdef USEMPI
      CALL aggregation_lc_worker_done ()
#endif
   end if

#ifdef CLMDEBUG
   call check_vector_data ('htop_patches ', tree_height_patches)
#endif 

   lndname = trim(dir_model_landdata)//'/htop_patches.nc'
   CALL ncio_create_file_vector (lndname, landpatch)
   CALL ncio_define_pixelset_dimension (lndname, landpatch)
   CALL ncio_write_vector (lndname, 'htop_patches', 'vector', landpatch, tree_height_patches, 1)

   if (p_is_worker) then
      deallocate ( tree_height_patches )
   end if
#endif


#ifdef IGBP_CLASSIFICATION
   IF (p_is_io) THEN
      CALL allocate_block_data (gland, htop)
   ENDIF

   IF (p_is_io) THEN
      dir_modis = trim(DEF_dir_rawdata) // '/srf_5x5' 
      CALL modis_read_data (dir_modis, 'HTOP', gland, htop)
#ifdef USEMPI
      CALL aggregation_lc_data_daemon (gland, htop)
#endif
   ENDIF

   IF (p_is_worker) THEN

      allocate (htop_patches (numpatch))

      DO ipatch = 1, numpatch

         IF (landpatch%ltyp(ipatch) /= 0) THEN
            CALL aggregation_lc_request_data (ipatch, gland, htop, htop_one, area_one)
            htop_patches(ipatch) = sum(htop_one * area_one) / sum(area_one)
         ENDIF

      ENDDO
      
#ifdef USEMPI
      CALL aggregation_lc_worker_done ()
#endif
   ENDIF

#ifdef CLMDEBUG
   CALL check_vector_data ('HTOP_patches ', htop_patches)
#endif

   lndname = trim(dir_model_landdata)//'/htop_patches.nc'
   CALL ncio_create_file_vector (lndname, landpatch)
   CALL ncio_define_pixelset_dimension (lndname, landpatch)
   CALL ncio_write_vector (lndname, 'htop_patches', 'vector', landpatch, htop_patches, 1)

   IF (p_is_worker) THEN
      IF (allocated(htop_patches)) deallocate (htop_patches)
      IF (allocated(htop_one))     deallocate (htop_one)
      IF (allocated(area_one))     deallocate (area_one)
   ENDIF
#endif

#ifdef PFT_CLASSIFICATION
   IF (p_is_io) THEN
      CALL allocate_block_data (gland, htop)
      CALL allocate_block_data (gland, pftPCT, N_PFT, lb1 = 0)
   ENDIF
     
   dir_modis = trim(DEF_dir_rawdata) // '/srf_5x5' 
      
   IF (p_is_io) THEN
      CALL modis_read_data     (dir_modis, 'HTOP',    gland, htop  )
      CALL modis_read_data_pft (dir_modis, 'PCT_PFT', gland, pftPCT)
#ifdef USEMPI
      CALL aggregation_pft_data_daemon (gland, pftPCT, data2 = htop)
#endif
   ENDIF

   IF (p_is_worker) THEN
      
      allocate (htop_patches (numpatch))
      allocate (htop_pfts    (numpft  ))

      DO ipatch = 1, numpatch

         CALL aggregation_pft_request_data (ipatch, gland, pftPCT, pct_one, &
            area = area_one, data2 = htop, dout2 = htop_one)

         htop_patches(ipatch) = sum(htop_one * area_one) / sum(area_one)

         IF (patchtypes(landpatch%ltyp(ipatch)) == 0) THEN
            DO ip = patch_pft_s(ipatch), patch_pft_e(ipatch)
               p = landpft%ltyp(ip)
               sumarea = sum(pct_one(p,:) * area_one)
               IF (sumarea > 0) THEN
                  htop_pfts(ip) = sum(htop_one * pct_one(p,:) * area_one) / sumarea
               ELSE
                  htop_pfts(ip) = htop_patches(ipatch)
               ENDIF
            ENDDO
         ENDIF
      ENDDO

#ifdef USEMPI
   CALL aggregation_pft_worker_done ()
#endif
   ENDIF

#ifdef CLMDEBUG
   CALL check_vector_data ('HTOP_patches ', htop_patches)
   CALL check_vector_data ('HTOP_pfts    ', htop_pfts   )
#endif

   lndname = trim(dir_model_landdata)//'/htop_patches.nc'
   CALL ncio_create_file_vector (lndname, landpatch)
   CALL ncio_define_pixelset_dimension (lndname, landpatch)
   CALL ncio_write_vector (lndname, 'htop_patches', 'vector', landpatch, htop_patches, 1)
   
   lndname = trim(dir_model_landdata)//'/htop_pfts.nc'
   CALL ncio_create_file_vector (lndname, landpft)
   CALL ncio_define_pixelset_dimension (lndname, landpft)
   CALL ncio_write_vector (lndname, 'htop_pfts', 'vector', landpft, htop_pfts, 1)
   
   IF (p_is_worker) THEN
      IF (allocated(htop_patches)) deallocate (htop_patches)
      IF (allocated(htop_pfts   )) deallocate (htop_pfts   )
      IF (allocated(htop_one)) deallocate (htop_one)
      IF (allocated(pct_one )) deallocate (pct_one )
      IF (allocated(area_one)) deallocate (area_one)
   ENDIF
#endif

#ifdef PC_CLASSIFICATION
   IF (p_is_io) THEN
      CALL allocate_block_data (gland, htop)
      CALL allocate_block_data (gland, pftPCT, N_PFT, lb1 = 0)
   ENDIF
     
   dir_modis = trim(DEF_dir_rawdata) // '/srf_5x5' 
      
   IF (p_is_io) THEN
      CALL modis_read_data     (dir_modis, 'HTOP',    gland, htop  )
      CALL modis_read_data_pft (dir_modis, 'PCT_PFT', gland, pftPCT)
#ifdef USEMPI
      CALL aggregation_pft_data_daemon (gland, pftPCT, data2 = htop)
#endif
   ENDIF

   IF (p_is_worker) THEN
      
      allocate (htop_patches (numpatch))
      allocate (htop_pcs (0:N_PFT-1, numpc))

      DO ipatch = 1, numpatch

         CALL aggregation_pft_request_data (ipatch, gland, pftPCT, pct_one, &
            area = area_one, data2 = htop, dout2 = htop_one)

         htop_patches(ipatch) = sum(htop_one * area_one) / sum(area_one)

         IF (patchtypes(landpatch%ltyp(ipatch)) == 0) THEN
            ip = patch2pc(ipatch)
            DO ipft = 0, N_PFT-1
               sumarea = sum(pct_one(ipft,:) * area_one)
               IF (sumarea > 0) THEN
                  htop_pcs(ipft,ip) = sum(htop_one * pct_one(ipft,:) * area_one) / sumarea
               ELSE
                  htop_pcs(ipft,ip) = htop_patches(ipatch)
               ENDIF
            ENDDO
         ENDIF
      ENDDO

#ifdef USEMPI
   CALL aggregation_pft_worker_done ()
#endif
   ENDIF

#ifdef CLMDEBUG
   CALL check_vector_data ('HTOP_patches ', htop_patches)
   CALL check_vector_data ('HTOP_pcs     ', htop_pcs    )
#endif

   lndname = trim(dir_model_landdata)//'/htop_patches.nc'
   CALL ncio_create_file_vector (lndname, landpatch)
   CALL ncio_define_pixelset_dimension (lndname, landpatch)
   CALL ncio_write_vector (lndname, 'htop_patches', 'vector', landpatch, htop_patches, 1)

   lndname = trim(dir_model_landdata)//'/htop_pcs.nc'
   CALL ncio_create_file_vector (lndname, landpc)
   CALL ncio_define_pixelset_dimension (lndname, landpc)
   CALL ncio_define_dimension_vector (lndname, 'pft', N_PFT)
   CALL ncio_write_vector (lndname, 'htop_pcs', 'pft', 'vector', landpc, N_PFT, htop_pcs, 1)
   
   IF (p_is_worker) THEN
      IF (allocated(htop_patches)) deallocate (htop_patches)
      IF (allocated(htop_pcs    )) deallocate (htop_pcs    )
      IF (allocated(htop_one)) deallocate (htop_one)
      IF (allocated(pct_one )) deallocate (pct_one )
      IF (allocated(area_one)) deallocate (area_one)
   ENDIF
#endif


END SUBROUTINE aggregation_forest_height
