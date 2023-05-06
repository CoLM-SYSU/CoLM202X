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
#ifdef CLMDEBUG 
   use mod_colm_debug
#endif
   use mod_aggregation
   USE mod_utils

   USE LC_Const
   USE mod_5x5_data
#ifdef PFT_CLASSIFICATION
   USE mod_landpft
#endif
#ifdef PC_CLASSIFICATION
   USE mod_landpc
#endif
#ifdef SinglePoint
   USE mod_single_srfdata
#endif

#ifdef SrfdataDiag
   USE mod_srfdata_diag
#endif

   IMPLICIT NONE

   ! arguments:
   type(grid_type),  intent(in) :: gland
   character(LEN=*), intent(in) :: dir_rawdata
   character(LEN=*), intent(in) :: dir_model_landdata

   ! local variables:
   ! ---------------------------------------------------------------
   character(len=256) :: landdir, lndname
   integer :: L, ipatch, p

   type (block_data_real8_2d) :: tree_height
   real(r8), allocatable :: tree_height_patches(:), tree_height_one(:)

   ! for IGBP data
   character(len=256) :: dir_5x5, suffix
   type (block_data_real8_2d) :: htop
   type (block_data_real8_3d) :: pftPCT
   real(r8), allocatable :: htop_patches(:), htop_pfts(:), htop_pcs(:,:)
   real(r8), allocatable :: htop_one(:), area_one(:), pct_one(:,:)
   INTEGER  :: ip, ipft
   REAL(r8) :: sumarea

   landdir = trim(dir_model_landdata) // '/htop/'

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif
   if (p_is_master) then
      write(*,'(/, A)') 'Aggregate forest height ...'
      CALL system('mkdir -p ' // trim(adjustl(landdir)))
   end if
#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef SinglePoint
   IF (USE_SITE_htop) THEN
      RETURN
   ENDIF
#endif

#ifdef USGS_CLASSIFICATION
   lndname = trim(dir_rawdata)//'/Forest_Height.nc' 

   if (p_is_io) then
      call allocate_block_data (gland, tree_height)
      call ncio_read_block (lndname, 'forest_height', gland, tree_height)

#ifdef USEMPI
      CALL aggregation_data_daemon (gland, data_r8_2d_in1 = tree_height)
#endif
   end if

   if (p_is_worker) then

      allocate (tree_height_patches (numpatch))
   
      do ipatch = 1, numpatch
         L = landpatch%settyp(ipatch)
         if(L/=0 .and. L/=1 .and. L/=16 .and. L/=24)then   
            ! NOT OCEAN(0)/URBAN and BUILT-UP(1)/WATER BODIES(16)/ICE(24)
            CALL aggregation_request_data (landpatch, ipatch, gland, &
               data_r8_2d_in1 = tree_height, data_r8_2d_out1 = tree_height_one)
            tree_height_patches (ipatch) = median (tree_height_one, size(tree_height_one))
         ELSE
            tree_height_patches (ipatch) = -1.0e36_r8
         ENDIF
      end do
      
#ifdef USEMPI
      CALL aggregation_worker_done ()
#endif
   end if

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CLMDEBUG
   call check_vector_data ('htop_patches ', tree_height_patches)
#endif 

#ifndef SinglePoint
   lndname = trim(landdir)//'/htop_patches.nc'
   CALL ncio_create_file_vector (lndname, landpatch)
   CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
   CALL ncio_write_vector (lndname, 'htop_patches', 'patch', landpatch, tree_height_patches, 1)
#else
   SITE_htop = tree_height_patches(1)
#endif 

   if (p_is_worker) then
      deallocate ( tree_height_patches )
   end if
#endif


#ifdef IGBP_CLASSIFICATION
   IF (p_is_io) THEN
      CALL allocate_block_data (gland, htop)
   ENDIF

   IF (p_is_io) THEN
      dir_5x5 = trim(dir_rawdata) // '/plant_15s_clim' 
      suffix  = 'MOD2005'
      CALL read_5x5_data (dir_5x5, suffix, gland, 'HTOP', htop)
#ifdef USEMPI
      CALL aggregation_data_daemon (gland, data_r8_2d_in1 = htop)
#endif
   ENDIF

   IF (p_is_worker) THEN

      allocate (htop_patches (numpatch))

      DO ipatch = 1, numpatch

         IF (landpatch%settyp(ipatch) /= 0) THEN
            CALL aggregation_request_data (landpatch, ipatch, gland, area = area_one, &
               data_r8_2d_in1 = htop, data_r8_2d_out1 = htop_one)
            htop_patches(ipatch) = sum(htop_one * area_one) / sum(area_one)
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
   CALL check_vector_data ('HTOP_patches ', htop_patches)
#endif

#ifndef SinglePoint
   lndname = trim(landdir)//'/htop_patches.nc'
   CALL ncio_create_file_vector (lndname, landpatch)
   CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
   CALL ncio_write_vector (lndname, 'htop_patches', 'patch', landpatch, htop_patches, 1)
#else
   SITE_htop = htop_patches(1)
#endif 

   IF (p_is_worker) THEN
      IF (allocated(htop_patches)) deallocate (htop_patches)
      IF (allocated(htop_one))     deallocate (htop_one)
      IF (allocated(area_one))     deallocate (area_one)
   ENDIF
#endif

#ifdef PFT_CLASSIFICATION
   IF (p_is_io) THEN
      CALL allocate_block_data (gland, htop)
      CALL allocate_block_data (gland, pftPCT, N_PFT_modis, lb1 = 0)
   ENDIF
     
   dir_5x5 = trim(dir_rawdata) // '/plant_15s_clim' 
   suffix  = 'MOD2005'
      
   IF (p_is_io) THEN
      CALL read_5x5_data     (dir_5x5, suffix, gland, 'HTOP',    htop  )
      CALL read_5x5_data_pft (dir_5x5, suffix, gland, 'PCT_PFT', pftPCT)
#ifdef USEMPI
      CALL aggregation_data_daemon (gland, &
         data_r8_2d_in1 = htop, data_r8_3d_in1 = pftPCT, n1_r8_3d_in1 = 16)
#endif
   ENDIF

   IF (p_is_worker) THEN
      
      allocate (htop_patches (numpatch))
      allocate (htop_pfts    (numpft  ))

      DO ipatch = 1, numpatch

         CALL aggregation_request_data (landpatch, ipatch, gland, area = area_one, &
            data_r8_2d_in1 = htop,   data_r8_2d_out1 = htop_one, &
            data_r8_3d_in1 = pftPCT, data_r8_3d_out1 = pct_one, n1_r8_3d_in1 = 16)

         htop_patches(ipatch) = sum(htop_one * area_one) / sum(area_one)

         IF (landpatch%settyp(ipatch) == 1) THEN
            DO ip = patch_pft_s(ipatch), patch_pft_e(ipatch)
               p = landpft%settyp(ip)
               sumarea = sum(pct_one(p,:) * area_one)
               IF (sumarea > 0) THEN
                  htop_pfts(ip) = sum(htop_one * pct_one(p,:) * area_one) / sumarea
               ELSE
                  htop_pfts(ip) = htop_patches(ipatch)
               ENDIF
            ENDDO
#ifdef CROP
         ELSEIF (landpatch%settyp(ipatch) == 12) THEN
            ip = patch_pft_s(ipatch)
            htop_pfts(ip) = htop_patches(ipatch)
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

#ifdef CLMDEBUG
   CALL check_vector_data ('HTOP_patches ', htop_patches)
   CALL check_vector_data ('HTOP_pfts    ', htop_pfts   )
#endif

#ifndef SinglePoint
   lndname = trim(landdir)//'/htop_patches.nc'
   CALL ncio_create_file_vector (lndname, landpatch)
   CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
   CALL ncio_write_vector (lndname, 'htop_patches', 'patch', landpatch, htop_patches, 1)
   
   lndname = trim(landdir)//'/htop_pfts.nc'
   CALL ncio_create_file_vector (lndname, landpft)
   CALL ncio_define_dimension_vector (lndname, landpft, 'pft')
   CALL ncio_write_vector (lndname, 'htop_pfts', 'pft', landpft, htop_pfts, 1)
#else
   allocate (SITE_htop_pfts(numpft))
   SITE_htop_pfts(:) = htop_pfts(:)
#endif
   
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
      CALL allocate_block_data (gland, pftPCT, N_PFT_modis, lb1 = 0)
   ENDIF
     
   dir_5x5 = trim(dir_rawdata) // '/plant_15s_clim' 
   suffix  = 'MOD2005'
      
   IF (p_is_io) THEN
      CALL read_5x5_data     (dir_5x5, suffix, gland, 'HTOP',    htop  )
      CALL read_5x5_data_pft (dir_5x5, suffix, gland, 'PCT_PFT', pftPCT)
#ifdef USEMPI
      CALL aggregation_data_daemon (gland, &
         data_r8_2d_in1 = htop, data_r8_3d_in1 = pftPCT, n1_r8_3d_in1 = 16)
#endif
   ENDIF

   IF (p_is_worker) THEN
      
      allocate (htop_patches (numpatch))
      allocate (htop_pcs (0:N_PFT-1, numpc))

      DO ipatch = 1, numpatch

         CALL aggregation_request_data (landpatch, ipatch, gland, area = area_one, &
            data_r8_2d_in1 = htop,   data_r8_2d_out1 = htop_one, &
            data_r8_3d_in1 = pftPCT, data_r8_3d_out1 = pct_one, n1_r8_3d_in1 = 16)

         htop_patches(ipatch) = sum(htop_one * area_one) / sum(area_one)

         IF (patchtypes(landpatch%settyp(ipatch)) == 0) THEN
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
   CALL aggregation_worker_done ()
#endif
   ENDIF

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CLMDEBUG
   CALL check_vector_data ('HTOP_patches ', htop_patches)
   CALL check_vector_data ('HTOP_pcs     ', htop_pcs    )
#endif

#ifndef SinglePoint
   lndname = trim(landdir)//'/htop_patches.nc'
   CALL ncio_create_file_vector (lndname, landpatch)
   CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
   CALL ncio_write_vector (lndname, 'htop_patches', 'patch', landpatch, htop_patches, 1)

   lndname = trim(landdir)//'/htop_pcs.nc'
   CALL ncio_create_file_vector (lndname, landpc)
   CALL ncio_define_dimension_vector (lndname, landpc, 'pc')
   CALL ncio_define_dimension_vector (lndname, landpc, 'pft', N_PFT)
   CALL ncio_write_vector (lndname, 'htop_pcs', 'pft', N_PFT, 'pc', landpc, htop_pcs, 1)
#else
   allocate (SITE_htop_pfts(N_PFT))
   SITE_htop_pfts(:) = htop_pcs(:,1)
#endif

   IF (p_is_worker) THEN
      IF (allocated(htop_patches)) deallocate (htop_patches)
      IF (allocated(htop_pcs    )) deallocate (htop_pcs    )
      IF (allocated(htop_one)) deallocate (htop_one)
      IF (allocated(pct_one )) deallocate (pct_one )
      IF (allocated(area_one)) deallocate (area_one)
   ENDIF
#endif


END SUBROUTINE aggregation_forest_height
