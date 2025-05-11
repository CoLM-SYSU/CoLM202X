#include <define.h>

SUBROUTINE Aggregation_ForestHeight ( &
      gland, dir_rawdata, dir_model_landdata, lc_year)

!-----------------------------------------------------------------------
!  Global Forest Height
!     (http://lidarradar.jpl.nasa.gov/)
!      Simard, M., N. Pinto, J. B. Fisher, and A. Baccini, 2011: Mapping
!      forest canopy height globally with spaceborne lidar.
!      J. Geophys. Res., 116, G04021.
!
!  Created by Yongjiu Dai, 02/2014
!
! !REVISIONS:
!  Hua Yuan,      ?/2020 : for land cover land use classifications
!  Shupeng Zhang, 01/2022: porting codes to MPI parallel version
!-----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_LandPatch
   USE MOD_NetCDFVector
   USE MOD_NetCDFBlock
#ifdef RangeCheck
   USE MOD_RangeCheck
#endif
   USE MOD_AggregationRequestData
   USE MOD_Utils

   USE MOD_Const_LC
   USE MOD_5x5DataReadin
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   USE MOD_LandPFT
#endif

#ifdef SrfdataDiag
   USE MOD_SrfdataDiag
#endif

   IMPLICIT NONE

   ! arguments:
   integer, intent(in) :: lc_year
   type(grid_type),  intent(in) :: gland
   character(len=*), intent(in) :: dir_rawdata
   character(len=*), intent(in) :: dir_model_landdata

   ! local variables:
   ! ---------------------------------------------------------------
   character(len=256) :: landdir, lndname, cyear
   integer :: L, ipatch, p

   type (block_data_real8_2d) :: tree_height
   real(r8), allocatable :: tree_height_patches(:), tree_height_one(:)

   ! for IGBP data
   character(len=256) :: dir_5x5, suffix
   type (block_data_real8_2d) :: htop
   type (block_data_real8_3d) :: pftPCT
   real(r8), allocatable :: htop_patches(:), htop_pfts(:), htop_pcs(:,:)
   real(r8), allocatable :: htop_one(:), area_one(:), pct_one(:,:)
   integer  :: ip, ipft
   real(r8) :: sumarea

#ifdef SrfdataDiag
   integer :: typpatch(N_land_classification+1), ityp
#ifndef CROP
   integer :: typpft  (N_PFT)
#else
   integer :: typpft  (N_PFT+N_CFT)
#endif
   integer :: typpc   (N_land_classification+1)
#endif

      write(cyear,'(i4.4)') lc_year
      landdir = trim(dir_model_landdata) // '/htop/' //trim(cyear)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         write(*,'(/, A)') 'Aggregate forest height ...'
         CALL system('mkdir -p ' // trim(adjustl(landdir)))
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef LULC_USGS
      lndname = trim(dir_rawdata)//'/Forest_Height.nc'

      IF (p_is_io) THEN
         CALL allocate_block_data (gland, tree_height)
         CALL ncio_read_block (lndname, 'forest_height', gland, tree_height)

#ifdef USEMPI
         CALL aggregation_data_daemon (gland, data_r8_2d_in1 = tree_height)
#endif
      ENDIF

      IF (p_is_worker) THEN

         allocate (tree_height_patches (numpatch))

         DO ipatch = 1, numpatch
            L = landpatch%settyp(ipatch)
            IF(L/=0 .and. L/=1 .and. L/=16 .and. L/=24)THEN
               ! NOT OCEAN(0)/URBAN and BUILT-UP(1)/WATER BODIES(16)/ICE(24)
               CALL aggregation_request_data (landpatch, ipatch, gland, zip = USE_zip_for_aggregation, &
                  data_r8_2d_in1 = tree_height, data_r8_2d_out1 = tree_height_one)
               tree_height_patches (ipatch) = median (tree_height_one, size(tree_height_one))
            ELSE
               tree_height_patches (ipatch) = -1.0e36_r8
            ENDIF
         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef RangeCheck
      CALL check_vector_data ('htop_patches ', tree_height_patches)
#endif

      lndname = trim(landdir)//'/htop_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'htop_patches', 'patch', landpatch, tree_height_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typpatch = (/(ityp, ityp = 0, N_land_classification)/)
      lndname  = trim(dir_model_landdata) // '/diag/htop_patch_' // trim(cyear) // '.nc'
      CALL srfdata_map_and_write (tree_height_patches, landpatch%settyp, typpatch, m_patch2diag, &
         -1.0e36_r8, lndname, 'htop', compress = 1, write_mode = 'one')
#endif

      IF (p_is_worker) THEN
         deallocate ( tree_height_patches )
      ENDIF
#endif


#ifdef LULC_IGBP
      IF (p_is_io) THEN
         CALL allocate_block_data (gland, htop)
      ENDIF

      IF (p_is_io) THEN
         dir_5x5 = trim(dir_rawdata) // '/plant_15s'
         suffix  = 'MOD'//trim(cyear)
         CALL read_5x5_data (dir_5x5, suffix, gland, 'HTOP', htop)
#ifdef USEMPI
         CALL aggregation_data_daemon (gland, data_r8_2d_in1 = htop)
#endif
      ENDIF

      IF (p_is_worker) THEN

         allocate (htop_patches (numpatch))

         DO ipatch = 1, numpatch

            IF (landpatch%settyp(ipatch) /= 0) THEN
               CALL aggregation_request_data (landpatch, ipatch, gland, zip = USE_zip_for_aggregation, &
                  area = area_one, data_r8_2d_in1 = htop, data_r8_2d_out1 = htop_one)
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

#ifdef RangeCheck
      CALL check_vector_data ('HTOP_patches ', htop_patches)
#endif

      lndname = trim(landdir)//'/htop_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'htop_patches', 'patch', landpatch, htop_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typpatch = (/(ityp, ityp = 0, N_land_classification)/)
      lndname  = trim(dir_model_landdata) // '/diag/htop_patch_' // trim(cyear) // '.nc'
      CALL srfdata_map_and_write (htop_patches, landpatch%settyp, typpatch, m_patch2diag, &
         -1.0e36_r8, lndname, 'htop', compress = 1, write_mode = 'one')
#endif

      IF (p_is_worker) THEN
         IF (allocated(htop_patches)) deallocate (htop_patches)
         IF (allocated(htop_one    )) deallocate (htop_one    )
         IF (allocated(area_one    )) deallocate (area_one    )
      ENDIF
#endif

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
      IF (p_is_io) THEN
         CALL allocate_block_data (gland, htop)
         CALL allocate_block_data (gland, pftPCT, N_PFT_modis, lb1 = 0)
      ENDIF

      dir_5x5 = trim(dir_rawdata) // '/plant_15s'
      suffix  = 'MOD'//trim(cyear)

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

            CALL aggregation_request_data (landpatch, ipatch, gland, zip = USE_zip_for_aggregation, &
               area = area_one, data_r8_2d_in1 = htop,   data_r8_2d_out1 = htop_one, &
               data_r8_3d_in1 = pftPCT, data_r8_3d_out1 = pct_one, n1_r8_3d_in1 = 16, lb1_r8_3d_in1 = 0)

            htop_patches(ipatch) = sum(htop_one * area_one) / sum(area_one)

#ifndef CROP
            IF (patchtypes(landpatch%settyp(ipatch)) == 0) THEN
#else
            IF (patchtypes(landpatch%settyp(ipatch)) == 0 .and. landpatch%settyp(ipatch)/=CROPLAND) THEN
#endif
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
            ELSEIF (landpatch%settyp(ipatch) == CROPLAND) THEN
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

#ifdef RangeCheck
      CALL check_vector_data ('HTOP_patches ', htop_patches)
      CALL check_vector_data ('HTOP_pfts    ', htop_pfts   )
#endif

      lndname = trim(landdir)//'/htop_patches.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'htop_patches', 'patch', landpatch, htop_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typpatch = (/(ityp, ityp = 0, N_land_classification)/)
      lndname  = trim(dir_model_landdata) // '/diag/htop_patch_' // trim(cyear) // '.nc'
      CALL srfdata_map_and_write (htop_patches, landpatch%settyp, typpatch, m_patch2diag, &
         -1.0e36_r8, lndname, 'htop', compress = 1, write_mode = 'one')
#endif

      lndname = trim(landdir)//'/htop_pfts.nc'
      CALL ncio_create_file_vector (lndname, landpft)
      CALL ncio_define_dimension_vector (lndname, landpft, 'pft')
      CALL ncio_write_vector (lndname, 'htop_pfts', 'pft', landpft, htop_pfts, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
#ifndef CROP
      typpft  = (/(ityp, ityp = 0, N_PFT-1)/)
#else
      typpft  = (/(ityp, ityp = 0, N_PFT+N_CFT-1)/)
#endif
      lndname = trim(dir_model_landdata) // '/diag/htop_pft_' // trim(cyear) // '.nc'
      CALL srfdata_map_and_write (htop_pfts, landpft%settyp, typpft, m_pft2diag, &
         -1.0e36_r8, lndname, 'htop', compress = 1, write_mode = 'one')
#endif

      IF (p_is_worker) THEN
         IF (allocated(htop_patches)) deallocate (htop_patches)
         IF (allocated(htop_pfts   )) deallocate (htop_pfts   )
         IF (allocated(htop_one    )) deallocate (htop_one    )
         IF (allocated(pct_one     )) deallocate (pct_one     )
         IF (allocated(area_one    )) deallocate (area_one    )
      ENDIF
#endif

END SUBROUTINE Aggregation_ForestHeight
