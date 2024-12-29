#include <define.h>

SUBROUTINE Aggregation_PercentagesPFT (gland, dir_rawdata, dir_model_landdata, lc_year)

! ----------------------------------------------------------------------
! Percentage of Plant Function Types
!
! Original from Hua Yuan's OpenMP version.
!
! REVISIONS:
! Hua Yuan,      ?/2020 : for land cover land use classifications
! Shupeng Zhang, 01/2022: porting codes to MPI parallel version
! ----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Vars_Global
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_LandPatch
#ifdef CROP
   USE MOD_LandCrop
#endif
   USE MOD_NetCDFBlock
   USE MOD_NetCDFVector
#ifdef RangeCheck
   USE MOD_RangeCheck
#endif
   USE MOD_AggregationRequestData

   USE MOD_Const_LC
   USE MOD_5x5DataReadin

#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   USE MOD_LandPFT
#endif
#ifdef SinglePoint
   USE MOD_SingleSrfdata
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
   ! ----------------------------------------------------------------------
   character(len=256) :: landdir, lndname

   ! for IGBP data
   character(len=256) :: dir_5x5, suffix, cyear
   ! for PFT
   type (block_data_real8_3d) :: pftPCT
   real(r8), allocatable :: pct_one(:), area_one(:)
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
   real(r8), allocatable :: pct_pft_one(:,:)
   real(r8), allocatable :: pct_pfts(:)
#endif
   integer  :: ipatch, ipc, ipft, p
   real(r8) :: sumarea
#ifdef SrfdataDiag
#ifdef CROP
   integer :: typcrop(N_CFT), ityp
   integer :: typpft(0:N_PFT+N_CFT-1)
#else
   integer :: typpft(0:N_PFT-1)
#endif
#endif

      write(cyear,'(i4.4)') lc_year
      landdir = trim(dir_model_landdata) // '/pctpft/' // trim(cyear)

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


#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)

#ifdef SinglePoint
      IF (USE_SITE_pctpfts) THEN
         RETURN
      ENDIF
#endif

      dir_5x5 = trim(dir_rawdata) // '/plant_15s'
      ! add parameter input for time year
      !write(cyear,'(i4.4)') lc_year
      suffix  = 'MOD'//trim(cyear)

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

#ifndef CROP
            IF (patchtypes(landpatch%settyp(ipatch)) == 0) THEN
#else
            IF (patchtypes(landpatch%settyp(ipatch)) == 0 .and. landpatch%settyp(ipatch)/=CROPLAND) THEN
#endif
               CALL aggregation_request_data (landpatch, ipatch, gland, zip = USE_zip_for_aggregation, area = area_one, &
                  data_r8_3d_in1 = pftPCT, data_r8_3d_out1 = pct_pft_one, n1_r8_3d_in1 = N_PFT_modis, lb1_r8_3d_in1 = 0)

#ifdef CROP
               pct_pft_one(N_PFT_modis-1,:) = 0.
#endif

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
            ELSEIF (landpatch%settyp(ipatch) == CROPLAND) THEN
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


#ifdef RangeCheck
      CALL check_vector_data ('PCT_PFTs ', pct_pfts)
#endif

#ifndef SinglePoint
      lndname = trim(landdir)//'/pct_pfts.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpft, 'pft')
      CALL ncio_write_vector (lndname, 'pct_pfts', 'pft', landpft, pct_pfts, DEF_Srfdata_CompressLevel)
#ifdef SrfdataDiag
#ifdef CROP
      typpft = (/(ipft, ipft = 0, N_PFT+N_CFT-1)/)
#else
      typpft = (/(ipft, ipft = 0, N_PFT-1)/)
#endif
      lndname = trim(dir_model_landdata)//'/diag/pct_pfts_'//trim(cyear)//'.nc'
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
      CALL ncio_write_vector (lndname, 'pct_crops', 'patch', landpatch, pctshrpch, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
      typcrop = (/(ityp, ityp = 1, N_CFT)/)
      lndname = trim(dir_model_landdata) // '/diag/pct_crop_patch_' // trim(cyear) // '.nc'
      CALL srfdata_map_and_write (pctshrpch, cropclass, typcrop, m_patch2diag, &
         -1.0e36_r8, lndname, 'pct_crop_patch', compress = 1, write_mode = 'one')
#endif
#else
      IF (.not. USE_SITE_pctcrop) THEN
         allocate (SITE_croptyp(numpatch))
         allocate (SITE_pctcrop(numpatch))
         SITE_croptyp = cropclass
         SITE_pctcrop = pctshrpch
      ENDIF
#endif
#endif

#endif

END SUBROUTINE Aggregation_PercentagesPFT
