#include <define.h>

SUBROUTINE Aggregation_LAI (gridlai, dir_rawdata, dir_model_landdata, lc_year)
! ----------------------------------------------------------------------
! 1. Global Plant Leaf Area Index
!    (http://globalchange.bnu.edu.cn)
!    Yuan H., et al., 2011:
!    Reprocessing the MODIS Leaf Area Index products for land surface
!    and climate modelling. Remote Sensing of Environment, 115: 1171-1187.
!
! Created by Yongjiu Dai, 02/2014
!
! REVISIONS:
! Hua Yuan,      ?/2020 : for land cover land use classifications
! Shupeng Zhang, 01/2022: porting codes to MPI parallel version
! Hua Yuan,      05/2023: TODO
! ----------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Vars_Global
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_TimeManager
   USE MOD_Grid
   USE MOD_LandPatch
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
   type(grid_type),  intent(in) :: gridlai
   character(len=*), intent(in) :: dir_rawdata
   character(len=*), intent(in) :: dir_model_landdata

   ! local variables:
   ! ----------------------------------------------------------------------
   character(len=256) :: landdir, lndname

   integer :: simulation_lai_year_start, simulation_lai_year_end
   integer :: idate(3)

   type (block_data_real8_2d) :: LAI          ! plant leaf area index (m2/m2)
   real(r8), allocatable :: LAI_patches(:), lai_one(:), area_one(:)
   integer :: itime, ntime, Julian_day, ipatch
   character(len=4) :: c2, c3, cyear
   integer :: start_year, end_year, iy

   ! for IGBP data
   character(len=256) :: dir_5x5, suffix
   integer :: month
   type (block_data_real8_2d) :: SAI          ! plant stem area index (m2/m2)
   real(r8), allocatable :: SAI_patches(:), sai_one(:)

   ! for PFT and PC
   type (block_data_real8_3d) :: pftLSAI, pftPCT
   real(r8), allocatable :: pct_one (:), pct_pft_one(:,:)
   real(r8), allocatable :: LAI_pfts(:), lai_pft_one(:,:)
   real(r8), allocatable :: SAI_pfts(:), sai_pft_one(:,:)
   integer  :: p, ip
   real(r8) :: sumarea

#ifdef SrfdataDiag
   integer :: typpatch(N_land_classification+1), ityp
#ifndef CROP
   integer :: typpft  (N_PFT)
#else
   integer :: typpft  (N_PFT+N_CFT)
#endif
   character(len=256) :: varname
#endif

      ! LAI data root directory->case/landdata/LAI
      landdir = trim(dir_model_landdata) // '/LAI/'

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif
      IF (p_is_master) THEN
         write(*,'(/, A)') 'Aggregate LAI ...'
         CALL system('mkdir -p ' // trim(adjustl(landdir)))
      ENDIF
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef SinglePoint
      IF (USE_SITE_LAI) THEN
         RETURN
      ENDIF
#endif

      idate(1) = DEF_simulation_time%start_year
      IF (.not. isgreenwich) THEN
         idate(3) = DEF_simulation_time%start_sec
         CALL monthday2julian (idate(1), &
            DEF_simulation_time%start_month, DEF_simulation_time%start_day, idate(2))
         CALL localtime2gmt(idate)
      ENDIF

      simulation_lai_year_start = idate(1)

      idate(1) = DEF_simulation_time%end_year
      IF (.not. isgreenwich) THEN
         idate(3) = DEF_simulation_time%end_sec
         CALL monthday2julian (idate(1), &
            DEF_simulation_time%end_month, DEF_simulation_time%end_day, idate(2))
         CALL localtime2gmt(idate)
      ENDIF

      simulation_lai_year_end = idate(1)

! ................................................
! ... global plant leaf area index
! ................................................

#if (defined LULC_USGS || defined LULC_IGBP)
      ! add time variation of LAI
      IF (DEF_LAI_MONTHLY) THEN
         ! monthly average LAI
         ! if use lai change, LAI data of simulation start year and end year will be made
         ! if do not use lai change, only make LAI data of defined lc year
#ifdef LULCC
         ! 07/2023, NOTE: if defined LULCC, only one year (lc_year) lai processed.
         start_year = lc_year
         end_year   = lc_year
         ntime      = 12
#else
         IF (DEF_LAI_CHANGE_YEARLY) THEN
            start_year = max(simulation_lai_year_start, DEF_LAI_START_YEAR)
            end_year   = min(simulation_lai_year_end,   DEF_LAI_END_YEAR  )
            ntime      = 12
         ELSE
            start_year = lc_year
            end_year   = lc_year
            ntime      = 12
         ENDIF
#endif
      ! 8-day LAI
      ELSE
         start_year = max(simulation_lai_year_start, DEF_LAI_START_YEAR)
         end_year   = min(simulation_lai_year_end,   DEF_LAI_END_YEAR  )
         ntime      = 46
      ENDIF

      ! ----- LAI -----
      IF (p_is_io) THEN
         CALL allocate_block_data (gridlai, LAI)
      ENDIF

      IF (p_is_worker) THEN
         allocate (LAI_patches (numpatch))
      ENDIF

#ifdef SinglePoint

      allocate (SITE_LAI_year (start_year:end_year))
      SITE_LAI_year = (/(iy, iy = start_year, end_year)/)

      IF (DEF_LAI_MONTHLY) THEN
         !TODO-yuan-done: for multiple years
         allocate (SITE_LAI_monthly (12,start_year:end_year))
      ELSE
         allocate (SITE_LAI_8day    (46,start_year:end_year))
      ENDIF

#endif

      IF(.not. DEF_USE_LAIFEEDBACK)THEN
         DO iy = start_year, end_year

            write(cyear,'(i4.4)') iy
            CALL system('mkdir -p ' // trim(landdir) // trim(cyear))

            ! loop for month or 8-day
            DO itime = 1, ntime
            ! -----------------------
            ! read in leaf area index
            ! -----------------------
               IF (DEF_LAI_MONTHLY) THEN
                  write(c3, '(i2.2)') itime
               ELSE
                  Julian_day = 1 + (itime-1)*8
                  write(c3, '(i3.3)') Julian_day
               ENDIF

               IF (p_is_master) THEN
                  write(*,'(A,I4,A1,I3,A1,I3)') 'Aggregate LAI :', iy, ':', itime, '/', ntime
               ENDIF

               IF (p_is_io) THEN
                  IF (DEF_LAI_MONTHLY) THEN
                     dir_5x5 = trim(dir_rawdata) // '/plant_15s'
                     suffix  = 'MOD'//trim(cyear)
                     CALL read_5x5_data_time (dir_5x5, suffix, gridlai, 'MONTHLY_LC_LAI', itime, LAI)
                  ELSE
                     lndname = trim(dir_rawdata)//'/lai_15s_8day/lai_8-day_15s_'//trim(cyear)//'.nc'
                     CALL ncio_read_block_time (lndname, 'lai', gridlai, itime, LAI)
                     CALL block_data_linear_transform (LAI, scl = 0.1)
                  ENDIF

#ifdef USEMPI
                  CALL aggregation_data_daemon (gridlai, data_r8_2d_in1 = LAI)
#endif
               ENDIF

! -------------------------------------------------------------------------------------------
! aggregate the plant leaf area index from the resolution of raw data to modelling resolution
! -------------------------------------------------------------------------------------------

               IF (p_is_worker) THEN
                  DO ipatch = 1, numpatch
                     CALL aggregation_request_data (landpatch, ipatch, gridlai, zip = USE_zip_for_aggregation, &
                        area = area_one, data_r8_2d_in1 = LAI, data_r8_2d_out1 = lai_one)
                     LAI_patches(ipatch) = sum(lai_one * area_one) / sum(area_one)
                  ENDDO

#ifdef USEMPI
                  CALL aggregation_worker_done ()
#endif
               ENDIF

#ifdef RangeCheck
               CALL check_vector_data ('LAI value '//trim(c3), LAI_patches)
#endif

#ifdef USEMPI
               CALL mpi_barrier (p_comm_glb, p_err)
#endif
! ---------------------------------------------------
! write out the plant leaf area index of grid patches
! ---------------------------------------------------
#ifndef SinglePoint
               IF (DEF_LAI_MONTHLY) THEN
                  lndname = trim(landdir) // trim(cyear) // '/LAI_patches' // trim(c3) // '.nc'
               ELSE
                  !TODO: rename filename of 8-day LAI
                  lndname = trim(landdir) // trim(cyear) // '/LAI_patches' // trim(c3) // '.nc'
               ENDIF

               CALL ncio_create_file_vector (lndname, landpatch)
               CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
               CALL ncio_write_vector (lndname, 'LAI_patches', 'patch', landpatch, LAI_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
               typpatch = (/(ityp, ityp = 0, N_land_classification)/)
               lndname  = trim(dir_model_landdata) // '/diag/LAI_patch_'// trim(cyear) // '.nc'
               IF (DEF_LAI_MONTHLY) THEN
                  varname = 'LAI'
               ELSE
               !TODO: rename file name of 8-day LAI
                  varname = 'LAI_8-day'
               ENDIF
               CALL srfdata_map_and_write (LAI_patches, landpatch%settyp, typpatch, m_patch2diag, &
                  -1.0e36_r8, lndname, trim(varname), compress = 0, write_mode = 'one', &
                  lastdimname = 'Itime', lastdimvalue = itime)
#endif
#else
               ! single point cases
               !TODO: parameter input for time year
               IF (DEF_LAI_MONTHLY) THEN
                  SITE_LAI_monthly(itime,iy) = LAI_patches(1)
               ELSE
                  SITE_LAI_8day(itime,iy) = LAI_patches(1)
               ENDIF
#endif
            ENDDO
         ENDDO
      ENDIF

      ! ----- SAI -----
      IF (DEF_LAI_MONTHLY) THEN

         IF (p_is_io) THEN
            CALL allocate_block_data (gridlai, SAI)
         ENDIF

         IF (p_is_worker) THEN
            allocate (SAI_patches (numpatch))
         ENDIF

#ifdef SinglePoint
         allocate (SITE_SAI_monthly (12,start_year:end_year))
#endif

         dir_5x5 = trim(dir_rawdata) // '/plant_15s'
         DO iy = start_year, end_year
            write(cyear,'(i4.4)') iy
            suffix  = 'MOD'//trim(cyear)

            DO itime = 1, 12
               write(c3, '(i2.2)') itime

               IF (p_is_master) THEN
                  write(*,'(A,I4,A1,I3,A1,I3)') 'Aggregate SAI :', iy, ':', itime, '/', ntime
               ENDIF

               IF (p_is_io) THEN
                  CALL read_5x5_data_time (dir_5x5, suffix, gridlai, 'MONTHLY_LC_SAI', itime, SAI)

#ifdef USEMPI
                  CALL aggregation_data_daemon (gridlai, data_r8_2d_in1 = SAI)
#endif
               ENDIF

! -------------------------------------------------------------------------------------------
! aggregate the plant stem area index from the resolution of raw data to modelling resolution
! -------------------------------------------------------------------------------------------

               IF (p_is_worker) THEN
                  DO ipatch = 1, numpatch

                     CALL aggregation_request_data (landpatch, ipatch, gridlai, zip = USE_zip_for_aggregation, &
                        area = area_one, data_r8_2d_in1 = SAI, data_r8_2d_out1 = sai_one)
                     SAI_patches(ipatch) = sum(sai_one * area_one) / sum(area_one)

                  ENDDO

#ifdef USEMPI
                  CALL aggregation_worker_done ()
#endif
               ENDIF

#ifdef RangeCheck
            CALL check_vector_data ('SAI value '//trim(c3), SAI_patches)
#endif

#ifdef USEMPI
               CALL mpi_barrier (p_comm_glb, p_err)
#endif
! ---------------------------------------------------
! write out the plant leaf area index of grid patches
! ---------------------------------------------------
#ifndef SinglePoint
               lndname = trim(landdir) // trim(cyear) // '/SAI_patches' // trim(c3) // '.nc'
               CALL ncio_create_file_vector (lndname, landpatch)
               CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
               CALL ncio_write_vector (lndname, 'SAI_patches', 'patch', landpatch, SAI_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
               typpatch = (/(ityp, ityp = 0, N_land_classification)/)
               lndname  = trim(dir_model_landdata) // '/diag/SAI_patch_'// trim(cyear) // '.nc'
               IF (DEF_LAI_MONTHLY) THEN
                  varname = 'SAI'
               ELSE
                  !TODO: rename varname
                  varname = 'SAI_8-day'
               ENDIF
               CALL srfdata_map_and_write (SAI_patches, landpatch%settyp, typpatch, m_patch2diag, &
                  -1.0e36_r8, lndname, trim(varname), compress = 0, write_mode = 'one', &
                  lastdimname = 'Itime', lastdimvalue = itime)
#endif
#else
               !TODO: single point case
               SITE_SAI_monthly(itime,iy) = SAI_patches(1)
#endif
            ENDDO
         ENDDO
      ENDIF
#endif

! For both PFT and PC run LAI!!!!!
#if (defined LULC_IGBP_PFT || defined LULC_IGBP_PC)
         ! add time variation of LAI
         ! monthly average LAI
         ! if use lai change, LAI data of simulation start year and end year will be made
         ! if not use lai change, only make LAI data of defined lc year
#ifdef LULCC
         ! 07/2023, NOTE: if defined LULCC, only one year (lc_year) lai processed.
         start_year = lc_year
         end_year   = lc_year
         ntime      = 12
#else
      IF (DEF_LAI_CHANGE_YEARLY) THEN
         start_year = max(simulation_lai_year_start, DEF_LAI_START_YEAR)
         end_year   = min(simulation_lai_year_end,   DEF_LAI_END_YEAR  )
         ntime      = 12
      ELSE
         start_year = lc_year
         end_year   = lc_year
         ntime      = 12
      ENDIF
#endif

      IF (p_is_io) THEN
         CALL allocate_block_data (gridlai, pftLSAI, N_PFT_modis, lb1 = 0)
         CALL allocate_block_data (gridlai, pftPCT,  N_PFT_modis, lb1 = 0)
      ENDIF

      IF (p_is_worker) THEN
         allocate(LAI_patches (numpatch))
         allocate(LAI_pfts    (numpft  ))
         allocate(SAI_patches (numpatch))
         allocate(SAI_pfts    (numpft  ))
      ENDIF

#ifdef SinglePoint
      allocate (SITE_LAI_year (start_year:end_year))
      SITE_LAI_year = (/(iy, iy = start_year, end_year)/)

      !TODO-yuan-done: for multiple years
      allocate (SITE_LAI_pfts_monthly (numpft,12,start_year:end_year))
      allocate (SITE_SAI_pfts_monthly (numpft,12,start_year:end_year))
#endif

      dir_5x5 = trim(dir_rawdata) // '/plant_15s'
      DO iy = start_year, end_year
         write(cyear,'(i4.4)') iy
         suffix  = 'MOD'//trim(cyear)
         CALL system('mkdir -p ' // trim(landdir) // trim(cyear))

         IF (p_is_io) THEN
            CALL read_5x5_data_pft (dir_5x5, suffix, gridlai, 'PCT_PFT', pftPCT)
         ENDIF

         IF(.not. DEF_USE_LAIFEEDBACK)THEN
            DO month = 1, 12
               IF (p_is_io) THEN
                  CALL read_5x5_data_pft_time (dir_5x5, suffix, gridlai, 'MONTHLY_PFT_LAI', month, pftLSAI)
#ifdef USEMPI
                  CALL aggregation_data_daemon (gridlai, &
                     data_r8_3d_in1 = pftPCT,  n1_r8_3d_in1 = 16, &
                     data_r8_3d_in2 = pftLSAI, n1_r8_3d_in2 = 16)
#endif
               ENDIF

! -------------------------------------------------------------------------------------------
! aggregate the plant leaf area index from the resolution of raw data to modelling resolution
! -------------------------------------------------------------------------------------------

               IF (p_is_worker) THEN
                  DO ipatch = 1, numpatch
                     CALL aggregation_request_data (landpatch, ipatch, gridlai, zip = USE_zip_for_aggregation, area = area_one, &
                        data_r8_3d_in1 = pftPCT,  data_r8_3d_out1 = pct_pft_one, n1_r8_3d_in1 = 16, lb1_r8_3d_in1 = 0, &
                        data_r8_3d_in2 = pftLSAI, data_r8_3d_out2 = lai_pft_one, n1_r8_3d_in2 = 16, lb1_r8_3d_in2 = 0)

                     IF (allocated(lai_one)) deallocate(lai_one)
                     allocate(lai_one(size(area_one)))

                     IF (allocated(pct_one)) deallocate(pct_one)
                     allocate(pct_one(size(area_one)))

                     pct_one = sum(pct_pft_one,dim=1)
                     pct_one = max(pct_one, 1.0e-6)

                     lai_one = sum(lai_pft_one * pct_pft_one, dim=1) / pct_one
                     LAI_patches(ipatch) = sum(lai_one * area_one) / sum(area_one)

#ifndef CROP
                     IF (patchtypes(landpatch%settyp(ipatch)) == 0) THEN
#else
                     IF (patchtypes(landpatch%settyp(ipatch)) == 0 .and. landpatch%settyp(ipatch)/=CROPLAND) THEN
#endif
                        DO ip = patch_pft_s(ipatch), patch_pft_e(ipatch)
                           p = landpft%settyp(ip)
                           sumarea = sum(pct_pft_one(p,:) * area_one)
                           IF (sumarea > 0) THEN
                              LAI_pfts(ip) = sum(lai_pft_one(p,:) * pct_pft_one(p,:) * area_one) / sumarea
                           ELSE
                           ! 07/2023, yuan: bug may exist below
                           !LAI_pfts(ip) = LAI_patches(ipatch)
                              LAI_pfts(ip) = 0.
                           ENDIF
                        ENDDO
#ifdef CROP
                     ELSEIF (landpatch%settyp(ipatch) == CROPLAND) THEN
                        ip = patch_pft_s(ipatch)
                        LAI_pfts(ip) = LAI_patches(ipatch)
#endif
                     ENDIF
                  ENDDO

#ifdef USEMPI
                  CALL aggregation_worker_done ()
#endif
               ENDIF

               write(c2,'(i2.2)') month
#ifdef RangeCheck
               CALL check_vector_data ('LAI_patches ' // trim(c2), LAI_patches)
               CALL check_vector_data ('LAI_pfts    ' // trim(c2), LAI_pfts   )
#endif
#ifdef USEMPI
               CALL mpi_barrier (p_comm_glb, p_err)
#endif

! ---------------------------------------------------
! write out the plant leaf area index of grid patches
! ---------------------------------------------------
#ifndef SinglePoint
               lndname = trim(landdir)//trim(cyear)//'/LAI_patches'//trim(c2)//'.nc'
               CALL ncio_create_file_vector (lndname, landpatch)
               CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
               CALL ncio_write_vector (lndname, 'LAI_patches', 'patch', landpatch, LAI_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
               typpatch = (/(ityp, ityp = 0, N_land_classification)/)
               lndname  = trim(dir_model_landdata) // '/diag/LAI_patch_'// trim(cyear) // '.nc'
               varname  = 'LAI'
               CALL srfdata_map_and_write (LAI_patches, landpatch%settyp, typpatch, m_patch2diag, &
                  -1.0e36_r8, lndname, trim(varname), compress = 0, write_mode = 'one', &
                  lastdimname = 'Itime', lastdimvalue = month)
#endif

               lndname = trim(landdir)//trim(cyear)//'/LAI_pfts'//trim(c2)//'.nc'
               CALL ncio_create_file_vector (lndname, landpft)
               CALL ncio_define_dimension_vector (lndname, landpft, 'pft')
               CALL ncio_write_vector (lndname, 'LAI_pfts', 'pft', landpft, LAI_pfts, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
#ifndef CROP
               typpft  = (/(ityp, ityp = 0, N_PFT-1)/)
#else
               typpft  = (/(ityp, ityp = 0, N_PFT+N_CFT-1)/)
#endif
               lndname = trim(dir_model_landdata) // '/diag/LAI_pft_'// trim(cyear) // '.nc'
               varname = 'LAI_pft'
               CALL srfdata_map_and_write (LAI_pfts, landpft%settyp, typpft, m_pft2diag, &
                  -1.0e36_r8, lndname, trim(varname), compress = 0, write_mode = 'one',  &
                  lastdimname = 'Itime', lastdimvalue = month)
#endif
#else
               !TODO: single point case
               SITE_LAI_pfts_monthly(:,month,iy) = LAI_pfts(:)
#endif
            ! loop end of month
            ENDDO

         ENDIF

         DO month = 1, 12
            IF (p_is_io) THEN
               CALL read_5x5_data_pft_time (dir_5x5, suffix, gridlai, 'MONTHLY_PFT_SAI', month, pftLSAI)
#ifdef USEMPI
               CALL aggregation_data_daemon (gridlai, &
                  data_r8_3d_in1 = pftPCT,  n1_r8_3d_in1 = 16, &
                  data_r8_3d_in2 = pftLSAI, n1_r8_3d_in2 = 16)
#endif
            ENDIF

! -------------------------------------------------------------------------------------------
! aggregate the plant leaf area index from the resolution of raw data to modelling resolution
! -------------------------------------------------------------------------------------------

            IF (p_is_worker) THEN
               DO ipatch = 1, numpatch

                  CALL aggregation_request_data (landpatch, ipatch, gridlai, zip = USE_zip_for_aggregation, area = area_one, &
                     data_r8_3d_in1 = pftPCT,  data_r8_3d_out1 = pct_pft_one, n1_r8_3d_in1 = 16, lb1_r8_3d_in1 = 0, &
                     data_r8_3d_in2 = pftLSAI, data_r8_3d_out2 = sai_pft_one, n1_r8_3d_in2 = 16, lb1_r8_3d_in2 = 0)

                  IF (allocated(sai_one)) deallocate(sai_one)
                  allocate(sai_one(size(area_one)))

                  IF (allocated(pct_one)) deallocate(pct_one)
                  allocate(pct_one(size(area_one)))

                  pct_one = sum(pct_pft_one,dim=1)
                  pct_one = max(pct_one, 1.0e-6)

                  sai_one = sum(sai_pft_one * pct_pft_one, dim=1) / pct_one
                  SAI_patches(ipatch) = sum(sai_one * area_one) / sum(area_one)

#ifndef CROP
                  IF (patchtypes(landpatch%settyp(ipatch)) == 0) THEN
#else
                  IF (patchtypes(landpatch%settyp(ipatch)) == 0 .and. landpatch%settyp(ipatch)/=CROPLAND) THEN
#endif
                     DO ip = patch_pft_s(ipatch), patch_pft_e(ipatch)
                        p = landpft%settyp(ip)
                        sumarea = sum(pct_pft_one(p,:) * area_one)
                        IF (sumarea > 0) THEN
                           SAI_pfts(ip) = sum(sai_pft_one(p,:) * pct_pft_one(p,:) * area_one) / sumarea
                        ELSE
                           ! 07/2023, yuan: bug may exist below
                           !SAI_pfts(ip) = SAI_patches(ipatch)
                           SAI_pfts(ip) = 0.
                        ENDIF
                     ENDDO
#ifdef CROP
                  ELSEIF (landpatch%settyp(ipatch) == CROPLAND) THEN
                     ip = patch_pft_s(ipatch)
                     SAI_pfts(ip) = SAI_patches(ipatch)
#endif
                  ENDIF
               ENDDO

#ifdef USEMPI
               CALL aggregation_worker_done ()
#endif
            ENDIF

         write(c2,'(i2.2)') month
#ifdef RangeCheck
         CALL check_vector_data ('SAI_patches ' // trim(c2), SAI_patches)
         CALL check_vector_data ('SAI_pfts    ' // trim(c2), SAI_pfts   )
#endif
#ifdef USEMPI
            CALL mpi_barrier (p_comm_glb, p_err)
#endif

! ---------------------------------------------------
! write out the plant stem area index of grid patches
! ---------------------------------------------------
#ifndef SinglePoint
            lndname = trim(landdir)//trim(cyear)//'/SAI_patches'//trim(c2)//'.nc'
            CALL ncio_create_file_vector (lndname, landpatch)
            CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
            CALL ncio_write_vector (lndname, 'SAI_patches', 'patch', landpatch, SAI_patches, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
            typpatch = (/(ityp, ityp = 0, N_land_classification)/)
            lndname  = trim(dir_model_landdata) // '/diag/SAI_patch_'// trim(cyear) // '.nc'
            varname  = 'SAI'
            CALL srfdata_map_and_write (SAI_patches, landpatch%settyp, typpatch, m_patch2diag, &
               -1.0e36_r8, lndname, trim(varname), compress = 0, write_mode = 'one', &
               lastdimname = 'Itime', lastdimvalue = month)
#endif

            lndname = trim(landdir)//trim(cyear)//'/SAI_pfts'//trim(c2)//'.nc'
            CALL ncio_create_file_vector (lndname, landpft)
            CALL ncio_define_dimension_vector (lndname, landpft, 'pft')
            CALL ncio_write_vector (lndname, 'SAI_pfts', 'pft', landpft, SAI_pfts, DEF_Srfdata_CompressLevel)

#ifdef SrfdataDiag
#ifndef CROP
            typpft  = (/(ityp, ityp = 0, N_PFT-1)/)
#else
            typpft  = (/(ityp, ityp = 0, N_PFT+N_CFT-1)/)
#endif
            lndname = trim(dir_model_landdata) // '/diag/SAI_pft_'// trim(cyear) // '.nc'
            varname = 'SAI_pft'
            CALL srfdata_map_and_write (SAI_pfts, landpft%settyp, typpft, m_pft2diag, &
               -1.0e36_r8, lndname, trim(varname), compress = 0, write_mode = 'one',  &
               lastdimname = 'Itime', lastdimvalue = month)
#endif
#else
            SITE_SAI_pfts_monthly(:,month,iy) = SAI_pfts(:)
#endif
         ! loop end of month
         ENDDO
      ! loop end of year
      ENDDO

      IF (p_is_worker) THEN

         IF (allocated(LAI_patches)) deallocate(LAI_patches)
         IF (allocated(LAI_pfts   )) deallocate(LAI_pfts   )
         IF (allocated(lai_one    )) deallocate(lai_one    )

         IF (allocated(SAI_patches)) deallocate(SAI_patches)
         IF (allocated(SAI_pfts   )) deallocate(SAI_pfts   )
         IF (allocated(sai_one    )) deallocate(sai_one    )
         IF (allocated(pct_one    )) deallocate(pct_one    )
         IF (allocated(pct_pft_one)) deallocate(pct_pft_one)
         IF (allocated(area_one   )) deallocate(area_one   )

      ENDIF
#endif

END SUBROUTINE Aggregation_LAI
