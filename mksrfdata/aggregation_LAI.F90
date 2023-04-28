#include <define.h>

SUBROUTINE aggregation_LAI (gridlai, dir_rawdata, dir_model_landdata)
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
#ifdef CLMDEBUG 
   USE mod_colm_debug
#endif
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
#ifdef SinglePoint
   USE mod_single_srfdata
#endif

   IMPLICIT NONE

   ! arguments:

   TYPE(grid_type),  intent(in) :: gridlai
   CHARACTER(LEN=*), intent(in) :: dir_rawdata
   CHARACTER(LEN=*), intent(in) :: dir_model_landdata

   ! local variables:
   ! ----------------------------------------------------------------------
   CHARACTER(len=256) :: landdir, lndname

   TYPE (block_data_real8_2d) :: LAI          ! plant leaf area index (m2/m2)
   REAL(r8), allocatable :: LAI_patches(:), lai_one(:), area_one(:)
   INTEGER :: itime, ntime, Julian_day, ipatch
   CHARACTER(LEN=4) :: c2, c3, cyear
   integer :: start_year, end_year, YY   

   ! for IGBP data
   CHARACTER(len=256) :: dir_modis
   INTEGER :: month
   TYPE (block_data_real8_2d) :: SAI          ! plant stem area index (m2/m2)
   REAL(r8), allocatable :: SAI_patches(:), sai_one(:)

   ! for PFT
   TYPE (block_data_real8_3d) :: pftLSAI, pftPCT
   REAL(r8), allocatable :: pct_one (:), pct_pft_one(:,:)
   REAL(r8), allocatable :: LAI_pfts(:), lai_pft_one(:,:)
   REAL(r8), allocatable :: SAI_pfts(:), sai_pft_one(:,:)
   INTEGER :: p, ip

   ! for PC
   REAL(r8), allocatable :: LAI_pcs(:,:), SAI_pcs(:,:)
   INTEGER :: ipc, ipft
   REAL(r8) :: sumarea
      
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

   ! ................................................
   ! ... global plant leaf area index
   ! ................................................

#if (defined USGS_CLASSIFICATION || defined IGBP_CLASSIFICATION)
   IF (DEF_LAI_CLIM) THEN
      start_year = 1
      end_year   = 1
      ntime = 12 
   ELSE
      start_year = DEF_simulation_time%start_year
      end_year   = DEF_simulation_time%end_year
      ntime = 46
   ENDIF

   ! ----- LAI -----
   IF (p_is_io) THEN
      CALL allocate_block_data (gridlai, LAI)
   ENDIF
   
   IF (p_is_worker) THEN
      allocate (LAI_patches (numpatch))
   ENDIF

#ifdef SinglePoint
   IF (DEF_LAI_CLIM) THEN
      allocate (SITE_LAI_clim (12))
   ELSE
      allocate (SITE_LAI_year (start_year:end_year))
      SITE_LAI_year = (/(YY, YY = start_year, end_year)/)

      allocate (SITE_LAI_modis (46,start_year:end_year))
   ENDIF
#endif

   DO YY = start_year, end_year

      IF (.not. DEF_LAI_CLIM) THEN
         write(cyear,'(i4.4)') YY
         CALL system('mkdir -p ' // trim(landdir) // '/' // trim(cyear))
      ENDIF

      DO itime = 1, ntime
         ! -----------------------
         ! read in leaf area index
         ! -----------------------
         IF (DEF_LAI_CLIM) THEN
            write(c3, '(i2.2)') itime
         ELSE
            Julian_day = 1 + (itime-1)*8
            write(c3, '(i3.3)') Julian_day
         ENDIF

         IF (p_is_master) THEN
            write(*,'(A,I4,A1,I3,A1,I3)') 'Aggregate LAI :', YY, ':', itime, '/', ntime 
         endif

         IF (p_is_io) THEN
            IF (DEF_LAI_CLIM) THEN
               dir_modis = trim(DEF_dir_rawdata) // '/plant_15s_clim' 
               CALL modis_read_data_time (dir_modis, 'MONTHLY_LC_LAI', gridlai, itime, LAI)
            ELSE
               lndname = trim(dir_rawdata)//'/lai_15s_8day/lai_8-day_15s_'//trim(cyear)//'.nc'
               CALL ncio_read_block_time (lndname, 'lai', gridlai, itime, LAI)
               CALL block_data_linear_transform (LAI, scl = 0.1)
            ENDIF

#ifdef USEMPI
            CALL aggregation_lc_data_daemon (gridlai, LAI)
#endif
         ENDIF


         ! ---------------------------------------------------------------
         ! aggregate the plant leaf area index from the resolution of raw data to modelling resolution
         ! ---------------------------------------------------------------

         IF (p_is_worker) THEN
            DO ipatch = 1, numpatch
               CALL aggregation_lc_request_data (ipatch, gridlai, LAI, lai_one, area_one)
               LAI_patches(ipatch) = sum(lai_one * area_one) / sum(area_one)
            ENDDO

#ifdef USEMPI
            CALL aggregation_lc_worker_done ()
#endif
         ENDIF

#ifdef CLMDEBUG 
         CALL check_vector_data ('LAI value '//trim(c3), LAI_patches)
#endif

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif
         ! ---------------------------------------------------
         ! write out the plant leaf area index of grid patches
         ! ---------------------------------------------------
#ifndef SinglePoint
         IF (DEF_LAI_CLIM) THEN
            lndname = trim(landdir) // '/LAI_patches' // trim(c3) // '.nc'
         ELSE
            lndname = trim(landdir) // '/' // trim(cyear) // '/LAI_patches' // trim(c3) // '.nc'
         ENDIF

         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'LAI_patches', 'patch', landpatch, LAI_patches, 1)
#else
         IF (DEF_LAI_CLIM) THEN
            SITE_LAI_clim(itime) = LAI_patches(1)
         ELSE
            SITE_LAI_modis(itime,YY) = LAI_patches(1)
         ENDIF
#endif
      ENDDO
   ENDDO
   
   ! ----- SAI -----
   IF (DEF_LAI_CLIM) THEN
     
      IF (p_is_io) THEN
         CALL allocate_block_data (gridlai, SAI)
      ENDIF
      IF (p_is_worker) THEN
         allocate (SAI_patches (numpatch))
      ENDIF
#ifdef SinglePoint
      allocate (SITE_SAI_clim (12))
#endif

      dir_modis = trim(DEF_dir_rawdata) // '/plant_15s_clim' 

      DO itime = 1, 12 
         write(c3, '(i2.2)') itime

         IF (p_is_master) THEN
            write(*,'(A,I3,A1,I3)') 'Aggregate SAI :', itime, '/', itime
         endif

         IF (p_is_io) THEN
            CALL modis_read_data_time (dir_modis, 'MONTHLY_LC_SAI', gridlai, itime, SAI)

#ifdef USEMPI
            CALL aggregation_lc_data_daemon (gridlai, SAI)
#endif
         ENDIF

         ! ---------------------------------------------------------------
         ! aggregate the plant stem area index from the resolution of raw data to modelling resolution
         ! ---------------------------------------------------------------

         IF (p_is_worker) THEN
            DO ipatch = 1, numpatch

               CALL aggregation_lc_request_data (ipatch, gridlai, SAI, sai_one, area_one)
               SAI_patches(ipatch) = sum(sai_one * area_one) / sum(area_one)

            ENDDO
      
#ifdef USEMPI
            CALL aggregation_lc_worker_done ()
#endif
         ENDIF

#ifdef CLMDEBUG 
         CALL check_vector_data ('SAI value '//trim(c3), SAI_patches)
#endif

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif
         ! ---------------------------------------------------
         ! write out the plant leaf area index of grid patches
         ! ---------------------------------------------------
#ifndef SinglePoint
         lndname = trim(landdir) // '/SAI_patches' // trim(c3) // '.nc'
         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'SAI_patches', 'patch', landpatch, SAI_patches, 1)
#else
         SITE_SAI_clim(itime) = SAI_patches(1)
#endif
      ENDDO
   ENDIF
#endif

#ifdef PFT_CLASSIFICATION
   IF (p_is_io) THEN
      CALL allocate_block_data (gridlai, pftLSAI, N_PFT_modis, lb1 = 0)
      CALL allocate_block_data (gridlai, pftPCT,  N_PFT_modis, lb1 = 0)
   ENDIF
     
   dir_modis = trim(DEF_dir_rawdata) // '/plant_15s_clim' 
      
   IF (p_is_io) THEN
      CALL modis_read_data_pft (dir_modis, 'PCT_PFT', gridlai, pftPCT)
   ENDIF

   IF (p_is_worker) THEN
      allocate(LAI_patches (numpatch))
      allocate(LAI_pfts    (numpft  ))
   ENDIF
      
#ifdef SinglePoint
   allocate (SITE_LAI_pfts_clim (numpft,12))  
#endif

   DO month = 1, 12
      IF (p_is_io) THEN
         CALL modis_read_data_pft_time (dir_modis, 'MONTHLY_LAI', gridlai, month, pftLSAI)
#ifdef USEMPI
         CALL aggregation_pft_data_daemon (gridlai, pftPCT, data3 = pftLSAI)
#endif
      ENDIF

      ! ---------------------------------------------------------------
      ! aggregate the plant leaf area index from the resolution of raw data to modelling resolution
      ! ---------------------------------------------------------------

      IF (p_is_worker) THEN
         DO ipatch = 1, numpatch

            CALL aggregation_pft_request_data (ipatch, gridlai, pftPCT, pct_pft_one, &
               area = area_one, data3 = pftLSAI, dout3 = lai_pft_one)
               
            IF (allocated(lai_one)) deallocate(lai_one)
            allocate(lai_one(size(area_one)))

            IF (allocated(pct_one)) deallocate(pct_one)
            allocate(pct_one(size(area_one)))

            pct_one = sum(pct_pft_one,dim=1)
            pct_one = max(pct_one, 1.0e-6)

            lai_one = sum(lai_pft_one * pct_pft_one, dim=1) / pct_one
            LAI_patches(ipatch) = sum(lai_one * area_one) / sum(area_one)
               
            IF (landpatch%settyp(ipatch) == 1) THEN
               DO ip = patch_pft_s(ipatch), patch_pft_e(ipatch)
                  p = landpft%settyp(ip)
                  sumarea = sum(pct_pft_one(p,:) * area_one)
                  IF (sumarea > 0) THEN
                     LAI_pfts(ip) = sum(lai_pft_one(p,:) * pct_pft_one(p,:) * area_one) / sumarea
                  ELSE
                     LAI_pfts(ip) = LAI_patches(ipatch)
                  ENDIF
               ENDDO
#ifdef CROP
            ELSEIF (landpatch%settyp(ipatch) == 12) THEN
               ip = patch_pft_s(ipatch)
               LAI_pfts(ip) = LAI_patches(ipatch)
#endif
            ENDIF
         ENDDO

#ifdef USEMPI
         CALL aggregation_pft_worker_done ()
#endif
      ENDIF

      write(c2,'(i2.2)') month
#ifdef CLMDEBUG 
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
      lndname = trim(landdir)//'/LAI_patches'//trim(c2)//'.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'LAI_patches', 'patch', landpatch, LAI_patches, 1)
      
      lndname = trim(landdir)//'/LAI_pfts'//trim(c2)//'.nc'
      CALL ncio_create_file_vector (lndname, landpft)
      CALL ncio_define_dimension_vector (lndname, landpft, 'pft')
      CALL ncio_write_vector (lndname, 'LAI_pfts', 'pft', landpft, LAI_pfts, 1)
#else
      SITE_LAI_pfts_clim(:,month) = LAI_pfts(:)
#endif

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

#ifdef SinglePoint
   allocate (SITE_SAI_pfts_clim (numpft,12))  
#endif
   
   DO month = 1, 12
      IF (p_is_io) THEN
         CALL modis_read_data_pft_time (dir_modis, 'MONTHLY_SAI', gridlai, month, pftLSAI)
#ifdef USEMPI
         CALL aggregation_pft_data_daemon (gridlai, pftPCT, data3 = pftLSAI)
#endif
      ENDIF

      ! ---------------------------------------------------------------
      ! aggregate the plant leaf area index from the resolution of raw data to modelling resolution
      ! ---------------------------------------------------------------

      IF (p_is_worker) THEN
         DO ipatch = 1, numpatch

            CALL aggregation_pft_request_data (ipatch, gridlai, pftPCT, pct_pft_one, &
               area = area_one, data3 = pftLSAI, dout3 = sai_pft_one)
               
            IF (allocated(sai_one)) deallocate(sai_one)
            allocate(sai_one(size(area_one)))
            
            IF (allocated(pct_one)) deallocate(pct_one)
            allocate(pct_one(size(area_one)))

            pct_one = sum(pct_pft_one,dim=1)
            pct_one = max(pct_one, 1.0e-6)

            sai_one = sum(sai_pft_one * pct_pft_one, dim=1) / pct_one
            SAI_patches(ipatch) = sum(sai_one * area_one) / sum(area_one)
               
            IF (landpatch%settyp(ipatch) == 1) THEN
               DO ip = patch_pft_s(ipatch), patch_pft_e(ipatch)
                  p = landpft%settyp(ip)
                  sumarea = sum(pct_pft_one(p,:) * area_one)
                  IF (sumarea > 0) THEN
                     SAI_pfts(ip) = sum(sai_pft_one(p,:) * pct_pft_one(p,:) * area_one) / sumarea
                  ELSE
                     SAI_pfts(ip) = SAI_patches(ipatch)
                  ENDIF
               ENDDO
#ifdef CROP
            ELSEIF (landpatch%settyp(ipatch) == 12) THEN
               ip = patch_pft_s(ipatch)
               SAI_pfts(ip) = SAI_patches(ipatch)
#endif
            ENDIF
         ENDDO

#ifdef USEMPI
         CALL aggregation_pft_worker_done ()
#endif
      ENDIF

      write(c2,'(i2.2)') month
#ifdef CLMDEBUG 
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
      lndname = trim(landdir)//'/SAI_patches'//trim(c2)//'.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'SAI_patches', 'patch', landpatch, SAI_patches, 1)
      
      lndname = trim(landdir)//'/SAI_pfts'//trim(c2)//'.nc'
      CALL ncio_create_file_vector (lndname, landpft)
      CALL ncio_define_dimension_vector (lndname, landpft, 'pft')
      CALL ncio_write_vector (lndname, 'SAI_pfts', 'pft', landpft, SAI_pfts, 1)
#else
      SITE_SAI_pfts_clim(:,month) = SAI_pfts(:)
#endif

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
      CALL allocate_block_data (gridlai, pftLSAI, N_PFT_modis, lb1 = 0)
      CALL allocate_block_data (gridlai, pftPCT,  N_PFT_modis, lb1 = 0)
   ENDIF
     
   dir_modis = trim(DEF_dir_rawdata) // '/plant_15s_clim' 
      
   IF (p_is_io) THEN
      CALL modis_read_data_pft (dir_modis, 'PCT_PFT', gridlai, pftPCT)
   ENDIF

   IF (p_is_worker) THEN
      allocate(LAI_patches (numpatch))
      allocate(LAI_pcs (0:N_PFT-1, numpc))
   ENDIF

#ifdef SinglePoint
   allocate (SITE_LAI_pfts_clim (0:N_PFT-1,12))  
#endif
      
   DO month = 1, 12
      IF (p_is_io) THEN
         CALL modis_read_data_pft_time (dir_modis, 'MONTHLY_LAI', gridlai, month, pftLSAI)
#ifdef USEMPI
         CALL aggregation_pft_data_daemon (gridlai, pftPCT, data3 = pftLSAI)
#endif
      ENDIF

      ! ---------------------------------------------------------------
      ! aggregate the plant leaf area index from the resolution of raw data to modelling resolution
      ! ---------------------------------------------------------------

      IF (p_is_worker) THEN
         DO ipatch = 1, numpatch

            CALL aggregation_pft_request_data (ipatch, gridlai, pftPCT, pct_pft_one, &
               area = area_one, data3 = pftLSAI, dout3 = lai_pft_one)
               
            IF (allocated(lai_one)) deallocate(lai_one)
            allocate(lai_one(size(area_one)))

            IF (allocated(pct_one)) deallocate(pct_one)
            allocate(pct_one(size(area_one)))

            pct_one = sum(pct_pft_one,dim=1)
            pct_one = max(pct_one, 1.0e-6)

            lai_one = sum(lai_pft_one * pct_pft_one, dim=1) / pct_one
            LAI_patches(ipatch) = sum(lai_one * area_one) / sum(area_one)
               
            IF (patchtypes(landpatch%settyp(ipatch)) == 0) THEN
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

      write(c2,'(i2.2)') month
#ifdef CLMDEBUG 
      CALL check_vector_data ('LAI_patches ' // trim(c2), LAI_patches)
      CALL check_vector_data ('LAI_pcs     ' // trim(c2), LAI_pcs   )
#endif
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      ! ---------------------------------------------------
      ! write out the plant leaf area index of grid patches
      ! ---------------------------------------------------
#ifndef SinglePoint
      lndname = trim(landdir)//'/LAI_patches'//trim(c2)//'.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'LAI_patches', 'patch', landpatch, LAI_patches, 1)

      lndname = trim(landdir)//'/LAI_pcs'//trim(c2)//'.nc'
      CALL ncio_create_file_vector (lndname, landpc)
      CALL ncio_define_dimension_vector (lndname, landpc, 'pc')
      CALL ncio_define_dimension_vector (lndname, landpc, 'pft', N_PFT)
      CALL ncio_write_vector (lndname, 'LAI_pcs', 'pft', N_PFT, 'pc', landpc, LAI_pcs, 1)
#else
      SITE_LAI_pfts_clim(:,month) = LAI_pcs(:,1)
#endif

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
#ifdef SinglePoint
   allocate (SITE_SAI_pfts_clim (0:N_PFT-1,12))  
#endif
      
   DO month = 1, 12
      IF (p_is_io) THEN
         CALL modis_read_data_pft_time (dir_modis, 'MONTHLY_SAI', gridlai, month, pftLSAI)
#ifdef USEMPI
         CALL aggregation_pft_data_daemon (gridlai, pftPCT, data3 = pftLSAI)
#endif
      ENDIF

      ! ---------------------------------------------------------------
      ! aggregate the plant leaf area index from the resolution of raw data to modelling resolution
      ! ---------------------------------------------------------------

      IF (p_is_worker) THEN
         DO ipatch = 1, numpatch

            CALL aggregation_pft_request_data (ipatch, gridlai, pftPCT, pct_pft_one, &
               area = area_one, data3 = pftLSAI, dout3 = sai_pft_one)
               
            IF (allocated(sai_one)) deallocate(sai_one)
            allocate(sai_one(size(area_one)))

            IF (allocated(pct_one)) deallocate(pct_one)
            allocate(pct_one(size(area_one)))

            pct_one = sum(pct_pft_one,dim=1)
            pct_one = max(pct_one, 1.0e-6)

            sai_one = sum(sai_pft_one * pct_pft_one, dim=1) / pct_one
            SAI_patches(ipatch) = sum(sai_one * area_one) / sum(area_one)
               
            IF (patchtypes(landpatch%settyp(ipatch)) == 0) THEN
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

      write(c2,'(i2.2)') month
#ifdef CLMDEBUG 
      CALL check_vector_data ('SAI_patches ' // trim(c2), SAI_patches)
      CALL check_vector_data ('SAI_pcs     ' // trim(c2), SAI_pcs   )
#endif
#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      ! ---------------------------------------------------
      ! write out the plant stem area index of grid patches
      ! ---------------------------------------------------
#ifndef SinglePoint
      lndname = trim(landdir)//'/SAI_patches'//trim(c2)//'.nc'
      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'SAI_patches', 'patch', landpatch, SAI_patches, 1)
      
      lndname = trim(landdir)//'/SAI_pcs'//trim(c2)//'.nc'
      CALL ncio_create_file_vector (lndname, landpc)
      CALL ncio_define_dimension_vector (lndname, landpc, 'pc')
      CALL ncio_define_dimension_vector (lndname, landpc, 'pft', N_PFT)
      CALL ncio_write_vector (lndname, 'SAI_pcs', 'pft', N_PFT, 'pc', landpc, SAI_pcs, 1)
#else
      SITE_SAI_pfts_clim(:,month) = SAI_pcs(:,1)
#endif

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
