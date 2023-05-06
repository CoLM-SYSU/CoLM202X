#include <define.h>
#ifdef CROP
SUBROUTINE aggregation_crop_parameters (gridcrop, dir_rawdata, dir_model_landdata)
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

   USE mod_aggregation

   USE LC_Const
#ifdef PFT_CLASSIFICATION
   USE mod_landpft
#endif
#ifdef PC_CLASSIFICATION
   USE mod_landpc
#endif

   IMPLICIT NONE

   ! arguments:

   TYPE(grid_type),  intent(in) :: gridcrop
   CHARACTER(LEN=*), intent(in) :: dir_rawdata
   CHARACTER(LEN=*), intent(in) :: dir_model_landdata

   ! local variables:
   ! ----------------------------------------------------------------------
   CHARACTER(len=256) :: landdir, lndname, lndname_out

   TYPE (block_data_real8_2d) :: plantdate_rice2
   TYPE (block_data_real8_2d) :: plantdate
   TYPE (block_data_real8_2d) :: fertnitro          ! plant leaf area index (m2/m2)
   TYPE (block_data_real8_2d) :: pct_cft            ! fraction of crop patch (unitless)
   REAL(r8), allocatable :: plantdate_one(:)
   REAL(r8), allocatable :: fertnitro_one(:)
   REAL(r8), allocatable :: pct_cft_one(:)
   REAL(r8), allocatable :: area_one(:)
   REAL(r8), allocatable :: plantdate_rice2_patches(:)
   REAL(r8), allocatable :: plantdate_pfts(:)
   REAL(r8), allocatable :: fertnitro_pfts(:)

   INTEGER :: itime, ntime, Julian_day, ipatch, ipft
   CHARACTER(LEN=4) ::cx, c2, c3, c4, cyear,c
   integer :: start_year, end_year, YY,nsl, cft

   landdir = trim(dir_model_landdata) // '/crop/'

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif
   IF (p_is_master) THEN
      write(*,'(/, A)') 'Aggregate planting date ...'
      CALL system('mkdir -p ' // trim(adjustl(landdir)))
   ENDIF
#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif



   ! ................................................
   ! ... global ****
   ! ................................................

!   IF (DEF_LAI_CLIM) THEN
      start_year = 1
      end_year   = 1
      ntime = 1
!   ELSE
!      start_year = DEF_simulation_time%start_year
!      end_year   = DEF_simulation_time%end_year
!      ntime = 46
!   ENDIF

   IF (p_is_io) THEN
      CALL allocate_block_data (gridcrop, plantdate_rice2)
   ENDIF

   IF (p_is_worker) THEN
      allocate (plantdate_rice2_patches     (numpatch))
   ENDIF


!   DO YY = start_year, end_year
   DO itime = 1, ntime
         ! -----------------------
         ! read in leaf area index
         ! -----------------------
      write(c3, '(i2.2)') itime
        ! IF (p_is_master) THEN
        !    write(*,'(A,I4,A1,I3,A1,I3)') 'Aggregate CONC_O2_UNSAT Level:',cx, ':', YY, ':', itime, '/', ntime
        ! endif

      IF (p_is_io) THEN
         lndname = trim(dir_rawdata)//'/crop/plantdt-colm-64cfts-rice2_fillcoast.nc'
         print *, lndname
         CALL ncio_read_block_time (lndname, 'pdrice2', gridcrop, itime, plantdate_rice2)

      ENDIF

!------------distribute plantdate_rice2 to worker---------------------
#ifdef USEMPI
      IF (p_is_io) THEN
         CALL aggregation_data_daemon (gridcrop, data_r8_2d_in1 = plantdate_rice2)
      ENDIF
#endif

      IF (p_is_worker) THEN
         DO ipatch = 1, numpatch
            CALL aggregation_request_data (landpatch, ipatch, gridcrop, area = area_one, &
               data_r8_2d_in1 = plantdate_rice2, data_r8_2d_out1 = plantdate_one)
            plantdate_rice2_patches(ipatch) = sum(plantdate_one * area_one) / sum(area_one)
         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CLMDEBUG
      CALL check_vector_data ('plant date value for rice2 '//trim(c3), plantdate_rice2_patches)
#endif

         ! ---------------------------------------------------
         ! write out the plant leaf area index of grid patches
         ! ---------------------------------------------------
      lndname = trim(landdir) // '/plantdate_patches.nc'

      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'plantdate_rice2_patches', 'patch', landpatch, plantdate_rice2_patches, 1)
   ENDDO


#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif
   IF (p_is_master) THEN
      write(*,'(/, A)') 'Aggregate plantdate ...'
      CALL system('mkdir -p ' // trim(adjustl(landdir)))
   ENDIF
#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

   ! ................................................
   ! ... global ****
   ! ................................................

!   IF (DEF_LAI_CLIM) THEN
      start_year = 1
      end_year   = 1
      ntime = 1
!   ELSE
!      start_year = DEF_simulation_time%start_year
!      end_year   = DEF_simulation_time%end_year
!      ntime = 46
!   ENDIF

   ! ----- CONC_O2_UNSAT -----
   IF (p_is_io) THEN
      CALL allocate_block_data (gridcrop, pct_cft)
      CALL allocate_block_data (gridcrop, plantdate)
   ENDIF

   IF (p_is_worker) THEN
      allocate (plantdate_pfts (numpft))
   ENDIF

   IF (p_is_io) THEN
      lndname = trim(dir_rawdata)//'/crop/plantdt-colm-64cfts-rice2_fillcoast.nc'
      print *, lndname
      call system('ls -l '//trim(lndname))
   ENDIF

   lndname_out = trim(landdir) // '/plantdate_pfts.nc'

   CALL ncio_create_file_vector (lndname_out, landpft)
   CALL ncio_define_dimension_vector (lndname_out, landpft, 'pft')

   DO cft = 15,78
      write(c4, '(i2.2)') cft
      IF (p_is_io) THEN
         CALL ncio_read_block_time (lndname, 'PLANTDATE_CFT_'//trim(c4), gridcrop, 1, plantdate)
      ENDIF

#ifdef USEMPI
      IF (p_is_io) THEN
         CALL aggregation_data_daemon (gridcrop, data_r8_2d_in1 = plantdate)
      ENDIF
#endif

      IF (p_is_worker) THEN
         DO ipatch = 1, numpatch
            CALL aggregation_request_data (landpatch, ipatch, gridcrop, area = area_one, &
               data_r8_2d_in1 = plantdate, data_r8_2d_out1 = plantdate_one)
            IF(landpatch%settyp(ipatch) .eq. 12)THEN
               DO ipft = patch_pft_s(ipatch),patch_pft_e(ipatch)
                  IF(landpft%settyp(ipft) .eq. cft)THEN
                     IF(sum(area_one) .ne. 0)THEN
                        plantdate_pfts(ipft) = sum(plantdate_one * area_one) / sum(area_one)
                     ELSE
                        write(*,*),'cft:'//trim(c4)//' crop distribution mismatch between model surface data and fertilization data,ipatch:'&
                                                   ,ipatch,'ipft',ipft,'p_iam_glb',p_iam_glb
                        write(*,*),'pct_cft_one',pct_cft_one
                        write(*,*),'plantdate_one',plantdate_one
                        call abort
                     END IF
                  ENDIF
               ENDDO
            ELSE
               IF (landpatch%settyp(ipatch) == 1) THEN
                  DO ipft = patch_pft_s(ipatch),patch_pft_e(ipatch)
                     plantdate_pfts(ipft) = -9999._r8
                  ENDDO
               ENDIF
            ENDIF
         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

!               CALL aggregation_lc_request_data (ipatch, gridcrop, plantdate_rice2, plantdate_one, area_one)
!               CALL aggregation_lc_request_data (ipatch, gridcrop, cftfrac_rice2, cftfrac_one)
!               DO m = patch_pft_s(ipatch), patch_pft_e(ipatch)
!                  IF(pftclass(m) .eq. 17 .or. pftclass(m) .eq. 18)then
!                     IF(sum(area_one * cftfrac_one * PCT_CROP_one) .ne. 0)then
!                        plantdate_pft(m) = sum(plantdate_one * area_one * cftfrac_one * PCT_CROP_one) / sum(area_one * cftfrac_one * PCT_CROP_one)
!                     ELSE
!                        write(*,*),'cft distribution data and plant date (rice2) data mismatch, ipatch=',ipatch,p_iam_glb
!                        call abort
!                     END IF
!                  END IF
!               END DO


         ! ---------------------------------------------------------------
         ! aggregate the plant leaf area index from the resolution of raw data to modelling resolution
         ! ---------------------------------------------------------------

         ! ---------------------------------------------------
         ! write out the plant leaf area index of grid patches
         ! ---------------------------------------------------
   ENDDO

#ifdef CLMDEBUG
   if(p_is_worker)then
      CALL check_vector_data ('plantdate_pfts value '//trim(c4), plantdate_pfts)
   endif
#endif

   CALL ncio_write_vector (lndname_out, 'plantdate_pfts', 'pft', landpft, plantdate_pfts, 1)


#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif
   IF (p_is_master) THEN
      write(*,'(/, A)') 'Aggregate fertnitro ...'
      CALL system('mkdir -p ' // trim(adjustl(landdir)))
   ENDIF
#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

   ! ................................................
   ! ... global ****
   ! ................................................

!   IF (DEF_LAI_CLIM) THEN
      start_year = 1
      end_year   = 1
      ntime = 1
!   ELSE
!      start_year = DEF_simulation_time%start_year
!      end_year   = DEF_simulation_time%end_year
!      ntime = 46
!   ENDIF

   ! ----- CONC_O2_UNSAT -----
   IF (p_is_io) THEN
      CALL allocate_block_data (gridcrop, pct_cft)
      CALL allocate_block_data (gridcrop, fertnitro)
   ENDIF

   IF (p_is_worker) THEN
      allocate (fertnitro_pfts (numpft))
   ENDIF

   IF (p_is_io) THEN
      lndname = trim(dir_rawdata)//'/crop/fertnitro_fillcoast.nc'
      print *, lndname
      call system('ls -l '//trim(lndname))
   ENDIF

   lndname_out = trim(landdir) // '/fertnitro_pfts.nc'

   CALL ncio_create_file_vector (lndname_out, landpft)
   CALL ncio_define_dimension_vector (lndname_out, landpft, 'pft')

   DO cft = 15,78
      write(c4, '(i2.2)') cft
      IF (p_is_io) THEN
         CALL ncio_read_block_time (lndname, 'CONST_FERTNITRO_CFT_'//trim(c4), gridcrop, 1, fertnitro)
         CALL ncio_read_block_time (lndname, 'PCT_CFT_'//trim(c4), gridcrop, 1, pct_cft)
      ENDIF

#ifdef USEMPI
      IF (p_is_io) THEN
         CALL aggregation_data_daemon (gridcrop, data_r8_2d_in1 = fertnitro, data_r8_2d_in2 = pct_cft)
      ENDIF
#endif

      IF (p_is_worker) THEN
         DO ipatch = 1, numpatch
            CALL aggregation_request_data (landpatch, ipatch, gridcrop, area = area_one, &
               data_r8_2d_in1 = fertnitro, data_r8_2d_out1 = fertnitro_one, &
               data_r8_2d_in2 = pct_cft  , data_r8_2d_out2 = pct_cft_one    )
            IF(landpatch%settyp(ipatch) .eq. 12)THEN
               DO ipft = patch_pft_s(ipatch),patch_pft_e(ipatch)
                  IF(landpft%settyp(ipft) .eq. cft)THEN
                     IF(sum(area_one) .ne. 0)THEN
                        fertnitro_pfts(ipft) = sum(fertnitro_one * area_one) / sum(area_one)
                     ELSE
                        write(*,*),'cft:'//trim(c4)//' crop distribution mismatch between model surface data and fertilization data,ipatch:'&
                                                   ,ipatch,'ipft',ipft,'p_iam_glb',p_iam_glb
                        write(*,*),'pct_cft_one',pct_cft_one
                        write(*,*),'fertnitro_one',fertnitro_one
                        call abort
                     END IF
                  ENDIF
               ENDDO
            ELSE
               IF (landpatch%settyp(ipatch) == 1) THEN
                  DO ipft = patch_pft_s(ipatch),patch_pft_e(ipatch)
                     fertnitro_pfts(ipft) = -9999._r8
                  ENDDO
               ENDIF
            ENDIF
         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

!               CALL aggregation_lc_request_data (ipatch, gridcrop, plantdate_rice2, plantdate_one, area_one)
!               CALL aggregation_lc_request_data (ipatch, gridcrop, cftfrac_rice2, cftfrac_one)
!               DO m = patch_pft_s(ipatch), patch_pft_e(ipatch)
!                  IF(pftclass(m) .eq. 17 .or. pftclass(m) .eq. 18)then
!                     IF(sum(area_one * cftfrac_one * PCT_CROP_one) .ne. 0)then
!                        plantdate_pft(m) = sum(plantdate_one * area_one * cftfrac_one * PCT_CROP_one) / sum(area_one * cftfrac_one * PCT_CROP_one)
!                     ELSE
!                        write(*,*),'cft distribution data and plant date (rice2) data mismatch, ipatch=',ipatch,p_iam_glb
!                        call abort
!                     END IF
!                  END IF
!               END DO


         ! ---------------------------------------------------------------
         ! aggregate the plant leaf area index from the resolution of raw data to modelling resolution
         ! ---------------------------------------------------------------

         ! ---------------------------------------------------
         ! write out the plant leaf area index of grid patches
         ! ---------------------------------------------------
   ENDDO

#ifdef CLMDEBUG
   if(p_is_worker)then
      CALL check_vector_data ('fert nitro value '//trim(c4), fertnitro_pfts)
   endif
#endif

   CALL ncio_write_vector (lndname_out, 'fertnitro_pfts', 'pft', landpft, fertnitro_pfts, 1)




END SUBROUTINE aggregation_crop_parameters
#endif
