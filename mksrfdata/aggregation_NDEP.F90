#include <define.h>

SUBROUTINE aggregation_NDEP (gridndep, dir_rawdata, dir_model_landdata)
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

   IMPLICIT NONE

   ! arguments:

   TYPE(grid_type),  intent(in) :: gridndep
   CHARACTER(LEN=*), intent(in) :: dir_rawdata
   CHARACTER(LEN=*), intent(in) :: dir_model_landdata

   ! local variables:
   ! ----------------------------------------------------------------------
   CHARACTER(len=256) :: landdir, lndname

   TYPE (block_data_real8_2d) :: NDEP          ! nitrogen deposition (gN/m2/s)
   REAL(r8), allocatable :: NDEP_patches(:), ndep_one(:), area_one(:)
   INTEGER :: itime, ipatch
   CHARACTER(LEN=4) :: cyear
   integer :: start_year, end_year, YY


   landdir = trim(dir_model_landdata) // '/NDEP/'

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif
   IF (p_is_master) THEN
      write(*,'(/, A)') 'Aggregate NDEP ...'
      CALL system('mkdir -p ' // trim(adjustl(landdir)))
   ENDIF
#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif

!ifdef SinglePoint
!   IF (USE_SITE_NDEP) THEN
!     RETURN
!  ENDIF
!#endif

   ! ................................................
   ! ... global plant leaf area index
   ! ................................................

   start_year = DEF_simulation_time%start_year
   end_year   = DEF_simulation_time%end_year

   ! ----- NDEP -----
   IF (p_is_io) THEN
      CALL allocate_block_data (gridndep, NDEP)
   ENDIF

   IF (p_is_worker) THEN
      allocate (NDEP_patches (numpatch))
   ENDIF

   DO YY = start_year, end_year

      write(cyear,'(i4.4)') YY
      itime = max(min(YY,2006),1849) - 1848
         !write(c1,'(i4.4)') YY
         ! ---------------------------
         ! read in nitrofen deposition
         ! ---------------------------
      IF (p_is_master) THEN
         write(*,'(A,I4,A9,I4.4,A1,I3)') 'Aggregate NDEP:', YY,'(data in:',itime+1848,')'
      endif

      IF (p_is_io) THEN
         lndname = trim(dir_rawdata)//'/ndep/fndep_colm_hist_simyr1849-2006_1.9x2.5_c100428.nc'
         CALL ncio_read_block_time (lndname, 'NDEP_year', gridndep, itime, NDEP)
      ENDIF

#ifdef USEMPI
         CALL aggregation_data_daemon (gridndep, data_r8_2d_in1 = NDEP)
#endif

         ! ---------------------------------------------------------------
         ! aggregate the plant leaf area index from the resolution of raw data to modelling resolution
         ! ---------------------------------------------------------------

      IF (p_is_worker) THEN
         DO ipatch = 1, numpatch
            CALL aggregation_request_data (landpatch, ipatch, gridndep, area = area_one, &
               data_r8_2d_in1 = NDEP, data_r8_2d_out1 = ndep_one)
            NDEP_patches(ipatch) = sum(ndep_one * area_one) / sum(area_one)
         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CLMDEBUG
      CALL check_vector_data ('NDEP value ', NDEP_patches)
#endif

         ! ---------------------------------------------------
         ! write out the plant leaf area index of grid patches
         ! ---------------------------------------------------
      lndname = trim(landdir) // '/NDEP_'//trim(cyear)//'_patches.nc'

      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'NDEP_patches', 'patch', landpatch, NDEP_patches, 1)
   ENDDO


   IF (p_is_worker) THEN
      IF (allocated(NDEP_patches)) deallocate(NDEP_patches)
      IF (allocated(ndep_one    )) deallocate(ndep_one    )
      IF (allocated(area_one    )) deallocate(area_one    )
   ENDIF

END SUBROUTINE aggregation_NDEP
