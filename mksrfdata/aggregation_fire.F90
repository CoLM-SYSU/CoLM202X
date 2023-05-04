#include <define.h>

SUBROUTINE aggregation_fire (gfire, dir_rawdata, dir_model_landdata)
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

   IMPLICIT NONE

   ! arguments:

   TYPE(grid_type),  intent(in) :: gfire
   CHARACTER(LEN=*), intent(in) :: dir_rawdata
   CHARACTER(LEN=*), intent(in) :: dir_model_landdata

   ! local variables:
   ! ----------------------------------------------------------------------
   CHARACTER(len=256) :: landdir, lndname

   TYPE (block_data_real8_2d) :: abm          ! global peak month of crop fire emissions (month)
   TYPE (block_data_real8_2d) :: hdm          ! Human population density
   TYPE (block_data_real8_2d) :: peatf        ! peatland fraction data
   TYPE (block_data_real8_2d) :: gdp          ! gdp data
   REAL(r8), allocatable :: abm_patches(:), abm_one(:), area_one(:)
   REAL(r8), allocatable :: hdm_patches(:), hdm_one(:)
   REAL(r8), allocatable :: peatf_patches(:), peatf_one(:)
   REAL(r8), allocatable :: gdp_patches(:), gdp_one(:)
   INTEGER :: itime, ipatch
   CHARACTER(LEN=4) :: c3, cyear
   integer :: start_year, end_year, YY

   landdir = trim(dir_model_landdata) // '/FIRE/'

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif
   IF (p_is_master) THEN
      write(*,'(/, A)') 'Aggregate FIRE ...'
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
      CALL allocate_block_data (gfire, abm)
      CALL allocate_block_data (gfire, hdm)
      CALL allocate_block_data (gfire, peatf)
      CALL allocate_block_data (gfire, gdp)
   ENDIF

   IF (p_is_worker) THEN
      allocate (abm_patches   (numpatch))
      allocate (hdm_patches   (numpatch))
      allocate (peatf_patches (numpatch))
      allocate (gdp_patches   (numpatch))
   ENDIF

   DO YY = start_year, end_year

      write(cyear,'(i4.4)') YY
      itime = max(1850,min(YY,2016)) - 1849
         !write(c1,'(i4.4)') YY
         ! ---------------------------
         ! read in hdm
         ! ---------------------------
      IF (p_is_master) THEN
         write(*,'(A,I4,A9,I4.4,A1,I3)') 'Aggregate Human population density (hdm):', YY,'(data in:',itime+1849,')'
      endif

      IF (p_is_io) THEN
              ! lndname = trim(dir_rawdata)//'/lai-true/'//trim(cyear)//'/NDEP_BNU_'//trim(cyear)//'_'//trim(c3)//'.h5'
         lndname = trim(dir_rawdata)//'/fire/colmforc.Li_2017_HYDEv3.2_CMIP6_hdm_0.5x0.5_AVHRR_simyr1850-2016_c180202.nc'
         CALL ncio_read_block_time (lndname, 'hdm', gfire, itime, hdm)
      ENDIF

#ifdef USEMPI
      CALL aggregation_data_daemon (gfire, data_r8_2d_in1 = hdm)
#endif

         ! ---------------------------------------------------------------
         ! aggregate the plant leaf area index from the resolution of raw data to modelling resolution
         ! ---------------------------------------------------------------

      IF (p_is_worker) THEN
         DO ipatch = 1, numpatch
            CALL aggregation_request_data (landpatch, ipatch, gfire, area = area_one, &
               data_r8_2d_in1 = hdm, data_r8_2d_out1 = hdm_one)
            hdm_patches(ipatch) = sum(hdm_one * area_one) / sum(area_one)
         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CLMDEBUG
      CALL check_vector_data ('hdm value ', hdm_patches)
#endif

         ! ---------------------------------------------------
         ! write out the plant leaf area index of grid patches
         ! ---------------------------------------------------
      lndname = trim(landdir) // '/hdm_'//trim(cyear)//'_patches.nc'

      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'hdm_patches', 'patch', landpatch, hdm_patches, 1)
   ENDDO

   IF (p_is_master) THEN
      write(*,'(A,I4,A1,I3,A1,I3)') 'Aggregate abm'
   endif

   itime = 1

   IF (p_is_io) THEN
          ! lndname = trim(dir_rawdata)//'/lai-true/'//trim(cyear)//'/NDEP_BNU_'//trim(cyear)//'_'//trim(c3)//'.h5'
      lndname = trim(dir_rawdata)//'/fire/abm_colm_double_fillcoast.nc'
      CALL ncio_read_block_time (lndname, 'abm', gfire, itime, abm)
   ENDIF

#ifdef USEMPI
      CALL aggregation_data_daemon (gfire, data_r8_2d_in1 = abm)
#endif

         ! ---------------------------------------------------------------
         ! aggregate the plant leaf area index from the resolution of raw data to modelling resolution
         ! ---------------------------------------------------------------

      IF (p_is_worker) THEN
         DO ipatch = 1, numpatch
            CALL aggregation_request_data (landpatch, ipatch, gfire, area = area_one, &
               data_r8_2d_in1 = abm, data_r8_2d_out1 = abm_one)
            abm_patches(ipatch) = sum(abm_one * area_one) / sum(area_one)
         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CLMDEBUG
      CALL check_vector_data ('abm value ', abm_patches)
#endif

         ! ---------------------------------------------------
         ! write out the plant leaf area index of grid patches
         ! ---------------------------------------------------
      lndname = trim(landdir) // '/abm_patches.nc'

      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'abm_patches', 'patch', landpatch, abm_patches, 1)


   IF (p_is_master) THEN
      write(*,'(A,I4,A1,I3,A1,I3)') 'Aggregate peatf'
   endif

   itime = 1

   IF (p_is_io) THEN
          ! lndname = trim(dir_rawdata)//'/lai-true/'//trim(cyear)//'/NDEP_BNU_'//trim(cyear)//'_'//trim(c3)//'.h5'
      lndname = trim(dir_rawdata)//'/fire/peatf_colm_360x720_c100428.nc'
      CALL ncio_read_block_time (lndname, 'peatf', gfire, itime, peatf)
   ENDIF

#ifdef USEMPI
      CALL aggregation_data_daemon (gfire, data_r8_2d_in1 = peatf)
#endif

         ! ---------------------------------------------------------------
         ! aggregate the plant leaf area index from the resolution of raw data to modelling resolution
         ! ---------------------------------------------------------------

      IF (p_is_worker) THEN
         DO ipatch = 1, numpatch
            CALL aggregation_request_data (landpatch, ipatch, gfire, area = area_one, &
               data_r8_2d_in1 = peatf, data_r8_2d_out1 = peatf_one)
            peatf_patches(ipatch) = sum(peatf_one * area_one) / sum(area_one)
         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CLMDEBUG
      CALL check_vector_data ('peatf value ', peatf_patches)
#endif

         ! ---------------------------------------------------
         ! write out the plant leaf area index of grid patches
         ! ---------------------------------------------------
      lndname = trim(landdir) // '/peatf_patches.nc'

      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'peatf_patches', 'patch', landpatch, peatf_patches, 1)

   IF (p_is_master) THEN
      write(*,'(A,I4,A1,I3,A1,I3)') 'Aggregate gdp'
   endif

   itime = 1

   IF (p_is_io) THEN
          ! lndname = trim(dir_rawdata)//'/lai-true/'//trim(cyear)//'/NDEP_BNU_'//trim(cyear)//'_'//trim(c3)//'.h5'
      lndname = trim(dir_rawdata)//'/fire/gdp_colm_360x720_c100428.nc'
      CALL ncio_read_block_time (lndname, 'gdp', gfire, itime, gdp)
   ENDIF

#ifdef USEMPI
      CALL aggregation_data_daemon (gfire, data_r8_2d_in1 = gdp)
#endif

         ! ---------------------------------------------------------------
         ! aggregate the plant leaf area index from the resolution of raw data to modelling resolution
         ! ---------------------------------------------------------------

      IF (p_is_worker) THEN
         DO ipatch = 1, numpatch
            CALL aggregation_request_data (landpatch, ipatch, gfire, area = area_one, &
               data_r8_2d_in1 = gdp, data_r8_2d_out1 = gdp_one)
            gdp_patches(ipatch) = sum(gdp_one * area_one) / sum(area_one)
         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CLMDEBUG
      CALL check_vector_data ('gdp value ', gdp_patches)
#endif

         ! ---------------------------------------------------
         ! write out the plant leaf area index of grid patches
         ! ---------------------------------------------------
      lndname = trim(landdir) // '/gdp_patches.nc'

      CALL ncio_create_file_vector (lndname, landpatch)
      CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
      CALL ncio_write_vector (lndname, 'gdp_patches', 'patch', landpatch, gdp_patches, 1)

   IF (p_is_worker) THEN
      IF (allocated(hdm_patches)) deallocate(hdm_patches)
      IF (allocated(hdm_one    )) deallocate(hdm_one    )
      IF (allocated(abm_patches)) deallocate(abm_patches)
      IF (allocated(abm_one    )) deallocate(abm_one    )
      IF (allocated(peatf_patches)) deallocate(peatf_patches)
      IF (allocated(peatf_one    )) deallocate(peatf_one    )
      IF (allocated(gdp_patches)) deallocate(gdp_patches)
      IF (allocated(gdp_one    )) deallocate(gdp_one    )
      IF (allocated(area_one   )) deallocate(area_one    )
   ENDIF

END SUBROUTINE aggregation_fire
