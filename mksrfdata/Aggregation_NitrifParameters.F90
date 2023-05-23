#include <define.h>

SUBROUTINE Aggregation_NitrifParameters (gridnitrif, dir_rawdata, dir_model_landdata)
   ! ----------------------------------------------------------------------
   ! !DESCRIPTION:
   ! Global nitrification data (climatological data of soil O2) from CLM5 simulation
   !
   ! Created by Zhongwang Wei and modified by Xingjie Lu, 09/2022
   !
   ! ----------------------------------------------------------------------
   USE precision
   USE GlobalVars
   USE mod_namelist
   USE spmd_task
   USE mod_grid
   USE MOD_LandPatch
   USE ncio_block
   USE ncio_vector
#ifdef CoLMDEBUG 
   USE mod_colm_debug
#endif

   USE MOD_AggregationRequestData

   USE LC_Const
#ifdef PFT_CLASSIFICATION
   USE MOD_LandPFT
#endif
#ifdef PC_CLASSIFICATION
   USE MOD_LandPC
#endif

#ifdef SrfdataDiag
   USE MOD_SrfdataDiag
#endif

   IMPLICIT NONE

   ! arguments:

   TYPE(grid_type),  intent(in) :: gridnitrif
   CHARACTER(LEN=*), intent(in) :: dir_rawdata
   CHARACTER(LEN=*), intent(in) :: dir_model_landdata

   ! local variables:
   ! ----------------------------------------------------------------------
   CHARACTER(len=256) :: landdir, lndname

   TYPE (block_data_real8_2d) :: CONC_O2_UNSAT, O2_DECOMP_DEPTH_UNSAT          ! plant leaf area index (m2/m2)
   REAL(r8), allocatable :: CONC_O2_UNSAT_patches(:), CONC_O2_UNSAT_one(:), area_one(:)
   REAL(r8), allocatable :: O2_DECOMP_DEPTH_UNSAT_patches(:), O2_DECOMP_DEPTH_UNSAT_one(:)

   INTEGER :: itime, ntime, Julian_day, ipatch
   CHARACTER(LEN=4) ::cx, c2, c3, cyear,c
   integer :: start_year, end_year, YY,nsl

#ifdef SrfdataDiag
   INTEGER :: typpatch(N_land_classification+1), ityp
   CHARACTER(len=256) :: varname
#endif

   landdir = trim(dir_model_landdata) // '/nitrif/'

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif
   IF (p_is_master) THEN
      write(*,'(/, A)') 'Aggregate CONC_O2_UNSAT ...'
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
      ntime = 12 
!   ELSE
!      start_year = DEF_simulation_time%start_year
!      end_year   = DEF_simulation_time%end_year
!      ntime = 46
!   ENDIF

   ! ----- CONC_O2_UNSAT -----
   IF (p_is_io) THEN
      CALL allocate_block_data (gridnitrif, CONC_O2_UNSAT)
   ENDIF
   
   IF (p_is_worker) THEN
      allocate (CONC_O2_UNSAT_patches (numpatch))
   ENDIF
   
DO nsl = 1, 20
   write(cx,'(i2.2)') nsl
   
   DO YY = start_year, end_year
      DO itime = 1, ntime
         ! -----------------------
         ! read in leaf area index
         ! -----------------------
         write(c3, '(i2.2)') itime
        ! IF (p_is_master) THEN
        !    write(*,'(A,I4,A1,I3,A1,I3)') 'Aggregate CONC_O2_UNSAT Level:',cx, ':', YY, ':', itime, '/', ntime 
        ! endif

         IF (p_is_io) THEN   
               lndname = trim(dir_rawdata)//'/nitrif/CONC_O2_UNSAT/CONC_O2_UNSAT_l'//trim(cx)//'.nc'
               print *, lndname
               CALL ncio_read_block_time (lndname, 'CONC_O2_UNSAT', gridnitrif, itime, CONC_O2_UNSAT)

#ifdef USEMPI
            CALL aggregation_data_daemon (gridnitrif, data_r8_2d_in1 = CONC_O2_UNSAT)
#endif
         ENDIF


         ! ---------------------------------------------------------------
         ! aggregate the plant leaf area index from the resolution of raw data to modelling resolution
         ! ---------------------------------------------------------------

         IF (p_is_worker) THEN
            DO ipatch = 1, numpatch
               CALL aggregation_request_data (landpatch, ipatch, gridnitrIF, area = area_one, &
                  data_r8_2d_in1 = CONC_O2_UNSAT, data_r8_2d_out1 = CONC_O2_UNSAT_one)
               CONC_O2_UNSAT_patches(ipatch) = sum(CONC_O2_UNSAT_one * area_one) / sum(area_one)
            ENDDO

#ifdef USEMPI
            CALL aggregation_worker_done ()
#endif
         ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CoLMDEBUG 
         CALL check_vector_data ('CONC_O2_UNSAT value '//trim(c3), CONC_O2_UNSAT_patches)
#endif

         ! ---------------------------------------------------
         ! write out the plant leaf area index of grid patches
         ! ---------------------------------------------------
         lndname = trim(landdir) // '/CONC_O2_UNSAT_patches_l' // trim(cx)//'_'// trim(c3) // '.nc'
       

         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'CONC_O2_UNSAT_patches', 'patch', landpatch, CONC_O2_UNSAT_patches, 1)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/CONC_O2_UNSAT_patch.nc'
         varname = 'CONC_O2_UNSAT_l' // trim(cx) // '_' // trim(c3)
         CALL srfdata_map_and_write (CONC_O2_UNSAT_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, trim(varname), compress = 1, write_mode = 'one')
#endif
      ENDDO
   ENDDO
   

enddo

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------

#ifdef USEMPI
   CALL mpi_barrier (p_comm_glb, p_err)
#endif
   IF (p_is_master) THEN
      write(*,'(/, A)') 'Aggregate O2_DECOMP_DEPTH_UNSAT ...'
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
      ntime = 12 
!   ELSE
!      start_year = DEF_simulation_time%start_year
!      end_year   = DEF_simulation_time%end_year
!      ntime = 46
!   ENDIF

   ! ----- CONC_O2_UNSAT -----
   IF (p_is_io) THEN
      CALL allocate_block_data (gridnitrif, O2_DECOMP_DEPTH_UNSAT)
   ENDIF
   
   IF (p_is_worker) THEN
      allocate (O2_DECOMP_DEPTH_UNSAT_patches (numpatch))
   ENDIF
   
DO nsl = 1, 25
   write(cx,'(i2.2)') nsl
   
   DO YY = start_year, end_year
      DO itime = 1, ntime
         ! -----------------------
         ! read in leaf area index
         ! -----------------------
         write(c3, '(i2.2)') itime
        ! IF (p_is_master) THEN
        !    write(*,'(A,I4,A1,I3,A1,I3)') 'Aggregate CONC_O2_UNSAT Level:',cx, ':', YY, ':', itime, '/', ntime 
        ! endif

         IF (p_is_io) THEN   
               lndname = trim(dir_rawdata)//'/nitrif/O2_DECOMP_DEPTH_UNSAT/O2_DECOMP_DEPTH_UNSAT_l'//trim(cx)//'.nc'
               print *, lndname
               CALL ncio_read_block_time (lndname, 'O2_DECOMP_DEPTH_UNSAT', gridnitrif, itime, O2_DECOMP_DEPTH_UNSAT)

#ifdef USEMPI
            CALL aggregation_data_daemon (gridnitrif, data_r8_2d_in1 = O2_DECOMP_DEPTH_UNSAT)
#endif
         ENDIF


         ! ---------------------------------------------------------------
         ! aggregate the plant leaf area index from the resolution of raw data to modelling resolution
         ! ---------------------------------------------------------------

         IF (p_is_worker) THEN
            DO ipatch = 1, numpatch
               CALL aggregation_request_data (landpatch, ipatch, gridnitrIF, area = area_one, &
                  data_r8_2d_in1 = O2_DECOMP_DEPTH_UNSAT, data_r8_2d_out1 = O2_DECOMP_DEPTH_UNSAT_one)
               O2_DECOMP_DEPTH_UNSAT_patches(ipatch) = sum(O2_DECOMP_DEPTH_UNSAT_one * area_one) / sum(area_one)
            ENDDO

#ifdef USEMPI
            CALL aggregation_worker_done ()
#endif
         ENDIF

#ifdef USEMPI
         CALL mpi_barrier (p_comm_glb, p_err)
#endif

#ifdef CoLMDEBUG 
         CALL check_vector_data ('O2_DECOMP_DEPTH_UNSAT value '//trim(c3), O2_DECOMP_DEPTH_UNSAT_patches)
#endif

         ! ---------------------------------------------------
         ! write out the plant leaf area index of grid patches
         ! ---------------------------------------------------
         lndname = trim(landdir) // '/O2_DECOMP_DEPTH_UNSAT_patches_l' // trim(cx)//'_'// trim(c3) // '.nc'
       

         CALL ncio_create_file_vector (lndname, landpatch)
         CALL ncio_define_dimension_vector (lndname, landpatch, 'patch')
         CALL ncio_write_vector (lndname, 'O2_DECOMP_DEPTH_UNSAT_patches', 'patch', landpatch, O2_DECOMP_DEPTH_UNSAT_patches, 1)

#ifdef SrfdataDiag
         typpatch = (/(ityp, ityp = 0, N_land_classification)/)
         lndname  = trim(dir_model_landdata) // '/diag/O2_DECOMP_DEPTH_UNSAT_patch.nc'
         varname = 'O2_DECOMP_DEPTH_UNSAT_l' // trim(cx) // '_' // trim(c3)
         CALL srfdata_map_and_write (O2_DECOMP_DEPTH_UNSAT_patches, landpatch%settyp, typpatch, m_patch2diag, &
            -1.0e36_r8, lndname, trim(varname), compress = 1, write_mode = 'one')
#endif
      ENDDO
   ENDDO
   

enddo


END SUBROUTINE Aggregation_NitrifParameters
