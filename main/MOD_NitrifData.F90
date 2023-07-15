#include <define.h>

#ifdef BGC
MODULE MOD_NitrifData
 !-----------------------------------------------------------------------
 ! !DESCRIPTION:
 ! This module read in nitrif data.
 !
 ! !ORIGINAL:
 ! Lu Xingjie and Zhang Shupeng, 2023, prepare the original version of the nitrif data module.

   USE MOD_Grid
   USE MOD_Mapping_Grid2Pset
   use MOD_BGC_Vars_TimeVariables, only : tCONC_O2_UNSAT, tO2_DECOMP_DEPTH_UNSAT
   IMPLICIT NONE

   TYPE(grid_type) :: grid_nitrif
   type(mapping_grid2pset_type) :: mg2p_nitrif

CONTAINS

   ! ----------
   SUBROUTINE init_nitrif_data (idate)

   !----------------------
   ! DESCTIPTION:
   ! open nitrif netcdf file from DEF_dir_runtime, read latitude and longitude info.
   ! Initialize nitrif data read in.

      use MOD_TimeManager
      USE MOD_Namelist
      USE MOD_Grid
      USE MOD_NetCDFSerial
      USE MOD_LandPatch
      IMPLICIT NONE

      integer, intent(in) :: idate(3)

      ! Local Variables
      CHARACTER(len=256) :: file_nitrif
      REAL(r8), allocatable :: lat(:), lon(:)
      integer :: month, mday

      file_nitrif = trim(DEF_dir_runtime)//'/nitrif/CONC_O2_UNSAT/CONC_O2_UNSAT_l01.nc'

      CALL ncio_read_bcast_serial (file_nitrif, 'lat', lat)
      CALL ncio_read_bcast_serial (file_nitrif, 'lon', lon)

      CALL grid_nitrif%define_by_center (lat, lon)

      call mg2p_nitrif%build (grid_nitrif, landpatch)

      IF (allocated(lon)) deallocate(lon)
      IF (allocated(lat)) deallocate(lat)

      CALL julian2monthday (idate(1), idate(2), month, mday)

      CALL update_nitrif_data (month)

   END SUBROUTINE init_nitrif_data

   ! ----------
   SUBROUTINE update_nitrif_data (month)

      use MOD_SPMD_Task
      use MOD_Namelist
      USE MOD_DataType
      USE MOD_Vars_Global, only : nl_soil
      USE MOD_NetCDFBlock
      use MOD_LandPatch
      use MOD_Vars_TimeInvariants
      USE MOD_RangeCheck
      IMPLICIT NONE

      integer,  intent(in) :: month

      ! Local Variables
      CHARACTER(len=256) :: file_nitrif
      TYPE(block_data_real8_2d) :: f_xy_nitrif
      REAL(r8), allocatable :: tCONC_O2_UNSAT_tmp(:)
      REAL(r8), allocatable :: tO2_DECOMP_DEPTH_UNSAT_tmp(:)
      character(len=2) :: cx
      integer :: nsl, npatch, m

      IF (p_is_worker) THEN
         allocate(tCONC_O2_UNSAT_tmp        (numpatch))
         allocate(tO2_DECOMP_DEPTH_UNSAT_tmp(numpatch))
      ENDIF 
      
      IF (p_is_io) THEN
         CALL allocate_block_data (grid_nitrif, f_xy_nitrif)
      ENDIF

      DO nsl = 1, nl_soil
      
         write(cx,'(i2.2)') nsl
         file_nitrif = trim(DEF_dir_runtime)//'/nitrif/CONC_O2_UNSAT/CONC_O2_UNSAT_l'//trim(cx)//'.nc'
         IF (p_is_io) THEN
            CALL ncio_read_block_time (file_nitrif, 'CONC_O2_UNSAT', grid_nitrif, month, f_xy_nitrif)
         ENDIF

         call mg2p_nitrif%map_aweighted (f_xy_nitrif, tCONC_O2_UNSAT_tmp)

         if (p_is_worker) then
            if (numpatch > 0) then
               do npatch = 1, numpatch
                  m = patchclass(npatch)
                  if( m == 0 )then
                     tCONC_O2_UNSAT(nsl,npatch)  = 0.
                  else
                     tCONC_O2_UNSAT(nsl,npatch)  = tCONC_O2_UNSAT_tmp(npatch)
                  endif
                  if (tCONC_O2_UNSAT(nsl,npatch) < 1E-10) then
                     tCONC_O2_UNSAT(nsl,npatch)=0.0
                  endif
               end do

            ENDIF
         ENDIF
      END do

#ifdef RangeCheck
      call check_vector_data ('CONC_O2_UNSAT', tCONC_O2_UNSAT)
#endif

      DO nsl = 1, nl_soil
         
         write(cx,'(i2.2)') nsl
         file_nitrif = trim(DEF_dir_runtime)//'/nitrif/O2_DECOMP_DEPTH_UNSAT_l'//trim(cx)//'.nc'
         IF (p_is_io) THEN
            CALL ncio_read_block_time (file_nitrif, 'O2_DECOMP_DEPTH_UNSAT', grid_nitrif, month, f_xy_nitrif)
         ENDIF 

         call mg2p_nitrif%map_aweighted (f_xy_nitrif, tO2_DECOMP_DEPTH_UNSAT_tmp)

         if (p_is_worker) then
            if (numpatch > 0) then
               do npatch = 1, numpatch
                  m = patchclass(npatch)
                  if( m == 0 )then
                     tO2_DECOMP_DEPTH_UNSAT(nsl,npatch)  = 0.
                  else
                     tO2_DECOMP_DEPTH_UNSAT(nsl,npatch)  = tO2_DECOMP_DEPTH_UNSAT_tmp(npatch)
                  endif
                  if (tO2_DECOMP_DEPTH_UNSAT(nsl,npatch) < 1E-10) then
                     tO2_DECOMP_DEPTH_UNSAT(nsl,npatch)=0.0
                  endif
               end do

            ENDIF
         ENDIF
      END do

#ifdef RangeCheck
      call check_vector_data ('O2_DECOMP_DEPTH_UNSAT', tO2_DECOMP_DEPTH_UNSAT)
#endif

      IF (p_is_worker) THEN
         deallocate (tCONC_O2_UNSAT_tmp)
         deallocate (tO2_DECOMP_DEPTH_UNSAT_tmp)
      ENDIF

   END SUBROUTINE update_nitrif_data

END MODULE MOD_NitrifData
#endif
