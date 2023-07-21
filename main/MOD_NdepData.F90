#include <define.h>

#ifdef BGC
MODULE MOD_NdepData
 !-----------------------------------------------------------------------
 ! !DESCRIPTION:
 ! This module read in ndep data.
 !
 ! !ORIGINAL:
 ! Lu Xingjie and Zhang Shupeng, 2023, prepare the original version of the ndep data module.

   USE MOD_Grid
   USE MOD_Mapping_Grid2Pset
   use MOD_BGC_Vars_TimeVariables, only : ndep
   use MOD_BGC_Vars_1DFluxes, only: ndep_to_sminn
   IMPLICIT NONE
      
   CHARACTER(len=256) :: file_ndep

   TYPE(grid_type) :: grid_ndep
   type(mapping_grid2pset_type) :: mg2p_ndep

CONTAINS

   ! ----------
   SUBROUTINE init_ndep_data (YY)

   !----------------------
   ! DESCTIPTION:
   ! open ndep netcdf file from DEF_dir_runtime, read latitude and longitude info.
   ! Initialize ndep data read in.

      use MOD_TimeManager
      USE MOD_Namelist
      USE MOD_Grid
      USE MOD_NetCDFSerial
      USE MOD_LandPatch
      IMPLICIT NONE

      integer, intent(in) :: YY

      ! Local Variables
      REAL(r8), allocatable :: lat(:), lon(:)

      file_ndep = trim(DEF_dir_runtime) // '/ndep/fndep_colm_hist_simyr1849-2006_1.9x2.5_c100428.nc'

      CALL ncio_read_bcast_serial (file_ndep, 'lat', lat)
      CALL ncio_read_bcast_serial (file_ndep, 'lon', lon)

      CALL grid_ndep%define_by_center (lat, lon)

      call mg2p_ndep%build (grid_ndep, landpatch)

      IF (allocated(lon)) deallocate(lon)
      IF (allocated(lat)) deallocate(lat)

      CALL update_ndep_data (YY, iswrite = .true.)

   END SUBROUTINE init_ndep_data

   ! ----------
   SUBROUTINE update_ndep_data (YY, iswrite)
! ===========================================================
!
! !DESCRIPTION:
! Read in the Nitrogen deposition data from CLM5.
!
! !REFERENCE:
! Galloway, J.N., et al. 2004. Nitrogen cycles: past, present, and future. Biogeochem. 70:153-226.
!
! !ORIGINAL:
! Created by Xingjie Lu and Shupeng Zhang, 2022
! ===========================================================

      use MOD_SPMD_Task
      USE MOD_Namelist, only : DEF_USE_PN
      USE MOD_DataType
      USE MOD_NetCDFBlock
      use MOD_LandPatch
      use MOD_Vars_TimeInvariants
      USE MOD_RangeCheck
      IMPLICIT NONE

      integer, intent(in) :: YY
      logical, INTENT(in) :: iswrite

      ! Local Variables
      TYPE(block_data_real8_2d) :: f_xy_ndep
      integer :: itime, npatch, m

      itime = max(min(YY,2006),1849) - 1848

      IF (p_is_io) THEN
         CALL allocate_block_data  (grid_ndep, f_xy_ndep)
         CALL ncio_read_block_time (file_ndep, 'NDEP_year', grid_ndep, itime, f_xy_ndep)
      ENDIF

      call mg2p_ndep%map_aweighted (f_xy_ndep, ndep)

      if (p_is_worker .and. iswrite) then
         if (numpatch > 0) then
            do npatch = 1, numpatch
               m = patchclass(npatch)
               if(m == 0)then
                  ndep_to_sminn(npatch) = 0.
               else
                  if(DEF_USE_PN)then
                     ndep_to_sminn(npatch)  = ndep(npatch) / 3600. / 365. / 24. * 5
                  else
                     ndep_to_sminn(npatch)  = ndep(npatch) / 3600. / 365. / 24.
                  end if
               end if
            end do

         ENDIF
      ENDIF

#ifdef RangeCheck
      call check_vector_data ('ndep', ndep)
#endif

   END SUBROUTINE update_ndep_data

END MODULE MOD_NdepData
#endif
