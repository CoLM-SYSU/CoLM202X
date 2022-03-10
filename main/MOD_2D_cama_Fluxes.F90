#include <define.h>

MODULE MOD_2D_cama_Fluxes
! ----------------------------------------------------------------------
! perfrom the grid average mapping: average a subgrid input 1d vector 
! of length numpatch to a output 2d array of length [lon_points,lat_points]
!
! Created by Yongjiu Dai, 03/2014
!---------------------------------------------------------------------
   use mod_data_type
   USE GlobalVars

   IMPLICIT NONE
   SAVE
   type(block_data_real8_2d) :: f_rnof_cama    ! total runoff [mm/s]
   type(block_data_real8_2d) :: IO_Effdepth    ! inundation to water depth [m]
   type(block_data_real8_2d) :: IO_Effarea
   ! PUBLIC MEMBER FUNCTIONS:
   public :: allocate_2D_cama_Fluxes

CONTAINS

   SUBROUTINE allocate_2D_cama_Fluxes (grid)
      ! --------------------------------------------------------------------
      ! Allocates memory for CLM 2d [lon_points,lat_points] variables
      ! --------------------------------------------------------------------

      use spmd_task
      use mod_grid
      use mod_data_type
      implicit none

      type(grid_type), intent(in) :: grid

      if (p_is_io) then
         call allocate_block_data (grid, f_rnof_cama)  ! total runoff [mm/s]
         call allocate_block_data (grid, IO_Effdepth)  ! inundation to water depth [m]
         call allocate_block_data (grid, IO_Effarea)  ! inundation to water depth [m]
      end if

   END SUBROUTINE allocate_2D_cama_Fluxes

END MODULE MOD_2D_cama_Fluxes
