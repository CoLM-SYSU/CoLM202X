#include <define.h>

MODULE MOD_Vars_2DForcing
   ! -------------------------------
   ! Meteorogical Forcing
   !
   ! Created by Yongjiu Dai, 03/2014
   ! -------------------------------

   use mod_data_type
   IMPLICIT NONE
   SAVE

   ! -----------------------------------------------------------------
   type(block_data_real8_2d) :: forc_xy_pco2m  ! CO2 concentration in atmos. (pascals)
   type(block_data_real8_2d) :: forc_xy_po2m   ! O2 concentration in atmos. (pascals)
   type(block_data_real8_2d) :: forc_xy_us     ! wind in eastward direction [m/s]
   type(block_data_real8_2d) :: forc_xy_vs     ! wind in northward direction [m/s]
   type(block_data_real8_2d) :: forc_xy_t      ! temperature at reference height [kelvin]
   type(block_data_real8_2d) :: forc_xy_q      ! specific humidity at reference height [kg/kg]
   type(block_data_real8_2d) :: forc_xy_prc    ! convective precipitation [mm/s]
   type(block_data_real8_2d) :: forc_xy_prl    ! large scale precipitation [mm/s]
   type(block_data_real8_2d) :: forc_xy_psrf   ! atmospheric pressure at the surface [pa]
   type(block_data_real8_2d) :: forc_xy_pbot   ! atm bottom level pressure (or reference height) (pa)
   type(block_data_real8_2d) :: forc_xy_sols   ! atm vis direct beam solar rad onto srf [W/m2]
   type(block_data_real8_2d) :: forc_xy_soll   ! atm nir direct beam solar rad onto srf [W/m2]
   type(block_data_real8_2d) :: forc_xy_solsd  ! atm vis diffuse solar rad onto srf [W/m2]
   type(block_data_real8_2d) :: forc_xy_solld  ! atm nir diffuse solar rad onto srf [W/m2]
   type(block_data_real8_2d) :: forc_xy_frl    ! atmospheric infrared (longwave) radiation [W/m2]
   type(block_data_real8_2d) :: forc_xy_hgt_u  ! observational height of wind [m]
   type(block_data_real8_2d) :: forc_xy_hgt_t  ! observational height of temperature [m]
   type(block_data_real8_2d) :: forc_xy_hgt_q  ! observational height of humidity [m]
   type(block_data_real8_2d) :: forc_xy_rhoair ! air density [kg/m3]
   type(block_data_real8_2d) :: forc_xy_hpbl   ! atmospheric boundary layer height [m]

   ! PUBLIC MEMBER FUNCTIONS:
   public :: allocate_2D_Forcing

CONTAINS

   !-----------------------------------------------------------------------

   SUBROUTINE allocate_2D_Forcing (grid)
      ! ------------------------------------------------
      ! Allocates memory for CoLM 2d [lon_points,lat_points] variables
      ! ------------------------------------------------
      use spmd_task
      use mod_grid
      use mod_data_type
      IMPLICIT NONE

      type(grid_type), intent(in) :: grid

      if (p_is_io) then

         call allocate_block_data (grid, forc_xy_pco2m ) ! CO2 concentration in atmos. (pascals)
         call allocate_block_data (grid, forc_xy_po2m  ) ! O2 concentration in atmos. (pascals)
         call allocate_block_data (grid, forc_xy_us    ) ! wind in eastward direction [m/s]
         call allocate_block_data (grid, forc_xy_vs    ) ! wind in northward direction [m/s]
         call allocate_block_data (grid, forc_xy_t     ) ! temperature at reference height [kelvin]
         call allocate_block_data (grid, forc_xy_q     ) ! specific humidity at reference height [kg/kg]
         call allocate_block_data (grid, forc_xy_prc   ) ! convective precipitation [mm/s]
         call allocate_block_data (grid, forc_xy_prl   ) ! large scale precipitation [mm/s]
         call allocate_block_data (grid, forc_xy_psrf  ) ! atmospheric pressure at the surface [pa]
         call allocate_block_data (grid, forc_xy_pbot  ) ! atm bottom level pressure (or reference height) (pa)
         call allocate_block_data (grid, forc_xy_sols  ) ! atm vis direct beam solar rad onto srf [W/m2]
         call allocate_block_data (grid, forc_xy_soll  ) ! atm nir direct beam solar rad onto srf [W/m2]
         call allocate_block_data (grid, forc_xy_solsd ) ! atm vis diffuse solar rad onto srf [W/m2]
         call allocate_block_data (grid, forc_xy_solld ) ! atm nir diffuse solar rad onto srf [W/m2]
         call allocate_block_data (grid, forc_xy_frl   ) ! atmospheric infrared (longwave) radiation [W/m2]
         call allocate_block_data (grid, forc_xy_hgt_u ) ! observational height of wind [m]
         call allocate_block_data (grid, forc_xy_hgt_t ) ! observational height of temperature [m]
         call allocate_block_data (grid, forc_xy_hgt_q ) ! observational height of humidity [m]
         call allocate_block_data (grid, forc_xy_rhoair) ! air density [kg/m3]
         call allocate_block_data (grid, forc_xy_hpbl  ) ! atmospheric boundary layer height [m]
      end if

   END SUBROUTINE allocate_2D_Forcing

END MODULE MOD_Vars_2DForcing
! ------ EOP --------
