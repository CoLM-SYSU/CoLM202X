#include <define.h>

MODULE MOD_Urban_Vars_2DFluxes

   use MOD_DataType

   type(block_data_real8_2d) :: f_t_room    ! temperature of inner building [K]
   type(block_data_real8_2d) :: f_tafu      ! temperature of outer building [K]
   type(block_data_real8_2d) :: f_fhac      ! sensible flux from heat or cool AC [W/m2]
   type(block_data_real8_2d) :: f_fwst      ! waste heat flux from heat or cool AC [W/m2]
   type(block_data_real8_2d) :: f_fach      ! flux from inner and outter air exchange [W/m2]
   type(block_data_real8_2d) :: f_fahe      ! flux from metabolism and vehicle [W/m2]
   type(block_data_real8_2d) :: f_fhah      ! sensible flux from heating [W/m2]
   type(block_data_real8_2d) :: f_fvehc     ! flux from vehicle [W/m2]
   type(block_data_real8_2d) :: f_fmeta     ! flux from metabolism [W/m2]

   type(block_data_real8_2d) :: f_senroof   ! sensible heat flux from roof [W/m2]
   type(block_data_real8_2d) :: f_senwsun   ! sensible heat flux from sunlit wall [W/m2]
   type(block_data_real8_2d) :: f_senwsha   ! sensible heat flux from shaded wall [W/m2]
   type(block_data_real8_2d) :: f_sengimp   ! sensible heat flux from impervious road [W/m2]
   type(block_data_real8_2d) :: f_sengper   ! sensible heat flux from pervious road [W/m2]
   type(block_data_real8_2d) :: f_senurbl   ! sensible heat flux from urban vegetation [W/m2]

   type(block_data_real8_2d) :: f_lfevproof ! latent heat flux from roof [W/m2]
   type(block_data_real8_2d) :: f_lfevpgimp ! latent heat flux from impervious road [W/m2]
   type(block_data_real8_2d) :: f_lfevpgper ! latent heat flux from pervious road [W/m2]
   type(block_data_real8_2d) :: f_lfevpurbl ! latent heat flux from urban vegetation [W/m2]

   type(block_data_real8_2d) :: f_troof     ! temperature of roof [K]
   type(block_data_real8_2d) :: f_twall     ! temperature of wall [K]

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: allocate_Urban_2DFluxes

CONTAINS

   SUBROUTINE allocate_Urban_2DFluxes(grid)

      use MOD_SPMD_Task
      use MOD_Grid
      use MOD_DataType
      implicit none

      type(grid_type), intent(in) :: grid

      IF (p_is_io) THEN
         call allocate_block_data (grid, f_t_room   ) ! temperature of inner building
         call allocate_block_data (grid, f_tafu     ) ! temperature of outer building
         call allocate_block_data (grid, f_fhac     ) ! sensible flux from heat
         call allocate_block_data (grid, f_fwst     ) ! waste heat flux from hea
         call allocate_block_data (grid, f_fach     ) ! flux from inner and outt
         call allocate_block_data (grid, f_fahe     ) ! flux from metabolism and
         call allocate_block_data (grid, f_fhah     ) ! sensible flux from heati
         call allocate_block_data (grid, f_fvehc    ) ! flux from vehicle [W/m2]
         call allocate_block_data (grid, f_fmeta    ) ! flux from metabolism [W/
         call allocate_block_data (grid, f_senroof  ) ! sensible heat flux from
         call allocate_block_data (grid, f_senwsun  ) ! sensible heat flux from
         call allocate_block_data (grid, f_senwsha  ) ! sensible heat flux from
         call allocate_block_data (grid, f_sengimp  ) ! sensible heat flux from
         call allocate_block_data (grid, f_sengper  ) ! sensible heat flux from
         call allocate_block_data (grid, f_senurbl  ) ! sensible heat flux from
         call allocate_block_data (grid, f_lfevproof) ! latent heat flux from ro
         call allocate_block_data (grid, f_lfevpgimp) ! latent heat flux from im
         call allocate_block_data (grid, f_lfevpgper) ! latent heat flux from pe
         call allocate_block_data (grid, f_lfevpurbl) ! latent heat flux from ur
         call allocate_block_data (grid, f_troof    ) ! temperature of roof [K]
         call allocate_block_data (grid, f_twall    ) ! temperature of wall [K]
      ENDIF

   END SUBROUTINE allocate_Urban_2DFluxes
END MODULE MOD_Urban_Vars_2DFluxes
