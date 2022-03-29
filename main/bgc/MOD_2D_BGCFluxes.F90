#include <define.h>

MODULE MOD_2D_BGCFluxes
! ----------------------------------------------------------------------
! perfrom the grid average mapping: average a subgrid input 1d vector 
! of length numpatch to a output 2d array of length [lon_points,lat_points]
!
! Created by Yongjiu Dai, 03/2014
!---------------------------------------------------------------------
#ifdef BGC

   use mod_data_type
   USE GlobalVars

   IMPLICIT NONE
   SAVE

   type(block_data_real8_2d) :: f_leafc              ! leaf carbon display pool  (gC/m2)
   type(block_data_real8_2d) :: f_leafc_storage      ! leaf carbon storage pool  (gC/m2)
   type(block_data_real8_2d) :: f_leafc_xfer         ! leaf carbon transfer pool (gC/m2)
   type(block_data_real8_2d) :: f_frootc             ! fine root carbon display pool  (gC/m2)
   type(block_data_real8_2d) :: f_frootc_storage     ! fine root carbon storage pool  (gC/m2)
   type(block_data_real8_2d) :: f_frootc_xfer        ! fine root carbon transfer pool (gC/m2)
   type(block_data_real8_2d) :: f_livestemc          ! live stem carbon display pool  (gC/m2)
   type(block_data_real8_2d) :: f_livestemc_storage  ! live stem carbon storage pool  (gC/m2)
   type(block_data_real8_2d) :: f_livestemc_xfer     ! live stem carbon transfer pool (gC/m2)
   type(block_data_real8_2d) :: f_deadstemc          ! dead stem carbon display pool  (gC/m2)
   type(block_data_real8_2d) :: f_deadstemc_storage  ! dead stem carbon storage pool  (gC/m2)
   type(block_data_real8_2d) :: f_deadstemc_xfer     ! dead stem carbon transfer pool (gC/m2)
   type(block_data_real8_2d) :: f_livecrootc         ! live coarse root carbon display pool  (gC/m2)
   type(block_data_real8_2d) :: f_livecrootc_storage ! live coarse root carbon storage pool  (gC/m2)
   type(block_data_real8_2d) :: f_livecrootc_xfer    ! live coarse root carbon transfer pool (gC/m2)
   type(block_data_real8_2d) :: f_deadcrootc         ! dead coarse root carbon display pool  (gC/m2)
   type(block_data_real8_2d) :: f_deadcrootc_storage ! dead coarse root carbon storage pool  (gC/m2)
   type(block_data_real8_2d) :: f_deadcrootc_xfer    ! dead coarse root carbon transfer pool (gC/m2)
   type(block_data_real8_2d) :: f_leafn              ! leaf nitrogen display pool  (gC/m2)
   type(block_data_real8_2d) :: f_leafn_storage      ! leaf nitrogen storage pool  (gC/m2)
   type(block_data_real8_2d) :: f_leafn_xfer         ! leaf nitrogen transfer pool (gC/m2)
   type(block_data_real8_2d) :: f_frootn             ! fine root nitrogen display pool  (gC/m2)
   type(block_data_real8_2d) :: f_frootn_storage     ! fine root nitrogen storage pool  (gC/m2)
   type(block_data_real8_2d) :: f_frootn_xfer        ! fine root nitrogen transfer pool (gC/m2)
   type(block_data_real8_2d) :: f_livestemn          ! live stem nitrogen display pool  (gC/m2)
   type(block_data_real8_2d) :: f_livestemn_storage  ! live stem nitrogen storage pool  (gC/m2)
   type(block_data_real8_2d) :: f_livestemn_xfer     ! live stem nitrogen transfer pool (gC/m2)
   type(block_data_real8_2d) :: f_deadstemn          ! dead stem nitrogen display pool  (gC/m2)
   type(block_data_real8_2d) :: f_deadstemn_storage  ! dead stem nitrogen storage pool  (gC/m2)
   type(block_data_real8_2d) :: f_deadstemn_xfer     ! dead stem nitrogen transfer pool (gC/m2)
   type(block_data_real8_2d) :: f_livecrootn         ! live coarse root nitrogen display pool  (gC/m2)
   type(block_data_real8_2d) :: f_livecrootn_storage ! live coarse root nitrogen storage pool  (gC/m2)
   type(block_data_real8_2d) :: f_livecrootn_xfer    ! live coarse root nitrogen transfer pool (gC/m2)
   type(block_data_real8_2d) :: f_deadcrootn         ! dead coarse root nitrogen display pool  (gC/m2)
   type(block_data_real8_2d) :: f_deadcrootn_storage ! dead coarse root nitrogen storage pool  (gC/m2)
   type(block_data_real8_2d) :: f_deadcrootn_xfer    ! dead coarse root nitrogen transfer pool (gC/m2)

   type(block_data_real8_2d) :: f_grainc             ! grain carbon display pool (gC/m2)
   type(block_data_real8_2d) :: f_grainc_storage     ! grain carbon storage pool (gC/m2)
   type(block_data_real8_2d) :: f_grainc_xfer        ! grain carbon transfer pool (gC/m2)

   type(block_data_real8_2d) :: f_cphase             ! crop phase
   type(block_data_real8_2d) :: f_cropprod1c         ! 1-yr crop production carbon
   type(block_data_real8_2d) :: f_cropprod1c_loss    ! loss of 1-yr crop production carbon
   type(block_data_real8_2d) :: f_cropseedc_deficit  ! crop seed carbon deficit
   type(block_data_real8_2d) :: f_grainc_to_cropprodc! grain to crop production
   type(block_data_real8_2d) :: f_grainc_to_seed     ! grain to crop seed

   type(block_data_real8_2d) :: f_gpp                ! net primary production (gC/m2)
   type(block_data_real8_2d) :: f_downreg            ! gpp downregulation due to N limitation

   type(block_data_real8_2d) :: f_litr1c             ! soil carbon pool [gC/m2]
   type(block_data_real8_2d) :: f_litr2c             ! soil carbon pool [gC/m2]
   type(block_data_real8_2d) :: f_litr3c             ! soil carbon pool [gC/m2]
   type(block_data_real8_2d) :: f_cwdc               ! soil carbon pool [gC/m2]
   type(block_data_real8_2d) :: f_soil1c             ! soil carbon pool [gC/m2]
   type(block_data_real8_2d) :: f_soil2c             ! soil carbon pool [gC/m2]
   type(block_data_real8_2d) :: f_soil3c             ! soil carbon pool [gC/m2]

   type(block_data_real8_2d) :: f_litr1n             ! soil nitrogen pool [gN/m2]
   type(block_data_real8_2d) :: f_litr2n             ! soil nitrogen pool [gN/m2]
   type(block_data_real8_2d) :: f_litr3n             ! soil nitrogen pool [gN/m2]
   type(block_data_real8_2d) :: f_cwdn               ! soil nitrogen pool [gN/m2]
   type(block_data_real8_2d) :: f_soil1n             ! soil nitrogen pool [gN/m2]
   type(block_data_real8_2d) :: f_soil2n             ! soil nitrogen pool [gN/m2]
   type(block_data_real8_2d) :: f_soil3n             ! soil nitrogen pool [gN/m2]
   type(block_data_real8_2d) :: f_sminn              ! soil mineral nitrogen pool [gN/m2]

   ! PUBLIC MEMBER FUNCTIONS:
   public :: allocate_2D_BGCFluxes

CONTAINS

   SUBROUTINE allocate_2D_BGCFluxes (grid)
      ! --------------------------------------------------------------------
      ! Allocates memory for CLM 2d [lon_points,lat_points] variables
      ! --------------------------------------------------------------------

      use spmd_task
      use mod_grid
      use mod_data_type
      implicit none

      type(grid_type), intent(in) :: grid

      if (p_is_io) then
         
         call allocate_block_data (grid, f_leafc              ) ! leaf carbon display pool  (gC/m2)
         call allocate_block_data (grid, f_leafc_storage      ) ! leaf carbon storage pool  (gC/m2)
         call allocate_block_data (grid, f_leafc_xfer         ) ! leaf carbon transfer pool (gC/m2)
         call allocate_block_data (grid, f_frootc             ) ! fine root carbon display pool  (gC/m2)
         call allocate_block_data (grid, f_frootc_storage     ) ! fine root carbon storage pool  (gC/m2)
         call allocate_block_data (grid, f_frootc_xfer        ) ! fine root carbon transfer pool (gC/m2)
         call allocate_block_data (grid, f_livestemc          ) ! live stem carbon display pool  (gC/m2)
         call allocate_block_data (grid, f_livestemc_storage  ) ! live stem carbon storage pool  (gC/m2)
         call allocate_block_data (grid, f_livestemc_xfer     ) ! live stem carbon transfer pool (gC/m2)
         call allocate_block_data (grid, f_deadstemc          ) ! dead stem carbon display pool  (gC/m2)
         call allocate_block_data (grid, f_deadstemc_storage  ) ! dead stem carbon storage pool  (gC/m2)
         call allocate_block_data (grid, f_deadstemc_xfer     ) ! dead stem carbon transfer pool (gC/m2)
         call allocate_block_data (grid, f_livecrootc         ) ! live coarse root carbon display pool  (gC/m2)
         call allocate_block_data (grid, f_livecrootc_storage ) ! live coarse root carbon storage pool  (gC/m2)
         call allocate_block_data (grid, f_livecrootc_xfer    ) ! live coarse root carbon transfer pool (gC/m2)
         call allocate_block_data (grid, f_deadcrootc         ) ! dead coarse root carbon display pool  (gC/m2)
         call allocate_block_data (grid, f_deadcrootc_storage ) ! dead coarse root carbon storage pool  (gC/m2)
         call allocate_block_data (grid, f_deadcrootc_xfer    ) ! dead coarse root carbon transfer pool (gC/m2)
         call allocate_block_data (grid, f_leafn              ) ! leaf nitrogen display pool  (gC/m2)
         call allocate_block_data (grid, f_leafn_storage      ) ! leaf nitrogen storage pool  (gC/m2)
         call allocate_block_data (grid, f_leafn_xfer         ) ! leaf nitrogen transfer pool (gC/m2)
         call allocate_block_data (grid, f_frootn             ) ! fine root nitrogen display pool  (gC/m2)
         call allocate_block_data (grid, f_frootn_storage     ) ! fine root nitrogen storage pool  (gC/m2)
         call allocate_block_data (grid, f_frootn_xfer        ) ! fine root nitrogen transfer pool (gC/m2)
         call allocate_block_data (grid, f_livestemn          ) ! live stem nitrogen display pool  (gC/m2)
         call allocate_block_data (grid, f_livestemn_storage  ) ! live stem nitrogen storage pool  (gC/m2)
         call allocate_block_data (grid, f_livestemn_xfer     ) ! live stem nitrogen transfer pool (gC/m2)
         call allocate_block_data (grid, f_deadstemn          ) ! dead stem nitrogen display pool  (gC/m2)
         call allocate_block_data (grid, f_deadstemn_storage  ) ! dead stem nitrogen storage pool  (gC/m2)
         call allocate_block_data (grid, f_deadstemn_xfer     ) ! dead stem nitrogen transfer pool (gC/m2)
         call allocate_block_data (grid, f_livecrootn         ) ! live coarse root nitrogen display pool  (gC/m2)
         call allocate_block_data (grid, f_livecrootn_storage ) ! live coarse root nitrogen storage pool  (gC/m2)
         call allocate_block_data (grid, f_livecrootn_xfer    ) ! live coarse root nitrogen transfer pool (gC/m2)
         call allocate_block_data (grid, f_deadcrootn         ) ! dead coarse root nitrogen display pool  (gC/m2)
         call allocate_block_data (grid, f_deadcrootn_storage ) ! dead coarse root nitrogen storage pool  (gC/m2)
         call allocate_block_data (grid, f_deadcrootn_xfer    ) ! dead coarse root nitrogen transfer pool (gC/m2)

         call allocate_block_data (grid, f_grainc             ) ! grain carbon display pool  (gC/m2)
         call allocate_block_data (grid, f_grainc_storage     ) ! grain carbon storage pool  (gC/m2)
         call allocate_block_data (grid, f_grainc_xfer        ) ! grain carbon transfer pool (gC/m2)

         call allocate_block_data (grid, f_cphase             )  ! crop phase
         call allocate_block_data (grid, f_cropprod1c         )  ! 1-yr crop production carbon
         call allocate_block_data (grid, f_cropprod1c_loss    )  ! loss of 1-yr crop production carbon
         call allocate_block_data (grid, f_cropseedc_deficit  )  ! crop seed carbon deficit
         call allocate_block_data (grid, f_grainc_to_cropprodc ) ! grain to crop production
         call allocate_block_data (grid, f_grainc_to_seed     )  ! grain to crop seed

         call allocate_block_data (grid, f_gpp                ) ! net primary production (gC/m2)
         call allocate_block_data (grid, f_downreg            ) ! gpp downregulation due to N limitation

         call allocate_block_data (grid, f_litr1c             )  ! soil carbon pool (gC/m2)
         call allocate_block_data (grid, f_litr2c             )  ! soil carbon pool (gC/m2)
         call allocate_block_data (grid, f_litr3c             )  ! soil carbon pool (gC/m2)
         call allocate_block_data (grid, f_cwdc               )  ! soil carbon pool (gC/m2)
         call allocate_block_data (grid, f_soil1c             )  ! soil carbon pool (gC/m2)
         call allocate_block_data (grid, f_soil2c             )  ! soil carbon pool (gC/m2)
         call allocate_block_data (grid, f_soil3c             )  ! soil carbon pool (gC/m2)

         call allocate_block_data (grid, f_litr1n             )  ! soil nitrogen pool (gN/m2)
         call allocate_block_data (grid, f_litr2n             )  ! soil nitrogen pool (gN/m2)
         call allocate_block_data (grid, f_litr3n             )  ! soil nitrogen pool (gN/m2)
         call allocate_block_data (grid, f_cwdn               )  ! soil nitrogen pool (gN/m2)
         call allocate_block_data (grid, f_soil1n             )  ! soil nitrogen pool (gN/m2)
         call allocate_block_data (grid, f_soil2n             )  ! soil nitrogen pool (gN/m2)
         call allocate_block_data (grid, f_soil3n             )  ! soil nitrogen pool (gN/m2)
         call allocate_block_data (grid, f_sminn              )  ! soil mineral nitrogen pool (gN/m2)

      end if


   END SUBROUTINE allocate_2D_BGCFluxes

#endif
END MODULE MOD_2D_BGCFluxes
