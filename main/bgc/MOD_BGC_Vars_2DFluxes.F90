#include <define.h>

MODULE MOD_BGC_Vars_2DFluxes
! ----------------------------------------------------------------------
! !DESCRIPTION:
! perfrom the grid average mapping: average a subgrid input 1d vector
! of length numpatch to a output 2d array of length [lon_points,lat_points] for biogeochemical flux variables
!
! !ORIGINAL:
! Xingjie Lu, 2022, created the original version.
!---------------------------------------------------------------------
#ifdef BGC

   use mod_data_type
   USE GlobalVars

   IMPLICIT NONE
   SAVE

   type(block_data_real8_2d) :: f_leafc              ! 2D grid: leaf carbon display pool  (gC m-2)
   type(block_data_real8_2d) :: f_leafc_storage      ! 2D grid: leaf carbon storage pool  (gC m-2)
   type(block_data_real8_2d) :: f_leafc_xfer         ! 2D grid: leaf carbon transfer pool (gC m-2)
   type(block_data_real8_2d) :: f_frootc             ! 2D grid: fine root carbon display pool  (gC m-2)
   type(block_data_real8_2d) :: f_frootc_storage     ! 2D grid: fine root carbon storage pool  (gC m-2)
   type(block_data_real8_2d) :: f_frootc_xfer        ! 2D grid: fine root carbon transfer pool (gC m-2)
   type(block_data_real8_2d) :: f_livestemc          ! 2D grid: live stem carbon display pool  (gC m-2)
   type(block_data_real8_2d) :: f_livestemc_storage  ! 2D grid: live stem carbon storage pool  (gC m-2)
   type(block_data_real8_2d) :: f_livestemc_xfer     ! 2D grid: live stem carbon transfer pool (gC m-2)
   type(block_data_real8_2d) :: f_deadstemc          ! 2D grid: dead stem carbon display pool  (gC m-2)
   type(block_data_real8_2d) :: f_deadstemc_storage  ! 2D grid: dead stem carbon storage pool  (gC m-2)
   type(block_data_real8_2d) :: f_deadstemc_xfer     ! 2D grid: dead stem carbon transfer pool (gC m-2)
   type(block_data_real8_2d) :: f_livecrootc         ! 2D grid: live coarse root carbon display pool  (gC m-2)
   type(block_data_real8_2d) :: f_livecrootc_storage ! 2D grid: live coarse root carbon storage pool  (gC m-2)
   type(block_data_real8_2d) :: f_livecrootc_xfer    ! 2D grid: live coarse root carbon transfer pool (gC m-2)
   type(block_data_real8_2d) :: f_deadcrootc         ! 2D grid: dead coarse root carbon display pool  (gC m-2)
   type(block_data_real8_2d) :: f_deadcrootc_storage ! 2D grid: dead coarse root carbon storage pool  (gC m-2)
   type(block_data_real8_2d) :: f_deadcrootc_xfer    ! 2D grid: dead coarse root carbon transfer pool (gC m-2)
#ifdef CROP
   type(block_data_real8_2d) :: f_grainc             ! 2D grid: grain carbon display pool (gC m-2)
   type(block_data_real8_2d) :: f_grainc_storage     ! 2D grid: grain carbon storage pool (gC m-2)
   type(block_data_real8_2d) :: f_grainc_xfer        ! 2D grid: grain carbon transfer pool (gC m-2)
#endif
   type(block_data_real8_2d) :: f_leafn              ! 2D grid: leaf nitrogen display pool  (gN m-2)
   type(block_data_real8_2d) :: f_leafn_storage      ! 2D grid: leaf nitrogen storage pool  (gN m-2)
   type(block_data_real8_2d) :: f_leafn_xfer         ! 2D grid: leaf nitrogen transfer pool (gN m-2)
   type(block_data_real8_2d) :: f_frootn             ! 2D grid: fine root nitrogen display pool  (gN m-2)
   type(block_data_real8_2d) :: f_frootn_storage     ! 2D grid: fine root nitrogen storage pool  (gN m-2)
   type(block_data_real8_2d) :: f_frootn_xfer        ! 2D grid: fine root nitrogen transfer pool (gN m-2)
   type(block_data_real8_2d) :: f_livestemn          ! 2D grid: live stem nitrogen display pool  (gN m-2)
   type(block_data_real8_2d) :: f_livestemn_storage  ! 2D grid: live stem nitrogen storage pool  (gN m-2)
   type(block_data_real8_2d) :: f_livestemn_xfer     ! 2D grid: live stem nitrogen transfer pool (gN m-2)
   type(block_data_real8_2d) :: f_deadstemn          ! 2D grid: dead stem nitrogen display pool  (gN m-2)
   type(block_data_real8_2d) :: f_deadstemn_storage  ! 2D grid: dead stem nitrogen storage pool  (gN m-2)
   type(block_data_real8_2d) :: f_deadstemn_xfer     ! 2D grid: dead stem nitrogen transfer pool (gN m-2)
   type(block_data_real8_2d) :: f_livecrootn         ! 2D grid: live coarse root nitrogen display pool  (gN m-2)
   type(block_data_real8_2d) :: f_livecrootn_storage ! 2D grid: live coarse root nitrogen storage pool  (gN m-2)
   type(block_data_real8_2d) :: f_livecrootn_xfer    ! 2D grid: live coarse root nitrogen transfer pool (gN m-2)
   type(block_data_real8_2d) :: f_deadcrootn         ! 2D grid: dead coarse root nitrogen display pool  (gN m-2)
   type(block_data_real8_2d) :: f_deadcrootn_storage ! 2D grid: dead coarse root nitrogen storage pool  (gN m-2)
   type(block_data_real8_2d) :: f_deadcrootn_xfer    ! 2D grid: dead coarse root nitrogen transfer pool (gN m-2)

#ifdef CROP
   type(block_data_real8_2d) :: f_grainn             ! 2D grid: grain nitrogen display pool (gN m-2)
   type(block_data_real8_2d) :: f_grainn_storage     ! 2D grid: grain nitrogen storage pool (gN m-2)
   type(block_data_real8_2d) :: f_grainn_xfer        ! 2D grid: grain nitrogen transfer pool (gN m-2)
#endif
   type(block_data_real8_2d) :: f_retransn           ! 2D grid: retranslocation nitrogen pool (gN m-2)

#ifdef CROP
   type(block_data_real8_2d) :: f_cphase             ! 2D grid: crop phase
   type(block_data_real8_2d) :: f_cropprod1c         ! 2D grid: 1-yr crop production carbon (gC m-2)
   type(block_data_real8_2d) :: f_cropprod1c_loss    ! 2D grid: loss of 1-yr crop production carbon (gC m-2 s-1)
   type(block_data_real8_2d) :: f_cropseedc_deficit  ! 2D grid: crop seed carbon deficit (gC m-2 s-1)
   type(block_data_real8_2d) :: f_grainc_to_cropprodc! 2D grid: grain to crop production carbon (gC m-2 s-1)
   type(block_data_real8_2d) :: f_grainc_to_seed     ! 2D grid: grain to crop seed carbon (gC m-2 s-1)
   type(block_data_real8_2d) :: f_fert_to_sminn      ! 2D grid: fertilization (gN m-2 s-1)
   type(block_data_real8_2d) :: f_plantdate          ! 2D grid: planting date
#endif
   type(block_data_real8_2d) :: f_ndep_to_sminn      ! 2D grid: nitrogen deposition (gN m-2 s-1)

   type(block_data_real8_2d) :: f_gpp                ! 2D grid: net primary production (gC m-2 s-1)
   type(block_data_real8_2d) :: f_downreg            ! 2D grid: gpp downregulation due to N limitation
   type(block_data_real8_2d) :: f_ar                 ! 2D grid: autotrophic respiration (gC m-2 s-1)
   type(block_data_real8_2d) :: f_cwdprod            ! 2D grid: CWD production (gC m-2 s-1)
   type(block_data_real8_2d) :: f_cwddecomp          ! 2D grid: CWD decomposition (gC m-2 s-1)
   type(block_data_real8_2d) :: f_hr                 ! 2D grid: heterotrophic respiration (gC m-2 s-1)
   type(block_data_real8_2d) :: f_fpg                ! 2D grid: fraction of gpp potential
   type(block_data_real8_2d) :: f_fpi                ! 2D grid: fraction of immobalization
   type(block_data_real8_2d) :: f_gpp_enftemp        ! 2D grid: gross primary productivity for needleleaf evergreen temperate tree (gC m-2 s-1)
   type(block_data_real8_2d) :: f_gpp_enfboreal      ! 2D grid: gross primary productivity for needleleaf evergreen boreal tree (gC m-2 s-1)
   type(block_data_real8_2d) :: f_gpp_dnfboreal      ! 2D grid: gross primary productivity for needleleaf deciduous boreal tree (gC m-2 s-1)
   type(block_data_real8_2d) :: f_gpp_ebftrop        ! 2D grid: gross primary productivity for broadleaf evergreen tropical tree (gC m-2 s-1)
   type(block_data_real8_2d) :: f_gpp_ebftemp        ! 2D grid: gross primary productivity for broadleaf evergreen temperate tree (gC m-2 s-1)
   type(block_data_real8_2d) :: f_gpp_dbftrop        ! 2D grid: gross primary productivity for broadleaf deciduous tropical tree (gC m-2 s-1)
   type(block_data_real8_2d) :: f_gpp_dbftemp        ! 2D grid: gross primary productivity for broadleaf deciduous temperate tree (gC m-2 s-1)
   type(block_data_real8_2d) :: f_gpp_dbfboreal      ! 2D grid: gross primary productivity for broadleaf deciduous boreal tree (gC m-2 s-1)
   type(block_data_real8_2d) :: f_gpp_ebstemp        ! 2D grid: gross primary productivity for broadleaf evergreen temperate shrub (gC m-2 s-1)
   type(block_data_real8_2d) :: f_gpp_dbstemp        ! 2D grid: gross primary productivity for broadleaf deciduous temperate shrub (gC m-2 s-1)
   type(block_data_real8_2d) :: f_gpp_dbsboreal      ! 2D grid: gross primary productivity for broadleaf deciduous boreal shrub (gC m-2 s-1)
   type(block_data_real8_2d) :: f_gpp_c3arcgrass     ! 2D grid: gross primary productivity for c3 arctic grass (gC m-2 s-1)
   type(block_data_real8_2d) :: f_gpp_c3grass        ! 2D grid: gross primary productivity for c3 grass (gC m-2 s-1)
   type(block_data_real8_2d) :: f_gpp_c4grass        ! 2D grid: gross primary productivity for c4 grass (gC m-2 s-1)
   type(block_data_real8_2d) :: f_leafc_enftemp      ! 2D grid: leaf carbon display pool for needleleaf evergreen temperate tree (gC m-2)
   type(block_data_real8_2d) :: f_leafc_enfboreal    ! 2D grid: leaf carbon display pool for needleleaf evergreen boreal tree (gC m-2)
   type(block_data_real8_2d) :: f_leafc_dnfboreal    ! 2D grid: leaf carbon display pool for needleleaf deciduous boreal tree (gC m-2)
   type(block_data_real8_2d) :: f_leafc_ebftrop      ! 2D grid: leaf carbon display pool for broadleaf evergreen tropical tree (gC m-2)
   type(block_data_real8_2d) :: f_leafc_ebftemp      ! 2D grid: leaf carbon display pool for broadleaf evergreen temperate tree (gC m-2)
   type(block_data_real8_2d) :: f_leafc_dbftrop      ! 2D grid: leaf carbon display pool for broadleaf deciduous tropical tree (gC m-2)
   type(block_data_real8_2d) :: f_leafc_dbftemp      ! 2D grid: leaf carbon display pool for broadleaf deciduous temperate tree (gC m-2)
   type(block_data_real8_2d) :: f_leafc_dbfboreal    ! 2D grid: leaf carbon display pool for broadleaf deciduous boreal tree (gC m-2)
   type(block_data_real8_2d) :: f_leafc_ebstemp      ! 2D grid: leaf carbon display pool for broadleaf evergreen temperate shrub (gC m-2)
   type(block_data_real8_2d) :: f_leafc_dbstemp      ! 2D grid: leaf carbon display pool for broadleaf deciduous temperate shrub (gC m-2)
   type(block_data_real8_2d) :: f_leafc_dbsboreal    ! 2D grid: leaf carbon display pool for broadleaf deciduous boreal shrub (gC m-2)
   type(block_data_real8_2d) :: f_leafc_c3arcgrass   ! 2D grid: leaf carbon display pool for c3 arctic grass (gC m-2)
   type(block_data_real8_2d) :: f_leafc_c3grass      ! 2D grid: leaf carbon display pool for c3 grass (gC m-2)
   type(block_data_real8_2d) :: f_leafc_c4grass      ! 2D grid: leaf carbon display pool for c4 grass (gC m-2)

   type(block_data_real8_3d) :: f_litr1c_vr          ! 3D grid (vertical resolved): litter 1 (metabolic litter) carbon density in soil layers (gC m-3)
   type(block_data_real8_3d) :: f_litr2c_vr          ! 3D grid (vertical resolved): litter 2 (cellulosic litte) carbon density in soil layers (gC m-3)
   type(block_data_real8_3d) :: f_litr3c_vr          ! 3D grid (vertical resolved): litter 3 (lignin litter) carbon density in soil layers (gC m-3)
   type(block_data_real8_3d) :: f_cwdc_vr            ! 3d grid (vertical resolved): coarse woody debris carbon density in soil layers (gC m-3)
   type(block_data_real8_3d) :: f_soil1c_vr          ! 3D grid (vertical resolved): soil 1 (active soil organic matter) carbon density in soil layers (gC m-3)
   type(block_data_real8_3d) :: f_soil2c_vr          ! 3D grid (vertical resolved): soil 2 (slow soil organic matter) carbon density in soil layers (gC m-3)
   type(block_data_real8_3d) :: f_soil3c_vr          ! 3D grid (vertical resolved): soil 3 (passive soil organic matter) carbon density in soil layers (gC m-3)

   type(block_data_real8_3d) :: f_litr1n_vr          ! 3D grid (vertical resolved): litter 1 (metabolic litter) carbon density in soil layers (gN m-3)
   type(block_data_real8_3d) :: f_litr2n_vr          ! 3D grid (vertical resolved): litter 2 (cellulosic litte) carbon density in soil layers (gN m-3)
   type(block_data_real8_3d) :: f_litr3n_vr          ! 3D grid (vertical resolved): litter 3 (lignin litter) carbon density in soil layers (gN m-3)
   type(block_data_real8_3d) :: f_cwdn_vr            ! 3d grid (vertical resolved): coarse woody debris carbon density in soil layers (gN m-3)
   type(block_data_real8_3d) :: f_soil1n_vr          ! 3D grid (vertical resolved): soil 1 (active soil organic matter) carbon density in soil layers (gN m-3)
   type(block_data_real8_3d) :: f_soil2n_vr          ! 3D grid (vertical resolved): soil 2 (slow soil organic matter) carbon density in soil layers (gN m-3)
   type(block_data_real8_3d) :: f_soil3n_vr          ! 3D grid (vertical resolved): soil 3 (passive soil organic matter) carbon density in soil layers (gN m-3)
   type(block_data_real8_3d) :: f_sminn_vr           ! 3D grid (vertical resolved): mineral nitrogen density in soil layers (gN m-3)

#ifdef NITRIF
   type(block_data_real8_3d) :: f_O2_DECOMP_DEPTH_UNSAT  ! 3D grid (vertical resolved): O2 consumption from heterotrophic respiration and autotrophic respiration for non-inundated area (mol m-3 s-1)
   type(block_data_real8_3d) :: f_CONC_O2_UNSAT          ! 3D grid (vertical resolved): O2 soil Concentration for non-inundated area (mol m-3)
#endif
#ifdef CROP
   type(block_data_real8_2d) :: f_hui                    ! 2D grid: gdd since planting (degree-days)
   type(block_data_real8_2d) :: f_vf                     ! 2D grid: vernalization factor for cereal
   type(block_data_real8_2d) :: f_gddplant               ! 2D grid: gdd since planting (degree-days)
   type(block_data_real8_2d) :: f_gddmaturity            ! 2D grid: gdd needed to harvest (degree-days)
   type(block_data_real8_2d) :: f_pdcorn                 ! 2D grid: planting date of corn (day of year)
   type(block_data_real8_2d) :: f_pdswheat               ! 2D grid: planting date of spring wheat (day of year)
   type(block_data_real8_2d) :: f_pdwwheat               ! 2D grid: planting date of winter wheat (day of year)
   type(block_data_real8_2d) :: f_pdsoybean              ! 2D grid: planting date of soybean (day of year)
   type(block_data_real8_2d) :: f_pdcotton               ! 2D grid: planting date of cotton (day of year)
   type(block_data_real8_2d) :: f_pdrice1                ! 2D grid: planting date of rice1 (day of year)
   type(block_data_real8_2d) :: f_pdrice2                ! 2D grid: planting date of rice2 (day of year)
   type(block_data_real8_2d) :: f_pdsugarcane            ! 2D grid: planting date of sugarcane (day of year)
   type(block_data_real8_2d) :: f_fertnitro_corn         ! 2D grid: nitrogen fertilizer for corn (gN m-2 year-1)
   type(block_data_real8_2d) :: f_fertnitro_swheat       ! 2D grid: nitrogen fertilizer for spring wheat (gN m-2 year-1)
   type(block_data_real8_2d) :: f_fertnitro_wwheat       ! 2D grid: nitrogen fertilizer for winter wheat (gN m-2 year-1)
   type(block_data_real8_2d) :: f_fertnitro_soybean      ! 2D grid: nitrogen fertilizer for soybean (gN m-2 year-1)
   type(block_data_real8_2d) :: f_fertnitro_cotton       ! 2D grid: nitrogen fertilizer for cotton (gN m-2 year-1)
   type(block_data_real8_2d) :: f_fertnitro_rice1        ! 2D grid: nitrogen fertilizer for rice1 (gN m-2 year-1)
   type(block_data_real8_2d) :: f_fertnitro_rice2        ! 2D grid: nitrogen fertilizer for rice2 (gN m-2 year-1)
   type(block_data_real8_2d) :: f_fertnitro_sugarcane    ! 2D grid: nitrogen fertilizer for sugarcane (gN m-2 year-1)
#endif
#ifdef Fire
   type(block_data_real8_2d) :: f_abm                    ! 2D grid: peak crop fire month
   type(block_data_real8_2d) :: f_gdp                    ! 2D grid: global gdp
   type(block_data_real8_2d) :: f_peatf                  ! 2D grid: global peatland fraction data (0-1)
   type(block_data_real8_2d) :: f_hdm                    ! 2D grid: human population density (counts km-2)
   type(block_data_real8_2d) :: f_lnfm                   ! 2D grid: lightning frequency (counts km-2 hr-1)
#endif

   ! PUBLIC MEMBER FUNCTIONS:
   public :: allocate_2D_BGCFluxes

CONTAINS

   SUBROUTINE allocate_2D_BGCFluxes (grid)
      ! --------------------------------------------------------------------
      ! Allocates memory for CoLM 2d [lon_points,lat_points] variables
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
#ifdef CROP
         call allocate_block_data (grid, f_grainc             ) ! grain carbon display pool  (gC/m2)
         call allocate_block_data (grid, f_grainc_storage     ) ! grain carbon storage pool  (gC/m2)
         call allocate_block_data (grid, f_grainc_xfer        ) ! grain carbon transfer pool (gC/m2)
#endif
         call allocate_block_data (grid, f_leafn              ) ! leaf nitrogen display pool  (gN/m2)
         call allocate_block_data (grid, f_leafn_storage      ) ! leaf nitrogen storage pool  (gN/m2)
         call allocate_block_data (grid, f_leafn_xfer         ) ! leaf nitrogen transfer pool (gN/m2)
         call allocate_block_data (grid, f_frootn             ) ! fine root nitrogen display pool  (gN/m2)
         call allocate_block_data (grid, f_frootn_storage     ) ! fine root nitrogen storage pool  (gN/m2)
         call allocate_block_data (grid, f_frootn_xfer        ) ! fine root nitrogen transfer pool (gN/m2)
         call allocate_block_data (grid, f_livestemn          ) ! live stem nitrogen display pool  (gN/m2)
         call allocate_block_data (grid, f_livestemn_storage  ) ! live stem nitrogen storage pool  (gN/m2)
         call allocate_block_data (grid, f_livestemn_xfer     ) ! live stem nitrogen transfer pool (gN/m2)
         call allocate_block_data (grid, f_deadstemn          ) ! dead stem nitrogen display pool  (gN/m2)
         call allocate_block_data (grid, f_deadstemn_storage  ) ! dead stem nitrogen storage pool  (gN/m2)
         call allocate_block_data (grid, f_deadstemn_xfer     ) ! dead stem nitrogen transfer pool (gN/m2)
         call allocate_block_data (grid, f_livecrootn         ) ! live coarse root nitrogen display pool  (gN/m2)
         call allocate_block_data (grid, f_livecrootn_storage ) ! live coarse root nitrogen storage pool  (gN/m2)
         call allocate_block_data (grid, f_livecrootn_xfer    ) ! live coarse root nitrogen transfer pool (gN/m2)
         call allocate_block_data (grid, f_deadcrootn         ) ! dead coarse root nitrogen display pool  (gN/m2)
         call allocate_block_data (grid, f_deadcrootn_storage ) ! dead coarse root nitrogen storage pool  (gN/m2)
         call allocate_block_data (grid, f_deadcrootn_xfer    ) ! dead coarse root nitrogen transfer pool (gN/m2)

#ifdef CROP
         call allocate_block_data (grid, f_grainn             ) ! grain nitrogen display pool  (gN/m2)
         call allocate_block_data (grid, f_grainn_storage     ) ! grain nitrogen storage pool  (gN/m2)
         call allocate_block_data (grid, f_grainn_xfer        ) ! grain nitrogen transfer pool (gN/m2)
#endif
         call allocate_block_data (grid, f_retransn           ) ! retranslocation nitrogen pool (gN/m2)

#ifdef CROP
         call allocate_block_data (grid, f_cphase             )  ! crop phase
         call allocate_block_data (grid, f_cropprod1c         )  ! 1-yr crop production carbon
         call allocate_block_data (grid, f_cropprod1c_loss    )  ! loss of 1-yr crop production carbon
         call allocate_block_data (grid, f_cropseedc_deficit  )  ! crop seed carbon deficit
         call allocate_block_data (grid, f_grainc_to_cropprodc ) ! grain to crop production
         call allocate_block_data (grid, f_grainc_to_seed     )  ! grain to crop seed
         call allocate_block_data (grid, f_fert_to_sminn      )  ! grain to crop seed
         call allocate_block_data (grid, f_plantdate          )  ! planting date
#endif
         call allocate_block_data (grid, f_ndep_to_sminn      )  ! grain to crop seed

         call allocate_block_data (grid, f_gpp                ) ! net primary production (gC/m2)
         call allocate_block_data (grid, f_downreg            ) ! gpp downregulation due to N limitation
         call allocate_block_data (grid, f_ar                 )
         call allocate_block_data (grid, f_cwdprod            )
         call allocate_block_data (grid, f_cwddecomp          )
         call allocate_block_data (grid, f_hr                 )
         call allocate_block_data (grid, f_fpg                ) ! fraction of potential gpp
         call allocate_block_data (grid, f_fpi                ) ! fraction of potential immobilization
         call allocate_block_data (grid, f_gpp_enftemp        ) !1
         call allocate_block_data (grid, f_gpp_enfboreal      ) !2
         call allocate_block_data (grid, f_gpp_dnfboreal      ) !3
         call allocate_block_data (grid, f_gpp_ebftrop        ) !4
         call allocate_block_data (grid, f_gpp_ebftemp        ) !5
         call allocate_block_data (grid, f_gpp_dbftrop        ) !6
         call allocate_block_data (grid, f_gpp_dbftemp        ) !7
         call allocate_block_data (grid, f_gpp_dbfboreal      ) !8
         call allocate_block_data (grid, f_gpp_ebstemp        ) !9
         call allocate_block_data (grid, f_gpp_dbstemp        ) !10
         call allocate_block_data (grid, f_gpp_dbsboreal      ) !11
         call allocate_block_data (grid, f_gpp_c3arcgrass     ) !12
         call allocate_block_data (grid, f_gpp_c3grass        ) !13
         call allocate_block_data (grid, f_gpp_c4grass        ) !14
         call allocate_block_data (grid, f_leafc_enftemp      ) !1
         call allocate_block_data (grid, f_leafc_enfboreal    ) !2
         call allocate_block_data (grid, f_leafc_dnfboreal    ) !3
         call allocate_block_data (grid, f_leafc_ebftrop      ) !4
         call allocate_block_data (grid, f_leafc_ebftemp      ) !5
         call allocate_block_data (grid, f_leafc_dbftrop      ) !6
         call allocate_block_data (grid, f_leafc_dbftemp      ) !7
         call allocate_block_data (grid, f_leafc_dbfboreal    ) !8
         call allocate_block_data (grid, f_leafc_ebstemp      ) !9
         call allocate_block_data (grid, f_leafc_dbstemp      ) !10
         call allocate_block_data (grid, f_leafc_dbsboreal    ) !11
         call allocate_block_data (grid, f_leafc_c3arcgrass   ) !12
         call allocate_block_data (grid, f_leafc_c3grass      ) !13
         call allocate_block_data (grid, f_leafc_c4grass      ) !14

         call allocate_block_data (grid, f_litr1c_vr  ,nl_soil)  ! soil carbon pool (gC/m2)
         call allocate_block_data (grid, f_litr2c_vr  ,nl_soil)  ! soil carbon pool (gC/m2)
         call allocate_block_data (grid, f_litr3c_vr  ,nl_soil)  ! soil carbon pool (gC/m2)
         call allocate_block_data (grid, f_cwdc_vr    ,nl_soil)  ! soil carbon pool (gC/m2)
         call allocate_block_data (grid, f_soil1c_vr  ,nl_soil)  ! soil carbon pool (gC/m2)
         call allocate_block_data (grid, f_soil2c_vr  ,nl_soil)  ! soil carbon pool (gC/m2)
         call allocate_block_data (grid, f_soil3c_vr  ,nl_soil)  ! soil carbon pool (gC/m2)

         call allocate_block_data (grid, f_litr1n_vr  ,nl_soil)  ! soil nitrogen pool (gN/m2)
         call allocate_block_data (grid, f_litr2n_vr  ,nl_soil)  ! soil nitrogen pool (gN/m2)
         call allocate_block_data (grid, f_litr3n_vr  ,nl_soil)  ! soil nitrogen pool (gN/m2)
         call allocate_block_data (grid, f_cwdn_vr    ,nl_soil)  ! soil nitrogen pool (gN/m2)
         call allocate_block_data (grid, f_soil1n_vr  ,nl_soil)  ! soil nitrogen pool (gN/m2)
         call allocate_block_data (grid, f_soil2n_vr  ,nl_soil)  ! soil nitrogen pool (gN/m2)
         call allocate_block_data (grid, f_soil3n_vr  ,nl_soil)  ! soil nitrogen pool (gN/m2)
         call allocate_block_data (grid, f_sminn_vr   ,nl_soil)  ! soil mineral nitrogen pool (gN/m2)
#ifdef NITRIF
         call allocate_block_data (grid, f_O2_DECOMP_DEPTH_UNSAT, nl_soil)
         call allocate_block_data (grid, f_CONC_O2_UNSAT        , nl_soil)
#endif
#ifdef CROP
         call allocate_block_data (grid, f_hui                 )
         call allocate_block_data (grid, f_vf                  )
         call allocate_block_data (grid, f_gddmaturity         )
         call allocate_block_data (grid, f_gddplant            )
         call allocate_block_data (grid, f_pdcorn              )
         call allocate_block_data (grid, f_pdswheat            )
         call allocate_block_data (grid, f_pdwwheat            )
         call allocate_block_data (grid, f_pdsoybean           )
         call allocate_block_data (grid, f_pdcotton            )
         call allocate_block_data (grid, f_pdrice1             )
         call allocate_block_data (grid, f_pdrice2             )
         call allocate_block_data (grid, f_pdsugarcane         )
         call allocate_block_data (grid, f_fertnitro_corn      )
         call allocate_block_data (grid, f_fertnitro_swheat    )
         call allocate_block_data (grid, f_fertnitro_wwheat    )
         call allocate_block_data (grid, f_fertnitro_soybean   )
         call allocate_block_data (grid, f_fertnitro_cotton    )
         call allocate_block_data (grid, f_fertnitro_rice1     )
         call allocate_block_data (grid, f_fertnitro_rice2     )
         call allocate_block_data (grid, f_fertnitro_sugarcane )
#endif
#ifdef Fire
         call allocate_block_data (grid, f_abm                 )
         call allocate_block_data (grid, f_gdp                 )
         call allocate_block_data (grid, f_peatf               )
         call allocate_block_data (grid, f_hdm                 )
         call allocate_block_data (grid, f_lnfm                )
#endif

      end if


   END SUBROUTINE allocate_2D_BGCFluxes

#endif
END MODULE MOD_BGC_Vars_2DFluxes
