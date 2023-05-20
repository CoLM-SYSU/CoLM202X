#include <define.h>

MODULE MOD_BGC_Vars_1DFluxes
#ifdef BGC
  !---------------------------------------------------------------------------------------------------------
  ! !DESCRIPTION
  ! Define, allocate, and deallocate biogeochmeical flux variables at patch level

  ! !ORIGINAL:
  ! Xingjie Lu, 2022, created the original version


  USE precision
  IMPLICIT NONE
  SAVE

!--------------------- BGC variables --------------------------------------
! ecosystem vegetation carbon/nitrogen flux
  REAL(r8), allocatable :: gpp                        (:) ! gross primary productivity (gC m-2 s-1)
  REAL(r8), allocatable :: gpp_enftemp                (:) ! gross primary productivity for needleleaf evergreen temperate tree (gC m-2 s-1)
  REAL(r8), allocatable :: gpp_enfboreal              (:) ! gross primary productivity for needleleaf evergreen boreal tree (gC m-2 s-1)
  REAL(r8), allocatable :: gpp_dnfboreal              (:) ! gross primary productivity for needleleaf deciduous boreal tree (gC m-2 s-1)
  REAL(r8), allocatable :: gpp_ebftrop                (:) ! gross primary productivity for broadleaf evergreen tropical tree (gC m-2 s-1)
  REAL(r8), allocatable :: gpp_ebftemp                (:) ! gross primary productivity for broadleaf evergreen temperate tree (gC m-2 s-1)
  REAL(r8), allocatable :: gpp_dbftrop                (:) ! gross primary productivity for broadleaf deciduous tropical tree (gC m-2 s-1) 
  REAL(r8), allocatable :: gpp_dbftemp                (:) ! gross primary productivity for broadleaf deciduous temperate tree (gC m-2 s-1)
  REAL(r8), allocatable :: gpp_dbfboreal              (:) ! gross primary productivity for broadleaf deciduous boreal tree (gC m-2 s-1)
  REAL(r8), allocatable :: gpp_ebstemp                (:) ! gross primary productivity for broadleaf evergreen temperate shrub (gC m-2 s-1)
  REAL(r8), allocatable :: gpp_dbstemp                (:) ! gross primary productivity for broadleaf deciduous temperate shrub (gC m-2 s-1)
  REAL(r8), allocatable :: gpp_dbsboreal              (:) ! gross primary productivity for broadleaf deciduous boreal shrub (gC m-2 s-1)
  REAL(r8), allocatable :: gpp_c3arcgrass             (:) ! gross primary productivity for c3 arctic grass (gC m-2 s-1)
  REAL(r8), allocatable :: gpp_c3grass                (:) ! gross primary productivity for c3 grass (gC m-2 s-1)
  REAL(r8), allocatable :: gpp_c4grass                (:) ! gross primary productivity for c4 grass (gC m-2 s-1)
  REAL(r8), allocatable :: leafc_enftemp              (:) ! leaf carbon display pool for needleleaf evergreen temperate tree (gC m-2)
  REAL(r8), allocatable :: leafc_enfboreal            (:) ! leaf carbon display pool for needleleaf evergreen boreal tree (gC m-2)
  REAL(r8), allocatable :: leafc_dnfboreal            (:) ! leaf carbon display pool for needleleaf deciduous boreal tree (gC m-2)
  REAL(r8), allocatable :: leafc_ebftrop              (:) ! leaf carbon display pool for broadleaf evergreen tropical tree (gC m-2)
  REAL(r8), allocatable :: leafc_ebftemp              (:) ! leaf carbon display pool for broadleaf evergreen temperate tree (gC m-2)
  REAL(r8), allocatable :: leafc_dbftrop              (:) ! leaf carbon display pool for broadleaf deciduous tropical tree (gC m-2)
  REAL(r8), allocatable :: leafc_dbftemp              (:) ! leaf carbon display pool for broadleaf deciduous temperate tree (gC m-2)
  REAL(r8), allocatable :: leafc_dbfboreal            (:) ! leaf carbon display pool for broadleaf deciduous boreal tree (gC m-2)
  REAL(r8), allocatable :: leafc_ebstemp              (:) ! leaf carbon display pool for broadleaf evergreen temperate shrub (gC m-2)
  REAL(r8), allocatable :: leafc_dbstemp              (:) ! leaf carbon display pool for broadleaf deciduous temperate shrub (gC m-2)
  REAL(r8), allocatable :: leafc_dbsboreal            (:) ! leaf carbon display pool for broadleaf deciduous boreal shrub (gC m-2)
  REAL(r8), allocatable :: leafc_c3arcgrass           (:) ! leaf carbon display pool for c3 arctic grass (gC m-2)
  REAL(r8), allocatable :: leafc_c3grass              (:) ! leaf carbon display pool for c3 grass (gC m-2)
  REAL(r8), allocatable :: leafc_c4grass              (:) ! leaf carbon display pool for c4 grass (gC m-2)
  REAL(r8), allocatable :: ar                         (:) ! autotrophic respiration (gC m-2 s-1)
  REAL(r8), allocatable :: cwdprod                    (:) ! CWD production (gC m-2 s-1)
  REAL(r8), allocatable :: cwddecomp                  (:) ! CWD decomposition (gC m-2 s-1)
  REAL(r8), allocatable :: hr                         (:) ! heterotrophic respiration (gC m-2 s-1)
  REAL(r8), allocatable :: er                         (:) ! total ecosystem respiration, autotrophic + heterotrophic (gC m-2 s-1)
  REAL(r8), allocatable :: fire_closs                 (:) ! total C emissions due to fire (gC m-2 s-1)
  REAL(r8), allocatable :: fire_nloss                 (:) ! total N emissions due to fire (gN m-2 s-1)
  REAL(r8), allocatable :: hrv_xsmrpool_to_atm        (:) ! maintenance respiration storage C to atmosphere due to harvest (gC m-2 s-1)
  REAL(r8), allocatable :: wood_harvestc              (:) ! harvested wood C (gC m-2 s-1)
  REAL(r8), allocatable :: wood_harvestn              (:) ! harvested wood N (gN m-2 s-1)
  REAL(r8), allocatable :: grainc_to_cropprodc        (:) ! grain to crop production carbon (gC m-2 s-1) 
  REAL(r8), allocatable :: grainc_to_seed             (:) ! grain to crop seed carbon (gC m-2 s-1) 
  REAL(r8), allocatable :: grainn_to_cropprodn        (:) ! grain to crop production nitrogen (gN m-2 s-1)
  REAL(r8), allocatable :: cropprod1c_loss            (:) ! loss rate of 1-yr crop production carbon (gC m-2 s-1)

! decomposition carbon fluxes
  REAL(r8), allocatable :: decomp_cpools_sourcesink   (:,:,:)       ! vertical resolved: the input of litter & soil carbon pools (donor or receiver) from phenology-associated litterfall and decomposition (gC m-3 timestep-1)
  REAL(r8), allocatable :: decomp_ctransfer_vr        (:,:,:)       ! vertical resolved: the non-respiratory portion of potential carbon transfer from one litter & soil carbon pool to another (gC m-3 s-1)
  REAL(r8), allocatable :: decomp_hr_vr               (:,:,:)       ! vertical resolved: the heterotrophic respiration portion of potential carbon loss from one litter & soil carbon pool to another (gC m-3 s-1)
  REAL(r8), allocatable :: decomp_hr                  (:)           ! the heterotrophic respiration portion of potential carbon loss from one decomposition carbon pool to another (gC m-3 s-1)
  REAL(r8), allocatable :: phr_vr                     (:,:)         ! vertical resolved: the potential heterotrophic respiration carbon flux (gC m-3 s-1)
  REAL(r8), allocatable :: m_decomp_cpools_to_fire_vr (:,:,:)       ! vertical resolved: the carbon from decomposition pools to fire emissions (gC m-3 s-1)
  REAL(r8), allocatable :: decomp_cpools_transport_tendency(:,:,:)  ! vertical resolved: the carbon tendency due to vertical transport in decomposition carbon pools (gC m-3 s-1)
  REAL(r8), allocatable :: som_c_leached              (:)           ! total soil organic matter C loss from vertical transport (gC m-2 s-1)

! vegetation to decomposition carbon fluxes
  REAL(r8), allocatable :: phenology_to_met_c       (:,:) ! phenology-associated plant C loss to metabolic litter C (gC m-3 s-1)
  REAL(r8), allocatable :: phenology_to_cel_c       (:,:) ! phenology-associated plant C loss to cellulosic litter C (gC m-3 s-1)
  REAL(r8), allocatable :: phenology_to_lig_c       (:,:) ! phenology-associated plant C loss to lignin litter C (gC m-3 s-1)
  REAL(r8), allocatable :: gap_mortality_to_met_c   (:,:) ! gap mortality-associated plant C loss to metabolic litter C (gC m-3 s-1)
  REAL(r8), allocatable :: gap_mortality_to_cel_c   (:,:) ! gap mortality-associated plant C loss to cellulosic litter C (gC m-3 s-1)
  REAL(r8), allocatable :: gap_mortality_to_lig_c   (:,:) ! gap mortality-associated plant C loss to lignin litter C (gC m-3 s-1)
  REAL(r8), allocatable :: gap_mortality_to_cwdc    (:,:) ! gap mortality-associated plant C loss to coarse woody debris C (gC m-3 s-1)
  REAL(r8), allocatable :: fire_mortality_to_met_c  (:,:) ! fire mortality-associated plant C loss to metabolic litter C (gC m-3 s-1)
  REAL(r8), allocatable :: fire_mortality_to_cel_c  (:,:) ! fire mortality-associated plant C loss to cellulosic litter C (gC m-3 s-1)
  REAL(r8), allocatable :: fire_mortality_to_lig_c  (:,:) ! fire mortality-associated plant C loss to lignin litter C (gC m-3 s-1)
  REAL(r8), allocatable :: fire_mortality_to_cwdc   (:,:) ! fire mortality-associated plant C loss to coarse woody debris C (gC m-3 s-1)

! decomposition nitrogen fluxes
  REAL(r8), allocatable :: decomp_npools_sourcesink   (:,:,:)       ! vertical resolved: the the input of litter & soil nitrogen pools (donor or receiver) (gN m-3 timestep)
  REAL(r8), allocatable :: decomp_ntransfer_vr        (:,:,:)       ! vertical resolved: the nitrogen flux transfer from one litter & soil nitrogen pool to another (gN m-3 s-1)
  REAL(r8), allocatable :: decomp_sminn_flux_vr       (:,:,:)       ! vertical resolved: the nitrogen mineralization flux from each nitrogen transfer between litter & soil pools (gN m-3 s-1)
  REAL(r8), allocatable :: sminn_to_denit_decomp_vr   (:,:,:)       ! vertical resolved: the nitrogen denitrification flux from each nitrogen transfer between litter & soil pools (gN m-3 s-1) 
  REAL(r8), allocatable :: m_decomp_npools_to_fire_vr (:,:,:)       ! vertical resolved: the litter & soil nitrogen loss associated to the fire (gN m-3 s-1)
  REAL(r8), allocatable :: decomp_npools_transport_tendency(:,:,:)  ! vertical resolved: the nitrogen tendency due to vertical transport in decomposition nitrogen pools (gN m-3 s-1)
  REAL(r8), allocatable :: som_n_leached              (:)           ! total soil organic matter N loss from vertical transport (gN m-2 s-1)

! vegetation to decomposition nitrogen fluxes
  REAL(r8), allocatable :: phenology_to_met_n       (:,:)   ! phenology-associated plant N loss to metabolic litter N (gN m-3 s-1)
  REAL(r8), allocatable :: phenology_to_cel_n       (:,:)   ! phenology-associated plant N loss to cellulosic litter N (gN m-3 s-1)
  REAL(r8), allocatable :: phenology_to_lig_n       (:,:)   ! phenology-associated plant N loss to lignin litter N (gN m-3 s-1)
  REAL(r8), allocatable :: gap_mortality_to_met_n   (:,:)   ! gap mortality-associated plant N loss to metabolic litter N (gN m-3 s-1)
  REAL(r8), allocatable :: gap_mortality_to_cel_n   (:,:)   ! gap mortality-associated plant N loss to cellulosic litter N (gN m-3 s-1)
  REAL(r8), allocatable :: gap_mortality_to_lig_n   (:,:)   ! gap mortality-associated plant N loss to lignin litter N (gN m-3 s-1)
  REAL(r8), allocatable :: gap_mortality_to_cwdn    (:,:)   ! gap mortality-associated plant N loss to coarse woody debris N (gN m-3 s-1)
  REAL(r8), allocatable :: fire_mortality_to_met_n  (:,:)   ! fire mortality-associated plant N loss to metabolic litter N (gN m-3 s-1)
  REAL(r8), allocatable :: fire_mortality_to_cel_n  (:,:)   ! fire mortality-associated plant N loss to cellulosic litter N (gN m-3 s-1)
  REAL(r8), allocatable :: fire_mortality_to_lig_n  (:,:)   ! fire mortality-associated plant N loss to lignin litter N (gN m-3 s-1)
  REAL(r8), allocatable :: fire_mortality_to_cwdn   (:,:)   ! fire mortality-associated plant N loss to coarse woody debris N (gN m-3 s-1)

  REAL(r8), allocatable :: sminn_leached_vr         (:,:)   ! vertical resolved: soil mineral N loss due to leaching (gN m-3 s-1)
  REAL(r8), allocatable :: smin_no3_leached_vr      (:,:)   ! vertical resolved: soil mineral NO3 loss due to leaching (gN m-3 s-1)
  REAL(r8), allocatable :: smin_no3_runoff_vr       (:,:)   ! vertical resolved: soil mineral NO3 loss due to runoff (gN m-3 s-1)
  REAL(r8), allocatable :: net_nmin_vr              (:,:)   ! vertical resolved: net N mineralization (gN m-3 s-1) 
  REAL(r8), allocatable :: gross_nmin_vr            (:,:)   ! vertical resolved: total N mineralization (gN m-3 s-1) 
  REAL(r8), allocatable :: net_nmin                 (:)     ! net N mineralization (gN m-2 s-1) 
  REAL(r8), allocatable :: gross_nmin               (:)     ! total N mineralization (gN m-2 s-1) 
  REAL(r8), allocatable :: plant_ndemand            (:)     ! potential plant N uptake (gN m-2 s-1) 
  REAL(r8), allocatable :: actual_immob_vr          (:,:)   ! vertical resolved: actual N immobilization (gN m-3 s-1)
  REAL(r8), allocatable :: actual_immob_nh4_vr      (:,:)   ! vertical resolved: actual NH4 immobilization (gN m-3 s-1)
  REAL(r8), allocatable :: actual_immob_no3_vr      (:,:)   ! vertical resolved: actual NO3 immobilization (gN m-3 s-1)
  REAL(r8), allocatable :: potential_immob_vr       (:,:)   ! vertical resolved: potential N immobilization (gN m-3 s-1)
  REAL(r8), allocatable :: pmnf_decomp              (:,:,:) ! vertical resolved: potential N mineralization flux of each transfer between litter & soil pools (gN m-3 s-1)
  REAL(r8), allocatable :: p_decomp_cpool_loss      (:,:,:) ! vertical resolved: potential C exit rate (transfer+hr) of the donor pool of each transfer between litter & soil pools (gC m-3 s-1)
  REAL(r8), allocatable :: sminn_to_plant           (:)     ! plant uptake N (gN m-2 s-1)
  REAL(r8), allocatable :: sminn_to_plant_vr        (:,:)   ! vertical resolved: plant uptake N (gN m-3 s-1)
  REAL(r8), allocatable :: smin_nh4_to_plant_vr     (:,:)   ! vertical resolved: plant uptake NH4 (gN m-3 s-1)
  REAL(r8), allocatable :: smin_no3_to_plant_vr     (:,:)   ! vertical resolved: plant uptake NO3 (gN m-3 s-1)
  REAL(r8), allocatable :: supplement_to_sminn_vr   (:,:)   ! vertical resolved: supplemental N supply to soil mineral N (gN m-3 s-1)
  REAL(r8), allocatable :: sminn_to_plant_fun_vr    (:,:)   ! vertical resolved: plant uptake N by FUN (gN m-3 s-1)
  REAL(r8), allocatable :: sminn_to_plant_fun_nh4_vr(:,:)   ! vertical resolved: plant uptake NH4 by FUN (gN m-3 s-1)
  REAL(r8), allocatable :: sminn_to_plant_fun_no3_vr(:,:)   ! vertical resolved: plant uptake NO3 by FUN (gN m-3 s-1)
  REAL(r8), allocatable :: sminn_to_denit_excess_vr (:,:)   ! vertical resolved: denitrification from excess mineral N (gN m-3 s-1)
  REAL(r8), allocatable :: f_nit_vr                 (:,:)   ! vertical resolved: nitrification (gN m-3 s-1)
  REAL(r8), allocatable :: f_denit_vr               (:,:)   ! vertical resolved: denitrification (gN m-3 s-1)
  REAL(r8), allocatable :: f_n2o_nit_vr             (:,:)   ! vertical resolved: N2O emission from N nitrification (gN m-3 s-1)
  REAL(r8), allocatable :: f_n2o_denit_vr           (:,:)   ! vertical resolved: N2O emission from N denitrification (gN m-3 s-1)
  REAL(r8), allocatable :: pot_f_nit_vr             (:,:)   ! vertical resolved: potential N nitrification (gN m-3 s-1)
  REAL(r8), allocatable :: pot_f_denit_vr           (:,:)   ! vertical resolved: potential N denitrification (gN m-3 s-1)
  REAL(r8), allocatable :: n2_n2o_ratio_denit_vr    (:,:)   ! vertical resolved: ratio of N2 to N2O production by denitrification (gN gN-1)
  REAL(r8), allocatable :: ndep_to_sminn            (:)     ! atmospheric N deposition to soil mineral N (gN m-2 s-1)
  REAL(r8), allocatable :: ffix_to_sminn            (:)     ! free living N fixation to soil mineral N (gN m-2 s-1)
  REAL(r8), allocatable :: nfix_to_sminn            (:)     ! N fixation to soil mineral N (gN m-2 s-1)
  REAL(r8), allocatable :: somc_fire                (:)     ! soil organic matters C to fire emissions (gC m-2 s-1)
  REAL(r8), allocatable :: supplement_to_sminn      (:)     ! supplemental N supply to soil mineral N (gN m-2 s-1)
  REAL(r8), allocatable :: fert_to_sminn            (:)     ! fertilizer N to soil mineral N (gN m-2 s-1)
  REAL(r8), allocatable :: soyfixn_to_sminn         (:)     ! soybean N fixation to soil mineral N (gN m-2 s-1)
  REAL(r8), allocatable :: denit                    (:)     ! total N denitrification (gN m-2 s-1)
  REAL(r8), allocatable :: sminn_leached            (:)     ! soil mineral N loss due to leaching (gN m-2 s-1)
  REAL(r8), allocatable :: f_n2o_nit                (:)     ! flux of N2O from N nitrification (gN m-2 s-1)
  REAL(r8), allocatable :: smin_no3_leached         (:)     ! soil mineral NO3 loss due to leaching (gN m-2 s-1)
  REAL(r8), allocatable :: smin_no3_runoff          (:)     ! soil mineral NO3 loss due to runoff (gN m-2 s-1)
!----------------- end BGC variables -----------------------------------

! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: allocate_1D_BGCFluxes
  PUBLIC :: deallocate_1D_BGCFluxes

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_1D_BGCFluxes
  ! --------------------------------------------------------------------
  ! Allocates memory for CLM 1d [numpatch] variables 
  ! --------------------------------------------------------------------
     USE precision
     USE GlobalVars
     USE spmd_task
     USE mod_landpatch
     IMPLICIT NONE


      if (p_is_worker) then

         if (numpatch > 0) then

! bgc variables
! ecosystem carbon flux
            allocate (gpp                        (numpatch))
            allocate (gpp_enftemp                (numpatch)) !1
            allocate (gpp_enfboreal              (numpatch)) !2
            allocate (gpp_dnfboreal              (numpatch)) !3
            allocate (gpp_ebftrop                (numpatch)) !4
            allocate (gpp_ebftemp                (numpatch)) !5
            allocate (gpp_dbftrop                (numpatch)) !6
            allocate (gpp_dbftemp                (numpatch)) !7
            allocate (gpp_dbfboreal              (numpatch)) !8
            allocate (gpp_ebstemp                (numpatch)) !9
            allocate (gpp_dbstemp                (numpatch)) !10
            allocate (gpp_dbsboreal              (numpatch)) !11
            allocate (gpp_c3arcgrass             (numpatch)) !12
            allocate (gpp_c3grass                (numpatch)) !13
            allocate (gpp_c4grass                (numpatch)) !14
            allocate (leafc_enftemp              (numpatch)) !1
            allocate (leafc_enfboreal            (numpatch)) !2
            allocate (leafc_dnfboreal            (numpatch)) !3
            allocate (leafc_ebftrop              (numpatch)) !4
            allocate (leafc_ebftemp              (numpatch)) !5
            allocate (leafc_dbftrop              (numpatch)) !6
            allocate (leafc_dbftemp              (numpatch)) !7
            allocate (leafc_dbfboreal            (numpatch)) !8
            allocate (leafc_ebstemp              (numpatch)) !9
            allocate (leafc_dbstemp              (numpatch)) !10
            allocate (leafc_dbsboreal            (numpatch)) !11
            allocate (leafc_c3arcgrass           (numpatch)) !12
            allocate (leafc_c3grass              (numpatch)) !13
            allocate (leafc_c4grass              (numpatch)) !14
            allocate (ar                         (numpatch))
            allocate (cwdprod                    (numpatch))
            allocate (cwddecomp                  (numpatch))
            allocate (hr                         (numpatch))
            allocate (er                         (numpatch))
            allocate (fire_closs                 (numpatch))
            allocate (fire_nloss                 (numpatch))
            allocate (hrv_xsmrpool_to_atm        (numpatch))
            allocate (wood_harvestc              (numpatch))
            allocate (wood_harvestn              (numpatch))
            allocate (grainc_to_cropprodc        (numpatch))
            allocate (grainc_to_seed             (numpatch))
            allocate (grainn_to_cropprodn        (numpatch))
            allocate (cropprod1c_loss            (numpatch))


! decomposition carbon fluxes
            allocate (decomp_cpools_sourcesink   (nl_soil_full,ndecomp_pools,numpatch))
            allocate (decomp_ctransfer_vr        (nl_soil_full,ndecomp_transitions,numpatch))
            allocate (decomp_hr_vr               (nl_soil_full,ndecomp_transitions,numpatch))
            allocate (decomp_hr                  (numpatch))
            allocate (phr_vr                     (nl_soil_full,numpatch))
            allocate (m_decomp_cpools_to_fire_vr (nl_soil_full,ndecomp_pools,numpatch))
            allocate (decomp_cpools_transport_tendency(nl_soil_full,ndecomp_pools,numpatch))
            allocate (som_c_leached              (numpatch))

! vegetation to decomposition carbon fluxes
            allocate (phenology_to_met_c       (nl_soil,numpatch))
            allocate (phenology_to_cel_c       (nl_soil,numpatch))
            allocate (phenology_to_lig_c       (nl_soil,numpatch))
            allocate (gap_mortality_to_met_c   (nl_soil,numpatch))
            allocate (gap_mortality_to_cel_c   (nl_soil,numpatch))
            allocate (gap_mortality_to_lig_c   (nl_soil,numpatch))
            allocate (gap_mortality_to_cwdc    (nl_soil,numpatch))
            allocate (fire_mortality_to_met_c  (nl_soil,numpatch))
            allocate (fire_mortality_to_cel_c  (nl_soil,numpatch))
            allocate (fire_mortality_to_lig_c  (nl_soil,numpatch))
            allocate (fire_mortality_to_cwdc   (nl_soil,numpatch))  

! decomposition nitrogen fluxes
            allocate (decomp_npools_sourcesink   (nl_soil_full,ndecomp_pools,numpatch))
            allocate (decomp_ntransfer_vr        (nl_soil_full,ndecomp_transitions,numpatch))
            allocate (decomp_sminn_flux_vr       (nl_soil_full,ndecomp_transitions,numpatch))
            allocate (sminn_to_denit_decomp_vr   (nl_soil_full,ndecomp_transitions,numpatch))
            allocate (m_decomp_npools_to_fire_vr (nl_soil_full,ndecomp_pools,numpatch))
            allocate (decomp_npools_transport_tendency(nl_soil_full,ndecomp_pools,numpatch))
            allocate (som_n_leached            (numpatch))

! vegetation to decomposition nitrogen fluxes
            allocate (phenology_to_met_n       (nl_soil,numpatch))
            allocate (phenology_to_cel_n       (nl_soil,numpatch))
            allocate (phenology_to_lig_n       (nl_soil,numpatch))
            allocate (gap_mortality_to_met_n   (nl_soil,numpatch))
            allocate (gap_mortality_to_cel_n   (nl_soil,numpatch))
            allocate (gap_mortality_to_lig_n   (nl_soil,numpatch))
            allocate (gap_mortality_to_cwdn    (nl_soil,numpatch))
            allocate (fire_mortality_to_met_n  (nl_soil,numpatch))
            allocate (fire_mortality_to_cel_n  (nl_soil,numpatch))
            allocate (fire_mortality_to_lig_n  (nl_soil,numpatch))
            allocate (fire_mortality_to_cwdn   (nl_soil,numpatch))

            allocate (sminn_leached_vr         (nl_soil,numpatch))
            allocate (smin_no3_leached_vr      (nl_soil,numpatch))
            allocate (smin_no3_runoff_vr       (nl_soil,numpatch))
            allocate (net_nmin_vr              (nl_soil,numpatch))
            allocate (gross_nmin_vr            (nl_soil,numpatch))
            allocate (net_nmin                 (numpatch))
            allocate (gross_nmin               (numpatch))
            allocate (plant_ndemand            (numpatch))
            allocate (actual_immob_vr          (nl_soil,numpatch))
            allocate (actual_immob_nh4_vr      (nl_soil,numpatch))
            allocate (actual_immob_no3_vr      (nl_soil,numpatch))
            allocate (potential_immob_vr       (nl_soil,numpatch))
            allocate (pmnf_decomp              (nl_soil,ndecomp_transitions,numpatch))
            allocate (p_decomp_cpool_loss      (nl_soil,ndecomp_transitions,numpatch))
            allocate (sminn_to_plant           (numpatch))
            allocate (sminn_to_plant_vr        (nl_soil,numpatch))
            allocate (smin_nh4_to_plant_vr     (nl_soil,numpatch))
            allocate (smin_no3_to_plant_vr     (nl_soil,numpatch))
            allocate (supplement_to_sminn_vr   (nl_soil,numpatch))
            allocate (sminn_to_plant_fun_vr    (nl_soil,numpatch))
            allocate (sminn_to_plant_fun_nh4_vr(nl_soil,numpatch))
            allocate (sminn_to_plant_fun_no3_vr(nl_soil,numpatch))
            allocate (sminn_to_denit_excess_vr (nl_soil,numpatch))
            allocate (f_nit_vr                 (nl_soil,numpatch))
            allocate (f_denit_vr               (nl_soil,numpatch))
            allocate (f_n2o_nit_vr             (nl_soil,numpatch))
            allocate (f_n2o_denit_vr           (nl_soil,numpatch))
            allocate (pot_f_nit_vr             (nl_soil,numpatch))
            allocate (pot_f_denit_vr           (nl_soil,numpatch))
            allocate (n2_n2o_ratio_denit_vr    (nl_soil,numpatch))
            allocate (ndep_to_sminn            (numpatch))
            allocate (ffix_to_sminn            (numpatch))
            allocate (nfix_to_sminn            (numpatch))
            allocate (somc_fire                (numpatch))
            allocate (supplement_to_sminn      (numpatch))
            allocate (fert_to_sminn            (numpatch))
            allocate (soyfixn_to_sminn         (numpatch))
            allocate (denit                    (numpatch))
            allocate (sminn_leached            (numpatch))
            allocate (f_n2o_nit                (numpatch))
            allocate (smin_no3_leached         (numpatch))
            allocate (smin_no3_runoff          (numpatch))
         end if
      end if


   END SUBROUTINE allocate_1D_BGCFluxes

   SUBROUTINE deallocate_1D_BGCFluxes ()
  ! --------------------------------------------------------------------
  ! deallocates memory for CLM 1d [numpatch] variables
  ! --------------------------------------------------------------------
     USE spmd_task
     USE mod_landpatch

     if (p_is_worker) then
        
        if (numpatch > 0) then

! bgc variables
! ecosystem carbon flux
           deallocate (gpp                        )
           deallocate (gpp_enftemp                ) !1
           deallocate (gpp_enfboreal              ) !2
           deallocate (gpp_dnfboreal              ) !3
           deallocate (gpp_ebftrop                ) !4
           deallocate (gpp_ebftemp                ) !5
           deallocate (gpp_dbftrop                ) !6
           deallocate (gpp_dbftemp                ) !7
           deallocate (gpp_dbfboreal              ) !8
           deallocate (gpp_ebstemp                ) !9
           deallocate (gpp_dbstemp                ) !10
           deallocate (gpp_dbsboreal              ) !11
           deallocate (gpp_c3arcgrass             ) !12
           deallocate (gpp_c3grass                ) !13
           deallocate (gpp_c4grass                ) !14
           deallocate (leafc_enftemp              ) !1
           deallocate (leafc_enfboreal            ) !2
           deallocate (leafc_dnfboreal            ) !3
           deallocate (leafc_ebftrop              ) !4
           deallocate (leafc_ebftemp              ) !5
           deallocate (leafc_dbftrop              ) !6
           deallocate (leafc_dbftemp              ) !7
           deallocate (leafc_dbfboreal            ) !8
           deallocate (leafc_ebstemp              ) !9
           deallocate (leafc_dbstemp              ) !10
           deallocate (leafc_dbsboreal            ) !11
           deallocate (leafc_c3arcgrass           ) !12
           deallocate (leafc_c3grass              ) !13
           deallocate (leafc_c4grass              ) !14
           deallocate (ar                         )
           deallocate (cwdprod                    )
           deallocate (cwddecomp                  )
           deallocate (hr                         )
           deallocate (er                         )
           deallocate (fire_closs                 )
           deallocate (fire_nloss                 )
           deallocate (hrv_xsmrpool_to_atm        )
           deallocate (wood_harvestc              )
           deallocate (wood_harvestn              )
           deallocate (grainc_to_cropprodc        )
           deallocate (grainc_to_seed             )
           deallocate (grainn_to_cropprodn        )
           deallocate (cropprod1c_loss            )
      
      
! decomposition carbon fluxes
           deallocate (decomp_cpools_sourcesink   )
           deallocate (decomp_ctransfer_vr        )
           deallocate (decomp_hr_vr               )
           deallocate (decomp_hr                  )
           deallocate (phr_vr                     )
           deallocate (m_decomp_cpools_to_fire_vr )
           deallocate (decomp_cpools_transport_tendency)
           deallocate (som_c_leached              )
      
! vegetation to decomposition carbon fluxes
           deallocate (phenology_to_met_c       )
           deallocate (phenology_to_cel_c       )
           deallocate (phenology_to_lig_c       )
           deallocate (gap_mortality_to_met_c   )
           deallocate (gap_mortality_to_cel_c   )
           deallocate (gap_mortality_to_lig_c   )
           deallocate (gap_mortality_to_cwdc    )
           deallocate (fire_mortality_to_met_c  )
           deallocate (fire_mortality_to_cel_c  )
           deallocate (fire_mortality_to_lig_c  )
           deallocate (fire_mortality_to_cwdc   )

! decomposition nitrogen fluxes
           deallocate (decomp_npools_sourcesink   )
           deallocate (decomp_ntransfer_vr        )
           deallocate (decomp_sminn_flux_vr       )
           deallocate (sminn_to_denit_decomp_vr   )
           deallocate (m_decomp_npools_to_fire_vr )
           deallocate (decomp_npools_transport_tendency)
           deallocate (som_n_leached              )

! vegetation to decomposition nitrogen fluxes
           deallocate (phenology_to_met_n       )
           deallocate (phenology_to_cel_n       )
           deallocate (phenology_to_lig_n       )
           deallocate (gap_mortality_to_met_n   )
           deallocate (gap_mortality_to_cel_n   )
           deallocate (gap_mortality_to_lig_n   )
           deallocate (gap_mortality_to_cwdn    )
           deallocate (fire_mortality_to_met_n  )
           deallocate (fire_mortality_to_cel_n  )
           deallocate (fire_mortality_to_lig_n  )
           deallocate (fire_mortality_to_cwdn   )

           deallocate (sminn_leached_vr         )
           deallocate (smin_no3_leached_vr      )
           deallocate (smin_no3_runoff_vr       )
           deallocate (net_nmin_vr              )
           deallocate (gross_nmin_vr            )
           deallocate (net_nmin                 )
           deallocate (gross_nmin               )
           deallocate (plant_ndemand            )
           deallocate (actual_immob_vr          )
           deallocate (actual_immob_nh4_vr      )
           deallocate (actual_immob_no3_vr      )
           deallocate (potential_immob_vr       )
           deallocate (pmnf_decomp              )
           deallocate (p_decomp_cpool_loss      )
           deallocate (sminn_to_plant           )
           deallocate (sminn_to_plant_vr        )
           deallocate (smin_nh4_to_plant_vr     )
           deallocate (smin_no3_to_plant_vr     )
           deallocate (supplement_to_sminn_vr   )
           deallocate (sminn_to_plant_fun_vr    )
           deallocate (sminn_to_plant_fun_nh4_vr)
           deallocate (sminn_to_plant_fun_no3_vr)
           deallocate (sminn_to_denit_excess_vr )
           deallocate (f_nit_vr                 )
           deallocate (f_denit_vr               )
           deallocate (f_n2o_nit_vr             )
           deallocate (f_n2o_denit_vr           )
           deallocate (pot_f_nit_vr             )
           deallocate (n2_n2o_ratio_denit_vr    )
           deallocate (ndep_to_sminn            )
           deallocate (ffix_to_sminn            )
           deallocate (nfix_to_sminn            )
           deallocate (somc_fire                )
           deallocate (supplement_to_sminn      )
           deallocate (fert_to_sminn            )
           deallocate (soyfixn_to_sminn         )
           deallocate (denit                    )
           deallocate (sminn_leached            )
           deallocate (f_n2o_nit                )
           deallocate (smin_no3_leached         )
           deallocate (smin_no3_runoff          )
        end if
     end if

   END SUBROUTINE deallocate_1D_BGCFluxes

#endif
END MODULE MOD_BGC_Vars_1DFluxes
! ---------- EOP ------------

