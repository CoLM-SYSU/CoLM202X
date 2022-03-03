#include <define.h> 

MODULE MOD_TimeVariables
! -------------------------------
! Created by Yongjiu Dai, 03/2014
! -------------------------------

  USE precision
  USE timemanager
  IMPLICIT NONE
  SAVE

! -----------------------------------------------------------------
! Time-varying state variables which reaquired by restart run
  REAL(r8), allocatable :: z_sno       (:,:) !node depth [m]
  REAL(r8), allocatable :: dz_sno      (:,:) !interface depth [m]
  REAL(r8), allocatable :: t_soisno    (:,:) !soil temperature [K]
  REAL(r8), allocatable :: wliq_soisno (:,:) !liquid water in layers [kg/m2]
  REAL(r8), allocatable :: wice_soisno (:,:) !ice lens in layers [kg/m2]
  REAL(r8), allocatable :: smp         (:,:) !soil matric potential [mm]
  REAL(r8), allocatable :: h2osoi      (:,:) !volumetric soil water in layers [m3/m3]
  REAL(r8), allocatable :: rstfac        (:) !factor of soil water stress 
  REAL(r8), allocatable :: t_grnd        (:) !ground surface temperature [K]

  REAL(r8), allocatable :: tleaf         (:) !leaf temperature [K]
  REAL(r8), allocatable :: ldew          (:) !depth of water on foliage [mm]
  REAL(r8), allocatable :: sag           (:) !non dimensional snow age [-]
  REAL(r8), allocatable :: scv           (:) !snow cover, water equivalent [mm]
  REAL(r8), allocatable :: snowdp        (:) !snow depth [meter]
  REAL(r8), allocatable :: fveg          (:) !fraction of vegetation cover
  REAL(r8), allocatable :: fsno          (:) !fraction of snow cover on ground
  REAL(r8), allocatable :: sigf          (:) !fraction of veg cover, excluding snow-covered veg [-]
  REAL(r8), allocatable :: green         (:) !leaf greenness
  REAL(r8), allocatable :: tlai          (:) !leaf area index
  REAL(r8), allocatable :: lai           (:) !leaf area index
  REAL(r8), allocatable :: laisun        (:) !leaf area index
  REAL(r8), allocatable :: laisha        (:) !leaf area index
  REAL(r8), allocatable :: tsai          (:) !stem area index
  REAL(r8), allocatable :: sai           (:) !stem area index
  REAL(r8), allocatable :: coszen        (:) !cosine of solar zenith angle
  REAL(r8), allocatable :: alb       (:,:,:) !averaged albedo [-]
  REAL(r8), allocatable :: ssun      (:,:,:) !sunlit canopy absorption for solar radiation (0-1)
  REAL(r8), allocatable :: ssha      (:,:,:) !shaded canopy absorption for solar radiation (0-1)
  REAL(r8), allocatable :: thermk        (:) !canopy gap fraction for tir radiation
  REAL(r8), allocatable :: extkb         (:) !(k, g(mu)/mu) direct solar extinction coefficient
  REAL(r8), allocatable :: extkd         (:) !diffuse and scattered diffuse PAR extinction coefficient
  REAL(r8), allocatable :: zwt           (:) !the depth to water table [m]
  REAL(r8), allocatable :: wa            (:) !water storage in aquifer [mm]
  REAL(r8), allocatable :: wat           (:) !total water storage [mm]

  REAL(r8), allocatable :: t_lake      (:,:) !lake layer teperature [K]
  REAL(r8), allocatable :: lake_icefrac(:,:) !lake mass fraction of lake layer that is frozen

  REAL(r8), allocatable :: trad          (:) !radiative temperature of surface [K]
  REAL(r8), allocatable :: tref          (:) !2 m height air temperature [kelvin]
  REAL(r8), allocatable :: qref          (:) !2 m height air specific humidity
  REAL(r8), allocatable :: rst           (:) !canopy stomatal resistance (s/m)
  REAL(r8), allocatable :: emis          (:) !averaged bulk surface emissivity
  REAL(r8), allocatable :: z0m           (:) !effective roughness [m]
  REAL(r8), allocatable :: displa        (:) !zero displacement height [m]
  REAL(r8), allocatable :: zol           (:) !dimensionless height (z/L) used in Monin-Obukhov theory
  REAL(r8), allocatable :: rib           (:) !bulk Richardson number in surface layer
  REAL(r8), allocatable :: ustar         (:) !u* in similarity theory [m/s]
  REAL(r8), allocatable :: qstar         (:) !q* in similarity theory [kg/kg]
  REAL(r8), allocatable :: tstar         (:) !t* in similarity theory [K]
  REAL(r8), allocatable :: fm            (:) !integral of profile function for momentum
  REAL(r8), allocatable :: fh            (:) !integral of profile function for heat
  REAL(r8), allocatable :: fq            (:) !integral of profile function for moisture

! bgc variables
  REAL(r8), allocatable :: decomp_cpools_vr  (:,:,:)
  REAL(r8), allocatable :: decomp_cpools     (:,:)
  REAL(r8), allocatable :: decomp_k          (:,:,:) ! soil decomposition rate [1/s]
  REAL(r8), allocatable :: ctrunc_vr         (:,:)
  REAL(r8), allocatable :: ctrunc_veg        (:)
  REAL(r8), allocatable :: ctrunc_soil       (:)

  REAL(r8), allocatable :: t_scalar          (:,:)   ! soil decomposition temperature scalars [unitless]
  REAL(r8), allocatable :: w_scalar          (:,:)   ! soil decomposition water scalars [unitless]
  REAL(r8), allocatable :: o_scalar          (:,:)   ! soil decomposition oxygen scalars [unitless]
  REAL(r8), allocatable :: depth_scalar      (:,:)   ! soil decomposition depth scalars [unitless]

! Soil CN diffusion and advection
  REAL(r8), allocatable :: som_adv_coef             (:,:)
  REAL(r8), allocatable :: som_diffus_coef          (:,:)

! Active Layer
  REAL(r8), allocatable :: altmax                   (:)
  REAL(r8), allocatable :: altmax_lastyear          (:)
  INTEGER , allocatable :: altmax_lastyear_indx     (:)

  REAL(r8), allocatable :: totlitc                  (:)
  REAL(r8), allocatable :: totvegc                  (:)
  REAL(r8), allocatable :: totsomc                  (:)
  REAL(r8), allocatable :: totcwdc                  (:)
  REAL(r8), allocatable :: totcolc                  (:)
  REAL(r8), allocatable :: col_begcb                (:)
  REAL(r8), allocatable :: col_endcb                (:)
  REAL(r8), allocatable :: col_vegbegcb             (:)
  REAL(r8), allocatable :: col_vegendcb             (:)
  REAL(r8), allocatable :: col_soilbegcb            (:)
  REAL(r8), allocatable :: col_soilendcb            (:)

  REAL(r8), allocatable :: totlitn                  (:)
  REAL(r8), allocatable :: totvegn                  (:)
  REAL(r8), allocatable :: totsomn                  (:)
  REAL(r8), allocatable :: totcwdn                  (:)
  REAL(r8), allocatable :: totcoln                  (:)
  REAL(r8), allocatable :: col_begnb                (:)
  REAL(r8), allocatable :: col_endnb                (:)
  REAL(r8), allocatable :: col_vegbegnb             (:)
  REAL(r8), allocatable :: col_vegendnb             (:)
  REAL(r8), allocatable :: col_soilbegnb            (:)
  REAL(r8), allocatable :: col_soilendnb            (:)
  REAL(r8), allocatable :: col_sminnbegnb           (:)
  REAL(r8), allocatable :: col_sminnendnb           (:)

  REAL(r8), allocatable :: leafc                    (:)
  REAL(r8), allocatable :: leafc_storage            (:)
  REAL(r8), allocatable :: leafc_xfer               (:)
  REAL(r8), allocatable :: frootc                   (:)
  REAL(r8), allocatable :: frootc_storage           (:)
  REAL(r8), allocatable :: frootc_xfer              (:)
  REAL(r8), allocatable :: livestemc                (:)
  REAL(r8), allocatable :: livestemc_storage        (:)
  REAL(r8), allocatable :: livestemc_xfer           (:)
  REAL(r8), allocatable :: deadstemc                (:)
  REAL(r8), allocatable :: deadstemc_storage        (:)
  REAL(r8), allocatable :: deadstemc_xfer           (:)
  REAL(r8), allocatable :: livecrootc               (:)
  REAL(r8), allocatable :: livecrootc_storage       (:)
  REAL(r8), allocatable :: livecrootc_xfer          (:)
  REAL(r8), allocatable :: deadcrootc               (:)
  REAL(r8), allocatable :: deadcrootc_storage       (:)
  REAL(r8), allocatable :: deadcrootc_xfer          (:)
  REAL(r8), allocatable :: grainc                   (:)
  REAL(r8), allocatable :: grainc_storage           (:)
  REAL(r8), allocatable :: grainc_xfer              (:)
  REAL(r8), allocatable :: xsmrpool                 (:)
  REAL(r8), allocatable :: downreg                  (:)
  REAL(r8), allocatable :: cropprod1c               (:)
  REAL(r8), allocatable :: cropseedc_deficit        (:)

  REAL(r8), allocatable :: leafn                    (:)
  REAL(r8), allocatable :: leafn_storage            (:)
  REAL(r8), allocatable :: leafn_xfer               (:)
  REAL(r8), allocatable :: frootn                   (:)
  REAL(r8), allocatable :: frootn_storage           (:)
  REAL(r8), allocatable :: frootn_xfer              (:)
  REAL(r8), allocatable :: livestemn                (:)
  REAL(r8), allocatable :: livestemn_storage        (:)
  REAL(r8), allocatable :: livestemn_xfer           (:)
  REAL(r8), allocatable :: deadstemn                (:)
  REAL(r8), allocatable :: deadstemn_storage        (:)
  REAL(r8), allocatable :: deadstemn_xfer           (:)
  REAL(r8), allocatable :: livecrootn               (:)
  REAL(r8), allocatable :: livecrootn_storage       (:)
  REAL(r8), allocatable :: livecrootn_xfer          (:)
  REAL(r8), allocatable :: deadcrootn               (:)
  REAL(r8), allocatable :: deadcrootn_storage       (:)
  REAL(r8), allocatable :: deadcrootn_xfer          (:)
  REAL(r8), allocatable :: grainn                   (:)
  REAL(r8), allocatable :: grainn_storage           (:)
  REAL(r8), allocatable :: grainn_xfer              (:)
  REAL(r8), allocatable :: retransn                 (:)

  REAL(r8), allocatable :: decomp_npools_vr         (:,:,:)
  REAL(r8), allocatable :: decomp_npools            (:,:)
  REAL(r8), allocatable :: ntrunc_vr                (:,:)
  REAL(r8), allocatable :: ntrunc_veg               (:)
  REAL(r8), allocatable :: ntrunc_soil              (:)

  REAL(r8), allocatable :: sminn_vr                 (:,:)
  REAL(r8), allocatable :: smin_no3_vr              (:,:)
  REAL(r8), allocatable :: smin_nh4_vr              (:,:)
  REAL(r8), allocatable :: sminn                    (:)

  REAL(r8), allocatable :: ndep_prof                (:,:)
  REAL(r8), allocatable :: nfixation_prof           (:,:)

  REAL(r8), allocatable :: cn_decomp_pools          (:,:,:) 
  REAL(r8), allocatable :: fpi_vr                   (:,:)
  REAL(r8), allocatable :: fpi                      (:)
  REAL(r8), allocatable :: fpg                      (:)

  REAL(r8), allocatable :: cropf                    (:)
  REAL(r8), allocatable :: lfwt                     (:)
  REAL(r8), allocatable :: fuelc                    (:)
  REAL(r8), allocatable :: fuelc_crop               (:)
  REAL(r8), allocatable :: fsr                      (:)
  REAL(r8), allocatable :: fd                       (:)
  REAL(r8), allocatable :: rootc                    (:)
  REAL(r8), allocatable :: lgdp                     (:)
  REAL(r8), allocatable :: lgdp1                    (:)
  REAL(r8), allocatable :: lpop                     (:)
  REAL(r8), allocatable :: wtlf                     (:)
  REAL(r8), allocatable :: trotr1                   (:)
  REAL(r8), allocatable :: trotr2                   (:)
  REAL(r8), allocatable :: hdmlf                    (:)
  REAL(r8), allocatable :: lnfm                     (:)
  REAL(r8), allocatable :: baf_crop                 (:)
  REAL(r8), allocatable :: baf_peatf                (:)
  REAL(r8), allocatable :: farea_burned             (:)
  REAL(r8), allocatable :: nfire                    (:)
  REAL(r8), allocatable :: fsat                     (:)
  REAL(r8), allocatable :: prec10                   (:) ! 10-day running mean of total      precipitation [mm/s]
  REAL(r8), allocatable :: prec60                   (:) ! 60-day running mean of total      precipitation [mm/s]
  REAL(r8), allocatable :: prec365                  (:) ! 365-day running mean of tota     l precipitation [mm/s]
  REAL(r8), allocatable :: prec_today               (:) ! today's daily precipitation      [mm/day]
  REAL(r8), allocatable :: prec_daily               (:,:) ! daily total precipitation      [mm/day]
  REAL(r8), allocatable :: wf2                      (:)
  REAL(r8), allocatable :: tsoi17                   (:)
  REAL(r8), allocatable :: rh30                     (:) ! 30-day running mean of relative humidity
  REAL(r8), allocatable :: accumnstep               (:) ! 30-day running mean of relative humidity
  REAL(r8), allocatable :: cphase                   (:) ! 30-day running mean of relative humidity

  REAL(r8), allocatable :: dayl                     (:)
  REAL(r8), allocatable :: prev_dayl                (:)

!--------------SASU variables---------------------------
  REAL(r8), allocatable :: decomp0_cpools_vr           (:,:,:)
  REAL(r8), allocatable :: I_met_c_vr_acc              (:,:)
  REAL(r8), allocatable :: I_cel_c_vr_acc              (:,:)
  REAL(r8), allocatable :: I_lig_c_vr_acc              (:,:)
  REAL(r8), allocatable :: I_cwd_c_vr_acc              (:,:)
  REAL(r8), allocatable :: AKX_met_to_soil1_c_vr_acc   (:,:)
  REAL(r8), allocatable :: AKX_cel_to_soil1_c_vr_acc   (:,:)
  REAL(r8), allocatable :: AKX_lig_to_soil2_c_vr_acc   (:,:)
  REAL(r8), allocatable :: AKX_soil1_to_soil2_c_vr_acc (:,:)
  REAL(r8), allocatable :: AKX_cwd_to_cel_c_vr_acc     (:,:)
  REAL(r8), allocatable :: AKX_cwd_to_lig_c_vr_acc     (:,:)
  REAL(r8), allocatable :: AKX_soil1_to_soil3_c_vr_acc (:,:)
  REAL(r8), allocatable :: AKX_soil2_to_soil1_c_vr_acc (:,:)
  REAL(r8), allocatable :: AKX_soil2_to_soil3_c_vr_acc (:,:)
  REAL(r8), allocatable :: AKX_soil3_to_soil1_c_vr_acc (:,:)
  REAL(r8), allocatable :: AKX_met_exit_c_vr_acc       (:,:)
  REAL(r8), allocatable :: AKX_cel_exit_c_vr_acc       (:,:)
  REAL(r8), allocatable :: AKX_lig_exit_c_vr_acc       (:,:)
  REAL(r8), allocatable :: AKX_cwd_exit_c_vr_acc       (:,:)
  REAL(r8), allocatable :: AKX_soil1_exit_c_vr_acc     (:,:)
  REAL(r8), allocatable :: AKX_soil2_exit_c_vr_acc     (:,:)
  REAL(r8), allocatable :: AKX_soil3_exit_c_vr_acc     (:,:)

  REAL(r8), allocatable :: decomp0_npools_vr           (:,:,:)
  REAL(r8), allocatable :: I_met_n_vr_acc              (:,:)
  REAL(r8), allocatable :: I_cel_n_vr_acc              (:,:)
  REAL(r8), allocatable :: I_lig_n_vr_acc              (:,:)
  REAL(r8), allocatable :: I_cwd_n_vr_acc              (:,:)
  REAL(r8), allocatable :: AKX_met_to_soil1_n_vr_acc   (:,:)
  REAL(r8), allocatable :: AKX_cel_to_soil1_n_vr_acc   (:,:)
  REAL(r8), allocatable :: AKX_lig_to_soil2_n_vr_acc   (:,:)
  REAL(r8), allocatable :: AKX_soil1_to_soil2_n_vr_acc (:,:)
  REAL(r8), allocatable :: AKX_cwd_to_cel_n_vr_acc     (:,:)
  REAL(r8), allocatable :: AKX_cwd_to_lig_n_vr_acc     (:,:)
  REAL(r8), allocatable :: AKX_soil1_to_soil3_n_vr_acc (:,:)
  REAL(r8), allocatable :: AKX_soil2_to_soil1_n_vr_acc (:,:)
  REAL(r8), allocatable :: AKX_soil2_to_soil3_n_vr_acc (:,:)
  REAL(r8), allocatable :: AKX_soil3_to_soil1_n_vr_acc (:,:)
  REAL(r8), allocatable :: AKX_met_exit_n_vr_acc       (:,:)
  REAL(r8), allocatable :: AKX_cel_exit_n_vr_acc       (:,:)
  REAL(r8), allocatable :: AKX_lig_exit_n_vr_acc       (:,:)
  REAL(r8), allocatable :: AKX_cwd_exit_n_vr_acc       (:,:)
  REAL(r8), allocatable :: AKX_soil1_exit_n_vr_acc     (:,:)
  REAL(r8), allocatable :: AKX_soil2_exit_n_vr_acc     (:,:)
  REAL(r8), allocatable :: AKX_soil3_exit_n_vr_acc     (:,:)

  REAL(r8), allocatable :: diagVX_c_vr_acc             (:,:,:)
  REAL(r8), allocatable :: upperVX_c_vr_acc            (:,:,:)
  REAL(r8), allocatable :: lowerVX_c_vr_acc            (:,:,:)
  REAL(r8), allocatable :: diagVX_n_vr_acc             (:,:,:)
  REAL(r8), allocatable :: upperVX_n_vr_acc            (:,:,:)
  REAL(r8), allocatable :: lowerVX_n_vr_acc            (:,:,:)
  LOGICAL , allocatable :: skip_balance_check          (:)
!------------------------------------------------------

! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: allocate_TimeVariables
  PUBLIC :: READ_TimeVariables
  PUBLIC :: WRITE_TimeVariables
  PUBLIC :: deallocate_TimeVariables

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_TimeVariables 
! --------------------------------------------------------------------
! Allocates memory for CLM 1d [numpatch] variables
! ------------------------------------------------------

     USE precision
     USE GlobalVars
     USE MOD_PFTimeVars
     USE MOD_PCTimeVars
     IMPLICIT NONE

     allocate (z_sno             (maxsnl+1:0,numpatch))
     allocate (dz_sno            (maxsnl+1:0,numpatch))
     allocate (t_soisno    (maxsnl+1:nl_soil,numpatch))
     allocate (wliq_soisno (maxsnl+1:nl_soil,numpatch))
     allocate (wice_soisno (maxsnl+1:nl_soil,numpatch))
     allocate (smp                (1:nl_soil,numpatch))
     allocate (h2osoi             (1:nl_soil,numpatch))
     allocate (rstfac                       (numpatch))
     allocate (t_grnd                       (numpatch))
     allocate (tleaf                        (numpatch))
     allocate (ldew                         (numpatch))
     allocate (sag                          (numpatch))
     allocate (scv                          (numpatch))
     allocate (snowdp                       (numpatch))
     allocate (fveg                         (numpatch))
     allocate (fsno                         (numpatch))
     allocate (sigf                         (numpatch))
     allocate (green                        (numpatch))
     allocate (tlai                         (numpatch))
     allocate (lai                          (numpatch))
     allocate (laisun                       (numpatch))
     allocate (laisha                       (numpatch))
     allocate (tsai                         (numpatch))
     allocate (sai                          (numpatch))
     allocate (coszen                       (numpatch))
     allocate (alb                      (2,2,numpatch))
     allocate (ssun                     (2,2,numpatch))
     allocate (ssha                     (2,2,numpatch))
     allocate (thermk                       (numpatch))
     allocate (extkb                        (numpatch))
     allocate (extkd                        (numpatch))
     allocate (zwt                          (numpatch))
     allocate (wa                           (numpatch))
     allocate (wat                          (numpatch))

     allocate (t_lake               (nl_lake,numpatch)) 
     allocate (lake_icefrac         (nl_lake,numpatch))

     allocate (trad                         (numpatch))
     allocate (tref                         (numpatch))
     allocate (qref                         (numpatch))
     allocate (rst                          (numpatch))
     allocate (emis                         (numpatch))
     allocate (z0m                          (numpatch))
     allocate (displa                       (numpatch))
     allocate (zol                          (numpatch))
     allocate (rib                          (numpatch))
     allocate (ustar                        (numpatch))
     allocate (qstar                        (numpatch))
     allocate (tstar                        (numpatch))
     allocate (fm                           (numpatch))
     allocate (fh                           (numpatch))
     allocate (fq                           (numpatch))

! bgc variables
     allocate (decomp_cpools_vr             (nl_soil_full,ndecomp_pools,numpatch))
     allocate (decomp_cpools                (ndecomp_pools,numpatch))
     allocate (ctrunc_vr                    (nl_soil_full,numpatch))
     allocate (ctrunc_veg                   (numpatch))
     allocate (ctrunc_soil                  (numpatch))
     allocate (decomp_k                     (nl_soil_full,ndecomp_pools,numpatch))

     allocate (t_scalar                     (nl_soil,numpatch))
     allocate (w_scalar                     (nl_soil,numpatch))
     allocate (o_scalar                     (nl_soil,numpatch))
     allocate (depth_scalar                 (nl_soil,numpatch))

     allocate (som_adv_coef                 (nl_soil_full,numpatch))
     allocate (som_diffus_coef              (nl_soil_full,numpatch))

     allocate (altmax                       (numpatch))
     allocate (altmax_lastyear              (numpatch))
     allocate (altmax_lastyear_indx         (numpatch))

     allocate (totlitc                      (numpatch))
     allocate (totvegc                      (numpatch))
     allocate (totsomc                      (numpatch))
     allocate (totcwdc                      (numpatch))
     allocate (totcolc                      (numpatch))
     allocate (col_begcb                    (numpatch))
     allocate (col_endcb                    (numpatch))
     allocate (col_vegbegcb                 (numpatch))
     allocate (col_vegendcb                 (numpatch))
     allocate (col_soilbegcb                (numpatch))
     allocate (col_soilendcb                (numpatch))

     allocate (totlitn                      (numpatch))
     allocate (totvegn                      (numpatch))
     allocate (totsomn                      (numpatch))
     allocate (totcwdn                      (numpatch))
     allocate (totcoln                      (numpatch))
     allocate (col_begnb                    (numpatch))
     allocate (col_endnb                    (numpatch))
     allocate (col_vegbegnb                 (numpatch))
     allocate (col_vegendnb                 (numpatch))
     allocate (col_soilbegnb                (numpatch))
     allocate (col_soilendnb                (numpatch))
     allocate (col_sminnbegnb               (numpatch))
     allocate (col_sminnendnb               (numpatch))

     allocate (leafc                        (numpatch))
     allocate (leafc_storage                (numpatch))
     allocate (leafc_xfer                   (numpatch))
     allocate (frootc                       (numpatch))
     allocate (frootc_storage               (numpatch))
     allocate (frootc_xfer                  (numpatch))
     allocate (livestemc                    (numpatch))
     allocate (livestemc_storage            (numpatch))
     allocate (livestemc_xfer               (numpatch))
     allocate (deadstemc                    (numpatch))
     allocate (deadstemc_storage            (numpatch))
     allocate (deadstemc_xfer               (numpatch))
     allocate (livecrootc                   (numpatch))
     allocate (livecrootc_storage           (numpatch))
     allocate (livecrootc_xfer              (numpatch))
     allocate (deadcrootc                   (numpatch))
     allocate (deadcrootc_storage           (numpatch))
     allocate (deadcrootc_xfer              (numpatch))
     allocate (grainc                       (numpatch))
     allocate (grainc_storage               (numpatch))
     allocate (grainc_xfer                  (numpatch))
     allocate (xsmrpool                     (numpatch))
     allocate (downreg                      (numpatch))
     allocate (cropprod1c                   (numpatch))
     allocate (cropseedc_deficit            (numpatch))

     allocate (leafn                        (numpatch))
     allocate (leafn_storage                (numpatch))
     allocate (leafn_xfer                   (numpatch))
     allocate (frootn                       (numpatch))
     allocate (frootn_storage               (numpatch))
     allocate (frootn_xfer                  (numpatch))
     allocate (livestemn                    (numpatch))
     allocate (livestemn_storage            (numpatch))
     allocate (livestemn_xfer               (numpatch))
     allocate (deadstemn                    (numpatch))
     allocate (deadstemn_storage            (numpatch))
     allocate (deadstemn_xfer               (numpatch))
     allocate (livecrootn                   (numpatch))
     allocate (livecrootn_storage           (numpatch))
     allocate (livecrootn_xfer              (numpatch))
     allocate (deadcrootn                   (numpatch))
     allocate (deadcrootn_storage           (numpatch))
     allocate (deadcrootn_xfer              (numpatch))
     allocate (grainn                       (numpatch))
     allocate (grainn_storage               (numpatch))
     allocate (grainn_xfer                  (numpatch))
     allocate (retransn                     (numpatch))

     allocate (decomp_npools_vr             (nl_soil_full,ndecomp_pools,numpatch))
     allocate (decomp_npools                (ndecomp_pools,numpatch))
     allocate (ntrunc_vr                    (nl_soil_full,numpatch))
     allocate (ntrunc_veg                   (numpatch))
     allocate (ntrunc_soil                  (numpatch))
     allocate (sminn_vr                     (nl_soil,numpatch))
     allocate (smin_no3_vr                  (nl_soil,numpatch))
     allocate (smin_nh4_vr                  (nl_soil,numpatch))
     allocate (sminn                        (numpatch))

     allocate (ndep_prof                    (nl_soil,numpatch))
     allocate (nfixation_prof               (nl_soil,numpatch))

     allocate (cn_decomp_pools              (nl_soil,ndecomp_pools,numpatch))
     allocate (fpi_vr                       (nl_soil,numpatch))
     allocate (fpi                          (numpatch))
     allocate (fpg                          (numpatch))

     allocate (cropf                        (numpatch))
     allocate (lfwt                         (numpatch))
     allocate (fuelc                        (numpatch))
     allocate (fuelc_crop                   (numpatch))
     allocate (fsr                          (numpatch))
     allocate (fd                           (numpatch))
     allocate (rootc                        (numpatch))
     allocate (lgdp                         (numpatch))
     allocate (lgdp1                        (numpatch))
     allocate (lpop                         (numpatch))
     allocate (wtlf                         (numpatch))
     allocate (trotr1                       (numpatch))
     allocate (trotr2                       (numpatch))
     allocate (hdmlf                        (numpatch))
     allocate (lnfm                         (numpatch))
     allocate (baf_crop                     (numpatch))
     allocate (baf_peatf                    (numpatch))
     allocate (farea_burned                 (numpatch))
     allocate (nfire                        (numpatch))
     allocate (fsat                         (numpatch))
     allocate (prec10                       (numpatch)) ! 10-day running mean of total      precipitation [mm/s]
     allocate (prec60                       (numpatch)) ! 60-day running mean of total      precipitation [mm/s]
     allocate (prec365                      (numpatch)) ! 365-day running mean of tota     l precipitation [mm/s]
     allocate (prec_today                   (numpatch)) ! today's daily precipitation      [mm/day]
     allocate (prec_daily               (365,numpatch)) ! daily total precipitation      [mm/day]
     allocate (wf2                          (numpatch))
     allocate (tsoi17                       (numpatch))
     allocate (rh30                         (numpatch)) ! 30-day running mean of relative humidity
     allocate (accumnstep                   (numpatch)) ! 30-day running mean of relative humidity
     allocate (cphase                       (numpatch)) ! 30-day running mean of relative humidity

     allocate (dayl                         (numpatch))
     allocate (prev_dayl                    (numpatch))
!---------------------------SASU variables--------------------------------------
     allocate (decomp0_cpools_vr            (nl_soil,ndecomp_pools,numpatch))
     allocate (I_met_c_vr_acc               (nl_soil,numpatch))
     allocate (I_cel_c_vr_acc               (nl_soil,numpatch))
     allocate (I_lig_c_vr_acc               (nl_soil,numpatch))
     allocate (I_cwd_c_vr_acc               (nl_soil,numpatch))
     allocate (AKX_met_to_soil1_c_vr_acc    (nl_soil,numpatch))
     allocate (AKX_cel_to_soil1_c_vr_acc    (nl_soil,numpatch))
     allocate (AKX_lig_to_soil2_c_vr_acc    (nl_soil,numpatch))
     allocate (AKX_soil1_to_soil2_c_vr_acc  (nl_soil,numpatch))
     allocate (AKX_cwd_to_cel_c_vr_acc      (nl_soil,numpatch))
     allocate (AKX_cwd_to_lig_c_vr_acc      (nl_soil,numpatch))
     allocate (AKX_soil1_to_soil3_c_vr_acc  (nl_soil,numpatch))
     allocate (AKX_soil2_to_soil1_c_vr_acc  (nl_soil,numpatch))
     allocate (AKX_soil2_to_soil3_c_vr_acc  (nl_soil,numpatch))
     allocate (AKX_soil3_to_soil1_c_vr_acc  (nl_soil,numpatch))
     allocate (AKX_met_exit_c_vr_acc        (nl_soil,numpatch))
     allocate (AKX_cel_exit_c_vr_acc        (nl_soil,numpatch))
     allocate (AKX_lig_exit_c_vr_acc        (nl_soil,numpatch))
     allocate (AKX_cwd_exit_c_vr_acc        (nl_soil,numpatch))
     allocate (AKX_soil1_exit_c_vr_acc      (nl_soil,numpatch))
     allocate (AKX_soil2_exit_c_vr_acc      (nl_soil,numpatch))
     allocate (AKX_soil3_exit_c_vr_acc      (nl_soil,numpatch))

     allocate (decomp0_npools_vr            (nl_soil,ndecomp_pools,numpatch))
     allocate (I_met_n_vr_acc               (nl_soil,numpatch))
     allocate (I_cel_n_vr_acc               (nl_soil,numpatch))
     allocate (I_lig_n_vr_acc               (nl_soil,numpatch))
     allocate (I_cwd_n_vr_acc               (nl_soil,numpatch))
     allocate (AKX_met_to_soil1_n_vr_acc    (nl_soil,numpatch))
     allocate (AKX_cel_to_soil1_n_vr_acc    (nl_soil,numpatch))
     allocate (AKX_lig_to_soil2_n_vr_acc    (nl_soil,numpatch))
     allocate (AKX_soil1_to_soil2_n_vr_acc  (nl_soil,numpatch))
     allocate (AKX_cwd_to_cel_n_vr_acc      (nl_soil,numpatch))
     allocate (AKX_cwd_to_lig_n_vr_acc      (nl_soil,numpatch))
     allocate (AKX_soil1_to_soil3_n_vr_acc  (nl_soil,numpatch))
     allocate (AKX_soil2_to_soil1_n_vr_acc  (nl_soil,numpatch))
     allocate (AKX_soil2_to_soil3_n_vr_acc  (nl_soil,numpatch))
     allocate (AKX_soil3_to_soil1_n_vr_acc  (nl_soil,numpatch))
     allocate (AKX_met_exit_n_vr_acc        (nl_soil,numpatch))
     allocate (AKX_cel_exit_n_vr_acc        (nl_soil,numpatch))
     allocate (AKX_lig_exit_n_vr_acc        (nl_soil,numpatch))
     allocate (AKX_cwd_exit_n_vr_acc        (nl_soil,numpatch))
     allocate (AKX_soil1_exit_n_vr_acc      (nl_soil,numpatch))
     allocate (AKX_soil2_exit_n_vr_acc      (nl_soil,numpatch))
     allocate (AKX_soil3_exit_n_vr_acc      (nl_soil,numpatch))

     allocate (diagVX_c_vr_acc              (nl_soil,ndecomp_pools,numpatch))
     allocate (upperVX_c_vr_acc             (nl_soil,ndecomp_pools,numpatch))
     allocate (lowerVX_c_vr_acc             (nl_soil,ndecomp_pools,numpatch))
     allocate (diagVX_n_vr_acc              (nl_soil,ndecomp_pools,numpatch))
     allocate (upperVX_n_vr_acc             (nl_soil,ndecomp_pools,numpatch))
     allocate (lowerVX_n_vr_acc             (nl_soil,ndecomp_pools,numpatch))

     allocate (skip_balance_check           (numpatch))
!---------------------------------------------------------------------------

#ifdef PFT_CLASSIFICATION
     CALL allocate_PFTimeVars
#endif

#ifdef PC_CLASSIFICATION
     CALL allocate_PCTimeVars
#endif
 
  END SUBROUTINE allocate_TimeVariables


  SUBROUTINE READ_TimeVariables (idate,dir_restart_hist,casename)
! --------------------------------------------------------------------
! Read the model variables for restart run [histTimeVar]
! ...............................................................

     USE precision
     USE MOD_PFTimeVars
     USE MOD_PCTimeVars
     IMPLICIT NONE

     CHARACTER(LEN=255), intent(in) :: dir_restart_hist
     CHARACTER(LEN=256), intent(in) :: casename
     INTEGER, intent(in) :: idate(3)     !calendar (year, julian day, seconds)

     INTEGER :: lhistTimeVar             !logical unit number of restart time-varying file
     INTEGER :: id(3)                    !calendar (year, julian day, seconds)
     CHARACTER(LEN=255) :: cdate         !character for date
     CHARACTER(LEN=256) :: fhistTimeVar  !file name of time-varying file

     ! the model variables for restart run
     write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') idate(1),idate(2),idate(3)

     lhistTimeVar = 100
     fhistTimeVar = trim(dir_restart_hist)//trim(casename)//'-'//'rstTimeVar'//'-'//trim(cdate)
     print*,trim(fhistTimeVar)
     open(unit=lhistTimeVar,file=trim(fhistTimeVar),status='unknown',&
                            form='unformatted',action='read')

     read (lhistTimeVar) id, & !
         ! Time-varying state variables which reaquired by restart run
           z_sno,           &! node depth [m]
           dz_sno,          &! interface depth [m]
           t_soisno,        &! soil temperature [K]
           wliq_soisno,     &! liquid water in layers [kg/m2]
           wice_soisno,     &! ice lens in layers [kg/m2]
           smp,             &! soil matric potential
           t_grnd,          &! ground surface temperature [K]
           tleaf,           &! leaf temperature [K]
           ldew,            &! depth of water on foliage [mm]
           sag,             &! non dimensional snow age [-]
           scv,             &! snow cover, water equivalent [mm]
           snowdp,          &! snow depth [meter]
           fveg,            &! fraction of vegetation cover
           fsno,            &! fraction of snow cover on ground
           sigf,            &! fraction of veg cover, excluding snow-covered veg [-]
           green,           &! leaf greenness
           lai,             &! leaf area index
           tlai,            &! leaf area index
           sai,             &! stem area index
           tsai,            &! stem area index
           coszen,          &! cosine of solar zenith angle
           alb,             &! averaged albedo [-]
           ssun,            &! sunlit canopy absorption for solar radiation (0-1)
           ssha,            &! shaded canopy absorption for solar radiation (0-1)
           thermk,          &! canopy gap fraction for tir radiation
           extkb,           &! (k, g(mu)/mu) direct solar extinction coefficient
           extkd,           &! diffuse and scattered diffuse PAR extinction coefficient
           zwt,             &! the depth to water table [m]
           wa,              &! water storage in aquifer [mm]
 
           t_lake,          &! lake layer teperature [K]
           lake_icefrac,    &! lake mass fraction of lake layer that is frozen

         ! Additional variables required by reginal model (such as WRF & RSM) 
           trad,            &! radiative temperature of surface [K]
           tref,            &! 2 m height air temperature [kelvin]
           qref,            &! 2 m height air specific humidity
           rst,             &! canopy stomatal resistance (s/m)
           emis,            &! averaged bulk surface emissivity
           z0m,             &! effective roughness [m]
           zol,             &! dimensionless height (z/L) used in Monin-Obukhov theory
           rib,             &! bulk Richardson number in surface layer
           ustar,           &! u* in similarity theory [m/s]
           qstar,           &! q* in similarity theory [kg/kg]
           tstar,           &! t* in similarity theory [K]
           fm,              &! integral of profile function for momentum
           fh,              &! integral of profile function for heat
           fq,              &! integral of profile function for moisture

! bgc variables
           totlitc,                &
           totvegc,                &
           totsomc,                &
           totcwdc,                &
           totcolc,                &

           totlitn,                &
           totvegn,                &
           totsomn,                &
           totcwdn,                &
           totcoln,                &
       
           sminn,                  &

           decomp_cpools_vr,    &
           ctrunc_vr,           &
!           decomp_k,            &

!           t_scalar,            &
!           w_scalar,            &
!           o_scalar,            &
!           depth_scalar,        &

!           som_adv_coef,        &
!           som_diffus_coef,     &

           altmax,              &
           altmax_lastyear,     &
           altmax_lastyear_indx,&

!           totlitc,             &
!           totvegc,             &
!           totsomc,             &

           decomp_npools_vr,    &
           ntrunc_vr,          &
           sminn_vr,            &
           smin_no3_vr,         &
           smin_nh4_vr,         &

!           ndep_prof,           &
!           nfixation_prof,      &

!           cn_decomp_pools,     &
!           fpi_vr,              &
!           fpi,                 &
!           fpg,                 &

!           cropf,               &
!           lfwt,                &
!           fuelc,               &
!           fuelc_crop,          &
!           fsr,                 &
!           fd,                  &
!           rootc,               &
!           lgdp,                &
!           lgdp1,               &
!           lpop,                &
!           wtlf,                &
!           trotr1,              &
!           trotr2,              &
!           hdmlf,               &
!           lnfm,                &
!           baf_crop,            &
!           baf_peatf,           &
!           farea_burned,        &
!           nfire,               &
!           fsat,                &
           prec10,              &
           prec60,              &
           prec365,             &
           prec_today,          &
           prec_daily,          &
!           wf2,                 &
           tsoi17,              &
           rh30,                &
           accumnstep,          &
           cphase,          &

!---------------SASU variables-----------------------
           decomp0_cpools_vr            , &
           I_met_c_vr_acc               , &
           I_cel_c_vr_acc               , &
           I_lig_c_vr_acc               , &
           I_cwd_c_vr_acc               , &
           AKX_met_to_soil1_c_vr_acc    , &
           AKX_cel_to_soil1_c_vr_acc    , &
           AKX_lig_to_soil2_c_vr_acc    , &
           AKX_soil1_to_soil2_c_vr_acc  , &
           AKX_cwd_to_cel_c_vr_acc      , &
           AKX_cwd_to_lig_c_vr_acc      , &
           AKX_soil1_to_soil3_c_vr_acc  , &
           AKX_soil2_to_soil1_c_vr_acc  , &
           AKX_soil2_to_soil3_c_vr_acc  , &
           AKX_soil3_to_soil1_c_vr_acc  , &
           AKX_met_exit_c_vr_acc        , &
           AKX_cel_exit_c_vr_acc        , &
           AKX_lig_exit_c_vr_acc        , &
           AKX_cwd_exit_c_vr_acc        , &
           AKX_soil1_exit_c_vr_acc      , &
           AKX_soil2_exit_c_vr_acc      , &
           AKX_soil3_exit_c_vr_acc      , &

           decomp0_npools_vr            , &
           I_met_n_vr_acc               , &
           I_cel_n_vr_acc               , &
           I_lig_n_vr_acc               , &
           I_cwd_n_vr_acc               , &
           AKX_met_to_soil1_n_vr_acc    , &
           AKX_cel_to_soil1_n_vr_acc    , &
           AKX_lig_to_soil2_n_vr_acc    , &
           AKX_soil1_to_soil2_n_vr_acc  , &
           AKX_cwd_to_cel_n_vr_acc      , &
           AKX_cwd_to_lig_n_vr_acc      , &
           AKX_soil1_to_soil3_n_vr_acc  , &
           AKX_soil2_to_soil1_n_vr_acc  , &
           AKX_soil2_to_soil3_n_vr_acc  , &
           AKX_soil3_to_soil1_n_vr_acc  , &
           AKX_met_exit_n_vr_acc        , &
           AKX_cel_exit_n_vr_acc        , &
           AKX_lig_exit_n_vr_acc        , &
           AKX_cwd_exit_n_vr_acc        , &
           AKX_soil1_exit_n_vr_acc      , &
           AKX_soil2_exit_n_vr_acc      , &
           AKX_soil3_exit_n_vr_acc      , &

           diagVX_c_vr_acc              , &
           upperVX_c_vr_acc             , &
           lowerVX_c_vr_acc             , &
           diagVX_n_vr_acc              , &
           upperVX_n_vr_acc             , &
           lowerVX_n_vr_acc             , &

           skip_balance_check           
!----------------------------------------------------
     ! PFT/PC time variabls
#ifdef PFT_CLASSIFICATION
     read (lhistTimeVar)    &!
           tleaf_p,         &! shaded leaf temperature [K]
           ldew_p,          &! depth of water on foliage [mm]
           sigf_p,          &! fraction of veg cover, excluding snow-covered veg [-]
           lai_p,           &! leaf area index
           tlai_p,          &! true leaf area index
           sai_p,           &! stem area index
           tsai_p,          &! true stem area index
           ssun_p,          &! sunlit canopy absorption for solar radiation (0-1)      
           ssha_p,          &! shaded canopy absorption for solar radiation (0-1)
           thermk_p,        &! canopy gap fraction for tir radiation
           extkb_p,         &! (k, g(mu)/mu) direct solar extinction coefficient
           extkd_p,         &! diffuse and scattered diffuse PAR extinction coefficient
           tref_p,          &! 2 m height air temperature [kelvin]
           qref_p,          &! 2 m height air specific humidity
           rst_p,           &! canopy stomatal resistance (s/m)
           z0m_p,           &! effective roughness [m]                                 

! bgc PFT variables
           leafc_p,                &
           leafc_storage_p,        &
           leafc_xfer_p,           &
           frootc_p,               &
           frootc_storage_p,       &
           frootc_xfer_p,          &
           livestemc_p,            &
           livestemc_storage_p,    &
           livestemc_xfer_p,       &
           deadstemc_p,            &
           deadstemc_storage_p,    &
           deadstemc_xfer_p,       &
           livecrootc_p,           &
           livecrootc_storage_p,   &
           livecrootc_xfer_p,      &
           deadcrootc_p,           &
           deadcrootc_storage_p,   &
           deadcrootc_xfer_p,      &
           grainc_p,               &
           grainc_storage_p,       &
           grainc_xfer_p,          &
           cropseedc_deficit_p,    &
           xsmrpool_p,             &
           gresp_storage_p,        &
           gresp_xfer_p,           &
           cpool_p,                &
           totvegc_p,              &
           cropprod1c_p,           &
           cropseedc_deficit,      &
      
           leafn_p,                &
           leafn_storage_p,        &
           leafn_xfer_p,           &
           frootn_p,               &
           frootn_storage_p,       &
           frootn_xfer_p,          &
           livestemn_p,            &
           livestemn_storage_p,    &
           livestemn_xfer_p,       &
           deadstemn_p,            &
           deadstemn_storage_p,    &
           deadstemn_xfer_p,       &
           livecrootn_p,           &
           livecrootn_storage_p,   &
           livecrootn_xfer_p,      &
           deadcrootn_p,           &
           deadcrootn_storage_p,   &
           deadcrootn_xfer_p,      &
           grainn_p,               &
           grainn_storage_p,       &
           grainn_xfer_p,          &
           cropseedn_deficit_p,    &
           retransn_p,             &
           totvegn_p,              &
      
           harvdate_p,             &
      
           tempsum_potential_gpp_p,&
           tempmax_retransn_p,     &
           tempavg_tref_p,         &
           tempsum_npp_p,          &
           tempsum_litfall_p,      &
           annsum_potential_gpp_p, &
           annmax_retransn_p,      &
           annavg_tref_p,          &
           annsum_npp_p,           &
           annsum_litfall_p,       &
      
           bglfr_p,                &
           bgtr_p,                 &
           lgsf_p,                 &
           gdd0_p,                 &
           gdd8_p,                 &
           gdd10_p,                &
           gdd020_p,               &
           gdd820_p,               &
           gdd1020_p,              &
           nyrs_crop_active_p,     &
      
           offset_flag_p,          &
           offset_counter_p,       &
           onset_flag_p,           &
           onset_counter_p,        &
           onset_gddflag_p,        &
           onset_gdd_p,            &
           onset_fdd_p,            &
           onset_swi_p,            &
           offset_fdd_p,           &
           offset_swi_p,           &
           dormant_flag_p,         &
           prev_leafc_to_litter_p, &
           prev_frootc_to_litter_p,&
           days_active_p,          &
      
           burndate_p,             &
           grain_flag_p,           &
           ctrunc_p,               &
           ntrunc_p,               &
           npool_p,                &

! crop variables 
           croplive_p,             &
           gddtsoi_p,              &
           huileaf_p,              &
           gddplant_p,             &
           huigrain_p,             &
           peaklai_p,              &
           aroot_p,                &
           astem_p,                &
           arepr_p,                &
           aleaf_p,                &
           astemi_p,               &
           aleafi_p,               &
           gddmaturity_p,          &

           cropplant_p,            &
           idop_p,                 &
           a5tmin_p,               &
           a10tmin_p,              &
           t10_p,                  &
           cumvd_p,                &
           hdidx_p,                &
           vf_p,                   &
           cphase_p,               &
           fert_counter_p,         &
           fert_p,                 &
           tref_min_p,             &
           tref_max_p,             &
           tref_min_inst_p,        &
           tref_max_inst_p,        &
           fertnitro_p,            &
           latbaset_p,             &

! SASU variables
           leafc0_p,                &
           leafc0_storage_p,        &
           leafc0_xfer_p,           &
           frootc0_p,               &
           frootc0_storage_p,       &
           frootc0_xfer_p,          &
           livestemc0_p,            &
           livestemc0_storage_p,    &
           livestemc0_xfer_p,       &
           deadstemc0_p,            &
           deadstemc0_storage_p,    &
           deadstemc0_xfer_p,       &
           livecrootc0_p,           &
           livecrootc0_storage_p,   &
           livecrootc0_xfer_p,      &
           deadcrootc0_p,           &
           deadcrootc0_storage_p,   &
           deadcrootc0_xfer_p,      &
           grainc0_p,               &
           grainc0_storage_p,       &
           grainc0_xfer_p,          &

           leafn0_p,                &
           leafn0_storage_p,        &
           leafn0_xfer_p,           &
           frootn0_p,               &
           frootn0_storage_p,       &
           frootn0_xfer_p,          &
           livestemn0_p,            &
           livestemn0_storage_p,    &
           livestemn0_xfer_p,       &
           deadstemn0_p,            &
           deadstemn0_storage_p,    &
           deadstemn0_xfer_p,       &
           livecrootn0_p,           &
           livecrootn0_storage_p,   &
           livecrootn0_xfer_p,      &
           deadcrootn0_p,           &
           deadcrootn0_storage_p,   &
           deadcrootn0_xfer_p,      &
           grainn0_p,               &
           grainn0_storage_p,       &
           grainn0_xfer_p,          &
           retransn0_p,             &

           I_leafc_p_acc         , &
           I_leafc_st_p_acc      , &
           I_frootc_p_acc        , &
           I_frootc_st_p_acc     , &
           I_livestemc_p_acc     , &
           I_livestemc_st_p_acc  , &
           I_deadstemc_p_acc     , &
           I_deadstemc_st_p_acc  , &
           I_livecrootc_p_acc    , &
           I_livecrootc_st_p_acc , &
           I_deadcrootc_p_acc    , &
           I_deadcrootc_st_p_acc , &
           I_grainc_p_acc        , &
           I_grainc_st_p_acc     , &
           I_leafn_p_acc         , &
           I_leafn_st_p_acc      , &
           I_frootn_p_acc        , &
           I_frootn_st_p_acc     , &
           I_livestemn_p_acc     , &
           I_livestemn_st_p_acc  , &
           I_deadstemn_p_acc     , &
           I_deadstemn_st_p_acc  , &
           I_livecrootn_p_acc    , &
           I_livecrootn_st_p_acc , &
           I_deadcrootn_p_acc    , &
           I_deadcrootn_st_p_acc , &
           I_grainn_p_acc        , &
           I_grainn_st_p_acc     , &

           AKX_leafc_xf_to_leafc_p_acc              , &
           AKX_frootc_xf_to_frootc_p_acc            , &
           AKX_livestemc_xf_to_livestemc_p_acc      , &
           AKX_deadstemc_xf_to_deadstemc_p_acc      , &
           AKX_livecrootc_xf_to_livecrootc_p_acc    , &
           AKX_deadcrootc_xf_to_deadcrootc_p_acc    , &
           AKX_grainc_xf_to_grainc_p_acc            , &
           AKX_livestemc_to_deadstemc_p_acc         , &
           AKX_livecrootc_to_deadcrootc_p_acc       , &
           
           AKX_leafc_st_to_leafc_xf_p_acc           , &
           AKX_frootc_st_to_frootc_xf_p_acc         , &
           AKX_livestemc_st_to_livestemc_xf_p_acc   , &
           AKX_deadstemc_st_to_deadstemc_xf_p_acc   , &
           AKX_livecrootc_st_to_livecrootc_xf_p_acc , &
           AKX_deadcrootc_st_to_deadcrootc_xf_p_acc , &
           AKX_grainc_st_to_grainc_xf_p_acc         , &

           AKX_leafc_exit_p_acc                     , &
           AKX_frootc_exit_p_acc                    , &
           AKX_livestemc_exit_p_acc                 , &
           AKX_deadstemc_exit_p_acc                 , &
           AKX_livecrootc_exit_p_acc                , &
           AKX_deadcrootc_exit_p_acc                , &
           AKX_grainc_exit_p_acc                    , &

           AKX_leafc_st_exit_p_acc                  , &
           AKX_frootc_st_exit_p_acc                 , &
           AKX_livestemc_st_exit_p_acc              , &
           AKX_deadstemc_st_exit_p_acc              , &
           AKX_livecrootc_st_exit_p_acc             , &
           AKX_deadcrootc_st_exit_p_acc             , &
           AKX_grainc_st_exit_p_acc                 , &

           AKX_leafc_xf_exit_p_acc                  , &
           AKX_frootc_xf_exit_p_acc                 , &
           AKX_livestemc_xf_exit_p_acc              , &
           AKX_deadstemc_xf_exit_p_acc              , &
           AKX_livecrootc_xf_exit_p_acc             , &
           AKX_deadcrootc_xf_exit_p_acc             , &
           AKX_grainc_xf_exit_p_acc                 , &
           
           AKX_leafn_xf_to_leafn_p_acc              , &
           AKX_frootn_xf_to_frootn_p_acc            , &
           AKX_livestemn_xf_to_livestemn_p_acc      , &
           AKX_deadstemn_xf_to_deadstemn_p_acc      , &
           AKX_livecrootn_xf_to_livecrootn_p_acc    , &
           AKX_deadcrootn_xf_to_deadcrootn_p_acc    , &
           AKX_grainn_xf_to_grainn_p_acc            , &
           AKX_livestemn_to_deadstemn_p_acc         , &
           AKX_livecrootn_to_deadcrootn_p_acc       , &

           AKX_leafn_st_to_leafn_xf_p_acc           , &
           AKX_frootn_st_to_frootn_xf_p_acc         , &
           AKX_livestemn_st_to_livestemn_xf_p_acc   , &
           AKX_deadstemn_st_to_deadstemn_xf_p_acc   , &
           AKX_livecrootn_st_to_livecrootn_xf_p_acc , &
           AKX_deadcrootn_st_to_deadcrootn_xf_p_acc , &
           AKX_grainn_st_to_grainn_xf_p_acc         , &

           AKX_leafn_to_retransn_p_acc              , &
           AKX_frootn_to_retransn_p_acc             , &
           AKX_livestemn_to_retransn_p_acc          , &
           AKX_livecrootn_to_retransn_p_acc         , &

           AKX_retransn_to_leafn_p_acc              , &
           AKX_retransn_to_frootn_p_acc             , &
           AKX_retransn_to_livestemn_p_acc          , &
           AKX_retransn_to_deadstemn_p_acc          , &
           AKX_retransn_to_livecrootn_p_acc         , &
           AKX_retransn_to_deadcrootn_p_acc         , &
           AKX_retransn_to_grainn_p_acc             , &

           AKX_retransn_to_leafn_st_p_acc           , &
           AKX_retransn_to_frootn_st_p_acc          , &
           AKX_retransn_to_livestemn_st_p_acc       , &
           AKX_retransn_to_deadstemn_st_p_acc       , &
           AKX_retransn_to_livecrootn_st_p_acc      , &
           AKX_retransn_to_deadcrootn_st_p_acc      , &
           AKX_retransn_to_grainn_st_p_acc          , &

           AKX_leafn_exit_p_acc                     , &
           AKX_frootn_exit_p_acc                    , &
           AKX_livestemn_exit_p_acc                 , &
           AKX_deadstemn_exit_p_acc                 , &
           AKX_livecrootn_exit_p_acc                , &
           AKX_deadcrootn_exit_p_acc                , &
           AKX_grainn_exit_p_acc                    , &
           AKX_retransn_exit_p_acc                  , &

           AKX_leafn_st_exit_p_acc                  , &
           AKX_frootn_st_exit_p_acc                 , &
           AKX_livestemn_st_exit_p_acc              , &
           AKX_deadstemn_st_exit_p_acc              , &
           AKX_livecrootn_st_exit_p_acc             , &
           AKX_deadcrootn_st_exit_p_acc             , &
           AKX_grainn_st_exit_p_acc                 , &

           AKX_leafn_xf_exit_p_acc                  , &
           AKX_frootn_xf_exit_p_acc                 , &
           AKX_livestemn_xf_exit_p_acc              , &
           AKX_deadstemn_xf_exit_p_acc              , &
           AKX_livecrootn_xf_exit_p_acc             , &
           AKX_deadcrootn_xf_exit_p_acc             , &
           AKX_grainn_xf_exit_p_acc                 
#endif

#ifdef PC_CLASSIFICATION
     read (lhistTimeVar)    &!
           tleaf_c,         &! leaf temperature [K]
           ldew_c,          &! depth of water on foliage [mm]
           sigf_c,          &! fraction of veg cover, excluding snow-covered veg [-]
           lai_c,           &! leaf area index
           tlai_c,          &! true leaf area index
           sai_c,           &! stem area index
           tsai_c,          &! true stem area index
           ssun_c,          &! sunlit canopy absorption for solar radiation (0-1)
           ssha_c,          &! shaded canopy absorption for solar radiation (0-1)
           thermk_c,        &! canopy gap fraction for tir radiation
           fshade_c,        &! canopy gap fraction for tir radiation
           extkb_c,         &! (k, g(mu)/mu) direct solar extinction coefficient
           extkd_c,         &! diffuse and scattered diffuse PAR extinction coefficient
           rst_c,           &! canopy stomatal resistance (s/m)
           z0m_c             ! effective roughness [m]                                 
#endif
      
           IF (id(1) /= idate(1) .or. id(2) /= idate(2) .or. id(3) /= idate(3)) THEN
              print*, 'id = ', id, 'idate = ', idate
              print*, 'The date of initial data is NOT IDENTICAL TO initial set-up'
!              CALL abort
           ENDIF

     close(lhistTimeVar)

  END SUBROUTINE READ_TimeVariables


  SUBROUTINE WRITE_TimeVariables (idate,dir_restart_hist,casename)
! --------------------------------------------------------------------
! Write out the model variables for restart run [histTimeVar]
! --------------------------------------------------------------------

     USE precision
     USE MOD_PFTimeVars
     USE MOD_PCTimeVars
     IMPLICIT NONE
     
     INTEGER, intent(in) :: idate(3)     !calendar (year, julian day, seconds)
     CHARACTER(LEN=255), intent(in) :: dir_restart_hist
     CHARACTER(LEN=256), intent(in) :: casename

     INTEGER :: lhistTimeVar             !logical unit number of restart time-varying file
     INTEGER :: id(3)                    !calendar (year, julian day, seconds), temporal
     CHARACTER(LEN=255) :: cdate         !character for date
     CHARACTER(LEN=256) :: fhistTimeVar  !file name of time-varying file

! ...............................................................

     id(:) = idate(:)
     CALL adj2begin(id)

     ! the model variables for restart run 
     write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') id(1), id(2), id(3)

     lhistTimeVar = 100
     fhistTimeVar = trim(dir_restart_hist)//trim(casename)//'-'//'rstTimeVar'//'-'//trim(cdate)
     print*,trim(fhistTimeVar)
     open(unit=lhistTimeVar,file=trim(fhistTimeVar),status='unknown',&
                            form='unformatted',action='write')

     write(lhistTimeVar) id, & !
         ! Time-varying state variables which reaquired by restart run
           z_sno,           &! node depth [m]
           dz_sno,          &! interface depth [m]
           t_soisno,        &! soil temperature [K]
           wliq_soisno,     &! liquid water in layers [kg/m2]
           wice_soisno,     &! ice lens in layers [kg/m2]
           smp,             &! soil matric potential
           t_grnd,          &! ground surface temperature [K]
           tleaf,           &! leaf temperature [K]
           ldew,            &! depth of water on foliage [mm]
           sag,             &! non dimensional snow age [-]
           scv,             &! snow cover, water equivalent [mm]
           snowdp,          &! snow depth [meter]
           fveg,            &! fraction of vegetation cover
           fsno,            &! fraction of snow cover on ground
           sigf,            &! fraction of veg cover, excluding snow-covered veg [-]
           green,           &! leaf greenness
           lai,             &! leaf area index
           tlai,            &! leaf area index
           sai,             &! stem area index
           tsai,            &! stem area index
           coszen,          &! cosine of solar zenith angle
           alb,             &! averaged albedo [-]
           ssun,            &! sunlit canopy absorption for solar radiation (0-1)
           ssha,            &! shaded canopy absorption for solar radiation (0-1)
           thermk,          &! canopy gap fraction for tir radiation
           extkb,           &! (k, g(mu)/mu) direct solar extinction coefficient
           extkd,           &! diffuse and scattered diffuse PAR extinction coefficient
           zwt,             &! the depth to water table [m]
           wa,              &! water storage in aquifer [mm]
 
           t_lake,          &! lake layer teperature [K]
           lake_icefrac,    &! lake mass fraction of lake layer that is frozen

         ! Additional variables required by reginal model (such as WRF & RSM) 
           trad,            &! radiative temperature of surface [K]
           tref,            &! 2 m height air temperature [kelvin]
           qref,            &! 2 m height air specific humidity
           rst,             &! canopy stomatal resistance (s/m)
           emis,            &! averaged bulk surface emissivity
           z0m,             &! effective roughness [m]
           zol,             &! dimensionless height (z/L) used in Monin-Obukhov theory
           rib,             &! bulk Richardson number in surface layer
           ustar,           &! u* in similarity theory [m/s]
           qstar,           &! q* in similarity theory [kg/kg]
           tstar,           &! t* in similarity theory [K]
           fm,              &! integral of profile function for momentum
           fh,              &! integral of profile function for heat
           fq,              &! integral of profile function for moisture

! bgc variables
           totlitc,                &
           totvegc,                &
           totsomc,                &
           totcwdc,                &
           totcolc,                &

           totlitn,                &
           totvegn,                &
           totsomn,                &
           totcwdn,                &
           totcoln,                &

           sminn,                  &

           decomp_cpools_vr,    &
           ctrunc_vr,           &
!           decomp_k,            &

!           t_scalar,            &
!           w_scalar,            &
!           o_scalar,            &
!           depth_scalar,        &

!           som_adv_coef,        &
!           som_diffus_coef,     &

           altmax,              &
           altmax_lastyear,     &
           altmax_lastyear_indx,&

!           totlitc,             &
!           totvegc,             &
!           totsomc,             &

           decomp_npools_vr,    &
           ntrunc_vr,           &
           sminn_vr,            &
           smin_no3_vr,         &
           smin_nh4_vr,         &

!           ndep_prof,           &
!           nfixation_prof,      &

!           cn_decomp_pools,     &
!           fpi_vr,              &
!           fpi,                 &
!           fpg,                 &

!           cropf,               &
!           lfwt,                &
!           fuelc,               &
!           fuelc_crop,          &
!           fsr,                 &
!           fd,                  &
!           rootc,               &
!           lgdp,                &
!           lgdp1,               &
!           lpop,                &
!           wtlf,                &
!           trotr1,              &
!           trotr2,              &
!           hdmlf,               &
!           lnfm,                &
!           baf_crop,            &
!           baf_peatf,           &
!           farea_burned,        &
!           nfire,               &
!           fsat,                &
           prec10,              &
           prec60,              &
           prec365,             &
           prec_today,          &
           prec_daily,          &
!           wf2,                 &
           tsoi17,              &
           rh30,                &
           accumnstep,             &
           cphase,          &

!---------------SASU variables-----------------------
           decomp0_cpools_vr            , &
           I_met_c_vr_acc               , &
           I_cel_c_vr_acc               , &
           I_lig_c_vr_acc               , &
           I_cwd_c_vr_acc               , &
           AKX_met_to_soil1_c_vr_acc    , &
           AKX_cel_to_soil1_c_vr_acc    , &
           AKX_lig_to_soil2_c_vr_acc    , &
           AKX_soil1_to_soil2_c_vr_acc  , &
           AKX_cwd_to_cel_c_vr_acc      , &
           AKX_cwd_to_lig_c_vr_acc      , &
           AKX_soil1_to_soil3_c_vr_acc  , &
           AKX_soil2_to_soil1_c_vr_acc  , &
           AKX_soil2_to_soil3_c_vr_acc  , &
           AKX_soil3_to_soil1_c_vr_acc  , &
           AKX_met_exit_c_vr_acc        , &
           AKX_cel_exit_c_vr_acc        , &
           AKX_lig_exit_c_vr_acc        , &
           AKX_cwd_exit_c_vr_acc        , &
           AKX_soil1_exit_c_vr_acc      , &
           AKX_soil2_exit_c_vr_acc      , &
           AKX_soil3_exit_c_vr_acc      , &

           decomp0_npools_vr            , &
           I_met_n_vr_acc               , &
           I_cel_n_vr_acc               , &
           I_lig_n_vr_acc               , &
           I_cwd_n_vr_acc               , &
           AKX_met_to_soil1_n_vr_acc    , &
           AKX_cel_to_soil1_n_vr_acc    , &
           AKX_lig_to_soil2_n_vr_acc    , &
           AKX_soil1_to_soil2_n_vr_acc  , &
           AKX_cwd_to_cel_n_vr_acc      , &
           AKX_cwd_to_lig_n_vr_acc      , &
           AKX_soil1_to_soil3_n_vr_acc  , &
           AKX_soil2_to_soil1_n_vr_acc  , &
           AKX_soil2_to_soil3_n_vr_acc  , &
           AKX_soil3_to_soil1_n_vr_acc  , &
           AKX_met_exit_n_vr_acc        , &
           AKX_cel_exit_n_vr_acc        , &
           AKX_lig_exit_n_vr_acc        , &
           AKX_cwd_exit_n_vr_acc        , &
           AKX_soil1_exit_n_vr_acc      , &
           AKX_soil2_exit_n_vr_acc      , &
           AKX_soil3_exit_n_vr_acc      , &

           diagVX_c_vr_acc              , &
           upperVX_c_vr_acc             , &
           lowerVX_c_vr_acc             , &
           diagVX_n_vr_acc              , &
           upperVX_n_vr_acc             , &
           lowerVX_n_vr_acc             , &

           skip_balance_check           
!----------------------------------------------------
      ! PFT/PC time variabls
#ifdef PFT_CLASSIFICATION
     write(lhistTimeVar)    &!
           tleaf_p,         &! shaded leaf temperature [K]
           ldew_p,          &! depth of water on foliage [mm]
           sigf_p,          &! fraction of veg cover, excluding snow-covered veg [-]
           lai_p,           &! leaf area index
           tlai_p,          &! true leaf area index
           sai_p,           &! stem area index
           tsai_p,          &! true stem area index
           ssun_p,          &! sunlit canopy absorption for solar radiation (0-1)       
           ssha_p,          &! shaded canopy absorption for solar radiation (0-1)
           thermk_p,        &! canopy gap fraction for tir radiation
           extkb_p,         &! (k, g(mu)/mu) direct solar extinction coefficient
           extkd_p,         &! diffuse and scattered diffuse PAR extinction coefficient
           tref_p,          &! 2 m height air temperature [kelvin]
           qref_p,          &! 2 m height air specific humidity
           rst_p,           &! canopy stomatal resistance (s/m)
           z0m_p,           &! effective roughness [m]                                 
! bgc PFT variables
           leafc_p,                &
           leafc_storage_p,        &
           leafc_xfer_p,           &
           frootc_p,               &
           frootc_storage_p,       &
           frootc_xfer_p,          &
           livestemc_p,            &
           livestemc_storage_p,    &
           livestemc_xfer_p,       &
           deadstemc_p,            &
           deadstemc_storage_p,    &
           deadstemc_xfer_p,       &
           livecrootc_p,           &
           livecrootc_storage_p,   &
           livecrootc_xfer_p,      &
           deadcrootc_p,           &
           deadcrootc_storage_p,   &
           deadcrootc_xfer_p,      &
           grainc_p,               &
           grainc_storage_p,       &
           grainc_xfer_p,          &
           cropseedc_deficit_p,    &
           xsmrpool_p,             &
           gresp_storage_p,        &
           gresp_xfer_p,           &
           cpool_p,                &
           totvegc_p,              &
           cropprod1c_p,           &
           cropseedc_deficit,      &
      
           leafn_p,                &
           leafn_storage_p,        &
           leafn_xfer_p,           &
           frootn_p,               &
           frootn_storage_p,       &
           frootn_xfer_p,          &
           livestemn_p,            &
           livestemn_storage_p,    &
           livestemn_xfer_p,       &
           deadstemn_p,            &
           deadstemn_storage_p,    &
           deadstemn_xfer_p,       &
           livecrootn_p,           &
           livecrootn_storage_p,   &
           livecrootn_xfer_p,      &
           deadcrootn_p,           &
           deadcrootn_storage_p,   &
           deadcrootn_xfer_p,      &
           grainn_p,               &
           grainn_storage_p,       &
           grainn_xfer_p,          &
           cropseedn_deficit_p,    &
           retransn_p,             &
           totvegn_p,              &
      
           harvdate_p,             &
      
           tempsum_potential_gpp_p,&
           tempmax_retransn_p,     &
           tempavg_tref_p,         &
           tempsum_npp_p,          &
           tempsum_litfall_p,      &
           annsum_potential_gpp_p, &
           annmax_retransn_p,      &
           annavg_tref_p,          &
           annsum_npp_p,           &
           annsum_litfall_p,       &
      
           bglfr_p,                &
           bgtr_p,                 &
           lgsf_p,                 &
           gdd0_p,                 &
           gdd8_p,                 &
           gdd10_p,                &
           gdd020_p,               &
           gdd820_p,               &
           gdd1020_p,              &
           nyrs_crop_active_p,     &
      
           offset_flag_p,          &
           offset_counter_p,       &
           onset_flag_p,           &
           onset_counter_p,        &
           onset_gddflag_p,        &
           onset_gdd_p,            &
           onset_fdd_p,            &
           onset_swi_p,            &
           offset_fdd_p,           &
           offset_swi_p,           &
           dormant_flag_p,         &
           prev_leafc_to_litter_p, &
           prev_frootc_to_litter_p,&
           days_active_p,          &
      
           burndate_p,             &
           grain_flag_p,           &
           ctrunc_p,               &
           ntrunc_p,               &
           npool_p,                &

! crop variables 
           croplive_p,             &
           gddtsoi_p,              &
           huileaf_p,              &
           gddplant_p,             &
           huigrain_p,             &
           peaklai_p,              &
           aroot_p,                &
           astem_p,                &
           arepr_p,                &
           aleaf_p,                &
           astemi_p,               &
           aleafi_p,               &
           gddmaturity_p,          &

           cropplant_p,            &
           idop_p,                 &
           a5tmin_p,               &
           a10tmin_p,              &
           t10_p,                  &
           cumvd_p,                &
           hdidx_p,                &
           vf_p,                   &
           cphase_p,               &
           fert_counter_p,         &
           fert_p,                 &
           tref_min_p,             &
           tref_max_p,             &
           tref_min_inst_p,        &
           tref_max_inst_p,        &
           fertnitro_p,            &
           latbaset_p,             &

! SASU variables
           leafc0_p,                &
           leafc0_storage_p,        &
           leafc0_xfer_p,           &
           frootc0_p,               &
           frootc0_storage_p,       &
           frootc0_xfer_p,          &
           livestemc0_p,            &
           livestemc0_storage_p,    &
           livestemc0_xfer_p,       &
           deadstemc0_p,            &
           deadstemc0_storage_p,    &
           deadstemc0_xfer_p,       &
           livecrootc0_p,           &
           livecrootc0_storage_p,   &
           livecrootc0_xfer_p,      &
           deadcrootc0_p,           &
           deadcrootc0_storage_p,   &
           deadcrootc0_xfer_p,      &
           grainc0_p,               &
           grainc0_storage_p,       &
           grainc0_xfer_p,          &

           leafn0_p,                &
           leafn0_storage_p,        &
           leafn0_xfer_p,           &
           frootn0_p,               &
           frootn0_storage_p,       &
           frootn0_xfer_p,          &
           livestemn0_p,            &
           livestemn0_storage_p,    &
           livestemn0_xfer_p,       &
           deadstemn0_p,            &
           deadstemn0_storage_p,    &
           deadstemn0_xfer_p,       &
           livecrootn0_p,           &
           livecrootn0_storage_p,   &
           livecrootn0_xfer_p,      &
           deadcrootn0_p,           &
           deadcrootn0_storage_p,   &
           deadcrootn0_xfer_p,      &
           grainn0_p,               &
           grainn0_storage_p,       &
           grainn0_xfer_p,          &
           retransn0_p,             &

           I_leafc_p_acc         , &
           I_leafc_st_p_acc      , &
           I_frootc_p_acc        , &
           I_frootc_st_p_acc     , &
           I_livestemc_p_acc     , &
           I_livestemc_st_p_acc  , &
           I_deadstemc_p_acc     , &
           I_deadstemc_st_p_acc  , &
           I_livecrootc_p_acc    , &
           I_livecrootc_st_p_acc , &
           I_deadcrootc_p_acc    , &
           I_deadcrootc_st_p_acc , &
           I_grainc_p_acc        , &
           I_grainc_st_p_acc     , &
           I_leafn_p_acc         , &
           I_leafn_st_p_acc      , &
           I_frootn_p_acc        , &
           I_frootn_st_p_acc     , &
           I_livestemn_p_acc     , &
           I_livestemn_st_p_acc  , &
           I_deadstemn_p_acc     , &
           I_deadstemn_st_p_acc  , &
           I_livecrootn_p_acc    , &
           I_livecrootn_st_p_acc , &
           I_deadcrootn_p_acc    , &
           I_deadcrootn_st_p_acc , &
           I_grainn_p_acc        , &
           I_grainn_st_p_acc     , &

           AKX_leafc_xf_to_leafc_p_acc              , &
           AKX_frootc_xf_to_frootc_p_acc            , &
           AKX_livestemc_xf_to_livestemc_p_acc      , &
           AKX_deadstemc_xf_to_deadstemc_p_acc      , &
           AKX_livecrootc_xf_to_livecrootc_p_acc    , &
           AKX_deadcrootc_xf_to_deadcrootc_p_acc    , &
           AKX_grainc_xf_to_grainc_p_acc            , &
           AKX_livestemc_to_deadstemc_p_acc         , &
           AKX_livecrootc_to_deadcrootc_p_acc       , &
           
           AKX_leafc_st_to_leafc_xf_p_acc           , &
           AKX_frootc_st_to_frootc_xf_p_acc         , &
           AKX_livestemc_st_to_livestemc_xf_p_acc   , &
           AKX_deadstemc_st_to_deadstemc_xf_p_acc   , &
           AKX_livecrootc_st_to_livecrootc_xf_p_acc , &
           AKX_deadcrootc_st_to_deadcrootc_xf_p_acc , &
           AKX_grainc_st_to_grainc_xf_p_acc         , &

           AKX_leafc_exit_p_acc                     , &
           AKX_frootc_exit_p_acc                    , &
           AKX_livestemc_exit_p_acc                 , &
           AKX_deadstemc_exit_p_acc                 , &
           AKX_livecrootc_exit_p_acc                , &
           AKX_deadcrootc_exit_p_acc                , &
           AKX_grainc_exit_p_acc                    , &

           AKX_leafc_st_exit_p_acc                  , &
           AKX_frootc_st_exit_p_acc                 , &
           AKX_livestemc_st_exit_p_acc              , &
           AKX_deadstemc_st_exit_p_acc              , &
           AKX_livecrootc_st_exit_p_acc             , &
           AKX_deadcrootc_st_exit_p_acc             , &
           AKX_grainc_st_exit_p_acc                 , &

           AKX_leafc_xf_exit_p_acc                  , &
           AKX_frootc_xf_exit_p_acc                 , &
           AKX_livestemc_xf_exit_p_acc              , &
           AKX_deadstemc_xf_exit_p_acc              , &
           AKX_livecrootc_xf_exit_p_acc             , &
           AKX_deadcrootc_xf_exit_p_acc             , &
           AKX_grainc_xf_exit_p_acc                 , &
           
           AKX_leafn_xf_to_leafn_p_acc              , &
           AKX_frootn_xf_to_frootn_p_acc            , &
           AKX_livestemn_xf_to_livestemn_p_acc      , &
           AKX_deadstemn_xf_to_deadstemn_p_acc      , &
           AKX_livecrootn_xf_to_livecrootn_p_acc    , &
           AKX_deadcrootn_xf_to_deadcrootn_p_acc    , &
           AKX_grainn_xf_to_grainn_p_acc            , &
           AKX_livestemn_to_deadstemn_p_acc         , &
           AKX_livecrootn_to_deadcrootn_p_acc       , &

           AKX_leafn_st_to_leafn_xf_p_acc           , &
           AKX_frootn_st_to_frootn_xf_p_acc         , &
           AKX_livestemn_st_to_livestemn_xf_p_acc   , &
           AKX_deadstemn_st_to_deadstemn_xf_p_acc   , &
           AKX_livecrootn_st_to_livecrootn_xf_p_acc , &
           AKX_deadcrootn_st_to_deadcrootn_xf_p_acc , &
           AKX_grainn_st_to_grainn_xf_p_acc         , &

           AKX_leafn_to_retransn_p_acc              , &
           AKX_frootn_to_retransn_p_acc             , &
           AKX_livestemn_to_retransn_p_acc          , &
           AKX_livecrootn_to_retransn_p_acc         , &

           AKX_retransn_to_leafn_p_acc              , &
           AKX_retransn_to_frootn_p_acc             , &
           AKX_retransn_to_livestemn_p_acc          , &
           AKX_retransn_to_deadstemn_p_acc          , &
           AKX_retransn_to_livecrootn_p_acc         , &
           AKX_retransn_to_deadcrootn_p_acc         , &
           AKX_retransn_to_grainn_p_acc             , &

           AKX_retransn_to_leafn_st_p_acc           , &
           AKX_retransn_to_frootn_st_p_acc          , &
           AKX_retransn_to_livestemn_st_p_acc       , &
           AKX_retransn_to_deadstemn_st_p_acc       , &
           AKX_retransn_to_livecrootn_st_p_acc      , &
           AKX_retransn_to_deadcrootn_st_p_acc      , &
           AKX_retransn_to_grainn_st_p_acc          , &

           AKX_leafn_exit_p_acc                     , &
           AKX_frootn_exit_p_acc                    , &
           AKX_livestemn_exit_p_acc                 , &
           AKX_deadstemn_exit_p_acc                 , &
           AKX_livecrootn_exit_p_acc                , &
           AKX_deadcrootn_exit_p_acc                , &
           AKX_grainn_exit_p_acc                    , &
           AKX_retransn_exit_p_acc                  , &

           AKX_leafn_st_exit_p_acc                  , &
           AKX_frootn_st_exit_p_acc                 , &
           AKX_livestemn_st_exit_p_acc              , &
           AKX_deadstemn_st_exit_p_acc              , &
           AKX_livecrootn_st_exit_p_acc             , &
           AKX_deadcrootn_st_exit_p_acc             , &
           AKX_grainn_st_exit_p_acc                 , &

           AKX_leafn_xf_exit_p_acc                  , &
           AKX_frootn_xf_exit_p_acc                 , &
           AKX_livestemn_xf_exit_p_acc              , &
           AKX_deadstemn_xf_exit_p_acc              , &
           AKX_livecrootn_xf_exit_p_acc             , &
           AKX_deadcrootn_xf_exit_p_acc             , &
           AKX_grainn_xf_exit_p_acc                 
#endif

#ifdef PC_CLASSIFICATION
     WRITE(lhistTimeVar)    &!
           tleaf_c,         &! leaf temperature [K]
           ldew_c,          &! depth of water on foliage [mm]
           sigf_c,          &! fraction of veg cover, excluding snow-covered veg [-]
           lai_c,           &! leaf area index
           tlai_c,          &! true leaf area index
           sai_c,           &! stem area index
           tsai_c,          &! true stem area index
           ssun_c,          &! sunlit canopy absorption for solar radiation (0-1)
           ssha_c,          &! shaded canopy absorption for solar radiation (0-1)
           thermk_c,        &! canopy gap fraction for tir radiation
           fshade_c,        &! canopy gap fraction for tir radiation
           extkb_c,         &! (k, g(mu)/mu) direct solar extinction coefficient
           extkd_c,         &! diffuse and scattered diffuse PAR extinction coefficient
           rst_c,           &! canopy stomatal resistance (s/m)
           z0m_c             ! effective roughness [m]                                 
#endif
 
     close(lhistTimeVar)

  END SUBROUTINE WRITE_TimeVariables

  SUBROUTINE deallocate_TimeVariables
! --------------------------------------------------
! Deallocates memory for CLM 1d [numpatch] variables
! --------------------------------------------------
     USE MOD_PFTimeVars
     USE MOD_PCTimeVars

     deallocate (z_sno        )
     deallocate (dz_sno       )
     deallocate (t_soisno     )
     deallocate (wliq_soisno  )
     deallocate (wice_soisno  )
     deallocate (smp          )
     deallocate (h2osoi       )
     deallocate (rstfac       )
     deallocate (t_grnd       )
     deallocate (tleaf        )
     deallocate (ldew         )
     deallocate (sag          )
     deallocate (scv          )
     deallocate (snowdp       )
     deallocate (fveg         )
     deallocate (fsno         )
     deallocate (sigf         )
     deallocate (green        )
     deallocate (tlai         )
     deallocate (lai          )
     deallocate (laisun       )
     deallocate (laisha       )
     deallocate (tsai         )
     deallocate (sai          )
     deallocate (coszen       )
     deallocate (alb          )
     deallocate (ssun         )
     deallocate (ssha         )
     deallocate (thermk       )
     deallocate (extkb        )
     deallocate (extkd        )
     deallocate (zwt          )
     deallocate (wa           )
     deallocate (wat          )

     deallocate (t_lake       ) 
     deallocate (lake_icefrac )

     deallocate (trad         )
     deallocate (tref         )
     deallocate (qref         )
     deallocate (rst          )
     deallocate (emis         )
     deallocate (z0m          )
     deallocate (displa       )
     deallocate (zol          )
     deallocate (rib          )
     deallocate (ustar        )
     deallocate (qstar        )
     deallocate (tstar        )
     deallocate (fm           )
     deallocate (fh           )
     deallocate (fq           )
     
! bgc variables
     deallocate (decomp_cpools_vr     )
     deallocate (decomp_cpools        )
     deallocate (ctrunc_vr            )
     deallocate (ctrunc_veg           )
     deallocate (ctrunc_soil          )
     deallocate (decomp_k             )

     deallocate (t_scalar             )
     deallocate (w_scalar             )
     deallocate (o_scalar             )
     deallocate (depth_scalar         )

     deallocate (som_adv_coef         )
     deallocate (som_diffus_coef      )

     deallocate (altmax               )
     deallocate (altmax_lastyear      )
     deallocate (altmax_lastyear_indx )

     deallocate (totlitc              )
     deallocate (totvegc              )
     deallocate (totsomc              )
     deallocate (totcwdc              )
     deallocate (totcolc              )
     deallocate (col_begcb            )
     deallocate (col_endcb            )
     deallocate (col_vegbegcb         )
     deallocate (col_vegendcb         )
     deallocate (col_soilbegcb        )
     deallocate (col_soilendcb        )

     deallocate (totlitn              )
     deallocate (totvegn              )
     deallocate (totsomn              )
     deallocate (totcwdn              )
     deallocate (totcoln              )
     deallocate (col_begnb            )
     deallocate (col_endnb            )
     deallocate (col_vegbegnb         )
     deallocate (col_vegendnb         )
     deallocate (col_soilbegnb        )
     deallocate (col_soilendnb        )
     deallocate (col_sminnbegnb       )
     deallocate (col_sminnendnb       )

     deallocate (leafc                )
     deallocate (leafc_storage        )
     deallocate (leafc_xfer           )
     deallocate (frootc               )
     deallocate (frootc_storage       )
     deallocate (frootc_xfer          )
     deallocate (livestemc            )
     deallocate (livestemc_storage    )
     deallocate (livestemc_xfer       )
     deallocate (deadstemc            )
     deallocate (deadstemc_storage    )
     deallocate (deadstemc_xfer       )
     deallocate (livecrootc           )
     deallocate (livecrootc_storage   )
     deallocate (livecrootc_xfer      )
     deallocate (deadcrootc           )
     deallocate (deadcrootc_storage   )
     deallocate (deadcrootc_xfer      )
     deallocate (grainc               )
     deallocate (grainc_storage       )
     deallocate (grainc_xfer          )
     deallocate (xsmrpool             )
     deallocate (downreg              )
     deallocate (cropprod1c           )
     deallocate (cropseedc_deficit    )

     deallocate (leafn                )
     deallocate (leafn_storage        )
     deallocate (leafn_xfer           )
     deallocate (frootn               )
     deallocate (frootn_storage       )
     deallocate (frootn_xfer          )
     deallocate (livestemn            )
     deallocate (livestemn_storage    )
     deallocate (livestemn_xfer       )
     deallocate (deadstemn            )
     deallocate (deadstemn_storage    )
     deallocate (deadstemn_xfer       )
     deallocate (livecrootn           )
     deallocate (livecrootn_storage   )
     deallocate (livecrootn_xfer      )
     deallocate (deadcrootn           )
     deallocate (deadcrootn_storage   )
     deallocate (deadcrootn_xfer      )
     deallocate (grainn               )
     deallocate (grainn_storage       )
     deallocate (grainn_xfer          )
     deallocate (retransn             )

     deallocate (decomp_npools_vr     )
     deallocate (decomp_npools        )
     deallocate (ntrunc_vr            )
     deallocate (ntrunc_veg           )
     deallocate (ntrunc_soil          )
     deallocate (sminn_vr             )
     deallocate (smin_no3_vr          )
     deallocate (smin_nh4_vr          )
     deallocate (sminn                )

     deallocate (ndep_prof            )
     deallocate (nfixation_prof       )

     deallocate (cn_decomp_pools      )
     deallocate (fpi_vr               )
     deallocate (fpi                  )
     deallocate (fpg                  )

     deallocate (cropf                )
     deallocate (lfwt                 )
     deallocate (fuelc                )
     deallocate (fuelc_crop           )
     deallocate (fsr                  )
     deallocate (fd                   )
     deallocate (rootc                )
     deallocate (lgdp                 )
     deallocate (lgdp1                )
     deallocate (lpop                 )
     deallocate (wtlf                 )
     deallocate (trotr1               )
     deallocate (trotr2               )
     deallocate (hdmlf                )
     deallocate (lnfm                 )
     deallocate (baf_crop             )
     deallocate (baf_peatf            )
     deallocate (farea_burned         )
     deallocate (nfire                )
     deallocate (fsat                 )
     deallocate (prec10               )
     deallocate (prec60               )
     deallocate (prec365              )
     deallocate (prec_today           )
     deallocate (prec_daily           )
     deallocate (wf2                  )
     deallocate (tsoi17               )
     deallocate (rh30                 )
     deallocate (accumnstep           )
     deallocate (cphase               )

     deallocate (dayl                 )
     deallocate (prev_dayl            )

!---------------------------SASU variables--------------------------------------
     deallocate (decomp0_cpools_vr           )
     deallocate (I_met_c_vr_acc              )
     deallocate (I_cel_c_vr_acc              )
     deallocate (I_lig_c_vr_acc              )
     deallocate (I_cwd_c_vr_acc              )
     deallocate (AKX_met_to_soil1_c_vr_acc   )
     deallocate (AKX_cel_to_soil1_c_vr_acc   )
     deallocate (AKX_lig_to_soil2_c_vr_acc   )
     deallocate (AKX_soil1_to_soil2_c_vr_acc )
     deallocate (AKX_cwd_to_cel_c_vr_acc     )
     deallocate (AKX_cwd_to_lig_c_vr_acc     )
     deallocate (AKX_soil1_to_soil3_c_vr_acc )
     deallocate (AKX_soil2_to_soil1_c_vr_acc )
     deallocate (AKX_soil2_to_soil3_c_vr_acc )
     deallocate (AKX_soil3_to_soil1_c_vr_acc )
     deallocate (AKX_met_exit_c_vr_acc       )
     deallocate (AKX_cel_exit_c_vr_acc       )
     deallocate (AKX_lig_exit_c_vr_acc       )
     deallocate (AKX_cwd_exit_c_vr_acc       )
     deallocate (AKX_soil1_exit_c_vr_acc     )
     deallocate (AKX_soil2_exit_c_vr_acc     )
     deallocate (AKX_soil3_exit_c_vr_acc     )

     deallocate (decomp0_npools_vr           )
     deallocate (I_met_n_vr_acc              )
     deallocate (I_cel_n_vr_acc              )
     deallocate (I_lig_n_vr_acc              )
     deallocate (I_cwd_n_vr_acc              )
     deallocate (AKX_met_to_soil1_n_vr_acc   )
     deallocate (AKX_cel_to_soil1_n_vr_acc   )
     deallocate (AKX_lig_to_soil2_n_vr_acc   )
     deallocate (AKX_soil1_to_soil2_n_vr_acc )
     deallocate (AKX_cwd_to_cel_n_vr_acc     )
     deallocate (AKX_cwd_to_lig_n_vr_acc     )
     deallocate (AKX_soil1_to_soil3_n_vr_acc )
     deallocate (AKX_soil2_to_soil1_n_vr_acc )
     deallocate (AKX_soil2_to_soil3_n_vr_acc )
     deallocate (AKX_soil3_to_soil1_n_vr_acc )
     deallocate (AKX_met_exit_n_vr_acc       )
     deallocate (AKX_cel_exit_n_vr_acc       )
     deallocate (AKX_lig_exit_n_vr_acc       )
     deallocate (AKX_cwd_exit_n_vr_acc       )
     deallocate (AKX_soil1_exit_n_vr_acc     )
     deallocate (AKX_soil2_exit_n_vr_acc     )
     deallocate (AKX_soil3_exit_n_vr_acc     )

     deallocate (diagVX_c_vr_acc             )
     deallocate (upperVX_c_vr_acc            )
     deallocate (lowerVX_c_vr_acc            )
     deallocate (diagVX_n_vr_acc             )
     deallocate (upperVX_n_vr_acc            )
     deallocate (lowerVX_n_vr_acc            )

     deallocate (skip_balance_check          )
!---------------------------------------------------------------------------
#ifdef PFT_CLASSIFICATION
     CALL deallocate_PFTimeVars
#endif

#ifdef PC_CLASSIFICATION
     CALL deallocate_PCTimeVars
#endif
  
  END SUBROUTINE deallocate_TimeVariables

END MODULE MOD_TimeVariables
! ---------- EOP ------------
