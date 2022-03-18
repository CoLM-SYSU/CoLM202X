#include <define.h>

MODULE MOD_TimeVariables
! -------------------------------
! Created by Yongjiu Dai, 03/2014
! -------------------------------

use precision
use timemanager
IMPLICIT NONE
SAVE
! -----------------------------------------------------------------
! Time-varying state variables which reaquired by restart run
      real(r8), allocatable :: z_sno    (:,:)   ! node depth [m]
      real(r8), allocatable :: dz_sno   (:,:)   ! interface depth [m]
      real(r8), allocatable :: t_soisno (:,:)   ! soil temperature [K]
      real(r8), allocatable :: wliq_soisno(:,:) ! liquid water in layers [kg/m2]
      real(r8), allocatable :: wice_soisno(:,:) ! ice lens in layers [kg/m2]
      real(r8), allocatable :: h2osoi (:,:)     ! volumetric soil water in layers [m3/m3]
      real(r8), allocatable :: smp(:,:)         ! soil matrix potential [mm]
      real(r8), allocatable :: hk (:,:)         ! hydraulic conductivity [mm h2o/s]
      real(r8), allocatable :: rootr(:,:)       ! water exchange between soil and root. Positive: soil->root [?]
#ifdef PLANT_HYDRAULIC_STRESS
      real(r8), allocatable :: vegwp(:,:)       ! vegetation water potential [mm]
      real(r8), allocatable :: gs0sun   (:)     ! working copy of sunlit stomata conductance
      real(r8), allocatable :: gs0sha   (:)     ! working copy of shalit stomata conductance
#endif
      real(r8), allocatable :: rstfacsun(:)     ! factor of soil water stress on sunlit leaf
      real(r8), allocatable :: rstfacsha(:)     ! factor of soil water stress on shaded leaf
      real(r8), allocatable :: t_grnd   (:)     ! ground surface temperature [K]

      real(r8), allocatable :: tleaf    (:)     ! leaf temperature [K]
      real(r8), allocatable :: ldew     (:)     ! depth of water on foliage [mm]
      real(r8), allocatable :: sag      (:)     ! non dimensional snow age [-]
      real(r8), allocatable :: scv      (:)     ! snow cover, water equivalent [mm]
      real(r8), allocatable :: snowdp   (:)     ! snow depth [meter]
      real(r8), allocatable :: fveg     (:)     ! fraction of vegetation cover
      real(r8), allocatable :: fsno     (:)     ! fraction of snow cover on ground
      real(r8), allocatable :: sigf     (:)     ! fraction of veg cover, excluding snow-covered veg [-]
      real(r8), allocatable :: green    (:)     ! leaf greenness
      real(r8), allocatable :: tlai     (:)     ! leaf area index
      real(r8), allocatable :: lai      (:)     ! leaf area index
      real(r8), allocatable :: laisun   (:)     ! leaf area index
      real(r8), allocatable :: laisha   (:)     ! leaf area index
      real(r8), allocatable :: tsai     (:)     ! stem area index
      real(r8), allocatable :: sai      (:)     ! stem area index
      real(r8), allocatable :: coszen   (:)     ! cosine of solar zenith angle
      real(r8), allocatable :: alb  (:,:,:)     ! averaged albedo [-]
      real(r8), allocatable :: ssun (:,:,:)     ! sunlit canopy absorption for solar radiation (0-1)
      real(r8), allocatable :: ssha (:,:,:)     ! shaded canopy absorption for solar radiation (0-1)
      real(r8), allocatable :: thermk   (:)     ! canopy gap fraction for tir radiation
      real(r8), allocatable :: extkb    (:)     ! (k, g(mu)/mu) direct solar extinction coefficient
      real(r8), allocatable :: extkd    (:)     ! diffuse and scattered diffuse PAR extinction coefficient
      real(r8), allocatable :: zwt      (:)     ! the depth to water table [m]
      real(r8), allocatable :: wa       (:)     ! water storage in aquifer [mm]
      real(r8), allocatable :: wat      (:)     ! total water storage [mm]
#ifdef VARIABLY_SATURATED_FLOW
      real(r8), allocatable :: dpond    (:)     ! depth of ponding water
#endif

      real(r8), allocatable :: t_lake(:,:)      ! lake layer teperature [K]
      real(r8), allocatable :: lake_icefrac(:,:)! lake mass fraction of lake layer that is frozen
      real(r8), allocatable :: savedtke1(:)     ! top level eddy conductivity (W/m K) 

      real(r8), allocatable :: trad     (:) ! radiative temperature of surface [K]
      real(r8), allocatable :: tref     (:) ! 2 m height air temperature [kelvin]
      real(r8), allocatable :: qref     (:) ! 2 m height air specific humidity
      real(r8), allocatable :: rst      (:) ! canopy stomatal resistance (s/m)
      real(r8), allocatable :: emis     (:) ! averaged bulk surface emissivity
      real(r8), allocatable :: z0m      (:) ! effective roughness [m]
      real(r8), allocatable :: displa   (:) ! zero displacement height [m]
      real(r8), allocatable :: zol      (:) ! dimensionless height (z/L) used in Monin-Obukhov theory
      real(r8), allocatable :: rib      (:) ! bulk Richardson number in surface layer
      real(r8), allocatable :: ustar    (:) ! u* in similarity theory [m/s]
      real(r8), allocatable :: qstar    (:) ! q* in similarity theory [kg/kg]
      real(r8), allocatable :: tstar    (:) ! t* in similarity theory [K]
      real(r8), allocatable :: fm       (:) ! integral of profile function for momentum
      real(r8), allocatable :: fh       (:) ! integral of profile function for heat
      real(r8), allocatable :: fq       (:) ! integral of profile function for moisture

!------------------------- BGC variables -------------------------------
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

!--------------BGC/SASU variables---------------------------
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
      public :: allocate_TimeVariables
      public :: deallocate_TimeVariables
      public :: READ_TimeVariables
      public :: WRITE_TimeVariables
#ifdef CLMDEBUG
      public :: check_TimeVariables
#endif

! PRIVATE MEMBER FUNCTIONS:


!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_TimeVariables 
! --------------------------------------------------------------------
! Allocates memory for CLM 1d [numpatch] variables
! ------------------------------------------------------

  use precision
  USE GlobalVars
#ifdef PFT_CLASSIFICATION
  USE MOD_PFTimeVars
#endif
#ifdef PC_CLASSIFICATION
  USE MOD_PCTimeVars
#endif
  use spmd_task
  use mod_landpatch, only : numpatch
  IMPLICIT NONE


  if (p_is_worker) then

     if (numpatch > 0) then

        allocate (z_sno      (maxsnl+1:0,      numpatch))
        allocate (dz_sno     (maxsnl+1:0,      numpatch))
        allocate (t_soisno   (maxsnl+1:nl_soil,numpatch))
        allocate (wliq_soisno(maxsnl+1:nl_soil,numpatch))
        allocate (wice_soisno(maxsnl+1:nl_soil,numpatch))
        allocate (smp        (1:nl_soil,numpatch))
        allocate (hk         (1:nl_soil,numpatch))
        allocate (h2osoi     (1:nl_soil,numpatch))
        allocate (rootr      (1:nl_soil,numpatch))
#ifdef PLANT_HYDRAULIC_STRESS
        allocate (vegwp      (1:nvegwcs,numpatch))
        allocate (gs0sun               (numpatch))
        allocate (gs0sha               (numpatch))
#endif
        allocate (rstfacsun            (numpatch))
        allocate (rstfacsha            (numpatch))
        allocate (t_grnd               (numpatch))
        allocate (tleaf                (numpatch))
        allocate (ldew                 (numpatch))
        allocate (sag                  (numpatch))
        allocate (scv                  (numpatch))
        allocate (snowdp               (numpatch))
        allocate (fveg                 (numpatch))
        allocate (fsno                 (numpatch))
        allocate (sigf                 (numpatch))
        allocate (green                (numpatch))
        allocate (tlai                 (numpatch))
        allocate (lai                  (numpatch))
        allocate (laisun               (numpatch))
        allocate (laisha               (numpatch))
        allocate (tsai                 (numpatch))
        allocate (sai                  (numpatch))
        allocate (coszen               (numpatch))
        allocate (alb              (2,2,numpatch))
        allocate (ssun             (2,2,numpatch))
        allocate (ssha             (2,2,numpatch))
        allocate (thermk               (numpatch))
        allocate (extkb                (numpatch))
        allocate (extkd                (numpatch))
        allocate (zwt                  (numpatch))
        allocate (wa                   (numpatch))
        allocate (wat                  (numpatch))
#ifdef VARIABLY_SATURATED_FLOW
        allocate (dpond                (numpatch))
#endif

        allocate (t_lake       (nl_lake,numpatch))    !new lake scheme
        allocate (lake_icefrac (nl_lake,numpatch))    !new lake scheme
        allocate (savedtke1            (numpatch))    !new lake scheme

        allocate (trad                 (numpatch))
        allocate (tref                 (numpatch))
        allocate (qref                 (numpatch))
        allocate (rst                  (numpatch))
        allocate (emis                 (numpatch))
        allocate (z0m                  (numpatch))
        allocate (displa               (numpatch))
        allocate (zol                  (numpatch))
        allocate (rib                  (numpatch))
        allocate (ustar                (numpatch))
        allocate (qstar                (numpatch))
        allocate (tstar                (numpatch))
        allocate (fm                   (numpatch))
        allocate (fh                   (numpatch))
        allocate (fq                   (numpatch))

#ifdef BGC
! bgc variables
        allocate (decomp_cpools_vr             (nl_soil_full,ndecomp_pools,numpatch))
        allocate (decomp_cpools                (ndecomp_pools,numpatch))
        allocate (ctrunc_vr                    (nl_soil,numpatch))
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
        allocate (ntrunc_vr                    (nl_soil,numpatch))
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

#ifdef SASU
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

!---------------------------------------------------------------------------
#endif
        allocate (skip_balance_check           (numpatch))
#endif
     end if
  end if


#ifdef PFT_CLASSIFICATION
     CALL allocate_PFTimeVars
#endif

#ifdef PC_CLASSIFICATION
     CALL allocate_PCTimeVars
#endif

  END SUBROUTINE allocate_TimeVariables



  SUBROUTINE deallocate_TimeVariables ()

     use spmd_task
     use mod_landpatch, only : numpatch
#if (defined PFT_CLASSIFICATION)
     USE MOD_PFTimeVars
#endif
#if (defined PC_CLASSIFICATION)
     USE MOD_PCTimeVars
#endif
     implicit none

     ! --------------------------------------------------
     ! Deallocates memory for CLM 1d [numpatch] variables
     ! --------------------------------------------------

     if (p_is_worker) then
        
        if (numpatch > 0) then

           deallocate (z_sno    )
           deallocate (dz_sno   )
           deallocate (t_soisno    )
           deallocate (wliq_soisno )
           deallocate (wice_soisno )
           deallocate (smp )
           deallocate (hk  )
           deallocate (h2osoi )
           deallocate (rootr  )
           deallocate (rstfacsun )
           deallocate (rstfacsha )
#ifdef PLANT_HYDRAULIC_STRESS 
           deallocate (vegwp  )
           deallocate (gs0sun )
           deallocate (gs0sha )
#endif
           deallocate (t_grnd )
           deallocate (tleaf  )
           deallocate (ldew   )
           deallocate (sag    )
           deallocate (scv    )
           deallocate (snowdp )
           deallocate (fveg   )
           deallocate (fsno   )
           deallocate (sigf   )
           deallocate (green  )
           deallocate (tlai   )
           deallocate (lai    )
           deallocate (laisun )
           deallocate (laisha )
           deallocate (tsai   )
           deallocate (sai    )
           deallocate (coszen )
           deallocate (alb    )
           deallocate (ssun   )
           deallocate (ssha   )
           deallocate (thermk )
           deallocate (extkb  )
           deallocate (extkd  )
           deallocate (zwt    )
           deallocate (wa     )
           deallocate (wat    )
#ifdef VARIABLY_SATURATED_FLOW
           deallocate (dpond  )
#endif

           deallocate (t_lake )      ! new lake scheme
           deallocate (lake_icefrac) ! new lake scheme
           deallocate (savedtke1)    ! new lake scheme

           deallocate (trad   )
           deallocate (tref   )
           deallocate (qref   )
           deallocate (rst    )
           deallocate (emis   )
           deallocate (z0m    )
           deallocate (displa )
           deallocate (zol    )
           deallocate (rib    )
           deallocate (ustar  )
           deallocate (qstar  )
           deallocate (tstar  )
           deallocate (fm     )
           deallocate (fh     )
           deallocate (fq     )

#ifdef BGC
! bgc variables
           deallocate (decomp_cpools_vr             )
           deallocate (decomp_cpools                )
           deallocate (ctrunc_vr                    )
           deallocate (ctrunc_veg                   )
           deallocate (ctrunc_soil                  )
           deallocate (decomp_k                     )
   
           deallocate (t_scalar                     )
           deallocate (w_scalar                     )
           deallocate (o_scalar                     )
           deallocate (depth_scalar                 )

           deallocate (som_adv_coef                 )
           deallocate (som_diffus_coef              )

           deallocate (altmax                       )
           deallocate (altmax_lastyear              )
           deallocate (altmax_lastyear_indx         )

           deallocate (totlitc                      )
           deallocate (totvegc                      )
           deallocate (totsomc                      )
           deallocate (totcwdc                      )
           deallocate (totcolc                      )
           deallocate (col_begcb                    )
           deallocate (col_endcb                    )
           deallocate (col_vegbegcb                 )
           deallocate (col_vegendcb                 )
           deallocate (col_soilbegcb                )
           deallocate (col_soilendcb                )

           deallocate (totlitn                      )
           deallocate (totvegn                      )
           deallocate (totsomn                      )
           deallocate (totcwdn                      )
           deallocate (totcoln                      )
           deallocate (col_begnb                    )
           deallocate (col_endnb                    )
           deallocate (col_vegbegnb                 )
           deallocate (col_vegendnb                 )
           deallocate (col_soilbegnb                )
           deallocate (col_soilendnb                )
           deallocate (col_sminnbegnb               )
           deallocate (col_sminnendnb               )

           deallocate (leafc                        )
           deallocate (leafc_storage                )
           deallocate (leafc_xfer                   )
           deallocate (frootc                       )
           deallocate (frootc_storage               )
           deallocate (frootc_xfer                  )
           deallocate (livestemc                    )
           deallocate (livestemc_storage            )
           deallocate (livestemc_xfer               )
           deallocate (deadstemc                    )
           deallocate (deadstemc_storage            )
           deallocate (deadstemc_xfer               )
           deallocate (livecrootc                   )
           deallocate (livecrootc_storage           )
           deallocate (livecrootc_xfer              )
           deallocate (deadcrootc                   )
           deallocate (deadcrootc_storage           )
           deallocate (deadcrootc_xfer              )
           deallocate (grainc                       )
           deallocate (grainc_storage               )
           deallocate (grainc_xfer                  )
           deallocate (xsmrpool                     )
           deallocate (downreg                      )
           deallocate (cropprod1c                   )
           deallocate (cropseedc_deficit            )

           deallocate (leafn                        )
           deallocate (leafn_storage                )
           deallocate (leafn_xfer                   )
           deallocate (frootn                       )
           deallocate (frootn_storage               )
           deallocate (frootn_xfer                  )
           deallocate (livestemn                    )
           deallocate (livestemn_storage            )
           deallocate (livestemn_xfer               )
           deallocate (deadstemn                    )
           deallocate (deadstemn_storage            )
           deallocate (deadstemn_xfer               )
           deallocate (livecrootn                   )
           deallocate (livecrootn_storage           )
           deallocate (livecrootn_xfer              )
           deallocate (deadcrootn                   )
           deallocate (deadcrootn_storage           )
           deallocate (deadcrootn_xfer              )
           deallocate (grainn                       )
           deallocate (grainn_storage               )
           deallocate (grainn_xfer                  )
           deallocate (retransn                     )

           deallocate (decomp_npools_vr             )
           deallocate (decomp_npools                )
           deallocate (ntrunc_vr                    )
           deallocate (ntrunc_veg                   )
           deallocate (ntrunc_soil                  )
           deallocate (sminn_vr                     )
           deallocate (smin_no3_vr                  )
           deallocate (smin_nh4_vr                  )
           deallocate (sminn                        )
   
           deallocate (ndep_prof                    )
           deallocate (nfixation_prof               )

           deallocate (cn_decomp_pools              )
           deallocate (fpi_vr                       )
           deallocate (fpi                          )
           deallocate (fpg                          )

           deallocate (cropf                        )
           deallocate (lfwt                         )
           deallocate (fuelc                        )
           deallocate (fuelc_crop                   )
           deallocate (fsr                          )
           deallocate (fd                           )
           deallocate (rootc                        )
           deallocate (lgdp                         )
           deallocate (lgdp1                        )
           deallocate (lpop                         )
           deallocate (wtlf                         )
           deallocate (trotr1                       )
           deallocate (trotr2                       )
           deallocate (hdmlf                        )
           deallocate (lnfm                         )
           deallocate (baf_crop                     )
           deallocate (baf_peatf                    )
           deallocate (farea_burned                 )
           deallocate (nfire                        )
           deallocate (fsat                         )
           deallocate (prec10                       )
           deallocate (prec60                       )
           deallocate (prec365                      )
           deallocate (prec_today                   )
           deallocate (prec_daily                   )
           deallocate (wf2                          )
           deallocate (tsoi17                       )
           deallocate (rh30                         )
           deallocate (accumnstep                   )
           deallocate (cphase                       )

           deallocate (dayl                         )
           deallocate (prev_dayl                    )

#ifdef SASU
!---------------------------SASU variables--------------------------------------
           deallocate (decomp0_cpools_vr            )
           deallocate (I_met_c_vr_acc               )
           deallocate (I_cel_c_vr_acc               )
           deallocate (I_lig_c_vr_acc               )
           deallocate (I_cwd_c_vr_acc               )
           deallocate (AKX_met_to_soil1_c_vr_acc    )
           deallocate (AKX_cel_to_soil1_c_vr_acc    )
           deallocate (AKX_lig_to_soil2_c_vr_acc    )
           deallocate (AKX_soil1_to_soil2_c_vr_acc  )
           deallocate (AKX_cwd_to_cel_c_vr_acc      )
           deallocate (AKX_cwd_to_lig_c_vr_acc      )
           deallocate (AKX_soil1_to_soil3_c_vr_acc  )
           deallocate (AKX_soil2_to_soil1_c_vr_acc  )
           deallocate (AKX_soil2_to_soil3_c_vr_acc  )
           deallocate (AKX_soil3_to_soil1_c_vr_acc  )
           deallocate (AKX_met_exit_c_vr_acc        )
           deallocate (AKX_cel_exit_c_vr_acc        )
           deallocate (AKX_lig_exit_c_vr_acc        )
           deallocate (AKX_cwd_exit_c_vr_acc        )
           deallocate (AKX_soil1_exit_c_vr_acc      )
           deallocate (AKX_soil2_exit_c_vr_acc      )
           deallocate (AKX_soil3_exit_c_vr_acc      )

           deallocate (decomp0_npools_vr            )
           deallocate (I_met_n_vr_acc               )
           deallocate (I_cel_n_vr_acc               )
           deallocate (I_lig_n_vr_acc               )
           deallocate (I_cwd_n_vr_acc               )
           deallocate (AKX_met_to_soil1_n_vr_acc    )
           deallocate (AKX_cel_to_soil1_n_vr_acc    )
           deallocate (AKX_lig_to_soil2_n_vr_acc    )
           deallocate (AKX_soil1_to_soil2_n_vr_acc  )
           deallocate (AKX_cwd_to_cel_n_vr_acc      )
           deallocate (AKX_cwd_to_lig_n_vr_acc      )
           deallocate (AKX_soil1_to_soil3_n_vr_acc  )
           deallocate (AKX_soil2_to_soil1_n_vr_acc  )
           deallocate (AKX_soil2_to_soil3_n_vr_acc  )
           deallocate (AKX_soil3_to_soil1_n_vr_acc  )
           deallocate (AKX_met_exit_n_vr_acc        )
           deallocate (AKX_cel_exit_n_vr_acc        )
           deallocate (AKX_lig_exit_n_vr_acc        )
           deallocate (AKX_cwd_exit_n_vr_acc        )
           deallocate (AKX_soil1_exit_n_vr_acc      )
           deallocate (AKX_soil2_exit_n_vr_acc      )
           deallocate (AKX_soil3_exit_n_vr_acc      )

           deallocate (diagVX_c_vr_acc              )
           deallocate (upperVX_c_vr_acc             )
           deallocate (lowerVX_c_vr_acc             )
           deallocate (diagVX_n_vr_acc              )
           deallocate (upperVX_n_vr_acc             )
           deallocate (lowerVX_n_vr_acc             )

!---------------------------------------------------------------------------
#endif
           deallocate (skip_balance_check           )
#endif
        end if
     end if

#if (defined PFT_CLASSIFICATION)
     CALL deallocate_PFTimeVars
#endif

#if (defined PC_CLASSIFICATION)
     CALL deallocate_PCTimeVars
#endif

  END SUBROUTINE deallocate_TimeVariables


  !---------------------------------------
  function save_to_restart (idate, deltim, itstamp, ptstamp) result(rwrite)

     use mod_namelist
     implicit none

     logical :: rwrite

     integer,  intent(in) :: idate(3)
     real(r8), intent(in) :: deltim
     type(timestamp), intent(in) :: itstamp, ptstamp

     ! added by yuan, 08/31/2014
     select case (trim(DEF_WRST_FREQ))
     case ('HOURLY')
        rwrite = isendofhour (idate, deltim)
     case ('DAILY')
        rwrite = isendofday(idate, deltim)
     case ('MONTHLY')
        rwrite = isendofmonth(idate, deltim)       
     case ('YEARLY')
        rwrite = isendofyear(idate, deltim)
     end select

     if (rwrite) then
        rwrite = (ptstamp < itstamp)
     end if

  end function save_to_restart

  !---------------------------------------
  SUBROUTINE WRITE_TimeVariables (idate, site, dir_restart)

     !=======================================================================
     ! Original version: Yongjiu Dai, September 15, 1999, 03/2014
     !=======================================================================

     use mod_namelist, only : DEF_REST_COMPRESS_LEVEL 
     USE mod_landpatch
     use ncio_vector
     USE GlobalVars
#if (defined PFT_CLASSIFICATION)
     USE MOD_PFTimeVars
#endif
#if (defined PC_CLASSIFICATION)
     USE MOD_PCTimeVars
#endif
     IMPLICIT NONE

     integer, INTENT(in) :: idate(3)
     character(LEN=*), intent(in) :: site
     character(LEN=*), intent(in) :: dir_restart
     
     ! Local variables
     character(LEN=256) :: file_restart
     character(len=14)  :: cdate
     integer :: compress

     compress = DEF_REST_COMPRESS_LEVEL 

     write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') idate(1), idate(2), idate(3)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_'//trim(cdate)//'.nc'

     call ncio_create_file_vector (file_restart, landpatch)
     CALL ncio_define_pixelset_dimension (file_restart, landpatch)
     
     CALL ncio_define_dimension_vector (file_restart, 'snow',     -maxsnl       )
     CALL ncio_define_dimension_vector (file_restart, 'soilsnow', nl_soil-maxsnl)
     CALL ncio_define_dimension_vector (file_restart, 'soil',     nl_soil)
     CALL ncio_define_dimension_vector (file_restart, 'lake',     nl_lake)

#ifdef PLANT_HYDRAULIC_STRESS
     CALL ncio_define_dimension_vector (file_restart, 'vegnodes', nvegwcs)
#endif
     
     CALL ncio_define_dimension_vector (file_restart, 'band',   2)
     CALL ncio_define_dimension_vector (file_restart, 'wetdry', 2)
#ifdef BGC
     CALL ncio_define_dimension_vector (file_restart, 'ndecomp_pools', ndecomp_pools)
     CALL ncio_define_dimension_vector (file_restart, 'doy' , 365)
#endif

     ! Time-varying state variables which reaquired by restart run
     call ncio_write_vector (file_restart, 'z_sno   '   , 'snow', -maxsnl, 'vector', landpatch, z_sno , compress) !  node depth [m]
     call ncio_write_vector (file_restart, 'dz_sno  '   , 'snow', -maxsnl, 'vector', landpatch, dz_sno, compress) !  interface depth [m]
     call ncio_write_vector (file_restart, 't_soisno'   , 'soilsnow', nl_soil-maxsnl, 'vector', landpatch, t_soisno   , compress) !  soil temperature [K]
     call ncio_write_vector (file_restart, 'wliq_soisno', 'soilsnow', nl_soil-maxsnl, 'vector', landpatch, wliq_soisno, compress) !  liquid water in layers [kg/m2]
     call ncio_write_vector (file_restart, 'wice_soisno', 'soilsnow', nl_soil-maxsnl, 'vector', landpatch, wice_soisno, compress) !  ice lens in layers [kg/m2]
     call ncio_write_vector (file_restart, 'smp',         'soil', nl_soil, 'vector', landpatch, smp, compress) !  soil matrix potential [mm]
     call ncio_write_vector (file_restart, 'hk',          'soil', nl_soil, 'vector', landpatch, hk, compress) !  hydraulic conductivity [mm h2o/s]
#ifdef PLANT_HYDRAULIC_STRESS
     call ncio_write_vector (file_restart, 'vegwp',   'vegnodes', nvegwcs, 'vector', landpatch, vegwp, compress) !  vegetation water potential [mm]
     call ncio_write_vector (file_restart, 'gs0sun  ',    'vector', landpatch, gs0sun, compress) !  working copy of sunlit stomata conductance
     call ncio_write_vector (file_restart, 'gs0sha  ',    'vector', landpatch, gs0sha, compress) !  working copy of shalit stomata conductance
#endif
     call ncio_write_vector (file_restart, 't_grnd  '   , 'vector', landpatch, t_grnd    , compress) !  ground surface temperature [K]
     call ncio_write_vector (file_restart, 'tleaf   '   , 'vector', landpatch, tleaf     , compress) !  leaf temperature [K]
     call ncio_write_vector (file_restart, 'ldew    '   , 'vector', landpatch, ldew      , compress) !  depth of water on foliage [mm]
     call ncio_write_vector (file_restart, 'sag     '   , 'vector', landpatch, sag       , compress) !  non dimensional snow age [-]
     call ncio_write_vector (file_restart, 'scv     '   , 'vector', landpatch, scv       , compress) !  snow cover, water equivalent [mm]
     call ncio_write_vector (file_restart, 'snowdp  '   , 'vector', landpatch, snowdp    , compress) !  snow depth [meter]
     call ncio_write_vector (file_restart, 'fveg    '   , 'vector', landpatch, fveg      , compress) !  fraction of vegetation cover
     call ncio_write_vector (file_restart, 'fsno    '   , 'vector', landpatch, fsno      , compress) !  fraction of snow cover on ground
     call ncio_write_vector (file_restart, 'sigf    '   , 'vector', landpatch, sigf      , compress) !  fraction of veg cover, excluding snow-covered veg [-]
     call ncio_write_vector (file_restart, 'green   '   , 'vector', landpatch, green     , compress) !  leaf greenness
     call ncio_write_vector (file_restart, 'lai     '   , 'vector', landpatch, lai       , compress) !  leaf area index
     call ncio_write_vector (file_restart, 'tlai    '   , 'vector', landpatch, tlai      , compress) !  leaf area index
     call ncio_write_vector (file_restart, 'sai     '   , 'vector', landpatch, sai       , compress) !  stem area index
     call ncio_write_vector (file_restart, 'tsai    '   , 'vector', landpatch, tsai      , compress) !  stem area index
     call ncio_write_vector (file_restart, 'coszen  '   , 'vector', landpatch, coszen    , compress) !  cosine of solar zenith angle
     call ncio_write_vector (file_restart, 'alb     '   , 'band', 2, 'wetdry', 2, 'vector', landpatch, alb , compress) !  averaged albedo [-]
     call ncio_write_vector (file_restart, 'ssun    '   , 'band', 2, 'wetdry', 2, 'vector', landpatch, ssun, compress) !  sunlit canopy absorption for solar radiation (0-1)
     call ncio_write_vector (file_restart, 'ssha    '   , 'band', 2, 'wetdry', 2, 'vector', landpatch, ssha, compress) !  shaded canopy absorption for solar radiation (0-1)
     call ncio_write_vector (file_restart, 'thermk  '   , 'vector', landpatch, thermk    , compress) !  canopy gap fraction for tir radiation
     call ncio_write_vector (file_restart, 'extkb   '   , 'vector', landpatch, extkb     , compress) !  (k, g(mu)/mu) direct solar extinction coefficient
     call ncio_write_vector (file_restart, 'extkd   '   , 'vector', landpatch, extkd     , compress) !  diffuse and scattered diffuse PAR extinction coefficient
     call ncio_write_vector (file_restart, 'zwt     '   , 'vector', landpatch, zwt       , compress) !  the depth to water table [m]
     call ncio_write_vector (file_restart, 'wa      '   , 'vector', landpatch, wa        , compress) !  water storage in aquifer [mm]
#ifdef VARIABLY_SATURATED_FLOW
     call ncio_write_vector (file_restart, 'dpond   '   , 'vector', landpatch, dpond     , compress) ! depth of ponding water
#endif

     call ncio_write_vector (file_restart, 't_lake  '   , 'lake', nl_lake, 'vector', landpatch, t_lake      , compress) !
     call ncio_write_vector (file_restart, 'lake_icefrc', 'lake', nl_lake, 'vector', landpatch, lake_icefrac, compress) !
     call ncio_write_vector (file_restart, 'savedtke1  ', 'vector', landpatch, savedtke1   , compress) !

     ! Additional va_vectorriables required by reginal model (such as WRF ) RSM) 
     call ncio_write_vector (file_restart, 'trad ', 'vector', landpatch, trad , compress) !     radiative temperature of surface [K]
     call ncio_write_vector (file_restart, 'tref ', 'vector', landpatch, tref , compress) !     2 m height air temperature [kelvin]
     call ncio_write_vector (file_restart, 'qref ', 'vector', landpatch, qref , compress) !     2 m height air specific humidity
     call ncio_write_vector (file_restart, 'rst  ', 'vector', landpatch, rst  , compress) !     canopy stomatal resistance (s/m)
     call ncio_write_vector (file_restart, 'emis ', 'vector', landpatch, emis , compress) !     averaged bulk surface emissivity
     call ncio_write_vector (file_restart, 'z0m  ', 'vector', landpatch, z0m  , compress) !     effective roughness [m]
     call ncio_write_vector (file_restart, 'zol  ', 'vector', landpatch, zol  , compress) !     dimensionless height (z/L) used in Monin-Obukhov theory
     call ncio_write_vector (file_restart, 'rib  ', 'vector', landpatch, rib  , compress) !     bulk Richardson number in surface layer
     call ncio_write_vector (file_restart, 'ustar', 'vector', landpatch, ustar, compress) !     u* in similarity theory [m/s]
     call ncio_write_vector (file_restart, 'qstar', 'vector', landpatch, qstar, compress) !     q* in similarity theory [kg/kg]
     call ncio_write_vector (file_restart, 'tstar', 'vector', landpatch, tstar, compress) !     t* in similarity theory [K]
     call ncio_write_vector (file_restart, 'fm   ', 'vector', landpatch, fm   , compress) !     integral of profile function for momentum
     call ncio_write_vector (file_restart, 'fh   ', 'vector', landpatch, fh   , compress) !     integral of profile function for heat
     call ncio_write_vector (file_restart, 'fq   ', 'vector', landpatch, fq   , compress) !     integral of profile function for moisture

#ifdef BGC
! bgc variables
     call ncio_write_vector (file_restart, 'totlitc              ', 'vector', landpatch, totlitc              )
     call ncio_write_vector (file_restart, 'totvegc              ', 'vector', landpatch, totvegc              )
     call ncio_write_vector (file_restart, 'totsomc              ', 'vector', landpatch, totsomc              )
     call ncio_write_vector (file_restart, 'totcwdc              ', 'vector', landpatch, totcwdc              )
     call ncio_write_vector (file_restart, 'totcolc              ', 'vector', landpatch, totcolc              )
     call ncio_write_vector (file_restart, 'totlitn              ', 'vector', landpatch, totlitn              )
     call ncio_write_vector (file_restart, 'totvegn              ', 'vector', landpatch, totvegn              )
     call ncio_write_vector (file_restart, 'totsomn              ', 'vector', landpatch, totsomn              )
     call ncio_write_vector (file_restart, 'totcwdn              ', 'vector', landpatch, totcwdn              )
     call ncio_write_vector (file_restart, 'totcoln              ', 'vector', landpatch, totcoln              )
       
     call ncio_write_vector (file_restart, 'sminn                ', 'vector', landpatch, sminn                )

     call ncio_write_vector (file_restart, 'decomp_cpools_vr     ', 'soil'  ,   nl_soil, 'ndecomp_pools', ndecomp_pools, &
                                                                    'vector', landpatch,     decomp_cpools_vr(1:nl_soil,:,:))
     call ncio_write_vector (file_restart, 'ctrunc_vr            ', 'soil'  ,   nl_soil, 'vector', landpatch, ctrunc_vr(1:nl_soil,:))

     call ncio_write_vector (file_restart, 'altmax               ', 'vector', landpatch, altmax               )
     call ncio_write_vector (file_restart, 'altmax_lastyear      ', 'vector', landpatch, altmax_lastyear      )
     call ncio_write_vector (file_restart, 'altmax_lastyear_indx ', 'vector', landpatch, altmax_lastyear_indx )

     call ncio_write_vector (file_restart, 'decomp_npools_vr     ', 'soil'  ,   nl_soil, 'ndecomp_pools', ndecomp_pools, &
                                                                    'vector', landpatch,     decomp_npools_vr(1:nl_soil,:,:))
     call ncio_write_vector (file_restart, 'ntrunc_vr            ', 'soil'  ,   nl_soil, 'vector', landpatch, ntrunc_vr(1:nl_soil,:))
     call ncio_write_vector (file_restart, 'sminn_vr             ', 'soil'  ,   nl_soil, 'vector', landpatch, sminn_vr    )
     call ncio_write_vector (file_restart, 'smin_no3_vr          ', 'soil'  ,   nl_soil, 'vector', landpatch, smin_no3_vr )
     call ncio_write_vector (file_restart, 'smin_nh4_vr          ', 'soil'  ,   nl_soil, 'vector', landpatch, smin_nh4_vr )

     call ncio_write_vector (file_restart, 'prec10               ', 'vector', landpatch, prec10               )
     call ncio_write_vector (file_restart, 'prec60               ', 'vector', landpatch, prec60               )
     call ncio_write_vector (file_restart, 'prec365              ', 'vector', landpatch, prec365              )
     call ncio_write_vector (file_restart, 'prec_today           ', 'vector', landpatch, prec_today           )
     call ncio_write_vector (file_restart, 'prec_daily           ', 'doy'   ,       365, 'vector', landpatch, prec_daily  )
     call ncio_write_vector (file_restart, 'tsoi17               ', 'vector', landpatch, tsoi17               )
     call ncio_write_vector (file_restart, 'rh30                 ', 'vector', landpatch, rh30                 )
     call ncio_write_vector (file_restart, 'accumnstep           ', 'vector', landpatch, accumnstep           )
     call ncio_write_vector (file_restart, 'cphase               ', 'vector', landpatch, cphase               )

#ifdef SASU
!---------------SASU variables-----------------------
     call ncio_write_vector (file_restart, 'decomp0_cpools_vr            ', 'soil'  ,   nl_soil, 'ndecomp_pools', ndecomp_pools, &
                                                                            'vector', landpatch, decomp0_cpools_vr            )
     call ncio_write_vector (file_restart, 'I_met_c_vr_acc               ', 'soil'  ,   nl_soil, 'vector', landpatch, I_met_c_vr_acc               )
     call ncio_write_vector (file_restart, 'I_cel_c_vr_acc               ', 'soil'  ,   nl_soil, 'vector', landpatch, I_cel_c_vr_acc               )
     call ncio_write_vector (file_restart, 'I_lig_c_vr_acc               ', 'soil'  ,   nl_soil, 'vector', landpatch, I_lig_c_vr_acc               )
     call ncio_write_vector (file_restart, 'I_cwd_c_vr_acc               ', 'soil'  ,   nl_soil, 'vector', landpatch, I_cwd_c_vr_acc               )
     call ncio_write_vector (file_restart, 'AKX_met_to_soil1_c_vr_acc    ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_met_to_soil1_c_vr_acc    )
     call ncio_write_vector (file_restart, 'AKX_cel_to_soil1_c_vr_acc    ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_cel_to_soil1_c_vr_acc    )
     call ncio_write_vector (file_restart, 'AKX_lig_to_soil2_c_vr_acc    ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_lig_to_soil2_c_vr_acc    )
     call ncio_write_vector (file_restart, 'AKX_soil1_to_soil2_c_vr_acc  ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_soil1_to_soil2_c_vr_acc  )
     call ncio_write_vector (file_restart, 'AKX_cwd_to_cel_c_vr_acc      ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_cwd_to_cel_c_vr_acc      )
     call ncio_write_vector (file_restart, 'AKX_cwd_to_lig_c_vr_acc      ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_cwd_to_lig_c_vr_acc      )
     call ncio_write_vector (file_restart, 'AKX_soil1_to_soil3_c_vr_acc  ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_soil1_to_soil3_c_vr_acc  )
     call ncio_write_vector (file_restart, 'AKX_soil2_to_soil1_c_vr_acc  ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_soil2_to_soil1_c_vr_acc  )
     call ncio_write_vector (file_restart, 'AKX_soil2_to_soil3_c_vr_acc  ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_soil2_to_soil3_c_vr_acc  )
     call ncio_write_vector (file_restart, 'AKX_soil3_to_soil1_c_vr_acc  ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_soil3_to_soil1_c_vr_acc  )
     call ncio_write_vector (file_restart, 'AKX_met_exit_c_vr_acc        ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_met_exit_c_vr_acc        )
     call ncio_write_vector (file_restart, 'AKX_cel_exit_c_vr_acc        ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_cel_exit_c_vr_acc        )
     call ncio_write_vector (file_restart, 'AKX_lig_exit_c_vr_acc        ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_lig_exit_c_vr_acc        )
     call ncio_write_vector (file_restart, 'AKX_cwd_exit_c_vr_acc        ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_cwd_exit_c_vr_acc        )
     call ncio_write_vector (file_restart, 'AKX_soil1_exit_c_vr_acc      ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_soil1_exit_c_vr_acc      )
     call ncio_write_vector (file_restart, 'AKX_soil2_exit_c_vr_acc      ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_soil2_exit_c_vr_acc      )
     call ncio_write_vector (file_restart, 'AKX_soil3_exit_c_vr_acc      ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_soil3_exit_c_vr_acc      )

     call ncio_write_vector (file_restart, 'decomp0_npools_vr            ', 'soil'  ,   nl_soil, 'ndecomp_pools', ndecomp_pools, &
                                                                            'vector', landpatch, decomp0_npools_vr            )
     call ncio_write_vector (file_restart, 'I_met_n_vr_acc               ', 'soil'  ,   nl_soil, 'vector', landpatch, I_met_n_vr_acc               )
     call ncio_write_vector (file_restart, 'I_cel_n_vr_acc               ', 'vector', landpatch, I_cel_n_vr_acc               )
     call ncio_write_vector (file_restart, 'I_lig_n_vr_acc               ', 'soil'  ,   nl_soil, 'vector', landpatch, I_lig_n_vr_acc               )
     call ncio_write_vector (file_restart, 'I_cwd_n_vr_acc               ', 'soil'  ,   nl_soil, 'vector', landpatch, I_cwd_n_vr_acc               )
     call ncio_write_vector (file_restart, 'AKX_met_to_soil1_n_vr_acc    ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_met_to_soil1_n_vr_acc    )
     call ncio_write_vector (file_restart, 'AKX_cel_to_soil1_n_vr_acc    ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_cel_to_soil1_n_vr_acc    )
     call ncio_write_vector (file_restart, 'AKX_lig_to_soil2_n_vr_acc    ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_lig_to_soil2_n_vr_acc    )
     call ncio_write_vector (file_restart, 'AKX_soil1_to_soil2_n_vr_acc  ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_soil1_to_soil2_n_vr_acc  )
     call ncio_write_vector (file_restart, 'AKX_cwd_to_cel_n_vr_acc      ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_cwd_to_cel_n_vr_acc      )
     call ncio_write_vector (file_restart, 'AKX_cwd_to_lig_n_vr_acc      ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_cwd_to_lig_n_vr_acc      )
     call ncio_write_vector (file_restart, 'AKX_soil1_to_soil3_n_vr_acc  ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_soil1_to_soil3_n_vr_acc  )
     call ncio_write_vector (file_restart, 'AKX_soil2_to_soil1_n_vr_acc  ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_soil2_to_soil1_n_vr_acc  )
     call ncio_write_vector (file_restart, 'AKX_soil2_to_soil3_n_vr_acc  ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_soil2_to_soil3_n_vr_acc  )
     call ncio_write_vector (file_restart, 'AKX_soil3_to_soil1_n_vr_acc  ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_soil3_to_soil1_n_vr_acc  )
     call ncio_write_vector (file_restart, 'AKX_met_exit_n_vr_acc        ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_met_exit_n_vr_acc        )
     call ncio_write_vector (file_restart, 'AKX_cel_exit_n_vr_acc        ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_cel_exit_n_vr_acc        )
     call ncio_write_vector (file_restart, 'AKX_lig_exit_n_vr_acc        ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_lig_exit_n_vr_acc        )
     call ncio_write_vector (file_restart, 'AKX_cwd_exit_n_vr_acc        ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_cwd_exit_n_vr_acc        )
     call ncio_write_vector (file_restart, 'AKX_soil1_exit_n_vr_acc      ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_soil1_exit_n_vr_acc      )
     call ncio_write_vector (file_restart, 'AKX_soil2_exit_n_vr_acc      ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_soil2_exit_n_vr_acc      )
     call ncio_write_vector (file_restart, 'AKX_soil3_exit_n_vr_acc      ', 'soil'  ,   nl_soil, 'vector', landpatch, AKX_soil3_exit_n_vr_acc      )

     call ncio_write_vector (file_restart, 'diagVX_c_vr_acc              ', 'soil'  ,   nl_soil, 'ndecomp_pools', ndecomp_pools, &
                                                                            'vector', landpatch, diagVX_c_vr_acc              )
     call ncio_write_vector (file_restart, 'upperVX_c_vr_acc             ', 'soil'  ,   nl_soil, 'ndecomp_pools', ndecomp_pools, &
                                                                            'vector', landpatch, upperVX_c_vr_acc             )
     call ncio_write_vector (file_restart, 'lowerVX_c_vr_acc             ', 'soil'  ,   nl_soil, 'ndecomp_pools', ndecomp_pools, &
                                                                            'vector', landpatch, lowerVX_c_vr_acc             )
     call ncio_write_vector (file_restart, 'diagVX_n_vr_acc              ', 'soil'  ,   nl_soil, 'ndecomp_pools', ndecomp_pools, &
                                                                            'vector', landpatch, diagVX_n_vr_acc              )
     call ncio_write_vector (file_restart, 'upperVX_n_vr_acc             ', 'soil'  ,   nl_soil, 'ndecomp_pools', ndecomp_pools, &
                                                                            'vector', landpatch, upperVX_n_vr_acc             )
     call ncio_write_vector (file_restart, 'lowerVX_n_vr_acc             ', 'soil'  ,   nl_soil, 'ndecomp_pools', ndecomp_pools, &
                                                                            'vector', landpatch, lowerVX_n_vr_acc             )

!----------------------------------------------------
#endif
     call ncio_write_vector (file_restart, 'skip_balance_check           ', 'vector', landpatch, skip_balance_check           )
#endif

#if (defined PFT_CLASSIFICATION)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_pft_'//trim(cdate)//'.nc'
     CALL WRITE_PFTimeVars (file_restart)
#endif

#if (defined PC_CLASSIFICATION)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_pc_'//trim(cdate)//'.nc'
     CALL WRITE_PCTimeVars (file_restart)
#endif

  end subroutine WRITE_TimeVariables

  !---------------------------------------
  SUBROUTINE READ_TimeVariables (idate, site, dir_restart)

     !=======================================================================
     ! Original version: Yongjiu Dai, September 15, 1999, 03/2014
     !=======================================================================

     use mod_namelist
     use spmd_task
     use ncio_vector
#ifdef CLMDEBUG 
   USE mod_colm_debug
#endif
     USE mod_landpatch
     USE GlobalVars
#ifdef PFT_CLASSIFICATION
     USE MOD_PFTimeVars
#endif
#ifdef PC_CLASSIFICATION
     USE MOD_PCTimeVars
#endif

     IMPLICIT NONE

     integer, INTENT(in) :: idate(3)
     character(LEN=*), intent(in) :: site
     character(LEN=*), intent(in) :: dir_restart

     ! Local variables
     character(LEN=256) :: file_restart
     character(len=14)  :: cdate
     
#ifdef USEMPI
     call mpi_barrier (p_comm_glb, p_err)
#endif

     if (p_is_master) then
        write(*,'(/,A26)') 'Loading Time Variables ...'
     end if

     write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') idate(1), idate(2), idate(3)
     file_restart = trim(dir_restart) // '/' // trim(site) //'_restart_'//trim(cdate)//'.nc'

     ! Time-varying state variables which reaquired by restart run
     call ncio_read_vector (file_restart, 'z_sno   '   , -maxsnl, landpatch, z_sno ) !  node depth [m]
     call ncio_read_vector (file_restart, 'dz_sno  '   , -maxsnl, landpatch, dz_sno) !  interface depth [m]
     call ncio_read_vector (file_restart, 't_soisno'   , nl_soil-maxsnl, landpatch, t_soisno   ) !  soil temperature [K]
     call ncio_read_vector (file_restart, 'wliq_soisno', nl_soil-maxsnl, landpatch, wliq_soisno) !  liquid water in layers [kg/m2]
     call ncio_read_vector (file_restart, 'wice_soisno', nl_soil-maxsnl, landpatch, wice_soisno) !  ice lens in layers [kg/m2]
     call ncio_read_vector (file_restart, 'smp',         nl_soil,        landpatch, smp        ) !  soil matrix potential [mm]
     call ncio_read_vector (file_restart, 'hk',          nl_soil,        landpatch, hk         ) !  hydraulic conductivity [mm h2o/s]
#ifdef PLANT_HYDRAULIC_STRESS
     call ncio_read_vector (file_restart, 'vegwp',       nvegwcs,        landpatch, vegwp      ) !  vegetation water potential [mm]
     call ncio_read_vector (file_restart, 'gs0sun  ',    landpatch, gs0sun     ) !  working copy of sunlit stomata conductance
     call ncio_read_vector (file_restart, 'gs0sha  ',    landpatch, gs0sha     ) !  working copy of shalit stomata conductance
#endif
     call ncio_read_vector (file_restart, 't_grnd  '   , landpatch, t_grnd     ) !  ground surface temperature [K]
     call ncio_read_vector (file_restart, 'tleaf   '   , landpatch, tleaf      ) !  leaf temperature [K]
     call ncio_read_vector (file_restart, 'ldew    '   , landpatch, ldew       ) !  depth of water on foliage [mm]
     call ncio_read_vector (file_restart, 'sag     '   , landpatch, sag        ) !  non dimensional snow age [-]
     call ncio_read_vector (file_restart, 'scv     '   , landpatch, scv        ) !  snow cover, water equivalent [mm]
     call ncio_read_vector (file_restart, 'snowdp  '   , landpatch, snowdp     ) !  snow depth [meter]
     call ncio_read_vector (file_restart, 'fveg    '   , landpatch, fveg       ) !  fraction of vegetation cover
     call ncio_read_vector (file_restart, 'fsno    '   , landpatch, fsno       ) !  fraction of snow cover on ground
     call ncio_read_vector (file_restart, 'sigf    '   , landpatch, sigf       ) !  fraction of veg cover, excluding snow-covered veg [-]
     call ncio_read_vector (file_restart, 'green   '   , landpatch, green      ) !  leaf greenness
     call ncio_read_vector (file_restart, 'lai     '   , landpatch, lai        ) !  leaf area index
     call ncio_read_vector (file_restart, 'tlai    '   , landpatch, tlai       ) !  leaf area index
     call ncio_read_vector (file_restart, 'sai     '   , landpatch, sai        ) !  stem area index
     call ncio_read_vector (file_restart, 'tsai    '   , landpatch, tsai       ) !  stem area index
     call ncio_read_vector (file_restart, 'coszen  '   , landpatch, coszen     ) !  cosine of solar zenith angle
     call ncio_read_vector (file_restart, 'alb     '   , 2, 2, landpatch, alb  ) !  averaged albedo [-]
     call ncio_read_vector (file_restart, 'ssun    '   , 2, 2, landpatch, ssun ) !  sunlit canopy absorption for solar radiation (0-1)
     call ncio_read_vector (file_restart, 'ssha    '   , 2, 2, landpatch, ssha ) !  shaded canopy absorption for solar radiation (0-1)
     call ncio_read_vector (file_restart, 'thermk  '   , landpatch, thermk     ) !  canopy gap fraction for tir radiation
     call ncio_read_vector (file_restart, 'extkb   '   , landpatch, extkb      ) !  (k, g(mu)/mu) direct solar extinction coefficient
     call ncio_read_vector (file_restart, 'extkd   '   , landpatch, extkd      ) !  diffuse and scattered diffuse PAR extinction coefficient
     call ncio_read_vector (file_restart, 'zwt     '   , landpatch, zwt        ) !  the depth to water table [m]
     call ncio_read_vector (file_restart, 'wa      '   , landpatch, wa         ) !  water storage in aquifer [mm]
#ifdef VARIABLY_SATURATED_FLOW
     call ncio_read_vector (file_restart, 'dpond   '   , landpatch, dpond      ) ! depth of ponding water
#endif

     call ncio_read_vector (file_restart, 't_lake  '   , nl_lake, landpatch, t_lake      ) !
     call ncio_read_vector (file_restart, 'lake_icefrc', nl_lake, landpatch, lake_icefrac) !
     call ncio_read_vector (file_restart, 'savedtke1', landpatch, savedtke1) !

     ! Additional variables required by reginal model (such as WRF ) RSM) 
     call ncio_read_vector (file_restart, 'trad ', landpatch, trad ) !     radiative temperature of surface [K]
     call ncio_read_vector (file_restart, 'tref ', landpatch, tref ) !     2 m height air temperature [kelvin]
     call ncio_read_vector (file_restart, 'qref ', landpatch, qref ) !     2 m height air specific humidity
     call ncio_read_vector (file_restart, 'rst  ', landpatch, rst  ) !     canopy stomatal resistance (s/m)
     call ncio_read_vector (file_restart, 'emis ', landpatch, emis ) !     averaged bulk surface emissivity
     call ncio_read_vector (file_restart, 'z0m  ', landpatch, z0m  ) !     effective roughness [m]
     call ncio_read_vector (file_restart, 'zol  ', landpatch, zol  ) !     dimensionless height (z/L) used in Monin-Obukhov theory
     call ncio_read_vector (file_restart, 'rib  ', landpatch, rib  ) !     bulk Richardson number in surface layer
     call ncio_read_vector (file_restart, 'ustar', landpatch, ustar) !     u* in similarity theory [m/s]
     call ncio_read_vector (file_restart, 'qstar', landpatch, qstar) !     q* in similarity theory [kg/kg]
     call ncio_read_vector (file_restart, 'tstar', landpatch, tstar) !     t* in similarity theory [K]
     call ncio_read_vector (file_restart, 'fm   ', landpatch, fm   ) !     integral of profile function for momentum
     call ncio_read_vector (file_restart, 'fh   ', landpatch, fh   ) !     integral of profile function for heat
     call ncio_read_vector (file_restart, 'fq   ', landpatch, fq   ) !     integral of profile function for moisture


#ifdef BGC
! bgc variables
     call ncio_read_vector (file_restart, 'totlitc              ', landpatch, totlitc              )
     call ncio_read_vector (file_restart, 'totvegc              ', landpatch, totvegc              )
     call ncio_read_vector (file_restart, 'totsomc              ', landpatch, totsomc              )
     call ncio_read_vector (file_restart, 'totcwdc              ', landpatch, totcwdc              )
     call ncio_read_vector (file_restart, 'totcolc              ', landpatch, totcolc              )
     call ncio_read_vector (file_restart, 'totlitn              ', landpatch, totlitn              )
     call ncio_read_vector (file_restart, 'totvegn              ', landpatch, totvegn              )
     call ncio_read_vector (file_restart, 'totsomn              ', landpatch, totsomn              )
     call ncio_read_vector (file_restart, 'totcwdn              ', landpatch, totcwdn              )
     call ncio_read_vector (file_restart, 'totcoln              ', landpatch, totcoln              )
       
     call ncio_read_vector (file_restart, 'sminn                ', landpatch, sminn                )

     call ncio_read_vector (file_restart, 'decomp_cpools_vr     ',   nl_soil, ndecomp_pools, landpatch, decomp_cpools_vr)
     call ncio_read_vector (file_restart, 'ctrunc_vr            ',   nl_soil, landpatch, ctrunc_vr            )

     call ncio_read_vector (file_restart, 'altmax               ', landpatch, altmax               )
     call ncio_read_vector (file_restart, 'altmax_lastyear      ', landpatch, altmax_lastyear      )
     call ncio_read_vector (file_restart, 'altmax_lastyear_indx ', landpatch, altmax_lastyear_indx )

     call ncio_read_vector (file_restart, 'decomp_npools_vr     ',   nl_soil, ndecomp_pools, landpatch, decomp_npools_vr)
     call ncio_read_vector (file_restart, 'ntrunc_vr            ',   nl_soil, landpatch, ntrunc_vr            )
     call ncio_read_vector (file_restart, 'sminn_vr             ',   nl_soil, landpatch, sminn_vr             )
     call ncio_read_vector (file_restart, 'smin_no3_vr          ',   nl_soil, landpatch, smin_no3_vr          )
     call ncio_read_vector (file_restart, 'smin_nh4_vr          ',   nl_soil, landpatch, smin_nh4_vr          )

     call ncio_read_vector (file_restart, 'prec10               ', landpatch, prec10               )
     call ncio_read_vector (file_restart, 'prec60               ', landpatch, prec60               )
     call ncio_read_vector (file_restart, 'prec365              ', landpatch, prec365              )
     call ncio_read_vector (file_restart, 'prec_today           ', landpatch, prec_today           )
     call ncio_read_vector (file_restart, 'prec_daily           ',       365, landpatch, prec_daily)
     call ncio_read_vector (file_restart, 'tsoi17               ', landpatch, tsoi17               )
     call ncio_read_vector (file_restart, 'rh30                 ', landpatch, rh30                 )
     call ncio_read_vector (file_restart, 'accumnstep           ', landpatch, accumnstep           )
     call ncio_read_vector (file_restart, 'cphase               ', landpatch, cphase               )

#ifdef SASU
!---------------SASU variables-----------------------
     call ncio_read_vector (file_restart, 'decomp0_cpools_vr            ',   nl_soil, ndecomp_pools, landpatch, decomp0_cpools_vr            )
     call ncio_read_vector (file_restart, 'I_met_c_vr_acc               ',   nl_soil, landpatch, I_met_c_vr_acc               )
     call ncio_read_vector (file_restart, 'I_cel_c_vr_acc               ',   nl_soil, landpatch, I_cel_c_vr_acc               )
     call ncio_read_vector (file_restart, 'I_lig_c_vr_acc               ',   nl_soil, landpatch, I_lig_c_vr_acc               )
     call ncio_read_vector (file_restart, 'I_cwd_c_vr_acc               ',   nl_soil, landpatch, I_cwd_c_vr_acc               )
     call ncio_read_vector (file_restart, 'AKX_met_to_soil1_c_vr_acc    ',   nl_soil, landpatch, AKX_met_to_soil1_c_vr_acc    )
     call ncio_read_vector (file_restart, 'AKX_cel_to_soil1_c_vr_acc    ',   nl_soil, landpatch, AKX_cel_to_soil1_c_vr_acc    )
     call ncio_read_vector (file_restart, 'AKX_lig_to_soil2_c_vr_acc    ',   nl_soil, landpatch, AKX_lig_to_soil2_c_vr_acc    )
     call ncio_read_vector (file_restart, 'AKX_soil1_to_soil2_c_vr_acc  ',   nl_soil, landpatch, AKX_soil1_to_soil2_c_vr_acc  )
     call ncio_read_vector (file_restart, 'AKX_cwd_to_cel_c_vr_acc      ',   nl_soil, landpatch, AKX_cwd_to_cel_c_vr_acc      )
     call ncio_read_vector (file_restart, 'AKX_cwd_to_lig_c_vr_acc      ',   nl_soil, landpatch, AKX_cwd_to_lig_c_vr_acc      )
     call ncio_read_vector (file_restart, 'AKX_soil1_to_soil3_c_vr_acc  ',   nl_soil, landpatch, AKX_soil1_to_soil3_c_vr_acc  )
     call ncio_read_vector (file_restart, 'AKX_soil2_to_soil1_c_vr_acc  ',   nl_soil, landpatch, AKX_soil2_to_soil1_c_vr_acc  )
     call ncio_read_vector (file_restart, 'AKX_soil2_to_soil3_c_vr_acc  ',   nl_soil, landpatch, AKX_soil2_to_soil3_c_vr_acc  )
     call ncio_read_vector (file_restart, 'AKX_soil3_to_soil1_c_vr_acc  ',   nl_soil, landpatch, AKX_soil3_to_soil1_c_vr_acc  )
     call ncio_read_vector (file_restart, 'AKX_met_exit_c_vr_acc        ',   nl_soil, landpatch, AKX_met_exit_c_vr_acc        )
     call ncio_read_vector (file_restart, 'AKX_cel_exit_c_vr_acc        ',   nl_soil, landpatch, AKX_cel_exit_c_vr_acc        )
     call ncio_read_vector (file_restart, 'AKX_lig_exit_c_vr_acc        ',   nl_soil, landpatch, AKX_lig_exit_c_vr_acc        )
     call ncio_read_vector (file_restart, 'AKX_cwd_exit_c_vr_acc        ',   nl_soil, landpatch, AKX_cwd_exit_c_vr_acc        )
     call ncio_read_vector (file_restart, 'AKX_soil1_exit_c_vr_acc      ',   nl_soil, landpatch, AKX_soil1_exit_c_vr_acc      )
     call ncio_read_vector (file_restart, 'AKX_soil2_exit_c_vr_acc      ',   nl_soil, landpatch, AKX_soil2_exit_c_vr_acc      )
     call ncio_read_vector (file_restart, 'AKX_soil3_exit_c_vr_acc      ',   nl_soil, landpatch, AKX_soil3_exit_c_vr_acc      )

     call ncio_read_vector (file_restart, 'decomp0_npools_vr            ',   nl_soil, ndecomp_pools, landpatch, decomp0_npools_vr            )
     call ncio_read_vector (file_restart, 'I_met_n_vr_acc               ',   nl_soil, landpatch, I_met_n_vr_acc               )
     call ncio_read_vector (file_restart, 'I_cel_n_vr_acc               ',   nl_soil, landpatch, I_cel_n_vr_acc               )
     call ncio_read_vector (file_restart, 'I_lig_n_vr_acc               ',   nl_soil, landpatch, I_lig_n_vr_acc               )
     call ncio_read_vector (file_restart, 'I_cwd_n_vr_acc               ',   nl_soil, landpatch, I_cwd_n_vr_acc               )
     call ncio_read_vector (file_restart, 'AKX_met_to_soil1_n_vr_acc    ',   nl_soil, landpatch, AKX_met_to_soil1_n_vr_acc    )
     call ncio_read_vector (file_restart, 'AKX_cel_to_soil1_n_vr_acc    ',   nl_soil, landpatch, AKX_cel_to_soil1_n_vr_acc    )
     call ncio_read_vector (file_restart, 'AKX_lig_to_soil2_n_vr_acc    ',   nl_soil, landpatch, AKX_lig_to_soil2_n_vr_acc    )
     call ncio_read_vector (file_restart, 'AKX_soil1_to_soil2_n_vr_acc  ',   nl_soil, landpatch, AKX_soil1_to_soil2_n_vr_acc  )
     call ncio_read_vector (file_restart, 'AKX_cwd_to_cel_n_vr_acc      ',   nl_soil, landpatch, AKX_cwd_to_cel_n_vr_acc      )
     call ncio_read_vector (file_restart, 'AKX_cwd_to_lig_n_vr_acc      ',   nl_soil, landpatch, AKX_cwd_to_lig_n_vr_acc      )
     call ncio_read_vector (file_restart, 'AKX_soil1_to_soil3_n_vr_acc  ',   nl_soil, landpatch, AKX_soil1_to_soil3_n_vr_acc  )
     call ncio_read_vector (file_restart, 'AKX_soil2_to_soil1_n_vr_acc  ',   nl_soil, landpatch, AKX_soil2_to_soil1_n_vr_acc  )
     call ncio_read_vector (file_restart, 'AKX_soil2_to_soil3_n_vr_acc  ',   nl_soil, landpatch, AKX_soil2_to_soil3_n_vr_acc  )
     call ncio_read_vector (file_restart, 'AKX_soil3_to_soil1_n_vr_acc  ',   nl_soil, landpatch, AKX_soil3_to_soil1_n_vr_acc  )
     call ncio_read_vector (file_restart, 'AKX_met_exit_n_vr_acc        ',   nl_soil, landpatch, AKX_met_exit_n_vr_acc        )
     call ncio_read_vector (file_restart, 'AKX_cel_exit_n_vr_acc        ',   nl_soil, landpatch, AKX_cel_exit_n_vr_acc        )
     call ncio_read_vector (file_restart, 'AKX_lig_exit_n_vr_acc        ',   nl_soil, landpatch, AKX_lig_exit_n_vr_acc        )
     call ncio_read_vector (file_restart, 'AKX_cwd_exit_n_vr_acc        ',   nl_soil, landpatch, AKX_cwd_exit_n_vr_acc        )
     call ncio_read_vector (file_restart, 'AKX_soil1_exit_n_vr_acc      ',   nl_soil, landpatch, AKX_soil1_exit_n_vr_acc      )
     call ncio_read_vector (file_restart, 'AKX_soil2_exit_n_vr_acc      ',   nl_soil, landpatch, AKX_soil2_exit_n_vr_acc      )
     call ncio_read_vector (file_restart, 'AKX_soil3_exit_n_vr_acc      ',   nl_soil, landpatch, AKX_soil3_exit_n_vr_acc      )

     call ncio_read_vector (file_restart, 'diagVX_c_vr_acc              ',   nl_soil, ndecomp_pools, landpatch, diagVX_c_vr_acc              )
     call ncio_read_vector (file_restart, 'upperVX_c_vr_acc             ',   nl_soil, ndecomp_pools, landpatch, upperVX_c_vr_acc             )
     call ncio_read_vector (file_restart, 'lowerVX_c_vr_acc             ',   nl_soil, ndecomp_pools, landpatch, lowerVX_c_vr_acc             )
     call ncio_read_vector (file_restart, 'diagVX_n_vr_acc              ',   nl_soil, ndecomp_pools, landpatch, diagVX_n_vr_acc              )
     call ncio_read_vector (file_restart, 'upperVX_n_vr_acc             ',   nl_soil, ndecomp_pools, landpatch, upperVX_n_vr_acc             )
     call ncio_read_vector (file_restart, 'lowerVX_n_vr_acc             ',   nl_soil, ndecomp_pools, landpatch, lowerVX_n_vr_acc             )

!----------------------------------------------------
#endif
     call ncio_read_vector (file_restart, 'skip_balance_check           ', landpatch, skip_balance_check           )
#endif

#if (defined PFT_CLASSIFICATION)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_pft_'//trim(cdate)//'.nc'
     CALL READ_PFTimeVars (file_restart)
#endif

#if (defined PC_CLASSIFICATION)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_pc_'//trim(cdate)//'.nc'
     CALL READ_PCTimeVars (file_restart)
#endif
     
#ifdef CLMDEBUG
     call check_TimeVariables
#endif
     
     if (p_is_master) then
        write(*,*) 'Loading Time Variables done.'
     end if

  end subroutine READ_TimeVariables

  !---------------------------------------
#ifdef CLMDEBUG
  SUBROUTINE check_TimeVariables ()

     use spmd_task
     use mod_colm_debug
#ifdef PFT_CLASSIFICATION
     USE MOD_PFTimeVars
#endif
#ifdef PC_CLASSIFICATION
     USE MOD_PCTimeVars
#endif

     IMPLICIT NONE

#ifdef USEMPI
     call mpi_barrier (p_comm_glb, p_err)
#endif
     if (p_is_master) then
        write(*,'(/,A27)') 'Checking Time Variables ...'
     end if

     call check_vector_data ('z_sno       ', z_sno )      !  node depth [m]
     call check_vector_data ('dz_sno      ', dz_sno)      !  interface depth [m]
     call check_vector_data ('t_soisno    ', t_soisno   ) !  soil temperature [K]
     call check_vector_data ('wliq_soisno ', wliq_soisno) !  liquid water in layers [kg/m2]
     call check_vector_data ('wice_soisno ', wice_soisno) !  ice lens in layers [kg/m2]
     call check_vector_data ('smp         ', smp        ) !  soil matrix potential [mm]
     call check_vector_data ('hk          ', hk         ) !  hydraulic conductivity [mm h2o/s]
#ifdef PLANT_HYDRAULIC_STRESS
     call check_vector_data ('vegwp       ', vegwp      ) !  vegetation water potential [mm]
     call check_vector_data ('gs0sun      ', gs0sun     ) !  working copy of sunlit stomata conductance
     call check_vector_data ('gs0sha      ', gs0sha     ) !  working copy of shalit stomata conductance
#endif
     call check_vector_data ('t_grnd      ', t_grnd     ) !  ground surface temperature [K]
     call check_vector_data ('tleaf       ', tleaf      ) !  leaf temperature [K]
     call check_vector_data ('ldew        ', ldew       ) !  depth of water on foliage [mm]
     call check_vector_data ('sag         ', sag        ) !  non dimensional snow age [-]
     call check_vector_data ('scv         ', scv        ) !  snow cover, water equivalent [mm]
     call check_vector_data ('snowdp      ', snowdp     ) !  snow depth [meter]
     call check_vector_data ('fveg        ', fveg       ) !  fraction of vegetation cover
     call check_vector_data ('fsno        ', fsno       ) !  fraction of snow cover on ground
     call check_vector_data ('sigf        ', sigf       ) !  fraction of veg cover, excluding snow-covered veg [-]
     call check_vector_data ('green       ', green      ) !  leaf greenness
     call check_vector_data ('lai         ', lai        ) !  leaf area index
     call check_vector_data ('tlai        ', tlai       ) !  leaf area index
     call check_vector_data ('sai         ', sai        ) !  stem area index
     call check_vector_data ('tsai        ', tsai       ) !  stem area index
     call check_vector_data ('coszen      ', coszen     ) !  cosine of solar zenith angle
     call check_vector_data ('alb         ', alb  ) !  averaged albedo [-]
     call check_vector_data ('ssun        ', ssun ) !  sunlit canopy absorption for solar radiation (0-1)
     call check_vector_data ('ssha        ', ssha ) !  shaded canopy absorption for solar radiation (0-1)
     call check_vector_data ('thermk      ', thermk     ) !  canopy gap fraction for tir radiation
     call check_vector_data ('extkb       ', extkb      ) !  (k, g(mu)/mu) direct solar extinction coefficient
     call check_vector_data ('extkd       ', extkd      ) !  diffuse and scattered diffuse PAR extinction coefficient
     call check_vector_data ('zwt         ', zwt        ) !  the depth to water table [m]
     call check_vector_data ('wa          ', wa         ) !  water storage in aquifer [mm]
#ifdef VARIABLY_SATURATED_FLOW
     call check_vector_data ('dpond       ', dpond      ) !  depth of ponding water
#endif
     call check_vector_data ('t_lake      ', t_lake      ) !
     call check_vector_data ('lake_icefrc ', lake_icefrac) !
     call check_vector_data ('savedtke1   ', savedtke1   ) !

#ifdef BGC
! bgc variables
     call check_vector_data ('decomp_cpools_vr  ', decomp_cpools_vr  ) 
     call check_vector_data ('decomp_cpools     ', decomp_cpools     ) 
     call check_vector_data ('decomp_k          ', decomp_k          ) 
     call check_vector_data ('ctrunc_vr         ', ctrunc_vr         ) 
     call check_vector_data ('ctrunc_veg        ', ctrunc_veg        ) 
     call check_vector_data ('ctrunc_soil       ', ctrunc_soil       ) 

     call check_vector_data ('t_scalar          ', t_scalar          ) 
     call check_vector_data ('w_scalar          ', w_scalar          ) 
     call check_vector_data ('o_scalar          ', o_scalar          ) 
     call check_vector_data ('depth_scalar      ', depth_scalar      ) 

! Soil CN diffusion and advection
     call check_vector_data ('som_adv_coef             ', som_adv_coef             ) 
     call check_vector_data ('som_diffus_coef          ', som_diffus_coef          ) 

! Active Layer
     call check_vector_data ('altmax                   ', altmax                   ) 
     call check_vector_data ('altmax_lastyear          ', altmax_lastyear          ) 
     !call check_vector_data ('altmax_lastyear_indx     ', altmax_lastyear_indx     ) 

     call check_vector_data ('totlitc                  ', totlitc                  ) 
     call check_vector_data ('totvegc                  ', totvegc                  ) 
     call check_vector_data ('totsomc                  ', totsomc                  ) 
     call check_vector_data ('totcwdc                  ', totcwdc                  ) 
     call check_vector_data ('totcolc                  ', totcolc                  ) 
     call check_vector_data ('col_begcb                ', col_begcb                ) 
     call check_vector_data ('col_endcb                ', col_endcb                ) 
     call check_vector_data ('col_vegbegcb             ', col_vegbegcb             ) 
     call check_vector_data ('col_vegendcb             ', col_vegendcb             ) 
     call check_vector_data ('col_soilbegcb            ', col_soilbegcb            ) 
     call check_vector_data ('col_soilendcb            ', col_soilendcb            ) 

     call check_vector_data ('totlitn                  ', totlitn                  ) 
     call check_vector_data ('totvegn                  ', totvegn                  ) 
     call check_vector_data ('totsomn                  ', totsomn                  ) 
     call check_vector_data ('totcwdn                  ', totcwdn                  ) 
     call check_vector_data ('totcoln                  ', totcoln                  ) 
     call check_vector_data ('col_begnb                ', col_begnb                ) 
     call check_vector_data ('col_endnb                ', col_endnb                ) 
     call check_vector_data ('col_vegbegnb             ', col_vegbegnb             ) 
     call check_vector_data ('col_vegendnb             ', col_vegendnb             ) 
     call check_vector_data ('col_soilbegnb            ', col_soilbegnb            ) 
     call check_vector_data ('col_soilendnb            ', col_soilendnb            ) 
     call check_vector_data ('col_sminnbegnb           ', col_sminnbegnb           ) 
     call check_vector_data ('col_sminnendnb           ', col_sminnendnb           ) 

     call check_vector_data ('leafc                    ', leafc                    ) 
     call check_vector_data ('leafc_storage            ', leafc_storage            ) 
     call check_vector_data ('leafc_xfer               ', leafc_xfer               ) 
     call check_vector_data ('frootc                   ', frootc                   ) 
     call check_vector_data ('frootc_storage           ', frootc_storage           ) 
     call check_vector_data ('frootc_xfer              ', frootc_xfer              ) 
     call check_vector_data ('livestemc                ', livestemc                ) 
     call check_vector_data ('livestemc_storage        ', livestemc_storage        ) 
     call check_vector_data ('livestemc_xfer           ', livestemc_xfer           ) 
     call check_vector_data ('deadstemc                ', deadstemc                ) 
     call check_vector_data ('deadstemc_storage        ', deadstemc_storage        ) 
     call check_vector_data ('deadstemc_xfer           ', deadstemc_xfer           ) 
     call check_vector_data ('livecrootc               ', livecrootc               ) 
     call check_vector_data ('livecrootc_storage       ', livecrootc_storage       ) 
     call check_vector_data ('livecrootc_xfer          ', livecrootc_xfer          ) 
     call check_vector_data ('deadcrootc               ', deadcrootc               ) 
     call check_vector_data ('deadcrootc_storage       ', deadcrootc_storage       ) 
     call check_vector_data ('deadcrootc_xfer          ', deadcrootc_xfer          ) 
     call check_vector_data ('grainc                   ', grainc                   ) 
     call check_vector_data ('grainc_storage           ', grainc_storage           ) 
     call check_vector_data ('grainc_xfer              ', grainc_xfer              ) 
     call check_vector_data ('xsmrpool                 ', xsmrpool                 ) 
     call check_vector_data ('downreg                  ', downreg                  ) 
     call check_vector_data ('cropprod1c               ', cropprod1c               ) 
     call check_vector_data ('cropseedc_deficit        ', cropseedc_deficit        ) 

     call check_vector_data ('leafn                    ', leafn                    ) 
     call check_vector_data ('leafn_storage            ', leafn_storage            ) 
     call check_vector_data ('leafn_xfer               ', leafn_xfer               ) 
     call check_vector_data ('frootn                   ', frootn                   ) 
     call check_vector_data ('frootn_storage           ', frootn_storage           ) 
     call check_vector_data ('frootn_xfer              ', frootn_xfer              ) 
     call check_vector_data ('livestemn                ', livestemn                ) 
     call check_vector_data ('livestemn_storage        ', livestemn_storage        ) 
     call check_vector_data ('livestemn_xfer           ', livestemn_xfer           ) 
     call check_vector_data ('deadstemn                ', deadstemn                ) 
     call check_vector_data ('deadstemn_storage        ', deadstemn_storage        ) 
     call check_vector_data ('deadstemn_xfer           ', deadstemn_xfer           ) 
     call check_vector_data ('livecrootn               ', livecrootn               ) 
     call check_vector_data ('livecrootn_storage       ', livecrootn_storage       ) 
     call check_vector_data ('livecrootn_xfer          ', livecrootn_xfer          ) 
     call check_vector_data ('deadcrootn               ', deadcrootn               ) 
     call check_vector_data ('deadcrootn_storage       ', deadcrootn_storage       ) 
     call check_vector_data ('deadcrootn_xfer          ', deadcrootn_xfer          ) 
     call check_vector_data ('grainn                   ', grainn                   ) 
     call check_vector_data ('grainn_storage           ', grainn_storage           ) 
     call check_vector_data ('grainn_xfer              ', grainn_xfer              ) 
     call check_vector_data ('retransn                 ', retransn                 ) 

     call check_vector_data ('decomp_npools_vr         ', decomp_npools_vr         ) 
     call check_vector_data ('decomp_npools            ', decomp_npools            ) 
     call check_vector_data ('ntrunc_vr                ', ntrunc_vr                ) 
     call check_vector_data ('ntrunc_veg               ', ntrunc_veg               ) 
     call check_vector_data ('ntrunc_soil              ', ntrunc_soil              ) 

     call check_vector_data ('sminn_vr                 ', sminn_vr                 ) 
     call check_vector_data ('smin_no3_vr              ', smin_no3_vr              ) 
     call check_vector_data ('smin_nh4_vr              ', smin_nh4_vr              ) 
     call check_vector_data ('sminn                    ', sminn                    ) 

     call check_vector_data ('ndep_prof                ', ndep_prof                ) 
     call check_vector_data ('nfixation_prof           ', nfixation_prof           ) 

     call check_vector_data ('cn_decomp_pools          ', cn_decomp_pools          ) 
     call check_vector_data ('fpi_vr                   ', fpi_vr                   ) 
     call check_vector_data ('fpi                      ', fpi                      ) 
     call check_vector_data ('fpg                      ', fpg                      ) 

     call check_vector_data ('cropf                    ', cropf                    ) 
     call check_vector_data ('lfwt                     ', lfwt                     ) 
     call check_vector_data ('fuelc                    ', fuelc                    ) 
     call check_vector_data ('fuelc_crop               ', fuelc_crop               ) 
     call check_vector_data ('fsr                      ', fsr                      ) 
     call check_vector_data ('fd                       ', fd                       ) 
     call check_vector_data ('rootc                    ', rootc                    ) 
     call check_vector_data ('lgdp                     ', lgdp                     ) 
     call check_vector_data ('lgdp1                    ', lgdp1                    ) 
     call check_vector_data ('lpop                     ', lpop                     ) 
     call check_vector_data ('wtlf                     ', wtlf                     ) 
     call check_vector_data ('trotr1                   ', trotr1                   ) 
     call check_vector_data ('trotr2                   ', trotr2                   ) 
     call check_vector_data ('hdmlf                    ', hdmlf                    ) 
     call check_vector_data ('lnfm                     ', lnfm                     ) 
     call check_vector_data ('baf_crop                 ', baf_crop                 ) 
     call check_vector_data ('baf_peatf                ', baf_peatf                ) 
     call check_vector_data ('farea_burned             ', farea_burned             ) 
     call check_vector_data ('nfire                    ', nfire                    ) 
     call check_vector_data ('fsat                     ', fsat                     ) 
     call check_vector_data ('prec10                   ', prec10                   ) 
     call check_vector_data ('prec60                   ', prec60                   ) 
     call check_vector_data ('prec365                  ', prec365                  ) 
     call check_vector_data ('prec_today               ', prec_today               ) 
     call check_vector_data ('prec_daily               ', prec_daily               ) 
     call check_vector_data ('wf2                      ', wf2                      ) 
     call check_vector_data ('tsoi17                   ', tsoi17                   ) 
     call check_vector_data ('rh30                     ', rh30                     ) 
     call check_vector_data ('accumnstep               ', accumnstep               ) 
     call check_vector_data ('cphase                   ', cphase                   ) 

     call check_vector_data ('dayl                     ', dayl                     ) 
     call check_vector_data ('prev_dayl                ', prev_dayl                ) 

#ifdef SASU
!--------------SASU variables---------------------------
     call check_vector_data ('decomp0_cpools_vr           ', decomp0_cpools_vr           ) 
     call check_vector_data ('I_met_c_vr_acc              ', I_met_c_vr_acc              ) 
     call check_vector_data ('I_cel_c_vr_acc              ', I_cel_c_vr_acc              ) 
     call check_vector_data ('I_lig_c_vr_acc              ', I_lig_c_vr_acc              ) 
     call check_vector_data ('I_cwd_c_vr_acc              ', I_cwd_c_vr_acc              ) 
     call check_vector_data ('AKX_met_to_soil1_c_vr_acc   ', AKX_met_to_soil1_c_vr_acc   ) 
     call check_vector_data ('AKX_cel_to_soil1_c_vr_acc   ', AKX_cel_to_soil1_c_vr_acc   ) 
     call check_vector_data ('AKX_lig_to_soil2_c_vr_acc   ', AKX_lig_to_soil2_c_vr_acc   ) 
     call check_vector_data ('AKX_soil1_to_soil2_c_vr_acc ', AKX_soil1_to_soil2_c_vr_acc ) 
     call check_vector_data ('AKX_cwd_to_cel_c_vr_acc     ', AKX_cwd_to_cel_c_vr_acc     ) 
     call check_vector_data ('AKX_cwd_to_lig_c_vr_acc     ', AKX_cwd_to_lig_c_vr_acc     ) 
     call check_vector_data ('AKX_soil1_to_soil3_c_vr_acc ', AKX_soil1_to_soil3_c_vr_acc ) 
     call check_vector_data ('AKX_soil2_to_soil1_c_vr_acc ', AKX_soil2_to_soil1_c_vr_acc ) 
     call check_vector_data ('AKX_soil2_to_soil3_c_vr_acc ', AKX_soil2_to_soil3_c_vr_acc ) 
     call check_vector_data ('AKX_soil3_to_soil1_c_vr_acc ', AKX_soil3_to_soil1_c_vr_acc ) 
     call check_vector_data ('AKX_met_exit_c_vr_acc       ', AKX_met_exit_c_vr_acc       ) 
     call check_vector_data ('AKX_cel_exit_c_vr_acc       ', AKX_cel_exit_c_vr_acc       ) 
     call check_vector_data ('AKX_lig_exit_c_vr_acc       ', AKX_lig_exit_c_vr_acc       ) 
     call check_vector_data ('AKX_cwd_exit_c_vr_acc       ', AKX_cwd_exit_c_vr_acc       ) 
     call check_vector_data ('AKX_soil1_exit_c_vr_acc     ', AKX_soil1_exit_c_vr_acc     ) 
     call check_vector_data ('AKX_soil2_exit_c_vr_acc     ', AKX_soil2_exit_c_vr_acc     ) 
     call check_vector_data ('AKX_soil3_exit_c_vr_acc     ', AKX_soil3_exit_c_vr_acc     ) 

     call check_vector_data ('decomp0_npools_vr           ', decomp0_npools_vr           ) 
     call check_vector_data ('I_met_n_vr_acc              ', I_met_n_vr_acc              ) 
     call check_vector_data ('I_cel_n_vr_acc              ', I_cel_n_vr_acc              ) 
     call check_vector_data ('I_lig_n_vr_acc              ', I_lig_n_vr_acc              ) 
     call check_vector_data ('I_cwd_n_vr_acc              ', I_cwd_n_vr_acc              ) 
     call check_vector_data ('AKX_met_to_soil1_n_vr_acc   ', AKX_met_to_soil1_n_vr_acc   ) 
     call check_vector_data ('AKX_cel_to_soil1_n_vr_acc   ', AKX_cel_to_soil1_n_vr_acc   ) 
     call check_vector_data ('AKX_lig_to_soil2_n_vr_acc   ', AKX_lig_to_soil2_n_vr_acc   ) 
     call check_vector_data ('AKX_soil1_to_soil2_n_vr_acc ', AKX_soil1_to_soil2_n_vr_acc ) 
     call check_vector_data ('AKX_cwd_to_cel_n_vr_acc     ', AKX_cwd_to_cel_n_vr_acc     ) 
     call check_vector_data ('AKX_cwd_to_lig_n_vr_acc     ', AKX_cwd_to_lig_n_vr_acc     ) 
     call check_vector_data ('AKX_soil1_to_soil3_n_vr_acc ', AKX_soil1_to_soil3_n_vr_acc ) 
     call check_vector_data ('AKX_soil2_to_soil1_n_vr_acc ', AKX_soil2_to_soil1_n_vr_acc ) 
     call check_vector_data ('AKX_soil2_to_soil3_n_vr_acc ', AKX_soil2_to_soil3_n_vr_acc ) 
     call check_vector_data ('AKX_soil3_to_soil1_n_vr_acc ', AKX_soil3_to_soil1_n_vr_acc ) 
     call check_vector_data ('AKX_met_exit_n_vr_acc       ', AKX_met_exit_n_vr_acc       ) 
     call check_vector_data ('AKX_cel_exit_n_vr_acc       ', AKX_cel_exit_n_vr_acc       ) 
     call check_vector_data ('AKX_lig_exit_n_vr_acc       ', AKX_lig_exit_n_vr_acc       ) 
     call check_vector_data ('AKX_cwd_exit_n_vr_acc       ', AKX_cwd_exit_n_vr_acc       ) 
     call check_vector_data ('AKX_soil1_exit_n_vr_acc     ', AKX_soil1_exit_n_vr_acc     ) 
     call check_vector_data ('AKX_soil2_exit_n_vr_acc     ', AKX_soil2_exit_n_vr_acc     ) 
     call check_vector_data ('AKX_soil3_exit_n_vr_acc     ', AKX_soil3_exit_n_vr_acc     ) 

     call check_vector_data ('diagVX_c_vr_acc             ', diagVX_c_vr_acc             ) 
     call check_vector_data ('upperVX_c_vr_acc            ', upperVX_c_vr_acc            ) 
     call check_vector_data ('lowerVX_c_vr_acc            ', lowerVX_c_vr_acc            ) 
     call check_vector_data ('diagVX_n_vr_acc             ', diagVX_n_vr_acc             ) 
     call check_vector_data ('upperVX_n_vr_acc            ', upperVX_n_vr_acc            ) 
     call check_vector_data ('lowerVX_n_vr_acc            ', lowerVX_n_vr_acc            ) 
!     call check_vector_data ('skip_balance_check          ', skip_balance_check          ) 
!------------------------------------------------------
#endif
#endif

#if (defined PFT_CLASSIFICATION)
     CALL check_PFTimeVars
#endif

#if (defined PC_CLASSIFICATION)
     CALL check_PCTimeVars
#endif

#ifdef USEMPI
     call mpi_barrier (p_comm_glb, p_err)
#endif
     
  end subroutine check_TimeVariables
#endif


END MODULE MOD_TimeVariables
! ------ EOP --------------
