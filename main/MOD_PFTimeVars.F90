#include <define.h> 

MODULE MOD_PFTimeVars
! -------------------------------
! Created by Hua Yuan, 08/2019
! -------------------------------

  USE precision
  USE timemanager

  IMPLICIT NONE
  SAVE
! -----------------------------------------------------------------
! Time-varying state variables which reaquired by restart run

  ! for PFT_CLASSIFICATION
  REAL(r8), allocatable :: tleaf_p   (:) !shaded leaf temperature [K]
  REAL(r8), allocatable :: ldew_p    (:) !depth of water on foliage [mm]
  REAL(r8), allocatable :: sigf_p    (:) !fraction of veg cover, excluding snow-covered veg [-]
  REAL(r8), allocatable :: tlai_p    (:) !leaf area index
  REAL(r8), allocatable :: lai_p     (:) !leaf area index
  REAL(r8), allocatable :: laisun_p  (:) !sunlit leaf area index
  REAL(r8), allocatable :: laisha_p  (:) !shaded leaf area index
  REAL(r8), allocatable :: tsai_p    (:) !stem area index
  REAL(r8), allocatable :: sai_p     (:) !stem area index                                      
  REAL(r8), allocatable :: ssun_p(:,:,:) !sunlit canopy absorption for solar radiation (0-1)
  REAL(r8), allocatable :: ssha_p(:,:,:) !shaded canopy absorption for solar radiation (0-1)
  REAL(r8), allocatable :: thermk_p  (:) !canopy gap fraction for tir radiation
  REAL(r8), allocatable :: extkb_p   (:) !(k, g(mu)/mu) direct solar extinction coefficient
  REAL(r8), allocatable :: extkd_p   (:) !diffuse and scattered diffuse PAR extinction coefficient
  REAL(r8), allocatable :: tref_p    (:) !2 m height air temperature [kelvin]
  REAL(r8), allocatable :: qref_p    (:) !2 m height air specific humidity
  REAL(r8), allocatable :: rst_p     (:) !canopy stomatal resistance (s/m)
  REAL(r8), allocatable :: z0m_p     (:) !effective roughness [m]

  ! bgc variables
  REAL(r8), allocatable :: leafc_p                  (:)
  REAL(r8), allocatable :: leafc_storage_p          (:)
  REAL(r8), allocatable :: leafc_xfer_p             (:)
  REAL(r8), allocatable :: frootc_p                 (:)
  REAL(r8), allocatable :: frootc_storage_p         (:)
  REAL(r8), allocatable :: frootc_xfer_p            (:)
  REAL(r8), allocatable :: livestemc_p              (:)
  REAL(r8), allocatable :: livestemc_storage_p      (:)
  REAL(r8), allocatable :: livestemc_xfer_p         (:)
  REAL(r8), allocatable :: deadstemc_p              (:)
  REAL(r8), allocatable :: deadstemc_storage_p      (:)
  REAL(r8), allocatable :: deadstemc_xfer_p         (:)
  REAL(r8), allocatable :: livecrootc_p             (:)
  REAL(r8), allocatable :: livecrootc_storage_p     (:)
  REAL(r8), allocatable :: livecrootc_xfer_p        (:)
  REAL(r8), allocatable :: deadcrootc_p             (:)
  REAL(r8), allocatable :: deadcrootc_storage_p     (:)
  REAL(r8), allocatable :: deadcrootc_xfer_p        (:)
  REAL(r8), allocatable :: grainc_p                 (:)
  REAL(r8), allocatable :: grainc_storage_p         (:)
  REAL(r8), allocatable :: grainc_xfer_p            (:)
  REAL(r8), allocatable :: cropseedc_deficit_p      (:)
  REAL(r8), allocatable :: xsmrpool_p               (:)
  REAL(r8), allocatable :: gresp_storage_p          (:)
  REAL(r8), allocatable :: gresp_xfer_p             (:)

  REAL(r8), allocatable :: leaf_prof_p              (:,:) !profile of leaves
  REAL(r8), allocatable :: stem_prof_p              (:,:) !profile of stem
  REAL(r8), allocatable :: froot_prof_p             (:,:) !profile of fine roots
  REAL(r8), allocatable :: croot_prof_p             (:,:) !profile of coarse roots
  REAL(r8), allocatable :: crootfr                  (:,:)

  REAL(r8), allocatable :: leafn_p                  (:)
  REAL(r8), allocatable :: leafn_storage_p          (:)
  REAL(r8), allocatable :: leafn_xfer_p             (:)
  REAL(r8), allocatable :: frootn_p                 (:)
  REAL(r8), allocatable :: frootn_storage_p         (:)
  REAL(r8), allocatable :: frootn_xfer_p            (:)
  REAL(r8), allocatable :: livestemn_p              (:)
  REAL(r8), allocatable :: livestemn_storage_p      (:)
  REAL(r8), allocatable :: livestemn_xfer_p         (:)
  REAL(r8), allocatable :: deadstemn_p              (:)
  REAL(r8), allocatable :: deadstemn_storage_p      (:)
  REAL(r8), allocatable :: deadstemn_xfer_p         (:)
  REAL(r8), allocatable :: livecrootn_p             (:)
  REAL(r8), allocatable :: livecrootn_storage_p     (:)
  REAL(r8), allocatable :: livecrootn_xfer_p        (:)
  REAL(r8), allocatable :: deadcrootn_p             (:)
  REAL(r8), allocatable :: deadcrootn_storage_p     (:)
  REAL(r8), allocatable :: deadcrootn_xfer_p        (:)
  REAL(r8), allocatable :: grainn_p                 (:)
  REAL(r8), allocatable :: grainn_storage_p         (:)
  REAL(r8), allocatable :: grainn_xfer_p            (:)
  REAL(r8), allocatable :: cropseedn_deficit_p      (:)
  REAL(r8), allocatable :: retransn_p               (:)

  REAL(r8), allocatable :: harvdate_p               (:)

  REAL(r8), allocatable :: tempsum_potential_gpp_p  (:)
  REAL(r8), allocatable :: tempmax_retransn_p       (:)
  REAL(r8), allocatable :: tempavg_tref_p           (:)
  REAL(r8), allocatable :: tempsum_npp_p            (:)
  REAL(r8), allocatable :: tempsum_litfall_p        (:)
  REAL(r8), allocatable :: annsum_potential_gpp_p   (:)
  REAL(r8), allocatable :: annmax_retransn_p        (:)
  REAL(r8), allocatable :: annavg_tref_p            (:)
  REAL(r8), allocatable :: annsum_npp_p             (:)
  REAL(r8), allocatable :: annsum_litfall_p         (:)

  REAL(r8), allocatable :: bglfr_p                  (:)
  REAL(r8), allocatable :: bgtr_p                   (:)
  REAL(r8), allocatable :: lgsf_p                   (:)
  REAL(r8), allocatable :: gdd0_p                   (:)
  REAL(r8), allocatable :: gdd8_p                   (:)
  REAL(r8), allocatable :: gdd10_p                  (:)
  REAL(r8), allocatable :: gdd020_p                 (:)
  REAL(r8), allocatable :: gdd820_p                 (:)
  REAL(r8), allocatable :: gdd1020_p                (:)
  INTEGER , allocatable :: nyrs_crop_active_p       (:)

  REAL(r8), allocatable :: offset_flag_p            (:)
  REAL(r8), allocatable :: offset_counter_p         (:)
  REAL(r8), allocatable :: onset_flag_p             (:)
  REAL(r8), allocatable :: onset_counter_p          (:)
  REAL(r8), allocatable :: onset_gddflag_p          (:)
  REAL(r8), allocatable :: onset_gdd_p              (:)
  REAL(r8), allocatable :: onset_fdd_p              (:)
  REAL(r8), allocatable :: onset_swi_p              (:)
  REAL(r8), allocatable :: offset_fdd_p             (:)
  REAL(r8), allocatable :: offset_swi_p             (:)
  REAL(r8), allocatable :: dormant_flag_p           (:)
  REAL(r8), allocatable :: prev_leafc_to_litter_p   (:)
  REAL(r8), allocatable :: prev_frootc_to_litter_p  (:)
  REAL(r8), allocatable :: days_active_p            (:)

  REAL(r8), allocatable :: burndate_p               (:)

  REAL(r8), allocatable :: c_allometry_p            (:)
  REAL(r8), allocatable :: n_allometry_p            (:)
  REAL(r8), allocatable :: downreg_p                (:)
  REAL(r8), allocatable :: grain_flag_p             (:)

!-----

! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: allocate_PFTimeVars
  PUBLIC :: deallocate_PFTimeVars

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE allocate_PFTimeVars ()
! ------------------------------------------------------
! Allocates memory for CLM 1d [numpft] variables
! ------------------------------------------------------
      USE precision
      USE GlobalVars
      IMPLICIT NONE

      allocate (tleaf_p      (numpft)) !leaf temperature [K]
      allocate (ldew_p       (numpft)) !depth of water on foliage [mm]
      allocate (sigf_p       (numpft)) !fraction of veg cover, excluding snow-covered veg [-]
      allocate (tlai_p       (numpft)) !leaf area index
      allocate (lai_p        (numpft)) !leaf area index
      allocate (laisun_p     (numpft)) !sunlit leaf area index
      allocate (laisha_p     (numpft)) !shaded leaf area index
      allocate (tsai_p       (numpft)) !stem area index
      allocate (sai_p        (numpft)) !stem area index
      allocate (ssun_p   (2,2,numpft)) !sunlit canopy absorption for solar radiation (0-1)
      allocate (ssha_p   (2,2,numpft)) !shaded canopy absorption for solar radiation (0-1)
      allocate (thermk_p     (numpft)) !canopy gap fraction for tir radiation
      allocate (extkb_p      (numpft)) !(k, g(mu)/mu) direct solar extinction coefficient
      allocate (extkd_p      (numpft)) !diffuse and scattered diffuse PAR extinction coefficient
      allocate (tref_p       (numpft)) !2 m height air temperature [kelvin]
      allocate (qref_p       (numpft)) !2 m height air specific humidity
      allocate (rst_p        (numpft)) !canopy stomatal resistance (s/m)
      allocate (z0m_p        (numpft)) !effective roughness [m]

  ! bgc variables
      allocate (leafc_p                  (numpft))
      allocate (leafc_storage_p          (numpft))
      allocate (leafc_xfer_p             (numpft))
      allocate (frootc_p                 (numpft))
      allocate (frootc_storage_p         (numpft))
      allocate (frootc_xfer_p            (numpft))
      allocate (livestemc_p              (numpft))
      allocate (livestemc_storage_p      (numpft))
      allocate (livestemc_xfer_p         (numpft))
      allocate (deadstemc_p              (numpft))
      allocate (deadstemc_storage_p      (numpft))
      allocate (deadstemc_xfer_p         (numpft))
      allocate (livecrootc_p             (numpft))
      allocate (livecrootc_storage_p     (numpft))
      allocate (livecrootc_xfer_p        (numpft))
      allocate (deadcrootc_p             (numpft))
      allocate (deadcrootc_storage_p     (numpft))
      allocate (deadcrootc_xfer_p        (numpft))
      allocate (grainc_p                 (numpft))
      allocate (grainc_storage_p         (numpft))
      allocate (grainc_xfer_p            (numpft))
      allocate (cropseedc_deficit_p      (numpft))
      allocate (xsmrpool_p               (numpft))
      allocate (gresp_storage_p          (numpft))
      allocate (gresp_xfer_p             (numpft))

      allocate (leaf_prof_p              (nl_soil,numpft))
      allocate (froot_prof_p             (nl_soil,numpft))
      allocate (croot_prof_p             (nl_soil,numpft))
      allocate (stem_prof_p              (nl_soil,numpft))
      allocate (crootfr                  (nl_soil,numpatch))

      allocate (leafn_p                  (numpft))
      allocate (leafn_storage_p          (numpft))
      allocate (leafn_xfer_p             (numpft))
      allocate (frootn_p                 (numpft))
      allocate (frootn_storage_p         (numpft))
      allocate (frootn_xfer_p            (numpft))
      allocate (livestemn_p              (numpft))
      allocate (livestemn_storage_p      (numpft))
      allocate (livestemn_xfer_p         (numpft))
      allocate (deadstemn_p              (numpft))
      allocate (deadstemn_storage_p      (numpft))
      allocate (deadstemn_xfer_p         (numpft))
      allocate (livecrootn_p             (numpft))
      allocate (livecrootn_storage_p     (numpft))
      allocate (livecrootn_xfer_p        (numpft))
      allocate (deadcrootn_p             (numpft))
      allocate (deadcrootn_storage_p     (numpft))
      allocate (deadcrootn_xfer_p        (numpft))
      allocate (grainn_p                 (numpft))
      allocate (grainn_storage_p         (numpft))
      allocate (grainn_xfer_p            (numpft))
      allocate (cropseedn_deficit_p      (numpft))
      allocate (retransn_p               (numpft))

      allocate (harvdate_p               (numpft))

      allocate (tempsum_potential_gpp_p  (numpft))
      allocate (tempmax_retransn_p       (numpft))
      allocate (tempavg_tref_p           (numpft))
      allocate (tempsum_npp_p            (numpft))
      allocate (tempsum_litfall_p        (numpft))
      allocate (annsum_potential_gpp_p   (numpft))
      allocate (annmax_retransn_p        (numpft))
      allocate (annavg_tref_p            (numpft))
      allocate (annsum_npp_p             (numpft))
      allocate (annsum_litfall_p         (numpft))

      allocate (bglfr_p                  (numpft))
      allocate (bgtr_p                   (numpft))
      allocate (lgsf_p                   (numpft))
      allocate (gdd0_p                   (numpft))
      allocate (gdd8_p                   (numpft))
      allocate (gdd10_p                  (numpft))
      allocate (gdd020_p                 (numpft))
      allocate (gdd820_p                 (numpft))
      allocate (gdd1020_p                (numpft))
      allocate (nyrs_crop_active_p       (numpft))

      allocate (offset_flag_p            (numpft))
      allocate (offset_counter_p         (numpft))
      allocate (onset_flag_p             (numpft))
      allocate (onset_counter_p          (numpft))
      allocate (onset_gddflag_p          (numpft))
      allocate (onset_gdd_p              (numpft))
      allocate (onset_fdd_p              (numpft))
      allocate (onset_swi_p              (numpft))
      allocate (offset_fdd_p             (numpft))
      allocate (offset_swi_p             (numpft))
      allocate (dormant_flag_p           (numpft))
      allocate (prev_leafc_to_litter_p   (numpft))
      allocate (prev_frootc_to_litter_p  (numpft))
      allocate (days_active_p            (numpft))

      allocate (burndate_p               (numpft))

      allocate (c_allometry_p            (numpft))
      allocate (n_allometry_p            (numpft))
      allocate (downreg_p                (numpft))
      allocate (grain_flag_p             (numpft))

   END SUBROUTINE allocate_PFTimeVars 
  
   SUBROUTINE deallocate_PFTimeVars ()
! --------------------------------------------------
! Deallocates memory for CLM 1d [numpft/numpc] variables
! --------------------------------------------------
      deallocate (tleaf_p  ) !leaf temperature [K]
      deallocate (ldew_p   ) !depth of water on foliage [mm]
      deallocate (sigf_p   ) !fraction of veg cover, excluding snow-covered veg [-]
      deallocate (tlai_p   ) !leaf area index
      deallocate (lai_p    ) !leaf area index
      deallocate (laisun_p ) !sunlit leaf area index
      deallocate (laisha_p ) !shaded leaf area index
      deallocate (tsai_p   ) !stem area index
      deallocate (sai_p    ) !stem area index                                   
      deallocate (ssun_p   ) !sunlit canopy absorption for solar radiation (0-1)
      deallocate (ssha_p   ) !shaded canopy absorption for solar radiation (0-1)
      deallocate (thermk_p ) !canopy gap fraction for tir radiation
      deallocate (extkb_p  ) !(k, g(mu)/mu) direct solar extinction coefficient
      deallocate (extkd_p  ) !diffuse and scattered diffuse PAR extinction coefficient
      deallocate (tref_p   ) !2 m height air temperature [kelvin]
      deallocate (qref_p   ) !2 m height air specific humidity
      deallocate (rst_p    ) !canopy stomatal resistance (s/m)
      deallocate (z0m_p    ) !effective roughness [m]                                 

! bgc variables
      deallocate (leafc_p                  )
      deallocate (leafc_storage_p          )
      deallocate (leafc_xfer_p             )
      deallocate (frootc_p                 )
      deallocate (frootc_storage_p         )
      deallocate (frootc_xfer_p            )
      deallocate (livestemc_p              )
      deallocate (livestemc_storage_p      )
      deallocate (livestemc_xfer_p         )
      deallocate (deadstemc_p              )
      deallocate (deadstemc_storage_p      )
      deallocate (deadstemc_xfer_p         )
      deallocate (livecrootc_p             )
      deallocate (livecrootc_storage_p     )
      deallocate (livecrootc_xfer_p        )
      deallocate (deadcrootc_p             )
      deallocate (deadcrootc_storage_p     )
      deallocate (deadcrootc_xfer_p        )
      deallocate (grainc_p                 )
      deallocate (grainc_storage_p         )
      deallocate (grainc_xfer_p            )
      deallocate (cropseedc_deficit_p      )
      deallocate (xsmrpool_p               )
      deallocate (gresp_storage_p          )
      deallocate (gresp_xfer_p             )

      deallocate (leaf_prof_p              )
      deallocate (froot_prof_p             )
      deallocate (croot_prof_p             )
      deallocate (stem_prof_p              )
      deallocate (crootfr                  )

      deallocate (leafn_p                  )
      deallocate (leafn_storage_p          )
      deallocate (leafn_xfer_p             )
      deallocate (frootn_p                 )
      deallocate (frootn_storage_p         )
      deallocate (frootn_xfer_p            )
      deallocate (livestemn_p              )
      deallocate (livestemn_storage_p      )
      deallocate (livestemn_xfer_p         )
      deallocate (deadstemn_p              )
      deallocate (deadstemn_storage_p      )
      deallocate (deadstemn_xfer_p         )
      deallocate (livecrootn_p             )
      deallocate (livecrootn_storage_p     )
      deallocate (livecrootn_xfer_p        )
      deallocate (deadcrootn_p             )
      deallocate (deadcrootn_storage_p     )
      deallocate (deadcrootn_xfer_p        )
      deallocate (grainn_p                 )
      deallocate (grainn_storage_p         )
      deallocate (grainn_xfer_p            )
      deallocate (cropseedn_deficit_p      )
      deallocate (retransn_p               )

      deallocate (harvdate_p               )

      deallocate (tempsum_potential_gpp_p  )
      deallocate (tempmax_retransn_p       )
      deallocate (tempavg_tref_p           )
      deallocate (tempsum_npp_p            )
      deallocate (tempsum_litfall_p        )
      deallocate (annsum_potential_gpp_p   )
      deallocate (annmax_retransn_p        )
      deallocate (annavg_tref_p            )
      deallocate (annsum_npp_p             )
      deallocate (annsum_litfall_p         )

      deallocate (bglfr_p                  )
      deallocate (bgtr_p                   )
      deallocate (lgsf_p                   )
      deallocate (gdd0_p                   )
      deallocate (gdd8_p                   )
      deallocate (gdd10_p                  )
      deallocate (gdd020_p                 )
      deallocate (gdd820_p                 )
      deallocate (gdd1020_p                )
      deallocate (nyrs_crop_active_p       )

      deallocate (offset_flag_p            )
      deallocate (offset_counter_p         )
      deallocate (onset_flag_p             )
      deallocate (onset_counter_p          )
      deallocate (onset_gddflag_p          )
      deallocate (onset_gdd_p              )
      deallocate (onset_fdd_p              )
      deallocate (onset_swi_p              )
      deallocate (offset_fdd_p             )
      deallocate (offset_swi_p             )
      deallocate (dormant_flag_p           )
      deallocate (prev_leafc_to_litter_p   )
      deallocate (prev_frootc_to_litter_p  )
      deallocate (days_active_p            )

      deallocate (burndate_p               )

      deallocate (c_allometry_p            )
      deallocate (n_allometry_p            )
      deallocate (downreg_p                )
      deallocate (grain_flag_p             )

   END SUBROUTINE deallocate_PFTimeVars

END MODULE MOD_PFTimeVars
! ---------- EOP ------------
