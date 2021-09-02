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

   END SUBROUTINE allocate_PFTimeVars
  
  
   SUBROUTINE deallocate_PFTimeVars
! --------------------------------------------------
! Deallocates memory for CLM 1d [numpft/numpc] variables
! --------------------------------------------------
      deallocate (tleaf_p  ) !leaf temperature [K]
      deallocate (ldew_p   ) !depth of water on foliage [mm]
      deallocate (sigf_p   ) !fraction of veg cover, excluding snow-covered veg [-]
      deallocate (tlai_p   ) !leaf area index
      deallocate (lai_p    ) !leaf area index
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
   END SUBROUTINE deallocate_PFTimeVars

END MODULE MOD_PFTimeVars
! ---------- EOP ------------
