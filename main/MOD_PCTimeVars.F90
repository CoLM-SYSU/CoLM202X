#include <define.h> 

MODULE MOD_PCTimeVars
! -------------------------------
! Created by Hua Yuan, 08/2019
! -------------------------------

  USE precision
  USE timemanager
  IMPLICIT NONE
  SAVE
! -----------------------------------------------------------------
! Time-varying state variables which reaquired by restart run

  ! for PC_CLASSIFICATION
  REAL(r8), allocatable :: tleaf_c    (:,:) !leaf temperature [K]
  REAL(r8), allocatable :: ldew_c     (:,:) !depth of water on foliage [mm]
  REAL(r8), allocatable :: sigf_c     (:,:) !fraction of veg cover, excluding snow-covered veg [-]
  REAL(r8), allocatable :: tlai_c     (:,:) !leaf area index
  REAL(r8), allocatable :: lai_c      (:,:) !leaf area index
  REAL(r8), allocatable :: tsai_c     (:,:) !stem area index
  REAL(r8), allocatable :: sai_c      (:,:) !stem area index
  REAL(r8), allocatable :: ssun_c (:,:,:,:) !sunlit canopy absorption for solar radiation (0-1)
  REAL(r8), allocatable :: ssha_c (:,:,:,:) !shaded canopy absorption for solar radiation (0-1)
  REAL(r8), allocatable :: thermk_c   (:,:) !canopy gap fraction for tir radiation
  REAL(r8), allocatable :: fshade_c   (:,:) !canopy gap fraction for tir radiation
  REAL(r8), allocatable :: extkb_c    (:,:) !(k, g(mu)/mu) direct solar extinction coefficient
  REAL(r8), allocatable :: extkd_c    (:,:) !diffuse and scattered diffuse PAR extinction coefficient
  REAL(r8), allocatable :: rst_c      (:,:) !canopy stomatal resistance (s/m)
  REAL(r8), allocatable :: z0m_c      (:,:) !effective roughness [m]

! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: allocate_PCTimeVars
  PUBLIC :: deallocate_PCTimeVars

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE allocate_PCTimeVars ()
! ------------------------------------------------------
! Allocates memory for CLM Plant Community (PC) 1D [numpc] variables
! ------------------------------------------------------
      USE precision
      USE GlobalVars
      IMPLICIT NONE

      allocate (tleaf_c    (0:N_PFT-1,numpc)) !leaf temperature [K]
      allocate (ldew_c     (0:N_PFT-1,numpc)) !depth of water on foliage [mm]
      allocate (sigf_c     (0:N_PFT-1,numpc)) !fraction of veg cover, excluding snow-covered veg [-]
      allocate (tlai_c     (0:N_PFT-1,numpc)) !leaf area index
      allocate (lai_c      (0:N_PFT-1,numpc)) !leaf area index
      allocate (tsai_c     (0:N_PFT-1,numpc)) !stem area index
      allocate (sai_c      (0:N_PFT-1,numpc)) !stem area index
      allocate (ssun_c (2,2,0:N_PFT-1,numpc)) !sunlit canopy absorption for solar radiation (0-1)
      allocate (ssha_c (2,2,0:N_PFT-1,numpc)) !shaded canopy absorption for solar radiation (0-1)
      allocate (thermk_c   (0:N_PFT-1,numpc)) !canopy gap fraction for tir radiation
      allocate (fshade_c   (0:N_PFT-1,numpc)) !canopy gap fraction for tir radiation
      allocate (extkb_c    (0:N_PFT-1,numpc)) !(k, g(mu)/mu) direct solar extinction coefficient
      allocate (extkd_c    (0:N_PFT-1,numpc)) !diffuse and scattered diffuse PAR extinction coefficient
      allocate (rst_c      (0:N_PFT-1,numpc)) !canopy stomatal resistance (s/m)
      allocate (z0m_c      (0:N_PFT-1,numpc)) !effective roughness [m]
 
   END SUBROUTINE allocate_PCTimeVars
  
  
   SUBROUTINE deallocate_PCTimeVars
! --------------------------------------------------
! Deallocates memory for CLM Plant Community (PC) 1D [numpc] variables
! --------------------------------------------------
      deallocate (tleaf_c  ) !leaf temperature [K]
      deallocate (ldew_c   ) !depth of water on foliage [mm]
      deallocate (sigf_c   ) !fraction of veg cover, excluding snow-covered veg [-]
      deallocate (tlai_c   ) !leaf area index
      deallocate (lai_c    ) !leaf area index
      deallocate (tsai_c   ) !stem area index
      deallocate (sai_c    ) !stem area index
      deallocate (ssun_c   ) !sunlit canopy absorption for solar radiation (0-1)
      deallocate (ssha_c   ) !shaded canopy absorption for solar radiation (0-1)
      deallocate (thermk_c ) !canopy gap fraction for tir radiation
      deallocate (fshade_c ) !canopy gap fraction for tir radiation
      deallocate (extkb_c  ) !(k, g(mu)/mu) direct solar extinction coefficient
      deallocate (extkd_c  ) !diffuse and scattered diffuse PAR extinction coefficient
      deallocate (rst_c    ) !canopy stomatal resistance (s/m)
      deallocate (z0m_c    ) !effective roughness [m]                                 
   END SUBROUTINE deallocate_PCTimeVars

END MODULE MOD_PCTimeVars
! ---------- EOP ------------
