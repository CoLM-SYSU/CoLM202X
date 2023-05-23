#include <define.h>

#if (defined PC_CLASSIFICATION)

MODULE MOD_Vars_PCTimeVars
! -----------------------------------------------------------------
! !DESCRIPTION:
! Define Plant Community time variables
!
! Created by Hua Yuan, 08/2019
! -----------------------------------------------------------------

  USE precision
  USE timemanager
  IMPLICIT NONE
  SAVE
! -----------------------------------------------------------------
! Time-varying state variables which reaquired by restart run

  ! for PC_CLASSIFICATION
  REAL(r8), allocatable :: tleaf_c    (:,:) !leaf temperature [K]
  REAL(r8), allocatable :: ldew_c     (:,:) !depth of water on foliage [mm]
!#ifdef CLM5_INTERCEPTION
  real(r8), allocatable :: ldew_rain_c     (:,:)     !depth of rain on foliage [mm]
  real(r8), allocatable :: ldew_snow_c     (:,:)     !depth of rain on foliage [mm]
!#endif
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
#ifdef PLANT_HYDRAULIC_STRESS
  real(r8), allocatable :: vegwp_c  (:,:,:) !vegetation water potential [mm]
  real(r8), allocatable :: gs0sun_c   (:,:) !working copy of sunlit stomata conductance
  real(r8), allocatable :: gs0sha_c   (:,:) !working copy of shalit stomata conductance
#endif

! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: allocate_PCTimeVars
  PUBLIC :: deallocate_PCTimeVars
  PUBLIC :: READ_PCTimeVars
  PUBLIC :: WRITE_PCTimeVars
#ifdef CoLMDEBUG
  PUBLIC :: check_PCTimeVars
#endif

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE allocate_PCTimeVars ()
! ------------------------------------------------------
! Allocates memory for CoLM Plant Community (PC) 1D [numpc] variables
! ------------------------------------------------------
      USE precision
      USE GlobalVars
      USE spmd_task
      USE MOD_LandPC
      IMPLICIT NONE

      IF (p_is_worker) THEN
         IF (numpc > 0) THEN
            allocate (tleaf_c    (0:N_PFT-1,numpc)) !leaf temperature [K]
            allocate (ldew_c     (0:N_PFT-1,numpc)) !depth of water on foliage [mm]
            allocate (ldew_rain_c(0:N_PFT-1,numpc)) !depth of rain on foliage [mm]
            allocate (ldew_snow_c(0:N_PFT-1,numpc)) !depth of snow on foliage [mm]
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
#ifdef PLANT_HYDRAULIC_STRESS
            allocate (vegwp_c    (1:nvegwcs,0:N_PFT-1,numpc))
            allocate (gs0sun_c   (0:N_PFT-1,numpc))
            allocate (gs0sha_c   (0:N_PFT-1,numpc))
#endif
         ENDIF
      ENDIF

   END SUBROUTINE allocate_PCTimeVars

   SUBROUTINE READ_PCTimeVars (file_restart)

      USE GlobalVars
      use mod_namelist
      use ncio_vector
      USE MOD_LandPC
      IMPLICIT NONE

      character(LEN=*), intent(in) :: file_restart

      call ncio_read_vector (file_restart, 'tleaf_c  ', N_PFT,     landpc, tleaf_c  ) !
      call ncio_read_vector (file_restart, 'ldew_c   ', N_PFT,     landpc, ldew_c   ) !
      call ncio_read_vector (file_restart, 'ldew_rain_c', N_PFT,   landpc, ldew_rain_c) !depth of rain on foliage [mm]
      call ncio_read_vector (file_restart, 'ldew_snow_c', N_PFT,   landpc, ldew_snow_c) !depth of snow on foliage [mm]
      call ncio_read_vector (file_restart, 'sigf_c   ', N_PFT,     landpc, sigf_c   ) !
      call ncio_read_vector (file_restart, 'tlai_c   ', N_PFT,     landpc, tlai_c   ) !
      call ncio_read_vector (file_restart, 'lai_c    ', N_PFT,     landpc, lai_c    ) !
      call ncio_read_vector (file_restart, 'tsai_c   ', N_PFT,     landpc, tsai_c   ) !
      call ncio_read_vector (file_restart, 'sai_c    ', N_PFT,     landpc, sai_c    ) !
      call ncio_read_vector (file_restart, 'ssun_c   ', 2,2,N_PFT, landpc, ssun_c   ) !
      call ncio_read_vector (file_restart, 'ssha_c   ', 2,2,N_PFT, landpc, ssha_c   ) !
      call ncio_read_vector (file_restart, 'thermk_c ', N_PFT,     landpc, thermk_c ) !
      call ncio_read_vector (file_restart, 'fshade_c ', N_PFT,     landpc, fshade_c ) !
      call ncio_read_vector (file_restart, 'extkb_c  ', N_PFT,     landpc, extkb_c  ) !
      call ncio_read_vector (file_restart, 'extkd_c  ', N_PFT,     landpc, extkd_c  ) !
      call ncio_read_vector (file_restart, 'rst_c    ', N_PFT,     landpc, rst_c    ) !
      call ncio_read_vector (file_restart, 'z0m_c    ', N_PFT,     landpc, z0m_c    ) !
#ifdef PLANT_HYDRAULIC_STRESS
      call ncio_read_vector (file_restart, 'vegwp_c  ', nvegwcs,   N_PFT,  landpc, vegwp_c ) !
      call ncio_read_vector (file_restart, 'gs0sun_c ', N_PFT,     landpc, gs0sun_c ) !
      call ncio_read_vector (file_restart, 'gs0sha_c ', N_PFT,     landpc, gs0sha_c ) !
#endif

   END SUBROUTINE READ_PCTimeVars

   SUBROUTINE WRITE_PCTimeVars (file_restart)

     USE GlobalVars
     use mod_namelist, only : DEF_REST_COMPRESS_LEVEL
     USE MOD_LandPC
     use ncio_vector
     IMPLICIT NONE

     character(LEN=*), intent(in) :: file_restart

     ! Local variables
     integer :: compress

     compress = DEF_REST_COMPRESS_LEVEL

     call ncio_create_file_vector      (file_restart, landpc               )
     CALL ncio_define_dimension_vector (file_restart, landpc, 'pc'         )
     CALL ncio_define_dimension_vector (file_restart, landpc, 'pft' , N_PFT)
     CALL ncio_define_dimension_vector (file_restart, landpc, 'band', 2    )
     CALL ncio_define_dimension_vector (file_restart, landpc, 'rtyp', 2    )
#ifdef PLANT_HYDRAULIC_STRESS
     CALL ncio_define_dimension_vector (file_restart, landpc, 'vegnodes', nvegwcs)
#endif

      call ncio_write_vector (file_restart, 'tleaf_c  ', 'pft', N_PFT, 'pc', landpc, tleaf_c  , compress) !
      call ncio_write_vector (file_restart, 'ldew_c   ', 'pft', N_PFT, 'pc', landpc, ldew_c   , compress) !
      call ncio_write_vector (file_restart, 'ldew_rain_c', 'pft', N_PFT, 'pc', landpc, ldew_rain_c, compress) ! depth of rain on foliage [mm]
      call ncio_write_vector (file_restart, 'ldew_snow_c', 'pft', N_PFT, 'pc', landpc, ldew_snow_c, compress) ! depth of snow on foliage [mm]

      call ncio_write_vector (file_restart, 'sigf_c   ', 'pft', N_PFT, 'pc', landpc, sigf_c   , compress) !
      call ncio_write_vector (file_restart, 'tlai_c   ', 'pft', N_PFT, 'pc', landpc, tlai_c   , compress) !
      call ncio_write_vector (file_restart, 'lai_c    ', 'pft', N_PFT, 'pc', landpc, lai_c    , compress) !
      call ncio_write_vector (file_restart, 'tsai_c   ', 'pft', N_PFT, 'pc', landpc, tsai_c   , compress) !
      call ncio_write_vector (file_restart, 'sai_c    ', 'pft', N_PFT, 'pc', landpc, sai_c    , compress) !
      call ncio_write_vector (file_restart, 'ssun_c   ', 'band', 2, 'rtyp', 2, 'pft', N_PFT, 'pc', landpc, ssun_c, compress) !
      call ncio_write_vector (file_restart, 'ssha_c   ', 'band', 2, 'rtyp', 2, 'pft', N_PFT, 'pc', landpc, ssha_c, compress) !
      call ncio_write_vector (file_restart, 'thermk_c ', 'pft', N_PFT, 'pc', landpc, thermk_c , compress) !
      call ncio_write_vector (file_restart, 'fshade_c ', 'pft', N_PFT, 'pc', landpc, fshade_c , compress) !
      call ncio_write_vector (file_restart, 'extkb_c  ', 'pft', N_PFT, 'pc', landpc, extkb_c  , compress) !
      call ncio_write_vector (file_restart, 'extkd_c  ', 'pft', N_PFT, 'pc', landpc, extkd_c  , compress) !
      call ncio_write_vector (file_restart, 'rst_c    ', 'pft', N_PFT, 'pc', landpc, rst_c    , compress) !
      call ncio_write_vector (file_restart, 'z0m_c    ', 'pft', N_PFT, 'pc', landpc, z0m_c    , compress) !
#ifdef PLANT_HYDRAULIC_STRESS
      call ncio_write_vector (file_restart, 'vegwp_c  ', 'vegnodes', nvegwcs, 'pft', N_PFT , 'pc'    , landpc, vegwp_c, compress)
      call ncio_write_vector (file_restart, 'gs0sun_c ', 'pft'     , N_PFT  , 'pc' , landpc, gs0sun_c, compress) !
      call ncio_write_vector (file_restart, 'gs0sha_c ', 'pft'     , N_PFT  , 'pc' , landpc, gs0sha_c, compress) !
#endif

   END SUBROUTINE WRITE_PCTimeVars


   SUBROUTINE deallocate_PCTimeVars
! --------------------------------------------------
! Deallocates memory for CoLM Plant Community (PC) 1D [numpc] variables
! --------------------------------------------------

      USE spmd_task
      USE MOD_LandPC

      IF (p_is_worker) THEN
         IF (numpc > 0) THEN
            deallocate (tleaf_c  ) !leaf temperature [K]
            deallocate (ldew_c   ) !depth of water on foliage [mm]
            deallocate (ldew_rain_c ) !depth of water on foliage [mm]
            deallocate (ldew_snow_c ) !depth of water on foliage [mm]
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
#ifdef PLANT_HYDRAULIC_STRESS
            deallocate (vegwp_c  ) !vegetation water potential [mm]
            deallocate (gs0sun_c ) !working copy of sunlit stomata conductance
            deallocate (gs0sha_c ) !working copy of shalit stomata conductance
#endif
         ENDIF
      ENDIF

   END SUBROUTINE deallocate_PCTimeVars

#ifdef CoLMDEBUG
   SUBROUTINE check_PCTimeVars

      use mod_colm_debug
      IMPLICIT NONE

      call check_vector_data ('tleaf_c  ', tleaf_c  )      !
      call check_vector_data ('ldew_c   ', ldew_c   )      !
      call check_vector_data ('ldew_rain_c', ldew_rain_c)  !  depth of rain on foliage [mm]
      call check_vector_data ('ldew_snow_c', ldew_snow_c)  !  depth of snow on foliage [mm]
      call check_vector_data ('sigf_c   ', sigf_c   )      !
      call check_vector_data ('tlai_c   ', tlai_c   )      !
      call check_vector_data ('lai_c    ', lai_c    )      !
      call check_vector_data ('tsai_c   ', tsai_c   )      !
      call check_vector_data ('sai_c    ', sai_c    )      !
      call check_vector_data ('ssun_c   ', ssun_c   )      !
      call check_vector_data ('ssha_c   ', ssha_c   )      !
      call check_vector_data ('thermk_c ', thermk_c )      !
      call check_vector_data ('fshade_c ', fshade_c )      !
      call check_vector_data ('extkb_c  ', extkb_c  )      !
      call check_vector_data ('extkd_c  ', extkd_c  )      !
      call check_vector_data ('rst_c    ', rst_c    )      !
      call check_vector_data ('z0m_c    ', z0m_c    )      !
#ifdef PLANT_HYDRAULIC_STRESS
      call check_vector_data ('vegwp_c  ', vegwp_c  )      !
      call check_vector_data ('gs0sun_c ', gs0sun_c )      !
      call check_vector_data ('gs0sha_c ', gs0sha_c )      !
#endif

   END SUBROUTINE check_PCTimeVars
#endif


END MODULE MOD_Vars_PCTimeVars

#endif
! ---------- EOP ------------
