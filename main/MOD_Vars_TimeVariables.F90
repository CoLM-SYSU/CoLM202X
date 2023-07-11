#include <define.h>

! -------------------------------
! Created by Yongjiu Dai, 03/2014
! -------------------------------

#if (defined LULC_IGBP_PFT)
MODULE MOD_Vars_PFTimeVariables
! -----------------------------------------------------------------
! !DESCRIPTION:
! Define PFT time variables
!
! Added by Hua Yuan, 08/2019
! -----------------------------------------------------------------

  USE MOD_Precision
  USE MOD_TimeManager
#ifdef BGC
  USE MOD_BGC_Vars_PFTimeVariables
#endif

  IMPLICIT NONE
  SAVE
! -----------------------------------------------------------------
! Time-varying state variables which reaquired by restart run

  ! for LULC_IGBP_PFT
  REAL(r8), allocatable :: tleaf_p      (:) !shaded leaf temperature [K]
  REAL(r8), allocatable :: ldew_p       (:) !depth of water on foliage [mm]
  real(r8), allocatable :: ldew_rain_p  (:)!depth of rain on foliage [mm]
  real(r8), allocatable :: ldew_snow_p  (:)!depth of snow on foliage [mm]
  REAL(r8), allocatable :: sigf_p       (:) !fraction of veg cover, excluding snow-covered veg [-]
  REAL(r8), allocatable :: tlai_p       (:) !leaf area index
  REAL(r8), allocatable :: lai_p        (:) !leaf area index
  REAL(r8), allocatable :: laisun_p     (:) !sunlit leaf area index
  REAL(r8), allocatable :: laisha_p     (:) !shaded leaf area index
  REAL(r8), allocatable :: tsai_p       (:) !stem area index
  REAL(r8), allocatable :: sai_p        (:) !stem area index
  REAL(r8), allocatable :: ssun_p   (:,:,:) !sunlit canopy absorption for solar radiation (0-1)
  REAL(r8), allocatable :: ssha_p   (:,:,:) !shaded canopy absorption for solar radiation (0-1)
  REAL(r8), allocatable :: thermk_p     (:) !canopy gap fraction for tir radiation
  REAL(r8), allocatable :: extkb_p      (:) !(k, g(mu)/mu) direct solar extinction coefficient
  REAL(r8), allocatable :: extkd_p      (:) !diffuse and scattered diffuse PAR extinction coefficient
  REAL(r8), allocatable :: tref_p       (:) !2 m height air temperature [kelvin]
  REAL(r8), allocatable :: qref_p       (:) !2 m height air specific humidity
  REAL(r8), allocatable :: rst_p        (:) !canopy stomatal resistance (s/m)
  REAL(r8), allocatable :: z0m_p        (:) !effective roughness [m]
! Plant Hydraulic variables
  real(r8), allocatable :: vegwp_p    (:,:) ! vegetation water potential [mm]
  real(r8), allocatable :: gs0sun_p     (:) ! working copy of sunlit stomata conductance
  real(r8), allocatable :: gs0sha_p     (:) ! working copy of shalit stomata conductance
! end plant hydraulic variables
!Ozone Stress Variables
  real(r8), allocatable :: o3coefv_sun_p(:) !Ozone stress factor for photosynthesis on sunlit leaf
  real(r8), allocatable :: o3coefv_sha_p(:) !Ozone stress factor for photosynthesis on shaded leaf
  real(r8), allocatable :: o3coefg_sun_p(:) !Ozone stress factor for stomata on sunlit leaf
  real(r8), allocatable :: o3coefg_sha_p(:) !Ozone stress factor for stomata on shaded leaf
  real(r8), allocatable :: lai_old_p    (:) !lai in last time step
  real(r8), allocatable :: o3uptakesun_p(:) !Ozone does, sunlit leaf (mmol O3/m^2)
  real(r8), allocatable :: o3uptakesha_p(:) !Ozone does, shaded leaf (mmol O3/m^2)
!End Ozone Stress Variables

! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: allocate_PFTimeVariables
  PUBLIC :: deallocate_PFTimeVariables
  PUBLIC :: READ_PFTimeVariables
  PUBLIC :: WRITE_PFTimeVariables
#ifdef CoLMDEBUG
  PUBLIC :: check_PFTimeVariables
#endif

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE allocate_PFTimeVariables ()
! ------------------------------------------------------
! Allocates memory for CoLM 1d [numpft] variables
! ------------------------------------------------------
      USE MOD_Precision
      USE MOD_SPMD_Task
      USE MOD_LandPFT
      USE MOD_Vars_Global
      IMPLICIT NONE

      IF (p_is_worker) THEN
         IF (numpft > 0) THEN
            allocate (tleaf_p      (numpft)) !leaf temperature [K]
            allocate (ldew_p       (numpft)) !depth of water on foliage [mm]
            allocate (ldew_rain_p  (numpft)) !depth of rain on foliage [mm]
            allocate (ldew_snow_p  (numpft)) !depth of snow on foliage [mm]
            allocate (sigf_p       (numpft)) !fraction of veg cover, excluding snow-covered veg [-]
            allocate (tlai_p       (numpft)) !leaf area index
            allocate (lai_p        (numpft)) !leaf area index
            allocate (laisun_p     (numpft)) !leaf area index
            allocate (laisha_p     (numpft)) !leaf area index
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
! Plant Hydraulic variables
            allocate (vegwp_p      (1:nvegwcs,numpft))
            allocate (gs0sun_p     (numpft))
            allocate (gs0sha_p     (numpft))
! end plant hydraulic variables
! Allocate Ozone Stress Variables
            allocate (o3coefv_sun_p(numpft)) !Ozone stress factor for photosynthesis on sunlit leaf
            allocate (o3coefv_sha_p(numpft)) !Ozone stress factor for photosynthesis on shaded leaf
            allocate (o3coefg_sun_p(numpft)) !Ozone stress factor for stomata on sunlit leaf
            allocate (o3coefg_sha_p(numpft)) !Ozone stress factor for stomata on shaded leaf
            allocate (lai_old_p    (numpft)) !lai in last time step
            allocate (o3uptakesun_p(numpft)) !Ozone does, sunlit leaf (mmol O3/m^2)
            allocate (o3uptakesha_p(numpft)) !Ozone does, shaded leaf (mmol O3/m^2)
! End allocate Ozone Stress Variables
         ENDIF
      ENDIF

#ifdef BGC
      CALL allocate_BGCPFTimeVariables
#endif

   END SUBROUTINE allocate_PFTimeVariables

   SUBROUTINE READ_PFTimeVariables (file_restart)

      USE MOD_Namelist, only: DEF_USE_PLANTHYDRAULICS, DEF_USE_OZONESTRESS
      use MOD_NetCDFVector
      USE MOD_LandPFT
      USE MOD_Vars_Global

      IMPLICIT NONE

      character(LEN=*), intent(in) :: file_restart

      call ncio_read_vector (file_restart, 'tleaf_p  ', landpft, tleaf_p    ) !
      call ncio_read_vector (file_restart, 'ldew_p   ', landpft, ldew_p     ) !
      call ncio_read_vector (file_restart, 'ldew_rain_p', landpft, ldew_rain_p) !depth of rain on foliage [mm]
      call ncio_read_vector (file_restart, 'ldew_snow_p', landpft, ldew_snow_p) !depth of snow on foliage [mm]
      call ncio_read_vector (file_restart, 'sigf_p   ', landpft, sigf_p     ) !
      call ncio_read_vector (file_restart, 'tlai_p   ', landpft, tlai_p     ) !
      call ncio_read_vector (file_restart, 'lai_p    ', landpft, lai_p      ) !
      call ncio_read_vector (file_restart, 'laisun_p ', landpft, laisun_p   ) !
      call ncio_read_vector (file_restart, 'laisha_p ', landpft, laisha_p   ) !
      call ncio_read_vector (file_restart, 'tsai_p   ', landpft, tsai_p     ) !
      call ncio_read_vector (file_restart, 'sai_p    ', landpft, sai_p      ) !
      call ncio_read_vector (file_restart, 'ssun_p   ', 2,2, landpft, ssun_p) !
      call ncio_read_vector (file_restart, 'ssha_p   ', 2,2, landpft, ssha_p) !
      call ncio_read_vector (file_restart, 'thermk_p ', landpft, thermk_p   ) !
      call ncio_read_vector (file_restart, 'extkb_p  ', landpft, extkb_p    ) !
      call ncio_read_vector (file_restart, 'extkd_p  ', landpft, extkd_p    ) !
      call ncio_read_vector (file_restart, 'tref_p   ', landpft, tref_p     ) !
      call ncio_read_vector (file_restart, 'qref_p   ', landpft, qref_p     ) !
      call ncio_read_vector (file_restart, 'rst_p    ', landpft, rst_p      ) !
      call ncio_read_vector (file_restart, 'z0m_p    ', landpft, z0m_p      ) !
      IF(DEF_USE_PLANTHYDRAULICS)THEN
         call ncio_read_vector (file_restart, 'vegwp_p  ', nvegwcs, landpft, vegwp_p ) !
         call ncio_read_vector (file_restart, 'gs0sun_p ', landpft, gs0sun_p   ) !
         call ncio_read_vector (file_restart, 'gs0sha_p ', landpft, gs0sha_p   ) !
      END IF
      IF(DEF_USE_OZONESTRESS)THEN
         call ncio_read_vector (file_restart, 'lai_old_p    ', landpft, lai_old_p    , defval = 0._r8)
         call ncio_read_vector (file_restart, 'o3uptakesun_p', landpft, o3uptakesun_p, defval = 0._r8)
         call ncio_read_vector (file_restart, 'o3uptakesha_p', landpft, o3uptakesha_p, defval = 0._r8)
      ENDIF

#ifdef BGC
      CALL read_BGCPFTimeVariables (file_restart)
#endif

   END SUBROUTINE READ_PFTimeVariables

   SUBROUTINE WRITE_PFTimeVariables (file_restart)

     use MOD_Namelist, only : DEF_REST_COMPRESS_LEVEL, DEF_USE_PLANTHYDRAULICS, DEF_USE_OZONESTRESS
     USE MOD_LandPFT
     use MOD_NetCDFVector
     USE MOD_Vars_Global
     IMPLICIT NONE

     character(LEN=*), intent(in) :: file_restart

     ! Local variables
     integer :: compress

     compress = DEF_REST_COMPRESS_LEVEL

     call ncio_create_file_vector (file_restart, landpft)
     CALL ncio_define_dimension_vector (file_restart, landpft, 'pft')
     CALL ncio_define_dimension_vector (file_restart, landpft, 'band', 2)
     CALL ncio_define_dimension_vector (file_restart, landpft, 'rtyp', 2)
     if(DEF_USE_PLANTHYDRAULICS)then
        CALL ncio_define_dimension_vector (file_restart, landpft, 'vegnodes', nvegwcs)
     end if

     call ncio_write_vector (file_restart, 'tleaf_p  ', 'pft', landpft, tleaf_p  , compress) !
     call ncio_write_vector (file_restart, 'ldew_p   ', 'pft', landpft, ldew_p   , compress) !
     call ncio_write_vector (file_restart, 'ldew_rain_p', 'pft', landpft, ldew_rain_p, compress) !depth of rain on foliage [mm]
     call ncio_write_vector (file_restart, 'ldew_snow_p', 'pft', landpft, ldew_snow_p, compress) !depth of snow on foliage [mm]
     call ncio_write_vector (file_restart, 'sigf_p   ', 'pft', landpft, sigf_p   , compress) !
     call ncio_write_vector (file_restart, 'tlai_p   ', 'pft', landpft, tlai_p   , compress) !
     call ncio_write_vector (file_restart, 'lai_p    ', 'pft', landpft, lai_p    , compress) !
     call ncio_write_vector (file_restart, 'laisun_p ', 'pft', landpft, laisun_p , compress) !
     call ncio_write_vector (file_restart, 'laisha_p ', 'pft', landpft, laisha_p , compress) !
     call ncio_write_vector (file_restart, 'tsai_p   ', 'pft', landpft, tsai_p   , compress) !
     call ncio_write_vector (file_restart, 'sai_p    ', 'pft', landpft, sai_p    , compress) !
     call ncio_write_vector (file_restart, 'ssun_p   ', 'band', 2, 'rtyp', 2, 'pft', landpft, ssun_p, compress) !
     call ncio_write_vector (file_restart, 'ssha_p   ', 'band', 2, 'rtyp', 2, 'pft', landpft, ssha_p, compress) !
     call ncio_write_vector (file_restart, 'thermk_p ', 'pft', landpft, thermk_p , compress) !
     call ncio_write_vector (file_restart, 'extkb_p  ', 'pft', landpft, extkb_p  , compress) !
     call ncio_write_vector (file_restart, 'extkd_p  ', 'pft', landpft, extkd_p  , compress) !
     call ncio_write_vector (file_restart, 'tref_p   ', 'pft', landpft, tref_p   , compress) !
     call ncio_write_vector (file_restart, 'qref_p   ', 'pft', landpft, qref_p   , compress) !
     call ncio_write_vector (file_restart, 'rst_p    ', 'pft', landpft, rst_p    , compress) !
     call ncio_write_vector (file_restart, 'z0m_p    ', 'pft', landpft, z0m_p    , compress) !
     IF(DEF_USE_PLANTHYDRAULICS)then
        call ncio_write_vector (file_restart, 'vegwp_p  '  , 'vegnodes', nvegwcs, 'pft',   landpft, vegwp_p, compress)
        call ncio_write_vector (file_restart, 'gs0sun_p '  , 'pft', landpft, gs0sun_p   , compress) !
        call ncio_write_vector (file_restart, 'gs0sha_p '  , 'pft', landpft, gs0sha_p   , compress) !
     END IF
     IF(DEF_USE_OZONESTRESS)THEN
        call ncio_write_vector (file_restart, 'lai_old_p    ', 'pft', landpft, lai_old_p    , compress)
        call ncio_write_vector (file_restart, 'o3uptakesun_p', 'pft', landpft, o3uptakesun_p, compress)
        call ncio_write_vector (file_restart, 'o3uptakesha_p', 'pft', landpft, o3uptakesha_p, compress)
     ENDIF

#ifdef BGC
      CALL WRITE_BGCPFTimeVariables (file_restart)
#endif

   END SUBROUTINE WRITE_PFTimeVariables


   SUBROUTINE deallocate_PFTimeVariables
! --------------------------------------------------
! Deallocates memory for CoLM 1d [numpft/numpc] variables
! --------------------------------------------------
      USE MOD_SPMD_Task
      USE MOD_LandPFT

      IF (p_is_worker) THEN
         IF (numpft > 0) THEN
            deallocate (tleaf_p  ) !leaf temperature [K]
            deallocate (ldew_p   ) !depth of water on foliage [mm]
            deallocate (ldew_rain_p)
            deallocate (ldew_snow_p)
            deallocate (sigf_p   ) !fraction of veg cover, excluding snow-covered veg [-]
            deallocate (tlai_p   ) !leaf area index
            deallocate (lai_p    ) !leaf area index
            deallocate (laisun_p ) !leaf area index
            deallocate (laisha_p ) !leaf area index
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
! Plant Hydraulic variables
            deallocate (vegwp_p  ) !vegetation water potential [mm]
            deallocate (gs0sun_p ) !working copy of sunlit stomata conductance
            deallocate (gs0sha_p ) !working copy of shalit stomata conductance
! end plant hydraulic variables
! Ozone Stress variables
            deallocate (o3coefv_sun_p ) !Ozone stress factor for photosynthesis on sunlit leaf
            deallocate (o3coefv_sha_p ) !Ozone stress factor for photosynthesis on shaded leaf
            deallocate (o3coefg_sun_p ) !Ozone stress factor for stomata on sunlit leaf
            deallocate (o3coefg_sha_p ) !Ozone stress factor for stomata on shaded leaf
            deallocate (lai_old_p     ) !lai in last time step
            deallocate (o3uptakesun_p ) !Ozone does, sunlit leaf (mmol O3/m^2)
            deallocate (o3uptakesha_p ) !Ozone does, shaded leaf (mmol O3/m^2)
! Ozone Stress variables
         ENDIF
      ENDIF

#ifdef BGC
      CALL deallocate_BGCPFTimeVariables
#endif

   END SUBROUTINE deallocate_PFTimeVariables

#ifdef CoLMDEBUG
   SUBROUTINE check_PFTimeVariables

      use MOD_CoLMDebug
      use MOD_Namelist, only : DEF_USE_PLANTHYDRAULICS, DEF_USE_OZONESTRESS

      IMPLICIT NONE

      call check_vector_data ('tleaf_p  ', tleaf_p  )      !
      call check_vector_data ('ldew_p   ', ldew_p   )      !
      call check_vector_data ('ldew_rain_p', ldew_rain_p ) !depth of rain on foliage [mm]
      call check_vector_data ('ldew_snow_p', ldew_snow_p ) !depth of snow on foliage [mm]
      call check_vector_data ('sigf_p   ', sigf_p   )      !
      call check_vector_data ('tlai_p   ', tlai_p   )      !
      call check_vector_data ('lai_p    ', lai_p    )      !
      call check_vector_data ('laisun_p ', lai_p    )      !
      call check_vector_data ('laisha_p ', lai_p    )      !
      call check_vector_data ('tsai_p   ', tsai_p   )      !
      call check_vector_data ('sai_p    ', sai_p    )      !
      call check_vector_data ('ssun_p   ', ssun_p   )      !
      call check_vector_data ('ssha_p   ', ssha_p   )      !
      call check_vector_data ('thermk_p ', thermk_p )      !
      call check_vector_data ('extkb_p  ', extkb_p  )      !
      call check_vector_data ('extkd_p  ', extkd_p  )      !
      call check_vector_data ('tref_p   ', tref_p   )      !
      call check_vector_data ('qref_p   ', qref_p   )      !
      call check_vector_data ('rst_p    ', rst_p    )      !
      call check_vector_data ('z0m_p    ', z0m_p    )      !
      IF(DEF_USE_PLANTHYDRAULICS)THEN
         call check_vector_data ('vegwp_p  ', vegwp_p  )      !
         call check_vector_data ('gs0sun_p ', gs0sun_p )      !
         call check_vector_data ('gs0sha_p ', gs0sha_p )      !
      ENDIF
      IF(DEF_USE_OZONESTRESS)THEN
         call check_vector_data ('o3coefv_sun_p', o3coefv_sun_p)
         call check_vector_data ('o3coefv_sha_p', o3coefv_sha_p)
         call check_vector_data ('o3coefg_sun_p', o3coefg_sun_p)
         call check_vector_data ('o3coefg_sha_p', o3coefg_sha_p)
         call check_vector_data ('lai_old_p    ', lai_old_p    )
         call check_vector_data ('o3uptakesun_p', o3uptakesun_p)
         call check_vector_data ('o3uptakesha_p', o3uptakesha_p)
      ENDIF

#ifdef BGC
      CALL check_BGCPFTimeVariables
#endif

   END SUBROUTINE check_PFTimeVariables
#endif

END MODULE MOD_Vars_PFTimeVariables
#endif



#if (defined LULC_IGBP_PC)
MODULE MOD_Vars_PCTimeVariables
! -----------------------------------------------------------------
! !DESCRIPTION:
! Define Plant Community time variables
!
! Added by Hua Yuan, 08/2019
! -----------------------------------------------------------------

  USE MOD_Precision
  USE MOD_TimeManager
  IMPLICIT NONE
  SAVE
! -----------------------------------------------------------------
! Time-varying state variables which reaquired by restart run

  ! for LULC_IGBP_PC
  REAL(r8), allocatable :: tleaf_c    (:,:) !leaf temperature [K]
  REAL(r8), allocatable :: ldew_c     (:,:) !depth of water on foliage [mm]
!#ifdef CLM5_INTERCEPTION
  real(r8), allocatable :: ldew_rain_c(:,:) !depth of rain on foliage [mm]
  real(r8), allocatable :: ldew_snow_c(:,:) !depth of rain on foliage [mm]
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
!Plant Hydraulic parameters
  real(r8), allocatable :: vegwp_c  (:,:,:) !vegetation water potential [mm]
  real(r8), allocatable :: gs0sun_c   (:,:) !working copy of sunlit stomata conductance
  real(r8), allocatable :: gs0sha_c   (:,:) !working copy of shalit stomata conductance
!end plant hydraulic parameters
!Ozone Stress Variables
  real(r8), allocatable :: o3coefv_sun_c(:,:) !Ozone stress factor for photosynthesis on sunlit leaf
  real(r8), allocatable :: o3coefv_sha_c(:,:) !Ozone stress factor for photosynthesis on shaded leaf
  real(r8), allocatable :: o3coefg_sun_c(:,:) !Ozone stress factor for stomata on sunlit leaf
  real(r8), allocatable :: o3coefg_sha_c(:,:) !Ozone stress factor for stomata on shaded leaf
  real(r8), allocatable :: lai_old_c    (:,:) !lai in last time step
  real(r8), allocatable :: o3uptakesun_c(:,:) !Ozone does, sunlit leaf (mmol O3/m^2)
  real(r8), allocatable :: o3uptakesha_c(:,:) !Ozone does, shaded leaf (mmol O3/m^2)
!End Ozone Stress Variables

! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: allocate_PCTimeVariables
  PUBLIC :: deallocate_PCTimeVariables
  PUBLIC :: READ_PCTimeVariables
  PUBLIC :: WRITE_PCTimeVariables
#ifdef CoLMDEBUG
  PUBLIC :: check_PCTimeVariables
#endif

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE allocate_PCTimeVariables ()
! ------------------------------------------------------
! Allocates memory for CoLM Plant Community (PC) 1D [numpc] variables
! ------------------------------------------------------
      USE MOD_Precision
      USE MOD_Vars_Global
      USE MOD_SPMD_Task
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
!Plant Hydraulic parameters
            allocate (vegwp_c    (1:nvegwcs,0:N_PFT-1,numpc))
            allocate (gs0sun_c   (0:N_PFT-1,numpc))
            allocate (gs0sha_c   (0:N_PFT-1,numpc))
!end plant hydraulic parameters
!Ozone Stress Variables
            allocate (o3coefv_sun_c(0:N_PFT-1,numpc)) !Ozone stress factor for photosynthesis on sunlit leaf
            allocate (o3coefv_sha_c(0:N_PFT-1,numpc)) !Ozone stress factor for photosynthesis on shaded leaf
            allocate (o3coefg_sun_c(0:N_PFT-1,numpc)) !Ozone stress factor for stomata on sunlit leaf
            allocate (o3coefg_sha_c(0:N_PFT-1,numpc)) !Ozone stress factor for stomata on shaded leaf
            allocate (lai_old_c    (0:N_PFT-1,numpc)) !lai in last time step
            allocate (o3uptakesun_c(0:N_PFT-1,numpc)) !Ozone does, sunlit leaf (mmol O3/m^2)
            allocate (o3uptakesha_c(0:N_PFT-1,numpc)) !Ozone does, shaded leaf (mmol O3/m^2)
!End Ozone Stress Variables
         ENDIF
      ENDIF

   END SUBROUTINE allocate_PCTimeVariables

   SUBROUTINE READ_PCTimeVariables (file_restart)

      USE MOD_Vars_Global
      use MOD_Namelist, only: DEF_USE_PLANTHYDRAULICS, DEF_USE_OZONESTRESS
      use MOD_NetCDFVector
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
      if(DEF_USE_PLANTHYDRAULICS)then
         call ncio_read_vector (file_restart, 'vegwp_c  ', nvegwcs,   N_PFT,  landpc, vegwp_c ) !
         call ncio_read_vector (file_restart, 'gs0sun_c ', N_PFT,     landpc, gs0sun_c ) !
         call ncio_read_vector (file_restart, 'gs0sha_c ', N_PFT,     landpc, gs0sha_c ) !
      end if
      IF(DEF_USE_OZONESTRESS)THEN
!Ozone Stress Variables
         call ncio_read_vector (file_restart, 'o3coefv_sun_c', N_PFT, landpc, o3coefv_sun_c)!Ozone stress factor for photosynthesis on sunlit leaf
         call ncio_read_vector (file_restart, 'o3coefv_sha_c', N_PFT, landpc, o3coefv_sha_c) !Ozone stress factor for photosynthesis on shaded leaf
         call ncio_read_vector (file_restart, 'o3coefg_sun_c', N_PFT, landpc, o3coefg_sun_c) !Ozone stress factor for stomata on sunlit leaf
         call ncio_read_vector (file_restart, 'o3coefg_sha_c', N_PFT, landpc, o3coefg_sha_c) !Ozone stress factor for stomata on shaded leaf
         call ncio_read_vector (file_restart, 'lai_old_c    ', N_PFT, landpc, lai_old_c    ) !lai in last time step
         call ncio_read_vector (file_restart, 'o3uptakesun_c', N_PFT, landpc, o3uptakesun_c) !Ozone does, sunlit leaf (mmol O3/m^2)
         call ncio_read_vector (file_restart, 'o3uptakesha_c', N_PFT, landpc, o3uptakesha_c) !Ozone does, shaded leaf (mmol O3/m^2)
!End Ozone Stress Variables
      ENDIF

   END SUBROUTINE READ_PCTimeVariables

   SUBROUTINE WRITE_PCTimeVariables (file_restart)

      USE MOD_Vars_Global
      use MOD_Namelist, only : DEF_REST_COMPRESS_LEVEL, DEF_USE_PLANTHYDRAULICS, DEF_USE_OZONESTRESS
      USE MOD_LandPC
      use MOD_NetCDFVector
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
      if(DEF_USE_PLANTHYDRAULICS)then
         CALL ncio_define_dimension_vector (file_restart, landpc, 'vegnodes', nvegwcs)
      end if

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
      if(DEF_USE_PLANTHYDRAULICS)then
         call ncio_write_vector (file_restart, 'vegwp_c  ', 'vegnodes', nvegwcs, 'pft', N_PFT , 'pc'    , landpc, vegwp_c, compress)
         call ncio_write_vector (file_restart, 'gs0sun_c ', 'pft'     , N_PFT  , 'pc' , landpc, gs0sun_c, compress) !
         call ncio_write_vector (file_restart, 'gs0sha_c ', 'pft'     , N_PFT  , 'pc' , landpc, gs0sha_c, compress) !
      end if
      IF(DEF_USE_OZONESTRESS)THEN
!Ozone Stress Variables
         call ncio_write_vector (file_restart, 'o3coefv_sun_c', 'pft' , N_PFT  , 'pc' , landpc, o3coefv_sun_c, compress)!Ozone stress factor for photosynthesis on sunlit leaf
         call ncio_write_vector (file_restart, 'o3coefv_sha_c', 'pft' , N_PFT  , 'pc' , landpc, o3coefv_sha_c, compress) !Ozone stress factor for photosynthesis on shaded leaf
         call ncio_write_vector (file_restart, 'o3coefg_sun_c', 'pft' , N_PFT  , 'pc' , landpc, o3coefg_sun_c, compress) !Ozone stress factor for stomata on sunlit leaf
         call ncio_write_vector (file_restart, 'o3coefg_sha_c', 'pft' , N_PFT  , 'pc' , landpc, o3coefg_sha_c, compress) !Ozone stress factor for stomata on shaded leaf
         call ncio_write_vector (file_restart, 'lai_old_c    ', 'pft' , N_PFT  , 'pc' , landpc, lai_old_c    , compress) !lai in last time step
         call ncio_write_vector (file_restart, 'o3uptakesun_c', 'pft' , N_PFT  , 'pc' , landpc, o3uptakesun_c, compress) !Ozone does, sunlit leaf (mmol O3/m^2)
         call ncio_write_vector (file_restart, 'o3uptakesha_c', 'pft' , N_PFT  , 'pc' , landpc, o3uptakesha_c, compress) !Ozone does, shaded leaf (mmol O3/m^2)
      ENDIF

   END SUBROUTINE WRITE_PCTimeVariables


   SUBROUTINE deallocate_PCTimeVariables
! --------------------------------------------------
! Deallocates memory for CoLM Plant Community (PC) 1D [numpc] variables
! --------------------------------------------------

      USE MOD_SPMD_Task
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
!Plant Hydraulic parameters
            deallocate (vegwp_c  ) !vegetation water potential [mm]
            deallocate (gs0sun_c ) !working copy of sunlit stomata conductance
            deallocate (gs0sha_c ) !working copy of shalit stomata conductance
!end plant hydraulic parameters
!Ozone Stress Variables
            deallocate (o3coefv_sun_c) !Ozone stress factor for photosynthesis on sunlit leaf
            deallocate (o3coefv_sha_c) !Ozone stress factor for photosynthesis on shaded leaf
            deallocate (o3coefg_sun_c) !Ozone stress factor for stomata on sunlit leaf
            deallocate (o3coefg_sha_c) !Ozone stress factor for stomata on shaded leaf
            deallocate (lai_old_c    ) !lai in last time step
            deallocate (o3uptakesun_c) !Ozone does, sunlit leaf (mmol O3/m^2)
            deallocate (o3uptakesha_c) !Ozone does, shaded leaf (mmol O3/m^2)
!End Ozone Stress Variables
         ENDIF
      ENDIF

   END SUBROUTINE deallocate_PCTimeVariables

#ifdef CoLMDEBUG
   SUBROUTINE check_PCTimeVariables

      use MOD_CoLMDebug
      use MOD_Namelist, only: DEF_USE_PLANTHYDRAULICS, DEF_USE_OZONESTRESS
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
      if(DEF_USE_PLANTHYDRAULICS)then
         call check_vector_data ('vegwp_c  ', vegwp_c  )      !
         call check_vector_data ('gs0sun_c ', gs0sun_c )      !
         call check_vector_data ('gs0sha_c ', gs0sha_c )      !
      end if
      IF(DEF_USE_OZONESTRESS)THEN
         call check_vector_data ('o3coefv_sun_c', o3coefv_sun_c) !Ozone stress factor for photosynthesis on sunlit leaf
         call check_vector_data ('o3coefv_sha_c', o3coefv_sha_c) !Ozone stress factor for photosynthesis on shaded leaf
         call check_vector_data ('o3coefg_sun_c', o3coefg_sun_c) !Ozone stress factor for stomata on sunlit leaf
         call check_vector_data ('o3coefg_sha_c', o3coefg_sha_c) !Ozone stress factor for stomata on shaded leaf
         call check_vector_data ('lai_old_c    ', lai_old_c    ) !lai in last time step
         call check_vector_data ('o3uptakesun_c', o3uptakesun_c) !Ozone does, sunlit leaf (mmol O3/m^2)
         call check_vector_data ('o3uptakesha_c', o3uptakesha_c) !Ozone does, shaded leaf (mmol O3/m^2)
      END IF

   END SUBROUTINE check_PCTimeVariables
#endif

END MODULE MOD_Vars_PCTimeVariables
#endif



MODULE MOD_Vars_TimeVariables
! -------------------------------
! Created by Yongjiu Dai, 03/2014
! -------------------------------

  use MOD_Precision
  use MOD_TimeManager
#ifdef LULC_IGBP_PFT
  USE MOD_Vars_PFTimeVariables
#endif
#ifdef LULC_IGBP_PC
  USE MOD_Vars_PCTimeVariables
#endif
#ifdef BGC
  USE MOD_BGC_Vars_TimeVariables
#endif
#ifdef LATERAL_FLOW
  USE MOD_Hydro_Vars_TimeVariables
#endif
#ifdef URBAN_MODEL
  USE MOD_Urban_Vars_TimeVariables
#endif

  IMPLICIT NONE
  SAVE
! -----------------------------------------------------------------
! Time-varying state variables which reaquired by restart run
     real(r8), allocatable :: z_sno      (:,:) ! node depth [m]
     real(r8), allocatable :: dz_sno     (:,:) ! interface depth [m]
     real(r8), allocatable :: t_soisno   (:,:) ! soil temperature [K]
     real(r8), allocatable :: wliq_soisno(:,:) ! liquid water in layers [kg/m2]
     real(r8), allocatable :: wice_soisno(:,:) ! ice lens in layers [kg/m2]
     real(r8), allocatable :: h2osoi     (:,:) ! volumetric soil water in layers [m3/m3]
     real(r8), allocatable :: cvsoil     (:,:) !heat capacity [J/(m2 K)]
     real(r8), allocatable :: smp        (:,:) ! soil matrix potential [mm]
     real(r8), allocatable :: hk         (:,:) ! hydraulic conductivity [mm h2o/s]
     real(r8), allocatable :: rootr(:,:)       ! water exchange between soil and root. Positive: soil->root [?]
!Plant Hydraulic variables
     real(r8), allocatable :: vegwp(:,:)       ! vegetation water potential [mm]
     real(r8), allocatable :: gs0sun       (:) ! working copy of sunlit stomata conductance
     real(r8), allocatable :: gs0sha       (:) ! working copy of shalit stomata conductance
!end plant hydraulic variables
!Ozone stress variables
     real(r8), allocatable :: o3coefv_sun  (:) ! Ozone stress factor for photosynthesis on sunlit leaf
     real(r8), allocatable :: o3coefv_sha  (:) ! Ozone stress factor for photosynthesis on shaded leaf
     real(r8), allocatable :: o3coefg_sun  (:) ! Ozone stress factor for stomata on sunlit leaf
     real(r8), allocatable :: o3coefg_sha  (:) ! Ozone stress factor for stomata on shaded leaf
     real(r8), allocatable :: lai_old      (:) ! lai in last time step
     real(r8), allocatable :: o3uptakesun  (:) ! Ozone does, sunlit leaf (mmol O3/m^2)
     real(r8), allocatable :: o3uptakesha  (:) ! Ozone does, shaded leaf (mmol O3/m^2)
!End ozone stress variables
     real(r8), allocatable :: rstfacsun_out(:) ! factor of soil water stress on sunlit leaf
     real(r8), allocatable :: rstfacsha_out(:) ! factor of soil water stress on shaded leaf
     real(r8), allocatable :: gssun_out    (:) ! stomata conductance on sunlit leaf
     real(r8), allocatable :: gssha_out    (:) ! stomata conductance on shaded leaf
     real(r8), allocatable :: t_grnd       (:) ! ground surface temperature [K]

     real(r8), allocatable :: assimsun_out         (:) !1
     real(r8), allocatable :: assimsha_out         (:) !1
     real(r8), allocatable :: etrsun_out           (:) !1
     real(r8), allocatable :: etrsha_out           (:) !1

     real(r8), allocatable :: tleaf        (:) ! leaf temperature [K]
     real(r8), allocatable :: ldew         (:) ! depth of water on foliage [mm]
     real(r8), allocatable :: ldew_rain    (:) ! depth of rain on foliage [mm]
     real(r8), allocatable :: ldew_snow    (:) ! depth of rain on foliage [mm]
     real(r8), allocatable :: sag          (:) ! non dimensional snow age [-]
     real(r8), allocatable :: scv          (:) ! snow cover, water equivalent [mm]
     real(r8), allocatable :: snowdp       (:) ! snow depth [meter]
     real(r8), allocatable :: fveg         (:) ! fraction of vegetation cover
     real(r8), allocatable :: fsno         (:) ! fraction of snow cover on ground
     real(r8), allocatable :: sigf         (:) ! fraction of veg cover, excluding snow-covered veg [-]
     real(r8), allocatable :: green        (:) ! leaf greenness
     real(r8), allocatable :: tlai         (:) ! leaf area index
     real(r8), allocatable :: lai          (:) ! leaf area index
     real(r8), allocatable :: laisun       (:) ! leaf area index for sunlit leaf
     real(r8), allocatable :: laisha       (:) ! leaf area index for shaded leaf
     real(r8), allocatable :: tsai         (:) ! stem area index
     real(r8), allocatable :: sai          (:) ! stem area index
     real(r8), allocatable :: coszen       (:) ! cosine of solar zenith angle
     real(r8), allocatable :: alb      (:,:,:) ! averaged albedo [-]
     real(r8), allocatable :: ssun     (:,:,:) ! sunlit canopy absorption for solar radiation (0-1)
     real(r8), allocatable :: ssha     (:,:,:) ! shaded canopy absorption for solar radiation (0-1)
     real(r8), allocatable :: thermk       (:) ! canopy gap fraction for tir radiation
     real(r8), allocatable :: extkb        (:) ! (k, g(mu)/mu) direct solar extinction coefficient
     real(r8), allocatable :: extkd        (:) ! diffuse and scattered diffuse PAR extinction coefficient
     real(r8), allocatable :: zwt          (:) ! the depth to water table [m]
     real(r8), allocatable :: wa           (:) ! water storage in aquifer [mm]
     real(r8), allocatable :: wat          (:) ! total water storage [mm]
     real(r8), allocatable :: wdsrf        (:) ! depth of surface water [mm]

     real(r8), allocatable :: t_lake     (:,:) ! lake layer teperature [K]
     real(r8), allocatable :: lake_icefrac(:,:)! lake mass fraction of lake layer that is frozen
     real(r8), allocatable :: savedtke1    (:) ! top level eddy conductivity (W/m K)

     REAL(r8), allocatable :: snw_rds    (:,:) ! effective grain radius (col,lyr) [microns, m-6]
     REAL(r8), allocatable :: mss_bcpho  (:,:) ! mass of hydrophobic BC in snow  (col,lyr) [kg]
     REAL(r8), allocatable :: mss_bcphi  (:,:) ! mass of hydrophillic BC in snow (col,lyr) [kg]
     REAL(r8), allocatable :: mss_ocpho  (:,:) ! mass of hydrophobic OC in snow  (col,lyr) [kg]
     REAL(r8), allocatable :: mss_ocphi  (:,:) ! mass of hydrophillic OC in snow (col,lyr) [kg]
     REAL(r8), allocatable :: mss_dst1   (:,:) ! mass of dust species 1 in snow  (col,lyr) [kg]
     REAL(r8), allocatable :: mss_dst2   (:,:) ! mass of dust species 2 in snow  (col,lyr) [kg]
     REAL(r8), allocatable :: mss_dst3   (:,:) ! mass of dust species 3 in snow  (col,lyr) [kg]
     REAL(r8), allocatable :: mss_dst4   (:,:) ! mass of dust species 4 in snow  (col,lyr) [kg]
     REAL(r8), allocatable :: ssno   (:,:,:,:) ! snow layer absorption [-]

     real(r8), allocatable :: trad         (:) ! radiative temperature of surface [K]
     real(r8), allocatable :: tref         (:) ! 2 m height air temperature [kelvin]
     real(r8), allocatable :: qref         (:) ! 2 m height air specific humidity
     real(r8), allocatable :: rst          (:) ! canopy stomatal resistance (s/m)
     real(r8), allocatable :: emis         (:) ! averaged bulk surface emissivity
     real(r8), allocatable :: z0m          (:) ! effective roughness [m]
     real(r8), allocatable :: displa       (:) ! zero displacement height [m]
     real(r8), allocatable :: zol          (:) ! dimensionless height (z/L) used in Monin-Obukhov theory
     real(r8), allocatable :: rib          (:) ! bulk Richardson number in surface layer
     real(r8), allocatable :: ustar        (:) ! u* in similarity theory [m/s]
     real(r8), allocatable :: qstar        (:) ! q* in similarity theory [kg/kg]
     real(r8), allocatable :: tstar        (:) ! t* in similarity theory [K]
     real(r8), allocatable :: fm           (:) ! integral of profile function for momentum
     real(r8), allocatable :: fh           (:) ! integral of profile function for heat
     real(r8), allocatable :: fq           (:) ! integral of profile function for moisture

     ! PUBLIC MEMBER FUNCTIONS:
     public :: allocate_TimeVariables
     public :: deallocate_TimeVariables
     public :: READ_TimeVariables
     public :: WRITE_TimeVariables
#ifdef CoLMDEBUG
     public :: check_TimeVariables
#endif


!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_TimeVariables
! --------------------------------------------------------------------
! Allocates memory for CoLM 1d [numpatch] variables
! ------------------------------------------------------

     use MOD_Precision
     USE MOD_Vars_Global
     use MOD_SPMD_Task
     use MOD_LandPatch, only: numpatch
     IMPLICIT NONE


     if (p_is_worker) then

        if (numpatch > 0) then

           allocate (z_sno      (maxsnl+1:0,      numpatch))
           allocate (dz_sno     (maxsnl+1:0,      numpatch))
           allocate (t_soisno   (maxsnl+1:nl_soil,numpatch))
           allocate (wliq_soisno(maxsnl+1:nl_soil,numpatch))
           allocate (wice_soisno(maxsnl+1:nl_soil,numpatch))
           allocate (smp               (1:nl_soil,numpatch))
           allocate (hk                (1:nl_soil,numpatch))
           allocate (h2osoi            (1:nl_soil,numpatch))
           allocate (cvsoil            (1:nl_soil,numpatch))
           allocate (rootr             (1:nl_soil,numpatch))
!Plant Hydraulic variables
           allocate (vegwp             (1:nvegwcs,numpatch))
           allocate (gs0sun                      (numpatch))
           allocate (gs0sha                      (numpatch))
!end plant hydraulic variables
!Ozone Stress variables
           allocate (o3coefv_sun                 (numpatch)) ! Ozone stress factor for photosynthesis on sunlit leaf
           allocate (o3coefv_sha                 (numpatch)) ! Ozone stress factor for photosynthesis on shaded leaf
           allocate (o3coefg_sun                 (numpatch)) ! Ozone stress factor for stomata on sunlit leaf
           allocate (o3coefg_sha                 (numpatch)) ! Ozone stress factor for stomata on shaded leaf
           allocate (lai_old                     (numpatch)) ! lai in last time step
           allocate (o3uptakesun                 (numpatch)) ! Ozone does, sunlit leaf (mmol O3/m^2)
           allocate (o3uptakesha                 (numpatch)) ! Ozone does, shaded leaf (mmol O3/m^2)
!End ozone stress variables
           allocate (rstfacsun_out               (numpatch))
           allocate (rstfacsha_out               (numpatch))
           allocate (gssun_out                   (numpatch))
           allocate (gssha_out                   (numpatch))
           allocate (assimsun_out                (numpatch))
           allocate (assimsha_out                (numpatch))
           allocate (etrsun_out                  (numpatch))
           allocate (etrsha_out                  (numpatch))

           allocate (t_grnd                      (numpatch))
           allocate (tleaf                       (numpatch))
           allocate (ldew                        (numpatch))
           allocate (ldew_rain                   (numpatch))
           allocate (ldew_snow                   (numpatch))
           allocate (sag                         (numpatch))
           allocate (scv                         (numpatch))
           allocate (snowdp                      (numpatch))
           allocate (fveg                        (numpatch))
           allocate (fsno                        (numpatch))
           allocate (sigf                        (numpatch))
           allocate (green                       (numpatch))
           allocate (tlai                        (numpatch))
           allocate (lai                         (numpatch))
           allocate (laisun                      (numpatch))
           allocate (laisha                      (numpatch))
           allocate (tsai                        (numpatch))
           allocate (sai                         (numpatch))
           allocate (coszen                      (numpatch))
           allocate (alb                     (2,2,numpatch))
           allocate (ssun                    (2,2,numpatch))
           allocate (ssha                    (2,2,numpatch))
           allocate (thermk                      (numpatch))
           allocate (extkb                       (numpatch))
           allocate (extkd                       (numpatch))
           allocate (zwt                         (numpatch))
           allocate (wa                          (numpatch))
           allocate (wat                         (numpatch))
           allocate (wdsrf                       (numpatch))

           allocate (t_lake              (nl_lake,numpatch))!new lake scheme
           allocate (lake_icefrac        (nl_lake,numpatch))!new lake scheme
           allocate (savedtke1                   (numpatch))!new lake scheme

           allocate (snw_rds          (maxsnl+1:0,numpatch))
           allocate (mss_bcpho        (maxsnl+1:0,numpatch))
           allocate (mss_bcphi        (maxsnl+1:0,numpatch))
           allocate (mss_ocpho        (maxsnl+1:0,numpatch))
           allocate (mss_ocphi        (maxsnl+1:0,numpatch))
           allocate (mss_dst1         (maxsnl+1:0,numpatch))
           allocate (mss_dst2         (maxsnl+1:0,numpatch))
           allocate (mss_dst3         (maxsnl+1:0,numpatch))
           allocate (mss_dst4         (maxsnl+1:0,numpatch))
           allocate (ssno         (2,2,maxsnl+1:1,numpatch))

           allocate (trad                        (numpatch))
           allocate (tref                        (numpatch))
           allocate (qref                        (numpatch))
           allocate (rst                         (numpatch))
           allocate (emis                        (numpatch))
           allocate (z0m                         (numpatch))
           allocate (displa                      (numpatch))
           allocate (zol                         (numpatch))
           allocate (rib                         (numpatch))
           allocate (ustar                       (numpatch))
           allocate (qstar                       (numpatch))
           allocate (tstar                       (numpatch))
           allocate (fm                          (numpatch))
           allocate (fh                          (numpatch))
           allocate (fq                          (numpatch))

        end if
  end if

#ifdef LULC_IGBP_PFT
     CALL allocate_PFTimeVariables
#endif

#ifdef LULC_IGBP_PC
     CALL allocate_PCTimeVariables
#endif

#ifdef BGC
     CALL allocate_BGCTimeVariables
#endif

#ifdef LATERAL_FLOW
     CALL allocate_HydroTimeVariables
#endif

#ifdef URBAN_MODEL
     CALL allocate_UrbanTimeVariables
#endif

  END SUBROUTINE allocate_TimeVariables



  SUBROUTINE deallocate_TimeVariables ()

     use MOD_SPMD_Task
     use MOD_LandPatch, only: numpatch
     implicit none

     ! --------------------------------------------------
     ! Deallocates memory for CoLM 1d [numpatch] variables
     ! --------------------------------------------------

     if (p_is_worker) then

        if (numpatch > 0) then

           deallocate (z_sno                  )
           deallocate (dz_sno                 )
           deallocate (t_soisno               )
           deallocate (wliq_soisno            )
           deallocate (wice_soisno            )
           deallocate (smp                    )
           deallocate (hk                     )
           deallocate (h2osoi                 )
           deallocate (cvsoil                 )
           deallocate (rootr                  )
!Plant Hydraulic variables
           deallocate (vegwp                  )
           deallocate (gs0sun                 )
           deallocate (gs0sha                 )
!End plant hydraulic variables
!Ozone stress variables
           deallocate (o3coefv_sun            ) ! Ozone stress factor for photosynthesis on sunlit leaf
           deallocate (o3coefv_sha            ) ! Ozone stress factor for photosynthesis on shaded leaf
           deallocate (o3coefg_sun            ) ! Ozone stress factor for stomata on sunlit leaf
           deallocate (o3coefg_sha            ) ! Ozone stress factor for stomata on shaded leaf
           deallocate (lai_old                ) ! lai in last time step
           deallocate (o3uptakesun            ) ! Ozone does, sunlit leaf (mmol O3/m^2)
           deallocate (o3uptakesha            ) ! Ozone does, shaded leaf (mmol O3/m^2)
!End Ozone stress variables
           deallocate (rstfacsun_out          )
           deallocate (rstfacsha_out          )
           deallocate (gssun_out              )
           deallocate (gssha_out              )
           deallocate ( assimsun_out          )
           deallocate ( assimsha_out          )
           deallocate ( etrsun_out            )
           deallocate ( etrsha_out            )

           deallocate (t_grnd                 )
           deallocate (tleaf                  )
           deallocate (ldew                   )
           deallocate (ldew_rain              )
           deallocate (ldew_snow              )
           deallocate (sag                    )
           deallocate (scv                    )
           deallocate (snowdp                 )
           deallocate (fveg                   )
           deallocate (fsno                   )
           deallocate (sigf                   )
           deallocate (green                  )
           deallocate (tlai                   )
           deallocate (lai                    )
           deallocate (laisun                 )
           deallocate (laisha                 )
           deallocate (tsai                   )
           deallocate (sai                    )
           deallocate (coszen                 )
           deallocate (alb                    )
           deallocate (ssun                   )
           deallocate (ssha                   )
           deallocate (thermk                 )
           deallocate (extkb                  )
           deallocate (extkd                  )
           deallocate (zwt                    )
           deallocate (wa                     )
           deallocate (wat                    )
           deallocate (wdsrf                  )

           deallocate (t_lake                 ) ! new lake scheme
           deallocate (lake_icefrac           ) ! new lake scheme
           deallocate (savedtke1              ) ! new lake scheme

           deallocate (snw_rds                )
           deallocate (mss_bcpho              )
           deallocate (mss_bcphi              )
           deallocate (mss_ocpho              )
           deallocate (mss_ocphi              )
           deallocate (mss_dst1               )
           deallocate (mss_dst2               )
           deallocate (mss_dst3               )
           deallocate (mss_dst4               )
           deallocate (ssno                   )

           deallocate (trad                   )
           deallocate (tref                   )
           deallocate (qref                   )
           deallocate (rst                    )
           deallocate (emis                   )
           deallocate (z0m                    )
           deallocate (displa                 )
           deallocate (zol                    )
           deallocate (rib                    )
           deallocate (ustar                  )
           deallocate (qstar                  )
           deallocate (tstar                  )
           deallocate (fm                     )
           deallocate (fh                     )
           deallocate (fq                     )

        end if
     end if

#if (defined LULC_IGBP_PFT)
     CALL deallocate_PFTimeVariables
#endif

#if (defined LULC_IGBP_PC)
     CALL deallocate_PCTimeVariables
#endif

#if (defined BGC)
     CALL deallocate_BGCTimeVariables
#endif

#ifdef LATERAL_FLOW
     CALL deallocate_HydroTimeVariables
#endif

#if (defined URBAN_MODEL)
     CALL deallocate_UrbanTimeVariables
#endif

  END SUBROUTINE deallocate_TimeVariables


  !---------------------------------------
  function save_to_restart (idate, deltim, itstamp, ptstamp) result(rwrite)

     use MOD_Namelist
     implicit none

     logical :: rwrite

     integer,  intent(in) :: idate(3)
     real(r8), intent(in) :: deltim
     type(timestamp), intent(in) :: itstamp, ptstamp


     ! added by yuan, 08/31/2014
     select case (trim(adjustl(DEF_WRST_FREQ)))
     case ('TIMESTEP')
        rwrite = .true.
     case ('HOURLY')
        rwrite = isendofhour (idate, deltim)
     case ('DAILY')
        rwrite = isendofday(idate, deltim)
     case ('MONTHLY')
        rwrite = isendofmonth(idate, deltim)
     case ('YEARLY')
        rwrite = isendofyear(idate, deltim)
     case default
        write(*,*) 'Warning: Please use one of TIMESTEP/HOURLY/DAILY/MONTHLY/YEARLY for restart frequency.'
     end select

     if (rwrite) then
        rwrite = (ptstamp < itstamp)
     end if

  end function save_to_restart

  !---------------------------------------
  SUBROUTINE WRITE_TimeVariables (idate, lc_year, site, dir_restart)

     !=======================================================================
     ! Original version: Yongjiu Dai, September 15, 1999, 03/2014
     !=======================================================================

     use MOD_Namelist, only : DEF_REST_COMPRESS_LEVEL, DEF_USE_PLANTHYDRAULICS, DEF_USE_OZONESTRESS
     USE MOD_LandPatch
     use MOD_NetCDFVector
     USE MOD_Vars_Global
     IMPLICIT NONE

     integer, INTENT(in) :: idate(3)
     INTEGER, intent(in) :: lc_year      !year of land cover TYPE data
     character(LEN=*), intent(in) :: site
     character(LEN=*), intent(in) :: dir_restart

     ! Local variables
     character(LEN=256) :: file_restart
     character(len=14)  :: cdate
     CHARACTER(len=256) :: cyear         !CHARACTER for lc_year
     integer :: compress

     compress = DEF_REST_COMPRESS_LEVEL

     ! land cover type year
     write(cyear,'(i4.4)') lc_year

     write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') idate(1), idate(2), idate(3)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_'//trim(cdate)//'_lc'//trim(cyear)//'.nc'


     call ncio_create_file_vector (file_restart, landpatch)
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'patch')

     CALL ncio_define_dimension_vector (file_restart, landpatch, 'snow',     -maxsnl       )
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'snowp1',   -maxsnl+1     )
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'soilsnow', nl_soil-maxsnl)
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'soil',     nl_soil)
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'lake',     nl_lake)

     if(DEF_USE_PLANTHYDRAULICS)then
        CALL ncio_define_dimension_vector (file_restart, landpatch, 'vegnodes', nvegwcs)
     end if

     CALL ncio_define_dimension_vector (file_restart, landpatch, 'band', 2)
     CALL ncio_define_dimension_vector (file_restart, landpatch, 'rtyp', 2)

     ! Time-varying state variables which reaquired by restart run
     call ncio_write_vector (file_restart, 'z_sno   '   , 'snow', -maxsnl, 'patch', landpatch, z_sno , compress)                 ! node depth [m]
     call ncio_write_vector (file_restart, 'dz_sno  '   , 'snow', -maxsnl, 'patch', landpatch, dz_sno, compress)                 ! interface depth [m]
     call ncio_write_vector (file_restart, 't_soisno'   , 'soilsnow', nl_soil-maxsnl, 'patch', landpatch, t_soisno   , compress) ! soil temperature [K]
     call ncio_write_vector (file_restart, 'wliq_soisno', 'soilsnow', nl_soil-maxsnl, 'patch', landpatch, wliq_soisno, compress) ! liquid water in layers [kg/m2]
     call ncio_write_vector (file_restart, 'wice_soisno', 'soilsnow', nl_soil-maxsnl, 'patch', landpatch, wice_soisno, compress) ! ice lens in layers [kg/m2]
     call ncio_write_vector (file_restart, 'smp',         'soil', nl_soil, 'patch', landpatch, smp, compress)                    ! soil matrix potential [mm]
     call ncio_write_vector (file_restart, 'hk',          'soil', nl_soil, 'patch', landpatch, hk, compress)                     ! hydraulic conductivity [mm h2o/s]
     IF(DEF_USE_PLANTHYDRAULICS)THEN
        call ncio_write_vector (file_restart, 'vegwp',   'vegnodes', nvegwcs, 'patch', landpatch, vegwp, compress)               ! vegetation water potential [mm]
        call ncio_write_vector (file_restart, 'gs0sun',    'patch', landpatch, gs0sun, compress)                                 ! working copy of sunlit stomata conductance
        call ncio_write_vector (file_restart, 'gs0sha',    'patch', landpatch, gs0sha, compress)                                 ! working copy of shalit stomata conductance
     ENDIF
     IF(DEF_USE_OZONESTRESS)THEN
        call ncio_write_vector (file_restart, 'lai_old    ', 'patch', landpatch, lai_old    , compress)
        call ncio_write_vector (file_restart, 'o3uptakesun', 'patch', landpatch, o3uptakesun, compress)
        call ncio_write_vector (file_restart, 'o3uptakesha', 'patch', landpatch, o3uptakesha, compress)
     ENDIF
     call ncio_write_vector (file_restart, 't_grnd  '   , 'patch', landpatch, t_grnd    , compress)                    ! ground surface temperature [K]
     call ncio_write_vector (file_restart, 'tleaf   '   , 'patch', landpatch, tleaf     , compress)                    ! leaf temperature [K]
     call ncio_write_vector (file_restart, 'ldew    '   , 'patch', landpatch, ldew      , compress)                    ! depth of water on foliage [mm]
     call ncio_write_vector (file_restart, 'ldew_rain'  , 'patch', landpatch, ldew_rain , compress)                    ! depth of water on foliage [mm]
     call ncio_write_vector (file_restart, 'ldew_snow'  , 'patch', landpatch, ldew_snow , compress)                    ! depth of water on foliage [mm]
     call ncio_write_vector (file_restart, 'sag     '   , 'patch', landpatch, sag       , compress)                    ! non dimensional snow age [-]
     call ncio_write_vector (file_restart, 'scv     '   , 'patch', landpatch, scv       , compress)                    ! snow cover, water equivalent [mm]
     call ncio_write_vector (file_restart, 'snowdp  '   , 'patch', landpatch, snowdp    , compress)                    ! snow depth [meter]
     call ncio_write_vector (file_restart, 'fveg    '   , 'patch', landpatch, fveg      , compress)                    ! fraction of vegetation cover
     call ncio_write_vector (file_restart, 'fsno    '   , 'patch', landpatch, fsno      , compress)                    ! fraction of snow cover on ground
     call ncio_write_vector (file_restart, 'sigf    '   , 'patch', landpatch, sigf      , compress)                    ! fraction of veg cover, excluding snow-covered veg [-]
     call ncio_write_vector (file_restart, 'green   '   , 'patch', landpatch, green     , compress)                    ! leaf greenness
     call ncio_write_vector (file_restart, 'lai     '   , 'patch', landpatch, lai       , compress)                    ! leaf area index
     call ncio_write_vector (file_restart, 'tlai    '   , 'patch', landpatch, tlai      , compress)                    ! leaf area index
     call ncio_write_vector (file_restart, 'sai     '   , 'patch', landpatch, sai       , compress)                    ! stem area index
     call ncio_write_vector (file_restart, 'tsai    '   , 'patch', landpatch, tsai      , compress)                    ! stem area index
     call ncio_write_vector (file_restart, 'coszen  '   , 'patch', landpatch, coszen    , compress)                    ! cosine of solar zenith angle
     call ncio_write_vector (file_restart, 'alb     '   , 'band', 2, 'rtyp', 2, 'patch', landpatch, alb , compress)    ! averaged albedo [-]
     call ncio_write_vector (file_restart, 'ssun    '   , 'band', 2, 'rtyp', 2, 'patch', landpatch, ssun, compress)    ! sunlit canopy absorption for solar radiation (0-1)
     call ncio_write_vector (file_restart, 'ssha    '   , 'band', 2, 'rtyp', 2, 'patch', landpatch, ssha, compress)    ! shaded canopy absorption for solar radiation (0-1)
     call ncio_write_vector (file_restart, 'thermk  '   , 'patch', landpatch, thermk    , compress)                    ! canopy gap fraction for tir radiation
     call ncio_write_vector (file_restart, 'extkb   '   , 'patch', landpatch, extkb     , compress)                    ! (k, g(mu)/mu) direct solar extinction coefficient
     call ncio_write_vector (file_restart, 'extkd   '   , 'patch', landpatch, extkd     , compress)                    ! diffuse and scattered diffuse PAR extinction coefficient
     call ncio_write_vector (file_restart, 'zwt     '   , 'patch', landpatch, zwt       , compress)                    ! the depth to water table [m]
     call ncio_write_vector (file_restart, 'wa      '   , 'patch', landpatch, wa        , compress)                    ! water storage in aquifer [mm]
     call ncio_write_vector (file_restart, 'wdsrf   '   , 'patch', landpatch, wdsrf     , compress)                    ! depth of surface water [mm]

     call ncio_write_vector (file_restart, 't_lake  '   , 'lake', nl_lake, 'patch', landpatch, t_lake      , compress) !
     call ncio_write_vector (file_restart, 'lake_icefrc', 'lake', nl_lake, 'patch', landpatch, lake_icefrac, compress) !
     call ncio_write_vector (file_restart, 'savedtke1  ', 'patch', landpatch, savedtke1   , compress)                  !
     call ncio_write_vector (file_restart, 'snw_rds  ', 'snow', -maxsnl, 'patch', landpatch, snw_rds  , compress)
     call ncio_write_vector (file_restart, 'mss_bcpho', 'snow', -maxsnl, 'patch', landpatch, mss_bcpho, compress)
     call ncio_write_vector (file_restart, 'mss_bcphi', 'snow', -maxsnl, 'patch', landpatch, mss_bcphi, compress)
     call ncio_write_vector (file_restart, 'mss_ocpho', 'snow', -maxsnl, 'patch', landpatch, mss_ocpho, compress)
     call ncio_write_vector (file_restart, 'mss_ocphi', 'snow', -maxsnl, 'patch', landpatch, mss_ocphi, compress)
     call ncio_write_vector (file_restart, 'mss_dst1 ', 'snow', -maxsnl, 'patch', landpatch, mss_dst1 , compress)
     call ncio_write_vector (file_restart, 'mss_dst2 ', 'snow', -maxsnl, 'patch', landpatch, mss_dst2 , compress)
     call ncio_write_vector (file_restart, 'mss_dst3 ', 'snow', -maxsnl, 'patch', landpatch, mss_dst3 , compress)
     call ncio_write_vector (file_restart, 'mss_dst4 ', 'snow', -maxsnl, 'patch', landpatch, mss_dst4 , compress)
     call ncio_write_vector (file_restart, 'ssno', 'band', 2, 'rtyp', 2, 'snowp1', -maxsnl+1, 'patch', landpatch, ssno, compress)

     ! Additional va_vectorriables required by reginal model (such as WRF ) RSM)
     call ncio_write_vector (file_restart, 'trad ', 'patch', landpatch, trad , compress) ! radiative temperature of surface [K]
     call ncio_write_vector (file_restart, 'tref ', 'patch', landpatch, tref , compress) ! 2 m height air temperature [kelvin]
     call ncio_write_vector (file_restart, 'qref ', 'patch', landpatch, qref , compress) ! 2 m height air specific humidity
     call ncio_write_vector (file_restart, 'rst  ', 'patch', landpatch, rst  , compress) ! canopy stomatal resistance (s/m)
     call ncio_write_vector (file_restart, 'emis ', 'patch', landpatch, emis , compress) ! averaged bulk surface emissivity
     call ncio_write_vector (file_restart, 'z0m  ', 'patch', landpatch, z0m  , compress) ! effective roughness [m]
     call ncio_write_vector (file_restart, 'zol  ', 'patch', landpatch, zol  , compress) ! dimensionless height (z/L) used in Monin-Obukhov theory
     call ncio_write_vector (file_restart, 'rib  ', 'patch', landpatch, rib  , compress) ! bulk Richardson number in surface layer
     call ncio_write_vector (file_restart, 'ustar', 'patch', landpatch, ustar, compress) ! u* in similarity theory [m/s]
     call ncio_write_vector (file_restart, 'qstar', 'patch', landpatch, qstar, compress) ! q* in similarity theory [kg/kg]
     call ncio_write_vector (file_restart, 'tstar', 'patch', landpatch, tstar, compress) ! t* in similarity theory [K]
     call ncio_write_vector (file_restart, 'fm   ', 'patch', landpatch, fm   , compress) ! integral of profile function for momentum
     call ncio_write_vector (file_restart, 'fh   ', 'patch', landpatch, fh   , compress) ! integral of profile function for heat
     call ncio_write_vector (file_restart, 'fq   ', 'patch', landpatch, fq   , compress) ! integral of profile function for moisture

#if (defined LULC_IGBP_PFT)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_pft_'//trim(cdate)//'_lc'//trim(cyear)//'.nc'
     CALL WRITE_PFTimeVariables (file_restart)
#endif

#if (defined LULC_IGBP_PC)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_pc_'//trim(cdate)//'_lc'//trim(cyear)//'.nc'
     CALL WRITE_PCTimeVariables (file_restart)
#endif

#if (defined BGC)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_bgc_'//trim(cdate)//'_lc'//trim(cyear)//'.nc'
     CALL WRITE_BGCTimeVariables (file_restart)
#endif

#if (defined LATERAL_FLOW)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_basin_'//trim(cdate)//'_lc'//trim(cyear)//'.nc'
     CALL WRITE_HydroTimeVariables (file_restart)
#endif

#if (defined URBAN_MODEL)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_urban_'//trim(cdate)//'_lc'//trim(cyear)//'.nc'
     CALL WRITE_UrbanTimeVariables (file_restart)
#endif
  end subroutine WRITE_TimeVariables

  !---------------------------------------
  SUBROUTINE READ_TimeVariables (idate, lc_year, site, dir_restart)

     !=======================================================================
     ! Original version: Yongjiu Dai, September 15, 1999, 03/2014
     !=======================================================================

     use MOD_Namelist
     use MOD_SPMD_Task
     use MOD_NetCDFVector
#ifdef CoLMDEBUG
     USE MOD_CoLMDebug
#endif
     USE MOD_LandPatch
     USE MOD_Vars_Global

     IMPLICIT NONE

     integer, INTENT(in) :: idate(3)
     INTEGER, intent(in) :: lc_year      !year of land cover type data
     character(LEN=*), intent(in) :: site
     character(LEN=*), intent(in) :: dir_restart

     ! Local variables
     character(LEN=256) :: file_restart
     character(len=14)  :: cdate, cyear

#ifdef USEMPI
     call mpi_barrier (p_comm_glb, p_err)
#endif

     if (p_is_master) then
        write(*,'(/,A26)') 'Loading Time Variables ...'
     end if

     ! land cover type year
     write(cyear,'(i4.4)') lc_year

     write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') idate(1), idate(2), idate(3)
     file_restart = trim(dir_restart) // '/' // trim(site) //'_restart_'//trim(cdate)//'_lc'//trim(cyear)//'.nc'

     ! Time-varying state variables which reaquired by restart run
     call ncio_read_vector (file_restart, 'z_sno   '   , -maxsnl, landpatch, z_sno )             ! node depth [m]
     call ncio_read_vector (file_restart, 'dz_sno  '   , -maxsnl, landpatch, dz_sno)             ! interface depth [m]
     call ncio_read_vector (file_restart, 't_soisno'   , nl_soil-maxsnl, landpatch, t_soisno   ) ! soil temperature [K]
     call ncio_read_vector (file_restart, 'wliq_soisno', nl_soil-maxsnl, landpatch, wliq_soisno) ! liquid water in layers [kg/m2]
     call ncio_read_vector (file_restart, 'wice_soisno', nl_soil-maxsnl, landpatch, wice_soisno) ! ice lens in layers [kg/m2]
     call ncio_read_vector (file_restart, 'smp',         nl_soil,        landpatch, smp        ) ! soil matrix potential [mm]
     call ncio_read_vector (file_restart, 'hk',          nl_soil,        landpatch, hk         ) ! hydraulic conductivity [mm h2o/s]
     if(DEF_USE_PLANTHYDRAULICS)then
        call ncio_read_vector (file_restart, 'vegwp',       nvegwcs,        landpatch, vegwp      ) ! vegetation water potential [mm]
        call ncio_read_vector (file_restart, 'gs0sun  ',    landpatch, gs0sun     ) ! working copy of sunlit stomata conductance
        call ncio_read_vector (file_restart, 'gs0sha  ',    landpatch, gs0sha     ) ! working copy of shalit stomata conductance
     end if
     IF(DEF_USE_OZONESTRESS)THEN
        call ncio_read_vector (file_restart, 'lai_old    ', landpatch, lai_old    )
        call ncio_read_vector (file_restart, 'o3uptakesun', landpatch, o3uptakesun)
        call ncio_read_vector (file_restart, 'o3uptakesha', landpatch, o3uptakesha)
     ENDIF
     call ncio_read_vector (file_restart, 't_grnd  '   , landpatch, t_grnd     ) ! ground surface temperature [K]
     call ncio_read_vector (file_restart, 'tleaf   '   , landpatch, tleaf      ) ! leaf temperature [K]
     call ncio_read_vector (file_restart, 'ldew    '   , landpatch, ldew       ) ! depth of water on foliage [mm]
     call ncio_read_vector (file_restart, 'ldew_rain'  , landpatch, ldew_rain  ) ! depth of rain on foliage [mm]
     call ncio_read_vector (file_restart, 'ldew_snow'  , landpatch, ldew_snow  ) ! depth of snow on foliage [mm]
     call ncio_read_vector (file_restart, 'sag     '   , landpatch, sag        ) ! non dimensional snow age [-]
     call ncio_read_vector (file_restart, 'scv     '   , landpatch, scv        ) ! snow cover, water equivalent [mm]
     call ncio_read_vector (file_restart, 'snowdp  '   , landpatch, snowdp     ) ! snow depth [meter]
     call ncio_read_vector (file_restart, 'fveg    '   , landpatch, fveg       ) ! fraction of vegetation cover
     call ncio_read_vector (file_restart, 'fsno    '   , landpatch, fsno       ) ! fraction of snow cover on ground
     call ncio_read_vector (file_restart, 'sigf    '   , landpatch, sigf       ) ! fraction of veg cover, excluding snow-covered veg [-]
     call ncio_read_vector (file_restart, 'green   '   , landpatch, green      ) ! leaf greenness
     call ncio_read_vector (file_restart, 'lai     '   , landpatch, lai        ) ! leaf area index
     call ncio_read_vector (file_restart, 'tlai    '   , landpatch, tlai       ) ! leaf area index
     call ncio_read_vector (file_restart, 'sai     '   , landpatch, sai        ) ! stem area index
     call ncio_read_vector (file_restart, 'tsai    '   , landpatch, tsai       ) ! stem area index
     call ncio_read_vector (file_restart, 'coszen  '   , landpatch, coszen     ) ! cosine of solar zenith angle
     call ncio_read_vector (file_restart, 'alb     '   , 2, 2, landpatch, alb  ) ! averaged albedo [-]
     call ncio_read_vector (file_restart, 'ssun    '   , 2, 2, landpatch, ssun ) ! sunlit canopy absorption for solar radiation (0-1)
     call ncio_read_vector (file_restart, 'ssha    '   , 2, 2, landpatch, ssha ) ! shaded canopy absorption for solar radiation (0-1)
     call ncio_read_vector (file_restart, 'thermk  '   , landpatch, thermk     ) ! canopy gap fraction for tir radiation
     call ncio_read_vector (file_restart, 'extkb   '   , landpatch, extkb      ) ! (k, g(mu)/mu) direct solar extinction coefficient
     call ncio_read_vector (file_restart, 'extkd   '   , landpatch, extkd      ) ! diffuse and scattered diffuse PAR extinction coefficient
     call ncio_read_vector (file_restart, 'zwt     '   , landpatch, zwt        ) ! the depth to water table [m]
     call ncio_read_vector (file_restart, 'wa      '   , landpatch, wa         ) ! water storage in aquifer [mm]
     call ncio_read_vector (file_restart, 'wdsrf   '   , landpatch, wdsrf      ) ! depth of surface water [mm]

     call ncio_read_vector (file_restart, 't_lake  '   , nl_lake, landpatch, t_lake      ) !
     call ncio_read_vector (file_restart, 'lake_icefrc', nl_lake, landpatch, lake_icefrac) !
     call ncio_read_vector (file_restart, 'savedtke1', landpatch, savedtke1) !

     call ncio_read_vector (file_restart, 'snw_rds  ', -maxsnl, landpatch, snw_rds  ) !
     call ncio_read_vector (file_restart, 'mss_bcpho', -maxsnl, landpatch, mss_bcpho) !
     call ncio_read_vector (file_restart, 'mss_bcphi', -maxsnl, landpatch, mss_bcphi) !
     call ncio_read_vector (file_restart, 'mss_ocpho', -maxsnl, landpatch, mss_ocpho) !
     call ncio_read_vector (file_restart, 'mss_ocphi', -maxsnl, landpatch, mss_ocphi) !
     call ncio_read_vector (file_restart, 'mss_dst1 ', -maxsnl, landpatch, mss_dst1 ) !
     call ncio_read_vector (file_restart, 'mss_dst2 ', -maxsnl, landpatch, mss_dst2 ) !
     call ncio_read_vector (file_restart, 'mss_dst3 ', -maxsnl, landpatch, mss_dst3 ) !
     call ncio_read_vector (file_restart, 'mss_dst4 ', -maxsnl, landpatch, mss_dst4 ) !
     call ncio_read_vector (file_restart, 'ssno', 2,2, -maxsnl+1, landpatch, ssno) !

     ! Additional variables required by reginal model (such as WRF ) RSM)
     call ncio_read_vector (file_restart, 'trad ', landpatch, trad ) ! radiative temperature of surface [K]
     call ncio_read_vector (file_restart, 'tref ', landpatch, tref ) ! 2 m height air temperature [kelvin]
     call ncio_read_vector (file_restart, 'qref ', landpatch, qref ) ! 2 m height air specific humidity
     call ncio_read_vector (file_restart, 'rst  ', landpatch, rst  ) ! canopy stomatal resistance (s/m)
     call ncio_read_vector (file_restart, 'emis ', landpatch, emis ) ! averaged bulk surface emissivity
     call ncio_read_vector (file_restart, 'z0m  ', landpatch, z0m  ) ! effective roughness [m]
     call ncio_read_vector (file_restart, 'zol  ', landpatch, zol  ) ! dimensionless height (z/L) used in Monin-Obukhov theory
     call ncio_read_vector (file_restart, 'rib  ', landpatch, rib  ) ! bulk Richardson number in surface layer
     call ncio_read_vector (file_restart, 'ustar', landpatch, ustar) ! u* in similarity theory [m/s]
     call ncio_read_vector (file_restart, 'qstar', landpatch, qstar) ! q* in similarity theory [kg/kg]
     call ncio_read_vector (file_restart, 'tstar', landpatch, tstar) ! t* in similarity theory [K]
     call ncio_read_vector (file_restart, 'fm   ', landpatch, fm   ) ! integral of profile function for momentum
     call ncio_read_vector (file_restart, 'fh   ', landpatch, fh   ) ! integral of profile function for heat
     call ncio_read_vector (file_restart, 'fq   ', landpatch, fq   ) ! integral of profile function for moisture

#if (defined LULC_IGBP_PFT)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_pft_'//trim(cdate)//'_lc'//trim(cyear)//'.nc'
     CALL READ_PFTimeVariables (file_restart)
#endif

#if (defined LULC_IGBP_PC)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_pc_'//trim(cdate)//'_lc'//trim(cyear)//'.nc'
     CALL READ_PCTimeVariables (file_restart)
#endif

#if (defined BGC)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_bgc_'//trim(cdate)//'_lc'//trim(cyear)//'.nc'
     CALL READ_BGCTimeVariables (file_restart)
#endif

#if (defined LATERAL_FLOW)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_basin_'//trim(cdate)//'_lc'//trim(cyear)//'.nc'
     CALL READ_HydroTimeVariables (file_restart)
#endif

#if (defined URBAN_MODEL)
     file_restart = trim(dir_restart)// '/' // trim(site) //'_restart_urban_'//trim(cdate)//'_lc'//trim(cyear)//'.nc'
     CALL READ_UrbanTimeVariables (file_restart)
#endif

#ifdef CoLMDEBUG
     call check_TimeVariables
#endif

     if (p_is_master) then
        write(*,*) 'Loading Time Variables done.'
     end if

  end subroutine READ_TimeVariables

  !---------------------------------------
#ifdef CoLMDEBUG
  SUBROUTINE check_TimeVariables ()

     use MOD_SPMD_Task
     use MOD_CoLMDebug
     use MOD_Namelist, only: DEF_USE_PLANTHYDRAULICS, DEF_USE_OZONESTRESS

     IMPLICIT NONE

#ifdef USEMPI
     call mpi_barrier (p_comm_glb, p_err)
#endif
     if (p_is_master) then
        write(*,'(/,A27)') 'Checking Time Variables ...'
     end if

     call check_vector_data ('z_sno       ', z_sno )      ! node depth [m]
     call check_vector_data ('dz_sno      ', dz_sno)      ! interface depth [m]
     call check_vector_data ('t_soisno    ', t_soisno   ) ! soil temperature [K]
     call check_vector_data ('wliq_soisno ', wliq_soisno) ! liquid water in layers [kg/m2]
     call check_vector_data ('wice_soisno ', wice_soisno) ! ice lens in layers [kg/m2]
     call check_vector_data ('smp         ', smp        ) ! soil matrix potential [mm]
     call check_vector_data ('hk          ', hk         ) ! hydraulic conductivity [mm h2o/s]
     if(DEF_USE_PLANTHYDRAULICS)then
        call check_vector_data ('vegwp       ', vegwp      ) ! vegetation water potential [mm]
        call check_vector_data ('gs0sun      ', gs0sun     ) ! working copy of sunlit stomata conductance
        call check_vector_data ('gs0sha      ', gs0sha     ) ! working copy of shalit stomata conductance
     end if
     IF(DEF_USE_OZONESTRESS)THEN
        call check_vector_data ('o3coefv_sun', o3coefv_sun)
        call check_vector_data ('o3coefv_sha', o3coefv_sha)
        call check_vector_data ('o3coefg_sun', o3coefg_sun)
        call check_vector_data ('o3coefg_sha', o3coefg_sha)
        call check_vector_data ('lai_old    ', lai_old    )
        call check_vector_data ('o3uptakesun', o3uptakesun)
        call check_vector_data ('o3uptakesha', o3uptakesha)
     ENDIF
     call check_vector_data ('t_grnd      ', t_grnd     ) ! ground surface temperature [K]
     call check_vector_data ('tleaf       ', tleaf      ) ! leaf temperature [K]
     call check_vector_data ('ldew        ', ldew       ) ! depth of water on foliage [mm]
     call check_vector_data ('ldew_rain   ', ldew_rain  ) ! depth of rain on foliage [mm]
     call check_vector_data ('ldew_snow   ', ldew_snow  ) ! depth of snow on foliage [mm]
     call check_vector_data ('sag         ', sag        ) ! non dimensional snow age [-]
     call check_vector_data ('scv         ', scv        ) ! snow cover, water equivalent [mm]
     call check_vector_data ('snowdp      ', snowdp     ) ! snow depth [meter]
     call check_vector_data ('fveg        ', fveg       ) ! fraction of vegetation cover
     call check_vector_data ('fsno        ', fsno       ) ! fraction of snow cover on ground
     call check_vector_data ('sigf        ', sigf       ) ! fraction of veg cover, excluding snow-covered veg [-]
     call check_vector_data ('green       ', green      ) ! leaf greenness
     call check_vector_data ('lai         ', lai        ) ! leaf area index
     call check_vector_data ('tlai        ', tlai       ) ! leaf area index
     call check_vector_data ('sai         ', sai        ) ! stem area index
     call check_vector_data ('tsai        ', tsai       ) ! stem area index
     call check_vector_data ('coszen      ', coszen     ) ! cosine of solar zenith angle
     call check_vector_data ('alb         ', alb        ) ! averaged albedo [-]
     call check_vector_data ('ssun        ', ssun       ) ! sunlit canopy absorption for solar radiation (0-1)
     call check_vector_data ('ssha        ', ssha       ) ! shaded canopy absorption for solar radiation (0-1)
     call check_vector_data ('thermk      ', thermk     ) ! canopy gap fraction for tir radiation
     call check_vector_data ('extkb       ', extkb      ) ! (k, g(mu)/mu) direct solar extinction coefficient
     call check_vector_data ('extkd       ', extkd      ) ! diffuse and scattered diffuse PAR extinction coefficient
     call check_vector_data ('zwt         ', zwt        ) ! the depth to water table [m]
     call check_vector_data ('wa          ', wa         ) ! water storage in aquifer [mm]
     call check_vector_data ('wdsrf       ', wdsrf      ) ! depth of surface water [mm]

     call check_vector_data ('t_lake      ', t_lake      )!
     call check_vector_data ('lake_icefrc ', lake_icefrac)!
     call check_vector_data ('savedtke1   ', savedtke1   )!

#if (defined LULC_IGBP_PFT)
     CALL check_PFTimeVariables
#endif

#if (defined LULC_IGBP_PC)
     CALL check_PCTimeVariables
#endif

#if (defined BGC)
     CALL check_BGCTimeVariables
#endif

#ifdef USEMPI
     call mpi_barrier (p_comm_glb, p_err)
#endif

  end subroutine check_TimeVariables
#endif


END MODULE MOD_Vars_TimeVariables
! ---------- EOP ------------
