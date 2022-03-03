#include <define.h> 

#if (defined PFT_CLASSIFICATION)

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
  PUBLIC :: READ_PFTimeVars
  PUBLIC :: WRITE_PFTimeVars
#ifdef CLMDEBUG
  PUBLIC :: check_PFTimeVars
#endif

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE allocate_PFTimeVars ()
! ------------------------------------------------------
! Allocates memory for CLM 1d [numpft] variables
! ------------------------------------------------------
      USE precision
      USE spmd_task
      USE mod_landpft
      IMPLICIT NONE

      IF (p_is_worker) THEN
         IF (numpft > 0) THEN
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
         ENDIF
      ENDIF 

   END SUBROUTINE allocate_PFTimeVars
  
   SUBROUTINE READ_PFTimeVars (file_restart)

      use ncio_vector
      USE mod_landpft
      IMPLICIT NONE

      character(LEN=*), intent(in) :: file_restart

      call ncio_read_vector (file_restart, 'tleaf_p  ', landpft, tleaf_p    ) !  
      call ncio_read_vector (file_restart, 'ldew_p   ', landpft, ldew_p     ) !  
      call ncio_read_vector (file_restart, 'sigf_p   ', landpft, sigf_p     ) !  
      call ncio_read_vector (file_restart, 'tlai_p   ', landpft, tlai_p     ) !  
      call ncio_read_vector (file_restart, 'lai_p    ', landpft, lai_p      ) !  
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

   END SUBROUTINE READ_PFTimeVars 
   
   SUBROUTINE WRITE_PFTimeVars (file_restart)

     use mod_namelist, only : DEF_REST_COMPRESS_LEVEL 
     USE mod_landpft
     use ncio_vector
     IMPLICIT NONE

     character(LEN=*), intent(in) :: file_restart
     
     ! Local variables
     integer :: compress

     compress = DEF_REST_COMPRESS_LEVEL 

     call ncio_create_file_vector (file_restart, landpft)
     CALL ncio_define_pixelset_dimension (file_restart, landpft)
     CALL ncio_define_dimension_vector (file_restart, 'band',   2)
     CALL ncio_define_dimension_vector (file_restart, 'wetdry', 2)

     call ncio_write_vector (file_restart, 'tleaf_p  ', 'vector', landpft, tleaf_p  , compress) !  
     call ncio_write_vector (file_restart, 'ldew_p   ', 'vector', landpft, ldew_p   , compress) !  
     call ncio_write_vector (file_restart, 'sigf_p   ', 'vector', landpft, sigf_p   , compress) !  
     call ncio_write_vector (file_restart, 'tlai_p   ', 'vector', landpft, tlai_p   , compress) !  
     call ncio_write_vector (file_restart, 'lai_p    ', 'vector', landpft, lai_p    , compress) !  
     call ncio_write_vector (file_restart, 'tsai_p   ', 'vector', landpft, tsai_p   , compress) !  
     call ncio_write_vector (file_restart, 'sai_p    ', 'vector', landpft, sai_p    , compress) !  
     call ncio_write_vector (file_restart, 'ssun_p   ', 'band', 2, 'wetdry', 2, 'vector', landpft, ssun_p, compress) !  
     call ncio_write_vector (file_restart, 'ssha_p   ', 'band', 2, 'wetdry', 2, 'vector', landpft, ssha_p, compress) !  
     call ncio_write_vector (file_restart, 'thermk_p ', 'vector', landpft, thermk_p , compress) !  
     call ncio_write_vector (file_restart, 'extkb_p  ', 'vector', landpft, extkb_p  , compress) !  
     call ncio_write_vector (file_restart, 'extkd_p  ', 'vector', landpft, extkd_p  , compress) !  
     call ncio_write_vector (file_restart, 'tref_p   ', 'vector', landpft, tref_p   , compress) !  
     call ncio_write_vector (file_restart, 'qref_p   ', 'vector', landpft, qref_p   , compress) !  
     call ncio_write_vector (file_restart, 'rst_p    ', 'vector', landpft, rst_p    , compress) !  
     call ncio_write_vector (file_restart, 'z0m_p    ', 'vector', landpft, z0m_p    , compress) !  

   END SUBROUTINE WRITE_PFTimeVars

  
   SUBROUTINE deallocate_PFTimeVars
! --------------------------------------------------
! Deallocates memory for CLM 1d [numpft/numpc] variables
! --------------------------------------------------
      USE spmd_task
      USE mod_landpft

      IF (p_is_worker) THEN
         IF (numpft > 0) THEN
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
         ENDIF
      ENDIF 

   END SUBROUTINE deallocate_PFTimeVars

#ifdef CLMDEBUG
   SUBROUTINE check_PFTimeVars

      use mod_colm_debug
      IMPLICIT NONE

      call check_vector_data ('tleaf_p  ', tleaf_p  )      !  
      call check_vector_data ('ldew_p   ', ldew_p   )      !  
      call check_vector_data ('sigf_p   ', sigf_p   )      !  
      call check_vector_data ('tlai_p   ', tlai_p   )      !  
      call check_vector_data ('lai_p    ', lai_p    )      !  
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

   END SUBROUTINE check_PFTimeVars
#endif

END MODULE MOD_PFTimeVars

#endif
! ---------- EOP ------------
