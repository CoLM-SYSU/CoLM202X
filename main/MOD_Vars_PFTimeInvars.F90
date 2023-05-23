#include <define.h>

#ifdef PFT_CLASSIFICATION

MODULE MOD_Vars_PFTimeInvars
! -----------------------------------------------------------------
! !DESCRIPTION:
! Define PFT time invariables
!
! Created by Hua Yuan, 08/2019
! -----------------------------------------------------------------

  USE precision
  USE MOD_Vars_Global
  IMPLICIT NONE
  SAVE

  ! for PFT_CLASSIFICATION
  INTEGER , allocatable :: pftclass    (:)    !PFT type
! INTEGER , allocatable :: patch_pft_s (:)    !start PFT index of a patch
! INTEGER , allocatable :: patch_pft_e (:)    !end PFT index of a patch
! INTEGER , allocatable :: pft2patch   (:)    !patch index of a PFT
  REAL(r8), allocatable :: pftfrac     (:)    !PFT fractional cover
  REAL(r8), allocatable :: htop_p      (:)    !canopy top height [m]
  REAL(r8), allocatable :: hbot_p      (:)    !canopy bottom height [m]

! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: allocate_PFTimeInvars
  PUBLIC :: READ_PFTimeInvars
  PUBLIC :: WRITE_PFTimeInvars
  PUBLIC :: deallocate_PFTimeInvars
#ifdef CoLMDEBUG
  PUBLIC :: check_PFTimeInvars
#endif

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_PFTimeInvars
  ! --------------------------------------------------------------------
  ! Allocates memory for CoLM PFT 1d [numpft] variables
  ! --------------------------------------------------------------------

     USE spmd_task
     USE mod_landpft,   only : numpft
     USE precision
     IMPLICIT NONE

     IF (p_is_worker) THEN
        IF (numpft > 0) THEN
           allocate (pftclass      (numpft))
           allocate (pftfrac       (numpft))
           allocate (htop_p        (numpft))
           allocate (hbot_p        (numpft))
        ENDIF
     ENDIF

  END SUBROUTINE allocate_PFTimeInvars

  SUBROUTINE READ_PFTimeInvars (file_restart)

     use ncio_vector
     USE mod_landpft
     IMPLICIT NONE

     character(LEN=*), intent(in) :: file_restart

     call ncio_read_vector (file_restart, 'pftclass', landpft, pftclass) !
     call ncio_read_vector (file_restart, 'pftfrac ', landpft, pftfrac ) !
     call ncio_read_vector (file_restart, 'htop_p  ', landpft, htop_p  ) !
     call ncio_read_vector (file_restart, 'hbot_p  ', landpft, hbot_p  ) !

  end subroutine READ_PFTimeInvars

  SUBROUTINE WRITE_PFTimeInvars (file_restart)

     use ncio_vector
     use mod_landpft
     USE mod_namelist
     USE MOD_Vars_Global
     IMPLICIT NONE

     ! Local variables
     character(len=*), intent(in) :: file_restart
     integer :: compress

     compress = DEF_REST_COMPRESS_LEVEL

     call ncio_create_file_vector (file_restart, landpft)
     CALL ncio_define_dimension_vector (file_restart, landpft, 'pft')

     call ncio_write_vector (file_restart, 'pftclass', 'pft', landpft, pftclass, compress) !
     call ncio_write_vector (file_restart, 'pftfrac ', 'pft', landpft, pftfrac , compress) !
     call ncio_write_vector (file_restart, 'htop_p  ', 'pft', landpft, htop_p  , compress) !
     call ncio_write_vector (file_restart, 'hbot_p  ', 'pft', landpft, hbot_p  , compress) !

  end subroutine WRITE_PFTimeInvars

  SUBROUTINE deallocate_PFTimeInvars
! --------------------------------------------------
! Deallocates memory for CoLM PFT 1d [numpft] variables
! --------------------------------------------------
     USE spmd_task
     USE mod_landpft

     IF (p_is_worker) THEN
        IF (numpft > 0) THEN
           deallocate (pftclass)
           deallocate (pftfrac )
           deallocate (htop_p  )
           deallocate (hbot_p  )
        ENDIF
     ENDIF

  END SUBROUTINE deallocate_PFTimeInvars

#ifdef CoLMDEBUG
  SUBROUTINE check_PFTimeInvars ()

     use mod_colm_debug
     IMPLICIT NONE

     call check_vector_data ('pftfrac', pftfrac) !
     call check_vector_data ('htop_p ', htop_p ) !
     call check_vector_data ('hbot_p ', hbot_p ) !

  end subroutine check_PFTimeInvars
#endif

END MODULE MOD_Vars_PFTimeInvars

#endif
! ---------- EOP ------------
