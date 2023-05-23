#include <define.h>

#ifdef PC_CLASSIFICATION

MODULE MOD_Vars_PCTimeInvars
! -----------------------------------------------------------------
! !DESCRIPTION:
! Define Plant Community time invariables
!
! Created by Hua Yuan, 08/2019
! -----------------------------------------------------------------

  USE MOD_Precision
  USE GlobalVars
  IMPLICIT NONE
  SAVE

  ! for PC_CLASSIFICATION
  REAL(r8), allocatable :: pcfrac(:,:)    !PC fractional cover
  REAL(r8), allocatable :: htop_c(:,:)    !canopy top height [m]
  REAL(r8), allocatable :: hbot_c(:,:)    !canopy bottom height [m]

! PUBLIC MEMBER FUNCTIONS:
  PUBLIC :: allocate_PCTimeInvars
  PUBLIC :: READ_PCTimeInvars
  PUBLIC :: WRITE_PCTimeInvars
  PUBLIC :: deallocate_PCTimeInvars
#ifdef CoLMDEBUG
  PUBLIC :: check_PCTimeInvars
#endif

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_PCTimeInvars
  ! --------------------------------------------------------------------
  ! Allocates memory for CoLM Plant Community (PC) [numpc] variables
  ! --------------------------------------------------------------------

     USE MOD_Precision
     USE MOD_SPMD_Task
     USE mod_landpc
     IMPLICIT NONE

     IF (p_is_worker) THEN

        IF (numpc > 0) THEN
           allocate (pcfrac   (0:N_PFT-1,numpc))
           allocate (htop_c   (0:N_PFT-1,numpc))
           allocate (hbot_c   (0:N_PFT-1,numpc))
        ENDIF
     ENDIF

  END SUBROUTINE allocate_PCTimeInvars

  SUBROUTINE READ_PCTimeInvars (file_restart)

     use MOD_NetCDFVector
     USE mod_landpc
     IMPLICIT NONE

     character(LEN=*), intent(in) :: file_restart

     call ncio_read_vector (file_restart, 'pcfrac', N_PFT, landpc, pcfrac) !
     call ncio_read_vector (file_restart, 'htop_c', N_PFT, landpc, htop_c) !
     call ncio_read_vector (file_restart, 'hbot_c', N_PFT, landpc, hbot_c) !

  end subroutine READ_PCTimeInvars

  SUBROUTINE WRITE_PCTimeInvars (file_restart)

     use MOD_NetCDFVector
     use mod_landpc
     USE MOD_Namelist
     USE GlobalVars
     IMPLICIT NONE

     ! Local variables
     character(len=*), intent(in) :: file_restart
     integer :: compress

     compress = DEF_REST_COMPRESS_LEVEL

     call ncio_create_file_vector (file_restart, landpc)
     CALL ncio_define_dimension_vector (file_restart, landpc, 'pc')
     CALL ncio_define_dimension_vector (file_restart, landpc, 'pft', N_PFT)

     call ncio_write_vector (file_restart, 'pcfrac', 'pft', N_PFT, 'pc', landpc, pcfrac, compress) !
     call ncio_write_vector (file_restart, 'htop_c', 'pft', N_PFT, 'pc', landpc, htop_c, compress) !
     call ncio_write_vector (file_restart, 'hbot_c', 'pft', N_PFT, 'pc', landpc, hbot_c, compress) !

  end subroutine WRITE_PCTimeInvars

  SUBROUTINE deallocate_PCTimeInvars
! --------------------------------------------------
! Deallocates memory for CoLM Plant Community (PC) variables
! --------------------------------------------------

     USE MOD_SPMD_Task
     USE mod_landpc

     IF (p_is_worker) THEN
        IF (numpc > 0) THEN
           deallocate (pcfrac   )
           deallocate (htop_c   )
           deallocate (hbot_c   )
        ENDIF
     ENDIF

  END SUBROUTINE deallocate_PCTimeInvars

#ifdef CoLMDEBUG
  SUBROUTINE check_PCTimeInvars ()

     use MOD_CoLMDebug
     IMPLICIT NONE

     call check_vector_data ('pcfrc ', pcfrac) !
     call check_vector_data ('htop_c', htop_c) !
     call check_vector_data ('hbot_c', hbot_c) !

  end subroutine check_PCTimeInvars
#endif

END MODULE MOD_Vars_PCTimeInvars

#endif
! ---------- EOP ------------
