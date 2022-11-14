#include <define.h> 

#ifdef PC_CLASSIFICATION

MODULE MOD_PCTimeInvars 
! -------------------------------
! Created by Hua Yuan, 08/2019
! -------------------------------

  USE precision
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
#ifdef CLMDEBUG
  PUBLIC :: check_PCTimeInvars
#endif

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

  CONTAINS

!-----------------------------------------------------------------------

  SUBROUTINE allocate_PCTimeInvars
  ! --------------------------------------------------------------------
  ! Allocates memory for CLM Plant Community (PC) [numpc] variables
  ! --------------------------------------------------------------------

     USE precision
     USE spmd_task
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

     use ncio_vector
     USE mod_landpc
     IMPLICIT NONE

     character(LEN=*), intent(in) :: file_restart
  
     call ncio_read_vector (file_restart, 'pcfrac', N_PFT, landpc, pcfrac) !
     call ncio_read_vector (file_restart, 'htop_c', N_PFT, landpc, htop_c) !
     call ncio_read_vector (file_restart, 'hbot_c', N_PFT, landpc, hbot_c) !

  end subroutine READ_PCTimeInvars

  SUBROUTINE WRITE_PCTimeInvars (file_restart)

     use ncio_vector
     use mod_landpc
     USE mod_namelist
     USE GlobalVars
     IMPLICIT NONE

     ! Local variables
     character(len=*), intent(in) :: file_restart
     integer :: compress
     
     compress = DEF_REST_COMPRESS_LEVEL 

     call ncio_create_file_vector (file_restart, landpc)
     CALL ncio_define_pixelset_dimension (file_restart, landpc)
     CALL ncio_define_dimension_vector (file_restart, 'pft', N_PFT)
     
     call ncio_write_vector (file_restart, 'pcfrac', 'pft', N_PFT, 'vector', landpc, pcfrac, compress) !
     call ncio_write_vector (file_restart, 'htop_c', 'pft', N_PFT, 'vector', landpc, htop_c, compress) !
     call ncio_write_vector (file_restart, 'hbot_c', 'pft', N_PFT, 'vector', landpc, hbot_c, compress) !
    
  end subroutine WRITE_PCTimeInvars

  SUBROUTINE deallocate_PCTimeInvars
! --------------------------------------------------
! Deallocates memory for CLM Plant Community (PC) variables
! --------------------------------------------------

     USE spmd_task
     USE mod_landpc

     IF (p_is_worker) THEN
        IF (numpc > 0) THEN
           deallocate (pcfrac   )  
           deallocate (htop_c   ) 
           deallocate (hbot_c   ) 
        ENDIF
     ENDIF 

  END SUBROUTINE deallocate_PCTimeInvars

#ifdef CLMDEBUG
  SUBROUTINE check_PCTimeInvars ()

     use mod_colm_debug
     IMPLICIT NONE

     call check_vector_data ('pcfrc ', pcfrac) !
     call check_vector_data ('htop_c', htop_c) !
     call check_vector_data ('hbot_c', hbot_c) !

  end subroutine check_PCTimeInvars 
#endif

END MODULE MOD_PCTimeInvars

#endif
! ---------- EOP ------------
