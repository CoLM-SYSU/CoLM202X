#include <define.h>

MODULE MOD_DA_Main
!-----------------------------------------------------------------------------
! DESCRIPTION:
!     Main procedures for data assimilation
!
! AUTHOR:
!     Lu Li, 12/2024
!     Zhilong Fan, Lu Li, 03/2024
!-----------------------------------------------------------------------------
   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_DA_GRACE
   USE MOD_DA_SMAP
   USE MOD_DA_FY3D
   IMPLICIT NONE
   SAVE


!-----------------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------------

   SUBROUTINE init_DA ()

!-----------------------------------------------------------------------------
   IMPLICIT NONE

      IF (DEF_DA_GRACE) THEN
         CALL init_DA_GRACE ()
      ENDIF

      IF (DEF_DA_SMAP) THEN
         CALL init_DA_SMAP  ()
      ENDIF

      IF (DEF_DA_FY3D) THEN
         CALL init_DA_FY3D  ()
      ENDIF

   END SUBROUTINE init_DA



!-----------------------------------------------------------------------------

   SUBROUTINE run_DA (idate, deltim)

!-----------------------------------------------------------------------------
   IMPLICIT NONE

!---------------------Dummy arguments-----------------------------------------
   integer,  intent(in) :: idate(3)
   real(r8), intent(in) :: deltim

!-----------------------------------------------------------------------------

      IF (DEF_DA_GRACE) THEN
         CALL run_DA_GRACE (idate, deltim)
      ENDIF

      IF (DEF_DA_SMAP) THEN
         CALL run_DA_SMAP  (idate, deltim)
      ENDIF

      IF (DEF_DA_FY3D) THEN
         CALL run_DA_SMAP  (idate, deltim)
      ENDIF

   END SUBROUTINE run_DA

!-----------------------------------------------------------------------------

   SUBROUTINE end_DA ()

!-----------------------------------------------------------------------------
   IMPLICIT NONE

      IF (DEF_DA_GRACE) THEN
         CALL end_DA_GRACE ()
      ENDIF

      IF (DEF_DA_SMAP) THEN
         CALL end_DA_SMAP  ()
      ENDIF

      IF (DEF_DA_FY3D) THEN
         CALL end_DA_FY3D  ()
      ENDIF

   END SUBROUTINE end_DA



!-----------------------------------------------------------------------------
END MODULE MOD_DA_Main
