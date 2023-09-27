#include <define.h>

#ifdef DataAssimilation
MODULE MOD_DataAssimilation

   USE MOD_Precision
   USE MOD_DA_GRACE
   IMPLICIT NONE

CONTAINS

   ! ----------
   SUBROUTINE init_DataAssimilation ()
      
      IMPLICIT NONE
      
      CALL init_DA_GRACE ()

   END SUBROUTINE init_DataAssimilation

   ! ----------
   SUBROUTINE do_DataAssimilation (idate, deltim)
      
      IMPLICIT NONE
      
      INTEGER,  INTENT(in) :: idate(3)
      REAL(r8), INTENT(in) :: deltim
      
      CALL do_DA_GRACE (idate, deltim)

   END SUBROUTINE do_DataAssimilation

   ! ---------
   SUBROUTINE final_DataAssimilation ()

      IMPLICIT NONE

      CALL final_DA_GRACE ()

   END SUBROUTINE final_DataAssimilation

END MODULE MOD_DataAssimilation
#endif
