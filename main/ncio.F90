#include <define.h>

MODULE ncio
!-------------------------------------------------------------------------------
! Purpose:
!       tool functions for netcdf file format
!-------------------------------------------------------------------------------
   USE precision
   USE netcdf
   
   IMPLICIT NONE
   SAVE
   
   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: nccheck

CONTAINS 

   SUBROUTINE nccheck(status)
      INTEGER, INTENT(IN) :: status

      IF (status /= nf90_noerr) THEN
         print *, trim(nf90_strerror(status))
         stop 2
      END IF 
   END SUBROUTINE nccheck

END MODULE ncio
