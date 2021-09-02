
#include <define.h>

SUBROUTINE lpwrite(idate,deltim,lwrite,rwrite)

  use timemanager

  implicit none
  integer, intent(in) :: idate(3)
  real,    intent(in) :: deltim
  logical, intent(inout) :: lwrite
  logical, intent(inout) :: rwrite
         
#if(defined WO_HOURLY)
  lwrite = .true.
#elif(defined WO_DAILY)
  lwrite = isendofday(idate, deltim)
#elif(defined WO_MONTHLY)
  lwrite = isendofmonth(idate, deltim)       
#elif(defined WO_YEARLY)
  lwrite = isendofyear(idate, deltim)
#endif

#if(defined WR_HOURLY)
  rwrite = .true.
#elif(defined WR_DAILY)
  rwrite = isendofday(idate, deltim)
#elif(defined WR_MONTHLY)
  rwrite = isendofmonth(idate, deltim)       
#elif(defined WR_YEARLY)
  rwrite = isendofyear(idate, deltim)
#endif

END SUBROUTINE lpwrite
! ------------------------------------------------------------------------
! EOP
