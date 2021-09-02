FUNCTION median(x, n)
! ------------------------
use precision
use array_sorting_mod
IMPLICIT NONE

      integer, intent(in) :: n
      real(r8), intent(in) :: x(n)
      integer :: i
      real(r8) a(n)
      real(r8) median
! ------------------------
      ! sort into ascending order
      call merge_sort(x, a, n)

      ! get median value
      if(mod(n,2) == 0) then
         median = (a(n/2) + a(n/2+1)) / 2.0
      else
         median = a(n/2+1)
      end if

END FUNCTION median
