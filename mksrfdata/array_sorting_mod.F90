MODULE array_sorting_mod

   use precision

   IMPLICIT none
!
!  PUBLIC: interfaces
!
!  Initial Author : Liu Li, 02/2014
!  --------------------------------

   public :: merge_sort
   public :: bubble_sort


contains

SUBROUTINE merge_sort(array_src, array_dst, array_size)
! ------------------------
   IMPLICIT NONE

   integer,  intent(in)    :: array_size
   real(r8), intent(in)    :: array_src(array_size)
   real(r8), intent(out)   :: array_dst(array_size)
   real(r8)                :: array_buf(array_size)
   integer                 :: i
! ------------------------

   array_dst(:) = array_src(:)
 
   if (array_size .le. 1) return

   call merge_sort_recursively(array_dst, array_buf, array_size, 1, array_size)

END SUBROUTINE merge_sort



RECURSIVE SUBROUTINE merge_sort_recursively(array_dst, array_buf, array_size, start_indx, end_indx)
! ------------------------
   IMPLICIT NONE
   integer,  intent(in)     :: array_size, start_indx, end_indx
   real(r8), intent(out )   :: array_buf(array_size)
   real(r8), intent(inout ) :: array_dst(array_size)
   integer                  :: mid_indx
   integer                  :: i, j, k, l

   mid_indx = (start_indx+end_indx) / 2

   if (mid_indx > start_indx) then
      call merge_sort_recursively(array_dst, array_buf, array_size, start_indx, mid_indx)
   endif
   if (mid_indx+1 < end_indx) then
      call merge_sort_recursively(array_dst, array_buf, array_size, mid_indx+1, end_indx)
   endif

   i = start_indx; j = mid_indx+1; k = start_indx
   do while (i .le. mid_indx .and. j .le. end_indx)
      if (array_dst(i) .le. array_dst(j)) then
         array_buf(k) = array_dst(i)
         i = i + 1
      else
         array_buf(k) = array_dst(j)
         j = j + 1
      endif
      k = k + 1
   enddo

   do l = i, mid_indx
      array_buf(k) = array_dst(l) 
      k = k + 1
   enddo
   do l = j, end_indx
      array_buf(k) = array_dst(l) 
      k = k + 1
   enddo

   array_dst(start_indx:end_indx) = array_buf(start_indx:end_indx) 

END SUBROUTINE merge_sort_recursively



SUBROUTINE bubble_sort(x, a, n ) 
! ------------------------
IMPLICIT NONE

   integer,  intent(in) :: n
   real(r8), intent(in) :: x(n)
   integer :: i
   integer :: Loca(1), Location
   real(r8), intent(out) :: a(n)
   real(r8) temp
! ------------------------

   a(:) = x(:)

   ! sort into ascending order
   do i = 1, n-1
      Loca = minloc(a(i:n))
      Location = Loca(1) + i - 1

      if(Location /= i)then
         temp = a(i)
         a(i) = a(Location)
         a(Location) = temp
      endif
   enddo
END SUBROUTINE bubble_sort 

END MODULE array_sorting_mod
