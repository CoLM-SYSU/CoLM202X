MODULE mod_utils

   ! ---- PUBLIC subroutines ----

   PUBLIC :: normalize_longitude

   interface expand_list
      MODULE procedure expand_list_int32
      MODULE procedure expand_list_real8
   END interface expand_list

   PUBLIC :: append_to_list

   PUBLIC :: insert_into_sorted_list1
   PUBLIC :: insert_into_sorted_list2 

   PUBLIC :: find_in_sorted_list1 
   PUBLIC :: find_in_sorted_list2 

   PUBLIC :: find_nearest_south
   PUBLIC :: find_nearest_north
   PUBLIC :: find_nearest_west
   PUBLIC :: find_nearest_east

   PUBLIC :: lon_between_floor
   PUBLIC :: lon_between_ceil

   interface quicksort
      MODULE procedure quicksort_int32
      MODULE procedure quicksort_int64
   END interface quicksort

   PUBLIC :: quickselect
   PUBLIC :: median

   PUBLIC :: areaquad

   interface unpack_inplace
      MODULE procedure unpack_inplace_int32
      MODULE procedure unpack_inplace_real8
      MODULE procedure unpack_inplace_lastdim_real8
   END interface unpack_inplace

CONTAINS

   !---------------------------------
   SUBROUTINE normalize_longitude (lon)

      USE precision
      IMPLICIT NONE

      REAL(r8), intent(inout) :: lon
      
      DO while (lon >= 180.0)
         lon = lon - 360.0
      ENDDO

      DO while (lon < -180.0)
         lon = lon + 360.0
      ENDDO

   END SUBROUTINE normalize_longitude

   !--------------------------------------------------
   SUBROUTINE expand_list_int32 (list, percent)

      USE precision
      IMPLICIT NONE

      INTEGER, allocatable, intent(inout) :: list (:)
      REAL(r8), intent(in) :: percent

      ! Local variables
      INTEGER :: n0, n1
      INTEGER, allocatable :: temp (:)

      n0 = size(list)

      allocate (temp(n0))
      temp = list
      
      n1 = ceiling(n0 * (1+percent))

      deallocate(list)
      allocate (list(n1))
      list(1:n0) = temp      

      deallocate (temp)

   END SUBROUTINE expand_list_int32

   !--------------------------------------------------
   SUBROUTINE expand_list_real8 (list, percent)

      USE precision
      IMPLICIT NONE

      REAL(r8), allocatable, intent(inout) :: list (:)
      REAL(r8), intent(in) :: percent

      ! Local variables
      INTEGER :: n0, n1
      REAL(r8), allocatable :: temp (:)

      n0 = size(list)

      allocate (temp(n0))
      temp = list
      
      n1 = ceiling(n0 * (1+percent))
      
      deallocate(list)
      allocate (list(n1))
      list(1:n0) = temp      

      deallocate (temp)

   END SUBROUTINE expand_list_real8

   !--------------------------------------------------
   SUBROUTINE append_to_list (list1, list2)

      IMPLICIT NONE

      INTEGER, allocatable, intent(inout) :: list1 (:)
      INTEGER, intent(in) :: list2 (:)

      ! Local variables
      INTEGER :: n1, n2
      INTEGER, allocatable :: temp (:)

      IF (.not. allocated(list1)) THEN
         n1 = 0
      ELSE
         n1 = size(list1)
      ENDIF

      n2 = size(list2)

      IF (n1 > 0) THEN
         allocate (temp(n1))
         temp = list1

         deallocate(list1)
         allocate (list1(n1+n2))
         list1(1:n1) = temp      

         deallocate (temp)
      ELSE
         IF (n2 > 0) allocate (list1(n2))
      ENDIF

      IF (n1 + n2 > 0) THEN
         list1(n1+1:n1+n2) = list2
      ENDIF

   END SUBROUTINE append_to_list

   !--------------------------------------------------
   SUBROUTINE insert_into_sorted_list1 (x, n, list, iloc, is_new_out)

      IMPLICIT NONE

      INTEGER, intent(in) :: x
      INTEGER, intent(inout) :: n
      INTEGER, intent(inout) :: list(:)
      INTEGER, intent(out)   :: iloc
      LOGICAL, intent(out), optional :: is_new_out

      ! Local variables
      LOGICAL :: is_new
      INTEGER :: ileft, iright

      IF (n == 0) THEN
         iloc = 1
         is_new = .true.
      ELSEIF (x <= list(1)) THEN
         iloc = 1
         is_new = (x /= list(1))
      ELSEIF (x > list(n)) THEN
         iloc = n + 1
         is_new = .true. 
      ELSEIF (x == list(n)) THEN
         iloc = n
         is_new = .false. 
      ELSE
         ileft  = 1
         iright = n

         DO while (.true.)
            IF (iright - ileft > 1) THEN
               iloc = (ileft + iright) / 2
               IF (x > list(iloc)) THEN
                  ileft = iloc
               ELSEIF (x < list(iloc)) THEN
                  iright = iloc
               ELSE
                  is_new = .false.
                  exit 
               ENDIF
            ELSE
               iloc = iright
               is_new = .true.
               exit
            ENDIF
         ENDDO
      ENDIF

      IF (is_new) THEN
         IF (iloc <= n) THEN
            list(iloc+1:n+1) = list(iloc:n)
         ENDIF

         list(iloc) = x
         n = n + 1
      ENDIF

      IF (present(is_new_out)) THEN
         is_new_out = is_new
      ENDIF

   END SUBROUTINE insert_into_sorted_list1

   !--------------------------------------------------
   SUBROUTINE insert_into_sorted_list2 (x, y, n, xlist, ylist, iloc, is_new_out)

      IMPLICIT NONE

      INTEGER, intent(in) :: x, y
      INTEGER, intent(inout) :: n
      INTEGER, intent(inout) :: xlist(:), ylist(:)
      INTEGER, intent(out)   :: iloc
      LOGICAL, intent(out), optional :: is_new_out

      ! Local variables
      LOGICAL :: is_new
      INTEGER :: ileft, iright

      IF (n == 0) THEN
         iloc = 1
         is_new = .true.
      ELSEIF ((y < ylist(1)) .or. ((y == ylist(1)) .and. (x <= xlist(1)))) THEN
         iloc = 1
         is_new = (x /= xlist(1)) .or. (y /= ylist(1))
      ELSEIF ((y > ylist(n)) .or. ((y == ylist(n)) .and. (x > xlist(n)))) THEN
         iloc = n + 1
         is_new = .true. 
      ELSEIF ((x == xlist(n)) .and. (y == ylist(n))) THEN
         iloc = n
         is_new = .false. 
      ELSE
         ileft  = 1
         iright = n

         DO while (.true.)
            IF (iright - ileft > 1) THEN
               iloc = (ileft + iright) / 2
               IF ((y > ylist(iloc)) .or. ((y == ylist(iloc)) .and. (x > xlist(iloc)))) THEN
                  ileft = iloc
               ELSEIF ((y < ylist(iloc)) .or. ((y == ylist(iloc)) .and. (x < xlist(iloc)))) THEN
                  iright = iloc
               ELSE
                  is_new = .false.
                  exit 
               ENDIF
            ELSE
               iloc = iright
               is_new = .true.
               exit
            ENDIF
         ENDDO
      ENDIF

      IF (is_new) THEN
         IF (iloc <= n) THEN
            xlist(iloc+1:n+1) = xlist(iloc:n)
            ylist(iloc+1:n+1) = ylist(iloc:n)
         ENDIF

         xlist(iloc) = x
         ylist(iloc) = y
         n = n + 1
      ENDIF

      IF (present(is_new_out)) THEN
         is_new_out = is_new
      ENDIF

   END SUBROUTINE insert_into_sorted_list2

   !--------------------------------------------------
   FUNCTION find_in_sorted_list1 (x, n, list) result(iloc)

      IMPLICIT NONE

      INTEGER :: iloc

      INTEGER, intent(in) :: x
      INTEGER, intent(in) :: n
      INTEGER, intent(in) :: list (n)

      ! Local variables
      INTEGER :: i, ileft, iright

      iloc = 0
      IF (n > 0) THEN
         IF ((x >= list(1)) .and. (x <= list(n))) THEN
            IF (x == list(1)) THEN
               iloc = 1
            ELSEIF (x == list(n)) THEN
               iloc = n
            ELSE
               ileft  = 1
               iright = n

               DO while (iright - ileft > 1)
                  i = (ileft + iright) / 2
                  IF (x == list(i)) THEN
                     iloc = i
                     exit
                  ELSEIF (x > list(i)) THEN
                     ileft = i
                  ELSEIF (x < list(i)) THEN
                     iright = i
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      ENDIF

   END FUNCTION find_in_sorted_list1

   !--------------------------------------------------
   FUNCTION find_in_sorted_list2 (x, y, n, xlist, ylist) result(iloc)

      IMPLICIT NONE

      INTEGER :: iloc

      INTEGER, intent(in) :: x, y
      INTEGER, intent(in) :: n
      INTEGER, intent(in) :: xlist(:), ylist(:)

      ! Local variables
      INTEGER :: i, ileft, iright

      iloc = 0
      IF (n < 1) RETURN

      IF ((y < ylist(1)) .or. ((y == ylist(1)) .and. (x < xlist(1)))) THEN
         iloc = 0
      ELSEIF ((y > ylist(n)) .or. ((y == ylist(n)) .and. (x > xlist(n)))) THEN
         iloc = 0
      ELSEIF ((x == xlist(1)) .and. (y == ylist(1))) THEN
         iloc = 1
      ELSEIF ((x == xlist(n)) .and. (y == ylist(n))) THEN
         iloc = n
      ELSE
         ileft  = 1
         iright = n

         DO while (.true.)
            IF (iright - ileft > 1) THEN
               i = (ileft + iright) / 2
               IF ((y == ylist(i)) .and. (x == xlist(i))) THEN
                  iloc = i
                  exit
               ELSEIF ((y > ylist(i)) .or. ((y == ylist(i)) .and. (x > xlist(i)))) THEN
                  ileft = i
               ELSEIF ((y < ylist(i)) .or. ((y == ylist(i)) .and. (x < xlist(i)))) THEN
                  iright = i
               ENDIF
            ELSE
               iloc = 0
               exit
            ENDIF
         ENDDO
      ENDIF

   END FUNCTION find_in_sorted_list2
   
   !-----------------------------------------------------
   FUNCTION find_nearest_south (y, n, lat) result(iloc)

      USE precision
      IMPLICIT NONE
      
      INTEGER :: iloc

      REAL(r8), intent(in) :: y
      INTEGER,  intent(in) :: n
      REAL(r8), intent(in) :: lat (n)

      ! Local variables
      INTEGER :: i, iright, ileft

      IF (lat(1) < lat(n))  THEN
         IF (y <= lat(1)) THEN
            iloc = 1
         ELSEIF (y >= lat(n)) THEN
            iloc = n
         ELSE
            ileft = 1;  iright = n

            DO while (iright - ileft > 1)
               i = (iright + ileft) / 2
               IF (y >= lat(i)) THEN
                  ileft = i
               ELSE
                  iright = i
               ENDIF
            ENDDO

            iloc = ileft
         ENDIF
      ELSE
         IF (y >= lat(1)) THEN
            iloc = 1
         ELSEIF (y <= lat(n)) THEN
            iloc = n
         ELSE
            ileft = 1;  iright = n

            DO while (iright - ileft > 1)
               i = (iright + ileft) / 2
               IF (y >= lat(i)) THEN
                  iright = i
               ELSE
                  ileft = i
               ENDIF
            ENDDO

            iloc = iright
         ENDIF
      ENDIF

   END FUNCTION find_nearest_south

   !-----------------------------------------------------
   FUNCTION find_nearest_north (y, n, lat)  result(iloc)

      USE precision
      IMPLICIT NONE

      REAL(r8), intent(in) :: y
      INTEGER,  intent(in) :: n
      REAL(r8), intent(in) :: lat (n)
      
      INTEGER :: iloc

      ! Local variables
      INTEGER :: i, iright, ileft

      IF (lat(1) < lat(n))  THEN
         IF (y <= lat(1)) THEN
            iloc = 1
         ELSEIF (y >= lat(n)) THEN
            iloc = n
         ELSE
            ileft = 1;  iright = n

            DO while (iright - ileft > 1)
               i = (iright + ileft) / 2
               IF (y > lat(i)) THEN
                  ileft = i
               ELSE
                  iright = i
               ENDIF
            ENDDO

            iloc = iright
         ENDIF
      ELSE
         IF (y >= lat(1)) THEN
            iloc = 1
         ELSEIF (y <= lat(n)) THEN
            iloc = n
         ELSE
            ileft = 1;  iright = n

            DO while (iright - ileft > 1)
               i = (iright + ileft) / 2
               IF (y > lat(i)) THEN
                  iright = i
               ELSE
                  ileft = i
               ENDIF
            ENDDO

            iloc = ileft
         ENDIF
      ENDIF

   END FUNCTION find_nearest_north

   !-----------------------------------------
   LOGICAL FUNCTION lon_between_floor (lon, west, east)

      USE precision
      IMPLICIT NONE

      REAL(r8), intent(in) :: lon, west, east ! [-180, 180)

      IF (west >= east) THEN
         lon_between_floor = (lon >= west) .or. (lon < east) 
      ELSE
         lon_between_floor = (lon >= west) .and. (lon < east)
      ENDIF

   END FUNCTION lon_between_floor 

   !-----------------------------------------
   LOGICAL FUNCTION lon_between_ceil (lon, west, east)

      USE precision
      IMPLICIT NONE

      REAL(r8), intent(in) :: lon, west, east ! [-180, 180)

      IF (west >= east) THEN
         lon_between_ceil = (lon > west) .or. (lon <= east) 
      ELSE
         lon_between_ceil = (lon > west) .and. (lon <= east)
      ENDIF

   END FUNCTION lon_between_ceil

   !-----------------------------------------------------
   FUNCTION find_nearest_west (x, n, lon) result(iloc)

      USE precision
      IMPLICIT NONE

      REAL(r8), intent(in) :: x
      INTEGER,  intent(in) :: n
      REAL(r8), intent(in) :: lon (n)
      
      INTEGER :: iloc

      ! Local variables
      INTEGER :: i, iright, ileft

      IF (n == 1) THEN
         iloc = 1
         RETURN 
      ENDIF

      IF (lon_between_floor (x, lon(n), lon(1))) THEN
         iloc = n
         RETURN
      ENDIF

      ileft = 1; iright = n
      DO while (iright - ileft > 1)
         i = (iright + ileft)/2
         IF (lon_between_floor(x,lon(i),lon(iright))) THEN
            ileft = i
         ELSE
            iright = i
         ENDIF
      ENDDO
      
      iloc = ileft

   END FUNCTION find_nearest_west

   !-----------------------------------------------------
   FUNCTION find_nearest_east (x, n, lon) result(iloc)

      USE precision
      IMPLICIT NONE

      REAL(r8), intent(in) :: x
      INTEGER,  intent(in) :: n
      REAL(r8), intent(in) :: lon (n)
      
      INTEGER :: iloc

      ! Local variables
      INTEGER :: i, iright, ileft

      IF (n == 1) THEN
         iloc = 1
         RETURN 
      ENDIF

      IF (lon_between_ceil (x, lon(n), lon(1))) THEN
         iloc = 1
         RETURN
      ENDIF

      ileft = 1; iright = n
      DO while (iright - ileft > 1)
         i = (iright + ileft)/2
         IF (lon_between_ceil(x,lon(i),lon(iright))) THEN
            ileft = i
         ELSE
            iright = i
         ENDIF
      ENDDO
      
      iloc = iright

   END FUNCTION find_nearest_east


   !-----------------------------------------------------
   recursive SUBROUTINE quicksort_int32 (nA, A, order)

      USE precision
      IMPLICIT NONE

      INTEGER, intent(in) :: nA
      INTEGER, intent(inout) :: A     (nA)
      INTEGER, intent(inout) :: order (nA)

      ! Local variables
      INTEGER :: left, right
      INTEGER :: pivot
      INTEGER :: marker
      INTEGER :: itemp

      IF (nA > 1) THEN

         pivot = A (nA/2)
         left  = 0
         right = nA + 1

         DO while (left < right)
            right = right - 1
            DO while (A(right) > pivot)
               right = right - 1
            ENDDO

            left = left + 1
            DO while (A(left) < pivot)
               left = left + 1
            ENDDO

            IF (left < right) THEN
               itemp    = A(left)
               A(left)  = A(right)
               A(right) = itemp

               itemp        = order(left)
               order(left)  = order(right)
               order(right) = itemp
            ENDIF
         ENDDO

         marker = right

         CALL quicksort_int32 (marker,    A(1:marker),    order(1:marker))
         CALL quicksort_int32 (nA-marker, A(marker+1:nA), order(marker+1:nA))

      ENDIF

   END SUBROUTINE quicksort_int32

   !-----------------------------------------------------
   recursive SUBROUTINE quicksort_int64 (nA, A, order)

      USE precision
      IMPLICIT NONE

      INTEGER*8, intent(in) :: nA
      INTEGER*8, intent(inout) :: A     (nA)
      INTEGER*8, intent(inout) :: order (nA)

      ! Local variables
      INTEGER*8 :: left, right
      INTEGER*8 :: pivot
      INTEGER*8 :: marker
      INTEGER*8 :: itemp

      IF (nA > 1) THEN

         pivot = A (nA/2)
         left  = 0
         right = nA + 1

         DO while (left < right)
            right = right - 1
            DO while (A(right) > pivot)
               right = right - 1
            ENDDO

            left = left + 1
            DO while (A(left) < pivot)
               left = left + 1
            ENDDO

            IF (left < right) THEN
               itemp    = A(left)
               A(left)  = A(right)
               A(right) = itemp

               itemp        = order(left)
               order(left)  = order(right)
               order(right) = itemp
            ENDIF
         ENDDO

         marker = right

         CALL quicksort_int64 (marker,    A(1:marker),    order(1:marker))
         CALL quicksort_int64 (nA-marker, A(marker+1:nA), order(marker+1:nA))

      ENDIF

   END SUBROUTINE quicksort_int64

   !-----------------------------------------------------
   recursive FUNCTION quickselect (nA, A, k) result(selected)

      USE precision
      IMPLICIT NONE

      REAL(r8) :: selected

      INTEGER , intent(in)    :: nA
      REAL(r8), intent(inout) :: A (nA)
      INTEGER,  intent(in)    :: k

      ! Local variables
      INTEGER  :: left, right
      REAL(r8) :: pivot
      INTEGER  :: marker
      REAL(r8) :: rtemp 

      IF (nA > 1) THEN

         pivot = A (nA/2)
         left  = 0
         right = nA + 1

         DO while (left < right)
            right = right - 1
            DO while (A(right) > pivot)
               right = right - 1
            ENDDO

            left = left + 1
            DO while (A(left) < pivot)
               left = left + 1
            ENDDO

            IF (left < right) THEN
               rtemp    = A(left)
               A(left)  = A(right)
               A(right) = rtemp
            ENDIF
         ENDDO

         marker = right

         IF (k <= marker) THEN
            selected = quickselect (marker, A(1:marker), k)
         ELSE
            selected = quickselect (nA-marker, A(marker+1:nA), k-marker)
         ENDIF

      ELSE
         selected = A(1) 
      ENDIF

   END FUNCTION quickselect
   
   
   ! ------------------------
   FUNCTION median(x, n, spval) result(mval)

      USE precision
      IMPLICIT NONE

      REAL(r8) :: mval

      INTEGER,  intent(in) :: n
      REAL(r8), intent(in) :: x(n)
      REAL(r8), intent(in), optional :: spval

      ! Local variables
      INTEGER  :: nc
      REAL(r8), allocatable :: xtemp(:)
      LOGICAL,  allocatable :: msk  (:)
      REAL(r8) :: right, left

      IF (present(spval)) THEN
         allocate (msk (n))
         msk = (x /= spval)
         nc  = count(msk)
         IF (nc /= 0) THEN
            
            allocate (xtemp(nc))
            xtemp = pack(x, msk)

            deallocate (msk)
         ELSE
           
            mval = spval

            deallocate(msk)
            RETURN
         ENDIF
      ELSE
         nc = n
         allocate (xtemp(nc))
         xtemp = x
      ENDIF
      
      IF (mod(nc,2) == 0) THEN
         left  = quickselect(nc,xtemp,nc/2)
         right = quickselect(nc,xtemp,nc/2+1)
         mval = (left + right) / 2.0_r8
      ELSE
         mval = quickselect(nc,xtemp,nc/2+1)
      ENDIF

      deallocate (xtemp)

   END FUNCTION median

   
   !-----------------------------------------------------
   FUNCTION areaquad (lats, latn, lonw, lone) result(area)

      USE precision
      USE MathConstants, only : deg2rad
      IMPLICIT NONE

      REAL(r8) :: area
      REAL(r8), parameter :: re = 6.37122e3 ! kilometer 
      REAL(r8), intent(in) :: lats, latn, lonw, lone

      ! Local variables
      REAL(r8) :: dx, dy

      IF (lone < lonw) THEN
         dx = (lone + 360 - lonw) * deg2rad
      ELSE
         dx = (lone - lonw) * deg2rad
      ENDIF

      dy = sin(latn * deg2rad) - sin(lats * deg2rad)

      area = dx * dy * re * re

   END FUNCTION areaquad

   ! --- spherical distance  ---
   FUNCTION arclen (lat1, lon1, lat2, lon2)

      USE precision
      IMPLICIT NONE

      REAL(r8) :: arclen
      REAL(r8), intent(in) :: lat1, lon1, lat2, lon2
      
      REAL(r8), parameter :: re = 6.37122e3 ! kilometer 
      
      arclen = re * acos (sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2) * cos(lon1-lon2))

   END FUNCTION arclen

   !-----------------------------------------------------
   SUBROUTINE unpack_inplace_int32 (din, msk, dout)

      IMPLICIT NONE

      INTEGER, intent(in) :: din (:)
      LOGICAL, intent(in) :: msk (:)
      INTEGER, intent(inout) :: dout (:)

      ! Local variables
      INTEGER :: n, i

      n = 0
      DO i = 1, size(msk)
         IF (msk(i)) THEN
            n = n + 1
            dout(i) = din(n)
         ENDIF
      ENDDO

   END SUBROUTINE unpack_inplace_int32

   !-----------------------------------------------------
   SUBROUTINE unpack_inplace_real8 (din, msk, dout)

      USE precision
      IMPLICIT NONE

      REAL(r8), intent(in) :: din (:)
      LOGICAL,  intent(in) :: msk (:)
      REAL(r8), intent(inout) :: dout (:)

      ! Local variables
      INTEGER :: n, i

      n = 0
      DO i = 1, size(msk)
         IF (msk(i)) THEN
            n = n + 1
            dout(i) = din(n)
         ENDIF
      ENDDO

   END SUBROUTINE unpack_inplace_real8
   
   !-----------------------------------------------------
   SUBROUTINE unpack_inplace_lastdim_real8 (din, msk, dout)

      USE precision
      IMPLICIT NONE

      REAL(r8), intent(in) :: din (:,:)
      LOGICAL,  intent(in) :: msk (:)
      REAL(r8), intent(inout) :: dout (:,:)

      ! Local variables
      INTEGER :: n, i

      n = 0
      DO i = 1, size(msk)
         IF (msk(i)) THEN
            n = n + 1
            dout(:,i) = din(:,n)
         ENDIF
      ENDDO

   END SUBROUTINE unpack_inplace_lastdim_real8

   !---------------------------------------------------
   INTEGER FUNCTION num_max_frequency (data_in)
      
      IMPLICIT NONE

      INTEGER, intent(in) :: data_in(:)

      ! Local Variables
      INTEGER, allocatable :: data_(:), cnts(:)
      INTEGER :: ndata, i, n, iloc
      LOGICAL :: is_new

      ndata = size(data_in)
      allocate (data_(ndata))
      allocate (cnts (ndata))

      n = 0
      cnts(:) = 0
      DO i = 1, ndata
         CALL insert_into_sorted_list1 (data_in(i), n, data_, iloc, is_new) 
         IF (is_new) THEN
            IF (iloc < n) cnts(iloc+1:ndata) = cnts(iloc:ndata-1)
            cnts(iloc) = 1
         ELSE
            cnts(iloc) = cnts(iloc) + 1
         ENDIF
      ENDDO

      num_max_frequency = data_(maxloc(cnts,dim=1))

      deallocate(data_)
      deallocate(cnts )

   END FUNCTION num_max_frequency 

END MODULE mod_utils
