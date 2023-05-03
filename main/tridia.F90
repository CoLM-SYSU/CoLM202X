      SUBROUTINE tridia (n, a, b, c, r, u)
      
      USE precision
      IMPLICIT NONE
      INTEGER,  intent(in) :: n       !length of diagonal element vector
      REAL(r8), intent(in) :: a(1:n)  !subdiagonal elements
      REAL(r8), intent(in) :: b(1:n)  !diagonal elements
      REAL(r8), intent(in) :: c(1:n)  !superdiagonal elements
      REAL(r8), intent(in) :: r(1:n)  !right hand side
      REAL(r8), intent(out) :: u(1:n) !solution vector

      INTEGER j
      REAL(r8) gam(1:n),bet
! -----------------------------------------------------------------

      bet = b(1)
      u(1) = r(1) / bet
      DO j = 2, n
            gam(j) = c(j-1) / bet
            bet = b(j) - a(j) * gam(j)
            u(j) = (r(j) - a(j)*u(j-1)) / bet
      ENDDO
      DO j = n-1, 1, -1
            u(j) = u(j) - gam(j+1) * u(j+1)
      ENDDO

      END SUBROUTINE tridia
