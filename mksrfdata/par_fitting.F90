Module par_fitting

implicit none
private

public :: lmder

contains

subroutine lmder ( fcn, m, n, x, fvec, fjac, ldfjac, ftol, xtol, gtol, maxfev, &
  diag, mode, factor, nprint, info, nfev, njev, ipvt, qtf, xdat, npoint, ydat, nptf, phi, isiter,&
  vf_om, vf_sand, vf_gravels)

!*******************************************************************************
!
!! LMDER minimizes M functions in N variables by the Levenberg-Marquardt method 
!  implemented for fitting the SW retention parameters of the Campbell/van Genuchten model
!  and Ke-Sr relationship in Balland V. and P. A. Arp (2005) model.
!
!  Discussion:
!
!    LMDER minimizes the sum of the squares of M nonlinear functions in
!    N variables by a modification of the Levenberg-Marquardt algorithm.
!    The user must provide a subroutine which calculates the functions
!    and the jacobian.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!    Modified by Nan Wei, 2019/06
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, external FCN, the name of the user-supplied subroutine which
!    calculates the functions and the jacobian.  FCN should have the form:
!      subroutine fcn ( m, n, x, fvec, fjac, ldfjac, iflag, xdat, npoint, ydat, nptf, phi, isiter, &
!                       vf_om, vf_sand, vf_gravels)
!      integer ( kind = 4 ) ldfjac
!      integer ( kind = 4 ) n
!      real ( kind = 8 ) fjac(ldfjac,n)
!      real ( kind = 8 ) fvec(m)
!      integer ( kind = 4 ) iflag
!      real ( kind = 8 ) x(n)
!      xdat, npoint, ydat, nptf, phi and isiter are transfered as the coefficients of the fitting functions.
!      vf_om, vf_sand and vf_gravels are only used for fitting Ke-Sr relationship.
!
!    If IFLAG = 0 on input, then FCN is only being called to allow the user
!    to print out the current iterate.
!    If IFLAG = 1 on input, FCN should calculate the functions at X and
!    return this vector in FVEC.
!    If IFLAG = 2 on input, FCN should calculate the jacobian at X and
!    return this matrix in FJAC.
!    To terminate the algorithm, FCN may set IFLAG negative on return.
!
!    Input, integer ( kind = 4 ) M, is the number of functions.
!
!    Input, integer ( kind = 4 ) N, is the number of variables.  
!    N must not exceed M.
!
!    Input/output, real ( kind = 8 ) X(N).  On input, X must contain an initial
!    estimate of the solution vector.  On output X contains the final
!    estimate of the solution vector.
!
!    Output, real ( kind = 8 ) FVEC(M), the functions evaluated at the output X.
!
!    Output, real ( kind = 8 ) FJAC(LDFJAC,N), an M by N array.  The upper
!    N by N submatrix of FJAC contains an upper triangular matrix R with
!    diagonal elements of nonincreasing magnitude such that
!      P' * ( JAC' * JAC ) * P = R' * R,
!    where P is a permutation matrix and JAC is the final calculated jacobian.
!    Column J of P is column IPVT(J) of the identity matrix.  The lower
!    trapezoidal part of FJAC contains information generated during
!    the computation of R.
!
!    Input, integer ( kind = 4 ) LDFJAC, the leading dimension of FJAC.
!    LDFJAC must be at least M.
!
!    Input, real ( kind = 8 ) FTOL.  Termination occurs when both the actual
!    and predicted relative reductions in the sum of squares are at most FTOL.
!    Therefore, FTOL measures the relative error desired in the sum of
!    squares.  FTOL should be nonnegative.
!
!    Input, real ( kind = 8 ) XTOL.  Termination occurs when the relative error
!    between two consecutive iterates is at most XTOL.  XTOL should be
!    nonnegative.
!
!    Input, real ( kind = 8 ) GTOL.  Termination occurs when the cosine of the
!    angle between FVEC and any column of the jacobian is at most GTOL in
!    absolute value.  Therefore, GTOL measures the orthogonality desired
!    between the function vector and the columns of the jacobian.  GTOL should
!    be nonnegative.
!
!    Input, integer ( kind = 4 ) MAXFEV.  Termination occurs when the number of
!    calls to FCN with IFLAG = 1 is at least MAXFEV by the end of an iteration.
!
!    Input/output, real ( kind = 8 ) DIAG(N).  If MODE = 1, then DIAG is set
!    internally.  If MODE = 2, then DIAG must contain positive entries that
!    serve as multiplicative scale factors for the variables.
!
!    Input, integer ( kind = 4 ) MODE, scaling option.
!    1, variables will be scaled internally.
!    2, scaling is specified by the input DIAG vector.
!
!    Input, real ( kind = 8 ) FACTOR, determines the initial step bound.  This
!    bound is set to the product of FACTOR and the euclidean norm of DIAG*X if
!    nonzero, or else to FACTOR itself.  In most cases, FACTOR should lie
!    in the interval (0.1, 100) with 100 the recommended value.
!
!    Input, integer ( kind = 4 ) NPRINT, enables controlled printing of iterates
!    if it is positive.  In this case, FCN is called with IFLAG = 0 at the
!    beginning of the first iteration and every NPRINT iterations thereafter
!    and immediately prior to return, with X and FVEC available
!    for printing.  If NPRINT is not positive, no special calls
!    of FCN with IFLAG = 0 are made.
!
!    Output, integer ( kind = 4 ) INFO, error flag.  If the user has terminated
!    execution, INFO is set to the (negative) value of IFLAG. See description
!    of FCN.  Otherwise, INFO is set as follows:
!    0, improper input parameters.
!    1, both actual and predicted relative reductions in the sum of
!       squares are at most FTOL.
!    2, relative error between two consecutive iterates is at most XTOL.
!    3, conditions for INFO = 1 and INFO = 2 both hold.
!    4, the cosine of the angle between FVEC and any column of the jacobian
!       is at most GTOL in absolute value.
!    5, number of calls to FCN with IFLAG = 1 has reached MAXFEV.
!    6, FTOL is too small.  No further reduction in the sum of squares
!       is possible.
!    7, XTOL is too small.  No further improvement in the approximate
!       solution X is possible.
!    8, GTOL is too small.  FVEC is orthogonal to the columns of the
!       jacobian to machine precision.
!
!    Output, integer ( kind = 4 ) NFEV, the number of calls to FCN with
!    IFLAG = 1.
!
!    Output, integer ( kind = 4 ) NJEV, the number of calls to FCN with
!    IFLAG = 2.
!
!    Output, integer ( kind = 4 ) IPVT(N), defines a permutation matrix P
!    such that JAC*P = Q*R, where JAC is the final calculated jacobian, Q is
!    orthogonal (not stored), and R is upper triangular with diagonal
!    elements of nonincreasing magnitude.  Column J of P is column
!    IPVT(J) of the identity matrix.
!
!    Output, real ( kind = 8 ) QTF(N), contains the first N elements of Q'*FVEC.
!
  implicit none

  integer ( kind = 4 ) ldfjac
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) actred
  real ( kind = 8 ) delta
  real ( kind = 8 ) diag(n)
  real ( kind = 8 ) dirder
  real ( kind = 8 ) epsmch
  real ( kind = 8 ) factor
  external fcn
  real ( kind = 8 ) fjac(ldfjac,n)
  real ( kind = 8 ) fnorm
  real ( kind = 8 ) fnorm1
  real ( kind = 8 ) ftol
  real ( kind = 8 ) fvec(m)
  real ( kind = 8 ) gnorm
  real ( kind = 8 ) gtol
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) maxfev
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) nfev
  integer ( kind = 4 ) njev
  integer ( kind = 4 ) nprint
  real ( kind = 8 ) par
  logical pivot
  real ( kind = 8 ) pnorm
  real ( kind = 8 ) prered
  real ( kind = 8 ) qtf(n)
  real ( kind = 8 ) ratio
  real ( kind = 8 ) sum2
  real ( kind = 8 ) temp
  real ( kind = 8 ) temp1
  real ( kind = 8 ) temp2
  real ( kind = 8 ) wa1(n)
  real ( kind = 8 ) wa2(n)
  real ( kind = 8 ) wa3(n)
  real ( kind = 8 ) wa4(m)
  real ( kind = 8 ) xnorm
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xtol
  real ( kind = 8 ) phi
  integer ( kind = 4 ) isiter
  integer ( kind = 4 ) npoint
  integer ( kind = 4 ) nptf
  real ( kind = 8 ) xdat(npoint)
  real ( kind = 8 ) ydat(nptf,npoint)
  real ( kind = 8 ), optional :: vf_om
  real ( kind = 8 ), optional :: vf_sand
  real ( kind = 8 ), optional :: vf_gravels

  epsmch = epsilon ( epsmch )

  info = 0
  iflag = 0
  nfev = 0
  njev = 0
!
!  Check the input parameters for errors.
!
  if ( n <= 0 ) then
    go to 300
  end if

  if ( m < n ) then
    go to 300
  end if

  if ( ldfjac < m &
    .or. ftol < 0.0D+00 .or. xtol < 0.0D+00 .or. gtol < 0.0D+00 &
     .or. maxfev <= 0 .or. factor <= 0.0D+00 ) then
    go to 300
  end if

  if ( mode == 2 ) then
    do j = 1, n
      if ( diag(j) <= 0.0D+00 ) then
        go to 300
      end if
    end do
  end if
!
!  Evaluate the function at the starting point and calculate its norm.
!
  iflag = 1
  if (present(vf_om) .and. present(vf_sand) .and. present(vf_gravels)) then
      call fcn ( m, n, x, fvec, fjac, ldfjac, iflag, xdat, npoint, ydat, nptf, phi, isiter,&
                 vf_om, vf_sand, vf_gravels)
  else
      call fcn ( m, n, x, fvec, fjac, ldfjac, iflag, xdat, npoint, ydat, nptf, phi, isiter )
  end if
  nfev = 1
  if ( iflag < 0 ) then
    go to 300
  end if

  fnorm = enorm ( m, fvec )
!
!  Initialize Levenberg-Marquardt parameter and iteration counter.
!
  par = 0.0D+00
  iter = 1
!
!  Beginning of the outer loop.
!
  do
!
!  Calculate the jacobian matrix.
!
    iflag = 2
    if (present(vf_om) .and. present(vf_sand) .and. present(vf_gravels)) then
        call fcn ( m, n, x, fvec, fjac, ldfjac, iflag, xdat, npoint, ydat, nptf, phi, isiter,&
                   vf_om, vf_sand, vf_gravels)
    else
        call fcn ( m, n, x, fvec, fjac, ldfjac, iflag, xdat, npoint, ydat, nptf, phi, isiter )
    end if

    njev = njev + 1

    if ( iflag < 0 ) then
      go to 300
    end if
!
!  If requested, call FCN to enable printing of iterates.
!
    if ( 0 < nprint ) then
      iflag = 0
      if ( mod ( iter - 1, nprint ) == 0 ) then
           if (present(vf_om) .and. present(vf_sand) .and. present(vf_gravels)) then
               call fcn ( m, n, x, fvec, fjac, ldfjac, iflag, xdat, npoint, ydat, nptf, phi, isiter,&
                          vf_om, vf_sand, vf_gravels)
           else
               call fcn ( m, n, x, fvec, fjac, ldfjac, iflag, xdat, npoint, ydat, nptf, phi, isiter )
           end if
      end if
      if ( iflag < 0 ) then
        go to 300
      end if
    end if
!
!  Compute the QR factorization of the jacobian.
!
    pivot = .true.
    call qrfac ( m, n, fjac, ldfjac, pivot, ipvt, n, wa1, wa2 )
!
!  On the first iteration and if mode is 1, scale according
!  to the norms of the columns of the initial jacobian.
!
    if ( iter == 1 ) then

      if ( mode /= 2 ) then
        diag(1:n) = wa2(1:n)
        do j = 1, n
          if ( wa2(j) == 0.0D+00 ) then
            diag(j) = 1.0D+00
          end if
        end do
      end if
!
!  On the first iteration, calculate the norm of the scaled X
!  and initialize the step bound DELTA.
!
      wa3(1:n) = diag(1:n) * x(1:n)

      xnorm = enorm ( n, wa3 )

      if ( xnorm == 0.0D+00 ) then
        delta = factor
      else
        delta = factor * xnorm
      end if

    end if
!
!  Form Q'*FVEC and store the first N components in QTF.
!
    wa4(1:m) = fvec(1:m)

    do j = 1, n

      if ( fjac(j,j) /= 0.0D+00 ) then
        sum2 = dot_product ( wa4(j:m), fjac(j:m,j) )
        temp = - sum2 / fjac(j,j)
        wa4(j:m) = wa4(j:m) + fjac(j:m,j) * temp
      end if

      fjac(j,j) = wa1(j)
      qtf(j) = wa4(j)

    end do
!
!  Compute the norm of the scaled gradient.
!
    gnorm = 0.0D+00

    if ( fnorm /= 0.0D+00 ) then

      do j = 1, n
        l = ipvt(j)
        if ( wa2(l) /= 0.0D+00 ) then
          sum2 = dot_product ( qtf(1:j), fjac(1:j,j) ) / fnorm
          gnorm = max ( gnorm, abs ( sum2 / wa2(l) ) )
        end if
      end do

    end if
!
!  Test for convergence of the gradient norm.
!
    if ( gnorm <= gtol ) then
      info = 4
      go to 300
    end if
!
!  Rescale if necessary.
!
    if ( mode /= 2 ) then
      do j = 1, n
        diag(j) = max ( diag(j), wa2(j) )
      end do
    end if
!
!  Beginning of the inner loop.
!
    do
!
!  Determine the Levenberg-Marquardt parameter.
!
      call lmpar ( n, fjac, ldfjac, ipvt, diag, qtf, delta, par, wa1, wa2 )
!
!  Store the direction p and x + p. calculate the norm of p.
!
      wa1(1:n) = - wa1(1:n)
      wa2(1:n) = x(1:n) + wa1(1:n)
      wa3(1:n) = diag(1:n) * wa1(1:n)

      pnorm = enorm ( n, wa3 )
!
!  On the first iteration, adjust the initial step bound.
!
      if ( iter == 1 ) then
        delta = min ( delta, pnorm )
      end if
!
!  Evaluate the function at x + p and calculate its norm.
!
      iflag = 1
      if (present(vf_om) .and. present(vf_sand) .and. present(vf_gravels)) then
          call fcn ( m, n, wa2, wa4, fjac, ldfjac, iflag, xdat, npoint, ydat, nptf, phi, isiter,&
                     vf_om, vf_sand, vf_gravels)
      else
          call fcn ( m, n, wa2, wa4, fjac, ldfjac, iflag, xdat, npoint, ydat, nptf, phi, isiter )
      end if

      nfev = nfev + 1

      if ( iflag < 0 ) then
        go to 300
      end if

      fnorm1 = enorm ( m, wa4 )
!
!  Compute the scaled actual reduction.
!
      if ( 0.1D+00 * fnorm1 < fnorm ) then
        actred = 1.0D+00 - ( fnorm1 / fnorm ) ** 2
      else
        actred = - 1.0D+00
      end if
!
!  Compute the scaled predicted reduction and
!  the scaled directional derivative.
!
      do j = 1, n
        wa3(j) = 0.0D+00
        l = ipvt(j)
        temp = wa1(l)
        wa3(1:j) = wa3(1:j) + fjac(1:j,j) * temp
      end do

      temp1 = enorm ( n, wa3 ) / fnorm
      temp2 = ( sqrt ( par ) * pnorm ) / fnorm
      prered = temp1 ** 2 + temp2 ** 2 / 0.5D+00
      dirder = - ( temp1 ** 2 + temp2 ** 2 )
!
!  Compute the ratio of the actual to the predicted reduction.
!
      if ( prered /= 0.0D+00 ) then
        ratio = actred / prered
      else
        ratio = 0.0D+00
      end if
!
!  Update the step bound.
!
      if ( ratio <= 0.25D+00 ) then

        if ( 0.0D+00 <= actred ) then
          temp = 0.5D+00
        end if

        if ( actred < 0.0D+00 ) then
          temp = 0.5D+00 * dirder / ( dirder + 0.5D+00 * actred )
        end if

        if ( 0.1D+00 * fnorm1 >= fnorm .or. temp < 0.1D+00 ) then
          temp = 0.1D+00
        end if

        delta = temp * min ( delta, pnorm / 0.1D+00 )
        par = par / temp

      else

        if ( par == 0.0D+00 .or. ratio >= 0.75D+00 ) then
          delta = 2.0D+00 * pnorm
          par = 0.5D+00 * par
        end if

      end if
!
!  Successful iteration.
!
!  Update X, FVEC, and their norms.
!
      if ( 0.0001D+00 <= ratio ) then
        x(1:n) = wa2(1:n)
        wa2(1:n) = diag(1:n) * x(1:n)
        fvec(1:m) = wa4(1:m)
        xnorm = enorm ( n, wa2 )
        fnorm = fnorm1
        iter = iter + 1
      end if
!
!  Tests for convergence.
!
      if ( abs ( actred) <= ftol .and. &
        prered <= ftol .and. &
        0.5D+00 * ratio <= 1.0D+00 ) then
        info = 1
      end if

      if ( delta <= xtol * xnorm ) then
        info = 2
      end if

      if ( abs ( actred) <= ftol .and. prered <= ftol &
        .and. 0.5D+00 * ratio <= 1.0D+00 .and. info == 2 ) then
        info = 3
      end if

      if ( info /= 0 ) then
        go to 300
      end if
!
!  Tests for termination and stringent tolerances.
!
      if ( nfev >= maxfev ) then
        info = 5
      end if

      if ( abs ( actred ) <= epsmch .and. prered <= epsmch &
        .and. 0.5D+00 * ratio <= 1.0D+00 ) then
        info = 6
      end if

      if ( delta <= epsmch * xnorm ) then
        info = 7
      end if

      if ( gnorm <= epsmch ) then
        info = 8
      end if

      if ( info /= 0 ) then
        go to 300
      end if
!
!  End of the inner loop. repeat if iteration unsuccessful.
!
      if ( 0.0001D+00 <= ratio ) then
        exit
      end if

    end do
!
!  End of the outer loop.
!
  end do

  300 continue
!
!  Termination, either normal or user imposed.
!
  if ( iflag < 0 ) then
    info = iflag
  end if

  iflag = 0

  if ( 0 < nprint ) then
       if (present(vf_om) .and. present(vf_sand) .and. present(vf_gravels)) then
           call fcn ( m, n, x, fvec, fjac, ldfjac, iflag, xdat, npoint, ydat, nptf, phi, isiter,&
                      vf_om, vf_sand, vf_gravels)
       else
           call fcn ( m, n, x, fvec, fjac, ldfjac, iflag, xdat, npoint, ydat, nptf, phi, isiter )
       end if
  end if

  return
end subroutine lmder

subroutine lmpar ( n, r, ldr, ipvt, diag, qtb, delta, par, x, sdiag )

!*****************************************************************************80
!
!! LMPAR computes a parameter for the Levenberg-Marquardt method.
!
!  Discussion:
!
!    Given an M by N matrix A, an N by N nonsingular diagonal
!    matrix D, an M-vector B, and a positive number DELTA,
!    the problem is to determine a value for the parameter
!    PAR such that if X solves the system
!
!      A*X = B,
!      sqrt ( PAR ) * D * X = 0,
!
!    in the least squares sense, and DXNORM is the euclidean
!    norm of D*X, then either PAR is zero and
!
!      ( DXNORM - DELTA ) <= 0.1 * DELTA,
!
!    or PAR is positive and
!
!      abs ( DXNORM - DELTA) <= 0.1 * DELTA.
!
!    This function completes the solution of the problem
!    if it is provided with the necessary information from the
!    QR factorization, with column pivoting, of A.  That is, if
!    A*P = Q*R, where P is a permutation matrix, Q has orthogonal
!    columns, and R is an upper triangular matrix with diagonal
!    elements of nonincreasing magnitude, then LMPAR expects
!    the full upper triangle of R, the permutation matrix P,
!    and the first N components of Q'*B.  On output
!    LMPAR also provides an upper triangular matrix S such that
!
!      P' * ( A' * A + PAR * D * D ) * P = S'* S.
!
!    S is employed within LMPAR and may be of separate interest.
!
!    Only a few iterations are generally needed for convergence
!    of the algorithm.  
!
!    If, however, the limit of 10 iterations is reached, then the output 
!    PAR will contain the best value obtained so far.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    24 January 2014
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of R.
!
!    Input/output, real ( kind = 8 ) R(LDR,N),the N by N matrix.  The full
!    upper triangle must contain the full upper triangle of the matrix R.
!    On output the full upper triangle is unaltered, and the strict lower
!    triangle contains the strict upper triangle (transposed) of the upper
!    triangular matrix S.
!
!    Input, integer ( kind = 4 ) LDR, the leading dimension of R.  LDR must be
!    no less than N.
!
!    Input, integer ( kind = 4 ) IPVT(N), defines the permutation matrix P 
!    such that A*P = Q*R.  Column J of P is column IPVT(J) of the 
!    identity matrix.
!
!    Input, real ( kind = 8 ) DIAG(N), the diagonal elements of the matrix D.
!
!    Input, real ( kind = 8 ) QTB(N), the first N elements of the vector Q'*B.
!
!    Input, real ( kind = 8 ) DELTA, an upper bound on the euclidean norm
!    of D*X.  DELTA should be positive.
!
!    Input/output, real ( kind = 8 ) PAR.  On input an initial estimate of the
!    Levenberg-Marquardt parameter.  On output the final estimate.
!    PAR should be nonnegative.
!
!    Output, real ( kind = 8 ) X(N), the least squares solution of the system
!    A*X = B, sqrt(PAR)*D*X = 0, for the output value of PAR.
!
!    Output, real ( kind = 8 ) SDIAG(N), the diagonal elements of the upper
!    triangular matrix S.
!
  implicit none

  integer ( kind = 4 ) ldr
  integer ( kind = 4 ) n

  real ( kind = 8 ) delta
  real ( kind = 8 ) diag(n)
  real ( kind = 8 ) dwarf
  real ( kind = 8 ) dxnorm
  real ( kind = 8 ) gnorm
  real ( kind = 8 ) fp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) nsing
  real ( kind = 8 ) par
  real ( kind = 8 ) parc
  real ( kind = 8 ) parl
  real ( kind = 8 ) paru
  real ( kind = 8 ) qnorm
  real ( kind = 8 ) qtb(n)
  real ( kind = 8 ) r(ldr,n)
  real ( kind = 8 ) sdiag(n)
  real ( kind = 8 ) sum2
  real ( kind = 8 ) temp
  real ( kind = 8 ) wa1(n)
  real ( kind = 8 ) wa2(n)
  real ( kind = 8 ) x(n)
!
!  DWARF is the smallest positive magnitude.
!
  dwarf = tiny ( dwarf )
!
!  Compute and store in X the Gauss-Newton direction.
!
!  If the jacobian is rank-deficient, obtain a least squares solution.
!
  nsing = n

  do j = 1, n
    wa1(j) = qtb(j)
    if ( r(j,j) == 0.0D+00 .and. nsing == n ) then
      nsing = j - 1
    end if
    if ( nsing < n ) then
      wa1(j) = 0.0D+00
    end if
  end do

  do k = 1, nsing
    j = nsing - k + 1
    wa1(j) = wa1(j) / r(j,j)
    temp = wa1(j)
    wa1(1:j-1) = wa1(1:j-1) - r(1:j-1,j) * temp
  end do

  do j = 1, n
    l = ipvt(j)
    x(l) = wa1(j)
  end do
!
!  Initialize the iteration counter.
!  Evaluate the function at the origin, and test
!  for acceptance of the Gauss-Newton direction.
!
  iter = 0
  wa2(1:n) = diag(1:n) * x(1:n)
  dxnorm = enorm ( n, wa2 )
  fp = dxnorm - delta

  if ( fp <= 0.1D+00 * delta ) then
    if ( iter == 0 ) then
      par = 0.0D+00
    end if
    return
  end if
!
!  If the jacobian is not rank deficient, the Newton
!  step provides a lower bound, PARL, for the zero of
!  the function.
!
!  Otherwise set this bound to zero.
!
  parl = 0.0D+00

  if ( n <= nsing ) then

    do j = 1, n
      l = ipvt(j)
      wa1(j) = diag(l) * ( wa2(l) / dxnorm )
    end do

    do j = 1, n
      sum2 = dot_product ( wa1(1:j-1), r(1:j-1,j) )
      wa1(j) = ( wa1(j) - sum2 ) / r(j,j)
    end do

    temp = enorm ( n, wa1 )
    parl = ( ( fp / delta ) / temp ) / temp

  end if
!
!  Calculate an upper bound, PARU, for the zero of the function.
!
  do j = 1, n
    sum2 = dot_product ( qtb(1:j), r(1:j,j) )
    l = ipvt(j)
    wa1(j) = sum2 / diag(l)
  end do

  gnorm = enorm ( n, wa1 )
  paru = gnorm / delta

  if ( paru == 0.0D+00 ) then
    paru = dwarf / min ( delta, 0.1D+00 )
  end if
!
!  If the input PAR lies outside of the interval (PARL, PARU),
!  set PAR to the closer endpoint.
!
  par = max ( par, parl )
  par = min ( par, paru )
  if ( par == 0.0D+00 ) then
    par = gnorm / dxnorm
  end if
!
!  Beginning of an iteration.
!
  do
 
    iter = iter + 1
!
!  Evaluate the function at the current value of PAR.
!
    if ( par == 0.0D+00 ) then
      par = max ( dwarf, 0.001D+00 * paru )
    end if

    wa1(1:n) = sqrt ( par ) * diag(1:n)

    call qrsolv ( n, r, ldr, ipvt, wa1, qtb, x, sdiag )

    wa2(1:n) = diag(1:n) * x(1:n)
    dxnorm = enorm ( n, wa2 )
    temp = fp
    fp = dxnorm - delta
!
!  If the function is small enough, accept the current value of PAR.
!
    if ( abs ( fp ) <= 0.1D+00 * delta ) then
      exit
    end if
!
!  Test for the exceptional cases where PARL
!  is zero or the number of iterations has reached 10.
!
    if ( parl == 0.0D+00 .and. fp <= temp .and. temp < 0.0D+00 ) then
      exit
    else if ( iter == 10 ) then
      exit
    end if
!
!  Compute the Newton correction.
!
    do j = 1, n
      l = ipvt(j)
      wa1(j) = diag(l) * ( wa2(l) / dxnorm )
    end do

    do j = 1, n
      wa1(j) = wa1(j) / sdiag(j)
      temp = wa1(j)
      wa1(j+1:n) = wa1(j+1:n) - r(j+1:n,j) * temp
    end do

    temp = enorm ( n, wa1 )
    parc = ( ( fp / delta ) / temp ) / temp
!
!  Depending on the sign of the function, update PARL or PARU.
!
    if ( 0.0D+00 < fp ) then
      parl = max ( parl, par )
    else if ( fp < 0.0D+00 ) then
      paru = min ( paru, par )
    end if
!
!  Compute an improved estimate for PAR.
!
    par = max ( parl, par + parc )
!
!  End of an iteration.
!
  end do
!
!  Termination.
!
  if ( iter == 0 ) then
    par = 0.0D+00
  end if

  return
end subroutine lmpar

subroutine qrfac ( m, n, a, lda, pivot, ipvt, lipvt, rdiag, acnorm )

!*****************************************************************************80
!
!! QRFAC computes a QR factorization using Householder transformations.
!
!  Discussion:
!
!    This function uses Householder transformations with optional column
!    pivoting to compute a QR factorization of the
!    M by N matrix A.  That is, QRFAC determines an orthogonal
!    matrix Q, a permutation matrix P, and an upper trapezoidal
!    matrix R with diagonal elements of nonincreasing magnitude,
!    such that A*P = Q*R.  
!
!    The Householder transformation for column K, K = 1,2,...,min(M,N), 
!    is of the form
!
!      I - ( 1 / U(K) ) * U * U'
!
!    where U has zeros in the first K-1 positions.  
!
!    The form of this transformation and the method of pivoting first
!    appeared in the corresponding LINPACK routine.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows of A.
!
!    Input, integer ( kind = 4 ) N, the number of columns of A.
!
!    Input/output, real ( kind = 8 ) A(LDA,N), the M by N array.
!    On input, A contains the matrix for which the QR factorization is to
!    be computed.  On output, the strict upper trapezoidal part of A contains
!    the strict upper trapezoidal part of R, and the lower trapezoidal
!    part of A contains a factored form of Q, the non-trivial elements of
!    the U vectors described above.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A, which must
!    be no less than M.
!
!    Input, logical PIVOT, is TRUE if column pivoting is to be carried out.
!
!    Output, integer ( kind = 4 ) IPVT(LIPVT), defines the permutation matrix P 
!    such that A*P = Q*R.  Column J of P is column IPVT(J) of the identity 
!    matrix.  If PIVOT is false, IPVT is not referenced.
!
!    Input, integer ( kind = 4 ) LIPVT, the dimension of IPVT, which should 
!    be N if pivoting is used.
!
!    Output, real ( kind = 8 ) RDIAG(N), contains the diagonal elements of R.
!
!    Output, real ( kind = 8 ) ACNORM(N), the norms of the corresponding
!    columns of the input matrix A.  If this information is not needed,
!    then ACNORM can coincide with RDIAG.
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) lipvt
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) acnorm(n)
  real ( kind = 8 ) ajnorm
  real ( kind = 8 ) epsmch
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_temp
  integer ( kind = 4 ) ipvt(lipvt)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ) minmn
  logical pivot
  real ( kind = 8 ) r8_temp(m)
  real ( kind = 8 ) rdiag(n)
  real ( kind = 8 ) temp
  real ( kind = 8 ) wa(n)

  epsmch = epsilon ( epsmch )
!
!  Compute the initial column norms and initialize several arrays.
!
  do j = 1, n
    acnorm(j) = enorm ( m, a(1:m,j) )
  end do

  rdiag(1:n) = acnorm(1:n)
  wa(1:n) = acnorm(1:n)

  if ( pivot ) then
    do j = 1, n
      ipvt(j) = j
    end do
  end if
!
!  Reduce A to R with Householder transformations.
!
  minmn = min ( m, n )

  do j = 1, minmn
!
!  Bring the column of largest norm into the pivot position.
!
    if ( pivot ) then

      kmax = j

      do k = j, n
        if ( rdiag(kmax) < rdiag(k) ) then
          kmax = k
        end if
      end do

      if ( kmax /= j ) then

        r8_temp(1:m) = a(1:m,j)
        a(1:m,j)     = a(1:m,kmax)
        a(1:m,kmax)  = r8_temp(1:m)

        rdiag(kmax) = rdiag(j)
        wa(kmax) = wa(j)

        i4_temp    = ipvt(j)
        ipvt(j)    = ipvt(kmax)
        ipvt(kmax) = i4_temp

      end if

    end if
!
!  Compute the Householder transformation to reduce the
!  J-th column of A to a multiple of the J-th unit vector.
!
    ajnorm = enorm ( m-j+1, a(j,j) )

    if ( ajnorm /= 0.0D+00 ) then

      if ( a(j,j) < 0.0D+00 ) then
        ajnorm = -ajnorm
      end if

      a(j:m,j) = a(j:m,j) / ajnorm
      a(j,j) = a(j,j) + 1.0D+00
!
!  Apply the transformation to the remaining columns and update the norms.
!
      do k = j + 1, n

        temp = dot_product ( a(j:m,j), a(j:m,k) ) / a(j,j)

        a(j:m,k) = a(j:m,k) - temp * a(j:m,j)

        if ( pivot .and. rdiag(k) /= 0.0D+00 ) then

          temp = a(j,k) / rdiag(k)
          rdiag(k) = rdiag(k) * sqrt ( max ( 0.0D+00, 1.0D+00-temp ** 2 ) )

          if ( 0.05D+00 * ( rdiag(k) / wa(k) ) ** 2 <= epsmch ) then
            rdiag(k) = enorm ( m-j, a(j+1,k) )
            wa(k) = rdiag(k)
          end if

        end if

      end do

    end if

    rdiag(j) = - ajnorm

  end do

  return
end subroutine qrfac

subroutine qrsolv ( n, r, ldr, ipvt, diag, qtb, x, sdiag )

!*****************************************************************************80
!
!! QRSOLV solves a rectangular linear system A*x=b in the least squares sense.
!
!  Discussion:
!
!    Given an M by N matrix A, an N by N diagonal matrix D,
!    and an M-vector B, the problem is to determine an X which
!    solves the system
!
!      A*X = B
!      D*X = 0
!
!    in the least squares sense.
!
!    This function completes the solution of the problem
!    if it is provided with the necessary information from the
!    QR factorization, with column pivoting, of A.  That is, if
!    A*P = Q*R, where P is a permutation matrix, Q has orthogonal
!    columns, and R is an upper triangular matrix with diagonal
!    elements of nonincreasing magnitude, then QRSOLV expects
!    the full upper triangle of R, the permutation matrix p,
!    and the first N components of Q'*B.
!
!    The system is then equivalent to
!
!      R*Z = Q'*B
!      P'*D*P*Z = 0
!
!    where X = P*Z.  If this system does not have full rank,
!    then a least squares solution is obtained.  On output QRSOLV
!    also provides an upper triangular matrix S such that
!
!      P'*(A'*A + D*D)*P = S'*S.
!
!    S is computed within QRSOLV and may be of separate interest.
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of R.
!
!    Input/output, real ( kind = 8 ) R(LDR,N), the N by N matrix.
!    On input the full upper triangle must contain the full upper triangle
!    of the matrix R.  On output the full upper triangle is unaltered, and
!    the strict lower triangle contains the strict upper triangle
!    (transposed) of the upper triangular matrix S.
!
!    Input, integer ( kind = 4 ) LDR, the leading dimension of R, which must be
!    at least N.
!
!    Input, integer ( kind = 4 ) IPVT(N), defines the permutation matrix P such 
!    that A*P = Q*R.  Column J of P is column IPVT(J) of the identity matrix.
!
!    Input, real ( kind = 8 ) DIAG(N), the diagonal elements of the matrix D.
!
!    Input, real ( kind = 8 ) QTB(N), the first N elements of the vector Q'*B.
!
!    Output, real ( kind = 8 ) X(N), the least squares solution.
!
!    Output, real ( kind = 8 ) SDIAG(N), the diagonal elements of the upper
!    triangular matrix S.
!
  implicit none

  integer ( kind = 4 ) ldr
  integer ( kind = 4 ) n

  real ( kind = 8 ) c
  real ( kind = 8 ) cotan
  real ( kind = 8 ) diag(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) nsing
  real ( kind = 8 ) qtb(n)
  real ( kind = 8 ) qtbpj
  real ( kind = 8 ) r(ldr,n)
  real ( kind = 8 ) s
  real ( kind = 8 ) sdiag(n)
  real ( kind = 8 ) sum2
  real ( kind = 8 ) t
  real ( kind = 8 ) temp
  real ( kind = 8 ) wa(n)
  real ( kind = 8 ) x(n)
!
!  Copy R and Q'*B to preserve input and initialize S.
!
!  In particular, save the diagonal elements of R in X.
!
  do j = 1, n
    r(j:n,j) = r(j,j:n)
    x(j) = r(j,j)
  end do

  wa(1:n) = qtb(1:n)
!
!  Eliminate the diagonal matrix D using a Givens rotation.
!
  do j = 1, n
!
!  Prepare the row of D to be eliminated, locating the
!  diagonal element using P from the QR factorization.
!
    l = ipvt(j)

    if ( diag(l) /= 0.0D+00 ) then

      sdiag(j:n) = 0.0D+00
      sdiag(j) = diag(l)
!
!  The transformations to eliminate the row of D
!  modify only a single element of Q'*B
!  beyond the first N, which is initially zero.
!
      qtbpj = 0.0D+00

      do k = j, n
!
!  Determine a Givens rotation which eliminates the
!  appropriate element in the current row of D.
!
        if ( sdiag(k) /= 0.0D+00 ) then

          if ( abs ( r(k,k) ) < abs ( sdiag(k) ) ) then
            cotan = r(k,k) / sdiag(k)
            s = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * cotan ** 2 )
            c = s * cotan
          else
            t = sdiag(k) / r(k,k)
            c = 0.5D+00 / sqrt ( 0.25D+00 + 0.25D+00 * t ** 2 )
            s = c * t
          end if
!
!  Compute the modified diagonal element of R and
!  the modified element of (Q'*B,0).
!
          r(k,k) = c * r(k,k) + s * sdiag(k)
          temp = c * wa(k) + s * qtbpj
          qtbpj = - s * wa(k) + c * qtbpj
          wa(k) = temp
!
!  Accumulate the tranformation in the row of S.
!
          do i = k + 1, n
            temp = c * r(i,k) + s * sdiag(i)
            sdiag(i) = - s * r(i,k) + c * sdiag(i)
            r(i,k) = temp
          end do

        end if

      end do

    end if
!
!  Store the diagonal element of S and restore
!  the corresponding diagonal element of R.
!
    sdiag(j) = r(j,j)
    r(j,j) = x(j)

  end do
!
!  Solve the triangular system for Z.  If the system is
!  singular, then obtain a least squares solution.
!
  nsing = n

  do j = 1, n

    if ( sdiag(j) == 0.0D+00 .and. nsing == n ) then
      nsing = j - 1
    end if

    if ( nsing < n ) then
      wa(j) = 0.0D+00
    end if

  end do

  do j = nsing, 1, -1
    sum2 = dot_product ( wa(j+1:nsing), r(j+1:nsing,j) )
    wa(j) = ( wa(j) - sum2 ) / sdiag(j)
  end do
!
!  Permute the components of Z back to components of X.
!
  do j = 1, n
    l = ipvt(j)
    x(l) = wa(j)
  end do

  return
end subroutine qrsolv

function enorm ( n, x )

!*****************************************************************************80
!
!! ENORM computes the Euclidean norm of a vector.
!
!  Discussion:
!
!    This is an extremely simplified version of the original ENORM
!    routine, which has been renamed to "ENORM2".
!
!  Licensing:
!
!    This code may freely be copied, modified, and used for any purpose.
!
!  Modified:
!
!    06 April 2010
!
!  Author:
!
!    Original FORTRAN77 version by Jorge More, Burton Garbow, Kenneth Hillstrom.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jorge More, Burton Garbow, Kenneth Hillstrom,
!    User Guide for MINPACK-1,
!    Technical Report ANL-80-74,
!    Argonne National Laboratory, 1980.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, is the length of the vector.
!
!    Input, real ( kind = 8 ) X(N), the vector whose norm is desired.
!
!    Output, real ( kind = 8 ) ENORM, the Euclidean norm of the vector.
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) enorm

  enorm = sqrt ( sum ( x(1:n) ** 2 ))

  return
end function enorm

end Module
