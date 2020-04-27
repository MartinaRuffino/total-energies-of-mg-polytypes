      subroutine l2a(m, n, m1, l, a, b, w, tol, mm, nn,
     * n1, ipivot, x, res, r, q, c, ifault)
c ** purpose **
c subroutine l2a computes least squares solutions to overdetermined
c systems of linear equations.  the method used is a modified
c gram-schmidt orthogonal decomposition with iterative refinement of
c the solution.  the solution may be subject to linear equality
c constraints.  output includes the least squares coefficients,
c residuals, unscaled covariance matrix, and information on the
c behavior of the iterative refinement procedure.
c matrix a is the given matrix of a system of m linear equations in n
c unknowns, and matrix w is a given diagonal matrix of weights with all
c diagonal elements nonnegative.  let h = (sqrt(w))*a.
c in the event that n1 (the computed rank of matrix h) is less than n
c (the number of unknown coefficients), the original matrix h (m by n)
c is replaced by a smaller matrix (m by n1) whose columns are linearly
c independent, and a solution is sought for the smaller system of
c equations.  thus n - n1 columns of the original matrix h are deleted,
c and coefficients corresponding to these n - n1 deleted columns will
c be set equal to zero.
c
c ** input variables **
c m      total number of equations.
c n      number of unknown coefficients.
c m1     number of linear constraints (0 <= m1 <= m and m1 <= n).
c l      number of right-hand sides (vectors of observations).
c a      two-dimensional array of size (mm,n).  on entry, the array a
c        contains the given matrix of a system of m linear equations
c        in n unknowns, where the first m1 equations are to be
c        satisfied exactly.  a is left intact on exit.
c b      two-dimensional array of size (mm,l).  on entry, b contains
c        the l given right-hand sides (vectors of observations).  b is
c        left intact on exit.
c w      vector of size m.  on entry, w contains the diagonal elements
c        of a given diagonal matrix of weights, all nonnegative.
c        (the first m1 elements of w are set equal to 1.0 by the
c        program when m1 is greater than zero.)  on exit, the original
c        elements of w have been replaced by their square roots.
c tol    parameter used in determining the rank of matrix h.
c        note --
c        (1) if tol equals zero, the tolerance used in subroutine
c            decom1 will be based on machine precision.
c        (2) if tol is greater than zero, this value of tol will be
c            used in setting an absolute tolerance for comparison with
c            diagonal elements of the triangular matrix obtained in
c            subroutine decom1.  the value of tol can be based on
c            knowledge concerning the accuracy of the data.
c mm     dimensioning parameter specifying maximum number of rows in
c        the arrays a, b, res and q.  mm must satisfy mm >= m.
c nn     dimensioning parameter specifying maximum number of rows in
c        the arrays x and r.  nn must satisfy nn >= n.
c
c ** output variables and internal variables **
c n1     computed rank of matrix h, where h = (sqrt(w))*a.
c ipivot vector of size n.  on exit, this array records the order
c        in which the columns of h were selected by the pivoting
c        scheme in the course of the orthogonal decomposition.
c        whenever n1 < n, the first n1 elements of ipivot indicate
c        which columns of h were found to be linearly independent.
c x      two-dimensional array of size (nn,l).  on exit, x contains
c        the solution vectors.
c res    two-dimensional array of size (mm,l).  on exit, res contains
c        the residual vectors.
c r      two-dimensional array of size (nn,n).  on exit, r contains
c        the lower triangular portion of the symmetric unscaled
c        covariance matrix.  (this array is used internally to store
c        results from subroutine decom1 which are destroyed in
c        computing the covariance matrix.)
c q      two-dimensional array of size (mm,n) used internally only.
c c      vector having at least 4*(m+n)+2*l elements used (1) for
c        internal work space and (2) for returning information on the
c        behavior of the iterative refinement procedure.
c        (a) numit is the number of iterations carried out during the
c            iterative refinement in attempting to obtain a solution
c            for the k-th right-hand side.
c            on exit, c(k) = +numit if the solution converged, and
c                     c(k) = -numit if the solution failed to converge.
c        (b) digitx gives an estimate of the number of correct digits
c            in the initial solution of the coefficients for the k-th
c            right-hand side.  on exit, c(k+l) = digitx.
c ifault fault indicator which is zero if no errors were encountered
c        and positive if errors were detected or if evidence of severe
c        ill-conditioning was found.  diagnostic messages are printed
c        from subroutine error.  if ifault is set equal to 1, 2, 3, 4,
c        5, 6 or 7, execution is terminated.  execution continues when
c        ifault is set equal to 8, 9 or 10 provided that a solution
c        was obtained for at least one right-hand side.  the value of
c        ifault is used to indicate the following --
c        0 = no errors encountered.
c        1 = bad input parameter (m, n or l).
c        2 = bad input parameter (m1).
c        3 = bad dimension.  either m > mm or n > nn.
c        4 = at least one weight is negative.
c        5 = either matrix h or matrix of constraints equals zero.
c        6 = constraints are linearly dependent.
c        7 = all solutions failed to converge.
c        8 = solution failed to converge for at least one right-hand
c            side.
c        9 = large number of iterations required for convergence.
c       10 = estimated number of digits in initial solution of
c            coefficients is small.
c       11 = diagonal element of covariance matrix was computed to be
c            negative owing to rounding error.
c
c ** subroutines required **
c subroutine decom1
c        uses modified gram-schmidt algorithm with pivoting to
c        obtain an orthogonal decomposition of the input matrix.
c subroutine solve1
c        computes coefficients and residuals.  iterative refinement is
c        used to improve the accuracy of the initial solution.
c subroutine covar
c        computes unscaled covariance matrix of the coefficients.
c subroutine error
c        prints error diagnostics when errors are detected or when
c        evidence of severe ill-conditioning is found.
c
c ** storage requirements **
c the storage required for the dimensioned arrays in subroutine l2a is
c   m*(2*n + 2*l + 5) + n*(n + l + 5) + 2*l
c locations.  all arrays required in subroutines called by l2a are
c declared herein and are transmitted only through parameter lists of
c call-sequences.
c
c ** precision of arithmetic calculations **
c single precision arithmetic is used for all calculations except the
c double precision accumulation of inner products.  (the variable sum
c is declared to be double precision in subroutines decom1, solve1 and
c covar.)  it is essential for the success of the iterative refinement
c procedure that inner products be accumulated in double precision.
c
c ** conversion of the program to double precision **
c *********************************************************************
c * on computers having short word length (as the ibm 360/370) it may *
c * be desirable to perform all calculations in double precision.  in *
c * this case, the iterative refinement presently included in solve1  *
c * should be omitted.                                                *
c * to convert the program to double precision, the following         *
c * approach is suggested.                                            *
c *                                                                   *
c * 1. variables presently declared to be real should be declared     *
c *    double precision.  those typed integer, double precision and   *
c *    logical should not be changed.                                 *
c * 2. the use of fail, numit and digitx should be omitted.           *
c * 3. description of variable c (at l2a 700-800) should read --      *
c *    c  vector having at least 4*(m+n) elements used only for       *
c *       internal work space.                                        *
c * 4. the value of eta (at l2a 1930) should be set so that it is the *
c *    smallest positive double precision number such that 1.0 + eta  *
c *    is greater than 1.0 in double precision arithmetic.            *
c *    for ibm computer type, eta = 16.**(-13)                        *
c *    for univac computer type, eta = 2.**(-59)                      *
c * 5. the following fortran functions should be changed --           *
c *      single precision name     double precision name              *
c *             dble(x)                x                              *
c *             float(n)               dble(float(n))                 *
c *             sqrt(x)                dsqrt(x)                       *
c *    dble(x) is used in subroutines decom1, solve1 and covar.       *
c *    float(n) is used in subroutine decom1.                         *
c *    sqrt(x) is used in subroutine l2a.                             *
c * 6. it may be necessary or desirable to change certain formats in  *
c *    subroutine error, replacing g specifications by d.             *
c * 7. replace statement l2a 2450 by a statement reading              *
c *          k3 = 1                                                   *
c * 8. further details are given in subroutine solve1 in connection   *
c *    with the omission of iterative refinement.                     *
c * 9. in subroutine l2a, statements l2a 960-1010, 1790-1800, 1990,   *
c *    2320-2330, 2430-2440, 2970, 3170-3460 and 3480-3510 should be  *
c *    omitted.                                                       *
c *    statement numbers given here refer to those in the right-hand  *
c *    margin.  certain comments in subroutine l2a do not apply to    *
c *    the double precision version.                                  *
c *                                                                   *
c *********************************************************************
c
      integer ipivot(n)
      double precision a(mm,n), b(mm,l), c(1), eta, q(mm,n), r(nn,n),
     * res(mm,l), tol, w(m), x(nn,l), z, d1mach
C      double precision digitx
C      logical fail
      logical sing
c
c set value of eta, a machine-dependent parameter.
c eta, the relative machine precision, is the smallest positive real
c number such that 1.0 + eta is greater than 1.0 in floating-point
c arithmetic.
c
c for ibm computer type, eta = 16.**(-5)
c for univac computer type, eta = 2.**(-26)
c for cdc computer type, eta = 2.**(-47)
c for honeywell computer type, eta = 2.**(-27)
c
      eta = d1mach(3)
c
c default value for tol is zero.
c
      if (tol < 0.0d0) tol = 0.0d0
      ifault = 0
C      ksum = 0
c
c perform initial checking of input parameters, dimensions and
c weights for possible errors.
c
      if (m > 0 .and. n > 0 .and. l > 0) go to 10
      ifault = 1
      call error(ifault, k, z)
      return
   10 if (m1 <= m .and. m1 <= n .and. m1 >= 0) go to 20
      ifault = 2
      call error(ifault, k, z)
      return
   20 if (m <= mm .and. n <= nn) go to 30
      ifault = 3
      k = 1
      call error(ifault, k, z)
      return
   30 do 40 i=1,m
        if (m1 > 0 .and. i <= m1) w(i) = 1.0d0
        if (w(i) >= 0.0d0) go to 40
        ifault = 4
        z = w(i)
        call error(ifault, i, z)
   40 continue
      if (ifault == 4) return
      do 50 i=1,m
        w(i) = dsqrt(w(i))
   50 continue
c
c set parameters which allocate vector c to contain certain final
c results and also to be used as work space.
c
Cc k1 is starting point for numit and fail, of length l.
Cc k2 is starting point for digitx, of length l.
c k3 is starting point for d, of length n.
c k4 is starting point for k-th column of b, of length m.
c k5 is starting point for k-th column of x, of length n.
c k6 is starting point for k-th column of res, of length m.
c k7 is starting point for work space of length m.
c k8 is starting point for work space of length m.
c k9 is starting point for work space of length n.
c k10 is starting point for work space of length n.
c
C      k1 = 1
C      k2 = k1 + l
c      k3 = k2 + l
      k3 = 1
      k4 = k3 + n
      k5 = k4 + m
      k6 = k5 + n
      k7 = k6 + m
      k8 = k7 + m
      k9 = k8 + m
      k10 = k9 + n
      k = k10 + n - 1
c
c multiply each row of matrix a (m by n) by its appropriate weight and
c store the result in array q.  set arrays c and r equal to zero.
c
      do 60 i=1,k
        c(i) = 0.0d0
   60 continue
      do 90 j=1,n
        do 70 i=1,m
          q(i,j) = a(i,j)*w(i)
   70   continue
        do 80 i=1,n
          r(i,j) = 0.0d0
   80   continue
   90 continue
c
c obtain an orthogonal decomposition of the matrix q and compute its
c rank.
c
      call decom1(m, n, m1, eta, tol, q, r, c(k3), n1, ipivot, sing,
     * mm, nn)
c
      if (.not.sing) go to 110
      if (n1 > 0) go to 100
      ifault = 5
      call error(ifault, k, z)
      return
  100 ifault = 6
      call error(ifault, k, z)
      return
c
c seek a solution (coefficients and residuals) for each of the l least
c squares problems whose right-hand sides are given in the array b.
c
  110 do 200 k=1,l
c k-th right-hand side.
        k0 = k4 - 1
        do 120 i=1,m
          k0 = k0 + 1
          c(k0) = b(i,k)
  120   continue
c
        call solve1(m, n, m1, a, c(k4), w, n1, ipivot, q, r, c(k3),
c     *   eta, fail, numit, digitx,
     *   c(k5), c(k6), c(k7), c(k8), c(k9), c(k10), mm, nn)
c
        k0 = k5 - 1
        do 130 j=1,n
          k0 = k0 + 1
          x(j,k) = c(k0)
  130   continue
        if (m1 == 0) go to 150
        do 140 i=1,m1
          res(i,k) = 0.0d0
  140   continue
  150   m1p1 = m1 + 1
        if (m1p1 > m) go to 170
        k0 = k6 + m1 - 1
        do 160 i=m1p1,m
          k0 = k0 + 1
          res(i,k) = c(k0)
  160   continue
  170   continue
c
Cc for right-hand sides where convergence of a solution is reported,
Cc a check is made for evidence of severe ill-conditioning.  such
Cc evidence is furnished by large values of numit (number of iterations
Cc before convergence was obtained) and small values of digitx
Cc (estimate of the number of correct digits in the initial solution
Cc of the coefficients).  if numit exceeds -alog10(eta) a diagnostic
Cc message is printed to warn of ill-conditioning.  if digitx is less
Cc than 0.5 (half a decimal digit) a similar diagnostic message is
Cc printed.
Cc
C        c(k) = dble(numit)
C        if (fail) c(k) = -c(k)
C        k0 = k2 + k - 1
C        c(k0) = digitx
C        if (.not.fail) go to 180
Cc ksum is a tally of solutions which failed to converge.
C        ksum = ksum + 1
C        ifault = 8
C        call error(ifault, k, z)
C        go to 200
C  180   z = -dlog10(eta)
C        if (dble(numit) <= z) go to 190
C        ifault = 9
C        z = dble(numit)
C        call error(ifault, k, z)
C  190   if (digitx >= 0.5d0) go to 200
C        ifault = 10
C        z = digitx
C        call error(ifault, k, z)
  200 continue
c      if (ksum < l) go to 210
c      ifault = 7
c      call error(ifault, k, z)
c      return
c
c compute the unscaled covariance matrix of the coefficients.
c
  210 call covar(n, m1, n1, ipivot, r, c(k3), c(k9), nn)
c
c in certain problems, some diagonal terms of the unscaled covariance
c matrix are equal to zero or to small positive numbers.  because of
c rounding errors, computed values for these terms may be small
c negative numbers.  an error diagnostic is printed if any diagonal
c term is negative.
c
      do 220 j=1,n
        if (r(j,j) >= 0.0d0) go to 220
        ifault = 11
        z = r(j,j)
        call error(ifault, j, z)
  220 continue
      return
      end
      subroutine l2b(m, n, m1, l, a, b, w, tol, mm, nn, mmpnn,
     * n1, ipivot, x, res, qr, c, ifault)
c ** purpose **
c subroutine l2b computes least squares solutions to overdetermined
c and underdetermined systems of linear equations.  the method used is
c a modified gram-schmidt orthogonal decomposition with iterative
c refinement of the solution.  the solution may be subject to linear
c equality constraints.  output includes the least squares
c coefficients, residuals, unscaled covariance matrix, and information
c on the behavior of the iterative refinement procedure.
c matrix a is the given matrix of a system of m linear equations in n
c unknowns, and matrix w is a given diagonal matrix of weights with all
c diagonal elements nonnegative.  let h = (sqrt(w))*a.
c in the event that n1 (the computed rank of matrix h) is less than n
c (the number of unknown coefficients), a unique solution vector having
c n elements can be obtained by imposing the condition that the
c solution be of minimal euclidean norm.  such a solution is sought in
c the case of underdetermined or rank-deficient problems.
c
c ** input variables **
c m      total number of equations.
c n      number of unknown coefficients.
c m1     number of linear constraints (0 <= m1 <= m and m1 <= n).
c l      number of right-hand sides (vectors of observations).
c a      two-dimensional array of size (mm,n).  on entry, the array a
c        contains the given matrix of a system of m linear equations
c        in n unknowns, where the first m1 equations are to be
c        satisfied exactly.  a is left intact on exit.
c b      two-dimensional array of size (mm,l).  on entry, b contains
c        the l given right-hand sides (vectors of observations).  b is
c        left intact on exit.
c w      vector of size m.  on entry, w contains the diagonal elements
c        of a given diagonal matrix of weights, all nonnegative.
c        (the first m1 elements of w are set equal to 1.0 by the
c        program when m1 is greater than zero.)  on exit, the original
c        elements of w have been replaced by their square roots.
c tol    parameter used in determining the rank of matrix h.
c        note --
c        (1) if tol equals zero, the tolerance used in subroutine
c            decom2 will be based on machine precision.
c        (2) if tol is greater than zero, this value of tol will be
c            used in setting an absolute tolerance for comparison with
c            diagonal elements of the triangular matrix obtained in
c            subroutine decom2.  the value of tol can be based on
c            knowledge concerning the accuracy of the data.
c mm     dimensioning parameter specifying maximum number of rows in
c        the arrays a, b and res.  mm must satisfy mm >= m.
c nn     dimensioning parameter specifying maximum number of rows in
c        the array x.  nn must satisfy nn >= n.
c mmpnn  dimensioning parameter specifying maximum number of rows in
c        the array qr.  mmpnn must satisfy mmpnn >= m+n.
c
c ** output variables and internal variables **
c n1     computed rank of matrix h, where h = (sqrt(w))*a.
c ipivot vector of size n.  on exit, this array records the order
c        in which the columns of h were selected by the pivoting
c        scheme in the course of the orthogonal decomposition.
c        whenever n1 < n, the first n1 elements of ipivot indicate
c        which columns of h were found to be linearly independent.
c x      two-dimensional array of size (nn,l).  on exit, x contains
c        the solution vectors.
c res    two-dimensional array of size (mm,l).  on exit, res contains
c        the residual vectors.
c qr     two-dimensional array of size (mmpnn,n).  on exit, qr
c        contains the lower triangular portion of the symmetric
c        unscaled covariance matrix.  (this array is used internally
c        to store results from subroutine decom2 which are
c        destroyed in computing the covariance matrix.)
c c      vector having at least 6*(m+n)+2*l elements used (1) for
c        internal work space and (2) for returning information on the
c        behavior of the iterative refinement procedure.
c        (a) numit is the number of iterations carried out during the
c            iterative refinement in attempting to obtain a solution
c            for the k-th right-hand side.
c            on exit, c(k) = +numit if the solution converged, and
c                     c(k) = -numit if the solution failed to converge.
c        (b) digitx gives an estimate of the number of correct digits
c            in the initial solution of the coefficients for the k-th
c            right-hand side.  on exit, c(k+l) = digitx.
c ifault fault indicator which is zero if no errors were encountered
c        and positive if errors were detected or if evidence of severe
c        ill-conditioning was found.  diagnostic messages are printed
c        from subroutine error.  if ifault is set equal to 1, 2, 3, 4,
c        5, 6 or 7, execution is terminated.  execution continues when
c        ifault is set equal to 8, 9 or 10 provided that a solution
c        was obtained for at least one right-hand side.  the value of
c        ifault is used to indicate the following --
c        0 = no errors encountered.
c        1 = bad input parameter (m, n or l).
c        2 = bad input parameter (m1).
c        3 = bad dimension.  either m > mm, n > nn or m+n > mmpnn.
c        4 = at least one weight is negative.
c        5 = either matrix h or matrix of constraints equals zero.
c        6 = constraints are linearly dependent.
c        7 = all solutions failed to converge.
c        8 = solution failed to converge for at least one right-hand
c            side.
c        9 = large number of iterations required for convergence.
c       10 = estimated number of digits in initial solution of
c            coefficients is small.
c       11 = diagonal element of covariance matrix was computed to be
c            negative owing to rounding error.
c
c ** subroutines required **
c subroutine decom2
c        uses modified gram-schmidt algorithm with pivoting to
c        obtain an orthogonal decomposition of the input matrix.
c subroutine solve2
c        computes coefficients and residuals.  iterative refinement is
c        used to improve the accuracy of the initial solution.
c subroutine solve3
c        called only by subroutine solve2.
c subroutine covar
c        computes unscaled covariance matrix of the coefficients.
c subroutine error
c        prints error diagnostics when errors are detected or when
c        evidence of severe ill-conditioning is found.
c
c ** storage requirements **
c the storage required for the dimensioned arrays in subroutine l2b is
c   m*(2*n + 2*l + 7) + n*(n + l + 7) + 2*l
c locations.  all arrays required in subroutines called by l2b are
c declared herein and are transmitted only through parameter lists of
c call-sequences.
c
c ** precision of arithmetic calculations **
c single precision arithmetic is used for all calculations except the
c double precision accumulation of inner products.  (the variable sum
c is declared to be double precision in subroutines decom2, solve2,
c solve3 and covar.)  it is essential for the success of the iterative
c refinement procedure that inner products be accumulated in double
c precision.
c
c ** conversion of the program to double precision **
c *********************************************************************
c * on computers having short word length (as the ibm 360/370) it may *
c * be desirable to perform all calculations in double precision.  in *
c * this case, the iterative refinement presently included in solve2  *
c * should be omitted.                                                *
c * to convert the program to double precision, the following         *
c * approach is suggested.                                            *
c *                                                                   *
c * 1. variables presently declared to be real should be declared     *
c *    double precision.  those typed integer, double precision and   *
c *    logical should not be changed.                                 *
c * 2. the use of fail, numit and digitx should be omitted.           *
c * 3. description of variable c (at l2b 690-790) should read --      *
c *    c  vector having at least 6*(m+n) elements used only for       *
c *       internal work space.                                        *
c * 4. the value of eta (at l2b 1960) should be set so that it is the *
c *    smallest positive double precision number such that 1.0 + eta  *
c *    is greater than 1.0 in double precision arithmetic.            *
c *    for ibm computer type, eta = 16.**(-13)                        *
c *    for univac computer type, eta = 2.**(-59)                      *
c * 5. the following fortran functions should be changed --           *
c *      single precision name     double precision name              *
c *             dble(x)                x                              *
c *             float(n)               dble(float(n))                 *
c *             sqrt(x)                dsqrt(x)                       *
c *    dble(x) is used in subroutines decom2, solve2, solve3 and      *
c *       covar.                                                      *
c *    float(n) is used in subroutine decom2.                         *
c *    sqrt(x) is used in subroutine l2b.                             *
c * 6. it may be necessary or desirable to change certain formats in  *
c *    subroutine error, replacing g specifications by d.             *
c * 7. replace statement l2b 2500 by a statement reading              *
c *          k3 = 1                                                   *
c * 8. further details are given in subroutine solve2 in connection   *
c *    with the omission of iterative refinement.                     *
c * 9. in subroutine l2b, statements l2b 950-1000, 1820-1830, 2020,   *
c *    2350-2360, 2480-2490, 3070, 3280-3570 and 3590-3620 should be  *
c *    omitted.                                                       *
c *    statement numbers given here refer to those in the right-hand  *
c *    margin.  certain comments in subroutine l2b do not apply to    *
c *    the double precision version.                                  *
c *                                                                   *
c *********************************************************************
c
      integer ipivot(n)
      double precision a(mm,n), b(mm,l), c(1), eta, qr(mmpnn,n),
     * res(mm,l), tol, w(m), x(nn,l), z
      double precision digitx
      logical fail
      logical sing
c
c set value of eta, a machine-dependent parameter.
c eta, the relative machine precision, is the smallest positive real
c number such that 1.0 + eta is greater than 1.0 in floating-point
c arithmetic.
c
c for ibm computer type, eta = 16.**(-5)
c for univac computer type, eta = 2.**(-26)
c for cdc computer type, eta = 2.**(-47)
c for honeywell computer type, eta = 2.**(-27)
c
      eta = 2.d0**(-26)
c
c default value for tol is zero.
c
      if (tol < 0.0d0) tol = 0.0d0
      ifault = 0
      ksum = 0
c
c perform initial checking of input parameters, dimensions and
c weights for possible errors.
c
      if (m > 0 .and. n > 0 .and. l > 0) go to 10
      ifault = 1
      call error(ifault, k, z)
      return
   10 if (m1 <= m .and. m1 <= n .and. m1 >= 0) go to 20
      ifault = 2
      call error(ifault, k, z)
      return
   20 if (m <= mm .and. n <= nn .and. m+n <= mmpnn) go to 30
      ifault = 3
      k = 2
      call error(ifault, k, z)
      return
   30 do 40 i=1,m
        if (m1 > 0 .and. i <= m1) w(i) = 1.0d0
        if (w(i) >= 0.0d0) go to 40
        ifault = 4
        z = w(i)
        call error(ifault, i, z)
   40 continue
      if (ifault == 4) return
      do 50 i=1,m
        w(i) = dsqrt(w(i))
   50 continue
c
c set parameters which allocate vector c to contain certain final
c results and also to be used as work space.
c
c k1 is starting point for numit and fail, of length l.
c k2 is starting point for digitx, of length l.
c k3 is starting point for d, of length n.
c k4 is starting point for k-th column of b, of length m.
c k5 is starting point for k-th column of x, of length n.
c k6 is starting point for k-th column of res, of length m.
c k7 is starting point for work space of length m.
c k8 is starting point for work space of length m.
c k9 is starting point for work space of length n.
c k10 is starting point for work space of length n.
c k11 is starting point for work space of length m + n.
c k12 is starting point for work space of length m + n.
c
      k1 = 1
      k2 = k1 + l
      k3 = k2 + l
      k4 = k3 + n
      k5 = k4 + m
      k6 = k5 + n
      k7 = k6 + m
      k8 = k7 + m
      k9 = k8 + m
      k10 = k9 + n
      k11 = k10 + n
      k12 = k11 + m + n
      k = k12 + m + n - 1
c
c multiply each row of matrix a (m by n) by its appropriate weight and
c store the result in the first m rows of array qr.  set array c and
c the last n rows of array qr equal to zero.
c
      do 60 i=1,k
        c(i) = 0.0d0
   60 continue
      mp1 = m + 1
      mpn = m + n
      do 90 j=1,n
        do 70 i=1,m
          qr(i,j) = a(i,j)*w(i)
   70   continue
        do 80 i=mp1,mpn
          qr(i,j) = 0.0d0
   80   continue
   90 continue
c
c obtain an orthogonal decomposition of the matrix stored in the first
c m rows of array qr and compute its rank.
c
      call decom2(m, n, m1, eta, tol, qr, c(k3), n1, ipivot, sing,
     * mmpnn)
c
      if (.not.sing) go to 110
      if (n1 > 0) go to 100
      ifault = 5
      call error(ifault, k, z)
      return
  100 ifault = 6
      call error(ifault, k, z)
      return
c
c seek a solution (coefficients and residuals) for each of the l least
c squares problems whose right-hand sides are given in the array b.
c
  110 do 200 k=1,l
c k-th right-hand side.
        k0 = k4 - 1
        do 120 i=1,m
          k0 = k0 + 1
          c(k0) = b(i,k)
  120   continue
c
        call solve2(m, n, m1, a, c(k4), w, n1, ipivot, qr, c(k3),
     *   eta, fail, numit, digitx,
     *   c(k5), c(k6), c(k7), c(k8), c(k9), c(k10), c(k11), c(k12),
     *   mm, mmpnn)
c
        k0 = k5 - 1
        do 130 j=1,n
          k0 = k0 + 1
          x(j,k) = c(k0)
  130   continue
        if (m1 == 0) go to 150
        do 140 i=1,m1
          res(i,k) = 0.0d0
  140   continue
  150   m1p1 = m1 + 1
        if (m1p1 > m) go to 170
        k0 = k6 + m1 - 1
        do 160 i=m1p1,m
          k0 = k0 + 1
          res(i,k) = c(k0)
  160   continue
  170   continue
c
c for right-hand sides where convergence of a solution is reported,
c a check is made for evidence of severe ill-conditioning.  such
c evidence is furnished by large values of numit (number of iterations
c before convergence was obtained) and small values of digitx
c (estimate of the number of correct digits in the initial solution
c of the coefficients).  if numit exceeds -alog10(eta) a diagnostic
c message is printed to warn of ill-conditioning.  if digitx is less
c than 0.5 (half a decimal digit) a similar diagnostic message is
c printed.
c
        c(k) = dble(numit)
        if (fail) c(k) = -c(k)
        k0 = k2 + k - 1
        c(k0) = digitx
        if (.not.fail) go to 180
c ksum is a tally of solutions which failed to converge.
        ksum = ksum + 1
        ifault = 8
        call error(ifault, k, z)
        go to 200
  180   z = -dlog10(eta)
        if (dble(numit) <= z) go to 190
        ifault = 9
        z = dble(numit)
        call error(ifault, k, z)
  190   if (digitx >= 0.5d0) go to 200
        ifault = 10
        z = digitx
        call error(ifault, k, z)
  200 continue
      if (ksum < l) go to 210
      ifault = 7
      call error(ifault, k, z)
      return
  210 if (n1 < n) return
      do 230 i=1,n
        mpi = m + i
        do 220 j=1,n
          qr(i,j) = qr(mpi,j)
  220   continue
        qr(i,i) = 0.0d0
  230 continue
c
c compute the unscaled covariance matrix of the coefficients.
c
      call covar(n, m1, n1, ipivot, qr, c(k3), c(k9), mmpnn)
c
c in certain problems, some diagonal terms of the unscaled covariance
c matrix are equal to zero or to small positive numbers.  because of
c rounding errors, computed values for these terms may be small
c negative numbers.  an error diagnostic is printed if any diagonal
c term is negative.
c
      do 240 j=1,n
        if (qr(j,j) >= 0.0d0) go to 240
        ifault = 11
        z = qr(j,j)
        call error(ifault, j, z)
  240 continue
      return
      end
c     subroutine decom1(...)
c subroutine decom1 uses a modified gram-schmidt algorithm with
c pivoting to obtain an orthogonal decomposition of the input matrix
c given in q.
c the input parameter tol (equal either to zero or to a positive
c number) is used in determining the rank of matrix q.
c note --
c     (1) if tol equals zero, the tolerance used at statement dc1 1080
c         will be based on machine precision.
c         under this approach, the tolerance (tol2) is set equal to
c         (float(n)*eta)**2*d(m1+1) at statement dc1 1070.
c         if desired, the user can obtain a more conservative
c         tolerance by replacing n in this statement by a larger
c         quantity.
c     (2) if tol is greater than zero, tol2 (equal to the square of
c         tol) will be used at statement dc1 1080 as an absolute
c         tolerance for comparison with diagonal elements of the
c         triangular matrix obtained in the decomposition.  under this
c         approach, the value of tol can be based on knowledge
c         concerning the accuracy of the data.
c on exit, the arrays q, r, d and ipivot contain the results of the
c decomposition which are needed for obtaining an initial solution
c and for iterative refinement.
c on exit, n1 is the computed rank of the input matrix q.
c on exit, sing is set equal to .true. whenever
c     (1) n1 = 0 (i.e., input matrix q equals zero or matrix of
c         constraints equals zero), or
c     (2) n1 is less than m1 (i.e., the m1 by n matrix of linear
c         constraints is singular).
c otherwise, on exit from decom1, sing = .false.
c on exit, the vector ipivot records the order in which the columns
c of q were selected by the pivoting scheme in the course of the
c orthogonal decomposition.
      subroutine decom1(m, n, m1, eta, tol, q, r, d, n1, ipivot, sing,
     * mm, nn)
      integer ipivot(n)
      double precision c, d(1), dm, ds, eta, q(mm,n), r(nn,n),
     .  rsj, tol, tol2
      double precision sum
      logical fsum, sing
      n1 = 0
      sing = .true.
      fsum = .true.
      mv = 1
      mh = m1
      if (tol > 0.0d0) tol2 = tol**2
      do 10 j=1,n
        d(j) = 0.0d0
        ipivot(j) = j
   10 continue
c step number ns of the decomposition.
      do 210 ns=1,n
        nsm1 = ns - 1
        nsp1 = ns + 1
        if (ns == m1+1) go to 20
        go to 30
   20   if (m1 == m) go to 150
        mv = m1 + 1
        mh = m
        fsum = .true.
c pivot search.
   30   ds = 0.0d0
        np = ns
        do 80 j=ns,n
          ik = ipivot(j)
          if (fsum) go to 40
          go to 60
   40     sum = 0.0d0
          do 50 l=mv,mh
            sum = sum + dble(q(l,ik))*dble(q(l,ik))
   50     continue
          d(j) = sum
   60     if (ds < d(j)) go to 70
          go to 80
   70     ds = d(j)
          np = j
   80   continue
c end pivot search.
        ik = ipivot(np)
        if (fsum) dm = ds
        if (ds < eta*dm) go to 90
        fsum = .false.
        go to 100
   90   fsum = .true.
  100   if (fsum) go to 30
        if (np /= ns) go to 110
        go to 130
c column interchange.
  110   ipivot(np) = ipivot(ns)
        ipivot(ns) = ik
        d(np) = d(ns)
        if (ns == 1) go to 130
        do 120 l=1,nsm1
          c = r(l,np)
          r(l,np) = r(l,ns)
          r(l,ns) = c
  120   continue
c end column interchange.
c return here if n1 = 0.  either input matrix q equals zero or matrix
c of constraints equals zero.
  130   if (ns == 1 .and. ds == 0.0d0) return
        sum = 0.0d0
        do 140 l=mv,mh
          sum = sum + dble(q(l,ik))*dble(q(l,ik))
  140   continue
        d(ns) = sum
        ds = d(ns)
        if (tol == 0.0d0) tol2 = (dble(n)*eta)**2*d(m1+1)
        if (ns > m1 .and. ds <= tol2) go to 150
        go to 160
  150   sing = .false.
c return here if n1 < n, n1 > 0 and n1 >= m1.
        return
c return here if matrix of constraints is found to be singular.
  160   if (ds == 0.0d0) return
        if (nsp1 > n) go to 200
c begin orthogonalization.
        do 190 j=nsp1,n
          np = ipivot(j)
          sum = 0.0d0
          do 170 l=mv,mh
            sum = sum + dble(q(l,np))*dble(q(l,ik))
  170     continue
          rsj = sum
          rsj = rsj/ds
          r(ns,j) = rsj
          do 180 l=mv,m
            q(l,np) = q(l,np) - rsj*q(l,ik)
  180     continue
          d(j) = d(j) - ds*rsj*rsj
  190   continue
c end orthogonalization.
  200   n1 = n1 + 1
  210 continue
c end step number ns.
      sing = .false.
c normal return.  n1 = n.
      return
      end
c     subroutine decom2(...)
c subroutine decom2 uses a modified gram-schmidt algorithm with
c pivoting to obtain an orthogonal decomposition of the input matrix
c given in qr.
c the input parameter tol (equal either to zero or to a positive
c number) is used in determining the rank of matrix qr.
c note --
c     (1) if tol equals zero, the tolerance used at statement dc2 1180
c         will be based on machine precision.
c         under this approach, the tolerance (tol2) is set equal to
c         (float(n)*eta)**2*d(m1+1) at statement dc2 1170.
c         if desired, the user can obtain a more conservative
c         tolerance by replacing n in this statement by a larger
c         quantity.
c     (2) if tol is greater than zero, tol2 (equal to the square of
c         tol) will be used at statement dc2 1180 as an absolute
c         tolerance for comparison with diagonal elements of the
c         triangular matrix obtained in the decomposition.  under this
c         approach, the value of tol can be based on knowledge
c         concerning the accuracy of the data.
c on exit, the arrays qr, d and ipivot contain the results of the
c decomposition which are needed for obtaining an initial solution
c and for iterative refinement.
c on exit, n1 is the computed rank of the input matrix qr.
c on exit, sing is set equal to .true. whenever
c     (1) n1 = 0 (i.e., input matrix qr equals zero or matrix of
c         constraints equals zero), or
c     (2) n1 is less than m1 (i.e., the m1 by n matrix of linear
c         constraints is singular).
c otherwise, on exit from decom2, sing = .false.
c on exit, the vector ipivot records the order in which the columns
c of qr were selected by the pivoting scheme in the course of the
c orthogonal decomposition.
      subroutine decom2(m, n, m1, eta, tol, qr, d, n1, ipivot, sing,
     * mmpnn)
      integer ipivot(n)
      double precision c, d(1), dm, ds, eta, qr(mmpnn,n), rsj, tol, tol2
      double precision sum
      logical finis, fsum, sing
      n1 = 0
      sing = .true.
      fsum = .true.
      mv = 1
      mh = m1
      ms = m
      mp1 = m + 1
      finis = .false.
      if (tol > 0.0d0) tol2 = tol**2
      do 10 j=1,n
        d(j) = 0.0d0
        ipivot(j) = j
   10 continue
c step number ns of the decomposition.
      do 350 ns=1,n
        k = m + ns
        if (ns == m1+1) go to 20
        go to 30
   20   if (m1 == m) go to 200
        mv = m1 + 1
        mh = m
        fsum = .true.
   30   if (.not.finis) go to 40
        go to 150
c pivot search.
   40   ds = 0.0d0
        np = ns
        do 90 j=ns,n
          if (fsum) go to 50
          go to 70
   50     sum = 0.0d0
          do 60 l=mv,mh
            sum = sum + dble(qr(l,j))*dble(qr(l,j))
   60     continue
          d(j) = sum
   70     if (ds < d(j)) go to 80
          go to 90
   80     ds = d(j)
          np = j
   90   continue
        if (fsum) dm = ds
        if (ds < eta*dm) go to 100
        fsum = .false.
        go to 110
  100   fsum = .true.
  110   if (fsum) go to 40
        if (np /= ns) go to 120
        go to 140
c column interchange.
  120   ik = ipivot(np)
        ipivot(np) = ipivot(ns)
        ipivot(ns) = ik
        d(np) = d(ns)
        km1 = k - 1
        do 130 l=1,km1
          c = qr(l,np)
          qr(l,np) = qr(l,ns)
          qr(l,ns) = c
  130   continue
c end column interchange.
c end pivot search.
c return here if n1 = 0.  either input matrix qr equals zero or
c matrix of constraints equals zero.
  140   if (ns == 1 .and. ds == 0.0d0) return
        go to 160
  150   ms = k - 1
        mh = k - 1
  160   if (finis) go to 170
        c = 0.0d0
        go to 180
  170   c = 1.0d0
  180   sum = dble(c)
        do 190 l=mv,mh
          sum = sum + dble(qr(l,ns))*dble(qr(l,ns))
  190   continue
        d(ns) = sum
        ds = d(ns)
        if (tol == 0.0d0) tol2 = (dble(n)*eta)**2*d(m1+1)
        if (.not.finis .and. ns > m1 .and. ds <= tol2) go to 200
        go to 290
  200   finis = .true.
        mv = m + 1
        do 280 np=ns,n
          if (1 > m1) go to 250
          do 210 l=1,m1
            qr(l,np) = 0.0d0
  210     continue
          do 240 j=1,m1
            sum = 0.0d0
            do 220 l=1,m
              sum = sum + dble(qr(l,j))*dble(qr(l,np))
  220       continue
            c = sum
            c = c/d(j)
            do 230 l=1,m1
              qr(l,np) = qr(l,np) - c*qr(l,j)
  230       continue
  240     continue
  250     mpn1 = m + n1
          do 270 jj=mp1,mpn1
            j = (m + 1) + (m + n1) - jj
            sum = 0.0d0
            do 260 l=j,mpn1
              lmm = l - m
              sum = sum + dble(qr(j,lmm))*dble(qr(l,np))
  260       continue
            qr(j,np) = -sum
  270     continue
  280   continue
        go to 150
c return here if matrix of constraints is found to be singular.
  290   if (ds == 0.0d0) return
        qr(k,ns) = -1.0d0
        nsp1 = ns + 1
        if (nsp1 > n) go to 340
c begin orthogonalization.
        do 330 j=nsp1,n
          sum = 0.0d0
          do 300 l=mv,mh
            sum = sum + dble(qr(l,j))*dble(qr(l,ns))
  300     continue
          rsj = sum
          rsj = rsj/ds
          qr(k,j) = rsj
          do 310 l=1,ms
            qr(l,j) = qr(l,j) - rsj*qr(l,ns)
  310     continue
          if (.not.finis) go to 320
          go to 330
  320     d(j) = d(j) - ds*rsj*rsj
  330   continue
c end orthogonalization.
  340   if (.not.finis) n1 = n1 + 1
  350 continue
c end step number ns.
      sing = .false.
c normal return.
      return
      end
c     subroutine solve1(...)
c subroutine solve1 uses the orthogonal decomposition stored in q, r,
c d and ipivot to compute the solution (coefficients and residuals)
c to the least squares problem whose right-hand side is given in b.
c in the event that n1 (the computed rank of matrix h) is less than n
c (the number of unknown coefficients), the original matrix h (m by n)
c is replaced by a smaller matrix (m by n1) whose columns are linearly
c independent, and a solution is sought for the smaller system of
c equations.  thus n - n1 columns of the original matrix h are deleted,
c and coefficients corresponding to these n - n1 deleted columns will
c be set equal to zero.
c in normal exits, the solution is contained in the vector x
c (coefficients) and the vector res (residuals).
c iterative refinement is used to improve the accuracy of the initial
c solution.
c on exit, fail is set equal to .true. if the solution fails to
c improve sufficiently.  otherwise, fail = .false.  information on the
c behavior of the iterative refinement procedure is given by numit and
c digitx.  numit is the number of iterations carried out in attempting
c to obtain a solution.  digitx is an estimate of the number of
c correct digits in the initial solution of the coefficients.
c
c ********* conversion of this subroutine to double precision *********
c * if the program is converted so that all calculations are done in  *
c * double precision arithmetic, the iterative refinement presently   *
c * included in solve1 should be omitted, since the success of this   *
c * procedure depends on computing inner products in greater          *
c * precision than other calculations.                                *
c * see comments in subroutine l2a regarding conversion to double     *
c * precision.  in addition, the following comments indicate how to   *
c * omit the iterative refinement from this subroutine.  statement    *
c * numbers given here refer to those in the right-hand margin.       *
c *                                                                   *
c * 1. in statement sv1 480 change real to double precision.          *
c * 2. replace statement sv1 810 by a statement reading               *
c *       30 do 50 i=1,m                                              *
c * 3. after statement sv1 1050 insert a statement reading            *
c *          return                                                   *
c * 4. omit statements sv1 140-210, 450, 500-510, 530-560, 610,       *
c *    680-780, 1600-1780 and 1800-1860.                              *
c *                                                                   *
c *********************************************************************
c
      subroutine solve1(m, n, m1, a, b, w, n1, ipivot, q, r, d,
C     * eta, fail, numit, digitx,
     * x, res, f, wres, g, y, mm, nn)
      integer ipivot(n)
      double precision a(mm,n), b(1), c, d(1), f(1), g(1), q(mm,n),
     * r(nn,n), res(1), w(m), wres(1), x(1), y(1)
C      double precision digitx, dxnorm, eta, eta2,
C     .  rdr1, rdr2, rdx1, rdx2, rnr, rnx, xnorm
      double precision sum
C      logical fail
C      numit = 0
C      kz = 0
C      eta2 = eta*eta
      do 10 i=1,m
        f(i) = b(i)*w(i)
        wres(i) = 0.0d0
        res(i) = 0.0d0
c        if (w(i) == 0.0d0) kz = kz + 1
   10 continue
      do 20 j=1,n
        x(j) = 0.0d0
        g(j) = 0.0d0
   20 continue
      k = 0
c      rdx2 = 0.0d0
c      rdr2 = 0.0d0
cc begin k-th iteration step.
c   30 if (k < 2) go to 40
c      if (((64.d0*rdx2 < rdx1) .and. (rdx2 > eta2*rnx)) .or.
c     * ((64.d0*rdr2 < rdr1) .and. (rdr2 > eta2*rnr))) go to 40
c      go to 300
c   40 rdx1 = rdx2
c      rdr1 = rdr2
c      rdx2 = 0.0d0
c      rdr2 = 0.0d0
      if (k == 0) go to 100
c new residuals.
   30 do 50 i=1,m
        wres(i) = wres(i) + f(i)*w(i)
        if (w(i) == 0.0d0) go to 50
        res(i) = res(i) + f(i)/w(i)
   50 continue
      do 70 ns=1,n1
        j = ipivot(ns)
        x(j) = x(j) + g(ns)
        sum = 0.0d0
        do 60 l=1,m
          sum = sum + dble(a(l,j))*dble(wres(l))
   60   continue
        g(ns) = -sum
   70 continue
      do 90 i=1,m
        sum = 0.0d0
        if (i > m1) sum = dble(res(i))
        do 80 l=1,n
          sum = sum + dble(a(i,l))*dble(x(l))
   80   continue
        sum = sum - dble(b(i))
        f(i) = -sum
        f(i) = f(i)*w(i)
        if (w(i) == 0.0d0) res(i) = dble(res(i)) - sum
   90 continue
      return
c end new residuals.
  100 mv = 1
      mh = m1
      do 160 ns=1,n1
        j = ipivot(ns)
        if (ns /= m1+1) go to 110
        mv = m1 + 1
        mh = m
  110   nsm1 = ns - 1
        sum = -dble(g(ns))
        if (1 > nsm1) go to 130
        do 120 l=1,nsm1
          sum = sum + dble(r(l,ns))*dble(y(l))
  120   continue
  130   y(ns) = -sum
        c = 0.0d0
        if (ns > m1) c = -y(ns)
        sum = dble(c)
        do 140 l=mv,mh
          sum = sum + dble(q(l,j))*dble(f(l))
  140   continue
        c = sum
        c = c/d(ns)
        g(ns) = c
        do 150 i=mv,m
          f(i) = f(i) - c*q(i,j)
  150   continue
  160 continue
      if (1 > m1) go to 210
      do 170 i=1,m1
        f(i) = 0.0d0
  170 continue
      do 200 ns=1,m1
        j = ipivot(ns)
        sum = -dble(y(ns))
        do 180 l=1,m
          sum = sum + dble(q(l,j))*dble(f(l))
  180   continue
        c = sum
        c = c/d(ns)
        do 190 i=1,m1
          f(i) = f(i) - c*q(i,j)
  190   continue
  200 continue
  210 do 240 i=1,n1
        ns = n1 + 1 - i
        nsp1 = ns + 1
        sum = -dble(g(ns))
        if (nsp1 > n1) go to 230
        do 220 l=nsp1,n1
          sum = sum + dble(r(ns,l))*dble(g(l))
  220   continue
  230   g(ns) = -sum
  240 continue
C      do 250 ns=1,n1
C        rdx2 = rdx2 + g(ns)*g(ns)
C  250 continue
C      do 260 i=1,m
C        rdr2 = rdr2 + f(i)*f(i)
C  260 continue
C      if (k /= 0) go to 270
C      rnx = rdx2
C      rnr = rdr2
C  270 if (k /= 1) go to 290
C      xnorm = dsqrt(rnx)
C      dxnorm = dsqrt(rdx2)
C      if (xnorm /= 0.0d0) go to 280
C      digitx = -dlog10(eta)
C      go to 290
C  280 digitx = -dlog10(dmax1(dxnorm/xnorm,eta))
Cc end k-th iteration step.
C  290 numit = k
C      k = k + 1
      go to 30
C  300 if ((m1+kz == m) .and. (rdx2 > 4.d0*eta2*rnx)) go to 310
C      if ((rdr2 > 4.d0*eta2*rnr) .and.
C     * (rdx2 > 4.d0*eta2*rnx)) go to 310
C      fail = .false.
C      return
C  310 fail = .true.
C      return
      end
c     subroutine solve2(...)
c subroutine solve2 uses the orthogonal decomposition stored in qr, d
c and ipivot to compute the solution (coefficients and residuals)
c to the least squares problem whose right-hand side is given in b.
c in the event that n1 (the computed rank of matrix h) is less than n
c (the number of unknown coefficients), a unique solution vector having
c n elements can be obtained by imposing the condition that the
c solution be of minimal euclidean norm.  such a solution is sought in
c the case of underdetermined or rank-deficient problems.
c in normal exits, the solution is contained in the vector x
c (coefficients) and the vector res (residuals).
c iterative refinement is used to improve the accuracy of the initial
c solution.
c on exit, fail is set equal to .true. if the solution fails to
c improve sufficiently.  otherwise, fail = .false.  information on the
c behavior of the iterative refinement procedure is given by numit and
c digitx.  numit is the number of iterations carried out in attempting
c to obtain a solution.  digitx is an estimate of the number of
c correct digits in the initial solution of the coefficients.
c this subroutine calls subroutine solve3.
c
c ********* conversion of this subroutine to double precision *********
c * if the program is converted so that all calculations are done in  *
c * double precision arithmetic, the iterative refinement presently   *
c * included in solve2 should be omitted, since the success of this   *
c * procedure depends on computing inner products in greater          *
c * precision than other calculations.                                *
c * see comments in subroutine l2b regarding conversion to double     *
c * precision.  in addition, the following comments indicate how to   *
c * omit the iterative refinement from this subroutine.  statement    *
c * numbers given here refer to those in the right-hand margin.       *
c *                                                                   *
c * 1. in statement sv2 470 change real to double precision.          *
c * 2. replace statement sv2 880 by a statement reading               *
c *       30 do 50 i=1,m                                              *
c * 3. replace statements sv2 1310-1400 by a statement reading        *
c *          return                                                   *
c * 4. omit statements sv2 120-190, 440, 490-500, 520-550, 650,       *
c *    750-850, 1650-1830 and 1850-1910.                              *
c *                                                                   *
c *********************************************************************
c
      subroutine solve2(m, n, m1, a, b, w, n1, ipivot, qr, d,
     * eta, fail, numit, digitx,
     * x, res, wres, y1, y2, y, f, g, mm, mmpnn)
      integer ipivot(n)
      double precision a(mm,n), b(1), c, d(1), f(1), g(1),
     * qr(mmpnn,n), res(1), w(m), wres(1), x(1), y(1), y1(1), y2(1)
      double precision digitx, dxnorm, eta, eta2, rdr1, rdr2,
     .  rdx1, rdx2, rnr, rnx, xnorm
      double precision sum
      logical fail
      numit = 0
      kz = 0
      eta2 = eta*eta
      mp1 = m + 1
      mpn = m + n
      n1p1 = n1 + 1
      do 10 i=1,m
        f(i) = b(i)*w(i)
        g(i) = 0.0d0
        wres(i) = 0.0d0
        res(i) = 0.0d0
        y1(i) = 0.0d0
        if (w(i) == 0.0d0) kz = kz + 1
   10 continue
      do 20 ns=1,n
        j = m + ns
        f(j) = 0.0d0
        g(j) = 0.0d0
        x(ns) = 0.0d0
        y2(ns) = 0.0d0
   20 continue
      k = 0
      rdx2 = 0.0d0
      rdr2 = 0.0d0
c begin k-th iteration step.
   30 if (k < 2) go to 40
      if (((64.d0*rdx2 < rdx1) .and. (rdx2 > eta2*rnx)) .or.
     * ((64.d0*rdr2 < rdr1) .and. (rdr2 > eta2*rnr))) go to 40
      go to 270
   40 rdx1 = rdx2
      rdr1 = rdr2
      rdx2 = 0.0d0
      rdr2 = 0.0d0
      if (k == 0) go to 160
c new residuals.
      do 50 i=1,m
        wres(i) = wres(i) + f(i)*w(i)
        if (w(i) == 0.0d0) go to 50
        res(i) = res(i) + f(i)/w(i)
        y1(i) = y1(i) + g(i)
   50 continue
      do 100 ns=1,n
        j = m + ns
        np = ipivot(ns)
        x(np) = x(np) + f(j)
        y2(np) = y2(np) + g(j)
        sum = -dble(x(np))
        do 60 l=1,m
          sum = sum + dble(a(l,np))*dble(y1(l))
   60   continue
        g(j) = -sum
        if (ns > n1) go to 70
        go to 80
   70   f(j) = 0.0d0
        go to 100
   80   sum = 0.0d0
        do 90 l=1,m
          sum = sum + dble(a(l,np))*dble(wres(l))
   90   continue
        f(j) = -sum
  100 continue
      do 130 i=1,m
        sum = 0.0d0
        if (i > m1) sum = dble(res(i))
        do 110 l=1,n
          sum = sum + dble(a(i,l))*dble(x(l))
  110   continue
        sum = sum - dble(b(i))
        f(i) = -sum
        f(i) = f(i)*w(i)
        if (w(i) == 0.0d0) res(i) = dble(res(i)) - sum
        sum = 0.0d0
        if (i > m1) sum = dble(y1(i))
        do 120 l=1,n
          sum = sum + dble(a(i,l))*dble(y2(l))
  120   continue
        g(i) = -sum
  130 continue
      if (n1p1 > n) go to 160
      do 150 i=n1p1,n
        ns = n + n1p1 - i
        j = m + ns
        sum = 0.0d0
        do 140 l=1,j
          sum = sum + dble(qr(l,ns))*dble(g(l))
  140   continue
        g(j) = sum
  150 continue
c end new residuals.
c
  160 call solve3(f, m1, m, n1, qr, d, y, mmpnn)
c
      if (n1p1 > n) go to 200
      do 190 ns=n1p1,n
        j = m + ns
        sum = dble(g(j))
        do 170 l=mp1,j
          sum = sum + dble(qr(l,ns))*dble(f(l))
  170   continue
        c = sum
        c = c/d(ns)
        do 180 i=1,j
          f(i) = f(i) - c*qr(i,ns)
  180   continue
  190 continue
  200 do 210 j=mp1,mpn
        g(j) = 0.0d0
        if (j <= m+n1) g(j) = g(j) + f(j)
  210 continue
c
      call solve3(g, m1, m, n1, qr, d, y, mmpnn)
c
      do 220 i=1,m
        rdr2 = rdr2 + f(i)*f(i)
  220 continue
      do 230 i=mp1,mpn
        rdx2 = rdx2 + f(i)*f(i)
  230 continue
      if (k /= 0) go to 240
      rnr = rdr2
      rnx = rdx2
  240 if (k /= 1) go to 260
      xnorm = dsqrt(rnx)
      dxnorm = dsqrt(rdx2)
      if (xnorm /= 0.0d0) go to 250
      digitx = -dlog10(eta)
      go to 260
  250 digitx = -dlog10(dmax1(dxnorm/xnorm,eta))
c end k-th iteration step.
  260 numit = k
      k = k + 1
      go to 30
  270 if ((m1+kz == m) .and. (rdx2 > 4.d0*eta2*rnx)) go to 280
      if ((rdr2 > 4.d0*eta2*rnr) .and.
     * (rdx2 > 4.d0*eta2*rnx)) go to 280
      fail = .false.
      return
  280 fail = .true.
      return
      end
      subroutine solve3(f, m1, m, n1, qr, d, y, mmpnn)
c subroutine solve3 is called only by subroutine solve2.
c this subroutine calculates new values of f.
      double precision c, d(1), f(1), qr(mmpnn,n1), y(1)
      double precision sum
      mv = 1
      mh = m1
      do 100 ns=1,n1
        j = m + ns
        if (ns == m1+1) go to 10
        go to 20
   10   mv = m1 + 1
        mh = m
   20   nsm1 = ns - 1
        sum = -dble(f(j))
        if (ns == 1) go to 40
        do 30 l=1,nsm1
          mpl = m + l
          sum = sum + dble(qr(mpl,ns))*dble(y(l))
   30   continue
   40   y(ns) = -sum
        if (ns > m1) go to 50
        go to 60
   50   c = -y(ns)
        go to 70
   60   c = 0.0d0
   70   sum = dble(c)
        do 80 l=mv,mh
          sum = sum + dble(qr(l,ns))*dble(f(l))
   80   continue
        c = sum
        c = c/d(ns)
        f(j) = c
        do 90 l=mv,m
          f(l) = f(l) - c*qr(l,ns)
   90   continue
  100 continue
      if (1 > m1) go to 150
      do 110 l=1,m1
        f(l) = 0.0d0
  110 continue
      do 140 ns=1,m1
        sum = -dble(y(ns))
        do 120 l=1,m
          sum = sum + dble(qr(l,ns))*dble(f(l))
  120   continue
        c = sum
        c = c/d(ns)
        do 130 l=1,m1
          f(l) = f(l) - c*qr(l,ns)
  130   continue
  140 continue
  150 do 170 ns=1,n1
        j = m + n1 + 1 - ns
        mpn1 = m + n1
        sum = 0.0d0
        do 160 l=j,mpn1
          lmm = l - m
          sum = sum + dble(qr(j,lmm))*dble(f(l))
  160   continue
        f(j) = -sum
  170 continue
      return
      end
      subroutine covar(n, m1, n1, ipivot, c, d, z, nn)
c subroutine covar uses results from the orthogonal decomposition
c stored in c, d and ipivot to compute the unscaled covariance matrix
c of the least squares coefficients.
c on entry, the first n rows and the first n columns of c contain the
c upper triangular matrix obtained from the decomposition.  this input
c matrix is destroyed in subsequent calculations.
c on exit, the lower triangular portion of the symmetric unscaled
c covariance matrix is contained in
c     c(1,1)
c     c(2,1) c(2,2)
c     . . .
c     c(n,1) c(n,2) ... c(n,n)
c if n1 is less than n, one or more columns of the matrix
c h = (sqrt(w))*a were rejected as being linearly dependent.  whenever
c the k-th column of h was so rejected, c(i,j) is set equal to zero,
c for i = k or j = k, i >= j.
      integer ipivot(n)
      double precision c(nn,n), d(1), z(1)
      double precision sum
      l = n1
      if (l > m1) c(l,l) = 1.0d0/d(l)
      if (l == 1) go to 60
   10 j = l - 1
      if (j > m1) c(j,j) = 1.0d0/d(j)
      do 20 k=l,n1
        z(k) = c(j,k)
   20 continue
      i = n1
      do 40 ka=j,n1
        sum = 0.0d0
        if (i == j) sum = dble(c(i,j))
        do 30 k=l,n1
          sum = sum - dble(z(k))*dble(c(k,i))
   30   continue
        c(i,j) = sum
        i = i - 1
   40 continue
      do 50 k=l,n1
        c(j,k) = c(k,j)
   50 continue
      l = l - 1
      if (l > 1) go to 10
   60 if (n1 == n) go to 90
      n1p1 = n1 + 1
      do 80 i=1,n
        do 70 j=n1p1,n
          c(i,j) = 0.0d0
   70   continue
   80 continue
c permute the columns and rows of matrix c to account for pivoting.
   90 do 120 i=1,n
        do 100 j=1,n
          k = ipivot(j)
          z(k) = c(i,j)
  100   continue
        do 110 j=1,n
          c(i,j) = z(j)
  110   continue
  120 continue
      do 150 i=1,n
        do 130 j=1,n
          k = ipivot(j)
          z(k) = c(j,i)
  130   continue
        do 140 j=i,n
          c(j,i) = z(j)
  140   continue
  150 continue
      return
      end
      subroutine error(ifault, k, z)
c subroutine error prints error diagnostics in the case of error
c failure.
c also printed are some informative diagnostics related to the
c iterative refinement of ill-conditioned systems of equations and to
c rounding error problems in computing the covariance matrix.
c in the data statement below, nw is the printer device number.
      double precision z
      data nw /6/
      go to (10,20,30,40,50,60,70,80,90,100,110), ifault
   10 write (nw,99999)
      return
   20 write (nw,99998)
      return
   30 write (nw,99997)
      if (k == 2) write (nw,99996)
      return
   40 write (nw,99995) k,z
      return
   50 write (nw,99994)
      return
   60 write (nw,99993)
      return
   70 write (nw,99992)
      return
   80 write (nw,99991) k
      return
   90 write (nw,99990) k,z
      return
  100 write (nw,99989) k,z
      return
  110 write (nw,99988) k,z
      return
c format statements.
99999 format (50h0***  parameter error.  m, n and l must be greater,
     * 11h than zero.)
99998 format (50h0***  parameter error.  m1 cannot exceed m or n, b,
     * 26hut m1 must be nonnegative.)
99997 format (50h0***  dimension error.  one or more of the followi,
     * 32hng error conditions was found --/7x,12hm exceeds mm/7x,
     * 12hn exceeds nn)
99996 format (5x,17hm+n exceeds mmpnn)
99995 format (43h0***  weights must be nonnegative.  for i =,i3,2x,
     * 9hweight = ,g15.8)
99994 format (50h0***  either matrix h equals zero or matrix of con,
     * 51hstraints equals zero.  no solution can be computed.)
99993 format (50h0***  since the constraints are linearly dependent,
     * 29h no solution can be computed.)
99992 format (39h0***  all solutions failed to converge.)
99991 format (22h0***  for b-vector no.,i3,21h solution failed to c,
     * 8honverge.)
99990 format (22h0***  for b-vector no.,i3,21h the number of iterat,
     * 45hions required for convergence of solution was,f4.0/4h ***,
     * 57h  this number is large, indicating the problem is ill-con,
     * 40hditioned and solution may be inaccurate.)
99989 format (22h0***  for b-vector no.,i3,21h estimated number of ,
     * 54hcorrect digits in initial solution of coefficients is ,
     * g15.8/51h ***  since this is small, the final solution may b,
     * 13he inaccurate.)
99988 format (26h0***  diagonal element no.,i3,17h of the unscaled ,
     * 57hcovariance matrix was computed to be negative owing to ro,
     * 13hunding error./29h ***  the computed value was ,g15.8)
      end
