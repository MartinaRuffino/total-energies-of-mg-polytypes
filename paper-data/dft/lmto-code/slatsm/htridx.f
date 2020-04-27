C#define UNROLL
      subroutine htridx(nm,n,ar,ai,d,e,e2,tau)
C- a RISC adaptation of htridi (also parallel for SGI)
C  NB: SGI-PARALLEL requires e2 be dimensioned e2(n,2) for extra wk spc
      implicit none
      integer i,j,k,l,n,nm,jj
      double precision ar(nm,n),ai(nm,n),d(n),e(n),e2(n,2),tau(n,2)
      double precision f,g,h,fi,gi,hh,si,scale,pythag,dr(4),di(4)
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure tred1, num. math. 11, 181-195(1968)
c     by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a complex hermitian matrix
c     to a real symmetric tridiagonal matrix using
c     unitary similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex hermitian input matrix.
c          only the lower triangle of the matrix need be supplied.
c
c     on output
c
c        ar and ai contain information about the unitary trans-
c          formations used in the reduction in their full lower
c          triangles.  their strict upper triangles and the
c          diagonal of ar are unaltered.
c
c        d contains the diagonal elements of the the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2 may coincide with e if the squares are not needed.
c
c        tau contains further information about the transformations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version adapted from august 1983 htridi.  Differences
C     with htridi are that indices of tau and a are permuted (uses
C     upper triangle of a)
c     ------------------------------------------------------------------
C#ifdefC SGI-PARALLEL
C      integer mp_numthreads,nproc
C      nproc = mp_numthreads()
C#endif

      call tcn('htridx')

      tau(n,1) = 1.0d0
      tau(n,2) = 0.0d0
c
      do  100  i = 1, n
  100 d(i) = ar(i,i)
      do  300  i = n, 1, -1
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l < 1) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do  120  k = 1, l
  120    scale = scale + dabs(ar(k,i)) + dabs(ai(k,i))
c
         if (scale /= 0.0d0) go to 140
         tau(l,1) = 1.0d0
         tau(l,2) = 0.0d0
  130    e(i) = 0.0d0
         e2(i,1) = 0.0d0
         go to 290
c
  140    do  150  k = 1, l
            ar(k,i) = ar(k,i) / scale
            ai(k,i) = ai(k,i) / scale
            h = h + ar(k,i) * ar(k,i) + ai(k,i) * ai(k,i)
  150    continue
c
         e2(i,1) = scale * scale * h
         g = dsqrt(h)
         e(i) = scale * g
         f = pythag(ar(l,i),ai(l,i))
c     .......... form next diagonal element of matrix t ..........
         if (f == 0.0d0) go to 160
         tau(l,1) = (- ai(l,i) * tau(i,2) - ar(l,i) * tau(i,1)) / f
         si = (ar(l,i) * tau(i,2) - ai(l,i) * tau(i,1)) / f
         h = h + f * g
         g = 1.0d0 + g / f
         ar(l,i) = g * ar(l,i)
         ai(l,i) = g * ai(l,i)
         if (l == 1) go to 270
         go to 170
  160    tau(l,1) = -tau(i,1)
         si = tau(i,2)
         ar(l,i) = g
  170    continue
         f = 0.0d0
C#ifdef UNROLL | SGI-PARALLEL
C#ifdefC SGI-PARALLEL
CC$DOACROSS local(dr,di,j,k,jj)
C#endif
         do  1240  j = 1, l-3, 4
c .......... form element of a*u ..........
           dr(1) = 0
           di(1) = 0
           dr(2) = 0
           di(2) = 0
           dr(3) = 0
           di(3) = 0
           dr(4) = 0
           di(4) = 0
           do  1180  k = 1, j
             dr(1) = dr(1) + ai(k,j+0) * ai(k,i) + ar(k,j+0) * ar(k,i)
             di(1) = di(1) - ai(k,j+0) * ar(k,i) + ar(k,j+0) * ai(k,i)
             dr(2) = dr(2) + ai(k,j+1) * ai(k,i) + ar(k,j+1) * ar(k,i)
             di(2) = di(2) - ai(k,j+1) * ar(k,i) + ar(k,j+1) * ai(k,i)
             dr(3) = dr(3) + ai(k,j+2) * ai(k,i) + ar(k,j+2) * ar(k,i)
             di(3) = di(3) - ai(k,j+2) * ar(k,i) + ar(k,j+2) * ai(k,i)
             dr(4) = dr(4) + ai(k,j+3) * ai(k,i) + ar(k,j+3) * ar(k,i)
             di(4) = di(4) - ai(k,j+3) * ar(k,i) + ar(k,j+3) * ai(k,i)
C             e(j) = e(j) + ar(k,j) * ar(k,i) + ai(k,j) * ai(k,i)
C             tau(j,2) = tau(j,2) + ar(k,j) * ai(k,i) - ai(k,j) * ar(k,i)
 1180      continue
           do  1200  k = j+4, l
             dr(1) = dr(1)  + ar(j+0,k) * ar(k,i) - ai(j+0,k) * ai(k,i)
             di(1) = di(1)  + ar(j+0,k) * ai(k,i) + ai(j+0,k) * ar(k,i)
             dr(2) = dr(2)  + ar(j+1,k) * ar(k,i) - ai(j+1,k) * ai(k,i)
             di(2) = di(2)  + ar(j+1,k) * ai(k,i) + ai(j+1,k) * ar(k,i)
             dr(3) = dr(3)  + ar(j+2,k) * ar(k,i) - ai(j+2,k) * ai(k,i)
             di(3) = di(3)  + ar(j+2,k) * ai(k,i) + ai(j+2,k) * ar(k,i)
             dr(4) = dr(4)  + ar(j+3,k) * ar(k,i) - ai(j+3,k) * ai(k,i)
             di(4) = di(4)  + ar(j+3,k) * ai(k,i) + ai(j+3,k) * ar(k,i)
 1200      continue
           k = j+1
           dr(1) = dr(1)  + ar(j+0,k) * ar(k,i) - ai(j+0,k) * ai(k,i)
           di(1) = di(1)  + ar(j+0,k) * ai(k,i) + ai(j+0,k) * ar(k,i)
           dr(2) = dr(2) + ai(k,j+1) * ai(k,i) + ar(k,j+1) * ar(k,i)
           di(2) = di(2) - ai(k,j+1) * ar(k,i) + ar(k,j+1) * ai(k,i)
           dr(3) = dr(3) + ai(k,j+2) * ai(k,i) + ar(k,j+2) * ar(k,i)
           di(3) = di(3) - ai(k,j+2) * ar(k,i) + ar(k,j+2) * ai(k,i)
           dr(4) = dr(4) + ai(k,j+3) * ai(k,i) + ar(k,j+3) * ar(k,i)
           di(4) = di(4) - ai(k,j+3) * ar(k,i) + ar(k,j+3) * ai(k,i)
           k = j+2
           dr(1) = dr(1)  + ar(j+0,k) * ar(k,i) - ai(j+0,k) * ai(k,i)
           di(1) = di(1)  + ar(j+0,k) * ai(k,i) + ai(j+0,k) * ar(k,i)
           dr(2) = dr(2)  + ar(j+1,k) * ar(k,i) - ai(j+1,k) * ai(k,i)
           di(2) = di(2)  + ar(j+1,k) * ai(k,i) + ai(j+1,k) * ar(k,i)
           dr(3) = dr(3) + ai(k,j+2) * ai(k,i) + ar(k,j+2) * ar(k,i)
           di(3) = di(3) - ai(k,j+2) * ar(k,i) + ar(k,j+2) * ai(k,i)
           dr(4) = dr(4) + ai(k,j+3) * ai(k,i) + ar(k,j+3) * ar(k,i)
           di(4) = di(4) - ai(k,j+3) * ar(k,i) + ar(k,j+3) * ai(k,i)
           k = j+3
           dr(1) = dr(1)  + ar(j+0,k) * ar(k,i) - ai(j+0,k) * ai(k,i)
           di(1) = di(1)  + ar(j+0,k) * ai(k,i) + ai(j+0,k) * ar(k,i)
           dr(2) = dr(2)  + ar(j+1,k) * ar(k,i) - ai(j+1,k) * ai(k,i)
           di(2) = di(2)  + ar(j+1,k) * ai(k,i) + ai(j+1,k) * ar(k,i)
           dr(3) = dr(3)  + ar(j+2,k) * ar(k,i) - ai(j+2,k) * ai(k,i)
           di(3) = di(3)  + ar(j+2,k) * ai(k,i) + ai(j+2,k) * ar(k,i)
           dr(4) = dr(4) + ai(k,j+3) * ai(k,i) + ar(k,j+3) * ar(k,i)
           di(4) = di(4) - ai(k,j+3) * ar(k,i) + ar(k,j+3) * ai(k,i)
c     .......... form element of p ..........
           do  1240  k = 1, 4
           jj = j+k-1
           e(jj)     = dr(k) / h
           tau(jj,2) = di(k) / h
C#ifdefC SGI-PARALLEL
C           e2(jj,2) = e(jj) * ar(jj,i) + tau(jj,2) * ai(jj,i)
C#else
           f = f + e(jj) * ar(jj,i) + tau(jj,2) * ai(jj,i)
C#endif
 1240    continue
         jj = (l/4)*4+1
C#elseC
C         jj = 1
C#endif
         do  240  j = jj, l
c     .......... form element of a*u ..........
            g = 0.d0
            gi = 0.d0
            do  180  k = 1, j
               g = g + ar(k,j) * ar(k,i) + ai(k,j) * ai(k,i)
               gi = gi + ar(k,j) * ai(k,i) - ai(k,j) * ar(k,i)
  180       continue
            do  200  k = j+1, l
               g = g + ar(j,k) * ar(k,i) - ai(j,k) * ai(k,i)
               gi = gi + ar(j,k) * ai(k,i) + ai(j,k) * ar(k,i)
  200       continue
c     .......... form element of p ..........
  220       e(j) = g / h
            tau(j,2) = gi / h
C#ifdefC SGI-PARALLEL
C            e2(j,2) = e(j) * ar(j,i) + tau(j,2) * ai(j,i)
C#else
            f = f + e(j) * ar(j,i) + tau(j,2) * ai(j,i)
C#endif
  240    continue

C#ifdefC SGI-PARALLEL
C         f = 0
C         do  242  j = 1, l
C  242    f = f + e2(j,2)
C#endif
         hh = f / (h + h)

c     .......... form reduced a ..........
C#ifdefC UNROLL2
C#ifdefC SGI-PARALLEL
CC$DOACROSS local(k)
C#endifC
C         do  1260  j = 1, l-1, 2
C            dr(1) = ar(j,i)
C            dr(2) = e(j) - hh * dr(1)
C            e(j) = dr(2)
C            di(1) = ai(j,i)
C            di(2) = tau(j,2) - hh * di(1)
C            tau(j,2) = -di(2)
C
C            jj = j+1
C            dr(3) = ar(jj,i)
C            dr(4) = e(jj) - hh * dr(3)
C            e(jj) = dr(4)
C            di(3) = ai(jj,i)
C            di(4) = tau(jj,2) - hh * di(3)
C            tau(jj,2) = -di(4)
C
C            do  1262  k = 1, j
C              ar(k,j) = ar(k,j) - dr(1) * e(k) - dr(2) * ar(k,i)
C     .                          + di(1) * tau(k,2) - di(2) * ai(k,i)
C              ai(k,j) = ai(k,j) + dr(1) * tau(k,2) - dr(2) * ai(k,i)
C     .                          + di(1) * e(k) + di(2) * ar(k,i)
C              ar(k,jj) = ar(k,jj) - dr(3) * e(k) - dr(4) * ar(k,i)
C     .                          + di(3) * tau(k,2) - di(4) * ai(k,i)
C              ai(k,jj) = ai(k,jj) + dr(3) * tau(k,2) - dr(4) * ai(k,i)
C     .                          + di(3) * e(k) + di(4) * ar(k,i)
C 1262       continue
C            ar(jj,jj) = ar(jj,jj) - dr(3) * e(jj) - dr(4) * ar(jj,i)
C     .                          + di(3) * tau(jj,2) - di(4) * ai(jj,i)
C            ai(jj,jj) = ai(jj,jj) + dr(3) * tau(jj,2) - dr(4) * ai(jj,i)
C     .                          + di(3) * e(jj) + di(4) * ar(jj,i)
C 1260    continue
C         jj = (l/2)*2+1
C#else
         jj = 1
C#endif
         do  260  j = jj, l
            f = ar(j,i)
            g = e(j) - hh * f
            e(j) = g
            fi = ai(j,i)
            gi = tau(j,2) - hh * fi
            tau(j,2) = -gi
c
C#ifdefC SGI-PARALLEL
CC$DOACROSS local(k)
C#endif
            do  262  k = 1, j
              ar(k,j) = ar(k,j) - f * e(k) - g * ar(k,i)
     .                          + fi * tau(k,2) - gi * ai(k,i)
              ai(k,j) = ai(k,j) + f * tau(k,2) - g * ai(k,i)
     .                          + fi * e(k) + gi * ar(k,i)
  262       continue
  260    continue
c
  270    do  280  k = 1, l
            ar(k,i) = scale * ar(k,i)
            ai(k,i) = scale * ai(k,i)
  280    continue
c
         tau(l,2) = -si
  290    hh = d(i)
         d(i) = ar(i,i)
         ar(i,i) = hh
         ai(i,i) = -scale * dsqrt(h)
  300 continue
c
      call tcx('htridx')
      end
