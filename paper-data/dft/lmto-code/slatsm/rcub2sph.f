      subroutine rcub2sph(forward,u)
C- Returns 2x2 rotation matrices for rotating cubic to spherical harmonics
C ----------------------------------------------------------------------
Ci Inputs
Ci  forward:T => Returns u as described in Remarks.     rotates Real -> Spherical.
Ci         :F => Returns u^-1 as described in Remarks.  rotates Spherical -> Real.
Co Outputs
Co  u     :u(:,:,1) rotation for m odd
Co        :u(:,:,2) rotation for m even
Co        :See Eq. (9) on the web page for definition and uses of u, and also Remarks
Cr Remarks
Cr   Standard definition (https://www.questaal.org/docs/numerics/spherical_harmonics)
Cr   Call R_l,m and Y_l,m the real and spherical harmonics respectively.
Cr   Rotate a function expressed in Y_lm to basis of R_lm.
Cr   Since l's do not mix, it is necessary to consider rotation only between (-m,m), or (-,+)
Cr
Cr   f(r) = Sum c_m Y_m(r) = Sum r_nu R_nu(r)
Cr
Cr   Relation is c_m = sum_nu u_m,nu r_nu (see web site). where for m odd
Cr
Cr                   1     (+i   1)                  R-mu = (i*Y-m + i*Ym)/sqrt(2) (r- = 1, r+ = 0)
Cr    (m odd)  u = -----   (      )      c=ur =>
Cr                 sqrt(2) (+i  -1)                  Rmu  = (Y-m - Ym)/sqrt(2)     (r- = 0, r += 1)
Cr
Cr             1     (-i  -i)   (r-)                 Y-m = (-i*R-mu + Rmu)/sqrt(2) (c -= 1, c += 0)
Cr    u^-1 = -----   (      ) = (  )     r=u^-1c =>
Cr           sqrt(2) ( i  -1)   (r+)                 Ym  = (-i*R-mu - Rmu)/sqrt(2) (c- = 0, c += 1)
Cr
Cr                   1     (+i   1)                  R-mu = (i*Y-m - i*Ym)/sqrt(2) (r- = 1, r+ = 0)
Cr   (m even)  u = -----   (      )      c=ur =>
Cr                 sqrt(2) (-i   1)                  Rmu  = (Y-m + Ym)/sqrt(2)     (r- = 0, r += 1)
Cr
Cu Updates
Cu   26 Mar 18 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical forward
      complex(8) :: u(2,2,2)
C ... Local parameters
      complex(8), parameter :: sr2 = dcmplx(1d0/dsqrt(2d0),0d0), iot = dcmplx(0d0,1d0)/dsqrt(2d0)

      if (forward) then  ! Return r^-1

        u(1,1:2,1) = [iot,  sr2] ! m odd
        u(2,1:2,1) = [iot, -sr2]
        u(1,1:2,2) = [iot,  sr2] ! m even
        u(2,1:2,2) = [-iot, sr2]


C        call zprm('R to Y, m odd',2,u(1,1,1),2,2,2)
C        call zprm('R to Y, m even',2,u(1,1,2),2,2,2)
      else ! Return r

        u(1,1:2,1) = [-iot, -iot] ! m odd
        u(2,1:2,1) = [ sr2, -sr2]
        u(1,1:2,2) = [-iot,  iot] ! m even
        u(2,1:2,2) = [ sr2,  sr2]

C        call zprm('Y to R, m odd',2,u(1,1,1),2,2,2)
C        call zprm('Y to R, m even',2,u(1,1,2),2,2,2)
      endif
      end

C     Test
C      subroutine fmain
C      double complex r2s(2,2,2),r2si(2,2,2),prd(2,2),r(2),y(2)
C
C      call rcub2sph(.true.,r2s)
C      call rcub2sph(.false.,r2si)
C      call zmpy22(r2s,r2si,prd)
C      call zprm('product m > 0',2,prd,2,2,2)
C      call zmpy22(r2s(1,1,2),r2si(1,1,2),prd)
C      call zprm('product m < 0',2,prd,2,2,2)
C
CC     Y11 = (-x - iy)/sqrt(2)
CC     u * R11 = Y11
C      r = [-(0d0,1d0), -(1d0,0d0)]/dsqrt(2d0)
C      call zgemm('N','N',2,1,2,(1d0,0d0),r2s,2,r,2,(0d0,0d0),y,2)
C      print *, y(1)
C      print *, y(2)
CC      r = [-(0d0,1d0), (1d0,0d0)]
CC      call zgemm('N','N',2,1,2,(1d0,0d0),r2s,2,r,2,(0d0,0d0),y,2)
CC      print *, y
C
C      stop
C
C      end
