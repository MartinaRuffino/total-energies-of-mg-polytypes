      subroutine dfphi(dw1,dwl1,dw2,dwl2,nr,g,gp,fivep)
C- Numerically differentiates g(:,w) using five point finite difference
C ----------------------------------------------------------------------
Ci Inputs
Ci   dw1   :One term in w change, 3 and 5 point formulas;  see Remarks
Ci   dwl1  :One term in w change, 3 and 5 point formulas;  see Remarks
Ci   dw2   :One term in w change, 5 point formulas only;  see Remarks
Ci   dwl2  :One term in w change, 5 point formulas only;  see Remarks
Ci   nr    :number of mesh points
Ci   g     :function on a mesh of points 1:nr, at w=0.
Ci   fivep :if true, five-point formula; else three-point formula
Cio Inputs/Outputs
Cio  gp    :On input, for each of nr points, gp(*,1:4) contains g(*) at
Cio        :  w = dw1+dwl1, -dw1+dwl1                      (three  point)
Cio        :  w = dw1+dwl1, -dw1+dwl1, dw2+dwl2, -dw2+dwl2 (five point)
Cio        :On output, gp(*,:) contains derivatives of g(*) wrt w :
Cio        :two derivatives or four derivatives, depending on fivep
Cr Remarks
Cr   A simple five-point estimate for the numerical differentiation for
Cr   the first four derivatives of phi is based on integration of phi
Cr   at four energies dw1+dwl1, -dw1+dwl1, dw2+dwl2, -dw2+dwl2.
Cr   Tthe deviations dwl1 and dwl2 are only nonzero to second order in
Cr   he increment de, and this can be exploited to obtain accurate five point
Cr   estimates without having to solve for the four simultaneous equations.
Cr
Cr   A three point for differences wrt about w=0, ep = de+dwl,
Cr   em = -de+dwl gives
Cr
Cr   (1)  gp =  2 dwl/(de^2-dwl^2) g(0) +
Cr              (de-dwl)/(de+dwl)/(2 de) g(ep) -
Cr              (de+dwl)/(de-dwl)/(2 de) g(em) +
Cr
Cr              1/6 gppp (de^2-dwl^2) + 1/12 gpppp dwl (de^2-dwl^2) + ...
Cr
Cr
Cr   (2) gpp = -2/(de^2-dwl^2) g(0) +
Cr              1/(de+dwl)/de g(ep) +
Cr              1/(de-dwl)/de g(em) -
Cr
Cr              2/3 gppp dwl + 1/12 gpppp dwl (de^2 + 3 dwl^2) + ...
Cr
Cr
Cr   The gppp term in (1) can be knocked out by taking two three point
Cr   formulas in linear combination (dw1^2-dwl1^2) and (dw2^2-dwl2^2)
Cr   leaving only fourth and higher order terms.  Because dwl is of
Cr   order de^2, the fourth order term is of the same order as the
Cr   fifth order term and there is no advantage in eliminating it.
Cr   Also from the difference between this more accurate (5-point)
Cr   estimate for gp, an estimate for gppp can be made as
Cr
Cr   (3) gppp = -6 ( gp(five point) - gp(three point) ) /(dw1^2-dwl1^2)
Cr
Cr             + 1/2 gpppp dwl1
Cr
Cr   which is again accurate to order de^2.  Once gppp is known to this
Cr   order the term proportional to gppp in three point estimate for
Cr   gpp can be subtracted out directly and the gpppp term can be
Cr   eliminated by taking three point formulas (with the gppp term
Cr   removed) in linear combinations (dw1^2+3*dwl1^2) and
Cr   (dw2^2+3*dwl2^2), and finally the fourth order derivative can be
Cr   estimated from the difference in the five-point estimate for gpp
Cr   and the three point estimate.
C  ----------------------------------------------------------------
      implicit none
C     Passed parameters
      integer nr
      logical fivep
      double precision dw1,dwl1,dw2,dwl2,g(nr),gp(nr,4)
C     Local parameters
      integer i
      double precision gp5p,gpp5p,gppp,gpp3p,gpp32,gp3p,gpppp
      double precision dx(4),w3(0:2),w3p(0:2),w5(0:3),w5p(0:3)

      real(8) :: dw1i, dw2i

C --- Constants common to 3-point and 5-point formulas ---
! -march=native on avx2 and avx512 cpus may give more accurate results because this algorithm is rather sensitive.
! To stabilise it a little and improve consistency between vectorised and scalar compilations 1/dw1 was exported and (dw1**2 - dwl1**2) replaced with (dw1+dwl1)*(dw1-dwl1).

      dw1i = 1/dw1

      w3(0) = 2*dwl1/((dw1+dwl1)*(dw1-dwl1))
      w3(1) = 0.5_8*dw1i*(dw1-dwl1)/(dw1+dwl1)
      w3(2) = 0.5_8*dw1i*(dw1+dwl1)/(dw1-dwl1)
      w3p(0) = -2/((dw1+dwl1)*(dw1-dwl1))
      w3p(1) = dw1i/(dw1+dwl1)
      w3p(2) = dw1i/(dw1-dwl1)

      if (.not. fivep) goto 20
      if (dabs(dwl1/dw1) > .1 .or.  dabs(dwl2/dw2) > .1) then
C       if (iprint() >= 20) print *, 'dfphi:  large dwl; use 3 point'
        goto 20
      endif

C --- Extra constants for 5-point formula ---
      dw2i = 1/dw2

      dx(1) = (dw1+dwl1)*(dw1-dwl1)
      dx(2) = (dw2+dwl2)*(dw2-dwl2)
      dx(3) = dw1**2 + 3*dwl1**2
      dx(4) = dw2**2 + 3*dwl2**2
      w5(0) = 2*dwl2/((dw2+dwl2)*(dw2-dwl2))
      w5(1) = 0.5_8*dw2i*(dw2-dwl2)/(dw2+dwl2)
      w5(2) = 0.5_8*dw2i*(dw2+dwl2)/(dw2-dwl2)
      w5(3) = 2d0*dwl1/3
      w5p(0) = -2/((dw2+dwl2)*(dw2-dwl2))
      w5p(1) = dw2i/(dw2+dwl2)
      w5p(2) = dw2i/(dw2-dwl2)
      w5p(3) = 2d0*dwl2/3

      do  i = 1, nr

C Three point formula for gp; store in temporary gp3p
        gp3p = w3(0)*g(i) + w3(1)*gp(i,1) - w3(2)*gp(i,2)

C Five point estimate for gp
        gp5p = (dx(2)*gp3p - dx(1)*(w5(0)*g(i) + w5(1)*gp(i,3) - w5(2)*gp(i,4)))/(dx(2)-dx(1))

C Difference between five point and three point gives estimate for gppp
        gppp = -6/dx(1)*(gp5p - gp3p)

C Three point estimates for gpp with correction for gppp
        gpp3p = w3p(0)*g(i)+w3p(1)*gp(i,1)+w3p(2)*gp(i,2)-w5(3)*gppp
        gpp32 = w5p(0)*g(i)+w5p(1)*gp(i,3)+w5p(2)*gp(i,4)-w5p(3)*gppp

C Five point estimate for gpp with correction for gppp
        gpp5p = (gpp3p*dx(4) - gpp32*dx(3)) / (dx(4) - dx(3))

C Difference between five point and three point gives est for gpppp
        gpppp = -12/dx(3)*(gpp5p - gpp3p)

        gp(i,1) = gp5p
        gp(i,2) = gpp5p
        gp(i,3) = gppp
        gp(i,4) = gpppp
      enddo
      return

C ... Three point formulae:  only gp, gpp calculated
   20 continue
      do  i = 1, nr
        gp3p    = w3(0)*g(i)  + w3(1)*gp(i,1)  - w3(2)*gp(i,2)
        gp(i,2) = w3p(0)*g(i) + w3p(1)*gp(i,1) + w3p(2)*gp(i,2)
        gp(i,1) = gp3p
      enddo

      end

C     Test
C      subroutine fmain
C      implicit none
C      double precision x0,x1,x2,f0,f1,f2,fp(2),dw1,dwl1,fpp
C
C      x0 = 1.5d0
C      x1 = 1d0
C      x2 = 0d0
C
C      fpp = 18
C
C      f0 = 3 + 4*x0 + fpp/2*x0**2
C      f1 = 3 + 4*x1 + fpp/2*x1**2
C      f2 = 3 + 4*x2 + fpp/2*x2**2
C
CC      f0 = 0.0003944900053837638d0
CC      f1 = 0.0015489281763727354d0
CC      f2 = 0.006138053346620949d0
C
CC     Get f', f'' around x0 from finite difference
C      dw1 = (x1 - x2)/2; dwl1 = (x1 + x2)/2 - x0
C      fp(1) = f1; fp(2) = f2
C      print *, 'x+, x-'
C      print *, x0 + dw1 + dwl1, x0 - dw1 + dwl1
C
C      call dfphi(dw1,dwl1,dw1,dwl1,1,f0,fp,.false.)
C      print *, fp
C      print *, 'should get'
C      print *, 4d0+fpp*x0, fpp
C      end
