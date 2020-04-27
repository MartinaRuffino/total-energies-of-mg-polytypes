      subroutine dfphiz(de1,del1,de2,del2,nr,g,gp,fivep)
C- Finite-difference energy differentiation of g, complex energy
C ----------------------------------------------------------------------
Ci Inputs
Ci   de1   :see Remarks
Ci   del1  :see Remarks
Ci   de2   :see Remarks
Ci   del2  :see Remarks
Ci   nr    :number of mesh points
Ci   g     :normalized wave function times r
Ci   fivep :if true, five-point formula; else three-point formula
Co Outputs
Co   gp(i) contains estimate for ith derivative of g
Cr Remarks
Cr   This is a complex analog of dfphi, which see.
C  ----------------------------------------------------------------
      implicit none
C     Passed parameters
      integer nr
      logical fivep
      double complex de1,del1,de2,del2,g(nr),gp(nr,4)
C     Local parameters
      integer i
      double complex gp5p,gpp5p,gppp,gpp3p,gpp32
      double complex xx1,xx2,xx3,xx4,gp3p,w01,w11,w21,w02,w12,w22,
     .  w01d,w11d,w21d,w02d,w12d,w22d,wp1d,wp2d,gpppp

C --- Constants common to 3-point and 5-point formulas ---
      w01 = 2*del1/(de1**2-del1**2)
      w11 = (de1-del1)/(de1+del1)/(2*de1)
      w21 = (de1+del1)/(de1-del1)/(2*de1)
      w01d = -2/(de1**2-del1**2)
      w11d = 1/(de1+del1)/de1
      w21d = 1/(de1-del1)/de1

      if (.not. fivep) goto 20
      if (abs(del1/de1) > .1 .or.  abs(del1/de1) > .1) then
C       if (iprint() >= 20) print *, 'dfphi:  large del; use 3 point'
        goto 20
      endif

C --- Extra constants for 5-point formula ---
      xx1 = de1**2 - del1**2
      xx2 = de2**2 - del2**2
      xx3 = de1**2 + 3*del1**2
      xx4 = de2**2 + 3*del2**2
      w02 = 2*del2/(de2**2-del2**2)
      w12 = (de2-del2)/(de2+del2)/(2*de2)
      w22 = (de2+del2)/(de2-del2)/(2*de2)
      wp1d = 2d0*del1/3
      w02d = -2/(de2**2-del2**2)
      w12d = 1/(de2+del2)/de2
      w22d = 1/(de2-del2)/de2
      wp2d = 2d0*del2/3

      do  10  i = 1, nr

C Three point formula for gp; store in temporary gp3p
        gp3p = w01*g(i) + w11*gp(i,1) - w21*gp(i,2)

C Five point estimate for gp
        gp5p = (xx2*gp3p - xx1*(w02*g(i) + w12*gp(i,3) - w22*gp(i,4)))
     .          /(xx2-xx1)

C Difference between five point and three point gives estimate for gppp
        gppp = -6/xx1*(gp5p - gp3p)

C Three point estimates for gpp with correction for gppp
        gpp3p = w01d * g(i) + w11d * gp(i,1) + w21d*gp(i,2)
     .         -wp1d * gppp
        gpp32 = w02d * g(i) + w12d * gp(i,3) + w22d*gp(i,4)
     .         -wp2d * gppp

C Five point estimate for gpp with correction for gppp
        gpp5p = (gpp3p*xx4 - gpp32*xx3) / (xx4 - xx3)

C Difference between five point and three point gives est for gpppp
        gpppp = -12/xx3*(gpp5p - gpp3p)

        gp(i,1) = gp5p
        gp(i,2) = gpp5p
        gp(i,3) = gppp
        gp(i,4) = gpppp
   10 continue
      return

C Three point formulae:  only gp, gpp calculated
   20 continue
      do  30  i = 1, nr
        gp3p    = w01*g(i)  + w11*gp(i,1)  - w21*gp(i,2)
        gp(i,2) = w01d*g(i) + w11d*gp(i,1) + w21d*gp(i,2)
        gp(i,1) = gp3p
   30 continue

      end
