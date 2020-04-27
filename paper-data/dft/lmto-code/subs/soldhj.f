      subroutine soldhj(r,e,loka,lmax,hl,bl)
C- Real solid hankel and bessel functions.
C ----------------------------------------------------------------
Ci Inputs
Ci   r     :radius (or radius/avw using OKA's conventions)
Ci   e     :energy (or energy*avw**2 using OKA's conventions)
Ci   loka  :conventions for Hankels, Bessels; see besslr
Ci   lmax  :maximum l for a given site
Co Outputs
Co   HL,BL: Hankel and Bessel functions:  calculated to (lmax+1)**2
Cr Remarks
Cr   Generates bessel function * real spherical harmonic
Cr   MSM's standard defs, notes IV-43.
Cu Updates
Cu   19 May 04 Changed loka from logical to integer
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer lmax,loka
      double precision e
      double precision bl(*),hl(*),r(3)
C Local parameters
      integer ilm,l,m,ll
      double precision rfac,r2
      double precision phi(0:10),psi(0:10)
      external ll

C     call sylm(r,hl,lmax,r2)
      call ropyln(1,r(1),r(2),r(3),lmax,1,hl,r2)
      call besslr(e*r2,loka,0,lmax,phi,psi)
      ilm = 0
      rfac = dsqrt(r2)
      if (r2 < 1.d-10) r2 = 1.d0
      do  l = 0, lmax
C       rfac = 1/r**(2l+1), or 1/(r/w)**(2l+1) using OKA conventions
        rfac = rfac/r2
        do  m = -l, l
        ilm = ilm+1
C       xx = cy(ilm)*hl(ilm)
C       bl(ilm) = phi(l)*xx
C       hl(ilm) = (rfac*psi(l))*xx
        bl(ilm) = phi(l)*hl(ilm)
        hl(ilm) = (rfac*psi(l))*hl(ilm)
        enddo
      enddo
      end
