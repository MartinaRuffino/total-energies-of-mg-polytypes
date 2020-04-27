      subroutine shorpsa(nbas,plat,mode,afmt,ng,istab,pin,pout)
C- Shift each basis vector by multiples of plat according to mode
C ----------------------------------------------------------------
Ci Inputs:  nbas,plat
Ci   nbas<0: sign used as a switch to make shorps return in
Ci         array pout the change in pin, in units of plat
Ci   plat  :primitive lattice vectors, in units of alat
Ci   mode  vector of length 3 governing shifts along selected axes.
Ci         0 suppresses shifts along plat(j)
Ci         1 shifts to unit cell at origin (pos in 1st quadrant)
Ci         2 shifts to minimize length of pos
Ci         10s digit of mode(1) (1st digit mode=2 only)
Ci         If nonzero, do not shorten lattice vector unless
Ci         change in length of vector is > 10^i, where i is
Ci         10s digit.
Ci   pin:  position (basis) vectors
Co Outputs:
Co   pout  (may point to the same address space as pin).
Cr Remarks
Cr   pos = f . plat, with f = integer + fraction for each plat.
Cr   Integer part according to mode.
Cu Updates
Cu   23 Feb 09 New 10s mode in order to suppress change in basis
Ci             vectors when they are not shortened.
Cu   09 Jan 09 (A. Lozovoi) Do not change initial values of basis
Cu             vectors if they fall on the cell boundary
Cu   10 Apr 02 Patch to handle mixed boundary conditions
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nbas,mode(3),ng
      integer istab(nbas,ng+1)
      double precision plat(3,3),pin(3,nbas),pout(3,nbas),afmt(3)
C ... Local parameters
      logical latvec
      integer ib,jb
      double precision dsqrt,ddot,ppa(3),ppb(3),tol,qlat(3,3),xx
      parameter (tol=1d-12)

      call shorps(nbas,plat,mode,pin,pout)
      if (ddot(3,afmt,1,afmt,1) == 0) return

C ... qlat = (plat^-1)^T so that qlat^T . plat = 1
      call mkqlat(plat,qlat,xx)
      call dpcopy(afmt,ppa,1,3,2d0)
      if (.not. latvec(1,tol,qlat,ppa)) then
        call rx('shorpsa: SYMGRP->AFM is not a half-translation vector')
      endif
      do  ib = 1, nbas
        jb = istab(ib,ng+1)
        ppa(:) = pout(:,ib) + afmt(:) - pout(:,jb)
C       Sanity check: require p(ib)-p(jb) = afmt
        if (.not. latvec(1,tol,qlat,ppa)) then
          call rx('shorpsa: improper use of SYMGRP->AFM')
        endif
        if (jb > ib) then
C         Require p(jb) = p(ib) + afmt
          ppb(:) = pout(:,ib) + afmt(:)
C         Shorten average of the two
          ppa(:) = pout(:,ib)/2 + ppb(:)/2
          call shorps(-1,plat,mode,ppa,ppa)
          call dgemm('N','N',3,1,3,1d0,plat,3,ppa,3,0d0,ppb,3)
          print 333, ib,jb,pout(:,ib),pout(:,jb),
     .      dsqrt(ddot(3,pout(1,ib),1,pout(1,ib),1)),
     .      dsqrt(ddot(3,pout(1,jb),1,pout(1,jb),1))
          print 333, ib,jb,ppa,ppb
          pout(:,ib) = pout(:,ib) + ppb(:)
          pout(:,jb) = pout(:,ib) + afmt(:)
          print 333, ib,jb,pout(:,ib),pout(:,jb),
     .      dsqrt(ddot(3,pout(1,ib),1,pout(1,ib),1)),
     .      dsqrt(ddot(3,pout(1,jb),1,pout(1,jb),1))
  333     format(2i4,3f12.6,3f12.6,2x,2f12.6)
        endif
      enddo

      end
