      subroutine hanr(rsq,lmin,lmax,nrx,nr,e,xi)
C- Vector of unsmoothed hankel functions for l=0...lmax, negative e.
C ---------------------------------------------------------------
Ci Inputs
Ci   rsq,nr  vector of points r**2, and number of points.
Ci   e       energy of Hankel.
Ci   lmin:lmax generate xi from lmin:lmax.  lmin must be -1 or 0.
Ci   nrx     leading dimension of xi.
Co Outputs
Co   xi      radial Hankel functions/r**l for: xi(1..nr, lmin:lmax)
Co           Solid hankel function is hl(ilm) = xi(l)*cy(ilm)*yl(ilm)
Co           where yl(ilm) are unnormalized sph. harm. polynomials
Co           Energy derivative is    hlp(ilm) = x(l-1)/2*cy(ilm)*yl(ilm)
C ---------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nrx,nr,lmin,lmax
      double precision rsq(nr),e,xi(nrx,lmin:lmax)
C ... Local parameters
      double precision sre,r2,r,h0,xx,akap
      integer l,ir

      if (lmin /= 0 .and. lmin /= -1)
     .  call rx('hanr: input lmin must be -1 or 0')
      if (lmax < lmin .or. nr <= 0) return
      akap = dsqrt(-e)

C --- Make xi(lmin), xi(lmin+1) ---
C ... xi(-1) for lmax=-1 only
      if (lmin == -1 .and. lmax == -1) then
        do  ir = 1, nr
          r = dsqrt(rsq(ir))
          sre = akap*r
          h0 = dexp(-sre)/r
          xi(ir,-1) = -h0*sre/e
        enddo
      elseif (lmin == -1) then
        do  ir = 1, nr
          r = dsqrt(rsq(ir))
          sre = akap*r
          h0 = dexp(-sre)/r
          xi(ir,0) = h0
          xi(ir,-1) = -h0*sre/e
        enddo
C ... xi(0) for lmax=0 only
      elseif (lmax == 0) then
        do  ir = 1, nr
          r2 = rsq(ir)
          r = dsqrt(r2)
          xi(ir,0) = dexp(-akap*r)/r
        enddo
      else
        do  ir = 1, nr
          r = dsqrt(rsq(ir))
          sre = akap*r
          h0 = dexp(-sre)/r
          xi(ir,0) = h0
          xi(ir,1) = h0*(1d0+sre)/rsq(ir)
        enddo
      endif

C --- xi(*,lmin+2:lmax) by upward recursion ---
      do  l = lmin+2, lmax
        xx = 2*l-1
        do  ir = 1, nr
          xi(ir,l) = (xx*xi(ir,l-1)-e*xi(ir,l-2))/rsq(ir)
        enddo
      enddo

      end
