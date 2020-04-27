      subroutine rophs1(lmax,c,nc0,nc,ww,u2,xi,ndim,n)
c  Makes vectors of smoothed Hankel functions, l=0..lmax
c  Here: Chebyshev polnomials for all l.
      implicit real*8 (a-h,p-z), integer (o)
      dimension c(nc0,1),xi(ndim,0:lmax),u2(n),ww(3*n)
      if(lmax < 0.or.n <= 0) return
      do 10 l=0,lmax
      call ropech(c(1,l+1),nc,n,u2,xi(1,l),ww)
  10  continue
      end
