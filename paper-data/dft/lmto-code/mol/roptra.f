
      subroutine roptra(r1,r2,a,n,rsq,u2)
c  makes vector of transformed distances for ropech.
      implicit real*8 (a-h,p-z), integer(o)
      dimension rsq(n),u2(n)
      if(n <= 0) return
      call ropg0(a,r1,y1)
      call ropg0(a,r2,y2)
      bma=0.5d0*(y2-y1)
      bpa=0.5d0*(y2+y1)
      a1=2d0/bma
      a2=-2d0*bpa/bma
      do 10 i=1,n
      xx=dexp((-2d0*a)*dsqrt(rsq(i)))
      yy=(1d0-xx)/(1d0+xx)
  10  u2(i)=a1*(yy*yy)+a2
      end
