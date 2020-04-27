c  transformation and back-transformation for cheyshev-fit.
      subroutine ropg0(a,x,y)
      implicit real*8 (a-h,p-z)
      dp=dexp(a*x)
      dm=1d0/dp
      y=(dp-dm)/(dp+dm)
      y=y*y
      end

      subroutine ropgi(a,y,x)
      implicit real*8 (a-h,p-z)
      u=dsqrt(y)
      x=dlog((1d0+u)/(1d0-u))/(2d0*a)
      end
