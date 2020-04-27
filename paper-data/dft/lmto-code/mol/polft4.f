      subroutine polft4(f,x,val,slo,sum)
c  makes 4-th order polynomial fit to five points,
c  returns value,slope at x and integral from zero to x.
      implicit real*8 (a-h,p-z), integer(o)
      dimension a(0:4),f(0:4)
      a(0)=f(0)
      a(1)=(-25d0*f(0)+48d0*f(1)-36d0*f(2)+16d0*f(3)-3d0*f(4))/12d0
      a(2)=(35d0*f(0)-104d0*f(1)+114d0*f(2)-56d0*f(3)+11d0*f(4))/24d0
      a(3)=(-5d0*f(0)+18d0*f(1)-24d0*f(2)+14d0*f(3)-3d0*f(4))/12d0
      a(4)=(f(0)-4d0*f(1)+6d0*f(2)-4d0*f(3)+f(4))/24d0
      val=a(0)
      slo=0.d0
      sum=a(0)*x
      xjm1=1.d0
      do 41 j=1,4
      val=val+a(j)*xjm1*x
      slo=slo+a(j)*xjm1*j
      sum=sum+a(j)*xjm1*x*x/(j+1)
  41  xjm1=xjm1*x


      return
      end
