      subroutine dpdist(a,b,n,d)
      implicit real*8 (a-h,p-z), integer(o)
      dimension a(n),b(n)
      sum=0d0
      do 10 i=1,n
  10  sum=sum+(a(i)-b(i))**2
      d=dsqrt(sum)
      end
