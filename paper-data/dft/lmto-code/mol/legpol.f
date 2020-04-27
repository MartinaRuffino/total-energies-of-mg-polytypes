      subroutine legpol(n,x,p)
C  Legendre polynomials p_0(x) ... p_n(x)
      implicit real*8 (a-h,p-z),integer (o)
      dimension p(0:n)
      p(0)=1d0
      if(n >= 1) p(1)=x
      do 1 m=2,n
  1   p(m)=((2*m-1)*x*p(m-1)-(m-1)*p(m-2))/m
      end
