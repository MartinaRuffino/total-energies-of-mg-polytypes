C#define BLAS
      subroutine dpdot(a,b,n,sum)
      implicit none

      integer n,i
      double precision a(1),b(1),sum
C#ifdefC APOLLO | HP
C      double precision vec_$ddot
C      sum = vec_$ddot(a,b,n)
C#elseifC CRAY
C      sum = sdot(n,a,1,b,1)
C#elseif BLAS
      double precision ddot
      sum = ddot(n,a,1,b,1)
C#elseC
C      sum = 0d0
C      do  10  i = 1, n
C   10 sum = sum + a(i)*b(i)
C#endif
      end
