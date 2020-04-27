C#define BLAS
      subroutine dpsadd(adest,asrc,nel,n1,n2,fac)
C- Shift and add. nel=number of elements, n1,n2= start in adest,asrc
      implicit none
      integer n1,n2,i,iadd,ntop,nel
      double precision asrc(1),adest(1),fac

C#ifdefC CRAY
C      call scopy(nel,asrc(n2),1,adest(n1),1)
C      call saxpy(nel,fac,asrc(n2),1,adest(n1),1)
C#elseif BLAS
      call daxpy(nel,fac,asrc(n2),1,adest(n1),1)
C#elseC   GENERIC
C      iadd = n1-n2
C      ntop = n2+nel-1
C      do  10  i = n2, ntop
C   10 adest(i+iadd) = adest(i+iadd) + fac*asrc(i)
C#endif
      end
C      program test
C      implicit none
C      double precision from(10), to(10)
C      integer i
C
C      do  10  i = 1, 10
C      to(i) = 100+i
C   10 from(i) = i
C
C      call dpsadd(to,from,4,3,5,2d0)
C      print 333, to
C  333 format(10f7.2)
C      end
