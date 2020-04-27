C#define BLAS
      subroutine dpscop(afrom,ato,nel,n1,n2,fac)
C- shift and copy.
Ci nel number of elements
Ci n1  offset in afrom
Ci n2  ofset in ato
      implicit none
      integer n1,n2,i,iadd,ntop,nel
      double precision afrom(1),ato(1),fac

      if (fac /= 1d0) goto 100
C#ifdefC CRAY
C      call scopy(nel,afrom(n1),1,ato(n2),1)
C#elseifC APOLLO | HP
C        call vec_$dcopy(afrom(n1),ato(n2),nel)
C#elseif BLAS
      call dcopy(nel,afrom(n1),1,ato(n2),1)
C#elseC   GENERIC
C      iadd = n2-n1
C      ntop = n1+nel-1
C      do  10  i = n1, ntop
C   10 ato(i+iadd) = afrom(i)
C#endif
      return

C --- fac not unity ---
  100 continue
      iadd = n2-n1
      ntop = n1+nel-1
      do  110  i = n1, ntop
  110 ato(i+iadd) = fac*afrom(i)
      end
C      program test
C      implicit none
C      double precision from(10), to(10)
C      integer i
C
C      do  10  i = 1, 10
C      to(i) = 0
C   10 from(i) = i
C
C      call dpscop(from,to,4,3,5,2d0)
C      print 333, to
C      call dpscop(from,to,4,3,5,1d0)
C      print 333, to
C  333 format(10f7.3)
C      end
