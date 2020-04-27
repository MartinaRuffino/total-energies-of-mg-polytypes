      subroutine dpzero(array,leng)
C- Initializes double precision array to zero
      implicit none
      integer leng,i
      double precision array(leng)
      forall (i=1:leng) array(i) = 0
      end
C
C      subroutine dpzero(a,n)
C      real*8 a(n)
C      do 1 i=1,n
C        a(i)=0.d0
C    1 continue
C      return
C      end
