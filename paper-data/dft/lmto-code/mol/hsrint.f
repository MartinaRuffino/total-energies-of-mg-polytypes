      subroutine hsrint(xp,yp,wp,wk,wk2,idx,n)
C- Sorts x,y pairs by decreasing length
      implicit none

      integer n,idx(n),i
      double precision xp(1),yp(1),wp(1),wk(2,n),wk2(n,3)


      do  10  i = 1, n
        wk(1,i) = xp(i)
        wk(2,i) = yp(i)
   10 continue

      call dvshel(2,n,wk,idx,11)

      do  20  i = 1, n
        wk2(n-i+1,1) = xp(idx(i)+1)
        wk2(n-i+1,2) = yp(idx(i)+1)
        wk2(n-i+1,3) = wp(idx(i)+1)
   20 continue

      do  30  i = 1, n
        xp(i) = wk2(i,1)
        yp(i) = wk2(i,2)
        wp(i) = wk2(i,3)
   30 continue

C      do  40  i = 1, n
C        print 345, i, xp(i), yp(i), xp(i)**2+yp(i)**2
C  345   format(i5,3f20.10)
C   40 continue
C      stop

      end
