c --------- gpfnd1 --------------------------------
      subroutine gpfnd1(g,ag,p,nbas,pos,d,jb)
c  finds atom jb at map of p under group operation (g,ag).
      implicit real*8 (a-h,p-z), integer(o)
      dimension g(3,3),ag(3),pos(3,1),d(3),p(3)
      do 4 m=1,3
      d(m)=ag(m)
      do 4 n=1,3
  4   d(m)=d(m)+g(m,n)*p(n)
      jb=0
      do 5 k=1,nbas
      ddd=(pos(1,k)-d(1))**2+(pos(2,k)-d(2))**2+(pos(3,k)-d(3))**2
  5   if(dsqrt(ddd) < 1d-5) jb=k
      end
c --------- gpfnd2 --------------------------------
      subroutine gpfnd2(g,ag,p,nbas,pos,jb)
c  finds atom jb which is transformed into pos by operation (g,ag).
      implicit real*8 (a-h,p-z), integer(o)
      dimension g(3,3),ag(3),pos(3,1),d(3),p(3)
      jb=0
      do 1 kb=1,nbas
      do 2 m=1,3
      d(m)=ag(m)-p(m)
      do 2 k=1,3
  2   d(m)=d(m)+g(m,k)*pos(k,kb)
      if(d(1)**2+d(2)**2+d(3)**2 < 1.d-8) then
        jb=kb
        return
        endif
  1   continue
      return
      end
