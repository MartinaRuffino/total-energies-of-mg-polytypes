      subroutine han1ci(e1,e2,r,lmax,s0)
c  Integrals of products of hankel fcts on same site,
c  heads subtracted out
      implicit real*8 (a-h,p-z), integer(o)
      dimension s0(0:1),fkj1(10),fkj2(10),fkk(10),fjj(10),fjk(10)

c --- for integrals outside spheres ---
      call wronkj(e1,e2,r,lmax,fkk,fkj1,fjk,fjj)
      do 3 l=0,lmax
    3 s0(l)=fkk(l+1)

c      print 333, (s0(l),l=0,lmax)
c  333 format(9f8.5)
      end
