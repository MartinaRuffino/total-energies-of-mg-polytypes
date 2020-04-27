      subroutine setylm(nlm,n,xyz,cy,yl)
c  Spheric harmonic polynomials for list of point
      implicit real*8 (a-h,p-z), integer(o)
      dimension xyz(3,1),yl(nlm,1),cy(1)
      if(n <= 0.or.nlm <= 0) return
      lmax=ll(nlm)
      do 2 i=1,n
      call sylm(xyz(1,i),yl(1,i),lmax,r2s)
      do 3 ilm=1,nlm
  3   yl(ilm,i)=yl(ilm,i)*cy(ilm)
  2   continue
      end
