      subroutine hyfsch(ndist,nalf,d0,adec,wx,dist,wist,f)
C- Sets Chebyshev polynomials and distances
      implicit real*8 (a-h,p-z), integer (o)
      dimension dist(1),wist(1),f(ndist,1)
      call getpr(ipr)
      pi=4d0*datan(1d0)
      if(ipr >= 20) write(6,200) adec,d0,ndist,nalf
  200 format(/' ftcsch:  adec=',f7.3,'   d0=',f7.3,'   ndist=',i3,
     .   '   nalf=',i3)
      if(ndist == 1) then
        dist(1)=d0
        wist(1)=1d0
        f(1,1)=1d0
        return
        endif
      if(ipr >= 30) write(6,782)
c ------ start loop over distances ---------
      do 10 k=1,ndist
      y=dcos(pi*(k-0.5d0)/ndist)
      x=0.5d0*(y+1d0)
      dist(k)=d0-dlog(x)/adec
      wist(k)=1d0
      if (wx /= 0d0) wist(k) = x**wx
      do 12 j=1,nalf
  12  f(k,j)=dcos((pi*(j-1))*((k-0.5d0)/ndist))
      lp=0
      if(ipr >= 30.and.(k == 1.or.k == ndist)) lp=1
      if(ipr >= 40) lp=1
      if(lp == 1) write(6,781) k,y,x,dist(k),wist(k)
  781 format(i5,4f12.6)
  782 format(/'   id',6x,'y',11x,'x',11x,'dist',8x,'wist')
  10  continue
c ------- printout cheby polynomials ------------
      if(ipr < 40) return
      write(6,334)
      do 30 k=1,ndist
  30  write(6,333) k,dist(k),(f(k,j),j=1,nalf)
  333 format(i5,f11.4,1x,5f11.5:/(17x,5f11.5))
  334 format(/'  idist     d        Chebyshev polynomials ..')
      end
