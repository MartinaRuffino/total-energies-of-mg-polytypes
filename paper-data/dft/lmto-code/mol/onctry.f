      subroutine onctry(r,lp1,ep1,lp2,ep2,nx,lx,ex,
     .   s,kp1,fp1,kp2,fp2,mx,kx,fx,lmax,jxi,lok,rdiff)
c  Compares parameters for which the 1-c-f is needed with
c  specifications of a table. returns lok=1 if table is ok.
c  Then returns pointers for the xi-fcts and radii difference.
      implicit real*8 (a-h,p-z), integer (o)
      dimension lx(1),ex(1),kx(10),fx(10),jxi(10)
      call getpr(ipr)
      fuzz=1.d-10
      lok=0
      rdiff=0
c ------- compare radii and phi's -------------
      if(s-r > fuzz) return
      if(dabs(ep1-fp1) > fuzz) return
      if(dabs(ep2-fp2) > fuzz) return
      if(kp1 < lp1.or.kp2 < lp2) return
c ------- look for table-xi's in those of calc -------
      do 21 i=1,mx
      j0=0
      do 22 j=1,nx
      lmux=min0(lmax,kx(i))
  22  if(lx(j) >= lmux.and.dabs(ex(j)-fx(i)) < fuzz) j0=j
      if(j0 == 0) return
  21  jxi(i)=j0
c ------- here if the table can be used ------------
      rdiff=r-s
      lok=1
      return
      end
