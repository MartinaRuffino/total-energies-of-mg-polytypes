      subroutine oncloc(r,lp1,ep1,lp2,ep2,nx,lx,ex,t,itbl,jxi)
c  locate the correct table for the 1-center case
      implicit real*8 (a-h,p-z), integer (o)
      dimension t(100,1),lx(1),ex(1),kx(10),fx(10),jxi(10),jx0(10)
      fuzz=1d-10
      call getigv(1,ipr)
      ntbl=idnint(t(26,1))
      itbl=0
c ------- start loop over tables ---------
      dr=1.d6
      do 10 it=1,ntbl
      lsym        =idnint( t(13,it) )
      if(lsym /= 3) goto 10
      s           =        t( 1,it)
      lmax        =idnint( t( 8,it) )
      kp1         =idnint( t(15,it) )
      fp1         =        t(16,it)
      kp2         =idnint( t(17,it) )
      fp2         =        t(18,it)
      mx          =idnint( t(19,it) )
      do 14 i=1,mx
      kx(i)       =idnint( t(40+i,it) )
  14  fx(i)       =        t(50+i,it)
c ------- compare to parameters -------------
      call onctry(r,lp1,ep1,lp2,ep2,nx,lx,ex,
     .   s,kp1,fp1,kp2,fp2,mx,kx,fx,lmax,jx0,lok,rdiff)
      if(lok == 1.and.rdiff < dr) then
        itbl=it
        dr=rdiff
        do 31 i=1,mx
  31    jxi(i)=jx0(i)
        endif
c ------- if table is not symmetric, try switching phi's -----
      if(kp1 /= kp2.or.dabs(fp1-fp2) > fuzz) then
      call onctry(r,lp1,ep1,lp2,ep2,nx,lx,ex,
     .   s,kp2,fp2,kp1,fp1,mx,kx,fx,lmax,jx0,lok,rdiff)
      if(lok == 1.and.rdiff < dr) then
        itbl=-it
        dr=rdiff
        do 32 i=1,mx
  32    jxi(i)=jx0(i)
        endif
      endif
  10  continue
      return
      end
