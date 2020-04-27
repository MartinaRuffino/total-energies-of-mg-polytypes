      subroutine hyfloc(r1,lp1,ep1,r2,lp2,ep2,nx1,lx1,ex1,
     .   nx2,lx2,ex2,t,itbl,jxi1,jxi2)
c  locate the correct table for 2-center case
      implicit real*8 (a-h,p-z), integer (o)
      dimension t(100,1),lx1(1),ex1(1),lx2(1),ex2(1),kx1(10),
     .   fx1(10),kx2(10),fx2(10),jxi1(10),jxi2(10),jx01(10),jx02(10)
      call getpr(ipr)
      ntbl=idnint(t(26,1))
      itbl=0
c ------- start loop over tables ---------
      dr=1.d6
      do 10 it=1,ntbl
      lsym        =idnint( t(13,it) )
      if(lsym /= 0.and.lsym /= 1) goto 10
      s1          =        t( 1,it)
      s2          =        t( 2,it)
      kp1         =idnint( t(15,it) )
      fp1         =        t(16,it)
      kp2         =idnint( t(17,it) )
      fp2         =        t(18,it)
      mx1         =idnint( t(19,it) )
      mx2         =idnint( t(20,it) )
      do 14 i=1,mx1
      kx1(i)      =idnint( t(40+i,it) )
  14  fx1(i)      =        t(50+i,it)
      do 15 i=1,mx2
      kx2(i)      =idnint( t(60+i,it) )
  15  fx2(i)      =        t(70+i,it)
c ------- compare to parameters -------------
      call hyftry(r1,lp1,ep1,nx1,lx1,ex1,r2,lp2,ep2,nx2,lx2,ex2,
     .   s1,kp1,fp1,mx1,kx1,fx1,s2,kp2,fp2,mx2,kx2,fx2,
     .   jx01,jx02,lok,rdiff)
      if(lok == 1.and.rdiff < dr) then
        itbl=it
        dr=rdiff
        do 32 i=1,mx1
  32    jxi1(i)=jx01(i)
        do 33 i=1,mx2
  33    jxi2(i)=jx02(i)
        endif
c ------- if table is not symmetric, try switching sites -----
      if(lsym == 0) then
      call hyftry(r1,lp1,ep1,nx1,lx1,ex1,r2,lp2,ep2,nx2,lx2,ex2,
     .   s2,kp2,fp2,mx2,kx2,fx2,s1,kp1,fp1,mx1,kx1,fx1,
     .   jx01,jx02,lok,rdiff)
      if(lok == 1.and.rdiff < dr) then
        itbl=-it
        dr=rdiff
        do 42 i=1,mx2
  42    jxi1(i)=jx01(i)
        do 43 i=1,mx1
  43    jxi2(i)=jx02(i)
        endif
      endif
  10  continue
      return
      end
