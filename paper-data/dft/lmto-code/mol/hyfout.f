      subroutine hyfout(r1,r2,rsm1,rsm2,d0,adec,a1,ndist,nalf,mmax,
     .   ncof,n1,n2,ri1,ri2,nxi,lp1,ep1,lp2,ep2,lsym,lx1,ex1,nx1,
     .   lx2,ex2,nx2,cof,err,ifi)
c  writes tcf specifications and data into file ifi
      implicit real*8 (a-h,p-z), integer (o)
      dimension cof(ncof,nalf),lx1(8),ex1(8),lx2(8),ex2(8),t(100),
     .   err(20)
      data t / 100*0d0 /
      if(max0(nx1,nx2) > 10) call rx('hyfout: nx too large''')
      call getpr(ipr)
c ------- put info into vector t ---------------
      t( 1)=r1
      t( 2)=r2
      t( 3)=d0
      t( 4)=adec
      t( 5)=a1
      t( 6)=ndist
      t( 7)=nalf
      t( 8)=mmax
      t( 9)=ncof
      t(10)=n1
      t(11)=n2
      t(12)=nxi
      t(13)=lsym
      t(15)=lp1
      t(16)=ep1
      t(17)=lp2
      t(18)=ep2
      t(19)=nx1
      t(20)=nx2
      t(21)=rsm1
      t(22)=rsm2
      t(23)=ri1
      t(24)=ri2
      do 10 i=1,nx1
      t(40+i)=lx1(i)
  10  t(50+i)=ex1(i)
      do 11 i=1,nx2
      t(60+i)=lx2(i)
  11  t(70+i)=ex2(i)
      nerr=(lp1+1)*(lp1+1)
      if(nerr > 20) call rx('hyfout:  array err too long''')
      do 12 i=1,nerr
  12  t(i+80)=err(i)
c ------- write into file ----------------------
      if(ipr >= 10) write(6,100) ncof,nalf,ifi
  100 format(/' hyfout:  ncof=',i6,'   nalf=',i3,'   file',i3)
      if(ipr >= 60) write(6,200) t
  200 format(1x,10f8.2)
      call dpdump(t,100,-ifi)
      call dpdump(cof,ncof*nalf,-ifi)

      return
      end
