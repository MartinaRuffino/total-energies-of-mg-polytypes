      subroutine hyfchk(nspec,r,rsm,nphi,lphi,ephi,nxi,lxi,exi,
     .   slab,n0,t)
c  checks if tables are available for all needed cases
      implicit real*8 (a-h,p-z), integer (o)
      dimension t(100,1),r(1),nphi(1),lphi(n0,1),ephi(n0,1),nxi(1),
     .   lxi(n0,1),exi(n0,1),jxi1(10),jxi2(10),fx1(10),fx2(10),
     .   lnot1(6),lnot2(6),lnot(6),kx1(10),kx2(10),rsm(1)
      character*4 slab(1)
      character*7 cnot1,cnot2
      character*15 ccc
      call getpr(ipr)
      do 2 k=1,6
  2   lnot(k)=0
      write(6,209)

c ------- printout data for species -----------------
      if(ipr >= 30) then
      write(6,448)
      do 3 i=1,nspec
      write(ccc,'(3f5.1)') (ephi(ie,i),ie=1,nphi(i))
      lx=0
      do 7 j=1,nxi(i)
  7   lx=10*lx+lxi(j,i)
      lp=0
      do 8 j=1,nphi(i)
      lj=lphi(j,i)
      if(lj < 0) lj=9
  8   lp=10*lp+lj
  3   write(6,449) slab(i),r(i),lp,ccc,lx,(exi(j,i),j=1,nxi(i))
  449 format(1x,a4,f7.3,i6,1x,a15,i7,6f6.1)
  448 format(/' spec    rmt  lphi   ephi',13x,'lxi   exi')
      write(6,*) ' '
      endif

c ------- start loop over cases for which tables are needed ----
      write(6,211)
      write(6,210)
      do 10 i=1,nspec
      do 10 j=1,i
      do 10 ie=1,nphi(i)
      jetop=nphi(j)
      if(j == i) jetop=ie
      do 10 je=1,jetop
      r1=r(i)
      r2=r(j)
      rm1=rsm(i)
      rm2=rsm(j)
      nx1=nxi(i)
      nx2=nxi(j)
      e1=ephi(ie,i)
      e2=ephi(je,j)
      lp1=lphi(ie,i)
      lp2=lphi(je,j)
      if(lp1 < 0.or.lp2 < 0) goto 10
      call hyfloc(r1,lp1,e1,r2,lp2,e2,nx1,lxi(1,i),exi(1,i),
     .   nx2,lxi(1,j),exi(1,j),t,itbl,jxi1,jxi2)

c ------- get table data ------------------------
      it=iabs(itbl)
      if(it > 0) then
      s1          =        t( 1,it)
      s2          =        t( 2,it)
      lsym        =idnint( t(13,it) )
      kp1         =idnint( t(15,it) )
      fp1         =        t(16,it)
      kp2         =idnint( t(17,it) )
      fp2         =        t(18,it)
      mx1         =idnint( t(19,it) )
      mx2         =idnint( t(20,it) )
      sm1         =        t(21,it)
      sm2         =        t(22,it)
      ll1=0
      do 4 k=1,mx1
      kx1(k)=    idnint( t(40+k,it) )
      ll1=10*ll1+idnint( t(40+k,it) )
  4   fx1(k)      =      t(50+k,it)
      ll2=0
      do 5 k=1,mx2
      kx2(k)=    idnint( t(60+k,it) )
      ll2=10*ll2+idnint( t(60+k,it) )
  5   fx2(k)      =      t(70+k,it)
      endif

c ------- printout ----------------------------
      if(itbl == 0) write(6,203) slab(i),e1,slab(j),e2
  203 format(/1x,a4,f7.2,14x,'---- no table found ----'
     .       /1x,a4,f7.2)
      if(itbl > 0) then
      call hyfnot(r1,rm1,lp1,nx1,lxi(1,i),s1,sm1,kp1,mx1,kx1,jxi1,
     .   cnot1,lnot)
      call hyfnot(r2,rm2,lp2,nx2,lxi(1,j),s2,sm2,kp2,mx2,kx2,jxi2,
     .   cnot2,lnot)
      write(6,200) slab(i),e1,itbl,cnot1,s1,kp1,ll1,(fx1(k),k=1,mx1)
      write(6,201) slab(j),e2,     cnot2,s2,kp2,ll2,(fx2(k),k=1,mx2)
      endif
      if(itbl < 0) then
      call hyfnot(r1,rm1,lp1,nx1,lxi(1,i),s2,sm2,kp2,mx2,kx2,jxi1,
     .   cnot1,lnot)
      call hyfnot(r2,rm2,lp2,nx2,lxi(1,j),s1,sm1,kp1,mx1,kx1,jxi2,
     .   cnot2,lnot)
      write(6,200) slab(i),e1,itbl,cnot1,s2,kp2,ll2,(fx2(k),k=1,mx2)
      write(6,201) slab(j),e2,     cnot2,s1,kp1,ll1,(fx1(k),k=1,mx1)
      endif
  200 format(/1x,a4,f7.2,i4,1x,a7,f7.3,i4,i6,5f6.1:
     .  /41x,5f6.1)
  201 format( 1x,a4,f7.2,4x,1x,a7,f7.3,i4,i6,5f6.1:
     .  /41x,5f6.1)
  209 format(/' hyfchk:')
  210 format(' spec   ephi  tbl  notes    rmt  lp   lxi   exi')
  211 format(' <= case ==>  ',
     .  '<============ data for selected table =================>')
  10  continue
      call hyfprn(lnot)
      return
      end
