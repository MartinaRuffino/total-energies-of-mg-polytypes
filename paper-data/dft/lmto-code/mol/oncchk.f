      subroutine oncchk(nspec,r,rsm,nphi,lphi,ephi,nxi,lxi,exi,
     .   slab,n0,t)
c  checks if tables are available for all needed 1c cases
      implicit real*8 (a-h,p-z), integer (o)
      dimension t(100,1),r(1),nphi(1),lphi(n0,1),ephi(n0,1),
     .   nxi(1),lxi(n0,1),exi(n0,1),jxi(10),fxi(10),
     .   lnot(6),kxi(10),rsm(1)
      character*8 slab(1)
      character*7 cnot
      character*15 ccc
      call getigv(1,ipr)
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
      if(lj == -1) lj=9
  8   lp=10*lp+lj
  3   write(6,449) slab(i),r(i),lp,ccc,lx,(exi(j,i),j=1,nxi(i))
  449 format(' ',a4,f7.3,i6,1x,a15,i7,6f6.1)
  448 format(/' spec    rmt  lphi   ephi',13x,'lxi   exi')
      write(6,*) ' '
      endif

c ------- start loop over cases for which tables are needed ----
      write(6,211)
      write(6,210)
      do 10 i=1,nspec
      do 10 ie=1,nphi(i)
      do 10 je=1,nphi(i)
      r1=r(i)
      rm1=rsm(i)
      nx1=nxi(i)
      e1=ephi(ie,i)
      e2=ephi(je,i)
      lp1=lphi(ie,i)
      lp2=lphi(je,i)
      if(lp1 < 0.or.lp2 < 0) goto 10
      call oncloc(r1,lp1,e1,lp2,e2,nx1,lxi(1,i),exi(1,i),t,itbl,jxi)

c ------- get table data ------------------------
      it=iabs(itbl)
      if(it > 0) then
      s1          =        t( 1,it)
      lsym        =idnint( t(13,it) )
      lmax        =idnint( t( 8,it) )
      kp1         =idnint( t(15,it) )
      fp1         =        t(16,it)
      kp2         =idnint( t(17,it) )
      fp2         =        t(18,it)
      mx1         =idnint( t(19,it) )
      sm1         =        t(21,it)
      ll1=0
      do 4 k=1,mx1
      kxi(k)=    idnint( t(40+k,it) )
      ll1=10*ll1+idnint( t(40+k,it) )
  4   fxi(k)      =      t(50+k,it)
      endif

c ------- printout ----------------------------
      if(itbl == 0) write(6,203) slab(i),e1,e2
  203 format(/1x,a4,f7.2,14x,'---- no table found ----'
     .       /5x,f7.2)
      if(itbl /= 0) then
      call oncnot(r1,rm1,lp1,lp2,nx1,lxi(1,i),s1,sm1,kp1,kp2,mx1,kxi,
     .   lmax,jxi,cnot,lnot)
      write(6,200) slab(i),e1,itbl,cnot,s1,kp1,ll1,(fxi(k),k=1,mx1)
      write(6,201) e2,kp2
      endif
  200 format(/1x,a4,f7.2,i4,1x,a7,f7.3,i4,i6,5f6.1:
     .  /41x,5f6.1)
  201 format(5x,f7.2,19x,i4)
  209 format(/' oncchk:')
  210 format(' spec   ephi  tbl  notes    rmt','  lp   lxi   exi')
  211 format(' <= case ==>  ',
     .  '<============ data for selected table ==============>')
  10  continue
      call oncprn(lnot)
      end
