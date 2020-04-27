      subroutine oncnot(r,rm,lp1,lp2,nx,lx,s,sm,kp1,kp2,mx,kx,lmax,
     .   jxi,cnot,lnot)
c  sub to set note switches, make note string
      implicit real*8 (a-h,p-z), integer (o)
      dimension jxi(1),lx(1),kx(1),lnot(6),lnt(6)
      character*7 cnot
      fuzz=1.d-10
      ldif=0
      do 7 i=1,mx
      j=jxi(i)
  7   ldif=max0(ldif,lx(j)-kx(i))
      do 9 k=1,6
  9   lnt(k)=0
      if(r-s > fuzz)                lnt(1)=1
      if(r-s > 0.3d0)               lnt(2)=1
      if(lp1 < kp1.or.lp2 < kp2)   lnt(3)=1
c|    if(...      )                  lnt(4)=1
      if(mx < nx)                   lnt(5)=1
      if(ldif > 0)                  lnt(6)=1
      notes=0
      do 8 k=1,6
      if(lnt(k) == 1) lnot(k)=1
  8   if(lnt(k) == 1) notes=notes*10+k
      cnot='      '
      if(notes > 0) write(cnot,'(i6,'')'')') notes
      return
      end
c -------- sub oncprn ----------------------------
      subroutine oncprn(lnot)
c  print out notes
      implicit real*8 (a-h,p-z), integer (o)
      dimension lnot(6)
      call getigv(1,ipr)
      if(ipr < 10) return
      write(6,*) ' '
      if(lnot(1) == 1) write(6,*)
     . ' 1) sphere radius in table is smaller than',
     . ' radius used in calculation.'
      if(lnot(2) == 1) write(6,*)
     . ' 2) radius difference is large ---',
     . ' accuracy could be bad.'
      if(lnot(3) == 1) write(6,*)
     . ' 3) table has higher phi angular momenta --',
     . ' smaller table is possible.'
c|    if(lnot(4) == 1) write(6,*)
c|   . ' 4) bla bla bla..'
      if(lnot(5) == 1) write(6,*)
     . ' 5) some xi-localizations not used,',
     . ' more accurate table is possible.'
      if(lnot(6) == 1) write(6,*)
     . ' 6) some xi-angular-momenta not used --',
     . ' more accurate table is possible.'
      return
      end