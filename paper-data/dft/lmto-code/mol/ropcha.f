      subroutine ropcha(spid,nspec,n0,nxi,lxi,exi,rsm,rcut,rint,radd,
     .   nc0,nc1,c1,nc2,c2,atr,tol1,tol2)
c  Sets up all chebyshev fits to the smooth hankel functions.
      implicit real*8 (a-h,p-z), integer (o)
      dimension nxi(nspec),lxi(n0,nspec),exi(n0,nspec),rsm(nspec),
     .   rcut(nspec),rint(nspec),nc2(nspec),c2(nc0,2,n0,nspec),
     .   c1(nc0,n0,n0,nspec),nc1(nspec)
      character*4 spid(1)
      call getpr(ipr)
      kpr=0
      if(ipr >= 70) kpr=3
      if(ipr >= 20) write(6,973) nc0,atr,radd,tol1,tol2
  973 format(/' ropcha:  n0=',i3,'   a=',f6.3,'   radd=',f6.3,
     .  '   tol=',1p,2d9.1)
c ------ start loop over classes ----------
      if(ipr >= 30) write(6,752)
      do 10 is=1,nspec
      rs=rsm(is)
      r0=0d0
      r1=rcut(is)
      r2=rint(is)+radd
      ne=nxi(is)
      nc1(is)=2
      nc2(is)=2
      n1mn=nc0
      n2mn=nc0
c ------ start loop over energies ---------
      do 11 ie=1,ne
      lmax=lxi(ie,is)
      lmux=1
      e=exi(ie,is)
      call ropchh(e,rs,lmax,r0,r1,atr,nc0,c1(1,1,ie,is),tol1,n1,kpr)
      call ropchh(e,rs,lmux,r1,r2,atr,nc0,c2(1,1,ie,is),tol2,n2,kpr)
      if(ipr >= 50) write(6,951) e,n1,n2
  951 format('     exi=',f11.4,'    n1,n2=',2i6)
      n1mn=min0(n1mn,n1)
      n2mn=min0(n2mn,n2)
      nc1(is)=max0(n1,nc1(is))
  11  nc2(is)=max0(n2,nc2(is))
      if(ipr >= 30) write(6,751) spid(is),rint(is),rcut(is),
     .   n1mn,nc1(is),n2mn,nc2(is)
  751 format(1x,a4,2f10.4,2x,2i6,2x,2i6)
  752 format(' spec',6x,'Rint',6x,'Rcut     n1min  n1     n2min  n2')
  10  if(max0(nc1(is),nc2(is)) >= nc0) call rx('ropcha: n gt nc0')
      end
