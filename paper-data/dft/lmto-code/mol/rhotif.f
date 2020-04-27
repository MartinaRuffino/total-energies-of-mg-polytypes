      subroutine rhotif(wgt,ewt,nstate,el,nel,nphi,lphi,nxi,lxi,exi,
     .  pos,r,n0,nbas,ips,alat,tspec,tdata,cg,jcg,indxcg,
     .  ioff,nhs,t,d,rhoi,zetp,f)
c  Adds together output intstl density and corresponding forces.
      parameter( nrot=81, nfx=300 )
      implicit real*8 (a-h,p-z), integer (o)
      dimension nxi(1),lxi(n0,1),exi(n0,1),pos(3,1),dr(3),r(1),
     .  nphi(1),lphi(n0,1),cg(1),jcg(1),indxcg(1),d(nhs,nhs),
     .  tspec(1),tdata(1),ips(1),t(nhs,nhs),el(1),ioff(1),rhoi(1),
     .  zetp(1),f(3,nbas),fp(3),ewt(1)
      real w(1)
      common /w/ w
      call getpr(ipr)
      call tcn('rhotif')
      istate = 1
      jstate = nstate

c --- make density matrix ----------
      call dpzero(d,   nhs*nhs)
      do 15 j=1,nhs
      do 15 nu=1,nstate
      wnu=wgt*ewt(nu)
      do 15 i=1,j
  15  d(i,j)=d(i,j)+ (wnu*t(j,nu)) * t(i,nu)
      do 16 j=1,nhs
      do 16 i=1,j
  16  d(j,i)=d(i,j)

c --- setup -----------------
      call defrr(ormat,   nrot*nrot)
      call defrr(owk,     nfx)
      call defrr(owr,     nfx)
      call defrr(ox,      nrot*nrot)
      call defrr(oy,      nrot*nrot)

c --- start loop over hamiltonian blocks -----
      iham=1
      do 10 ie=1,nel
      e1=el(ie)
      do 10 ib=1,nbas
      is=ips(ib)
      r1=r(is)
      jham=1
      do 12 je=1,nel
      e2=el(je)
      do 12 jb=1,nbas
      js=ips(jb)
      r2=r(js)
      do 14 i=1,3
   14 dr(i)=alat*(pos(i,jb)-pos(i,ib))
      l1=lphi(ie,is)
      l2=lphi(je,js)
      nlm1=(l1+1)**2
      nlm2=(l2+1)**2
      nf1=ioff(ib+1)-ioff(ib)
      nf2=ioff(jb+1)-ioff(jb)
      i1=ioff(ib)+1
      j1=ioff(jb)+1
      if (l1 < 0 .or. l2 < 0) goto 12
      if(jham > iham) goto 12
      if(nf1 > nfx.or.nf2 > nfx) call rx('rhotif increase nfx')
      wgt1=1d0
      if(iham /= jham) wgt1=2d0
      call defrr(oc1,     nf1*nlm1*nlm2)
      call defrr(oc2,     nf2*nlm1*nlm2)
      call defrr(og1,     nf1*nlm1*nlm2)
      call defrr(og2,     nf2*nlm1*nlm2)
c --- case ib ne jb ---
      if (ib /= jb) then
        call hyfevg(dr,r1,e1,nlm1,r2,e2,nlm2,w(oc1),w(oc2),
     .    w(og1),w(og2),w(ormat),nrot,nf1,nf2,nxi(is),lxi(1,is),
     .    exi(1,is),nxi(js),lxi(1,js),exi(1,js),tspec,tdata,itab)
c ... rotate subblock of density matrix
        call defrr(ob,    nlm1*nlm2)
        call defrr(owk2,  nlm1*nlm2)
        call prerom(nrot,w(ormat),nlm1,nlm2,nhs,d(iham,jham),
     .    w(owk2),w(ob))
c ... add to output density
        call xrhi42(nlm1,nlm2,w(ob),nf1,nf2,w(oc1),w(oc2),
     .    wgt1,w(owk),w(owr),rhoi(i1),
     .    rhoi(j1),w(ormat),nrot,nxi(is),lxi(1,is),nxi(js),lxi(1,js))
c  the next three calls make the force in the rotated system
        call yyhi42(nlm1,nlm2,nf1,nf2,w(og1),w(og2),w(ob),
     .    w(owk),zetp(i1),zetp(j1),w(ormat),
     .    nrot,nxi(is),lxi(1,is),nxi(js),lxi(1,js),w(ox),fz)
        call zzhi42(nlm1,nlm2,nf1,nf2,w(oc1),w(oc2),w(ob),
     .    zetp(i1),zetp(j1),w(ormat),nrot,nxi(is),lxi(1,is),
     .    nxi(js),lxi(1,js),w(ox),w(oy),w(owk),1,fx)
        call zzhi42(nlm1,nlm2,nf1,nf2,w(oc1),w(oc2),w(ob),
     .    zetp(i1),zetp(j1),w(ormat),nrot,nxi(is),lxi(1,is),
     .    nxi(js),lxi(1,js),w(ox),w(oy),w(owk),2,fy)
c --- add to forces ------------
      dd=dsqrt(dr(1)**2+dr(2)**2+dr(3)**2)
      fac=-wgt1
      if(itab < 0) fac=-fac
      fx=-fx*fac/dd
      fy=-fy*fac/dd
      fz=fz*fac
      call zzrotf(w(ormat),nrot,fx,fy,fz,fp(1),fp(2),fp(3))
      do 18 m=1,3
      f(m,ib)=f(m,ib)-fp(m)
  18  f(m,jb)=f(m,jb)+fp(m)

c --- case ib eq jb ---
      else
        call oncget(r1,e1,nlm1,e2,nlm2,w(oc1),nf1,nxi(is),
     .    lxi(1,is),exi(1,is),tspec,tdata,cg,jcg,indxcg)
        call xrhi41(nlm1,nlm2,nf1,w(oc1),nhs,d(iham,jham),
     .    wgt1,rhoi(i1))
      endif

      call rlse(oc1)
   12 jham=jham+nlm2
   10 iham=iham+nlm1
      call rlse(ormat)

c --- printout forces ----------
      if(ipr >= 30) then
      write(6,289)
      do 55 ib=1,nbas
  55  write(6,288) ib,(f(m,ib),m=1,3)
  288 format(i6,3f13.6)
  289 format(/'    ib     total eigenvalue force ')
      endif
      call tcx('rhotif')
      end
c -------- xrhi41: 1c contribution to density -----
      subroutine xrhi41(nlm1,nlm2,nf,c,nhs,d,wgt1,rhoi)
      implicit real*8 (a-h,p-z), integer (o)
      dimension c(nf,nlm1,nlm2),d(nhs,nlm2),rhoi(1)
      do 10 i1=1,nlm1
      do 10 i2=1,nlm2
      do 10 m=1,nf
  10  rhoi(m)=rhoi(m)+c(m,i1,i2)*(d(i1,i2)*wgt1)
      end
c -------- xrhi42: 2c contribution to density -----
      subroutine xrhi42(nlm1,nlm2,b,nf1,nf2,c1,c2,
     .  wgt1,wk,wr,rhoi1,rhoi2,rmat,nrot,nx1,lx1,nx2,lx2)
      implicit real*8 (a-h,p-z), integer (o)
      dimension c1(nf1,nlm1,nlm2),c2(nf2,nlm2,nlm1),b(nlm1,nlm2),
     .  rhoi1(1),rhoi2(1),wk(1),wr(1),lx1(1),lx2(1),rmat(nrot,nrot)
c  add together unrotated density on site 1
      call dpzero(wk,nf1)
      do 10 i1=1,nlm1
      do 10 i2=1,nlm2
      do 10 m=1,nf1
  10  wk(m)=wk(m)+c1(m,i1,i2)*b(i1,i2)
c  rotate density and add into rho1
      i1=1
      do 15 ie=1,nx1
      nlm=(lx1(ie)+1)**2
      call posrot(wk(i1),wr(i1),nlm,rmat,nrot)
  15  i1=i1+nlm
      do 11 m=1,nf1
  11  rhoi1(m)=rhoi1(m)+wr(m)*wgt1
c  add together unrotated density on site 2
      call dpzero(wk,nf2)
      do 20 i1=1,nlm1
      do 20 i2=1,nlm2
      do 20 m=1,nf2
  20  wk(m)=wk(m)+c2(m,i2,i1)*b(i1,i2)
c  rotate density and add into rhoi2
      i1=1
      do 16 ie=1,nx2
      nlm=(lx2(ie)+1)**2
      call posrot(wk(i1),wr(i1),nlm,rmat,nrot)
  16  i1=i1+nlm
      do 22 m=1,nf2
  22  rhoi2(m)=rhoi2(m)+wr(m)*wgt1
      end
c -------- yyhi42: force parallel to connecting vector -----
      subroutine yyhi42(nlm1,nlm2,nf1,nf2,g1,g2,b,
     .  u,zet1,zet2,rmat,nrot,nx1,lx1,nx2,lx2,x,sum)
      implicit real*8 (a-h,p-z), integer (o)
      dimension g1(nf1,nlm1,nlm2),g2(nf2,nlm2,nlm1),x(1),
     .  b(nlm1,nlm2),zet1(1),zet2(1),
     .  w(1),u(1),lx1(1),lx2(1),rmat(nrot,nrot)
c  rotate zet1, make help matrix x, sp with density matrix
      i1=1
      do 15 ie=1,nx1
      nlm=(lx1(ie)+1)**2
      call prerot(zet1(i1),u(i1),nlm,rmat,nrot)
  15  i1=i1+nlm
      call hhrif1(nlm1,nlm2,nf1,g1,u,x)
      call dpdot(x,b,nlm1*nlm2,sum1)
c  rotate zet2, make help matrix x, sp with density matrix
      i1=1
      do 16 ie=1,nx2
      nlm=(lx2(ie)+1)**2
      call prerot(zet2(i1),u(i1),nlm,rmat,nrot)
  16  i1=i1+nlm
      call hhrif2(nlm1,nlm2,nf2,g2,u,x)
      call dpdot(x,b,nlm1*nlm2,sum2)
      sum=sum1+sum2
      end
c -------- zzhi42: force orthogonal to connecting vector -----
      subroutine zzhi42(nlm1,nlm2,nf1,nf2,c1,c2,b,zet1,zet2,
     .   rmat,nrot,nx1,lx1,nx2,lx2,x,y,u,ixy,sum)
      implicit real*8 (a-h,p-z), integer (o)
      parameter( nfx=200 )
      dimension c1(nf1,nlm1,nlm2),c2(nf2,nlm2,nlm1),x(1),y(1),
     .  b(nlm1,nlm2),zet1(1),zet2(1),u(nlm1,nlm2),
     .  lx1(1),lx2(1),rmat(nrot,nrot),f(nfx),g(nfx)
      sum=0d0
c  make S*B + B*St where S = diff rot matrix
      call dpzero(u,   nlm1*nlm2)
      call premxy(b,u,nlm1,nlm2,ixy)
      call pretxy(b,u,nlm1,nlm2,ixy)
C|    call zzhi40(nlm1,nlm2,ixy,b,u)
c  rotate zet1
      i1=1
      do 15 ie=1,nx1
      nlm=(lx1(ie)+1)**2
      call prerot(zet1(i1),f(i1),nlm,rmat,nrot)
      call prerxy(f(i1),g(i1),nlm,ixy)
  15  i1=i1+nlm
c  make help matrices x,y for first site
      call hhrif1(nlm1,nlm2,nf1,c1,g,x)
      call hhrif1(nlm1,nlm2,nf1,c1,f,y)
c  contract with B and S*B + B*St
      call dpdot(x,b,nlm1*nlm2,sum1)
      call dpdot(y,u,nlm1*nlm2,sum2)
c  rotate zet2
      i1=1
      do 16 ie=1,nx2
      nlm=(lx2(ie)+1)**2
      call prerot(zet2(i1),f(i1),nlm,rmat,nrot)
      call prerxy(f(i1),g(i1),nlm,ixy)
  16  i1=i1+nlm
c  make help matrices x,y for second site
      call hhrif2(nlm1,nlm2,nf2,c2,g,x)
      call hhrif2(nlm1,nlm2,nf2,c2,f,y)
c  contract with B and S*B + B*St
      call dpdot(x,b,nlm1*nlm2,sam1)
      call dpdot(y,u,nlm1*nlm2,sam2)
      sum=sum1+sum2+sam1+sam2
      end
c -------- little help routines ----------------
      subroutine hhrif1(nlm1,nlm2,nf,c,w,x)
      implicit real*8 (a-h,p-z)
      dimension c(nf,nlm1,nlm2),w(nf),x(nlm1,nlm2)
      call dpzero(x,  nlm1*nlm2)
      do 26 i=1,nlm1
      do 26 j=1,nlm2
      do 26 m=1,nf
  26  x(i,j)=x(i,j)+c(m,i,j)*w(m)
         end
      subroutine hhrif2(nlm1,nlm2,nf,c,w,x)
      implicit real*8 (a-h,p-z)
      dimension c(nf,nlm2,nlm1),w(nf),x(nlm1,nlm2)
      call dpzero(x,  nlm1*nlm2)
      do 26 i=1,nlm1
      do 26 j=1,nlm2
      do 26 m=1,nf
  26  x(i,j)=x(i,j)+c(m,j,i)*w(m)
         end
c -------- zzrotf: rotate summed force for atom pair -----
      subroutine zzrotf(r,nrot,fx,fy,fz,gx,gy,gz)
      implicit real*8 (a-h,p-z), integer (o)
      dimension r(nrot,nrot)
      gy= r(2,2)*fy + r(3,2)*fz + r(4,2)*fx
      gz= r(2,3)*fy + r(3,3)*fz + r(4,3)*fx
      gx= r(2,4)*fy + r(3,4)*fz + r(4,4)*fx
      end
c -------- prerom: pre-rotate matrix -----------
      subroutine prerom(nrot,rmat,nlm1,nlm2,nhs,d,wk,b)
      implicit real*8 (a-h,p-z), integer (o)
      dimension rmat(nrot,nrot),d(nhs,nlm2),wk(nlm1,nlm2),b(nlm1,nlm2)
      call dpzero(wk,nlm1*nlm2)
      do 1 j=1,nlm2
      do 1 k=1,nlm1
      do 1 i=1,nlm1
  1   wk(i,j)=wk(i,j)+rmat(i,k)*d(k,j)
      call dpzero(b,nlm1*nlm2)
      do 2 j=1,nlm2
      do 2 k=1,nlm2
      do 2 i=1,nlm1
  2   b(i,j)=b(i,j)+wk(i,k)*rmat(j,k)
      end
