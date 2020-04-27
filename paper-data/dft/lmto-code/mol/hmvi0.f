      subroutine hmvi0(el,nel,nphi,lphi,nxi,lxi,exi,n0,r,ips,nbas,isp,
     .  tspec,tdata,cg,jcg,indxcg,ioff,nhtab,iax,rtab,ioffb,zetp,nhs,h)
C- Hamiltonian matrix from interstitial potential, q=0 (real H)
      implicit real*8 (a-h,p-z), integer (o)
      dimension r(1),nxi(1),lxi(n0,1),exi(n0,1),
     .  nphi(1),lphi(n0,1),cg(1),jcg(1),indxcg(1),
     .  tspec(1),tdata(1),zetp(1),ips(1),h(nhs,nhs),el(1),ioff(1),
     .  ioffb(nbas,nel),iax(10,1),rtab(3,1)
      real w(1)
      common /w/ w

C     call dpzero(h,nhs*nhs)
      call getpr(ipr)
      call tcn('hmvi0')
      nrot=81
      nsp=lsp()+1
      nri=ioff(nbas+1)
      call defrr(ormat,    nrot*nrot)
c ------- start loop over matrix sub-blocks --------
      do 10 it=1,nhtab
      ib=iax(1,it)
      jb=iax(2,it)
      is=ips(ib)
      js=ips(jb)
      do 10 ie=1,nel
      do 10 je=1,nel
      e1=el(ie)
      e2=el(je)
      r1=r(is)
      r2=r(js)
      l1=lphi(ie,is)
      l2=lphi(je,js)
      nlm1=(l1+1)**2
      nlm2=(l2+1)**2
      nf1=ioff(ib+1)-ioff(ib)
      nf2=ioff(jb+1)-ioff(jb)
      i1=1+ioff(ib)+nri*(isp-1)
      j1=1+ioff(jb)+nri*(isp-1)
      iham=(ioffb(ib,ie)+1)
      jham=(ioffb(jb,je)+1)
      if (l1 < 0 .or. l2 < 0 ) goto 10
      if(jham > iham) goto 10
c --- case dr nonzero ---
      if(dabs(rtab(1,it))+dabs(rtab(2,it))+dabs(rtab(3,it)) > 1d-8)then
        call defrr(oc1,     nf1*nlm1*nlm2)
        call defrr(oc2,     nf2*nlm1*nlm2)
        call defrr(ot,      nlm1*nlm2)
        call defrr(ou,      nlm1*nlm2)
        call defrr(owr,     max0(nf1,nf2))
        call hyfevl(rtab(1,it),r1,e1,nlm1,r2,e2,nlm2,w(oc1),w(oc2),
     .    w(ormat),nrot,nf1,nf2,nxi(is),lxi(1,is),exi(1,is),
     .    nxi(js),lxi(1,js),exi(1,js),tspec,tdata)
        call xhmv42(nlm1,nlm2,nf1,nf2,w(oc1),w(oc2),zetp(i1),
     .    zetp(j1),nhs,h(iham,jham),w(ormat),nrot,nxi(is),lxi(1,is),
     .    nxi(js),lxi(1,js),w(ot),w(ou),w(owr))
c --- case dr=0 ---
      else
        call defrr(oc1,     nf1*nlm1*nlm2)
        call oncget(r1,e1,nlm1,e2,nlm2,w(oc1),nf1,nxi(is),
     .    lxi(1,is),exi(1,is),tspec,tdata,cg,jcg,indxcg)
        call xhmv41(nlm1,nlm2,nf1,w(oc1),zetp(i1),nhs,h(iham,jham))
      endif
      call rlse(oc1)
   10 continue
      call rlse(ormat)

c --- copy to other triangle -------
      do 30 i=1,nhs
      do 30 j=1,i
  30  h(j,i)=h(i,j)

C     call prmx('h in hmvi0',h,nhs,nhs,nhs)
      call tcx('hmvi0')
      end
C --------- sub xhmv41 ----------------------------
      subroutine xhmv41(nlm1,nlm2,nf1,c1,zet1,nhs,h)
      implicit real*8 (a-h,p-z), integer (o)
      dimension h(nhs,nhs),c1(nf1,nlm1,nlm2),zet1(1)
      do 10 i2=1,nlm2
      do 10 i1=1,nlm1
      do 10 m=1,nf1
   10 h(i1,i2)=h(i1,i2)+zet1(m)*c1(m,i1,i2)
      end
C --------- sub xhmv42 ----------------------------
      subroutine xhmv42(nlm1,nlm2,nf1,nf2,c1,c2,zet1,zet2,
     .   nhs,h,rmat,nrot,nx1,lx1,nx2,lx2,t,u,wr)
      implicit real*8 (a-h,p-z), integer (o)
      dimension h(nhs,nhs),c1(nf1,nlm1,nlm2),c2(nf2,nlm2,nlm1),
     .  zet1(1),zet2(1),lx1(1),lx2(1),t(nlm1,nlm2),u(nlm1,nlm2),
     .  wr(1),rmat(nrot,nrot)
      call dpzero(t,   nlm1*nlm2)
c  rotate zet1 into wr
      i1=1
      do 12 ie=1,nx1
      nlm=(lx1(ie)+1)**2
      call prerot(zet1(i1),wr(i1),nlm,rmat,nrot)
C     call dgemm('n','n',nlm,1,nlm,1d0,rmat,nrot,zet1(i1),nlm,
C    .  0d0,wr(i1),nlm)
  12  i1=i1+nlm
c  add to unrotated hamiltonian
      do 10 i2=1,nlm2
      do 10 i1=1,nlm1
      do 10 m=1,nf1
   10 t(i1,i2)=t(i1,i2)+wr(m)*c1(m,i1,i2)
c  rotate zet2 into wr
      i1=1
      do 22 ie=1,nx2
      nlm=(lx2(ie)+1)**2
      call prerot(zet2(i1),wr(i1),nlm,rmat,nrot)
C     call dgemm('n','n',nlm,1,nlm,1d0,rmat,nrot,zet2(i1),nlm,
C    .  0d0,wr(i1),nlm)
  22  i1=i1+nlm
c  add to unrotated hamiltonian
      do 20 i2=1,nlm2
      do 20 i1=1,nlm1
      do 20 m=1,nf2
   20 t(i1,i2)=t(i1,i2)+wr(m)*c2(m,i2,i1)
c  apply 2 more rotations, add into hamiltonian matrix
C     call dpzero(u,nlm1*nlm2)
C     do 52 ilm2=1,nlm2
C     do 52 jlm2=1,nlm2
C     do 52 ilm1=1,nlm1
C 52  u(ilm1,ilm2)=u(ilm1,ilm2)+t(ilm1,jlm2)*rmat(jlm2,ilm2)
      call dgemm('N','N',nlm1,nlm2,nlm2,1d0,t,nlm1,rmat,nrot,
     .           0d0,u,nlm1)
C     do 54 ilm2=1,nlm2
C     do 54 ilm1=1,nlm1
C     do 54 jlm1=1,nlm1
C 54  h(ilm1,ilm2)=h(ilm1,ilm2)+rmat(jlm1,ilm1)*u(jlm1,ilm2)
      call dgemm('T','N',nlm1,nlm2,nlm1,1d0,rmat,nrot,u,nlm1,
     .           1d0,h,nhs)
      end
