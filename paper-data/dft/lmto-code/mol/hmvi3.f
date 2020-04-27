      subroutine hmvi3(el,nel,nphi,lphi,nxi,lxi,exi,pos,r,
     .  n0,nri,nvi,nbas,ips,alat,tspec,tdata,
     .  cg,jcg,indxcg,ioff,ioffv0,zetp,zets,zet0,zet0s,nhs,h)
C- Adds contribution from interstitial potential to hamiltonian matrix
C  using hyfget.
      implicit real*8 (a-h,p-z), integer (o)
      dimension r(1),nxi(1),lxi(n0,1),exi(n0,1),pos(3,1),dr(3),
     .  nphi(1),lphi(n0,1),cg(1),jcg(1),indxcg(1),zets(nri,1),
     .  tspec(1),tdata(1),zetp(nri),ips(1),h(nhs,nhs),el(1),ioff(1),
     .  zet0(nvi),zet0s(nvi,nbas),ioffv0(1)
      real w(1)
      common /w/ w
      call getpr(ipr)
      call tcn('hmvi3')
      nrot=81
      call defrr(ormat,    nrot*nrot)
c ------- start loop over matrix sub-blocks --------
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
      nx1=nxi(is)
      nx2=nxi(js)
      nlm1=(l1+1)**2
      nlm2=(l2+1)**2
      nf1=ioff(ib+1)-ioff(ib)
      nf2=ioff(jb+1)-ioff(jb)
      nf1g=nf1+(lxi(nx1+1,is)+1)**2
      nf2g=nf2+(lxi(nx2+1,js)+1)**2
      i1=ioff(ib)+1
      j1=ioff(jb)+1
C|    if (l1 < 0 .or. l2 < 0 .or. jham > iham) goto 12
      if (l1 < 0 .or. l2 < 0 ) goto 12
c      if(jham > iham) goto 12
c --- case ib ne jb ---
      if (ib /= jb) then
        call defrr(oc1,     nf1*nlm1*nlm2)
        call defrr(oc2,     nf2*nlm1*nlm2)
        call defrr(ot,      nlm1*nlm2)
        call defrr(ou,      nlm1*nlm2)
        call defrr(owr,     max0(nf1,nf2))
        call hyfevl(dr,r1,e1,nlm1,r2,e2,nlm2,w(oc1),w(oc2),
     .    w(ormat),nrot,nf1,nf2,nxi(is),lxi(1,is),exi(1,is),
     .    nxi(js),lxi(1,js),exi(1,js),tspec,tdata)
        call hmvi2(nlm1,nlm2,nf1,nf2,w(oc1),w(oc2),zetp(i1),zetp(j1),
     .    zets(i1,ib),zets(i1,jb),zets(j1,ib),zets(j1,jb),
     .    nhs,w(ormat),nrot,nxi(is),lxi(1,is),
     .    nxi(js),lxi(1,js),w(ot),w(ou),w(owr),h(iham,jham))

c --- case ib eq jb ---
      else
        call defrr(oc1,     nf1g*nlm1*nlm2)
        call oncget(r1,e1,nlm1,e2,nlm2,w(oc1),nf1g,nxi(is)+1,
     .    lxi(1,is),exi(1,is),tspec,tdata,cg,jcg,indxcg)
        i1g=ioffv0(ib)+1
        call hmvi1g(nlm1,nlm2,nf1,nf1g,w(oc1),zetp(i1),zets(i1,ib),
     .    zet0(i1g),zet0s(i1g,ib),nhs,h(iham,jham))
      endif
      call rlse(oc1)
   12 jham=jham+nlm2
   10 iham=iham+nlm1
      call rlse(ormat)

c --- copy to other triangle -------
      do 30 i=1,nhs
      do 30 j=1,i
  30  h(j,i)=h(i,j)

      call tcx('hmvi3')
      end

C --------- sub hmvi1g ---------------------------
      subroutine hmvi1g(nlm1,nlm2,nf1,nf1g,c1,zet1,zets,
     .   zet0,zet0s,nhs,h)
      implicit real*8 (a-h,p-z), integer (o)
      dimension h(nhs,nhs),c1(nf1g,nlm1,nlm2),zet1(1),zets(1),
     .   zet0(1),zet0s(1)
      do 10 i2=1,nlm2
      do 10 m=1,nf1
      do 10 i1=1,nlm1
   10 h(i1,i2)=h(i1,i2) + (zet1(m)-zets(m))*c1(m,i1,i2)
      do 11 i2=1,nlm2
      do 11 m=1,nf1g-nf1
      do 11 i1=1,nlm1
   11 h(i1,i2)=h(i1,i2) + (zet0(m)-zet0s(m))*c1(m+nf1,i1,i2)
      end
C --------- sub hmvi2 ----------------------------
      subroutine hmvi2(nlm1,nlm2,nf1,nf2,c1,c2,zet1,zet2,zets11,zets21,
     .  zets12,zets22,nhs,rmat,nrot,nx1,lx1,nx2,lx2,t,u,wr,h)
      implicit real*8 (a-h,p-z), integer (o)
      parameter (nlx=49)
      dimension h(nhs,nhs),c1(nf1,nlm1,nlm2),c2(nf2,nlm2,nlm1),
     .  zet1(1),zet2(1),lx1(1),lx2(1),t(nlm1,nlm2),u(nlm1,nlm2),
     .  wr(1),rmat(nrot,nrot),wk(nlx),
     .  zets11(1),zets12(1),zets21(1),zets22(1)
      call dpzero(t,   nlm1*nlm2)
C ... Subtract spherical contributions to zetl, rotate into wr
      i1=0
      do 12 ie=1,nx1
      nlm=(lx1(ie)+1)**2
      if (nlm > nlx) call rx('hmvi: nlm gt nlx')
      do  11  ilm = 1, nlm
   11 wk(ilm) = zet1(i1+ilm) - zets11(i1+ilm) - zets21(i1+ilm)
      call prerot(wk,wr(i1+1),nlm,rmat,nrot)
   12 i1=i1+nlm
C ... Add to unrotated hamiltonian
      do 10 i2=1,nlm2
      do 10 m=1,nf1
      do 10 i1=1,nlm1
   10 t(i1,i2)=t(i1,i2)+wr(m)*c1(m,i1,i2)
C ... Subtract spherical contributions to zet2, rotate into wr
      i2=0
      do 22 ie=1,nx2
      nlm=(lx2(ie)+1)**2
      if (nlm > nlx) call rx('hmvi: nlm gt nlx')
      do  21  ilm = 1, nlm
   21 wk(ilm) = zet2(i2+ilm) - zets12(i2+ilm) - zets22(i2+ilm)
      call prerot(wk,wr(i2+1),nlm,rmat,nrot)
   22 i2=i2+nlm
C ... Add to unrotated hamiltonian
      do 20 i2=1,nlm2
      do 20 m=1,nf2
      do 20 i1=1,nlm1
   20 t(i1,i2)=t(i1,i2)+wr(m)*c2(m,i2,i1)
C ... Apply 2 more rotations, add into hamiltonian matrix
      call dpzero(u,nlm1*nlm2)
      do 52 ilm2=1,nlm2
      do 52 jlm2=1,nlm2
      do 52 ilm1=1,nlm1
   52 u(ilm1,ilm2)=u(ilm1,ilm2)+rmat(jlm2,ilm2)*t(ilm1,jlm2)
      do 54 ilm2=1,nlm2
      do 54 jlm1=1,nlm1
      do 54 ilm1=1,nlm1
   54 h(ilm1,ilm2)=h(ilm1,ilm2)+rmat(jlm1,ilm1)*u(jlm1,ilm2)
      end
      subroutine hprt(h,ndim,n)
      implicit real*8 (a-h,p-z), integer (o)
      double precision h(ndim,n)

      do  10  i = 1, n
      do  10  j = 1, i
        top = dmax1(dabs(h(i,j)), dabs(h(j,i)))
        if (top > 1d-5) print 333, i,j, h(i,j), h(j,i), h(i,j)-h(j,i)
  333   format(2i5,3f12.6)
   10 continue
      end
