      subroutine vesm2f(nxi,lxi,exi,alat,pos,rmt,n0,pot0,rhoi,
     .   nbas,ips,cg,jcg,indxcg,cy,ioff,zeta,f)
C- Add part from interstitial phi0 to zeta=integral phi * xi
      implicit real*8 (a-h,p-z), integer (o)
      dimension rmt(1),nxi(1),lxi(n0,1),exi(n0,1),pos(3,1),
     .  cg(1),jcg(1),indxcg(1),cy(1),pot0(1),ips(1),
     .  ioff(1),zeta(1),taui(3),tauj(3),rhoi(1),f(3,1),fj(3)
      real w(1)
      common /w/ w
      call getpr(ipr)
      call tcn('vesm2f')
      ndim=49
      call defrr(os,     ndim*ndim)
      call defrr(ogs,    ndim*ndim*3)
      call defrr(ob1,    ndim*ndim)
      call defrr(ogb1,   ndim*ndim*3)
      call defrr(ob2,    ndim*ndim)
      call defrr(ogb2,   ndim*ndim*3)
      e1=0d0
      i0=0
c ------ start loop over sites for functions with e=0 -- ---
      do 10 ib=1,nbas
      is=ips(ib)
      l1=0
      do 3 m=1,3
  3   taui(m)=alat*pos(m,ib)
      do 2 ie=1,nxi(is)
  2   l1=max(l1,lxi(ie,is))
      do 12 jb=1,nbas
      js=ips(jb)
      j0=ioff(jb)
      do 4 m=1,3
  4   tauj(m)=alat*pos(m,jb)
      do 12 je=1,nxi(js)
      r1=rmt(is)
      r2=rmt(js)
      l2=lxi(je,js)
      e2=exi(je,js)
      nlm1=(l1+1)**2
      nlm2=(l2+1)**2
      if(nlm1*nlm2 > ndim*ndim) call rx('vesum2: increase ndim')
c ------ case ib ne jb ---
      if (ib /= jb) then
        do 1 m=1,3
  1     fj(m)=0d0
        call han2cg(taui,tauj,e1,e2,r1,r2,nlm1,nlm2,w(ob1),w(ob2),
     .     w(ogb1),w(ogb2),cg,indxcg,jcg,cy,w(os),w(ogs))
        call vs0x2f(nlm1,nlm2,i0,j0,pot0,rhoi,w(os),zeta,w(ogs),fj)
        do 7 m=1,3
        f(m,ib)=f(m,ib)-fj(m)*0.5d0
  7     f(m,jb)=f(m,jb)+fj(m)*0.5d0

c ------ case ib eq jb ---
      else
        lmax=min(l1,l2)
        call han1ci(e1,e2,r1,lmax,w(os))
        call ves0x1(lmax,i0,j0,pot0,w(os),zeta)
      endif
   12 j0=j0+nlm2
   10 i0=i0+nlm1
      call rlse(os)
      if(ipr >= 60) then
      do 33 ib=1,nbas
  33  write(6,339) ib,(f(m,ib),m=1,3)
  339 format(i6,'  f=',3f13.6)
      endif

c ------ printout ---
      if (ipr < 80) return
      write(6,300)
      i=0
      do 30 ib=1,nbas
      is=ips(ib)
      do 30 ie=1,nxi(is)
      do 30 ilm=1,(lxi(ie,is)+1)**2
      i=i+1
      top=dabs(zeta(i))
      if(top >= 1d-5) write(6,301) ib,exi(ie,is),ilm,zeta(i)
  301  format(i4,f8.3,i5,f12.6,1x,f14.6)
  300  format(/' vesum2:'/'  ib    exi   ilm        zeta')
   30 continue
      call tcx('vesm2f')
      end

c -------- vs0x2f: add 2c interstitial contributions -----
      subroutine vs0x2f(nlm1,nlm2,i0,j0,pot0,rhoi,s,zeta,gs,fj)
      implicit real*8 (a-h,p-z), integer (o)
      dimension pot0(1),s(nlm1,1),zeta(1),rhoi(1),gs(nlm1,nlm2,3),
     .   fj(3)
      do 20 jlm=1,nlm2
      do 20 ilm=1,nlm1
      do 21 m=1,3
  21  fj(m)=fj(m)-pot0(i0+ilm)*gs(ilm,jlm,m)*rhoi(j0+jlm)
  20  zeta(j0+jlm)=zeta(j0+jlm)+pot0(i0+ilm)*s(ilm,jlm)
      end

c -------- ves0x1: add 1c interstitial contributions -----
      subroutine ves0x1(lmax,i0,j0,pot0,s0,zeta)
      implicit real*8 (a-h,p-z), integer (o)
      dimension pot0(1),s0(0:lmax),zeta(1)
      ilm=0
      do 20 l=0,lmax
      do 20 m=-l,l
      ilm=ilm+1
  20  zeta(j0+ilm)=zeta(j0+ilm)+pot0(i0+ilm)*s0(l)
      end
