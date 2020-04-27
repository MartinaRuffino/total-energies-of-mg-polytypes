      subroutine vesm1f(nxi,lxi,exi,alat,pos,rmt,n0,poti,rhoi,
     .  nbas,ips,cg,jcg,indxcg,cy,ioff,zeta,f)
c  Make zeta, which is the integral of smooth-ves * xi_m
c  over interstitial. Here: add part from hankels with e /= 0.
      implicit real*8 (a-h,p-z), integer (o)
      dimension rmt(1),nxi(1),lxi(n0,1),exi(n0,1),pos(3,1),
     .  cg(1),jcg(1),indxcg(1),cy(1),ips(1),taui(3),tauj(3),
     .  ioff(1),zeta(1),poti(1),f(3,1),rhoi(1),fi(3),fj(3)
      real w(1)
      common /w/ w
      call getpr(ipr)
      call tcn('vesm1f')
      ndim=49
      call defrr(os,     ndim*ndim)
      call defrr(ogs,    ndim*ndim*3)
      call defrr(ob1,    ndim*ndim)
      call defrr(ogb1,   ndim*ndim*3)
      call defrr(ob2,    ndim*ndim)
      call defrr(ogb2,   ndim*ndim*3)
c --- start loop over sites ---
      do 10 ib=1,nbas
      is=ips(ib)
      r1=rmt(is)
      do 3 m=1,3
  3   taui(m)=alat*pos(m,ib)
      do 10 jb=1,ib
      js=ips(jb)
      r2=rmt(js)
      do 4 m=1,3
  4   tauj(m)=alat*pos(m,jb)
      i0=ioff(ib)
      do 11 ie=1,nxi(is)
      l1=lxi(ie,is)
      e1=exi(ie,is)
      nlm1=(l1+1)**2
      j0=ioff(jb)
      do 12 je=1,nxi(js)
      l2=lxi(je,js)
      e2=exi(je,js)
      nlm2=(l2+1)**2
      if(nlm1*nlm2 > ndim*ndim) call rx('vesum1: increase ndim')
c --- case ib ne jb ---
      if (ib /= jb) then
        do 1 m=1,3
        fi(m)=0d0
  1     fj(m)=0d0
        call han2cg(taui,tauj,e1,e2,r1,r2,nlm1,nlm2,w(ob1),w(ob2),
     .     w(ogb1),w(ogb2),cg,indxcg,jcg,cy,w(os),w(ogs))
        call vesx1f(nlm1,nlm2,i0,j0,poti,rhoi,w(os),zeta,w(ogs),fi,fj)
       do 2 m=1,3
       f(m,ib)=f(m,ib)+(fi(m)-fj(m))*0.5d0
  2    f(m,jb)=f(m,jb)+(fj(m)-fi(m))*0.5d0

c --- case ib eq jb ---
      else
        lmax=min(l1,l2)
        call han1ci(e1,e2,r1,lmax,w(os))
        call vesix1(lmax,i0,j0,poti,w(os),zeta)
      endif
   12 j0=j0+nlm2
   11 i0=i0+nlm1
   10 continue
      call rlse(os)
      if(ipr >= 60) then
      do 33 ib=1,nbas
  33  write(6,339) ib,(f(m,ib),m=1,3)
  339 format(i6,'  f=',3f13.6)
      endif

c --- printout ---
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
  301  format(i4,f8.3,i5,f12.6,1x,f12.6)
  300  format(/' vesum1:'/'  ib    exi   ilm      zeta')
   30 continue
      call tcx('vesm1f')
      end

c ----- vesx1f: add 2c interstitial contributions ----
      subroutine vesx1f(nlm1,nlm2,i0,j0,poti,rhoi,s,zeta,gs,fi,fj)
      implicit real*8 (a-h,p-z), integer (o)
      dimension poti(1),s(nlm1,1),zeta(1),gs(nlm1,nlm2,3),
     .   fi(3),fj(3),rhoi(1)
      do 10 jlm=1,nlm2
      do 10 ilm=1,nlm1
      do 11 m=1,3
   11 fi(m)=fi(m)+poti(j0+jlm)*gs(ilm,jlm,m)*rhoi(i0+ilm)
   10 zeta(i0+ilm)=zeta(i0+ilm)+poti(j0+jlm)*s(ilm,jlm)
      do 20 ilm=1,nlm1
      do 20 jlm=1,nlm2
      do 21 m=1,3
   21 fj(m)=fj(m)-poti(i0+ilm)*gs(ilm,jlm,m)*rhoi(j0+jlm)
   20 zeta(j0+jlm)=zeta(j0+jlm)+poti(i0+ilm)*s(ilm,jlm)
      end

c ----- vesix1: add 1c interstitial contributions -----
      subroutine vesix1(lmax,i0,j0,poti,s0,zeta)
      implicit real*8 (a-h,p-z), integer (o)
      dimension poti(1),s0(0:lmax),zeta(1)
      ilm=0
      do 20 l=0,lmax
      do 20 m=-l,l
      ilm=ilm+1
  20  zeta(i0+ilm)=zeta(i0+ilm)+poti(j0+ilm)*s0(l)
      end
