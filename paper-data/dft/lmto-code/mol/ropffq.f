      subroutine ropffq(nbas,ips,n0,nxi,lxi,exi,nri,rsm,alat,pos,fac,
     .  nd,n,G,rhoi,rho0,rhoq,lfrc,grhoq,grhoi,job,jobg)
C-  Tabulation of Fourier transforms of smooth H's and G's.
C   Job=1: input is rhoi,rho0, output is FT of total rho in rhoq.
C          For this case, fac should be (cell vol)^-1
C   Job=2: input is FT of vxc in rhoq, output is scalar product of rhoq
C          with transforms of smooth H's and G's in rhoi,rho0.
C          For this case, fac should be 1
C   Jobg=0:no gaussians
C   Jobg=1:one gaussian with rsm same as exi
      implicit real*8 (a-h,p-z), integer (o)
      logical lfrc
      dimension ips(1),lxi(n0,1),exi(n0,1),nxi(1),rhoi(1),rho0(1),
     .  G(nd,3),rsm(1),pos(3,1),grhoi(1)
      complex*16 rhoq(n+1),grhoq(nd+1,3)
      real w(1)
      common /w/ w
      if(job /= 1.and.job /= 2) call rx('ropffq: job not 1 or 2')
c ------ get maximal angular momentum -------------
      ltop=-1
      do 1 ib=1,nbas
      is=ips(ib)
      do 1 ie=1,nxi(is)
  1   ltop=max0(ltop,lxi(ie,is))
      nlmtop=(ltop+1)**2
c ------ combine with FT's of H's and G's ---------
      call defrr(oyl,     n*nlmtop)
      call defrr(oqsq,    n)
      call ropyln(n,G,G(1,2),G(1,3),ltop,n,w(oyl),w(oqsq))

      call defcc(oh1,     n)
      call defcc(oh2,     n)
      rhoq(1) = (0d0,0d0)
      call hhrpfq(nbas,ips,n0,nxi,lxi,exi,nri,rsm,alat,pos,fac,n,
     .  G,G(1,2),G(1,3),w(oqsq),w(oyl),w(oh1),w(oh2),rhoi,rho0,
     .  rhoq(2),lfrc,grhoq(2,1),grhoi,nd+1,job,jobg)
      call rlse(oyl)
      end
c -------- sub hhrpfq -----------
      subroutine hhrpfq(nbas,ips,n0,nxi,lxi,exi,nri,rsm,alat,pos,fac,
     .  n,Gx,Gy,Gz,qsq,yl,h1,h2,rhoi,rho0,rhoq,lfrc,grhoq,grhoi,
     .  ndg,job,jobg)
      implicit none
      logical lfrc
      integer i,ib,ie,ilm,iri,ivi,is,job,l,lmax,lx,m,n,n0,nbas,nri,jobg
      double precision alat,cof,e,fac,gam,pi,sp,sum,tpi
      integer ips(1),lxi(n0,1),nxi(1),ip2,ndg,ix,nlm
      double precision exi(n0,1),rhoi(1),rho0(1),
     .  rsm(1),pos(3,1),Gx(n),Gy(n),Gz(n),qsq(n),yl(n,1),
     .  grhoi(nri,3)
      complex*16 rhoq(n),grhoq(ndg,3),h1(n),h2(n),sm1
      pi=4d0*datan(1d0)
      tpi=2d0*pi
      if(job == 1) call dpzero(rhoq,   2*n)

c ----- start loop over atoms ---------------
      iri=0
      ivi=0
      do 10 ib=1,nbas
      is=ips(ib)
      gam=0.25d0*rsm(is)**2
      do 14 i=1,n
      sp=alat*(Gx(i)*pos(1,ib)+Gy(i)*pos(2,ib)+Gz(i)*pos(3,ib))
   14 h1(i)=cdexp( dcmplx(-gam*qsq(i),-sp) )
      lmax=-1
c ----- loop over energies; ft of smooth hankels ----
      do 12 ie=1,nxi(is)
      e=exi(ie,is)
      lx=lxi(ie,is)
      nlm=(lx+1)**2
      lmax=max0(lmax,lx)
      do 16 i=1,n
   16 h2(i)=h1(i)/(e-qsq(i))
      sm1 = (0d0,1d0)
      ip2 = -1
      ilm = 0
      do 12 l = 0, lx
      sm1 = sm1*(0d0,-1d0)
      ip2 = mod(ip2+1,4)
      do 12 m = -l, l
      ilm = ilm+1
      iri=iri+1
      if(job == 1) then
        cof = -4d0*pi*dexp(gam*e)*rhoi(iri)*fac
        goto (40,41,42,43) ip2+1
   40   do 46 i = 1, n
   46   rhoq(i)=rhoq(i) + (cof*yl(i,ilm)) * h2(i)
        goto 44
   41   do 47 i = 1, n
   47   rhoq(i)=rhoq(i) + (cof*yl(i,ilm)) *
     .      dcmplx(dimag(h2(i)),-dble(h2(i)))
        goto 44
   42   do 48 i = 1, n
   48   rhoq(i)=rhoq(i) - (cof*yl(i,ilm)) * h2(i)
        goto 44
   43   do 49 i = 1, n
   49   rhoq(i)=rhoq(i) - (cof*yl(i,ilm)) *
     .      dcmplx(dimag(h2(i)),-dble(h2(i)))
   44   continue
C       do 19 i=1,n
C  19   rhoq(i)=rhoq(i) + (cof*yl(i,ilm)) * h2(i)*sm1
      else
        sum = 0d0
C       do 29 i=1,n
C  29   sum = sum +  yl(i,ilm) * dble(dconjg(rhoq(i))*h2(i)*sm1)
        goto (50,51,52,53) ip2+1
   50   do 56 i = 1, n
   56   sum = sum +  yl(i,ilm) * dble(dconjg(rhoq(i))*h2(i))
        goto 54
   51   do 57 i = 1, n
   57   sum = sum +  yl(i,ilm) *
     .      dble(dconjg(rhoq(i))*dcmplx(dimag(h2(i)),-dble(h2(i))))
        goto 54
   52   do 58 i = 1, n
   58   sum = sum -  yl(i,ilm) * dble(dconjg(rhoq(i))*h2(i))
        goto 54
   53   do 59 i = 1, n
   59   sum = sum -  yl(i,ilm) *
     .      dble(dconjg(rhoq(i))*dcmplx(dimag(h2(i)),-dble(h2(i))))
   54   continue
        rhoi(iri) = rhoi(iri) - 4d0*pi*dexp(gam*e)*sum*fac
        if (lfrc) then
          do  120  ix = 1, 3
          sum = 0d0
C         do  129  i = 1, n
C 129     sum = sum +  yl(i,ilm) * dble(dconjg(grhoq(i,ix))*h2(i)*sm1)
          goto (150,151,152,153) ip2+1
  150     do 156 i = 1, n
  156     sum = sum +  yl(i,ilm) * dble(dconjg(grhoq(i,ix))*h2(i))
          goto 154
  151     do 157 i = 1, n
  157     sum = sum +  yl(i,ilm) *
     .       dble(dconjg(grhoq(i,ix))*dcmplx(dimag(h2(i)),-dble(h2(i))))
          goto 154
  152     do 158 i = 1, n
  158     sum = sum -  yl(i,ilm) * dble(dconjg(grhoq(i,ix))*h2(i))
          goto 154
  153     do 159 i = 1, n
  159     sum = sum -  yl(i,ilm) *
     .       dble(dconjg(grhoq(i,ix))*dcmplx(dimag(h2(i)),-dble(h2(i))))
  154     continue
          grhoi(iri,ix)=grhoi(iri,ix)-4d0*pi*dexp(gam*e)*sum*fac
  120     continue
        endif
      endif
   12 continue

c ----- ft of gaussians -------
      if (jobg == 0) goto 22
      if (lfrc) call rx('hhrpfq not set up for gaussians and forces')
      sm1 = (0d0,1d0)
      ip2 = -1
      ilm = 0
      do 20 l = 0, lx
      ip2 = mod(ip2+1,4)
      sm1 = sm1*(0d0,-1d0)
      do 20 m = -l, l
      ilm = ilm+1
      ivi=ivi+1
      if(job == 1) then
        cof = rho0(ivi)*fac
C       do 21 i=1,n
C  21   rhoq(i)=rhoq(i) + (cof*yl(i,ilm)) * h1(i)*sm1
        goto (60,61,62,63) ip2+1
   60   do 66 i = 1, n
   66   rhoq(i)=rhoq(i) + (cof*yl(i,ilm)) * h1(i)
        goto 64
   61   do 67 i = 1, n
   67   rhoq(i)=rhoq(i) + (cof*yl(i,ilm)) *
     .      dcmplx(dimag(h1(i)),-dble(h1(i)))
        goto 64
   62   do 68 i = 1, n
   68   rhoq(i)=rhoq(i) - (cof*yl(i,ilm)) * h1(i)
        goto 64
   63   do 69 i = 1, n
   69   rhoq(i)=rhoq(i) - (cof*yl(i,ilm)) *
     .      dcmplx(dimag(h1(i)),-dble(h1(i)))
   64   continue
      else
        sum = 0d0
c       do 24 i=1,n
c  24   sum = sum +  yl(i,ilm) * dble(dconjg(rhoq(i))*h1(i)*sm1)
        goto (70,71,72,73) ip2+1
   70   do 76 i = 1, n
   76   sum = sum +  yl(i,ilm) * dble(dconjg(rhoq(i))*h1(i))
        goto 74
   71   do 77 i = 1, n
   77   sum = sum +  yl(i,ilm) *
     .      dble(dconjg(rhoq(i))*dcmplx(dimag(h1(i)),-dble(h1(i))))
        goto 74
   72   do 78 i = 1, n
   78   sum = sum -  yl(i,ilm) * dble(dconjg(rhoq(i))*h1(i))
        goto 74
   73   do 79 i = 1, n
   79   sum = sum -  yl(i,ilm) *
     .      dble(dconjg(rhoq(i))*dcmplx(dimag(h1(i)),-dble(h1(i))))
   74   continue
        rho0(ivi) = rho0(ivi) + sum*fac
      endif
   20 continue
   22 continue

   10 continue
      end
