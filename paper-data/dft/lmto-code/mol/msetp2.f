      subroutine msetp2(spid,z,rmt,rsm,rint,rcut,nphi,lphi,ephi,nxi,
     .  lxi,exi,nsp,lmxl,lmxa,pin,p,piz,pz,lsc,
     .  nspec,n0,nel,el,atid,pos0,pol,pos,ips,nbas,nbasp,lmaxp,amp,
     .  ioff,ioffv0,ioffp,nri,nvi,ndim,nhs,nla,nll,qval,qsmc)
c  all kinds of setup for m-programs (two-panel)
      implicit real*8 (a-h,p-z), integer (o)
      dimension rmt(1),rsm(1),z(1),nxi(1),lxi(n0,1),exi(n0,1),
C     .  pin(n0,nsp,1),p(n0,nsp,1),pz(n0,nsp,1),piz(n0,nsp,1),
     .  pin(n0,2,1),p(n0,2,1),pz(n0,2,1),piz(n0,2,1),
     .  nphi(1),ephi(n0,1),lphi(n0,1),lmxa(1),lmxl(1),rint(1),
     .  ioffv0(1),pos(3,1),pol(3,1),ips(1),ioff(1),el(nel),
     .  ndim(nel),konf(10),ioffp(1),pos0(3,1),rcut(1)
      character*8 spid(1),atid(1)
      character*8 t,tt
      character*1 dig(-1:9)
      data dig /'-','0','1','2','3','4','5','6','7','8','9'/
      call getpr(ipr)
      if(ipr >= 10) write(6,330) nbas,nbasp,nspec,amp
  330 format(/' msetp2:   nbas=',i3,'  nbasp=',i3,'   nspec=',i3,
     .        '    amp=',f10.6)

c ------ put lmto energies into ephi ----------
      do 7 is=1,nspec
      if(nphi(is) > nel)
     .  call rx(' more phis than in elmto line,  sp='//spid(is))
      do 8 ie=1,nel
    8 ephi(ie,is)=el(ie)
      do 9 ie=nphi(is)+1,nel
    9 lphi(ie,is)=-1
    7 continue
CL      write(71,766) (el(ie),ie=1,nel)
  766 format(' mce  ephi',6f9.3)

c ------ add energy "zero" to exi lists ----------
      do 18 is=1,nspec
      nx=nxi(is)
      exi(nx+1,is)=0d0
      lmaxv=-1
      do 19 ie=1,nx
  19  lmaxv=max0(lmaxv,lxi(ie,is))
  18  lxi(nx+1,is)=lmaxv

c ------ print out species data ----------------
      if(ipr >= 1) write(6,261)
  261 format(/' spec    Z    conf    rmt',
     .    '     atoms   lmto-basis')
      do 26 i=1,nspec
      n=0
      do 27 j=1,nbas
   27 if(atid(j) == spid(i)) n=n+1
      write(t,'(8a1)') ' ','<',(dig(lphi(j,i)),j=1,nphi(i)),'>'
      nl=lmxa(i)+1
      do 23 k=1,nl
   23 konf(k)=pin(k,1,i)
      write(tt,'(i3,6i1)') (konf(k),k=1,nl)
      iz=z(i)+.1d0
CL      write(71,263) spid(i),n,rmt(i),rsm(i),rint(i),rcut(i),lmxa(i),t
   26 if(ipr >= 10) write(6,260) spid(i),iz,tt,
     .  rmt(i),n,t,(ephi(j,i),j=1,nphi(i))
  260 format(1x,a4,i5,2x,a,f9.4,i5,4x,a,5f7.3)
  263 format(' mce <',a4,'>  na',i3,'   r ',4f7.3,'   lxa',i2,'   bs',a)
      if(ipr >= 10) write(6,271)
  271 format(/' spec    rint    rsm     lmxl',
     .  '    interstitial chden basis')
      do 28 i=1,nspec
      if (nxi(i) == 7)
     .write(t,'(8a1)') '<',(dig(lxi(j,i)),j=1,nxi(i))
      if (nxi(i) < 7)
     .write(t,'(8a1)') '<',(dig(lxi(j,i)),j=1,nxi(i)),'>'
CL      write(71,273) spid(i),lmxl(i),t,(exi(j,i),j=1,nxi(i))
   28 if(ipr >= 10) write(6,270) spid(i),rint(i),rsm(i),lmxl(i),t,
     .  (exi(j,i),j=1,nxi(i))
  270 format(1x,a4,1x,f7.2,f8.3,i6,6x,a,6f6.1)
  273 format(' mce  <',a4,'>  lxl',i2,'   cbs ',a,'  exi',6f7.2)

c ------ atom pointers, list of atoms ------------------
      if(ipr >= 20) write(6,401)
      do 10 ib=1,nbas
      j=0
      do 11 is=1,nspec
  11  if(spid(is) == atid(ib)) j=is
      if(j == 0) call rx(' species <'//atid(ib)//'> not defined')
      ips(ib)=j
      do 12 m=1,3
  12  pos(m,ib)=pos0(m,ib)
C  12  pos(m,ib)=pos0(m,ib)+amp*pol(m,ib)
      if(ipr >= 20) write(6,400) ib,atid(ib),(pos(m,ib),m=1,3)
CL      write(71,404) ib,atid(ib),(pos(m,ib),m=1,3)
  400 format(i4,3x,a4,3f9.5)
  404 format(' mce  ',i4,'  <',a4,'>  ',3f9.5)
  401 format(/'  ib   spec',16x,'position')
  10  continue

c ------ offsets -------
      npjj=0
      do  42 ib = 1, nbas
      is=ips(ib)
      ioffp(ib)=npjj
      npjj = npjj + (lmxa(is)+1)**4
   42 continue
      ioffp(nbas+1)=npjj
c ... ioff, ioffv0
      nri=0
      nvi=0
      do 20 ib=1,nbas
      ioff(ib)=nri
      is=ips(ib)
      lmax=0
      do 22 j=1,nxi(is)
      lmax=max(lmax,lxi(j,is))
   22 nri=nri+(lxi(j,is)+1)**2
      if (lmax < lmxa(is)) call rx('msetp2: lxi lt lmxa')
      ioffv0(ib)=nvi
      nvi=nvi+(lmax+1)**2
   20 continue
      ioff(nbas+1)=nri
      ioffv0(nbas+1)=nvi
      if (nbasp > nbas) then
        do ib = nbas+1, nbasp
          ioffv0(ib)=nvi
          nvi=nvi+(lmaxp+1)**2
        enddo
        ioffv0(nbasp+1)=nvi
      endif
      if(ipr >= 50) then
        write(6,660) 'ioff  ',(ioff(ib),ib=1,nbas+1)
        write(6,660) 'ioffv0',(ioffv0(ib),ib=1,nbas+1)
        write(6,660) 'ioffp ',(ioffp(ib),ib=1,nbas+1)
  660   format(1x,a6,10i7:/(7x,10i7))
      endif

c ------ set pnu's for each atom, get valence and sc charge -----
      if(ipr >= 10) write(6,391)
      qval=0d0
      qsmc=0d0
      do 30 is=1,nspec
      qc=0d0
      qs=0d0
      lxa=lmxa(is)
      do 31 l=0,lxa
      konfig=pin(l+1,1,is)
      if(lsc == 1) then
        konfiz=piz(l+1,1,is)
        if(konfiz < konfig) qs=qs+2*(2*l+1)
        endif
   31 qc=qc+(konfig-l-1)*2*(2*l+1)
      qv=z(is)-qc
      do 32 ib=1,nbas
      if(ips(ib) == is) then
        qval=qval+qv
        qsmc=qsmc+qs
        do 33 l=0,lxa
        if(lsc == 1) pz(l+1,1,ib)=piz(l+1,1,is)
        if(lsc == 1) pz(l+1,nsp,ib)=piz(l+1,nsp,is)
        p(l+1,1,ib)=pin(l+1,1,is)
   33   p(l+1,nsp,ib)=pin(l+1,nsp,is)
      endif
   32 continue
      if(ipr >= 10) write(6,390) spid(is),z(is),qc,qv,
     .  (pin(k,1,is),k=1,lxa+1)
      if(ipr >= 10 .and. nsp == 2) write(6,392) (pin(k,2,is),k=1,lxa+1)
   30 continue
  390 format(1x,a4,f5.1,2f6.1,2x,6f9.3)
  392 format(17x,'spin 2:',6f9.3)
  391 format(/' spec    Z    qc    qv       pnu ...')
      if(ipr >= 20) write(6,870) qval,qsmc
  870 format(' qval=',f12.6,'   qsmc=',f12.6)

c ------ dimension for eigenvalue problem and augmentation -----
      nla=0
      nll=0
      do 14 i=1,nel
   14 ndim(i)=0
      do 15 ib=1,nbas
      is=ips(ib)
      nla=nla+(lmxa(is)+1)**2
      nll=nll+(lmxl(is)+1)**2
      do 15 i=1,nel
   15 ndim(i)=ndim(i)+(lphi(i,is)+1)**2
      nhs=0
      do 16 i=1,nel
   16 nhs=nhs+ndim(i)

      if(ipr >= 10) write(6,333) nri,nvi,nla,nll,nhs,ndim
  333 format(/' nri=',i5,'     nvi=',i5,'     nla=',i5,'     nll=',i5,
     .   /' secmat dimension=',i5,'      ndim(1..nel)=',6i5)
      if(ipr >= 10) write(6,334) qval,amp
  334 format(' total valence charge=',f8.2
     .     /' Positions calculated for   amp=',f10.6)

c ------ write into log -------------------
CL      write(71,710) nri,nvi,nla,nhs,(ndim(i),i=1,nel)
CL      write(71,713) (el(i),i=1,nel)
  710 format(' mce  nri',i5,'  nvi',i5,'  nla',i5,'  nhs',i5,
     .   ' :',6i5)
  713 format(' mce  el',6f8.3)
      end
