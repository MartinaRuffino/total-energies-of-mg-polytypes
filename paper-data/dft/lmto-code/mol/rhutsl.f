      subroutine rhutsl(wgt,ewt,evl,sumev,nstate,nbas,isp,lmxa,ips,
     .  n0,nhs,nla,ioffa,ioffb,ioffl,ioffp,llink,b,bl,ceadd,t,nel,gh,gj,
     .  v,s,quu,qus,qss,qus0)
C- Accumulate output density inside spheres, linked basis
c  v,s are workspaces of length nla, d of length nla*nel.
C  6 Jan 95 accumulates quu etc for spin isp (MvS)
      implicit real*8 (a-h,p-z), integer(o)
      dimension lmxa(1),gh(1),gj(1),ips(1),quu(1),qus(1),qss(1),
     .  t(nhs,1),b(nla,1),v(nla),s(nla),ewt(nhs,isp),
     .  evl(nhs,1),qus0(3,n0,nbas,2),ioffa(nbas),ioffb(nbas,1),ioffl(1),
     .  ioffp(1)
      real w(1)
      common /w/ w
      call getpr(ipr)
      call tcn('rhutsl')
      nsp=lsp()+1
      istate = 1
      jstate = nstate
C ... Third column of wk must be at least max(ndim)
      call defrr(owk,      nla*3)

c ------ start loop over states ------
      sumev=0d0
      do 80 is=istate,jstate
        wi=wgt*ewt(is,isp)
        sumev=sumev+wi*evl(is,isp)
        call vsstal(nbas,ips,w(owk),nla,gh,gj,n0,
     .    nel,ioffa,ioffb,ioffl,llink,b,bl,ceadd,t(1,is),v,s)
        k1=1
        do 50 ib=1,nbas
          lmaxa=lmxa(ips(ib))
          nlma=(lmaxa+1)**2
C i1 is index to quu(ib,isp)
          i1 = 1+(isp-1)*ioffp(nbas+1)+ioffp(ib)
          call hhadq0(wi,lmaxa,v(k1),s(k1),n0,qus0(1,1,ib,isp))
          call hhadq1(wi,nlma,v(k1),s(k1),quu(i1),qus(i1),qss(i1))
  50      k1=k1+nlma
  80    continue
      call rlse(owk)
CL      write(71,710) (evl(is,isp),is=istate,jstate)
  710 format(' mce  e=',8f8.4)
      call tcx('rhutsl')

c ------ printout quu, qus, qss -----------------
      if(ipr <= 60) return
      k=(isp-1)*ioffp(nbas+1)
      do 40 ib=1,nbas
      if(nsp == 1) write(6,220) ib
      if(nsp == 2) write(6,220) ib,isp
      nlma=(lmxa(ips(ib))+1)**2
      do 40 j=1,nlma
      do 40 i=1,nlma
      k=k+1
      top=dmax1(dabs(quu(k)),dabs(qus(k)),dabs(qss(k)))
      if(top >= 1d-5) write(6,221) i,j,quu(k),qus(k),qss(k)
  221 format(2i5,3f12.6)
  220 format(/'  ilm  jlm',7x,'quu',9x,'qus',9x,'qss  for ib=',i3:
     .  '  spin',i2)
  40  continue

      end
