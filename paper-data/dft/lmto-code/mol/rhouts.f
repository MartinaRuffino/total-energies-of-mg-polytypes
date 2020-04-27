      subroutine rhouts(wgt,nbas,lmxa,ips,lphi,n0,nhs,nla,
     .   b,t,evl,qval,nel,ndim,gh,gj,sumev,v,s,
     .   quu,qus,qss,qus0,d,njj,qjjs,nv,qjj0)
c  Accumulate output density inside spheres.
c  v,s are workspaces of length nla, d of length nla*nel.
      implicit real*8 (a-h,p-z), integer(o)
      dimension lmxa(1),evl(1),gh(1),gj(1),ips(1),lphi(n0,1),
     .  quu(1),qus(1),qss(1),t(nhs,1),b(nla,1),ndim(1),v(nla),s(nla),
     .  qus0(3,n0,nbas),qjjs(njj,1),d(nla,nel),qjj0(nv,n0,1)
      real w(1)
      common /w/ w
      call getpr(ipr)
      call tcn('rhouts')
      nqval=qval+0.01
      nstate=nqval/2
      istate = 1
      jstate = nstate
      call defrr(owk,      nla)

c ------ start loop over states ------
      do 80 is=istate,jstate
        ei=evl(is)
        sumev=sumev+wgt*ei
        call vsstat(nbas,ips,lmxa,w(owk),nla,gh,gj,lphi,n0,ndim,
     .     nel,b,t(1,is),v,s)
        call tlstat(nla,ndim,nel,b,t(1,is),d)
        i1=1
        k1=1
        do 50 ib=1,nbas
          lmaxa=lmxa(ips(ib))
          nlma=(lmaxa+1)**2
          call hhadq0(wgt,lmaxa,v(k1),s(k1),n0,qus0(1,1,ib))
          call hhadj0(wgt,lmaxa,nla,nel,d(k1,1),n0,nv,qjj0(1,1,ib))
          call hhadq1(wgt,nlma,v(k1),s(k1),quu(i1),qus(i1),qss(i1))
          call hhadqj(wgt,nlma,nla,nel,d(k1,1),njj,qjjs(i1,1))
          i1=i1+nlma*nlma
  50      k1=k1+nlma
  80    continue

      call rlse(owk)
CL      write(71,710) (evl(is),is=istate,jstate)
  710 format(' evl=',9f8.4)
      call rlse(owk)
c ------ printout quu etc -----------------------
      if(ipr <= 60) return
      k=0
      do 40 ib=1,nbas
      write(6,220) ib
      nlma=(lmxa(ips(ib))+1)**2
      do 40 j=1,nlma
      do 40 i=1,nlma
      k=k+1
      top=dmax1(dabs(quu(k)),dabs(qus(k)),dabs(qss(k)))
      if(top >= 1d-5) write(6,221) i,j,quu(k),qus(k),qss(k)
  221 format(2i5,3f12.6)
  220 format(/'  ilm  jlm',7x,'quu',9x,'qus',9x,'qss     for  ib=',i3)
  40  continue
      nv=(nel*(nel+1))/2
      k=0
      do 42 ib=1,nbas
      write(6,320) ib
      nlma=(lmxa(ips(ib))+1)**2
      do 42 j=1,nlma
      do 42 i=1,nlma
      k=k+1
      top=0d0
      do 43 iv=1,nv
  43  top=dmax1(top,dabs(qjjs(k,iv)))
      if(top >= 1d-5) write(6,221) i,j,(qjjs(k,iv),iv=1,nv)
  321 format(2i5,6f12.6)
  320 format(/'  ilm  jlm',7x,'qjjs(1..nv)    for ib=',i3)
  42  continue

      call tcx('rhouts')
      end
