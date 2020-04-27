      subroutine vsstal(nbas,ips,w,nla,gh,gj,n0,nel,
     .  ioffa,ioffb,ioffl,llink,b,bl,ceadd,t,v,s)
C- Values and slopes for one eigenvector on all spheres, linked basis
      implicit real*8 (a-h,p-z), integer(o)
      dimension gh(nla*2,nel),gj(nla*2,nel),ips(1),
     .   ioffa(nbas),ioffb(nbas,nel),ioffl(1),ceadd(25,5,1)
      double precision b(nla,1),bl(nla,1),t(1),s(nla),v(nla),w(nla,3)

      do 10 ila=1,nla
      v(ila)=0d0
   10 s(ila)=0d0

      if (llink == 0) then
        ne=1
        do 20 ie=1,nel
          ndim = ioffb(nbas+1,ie)-ioffb(1,ie)
          call hhvst(nbas,ioffa,ioffb(1,ie),w,nla,n0,
     .      t(ne),b(1,ne),gh(1,ie),gj(1,ie),ndim,v,s)
          ne=ne+ndim
   20   continue
      else
        ne=1
        do 30 ie=1,nel
          ndim = ioffb(nbas+1,ie)-ioffb(1,ie)
          call hmcadd(ips,nbas,ioffb,ie,ceadd,w(1,3))
          call hhvsl(nbas,w,w(1,2),nla,n0,t(ne),b(1,ne),bl,
     .      gh(1,ie),gj(1,ie),gh(1,nel+1),gj(1,nel+1),
     .      w(1,3),ioffa,ioffb(1,ie),ioffl,ndim,v,s)
          ne=ne+ndim
   30   continue
      endif
      end
      subroutine hhvsl(nbas,w1,wl,nla,n0,t,
     .  b,bl,gh,gj,ghl,gjl,add,ioffa,ioffb,ioffl,ndim,v,s)
C- Adds to v,s for orbitals with one energy, linked basis
      implicit real*8 (a-h,p-z)
      dimension ioffb(1),ioffl(1),ioffa(1)
      dimension t(ndim),v(nla),s(nla),w1(nla),wl(nla),
     .  gh(nla,2),gj(nla,2),ghl(nla,2),gjl(nla,2),
     .  b(nla,ndim),bl(nla,ndim),add(ndim)
      if(ndim <= 0) return
c ------------- add part from hankels --------
      koff=0
      ioff=0
      do 10 kb=1,nbas
      nlmb=ioffb(kb+1)-ioffb(kb)
      nlma=ioffa(kb+1)-ioffa(kb)
      do 11 ilm=1,nlmb
      k=ilm+koff
      i=ilm+ioff
      v(k)=v(k)+t(i)*(gh(k,1)+add(i)*ghl(k,1))
   11 s(k)=s(k)+t(i)*(gh(k,2)+add(i)*ghl(k,2))
      koff=koff+nlma
   10 ioff=ioff+nlmb
c ------------- add part from bessels --------------
      call dpzero(w1,nla)
      call dpzero(wl,nla)
      do 20 ib=1,nbas
      nlmb=ioffb(ib+1)-ioffb(ib)
      do 20 ilm=1,nlmb
      i=ioffb(ib)-ioffb(1)+ilm
      l=ioffl(ib)-ioffl(1)+ilm
      xxi = add(i)*t(i)
      do 20 k=1,nla
      wl(k)=wl(k)+bl(k,l)*xxi
      w1(k)=w1(k)+b(k,i)*t(i)
   20 continue
      do 30 k=1,nla
      v(k)=v(k)+gj(k,1)*w1(k)+gjl(k,1)*wl(k)
   30 s(k)=s(k)+gj(k,2)*w1(k)+gjl(k,2)*wl(k)
      end
      subroutine hhvst(nbas,ioffa,ioffb,w,nla,n0,t,b,gh,gj,ndim,
     .  v,s)
C- adds to v,s for orbitals with one energy
      implicit real*8 (a-h,p-z)
      dimension gh(nla,2),gj(nla,2),ioffa(1),ioffb(1)
      dimension t(ndim),v(nla),s(nla),w(nla),b(nla,ndim)
      if(ndim <= 0) return
c ------------- add part from hankels --------
      koff=0
      ioff=0
      do 10 kb=1,nbas
      nlmb=ioffb(kb+1)-ioffb(kb)
      nlma=ioffa(kb+1)-ioffa(kb)
      do 11 ilm=1,nlmb
      k=ilm+koff
      i=ilm+ioff
      v(k)=v(k)+t(i)*gh(k,1)
   11 s(k)=s(k)+t(i)*gh(k,2)
      koff=koff+nlma
   10 ioff=ioff+nlmb
c ------------- add part from bessels --------------
      do 18 k=1,nla
   18 w(k)=0d0
      do 20 i=1,ndim
      do 20 k=1,nla
   20 w(k)=w(k)+b(k,i)*t(i)
      do 21 k=1,nla
      v(k)=v(k)+gj(k,1)*w(k)
   21 s(k)=s(k)+gj(k,2)*w(k)
      end
c --------- sub hhadq1 ---------------------------------
      subroutine hhadq1(wgt,nlma,v,s,quu,qus,qss)
      implicit real*8 (a-h,p-z)
      dimension v(1),s(1),quu(nlma,1),qus(nlma,1),qss(nlma,1)
      do 12 jlm=1,nlma
      vw  = v(jlm)*wgt
      sw  = s(jlm)*wgt
      sw2 = sw*2
      do 12 ilm=1,nlma
      quu(ilm,jlm)=quu(ilm,jlm)+v(ilm)*vw
      qus(ilm,jlm)=qus(ilm,jlm)+v(ilm)*sw2
   12 qss(ilm,jlm)=qss(ilm,jlm)+s(ilm)*sw
      end
c --------- sub hhadq0 ---------------------------------
      subroutine hhadq0(wgt,lmaxa,v,s,n0,qus0)
      implicit real*8 (a-h,p-z)
      dimension v(1),s(1),qus0(3,0:1)
      i=0
      do 12 l=0,lmaxa
      do 12 m=-l,l
      i=i+1
      qus0(1,l)=qus0(1,l)+v(i)*v(i)*wgt
      qus0(2,l)=qus0(2,l)+v(i)*s(i)*wgt*2d0
   12 qus0(3,l)=qus0(3,l)+s(i)*s(i)*wgt
      end
