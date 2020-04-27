      subroutine smatm(el,nel,b,bb,svecs,si,nbas,lphi,lmxa,ips,
     .  n0,nhs,ndim,nla,s)
c  lmto-asa matrix (overlap or hamiltonian)
      implicit real*8 (a-h,p-z), integer(o)
      dimension svecs(4*nla,1),si(6),el(1),
     .  lmxa(1),ips(1),lphi(n0,1),ndim(1)
      real*8 b(nla,nhs),bb(nla,nhs),s(nhs,nhs)
      real w(1)
      common /w/ w

      iv=0
      ne1=1
      do 20 ie=1,nel
        ne2=ne1
        do 30 je=ie,nel
        iv=iv+1
        ldiag=0
        if (ie == je) ldiag=1
        call smatm1(lmxa,ips,nbas,b(1,ne1),b(1,ne2),n0,nla,
     .    ndim(ie),ndim(je),lphi(ie,1),lphi(je,1),el(ie),el(je),
     .    svecs(1,iv),si(iv),bb(1,ne1),s(ne1,ne2),nhs,ldiag)
        ne2=ne2+ndim(je)
   30 continue
      ne1=ne1+ndim(ie)
   20 continue

c --- copy into other triangle ---
      do 11 j=1,nhs
      do 11 i=1,j
   11 s(j,i)=s(i,j)
      end
c ---------- smatm1 -----------------------
      subroutine smatm1(lmxa,ips,nbas,b1,b2,n0,nla,nlb1,nlb2,
     .  lphi1,lphi2,e1,e2,fvec,si,bb,s,nh,ldiag)
c  makes one (e1,e2)-block of overlap matrix
c  ldiag=1 if a diagonal block is to be made
      implicit real*8 (a-h,p-z), integer(o)
      dimension lmxa(1),ips(1),lphi1(n0,1),lphi2(n0,1),fvec(nla,4)
      real*8 b1(nla,nlb1),b2(nla,nlb2),s(nh,1),cx,bb(nla,1)
      if(nlb1 <= 0.or.nlb2 <= 0) return
      xxx=0.d0
      if(ldiag == 0) xxx=si/(e2-e1)
      do 1 j=1,nlb2
      itop=ldiag*j+(1-ldiag)*nlb1
      do 1 i=1,itop
  1   s(i,j)=0.d0
c ----------- b1(dagger) * fjj * b2 ----------------
      koff=0
      do 10 kb=1,nbas
        ks=ips(kb)
        nlma=(lmxa(ks)+1)**2
        do 11 j=1,nlb2
        itop=ldiag*j+(1-ldiag)*nlb1
        do 11 klm=1,nlma
          k=klm+koff
          cx=fvec(k,4)*b2(k,j)
          do 14 i=1,itop
   14     s(i,j)=s(i,j)+b1(k,i)*cx
   11   continue
   10 koff=koff+nlma
c ----------- add fhh ----------------------------
      iof1=0
      iof2=0
      koff=0
      do 20 kb=1,nbas
        ks=ips(kb)
        nlm1=(lphi1(1,ks)+1)**2
        nlm2=(lphi2(1,ks)+1)**2
        nlma=(lmxa(ks) +1)**2
        do 21 ilm=1,min0(nlm1,nlm2)
   21   s(ilm+iof1,ilm+iof2)=s(ilm+iof1,ilm+iof2)+fvec(ilm+koff,1)
        iof1=iof1+nlm1
        iof2=iof2+nlm2
   20 koff=koff+nlma
c ----------- add  fhj * b2 -----------------------
      iof1=0
      koff=0
      do 30 kb=1,nbas
        ks=ips(kb)
        nlm1=(lphi1(1,ks)+1)**2
        nlma=(lmxa(ks) +1)**2
        do 31 ilm=1,nlm1
        i=iof1+ilm
        k=koff+ilm
        jbot=ldiag*i+(1-ldiag)*1
        do 31 j=jbot,nlb2
   31   s(i,j)=s(i,j)+(fvec(k,2)+xxx)*b2(k,j)
        iof1=iof1+nlm1
   30 koff=koff+nlma
c ----------- add  b1(dagger) * wjh  --------------
      iof2=0
      koff=0
      do 40 kb=1,nbas
        ks=ips(kb)
        nlm2=(lphi2(1,ks)+1)**2
        nlma=(lmxa(ks) +1)**2
        do 41 jlm=1,nlm2
        j=iof2+jlm
        k=koff+jlm
        itop=ldiag*j+(1-ldiag)*nlb1
        do 41 i=1,itop
   41   s(i,j)=s(i,j)+(fvec(k,3)-xxx)*b1(k,i)
        iof2=iof2+nlm2
   40 koff=koff+nlma
c ----------- add global terms for ldiag=1 ---------
      if(ldiag == 1) then
      koff=0
      iof1=0
      do 36 ib=1,nbas
      is=ips(ib)
      nlm1=(lphi1(1,is)+1)**2
      nlma=(lmxa(is) +1)**2
      do 37 ilm=1,nlm1
      i=iof1+ilm
      k=koff+ilm
      do 37 j=i,nlb2
   37 s(i,j)=s(i,j)+si*bb(k,j)
      iof1=iof1+nlm1
   36 koff=koff+nlma
      endif
      end
