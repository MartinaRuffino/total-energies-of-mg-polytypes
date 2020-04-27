      subroutine hmat2(el,nel,b,bb,hvecs,hi,nbas,lphi,lmxa,ips,
     .  gh,gj,puu,pus,pss,ioffp,pjjs,njj,n0,nhs,ndim,nla,h)
C- Hamiltonian matrix (real h)
      implicit none
      integer n0,nbas,nel,nhs,nla,njj,
     .  lmxa(1),ips(1),lphi(n0,1),ndim(1),ioffp(1)
      double precision hvecs(nla,4,1),hi(6),el(1),
     .  gh(nla,2,nel),gj(nla,2,nel),pjjs(njj,1),puu(1),pus(1),pss(1)
      double precision b(nla,nhs),bb(nla,nhs),h(nhs,nhs),xxx
      integer i,ie,iv,jv,j,je,ldiag,ne1,ne2,owk,
     .  koff,iof1,iof2,kb,ks,nlm1,nlm2,nlma,nx,k0,io
      parameter (nx=49)
      double precision pkk(nx*nx),pkj(nx*nx),pjk(nx*nx),pjj(nx*nx)
      real w(1)
      common /w/ w

C --- Loop over energy pairs and basis ---
      iv=0
      ne1=1
      do  20  ie = 1, nel
        ne2 = ne1
        do  30  je = ie, nel
        iv = iv+1
        ldiag = 0
        if (ie == je) ldiag = 1
        koff = 0
        iof1 = 0
        iof2 = 0
        do  10  kb = 1, nbas
        ks = ips(kb)
        nlm1 = (lphi(ie,ks)+1)**2
        nlm2 = (lphi(je,ks)+1)**2
        nlma = (lmxa(ks)+1)**2

C ---   Make perturbation matrices ---
        if (nlma > nx) call rx('hmat2: nlma gt nx')
        io = ioffp(kb)+1
        k0 = koff+1
        xxx = 0d0
        if (ldiag == 0) xxx = hi(iv)/(el(je)-el(ie))
        call pertmb(puu(io),pus(io),pss(io),hvecs(k0,1,iv),nla,xxx,
     .    gh(k0,1,ie),gh(k0,2,ie),gh(k0,1,je),gh(k0,2,je),
     .    gj(k0,1,ie),gj(k0,2,ie),gj(k0,1,je),gj(k0,2,je),
     .    nlm1,nlm2,nlma,pkk,pkj,pjk,pjj)
        call daxpy(nlma**2,-1d0,pjjs(io,iv),1,pjj,1)

C ---   Add to Hamiltonian this (ie,je,kb) block ---
        call defrr(owk,   nlma*ndim(ie))
        if (ndim(ie) > 0 .and. ndim(je) > 0)
     .    call hmatm1(ips,nbas,b(1+koff,ne1),b(1+koff,ne2),
     .    bb(1+koff,ne1),nla,nla,ndim(ie),ndim(je),
     .    nlm1,nlm2,iof1,iof2,nlma,hi(iv),
     .    pkk,pkj,pjk,pjj,w(owk),ldiag,nhs,h(ne1,ne2))
        call rlse(owk)
        iof1 = iof1+nlm1
        iof2 = iof2+nlm2
        koff = koff+nlma
   10   continue
        ne2 = ne2+ndim(je)
   30 continue
      ne1 = ne1+ndim(ie)
   20 continue

C --- copy into other triangle ---
      do  11  j = 1, nhs
      do  11  i = 1, j
   11 h(j,i) = h(i,j)
      end
      subroutine hmatm1(ips,nbas,b1,b2,bb,nla1,nla2,nlb1,nlb2,
     .  nlm1,nlm2,iof1,iof2,nlma,hi,pkk,pkj,pjk,pjj,wk,ldiag,nhs,h)
C- Adds contribution of one atom to an (e1,e2)-block of hamiltonian.
C  ldiag=1 if diagonal block
      implicit none
      integer nbas,nla1,nla2,nlb1,nlb2,nhs,ldiag,ips(1),
     .  nlm1,nlm2,nlma,iof1,iof2
      double precision pkk(nlm1,nlm2),pkj(nlm1,nlma),
     .  pjk(nlma,nlm2),pjj(nlma,nlma),hi
      double precision b1(nla1,nlb1),b2(nla2,nlb2),bb(nla1,1),
     .  h(nhs,1),wk(nlb1,1)
      integer i,itop,j,ilm1,ilm2,klm

C --- b1(+) * pjj * b2, using wk for b1(+) pjj ---
        call dgemm('T','N',nlb1,nlma,nlma,1d0,b1,nla1,
     .    pjj,nlma,0d0,wk,nlb1)
C        call dgemm('N','N',nlb1,nlb2,nlma,1d0,wk,nlb1,
C     .    b2,nla2,1d0,h,nhs)
        do  11  j = 1, nlb2
          itop = ldiag*j + (1-ldiag)*nlb1
          do  11  klm = 1, nlma
          do  14  i = 1, itop
   14     h(i,j) = h(i,j) + wk(i,klm)*b2(klm,j)
   11   continue

C --- Add  pkj * b2 ---
        if (nlm1 > 0)
     .    call dgemm('N','N',nlm1,nlb2,nlma,1d0,pkj,nlm1,
     .    b2,nla2,1d0,h(1+iof1,1),nhs)

C --- Add  b1(+) pjk ---
        if (nlm2 > 0)
     .    call dgemm('T','N',nlb1,nlm2,nlma,1d0,b1,nla1,
     .    pjk,nlma,1d0,h(1,1+iof2),nhs)

C --- Add pkk ---
        do  42  ilm2 = 1, nlm2
        itop = ldiag*ilm2 + (1-ldiag)*nlm1
        do  41  ilm1 = 1, itop
   41   h(ilm1+iof1,ilm2+iof2) = h(ilm1+iof1,ilm2+iof2) + pkk(ilm1,ilm2)
   42   continue

C --- Add global terms for ldiag=1 ---
        if (ldiag == 1) then
          do  57  ilm1 = 1, nlm1
          i = iof1+ilm1
          do  57  j = i, nlb2
   57     h(i,j) = h(i,j) + hi*bb(ilm1,j)
        endif

      end
