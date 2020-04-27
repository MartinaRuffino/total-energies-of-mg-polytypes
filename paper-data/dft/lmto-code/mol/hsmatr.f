C#define BLAS3
      subroutine hsmatr(lhs,el,nel,b,bd,bl,bld,ceadd,hvecs,hi,svecs,si,
     .  nbas,isp,n0,nphi,lphi,lmxa,cg,cy,alat,indxcg,jcg,pos,
     .  iprs,ntab,iax,rtab,ntab3,iax3,rtab3,iofh,ixp,dxp,
     .  ips,gh,gj,puu,pus,pss,ioffa,ioffb,ioffl,iofc,iofcl,iofc3,ioffp,
     .  nhs,llink,nla,h,s,lstore)
C- Real-space ASA+warp hamiltonian and overlap matrix, linked basis.
C  Storage of h (and s): h stored as a sequence of (nlmi,nlmj) blocks.
C    Only the upper triangle of blocks is stored.  iofh(iv,pair) holds
C    offset to start of (ie,je,pair) block, with iv a combined (ie,je)
C    index, given by iv=je+(ie-1)*(2*nel-ie)/2 for je>=ie
C  iofc,iofcl: local work arrays of dimension at least (mxcsiz,nel)
C  lhs ones digit nonzero: calculate s; tens digit nonzero: calculate h
C  b,bd must be allocated max(nlma)*nhs !?
C  lstore: 0 stores h has whole array, 1 in sparse form.
      implicit none
      integer lhs,nbas,isp,nel,n0,nhs,llink,nla,ips(1),nphi(1),
     .  lphi(n0,1),lmxa(1),ixp(1),ioffp(1),iofc(1),iofc3(1),lstore,
     .  iofcl(nbas+1),ioffa(1),ioffb(nbas,nel),ioffl(nbas,1),iofh(1),
     .  jcg(1),indxcg(1),iprs(1),ntab(1),iax(10,1),ntab3(1),iax3(40)
      double precision hvecs(nla,4,1),hi(6),el(1),ceadd(25,5,1),alat,
     .  cg(1),cy(1),gh(nla,2,nel),gj(nla,2,nel),puu(1),pus(1),pss(1),
     .  pos(3,nbas),dxp(1),svecs(nla,4,1),si(6),rtab(3,1),rtab3(3,1)
      double precision b(1),bl(1),bd(1),bld(1),h(nhs,nhs),s(nhs,nhs),
     .  xxh,xxs,xx2,xx3,xxt,alatl
      integer i,ie,iv,ivs,j,je,k,ic,ic3,ldiag,ofe1,ofe2,owk,owk2,iofi,
     .  kb,ks,nlm1,nlm2,nx,k0,io,je1,ivh,ivl,jvl,ivls,jvls,il,ilfi,
     .  nlai,nlbi,nlmai,nlmbs,nlmak,nclus,nclus3,nlbx,nvl,nvls,lsp,nsp,
     .  oadd1,oadd2,obl1,obl2,obl3i,obl3j,ndimi,ndimj,ndim3i,ndim3j,
     .  opkk,opkj,opjk,opjj,opkk2,opkj2,opjk2,opjj2,
     .  opkkt,opkjt,opjkt,opjjt,opkk3,opkj3,opjk3,opjj3,
     .  ojkl(0:20),oba(10),obh(10),obp,obl,oblp,oiofa,
     .  ob3i,ob3j,oips
      logical lh,ls
c     character dnam*4, outs*80

      parameter (nx=49)
      real w(1)
      common /w/ w

c     dnam = '.'

C --- Setup  ---
Cshould be no kb entering at all after hoffr, except references of hvec
      call tcn('hsmatr')
c     call wkprnt(1)
      nsp = lsp()+1
      alatl = 1d0
      lh = mod(lhs/10,10) /= 0
      ls = mod(lhs,10) /= 0
      call defrr(opkk, nx*nx)
      call defrr(opkj, nx*nx)
      call defrr(opjk, nx*nx)
      call defrr(opjj, nx*nx)
      if (llink > 0) then
        call defrr(opkk2, nx*nx)
        call defrr(opkj2, nx*nx)
        call defrr(opjk2, nx*nx)
        call defrr(opjj2, nx*nx)
        call defrr(opkkt, nx*nx)
        call defrr(opkjt, nx*nx)
        call defrr(opjkt, nx*nx)
        call defrr(opjjt, nx*nx)
        call defrr(opkk3, nx*nx)
        call defrr(opkj3, nx*nx)
        call defrr(opjk3, nx*nx)
        call defrr(opjj3, nx*nx)
        il = nel+1
        nvl = ((nel+1)*(nel+2))/2
        nvls= nsp*(nvl-1)+isp
      endif

C --- Loop over basis ---
      do  10  kb = 1, nbas
C ...   iofc,iofl: table of offsets for cluster strux
        ic  = 1+ntab(kb)
        ic3 = 1+ntab3(kb)
        do  11  j = ntab(kb)+1, ntab(kb+1)
   11   iprs(j-ntab(kb)) = ips(iax(2,j))
        nclus = ntab(kb+1) - ntab(kb)
        nclus3 = ntab3(kb+1)-ntab3(kb)
        call defrr(oadd1,25*nclus)
        call defrr(oadd2,25*nclus)
        call defrr(oips, nclus3**2)
        call pairmc(kb,nclus3,iax,iax3(ic3),ntab,rtab,w(oips))
        call defrr(oiofa, nclus+1)
        call hoffs(nclus,lphi,lmxa,n0,iprs,nel,iofc,w(oiofa),iofcl,nlbx)
        call hoff3(nclus3,iax3(ic3),iax(1,ic),lphi,lmxa,
     .    n0,ips,nel,iofc3)

        call rlse(oiofa)
C ...   Dimensions of b,bd are locally b(nlmak,*)
        nlmak = ioffa(kb+1) - ioffa(kb)
        k0   = 1+ioffa(kb)
        if (nlmak > nx) call rx('hsmatr: nlmak gt nx')
        call astrxj(nphi,lphi,nel,el,nlmak,rtab(1,ic),n0,nclus,1,iprs,
     .    alatl,cg,jcg,indxcg,cy,llink,iofc,iofcl,0,dxp,b,bd,bl,bld)
C       call prmx('b',b,nlmak,nlmak,nhs)

C ---   Loop over energy pairs ---
C ...  compound (ie,je) indices:
C iv:  index to (ie,je) pair for hi,si,hvecs,svecs
C ivh: index to (ie,je) for hamiltonian.
C      They differ in case of linked basis to make room for extra
C      potential needed from linking part of matrix elements.
C ivs: index to (isp,ie,je) for hvecs,svecs
C ivl,ivls:  position of (ie,linked) or (isp,ie,linked) block
C jvl,jvls:  position of (je,linked) or (isp,je,linked) block
        do  20  ie = 1, nel
C ...   Local cluster dimensions for this (ie,kb)
        ivl = (ie*(il+il-ie+1))/2
        ivls= nsp*(ivl-1)+isp
        ndimi = iofc(1+ie*nclus)   - iofc(1+(ie-1)*nclus)
        ndim3i= iofc3(1+ie*nclus3) - iofc3(1+(ie-1)*nclus3)
        ofe1 = iofc(1+(ie-1)*nclus)
        nlm1 = iofc(2+(ie-1)*nclus) - iofc(1+(ie-1)*nclus)

C ...   3C parts of strux for this ie,kb
        ofe1 = iofc(1+(ie-1)*nclus)
        call defrr(ob3i,    nlmak*ndim3i)
        call pmblk(nclus3,iax3(ic3),iofc(1+(ie-1)*nclus),
     .    iofc3(1+(ie-1)*nclus3),nlmak,b(1+nlmak*ofe1),111,
     .    0d0,nlmak,nlmak,w(ob3i),k)
        if (k /= ndim3i) call rx('bug in hsmatr')

C ...   Linked strux for this ie,kb
        if (llink > 0) then
          call hmcadd(iprs,nclus,iofc,ie,ceadd,w(oadd1))
          call defrr(obl1,   nlmak*ndimi)
          call pmblk(nclus,0,iofcl,iofc(1+(ie-1)*nclus),
     .      nlmak,bl,10011,w(oadd1),nlmak,nlmak,w(obl1),k)
          if (k /= ndimi) call rx('bug in hsmatr')

C         call prmx('bl1 in hsmatr',w(obl1),nlmak,nlmak,ndimi)

          call defrr(obl3i,    nlmak*ndim3i)
          call pmblk(nclus3,iax3(ic3),iofc(1+(ie-1)*nclus),
     .      iofc3(1+(ie-1)*nclus3),nlmak,w(obl1),111,0d0,nlmak,nlmak,
     .      w(obl3i),k)
          if (k /= ndim3i) call rx('bug in hsmatr')

C          call prmx('bl3i in hsmatr',w(obl3i),nlmak,nlmak,ndim3i)

        endif

        do  30  je = ie, nel
C ...   Local cluster dimensions for this (je,kb)
        jvl = (je*(il+il-je+1))/2
        jvls= nsp*(jvl-1)+isp
        ndimj = iofc(1+je*nclus)   - iofc(1+(je-1)*nclus)
        ndim3j= iofc3(1+je*nclus3) - iofc3(1+(je-1)*nclus3)
        ofe2 = iofc(1+(je-1)*nclus)
        nlm2 = iofc(2+(je-1)*nclus) - iofc(1+(je-1)*nclus)
        iv = je+(ie-1)*(2*nel-ie)/2
        ivh = 1+(iv-1)*ntab(nbas+1)
        if (llink > 0) iv = iv + ie-1
        ivs = nsp*(iv-1)+isp
        ldiag = 0
        if (ie == je) ldiag = 1

C        call awrit5('hsm ie,je,kb %i %i %,2i ndim3i,j %,3i %,3i',
C     .    ' ',80,6,ie,je,kb,ndim3i,ndim3j)

C ...   3C parts of strux for this je,kb
        call defrr(ob3j,    nlmak*ndim3j)
        call pmblk(nclus3,iax3(ic3),iofc(1+(je-1)*nclus),
     .    iofc3(1+(je-1)*nclus3),nlmak,b(1+nlmak*ofe2),111,
     .    0d0,nlmak,nlmak,w(ob3j),k)
        if (k /= ndim3j) call rx('bug in hsmatr')
C ...   Linked strux for this je,kb
        if (llink > 0) then
          call hmcadd(iprs,nclus,iofc,je,ceadd,w(oadd2))
          call defrr(obl2,   nlmak*ndimj)
C         call xxmcsx(nclus,nlmak,nlmak,iofc(1+(je-1)*nclus),iofcl,bl,
C     .     w(obl2),w(oadd2))
          call pmblk(nclus,0,iofcl,iofc(1+(je-1)*nclus),
     .      nlmak,bl,10011,w(oadd2),nlmak,nlmak,w(obl2),k)
          if (k /= ndimj) call rx('bug in hsmatr')
          call defrr(obl3j,    nlmak*ndim3j)
          call pmblk(nclus3,iax3(ic3),iofc(1+(je-1)*nclus),
     .      iofc3(1+(je-1)*nclus3),nlmak,w(obl2),111,
     .      0d0,nlmak,nlmak,w(obl3j),k)
          if (k /= ndim3j) call rx('bug in hsmatr')
        endif

C ---   Perturbation matrices for H for this kb ---
        if (lh) then
          io = nsp*ioffp(kb)+1
          if (isp == 2) io = io + nlmak**2
          xxh = 0d0
          if (ldiag == 0) xxh = hi(iv)/(el(je)-el(ie))
          call pertmb(puu(io),pus(io),pss(io),hvecs(k0,1,ivs),nla,xxh,
     .      gh(k0,1,ie),gh(k0,2,ie),gh(k0,1,je),gh(k0,2,je),
     .      gj(k0,1,ie),gj(k0,2,ie),gj(k0,1,je),gj(k0,2,je),
     .      nlm1,nlm2,nlmak,w(opkk),w(opkj),w(opjk),w(opjj))
C          call awrit3('%x pkk.%i%i%i',outs,len(outs),0,ie,je,kb)
C          call prmx(outs,w(opkk),nlm1,nlm1,nlm2)
C          call awrit3('diff out.o2 xx/pkk.%i%i%i',
C     .      outs,len(outs),0,ie,je,kb)
C          call fsystm(outs,i)
C          call awrit1('shell command '//outs//'%a returned %i',
C     .      outs,len(outs),-6,i)
        endif
        if (ls) then
          xxs = 0d0
          if (ldiag == 0) xxs = si(iv)/(el(je)-el(ie))
        endif

C ---   Add to H,S this (ie,je,kb) block (no linked basis) ---
        if (llink == 0) then
          call defrr(owk,   nlmak*max(ndimi,ndimj))
          call defrr(owk2,  nlmak*max(ndimi,ndimj))
          if (lh .and. nlm1*ndimj+nlm2*ndimi > 0) then
            call xxhm2c(b(1+nlmak*ofe1),b(1+nlmak*ofe2),
     .        bd(1+nlmak*ofe1),nlmak,nlmak,ndimi,ndimj,nlm1,nlm2,0,
     .        nlmak,hi(iv),w(opkk),w(opkj),w(opjk),ldiag,w(owk),w(owk2))
            call hscop(ic,nclus,ie,je,iax,w(oips),0,2,
     .        iofc,iofh(ivh),nlm1,w(owk),nlm2,w(owk2),ndimi,
     .        w(owk),nhs,h)
          endif
          if (ls .and. nlm1*ndimj+nlm2*ndimi > 0) then
            call xxsm2c(b(1+nlmak*ofe1),b(1+nlmak*ofe2),
     .        bd(1+nlmak*ofe1),nlmak,nlmak,ndimi,ndimj,nlm1,nlm2,0,
     .        nlmak,ldiag,si(iv),nla,svecs(k0,1,ivs),xxs,w(owk),w(owk2))
            call hscop(ic,nclus,ie,je,iax,w(oips),0,2,
     .        iofc,iofh(ivh),nlm1,w(owk),nlm2,w(owk2),ndimi,
     .        w(owk),nhs,s)
          endif
          call rlse(owk)
          call defrr(owk,   nlmak*ndim3i*ndim3j)
          call defrr(owk2,  nlmak*max(ndim3i,ndim3j))
          if (lh .and. ndimi*ndimj > 0) then
            call xxhm3c(w(ob3i),w(ob3j),nlmak,nlmak,
     .        ndim3i,ndim3j,nlmak,w(opjj),ldiag,w(owk2),w(owk))
            call hscop(ic,nclus3,ie,je,iax,w(oips),
     .        iax3(ic3),3,iofc3,iofh(ivh),nlm1,w(owk),nlm2,
     .        w(owk2),ndim3i,w(owk),nhs,h)
          endif
          if (ls .and. ndimi*ndimj > 0) then
            call xxsm3c(w(ob3i),w(ob3j),nlmak,nlmak,ndim3i,ndim3j,
     .        nlmak,nla,svecs(k0,1,ivs),ldiag,w(owk2),w(owk))
            call hscop(ic,nclus3,ie,je,iax,w(oips),
     .        iax3(ic3),3,iofc3,iofh(ivh),nlm1,w(owk),nlm2,
     .        w(owk2),ndim3i,w(owk),nhs,s)
          endif
          call rlse(owk)
C ---   Add to H,S this (ie,je,kb) block (linked basis) ---
        else
          call defrr(owk,   nlmak*max(ndimi,ndimj))
          call defrr(owk2,  nlmak*max(ndimi,ndimj))
          if (lh .and. nlm1*ndimj+nlm2*ndimi > 0) then
C ...     pert matrices for (ie,linked) block
            xxh = hi(ivl)/(el(il)-el(ie))
            call pertmb(puu(io),pus(io),pss(io),hvecs(k0,1,ivls),nla,
     .      xxh,gh(k0,1,ie),gh(k0,2,ie),gh(k0,1,il),gh(k0,2,il),
     .      gj(k0,1,ie),gj(k0,2,ie),gj(k0,1,il),gj(k0,2,il),
     .      nlm1,nlm2,nlmak,w(opkk2),w(opkj2),w(opjk2),w(opjj2))
C ...     pert matrices for (je,linked) block
            xxh = hi(jvl)/(el(il)-el(je))
            call pertmb(puu(io),pus(io),pss(io),hvecs(k0,1,jvls),nla,
     .      xxh,gh(k0,1,je),gh(k0,2,je),gh(k0,1,il),gh(k0,2,il),
     .      gj(k0,1,je),gj(k0,2,je),gj(k0,1,il),gj(k0,2,il),
     .      nlm2,nlm1,nlmak,w(opkkt),w(opkjt),w(opjkt),w(opjjt))
C ...     pert matrices for (linked,linked) block
            xxh = 0d0
            call pertmb(puu(io),pus(io),pss(io),hvecs(k0,1,nvls),nla,
     .      xxh,gh(k0,1,il),gh(k0,2,il),gh(k0,1,il),gh(k0,2,il),
     .      gj(k0,1,il),gj(k0,2,il),gj(k0,1,il),gj(k0,2,il),
     .      nlm1,nlm2,nlmak,w(opkk3),w(opkj3),w(opjk3),w(opjj3))
C ...     One and two-center contribution to h from this kb
            call xxhl2c(b(1+nlmak*ofe1),b(1+nlmak*ofe2),
     .        bd(1+nlmak*ofe1),w(obl1),w(obl2),nlmak,nlmak,
     .        ndimi,ndimj,nlm1,nlm2,0,0,nlmak,hi(iv),
     .        w(oadd1),w(oadd2),w(opkk),w(opkj),w(opjk),
     .        w(opkk2),w(opkj2),w(opjk2),
     .        w(opkk3),w(opkj3),w(opjk3),
     .        w(opkkt),w(opkjt),w(opjkt),
     .        ldiag,w(owk),w(owk2))
            call xxhmgl(nclus,hi(nvl),nlmak,bld,w(oadd1),
     .        w(oadd2),iofc(1+(je-1)*nclus),iofcl,nlm1,w(owk))
            call hscop(ic,nclus,ie,je,iax,w(oips),0,2,
     .        iofc,iofh(ivh),nlm1,w(owk),nlm2,w(owk2),ndimi,
     .        w(owk),nhs,h)
          endif
C ...     One and two-center contribution to s from this kb
          if (ls .and. nlm1*ndimj+nlm2*ndimi > 0) then
            xx2 = si(ivl)/(el(il)-el(ie))
            xxt = si(jvl)/(el(il)-el(je))
            xx3 = 0d0
            call xxsl2c(b(1+nlmak*ofe1),b(1+nlmak*ofe2),
     .        bd(1+nlmak*ofe1),w(obl1),w(obl2),nlmak,nlmak,
     .        ndimi,ndimj,nlm1,nlm2,0,0,nlmak,si(iv),
     .        w(oadd1),w(oadd2),nla,svecs(k0,1,ivs),xxs,
     .        svecs(k0,1,ivls),xx2,svecs(k0,1,nvls),xx3,
     .        svecs(k0,1,jvls),xxt,ldiag,w(owk),w(owk2))
            call xxhmgl(nclus,si(nvl),nlmak,bld,w(oadd1),
     .        w(oadd2),iofc(1+(je-1)*nclus),iofcl,nlm1,w(owk))
            call hscop(ic,nclus,ie,je,iax,w(oips),0,2,
     .        iofc,iofh(ivh),nlm1,w(owk),nlm2,w(owk2),ndimi,
     .        w(owk),nhs,s)
          endif
          call rlse(owk)
          call defrr(owk,   nlmak*ndim3i*ndim3j)
          call defrr(owk2,  nlmak*max(ndim3i,ndim3j))
C ...     Three-center contribution to h from this kb
          if (lh .and. ndimi*ndimj > 0) then
            call xxhl3c(w(ob3i),w(ob3j),w(obl3i),w(obl3j),nlmak,nlmak,
     .        ndim3i,ndim3j,nlmak,w(opjj),w(opjj2),w(opjj3),w(opjjt),
     .        ldiag,w(owk2),w(owk))
            call hscop(ic,nclus3,ie,je,iax,w(oips),
     .        iax3(ic3),3,iofc3,iofh(ivh),nlm1,w(owk),nlm2,
     .        w(owk2),ndim3i,w(owk),nhs,h)
          endif
C ...     Three-center contribution to s from this kb
          if (ls .and. ndimi*ndimj > 0) then
            call xxsl3c(w(ob3i),w(ob3j),w(obl3i),w(obl3j),nlmak,nlmak,
     .        ndim3i,ndim3j,nlmak,nla,svecs(k0,1,ivs),xxs,
     .        svecs(k0,1,ivls),xx2,svecs(k0,1,nvls),xx3,
     .        svecs(k0,1,jvls),xxt,ldiag,w(owk2),w(owk))
            call hscop(ic,nclus3,ie,je,iax,w(oips),
     .        iax3(ic3),3,iofc3,iofh(ivh),nlm1,w(owk),nlm2,
     .        w(owk2),ndim3i,w(owk),nhs,s)
          endif
          call rlse(owk)
        endif
        call rlse(ob3j)
   30 continue
      call rlse(ob3i)
   20 continue
      call rlse(oadd1)
   10 continue
      call rlse(opkk)

C --- Copy into other triangle ---
      if (lh .and. lstore == 1) then
        do  41  j = 1, nhs
        do  41  i = 1, j
   41   h(j,i) = h(i,j)
      endif

      if (ls .and. lstore == 1) then
        do  42  j = 1, nhs
        do  42  i = 1, j
   42   s(j,i) = s(i,j)
      endif

C      if (lh) call prmx('h in hsmatr',h,nhs,nhs,nhs)
C      if (ls) call prmx('s in hsmatr',s,nhs,nhs,nhs)
C     rewind 81
C      write(81) nhs,nhs,11
C      call dpdump(s,nhs*nhs,-81)
C      close(81)
C      call rx('done')

      call tcx('hsmatr')
      end
      subroutine xxhmgl(nclus,hil,nlma,bld,add1,add2,iofc,iofcl,nlm1,h1)
C- Adds bdot contribution from linked basis for one site
C  NB:  alternative would be to extract bld for this e using xxmcsx
      implicit none
      integer nhs,nlm1,nclus,nlma,iofc(nclus),iofcl(nclus)
      double precision h1(nlm1,1),hil,bld(nlma,1),add1(1),add2(1),xx
      integer nlmj,iofj,ilfj,jb,ilm,jlm

      if (nlm1 <= 0) return
      do  20  jb = 1, nclus
        nlmj = iofc(jb+1)-iofc(jb)
        iofj = iofc(jb)-iofc(1)
        ilfj = iofcl(jb)-iofcl(1)
        do  20  jlm = 1, nlmj
        xx = add2(jlm+iofj)*hil
        do  20  ilm = 1, nlm1
          h1(ilm,jlm+iofj) = h1(ilm,jlm+iofj) +
     .                       add1(ilm)*xx*bld(ilm,jlm+ilfj)
   20 continue
      end

      subroutine xxsm2c(b1,b2,bd,nla1,nla2,nlb1,nlb2,
     .  nlm1,nlm2,iof2,nlma,ldiag,si,nlf,fvec,xxx,s1,s2)
C- Adds 1C and 2C ASA contribution to S for all strux around one site.
C  ldiag=1 if diagonal block
      implicit none
      integer nla1,nla2,nlb1,nlb2,ldiag,nlm1,nlm2,nlma,iof2,nlf
      double precision si,fvec(nlf,4),xxx
      double precision b1(nla1,nlb1),b2(nla2,nlb2),bd(nlma,1),
     .  s1(nlm1,nlb2),s2(nlm2,nlb1)
      integer i,j,ilm,itop

C --- fv2 b(e2) (2C terms <ib,e1 | jb,e2>, augmentation at ib) ---
      do  11  j = 1, nlb2
      do  11  ilm = 1, nlm1
   11 s1(ilm,j) = (fvec(ilm,2)+xxx)*b2(ilm,j)

C --- fv3 b(e1) (2C terms <ib,e2 | jb,e1>, augmentation at ib) ---
      do  31  j = 1, nlb1
      do  31  ilm = 1, nlm2
   31 s2(ilm,j) = (fvec(ilm,3)-xxx)*b1(ilm,j)

C --- Add fv(1) ---
      do  42  ilm = 1, min(nlm1,nlm2)
   42 s1(ilm,ilm+iof2) = s1(ilm,ilm+iof2) + fvec(ilm,1)

C --- Add bdot (global) term associated with this strux
      if (ldiag == 1) then
        do  50  j = 1, nlb1
        itop = min(nlm1,j)
        do  50  i = 1, itop
   50   s1(i,j) = s1(i,j) + si*bd(i,j)
      endif

      end
      subroutine xxsm3c(b1,b2,nla1,nla2,nlb1,nlb2,nlma,nlf,fvec,
     .  ldiag,wk,s3)
C- Adds 3C ASA contribution to H for all strux around one site.
C  ldiag=1 if diagonal block
      implicit none
      integer nla1,nla2,nlb1,nlb2,ldiag,nlma,nlf
      double precision b1(nla1,nlb1),b2(nla2,nlb2),
     .  s3(nlb1,nlb2),wk(nlb1,nlma),fvec(nlf,4)
      integer i,itop,j,klm,ncol,irow,mrow,nrow
      parameter (mrow=16)

C --- b1(+) * fv(4) * b2, using wk for b1(+) fv(4) ---
      do  10  klm = 1, nlma
      do  10  j = 1, nlb1
   10 wk(j,klm) = b1(klm,j)*fvec(klm,4)
      if (ldiag == 1) then
C#ifdef BLAS3
        call dpzero(s3,nlb1*nlb2)
        do  11  irow = 1, nlb2, mrow
          nrow = min(nlb2-irow+1,mrow)
          ncol = nlb2-irow+1
          call dgemm('N','N',nrow,ncol,nlma,1d0,wk(irow,1),nlb1,
     .      b2(1,irow),nla2,0d0,s3(irow,irow),nlb1)
   11   continue
C#elseC
C        do  11  j = 1, nlb2
C        itop = ldiag*j + (1-ldiag)*nlb1
C        do  12  i = 1, itop
C   12   h3(i,j) = 0
C        do  11  klm = 1, nlma
C          do  14  i = 1, itop
C   14     h3(i,j) = h3(i,j) + wk(i,klm)*b2(klm,j)
C   11   continue
C#endif
      else
        call dgemm('N','N',nlb1,nlb2,nlma,1d0,wk,nlb1,
     .    b2,nla2,0d0,s3,nlb1)
      endif

      end
      subroutine xxsl2c(b1,b2,bd,bl1,bl2,nla1,nla2,nlb1,nlb2,
     .  nlm1,nlm2,iof1,iof2,nlma,si,add1,add2,
     .  nlf,fvec,xxx,fvc2,xx2,fvc3,xx3,fvct,xxt,ldiag,s1,s2)
C- Adds 1C and 2C ASA contribution to S for all strux around one site.
      implicit none
      integer nla1,nla2,nlb1,nlb2,ldiag,nlm1,nlm2,nlma,iof1,iof2,nlf
      double precision si,add1(nlb1),add2(nlb2),
     .  fvec(nlf,4),fvc2(nlf,4),fvc3(nlf,4),fvct(nlf,4),
     .  xxx,xx2,xx3,xxt,b1(nla1,nlb1),b2(nla2,nlb2),bd(nlma,1),
     .  bl1(nlma,nlb1),bl2(nlma,nlb2),
     .  s1(nlm1,nlb2),s2(nlm2,nlb1),
     .  x1,x2
      integer i,itop,j,ilm1,ilm2,ilm,ila,jlm

C --- fv(2) b(e2) (2C terms <ib,e1+el | jb,e2>, augmentation at ib) ---
      do  21  ilm = 1, nlm1
      x1 = (fvec(ilm,2)+xxx + add1(ilm+iof1)*(fvct(ilm,3)-xxt))
      x2 = (fvc2(ilm,2)+xx2 + add1(ilm+iof1)*(fvc3(ilm,2)+xx3))
      do  21  j = 1, nlb2
   21 s1(ilm,j) = x1*b2(ilm,j) + x2*bl2(ilm,j)

C --- fv3 b(e1) (2C terms <ib,e2+el | jb,e1>, augmentation at ib) ---
      do  31  jlm = 1, nlm2
      x1 = (fvec(jlm,3)-xxx + add2(iof2+jlm)*(fvc2(jlm,3)-xx2))
      x2 = (fvct(jlm,2)+xxt + add2(iof2+jlm)*(fvc3(jlm,3)-xx3))
      do  31  j = 1, nlb1
   31 s2(jlm,j) = x1*b1(jlm,j) + x2*bl1(jlm,j)

C --- Add fv(1) ---
      do  41  ilm = 1, min(nlm1,nlm2)
   41 s1(ilm,ilm+iof2) = s1(ilm,ilm+iof2) + fvec(ilm,1)
     .    + fvc2(ilm,1)*add2(ilm+iof2) + add1(ilm+iof1)*fvct(ilm,1)
     .    + add1(ilm+iof1)*fvc3(ilm,1)*add2(ilm+iof2)

C --- Add bdot (global) term associated with this strux
      if (ldiag == 1) then
        do  50  j = 1, nlb1
        itop = min(nlm1,j)
        do  50  i = 1, itop
   50   s1(i,j) = s1(i,j) + si*bd(i,j)
      endif

      end
      subroutine xxsl3c(b1,b2,bl1,bl2,nla1,nla2,nlb1,nlb2,
     .  nlma,nlf,fvec,xxx,fvc2,xx2,fvc3,xx3,fvct,xxt,ldiag,wk,s3)
C- Adds 3C ASA contribution to S for all strux around one site.
C  ldiag=1 if diagonal block
      implicit none
      integer nla1,nla2,nlb1,nlb2,ldiag,nlma,nlf
      double precision fvec(nlf,4),fvc2(nlf,4),fvc3(nlf,4),fvct(nlf,4),
     .  xxx,xx2,xx3,xxt,b1(nla1,nlb1),b2(nla2,nlb2),
     .  bl1(nlma,nlb1),bl2(nlma,nlb2),s3(nlb1,nlb2),wk(nlb1,nlma)
      integer i,itop,j,klm,nrow,ncol,mrow,irow
      parameter (mrow=16)

C --- (fv4 * b2 + fvct bl)+ * b2 ---
      do  10  klm = 1, nlma
      do  10  j = 1, nlb1
   10 wk(j,klm) = b1(klm,j)*fvec(klm,4) + bl1(klm,j)*fvct(klm,4)
      if (ldiag == 1) then
C#ifdef BLAS3
        call dpzero(s3,nlb1*nlb2)
        do  11  irow = 1, nlb2, mrow
          nrow = min(nlb2-irow+1,mrow)
          ncol = nlb2-irow+1
          call dgemm('N','N',nrow,ncol,nlma,1d0,wk(irow,1),nlb1,
     .      b2(1,irow),nla2,0d0,s3(irow,irow),nlb1)
   11   continue
C#elseC
C        do  11  j = 1, nlb2
C        itop = ldiag*j + (1-ldiag)*nlb1
C        do  12  i = 1, itop
C   12   s3(i,j) = 0
C        do  11  klm = 1, nlma
C          do  14  i = 1, itop
C   14     s3(i,j) = s3(i,j) + wk(i,klm)*b2(klm,j)
C   11   continue
C#endif
      else
        call dgemm('N','N',nlb1,nlb2,nlma,1d0,wk,nlb1,
     .    b2,nla2,0d0,s3,nlb1)
      endif

C --- (fv_2(4) * b1 + fv_3(4) bl1)+ * bl2  ---
      do  18  klm = 1, nlma
      do  18  j = 1, nlb1
   18 wk(j,klm) = b1(klm,j)*fvc2(klm,4) + bl1(klm,j)*fvc3(klm,4)
      if (ldiag == 1) then
C#ifdef BLAS3
        do  15  irow = 1, nlb2, mrow
          nrow = min(nlb2-irow+1,mrow)
          ncol = nlb2-irow+1
          call dgemm('N','N',nrow,ncol,nlma,1d0,wk(irow,1),nlb1,
     .      bl2(1,irow),nlma,1d0,s3(irow,irow),nlb1)
   15   continue
C#elseC
C        do  15  j = 1, nlb2
C        itop = ldiag*j + (1-ldiag)*nlb1
C        do  15  klm = 1, nlma
C          do  16  i = 1, itop
C   16     s3(i,j) = s3(i,j) + wk(i,klm)*bl2(klm,j)
C   15   continue
C#endif
      else
        call dgemm('N','N',nlb1,nlb2,nlma,1d0,wk,nlb1,
     .    bl2,nlma,1d0,s3,nlb1)
      endif

      end
      subroutine xxhm2c(b1,b2,bd,nla1,nla2,nlb1,nlb2,
     .  nlm1,nlm2,iof2,nlma,hi,pkk,pkj,pjk,ldiag,h1,h2)
C- Adds 1C and 2C ASA contribution to H for all strux around one site.
C  ldiag=1 if diagonal block
      implicit none
      integer nla1,nla2,nlb1,nlb2,ldiag,nlm1,nlm2,nlma,iof2
      double precision hi,pkk(nlm1,nlm2),pkj(nlm1,nlma),pjk(nlma,nlm2)
      double precision b1(nla1,nlb1),b2(nla2,nlb2),bd(nlma,1),
     .  h1(nlm1,nlb2),h2(nlm2,nlb1)
      integer i,itop,j,ilm1,ilm2

C --- pkj b(e2) (2C terms <ib,e1 | jb,e2>, augmentation at ib) ---
      if (nlm1 > 0)
     .  call dgemm('N','N',nlm1,nlb2,nlma,1d0,pkj,nlm1,
     .  b2,nla2,0d0,h1,nlm1)

C --- pjk+ b(e1) (2C terms <ib,e2 | jb,e1>, augmentation at ib) ---
      if (nlm2 > 0)
     .  call dgemm('T','N',nlm2,nlb1,nlma,1d0,pjk,nlma,
     .  b1,nla1,0d0,h2,nlm2)

C --- Add pkk (upper triangle only if ldiag=1) ---
      do  42  ilm2 = 1, nlm2
        itop = ldiag*ilm2 + (1-ldiag)*nlm1
        do  41  ilm1 = 1, itop
   41   h1(ilm1,ilm2+iof2) = h1(ilm1,ilm2+iof2) + pkk(ilm1,ilm2)
   42 continue

C --- Add bdot (global) term associated with this strux
      if (ldiag == 1) then
        do  50  j = 1, nlb1
        itop = min(nlm1,j)
        do  50  i = 1, itop
   50   h1(i,j) = h1(i,j) + hi*bd(i,j)
      endif

      end

      subroutine xxhm3c(b1,b2,nla1,nla2,nlb1,nlb2,nlma,pjj,ldiag,wk,h3)
C- Adds 3C ASA contribution to H for all strux around one site.
C  ldiag=1 if diagonal block
      implicit none
      integer nla1,nla2,nlb1,nlb2,ldiag,nlma
      double precision pjj(nlma,nlma)
      double precision b1(nla1,nlb1),b2(nla2,nlb2),
     .  h3(nlb1,nlb2),wk(nlb1,nlma)
      integer i,itop,j,klm,ncol,irow,mrow,nrow
      parameter (mrow=16)

C --- b1(+) * pjj * b2, using wk for b1(+) pjj ---
      call dgemm('T','N',nlb1,nlma,nlma,1d0,b1,nla1,
     .  pjj,nlma,0d0,wk,nlb1)
      if (ldiag == 1) then
C#ifdef BLAS3
        call dpzero(h3,nlb1*nlb2)
        do  11  irow = 1, nlb2, mrow
          nrow = min(nlb2-irow+1,mrow)
          ncol = nlb2-irow+1
          call dgemm('N','N',nrow,ncol,nlma,1d0,wk(irow,1),nlb1,
     .      b2(1,irow),nla2,0d0,h3(irow,irow),nlb1)
   11   continue
C#elseC
C        do  11  j = 1, nlb2
C        itop = ldiag*j + (1-ldiag)*nlb1
C        do  12  i = 1, itop
C   12   h3(i,j) = 0
C        do  11  klm = 1, nlma
C          do  14  i = 1, itop
C   14     h3(i,j) = h3(i,j) + wk(i,klm)*b2(klm,j)
C   11   continue
C#endif
      else
        call dgemm('N','N',nlb1,nlb2,nlma,1d0,wk,nlb1,
     .    b2,nla2,0d0,h3,nlb1)
      endif

C      call prmx('pjj in xxhm3c',pjj,nlma,nlma,nlma)
C      call prmx('h3  in xxhm3c',h3,nlb1,nlb1,nlb2)

      end

      subroutine xxhl2c(b1,b2,bd,bl1,bl2,nla1,nla2,nlb1,nlb2,
     .  nlm1,nlm2,iof1,iof2,nlma,hi,add1,add2,
     .  pkk,pkj,pjk,pkk2,pkj2,pjk2,pkk3,pkj3,pjk3,
     .  pkkt,pkjt,pjkt,ldiag,h1,h2)
C- Adds 1C and 2C ASA contribution to H for all strux around one site.
C  ldiag=1 if diagonal block
      implicit none
      integer nla1,nla2,nlb1,nlb2,ldiag,
     .  nlm1,nlm2,nlma,iof1,iof2
      double precision hi,add1(nlb1),add2(nlb2),
     .  pkk (nlm1,nlm2),pkj (nlm1,nlma),pjk (nlma,nlm2),
     .  pkk2(nlm1,nlm2),pkj2(nlm1,nlma),pjk2(nlma,nlm2),
     .  pkkt(nlm2,nlm2),pkjt(nlm2,nlma),pjkt(nlma,nlm2),
     .  pkk3(nlm1,nlm2),pkj3(nlm1,nlma),pjk3(nlma,nlm2)
      double precision b1(nla1,nlb1),b2(nla2,nlb2),bd(nlma,1),
     .  bl1(nlma,nlb1),bl2(nlma,nlb2),
     .  h1(nlm1,nlb1),h2(nlm2,nlb2)
      integer i,itop,j,ilm1,ilm2,ilm,ila

C --- pkj b(e2) (2C terms <ib,e1+el | jb,e2>, augmentation at ib) ---
      if (nlm1 > 0) then
        do  21  ila = 1, nlma
        do  21  ilm = 1, nlm1
   21   pkj(ilm,ila) = pkj(ilm,ila) + add1(ilm+iof1)*pjkt(ila,ilm)
        call dgemm('N','N',nlm1,nlb2,nlma,1d0,pkj,nlm1,
     .    b2,nla2,0d0,h1,nlm1)
C         call prmx('h in xxhl2c',h1,nlm1,nlm1,nlb2)
        do  22  ila = 1, nlma
        do  22  ilm = 1, nlm1
   22   pkj3(ilm,ila) = add1(ilm+iof1)*pkj3(ilm,ila) + pkj2(ilm,ila)
        call dgemm('N','N',nlm1,nlb2,nlma,1d0,pkj3,nlm1,
     .    bl2,nlma,1d0,h1,nlm1)
C         call prmx('h in xxhl2c',h1,nlm1,nlm1,nlb2)
      endif

C --- pjk+ b(e1) (2C terms <ib,e2+el | jb,e1>, augmentation at ib) ---
      if (nlm2 > 0) then
        do  31  ilm = 1, nlm2
        do  31  ila = 1, nlma
   31   pjk(ila,ilm) = pjk(ila,ilm) + pjk2(ila,ilm)*add2(ilm+iof2)
        call dgemm('T','N',nlm2,nlb1,nlma,1d0,pjk,nlma,
     .    b1,nla1,0d0,h2,nlm2)
        do  33  ilm = 1, nlm2
        do  33  ila = 1, nlma
   33   pjk3(ila,ilm) = pjk3(ila,ilm)*add2(ilm+iof2) + pkjt(ilm,ila)
        call dgemm('T','N',nlm2,nlb1,nlma,1d0,pjk3,nlma,
     .    bl1,nlma,1d0,h2,nlm2)
      endif

C --- Add pkk (upper triangle only if ldiag=1) ---
      do  42  ilm2 = 1, nlm2
        itop = ldiag*ilm2 + (1-ldiag)*nlm1
        do  41  ilm1 = 1, itop
   41   h1(ilm1,ilm2+iof2) = h1(ilm1,ilm2+iof2) + pkk(ilm1,ilm2)
     .      + pkk2(ilm1,ilm2)*add2(ilm2+iof2)
     .      + add1(ilm1+iof1)*pkkt(ilm2,ilm1)
     .      + add1(ilm1+iof1)*pkk3(ilm1,ilm2)*add2(ilm2+iof2)
   42 continue

C --- Add bdot (global) term associated with this strux
      if (ldiag == 1) then
        do  50  j = 1, nlb1
        itop = min(nlm1,j)
        do  50  i = 1, itop
   50   h1(i,j) = h1(i,j) + hi*bd(i,j)
      endif

      end
      subroutine xxhl3c(b1,b2,bl1,bl2,nla1,nla2,nlb1,nlb2,
     .  nlma,pjj,pjj2,pjj3,pjjt,ldiag,wk,h3)
C- Adds 3C ASA contribution to H for all strux around one site.
C  ldiag=1 if diagonal block
      implicit none
      integer nla1,nla2,nlb1,nlb2,ldiag,nlma
      double precision
     .  pjj(nlma,nlma),pjj2(nlma,nlma),pjjt(nlma,nlma),pjj3(nlma,nlma)
      double precision b1(nla1,nlb1),b2(nla2,nlb2),
     .  bl1(nlma,nlb1),bl2(nlma,nlb2),
     .  h3(nlb1,nlb2),wk(nlb1,nlma)
      integer i,itop,j,klm,nrow,ncol,mrow,irow
      parameter (mrow=16)

C --- b1(+) * pjj * b2, using wk for b1(+) pjj ---
      call dgemm('T','N',nlb1,nlma,nlma,1d0,b1,nla1,
     .  pjj,nlma,0d0,wk,nlb1)
      call dgemm('T','T',nlb1,nlma,nlma,1d0,bl1,nlma,
     .  pjjt,nlma,1d0,wk,nlb1)
      if (ldiag == 1) then
C#ifdef BLAS3
        call dpzero(h3,nlb1*nlb2)
        do  11  irow = 1, nlb2, mrow
          nrow = min(nlb2-irow+1,mrow)
          ncol = nlb2-irow+1
          call dgemm('N','N',nrow,ncol,nlma,1d0,wk(irow,1),nlb1,
     .      b2(1,irow),nla2,0d0,h3(irow,irow),nlb1)
   11   continue
C#elseC
C        do  11  j = 1, nlb2
C        itop = ldiag*j + (1-ldiag)*nlb1
C        do  12  i = 1, itop
C   12   h3(i,j) = 0
C        do  11  klm = 1, nlma
C          do  14  i = 1, itop
C   14     h3(i,j) = h3(i,j) + wk(i,klm)*b2(klm,j)
C   11   continue
C#endif
      else
        call dgemm('N','N',nlb1,nlb2,nlma,1d0,wk,nlb1,
     .    b2,nla2,0d0,h3,nlb1)
      endif

      call dgemm('T','N',nlb1,nlma,nlma,1d0,b1,nla1,
     .  pjj2,nlma,0d0,wk,nlb1)
      call dgemm('T','N',nlb1,nlma,nlma,1d0,bl1,nlma,
     .  pjj3,nlma,1d0,wk,nlb1)
      if (ldiag == 1) then
C#ifdef BLAS3
        do  15  irow = 1, nlb2, mrow
          nrow = min(nlb2-irow+1,mrow)
          ncol = nlb2-irow+1
          call dgemm('N','N',nrow,ncol,nlma,1d0,wk(irow,1),nlb1,
     .      bl2(1,irow),nlma,1d0,h3(irow,irow),nlb1)
   15   continue
C#elseC
C        do  15  j = 1, nlb2
C        itop = ldiag*j + (1-ldiag)*nlb1
C        do  15  klm = 1, nlma
C          do  16  i = 1, itop
C   16     h3(i,j) = h3(i,j) + wk(i,klm)*bl2(klm,j)
C   15   continue
C#endif
      else
        call dgemm('N','N',nlb1,nlb2,nlma,1d0,wk,nlb1,
     .    bl2,nlma,1d0,h3,nlb1)
      endif

      end
