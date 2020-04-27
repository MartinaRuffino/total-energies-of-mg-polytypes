C#define BLAS3
      subroutine hmatl(el,nel,b,bd,bl,bld,ceadd,hvecs,hi,nbas,isp,
     .  ips,gh,gj,puu,pus,pss,ioffa,ioffb,ioffl,ioffp,nhs,nhl,nla,hbf,h)
C- Hamiltonian matrix, linked basis (real h).
C  Compressed bd,bld arrays.
      implicit none
      integer nbas,nel,nhs,nla,ips(1),isp,
     .  ioffp(1),ioffb(nbas,nel),ioffl(nbas+1),ioffa(1),nhl
      double precision hvecs(nla,4,1),hi(6),el(1),ceadd(25,5,1),
     .  gh(nla,2,nel),gj(nla,2,nel),puu(1),pus(1),pss(1)
      double precision b(nla,nhs),bl(nla,1),bd(1),bld(1),h(nhs,nhs),
     .  hbf(nhs,nhs),xxx
      integer i,ie,iv,ivs,j,je,ldiag,ne1,ne2,owk,koff,iof1,iof2,iofi,
     .  kb,nlm1,nlm2,nlma,nx,k0,io,je1,ivl,ivls,jvl,jvls,il,ilfi,
     .  nvl,nvls,oadd1,oadd2,obl1,obl2,ndimi,ndimj,lsp,nsp,
     .  opkk,opkj,opjk,opjj,opkk2,opkj2,opjk2,opjj2,
     .  opkkt,opkjt,opjkt,opjjt,opkk3,opkj3,opjk3,opjj3
      parameter (nx=49)
      integer , dimension(:), allocatable :: kkoff
      integer , dimension(:,:), allocatable :: iiv
      real w(1)
      common /w/ w

C For MPI ...
      integer, dimension(:), allocatable :: bproc
      integer mpipid,procid,master,numprocs,length,ierr
      integer ib1,ib2,iprint,lgunit
      logical mlog,cmdopt,MPI
      character outs*80, shortname*10, datim*26
      procid = mpipid(1)
      numprocs = mpipid(0)
      MPI = numprocs > 1
      call mpiprn(shortname,length)
      master = 0
      mlog = cmdopt('--mlog',6,0,outs)

      if (nhl /= 0 .and. MPI) then
        if (iprint() >= 30) then
          write (lgunit(1),1)
    1     format (' HMATL: MPI not checked for linked basis.'//
     .      ' Using serial implementation ...')
        endif
        MPI = .false.
      endif

      call tcn('hmatl')
      il = 0
      nsp = lsp()+1
      call defrr(opkk, nx*nx)
      call defrr(opkj, nx*nx)
      call defrr(opjk, nx*nx)
      call defrr(opjj, nx*nx)
      if (nhl > 0) then
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
        call defrr(oadd1,25*nbas)
        call defrr(oadd2,25*nbas)
        il = nel+1
        nvl = ((nel+1)*(nel+2))/2
        nvls= nsp*(nvl-1)+isp
      endif

C --- Pointer arrays ---
      allocate (kkoff(1:nbas))
      kkoff(1) = 0
      do  kb = 2, nbas
        kkoff(kb) = kkoff(kb-1) + ioffa(kb+1)-ioffa(kb)
      enddo
      allocate (iiv(1:nel,1:nel))
      iv = 0
      do  ie = 1, nel
        do  je = ie, nel
          iv = iv + 1
          iiv(ie,je) = iv
        enddo
        if (nhl /= 0) iv = iv + 1
      enddo

C --- buffer input hamiltonian and zero out ---
      if (MPI) then
        call dcopy(nhs*nhs,h,1,hbf,1)
        call dcopy(nhs*nhs,0d0,0,h,1)
      endif

C --- Loop over energy pairs and basis ---
      if (MPI) then
        allocate (bproc(0:numprocs), stat=ierr)
        call dstrbp(nbas,numprocs,1,bproc(0))
        ib1 = bproc(procid)
        ib2 = bproc(procid+1)-1
      else
        ib1 = 1
        ib2 = nbas
      endif
      do   kb = ib1, ib2
        if (MPI .and. mlog) then
        if (kb == bproc(procid)) then
          call gettime(datim)
          call awrit4(' hmatl '//datim//' Process %i of %i on '
     .       //shortname(1:length)//
     .       ' starting atoms %i to %i',' ',256,lgunit(3),
     .        procid,numprocs,bproc(procid),bproc(procid+1)-1)
        endif
        endif
        do  ie = 1, nel
c...sl
            ndimi = ioffb(nbas+1,ie)-ioffb(1,ie)
c         je1 = 1 + (ie-1)*ndimi**2
c...sl
          if (nhl > 0) then
            call hmcadd(ips,nbas,ioffb,ie,ceadd,w(oadd1))
          endif
          do  je = ie, nel
            ivl = (ie*(il+il-ie+1))/2
            ivls= nsp*(ivl-1)+isp
c           ndimi = ioffb(nbas+1,ie)-ioffb(1,ie)
            jvl = (je*(il+il-je+1))/2
            jvls= nsp*(jvl-1)+isp
            ndimj = ioffb(nbas+1,je)-ioffb(1,je)
            if (nhl > 0) then
              call hmcadd(ips,nbas,ioffb,je,ceadd,w(oadd2))
            endif
C       iv: index to (ie,je) pair;  ivs: index to (isp,ie,je)
            iv = iiv(ie,je)
            ivs = nsp*(iv-1)+isp
            ldiag = 0
            if (ie == je) ldiag = 1
            koff  = kkoff(kb)
            iof1 = ioffb(kb,ie)-ioffb(1,ie)
            iof2 = ioffb(kb,je)-ioffb(1,je)
            ne1  = 1+ioffb(1,ie)
            ne2  = 1+ioffb(1,je)
            nlm1 = ioffb(kb+1,ie)-ioffb(kb,ie)
            nlm2 = ioffb(kb+1,je)-ioffb(kb,je)
            nlma = ioffa(kb+1)-ioffa(kb)

C --- Make perturbation matrices ---
            if (nlma > nx) call rx('hmatl: nlma gt nx')
            io = nsp*ioffp(kb)+1
            if (isp == 2) io = io + nlma**2
            k0 = koff+1
            xxx = 0d0
            if (ldiag == 0) xxx = hi(iv)/(el(je)-el(ie))
            call pertmb(puu(io),pus(io),pss(io),hvecs(k0,1,ivs),nla,xxx,
     .                  gh(k0,1,ie),gh(k0,2,ie),gh(k0,1,je),gh(k0,2,je),
     .                  gj(k0,1,ie),gj(k0,2,ie),gj(k0,1,je),gj(k0,2,je),
     .                  nlm1,nlm2,nlma,w(opkk),w(opkj),w(opjk),w(opjj))
C --- Add to Hamiltonian this (ie,je,kb) block ---
            call defrr(owk,   nlma*max(ndimi,ndimj))
C ...   Case no linked basis
            if (nhl == 0) then
              if (ndimi > 0 .and. ndimj > 0) then
                call xxhmt5(b(1+koff,ne1),b(1+koff,ne2),nla,nla,
     .            ndimi,ndimj,nlm1,nlm2,iof1,iof2,nlma,
     .            w(opkk),w(opkj),w(opjk),w(opjj),w(owk),ldiag,nhs,
     .            h(ne1,ne2))
              endif
            else
C ...     pert matrices for (ie,linked) block
              xxx = hi(ivl)/(el(il)-el(ie))
              call pertmb(puu(io),pus(io),pss(io),hvecs(k0,1,ivls),
     .           nla,xxx,gh(k0,1,ie),gh(k0,2,ie),gh(k0,1,il),
     .           gh(k0,2,il),gj(k0,1,ie),gj(k0,2,ie),gj(k0,1,il),
     .           gj(k0,2,il),nlm1,nlm2,nlma,w(opkk2),w(opkj2),w(opjk2),
     .           w(opjj2))
C ...     pert matrices for (je,linked) block
              xxx = hi(jvl)/(el(il)-el(je))
              call pertmb(puu(io),pus(io),pss(io),hvecs(k0,1,jvls),
     .           nla,xxx,gh(k0,1,je),gh(k0,2,je),gh(k0,1,il),
     .           gh(k0,2,il),gj(k0,1,je),gj(k0,2,je),gj(k0,1,il),
     .           gj(k0,2,il),nlm2,nlm1,nlma,w(opkkt),w(opkjt),w(opjkt),
     .           w(opjjt))
C ...     pert matrices for (linked,linked) block
              xxx = 0d0
              call pertmb(puu(io),pus(io),pss(io),hvecs(k0,1,nvls),
     .           nla,xxx,gh(k0,1,il),gh(k0,2,il),gh(k0,1,il),
     .           gh(k0,2,il),gj(k0,1,il),gj(k0,2,il),gj(k0,1,il),
     .           gj(k0,2,il),nlm1,nlm2,nlma,w(opkk3),w(opkj3),w(opjk3),
     .           w(opjj3))
C ...     Scale and copy linked strux to work arrays
              call defrr(obl1,   nlma*ndimi)
              call defrr(obl2,   nlma*ndimj)
              call xxmcsx(nbas,nla,nlma,ioffb(1,ie),ioffl,bl(1+koff,1),
     .          w(obl1),w(oadd1))
              call xxmcsx(nbas,nla,nlma,ioffb(1,je),ioffl,bl(1+koff,1),
     .          w(obl2),w(oadd2))
              if (ndimi > 0 .and. ndimj > 0) then
                call xxhmtl(b(1+koff,ne1),b(1+koff,ne2),w(obl1),w(obl2),
     .                     nla,nla,ndimi,ndimj,nlm1,nlm2,iof1,iof2,nlma,
     .                     w(owk),w(oadd1),w(oadd2),ldiag,nhs,
     .                     w(opkk),w(opkj),w(opjk),w(opjj),
     .                     w(opkk2),w(opkj2),w(opjk2),w(opjj2),
     .                     w(opkk3),w(opkj3),w(opjk3),w(opjj3),
     .                     w(opkkt),w(opkjt),w(opjkt),w(opjjt),
     .                     h(ne1,ne2))
              endif
            endif
            call rlse(owk)

C --- Global (bdot) terms from linked basis ---
            if (nhl /= 0) then
              iofi = ioffb(kb,ie)-ioffb(1,ie)
              ilfi = ioffl(kb)-ioffl(1)
              call yyhmtl(nbas,iofi,ilfi,hi(nvl),nhl,bld,
     .          w(oadd1),w(oadd2),ioffb(1,je),ioffl,nlm1,nhs,h(ne1,ne2))
            endif
          enddo ! loop over ie
        enddo   ! loop over je
      enddo     ! loop over kb

      if (MPI) then
        call mpibc2(h,nhs*nhs,4,3,mlog,'hmatl','h')
C --- restore buffered components of h --
        ierr = mpipid(2)
        call daxpy(nhs*nhs,1d0,hbf,1,h,1)
        deallocate(bproc)
      endif

C --- Global (bdot) diagonal terms  ---
c...sl
      je1 = 1
      do  ie = 1, nel
c       je1 = 1 + (ie-1)*ndimi**2
        iv = iiv(ie,ie)
        ndimi = ioffb(nbas+1,ie)-ioffb(1,ie)
        ne1  = 1+ioffb(1,ie)
        call yyhmt5(hi(iv),ndimi,bd(je1),nhs,h(ne1,ne1))
        je1 = je1 + ndimi**2
      enddo
c...sl

      call rlse(opkk)

C --- Copy into other triangle ---
      do  j = 1, nhs
        do  i = 1, j
          h(j,i) = h(i,j)
        enddo
      enddo

C     call prmx('h in hmatl',h,nhs,nhs,nhs)
      deallocate(kkoff)
      deallocate(iiv)
      call tcx('hmatl')
      end
      subroutine xxmcsx(nbas,nla,nlma,ioffb,ioffl,bl,blc,alfa)
C- Contracts strux for linked energy expanded about one site
      implicit none
      integer nbas,nla,nlma,ioffb(nbas),ioffl(nbas)
      double precision blc(nlma,1),bl(nla,1),alfa(1)
      integer ib,ilma,ilfi,ilmb,iofi,nlmb

      do  10  ib = 1, nbas
        nlmb = ioffb(ib+1)-ioffb(ib)
        iofi = ioffb(ib)-ioffb(1)
        ilfi = ioffl(ib)-ioffl(1)
        do  20  ilmb = 1, nlmb
        do  20  ilma = 1, nlma
   20   blc(ilma,ilmb+iofi) = bl(ilma,ilmb+ilfi)*alfa(ilmb+iofi)
   10 continue
      end
      subroutine hmcadd(ips,nbas,ioffb,ie,ceadd,add)
C- Expands linked coefficient array for one energy
      implicit none
      integer ips(1),nbas,ioffb(nbas,1),ie
      double precision ceadd(25,5,1),add(1)
      integer kb,ioff,ks,ilm,nlm
      do  10  kb = 1, nbas
      ks = ips(kb)
      ioff = ioffb(kb,ie)-ioffb(1,ie)
      nlm = ioffb(kb+1,ie)-ioffb(kb,ie)
      do  10  ilm = 1, nlm
   10 add(ilm+ioff) = ceadd(ilm,ie,ks)
      end
      subroutine yyhmtl(nbas,iofi,ilfi,hil,nhl,bld,add1,add2,
     .  ioffb,ioffl,nlmi,nhs,h)
C- Adds bdot (global) terms from from linked basis
      implicit none
      integer nhs,nlmi,nbas,nhl,ioffb(nbas),ioffl(nbas),ilfi,iofi
      double precision h(nhs,nhs),hil,bld(nhl,nhl),add1(1),add2(1),xx
      integer nlmj,iofj,ilfj,jb,ilm,jlm

      if (nlmi <= 0) return
      do  20  jb = 1, nbas
      nlmj = ioffb(jb+1)-ioffb(jb)
      iofj = ioffb(jb)-ioffb(1)
      ilfj = ioffl(jb)-ioffl(1)
      do  20  jlm = 1, nlmj
      xx = add2(jlm+iofj)*hil
      do  20  ilm = 1, nlmi
   20 h(ilm+iofi,jlm+iofj) = h(ilm+iofi,jlm+iofj) +
     .    add1(ilm+iofi)*xx*bld(ilm+ilfi,jlm+ilfj)
      end
      subroutine yyhmt5(hi,nlb,bd,nh,h)
C- Adds bdot (global) terms
      implicit none
C ... Passed parameters
      integer nlb,nh
      double precision hi,h(nh,1),bd(nlb,nlb)
C ... Local parameters
      integer i,j
      if (nlb <= 0) return
      do  37  i = 1, nlb
      do  37  j = i, nlb
   37 h(i,j) = h(i,j) + hi*bd(i,j)
      end
      subroutine xxhmt5(b1,b2,nla1,nla2,nlb1,nlb2,
     .  nlm1,nlm2,iof1,iof2,nlma,pkk,pkj,pjk,pjj,wk,ldiag,nhs,h)
C- Adds contribution of one atom to an (e1,e2)-block of hamiltonian.
C  ldiag=1 if diagonal block
      implicit none
      integer nla1,nla2,nlb1,nlb2,nhs,ldiag,
     .  nlm1,nlm2,nlma,iof1,iof2
      double precision pkk(nlm1,nlm2),pkj(nlm1,nlma),
     .  pjk(nlma,nlm2),pjj(nlma,nlma)
      double precision b1(nla1,nlb1),b2(nla2,nlb2),
     .  h(nhs,1),wk(nlb1,1)
      integer itop,ilm1,ilm2,ncol,irow,mrow,nrow
      parameter (mrow=32)

C --- b1(+) * pjj * b2, using wk for b1(+) pjj ---
      call dgemm('T','N',nlb1,nlma,nlma,1d0,b1,nla1,
     .  pjj,nlma,0d0,wk,nlb1)
      if (ldiag == 1) then
C#ifdef BLAS3
        do  11  irow = 1, nlb2, mrow
          nrow = min(nlb2-irow+1,mrow)
          ncol = nlb2-irow+1
          call dgemm('N','N',nrow,ncol,nlma,1d0,wk(irow,1),nlb1,
     .      b2(1,irow),nla2,1d0,h(irow,irow),nhs)
   11   continue
C#elseC
C        do  11  j = 1, nlb2
C        itop = ldiag*j + (1-ldiag)*nlb1
C        do  11  klm = 1, nlma
C          do  14  i = 1, itop
C   14     h(i,j) = h(i,j) + wk(i,klm)*b2(klm,j)
C   11   continue
C#endif
      else
        call dgemm('N','N',nlb1,nlb2,nlma,1d0,wk,nlb1,
     .    b2,nla2,1d0,h,nhs)
      endif

C --- Add  pkj * b2 ---
      if (nlm1 > 0)
     .  call dgemm('N','N',nlm1,nlb2,nlma,1d0,pkj,nlm1,
     .  b2,nla2,1d0,h(1+iof1,1),nhs)

C --- Add  b1(+) pjk ---
      if (nlm2 > 0)
     .  call dgemm('T','N',nlb1,nlm2,nlma,1d0,b1,nla1,
     .  pjk,nlma,1d0,h(1,1+iof2),nhs)

C --- Add pkk ---
      do  42  ilm2 = 1, nlm2
        itop = ldiag*ilm2 + (1-ldiag)*nlm1
        do  41  ilm1 = 1, itop
   41   h(ilm1+iof1,ilm2+iof2) = h(ilm1+iof1,ilm2+iof2) + pkk(ilm1,ilm2)
   42 continue

      end
      subroutine xxhmtl(b1,b2,bl1,bl2,nla1,nla2,nlb1,nlb2,
     .  nlm1,nlm2,iof1,iof2,nlma,wk,add1,add2,ldiag,nhs,
     .  pkk,pkj,pjk,pjj,pkk2,pkj2,pjk2,pjj2,pkk3,pkj3,pjk3,pjj3,
     .  pkkt,pkjt,pjkt,pjjt,h)
C- Contribution of one atom to an (e1,e2)-block of linked hamiltonian
C  from strux expanded at one site.  ldiag=1 if diagonal block.
      implicit none
      integer nla1,nla2,nlb1,nlb2,nhs,ldiag,nlm1,nlm2,nlma,iof1,iof2
      double precision add1(nlb1),add2(nlb2),
     .  pkk (nlm1,nlm2),pkj (nlm1,nlma),pjk (nlma,nlm2),pjj (nlma,nlma),
     .  pkk2(nlm1,nlm2),pkj2(nlm1,nlma),pjk2(nlma,nlm2),pjj2(nlma,nlma),
     .  pkkt(nlm2,nlm2),pkjt(nlm2,nlma),pjkt(nlma,nlm2),pjjt(nlma,nlma),
     .  pkk3(nlm1,nlm2),pkj3(nlm1,nlma),pjk3(nlma,nlm2),pjj3(nlma,nlma)
      double precision b1(nla1,nlb1), b2(nla2,nlb2),h(nhs,1),
     .                 bl1(nlma,nlb1),bl2(nlma,nlb2),wk(nlb1,1)
      integer itop,ilm1,ilm2,ilm,ila,ncol,irow,mrow,nrow
      parameter (mrow=32)

C --- b1(+) * pjj * b2, using wk for b1(+) pjj ---
      call dgemm('T','N',nlb1,nlma,nlma,1d0,b1,nla1,
     .  pjj,nlma,0d0,wk,nlb1)
      call dgemm('T','T',nlb1,nlma,nlma,1d0,bl1,nlma,
     .  pjjt,nlma,1d0,wk,nlb1)
      if (ldiag == 1) then
C#ifdef BLAS3
        do  11  irow = 1, nlb2, mrow
          nrow = min(nlb2-irow+1,mrow)
          ncol = nlb2-irow+1
          call dgemm('N','N',nrow,ncol,nlma,1d0,wk(irow,1),nlb1,
     .      b2(1,irow),nla2,1d0,h(irow,irow),nhs)
   11   continue
C#elseC
C        do  11  j = 1, nlb2
C        itop = ldiag*j + (1-ldiag)*nlb1
C        do  11  klm = 1, nlma
C          do  12  i = 1, itop
C   12     h(i,j) = h(i,j) + wk(i,klm)*b2(klm,j)
C   11   continue
C#endif
      else
        call dgemm('N','N',nlb1,nlb2,nlma,1d0,wk,nlb1,
     .    b2,nla2,1d0,h,nhs)
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
     .      bl2(1,irow),nlma,1d0,h(irow,irow),nhs)
   15   continue
C#elseC
C        do  15  j = 1, nlb2
C          itop = ldiag*j + (1-ldiag)*nlb1
C          do  15  klm = 1, nlma
C          do  16  i = 1, itop
C   16     h(i,j) = h(i,j) + wk(i,klm)*bl2(klm,j)
C   15   continue
C#endif
      else
        call dgemm('N','N',nlb1,nlb2,nlma,1d0,wk,nlb1,
     .    bl2,nlma,1d0,h,nhs)
      endif

C --- pkj * b(e2), (2C terms <ib,e1+el | jb,e2>, augmentation at ib) ---
      if (nlm1 > 0) then
        do  21  ila = 1, nlma
        do  21  ilm = 1, nlm1
   21   pkj(ilm,ila) = pkj(ilm,ila) + add1(ilm+iof1)*pjkt(ila,ilm)
        call dgemm('N','N',nlm1,nlb2,nlma,1d0,pkj,nlm1,
     .    b2,nla2,1d0,h(1+iof1,1),nhs)
        do  22  ila = 1, nlma
        do  22  ilm = 1, nlm1
   22   pkj3(ilm,ila) = add1(ilm+iof1)*pkj3(ilm,ila) + pkj2(ilm,ila)
C       call prmx('h in xxhmtl',h,nhs,nlm1,nlb2)
        call dgemm('N','N',nlm1,nlb2,nlma,1d0,pkj3,nlm1,
     .    bl2,nlma,1d0,h(1+iof1,1),nhs)
C       call prmx('h in xxhmtl',h,nhs,nlm1,nlb2)
      endif

C --- b+(e2) pjk (2C terms <ib,e2+el | jb,e1>, augmentation at ib) ---
      if (nlm2 > 0) then
        do  31  ilm = 1, nlm2
        do  31  ila = 1, nlma
   31   pjk(ila,ilm) = pjk(ila,ilm) + pjk2(ila,ilm)*add2(ilm+iof2)
        call dgemm('T','N',nlb1,nlm2,nlma,1d0,b1,nla1,
     .    pjk,nlma,1d0,h(1,1+iof2),nhs)
        do  33  ilm = 1, nlm2
        do  33  ila = 1, nlma
   33   pjk3(ila,ilm) = pjk3(ila,ilm)*add2(ilm+iof2) + pkjt(ilm,ila)
        call dgemm('T','N',nlb1,nlm2,nlma,1d0,bl1,nlma,
     .    pjk3,nlma,1d0,h(1,1+iof2),nhs)
      endif

C --- Add pkk ---
      do  42  ilm2 = 1, nlm2
        itop = ldiag*ilm2 + (1-ldiag)*nlm1
        do  41  ilm1 = 1, itop
   41   h(ilm1+iof1,ilm2+iof2) = h(ilm1+iof1,ilm2+iof2) + pkk(ilm1,ilm2)
     .      + pkk2(ilm1,ilm2)*add2(ilm2+iof2)
     .      + add1(ilm1+iof1)*pkkt(ilm2,ilm1)
     .      + add1(ilm1+iof1)*pkk3(ilm1,ilm2)*add2(ilm2+iof2)
   42 continue

      end
