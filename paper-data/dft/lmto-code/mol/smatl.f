C#define BLAS3
      subroutine smatl(el,nel,b,bd,bl,bld,ceadd,svecs,si,nbas,isp,
     .  ips,ioffa,ioffb,ioffl,nhs,nhl,nla,s)
C- Overlap matrix, linked basis (real s).
C  Compressed bd,bld arrays.
C  5 Jan makes s for spin isp (MvS)
      implicit none
      integer nbas,isp,nel,nhs,nla,ips(1),ioffb(nbas,nel),
     .  ioffa(1),ioffl(nbas+1),nhl
      double precision b(nla,nhs),bl(nla,1),bd(1),bld(1),s(nhs,nhs),
     .  svecs(nla,4,1),si(6),el(1),ceadd(25,5,1),xxx,xx2,xx3,xxt
      integer i,ie,iv,ivs,j,je,ldiag,ne1,ne2,owk,koff,iof1,iof2,iofi,
     .  kb,ks,nlm1,nlm2,nlma,k0,je1,ivl,jvl,ivls,jvls,il,ilfi,
     .  nvl,nvls,oadd1,oadd2,obl1,obl2,ndimi,ndimj,nsp,lsp,lgunit
      integer , dimension(:), allocatable :: kkoff
      integer , dimension(:,:), allocatable :: iiv
      real w(1)
      common /w/ w

C For MPI ...
      integer, dimension(:), allocatable :: bproc
      integer mpipid,procid,master,numprocs,length,ierr
      integer ib1,ib2,iprint
      logical mlog,cmdopt,MPI
      character outs*80, shortname*10, datim*26
      procid = mpipid(1)
      numprocs = mpipid(0)
      MPI = numprocs > 1
      call mpiprn(shortname,length)
      master = 0
      mlog = cmdopt('--mlog',6,0,outs)

      call tcn('smatl')
      if (nhl /= 0 .and. MPI) then
        if (iprint() >= 30) then
          write (lgunit(1),1)
    1     format (' SMATL: MPI not checked for linked basis.'//
     .      ' Using serial implementation ...')
        endif
        MPI = .false.
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

      il = 0
      nsp = lsp()+1
      if (nhl > 0) then
        call defrr(oadd1,25*nbas)
        call defrr(oadd2,25*nbas)
        il = nel+1
        nvl = ((nel+1)*(nel+2))/2
        nvls= nsp*(nvl-1)+isp
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
          call awrit4(' smatl '//datim//' Process %i of %i on '
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
            koff  = kkoff(kb)
            ivl   = (ie*(il+il-ie+1))/2
            ivls  = nsp*(ivl-1)+isp
c           ndimi = ioffb(nbas+1,ie)-ioffb(1,ie)
            jvl   = (je*(il+il-je+1))/2
            jvls  = nsp*(jvl-1)+isp
            ndimj = ioffb(nbas+1,je)-ioffb(1,je)

            if (nhl > 0) then
              call hmcadd(ips,nbas,ioffb,je,ceadd,w(oadd2))
            endif
C       iv: index to (ie,je) pair;  ivs: index to (isp,ie,je)
            iv = iiv(ie,je)
            ivs=nsp*(iv-1)+isp
            ldiag = 0
            if (ie == je) ldiag = 1

            ks = ips(kb)
            iof1 = ioffb(kb,ie)-ioffb(1,ie)
            iof2 = ioffb(kb,je)-ioffb(1,je)
            ne1  = 1+ioffb(1,ie)
            ne2  = 1+ioffb(1,je)
            nlm1 = ioffb(kb+1,ie)-ioffb(kb,ie)
            nlm2 = ioffb(kb+1,je)-ioffb(kb,je)
            nlma = ioffa(kb+1)-ioffa(kb)
            k0 = koff+1
            xxx = 0d0
            if (ldiag == 0) xxx = si(iv)/(el(je)-el(ie))

C --- Add to Hamiltonian this (ie,je,kb) block ---
            call defrr(owk,   nlma*max(ndimi,ndimj))
C ...   Case no linked basis
            if (nhl == 0) then
              if (ndimi > 0 .and. ndimj > 0)
     .          call xxsmt5(b(1+koff,ne1),b(1+koff,ne2),nla,nla,
     .          ndimi,ndimj,nlm1,nlm2,iof1,iof2,nlma,
     .          w(owk),ldiag,nhs,nla,svecs(k0,1,ivs),xxx,s(ne1,ne2))
            else
              xx2 = si(ivl)/(el(il)-el(ie))
              xxt = si(jvl)/(el(il)-el(je))
              xx3 = 0d0
C ...     Scale and copy linked strux to work arrays
              call defrr(obl1,   nlma*ndimi)
              call defrr(obl2,   nlma*ndimj)
              call xxmcsx(nbas,nla,nlma,ioffb(1,ie),ioffl,bl(1+koff,1),
     .          w(obl1),w(oadd1))
              call xxmcsx(nbas,nla,nlma,ioffb(1,je),ioffl,bl(1+koff,1),
     .          w(obl2),w(oadd2))
              if (ndimi > 0 .and. ndimj > 0) then
                call xxsmtl(b(1+koff,ne1),b(1+koff,ne2),w(obl1),w(obl2),
     .            nla,nla,ndimi,ndimj,nlm1,nlm2,iof1,iof2,nlma,
     .            w(owk),w(oadd1),w(oadd2),ldiag,nhs,nla,
     .            svecs(k0,1,ivs),xxx,svecs(k0,1,ivls),xx2,
     .            svecs(k0,1,nvls),xx3,svecs(k0,1,jvls),xxt,s(ne1,ne2))
              endif
            endif

            call rlse(owk)

C --- Global (bdot) terms from linked basis ---
            if (nhl /= 0) then
              iofi = ioffb(kb,ie)-ioffb(1,ie)
              ilfi = ioffl(kb)-ioffl(1)
              call yyhmtl(nbas,iofi,ilfi,si(nvl),nhl,bld,
     .          w(oadd1),w(oadd2),ioffb(1,je),ioffl,nlm1,nhs,s(ne1,ne2))
            endif
          enddo ! loop over ie
        enddo   ! loop over je
      enddo     ! loop over kb

      if (MPI) then
        call mpibc2(s,nhs*nhs,4,3,mlog,'smatl','s')
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
        call yyhmt5(si(iv),ndimi,bd(je1),nhs,s(ne1,ne1))
        je1 = je1 + ndimi**2
      enddo
c...sl

C --- Copy into other triangle ---
      do   j = 1, nhs
        do   i = 1, j
          s(j,i) = s(i,j)
        enddo
      enddo

      if (nhl > 0) call rlse(oadd1)
C     call prmx('s in smatl',s,nhs,nhs,nhs)
      deallocate(kkoff)
      deallocate(iiv)

      call tcx('smatl')
      end
      subroutine yysmtl(nbas,iofi,ilfi,sil,nhl,bld,add1,add2,
     .  ioffb,ioffl,nlmi,nhs,s)
C- Adds bdot (global) terms from from linked basis
      implicit none
      integer nhs,nlmi,nbas,nhl,ioffb(nbas),ioffl(nbas),ilfi,iofi
      double precision s(nhs,nhs),sil,bld(nhl,nhl),add1(1),add2(1),xx
      integer nlmj,iofj,ilfj,jb,ilm,jlm

      if (nlmi <= 0) return
      do  20  jb = 1, nbas
      nlmj = ioffb(jb+1)-ioffb(jb)
      iofj = ioffb(jb)-ioffb(1)
      ilfj = ioffl(jb)-ioffl(1)
      do  20  jlm = 1, nlmj
      xx = add2(jlm+iofj)*sil
      do  20  ilm = 1, nlmi
   20 s(ilm+iofi,jlm+iofj) = s(ilm+iofi,jlm+iofj) +
     .    add1(ilm+iofi)*xx*bld(ilm+ilfi,jlm+ilfj)
      end
      subroutine xxsmt5(b1,b2,nla1,nla2,nlb1,nlb2,nlm1,nlm2,iof1,iof2,
     .  nlma,wk,ldiag,nhs,nlf,fvec,xxx,s)
C- Contribution of one atom to an (e1,e2)-block of overlap matrix
C  ldiag=1 if diagonal block.  nla,nlf dimension b and fvec
      implicit none
      integer nla1,nla2,nlb1,nlb2,nhs,ldiag,nlm1,nlm2,nlma,iof1,iof2,nlf
      double precision xxx,fvec(nlf,4)
      double precision b1(nla1,nlb1), b2(nla2,nlb2),s(nhs,1),wk(nlb1,1)
      integer i,itop,j,klm,ilm,jlm,jbot,ncol,irow,mrow,nrow
      parameter (mrow=32)

C --- b1(+) * fv(4) * b2, using wk for b1(+) fv(4) ---
      do  10  klm = 1, nlma
      do  10  j = 1, nlb1
   10 wk(j,klm) = b1(klm,j)*fvec(klm,4)
      if (ldiag == 1) then
C#ifdef BLAS3
        do  11  irow = 1, nlb2, mrow
          nrow = min(nlb2-irow+1,mrow)
          ncol = nlb2-irow+1
          call dgemm('N','N',nrow,ncol,nlma,1d0,wk(irow,1),nlb1,
     .      b2(1,irow),nla2,1d0,s(irow,irow),nhs)
   11   continue
C#elseC
C        do  11  j = 1, nlb2
C          itop = ldiag*j + (1-ldiag)*nlb1
C          do  11  klm = 1, nlma
C            do  12  i = 1, itop
C   12       s(i,j) = s(i,j) + wk(i,klm)*b2(klm,j)
C   11   continue
C#endif
      else
        call dgemm('N','N',nlb1,nlb2,nlma,1d0,wk,nlb1,
     .    b2,nla2,1d0,s,nhs)
      endif

C --- fv(2) b(e2) (2C terms <ib,e1 | jb,e2>, augmentation at ib) ---
      do  21  ilm = 1, nlm1
      i = iof1+ilm
      jbot = ldiag*i + (1-ldiag)
      do  21  j = jbot, nlb2
   21 s(ilm+iof1,j) = s(ilm+iof1,j) + (fvec(ilm,2)+xxx)*b2(ilm,j)

C --- fv(3) b(e1) (2C terms <ib,e2 | jb,e1>+, augmentation at ib) ---
      do  31  jlm = 1, nlm2
      j = iof2+jlm
      itop = ldiag*j + (1-ldiag)*nlb1
      do  31  i = 1, itop
   31 s(i,jlm+iof2) = s(i,jlm+iof2) + (fvec(jlm,3)-xxx)*b1(jlm,i)

C --- Add fv(1) ---
      do  41  ilm = 1, min(nlm1,nlm2)
   41 s(ilm+iof1,ilm+iof2) = s(ilm+iof1,ilm+iof2) + fvec(ilm,1)

      end
      subroutine xxsmtl(b1,b2,bl1,bl2,nla1,nla2,nlb1,nlb2,
     .  nlm1,nlm2,iof1,iof2,nlma,wk,add1,add2,ldiag,nhs,nlf,
     .  fvec,xxx,fvc2,xx2,fvc3,xx3,fvct,xxt,s)
C- Contribution of one atom to an (e1,e2)-block of linked overlap
C  ldiag=1 if diagonal block
      implicit none
      integer nla1,nla2,nlb1,nlb2,nhs,ldiag,nlm1,nlm2,nlma,iof1,iof2,nlf
      double precision add1(nlb1),add2(nlb2),xxx,xx2,xx3,xxt,
     .  fvec(nlf,4),fvc2(nlf,4),fvc3(nlf,4),fvct(nlf,4)
      double precision b1(nla1,nlb1), b2(nla2,nlb2),s(nhs,1),
     .                 bl1(nlma,nlb1),bl2(nlma,nlb2),wk(nlb1,1)
      integer i,itop,j,klm,ilm,jlm,jbot,ncol,irow,mrow,nrow
      parameter (mrow=32)

C --- b1(+) * fv(4) * b2, using wk for b1(+) fv(4) ---
      do  10  klm = 1, nlma
      do  10  j = 1, nlb1
   10 wk(j,klm) = b1(klm,j)*fvec(klm,4) + bl1(klm,j)*fvct(klm,4)
      if (ldiag == 1) then
C#ifdef BLAS3
        do  11  irow = 1, nlb2, mrow
          nrow = min(nlb2-irow+1,mrow)
          ncol = nlb2-irow+1
          call dgemm('N','N',nrow,ncol,nlma,1d0,wk(irow,1),nlb1,
     .      b2(1,irow),nla2,1d0,s(irow,irow),nhs)
   11   continue
C#elseC
C        do  11  j = 1, nlb2
C          itop = ldiag*j + (1-ldiag)*nlb1
C          do  11  klm = 1, nlma
C            do  12  i = 1, itop
C   12       s(i,j) = s(i,j) + wk(i,klm)*b2(klm,j)
C   11   continue
C#endif
      else
        call dgemm('N','N',nlb1,nlb2,nlma,1d0,wk,nlb1,
     .    b2,nla2,1d0,s,nhs)
      endif

      do  18  klm = 1, nlma
      do  18  j = 1, nlb1
   18 wk(j,klm) = b1(klm,j)*fvc2(klm,4) + bl1(klm,j)*fvc3(klm,4)
      if (ldiag == 1) then
C#ifdef BLAS3
        do  15  irow = 1, nlb2, mrow
          nrow = min(nlb2-irow+1,mrow)
          ncol = nlb2-irow+1
          call dgemm('N','N',nrow,ncol,nlma,1d0,wk(irow,1),nlb1,
     .      bl2(1,irow),nlma,1d0,s(irow,irow),nhs)
   15   continue
C#elseC
C        do  15  j = 1, nlb2
C          itop = ldiag*j + (1-ldiag)*nlb1
C          do  15  klm = 1, nlma
C            do  16  i = 1, itop
C   16       s(i,j) = s(i,j) + wk(i,klm)*bl2(klm,j)
C   15   continue
C#endif
      else
       call dgemm('N','N',nlb1,nlb2,nlma,1d0,wk,nlb1,bl2,nlma,1d0,s,nhs)
      endif

C --- fv(2) b(e2) (2C terms <ib,e1 | jb,e2>, augmentation at ib) ---
      do  21  ilm = 1, nlm1
      i = iof1+ilm
      jbot = ldiag*i + (1-ldiag)
      do  21  j = jbot, nlb2
   21 s(ilm+iof1,j) = s(ilm+iof1,j)
     . + (fvec(ilm,2)+xxx + add1(ilm+iof1)*(fvct(ilm,3)-xxt))*b2(ilm,j)
     . + (fvc2(ilm,2)+xx2 + add1(ilm+iof1)*(fvc3(ilm,2)+xx3))*bl2(ilm,j)

C --- fv(3) b(e1) (2C terms <ib,e2 | jb,e1>+, augmentation at ib) ---
      do  31  jlm = 1, nlm2
      j = iof2+jlm
      itop = ldiag*j + (1-ldiag)*nlb1
      do  31  i = 1, itop
   31 s(i,jlm+iof2) = s(i,jlm+iof2)
     . + (fvec(jlm,3)-xxx + add2(jlm+iof2)*(fvc2(jlm,3)-xx2))*b1(jlm,i)
     . + (fvct(jlm,2)+xxt + add2(jlm+iof2)*(fvc3(jlm,3)-xx3))*bl1(jlm,i)

C --- Add fv(1) ---
        do  41  ilm = 1, min(nlm1,nlm2)
   41   s(ilm+iof1,ilm+iof2) = s(ilm+iof1,ilm+iof2) + fvec(ilm,1)
     .      + fvc2(ilm,1)*add2(ilm+iof2) + add1(ilm+iof1)*fvct(ilm,1)
     .      + add1(ilm+iof1)*fvc3(ilm,1)*add2(ilm+iof2)
   42   continue

      end
