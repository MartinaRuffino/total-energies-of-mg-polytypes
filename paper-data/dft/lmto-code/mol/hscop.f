      subroutine hscop(ic,nclus,ie,je,iax,ipr,iax3,
     .  mode,ioffc,iofh,nlm1,h2ci,nlm2,h2cj,ndim1,h3c,nhs,h)
C- Add local cluster hamiltonian into the hamiltonian matrix
C  Inputs:
C    ic,nclus: index to start of this cluster in iax, and cluster size
C    ie,je and various tables of pointers:
C    ioffc: table of offsets to local cluster hamiltonian
C    iax (cf hpair.f) is a table information about cluster pairs
C    iax3 (cf hpair3.f) contains which elts in iax belong to 3C cluster
C    ipr (cf pair3.f) holds pos'n in iax for each pair in 3C cluster
C    iofh (cf hoffc.f) holds offset in relative to start of h for each
C                     (ie,je) (compound index) and cluster pair.
C    mode, one's digit: 2, copy 2C h; 3, copy 3C h
C          ten's digit: 1, copy h in unpacked form (molecules only)
C                       2, copy h in packed form.
C    h2ci: 2c h connecting <ib,ie|jb,je> (augmentation at ib) (mode=2)
C    h2cj: 2c h connecting <ib,je|jb,ie> (augmentation at ib) (mode=2)
C    h3c:  3c h connecting all pairs, with augmentation at ic (mode=3)
C  Output:  local cluster 2C h2ci,h2cj or 3C h3c is added into h
C  Remarks:
C    hscop copies the 2- and 3- center parts of a local cluster into a
C    hamiltonian or overlap matrix h, with minimal assumptions about the
C    structure of h.  Much of the needed information is passed in the
C    pointer arrays iax,iax3,iofh, and (for the 3-center terms) ipr.
C    Only two modes are allowed here (though minimal changes are needed
C    to make it completely general), but for now that seems like
C    overkill.  One mode is the normal molecule hamiltonian.  The other
C    stores h as a collection of (nlmi,nlmj) blocks, whose offset is
C    given by iofh.
      implicit none
      integer ic,nclus,nhc,nhs,mode,ie,je,iax(10,9),nlm1,nlm2,
     .  ndim1,ioffc(nclus,1),iax3(nclus),iofh(ic+nclus),ipr(nclus,nclus)
      double precision h(1),h2ci(nlm1,1),h2cj(nlm2,1),h3c(ndim1,1)
      integer ofhi,ofhj,ofci,ofcj,i,ib,j,jb,l1,l2,nlmi,nlmj,jbot,is,
     .  ofh,i1mach
      logical lxx

C --- Accumulate 2C terms to Hamiltonian ---
      if (mod(mode,10) == 2) then
C   ... Accumulate h2ci = h(ib,ie,jb,je) into h
        nlmi = ioffc(2,ie)-ioffc(1,ie)
Ch      ib = iax(1,ic)
Ch      nlmi = ioffb(ib+1,ie)-ioffb(ib,ie)
        do  10  j = 1, nclus
        nlmj = ioffc(j+1,je)-ioffc(j,je)
Ch      jb = iax(2,ic-1+j)
Ch      nlmj = ioffb(jb+1,je)-ioffb(jb,je)
        ofcj = ioffc(j,je)-ioffc(1,je)
C ...   skip blocks in lower triangle
        lxx = (iofh(ic-1+j) >= iofh(iax(6,ic-1+j))) .or. ie /= je
Ch      ofhi = ioffb(ib,ie)
Ch      ofhj = ioffb(jb,je)
Ch      lxx = (ofhi <= ofhj)
Ch      if (lxx .xor. (ofhi <= ofhj)) call rx('problem here')

c       if (iofh(ic-1+j) == 0)
c    .    call awrit6('2c j= %1,2i  ie,je= %i %i  ofcj%1,5i is=%i '//
c    .  '%,12;12d',' ',80,i1mach(2), j,ie,je,ofcj,ic-1+j,h2ci(1,1+ofcj))

        if (lxx) then
          do  12  l2 = 1, nlmj
Ch        ofh = (l2-1+ofhj)*nhs + ofhi
          ofh = (l2-1)*nhs + iofh(ic-1+j)
          do  12  l1 = 1, nlmi
   12     h(l1+ofh) = h(l1+ofh) + h2ci(l1,l2+ofcj)
        endif
   10   continue
C   ... Accumulate h2cj = h(ib,je,jb,ie) into h+(jb,ie,ib,je)
Ch      jb = ib
        nlmj = ioffc(2,je)-ioffc(1,je)
Ch      nlmj = ioffb(jb+1,je)-ioffb(jb,je)
        do  20  i = 1, nclus
        j = iax(6,ic-1+i)
        nlmi = ioffc(i+1,ie)-ioffc(i,ie)
Ch      ib = iax(2,ic-1+i)
Ch      nlmi = ioffb(ib+1,ie)-ioffb(ib,ie)
        ofci = ioffc(i,ie)-ioffc(1,ie)
C   ... skip blocks in lower triangle
        lxx = (iofh(j) >= iofh(iax(6,j))) .or. ie /= je
Ch      ofhi = ioffb(ib,ie)
Ch      ofhj = ioffb(jb,je)
Ch      lxx = (ofhi <= ofhj)
Ch      if (lxx .xor. (ofhi <= ofhj)) call rx('problem here')
        if (lxx) then
          do  22  l2 = 1, nlmj
          ofh = (l2-1)*nhs + iofh(j)
Ch        ofh = (l2-1+ofhj)*nhs + ofhi
          do  22  l1 = 1, nlmi
   22     h(l1+ofh) = h(l1+ofh) + h2cj(l2,l1+ofci)
        endif
   20   continue

C --- Accumulate 3C terms to hamiltonian ---
      elseif (mod(mode,10) == 3) then
        do  30  i = 1, nclus
        nlmi = ioffc(i+1,ie)-ioffc(i,ie)
Ch      ib = iax(2,ic-1+iax3(i))
Ch      nlmi = ioffb(ib+1,ie)-ioffb(ib,ie)
        ofci = ioffc(i,ie)-ioffc(1,ie)
        do  30  j = 1, nclus
          nlmj = ioffc(j+1,je)-ioffc(j,je)
Ch        jb = iax(2,ic-1+iax3(j))
Ch        nlmj = ioffb(jb+1,je)-ioffb(jb,je)
          ofcj = ioffc(j,je)-ioffc(1,je)
          is = ipr(i,j)
C     ... skip blocks in lower triangle or if 2C pair is missing
          if (is == 0) goto 30
          lxx = (iofh(is) >= iofh(iax(6,is))) .or. ie /= je
Ch        ofhi = ioffb(ib,ie)
Ch        ofhj = ioffb(jb,je)
Ch        if (lxx .xor. (ofhi <= ofhj)) call rx('problem here')
Ch        lxx = (ofhi <= ofhj)
          if (.not. lxx) then
          elseif (ofci > ofcj .and. ie == je) then
            do  132  l2 = 1, nlmj
            ofh = (l2-1)*nhs + iofh(is)
            do  132  l1 = 1, nlmi
  132       h(l1+ofh) = h(l1+ofh) + h3c(l2+ofcj,l1+ofci)
          else
            do  32  l2 = 1, nlmj
            ofh = (l2-1)*nhs + iofh(is)
            do  32  l1 = 1, nlmi
   32       h(l1+ofh) = h(l1+ofh) + h3c(l1+ofci,l2+ofcj)
          endif
   30   continue
      else
        call rx('hscop: bad mode, man')
      endif
C     call prmx('h in hscop',h,nhs,nhs,nhs)

      end
