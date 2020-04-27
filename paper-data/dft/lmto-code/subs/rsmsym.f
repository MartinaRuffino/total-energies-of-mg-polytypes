      subroutine rsmsym(lbloch,plat,mxorb,iprmb,ldim,nbas,pos,nl,nsp,
     .  is1,is2,ntab,iax,g,istab,ng,nsafm,nds,s,syms)
C- Symmetrize a real-space matrix according to given group operations
C ----------------------------------------------------------------------
Ci Inputs:
Ci   lbloch    :pertains to how s and syms are represented
Ci             :1s digit:
Ci               0 if s, syms are real
Ci               1 if s, syms are complex
Ci          10s digit:
Ci               0 no additional operations
Ci               1 replace S_RL,R'L' and S_R'L',RL with their average
Ci                 (symmetric complex matrices)
Ci               2 replace S_RL,R'L' and S*_R'L',RL with their average
Ci                 (hermitian matrices)
Ci         100s digit distinguishes how complex arithmetic is handled
Ci              (parameter kcplx below)
Ci           0: sll or s has real, imaginary separated
Ci              sll = sll(ldl,ldl2,2), with sll(*,*,1..2) = real..imag
Ci           1: sll or s are in complex*16 format:
Ci              sll = sll(2,ldl,ldl2), with sll(1..2,*,*) = real..imag
Ci           2: sll or s have real, imaginary separated by columns
Ci              sll = sll(ldl,2,ldl2), with sll(*,1..2,*) = real..imag
Ci       1000s digit not used
Ci           1: permits missing symmetry-equivalent sites
Ci      10000s digit not used
Ci
Ci     100000s digit set when s is not structured in a
Ci             simple s,p,d,.. form (no permutations, 1 l channel/site)
Ci           0 if s has normal order
Ci           1 if s has a permuted orbital order
Ci
Ci   plat,nbas,pos : site information
Ci   mxorb     :leading dimension of iprmb, and the maximum number
Ci             :of orbitals connected with a single site.
Ci             :Special case: mxorb=1
Ci             :The "matrix" s is just a scalar.  In that case:
Ci             :iprmb,ldim are not used
Ci             :nl,nsp  must be 1  (In future, may allow nsp=2)
Ci             :nds must be 1
Ci   iprmb     :permutation indices ordering orbitals in downfolding
Ci             :order. iprmb(1:mxorb,ib) are the permutation indices
Ci             :for site ib.
Ci   ldim      :dimension of hamiltonian. Orbitals for which iprmb>ldim
Ci             :are not included in the basis.
Ci             :Used only when s has a permuted orbital order
Ci   nl        :at least 1+ maximum l for any orbital in the hamiltonian
Ci             :nl does not need to be synchronized with mxorb
Ci   nsp       :2 for spin-polarized case, otherwise 1
Ci   is1,is2   :symmetrize only pairs between is1,is2.  NB: g MUST NOT
Ci              map any pair in (is1,is2) out of that subset!
Ci   ntab,iax  :neighbor table information; see pairc
Ci   iax       :neighbor table containing pair information (pairc.f)
Ci              For each pair i, the following portion is used here:
Ci              iax(1,i): basis atom for source (column) index
Ci                        If <= 0, pair excluded from symmetrization
Ci              iax(2,i): basis atom for augmentation (row) index
Ci                        If <= 0, pair excluded from symmetrization
Ci              iax(3..5,i): lattice vectors separating the two sites
Ci                            as multiples of plat
Ci   g,ng      :point group operations (see spcgrp) and number
Ci             :NOTE: this must NOT include special operations which
Ci             :are present only because of time-reversal symmetry
Ci   istab     :site i is transformed into istab(i,ig) by grp op ig
Ci   nds       :leading dimensions of s and syms
Ci   s         :un-symmetrized matrix
Co Outputs:
Co   syms      :symmetrized matrix
Cl Local variables
Cl   lmiss     :if true, permits missing symmetry-equivalent sites
Cl   lsclr     :if true, special case mxorb=1
Cu Updates
Cu   07 Aug 17 Revise to include AFM symmetrization
Cu   28 Feb 07 Simplified treatment when s is a scalar
Cu   06 Aug 06 change call to symstr
Cu   25 Sep 04 permits missing symmetry-equivalent sites
Cu   26 Sep 03 Bug fix: use grpfnd instead of gpfndx
Cu   09 May 03 Make spin-polarized
Cu   08 May 03 Add 20's digit lbloch (average S and S+)
Cu   10 Jan 03 Extensively revised
C ---------------------------------------------------------------
      use mpi
      implicit none
C ... Passed parameters
      integer, parameter :: niax=10, n0=10, nkap0=4
      integer nbas,nl,ntab(*),nsp,ng,nsafm,is1,is2,lbloch,nds,mxorb,ldim
      integer istab(nbas,*)
      integer iax(niax,*),iprmb(mxorb,*)
      double precision pos(3,*),g(9,*),plat(3,3),s(nds,nds,nsp,2*is2),syms(nds,nds,nsp,2*is2)
C     double precision ag(3,*)
C     double complex symz(nds,nds,nsp,*)
C ... Dynamically allocated local arrays
      integer, dimension(:),allocatable :: igproc
C ... Local parameters
      character strn*120,dc*1
      logical lmiss,lmissi,lsclr
C#ifdefC DEBUG
C      logical :: ldebug = .false.
C#endif
      integer ib,jb,isite,ig,kg,nlm,ibp,jbp,ksite,is,ks,ipr,isp,
     .  ksp,ksps,offz,ix,i,kcplx,scplx,lsympr,lprmr,iprmbl(mxorb),nmiss,sflg(is2)
      double precision qlat(3,3),v(3),rv(3),vk(3),dx,fuzz
      double precision rmat(nl**2,nl**2),wk(mxorb,mxorb),dum
      double complex wkz1(mxorb,mxorb),wkz2(mxorb,mxorb)
C     For orbital permutations
      integer norbi,ltabi(n0*nkap0),ktabi(n0*nkap0),offli(n0*nkap0),ndimi,offri
      integer norbj,ltabj(n0*nkap0),ktabj(n0*nkap0),offlj(n0*nkap0),ndimj,offrj
      integer nlmi,nlmj,iorb,jorb,li,lj,offr,offi,offj,k,j
C     Convert pair index to effective offset
C     offz(i) = i if real; offz(i) = 2*i-1 if complex
      offz(scplx,i) = (scplx+1)*(i-1) + 1
C     ix-th component of connecting vector of pair i
      dx(ix,i) = pos(ix,iax(2,i)) - pos(ix,iax(1,i)) +
     .  plat(ix,1)*iax(3,i) + plat(ix,2)*iax(4,i) + plat(ix,3)*iax(5,i)
C     For MPI
      integer procid, master, numprocs, ierr
      logical mlog
      procedure(logical) :: isanrg,cmdopt
      procedure(integer) :: wordsw,a2vec

C      is = 1155
C      print 333, s(1,1,:,is),s(1,1,:,is+1),s(1,1,:,is+2)


C --- Setup ---
      if (ng <= 0) return
      if (nsafm > 0) lsclr = isanrg(ng,2,2,'rsmsym:','ng, afm case',.true.)
      fuzz = 1d-5
      if (cmdopt('--fixpos',8,0,strn)) then
        dc = strn(9:9)
        ib = wordsw(strn,dc,'tol','= ',jb)
        if (ib /= 0) then
          call rxx(a2vec(strn,len_trim(strn),jb,4,', '//dc,3,3,1,is,fuzz) /= 1,
     .      'failed to parse '//trim(strn))
        endif
      endif

      call tcn('rsmsym')
      call mpi_comm_rank( mpi_comm_world, procid, ierr )
      call mpi_comm_size( mpi_comm_world, numprocs, ierr )
      master = 0
      nmiss = 0
      mlog = cmdopt('--mlog',6,0,strn)
      kcplx = mod(lbloch/100,10)
      scplx = mod(lbloch,10)
      lsympr = mod(lbloch/10,10)
      lprmr = mod(lbloch/100000,10)
      lmiss = mod(lbloch/1000,10) == 1
      lsclr = mxorb <= 1
      call getpr(ipr)
      nlm = nl*nl
      call dinv33(plat,1,qlat,v)
      call info5(20,0,0,' rsmsym: symmetrizing %?#n#complex#real#'//
     .  ' s(%i..%i) using %i group operations %?#n#(AFM)##',scplx,is1,is2,ng,nsafm)
      if (scplx == 0 .and. kcplx /= 0)
     .  call rxi('rsmsym: incompatible switches in lbloch:',lbloch)
C ... Setup norbi,ltabi,ktabi,offri,offli,ndimi when no orbital permutation
      if (lprmr == 0 .and. .not. lsclr) then
        forall (i = 1:mxorb) iprmbl(i) = i
        call orbl(1,0,mxorb,iprmbl,norbi,ltabi,ktabi,offri,offli,ndimi)
        call orbl(1,0,mxorb,iprmbl,norbj,ltabj,ktabj,offrj,offlj,ndimj)
      endif

C --- Contribution from identity symmetry operation ---
      if (procid == 0) then
        call dcopy(nds*nds*(is2-is1+1)*(1+scplx)*nsp,
     .    s(1,1,1,offz(scplx,is1)),1,syms(1,1,1,offz(scplx,is1)),1)
      else
        call dpzero(syms(1,1,1,offz(scplx,is1)),nds*nds*(is2-is1+1)*(1+scplx)*nsp)
      endif
C     if (lsympr /= 0) then
C       do  i = 1, is2
C         sflg(i) = 0
C       enddo
C       call symstr(scplx,nds,is2,iax,nsp,sflg,syms,syms,dum)
C       call info(30,0,0,'         max asymmetry=%d from symstr',dum,0)
C     endif
      if (ng == 1) goto 999

C      isite = 2
C      is = offz(scplx,isite)
C      ksp = 1
C      print *, 'start ... isite=',isite
C      call pvtro9(scplx,1,2,nds,syms(1,1,ksp,is),syms(1,1,ksp,is))
C      if (scplx == 0) call prmx('sr',syms(1,1,ksp,is),nds,nds,nds)
C      if (scplx == 1) call zprm('sz',2,syms(1,1,ksp,is),nds,nds,nds)

C --- For each symmetry operation, do ---
      allocate (igproc(0:numprocs))
      if (numprocs > 1) then
C      call pshpr(ipr-10)
      call dstrbp(ng,numprocs,1,igproc)
C     Alternatively, distribute loops over pairs.  Doesn't improve speed.
C     Probably need to reverse order of ng, isite
C      call dstrbp(is2-is1+1,numprocs,1,igproc)
C      do  isite = 0, numprocs
C        igproc(isite) = igproc(isite) + is1-1
C      enddo
C      call poppr
      else
        igproc(0:1) = [2,ng+1]
      endif

      ! serial:  do  ig = 2, ng
      do  ig = max(igproc(procid),2), igproc(procid+1)-1

        kg = ig ; if (nsafm > 0) kg = nsafm ! AFM case: substitute kg for ig=2

C   ... Rotation matrix for this ig
        call ylmrtg(nlm,g(1,kg),rmat)
C#ifdefC DEBUG
C        if (ldebug) then
C          call prmx('rmat',rmat,nlm,nlm,nlm)
C        endif
C#endif

C   --- For each site, do ---
CC#ifdef MPI
C        do isite = igproc(procid), igproc(procid+1)-1
CC#elseC
CC        do  isite = is1, is2
CC#endif
        do  isite = is1, is2

          ib = iax(1,isite)
          jb = iax(2,isite)
C         Any nonpositive site indices are excluded from symmetrization
          if (jb <= 0 .or. ib <= 0) goto 20

C     ... Find ibp and jbp = sites that ib,jb are rotated into
C         call gpfndx(g(1,kg),ag(1,kg),ib,ibp,pos,nbas,plat,qlat)
C         call gpfndx(g(1,kg),ag(1,kg),jb,jbp,pos,nbas,plat,qlat)

C         This is apparently correct
C         call grpfnd(fuzz,g,ag,kg,pos,nbas,qlat,ib,ibp,dlat)
C         call grpfnd(fuzz,g,ag,kg,pos,nbas,qlat,jb,jbp,dlat)
C         This should accomplish the same result provided
C         site ib is transformed into istab(ib,kg) by grp op kg
          ibp = istab(ib,kg)
          jbp = istab(jb,kg)

C     ... Case permute orbital order in s: offr[ab] = offset to orbitals
          if (lprmr /= 0 .and. .not. lsclr) then
C           korb = 0
            call orbl(ib,0,ldim,iprmb,norbi,ltabi,ktabi,offri,offli,ndimi)
            call orbl(jb,0,ldim,iprmb,norbj,ltabj,ktabj,offrj,offlj,ndimj)
            if (ndimi == 0) goto 20
            if (ndimj == 0) goto 20
            if (ndimi > nds .or. ndimj > nds)
     .        call rx('rsmsym: s dimension incompatible with iprmb')
          else
          endif

C     ... Original connecting vector v, and the rotated one rv
          v(1) = dx(1,isite)
          v(2) = dx(2,isite)
          v(3) = dx(3,isite)
          call rotpnt(v,rv,g(1,kg))

C     ... Find ksite = pair corresponding to ibp,jbp and rv
          lmissi = .false.
          do  ks = ntab(ibp)+1, ntab(ibp+1)
            ksite = ks
C       ... This is the corresponding pair only if jb of this pair = jbp
            if (iax(2,ksite) /= jbp) goto 40
C       ... Also the connecting vector of ksite must be rv
            vk(1) = dx(1,ksite)
            vk(2) = dx(2,ksite)
            vk(3) = dx(3,ksite)
            if (abs(rv(1)-vk(1)) + abs(rv(2)-vk(2)) + abs(rv(3)-vk(3)) < fuzz) goto 42
   40     continue
          enddo
C         Case missing symmetry-equivalent sites permitted
          if (lmiss) then
            isp = 35; if (nmiss > 0) isp = 50
            call info5(isp,0,0,' rsmsym (warning) no connecting vector '//
     .        'for is,ib,jb= %i %i %i  v=(%s,%3;5,5d)',isite,ib,jb,rv,0)
            lmissi = .true.
            ksite = isite
            nmiss = nmiss+1
          else
            call info5(10,0,0,' rsmsym (aborting) no connecting vector '//
     .        'for is,ib,jb= %i %i %i',isite,ib,jb,0,0)
            call rx('RSMSYM: bad or incomplete neighbor table.  Consider --fixpos:tol=...')
          endif
   42     continue

C          if (isite == 2)
C     .      call info8(0,0,0,' ig=%i isite=%i, ib,jb=%i %i  maps to'
C     .      //' site=%i, jb,jbp=%i %i',kg,isite,ib,jb,ksite,ibp,jbp,0)

          if (ipr >= 110 .and. .not. lmissi)
     .      call info8(0,0,0,' rsmsym ig=%i isite=%i, ib,jb=%i %i  maps'
     .     //' to site=%i, jb,jbp=%i %i',kg,isite,ib,jb,ksite,ibp,jbp,0)

C     --- Rotate s(isite) into syms(ksite) ---
          is = offz(scplx,isite)
          ks = offz(scplx,ksite)

          do  isp = 1, nsp
          ksp = offz(scplx,isp)
          ksps = ksp
          if (nsafm /= 0) ksps = offz(scplx,nsp+1-isp) ! Opposite spin if nsp=2
          if (lmissi .and. nsafm /= 0) then
C            print 333, s(1,1,:,is),s(1,1,:,is+1)
C            print 333, syms(1,1,:,is),syms(1,1,:,is+1)
C            print 333, syms(1,1,:,offz(scplx,isite):offz(scplx,isite)+1)
C  333       format(6f15.10)
C            print *, isite,is,offz(scplx,isite),(1+scplx)*nsp
            call dcopy(nds*nds*(1+scplx)*nsp,s(1,1,1,offz(scplx,isite)),1,
     .        syms(1,1,1,offz(scplx,isite)),1)
C           print 333, syms(1,1,:,is),syms(1,1,:,is+1),syms(1,1,:,is+2)
            call dscal(nds*nds*(1+scplx)*nsp,dble(ng),syms(1,1,1,offz(scplx,isite)),1)
C           print 333, syms(1,1,:,is),syms(1,1,:,is+1),syms(1,1,:,is+2)
            cycle
          endif

          if (lsclr) then
            syms(1,1,ksps,ks) = syms(1,1,ksps,ks) + s(1,1,ksp,is)
          else

C         Complex case: copy s to work array wkz1
          if (scplx /= 0) then
C            print *, 'ks is', ksite,ks
C            call zprm('s',2,s(1,1,ksp,ks),nds,nds,nds)
            call zmscop(0,nds,nds,nds,mxorb,0,0,0,0,s(1,1,ksp,is),wkz1)
            if (kcplx /= 1) call ztoyy(wkz1,mxorb,mxorb,nds,nds,kcplx,1)
C#ifdefC DEBUG
C            if (ldebug) then
C              call yprmi('wkz1 isite,ksite=%s,%2i ig=%i',[isite,ksite],ig,3,wkz1,0,mxorb,nds,nds)
C            endif
C#endif
          endif

          if (.not. lmissi) then
C     ... For each iorb, wk(iorb) =  rmat * s(is).  Only rows mix
C         Index jb corresponds to row (augmentation) index
C         Index ib corresponds to col (source) index
C         Example: norbj=3 with lj=0,1,0: product has structure:
C               r                  s
C                          <------ ndimi ------->
C           (x|   | )     (______________________)
C           ( |xxx| )     ( |                    )
C           ( |xxx| )     ( |                    )
C           ( |xxx| )     (_|____________________)
C           ( |   |x)     (_|____________________)
          if (scplx /= 0) call dpzero(wkz2,mxorb**2*2)
          do  jorb = 1, norbj
            lj = ltabj(jorb)
            nlmj = 2*lj + 1
            offr = lj**2
            offj = offlj(jorb) - offrj
            if (scplx == 0) then
            do  k = 1, ndimi
            do  i = 1, nlmj
            wk(i+offj,k) = 0
            do  j = 1, nlmj
              wk(i+offj,k) = wk(i+offj,k) +
     .                       rmat(i+offr,j+offr)*s(j+offj,k,ksp,is)
            enddo
            enddo
            enddo
            else
              do  k = 1, ndimi
              do  i = 1, nlmj
              do  j = 1, nlmj
              wkz2(i+offj,k) = wkz2(i+offj,k) +
     .                         rmat(i+offr,j+offr)*wkz1(j+offj,k)
              enddo
              enddo
              enddo
            endif
          enddo

C     ... For each iorb, syms += wk(iorb) * rmat+.  Only columns mix
C         Index jb corresponds to row (augmentation) index
C         Index ib corresponds to col (source) index
C         Example: norbj=3 with lj=0,1,0: product has structure:
C                      wk                      r+
C
C       ^     (______________________)     (x|   | )
C       |     ( |                    )     ( |xxx| )
C     ndimj   ( |                    )     ( |xxx| )
C       |     (_|____________________)     ( |xxx| )
C       |     (_|____________________)     ( |   |x)
          if (scplx /= 0) call dpzero(wkz1,mxorb**2*2)
          do  iorb = 1, norbi
            li = ltabi(iorb)
            nlmi = 2*li + 1
            offr = li**2
            offi = offli(iorb) - offri
            if (scplx == 0) then
              do  i = 1, nlmi
              do  j = 1, nlmi
              do  k = 1, ndimj
                syms(k,i+offi,ksps,ks) = syms(k,i+offi,ksps,ks) +
     .                                   wk(k,j+offi)*rmat(i+offr,j+offr)
              enddo
              enddo
              enddo
            else
              do  i = 1, nlmi
              do  j = 1, nlmi
              do  k = 1, ndimj
                wkz1(k,i+offi) = wkz1(k,i+offi) +
     .                           wkz2(k,j+offi)*rmat(i+offr,j+offr)
              enddo
              enddo
              enddo
            endif
          enddo

C#ifdefC DEBUG
C          if (ldebug) then
C            call yprmi('wkz1 after rotation isite,ksite=%s,%2i ig=%i',[isite,ksite],ig,3,wkz1,0,mxorb,nds,nds)
C          endif
C#endif

C         Case no matching vector: just add s(ksp,is) to syms(ksps,is)
C         so norm correct when scaling by 1/ng
          else if (scplx == 0) then ! lmiss = T
            call dmscop(1,nds,nds,mxorb,nds,0,0,0,0,syms(1,1,ksps,ks),syms(1,1,ksps,ks))

          endif  ! Finished rotation

C         Complex case: add wkz1 to syms(ksps,ks)
          if (scplx /= 0) then
            if (kcplx /= 1) call ztoyy(wkz1,mxorb,mxorb,nds,nds,1,kcplx)

C            if (ksp == 1 .and. ks == 1) then
C              print 338, isite, kg, wkz1(1,1)
C  338         format(' adding from isite',i4,'  for ig',i4,2f12.6)
C            endif

            call zmscop(1,nds,nds,mxorb,nds,0,0,0,0,wkz1,syms(1,1,ksps,ks))
          endif
        endif

C          i = offz(scplx,ksite)
C          if (ksite == 2)
C     .      call pvtro9(scplx,1,2,nds,syms(1,1,ksps,i),syms(1,1,ksps,i))

C          if (ksite == 3) then
C            print *, 'ig,isite,ksite=',kg,isite,ksite
C            call pvtro9(scplx,5,5,nds,syms(1,1,ksps,ks),syms(1,1,ksps,ks))
C            if (scplx == 0) call prmx('syms',syms(1,1,ksps,ks),nds,nds,nds)
C            if (scplx == 1) call zprm('syms',2,syms(1,1,ksps,ks),nds,nds,nds)
C            if (scplx == 0) call prmx('s(is)',s(1,1,ksp,is),nds,nds,nds)
C            if (scplx == 1) call zprm('s(is)',2,s(1,1,ksp,is),nds,nds,nds)
C            if (scplx == 0) call prmx('s(ks)',s(1,1,ksp,ks),nds,nds,nds)
C            if (scplx == 1) call zprm('s(ks)',2,s(1,1,ksp,ks),nds,nds,nds)
C          endif

C     --- s_rot = rmat s rmat^T: add directly into syms ---
C          is = offz(scplx,isite)
C          ks = offz(scplx,ksite)
C          call dgemm('N','T',nlm,nlm,nlm,1d0,s(1,1,ksp,is),nlm,
C     .      rmat,nlm,0d0,wk,nlm)
C          call dgemm('N','N',nlm,nlm,nlm,1d0,rmat,nlm,
C     .      wk,nlm,1d0,syms(1,1,ksps,ks),nlm)
C
C*         if (ksite == 1) print *, kg,isite,is,ks,s(1,1,ksp,is)
C
C          if (scplx == 1) then
C            call dgemm('N','T',nlm,nlm,nlm,1d0,s(1,1,ksp,is+1),nlm,
C     .        rmat,nlm,0d0,wk,nlm)
C            call dgemm('N','N',nlm,nlm,nlm,1d0,rmat,nlm,
C     .        wk,nlm,1d0,syms(1,1,ksps,ks+1),nlm)
C          endif

C ... debugging
C          if (ksite == 1) then
C          print *, 'ig,isite,ksite=',kg,isite,ksite,s(1,1,ksp,is)
C          call prmx('start',s(1,1,ksp,is),nlm,nlm,nlm)
C          call dcopy(81,rmat,1,dwk,1)
C          call dcopy(81,wk,1,wk,1)
C          call prmx('dest',s(1,1,ksp,ks),nlm,nlm,nlm)
C          call dgemm('N','N',nlm,nlm,nlm,1d0,rmat,nlm,
C     .      wk,nlm,0d0,wk,nlm)
C          call prmx('rotated',wk,nlm,nlm,nlm)
C          call prmx('dest now',syms(1,1,ksps,ks),nlm,nlm,nlm)
C          endif
        enddo                     ! isp loop
   20   continue

C      is = 1155
C      print 334, isite,syms(1,1,:,is),syms(1,1,:,is+1),syms(1,1,:,is+2)
C  334 format('isite',i6,6f15.10)


        enddo    ! isite loop
      enddo      ! ig loop

C ... Combine contributions to syms from separate threads
      call mpibc2(syms(1,1,1,offz(scplx,is1)),
     .  nds*nds*(is2-is1+1)*(1+scplx)*nsp,4,3,mlog,'rsmsym','syms')

C --- Scale by 1/ng ---
      call dscal(nds*nds*(is2-is1+1)*(1+scplx)*nsp,1d0/ng,
     .  syms(1,1,1,offz(scplx,is1)),1)

C     Debugging printout
C      if (procid == master) then
C        v(1) = 0
C        v(2) = 0
C        do  i = is1, is2
C          print 543, i, iax(1,i), iax(2,i),
C     .      s(1,1,1,i),syms(1,1,1,i),(syms(1,1,1,i)-s(1,1,1,i))
C  543     format(3i6,3f12.6)
C          v(1) = v(1) + s(1,1,1,i)
C          v(2) = v(2) + syms(1,1,1,i)
C        enddo
C        call info5(10,0,0,'sum s = %d    sum syms = %d  diff = %d',
C     .    v(1),v(2),v(2)-v(1),0,0)
C      endif
C      call rx0('done')

C      isite = 2
C      ksp = 1
C      is = offz(scplx,isite)
C      print *, 'end ... isite=',isite
C      call pvtro9(scplx,1,2,nds,syms(1,1,ksp,is),syms(1,1,ksp,is))
C      if (scplx == 0) call prmx('sr',syms(1,1,ksp,is),nds,nds,nds)
C      if (scplx == 1) call zprm('sz',2,syms(1,1,ksp,is),nds,nds,nds)

C      print *, 'rotated strux'
C  998 print *, 'is (0 to exit)?'
C      read(*,*) is
C      if (is == 0) goto 999
C      do  i = 1, min(nds,9)
C        print 393, (dble(symz(i,j,is)), j=1,min(nds,9))
C      enddo
C      print 393
C      do  i = 1, min(nds,9)
C        print 393, (dimag(symz(i,j,is)), j=1,min(nds,9))
C      enddo
C  393 format(9f12.6)
C      goto 998

  999 continue
      if (lsympr /= 0) then
        call iinit(sflg,is2)
C       i = 0 for real matrix, 1 to avg (sij,sji), 2 to avg (sij,sji+)
        i = scplx*lsympr
        call symstr(i,nds,is2,iax,nsp,1,sflg,syms,syms,dum)
        call info(30,0,0,' symstr: max asymmetry = %;3g',dum,0)
C       if (dum /= 0) goto 998
      endif
      if (nmiss > 0) then
        call info2(35,0,0,' rsmsym (warning) missed %i connecting vectors',nmiss,2)
      endif

C      is = 1155
C      print 333, syms(1,1,:,is),syms(1,1,:,is+1),syms(1,1,:,is+2)
C  333 format(6f15.10)
C      pause
      call tcx('rsmsym')
C     call rx('done')
      end
