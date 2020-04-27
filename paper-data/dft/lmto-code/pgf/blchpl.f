      subroutine blchpl(lbloch,qp,nl,plat,iprm,ldham,ipl,npl,pgplp,
     .  iax,ntab,s,nds,ldl,ldi,ldl2,klu,sll,sil,sii)
C- Bloch transform strux connecting to a single principal layer
C ----------------------------------------------------------------------
Ci Inputs
Ci   lbloch:1s digit pertains to storage of Bloch summed hamiltonian
Ci           0: s is stored in unpacked form
Ci           1: s is stored in banded form (not implemented;see Remarks)
Ci          10s digit distinguishes how complex arithmetic is handled
Ci           0: sk has real, imaginary separated
Ci              sk = sk(ldl,ldl2,2), with sk(*,*,1..2) = real..imag
Ci           1: sk is returned complex*16 format:
Ci              sk = sk(2,ldl,ldl2), with sk(1..2,*,*) = real..imag
Ci           2: sk has real, imaginary separated by columns
Ci              sk = sk(ldl,2,ldl2), with sk(*,1..2,*) = real..imag
Ci
Ci        100s digit
Ci           1 if to add to s (s not initialized to zero)
Ci           2 subtract from s
Ci           3 combination of 1+2
Ci
Ci       1000s digit 1 if to convert s to spherical harmonics
Ci
Ci   qp    :k-point (two-dimensional)
Ci
Ci   nl    :(global maximum l) + 1
Ci
Ci   plat  :primitive lattice vectors, in units of alat
Ci
Ci   iprm  :permutation indices ordering orbitals in downfolding order.
Ci          Unpermuted orbitals run from 1 ... n, with indices to
Ci          orbitals corresponding to site ib starting at 1+nl*nl*(ib-1).
Ci          Orbital indices i for which 0<iprm(i)<=ldim are accumulated
Ci          into the lower set; those for which ldim<iprm(i)<=ldim+idim
Ci          are accumulated the intermediate set.  Indices for which
Ci          iprm(i) lie outside this range are not accumulated.
Ci          ldim is specified by pgplp.
Ci
Ci   ldham :vector describing hamiltonian dimensions:
Ci         :used for shifting subblocks in a PL
Ci         :1: ldim   = dimension of lmto basis
Ci         :2: lidim  = ldim + size of downfolding block
Ci         :3: lidhim = lidim + size of higher block
Ci
Ci   ipl   :principal layer.
Ci
Ci   pgplp :index and dimensioning information for each PL (pgfset.f)
Ci
Ci   iax   :neighbor table containing pair information for a padded basis.
Ci          Pairs in the neighbor table coupling layer ipl to itself
Ci          and to ipl-1 and ipl+1 are used.
Ci          See pgfset.f for description of a padded basis;
Ci          see pairc.f for generation of iax table.
Ci
Ci   ntab  :ntab(ib)=offset to neighbor table for site ib (pairc.f)
Ci
Ci   s     :real-space matrix to be Bloch summed
Ci
Ci   nds   :leading dimension of s
Ci
Ci   ldl   :leading dimension of sll
Ci
Ci   ldi   :leading and second dimension of sii
Ci
Ci   ldl2  :second dimension of sll and sil
Ci
Ci   klu   :size of sub- and super-diagonal, if s stored banded form
Ci
Cio Outputs: :lower-lower block of Bloch summed matrix
Co
Co   sil   :lower-intermediate block of Bloch summed matrix
Co
Co   sii   :intermediate-intermediate block of Bloch summed matrix
Co
Cr Remarks
Cr  *This routine assembles a bloch sum of a real-space matrix coupling
Cr   a principal layer to itself and its neighbors, viz
Cr      s(k;r1,l1,r2,l2) = sum_T s(r1,l1,T+r2,l2) * exp(i k . T)
Cr   where r1 and r2 are basis vectors and T = t2-t1 is the difference
Cr   in the (two-dimensional) primitive lattice translation vectors.
Cr
Cl Local variables
Cu Updates
Cu   13 Dec 01 no longer relies on orbitals being ordered by site
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
C     integer iipl
      integer lbloch,nds,nl,ipl,npl,ldl,ldi,ldl2,niax,klu,ldham(3),
     .  ntab(*),iprm(*),pgplp(6,-1:*)
      parameter (niax=10)
      integer iax(niax,1)
      double precision qp(3),plat(3,3)
      double precision s(nds,nds,*)
C     real + imaginary storage mode
      double precision sll(ldl,ldl2,2),sil(ldi,ldl2,2),sii(ldi,ldi,2)

C Local parameters
      integer lblch0,ib1,ib2,is1,is2,nl2,kcplx,ldrs,nglob,
     .  ibL1,ibL2,ibR1,ibR2,offi,offL,offR,idim,ldima,ldiml,ldimr,mxorb
C     integer i,fopna

C --- Setup ---
      nl2 = nl*nl
      mxorb = nglob('mxorb')
      lblch0 = mod(lbloch,10)
      call sanrg(.true.,lblch0,0,0,'blchpl:','1s digit lbloch')
C     Parameters describing treatment of complex arithmetic
      kcplx = mod(lbloch/10,10)
      ldrs = ldl
      if (kcplx > 0) ldrs = ldl*2

      ldima = pgplp(4,ipl)
      ldiml = pgplp(4,max(ipl-1,-1))
      ldimr = pgplp(4,min(ipl+1,npl))
C     No downfolding for now
      if (ldima /= pgplp(5,ipl))
     .  call rx('blchpl not set up for downfolding')
      idim = 0
C ... Range of sites within PL, to to right- and left- connecting PL
      call gtibpl(ipl,npl,pgplp,ib1,ib2)
      call gtibpl(max(ipl-1,-1),npl,pgplp,ibL1,ibL2)
      call gtibpl(min(ipl+1,npl),npl,pgplp,ibR1,ibR2)
C ... Subtract hamiltonian offsets in iprm to ipl, ipl-1, ipl+1
C     because we want to bloch to copy array into strx(1)
      call pghoff(ipl,npl,pgplp,iprm,mxorb,offL,offi,offR)

      call pblch3(0,1,offi,ib1,ib2,mxorb,ldham,iprm)
C     call pblch2(-offi,ib1,ib2,nl2,iprm)
      if (offL >= 0 .and. offL /= offi)
     .  call pblch3(0,1,offL,ibL1,ibL2,mxorb,ldham,iprm)
C    .  call pblch2(-offL,ibL1,ibL2,nl2,iprm)
      if (offR >= 0 .and. offR /= offi)
     .  call pblch3(0,1,offR,ibR1,ibR2,mxorb,ldham,iprm)
C    .  call pblch2(-offR,ibR1,ibR2,nl2,iprm)

C --- Strux within PL ---
      is1 = 1+ntab(ib1)
      is2 = ntab(ib2+1)
C ... Mask pairs outside PL by changing sign of iax(1),iax(2)
      call pblchp(iax,is1,is2,ib1,ib2,ib1,ib2)
      call bloch(lbloch,qp,nl,plat,nl**2,iprm,is1,is2,iax,s,nds,1,1,
     .  ldima,ldima,idim,ldl,ldi,ldl2,klu,sll,sil,sii)
      call pblchp(iax,is1,is2,ib1,ib2,ib1,ib2)

C --- Strux connecting PL to PL-1 ---
      if (ipl >= 0) then
        is1 = 1+ntab(ibL1)
        is2 = ntab(ibL2+1)
        call pblchp(iax,is1,is2,ib1,ib2,ibL1,ibL2)
      else
        is1 = 1+ntab(ib1)
        is2 = ntab(ib2+1)
        call pblchp(iax,is1,is2,ibR1,ibR2,ib1,ib2)
      endif
      call bloch(lbloch+100,qp,nl,plat,nl**2,iprm,is1,is2,iax,s,nds,1,1,
     .  ldima,ldiml,idim,ldl,ldi,ldl2,klu,sll(1+ldrs*ldima,1,1),sil,sii)
      if (ipl >= 0) then
        is1 = 1+ntab(ibL1)
        is2 = ntab(ibL2+1)
        call pblchp(iax,is1,is2,ib1,ib2,ibL1,ibL2)
      else
        is1 = 1+ntab(ib1)
        is2 = ntab(ib2+1)
        call pblchp(iax,is1,is2,ibR1,ibR2,ib1,ib2)
      endif

C --- Strux connecting PL to PL+1 ---
      if (ipl <= npl-1) then
        is1 = 1+ntab(ibR1)
        is2 = ntab(ibR2+1)
        call pblchp(iax,is1,is2,ib1,ib2,ibR1,ibR2)
      else
        is1 = 1+ntab(ib1)
        is2 = ntab(ib2+1)
        call pblchp(iax,is1,is2,ibL1,ibL2,ib1,ib2)
      endif
      call bloch(lbloch+100,qp,nl,plat,nl**2,iprm,is1,is2,iax,s,nds,1,1,
     .  ldima,ldimr,idim,ldl,ldi,ldl2,klu,sll(1+ldrs*(ldima+ldiml),1,1),
     .  sil,sii)
      if (ipl <= npl-1) then
        is1 = 1+ntab(ibR1)
        is2 = ntab(ibR2+1)
        call pblchp(iax,is1,is2,ib1,ib2,ibR1,ibR2)
      else
        is1 = 1+ntab(ib1)
        is2 = ntab(ib2+1)
        call pblchp(iax,is1,is2,ibL1,ibL2,ib1,ib2)
      endif

C      i = fopna('out',-1,0)
C      call ywrm(0,' ',2,i,'(9f15.6)',sll,ldl*ldl2,ldl,ldl,ldl2)
C      call fclose(i)
C      stop

C ... Undo changes to iprm table
      call pblch3(11,1,offi,ib1,ib2,mxorb,ldham,iprm)
C     call pblch2(-offi,ib1,ib2,nl2,iprm)
      if (offL >= 0 .and. offL /= offi)
     .  call pblch3(11,1,offL,ibL1,ibL2,mxorb,ldham,iprm)
C    .  call pblch2(-offL,ibL1,ibL2,nl2,iprm)
      if (offR >= 0 .and. offR /= offi)
     .  call pblch3(11,1,offR,ibR1,ibR2,mxorb,ldham,iprm)
C    .  call pblch2(-offR,ibR1,ibR2,nl2,iprm)

C      call priprm(' ham offsets for current PL',
C     .  1,24,mxorb,ldham,ldham,ldham,iprm)
C      stop

      end
      subroutine pblchp(iax,is1,is2,ia1,ia2,ib1,ib2)
C- Flag pairs in iax table outside specified range
C ----------------------------------------------------------------------
Ci Inputs
Ci  is1,is2:range of pairs to loop over
Ci  ia1,ia2:range of site index for augmentation channel; see Remarks
Ci  ib1,ib2:range of site index for basis channel; see Remarks
Cio Inputs/Outputs
Cio  iax   :neighbor table containing pair information (pairc.f)
Cio        :basis site index iax(1) for pairs is1..is2 is masked
Cio        :for pairs outsite range ib1..ib2 by changing sign of iax(1)
Cio        :augmentation site index iax(1) for pairs is1..is2 is masked
Cio        :for pairs outsite range ia1..ia2 by changing sign of iax(1)
Cr Remarks
Cr   A pair i in the iax table is flagged by changing the sign of
Cr   either iax(1,i) or iax(2,i) or both.  For example,
Cr     iax(1,i) <= 0, bloch excludes this pair from the sum
Cr     iax(2,i) <= 0, bloch excludes this pair from the sum
C ----------------------------------------------------------------------
      implicit none
      integer niax
      parameter (niax=10)
      integer iax(niax,1),is1,is2,ia1,ia2,ib1,ib2
      integer is

      do  is = is1, is2
        if (iax(1,is) < ib1 .or. iax(1,is) > ib2) iax(1,is)=-iax(1,is)
        if (iax(2,is) < ia1 .or. iax(2,is) > ia2) iax(2,is)=-iax(2,is)
      enddo

      end
      subroutine pblch2(offi,ib1,ib2,mxorb,iprm)
C- Add a constant offset to iprm for sites ib1..ib2
C ----------------------------------------------------------------------
Ci Inputs
Ci   offi  :offset to add
Ci   ib1   :add offsets for orbitals in (ib1..ib2)
Ci   ib2   :add offsets for orbitals in (ib1..ib2)
Ci   mxorb :iprm is effectively dimensioned iprm(mxorb,ib2)
Co Outputs
Cio  iprm  :hamiltonian offsets.  On exit iprm for sites ib1:ib2
Cio        :is shifted by offi
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
C ----------------------------------------------------------------------
      implicit none
      integer offi,ib1,ib2,mxorb,iprm(486)
      integer ib,iorb,lmi

      lmi = (ib1-1)*mxorb
      do  ib = ib1, ib2
      do  iorb = 1, mxorb
        lmi = lmi+1
        iprm(lmi) = iprm(lmi) + offi
      enddo
      enddo
      end

      subroutine pblch3(opt,nblk,offi,ib1,ib2,mxorb,ldham,iprmb)
C- Add a block-dependent offset to iprmb for sites ib1..ib2
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :1s digit
Ci         :0 subtract offi from iprmb
Ci         :1 add offi to iprmb
Ci         :10s digit
Ci         :0 initial shift; see Remarks
Ci         :1 Restore initial shift; see Remarks.
Ci   nblk  :number of blocks hamiltonian is partitioned into
Ci   offi  :block-dependent offset to subtract (see Remarks)
Ci         :offi(i) = block i, i=1..nblk
Ci         :Note: it is assumed that offi(i+1) >= offi(i)
Ci   ib1   :add offsets for orbitals in (ib1..ib2)
Ci   ib2   :add offsets for orbitals in (ib1..ib2)
Ci   ldham :ldham(i) = upper limit to offset for block i
Ci         :An orbital is in the ith block if it is not in a lower block,
Ci         :and it falls below the upper limit; but see Remarks
Ci   mxorb :iprmb is effectively dimensioned iprmb(mxorb,ib2)
Co Outputs
Cio  iprmb :hamiltonian offsets.  On exit iprmb for sites ib1:ib2
Cio        :is shifted by offi
Cl Local variables
Cr Remarks
Cr   Orbital iorb belongs to block i if it does not belong to a lower
Cr   block, and:
Cr     if (iprmb(iorb) < ldham(i))         10s digit opt=0 (initial shift)
Cr     if (iprmb(iorb)+offi(i) < ldham(i)) 10s digit opt=1 (restore shift)
Cr   If the orbital does not belong to any block, it is not shifted
Cu Updates
Cu   21 Jan 10 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer opt,nblk,ib1,ib2,mxorb
      integer offi(nblk),ldham(nblk),iprmb(mxorb,*)
C ... Local parameters
      integer ib,iorb,offk,opt0,opt1,i

      opt0 = mod(opt,10)
      opt1 = mod(opt/10,10)
      do  ib = ib1, ib2
      do  iorb = 1, mxorb
        offk = 0
        if (opt1 == 0) then
          do  i = 1, nblk
            if (iprmb(iorb,ib) <= ldham(i)) then
              offk = offi(i)
              exit
            endif
          enddo
        else
          do  i = 1, nblk
            if (iprmb(iorb,ib)+offi(i) <= ldham(i)) then
              offk = offi(i)
              exit
            endif
          enddo
        endif
        if (opt0 == 0) offk = -offk
        iprmb(iorb,ib) = iprmb(iorb,ib) + offk
      enddo
      enddo
      end

      subroutine pghoff(ipl,npl,pgplp,iprmb,mxorb,offL,offi,offR)
C- Generate offsets to full hamiltonian subblock for PL
C ----------------------------------------------------------------------
Ci Inputs
Ci   ipl   : PL for which to calculated offsets
Ci   npl   : total number of PL, needed for pgplp
Ci   pgplp : index and dimensioning information for each PL (pgfset.f)
Ci   iprmb : orbital permutation table
Co Outputs
Co   offL  : offset to ipl-1 (not defined if ipl=-1); returns -1
Co   offi  : offset to ipl
Co   offR  : offset to ipl+1 (not defined if ipl=npl); returns -1
Cu Updates
Cu   13 Dec 01 pghoff no longer relies on orbitals being ordered by site
Cu             New argument list
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ipl,npl,mxorb,pgplp(6,-1:npl),offL,offi,offR,iprmb(mxorb,*)
C ... Local parameters
C     integer isum,ibL1,ibL2,ibR1,ibR2,ldim
      integer k,ib1,ib2

      call gtibpl(ipl,npl,pgplp,ib1,ib2)
      call imxmn((ib2-ib1+1)*mxorb,iprmb(1,ib1),1,offi,k)
      offi = offi-1
      offL = -1
      offR = -1

      if (ipl < npl) then
        call gtibpl(min(ipl+1,npl),npl,pgplp,ib1,ib2)
        call imxmn((ib2-ib1+1)*mxorb,iprmb(1,ib1),1,offR,k)
        offR = offR-1
      endif

      if (ipl > -1) then
        call gtibpl(max(ipl-1,-1),npl,pgplp,ib1,ib2)
        call imxmn((ib2-ib1+1)*mxorb,iprmb(1,ib1),1,offL,k)
        offL = offL-1
      endif

c     debugging ... old way (assumes site ordering of orbitals)
C      if (ipl == -1) then
C        k = offi
C        offi = isum(npl,pgplp(4,0),6)
C        if (offi /= k) call rxi('oops',1)
C        if (offR /= 0) call rxi('oops',11)
C      elseif (ipl == npl) then
C        k = isum(npl,pgplp(4,0),6)
C        if (offi /= k + pgplp(4,0)) call rxi('oops',2)
C        if (offl /= k - pgplp(4,npl)) call rxi('oops',2)
C        if (offR /= -1) call rxi('oops',2)
C      else
C        k = offi
C        offi = isum(ipl,pgplp(4,0),6)
C        if (offi /= k) call rxi('oops',1)
C        k = offL
C        offL = offi - pgplp(4,ipl-1)
C        if (ipl == 0) offL = isum(npl,pgplp(4,0),6)
C        if (offL /= k) call rxi('oops',3)
C        k = offR
C        offR = offi + pgplp(4,ipl)
C        if (ipl == npl-1) offR = offR + pgplp(4,0)
C        if (offR /= k) call rxi('oops',3)
C      endif

      end
