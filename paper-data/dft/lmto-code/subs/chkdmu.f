      subroutine chkdmu(nbas,nsp,nlibu,lmaxu,s_site,s_spec,s_ham,idvsh,
     .  dmatu,dmatuo,vorb,tolu,umix,lldau,ng,g,istab)
C- LDA+U total energy and mixing of lda+U density matrix
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read: spec
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: praldm rotycs symdmu
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: lmxa idu uh jh
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: praldm rotycs symdmu
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read: eterms
Co     Stored:    ehk eterms
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: *
Ci Inputs
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nlibu
Ci   lmaxu :dimensioning parameter for U matrix
Ci   idvsh :0 dmatu, dmatuo, vorb input/output in real harmonics
Ci         :1 dmatu, dmatuo, vorb input/output in spherical harmonics
Ci   dmatu : dmatu produced in current iteration
Ci         : dmatu is passed in real harmonics
Ci   dmatuo: dmatu produced in prior iteration
Ci         : dmatuo is passed in real harmonics
Ci   tolu  :convergence tolerance density-matrix
Ci   umix  :linear mixing parameter for density matrix
Ci   lldau :lldau(ib)=0 => no U on this site otherwise
Ci         :U on site ib with dmat in dmats(*,lldau(ib))
Ci   ng    :number of group operations
Ci   g     :point group operations
Ci   istab :site istab(i,ig) is transformed into site i by grp op ig
Cio Inputs/Outputs
Cio  vorb  :orbital dependent potential matrices
Cio        :vorb is updated on output
Cs Command-line switches
Cs   --nosymdm : Not documented
Cl Local variables
Cl   eorb  : U contribution to LDA+U total energy
Cb Bugs
Cb   This routine should not update vorb
Cu Updates
Cu   17 Jun 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
Cu   12 Nov 07 Works with dmatu and vorb in either real or spher. harmonics
Cu   31 Jan 06 Printouts in spherical harmonics
Cu   09 Nov 05 Convert dmat to complex form
Cu   29 Oct 05 doesn't update vorb in tot. E eval; restores dmatu if conv.
Cu    2 Jun 05 Evaluates total energy contribution from output dmatu
Cu   27 Apr 05 Lambrecht first created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,nsp,nlibu,lmaxu,ng,idvsh
      integer lldau(nbas),istab(nbas,ng)
      double precision tolu,umix,g(9,ng)
      double complex dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      double complex dmatuo(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      double precision eorb
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_ham)::   s_ham
C ... Local parameters
      integer l,idu(4),lmxa,ib,is,iblu,idmat,fopna,ivsiz
      integer nglob,stdl,iprint,ipl,havesh
      double precision ddmat,uh(4),jh(4),eorbi,eterms(22),eks,ddot,xx
      equivalence (eterms(2),eks)
C ... Dynamically allocated arrays
      complex(8), pointer :: dmwk(:)
C ... MPI
      integer procid,master,mpipid
C     logical mlog,cmdopt

C     print *, 'dmatu',dmatu(-lmaxu,-lmaxu,1,1)
C     call rx0('done')
C     call prmx('dmatu',dmatu,2*lmaxu+1,2*lmaxu+1,2*lmaxu+1)

      if (nlibu == 0) return
C     stdo = nglob('stdo')
      havesh = idvsh
      stdl = nglob('stdl')
      ipl = 1
      ivsiz = nsp*nlibu*(lmaxu*2+1)**2
      call info0(20,1,0,' chkdmu:  '//
     .  'check LDA+U density-matrix for convergence and update ...')
C ... MPI
      procid = mpipid(1)
      master = 0
C     mlog = cmdopt('--mlog',6,0,strn)

C --- Symmetrize output dmatu (req. real harmonics); compare diff ---
      call praldm(0,60,60,havesh,nbas,nsp,lmaxu,lldau,
     .  s_spec,s_site,' Unsymmetrized output dmats',dmatu)
      if (havesh == 1) then
        call rotycs(-1,dmatu,nbas,nsp,lmaxu,s_spec,s_site,lldau)
        havesh = 0
      endif
      allocate(dmwk(ivsiz)); call dpzero(dmwk,2*ivsiz)
      call symdmu(dmatu,dmwk,nbas,nsp,lmaxu,s_spec,s_site,ng,g,istab,
     .  lldau,xx)
      deallocate(dmwk)
      if (havesh /= idvsh) then
        call rotycs(2*idvsh-1,dmatu,nbas,nsp,lmaxu,
     .    s_spec,s_site,lldau)
        havesh = idvsh
      endif
C     RMS change : dmatu-dmatuo; restore dmatuo
      call daxpy(2*ivsiz,-1d0,dmatu,1,dmatuo,1)
      ddmat = dsqrt(ddot(2*ivsiz,dmatuo,1,dmatuo,1)/(2*ivsiz))
      call daxpy(2*ivsiz,1d0,dmatu,1,dmatuo,1)

C --- Printout dmatu in real or spherical harmonics, fixed by idvsh ---
      call info2(30,0,0,' chkdmu:  RMS change in dmat'//
     .  ' from symmetrization = %,6d',xx,0)

C --- Compute U contribution to total energy; make vorb ---
C     This block requires dmatu to be in spherical harmonics
      if (havesh /= 1) then
        call rotycs(1,dmatu,nbas,nsp,lmaxu,s_spec,s_site,lldau)
        havesh = 1
      endif
      call info0(20,1,0,'%9pLDA+U total energy ...')
      eorb = 0
      iblu = 0
      do  ib = 1, nbas
      if (lldau(ib) /= 0) then
        is = s_site(ib)%spec
        lmxa = s_spec(is)%lmxa
        idu = s_spec(is)%idu
        uh = s_spec(is)%uh
        jh = s_spec(is)%jh
        do  l = 0, min(lmxa,3)
          if (idu(l+1) /= 0) then
            iblu = iblu+1
            eorbi = 999
            call ldau(100+mod(idu(l+1),10),l,iblu,uh(l+1),jh(l+1),dmatu,
     .        nsp,lmaxu,vorb,eorbi)
            eorb = eorb + eorbi
          endif
        enddo
      endif
      enddo

C --- LDA total energy terms ---
      eterms = s_ham%eterms
      eks = eks + eorb
      s_ham%ehk = eks
      s_ham%eterms = eterms
      call info5(20,0,0,'%9peks = %,6;6d  '//
     .  'e[U] = %,6;6d  Etot(LDA+U) = %,6;6d',eks-eorb,eorb,eks,0,0)
      if (mpipid(1) == 0 .and. ipl > 1)
     .  write (stdl,720) eks-eorb,eorb,eks
  720 format('ldau EHK ',f14.6,'  U',f12.6,'  ELDA+U ',f14.6)

C --- Restore dmatu, vorb to harmonics specified by idvsh ---
      if (havesh /= idvsh) then
        call rotycs(2*idvsh-1,dmatu,nbas,nsp,lmaxu,
     .    s_spec,s_site,lldau)
        call rotycs(2*idvsh-1,vorb,nbas,nsp,lmaxu,
     .    s_spec,s_site,lldau)
        havesh = idvsh
      endif

C --- Case self-consistent within tolerance tolu ---
      if (ddmat < tolu) then
        call info5(20,0,0,' LDA+U dmatu converged'//
     .    '  RMS diff (%;3g) < tolu (%;3g)',ddmat,tolu,0,0,0)
C       Restore dmatuo to dmatu
        call dcopy(2*ivsiz,dmatuo,1,dmatu,1)
        return

C --- Case no mixing ---
      elseif (umix == 0) then
        call info2(20,0,0,'%9fRMS diff in dens mat(%;3g) (no mixing)',
     .    ddmat,tolu)
C       Restore dmatuo to dmatu
        call dcopy(2*ivsiz,dmatuo,1,dmatu,1)
        return

C --- Case not self-consistent ---
      else
        call info0(20,0,0,'%9pLDA+U update density matrix ...')
        call info5(20,0,0,'%9fRMS diff in dens mat(%;3g) > tolu (%;3g)'
     .    //' Linear mix with beta=%;3g',ddmat,tolu,umix,0,0)

C   ... Make new dmatu by mixing    new*umix + old*(1-umix)
        call dscal(2*ivsiz,umix,dmatu,1)
        call daxpy(2*ivsiz,1-umix,dmatuo,1,dmatu,1)

C   --- Make Vorb from mixed dmatu ---
C       This block requires dmatu to be in spherical harmonics
        if (havesh /= 1) then
          call rotycs(1,dmatu,nbas,nsp,lmaxu,s_spec,s_site,lldau)
          havesh = 1
        endif
        iblu = 0
        do  ib = 1, nbas
        if (lldau(ib) /= 0) then
          is = s_site(ib)%spec
          lmxa = s_spec(is)%lmxa
          idu = s_spec(is)%idu
          uh = s_spec(is)%uh
          jh = s_spec(is)%jh
          do  l = 0, min(lmxa,3)
            if (idu(l+1) /= 0) then
              iblu = iblu+1
              call pshpr(iprint()-20)
              call ldau(mod(idu(l+1),10),l,iblu,uh(l+1),jh(l+1),dmatu,
     .          nsp,lmaxu,vorb,eorb)
              call poppr
            endif
          enddo
        endif
        enddo
C       At this point, dmatu and vorb are in spherical harmonics

C   ... Symmetrize vorb to check (symdmu requires real harmonics)
        allocate(dmwk(ivsiz)); call dpzero(dmwk,2*ivsiz)
        call rotycs(-1,vorb,nbas,nsp,lmaxu,s_spec,s_site,lldau)
        call symdmu(vorb,dmwk,nbas,nsp,lmaxu,s_spec,s_site,ng,g,istab,
     .    lldau,xx)
        call rotycs(1,vorb,nbas,nsp,lmaxu,s_spec,s_site,lldau)
        deallocate(dmwk)
C       At this point, dmatu and vorb are in spherical harmonics

C   ... Printout
        call info2(20,0,0,'         RMS change in vorb '//
     .    'from symmetrization = %,6d',xx,0)
        if (xx > .0001d0) call info0(30,0,0,
     .    '         (warning) RMS change unexpectely large')
        call praldm(0,30,30,havesh,nbas,nsp,lmaxu,lldau,
     .    s_spec,s_site,
     .    ' Mixed dmats',dmatu)
        call praldm(0,30,30,havesh,nbas,nsp,lmaxu,lldau,
     .    s_spec,s_site,
     .    ' New vorb',vorb)

C   ... Write dmatu to file
        if (procid == master) then
          idmat = fopna('dmats',-1,0)
          rewind idmat
          call praldm(idmat,0,0,havesh,nbas,nsp,lmaxu,lldau,
     .      s_spec,s_site,' mixed dmats',dmatu)
          call fclose(idmat)
        endif

C   ... Exit with dmatu, vorb in real harmonics, depending on idvsh
        if (idvsh == 0) then
          call rotycs(-1,dmatu,nbas,nsp,lmaxu,s_spec,s_site,lldau)
          call rotycs(-1,vorb,nbas,nsp,lmaxu,s_spec,s_site,lldau)
          call praldm(0,60,60,havesh,nbas,nsp,lmaxu,lldau,
     .      s_spec,s_site,' Mixed dmats',dmatu)
          call praldm(0,60,60,havesh,nbas,nsp,lmaxu,lldau,
     .      s_spec,s_site,' New vorb',vorb)
        endif

      endif

      end
