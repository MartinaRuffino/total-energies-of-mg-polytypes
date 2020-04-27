      subroutine sudmtu(nbas,nsp,nlibu,lmaxu,s_site,s_spec,
     .  idvsh,lldau,ng,g,istab,dmatu,vorb)
C- Initialize site density matrix and vorb for LDA+U
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  praldm rotycs symdmu
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa idu uh jh name
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  praldm rotycs symdmu
Ci Inputs
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nlibu : nlibu total number of U blocks
Ci   lmaxu :dimensioning parameter for U matrix
Ci   idvsh :0 dmatu and vorb returned in real harmonics
Ci         :1 dmatu and vorb returned in spherical harmonics
Ci   lldau :lldau(ib)=0 => no U on this site otherwise
Ci         :U on site ib with dmat in dmats(*,lldau(ib))
Ci   ng    :number of group operations
Ci   g     :point group operations
Ci   istab :site istab(i,ig) is transformed into site i by grp op ig
Co Outputs
Co   dmatu :density matrix for LDA+U orbitals
Co         :in real spherical harmonics basis
Co   vorb  :orbital dependent potential matrices
Co         :in real spherical harmonics basis
Cs Command-line switches
Cs   --nosymdm : Not documented
Cs   --wdmats  : write dmats to file (useful when starting from occnum)
Cl Local variables
Cl   eorb  : U contribution to LDA+U total energy
Cr Remarks
Cr   Reads in diagonal occupation numbers from file occnum.ext or dmatu
Cr   given in order of m=-l,l, isp=1,2, and constructs initial vorb
Cu Updates
Cu   17 Jun 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
Cu   12 Nov 07 Generate dmatu and vorb in either real or spher. harmonics
Cu   07 May 07 Bug fix MPI mode, when reading occnum instead of dmats
Cu   31 Jan 06 Setup and printouts in spherical harmonics
Cu   09 Nov 05 (wrl) Convert dmat to complex form
Cu   27 Apr 05 Lambrecht first created
C-------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,nsp,nlibu,lmaxu,ng,idvsh
      integer lldau(nbas),istab(nbas,ng)
      double precision g(9,ng)
      double complex dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      double complex Vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
C ... Local parameters
C     logical haveocc
      integer i,isp,ib,l,lmxa,m,m2,foccn,havesh,ivsiz,ipr
C     integer nglob,stdo
      double precision nocc(-3:3,2),iv(7)
      integer idu(4),is,idmat,iblu,nlm,nlmu,a2vec
      double precision uh(4),jh(4),eorb,tmp(2,7,7),xx
      character str*80,spid*8
      complex(8), allocatable :: dmwk(:)
      procedure(logical) rdstrn,parstr,cmdopt
      procedure(integer) :: fopna,fxst,rdm
C ... MPI
      integer procid,master,mpipid
      logical mlog
C ... External calls
      external awrit2,daxpy,dcopy,dpzero,fclose,getpr,info0,info2,info8,
     .         ldau,mpibc1,praldm,rotycs,rxi,rxx,symdmu,zmscop


C     stdo = nglob('stdo')
      call rxx(nsp /= 2,'LDA+U must be spin-polarized!')
      procid = mpipid(1)
      master = 0
      mlog = .false.
      call getpr(ipr)
C     ipr = 55

C --- Read in dmatu if file  dmats.ext  exists ---
      if (procid == master) then
      if (fxst('dmats') == 1) then
        idmat = fopna('dmats',-1,1)
        rewind idmat
  825   continue
        if (rdstrn(idmat,str,len(str),.false.)) then
          if (str(1:1) == '#') goto 825
          i = 0
          if (parstr(str,'sharm ',len(str)-6,5,' ',i,m)) then
            havesh = 1
          else
            havesh = 0
          endif
        endif
        call info2(20,0,0,' sudmtu:  reading density matrix from file'//
     .    ' dmats in %?#n#spherical#real# harmonics',havesh,0)
        backspace idmat
        iblu = 0
        do  ib = 1, nbas
        if (lldau(ib) /= 0) then
          is = s_site(ib)%spec
          lmxa = s_spec(is)%lmxa
          idu = s_spec(is)%idu
          do  l = 0, min(lmxa,3)
          if (idu(l+1) /= 0) then
            iblu = iblu+1
            nlm = 2*l+1
            nlmu = 2*lmaxu+1
            do  isp = 1, 2
              if (rdm(idmat,40,2*nlm**2,' ',tmp,nlm,nlm) /= 2)
     .          call rxi('sudmtu failed to read dmats for site',ib)
C              print *, 'dmat for ib=',ib,' spin',isp
C              call zprm('sudmtu, dmatu',2,tmp,nlm,nlm,nlm)
C              call dmscop(dmatu(-l,-l,isp,iblu),
C     .          nlmu,tmp,nlm,1,nlm,1,nlm,1,1,1d0)
              call zmscop(0,nlm,nlm,nlm,nlmu,0,0,0,0,tmp,
     .          dmatu(-l,-l,isp,iblu))
            enddo
          endif
          enddo
        endif
        enddo

C --- Otherwise read in occnumbers and construct diagonal dmatu ---
      else
        call info0(20,1,0,
     .  ' sudmtu:  no file dmats ... read (diagonal) density-matrix from occnum file')
        if (fxst('occnum') == 0) then
        call info0(20,0,0,
     .    '%10fno file occnum ... initial density matrix will be 1/2 !')
          foccn = 0
          nocc = 0.5d0
        else
          foccn = fopna('occnum',-1,1)
          rewind foccn
        endif
        havesh = 1
   12   continue
        if (foccn /= 0) then
        if (.not. rdstrn(foccn,str,len(str),.false.)) goto 99
        if (str(1:1) == '#') goto 12
        if (str(1:1) == '%') then
          i = 0
          if (parstr(str,'real ',len(str)-5,4,' ',i,m)) havesh = 0
        else
          backspace foccn
        endif
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
            if (foccn /= 0) then
            do  isp = 1, 2
   11         continue
              if (.not. rdstrn(foccn,str,len(str),.false.)) goto 99
C             Skip comment lines
              if (str(1:1) == '#') goto 11
              i = 0
              m = a2vec(str,len(str),i,4,', ',2,3,2*l+1,iv,nocc(-l,isp))
              if (m < 0) goto 99
            enddo
            endif
            call info8(20,0,0,' occ numbers, site %i l=%i:  '//
     .        '%n:1d (spin 1)  %n:1d (spin 2)',
     .        ib,l,2*l+1,nocc(-l,1),2*l+1,nocc(-l,2),0,0)
            do  isp = 1, 2
              do  m = -l, l
                do  m2 = -l, l
                  dmatu(m,m2,isp,iblu) = dcmplx(0d0,0d0)
                enddo
                dmatu(m,m,isp,iblu) = dcmplx(nocc(m,isp),0d0)
C                In order do this, assign new modes 6,7.
C                Such a scheme cannot be used self-consistently
C                if (idu(l+1) == 4) then
C                  dmatu(m,m,1,iblu) = nocc(m,isp)*uh(l+1)
C                  dmatu(m,m,2,iblu) = nocc(m,isp)*uh(l+1)
C                elseif (idu(l+1) == 5) then
C                  dmatu(m,m,1,iblu) = nocc(m,isp)*uh(l+1)
C                  dmatu(m,m,2,iblu) = nocc(m,isp)*jh(l+1)
C                endif
              enddo
            enddo
          endif
          enddo
        endif
      enddo
      if (foccn /= 0) call fclose(foccn)
      endif
C ... Initial printout
      call praldm(0,51,51,havesh,nbas,nsp,lmaxu,lldau,s_spec,s_site,' dmats read from disk',dmatu)
      if (cmdopt('--wdmats',8,0,str)) then
        call info0(30,0,1,' sudmtu:  writing initial dmats to file')
        i = fopna('dmats',-1,0)
        rewind i
        call praldm(i,0,0,0,nbas,nsp,lmaxu,lldau,s_spec,s_site,' initial dmats',dmatu)
        call fclose(i)
      endif
      endif ! procid == master

      ivsiz = nsp*nlibu*(lmaxu*2+1)**2
      call mpibc1(dmatu,2*ivsiz,4,mlog,'sudmtu','dmatu')
      call mpibc1(havesh,1,2,.false.,' ',' ')

C ... Density matrix in real or spherical harmonics (fixed by idvsh)
      if (havesh /= idvsh) then
        call rotycs(2*idvsh-1,dmatu,nbas,nsp,lmaxu,s_spec,s_site,lldau)
        havesh = idvsh
      endif

C ... Symmetrize dmatu (symdmu requires real harmonics)
      allocate(dmwk(ivsiz)); call dpzero(dmwk,2*ivsiz)
      if (havesh == 1) then
        call rotycs(-1,dmatu,nbas,nsp,lmaxu,s_spec,s_site,lldau)
        havesh = 0
      endif
      call symdmu(dmatu,dmwk,nbas,nsp,lmaxu,s_spec,s_site,ng,g,istab,lldau,xx)
      if (havesh /= idvsh) then
        call rotycs(2*idvsh-1,dmatu,nbas,nsp,lmaxu,s_spec,s_site,lldau)
        call rotycs(2*idvsh-1,dmwk,nbas,nsp,lmaxu,s_spec,s_site,lldau)
        havesh = idvsh
      endif
      if (ng /= 0) then
      call info2(30,0,0,' sudmtu:  RMS change in dmats from symmetrization = %,6d',xx,0)
      if (xx > .01d0) call info0(30,0,0,'          (warning) RMS change unexpectely large')
      call daxpy(ivsiz*2,-1d0,dmatu,1,dmwk,1)
      call info0(60,0,0,' change in dmat wrought by symmetrization')
      call praldm(0,60,60,0,nbas,nsp,lmaxu,lldau,s_spec,s_site,' ',dmwk)
      endif

C     Print dmats in specified harmonics
      call dpzero(dmwk,2*ivsiz)
      call dcopy(ivsiz*2,dmatu,1,dmwk,1)
      if (havesh /= idvsh) then
        call rotycs(2*idvsh-1,dmwk,nbas,nsp,lmaxu,s_spec,s_site,lldau)
      endif
      call info0(30,0,0,' ')
      call praldm(0,30,30,idvsh,nbas,nsp,lmaxu,lldau,
     .  s_spec,s_site,' Symmetrized dmats',dmwk)

C     Print dmats in complementary harmonics
      if (ipr >= 45 .or. ipr >= 40 .and. idvsh == 0) then
        i = 1-idvsh
        call rotycs(2*i-1,dmwk,nbas,nsp,lmaxu,s_spec,s_site,lldau)
        call info0(30,0,0,' ')
        call praldm(0,30,30,i,nbas,nsp,lmaxu,lldau,s_spec,s_site,
     .    ' Symmetrized dmats',dmwk)
      endif
      deallocate(dmwk)

C ... Make Vorb (ldau requires spherical harmonics)
      if (havesh /= 1) then
        call rotycs(1,dmatu,nbas,nsp,lmaxu,s_spec,s_site,lldau)
        havesh = 1
      endif
      iblu = 0
      do  20  ib = 1, nbas
        if (lldau(ib) == 0) goto 20
        is = s_site(ib)%spec
        lmxa = s_spec(is)%lmxa
        idu = s_spec(is)%idu
        uh = s_spec(is)%uh
        jh = s_spec(is)%jh
        spid = s_spec(is)%name
        i = min(lmxa,3)
        call info8(30,1,0,' Species '//spid//
     .    '%a: mode=%n:1i    U=%n:1d    J=%n:1d',
     .    i+1,idu,i+1,uh,i+1,jh,0,0)
        do  22  l = 0, i
          if (idu(l+1) /= 0) then
            iblu = iblu+1
            call ldau(mod(idu(l+1),10),l,iblu,uh(l+1),jh(l+1),dmatu,nsp,
     .                lmaxu,vorb,eorb)
C            call prdmts(0,30,30,0,' sudmtu:  Vorb in spherical harm',ib,
C     .        l,lmaxu,iblu,vorb,nsp,1)

          endif
   22   continue
   20 continue
      call praldm(0,60,60,havesh,nbas,nsp,lmaxu,lldau,
     .  s_spec,s_site,' Unsymmetrized vorb',vorb)
C     At this point, dmatu and vorb are in spherical harmonics

C ... Symmetrize vorb to check (symdmu requires real harmonics)
      allocate(dmwk(ivsiz)); call dpzero(dmwk,2*ivsiz)
      call rotycs(-1,vorb,nbas,nsp,lmaxu,s_spec,s_site,lldau)
      call symdmu(vorb,dmwk,nbas,nsp,lmaxu,s_spec,s_site,ng,g,istab,
     .  lldau,xx)
      deallocate(dmwk)

C ... Exit with vorb,dmatu in spherical or real harmonics
C     EITHER: vorb =>  spherical harmonics OR dmatu => real harmonics
      if (idvsh == 1) then
        call rotycs(1,vorb,nbas,nsp,lmaxu,s_spec,s_site,lldau)
        havesh = 1
      else
        call rotycs(-1,dmatu,nbas,nsp,lmaxu,s_spec,s_site,lldau)
        havesh = 0
      endif
      if (ng /= 0) then
      call info2(30,1,0,' sudmtu:  RMS change in vorb '//
     .  'from symmetrization = %,6d',xx,0)
      if (xx > .01d0) call info0(30,0,0,
     .  '          (warning) RMS change unexpectedly large')
      endif

C     Print vorb in specified harmonics
      call praldm(0,30,30,havesh,nbas,nsp,lmaxu,lldau,
     .  s_spec,s_site,' Symmetrized vorb',vorb)

C     Print vorb in complementary harmonics
      if (ipr >= 45 .or. ipr >= 40 .and. idvsh == 0) then
        i = 1-idvsh
        allocate(dmwk(ivsiz)); call dpzero(dmwk,2*ivsiz)
        call dcopy(ivsiz*2,vorb,1,dmwk,1)
        call rotycs(2*i-1,dmwk,nbas,nsp,lmaxu,s_spec,s_site,lldau)
        call info0(30,0,0,' ')
        call praldm(0,30,30,i,nbas,nsp,lmaxu,lldau,
     .    s_spec,s_site,' Vorb',dmwk)
        deallocate(dmwk)
      endif

      eorb = 0

C     call prmx('vorb',vorb(1,1,1,1),2*lmaxu+1,2*lmaxu+1,2*lmaxu+1)

C --- Error exit ---
      return
   99 continue
      call awrit2('bad occnum file, site %i, l=%i',str,len(str),0,ib,l)
      call rx(str)

      end
