C#define v2
      subroutine asasx(s_ctrl,s_spec,s_site,s_lat,s_ham,s_pot,s_bz,
     .  s_str,s_strn,sxopts,aintra,nl,nsp,efermi,eband,nbmx,zdel,vintra)
C- Calculate ASA static screened response
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbas nclass lasa lncol lham nbasp nspec nl lcgf
Ci                 lgen3 nspin lpgf ldlm nclasp ipc loptc
Co     Stored:     lasa
Co     Allocated:  *
Cio    Elts passed:ipc rmax lfp lgen3 ips ncomp idcc lqp lasa
Cio    Passed to:  subasi suham pp2alp secmat secmtn secmt2
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa lmxl p pz idxdn rmt lmxb ncomp name ngcut orbp
Ci                 hcr
Co     Stored:     idxdn ngcut orbp
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  subasi suham atfold makidx nscpa showbs sugcut
Cio                uspecb secmat secmtn
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec clabel
Co     Stored:     pnu pz norb
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  suham setnorb showbs pvioeu
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  avw alat plat nsgrp awald nkd nkq ng tolft qlat kv
Ci                 kv2 igv igv2 nabc gmax vol
Co     Stored:     igv igv2 ng gmax nabc
Co     Allocated:  gv kv igv kv2 igv2
Cio    Elts passed: pos istab symgr ag dlv qlv gv qlat igv igv2 cg jcg
Cio                indxcg
Cio    Passed to:  suham sugvec sugvec0 secmat secmtn sstrxq
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  ldham iprmb pwmode pwemin pwemax npwpad lncol neula
Ci                 qss nbf ndhrs nmto kmto
Co     Stored:     ndham ndofH ldham lmxax npwmin npwpad hord
Co     Allocated:  offH iprmb bdots
Cio    Elts passed:offH iprmb eula magf bdots nprs iaxs hrs
Cio    Passed to:  subasi suham secmat secmtn sstrxq
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  vmtz
Co     Stored:     nlma nlml
Co     Allocated:  pti
Cio    Elts passed:qnu pp pti ppn sop
Cio    Passed to:  suham secmat secmtn secmt2
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  nkp nkabc lshft zval
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:ipq star qp wtkp
Cio    Passed to:  *
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  npr
Co     Stored:     nkaps nitab nttab lmaxw nds
Co     Allocated:  iax npr alp s adot sdot
Cio    Elts passed:kaps alph nkaps iax s sdot nitab adot
Cio    Passed to:  suham rdstrx pp2alp secmat secmtn secmt2
Cio  s_strn :struct for global strings; see structures.h
Ci     Elts read:  strn
Co     Stored:     strn
Co     Allocated:  *
Cio    Elts passed:strn
Cio    Passed to:  suham str_pack
Ci Inputs
Ci   sxopts
Ci     Passed to: len partok
Ci   aintra:T, include onsite W in SX
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   efermi:Fermi energy
Ci   eband :energy bands; alias eb (sec*.f)
Ci   nbmx
Ci   zdel
Ci   vintra
Co Outputs
Co   Files 'sigm' contains Sigma=SX-LDA potential
Cl Local variables
Cl   ldim  :dimension of hamiltonian
Cl   lpdim :number of channels in response matrix
Cr Remarks
Cr
Cu Updates
Cu   08 May 13 Complete migration to f90 structures; eliminate s_array
Cu   12 Nov 12 migrated structures to f90 pointers
Cu   10 Nov 11 Begin migration to f90 structures
Cu   29 Sep 04 v3 format --- revert to old subtr sigma^loc
Cu   23 Sep 03 SX always writes file sigm, in format compatible
Cu             with GW.  Symmetrization now done at read time
Cu   03 Oct 01 Made to work for neglected orbitals (idxdn=4)
Cu   12 Jun 01 New call to wsloca : better mimic DF subtraction, on-site
C ----------------------------------------------------------------------
      use iolib_cb, only : optio
      use structures
      implicit none
C ... Passed parameters
      character sxopts*(*)
      logical aintra
      integer nl,nsp,nbmx
      double precision eband(nbmx,nsp,*),vintra(nl,nl,nsp,*),
     .  efermi,zdel(2)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
      type(str_lat)::   s_lat
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_bz)::    s_bz
      type(str_str)::   s_str
      type(str_strn) :: s_strn(*)
C ... Dynamically allocated local arrays
C#ifndef LINUXF
      integer, allocatable, target :: iax(:)
C#elseC
C      integer, pointer :: iax(:) => null()
C#endif
      type(str_str0):: s_sts
      integer, allocatable :: lmx(:)
      complex(8), allocatable :: v0(:)
      complex(8), allocatable :: v01(:)
      complex(8), allocatable :: v00(:)
      complex(8), allocatable :: wscr(:)
      complex(8), allocatable :: wscrk(:)
      complex(8), allocatable :: p0(:)
      complex(8), allocatable :: mad(:)
      complex(8), allocatable :: wscrl(:)
      complex(8), allocatable :: s(:)
      complex(8), allocatable :: evec(:)
      complex(8), allocatable :: h(:)
      complex(8), allocatable :: h2(:)
      complex(8), allocatable :: h3(:)
      real(8), allocatable :: pa(:)
      complex(8), allocatable :: v02(:)
      real(8), allocatable :: dedv(:)
      real(8), allocatable :: vsxl(:)
      complex(8), allocatable :: h4(:)
C     real(8), allocatable :: rmat(:)
      complex(8), allocatable :: ecorr(:)
C#ifdefC v1
C      complex(8), allocatable :: p1(:)
C      complex(8), allocatable :: p0l(:)
C      real(8), allocatable :: mad0(:)
C      real(8), allocatable :: mad1(:)
C      real(8), allocatable :: p00(:)
C      real(8), allocatable :: p01(:)
C      real(8), allocatable :: p02(:)
C      real(8), allocatable :: w0(:)
C      real(8), allocatable :: w1(:)
C      real(8), allocatable :: w2(:)
C      real(8), allocatable :: w0m(:)
C      real(8), allocatable :: w1m(:)
C      real(8), allocatable :: w2m(:)
C#endif

C ... Local parameters
      type(dp_wk_vec) :: s_wk(3)
      logical ltmp,lp0,bittst
      character*80 outs
C#ifdefC v1
C      double precision v2
C      integer op00,op01,op02,op0l,op1,omad0,omad1,
C     .  ow0,ow0m,ow1,ow1m,ow2,ow2m
C#endif
C#ifdefC v3
C      integer ov02
C#endif
      logical llshft(3)
      double precision qlat(3,3),vol,dum,rsrnge,avw,alat,plat(3,3),
     .  ddot,qp(3),qp2(3),awald,g(3,3),ag(3),ckbas,cksumf,zval,evtop,
     .  ecbot,qpoff(3),xx
      integer i,idim,ifi,ifisq,ifit,ig,ipr,
     .  iprint,iq,isp,isw,iter,ix,j,l,lasa,ld2,ldham(16),lgunit,ldim,
     .  lidim,lihdim,lncol,lham,lp2,lpdim,lrssig,lshft(3),lsharm,nbas,
     .  nclass,nev,nttab,niter,nkd,nkp,nkr,nkxyz(3),nsgrp,izval,
     .  stdo,vsn,xxia(1)
      integer, parameter :: niax = 10
      equivalence (ldim,ldham(1)),(lidim,ldham(2)),(lihdim,ldham(3))
      integer,allocatable:: ntabs(:)
C     dummy entries for LDA+U
      integer nlibu,lmaxu,vorb(1),idu(1),ludiag
      parameter (nlibu=0,lmaxu=0,ludiag=0)
      procedure(integer) :: bitand,fopn,fopna,fopnT,partok

C#ifdefC v1
C     vsn = 1
C#elseif v2
      vsn = 2
C#elseifC v3
C      vsn = 3
C#elseifC v4
C      vsn = 4
C#endif

C --- Setup, unpack SX options ---
      call tcn('asasx')
      nbas = s_ctrl%nbas
      nclass = s_ctrl%nclass
      lasa = s_ctrl%lasa
      lncol = s_ctrl%lncol
      lham = s_ctrl%lham
      avw = s_lat%avw
      alat = s_lat%alat
      plat = s_lat%plat
      nsgrp = s_lat%nsgrp
      nkp = s_bz%nkp
      nkxyz = s_bz%nkabc
      lshft = s_bz%lshft
      zval = s_bz%zval
      awald = s_lat%awald
      nkd = s_lat%nkd
      nkr = s_lat%nkq
      ckbas = cksumf(s_lat%pos,3*nbas)
      lsharm = 0
      if (bittst(lham,256)) lsharm = 20
      call dpzero(qpoff,3)

      call getpr(ipr)
      stdo = lgunit(1)
C ... Options can be one of the following:
C     'pstat' generates the q=0 response function and exits
C     rssig[=range] averages sigma over k
C     nit=#         iterates for sigma # times.
      call partk0(0,len(sxopts),1,-1,0,len(sxopts),-1,31,.false.)
!       call stlibv(12,11,.false.)
      optio = 11
      i = partok(sxopts,'pstat',' ;',lp0,' ',0,0,0,0)
      lrssig = 0
      niter = 1

      if (lp0) then
        call info0(10,1,0,' --- ASASX : make pstat(q=0) ---')
      else
        i = partok(sxopts,'rssig',' ;',ltmp,' ',0,0,0,0)
        if (ltmp) then
          lrssig = 1
          rsrnge = 6
          i = partok(sxopts,'rssig=','=',rsrnge,' ',-1,4,0,0)
        endif
        i = partok(sxopts,'nit=','= ;',niter,' ',-1,2,0,0)
        call info5(10,1,0,' --- ASASX(v%i): make P0,W, Sig (nit=%i)'//
     .  '%?#n# (rs sigma, r=%d)#%j#%?#n# vintra=T## ---',
     .    vsn,niter,lrssig,rsrnge,
     .    isw(aintra.and.ddot(nl*nl*nsp*nclass,vintra,1,vintra,1) > 0))
      endif
      if (lrssig == 1) call rx('no lrssig')

      isp = 1
      call dinv33(plat,1,qlat,vol)
      vol = dabs(vol)*(alat**3)
      do  8  i = 1, 3
    8 llshft(i) = lshft(i) /= 0
      call subasi(s_ctrl,s_spec,s_ham)
      call suham(s_ctrl,s_lat,s_spec,s_site,s_ham,s_pot,s_str,s_strn)
      ldham = s_ham%ldham

C ... debugging
*     call query('continue?',-1,0)
*     print *, 'zero vintra, aintra=t'
*     call dpzero(vintra,nl**2*nsp*nclass)
*     aintra = .true.

C ... Decide on size for screening matrix
      if (aintra) then
        ix = 0
        lpdim = 0
        do  10  i = 1, nbas
        do  10  l = 1, nl
          if (s_ham%iprmb(1+ix) <= ldim) lpdim = lpdim+1
          ix = ix+2*l-1
   10   continue
      else
         lpdim = nbas
      endif

C ... Big memory allocation (needs be reorganized)
      ld2 = ldim**2
      lp2 = lpdim**2
      allocate(v0(lp2))
      allocate(v01(lp2))
      allocate(v00(lp2))
      allocate(wscr(lp2))
      allocate(wscrk(lp2*nkp))
      allocate(p0(lp2*2)); call dpzero(p0,2*lp2*2)
      allocate(mad(nbas*nbas))
      allocate(lmx(nbas))
      allocate(wscrl(lp2)); call dpzero(wscrl,2*lp2)
      call icopy(nbas,nl-1,0,lmx,1)
C#ifdefC v1
C      deallocate(wscrl)
C      allocate(wscrl(lp2*nkp)); call dpzero(wscrl,2*lp2*nkp)
C      allocate(p1(lp2))
C      allocate(p0l(lp2)); call dpzero(p0l,2*lp2)
C      allocate(mad0(nbas*nbas))
C      allocate(mad1(nbas*nbas))
C      allocate(p00(lp2))
C      allocate(p01(lp2))
C      allocate(p02(lp2))
C      allocate(w0(lp2))
C      allocate(w1(lp2))
C      allocate(w2(lp2))
C      allocate(w0m(lp2))
C      allocate(w1m(lp2))
C      allocate(w2m(lp2))
C#endif
C#ifdefC v3
C#elseif v2 | v4
      allocate(pa(nbas*nbas))
      allocate(v02(lp2))
C#endif

C ... Work array to hold eigenvectors
      allocate(evec(ld2*nkp)); call dpzero(evec,2*ld2*nkp)
C ... Work arrays to make P0 and Sigma
      allocate(h(ld2))
      allocate(h2(ld2))
      allocate(h3(ld2))

      do  isp = 1, nsp

C --- Read in all the evecs for one spin for now ---
      ifi = fopna('EVEC',-1,4)
      rewind ifi
C     Read header nkp,nsp
      read(ifi,err=99,end=99) i,j
      if (i /= nkp) call rx1(' ASASX: evec file has nkp=%i',i)
      if (j /= nsp) call rx1(' ASASX: evec file has nsp=%i',j)
      do  iq = 1, nkp
C       Record: eb header
        read(ifi,err=99,end=99) i,nev
C       Record: eband
C       call dpdump(eband(1,1,iq),ldim,ifi)
        read(ifi,err=99,end=99) dum
C       Record: evec header
        read(ifi,err=99,end=99) i,nev
        if (isp == 2) then
C         Record: evecs
          read(ifi,err=99,end=99) dum
C         Record: eb header
          read(ifi,err=99,end=99) i,nev
C         Record: eband
          read(ifi,err=99,end=99) dum
C         Record: evec header
          read(ifi,err=99,end=99) i,nev
        endif
        if (i /= ldim .or. nev /= ldim) call fexit2
     .    (-1,111,' ASASX: expected nev=%i but read nev=%i',ldim,nev)
C       efermi = eband(4,1,1) + 1d-12
        call dpsdmp(evec,1+(iq-1)*ld2*2,iq*ld2*2,ifi)
        if (isp == 1 .and. nsp == 2) then
          read(ifi,err=99,end=99) i
          read(ifi,err=99,end=99) dum
          read(ifi,err=99,end=99) i
          read(ifi,err=99,end=99) dum
        endif
      enddo
C     call yprm('eband',1,eband,0,ldim,ldim,nkp)
      call fclose(ifi)

C --- Nonlocal, static unscreened response P0 for q=0 ---
      if (lp0) then
        call asxp0(lsharm,aintra,1,isp,plat,nl,nsp,nbas,s_ham%offH,
     .    s_ham%iprmb,nkp,llshft,s_bz%ipq,s_bz%star,ldim,lpdim,eband,
     .    nbmx,evec,zdel,s_lat%pos,s_lat%istab,s_lat%symgr,
     .    s_lat%ag,efermi,h,h2,h3,p0)

C   ... Save p0(q=0) and exit
        i = fopn('PSTA')
        call awrit1('%x nsp=%i',outs,80,0,nsp)
        call dscal(lpdim**2,1d0/nsp,p0,1)
        call ywrm(0,outs,1,i,'(5f15.10)',p0,
     .    lpdim*lpdim,lpdim,lpdim,lpdim)
C   ... Make psta for second spin before exiting
        if (isp < nsp) cycle
        call fclose(i)
        call tcx('asasx')
        return
      endif

      if (nsp == 2) call rx('asasx not ready for nsp=2')

C --- Nonlocal, screened potential for all q ---
      do  iq = 1, nkp

        call dpscop(s_bz%qp,qp,3,3*iq-2,1,1d0)
C   ... q-dependent bare Coulomb interaction
        call madmtq(0,qp,nbas,s_lat%pos,awald,alat,vol,
     .    s_lat%dlv,nkd,s_lat%qlv,nkr,mad,xx,xx)
        call vbare(0,nsp,aintra,vintra,mad,xx,nbas,nl,lpdim,
     .    s_ham%iprmb,s_ctrl%ipc,ldim,s_ctrl%rmax,v0,xx)
C#ifdef v2 | v3
        if (iq == 1) call dcopy(lp2*2,v0,1,v01,1)
        if (iq == 2) call dcopy(lp2*2,v0,1,v02,1)
C#endif

C   ... Nonlocal, static response P0 for non-interacting particles
        call asxp0(lsharm,aintra,iq,isp,plat,nl,nsp,nbas,s_ham%offH,
     .    s_ham%iprmb,nkp,llshft,s_bz%ipq,s_bz%star,ldim,lpdim,eband,
     .    nbmx,evec,zdel,s_lat%pos,s_lat%istab,s_lat%symgr,
     .    s_lat%ag,efermi,h,h2,h3,p0)
C       print *, 'iq=',iq
C       call yprm('P0(q)',2,p0,lpdim*lpdim,lpdim,lpdim,lpdim)

C   ... Accumulate local part of P0
C#ifdefC v1
C        call p0loc(iq,s_bz%wtkp,p0,p0l,lpdim)
C        if (iq == 2) then
C          call dcopy(lpdim*lpdim*2,p0,1,p1,1)
C        endif
C#endif
C#ifdef v2 | v4
        if (iq == 1) then
          j = 1
          if (aintra) j = 11
          call plm2pa(p0,nbas,lpdim,ldim,s_ham%iprmb,j,1,pa)
        endif
C#endif

C   ... Screened potential for this q
        call wstat(lpdim,p0,v0,wscr)
        call dpscop(wscr,wscrk,lp2*2,1,1+(iq-1)*lp2*2,1d0)

      enddo

C --- w(q = lim_x->0 q1*x) ---
      call dpscop(s_bz%qp,qp2,3,4,1,1d0)
C#ifdef v2 | v3
      call wsloc(wscrk,lpdim,v01,v02,alat**3/vol,
     .  s_bz%qp,nkp,nkxyz,s_bz%wtkp,wscrl)
      call wsloca(nsp,s_ham%iprmb,lihdim,nbas,nl,s_ctrl%ipc,s_pot%qnu,
     .  lpdim,wscrl)
C#endif

C#ifdefC v1
CC --- Local part of screened potential  ---
C      do  40  iq = 1, nkp
C        call dpscop(s_bz%qp,qp,3,3*iq-2,1,1d0)
CC   ... Remake q-dependent bare Coulomb interaction
C        call madmtq(0,qp,nbas,s_lat%pos,awald,alat,vol,
C     .       s_lat%dlv,nkd,s_lat%qlv,nkr,mad,xx,xx)
C        call vbare(0,nsp,aintra,vintra,mad,xx,nbas,nl,lpdim,
C     .    s_ham%iprmb,s_ctrl%ipc,ldim,s_ctrl%rmax,v0,xx)
CC   ... Screened response for this q, local part of P0 ---
C        call wstat(lpdim,p0l,v0,wscr)
C        call dpscop(wscr,wscrl,lp2*2,1,1+(iq-1)*lp2*2,1d0)
C   40 continue
C
CC --- w(q = lim_eps->0 q1*eps) ---
C      call madmtq(1,qp2,nbas,s_lat%pos,awald,alat,vol,
C     .  s_lat%dlv,nkd,s_lat%qlv,nkr,mad0,mad1,v2)
C      call vbare(1,nsp,aintra,vintra,mad0,mad1,nbas,nl,lpdim,
C     .  s_ham%iprmb,s_ctrl%ipc,ldim,s_ctrl%rmax,v00,v01)
Cc ... Nonlocal, static response P0 for q1
C      call asxp0(lsharm,aintra,1,isp,plat,nl,nsp,nbas,s_ham%offH,
C     .  s_ham%iprmb,nkp,llshft,s_bz%ipq,s_bz%star,ldim,lpdim,eband,nbmx,
C     .  evec,zdel,s_lat%pos,s_lat%istab,s_lat%symgr,s_lat%ag,efermi,
c     .h,h2
C     .  ,h3,p0)
CC ... local part of q->0 limit of P0
C      call p0q0(qp2,lpdim,p0,p1,p00,p01,p02)
CC ... q->0 limit of w
C      call wq0(lpdim,p00,p01,p02,v00,v01,v2,
C     .     w0,w1,w2)
CC ... q->0 limit of w from local part of P
C      call p0q0(qp2,lpdim,p0l,p0l,p00,p01,p02)
C      call wq0(lpdim,p00,p01,p02,v00,v01,v2,
C     .     w0m,w1m,w2m)
CC --- Difference w(P0) and w(local-P0) ---
C      call wdiff(lpdim,nkp,wscrk,wscrl,w0,w1,w2,
C     .           w0m,w1m,w2m)
C#endif
      enddo

C --- Make Sigma, new evals and evecs for niter iterations ---
      do  iter = 1, niter

      allocate(ecorr(ldim))
C ... Write sigma header
      ifisq = fopna('SIGM',-1,4)
*     ifiqp = fopnn('QPBS')
      rewind ifisq
C     i = iosig(lrssig,-ifisq,0,nkp,ldim,xx)
      call iosigh(3,0,nsp,1,ldim,ldim,nkxyz(1),nkxyz(2),nkxyz(3),nkp,nkp,
     .  lshft(1),lshft(2),lshft(3),-ifisq,0d0)

C#ifdef v2 | v4
      allocate(dedv(ldim)); call dpzero(dedv,ldim)
      allocate(vsxl(ldim)); call dpzero(vsxl,ldim)
C#endif

C --- Make nonlocal sigma (v2 missing local part), write to disk ---
C     ifit = fopn('TMP')
      ifit = fopnT('TMP',-1,4,0)
      rewind ifit
      do  iq = 1, nkp

C   ... Make Sigma (h2), dump to disk
C#ifdefC v1
C        call asxsig(aintra,iq,plat,nl,nsp,nbas,s_ham%offH,s_ham%iprmb,nkp,
C     .    s_bz%wtkp,llshft,s_bz%ipq,s_bz%star,ldim,lpdim,eband,nbmx,
C     .    evec,s_lat%pos,s_lat%istab,s_lat%symgr,s_lat%ag,efermi,
C     .    wscrk,w0,w1,w2,h,h3,p0,h2)
C#elseif v2 | v4 | v3
        call asxsgm(lsharm,aintra,iq,plat,nl,nbas,s_ham%offH,
     .    s_ham%iprmb,nkp,s_bz%wtkp,llshft,s_bz%ipq,s_bz%star,ldim,
     .    lpdim,eband,nbmx,evec,s_lat%pos,s_lat%istab,s_lat%symgr,
     .    s_lat%ag,efermi,wscrk,wscrl,h,h3,p0,h2)
C       call yprm('sigma',2,h2,ldim*ldim,ldim,ldim,ldim)
C   ... Accumulate linear response matrix dE_sx/dv
C#ifdef v2 | v4
        call desxdv(ldim,nsp,nbas,h2,evec,h3,
     .    s_bz%wtkp,efermi,zdel,eband,nbmx,iq,dedv)
C#endif
C  ...  debugging
C       call plm2pa(dedv,nbas,lpdim,ldim,s_ham%iprmb,20,1,h3)
C       call awrit2('dedv=%n:1,14;14d',' ',80,6,nbas,h3)
C#endif
        call dpdump(h2,ld2*2,-ifit)
      enddo

C --- Make local sigma (v2 | v4) ---
C#ifdef v2 | v4
      call plm2pa(dedv,nbas,ldim,ldim,s_ham%iprmb,20,1,vsxl)
      call awrit2('dedv=%n:1,8;8d',' ',80,6,nbas,vsxl)
      call desxdn(nbas,pa,vsxl,vsxl)
      call awrit2('vsxl=%n:1,8;8d',' ',80,6,nbas,vsxl)
C#endif

C --- Make and save sigma-sigma(loc), pert. corr of evals ---
      rewind ifit
      do  iq = 1, nkp
        call dpscop(s_bz%qp,qp,3,3*iq-2,1,1d0)
        if (ipr >= 20) then
          do  j = 1, 2
          call awrit3('  QP pert. corr. (kpt %i of %i), '//
     .      'q =%3:1,4;4d',' ',80,lgunit(j),iq,nkp,qp)
*         write(ifiqp,'(3f10.5)') (qp(j),j=1,3)
        enddo
        endif

C   ... Read Sigma, subtract local part (v2 | v4)
        call dpdump(h2,ld2*2,ifit)
C#ifdef v2 | v4
        call dvsxl(nbas,ldim,aintra,s_ham%offH,s_ham%iprmb,vsxl,h2)
C#endif

C   ... For Gamma point, symmetrize according to group ops
        if (iq == 1) then
C     ... Allocate space for sig(q=0), and hold it
          allocate(h4(ld2)); call dpzero(h4,2*ld2)
C         allocate(rmat(nl**4))
C         call yprm('Sig (q=0)',2,h2,ldim*ldim,ldim,ldim,ldim)
          do  ig = 1, nsgrp
            call dcopy(ld2*2,h2,1,h3,1)
            call dpscop(s_lat%symgr,g,9,9*ig-8,1,1d0)
            call dpscop(s_lat%ag,ag,3,3*ig-2,1,1d0)
            call roth(0,nl,nbas,s_lat%pos,xx,s_ham%offH,s_ham%iprmb,
     .        s_lat%istab(1,ig),g,ag,s_bz%qp,ldim,ldim,h,h3)
            call daxpy(ld2*2,1d0/nsgrp,h3,1,h4,1)
*           call daxpy(ld2*2,-1d0,h2,1,h3,1)
*           print *, ig, ddot(ld2*2,h3,1,h3,1)
          enddo
*         call yprm('Symm. sig (q=0)',2,h4,ldim*ldim,ldim,ldim,ldim)
          call dcopy(ld2*2,h4,1,h2,1)
          deallocate(h4)
        endif

C   ... Perturbation correction
        call asxdev(iq,h2,ldim,evec,ecorr)

C   ... Accumulate inverse Bloch transform of Sigma
        if (lrssig == 1) then
C     ... Make pair table and allocate memory for rs sigma
          if (iq == 1) then
            call pshpr(ipr-10)
            allocate(s_sts%n(nbas+1))
            i = 0
            call pairs(nbas,nbas,1d0,plat,[rsrnge*avw/alat/2],s_lat%pos,
     .        [-1],3,-1,xxia,nttab,s_sts%n,iax,i)
            s_sts%i => iax
            allocate(ntabs(nbas+1))
            call icopy(nbas+1,s_sts%n,1,ntabs,1)
            call poppr
            i = ntabs(nbas+1)*nl**4
            allocate(s_sts%s(i)); call dpzero(s_sts%s,i)
            allocate(s_sts%a(nl**2*nbas))
          endif
C          print *, 'iq=',iq
C          call yprm('Sig(iq)',2,h2,ldim*ldim,ldim,ldim,ldim)
          call rx('fix call to ibloch')
C          call ibloch(0,isp,nsp,s_bz%qp,s_bz%wtkp,iq,plat,s_ham%offH,
C     .      s_ham%iprmb,1,nttab,iax,ldim,ldim,h2,s_sts%s,nl**2)
        endif

C   ... Write sigma(q)
        if (lrssig == 0 .or. iq == 1) then
          write(ifisq) qp
C          call yprm('sig-q',2,h2,ldim*ldim,ldim,ldim,ldim)
          call dpdump(h2,ld2*2,-ifisq)
        endif
      enddo

C --- Symmetrize sigma(r) ---

      if (lrssig == 1) then
        call rx('update call to rsmsym')
        i = ntabs(nbas+1)*nl**4
C       call rx('need to disginguish sigrs,sigr')
        allocate(s_sts%s(i)); call dpzero(s_sts%s,2*i)
C        call rsmsym(0,plat,nbas,s_lat%pos,nl,nsp,1,nttab,
C     .    ntabs,iax,s_lat%symgr,s_lat%ag,nsgrp,0,1,s_sts%s,s_sts%s)
        call rx('asasx 555 not ready')
C       call defi(ontab, nbas+1)
        allocate(s_sts%n(nbas+1))
        call iostr(1,'SIGR',nl,nbas,1,0d0,0,ckbas,-1,nttab,s_sts)
C    .    oalph,oiaxs,ontab,osigrs)
        call fclose(fopn('SIGR'))

C   ... debugging printout
        call query('V>=60 to display rs sigma',-1,0)
        if (iprint() >= 60) then
          call shostr(nl,nttab,nbas,1,plat,s_lat%pos,01000,s_sts%a,
     .      iax,ntabs,s_sts%s,1,dum,1,1d0)
        endif
C    ... Debugging: see how pert corr changes when make from Bloch
C        if (iosig(lrssig,ifisq,1,nkp,ldim,h3) < 0) stop 'oops'
C        do  160  iq = 1, nkp
C          call dpscop(s_bz%qp,qp,3,3*iq-2,1,1d0)
C          call dpzero(qp0,3)
C          call dcopy(ld2*2,h3,1,h2,1)
C          call dscal(ld2*2,-1d0,h2,1)
C          call blcho(10,isp,nsp,qp0,nbas,plat,s_ham%offH,1,nttab,
C     .      iax,s_sts%s,nl**2,ldim,0,ldim,h2,xx,xx)
C          call dscal(ld2*2,-1d0,h2,1)
C          call blcho(10,isp,nsp,qp,nbas,plat,s_ham%offH,1,nttab,
C     .      iax,s_sts%s,nl**2,ldim,0,ldim,h2,xx,xx)
C          call asxdev(iq,h2,ldim,evec,ecorr)
CC          call blcho(0,isp,nsp,qp,nbas,plat,s_ham%offH,1,nttab,
CC     .      iax,s_sts%s,nl**2,ldim,0,ldim,h2,xx,xx)
C        call yprm('sigx',2,h2,ldim*ldim,ldim,ldim,ldim)
C 160   continue
      endif

C --- Band pass, with SX correction ---

C ... Setup for new band pass
      call rxx(lncol /= 0,'asasx not implemented for noncol.')
      rewind ifisq
      call iosigh(3,0,nsp,1,ldim,ldim,nkxyz(1),nkxyz(2),nkxyz(3),nkp,iq,
     .  lshft(1),lshft(2),lshft(3),ifisq,xx)

C      if (iosig(lrssig,ifisq,0,nkp,ldim,xx) < 0)
C     .  call rx('ASASX: bad sigma file')

C ... For each qp, make new evals and evecs ...
C ... Suppress combined correction until last iteration
      s_ctrl%lasa = lasa - bitand(lasa,4)
      isp = 1
      izval = zval/2
      evtop = -9999
      ecbot = -evtop
      do  iq = 1, nkp
        call dpscop(s_bz%qp,qp,3,3*iq-2,1,1d0)
C   ... Evecs in h, make without ccor except for last iter
        if (iter == niter) s_ctrl%lasa = lasa
        allocate(s(ld2))
        idim = lidim - ldim
        allocate(s_wk(1)%p(max(ldim*idim*2,1)))
        allocate(s_wk(3)%p(max(ldim*ldim*2,1)))
        if (lncol /= 0.and.idim > 0) then
          allocate(s_wk(2)%p(8*ldim**2))
        endif
        call secmat(s_ctrl,s_spec,s_ham,s_pot,s_lat,s_str,iq,
     .    nkp,qp,xx,isp,-1,ldim,99d0,nlibu,lmaxu,idu,ludiag,vorb,nev,
     .    h,s,s_wk,eband(1,isp,iq),xx,xx,xx)
        if (lncol /= 0.and.idim > 0) deallocate(s_wk(2)%p)
        deallocate(s_wk(1)%p,s_wk(3)%p)
        deallocate(s)
        call dpscop(h,evec,ld2*2,1,1+(iq-1)*ld2*2,1d0)
        evtop = max(evtop,eband(izval,isp,iq))
        ecbot = min(ecbot,eband(izval+1,isp,iq))
      enddo
      deallocate(ecorr)
      if (ipr >= 20)
     .  call awrit3(' %N Highest occ. level = %,5;5d '//
     .  ' Lowest unocc. = %,5;5d  diff =  %,5;5d',' ',80,stdo,
     .  evtop,ecbot,ecbot-evtop)
      enddo
      deallocate(v0,v01,v00,wscr,wscrk,p0,lmx,wscrl)
      deallocate(evec,h,h2,h3)
C#ifdef v2 | v4
      deallocate(pa,v02,dedv,vsxl)
C#endif

      call fclose(ifisq)
      call fclose(ifit)
      call tcx('asasx')
      return

C --- Error exit ---
   99 call rx('asasx: failed to read evecs')

      end
