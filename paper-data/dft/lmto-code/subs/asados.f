      subroutine asados(prgnam,s_bz,s_ctrl,s_ham,s_lat,s_spec,s_site)
C- Partial density-of-states generator
C ----------------------------------------------------------------------
Cio Structures
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  n w efmax lcond dosw nkabc nkp ntet lshft
Co     Stored:     nkp nkabc ntet
Co     Allocated:  qp wtkp idtet star ipq
Cio    Elts passed:idtet qp ipq star wtkp lopt lio
Cio    Passed to:  mkqp
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbas nclass nl nspec lncol lham ldos nclasp ipc lpgf
Ci                 lmet
Co     Stored:     lmet
Co     Allocated:  *
Cio    Elts passed:ics lmet dclabl ipc lsx lscr
Cio    Passed to:  subasi mkqp
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  subasi
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  avw alat plat nsgrp npgrp
Co     Stored:     npgrp
Co     Allocated:  *
Cio    Elts passed:symgr
Cio    Passed to:  mkqp
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa lmxb z
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  subasi sumlst
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  clabel spec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  sumlst
Ci Inputs
Ci   prgnam:string labelling name of caller (for printout)
Ci   Optional command-line switches:
Ci     --dos:option1:option2...
Ci       Options are:
Ci       wtfn=name  read dos weights from file `name'
Ci       cls        is equivalent to wtfn=CLS
Ci       tbdos      uses moments file with tbe conventions; also moments
Ci                  file is named 'BAND'
Ci       mode=#     makes dos-related quantities:
Ci                  #=0 makes dos
Ci                  #=1 (ballistic) conductivity,
Ci                      1/2 | grad_k E(k) . vec
Ci                      In this mode, vec must also be specified.
Ci                  #=2 (diffusive) conductivity, or actually
Ci                      grad_1 E(k) . grad_2 E(k) where directions
Ci                      1 and 2 must be specifed by vec and vec2 below
Ci       idos       generate energy-integrated dos
Ci       rdm        Write DOS in rdm-compatible format
Ci       ldig       Write DOS with extra digits
Ci       fpdos      tells asados to run in 'fp' mode.  It only
Ci                  affects printout associating dos channels with
Ci                  site and class names
Ci       npts=#     number of energy points.  If not specified with
Ci                  command-line --dos, user is prompted for input.
Ci       window=#,# energy window over which data is to be computed
Ci                  #,# are minimum, maximum energies.  If not
Ci                  specified with command-line --dos, user is
Ci                  prompted for input.
Ci       vec=#,#,#  direction vector for conductivity; see mode=#
Ci       vec2=#,#,# second direction vector for `diffusive'
Ci                  conductivity; see mode= above
Ci       totdos     compute total dos by adding weights from all
Ci                  partial dos
Ci       bands=list compute contribution to dos from a prescribed
Ci                  list of bands
Ci       classes=list generate dos for a specified list of classes.
Ci                  list syntax follows standard class-list syntax
Ci                  as described in lmto.html.
Co Outputs
Co    DOS is generated and written to file 'dos'
Cl Local variables
Cb Bugs
Cb   Sampling DOS doesn't generate spin 2 dos in a case with SO coupling
Cb   Test case: nc/test/ctrl.cdte
Cb   lm --noinv  --dos:npts=2001:window=-1.2,.5 cdte --soc -vso=t --iactiv '--pdos'
Cb   lmdos --noinv  --dos:npts=2001:window=-1.2,.5 cdte --soc -vso=t --iactiv '--pdos'
Cb   echo 8 5 -1.0 .4 | pldos -fplot -lst='1;3;5;7;9;11;13;15;17' -lst2 dos.cdte
Cu Updates
Cu   28 Oct 15  Repackaged parsing of --dos options into a subroutine
Cu   08 May 13  Complete migration to f90 structures; eliminate s_array
Cu   06 Sep 11  Begin migration to f90 structures
Cu   27 Mar 10  Eliminate bzmap (no longer used)
Cu   05 Mar 10  Better printout when channels generated through sumlst
Cu              Mulliken and pdos treated in a more consistent manner
Cu   08 Jun 06  Bug fix for noncollinear case; new --mull argument
Cu    8 Jun 05  (ATP) handles empirical TB case
Cu   13 Sep 03  Revised printout of channels
Cu   18 Jan 02  Revised printout of identification of channels
Cu              with site orbitals
Cu   07 Feb 01 *several revisions to deal with the following cases:
Cu              different modes for DOS-related quantities
Cu              integration over a subset of bands
Cu              possibility to integrate quantities for total DOS
Cu              unification of nfp,etb,asa
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      character*(*) prgnam*8
C ... For structures
!      include 'structures.h'
      type(str_bz)::    s_bz
      type(str_ctrl)::  s_ctrl
      type(str_ham)::   s_ham
      type(str_lat)::   s_lat
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
C ... Dynamically allocated local arrays
      integer, allocatable :: iblst(:)
      integer, allocatable :: chan(:)
      integer, allocatable :: lmx(:)
      integer,allocatable:: lchan(:,:)
      real(8), allocatable :: z(:)
      real(8), allocatable :: eband(:)
      real(8), allocatable :: doswt(:)
      real(8), allocatable :: dos(:)
      real(8), allocatable :: wk(:)
C ... Local variables
      character*120 outs,sopts,wtfn
      logical T,F,ltet,cmdopt,tbdos,TBU,totdos,lcond,lrange,lidos,lfp,ldig,ltmp,lkres
      integer fopnn,fopna,i,ib,ic,ichds,iclbas,icond,ifmt,il,ief,
     .  imax,iprint,is,iv(10),j,j1,j2,l,ldos,lham,lmax,lmdim,lmxch,
     .  lncol,lstyle,m,moddos,nband,nbandx,nbas,nblst,nchan,nchan2,
     .  nchds,nchmx,nclasp,nclass,ndum,nevmx,nfilem,nfstg,nkp,
     .  nkxyz(3),nl,nlf,nll,nlo,nlstc,norder,npgrp,npts,nsite,
     .  nsp,nspc,nspd,nspec,nspx,ntet,parg,scrwid,stdo
      integer, parameter :: nsitmx=1024, NULLI=-99999
      double precision ef
C     DOS switches set by sudossw
      logical doslsw(9)
      integer dosisw(7)
      double precision dosrsw(8)
      equivalence (lidos,doslsw(1)), (totdos,doslsw(2)), (lrange,doslsw(3))
      equivalence (ldig,doslsw(4)), (lfp,doslsw(5)), (tbdos,doslsw(6)), (tbu,doslsw(7))
      equivalence (npts,dosisw(1)), (ifmt,dosisw(2)), (icond,dosisw(3))

C     Site or class list
      integer lsite(nsitmx)
      double precision vmtz(2),efermi(2),efmax,emin,emax,drange(2),
     .  width,cvec(0:6),ddot,avw,alat,plat(3,3),xx
      character*40 modstr(0:4),strn*120,strn2*120
      character*8 clabl
      parameter (T=.true., F=.false., scrwid=80)

      procedure(logical) isanrg
      procedure(integer) :: nglob,mkilsd,a2vec,isw
      procedure(real(8)) :: dlength

      data modstr /'dos','sigma(ballistic)','sigma(v*v)','|v|','*'/

C --- Setup ---
      stdo = nglob('stdo')
      call lvset(doslsw,1,size(doslsw),.false.)
      call ivset(dosisw,1,size(dosisw),0)
      call dvset(dosrsw,1,size(dosrsw),0d0)
      nbas = s_ctrl%nbas
      nclass = s_ctrl%nclass
      nl = s_ctrl%nl
      nspec = s_ctrl%nspec
      lncol = s_ctrl%lncol
      lham = s_ctrl%lham
      ldos = s_ctrl%ldos
      nclasp = s_ctrl%nclasp
      avw = s_lat%avw
      alat = s_lat%alat
      plat = s_lat%plat
      norder = s_bz%n
      width = s_bz%w
      efmax = s_bz%efmax
      cvec(0:3) = s_bz%lcond(1:4)
      drange = s_bz%dosw
      norder = isign(1,norder) * mod(iabs(norder),100)

C     For now... until input in rdccat is patched
      call dcopy(3,cvec(1),1,cvec(4),1)
C ... Class-dependent arrays
      allocate(z(nclasp))
      call spec2class(s_spec,nclasp,s_ctrl%ics,'z',1,xx,z)
      allocate(lmx(nclasp))
      call spec2class(s_spec,nclasp,s_ctrl%ics,'lmxa',1,lmx,xx)
C ... Other initialization
CC    lmet = lgors('ctrl lmet,1',sctrl)
C     lmet = IAND(s_ctrl%lmet,1)
C     ltet = lgors('ctrl lmet,2',sctrl)
      ltet = IAND(s_ctrl%lmet,2) /= 0
      icond = nint(cvec(0))
      lkres  = .false.
      call subasi(s_ctrl,s_spec,s_ham)
C     lfp = .false.
      ifmt = 1
C --- Parse possible --kres switch
      i = 7; ef=0
      if (cmdopt('--kres=',i,0,strn)) then
        is = a2vec(strn,len(strn),i,4,', ',2,2,1,iv,ef)
        if (is /= 1) call rx('asados failed to parse '//trim(strn))
        lkres  = .true.
      endif

C --- Check for --pdos or --mull command-line argument ---
C     Tells lmdos that channels are oriented by site, not class.
C     It does not affect calculation, only the identification
C     of the channels printed at the end
C     ldos is not set properly when asados is called, so make assumptions
      nsite = 0
      lsite(1) = -1
      nll = nl
      if (ldos == 0) ldos = 2
      if (cmdopt('--pdos',6,0,strn).or.cmdopt('--mull',6,0,strn2)) then
        nchmx = nbas*nl**2
        allocate(chan(nchmx))
        j = 0
        if (cmdopt('--mull',6,0,strn)) j = 1
        call sumlst(j,nchmx,nbas,1,s_site,s_spec,strn(7:),moddos,
     .    nsite,lsite,lmxch,nchan,lmdim,chan,nll)
        allocate(lchan(lmdim,nsite))
        call icopy(lmdim*nsite,chan,1,lchan,1)
        deallocate(chan)
        if (moddos == 0 .or. moddos == 3) ldos = 0
        if (moddos == 1 .or. moddos == 4) ldos = 2
        if (moddos == 2 .or. moddos == 5) ldos = 4
      endif
      nlo = nl
      if (ldos == 0) nlo = 1
      if (ldos == 2) nlo = nll
      if (ldos == 4) nlo = nll*nll

C --- Parse DOS options ---
C ... Defaults
      nlstc = 0
      nchds = -1
      lstyle = 1
      wtfn = ' '
      allocate(iblst(1)); iblst(1) = 0
      nblst = 0
      npts = 0
C     If .true. calculate conductivity rather than DOS
C     lcond = .false.
      if (cmdopt('--dos',5,0,sopts)) then
        icond = 0; dosrsw(1:2) = drange(1:2); dosrsw(3:8) = cvec(1:6)
        call sudossw(sopts,0,doslsw,dosisw,dosrsw,wtfn)
        drange(1:2) = dosrsw(1:2); cvec(1:6) = dosrsw(3:8)
C       Parse bands list, if specified
        if (dosisw(5) > dosisw(4)) then
          j1 = dosisw(4); j2 = dosisw(5); i = 6
          nblst = mkilsd(sopts(j1+i:j2),-1,iblst) ! count number
          call rxx(nblst<=0,'no bands specified: '//sopts(j1:j2))
          deallocate(iblst); allocate(iblst(nblst))
          nblst = mkilsd(sopts(j1+i:j2),nblst,iblst) ! get list
        endif
C       Additional parsing for classes
        if (dosisw(7) > dosisw(6)) then
          j1 = dosisw(6); j2 = dosisw(7)
          if (nsite /= 0) call rx('ASADOS: specification of class list not compatible with --pdos')
          j1 = j1+7
          m = j1
          i = parg('style=',2,sopts,m,len(sopts),sopts(j1:j1)//' ',1,1,iv,lstyle)
          if (i < 0) then
            goto 999
          elseif (i == 0) then
            j1 = j1+1
          else
            j1 = m+2
          endif
          call clist(lstyle,sopts(j1:j2),s_ctrl%dclabl,z,nclass,nlstc,lsite)
        endif
      endif

C      if (cmdopt('--dos',5,0,sopts)) then
C      dc = sopts(6:6)
C      if (dc /= ' ') then
CC   ... Return here to resume parsing for arguments
C        j2 = 5
C   10   continue
C        j2 = j2+1
C        if (sopts(j2:j2) == dc) goto 10
C        j1 = min(len(sopts),j2)
C        call nwordg(sopts,0,dc//' ',1,j1,j2)
C        if (j2 >= j1) then
C          if (.false.) then
C          elseif (sopts(j1:j1+4) == 'wtfn=')  then
C            if (j1+5 > j2) call rx('ASADOS: bad file name')
C            wtfn = sopts(j1+5:j2)
CC         cls switch just sets file name to 'cls'
C          elseif (sopts(j1:j2) == 'cls')  then
C            wtfn = 'CLS'
C            lfp = .true.
CC         generate integrated DOS
C          elseif (sopts(j1:j2) == 'idos')  then
C            lidos = .true.
CC         Write does with extra digits
C          elseif (sopts(j1:j2) == 'ldig')  then
C            ldig = .true.
CC         tbdos has different dimensions in moments file
C          elseif (sopts(j1:j2) == 'tbdos')  then
C            tbdos = .true.
CC         tbdos has different dimensions in moments file
C          elseif (sopts(j1:j2) == 'tbu')  then
C            tbdos = .true.
C            TBU = .true.
CC         fp mode has site-dependent dos channels (for printout)
C          elseif (sopts(j1:j2) == 'fpdos')  then
C            lfp = .true.
CC         Generate total DOS (no partial decomposition into channels)
C          elseif (sopts(j1:j2) == 'totdos')  then
C            totdos = .true.
CC         Specify dos file format to be rdm-compatible
C          elseif (sopts(j1:j2) == 'rdm')  then
C            ifmt = 0
CC         Calculate conductivity or some quantity other than dos
C          elseif (sopts(j1:j1+4) == 'mode=')  then
C            m = 0
C            i = parg('mode=',2,sopts(j1:),m,len(sopts(j1:)),
C     .        dc//' ',1,1,iv,icond)
C            if (i <= 0) goto 999
C            call sanrg(.true.,icond,0,2,' ','dos mode')
CC           lcond = .true.
CC         Number of energy points
C          elseif (sopts(j1:j1+4) == 'npts=')  then
C            m = 0
C            i = parg('npts=',2,sopts(j1:),m,len(sopts(j1:)),
C     .        dc//' ',1,1,iv,npts)
C            if (i <= 0) goto 999
CC         DOS window
C          elseif (sopts(j1:j1+6) == 'window=')  then
C            m = 0
C            i = parg('window=',4,sopts(j1:),m,len(sopts(j1:)),
C     .        ', '//dc,2,2,iv,drange)
C            if (i <= 1) goto 999
C            lrange = .true.
CC         First direction vector (icond=1,2)
C          elseif (sopts(j1:j1+3) == 'vec=')  then
C            m = 0
C            i = parg('vec=',4,sopts(j1:),m,len(sopts(j1:)),
C     .        ', '//dc,2,3,iv,cvec(1))
C            if (i <= 0) goto 999
CC         Second direction vector (icond=2)
C          elseif (sopts(j1:j1+4) == 'vec2=')  then
C            m = 0
C            i = parg('vec2=',4,sopts(j1:),m,len(sopts(j1:)),
C     .        ', '//dc,2,3,iv,cvec(4))
C            if (i <= 0) goto 999
CC         DOS for a prescribed list of bands only
C          elseif (sopts(j1:j1+5) == 'bands=') then
C            if (j1+6 > j2) call rx('ASADOS: bad list')
C            call mkils0(sopts(j1+6:j2),nblst,iblst)
C            if (nblst > 100) call rx('increase size of iblst')
C            call mkilst(sopts(j1+6:j2),nblst,iblst)
CC         DOS for a prescribed list of classes only
C          elseif (sopts(j1:j1+6) == 'classes') then
C            if (nsite /= 0) call rx('ASADOS: specification of class '
C     .        //'list not compatible with --pdos')
C            j1 = j1+7
C            m = j1
C            i = parg('style=',2,sopts,m,len(sopts),sopts(j1:j1)//' ',1,1,iv,lstyle)
C            if (i < 0) then
C              goto 999
C            elseif (i == 0) then
C              j1 = j1+1
C            else
C              j1 = m+2
C            endif
C            call clist(lstyle,sopts(j1:j2),s_ctrl%dclabl,z,nclass,nlstc,lsite)
C          else
C            goto 999
C          endif
C          goto 10
C        endif
C      endif
C      endif

C --- Post-options initialization ---
C     lcond=T when contet branch is needed
      lcond = icond /= 0 .or. totdos .or. nblst /= 0 .or. lkres

C     List of irreducible qp, tetrahedra, ipq
      call mkqp(s_ctrl,s_bz,s_lat,ltet,lcond,1,0)
      nkxyz = s_bz%nkabc
      nkp = s_bz%nkp
      ntet = s_bz%ntet

      if (.not. ltet) call info2(20,1,1,' Integration by'//
     .  ' Methfessel-Paxton sampling,  order=%i  width=%d',norder,width)

C ... Open file containing bands and dos weights
      if (wtfn == ' ') then
        wtfn = 'MOMS'
        if (tbdos) wtfn = 'BAND'
      endif
      call info0(20,1,0,' ASADOS: reading weights from file '//wtfn//'%a')
      call info2(20,0,0,
     .  '%9fexpecting file to be resolved by '//
     .  '%?#n==0#site##%-1j'//
     .  '%?#(n==1|n==2)#l##%-1j'//
     .  '%?#n==4#l and m##',ldos,0)
C      nfilem = fopno(wtfn)
      nfilem = fopna(wtfn,-1,4+1)

C ... Get dimensions (nsp,nspc,nband,nfstg); allocate memory for eband
      call iomomq(nfilem,0,nlf,nsp,nspc,ndum,nband,nfstg,j,ndum,
     .  nchan,nchan,i,xx,xx,xx,xx,efermi,vmtz)
      ltmp = isanrg(nlf,nl,nl,' asados (warning)','file nl',.false.)
      nspd = nsp
C     Empirical TB stores stuff a little differently
      if (tbdos) then
        if (TBU) then
          nsp = 2
          nspc = 1
        else
        nband = nband*nsp
        nsp = 1
        nspc = 1
        endif
      endif
      nspx = nsp / nspc
      nbandx = nband*nspc
      allocate(eband(nband*nsp*nkp))

C ... Determine number of DOS channels (nchan=0 -> total dos only)
      read (nfilem) nchan,ndum
      call info2(20,0,0,'%9ffile has %i channel(s)%?#n==0# => generate total DOS##',nchan,nchan)
      if (nchan /= 0) then
        if (mod(nfstg/10,10) == 0) call rx(prgnam//'%a:  moments file missing DOS wts')
      else
        totdos = .true.
      endif
      if (nlstc /= 0) then
        nchds = nlstc*nl
      else
        nchds = nchan
      endif
C     No decomposition of dos; nl is no longer relevant
      if (totdos) then
        nchan = 1
        nchds = 1
        nl = nlf
      endif

C ... Read npts, emin, emax if not already available
      emin = drange(1)
      emax = drange(2)
      if (npts == 0 .or. .not. lrange) then
        if (npts == 0) npts = 501
        call info5(2,0,-1,' Enter npts (def=%i), emin and emax (def='//
     .    '%;3d,%;3d): ',npts,emin,emax,0,0)
        read (*,*) npts, emin, emax
      else
        call info5(2,0,0,'%9fUsing npts=%i  emin=%;3d  emax=%;3d',npts,emin,emax,0,0)
      endif
      if (emin == emax) then
        call rx('no energy window: emin=emax')
      elseif (emin > emax) then
        call info0(2,0,0,'%9fwarning! emin > emax ... swapping limits')
      endif
C     Adjust mesh so one point coincides with ef
      if (lkres) then
C       Point of closest approach
        xx = (npts-1)*(ef-emin)/(emax-emin)
        xx = (xx - nint(xx))*(emax-emin)/(npts-1)
        emax = emax + xx
        emin = emin + xx
        ief = nint((npts-1)*(ef-emin)/(emax-emin))
        call info(30,0,0,'%9fadjust emin,emax: point %i coincides with ef=%;6d',ief,ef)
      endif

C ... Allocate arrays for dos and moms if tetrahedron
      if (ltet) then
        allocate(doswt(nchds*nbandx*nsp*(nkp+1)))
        call dpzero(doswt,nchds*nbandx*nsp*(nkp+1))
        nfstg = 11
        if (nspc == 2) nfstg = 21
C       Suppress projection of unity into parts for total dos
        if (totdos) then
          nfstg = 1
          call dvset(doswt,1,nchds*nbandx*nsp*(nkp+1),1d0)
        endif
      else
        allocate(doswt(nchds*nbandx*nspc))
        nfstg = 1
      endif
      allocate(dos(npts*nsp*nchds),wk(npts))

C ... Make sure missing bands are high enough
      call dvset(eband,1,nband*nsp*nkp,9d9)

C ... Ensure nl,nsp,nspc,nkp,nband match; read in bands (& moms if tet.)
      if (nchds == nchan) then
        call iomomq(nfilem,32,nl,nsp,nspc,nkp,nband,nfstg,j,nband,nchan,
     .    nchan,nevmx,eband,doswt,doswt,xx,efermi,vmtz)
      else
        call iomomx(nfilem,32,nl,nsp,nspc,nkp,nband,nfstg,j,nband,nchan,
     .    nchan2,nchds,nevmx,nlstc,lsite,eband,xx,doswt,xx,efermi,vmtz)
      endif
C     call dfdump(eband,nband*nsp*nkp,-66)
      if (j /= nkp) call rx(prgnam//': moments file missing qpts')
      if (efmax-.2d0 < emax) call info2(10,1,0,' ASADOS'//
     .  '%a (warning) efmax (=%;3d) is not >> emax (=%;3d)',efmax,emax)

C ... Printout switches
      call awrit4('%N ASADOS:  make dos for %i points from %i bands in'
     .  //' window (%;3d,%;3d)',' ',80,stdo,npts,nevmx,emin,emax)
      outs = ' '
      call awrit8(' options '//
     .  '%?#n# mode: '//modstr(icond)//'%a;##'//
     .  '%?#(n==1|n==2)#%b,v=(%3;3d);#%j#'//
     .  '%?#n# kres@ef=%d#%j#'//
     .  '%?#n# totdos;##'//
     .  '%?#n# idos;##'//
     .  '%?#n# tbdos;##'//
     .  '%b %b',outs,80,0,
     .  icond,icond,cvec(1),isw(lkres),ef,isw(totdos),lidos,tbdos)
      if (outs /= ' options') call awrit0('%a',outs,-80,-stdo)
      if (nlstc /= 0) call awrit2('%10fclasses: %n:1i',' ',
     .  scrwid,stdo,nlstc,lsite)
      if (nblst /= 0) call awrit2('%10fbands: %n:1i',' ',
     .  scrwid,stdo,nblst,iblst)

C ... Sanity checks
      if (nlstc /= 0 .and. totdos) call rx(' class list is incompatible with totdos')
      do  i = 1, nblst
        if (iblst(i) > nband .or. iblst(i) <= 0)
     .    call fexit2(-1,111,' Exit -1  band %i '
     .    //'exceeds number of available bands (%i)',iblst(i),nband)
      enddo
C ... cvec(2) reverts to cvec(1) if it wasn't defined.
      if (ddot(3,cvec(4),1,cvec(4),1) == 0)
     .  call dcopy(3,cvec(1),1,cvec(4),1)
C ... Normalize cvec(1),cvec(2)
      if (icond == 1 .or. icond == 2) then
        if (ddot(3,cvec(1),1,cvec(1),1) == 0) call
     .    rx(' direction vector required for this mode but not defined')
        call dscal(3,1/sqrt(ddot(3,cvec(1),1,cvec(1),1)),cvec(1),1)
        call dscal(3,1/sqrt(ddot(3,cvec(4),1,cvec(4),1)),cvec(4),1)
      endif

C --- Ballistic conductivity ---
      if (lcond) then
        call rxx(nchds /= nchan,'cond not set for subchannels')
        if (.not. ltet) call
     .    rx('LMDOS: conductivity only with tetrahedron integration')
C   ... To save memory, allocate dos later
        deallocate(dos)
C   ... for debugging
C         nkp = s_bz%nkp
C        iblst(1) = 0
C        nblst = 1
C        nevmx = 1
C        nband = 1
C        nchan = 1
C        nbandx = 1
C        nsp = 1
C        nspx = 1
C        call snot(s_bz%qp,nkp,eband,plat)
C   ... Remake qp with reduced symmetry; assume no grp ops for now
        npgrp = 1
        s_lat%npgrp = npgrp
        call info2(2,0,0,'%10fremake qp with %i symop%?#n>1#s##',npgrp,npgrp)
        call mkqp(s_ctrl,s_bz,s_lat,ltet,F,0,2)
        nkp = s_bz%nkp
        nkxyz = s_bz%nkabc
        ntet = s_bz%ntet
        allocate(dos(npts*nsp*nchan))
        if (nblst == 0) nblst = nevmx
        if (lkres) then
          icond = icond+10
        endif

        call contet(icond,nbandx,nsp,nspx,nblst,nchan,nkxyz(1),
     .    nkxyz(2),nkxyz(3),ntet,s_bz%idtet,alat,plat,s_bz%qp,s_bz%ipq,s_bz%star,
     .    iblst,eband,cvec(1),doswt,npts,emin,emax,lidos,ief,dos)

C --- make DOS or NOS ---
      else
        if (ltet) then
          call dostet(nbandx,nsp,nspx,nevmx,nchds,nkxyz(1),nkxyz(2),nkxyz(3),
     .      ntet,s_bz%idtet,eband,doswt,npts,emin,emax,lidos,wk,dos)

C   ... dosspl reads bands and moms on the fly to save work space
        else
          call rxx(nchds /= nchan,'sampling not set for subchannels')
          call dosspl(nfilem,nbandx,nsp,nspc,nchan,norder,width,nkp,
     .      s_bz%wtkp,eband,doswt,npts,emin,emax,lidos,wk,dos)
        endif
        if (nspd > nsp) call dscal(npts*nsp*nchan,.5d0,dos,1)
      endif
      i = nsp
      if (lncol > 0) i = 1
C     rewind fopnn('DOS')
      j = 100*isw(ldig) + 3
C     call iodos(j,-fopnn('DOS'),dos,npts*(1+nsp-i),nchds,npts,nchds,emin,emax,i,efermi,ifmt)
      call iodos(j,-fopnn('DOS'),dos,npts,nchds,npts,nchds,emin,emax,nsp,efermi,ifmt) ! Change May 2017

      deallocate(doswt,dos)

C --- List channels in DOS file ---
      if (iprint() >= 30 .and. .not. totdos .and. .not. cmdopt('--cls',5,0,strn)) then
C      if (nsp == 1) write(stdo,1)
C      if (nsp == 2) write(stdo,1) 'spin-2'
C    1 format(/' Channels in dos file generated by LMDOS:'
C     .       /' site class label   spin-1':23x,a)
      call info2(1,1,0,' %?#n==0#Columns#Channels# in dos file generated by LMDOS:%N'//
     .  ' site class label   spin-1%?#n==2#%23fspin-2##',ifmt,nsp)

      if (nsite /= 0) then
        imax = nbas
      else
        imax = nclass
      endif

C ... If the channel list made by sumlst, use it
      deallocate(iblst); allocate(iblst(1000))
      if (allocated(lchan)) then
        do  i = 1, nsite
          ib = lsite(i)
          clabl = s_site(ib)%clabel
          ic = s_ctrl%ipc(ib)
          outs = ' '
            write (outs,2) ib,ic,clabl
          j = 0
          do  l = 1, lmdim
            iblst(l) = lchan(l,i)
            if (iblst(l) /= 0) then
              j = l
              if (nsp == 2) iblst(l) = 2*lchan(l,i)-1
            endif
          enddo
          if (ifmt == 0) iblst(1:j) = iblst(1:j)+1
          call ilst2a(iblst,j,outs(21:))
          if (nsp == 2) then
            j = 0
            do  l = 1, lmdim
              iblst(l) = lchan(l,i)
              if (iblst(l) /= 0) then
                j = l
                iblst(l) = 2*lchan(l,i)
              endif
            enddo
            if (ifmt == 0) iblst(1:j) = iblst(1:j)+1
            call ilst2a(iblst,j,outs(50:))
          endif
          call awrit0('%a',outs,80,-stdo)
        enddo
      else
        il = 1
        ichds = 0
C       Loop over all sites (or classes)
        do i = 1, imax
          if (lsite(1) /= -1) then
            if (i /= lsite(il)) cycle
          endif
          if (nsite /= 0) then
            ib = i
            ic = s_ctrl%ipc(ib)
          else
            ic = i
            ib = iclbas(i,s_ctrl%ipc,size(s_ctrl%ipc))
            if (ib == 0) then
              call info2(10,0,0,'no site corresponding to class',ic,0)
              cycle
            endif
          endif

C         Make a new nlo if channels are site-dependent
          if (lfp) then
C           Global upper bound to lmax
            nl = nll
            is = s_site(ib)%spec
            lmax = s_spec(is)%lmxa
            if (cmdopt('--mull',6,0,strn)) lmax = s_spec(is)%lmxb
            nll = min(lmax+1,nl)
            nlo = nll
            if (ldos == 0) nlo = 1
            if (ldos == 2) nlo = nll
            if (ldos == 4) nlo = nll*nll
          endif

C         call spacks(0,'site clabel',ssite,clabl,ib,ib)
          clabl = s_site(ib)%clabel
          outs = ' '
          write (outs,2) ib,ic,clabl
          j = 0
          do l = 0, nlo-1
            j = j+1
            iblst(j) = ichds+nsp*l+1
          enddo
          if (ifmt == 0) iblst(1:j) = iblst(1:j)+1
          call ilst2a(iblst,j,outs(21:))
C          do  44  l = 0, nlo-1
C   44     call awrit1('%a%i,',outs(21:),80,0,ichds+nsp*l+1)
C          call awrit0('%a%b ',outs(21:),80,0)
          if (nsp == 2) then
            j = 0
            do l = 0,nlo-1
              j = j+1
              iblst(j) = ichds+nsp*l+2
            enddo
            if (ifmt == 0) iblst(1:j) = iblst(1:j)+1
            call ilst2a(iblst,j,outs(50:))
C            do  46  l = 0, nlo-1
C   46       call awrit1('%a%i,',outs(50:),80,0,ichds+nsp*l+2)
C            call awrit0('%a%b ',outs(50:),80,0)
          endif
          call awrit0('%a',outs,80,-stdo)
          il = il+1

          ichds = ichds + nlo*nsp
        enddo
      endif
      endif
      return
    2 format(2i5,3x,a8)

  999 continue
      call rxs('asados: failed to parse dos options:',sopts)
      end
C      subroutine snot(qp,nkp,eband,plat)
C      implicit none
C      integer nkp
C      double precision qp(3,nkp),eband(nkp),plat(3,3)
C      integer iq
C      double precision q(3),qlat(3,3),xx
C
C      call dinv33(plat,1,qlat,xx)
C      do  iq = 1, nkp
C        call dcopy(3,qp(1,iq),1,q,1)
C        call shorbz(qp(1,iq),q,qlat,plat)
C        eband(iq) = q(1)**2+q(2)**2+q(3)**2
CC        print 333, iq, qp(1,iq),qp(2,iq),qp(3,iq), q, eband(iq)
CC  333   format(i5,6f12.6,2x,f12.6)
C      enddo
C      end

