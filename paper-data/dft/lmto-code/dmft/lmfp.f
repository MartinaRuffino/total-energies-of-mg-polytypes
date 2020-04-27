C#define DMFT
      subroutine lmfp(prgnam,s_bz,s_ctrl,s_ham,s_lat,s_mix,
     .  s_move,s_pot,s_str,s_spec,s_site,s_optic,s_gw,s_dmft,s_strn)
C- LM-FP self-consistency loop
C ----------------------------------------------------------------------
Cio Structures
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  lmet nevmx nkabc nkp ef def w n dos idtet ipq pdos
Ci                 qp star wtkp wtkb swtk ntet efmax fsmom ndos dosw
Ci                 numq lshft egap nef sopertp range
Co     Stored:     nevmx ef def w n dos idtet ipq pdos qp star wtkp
Co                 wtkb swtk nef egap numq ndos dosw sopertp nkp nkabc
Co                 ntet
Co     Allocated:  qp wtkb swtk wtkp idtet star ipq
Cio    Elts passed: wtkp qp ipq lio sopertp nef numq swtk wtkb idtet
Cio                egap n w def lopt star
Cio    Passed to:  popted rsedit iorsf bcast_strx iinit mpibc1 chimedit
Cio                bndfp rdsigm subzi addrbl sosite optinq optin2
Cio                optint mkqp
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  lrs nbas nl nspec nspin npadl npadr lpgf lcgf lscr
Ci                 lsx zbak maxit lrel lfrce nitmv mdprm ltb plbnd quit
Ci                 tol nvario defm nbasp lham lgen3 lncol ldlm nclasp
Ci                 ipc cllst clp clssl group initc ics idcc ipcp ips
Ci                 mxcst nrc ncomp pgfsl pgfvl pgord pgplp spid dclabl
Ci                 pos rmax pfloat ldos lmet lwsig lfp loptc
Co     Stored:     lfrce mdprm lrs maxit lwsig ltb cllst clp clssl
Co                 group initc ics idcc ipc ipcp ips mxcst nrc ncomp
Co                 pgfsl pgfvl pgord pgplp spid dclabl pos rmax lmet
Co     Allocated:  *
Cio    Elts passed: lrs lcd lbas lham ips lbxc nvario lmet ldos lfp
Cio                lgen3 ipc ncomp idcc spid lncol maxmem lwsig lsx lscr
Cio    Passed to:  popted rlxstp supot suham pp2alp rsedit iorsf bcast_
Cio                strx iinit mpibc1 lattic chimedit smshft bndfp
Cio                grp2av rdsigm siged suham2 fpopm optinq vcdmel rxes
Cio                dfrce relax mkqp
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  lsig ndham ehf ehk seref eterms pwmode pwemin pwemax
Ci                 npwpad lncol neula qss nbf lrsa eula elind ldham
Ci                 oveps ovncut evals nqsig eseavr sigp pmin pmax
Ci                 pnudef rsrnge rsstol nprs qsig offH nlmto ndhrs
Co     Stored:     nlibu lmaxu nlmto ndham ndofH ldham lmxax npwmin
Co                 npwpad hord lsig eterms sigp ehf ehk eula ndhrs
Co                 eseavr nqsig iaxs
Co     Allocated:  offH iprmb bdots nprs iaxs hrs qsig
Cio    Elts passed: seref offH iprmb eula magf bdots lrsa lncol nlmto
Cio                qsig evals nprs iaxs hrs rsrnge
Cio    Passed to:  suham rsedit rdovfa ovlpfa ovlocr smshft bndfp mkpot
Cio                locpot rdsigm hft2rs siged suham2 hambls hambl
Cio                augmbl sopert3 mkhso blsig makusq addrbl sosite
Cio                mkrout mkehkf mkekin chkdmu
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  plat gam plat0 alat nsgrp qlat gmax vol awald nkd
Ci                 nkq nabc ng kv kv2 igv igv2 tolft tol ag bgv cg cy
Ci                 dlv gv gvq indxcg ips0 istab jcg pos qlv symgr s_sym
Ci                 afmt npgrp napw as rpad nkdmx nkqmx platl platr
Ci                 ldist dist
Co     Stored:     gam ldist plat dist gmax nabc ng igv igv2 alat ag
Co                 bgv cg cy dlv gv gvq indxcg ips0 istab jcg kv kv2
Co                 pos qlv symgr s_sym napw vol plat0 qlat platl platr
Co                 awald nkd nkq
Co     Allocated:  ips0 bgv gv kv igv kv2 igv2 dlv qlv
Cio    Elts passed: pos symgr istab qlat alat ag dlv qlv nsgrp gv ips0
Cio                bgv igv igv2 gmax kv cg indxcg jcg cy vol plat napw
Cio    Passed to:  popted lmfopb supot sugvec0 sugvec suham rsedit
Cio                rdovfa ovlpfa ovlocr hxpbl ghibl hklbl gklbl iorsf
Cio                bcast_strx iinit mpibc1 prsed2 prsed4 chimedit
Cio                smshft pvsms1 rhgcmp symsmr bndfp rhopos mkpot smves
Cio                vesft vesgcm mshvmt symvvl ugcomp ggugbl gfigbl
Cio                fklbl hhugbl hhigbl phhigb hsmbl hgugbl smvxc2
Cio                vxcnlm smvxcm smcorm locpot rdsigm siged suham2
Cio                hambls hambl augmbl bstrux hxpgbl ghigbl hklgbl
Cio                smhsbl hhibl phhibl hsibq mkhso makusq pusq1 fpopm
Cio                addrbl fsmbl rsibl rlocbl sosite vcdmel rxes mkrout
Cio                symrat symrho dfrce pvdf4 pvdf2 pvdf1 mkekin totfrc
Cio                mixrho ioden lattic mkqp
Cio  s_mix :contains parameters for charge mixing; see structures.h
Ci     Elts read:  b bv w wc nsave mmix umix tolu
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
C
Cio  s_move :struct for dynamics parameters; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  rho vmtz nlml nlma ves aamom bxc cp ddpf dddpf ddpfr
Ci                 dlmwt dmatk dpf dpfr gibbs gma gmar grrme mad mxy
Ci                 palp papg pf pfnc pfr pmpol pnu pp ppn pprel pti qc
Ci                 qcorr qnu qnur qpp qt rhat rnew rhos rhrmx socscl
Ci                 sop shfac thetcl vdif vintr vrmax vshft smpot smrho
Ci                 smrout GFr bxcscali
Co     Stored:     bxcscali nlma nlml ves aamom bxc cp ddpf dddpf ddpfr
Co                 dlmwt dmatk dpf dpfr gibbs gma gmar grrme mad mxy
Co                 palp papg pf pfnc pfr pmpol pnu pp ppn pprel pti qc
Co                 qcorr qnu qnur qpp qt rhat rnew rhos rhrmx socscl
Co                 sop shfac thetcl vdif vintr vrmax vshft smpot smrho
Co                 smrout GFr
Co     Allocated:  shfac mad smrho smpot rhat rnew pti smrout
Cio    Elts passed: ppn shfac bxcscali rhat smrho mad pp pti smpot rnew
Cio                smrout
Cio    Passed to:  supot suham rsedit rdovfa iorsf bcast_strx iinit
Cio                mpibc1 chimedit smshft bndfp mkpot locpot suham2
Cio                mixrho
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  *
Co     Stored:     nkaps nitab nttab lmaxw nds
Co     Allocated:  iax npr alp s adot sdot
Cio    Elts passed:kaps alph
Cio    Passed to:  suham rdstrx pp2alp
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  name shfac lmxb z rmt lmxa pz orbp idxdn p norp
Ci                 ntorb rsma kmxt lmxl rg rsmv kmxv lfoca rfoca ncomp
Ci                 ngcut idu uh jh a nr qc rhoc coreh coreq ctail etail
Ci                 stc idmod nxi exi chfa rsmfa grp2 rs3 eh3 vmtz mxcst
Co     Stored:     orbp norp ntorb p pz idxdn lmxb ngcut lmxa a nr qc
Co                 nxi exi chfa rsmfa ctail etail stc rhoc name rmt z
Co                 lmxl kmxt lfoca rfoca coreh pb1 pb2
Co     Allocated:  rhoc
Cio    Elts passed:pz name idu rhoc z
Cio    Passed to:  lmfopb uspecb lmfop2 ioorbp praugm suham atfold
Cio                makidx nscpa showbs sugcut suldau sudmtu praldm
Cio                rotycs symdmu rsedit chimedit rdovfa bcast_strx
Cio                dfratm gtpcor ovlocr corprm adbkql iorsf pvsms2
Cio                smshft pvsms1 rhgcmp rhogkl bndfp grp2av pgrp2a
Cio                prrhat rhopos siterhopos dfaugm lkdim mkpot rhomom
Cio                smves vesgcm mshvmt symvvl ugcomp smvxcm smcorm
Cio                smvxc4 elocp locpot rdsigm siged dfqkkl suham2
Cio                suclst surho sumlst hambls hambl augmbl bstrux
Cio                smhsbl hsibq tbhsi hsubblock hsibq2 hsibq4 mkhso
Cio                makusq pusq1 mkorbm mullmf mkpdos fpopm mkdmtu
Cio                addrbl fsmbl fsmbpw rsibl rsibl1 rlocbl sosite
Cio                psymrq1 iorbtm mshn3p mchan vcdmel rxes mkrout
Cio                symrat symrho pnunew dfrce pvdf4 pvdf2 pvdf1 mkekin
Cio                mixrho ftlxp pvmix5 pvmix3 pvmix7 chkdmu ioden
Cio                iosits relax
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos pos0 relax clabel rho1 rho2 rhoc rho1x
Ci                 rho2x rhocx class force vel pnu pz v0 v1 bxc cpawt
Ci                 omg omgn domg dmat gc gcu gcorr gii sfvrtx j0 pdos
Ci                 qhhl qhkl qkkl eqhhl eqhkl eqkkl sighh sighk sigkk
Ci                 tauhh tauhk taukk pihh pihk pikk sohh sohk sokk
Ci                 sighhx sighkx sigkkx tauhhx tauhkx taukkx pihhx
Ci                 pihkx pikkx thet eula pl vshft ndelta delta
Co     Stored:     pos pnu pz norb vel rhoc pos0 force spec clabel bxc
Co                 cpawt omg omgn domg dmat gc gcu gcorr gii sfvrtx j0
Co                 pdos rho1 rho2 rho1x rho2x rhocx qhhl qhkl qkkl
Co                 eqhhl eqhkl eqkkl sighh sighk sigkk tauhh tauhk
Co                 taukk pihh pihk pikk sohh sohk sokk sighhx sighkx
Co                 sigkkx tauhhx tauhkx taukkx pihhx pihkx pikkx thet
Co                 v0 v1 saxis eula pl relax vshft
Co     Allocated:  v0 v1 rho1 rho2 rhoc rho1x rho2x rhocx sigkk taukk
Co                 sigkkx taukkx pikk sokk pikkx sighk tauhk sighkx
Co                 tauhkx pihk sohk pihkx sighh tauhh sighhx tauhhx
Co                 pihh sohh pihhx qkkl eqkkl qhkl eqhkl qhhl eqhhl
Cio    Elts passed: rhoc rho1 rho2 rho1x rho2x rhocx v0 qkkl qhkl qhhl
Cio                eqkkl eqhkl eqhhl
Cio    Passed to:  rlxstp suham setnorb showbs pvioeu suldau sudmtu
Cio                praldm rotycs symdmu rsedit lattic ovlpfa chimedit
Cio                rdovfa dfratm ovlocr adbkql iorsf pvsms2 bcast_strx
Cio                smshft pvsms1 rhgcmp rhogkl iopos bndfp grp2av
Cio                pgrp2a prrhat rhopos siterhopos dfaugm lkdim mkpot
Cio                rhomom smves vesgcm mshvmt symvvl ugcomp smvxcm
Cio                smcorm smvxc4 elocp locpot dfqkkl suham2 suclst
Cio                surho sumlst hambls hambl augmbl bstrux smhsbl hsibq
Cio                hsubblock hsibq2 hsibq4 mkhso makusq pusq1 mkorbm
Cio                mullmf mkpdos mkdmtu addrbl fsmbl fsmbpw rsibl
Cio                rsibl1 rlocbl sosite psymrq1 mshn3p mchan vcdmel
Cio                rxes mkrout symrat symrho pnunew dfrce pvdf4 pvdf2
Cio                pvdf1 mkekin totfrc mixrho ftlxp pvmix5 pvmix3
Cio                pvmix7 chkdmu ioden iosits cppos relax
Cio  s_optic:struct for parameters for optical ME; see structures.h
Ci     Elts read:  loptic kmxax ltet ocrng unrng dw window mefac nlg
Ci                 optme lpart esciss iq
Co     Stored:     nlg ocrng unrng nfilm nempm optme
Co     Allocated:  *
Cio    Elts passed:rgrad rgrade kmxax esmr loptic ocrng unrng optme
Cio    Passed to:  bndfp mkpot locpot ocnock opt_nstate fpopm addrbl
Cio                rsibl optinq optin2 optint
Cio  s_gw  :struct for gw-related parameters
Ci     Elts read:  nkabc eseavr lshft
Co     Stored:     eseavr
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  bndfp rdsigm siged
Cio  s_strn :struct for global strings; see structures.h
Ci     Elts read:  strn
Co     Stored:     strn
Co     Allocated:  *
Cio    Elts passed:strn
Cio    Passed to:  str_pack suham bndfp suham2
Co Outputs
Cs Command-line switches
Cs   --band        : Tabulate energy bands; see doc/Command-line-options.html
Cs   --bonly       : Use in conjuction with --lxb switch
Cs   --chimedit    : Invoke magnetic susceptibility editor
Cs   --chklim      : Check band limits in optics calculations
Cs   --cls         : For core-level spectroscopy
Cs   --cv:         : Calculate electronic specific heat, eV
Cs   --cvK:        : Calculate electronic specific heat, Ry
Cs    -ef=         : Override file Fermi level; use with --band
Cs   --etot        : 1-shot total energy mode.  No forces, dynamics, output density
Cs   --evec        : Writes eigenvectors to disk
Cs   --eveca       : Writes eigenvectors to disk in ASCII form
Cs   --fixpos      : Adjust positions slightly, rendering them consistent with the symmetry group.
Cs                 : Has effect when used in lattice relaxation
Cs   --grfac=      : Used in line minimization algorithms that relax the lattice
Cs   --jdosw       : Channels for optical properties; See doc/optics.html
Cs   --jdosw2      : Channels for optical properties; See doc/optics.html
Cs   --lxb         : Check band limits in optics calculations
Cs   --mixsig=     : For self-energy; see Command-line-options.html
Cs   --mlog        : (MPI) write MPI commands to log file
Cs   --mull        : Mulliken analysis; see doc/Command-line-options.html
Cs   --no-fixef0   : Do not adjust estimate of Fermi level after 1st band pass
Cs   --nosymdm     : Not documented
Cs   --oldbz       : Not documented
Cs   --oldvc       : Reference potential defined so average cell potential is zero
Cs   --onesp       : Generate bands for one spin.  Better: use --band~spin1...
Cs   --opt:read    : Read optical matrix elements; See doc/optics.html
Cs   --opt:write   : Write optical matrix elements; See doc/optics.html
Cs   --optbas      : Special branch to optimise basis
Cs   --pdiag       : (MPI) parallel diagonaliser
Cs   --pdos        : Partial DOS; see doc/Command-line-options.html
Cs   --phase       : Write to disk lattice shifts caused by iorsf or shorps
Cs   --popted      : Invoke the optics editor
Cs   --quit=XX     : quit after execution of certain blocks, e.g. --quit=ham
Cs   --rsedit      : Invoke the restart file editor
Cs   --rsig        : For reading the self-energy; see Command-line-options.html
Cs   --rxes        : RXES spectroscopy
Cs   --shorps      : Shorten basis vectors after final input (e.g. from iorsf)
Cs   --shorten=    : Suppress all shortening of basis vectors
Cs   --symsig      : Symmetrize sigma (if read) overriding default behavior
Cs   --wden        : Write density to file
Cs   --wforce=     :
Cs   --window=     : Generate output density in a restricted energy window
Cs   --wpotmt      : Write augmentation potentials to atom files
Cs   --wrhoat      : Write partial atomic densities to atom files
Cs   --wrhomt      : Write augmentation densities to files
Cs   --wratrho     : Write sphere rhov, rhoc, v0 to atomic file atmx
Cs   --rdatrho     : Read sphere rhov, rhoc, v0 from atomic file atmx
Cs                 : Options --rdatrho~rhov~rhoc~v0
Cs   --wsig        : For writing the self-energy; see Command-line-options.html
Cs   --wsmpot      : Write smooth potential to file
Cs   -ef=          : Overwrite Fermi level; use with --band
Cl Local variables
Cl   lmaxu :max l for a U (used for dimensioning)
Cl   nlibu :total number of U blocks
Cl   irlxsh:counter for shear relaxations.  irlxsh=0 => new step
Cr Remarks
Cr  This is the entry point for the self-consistent FP band program
Cu Updates
Cu   21 Sep 17 Extra tolerance added to nwit (correction to HF forces)
Cu   23 Dec 16 New switches --wratrho and --rdatrho
Cu   07 Oct 15 Bug fix in handling --wsite and --wpos switches
Cu   19 Nov 14 New option to write rst file only when RMS DQ has dropped
Cu   24 Jun 13 shfac is created, and read into s_pot%shfac
Cu   08 May 13 Complete migration to f90 structures; eliminate s_array
Cu   12 Nov 12 migrated structures to f90 pointers
Cu   30 Aug 12 Modifications to enable call from routine other than main
Cu   06 Sep 11 Started migration to f90 structures
Cu   10 May 09 Some improvements to basis optimization
Cu   05 Jul 08 Setup for new PW addition to basis
Cu   04 Jul 08 New restart file editor
Cu   20 Jun 06 Repackaged MPI
Cu   21 Mar 06 First cut at shear relaxations
Cu   08 Mar 06 Relaxation restores pos at minimum g when not convgd
Cu   08 Jan 06 can write to restart file rst.#, where # = iter-no.
Cu   09 Nov 05 (wrl) Convert dmat to complex form
Cu   07 Jul 05 rst file version 1.04
Cu   27 Apr 05 LDA+U added (Lambrecht)
Cu   26 Mar 05 Added switch --shorten=no to suppress pos shortening
Cu   23 Feb 05 Bug fix: forces correspondence betw/ pos and site->pos
Cu             after file read of positions.
Cu   11 Jan 05 energy convergence set to ehk when sigma included
Cu   21 Dec 04 Add option to rotate local density on file read
Cu             and to shorten basis vectors after file read
Cu   06 Sep 03 1st cut at automatic optimization of wave function
Cu    9 Jan 03 Undoes lattice shear when writing basis vectors
Cu   21 May 02 Writes restart file after smshft when moving atoms
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   08 Jun 01 Revised call to nwit
Cu   15 Feb 01 added density I/O; arguments to bndfp changed.
Cu   17 Jun 00 alpha version.  No relaxations yet.
C ----------------------------------------------------------------------
C#ifdefC LMFRS
C      use mod_str, only : str_b,str_bp,str_g,str_x
C      use mod_iface
C#endif
      use structures
      implicit none
C ... Passed parameters:
      character prgnam*8
C ... For structures
!      include 'structures.h'
      type(str_bz)::    s_bz
      type(str_ctrl)::  s_ctrl
      type(str_ham)::   s_ham
      type(str_lat)::   s_lat
      type(str_mix)::   s_mix
      type(str_move)::  s_move
      type(str_pot)::   s_pot
      type(str_str)::   s_str
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
      type(str_optic):: s_optic
      type(str_gw)::    s_gw
      type(str_dmft)::  s_dmft
      type(str_strn) :: s_strn(*)

C#ifdefC LMFRS
Cc     include 'par.h'
Cc     include 'str.h'
C      type(str_g) :: sstr_gen
C      type(str_b), allocatable, save :: sstr_bas(:)
C      type(str_bp), pointer          :: sstr_mesh(:)
C      type(str_x), pointer           :: sstr_x
C!     type(str_s), allocatable, save :: sstr_spec(:)
C#endif
C ... Dynamical arrays
      integer, allocatable :: indrlx(:,:)
      real(8), pointer :: pos2(:),wk(:)
      real(8), allocatable :: frc(:),grwk(:),hessw(:)
      real(8), allocatable :: shfac(:,:)
      character, allocatable:: slabl(:)*8
C ... Local variables
      integer procid,master,mpipid,nproc
      logical mlog
      integer i,ifi,ipr,irs(5),iscr,ix(5),j,
     .  k,lcgf,leks,lfrce,lgunit,lpgf,lpnu,lrel,lrout,lsx,nbas,nat,
     .  nbaspp,ndham,nevmx,nglob,nit1,nl,npadl,npadr,nsp,nspec,numq,
     .  stdo,pdim,lsc
      double precision plat(3,3),qlat(3,3),qbg,xv(10),fptol,umix,xx
      character strn*120, fileid*68, alabl*8, flg*3
      integer, parameter :: NULLI=-99999
      procedure(logical) :: cmdopt,bittst
      procedure(integer) :: fopn,iorsf,iosits,isw,str_pack,fopna,fopng,fxst
      procedure(real(8)) :: dlength,dsum
C For mixing.  Default parameters dmxp:
C 1(I) mixing scheme; 2 beta; 3 wc; 4,5 wgt; 6(I)nsave 7(I) mmix;
C 8(I) nkill; 9 betv; 10 rmscst.  11-20 are outputs:
C 11 rmsdel 12 rms2 13 nmix 14 actual broy 15 actual beta 16-24 tj
C 25 1 if wt ne 0, 10 if wa ne 0, 11 if all nonzero
C 27..29: hold parms for static parms block regular mixing
C 30..32: hold parms for static parms block Euler angle mixing
C 33 : Lindhard screening parameter
      double precision dmxp(33)
C ... for iterations
      logical lhf,lbin,a2bin
      integer maxit,iter
      double precision seref,etot(2),amom,qdiff,qtol,etol,alat
      equivalence (qdiff,dmxp(11))
C ... for relaxation
      logical xyzfrz(3),lshr,ltmp
      integer icom,natrlx,nvrelx,ltb,itrlx,nm,irlxsh,nitrlx,bitor,ng,init_obp
      double precision mdprm(7),gam(4),gam1,bstim,rhosig,pletot(6,2),
     .  plat0(3,3),dist0(9),dist(9)
      parameter (nm=3)
C ... for LDA+U
      integer nlibu,lmaxu
      double precision tolu
C ... for LMFRS
      integer lmfrs
C#ifdefC LMFRS
C      integer imjob,ierr,ib
C#endif
C     For LDA+U
      integer, allocatable :: lldau(:)
      complex(8), allocatable :: vorb(:),dmatu(:),dmato(:)

      data irlxsh /0/ dist0 /9*0d0/
C ... for debugging
C      double precision vol,repnl(2),rmunl(2)

C     parameter (T=.true., F=.false.)

      call tcn('lmfp')

      nvrelx = 0
      etot(1) = 0
      allocate(vorb(1),dmatu(1),dmato(1))
      nullify(s_pot%ppn)

C      if (cmdopt('--rdbasp',8,0,strn)) then
C        fileid = 'basp'
C        if (strn(9:12) == ':fn=') then
C          fileid = strn(13:)
C        else
C        endif
C        call strip(fileid,i,j)
C        ifi = fopna(fileid(1:j),-1,0)
C        rewind ifi
C         nspec = s_ctrl%nspec
C        if (.not. ioorbp(111,2,1,nspec,sspec,s_spec,k,ifi))
C     .    call rxs2('lmfp: failed to find BASIS: token in file "',
C     .    fileid(1:j),'"')
C        call fclr(' ',ifi)
C      endif

C ... MPI-specific
      nproc  = mpipid(0)
      procid = mpipid(1)
      master = 0
      mlog = cmdopt('--mlog',6,0,strn)

C ... popt file editor
      if (cmdopt('--popted',8,0,strn)) then
        call popted(strn(9:),s_ctrl,s_lat,s_bz)
        call rx0('exit popted')
      endif

C -------------------------- Basis optimization -------------------
    2 continue
      if (cmdopt('--optbas',8,0,strn)) then
C       call wkdbg2
C       No forces or dynamics
        call dpzero(mdprm,7)
        s_ctrl%lfrce = 0
        s_ctrl%mdprm = mdprm
        call lmfopb(strn(9:),s_lat,s_spec,etot(1)-s_ham%seref,init_obp)
C       No self-consistency
        if (init_obp == 0) s_bz%nevmx = -1
C       s_ctrl%quit = struc
C       call wkdbg2

C -------------------------- Total energy mode -------------------
      elseif (cmdopt('--etot',6,0,strn)) then
C       No forces or dynamics
        call dpzero(mdprm,7)
        s_ctrl%lfrce = 0
        s_ctrl%mdprm = mdprm
C       Suppress writing output density
        s_ctrl%lrs = s_ctrl%lrs - IAND(s_ctrl%lrs,16+8)
C       Exactly one iteration
        s_ctrl%maxit = 1
C       Suppress mixing of output density
        j = str_pack('mix',1,s_strn,'none')
      endif

C -------------------------- Unpack & initialization -------------------
      call getpr(ipr)
      nbas = s_ctrl%nbas
      nl = s_ctrl%nl
      nspec = s_ctrl%nspec
      nsp = s_ctrl%nspin
      nat = nglob('nat')
      nsp = s_ctrl%nspin
      npadl = s_ctrl%npadl
      npadr = s_ctrl%npadr
      lpgf = s_ctrl%lpgf(1)
      lcgf = s_ctrl%lcgf
      iscr = s_ctrl%lscr
      lsx = s_ctrl%lsx

      if (s_ctrl%lwsig == 8) then
      elseif (cmdopt('--eveca',7,0,strn)) then
        s_ctrl%lwsig=7
        call info0(30,1,0,
     .    ' LMFP:   eigenvectors will be written to ASCII file eveca')
      elseif (cmdopt('--evecl',6,0,strn)) then
      elseif (cmdopt('--evec',6,0,strn)) then
        ifi = fopna('eveca',-1,0); rewind ifi
        s_ctrl%lwsig=5
        call info0(30,1,0,
     .    ' LMFP:   eigenvectors will be written to file evec')
      endif

C     lncol = s_ctrl%lncol

      qbg    = s_ctrl%zbak(1)
      maxit  = s_ctrl%maxit
      lrel   = isw(mod(s_ctrl%lrel,10) /= 0)
      lhf    = IAND(s_ctrl%lcd,2) /= 0

      if (lhf) maxit = min(maxit,1)
C     nbasp  = nbas +    npadl + npadr
      nbaspp = nbas + 2*(npadl + npadr)
      stdo   = lgunit(1)
C     stdl   = lgunit(2)
      call setcc(lrel,0d0)

      irs(1) =     IAND(s_ctrl%lrs,7)
     .       +     8*isw(IAND(s_ctrl%lrs,256) /= 0)
      irs(2) =     IAND(s_ctrl%lrs,8+16)/8 + 10*isw(IAND(s_ctrl%lrs,512) /= 0)
      irs(3) = isw(IAND(s_ctrl%lrs,32) /= 0)
      irs(4) = isw(IAND(s_ctrl%lrs,64) /= 0)
      irs(5) = isw(IAND(s_ctrl%lrs,128) /= 0)

C ... --rs=3 => always read from atom file
      if (IAND(s_ctrl%lrs,3) == 3) irs(1) = 0

C     Sanity checks: most ASA "extras" are not implemented here
      call sanrg(.true.,lcgf, 0,0,'lmfp:','lcgf')
      call sanrg(.true.,lpgf, 0,0,'lmfp:','lpgf')
      call sanrg(.true.,iscr, 0,0,'lmfp:','lscr')
      call sanrg(.true.,lsx,  0,0,'lmfp:','lsx')
C     call sanrg(.true.,lncol,0,0,'lmfp:','lncol')

C ... Printout properties of species
      if (ipr >= 30) then
C       call pr_basis (sspec,0)
        call praugm(s_spec,0)
      endif

C ... Setup for no screening transformation
      plat = s_lat%plat
      if (IAND(s_ctrl%lbas,1) /= 0 .and. s_ham%lsig == 0) then

C       Shorten site positions
        if (.not. cmdopt('--shorten=no',12,0,strn)) then
C       unpack from site structure
        allocate(pos2(3*nbas))
        call sitepack(s_site,1,nbas,'pos',3,xx,pos2)
        ix(1) = 2
        ix(2) = 2
        ix(3) = 2
C       call shorpsa(nbas,plat,ix,afmt,npgrp,s_lat%istab,pos2,s_lat%pos)
        call shorps(nbas,plat,ix,pos2,s_lat%pos)
        call sitepack(s_site,1,nbaspp,'-pos',3,xx,s_lat%pos)
C       Debugging printout
C       call prmx('starting basis vectors',pos2,3,3,nbas)
C       call prmx('shortened basis vectors',s_lat%pos,3,3,nbas)
        call daxpy(3*nbas,-1d0,s_lat%pos,1,pos2,1)
        if (dlength(3*nbas,pos2,1) > 1d-6)
     .    call info0(30,1,0,' lmfp : shortened basis vectors ... ')
        deallocate(pos2)
        endif
      endif

C ... Setup for charge mixing
      call dpzero(dmxp,33)
      dmxp(2) = s_mix%b
      dmxp(9) = s_mix%bv
      dmxp(4:6) = s_mix%w
      dmxp(3) = s_mix%wc
      dmxp(6) = s_mix%nsave
      dmxp(7) = s_mix%mmix
      call parms0(0,0,0d0,0)

C ... Allocate memory for forces
      lfrce = s_ctrl%lfrce
      if (lfrce /= 0) then
        numq = 1
        if (s_bz%lmet == 4) numq = 3
        if (.not. allocated(frc)) allocate(frc(3*nbas*numq))
      endif
        if (.not. allocated(frc)) allocate(frc(1)) ! Enables lmf to pass bounds check

C ... Relaxation setup
      itrlx = 1
C     nstack = 0
C     Initial shear was already folded into plat
      gam = s_lat%gam
      gam1 = gam(4)
      gam(4) = 1
      s_lat%gam = gam
      s_lat%ldist = 0

      nitrlx = s_ctrl%nitmv
      mdprm = s_ctrl%mdprm
      ltb = s_ctrl%ltb
      lshr = nint(mdprm(1)) > 100
      if (nint(mdprm(1)) == 0) nitrlx = 0
      if (nint(mdprm(1)) > 0 .and. nint(mdprm(1)) < 4) then
        call rx('lmf not set up for MD yet')
      endif
      if (nitrlx > 0) then
        s_ctrl%mdprm = mdprm
        s_ctrl%ltb = bitor(ltb,16)
        if (.not. allocated(indrlx)) allocate(indrlx(2,3*(nbas+6)))
C       Next lines in case lattice relaxation
        if (lshr) then
          if (abs(gam(4)-1) > 1d-10) call rx('lmfp: '//
     .        'use of SHEAR= incompatible w/ lattice relaxation')
          plat0 = s_lat%plat0
        endif
        call rlxstp(s_ctrl,s_site,natrlx,nvrelx,indrlx,xyzfrz,pdim)
        icom = 0
        if (nvrelx /= 0) then
          if (natrlx > 3*(nbas+6)) call rx('reallocate indrlx')
          if (.not. allocated(grwk)) then
            allocate(grwk(pdim),hessw(nvrelx*nvrelx))
          endif
        endif
        alat = s_lat%alat
        if (procid == master) then
        ifi = fopna('bsmv',-1,0)
        allocate(pos2(3*nbas))
        j = 1
        call ivset(ix,1,3,j)
        call shorps(nbas,plat,ix,s_lat%pos,pos2)
        call iobsm0(0,bstim,0d0,0d0,nbas,alat,pos2,ifi)
        deallocate(pos2)
        endif
      endif

C ... Re-entry for shear distortion
    4 continue

C ... Potential setup
      call info0(50,0,0,' lmfp : potential setup ... ')
      call supot(0,s_ctrl,s_lat,s_pot)

C#ifdefC LMFRS
CCC ... allocate str structures if has not been done before
CC      if (.not. allocated(sstr_bas)) then
CC        allocate(sstr_bas(nbas), stat=ierr)
CC        if(ierr /= 0)
CC     .    call rx0(' lmfp: allocation of sstr_bas failed')
CC      endif
CC ... initialize pointers into disassociated state
Cc     nullify (imesh,fmesh,gmesh)
Cc ... initialize and allocate parent pointer arrays and nullify subpointers
C      nullify (sstr_mesh)
C      allocate(sstr_mesh(nbas), stat=ierr)
C      if(ierr == 0) then
C        do ib = 1, nbas
C          nullify (sstr_mesh(ib)%ind)
C          nullify (sstr_mesh(ib)%han)
C          nullify (sstr_mesh(ib)%gs)
C        enddo
C      else
C        call rx0(' lmfp: allocation of sstr_mesh failed')
C      endif
C      nullify (sstr_x)
C      allocate(sstr_x, stat=ierr)
C      if(ierr == 0) then
C        nullify (sstr_x%npr,sstr_x%nprg,sstr_x%iax,
C     .           sstr_x%s,sstr_x%sg,sstr_x%hsvx,sstr_x%hslx)
C      else
C        call rx0(' lmfp: allocation of sstr_x failed')
C      endif
CC ... allocate str structures if has not been done before
C      if (.not. allocated(sstr_bas)) then
C        allocate(sstr_bas(nbas), stat=ierr)
C        if(ierr /= 0)
C     .    call rx0(' lmfp: allocation of sstr_bas failed')
C      endif
Cc...deb
C        print *,' lmfp: size(sstr_bas) = ',size(sstr_bas)
C        print *,' lmfp: size(sstr_mesh) = ',size(sstr_mesh)
Cc...deb
C
Cc...deb
Cc lmfrs, rmaxm, etc should have been set earlier. Before this is programmed, do it here
Cc     sstr_gen%lmfrs = 0
C      sstr_gen%lmfrs = 1
C      alat = s_lat%alat
Cc     sstr_bas(1:nbas)%rmaxm = 1.5d0*alat
C      sstr_bas(1:nbas)%rmaxm = 2.0d0*alat
Cc     sstr_bas(3)%rmaxm = 1.8d0*alat
Cc...deb
Cc ... copy lmfrs from structure
C      lmfrs = sstr_gen%lmfrs
C        print *,' lmfp: lmfrs = ',lmfrs
C        print *,' lmfp: alat = ',alat
C        print *,' lmfp: rmaxm = ',sstr_bas(1:nbas)%rmaxm
C
C
CC ... arrange extended supercell for real-space mode
C      if (lmfrs /= 0) then
C        call info0(50,0,0,' lmfp : extended supercell setup ... ')
C        call subsc(s_lat,s_ctrl,sstr_gen,sstr_bas)
C      endif
C        print *,' lmfp: nbsc  = ',sstr_gen%nbsc(1:3)
C        print *,' lmfp: rmaxm = ',sstr_bas(1:nbas)%rmaxm
C#else
      lmfrs = 0
C#endif

C ... Setup of hamiltonian, augmentation
      if (ipr >= 50) then
        call info0(50,0,0,' lmfp : basis setup ... ')
      else
        call info0(30,0,0,' ')
      endif
C     call subasi(s_ctrl,s_spec,s_ham)
      call suham(s_ctrl,s_lat,s_spec,s_site,s_ham,s_pot,s_str,s_strn)

C --- Setup for iterations in a self-consistency cycle ---
C ... Unpack or allocate some permanent arrays
      ndham = s_ham%ndham

C ... Set various switches
C     Whether forces, and how to calculate non Helman-Feynman corr.
      lfrce = s_ctrl%lfrce
C     Maximum number of eigenvalues
      nevmx = s_bz%nevmx
C     Whether to evaluate output density and/or KS energy
      lrout = 1
      leks = 1
      if (s_ctrl%plbnd == 2 .or. s_ctrl%plbnd == 3) then
        lrout = 0
        leks = 0
        s_ctrl%lfrce = 0
        lfrce = 0
        nevmx = -1
      elseif (cmdopt('-leks=',6,0,strn)) then
        j = 6
        if (.not. a2bin(strn,leks,2,0,' ',j,len(strn)))
     .    call rxs('failed to parse',strn)
      endif
      if (nevmx == -1) then
        lrout = 0
        leks = 0
        s_bz%nevmx = nevmx
      endif
C     Whether to float pnu's
      lpnu = 1
C     Sanity checks
      if (lrout == 0 .and. lfrce /= 0) then
        write(stdo,333) 'when forces sought'
  333   format('lmfp (fatal): output density required ',a/
     .    '      To make output density turn off HF=t and/or NEVMX<0')
        call rx('incompatible input')
      endif
C     Sanity checks
      if (lrout == 0 .and. cmdopt('--etot',6,0,strn)) then
        write(stdo,333) 'with --etot switch.'
        call rx('incompatible input')
      endif
      if (lrout == 0 .and. maxit > 1 .and. ipr >= 20) then
        call awrit1('%N lmfp (warning): %i iterations sought but no'//
     .    ' output rho ... do 1 iteration',' ',80,stdo,maxit)
        maxit = 1
      endif

C ... LDA+U initialization
      if  (.not. allocated(lldau)) allocate(lldau(nbas))
      call iinit(lldau,nbas)
C     Check for LDA+U ... return nlibu > 0 if any U blocks.
      call suldau(nbas,s_spec,s_site,nlibu,lmaxu,lldau)
      s_ham%nlibu = nlibu
      s_ham%lmaxu = lmaxu
C ... Read LDA+U Hamiltonian
      if (nlibu > 0) then
        deallocate(vorb,dmatu,dmato)
        i = nsp*nlibu*(lmaxu*2+1)**2
        allocate(vorb(i)); call dpzero(vorb,2*i)
        allocate(dmatu(i)); call dpzero(dmatu,2*i)
        allocate(dmato(i)); call dpzero(dmato,2*i)
C       need group info to symmetrize site density matrix
        ng = s_lat%nsgrp
C       defaults
        umix = s_mix%umix
        tolu = s_mix%tolu
C       if (umix == 0) umix = 1  ! Comment out May 2017
C       if (tolu == 0d0) tolu = 1d-4
C       initialize density matrix for LDA+U
        call sudmtu(nbas,nsp,nlibu,lmaxu,s_site,s_spec,0,lldau,
     .    ng,s_lat%symgr,s_lat%istab,dmatu,vorb)
      endif
C     end LDA+U  initialization section

C ... Invoke the restart editor
      if (cmdopt('--rsedit',8,0,strn)) then
        call rsedit(strn(9:),1,s_ctrl,s_site,s_spec,s_lat,s_pot,s_ham,
     .    s_bz,nbas,nat,nspec)
        call rx0('lmfp from rsedit')
      endif

C ... Invoke the response function editor
      if (cmdopt('--chimedit',10,0,strn)) then
        call chimedit(strn(11:),1,
     .    s_ctrl,s_site,s_spec,s_lat,s_pot,s_bz,nbas,nat,nspec)
        call rx0('lmfp from chimedit')
      endif

C --- Create initial file for scaling parameters ---
      if (procid == master .and. fxst('shfac') == 0
     .   .and. IAND(s_ctrl%lham,512) /= 0) then
        allocate(shfac(3,nbas))
        call spec2class(s_spec,nbas,s_ctrl%ips,'shfac',3,xx,shfac)
        if (dlength(3*nbas,shfac,1) /= 0) then
          ifi = fopn('shfac')
          rewind ifi
          call info0(20,1,0,' lmfp : writing initial shfac file')
C         call iomagf(2,nbas,1,shfac,shfac,1,-ifi)
          call ioextf(2,'shfac file',nbas,1,shfac,shfac,3,1,-ifi)
          call fclr('shfac',ifi)
        endif
        deallocate(shfac)
      endif

C --- Read shfac file for scaling parameters, incl XC field ---
      s_pot%bxcscali = 1
      call ptr_pot(s_pot,8+1,'shfac',3,nbas,xx)
      ifi = 0
      if (procid == master .and. IAND(s_ctrl%lham,512) /= 0) then
        if (fxst('shfac') == 1) then
          ifi = fopn('shfac'); rewind ifi
C         call info0(20,0,0,' lmfp :  reading initial shfac file')
C         call iomagf(2,nbas,1,s_pot%shfac,s_pot%shfac,1,ifi)
          call ioextf(2,'shfac ',nbas,1,s_pot%shfac,s_pot%shfac,3,1,ifi)
          call fclr('shfac',ifi)
        else
         call info0(20,1,0,' LMFP (warning) : HAM_SHFAC set but no shfac file')
        endif
      endif
      call mpibc1(ifi,1,2,.false.,'lmfp','ifi')
      if (ifi /= 0) then
        call mpibc1(s_pot%shfac,nbas*3,4,.false.,'lmfp','shfac')
        if (IAND(s_ctrl%lbxc,8) == 8) then
          s_pot%bxcscali = 1 + dsum(nbas,s_pot%shfac,3)/nbas
          call info2(20,0,0,'%9fscale interstitial Bxc by: %,6;6d',s_pot%bxcscali,0)
        endif
      endif

C --- Read socscl file for spin-orbit scaling parameters ---
      call ptr_pot(s_pot,8+1,'socscl',nbas,0,xx)
      call dvset(s_pot%socscl,1,size(s_pot%socscl),1d0)
      if (bittst(s_ctrl%lncol,4)) then  ! If SOC turned on
        if (procid == master) then
        if (fxst('socscl') == 1) then
          ifi = fopn('socscl')
          rewind ifi
          call ioextf(2,'SO scaling from socscl',nbas,1,s_pot%socscl,s_pot%socscl,1,1,ifi)
          call fclr('socscl',ifi)
          write(stdo,901) ! (s_pot%socscl(i),i=1,nbas)
  901     format(' SOC coefficients scaled by: ':/16F6.2)
          call arrprt(' Site  scale','%,5i%;7,3D','Id',nbas,stdo,7,0,' | ',
     .      xx,s_pot%socscl,xx,xx,xx,xx,xx,xx)
        endif
        endif
        call mpibc1(s_pot%socscl,nbas,4,.false.,'gfasa','socscl')
      endif

C ... Quit if --quit=ham given
      if (s_ctrl%quit == 8) call rx0('quit = ham')

C ---------------- Re-entry point for a new iteration ---------------
      iter = 1
    5 continue

C --- Read restart file or overlap free atom densities ---
C     irs(1) tells what to read and whether to invoke smshft.
C     4s' bit of irs(1) -> invoke smshft after file read.
C     8s' bit of irs(1) -> rotate local density after file read
C     0+1's bits irs(1)     action
C           0              read from atom file
C           1              read from binary rst file
C           2              read from ascii rsta file
C           3              read nothing (data already input)
   10 continue
C     Harris-Foulkes -> always overlap free-atom densities
      if (irs(1) == 0) then
        call rdovfa(nbas,nspec,s_site,s_spec,s_lat,s_pot,s_ham,qbg)
        nit1 = 0
      elseif (mod(irs(1),4) >= 1 .and. mod(irs(1),4) <= 2) then
        lbin = .not. bittst(irs(1),2)
        ifi = -1
        if (procid == master) then
          if (lbin) then
            if (fxst('rst') /= 1) then
              call info0(10,1,-1,' lmfp  : '//
     .          'no rst file ... try to overlap atomic densities')
              irs(1) = 0
C             goto 10  ! defer jumping so MPI jobs do not hang
            else
              ifi = fopna('rst',-1,4)
            endif
          else
            if (fxst('rsta') /= 1) then
              call info0(10,1,-1,' lmfp  : '//
     .          'no rsta file ... try to overlap atomic densities')
              irs(1) = 0
C             goto 10  ! defer jumping so MPI jobs do not hang
            else
              ifi = fopna('rsta',-1,0)
            endif
          endif
        endif
        call mpibc1(irs(1),1,2,mlog,'lmfp','ifi')
        call mpibc1(ifi,1,2,mlog,'lmfp','ifi')
        if (ifi == -1) goto 10  ! No file existed; attempt to read overlap matrix
        k = iorsf(1,s_ctrl,s_site,s_spec,s_lat,s_pot,s_bz,
     .    fileid,nbas,nat,nspec,xx,nit1,lbin,ifi)
        if (procid == master) then
          call fclose(ifi)
        endif
        call mpibc1(k,1,2,mlog,'lmfp','k')
        if (k < 0) then
C         irs(1) = irs(1) - mod(irs(1),4)
C          call rx('MPI: rst read failed. Restart with --rs=0')
          irs(1) = 0
          goto 10
        endif

C   ... Write positions array from site structure
        call sitepack(s_site,1,nbas,'pos',3,xx,s_lat%pos)
        call bcast_strx(2**6,xx,xx,xx,xx,s_lat,xx,xx,xx,xx,0,0)
C        if (k < 0) then
C          irs(1) = 0
C          goto 10
C        endif
        if (mod(irs(1),8) >= 4) then
C         If no force switch set, use default
          k = s_ctrl%lfrce
          if (k == 0) s_ctrl%lfrce = 1
          call dfratm(s_site,s_spec,8,1,nbas,s_pot%rhat) ! link local rho to s_pot%rho
          call smshft(s_site,s_spec,s_lat,s_ctrl,s_ham,s_pot)
C         Restore force switch
          if (k == 0) s_ctrl%lfrce = k
        endif
        if (mod(irs(1),16) >= 8) then
C          xv(1:9) = s_lat%dist(1:9,1)
C          call pvsms2(s_site,s_spec,xv,nbas,nsp)
          irs(1) = irs(1)-8
C          i = i-256
          i = s_ctrl%lrs - IAND(s_ctrl%lrs,256)
          s_ctrl%lrs = i
        endif
      endif

C ... Write atmx file
      if (cmdopt('--wratrho',9,0,strn)) then
        i = fopna('atmx',-1,0)
        call ioatmx(0,1,nbas,nspec,s_site,s_spec,s_ham,-i)
      endif

      if (cmdopt('--rdatrho',9,0,strn) .and. iter == 1) then
        j = 16+8+4+2
        i = fopna('atmx',-1,0)
        if (strn(10:10) /= ' ') then
          j = 0
          if (index(strn,strn(10:10)//'rhoc') > 0) j = j+2
          if (index(strn,strn(10:10)//'rhov') > 0) j = j+8
          if (index(strn,strn(10:10)//'v0') > 0) j = j+16
        endif
        call ioatmx(j,1,nbas,nspec,s_site,s_spec,s_ham,i)
      endif

C ... Write positions after file read, and repack
      if (ipr >= 50) then
      write(stdo,357) 'Basis, after reading restart file'
  357 format(/1x,a/' site spec',8x,'pos (Cartesian coordinates)',9x,
     .  'pos (multiples of plat)')
      call dinv33(plat,1,qlat,xv)
      do  i = 1, nbas
        j = s_site(i)%spec
        xv(1:3) = s_site(i)%pos
C       call spacks(0,'spec name',sspec,alabl,j,j)
        alabl = s_spec(j)%name
        call dgemm('T','N',3,1,3,1d0,qlat,3,xv,3,0d0,xv(4),3)
        write(stdo,345) i, alabl, (xv(j),j=1,3), (xv(3+j),j=1,3)
  345   format(i4,2x,a8,f10.6,2f11.6,1x,3f11.6)
      enddo
      endif

C --- Optionally re-shorten basis vectors ---
      if (cmdopt('--shorps',8,0,strn)) then
        allocate(pos2(3*nbas))
        call sitepack(s_site,1,nbas,'pos',3,xx,pos2)
        ix(1) = 2
        ix(2) = 2
        ix(3) = 2
        call shorps(-nbas,plat,ix,pos2,s_lat%pos)
        call info0(30,1,-1,
     .    ' lmfp  : write shortened vectors to file shorps ...')
        call iopos(.true.,1,'shorps',nbas,s_lat%pos,s_site)
        call shorps(nbas,plat,ix,pos2,s_lat%pos)
        call sitepack(s_site,1,nbas,'-pos',3,xx,s_lat%pos)
C       Debugging printout
C       call prmx('starting basis vectors',pos2,3,3,nbas)
C       call prmx('shortened basis vectors',s_lat%pos,3,3,nbas)
c       print *,' lmfp: nbas, nbaspp = ',nbas,nbaspp
c       print *,' lmfp: starting lattice vectors:'
c       call xxr(pos2,3*nbas)
c       print *,' lmfp: shortened lattice vectors:'
C       call sitepack(s_site,1,nbas,'pos',3,xx,s_lat%pos)
c       call xxr(s_lat%pos,3*nbaspp)
        deallocate(pos2)
      endif

C --- Write to disk lattice shifts caused by iorsf or shorps ---
      if (cmdopt('--phase',7,0,strn)) then
        allocate(pos2(3*nbas),wk(3*nbas))
        call sitepack(s_site,1,nbas,'pos',3,xx,pos2)
        call sitepack(s_site,1,nbas,'pos0',3,xx,wk)
        call daxpy(3*nbas,-1d0,wk,1,pos2,1)
C        ix(1) = 0
C        ix(2) = 0
C        ix(3) = 0
        call dgemm('T','N',3,nbas,3,1d0,s_lat%qlat,3,pos2,3,0d0,wk,3)
        call iopos(.true.,1,'phase',nbas,wk,s_site)
        deallocate(pos2,wk)
        call rx0('wrote phase file')
      endif


C#ifdefC LMFRS
Cc     if (lmfrs /= 0 .and. irs(1) /= 3) then
C      if (irs(1) /= 3) then
C        imjob = 1                                               ! imjob = 1 make functions on the mesh
Cc       imjob = 0                                               ! imjob = 0 don't make functions on the mesh
C          print *,' lmfp: nl = ',nl
C          print *,' lmfp: rmaxm = ',sstr_bas(1:nbas)%rmaxm
CC ...   update the pair table
C         call mkiaxg(1,ssite,sstr,
C     .  s_ctrl,s_lat,s_site,s_spec,s_str,
C     .  sstr_x)
CC ...   make structure constants
C        call info0(30,0,0,' lmfp : make structure constants ...')
C        call mkstr(10,ssite,sstr,sstr_x,nbas,nl)
CC ...   make basis functions on the mesh
C        if (imjob /= 0) call info0(30,0,0,
C     .  ' lmfp : set real space clusters on the mesh ...')
C        call sumshv(imjob,ssite,
C     .  s_ctrl,s_lat,s_site,s_spec,sstr,sstr_gen,sstr_bas,
C     .  sstr_x,sstr_mesh,nbas,nl)
Cc...deb
CC stop for now
C        return
Cc...deb
C      endif
C#endif

C#ifdef DMFT
      call sudmft(s_bz,s_ctrl,s_ham,s_lat,s_pot,s_spec,s_site,s_optic,s_gw,s_dmft,dmxp,iter,s_strn)
      call rx0('done sudmft')
C#endif

C     Hang on to previous site density matrix for this iteration
      if (nlibu > 0) then
        i = nsp*nlibu*(lmaxu*2+1)**2
        call zcopy(i,dmatu,1,dmato,1)
        call dpzero(dmatu,2*i)
      endif

C --- Make and diagonalize hamiltonian, make new charge density ---
      if (maxit == 0) call info0(20,1,0,
     .  ' lmfp  : zero iterations sought ... no band pass')

      call bndfp(s_strn,s_site,s_spec,s_lat,s_ctrl,s_ham,s_pot,s_bz,
     .  s_optic,s_gw,nbas,nsp,nlibu,lmaxu,lldau,leks,lrout,
     .  lfrce,lpnu,dmxp,iter,maxit,frc,dmatu,vorb)

C ... Quit if --quit=rho given
      if (s_ctrl%quit == 32) call rx0('quit = rho')

C ... Cases that require no further processing
      if (s_ctrl%plbnd == 2 .or. s_ctrl%plbnd == 3) then
        return
      endif

C ... Basis optimization: just extract etot(1) and return to opt.
      if (cmdopt('--optbas',8,0,strn)) then
        etot(1) = s_ham%ehf
        etot(2) = s_ham%ehk
        if (init_obp == 0) goto 2
      endif

C ... check convergence of dmatu and update it and vorb if necessary
      if (nlibu > 0 .and. maxit > 0 .and. lrout > 0) then
        call chkdmu(nbas,nsp,nlibu,lmaxu,s_site,s_spec,s_ham,0,
     .    dmatu,dmato,vorb,tolu,umix,lldau,ng,s_lat%symgr,s_lat%istab)
      endif

C ... Write smoothed charge density for contour plotting
      if (cmdopt('--wden',6,0,strn)) then
        if (procid == master) then
        if (maxit == 0) call dfratm(s_site,s_spec,8,1,nbas,s_pot%rhat) ! In case not already done
        call ioden(strn(7:),s_lat,s_site,s_spec,s_pot%rhat,s_pot%smrho)
        endif
        call rx0('density written to disk')
      endif

C --- Write restart file (skip if --quit=band) ---
      if (procid == master .and. s_ctrl%quit /= 4) then
      ltmp = irs(2) < 10 .or. iter > 1 .and. iter == dmxp(26)
      if (ltmp .and. irs(2) > 10) then
        call info2(10,0,0,' Writing rst file ... RMS DQ is now %1,3;3e',dmxp(27),0)
      elseif (.not. ltmp .and. irs(2) > 10) then
        call info2(10,0,0,' Suppress writing rst file ... RMS DQ=%1,3;3e  min=%1,3;3e',dmxp(11),dmxp(27))
      endif
C     Suppress saving rst file in the middle of a shear (irlxsh > 0)
      ltmp = ltmp .and. irlxsh == 0
      if (ltmp .and. irs(2) > 0 .and. (lrout > 0 .or. maxit == 0)) then
        lbin = mod(irs(2),10) /= 2
        if (lbin) fileid = 'rst'
        if (.not. lbin) fileid = 'rsta'
        if (mod(irs(2),10) == 3) then
          call word(fileid,1,i,j)
          j = j+1
          fileid(j:j) = '.'
          call bin2a(' ',0,0,iter,2,0,len(fileid),fileid,j)
          if (lbin) ifi = fopng(fileid,-1,8+4)
          if (.not. lbin) ifi = fopng(fileid,-1,8)
          call info0(10,1,-1,' lmfp:  writing to restart file '//fileid)
        else
          if (lbin) ifi = fopna(fileid,-1,4)
          if (.not. lbin) ifi = fopna(fileid,-1,0)
        endif
        i = str_pack('jobid',-1,s_strn,strn)
        fileid = 'lmfp:  ' // strn
        k = iorsf(1,s_ctrl,s_site,s_spec,s_lat,s_pot,s_bz,
     .    fileid,nbas,nat,nspec,xx,iter,lbin,-ifi)
        call fclose(ifi)
      endif
      endif ! write restart file

      if (cmdopt('--window=',9,0,strn))
     .  call rx0('lmf : early exit (--window option)')

C ... Continue basis optimization; now read from rst file
      if (cmdopt('--optbas',8,0,strn)) then
        i = s_ctrl%lrs
        if (mod(i,2) == 0) i = i+1
        s_ctrl%lrs = i
        goto 2
      endif

C --- Add to save file; decide on next iteration ---
      if (maxit > 0) then
      etot(1) = s_ham%ehf
      etot(2) = s_ham%ehk
C ... Subtract reference energy
      seref = s_ham%seref
      etot(1) = etot(1) - seref
      if (etot(2) /= 0) etot(2) = etot(2) - seref
      amom = s_ham%eterms(15)
C     The desired tolerances in q,e
      qtol = s_ctrl%tol(1)
      etol = s_ctrl%tol(3)
      if (procid == master) then
        rhosig = s_ham%eterms(19)
        i = 0
        if (rhosig /= NULLI .and. rhosig /= 0) i = 10
C     .    lhf.or.irs(1) == 0.and.iter == 1,leks+i,etol,qtol,qdiff,
C     .    'cxhi',amom,etot,lsc)
        call nwit(s_ctrl%nvario,iter,maxit,lhf.or.irs(1)==0.and.iter==1,
     .    leks+i,etol,qtol,qdiff,s_ctrl%tol(4),s_ctrl%tol(5),'cxhi',amom,etot,lsc)
      endif
      call mpibc1(lsc,1,2,mlog,'lmfp','lsc')
      if (lsc == 2 .and. .not. lhf .and. maxit > 1) lsc = 3
      if (lsc == 1 .and. lrout > 0  .or. lsc == 3) then
        call querym('maxit',2,maxit)
        if (iter >= maxit) lsc = 1
        if (iter < maxit) lsc = 3
      endif

      if (s_ctrl%quit == 4) call rx0('lmf : exit (--quit=band)')

C ... Write forces to file
      if (cmdopt('--wforce=',9,0,strn) .and. procid == master)
     .    call iopos(.true.,0,strn(10:),nbas,frc,s_site)

C ... Continue iterations toward self-consistency
      iter = iter+1
      if (lsc > 2) then
        irs(1) = 3
        goto 5
      endif

C ... Reset quantities for iterations towards self-consistency
      if (nvrelx > 0 .and. nitrlx > 0) then
        iter = 1
        dmxp(11) = 0
      endif
      endif

C ... Write site file
      if (cmdopt('--wsite',7,0,strn)) then
        i = 1000*(2+4+8+16+32) + 1 ! for --wsite
        k = 7
        if (strn(1:8) == '--wsitex') then
          i = 1000*(2+4+8+16+32) + 1 + 10
          k = 8
        endif
        allocate(slabl(nspec))
        slabl(1:nspec) = s_spec(1:nspec)%name
        if (iosits(i,3d0,0,strn(k+2:),ifi,slabl,s_lat%alat,plat,nbas,
     .    nspec,s_spec,s_site) < 0) call rx('failed to write site')
        deallocate(slabl)
      endif

C ... Write positions to file
      if (cmdopt('--wpos=',7,0,strn) .or. cmdopt('--wpos:mode1:',13,0,strn)) then
        gam = s_lat%gam
        allocate(wk(3*nbas))
        call rdistn(s_lat%pos,wk,nbas,gam(1),gam(2),gam(3),1/gam1)
        call iopos(.true.,0,strn(8:),nbas,wk,s_site)
        deallocate(wk)
      endif

C --- Molecular statics ---
      if (maxit > 0) then
      if (nitrlx > 0 .and. lsc <= 2) then
        call cppos(1,nbas,s_site)
        call sitepack(s_site,1,nbas,'pos',3,xx,s_lat%pos)
        mdprm = s_ctrl%mdprm
        if (lshr) then
          call grdepl(nvrelx,indrlx,0.01d0,etot,irlxsh,pletot,dist)
          if (irlxsh /= 0) then
            call grdep2(1,nvrelx,indrlx,dist0,dist)
            goto 98
          else
            call rx('lmfp: update shear relaxation step (pos)')
            call relax(prgnam,s_ctrl,s_site,s_spec,itrlx,indrlx,
     .        natrlx,nvrelx,pletot,grwk,hessw,0,0d0,dist0,icom)
            call dpzero(dist,6)
            call grdep2(1,nvrelx,indrlx,dist0,dist)
            dist(7) = 1
          endif
        else
          call relax(prgnam,s_ctrl,s_site,s_spec,itrlx,indrlx,natrlx,
     .      nvrelx,frc,grwk,hessw,0,0d0,s_lat%pos,icom)
          call sitepack(s_site,1,nbas,'-pos',3,xx,s_lat%pos)
        endif

C       Restore lattice symmetry to machine precision
        if (cmdopt('--fixpos',8,0,strn)) then
          ng = s_lat%nsgrp
C         call shoist(0,s_lat%istab,nbas,s_lat%ag,s_lat%symgr,ng,0)
          j = 8+1
          if (strn(9:13) == ':tol=') then
            j = 13
          endif
          if (strn(9:9) /= ':' .or.
     .      .not. a2bin(strn,fptol,4,0,' ',j,len(strn))) fptol = 1d-5
          call fixpos(s_lat%pos,nbas,fptol,ng,plat,s_lat%symgr,s_lat%ag,
     .      s_lat%istab)
        endif

C       Write relaxed positions to file
        if (cmdopt('--wpos=',7,0,strn)) then
          gam = s_lat%gam
          allocate(wk(3*nbas))
          call rdistn(s_lat%pos,wk,nbas,gam(1),gam(2),gam(3),1/gam1)
          call iopos(.true.,0,strn(8:),nbas,wk,s_site)
          deallocate(wk)
        endif

C       Write updated positions to bsmv file
        if (procid == master .and. .not. lshr) then
          ifi = fopna('bsmv',-1,0)
          call poseof(ifi)
          bstim = bstim+1
          allocate(pos2(3*nbas))
          j = 1
          call ivset(ix,1,3,j)
          call shorps(nbas,plat,ix,s_lat%pos,pos2)
          call iobsmv(0,bstim,0d0,0d0,nbas,alat,pos2,-ifi)
          call fclose(ifi)
          deallocate(pos2)
        endif
C       repack updated positions in site structure
        call sitepack(s_site,1,nbas,'-pos',3,xx,s_lat%pos)

C        if (icom == 1) then
C          print *, '!!'
C        endif

C   ... Exit when relaxation converged or maximum number of iterations
        if (icom == 1) then
          if (procid == master) then
            k = s_ctrl%nvario
            flg = 'C67'
            call nwitsv(1+2,k,flg,nsp,amom,etot)
          endif
          call tcx('lmfp')
          call fexit(0,111,
     .    ' LMFP: relaxation converged after %i iteration(s)',itrlx)
        else
          call query('proceed with next relaxation step',-1,0)
        endif

C   ... Restore minimum gradient positions if this is last step
C        if (itrlx == nitrlx .and. icom == -1) then
        if (itrlx == nitrlx) then

          if (.not. lshr) then
            call info0(20,1,0,' lmfp: restore positions for minimum g')

C           call prmx('initial positions',s_lat%pos,3,3,nbas)
            call prelx1(1,nm,lshr,natrlx,nvrelx,indrlx,grwk,s_lat%pos)
C           call prmx('minimum-g positions',s_lat%pos,3,3,nbas)
          else
            call prelx1(1,nm,lshr,natrlx,nvrelx,indrlx,grwk,dist0)
            call dpzero(dist,6)
            call grdep2(1,nvrelx,indrlx,dist0,dist)
            call info2(20,0,0,
     .        ' lmfp : strain of minimum gradient:'//
     .        '%N   PDEF=%6;8,4D'//
     .        '%N STRAIN=%6;8,4D',
     .        dist0,dist)
          endif

C         Repack updated positions in site structure
          call sitepack(s_site,1,nbas,'-pos',3,xx,s_lat%pos)

        endif

C   ... New density after atom shifts
C       If explicitly told to read from atom files after atom movmment
        if (IAND(s_ctrl%lrs,3) == 3) then
          irs(1) = 0

C       Else, use self-consistent
        else if (.not. lshr) then
          irs(1) = 3
          call smshft(s_site,s_spec,s_lat,s_ctrl,s_ham,s_pot)
        endif

C   ... Write restart file (to include new positions)
        if (procid == master .and. .not. lshr) then
          ifi = fopna('rst',-1,4)
          i = str_pack('jobid',-1,s_strn,strn)
          fileid = 'lmfp:  ' // strn
          k = iorsf(1,s_ctrl,s_site,s_spec,s_lat,s_pot,s_bz,
     .      fileid,nbas,nat,nspec,xx,iter,.true.,-ifi)
          call fclose(ifi)
        endif

C   ... Remove mixing file
        if (procid == master) then
          call info0(20,0,0,' Delete mixing and band weights files ...')
        ifi = fopna('mixm',-1,4)
        call dfclos(ifi)
        ifi = fopna('wkp',-1,4)
        call dfclos(ifi)
        endif
C       reset mixing block
        call parms0(0,0,0d0,0)

C   ... Exit when maximum number of iterations encountered
        if (itrlx == nitrlx) then
          if (procid == master) then
            call tcx('lmfp')
            call fexit(1,111,
     .    ' LMFP: relaxation incomplete after %i iteration(s)',nitrlx)
          else
            call tcx('lmfp')
            call fexit(1,111,' ',0)
          endif
        endif
        itrlx = itrlx+1

        if (lshr) then
          goto 98
        else
          goto 5
        endif
      endif
      endif

      call tcx('lmfp')
C     if (maxit == 0) call rx0(' zero iterations sought ... quitting')
      return

C --- Setup to start calculation at new shear ---
   98 continue

      if (procid == master) then
        call info0(20,0,0,' Delete mixing and band weights files ...')
        ifi = fopna('mixm',-1,4)
        call dfclos(ifi)
        ifi = fopna('wkp',-1,4)
        call dfclos(ifi)
      endif

C     Restore plat, pos to their undistorted state:
C     undo original transformation = P P_0^-1
      allocate(pos2(3*nbas))
      call dinv33(plat,0,xv,xv(10))
      call dgemm('N','N',3,nbas,3,1d0,xv,3,s_lat%pos,3,0d0,pos2,
     .  3)
C     Simultaneously pack in lat->pos and site->pos
      call dgemm('N','N',3,nbas,3,1d0,plat0,3,pos2,3,0d0,
     .  s_lat%pos,3)
      deallocate(pos2)
C     call prmx('positions for plat0',s_lat%pos,3,3,nbas)
      call sitepack(s_site,1,nbas,'-pos',3,xx,s_lat%pos)
      call cppos(1,nbas,s_site)
C     New shear
      s_lat%plat = plat0
      s_lat%ldist = 3
      s_lat%dist(1:3,1:3) = reshape(dist(1:9),[3,3])
C     A little memory leakage rel to 1st pass, but not so serious
      call lattic(s_lat,s_ctrl,s_site)
      plat = s_lat%plat
C     Remake qp
      ltmp = IAND(s_ctrl%lmet,1) /= 0 .or. IAND(s_ctrl%ldos,4+2+1) /= 0
      call mkqp(s_ctrl,s_bz,s_lat,ltmp,.false.,1,-2)

C ... Write restart file (to include new positions)
      if (procid == master .and. irlxsh == 0) then
        ifi = fopna('rst',-1,4)
        i = str_pack('jobid',-1,s_strn,strn)
        fileid = 'lmfp:  ' // strn
        k = iorsf(1,s_ctrl,s_site,s_spec,s_lat,s_pot,s_bz,
     .    fileid,nbas,nat,nspec,xx,iter,.true.,-ifi)
        call fclose(ifi)
      endif

C     Decide on what density to use
      if (IAND(s_ctrl%lrs,3) == 3) then
        irs(1) = 0
C     Else, use file density
      else
        irs(1) = IAND(s_ctrl%lrs,7)
      endif
      goto 4

c...deb
c ... use contains statement
c     contains
c       include 'sumesh.f'
c...deb
      end

      subroutine cppos(ib1,ib2,s_site)
C- Copy site positions to p0 for a range of sites
      use structures
      implicit none
      integer ib1,ib2
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      double precision xx
      real(8),pointer :: p_pos(:,:)

C      do  ib = ib1, ib2
C      enddo
      allocate(p_pos(3,ib2))
      call sitepack(s_site,ib1,ib2,'pos',3,xx,p_pos)
      call sitepack(s_site,ib1,ib2,'-pos0',3,xx,p_pos)
      deallocate(p_pos)
      end

      subroutine grdepl(nvrelx,indrlx,alpha,etot,irlxsh,pletot,dist)
C-
      implicit none
      integer irlxsh,nvrelx,indrlx(nvrelx)
      double precision pletot(6,2),etot,dist(9),alpha
      double precision grad,vec1(6)
      integer iv,ipm,ipv

      call info5(30,1,0,' GRDEPL: point %i of %i for grad shear: '//
     .  'etot=%d',irlxsh,2*nvrelx,etot,0,0)

C   3 continue

C     Get index for current shear and store energy for that shear
      if (irlxsh > 0) then
        iv = (irlxsh-1)/2 + 1
        ipv = iv
C       ipv = indrlx(iv)
        ipm = mod((irlxsh-1),2) + 1
        pletot(ipv,ipm) = etot
      endif

C     If this is last point, form gradient and exit
      if (irlxsh == 2*nvrelx) then
        do  iv = 1, nvrelx
C         ipv = indrlx(iv)
          ipv = iv
          grad = (pletot(ipv,1) - pletot(ipv,2))/(2*alpha)
          pletot(ipv,1) = grad
        enddo
        irlxsh = 0
        return
      endif

C     Get shear index for next shear and whether + or -
      irlxsh = irlxsh+1
      iv = (irlxsh-1)/2 + 1
      ipv = indrlx(iv)
      ipm = mod((irlxsh-1),2) + 1
      if (ipv < 1 .or. ipv > 6)
     .  call rx('grdepl: something wrong with indrlx')
C     Make new shear
      call dvset(vec1,1,6,alpha)
      if (ipm == 2) call dvset(vec1,1,6,-alpha)
      call dpzero(dist,6)
      call grdep2(iv,iv,indrlx,vec1,dist)
      dist(7) = 1

C     goto 3
      end
