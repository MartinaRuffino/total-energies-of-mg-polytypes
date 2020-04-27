      subroutine sugws(prgnam,s_bz,s_ctrl,s_ham,s_lat,s_mix,
     .  s_move,s_pot,s_str,s_spec,s_site,s_optic,s_gw,s_strn)
C- Entry point for self-energy manipulations
C ----------------------------------------------------------------------
Cio Structures
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  nkp nkabc lshft ef egap lmet nevmx def w dos idtet
Ci                 ipq pdos qp star wtkp wtkb swtk ntet n efmax fsmom
Ci                 ndos dosw numq range
Co     Stored:     nkabc lshft star nkp ntet nevmx ef def w dos idtet
Co                 ipq pdos qp wtkp wtkb swtk numq ndos dosw
Co     Allocated:  qp wtkp idtet star ipq wtkb swtk
Cio    Elts passed: ipq wtkp qp lopt lio idtet star numq swtk wtkb egap
Cio                n w def
Cio    Passed to:  sfuned mkqp evalqp lmfp popted rsedit iorsf bcast_
Cio                strx iinit mpibc1 chimedit bndfp rdsigm subzi addrbl
Cio                optinq optin2 optint
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  lpgf lmet nbas nbasp nspin nspec nl lham lcgf lgen3
Ci                 lncol ldlm plbnd lrs lscr lsx zbak maxit lrel lfrce
Ci                 nitmv mdprm ltb quit tol nvario defm loptc pfloat
Ci                 ldos lfp
Co     Stored:     lmet plbnd lfrce mdprm lrs maxit ltb
Co     Allocated:  *
Cio    Elts passed: lmet ldos lsx lscr lcd lfp lgen3 lrs lbas nvario
Cio                lncol
Cio    Passed to:  sfuned mkqp supot suham evalqp lmfp popted rlxstp
Cio                rsedit iorsf bcast_strx iinit mpibc1 chimedit smshft
Cio                bndfp suham2 optinq optin2 optint vcdmel rxes dfrce
Cio                relax lattic
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  ndham pwmode pwemin pwemax npwpad lncol neula qss
Ci                 nbf evals ehf ehk seref eterms elind ldham lsig
Ci                 oveps ovncut rsrnge nqsig eseavr sigp pmin pmax
Ci                 pnudef rsstol nprs qsig ndhrs
Co     Stored:     ndham ndofH ldham lmxax npwmin npwpad hord nlibu
Co                 lmaxu lsig eterms sigp ehf ehk nqsig ndhrs eseavr
Co     Allocated:  offH iprmb bdots nprs hrs qsig iaxs
Cio    Elts passed: offH iprmb eula magf bdots evals seref qsig lncol
Cio                hrs iaxs nprs
Cio    Passed to:  sfuned suham evalqp lmfp smshft bndfp mkpot rdsigm
Cio                hft2rs suham2 hambls hambl sopert3 blsig makusq
Cio                addrbl mkrout mkehkf mkekin chkdmu
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat qlat npgrp nsgrp vol awald nkd nkq nabc
Ci                 gmax ng kv kv2 igv igv2 tolft gam plat0 tol dist ag
Ci                 bgv cg cy dlv gv gvq indxcg ips0 istab jcg pos qlv
Ci                 symgr afmt napw as rpad nkdmx nkqmx platl platr ldist
Co     Stored:     nabc gmax ng igv igv2 gam ldist plat dist alat ag
Co                 bgv cg cy dlv gv gvq indxcg ips0 istab jcg kv kv2
Co                 pos qlv symgr napw vol plat0 qlat platl platr awald
Co                 nkd nkq
Co     Allocated:  ips0 bgv gv kv igv kv2 igv2 dlv qlv
Cio    Elts passed: symgr pos dlv qlv nsgrp ag gv ips0 bgv qlat plat
Cio                igv igv2 istab kv cg indxcg jcg cy vol
Cio    Passed to:  sfuned mkqp supot sugvec0 sugvec suham evalqp lmfp
Cio                popted lmfopb rsedit rdovfa ovlpfa ovlocr hxpbl
Cio                ghibl hklbl gklbl iorsf bcast_strx iinit mpibc1
Cio                prsed2 prsed4 chimedit smshft pvsms1 rhgcmp symsmr
Cio                bndfp mkpot smves vesft vesgcm mshvmt symvvl ugcomp
Cio                ggugbl gfigbl fklbl hhugbl hhigbl phhigb hsmbl
Cio                hgugbl smvxc2 vxcnlm smvxcm smcorm locpot rdsigm
Cio                suham2 hambls hambl augmbl bstrux hxpgbl ghigbl
Cio                hklgbl smhsbl hhibl phhibl hsibq makusq pusq1 addrbl
Cio                fsmbl rsibl rlocbl vcdmel rxes mkrout symrat symrho
Cio                dfrce pvdf4 pvdf2 pvdf1 mkekin totfrc mixrho ioden
Cio                lattic
Cio  s_mix :contains parameters for charge mixing; see structures.h
Ci     Elts read:  b bv w wc nsave mmix umix tolu
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  sfuned evalqp lmfp
C
Cio  s_move :struct for dynamics parameters; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  sfuned evalqp lmfp
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  vmtz rho nlml nlma ves aamom bxc cp ddpf dddpf ddpfr
Ci                 dlmwt dpf dpfr gibbs gma gmar grrme mad mxy palp
Ci                 papg pf pfnc pfr pmpol pnu pp ppn pprel pti qc qcorr
Ci                 qnu qpp qt rhat rnew rhos rhrmx sop thetcl vdif
Ci                 vintr vrmax vshft smpot smrho smrout
Co     Stored:     nlma nlml ves aamom bxc cp ddpf dddpf ddpfr dlmwt
Co                 dpf dpfr gibbs gma gmar grrme mad mxy palp papg pf
Co                 pfnc pfr pmpol pnu pp ppn pprel pti qc qcorr qnu qpp
Co                 qt rhat rnew rhos rhrmx sop thetcl vdif vintr vrmax
Co                 vshft smpot smrho smrout
Co     Allocated:  mad smrho smpot rhat rnew pti smrout
Cio    Elts passed:mad pp pti rhat smrho smpot rnew smrout
Cio    Passed to:  sfuned supot suham evalqp lmfp rsedit rdovfa iorsf
Cio                bcast_strx iinit mpibc1 chimedit smshft bndfp mkpot
Cio                suham2 mixrho
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  *
Co     Stored:     nkaps nitab nttab lmaxw nds
Co     Allocated:  iax npr alp s adot sdot
Cio    Elts passed:kaps alph
Cio    Passed to:  sfuned suham rdstrx pp2alp evalqp lmfp
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa lmxl p pz idxdn rmt lmxb ncomp name ngcut orbp
Ci                 z norp ntorb rsma kmxt rg rsmv kmxv lfoca rfoca idu
Ci                 uh jh nr a qc rhoc coreh coreq ctail etail stc idmod
Ci                 nxi exi chfa rsmfa rs3 eh3 vmtz mxcst
Co     Stored:     idxdn ngcut orbp norp ntorb p pz lmxb lmxa a nr qc
Co                 nxi exi chfa rsmfa ctail etail stc rhoc name rmt z
Co                 lmxl kmxt lfoca rfoca coreh pb1 pb2
Co     Allocated:  rhoc
Cio    Elts passed:pz name rhoc
Cio    Passed to:  sfuned suham atfold makidx nscpa showbs sugcut
Cio                uspecb evalqp lmfp lmfopb lmfop2 ioorbp praugm
Cio                suldau sudmtu praldm rotycs symdmu rsedit dfratm
Cio                chimedit rdovfa bcast_strx gtpcor ovlocr corprm
Cio                adbkql iorsf pvsms2 smshft pvsms1 rhgcmp rhogkl bndfp
Cio                dfaugm mkpot rhomom smves vesgcm mshvmt symvvl
Cio                ugcomp smvxcm smcorm smvxc4 elocp locpot dfqkkl
Cio                suham2 suclst surho sumlst hambls hambl augmbl
Cio                bstrux smhsbl hsibq tbhsi hsubblock hsibq2 hsibq4
Cio                makusq pusq1 mkorbm mullmf mkpdos mkdmtu addrbl
Cio                fsmbl fsmbpw rsibl rsibl1 rlocbl iorbtm mshn3p mchan
Cio                vcdmel rxes mkrout symrat symrho prrhat pnunew dfrce
Cio                pvdf4 pvdf2 pvdf1 mkekin mixrho ftlxp pvmix5 pvmix3
Cio                pvmix7 chkdmu ioden relax
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec clabel pos relax rho1 rho2 rhoc rho1x rho2x
Ci                 rhocx class force vel pnu pz v0 v1 bxc cpawt omg
Ci                 omgn domg gc gcu gcorr sfvrtx j0 pdos qhhl qhkl qkkl
Ci                 eqhhl eqhkl eqkkl sighh sighk sigkk tauhh tauhk
Ci                 taukk pihh pihk pikk sighhx sighkx sigkkx tauhhx
Ci                 tauhkx taukkx pihhx pihkx pikkx thet pos0
Co     Stored:     pnu pz norb pos vel pos0 force spec clabel bxc cpawt
Co                 omg omgn domg gc gcu gcorr sfvrtx j0 pdos rho1 rho2
Co                 rhoc rho1x rho2x rhocx qhhl qhkl qkkl eqhhl eqhkl
Co                 eqkkl sighh sighk sigkk tauhh tauhk taukk pihh pihk
Co                 pikk sighhx sighkx sigkkx tauhhx tauhkx taukkx pihhx
Co                 pihkx pikkx thet v0 v1
Co     Allocated:  rho1 rho2 rhoc rho1x rho2x rhocx v0 v1 sigkk taukk
Co                 sigkkx taukkx pikk pikkx sighk tauhk sighkx tauhkx
Co                 pihk pihkx sighh tauhh sighhx tauhhx pihh pihhx qkkl
Co                 eqkkl qhkl eqhkl qhhl eqhhl
Cio    Elts passed: rho1 rho2 rhoc rho1x rho2x rhocx v0 qkkl qhkl qhhl
Cio                eqkkl eqhkl eqhhl
Cio    Passed to:  sfuned suham setnorb showbs pvioeu evalqp lmfp
Cio                rlxstp suldau sudmtu praldm rotycs symdmu rsedit
Cio                dfratm chimedit rdovfa ovlpfa ovlocr adbkql iorsf
Cio                pvsms2 bcast_strx smshft pvsms1 rhgcmp rhogkl iopos
Cio                bndfp dfaugm mkpot rhomom smves vesgcm mshvmt symvvl
Cio                ugcomp smvxcm smcorm smvxc4 elocp locpot dfqkkl
Cio                suham2 suclst surho sumlst hambls hambl augmbl
Cio                bstrux smhsbl hsibq hsubblock hsibq2 hsibq4 makusq
Cio                pusq1 mkorbm mullmf mkpdos mkdmtu addrbl fsmbl
Cio                fsmbpw rsibl rsibl1 rlocbl mshn3p mchan vcdmel rxes
Cio                mkrout symrat symrho prrhat pnunew dfrce pvdf4 pvdf2
Cio                pvdf1 mkekin totfrc mixrho ftlxp pvmix5 pvmix3
Cio                pvmix7 chkdmu ioden cppos relax lattic
Cio  s_optic:struct for parameters for optical ME; see structures.h
Ci     Elts read:  ltet ocrng unrng dw window lpart esciss iq
Co     Stored:     ocrng unrng
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  sfuned evalqp lmfp bndfp optinq optin2 optint
Cio  s_gw  :struct for gw-related parameters
Ci     Elts read:  nkabc lshft
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  sfuned evalqp lmfp bndfp
Ci Inputs
Ci   prgnam:name of main program
Co Outputs
Cs Command-line switches
Cs   --sfuned :Invoke the self-energy editor
Cs   --chksig :check causality of sigma
Cl Local variables
Cr Remarks
Cr   At present this code only serves as a front end to sfuned
Cu Updates
Cu   12 Jun 17 Corrections to optics and joint DOS
Cu   05 Jun 17 Optics and joint DOS
Cu   21 Apr 17 Extend this routine for self-energies not made by QSGW, eg DMFT
Cu   30 Nov 16 New 'merges12' to combine two spins in seq file
Cu   13 Nov 16 Write interpolated spectral functions, self-energies to disk
Cu   22 Feb 13 Incorporated Julien Vidal's additions for modeling PES
Cu   18 Jul 12 First created
C ----------------------------------------------------------------------
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
      type(str_dmft)::  s_dmft  ! Just a placeholder; not used
      type(str_strn) :: s_strn(*)
C ... Local variables
      integer procid,master,mpipid,nproc
      logical cmdopt
      character strn*120

      call tcn('sugws')

C ... MPI-specific
      nproc  = mpipid(0)
      procid = mpipid(1)
      master = 0

C ... Electron-phonon self-energy editor
      if (cmdopt('--ephed',7,0,strn)) then
        call ephed(strn(7+1:),s_bz,s_ctrl,s_ham,s_lat,s_mix,
     .  s_move,s_pot,s_str,s_spec,s_site,s_optic,s_gw,s_dmft,s_strn)
        call rx0('exit sfuned')
      endif

C ... Self-energy file editor
      if (cmdopt('--sfuned',8,0,strn)) then
        call sfuned(strn(8+1:),s_bz,s_ctrl,s_ham,s_lat,s_mix,
     .  s_move,s_pot,s_str,s_spec,s_site,s_optic,s_gw,s_dmft,s_strn)
        call rx0('exit sfuned')
      endif

      end

      subroutine sfuned(sopts,s_bz,s_ctrl,s_ham,s_lat,s_mix,
     .  s_move,s_pot,s_str,s_spec,s_site,s_optic,s_gw,s_dmft,s_strn)
C- Frequency-dependent sigma file editor
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read: *
Co     Stored:    *
Cio    Passed to: *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read: alat plat qlat npgrp nsgrp osymgr
Co     Stored:    *
Cio    Passed to: *
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read: *
Co     Stored:    *
Cio    Passed to: *
Cio  s_gw  :struct for gw-related parameters
Ci     Elts read: nkabc
Co     Stored:    *
Cio    Passed to: *
Ci Inputs
Ci   sopts :command options performed automatically, before reading
Ci         :from standard input
Cr Remarks
Cr   sfuned never returns.
Cl Local variables
Cl   agb   :holds spectral functions for multiple bands (dos branch)
Cl   lunits:0 eig,sigma,omega in Ry units
Cl         :1 eig,sigma,omega in eV units
Cl   limesh:0 => regular mesh is the original regular mesh containing SE
Cl         :1 => New regular mesh is to be constructed, interpolating sigma and maybe QP levels
Cl         :2 => regular mesh is the original regular mesh containing SE, and also
Cl         :     any of nqbnd extra points are on this mesh
Cl         :3 => no regular mesh; calculations restricted to qp at points in file
Cl         :     limesh=3 and irr<0 specify the same thing
Cl         :4 => like 2, but data may be interpolated to new mesh
Cl linputef:T when user supplies Ef
Cl   haves :0 no self-energy ready
Cl         :1 diagonal self-energy read
Cl         :2 quasiparticle hamiltonians only
Cl   havesi:0 no intepolated self-energy
Cl         :1 interpolated self-energy generated
Cl         :2 q-interpolated DOS (summed over BZ)
Cl         :11 interpolated self-energy generated on Matsubara frequencies
Cl  loptics:1 => calculate Im eps
Cl         :2 => calculate Im eps for noninteracting sigma
Cl   havepe:0 no PE spectrum
Cl         :1 PE spectrum generated
Cl   irr   :1  => qp match a regular mesh of points
Cl         :0  => qp do not match a regular mesh
Cl         :-1 => qp do not match a regular mesh; restrict calc to file q (no interpolation).
Cl         :      irr<0 specifies the same thing as limesh=3
Cl   ipqbnd:maps list of qp where spectra are calculated to qp supplied from input file.
Cl         :used when limesh is 2 or 3
Cl         :qp(iq) = qpbnd(ipqbnd(iq))
Cl   ldos  :.false. => generating information for one qp
Cl         :.true.  => generating information for multiple qp
Cl   lresib:.false. => Sum over bands.  Only spectral function or DOS may be written
Cl         :.true.  => Resolve by band.  Both sigma and spectral function may be written
Cl   lband :.true.  => generate spectral function along symmetry lines for bands
Cl   ef0   :Fermi level.  ef0=NULLI if not available.
Cl         :Note: ef0 is always in Ry units
Cl  chempot:chemical potential read from SE file.  Not used so far.  chempot=NULLI if not available.
Cl   nband :Number of one-particle states read from disk, or pared from disk read
Cl   nqbnd :number of qp where sigma, spectral functions evaluated off the regular mesh
Cl         :-1 not been assigned
Cl         : 0 no extra qp
Cl         :>0 some extra qp
Cl   qpbnd :list of nqbnd points
Cl   nqf   :number of qp in se file
Cl   qpf   :qp read from se file
Cl   eigf  :file QP levels, read with self-energies
Cl         :eigf are relative to the Fermi level.  Units depend on lunits; see above
Cl   sexf  :file Fock exchange potential, if read
Cl   vxcf  :file LDA exchange-correlation potential, if read
Cl   sigvf :File QSGW XC potential for each band, q, spin
Cl         :To satisfy the QSGW condition, sigwq(w=)=sigvf
Cl   sigwf :File self-energy(omega) for each band, q, spin.  May be shifted by schkpole?
Cl   nksabc,nqsibz,nqsfbz,ipqs,wsibz:
Cl         :analogs of nkabc,nqibz,nqfbz,ipq,qibz,wibz for self-energy interpolation
Cl   nqeibz,nqefbz,ipqe,qeibz,weibz,qefbz:
Cl         :analogs of nqibz,nqfbz,ipq,qibz,wibz for eval interpolation
Cl   nfinal:Number of final q-points to be summed over.
Cl         :Usually 1, but some number when calculating PE spectra
Cl   eigi  :eigf  interpolated to different q-mesh
Cl   sigi  :sigvf interpolated to different q-mesh
Cl   sexi  :sexf  interpolated to different q-mesh
Cl   vxci  :vxcf  interpolated to different q-mesh
Cl   sigewq:sigwf-sigi
Cl   lag0  :T => compute joint dos or optics for noninteracting g
Cl   lgetev:Calculate qp levels and use to improve interpolation
Cl   lqsgws:self-energy file has QSGW sigma
Cl   eigib :eQP for specified bands, on list of qp
Cu Updates
Cu   04 May 17 Redesigned to read DMFT-generated self-energy
Cu   30 Aug 12 Redesigned so QP eigenvalues gen by lmf improve results
Cu   17 Jul 12 First created.
C  ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      character sopts*(*)
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
C ... Dynamically allocated arrays
C     parameter (n0=10)
      integer,allocatable,target :: iblst(:),ipqbnd(:)
      integer,pointer :: iblstf(:)
      real(8),allocatable:: omgn(:),omgx(:)
      real(8),pointer:: qpf(:,:),qpbnd(:,:)
      real(8),allocatable:: sumag(:),sumag0(:),ag(:),ag0(:),dag(:),agb(:,:)
      real(8),allocatable:: dos(:,:),sag(:),sag0(:),sagkb(:,:),sag0kb(:,:),wk3(:,:,:)
      real(8),pointer:: eigiq0(:),eigi(:,:,:),sigi(:,:,:),sexi(:,:,:),vxci(:,:,:)
      complex(8),allocatable,target:: seitp(:),seb(:,:)
      complex(8),allocatable,target:: sigwf(:,:,:,:)
      complex(8),pointer:: sigewq(:,:,:,:),seitpw(:),seitpwkb(:,:)
      real(8),allocatable,target:: sigvf(:,:,:),eigf(:,:,:),sexf(:,:,:),vxcf(:,:,:),eigkb(:),cwt(:)
      real(8),pointer:: wibz(:),weibz(:)
      real(8),pointer:: qfbz(:,:),qibz(:,:),qeibz(:,:)
      complex(8),allocatable:: omgmats(:),seimats(:,:)
      real(8),pointer:: qefbz(:,:),eige(:,:,:),sige(:,:,:)
      integer,pointer :: ipq(:,:,:),ipqe(:,:,:)
C     For Pade interpolation
      integer, allocatable :: iwk(:)
      integer, allocatable :: iprm(:)  ! Used in permutation of mesh points
C     integer,pointer :: ipqs(:,:,:)
C     real(8),pointer:: wsibz(:)
      complex(8),allocatable :: znoso(:,:,:,:,:),zso(:,:,:,:,:),zwk2(:,:)
      integer, allocatable :: bmap(:)
C ... Local parameters
      logical lnsave,lband,lbin,lsopts,lse,lsem,lpe,lpe3d,lpeqp,ldos,ltet,
     .  linputef,lresib,lhf,llxc,lqsgws,lgetev
      integer ljdos,loptics
      integer i,j,i1,i2,i3,k,iw0,j1,j2,js1,js2,lascii,fmode,lunits,multi,nbas,ndham,ifiz,stdo,irr,limesh
      integer iq,ib,kb,ix(10),nsp,nspc,nspse,nomg,nqeibz,nqefbz,nblst,nblstf,nmatsubara,npade,nqbnd,onesp
      integer nband,nbandx,ndimhx,nlmto,ierr,nspx,lso,ifi,jfi,nkpan(100)
C      integer nksabc(3),nqsfbz,nqsibz
      integer npgrp,nsgrp,nqibz,isw,nqfbz,isp,ispx
      integer iwx,nomgx,npoly,jx1,jx2,iw1,iw2,ioff,nqe,nw12,nqp
      integer haves,havesi,havepe,nqf
      integer n1,n2,n3,nkabc(3),lshft(3)
      integer n1e,n2e,n3e,nkeabc(3)
      equivalence (n1,nkabc(1)),(n2,nkabc(2)),(n3,nkabc(3))
      equivalence (n1e,nkeabc(1)),(n2e,nkeabc(2)),(n3e,nkeabc(3))
      double precision qpi(3),ommin,ommax,ommaxs(3),ommaxp(3),eigitp,sigitp,Zfac,TeV,chempot
      double precision qb(3,3),alat,vol,plat(3,3),qlat(3,3),q(3),qtarg(3)
      double precision rytoev,dymx,pi,range(2),shft,ef0,dy,sdy,sdy2,ssdy,evshft,beta,sigdc
      double precision qsave(3),eigitp_save,xx(4),cwitp,cwt_save
      integer NULLI,LW5,nbpw,i1mach
      parameter (rytoev=13.60569193d0,npoly=6,dymx=1d-5,NULLI=-99999,LW5=5)
      real(8),parameter:: NULLR =-99999, a0=0.529177d0
      real(8):: tolq=1d-7,eps=5d-3/rytoev
C     Parameters needed to call bndfp
      logical lcnvtso,ltmp
      logical :: nofilter = .false.    ! if T do not filter A(w) by integrating, differentiating (see sfitrp)
      logical :: lupdateqp = .false.   ! if T do update eigqp from dominant peak in A(w)
      logical :: lchkimg               ! if T check properties of Im Sigma (e.g. causality)

      logical, parameter :: T=.true., F=.false.
      character dc*1, dc2*1, fn*120, fn2*120, outs*256, strn*120, tag*6
C     Parameters for optics
      integer nwft,npad,k1,k2,k3,ndos
C     integer cyclic
      double precision kbT,jacobian
      real(8),allocatable :: convolve(:,:,:),wk(:),fermifn(:),optmt(:,:,:,:,:)
      complex(8),allocatable :: zagb(:,:,:),zwk(:),zfermi(:) !,zconvolve(:)
!     double complex zdotc
C     Parameters for final state
      integer ifinal,nfinal
      double precision fac,ek0,sqrk0,qpar(3),qperp(3),aqpbar,radk,kfbrd,dqperp,knorm(3)
C     double precision temp_1
C     double precision Epho, Phi_s, V_0
      real(8),allocatable :: wtqf(:), qfinal(:,:), pes(:), pesi(:) !, pesk(:,:)
      double complex zer,one
      parameter (zer=(0d0,0d0),one=(1d0,0d0))

      procedure(logical) :: latvec,isanrg,cmdopt
      procedure(integer) :: fopnx,fopna,fopng,fxst,nglob,a2vec,iprint,seqitp,wordsw,mkilsd
      procedure(real(8)) :: ddot,dglob,dlength

C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
c      double precision qk
c      integer jj1,jj2,jj3,k
c      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
c     .                    (jj2*ifac(2)-1)*qb(k,2) +
c     .                    (jj3*ifac(3)-1)*qb(k,3)

C     cyclic(i,j) = i-floor(dble(i)/j)*j  ! Cyclic modulus: returns n>=0 and n<j


      pi = 4*datan(1d0)
      kbT = NULLI
      nbas = nglob('nbas')
      nsp = nglob('nsp')
      stdo = nglob('stdo')
      haves = 0  ! 0 if no self energy
                 ! 1 if file self-energy has been read
      havesi = 0 ! 0 = nothing
                 ! 1 = interpolated self-energy
                 ! 2 = dos
                 !11 = interpolated self-energy on Matsubara frequencies
      havepe = 0 ! 0 = no PE spectrum
      lunits = 0 ! Assume Ry units
      lnsave = .false.          ! .true. when local SE generated
      linputef = .false.        ! .true. when user supplies Ef
      lband = .false.           ! .true. in band generation mode
      ef0 = NULLI               ! If not NULLI, ef0 is Fermi level got by running evsync
      sigdc = NULLI             ! If not NULLI, double counting to be added to sigma
      ioff = 0                  ! Offset to 1st band in QP levels corresponding to those in se
      nullify(qfbz,qibz,wibz,ipq,sigewq,seitpw,qpf,iblstf)
      nblstf = NULLI; nblst = NULLI
      lcnvtso = F
      nband = -1
      nspc = nglob('nspc')    ! 2 for noncollinear case
      lso = isw(IAND(s_ctrl%lncol,4) /= 0)
      lgetev = .false.
      lqsgws = .false.
      lhf = .false.; llxc = .false.

C ... Lattice information and symmetry
      alat = s_lat%alat
      plat = s_lat%plat
      qlat = s_lat%qlat
      npgrp = s_lat%npgrp
      nsgrp = s_lat%nsgrp
      vol = s_lat%vol
      qpi = 0

C ... qp in full and irreducible BZ
      call info0(2,1,0,' SFUNED: generating qp for irr and full BZ from gw nkabc')
      nkabc = s_gw%nkabc
      lshft = s_gw%lshft
      s_bz%nkabc = nkabc
      s_bz%lshft = lshft
      ltet = IAND(s_ctrl%lmet,3) == 3 .or.
     .       IAND(s_ctrl%ldos,4+2+1) /= 0
      nqfbz = n1*n2*n3
      if (nqfbz /= 0) then
        if (associated(qfbz)) deallocate(qfbz,wibz,ipq,qibz)
        allocate(qfbz(3,nqfbz),wibz(nqfbz),ipq(n1,n2,n3))
        wibz(1) = 0
        call bzmesh(plat,qb,n1,n2,n3,lshft,s_lat%symgr,1,ipq,
     .    qfbz,wibz,nqibz,nqfbz,0,0)
        deallocate(wibz)
        call mkqp(s_ctrl,s_bz,s_lat,ltet,F,1,-2)
        nqibz = s_bz%nkp
        allocate(qibz(3,nqibz),wibz(nqibz))
        call icopy(nqfbz,s_bz%ipq,1,ipq,1)
        call dcopy(nqibz,s_bz%wtkp,1,wibz,1)
        call dcopy(nqibz*3,s_bz%qp,1,qibz,1)
        irr = 1
        nqefbz = nqfbz
        nqeibz = nqibz
        qefbz => qfbz
        qeibz => qibz
        weibz => wibz
        ipqe => ipq
      endif

      call pshpr(1)
      call supot(0,s_ctrl,s_lat,s_pot)
      call suham(s_ctrl,s_lat,s_spec,s_site,s_ham,s_pot,s_str,s_strn)
      ndham = s_ham%ndham
      call poppr

C      s_ctrl%plbnd = 2
C      allocate(s_ham%evals(ndham,nsp,nqibz))
C      call lmfp('sugws',s_bz,s_ctrl,s_ham,s_lat,s_mix,
C     .  s_move,s_pot,s_str,s_spec,s_site,s_optic,s_gw,s_dmft,s_strn)
C      stop

C      ifi = fopna('wkp',-1,4)
C      call getef(ifi,0,ef0)  ! ef0=999 if read fails

      dc = sopts(1:1)
      if (dc /= ' ') then
        print 301
  301   format(//' Entering the spectral function file editor. ',
     .    'Parsing command-line options ...')
        lsopts = .true.
        js2 = 0
      else
        print 302
  302   format(//' Welcome to the spectral function file editor.  ',
     .    'Enter ''?'' to see options.')
        lsopts = .false.
      endif

C ... Return here to resume parsing for arguments
   10 continue
      if (lsopts) then
        js2 = js2+1
        if (js2 > len(sopts)) then
          lsopts = .false.
          goto 10
        endif
        if (sopts(js2:js2) == dc) goto 10
        js1 = min(len(sopts),js2)
        call nwordg(sopts,0,dc,1,js1,js2)
        if (js2 < js1) lsopts = .false.
      endif

      write(stdo,"(/' Option : ')",advance='no')
      outs = ' '
      if (lsopts) then
        print '(a)', trim(sopts(js1:js2))
        outs = sopts(js1:js2)
      else
        read(*,'(a150)') outs
      endif
      if (havepe==0 .and. allocated(qfinal)) deallocate(qfinal)

      call locase(outs)

C ... Parse and execute the next command
      if (.false.) then

C --- Null input ---
      elseif (outs == ' ') then
        print 304
  304   format(' Enter ''q'' to exit, ''a'' to abort',
     .    ' ''?'' to see menu')
        goto 10

C --- Read self-energy from file ---
      elseif (outs(1:7)=='readsek' .or. outs(1:8)=='readsekb') then

C       Set up defaults
        lbin = outs(1:8)=='readsekb'
        k = 7 ; if (lbin) k = 8
        fn = 'se'; if (lbin) fn = 'seb'
        dc2 = outs(k+1:k+1)

C       See if optional switch has fn
        i = wordsw(outs,dc2,'fn=','',j1)
        if (i /= 0) then
          call nwordg(outs,0,dc2//' ',1,j1,j2)
          fn = outs(j1:j2)
        endif
        if (fxst(fn) /= 1) then
          call info0(10,0,0,"missing file '"//trim(fn)//"' ... nothing done")
          goto 10
        endif

C       Open and read file
        i = 1; if (lbin) i = 4
        ifi = fopna(fn,-1,i)
C       if (associated(qpf)) deallocate(qpf)
        if (allocated(eigf)) deallocate(eigf)
        if (allocated(sigvf)) deallocate(sigvf)
        if (allocated(sigwf)) deallocate(sigwf)
        lascii = 20; if (lbin) lascii = 0
        nblst = 0
        i = 1000 + lunits*100+lascii+4
        call ioseh(i,fmode,j,nspse,nband,nqf,nomg,ommin,ommax,chempot,nblstf,ix,ix,ifi)
        if (nspc == 2 .and. nspse /= 1) then
          call info0(10,0,0,' SE not set up for noncollinear, sorry ... nothing done')
          goto 10
        endif

        if (j /= nsp) call rx('sugws: spin mismatch in se file')
        if (nomg < 0) call rx('sugws: not ready for sigma on Matsubara frequencies')
        if (nblstf > 0) then ! Pick up iblstf
          rewind ifi
          allocate(iblstf(nblstf))
          call ioseh(4000+i,fmode,nsp,nspse,nband,nqf,nomg,ommin,ommax,chempot,nblstf,iblstf,ix,ifi)
        else
          allocate(iblstf(nband))
          forall (ib=1:nband) iblstf(ib) = ib
          nblstf = nband
        endif
        ioff = iblstf(1)-1
        if (chempot /= NULLI .and. lunits == 1) chempot = chempot/rytoev ! chempot always in Ry units
        lqsgws = mod(fmode/10,2) >= 1 ! T if have static QSGW sigma
        lhf = mod(fmode/10,4) >= 2    ! T if have Fock Exchange
        llxc = mod(fmode/10,10) >= 4  ! T if have LDA exchange-correlation potential
        nbandx = nband ! for now
C       qp, evals and QP potential at irr points given in se file
        allocate(qpf(3,nqf),eigf(nband,nqf,nspse),sigvf(nband,nqf,nspse),sigwf(nomg,nband,nqf,nspse))
        call dpzero(sigvf,size(sigvf))
        if (lhf) allocate(sexf(nband,nqf,nsp))
        if (.not. lhf) allocate(sexf(1,1,2))
        if (llxc) allocate(vxcf(nband,nqf,nsp))
        if (.not. llxc) allocate(vxcf(1,1,2))
        call iose2(lascii,fmode,fmode,nspse,nband,nqf,nomg,qpf,eigf,sigvf,sigwf,sexf,vxcf,ifi)
        call info8(10,0,-1,
     .    '%1fData from file '//trim(fn)//'%?#n# incl Vx##'//
     .    ':  nq=%i  nband=%i  nsp=%i  omega interval (%;6d,%;6d) '//
     .    '%?#n#eV#Ry# with %i points',
     .    isw(lhf),nqf,nband,nsp,ommin,ommax,lunits,nomg)
        call info2(10,0,0,'  Ef=%d',chempot,2)
        haves = 1
C       call zprm('se',2,sigwf,nomg,nomg,nband*nqf*nspse)
        lchkimg = cmdopt('--chksig',8,0,strn)
        if (lchkimg) call modsig(2,nomg,nband,nqf,nspse,ommin,ommax,sigwf)

C   ... Copy file chemical potential to ef0
        if (wordsw(outs,dc2,'useef','',j1) /= 0) then
          if (chempot == NULLI) then
            print *, 'Chemical potential not in file, sorry'
            haves = 0
            goto 10
          endif
          linputef = .true.
          ef0 = chempot
          call info2(10,0,0,' using file chemical potential, Ef=%;6d Ry = %;6d eV',ef0,ef0*rytoev)
        elseif (linputef .and. chempot /= NULLI) then ! Realign energies
          xx(1) = chempot-ef0; if (lunits == 1) xx(1) = xx(1)*rytoev
          call dvadd(eigf,1,size(eigf),xx(1))
        elseif (chempot /= NULLI) then ! Set this as the default Fermi level
          ef0 = chempot
        endif

C   ... Flag that points are not on mesh (no interpolation in this case)
        if (wordsw(outs,dc2,'irrmesh','',j1) /= 0) irr = -1

C   ... Pare list of eigenstates
        if (wordsw(outs,dc2,'ib=','',j1) /= 0) then
          call nwordg(outs,0,dc2//' ',1,j1,j2)
          nblst = mkilsd(outs(j1:j2),-1,i)
          if (nblst <= 0) call rxs(' band or null band list,',outs)
          if (allocated(iblst)) deallocate(iblst)
          allocate(iblst(nblst),iprm(nblst))
          call mkilssr(11,outs(j1:j2),nblst,iblst,[1,nband])
          if (nblst <= 0 .or. nblst > nband) goto 98
          call ilst2a(iblst,nblst,strn)
          call info2(10,1,0,' Reducing sigma to %i bands : '//trim(strn),nblst,2)

C         These are temporary work arrays
          allocate(eigi(nband,nqf,nsp),sigi(nband,nqf,nsp),sigewq(nomg,nband,nqf,nsp),
     .             sexi(nband,nqf,nsp),vxci(nband,nqf,nsp))

          call dcopy(2*size(sigwf),sigwf,1,sigewq,1)
          call dcopy(1*size(sigvf),sigvf,1,sigi,1)
          call dcopy(1*size(eigf),eigf,1,eigi,1)
          if (lhf) call dcopy(1*size(sexf),sexf,1,sexi,1)
          if (llxc) call dcopy(1*size(vxcf),vxcf,1,vxci,1)
          deallocate(eigf,sigvf,sigwf)
          allocate(eigf(nblst,nqf,nsp),sigvf(nblst,nqf,nsp),sigwf(nomg,nblst,nqf,nspse))
          if (lhf) then
            deallocate(sexf); allocate(sexf(nblst,nqf,nsp))
          endif
          if (llxc) then
            deallocate(vxcf); allocate(vxcf(nblst,nqf,nsp))
          endif
          do  isp = 1, nspse
          do  iq = 1, nqf
            do  ib = 1, nblst
C             ib = band in reduced or permuted list; kb = corresponding index from se file
              call hunti(iblstf,nblstf,iblst(ib),[0],kb); kb = kb+1
              if (iblstf(kb) /= iblst(ib)) then
                call info0(10,0,0,' band %i not present ... skipping read')
                haves = 0
                goto 10
              endif
              iprm(ib) = iblstf(kb)
              eigf(ib,iq,isp) = eigi(kb,iq,isp)
              sigvf(ib,iq,isp) = sigi(kb,iq,isp)
              sigwf(:,ib,iq,isp) = sigewq(:,kb,iq,isp)
              if (lhf) sexf(ib,iq,isp) = sexi(kb,iq,isp)
              if (llxc) vxcf(ib,iq,isp) = vxci(kb,iq,isp)
            enddo
          enddo
          enddo
          if (nblstf > 0) then
            deallocate(iblstf); allocate(iblstf(nblst))
            call icopy(nblst,iprm,1,iblstf,1)
            nblstf = nblst
            ioff = iblstf(1)-1
          endif
C         Until further checks are done, require iblst to be contiguous (note evsync!)
          do  ib = 2, nblst
            if (iblst(ib) /= iblst(ib-1)+1) call rx('ib list is not contiguous')
          enddo
          forall (ib=1:nblst) iblst(ib) = ib  ! Reorder internal band index 1,2,3,...
          deallocate(eigi,sigi,sigewq,sexi,vxci)
          nband = nblst
          nbandx = nblst
C          xx(1) = minval(eigf(1,:,:))
C          xx(2) = maxval(eigf(1,:,:))
C          xx(3) = minval(eigf(nblst,:,:))
C          xx(4) = maxval(eigf(nblst,:,:))
C          call info5(10,0,0,' min eval, band 1 = :string,a1,a2)
        endif

C   ... Printout maxmin of evals.  dos here is used as a work array
        if (wordsw(outs,dc2,'minmax','',j1) > 0) then
          allocate(dos(nband,2))
          do  isp = 1, nspse
            call info2(10,0,0,' min, max evals for all qp, spin %i',isp,2)
            forall (ib=1:nband) dos(ib,1) = minval(eigf(ib,:,isp))
            forall (ib=1:nband) dos(ib,2) = maxval(eigf(ib,:,isp))
            strn(1:3) = 'Idd'; if (nblstf > 0) strn(1:3) = 'idd'
            call arrprt('  ib     low     high','%,4i%,3;9D%,3;9D',strn(1:3),nband,0,4,
     .        0,'  | ',iblstf,dos,dos(1,2),xx,xx,xx,xx,xx)
          enddo
          deallocate(dos)
C         call rx0('done printing evals range')
        endif

        if (wordsw(outs,dc2,'dc=','',j1) > 0) then
          call nwordg(outs,0,dc2//' ',1,j1,j2)
          if (j2 <= j1) goto 99
          j = 0; j = a2vec(outs(j1:j2),len(outs(j1:j2)),j,4,' ',1,1,1,ix,sigdc)
          if (j <= 0) goto 99
          call info2(10,0,0,"%1fSubtract dc=%g from sigma",sigdc,0)
          call daxpy(size(sigwf),1d0,[-sigdc],0,sigwf,2)
        endif

C   ... Handle optional quasiparticlization of self-energy
        if (wordsw(outs,dc2,'makeqpse','',j1) > 0) then
          call qpse(1,nomg,nomg,ommin,ommax,nband,nqp,nsp,eigf,sigwf,sigvf)
        endif

        if (irr == 1) then
          if (nqf /= nqibz) then
            call info2(10,0,0,'%1f(warning) mismatch file nq (%i) '//
     .        'with bzmesh nq (%i) ... no connection to BZ established',
     .        nqf,nqibz)
            irr = 0
          endif
        endif

C   ... Match file qp to qp given on qp mesh
        if (irr == 1) then
          do  i = 1, nqf
            call findqbz(nqibz,qibz,plat,1d-5,qpf(1,i),iq)
            if (i /= iq) then
              call info2(10,0,0,'%1f(warning) mismatch file qp iq=%i '//
     .        'with bzmesh ... no connection to BZ established',i,iq)
              irr = 0
              exit
            endif
          enddo
        endif
        if (irr >= 0) then
          call info2(10,0,0,'%1ffile q-points match standard mesh? ... '
     .      //'%?#(n==1)#yes#no#',irr,0)
C        else
C          call info0(10,0,0,'%1ffile q-points not on mesh ... ')
        endif

        if (irr == 0) then
          haves = 0
          goto 10
        endif

C       Print out data showing quality of QSGW condition
        call schkpole(21,nband,nqf,nspse,nomg,ommin,ommax,eigf,sigwf,sigvf)

        goto 10

CC --- Construct quasiparticle self-energy from dynamical one---
C      elseif (outs(1:5) == 'qpse ') then

C --- Construct quasiparticle hamiltonian using lmf ---
      elseif (outs(1:6) == 'qsgwh ') then

        if (allocated(eigf)) deallocate(eigf)
        if (allocated(sigvf)) deallocate(sigvf)
        if (allocated(sigwf)) deallocate(sigwf)
        nband = ndham
        nbandx = ndham*nspc
        nspx = nsp; if (nspc == 2) nspx = 1
        nqf = nqibz
        nomg = 401
        ommin = -2*rytoev
        ommax =  2*rytoev
        allocate(eigf(nbandx,nqf,nsp))
        allocate(sigwf(nomg,nbandx,nqf,nspx),sigvf(nbandx,nqf,nspx))
        if (.not. allocated(sexf)) allocate(sexf(1,1,2))
        if (.not. allocated(vxcf)) allocate(vxcf(1,1,2))
        if (irr < 0) call rx('qsgwh not ready for irr<0')
        qpf => qibz
        sigwf = dcmplx(0d0,eps)
        call info8(10,0,0,
     .    '%1fUse data from ctrl file '//
     .    ':  nq=%i  nband=%i  nsp=%i  omega interval (%;6d,%;6d) '//
     .    '%?#n#eV#Ry# with %i points',
     .    nqf,nbandx,nsp,ommin,ommax,lunits,nomg,0)

        i = 5701+2 ! Generate and return eigf and ef0
        ioff = 0
        call evalqp(i,s_bz,s_ctrl,s_ham,s_lat,s_mix,s_move,s_pot,
     .    s_str,s_spec,s_site,s_optic,s_gw,s_dmft,s_strn,nbandx,
     .    ioff,nsp,nqf,nkabc,nqf,0,lunits,ef0,eigf,qibz,wibz,ipq,evshft)

        allocate(iblstf(nbandx))
        nblstf = nbandx
        forall (ib=1:nbandx) iblstf(ib) = ib

        haves = 2
        goto 10

C --- Recalculate QP levels ---
      elseif (outs(1:6) == 'evsync' .and. (outs(7:7) == ' ' .or. outs(7:7) == '=')) then

        call info0(10,0,0,' Recalculate QP levels to adjust sigw and sigqp on irr mesh')

        if (haves == 0) goto 97
        if (haves == 2) then
          call info0(10,0,0,' evsync makes no sense in QP context ... nothing done')
          goto 10
        endif

        if (irr == 0) then
          call info0(10,0,0,"%1fNo connection established with BZ."//
     .      "%N Modify GW_NKABC and restart program.")
          goto 10
        endif
        if (nsp /= nglob('nsp')) then
         call info2(10,0,0," Number of spins does not match se file "//
     .    "value (%i)"//"%N Reconcile spins and restart program.",nsp,0)
         goto 10
        endif

C       Find the Fermi level for a user-specified mesh
        if (outs(7:7) == '=') then
          j = 7
          j = a2vec(outs,len(outs),j,2,', ',2,-2,-3,ix,nkeabc)
          if (j < 1) goto 99
          call fill3in(j,nkeabc)
          call info2(10,0,1,' Find Fermi level from internal mesh %s,(%3i)',nkeabc,0)
          i = 8001
          call evalqp(i,s_bz,s_ctrl,s_ham,s_lat,s_mix,
     .      s_move,s_pot,s_str,s_spec,s_site,s_optic,s_gw,s_dmft,s_strn,nbandx,
     .      ioff,nsp,nqf,nkeabc,nqf,0,lunits,ef0,eigf,qibz,wibz,ipq,evshft)
          deallocate(s_ham%evals)
          linputef = .true.
        elseif (linputef) then
          call info2(10,0,1,' Use Fermi level %d',ef0,0)
        else
          call info0(10,0,1,' Find Fermi level from mesh given by se file')
        endif

        i = 7701
        if (linputef) i = i+2 ! use given ef0
        call evalqp(i,s_bz,s_ctrl,s_ham,s_lat,s_mix,
     .    s_move,s_pot,s_str,s_spec,s_site,s_optic,s_gw,s_dmft,s_strn,nbandx,
     .    ioff,nsp,nqf,nkabc,nqf,0,lunits,ef0,eigf,qibz,wibz,ipq,evshft)
        ommin = ommin - evshft; ommax = ommax - evshft
        call info5(10,1,0,' sugws: shifting omega interval by %d to '//
     .    '(%;6d,%;6d)',evshft,ommin,ommax,0,0)
C       Enforce QSGW condition
        if (lqsgws) then
          i=31
          if (iprint() >= 45) i=32
          call schkpole(i,nbandx,nqf,nsp,nomg,ommin,ommax,eigf,sigwf,sigvf)
C         call schkpole(21,nbandx,nqf,nsp,nomg,ommin,ommax,eigf,sigwf,sigvf)
        endif

C   ... Convert sig into relativistic form
        call word(outs,2,j1,j2)
        if (j2 < j1) goto 10
        if (outs(j1:j2) /= 'so') goto 10
        lcnvtso = T
        if (lcnvtso) then
          ndimhx = ndham*2
C         Read nonrel evecs
          ifiz = fopna('evec',-1,4)
          rewind ifiz
          call iosigh(3,LW5,nsp,nspc,ndham,nlmto,n1,n2,n3,nqf,iq,lshft(1),lshft(2),lshft(3),ifiz,xx)
          allocate(znoso(ndham,2,ndham,2,nqf),zwk2(ndham,ndham))
          call dpzero(znoso,2*ndham*2*ndham*2*nqf)
          allocate(s_ham%evals(ndham,nsp,nqf))
          do  iq = 1, nqf
            do  isp = 1, nsp
C             Sanity check: qp must match
              read(ifiz) q
              if (.not. latvec(1,tolq,qlat,q-qibz(1:3,iq))) then
                write(stdo,*) iq
                write(stdo,*) sngl(q)
                write(stdo,*) sngl(qibz(:,iq))
                call rx('incompatible q-mesh')
              endif
              call dpdump(s_ham%evals(1,isp,iq),ndham,ifiz)
              call dpdump(zwk2,ndham**2*2,ifiz)
              znoso(1:ndham,isp,1:ndham,isp,iq) = zwk2(1:ndham,1:ndham)
            enddo
          enddo
C         Permute the eigenvectors
          nbpw = int(dlog(dble(i1mach(9))+1d0)/dlog(2d0))
          allocate(bmap((ndham*nsp*nqf/nbpw+1)))
          call iinit(bmap,(ndham*nsp*nqf/nbpw+1))
          call ebcpl(0,ndham,ndham,nsp,1,nqf,nbpw,bmap,zwk2,s_ham%evals)
C         call prmx('permuted evl',s_ham%evals,ndimhx,ndimhx,nqf)
C          iq = 1
C          call zprm('evec, noso',2,znoso(1,1,1,1,iq),
C     .        ndimhx,ndimhx,ndimhx)
          call ebcpl2(10,ndimhx,ndham,ndham,nsp,1,nqf,nbpw,bmap,xx,znoso)
C          call zprm('evec, noso',2,znoso(1,1,1,1,iq),
C     .        ndimhx,ndimhx,ndimhx)
          deallocate(zwk2)
          call fclose(ifiz)

C         Generate relativistic evecs and read from disk
          i = 5701
          if (linputef) i = i+2 ! use given ef0
          xx(1) = dglob('nspc',2d0,1)
          xx(1) = dglob('nspc3',2d0,1)
          s_ctrl%lncol = 4
          s_ham%lncol = 4
          call evalqp(i,s_bz,s_ctrl,s_ham,s_lat,s_mix,s_move,s_pot,
     .      s_str,s_spec,s_site,s_optic,s_gw,s_dmft,s_strn,nbandx,ioff,nsp,nqf,
     .      nkabc,nqf,0,lunits,ef0,eigf,qibz,wibz,ipq,evshft)
          ifiz = fopna('evec',-1,4)
          rewind ifiz
          call iosigh(3,LW5,nsp,nspc,ndham,nlmto,n1,n2,n3,nqf,iq,lshft(1),lshft(2),lshft(3),ifiz,xx)
          allocate(zso(ndham,2,ndham,2,nqf))
          do  iq = 1, nqf
C           Sanity check: qp must match
            read(ifiz) q
            if (.not. latvec(1,tolq,plat,q-qibz(1:3,iq))) then
              write(stdo,*) iq
              write(stdo,*) sngl(q)
              write(stdo,*) sngl(qibz(:,iq))
              call rx('incompatible q-mesh')
            endif
            call dpdump(xx,1,ifiz) ! Skip over evals
            call dpdump(zso(1,1,1,1,iq),(ndham*2)**2*2,ifiz)
          enddo

          allocate(zwk2(ndimhx,ndimhx+1))
          do  iq = 1, nqf
            call zprm('evec, noso',2,znoso(1,1,1,1,iq),
     .        ndimhx,ndimhx,ndimhx)
            call zprm('evec, so',2,zso(1,1,1,1,iq),
     .        ndimhx,ndimhx,ndimhx)
            call zqinv('N',zso(1,1,1,1,iq),ndimhx,0,ndimhx,zwk2,
     .        ndimhx,ierr)
            call rxx(ierr<0,'rotevs: failed to invert zso')
            call zgemm('N','N',ndimhx,ndimhx,ndimhx,one,
     .        znoso(1,1,1,1,iq),ndimhx,zso(1,1,1,1,iq),ndimhx,
     .        zer,zwk2,ndimhx)
            call zprm('U',2,zwk2,ndimhx,ndimhx,ndimhx)
          enddo

          stop 'not ready'
        endif

        goto 10

C --- evsync[flag] extract QP levels from peaks in spectral function ---
      elseif (outs(1:6) == 'evsync') then
        dc2 = outs(7:7)
        if (wordsw(outs,dc2,'sf',dc//dc2//' ',j1) > 0) then
          if (haves /= 1) goto 97
          if (allocated(iblst)) deallocate(iblst)
          if (wordsw(outs,dc2,'ib=','',j1) > 0) then
            j2 = scan(outs(j1:),dc//dc2//' ') + j1-2
            nblst = mkilsd(outs(j1:j2),-1,[i])
            if (nblst <= 0) goto 99
            allocate(iblst(nblst))
            if (mkilsd(outs(j1:j2),nblst,iblst) < 0) call rx('sugws bug in mkilsd')
          else
            allocate(iblst(nband))
            nblst = nband
            forall (ib=1:nband) iblst(ib) = ib
          endif
          dy = -NULLI
          if (wordsw(outs,dc2,'demax=','',j1) > 0) then
            j = j1-1
            if (a2vec(outs,len(outs),j,4,dc//dc2//' ',3,1,1,ix,dy) /= 1) goto 99
          endif
          sdy = 1
          if (wordsw(outs,dc2,'acut=','',j1) > 0) then
            j = j1-1
            if (a2vec(outs,len(outs),j,4,dc//dc2//' ',3,1,1,ix,sdy) /= 1) goto 99
          endif
          range(1) = NULLR; range(2) = -NULLR
          if (wordsw(outs,dc2,'range=','',j1) > 0) then
            j = j1-1
            if (a2vec(outs,len(outs),j,4,', '//dc//dc2,4,-2,-3,ix,range) /= 2) goto 99
          endif


          ltmp = wordsw(outs,dc2,'warn','',j1) > 0

          if (allocated(omgn)) deallocate(omgn,omgx); allocate(omgn(nomg),omgx(nomg))
          call sfitrp(1,nomg,nomg,ommin,ommax,1,xx,xx,xx,omgn,omgx,xx,xx,xx,xx,xx,xx,ix,ix)

C         allocate(seitp(nomg),sag(nomg),sag0(nomg),dag(nomg),seitp(nomg),seitpw(nomg))
          allocate(ag(nomg),ag0(nomg),sumag(nomg),sumag0(nomg),seitp(nomg),seitpw(nomg),dag(nomg))

          call info0(10,0,0,
     .    '  ib Npeak   eQP       w(Amax)     shift%6fZ%7fFWHM     w(A2)%7fshift   A2/Amax     w*    rank')

          ssdy = 0
          do  isp = 1, nsp
          do  iq = 1, nqf
          do  kb = 1, nblst

            ib = iblst(kb)

            if (eigf(ib,iq,isp) < range(1) .or. eigf(ib,iq,isp) > range(2)) cycle

            forall (i=1:nomg) seitp(i) = sigwf(i,ib,iq,isp) - sigvf(ib,iq,isp)

            call sfitrp(2,nomg,nomg,ommin,ommax,1,eigf(ib,iq,isp),seitp,
     .        eps,omgn,omgn,seitpw,sumag,sumag0,ag,ag0,Zfac,jx1,jx2)

C       ... Find and print to stdout largest peaks in A
            shft = 0
            call spec_peaks(1,nomg,nomg,omgn,ifi,ib,shft,eigf(ib,iq,isp),ag,dag,Zfac,ommaxs)
            if (abs(ommaxs(1)-eigf(ib,iq,isp)) < dy) cycle ! Main peak within tolerance
            if (abs(ommaxs(2)-eigf(ib,iq,isp)) < dy .and. ommaxs(3)>=sdy) cycle ! Nearest peak large enough

            ssdy = max(ssdy,abs(eigf(ib,iq,isp)-ommaxs(2)))

            call awrit7(' shift %;7F in eig %;7F for ib,iq,isp = %,3i%,3i%,3i exceeds tolerance %;7F.  A(nn)/A(max)=%;3g',
     .        strn,len(strn),0,ommaxs(1)-eigf(ib,iq,isp),eigf(ib,iq,isp),ib,iq,isp,dy,ommaxs(3))
            if (.not. ltmp) call rx('sugws: QP '//trim(strn))
            call info0(20,0,0,strn)

C           Debugging
            allocate(dos(nomg,2))
            call dcopy(nomg,omgn,1,dos,1)
            call dcopy(nomg,ag,1,dos(1,2),1)
            call prmx('A(w)',dos,nomg,nomg,2)
            deallocate(dos)

          enddo
          enddo
          enddo

          call info2(20,1,0,' '//trim(outs)//': maximum shift = %;7F',ssdy,0)
          stop 'here'

        else
          goto 98
        endif

C --- Set smearing ---
      elseif (outs(1:4) == 'eps ' .or. outs(1:4) == 'eps=') then

        if (outs(1:4) == 'eps ') then
          call word(outs,2,j1,j2)
        else
          call word(outs,1,j1,j2)
          j1 = j1+4
        endif

        if (j2 < j1) goto 99
        j = 0; j = a2vec(outs(j1:j2),len(outs(j1:j2)),j,4,' ',1,1,1,ix,eps)
        if (j <= 0) goto 99
        call info2(10,0,0,'%1fSmearing set to %g %?#n#eV#Ry#',eps,lunits)
        goto 10

C --- Set Fermi ---
      elseif (outs(1:3) == 'ef ' .or. outs(1:3) == 'ef=') then

        if (outs(1:3) == 'ef ') then
          call word(outs,2,j1,j2)
        else
          call word(outs,1,j1,j2)
          j1 = j1+3
        endif

        if (j2 <= j1) goto 99
        j = 0
        j = a2vec(outs(j1:j2),len(outs(j1:j2)),j,4,' ',1,1,1,ix,ef0)
        if (j <= 0) then
          ef0 = NULLI
          goto 99
        endif
        call info2(10,0,0,"%1fFermi level set to %g Ry",ef0,0)
        linputef = .true.
        goto 10

C --- Set units ---
      elseif (outs(1:6) == 'units' .or. outs(1:6) == 'units=') then

        if (outs(1:6) == 'units') then
          call word(outs,2,j1,j2)
        else
          call word(outs,1,j1,j2)
          j1 = j1+6
        endif

        if (j2 <= j1) then
          goto 99
        elseif (outs(j1:j2) == 'ry' .or. outs(j1:j2) == 'Ry') then
          lunits = 0
          call info0(10,0,0,' Units set to Ry')
        elseif (outs(j1:j2) == 'ev' .or. outs(j1:j2) == 'eV') then
          lunits = 1
          call info0(10,0,0,' Units set to eV')
        else
          goto 99
        endif

C --- Interpolate self-energy, spectral functions to arbitrary q, or make DOS ---
      elseif (outs(1:2) == 'se' .or.
     .        outs(1:2) == 'pe' .or.
     .        outs(1:3) == 'sem' .or.
     .        outs(1:3) == 'dos' .or.
     .        outs(1:4) == 'pe3d' .or.
     .        outs(1:4) == 'peqp' .or.
     .        outs(1:4) == 'jdos' .or.
     .        outs(1:5) == 'imeps') then
        if (haves /= 1 .and. outs(1:4) /= 'peqp') then
          call info0(10,0,0,"%1fSelf-energy must be read first.")
          goto 10
        endif

        dc2 = outs(3:3)
        lpe  = outs(1:2) == 'pe'   ! T if seeking PE
        lse  = outs(1:2) == 'se'   ! T if seeking SE
        ldos = outs(1:3) == 'dos'  ! T if seeking DOS.  More generally ldos set when looping over qp
        lsem = outs(1:3) == 'sem'  ! T if seeking SE on Matsubara
        lpe3d= outs(1:4) == 'pe3d' ! T if seeking PEqp (qp spectral function, e.g. with SO)
        lpeqp= outs(1:4) == 'peqp' ! T if seeking PEqp (qp spectral function, e.g. with SO)
        ljdos   = isw(outs(1:4)== 'jdos') !  1 if seeking joint DOS
        loptics = isw(outs(1:5)== 'imeps') ! 1 if seeking optics
        if (ldos .or. lsem) dc2 = outs(4:4)
        if (lpe3d .or. lpeqp .or. ljdos>0) dc2 = outs(5:5)
        if (loptics>0) dc2 = outs(6:6)
        if (lpe3d) lpe = .false.
        if (lpe3d .or. lpe .or. lpeqp) lse = .true.
        nfinal = 0  ! Default
        TeV = .0259d0
        kfbrd = 0.19

C   ... Simulate PE by integrating along a line
        if (lpe .or. lpeqp) then
C         Define constant for final state
          nfinal = 200
C         dqperp = 0.003 * (2*pi)/(alat*a0)  ! Julien's convention for now
          dqperp = 0.003 * (alat*a0)/(2*pi) ! Julien's convention for now
C         Define the perpendicular direction
          knorm(1)=1
          knorm(2)=1
          knorm(3)=0
C   ...   Define the k-point path for the band structure calculation
C         Epho = kinetic energy of the photon
C         Phi_s = work function of the sample
C         V_0 = inner potential
C          Epho=139.0
C          Phi_s=5.0
C          V_0=14.8
C         temp_1=0.5123*SQRT(Epho-Phi_s+V_0) ! 0.5123^2 = 2m/hbar^2 eV^-1 A^-2
C         Ek0 = Ekin + V0 = Ephoton-Phi_s+V_0; see Notes
C          Ek0=139d0-5d0+14.8d0
C          if (lunits == 1) then
C            sqrk0 = Ek0/rytoev  ! k0^2 in atomic units
C            sqrk0 = sqrk0*(alat/2/pi)**2 ! In units of 2*pi/alat
C          endif
          Ek0 = NULLI
        endif

C   ... Construct nqbnd and qp list of points to evaluate
        nqbnd = -1   ! number of non-mesh qp where sigma, A made.  Set here when it can be determined
        nkpan(1) = 0 ! number of panels for special band mode
        if (ldos .or. loptics>0 .or. ljdos>0) then
          nqbnd = 0
          ldos = .true.         ! use 'dos' branch
        elseif (lse .or. lsem .or. lpeqp) then ! SE specific
          if ((lse .or. lsem) .and. wordsw(outs,dc2,'allq',dc//dc2//' ',j1) > 0) then
            if (lpe3d) goto 98
            ldos = .true.       ! use 'dos' branch, points on a mesh to be determined later
            nqbnd = 0           ! points will be all be on a mesh
          elseif ((lse .or. lsem) .and. wordsw(outs,dc2,'band','',j1) > 0) then
            if (lpe3d) goto 98
            ldos = .true.       ! use 'dos' branch
            lband = .true.      ! in band generation mode
            call info0(30,1,0,' ... make A(w) and sigma(w) for qpts read from file')
            j2 = scan(outs(j1:),dc//dc2//' ') + j1-2
            call rdqlst(0,s_lat,outs(j1:j2),onesp,nqbnd,xx,xx)
            allocate(qpbnd(3,nqbnd))
            call rdqlst(3,s_lat,outs(j1:j2),onesp,nqbnd,qpbnd,nkpan)
            if (nqbnd == 0) call rx('sugws: no qp specified with band mode')
          elseif (wordsw(outs,dc2,'iq=','',j1) > 0) then
            nqbnd = 1           ! just one qp
            j = j1-1
            if (a2vec(outs,len(outs),j,2,dc//dc2//' ',3,1,1,ix,iq) /= 1) goto 99
            if (iq > nqf) goto 98
            allocate(qpbnd(3,1))
            ltmp = isanrg(iq,1,nqf,' sfuned: ','iq',.true.)
            qpi(1:3) = qpf(1:3,iq) ! keep since qpbnd may be reallocated
            qpbnd(1:3,1) = qpf(1:3,iq)
          elseif (wordsw(outs,dc2,'q=','',j1) > 0) then
C            if (lpe3d) call info0(30,1,0,' ... for pe3d specify qp through iq')
C            if (lpe3d) goto 98  ! No interpolation for now
            nqbnd = 1           ! just one qp
            j = j1-1
            if (a2vec(outs,len(outs),j,4,', '//dc//dc2,4,2,3,ix,qpi) /= 3) goto 99
            allocate(qpbnd(3,1))
            qpbnd(1:3,1) = qpi(1:3)
          else
            goto 99
          endif
        endif
        qtarg = qpi

C   ... Band list
        lresib = .false.
        ltmp = wordsw(outs,dc2,'ib=','',j1) > 0
        if (.not. ltmp) ltmp = wordsw(outs,dc2,'ibx=','',j1) > 0
        if (ltmp) then
          lresib = outs(j1-2:j1-2) == 'x'  ! Resolve by band
          if (.not. lresib .or. lse .and. .not. (lpe .or. lpeqp)) then !lresib only for SE
          else
            call info0(0,0,0,' ibx not allowed with this mode, sorry')
            goto 98
          endif
          j2 = scan(outs(j1:),dc//dc2//' ') + j1-2
          nblst = mkilsd(outs(j1:j2),-1,[i])
          if (nblst <= 0) goto 99
          if (allocated(iblst)) deallocate(iblst); allocate(iblst(nblst))
          if (mkilsd(outs(j1:j2),nblst,iblst) < 0) call rx('sugws bug in mkilsd')
          if (.not. associated(iblstf)) call rx('iblstf not allocated')
          call imxmn(nblstf,iblstf,1,j1,j2)
          call imxmn(nblst,iblst,1,i,j)
          if (i < j1 .or. j > j2) then
            call info2(0,0,0,' oops ... elements in ib must be in range (%i,%i)',j1,j2)
            goto 98
          endif
          if (iblstf(1) /= 1) call rx('se bands not yet synchronized with selected ones')
          deallocate(iblstf); nblstf = nblst; iblstf => iblst
        endif

        if (lsem) then  ! parameters specific to SEM
          tag = 'beta='
          if (.not. wordsw(outs,dc2,trim(tag),'',j1) > 0) goto 96
          j = j1-1
          if (a2vec(outs,len(outs),j,4,' '//dc//dc2,3,1,1,ix,beta) /= 1) goto 99

          tag = 'nm='
          if (.not. wordsw(outs,dc2,trim(tag),'',j1) > 0) goto 96
          j = j1-1
          if (a2vec(outs,len(outs),j,2,' '//dc//dc2,3,1,1,ix,nmatsubara) /= 1) goto 99
        endif

C   ... Optional switches
        multi = 1; isp = 1
        range(1) = NULLI; range(2) = -range(1)
C       nksabc = nkabc
        nkeabc = 0              ! flag that no interpolating mesh params nkeabc read

        if (wordsw(outs,dc2,'getev','',j1) > 0) then
          lgetev = .true.
          if (outs(j1:j1) == '=') then
            j = j1
            j = a2vec(outs,len(outs),j,2,', '//dc//dc2,4,2,3,ix,nkeabc)
            if (j < 1) goto 99
            call fill3in(j,nkeabc)
          endif
        endif

        lupdateqp = wordsw(outs,dc2,'a2qp','',j1) > 0

        nofilter = wordsw(outs,dc2,'rawnw=','',j1) > 0
        if (.not. nofilter) ltmp = wordsw(outs,dc2,'nw=','',j1) > 0
        if (nofilter .or. ltmp) then
          j = j1-1
          if (a2vec(outs,len(outs),j,2,' '//dc//dc2,3,1,1,ix,multi) /= 1) goto 99
        endif

        if (wordsw(outs,dc2,'a0','',j1) > 0) then
          if (loptics>0) then
            loptics = 2
          elseif (ljdos>0) then
            ljdos = 2
          else
            call info0(1,0,0,' ignoring a0 ...')
          endif
        endif

        if (wordsw(outs,dc2,'kt=','',j1)>0) then
          if (loptics>0 .or. ljdos>0) then
            j = j1-1
            if (a2vec(outs,len(outs),j,4,' '//dc//dc2,3,1,1,ix,kbT) /= 1) goto 99
          else
            call info0(1,0,0,' ignoring kt= ...')
          endif
        endif

        if (wordsw(outs,dc2,'domg=','',j1) > 0) then
          print *, 'SE not ready for domg= option ...'; goto 99
        endif

        if (wordsw(outs,dc2,'nqf=','',j1) > 0) then
          if (lpe .or. lpeqp) then
            j = j1-1
            if (a2vec(outs,len(outs),j,2,' '//dc//dc2,3,1,1,ix,nfinal) /= 1) goto 99
          else
            call info0(1,0,0,' ignoring nqf= ...')
          endif
        endif

        if (wordsw(outs,dc2,'ke0=','',j1) > 0) then
          if (lpe3d .or. lpe .or. lpeqp) then
            j = j1-1
            if (a2vec(outs,len(outs),j,4,' '//dc//dc2,3,1,1,ix,Ek0) /= 1) goto 99
            if (lunits == 1) then
              sqrk0 = Ek0/rytoev ! k0^2 in atomic units
              sqrk0 = sqrk0*(alat/2/pi)**2 ! In units of 2*pi/alat
            endif
          else
            call info0(1,0,0,' ignoring ke0= ...')
          endif
        endif

        if (wordsw(outs,dc2,'delk=','',j1) > 0) then
          if (lpe3d .or. lpe .or. lpeqp) then
            j = j1-1
            if (a2vec(outs,len(outs),j,4,' '//dc//dc2,3,1,1,ix,kfbrd) /= 1) goto 99
          else
            call info0(1,0,0,' ignoring kfbrd= ...')
          endif
        endif

        if (wordsw(outs,dc2,'isp=','',j1) > 0) then
          if (loptics>0 .or. ljdos>0) then
            call info0(1,0,0,' ignoring isp= ...')
          else
            j = j1-1
            if (a2vec(outs,len(outs),j,2,' '//dc//dc2,3,1,1,ix,isp) /= 1) goto 99
            if (isanrg(isp,1,nsp,' sfuned: ','isp',.false.)) goto 99
            if (isp < 1 .or. isp > 2 .or. isp == 2 .and. nspc == 2) then
              call info0(1,0,0,'isp=2 but noncollinear ...')
              goto 98
            endif
          endif
        endif

        if (wordsw(outs,dc2,'range=','',j1) > 0) then
          j = j1-1
          if (a2vec(outs,len(outs),j,4,', '//dc//dc2,4,-2,-3,ix,range) /= 2) goto 99
        endif

        if (ldos .and. wordsw(outs,dc2,'nq=','',j1) > 0) then
          j = j1-1
          j = a2vec(outs,len(outs),j,2,', '//dc//dc2,4,2,3,ix,nkeabc)
          if (j < 1) goto 99
          call fill3in(j,nkeabc)
        endif

        if (wordsw(outs,dc2,'radk=','',j1) > 0) then
          if (.not. lpe3d) then
            call info0(1,0,0,' ignoring radk= ...')
          else
            j = j1-1
            if (a2vec(outs,len(outs),j,4,' '//dc//dc2,3,1,1,ix,radk) /= 1) goto 99
          endif
        elseif (lpe3d) then
          call info0(1,0,0,' oops! pe3d requires radk= ...'); goto 98
        endif

C   ... This ends parsing of optional switches.  Sanity checks
        if (irr == 0) then
          call info0(10,0,0,"%1fNo connection established with BZ."//
     .      "%N Modify GW_NKABC and restart program.")
          goto 10
        endif
C       Redo kfbrd in units of 2*pi/alat
        kfbrd = kfbrd * (alat*a0)/(2*pi)

C   ... Distinguish new interpolating mesh, and mesh corresponding to given spectrum mesh
        if (irr < 0) then  ! No regular mesh
          limesh = 3
        elseif (nkeabc(1) == 0) then ! points on the original regular mesh
          limesh = 0
          nkeabc = nkabc
        else ! A new regular mesh is to be constructed, interpolating sigma
          limesh = 1
        endif

C       Fill out band list iblst
        if (ldos) then
C          nqbnd = 0       ! No extra points in usual dos mode
C          if (lband) then
C            nqbnd = nqbnd ! But may be extra points in band-plotting mode
C            allocate(qfinal(3,nqbnd))
C            call dcopy(3*nqbnd,qpbnd,1,qfinal,1)
C          endif
          if (nblst <= 0) then
            if (allocated(iblst)) deallocate(iblst)
            allocate(iblst(nbandx))
            nblst = nbandx
            forall (ib=1:nband) iblst(ib) = ib
          endif
        else
C          stop '1059 next 3 lines should not be needed'
C          nqbnd = 1
C          allocate(qfinal(3,1))
C          qfinal(:,1) = qpi
          if (range(1) /= NULLI) then
            call info5(10,1,0,' Spectral function at %s, q=(%3;6d)'//
     .        ' window (%;3,3d,%;3,3d) dw=%;5,5d eps=%;5G',
     .        qpi,range(1),range(2),(ommax-ommin)/((nomg-1)*multi),eps)
          else
            call info5(10,1,0,' Spectral function at %s, q=(%3;6d)'//
     .        ' window (%;3,3d,%;3,3d) dw=%;5,5d eps=%;5G',
     .        qpi,ommin,ommax,(ommax-ommin)/((nomg-1)*multi),eps)
          endif
        endif

C       PE simulation (1D): make list of final states -> qpbnd
        if (lpe .or. lpeqp) then
C          if (ek0 == NULLI) then
C            call info0(0,0,0,' sfuned:  '//
C     .        'must supply ke0 for PE spectra ... nothing done')
C            goto 10
C          endif
          qsave = qpi
C         KE correction for q
          if (ek0 /= NULLI) then
            dy = ddot(3,knorm,1,qpi,1)/ddot(3,knorm,1,knorm,1)
            qpar(:) = qpi(:) - dy*knorm(:) ! projection || surface
            qperp(:) = qpi(:) - qpar(:) ! projection perp surface
            aqpbar = dsqrt(sqrk0 - ddot(3,qpar,1,qpar,1)) ! |bar kperp|
            qpi(:) = qpar(:) - knorm(:)* ! I think should be minus (+ for Julien)
     .        (aqpbar - dlength(3,qperp,1))/dlength(3,knorm,1)
            call shorps(1,qlat,[72,2,2],qpi,qpi)
            call info8(10,1,0,' PES:  KE+V0=%;3d  kperp=%3:-1;6,6d'//
     .        '  q(internal)=%3:1;6,6d%N'//
     .        ' kf:  broadening=%;3d=%;3dA^-1  range=+/-%;2gA^-1 (%i pts)'//
     .        '  kT=%;3geV',
     .        Ek0,knorm,qpi,kfbrd,kfbrd*(2*pi)/(alat*a0),
     .        dlength(3,nfinal/2*dqperp*knorm(:),1),nfinal,TeV)
          else
            call info5(10,1,0,' PES:  no KE correction of q.'//
     .        '  kf:  broadening=%;3dA^-1  range=+/-%;2gA^-1 (%i pts)'//
     .        '  kT=%;3geV',
     .        kfbrd*(2*pi)/(alat*a0),
     .        dlength(3,nfinal/2*dqperp*knorm(:),1),nfinal,TeV,0)
          endif

          nqbnd = nfinal+1
          allocate(qpbnd(3,nfinal+1),wtqf(nfinal+1))
          i = 1
          qpbnd(:,1) = qpi
C         Integration over final state broadened by kfbrd:
C         The following integral
C         (kfbrd/2)/pi * int_-infty^infty A(qf) 1/((kfbrd/2)^2 + (qf-qf0)^2) d(qf)
C         evaluates to 1 when A(qf) = 1
C         If quadrature is a uniform mesh of points spaced dqperp*|knorm|:
C         wtqf(i) = (kfbrd/2)/pi * dqperp*|knorm| /((kfbrd/2)^2 + qf^2)
          fac = (kfbrd/2)/pi*dqperp*dlength(3,knorm,1)
          wtqf(1) = fac/(kfbrd/2)**2
          dy = wtqf(1)
          do  ifinal = -nfinal/2, nfinal/2
            if (ifinal == 0) cycle
            i = i+1
            qpbnd(:,i) = ifinal*dqperp*knorm(:) + qpbnd(:,1)
            wtqf(i) = fac / ((kfbrd/2)**2 + dlength(3,ifinal*dqperp*knorm(:),1)**2)
            dy = dy + wtqf(i)
          enddo

          call info8(30,1,0,' Range of integration over kfinal = '//
     .      '(%,3;3d,%,3;3d)*(%;3d,%;3d,%;3d).%N'//
     .      ' Quadrature of %i+1 points sums to %d ...'//
     .      ' scale to normalize',
     .      -nfinal/2*dqperp,nfinal/2*dqperp,knorm(1),knorm(2),
     .      knorm(3),nfinal,dy,0)
          call dscal(nfinal+1,1/dy,wtqf,1)

C          ifi = fopna('qpbnd',-1,0) ! Dump list of qpfinal
C          do  i = 1, nfinal
C            write(ifi,345) qpbnd(:,i), wtqf(i)
C  345       format(3f12.6,2x,f12.6)
C          enddo
C          call rx0('done')
        endif
C       End of assembling qpbnd.
C       irr>0 : Interpolate sigma and/or A at nqbnd points
C       irr<0 : Interpolate sigma in frequency but not k
        if (nqbnd < 0) call rx('no nqbnd')

C   ... Set limesh=2 and ipqbnd if 'extra' points all coincide with points read from file
        if (nqbnd > 0 .and. limesh == 0) then
C          allocate(ipqbnd(nqbnd)); call iinit(ipqbnd,size(ipqbnd))
C          do  iq = 1, nqbnd
C            do  i = 1, nqibz
C              if (latvec(1,tolq,plat,qpbnd(1:3,iq)-qibz(1:3,i))) then
C                ipqbnd(iq) = i
C                exit
C              endif
C            enddo
C            if (ipqbnd(iq) == 0) exit
C            if (iq == nqbnd) limesh = 2 ! Extra points are all on regular mesh
C          enddo
          allocate(ipqbnd(nqbnd))
          call alignqp(nqbnd,qpbnd,nqibz,qibz,plat,tolq,ipqbnd,i)
          if (i == 0) limesh = 2
        endif

C   ... Make vector of qp and corresponding evals
C      *Case irr > 0:
C   ... Make (1) special points 1:nqbnd on a mesh, from qpbnd
C            (2) points nqbnd+1:nqbnd+neibz interpolated q mesh
C            (3) evals for all these points if (lgetev)
C       Eventually need either (1) only or (2) only.
C       But when interpolating a mesh we always need (2) initially to refine the fermi level
C       After refinement, merge (2) into (1)
C       This block makes:
C       nqeibz = number of points in the interpolated mesh of the irr BZ
C       nqefbz = number of points in the interpolated mesh of the full BZ + nqbnd
C                Note nqbnd added here because evals stored at these points for eval "interpolation"
C       qefbz  = corresponding qp
C       nqe    = nqeibz + nqbnd
C       eige   = (lgetev=T) evals at the nqe points.
C       qeibz  = corresponding qp
C       weibz  = weights for irr zone, interpolated mesh.  Note weibz(1:nqbnd) = 0
C       ipqe   = map full BZ + nqbnd to nqe points in irr zone + nqbnd
C      *Case irr < 0:
C       nqeibz = nqe = nqefbz = nqbnd
C       qeibz => qpbnd
C       eige, iqpe, qefbz, weibz are no longer defined

        call info0(10,0,0,' ')
C   ... Setup for simulations when file self-energy lies on a regular mesh
        if (lpe3d .or. (irr >= 0 .and. (limesh == 1 .or. nqbnd /= 0 .or. lgetev))) then
          if (irr < 0) call rx('irr<0 not ready for getev')

C         Make qeibz, nqeibz = irr qp and number on new regular mesh, and ipqe
          call pshpr(1)
          i = nkeabc(1)*nkeabc(2)*nkeabc(3)
          allocate(qeibz(3,i),weibz(i),ipqe(nkeabc(1),nkeabc(2),nkeabc(3)))
          call bzmshp('sugws',3,nkeabc,lshft,plat,s_lat%symgr,nsgrp,
     .      nsgrp,.false.,0,0,qb,nqeibz,qeibz,weibz,ipqe,xx)
          deallocate(weibz)
          call poppr

C         Set limesh=4 and ipqbnd if 'extra' points all coincide with points on new regular mesh
          if (limesh /= 2 .and. limesh /= 3 .and. nqbnd > 0) then
            if (allocated(ipqbnd)) deallocate(ipqbnd); allocate(ipqbnd(nqbnd))
            call alignqp(nqbnd,qpbnd,nqeibz,qeibz,plat,tolq,ipqbnd,i)
            if (i == 0) limesh = 4
          endif
          deallocate(qeibz)

C         limesh == 2 or 4 implies "extra points" are on regular mesh read from disk.
C         Extra points will be given by nfinal and ipqbnd which references the mesh points
          if (limesh == 2 .or. limesh == 4) then
            nqbnd = 0
          endif
          nqe = nqeibz+nqbnd  ! This is final value for nqe

C         Make qefbz = (qp on new mesh in full BZ) + nband extra qp
          nqefbz = n1e*n2e*n3e + nqbnd
          allocate(qefbz(3,nqefbz),weibz(nqefbz),iwk(nqefbz))
          if (nqbnd > 0) call dcopy(3*nqbnd,qpbnd,1,qefbz,1)
          weibz(1) = 0
          call pshpr(1)
          call bzmesh(plat,qb,n1e,n2e,n3e,lshft,s_lat%symgr,1,iwk,qefbz(1,1+nqbnd),weibz,i,n1e*n2e*n3e,0,0)
          call poppr
          deallocate(weibz,iwk)
          if (i /= n1e*n2e*n3e) call rx('bug in bzmesh')

C     ... PES, 3D: find all qp in full BZ inside radius of qpi.
C         This block requires nqbnd = 0 for qefbz, ipqe to be aligned
          if (lpe3d) then
            if (limesh /= 2 .and. limesh /= 4) call rx('bug in sugws, pes3d')
            nfinal = 1
            do  i = 1, 2
              deallocate(ipqbnd); allocate(ipqbnd(nfinal)); nfinal = 1
              do  iq = 1, nqefbz
                call shorps(1,qlat,[72,2,2],qefbz(1:3,iq)-qpi(:),q)
C               call shorbz(qefbz(1:3,iq)-qpi(:),q,qlat,plat)
                if (dlength(3,q,1) == 0) then
                  ipqbnd(1) = iq  ! ipqbnd is always allocated at least length 1
                elseif (dlength(3,q,1) <= radk) then
                  nfinal = nfinal+1
                  if (i == 2) ipqbnd(nfinal) = iq
C                  print 864, nfinal,iq,qefbz(:,iq),q,dlength(3,q,1)
C  864             format(2i5,3f9.4,2x,3f9.4,2x,f9.4)
                endif
              enddo
            enddo

C           Compute weights; store for now in wk
            allocate(wk(nfinal))
            dy = 0
            do  ifinal = 1, nfinal
              i = i+1
              call shorps(1,qlat,[72,2,2],qefbz(1:3,ipqbnd(ifinal))-qpi(:),q)
              wk(ifinal) = 1 / ((kfbrd/2)**2 + ddot(3,q,1,q,1))
              dy = dy + wk(ifinal)
            enddo
            call dscal(nfinal,1/dy,wk,1)

C           Assemble list equivalent points on irreducible wedge
            allocate(iwk(nfinal))
            do  ifinal = 1, nfinal
              iq = ipqbnd(ifinal)
              i1 = mod(iq-1,n1e)+1
              i2 = mod((iq-1)/n1e,n2e)+1
              i3 = mod((iq-1)/(n1e*n2e),n3e)+1
              i  = ipqe(i1,i2,i3) ! Corresponding irr point
              iwk(ifinal) = i     ! This is irr point
C              print 853, ifinal,i1,i2,i3,iq,i,qefbz(:,iq),wk(ifinal)
C  853         format(6i5,3f9.4,2x,f12.6)
            enddo

C           Contract list in full BZ to unique points in irr wedge
            allocate(wtqf(nfinal),iprm(nfinal))
            iq = iwk(1); ipqbnd(1) = iq ! First point is center point
            call dpzero(wtqf,nfinal); wtqf(1) = wk(1) ! Weight from center point
            call ivheap(1,nfinal-1,iwk(2),iprm(2),1) ! Sort remaining points
            forall (i = 2:nfinal) iprm(i) = iprm(i) + 1
            k = 1  ! Counter for current irr qp (1st point reserved)
            i = iwk(iprm(nfinal)) ! Will become counter for prior irr qp
            do  ifinal = 2, nfinal
              if (iwk(iprm(ifinal)) == iq) then ! Duplication of central point
                wtqf(1) = wtqf(1) + wk(ifinal)
                cycle
              endif
              if (i /= iwk(iprm(ifinal))) k = k+1 ! A new entry in the list
              i = iwk(iprm(ifinal))
              ipqbnd(k) = i; wtqf(k) = wtqf(k) + wk(ifinal)
            enddo
            deallocate(iwk)

            call info8(10,1,0,' PES(3d) at q=%s,(%3;6d)  KE+V0=%;3d  kT=%;3geV'//
     .        '  kf broadening %;3d(%;3dA^-1) over %i qp (%i irr)',
     .        qpi,Ek0,TeV,kfbrd,kfbrd*(2*pi)/(alat*a0),nfinal,k,8)

            nfinal = k - 1 ! nfinal = points other than given qp (consistent with lqpe def)

C            call info5(10,1,0,' PES:  no KE correction of q.'//
C     .        '  kf:  broadening=%;3dA^-1  range=+/-%;2gA^-1 (%i pts)'//
C     .        '  kT=%;3geV',
C     .        kfbrd*(2*pi)/(alat*a0),
C     .        dlength(3,nfinal/2*dqperp*knorm(:),1),nfinal,TeV,0)

          elseif (limesh == 0) then
            call info2(10,0,0,' Use as-given %sx%3i k-mesh'//
     .        '%?#n# ... include %-1j%i extra q not on mesh##',nkeabc,nqbnd)
          else if (limesh == 1 .or. limesh == 4) then
            call info2(10,0,0,' Interpolate to %sx%3i k mesh'//
     .        '%?#n# ... include %-1j%i extra q not on mesh##',nkeabc,nqbnd)
          else
            call info2(10,0,0,' Use as-given %sx%3i k-mesh'//
     .        '%?#n# ... %-1j%i extra qp also on the mesh##',nkeabc,nqbnd)
          endif
          deallocate(ipqe)  ! It will be remade, this time with nqe points

          if (nqbnd > 0 .and. .not. lgetev) then
            print *, 'you must use getev in this mode, sorry'
            goto 10
          endif

C         QP levels for the nqe points
          allocate(eige(nbandx,nqe,nsp),sige(nbandx,nqe,nsp),qeibz(3,nqe),weibz(nqe))
          allocate(ipqe(nqefbz,1,1))
          i = 71                 ! Make qeibz(:,nqbnd+1:), weibz(nqbnd+1:), ipqe(nqbnd+1:,1,1))
          if (lgetev) i = i+1000 ! Make eige also
          if (linputef) i = i+2  ! use given ef0
          i = i+4000             ! Make spin weights if noncollinear
          call evalqp(i,s_bz,s_ctrl,s_ham,s_lat,s_mix,s_move,s_pot,
     .      s_str,s_spec,s_site,s_optic,s_gw,s_dmft,s_strn,nbandx,ioff,nsp,nqe,
     .      nkeabc,nqe-nqbnd,nqbnd,lunits,ef0,eige,qeibz,weibz,ipqe,evshft)

C         Fill in ipqe, possibly eige for first nqbnd points
          if (nqbnd > 0) then
            call dcopy(3*nqbnd,qpbnd,1,qeibz,1)
            call dcopy(3*nqbnd,qpbnd,1,qefbz,1)
            i = 50               ! Make qeibz(:,1:nqbnd), ipqe(1:nqbnd,1,1))
                                 ! qeibz is not made, merely written back to itself
            if (lgetev) i = 1050 ! eige also
            i = i+4000           ! spin weights if noncollinear
            weibz(1:nqbnd) = 0   ! No weights --- offset points not part of a mesh
            call evalqp(i,s_bz,s_ctrl,s_ham,s_lat,s_mix,s_move,
     .        s_pot,s_str,s_spec,s_site,s_optic,s_gw,s_dmft,s_strn,nbandx,ioff,
     .        nsp,nqe,nkeabc,nqbnd,0,lunits,ef0,eige,qeibz,xx,ipqe,evshft)
          endif

C   ... irr < 0 : calculations limited to supplied (irregular) k mesh
        elseif (irr < 0) then
          if (lgetev) call info0(20,0,0,' Note: --getev is ignored in irr mesh for now')

          nqeibz = nqbnd
          nqe    = nqbnd
          nqefbz = nqbnd
          qeibz => qpbnd
!         nullify(qefbz,weibz,ipqe)
          deallocate(qefbz,weibz); allocate(qefbz(1,1),weibz(1))
          allocate(ipqbnd(nqbnd))
          do  iq = 1, nqbnd
            ipqbnd(iq) = 0
            do  i = 1, nqf
              q = qpbnd(1:3,iq)-qpf(1:3,i)
              if (latvec(1,tolq,plat,q)) then
                ipqbnd(iq) = i
                exit
              endif
            enddo
            if (ipqbnd(iq) == 0) call rx1('qp = %s,(%3;5,5d) not on supplied irregular mesh',qpbnd(1,iq))
          enddo

C       irr=0.  Supposed to be when qp do not match a regular mesh.  Which case does this handle?
        else
          if (irr == 0) call rx('sugws: check irr=0 case')
          nqeibz = nqibz
          nqe    = nqibz
          nqefbz = nqfbz
          qeibz => qibz
          qefbz => qfbz
          weibz => wibz
          ipqe => ipq
C         sigewq => sigwf
          allocate(eige(nbandx,nqe,nsp))
        endif

C   ... Make sigm, evals at interpolated k-points.  This block does the following:
C       Subtracts vxc from sigm; interpolate sigm-vxc to list of (new irr mesh + nqbnd points)
C       eigi from eigf (eigf=evals at qp given in se file)
C       sigi from sigvf (sigvf=QSGW potential at qp given in se file)
C       sigewq = sigwf-sigi (dynamical self-energy - static potential at qp given in se file)
C       If (lgetev), eigenvalues are be calculated at the interpolated points.
C       Used for exact interpolation of eigi to improve interpolation for sigi.
        ispx = isp; if (nspc == 2) ispx = 1
   18   continue
        if (irr >= 0) then
        if (ispx == isp) then
          if (associated(sigewq)) deallocate(sigewq)
          allocate(sigewq(nomg,nblst,nqe,nsp),eigkb(nblst),cwt(nblst))
          eigkb = 0
          allocate(sigi(nblst,nqe,nsp),eigi(nblst,nqe,nsp),eigiq0(nblst))
          allocate(sexi(nblst,nqe,nsp),vxci(nblst,nqe,nsp))
        endif
        call info8(10,1,0,' SUGWS:  interpolate self-energies, '//
     .    '%i band(s), %i%?#n#%-1j+%i## qp%?#(n==2)#, spin %i#%j#'//
     .    '%?#n#, independent eQP calc##',
     .    nblst,nqe-nqbnd,nqbnd,nsp,ispx,isw(lgetev),7,8)

        ssdy = 0
        do  iq = 1, nqe
          call shorps(1,qlat,[72,2,2],qeibz(1:3,iq),q)
          if (sum(abs(qeibz(1:3,iq)-q)) > tolq) then
            call info2(41,0,0,'%1fq translated from q =%3;12,6D '//
     .        'to q =%3:1,6;6d',qeibz(1:3,iq),q)
          endif
C          call info5(10,0,0,'%1fSE spectra for q=%3;12,6D  ib=%n:1i',
C     .      q,nblst,iblst,0,0)
C         q-interpolation from data read from se file (nqfbz,ipq,eigf,sigvf,sigwf).
C         To identify neighboring points the full BZ must be used.
C         But evals are given for irreducible points only.
C         ipqe maps point in full BZ to irr point.
          i = 7                 ! interpolate eigi<-eigf,sigi<-sigvf,sigewq<-sigwf-sigi
          if (lpeqp) i = 1      ! interpolate eigi<-eigf
          if (lhf) i = i+16     ! interpolate sexi<-sexf
          if (llxc) i = i+32    ! interpolate vxci<-vxcf
C         eigi = 0; print *, '!!'
          if (seqitp(i,q,plat,nqfbz,qfbz,ipq,nomg,nbandx,nqf,ispx,iblst,nblst,eigf,sigvf,sexf,vxcf,
     .      sigwf,sigi(1:nblst,iq,ispx),eigi(1:nblst,iq,ispx),sigewq(1:nomg,1:nblst,iq,ispx),
     .      sexi(1:nblst,iq,ispx),vxci(1:nblst,iq,ispx)) < 0) goto 10
          if (lpeqp .and. lso == 1) then ! the second spin
            if (seqitp(1,q,plat,nqfbz,qfbz,ipq,nomg,nbandx,nqf,ispx,iblst,nblst,
     .        eigf(1,1,2),sigvf,sexf(1,1,2),vxcf(1,1,2),sigwf,sigi(1:nblst,iq,2),eigi(1:nblst,iq,2),
     .        sigewq,sexi(1:nblst,iq,2),vxci(1:nblst,iq,2)) < 0) goto 10
          endif
C         Improve eigenvalue interpolation using eige, which is made by evalqp
C         eige and eigf may differ because of a different Fermi energy or
C         they may originate from a different QSGW potential.
C         Improvement calculable for evals only, but ...
C         assume that shift in sigi tracks shift in eigi.
C         No shift for sigwf because sigi was already subtracted.
C         Note: sigi is shifted only if QSGW sigma is available
C         Maybe should consider shifting sigwf in that case.
          if (lgetev) then
            call dcopy(nblst,eigi(1,iq,ispx),1,eigiq0,1)
            sdy=0; sdy2=0
            do  kb = 1, nblst
              ib = iblst(kb)
              dy = eige(ib,iq,ispx) - eigi(kb,iq,ispx)
              sdy = sdy + dy/nblst    ! Average
              sdy2 = sdy2 + dy*dy/nblst ! for root mean square
              eigi(kb,iq,ispx) = eigi(kb,iq,ispx) + dy
              if (lqsgws) sigi(kb,iq,ispx) = sigi(kb,iq,ispx) + dy
            enddo
            ssdy = ssdy + sdy
            sdy2 = sdy2-sdy**2+1d-12
C           call pshpr(45)
            call info5(45,0,0,' eval adjusted iq=%i:%N old %n;12,6D',iq,nblst,eigiq0,0,0)
            call info2(45,0,0,' new %n;12,6D',nblst,eigi(1,iq,ispx))
            call info5(45,0,0,' shf %n;12,6D  mean %;4d  rms %;3g',
     .        nblst,eigi(1:nblst,iq,ispx)-eigiq0,sdy,sqrt(sdy2),0)
C           call poppr
          endif
        enddo
        call info5(31,0,0,' mean eval shift (%i qp, %i bands) = %;6,6d',nqe,nblst,ssdy/nqe,4,5)
C       Shift sigvf or sigwf to enforce QSGW condition, (sig-vxc)_w=eQP=0
C       If: (1) evsync had been invoked previously and (2) Fermi level unchanged,
C       should not affect states on original mesh.  It will affect states on interpolated mesh
        if (lgetev .and. lqsgws .and. .not. lpeqp) then
          i = 11
          if (iprint() >= 41) i = 12
C         Bug fix 18 Nov 2016.  nqf -> nqe
          call schkpole(i,nblst,nqe,1,nomg,ommin,ommax,eigi(1:nblst,1:nqe,ispx),
     &                                   sigewq(1:nomg,1:nblst,1:nqe,ispx),shft) ! passing only eigi(1,1,ispx) is formally not correct and gfortran cannot be silenced about it, eigi(1:nblst,1:nqe,ispx) is contiguous and will not produce temporaries.
C debug check; shouldn't be needed
C          call schkpole(01,nblst,nqe,1,nomg,ommin,ommax,eigi(1,1,ispx),
C     .      sigewq(1,1,1,ispx),shft)

        endif

C       irr < 0  requires iblst, ipqbnd
        else ! irr < 0
          allocate(sigewq(nomg,nblst,nqe,nsp),eigkb(nblst),cwt(nblst))
          allocate(sigi(nblst,nqe,nsp),eigi(nblst,nqe,nsp),eigiq0(nblst))
          allocate(sexi(nblst,nqe,nsp),vxci(nblst,nqe,nsp))
          eigkb = 0
          do  isp = 1, nspse
          do  iq = 1, nqe
            do  kb = 1, nblst
              eigi(kb,iq,isp) = eigf(iblst(kb),ipqbnd(iq),isp)
              if (.not. lpeqp) then
                sigi(kb,iq,isp) = sigvf(iblst(kb),ipqbnd(iq),isp)
                forall (i=1:nomg) sigewq(i,kb,iq,isp) = sigwf(i,iblst(kb),ipqbnd(iq),isp)
              endif
              if (lhf) sexi(kb,iq,isp) = sexf(iblst(kb),ipqbnd(iq),isp)
              if (llxc) vxci(kb,iq,isp) = vxcf(iblst(kb),ipqbnd(iq),isp)
            enddo
          enddo
          enddo
        endif ! irr blocks

        if (loptics>0 .or. ljdos>0) then
          if (ispx == 1 .and. nsp > nspc) then
            ispx = 2
            goto 18
          else
            ispx = 1
          endif
        endif

C   --- Setup for interpolated frequency mesh ---
C       omgn = given  mesh; omgx = interpolated mesh
        nomgx = (nomg-1)*multi + 1
        if (allocated(omgn)) deallocate(omgn,omgx)
        allocate(omgn(nomg),omgx(nomgx))
        call sfitrp(1,nomg,nomgx,ommin,ommax,multi,xx,xx,eps,
     .    omgn,omgx,xx,xx,xx,xx,xx,xx,jx1,jx2)
        iw1 = 1
        iw2 = nomgx
        if (range(1) /= NULLI) then
          iw2 = 1
          do  iwx = 1, nomgx
            if (omgx(iwx) < range(1)) cycle
            if (omgx(iwx) > range(2)) cycle
            if (iw1 == 1) iw1 = iwx
            iw2 = max(iw2,iwx)
          enddo
        endif

        if (lsem) then
          i = nblst ; if (ldos) i = 1
          allocate(omgmats(nmatsubara),seimats(nmatsubara,i))
          call info2(10,0,0,'%1fSelf-energy on Matsubara axis, beta=%d, nmats=%i',
     .      beta,nmatsubara)
          forall (i=1:nmatsubara) omgmats(i) = dcmplx(0d0,pi*(2*i-1)/beta)
          allocate(iprm(1))
        endif

C   ... At this point all spectral information is folded into eigi,sigewq
        if (associated(seitpw)) deallocate(seitpw)
        if (ldos) goto 14

C   --- Spectral function, PE spectrum, or self-energy at one k-point qpi ---
        if (allocated(sumag)) deallocate(sumag,sumag0,ag,ag0,sag,sag0)
C                  int A(ib)   int A0(ib)     A(ib)     A0(ib)
        allocate(sumag(nomgx),sumag0(nomgx),ag(nomgx),ag0(nomgx))
C                sum_ib A   sum_ib A0     dA/dw     sigma-vxc   interp sig
        allocate(sag(nomgx),sag0(nomgx),dag(nomgx),seitp(nomg),seitpw(nomgx))
        if (lpe3d .or. lpe .or. lpeqp) then
          if (.not. allocated(pes)) allocate(pes(nomgx),pesi(nomgx))
          call dpzero(pes,nomgx)
          if (.not. lpe3d) ifi = fopna('ebqp',-1,0) ! Dump bands along qperp to file
          if (lresib) call rx('For PE simulation resolved by band, do one band at a time. Sorry')
        endif
        if (lresib) then
          allocate(sagkb(nomgx,nblst),sag0kb(nomgx,nblst),seitpwkb(nomgx,nblst))
          call dpzero(sagkb,nomgx*nblst); call dpzero(sag0kb,nomgx*nblst)
        endif

C   --- Loop over final states (PE) ... just 1 final state if SE only
        if (lpe3d.or.lpe.or.lpeqp) then
          call info5(10,1,0,' PES at q=%s,(%3;6d)  ib=%ni  broadened over 1+%i irr qp',
     .      qpi,nblst,iblst,nfinal,5)
          if (lpe3d) qsave = qpi
        endif

        call dpzero(sag,nomgx); call dpzero(sag0,nomgx)
        do  ifinal = 1, nfinal+1

        if (limesh == 2 .or. limesh == 4) then
          qpi(:) = qeibz(:,ipqbnd(ifinal))
        else
          qpi = qeibz(:,ifinal)
        endif

        call shorps(1,qlat,[72,2,2],qpi-qtarg,q); fac = dlength(3,q,1) ! Length

C       Shorten qpi
        call shorps(1,qlat,[72,2,2],qpi,q)
        if (sum(abs(qpi-q)) > tolq) then
          call info2(40,0,0,'%1fq translated from %s,q=(%3;6d) to q=(%3;6d)',qpi,q)
          qpi = q
        endif

        call info8(10,0,0,' Spectra for q=%3:1,6;6d  ib=%n:1i%?#n>1#  d=%;5F  weight=%;5F##',
     .    qpi,nblst,iblst,nfinal,fac,wtqf(ifinal),7,8)
        call info0(10,0,0,
     .    '  ib Npeak   eQP       w(Amax)     shift%6fZ%7fFWHM     w(A2)%7fshift   A2/Amax     w*    rank')

        if (lpe3d .or. lpe .or. lpeqp) call dpzero(pesi,nomgx)

        do  kb = 1, nblst

C     ... Interpolate sigma to qpi ... point should already be in nqefbz
C         call pshpr(45)
C         Extract from given data at qpi.
C         eigitp and sigewq were mapped from original file data
C         See branch 'irr > 0' above. eigi, sigewq on regular mesh, kb mapped; not iq,isp
          if (limesh == 2 .or. limesh == 4) then
C           eigitp = eigf(iblst(kb),ipqbnd(iq),ispx)
            eigitp = eigi(kb,ipqbnd(ifinal),ispx)
            forall (i=1:nomg) seitp(i) = sigewq(i,kb,ipqbnd(ifinal),ispx)

C         See branch 'irr < 0' above. kb,iq,ispx all mapped
          elseif (limesh == 3) then
            eigitp = eigi(kb,ifinal,1)
            forall (i=1:nomg) seitp(i) = sigewq(i,kb,ifinal,1)

C         Interpolate to desired qp.  eigi, sigewq range over bands in nblst
          else
            i = 9 ; if (lpeqp) i = 1
            if (.not. associated(qefbz)) call rx('mode not accessible ... no BZ mesh defined')
            if (seqitp(i,qpi,plat,nqefbz,qefbz,ipqe,nomg,nblst,nqe,
     .      ispx,kb,1,eigi,xx,xx,xx,sigewq,xx,eigitp,seitp,xx,xx) < 0)
     .      goto 10
          endif
          eigkb(kb) = eigitp

          cwitp = 1; cwt(kb) = 1
          if (lpeqp .and. lso == 1) then
            if (limesh == 2 .or. limesh == 4 .or. limesh == 3) call rx('not ready for this branch')
            if (seqitp(1,qpi,plat,nqefbz,qefbz,ipqe,nomg,nblst,nqe,
     .        ispx,kb,1,eigi(1,1,2),xx,xx,xx,sigewq,xx,cwitp,seitp,xx,xx) < 0)
     .      goto 10
            if (isp == 2) cwitp = 1-cwitp
            cwt(kb) = cwitp
          endif

          if (kb == nblst .and. .not. lpe3d) then
            ifi = fopna('ebqp',-1,0)
            write(ifi,"(i6,3f12.6,2x,20f12.7)") ifinal, qpi, eigkb
          endif

C     ... sigma(qpi) on Matsubara frequencies
          if (lsem) then
            call susepade(0,nomg,omgn,-2d0,2d0,-10d0,10d0,npade,[0])
            deallocate(iprm); allocate(iprm(npade))
            call susepade(1,nomg,omgn,-2d0,2d0,-10d0,10d0,npade,iprm)
            call sig2mats(0,nomg,omgn,eigitp,seitp,npade,iprm,nmatsubara,omgmats,seimats(1,kb))
            cycle
          endif

C     ... Interpolate sigma(qpi) to finer energy mesh; generate ag
          i = 2 ; if (lpeqp) i = 4
          call sfitrp(i,nomg,nomgx,ommin,ommax,multi,eigitp,seitp,
     .      eps,omgn,omgx,seitpw,sumag,sumag0,ag,ag0,Zfac,jx1,jx2)

          call info5(50,0,0,' ib %,4i int A(w)dw = %d  int A0(w)dw = %d'
     .      //'%?#n#  spin projection %d##',iblst(kb),
     .      sumag(nomgx),sumag0(nomgx),isw(cwt(kb) /= 1),cwt(kb))

C          call prrmsh('A[G]',omgx,sumag,nomgx,nomgx,1)
C          call prrmsh('A[G0]',omgx,sumag0,nomgx,nomgx,1)
C          call prrmsh('A[G]',omgx,ag,nomgx,nomgx,1)
C          call prrmsh('A[G0]',omgx,ag0,nomgx,nomgx,1)

          if (lpeqp) ag = ag0

C   ...   Find and print to stdout all peaks in A
          call spec_peaks(iw1,iw2,nomgx,omgx,ifi,iblst(kb),shft,eigitp,ag,dag,Zfac,ommaxs)

C     ... Add spectral function from this QP to spectral function
          if (ifinal == 1) then
            eigitp_save = eigitp
            cwt_save = cwitp
            call daxpy(nomgx,cwitp,ag,1,sag,1)
            call daxpy(nomgx,cwitp,ag0,1,sag0,1)
            if (lresib) then
              call daxpy(nomgx,cwitp,ag,1,sagkb(1,kb),1)
              call daxpy(nomgx,cwitp,ag0,1,sag0kb(1,kb),1)
            endif
          endif

C     ... Add to PE spectra
C         if (lpe3d .or. lpe .or. lpeqp) then
          if (lpe .or. lpeqp) then
C           PE for this k: contribution weighted PE by Fermi function
            call daxpy(nomgx,cwitp/(1+exp(eigitp/TeV)),ag,1,pesi,1)
          endif
C         A(w) should be weighted by Ef for each energy instead of 1 weight for whole QP
          if (lpe3d) then
C           PE for this k: contribution weighted PE by Fermi function
            call daxpy(nomgx,cwitp/(1+exp(eigitp/TeV)),ag,1,pesi,1)
            do  i = 1, nomgx
              fac = cwitp/(1+dexp(min(omgx(i)/TeV,100d0)))
              pesi(i) = pesi(i) + fac*ag(i)
            enddo
          endif

        enddo ! loop over QP levels

        if (lpe3d .or. lpe .or. lpeqp) then
          call daxpy(nomgx,wtqf(ifinal),pesi,1,pes,1)
          havepe = 1

C         Debugging
C          allocate(wk3(nomgx,4,1))
C          forall (i=1:nomgx) wk3(i,1,1) = omgx(i)
C          forall (i=1:nomgx) wk3(i,2,1) = pes(i)
C          forall (i=1:nomgx) wk3(i,3,1) = pesi(i)
C          forall (i=1:nomgx) wk3(i,4,1) = ag(i)
C          call prmx('pes',wk3,nomgx,nomgx,4)
C          deallocate(wk3)

        endif

        enddo ! loop over final states

C   ... Printout peaks in cumulative spectral function and PES
        if (lpe3d .or. lpe .or. lpeqp) then
          deallocate(pesi)
C   ...   Find and print to stdout all peaks in Spectral function, PES
          shft = eigitp_save
          call info0(10,1,0,' Search for peaks in cumulative spectral function')
          call spec_peaks(iw1,iw2,nomgx,omgx,0,1,shft,eigitp_save,sag,dag,NULLR,ommaxs)
          call info0(10,1,0,' Search for peaks in PES')
          call spec_peaks(iw1,iw2,nomgx,omgx,0,1,shft,eigitp_save,pes,dag,NULLR,ommaxp)
          call info5(10,1,0,' qp=%3:1;6,6d  w(max A)=%;6,6d   w(max P)=%;6,6d',
     .      qsave,ommaxs,ommaxp,0,0)
        else

        endif
        deallocate(dag,seitp)
        havesi = 1
        if (lsem) havesi = 11
        lnsave = .true.
        goto 10

C   --- Make spectrum DOS or optics or k-resolved self-energies and/or spectral functions ---
C       This branch handles these modes:
C       ljdos>0                        q->0 optics, no vertex or LF
C       loptics>0                      q->0 optics, no vertex or LF
C       lse=.false. and lsem=.false. : integrated DOS.  Not written to disk.  No band mode
C       lse=.true.                     k-resolved spectral function.  Written to disk and program exits.
C                                      If no band mode, also write self-energy seq2
C       lsem=.true.                    k-resolved spectral function on Matsubara frequencies.  No band mode

   14   continue                ! Entry point when ldos=.true.
        if ((loptics>0 .or. ljdos>0) .and. (lband .or. lse .or. lsem))
     .    call rx('sugws: switches incompatible with optics')

C       Generate refined k mesh for DOS, weights
        nw12 = iw2-iw1+1        ! number of points within (iw1,iw2) window

        if (allocated(dos)) deallocate(dos)
        allocate(seitpw(nomgx),dos(nomgx,2))
        call dvset(dos,1,size(dos),dble(NULLI)) ! For debugging

        nqp = nqe               ! number of qp to loop over
        if (lband) nqp = nqbnd  ! number of qp in band mode

        if (lse) then ! write on the fly
          call info8(10,0,0,' Spectral function%?#(n==0)# and self-energy## for '//
     .      '%s,%i qp  energy window (%;4d,%;4d), %i points (spacing=%;4d)  eps=%;4d',
     .      isw(lband),nqp,omgx(iw1),omgx(iw2),nw12,omgx(2)-omgx(1),eps,8)

          fn = 'spq'; if (isp == 2) fn = 'spq2';
          lascii = 20
          if (lascii /= 0) then
            ifi = fopna(trim(fn),-1,0)
            if (lupdateqp) then
              jfi = ifi
              ifi = fopna('sqpx',-1,0)
            endif
          else
            ifi = fopna(trim(fn),-1,4)
            if (lupdateqp) then
              jfi = ifi
              ifi = fopna('sqpx',-1,4)
            endif
          endif
          i = lunits*100 + lascii + 4
          j = 112; if (lhf) j = j+20  ! Write spectral function later
          xx(1) = ef0; if (ef0 /= NULLI .and. lunits == 1) xx(1) = ef0*rytoev
          call ioseh(i,j,1,1,nblst,nqp,nw12,omgx(iw1),omgx(iw2),xx,nblst,iblst,nkpan,-ifi)
          call iose2(lascii,j,j,1,nblst,nqp,nw12,qeibz,eigi(1,1,ispx),sigi(1,1,ispx),xx,sexi(1,1,ispx),vxci(1,1,ispx),-ifi)
          if (.not. lband) then ! Also write header for the self-energy
            fn2 = 'seq'; if (isp == 2) fn2 = 'seq2'
            jfi = fopna(trim(fn2),-1,0)
            j = 12; if (lhf) j = j+20
            call ioseh(i,j,1,1,nblst,nqp,nw12,omgx(iw1),omgx(iw2),chempot,nblst,iblst,nkpan,-jfi)
            call iose2(lascii,j,j,1,nblst,nqp,nw12,qeibz,eigi(1,1,ispx),sigi(1,1,ispx),xx,sexi(1,1,ispx),vxci(1,1,ispx),-jfi)
            allocate(seb(nw12,nblst))  ! Self-energy for each band
          endif
        elseif (lsem) then ! write on the fly
          fn = 'seim'; if (isp == 2) fn = 'seim2'
          ifi = fopna(trim(fn),-1,4)
          rewind ifi
          if (nomgx == nomg) deallocate(seitpw) ! Point to siqewq instead
          call info5(10,0,0,' Self-energy on Matsubara axis for %i qp, beta=%d, nmats=%i',
     .      nqp,beta,nmatsubara,0,0)
          write(ifi) 0, 1, nkeabc, beta, nmatsubara, eps
        elseif (loptics>0 .or. ljdos>0) then      ! optics
C         Find point closest to omega=0
          iw0 = 0
          call huntx(omgx(iw1),nw12,0d0,[0],iw0)
          if (iw0 == nw12) call rx('band energy window for optics')
          if (abs(omgx(iw1-1+iw0)) < abs(omgx(iw1-1+iw0+1))) iw0 = iw0-1
          if (kbT == NULLI .and. lunits == 1) kbT = .0259d0
          if (kbT == NULLI .and. lunits == 0) kbT = .0259d0/rytoev
          call info8(10,1,0,' Optics on  %s,(%3i) qp mesh'//
     .      '%N energy window (%2;4d) %i points (spacing=%;4d)  omg(%i)=%d  kT=%;4d  eps=%;4d',
     .      nkeabc,[omgx(iw1),omgx(iw2)],nw12,omgx(2)-omgx(1),iw0+1,omgx(iw1+iw0),kbT,eps)
C         Fermi function
          allocate(fermifn(nw12))
          do  i = 1, nw12
            call delstp(-1,omgx(iw1-1+i)/kbT,xx,fermifn(i),xx)
          enddo
        else
          call info8(10,0,0,' DOS on  '//'%s,(%3i) qp mesh'//
     .      '  energy window (%;4d,%;4d), %i points (spacing=%;4d)  eps=%;4d',
     .      nkeabc,omgx(iw1),omgx(iw2),nw12,omgx(2)-omgx(1),eps,7,8)
        endif

        allocate(agb(nw12,nblst))  ! Spectral function for each band
        if (allocated(sumag)) deallocate(sumag,sumag0,ag,ag0)
        allocate(sumag(nomgx),sumag0(nomgx),ag(nomgx),ag0(nomgx))

        if (.not. lse .and. .not. lsem) then
          call dpzero(dos,nomgx*2)
          call dscal(nqe,0.5d0,weibz,1) ! so that sum of weibz is 1
        endif
        if (ljdos>0) then
          deallocate(dos)
          ndos = nw12-iw0
          allocate(dos(ndos,0:nsp))
C       Read and synchronize velocity matrix elements
        elseif (loptics>0) then
          deallocate(dos)
          ndos = nw12-iw0
          allocate(dos(ndos,0:3*nsp)); call dpzero(dos,size(dos))
          call info0(10,1,0,' reading optics data from file optdatac ...')
          ifi = fopna('optdatac',-1,4); rewind ifi
          read(ifi) xx(1), i,j,k   ! Fermi level, ndham, nspx, nkp
          ltmp = isanrg(i,ndham,ndham,'sufuned:','optdatac file''s ndham',.true.)
          ltmp = isanrg(j,nsp,nsp,'sufuned:','optdatac file''s nsp',.true.)
          ltmp = isanrg(k,nqe,nqe,'sufuned:','optdatac file''s nkp',.true.)
          allocate(wk3(ndham,nsp,nqe))
          call dpdump(wk3,ndham*nsp*nqe,ifi)
          forall (i=1:ndham, j=1:nsp, k=1:nqe) wk3(i,j,k) = wk3(i,j,k) - xx(1)
          if (lunits == 1) call dscal(ndham*nsp*nqe,rytoev,wk3,1)
          read(ifi) ix(1), ix(2),ix(3),ix(4),ix(5),ix(6) !, qp ! nfilo,nfiup,nemlo,nemup,nspx,nkp,s_bz%qp
          i = minval(iblst)
          if (ix(1)>i .or. ix(3)>i)
     .      call rx2('sfuned: velocity matrix elements start at band %i but expected band %i',
     .      max(ix(1),ix(3)), i)
          i = maxval(iblst)
          if (ix(2)<i .or. ix(4)<i)
     .      call rx2('sfuned: velocity matrix elements extend to band %i but expected band %i',
     .      min(ix(1),ix(3)), i)
          if (ix(5) /= nsp .or. ix(6) /= nqe) call rx('mismatch nsp or nkp')
C          allocate(wk5(3,ix(1):ix(2),ix(3):ix(4),nsp,nqe))
C          allocate(optmt(3,1:nblst,1:nblst,nsp,nqe))
          allocate(optmt(3,ix(1):ix(2),ix(3):ix(4),nsp,nqe))
          call dpdump(optmt,3*(ix(2)-ix(1)+1)*(ix(4)-ix(3)+1)*nsp*nqe,ifi)
          sdy = 0 ; sdy2 = 0; k = 0
          do  isp = 1, nsp
            do  iq = 1, nqe
              do  kb = 1, nblst
                j = j+1
                ib = iblstf(kb)
                dy = wk3(ib,isp,iq) - eigi(kb,iq,isp)
                sdy = sdy + dy      ! for average
                sdy2 = sdy2 + dy*dy ! for root mean square
              enddo

CC             Copy velocity matrix elements for subbands
C              do  ib = 1, nblst
C              do  kb = 1, nblst
C                i = iblst(ib); k = iblst(kb)
C                optmt(1:3,ib,kb,isp,iq) = optmt(1:3,i,k,isp,iq)
C              enddo
C              enddo

            enddo
          enddo
          isp = 1
          sdy = sdy/j; sdy2 = sdy2/j - sdy**2 + 1d-12
          call info5(20,0,0,' optdatac eigenvalues shifted by mean %,4;4d  rms %;3g (%i pts)',
     .      sdy,sqrt(sdy2),j,4,5)
        endif
        if (ljdos>0 .or. loptics>0) then
          forall (i=1:ndos) dos(i,0) = omgx(iw0+iw1+i-1) ! First column = frequency
        endif

C   ... Re-entry point spin 2 loptics or jdos
        if (lupdateqp) allocate(dag(nomgx))
        if (nomgx == 1) nofilter = .true.
   16   continue
        call pshpr(iprint()-30)
        do  iq = 1, nqp
!          print *, '!!'; if (iq /= 7) cycle
          do  kb = 1, nblst
            ib = iblst(kb)

C       ... Interpolate sigma to finer frequency mesh; generate spectral function
            eigitp = eigi(kb,iq,ispx)
            sigitp = sigi(kb,iq,ispx)
            if (lsem .and. nomgx == nomg) then
              seitpw => sigewq(:,kb,iq,ispx) ! Special Matsubara SE mode with no interpolation
            else
              i = 2; if (nofilter) i = i+100
              call sfitrp(i,nomg,nomgx,ommin,ommax,multi,eigitp,sigewq(1:nomg,kb,iq,ispx),
     .          eps,omgn,omgx,seitpw,sumag,sumag0,ag,ag0,Zfac,jx1,jx2)
              if (.not. nofilter .and. (jx1 /= 0 .or. jx2 /= 0)) then
                call info5(0,0,0,' (warning) rational function '//
     .            'interpolation failed,  i2=%i  ib=%i  kb=%i',iq,ib,kb,0,0)
              endif
            endif

C       ... Self-energy on Matsubara frequency for this iq,ib
            if (lsem) then
              call susepade(0,nomgx,omgx,-2d0,2d0,-10d0,10d0,npade,[0])
              deallocate(iprm); allocate(iprm(npade))
              call susepade(1,nomgx,omgx,-2d0,2d0,-10d0,10d0,npade,iprm)
              call sig2mats(0,nomgx,omgx,eigitp,seitpw,npade,iprm,nmatsubara,omgmats,seimats)
              write(ifi) iq,qeibz(1:3,iq),ib
              write(ifi) seimats
              cycle
            endif

            if (cmdopt('--table',7,0,strn)) call dcopy(nw12,ag(iw1),1,agb(1,kb),1)

            if (lupdateqp) then
              shft = eigitp
              call togpr
              call spec_peaks(iw1,iw2,nomgx,omgx,0,1,shft,eigitp,ag,dag,Zfac,ommaxp)
              call togpr
              eigi(kb,iq,ispx) = ommaxp(1)
            endif

C           Retain spectral function for this band; exit loop
            if (lse .or. loptics>0 .or. ljdos>0) then
              call dcopy(nw12,ag(iw1),1,agb(1,kb),1)
              if (loptics>1 .or. ljdos>1) then
                call dcopy(nw12,ag0(iw1),1,agb(1,kb),1)
              endif
C             Next line for debugging : find index near position of QP peak
C             call huntx(omgx(iw1),nomgx-iw1,eigitp,[0],iwx); print *, kb, iwx
              if (lse .and. .not. lband) then ! and self-energy for this band
                call dcopy(2*nw12,seitpw(iw1),1,seb(1,kb),1)
                seb(:,kb) = seb(:,kb) + sigitp
              endif
              cycle
            endif

C           Accumulate dos
C           print *, kb,iq,ag(10)
            do  iwx = 1, nomgx
              dos(iwx,1) = dos(iwx,1) +  weibz(iq)*2/dble(nsp) * ag(iwx)
              dos(iwx,2) = dos(iwx,2) +  weibz(iq)*2/dble(nsp) * ag0(iwx)
            enddo
          enddo                 ! loop over bands

          if (cmdopt('--table',7,0,strn)) then
            write(strn,"('debugging agb for iq=',i4,' q=',3f12.6)") iq, qeibz(1:3,iq)
            call prrmsh(trim(strn),omgx(iw1),agb,nw12,nw12,nblst)
          endif

C     ... RPA optics or joint DOS
          if (loptics>0 .or. ljdos>0) then

CC           Debugging test: compare cycling and noncyclic convolution, direct and FFT.
CC           Noncyclic FFT done by padding array
CC           For cyclic test, set npad = 0, uncomment noncyclic lines for wk(i) in brute force calculation
CC           For noncyclic test, set npad = 0, uncomment cyclic line for wk(i) in brute force calculation
C            npad = nw12                ! number of padding points.  Use npad=0 for cyclic
C            nwft = nw12+npad
C            allocate(zagb(nwft,nblst)) ! Map spectral function for FFT; shift to origin
C            allocate(convolve(nw12,nblst,nblst)) ! Convolutions of zagb
C            call fftz30(nwft,1,1,k1,k2,k3)
CC           Copy to complex, padded zagb. Do FFT's in advance so they are done only once
C            do  kb = 1, nblst
C              call suconvolve1D(0,nw12,iw0,nw12-iw0,npad,1,agb(1,kb),zagb(1,kb))
C              call fftz3(zagb(1,kb),nwft,1,1,k1,k2,k3,1,0,-1)
C            enddo
C
CC           Test: establish convolution by FFT matches brute force calculation
C            allocate(wk(nw12),zwk(nwft))
C            do  ib = 3, 5
C              do  kb = 4, 6
C
CC               Convolution by brute force
C                do  i = 0, nw12-1
C                  wk(1+i) = 0
C                  do  j = 0, nw12-1
CC                   Next two lines for noncyclic
C                    if (i+j >= nw12) cycle
C                    wk(1+i) = wk(1+i) + agb(1+j,ib)*agb(1+i+j,kb)/nw12
CC                   Next line for cyclic
CC                   wk(1+i) = wk(1+i) + agb(1+j,ib)*agb(1+cyclic(i+j,nw12),kb)/nw12
C                  enddo
C                enddo
C
CC               Convolution by FFT
C                call dcopy(2*nwft,zagb(1,ib),1,zwk,1)
C                call fftz3c(zwk,zagb(1,kb),nwft,1,1,k1,k2,k3,13,-1) ! data given as FT: invert zwk only on return
CC               Copy to real array, reversing order and taking last elements (corresponding to w>0)
C                j = 0
C                convolve(1,ib,kb) = zwk(1) * dble(nwft)/nw12  ! Point at origin
C                do  i = 1, nw12-1
C                  convolve(1+i,ib,kb) = zwk(nwft+1-i) * dble(nwft)/nw12
C                enddo
C                wk(:) = wk(:) - convolve(:,ib,kb)
C                xx(1) = dlength(nw12,convolve(1,ib,kb),1)
C                xx(2) = dlength(nw12,wk,1)
C                print *, ib,kb,xx
C              enddo
C            enddo
C            deallocate(zwk,wk)

            npad = nw12                ! number of padding points.  Use npad=0 for cyclic
            nwft = nw12+npad
            allocate(zagb(nwft,nblst,2)) ! Map spectral function for FFT; shift to origin
            allocate(convolve(nw12,nblst,nblst)) ! Convolutions of zagb
            allocate(wk(nw12),zwk(nwft),zwk2(nwft,0:2))
            allocate(zfermi(nwft)) ! shift and pad Fermi function
            call suconvolve1D(0,nw12,iw0,nw12-iw0,npad,1,fermifn,zfermi)
            call fftz30(nwft,1,1,k1,k2,k3)
C           Copy to complex, padded zagb. Do FFT's in advance so they are done only once
            do  ib = 1, nblst
              call suconvolve1D(0,nw12,iw0,nw12-iw0,npad,1,agb(1,ib),zagb(1,ib,1))
              forall (i=1:nwft) zagb(i,ib,2) = zfermi(i)*zagb(i,ib,1) ! fermifn * zagb
              call fftz3(zagb(1,ib,1),nwft,1,1,k1,k2,k3,1,0,-1) ! FFT of zagb
              call fftz3(zagb(1,ib,2),nwft,1,1,k1,k2,k3,1,0,-1) ! FFT of fermifn * zagb
            enddo

C           Convolution by FFT [f(w') A_ib(w')] A_kb(w+w') - A_ib(w') [f(w+w') A_kb(w+w')]
!                print *, '!!';
            jacobian = (omgx(2)-omgx(1)) * dble(nwft)
            do  ib = 1, nblst
              do  kb = 1, nblst
!                if (ib /= 32 .or. kb /= 33) cycle
                call dcopy(2*nwft,zagb(1,ib,1),1,zwk2(1,0),1)
                do  i = 1, nwft
                  zwk2(i,1) = zfermi(i)*zagb(i,ib,1)
                  zwk2(i,2) = zfermi(i)*zagb(i,kb,1)
C                 print *, i, zagb(i,ib,1)-zagb(i,kb,1)
                enddo
C                print *, ib,kb,sum(zagb(:,ib,1)-zagb(:,kb,1))
C                print *, dlength(nwft*2,zwk2(:,1)-zwk2(:,2),1)
C                print *, dlength(nwft*2,zagb(:,ib,1)-zagb(:,kb,1),1)

!                 call huntx(omgx(iw1),nomgx-iw1,eigi(ib,iq,ispx),[0],iwx); print *, ib, iwx; i = iwx
!                 call huntx(omgx(iw1),nomgx-iw1,eigi(kb,iq,ispx),[0],iwx); print *, kb, iwx
!                 print 642, ib,kb,i,iwx,iwx-i
!                 call prmx('agb(ib)',agb(1,ib),nw12,nw12,1)
!                 call prmx('agb(kb)',agb(1,kb),nw12,nw12,1)
C                 call zprm('zagb(ib)',1,zagb(1,ib,1),nwft,nwft,1)
C                 call zprm('zagb(kb)',1,zagb(1,kb,1),nwft,nwft,1)
C                 call zprm('Fermi * zagb(ib)',1,zwk2(1,1),nwft,nwft,1)
C                 call zprm('Fermi * zagb(kb)',1,zwk2(1,2),nwft,nwft,1)

C               Convolution [f(w') A_ib(w')] A_kb(w+w')
                call dcopy(2*nwft,zagb(1,ib,2),1,zwk2(1,1),1) ! FFT of fermifn * zagb(ib)
                call fftz3c(zwk2(1,1),zagb(1,kb,1),nwft,1,1,k1,k2,k3,13,-1) ! zwk2(1) -> convolution in w repsn
C               Convolution A_ib(w') [f(w+w') A_kb(w+w')]
                call dcopy(2*nwft,zagb(1,ib,1),1,zwk2(1,0),1) ! FFT of zagb(ib)
                call fftz3c(zwk2(1,0),zagb(1,kb,2),nwft,1,1,k1,k2,k3,13,-1) ! zwk2(0) -> convolution in w repsn
C               Copy to real array, reversing order and taking last elements (corresponding to w>0)
                j = 0
                convolve(1,ib,kb) = (zwk2(1,1)-zwk2(1,0)) * jacobian  ! Point at origin
                do  i = 1, nw12-1
                  convolve(1+i,ib,kb) = (zwk2(nwft+1-i,1)-zwk2(nwft+1-i,0)) * jacobian
                  if (convolve(1+i,ib,kb) < 0) then
                    call logwarn(2,'%26fwarning!  negative DOS')
                  endif
                enddo
!                call huntx(omgx(iw1),nomgx-iw1,eigi(ib,iq,ispx),[0],iwx); print *, ib, iwx; i = iwx
!                call huntx(omgx(iw1),nomgx-iw1,eigi(kb,iq,ispx),[0],iwx); print *, kb, iwx
!                print 642, ib,kb,i,iwx,iwx-i,omgx(iwx+iw1)-omgx(i+iw1)
!                print *, ib,kb; call zprm('convolve 1',1,zwk2(1,1),nwft,nwft,1)
!                print *, ib,kb; call zprm('convolve 2',1,zwk2(1,0),nwft,nwft,1)
!                print *, ib,kb; call prmx('convolve',convolve(1,ib,kb),nw12,nw12,1)
              enddo
            enddo

            call togpr
            if (iprint() >= 60) then
              call info5(10,0,0,' Sum of joint DOS decomposed by pairs, q=%s,(%3;4d)',qeibz(1,iq),2,3,4,5)
              call info2(10,0,0,' eqp=%n:;11,6D',nblst,eigi(iblst,iq,ispx))
              call info2(10,0,0,'%n:,11i',nblst,iblst)
              do  ib = 1, nblst
                write (stdo, 366) ib, (sum(convolve(:,ib,kb)), kb = 1, nblst)
  366           format(i5,100f12.6)
              enddo
C             call prmx('agb',agb,nw12,nw12,nblst)

C             deugging : pictures
C              do  ib = 1, nblst
C              do  kb = 1, nblst
C                if (sum(convolve(:,ib,kb)) < 0.1d0) cycle
C                call huntx(omgx(iw1),nomgx-iw1,eigi(ib,iq,ispx),[0],iwx); print *, ib, iwx; i = iwx
C                call huntx(omgx(iw1),nomgx-iw1,eigi(kb,iq,ispx),[0],iwx); print *, kb, iwx
C                print 642, ib,kb,i,iwx,iwx-i
C 642           format(' ib,kb,iw(ib),iw(kb)',4i5,' expect peak at',i4:' dw=',f12.6)
C                call prmx('convolve',convolve(1,ib,kb),nw12,nw12,1)
C              enddo
C              enddo

            endif
            if (mod(iq+1,100) == 1) call info2(30,0,0,' completed iq=%i of %i',iq,nqp)
            call togpr

C
            do  ib = 1, nblst
              do  kb = 1, nblst
!               print *, '!!!'; if (ib /= 32 .or. kb /= 33) cycle
!               call huntx(omgx(iw1),nomgx-iw1,eigi(ib,iq,ispx),[0],iwx); print *, ib, iwx; i = iwx
!               call huntx(omgx(iw1),nomgx-iw1,eigi(kb,iq,ispx),[0],iwx); print *, kb, iwx

C               Joint DOS
                if (ljdos>0) then
                  call daxpy(ndos,weibz(iq),convolve(1,ib,kb),1,dos(1,isp),1)
C               Joint DOS * |<v>|2
                else
                  j = 3*(isp-1)
                  do  i = 1, 3
                    fac = weibz(iq)*optmt(i,iblstf(ib),iblstf(kb),isp,iq)*2/dble(nsp)
                    call daxpy(ndos,fac,convolve(1,ib,kb),1,dos(1,i+j),1)
                  enddo
                endif
              enddo
            enddo
C            do  i = 1, nw12
C              if (dos(i,0) < .2d0) print *, dos(i,0), dos(i,1)
C            enddo
C            stop

            deallocate(wk,zwk,zwk2,zagb,convolve,zfermi)
          endif

          if (lse) then
            call dfdump(agb,nw12*nblst,-ifi)
            if (.not. lband) call dfdump(seb,2*nw12*nblst,-jfi)
          endif
        enddo                   ! Loop over qp
        call poppr

        if (lse .and. lupdateqp) then ! Remake file with updated header
          rewind jfi
          i = lunits*100 + lascii + 4
          j = 112; if (lhf) j = j+20 ! Write spectral function later
          xx(1) = ef0; if (ef0 /= NULLI .and. lunits == 1) xx(1) = ef0*rytoev
          call ioseh(i,j,1,nspse,nblst,nqp,nw12,omgx(iw1),omgx(iw2),xx,nblst,iblst,nkpan,-jfi)
          call iose2(lascii,j,j,1,nblst,nqp,nw12,qeibz,eigi(1,1,ispx),sigi(1,1,ispx),xx,sexi(1,1,ispx),vxci(1,1,ispx),-jfi)
          rewind ifi
          call ioseh(i-4+1,j,1,nspse,nblst,nqp,nw12,omgx(iw1),omgx(iw2),xx,nblst,iblst,nkpan,ifi)
          call iose2(lascii,j,j,1,nblst,nqp,nw12,qeibz,eigi(1,1,ispx),sigi(1,1,ispx),xx,sexi(1,1,ispx),vxci(1,1,ispx),ifi)
          do  iq = 1, nqp
            call dfdump(agb,nw12*nblst,ifi)
            call dfdump(agb,nw12*nblst,-jfi)
          enddo                 ! Loop over qp
          call dfclos(ifi)
        endif

        if (lse .or. lsem) then
          if (.not. lband) call info0(1,0,0,'%8fwrote self-energy to file '//trim(fn2))
          call rx0('wrote spectral function to file '//trim(fn))
        endif

        if (ljdos>0 .or. loptics>0) then
C         deallocate(fermifn,agb,sumag,sumag0,ag,ag0)
          if (isp == 1 .and. nsp == 2 .and. nspc == 1) then
            ispx = 2
            isp = 2
            call info0(20,0,0,' begin spin 2')
            goto 16
          endif
          if (loptics>0) then
C           Extra factors for Im eps
            fac = .001; if (lunits == 1) fac = fac*rytoev
            do  iwx = 1, ndos
              j = iwx
              if (dos(iwx,0) > fac) exit
            enddo
            fac = 32d0*pi**2/vol
            if (lunits == 1) fac = fac*rytoev**3
            k = 3*nsp
            forall (i=j:ndos) dos(i,1:k) = dos(i,1:k) * fac/dos(i,0)**2
            fn = 'opt'; if (loptics>1) fn = 'optni'
            strn = 'optics'
          else
            fn = 'jdos'; if (ljdos>1) fn = 'jdosni'
            j = 0
            k = nsp
            strn = 'joint DOS'
          endif
          ifi = fopna(trim(fn),-1,0)
          call ywrm(0,' ',1,ifi,'(1p15d15.7)',dos(1+j,0),0,ndos,ndos-j,1+k)
          call rx0('wrote '//trim(strn)//' to file '//trim(fn))
        endif

        if (.not. associated(qefbz,qfbz)) then
          deallocate(qefbz,qeibz,weibz,ipqe)
        endif

        havesi = 2
        lnsave = .true.
        goto 10

C ... merge
      elseif (outs(1:8) == 'merges12') then

        ifi = fopna('seq1',-1,1); jfi = fopna('seq2',-1,1); i = fopna('seq',-1,1);
        call mergesegw12(ifi,jfi,i)  ! Does not return

C ... show
      elseif (outs(1:5) == 'show ') then

C  --- Save interpolated SE ---
      elseif (outs(1:7) == 'savese ' .or.
     .        outs(1:8) == 'savesea ') then

        if (.not. lnsave) then
          call info0(10,0,0,' nothing to save, nothing written')
          goto 10
        endif

        lbin = outs(1:7) == 'savese '
        call word(outs,2,j1,j2)
        if (j2 >= j1) fn = outs(j1:j2)

C       File name for self-energy or DOS
        select case (havesi)
        case (1)
          if (j2 < j1) fn = 'seia'
          if (j2 < j1 .and. isp == 2) fn = 'seia2'
        case (2)
          if (j2 < j1) fn = 'sdos'
          if (j2 < j1 .and. isp == 2) fn = 'sdos2'
        case (11)
          if (j2 < j1) fn = 'seim'
          if (j2 < j1 .and. isp == 2) fn = 'seim2'
        case default
          call info0(2,0,0,' oops ... nothing to save')
          goto 10
        end select

C       File name for PES
        select case (havepe)
        case (1)
          fn = 'pes'
          if (lpeqp) fn = 'pesqp'
          if (isp == 2) then
            fn = 'pes2'
            if (lpeqp) fn = 'pesqp2'
          endif
        case default
        end select

C   ... iblstf not allocated:  substitute with iblst
        if (nblstf < 0) then
          call ilst2a(iblst,nblst,strn)
        else
          call ilst2a(iblstf,nblstf,strn)
        endif

        if (lbin) then
          print *, 'binary write not implemented, sorry'
          goto 10
C          if (j2 < j1) fn = 'seib'
C          ifi = fopna(fn,-1,4)
C          rewind ifi
C          call ywrm(1,' ',1,ifi,' ',seitp,0,nomg,nomg,nsei)
C          call fclose(ifi)
C          call info0(2,0,0,' wrote sei to binary file '//trim(fn))
        else
          ifi = fopna(fn,-1,0)
          rewind ifi
          kb = 0
          do  while (kb < nblst)
C         Write ib-summed PE spectra
          if (havepe == 1) then
            kb = nblst
            call awrit4('# ib='//trim(strn)//' qp(int)=%3:1;6,6d qp(meas)=%3:1;6,6d'//
     .        '%?#(n==1)# eqp=%d##',outs,len(outs),ifi,qpbnd,qsave,nblst,eigitp_save)
              write(ifi,324)
  324         format('#',5x,'omega',10x,'PES(w)',9x,'A(w)',11x,'A0(w)')

C         Write ib-summed spectral function or DOS
          elseif (havesi <= 2 .and. .not. lresib) then
            kb = nblst
            call awrit5('# ib='//trim(strn)//
     .        ' qp =%3:1;6,6d  eps=%;5d  %?#(n==2)##eqp=%n:1;6;6d#',
     .        outs,len(outs),ifi,qpi,eps,havesi,nblst,eigkb)
            write(ifi,323)
  323       format('#',5x,'omega',9x,'A(w)',11x,'A0(w)')

C         Write ib-resolved spectral function or dos
          elseif (havesi == 1) then
            kb = kb+1
            call awrit4('# ib=%i qp =%3:1;6,6d  eps=%;5d  eqp=%d',
     .        outs,len(outs),ifi,iblst(kb),qpi,eps,eigkb(kb))
            write(ifi,321)
            call awrit2('%% rows %i cols %i',' ',80,ifi,iw2-iw1+1,7)
C           endif
  321       format('#',5x,'omega',9x,'Re sigm-vxc',4x,'Im sigm-vxc',6x,
     .        'int A(w)',6x,'int A0(w)',7x,'A(w)',11x,'A0(w)')

C         Write sigma on Matsubara frequencies
          elseif (havesi == 11) then
            call ilst2a(iblst,nblst,strn)
            call awrit4('# beta=%d  nfreq=%i  ib='//trim(strn)//
     .        '  qp =%3:1;6,6d  eps=%;5,5d',
     .        outs,len(outs),ifi,beta,nmatsubara,qpi,eps)
            do  iwx = 1, nmatsubara
              write(ifi,"(' ',f12.7,1p,100(1x,2d15.7))") dimag(omgmats(iwx)),
     .          seimats(iwx,1:nblst)
            enddo

C         Missed something
          else
            call rx('sugws bug : mode not implemented')
          endif
          do  iwx = 1, nomgx
            if (havesi == 11) cycle
            if (iwx >= iw1 .and. iwx <= iw2) then
              if (havepe == 1) then
                write(ifi,"(' ',4d15.7)")
     .            omgx(iwx),pes(iwx),sag(iwx),sag0(iwx)
              elseif (havesi == 1 .and. .not. lresib) then
                write(ifi,"(' ',15d15.7)") omgx(iwx),sag(iwx),sag0(iwx)
              elseif (havesi == 1) then
                write(ifi,"(' ',15d15.7)") omgx(iwx),seitpw(iwx),
     .            sumag(iwx),sumag0(iwx),sagkb(iwx,kb),sag0kb(iwx,kb)
              elseif (havesi == 2) then
                write(ifi,"(' ',15d15.7)") omgx(iwx),dos(iwx,1),dos(iwx,2)
              endif
            endif
          enddo  ! frequency loop
          enddo  ! while kb<nblst

          call fclose(ifi)
          call info0(2,0,0,' wrote data to ascii file '//trim(fn))
        endif

        lnsave = .false.

C --- Abort ---
      elseif (outs(1:2) == 'a ') then
        call rx0('aborting sigfun editor ... no file written')

C --- quit ---
      elseif (outs(1:2) == 'q '. or. outs(1:5) == 'quit ') then
        if (lnsave) then
          print '('' se file not saved ... really quit?'')'
          read(*,'(a150)') outs
          call locase(outs)
          if (.not. (outs(1:1) == 'y' .or. outs(1:1) == 'q'))
     .      goto 10
        endif
        call rx0('exit sigfun editor')

C --- Help ---
      elseif (outs == '?') then
        print 310
        print 311
        print 312
        print 313
        print 314
        print 315
        print 316
  310   format(
     .    ' Select one of these options:'/
     .    t4,'readsek[flags]'/
     .    t4,'readsekb[flags]'/t12,
     .    'read spectral function from file, ASCII or binary format.'/
     .    t12,"File name is 'se.ext' or 'seb.ext' unless specified."/
     .    t12,"Optional flags are separated by a delimiter, which is the first character, e.g. '@'"/
     .    t14,"@fn=nam",t25,"use nam file file name"/
     .    t14,"@useef",t25,"specify that the file chemical potential be the Fermi level"/
     .    t14,"@irrmesh",t25,"points are not on a regular q mesh -- no q integrations"/
     .    t14,"@ib=list",t25,"restrict bands read from file to those in list"/
     .    t14,"@minmax",t25,"print minimum and maximum QP levels for each band"/
     .    t14,"@dc=#",t25,"subtract double counting # from Re sigma(omega) after reading"/
     .    t14,"@makeqpse",t25,"Not documented"//
     .
     .    t4,'evsync|evsync=#,#,#',t12,
     .    "replace QP levels read from 'se' by recalculating them"//
     .    t4,'units Ry|eV'/t12,
     .    'select units data for which data is to be returned.  Default is Ry.'//
     .    t4,'eps #',t12,'select smearing for spectral functions.'//
     .    t4,'ef #', t12,'Use # for Fermi level, overriding internally'
     .    ' calculated value.')

  311   format(/t4,
     .    'dos|jdos|imeps [ib=list] [nq=#1,#2,#3] ',
     .    '[nw=#|domg=#] [range=#1,#2] [a0] [kT=#] [isp=#]'/t12,
     .    'Integrate spectral function to make QP and spectrum DOS'/t12
     .    'Options:'/t12
     .    'ib    restrict contribution to A(w) from QP states in list'/t12,
     .    'nq    interpolate A(w) new mesh for q-integration'/t12,
     .    'nw    refine DOS to multiple of given energy mesh'/t12,
     .    'range generate dos in specified energy window'/t12,
     .    'kT    Temperature, units of omega (applies only to jdos and imeps)'/t12,
     .    'a0    Use noninteracting spectral function (only for jdos and imeps)'/t12,
     .    'isp   dos for given spin (default=1).  Applies to dos only.'
     .    )

  312   format(/t4,
     .    'se  iq=#|q=#1,#2,#3|allq|band[~args] ib=list|ibx=list [getev[=#1,#2,#3]] [a2qp] [nw=#|domg=#]',
     .    ' [isp=#] [range=#1,#2]'/t12,
     .    'Make sigma(omega) and A(omega) for given q and band.  Writes to spq.ext or spq2.ext'/t12
     .    'iq|q  single q : qp index, or q in Cartesian coordinates ','(units of 2pi/alat)'/t12,
     .    'allq  sigma is made for all q in the irreducible BZ and written to disk'/t12,
     .    'band  A(w),sigma(w) are made for qp along symmetry lines and written to disk'/t12,
     .    '      Optional ~args are parsed like options of the --band switch'/t12,
     .    '      In the multiple-q (latter) cases, data is written to disk as it is generated'/t12,
     .    'ib=   include contribution to A(w) from QP states in list'/t12,
     .    'ibx=  same as ib=, but resolve contribution to A(w) by band'/t12,
     .    'Options:'/t12
     .    'getev Do not interpolate energy but calculate it at q'/t12,
     .    'getev=#1,#2,#3 generates evals on independent mesh with nq=#1,#2,#3 divisions'/t12,
     .    'a2ap  Extract the QP energy from the peak in A and write it for the one-particle energy in spq.ext'/t12
     .    'nw    interpolate DOS to multiple of given energy mesh'/t12,
     .    'isp   functions for given spin (default=1)'/t12,
     .    'range generate sigma, A in specified energy window'
     .    )

  313   format(/t4,
     .    'pe|peqp|pe3d  iq=#|q=#1,#2,#3 ib=# [getev[=#1,#2,#3]]',
     .    ' [nw=#|domg=#] [nqf=#] [ek0=#] [isp=#] [range=#1,#2]'/t12,
     .    'Model ARPES for given q and band(s)'/t12
     .    'iq|q  qp index, or q in Cartesian coordinates ',
     .    '(units of 2pi/alat)'/t12,
     .    'ib    Include contribution to ARPES from QP-derived states in list'/t12,
     .    'Options:'/t12
     .    'getev Do not interpolate energy but calculate it at q'/t12,
     .    'getev=#1,#2,#3 generates evals on independent mesh'/t12
     .    '      with nq=#1,#2,#3 divisions'/t12,
     .    'nw    interpolate DOS to multiple of given energy mesh'/t12,
     .    'isp   functions for given spin (default=1)'/t12,
     .    'nqf   no. mesh points for final state intg (def=200)'/t12,
     .    'ke0   KE+V0=hbar*omega-phis+V0  KE of emitted electron'/t12,
     .    'range generate sigma, A in specified energy window'
     .    )

  314   format(/t4,
     .    'sem  allq|iq=#|q=#1,#2,#3|band=args ib=list beta=# nm=# [getev[=#1,#2,#3]] [nw=#|domg=#]',
     .    ' [isp=#] [range=#1,#2]'/t12,
     .    'Continue sigma to the Matsubara axis for given q'/t12,
     .    'iq|q  qp index, or q in Cartesian coordinates '/t12,
     .    'allq  special mode: sigma is made for all q in the irreducible BZ and written to disk'/t12,
     .    '(units of 2pi/alat)'/t12,
     .    'ib    include contribution to sigma(w) from QP states in list'/t12,
     .    'beta  Inverse temperature'/t12
     .    'nm=#  Number of Matsubara frequencies'/t12
     .    'Options:'/t12
     .    'getev Do not interpolate energy but calculate it at q'/t12,
     .    'getev=#1,#2,#3 generates evals on independent mesh'/t12
     .    '      with nq=#1,#2,#3 divisions'/t12,
     .    'nw    interpolate DOS to multiple of given energy mesh'/t12,
     .    'isp   functions for given spin (default=1)'/t12,
     .    'range generate sigma, A in specified energy window'
     .    )

  315   format(/
     .    t4,'qsgwh',t12,'Construct quasiparticle hamiltonian from lmf evals')

  316   format(/
     .    t4,'merges12'/t12,'merge spin1 and spin2 self-energies.'/
     .    t12,'Requires seq1.ext and seq2.ext reside on disk.'//
     .    t4,'savesea [fn]'/
     .    t4,'savese  [fn]'/t12,'saves q-interp self-energy, ',
     .    'ASCII or binary format.'/t12,
     .    "File name is 'seia.ext' or 'seib.ext' unless specified."//

     .    t4,'q',t12,'to quit the editor'/
     .    t4,'a',t12,'to abort')
      else
        call word(outs,1,j1,j2)
        call info0(0,0,0,' sfuned:  "'//trim(outs(j1:j2))//
     .    '" not recognized ... enter "?" for a list of options')
      endif
      goto 10

C   95 call info0(0,0,0,' sfuned:  illegal value at ...'//trim(outs(j1:))//
C     .  ' ... nothing done')
C      goto 10
   96 call info0(0,0,0,' sfuned:  missing tag '//trim(tag)//' in '//trim(outs)//
     .  ' ... nothing done')
      goto 10
   97 call info0(0,0,0,'%10pself-energy must be read before invoking this command ')
   98 call info0(0,0,0,' sfuned:  improper usage of '//trim(outs)//
     .  ' ... nothing done')
      goto 10
   99 call info0(0,0,0,' sfuned:  failed to parse arg starting at: '//
     .  trim(outs(j1:))//' ... nothing done')
      goto 10

      end

      subroutine spec_peaks(iw1,iw2,nomgx,omgx,ifi,ib,wpeak,eigitp,ag,dag,Zfac,ommax)
C- Return frequencies corresponding to peaks in the spectral function
C ----------------------------------------------------------------------
Ci Inputs
Ci  iw1,iw2:omgx(iw1:iw2) is range of omgx used by this routine
Ci   nomgx :dimensions omgx
Ci   omgx  :frequency mesh
Ci   ifi   :file logical unit where anomalies in peaks are written
Ci   ib    :band index (printout only)
Ci   wpeak :Input:  Initial value of omega for which to seek peak
Ci         :DESTROYED on output
Ci   eigitp:QP level
Ci   ag    :spectral function
Ci   Zfac  :Z factor (used for printout)
Co Outputs
Co   dag   :energy derivative of ag
Co   ommax :ommax(1) = energy of dominant peak
Co         :ommax(2) = energy of closest peak
Co         :ommax(3) = ratio of A(closest)/A(dominant)
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   16 Feb 13 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer iw1,iw2,nomgx,ifi,ib
      double precision wpeak,eigitp,Zfac,ommax(3)
      double precision omgx(nomgx),ag(nomgx),dag(nomgx)
C ... Local parameters
      logical T,F
      parameter (T=.true., F=.false.)
C     integer iwx,nomgx,npoly,jx1,jx2,iw1,iw2,ioff,nqe,nqbnd
      integer,allocatable :: iprm(:)
      real(8),allocatable:: zerom(:,:)
      integer jx1,nz,nbrack,iw,ilo,izeig,isw
      double precision dqp,dy,dwm,dwp
      character outs*128
      procedure(integer) :: iinear
      real(8),parameter:: NULLR =-99999

      outs = ' '

C ... Make A'(omega)
      call poldvm(omgx,ag,nomgx,6,F,1d-6,jx1,dag)

C ... Find all omega where A'(omega) vanishes
      if (jx1 /= 0) then
        nz = 0
      else
        nz = nbrack(dag(iw1),iw2-iw1+1,0d0)
      endif
      if (nz == 0) then
        call info2(20,1,0,' No peaks in A found in range (%d,%d)',omgx(iw1),omgx(iw2))
      else
        if (allocated(zerom)) deallocate(zerom,iprm)
        allocate(zerom(nz,2),iprm(nz))
        nz = 0
        iw = iw1
        dqp = 999 ! distance to QP level; find point closest
C       Store zeros in slope zerom(:,1); fn vals in zerom(:,2)
        do while (iw < iw2)
          call pshpr(1)
          call ytab0(0,omgx,dag,iw2,0d0,6,.false.,1d-8,wpeak,iw)
          call poppr
          if (iw < 0) exit
          nz = nz+1
          zerom(nz,1) = wpeak
          ilo = min(max(1,iw-6/2),iw2-6+1)
C         Value of A at this extremum
          call polinx(omgx(ilo),ag(ilo),6,zerom(nz,1),ag(iw)/1d6,zerom(nz,2),dy)
          iw = iw+1
C         Find extremum closest to eigitp
          if (abs(zerom(nz,1)-eigitp) < dqp) then
            dqp = abs(zerom(nz,1)-eigitp)
            izeig = nz
          endif
        enddo
        call dvheap(1,nz,zerom(1,2),iprm,0d0,1)
        ommax(1) = zerom(iprm(nz),1)
        ommax(2) = zerom(izeig,1)
        ommax(3) = zerom(izeig,2)/zerom(iprm(nz),2)
        call afwhm(omgx,nomgx,ag,ommax(1),zerom(iprm(nz),2),dwm,dwp)
        call awrit8('%x%2,4i%;11,6D%;11,6D%;11,6D%?#n#%j%5f----#%;9,4D#%?#n#%;9,4D#%j%5f----#',
     .    outs,len(outs),0,[ib,nz],eigitp,ommax(1),ommax(1)-eigitp,
     .    isw(Zfac==NULLR),Zfac,isw(dwp /= 0.and.dwm /= 0),dwp-dwm)
        if (nz > 1) then
          iw = iinear(nz,izeig,iprm,1)
          call awrit5('%a%;11,6D%;11,6D%:2;7F%;11,6D (%i)',outs,len(outs),0,
     .      zerom(iprm(nz-1),1),
     .      zerom(iprm(nz-1),1)-eigitp,
     .      zerom(iprm(nz-1),2)/zerom(iprm(nz),2),zerom(izeig,1),nz-iw+1)
        endif
        call info0(20,0,0,trim(outs))
      endif

      end

      integer function seqitp(mode,qpi,plat,nqfbz,qfbz,ipq,nomg,nband,
     .  nq,isp,iblst,nblst,eig,sigvq,sex,vxc,sigwq,sigi,eigi,sei,sexi,vxci)
C- Linearly interpolate self-energy to specified q-point from mesh points
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 do nothing
Ci         :bit 1  Interpolate eigi
Ci         :bit 2  Interpolate sigi
Ci         :bit 4  Interpolate sei, subtract sigi
Ci         :bit 8  Interpolate sei, do not subtract sigi
Ci         :bit 16 Interpolate sexi
Ci         :bit 32 Interpolate vxcdla
Ci         :bits 1,2 may be taken in combination with 4 or 8
Ci   qpi   :interpolate SE to this k-point
Ci   plat  :primitive lattice vectors, in units of alat
Ci   nqfbz :number of k-points in full BZ
Ci   qfbz  :k-points in full BZ
Ci   ipq   :ipq(i1,i2,i3) points to the irreducible qp into which
Ci         :mesh point (i1,i2,i3) is mapped (bzmesh.f)
Ci         :ipq(1)=0 => not used
Ci   nomg  :number of energy mesh points on which sigwq is given
Ci   nband :number of bands
Ci   nq    :number of k-points in the irreducible BZ.
Ci         :NB: if ipq(1) is zero, nq should be nqfbz!
Ci   isp   :current spin channel (1 or 2)
Ci   iblst :interpolate self-energy for bands in iblst
Ci   nblst :size of iblst
Ci   eig   :QP levels for each band, ib, spin
Ci   sigvq :QSGW VXC for each band, ib, spin
Ci   sigwq :self-energy for each band, ib, spin
Co Outputs
Co   sigi  :QSGW VXC interpolated to qpi
Co   eigi  :eig interpolated to qpi
Co   sei   :sigwq-vxc interpolated to qpi
Co   seqitp: 0 if interpolation is normal
Co         :-1 if interpolation failed
Co         : 1 if any of the eig(ib,iq,isp) is NULLI
Cl Local variables
Cl         :
Cr Remarks
Cu Updates
Cu   17 Jul 12
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nomg,nband,nq,isp,nqfbz,nblst,iblst(nblst),ipq(nqfbz),mode
      double precision qpi(3),qfbz(3,nqfbz),sigi(nblst),eigi(nblst),sexi(nblst),vxci(nblst)
      double precision eig(nband,nq,isp),sigvq(nband,nq,isp),plat(3,3)
      double precision sex(nband,nq,isp),vxc(nband,nq,isp)
      complex(8):: sigwq(nomg,nband,nq,isp),sei(nomg,nblst)
C ... Local parameters
      integer iprint
      integer i,j,k,iw,qplin,iqfbz(4,2),mode0,mode1,mode2,mode3,mode4,mode5,ib,kb
      double precision ddot,wt(4),xx,qp(3),qlat(3,3)
      complex(8):: cxx
      integer NULLI
      parameter (NULLI=-99999)

      mode0 = mod(mode,2)
      mode1 = mod(mode/2,2)
      mode2 = mod(mode/4,2)
      mode3 = mod(mode/8,2)
      mode4 = mod(mode/16,2)
      mode5 = mod(mode/32,2)

      if (mode == 0) return

        call pshpr(1)
C       print *, '!!'; qfbz(1:3,1) = (/.4d0,.5d0,.6d0/)
        seqitp = qplin(qpi,plat,nqfbz,qfbz,ipq,iqfbz,wt)
C       If no mapping to IBZ, substitute full BZ for irr
        if (ipq(1) == 0) then
          call icopy(4,iqfbz(1,1),1,iqfbz(1,2),1)
          if (nq /= nqfbz)
     .      call rx('seqitp: inconsistent use of mapping to irr mesh')
        endif
        call poppr
        call info5(41,0,0,' Interp SE for q=%3:1,5;5d, wt=%4:1,5;5d',qpi,wt,0,0,0)
        if (iprint() >= 45) then
          call mkqlat(plat,qlat,xx)
          do  j = 1, 4
            k = iqfbz(j,1)
            i = iqfbz(j,2)
            call shorbz(qfbz(:,k)-qpi(:),qp,qlat,plat)
            xx = dsqrt(ddot(3,qp,1,qp,1))
            call info8(40,0,0,' qi(%i) =%3;11,6D  |qi-qpi|=%,6;6d  iq=%i irr=%i wt=%,6;6d',
     .        j,qfbz(1,k),xx,iqfbz(j,1),i,wt(j),0,0)
          enddo
        endif
        if (seqitp < 0) then
          call info2(20,0,0,' wt gen failed for qp=%3:1,5;5d, wt =%4:1;9F',qpi,wt)
          return
      endif

      do  kb = 1, nblst
      ib = iblst(kb)

C ... List of corresponding irr QP and interpolated eval
      if (mode0 /= 0) eigi(kb) = 0
      if (mode1 /= 0) sigi(kb) = 0
      if (mode4 /= 0) sexi(kb) = 0
      if (mode5 /= 0) vxci(kb) = 0
      do  j = 1, 4
        if (eig(ib,iqfbz(j,2),isp) == NULLI) then
          seqitp = 1
          return
        endif
        if (mode0 /= 0)
     .  eigi(kb) = eigi(kb) + wt(j)*eig(ib,iqfbz(j,2),isp)
        if (mode1 /= 0)
     .  sigi(kb) = sigi(kb) + wt(j)*sigvq(ib,iqfbz(j,2),isp)
        if (mode4 /= 0)
     .  sexi(kb) = sexi(kb) + wt(j)*sex(ib,iqfbz(j,2),isp)
        if (mode5 /= 0)
     .  vxci(kb) = vxci(kb) + wt(j)*vxc(ib,iqfbz(j,2),isp)

C        if (mode0 /= 0)
C     .    print *, "!!e", ib, iqfbz(j,2), eig(ib,iqfbz(j,2),1:2)
C        if (mode1 /= 0)
C     .    print *, "!!s", ib, iqfbz(j,2), sigvq(ib,iqfbz(j,2),1:2)
      enddo

C ... Make interpolated self-energy; subtract vxc
      if (mode2 /= 0 .or. mode3 /= 0) then
        do  iw = 1, nomg
          cxx = 0
          do  j = 1, 4
            cxx = cxx + wt(j)*sigwq(iw,ib,iqfbz(j,2),isp)
          enddo
          if (mode2 /= 0) then
            sei(iw,kb) = cxx - sigi(kb)
          else
            sei(iw,kb) = cxx
          endif
        enddo
      endif

      enddo

      end

      subroutine sfitrp(mode,nomg,nomgx,ommin,ommax,multi,eigq,seiq,
     .  eps,omgn,omgx,seiqw,sumag,sumag0,ag,ag0,Zfac,jx1,jx2)
C- Interpolate one self-energy to finer energy mesh; make spectral function
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 do nothing
Ci         :1 generate omgn,omgx
Ci         :2 make spectral functions
Ci         :3 Both 1 and 2
Ci         :4 make spectral ag0 only
Ci         :5 Both 1 and 4
Ci         :10s digit:
Ci         :0 Integrate spectral function by trapezoidal rule
Ci         :1 Integrate spectral function using Simpson's rule
Ci         :100s digit:
Ci         :0 Make spectral function from [omega - (E_qp + sigm(omega) - vxc + i*eps)]^-1;
Ci         :  massage it by integration, then differentiation.
Ci         :1 Make spectral function directly
Ci   nomg  :Number of frequencies on which self-energy is given
Ci   nomgx :Number of frequencies on which to interpolate sigma, A
Ci   ommin :energy of seiq(1)
Ci   ommax :energy of seiq(nomg)
Ci   multi :Divide seiq into this number of divisions
Ci   eigq  :QP energy
Ci   seiq  :Given self-energy - vxc
Ci   eps   :Smearing (added to Im Sigma for smearing)
Cio Inputs/Outputs
Cio  omgn  :mesh of frequencies on which input seiq is defined
Cio        :Output if 1s digit mode is set; else input
Cio  omgx  :mesh of frequencies on which output seiqw,ag,ag0 are defined
Cio        :Output if 1s digit mode is set; else input
Co Outputs
Co   The following are generated only for mode>1:
Co   seiqw :Self-energy on interpolated mesh, QP part subtracted so that
Co         :G has a pole where seiqw=0
Co   sumag :Integral of spectral function on interpolated mesh
Co   sumag0:Integral of QP spectral function on interpolated mesh
Co   ag    :Spectral function on interpolated mesh
Co   ag0   :QP Spectral function on interpolated mesh
Co   jx1   :Case mode>1:
Co         :0 if differentiation to make ag with rational functions
Co         :1 if differentiation to make ag with polynomial
Co   jx2   :0 if differentiation to make ag0 with rational functions
Co         :1 if differentiation to make ag0 with polynomial
Cs Command-line switches
Cl Local variables
Cl  awi(iw,i,j) : intermediate quantities for spectral function
Cl          iw = frequency point
Cl          i = 1 for A at freq point; i = 2 for A at midpoint to next pt
Cl          j = 1 for A, 2 for A0
Cr Remarks
Cr
Cu Updates
Cu   23 Jul 12 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nomg,nomgx,multi,jx1,jx2
      double precision ommin,ommax,eigq,eps,Zfac
      complex(8):: seiq(nomg),seiqw(nomgx)
      double precision omgn(nomg),omgx(nomgx)
      double precision sumag(nomgx),sumag0(nomgx),ag(nomgx),ag0(nomgx)
C ... Local parameters
      logical nofilter
      double precision dw0,dw,dymx,sigwr,sigwi,dy,pi,sumchk,sumchk2,aw,exx
      double precision seiqr(nomg),seiqi(nomg)
      real(8),allocatable:: awi(:,:,:)
      complex(8):: den,cxx
      integer iwx,iw,im,npoly,mode0
      parameter (npoly=6,dymx=1d-5)

      allocate(awi(nomgx,2,2))

      mode0 = mod(mode,10)
!     mode2 = mod(mode/100,10)
      nofilter = mod(mode/100,10) == 1
      if (mode0 == 0) return
      if (nomg < 0 .or. multi < 1) call rx('sfitrp: bad frequency mesh')
      dw0 = 0; if (nomg > 1) dw0 = (ommax-ommin)/max(nomg-1,1)
      dw = dw0/multi
      pi = 4*datan(1d0)

C ... Make energy meshes
      if (mod(mode0,2) == 1) then
      iwx = 0
      do  iw = 1, max(nomg-1,1)
        omgn(iw) = ommin + dw0*(iw-1)
        do  im = 0, multi
          if (im == multi .and. iw /= nomg-1) cycle
          iwx = iwx+1
          omgx(iwx) = ommin + dw0*(iw-1 + im/dble(multi))
C          if (exx>=range(1) .and. exx<=range(2)) then
C            iws = iws+1
C          endif
        enddo
      enddo
      omgn(nomg) = ommax
      if (iwx > nomgx) call rx('sfitrp: dimensioning error')
      endif
      if (mode0 <= 1) return

C     call brack0(0,omgx,nomgx,eigq,iwx); print *,iwx; stop

C ... Interpolate self-energy to specified omega-mesh
      sumchk =  0
      sumchk2 = 0
      jx1 = 0
      seiqr = dble(seiq); seiqi = dimag(seiq)
      do  iwx = 1, nomgx

        exx = (omgx(iwx) + omgx(min(iwx+1,nomgx)))/2

        if (mode0 <= 3) then
C       Interpolate sigma to omgx(iwx)
        call polint(omgn,seiqr,nomg,npoly,omgx(iwx),dymx,3,jx1,sigwr,dy)
        call polint(omgn,seiqi,nomg,npoly,omgx(iwx),dymx,3,jx1,sigwi,dy)
        seiqw(iwx) = dcmplx(sigwr,sigwi)

C   ... Denominator = omega - (E_qp + sigm(omega) - vxc + i*eps)
        den = dcmplx(dble(omgx(iwx) - (Eigq + seiqw(iwx))),
     .               abs(dimag(seiqw(iwx)))+eps) ! eps is smearing

C   ... Integrate by trapezoidal rule A(omgx) = 1/pi abs(dimag(1d0/den))
C       and accumulate spectrum DOS
        awi(iwx,1,1) = abs(dimag(1/pi/den)) !A[omgx(iwx)]
        aw = abs(dimag(1/pi/den))
        sumchk = sumchk + dw*aw
        sumag(iwx) = sumchk
C       print 333, iwx,omgx(iwx),den,-1/den/pi

        call polint(omgn,seiqr,nomg,npoly,exx,dymx,3,jx1,sigwr,dy) ! Midpoint
        call polint(omgn,seiqi,nomg,npoly,exx,dymx,3,jx1,sigwi,dy)
        cxx = dcmplx(sigwr,sigwi)
        den = dcmplx(dble(exx - (Eigq + cxx)),
     .               abs(dimag(cxx))+eps) ! eps is smearing
        awi(iwx,2,1) = abs(dimag(1/pi/den))  !A[(omgx(iwx)+omgx(iwx+1))/2]
        else
        awi(iwx,1,1) = 0
        awi(iwx,2,1) = 0
        endif

C   ... Integrate by trapezoidal rule QP DOS
        den = dcmplx(dble(omgx(iwx)-eigq), eps) !eps = smearing
        awi(iwx,1,2) = abs(dimag(1/pi/den))  !A0[omgx(iwx)]
        aw = abs(dimag(1/pi/den))

        den = dcmplx(dble(exx-eigq), eps) !eps = smearing
        awi(iwx,2,2) = abs(dimag(1/pi/den)) !A0[(omgx(iwx)+omgx(iwx+1))/2]

C       dosqp(iwx) = dosqp(iwx) + wibz(iq) * aw
        sumchk2 = sumchk2 + dw*aw
        sumag0(iwx) = sumchk2

C       dos(iwx) = dos(iwx) +  wibz(iq) * aw
C       if (iwx == 1) print *, iq,wibz(iq)

      enddo

C ... Integrate by trapezoidal rule A(omgx) = 1/pi*abs(dimag(1d0/den))
      sumchk = 0; sumchk2 = 0; sumag(1) = 0; sumag0(1) = 0
C     sumchk = dw*awi(1,1,1); sumchk2 = dw*awi(1,1,2)
      do  iwx = 2, nomgx
C       Constant rule
C        sumchk = sumchk + dw*awi(iwx,1,1)
C        sumchk2 = sumchk2 + dw*awi(iwx,1,2)
C       Trapezoidal rule
C        sumchk = sumchk + dw/2*(awi(iwx,1,1)+awi(iwx-1,1,1))
C        sumchk2 = sumchk2 + dw/2*(awi(iwx,1,2)+awi(iwx-1,1,2))
C       Simpson's rule
        if (mode0 <= 3) then
          sumchk = sumchk + dw/6*(awi(iwx,1,1)+awi(iwx-1,1,1)+4*awi(iwx-1,2,1))
        endif
        sumchk2 = sumchk2 + dw/6*(awi(iwx,1,2)+awi(iwx-1,1,2)+4*awi(iwx-1,2,2))

C       print *, sumag(iwx) - sumchk
C       print *, sumag0(iwx) - sumchk2
        sumag(iwx) = sumchk
        sumag0(iwx) = sumchk2
      enddo

C ... Z factor evaluated at eigq
C     call prrmsh('Re sig',omgn,seiqr,nomg,nomg,1)
      if (mode0 <= 3 .and. nomg > 1) then
        call poldvm(omgn,seiqr,nomg,6,.false.,1d-6,iw,ag)
        call polint(omgn,ag,nomg,npoly,eigq,dymx,3,jx1,Zfac,dy)
        Zfac = 1/(1-Zfac)
      else
        Zfac = 1
      endif
C     print *, eigq,Zfac,dy; stop

C ... Generate A[G(w)] and A[G0(w)] by numerical differentiation
C     jx=0 => rational interpolation (attempt first), else polynomial
      if (mode0 <= 3 .and. nofilter) then
        forall (iw=1:nomgx) ag(iw) = awi(iw,1,1)
      else if (mode0 <= 3) then
        call poldvm(omgx,sumag,nomgx,6,.true.,1d-6,jx1,ag)
        if (jx1 /= 0) then
          call poldvm(omgx,sumag,nomgx,6,.false.,1d-6,iw,ag)
        endif
      endif
      if (nofilter) then
        forall (iw=1:nomgx) ag0(iw) = awi(iw,1,2)
      else
        call poldvm(omgx,sumag0,nomgx,6,.true.,1d-6,jx2,ag0)
        if (jx2 /= 0) then
          call poldvm(omgx,sumag0,nomgx,6,.false.,1d-6,iw,ag0)
        endif
      endif

C      do  iwx = 1, nomgx
C        print 333, iwx,omgx(iwx),awi(iwx,1,1),ag(iwx)
C      enddo
  333 format(i4,f12.6,2x,2f12.6,2x,2f12.6)

      end
      subroutine schkpole(opt,nband,nq,nspse,nomg,ommin,ommax,eigs,sigwq,sigvq)
C- Shift sigvq or sigwq to enforce QSGW condition, (sig-vxc)_w=eQP=0
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :1s digit printout
Ci         :1 print out avg shift in sigvq, and RMS deviation
Ci         :2 print out shift in sigvq for each state
Ci         :10s digit for sigvq
Ci         :0 sigvq has been folded into sigwq; sigvq not used
Ci         :1 Same as 0, except that sigwq is adjusted to reconcile sigwq=eigs
Ci         :2 sigwq does not contain static shift.
Ci         :  sigvq is subtracted from sigwq to find where it crosses zero
Ci            but sigvq left unchanged
Ci         :3 Same as 2, but sigvq is adjusted to reconcile sigwq-sigvq=eigs
Ci   nband :number of bands
Ci   nq    :number of k-points on a single symmetry line for bands plotting
Ci   nspse :2 for spin-polarized, collinear case, otherwise 1
Ci   nomg  :number of frequency points
Ci   ommin :lowest frequency
Ci   ommax :highest frequency
Ci   eigs  :QP levels
Cio Inputs/Outputs
Cio  sigwq :frequency dependent self-energy, or SE-vxc
Cio        :If 10s digit opt is 1, adjusted on output
Cio  sigvq :QSGW vxc
Cio        :Only used if 10s digit opt >= 2
Cio        :If 10s digit opt is 3, adjusted on output
Cl Local variables
Cl   seqp  : Value of sig(w)-vxc at eqp
Cl   shft  : change in w to make sig(w)-vxc=0
Cr Remarks
Cr
Cu Updates
Cu   25 Aug 12 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer opt,nomg,nband,nq,nspse
      double precision ommin,ommax
      double precision eigs(nband,nq,nspse),sigvq(nband,nq,nspse)
      double complex sigwq(nomg,nband,nq,nspse)
C ... Local parameters
      integer isp,iq,ib,iw,nglob,stdo,opt0,opt1,jx,nincl
      double precision xx,w0,eqp,seqp,dy,shft,sshft,sshft2,sseqp,sseqp2
      real(8),allocatable:: omgn(:),omgx(:),dos(:)

      stdo = nglob('stdo')
      opt0 = mod(opt,10)
      opt1 = mod(opt/10,10)
C      opt2 = mod(opt/100,10)
C      opt3 = mod(opt/1000,10)
C      opt4 = mod(opt/10000,10)


      allocate(omgn(nomg),omgx(nomg),dos(nomg))
      call sfitrp(1,nomg,nomg,ommin,ommax,1,xx,xx,xx,
     .  omgn,omgx,xx,xx,xx,xx,xx,xx,iw,iw)

      if (opt0 >= 2) then
        write(stdo,333)
  333   format('  ib  iq isp  iw     w(sig=0)      eQP',9x,'diff    Sig-vxc(eQP)')
      endif

      sshft=0; sshft2=0; sseqp=0; sseqp2=0; nincl=0
      do  isp = 1, nspse
        do  iq = 1, nq
          jx = 1
          do  ib = 1, nband
            eqp = eigs(ib,iq,isp)
            if (eqp < ommin .or. eqp > ommax) cycle
            nincl = nincl+1
C           dos <- sigm(w)
            call dcopy(nomg,sigwq(1,ib,iq,isp),2,dos,1)
C           Subtract sigvq from dos
            if (opt1 >= 2) call dvadd(dos,1,nomg,-sigvq(ib,iq,isp))
C           Closest energy w0 where Re(sig-vxc) crosses zero
            call pshpr(1)
            w0 = eqp
            call ytab0(2,omgn,dos,nomg,0d0,6,.false.,1d-10,w0,iw)
C           Value of sig(w)-vxc at eqp
            call polint(omgn,dos,nomg,6,eqp,0d0,0,jx,seqp,dy)
C           Get energy shift
C           call ytab0(2,omgn,dos,nomg,seqp,6,.false.,1d-10,evshft,iw)
            call poppr
            shft = eqp-w0
            sshft = sshft + shft         ! Average
            sshft2 = sshft2 + shft*shft  ! for root mean square
            sseqp = sseqp + seqp
            sseqp2 = sseqp2 + seqp*seqp
            if (opt0 >= 2) write(stdo,"(3i4,i5,4f12.6)")
     .        ib,iq,isp,iw,w0,eigs(ib,iq,isp),shft,seqp
            if (opt1 == 1) then
              call dvadd(dos,1,nomg,-seqp)
              call dcopy(nomg,dos,1,sigwq(1,ib,iq,isp),2)
            elseif (opt1 == 3) then
              sigvq(ib,iq,isp)= sigvq(ib,iq,isp)+seqp
            endif
          enddo
        enddo
      enddo
      if (nincl == 0) return
      sshft2 = sshft2/nincl; sshft = sshft/nincl
      sshft2 = dsqrt(sshft2-sshft**2+1d-16)
      sseqp2 = sseqp2/nincl; sseqp = sseqp/nincl
      sseqp2 = dsqrt(sseqp2-sseqp**2+1d-16)

C --- Printout summary ---
      if (opt0 == 0 .or. nincl == 0) return
      if (nincl == 1) then
        call info5(20,0,0,
     .    ' QSGW condition for %i pts  (sig-vxc)@eQP>=%d  shift <dw>=%d  adjust=%l',
     .    nq,-sseqp,sshft,mod(opt1,2) /= 0,0)
      else
        call info8(20,0,0,
     .    ' QSGW condition for %i pts  (sig-vxc)@eQP>=%d RMS=%d  shift <dw>=%d RMS=%d  Adjust=%l',
     .    nq,-sseqp,sseqp2,sshft,sshft2,mod(opt1,2) /= 0,0,0)
      endif
      end

      subroutine afwhm(omgx,nomgx,ag,wtop,atop,dwm,dwp)
C- Points of half maximum
C ----------------------------------------------------------------------
Ci Inputs
Ci   omgx
Ci   nomgx
Ci   ag    :translation part of space group
Ci   wtop
Ci   atop
Co Outputs
Ci   dwm
Ci   dwp
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   30 Aug 12 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nomgx
      double precision wtop,atop,omgx(nomgx),ag(nomgx),dwm,dwp
C ... Local parameters
      integer iwx,iw
C     double precision

      call pshpr(1)
      iwx = -1
      call brack0(0,omgx,nomgx,wtop,iwx)
      iw = iwx
      dwp = wtop
      call ytab0(0,omgx,ag,nomgx,atop/2,6,.false.,1d-10,dwp,iw)
      if (iw > 0) then
        dwp = dwp - wtop
      else
        dwp = 0
      endif
      iw = iwx+1
      dwm = wtop
      call ytab0(1,omgx,ag,nomgx,atop/2,6,.false.,1d-10,dwm,iw)
      if (iw > 0) then
        dwm = dwm - wtop
      else
        dwm = 0
      endif
      call poppr
      end
      subroutine susepade(mode,nomg,omg,ommin1,ommax1,ommin2,ommax2,npade,iprm)
C- Generate a list of points on which to make Pade interpolation
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 return npade
Ci         :1 return iprm
Ci   nomg  :total number of points for which sigma is available on the real axis
Ci   omg   :points for which sigma is available on the real axis
Ci   ommin2:Select geometric points between [ommin2,ommin1]; see Remarks
Ci   ommin1:Select all the points between [ommin1,ommax1];   see Remarks
Ci   ommax1:ditto
Ci   ommax2:Select geometric points between [ommax1,ommax2]; see Remarks
Cio Inputs/Outputs
Cio  npade :number of points to include in Pade interpolation
Cio        :(mode=0): output   (mode=1): input
Co Outputs
Co   iprm  :(mode=1) list of points, sorted
Cr Remarks
Cr   Algorithm:
Cr   interval  [ommin1,ommax1]  include up all points in
Cr   interval  [ommin1,ommin2]  points 1,3,7,15,... decreasing from ommin1 to ommin2
Cr   interval  [ommax1,ommax2]  points 1,3,7,15,... increasing from ommax1 to ommax2
Cu Updates
Cu   12 Nov 16 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nomg,npade,iprm(npade)
      double precision omg(nomg) !,omgp(npade)
      double precision ommin1,ommax1,ommin2,ommax2
C ... Local parameters
      integer,allocatable:: iwk(:)
      integer iw,deltai,npad0
      integer low1,low2,hi1,hi2

C     Find points closest to omega=0
      call huntx(omg,nomg,ommin1,[0],low1)
      call huntx(omg,nomg,ommin2,[0],low2)
      call huntx(omg,nomg,ommax1,[0],hi1)
      call huntx(omg,nomg,ommax2,[0],hi2)
C     Rationalize the limits
      low2 = max(min(low1,low2),1)
      hi1 = min(max(low1+1,hi1),nomg)
      hi2 = min(max(hi1-1,hi2),nomg)

      npad0 = 0
C     if (mode > 0) allocate(omgp(npade))

C ... All the points between ommin1 and ommax1
      do  iw = low1+1, hi1-1
        npad0 = npad0+1
        if (mode > 0) iprm(npad0) = iw
      enddo
C ... Geometric points between ommin2 and ommin1
      iw = low1+1; deltai = 1
C     do  iw = low2, low1
      do while (iw-deltai >= low2)
        iw = iw - deltai
        npad0 = npad0+1
        if (mode > 0) iprm(npad0) = iw
        deltai = 2*deltai       !points 1,3,7,15,...
      enddo
C ... Geometric points between ommax1 and ommax2
      iw = hi1-1; deltai = 1
C     do  iw = hi1, hi2
      do while (iw+deltai <= hi2)
        iw = iw + deltai
        npad0 = npad0+1
        if (mode > 0) iprm(npad0) = iw
        deltai = 2*deltai       !points 1,3,7,15,...
      enddo

      if (mode == 0) then
        npade = npad0
        call info5(10,0,0,
     .    ' susepade : limits (%d,%d,%d,%d) yield %i points',
     .    ommin1,ommax1,ommin2,ommax2,npade)
      else
        if (npade /= npad0) call rx('susepade: dimensioning error')
        allocate(iwk(npade))
        call ivheap(1,npade,iprm,iwk,0)
        deallocate(iwk)
      endif

      end

      subroutine sig2mats(mode,nomg,omgn,eigtp,seitp,npade,iprm,nmats,omgm,
     .  semats)
C- Analytically continue self-energy to complex frequencies
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :not used now
Ci   nomg  :Number of frequencies on which self-energy is given
Ci   omgn  :mesh of frequencies on which input seitp is defined
Ci   eigtp :QP energy
Ci   seitp :self-energy - vxc
Ci   npade :number of points to include in Pade interpolation (susepade)
Ci   iprm  :selection of npade points from omgn
Ci   nmats :number of Matsubara points (or generally, points where to continue seitp)
Ci   omgm  :omgm(nmats) = frequencies where to make semats
Co Outputs
Co   semats:Self-energies at omgm
Cr Remarks
Cu Updates
Cu   12 Nov 16 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nomg,nmats,npade,iprm(npade)
      double precision eigtp,omgn(nomg)
      complex(8):: omgm(nmats),seitp(nomg),semats(nmats)
C ... Dynamically allocated arrays
      complex(8), allocatable :: seip(:)        ! Self-energy at Pade points
      complex(8), allocatable :: zomgp(:)       ! Points from omg for Pade
      complex(8), allocatable :: padcoff(:,:,:) ! Pade coefficients
      complex(8), allocatable :: wk(:)          ! work array
C ... Local parameters
      logical :: ltest=.true.
      integer iw,mode0,i
      double precision rmserr

      mode0 = mod(mode,10); i = 1
      allocate(zomgp(npade),seip(npade),padcoff(npade,npade+2,i))
      do  iw = 1, npade
        zomgp(iw) = omgn(iprm(iw))
        seip(iw)  = seitp(iprm(iw))
      enddo
C      call prrmsh('Re',dble(zomgp),dble(seip),npade,npade,1)
C      call prrmsh('Im',dble(zomgp),dimag(seip),npade,npade,1)

      call padcof(npade,zomgp,seip,npade,padcoff)

C ... test the reliability of the Pade fit
      if (ltest) then
        rmserr = 0
        allocate(wk(npade))
        call pade(npade,zomgp,npade,zomgp,npade,padcoff(1,1,i),wk)
        do  iw = 1, npade
          rmserr = rmserr + (wk(iw)-seip(iw))*dconjg(wk(iw)-seip(iw))
        enddo
        rmserr = sqrt(rmserr/npade)
        deallocate(wk)
        call info2(10,0,0,' sig2mats: check Pade fit (%i points) rms err=%g',npade,rmserr)
      endif

C ... Analytically continue to Matsubara (or any complex) frequency
      call pade(nmats,omgm,npade,zomgp,npade,padcoff(1,1,i),semats)

      deallocate(zomgp,seip,padcoff)
      end

      subroutine mergesegw12(ifi,jfi,kfi)
C- Merge spin1 and spin2 spectral functions, contained in the separate se file formats
C ----------------------------------------------------------------------
Ci Inputs
Ci   ifi   : spin 1 file logical unit
Ci   jfi   : spin 2 file logical unit
Ci   kfi   : merged file logical unit
Co Outputs
Co   seq1 (file ifi) and seq2 (file jfi) are merged into seq (file kfi)
Cr Remarks
Cr  The se file format is documented here: /docs/input/data_format/#the-se-file
Cu Updates
Cu   20 Nov 16  First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ifi,jfi,kfi
C     integer mode,isp,nsp,nblst,iblst(nblst),nkp,ifi
C ... Dynamically allocated arrays
      integer, allocatable :: nqlst(:),ifblst(:)
      real(8),allocatable:: eigiq(:,:,:),sexiq(:,:,:),vxciq(:,:,:),sigiq(:,:,:),qp(:,:)
      complex(8),allocatable:: seq(:,:,:)
C ... Local parameters
      logical ltmp
      real(8), parameter:: tolq=1d-7
      integer fmode,i,lascii,i1,i2,i3,i4,i5,i6,i7,iq,nomg,nfblst,nband,nqp,nsp,npan,nspse
      double precision xv(10),ommin,ommax,chempot
      procedure(logical) :: isanrg
      procedure(integer) :: iinear

C --- Read header ---
      lascii = 20
      i = 3000 + 200 + lascii + 4
      rewind ifi
      call ioseh(i,fmode,nsp,nspse,nband,nqp,nomg,ommin,ommax,chempot,nfblst,[0],npan,ifi)
      ltmp = isanrg(nsp,1,1,'iosegw:','file''s nsp',.true.)
      rewind jfi
      call ioseh(i,i1,i2,nspse,i3,i4,i5,xv(1),xv(2),i6,chempot,[0],i7,jfi)
C     Header info much match
      if (fmode /= i1  .or. nsp /= i2 .or. nband /= i3 .or. nqp /= i4 .or. nomg /= i5 .or.
     .  ommin /= xv(1) .or. ommax /= xv(2) .or. nfblst /= i6 .or. npan /= i7)
     .  call rx('mergesegw12 : file mismatch')

      if (fmode /= 32) call rx('mergesegw12 can only handle fmode=32 for now')

      allocate(ifblst(nfblst),nqlst(0:npan))
      i = i + 4000
      rewind ifi; rewind jfi
      call ioseh(i,fmode,nsp,nspse,nband,nqp,nomg,ommin,ommax,chempot,nfblst,ifblst,nqlst,ifi)
      call ioseh(i,i1,i2,nspse,i3,i4,i5,xv(1),xv(2),chempot,i6,ifblst,nqlst,jfi) ! A little buglet here ... should check ifblst, nqlst

      rewind kfi
      call ioseh(i,fmode,2,nspse,nband,nqp,nomg,ommin,ommax,chempot,nfblst,ifblst,nqlst,-kfi)

C ... Read, merge and write q, eqp, sigqp.  Buglet: no check is made on the consistency of qp
      allocate(qp(3,nqp),sigiq(nfblst,nqp,2),eigiq(nfblst,nqp,2),sexiq(nfblst,nqp,2),vxciq(nfblst,nqp,2))
      call iose2(lascii,fmode,fmode,1,nfblst,nqp,nomg,qp,eigiq(1,1,1),sigiq(1,1,1),xv,sexiq(1,1,1),vxciq(1,1,1),ifi)
      call iose2(lascii,fmode,fmode,1,nfblst,nqp,nomg,qp,eigiq(1,1,2),sigiq(1,1,2),xv,sexiq(1,1,2),vxciq(1,1,2),jfi)
      call iose2(lascii,fmode,fmode,2,nfblst,nqp,nomg,qp,eigiq(1,1,1),sigiq(1,1,1),xv,sexiq(1,1,1),vxciq(1,1,1),-kfi)

C --- Read frequency dependent self-energy ---
      allocate(seq(nomg,nfblst,2))
      do  iq = 1, nqp
        call dfdump(seq(1,1,1),2*nomg*nfblst,ifi)
        call dfdump(seq(1,1,1),2*nomg*nfblst,-kfi)
        call dfdump(seq(1,1,2),2*nomg*nfblst,jfi)
        call dfdump(seq(1,1,2),2*nomg*nfblst,-kfi)
      enddo
      deallocate(seq)
      call rx0('completed spin merge')

      end

      subroutine qpse(iw1,iw2,nomg,ommin,ommax,nband,nqp,nsp,eigs,sigwq,sgq)
C- Find the 'quasiparticle self-energy' for one band and k-point
Cu Updates
Cu   15 Apr 17 First attempt ... still under construction
      implicit none
C ... Passed parameters
      integer iw1,iw2,nomg,nband,nqp,nsp
      double precision ommin,ommax
      double precision eigs(nband,nqp,nsp),sgq(nband,nqp,nsp)
      complex(8):: sigwq(nomg,nband,nqp,nsp)
C ... Dynamically allocated arrays
      real(8),allocatable:: omgn(:),omgx(:)
      real(8),allocatable:: sigom(:,:)
C ... Local parameters
      integer isp,iq,ib,iw
      double precision xx,z,omgi(2),wlo,whi
      equivalence (wlo,omgi(1)),(whi,omgi(2))
      integer, parameter :: NULLI=-99999

      allocate(omgn(nomg),omgx(nomg))
      call sfitrp(1,nomg,nomg,ommin,ommax,1,xx,xx,xx,
     .  omgn,omgx,xx,xx,xx,xx,xx,xx,iw,iw)

C     Debugging : dump (omg,sig) to file for one qp and spin
      allocate(sigom(nomg,0:nband))
      iq = 1; isp = 1
      do  iw = 1, nomg
        sigom(iw,0) = omgn(iw)
        forall (ib=1:nband) sigom(iw,ib) = sigwq(iw,ib,iq,isp)
      enddo
      call prmx('sigw',sigom,nomg,nomg,nband+1)

C      call rx('qpse not ready')

      do  isp = 1, nsp
        do  iq = 1, nqp
          do  ib = 1, nband
C         do  ib = 32, 35
            call qpsei(iw1,iw2,omgn,eigs(ib,iq,isp),sigwq(1,ib,iq,isp),omgi,z,sgq(ib,iq,isp))
            if (wlo == NULLI .and. whi == NULLI) then
              print 333, ib,iq,'none',eigs(ib,iq,isp)
            elseif (wlo == NULLI) then
              print 333, ib,iq,'whi', eigs(ib,iq,isp),whi-eigs(ib,iq,isp)
            elseif (whi == NULLI) then
              print 333, ib,iq,'wlo', eigs(ib,iq,isp),wlo-eigs(ib,iq,isp)
            else
              print 333, ib,iq,'both', eigs(ib,iq,isp),wlo-eigs(ib,iq,isp), whi-eigs(ib,iq,isp)
            endif
          enddo
        enddo
      enddo
  333 format(2i4,2x,a4,1x,3f12.6)

      stop
      end
      subroutine qpsei(iw1,iw2,omg,eig,sigw,omgi,z,sgq)
C- Find the 'quasiparticle self-energy' for one band and k-point
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :
Ci  iw1,iw2:omgx(iw1:iw2) is range of omgx used by this routine
Ci   nomgx :dimensions omgx
Ci   omgx  :frequency mesh
Ci   ifi   :file logical unit where anomalies in peaks are written
Ci   ib    :band index (printout only)
Ci   wpeak :Input:  Initial value of omega for which to seek peak
Ci         :DESTROYED on output
Ci   eigitp:QP level
Ci   ag    :spectral function
Co Outputs
Co   sigqp :potential where (omega - e - sigw(omega)) crosses zero
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr   Should return z vactor, sgq
Cu Updates
Cu   15 Apr 17 First attempt
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer iw1,iw2
      double precision eig,omgi(2),sgq,z
      double precision omg(iw2)
      complex(8) :: sigw(iw2)
C ... Local parameters
      integer, parameter :: NULLI=-99999
      logical, parameter :: T=.true., F=.false.
C     integer iwx,nomgx,npoly,jx2,iw1,iw2,ioff,nqe,nqbnd
      integer iw,low1
      double precision linri,wlo,whi
C     Find x that causes linear interpolation of f to vanish
      double precision x1,x2,f1,f2
C     Linearly interpolate a function from data at 2 points
C     linr(x,x1,x2,f1,f2) = ((x-x1)*f2-(x-x2)*f1)/(x2-x1)
C     Find a zero of a linearly interpolated function
      linri(x1,x2,f1,f2) = (x1*f2-x2*f1)/(f2-f1)

C      print *, linr(2d0,1d0,2d0,2*1d0-.5d0,2*2d0-.5d0)
C      print *, linr(1d0,1d0,2d0,2*1d0-.5d0,2*2d0-.5d0)
C      print *, linr(.0d0,1d0,2d0,2*1d0-.5d0,2*2d0-.5d0)
C      print *, linr(.25d0,1d0,2d0,2*1d0-.5d0,2*2d0-.5d0)
C      print *, linri(1d0,2d0,2*1d0-.5d0,2*2d0-.5d0)
C      stop

C     Find point closest to omega=0
      low1 = 1
      call huntx(omg,iw2,eig,[0],low1)

C     Look for crossing in omega - e - sigw(omega) for interval [low1,nomg]
      whi = NULLI
      do  iw = low1, iw2-1
        if ((omg(iw) - eig - dble(sigw(iw)))*(omg(iw+1) - eig - dble(sigw(iw+1))) < 0) then
          whi = linri(omg(iw),omg(iw+1),omg(iw) - eig - dble(sigw(iw)),omg(iw+1) - eig - dble(sigw(iw+1)))
!         print *, whi,linr(whi,omg(iw),omg(iw+1),omg(iw) - eig - dble(sigw(iw)),omg(iw+1) - eig - dble(sigw(iw+1)))
          exit
        endif
      enddo

C     Look for crossing in omega - e - sigw(omega) for interval [1,low1]
      wlo = NULLI
      do  iw = low1, iw1+1, -1
C       print *, iw, (omg(iw) - eig - dble(sigw(iw))),(omg(iw-1) - eig - dble(sigw(iw-1)))
        if ((omg(iw) - eig - dble(sigw(iw)))*(omg(iw-1) - eig - dble(sigw(iw-1))) < 0) then
          wlo = linri(omg(iw),omg(iw-1),omg(iw) - eig - dble(sigw(iw)),omg(iw-1) - eig - dble(sigw(iw-1)))
!         print *, wlo,linr(wlo,omg(iw),omg(iw-1),omg(iw) - eig - dble(sigw(iw)),omg(iw-1) - eig - dble(sigw(iw-1)))
          exit
        endif
      enddo

      omgi(1) = wlo; omgi(2) = whi

      end

