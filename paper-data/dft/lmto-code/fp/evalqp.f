      subroutine evalqp(mode,s_bz,s_ctrl,s_ham,s_lat,s_mix,s_move,s_pot,
     .  s_str,s_spec,s_site,s_optic,s_gw,s_dmft,s_strn,nband,ioff,nsp,ldq,
     .  nkabc,nq,nqoff,lunits,ef0,eig,qp,wbz,ipq,efshft)
C- qp, and optionally eigenvalues, Fermi level for a k-mesh or list of k
C ----------------------------------------------------------------------
Cio Structures
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  nkp nkabc ef egap lshft lmet nevmx def w dos idtet
Ci                 ipq pdos qp star wtkp wtkb swtk ntet n efmax fsmom
Ci                 ndos dosw numq range
Co     Stored:     nkp nkabc star ntet nevmx ef def w dos idtet ipq
Co                 pdos qp wtkp wtkb swtk numq ndos dosw
Co     Allocated:  qp ipq wtkp idtet star wtkb swtk
Cio    Elts passed: ipq qp wtkp lopt lio idtet star numq swtk wtkb egap
Cio                n w def
Cio    Passed to:  mkqp lmfp popted rsedit iorsf bcast_strx iinit mpibc1
Cio                chimedit bndfp rdsigm subzi addrbl optinq optin2
Cio                optint
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  plbnd lpgf lmet lrs nbas nl nspec nspin npadl npadr
Ci                 lcgf lscr lsx zbak maxit lrel lfrce nitmv mdprm ltb
Ci                 quit tol nvario defm nbasp lham lgen3 lncol ldlm
Ci                 nclasp ipc loptc pfloat ldos lfp
Co     Stored:     plbnd lmet lfrce mdprm lrs maxit ltb
Co     Allocated:  *
Cio    Elts passed: lmet ldos lsx lscr lrs lcd lbas nvario lfp lgen3
Cio                ips ipc ncomp idcc lncol
Cio    Passed to:  mkqp lmfp popted rlxstp supot suham pp2alp rsedit
Cio                iorsf bcast_strx iinit mpibc1 chimedit smshft bndfp
Cio                suham2 optinq optin2 optint vcdmel rxes dfrce relax
Cio                lattic
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  ndham evals ehf ehk seref eterms pwmode pwemin
Ci                 pwemax npwpad lncol neula qss nbf elind ldham lsig
Ci                 oveps ovncut rsrnge nqsig eseavr sigp pmin pmax
Ci                 pnudef rsstol nprs qsig ndhrs
Co     Stored:     nlibu lmaxu ndham ndofH ldham lmxax npwmin npwpad
Co                 hord lsig eterms sigp ehf ehk nqsig ndhrs eseavr
Co     Allocated:  offH iprmb bdots nprs hrs qsig iaxs
Cio    Elts passed: evals seref offH iprmb eula magf bdots qsig lncol
Cio                hrs iaxs nprs
Cio    Passed to:  lmfp suham smshft bndfp mkpot rdsigm hft2rs suham2
Cio                hambls hambl sopert3 blsig makusq addrbl mkrout
Cio                mkehkf mkekin chkdmu
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  plat nsgrp npgrp gam plat0 alat qlat gmax vol awald
Ci                 nkd nkq nabc ng kv kv2 igv igv2 tolft tol dist ag
Ci                 bgv cg cy dlv gv gvq indxcg ips0 istab jcg pos qlv
Ci                 symgr afmt napw as rpad nkdmx nkqmx platl platr ldist
Co     Stored:     gam ldist plat dist gmax nabc ng igv igv2 alat ag
Co                 bgv cg cy dlv gv gvq indxcg ips0 istab jcg kv kv2
Co                 pos qlv symgr napw vol plat0 qlat platl platr awald
Co                 nkd nkq
Co     Allocated:  ips0 bgv gv kv igv kv2 igv2 dlv qlv
Cio    Elts passed: symgr pos istab ag dlv qlv nsgrp gv ips0 bgv qlat
Cio                igv igv2 kv cg indxcg jcg cy vol plat
Cio    Passed to:  mkqp lmfp popted lmfopb supot sugvec0 sugvec suham
Cio                rsedit rdovfa ovlpfa ovlocr hxpbl ghibl hklbl gklbl
Cio                iorsf bcast_strx iinit mpibc1 prsed2 prsed4 chimedit
Cio                smshft pvsms1 rhgcmp symsmr bndfp mkpot smves vesft
Cio                vesgcm mshvmt symvvl ugcomp ggugbl gfigbl fklbl
Cio                hhugbl hhigbl phhigb hsmbl hgugbl smvxc2 vxcnlm
Cio                smvxcm smcorm locpot rdsigm suham2 hambls hambl
Cio                augmbl bstrux hxpgbl ghigbl hklgbl smhsbl hhibl
Cio                phhibl hsibq makusq pusq1 addrbl fsmbl rsibl rlocbl
Cio                vcdmel rxes mkrout symrat symrho dfrce pvdf4 pvdf2
Cio                pvdf1 mkekin totfrc mixrho ioden lattic
Cio  s_mix :contains parameters for charge mixing; see structures.h
Ci     Elts read:  b bv w wc nsave mmix umix tolu
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  lmfp
Cio  s_move :struct for dynamics parameters; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  lmfp
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  rho vmtz nlml nlma ves aamom bxc cp ddpf dddpf ddpfr
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
Cio    Elts passed:rhat smrho mad pp pti smpot rnew smrout
Cio    Passed to:  lmfp supot suham rsedit rdovfa iorsf bcast_strx iinit
Cio                mpibc1 chimedit smshft bndfp mkpot suham2 mixrho
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  *
Co     Stored:     nkaps nitab nttab lmaxw nds
Co     Allocated:  iax npr alp s adot sdot
Cio    Elts passed:kaps alph
Cio    Passed to:  lmfp suham rdstrx pp2alp
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  name lmxb z rmt lmxa pz orbp idxdn p norp ntorb rsma
Ci                 kmxt lmxl rg rsmv kmxv lfoca rfoca ncomp ngcut idu
Ci                 uh jh nr a qc rhoc coreh coreq ctail etail stc idmod
Ci                 nxi exi chfa rsmfa rs3 eh3 vmtz mxcst
Co     Stored:     orbp norp ntorb p pz idxdn lmxb ngcut lmxa a nr qc
Co                 nxi exi chfa rsmfa ctail etail stc rhoc name rmt z
Co                 lmxl kmxt lfoca rfoca coreh pb1 pb2
Co     Allocated:  rhoc
Cio    Elts passed:pz name rhoc
Cio    Passed to:  lmfp lmfopb uspecb lmfop2 ioorbp praugm suham atfold
Cio                makidx nscpa showbs sugcut suldau sudmtu praldm
Cio                rotycs symdmu rsedit dfratm chimedit rdovfa bcast_
Cio                strx gtpcor ovlocr corprm adbkql iorsf pvsms2 smshft
Cio                pvsms1 rhgcmp rhogkl bndfp dfaugm mkpot rhomom smves
Cio                vesgcm mshvmt symvvl ugcomp smvxcm smcorm smvxc4
Cio                elocp locpot dfqkkl suham2 suclst surho sumlst
Cio                hambls hambl augmbl bstrux smhsbl hsibq tbhsi
Cio                hsubblock hsibq2 hsibq4 makusq pusq1 mkorbm mullmf
Cio                mkpdos mkdmtu addrbl fsmbl fsmbpw rsibl rsibl1
Cio                rlocbl iorbtm mshn3p mchan vcdmel rxes mkrout symrat
Cio                symrho prrhat pnunew dfrce pvdf4 pvdf2 pvdf1 mkekin
Cio                mixrho ftlxp pvmix5 pvmix3 pvmix7 chkdmu ioden relax
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos relax clabel rho1 rho2 rhoc rho1x rho2x
Ci                 rhocx class force vel pnu pz v0 v1 bxc cpawt omg
Ci                 omgn domg gc gcu gcorr sfvrtx j0 pdos qhhl qhkl qkkl
Ci                 eqhhl eqhkl eqkkl sighh sighk sigkk tauhh tauhk
Ci                 taukk pihh pihk pikk sighhx sighkx sigkkx tauhhx
Ci                 tauhkx taukkx pihhx pihkx pikkx thet pos0
Co     Stored:     pos pnu pz norb vel pos0 force spec clabel bxc cpawt
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
Cio    Passed to:  lmfp rlxstp suham setnorb showbs pvioeu suldau
Cio                sudmtu praldm rotycs symdmu rsedit dfratm chimedit
Cio                rdovfa ovlpfa ovlocr adbkql iorsf pvsms2 bcast_strx
Cio                smshft pvsms1 rhgcmp rhogkl iopos bndfp dfaugm mkpot
Cio                rhomom smves vesgcm mshvmt symvvl ugcomp smvxcm
Cio                smcorm smvxc4 elocp locpot dfqkkl suham2 suclst
Cio                surho sumlst hambls hambl augmbl bstrux smhsbl hsibq
Cio                hsubblock hsibq2 hsibq4 makusq pusq1 mkorbm mullmf
Cio                mkpdos mkdmtu addrbl fsmbl fsmbpw rsibl rsibl1
Cio                rlocbl mshn3p mchan vcdmel rxes mkrout symrat symrho
Cio                prrhat pnunew dfrce pvdf4 pvdf2 pvdf1 mkekin totfrc
Cio                mixrho ftlxp pvmix5 pvmix3 pvmix7 chkdmu ioden cppos
Cio                relax lattic
Cio  s_optic:struct for parameters for optical ME; see structures.h
Ci     Elts read:  ltet ocrng unrng dw window lpart esciss iq
Co     Stored:     ocrng unrng
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  lmfp bndfp optinq optin2 optint
Cio  s_gw  :struct for gw-related parameters
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  lmfp bndfp
Ci Inputs
Ci   mode  :controls what is generated:
Ci         :1s digit is for k-mesh
Ci         :0 qk data supplied by qp,nq.  No BZ integration.
Ci         :This branch requires nq as input; nkabc is not used
Ci         :1 Generate k-mesh internally; BZ integration of some quantities
Ci         :2 Generate k-mesh, but no BZ-integrated properties
Ci         :  Outputs: qibz(3,0:nqfbz) or qibz(3,1:nqfbz),
Ci         :           ipq(0:nqfbz) or ipq(1:nqfbz)
Ci         :3 Same as 1, except special treatment of ef0 if 1000s digit also set
Ci         :Branches 1,2,3 require nkabc as input; nq is not used
Ci         :10s digit (only used if 1s digit is nonzero)
Ci         :1 copy internally generated k-points qpl to qp
Ci         :2 copy internally generated k-point weights wbzl to wbz
Ci         :4 copy internally generated k-mapping array ipql to ipq
Ci         :Bits 1,2,4 may be combined
Ci         :100s digit like 10s digit, but internal qpl,wbzl,ipql are not
Ci               copied but checked against corresponding passed arrays.
Ci         :1 check internally generated k-points array qpl to qp
Ci         :2 check internally generated k-point weights wbzl to wbz
Ci         :4 check internally generated ipql to ipq
Ci         :Bits 1,2,4 may be combined
Ci         :1000s digit.  Evals on k-mesh, either supplied or generated
Ci         :0 No evals generated
Ci         :>0 get evals on k-mesh, store in s_ham%evals.
Ci         :   and if 1s digit has 1s bit set, copy calculated Ef to ef0
Ci         :1 s_ham%evals is shifted to so eval=0 falls at Ef, 
Ci         :2 Compare calculated evals to those given in eig,
Ci         :  and return the mean shift in efshft.
Ci         :  It may often happen that the Fermi level changes because of a change
Ci         :  in k-mesh; efshft may be used to make this adjustment.
Ci         :4 (if noncollinear case) return in eig(:,:,2) projection
Ci         :  of eigenvector onto spin 1
Ci         :  If this digit < 8 s_ham%evals are copied to eig.
Ci         :8 Do not copy or compare to eig; evals are returned in s_ham%evals
Ci   nband :number of QP levels for spectral function; also dimensions eig
Ci   ioff  :Offset to 1st band in QP levels corresponding to those in eig
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   ldq   :dimensions eig,qp,wbz
Ci   nkabc :number of divisions in k-mesh.
Ci         :Used only when 1s digit mode>0; otherwise not used
Ci   nq    :number of k-points to calculate
Ci         :Used only when 1s digit mode=0; otherwise not used
Ci   lunits:0 eig in Ry units
Ci         :1 eig in eV units
Ci   nqoff :Offset for qp-related arrays: poke qp and band info in iq+nqoff
Cio Inputs/Outputs
Cio  ef0   :Fermi level, subtracted from eigenvalues.
Cio        :Not used unless eigenvalues are calculated, in which case:
Cio        :Output if 1s digit mode = 1 and 1000s digit > 0
Cio        :Otherwise, ef0 is input
Cio  qp    :List of qp.
Cio        :used as input qp list if 1s digit mode is 0
Cio        :Output if 1s bit set, 10s digit mode
Cio        :Used for check if 1s bit set, 100s digit mode
Cio        :Otherwise, not used.
Cio  wbz   :qp weights.
Cio        :Output if 2s bit set, 10s digit mode
Cio        :Used for check if 4s bit set, 100s digit mode
Cio        :Otherwise, not used.
Cio  ipq   :mapIBZ->FBZ
Cio        :Output if 4s bit set, 10s digit mode
Cio        :Used for check if 4s bit set, 100s digit mode
Cio        :Otherwise, not used.
Cio  eig   :eigenvalues generated by lmfp
Cio        :Output if 1s bit set, 1000s digit mode
Cio        :Used for check if 2s bit set, 1000s digit mode
Cio  efshft:shift in Fermi level, calculated by lmf compared to GW
Cio        :If 2's bit set in 10000's digit, calculated from
Cio        :mean shift in lmf generated QP levels relative to given values
Cio        :This is only legitimate if the RMS deviation is small
Cio        :(the RMS change is printed out when efshft is calculated)
Cl Local variables
Cl   s_ctrl%plbnd 2 => BZ integration for Ef.  Mesh must be regular mesh
Cl   s_ctrl%plbnd 3 => No BZ integration.  No requirement for mesh.
Cl   nqfbz : number of points in full BZ.  Reduces to nq if 1s digit mode=0
Cb Bugs
Cb   Leading dimension of ipq isn't supplied; check cannot be made
Cb   to ensure that it is at least as large as nqfbz+nqoff
Cu Updates
Cu   02 Jun 13 completed replacement of f77 pointers
Cu             Extended to noncollinear case
Cu   20 Aug 12 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nband,ioff,nsp,ldq,nq,nqoff,lunits,nkabc(3)
      integer ipq(*)
      double precision ef0,eig(nband,ldq,nsp),qp(3,ldq),wbz(ldq)
      double precision efshft
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
C ... Local parameters
      real(8), pointer :: evall(:,:,:)
      real(8), allocatable :: qpl(:,:),wbzl(:)
      real(8), allocatable :: qph(:),wtkph(:)
      integer, allocatable :: ipqh(:),ipql(:)
      integer NULLI,iprint,fopna,LW5
      logical T,F,ltet,lcw,latvec
      double precision rytoev,sshft2,shft,xx(2),ql(3)
      parameter (rytoev=13.6058d0,T=.true.,F=.false.,NULLI=-99999,LW5=5)
      double precision tolq
      parameter (tolq=1d-6)
      integer isp,iq,ib,nglob,stdo,j,nql,mode0,mode1,mode2,mode3,
     .  isw,nqfbz,nkabch(3),nkph,ndham,nspc,ndhamx,nspx,lso,ifiz
      complex(8),allocatable :: zso(:,:,:,:)
      integer, allocatable :: ifblst(:) ! List of orbitals for color weights
      real(8), allocatable :: cwt(:)
      integer nk1,nk2,nk3,lshft1,lshft2,lshft3,lwsig_save

C ... Setup and initial printout
      stdo = nglob('stdo')
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      mode2 = mod(mode/100,10)
      mode3 = mod(mode/1000,10)
      nkabch = 0
      ndham = s_ham%ndham
      nspc = nglob('nspc')    ! 2 for noncollinear case
      ndhamx = ndham*nspc
      nspx = nsp; if (nspc == 2) nspx = 1
      lso = isw(IAND(s_ctrl%lncol,4) /= 0)
      lcw = mode3 >= 4 .and. mode3 /= 8 .and. lso == 1
      if (mode3 < 8) mode3 = mod(mode3,4)
      nullify(evall)
      select case (mode0)
        case (0,2); s_ctrl%plbnd = 3
        case (1,3); s_ctrl%plbnd = 2
        case DEFAULT; call rxi('illegal 1s digit mode',mode0)
      end select
      call info5(10,0,0,' evalqp: Get  '//
     .  '%?#n#k-mesh, ##'//
     .  '%?#n#bands, ##'//
     .  '%?#n#Fermi level, ##%b%b ...',isw(mode0 /= 0),
     .  isw(mode3 /= 0),isw(mode3 /= 0.and.s_ctrl%plbnd == 2),0,0)

C --- Generate or copy k-points ---
C     Hang onto these pointers so they can be restored on exit
      j = size(s_bz%ipq)
      allocate(ipqh(j)); call icopy(j,s_bz%ipq,1,ipqh,1)
      j = size(s_bz%qp); allocate(qph(j)); call dcopy(j,s_bz%qp,1,qph,1)
      j = size(s_bz%wtkp)
      allocate(wtkph(j)); call dcopy(j,s_bz%wtkp,1,wtkph,1)

      nkph = s_bz%nkp
      if (mode0 == 0) then
        call ptr_bz(s_bz,1,'qp',3,nq,xx)
        call dcopy(3*nq,qp,1,s_bz%qp,1)
        s_bz%nkp  = nq
        nql = nq
        nqfbz = nq
      else
        nkabch = s_bz%nkabc
        s_bz%nkabc = nkabc
        ltet = IAND(s_ctrl%lmet,3) == 3 .or.
     .         IAND(s_ctrl%ldos,4+2+1) /= 0
        call mkqp(s_ctrl,s_bz,s_lat,ltet,F,1,-2)
        nql = s_bz%nkp
        nqfbz = nkabc(1)*nkabc(2)*nkabc(3)
      endif

C --- Generate evals ---
      if (mode3 /= 0) then
        allocate(s_ham%evals(ndham,nsp,nql))
        call dpzero(s_ham%evals,ndham*nsp*nql)
        lwsig_save = s_ctrl%lwsig
        if (lcw) s_ctrl%lwsig = LW5
        call lmfp('evalqp:',s_bz,s_ctrl,s_ham,s_lat,s_mix,
     .    s_move,s_pot,s_str,s_spec,s_site,s_optic,s_gw,s_dmft,s_strn)
        s_ctrl%lwsig = lwsig_save
        shft = ef0  ! Save input ef0
        if (s_ctrl%plbnd == 2) then
          ef0 = s_bz%ef
          if (ef0 == 999) then
            ef0 = NULLI
          else
            if (s_bz%egap /= NULLI) then
              ef0 = s_bz%ef + s_bz%egap/2
            endif
          endif
        endif
        if (s_ctrl%plbnd == 2 .and. shft /= NULLI) then
          call info5(30,1,-1,' evalqp:  Ef=%d (calc)  Ef=%d (given)'//
     .      '  shift=%d ... Use %?#(n==1)#calculated#given# value',
     .      ef0,shft,ef0-shft,mode0,0)
          if (mode0 == 3) ef0 = shft
        endif
        call info5(30,1,0,' evalqp:  Ef=%d Ry'//
     .    '%?#n# (incl Egap/2=%d)#%j#.'//
     .    '  Shift bands to Ef=0'//
     .    '%?#n# and scale to eV ##...',ef0,
     .    isw(s_bz%egap /= NULLI),s_bz%egap/2,lunits,0)
        call dvadd(s_ham%evals,1,ndham*nsp*nql,-ef0)
        if (lunits == 1) then
          call dscal(ndham*nsp*nql,rytoev,s_ham%evals,1)
        endif
      endif

C --- Copy or check generated BZ arrays to passed arrays ---
      allocate(qpl(3,nql),wbzl(nql),ipql(nqfbz))
C     Straight copy
      call dcopy(nql*3,s_bz%qp,1,qpl,1)
      if (mod(mode1/2,2) /= 0 .or. mod(mode2/2,2) /= 0) then
        call dcopy(nql,s_bz%wtkp,1,wbzl,1)
      endif
      if (mode0 == 0) then
        do  iq = 1, nq
          ipql(iq) = iq
        enddo
      else
        call icopy(nqfbz,s_bz%ipq,1,ipql,1)
        if (nqoff > 0) then
          do  iq = 1, nqfbz
            ipql(iq) = ipql(iq)+nqoff
          enddo
        endif
      endif

      if (mode1 /= 0 .or. mode2 /= 0) then
        if (nql+nqoff > ldq)
     .    call rxi('ldq is dimensioned too small, need',nql+nqoff)
        if (mod(mode1,2) /= 0) then
          call dcopy(nql*3,qpl,1,qp(1,1+nqoff),1)
        elseif (mod(mode2,2) /= 0) then
          call daxpy(nql*3,-1d0,qp(1,1+nqoff),1,qpl,1)
          if (max(maxval(qpl),-minval(qpl)) > tolq)
     .      call rx('generated and given qp are not equal')
        endif
        if (mod(mode1/2,2) /= 0) then
          call dcopy(nql,wbzl,1,wbz(1+nqoff),1)
        elseif (mod(mode2/2,2) /= 0) then
          call daxpy(nql,-1d0,wbz(1+nqoff),1,wbzl,1)
          if (max(maxval(wbzl),-minval(wbzl)) > tolq)
     .      call rx('generated and given wbz are not equal')
        endif
        if (mod(mode1/4,2) /= 0) then
          call icopy(nqfbz,ipql,1,ipq(1+nqoff),1)
        elseif (mod(mode2/4,2) /= 0) then
          call iaxpy(nqfbz,-1,ipq,1,ipql,1)
          if (max(maxval(ipql),-minval(ipql)) > 0)
     .      call rx('generated and given ipq are not equal')
        endif
      endif

C --- Compare evals to eig and/or copy evals to eig ---
      if (mode3 /= 0 .and. mode3 < 8) then
        if (.not. associated(evall)) evall => s_ham%evals
        if (ef0 == NULLI)
     .    call rx('evalqp: cannot assign evals w/out Fermi level')
        call info0(43,0,0,'   ib   iq   isp%5forig%7frecalc%6fshift')
        if (mode3 >= 2) then
          efshft = 0
          sshft2 = 0
        endif

C       Setup for suqlsz call to get spin weights in SO coupled case
        if (lcw) then
          allocate(zso(ndham,2,ndham,2))
          allocate(ifblst(ndham),cwt(ndhamx))
          do  j = 1, ndham
            ifblst(j) = j
          enddo
          ifiz = fopna('evec',-1,4)
          rewind ifiz
          call iosigh(2,LW5,nsp,nspc,j,j,nk1,nk2,nk3,ib,iq,lshft1,lshft2,lshft3,ifiz,xx)
          if (j /= ndham) call rx('evalqp: evec file mismatch, hamiltonian dimension')
          if (ib /= nq) call rx('evalqp: evec file mismatch, no k-points')
C         Make copy of qp list actually used
          if (mode0 == 0) then
            call dcopy(3*nq,qp,1,qpl,1)
          else
            call dcopy(nql*3,s_bz%qp,1,qpl,1)
          endif
        endif

        do  isp = 1, nspx
          do  iq = 1, nql

            if (lcw) then  ! Projection of eigenfunction onto spin 1
              read(ifiz) ql
              if (.not. latvec(1,tolq,s_lat%qlat,ql-qpl(1:3,iq))) then
                write(stdo,*) iq
                write(stdo,*) sngl(ql)
                write(stdo,*) sngl(qp(:,iq))
                call rx('incompatible q-mesh')
              endif
              call dpdump(xx,1,ifiz) ! Skip over evals
              call dpdump(zso,ndhamx**2*2,ifiz)
              call suqlsz(ndhamx,1,ndham,ifblst,ndham,zso,cwt)
            endif

            do  ib = 1, nband
              if (mode3 >= 2) then
                shft = eig(ib,iq+nqoff,isp) - evall(ib+ioff,isp,iq)
                efshft = efshft + shft
                sshft2 = sshft2 + shft*shft
                if (iprint() >= 43) then
                  write(stdo,"(3i5,3f12.6)") ib,iq,isp,
     .              eig(ib,iq+nqoff,isp),evall(ib+ioff,isp,iq),shft
                endif
              endif
              if (mod(mode3,2) /= 0 .and. mod(mode3,2) < 8) then
                eig(ib,iq+nqoff,isp) = evall(ib+ioff,isp,iq)
                if (lcw) then
                  eig(ib,iq+nqoff,2) = cwt(ib+ioff)
                endif
              endif
            enddo
          enddo
        enddo
        if (mode3 >= 2 .and. mod(mode3,2) < 8) then
          j = nband*nql*nsp; sshft2 = sshft2/j; efshft = efshft/j
          sshft2 = dsqrt(sshft2-efshft**2)
          call info2(20,0,0,'%10fAlign given evals with calculated: Average shift = %d  '//
     .      'RMS deviation = %d',efshft,sshft2)
        endif
        if (.not. associated(evall,s_ham%evals)) deallocate(evall)
        deallocate(s_ham%evals)
      endif
      if (lcw) then
        call fclose(ifiz)
      endif
      deallocate(qpl,wbzl,ipql)

C ... Restore s_bz parameters
      j = size(ipqh); call ptr_bz(s_bz,1,'ipq',j,0,xx)
      call icopy(j,ipqh,1,s_bz%ipq,1); deallocate(ipqh)

      j = size(qph); call ptr_bz(s_bz,1,'qp',3,j/3,xx)
      call dcopy(j,qph,1,s_bz%qp,1); deallocate(qph)

      j = size(wtkph); call ptr_bz(s_bz,1,'wtkp',j,0,xx)
      call dcopy(j,wtkph,1,s_bz%wtkp,1); deallocate(wtkph)

      s_bz%nkabc = nkabch
      s_bz%nkp = nkph

      end
