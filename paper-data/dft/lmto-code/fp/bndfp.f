      subroutine bndfp(s_strn,s_site,s_spec,s_lat,s_ctrl,s_ham,s_pot,
     .  s_bz,s_optic,s_gw,nbas,nsp,nlibu,lmaxu,lldau,leks,lrout,lfrce,
     .  lpnu,dmxp,iter,maxit,frc,dmatu,vorb)
C- One band pass, full-potential hamiltonian
C ----------------------------------------------------------------------
Cio Structures
C    ... lmf branch
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  pos spec rho1 rho2 rhoc rho1x rho2x rhocx class pnu
Ci                 pz v0 v1 sighh sighk sigkk pihh pihk pikk sighhx
Ci                 sighkx sigkkx pihhx pihkx pikkx sohh sohk sokk qkkl
Ci                 qhkl qhhl eqkkl eqhkl eqhhl tauhh tauhk taukk
Co     Stored:     saxis pnu pz force
Co     Allocated:  rho1 rho2 rhoc rho1x rho2x rhocx sigkk taukk sigkkx
Co                 taukkx pikk sokk pikkx sighk tauhk sighkx tauhkx
Co                 pihk sohk pihkx sighh tauhh sighhx tauhhx pihh sohh
Co                 pihhx qkkl eqkkl qhkl eqhkl qhhl eqhhl
Cio    Elts passed: rho1 rho2 rhoc rho1x rho2x rhocx v0 qkkl qhkl qhhl
Cio                eqkkl eqhkl eqhhl pihhx tauhh pikkx taukk pihkx tauhk
Cio    Passed to:  dfratm grp2av pgrp2a prrhat rhopos siterhopos dfaugm
Cio                lkdim mkpot rhomom smves vesgcm mshvmt symvvl ugcomp
Cio                smvxt smvxcm smcorm smvxc4 elocp locpot msh21c
Cio                suham2 dfqkkl suclst surho sumlst hambls hambl
Cio                augmbl bstrux smhsbl hsibq hsubblock hsibq2 hsibq4
Cio                mkhso makusq pusq1 mkorbm mullmf mkpdos mkdmtu
Cio                addrbl fsmbl fsmbpw rsibl rsibl1 rlocbl sosite
Cio                psymrq1 mshn3p mchan vcdmel rxes mkrout symrat
Cio                symrho pnunew dfrce pvdf4 pvdf2 pvdf1 mkekin totfrc
Cio                mixrho rhgcmp rhogkl ftlxp pvmix5 pvmix3 pvmix7
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  nr lmxl z rmt a grp2 lmxa lmxb kmxt qc rg lfoca
Ci                 rfoca ctail etail stc p pz name rs3 eh3 vmtz orbp
Ci                 rsma idu mxcst coreh coreq ngcut idxdn ncomp idmod
Ci                 nxi exi chfa rsmfa rsmv kmxv
Co     Stored:     orbp ngcut idxdn
Co     Allocated:  *
Cio    Elts passed:z rhoc
Cio    Passed to:  dfratm grp2av pgrp2a prrhat rhopos siterhopos dfaugm
Cio                lkdim mkpot rhomom corprm smves vesgcm mshvmt symvvl
Cio                ugcomp smvxt smvxcm smcorm smvxc4 elocp uspecb
Cio                locpot gtpcor msh21c suham2 sugcut rdsigm siged
Cio                makidx nscpa dfqkkl suclst surho sumlst hambls hambl
Cio                augmbl bstrux smhsbl hsibq tbhsi hsubblock hsibq2
Cio                hsibq4 mkhso makusq pusq1 mkorbm mullmf mkpdos fpopm
Cio                mkdmtu addrbl fsmbl fsmbpw rsibl rsibl1 rlocbl
Cio                sosite psymrq1 iorbtm mshn3p mchan vcdmel rxes
Cio                mkrout symrat symrho pnunew dfrce pvdf4 pvdf2 pvdf1
Cio                mkekin mixrho rhgcmp rhogkl ftlxp pvmix5 pvmix3
Cio                pvmix7
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat qlat vol nabc afmt nsgrp npgrp ng pos kv
Ci                 gv awald tol nkd nkq tolft kv2 igv igv2 napw gmax
Co     Stored:     napw igv igv2 ng gmax nabc
Co     Allocated:  gv kv igv kv2 igv2
Cio    Elts passed: igv2 plat pos istab symgr ag nsgrp vol nabc ng kv
Cio                gv ips0 bgv cy cg indxcg jcg qlv dlv alat qlat igv
Cio                napw
Cio    Passed to:  rhopos mkpot ioden2 symsmr smves vesft vesgcm mshvmt
Cio                symvvl ugcomp ggugbl gfigbl fklbl gklbl hhugbl
Cio                hhigbl phhigb hklbl hsmbl hgugbl smvxt smvxc2 vxcnlm
Cio                smvxcm smcorm locpot msh21c suham2 rdsigm siged
Cio                suqlst sugvec hambls hambl augmbl bstrux hxpbl ghibl
Cio                hxpgbl ghigbl hklgbl smhsbl hhibl phhibl hsibq mkhso
Cio                makusq pusq1 fpopm addrbl fsmbl rsibl rlocbl sosite
Cio                sugvec0 vcdmel rxes mkrout symrat symrho dfrce pvdf4
Cio                pvdf2 pvdf1 mkekin totfrc mixrho rhgcmp
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  lekkl ldos lmet zbak plbnd lwsig nl lfp ips lham
Ci                 loptc lfrce
Co     Stored:     lwsig
Co     Allocated:  *
Cio    Elts passed:lncol lbas lcd lbxc maxmem lwsig lfp ldos spid ips
Cio    Passed to:  grp2av suham2 rdsigm siged fpopm optinq vcdmel rxes
Cio                dfrce
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  ldham pwmode pwemin pwemax nbf lsig oveps ovncut
Ci                 ndham lncol elind evals nqsig eseavr sigp pmin pmax
Ci                 pnudef eterms eula lrsa rsrnge rsstol nprs qsig offH
Ci                 nlmto ndhrs
Co     Stored:     lsig eterms lncol elind sigp ehf ehk eula ndhrs
Co                 eseavr nqsig iaxs
Co     Allocated:  nprs iaxs hrs qsig
Cio    Elts passed: lncol lrsa magf nlmto qsig iprmb offH evals eula
Cio                nprs iaxs hrs rsrnge
Cio    Passed to:  mkpot locpot suham2 rdsigm hft2rs siged hambls hambl
Cio                augmbl sopert3 mkhso blsig makusq addrbl sosite
Cio                mkrout mkehkf mkekin
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  nlml nlma hab sab ppn smrout lcplxp vesrmt shfac
Ci                 bxcscali socscl rnew rhat
Co     Stored:     lcplxp qval
Co     Allocated:  hab sab smrout
Cio    Elts passed: vesrmt rhat smrho ppn lcplxp smpot smvcnst rnew
Cio                smrout hab sab smvextc bxcscali v0beta nlml
Cio    Passed to:  mkpot smvxt locpot suham2 mixrho
C
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  lmet nkabc nkp ntet n w nevmx efmax fsmom ndos dosw
Ci                 ef def numq lshft egap star nef qp wtkp sopertp range
Co     Stored:     nef egap numq wtkp wtkb ndos dosw ef def sopertp
Co     Allocated:  qp wtkb swtk
Cio    Elts passed: lio sopertp qp nef numq swtk wtkp wtkb lshft nkabc
Cio                ipq idtet egap n w def
Cio    Passed to:  rdsigm subzi addrbl sosite optinq optin2 optint
Cio  s_optic:struct for parameters for optical ME; see structures.h
Ci     Elts read:  loptic kmxax ltet ocrng unrng dw window mefac nlg
Ci                 alltrans optme lpart esciss iq ffmt imref
Co     Stored:     nlg ocrng unrng nfilm nempm optme imref
Co     Allocated:  *
Cio    Elts passed: rgrad rgrade kmxax esmr loptic ocrng unrng optme
Cio                imref kt
Cio    Passed to:  mkpot locpot ocnock opt_nstate fpopm addrbl rsibl
Cio                optinq optin2 optint
Cio  s_gw  :struct for gw-related parameters (not used)
C      ... GW driver branch
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec rho1 rho2 rhoc rho1x rho2x rhocx pos class pnu
Ci                 pz v0 v1 clabel sighh sighk sigkk pihh pihk pikk
Ci                 sighhx sighkx sigkkx pihhx pihkx pikkx
Co     Stored:     *
Co     Allocated:  rho1 rho2 rhoc rho1x rho2x rhocx sigkk taukk sigkkx
Co                 taukkx pikk pikkx sighk tauhk sighkx tauhkx pihk
Co                 pihkx sighh tauhh sighhx tauhhx pihh pihhx qkkl
Co                 eqkkl qhkl eqhkl qhhl eqhhl
Cio    Elts passed:rho1 rho2 rhoc rho1x rho2x rhocx v0
Cio    Passed to:  dfratm dfaugm mkpot rhomom smves vesgcm mshvmt
Cio                symvvl ugcomp smvxcm smcorm smvxc4 elocp locpot
Cio                dfqkkl suham2 suclst chkgwin sugwin sugw wlattc mk_
Cio                hamindex hambl augmbl bstrux smhsbl hsibq hsubblock
Cio                hsibq2 hsibq4 hambls makusq pusq1 gwcphi pwmat pwmat2
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  nr lmxl lmxa lmxb kmxt z qc a rmt rg lfoca rfoca
Ci                 ctail etail stc p pz name rs3 eh3 vmtz orbp rsma idu
Ci                 mxcst coreh coreq ngcut idxdn pb1 pb2
Co     Stored:     orbp ngcut
Co     Allocated:  *
Cio    Elts passed:name
Cio    Passed to:  dfratm dfaugm mkpot rhomom corprm smves vesgcm
Cio                mshvmt symvvl ugcomp smvxcm smcorm smvxc4 elocp
Cio                uspecb locpot gtpcor dfqkkl suham2 sugcut suclst
Cio                chkgwin sugwin sugw wlattc mk_hamindex hambl augmbl
Cio                bstrux smhsbl hsibq tbhsi hsubblock hsibq2 hsibq4
Cio                hambls makusq pusq1 gwcphi pwmat pwmat2
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat qlat nabc afmt vol ng nsgrp awald tol nkd
Ci                 nkq tolft npgrp igv2 kv2 kv igv gmax napw
Co     Stored:     igv2 kv2 igv napw
Co     Allocated:  gv kv igv kv2 igv2
Cio    Elts passed: gv kv cy symgr ag cg indxcg jcg qlv dlv pos istab
Cio                igv2 kv2 qlat igv s_sym vol plat
Cio    Passed to:  mkpot smves vesft vesgcm mshvmt symvvl ugcomp ggugbl
Cio                gfigbl fklbl gklbl hhugbl hhigbl phhigb hklbl hsmbl
Cio                hgugbl smvxc2 vxcnlm smvxcm smcorm locpot rdsigm
Cio                suham2 chkgwin sugwin sugvec sugw mk_hamindex hambl
Cio                augmbl bstrux hxpbl ghibl hxpgbl ghigbl hklgbl
Cio                smhsbl hhibl phhibl hsibq hambls makusq pusq1 pwmat
Cio                pwmat2
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  loptc pfloat ldos zbak plbnd lfp
Co     Stored:     lfp
Co     Allocated:  *
Cio    Elts passed:lncol lbas lfp lcd
Cio    Passed to:  suham2 sugw
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  ldham pwmode pwemin pwemax nbf lsig oveps ovncut
Ci                 ndham rsrnge elind evals eterms sigp rsstol nprs
Ci                 qsig nqsig ndhrs eseavr
Co     Stored:     lsig eterms nqsig ndhrs eseavr
Co     Allocated:  nprs hrs qsig iaxs
Cio    Elts passed:magf lncol iprmb qsig offH hrs iaxs nprs
Cio    Passed to:  mkpot rdsigm hft2rs suham2 sugwin sugw mk_hamindex
Cio                hambl hambls sopert3 makusq
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  nlml nlma smrout
Co     Stored:     *
Co     Allocated:  smrout
Cio    Elts passed:rhat smrho smpot
Cio    Passed to:  mkpot suham2
C
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  nkabc nkp ntet lmet n w nevmx efmax fsmom ndos dosw
Ci                 ef def numq lshft qp
Co     Stored:     numq
Co     Allocated:  wtkb swtk
Cio    Elts passed:lio numq n wtkb w def
Cio    Passed to:  rdsigm subzi chkgwin sugwin
Cio  s_optic:struct for parameters for optical ME; see structures.h
Ci     Elts read:  ltet ocrng unrng dw window
Co     Stored:     ocrng unrng
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_gw  :struct for gw-related parameters
Ci     Elts read:  code mksig lgw nkabc nband qoffp gcutb gcutx ecuts
Ci                 nime delre deltax deltaw pbtol gsmear
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  chkgwin sugwin sugw
Ci Inputs
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nlibu :number of U blocks (LDA+U) (used to dimension dmatu and vorb)
Ci   lmaxu :dimensioning parameter for U matrix (used to dimension dmatu and vorb)
Ci   lldau :lldau(ib)=0 => no U on this site otherwise
Ci          U on site ib with dmat beginning at dmats(*,lldau(ib))
Ci   leks  :>0 make the Hohnberg-Kohn-Sham energy
Ci          >1 use the HKS forces
Ci   lrout :>0 generate output density and attendant quantities
Ci         :Note: if bndfp operates in a special-purpose mode
Ci         :(s_ctrl%plbnd > 0 or (lwsig > 0 .and. lwsig <= LW4), bndfp sets lrout=0
Ci   lfrce : 0 suppress generation of forces
Ci   lpnu  : 1 make new pnu's
Ci   dmxp  :vector of mixing parameters; see mixrho.f for dmxp(1..25)
Ci         :Additionally:
Ci         :dmxp(26)  is the iteration where the minimum RMS error was found
Ci         :dmxp(27)  is the minimum RMS error.
Ci         :dmxp(33)  is the Lindhard parameter
Ci   iter  :current iteration number
Ci   maxit :maximum number of iterations (for printout only)
Cio LDA+U inputs and outputs
Cio  dmatu :density matrix for LDA+U (changed upon output)
Cio  vorb  :orbital dependent LDA+U potential
Co Outputs
Co   frc   :forces
Cs Command-line switches
Cs   --asars       : Write ASA restart file for input potential
Cs   --asars2      : Write ASA restart file for output potential, including charges
Cs   --band        : Tabulate energy bands; see doc/Command-line-options.html
Cs   --bonly       : Use in conjuction with --lxb switch
Cs   --chklim      : Check band limits in optics calculations
Cs   --cls         : For core-level spectroscopy
Cs   --cv:         : Calculate electronic specific heat, eV
Cs   --cvK:        : Calculate electronic specific heat, Ry
Cs   --dos[options]: Calculate total DOS, and qualifying switches
Cs   --ef=         : Override file Fermi level; use with --band
Cs   --efrnge      : Print out indices to bands that bracket the Fermi level
Cs   --jdosw       : Channels for optical properties; See doc/optics.html
Cs   --jdosw2      : Channels for optical properties; See doc/optics.html
Cs   --lxb         : Check band limits in optics calculations
Cs   --mixsig=     : For self-energy; see Command-line-options.html
Cs   --mlog        : (MPI) write MPI commands to log file
Cs   --mull        : Mulliken analysis; see doc/Command-line-options.html
Cs   --no-fixef0   : Do not adjust estimate of Fermi level after 1st band pass
Cs   --oldbz       : Not documented
Cs   --oldvc       : Reference potential defined so average cell potential is zero
Cs   --onesp       : Generate bands for one spin.  Better: use --band~spin1...
Cs   --opt:read    : Read optical matrix elements; See doc/optics.html
Cs   --opt:write   : Write square of optical matrix elements; See doc/optics.html
Cs   --opt:woptmc  : Write complex optical matrix elements
Cs   --optbas      : Special branch to optimise basis
Cs   --pdiag       : (MPI) parallel diagonaliser
Cs   --pdos        : Partial DOS; see doc/Command-line-options.html
Cs   --quit=       : quit after execution of certain blocks, e.g. --quit=ham
Cs   --rsig        : For reading the self-energy; see Command-line-options.html
Cs   --rxes        : RXES spectroscopy
Cs   --shorten=    : Suppress all shortening of basis vectors
Cs   --symsig      : Symmetrize sigma (if read) overriding default behavior
Cs   --SOefield    : Makes SO-weighted average of electric field at each site
Cs   --vext        : Add external potential
Cs   --wden[..]    : Generate output density to a file, in various ways
Cs   --window=     : Generate output density in a restricted energy window
Cs   --wpotmt      : Write augmentation potentials to atom files
Cs   --wrhoat      : Write partial atomic densities to atom files
Cs   --wrhomt      : Write augmentation densities to files
Cs   --wsig        : For writing the self-energy; see Command-line-options.html
Cs   --wsmpot      : Write smooth potential to file
Cs   --wsmrho      : Write smooth density to file
Cs  GW-specific
Cs   --make-Q0P  :Create Q0P file
Cs   --no-GWinput:Do not create GWinput file
Cs   --novxc     :Do not write matrix elements of vxc to disk
Cs   --sc        :Tell setup to set flags for QSGW
Cs   --vxcsig    :overwrite vxc with vxc + (sigm-vxc) (1-shot on QSGW)
Cs   --jobgw=    :set GW job through this switch (if missing, job read from stdin)
Cl Local variables
Ci   evlk  :band eigenvalues for current spin, qp
Cl   k1,k2,k3: dimensions smrho,smpot
Cl   lpdiag:0 use standard diagonalization (zhev)
Cl         :1 use parallel diagonalization (pzhev)
Cl         :2 diagonalization done internally (hambls)
Cl   lwndow:T if to make density in a specified energy window
Cl   jsp   :current spin index.
Cl         :In the collinear case, jsp and isp are equivalent
Cl         :In the noncollinear case, isp loops 1..2 for the
Cl         :purpose of assembling the hamiltonian.
Cl         :Once assembled, isp should not be used; and jsp=1
Cl   ispc  :2 when working on (2,2) block of noncollinear hamiltonian;
Cl         :otherwise 1
Cl   nspx: number of independent spin channels containing eigenvalues
Cl         and total DOS; nspx=nsp unless nspc=2, in which case nspx=1
Cl   onesp :do only one spin branch of isp loop (spec'd by onesp)
Cl         :also used when usual loop order (iq=1..nq, isp=1..2)
Cl         :needs to be reversed, as it does, e.g. when transforming
Cl         :sigma matrix.  Then onesp plays the role of spin index
Cl   llmet :0 nonmetal
Cl         :1 metal
Cl         :2 metal, efermi specified externally with --ef=#
Cl   lswtk :Flags whether to make 'spin weights' swtk
Cl         :-2 do not make spin weights
Ci         :1 given a set of weights, make 'spin weights' swtk
Cl   lwtkb :0 weights are neither required nor available a priori
Cl         :1 weights are required a priori, and are read from disk
Cl         :-1 weights are required a priori, but were not read
Cl         :2 weights were generated with constrained global moment
Cl   lekkl :0 do not accumulate eqkkl; 1 do accumulate eqkkl
Cl   lwsig :special modes to handling reading/writing of sigma or evecs
Cl         :LW1  Rotates sigm to LDA basis; saves in file 'sigm2'.
Cl         :LW2  Similar to lwsig=1, except
Cl         :     low- and high- energy blocks replaced by diagonal parts
Cl         :LW3  Proceeds like LW2, but sigm(orbital) -> sigm(orbital) before writing
Cl         :LW4  Writes evals,evecs of LDA hamiltonian to file 'evec'
Cl         :LW5  Writes evals,evecs of hamiltonian to file 'evec'
Cl         :LW6  Similar to LW5, but occupation numbers are also saved.
Cl         :LW7  Writes evals,evecs of hamiltonian to file 'eveca' in ASCII format
Cl         :LW8  Reads evals,evecs of hamiltonian from file 'evec' (compatibility with lmfgwd)
Cl  ldmatk :1 make and store density n_ij(k), eigenfunction basis
Cl   ltso  :1 Make site-resolved SO contribution to band energy
Cl   ndham :dimensioning parameter, at least as large as largest
Cl         :hamiltonian dimension
Cl   ndimh :Hamiltonian dimension for one spin and current k.
Cl         :In the APW case with pwmode>10 it may be less than ndham
Cl   nsmidb:smallest value of nmax encountered in truncating sigma
Cl         :   (only used for printout)
Cl   nevl  :By default, nevl = hamiltonian dimension.
Cl         :If the hamiltonian is mapped to a reduced hilbert space,
Cl         :nevl is the dimension of the reduced space (see ham->oveps)
Cl   nfbn  :number of channels in decomposition of DOS like objects
Cl         :(1) = number of channels, 1st decomposition
Cl         :(2) = number of channels, 2nd decomposition
Cl         :(3) = number of unocc channels, 1st decomp (joint DOS)
Cl         :(4) = number of unocc channels, 2nd decomp (joint DOS)
Cl   nphimx:max number of partial waves of a given l channel at any site
Cl   plbnd :Controls what FP code lmf makes in band pass (rdctrl2,evalqp,sugws)
Cl          0 band pass, integrates over BZ for density, total energy
Cl          1 band pass makes bands at qp tabulated in file, written to disk
Cl          2 like 0, but no output density.  Bands stored in s_ham%evals
Cl          3 like 2, but no attempt to integrate over BZ
Cl            Caller must allocate s_ham%evals
Cl         10 make potential and exit
Cl          plbnd>0 suppresses parsing of --dos switch
Cl   ldsdk :T => calculate dSigma/dk from QSGW for optics velocity operator
Cl   lmefac:Used for ME correction when nonlocal sigma added to LDA pot.
Cl          0 default operation
Cl          1 Flag for special LDA band pass before optics
Cl         -1 Regular pass after special LDA band pass
Cl   lwden :T => Makes density for only 1 k-point
Cr Remarks
Cr   Band pass consists of:
Cr   (1) make the effective potential,
Cr   (2) generate eigenvalues (and eigenvectors if lrout>0)
Cr   (3) if lrout>0, assemble the output density by BZ integration
Cr   (4) evaluate hf (and KS, if leks) energy by BZ integration
Cr   (5) mix the output density to make a new input density.
Cr
Cr  Starting with version 7.9b, q-vectors need not be shortened.
Cr  Algorithms were redesigned so that q-vectors can be shortened
Cr  internally where appropriate, e.g. in hambls.f, sugvec.f, makusq.f
Cr
Cr  Algorithms were constructed with the property that the result be
Cr  independent of shortening.  This introduces some complication
Cr  when the APW basis is used: q-independent and q-dependent APWs
Cr  must be treated slightly differently; see e.g. hambls.f, and
Cr  the code prior to pwmat in sugw.f
Cu Updates
Cu   09 Sep 18 Bug fix switch --dos:ef0...
Cu   05 Apr 18 s_ctrl%plbnd=10 causes bndfp to return after making potential
Cu   18 Mar 18 Enable reading of evals, evecs from contents of evec file (--revec:gw)
Cu             enable k-reordering of optical matrix element to follow qpts file (--opt:write:permqp)
Cu             Useful for exact compatibility with, e.g. lmfgwd
Cu   15 Mar 18 New lmfgwd compatibility mode (LW8)
Cu   28 Oct 17 bndfp can write asa-style rsta files (--asars or --asars2)
Cu   07 Aug 17 Revise AFM symmetrization
Cu   28 Jul 17 --pdos, --mull and band plotting with color weights works with MPIK
Cu   28 Jul 17 Extend --pdos to the noncollinear case
Cu   14 May 17 LDA+U can work with AFM symmetry
Cu   26 Mar 17 New switches --opt:woptmc and --opt:woptmc,rdqp
Cu      Jan 17 (B. Cunningham) write out complex optical matrix elements
Cu   28 Nov 16 Print out lowest, highest level for each band if ipr>50
Cu   10 Aug 16 Some adjustments enabling sudmft to generate output density
Cu   23 Apr 16 --wden~qp=~... saves density for a particular q
Cu   29 Mar 16 (P. Azarhoosh) restructure for reduced dimensions in evnl, evloc when doing optics
Cu   29 Oct 15 Added --dos switch
Cu   10 Oct 14 Some changes to adapt to the DMFT interface
Cu   21 Aug 14 Optics reworked to reflect the proper three-fold representation
Cu             of optical matrix elements
Cu   04 Jun 14 (P. Azarhoosh) optical ME correction factor with NL potential
Cu             MvS made to work for spin pol case
Cu   15 Nov 13 (Ben Kaube) first cut at optical matrix elements
Cu             (MvS) first cut at the noncollinear case
Cu   11 Oct 13 New --ef switch
Cu   09 Aug 13 Option to calculate bands, sumev by tetrahedron,
Cu             weights for density by sampling with lmet modes 2,3
Cu             Option to resolve SO contribution to energy by site
Cu             in 2nd order perturbation theory
Cu   02 Aug 13 New averaging l=0 rho by group, and rhopos
Cu   30 Jun 13 (WRL) bndfp writes evecs to ascii file if s_ctrl%lwsig=7
Cu   01 Jun 13 Implement lwsig=5 properly
Cu   08 May 13 Eliminate s_array
Cu   12 Apr 13 First cut at generating density matrix
Cu   10 Feb 13 New ability to shorten qp in local routines
Cu             to avoid global shortening.  For compatibility with new GW.
Cu             To test, set default for OPTIONS_SHORBZ to F.
Cu             All tests should be unaffected.
Cu   12 Nov 12 migrated structures to f90 pointers
Cu   09 Aug 12 New SO=3 option; new s_ctrl%plbnd=3
Cu   25 Oct 11 Started migration to f90 structures
Cu   18 Mar 11 New --quit=dos --- when only total DOS sought
Cu   11 Mar 11 Reworked --rhopos switch to force local rho > 0
Cu   19 Jul 10 Reworked fixed-spin-moment method
Cu   01 Jun 10 Option to render blocks of hamiltonian diagonal
Cu   01 Apr 10 Add external B field to potential.  New argument list
Cu   08 Mar 10 Bugs fixes for l-dependent cutoffs in Mulliken analysis;
Cu             pdos and Mulliken analysis on a consistent footing.
Cu   08 Feb 10 (D. Pashov) new METAL=5, a faster version of METAL=3, in F90
Cu   30 Jan 10 Handles epsovl with self-energy
Cu   05 Jan 10 New AFM symmetry
Cu   19 Jun 09 (W. Lambrecht) Calculate resonant x-ray emission spectra
Cu   04 Sep 09 First cut at optics: Make joint DOS
Cu   02 Jun 09 some fixes for MPI, v7; new call to pnunew
Cu   05 Jul 08 (T. Kotani) new PW basis
Cu             Option to accumulate energy-weighted output density
Cu   27 Jun 08 Redesigned transformation of sigma to new basis
Cu   09 Jul 07 MPIK enabled to plot bands
Cu   05 Jul 07 Enable onesp to be set as switch in --band:spin1
Cu   09 Jun 07 Fixed-spin-moment, noncollinear case
Cu   16 Jan 07 First cut at I/O of sigm transformed to,from LDA basis
Cu   26 May 07 Some preparation for rotation betw/ LDA, GW basis
Cu   17 Jul 06 Some MPI changes, fixes SO case
Cu   27 Jun 06 New constraints for floating pnu
Cu   08 Jun 06 Bug fix total DOS, noncollinear case;
Cu             Mulliken works for noncollinear case
Cu   02 Jan 06 better printout of magnetic moment
Cu   27 Nov 05 LDA+U => complex potential
Cu   09 Nov 05 (wrl) Convert dmat to complex form
Cu   06 Oct 05 (A. Chantis) bug fix dos when nspc=2
Cu   25 Jul 05 bug fix partial dos combined with LDA+U
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   29 Jun 05 (MvS) extended LDA+U to local orbitals
Cu   27 Apr 05 LDA+U (Lambrecht)
Cu   14 Feb 05 fixes for band plot, contour mode
Cu   03 Feb 05 (A. Chantis) implemented spin-orbit coupling by L.S
Cu   11 Jan 05 double-counting term rho*sig subtracted from ehks.
Cu   23 Dec 04 Extended to spin-coupled case
Cu   18 Nov 04 Sampling integration properly handles Fermi distribtion
Cu   25 Sep 04 (ATP) some patches for MPI parallelization
Cu    1 Sep 04 Adapted to handle complex ppi. Spin-orbit folded into ppi
Cu   12 Aug 04 First implementation of extended local orbitals
Cu   29 Jun 04 (A Chantis) First implementation of spin-orbit coupling
Cu             (Lz.Sz only)
Cu   19 Sep 03 (ATP) Modifications of CLS spectroscopy
Cu   24 May 03 New --window switch for accumulating density
Cu             in a specific energy window
Cu   24 May 03 New interpolation mode for sigma
Cu   14 Aug 02 Added file input of self-energy addition to LDA
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   24 Oct 01 Updated mixrho, dfrce
Cu   24 Aug 01 Extended to local orbitals.
Cu   22 Apr 01 Added driver for Kotani's GW
Cu   21 Mar 01 bug fix in call to makdos
Cu   20 Mar 01 (ATP) Added Mulliken analysis, CLS
Cu   15 Feb 01 eliminated smrho, smpot from passed arguments
Cu   23 Jan 01 bug fixes connected with lrout=0
C ----------------------------------------------------------------------
      use mpi
      use structures
      implicit none
C ... Passed parameters
C#ifdefC MPE
C      include "mpef.h"
C#endif
      integer numprocs, ierr
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
C#ifdefC MPI
C      integer dims(2)
C#endif
C#ifdefC MPE
CC Event numbers:
C      include "events.ins"
C#endif
      logical mlog
      integer procid,master
      integer nbx,nbas,n0,nab,nppn,iter,maxit,NULLI
      parameter ( nbx=512, n0=10, nppn=12, nab=9, NULLI=-99999)
      double precision NULLR
      parameter (NULLR=-99999d0)
      integer k1,k2,k3,leks,lrout,lfrce,lpnu
      integer nlibu,lmaxu,lldau(nbas)
      double precision dmxp(33),frc(3,nbas,*)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_ctrl)::  s_ctrl
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_bz)::    s_bz
      type(str_optic):: s_optic
      type(str_gw)::    s_gw
      type(str_strn) :: s_strn(*)

C ... Dynamically allocated arrays
      integer, allocatable :: ips(:)   ! species index
      integer, allocatable :: ifbls(:) ! List of orbitals entering into color weights
      integer, allocatable :: idoschan(:) ! Channels for partial DOS, Mulliken
      integer, allocatable :: idoschanl(:)

      real(8), allocatable :: dos(:)   ! total or partial DOS
      real(8), allocatable :: doswt(:,:,:) ! Weights for Mulliken or partial DOS
      real(8), allocatable :: doswtmpi(:,:,:,:) ! To gather partial dos wts for MPI
      real(8), allocatable :: fh(:)    ! Harris-Foulkes forces
      real(8), allocatable :: fes1(:)  ! partial contribution to HF forces (mkpot)
      real(8), allocatable :: fes2(:)  ! partial contribution to HK forces (mkpot)
      real(8), allocatable :: qmom(:)  ! multipole moments of valence sphere densities
      real(8), allocatable :: gpot0(:) ! integrals of gaussians times electrostatic potential
      real(8), allocatable :: vval(:)  ! coffs to YL expansion of es potential at MT boundary
      real(8), pointer     :: hab(:,:) ! augmentation matrices for local parts of hamiltonian
      real(8), pointer     :: vab(:,:) ! augmentation matrices for local parts of potential
      real(8), pointer     :: sab(:,:) ! augmentation matrices for local parts of overlap
      real(8), pointer     :: ppnl(:)  ! NMTO-like potential parameters
      real(8), pointer     :: qbyl(:,:)! l-decomposed site charge
      real(8), allocatable :: hbyl(:,:)! l-decomposed site eigenvalue sum
      real(8), allocatable :: wk(:),ww(:) ! work arrays
      real(8), allocatable :: orbtm(:) ! Orbital moments
      real(8), allocatable :: epsi(:)  ! Imaginary part of epsilon
      real(8), allocatable :: tso(:,:,:,:) ! Resolution of SO contribution to band energy by site
      real(8), allocatable :: rmat(:)  ! Rotation matrix for rotating eigenfunctions

      complex(8), pointer :: z(:,:,:), zr(:,:,:)
      complex(8), allocatable,target :: h(:,:,:),s(:,:,:),hso(:,:,:),dsdk(:,:,:,:)
      complex(8), pointer :: srout(:,:) ! interstitial output density
      complex(8), allocatable :: ausc(:,:,:,:,:)  ! val,slo of w.f. at MT sphere surface
      complex(8), allocatable :: ausu(:,:,:,:,:)  ! val,slo of w.f. at MT sphere surface for LDA+U
      complex(8), allocatable :: ausp(:,:,:,:,:)  ! val,slo of w.f. at MT sphere surface for pdos
      complex(8), allocatable :: cPkLq(:) ! Coefficients to Pkl expansion
      complex(8), allocatable :: zwk(:)   ! Work array
      complex(8), allocatable :: dmatk(:,:) !Density matrix

C     For optics
      integer,allocatable :: iwk1(:),iwk2(:)
      real(8),allocatable :: jdosw(:,:,:,:),wk3(:)
      real(8), allocatable, target:: optmt(:,:,:,:,:)
      complex(8), allocatable, target:: optmc(:,:,:,:,:)
      complex(8), allocatable, target:: optme(:,:,:)
      real(8), allocatable:: evlda(:,:,:),evnl(:,:,:)

C     For energy bands
      real(8),pointer:: evals(:,:,:),colwt(:,:,:,:)
      real(8),allocatable:: evlk(:,:),elohi(:,:,:)
      integer,allocatable:: nevls(:,:)
      type keigs
        integer :: nevk(2),nevlk(2)
        real(8) :: eomink(2)
        real(8),pointer :: evals(:,:)
        complex(8),pointer :: evecs(:,:)
      end type
      type(keigs),allocatable :: karchv(:)

C ... Local parameters
      logical :: lwndow,ltmp
      logical :: spintexture = .false. ! .true. when color weights made by spin texture
      logical :: lwden = .false.   ! .true. when making density for 1 k-point
      integer :: lwoptmc = 0       !  write complex optical matrix elements
      character strn*256,plbopt*256,clsopt*120,strn2*160,sjdosw(2)*160,dc*1
      integer i,ig,i1,i2,ifi,iopq,ipl,ipr,iq,ismidb,modehambls,isp,ispc,isqp,isum,
     .  j,jobgw,jsp,k,lbf,ldim,ldos,lekkl,lfrzw,llmet,lmet,loptic,ltet,lpdiag,
     .  lrep,lrsig,lso,lswtk,lwsig,lwtkb,ldmatk,ltso,mnevl,mpipid,mpsord,mxevcut,
     .  n1,n2,n3,ndham,ndhamx,ndimh,ndimhx,nkaph,nphimx,ndos,nev,nevl,nevmx,nfilem,
     .  nglob,njdosw,nk1,nk2,nk3,nkabc(3),nkp,nl,nsmidb,nsp,nspc,nspc3,nspx,ntet,
     .  nvl,onesp,ovncut,parg,plbnd,stdl,stdo,lmefac,nfstg,morder
      integer lshft(3),ngabc(3),iv(10),nfbn(4),nev_from_evec_file
      integer, parameter :: nkap0=4
      procedure(logical) cmdopt,a2bin,memcheck,latvec
      procedure(integer) :: cmdoptswx,cmdoptsw,iprint,isw,fopn,fopna,fopno,fxst,fadd,
     .  idalloc,allocvb,iobzwt,lgunit,str_pack,iorsa,wordsw
      procedure(real(8)) :: dval
C#ifdefC LMFGWD | LMFGW
C      complex(8),allocatable :: smrox(:)
C#else
      double precision sigp(10)
C#endif
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      equivalence (nk1,nkabc(1)), (nk2,nkabc(2)), (nk3,nkabc(3))
      double precision alat,alfsi,def,del,dosrng,dum,dlength,ebot,
     .  ecbot,ef0,ef00,eferm,efmax,ehar,eks,elind,emax,emin,eomin,
     .  epsovl,esmear,evtop,pi,qbg,qsc,qval,zval,sev,sev00,sev1,sumtv,vol
      double precision plat(3,3),qlat(3,3),sumev(2,3),sumqv(3,2),qp(3),
     .  eterms(22),fsmom(3),xv(20)

C     For dos
      logical doslsw(9),lidos,ldig,llohi,ldsdk,lquitdos
      integer dosisw(5),idfmt
      double precision dosw(2)
      equivalence (lidos,doslsw(1)),(ldig,doslsw(4)),(idfmt,dosisw(2))

      equivalence (emin,dosw(1)),(emax,dosw(2))
      integer, parameter :: LW1=1,LW2=2,LW3=3,LW4=4,LW5=5,LW6=6,LW7=7,LW8=8,LWU=10
      integer, parameter :: lBextz=256,lSOsite=512, lSOEfield=1024
      integer, parameter :: lSOf=4,lSzLz=32,lSzLzp=64,lSzLzp0=128
      logical, parameter :: T=.true., F=.false.
C     For CLS, DOS, Mulliken
      integer, parameter :: nsitmx=256
      real(8),parameter:: tolq=1d-7
      integer nlmax,icls,isite(nsitmx),iclsl(nsitmx),
     .        iclsn(nsitmx),irxes
      integer moddos,nsites,lsites(nbx),nchan,nschan,ng,iomoms,
     .        nchmx,lmdim,lmxch
      integer ifac(3)
      double precision qb(9),qpr(3)

C     For pzhev
      integer nblk,nprow,npcol,nmx
C     For LDA+U, external B field
      integer nbf
      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu),
     .              dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
C     For optics
C     nfilo,nfiup,nemlo,nemup dimension optmt (used at all k)
C     nfiloe,nfiupe,nemloe,nemupe dimension optme (used at one k)
      integer lwtkbsav,lmetsav,lrsigsav,lroutsav
      integer nfilo,nfiup,nemlo,nemup,nfilm,nempm,npopt
      integer nfiloe,nfiupe,nemloe,nemupe ! These dimension optdme
      integer ocrng(2),unrng(2),npol,lteto,rlto
      double precision optrng(4),dwopt(3)
      equivalence (nfilo,ocrng(1)),(nfiup,ocrng(2))
      equivalence (nemlo,unrng(1)),(nemup,unrng(2))
C ... PW basis
      integer pwmode,napw
      double precision pwemin,pwemax,pwgmin,pwgmax
C ... Special AFM symmetry operation
      logical lafms
      double precision afmt(3)
      integer npgrp
C ... For self-energy
      integer nqsig
      double precision eseavr(2)
      real(8), target :: xx(1)

      integer ib,is,lmxa,lmxh,nlma,nlmh,nelt(3),kmax
      integer, allocatable :: kpproc(:)
      double precision sttime,entime
      real(8), allocatable :: eomink(:)  ! Future: merge with karchv%eomink

C      for debugging and testing
C      call thxpbl(s_lat%cg,s_lat%indxcg,s_lat%jcg,s_lat%cy)

      xv = 0
      sumqv = 0
      sumev = 0
      dosisw = 0
      dosw = 0
      doslsw = .false.
      nfstg = 0

      if (maxit == 0) return

C ... MPI setup
      procid = mpipid(1)
      master = 0
      mlog = cmdopt('--mlog',6,0,strn)
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      if (numprocs > 1) then
        call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
        call strcop(shortname(procid),name,10,'.',i)
        namelen(procid) = i-1
      end if

C#ifdefC MPI
C      if (numprocs > nbas) call
C     .  rxi('MPI job cannot allocate more processors than nbas =',nbas)
C#endif

C --- Setup ---
      call tcn ('bndfp')
      napw = 0
      ipr  = iprint()
      ipl  = ipr
      nsp  = nglob('nsp')
      nspc = nglob('nspc')    ! 2 for noncollinear case
      nspc3 = nglob('nspc3')  ! 2 if lso=3 or 4 (spins coupled perturbatively)
      nkaph = nglob('nkaph')  ! number of envelope function types joined to (u,s,phiz)
      nphimx = nglob('nphimx')  ! number of envelope function types joined to (u,s,phiz)
      lso =   isw(IAND(s_ctrl%lncol,4) /= 0)
     .    + 2*isw(IAND(s_ctrl%lncol,lSzLz) /= 0)
     .    + 3*isw(IAND(s_ctrl%lncol,lSzLzp) /= 0)
     .    + 4*isw(IAND(s_ham%lncol,lSzLzp0) /= 0)
      lbf =   isw(IAND(s_ctrl%lncol,8) /= 0)
     .    + 2*isw(IAND(s_ctrl%lncol,lBextz) /= 0)
      ltso = 0; if (lso /= 0) ltso = IAND(s_ctrl%lncol,lSOsite)/lSOsite
      if (s_bz%lmet /= 3 .and. s_bz%lmet /= 5 .and. ltso /= 0) then
        ltso = 0  ! Since FS contributions not properly calculated
        call info0(2,1,0,' bndfp (warning): '//
     .    'SO pert turned off ... requires METAL=3 or 5')
      endif
      s_bz%nef = 1; if (ltso /= 0 .and. nspc == 2) s_bz%nef = 3
      s_bz%egap = NULLI
      loptic = s_optic%loptic
      s_pot%lcplxp = 0
      lekkl = s_ctrl%lekkl
      lmefac = 0; rlto = 0; ldsdk = .false.; allocate(dsdk(1,1,1,1))
      allocate(s_pot%vesrmt(nbas))
      llohi = ipr > 50       ! If T, print out lowest, highest level for each band
      if (cmdopt('--minmax',8,0,strn)) llohi = .true.
      call mpibc1(llohi,1,2,0,'','')
      call ptr_pot(s_pot,8+1,'qbyl',n0*nsp,nbas,[0d0])

C     Conditions that require complex local potential
      if (lso /= 0 .or. lbf /= 0) s_pot%lcplxp = 1  ! SO couplng or external B
      if (isum(nbas,lldau,1) /= 0) s_pot%lcplxp = 1   ! LDA+U
      if (nspc==2 .and. mod(s_ham%lrsa,10) /= 0) s_pot%lcplxp = 1 ! Noncollinear
      nspx = nsp / nspc         ! 1 if nsp=1 or if noncollinear
      ldim   = s_ham%ldham(1)
      pwmode = s_ham%pwmode
      pwemin = s_ham%pwemin
      pwemax = s_ham%pwemax
      nbf = s_ham%nbf
      onesp = 0
      call ivset(nfbn,1,4,0)
      stdo = lgunit(1)
      stdl = lgunit(2)
      nullify(ppnl,evals,colwt)
      morder = 0
      if (cmdopt('--sharm',7,0,strn)) then
        morder = 1
        if (cmdopt('--sharm=0',9,0,strn)) morder = 2
      endif

C     if (stdl < 0) ipl = 0
      ldos = s_ctrl%ldos
      lfrzw= isw(IAND(s_ctrl%lbas,16) /= 0)
      lrsig= s_ham%lsig
      if (procid == master) then
      if (lrsig /= 0
     .    .and. fxst('sigm') /= 1 .and. fxst('sigma') /= 1) then
        call info0(2,1,0,' bndfp (warning): no sigm file found ... DFT calculation only')
        lrsig = 0
      endif
      endif
      call mpibc1(lrsig,1,2,mlog,'rdsigm','lrsig')
      s_ham%lsig = lrsig
      jobgw= -999
      epsovl = s_ham%oveps
      ovncut = s_ham%ovncut
      lpdiag = isw(cmdopt('--pdiag',7,0,strn))
      eomin = 0d0
      mxevcut = NULLI
      mnevl = NULLI
      ndham = s_ham%ndham
      ndhamx = ndham*nspc
      allocate(evlk(ndham,nsp))
      if (llohi) then
        allocate(elohi(ndham,nsp,2))
        call dvset(elohi(1,1,1),1,ndham*nsp,-NULLR)
        call dvset(elohi(1,1,2),1,ndham*nsp,NULLR)
      endif
      ldmatk = min(IAND(s_bz%lio,4+8)/4,2) ! One of 0,1,2
C     if (ldmatk /= 0) call rx('bndfp not ready for dmatk')
C     allocate(s_pot%rhat(nbas),s_pot%rnew(nbas))
      call dfratm(s_site,s_spec,8,1,nbas,s_pot%rhat) ! In case not already done
      s_optic%nlg = 1
      allocate(s_optic%rgrad(1,1,1,1,1)) ! Empty array by default
      allocate(s_optic%rgrade(1,1,1,1,1)) ! Empty array by default

      if (numprocs > 1) lpdiag = 0
C#ifdefC LMFGWD | LMFGW
C      if (cmdopt('-jobgw=',7,0,strn)) then
C        i = 7
C      elseif (cmdopt('--jobgw=',8,0,strn)) then
C        i = 8
C      elseif (cmdopt('--job=',6,0,strn)) then
C        i = 6
C      else
C        i = 0
C      endif
C      if (i /= 0) then
C        if (.not. a2bin(strn,jobgw,2,0,' ',i,-1)) call
C     .    rxs2('BNDFP: failed to parse "',strn(1:30),'%a"')
C      endif
C      if (jobgw /= 1 .and. mpipid(0) > 1)
C     .  call rx('lmfgwd must be run with --job=1 when in multiprocess mode')
C#endif
      pi = 4d0*datan(1d0)
      if (iprint() >= 20) call awrit2('%N --- BNDFP:  '//
     .  'begin iteration %i of %i ---',' ',80,stdo,iter,maxit)
      if (nbas>nbx) call rxi('bndfp: nbx exceeded, need',nbas)
      if (ndham<=0)
     .  call rx('bndfp: hamiltonian matrix has zero dimension')

C ... Average density by groups, render rho positive
      if (iand(s_ctrl%lcd,32) /= 0) then
        call grp2av(s_ctrl,s_site,s_spec,s_pot%rhat)
      endif
      if (iand(s_ctrl%lcd,128) /= 0) then
        i = 1
        call rhopos(i,s_site,s_spec,s_lat,s_pot%rhat,s_pot%smrho)
      endif

C#ifdefC MPI
CC MPI Process configuration
C      if (lpdiag == 1) then
C      nblk = 16
C      dims(1) = 0
C      dims(2) = 0
C      call MPI_DIMS_CREATE(numprocs,2,dims,ierr)
C      npcol = dims(1)
C      nprow = dims(2)
C      if (iprint() >= 30) then
C        call awrit3(
C     .     ' MPI creating process configuration .. nprow=%i npcol=%i,'//
C     .     ' blocking factor %i',
C     .     ' ',256,lgunit(1),nprow,npcol,nblk)
C        endif
C      endif
C#endif

      call dvset(eterms,1,22,NULLR)
      eterms(19) = 0d0
      s_ham%eterms = eterms
      nvl = s_pot%nlml
      nchan = s_pot%nlma
      eks = 0
      alat = s_lat%alat; plat = s_lat%plat; qlat = s_lat%qlat; vol = s_lat%vol
      ngabc = s_lat%nabc
      lafms = s_lat%nsafm < 0
      if (nspc == 2 .and. lafms)
     .  call rx('bndfp: AFM symop not allowed in conjunction with noncol')

C ... Setup and optimize parameters for FFT
      call fftz30(n1,n2,n3,k1,k2,k3)
      allocate(zwk(n1*n2*n3))
      call fftzv(zwk,n1,n2,n3,1,1,21,-1)
      deallocate(zwk)

C ... Switch to generate matrix elements for velocity operator
      lwoptmc = max(cmdoptsw('--opt',',','woptmc',''),0)
C      if (lwoptmc > 0 .and. cmdoptsw('--opt',',','rdqp','') > 0) then
C        call info0(20,0,0,' ... Read qp for optical matrix elements from qpts file')
C        nkp = 0
C        ifi = fopno('qpts')
C        call getqp(0,ifi,nkp,s_bz%nkabc,s_bz%lshft,s_bz%ntet,xx,xx,xx)
C        if (nkp <= 0) call rx('improper or missing contents of qpts file')
C        if (associated(s_bz%star)) deallocate(s_bz%star)
C        s_bz%nkp = nkp
C        call ptr_bz(s_bz,8+1,'qp',3,nkp,xx)
C        call ptr_bz(s_bz,8+1,'wtkp',nkp,0,xx)
C        call ptr_bz(s_bz,1,'idtet',5,max(s_bz%ntet,1),xx)
C        call getqp(2,ifi,nkp,s_bz%nkabc,s_bz%lshft,s_bz%ntet,s_bz%qp,s_bz%wtkp,s_bz%idtet)
C        call fclose(ifi)
C      endif

C ... for BZ integration
      nkabc = s_bz%nkabc
      nkp = s_bz%nkp
      ntet = s_bz%ntet
      lmet = s_bz%lmet
      llmet = 0
      mpsord = s_bz%n
      esmear = s_bz%w
      nevmx = s_bz%nevmx
      efmax = s_bz%efmax
      fsmom(1:2) = s_bz%fsmom
      esmear = esmear + iabs(mpsord)
      if (mpsord < 0) esmear = -esmear
      mpsord = int(dabs(esmear))
      if (esmear < 0) mpsord = -mpsord

      ndos = s_bz%ndos
      dosw = s_bz%dosw
      ef0 = s_bz%ef
      def = s_bz%def
      ltet = mod(s_ctrl%lmet/2,2) + 10*mod(s_ctrl%lmet/16,2)
C#ifndef LMFGWD | LMFGW
      call rxx(ltet /= 0.and.ntet==0,
     .  'bndfp:  tetrahedron integration set but tetrahedra not available')
C#endif
      lwndow = cmdopt('--window=',9,0,strn)
      if (lwndow) then
        if (lmet == 0) then
          call rx(' bndfp: restart with METAL=2 for --window')
        endif
        if (ltet == 0) then
          call rx(' bndfp: restart with TETRA=T for --window')
        endif
        iq = 0
        i = parg('window=',4,strn,iq,len(strn),', ',2,2,iv,dosw)
        call info2(20,0,0,
     .    ' BNDFP: generating density in energy window %2:1d',dosw,0)
        lfrce = 0
        lpnu = 0
        efmax = 1d3
        nevmx = ndham
        if (lrout == 0) call rx('--window incompatible with no output density')
        call info0(20,0,0,' Delete band weights file ...')
        ifi = fopna('wkp',-1,4)
        call dfclos(ifi)
      endif
C     Always need DOS if just sampling
      if (lmet /= 0 .and. ltet /= 1) ldos=1
      qbg = s_ctrl%zbak(1)
      alfsi = s_ham%alfsi
C     Switch to plot bands at specified qp
      strn = ' '
C     Band pass for plotting
      if (cmdopt('--band',6,0,strn) .or. s_ctrl%plbnd == 1) then
        plbnd = 1; plbopt = strn(7:)
        lrout = 0; lfrce = 0
        nkp = 0; s_bz%numq = 1
        allocate(ifbls(ndhamx*3)); call iinit(ifbls,size(ifbls))
        if (lmet == 5) lmet = 3
C     Band pass just to make bands, integrated over BZ
      elseif (s_ctrl%plbnd == 2 .or. s_ctrl%plbnd == 3) then
        plbnd = s_ctrl%plbnd; plbopt = ' '
        lrout = 0; lfrce = 0
        s_bz%numq = 1; if (lmet == 5) lmet = 3
C     Band pass
      else
        plbnd = 0; s_bz%numq = 1
      endif
      if (cmdopt('--optbas',8,0,strn)) then
        if (lmet == 5) lmet = 3
      endif

C     Sanity checks
      if (lmet == 1) call rx('bndfp: lmet=1 not implemented')
      call sanrg(.true.,lmet,0,5,'bndfp:','lmet')
      lmetsav = NULLI

C --- Define local arrays used in the generation of the potential ---
      allocate(qmom(nvl),gpot0(nvl),vval(nchan),fes1(3*nbas))
      i = n0*nsp*nbas
C      allocate(hab(nab*i))
C      allocate(sab(nab*i))
      call ptr_pot(s_pot,1,'hab',nab*n0*nsp,nbas,xv)
      call ptr_pot(s_pot,1,'sab',nab*n0*nsp,nbas,xv)
      hab => s_pot%hab; sab => s_pot%sab
      allocate(vab(nab*n0*nsp,nbas))
      if (associated(s_pot%ppn)) then  ! If T, allocate ppn as s_pot%ppn to preserve
        deallocate(s_pot%ppn)
        allocate(s_pot%ppn(nppn*i))
        ppnl => s_pot%ppn
      else
        allocate(ppnl(nppn*i))
      endif
C     i = idalloc('site',allocvb()+2,3*nab+nppn,i)
      call dfaugm(nbas,s_pot%lcplxp,lso,0,s_site,s_spec)

C#ifndef LMFGWD | LMFGW
      if (mod(loptic,10)>0) then  ! Calculate epsilon_2
        call lkdim(1,nbas,s_site,s_spec,s_optic%kmxax,nlmax)
C       s_optic%nlg = s_ctrl%nl
        s_optic%nlg = nlmax+1
C       s_optic%nlg = s_optic%nlg + 1     !debugging check
C       s_optic%kmxax = s_optic%kmxax + 1 !debugging check
        deallocate(s_optic%rgrad,s_optic%rgrade)
        allocate(s_optic%rgrad(nkaph**2,2,s_optic%nlg,nsp*nspc,nbas))
        i = (s_optic%kmxax+1+nkaph)**2
        allocate(s_optic%rgrade(i,2,s_optic%nlg,nsp*nspc,nbas))
      endif
C#endif

C --- Make the potential sans XC part ---
C#ifdefC LMFGWD | LMFGW
C      if (jobgw == 1 .or. jobgw == 2 .or. jobgw == -999) then
C        call info(20,1,0,' Make potential without XC part ...',0,0)
C        call togpr()
C        i = 1 + 10*lfrzw + 100
C        call dfaugm(nbas,s_pot%lcplxp,lso,1,s_site,s_spec)
C        allocate(s_pot%smpoth(k1*k2*k3,nsp)); call dpzero(s_pot%smpoth,2*k1*k2*k3*nsp)
CC        allocate(smrox(k1*k2*k3*nsp)); call dpzero(smrox,2*k1*k2*k3*nsp)
CC        call dcopy(2*k1*k2*k3*nsp,s_pot%smrho,1,smrox,1)
C        call mkpot(s_site,s_spec,s_lat,s_ham,s_pot,s_optic,nbas,0,k1,k2,k3,s_pot%smrho,
C     .    s_pot%rhat,qbg,s_pot%smpoth,qmom,s_pot%smvcnst,ppnl,hab,vab,sab,
C     .    qval,qsc,gpot0,vval,fes1,i,vorb,nlibu,lmaxu,lldau,nbf,s_ham%magf)
C        call togpr()
CC       deallocate(smrox)
C      endif
C#endif

C --- Make the potential and total energy terms ---
      i = 1 + 10*lfrzw
      if (cmdopt('--wrhomt',8,0,strn)) then
        i = i + 10000
      else if (cmdopt('--wpotmt',8,0,strn)) then
        i = i + 20000
      endif
C#ifdefC LMFGWD | LMFGW
C      if (s_gw%code == 0 .or. s_gw%code == 2) i = i + 10000
C#endif
      if (iand(s_ctrl%lcd,128) /= 0) i = i + 300000
      if (iand(s_ctrl%lbxc,16) /= 0) i = i + 200     ! Use spin-avgd rho => Bxc=0
      call togpr()
      if (cmdopt('--SOefield',10,0,strn)) then  ! Switch to print out SO-weighted average E-field
        print *, iand(s_ham%lncol,lSOEfield)
        s_ham%lncol = s_ham%lncol + lSOEfield - iand(s_ham%lncol,lSOEfield)
      endif
      if (.not. associated(s_ham%magf)) s_ham%magf => xx
      call mkpot(s_site,s_spec,s_lat,s_ham,s_pot,s_optic,nbas,lfrce,k1,k2,k3,s_pot%smrho,
     .  s_pot%rhat,qbg,s_pot%smpot,qmom,s_pot%smvcnst,ppnl,hab,vab,sab,
     .  qval,qsc,gpot0,vval,fes1,i,vorb,nlibu,lmaxu,lldau,nbf,s_ham%magf)
      s_pot%qval = qval
      zval = qval
      if (s_bz%zval /= NULLI) zval = s_bz%zval

      if (.not. cmdopt('--asars2',8,0,strn) .and. cmdopt('--asars',7,0,strn)) then
        i = fopna('rsta',-1,0)
        strn = 'ASA format, gen by lmf'
        i = iorsa(1,s_ctrl,s_site,s_spec,s_lat,s_pot,s_bz,strn,nbas,nbas,nglob('nspec'),
     .    xx,iter,.false.,-i)
        call rx0('ASA style rsta file written')
      endif

      if (cmdopt('--SOefield',10,0,strn)) then
        s_ham%lncol = s_ham%lncol - iand(s_ham%lncol,lSOEfield)
      endif
      call togpr()
      if (cmdopt('--quit=pot',10,0,strn)) then
        call rx0('quit = pot')
      endif
      if (s_ctrl%plbnd == 10) return

C --- Setup for optics ---
C#ifndef LMFGWD | LMFGW
      allocate(optmt(1,1,1,1,1))
      if (loptic /= 0) then
        lteto = s_optic%ltet
        ocrng = s_optic%ocrng(1:2)
        unrng = s_optic%unrng(1:2)
        dwopt = s_optic%dw
        optrng = s_optic%window
        i = 1
        if (dwopt(2) /= 0) i = 11
        call dosmsh(-i,optrng,dwopt,0,npopt,xv)
        if (lmet > 0) then
          if (nfilo <= 0) nfilo = 1
          if (nfiup <= 0) nfiup = ndhamx
          if (nemlo <= 0) nemlo = 1
          if (nemup <= 0) nemup = ndhamx
        else
          if (nfilo <= 0) nfilo = 1
          if (nfiup <= 0) nfiup = (nspc*nint(zval)+1)/2
          if (nemlo <= 0) nemlo = nfiup+1
          if (nemup <= 0) nemup = ndhamx
        endif
        nemup = min(nemup,ndhamx)
        nfiup = min(nfiup,ndhamx)
        if (loptic == -5 .or. loptic == -6) then
          nemlo = nfilo
          nemup = nfiup
        else if (nfiup > nemup) then
          call info2(30,0,0,' (warning, optics) last filled '//
     .      'state (#%i) above last empty state (#%i)',nfiup,nemup)
C         nfiup = nemup
        endif
        if (nemlo < nfilo) then
          call info2(30,0,0,' (warning, optics) first empty '//
     .      'state (#%i) below first filled state (#%i)',nemlo,nfilo)
C         nemlo = nfilo
        endif
        s_optic%ocrng(1:2) = ocrng
        s_optic%unrng(1:2) = unrng
        nfilm = nfiup-nfilo+1
        nempm = nemup-nemlo+1
        s_optic%nfilm = nfilm
        s_optic%nempm = nempm
        call ivset(nfbn,1,4,0)

C       Allocate memory for optical matrix elements
        if (mod(loptic,10)>0) then
          i = mpipid(0)
          i = idalloc('optmt',allocvb()+2,2*3*nfilm*nempm,nkp*i)
          ltmp = memcheck('bndfp','optical matrix elements',s_ctrl%maxmem,T)
          deallocate(optmt)
          allocate(optmt(3,nfilo:nfiup,nemlo:nemup,nspx,nkp))
          call dpzero(optmt,size(optmt))
          if (lwoptmc > 0) then
            allocate(optmc(3,nfilo:nfiup,nemlo:nemup,nspx,nkp))
            call dpzero(optmc,2*size(optmc))
          else
            allocate(optmc(1,nfilo:nfilo,nemlo:nemlo,nspx,nkp)) ! Dummy array
          endif
        endif
      endif

C ... Parse switch '--jdosw': get list of projections
      sjdosw = ' '
      njdosw = 0
      if (loptic < 0 .and. cmdopt('--jdosw',7,0,strn)) then
        if (strn(1:8) == '--jdosw2')
     .    call rx('BNDFP: --jdosw must precede --jdosw2')
        sjdosw(1) = strn(8:)
        call getjdosw(1,8,strn,nfbn(1),nfbn(3),i,i)
        allocate(iwk1(maxval(nfbn)))
        allocate(iwk2(maxval(nfbn)))
        call getjdosw(3,8,strn,nfbn(1),nfbn(3),iwk1,iwk2)
        call ilst2a(iwk1,nfbn(1),strn2)
        call ilst2a(iwk2,nfbn(3),strn)
        call info2(20,1,0,
     .    ' BNDFP:  Mulliken projection of %?#(n==-5|n==-6)#DOS#JDOS#',loptic,0)
        call info2(20,0,0,'%9fgroup:  '//
     .    '%?#n#%-1j%i channel(s): '//strn2//'%a (occ)  ##' //
     .    '%?#n#%-1j%i channel(s): '//strn//'%a (unocc)  ##',
     .    nfbn(1),nfbn(3))
        if (loptic == -5 .or. loptic == -6) then
          if (nfbn(3) /= 0)
     .      call rxi(' unocc weights not allowed for loptic =',loptic)
        endif
        deallocate(iwk1,iwk2)
        allocate(jdosw(1,nfilm+nempm,nkp,nspx))
        njdosw = 1
      endif
      if (loptic < 0 .and. cmdopt('--jdosw2',8,0,strn)) then
        sjdosw(2) = strn(9:)
        call getjdosw(1,9,strn,nfbn(2),nfbn(4),i,i)
        allocate(iwk1(maxval(nfbn)))
        allocate(iwk2(maxval(nfbn)))
        call getjdosw(3,9,strn,nfbn(2),nfbn(4),iwk1,iwk2)
        call ilst2a(iwk1,nfbn(2),strn2)
        call ilst2a(iwk2,nfbn(4),strn)
        call info2(20,0,0,'%9fgrp 2:  '//
     .    '%?#n#%-1j%i channel(s): '//strn2//'%a (occ)  ##' //
     .    '%?#n#%-1j%i channel(s): '//strn//'%a (unocc)  ##',
     .    nfbn(2),nfbn(4))
        if (loptic == -5 .or. loptic == -6) then
          if (nfbn(4) /= 0)
     .      call rxi(' unocc weights not allowed for loptic =',loptic)
        endif
        deallocate(iwk1,iwk2)
        deallocate(jdosw)
        allocate(jdosw(2,nfilm+nempm,nkp,nspx))
        call dpzero(jdosw,size(jdosw))
        njdosw = 2
      endif
      if (.not. allocated(jdosw)) allocate(jdosw(1,1,1,1))
      call dpzero(jdosw,size(jdosw))
C#endif

C ... Other hamiltonian setup
      if (ltso == 1) then
        allocate(tso(2,2,2,0:nbas)); call dpzero(tso,8*(nbas+1))
      else
        allocate(tso(1,1,1,0:0))
      endif
      call dpzero(s_bz%sopertp,12)
      call suham2(s_ctrl,s_lat,s_spec,s_site,s_ham,s_pot,s_strn)

C ... Read and store self-energy addition to LDA potential
      lwsig = 0
      if (mod(lrsig,10) /= 0) then
C       Real-space range
        ifi = 1
C       if (procid == master) ifi = fopna('sigm',-1,4)
        if (procid == master) call sigfileunit(strn,ifi,ltmp)
        call mpibc1(ifi,1,2,.false.,'','')
        call rdsigm(lrsig,nbas,nsp,s_ham%nlmto,s_ctrl,s_spec,s_lat,s_ham,s_bz,s_gw,ifi,lwsig)
        call rxx(lwsig /= 0 .and. plbnd /= 0,'incompatible options, lwsig and plbnd')
        if (lwsig > 0 .and. lwsig <= LW4) then
          lwtkb = 0; lrout = 0; lfrce = 0
        endif
        call fclose(ifi)
        call phmbl3(1,0,0,0,0,0,0,0)
      endif
C     Check for other special modes that read/write h, sig, or evecs to disk
      if (s_ctrl%lwsig /= 0) then
        if (lwsig /= 0 .and. lwsig /= s_ctrl%lwsig)
     .    call rx2('incompatible options lwsig=%i and lwsig=%i',s_ctrl%lwsig,lwsig)
        lwsig = s_ctrl%lwsig
      endif
C     if (lwsig == LW3 .and. lmet == 5) lmet = 3
      if (lwsig == LW3 .or. lwsig == LW5 .or. lwsig == LW6 .or. lwsig == LW7) then
        if (lmet == 5) call rx('For trans restart with METAL=2 or 3')
        nevmx = ndhamx
        if (lwsig == LW5 .or. lwsig == LW6) then
          lwtkb = 0; lfrce = 0
        endif
        if (lwsig == LW6) then
C         nkp = s_bz%nkp  ? maybe after call to rdsigm?
          ldmatk = 1
          if (lrout == 0) call rxx(lrout==0,'incompatible options lwsig+lrout=0')
        elseif (lwsig == LW5) then
          lrout = 0
        endif
      endif
      s_ctrl%lwsig = lwsig

C     Transformation modes: read qp and skip over BZ setup
C#ifndef LMFGWD | LMFGW
      if (lwsig > 0 .and. lwsig <= LW4) then
        call sanrg(.true.,lwsig,0,4,'BNDFP:','lwsig')
        onesp = 0; lwtkb = 0; nsmidb = ndham
        ifi = fopna('qpts',-1,0)
        call getqp(0,ifi,nkp,nkabc,lshft,i,xv,xv,xv)
        call ptr_bz(s_bz,1,'qp',3,nkp,xv)
        call getqp(1,ifi,nkp,nkabc,lshft,i,s_bz%qp,xv,xv)
        goto 50
      endif
C#endif

C#ifdefC LMFGWD | LMFGW
C      if (s_gw%mksig > 0 .and. epsovl /= 0) then
C        call logwarn(2,'%N bndfp (warning!) using OVEPS with QSGW is not advised at present')
C      endif
C#endif

      elind = s_ham%elind
      if (elind < 0d0) elind=-(3*pi**2*(zval-qsc-qbg)/vol)**.66666d0*elind
      if (lrout /= 0) dmxp(33) = elind
      if (iand(s_ctrl%lfp,8) /= 0) s_ham%elind = elind

C --- Setup for BZ integration ---
      if (plbnd == 0 .or. plbnd == 2) then
        if (plbnd == 0) then
          i = idalloc('evals',allocvb()+2,ndham,nsp*nkp)
          allocate(evals(ndham,nsp,nkp))
        else
          evals => s_ham%evals
          iv(1:3) = shape(evals)
          if (iv(1) /= ndham .or. iv(2) /= nsp .or. iv(3)<nkp)
     .      call rx('bndfp: s_ham=>evals improperly dimensioned')
        endif
        if (cmdopt('--dos',5,0,strn)) nevmx = ndhamx
        call subzi(s_bz,lmet,ltet,lrout>0,s_bz%nef,ndham,nsp,nspc,
     .    nkp,zval-qbg,nevmx,lwtkb,eferm,lswtk,ef0)
C       ? not sure ... don't need weights if no
        if (lrout==0 .and. lwtkb==-1) lwtkb=0
C#ifndef LMFGWD | LMFGW
        if (lwtkb == -1) call info(20,1,0,' Start first of two band passes ...',0,0)
C#endif
        if (lwtkb == 1) then
          if (ef0 /= eferm) call info(20,0,1,'%8pReplace ef0 with file ef=%;6d',eferm,0)
          ef0 = eferm
        endif
C Future : merge with karchv%eomink ... information is the same
        if (numprocs > 1) then
        if (epsovl /= 0 .or. ovncut /= 0) then
          allocate(eomink(nsp*nkp))
        endif
        end if
      else
        if (plbnd == 3) then
          evals => s_ham%evals
          if (size(evals) /= ndham*nsp*nkp)
     .      call rx('bndfp: s_ham=>evals improperly dimensioned')
        else if (plbnd == 1) then
          nkp = 0
        endif
        ldos = 0; lwtkb = -1; icls = 0; lswtk = 0
      endif
      lquitdos = cmdopt('--quit=dos',10,0,strn) .and.
     .   .not. (cmdopt('--mull',6,0,strn).or.cmdopt('--pdos',6,0,strn))
      if (lquitdos) lwtkb = -1

C ... Special mode to generate density at one k-point
      qpr(1) = NULLI
      if (cmdopt('--wden',6,0,strn)) then
        if (plbnd /= 0) call rx('modes band and wden are incompatible')
        k = -ndham; allocate(iwk1(iabs(k)))
        call ioden2(1,s_lat,strn(7:),xv,xv,xv,xv,xv,xv,xv,xv,xv,xv,xv,qpr,k,iwk1,xv)
        lwden = qpr(1) /= NULLI
        if (lwden) then
          call info5(20,0,1,' ioden : make density for qp =%3:1,6;6d'//
     .      ' bands%n:1i',qpr,k,iwk1,0,0)
          if (k < 0) call rx('no bands specified')
          nkp = 1; lwtkb = 1; lfrce = 0; lmet = 0
          lmet = 0; loptic = 0; s_bz%numq = 1; s_bz%nef = 1; s_bz%wtkp(1) = 2d0
C         Assign band weights
          call ptr_bz(s_bz,8+1,'wtkb',ndham*nsp,0,xv)
          do  i = 1, k
            s_bz%wtkb(iwk1(i)) = 1; if (nspx == 2) s_bz%wtkb(iwk1(i)+ndham) = 1
          enddo
          call imxmn(k,iwk1,1,i,nevmx)
          deallocate(iwk1)
        endif
      endif

C ... Set switch to write file sigii
C#ifndef LMFGWD | LMFGW
      if (lrsig /= 0 .and. plbnd == 0 .and. procid == master) then
        i = fopna('sigii',-1,0)
        rewind i
        s_ham%sigp(9) = 1d0
      endif
C#endif

C ... More k-point independent local arrays
      if (ldos == 0 .and. lmet == 0) ndos = 1
      if (ndos < 0) ndos = -ndos
      allocate(dos(ndos*2*nsp))
      dos = 0
      if (lrout /= 0) then
        call dfqkkl(nbas,lekkl,s_site,s_spec,s_bz%numq)
        call ptr_pot(s_pot,8+1,'smrout',k1*k2*k3,s_bz%numq*nsp*nspc,xv)
        srout => s_pot%smrout
C       allocate(srout(k1*k2*k3*s_bz%numq*nsp*nspc))
        allocate(fh(3*nbas))
        allocate(fes2(3*nbas))
      else
        allocate(srout(1,1))  ! DEC compiler complains if not allocated
      endif

C ... Options for core level specta (CLS)
      if (cmdopt('--cls',5,0,strn)) then
        if (lmet > 0 .and. (lmet /= 2 .and. lmet /= 3))
     .    call rx('For CLS restart with METAL=2 or 3')
        icls = 1
        efmax = 1d3
        nevmx = ndhamx
        if (lrout == 0) call rx('bndfp: need output density for cls')
        call suclst(nsitmx,nbas,nsp,s_site,s_spec,strn(6:),isite,iclsl,iclsn,nsites)
        efmax = 1d3
        nevmx = ndhamx
        nlmax = nglob('mxorb') / nkaph
        i = idalloc('ausc',allocvb()+2,nlmax*ndhamx*nphimx,nsites*nsp*nkp*2)
        allocate(ausc(nlmax,ndhamx,nphimx,nsp,nsites*nkp)); call dpzero(ausc,2*size(ausc))
        if (s_optic%loptic /= 0) call rx('optics and CLS cannot be calculated together')
      else
        icls = 0
      endif
C ... Options for RXES
      if (cmdopt('--rxes',6,0,strn)) then
        if (lmet > 0 .and. (lmet /= 2 .and. lmet /= 3))
     .    call rx('For RXES restart with METAL=2 or 3')
        irxes = 1
        clsopt = strn(7:)
        efmax = 1d3
        nevmx = ndhamx
        if (lrout == 0) call rx('bndfp: need output density for cls')
        call suclst(nsitmx,nbas,nsp,s_site,s_spec,clsopt,isite,iclsl,iclsn,nsites)
        efmax = 1d3
        nevmx = ndhamx
        nlmax = nglob('mxorb') / nkaph
        i = idalloc('ausc',allocvb()+2,nlmax*ndhamx*nphimx,nsites*nsp*nkp*2)
        allocate(ausc(nlmax,ndhamx,nphimx,nsp,nsites*nkp)); call dpzero(ausc,2*size(ausc))
      else
        irxes = 0
      endif

C#ifdefC LMFGW
CC --- GW code ---
C      strn = ' '
C      call makegw(s_ctrl,s_site,s_spec,s_lat,s_ham,s_pot,s_bz,s_gw,nbas,ppnl)
C      call tcx('bndfp')
C      call rx0('bndfp')
C#endif

C#ifdefC LMFGWD
CC --- GW driver ---
C      strn = ' '
C      call sugw(s_ctrl,s_site,s_spec,s_lat,s_ham,s_pot,s_bz,s_gw,strn,nbas,ppnl,jobgw)
C      call tcx('bndfp')
C      call rx0('bndfp')
C#endif

C#ifndef LMFGWD | LMFGW

C ... Skip to optics, reading data from file
      if (cmdopt('--opt:read',10,0,strn)) then
        goto 499
      endif

C --- Start loop over k points; also, re-entry for second band pass ---
      nsmidb = ndham
   99 continue
      ebot = 1000d0
      call surho(nbas,s_site,s_spec,lmet,ldos,lrout,lekkl,s_bz%numq,
     .  k1,k2,k3,srout,ndos,dos,sumev,sumqv)
      if (lfrce > 0) then
        call dpzero(frc, 3*nbas*s_bz%numq); call dpzero(fh,3*nbas)
      endif

C ... Optics: initial printout and check on nevmx
      if (loptic /= 0) then  ! Extra unocc bands might be needed
        i = nint((zval-qbg)*nspc/2)
        if (loptic == -5 .or. loptic == -6) then
          call info5(30,1,0,' optics: DOS for states (%i,%i):'//
     .      '  %i states below Ef',ocrng,ocrng(2),i,4,5)
        else
          call info5(30,1,0,' optics: transitions for occ=(%i,%i)  '//
     .      'unocc=(%i,%i):  %i states below Ef',ocrng,ocrng(2),unrng,unrng(2),i)
        endif
        if (nfilo >= i+1) then
          call info0(30,0,0,'%9f(warning) first filled state above Ef')
        endif
        if (loptic /= -5 .and. loptic /= -6 .and. nemup <= i) then
          call info0(30,0,0,'%9f(warning) last empty state below EF')
        endif
        if (efmax < s_optic%window(2)+1d0) then
          call info2(20,0,1,'%9f(warning) increasing EFmax from %,1d Ry to ommax+1 (%,1d Ry)',
     .      efmax,s_optic%window(2)+1d0)
          efmax = s_optic%window(2)+1d0
        endif
        i = nevmx
        nevmx = max(nevmx,unrng(2))
        if (nevmx /= i) then
         call info2(30,0,0,'%9fincrease no. calc. evals to %i ',nevmx,0)
        else
          if (iprint() >= 20) write(stdo,"(1x)")
        endif
        if (mod(loptic,10)>0 .and. lrsig /= 0 .and. plbnd==0 .and.
     .    (lwsig==0.or.lwsig==8) .and. s_optic%mefac==1) then
          ldsdk = .true.
          call info0(30,0,0,'%9fCalculate dSigma/dk for velocity operator')
          if (lmet == 5) call rx('dSigma/dk does not yet work with metal=5, sorry')
        endif
        call info0(30,0,0,' ')
      endif
C ... Setup for approximate correction for NL sigma in optics
   40 continue
      if (mod(loptic,10)>0 .and. lrsig /= 0 .and. plbnd==0 .and.
     .        lmefac==0 .and. (lwsig==0.or.lwsig==8) .and. s_optic%mefac==2) then
        lmefac = 1              ! Flag for special LDA band pass before optics
        lwtkbsav=lwtkb; lmetsav=lmet; lrsigsav=lrsig; lroutsav=lrout
        lwtkb=-1; lmet=0; lrsig=0; lrout=0; loptic=0
C        if ( .not. allocated(evlda)) allocate(evlda(ndham,nsp,nkp))
C        if ( .not. allocated(evnl)) allocate(evnl(ndham,nsp,nkp))
        allocate(evlda(nfilo:nemup,nspx,nkp))
        allocate(evnl(nfilo:nemup,nspx,nkp))
        call info0(30,0,1,' bndfp:  make initial LDA band pass for optics correction factor')
      endif

      if (lswtk == 1) then  ! Need the entire evec matrix to invert it
        efmax = 1d3
        nevmx = ndhamx
        call dpzero(s_bz%swtk,ndhamx*nkp)
        if (lwtkb == -1 .and. lmet==5) lswtk=-2 ! need evecs for spin weights
      endif

C --- Setup moments file : write header ---
      nl = s_ctrl%nl
      nfilem = 0
C     MPI case: separate moms files for each procid
      if (procid == master .or. cmdopt('--mull',6,0,strn).or.cmdopt('--pdos',6,0,strn)) then
        nfilem = fopna('moms',-1,4)
      endif
      nlmax = nglob('nlmax') ! All nodes require this for partial DOS
      if (nfilem > 0) then
        nfstg = 1 ! Default nfstg : bands only, no DOS weights
        rewind nfilem
        i = 1
C   ... If switch '--mull'; get sites, number of channels
        if (cmdopt('--mull',6,0,strn).or.cmdopt('--pdos',6,0,strn)) then
          ng = s_lat%nsgrp
          efmax = 1d3
          nevmx = ndhamx
          nchmx = min(1024,nbas*nl**2)
          allocate(idoschanl(nchmx))
          j = 0
          if (cmdopt('--mull',6,0,strn)) j = 1
          i = 0
          call sumlst(j,nchmx,nbas,ng,s_site,s_spec,strn(7:),moddos,
     .      nsites,lsites,lmxch,nchan,lmdim,idoschanl,i)
          if (lwtkb >= 0) then
            if (allocated(idoschan)) deallocate(idoschan)
            allocate(idoschan(lmdim*nsites))
            call icopy(lmdim*nsites,idoschanl,1,idoschan,1)
          endif
          deallocate(idoschanl)
          if (cmdopt('--pdos',6,0,strn)) then
            nlmax = nglob('nlmax')
C            i = nlmax*ndham*3*nsp*nbas
C            if (16*i*nkp/1000000 > 24 .and. procid == master)
C     .        call info(20,0,0,' PDOS: %iMb memory for aus: nlmax=%i',i/1000000,nlmax)
            if (cmdopt('--mull',6,0,strn)) call rx('--pdos and --mull not allowed in conjunction')
            if (cmdopt('--cls',5,0,strn)) call rx('--pdos and --cls not allowed in conjunction')
          endif
          call iomomn(.true.,2,.false.,1,nspc,1,1,nfstg)
        endif
        i = iomoms(-nfilem,nl,nsp,nspc,nkp,
     .    ndham,nfstg,1,0,1,0,0,0,0,0,0d0,0d0,0d0,0d0,0d0,0d0)
      endif

C ... Case only generate bands at supplied qp: setup
      if (plbnd == 1) then
C       Try and read Fermi level from file
        if (procid == master) then
          ifi = fopna('wkp',-1,4); call getef(ifi,0,ef0); call fclr('wkp',ifi)
        endif

        call mpibc1(ef0,1,4,.false.,'','')
        iopq = 0
C       suqlst in MPIK mode; returns cumulative number of k-points
        if (numprocs > 1) iopq = 2
        if (cmdopt('--onesp',7,0,strn) .and. nspc == 1) onesp = 1
        i = nsp
        if (onesp /= 0 .or. nspc == 2 .or. lafms) i = 1
C       In parallel mode, suqlst call only serves to generate nkp
        if (procid == master) then
          if (pwmode > 10) iopq = 10+iopq
          call suqlstst(s_lat,plbopt,iopq,ndhamx,ef0,i,xv,nfbn,ifbls,nkp,qp,onesp,spintexture)
        endif
        call mpibc1(nkp,1,2,.false.,'','')
        if (nkp <= 0) call rx0('bndfp')
        call mpibc1(spintexture,1,1,.false.,'','')
        call mpibc1(nfbn,size(nfbn),2,.false.,'','')
        call mpibc1(ifbls,size(ifbls),2,.false.,'','')
        call mpibc1(onesp,1,2,.false.,'','')

C   ... MPIK: Setup to assemble all k-points into single list with qp table
        if (numprocs > 1) then

          call suqlsm(i)
          call info2(20,1,1,' bndfp:  MPIK band plotting mode %i:  %i q-points',i,nkp)
C         Re-allocate qp and evl arrays
          call ptr_bz(s_bz,1,'qp',3,nkp,xv)
          if (associated(evals)) deallocate(evals)
          i = idalloc('colwt',allocvb()+2,ndhamx*nspx,4*nkp)
          allocate(evals(ndham,nsp,nkp),colwt(ndhamx,nspx,4,nkp))
          if (epsovl /= 0 .or. ovncut /= 0) then
            if (allocated(eomink)) deallocate(eomink)
            allocate(eomink(nsp*nkp))
          endif
C         Loop through all qp; accumulate vector of qp.
C         Use i2 in place of nkp to preserve nkp
          if (procid == master) then
C         iq = running index to big qp list, i1 = index to current line
            iq = 0
            i2 = 0
  199       continue
            i = 1
            call pshpr(1)
            iopq = 1; if (pwmode > 10) iopq = 10+iopq
            call suqlstst(s_lat,plbopt,iopq,ndhamx,ef0,i,xv,nfbn,ifbls,i2,qp,onesp,spintexture)
            call poppr
            if (i2 > 0) then
              call pshpr(1)
              do  i1 = 1, i2
                iq = iq+1
                iopq = 1; if (pwmode > 10) iopq = 10+iopq
                call suqlst(s_lat,plbopt,iopq,ndhamx,ef0,i,xv,nfbn,ifbls,i2,qp,onesp)
                call dpscop(qp,s_bz%qp,3,1,3*iq-2,1d0)
              enddo
              call poppr
              call suqlsm(i)
              if (i /= 2 .and. i /= 3) goto 199
            endif
          endif
          call mpibc1(s_bz%qp,3*nkp,4,.false.,'bndfp','qp')
        endif
      endif

C ... For insulator, valence band top and conduction band bottom
      evtop = -99999
      ecbot = -evtop

C ... Setup for spin-orbit coupling
      if (lso /= 0) then
C        i = 0
C        call sumlst(10,nchmx,nbas,ng,ssite,s_site,s_spec,
C     .    mulopt,moddos,nsites,lsites,lmxch,nchan,j,xv,i)
        nlmax = nglob('nlmax')
        if (.not. allocated(orbtm)) allocate(orbtm(nl*nsp*nbas))
        call dpzero(orbtm,nl*nsp*nbas)
      endif

C ... Setup for case sigma or evecs written to disk.
C     No integrated quantities accumulated.  qp read from file
C     Also sigma must be written in (iq,isp) order, opposite to the
C     order in which they would be generated here.
C     Requires second loop over (iq,isp) and filtering of isp
C     First pass should have onesp=0 and lwtkb=0
C     Next line is re-entry point for 2nd spin when writing sigma
   50 continue
      if (lwsig /= 0 .and. lwsig < 10) then
        if (lwsig /= LW5 .and. lwsig /= LW7 .and. lwsig /= LW8) onesp = onesp + 1
        call info5(20,1,0,' BNDFP:  '//
     .    '%?#(n==1)#Write sigm(LDA)##%-1j'//
     .    '%?#(n==2)#Write modified sigm(LDA)##%-1j'//
     .    '%?#(n==3)#Write modified sigm(orb)##%-1j'//
     .    '%?#(n==4)#Write LDA evals,evecs to file##%-1j'//
     .    '%?#(n==5)#Write evals,evecs to file##%-1j'//
     .    '%?#(n==6)#Write evals,evecs,dmat to file##%-1j'//
     .    '%?#(n==7)#Write evecs to ASCII file##%-1j'//
     .    '%?#(n==8)#Read  evecs from evec file##%-1j'//
     .    '%j for %i qp%?#n==2#, spin 2##',
     .    lwsig,nkp,onesp,0,0)
      endif

C   ... debugging : test roth and rotwf
C#ifdefC TESTRWF
C        call trothf(s_bz,s_lat,s_site,s_spec,s_ham,nl,nbas,jsp,nsp,s_ham%offH,
C     .  s_lat%istab,s_bz%qp,ldim,s_pot%smpot,s_pot%smvcnst)
C        call rx('done')
C#endif

C ... Initialize for lmfgwd compatibility mode
      if (lwsig == LW8) then
        ifi = fopna('evec',-1,4); rewind ifi
        call iosigh(3,5,nsp,1,ndham,s_ham%nlmto,nkabc(1),nkabc(2),nkabc(3),nkp,iq,0,0,0,ifi,xv)
      endif

C ... Initialize for k and isp loop (first loop in parallel mode)
C     h,s,evecs are dimensioned ndimh in this loop (may be q dependent)
C     evlk is dimensioned evlk(ndham,2)
C     evals, s_bz%wtkb, s_bz%swtk are dimensioned (ndham,nsp,nkp)
      if (lwtkb==-1 .and. lmet==5) then
        allocate( karchv(nkp) )
        do  iq = 1, nkp
          nullify(karchv(iq)%evals,karchv(iq)%evecs)
        enddo
      endif
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_START_KLOOP,procid,"k-loop")
C#endif
      if (.not. allocated(kpproc)) allocate (kpproc(0:numprocs))
      if (numprocs > 1) then
      if (.not. allocated(karchv)) allocate( karchv(nkp) )
      call info0(30,1,0, ' ... Start MPI k-loop')
      if (numprocs > nkp) call
     .  rxi('MPIK job cannot allocate more processors than nkp =',nkp)
C      Keep track of k-dependent number of evals, for k-dependent APW basis (TK)
C      Future: merge with karchv%nevlk ... information is the same
      if (allocated(nevls)) deallocate(nevls)
      allocate(nevls(nkp,nsp))
      nevls = 0
      sttime = MPI_WTIME()
      call dstrbp(nkp,numprocs,1,kpproc(0))
C      if (lwtkb == -1 .and. lmet==5) then
C        allocate( karchv(kpproc(procid+1)-kpproc(procid)) )
C      endif
      else
        kpproc(0:1) = [1,nkp+1]
      end if
      do  iq = kpproc(procid), kpproc(procid+1)-1
        if (numprocs > 1) then
        if (iq == kpproc(procid)) then
          if (mlog) then
            call gettime(datim)
            call awrit4(' bndfp '//datim//' Process %i of %i on '
     .        //shortname(procid)(1:namelen(procid))//
     .        ' starting k-points %i to %i',' ',256,lgunit(3),
     .        procid,numprocs,kpproc(procid),kpproc(procid+1)-1)
          endif
        endif
        end if

        isqp = nsp*(iq-1)
C       Get qp either from qp list or read from suqlst
        if (numprocs > 1) then
          call dpscop(s_bz%qp,qp,3,3*iq-2,1,1d0)
        else
          if (plbnd == 1) then
            i = nsp; if (onesp /= 0 .or. nspc == 2 .or. lafms) i = 1
            call suqlst(s_lat,plbopt,0,ndhamx,ef0,i,xv,nfbn,ifbls,nkp,qp,onesp)
          elseif (lwden) then
            qp = qpr
          else
            call dpscop(s_bz%qp,qp,3,3*iq-2,1,1d0)
C           qp = [.1d0,.2d0,.3d0]; print *, 'set qp',qp
          endif
        endif

        if (mod(s_ctrl%lfp/2,2) /= 0 .and. .not. cmdopt('--shorten=no',12,0,strn)) then
          xv(1:3) = qp
C         call shorbz(qp,qp,qlat,plat)
          call shorps(1,qlat,(/72,2,2/),xv,qp)
          if (ipr > 50 .and. dlength(3,xv(1:3)-qp,1) > 1d-8)
     .      write (stdo,578) iq,xv(1:3),qp
  578     format(' iq=',i4,'  q=',3f7.3,'  shortened to',3f7.3)
        endif
C       print *, '!!'; qp(1:3) = qp(1:3) + s_lat%qlat(:,2)

C   --- For this qp, G vectors for PW basis and hamiltonian dimension ---
        if (pwemax > 0 .and. mod(pwmode,10) > 0) then
          pwgmin = dsqrt(pwemin); pwgmax = dsqrt(pwemax)
          call dpzero(xv,3); if (mod(pwmode/10,10) == 1) call dpcopy(qp,xv,1,3,1d0)
          call pshpr(iprint()-10)
          call sugvec(s_lat,100000+64+16+2,xv,pwgmin,pwgmax,0,0,napw)
          s_lat%napw = napw
          call poppr
          ndimh = ldim + napw
          if (mod(pwmode,10) == 2) ndimh = napw
          if (ndimh > ndham) then
            call fexit2(-1,111,'%N Exit -1 : BNDFP: ndimh=%i exceeds ndham=%i.'//
     .        '  Try increasing input NPWPAD',ndimh,ndham)
          endif
        else
          ndimh = ldim
        endif
        ndimhx = ndimh*nspc

C   ... Allocate hamiltonian, overlap.  In the SO case, still a loop isp=1..2
        if (iq == 1) then
          i = mpipid(0)
          i = idalloc('h+s',allocvb()+2,2*ndimh,ndimh*3*2*i)
          ltmp = memcheck('bndfp','hamiltonian',s_ctrl%maxmem,T)
        endif
        if (nspc3 > nspc) then
          allocate(h(ndimh,ndimh,3),s(ndimh,ndimh,3))
          call zinit(h,ndimh**2*3); call zinit(s,ndimh**2*3) ! isp=3 holds L+S-
        else
          allocate(h(ndimhx,ndimhx,1),s(ndimhx,ndimhx,1))
          call zinit(h,ndimhx**2); call zinit(s,ndimhx**2)
        endif
        if (ldsdk .and. nspc3 > nspc) then
          call rx('not ready for dsdk this nspc3')
        elseif (ldsdk) then
C          if (nspc == 2) stop 'fix dsdk'
          deallocate(dsdk)
          allocate(dsdk(ndimhx,ndimhx,3,nsp))
          call zinit(dsdk,size(dsdk))
        endif

C   --- Loop over spins ---
        do  isp = 1, nsp
        if (isp == 2 .and. lafms) cycle
        if (onesp == 0 .or. isp == onesp) then
C       if (isp == 1 .or. isp == onesp) call shorbz(qp,qp,qlat,plat)

C   ... Make Hamiltonian and overlap matrices
        nqsig = s_ham%nqsig
C       if (oqsig == 0) oqsig = 1

        ispc = min(isp,nspc3)    ! index to where h,s are stored (hambls)
        i = lrsig*10
        if (ldsdk) i = i+4

C       Special modes:
C       lwsig= LW1:  sigm orbital -> LDA basis
C       lwsig= LW2:  sigm orbital -> LDA basis, high energy parts replaced
C       lwsig= LW3:  sigm orbital -> orbital, high energy parts replaced
C       lwsig= LW4:  Write evecs of LDA hamiltonian to file
C       lwsig= LW5:  Write evecs of hamiltonian to file
C       lwsig= LW6:  Write evecs+dmatk of hamiltonian to file
C       lwsig= LW7:  Write evecs of hamiltonian to ASCII file
C       Note: hambls called in special modes for lwsig=LW1,LW2,LW3
C       hambls returns with evecs,evals of the LDA hamiltonian in s,h
C       Transform sigm to LDA basis: hambls returns sigm(LDA) in h
        if (lwsig == LW1) then
          i = i + 3000
        elseif (lwsig == LW2) then
          i = i + 4000
        elseif (lwsig == LW3) then
          i = i + 5000
        elseif (lwsig == LW4) then
          i = i + 1000
        endif

        if (lwsig==LW8 .and. lmefac<=0) then
          ifi = fopna('evec',-1,4)
          read(ifi) qp, j, nev_from_evec_file
          call alignqp(1,qp,nkp,s_bz%qp,plat,tolq,iv(2),iv(3))
          if (iv(3) /= 0) call rx('mismatch file qp (evec) and given qp')
          if (iv(2) /= iq) call rx('problem with evec qp list')
          call sanrg(.true.,j,ndimhx,ndimhx,'bndfp:','file ham dimsn')
        endif

C ...   First (maybe only) band pass all modes: make evals, evecs
        if (lwtkb == -1 .or. lmet /= 5) then
C       if (iprint()>=55) call pshpr(110)
        modehambls = i
        if (ldsdk) i = i-4
        call hambls(i,s_site,s_spec,s_lat,s_ham,nbas,isp,ispc,qp,
     .      k1,k2,k3,s_ham%qsig,nqsig,s_pot%smpot,s_pot%smvcnst,s_pot%lcplxp,lso,
     .      alfsi,ndimh,h,s,dsdk,ismidb,nevl)
        nsmidb = min(nsmidb,ismidb)

C   ... Mode to write sigma or evecs to file
        if (lwsig > 0 .and. lwsig <= LW6) then
        if (procid == master) then ! Until parallel write is implemented, save only master
C         Write header information
          lshft = s_bz%lshft
          if (iq == 1 .and. isp == 1) then
            if (lwsig <= LW3) then
              ifi = fopna('sigm2',-1,4)
            elseif (lwsig <= LW6) then
              ifi = fopna('evec',-1,4)
C            elseif (ldmatk /= 0) then
C              lwsig = 10
C              if (ldmatk == 2) lwsig = 11
C              ifi = fopna('dmatk',-1,4)
            endif
            rewind ifi
            eseavr = s_ham%eseavr
            call iosigh(0,lwsig,nsp,nspc,ndham,s_ham%nlmto,nk1,nk2,nk3,nkp,nkp,
     .        lshft(1),lshft(2),lshft(3),-ifi,eseavr)
          endif

C         lwsig=LW1,LW2 : dump sigma(LDA basis) into file sigm2
C         lwsig=LW3 : dump sigma(orb basis) into file sigm2
C         lwsig=LW4 : dump LDA eigenvectors into file evec (written after addrbl)
C         lwsig=LW5,LW6 file written after occ number array generated
          if (procid == master) then ! Until parallel write is implemented, save only master
          if (lwsig <= LW6) then
            if (lwsig >= LW4) then
              ifi = fopna('evec',-1,4)
            else
              ifi = fopna('sigm2',-1,4)
            endif
            if (lwsig <= LW4) then
              write(ifi) qp, ndimh, ndimh
              call dpdump(h,ndimh,-ifi)
              call dpdump(s,ndimh**2*2,-ifi)
              cycle
            endif
          endif
          endif
        endif
        endif  ! lwsig modes
        endif  ! lwtkb == -1 .or. lmet /= 5

C   ... Hamiltonian, or evals,evecs now generated for both spins
C       Noncollinear case or lso=3,4 => H made in 2 band passes.
C       For the rest of the isp loop:
C         jsp=isp in the collinear case
C         jsp=1   in the full noncollinear, nspc=2
C         jsp=1,2 in the pert noncollinear, nspc3=2 and nspc=1
C       Thus jsp should be used in place of isp.
C       isp serves as a flag for the noncollinear case
        if (lso == 3 .or. lso == 4) then
          if (isp == 1) goto 30 ! require both spin channels
        elseif (ispc /= nspc) then
          goto 30                 ! require both spin channels
        endif

        if (lwtkb == -1 .or. lmet /= 5) then
C       lso = 1 should also work ... substitutes for sopert branch below
C       if ((lso==3.or.lso==4) .and. i /= -1) then
        if ((lso==1 .or. lso==3 .or. lso==4) .and. i /= -1) then
          if (lso==4) then
            call shorps(1,qlat,(/72,2,2/),qp,xv)
            allocate(hso(ndimh,ndimh,4))
            call mkhso(0,s_site,s_spec,s_lat,s_ham,xv,ndimh,0,hso)
          else
            allocate(hso(1,1,1))
          endif
          call sopert3(10*lso,s_ham,ndimh,h,s,hso)
          deallocate(hso)
          i = -1
          nevl = ndimh
        endif
        endif

C       Re-entry point for loop over jsp=1,2 (nspc3=2 and nspc=1)
        jsp = isp                 ! Evals,etc go into channel jsp
        if (nspc3 == 2) jsp = 1
   20   continue

C   ... Special handling of evals made by hambls
        if (i == -1) then
          lpdiag = 2 ! evals made by hambls
          if (.not. allocated(nevls)) then
            allocate(nevls(nkp,nsp))
            call iinit(nevls,nkp*nsp)
          endif
          nevls(iq,jsp) = nevl
          nevl = ndimh
        endif

C   --- Diagonalize and add to density ---
C       Allocate evecs here, to handle lmet==5
        if (lpdiag == 2) then
          z => s
        else
          allocate(z(ndimhx,ndimhx,1))
        endif

        if (mod(iq,10) /= 1) call pshpr(iprint()-6)
        if (numprocs == 1) call info5(30,0,0,' bndfp:  kpt %i of %i, k=%3:2,5;5d'//
     .    '%?#n#   ndimh = %i##',iq,nkp,qp,mod(pwmode/10,10),ndimh)
        nmx = min(nevmx,ndimhx)
        if (lwtkb == -1 .and. lmet /= 5 .and. plbnd /= 3) nmx = -1

C       Need all eigenvalues if 'fat bands' plotting mode
        if (nfbn(1) > 0) then
          nmx = ndimhx; efmax = 99999
        endif

C   ... Make evals, evecs from hamiltonian
        if (lwtkb == -1 .or. lmet /= 5) then
        if (lpdiag == 1) then
          call rxx(nspc /= 1,'parallel diag not implemented for noncol')
C#ifdefC MPE & MPI
C          ierr = MPE_LOG_EVENT(EVENT_START_PZHEV,procid,"pzhev")
C#endif
          if (iq == 1) call info0(20,0,0,' bndfp:  diagonalise with SCALALPACK ..')
          call rx('update call to pzhev')
C          call pzhev(F,T,T,ndimh,oh,os,nblk,nprow,npcol,efmax,nmx,nev,
C     .               evlk(1,jsp),ot)
          nevl = ndimh
C#ifdefC MPE & MPI
C          ierr = MPE_LOG_EVENT(EVENT_END_PZHEV,procid,"pzhev")
C#endif
        elseif (lpdiag == 2) then
          nev = ndimhx; nmx = ndimhx; nevl = ndimhx
          if (allocated(nevls)) nevls(iq,jsp) = nevl
          call phmbls(2,ndimhx,ndimhx,evlk(1,jsp),xv,xv,xv,xv,xv,h)
          if (lwtkb /= -1) then
            allocate(ww(ndimhx**2*2))
            call blsig(1+lrsig*10,nbas,s_ham,jsp,nsp,nspc,plat,qp,
     .        lwtkb,zval-qbg,iq,s_bz%wtkp,s_bz%wtkb,ndimh,z,ww)
            deallocate(ww)
          endif
        else
          if (nspc == 2) then
            allocate(ww(2*ndimhx**2))
            call sopert(0,ndimh,nspc,ww,h,h)
            call sopert(0,ndimh,nspc,ww,s,s)
            deallocate(ww)
          endif
          allocate(ww(11*ndimhx))
C          if (iprint() >= 100) then
C            call zprm('h',12,h,ndimhx,ndimhx,ndimhx)
C            call zprm('s',12,s,ndimhx,ndimhx,ndimhx)
C          endif
          if (epsovl == 0 .and. ovncut == 0) then
C           call zhev(ndimhx,h,s,T,T,nmx,efmax,nev,ww,F,-1,evlk(1,jsp),z)
            call zhevx(ndimhx,ndimhx,h,s,1,T,nmx,efmax,nev,ww,F,evlk(1,jsp),ndimhx,z)
            nevl = ndimhx
            if (nmx > 0) nevl = nev  ! if evecs also calculated
C           call zprm('z',12,z,ndimhx,ndimhx,nev)
          else
            nevl = -1
            call dvset(ww,1,1,99999d0)
            call zhevo(ndimhx,ndimhx,h,s,nmx,efmax,epsovl,ovncut,nevl,nev,evlk(1,jsp),ww,ndimhx,z)
            eomin = ww(1)
            mxevcut = max(mxevcut,ndimhx-nevl)
            if (mnevl == NULLI) mnevl = nevl
            mnevl = min(nevl,mnevl)
          endif
          deallocate(ww)
        endif

C       Limit size of very large evals (PMT method can make nearly singular overlap)
        if (ndham > s_ham%nlmto) then
          do  i  = 1, ndimhx
C           evlk(i,jsp) = min(199d0,evlk(i,jsp))  ! may violate compiler bounds check
            call dvset(evlk(1,jsp),i,i,min(999d0,dval(evlk(1,jsp),i)))
          enddo
        endif

C   ... Read evals, evecs from evec
        if (lwsig==LW8 .and. lmefac<=0) then
          ifi = fopna('evec',-1,4)
          call dpdump(evals(1,jsp,iq),nev_from_evec_file,ifi)
          call dpdump(z,ndimhx**2*2,ifi)
          xv(1) = dlength(nevl,evals(1:nevl,jsp,iq)-evlk(1:nevl,jsp),1)
          call info2(10,0,0,' bndfp: difference in lmf and file(evec) evals for iq=%i : %;3g',iq,xv)
        endif

C       Keep track of lowest, highest eval for each ib
        if (llohi) then
          forall (i = 1:ndimhx) elohi(i,jsp,1) = min(elohi(i,jsp,1),evlk(i,jsp))
          forall (i = 1:ndimhx) elohi(i,jsp,2) = max(elohi(i,jsp,2),evlk(i,jsp))
        endif

C   ... lmet==5 branch, first pass
        if (lmet==5) then
          if (jsp == 1) then
C            if (associated(karchv(iq)%evals)) then
C              deallocate(karchv(iq)%evals,karchv(iq)%evecs)
C            endif
            allocate( karchv(iq)%evals(ndham,2) )
            allocate( karchv(iq)%evecs(ndimhx**2,nspx) )
          endif

          karchv(iq)%nevk(jsp) = nev
          karchv(iq)%nevlk(jsp) = nevl
          karchv(iq)%eomink(jsp) = eomin
          call dcopy(ndhamx,evlk(1,jsp),1,karchv(iq)%evals(1,jsp),1)
          call zcopy(ndimhx**2,z,1,karchv(iq)%evecs(1,jsp),1)
        endif

C   ... lmet==5 branch, second pass
        else
          if (.not. allocated(karchv)) call rx
     .      ('metal=5 not allowed in this branch ... arrays not allocated')
          nev = karchv(iq)%nevk(jsp)
          nevl = karchv(iq)%nevlk(jsp)
          eomin = karchv(iq)%eomink(jsp)
          call dcopy(ndhamx,karchv(iq)%evals(1,jsp),1,evlk(1,jsp),1)
          call zcopy(ndimhx**2,karchv(iq)%evecs(1,jsp),1,z,1)

          if (lswtk /= 1 .or. lwtkb == 2) then  ! Skip if preserve for spin wts
            if (jsp == nspx) then
              deallocate(karchv(iq)%evals,karchv(iq)%evecs)
            endif
          endif

C     ... Maintain evecs file pointer
          if (lwsig==LW8 .and. lmefac<=0) then
            ifi = fopna('evec',-1,4)
            read(ifi); read(ifi)
          endif

        endif

C   ... Fixed spin moment : shift potential by fsmom(2)
        if (fsmom(2) /= 0) then
C         Not ready for noncoll case
          if (nspc /= 2) then
            if (.not. associated(s_bz%swtk)) s_bz%swtk => xx
            call bzbfde(lswtk==1,1,1,ndham,ndham,(3-2*jsp)*fsmom(2),s_bz%swtk,evlk(1,jsp))
          endif
        endif

        if (ldsdk .and. lwtkb >= 0) then
          call hambldsdk(modehambls,s_site,s_spec,s_lat,s_ham,nbas,isp,ispc,qp,
     .    k1,k2,k3,s_ham%qsig,nqsig,s_pot%smpot,s_pot%smvcnst,s_pot%lcplxp,lso,
     .    alfsi,ndimh,z,nev,dsdk)
        endif

C       Pad evals between ndimh and ndham with a large positive number
C       to avoid mixing up integration routines
        if (ndhamx > nevl .and. nspc == 2) then
          call dvset(evlk,1+nevl,ndhamx,99999d0)
        elseif (ndhamx > nevl) then
          call dvset(evlk(1,jsp),1+nevl,ndham,99999d0)
        endif
        if (allocated(nevls) .and. lpdiag /= 2) then
          nevls(iq,jsp) = nevl
        endif

C       AFM condition: eigenvalues duplicated
        if (lafms) then
          if (fsmom(1) /= NULLR .and. fsmom(1) /= 0)
     .      call rx('Special AFM symmetry incompatible w/ FSMOM')
          evlk(:,2) = evlk(:,1)
          if (allocated(nevls)) nevls(iq,2) = nevls(iq,1)
        endif

        if (numprocs == 1) then
        if (epsovl /= 0 .or. ovncut /= 0) then
          if (lpdiag /= 2) then
            call info5(30,0,0,' Overlap''s smallest eigenvalue: %;3g.  '//
     .        '%?#(n>0)#H dim reduced from %i to %i#H dim not reduced#',
     .        eomin,ndimhx-nevl,ndimhx,nevl,0)
          elseif (nevls(iq,jsp) /= ndimhx) then
            call info2(30,0,0,' ham dimension reduced to %i',nevls(iq,jsp),0)
          endif
        endif
        call prtev(z,nevl,evlk(1,jsp),nmx,efmax,nev)
        if (iprint() >= 110 .or. .false.) then
          call yprm('evals',1,evlk(1,jsp),1,ndhamx,nevl,1)
          if (nev > 0) call zprm('evecs',2,z,ndimhx,ndimhx,nev)
        endif
        end if

        if (mod(iq,10) /= 1) call poppr

C ...  In k-parallel mode, defer this section until all qp available
        if (numprocs > 1) then
        if (plbnd == 0 .or. plbnd == 2) then
          if (lwtkb /= -1 .and. .not.lwndow) then
            ef00 = ef0
            if (iq == kpproc(procid) .and. jsp == nsp .and.
     .        .not. cmdopt('--no-fixef0',11,0,strn)) then
              if (procid == master) then
                call fixef0(zval-qbg,jsp,1,nevl,ndham,evlk,dosw,ef0)
              endif
              call mpibc1( ef0, 1, 4, mlog, 'bndfp', 'ef0' )
              call mpibc1( dosw, 2, 4, mlog, 'bndfp', 'dosw' )
              if (jsp == 2 .and. ef00 /= ef0 .and.
     .            lwtkb == 0 .and. lmet > 0 .and. lrout /= 0) then
                call info0(10,1,1,' ... Fermi level reset in second'//
     .            '  spin channel ... restart band pass')
                deallocate(h,s)
                if (iq == 1)
     .            i = idalloc('h+s',allocvb()+4,0,0)
                goto 99
              endif
            endif
          endif
        endif
        else
        ebot = dmin1(ebot,evlk(1,jsp))
        i = max(1,nint(zval-qbg)/(3-nspc))
        evtop = max(evtop,evlk(i,jsp))
        ecbot = min(ecbot,dval(evlk(1,jsp),i+1))
        if (lmetsav == NULLI) lmetsav = lmet
        if (lmetsav == 0 .and. iq == 1 .and. jsp == 1) ef0 = evtop
        if (plbnd == 0 .or. plbnd == 2) then
        if (ipr>=10 .and. iq==1 .and. ipl>1)
     .  write (stdl,"('fp evl',9f8.4)") ((evlk(i,j),i=1,nev/nspc),j=1,nspc)
        if (lwtkb /= -1 .and. .not.lwndow .and. .not.lwden) then
          if (iq == 1 .and. jsp == nsp .and. .not. cmdopt('--no-fixef0',11,0,strn)) then
            ef00 = ef0
            call fixef0(zval-qbg,jsp,1,nevl,ndham,evlk,dosw,ef0)
            if (jsp == 2 .and. ef00 /= ef0 .and.
     .        lwtkb == 0 .and. lmet > 0 .and. lrout /= 0) then
              if (procid == master) call info0(10,1,1,
     .          ' ... Fermi level reset in second spin channel ... restart band pass')
              if (iq == 1) i = idalloc('h+s',allocvb()+4,0,0)
              deallocate(h,s)
              goto 99
            endif
          endif
C         Check for cases when nevmx is too small : i=2 => fatal error
          i = 0
          if (nevmx>=0 .and. lmet /= 0) then
            dum = evlk(max(nev/nspc,1),jsp*nspc)
C           if (ef0 >= dum) i = 2
            if (ltet /= 1 .and. ef0+5*dabs(esmear-mpsord) > dum) i=2
            if (lmet==4 .and. ef0+def+5*dabs(esmear-mpsord)>dum)i=2
          endif
          if (i == 2) then
            call info5(2,1,0,' evl(nev=%i)=%;3d but ef0=%;3d ... restart with '//
     .        'larger efmax or nevmx',nev,evlk(max(nev,1),jsp),ef0,0,0)
            call rx('bndfp')
          endif
        endif
        endif
        end if

C   ... Remake local ocrng,unrng for this k.  Only if ef0 is known
        if (loptic /= 0 .and. lwtkb > 0) then
          xv(1) = ef0
          if (s_bz%egap /= NULLI) xv(1) = ef0 + s_bz%egap/2
          call ocnock(s_optic,iq,evlk(1,jsp),nevl,xv(1),s_optic%esmr)
        endif

C   ... Save evals for this qp
        if (plbnd /= 1) then
          if (lwtkb == -1 .or. lmet /= 5) then
C     ... Copy eigenvalues into array containing ev for all qp
          call dpscop(evlk(1,jsp),evals,ndhamx,1,1+ndham*(jsp-1+isqp),1d0)
          if (lafms) call dpscop(evlk(1,2),evals,ndhamx,1,1+ndham*(2-1+isqp),1d0)

C     ... Save evecs for this qp
          if (lwsig == LW5) then
C#ifdefC MPI
C            if (procid == master) then
C#endif
            ifi = fopna('evec',-1,4)
C            print 777, 'write evec q,evlk,z',procid,iq,0,qp,evlk(21,1),z(1,1,1)
C  777       format(a,3i5,10f14.7)
            write(ifi) qp, ndimhx, nev
            call dpdump(evlk(1,jsp),nev,-ifi)
            call dpdump(z,ndimhx**2*2,-ifi)
C#ifdefC MPI
C            endif
C#endif
          else if (lwsig == LW7) then
            ifi = fopna('eveca',-1,0)
            call ywrm(0,'evecs',3,ifi,'(9f20.10)',z,1,ndimhx,ndimhx,nev)
          endif

          if (numprocs > 1) then
          if (epsovl /= 0 .or. ovncut /= 0) then
            call dpscop(eomin,eomink,1,1,1+(jsp-1+isqp),1d0)
            if (lafms) then
              call dpscop(eomin,eomink,1,1,1+(2-1+isqp),1d0)
            endif
          endif
C          karchv(iq)%nevk(jsp) = nev
C          karchv(iq)%nevlk(jsp) = nevl
C          karchv(iq)%eomink(jsp) = eomin
          if (lafms) karchv(iq)%eomink(2) = eomin
          end if
          endif
C        call prmx('ev',evals,ndham,ndham,isqp+2)
C        call prmx('ev',evals,ndham,ndham,nkp*nsp)

C     --- Orbital magnetic moment (requires lso) ---
          if (lso /= 0 .and. lwtkb /= -1 .and. .not. lafms) then
            if (lwtkb == 0 .and. iq == 1) then
              call info0(20,0,0,' (warning) '//
     .          'skip orb. moment calculation ... metal weights required')
            else
              allocate(ausp(nlmax,ndhamx,nphimx,nsp,nbas)); call dpzero(ausp,2*size(ausp))
              call makusq(10,s_site,s_spec,s_lat,s_ham,nbas,nbas,0,nlmax,1,ndham,nphimx,
     .          ndimh,napw,s_lat%igv2,nev,nsp,nspc,nsp,jsp,1,qp,z,ppnl,ausp,xv)
              call mkorbm(s_site,s_spec,jsp,nsp,nspc,nlmax,ndham,nphimx,nev,
     .          s_bz%wtkb,iq,nbas,ppnl,ausp,nl,nkp,orbtm)
              deallocate(ausp)
            endif
          endif

C     --- Mulliken analysis and partial DOS for 1 irreducible qp and star of qp ---
          if (lwtkb /= -1) then
          if (procid == master .or. cmdopt('--mull',6,0,strn) .or. cmdopt('--pdos',6,0,strn)) then
C           call zprm('evecs',2,z,ndimhx,ndimhx,nev)
            call bzmsh00(s_lat%plat,s_bz%lshft,0,s_bz%nkabc,ifac,qb) ! Needed for symmetrization
            if (cmdopt('--mull',6,0,strn) .or. cmdopt('--pdos',6,0,strn)) then
              call rxx(lafms,'pdos not allowed in conjunction w/ AFM')
              ltmp = cmdopt('--pdos',6,0,strn)
              if (ltmp) then  ! --pdos
                i1 = idalloc('aus',allocvb()+2,nlmax*ndhamx*nkaph*nsp*nbas,2)
                ltmp = memcheck('bndfp','aus',s_ctrl%maxmem,T)
                allocate(ausp(nlmax,ndhamx,nphimx,nsp,nbas)); call dpzero(ausp,2*size(ausp))
              endif
              allocate(doswt(nchan,ndimhx,nspc)); call dpzero(doswt,size(doswt))

C             Debugging : choose phase of each so that phase of 1st orbital is zero
C              do  i = 1, nev
C                xv(1) = datan2(dimag(z(1,i,1)),dble(z(1,i,1)))
C                z(:,i,1) = z(:,i,1) * dcmplx(dcos(xv(1)),-dsin(xv(1)))
C              enddo

C         ... Setup for loop over star of k
              i1 = 0            ! Flags iqstar to start a new loop over the star of qp
              dum = 1           ! reduce weighting by the number of points in star
              qpr = qp          ! the rotated qp is the irreducible qp, unless it gets rotated
              do  while (.true.) ! loop over star of qp
                zr => z         !The rotated eigenvector is initially the unrotated one
                if (moddos == 2 .or. moddos == 5) then ! Unless m-resolved DOS, no need to rotate
                  if (i1 == 0) then ! First point: count number of points in star
                    call iqstar(-1,iq,nk1,nk2,nk3,s_bz%ipq,ifac,qb,i,qpr,[0])
                    dum = 1/dble(i) ! reduce doswt by number of points in star of iq
C                   dum = 1       ! Debugging use just 1st point in star

                  endif
                  call iqstar(1,iq,nk1,nk2,nk3,s_bz%ipq,ifac,qb,i1,qpr,[0]) ! Get next (i1,qpr) in star
                  if (i1 == 0) exit  ! No more qp in star
                  ig = s_bz%star(i1+1) ! group operation mapping qp to qpr
                  if (ig > 1) then  ! If not the identity operation, rotate z to zr
C                   if (ig /= 1) exit ! Debugging use just 1st point in star
                    zr => h            ! h is used as a help array here
                    allocate(rmat(nl**4*2))
                    i = 140
C                   Debugging : return rotation matrix R in zr rather than zr=R*z
C                    if (.false.) then
C                      i = 1140 ; zr = 0
C                    endif
                    call rotwf(i,nl,nbas,nspc,s_lat%pos,s_ham%offH,s_ham%iprmb,s_lat%istab(1,ig),
     .                s_lat%symgr(1,ig),s_lat%ag(1,ig),qpr,rmat,0,ldim,ldim,nev,z,zr)
                    deallocate(rmat)
                  endif
C                  i = 0
C                  if (i /= 0) call zprm('evec(r)',2,zr,ndimhx,ndimhx,nev)
                endif
                if (ltmp) then  ! --pdos
C                 Make (u,s) for (qpr, zpr)
                  call dpzero(ausp,2*size(ausp))
                  call makusq(1,s_site,s_spec,s_lat,s_ham,nbas,nbas,0,nlmax,1,ndham,nphimx,
     .              ndimh,napw,s_lat%igv2,nev,nsp,nspc,nsp,jsp,1,qpr,zr,ppnl,ausp,qpr)
C                 Accumulate partial DOS for (qpr, zpr)
                  call mkpdos(moddos,s_site,s_spec,jsp,nsp,nspc,nlmax,ndham,nphimx,nev,
     .              nchan,idoschan,lmdim,1,lsites,nsites,ppnl,ausp,doswt)
                else  ! --mull
                  call mullmf(nbas,s_site,s_spec,s_ham%iprmb,zr,ndimh,nspc,iq,jsp,
     .              moddos,nsites,lsites,nchan,idoschan,lmdim,ndham,doswt)
                endif
                if (.not. (moddos == 2 .or. moddos == 5)) exit ! No star unless m-resolved DOS
              enddo
              if (dum /= 1) call dscal(size(doswt),dum,doswt,1)
              nschan = mod(nfstg/10,10)
              if (ltmp) then  ! --pdos
                i1 = iomoms(-nfilem,nl,nsp,nspc,2,ndimh,nfstg,nschan,1,1,ndhamx,ndimhx,
     .            nchan,nchan,nev,evlk(1,jsp),0d0,doswt,0d0,0d0,0d0)
                deallocate(ausp)
                i1 = idalloc('aus',allocvb()+4,0,0)
              else  ! --mull
                i1 = iomoms(-nfilem,nl,nsp,nspc,2,ndimh,nfstg,nspc,1,1,ndhamx,
     .            ndimhx,nchan,nchan,nev,evlk(1,jsp),0d0,doswt,0d0,0d0,0d0)
              endif
              if (allocated(doswt)) deallocate(doswt)
            elseif (nfilem > 0) then
              write (nfilem) 0, ndimhx
              call dpdump(evlk(1,jsp),ndimhx,-nfilem)
              if (lafms) then
                write (nfilem) 0, ndimhx
                call dpdump(evlk(1,2),ndimhx,-nfilem)
              endif
            endif
           endif ! procid == master

C      .. Calculate matrix elements of gradient operator for eps
          if (mod(loptic,10)>0) then
C           call lkdim(1,nbas,s_site,s_spec,i,nlmax)
C           nlmax = (nlmax+1)**2
            nlmax = s_optic%nlg**2
C           nlmax = nglob('mxorb') / nkaph
            allocate(ausp(nlmax,ndhamx,nphimx,nsp,nbas)); call dpzero(ausp,2*size(ausp))
            k = s_optic%kmxax+1+nkaph
            i = nlmax*ndham*nspc*k*nbas*nspc
            allocate(cPkLq(i)); call dpzero(cPkLq,2*i) ! Coff to PkL exp
            ! Value, slope of phi,phidot at MT boundary
            call makusq(101+10*morder,s_site,s_spec,s_lat,s_ham,nbas,nbas,0,nlmax,k,ndham,
     .        nphimx,ndimh,napw,s_lat%igv2,nev,nsp,nspc,nspc,jsp,1,qp,z,ppnl,ausp,cPkLq)
            ! Matrix elements of gradient operator
            call opt_nstate(s_optic,3,nevl,nfiloe,nfiupe,nemloe,nemupe,i1,i2)
            allocate(optme(nfiloe:nfiupe,nemloe:nemupe,3)); call dpzero(optme,2*size(optme))
            s_optic%optme => optme
            call fpopm(s_ctrl,s_spec,s_lat,s_optic,nlmax,k,ndham,nphimx,jsp,
     .        nsp,nspc,nbas,nevl,ausp,cPkLq,10*isw(ldsdk)+morder,dsdk)
C           call zprm('optme(1)-optme(2)',2,optme,nfilm,nfilm,nempm)
            deallocate(ausp,cPkLq)
          endif

C   ...   Make new density matrix dmatu for LDA+U
          if (nlibu > 0 .and. nev > 0 .and. .not.lwden) then
            if (lwtkb == 0) call rx('metal weights required for LDA+U calculation')
            nl = s_ctrl%nl; nlmax = nl*nl
            allocate(ausu(nlmax,ndhamx,nphimx,nsp,nbas)); call dpzero(ausu,2*size(ausu))
            call makusq(0,s_site,s_spec,s_lat,s_ham,nbas,nbas,0,nlmax,1,ndham,nphimx,
     .        ndimh,napw,s_lat%igv2,nev,nsp,nspc,nsp,jsp,1,qp,z,ppnl,ausu,xv)
            call mkdmtu(s_site,s_spec,s_bz%wtkb,jsp,iq,nsp,nspc,nlmax,ndham,nphimx,
     .        nbas,nev,ppnl,ausu,dmatu,nlibu,lmaxu,lldau)
            deallocate(ausu)
          endif
          endif ! post processing lwtkb /= -1

C     ... Accumulate ausc for core-level spectroscopy
          if ((icls /= 0 .or. irxes /= 0) .and. lwtkb /= -1) then
            call rxx(nspc /= 1,'CLS not implemented in noncoll case')
C#ifdefC MPI
C            call rx('CLS only k-parallel')
C#endif
            call makusq(0,s_site,s_spec,s_lat,s_ham,nbas,nsites,isite,nlmax,1,ndham,
     .        nphimx,ndimh,napw,s_lat%igv2,nev,nsp,nspc,nsp,jsp,iq,qp,z,ppnl,ausc,xv)

C     --- Accumulate output density ---
C         Even if no output rho, still call addrbl to make DOS when lmet=4
          elseif (lwtkb /= -1 .and. (lrout /= 0 .or. lmet==4)) then
C           dum=frc(3,1,1)
            if (ldmatk /= 0) then
              allocate(dmatk(ndhamx,1))
            else
              allocate(dmatk(1,1))
            endif
            if (.not. associated(s_bz%swtk))  s_bz%swtk => xx
            call addrbl(s_site,s_spec,s_lat,s_ham,s_bz,s_optic,jsp,nsp,nspc,qp,
     .        ndham,ndimh,lmet,lrout,lwtkb,ldmatk,ltso,s_bz%wtkb,lswtk,s_bz%swtk,
     .        iq,lfrce,ldos,lekkl,k1,k2,k3,s_pot%smpot,s_pot%smvcnst,s_pot%lcplxp,
     .        s_bz%numq,zval-qbg,z,evlk,nevl,ef0,def,esmear,emin,emax,ndos,dos,
     .        srout,sumqv,sumev,frc,dmatk,tso)
            if (mod(iq,10) /= 1) call pshpr(iprint()-6)
            i = 1; if (iprint() >= 30) i = 3
            call fpopm_sfac(i,s_lat,s_optic,ldsdk,dsdk,iq,ndham,isp,nsp,nspc,nbas,nev)
            if (mod(iq,10) /= 1) call poppr
C           Replace evec with translated form, for AFM.  Destroys h
            if (lafms) then
              call dcopy(ndhamx**2*2,z,1,h,1)  ! copy evecs to work array
              npgrp = s_lat%npgrp
              allocate(wk(3*nbas))
              call sitepack(s_site,1,nbas,'pos',3,xv,wk)
              call afmevc(1,nbas,npgrp,s_lat%istab,s_ham%iprmb,qp,plat,
     .          wk,s_lat%igv2,afmt,ndimh,napw,nev,h,z)
              deallocate(wk)
              if (nlibu > 0 .and. nev > 0 .and. .not.lwden) then
                if (lwtkb == 0) call rx('metal weights required for LDA+U calculation')
                nl = s_ctrl%nl; nlmax = nl*nl
                call rx('bndfp 2348')
                allocate(ausu(nlmax,ndhamx,nphimx,nsp,nbas)); call dpzero(ausu,2*size(ausu))
                call makusq(0,s_site,s_spec,s_lat,s_ham,nbas,nbas,0,nlmax,1,ndham,nphimx,
     .            ndimh,napw,s_lat%igv2,nev,nsp,nspc,nsp,2,1,qp,z,ppnl,ausu,xv)
                call rx('bndfp 2347')
                call mkdmtu(s_site,s_spec,s_bz%wtkb,2,iq,nsp,nspc,nlmax,ndham,nphimx,
     .            nbas,nev,ppnl,ausu,dmatu,nlibu,lmaxu,lldau)
                deallocate(ausu)
              endif
              call addrbl(s_site,s_spec,s_lat,s_ham,s_bz,s_optic,2,nsp,nspc,qp,
     .          ndham,ndimh,lmet,lrout,lwtkb,ldmatk,ltso,s_bz%wtkb,lswtk,s_bz%swtk,
     .          iq,lfrce,ldos,lekkl,k1,k2,k3,s_pot%smpot,s_pot%smvcnst,s_pot%lcplxp,
     .          s_bz%numq,zval-qbg,z,evlk,nevl,ef0,def,esmear,emin,emax,ndos,dos,
     .          srout,sumqv,sumev,frc,dmatk,tso)
            endif
C           Return symmetrized |optme|^2 in optmt and possibly optme in optmc
            if (mod(loptic,10)>0) then
              call opt_nstate(s_optic,2,nevl,i,i,i,i,i1,i2)
              i = 1; if (lwoptmc /= 0) i = 5
              call symdme(i,s_lat%symgr,s_lat%nsgrp,nfiloe,nfiupe,nemloe,nemupe,
     .          nfilo,nfiup,nemlo,nemup,1d0,optme,optmt(1,nfilo,nemlo,jsp,iq),
     .          optmc(1,nfilo,nemlo,jsp,iq))
              deallocate(optme)
            endif
            if (lwsig == LW6) then
              ifi = fopna('evec',-1,4)
              write(ifi) qp
              call dpdump(evlk(1,jsp),ndimh,-ifi)
              call dpdump(dmatk(1,jsp),ndimh,-ifi)
              call dpdump(z,ndimh**2*2,-ifi)
            endif ! lwsig = LW6
            deallocate(dmatk)

C            print 399,iq,frc(3,1,1),frc(3,1,1)-dum
C  399       format(' after addrbl: frc(3,1)=',i4,2f12.6)
          endif

        elseif (plbnd == 1) then
          if (ndimhx /= nevl) call rx('color weights not implemented when nevl < ham. dimension')
          if (numprocs > 1) then
            call dpscop(evlk(1,jsp),evals,ndhamx,1,1+ndham*(jsp-1+isqp),1d0) ! keep evals for this qp
            if (nfbn(1) /= 0) then
              i = nsp; if (onesp /= 0 .or. nspc == 2 .or. lafms) i = 1
              do  j = 1, 3
                k = j ; if (spintexture) k = 0
                if (nfbn(j) > 0)
     .            call suqlse(1,ndimhx,jsp,i,ndimhx,k,nfbn,ifbls,ndhamx,ndham*nsp,z,colwt(1,jsp,j,iq))
C               print *, 'procid',procid,jsp,j,iq,colwt(1,jsp,j,iq)
              enddo
            endif
            if (epsovl /= 0 .or. ovncut /= 0) then
              call dpscop(eomin,eomink,1,1,1+(jsp-1+isqp),1d0)
            endif
          else
            if (lwsig == LW7) then
              ifi = fopna('eveca',-1,0)
              call ywrm(0,'evecs',3,ifi,'(9f20.10)',z,1,ndimhx,ndimhx,nev)
            endif
            i = nsp; if (onesp /= 0 .or. nspc == 2 .or. lafms) i = 1
            k = jsp; if (onesp /= 0) k = 1
            call suqlsw(nevl,k,i,evlk(1,jsp)) ! Write the energy bands to file for this qp
            if (nfbn(1) /= 0) then
              if (ndimhx /= nevl) call rx('color weights not implemented when nevl < ham. dimension')
              do  j = 1, 3
                k = j ; if (spintexture) k = 0
                if (nfbn(j) > 0) call suqlse(0,ndimhx,jsp,i,ndimhx,k,nfbn,ifbls,ndhamx,1,z,z)
              enddo
            endif
          endif
        endif

C   ... Color weights to project jdos
        if (loptic<0 .and. lwtkb /= -1 .and. maxval(nfbn)>0) then
          allocate(wk3(ndimhx))
          do  i = 1, njdosw
            allocate(iwk1(maxval(nfbn)))
            allocate(iwk2(maxval(nfbn)))
            call getjdosw(3,1,sjdosw(i),nfbn(i),nfbn(i+2),iwk1,iwk2)
            if (nfbn(i) /= 0) then    ! initial states, or just states
C             call info2(0,0,0,' weight list %n:1i',nfbn(i),iwk1)
              call suqlsz(ndimhx,1,nfbn(i),iwk1,nfbn(i),z,wk3)
              call dcopy(nfilm,wk3(nfilo),1,jdosw(i,1,iq,jsp),njdosw)
            else
              jdosw(i,1:nfilm,iq,jsp) = 1
            endif
            if (nfbn(i+2) /= 0) then  ! final states, if they exist
C             call info2(0,0,0,' weight list %n:1i',nfbn(i+2),iwk2)
              call suqlsz(ndimhx,1,nfbn(i+2),iwk2,1,z,wk3)
              call dcopy(nempm,wk3(nemlo),1,jdosw(i,1+nfilm,iq,jsp),njdosw)
            else
              jdosw(i,1+nfilm:nfilm+nempm,iq,jsp) = 1
            endif
            deallocate(iwk1,iwk2)
          enddo
          deallocate(wk3)
        endif

        if (.not. associated(z,target=s)) deallocate(z)

        if (jsp == 1 .and. nspc == 1 .and. nspc3 == 2) then
          jsp = 2
          call zcopy(ndimh**2,h(1,1,2),1,h,1)
          call zcopy(ndimh**2,s(1,1,2),1,s,1)
          goto 20
        endif

   30   continue
        endif
        enddo  !end loop over isp (main loop in parallel mode)
        if (iq == 1) i = idalloc('h+s',allocvb()+4,0,0)
        deallocate(h,s)

      enddo  ! end loop over iq (main loop in parallel mode)

C ... Second pass over iq for second spin when writing sigma
      if (lwsig /= 0) then
        if (onesp > 0 .and. onesp < nsp) goto 50
        if (lwsig == LW1 .or. lwsig == LW2 .or. lwsig == LW3) then
          if (procid == master) then ! Until parallel write is implemented, save only master
          ifi = fopna('sigm2',-1,4)
          call fclose(ifi)
          endif
          if (lwsig == LW3)
     .    call rx0('BNDFP:  sigm(orb basis) saved in file sigm2')
          call rx0('BNDFP:  sigm(LDA basis) saved in file sigm2')
        else if (lwsig >= LW4 .and. lwsig <= LW6) then
          if (procid == master) then ! Until parallel write is implemented, save only master
          ifi = fopna('evec',-1,4)
          call fclose(ifi)
          endif
          if (lwsig == 5) then
            call info0(20,0,0,' BNDFP:  evecs saved in file evec')
          else
            call rx0('BNDFP:  evecs saved in file evec')
          endif
        endif
C        if (lwsig == -1) then
C        ifi = fopna('sigm3',-1,4)
C        call fclose(ifi)
C        call rx0('BNDFP:  sigm(orbital basis) saved in file sigm3')
C        endif
C        if (lwsig == LW4) then
C        call rx0('BNDFP:  U(LDA-QP) saved in file sigm2')
C        endif
      endif

C ... Restore q=0 lattice vectors
      call pshpr(10)
      call sugvec0(s_lat)
      call poppr

C ... Set switch to write file sigii
C#ifndef LMFGWD | LMFGW
      if (lrsig /= 0 .and. plbnd == 0 .and. procid == master) then
        sigp = s_ham%sigp
        if (mpipid(0) == 1) call phmbl3(12,0,nsmidb,0,0,sigp,qp,qp(2))
        call fclose(fopna('sigii',-1,0))
      endif
C#endif

C --- Second loop over qp, needed to make k-parallelisation possible: ---
C     You can't do this until you have all the evals.
      if (numprocs > 1) then
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_END_KLOOP,procid,"k-loop")
C#endif
      entime = mpi_wtime()
      call info2(30,0,-1, ' ... Done MPI k-loop: %;1ds elapsed at root...',entime-sttime,0)
!     flush (6)   ! Generates error on intel v11 compiler
      call mpi_barrier(mpi_comm_world, ierr)
      call info2(30,0,0, ' %;1ds in all',mpi_wtime()-sttime,0)

      call info0(20,0,-1,' ... Sharing data between processes...')
      sttime = mpi_wtime()

C     Share bands among all processors
      call xmpbnd(kpproc,ndham,nkp,nsp,evals)
C     Share color weights among all processors
      if (plbnd == 1 .and. nfbn(1) /= 0) then
        call xmpbnd(kpproc,ndham*nsp,nkp,4,colwt)
      endif

      if (epsovl /= 0 .or. ovncut /= 0) then
        call xmpbnd(kpproc,1,nkp,nsp,eomink)
      endif
      if (icls /= 0 .and. lwtkb /= -1) then
        if (nsites /= 1) call rx('CLS: need modify xmpbnd for MPI if nsite>1')
        call rxx(nspc /= 1,'CLS not implemented in noncoll case')
        call xmpbnd(kpproc,(2*size(ausc))/(nkp*nsp),nkp,nsp,ausc)
      endif
C     Share number of eigenvalues for each band
      call mpibc2(nevls,nkp*nsp,2,3,mlog,'bndfp','nevls')
      deallocate(kpproc, stat=ierr)
      if (llohi) then
        do  isp = 1, nspx
          call mpibc2(elohi(1,isp,1),ndimhx,4,2,.false.,'','')
          call mpibc2(elohi(1,isp,2),ndimhx,4,1,.false.,'','')
        enddo
      endif

C     Allreduce density-related quantities
      if (lrout /= 0) then
      call mpibc2(sumqv,6,4,3,mlog,'bndfp','sumqv')
      call mpibc2(sumev,6,4,3,mlog,'bndfp','sumev')
      call mpibc2(srout,k1*k2*k3*nsp*s_bz%numq,6,3,mlog,'bndfp','smrho')
      if (lswtk == 1) then
        call mpibc2(s_bz%swtk,ndhamx*nkp,4,3,mlog,'bndfp','swtk')
      endif
C     Allreduce qkkl
      do  ib = 1, nbas
        is = s_site(ib)%spec
        lmxa = s_spec(is)%lmxa
        lmxh = s_spec(is)%lmxb
        kmax = s_spec(is)%kmxt
        nlma = (lmxa+1)**2
        nlmh = (lmxh+1)**2
        if (lmxa > -1) then
          nelt(1) = (kmax+1)*(kmax+1)*nlma*nlma
          nelt(2) = (kmax+1)*nkaph*nlma*nlmh
          nelt(3) = nkaph*nkaph*nlmh*nlmh
          i = s_bz%numq*nsp*nspc
          j = s_bz%numq*nsp*(2*nspc-1)
          call mpibc2(s_site(ib)%qkkl,nelt(1)*i,4,3,F,' ',' ')
          call mpibc2(s_site(ib)%qhkl,nelt(2)*j,4,3,F,' ',' ')
          call mpibc2(s_site(ib)%qhhl,nelt(3)*i,4,3,F,' ',' ')
          if (lekkl == 1) then
            call mpibc2(s_site(ib)%eqkkl,nelt(1)*i,4,3,F,' ',' ')
            call mpibc2(s_site(ib)%eqhkl,nelt(2)*j,4,3,F,' ',' ')
            call mpibc2(s_site(ib)%eqhhl,nelt(3)*i,4,3,F,' ',' ')
          endif
        endif
      enddo
C     Allreduce DOS, forces, dmatu, orbtm
      if (ndos > 0) call mpibc2(dos,ndos*2*nsp,4,3,mlog,'bndfp','dos')
      if (lfrce /= 0) call mpibc2(frc,3*nbas*s_bz%numq,4,3,mlog,'bndfp','frc')
      if (nlibu > 0) call mpibc2(dmatu,nsp*nlibu*(lmaxu*2+1)**2,6,3,mlog,'bndfp','dmatu')
      if (lso /= 0 .and. lwtkb /= -1) call mpibc2(orbtm,nl*nsp*nbas,4,3,mlog,'bndfp','orbtm')
      endif
      if (lpdiag == 2 .and. plbnd == 0) then
        eterms = s_ham%eterms
        call mpibc2(eterms(19),1,4,3,mlog,'bndfp','rhosig')
        s_ham%eterms = eterms
      endif
      if (loptic<0 .and. lwtkb /= -1 .and. njdosw>0) then
        call mpibc2(jdosw,size(jdosw),4,3,mlog,'bndfp','jdosw')
      endif
      if (mod(loptic,10)>0 .and. lwtkb /= -1) then
      call mpibc2(optmt,3*nfilm*nempm*nspx*nkp,4,3,mlog,'bndfp','optmt')
      endif

      entime = MPI_WTIME()
      call info2(20,0,0,' %;1d s',(entime-sttime),0)

C     Write bands in bands-plotting case: loop over qp getting evals from array
      if (plbnd == 1) then
      call info0(20,1,0,' Writing bands to bands file (MPIK) ...')
      if (procid == master) then
C     iq = running index to big qp list, i1 = index to current line
      iq = 0
      i2 = 0
  299 continue
      i = nsp
      if (onesp /= 0 .or. nspc == 2 .or. lafms) i = 1
      iopq = 0; if (pwmode > 10) iopq = 10+iopq
      call suqlstst(s_lat,plbopt,iopq,ndhamx,ef0,i,xv,nfbn,ifbls,i2,qp,onesp,spintexture)
      if (i2 <= 0) call rx0('bndfp')
      do  i1 = 1, i2
        iq = iq+1
        isqp = nsp*(iq-1)
        call suqlst(s_lat,plbopt,1,ndhamx,ef0,i,xv,nfbn,ifbls,i2,qp,onesp)
        call dpscop(qp,s_bz%qp,3,1,3*iq-2,1d0)
        do  isp = 1, nsp
          if (lafms .and. isp == 2) cycle
          ispc = min(isp,nspc)
C         jsp=isp in the collinear case; jsp=1 in the noncollinear
C         Thus jsp should be used in place of isp
C         isp serves as a flag for the noncollinear case
          if ((onesp==0 .or. isp==onesp) .and. (ispc==nspc)) then
            jsp = isp
            if (ispc == 2) jsp = 1
            call dpscop(evals,evlk(1,jsp),ndhamx,1+ndham*(jsp-1+isqp),1,1d0)
            if (mod(i1,10) /= 1) call pshpr(iprint()-6)
            call info5(30,0,0,' bndfp:  kpt %i of %i, k=%3:2,5;5d',i1,nkp,qp,0,0)
            if (mod(i1,10) /= 1) call poppr
            i = nsp; if (onesp /= 0 .or. nspc == 2 .or. lafms) i = 1
            k = jsp; if (onesp /= 0) k = 1
            call suqlsw(nevls(iq,jsp),k,i,evlk(1,jsp))
            if (nfbn(1) /= 0) then
              ndimhx = nevls(iq,jsp)  ! this should have been checked at generation time
              do  j = 1, 3
                k = j ; if (spintexture) k = 0
                if (nfbn(j) > 0)
     .            call suqlse(2,ndimhx,jsp,i,ndimhx,k,nfbn,ifbls,ndhamx,ndham*nsp,z,colwt(1,jsp,j,iq))
              enddo
            endif
          endif
        enddo
      enddo
      if (i /= 3) goto 299
      endif
      call rx0('done')
      endif

C     Repeat loop for printout.  Put evals back into local array
      do  iq = 1, nkp
        isqp = nsp*(iq-1)
        call dpscop(s_bz%qp,qp,3,3*iq-2,1,1d0)
C        call shorbz(qp,qp,qlat,plat)
        if (s_ctrl%lfp/2 /= 0 .and. .not. cmdopt('--shorbz=no',11,0,strn)) then
          xv(1:3) = qp
          call shorps(1,qlat,(/72,2,2/),xv,qp)
          if (ipr > 50 .and. dlength(3,xv(1:3)-qp,1) > 1d-8) write (stdo,578) iq,xv(1:3),qp
        endif
        do  isp = 1, nsp
        ispc = min(isp,nspc)
C       jsp=isp in the collinear case; jsp=1 in the noncollinear
C       Thus jsp should be used in place of isp
C       isp serves as a flag for the noncollinear case
        if ((onesp==0 .or. isp==onesp) .and. (ispc==nspc)) then
        jsp = isp
        if (ispc == 2) jsp = 1
        call dpscop(evals,evlk(1,jsp),ndhamx,1+ndham*(jsp-1+isqp),1,1d0)
        if (epsovl /= 0 .or. ovncut /= 0) then
          call dpscop(eomink,eomin,1,1+(jsp-1+isqp),1,1d0)
        endif
        if (mod(iq,10) /= 1) call pshpr(iprint()-6)

C       ndimhk must be MPI broadcast before using ...
C        call info5(30,0,0,' bndfp:  kpt %i of %i, k=%3:2,5;5d'//
C     .    '%?#n#   ndimh = %i##',
C     .    iq,nkp,qp,mod(pwmode/10,10),karchv(iq)%ndimhk)
        call info5(30,0,0,' bndfp:  kpt %i of %i, k=%3:2,5;5d',iq,nkp,qp,0,0)

        if (epsovl /= 0 .or. ovncut /= 0) then
          if (lpdiag /= 2) then
          call info5(30,0,0,' Overlap''s smallest eigenvalue: %;3g.  '//
     .      'H dim reduced to %i',eomin,nevls(iq,jsp),0,0,0)
          elseif (nevls(iq,jsp) /= ndimhx) then
            call info2(30,0,0,' ham dimension reduced to %i',nevls(iq,jsp),0)
          endif
        endif
        call prtev(xv,nevl,evlk(1,jsp),nmx,efmax,nev)
        if (mod(iq,10) /= 1) call poppr
        ebot = dmin1(ebot,evlk(1,jsp))
        i = max(1,nint(zval-qbg)/(3-nspc))
        evtop = max(evtop,evlk(i,jsp))
        ecbot = min(ecbot,evlk(i+1,jsp))
        if (lmetsav == NULLI) lmetsav = lmet
        if (lmetsav == 0 .and. iq == 1 .and. jsp == 1) ef0 = evtop
        if (plbnd == 0 .or. plbnd == 2) then
        if (ipr>=10 .and. iq==1 .and. ipl>1) write (stdl,712) (evlk(i,jsp),i=1,nev)
  712   format('fp evl',8f8.4)
        if (lwtkb /= -1 .and. .not. lwndow) then
          if (iq == 1 .and. jsp == nsp .and. .not. cmdopt('--no-fixef0',11,0,strn)) then
            ef00 = ef0
            call fixef0(zval-qbg,jsp,1,ndimh,ndham,evlk,dosw,ef0)
            if (jsp == 2 .and. ef00 /= ef0 .and.
     .        lwtkb == 0 .and. lmet > 0 .and. lrout /= 0) then
              if (procid == master) call info0(10,1,1,
     .          ' ... Fermi level reset in second spin channel ... restart band pass')
              goto 99
            endif
          endif
C         Check for cases when nevmx is too small : i=2 => fatal error
          i = 0
          if (nevmx>=0 .and. lmet /= 0) then
            dum = evlk(max(nev,1),jsp)
C           if (ef0 >= dum) i = 2
            if (ltet /= 1 .and. ef0+5*dabs(esmear-mpsord) > dum) i=2
            if (lmet==4 .and. ef0+def+5*dabs(esmear-mpsord)>dum)i=2
          endif
          if (i == 2) then
            call awrit3('%N evl(nev=%i)=%;3d but '//
     .        'ef0=%;3d ... restart with larger efmax or nevmx',
     .        ' ',80,stdo,nev,evlk(max(nev,1),jsp),ef0)
            call rx('bndfp')
          endif
        endif
        endif
C end second loop over isp
        endif
        enddo
C end second loop over iq
      enddo
      endif !second loop over qp (parallel k-points mode)

C ... Toggle s_optic%mefac and make 2nd pass if s_optic%mefac > 0
      if (lmefac > 0) then
        call dmcpy(evals(nfilo,1,1),ndhamx,1,evlda,nemup-nfilo+1,1,nemup-nfilo+1,nspx*nkp)
C        call prmx('original eval',evals(nfilo,1,1),ndhamx,nemup-nfilo+1,nspx*nkp)
C        call prmx('copied eval',evlda,nemup-nfilo+1,nemup-nfilo+1,nspx*nkp)
        if (allocated(karchv)) deallocate(karchv) ! Happens in MPIK case
C        if (lmet==5) then  ! should never happen
C          do  iq = 1, size(karchv)
C            if (associated(karchv(iq)%evals)) then
C              deallocate(karchv(iq)%evals,karchv(iq)%evecs)
C            endif
C          enddo
C          deallocate(karchv)
C        endif
        lmefac = -1             ! Toggle lmefac to do normal pass
        lwtkb=lwtkbsav; lmet=lmetsav; lrsig=lrsigsav; lrout=lroutsav
        loptic = s_optic%loptic
        call info0(30,1,1,' bndfp:  begin normal band pass')
        goto 40
      elseif (lmefac < 0) then
        call dmcpy(evals(nfilo,1,1),ndhamx,1,evnl,nemup-nfilo+1,1,nemup-nfilo+1,nspx*nkp)
c       if (lwtkb /= -1) lmefac = 1 ! Toggle lmefac back to starting point
      endif
      if (ldsdk .and. lwtkb >= 0) then
        call fpopm_sfac(2,s_lat,s_optic,ldsdk,dsdk,0,ndham,isp,nsp,nspc,nbas,nev)
      endif

C ... Print parameters for hilbert space reduction (epsovl or ovncut nonzero)
      if (mxevcut /= NULLI) then
        call mpibc2(mxevcut,1,2,1,mlog,'bndfp','mxevcut')
        call mpibc2(mnevl,1,2,2,mlog,'bndfp','mnevl')
        call info2(20,1,0,' BNDFP: max reduction in hilbert space = %i;'//
     .                    ' smallest dimension = %i',mxevcut,mnevl)
      endif

C ... Print upper and lower eval for each band
      if (llohi .and. procid == master) then
        do  isp = 1, nspx
          call info2(10,1,0,' Lowest, highest evals by band, spin %i',isp,0)
          call arrprt(' ','%,4i%,4;9D%,4;8D','Idd',ndimhx,0,3,
     .      0,'  | ',xv,elohi(1,isp,1),elohi(1,isp,2),xv,xv,xv,xv,xv)
        enddo
      endif

C ... Print orbital moment. SO by site, not by class
      if (lwtkb == 1 .and. lso /= 0) then
        call psymrq1(0,s_spec,s_site,nbas,nl,nsp,orbtm) ! Symmetrize orbital moment
        if (ltso == 1) then
          call mpibc2(tso,8*(1+nbas),4,3,0,' ',' ')
          call mpibc2(s_bz%sopertp(4),9,4,3,0,' ',' ')
          call psymrq1(1,s_spec,s_site,nbas,8,1,tso(1,1,1,1))
        endif
        if (procid == master) then
        allocate(ips(nbas))
        call sitepack(s_site,1,nbas,'spec',1,ips,xv)
        call iorbtm(s_spec,ltso,1,ips,nl,nl,nbas,nsp,orbtm,s_bz%sopertp,tso)
        deallocate(ips)
        endif
      endif

C ... Case generating bands: find next block of qp
      if (plbnd == 1) then
        if (allocated(nevls)) deallocate(nevls)
        goto 99
      endif
      if (plbnd == 3) goto 999
      if (lwden) goto 123     ! Integration weights already made
      if (ipl>1) write (stdl,715) nkp,ebot,zval,qbg,esmear
  715 format('nv nkp',i5,'  ebot',f9.4,'   qval',f10.4,'  qbg',f8.4,'  smr',f10.4)
C ... End of k point loop, all cases

C     call zprm3('smrho after k-point loop',0,srout,k1,k2,k3)

C --- Interpolate density to Fermi energy ---
      sev = sumev(1,1)
      if (lmet == 4) then
        call mshn3p(nbas,s_site,s_spec,lmet,lrout,lfrce,zval-qbg,ef0,
     .    def,sumqv,sumev,n1,n2,n3,k1,k2,k3,srout,frc,lrep)
C   ... Store val q & magnetic moment in sumqv(1) and sumqv(2)
        sumqv(2,1) = sumqv(1,1) - sumqv(1,2)
        sumqv(1,1) = sumqv(1,1) + sumqv(1,2)
C   ... Eigenvalue sum including entropy term
        sev = sumev(1,1) + sumev(2,1)
C   ... Remake sev,ef linearly interpolating tabulated sampling DOS
        sev00 = sev
        ef00  = ef0
        if (ldos /= 0) then
          call efldos(zval,nsp,emin,emax,ndos,dos,eferm,sev1)
          sev = sev1; ef0 = eferm
        endif
        if (ipr > 30 .and. ldos /= 0) then
          write(stdo,388) sev00,ef00,sev1,eferm,sev,ef0
  388     format(' ipol:  sev=',f12.6,'   ef=',f12.6:/
     .           ' dos:   sev=',f12.6,'   ef=',f12.6/
     .           ' use:   sev=',f12.6,'   ef=',f12.6)
           if (ipl>1) write (stdl,733) ef00,eferm,ef0,sev00,sev1,sev
  733      format('nf EF:',3f9.5,'    sev:',3f12.5)
        endif
        s_bz%ndos = ndos; s_bz%dosw = dosw; s_bz%ef = ef0; s_bz%def = def
        if (lrep == 1) then
          ef0 = -1
          call info0(2,0,0,'Input Fermi energy was too far off, repeat band pass')
          goto 99
        endif
      endif

C --- BZ integration for fermi level, band sum and qp weights ---
      if (lmet >= 0 .and. (lmet /= 4 .or. ltet /= 0)) then
        if (lwndow) then
          allocate(ww(ndham*nsp*nkp))
          eferm = min(dosw(1),dosw(2))
          call bzints(nkabc(1),nkabc(2),nkabc(3),evals,ww,nkp,
     .      ndham,ndham,nsp,xv,xv,xv,0,eferm,2,ntet,s_bz%idtet,sev,dum)
C         call prmx('w(min)',ww,ndham,ndham,nkp*nsp)
          eferm = max(dosw(1),dosw(2))
          call bzints(nkabc(1),nkabc(2),nkabc(3),evals,s_bz%wtkb,nkp,
     .      ndham,ndham,nsp,xv,xv,xv,0,eferm,2,ntet,s_bz%idtet,sev,dum)
C         call prmx('w(max)',s_bz%wtkb,ndham,ndham,nkp*nsp)
          call daxpy(ndham*nsp*nkp,-1d0,ww,1,s_bz%wtkb,1)
C         call prmx('w',s_bz%wtkb,ndham,ndham,nkp*nsp)
          deallocate(ww)
        else
          dosrng = 8
          if (mpsord < 0) dosrng = 16
          llmet = isw(lmet /= 0)
          if (lmet /= 0 .and. cmdopt('--ef=',5,0,strn)) then
            i = 5
            if (.not. a2bin(strn,eferm,4,0,' ',i,-1)) call
     .        rxs2('BNDFP: failed to parse "',strn(1:30),'%a"')
            llmet = 11
          endif
          if (.not. associated(s_bz%swtk)) s_bz%swtk => xx

          call bzwtsf(ndham,ndham,nsp,nspc,nkabc(1),nkabc(2),nkabc(3),
     .      nkp,ntet,s_bz%idtet,zval-qbg,fsmom,llmet,ltet,mpsord,
     .      ndos,dabs(esmear-mpsord),dosrng,s_bz%wtkp,evals,efmax,lswtk,
     .      s_bz%swtk,eferm,sev,xv,s_bz%wtkb,sumqv(1,2),s_bz%egap,lwtkb)
          if (fsmom(2) /= 0) s_ham%eterms(21) = fsmom(3)
C         Store val charge & magnetic moment in sumqv(1..2)
          if (lmet /= 4) then
            sumqv(1,1) = sumqv(1,2); sumqv(2,1) = sumqv(2,2)
          endif
        endif

        if (lmet /= 4) ef0 = eferm
        if (lmet /= 4) s_bz%ef = ef0
        if (lmet > 0) then
          if (procid == master) then
            ifi = fopna('wkp',-1,4)
            i = iobzwt(0,ndhamx,nkp,nspx,eferm,s_bz%wtkb,-ifi)
            call fclr('wkp',ifi)
          endif
        endif
C       skip second pass if only DOS are sought
        if (.not. lwndow .and. lquitdos) goto 123

        if (ltso /= 0 .and. s_bz%nef==3 .and. iabs(lwtkb)==1) then

          call info2(20,0,0,
     .      ' Sampling Fermi level and weights for LS pert',0,0)
          xv(1) = 0
          s_bz%sopertp(1) = eferm  ! In case Fermi level given in advance
          call bzwtsf(ndham,ndham,nsp,nspc,nkabc(1),nkabc(2),nkabc(3),
     .      nkp,ntet,s_bz%idtet,zval-qbg,xv,llmet,0,mpsord,ndos,
     .      dabs(esmear-mpsord),dosrng,s_bz%wtkp,evals,efmax,lswtk,
     .      s_bz%swtk,s_bz%sopertp(1),s_bz%sopertp(2),s_bz%sopertp(3),
     .      s_bz%wtkb(1+ndham*nkp*nsp),xv(6),xv(8),lwtkb)
C          call prmx('tetrahedron weights',s_bz%wtkb,ndhamx,12,nkp)
C          call prmx('sampling weights',s_bz%wtkb(1+ndham*nkp*nsp),
C     .      ndhamx,12,nkp)
C          stop

C          call bzwtsf(ndham,ndham,nsp,nspc,nkabc(1),nkabc(2),nkabc(3),
C     .      nkp,ntet,s_bz%idtet,zval-qbg,xv,11,ltet,mpsord,ndos,
C     .      dabs(esmear-mpsord),dosrng,s_bz%wtkp,evals,efmax,lswtk,
C     .      s_bz%swtk,eferm+def,sev,xv,s_bz%wtkb(1+ndham*nkp*nsp),
C     .      sumqv(1,2),xv(4),lwtkb)
C
C          call bzwtsf(ndham,ndham,nsp,nspc,nkabc(1),nkabc(2),nkabc(3),
C     .      nkp,ntet,s_bz%idtet,zval-qbg,xv,11,ltet,mpsord,ndos,
C     .      dabs(esmear-mpsord),dosrng,s_bz%wtkp,evals,efmax,lswtk,
C     .      s_bz%swtk,eferm-def,sev,xv,s_bz%wtkb(1+2*ndham*nkp*nsp),
C     .      sumqv(1,2),xv(4),lwtkb)

        endif

        if (lwtkb == -1 .and. lrout > 0) then
          call info0(20,0,0,' Start second band pass ...')
          lwtkb = mpipid(2) ! MPI barrier
          lwtkb = 1
          if (nspc == 2) lswtk = 1
          goto 99
        endif
        if (lwtkb == 2 .and. lrout > 0) then
          call info0(20,0,0,' New pass with constrained weights ...')
          goto 99
        endif
      endif

C ... Generate DOS on disk
  123 continue
      deallocate(dos)
      if (.not. lwndow .and. procid == master) then
      call iinit(dosisw,size(dosisw)); dosisw(1) = ndos; idfmt = 1; lidos = s_bz%ndos < 0
C     Parse --dos options if present
      if (cmdopt('--dos',5,0,strn)) then
        call sudossw(strn,1+2+4+16,doslsw,dosisw,dosw,' ')
        ndos = dosisw(1)
      endif
      if (IAND(s_ctrl%ldos,1) /= 0 .or. plbnd==0.and.cmdopt('--dos',5,0,strn)) then
        allocate(dos(3*ndos))
        call info2(30,1,0,' ... Generating %?#n>0# integrated#total# DOS',isw(lidos),0)
        ef0 = 0
        if (doslsw(8)) ef0 = s_bz%ef
        dosw(1) = dosw(1) + ef0; dosw(2) = dosw(2) + ef0
        if (ltet /= 0) then
          call bzints(nkabc(1),nkabc(2),nkabc(3),evals,dum,nkp,ndhamx,ndhamx,
     .      nspx,dosw(1),dosw(2),dos,ndos,eferm,1,ntet,s_bz%idtet,dum,dum)
          if (.not. lidos) call xxxdif(dosw(1),dosw(2),ndos,nspx,0,dos)
          del = 0d0
        else
          if (mpsord >= 100) mpsord = mod(mpsord,100)
          if (lidos) then
            call maknos(nkp,ndhamx,ndhamx,nspx,s_bz%wtkp,evals,mpsord,
     .      dabs(esmear-mpsord),-6d0,dosw(1),dosw(2),ndos,dos)
          else
            call makdos(nkp,ndhamx,ndhamx,nspx,s_bz%wtkp,evals,mpsord,
     .      dabs(esmear-mpsord),-6d0,dosw(1),dosw(2),ndos,dos)
          endif
          del = mpsord+dabs(esmear-mpsord)
        endif
        if (nspc == 2) call dscal(ndos,.5d0,dos,1)
        xv(1) = eferm
        i = 3; if (ldig) i = i+100
        call iodos(i,-fopn('DOS'),dos,ndos,nspx,ndos,1,dosw(1)-ef0,dosw(2)-ef0,nspx,eferm-ef0,idfmt)
        call fclose(fopn('DOS'))
        deallocate(dos)
      endif
      endif

C ... Dump weights to moms file, readable by lmdos
      if (lquitdos .and. procid==master) then
        do  iq = 1, nkp
          isqp = nsp*(iq-1)
          do  isp = 1, nsp
            if (isp == 2 .and. nspc > 1) cycle
            call dpscop(evals,evlk(1,isp),ndhamx,1+ndham*(isp-1+isqp),1,1d0)
            write (nfilem) 0, ndimhx
            call dpdump(evlk(1,isp),ndimhx,-nfilem)
            if (lafms) then
              write (nfilem) 0, ndimhx
              call dpdump(evlk(1,isp),ndimhx,-nfilem)
            endif
          enddo
        enddo
      endif

C ... Save Fermi level, nonmetal or sampling integration
      if ((lmet == 0 .or. .not. (lmet /= 4 .or. ltet == 1)) .and. .not.lwden) then
        if (lmet == 0) then
          ef0 = (evtop+ecbot)/2; eferm = ef0; s_bz%ef = ef0; s_bz%egap = ecbot-evtop
          call info5(20,0,0,' Highest occ. level = %,5;5d  Lowest unocc. = %,5;5d  '//
     .      'average = %,5;5d',evtop,ecbot,ef0,0,0)
C          if (ipr >= 20) call awrit3(' Highest occ. level = %,5;5d '//
C     .      ' Lowest unocc. = %,5;5d  average = %,5;5d',' ',80,
C     .      stdo,evtop,ecbot,ef0)
        endif
        if (procid == master) then
          ifi = fopna('wkp',-1,4)
          i = iobzwt(1,ndham,nkp,nsp,ef0,s_bz%wtkb,-ifi)
          call fclr('wkp',ifi)
        endif
      endif

C ... Collect partial dos from different disk files, write as a single file
C     if (nfilem > 0 .and. numprocs > 1) then
      if (nfilem > 0 .and. numprocs > 1 .and. (cmdopt('--mull',6,0,strn).or.cmdopt('--pdos',6,0,strn))) then
        if (allocated(kpproc)) deallocate(kpproc)
        allocate (kpproc(0:numprocs)); call dstrbp(nkp,numprocs,1,kpproc(0))
        allocate(doswtmpi(nchan,ndimhx,nsp,nkp)) ! evals is still allocated
        call dpzero(doswtmpi,size(doswtmpi))
        call iomomn(.true.,2,.false.,1,nspc,1,1,nfstg)
        rewind nfilem
        i = iomoms(nfilem,nl,nsp,nspc,nkp,
     .    ndham,nfstg,1,0,1,0,0,0,0,0,0d0,0d0,0d0,0d0,0d0,0d0)
        if (i < 0) call rx('failed to read moms file')
        do  iq = kpproc(procid), kpproc(procid+1)-1
          do  isp = 1, nsp
            nschan = mod(nfstg/10,10)
            i1 = iomoms(nfilem,nl,nsp,nspc,2,ndimh,nfstg,nschan,1,1,ndhamx,ndimhx,
     .        nchan,nchan,nev,evals(1,isp,iq),0d0,doswtmpi(1,1,isp,iq),0d0,0d0,0d0)
C            print 842, '2964',procid,iq,isp,doswtmpi(1,1,isp,iq)
C  842       format('bndfp ',a,3i4,10f12.6)
            if (nspc > 1) exit ! loop over isp is done through nfstg
          enddo                 ! loop over isp (main loop in parallel mode)
        enddo                   ! loop over iq (main loop in parallel mode)
        call dfclos(nfilem)     ! close and delete disk files
        call info2(20,0,0,' broadcasting partial DOS ...',entime-sttime,0)
        if (numprocs > 1) call mpibc2(doswtmpi,size(doswtmpi),4,3,.false.,'','')
C       print *, '2995', procid, sum(evals)

C       Master writes a single file with data for all qp
        if (procid == master) then
          nfilem = fopna('moms',-1,4)
C         Header
          i = iomoms(-nfilem,nl,nsp,nspc,nkp,
     .      ndham,nfstg,1,0,1,0,0,0,0,0,0d0,0d0,0d0,0d0,0d0,0d0)
          do  iq = 1, nkp
            do  isp = 1, nsp
              if (isp == 2 .and. nspc > 1) cycle
              i1 = iomoms(-nfilem,nl,nsp,nspc,2,ndimh,nfstg,nschan,1,1,ndhamx,ndimhx,
     .          nchan,nchan,nev,evals(1,isp,iq),0d0,doswtmpi(1,1,isp,iq),0d0,0d0,0d0)
            enddo
          enddo
        endif
      endif

C ... Cleanup ASA-style moments file, print table of DOS channels
      if (nfilem > 0 .and. procid == master) then
        i = iomoms(-nfilem,nl,nsp,nspc,nkp,ndimh,11,1,nkp*nsp+1,1,
     .    ndham,ndimh,nchan,nchan,ndimh,0d0,0d0,0d0,0d0,eferm,0d0)
        call fclose(nfilem)
        if (iprint() >= 10 .and. cmdopt('--mull',6,0,strn)) then
          call mchan(lmdim,s_site,s_spec,nsp,nsites,lsites,0,0,0,0,idoschan)
        endif
      endif

      if (allocated(nevls)) deallocate(nevls)
      if (allocated(idoschan)) deallocate(idoschan)

      if (cmdopt('--quit=dos',10,0,strn).and.plbnd==0 .or. llmet==11) then
        call rx0('completed band pass')
      endif

C ... Average forces so net force on system is zero (APW case)
      if (lfrce /= 0 .and. napw /= 0) then
        call dpzero(xv,3)
        do  i1 = 1, nbas
        do  i = 1, 3
          xv(i) = xv(i) + frc(i,i1,1)/nbas
        enddo
        enddo
        do  i1 = 1, nbas
        do  i = 1, 3
          frc(i,i1,1) = frc(i,i1,1) - xv(i)
        enddo
        enddo
      endif

C ... Approximate correction for NL potential
      if (lmefac /= 0) then
        xv(1) = eferm
        if (s_bz%egap /= NULLI) xv(1) = eferm + s_bz%egap/2
        call optsf(1+2*lwoptmc,nspx,nkp,xv(1),evnl,
     .    evlda,nfilo,nfiup,nemlo,nemup,optmt,optmc)
      endif

  499 continue  ! Entry point for --opt:read
      if (cmdopt('--opt',5,0,strn)) then
        dc = strn(6:6)
        ltmp = .false.
        if (procid == master) then
          if (wordsw(strn,dc,'read','',j) > 0) then
            ifi = fopna('optdata',-1,4)
            rewind ifi
            call info0(10,1,0,' reading optics data ...')
            read(ifi) eferm, i,j,k
            call sanrg(.true.,i,ndham,ndham,'bndfp:','file ham dimsn')
            call sanrg(.true.,j,nsp,nsp,'bndfp:','file nsp')
            call sanrg(.true.,k,nkp,nkp,'bndfp:','file nkp')
            call dpdump(evals,ndham*nsp*nkp,ifi)
            read(ifi) i,j
            if (loptic < 0) then
              call sanrg(.true.,i,njdosw,njdosw,'bndfp:','file njdosw')
              call sanrg(.true.,j,nfilm+nempm,nfilm+nempm,'bndfp:','file nfilm+nempm')
              call dpdump(jdosw,njdosw*(nfilm+nempm)*nspx*nkp,ifi)
            elseif (mod(loptic,10)>0) then
              call sanrg(.true.,i,nfilm,nfilm,'bndfp:','nfilm')
              call sanrg(.true.,j,nempm,nempm,'bndfp:','nempm')
              call dpdump(optmt,3*nfilm*nempm*nspx*nkp,ifi)
            endif

          elseif (wordsw(strn,dc,'write','',j) > 0 .or. wordsw(strn,dc,'woptmc','',j) > 0) then
            call info0(10,1,0,' writing optics data ...')
            if (mod(loptic,10)==0) call rx('no optical matrix elements calculated ... set OPTICS_MODE=1')
            if (wordsw(strn,dc//',','permqp','',j) > 0) then
              call info0(10,0,0,' reorder qp according to file qpts ...')
              ifi = fopna('QPTS',-1,1)
              call getqp(0,ifi,iq,iv,lshft,ntet,xx,xx,xx)
              call sanrg(.true.,iq,nkp,nkp,'bndfp:','qp list')
              allocate(ww(3*nkp),iwk1(nkp),wk(nkp))
              call getqp(2,ifi,iq,iv,lshft,0,ww,wk,xx)
              call alignqp(nkp,ww,nkp,s_bz%qp,plat,tolq,iwk1,i) ! iwk1 = permutation list
              call rxx(i /= 0,'qpts file is not compatible with available qp')
              call optpmme(3+4*lwoptmc,ndhamx,nspx,nkp,iwk1,evals,nfilo,nfiup,nemlo,nemup,optmt,optmc)
              call dcopy(3*nkp,ww,1,s_bz%qp,1)
              deallocate(ww,iwk1,wk)
            endif

            if (wordsw(strn,dc//',','write','',j)>0) then
              ifi = fopna('optdata',-1,4)
              rewind ifi
              write(ifi) eferm,ndham,nsp,nkp
              call dpdump(evals,ndham*nsp*nkp,-ifi)
              if (loptic < 0) then
                write(ifi) njdosw,nfilm+nempm
                call dpdump(jdosw,njdosw*(nfilm+nempm)*nspx*nkp,-ifi)
              endif
              if (mod(loptic,10)>0) then
                write(ifi) nfilm,nempm,nspx,nkp
                call dpdump(optmt,size(optmt),-ifi)
              endif
            endif
            if (wordsw(strn,dc//',','woptmc','',j)>0) then
              ifi = fopna('optdatac',-1,4)
              rewind ifi
              write(ifi) eferm,ndham,nsp,nkp
              call dpdump(evals,ndham*nsp*nkp,-ifi)
              write(ifi) nfilo,nfiup,nemlo,nemup,nspx,nkp,s_bz%qp
              call dpdump(optmt,size(optmt),-ifi)
              call dpdump(optmc,2*size(optmc),-ifi)
            endif
            ltmp = .true.
          endif
        endif
        call mpibc1(ltmp,1,1,.false.,'','')
        if (ltmp) call rx0('optics data written')
      endif
      if (nevmx > 0 .and. mod(iabs(s_optic%loptic),10) /= 0
     .  .and. lteto /= -1 .and. procid == master) then
        if (procid == master) then
        npol = 3
        if (loptic < 0) npol = 1
        if (loptic < 0 .and. njdosw > 0) npol = njdosw+1
        allocate(epsi(npopt*(1+nspx*npol)))
        call dpzero(epsi,npopt*(1+nspx*npol))
        if (lteto == 3) then
C         call prmx('jdosw',jdosw(1,1,1,1),njdosw,njdosw,(nfilm+nempm)*nkp)
          if (.not. associated(s_bz%ipq)) call rx('bndfp: ipq array never allocated')
          call optinq(s_ctrl,s_optic,s_bz,evals,ndham*nspc,nspx,nspc,
     .      eferm,s_bz%ipq,alat,plat,nfilm,nempm,npopt,npol,njdosw,
     .      optmt,jdosw,epsi)
        else
          call optint(s_optic,s_bz,lteto,evals,ndham*nspc,1,nspx,
     .      nspc,-1,eferm,s_bz%idtet,s_bz%wtkp,s_lat%vol,
     .      nfilm,nempm,npopt,npol,optmt,epsi)
        endif
        deallocate(epsi)
        else
        deallocate(optmt)
        allocate(optmt(1,1,1,1,1))
        endif
      endif
      if (cmdopt('--opt:read',10,0,strn)) call rx0('finished optics')

C --- Quantities involving output density ---
      if (lrout /= 0) then

C ... Core-level spectroscopy
      if (icls /= 0) then
        if (procid == master) then
          eferm = s_bz%ef
          call vcdmel(s_ctrl,s_site,s_spec,s_lat,nlmax,ndham,nphimx,ndimh,nkp,
     .      nsp,nspc,eferm,evals,ausc,nsites,isite,iclsl,iclsn)
        endif
        deallocate(ausc)
        call rx0('done generating core level spectra')
      endif

C ... Resonant X-ray emission spectroscopy
      if (irxes /= 0) then
        call info0(10,0,0,' doing RXES ...')
        if (procid == master) then
          call rxes(s_ctrl,s_site,s_spec,s_lat,nlmax,ndham,nphimx,ndimh,nkp,
     .      nsp,nspc,eferm,evals,ausc,nsites,isite,iclsl,iclsn,s_bz%wtkp)
        endif
        deallocate(ausc)
        call rx0('done generating rxes')
      endif

C ... Assemble output density
      call dfratm(s_site,s_spec,16+2,1,nbas,s_pot%rnew)
      qbyl => s_pot%qbyl
      allocate(hbyl(n0*nsp,nbas))
C     --window: Put output density into orhat and smrho, and exit
      if (lwndow .or. lwden) then
        call mkrout(' ',s_site,s_spec,s_lat,s_ham,nbas,nsp,ldim,lekkl,
     .    s_pot%rhat,hab,sab,qbyl,hbyl,lrout)
        call zcopy(k1*k2*k3*nsp,srout,1,s_pot%smrho,1)
        if (.not. lwden) then
        call symrho(s_site,s_spec,s_lat,7+10*lfrce,s_pot%smrho,s_pot%rhat,qbyl,hbyl,frc)
        endif
        goto 999
      endif
C     Dump [partial] atomic densities in atom files
      if (cmdopt('--wrhoat[:...]',8,0,strn)) then
        call mkrout(strn,s_site,s_spec,s_lat,s_ham,nbas,nsp,ldim,lekkl,
     .    s_pot%rnew,hab,sab,qbyl,hbyl,lrout)
      endif
      call mkrout(' ',s_site,s_spec,s_lat,s_ham,nbas,nsp,ldim,lekkl,
     .  s_pot%rnew,hab,sab,qbyl,hbyl,lrout)

C ... Symmetrize output density and forces
      call symrho(s_site,s_spec,s_lat,7+10*lfrce,s_pot%smrout,s_pot%rnew,qbyl,hbyl,frc)
C     call zprm3('smrho after symrho',0,srout,k1,k2,k3*nsp)
C     Average density by groups, render rho positive
      if (iand(s_ctrl%lcd,64) /= 0) then
        call grp2av(s_ctrl,s_site,s_spec,s_pot%rnew)
      endif
      if (iand(s_ctrl%lcd,256) /= 0) then
        i = 1
        call rhopos(i,s_site,s_spec,s_lat,s_pot%rnew,s_pot%smrout)
      endif

C ... New boundary conditions pnu for phi and phidot
C     call pshpr(iprint()-12)
      if (lpnu > 0) then
        xv(1:10) = s_ham%pmin; xv(11:20) = s_ham%pmax
        i = s_ham%pnudef
        if (i > 40) i = i-40  ! Strip off flag that retains lvs11 bug
        call pnunew(nbas,nsp,s_site,s_spec,xv,xv(11),10*i+lfrzw,hab,sab,qbyl,hbyl)
      endif
      endif                     ! output density branch

      if (cmdopt('--asars2',8,0,strn)) then
        i = fopna('rsta',-1,0)
        strn = 'ASA format, gen by lmf'
        i = iorsa(1,s_ctrl,s_site,s_spec,s_lat,s_pot,s_bz,strn,nbas,nbas,nglob('nspec'),
     .    xx,iter,.false.,-i)
        call rx0('ASA style rsta file written')
      endif

C --- Evaluate Harris energy ---
      call mkehkf(1,s_ham,sev,sumqv(2,1),sumtv,ehar)

C --- Evaluate KS total energy, correction to Harris force ---
      if (lrout /= 0) then

C   ... Correction to Harris force, nonzero when density not self-consistent
        call dfrce(s_site,s_spec,s_lat,s_ctrl,k1,k2,k3,nvl,s_pot%rhat,
     .    s_pot%rnew,elind,qmom,s_pot%smrho,srout,fh)

C   ... Evaluate KS total energy and output magnetic moment
        eks = 0d0
        if (leks >= 1) then
          call togpr()
          call mkekin(nbas,ldim,s_site,s_spec,s_lat,s_ham,
     .      s_pot%lcplxp,k1,k2,k3,s_pot%smvcnst,s_pot%smpot,srout,sev,sumtv)
          call pshpr(ipr-20)
          call mkpot(s_site,s_spec,s_lat,s_ham,s_pot,s_optic,nbas,lfrce,k1,k2,k3,srout,s_pot%rnew,
     .      qbg,s_pot%smpot,qmom,s_pot%smvcnst,ppnl,hab,vab,sab,
     .      qval,qsc,gpot0,vval,fes2,0,vorb,nlibu,lmaxu,lldau,nbf,s_ham%magf)
          call poppr()
          call mkehkf(2,s_ham,sev,sumqv(2,1),sumtv,eks)
          call togpr()
        endif

C   --- Add together force terms ---
        if (lfrce > 0)
     .    call totfrc(nbas,s_site,s_lat,leks,fes1,fes2,fh,frc,xv(1),s_ctrl%tol(5))

C   --- Mix input and output densities ---
        i1 = str_pack('mix',-2,s_strn,strn)
C#ifdefC MPE
C        ierr = MPE_LOG_EVENT(EVENT_START_MIXRHO,procid,"mixrho")
C#endif
        call mixrho(s_site,s_spec,s_lat,s_pot,srout,s_pot%smrho,
     .    nsp,iter,strn,zval-qbg,elind,k1,k2,k3,dmxp)
        if (dmxp(27)== 0 .or. dmxp(11)<dmxp(27)) then
          dmxp(26) = iter
          dmxp(27) = dmxp(11)
        endif
C#ifdefC MPE
C        ierr = MPE_LOG_EVENT(EVENT_END_MIXRHO,procid,"mixrho")
C#endif
      else
        eks = 0
      endif

      s_ham%ehf = ehar
      s_ham%ehk = eks
C#endif

C --- Cleanup ---
  999 continue
C#ifndef LMFGWD | LMFGW
      deallocate(s_optic%rgrad,s_optic%rgrade)
      if (size(optmt) /= 0) i = idalloc('optmt',allocvb()+3,0,0)
      i = idalloc('wtkb',allocvb()+3,0,0)
      i = idalloc('swtk',allocvb()+3,0,0)
      deallocate(optmt)
C#endif
      if (lrout /= 0) then
        deallocate(fh,fes2)
        if (.not. lwden) then
          deallocate(s_pot%qbyl,hbyl); nullify(s_pot%qbyl)
        endif
      else
        deallocate(srout)
      endif
      deallocate(qmom,gpot0,vval,fes1,vab)
      if (iand(s_ctrl%lfp,4) == 0) deallocate(s_pot%hab,s_pot%sab)
      if (.not. associated(s_pot%ppn)) then ! Deallocate ppnl unless it is to be preserved
        deallocate(ppnl)
      endif
      if (associated(evals) .and. .not.associated(evals,s_ham%evals)) then
        i = idalloc('evals',allocvb()+4,ndham,nsp*nkp)
        deallocate(evals)
      endif
      if (associated(colwt)) then
        i = idalloc('colwt',allocvb()+4,ndhamx*nspx,4*nkp)
        deallocate(colwt)
      endif
      deallocate(evlk)
      if (numprocs > 1) then
      if (epsovl /= 0 .or. ovncut /= 0) then
        deallocate(eomink)
      endif
      endif

      if (lmet==5) then
        do  iq = 1, size(karchv)
          if (associated(karchv(iq)%evals)) then
            deallocate(karchv(iq)%evals,karchv(iq)%evecs)
          endif
        enddo
        deallocate(karchv)
      endif
      call dfratm(s_site,s_spec,4+2,1,nbas,xv) ! Deallocate rho_out
      if (allocated(orbtm)) deallocate(orbtm)
      if (allocated(ausc)) deallocate(ausc)
      if (allocated(tso)) deallocate(tso)

C Patch to avoid PGI compiler bug on AMD processor
      call xxxbfp

      call tcx('bndfp')

      end

      subroutine xmpbnd(kpproc,ndham,nkp,nsp,eb)
C- Collect eb from various processors (MPI)
C ----------------------------------------------------------------------
Ci Inputs
Ci   kpproc
Ci   ndham :leading dimension of eb
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   eb    :energy bands; alias eband
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   05 Jul 13 Replace f77 pointers with f90 ones
Cu   10 Jul 06 First created
C ----------------------------------------------------------------------
      use mpi
      implicit none
C ... Passed parameters
      integer kpproc(0:*),ndham,nkp,nsp
      double precision eb(ndham,*)
C ... Dynamically allocated local arrays
      real(8), allocatable :: buf(:)
C ... Local parameters
      integer i,ista,iend
      integer procid
      integer numprocs, ierr
      integer, dimension(:),allocatable :: offset,length

      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

      allocate (offset(0:numprocs), stat=ierr)
      allocate (length(0:numprocs), stat=ierr)
      offset(0) = 0

      do  i = 0, numprocs-1
        ista = kpproc(i)
        iend = kpproc(i+1)-1
        length(i) = (iend - ista + 1)*nsp*ndham
        offset(i+1) = offset(i) + length(i)
      enddo
      ista = kpproc(procid)
      allocate(buf(ndham*nkp*nsp))
      call MPI_ALLGATHERV(eb(1,1+nsp*(ista-1)),
     ,  length(procid),mpi_real8,buf,length,
     ,  offset,mpi_real8,MPI_COMM_WORLD,ierr)
      call dcopy(ndham*nsp*nkp,buf,1,eb,1)
      deallocate(buf)

      deallocate(offset, stat=ierr)
      deallocate(length, stat=ierr)

C      if (procid == master) then
C        print *, procid,eb(1,1),eb(1,nsp*nkp)
C      endif
C      call rx('done')
      end
C Patch to avoid PGI compiler bug on AMD processor
      subroutine xxxbfp
      end
