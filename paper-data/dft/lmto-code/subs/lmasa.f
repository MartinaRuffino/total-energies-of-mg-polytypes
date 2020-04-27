      subroutine lmasa(prgnam,s_bz,s_ctrl,s_ham,s_lat,s_mix,
     .  s_pot,s_str,s_spec,s_site,s_tb,s_optic,s_gw,s_strn)
C- ASA self-consistency loop
C ----------------------------------------------------------------------
Cio Structures
C      ... Band code branch
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  nkabc nkp ef def w dos idtet ipq pdos qp star wtkp
Ci                 wtkb swtk zval lshft ntet lmet n range nevmx efmax
Ci                 lio ndos dosw fsmom
Co     Stored:     nevmx ef def w dos idtet ipq pdos qp star wtkp wtkb
Co                 swtk zval numq
Co     Allocated:  qp wtkb swtk
Cio    Elts passed:nkp wtkp qp ipq swtk wtkb idtet n w def star
Cio    Passed to:  popted rsedita iorsa bcast_strx iinit mpibc1 asars
Cio                clsprp getzv rdsigm bndasa subzi optint optin2
Cio                optinq intdrv gtpdss lmgdos asasx asalsq scrmom
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbas nclass nl nspec nspin loptc lasa lsx lncol lham
Ci                 lpgf sdmod sdprm maxit zbak lcgf npadl npadr nclasp
Ci                 ldlm tdlm nccomp ndlmcl lscr tol quit nvario nbasp
Ci                 lrs lmet ipcp nrc ips ipc rmax idcc ncomp ics initc
Ci                 lgen3 smalit lrel
Co     Stored:     lves lscr lasa lrs lham
Co     Allocated:  *
Cio    Elts passed: lscr lasa lrs lves ics nrc nvario ipc rmax ncomp
Cio                dclabl ipcp idcc ips lbxc initc lncol lrel lcd lham
Cio                lfp lgen3 ldos lstonr lqp
Cio    Passed to:  supot popted rsedita iorsa bcast_strx iinit mpibc1
Cio                asars clsprp subasi getzv asamad asavqm asaqmp
Cio                asvsph atscpp pp2alp suham bndasa nmpot secmat
Cio                secmtn secmt2 asaddq asaopm optint optin2 optinq
Cio                intdrv lmgdos dellmp asasx shoctl asalsq ioqpp
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  neula lsig eterms rdvext ehk thrpv seref ldham
Ci                 rsrnge pmin nmto kmto bdots etrms eula hrs iaxs
Ci                 iprmb lmxa magf nprs offH qsig pwmode pwemin pwemax
Ci                 npwpad lncol qss nbf sigp rsstol nlibu lmaxu udiag
Ci                 ndhrs elind
Co     Stored:     lsig nlibu lmaxu eterms ehk thrpv seref amgm bdots
Co                 etrms eula hrs iaxs iprmb lmxa magf nprs offH qsig
Co                 ndham ndofH ldham lmxax npwmin npwpad hord nqsig
Co                 ndhrs eseavr lham
Co     Allocated:  offH iprmb bdots nprs hrs qsig iaxs
Cio    Elts passed: etrms eterms iprmb entropy eula offH magf bdots
Cio                qsig hrs iaxs nprs
Cio    Passed to:  clsprp subasi asvsph bcast_strx iinit mpibc1 suham
Cio                rdsigm hft2rs bndasa nmpot secmat secmtn sstrxq
Cio                asaddq lmgdos chkdmu asasx asetot asalsq scrmom
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  avw nsgrp alat vol awald nkd nkq plat qlat nabc ag
Ci                 bgv cg cy dlv gv gvq indxcg ips0 istab jcg kv igv
Ci                 igv2 kv2 pos qlv symgr as tol nkdmx nkqmx platl
Ci                 platr ng tolft gmax
Co     Stored:     alat plat ag bgv cg cy dlv gv gvq indxcg ips0 istab
Co                 jcg kv igv igv2 kv2 pos qlv symgr ng gmax nabc
Co     Allocated:  gv kv igv kv2 igv2
Cio    Elts passed: symgr istab pos dlv qlv vol cg jcg indxcg ag gv
Cio                qlat igv igv2
Cio    Passed to:  supot popted rsedita iorsa bcast_strx iinit mpibc1
Cio                asars clsprp asamad asavqm asaqmp asvsph suham
Cio                sugvec sugvec0 rdsigm bndasa secmat secmtn sstrxq
Cio                asaddq asaopm lmgdos asasx asalsq scrmom
Cio  s_mix :contains parameters for charge mixing; see structures.h
Ci     Elts read:  b bv w wc nsave mmix umix tolu lxpot
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  nrhos pnu qnu nlml nlma ves aamom bxc cp ddpf dddpf
Ci                 ddpfr dlmwt dpf dpfr gibbs gma gmar grrme mad mxy
Ci                 palp papg pf pfnc pfr pmpol pp ppn pprel pti qc
Ci                 qcorr qpp qt rhat rnew rhos rhrmx sop thetcl vdif
Ci                 vintr vrmax vshft smpot smrho smrout vconst vmtz0
Ci                 vmtz
Co     Stored:     vconst ves aamom bxc cp ddpf dddpf ddpfr dlmwt dpf
Co                 dpfr gibbs gma gmar grrme mad mxy palp papg pf pfnc
Co                 pfr pmpol pnu pp ppn pprel pti qc qcorr qnu qpp qt
Co                 rhat rnew rhos rhrmx sop thetcl vdif vintr vrmax
Co                 vshft smpot smrho smrout vmtz nlma nlml
Co     Allocated:  mad pti ppn
Cio    Elts passed: pnu qnu vdif qc qt pp ves qpp rhos aamom vintr mad
Cio                vrmax rhrmx gibbs pmpol thetcl bxc pprel sop grrme
Cio                pti ppn
Cio    Passed to:  supot rsedita iorsa bcast_strx iinit mpibc1 asars
Cio                clsprp asamad asavqm asaqmp asvsph suham bndasa
Cio                nmpot secmat secmtn secmt2 asaddq asaopm lmgdos
Cio                asasx shoctl asalsq ioqpp
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  npr
Co     Stored:     nkaps nitab nttab lmaxw nds
Co     Allocated:  iax npr alp s adot sdot
Cio    Elts passed:alph kaps nkaps iax s sdot nitab adot
Cio    Passed to:  pp2alp suham rdstrx bndasa secmat secmtn secmt2
Cio                lmgdos dellmp asasx
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  z lmxa idmod grp2 name a nr rmt lmxl kmxt p pz lfoca
Ci                 qc rsma lmxb kmxv rsmv rfoca ctail etail stc nxi exi
Ci                 chfa rsmfa rhoc coreh coreq lmxf eref iq1 iq2 idxdn
Ci                 ncomp ngcut orbp idu uh jh hcr mxcst
Co     Stored:     name lmxa a nr rmt z lmxl kmxt p pz qc lfoca rfoca
Co                 coreh pb1 pb2 ctail etail stc nxi exi chfa rsmfa
Co                 rhoc idxdn ngcut orbp
Co     Allocated:  rhoc
Cio    Elts passed:rhoc
Cio    Passed to:  rsedita asars iorsa bcast_strx clsprp subasi getq
Cio                gtpcor asamad dlmq asavqm asvsph atscpp lmfitmr1
Cio                suham atfold makidx nscpa showbs sugcut uspecb
Cio                suldau sudmtu praldm rotycs symdmu sumlst bndasa
Cio                getidu nmpot secmat secmtn mullmf iorbtm lmgdos
Cio                lmrqit chkdmu asasx shoctl magtrq asalsq pvpqm1
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  class pnu pos spec force vel pz v0 v1 bxc cpawt omg
Ci                 omgn domg gc gcu gcorr sfvrtx j0 pdos rho1 rho2 rhoc
Ci                 rho1x rho2x rhocx qhhl qhkl qkkl eqhhl eqhkl eqkkl
Ci                 sighh sighk sigkk tauhh tauhk taukk pihh pihk pikk
Ci                 sighhx sighkx sigkkx tauhhx tauhkx taukkx pihhx
Ci                 pihkx pikkx thet clabel relax
Co     Stored:     pnu pos pos0 force vel spec clabel pz bxc cpawt omg
Co                 omgn domg gc gcu gcorr sfvrtx j0 pdos rho1 rho2 rhoc
Co                 rho1x rho2x rhocx qhhl qhkl qkkl eqhhl eqhkl eqkkl
Co                 sighh sighk sigkk tauhh tauhk taukk pihh pihk pikk
Co                 sighhx sighkx sigkkx tauhhx tauhkx taukkx pihhx
Co                 pihkx pikkx thet v0 v1 norb
Co     Allocated:  v0 v1
Cio    Elts passed:rhoc rho1 rho2
Cio    Passed to:  rsedita asars asars1 iorsa bcast_strx lmfitmr1 suham
Cio                setnorb showbs pvioeu suldau sudmtu praldm rotycs
Cio                symdmu sumlst bndasa getidu nmpot lmfitmr6 mullmf
Cio                lmgdos lmrqit chkdmu asasx magtrq
Cio  s_tb
Ci     Elts read:  fmode nbfit nbfitf wt alam alsc shft wg ebfit ndos
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:wg
Cio    Passed to:  *
Cio  s_optic:struct for parameters for optical ME; see structures.h
Ci     Elts read:  ltet ocrng unrng dw window nchi2 axes esciss lpart iq
Co     Stored:     ocrng unrng
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  bndasa asaopm optint optin2 optinq intdrv lmgdos
Cio  s_strn :struct for global strings; see structures.h
Ci     Elts read:  strn
Co     Stored:     strn
Co     Allocated:  *
Cio    Elts passed:strn
Cio    Passed to:  suham str_pack asasx
C      ... Crystal GF branch
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  semsh nkabc nkp ef def w dos idtet ipq pdos qp star
Ci                 wtkp wtkb swtk zval lshft dosw ndos lmet
Co     Stored:     semsh nevmx ef def w dos idtet ipq pdos qp star wtkp
Co                 wtkb swtk zval
Co     Allocated:  *
Cio    Elts passed:nkp wtkp qp ipq star
Cio    Passed to:  popted rsedita iorsa bcast_strx iinit mpibc1 asars
Cio                clsprp setne getzv rdsigm gfasa gfgw gfsigma asalsq
Cio                scrmom
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbas nclass nl nspec nspin loptc lasa lsx lncol lham
Ci                 lpgf sdmod sdprm maxit zbak lcgf npadl npadr nclasp
Ci                 ldlm tdlm nccomp ndlmcl lbxc npl lscr tol quit
Ci                 nvario nbasp lrs lmet ipcp nrc ips ipc rmax idcc
Ci                 ncomp ics initc lgen3 smalit lordn lrel ncl clssl
Ci                 clp ldomg
Co     Stored:     lves lscr lrs ldomg
Co     Allocated:  *
Cio    Elts passed: lscr lasa lrs lves ics ncomp idcc nrc nvario pgplp
Cio                ipc dclabl rmax ipcp ips lbxc initc lncol lrel lcd
Cio                lham lfp lgen3 ldos clssl clp
Cio    Passed to:  supot popted rsedita iorsa bcast_strx iinit mpibc1
Cio                asars clsprp subasi setne dlminit getzv asamad
Cio                asavqm asaqmp asvsph atscpp pp2alp suham gfasa
Cio                suordn suiord sugemb suhamz mkpotf makpfz mkpotf
Cio                mkcpa sugfrs gfordn gfibz gfomg gfp0ft hamfbz
Cio                gfenint mixomg mkgint shoctl gfgw gfsigma asalsq
Cio                ioqpp
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  neula lsig eterms rdvext ehk thrpv seref ldham
Ci                 rsrnge pmin nmto kmto bdots etrms eula hrs iaxs
Ci                 iprmb lmxa magf nprs offH qsig pwmode pwemin pwemax
Ci                 npwpad lncol qss nbf sigp rsstol nlibu lmaxu udiag
Ci                 hord lgen3 bandw ndhrs elind
Co     Stored:     lsig nlibu lmaxu eterms ehk thrpv seref amgm bdots
Co                 etrms eula hrs iaxs iprmb lmxa magf nprs offH qsig
Co                 ndham ndofH ldham lmxax npwmin npwpad hord nqsig
Co                 ndhrs eseavr
Co     Allocated:  etrms offH iprmb bdots nprs hrs qsig iaxs
Cio    Elts passed: etrms eterms iprmb entropy eula offH magf bdots
Cio                qsig hrs iaxs nprs
Cio    Passed to:  clsprp subasi dlminit asvsph bcast_strx iinit mpibc1
Cio                suham rdsigm hft2rs gfasa suhamz mkpotf makpfz
Cio                mkpotf mkcpa gfibz gf1kp gfg2g gfomg gfp0ft gfp0f1
Cio                pdglrs gfdpp gfg0g gfenint mkgint gfgw gfsigma
Cio                chkdmu asetot asalsq scrmom
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  avw nsgrp alat vol awald nkd nkq plat qlat nabc ag
Ci                 bgv cg cy dlv gv gvq indxcg ips0 istab jcg kv igv
Ci                 igv2 kv2 pos qlv symgr as tol nkdmx nkqmx platl
Ci                 platr ng tolft gmax
Co     Stored:     alat plat ag bgv cg cy dlv gv gvq indxcg ips0 istab
Co                 jcg kv igv igv2 kv2 pos qlv symgr ng gmax nabc
Co     Allocated:  gv kv igv kv2 igv2
Cio    Elts passed: symgr istab pos dlv qlv vol cg jcg indxcg ag gv
Cio                qlat igv igv2 cy
Cio    Passed to:  supot popted rsedita iorsa bcast_strx iinit mpibc1
Cio                asars clsprp asamad asavqm asaqmp asvsph suham
Cio                sugvec sugvec0 rdsigm gfasa suordn mkcpa sugfrs
Cio                gfordn gfp0ft gfgw gfsigma asalsq scrmom
Cio  s_mix :contains parameters for charge mixing; see structures.h
Ci     Elts read:  b bv w wc nsave mmix umix tolu lxpot
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  nrhos dlmwt pnu qnu nlml nlma ves aamom bxc cp ddpf
Ci                 dddpf ddpfr dpf dpfr gibbs gma gmar grrme mad mxy
Ci                 palp papg pf pfnc pfr pmpol pp ppn pprel pti qc
Ci                 qcorr qpp qt rhat rnew rhos rhrmx sop thetcl vdif
Ci                 vintr vrmax vshft smpot smrho smrout vconst vmtz0
Ci                 vmtz
Co     Stored:     gibbs vconst ves aamom bxc cp ddpf dddpf ddpfr dlmwt
Co                 dpf dpfr gma gmar grrme mad mxy palp papg pf pfnc
Co                 pfr pmpol pnu pp ppn pprel pti qc qcorr qnu qpp qt
Co                 rhat rnew rhos rhrmx sop thetcl vdif vintr vrmax
Co                 vshft smpot smrho smrout vmtz nlma nlml
Co     Allocated:  mad pti qcorr cp pf dpf ddpf dddpf pfr dpfr ddpfr
Co                 gmar papg palp gma
Cio    Elts passed: pnu qnu vdif bxc dlmwt vshft qc qt gibbs pp ves qpp
Cio                rhos aamom mxy thetcl mad vrmax rhrmx pmpol pprel
Cio                sop grrme vintr pti qcorr palp pf dpf ddpf dddpf gma
Cio                pfr ddpfr dpfr papg gmar cp
Cio    Passed to:  supot rsedita iorsa bcast_strx iinit mpibc1 asars
Cio                clsprp dlminit asamad asavqm asaqmp asvsph suham
Cio                gfasa suhamz mkpotf makpfz mkpotf mkcpa gfibz gf1kp
Cio                gfg2g gfomg gfp0ft gfp0f1 pdglrs gfdpp gfg0g gfenint
Cio                mkgint shoctl gfgw gfsigma asalsq ioqpp
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  npr
Co     Stored:     nkaps nitab nttab lmaxw nds
Co     Allocated:  iax npr alp s adot sdot
Cio    Elts passed:alph kaps iax npr s
Cio    Passed to:  pp2alp suham rdstrx gfasa
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  z lmxa idmod grp2 name a nr rmt lmxl kmxt p pz lfoca
Ci                 qc rsma lmxb kmxv rsmv rfoca ctail etail stc nxi exi
Ci                 chfa rsmfa rhoc ncomp ncpa nthet iscpa beff coreh
Ci                 coreq lmxf eref iq1 iq2 idxdn ngcut orbp idu uh jh
Ci                 hcr mxcst
Co     Stored:     name lmxa a nr rmt z lmxl kmxt p pz qc lfoca rfoca
Co                 coreh pb1 pb2 ctail etail stc nxi exi chfa rsmfa
Co                 rhoc idxdn ngcut orbp
Co     Allocated:  rhoc
Cio    Elts passed:rhoc xcpa beff
Cio    Passed to:  rsedita asars iorsa bcast_strx clsprp subasi dlminit
Cio                dlmwgts getq gtpcor asamad dlmq asavqm asvsph atscpp
Cio                lmfitmr1 suham atfold makidx nscpa showbs sugcut
Cio                uspecb suldau sudmtu praldm rotycs symdmu sumlst
Cio                gfasa getidu suordn suhamz mkpotf makpfz mkpotf
Cio                vorbydl shoctl iorbtm chkdmu magtrq asalsq pvpqm1
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  ncomp class pnu pos spec force vel pz v0 v1 bxc
Ci                 cpawt omg omgn domg gc gcu gcorr sfvrtx j0 pdos rho1
Ci                 rho2 rhoc rho1x rho2x rhocx qhhl qhkl qkkl eqhhl
Ci                 eqhkl eqkkl sighh sighk sigkk tauhh tauhk taukk pihh
Ci                 pihk pikk sighhx sighkx sigkkx tauhhx tauhkx taukkx
Ci                 pihhx pihkx pikkx thet dlmcl clabel norb relax
Co     Stored:     pnu pos pos0 force vel spec clabel pz bxc cpawt omg
Co                 omgn domg gc gcu gcorr sfvrtx j0 pdos rho1 rho2 rhoc
Co                 rho1x rho2x rhocx qhhl qhkl qkkl eqhhl eqhkl eqkkl
Co                 sighh sighk sigkk tauhh tauhk taukk pihh pihk pikk
Co                 sighhx sighkx sigkkx tauhhx tauhkx taukkx pihhx
Co                 pihkx pikkx thet v0 v1 ncomp norb
Co     Allocated:  v0 v1 cpawt thet j0 bxc omg omgn domg gc gcorr gcu
Co                 sfvrtx pdos
Cio    Elts passed: j0 ncomp rhoc rho1 rho2 thet cpawt bxc pdos dlmcl
Cio                gc gcorr omgn gcu domg sfvrtx omg
Cio    Passed to:  rsedita asars asars1 iorsa bcast_strx dlminit pkgdlm
Cio                lmfitmr1 suham setnorb showbs pvioeu suldau sudmtu
Cio                praldm rotycs symdmu sumlst rdomg gfasa getidu
Cio                suordn sugemb suhamz mkpotf makpfz mkpotf mkcpa
Cio                vorbydl gfomg mixomg womg mkgint cpomg chkdmu
Cio                dlmentrp magtrq
Cio  s_tb
Ci     Elts read:  fmode nbfit nbfitf wt alam alsc shft wg ebfit ndos
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_optic (not used)
C      ... Layer GF branch
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbas nbasp nl nspin lncol lpgf lham nclasp nofgl
Ci                 nofgr nclass nspec lgen3 npl lrel lasa ldlm nccomp
Co     Stored:     lham
Co     Allocated:  *
Cio    Elts passed:ldos ics ips ipc nrc dclabl initc
Cio    Passed to:  setne supgh supgbc supghz suhamz mkpotf makpfz
Cio                mkpotf pgsif pgflu shoctl
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pos
Cio    Passed to:  supgbc pgsif pgflu
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  z idxdn lmxb lmxa ncomp idu a nr rmt hcr
Co     Stored:     idxdn lmxa
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  pgflu7 makidx nscpa supgh getidu supgbc supghz
Cio                suhamz mkpotf makpfz mkpotf pgsif shoctl
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  ncomp dlmcl spec class pnu clabel v0
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  getidu supgbc supghz suhamz mkpotf makpfz mkpotf
Cio                pgsif
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  lmet nkp semsh dosw ndos
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:qp wtkp
Cio    Passed to:  setne
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  rmax npr
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:npr iax alph s nds
Cio    Passed to:  supgh supgbc pgsif plham pgbevl pgkmap pgflu sblham
Cio                pgemb pvemb pgdysn pgles pgcurr
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  nlibu lmaxu udiag eterms iprmb offH ldham lncol hord
Ci                 lgen3 lham neula bandw
Co     Stored:     offH iprmb eterms bandw
Co     Allocated:  offH iprmb
Cio    Elts passed:iprmb offH lncol eula
Cio    Passed to:  supgh pglusu supgbc supghz suhamz mkpotf makpfz
Cio                mkpotf pgsif plham pgbevl pgkmap pgflu sblham pgemb
Cio                pvemb pgdysn pgles pgcurr gfg2g
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  vmtz ves palp
Co     Stored:     *
Co     Allocated:  cp pf dpf ddpf dddpf pfr dpfr ddpfr gmar papg palp
Co                 gma pfnc
Cio    Elts passed: pp pf qc gma palp dpf ddpf dddpf pprel pfr ddpfr
Cio                dpfr papg gmar sop pfnc pnu qnu ves
Cio    Passed to:  supgh supgbc supghz suhamz mkpotf makpfz mkpotf
Cio                pgsif plham pgbevl pgkmap pgflu sblham pgemb pvemb
Cio                pgdysn pgles pgcurr gfg2g shoctl
Co Outputs
Cs Command-line switches
Cs   --band      : Tabulate energy bands; see doc/Command-line-options.html
Cs   --bonly     : Use in conjuction with --lxb switch
Cs   --chklim    : Check band limits in optics calculations
Cs   --cv:       : Calculate electronic specific heat, eV
Cs   --cvK:      : Calculate electronic specific heat, Ry
Cs   --dumprho   : Not documented
Cs   --invbl     : Not documented
Cs   --jdosw     : Channels for optical properties; See doc/optics.html
Cs   --jdosw2    : Channels for optical properties; See doc/optics.html
Cs   --lxb       : Check band limits in optics calculations
Cs   --mix=      : Not documented
Cs   --mixsig=   : For self-energy; see Command-line-options.html
Cs   --mlog      : (MPI) write MPI commands to log file
Cs   --mull      : Mulliken analysis; see doc/Command-line-options.html
Cs   --nosymdm   : Not documented
Cs   --oldbz     : Not documented
Cs   --onesp     : Generate bands for one spin; use with --band
Cs   --opt:read  : Read optical matrix elements; See doc/optics.html
Cs   --opt:write : Write optical matrix elements; See doc/optics.html
Cs   --pdos      : Partial DOS; see doc/Command-line-options.html
Cs   --popted    : Invoke the optics editor
Cs   --rs=       : Controls I/O with rst file; see Command-line-options.html
Cs   --rsedit    : Invoke the restart file editor
Cs   --rsig      : For reading the self-energy; see Command-line-options.html
Cs   --sh=       : Not documented
Cs   --wsig      : For writing the self-energy; see Command-line-options.html
Cs   -ef=        : Overwrite Fermi level; use with --band
Cs   -enu=       : Not documented
Cs   -noves      : Use Madelung potential from disk; do not remake from density
Cs   -novintr    : Exclude intra-atomic W when making screened exchange sigma
Cr Remarks
Cr  This is the entry point for the self-consistent ASA band program
Cr
Cr  lmasa contains two main blocks.  One part generates the potential
Cr  from input moments and generates potential parameters as output.
Cr  The ASA interatomic electrostatic potential is generated by asamad,
Cr  while asvsph handles the sphere part.
Cr
Cr  The other block generates eigenvalues, Green's functions, and
Cr  the like from input potential parameters.  One output of
Cr  these blocks are moments, which complete the cycle.  The
Cr  routines generating outputs of the hamiltonian are:
Cr     bndasa:  Eigenvectors and output moms from the ASA hamiltonian.
Cr     gfasa:   crystal Green's function
Cr     pgfasa:  layer Green's function
Cr  Also:
Cr     asasx:   Rucker's screened exchange potential is implemented to
Cr              date and requires output of bndasa.
Cr     magtrq:  generates magnetic forces, rotates magnetic spins
Cr              based on output from noncollinear bndasa, gfasa.
Cr
Cr  The first band iteration is nit=1.  It is possible to do
Cr  a band calculation without moments as input, provided pot. pars.
Cr  are supplied.
Cr
Cl Local variables
Cl   ehterm: parameters needed to make Harris energy
Cl            (1) etot for VH(rmt)=0 (2) sumeV(VH=0)
Cl            (3) sum Q_R V_R        (4) emad
Cl            (5) int input applied field * input density
Cl            (6) int input applied field * output density
Cl          NOTE: this array is archaic.  Will be superseded
Cl                by s_ham%eterms when the latter is stabilized.
Cl   latok :1s digit
Cl           0 input potential not read or created
Cl           1 input potential is available
Cl         10s digit
Cl          0 parameters to gen. total energy available
Cl          1 parameters to gen. total energy not available
Cl   makepp: F do not update potential or make ppars
Cl           T update potential and make ppars from moments
Cl   etrmss: backup double-counting terms s_ham%eterms, used
Cl         : when potential, P,Q are kept frozen (see makepp=F)
Cl         : in e.g. frozen-potential SD
Cl   nclspd: number of CPA classes
Cl   nwmoms: T output moments have been generated
Cl   sdmod : see doc/nc.txt.
Cl           1s digit:
Cl            0 Output Euler angles from density matrix
Cl            1 relax along force
Cl            2 dynamics
Cl           10s digit (for spin statics only)
Cl            0 mix Euler angles with charge moments P,Q
Cl            1 independently mix Euler angles
Cl           1000s digit
Cl            causes lm to prevent updating of atomic P and Q,
Cl            and the potential parameters.
Cl  vconst  :(gf,pgf) Constant estat potential shifts, used where
Cl          : the Fermi level is specified and the potential
Cl          : adjusts to it.
Cl          : vconst(1) = potential shift
Cl          : vconst(2) = potential shift of L end region (PGF)
Cl          : vconst(3) = potential shift of R end region (PGF)
Cl   vrl    : Difference in estat potential in the left- and right- electrodes
Cl          : Consider the equilibrium case first.
Cl          : Right now, the code solves the electrostatics by Ewald
Cl          : summation using periodic boundary conditions.  To make the
Cl          : potential continuous everywhere, an additional linear
Cl          : solution V(z)=[vconst(3)-vconst(2)]*z/Lz is required.  This
Cl          : drop in potential will be distributed properly by the
Cl          : self-consistency condition, satisfying boundary conditons
Cl          : that the potential difference across the leads is
Cl          : [vconst(3)-vconst(2)].  Thus, the potential is continuous.
Cl
Cl          : In the nonequilibrium case, the bias across the device
Cl          : creates a difference in Fermi level across the left and right
Cl          : leads.  Combining the two effects,
Cl          : vrl = [vconst(3)-vconst(2)] + ef(R)-efermi, where
Cl          : ef(R) and efermi are fermi energies of right and left leads.
Cl   nangl  : total number of DLM angles summed over all sites
Cb Bugs
C    For DLM pvpqm1 needs changes if something must be mixed
Cu Updates
Cu   28 Oct 17 New functionality with --rs= : irs(5) passed to asvsph as 10s digit imake
Cu    5 Oct 12 (Belashchenko) emesh mode 5 (imaginary axis)
Cu   17 Sep 17 Pnu not floated when s_ctrl%lbas,16 is set
Cu   21 Apr 15 Reworking of electrostatics for fully relativistic case
Cu   21 Sep 13 New shfac, vshft for the band program
Cu   08 May 13 Complete migration to f90 structures; eliminate s_array
Cu   25 Apr 12 (Belashchenko) CPA extended to treat chemical disorder
Cu   26 Dec 11 (Belashchenko) Major cleanup of the DLM part
Cu   06 Sep 11 Started migration to f90 structures
Cu   12 Apr 11 Improved treatment of orbital-dependent noncollinearity
Cu   29 Nov 10 (P. Larson) Added corrections for DLM
Cu   17 Nov 10 L-M fitting with C,enu in tandem
Cu   15 Aug 10 L-M fitting with MPI
Cu   01 Jul 10 Implement L-M fitting for spin-pol case
Cu   12 May 10 New features in Levenberg-Marquardt fit
Cu   16 Apr 10 Levenberg-Marquardt fit of ASA hamiltonian to DOS
Cu   11 Aug 09 Start of Levenberg-Marquardt fit of ASA hamiltonian to bands
Cu   08 Nov 07 (J. Xu) New LDA+U ASA hamiltonian
Cu   19 Jul 07 (S.Faleev) corrected electrostatics in pgf case, when
Cu             left- and right- leads are different
Cu   17 Mar 05 (T.Sandu) GF can generate SX sigma
Cu   10 Jul 04 (S.Faleev) Changes to handle non-equilibrium mode
Cu   21 Apr 04 Additions for an m-dependent spin-density matrix
Cu   23 Sep 03 SX patterned after GW.  sigm(irr k) generated
Cu             and stored in file sigm.  sigm(rs) obtained
Cu             by call to rdsigh.
Cu   04 Sep 03 Expand dxmprm so wa may have independent value
Cu   30 Apr 03 Added MPI read/write for rsta.
Cu   20 Mar 03 Added etol as a tolerance for self-consistency
Cu   10 Mar 03 First cut at making proper Kohn-Sham energy.
Cu   17 Feb 03 Added double-counting terms for applied field
Cu   27 Jan 03 Changes to noncollinear code
Cu   24 May 02 small changes to better decouple sx and scr modes
Cu   01 Mar 02 revisions to accomodate changes to PGF codes
Cu   19 Feb 02 revisions to accomodate changes to GF codes
Cu   07 Jun 01 use call to nwit to check convergence
Cu   17 May 01 First attempt for I/O from rst file
Cu   17 Sep 99 start streamlining with sham,spot,fewer passed arguments
Cu   28 Apr 98 add code for radial matrix elements of grad
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters:
      character*(*)  prgnam*8
C ... For structures
!      include 'structures.h'
      type(str_bz)::    s_bz
      type(str_ctrl)::  s_ctrl
      type(str_ham)::   s_ham
      type(str_lat)::   s_lat
      type(str_mix)::   s_mix
      type(str_pot)::   s_pot
      type(str_str)::   s_str
      type(str_tb)::    s_tb
      type(str_gw)::    s_gw
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
      type(str_optic):: s_optic
      type(str_strn) :: s_strn(*)
C ... Dynamically allocated local arrays
      integer, allocatable :: lmx(:),idmod(:,:),grp2(:),llmx(:)
      real(8), allocatable :: eband(:,:),frc(:),dq(:),z(:)
      real(8), allocatable :: pold(:,:),qold(:,:),vold(:),eulo(:),qrold(:)
      real(8), allocatable :: xnew(:),xold(:)
      real(8), allocatable :: wk(:),wk2(:),wkr(:),temp(:),shfac(:,:)
      complex(8), allocatable :: qppo(:)
      real(8), pointer     :: pnus(:,:),qnus(:,:)
C     real(8), pointer :: vshfs(:)

C     For LDA+U
      integer, allocatable :: lldau(:)
      complex(8), allocatable :: vorb(:),dmatu(:),dmato(:)

C ... Local variables
      character outs*120, sxopts*120, strn*160
C     character ch*1
      logical mlog,T,F,makepp,nwmoms,lmfrce,aintra,sw,lfree,lso
      parameter (T=.true., F=.false.)
      integer i,i1,i2,ifi,imake,stdo,lham,lves,iprt,iv(10),
     .  lasa,lpgf,lcgf,lncol,lrel,lsx,lnsph,loptc,nband,nbas,nclasp,
     .  nclspp,nclass,nevmx,nl,nlspcp,nsp,nspx,nspc,nspec,nlspcr,
     .  sdmod,iscr,neul,nrhos,latok,iscsav,str_pack
      integer,parameter:: NULLI = -99999
      integer mpipid,procid,master
      double precision amgm,ehf,ehk,ehk0,emad,etol,seref,sevat,sumev,
     .  thrpv,trumad,xx,elind
      real(8), target  :: xxp(1)
      double precision ehterm(6),vmtz(2),zval,zbak(2),amag(3),efermi(2),
     .  avw,xv(10),zdel(2),pmin(10),etot(2),etrmss(22),vconst(3),vcnsts(3)
      equivalence (etot(1),ehf),(etot(2),ehk)
C ... Parameters for iterations
      integer lsc,irs(5),nit,maxit
C ... Parameters for nonspherical density
      integer nqpp
C ... Parameters for Levenberg-Marquardt fitting
      integer,allocatable:: ivcof(:,:),lsite(:)
      real(8),allocatable:: ebf(:,:),sigmrq(:,:)
      real(8),allocatable:: gband(:,:,:),gdos(:,:,:)
      real(8),allocatable:: ewk(:,:),gwk(:,:,:),ppsave(:,:,:,:)
      real(8),allocatable:: ftcof(:),ftalp(:,:),ftcov(:,:),ftwk(:,:)
      real(8),allocatable:: refdos(:,:),fitdos(:,:),refdos2(:,:)
      integer fmode,fmode0,fmode1,nfit,nfcof,nfvar,nbf,nqf,nbfit(2),
     .  nbfitf,nchan,nsgrp,moddos,nsite,nfilem
      double precision fitmrw(3),ftchi(4),ftalam,ftalsc
C     1 = lower bound of ref dos
C     2 = upper bound of ref dos
C     3 = number of reference DOS points (integer)
C     4 = lower bound of fit dos
C     5 = upper bound of fit dos
C     6 = number of fit DOS points (integer)
C     7 = lower bound of ref dos, shifted by ref dos EF
C     8 = upper bound of ref dos, shifted by ref dos EF
C     9 = energy shift added ref dos to align Ef with fit dos
C    10 = lower bound of fit dos, shifted by fit dos EF
C    11 = upper bound of fit dos, shifted by fit dos EF
C    12 = number of confined fit DOS points (integer)
      double precision dosprm(12)
C ... Parameters for spin dynamics
      double precision sdprm(6)
C ... Parameters for mixing.  Default parameters dmxprm:
C 1(I) mixing scheme; 2 beta; 3 wc; 4,5 wgt; 6(I)nsave 7(I) mmix;
C 8(I) nkill; 9 betv; 10 rmscst.  11-20 are outputs:
C 11 rmsdel 12 rms2 13 nmix 14 actual broy 15 actual beta 16-24 tj
C 25 1 if wt ne 0, 10 if wa ne 0, 11 if all nonzero
C 27..29: hold parms for static parms block regular mixing
C 30..32: hold parms for static parms block Euler angle mixing
C 33 : Lindhard screening parameter
C 34 : wa
      integer, allocatable :: mixcst(:)
      double precision dmxprm(34),dmxeu(34),rmsdel,rms2,betsv,rmscst,betv,qtol
      equivalence (dmxprm(9),betv),(dmxprm(10),rmscst),
     .  (dmxprm(11),rmsdel),(dmxprm(12),rms2),(dmxprm(15),betsv)

C ... Parameters for GF
C#ifdefC GF
C      logical a2bin
C      integer nzp,npl
C      integer oomega2
CC ... Parameters for non-equilibrium mode
C      integer nzne,j
C      logical lnoneq
C
C      complex(8), pointer :: zp(:),wz(:)
C#endif
      double precision vne,vrl
      integer npadl,npadr,nbasp
C ... Parameters for DLM
      integer ldlm,ldlm1,nangl,ndlmcl,nclspd,nclsppd
      double precision omdiff,tdlm !,wk(2)
      logical omonly,lweiss
C#ifdefC GF
C      integer maxth
C      parameter (maxth=100)
C      integer ib,ifac(3),lshft(3)
C      double precision rb(3,3),qb(3,3),entropy
C      double precision plat(9)
C      logical llshft(3)
C#endif
C     for SX potential
      integer lrsig,lrsigl,ldham(16),lwsig
C..   new param for SX
      logical sxad
C.... Local parameters for LDA+U
      integer lmaxu,ng,ngi,nlibu
      double precision umix,tolu
      procedure(logical) :: bittst,cmdopt,parmxp
      procedure(integer) :: isw,nglob,iprint,a2vec,bitand,getdig,lgunit,
     .  fopn,fopno,fopna,fxst,fopng,lmfitmr1
      procedure(real(8)) :: dlength,dval

C ---------------------------- Setup ---------------------
C     call tcn('lmasa')
C ... Default values for mixing
      call dpzero(ehterm,6)
      call dpzero(dmxprm,34)
      call dpzero(vconst,3)
      allocate(vorb(1),dmatu(1),dmato(1))
      dmxprm(2) = s_mix%b
      dmxprm(9) = s_mix%bv
      dmxprm(4:6) = s_mix%w
      dmxprm(3) = s_mix%wc
      dmxprm(6) = s_mix%nsave
      dmxprm(7) = s_mix%mmix
      nbas = s_ctrl%nbas
      nclass = s_ctrl%nclass
      nl = s_ctrl%nl
      nspec = s_ctrl%nspec
      nsp = s_ctrl%nspin
      loptc = s_ctrl%loptc
      lasa = s_ctrl%lasa
      lsx = s_ctrl%lsx
      lrel = mod(s_ctrl%lrel,10)
C#ifdefC GF
CC ... If Dirac, unset the spin-orbit flag
C      if (lrel == 2 .and. bittst(s_ctrl%lncol,4)) then
C        s_ctrl%lncol = s_ctrl%lncol - 4
C        s_ham%lncol = s_ctrl%lncol
C      endif
C#endif
      lncol = s_ctrl%lncol
      lso = bittst(lncol,4)

      lham = s_ctrl%lham
      lpgf = s_ctrl%lpgf(1)
      neul = s_ham%neula
      sdmod = s_ctrl%sdmod
      sdprm(1:5) = s_ctrl%sdprm
      maxit = s_ctrl%maxit
      zbak = s_ctrl%zbak
      lcgf = s_ctrl%lcgf
      iscr = mod(s_ctrl%lscr,10)
      lrsig = s_ham%lsig
      npadl = s_ctrl%npadl
      npadr = s_ctrl%npadr
      nclasp = s_ctrl%nclasp
      nrhos = s_pot%nrhos
      avw = s_lat%avw
      lnsph = isw(IAND(s_ctrl%lasa,32) /= 0)
      lfree = IAND(s_ctrl%lasa,8) /= 0
      call supot(1,s_ctrl,s_lat,s_pot)
      vne = 0
      vrl = 0
      stdo = lgunit(1)
      procid = mpipid(1)
      master = 0
      mlog = cmdopt('--mlog',6,0,outs)
      ldlm = s_ctrl%ldlm
C#ifdefC GF
C      if ((lso .or. lrel == 2) .and. ldlm == 0)
C     .  call rx('Turn on CPA for spin-orbit coupling')
C#endif
      tdlm = s_ctrl%tdlm
      omonly = .false.

C ... popt file editor
      if (cmdopt('--popted',8,0,outs)) then
        call popted(outs(9:),s_ctrl,s_lat,s_bz)
        call rx0('lm from popted')
      endif

C ... Invoke the restart editor
      if (cmdopt('--rsedit',8,0,outs)) then
        nbasp  = nbas + npadl + npadr
        call rsedita(outs(9:),1,
     .    s_ctrl,s_site,s_spec,s_lat,s_pot,s_bz,nbasp,nbasp,nspec)
        call rx0('lmasa from rsedit')
      endif

C ... Input from restart file
      irs = 0; call readrssw(strn,i,irs)
      if (irs(1) /= 0) then
!     irs(1) = IAND(s_ctrl%lrs,7)
        if (lrel == 2) then
          call info0(10,0,0,' lmasa (warning):  --rs not fully implemented with lrel=2')
        endif
        if (procid == master) ifi = fopna('rsta',-1,0)
        call mpibc1(ifi,1,2,.false.,'lmasa','ifi')
        call asars(1,s_ctrl,s_site,s_spec,s_lat,s_pot,s_bz,
     .    s_pot%pnu,s_pot%qnu,.false.,ifi)
        call clsprp(1,s_ctrl,s_ham,s_pot,s_spec,s_lat,s_bz)
        if (procid == master) call fclr('rsta',ifi)
C       call shoctl(s_ctrl,s_spec,s_pot,0,stdo)
C       call rx('done')
      endif

C#ifdefC PGF
CC     for now
C      iv(3) = 0
C      iv(4) = 0
C      if (cmdopt('-nofg=',6,0,outs)) then
C        i = 6
C        i = a2vec(outs,len(outs),i,2,', ',2,2,2,iv,iv(3))
C      endif
C      s_ctrl%nofgl = iv(3)
C      s_ctrl%nofgr = iv(4)
C
C      if (lpgf == 2) iscr = 0
C#endif

      nangl = s_ctrl%nccomp
      ndlmcl = s_ctrl%ndlmcl
      if (ldlm == 0) nangl = 0
      omdiff = 0d0
      ldlm1  = mod(ldlm/10,10)
      lweiss = mod(ldlm1,2) == 1
      omonly = ldlm1 >= 2

      nbasp  = nbas + npadl + npadr
      nclspp = 2*nclasp-nclass
      nclspd = nclasp + nangl
      nclsppd = nclspp + nangl
      nlspcp = nl*nsp*max(nclsppd,nspec)
      nlspcr = 8*nl*nl*max(nclsppd,nspec)
C     nlspcd = nl*nsp*max(nclspp+nangl,nspec)

      nspc = 1
      if (bitand(lncol,1+2+4+8) /= 0) nspc = 2
      if (nspc == 2 .and. ldlm /= 0 .and. procid == master)
     .  print *,'DLM with coupled spins: IN DEVELOPMENT'
      nqpp   = (nl**2*(nl**2+1))/2
      lmfrce = bittst(lncol,16) .and. bittst(lncol,1)
      if (cmdopt('--pdos',6,0,outs)) then
        if (lmfrce) call info0(2,1,1,
     .    ' LMASA (warning): --pdos set; magnetic forces turned off')
        lmfrce = .false.
      endif
      if (cmdopt('-noves',6,0,outs)) s_ctrl%lves = IOR(s_ctrl%lves,1)
      lves   = 2*IAND(s_ctrl%lves,1)
      emad = 0
      trumad = 0
      seref = 0
      latok = 0
      call dpzero(s_pot%vdif,nclsppd)
C ... class-based arrays
      allocate(z(nclspd))
      call spec2class(s_spec,nclspd,s_ctrl%ics,'z',1,xx,z)
      allocate(lmx(nclspd))
      call spec2class(s_spec,nclspd,s_ctrl%ics,'lmxa',1,lmx,xx)
      allocate(idmod(nl,nclspd))
      call spec2class(s_spec,nclspd,s_ctrl%ics,'idmod',nl,idmod,xx)
      allocate(grp2(nclspd))
      call spec2class(s_spec,nclspd,s_ctrl%ics,'grp2',1,grp2,xx)

C --- Create initial file for scaling parameters ---
      if (procid == master .and. fxst('shfac') == 0 .and. IAND(s_ctrl%lham,512) /= 0) then
        allocate(shfac(3,nclsppd))
        call spec2class(s_spec,nclsppd,s_ctrl%ics,'shfac',3,xx,shfac)
        if (dlength(3*nclsppd,shfac,1) /= 0) then
          ifi = fopn('shfac')
          rewind ifi
          call info0(20,1,0,' LMASA:  writing initial shfac file')
C         call iomagf(2,nclsppd,1,shfac,shfac,1,-ifi)
          call ioextf(2,'shfac',nclsppd,1,shfac,shfac,3,1,-ifi)
          call fclr('shfac',ifi)
        endif
      endif

C --- Create initial file for vshft ---
C     2s bit of s_ctrl%lves should not be set by lmgf or lmpg
      if (IAND(s_ctrl%lves,2) /= 0 .and. procid == master) then
        if (fxst('vshft') == 0) then
          ifi = fopn('vshft')
          call iovshf(nbasp,-5,' ',0d0,0d0,0d0,s_pot%vshft,-ifi)
          call rx0('wrote vshft file')
        endif
      endif

C --- Read socscl file for spin-orbit scaling parameters ---
      call ptr_pot(s_pot,8+1,'socscl',nclasp+nangl,0,xx)
      if (lso) then
        call dvset(s_pot%socscl,1,size(s_pot%socscl),1d0)
        if (procid == master) then
        if (fxst('socscl') == 1) then
          ifi = fopn('socscl')
          rewind ifi
          call ioextf(2,'SO scaling from socscl',nclasp+nangl,1,s_pot%socscl,s_pot%socscl,1,1,ifi)
          call fclr('socscl',ifi)
          write(stdo,901) ! (s_pot%socscl(i),i=1,nclasp+nangl)
  901     format(' SOC coefficients scaled by: ':/16F6.2)
          call arrprt(' Class scale','%,5i%;7,3D','Id',nclasp+nangl,stdo,7,0,' | ',
     .      xx,s_pot%socscl,xx,xx,xx,xx,xx,xx)
        endif
        endif
        call mpibc1(s_pot%socscl,nclasp+nangl,4,.false.,'gfasa','socscl')
      endif

C --- Read shfac file for scaling parameters, incl XC field ---
      ifi = 0
      if (procid == master) then
        if (fxst('shfac') == 1) then
          ifi = fopn('shfac')
          rewind ifi
C         call info0(20,0,0,' GFASA:  reading initial shfac file')
C         call iomagf(2,nclass+nangl,1,s_pot%shfac,s_pot%shfac,1,ifi)
c         call ioextf(2,'constraining fields from shfac',nclass+nangl,1,s_pot%shfac,s_pot%shfac,3,1,ifi)
          call ioextf(2,'constraining fields from shfac',nclsppd,1,s_pot%shfac,s_pot%shfac,3,1,ifi)
          call fclr('shfac',ifi)
        endif
      endif
      call mpibc1(ifi,1,2,.false.,'lmasa','ifi')
c     if (ifi /= 0) call mpibc1(s_pot%shfac,(nclass+nangl)*3,4,.false.,'lmasa-gf','shfac')
      if (ifi /= 0) call mpibc1(s_pot%shfac,nclsppd*3,4,.false.,'lmasa-gf','shfac')

C     2s bit of s_ctrl%lves should not be set by lmgf or lmpg
      if (IAND(s_ctrl%lves,2) /= 0) then
        if (procid == master) then
          if (fxst('vshft') == 0) call rx('no file vshft')
          ifi = fopn('vshft')
          call iovshf(nbasp,4,'read potential shifts from file vshft',
     .      efermi,xv,vconst,s_pot%vshft,ifi)
        endif
        call mpibc1(s_pot%vshft,nbas,4,.false.,'lmasa','vshft')
      endif

C ... Other initialization
      call subasi(s_ctrl,s_spec,s_ham)
      nband = nbas*nglob('mxorb')
      allocate(eband(nband,nsp*s_bz%nkp))

C ... Setup for charge mixing
      call parms0(dmxprm(30),dmxprm(31),dmxprm(32),1)
      if (lmfrce) then
        call dpzero(dmxeu,34)
        call dcopy(25,dmxprm,1,dmxeu,1)
      endif
      if (cmdopt('--mix=',6,0,outs)) then
        i = 6
        iv(4) = 1
        i = a2vec(outs,len(outs),i,2,', ',2,2,2,iv,iv(3))
        iv(4) = max(iv(4),1)
        call parmx0(iv(3),iv(4),0d0)
      endif

C ... Green's function specific setup
C#ifdefC GF
CC ... setup for non-equilibrium mode
C      call setne(s_bz,s_ctrl,nzne,vne,lnoneq)
CC     if (lnoneq .and. (iscr /= 0 .or. lves /= 0) )
C      if (lnoneq .and. (lves /= 0) )
C     . call rx ('lmasa not ready for iscr!=0 or lves!=0 in non-eq mode')
C
C      i = 4
C      if (cmdopt('-ef=',i,0,outs)) then
C        if (.not. a2bin(outs,efermi,4,0,' ',i,-1)) call
C     .    rxs2('LMASA: failed to parse "',outs(1:30),'%a"')
C        xx = efermi(1) - s_bz%semsh(4)
C        call info2(10,1,1,' Override file fermi level, use ef= %,6d, '//
C     .    'delef= %,6d',efermi,xx)
C        i = mod(nint(s_bz%semsh(2)),100)
C        if (i /= 2) then
C          s_bz%semsh(3) = s_bz%semsh(3) + xx
C          s_bz%semsh(4) = s_bz%semsh(4) + xx
C        endif
C      endif
C
C      if (procid == master) then
C      i = 1+4
C#ifdefC PGF
C      i= 2+4
C#endifC
C        ifi = fopn('VSHFT')
C        xx = NULLI
C        call iovshf(nbasp,i,'read file vshft: ef=%d  vconst=%d',
C     .    xx,xx,vconst,s_pot%vshft,ifi)
C        s_pot%vconst = vconst
C        call fclr('VSHFT',-1)
C      endif
C      call mpibc1(vconst,3,4,.false.,'lmasa-gf','vconst')
C      call mpibc1(s_pot%vshft,nbas,4,.false.,'lmasa-gf','vshft')
C      call mpibc1(xx,1,4,.false.,'lmasa-gf','ef')
C#ifdefC PGF
C      vrl = vne + vconst(3)-vconst(2)
C#endifC
C
C      if (cmdopt('--fileef',8,0,outs) .and. xx /= NULLI) then
C        xx = xx - s_bz%semsh(4)
C        call info2(10,1,1,' Use file ef= %,6d, delef= %,6d',xx+s_bz%semsh(4),xx)
C        s_bz%semsh(3) = s_bz%semsh(3) + xx
C        s_bz%semsh(4) = s_bz%semsh(4) + xx
C      endif
C
CC ... Energy mesh
C      nzp = nint(s_bz%semsh(1))
C      allocate(zp(nzp+nzne))
C      allocate(wz(nzp+nzne))
C      call emesh(s_bz%semsh,zp,wz)
C
CC     debugging: read/write emesh on disk
C      if (cmdopt('--read-emesh',12,0,outs)) then
C        if (procid == master) then
C          ifi =  fopn('emesh2')
C          rewind ifi
C          read(ifi,*)nzp,nzne
C          print *,nzp,nzne
C          allocate(wk(2))
C          do  i = 1, nzp+nzne
C            read(ifi,822)wk(1),wk(2)
C            call dpscop(wk,zp,2,1,2*i-1,1d0)
C            print *,zp(i)
C          enddo
C          do  i = 1, nzp+nzne
C            read(ifi,822)wk(1),wk(2)
C            call dpscop(wk,wz,2,1,2*i-1,1d0)
C            print *,wz(i)
C          enddo
C          deallocate(wk)
C          call fclr('emesh2',ifi)
C        endif
C        call mpibc1(nzp,1,2,.false.,'lmasa-gf','nzp')
C        call mpibc1(nzne,1,2,.false.,'lmasa-gf','nzne')
C        s_bz%semsh(1) = nzp
C        s_bz%semsh(7) = nzne
C        call mpibc1(zp,(nzp+nzne)*2,4,.false.,'lmasa-gf','zp')
C        call mpibc1(wz,(nzp+nzne)*2,4,.false.,'lmasa-gf','wz')
C      endif
C
C      if (cmdopt('--write-emesh',13,0,outs)) then
C        if (procid == master) then
C          ifi = fopn('emesh')
C          rewind ifi
C          write(ifi,*)nzp,nzne
C          write(ifi,822)(dval(zp,i),i=1,(nzp+nzne)*2)
C          write(ifi,822)(dval(wz,i),i=1,(nzp+nzne)*2)
C  822     format(2F12.8)
C          call fclr('emesh',ifi)
C        endif
C      endif
C      npl = s_ctrl%npl
C
CC ... CPA and other site-specific setup
C      call dlminit(s_ctrl,s_ham,s_pot,s_spec,s_site,s_ctrl%ics,nbasp,nzp)
C      if (ldlm /= 0) then
Cc       tdlm = s_ctrl%tdlm
Cc       if (procid == master) then
Cc         if (lweiss) then
Cc           if (tdlm > 0) then
Cc            call awrit1('%N DLM: p(x) = -%,3d cos(t)',' ',80,stdo,tdlm)
Cc           else
Cc            call awrit1('%N DLM: p(x) = %,3d cos(t)',' ',80,stdo,-tdlm)
Cc           endif
Cc         elseif (tdlm /= 0d0) then
Cc           call awrit1('%N DLM: Spin T = %,0d K',' ',80,stdo,tdlm)
Cc         else
Cc           call awrit0('%N DLM: Infinite spin temperature',' ',80,stdo)
Cc         endif
Cc       endif
C
CC   ... Read/assign Gibbs weights for DLM
CC       ogibbs = odlmwt
C        s_pot%gibbs => s_pot%dlmwt
C        call pkgdlm(nbas,s_site,s_pot%dlmwt,s_pot%shfac)
CC   ... Assign class groups for DLM if necessary
C        call dlmgrp2(nclasp,s_ctrl%ncomp,s_ctrl%idcc,grp2)
C      else
C        allocate(s_ham%etrms(1))
C      endif
C
C#endif

C ------------------- Start of iterative loop ------------------
      makepp = mod(lasa,4) /= 0 .and. getdig(sdmod,3,10) == 0
      nwmoms = .false.
      nit = 0
   10 continue

C --- Madelung energy, potential, number of valence electrons ---
      call pshpr(1)
      allocate(dq(nclspd))
      call getq(nsp,nl,lmx,nclspd,z,s_pot%pnu,s_pot%qnu,
     .  s_ctrl%ics,s_spec,s_pot%qc,s_pot%qt,dq)
      call togpr
C#ifdefC GF
C      if (ldlm /= 0) then
C        call cpaz(nclasp,z,s_ctrl%ncomp,s_pot%gibbs)
C        call cpaz(nclasp,s_pot%qc,s_ctrl%ncomp,s_pot%gibbs)
C      endif
C#endif
      call getzv(nclasp,s_ctrl%nrc,z,s_pot%qc,s_bz,s_ctrl,zval)
      call togpr
      deallocate(dq)
      if (.not. associated(s_ham%etrms)) allocate(s_ham%etrms(1))

C ... Undo shift that was added to pp's when they were made
      call shftppr(lrel,nclspd,nl,nsp,s_pot%pp,s_pot%pprel,s_pot%ves,s_pot%ves,T,F)
      i = lves ; if (mod(irs(5),10) >= 4) i = 2 ! flags asamad to stick with given ves
      call asamad(s_ctrl,s_pot,s_lat,s_spec,i+100,s_pot%pnu,
     .  s_pot%qnu,vrl,s_pot%ves,emad,trumad,vmtz,s_ham%etrms)
C ... add shift for v as have now
      call shftppr(lrel,nclspd,nl,nsp,s_pot%pp,s_pot%pprel,s_pot%ves,s_pot%ves,F,T)
      call poppr

C --- Sphere program.  Make potential parameters ---
      imake = 0
      if (makepp) imake = 3
      if (bitand(lasa,3) >= 2 .and. nit == 0) imake = 2
      if (imake == 3 .and. lfree) imake = 1
      if (bitand(lasa,3) == 3) call rx('lmasa: not ready for BEGMOM=3')
C     Temporarily set lscr when making vintra only selected iterations
      iscsav = s_ctrl%lscr
      if (mod(iscsav,100) >= 20) then
        if (mod(nit,mod(iscsav,100)/10) /= 0) s_ctrl%lscr = 0
      endif
      if (irs(5) > 0) imake = imake + 10*irs(5)
      call asvsph(s_ctrl,s_lat,s_spec,s_ham,s_pot,vrl,imake,ehterm,latok)
      irs(5) = 0

      etrmss = s_ham%eterms
      if (nit == 0) then
        amgm = s_ham%eterms(15)
C        if (neul > 0) then
C          ifi = fopn('bxc')
C          rewind ifi
C          call info0(20,0,0,' LMASA:  writing initial bxc file')
C          call iomagf(1,nclass,1,s_pot%bxc,s_pot%bxc,1,-ifi)
C          call fclose(ifi)
C        endif
      else
        if (nsp == 2) call info2(30,0,0,' Magnetization from '//
     .    'sphere program: %1,6d  from output rho: %1,6d',
     .    s_ham%eterms(15),amgm)
      endif
      sevat = s_ham%eterms(16)
C     Reset lscr
      s_ctrl%lscr = iscsav

      fmode = s_tb%fmode
      if (fmode > 0) then
        if (ldlm /= 0) call rx0('LMASA-GF: fmode and DLM')
        if (.not. allocated(ppsave))
     .  allocate(ppsave(6,nl,nsp,max(nclspp,nspec)))
        call dcopy(6*nlspcp,s_pot%pp,1,ppsave,1)
      endif

C --- Remake empirical shifts in PPARS ---
      if (s_ham%rdvext /= 0) then
C     No shift unless new PP's were made
      if (imake == 2 .or. imake == 3 .or. nit == 0) then
C       No shift unless file present
        if (fxst('vext') == 1 .and. procid == master) then
C         sw = T
          ifi = fopn('vext')
          call pp2alp(1,s_ctrl,s_str,lham*0,s_pot%pp)
C         call prtcls(4,s_ctrl,s_pot,s_spec,1,nclass,stdo)
          i = s_ham%rdvext
          i = lmfitmr1(10*i+4,s_ctrl,s_site,s_spec,nl,nbasp,nsp,nclasp,
     .      z,iv,ldham,nfcof,nfvar,xv,s_pot%pp,xv)
C          call prtcls(4,s_ctrl,s_pot,s_spec,1,nclass,stdo)
C         call pp2alp(0,s_ctrl,s_str,lham*0,s_pot%pp)
C         call prtcls(4,s_ctrl,s_pot,s_spec,1,nclass,stdo)
        else
C         sw = F
          call info0(30,1,0,' LMASA (warning) missing file '//
     .      '"vext" required for HAM_RDVEXT ... skipping read')
        endif
      endif
C     call mpibc1(sw,1,1,.false.,'lmasa','sw_vext')
C     if (sw) call mpibc1(s_pot%pp,6*nlspcp,4,mlog,'lmasa','shifted pp')
      call mpibc1(s_pot%pp,6*nlspcp,4,mlog,'lmasa','shifted pp')
      endif

C ... Update second-padded pnu,qnu, in case they changed
C#ifdefC PGF
C       call clsprp(1,s_ctrl,s_ham,s_pot,s_spec,s_lat,s_bz)
C#endif

C ... Save convergence data in save file
      if (nit > 0 .and. procid == master) then
        call poseof(fopn('SV'))
        do i = 1, 3
          if (i /= 1 .or. iprint() >= 10) then
          ifi = lgunit(i)
          if (i == 3) ifi = fopn('SV')
          call iosv(-ifi,nit,dmxprm,ehk,amgm,.true.)
          endif
        enddo
        call fclose(fopn('SV'))
      endif

C --- Re-entry point for non-self-consistent calculations ---
   58 continue

C ... Output to restart file (eventually only for nit>0)
      if (cmdopt('--rs=',5,0,outs)) then
        irs(2) = IAND(s_ctrl%lrs,8+16)/8
        if (irs(2) > 0) then
          if (procid == master) then
            if (mod(irs(2),10) == 3) then
              outs = 'rsta.'
              i = 5
              call bin2a(' ',0,0,nit,2,0,len(outs),outs,i)
              ifi = fopng(outs,-1,0)
            else
              ifi = fopna('rsta',-1,0)
            endif
          endif
          call mpibc1(ifi,1,2,.false.,'lmasa','ifi')
          call asars(1,s_ctrl,s_site,s_spec,s_lat,s_pot,s_bz,
     .      s_pot%pnu,s_pot%qnu,.false.,-ifi)
          if (procid == master) call fclr('rsta',ifi)
        endif
      endif

C --- Print out summary of this iteration, check for convergence ---
      lsc = 3
      if (nit > 0) then
        xv(3) = ehf
        xv(4) = ehk
        qtol = s_ctrl%tol(1)
        etol = s_ctrl%tol(3)
        xv(1) = rmsdel
C       Independent mixing of Euler angles: combine rms dq, deuler
        if (getdig(sdmod,1,10) == 1) xv(1)=sqrt(rmsdel**2+dmxprm(32)**2)
        if (procid == master) then
          call nwit(s_ctrl%nvario,nit,maxit,.false.,1,etol,
     .      qtol,xv,0d0,0d0,'cxhi',amgm,xv(3),lsc)
        endif
        call mpibc1(lsc,1,2,.false.,'lmasa','lsc')
      elseif (iprint() >= 10 .and. makepp) then
C       Extract calculated energy parameters
        ehk0 = s_ham%ehk
        thrpv = s_ham%thrpv
        seref = s_ham%seref
        do i = 1, 2
          call awrit6('%x '//prgnam//'%a: it %i of %i  '//
     .      'ehk0=%1,6d  pv=%1;4d  mmom=%1;4d'//'  seref=%1;6d',outs,
     .      80,lgunit(i),nit,maxit,ehk0,thrpv/3,amgm,seref)
          if (ldlm /= 0) then
            call awrit2(' S = %1,6d  free energy = %1;6d',outs,
     .        80,lgunit(i),s_ham%entropy,
     .        ehk0 - (tdlm/11604.5d0/13.6057d0) * s_ham%entropy)
          endif
        enddo
      endif
      if (iprint() > 30) call cpudel(stdo,'...   Time this iter:',xv)

C ... Invoke shell command if supplied
      outs = ' '
      if (nit > 0 .and. cmdopt('--sh=',5,0,outs)) then
        call awrit0(' LM: invoke sh="'//outs(6:len(outs))//'%a"',
     .    ' ',len(outs),stdo)
        call fsystm(outs(6:len(outs)),i)
        call awrit1(' LM: shell returned value %i',' ',80,stdo,i)
      endif

C ... When quit is specified or insufficient information to continue
      if (s_ctrl%quit == 2 .or. getdig(latok,0,10) == 0 .or. bittst(lasa,8)) then
        call info0(0,0,0,' '//prgnam//'%a:  Q=ATOM encountered or missing input')
        goto 99
      endif
C ... Spin statics or dynamics: charge self-consistency not a criterion
      if (lsc == 0 .and. lmfrce .and. getdig(sdmod,0,10) > 1) lsc=3
C ... DLM case: omega must be converged, not just charge
      if (ldlm /= 0 .and. omdiff > 1d-3) lsc = 3
C ... Interactive query for number of iterations
      call ftflsh(stdo)
      if (lsc == 0) then
        if (procid == master) write(stdo,302) rms2
  302   format(/' Jolly good show !  You converged to rms DQ=',f10.6)
        goto 99
      endif
      call querym('maxit',2,maxit)
      if (nit >= maxit) goto 99
      nit = nit+1
C ... If subsequent iterations have moments in hand
      makepp = makepp .or. nwmoms
      makepp = makepp .and. getdig(sdmod,3,10) == 0

C ------------------------- Band pass --------------------------
C ... Keep a copy of input P,Q and related parameters
      if (.not. allocated(pold)) then
        allocate(pold(nl*nsp,max(nclsppd,nspec)),qold(3*nl*nsp,max(nclsppd,nspec)),vold(nclsppd))
        allocate(eulo(3*neul*nbasp))
        if (lrel == 2) then
          allocate(qrold(4*nlspcr))
        else
          allocate(qrold(1)); qrold(1) = NULLI; s_pot%qnur(1,1) = NULLI ! => arrays have no content
        endif
      endif
      call dcopy(nlspcp,s_pot%pnu,1,pold,1)
      call dcopy(3*nlspcp,s_pot%qnu,1,qold,1)
      call dcopy(nclsppd,s_pot%ves,1,vold,1)
      if (lrel == 2) call dcopy(4*nlspcr,s_pot%qnur,1,qrold,1)

      if (lnsph /= 0) then
        if (ldlm > 0) call rx0('DLM not ready for lnsph')
        if (.not. allocated(qppo)) allocate(qppo(nqpp*4*nsp*nbas))
        call dcopy(2*nqpp*4*nsp*nbas,s_pot%qpp,1,qppo,1)
      else
        if  (.not. allocated(qppo)) allocate(qppo(1))
      endif

      call suham(s_ctrl,s_lat,s_spec,s_site,s_ham,s_pot,s_str,s_strn)

C ... Read and store self-energy addition to LDA potential
      ldham = s_ham%ldham
      lrsigl = 0
      if (lrsig /= 0 .and. procid == master) then
        lrsigl = fxst('sigm')
      endif
      call mpibc1(lrsigl,1,2,.false.,'lmasa','lrsigl')
      if (lrsig /= 0 .and. lrsigl == 0) then
        call info0(10,1,0,' lmasa (warning): no sigm file found ... LDA calculation only')
        s_ham%lsig = 0
        lrsigl = 0
      else
        lrsigl = lrsig
      endif
      sxad = F

      if (mod(lrsigl,10) /= 0) then
C       Real-space range
        ifi = fopna('sigm',-1,4)
        call rdsigm(1000+lrsigl,nbas,nsp,ldham(1),s_ctrl,s_spec,s_lat,
     .    s_ham,s_bz,s_gw,ifi,lwsig)
        call fclose(ifi)
        sxad = T
      endif

C ... LDA+U initialization
      if  (.not. allocated(lldau)) allocate(lldau(nbasp))
      call iinit(lldau,nbasp)
C     Check for LDA+U ... return nlibu > 0 if any U blocks.
      call suldau(nbasp,s_spec,s_site,nlibu,lmaxu,lldau)
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
C       layer GF code doesn't have istab for padding layers.
        ngi = ng
        if (lpgf /= 0) then
          ngi = 0
        endif
        umix = s_mix%umix
        tolu = s_mix%tolu

C       initialize density matrix for LDA+U
        i = 0
        if (bittst(lham,256)) i = 1
        call sudmtu(nbasp,nsp,nlibu,lmaxu,s_site,s_spec,i,lldau,
     .    ngi,s_lat%symgr,s_lat%istab,dmatu,vorb)

C       Hang on to previous site density matrix for this iteration
        i = nsp*nlibu*(lmaxu*2+1)**2
        call dcopy(2*i,dmatu,1,dmato,1)
        call dpzero(dmatu,2*i)
      endif

C --- Levenberg-Marquardt Band fitting mode: setup ---
      fmode = s_tb%fmode
      fmode0 = mod(fmode,10)
      fmode1 = mod(fmode/10,10)
C ... Standard mode, no LM fitting
      if (fmode0 <= 0) then
        nfvar = 0
        if (.not. allocated(gband)) then
          allocate(gband(1,1,1))
        endif
        if (.not. allocated(ivcof)) allocate(ivcof(1,1))
C ... Levenberg-Marquardt fitting mode
      else
      nbfit = s_tb%nbfit
      nbfitf = s_tb%nbfitf
      fitmrw = s_tb%wt
      ftalam = s_tb%alam
      ftalsc = s_tb%alsc
      if (fmode0 == 1 .or. fmode0 == 2) then
        call info2(10,1,0,' Levenberg-Marquardt mode %i:'//
     .    '  fit ASA C%?#(n>1)#+enu##,Delta to data',fmode0,fmode0)
      else
        call rx1('LMASA: fitting mode %i not implemented',fmode)
      endif

C     Count the number of fitting parameters
      ldham = s_ham%ldham
      call pshpr(1)
      i2 = 0
      if (s_ham%rdvext /= 0 .and. fxst('vext') == 1) i2=5
      if (procid == master) then
        i = lmfitmr1(10*fmode0+i2,s_ctrl,s_site,s_spec,nl,nbasp,nsp,nclasp,
     .    z,s_ham%iprmb,ldham,nfcof,nfvar,xv,xv,xv)
      endif
      call mpibc1(nfcof,1,2,.false.,'lmasa','i2')
      call mpibc1(nfvar,1,2,.false.,'lmasa','i2')
      call poppr

C ... Setup for fitting mode 1 (ASA C and Delta)
      if (fmode0 == 1 .or. fmode0 == 2) then
        call info5(10,0,0,' Fit %i parameters in %i Rl channels.',
     .    nfvar,nfcof/2,0,0,0)

        allocate(ivcof(nfcof,2),ftcof(nfcof))
        call iinit(ivcof,nfcof*2)
        call pp2alp(1,s_ctrl,s_str,lham*0,s_pot%pp) ! pp's in orthogonal rep
C       Read ivcof,ftcof
        if (procid == master) then
        i = lmfitmr1(10*fmode0+1+i2,s_ctrl,s_site,s_spec,nl,nbasp,nsp,nclasp,
     .      z,s_ham%iprmb,ldham,nfcof,nfvar,ivcof,s_pot%pp,ftcof)
        endif
        call pp2alp(0,s_ctrl,s_str,lham*0,s_pot%pp) ! pp's in tb rep
        call mpibc1(ivcof,nfcof*2,2,.false.,'lmasa','ivcof')
        call mpibc1(ftcof,nfcof,4,.false.,'lmasa','ftcof')
        allocate(ftalp(nfcof,nfcof),ftcov(nfcof,nfcof),ftwk(nfcof,5))

        call info0(30,0,0,' LMASA:  saving initial vext0')
        call pshpr(10)
        if (procid == master) then
        i = lmfitmr1(10*fmode0+3,s_ctrl,s_site,s_spec,nl,nbasp,nsp,nclasp,
     .      z,s_ham%iprmb,ldham,nfcof,nfvar,ivcof,ppsave,ftcof)
        endif
        call poppr
C       call prmx('pp',s_pot%pp,6,6,nl*nsp*nclasp)
      endif
      if (nfvar <= 0) call rx('LMFITMR: no parameters to vary')

C ... Set up input for reference data
C ... Fit bands: check whether bands file can be read; get dimensions
      if (fmode1 == 0) then
        call info2(10,1,0,' Data mode %i: Fit to bands',fmode1,0)
        if (nbfit(2) == 0) nbfit(2) = ldham(1)
        if (procid == master) then
        ifi = fopno('refbnds')
        call pshpr(1)
C       Determine the number of bands and qp in reference file
        call suqlsr(21,ifi,max(nspc,nsp),nbf,i2,1,1,1,.false.,
     .    .false.,nqf,efermi(2),xv,xv)
        call poppr
        endif
        call mpibc1(nqf,1,2,.false.,'lmasa','nqf')
        call mpibc1(nbf,1,2,.false.,'lmasa','nbf')
        call mpibc1(efermi(2),1,4,.false.,'lmasa','efermi')

        call info8(10,0,0,
     .    ' bands file file contains nq=%i, nband=%i  efermi=%d'//
     .    '%N Fit to %i bands'//
     .    '%?#(n>1)#%-1j (%i to %i)#%j#.'//
     .    '%?#(n>1)#%-1j  Fit to file data, bands %i to %i.#%j#',
     .    nqf,nbf,efermi(2),
     .    nbfit(2)-nbfit(1)+1,nbfit(1),nbfit(2),
     .    nbfitf,nbfitf+nbfit(2)-nbfit(1))

C       Last band cannot exceed hamiltonian dimension
        call sanrg(.true.,nbfit(2),nbfit(1),ldham,'LMASA','2nd arg ITER_FIT_NBFIT')
C       Last band cannot exceed number of bands available
        if (nbfit(2)-nbfit(1) > nbf-nbfitf)
     .    call rx('LMASA: file does not have sufficient number of '//
     .    'bands to fit')

C       Read file energy bands for fitting, and file Fermi level
        nbf = nbfit(2)-nbfit(1)+1
        allocate(ebf(nbf,nqf*nsp),sigmrq(nbf,nqf*nsp))
        call dpzero(ebf,nbf*nqf)
        if (procid == master) then
        call suqlsr(25,ifi,max(nspc,nsp),i1,i2,nbf,nbfitf,
     .    nbfitf+nbfit(2)-nbfit(1),.false.,.false.,nqf,efermi(2),xv,ebf)
        endif
        call mpibc1(ebf,nbf*nqf,4,.false.,'lmasa','ebf')

C       If Fermi level is required, error if it was not read
        if (s_tb%shft(1) == 1) then
          if (efermi(2) == NULLI) call
     .      rx('LMASA: FIT_SHFTT=1 requires Fermi level in bands file')
        endif
C       Number of points to fit
        nfit = nbf*nqf*nsp

C       Fitting weights
        if (nint(fitmrw(1)) == 0) then
          call dvset(sigmrq,1,nfit,1d0)
        else
          if (efermi(2) == NULLI) call
     .      rx('LMASA: FIT_WT>0 requires Fermi level in bands file')
          do  i = 1, nfit
            iv(1) = mod(i-1,nbf)+1; iv(2) = (i-1)/nbf+1
            xv(1) = (ebf(iv(1),iv(2))-fitmrw(2)-efermi(2))/fitmrw(3)
            if (nint(fitmrw(1)) == 1) then
              sigmrq(iv(1),iv(2)) = dsqrt((exp(xv(1))+1))
            elseif (nint(fitmrw(1)) == 2) then
              sigmrq(iv(1),iv(2)) = dsqrt((exp(xv(1)**2)+1))
            elseif (nint(fitmrw(1)) == 3) then
              sigmrq(iv(1),iv(2)) = dsqrt(dexp(dabs(xv(1))))
            elseif (nint(fitmrw(1)) == 4) then
              sigmrq(iv(1),iv(2)) = dsqrt(dexp(dabs(xv(1)**2)/2))
            else
              call rx1('LMASA: fitting weights mode %d not implemented',
     .          fitmrw)
            endif
          enddo
        endif

C       Allocate gband
        nband = nbas*nglob('mxorb')
        allocate(gband(nfvar,nband,nqf*nsp))

C ... Fit partial DOS: get channels, read reference dos
      else if (fmode1 == 1) then
        xv(1) = s_tb%wg(1)
        xv(2) = s_tb%wg(2)
        call info5(10,1,0,' Data mode %i: Fit to partial DOS'//
     .  '%?#n#, smooth with wg=%g#, not smoothed%j#'//
     .  '%?#n#, ref wg=%g#%j#',
     .    fmode1,isw(xv(1) /= 0),xv,isw(xv(1) /= xv(2)),xv(2))

        if (cmdopt('--pdos',6,0,outs)) then
          i1 = min(1024,nbas*nl**2)
          allocate(lsite(nbas))
          nsgrp = s_lat%nsgrp
          i = 0
          call sumlst(0,0,nbas,nsgrp,s_site,s_spec,outs(7:),
     .      moddos,nsite,lsite,i1,nchan,i2,xx,i)
          deallocate(lsite)
        else
          call rx('LMASA: LM fit to dos requires --pdos switch')
        endif

C       First pass just gets number of energies
        if (procid == master) then
          ifi = fopno('refdos')
          call iodos(1,ifi,xv,xv,xv,i1,iv(2),xv,xv(2),iv(3),xv(3),1)
          call dvset(dosprm,1,6,dble(NULLI))
          dosprm(3) = i1
        endif
        call mpibc1(i1,1,2,.false.,'lmasa','dosprm')
        allocate(refdos(i1,nchan*nsp))
C       Second pass reads dos; requires match in ndos,nchan,nspin
        if (procid == master) then
          rewind ifi
          call iodos(13,ifi,refdos,i1,nchan*nsp,i1,nchan,dosprm(1),
     .      dosprm(2),nsp,efermi(2),1)
          call info5(10,0,0,' Reference DOS contains %i channels, '//
     .      '%i points in (%d,%d)  Ef=%d',
     .      nchan,i1,dosprm(1),dosprm(2),efermi(2))
          call fclr('refdos',ifi)
        endif
        call mpibc1(refdos,i1*nchan*nsp,4,mlog,'lmasa','refdos')
        call mpibc1(dosprm,6,4,.false.,'lmasa','dosprm')
        call mpibc1(efermi,2,4,.false.,'lmasa','efermi')
        dosprm(4:5) = s_tb%ebfit
        i1 = s_tb%ndos
        dosprm(6) = i1
        if (dosprm(4) == NULLI) call dcopy(2,dosprm,1,dosprm(4),1)
        if (dosprm(6) == NULLI) dosprm(6) = dosprm(3)
        i1 = int(dosprm(6))
        call info5(10,0,0,' Generate fit DOS on '//
     .    '%i points in interval (%d,%d)',i1,dosprm(4),dosprm(5),0,0)

C       Number of points to fit
        nspx = nsp / nspc
        nfit = i1*nchan*nspx
        allocate(sigmrq(i1,nchan*nspx))

C       Allocate gdos
        allocate(gdos(nfvar,i1,nchan*nspx))
        allocate(gband(1,1,1))  ! So compilers do not complain

      else
        call rx1('LMASA: fitting mode %i not implemented',fmode)
      endif

      call info5(10,0,0,' '//
     .  '%?#(n==0)#Unit##%-1j'//
     .  '%?#(n==1)#Fermi##%-1j'//
     .  '%?#(n==2)#Gaussian-like fermi##%-1j'//
     .  '%?#(n==3)#Exponential##%-1j'//
     .  '%?#(n==4)#Gaussian##%-1j'//
     .  '%?#(n>4)#Unknown## weights%-1j'//
     .  '%?#(n>0&n<3)#:  cutoff=Ef+%d width=%d%-2j##%-1j'//
     .  '%?#(n>2&n<5)#:  center=Ef+%d width=%d%-2j##%-1j'//
     .  ' ',
     .  nint(fitmrw(1)),fitmrw(2),fitmrw(3),0,0)

      if (fmode1 == 0) s_bz%nevmx = -1
      nit = 0 ! zeroth iteration for L-M fitting

      endif  ! End of Levenberg-Marquardt fitting setup

C --- Green's function approaches ---
C#ifdefC GF
C      if (lpgf /= 0) then
C#ifdefC PGF
CCC       Printout vshft
CC        call iozshf(nbasp,1,s_pot%vshft,s_ctrl%ipc,grp2,-stdo)
CCC       Replicate vshft->vshfs; distribute averages for pgfasa
CCC       Now done in pgfasa
CCC       call vshfav(11,npl,s_ctrl%pgplp,vshfs,s_pot%vshft)
C        call pgfasa(s_ctrl,s_lat,s_spec,s_site,s_bz,s_str,s_ham,
C     .    s_pot,s_strn,lldau,vorb,dmatu,vconst,s_pot%vshft,
C     .    zp,wz,s_ctrl%pgfsl,s_ctrl%pgplp,efermi,moddos,sumev,
C     .    s_pot%qnu,s_pot%rhos,amag,s_pot%aamom)
C        if (mod(mod(moddos,10),2) /= 1) goto 99
C        nevmx = 1
CC       call shoctl(s_ctrl,s_spec,s_pot,0,stdo)
C
CC       Add vshft(ib) to vdif to make etot when Ves ne Ves[n]
C        if (lves == 2) call
C     .    vsh2es(1,nbasp,s_ctrl%ipc,s_ctrl%nrc,s_pot%vshft,s_pot%vdif)
CC       Update vshft file
C        if (lpgf <= 2) then
C          ifi = fopn('VSHFT')
C          call iovshf(nbasp,6,' ',efermi,efermi,vconst,s_pot%vshft,-ifi)
C          call fclr('VSHFT',-1)
C        endif
C
C#endifC
C#ifdefC CGF
C      elseif (lcgf /= 0) then
C
CC       Temporarily set lscr when making psta only selected iterations
C        iscsav = s_ctrl%lscr
C        if (mod(iscsav,1000) >= 200) then
CC          if (mod(nit,mod(iscsav,1000)/100) /= 0)
C          if (mod(nit,mod(iscsav,1000)/100) /= 0)
C     .      s_ctrl%lscr = 0
C        endif
C        call dpcopy(vconst,vcnsts,1,3,1d0)
CC   ... Read Omega from disk for DLM
C        if (ldlm > 0 .and. nit == 1) then
C          call initomg(s_ctrl,s_site,nbas,nspc,nzp)
C          call rdomg(100*(nspc-1),s_site,nbas,nzp)
C        endif
CC       if (mod(ldlm/100,10) == 1) call rdomg(1,s_site,nbas,nzp)
C        if (.not. associated(s_pot%aamom)) s_pot%aamom => xxp
C        call gfasa(s_ctrl,s_lat,s_spec,s_site,s_bz,s_ham,s_pot,s_str,s_strn,lldau,vorb,dmatu,
C     .    vconst,s_pot%vshft,zp,wz,s_ctrl%pgplp,efermi,moddos,sumev,ehterm,
C     .    amag,s_pot%aamom,sxad,omdiff)
C
C        if (omonly) goto 99
C
CC       Reset lscr
C        s_ctrl%lscr = iscsav
C
C        if (mod(mod(moddos,10),2) /= 1) goto 99
C        nevmx = 1
CC       call shoctl(s_ctrl,s_spec,s_pot,0,stdo)
C
CC       Add vshft(ib) to vdif to make etot when Ves ne Ves[n]
C        if (lves == 2) call vsh2es(1,nbasp,s_ctrl%ipc,s_ctrl%nrc,s_pot%vshft,s_pot%vdif)
C
CC       Update vshft file
C        if (procid == master) then
C          ifi = fopn('VSHFT')
C          call iovshf(nbas,5,' ',efermi,efermi,vconst,s_pot%vshft,-ifi)
C          call fclr('VSHFT',-1)
C        endif
C
C
CC  ... Check convergence of dmatu and update it and vorb if necessary
C        if (nlibu > 0 .and. maxit > 0) then
C        i = 0
C        if (bittst(lham,256)) i = 1
C        call chkdmu(nbas,nsp,nlibu,lmaxu,s_site,s_spec,s_ham,i,dmatu,dmato,
C     .    vorb,tolu,umix,lldau,ng,s_lat%symgr,s_lat%istab)
C        endif
C
C#endifC
C      else
CC        if (mordrn == EMCLUS) then
CC          call rx('emc not ready')
CC          call emc2C(s_ordn,nl,nlo,nsp,nbas,lmx,nclass,s_ctrl%ipc,
CC     .      s_ham%eula,switch,s_pot%pp,vmtz,elin,wsr,
CC     .      efmax,nevmx,z,s_pot%qc,eband,nband,
CC     .      norder,width,range,npts,drange,
CC     .      avw,zval,plat,w(obas),s_ctrl%group,efermi,sumev,s_pot%rhos)
CC        endif
C        call rx('no technique specified to generate density')
C      endif
C
C#else
C --- Energy bands by direct diagonalization ---
C     Suppress ccor if generating evecs for response function
C     if (mod(lsx,2) == 1 .or. mod(iscr,2) == 1)
      if (mod(lsx,2) == 1 .or. mod(iscr,2) == 1)
     .  s_ctrl%lasa = lasa-bitand(lasa,4)

C ... Quit if --quit=ham given
      if (s_ctrl%quit == 8) call rx0('quit = ham')

C --- L-M fitting setup ---
   40 continue
      if (fmode0 > 0) then
        if (cmdopt('--band',6,0,outs)) then
          call rx('LMASA: cannot use L-M fitting with --band')
        endif
C   ... Get Fermi level if required
        if (s_tb%shft(1) == 1) then
          call info0(20,1,0,' LMASA: Find Ef for given parameters ...')
     .
          deallocate(eband)
          allocate(eband(nband,nsp*s_bz%nkp))
C         Delete weight file, if it exists
          ifi = fopna('wkp',-1,4)
          call dfclos(ifi)
          if (.not. associated(s_pot%aamom)) s_pot%aamom => xxp
          call bndasa(s_ctrl,s_lat,s_spec,s_site,s_bz,s_ham,
     .      s_pot,s_str,s_optic,vorb,dmatu,efermi,-1,0,0,0,xx,
     .      gband,eband,nband,nevmx,s_pot%qnu,sumev,s_pot%rhos,amag,
     .      s_pot%aamom)
        endif

C   ... Levenberg-Marquardt band fitting mode
        if (fmode1 == 0) then
          call info0(20,1,0,
     .      ' LMASA: Make gradient of bands wrt parameters ...')
          deallocate(eband)
          allocate(eband(nband,nsp*nqf))

C         call prmx('before pp',s_pot%pp,6,6,nl*nsp*nclasp)
          call bndasa(s_ctrl,s_lat,s_spec,s_site,s_bz,s_ham,
     .      s_pot,s_str,s_optic,vorb,dmatu,xv,-1,fmode0,nfvar,
     .      nfcof,ivcof,gband,eband,nband,nevmx,s_pot%qnu,sumev,
     .      s_pot%rhos,amag,s_pot%aamom)
C         call prmx('after pp',s_pot%pp,6,6,nl*nsp*nclasp)

          allocate(ewk(nbf,nqf*nsp),gwk(nfvar,nbf,nqf*nsp))
C         ewk(1:nbf,1:nqf) = eband(nbfit(1):nbfit(2),:)
          call dmcpy(eband(nbfit(1),1),nband,1,ewk,nbf,1,nbf,nqf*nsp)
C         gwk(:,1:nbf,:) = gband(:,nbfit(1):nbfit(2),:)
          call dmcpy(gband(1,nbfit(1),1),nfvar*nband,1,gwk,nfvar*nbf,1,
     .      nfvar*nbf,nqf*nsp)

C     ... Shift fit energy bands by constant to align Fermi levels
          if (s_tb%shft(1) == 1) then
            xv(1) = efermi(2) - efermi(1)
C           do  i2 = 1, nqf
C           do  i1 = 1, nbf
C             ewk(i1,i2) = ewk(i1,i2) + xv(1)
C           enddo
C           enddo
            ewk = ewk + xv(1)

            call info2(30,1,0,' LM fit:  '//
     .        'Shift bands by %d Ry to align Fermi levels ...',xv,0)
            if (iprint() > 30) then
              call info0(30,0,0,' Bands at first qp:%N ref')
              write(stdo,444) ebf(:,1)
              call info0(30,0,0,' fit')
              write(stdo,444) ewk(:,1)
              call info0(30,0,0,' sig')
              write(stdo,444) sigmrq(:,1)
  444         format(9f11.5)
              if (nsp == 2) then
              call info0(30,0,0,' Spin 2:%N ref')
              write(stdo,444) ebf(:,2)
              call info0(30,0,0,' fit')
              write(stdo,444) ewk(:,2)
              call info0(30,0,0,' sig')
              write(stdo,444) sigmrq(:,2)
              endif
            endif
          endif

C     ... chi^2 for current values of coffs
          call mrqcof(nfit,nfvar,nfcof*0,sigmrq,ebf,ewk,gwk,
     .      ftchi(3),xv,xv)
        endif

C   ... Levenberg-Marquardt DOS fitting mode
        if (fmode1 == 1) then
          call info0(20,1,0,
     .      ' LMASA: Make gradient of DOS wrt parameters ...')

C     ... Make the DOS
          if (procid == master) nfilem = fopna('moms',-1,4)
          i1 = int(dosprm(6))
          allocate(fitdos(i1,nchan*nspx))
          call gtpdss(1,nfilem,s_bz,.false.,nchan,nsp,
     .      i1,dosprm(4),dosprm(5),s_tb%wg(1),fitdos)

C     ... Setup shifted endpoints, number for reference and fit dos
          call dcopy(6,dosprm,1,dosprm(7),1)
          dosprm(9) = 0
C         Shift reference dos by constant to align Fermi levels
          if (s_tb%shft(1) == 1)
     .      dosprm(9) = efermi(1)-efermi(2)

C     ... Make shifted and scaled reference dos
          allocate(refdos2(i1,nchan*nsp))
          call lmrefdos(16,int(dosprm(3)),dosprm(1),dosprm(2),
     .      int(dosprm(6)),dosprm(4),dosprm(5),nchan*nspx,efermi,
     .      dosprm(9),fitmrw,s_tb%wg(2),
     .      refdos,fitdos,refdos2,sigmrq)

C         Debugging printout
          if (procid == master) then
          call info0(30,0,0,
     .      ' LMASA: saving fit and reference DOS, files dosf, dosr')
          i1 = int(dosprm(6))
          call iodosw(-fopn('dosr'),refdos2,sigmrq,i1,i1,nchan,
     .      dosprm(4),dosprm(5),nspx,efermi)
          call fclose(fopn('dosr'))
          i1 = int(dosprm(6))
          call iodosw(-fopn('dosf'),fitdos,sigmrq,i1,i1,nchan,dosprm(4),
     .      dosprm(5),nspx,efermi)
          call fclose(fopn('dosf'))
          endif

C     ... Make gradient of DOS wrt parameters
          call lmgdos(s_ctrl,s_lat,s_spec,s_site,s_bz,s_ham,
     .      s_pot,s_str,s_optic,vorb,dmatu,fmode0,
     .      nfvar,nfcof,ivcof,dosprm,s_tb%wg(2),nchan,nspc,
     .      eband,nband,s_pot%qnu,s_pot%rhos,gdos)

C     ... chi^2 for current values of coffs
          call mrqcof(nfit,nfvar,nfcof*0,sigmrq,refdos2,fitdos,gdos,
     .      ftchi(3),xv,xv)
        endif

        if (dsqrt(ftchi(3)/nfit) < 1d-10) ftalam = 0
        if (maxit > 0 .and. nit == maxit) ftalam = 0

C       Printout, Zeroth iteration
        if (nit == 0) then
          ftchi(4) = ftchi(3)
          call dpzero(ftchi,2)
          call info2(10,1,0,' LM fit, starting RMS error=%1,3;3g'//
     .      '  lambda=%;3g',dsqrt(ftchi(3)/nfit),ftalam)
        else if (nit > 0) then
          call info5(10,1,0,'  LMFITMR:  iter %i'//
     .      '  old RMS=%1,3;3g  new RMS=%1,3;3g  lambda=%;3g',
     .      nit,dsqrt(ftchi(1)/nfit),dsqrt(ftchi(3)/nfit),ftalam,0)
        endif

        if (ftalam > 0) then
        if (fmode1 == 0) then
          call lmrqit(10*fmode0+2,s_ctrl,s_site,s_spec,nl,nbasp,
     .      nsp,nclasp,z,s_ham%iprmb,ldham,s_pot%pp,nfit,nfcof,
     .      nfvar,ivcof,ebf,ewk,gwk,sigmrq,ftwk,ftcof,ftalsc,ftalp,
     .      ftalam,ftcov,ftchi,nit)
        elseif (fmode1 == 1) then
          call lmrqit(10*fmode0+2,s_ctrl,s_site,s_spec,nl,nbasp,
     .      nsp,nclasp,z,s_ham%iprmb,ldham,s_pot%pp,nfit,nfcof,
     .      nfvar,ivcof,refdos2,fitdos,gdos,sigmrq,ftwk,ftcof,ftalsc,
     .      ftalp,ftalam,ftcov,ftchi,nit)
        endif
        endif

C       Printout, final iteration
        if (ftalam <= 0) then
          call info5(10,1,0,'  LMFITMR:  iter %i'//
     .      '  starting RMS=%1,3;3g  final RMS=%1,3;3g  lambda=%;3g',
     .      nit,dsqrt(ftchi(4)/nfit),dsqrt(ftchi(3)/nfit),ftalam,0)
          call pp2alp(1,s_ctrl,s_str,lham*0,ppsave) ! save pp's in orthogonal rep
          if (procid == master) then
          call info0(10,0,0,'%12psaving file vext0')
          i = lmfitmr1(10*fmode0+3,s_ctrl,s_site,s_spec,nl,nbasp,nsp,nclasp,
     .      z,s_ham%iprmb,ldham,nfcof,nfvar,ivcof,ppsave,ftcof)
          endif
C         call pp2alp(0,s_ctrl,s_str,lham*0,ppsave) ! pp's in alpha rep
          call rx0('finished fit')
        endif

        if (allocated(ewk)) then
          deallocate(ewk,gwk)
        endif
        if (allocated(fitdos)) then
          deallocate(fitdos,refdos2)
        endif
        goto 40

      endif

C ... Band pass, normal mode
      if (.not. associated(s_pot%aamom)) s_pot%aamom => xxp

      call bndasa(s_ctrl,s_lat,s_spec,s_site,s_bz,s_ham,s_pot,
     .  s_str,s_optic,vorb,dmatu,efermi,nit,0,0,nfcof,ivcof,
     .  gband,eband,nband,nevmx,s_pot%qnu,sumev,s_pot%rhos,amag,
     .  s_pot%aamom)

C ... Check convergence of dmatu and update it and vorb if necessary
      if (nlibu > 0 .and. maxit > 0 .and. nevmx > 0) then
        i = 0
        if (bittst(lham,256)) i = 1
        call chkdmu(nbas,nsp,nlibu,lmaxu,s_site,s_spec,s_ham,i,
     .    dmatu,dmato,vorb,tolu,umix,lldau,ng,
     .    s_lat%symgr,s_lat%istab)
      endif

      s_ctrl%lasa = lasa

C --- Rucker's screened exchange ---
C#ifdefC SX
C      if (mod(lsx,2) == 1 .or. mod(iscr,2) == 1) then
C
C        aintra = T
C        zdel(1) = 0
C        zdel(2) = 0
CC       sxopts = ' '
CC       if (i2 >= i1) sxopts = sstrn(i1:i2)
C        i1 = str_pack('sxopt',-2,s_strn,sxopts)
C        if (mod(lsx,2) == 0) sxopts = 'pstat'
C        if (cmdopt('-novintr',4,0,outs))
C     .    call dpzero(s_pot%vintr,nclasp*nl*nl*nsp*nsp)
C        call asasx(s_ctrl,s_spec,s_site,s_lat,s_ham,s_pot,s_bz,
C     .    s_str,s_strn,sxopts,aintra,nl,nsp,efermi,eband,nband,zdel,
C     .    s_pot%vintr)
CC        call rx0('finished calling asasx ('//sxopts//'%a)')
C        goto 99
C      endif
C#endif (SX branch)

      if (nevmx > 0 .and. procid == master) then
        call shoctl(s_ctrl,s_spec,s_pot,0,fopn('LOG'))
      endif

C#ifdefC STONER
C      if (IAND(s_ctrl%lstonr(1),1) /= 0) then
C        real(8), allocatable :: nbar(:)
C        real(8), allocatable :: ewk(:)
C        real(8), allocatable :: mwk(:)
C        real(8), allocatable :: amom(:)
C        real(8), allocatable :: emag(:)
C        allocate(nbar(mnpts))
C        allocate(ewk(mnpts))
C        allocate(mwk(mnpts))
C        allocate(amom(nbas)); call dpzero(amom,nbas)
C        allocate(emag(nbas)); call dpzero(emag,nbas)
C        call rx('fix clabl for call to stoner')
C        call stoner(nl,nclass,clabl,w(oidxdn),efermi,
C     .    drange,iabs(npts),w(ozos),w(oindex),w(ostni),w(oammx),
C     .    mnpts,switch(44),w(onbar),w(oewk),w(omwk),w(oamom),
C     .    w(oemag))
C      endif
C#endif

C#endif    (hamiltonian or GF branches)

C --- Harris and HK energies if terms are available ---
C ... HF energy
      ehf = 0
      if (getdig(latok,1,10) /= 0) call asetot(1,s_ham,sumev,ehf)
C ... HK energy
      ehk = 0
      if (nevmx > 0) then
        call pshpr(max(iprt(1)-100,0))
        call pshpr(max(iprt(1)-20,0))
        call asvsph(s_ctrl,s_lat,s_spec,s_ham,s_pot,vrl,4,ehterm,latok)
C       Moment from output density
        amgm = s_ham%eterms(15)
C PGF doesn't make amag yet
C#ifndef PGF
        if (neul > 0 .and. nit > 0) amgm = dsqrt(amag(1)**2 + amag(2)**2 + amag(3)**2)*nbasp
C#endif
        call poppr
        call poppr
        if (getdig(latok,1,10) /= 0) call asetot(2,s_ham,sumev,ehk)
      endif

      if (ehf /= 0) then
        do  i = 1, 2
            if (iprint() >= 10 .or. i == 2.and.procid == master) then
            call awrit6('%N '//prgnam//'%a: '//
     .          'ehf=%1,7;7d  %?#n#ehk=%1,7;7d  #%j#sumev=%1,7;7d'//'%?#n#  delsev=%1,7;7d',
     .          outs,100,lgunit(i),ehf,isw(ehk /= 0),ehk,sumev,isw(sevat /= 0),sumev-sevat)
            if (ldlm /= 0) then
              call awrit2(' S = %1,6d  free energy = %1,7;7d',outs,80,lgunit(i),s_ham%entropy,
     .          ehk - (tdlm/11604.5d0/13.6057d0) * s_ham%entropy)
            endif
          endif
        enddo
      endif

C --- DLM single-site energies and distribution function-related quantities
C#ifdefC GF
C      if (ldlm /= 0 .and. nangl > 0 .and. .not. omonly) then
C        if (procid == master) then
C          call awrit0('%N Exchange constants J0 (mRy)',outs,80,stdo)
C          do  ib = 1, nbas
C            if (s_site(ib)%ncomp < 2) cycle
C            write(stdo,801) ib, (1d3*s_site(ib)%j0(j,1),j=1,s_site(ib)%ncomp)
C  801       format(' Site ',i3,':'/100(8F9.3/))
C          enddo
C        endif
C        call dlmmxy(s_ctrl%lbxc/4,nl,nspc,nclasp,nangl,
C     .    s_ctrl%dclabl,s_pot%mxy,s_pot%qnu,s_pot%thetcl,s_pot%shfac)
C        call dlmenrg(nclasp+1,nclasp+nangl,s_ham%etrms,efermi)
C        call dlmentrp(s_site,nbas,entropy)
C        if (procid == master) then
C          call awrit1('%N Entropy=%1,7;7d',outs,80,stdo,entropy)
C        endif
C        call pkgdlm(nbas,s_site,s_pot%gibbs,s_pot%shfac)
C        if (procid == master .and. bittst(s_ctrl%lbxc,4)) then
C          ifi = fopn('shfac')
C          rewind ifi
CC         call info0(20,0,0,' GFASA:  writing constraining B fields to shfac file ...')
CC         call iomagf(2,nclass+nangl,1,s_pot%shfac,s_pot%shfac,1,-ifi)
C          call ioextf(2,'constraining B fields to shfac',nclass+nangl,1,s_pot%shfac,s_pot%shfac,3,1,-ifi)
C          call fclr('shfac',ifi)
C        endif
C      endif
Cc     if (.not. lweiss) then
Cc       real(8), allocatable :: gbold(:)
Cc       allocate(gbold(nclasp+nangl))
Cc       call dcopy(nclspd,s_pot%gibbs,1,gbold,1)
Cc       call gibbswts(nclasp,s_ctrl%ncomp,maxth,tdlm,s_pot%thetcl,s_pot%dlmwt,
Cc    .    s_ham%etrms,s_pot%gibbs)
Cc       call vmix2(nclspd,s_ctrl%dclabl,0.2d0,gbold,s_pot%gibbs)
Cc       deallocate(gbold)
Cc     endif
C      if (lweiss .and. dabs(tdlm) > 1d-6) then
C        allocate(temp(ndlmcl))
C        call dlmtemp(0,nclasp,s_ctrl%ncomp,s_pot%thetcl,s_pot%gibbs,tdlm,s_ham%etrms,temp)
C        call dlmtemp(1,nclasp,s_ctrl%ncomp,s_pot%thetcl,s_pot%gibbs,tdlm,s_ham%etrms,temp)
C        deallocate(temp)
C      endif
C#endif

C --- Early program exit ---
      if (s_ctrl%quit == 4 .or. nevmx <= 0) then
        i = s_ctrl%nvario
        if (procid == master) then
          call nwitsv(1+2,i,'h67',nsp,amgm,etot)
        endif

C   ... Invoke shell command if supplied
        outs = ' '
        if (nit > 0 .and. cmdopt('--sh=',5,0,outs)) then
          call awrit0(' LM: invoke sh="'//outs(6:len(outs))//'%a"',' ',len(outs),stdo)
          call fsystm(outs(6:len(outs)),i)
          call awrit1(' LM: shell returned value %i',' ',80,stdo,i)
        endif

        call rx0(prgnam//': quit after bands')
      endif

C --- Magnetic torque ---
C#ifdefC NC
C      if (lmfrce) then
C        call dcopy(3*neul*nbasp,s_ham%eula,1,eulo,1)
C        allocate(frc(3*nbas+3))
C        call magtrq(nbasp,nl,nclasp,s_ctrl%ipc,sdmod,sdprm,s_site,
C     .    s_spec,s_pot%pp,s_pot%rhos,nrhos,ehf,s_ham%eula,neul,frc,
C     .    s_pot%aamom)
C        if (mod(sdmod,10) == 3) then
C          call rx('lmasa not ready for sdmod=3')
CCC     ... for now, copy eula to work array with one eula/atom
CC          stop 'nl->neula'
CC          real(8), allocatable :: wk2(:)
CC          allocate(wk2(2*nbas+3))
CC          call dpscop(s_ham%eula,wk2,nbasp,nbasp*nl+1,nbas+1,1d0)
CC          call dpscop(s_ham%eula,wk2,nbasp,1,1,1d0)
CC          real(8), allocatable :: wk(:)
CC          allocate(wk(3*nbas+3))
CC          print *, 'put irmmd in bsi,bswk'
CC          call u_bsi(s_sdyn,i,j,iv,xv,xx1,xx2,obswk,k,-1)
CCC     ... pull out second number from verbosity stack
CC          call togprt
CC          j = iprint()
CC          call poppr
CC          ifi = fopn('SAVE')
CC          call poseof(ifi)
CC          call mm_dyn(nbas,xsi,wk2,xx,wk,xx,xx,nstep,irmmd,xx,0,
CC     .      w(obswk),frc,amag,s_pot%aamom,ehf,s_sdyn,nvarms,ifi)
CC          call fclose(ifi)
CC          call pshpr(j)
CC          call togprt
CC          do  36  j = 1, nl
CC          call dpscop(wk2,s_ham%eula,nbasp,1,nbasp*(j-1)+1,1d0)
CC   36     call dpscop(wk2,s_ham%eula,nbasp,nbasp+1,nbasp*(nl+j-1)+1,
CC     .      1d0)
CC          deallocate(wk2,wk)
C        endif
C        deallocate(frc)
C      endif
C#endif

C --- Estimate self-consistent moments from dielectric response ---
      elind = s_ham%elind; i1 = str_pack('mix',-2,s_strn,strn)
      call pshpr(1)
      call parms0(i1,i2,xxp,1)
      xv = 0
      if (.not. parmxp(nit,strn,len(strn),xv,xv,xv,xv,xv,elind,xv,xv,
     .  outs,xv,xv,xv,xv)) call rx('LMASA: parse in parmxp failed')
      call parms0(i1,i2,xxp,-1)
      call poppr
      call asalsq(iscr,s_ctrl,s_lat,s_spec,s_ham,s_bz,s_pot,elind,vrl,qold,s_pot%qnu)

C --- Enforce charge neutrality ---
      if (cmdopt('--zerq',6,0,outs)) then
        call asazerq(outs(7:),2,s_ctrl,s_lat,s_spec,s_pot)
      endif

C ... Update changes in vconst
      s_pot%vconst = vconst

C --- Shift moments to the center of gravity of band ---
      makepp = getdig(sdmod,3,10) == 0
      if (makepp) then
C       Replace arrays with DLM-size arrays
        call mixpqc(nclspd,s_ctrl%nrc,nl,nsp,grp2,s_pot%pnu,s_pot%qnu)
        pmin = s_ham%pmin
        if (IAND(s_ctrl%lbas,16) == 0) then
          allocate(llmx(nclspd)); call icopy(nclspd,lmx,1,llmx,1)
C#ifdefC PGF
C          if (lpgf == 2) then
C            call ivset(llmx,1,nclspd,-1)
C            do  i = nbas+1, nbasp
C              j = s_ctrl%ipc(i)
C              llmx(j) = lmx(j)
C            enddo
C          endif
C#endif
          call shftpq(lrel,nclspd,s_ctrl%nrc,nsp,nl,llmx,s_ctrl%rmax,avw,
     .      s_pot%pp,s_pot%qnu,s_pot%qnur,idmod,.false.,pmin,s_pot%pnu,
     .      s_pot%qnu,s_pot%qnur,xv)
          deallocate(llmx)
        endif
      endif

C#ifdefC LMCNST
CC --- Make output ves ---
C      call asamad(s_ctrl,s_pot,s_lat,s_spec,
C     .  0,s_pot%pnu,s_pot%qnu,vrl,s_pot%ves,emad,trumad,vmtz,s_ham%etrms)
CC --- Map ves into vtil ---
C      real(8), allocatable :: vtil(:)
C      real(8), allocatable :: mad2(:)
C      real(8), allocatable :: umad(:)
C      real(8), allocatable :: evmad(:)
C      real(8), allocatable :: votil(:)
C      allocate(vtil(nbas))
C      allocate(votil(nbas))
C      allocate(mad2(nbas**2))
C      allocate(umad(nbas**2))
C      allocate(evmad(nbas))
C      call rotmad(nbas,nclasp,s_ctrl%ipc,wsr,mad,vold,s_pot%ves,
C     .  s_ctrl%nrc,emad0,modcst,mad2,s_ctrl%group,
C     .  evmad,umad,votil,vtil,nvmix)
CC --- Mix vtil; back transform to ves ---
C      call vmix(nclasp,nl,nsp,nbas,s_pot%ves,vold,s_pot%pp,
C     .  w(oiclas),s_ctrl%nrc,vtil,votil,umad,evmad,
C     .  modcst,s_ctrl%group,wsr,mad,nvmix,switch(4))
C      deallocate(vtil)
C      deallocate(votil)
C      deallocate(mad2)
C      deallocate(umad)
C      deallocate(evmad)
C#endif

C --- Mix moments and/or Euler angles ---
      nwmoms = .false.
      if (makepp) then
        dmxprm(9) = betv
C   ... Get mixing string; set up mixing constraint for pqmix
        i1 = str_pack('mix',-2,s_strn,strn)
C   ... Save qnu for independent potential mixing
        i = s_mix%lxpot
        pnus => s_pot%pnu
        qnus => s_pot%qnu
        if (i == 2) then
          if (lrel == 2) call rx('lrel=2 and indep mixing')
          allocate(pnus(nl*nsp,max(nclsppd,nspec)),qnus(3*nl*nsp,max(nclsppd,nspec)))
          call dcopy(nlspcp,s_pot%pnu,1,pnus,1)
          call dcopy(3*nlspcp,s_pot%qnu,1,qnus,1)
        endif
C   ... Extra degrees of freedom
        allocate(mixcst(0:nclasp+nangl))
        call iinit(mixcst,(nclasp+1+nangl))
        i = 3*nl*nl*nbas
        if (lnsph /= 0) i = i + nqpp*4*nsp*nbas
        allocate(xold(i))
        allocate(xnew(i))
        sw = lmfrce .and. getdig(sdmod,1,10) == 0 .and. neul > 0
        call pvpqm1(0,s_spec,nclass,nclspd,nsp,nbas,lpgf,F,lnsph /= 0,
     .    sw,neul,s_ctrl%nrc,s_ctrl%ics,s_ctrl%ipc,mixcst,nqpp,
     .    qppo,s_pot%qpp,xx,xx,eulo,s_ham%eula,vcnsts,vconst,i,
     .    xold,xnew)
        dmxprm(34) = dmxprm(4)
        if (.not. omonly) then
C  ...    Dummy DLM sites should not be mixed (would screw up Broyden)
C  ...    Add nthet to mixcst to cancel mixing of dummy DLM sites
C  ...    (s_ctrl%ncomp = 0 for normal sites)
C#ifdefC GF
C          if (ldlm /= 0) call iaxpy(nclasp,1,s_ctrl%ncomp,1,mixcst(1),1)
C#endif
          call pqmix(nclspd,nl,lmx,nsp,i,nit,strn,dmxprm,mixcst,
     .      pold,qold,qrold,xold,s_pot%pnu,s_pot%qnu,s_pot%qnur,xnew)
C  ...    Subtract nthet back from omxcst
          if (ldlm > 0) call iaxpy(nclasp,-1,s_ctrl%ncomp,1,mixcst(1),1)
        endif
        call pvpqm1(1,s_spec,nclass,nclasp,nsp,nbas,lpgf,F,lnsph /= 0,
     .    sw,neul,s_ctrl%nrc,s_ctrl%ics,s_ctrl%ipc,mixcst,nqpp,
     .    qppo,s_pot%qpp,xx,xx,eulo,s_ham%eula,vcnsts,vconst,i,
     .    xold,xnew)
        deallocate(mixcst,xold,xnew)

        betv = dmxprm(9)
C   ... Independent potential mixing
        if (lves == 2) then
C     ... ves [mixed rho or output rho], depending on mix->lxpot
C         NB: Not needed if screened ves available.
C         But this makes ves[screened output rho, so ok.
C         if (mod(iscr/2,2) == 0) then
          call pshpr(iprint()-20)
          allocate(wk(22*nclspd))
          call asamad(s_ctrl,s_pot,s_lat,s_spec,
     .      0,pnus,qnus,vrl,s_pot%ves,emad,trumad,vmtz,wk)
          deallocate(wk)
          call poppr
C         endif
          call vmix2(nclspd,s_ctrl%dclabl,betv,vold,s_pot%ves)
          if (s_mix%lxpot == 2) deallocate(pnus,qnus)
        endif
        if (procid == master) then
          call shoctl(s_ctrl,s_spec,s_pot,1,fopn('LOG'))
        endif
        if (mod(nint(dmxprm(25)),10) == 1) nwmoms = .true.

C     Potential is frozen: keep starting P,Q, d.c. terms
      else
C       call shoctl(s_ctrl,s_spec,s_pot,0,stdo)
        call dcopy(nlspcp,pold,1,s_pot%pnu,1)
        call dcopy(3*nlspcp,qold,1,s_pot%qnu,1)
        if (lrel == 2) call dcopy(4*nlspcr,qrold,1,s_pot%qnur,1)
C       call shoctl(s_ctrl,s_spec,s_pot,0,stdo)
        s_ham%eterms = etrmss
      endif

C ... Save mixed qpp array
      if (lnsph /= 0) call ioqpp(T,s_ctrl,s_pot)

C ... Independent mixing of Euler angles (moms into scratch)
      if (lmfrce .and. getdig(sdmod,1,10) == 1 .and. neul > 0) then
C       if (i2 <= i1) then
C         outs = 'A0,b=1,w=0,0,wa=1,fn=ma'
C       else
C         outs = sstrn(i1:i2)
C       endif
        i1 = str_pack('amix',-2,s_strn,outs)
        if (len_trim(outs) == 0) outs = 'A0,b=1,w=0,0,wa=1,fn=ma'
        call info0(31,1,0,' Independently mix Euler angles')
        if (lrel == 2) call rx('Incompatible with lrel=2')
        allocate(wk(nlspcp))
        allocate(wk2(3*nlspcp))
        i1 = size(s_pot%qnur)
        allocate(wkr(i1))
        call dcopy(nlspcp,s_pot%pnu,1,wk,1)
        call dcopy(3*nlspcp,s_pot%qnu,1,wk2,1)
        call dcopy(i1,s_pot%qnur,1,wkr,1)
C   ... Pull out regular mixing variables, poke in Euler-specific
        call parms0(dmxprm(27),dmxprm(28),dmxprm(29),1)
        call parms0(dmxprm(30),dmxprm(31),dmxprm(32),-1)

        allocate(mixcst(0:nclasp))
        call iinit(mixcst,(nclasp+1))
        allocate(xold(3*nl*nl*nbas))
        allocate(xnew(3*nl*nl*nbas))
        call pvpqm1(0,s_spec,nclass,nclasp,nsp,nbas,0,F,F,T,neul,
     .    s_ctrl%nrc,s_ctrl%ics,s_ctrl%ipc,mixcst,nqpp,xx,xx,xx,xx,
     .    eulo,s_ham%eula,vcnsts,vconst,i,xold,xnew)
        dmxeu(4) = 0
        dmxeu(5) = 0
        dmxeu(34) = 1
        call pqmix(nclasp,nl,lmx,nsp,i,nit,outs,dmxeu,0,
     .    s_pot%pnu,s_pot%qnu,s_pot%qnur,xold,wk,wk2,wkr,xnew)
        call pvpqm1(1,s_spec,nclass,nclasp,nsp,nbas,0,F,F,T,neul,
     .    s_ctrl%nrc,s_ctrl%ics,s_ctrl%ipc,mixcst,nqpp,xx,xx,xx,xx,
     .    eulo,s_ham%eula,vcnsts,vconst,i,xold,xnew)
        call parms0(dmxprm(30),dmxprm(31),dmxprm(32),1)
        call parms0(dmxprm(27),dmxprm(28),dmxprm(29),-1)
        deallocate(mixcst,wk,wk2,wkr,xold,xnew)
      endif

      if (s_ctrl%quit == 32) then
        i = s_ctrl%nvario
        if (procid == master) then
          call nwitsv(1+2,i,'h67',nsp,amgm,etot)
        endif
        call rx0(prgnam//'%a:  quit = rho')
      endif

C ... Write updated spin quantization axis and Euler angles to disk
      if (bittst(lncol,1) .and. getdig(sdmod,0,10) /= 3) then
        if (procid == master) then
C          ifi = fopn('bxc')
C          rewind ifi
C          call iomagf(1,nclass,1,s_pot%bxc,s_pot%bxc,1,-ifi)
C          call fclose(ifi)
          if (getdig(sdmod,0,10) /= 3) then
            ifi = fopn('EULA')
            rewind ifi
            call ioeula(nbasp,nl,s_ham%eula,neul,0d0,-ifi)
            call fclose(ifi)
          endif
        endif
      endif

C ... Skip over self-consistency when Euler angles are l-dependent
      if (.not. nwmoms) goto 58
      goto 10

   99 continue
C     call tcx('lmasa')

      end
