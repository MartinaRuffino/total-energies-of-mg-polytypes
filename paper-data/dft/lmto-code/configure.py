#!/usr/bin/env python

import os, sys, re
import subprocess
from nj.nj import loadfile, writefile, laddprefix, laddsuffix, lsandwich, nj_t


nj = nj_t()


ofln = list(filter(lambda o: not o.startswith('-'), sys.argv[1:]))
ofln = 'build.ninja' if len(ofln) == 0 else ofln[0]
nj.bld(['$srcpath/configure.py'], 'nj', ['$bldpath/'+ofln], si='flags.mk $srcpath/nj/nj.py $srcpath/testable.md'.split())
writefile(ofln, nj.rules)

#gwbuild = os.path.dirname(sys.argv[0]) + '/gw/configure.py'
#gwavail = os.path.isfile(gwbuild)
#if gwavail:
    #subprocess.call((gwbuild + ' build-gw.ninja -subninja').split())
    #nj.direct += '\nsubninja build-gw.ninja'

cflags  = nj.flags['cflags']
fflags  = nj.flags['fflags']
ppflags = nj.flags['ppflags']

if re.search(r'-dMPI\b', ppflags) != None:
    print('Atom based MPI parallelism is disabled untill improved version arrives! Replacing -dMPI with -dMPIK ...')
    ppflags.replace('-dMPI','-dMPIK')

mpiavail = '-dMPIK' in ppflags
ppflags = ppflags.replace('-dMPIK', '') # do not preprocess for mpi, the flags is only used to decide whether to use nullmpi.
cudavail = '-dCUDA' in ppflags
scaavail = '-dSCALAPACK' in ppflags
ppflags = ppflags.replace('-dSCALAPACK', '') # the flag is only used to decide whether to use nullsca
elpavail = '-dELPA' in ppflags
lsxavail = '-dLS++' in ppflags
ncavail  = '-dNC' in ppflags
mklavail = '-DMKL' in nj.flags['gwflags'] or '-dMKL' in ppflags
h5avail  = '-DUSE_HDF5' in nj.flags['gwflags']
makedso  = '-fpic' in fflags.lower() or '-fpie' in fflags.lower()
if makedso: nj.libsfx = '.so'

if scaavail and not mpiavail:
   sys.stderr.write('Disabling ScaLAPACK & friends (missing MPI)\n')
   ppflags = ppflags.replace('-dSCALAPACK', '')
   scaavail = False

if elpavail and not (scaavail and mpiavail):
   sys.stderr.write('Disabling ELPA (missing MPI or ScaLAPACK)\n')
   ppflags = ppflags.replace('-dELPA', '')
   elpavail = False

if cudavail and (env['NVCC'] == '' or env['CXX'] == ''):
   sys.stderr.write('Disabling CUDA (missing NVCC or CXX)\n')
   ppflags = re.sub(r'-dCUDA\w*', '', ppflags)

#print 'mpiavail:', mpiavail
#print 'cudavail:', cudavail
#print 'scaavail:', scaavail
#print 'elpavail:', elpavail


nj.emods.update(set('xc_f90_lib_m mpi mpi_f08 elpa1 elpa2'.split()))
if not mklavail: nj.emods.update(set(['mkl_service']))
nj.gincs.update(set(['gitref.h']))
nj.eincs.update(set(['hdf5.h']))

if not mpiavail: nj.emods -= set(['mpi'])

nj.bld(['$srcpath/utils/ccomp.c'],'cl',['$bldpath/ccomp'],vs=dict(cc=nj.flags['knrcc'], cflags = '-Dunix', ldflags= '', implicit_cflags='', implicit_ldflags='', libs=''))


#git = r'git --no-pager -C $srcpath' # the -C option is not available on archaic systems like rhel6/7
git = r'git --no-pager --work-tree=$srcpath --git-dir=$srcpath/.git' # --work-tree is enough on rhel7; 6 needs --git-dir too
tag = r'{git} tag --points-at HEAD'.format(git=git) # multiple tags will confuse this, HEAD added for compatibility with rhel7's git version (1.8.3.1). rhel6's git (1.7.1) does not support --points-at so only references are used there.
cmt = r'{git} log -n 1 --format=%h'.format(git=git)
vref_fl = '$bldpath/gitref.h'

nj.bld([],'phony',['refcheck']) # ensures gitref is always out of date.
nj.bld(['refcheck'],'cmd',[vref_fl],vs=dict(
    # this mess is necessary to avoid changing the timestamp of gitref.h when the content will not change. it may be neutralised by specifying it as order only dependency too but then an update will not trigger a rebuild resulting in misleading version being printed from executables.
    # compatible with /bin/sh -> dash
    command = r'printf "git version 1.7.10\n$$(git --version)" | sort -VC && vref="$$({tag})"; vref="      gitref = \"$${{vref:-"$$({cmt})"}}\""; [ -r {vref_fl} ] && [ "$$vref" = "$$(cat {vref_fl})" ] || echo "$$vref" > {vref_fl}'.format(tag=tag, cmt=cmt, vref_fl=vref_fl),
    restat = 1,
    description = "refcheck"
))


sl_src = laddprefix('slatsm/','''
    fmain.F90 fcuda.F90 fsubs.f90 mod_ctx.f90 vbdist.f90 darrayops.f90 posix.f90
    a2bin.f  initldev.f90 unievp.f90 elpaglue.F90 idalloc.f omppid.f
    a2bina.f  a2rotm.f  a2vec.f  afcmd.f  alloc.f  amix.f
    arrprt.f  awrite.f  bin2a.f  bitops.f  brentm.f  brmin.f
    broydn.f  broyj.f  cbabk2.f  cbal.f  cdasin.f  cdiv.f
    cg.f  chcase.f  chkhss.f  cksumf.f  cmdopt.f  cmdstr.f cmpi123.f
    comqr.f  comqr2.f  convolve.f  corth.f  cplxdm.f  cpstr.f
    cpudel.f  cpvprm.f  crdbf.f  cross.f  csmaln.f  cspline.f comfac.f
    csroot.f  d1mach.f  da2spr.f  dampy.f  daxpby.f  dbesj.f
    dbesnu.f  dbsstp.f  dcabsr.f  dcsqrt.f  decrnd.f
    derfc.f  dfdump.f  dfftqp.f  dgamma.f  dgedi.f  dgefa.f dfphi.f
    dgesl.f  diagno.f  dibesl.f  dinv.f  dinv33.f  discop.f
    djbesl.f  dkbesl.f  dlnref.f  dmadd.f  dmcpy.f  dmpy.f
    dmpy22.f  dmpyt.f  dmsadd.f  dmscop.f  dmx22.f  dot3.f
    dotprd.f  dpack.f  dpadd.f  dpcopy.f  dpdbyl.f  dpdot.f
    dpdump.f  dpmpy.f  dpsadd.f  dpscop.f  dpsmpy.f  dqinv.f
    dqinvb.f  dqsmpy.f  dsbak.f  dschd.f  dsev1.f  dshell.f
    dsidi.f  dsifa.f  dsisl.f  dsmpy.f  dspbmm.f
    dspdi.f  dspfa.f  dsprmv.f  dspsl.f  dsred.f  dsredx.f dstrbpx.f
    dsum.f  dsumdf.f  dtribx.f  dtridx.f  dupack.f  dvheap.f
    dvshel.f  dybesl.f  eigch.f  epslon.f  errmsg.f  extstr.f
    falsi.f  fftz3.f  fftzv.f  fgtcat.f  finits.f  fopna.f
    fpiint.f  garg.f  getcat.f  gettime.f  gradfl.f  gradzr.f
    gvgetf.f  hbak.f  hchd.f  hcinv.f  hred.f  hrmint.f
    htrib3.f  htribk.f  htribx.f  htrid3.f  htridi.f  htridx.f
    hunti.f  huntst.f  huntx.f  i1mach.f  iasum.f  iaxpy.f
    ichksm.f  idmxmn.f  idnear.f  iiamax.f  iinear.f
    imtql2.f  imtqlv.f  imxmn.f  info.f  ipdump.f  iprint.f
    ipsw.f  isanrg.f  ishell.f  isum.f  isw.f  iswap.f
    ivheap.f ivshel.f  iyamax.f  ll.f  lstshr.f linemin.f logwarn.f
    macro.f  minfit.f  mkdlst.f  mkilst.f  mklagw.f mklegw.f makdkron.f mpi_append_files.f
    mpibc.f  mrqstp.f  mstmin.f  mtch3f.f  mxint.f  nada.f
    nglob.f  nintsp.f  nxtarg.f  nxtgrp.f  packi.f  packv.f
    pade.f  parchv.f  parrng.f  pars1v.f  parsvc.f  pertev.f
    platsl.f  pmblk.f  polcof.f  poldif.f  poldvm.f  polint.f
    politg.f  poseof.f  ppfac.f  pralloc.f  praxis.f  pretty.f
    pythag.f  qdiv.f  query.f  r1mach.f  ran1.f
    rdfiln.f  rdm.f  rdm2.f  rdtok.f  rfalsi.f  ropcsm.f
    ropylg.f  ropyln.f  ropyln1.f  roundp.f  rs.f  ry2ev.f s8tor8.f
    scglp1.f  setfac.f  spchar.f  stkops.f  strings.f  strip.f
    strxls.f  svbksb.f  svd.f  symvar.f  symvec.f  syscall.f
    tc.f  tinvit.f  tmix.f  tql2.f  tqlrat.f  tred1.f
    tred2.f  verfc.f  vmem.f  wordg.f  words.f  ygemm.f
    ymcpy.f  ymsadd.f  ymscop.f  ymtrns.f  yqinv.f  yqinvb.f
    ysbnv.f  ysbnvb.f  yscal.f  yspbmm.f  yspr2a.f  ywrm.f
    yyaxpy.f  yydot.f  yydotc.f  yygefa.f  yygemm.f  yygesl.f
    yyhbak.f  yyhchd.f  yyhifa.f  yyhmpy.f  yyhrdx.f  yyhred.f
    yympy.f  yympyt.f  yyqinv.f  yyqnvb.f  yysbib.f  yysbmm.f
    yysbnv.f  yysp2a.f  yyxcpy.f  z1btf2.f  z1btrf.f  z1btrs.f
    z1tbsv.f  z2herm.f  zampy.f  zerfc.f  zfftqp.f  zgeevs.f
    zgefa.f  zhdmpy.f  zhev.f  zhisl.f  zhmpy.f  zmcpy.f zinv.f
    zmpy.f  zmpy22.f  zmsadd.f  zmscop.f  zpack.f  zpmblk.f
    zqinv.f  zqinvb.f  zqsmpy.f  zsum.f  ztoy.f  ztoyy.f
    dlength.f icopy.f  csubs.c fmagma.f90 ytab0.f init.f det33.f
    rotspa.f ymscpk.f zdia22.f meshes.f90 testcov.f90 mod_utils.f90
''')
# rotspv.f
# dci.f dsi.f lgand.f
#mklagw.f pzhev.f

if not mpiavail: sl_src.append('slatsm/nullmpi.f90')
if not scaavail: sl_src.append('slatsm/nullsca.f90')
if not cudavail:
   sl_src.extend(laddprefix('slatsm/', 'nullcuda.c nullmagma.f90'))
else:
   sl_src.extend(laddprefix('slatsm/', 'cuda_auxiliary.cu devremapshim.c'))
   if mpiavail:
      sl_src.append('slatsm/pcudgemm.cxx')

nj.mks(sl_src,vs=dict(ppflags = ppflags.replace('-dMPIK','-dMPI')))

nj.fppstyle = 'cpp'
h5_vs = dict()
if mpiavail: h5_vs['cflags'] = cflags + ' -DUSE_MPI'
if mpiavail: h5_vs['fflags'] = fflags + ' -DUSE_MPI'
if h5avail: nj.mks(laddprefix('slatsm/','h5.F90 h5.c'), vs=h5_vs)
nj.fppstyle = 'ccomp'

if mklavail:
    incs = re.findall(r'(?x)-I \s* (.*?) (?:\s+|$)',nj.flags['fflags'])
    mkl_service = [(inc+'/mkl_service.f90') for inc in incs if os.path.exists(inc+'/mkl_service.f90')]
    if len(mkl_service) == 0:
        sys.stderr.write('\nError: could not find mkl_service.f90 in the include paths, maybe the environment changed since generating the flags.mk file? Possible solutions:\n    (1) Restore the original environment:\n      a) reload the appropriate modules providing the MKL\n      b) source env scripts which may come with the intel compiler\n    (2) Regenerate flags.mk with genflags.py in the present environment\n')
        sys.exit(1)
    nj.bld(mkl_service[0:1],'fc',['$bldpath/slatsm/mkl_service.f90.o'],oi=['$bldpath/$modpath/mkl_service.mod'])

in_src = laddprefix('v7input/','rdctrl.f rdctrlchk.f m_rdctrl.f m_rdctrlchk.f m_gtv.f m_toksw.f')

nj.mks(in_src)

lm_src = laddprefix('subs/','''
    pointers.f
    acoord.f addbas.f addtog.f addtos.f afmevc.f aiocls.f aiocor.f
    aiogen.f aiomom.f aiomp.f aiopar.f aiopot.f aiorme.f aiosop.f aiova.f
    aioxtn.f alwats.f amagnc.f angtab.f asadc.f asaddq.f asadmu.f asados.f
    asalsq.f asamad.f asaopm.f asaqmp.f asars.f asastr.f asavqm.f asetot.f
    asldau.f asvsph.f asymop.f atfold.f atomsr.f atsev.f
    avwsr.f baslst.f bdotsr.f besslr.f bloch.f blochi.f bndasa.f bsortl.f
    bsrhos.f bzints.f bzjnos.f bzmesh.f bzmio.f bzwts.f bzwtsf.f cb2sph.f
    ccutup.f chgmsh.f chkdmu.f chksg.f chkstr.f clcoin.f clabel.f clist.f clsprm.f
    clsset.f cmplat.f config.f contet.f cpaidx.f cpaz.f cshear.f decmp2.f deflx.f defpq.f
    defspc.f defwsr.f delstd.f delstp.f dlmidx.f dlmn.f dlmq.f
    dmatio.f dosio.f dosmsh.f dosspl.f dostet.f drr2.f dsolv2.f dstrbp.f
    dvdos.f e0parm.f easypbe.f efrng2.f ehcor.f enutod.f evxc.f evxcp.f epsdipole.f
    evxcv.f fadflb.f fermi.f findes.f findqbz.f fixef0.f fixpos.f fixplat.f
    fovlp.f freeat.f fswgts.f gengrp.f gensym.f getbzp.f getcls.f getcor.f
    getef.f getidu.f getjdosw.f getmom.f getsclo.f getq.f getqp.f getzv.f ggiml.f
    gh1c.f ghdiml.f ghgiml.f ghiml.f gklml.f gprdsr.f groupg.f grpfnd.f
    grpgen.f grpop.f gskham.f gtbsl1.f gtbvec.f gtpcor.f gtpdos.f gvctof.f
    gvlist.f gvlst2.f gvmtch.f gxpml.f h1c.f h2grad.f hamfb3.f hanb2s.f hanr.f
    hansmd.f hansmr.f hansmz.f hansr4.f hansr5.f hansr.f hansrz.f hanszd.f
    hansz.f hblstr.f hcr2a.f hdxpml.f hft2rs.f
    hml2nc.f hmladb.f hmlt2c.f hmltnc.f hmltns.f hmltso.f hnsmft.f hoffs.f
    hpairm.f hs1cx.f hsmml.f hsmq.f hstrux.f hxpml.f i3cntr.f iax2fd.f
    ibloch.f iclbas.f intdrv.f intnos.f invalp.f invbl.f iobsmv.f iodos.f
    iodosw.f ioeula.f iofa.f ioka.f iomagf.f ioextf.f iomomq.f iomoms.f iomomx.f
    ioorbp.f iopos.f ioqpp.f ioatmx.f iorbtm.f iores.f iors.f iorsj.f iosave.f ioseh.f iosg.f
    iosig.f iosite.f iositp.f iosits.f iostr.f iostrt.f iosv.f iovext.f iovshf.f
    ioxbs.f ioxsf.f kink2g.f latlim.f lattc.f lattdf.f lattic.f latvec.f
    ldau.f lfrog.f lgen.f lmasa.f lmaux.f lmfitmr.f lmgdos.f lmrefdos.f
    lx2vec.f maadot.f madmat.f madpot.f magtrq.f makalp.f makbet.f
    makbzm.f makdia.f makdla.f makdos.f makdsd.f makidx.f makipt.f maknos.f
    maknrs.f makpph.f makptf.f makrcz.f makrm0.f makrot.f makrwf.f maks0.f
    makwts.f mcasim.f mchan.f mdinit.f mdout.f meptr.f mixpqb.f mixpqc.f
    mixpq.f mkalpv.f mkbdot.f mkbfld.f mkcond.f mkiaxd.f mkjdos.f mkplat.f
    mkppar.f mkqp.f mkrtab.f mkrwat.f mksbet.f mksodn.f mksoph.f mks_sym.f
    mksym.f mksyml.f mktra2.f mktral.f moment.f mpint.f mshsiz.f mstrx2.f mstrx3.f mstokm.f lmorder.f
    mtchae.f mullmf.f mxbnds.f mxmymz.f newalp.f nghbor.f nglob.f nmefac.f
    nm.f nmham.f nmpot.f nnrl.f norm2g.f nosbzm.f npr2tb.f nwit.f ogvlst.f
    optinq.f optint.f oraddp.f orbl.f ovlchk.f ovmin.f pa2plm.f pairc.f
    pairg.f pairp6.f parmxp.f pgbasp.f pgfp.f pgfset.f phidx.f phmbls.f
    phvsfp.f pkli.f plana.f plm2pa.f popted.f pp2enu.f pp2hvs.f pptrns.f
    pqedit.f pqmix.f prdmts.f prjpos.f prjrsp.f prmsk2.f prmx.f projql.f
    prrmsh.f prtbas.f prterr.f psymop.f pvpqm1.f qdist.f qplin.f radgkg.f radgkl.f
    radgkv.f radkj.f radmsh.f radmwt.f radpkl.f rdeq.f rdsigm.f readrssw.f relax.f rdeq_old.f
    remhyb.f rgrme.f rhocor.f rlxstp.f rmesh.f rmeshprm.f rnatm.f ropbes.f roplat.f rotevs.f
    roth.f rothrm.f rotmad.f rotpnt.f rotspn.f rotspu.f rotspv.f rotwf.f rotycs.f
    rsedita.f rseq.f rsmsym.f rsort.f rsqnri.f rstr0.f s2oka.f s2sph.f
    saiwav.f salias.f salph1.f salph2.f salphg.f savvxc.f sblham.f scalpv.f
    scalsr.f scg.f sclwsr.f scrmom.f secm2c.f secmat.f secmtn.f setcg.f
    setnorb.f sfill.f shear.f shffle.f shftpp.f shftpq.f shoblk.f shoctl.f
    shoist.f shopol.f shorbz.f shorigv.f shorpsa.f shorps.f shortn.f90 shosg.f shoshl.f
    shostr.f shosym.f siged.f siteid.f skhso.f slinz.f sogrp.f slist.f soldhj.f sole0g.f solgsg.f
    solgsv.f solhjg.f solhsg.f soprm.f spbndw.f spcgrp.f spinav.f
    splcls.f splwts.f sstrxq.f stoner.f strck.f strdif.f streqv.f strg1c.f
    strgck.f strg.f strrs.f strscr.f strx00.f strxops.f strxq.f strxsu.f
    subasi.f subzi.f suclst.f suctrl.f sudmtu.f suemph.f sugcut.f sugvec.f sudossw.f
    suham2.f suham.f suidx.f suldau.f sumlst.f sumsro.f supcel.f stackcel.f supot.f
    suqlst.f susite.f swapdv.f90 sylm.f sylmnc.f symcry.f symdmu.f symfor.f
    symiax.f symlat.f sympad.f symprj.f symq.f symstr.f symtab.f symtbl.f
    tailsm.f tbaddh.f tetfbz.f tetirr.f tetwtq.f trsdot.f trs.f trysop.f upd_nthet.f
    uspecb.f v0intr.f vecpkl.f verlet.f vldau.f vmix2.f vmix.f veecomp.f
    volsph.f vrbasp.f vsh2es.f vsl0.f vsl.f vtoq.f vxc0sp.f vxcgga.f vxcgr2.f vxnloc.f
    wkdbg.f wronjjv.f xc_libxc.f xlgen.f xxxbnd.f xxxdif.f ylmrtg.f yyhmul.f zercmv.f zhblock.f zslabl.f
    structures.f
''')
# partks.f

nj.mks(lm_src)

fp_src = laddprefix('fp/','''
    addbkg.f addrbl.f atqval.f augmat.f augmbl.f avsoef.f bndfp.f bstrux.f chimedit.f
    corprm.f dfaugm.f dfqkkl.f dfratm.f dfrce.f diagcv.f dpdbyl.f dpdftr.f
    efldos.f ephed.f elocp.f evalqp.f fklbl.f flocbl.f fpchk.f fradhd.f fradpk.f fsmbl.f
    fsmbpw.f ftlxp.f gaugm.f gfigbl.f ggugbl.f ghibl.f ghigbl.f ghios.f
    gklbl.f gklft.f gklq.f gradv1.f grfmsh.f grp2av.f gtbsl2.f gvsym.f hambl.f hambls.f hambldsdk.f
    hgugbl.f hhibl.f hhigbl.f hhugbl.f hklbl.f hklft.f hklgbl.f hklos.f
    hsibl.f hsibq.f hsmbl.f hsubblock.f hxpbl.f hxpgbl.f hxpos.f
    ioden.f iors.f jxpos.f lgstar.f lindsc.f lmfopb.f lmfp.f locpot.f locpt2.f
    loctsh.f makusp.f makusp2.f makusq.f mixrho.f mkdmtu.f mkehkf.f mkekin.f mkewgt.f
    mkorbm.f mkpdos.f mkpot.f mkrout.f modsig.f momusl.f msh21c.f mshdot.f mshint.f
    mshn3p.f mshvmt.f ncutcorrect.f ovlocr.f ovlpfa.f paug1.f paugb.f pkl2ro.f
    pnunew.f poinsp.f potpus.f praugm.f prrhat.f prtev.f rdovfa.f readevec.f rhgcmp.f
    rhogkl.f rhomom.f rhopos.f rlocbl.f rsedit.f rsibl.f rxes.f setofl.f sgvsym.f
    shorho.f smcorm.f smhsbl.f smshft.f smves.f smvxt.f smvxte.f smvxcm.f sopert.f
    splrho.f stonerpb.f stonerrsa.f sugws.f suphas.f supot.f surho.f suylg.f symrho.f
    tbhsi.f totfrc.f ugcomp.f vcdmel.f vesft.f vesgcm.f vlm2us.f vxcnlm.f vxcnls.f
    vxcnsp.f vxnlcc.f vxtrap.f wrhomt.f
    atwf.f addrwf.f
    hhiml.f hklml.f
''')

nj.mks(fp_src) #,vs=dict(ppflags = ppflags + ' -dLMF')

op_src = ['optics/asaopm.f','subs/bndasa.f'] + laddprefix('optics/','dosmsh.f fpopm.f gradme.f intdrv.f lindhard.f')  \
    + ['subs/lmasa.f'] + laddprefix('optics/','''ocnock.f optdme.f ooptdme.f optinq.f optint.f nesjdos.f optsf.f
    optshd.f optshs.f shdint.f shsint.f sliny.f tetsubs.f tetwtq.f tetwtt.f tetwtz.f veltet.f''')

nj.mks(op_src)


tb_src = laddprefix('tb/', '''mkstrx.F90 tbprl.f90 bldham.f90 tbeseld.f90 driver.f90
    secmtb.f90 bndtb.f90 tbfrce.f90 ztbloch.f90 mkstrxidx.f90 ztblochd.f90 mk2didc.f90
    findmolidc.f90 bndtb_alt.f90 rho2frce.f90 sockets.c fsockets.f90

    addes.f addoff.f90 addofff.f bangs.f bndtbf.f cffor.f chkme.f derivp.f
    dhmix.f dskham.f efgrad.f getujh.f hstr.f hstr0.f hstra.f invtb1.f
    iodell.f iomv.f ioqm.f iowts.f itrpck.f lumo.f makcg9.f makcgn.f
    makv0.f makvme.f makvpp.f mkrho0.f  mkwttb.f mrqmintb.f
    nlme.f pcut45.f qmix.f qmix2.f rcnsl0.f rcut.f rdtbh.f rhomix.f
    shotbm.f skham.f sulmxl.f swapv.f symre.f symrtb.f
    symtbd.f tbadh1.f tbadh2.f tbdiag.f  tbfdia.f tbfitmrq.f tbfrc2.f
    tbfrc3.f tbloch.f tbmpol.f tbshfl.f tbtote.f tbxbs.f  tbham.f90
    tbzint.f90 tiozll.f ume.f vmeder.f vppder.f tbesel.f secmtbf.f
    tbdos.f90
    ''')

if cudavail:
   tb_src.extend(laddprefix('tb/', 'mkstrxd_fcu.f90 mkstrxd_cu.cu sylm_c.cu h2rho.cu'))
else:
   tb_src.append('tb/mkstrxd.f90')

if not lsxavail: tb_src.append('tb/nulllsl++.f90')

nj.mks(tb_src,vs=dict(ppflags = ppflags.replace('-dMPIK','-dMPI')))



dt_src = laddprefix('dmft/', '''
    diagdenmat.f cmp_denmat.f cmp_overlap.f cmp_overlap_inv.f cmp_overlaptest.f cmpvalcharg_matsub4.f
    chk_selfcons.f chk_sigbar.f diagoham.f dmftevec.f
    embed_sigma.f getVnew.f prtchargneut.f iodsig.f iosegw.f iodmftdc.f
    dmftprojection.f makegloc2.f lmlproj.f makeeimp2.f makesiggwloc.f
    makeproj.f print_dmftu.f print_wtkb.f
    renorm_proj.f renorm_proj_cholesky.f sudmft.f readindmfl.f readsiginp.f zwrmdmft.f symdmftsig.f
    compressproj.f addsigmatog.f trivialsigqp.f makedc.f
    makehbar2.f mksigpade.f mkdmftstates.f
    makedelta3.f iodmftu.f print_sigf0.f print_gloc.f sigp2sigij.f testdenmat.f read_sigf0.f makebarechi.f printg.f vertex.f90
    ''')
nj.mks(dt_src) #, vs=dict(ppflags = ppflags + ' -dDMFT')
nj.mk('fp/lmfp.f',o='dmft/lmfp.f.o',vs=dict(ppflags = ppflags + ' -dDMFT'))


sx_src = laddprefix('sx/', '''
     asasx.f asxdev.f asxp0.f asxsgm.f asxsig.f desxdn.f desxdv.f dvsxl.f ibloch.f
     madmtq.f vbare.f wdiff.f wq0.f wsloc.f wsloca.f wstat.f wtrans.f
     ''')
nj.mks(sx_src) # , vs=dict(ppflags = ppflags + ' -dSX -dNC')
nj.mk('subs/lmasa.f',o='sx/lmasa.f.o',vs=dict(ppflags = ppflags + ' -dSX -dNC'))


nc_src = laddprefix('nc/', '''
 amagnc.f bdotsr.f bsrhos.f ham2nc.f hamss2.f hmfr2c.f hmfr3c.f
 hml2nc.f hmladb.f hmltnc.f hmltso.f magtrq.f mkbdia.f
 mksodn.f mksoph.f mmag.f mmdyn.f mmpair.f rotheu.f rotspn.f yhmlso.f
 rothnc.f
''')
nj.mks(nc_src) # , vs=dict(ppflags = ppflags + ' -dSX -dNC')
nj.mk('subs/bndasa.f', o='nc/bndasa.f.o', vs = dict(ppflags = ppflags + ' -dSX -dNC'))
nj.mk('subs/lmasa.f' , o='nc/lmasa.f.o' , vs = dict(ppflags = ppflags + ' -dSX -dNC'))
nj.mk('subs/lmaux.f' , o='nc/lmaux.f.o' , vs = dict(ppflags = ppflags + ' -dSX -dNC'))
nj.mk('subs/secm2c.f', o='nc/secm2c.f.o', vs = dict(ppflags = ppflags + ' -dSX -dNC'))
nj.mk('subs/secmat.f', o='nc/secmat.f.o', vs = dict(ppflags = ppflags + ' -dSX -dNC'))
nj.mk('subs/asaddq.f', o='nc/asaddq.f.o', vs = dict(ppflags = ppflags + ' -dSX -dNC'))


gf_src = laddprefix('gf/', '''
     gcpa.f sclmat.f adrotp.f advshp.f asagqr.f asajft.f besslz.f cpadlm.f cpagamma.f
     dfphiz.f dglrft.f dlmenrg.f dlmentrp.f dlmgrp2.f dlminit.f dlmmxy.f dlmsumev.f dlmtemp.f
     dlmwgts.f emesh.f exasa.f getpfi.f gf1kp.f gfasa.f gfasdc.f gfbz2r.f
     gfdmatu.f gfdpp.f gfefbz.f gfenint.f gffbz.f gfg0g.f gfg2g.f gfg2gr.f gfg2gnc.f
     gfgii.f mkpotf.f mkptfp.f pokeg0.f pokepf.f soscale.f
     gfgw.f gfibz.f gfidos.f gfidos2.f gfidosr.f gfomg.f gfp0ft.f gfp0io.f gfpdos.f
     gfqdys.f gfradp.f gfradp1.f gfradpl.f gfree.f gfrscp.f gfsbar.f
     gfsigma.f gfsubs.f gfwdiff.f gfwscrbz.f gfzerq.f gfzkin.f gibbswts.f
     gintz.f gippad.f gtibpl.f gtqval.f gvbma.f gvsym.f hamfbz.f hoffcl.f iogf.f
     iogibbs.f ioomg.f kmu2lms.f magcpa.f makdlz.f makpfz.f makpzf.f
     mixomg.f mkcpa.f mkfrpf.f mkgint.f mkpotj.f mksopf.f mstrxz.f
     nrm2gz.f p0herm.f p0k2p0.f packdmat.f pfrrot.f pgnei.f pgpdos.f pgqtot.f
     phidz.f pkgdlm.f qpitrp.f radhjz.f rdarry.f rotm.f setne.f setnth.f
     sgvsym.f specfun.f srvsym.f subdiag.f suhamz.f sxsgmtrns.f symcp.f vorbydl.f
     wrnhjz.f zseq.f pvgfevcp.f
     ''')
# gfrcls.f ioordn.f

nj.mks(gf_src) #, vs=dict(ppflags = ppflags + ' -dGF -dCGF -dSX -dNC')
nj.mk('subs/lmasa.f' , o='gf/lmasa.f.o' , vs = dict(ppflags = ppflags + ' -dGF -dCGF -dSX -dNC'))

pg_src = laddprefix('pgf/', '''
 blchpl.f lgupac.f pgbevl.f pgcops.f pgcurr.f pgdim.f pgdysn.f
 pgemb.f pgfasa.f pgfg2g.f pgfgii.f pgflu.f pgglst.f pggmal.f pginns.f pgkmap.f
 pgles.f pglusu.f pgsetg.f pgsif.f pgsoscl.f pgvpl.f plham.f rdsurfg.f supgbc.f supgh.f supghz.f
''')

nj.mks(pg_src) #, vs=dict(ppflags = ppflags + ' -dGF -uCGF -dPGF -dSX -dNC')
nj.mk('subs/lmasa.f' , o='pgf/lmasa.f.o' , vs = dict(ppflags = ppflags + ' -dGF -uCGF -dPGF -dSX -dNC'))
nj.mk('subs/atomsr.f' , o='pgf/atomsr.f.o' , vs = dict(ppflags = ppflags + ' -uALENA'))

gd_src = laddprefix('gwd/', '''bndconn_v2.f chimedit.f  chkgwin.f genqbz.f gwcphi.f mk_hamindex.f
    polinta.f pwmat.f read_bzdata.f stonerpb.f stonerrsa.f sugw.f sugwin.f rdqgcou.f gwdimpar.f suq0x.f''')

gx_src = laddprefix('gw/', '''makegw.f gwinputprodbas.f gworthphi.f iopbme.f makecgr.f mtoindex.f prodbas.f prodbaschk.f prodbasme.f
    prodphi.f rprodbas.f vcoulq.f rprodbasi.f''')

nj.mks(gd_src) #, vs=dict(ppflags = ppflags + ' -dLMFGWD')
nj.mks(gx_src) #, vs=dict(ppflags = ppflags + ' -dLMFGW')
nj.mk('fp/bndfp.f', o='gwd/bndfp.f.o', vs = dict(ppflags = ppflags + ' -dLMFGWD'))
nj.mk('fp/lmfp.f' , o='gwd/lmfp.f.o' , vs = dict(ppflags = ppflags + ' -dLMFGWD'))
nj.mk('fp/bndfp.f', o='gw/bndfp.f.o', vs = dict(ppflags = ppflags + ' -dLMFGW'))
nj.mk('fp/lmfp.f' , o='gw/lmfp.f.o' , vs = dict(ppflags = ppflags + ' -dLMFGW'))


mc_src = laddprefix('mol/', '''
    amsbb1.f amstrg.f amstrp.f astrg0.f astrgc.f astrgj.f astrgl.f astrx0.f
    astrxj.f augskj.f bisinl.f bonds.f chlr2f.f chlr2s.f cirint.f cstrx0.f
    cstrxq.f debprt.f diagrr.f distab.f djohn.f dlmtor.f dpdist.f dpmucp.f
    etacc1.f etsums.f evlwgt.f evlwsm.f evlwss.f fixcor.f ftcass.f ftcchk.f
    ftccmp.f ftcdim.f ftcdis.f ftcerr.f ftcgen.f ftchad.f ftchyf.f ftcinc.f
    ftcmak.f ftcmch.f ftcnrm.f ftcsch.f ftcseq.f ftcsev.f ftcslv.f ftcsrt.f
    ftcusc.f ftcxip.f gauskl.f gausr.f gauxpn.f getcor.f gklv.f gmeshx.f
    gpfnd1.f han1ci.f han2cg.f han2ci.f hansrf.f hansrg.f hdercg.f headln.f
    hmat2.f hmatl.f hmvi0.f hmvi3.f hoffs.f holint.f hpair3.f hrmsav.f
    hscop.f hsg2ci.f hsm1ci.f hsm2ci.f hsmatr.f hsmder.f hsmdvl.f hsmmom.f
    hsmxpn.f hsrint.f hsvecs.f hy2dm.f hy2gen.f hy2int.f hy2mak.f hy2msh.f
    hy2nrm.f hy2plt.f hy2sxi.f hyfchk.f hyfdim.f hyfevg.f hyfevl.f hyfgen.f
    hyfget.f hyfin0.f hyfin2.f hyfinp.f hyfip1.f hyfipd.f hyfipl.f hyfipt.f
    hyfloc.f hyfmak.f hyfnot.f hyfnrm.f hyfout.f hyfplt.f hyfrd.f hyfsbb.f
    hyfsbx.f hyfsch.f hyfsrt.f hyfstf.f hyftry.f hyfusc.f hyfxsc.f
    inline.f ioatom.f iojob.f iomc3.f iomc4.f iomc6.f iorodc.f iostat.f
    ious.f iunpac.f legpol.f legprj.f leqr2f.f lmca.f lmce.f lmcfit.f
    lmrho.f lmxbs.f lpfrog.f lsqfc1.f lxcf.f makusp.f maprho.f mcatm.f
    mcmixr.f mcmull.f mcrelx.f mdinf.f mdinfo.f mixr3x.f mkbs.f mopen.f
    move.f msadd.f msetp2.f msetp3.f msetpc.f msetpl.f mstrgg.f mstrpg.f
    mstrug.f mstrux.f mstrxg.f mstrxp.f nstru0.f on2err.f on2gen.f on2mak.f
    on2plt.f oncchk.f oncerr.f oncget.f oncloc.f oncmak.f oncnot.f onctry.f
    oncxpn.f outmc6.f ovldum.f parsgn.f parvls.f pauskl.f perco1.f perco2.f
    percor.f percos.f pertmb.f pertms.f pertmw.f phdtnr.f plch1.f plch2.f
    plch3.f plrho.f plwv2.f pmblk.f pnu2pl.f pnunew.f poinsp.f polft4.f
    potpm.f potpus.f potxcl.f ppot3.f premxy.f prerot.f prerxy.f prpot.f
    qifix.f qifix1.f qifix2.f qsmuz.f qsmuz3.f radint.f rdfa.f reduc1.f
    rfft3.f rhdma0.f rhdmp0.f rhdmpa.f rhdump.f rhodif.f rhoix7.f rhosum.f
    rhoti0.f rhotif.f rhouts.f rhovin.f rhoxci.f rhuti3.f rhutsl.f rodecz.f
    roisym.f ropcha.f ropchh.f ropecd.f ropech.f ropecs.f ropevx.f ropexc.f
    ropffq.f ropfxc.f ropg0.f ropgau.f rophs0.f rophs1.f rophsm.f ropj.f
    ropmro.f ropqlm.f roprfg.f roprho.f roptra.f ropvh7.f ropvho.f ropvr7.f
    ropvro.f ropxc0.f ropxcc.f ropxf7.f ropxfg.f ropxfv.f ropxvf.f ropy0m.f
    ropylm.f rovlas.f rovlaz.f rxcmro.f rxcrh1.f rxcrho.f rxcvhf.f rxcvrf.f
    setgrp.f setigv.f setrel.f setsp.f setylm.f shopos.f slhhg.f slhhpg.f
    sll0sr.f smatl.f smatm.f smrohd.f smrohx.f sofp.f sol0sr.f solhg.f
    solhpg.f solhsr.f spcfnd.f sphed3.f splpol.f sprho2.f spsm3.f spsmt.f
    spvab0.f spval3.f spvls2.f st2ary.f stewlm.f sthsv1.f strxsu.f subtls.f
    subtlz.f sumtlz.f svsvec.f sxlmnc.f symchd.f symcls.f symprj.f tcfchk.f
    tcfdim.f tcfdm.f tcferr.f tcfint.f tcfnrm.f tcfrt1.f tcfrt2.f tcfslv.f
    tcfsrt.f tcfusr.f tcfxpn.f tlstat.f tm.f trpos.f us2h.f us2hj.f
    us2sm.f vesint.f vesm1f.f vesm2f.f vesmkg.f vesmom.f vesmoz.f vprds.f
    vrobyl.f vsphrs.f vsstal.f vtailz.f vxc0gc.f vxcnlp.f vxcnls.f vxcnsp.f
    vxnlcc.f vxtrap.f writpd.f wronkj.f xcecut.f xcpmsh.f xcpmsx.f xczet0.f
    xczetx.f ylmrot.f ylmrt0.f z0head.f zthead.f
''')
nj.mks(mc_src)
#hyfzw.f



if '-time' in sys.argv: nj.printtimes()

sl_obj = laddsuffix('.o', sl_src)
if mklavail: sl_obj.append('slatsm/mkl_service.f90.o')
if h5avail: sl_obj.extend(laddsuffix('.o',laddprefix('slatsm/','h5.F90 h5.c')))
nj.mkl(sl_obj, 'sl')
nj.mkl(laddsuffix('.o', in_src), 'inp')
nj.mkl(laddsuffix('.o', lm_src), 'lm')
nj.mkl(laddsuffix('.o', sx_src + ['sx/lmasa.f']), 'sx')
nj.mkl(laddsuffix('.o', nc_src + 'nc/bndasa.f nc/lmasa.f nc/lmaux.f nc/secm2c.f nc/secmat.f nc/asaddq.f'.split()), 'nc')
nj.mkl(laddsuffix('.o', gf_src + ['gf/lmasa.f']), 'gf')
nj.mkl(laddsuffix('.o', pg_src + ['pgf/lmasa.f','pgf/atomsr.f']), 'pg')
nj.mkl(laddsuffix('.o', fp_src), 'fp')
nj.mkl(laddsuffix('.o', op_src), 'op')
nj.mkl(laddsuffix('.o', tb_src), 'tb')
nj.mkl(laddsuffix('.o', dt_src + ['dmft/lmfp.f']), 'dmft')
nj.mkl(laddsuffix('.o', gd_src + ['gwd/bndfp.f', 'gwd/lmfp.f']), 'gwd')
nj.mkl(laddsuffix('.o', gx_src + ['gw/bndfp.f', 'gw/lmfp.f']), 'gw')
nj.mkl(laddsuffix('.o', mc_src), 'mol')

nj.mk('lmv7.f', o='lm.o'     , vs=dict(ppflags = ppflags + ' -uLM -dLM -dSX -dNC'))
nj.mk('lmv7.f', o='lmstr.o'  , vs=dict(ppflags = ppflags + ' -uLM -dLMSTR'))
nj.mk('lmv7.f', o='lmf.o'    , vs=dict(ppflags = ppflags + ' -uLM -dLMF'))
nj.mk('lmv7.f', o='lmfa.o'   , vs=dict(ppflags = ppflags + ' -uLM -dLMFA'))
nj.mk('lmv7.f', o='tbe.o'    , vs=dict(ppflags = ppflags + ' -uLM -dTBE'))
nj.mk('lmv7.f', o='lmfdmft.o', vs=dict(ppflags = ppflags + ' -uLM -dLMF -dLMFDMFT'))
nj.mk('lmv7.f', o='lmgf.o'   , vs=dict(ppflags = ppflags + ' -uLM -dLMGF -dGF -dCGF -dSX -dNC'))
nj.mk('lmv7.f', o='lmpg.o'   , vs=dict(ppflags = ppflags + ' -uLM -dLM -dLMPG -dSX -dNC'))
nj.mk('lmv7.f', o='blm.o'    , vs=dict(ppflags = ppflags + ' -uLM -dBLM'))
nj.mk('lmv7.f', o='lmchk.o'  , vs=dict(ppflags = ppflags + ' -uLM -dLMCHK'))
nj.mk('lmv7.f', o='lmdos.o'  , vs=dict(ppflags = ppflags + ' -uLM -dLMDOS'))
nj.mk('lmv7.f', o='lmctl.o'  , vs=dict(ppflags = ppflags + ' -uLM -dLMCTL'))
nj.mk('lmv7.f', o='lmxbs.o'  , vs=dict(ppflags = ppflags + ' -uLM -dLMXBS'))
nj.mk('lmv7.f', o='lmplan.o' , vs=dict(ppflags = ppflags + ' -uLM -dLMPLAN'))
nj.mk('lmv7.f', o='lmscell.o', vs=dict(ppflags = ppflags + ' -uLM -dLMSCELL'))
nj.mk('lmv7.f', o='lmfgwd.o' , vs=dict(ppflags = ppflags + ' -uLM -dLMF -dLMFGWD'))
nj.mk('lmv7.f', o='lmfgw.o'  , vs=dict(ppflags = ppflags + ' -uLM -dLMF -dLMFGW'))
nj.mk('lmv7.f', o='lmfgws.o' , vs=dict(ppflags = ppflags + ' -uLM -dLMF -dLMFGWS'))
nj.mk('lmv7.f', o='lmmc.o'   , vs=dict(ppflags = ppflags + ' -uLM -dLMMC'))

nj.mks('gwd/lmf2gw_0.f gwd/lmf2gw_2.f')

nj.mks(laddprefix('utils/','fmin.f korr.f')) # , vs = dict(ppflags = ' -dMULTIMF')
nj.mks(laddprefix('plot/', 'fplot.f pldos.f plsub.f fpsub.f contour.f plbnds.f'))
nj.mks(laddprefix('utils/fit/', 'p2.f powell.f'))
nj.mk('utils/fit/pfit.f') #, vs=dict(ppflags = '-dDEBUG')
nj.mks(laddprefix('utils/', '''m_gradzr.f lm67.f spectral.f rdcmd.f
                  poscar2init.f poscar2site.f cif2init.f cif2site.f site2init.f dval.f mlist.f rdfile.f catdos.f''')
    + laddprefix('utils/mcsrc/', 'mc.f mcsubs.f rdms.f')) #, vs = dict(ppflags = '')


nj.mk('lm.o -linp -lnc -lsx -lop -llm -lsl', o='lm')
nj.mk('lmstr.o -linp -llm -lsl', o='lmstr')
nj.mk('lmf.o -linp -lfp -lop -llm -lsl', o='lmf')
if makedso:
    nj.mk('lmfa.o -linp -lfp -lop -llm -lsl', o='lmfa')
    nj.mk('lmfgwd.o -linp -lgwd -lgw -lfp -lop -llm -lsl', o ='lmfgwd')
    nj.mk('lmfgw.o -linp -lgw -lgwd -lfp -lop -llm -lsl', o ='lmfgw')
else:
    nj.mk('lmfa.o -linp -lfp -llm -lsl', o='lmfa')
    nj.mk('lmfgwd.o -linp -lgwd -lgw -lfp -llm -lsl', o ='lmfgwd')
    nj.mk('lmfgw.o -linp -lgw -lgwd -lfp -llm -lsl', o ='lmfgw')
nj.mk('lmfgws.o -linp -lfp -lop -llm -lsl', o ='lmfgws')
nj.mk('lmmc.o -linp -lmol -llm -lsl', o ='lmmc')
nj.mk('gwd/lmf2gw_0.f.o -lsl', o ='lmf2gw_0')
nj.mk('gwd/lmf2gw_2.f.o -lsl', o ='lmf2gw_2')

nj.mk('lmgf.o -linp -lgf -lnc -lsx -lop -llm -lsl', o='lmgf')
nj.mk('lmpg.o -linp -lpg -lgf -lnc -lsx -lop -lpg -llm -lsl', o='lmpg')
nj.mk('lmfdmft.o -linp -ldmft -lfp -lop -llm -lsl', o='lmfdmft')
nj.mk('tbe.o -linp -ltb -llm -lsl', o='tbe')

nj.mk('blm.o -linp -llm -lsl', o='blm')
nj.mk('lmchk.o -linp -lnc -lsx -llm -lsl', o='lmchk')
nj.mk('lmdos.o -linp -llm -lsl', o='lmdos')
nj.mk('lmctl.o -linp -llm -lsl', o='lmctl')
nj.mk('lmxbs.o -linp -llm -lsl', o='lmxbs')
nj.mk('lmplan.o -linp -llm -lsl', o='lmplan')
nj.mk('lmscell.o -linp -llm -lsl', o='lmscell')
nj.mk('utils/lm67.f.o -linp -llm -lsl', o='lm67')
nj.mk('utils/m_gradzr.f.o -llm -lsl', o='m_gradzr')
nj.mk('utils/spectral.f.o -llm -lsl', o='spectral')
nj.mk('utils/fmin.f.o utils/korr.f.o -lsl', o='fmin') # -linp
nj.mk('utils/site2init.f.o -llm -lsl', o='site2init')
nj.mk('utils/rdcmd.f.o -lsl', o='xrdcmd') # -linp
nj.mk('utils/poscar2init.f.o -lsl', o='poscar2init')
nj.mk('utils/poscar2site.f.o -lsl', o='poscar2site')
nj.mk('utils/cif2init.f.o -lsl', o='cif2init')
nj.mk('utils/cif2site.f.o -lsl', o='cif2site')
nj.mk('utils/dval.f.o -lsl', o='dval')
nj.mk('utils/mlist.f.o -lsl', o='mlist')
nj.mk('utils/catdos.f.o -lsl', o='catdos')
nj.mk('utils/rdfile.f.o -lsl', o='rdfile')
nj.mk(laddprefix('utils/mcsrc/', 'mc.f.o mcsubs.f.o rdms.f.o') + '-llm -lsl'.split(), o='mcx')
nj.mk(laddprefix('plot/', 'fplot.f.o  plsub.f.o fpsub.f.o contour.f.o') + ['-lsl'], o='fplot')
nj.mk(laddprefix('plot/', 'pldos.f.o  plsub.f.o fpsub.f.o') + ['-lsl'], o='pldos' )
nj.mk(laddprefix('plot/', 'plbnds.f.o plsub.f.o fpsub.f.o') + 'subs/ioseh.f.o -lsl'.split(), o='plbnds')
nj.mk(laddprefix('utils/fit/', 'pfit.f.o p2.f.o powell.f.o') + ['-lsl'], o='pfit')

nj.mks('dmft/broad_sig.f90'.split())
nj.mk('dmft/broad_sig.f90.o'.split(), o='broad_sig')
if h5avail:
    nj.mks('dmft/readvertex.f90'.split())
    nj.mk(['dmft/readvertex.f90.o  ']+laddprefix('slatsm/','h5.F90.o h5.c.o init.f.o ')+laddprefix('dmft/','vertex.f90.o') , o='readvertex')


# gw stuff

nj.fppstyle = 'cpp'
nj.eincs.update(set(['mpif.h']))

#fflags = nj.flags['fflags'] + ' ' + nj.flags['omp'] + ' ' + nj.flags['gwflags']
gwfflags = '$fflags $omp $gwflags'

cudavail = '-DCUD' in nj.flags['gwflags']


derfc = laddprefix('fpgw/nfpsrc/', 'derfc.F d1mach.F i1mach.F')


chk = laddprefix('fpgw/', 'main/hchknw.m.F gwsrc/genallcf_dump.F')

phig = laddprefix('fpgw/', '''
    Miyake/maxloc/hphig.F gwsrc/wcf.F gwsrc/tetwt4.F gwsrc/tetwt5.F gwsrc/gintxx.F
    gwsrc/pplmat.F gwsrc/getgv2.F gwsrc/x0kf_v2h.F gwsrc/ppbafp.fal_oldroutines.F
    gwsrc/rs.F nfpsrc/u_lat_0.F nfpsrc/wronkj.F nfpsrc/besslr.F nfpsrc/dbesnu.F
    nfpsrc/dgamma.F nfpsrc/dibesl.F nfpsrc/djbesl.F nfpsrc/dkbesl.F nfpsrc/dybesl.F''')
#fpgw/gwsrc/x0kf.F


qpe = laddprefix('fpgw/', '''gwsrc/switches.F gwsrc/keyvalue.F main/hqpe.m.F gwsrc/qpe1.F
    gwsrc/icompvv2.F gwsrc/iopenxx.F gwsrc/iopen.F gwsrc/rw.F gwsrc/rydberg.F''')


eo = laddprefix('fpgw/', 'tote/eout.F gwsrc/rydberg.F')


psig = laddprefix('fpgw/', '''Miyake/maxloc/hpsig.F gwsrc/wcf.F gwsrc/tetwt4.F gwsrc/x0kf.F
    gwsrc/tetwt5.F gwsrc/gintxx.F gwsrc/pplmat.F gwsrc/getgv2.F gwsrc/x0kf_v2h.F
    gwsrc/ppbafp.fal_oldroutines.F gwsrc/rs.F nfpsrc/u_lat_0.F nfpsrc/wronkj.F
    nfpsrc/besslr.F nfpsrc/dbesnu.F nfpsrc/dgamma.F nfpsrc/dibesl.F nfpsrc/djbesl.F nfpsrc/dkbesl.F nfpsrc/dybesl.F''')


convg = 'fpgw/main/convgwin.F'.split()


sxc = laddprefix('fpgw/', '''main/hsfp0.m.F gwsrc/wse.F gwsrc/wintzsg.F
    gwsrc/sxcf_fal2.F gwsrc/mkmelt.F gwsrc/bzints2.F gwsrc/genallcf_dump.F
    gwsrc/ppbafp.fal_oldroutines.F nfpsrc/diagcv2.F''')


uu = laddprefix('fpgw/', '''main/h_uumatrix.m.F gwsrc/wcf.F gwsrc/tetwt4.F gwsrc/tetwt5.F
    gwsrc/gintxx.F gwsrc/pplmat.F gwsrc/getgv2.F gwsrc/x0kf_v2h.F gwsrc/psi2b_v4.F
    gwsrc/hilbertmat.F gwsrc/x0kf.F gwsrc/mkzxq.F gwsrc/rs.F nfpsrc/u_lat_0.F
    nfpsrc/wronkj.F nfpsrc/besslr.F nfpsrc/dbesnu.F nfpsrc/dgamma.F nfpsrc/dibesl.F
    nfpsrc/djbesl.F nfpsrc/dkbesl.F nfpsrc/dybesl.F gwsrc/ppbafp.fal_oldroutines.F
    nfpsrc/diagcv2.F''')
#fpgw/gwsrc/x0kf_v3h.F gwsrc/cross.F

gw0 = laddprefix('fpgw/gwsrc/', '''MPI_fpgw.F m_hamindex.F readpomat.F keyvalue.F rppovl.F
    nocctotg.F ppbafp.fal.F psi2b_v2.F psi2b_v3.F meltpb.F wfacx.F sortea.F rydberg.F
    polinta.F efsimplef.F extension.F rangedq.F nword_in_dble.F scg.F matm.F rdpp.F mptauof.F w0w0i.F
    genallcf_mod.F rgwinf_mod.F rotdlmm.F iopen.F cputid.F rw.F ext.F ext2.F cross.F
    mate.F mate1.F icopy.F bib1.F index.F idxk.F maxnn.F reindx.F

    iprint.F bz.F bzmesh.F genqbz.F linpackdummy.F switches.F rwbzdata.F llnew.F readeigen.F locqp.F
    readqg.F zprm.F iqindx.F idalloc.F

    thdist.F mpi_mod.F90''') + laddprefix('fpgw/',
    'nfpsrc/mklegw.F nfpsrc/dpzero.F slatsmlib/dlength.F slatsmlib/dvheap.F slatsmlib/huntx.F slatsmlib/ylmrtg.F slatsmlib/ropyln.F slatsmlib/ropyln1.F slatsmlib/ropcsm.F')
#fpgw/gwsrc/pointops.F
#fpgw/gwsrc/iolib.F
#fpgw/slatsmlib/dlength.F
#fpgw/nfpsrc/rxx.F
#fpgw/gwsrc/mopen.F

qg = laddprefix('fpgw/','''gwsrc/conv2gwinput.F main/qg4gw.m.F gwsrc/getbzdata1.F gwsrc/tetfbz.F
 gwsrc/m_q0p.F gwsrc/mkqg.F gwsrc/q0irre.F gwsrc/getgv2.F gwsrc/tetwt4.F gwsrc/tetwt5.F gwsrc/zsvd.F
 nfpsrc/latvec.F''')


bas = ['fpgw/main/hbasfp0.m.F'] + laddprefix('fpgw/gwsrc/','''reindx.F maxnn.F
    icopy.F basnfp.F rgwinf_mod.F keyvalue.F switches.F gintxx.F rs.F ext.F
    iopen.F excore.F rydberg.F extension.F rangedq.F polinta.F llnew.F''') + laddprefix('fpgw/nfpsrc/','''besslr.F dbesnu.F dgamma.F dibesl.F djbesl.F dkbesl.F dybesl.F hansmr.F''')


vcc = ['fpgw/main/hvccfp0.m.F'] + laddprefix('fpgw/gwsrc/','''vcoulq.F radmwt.F
    radmsh.F ropbes.F wronjjv.F mkjp.F gintxx.F extension.F rangedq.F keyvalue.F
    switches.F strxq.F iopen.F pplmat.F matm.F getgv2.F cross.F llnew.F readqg.F
    iqindx.F cputid.F idalloc.F thdist.F dqsmpy.F zprm.F mpi_mod.F90''') + laddprefix('fpgw/nfpsrc/','''besslr.F dbesnu.F dgamma.F dibesl.F djbesl.F dkbesl.F dybesl.F''')

wmat = laddprefix('fpgw/','''Miyake/maxloc/hwmat.F Miyake/maxloc/maxloc0.F
 gwsrc/wse.F Miyake/maxloc/wmat.F gwsrc/ppbafp.fal_oldroutines.F gwsrc/genallcf_dump.F''')

sigmconv = laddprefix('fpgw/gwsrc/','switches.F keyvalue.F iopen.F') + ['fpgw/main/hsigmconv.m.F']

sxc_sc = ['fpgw/main/hsfp0.sc.m.F'] + laddprefix('fpgw/gwsrc/','''wse.F sxcf.sc.F iosigh.F
 mkmelt.F bzints2.F wintzsg.F mkseci.F mksecp.F ppbafp.fal_oldroutines.F zqsmpy.F''') + ['fpgw/nfpsrc/diagcv2.F'] + ['fpgw/slatsmlib/dpdump.F']
#fpgw/gwsrc/sxcf_fal2.F

parainfo = laddprefix('fpgw/', 'main/hparainfo.m.F gwsrc/charext.F')

merge = ['fpgw/main/hmergewv.m.F'] + laddprefix('fpgw/gwsrc/','switches.F keyvalue.F iopen.F')


uu2 = ['fpgw/Miyake/maxloc/huumat.F'] + laddprefix('fpgw/gwsrc/','''wcf.F tetwt4.F tetwt5.F
 gintxx.F pplmat.F getgv2.F x0kf_v2h.F psi2b_v4.F rs.F ppbafp.fal_oldroutines.F''') + laddprefix('fpgw/nfpsrc/','u_lat_0.F wronkj.F besslr.F dbesnu.F dgamma.F dibesl.F djbesl.F dkbesl.F dybesl.F''')


#gw0tot = laddprefix('fpgw/gwsrc/','''rwbzdata.F keyvalue.F genallcf_mod.F rgwinf_mod.F
 #nocctotg.F ppbafp.fal.F psi2b_v2.F psi2b_v3.F wfacx.F sortea.F rydberg.F polinta.F
 #efsimplef.F extension.F rangedq.F nword_in_dble.F scg.F matm.F rdpp.F mptauof.F rotdlmm.F
 #iopen.F cputid.F rw.F ext.F ext2.F cross.F mate.F mate1.F icopy.F bib1.F index.F
 #idxk.F maxnn.F reindx.F pointops.F iolib.F iprint.F bz.F bzmesh.F genqbz.F switches.F
 #linpackdummy.F rppovl.F thdist.F zprm.F llnew.F''')


kino_input_test = ['fpgw/main/kino_input_test.F']


bndout = ['fpgw/main/hbndout.m.F'] + laddprefix('fpgw/gwsrc/','''iqagree.F iopenxx.F iopen.F polinta.F
 rydberg.F extension.F rangedq.F switches.F keyvalue.F''')


mloc = laddprefix('fpgw/Miyake/maxloc/', 'hmaxloc.F maxloc0.F maxloc1.F maxloc2.F maxloc3.F') + laddprefix('fpgw/gwsrc/','wse.F genallcf_dump.F')



nfpl = laddprefix('fpgw/nfpsrc/', '''wronkj.F sylm.F sylmnc.F u_lat_0.F mklegw.F setpr.F

 hsmq.F lgen.F hansr5.F hansr4.F lattc.F ll.F dpcopy.F dpadd.F qdist.F dlmtor.F dpzero.F ropyln.F ropcsm.F
 dsisl.F dsifa.F diagcv2.F''') + ['fpgw/gwsrc/scg.F']

rdat_v2 = ['fpgw/main/rdata4gw_v2.m.F'] + laddprefix('fpgw/gwsrc/', '''
 keyvalue.F switches.F rwbzdata.F gintxx.F cinvrx.F idxk.F nword_in_dble.F gwinput_v2.F
 matm.F getgv2.F iopen.F pplmat.F bzmesh.F ext.F ext2.F cross.F rs.F extension.F iofa.F
 rangedq.F llnew.F iqindx.F polinta.F mpi_mod.F90''') + laddprefix('fpgw/nfpsrc/','''besslr.F dbesnu.F dgamma.F dibesl.F djbesl.F dkbesl.F dybesl.F''') \
 + laddprefix('fpgw/slatsmlib/', 'dfdump.F')

qpe_sc = ['fpgw/main/hqpe.sc.m.F'] + laddprefix('fpgw/gwsrc/','''switches.F keyvalue.F
 qpe1.sc.F icompvv2.F iopenxx.F iopen.F iosigh.F rw.F rydberg.F iprint.F rwbzdata.F zprm.F mpi_mod.F90''') \
 + laddprefix('fpgw/slatsmlib/', 'dmcpy.F dsifa.F dsisl.F dsidi.F amix.F')

heftet = laddprefix('fpgw/', 'main/heftet.m.F gwsrc/bzints2.F')

nfpltot = ['fpgw/nfpsrc/diagcv2.F']

comm = laddprefix('fpgw/', 'nfpsrc/rxx.F gwsrc/mopen.F slatsmlib/words.F slatsmlib/zmscop.F slatsmlib/zmsadd.F nfpsrc/syscalls2.F')

ecor = laddprefix('fpgw/tote/', 'hecor.F rpaq.F')

x0 = laddprefix('fpgw/','main/hx0fp0.m.F nfpsrc/diagcv2.F tote/rpaq.F') + laddprefix('fpgw/gwsrc/', '''
 wcf.F tetwt4.F x0kf.F tetwt5.F hilbertmat.F x0kr.F x0k_sym.F mkzxq.F cinvrx.F zsvd.F mkmelt.F''')
    #fpgw/gwsrc/x0kf_v2h.F
    #fpgw/gwsrc/x0kf_v3h.F

eo2 = laddprefix('fpgw/','tote/eout2.F gwsrc/rydberg.F')

x0mlw = laddprefix('fpgw/','''Miyake/maxloc/hx0fp0.m.F Miyake/maxloc/wcf.F gwsrc/tetwt4.F
 gwsrc/x0kf.F gwsrc/tetwt5.F gwsrc/x0kf_v2h.F Miyake/maxloc/x0kf_v3h.F nfpsrc/diagcv2.F
 tote/rpaq.F gwsrc/ppbafp.fal_oldroutines.F gwsrc/cinvrx.F''')


#mloc1d = '''
    #fpgw/Miyake/maxloc/maxloc0.F
    #fpgw/Miyake/maxloc/maxloc1.F
    #fpgw/Miyake/maxloc/maxloc2.F
    #fpgw/Miyake/maxloc/maxloc3.F
    #fpgw/gwsrc/wse.F
    #fpgw/gwsrc/genallcf_dump.F
#'''.split()
##fpgw/Miyake/maxloc/hmaxloc1D.F

hnocc_mlw = laddprefix('fpgw/','Miyake/maxloc/hnocc_mlw.F gwsrc/bzints2.F')

x0_sc = laddprefix('fpgw/','main/hx0fp0.sc.m.F nfpsrc/diagcv2.F') + laddprefix('fpgw/gwsrc/','''wcf.F tetwt4.F mkmelt.F
 tetwt5.F x0kf.F hilbertmat.F x0kr.F x0k_sym.F mkzxq.F''')
    #fpgw/gwsrc/x0k.F
    #fpgw/gwsrc/x0kf.F
    #fpgw/gwsrc/x0kf_v2h.F
    #fpgw/gwsrc/x0kf_v3h.F

hef = laddprefix('fpgw/','main/hef.m.F gwsrc/wse.F')

qpwf = ['fpgw/Miyake/maxloc/qpwf.F']


bse = laddprefix('fpgw/gwsrc/','''bib1.F bz.F cputid.F cross.F ext.F idalloc.F iopenxx.F
 iprint.F iqindx.F maxnn.F mopen.F rangedq.F reindx.F rw.F sortea.F thdist.F efsimplef.F
 extension.F genallcf_mod.F genqbz.F idxk.F index.F iopen.F keyvalue.F llnew.F mate.F
 matm.F meltpb.F m_hamindex.F mkmelt.F mptauof.F w0w0i.F ppbafp.fal.F rdpp.F readeigen.F locqp.F readqg.F
 rgwinf_mod.F rotdlmm.F rppovl.F rwbzdata.F rydberg.F scg.F switches.F mpi_mod.F90''')
bse += laddprefix('fpgw/bse/', 'nocctotg_bc.F kernel.F inverse_bc.F diag_bc.F bethesalpeter.F')
bse += laddprefix('fpgw/nfpsrc/','diagcv2.F syscalls2.F mklegw.F dpzero.F rxx.F')
bse += laddprefix('fpgw/slatsmlib/','ylmrtg.F ropyln.F ropyln1.F ropcsm.F dlength.F dvheap.F huntx.F')

cu = '''
    fpgw/gwsrc/fcuzwz.F90
'''.split()

if cudavail:
    cu.append('fpgw/gwsrc/zwz.cxx')


nj.mks(derfc + chk + phig + qpe + eo + psig + convg + sxc + uu + gw0 + qg + bas + vcc + wmat + sigmconv + sxc_sc + parainfo + merge + uu2 + kino_input_test + bndout + mloc + nfpl +  rdat_v2 + qpe_sc + heftet + nfpltot + comm + ecor + x0 + eo2 + x0mlw  + hnocc_mlw + x0_sc + hef + qpwf + bse + cu, vs=dict(fflags=gwfflags)) # + gw0tot + mloc1d


if '-time' in sys.argv: nj.printtimes()

for i in 'derfc gw0 nfpl nfpltot comm'.split(): # gw0tot mloc1d
    nj.mk(laddsuffix('.o', eval(i)), o='lib/liblmgw_'+i+'.a')


def njmk(name, group, libs):
    nj.mk(laddsuffix('.o', group) + list(map(lambda i: 'lib/liblmgw_'+i+'.a', libs.split())), o=name, vs=dict(fflags=gwfflags))

njmk('code2/hecor'          , ecor           , 'nfpltot gw0 comm')
njmk('code2/eout'           , eo             , 'comm            ')
njmk('code2/eout2'          , eo2            , 'comm            ')
njmk('code2/hsigmconv'      , sigmconv       , 'comm            ')
njmk('code2/qpwf'           , qpwf           , 'gw0 comm        ')
njmk('code2/qg4gw'          , qg             , 'gw0 comm        ')
njmk('code2/rdata4gw_v2'    , rdat_v2        , 'nfpl comm       ')
njmk('code2/hbasfp0'        , bas            , 'comm            ')
#njmk('code2/hvccfp0'        , vcc            , 'nfpl derfc comm ')
njmk('code2/hx0fp0'         , x0             , 'gw0 comm        ')
njmk('code2/hx0fp0_mlw'     , x0mlw          , 'gw0 comm        ')
njmk('code2/h_uumatrix'     , uu             , 'gw0 comm        ')
njmk('code2/huumat'         , uu2            , 'gw0 comm        ')
njmk('code2/hphig'          , phig           , 'gw0 comm     ') #lsl
njmk('code2/hpsig'          , psig           , 'gw0 comm        ')
njmk('code2/hx0fp0_sc'      , x0_sc          , 'gw0 comm        ')
njmk('code2/hwmat'          , wmat           , 'gw0 comm        ')
njmk('code2/hmaxloc'        , mloc           , 'nfpltot gw0 comm')
njmk('code2/hsfp0'          , sxc            , 'gw0 comm        ')
njmk('code2/hsfp0_sc'       , sxc_sc+cu      , 'gw0 comm        ')
njmk('code2/hnocc_mlw'      , hnocc_mlw      , 'gw0 comm        ')
njmk('code2/heftet'         , heftet         , 'gw0 comm        ')
njmk('code2/hef'            , hef            , 'gw0 comm        ')
njmk('code2/hchknw'         , chk            , 'gw0 comm        ')
njmk('code2/hqpe'           , qpe            , 'comm            ')
njmk('code2/hqpe_sc'        , qpe_sc         , 'comm            ')
njmk('code2/hmergewv'       , merge          , 'comm            ')
njmk('code2/hparainfo'      , parainfo       , 'gw0 comm        ')
njmk('code2/hbndout'        , bndout         , 'comm            ')
njmk('code2/convgwin'       , convg          , 'comm            ')
njmk('code2/kino_input_test', kino_input_test, 'gw0 comm        ')
njmk('code2/bethesalpeter'  , bse            , 'comm            ')
#njmk('code2/hmaxloc1D'      , mloc1d         , 'nfpltot gw0 comm')



omp_src = laddprefix('fpgw/main/','hx0fp0.m.F hx0fp0.sc.m.F hsfp0.m.F hsfp0.sc.m.F h_uumatrix.m.F hvccfp0.m.F') + laddprefix('fpgw/gwsrc/', 'ppbafp.fal.F mkzxq.F mkseci.F mksecp.F wintzsg.F radmsh.F radmwt.F vcoulq.F ropbes.F wronjjv.F sxcf.sc.F meltpb.F mkjp.F') + laddprefix('fpgw/nfpsrc/','syscalls2.F sylm.F wronkj.F ropyln.F')

ompflags = gwfflags+' -DOPENMP'

for i in omp_src:
    nj.mk(i, o=i+'_om.o', vs=dict(fflags=ompflags))

def njmk_om(name, srcs):
    ssrcs = ' '.join(srcs)
    for i in omp_src:
        ssrcs = ssrcs.replace(i,i+'_om')
    nj.mk(laddsuffix('.o', ssrcs), o=name, vs=dict(fflags=ompflags))

njmk_om('code2/hx0fp0_om'    , x0     + gw0 + comm)
njmk_om('code2/h_uumatrix_om', uu     + gw0 + comm)
njmk_om('code2/hx0fp0_sc_om' , x0_sc  + gw0 + comm)
njmk_om('code2/hsfp0_om'     , sxc    + gw0 + comm)
njmk_om('code2/hsfp0_sc_om'  , sxc_sc + cu + gw0 + comm)
njmk_om('code2/hvccfp0'      , vcc    + nfpl + derfc + comm)

nj.bld(['$srcpath/fpgw/mark/lmgw'],'cp',['$bldpath/lmgw'],si=['$bldpath/mpicmd'])
nj.bld([],'cmd',['$bldpath/mpicmd'],vs=dict(command='echo set mpirun = \\"$mpirun\\" > $bldpath/mpicmd'))

scripts = '''lmgwsc lmgw1-shot lmgwclear infgw savegwfiles
            analyze_lx0 gw-extract-prodbas-and-time-from-output
            gw-extract-qp-from-spex-output'''.split()

for script in scripts:
    nj.bld(['$srcpath/fpgw/mark/' + script], 'cp', ['$bldpath/'+script])

nj.bld(['$srcpath/fpgw/mark/hqpemetal'],'cp', ['$bldpath/code2/hqpemetal'])

cm_src = laddprefix('fpgw/cumulant/','''m_par.f90 par.f90 polint.f90 terp.f90 terpc.f90 spfcn2.f90''')
nj.mks(cm_src)
nj.mk(laddsuffix('.o', cm_src), o='spfcn')


if (nj.srcpath != nj.bldpath): nj.mkincludes()

scriptutils = 'rdcmd vextract band-edge absolute-path catf poszer add0 extract-lines SpectralFunction.sh band-plot'.split()
for i in scriptutils:
    nj.bld(['$srcpath/utils/'+i],'cp',['$bldpath/'+i])


nj.bld(['$srcpath/dmft/interface_solver/lmtriqs'],'cp',['$bldpath/lmtriqs'])
nj.bld(['$srcpath/dmft/interface_solver/average_sig'],'cp',['$bldpath/average_sig'])
nj.bld(['$srcpath/dmft/interface_solver/impurity.py'],'cp',['$bldpath/impurity.py'])
nj.bld(['$srcpath/dmft/sus/susceptibility'],'cp',['$bldpath/susceptibility'])
utils = 'm_gradzr fmin poscar2init site2init poscar2site catdos cif2init cif2site dval mlist rdfile mcx fplot pldos plbnds pfit broad_sig lmtriqs susceptibility average_sig'.split() + scriptutils


allgw = laddprefix('code2/','hecor eout eout2 hsigmconv qpwf qg4gw rdata4gw_v2 hbasfp0 hvccfp0 hx0fp0 hx0fp0_mlw h_uumatrix huumat hphig hpsig hx0fp0_sc hwmat hmaxloc hsfp0 hsfp0_sc hnocc_mlw heftet hef hchknw hqpe hqpe_sc hmergewv hparainfo hbndout convgwin kino_input_test hx0fp0_om h_uumatrix_om hx0fp0_sc_om hsfp0_om hsfp0_sc_om hqpemetal bethesalpeter') + scripts + 'spfcn lmgw'.split()

nj.bld(laddprefix('$bldpath/',allgw),'phony',['all-gw'])

progs = '''lm lmstr lmf lmfa lmfdmft tbe lmgf lmpg blm lmchk lmdos lmctl lmscell lm67 xrdcmd lmxbs lmplan lmfgwd lmfgw lmfgws spectral lmf2gw_0 lmf2gw_2'''.split() # lmmc
nj.bld(laddprefix('$bldpath/',utils),'phony',['utils'])
nj.bld(laddprefix('$bldpath/',progs)+['utils','all-gw'],'phony',['all']) #set(progs)-set(['lmmc'])

conf = 'mpicmd impurity.py'.split()
for i in conf:
    nj.bld(['$bldpath/'+i],'cp',['$prefix/bin/'+i])

allprogs = utils+progs+allgw
for i in allprogs:
    nj.bld(['$bldpath/'+i],'install',['$prefix/bin/'+i])
    #nj.bld(['$bldpath/'+i],'cp',['$prefix/bin/'+i])

allfiles = laddprefix('$prefix/bin/', allprogs + conf)

if makedso:
    ilibs = lsandwich(nj.libpath+'/lib','sl inp lm sx nc gf pg fp op tb dmft gwd gw',nj.libsfx) # mol
    for i in ilibs: nj.bld(['$bldpath/'+i],'install',['$prefix/'+i])
    allfiles.extend(laddprefix('$prefix/',ilibs))

nj.bld(allfiles,'phony',['install'])
nj.bld([],'cmd',['uninstall'],vs=dict(
    command='rm -f '+' '.join(allfiles)
))

# checks
for i in 'tm testrun.sh stat-tests.py'.split():
    nj.bld(['$srcpath/utils/'+i],'cp',['$bldpath/'+i])

nj.bld(['$srcpath/testing/eref.dat'],'cp',['$bldpath/eref.dat'])
nj.bld(laddprefix('$bldpath/','eref.dat stat-tests.py testrun.sh tm'),'phony',['testutils'])

nj.bld(['$srcpath/testing'], 'cmd', ['$bldpath/testing'], vs=dict(command='cp -r $srcpath/testing $bldpath/'))
for m in 'fp gf pgf sx optics tb mol dmft nc gwd gw'.split():
    nj.bld(['$srcpath/'+m+'/test'], 'cmd', ['$bldpath/'+m+'/test'], vs=dict(command='cp -r $srcpath/'+m+'/test $bldpath/'+m+'/'))

nj.bld(['testutils']+laddprefix('$bldpath/',['testing']+laddsuffix('/test','fp gf pgf sx optics tb mol dmft nc gwd gw')),'phony',['testpreps'])
nj.bld(['testpreps'],'phony',['prep-tests'])


tests_s = loadfile(os.path.dirname(sys.argv[0])+'/testable.md')
testtable_start = re.search(r'(?xm)^\s*\|(-+\|){4}\s*$',tests_s).end()
ts = re.findall(r'(?xm)^\s*\|([^|]+?)\|([^|]+?)\|([^|]*?)\|([^|]*?)\|\s*$',tests_s[testtable_start:])
tsd = []

for t in ts:
    cmd, path, np, tags = (i.strip() for i in t)
    if cmd != '' and 'disable' not in tags and (mpiavail or not (mpiavail or 'mpik' in tags)):
        if np == '': np = 1
        tsd.append(dict(
            cmd = '$bldpath/'+cmd,
            rd  = '$bldpath/checks/'+path+'.d',
            log  = '$bldpath/checks/'+path+'.log',
            np  = int(np),
            path= path,
            tags= tags.split()))

del ts

upaths = {}
for t in tsd:
    path, cmd = t['path'], t['cmd']
    if path not in upaths:
        upaths[path] = cmd
    else:
        sys.stderr.write('\nWarning: Tests "{cmd1}" and "{cmd2}" run in the same path {path}, expect trouble!\n'.format(cmd1=upaths[path],cmd2=cmd,path=path))
del upaths

ncore = 1
if sys.platform.startswith('linux'):
    #ncore = int(subprocess.check_output('lscpu -p=core'.split()).rsplit(None,1)[1])+1
    lscpu_p = subprocess.Popen('lscpu -p=core'.split(), shell=False, stdout=subprocess.PIPE)
    #ncore = int(lscpu_p.communicate()[0].rsplit(None,1)[1])+1
    ncore = int(lscpu_p.communicate()[0].rsplit(None,1)[1].split(b',',1)[0])+1
elif sys.platform.startswith('darwin'):
    sysctl_p = subprocess.Popen('sysctl -n hw.physicalcpu'.split(), shell=False, stdout=subprocess.PIPE)
    ncore = int(sysctl_p.communicate()[0])
else:
    print ('We have not tested on your OS!')
    print ('Cannot find the number of cores on your system and will stay on the safe side by assuming %d'%ncore)


nj.js.append('pool pool-tests\n    depth = %d' % (ncore))
#nj.js.append('pool pool-tests\n    depth = 4')

#npmx = max(t['np'] for t in tsd)
#if npmx > ncore:
    #sys.stderr.write('Warning: Parallel tests requesting more than %i cores may not run well on this machine.\n' % ncore)

testcmd = '$bldpath/testrun.sh {rd} {np:d} {cmd} &> /dev/null' # the redirection is a workaround to some glitch in ninja, without it ^C does not kill processes
#testcmd = '$bldpath/testrun.sh {rd} {np:d} {cmd} &>> oo'

for t in tsd:
    cmd = testcmd.format(**t)

    nj.bld([],'cmd',['test-'+t['path']+'-show'],vs=dict(
        command = 'echo "' + cmd.replace('"',r'\"') + '"',
        description = 'test-' + t['path'],
        #description = cmd,
        pool    = 'console',
    ))

    vs = dict(
        command = cmd,
        pool    = 'pool-tests',
        #pool    = 'console',
    )

    #if t['np'] > 1: vs['weight'] = '%d' % t['np']
    nj.bld([],'cmd',[t['rd']],so=['testpreps'],vs=vs)
    nj.bld([t['rd']],'phony',['test-'+t['path']])


def test_aliases(tag, tsd):
    if tag != '': tag = '-'+tag

    tests_all = set(t['path'] for t in tsd)
    tests_mpi = set(t['path'] for t in tsd if 'mpik' in t['tags'])
    tests_xtr = set(t['path'] for t in tsd if 'extra' in t['tags'])

    simple = tests_all - tests_xtr - tests_mpi
    extra  = tests_xtr - tests_mpi
    simple_mpi = tests_mpi - tests_xtr
    mpi_extra = tests_mpi & tests_xtr

    if tests_all  != set(): nj.bld(laddprefix('test-',tests_all ),'phony',['test'+tag+'-all'])
    if simple     != set(): nj.bld(laddprefix('test-',simple    ),'phony',['test'+tag])
    if extra      != set(): nj.bld(laddprefix('test-',extra     ),'phony',['test'+tag+'-extra'])
    if simple_mpi != set(): nj.bld(laddprefix('test-',simple_mpi),'phony',['test'+tag+'-mpi'])
    if mpi_extra  != set(): nj.bld(laddprefix('test-',mpi_extra ),'phony',['test'+tag+'-mpi-extra'])

    if tests_all  != set(): nj.bld(laddsuffix('-show',laddprefix('test-',tests_all )),'phony',['test'+tag+'-all-show'])
    if simple     != set(): nj.bld(laddsuffix('-show',laddprefix('test-',simple    )),'phony',['test'+tag+'-show'])
    if extra      != set(): nj.bld(laddsuffix('-show',laddprefix('test-',extra     )),'phony',['test'+tag+'-extra-show'])
    if simple_mpi != set(): nj.bld(laddsuffix('-show',laddprefix('test-',simple_mpi)),'phony',['test'+tag+'-mpi-show'])
    if mpi_extra  != set(): nj.bld(laddsuffix('-show',laddprefix('test-',mpi_extra )),'phony',['test'+tag+'-mpi-extra-show'])

tsda = [t for t in tsd if not t['path'].startswith('mol_')] # filter out mol tests from the aliases
tsde = ''
if nj.flags['qcmd'] == 'env' or nj.flags['qcmd'].startswith('tm -n'):
    tsde = ', '.join([t['path'] for t in tsda if t['np'] > ncore])
    tsda = [t for t in tsda if t['np'] <= ncore]


if tsde != '':
    tsde = tsde[::-1].replace(',','dna ',1)[::-1]
    sys.stderr.write('\nWarning: The following parallel tests request more cores than available on this machine and may not run well. They will be excluded from the standard test groups but may still be run with if explicitly specified on the command line and if qcmd = env is set in the flags.mk file. %s\n\n' % (tsde) )

test_aliases('',tsda)
prog_tags = set(t['path'].split('_',1)[0] for t in tsda)

for ptag in prog_tags:
    ptag_subset = [t for t in tsda if t['path'].startswith(ptag+'_')]
    test_aliases(ptag, ptag_subset)

nj.bld(['$bldpath/stat-tests.py'],'cmd',['stat-tests'],vs=dict(
        command = '$bldpath/stat-tests.py',
        description = '',
        pool    = 'console',
    ))
nj.bld(['$bldpath/stat-tests.py'],'cmd',['clear-tests'],vs=dict(
        command = "r=$$($bldpath/stat-tests.py | grep running | awk '{print $$1}' | sed s:'\.log$$':'.*':); [ -z \"$$r\" ] || rm -rf $$r",
        description = '',
        pool    = 'console',
    ))


if nj.flags['qcmd'].startswith('tm -n'):
    nj.bld(['$bldpath/tm'],'cmd',['stop-tests'],vs=dict(
            command = '$bldpath/tm -n 0 q',
            description = '',
            pool    = 'console',
        ))


nj.bld([],'cmd',['clean'],vs=dict(command='ninja -t clean'))

nj.direct += '\ndefault all\n'

nj.write(ofln,lw=96)
