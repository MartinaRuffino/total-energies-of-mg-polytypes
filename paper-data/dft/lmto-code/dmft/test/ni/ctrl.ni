TDMRHO
 rm -f mixm.ni
 lmf ni -vnk=12 --rs=1,0,0,1 -vnsp=2 -vtpd1=4 -vbigbas=3 -vehmax=-.4 -vrwa=0.3359 -vemaxs=2.5 -vconvc=1e-6 -vbeta=.3 --ldadc=82.2 -job=1 --quit=rho --pr40,20
 rm -f mixm.ni
 lmf ni -vnk=12 --rs=1,0,0,1 -vnsp=2 -vtpd1=4 -vbigbas=3 -vehmax=-.4 -vrwa=0.3359 -vemaxs=2.5 -vconvc=1e-6 -vbeta=.3 --ldadc=82.2 -job=1 -vkt=1/50 --quit=rho --pr40,20
 rm -f mixm.ni
 lmfdmft ni -vnk=12 --rs=1,0 -vnsp=2 -vtpd1=4 -vbigbas=3 -vehmax=-.4 -vrwa=0.3359 -vemaxs=2.5 -vconvc=1e-6 -vbeta=.3 --ldadc=82.2 -job=1 --pr40,20
 lmfdmft ni -vnk=12 --rs=1,2,0,1 -vnsp=2 -vtpd1=4 -vbigbas=3 -vehmax=-.4 -vrwa=0.3359 -vemaxs=2.5 -vconvc=1e-6 -vbeta=.3 --ldadc=82.2 -job=1 -vkt=1/50 --udrs --pr40,20
TMDMRHO
 rm -f mixm.ni
 mpirun -n 6 lmf ni -vnk=12 --rs=1,0,0,1 -vnsp=2 -vtpd1=4 -vbigbas=3 -vehmax=-.4 -vrwa=0.3359 -vemaxs=2.5 -vconvc=1e-6 -vbeta=.3 --ldadc=82.2 -job=1 --quit=rho --pr40,20
 rm -f mixm.ni
 mpirun -n 6 lmf ni -vnk=12 --rs=1,0,0,1 -vnsp=2 -vtpd1=4 -vbigbas=3 -vehmax=-.4 -vrwa=0.3359 -vemaxs=2.5 -vconvc=1e-6 -vbeta=.3 --ldadc=82.2 -job=1 -vkt=1/50 --quit=rho --pr40,20
 rm -f mixm.ni
 lmfdmft ni -vnk=12 --rs=1,0 -vnsp=2 -vtpd1=4 -vbigbas=3 -vehmax=-.4 -vrwa=0.3359 -vemaxs=2.5 -vconvc=1e-6 -vbeta=.3 --ldadc=82.2 -job=1 --pr40,20
 mpirun -n 6 lmfdmft ni -vnk=12 --rs=1,2,0,1 -vnsp=2 -vtpd1=4 -vbigbas=3 -vehmax=-.4 -vrwa=0.3359 -vemaxs=2.5 -vconvc=1e-6 -vbeta=.3 --ldadc=82.2 -job=1 -vkt=1/50 --udrs --pr40,20
TSTDMFT
 lmfdmft ni -vnk=12 --rs=1,0 -vnsp=2 -vtpd1=4 -vbigbas=3 -vehmax=-.4 -vrwa=0.3359 -vemaxs=2.5 -vconvc=1e-6 -vbeta=.3 --ldadc=82.2 -job=1 -vindmfl=f
 lmfdmft ni -vnk=12 --rs=1,0 -vnsp=2 -vtpd1=4 -vbigbas=3 -vehmax=-.4 -vrwa=0.3359 -vemaxs=2.5 -vconvc=1e-6 -vbeta=.3 --ldadc=82.2 -job=1 -vindmfl=f
 lmfdmft ni -vnk=12 --rs=1,0 -vnsp=2 -vtpd1=4 -vbigbas=3 -vehmax=-.4 -vrwa=0.3359 -vemaxs=2.5 -vconvc=1e-6 -vbeta=.3 --ldadc=82.2 -job=1 --gprt:mode=3:nom=20
TMPDMFT
 lmfdmft ni -vnk=12 --rs=1,0 -vnsp=2 -vtpd1=4 -vbigbas=3 -vehmax=-.4 -vrwa=0.3359 -vemaxs=2.5 -vconvc=1e-6 -vbeta=.3 --ldadc=82.2 -job=1 -vindmfl=f
 mpirun -n 6 lmfdmft ni -vnk=12 --rs=1,0 -vnsp=2 -vtpd1=4 -vbigbas=3 -vehmax=-.4 -vrwa=0.3359 -vemaxs=2.5 -vconvc=1e-6 -vbeta=.3 --ldadc=82.2 -job=1 -vindmfl=f
 mpirun -n 6 lmfdmft ni -vnk=12 --rs=1,0 -vnsp=2 -vtpd1=4 -vbigbas=3 -vehmax=-.4 -vrwa=0.3359 -vemaxs=2.5 -vconvc=1e-6 -vbeta=.3 --ldadc=82.2 -job=1 --gprt:mode=3:nom=20
TSTDDOS
 lmf ni -vnk=24 -vnsp=2 -vtpd1=4 -vbigbas=3 -vehmax=-.4 -vrwa=0.3359 -vemaxs=2.5 -vconvc=1e-6 -vbeta=.3 -job=1 --dos:npts=201:ef0:rdm:window=-.2,.2 --quit=dos
 cp dos.ni doslmf.ni
 lmfdmft ni -vnk=24 -vnsp=2 -vtpd1=4 -vbigbas=3 -vehmax=-.4 -vrwa=0.3359 -vemaxs=2.5 -vconvc=1e-6 -vbeta=.3 --ldadc=82.2 -job=1
 lmfdmft ni -vnk=24 -vnsp=2 -vtpd1=4 -vbigbas=3 -vehmax=-.4 -vrwa=0.3359 -vemaxs=2.5 -vconvc=1e-6 -vbeta=.3 --ldadc=82.2 -job=1 --dos:npts=201:ef0:rdm:window=-.2,.2
 cp dos.ni doslmfdmft.ni
TMPDOS
 lmf ni -vnk=24 -vnsp=2 -vtpd1=4 -vbigbas=3 -vehmax=-.4 -vrwa=0.3359 -vemaxs=2.5 -vconvc=1e-6 -vbeta=.3 -job=1 --dos:npts=201:ef0:rdm:window=-.2,.2 --quit=dos
 cp dos.ni doslmf.ni
 lmfdmft ni -vnk=24 -vnsp=2 -vtpd1=4 -vbigbas=3 -vehmax=-.4 -vrwa=0.3359 -vemaxs=2.5 -vconvc=1e-6 -vbeta=.3 --ldadc=82.2 -job=1
 mpirun -n 6 lmfdmft ni -vnk=24 -vnsp=2 -vtpd1=4 -vbigbas=3 -vehmax=-.4 -vrwa=0.3359 -vemaxs=2.5 -vconvc=1e-6 -vbeta=.3 --ldadc=82.2 -job=1 --dos:npts=201:ef0:rdm:window=-.2,.2
 cp dos.ni doslmfdmft.ni
# fplot dos.ni -ord x2+x3 dmft/test/ni/doslmf.ni
# SYMGRP i*r3(1,1,-1) r4x SOC=t
% const asa=f lmf=t nit=20 so=0 nsp=so?2:1 soscl=1
% const kt=0 # If nonzero, use finite temperature, kT, kT in eV.
# --- Structural and site information ---
# Structure-specific data must be entered in the site file.
# For a particular compound,  you must:
#   1. Create a site file with lattice constant, lattice vectors, basis vectors
#   2. Enter below number of species and labels elt1,elt2,... for each element
#   3. Tailor the const rwa (MT radius in units of the lattice constant)
# ... elemental fcc: 1 species, rwa set for -3% overlap
# % const nspec=1 rwa=sqrt(2)/4*.97
# This for Pt,Co,etc
# % char0 elt1=Al
# % char0 elt1=Co
# % char0 elt1=Cr
% char0 elt1=Ni
# % char0 elt1=Pd
# % char0 elt1=Pt
# % char0 elt1=Cu
# % char0 elt1=Ag
# % char0 elt1=Cd
# % char0 elt1=Au
# % char0 elt1=Ce
# % char0 elt1=Er
# % char0 elt1=U
# % char0 elt1=K
# % char0 elt1=Fe
# ... elemental bcc: 1 species, rwa set for -3% overlap
#% const nspec=1 rwa=sqrt(3)/4*.97
# ... elemental fcc: 1 species, rwa set for -3% overlap
% const nspec=1 rwa=sqrt(2)/4*.97
# ... elemental hcp: 1 species, rwa set for -3% overlap
# % const nspec=1 rwa=1/2*.97
# ... for tetragonal binary compounds (a,0,0) (0,a,0) (0,0,c)
# % const nspec=2 rwa=2/4*.97
# % char elt1=Co elt2=Pt
# % char elt1=Fe elt2=Pt
# ... for bcc based binary compounds
# % const nspec=2 rwa=sqrt(3)/4*.97
# % char elt1=Ni elt2=Al
# % char elt1=Cr elt2=Cr2
#
# --- Macros and default values ---
# ... Macro pval(tp,default-p,low-p,high-p) select which P valence
#     Macro ploc(tp,default-p,low-p,high-p) select which PZ local
# macro   makes  tp=0    tp=1    tp=2    tp=3    tp=4
#  pval     P     p0      plo     phi     phi     plo
#  ploc     PZ     0      0       plo  10+plo     phi
% macro pval(tp,p0,plo,phi) tp==0?p0:((tp==1|tp==4)?plo:((tp==2|tp==3)?phi:0))
% macro ploc(tp,p0,plo,phi) tp<2?0:(tp==2?plo:(tp==3?10+plo:(tp==4?phi:0)))
% macro pqnsc(z) z<37?3:(z<55?4:5)
% macro pqnd(z)  z<37?3:(z<55?4:5)
# ... default switches for local orbitals: no orbitals
% repeat i=1,nspec
%  const tps{i}=0 tpp{i}=0 tpd{i}=0
% end
# ... conventions for lmf basis
#             bigbas basis
#               0    spd
#               1    spd+spd
#               2    spd+spdf
#               3    spd+spdfg
# ... default value for ehmax
% ifndef ehmax<0
%   ifdef ehmax | ehmax==0
%     echo (warning) illegal ehmax = {ehmax} ... resetting to -.2
%   endif
%   var ehmax=-.2
% endif
# ... other default parameters
% const hf=f xcn=0 lmh=f nl={lmf?5:{asa?4:5}} ef0=0 da=0 tet=0 trig=0
# ASA specific
% const ccor=t twoc=f gamma=twoc adnf=0 gfmod=0
# LMF-specific
% const bigbas=2 gmax=9 ngd=15 elind=-.7 lmxa=4 rsma=0 kmx=4 lfoca=1 met=lmf?5:1
% const cbya=sqrt(8/3)
#
# ... Sample Levenberg-Marquardt fitting. Generate qpts and reference bands with:
# lmf -vnk=20 --rs=0,0 -vtet=0 --rs=0,0 --putqp --quit=ham -vsig=0 fe `cat switches-for-lm`
# lmf -vsig=0 fe `cat switches-for-lm` --band:qp:fn=qpts
# ASA-LDA; use lmfit=1 for Levenberg-Marquardt fit
# lm -vnsp=2 -vnk=8 --rs=0,0 -vsclwsr=1 -vlmf=0 -vasa=1 -vnit=0 -vnl=4 fe -vrdvext=0 -vlmfit=0
# lm -vnsp=2 -vnk=8 --rs=0,0 -vsclwsr=1 -vlmf=0 -vasa=1 -vnit=20 -vnl=4 fe -vrdvext=1 -vlmfit=1 --iactiv
VERS    LMF-6 LMASA-6 LM:7 FP:7 ASA:7
% const pwmode=0 pwemin=0 pwemax=3 oveps=0 lxcf=2 lxcg=0 pfloat=0 rdvext=0 lmfit=0 bf=0
HAM     GMAX={gmax} FTMESH={ngd} FORCES={so==0&hf==f} ELIND={elind} SHFAC=0 # BXC0=t
        NSPIN={nsp} REL=t SO={so} RDVEXT={rdvext} BFIELD={bf}
        AUTOBAS[PFLOAT=0,{pfloat} MTO=0] PWMODE={pwmode} PWEMIN={pwemin} PWEMAX={pwemax} OVEPS={oveps}
%iffile specialh
%echo including file specialh in HAM category
%include specialh
%endif
%ifdef  lxcf & lxcf2+1
        XCFUN=0,{lxcf},{lxcf2}
%endif
        XCFUN={lxcf} GGA={lxcg}
# the following line for GW
% const sig=12 modsgp=4 nmins=0 emins=0 nmaxs=0 emaxs=2 asig=0.0 bsig=0.09 efits=0
# for v7
        SIGP[MODE={modsgp} NMIN={nmins} EMIN={emins} NMAX={nmaxs} EMAX={emaxs} A={asig} B={bsig} EFIT={efits}]
        RDSIG={sig} SIGP:{modsgp},{nmins},{emins},{nmaxs},{emaxs},{asig},{bsig},{efits}
% const modsig=3 ecuts=2 pwb=2.7 pwc=2.2
        RSRNGE=8
GW      NKABC=nkgw nkgw nkgw2 GCUTB={pwb} GCUTX={pwc} MKSIG={modsig} ECUTS={ecuts}
IO      SHOW=f HELP=F VERBOS=31 20 WKP=F IACTIV={nit==0}
% cchar strn twoc==0 p3
GF      MODE={gfmod} GFOPTS={strn}
STR     RMAX=3.5
CONST   rc=0 nk=16 nk2=nk bzj=0 cbya={cbya} nkgw=8 nkgw2=nk2
HEADER  Master file for multinary TM compounds, suitable for GW
% includo atparms
STRUC   NL={nl} NSPEC={nspec} DALAT={da}
%ifdef file==2
        FILE=site{file}
%elseifd file
        FILE=site
%endif
        FILE=site  # Strux data always read from a SITE file
# to rotate the symmetry-line file, put rotation matrix into file r, and then:
#  cat syml.ext | grep -v ^# | awk '{print $2,$3,$4}' | mcx r . -t -x -t  > 1
#  cat syml.ext | grep -v ^# | awk '{print $5,$6,$7}' | mcx r . -t -x -t  > 2
#  cat syml.ext | awk '{print $1}' | mcx -ff5.0,3f12.7,1x,3f12.7 . 1 2 -ccat -ccat | grep -v rows | sed s/\\.//
%ifdef rot==111     # Rotates (1,1,1) to z, i.e. bcc P1+P2+P3 = (1,1,1) -> (0,0,1)
        ROT=z:pi/4,y:acos(1/sqrt(3))
%elseifd rot==110   # Rotates (1,-1,0) to z, i.e. bcc P2-P1 = (1,-1,0) -> (0,0,1)
        ROT=z:-pi/4,y:pi/2
%elseifd rot==001  # Rotates (1,-1,0) to z, i.e. bcc P2-P1 = (1,-1,0) -> (0,0,1)
        ROT=y:pi/2
%endif
% ifdef tet
        SHEAR=0 0 1 {tet}
SYMGRP  r4z i r2(1,1,0)
% elseifd trig
        SHEAR=1 1 1 {trig}
SYMGRP  i*r3d r2(1,0,-1)
% endif
SYMGRP  find SOC=2
%ifdef file==2
SITE    FILE=site{file}
# %elseifd file
# SITE    FILE=site
%endif
SITE    FILE=site  # Strux data always read from a SITE file
SPEC
        HFAC[L={soscl}]
% ifdef sclwsr
        SCLWSR=1 OMAX1={asa?.16:.04} WSRMAX=3.30
% endif
% repeat i=1:nspec
        ATOM={elt{i}} Z={Z{i}} R/A={rwa} LMXA={nl-1} EREF={eref{i}} A={lmh?.015:.025}
% iffile specialspec{i}
%   echo including file specialspec{i} in SPEC category
%   include specialspec{i}
% endif
% ifdef mom{i}
        MMOM=0 0 {mom{i}}
% endif
        LMXA={lmxa} RSMA={rsma} KMXA={kmx} LFOCA={lfoca}
        P=0,0,0,{{Z{i}}>71?5.15:4.15},5.12 IDMOD=0,0,0,1,1
# build up PZ from tps{i}, tpp{i}, tpd{i}
# % trace 4
% char strn={ploc({tps{i}},0,{{pqnsc(Z{i})}+.94},{{pqnsc(Z{i})}+2.5})}
% char strn="{strn},{ploc({tpp{i}},0,{{pqnsc(Z{i})}+.93},{{pqnsc(Z{i})}+2.5})}
% char strn="{strn},{ploc({tpd{i}},0,0,{{pqnd(Z{i})}+1.5})}
        PZ={strn} RS3=0.95
% char strn="{rsm{i}},{rsm{i}},{rsmd{i}},{bigbas>=2?rsm{i}:0},{bigbas>=3?rsm{i}:0}"
        RSMH={strn} EH={ehmax},{ehmax},{ehmax},{ehmax},{ehmax}
% ifdef bigbas
        RSMH2={rsm{i}},{tpp{i}==2?1.0:rsm{i}},{rsmd{i}} EH2={ehmax-.8},{tpp{i}==2?-2:ehmax-.8},{ehmax-.8}
% endif
% end

OPTIONS NSPIN={nsp} REL=t XCFUN={lxcf} XCN={xcn} HF={hf} SO={so}
        Q =BAND INVIT=T GRCOR={lmh} LMH={lmh} RQUAD=0
        ASA[ CCOR={ccor} ADNF={adnf} TWOC={twoc} GAMMA={gamma} ]
        SAVVEC={gfmod>=10}
% ifdef kT  # If finite temperature
% const MPN=-1 MPW=kT/13.606 tetra=0
% endif
% const MPN=1 MPW=.002 tetra=1
BZ      NKABC=nk nk nk2 N={MPN} W={MPW} NPTS=1001 NKABC2=6 JOB2=bzj SAVDOS=F
%ifdef fsmom
        FSMOM={fsmom}
%endif
%if hcp
        BZJOB=0 0 bzj
%else
        BZJOB=bzj
%endif
        EF0={ef0} DELEF=.1 TETRA={tetra} DOS={ef0}-1 {ef0}+.5 METAL={met}
% ifdef hf
        NEVMX=-1
% endif
%ifdef gfmod>=10 | mdos
        NEVMX=9999 EFMAX=10
%endif
        EMESH=12 10 -1 0 .5 .5 INVIT=f NOINV=T
EWALD   AS=2.0 TOL=1D-8 ALAT0=a NKRMX=600 NKDMX=600
% const beta=.5 convc=.005/50
# ITER    MIX=A0,b={beta} NIT={nit} CONV=.00001 CONVC={convc} UMIX={beta/10}
ITER    MIX=B3,b={beta} NIT={nit} CONV=.00001 CONVC={convc} UMIX={beta/3}
ITER    MIX=A3,b={beta} NIT={nit} CONV=.00001 CONVC={convc}
# For Levenberg-Marquardt fittting
#       FIT[MODE={lmfit} NBFIT=1 8 NBFITF=1 LAM=1 SCL=5 SHFT=1 WT=1,.1,.1]
MIX     MODE=A3,b={beta}
        NMIX=2 AMIX=T BETA=.8 CONV=.00001 CONVC={convc}
START   CNTROL={nit==0} BEGMOM={nit==0} NIT={nit} FREE=f
DMFT    PROJ=2  NLOHI=1,8  BETA=50  NOMEGA=1999  KNORM=1
% const indmfl=f
%ifndef indmfl
        BLOCK:
        SITES=1 L= 2
        SIDXD= 1 2 3 4 5
        UMODE= full~static
%endif
