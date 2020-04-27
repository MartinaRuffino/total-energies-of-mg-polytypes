% const lmh=f nsp=2 nl=3 ccor=f nit=1 nk=8 novmza=f gfmod=10
% const twoc=t levec=f mode=twoc?10:0 tetra=t bzj=0 cond=0 qss=0 theta=0
HEADER  Study of Ni
VERS    LMASA-6 LM:7 ASA:7
IO      SHOW=F HELP=f VERBOS=31 20 WKP=F IACTIV=t TIM=f
HAM     NSPIN={nsp} NONCOL=f QASA=0
OPTIONS NSPIN={nsp} LMH={lmh} INVIT=F NONCOL=f RQUAD=0
        ASA[ ADNF=T CCOR={ccor} TWOC={twoc} GAMMA=t ]
SYMGRP  R4X MX R3D
BZ      NKABC={nk} SAVDOS=F NPTS=501 NEVMX=18 EFMAX=9 PUTQP=f
        TETRA={tetra} COND:{cond},0,0,1 BZJOB={bzj} DOS=-1 .5
        EMESH=12 10 -1 0 .5 .5 INVIT=f
GF      MODE={gfmod} GFOPTS={?~twoc==0~p3;~}{?~levec==1~evec;~}
STR     RMAX=3.5 MODE={mode}
% ifdef file
STRUC   FILE=site NSPEC=1
SITE    FILE=site
SPEC    ATOM=A Z=28 R/W=1 IDXDN=0 0 0 2 MMOM=0,0,1 A=0.03
% endif

SPEC    ATOM=A Z=28 R/W=1 IDXDN=0 0 0 2 MMOM=0,0,1 A=0.03
STRUC   NBAS=1 NSPEC=1 NL={nl}
        ALAT=6.621 PLAT=  0 .5 .5  .5 0 .5  .5 .5 0
#       SLAT=  0 .5*2 .5*2  .5*2 0 .5*2  .5*2 .5*2 0
SITE    ATOM=A POS= 0 0 0
ITER    MIX=B NIT={nit} CONVC=1D-7 CONV=0
MIX     MODE=B
START   CNTROL={nit==0} BEGMOM={nit==0} NIT={nit} CNVG=1D-7
