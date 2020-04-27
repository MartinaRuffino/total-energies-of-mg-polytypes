% const lmh=f nsp=2 nl=3 ccor=f nit=1 nk=16 novmza=f gfmod=10
% const twoc=t levec=f mode=twoc?10:0 tetra=t bzj=0 cond=0 qss=0 theta=0
HEADER  Study of Ni
VERS    LMASA-6 LM:7 ASA:7
IO      SHOW=F HELP=f VERBOS=31 20 WKP=F IACTIV=f TIM=f
HAM     NSPIN={nsp} REL={rel} NONCOL=f SO=T QASA=0
OPTIONS NSPIN={nsp} LMH={lmh} REL={rel} SO=T INVIT=F
        ASA[ ADNF=T CCOR={ccor} TWOC=F GAMMA={gamma} ]
        NONCOL=f RQUAD=2
%ifdef sym
SYMGRP  I R4Z
%endif
SYMGRP  E
BZ      NKABC={nk} SAVDOS=F NPTS=501 NEVMX=18 EFMAX=9 PUTQP=f
        TETRA={tetra} COND:{cond},0,0,1 BZJOB={bzj} DOS=-1 .5
% const  ef=-0.051492
        EMESH=12 10 {ef-0.9} {ef} .5 0 INVIT=f
GF      MODE=1 GFOPTS=emom;p3;padtol=1d-5
STR     RMAX=3.5 MODE={mode}
% const idmod=0
% ifdef file
STRUC   FILE=site NSPEC=1
SITE    FILE=site
SPEC    ATOM=A Z=28 R/W=1 IDXDN=1 1 1 1 IDMOD={idmod} {idmod} {idmod} {idmod} MMOM=0,0,1 A=0.01
% endif
SPEC    ATOM=A Z=28 LMX=2 R/W=1 IDXDN=1 1 1 1 IDMOD={idmod} {idmod} {idmod} {idmod} MMOM=0,0,1 A=0.01
        ATOM=B Z=28 LMX=2 R/W=1 IDXDN=1 1 1 1 IDMOD={idmod} {idmod} {idmod} {idmod} MMOM=0,0,1 A=0.01 # artificial
STRUC   NBAS=1 NSPEC=2 NL={nl}
        ALAT=6.621 PLAT=  0 .5 .5  .5 0 .5  .5 .5 0
#       ALAT=6.721 PLAT=  0 .5 .5  .5 0 .5  .5 .5 0
#       SLAT=  0 .5*2 .5*2  .5*2 0 .5*2  .5*2 .5*2 0
SITE    ATOM=A POS= 0 0 0
% const beta=1
ITER    MIX=B,b={beta} NIT={nit} CONVC=1D-5 CONV=0
MIX     MODE=B
% const begmom=1
START   CNTROL={nit==0} BEGMOM={begmom}
#       Float P to band CG
        ATOM=A        P=  4.6731692  4.4099332  3.8738842
                          4.6749678  4.4204305  3.8385854
                      Q=  0.3275075  0.0000000  0.0063744
                          0.3683053  0.0000000  0.0058928
                          4.6164026  0.0000000  0.0423487
                          0.3328032  0.0000000  0.0061967
                          0.3983901  0.0000000  0.0067149
                          3.9565913  0.0000000  0.0357336
# Self-consistent moments for IDMOD=4
% ifdef idmod==4
        ATOM=A        P=  4.6753021  4.4156763  3.8912489
                          4.6728490  4.4146378  3.8175530
                      Q=  0.3273776 -0.0012584  0.0063719
                          0.3686497 -0.0042616  0.0059484
                          4.6148944 -0.0552343  0.0430414
                          0.3327559  0.0012198  0.0061856
                          0.3976535  0.0045079  0.0067385
                          3.9586690  0.0474143  0.0362463
% endif
