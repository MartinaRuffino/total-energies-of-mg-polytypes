%const  nit=30 begmom=nit<=0
% const lmh=f nsp=2 nl=3 ccor=f gamma=0 nit=1 nk=16 novmza=f gfmod=1 so=t
% const twoc=t levec=f mode=twoc?10:0 tetra=t bzj=0 cond=0 qss=0 theta=0
HEADER  Study of Ni
VERS    LMASA-6 LM:7 ASA:7
IO      SHOW=F HELP=f VERBOS=31 20 WKP=F IACTIV=f TIM=f
HAM     NSPIN={nsp} REL={rel} NONCOL=f SO={so} QASA=0
% const sharm=1
OPTIONS NSPIN={nsp} LMH={lmh} REL={rel} SO={so} INVIT=F 
        ASA[ ADNF=T CCOR={ccor} TWOC={twoc} GAMMA={gamma} ]
        NONCOL=f RQUAD=2
	SHARM={sharm} 
%ifdef sym
SYMGRP  I R4Z MX SOC=1
%endif
SYMGRP  E
BZ      NKABC={nk} SAVDOS=F NPTS=501 NEVMX=18 EFMAX=9 PUTQP=f
        TETRA={tetra} COND:{cond},0,0,1 BZJOB={bzj} DOS=-1 .5
% const  ef=-0.055806
        EMESH=12 10 {ef-0.9} {ef} .5 0 INVIT=f
GF      MODE={gfmod} GFOPTS=emom;padtol=3e-6  # ;shftef
STR     RMAX=3.5 MODE={mode}
% ifdef file
STRUC   FILE=site NSPEC=1
SITE    FILE=site
% endif
% const idmod=0
SPEC    ATOM=A Z=28 R/W=1 IDXDN=1 1 1 1 IDMOD={idmod} {idmod} {idmod} {idmod} MMOM=0,0,1 A=0.015 # NR=1501
STRUC   NBAS=1 NSPEC=1 NL={nl}
        ALAT=6.621 PLAT=  0 .5 .5  .5 0 .5  .5 .5 0
#	ALAT=6.721 PLAT=  0 .5 .5  .5 0 .5  .5 .5 0
#       SLAT=  0 .5*2 .5*2  .5*2 0 .5*2  .5*2 .5*2 0
SITE    ATOM=A POS= 0 0 0
% const beta=1
ITER    MIX=A3,k=8 NIT={abs(nit)} CONVC=3D-6 CONV=0
MIX     MODE=B
START   CNTROL={nit==0} BEGMOM={begmom}
#       Self-consistent with SO coupling included
        ATOM=A        P=  4.6725998  4.4175716  3.8417047
                          4.6702928  4.4076540  3.8745404
                      Q=  0.3443967  0.0000000  0.0061619
                          0.4007489  0.0000000  0.0067936
                          3.9587553  0.0000000  0.0362773
                          0.3395225  0.0000000  0.0062606
                          0.3733664  0.0000000  0.0060445
                          4.5832097  0.0000000  0.0430491
# Self-consistent moments for IDMOD=4
# lmgf -vsym=t -vrel=1 -vso=f -vnit=1 -vidmod=4 -vgamma=0 ni --iactiv --pr41,41
% ifdef idmod==4
        ATOM=A        P=  4.6701022  4.4117171  3.8212065
                          4.6731941  4.4134306  3.8911086
                      Q=  0.3425507  0.0015420  0.0061260
                          0.3984117  0.0046143  0.0067695
                          3.9601512  0.0470842  0.0364677
                          0.3411120 -0.0019677  0.0062659
                          0.3750908 -0.0045731  0.0061301
                          4.5826837 -0.0540647  0.0436888
% endif
