% const q=1.626952235 vol=155.316320
% const alat=({4*vol}/({q}*sqrt(3)))^(1/3)

HEADER  HCP Mg
VERS    LM:7 FP:7
BZ      NKABC=30 30 30 METAL=5 TETRA=0 N=-1 W=0.03
HAM     AUTOBAS[LMTO=4 MTO=4 PNU=1 LOC=1 GW=0]
           PWMODE=0 FORCES=1 ELIND=-1.0 GGA=3 XCFUN=4 GMAX=8.2
IO      SHOW=T HELP=f WKP=F IACTIV=f VERBOS=31
ITER    MIX=B2,b=.1,k=10,wc=1 NIT=1000 CONVC=1d-6 CONV=1d-6
SYMGRP  find
SPEC    ATOM=Mg Z=12 R=3.012678
        LMX=2 LMXA=3
STRUC   NBAS=3 NL=3 NSPEC=1
        ALAT={alat}
        PLAT=    1.0    0.0	0.0
	     	-0.5 sqrt(3)/2  0.0
		 0.0    0.0  	3*{q}/2
SITE    ATOM=Mg POS=0.0 0.0 0.0
        ATOM=Mg POS=0.5 sqrt(3)/6 {q}/2
	ATOM=Mg POS=0.0 sqrt(3)/3 {q}
