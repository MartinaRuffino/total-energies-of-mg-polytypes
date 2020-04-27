% const Harrison=1
% const io=0 nitq=30 qtol=1d-6 nx=5 beta=0.75
% const nk=4 diamond=0 azulene=0 graphene=0
% const verb=31 so=0 nsp=1 tetra=0 metal=1 ovlp=0 ul={azulene?1:0}
% const au=0.529177 Ry=13.61 Nd=6
% const dyn=0 relax=0 temp=300 taup=10 taub=100 time=100000 tstep=0.5
% const hess=1 gtol=1d-4 xtol=1d-4 step=0.1 nitf=100 nkill=100 relax=0
% const fs=0.048377 K=1/0.6333328d-5 amass=1.09716d-3
HEADER  Graphene, Azulene and Diamond
VERS    TB=9 LM=7
IO      SHOW=F HELP=F VERBOS={verb} WKP=F
CONST   ftmesh=24 nit=100 nmix=3 beta=1 conv=1d-4 convc=1d-4
# Use -vsi=1 for Si, rather than C
%  ifdef si==1
        d0d=2.35/{au} a0d=d0d*4/sqrt(3) V0=a0d^3 vfrac=1 V=V0*vfrac
%  else
        d0d=1.54/{au} a0d=d0d*4/sqrt(3) V0=a0d^3 vfrac=1 V=V0*vfrac
%  endif
        ad=v^(1/3)
        a0g=2.68525 ag=a0g q=8 r3=sqrt(3)
        nka={nk} nkb={nk} nkc={diamond?nk:1} mull=-1
        sss=-5/{Ry} sps=4.7/{Ry} pps=5.5/{Ry} ppp=-1.55/{Ry}
        n=2 nc=6.5 r0=1.5363/{au} rc=2.18/{au}
        q0sCnu=2 q0pCnu=2  q0sCnd=0 q0pCnd=0 q0Hu=1 q0Hd=0
        UC=1.2 UH=1.2 JC=0.0 JH=0.0
%  if Harrison==1
        esC=-19.37/{Ry} epC=-11.07/{Ry} esH=-1
%  else
        esC=-0.220 epC=0.273 esH=-0.349
%  endif
        esC=-0.220 epC=0.273 esH=-0.349
        CHsss=-0.4796324  CHsps=0.5008088
        CHn=0.5663000  CHnc=3.1955000  CHr0=2.0491493 CHrc=2.2686389
        CCA=8.18555/{Ry} CCnp=3.303 CCncp=8.6655
        CCr0p=1.64/{au} CCrcp=2.1052/{au} p=10
        CHA=11.4813/{Ry} CHnp=1.408 CHncp=3.5077
        CHr0p=1.084/{au} CHrcp=1.5474/{au} p=10
        rmaxh={diamond?0.5:(azulene?3:1.01)} mol={azulene?1:0}
STRUC   NBAS={azulene?18:2} NSPEC=2
        NL=2
% ifdef diamond
        ALAT=ad PLAT=0 1/2 1/2
                     1/2 0 1/2
                     1/2 1/2 0
% elseifd graphene
        ALAT=ag PLAT= 3/2  r3/2  0
                     -3/2  r3/2  0
                      0     0   q
% elseifd azulene
                ALAT=1 PLAT=p 0 0 0 p 0  0 0 p
% endif
SPEC    ATOM=C Z=6 R/W=1 IDXDN=1 1 COLOUR=0.8 0.8 0.8 RDAIUS=0.5 AMASS=12/{amass}
        ATOM=H Z=1 R/W=1 IDXDN=1 3 COLOUR=0.8 0.1 0.1 RADIUS=0.2 AMASS=1.00794/{amass}
SITE
% ifdef diamond
        ATOM=C POS=0 0 0
        ATOM=C POS=1/4 1/4 1/4
% elseifd graphene
        ATOM=C POS= 1/2  r3/2  0
        ATOM=C POS=-1/2  r3/2  0
% elseifd azulene
        ATOM=C POS= -1.382692 -0.015105 0
        ATOM=C POS=  1.378829 -0.012768 0
        ATOM=C POS=  2.922391 -2.041688 0
        ATOM=C POS=  2.317236 -4.528024 0
        ATOM=C POS= -0.001756 -5.606962 0
        ATOM=C POS= -2.316937 -4.526679 0
        ATOM=C POS= -2.921573 -2.044488 0
        ATOM=C POS=  2.098442  2.453647 0
        ATOM=C POS= -0.004822  3.922413 0
        ATOM=C POS= -2.107111  2.450709 0
        ATOM=H POS=  4.908613 -1.614889 0
        ATOM=H POS=  3.881018 -5.815217 0
        ATOM=H POS= -0.002395 -7.634188 0
        ATOM=H POS= -3.880021 -5.812731 0
        ATOM=H POS= -4.907183 -1.619900 0
        ATOM=H POS=  4.004513  3.096787 0
        ATOM=H POS= -0.004927  5.939342 0
        ATOM=H POS= -4.013888  3.090810 0
% endif
OPTIONS NSPIN={nsp} REL=T SO={so} XCFUN=2
SYMGRP  find
BZ      NKABC=nka nkb nkc TETRA={tetra} METAL={metal}
        EF0=0.95 DELEF=0.1 N=1 W=0.02
        NPTS=1001 BZJOB=1 SAVDOS=1 EFMAX=5 NEVMX=36 NOINV=T
        INVIT=F MULL=mull DOS=-0.5 1
STR     MXNBR=200
ME
        5
       C C MEMODE=5 PPMODE=30 POLY=5 CUTMOD=0 CUTPP=0 0
            | sss n nc r0 rc
              sps n nc r0 rc
              pps n nc r0 rc
              ppp n nc r0 rc
            ! CCA 1 -1 CCnp CCncp CCr0p CCrcp 0 0
       C H  MEMODE=5 PPMODE=30 POLY=5 CUTMOD=0 CUTPP=0 0
            | CHsss CHn CHnc CHr0 CHrc
              CHsps CHn CHnc CHr0 CHrc
              0 0 0 0 0
              0 0 0 0 0
            ! CHA 1 -1 CHnp CHncp CHr0p CHrcp 0 0
       H C  MEMODE=5 PPMODE=30 POLY=5 CUTMOD=0 CUTPP=0 0
            | CHsss CHn CHnc CHr0 CHrc
              CHsps CHn CHnc CHr0 CHrc
              0 0 0 0 0
              0 0 0 0 0
            ! CHA 1 -1 CHnp CHncp CHr0p CHrcp 0 0
       H H  MEMODE=5 PPMODE=0 POLY=5 CUTMOD=0 CUTPP=0 0
            | 0 0 0 0 0
              0 0 0 0 0
              0 0 0 0 0
              0 0 0 0 0
            ! 0 0 0 0 0 0 0 0 0
TB      FORCES=1 EVDISC=T RMAXH=rmaxh TRH=F RHO=F
        UL={UL} IODEL={io} OVLP={ovlp} 3PV=T MOL=mol
EWALD   TOL=1D-8 NKDMX=1999 NKRMX=1999
HAM     FTMESH=ftmesh ftmesh q*ftmesh
ITER    CONV=conv CONVC=qtol NIT={nitq} MIX=A{nx},b={beta}
DYN
% if dyn==1|dyn==2|dyn==3
        MD[MODE={dyn} TSTEP={tstep/fs} TEMP={temp/K} TAUP={taup/fs}
           TIME={time/fs} TAUB={taub/fs}]
% elseif relax>0
        MSTAT[MODE={relax} HESS={hess} XTOL={xtol} GTOL={gtol}
              STEP={step} NKILL={nkill}] NIT={nitf}
% endif
START   CNTROL=T
        ATOM=C P=2 2   2 2
               Q=q0sCnu   esC   UC
                 q0pCnu   epC   UC
                 q0sCnd   esC   UC
                 q0pCnd   epC   UC
        ATOM=H P=1 1   1 1
               Q=q0Hu   esH   UH
                    0    0     0
                 q0Hd   esH   UH
                    0   0      0
