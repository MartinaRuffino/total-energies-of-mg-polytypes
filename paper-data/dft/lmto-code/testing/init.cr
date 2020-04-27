# Init file for Cr (Wein2K equilibrium volume 11.773A, http://molmod.ugent.be/sites/default/files/deltadftcodes/supplmat/SupplMat-WIEN2k.pdf)
LATTICE
    ALAT=5.41632 PLAT= 1 0 0   0 1 0   0 0 1
SPEC
    ATOM=Cr   MMOM=0,0,1  TOKEN="IDMOD=0,\{idp}"
    ATOM=Cr2  MMOM=0,0,-1 TOKEN=IDMOD=0,\{idp}
SITE
    ATOM=Cr   X=0 0 0
    ATOM=Cr2  X=.5 .5 .5
BZ  FSMOM=1d-6
HAM REL=11