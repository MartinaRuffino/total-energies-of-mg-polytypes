# make self-consistent with:
lmstr ni -vnit=0 -vnk=32 -vbzj=0 -vtwoc=t -vccor=f -vexc=10 -vni=t -vnit=0
# make evecs with
lm ni -vnit=0 -vnk=32 -vbzj=0 -vtwoc=t -vccor=f -vexc=10 -vni=t -vnit=0
# make jr with
lmgfe ni -vnit=0 -vnk=32 -vbzj=0 -vtwoc=t -vccor=f -vexc=10 -vni=t -vnit=0
# make resultant quantities
lmgfe ni -vnit=0 -vnk=32 -vbzj=0 -vtwoc=t -vccor=f -vexc=11 -vni=t -vnit=0
