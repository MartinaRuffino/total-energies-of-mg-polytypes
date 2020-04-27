rdcmd:  echo -1 | lmfgwd cu -vnk=8 -vbigbas=t -vpwmode=11 --make-Q0P
 -----------------------  START LMFGWD -----------------------

 LMFGWD:   nbas = 1  nspec = 1  vn 7.11.b  verb 31,30,31
 special:  APW basis(q)
 pot:      XC:BH
 auto:     Pfloat:v6(CG v6)  Autoread:  none
 bz:       metal(3), tetra, invit 

                Plat                                  Qlat
   0.000000   0.500000   0.500000       -1.000000   1.000000   1.000000
   0.500000   0.000000   0.500000        1.000000  -1.000000   1.000000
   0.500000   0.500000   0.000000        1.000000   1.000000  -1.000000
   alat = 6.798  Cell vol = 78.538660

 LATTC:  as= 2.000   tol=1.00E-12   alat= 6.79800   awald= 0.467
         r1=  2.252   nkd= 201      q1=  6.871   nkg= 331

 SGROUP: 1 symmetry operations from 0 generators
 SYMLAT: Bravais system is cubic with 48 symmetry operations.
 SYMCRY: crystal invariant under 48 symmetry operations for tol=1e-5
 GROUPG: the following are sufficient to generate the space group:
         i*r3(1,1,-1) r4x
         i*r3(1,1,-1) r4x
 MKSYM:  found 48 space group operations
         no attempt to add inversion symmetry

 species data:  augmentation                           density
 spec       rmt   rsma lmxa kmxa      lmxl     rg   rsmv  kmxv foca   rfoca
 A        2.311  0.925    4    4         4  0.578  1.156    15    1   0.925

 MSHSIZ: mesh has 11 x 11 x 11 divisions; length 0.437, 0.437, 0.437
         generated from gmax = 9.0 a.u. : 941 vectors of 1331 (70%)

 GVLIST: gmax = 9 a.u. created 941 vectors of 1331 (70%)
         mesh has 11 x 11 x 11 divisions; length 0.437, 0.437, 0.437
 SGVSYM: 41 symmetry stars found for 941 reciprocal lattice vectors
 

 Makidx:  hamiltonian dimensions Low, Int, High, Negl: 31 0 28 16
 kappa   Low   Int   High  L+I  L+I+H  Neglected
   1       9     0    16     9    25       0
   2      13     0    12    13    25       0
   3       9     0     0     9     9      16
  all     31     0    28    31    59      16
 suham :  25 augmentation channels, 25 local potential channels  Maximum lmxa=4

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol=1.0e-6
 spec      l    rsm    eh     gmax    last term   cutoff
  A        0*   2.50  -0.01   2.974    2.30E-05      27 
  A        1    2.50  -0.01   3.093    1.29E-06      51 
  A        2    1.00  -0.01   8.508    1.16E-06     869 
  A        0    1.30  -1.00   5.718    2.28E-06     259 
  A        2    1.00  -1.00   8.508    1.16E-06     869 
  A        3    1.30  -0.01   6.806    2.09E-06     411 
 

 suham:  q-dependent PW basis with  Emin = 1 < E < 3.
         Est. min,max PW dimension = 3,10.  Use npwpad = 3 => ndham = 44

 iors  : read restart file (binary, mesh density) 
         use from  restart file: ef window, positions, pnu 
         ignore in restart file: *

 --- BNDFP:  begin iteration 1 of 1 ---

 Make potential without XC part ...

 avg es pot at rmt= 0.582618  avg sphere pot= 0.651618  vconst=-0.582618

 smooth rhoves     10.174210   charge     3.709217

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000135

 Energy terms:             smooth           local           total
   rhoval*vef             -2.771213      -184.378311      -187.149524
   rhoval*ves            -46.687198      -116.549186      -163.236383
   psnuc*ves              67.035618    -12972.035039    -12904.999421
   utot                   10.174210     -6544.292112     -6534.117902
   rho*exc                 0.000000      -127.271001      -127.271001
   rho*vxc                 0.000000      -168.525834      -168.525834
   valence chg             3.709217         7.290783        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 avg es pot at rmt= 0.582618  avg sphere pot= 0.651618  vconst=-0.582618

 smooth rhoves     10.174210   charge     3.709217
 smooth rhoeps =   -2.662560   rhomu =   -3.464668  avg vxc =   -0.832977 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 writing sphere density to file rhoMT.1 ...
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000135

 Energy terms:             smooth           local           total
   rhoval*vef             -5.980485      -184.378311      -190.358796
   rhoval*ves            -46.687198      -116.549186      -163.236383
   psnuc*ves              67.035618    -12972.035039    -12904.999421
   utot                   10.174210     -6544.292112     -6534.117902
   rho*exc                -2.662560      -127.271001      -129.933562
   rho*vxc                -3.464668      -168.525834      -171.990502
   valence chg             3.709217         7.290783        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 Start first of two band passes ...
  lmfgw: input one of the following jobs:
   -1 : creates files GWinput, QPNT, QIBZ, Q0P, QGpsi, QGcou, KPTin1BZ
   -2 : Performs some sanity checks on GWinput
    0 : init mode; creates files SYMOPS, LATTC, CLASS, NLAindx, ldima
    1 : GW driver mode; creates files gwb,gw1,gw2,gwa,vxc,evec,rhoMT.*,normchk
    4 : band mode 
    5 : eigenvalue-only mode 
  job?
 ... Creating files GWinput, QIBZ, KPTin1BZ
 BZMESH:  29 irreducible QP from 512 ( 8 8 8 )  shift= F F F

 ... Creating in GWinput the following product basis:
 site  phi(o)     phi(u)     dot(o)     dot(u)     loc(o)     loc(u)
   1   11100      11110      00000      11100      000--      111--   

 ... Creating files Q0P, QGpsi, QGcou
 suq0x : 6 special gamma-points  auxf = 3.526133 (error = 0.0)
         qa = 0.065873  q0x(1) = -0.03803162 0.03803162 0.03803162
 q0irre: 1 irreducible special gamma-point(s)

 Maximum number of PW = 29 for basis,  18 for coulomb
 Extended mesh has 149 irreducible qp out of 1024

 Exit 0 lmfgw, job -1 
 CPU time:    0.271s   Wall clock    0.319s   08:49:47 30.06.2014        on phpdl1.ph.kcl.ac.uk
rdcmd:  echo  0 | lmfgwd cu -vnk=8 -vbigbas=t -vpwmode=11
 -----------------------  START LMFGWD -----------------------

 LMFGWD:   nbas = 1  nspec = 1  vn 7.11.b  verb 31,30,31
 special:  APW basis(q)
 pot:      XC:BH
 auto:     Pfloat:v6(CG v6)  Autoread:  none
 bz:       metal(3), tetra, invit 

                Plat                                  Qlat
   0.000000   0.500000   0.500000       -1.000000   1.000000   1.000000
   0.500000   0.000000   0.500000        1.000000  -1.000000   1.000000
   0.500000   0.500000   0.000000        1.000000   1.000000  -1.000000
   alat = 6.798  Cell vol = 78.538660

 LATTC:  as= 2.000   tol=1.00E-12   alat= 6.79800   awald= 0.467
         r1=  2.252   nkd= 201      q1=  6.871   nkg= 331

 SGROUP: 1 symmetry operations from 0 generators
 SYMLAT: Bravais system is cubic with 48 symmetry operations.
 SYMCRY: crystal invariant under 48 symmetry operations for tol=1e-5
 GROUPG: the following are sufficient to generate the space group:
         i*r3(1,1,-1) r4x
         i*r3(1,1,-1) r4x
 MKSYM:  found 48 space group operations
         no attempt to add inversion symmetry

 species data:  augmentation                           density
 spec       rmt   rsma lmxa kmxa      lmxl     rg   rsmv  kmxv foca   rfoca
 A        2.311  0.925    4    4         4  0.578  1.156    15    1   0.925

 MSHSIZ: mesh has 11 x 11 x 11 divisions; length 0.437, 0.437, 0.437
         generated from gmax = 9.0 a.u. : 941 vectors of 1331 (70%)

 GVLIST: gmax = 9 a.u. created 941 vectors of 1331 (70%)
         mesh has 11 x 11 x 11 divisions; length 0.437, 0.437, 0.437
 SGVSYM: 41 symmetry stars found for 941 reciprocal lattice vectors
 

 Makidx:  hamiltonian dimensions Low, Int, High, Negl: 31 0 28 16
 kappa   Low   Int   High  L+I  L+I+H  Neglected
   1       9     0    16     9    25       0
   2      13     0    12    13    25       0
   3       9     0     0     9     9      16
  all     31     0    28    31    59      16
 suham :  25 augmentation channels, 25 local potential channels  Maximum lmxa=4

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol=1.0e-6
 spec      l    rsm    eh     gmax    last term   cutoff
  A        0*   2.50  -0.01   2.974    2.30E-05      27 
  A        1    2.50  -0.01   3.093    1.29E-06      51 
  A        2    1.00  -0.01   8.508    1.16E-06     869 
  A        0    1.30  -1.00   5.718    2.28E-06     259 
  A        2    1.00  -1.00   8.508    1.16E-06     869 
  A        3    1.30  -0.01   6.806    2.09E-06     411 
 

 suham:  q-dependent PW basis with  Emin = 1 < E < 3.
         Est. min,max PW dimension = 3,10.  Use npwpad = 3 => ndham = 44

 iors  : read restart file (binary, mesh density) 
         use from  restart file: ef window, positions, pnu 
         ignore in restart file: *

 --- BNDFP:  begin iteration 1 of 1 ---

 Make potential without XC part ...

 avg es pot at rmt= 0.582618  avg sphere pot= 0.651618  vconst=-0.582618

 smooth rhoves     10.174210   charge     3.709217

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000135

 Energy terms:             smooth           local           total
   rhoval*vef             -2.771213      -184.378311      -187.149524
   rhoval*ves            -46.687198      -116.549186      -163.236383
   psnuc*ves              67.035618    -12972.035039    -12904.999421
   utot                   10.174210     -6544.292112     -6534.117902
   rho*exc                 0.000000      -127.271001      -127.271001
   rho*vxc                 0.000000      -168.525834      -168.525834
   valence chg             3.709217         7.290783        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 avg es pot at rmt= 0.582618  avg sphere pot= 0.651618  vconst=-0.582618

 smooth rhoves     10.174210   charge     3.709217
 smooth rhoeps =   -2.662560   rhomu =   -3.464668  avg vxc =   -0.832977 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 writing sphere density to file rhoMT.1 ...
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000135

 Energy terms:             smooth           local           total
   rhoval*vef             -5.980485      -184.378311      -190.358796
   rhoval*ves            -46.687198      -116.549186      -163.236383
   psnuc*ves              67.035618    -12972.035039    -12904.999421
   utot                   10.174210     -6544.292112     -6534.117902
   rho*exc                -2.662560      -127.271001      -129.933562
   rho*vxc                -3.464668      -168.525834      -171.990502
   valence chg             3.709217         7.290783        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 Start first of two band passes ...
  lmfgw: input one of the following jobs:
   -1 : creates files GWinput, QPNT, QIBZ, Q0P, QGpsi, QGcou, KPTin1BZ
   -2 : Performs some sanity checks on GWinput
    0 : init mode; creates files SYMOPS, LATTC, CLASS, NLAindx, ldima
    1 : GW driver mode; creates files gwb,gw1,gw2,gwa,vxc,evec,rhoMT.*,normchk
    4 : band mode 
    5 : eigenvalue-only mode 
  job?

 gw setup, job 0 code  0

 Creating files SYMOPS, LATTC, CLASS, NLAindx, ldima

 OK! lmfgw generated files, job 0
rdcmd:  echo  1 | lmfgwd cu -vnk=8 -vbigbas=t -vpwmode=11 --wrange
 -----------------------  START LMFGWD -----------------------

 LMFGWD:   nbas = 1  nspec = 1  vn 7.11.b  verb 31,30,31
 special:  APW basis(q)
 pot:      XC:BH
 auto:     Pfloat:v6(CG v6)  Autoread:  none
 bz:       metal(3), tetra, invit 

                Plat                                  Qlat
   0.000000   0.500000   0.500000       -1.000000   1.000000   1.000000
   0.500000   0.000000   0.500000        1.000000  -1.000000   1.000000
   0.500000   0.500000   0.000000        1.000000   1.000000  -1.000000
   alat = 6.798  Cell vol = 78.538660

 LATTC:  as= 2.000   tol=1.00E-12   alat= 6.79800   awald= 0.467
         r1=  2.252   nkd= 201      q1=  6.871   nkg= 331

 SGROUP: 1 symmetry operations from 0 generators
 SYMLAT: Bravais system is cubic with 48 symmetry operations.
 SYMCRY: crystal invariant under 48 symmetry operations for tol=1e-5
 GROUPG: the following are sufficient to generate the space group:
         i*r3(1,1,-1) r4x
         i*r3(1,1,-1) r4x
 MKSYM:  found 48 space group operations
         no attempt to add inversion symmetry

 species data:  augmentation                           density
 spec       rmt   rsma lmxa kmxa      lmxl     rg   rsmv  kmxv foca   rfoca
 A        2.311  0.925    4    4         4  0.578  1.156    15    1   0.925

 MSHSIZ: mesh has 11 x 11 x 11 divisions; length 0.437, 0.437, 0.437
         generated from gmax = 9.0 a.u. : 941 vectors of 1331 (70%)

 GVLIST: gmax = 9 a.u. created 941 vectors of 1331 (70%)
         mesh has 11 x 11 x 11 divisions; length 0.437, 0.437, 0.437
 SGVSYM: 41 symmetry stars found for 941 reciprocal lattice vectors
 

 Makidx:  hamiltonian dimensions Low, Int, High, Negl: 31 0 28 16
 kappa   Low   Int   High  L+I  L+I+H  Neglected
   1       9     0    16     9    25       0
   2      13     0    12    13    25       0
   3       9     0     0     9     9      16
  all     31     0    28    31    59      16
 suham :  25 augmentation channels, 25 local potential channels  Maximum lmxa=4

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol=1.0e-6
 spec      l    rsm    eh     gmax    last term   cutoff
  A        0*   2.50  -0.01   2.974    2.30E-05      27 
  A        1    2.50  -0.01   3.093    1.29E-06      51 
  A        2    1.00  -0.01   8.508    1.16E-06     869 
  A        0    1.30  -1.00   5.718    2.28E-06     259 
  A        2    1.00  -1.00   8.508    1.16E-06     869 
  A        3    1.30  -0.01   6.806    2.09E-06     411 
 

 suham:  q-dependent PW basis with  Emin = 1 < E < 3.
         Est. min,max PW dimension = 3,10.  Use npwpad = 3 => ndham = 44

 iors  : read restart file (binary, mesh density) 
         use from  restart file: ef window, positions, pnu 
         ignore in restart file: *

 --- BNDFP:  begin iteration 1 of 1 ---

 Make potential without XC part ...

 avg es pot at rmt= 0.582618  avg sphere pot= 0.651618  vconst=-0.582618

 smooth rhoves     10.174210   charge     3.709217

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000135

 Energy terms:             smooth           local           total
   rhoval*vef             -2.771213      -184.378311      -187.149524
   rhoval*ves            -46.687198      -116.549186      -163.236383
   psnuc*ves              67.035618    -12972.035039    -12904.999421
   utot                   10.174210     -6544.292112     -6534.117902
   rho*exc                 0.000000      -127.271001      -127.271001
   rho*vxc                 0.000000      -168.525834      -168.525834
   valence chg             3.709217         7.290783        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 avg es pot at rmt= 0.582618  avg sphere pot= 0.651618  vconst=-0.582618

 smooth rhoves     10.174210   charge     3.709217
 smooth rhoeps =   -2.662560   rhomu =   -3.464668  avg vxc =   -0.832977 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 writing sphere density to file rhoMT.1 ...
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000135

 Energy terms:             smooth           local           total
   rhoval*vef             -5.980485      -184.378311      -190.358796
   rhoval*ves            -46.687198      -116.549186      -163.236383
   psnuc*ves              67.035618    -12972.035039    -12904.999421
   utot                   10.174210     -6544.292112     -6534.117902
   rho*exc                -2.662560      -127.271001      -129.933562
   rho*vxc                -3.464668      -168.525834      -171.990502
   valence chg             3.709217         7.290783        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 Start first of two band passes ...
  lmfgw: input one of the following jobs:
   -1 : creates files GWinput, QPNT, QIBZ, Q0P, QGpsi, QGcou, KPTin1BZ
   -2 : Performs some sanity checks on GWinput
    0 : init mode; creates files SYMOPS, LATTC, CLASS, NLAindx, ldima
    1 : GW driver mode; creates files gwb,gw1,gw2,gwa,vxc,evec,rhoMT.*,normchk
    4 : band mode 
    5 : eigenvalue-only mode 
  job?

 gw setup, job 1 code  0

 ... Generate radial wave functions (file gwa)

 getcor:  qcore= 18.00  qsc= 18.00  konf = 4  4  3  4  5 
 sum q=18.00  sum ec= -1904.12818394  sum tc=  3173.29726688  rho(rmax) 0.00031

 ... Make w.f. and matrix elements, 1024 qp: file(s) gwb, vxc, evec

 sugw:  kpt 1 of 1024, k=  0.00000  0.00000  0.00000
 -0.7285 -0.2571 -0.2571 -0.2571 -0.1934 -0.1934  1.6532  1.8365  1.8365
 sugw:  kpt 11 of 1024, k=  0.37500  0.12500  -0.12500
 -0.5816 -0.2896 -0.2407 -0.2323 -0.2030 -0.1782  1.0731  1.4084  1.5119
 sugw:  kpt 21 of 1024, k=  0.75000  0.25000  -0.25000
 -0.3634 -0.3200 -0.2554 -0.2079 -0.1473  0.0397  0.4660  0.9012  1.1572
 sugw:  kpt 31 of 1024, k=  0.12500  -0.62500  0.62500
 -0.3453 -0.3056 -0.2682 -0.2152 -0.1647  0.1284  0.3587  0.7828  1.2762
 sugw:  kpt 41 of 1024, k=  -0.37500  0.37500  -0.37500
 -0.4453 -0.2539 -0.2539 -0.2259 -0.1628 -0.1628  0.4805  1.5735  1.6080
 sugw:  kpt 51 of 1024, k=  0.00000  0.50000  -0.50000
 -0.3831 -0.2952 -0.2790 -0.2183 -0.1943 -0.0519  0.6898  0.7768  1.4825
 sugw:  kpt 61 of 1024, k=  -0.62500  -0.37500  0.37500
 -0.3931 -0.2739 -0.2449 -0.2305 -0.1495 -0.0204  0.2987  1.2633  1.3915
 sugw:  kpt 71 of 1024, k=  -0.37500  -0.12500  0.37500
 -0.4922 -0.2866 -0.2532 -0.2352 -0.1983 -0.1439  0.7755  1.1719  1.7128
 sugw:  kpt 81 of 1024, k=  0.12500  -0.12500  0.37500
 -0.5816 -0.2896 -0.2407 -0.2323 -0.2030 -0.1782  1.0731  1.4084  1.5119
 sugw:  kpt 91 of 1024, k=  0.50000  0.00000  0.25000
 -0.4796 -0.3055 -0.2572 -0.2060 -0.2022 -0.1416  0.9564  1.0047  1.3968
 sugw:  kpt 101 of 1024, k=  0.87500  0.12500  0.12500
 -0.3835 -0.3516 -0.1882 -0.1669 -0.1442  0.0654  0.6193  0.7257  0.9928
 sugw:  kpt 111 of 1024, k=  -0.75000  0.25000  0.00000
 -0.3750 -0.3345 -0.2451 -0.1726 -0.1619 -0.0136  0.7779  0.8329  0.8619
 sugw:  kpt 121 of 1024, k=  -0.25000  0.25000  0.00000
 -0.6196 -0.2789 -0.2437 -0.2380 -0.2093 -0.1806  1.2107  1.2629  1.8587
 sugw:  kpt 131 of 1024, k=  0.00000  0.50000  0.00000
 -0.5233 -0.3116 -0.2288 -0.2066 -0.2066 -0.1726  1.2319  1.2319  1.3327
 sugw:  kpt 141 of 1024, k=  0.37500  0.62500  -0.12500
 -0.3774 -0.3009 -0.2776 -0.2148 -0.1681 -0.0405  0.5621  0.9663  1.1821
 sugw:  kpt 151 of 1024, k=  -0.25000  -0.25000  0.75000
 -0.3634 -0.3200 -0.2554 -0.2079 -0.1473  0.0397  0.4660  0.9012  1.1572
 sugw:  kpt 161 of 1024, k=  0.25000  -0.25000  0.75000
 -0.3634 -0.3200 -0.2554 -0.2079 -0.1473  0.0397  0.4660  0.9012  1.1572
 sugw:  kpt 171 of 1024, k=  0.62500  -0.12500  0.62500
 -0.3453 -0.3056 -0.2682 -0.2152 -0.1647  0.1284  0.3587  0.7828  1.2762
 sugw:  kpt 181 of 1024, k=  1.00000  0.00000  0.50000
 -0.3404 -0.2955 -0.2955 -0.1968 -0.1396  0.4285  0.4285  0.5902  0.6155
 sugw:  kpt 191 of 1024, k=  -0.62500  0.12500  0.37500
 -0.3774 -0.3009 -0.2776 -0.2148 -0.1681 -0.0405  0.5621  0.9663  1.1821
 sugw:  kpt 201 of 1024, k=  -0.25000  0.25000  0.50000
 -0.4511 -0.2863 -0.2587 -0.2366 -0.1582 -0.1470  0.6352  1.3547  1.4103
 sugw:  kpt 211 of 1024, k=  0.12500  0.37500  0.37500
 -0.4922 -0.2866 -0.2532 -0.2352 -0.1983 -0.1439  0.7755  1.1719  1.7128
 sugw:  kpt 221 of 1024, k=  0.50000  0.50000  0.25000
 -0.4005 -0.2686 -0.2449 -0.2413 -0.1662 -0.0619  0.3921  1.1777  1.5248
 sugw:  kpt 231 of 1024, k=  -0.12500  -0.37500  -0.87500
 -0.3439 -0.3246 -0.2647 -0.1916 -0.1488  0.2169  0.4414  0.6938  0.8573
 sugw:  kpt 241 of 1024, k=  -0.62500  0.62500  0.12500
 -0.3453 -0.3056 -0.2682 -0.2152 -0.1647  0.1284  0.3587  0.7828  1.2762
 sugw:  kpt 251 of 1024, k=  -0.25000  0.75000  0.00000
 -0.3750 -0.3345 -0.2451 -0.1726 -0.1619 -0.0136  0.7779  0.8329  0.8619
 sugw:  kpt 261 of 1024, k=  0.00000  1.00000  0.00000
 -0.3965 -0.3629 -0.1513 -0.1397 -0.1397  0.0771  0.4993  0.9323  0.9323
 sugw:  kpt 271 of 1024, k=  0.37500  -0.87500  -0.12500
 -0.3439 -0.3246 -0.2647 -0.1916 -0.1488  0.2169  0.4414  0.6938  0.8573
 sugw:  kpt 281 of 1024, k=  -0.12500  0.12500  0.87500
 -0.3835 -0.3516 -0.1882 -0.1669 -0.1442  0.0654  0.6193  0.7257  0.9928
 sugw:  kpt 291 of 1024, k=  0.25000  0.25000  0.75000
 -0.3634 -0.3200 -0.2554 -0.2079 -0.1473  0.0397  0.4660  0.9012  1.1572
 sugw:  kpt 301 of 1024, k=  -0.37500  -0.62500  -0.37500
 -0.3931 -0.2739 -0.2449 -0.2305 -0.1495 -0.0204  0.2987  1.2633  1.3915
 sugw:  kpt 311 of 1024, k=  0.00000  -0.50000  -0.50000
 -0.3831 -0.2952 -0.2790 -0.2183 -0.1943 -0.0519  0.6898  0.7768  1.4825
 sugw:  kpt 321 of 1024, k=  0.37500  -0.37500  -0.37500
 -0.4453 -0.2539 -0.2539 -0.2259 -0.1628 -0.1628  0.4805  1.5735  1.6080
 sugw:  kpt 331 of 1024, k=  -0.25000  0.75000  0.50000
 -0.3548 -0.3010 -0.2619 -0.2225 -0.1547  0.1269  0.2843  1.0228  1.0453
 sugw:  kpt 341 of 1024, k=  0.12500  0.87500  0.37500
 -0.3439 -0.3246 -0.2647 -0.1916 -0.1488  0.2169  0.4414  0.6938  0.8573
 sugw:  kpt 351 of 1024, k=  -0.50000  0.00000  -0.75000
 -0.3255 -0.3119 -0.2938 -0.2001 -0.1616  0.1510  0.5016  0.6457  0.9854
 sugw:  kpt 361 of 1024, k=  0.00000  0.00000  -0.75000
 -0.4052 -0.3482 -0.1616 -0.1616 -0.1576 -0.1237  0.8170  1.0137  1.0137
 sugw:  kpt 371 of 1024, k=  0.37500  0.12500  -0.87500
 -0.3439 -0.3246 -0.2647 -0.1916 -0.1488  0.2169  0.4414  0.6938  0.8573
 sugw:  kpt 381 of 1024, k=  -0.25000  -0.75000  0.00000
 -0.3750 -0.3345 -0.2451 -0.1726 -0.1619 -0.0136  0.7779  0.8329  0.8619
 sugw:  kpt 391 of 1024, k=  0.00000  -0.50000  0.00000
 -0.5233 -0.3116 -0.2288 -0.2066 -0.2066 -0.1726  1.2319  1.2319  1.3327
 sugw:  kpt 401 of 1024, k=  0.50000  -0.50000  0.00000
 -0.3831 -0.2952 -0.2790 -0.2183 -0.1943 -0.0519  0.6898  0.7768  1.4825
 sugw:  kpt 411 of 1024, k=  0.87500  -0.37500  -0.12500
 -0.3439 -0.3246 -0.2647 -0.1916 -0.1488  0.2169  0.4414  0.6938  0.8573
 sugw:  kpt 421 of 1024, k=  -0.75000  -0.25000  -0.25000
 -0.3634 -0.3200 -0.2554 -0.2079 -0.1473  0.0397  0.4660  0.9012  1.1572
 sugw:  kpt 431 of 1024, k=  -0.37500  -0.12500  -0.37500
 -0.4922 -0.2866 -0.2532 -0.2352 -0.1983 -0.1439  0.7755  1.1719  1.7128
 sugw:  kpt 441 of 1024, k=  0.12500  -0.12500  -0.37500
 -0.5816 -0.2896 -0.2407 -0.2323 -0.2030 -0.1782  1.0731  1.4084  1.5119
 sugw:  kpt 451 of 1024, k=  0.37500  0.12500  -0.37500
 -0.4922 -0.2866 -0.2532 -0.2352 -0.1983 -0.1439  0.7755  1.1719  1.7128
 sugw:  kpt 461 of 1024, k=  0.75000  0.25000  -0.50000
 -0.3548 -0.3010 -0.2619 -0.2225 -0.1547  0.1269  0.2843  1.0228  1.0453
 sugw:  kpt 471 of 1024, k=  0.12500  -0.62500  0.37500
 -0.3774 -0.3009 -0.2776 -0.2148 -0.1681 -0.0405  0.5621  0.9663  1.1821
 sugw:  kpt 481 of 1024, k=  -0.37500  0.37500  -0.62500
 -0.3931 -0.2739 -0.2449 -0.2305 -0.1495 -0.0204  0.2987  1.2633  1.3915
 sugw:  kpt 491 of 1024, k=  0.00000  0.50000  -0.75000
 -0.3255 -0.3119 -0.2938 -0.2001 -0.1616  0.1510  0.5016  0.6457  0.9854
 sugw:  kpt 501 of 1024, k=  -0.62500  -0.37500  0.12500
 -0.3774 -0.3009 -0.2776 -0.2148 -0.1681 -0.0405  0.5621  0.9663  1.1821
 sugw:  kpt 511 of 1024, k=  -0.25000  -0.25000  0.00000
 -0.6196 -0.2789 -0.2437 -0.2380 -0.2093 -0.1806  1.2107  1.2629  1.8587
 sugw:  kpt 521 of 1024, k=  0.08697  -0.08697  0.16303
 -0.6912 -0.2655 -0.2502 -0.2499 -0.1990 -0.1924  1.4145  1.6676  1.6767
 sugw:  kpt 531 of 1024, k=  0.46197  0.03803  0.03803
 -0.5474 -0.3050 -0.2317 -0.2143 -0.2092 -0.1746  1.2167  1.2771  1.4068
 sugw:  kpt 541 of 1024, k=  0.83697  0.16303  -0.08697
 -0.3827 -0.3484 -0.1995 -0.1687 -0.1461  0.0319  0.6811  0.7404  0.9726
 sugw:  kpt 551 of 1024, k=  -0.78803  0.28803  -0.21197
 -0.3573 -0.3241 -0.2540 -0.2016 -0.1485  0.0846  0.4508  0.8302  1.0958
 sugw:  kpt 561 of 1024, k=  -0.28803  0.28803  -0.21197
 -0.5536 -0.2774 -0.2457 -0.2443 -0.1908 -0.1782  0.8497  1.5038  1.6885
 sugw:  kpt 571 of 1024, k=  0.08697  0.41303  -0.33697
 -0.4942 -0.2915 -0.2544 -0.2290 -0.2025 -0.1422  0.8305  1.1092  1.6352
 sugw:  kpt 581 of 1024, k=  0.33697  0.66303  -0.33697
 -0.3816 -0.2853 -0.2494 -0.2355 -0.1489  0.0034  0.3370  1.1481  1.3151
 sugw:  kpt 591 of 1024, k=  -0.28803  -0.21197  0.53803
 -0.4303 -0.2906 -0.2624 -0.2322 -0.1600 -0.1232  0.6060  1.2562  1.3497
 sugw:  kpt 601 of 1024, k=  0.21197  -0.21197  0.53803
 -0.4453 -0.3011 -0.2594 -0.2267 -0.1572 -0.1414  0.6959  1.2834  1.3158
 sugw:  kpt 611 of 1024, k=  0.58697  -0.08697  0.41303
 -0.3802 -0.2968 -0.2808 -0.2160 -0.1801 -0.0475  0.6038  0.8932  1.2796
 sugw:  kpt 621 of 1024, k=  0.96197  0.03803  0.28803
 -0.3649 -0.3406 -0.2264 -0.1757 -0.1408  0.2191  0.5427  0.6357  0.7781
 sugw:  kpt 631 of 1024, k=  -0.66303  0.16303  0.16303
 -0.4009 -0.3265 -0.2479 -0.2013 -0.1510 -0.0948  0.7282  1.0276  1.1360
 sugw:  kpt 641 of 1024, k=  -0.28803  0.28803  0.28803
 -0.5278 -0.2728 -0.2470 -0.2470 -0.1803 -0.1803  0.7584  1.6224  1.6362
 sugw:  kpt 651 of 1024, k=  0.08697  0.41303  0.16303
 -0.5560 -0.2949 -0.2434 -0.2268 -0.2005 -0.1726  1.0368  1.2984  1.5315
 sugw:  kpt 661 of 1024, k=  0.46197  0.53803  0.03803
 -0.3825 -0.2953 -0.2798 -0.2177 -0.1914 -0.0512  0.6617  0.8108  1.4101
 sugw:  kpt 671 of 1024, k=  -0.16303  -0.33697  0.91303
 -0.3484 -0.3302 -0.2517 -0.1886 -0.1497  0.2452  0.4203  0.6268  0.9460
 sugw:  kpt 681 of 1024, k=  -0.66303  0.66303  -0.08697
 -0.3385 -0.3173 -0.2658 -0.2045 -0.1628  0.2007  0.3530  0.6831  1.2107
 sugw:  kpt 691 of 1024, k=  -0.28803  0.78803  -0.21197
 -0.3573 -0.3241 -0.2540 -0.2016 -0.1485  0.0846  0.4508  0.8302  1.0958
 sugw:  kpt 701 of 1024, k=  0.08697  0.91303  -0.33697
 -0.3535 -0.3313 -0.2489 -0.1839 -0.1449  0.2241  0.5218  0.6152  0.8169
 sugw:  kpt 711 of 1024, k=  -0.66303  0.16303  0.66303
 -0.3477 -0.3107 -0.2617 -0.2135 -0.1570  0.1742  0.2980  0.8103  1.2296
 sugw:  kpt 721 of 1024, k=  -0.16303  0.16303  0.66303
 -0.4009 -0.3265 -0.2479 -0.2013 -0.1510 -0.0948  0.7282  1.0276  1.1360
 sugw:  kpt 731 of 1024, k=  0.21197  0.28803  0.53803
 -0.4303 -0.2906 -0.2624 -0.2322 -0.1600 -0.1232  0.6060  1.2562  1.3497
 sugw:  kpt 741 of 1024, k=  0.58697  0.41303  0.41303
 -0.4028 -0.2662 -0.2521 -0.2050 -0.1499 -0.0487  0.2719  1.3786  1.4696
 sugw:  kpt 751 of 1024, k=  -0.03803  -0.46197  -0.71197
 -0.3337 -0.3064 -0.2963 -0.2014 -0.1668  0.0829  0.5380  0.7100  1.0386
 sugw:  kpt 761 of 1024, k=  -0.53803  0.53803  0.28803
 -0.3930 -0.2718 -0.2461 -0.2278 -0.1586 -0.0210  0.3061  1.1952  1.4742
 sugw:  kpt 771 of 1024, k=  -0.28803  0.78803  0.28803
 -0.3548 -0.3160 -0.2574 -0.2106 -0.1499  0.1116  0.3680  0.8568  1.1828
 sugw:  kpt 781 of 1024, k=  0.08697  0.91303  0.16303
 -0.3823 -0.3526 -0.1878 -0.1636 -0.1436  0.1059  0.5691  0.7192  0.9497
 sugw:  kpt 791 of 1024, k=  0.46197  -0.96197  0.03803
 -0.3406 -0.3060 -0.2852 -0.1958 -0.1408  0.3736  0.4494  0.6039  0.6497
 sugw:  kpt 801 of 1024, k=  -0.03803  0.03803  -0.96197
 -0.3952 -0.3618 -0.1524 -0.1450 -0.1405  0.0758  0.5115  0.9025  0.9380
 sugw:  kpt 811 of 1024, k=  0.33697  0.16303  0.91303
 -0.3484 -0.3302 -0.2517 -0.1886 -0.1497  0.2452  0.4203  0.6268  0.9460
 sugw:  kpt 821 of 1024, k=  -0.28803  -0.71197  -0.21197
 -0.3686 -0.3146 -0.2611 -0.2117 -0.1487  0.0040  0.4880  0.9785  1.1348
 sugw:  kpt 831 of 1024, k=  0.08697  -0.58697  -0.33697
 -0.4008 -0.3058 -0.2735 -0.2113 -0.1778 -0.0819  0.6807  0.9633  1.2409
 sugw:  kpt 841 of 1024, k=  0.46197  -0.46197  -0.21197
 -0.4164 -0.2645 -0.2603 -0.2427 -0.1753 -0.0962  0.4970  1.1676  1.5836
 sugw:  kpt 851 of 1024, k=  -0.16303  0.66303  0.66303
 -0.3477 -0.3107 -0.2617 -0.2135 -0.1570  0.1742  0.2980  0.8103  1.2296
 sugw:  kpt 861 of 1024, k=  -0.78803  -0.21197  -0.46197
 -0.3418 -0.3114 -0.2709 -0.2129 -0.1545  0.1688  0.3077  0.9200  0.9751
 sugw:  kpt 871 of 1024, k=  -0.41303  -0.08697  -0.58697
 -0.3802 -0.2968 -0.2808 -0.2160 -0.1801 -0.0475  0.6038  0.8932  1.2796
 sugw:  kpt 881 of 1024, k=  0.08697  -0.08697  -0.58697
 -0.4563 -0.3234 -0.2371 -0.1985 -0.1676 -0.1622  0.9811  1.1597  1.1631
 sugw:  kpt 891 of 1024, k=  0.46197  0.03803  -0.71197
 -0.3337 -0.3064 -0.2963 -0.2014 -0.1668  0.0829  0.5380  0.7100  1.0386
 sugw:  kpt 901 of 1024, k=  -0.28803  -0.71197  0.28803
 -0.3681 -0.3049 -0.2596 -0.2207 -0.1480  0.0266  0.4031  1.0054  1.2228
 sugw:  kpt 911 of 1024, k=  0.08697  -0.58697  0.16303
 -0.4448 -0.3203 -0.2488 -0.2029 -0.1673 -0.1440  0.8836  1.1288  1.1809
 sugw:  kpt 921 of 1024, k=  0.58697  -0.58697  0.16303
 -0.3595 -0.2932 -0.2638 -0.2256 -0.1659  0.0629  0.3645  0.8945  1.3474
 sugw:  kpt 931 of 1024, k=  0.96197  -0.46197  0.03803
 -0.3406 -0.3060 -0.2852 -0.1958 -0.1408  0.3736  0.4494  0.6039  0.6497
 sugw:  kpt 941 of 1024, k=  -0.66303  -0.33697  -0.08697
 -0.3717 -0.3136 -0.2764 -0.2012 -0.1670 -0.0335  0.6284  0.9132  1.0767
 sugw:  kpt 951 of 1024, k=  -0.28803  -0.21197  -0.21197
 -0.5820 -0.2779 -0.2443 -0.2434 -0.1937 -0.1841  0.9405  1.5738  1.6046
 sugw:  kpt 961 of 1024, k=  0.08697  -0.08697  -0.08697
 -0.7080 -0.2614 -0.2533 -0.2533 -0.1948 -0.1948  1.5058  1.7114  1.7810
 sugw:  kpt 971 of 1024, k=  0.46197  0.03803  -0.21197
 -0.5153 -0.3015 -0.2503 -0.2158 -0.2009 -0.1580  0.9975  1.1220  1.4709
 sugw:  kpt 981 of 1024, k=  0.83697  0.16303  -0.33697
 -0.3490 -0.3266 -0.2582 -0.1938 -0.1502  0.1546  0.4416  0.7487  0.9641
 sugw:  kpt 991 of 1024, k=  0.21197  -0.71197  0.53803
 -0.3539 -0.3013 -0.2623 -0.2219 -0.1560  0.1268  0.2918  0.9498  1.1227
 sugw:  kpt 1001 of 1024, k=  -0.28803  0.28803  -0.46197
 -0.4538 -0.2687 -0.2564 -0.2445 -0.1595 -0.1538  0.5810  1.4309  1.4866
 sugw:  kpt 1011 of 1024, k=  0.08697  0.41303  -0.58697
 -0.3802 -0.2968 -0.2808 -0.2160 -0.1801 -0.0475  0.6038  0.8932  1.2796
 sugw:  kpt 1021 of 1024, k=  -0.53803  -0.46197  0.28803
 -0.4022 -0.2677 -0.2476 -0.2300 -0.1599 -0.0627  0.3558  1.2509  1.4584
 
 sugw: minimum reduced hilbert space = 34 ... check nband_chi0 in GWinput
 sugw: writing file erange:  ebot=-0.728  etop=16.664
 Exit 0 bndfp 
 CPU time:   30.109s   Wall clock   30.145s   08:50:17 30.06.2014        on phpdl1.ph.kcl.ac.uk
rdcmd:  echo cu | lmf2gw_0
 lmf2gw: read files from ext=.cu        1024  qp made by           1  threads
 File CLASS:  nclass =   1
 nqbze nqpz ndimh(max)= 1024  512   44
ikp isp nband q=     1     1    39      0.000000   0.000000   0.000000 geig sumcheck=  0.437603D+01 0.328353D-04
ikp isp nband q=     2     1    38      0.125000   0.125000  -0.125000 geig sumcheck= -0.467707D+01 0.398579D-04
ikp isp nband q=     3     1    38      0.250000   0.250000  -0.250000 geig sumcheck= -0.643107D+01-0.626321D-04
ikp isp nband q=     4     1    38      0.375000   0.375000  -0.375000 geig sumcheck= -0.324333D+01 0.651122D-04
ikp isp nband q=     5     1    37      0.500000   0.500000  -0.500000 geig sumcheck= -0.646356D+01 0.248575D-04
ikp isp nband q=     6     1    38      0.625000   0.625000  -0.625000 geig sumcheck= -0.436388D+00-0.484888D-04
ikp isp nband q=     7     1    38      0.750000   0.750000  -0.750000 geig sumcheck= -0.795778D+01-0.244820D-04
ikp isp nband q=     8     1    38      0.875000   0.875000  -0.875000 geig sumcheck= -0.543301D+01 0.242571D-05
ikp isp nband q=     9     1    38      0.125000  -0.125000   0.125000 geig sumcheck= -0.571383D+01-0.194364D-04
ikp isp nband q=    10     1    36      0.250000   0.000000   0.000000 geig sumcheck=  0.154644D+01-0.124492D-04
ikp isp nband q=    11     1    37      0.375000   0.125000  -0.125000 geig sumcheck=  0.273998D+01 0.125028D-04
ikp isp nband q=    12     1    39      0.500000   0.250000  -0.250000 geig sumcheck=  0.404149D+01 0.121328D-04
ikp isp nband q=    13     1    37      0.625000   0.375000  -0.375000 geig sumcheck= -0.102181D+02-0.205008D-04
ikp isp nband q=    14     1    36      0.750000   0.500000  -0.500000 geig sumcheck=  0.298401D+01-0.182588D-03
ikp isp nband q=    15     1    37      0.875000   0.625000  -0.625000 geig sumcheck= -0.608078D+01-0.189794D-03
ikp isp nband q=    16     1    39      1.000000   0.750000  -0.750000 geig sumcheck= -0.416399D+01-0.730183D-05
ikp isp nband q=    17     1    38      0.250000  -0.250000   0.250000 geig sumcheck= -0.347244D+01 0.491672D-06
ikp isp nband q=    18     1    37      0.375000  -0.125000   0.125000 geig sumcheck=  0.274048D+01 0.650499D-05
ikp isp nband q=    19     1    36      0.500000   0.000000   0.000000 geig sumcheck= -0.309273D+01-0.157348D-04
ikp isp nband q=    20     1    36      0.625000   0.125000  -0.125000 geig sumcheck= -0.565621D+01-0.182922D-04
ikp isp nband q=    21     1    36      0.750000   0.250000  -0.250000 geig sumcheck= -0.747627D+01-0.352990D-04
ikp isp nband q=    22     1    34      0.875000   0.375000  -0.375000 geig sumcheck=  0.490242D+01-0.311569D-05
ikp isp nband q=    23     1    39      1.000000   0.500000  -0.500000 geig sumcheck= -0.196085D+01-0.450140D-04
ikp isp nband q=    24     1    37      1.125000   0.625000  -0.625000 geig sumcheck=  0.223884D+01 0.152938D-05
ikp isp nband q=    25     1    38      0.375000  -0.375000   0.375000 geig sumcheck= -0.361946D+01 0.159256D-04
ikp isp nband q=    26     1    39      0.500000  -0.250000   0.250000 geig sumcheck=  0.436324D+01 0.240444D-05
ikp isp nband q=    27     1    36      0.625000  -0.125000   0.125000 geig sumcheck= -0.619456D+01-0.358113D-05
ikp isp nband q=    28     1    36      0.750000   0.000000   0.000000 geig sumcheck= -0.145743D+01-0.141001D-04
ikp isp nband q=    29     1    36      0.875000   0.125000  -0.125000 geig sumcheck=  0.282745D+01 0.106520D-04
ikp isp nband q=    30     1    34      1.000000   0.250000  -0.250000 geig sumcheck=  0.892230D+00 0.304149D-04
ikp isp nband q=    31     1    34      1.125000   0.375000  -0.375000 geig sumcheck=  0.490240D+01 0.145670D-04
ikp isp nband q=    32     1    36      1.250000   0.500000  -0.500000 geig sumcheck=  0.541860D+01-0.112842D-04
ikp isp nband q=    33     1    37      0.500000  -0.500000   0.500000 geig sumcheck= -0.634193D+01-0.258643D-04
ikp isp nband q=    34     1    37      0.625000  -0.375000   0.375000 geig sumcheck= -0.102178D+02 0.178037D-04
ikp isp nband q=    35     1    36      0.750000  -0.250000   0.250000 geig sumcheck= -0.692362D+01 0.200348D-04
ikp isp nband q=    36     1    36      0.875000  -0.125000   0.125000 geig sumcheck=  0.425525D+01-0.537942D-04
ikp isp nband q=    37     1    35      1.000000   0.000000   0.000000 geig sumcheck= -0.391185D+01 0.203017D-04
ikp isp nband q=    38     1    36      1.125000   0.125000  -0.125000 geig sumcheck=  0.379702D+01-0.224497D-04
ikp isp nband q=    39     1    36      1.250000   0.250000  -0.250000 geig sumcheck= -0.747608D+01 0.191332D-04
ikp isp nband q=    40     1    37      1.375000   0.375000  -0.375000 geig sumcheck=  0.323638D+01-0.356146D-04
ikp isp nband q=    41     1    38      0.625000  -0.625000   0.625000 geig sumcheck=  0.306268D+01 0.218516D-03
ikp isp nband q=    42     1    36      0.750000  -0.500000   0.500000 geig sumcheck=  0.298358D+01 0.284798D-03
ikp isp nband q=    43     1    34      0.875000  -0.375000   0.375000 geig sumcheck=  0.490249D+01 0.246904D-04
ikp isp nband q=    44     1    34      1.000000  -0.250000   0.250000 geig sumcheck=  0.694862D+00 0.475951D-04
ikp isp nband q=    45     1    36      1.125000  -0.125000   0.125000 geig sumcheck=  0.522483D+01-0.230953D-04
ikp isp nband q=    46     1    36      1.250000   0.000000   0.000000 geig sumcheck= -0.953578D+00-0.150916D-04
ikp isp nband q=    47     1    36      1.375000   0.125000  -0.125000 geig sumcheck= -0.619466D+01-0.831379D-06
ikp isp nband q=    48     1    39      1.500000   0.250000  -0.250000 geig sumcheck= -0.162047D+01 0.303205D-03
ikp isp nband q=    49     1    38      0.750000  -0.750000   0.750000 geig sumcheck= -0.703847D+01 0.276514D-04
ikp isp nband q=    50     1    37      0.875000  -0.625000   0.625000 geig sumcheck= -0.605094D+01-0.145258D-03
ikp isp nband q=    51     1    39      1.000000  -0.500000   0.500000 geig sumcheck= -0.196131D+01 0.154629D-03
ikp isp nband q=    52     1    34      1.125000  -0.375000   0.375000 geig sumcheck=  0.490250D+01-0.132995D-04
ikp isp nband q=    53     1    36      1.250000  -0.250000   0.250000 geig sumcheck= -0.747650D+01-0.129060D-04
ikp isp nband q=    54     1    36      1.375000  -0.125000   0.125000 geig sumcheck= -0.619458D+01 0.152807D-04
ikp isp nband q=    55     1    36      1.500000   0.000000   0.000000 geig sumcheck= -0.205003D+01-0.125424D-04
ikp isp nband q=    56     1    37      1.625000   0.125000  -0.125000 geig sumcheck=  0.274007D+01-0.103888D-04
ikp isp nband q=    57     1    38      0.875000  -0.875000   0.875000 geig sumcheck= -0.434392D+01-0.129573D-03
ikp isp nband q=    58     1    39      1.000000  -0.750000   0.750000 geig sumcheck= -0.422021D+01-0.239454D-04
ikp isp nband q=    59     1    37      1.125000  -0.625000   0.625000 geig sumcheck=  0.223914D+01-0.569322D-05
ikp isp nband q=    60     1    36      1.250000  -0.500000   0.500000 geig sumcheck= -0.516849D+01-0.568845D-04
ikp isp nband q=    61     1    37      1.375000  -0.375000   0.375000 geig sumcheck=  0.323668D+01 0.104838D-04
ikp isp nband q=    62     1    39      1.500000  -0.250000   0.250000 geig sumcheck= -0.472223D+01 0.847834D-03
ikp isp nband q=    63     1    37      1.625000  -0.125000   0.125000 geig sumcheck=  0.274034D+01 0.110982D-04
ikp isp nband q=    64     1    36      1.750000   0.000000   0.000000 geig sumcheck= -0.686944D+00-0.215106D-04
ikp isp nband q=    65     1    38     -0.125000   0.125000   0.125000 geig sumcheck= -0.354128D+01 0.899444D-05
ikp isp nband q=    66     1    36      0.000000   0.250000   0.000000 geig sumcheck=  0.148581D+01-0.331180D-05
ikp isp nband q=    67     1    37      0.125000   0.375000  -0.125000 geig sumcheck= -0.300031D+01 0.182671D-04
ikp isp nband q=    68     1    39      0.250000   0.500000  -0.250000 geig sumcheck=  0.213749D+01-0.377972D-04
ikp isp nband q=    69     1    37      0.375000   0.625000  -0.375000 geig sumcheck=  0.671653D+01 0.561917D-05
ikp isp nband q=    70     1    36      0.500000   0.750000  -0.500000 geig sumcheck=  0.298329D+01-0.508455D-03
ikp isp nband q=    71     1    37      0.625000   0.875000  -0.625000 geig sumcheck= -0.608150D+01 0.204575D-03
ikp isp nband q=    72     1    39      0.750000   1.000000  -0.750000 geig sumcheck= -0.422041D+01-0.926955D-04
ikp isp nband q=    73     1    36      0.000000   0.000000   0.250000 geig sumcheck=  0.156974D+01 0.737930D-05
ikp isp nband q=    74     1    38      0.125000   0.125000   0.125000 geig sumcheck= -0.623633D+01 0.438563D-04
ikp isp nband q=    75     1    39      0.250000   0.250000   0.000000 geig sumcheck= -0.422017D+01 0.237712D-04
ikp isp nband q=    76     1    37      0.375000   0.375000  -0.125000 geig sumcheck= -0.254304D+01-0.270996D-04
ikp isp nband q=    77     1    36      0.500000   0.500000  -0.250000 geig sumcheck=  0.266447D+01 0.434290D-03
ikp isp nband q=    78     1    37      0.625000   0.625000  -0.375000 geig sumcheck= -0.552319D+01-0.227132D-04
ikp isp nband q=    79     1    39      0.750000   0.750000  -0.500000 geig sumcheck=  0.363624D+01-0.280796D-04
ikp isp nband q=    80     1    37      0.875000   0.875000  -0.625000 geig sumcheck= -0.353534D+01-0.122887D-03
ikp isp nband q=    81     1    37      0.125000  -0.125000   0.375000 geig sumcheck= -0.149777D+01 0.200275D-03
ikp isp nband q=    82     1    39      0.250000   0.000000   0.250000 geig sumcheck= -0.416393D+01 0.140429D-04
ikp isp nband q=    83     1    37      0.375000   0.125000   0.125000 geig sumcheck= -0.198934D+01 0.832063D-05
ikp isp nband q=    84     1    37      0.500000   0.250000   0.000000 geig sumcheck= -0.340663D+01-0.539821D-04
ikp isp nband q=    85     1    37      0.625000   0.375000  -0.125000 geig sumcheck= -0.287293D+01 0.122883D-04
ikp isp nband q=    86     1    35      0.750000   0.500000  -0.250000 geig sumcheck= -0.329978D+01 0.592064D-03
ikp isp nband q=    87     1    37      0.875000   0.625000  -0.375000 geig sumcheck=  0.260042D+01 0.266376D-04
ikp isp nband q=    88     1    37      1.000000   0.750000  -0.500000 geig sumcheck=  0.298617D+01-0.750315D-05
ikp isp nband q=    89     1    39      0.250000  -0.250000   0.500000 geig sumcheck=  0.245876D+01 0.180608D-04
ikp isp nband q=    90     1    37      0.375000  -0.125000   0.375000 geig sumcheck= -0.109538D+01-0.256677D-04
ikp isp nband q=    91     1    37      0.500000   0.000000   0.250000 geig sumcheck= -0.184987D+01 0.188460D-04
ikp isp nband q=    92     1    36      0.625000   0.125000   0.125000 geig sumcheck= -0.701689D+01 0.514667D-04
ikp isp nband q=    93     1    36      0.750000   0.250000   0.000000 geig sumcheck= -0.772623D+00 0.996424D-07
ikp isp nband q=    94     1    37      0.875000   0.375000  -0.125000 geig sumcheck= -0.375627D+01 0.810138D-05
ikp isp nband q=    95     1    37      1.000000   0.500000  -0.250000 geig sumcheck=  0.401687D+01 0.102341D-04
ikp isp nband q=    96     1    37      1.125000   0.625000  -0.375000 geig sumcheck=  0.194612D+01 0.211086D-04
ikp isp nband q=    97     1    37      0.375000  -0.375000   0.625000 geig sumcheck=  0.591597D+01-0.632171D-05
ikp isp nband q=    98     1    36      0.500000  -0.250000   0.500000 geig sumcheck=  0.266446D+01-0.441054D-03
ikp isp nband q=    99     1    37      0.625000  -0.125000   0.375000 geig sumcheck= -0.477874D+01-0.183018D-04
ikp isp nband q=   100     1    36      0.750000   0.000000   0.250000 geig sumcheck= -0.647946D+01-0.538326D-05
ikp isp nband q=   101     1    36      0.875000   0.125000   0.125000 geig sumcheck= -0.418986D+01-0.169988D-03
ikp isp nband q=   102     1    35      1.000000   0.250000   0.000000 geig sumcheck= -0.436773D+01-0.117242D-04
ikp isp nband q=   103     1    37      1.125000   0.375000  -0.125000 geig sumcheck= -0.420128D+01 0.144752D-04
ikp isp nband q=   104     1    35      1.250000   0.500000  -0.250000 geig sumcheck= -0.329973D+01 0.324623D-03
ikp isp nband q=   105     1    36      0.500000  -0.500000   0.750000 geig sumcheck=  0.298144D+01 0.468117D-03
ikp isp nband q=   106     1    37      0.625000  -0.375000   0.625000 geig sumcheck= -0.626955D+01-0.308249D-04
ikp isp nband q=   107     1    35      0.750000  -0.250000   0.500000 geig sumcheck= -0.329923D+01-0.737254D-03
ikp isp nband q=   108     1    37      0.875000  -0.125000   0.375000 geig sumcheck= -0.375614D+01-0.229732D-04
ikp isp nband q=   109     1    35      1.000000   0.000000   0.250000 geig sumcheck= -0.690909D+01-0.167040D-04
ikp isp nband q=   110     1    36      1.125000   0.125000   0.125000 geig sumcheck= -0.772028D+00 0.853710D-04
ikp isp nband q=   111     1    36      1.250000   0.250000   0.000000 geig sumcheck= -0.772627D+00-0.392275D-05
ikp isp nband q=   112     1    37      1.375000   0.375000  -0.125000 geig sumcheck= -0.184136D+01 0.203270D-04
ikp isp nband q=   113     1    37      0.625000  -0.625000   0.875000 geig sumcheck= -0.463367D+01-0.109316D-03
ikp isp nband q=   114     1    39      0.750000  -0.500000   0.750000 geig sumcheck= -0.145151D+01 0.345194D-04
ikp isp nband q=   115     1    37      0.875000  -0.375000   0.625000 geig sumcheck=  0.194705D+01-0.386311D-04
ikp isp nband q=   116     1    37      1.000000  -0.250000   0.500000 geig sumcheck=  0.421530D+01 0.113941D-04
ikp isp nband q=   117     1    37      1.125000  -0.125000   0.375000 geig sumcheck= -0.420161D+01-0.991773D-05
ikp isp nband q=   118     1    36      1.250000   0.000000   0.250000 geig sumcheck= -0.772607D+00 0.745148D-05
ikp isp nband q=   119     1    36      1.375000   0.125000   0.125000 geig sumcheck= -0.662593D+01-0.217798D-04
ikp isp nband q=   120     1    37      1.500000   0.250000   0.000000 geig sumcheck= -0.391386D+01 0.211923D-04
ikp isp nband q=   121     1    39      0.750000  -0.750000   1.000000 geig sumcheck= -0.422019D+01 0.349877D-05
ikp isp nband q=   122     1    37      0.875000  -0.625000   0.875000 geig sumcheck=  0.188877D+01 0.152870D-03
ikp isp nband q=   123     1    37      1.000000  -0.500000   0.750000 geig sumcheck= -0.269995D-01-0.166725D-04
ikp isp nband q=   124     1    37      1.125000  -0.375000   0.625000 geig sumcheck=  0.260128D+01-0.126064D-04
ikp isp nband q=   125     1    35      1.250000  -0.250000   0.500000 geig sumcheck= -0.329933D+01-0.290273D-03
ikp isp nband q=   126     1    37      1.375000  -0.125000   0.375000 geig sumcheck= -0.287365D+01 0.592313D-05
ikp isp nband q=   127     1    37      1.500000   0.000000   0.250000 geig sumcheck= -0.391436D+01 0.407133D-04
ikp isp nband q=   128     1    37      1.625000   0.125000   0.125000 geig sumcheck= -0.280024D+01-0.763819D-05
ikp isp nband q=   129     1    38     -0.250000   0.250000   0.250000 geig sumcheck= -0.668697D+01-0.229930D-04
ikp isp nband q=   130     1    37     -0.125000   0.375000   0.125000 geig sumcheck=  0.274424D+01-0.110440D-04
ikp isp nband q=   131     1    36      0.000000   0.500000   0.000000 geig sumcheck= -0.577335D+00-0.289341D-04
ikp isp nband q=   132     1    36      0.125000   0.625000  -0.125000 geig sumcheck= -0.565845D+01-0.452744D-05
ikp isp nband q=   133     1    36      0.250000   0.750000  -0.250000 geig sumcheck= -0.692318D+01-0.534252D-05
ikp isp nband q=   134     1    34      0.375000   0.875000  -0.375000 geig sumcheck=  0.490242D+01-0.183746D-04
ikp isp nband q=   135     1    39      0.500000   1.000000  -0.500000 geig sumcheck= -0.196110D+01-0.177816D-03
ikp isp nband q=   136     1    37      0.625000   1.125000  -0.625000 geig sumcheck=  0.223976D+01 0.347647D-05
ikp isp nband q=   137     1    37     -0.125000   0.125000   0.375000 geig sumcheck= -0.248621D+01 0.309708D-03
ikp isp nband q=   138     1    39      0.000000   0.250000   0.250000 geig sumcheck= -0.422036D+01 0.610990D-04
ikp isp nband q=   139     1    37      0.125000   0.375000   0.125000 geig sumcheck= -0.230127D+01 0.670389D-05
ikp isp nband q=   140     1    37      0.250000   0.500000   0.000000 geig sumcheck= -0.391547D+01-0.555815D-04
ikp isp nband q=   141     1    37      0.375000   0.625000  -0.125000 geig sumcheck= -0.271053D+01 0.345972D-04
ikp isp nband q=   142     1    35      0.500000   0.750000  -0.250000 geig sumcheck= -0.405795D+01 0.191038D-03
ikp isp nband q=   143     1    37      0.625000   0.875000  -0.375000 geig sumcheck=  0.294454D+01 0.475986D-05
ikp isp nband q=   144     1    37      0.750000   1.000000  -0.500000 geig sumcheck=  0.295759D+01 0.905715D-05
ikp isp nband q=   145     1    36      0.000000   0.000000   0.500000 geig sumcheck= -0.621640D+00 0.410262D-04
ikp isp nband q=   146     1    37      0.125000   0.125000   0.375000 geig sumcheck= -0.718193D+01-0.505424D-05
ikp isp nband q=   147     1    38      0.250000   0.250000   0.250000 geig sumcheck= -0.528695D+01 0.494289D-05
ikp isp nband q=   148     1    37      0.375000   0.375000   0.125000 geig sumcheck=  0.572810D+01 0.261272D-04
ikp isp nband q=   149     1    39      0.500000   0.500000   0.000000 geig sumcheck= -0.110128D+01 0.804429D-04
ikp isp nband q=   150     1    34      0.625000   0.625000  -0.125000 geig sumcheck=  0.490242D+01-0.219356D-04
ikp isp nband q=   151     1    36      0.750000   0.750000  -0.250000 geig sumcheck=  0.320514D+01 0.108394D-04
ikp isp nband q=   152     1    36      0.875000   0.875000  -0.375000 geig sumcheck= -0.315205D+01 0.863581D-05
ikp isp nband q=   153     1    36      0.125000  -0.125000   0.625000 geig sumcheck= -0.565746D+01-0.142825D-04
ikp isp nband q=   154     1    37      0.250000   0.000000   0.500000 geig sumcheck= -0.340688D+01-0.536549D-05
ikp isp nband q=   155     1    37      0.375000   0.125000   0.375000 geig sumcheck= -0.161439D+01-0.173219D-04
ikp isp nband q=   156     1    39      0.500000   0.250000   0.250000 geig sumcheck= -0.881122D+01 0.395250D-04
ikp isp nband q=   157     1    37      0.625000   0.375000   0.125000 geig sumcheck= -0.259316D+01 0.382346D-03
ikp isp nband q=   158     1    37      0.750000   0.500000   0.000000 geig sumcheck= -0.274311D+01-0.846503D-05
ikp isp nband q=   159     1    37      0.875000   0.625000  -0.125000 geig sumcheck=  0.807124D+01 0.306133D-04
ikp isp nband q=   160     1    36      1.000000   0.750000  -0.250000 geig sumcheck=  0.320918D+01-0.952523D-05
ikp isp nband q=   161     1    36      0.250000  -0.250000   0.750000 geig sumcheck= -0.692327D+01 0.115040D-04
ikp isp nband q=   162     1    37      0.375000  -0.125000   0.625000 geig sumcheck= -0.287633D+01-0.371647D-04
ikp isp nband q=   163     1    39      0.500000   0.000000   0.500000 geig sumcheck= -0.387484D+01 0.835684D-04
ikp isp nband q=   164     1    37      0.625000   0.125000   0.375000 geig sumcheck= -0.469562D+01-0.261191D-05
ikp isp nband q=   165     1    36      0.750000   0.250000   0.250000 geig sumcheck= -0.806773D+01 0.213886D-03
ikp isp nband q=   166     1    37      0.875000   0.375000   0.125000 geig sumcheck=  0.142953D+00-0.205915D-06
ikp isp nband q=   167     1    39      1.000000   0.500000   0.000000 geig sumcheck=  0.543690D+01 0.404123D-04
ikp isp nband q=   168     1    37      1.125000   0.625000  -0.125000 geig sumcheck=  0.807112D+01 0.459480D-04
ikp isp nband q=   169     1    34      0.375000  -0.375000   0.875000 geig sumcheck=  0.490248D+01 0.295284D-04
ikp isp nband q=   170     1    35      0.500000  -0.250000   0.750000 geig sumcheck= -0.389938D+01 0.221178D-04
ikp isp nband q=   171     1    34      0.625000  -0.125000   0.625000 geig sumcheck=  0.490244D+01 0.111713D-04
ikp isp nband q=   172     1    37      0.750000   0.000000   0.500000 geig sumcheck= -0.274332D+01 0.164842D-05
ikp isp nband q=   173     1    37      0.875000   0.125000   0.375000 geig sumcheck= -0.833623D+01-0.144040D-05
ikp isp nband q=   174     1    34      1.000000   0.250000   0.250000 geig sumcheck= -0.628311D+00-0.378033D-04
ikp isp nband q=   175     1    37      1.125000   0.375000   0.125000 geig sumcheck= -0.385127D+00-0.133023D-03
ikp isp nband q=   176     1    37      1.250000   0.500000   0.000000 geig sumcheck= -0.274301D+01-0.580712D-05
ikp isp nband q=   177     1    39      0.500000  -0.500000   1.000000 geig sumcheck= -0.196128D+01 0.205344D-03
ikp isp nband q=   178     1    37      0.625000  -0.375000   0.875000 geig sumcheck=  0.294483D+01-0.373800D-04
ikp isp nband q=   179     1    36      0.750000  -0.250000   0.750000 geig sumcheck=  0.530550D+01-0.855490D-05
ikp isp nband q=   180     1    37      0.875000  -0.125000   0.625000 geig sumcheck=  0.754330D+01-0.412703D-04
ikp isp nband q=   181     1    39      1.000000   0.000000   0.500000 geig sumcheck= -0.518870D+01-0.342274D-04
ikp isp nband q=   182     1    37      1.125000   0.125000   0.375000 geig sumcheck= -0.151178D-01 0.891683D-04
ikp isp nband q=   183     1    36      1.250000   0.250000   0.250000 geig sumcheck= -0.902015D+01 0.565334D-04
ikp isp nband q=   184     1    37      1.375000   0.375000   0.125000 geig sumcheck= -0.212803D+01 0.111724D-02
ikp isp nband q=   185     1    37      0.625000  -0.625000   1.125000 geig sumcheck=  0.198065D+01-0.878594D-05
ikp isp nband q=   186     1    37      0.750000  -0.500000   1.000000 geig sumcheck=  0.268146D+01-0.143913D-04
ikp isp nband q=   187     1    36      0.875000  -0.375000   0.875000 geig sumcheck= -0.228118D+01-0.108873D-05
ikp isp nband q=   188     1    36      1.000000  -0.250000   0.750000 geig sumcheck=  0.320917D+01 0.610325D-05
ikp isp nband q=   189     1    37      1.125000  -0.125000   0.625000 geig sumcheck=  0.807153D+01-0.297701D-04
ikp isp nband q=   190     1    37      1.250000   0.000000   0.500000 geig sumcheck= -0.274351D+01-0.121053D-04
ikp isp nband q=   191     1    37      1.375000   0.125000   0.375000 geig sumcheck= -0.332926D+01-0.986538D-03
ikp isp nband q=   192     1    39      1.500000   0.250000   0.250000 geig sumcheck= -0.687447D+01 0.146036D-04
ikp isp nband q=   193     1    38     -0.375000   0.375000   0.375000 geig sumcheck= -0.686068D-01 0.826263D-04
ikp isp nband q=   194     1    39     -0.250000   0.500000   0.250000 geig sumcheck=  0.436249D+01 0.256015D-04
ikp isp nband q=   195     1    36     -0.125000   0.625000   0.125000 geig sumcheck= -0.619673D+01 0.200115D-05
ikp isp nband q=   196     1    36      0.000000   0.750000   0.000000 geig sumcheck= -0.763579D+00-0.132482D-04
ikp isp nband q=   197     1    36      0.125000   0.875000  -0.125000 geig sumcheck=  0.425486D+01-0.787214D-04
ikp isp nband q=   198     1    34      0.250000   1.000000  -0.250000 geig sumcheck=  0.146652D+00-0.235805D-03
ikp isp nband q=   199     1    34      0.375000   1.125000  -0.375000 geig sumcheck=  0.475742D+01-0.108722D-04
ikp isp nband q=   200     1    36      0.500000   1.250000  -0.500000 geig sumcheck= -0.516750D+01 0.104355D-03
ikp isp nband q=   201     1    39     -0.250000   0.250000   0.500000 geig sumcheck=  0.213624D+01-0.180455D-04
ikp isp nband q=   202     1    37     -0.125000   0.375000   0.375000 geig sumcheck= -0.112541D+01-0.269662D-03
ikp isp nband q=   203     1    37      0.000000   0.500000   0.250000 geig sumcheck= -0.184858D+01 0.344399D-04
ikp isp nband q=   204     1    36      0.125000   0.625000   0.125000 geig sumcheck= -0.706815D+01-0.339301D-04
ikp isp nband q=   205     1    36      0.250000   0.750000   0.000000 geig sumcheck= -0.304870D+00-0.939062D-05
ikp isp nband q=   206     1    37      0.375000   0.875000  -0.125000 geig sumcheck= -0.375504D+01 0.298785D-04
ikp isp nband q=   207     1    37      0.500000   1.000000  -0.250000 geig sumcheck= -0.649937D+00-0.808430D-05
ikp isp nband q=   208     1    37      0.625000   1.125000  -0.375000 geig sumcheck=  0.194734D+01 0.161395D-04
ikp isp nband q=   209     1    36     -0.125000   0.125000   0.625000 geig sumcheck= -0.619538D+01-0.180700D-04
ikp isp nband q=   210     1    37      0.000000   0.250000   0.500000 geig sumcheck= -0.184910D+01 0.457612D-04
ikp isp nband q=   211     1    37      0.125000   0.375000   0.375000 geig sumcheck= -0.112554D+01-0.270717D-03
ikp isp nband q=   212     1    39      0.250000   0.500000   0.250000 geig sumcheck= -0.742415D+01 0.102184D-04
ikp isp nband q=   213     1    37      0.375000   0.625000   0.125000 geig sumcheck= -0.334053D+01 0.303213D-03
ikp isp nband q=   214     1    37      0.500000   0.750000   0.000000 geig sumcheck= -0.274356D+01 0.110117D-04
ikp isp nband q=   215     1    37      0.625000   0.875000  -0.125000 geig sumcheck=  0.807068D+01 0.168360D-04
ikp isp nband q=   216     1    36      0.750000   1.000000  -0.250000 geig sumcheck=  0.357502D+01 0.559777D-06
ikp isp nband q=   217     1    36      0.000000   0.000000   0.750000 geig sumcheck= -0.204238D+01 0.764846D-05
ikp isp nband q=   218     1    36      0.125000   0.125000   0.625000 geig sumcheck= -0.321213D+01 0.307996D-04
ikp isp nband q=   219     1    39      0.250000   0.250000   0.500000 geig sumcheck= -0.449415D+01-0.268465D-04
ikp isp nband q=   220     1    38      0.375000   0.375000   0.375000 geig sumcheck=  0.164066D+01 0.114080D-02
ikp isp nband q=   221     1    36      0.500000   0.500000   0.250000 geig sumcheck= -0.905127D+01-0.323986D-04
ikp isp nband q=   222     1    34      0.625000   0.625000   0.125000 geig sumcheck=  0.339233D+01 0.239289D-04
ikp isp nband q=   223     1    34      0.750000   0.750000   0.000000 geig sumcheck=  0.371391D+01-0.552508D-04
ikp isp nband q=   224     1    36      0.875000   0.875000  -0.125000 geig sumcheck=  0.123720D+01-0.246967D-05
ikp isp nband q=   225     1    36      0.125000  -0.125000   0.875000 geig sumcheck=  0.531383D+01 0.129516D-03
ikp isp nband q=   226     1    36      0.250000   0.000000   0.750000 geig sumcheck= -0.772411D+00 0.749817D-05
ikp isp nband q=   227     1    37      0.375000   0.125000   0.625000 geig sumcheck= -0.334073D+01 0.381357D-04
ikp isp nband q=   228     1    36      0.500000   0.250000   0.500000 geig sumcheck=  0.298234D+01 0.179055D-03
ikp isp nband q=   229     1    37      0.625000   0.375000   0.375000 geig sumcheck= -0.766939D+01-0.202355D-04
ikp isp nband q=   230     1    35      0.750000   0.500000   0.250000 geig sumcheck= -0.383021D+01-0.278457D-04
ikp isp nband q=   231     1    37      0.875000   0.625000   0.125000 geig sumcheck=  0.468691D+01 0.336302D-04
ikp isp nband q=   232     1    35      1.000000   0.750000   0.000000 geig sumcheck=  0.768763D+01-0.129407D-04
ikp isp nband q=   233     1    34      0.250000  -0.250000   1.000000 geig sumcheck=  0.315604D+01 0.110291D-03
ikp isp nband q=   234     1    37      0.375000  -0.125000   0.875000 geig sumcheck= -0.375509D+01-0.330988D-04
ikp isp nband q=   235     1    37      0.500000   0.000000   0.750000 geig sumcheck= -0.274417D+01-0.753592D-05
ikp isp nband q=   236     1    34      0.625000   0.125000   0.625000 geig sumcheck=  0.332134D+01-0.104576D-03
ikp isp nband q=   237     1    35      0.750000   0.250000   0.500000 geig sumcheck= -0.390067D+01 0.548737D-04
ikp isp nband q=   238     1    34      0.875000   0.375000   0.375000 geig sumcheck=  0.475745D+01 0.421570D-04
ikp isp nband q=   239     1    37      1.000000   0.500000   0.250000 geig sumcheck=  0.163121D+01 0.963475D-04
ikp isp nband q=   240     1    37      1.125000   0.625000   0.125000 geig sumcheck=  0.414013D+01-0.401287D-04
ikp isp nband q=   241     1    34      0.375000  -0.375000   1.125000 geig sumcheck=  0.475737D+01-0.799235D-05
ikp isp nband q=   242     1    37      0.500000  -0.250000   1.000000 geig sumcheck=  0.460225D+01-0.151229D-04
ikp isp nband q=   243     1    37      0.625000  -0.125000   0.875000 geig sumcheck=  0.754237D+01-0.382139D-04
ikp isp nband q=   244     1    34      0.750000   0.000000   0.750000 geig sumcheck=  0.298331D+01 0.624840D-05
ikp isp nband q=   245     1    37      0.875000   0.125000   0.625000 geig sumcheck=  0.468667D+01 0.315780D-04
ikp isp nband q=   246     1    37      1.000000   0.250000   0.500000 geig sumcheck=  0.163133D+01 0.107459D-03
ikp isp nband q=   247     1    34      1.125000   0.375000   0.375000 geig sumcheck=  0.271487D+01-0.237150D-04
ikp isp nband q=   248     1    35      1.250000   0.500000   0.250000 geig sumcheck= -0.383027D+01 0.588319D-04
ikp isp nband q=   249     1    36      0.500000  -0.500000   1.250000 geig sumcheck= -0.991536D+01-0.496349D-04
ikp isp nband q=   250     1    37      0.625000  -0.375000   1.125000 geig sumcheck=  0.194770D+01-0.164172D-04
ikp isp nband q=   251     1    36      0.750000  -0.250000   1.000000 geig sumcheck=  0.357506D+01 0.450329D-05
ikp isp nband q=   252     1    36      0.875000  -0.125000   0.875000 geig sumcheck=  0.113949D+01-0.134719D-04
ikp isp nband q=   253     1    35      1.000000   0.000000   0.750000 geig sumcheck=  0.486289D+01 0.414197D-05
ikp isp nband q=   254     1    37      1.125000   0.125000   0.625000 geig sumcheck=  0.437002D+01-0.271195D-04
ikp isp nband q=   255     1    35      1.250000   0.250000   0.500000 geig sumcheck= -0.352847D+01-0.233091D-04
ikp isp nband q=   256     1    37      1.375000   0.375000   0.375000 geig sumcheck= -0.510814D+00 0.599737D-07
ikp isp nband q=   257     1    37     -0.500000   0.500000   0.500000 geig sumcheck= -0.276958D+01 0.755280D-05
ikp isp nband q=   258     1    37     -0.375000   0.625000   0.375000 geig sumcheck= -0.111485D+02 0.381597D-06
ikp isp nband q=   259     1    36     -0.250000   0.750000   0.250000 geig sumcheck= -0.692373D+01 0.184546D-05
ikp isp nband q=   260     1    36     -0.125000   0.875000   0.125000 geig sumcheck=  0.425467D+01-0.101265D-03
ikp isp nband q=   261     1    35      0.000000   1.000000   0.000000 geig sumcheck= -0.375699D+01 0.441607D-04
ikp isp nband q=   262     1    36      0.125000   1.125000  -0.125000 geig sumcheck=  0.379666D+01 0.125913D-03
ikp isp nband q=   263     1    36      0.250000   1.250000  -0.250000 geig sumcheck= -0.692316D+01 0.659955D-05
ikp isp nband q=   264     1    37      0.375000   1.375000  -0.375000 geig sumcheck=  0.243199D+01 0.800717D-04
ikp isp nband q=   265     1    37     -0.375000   0.375000   0.625000 geig sumcheck=  0.671465D+01 0.180652D-04
ikp isp nband q=   266     1    36     -0.250000   0.500000   0.500000 geig sumcheck=  0.266308D+01-0.477796D-04
ikp isp nband q=   267     1    37     -0.125000   0.625000   0.375000 geig sumcheck= -0.477841D+01 0.244917D-04
ikp isp nband q=   268     1    36      0.000000   0.750000   0.250000 geig sumcheck= -0.647966D+01 0.248401D-05
ikp isp nband q=   269     1    36      0.125000   0.875000   0.125000 geig sumcheck= -0.328336D+01-0.128358D-03
ikp isp nband q=   270     1    35      0.250000   1.000000   0.000000 geig sumcheck= -0.896287D+00-0.829763D-06
ikp isp nband q=   271     1    37      0.375000   1.125000  -0.125000 geig sumcheck= -0.420035D+01-0.785136D-05
ikp isp nband q=   272     1    35      0.500000   1.250000  -0.250000 geig sumcheck= -0.345894D+01 0.127532D-03
ikp isp nband q=   273     1    36     -0.250000   0.250000   0.750000 geig sumcheck= -0.747657D+01-0.291963D-04
ikp isp nband q=   274     1    37     -0.125000   0.375000   0.625000 geig sumcheck= -0.494312D+01-0.884113D-06
ikp isp nband q=   275     1    39      0.000000   0.500000   0.500000 geig sumcheck= -0.473469D+01 0.615927D-04
ikp isp nband q=   276     1    37      0.125000   0.625000   0.375000 geig sumcheck= -0.467849D+01 0.116863D-03
ikp isp nband q=   277     1    36      0.250000   0.750000   0.250000 geig sumcheck= -0.162075D+01-0.270341D-03
ikp isp nband q=   278     1    37      0.375000   0.875000   0.125000 geig sumcheck= -0.543069D+00-0.173520D-04
ikp isp nband q=   279     1    39      0.500000   1.000000   0.000000 geig sumcheck= -0.308329D+01-0.148038D-04
ikp isp nband q=   280     1    37      0.625000   1.125000  -0.125000 geig sumcheck=  0.754257D+01 0.253352D-06
ikp isp nband q=   281     1    36     -0.125000   0.125000   0.875000 geig sumcheck=  0.425521D+01 0.879272D-05
ikp isp nband q=   282     1    36      0.000000   0.250000   0.750000 geig sumcheck= -0.647945D+01 0.610196D-06
ikp isp nband q=   283     1    37      0.125000   0.375000   0.625000 geig sumcheck= -0.543124D+01 0.142196D-03
ikp isp nband q=   284     1    36      0.250000   0.500000   0.500000 geig sumcheck=  0.266394D+01-0.676357D-04
ikp isp nband q=   285     1    37      0.375000   0.625000   0.375000 geig sumcheck= -0.822368D+01-0.100099D-03
ikp isp nband q=   286     1    35      0.500000   0.750000   0.250000 geig sumcheck= -0.352739D+01-0.103752D-03
ikp isp nband q=   287     1    37      0.625000   0.875000   0.125000 geig sumcheck=  0.513233D+01-0.225928D-04
ikp isp nband q=   288     1    35      0.750000   1.000000   0.000000 geig sumcheck=  0.458358D+00 0.118815D-04
ikp isp nband q=   289     1    35      0.000000   0.000000   1.000000 geig sumcheck= -0.260371D+01-0.933517D-04
ikp isp nband q=   290     1    36      0.125000   0.125000   0.875000 geig sumcheck= -0.349165D+01 0.193068D-03
ikp isp nband q=   291     1    36      0.250000   0.250000   0.750000 geig sumcheck= -0.246172D+01 0.328778D-03
ikp isp nband q=   292     1    37      0.375000   0.375000   0.625000 geig sumcheck= -0.262397D+01 0.132587D-04
ikp isp nband q=   293     1    37      0.500000   0.500000   0.500000 geig sumcheck= -0.233010D+01 0.112586D-04
ikp isp nband q=   294     1    37      0.625000   0.625000   0.375000 geig sumcheck=  0.144080D+01-0.461315D-04
ikp isp nband q=   295     1    36      0.750000   0.750000   0.250000 geig sumcheck= -0.173267D+01 0.976104D-05
ikp isp nband q=   296     1    36      0.875000   0.875000   0.125000 geig sumcheck=  0.252636D+01-0.222498D-05
ikp isp nband q=   297     1    36      0.125000  -0.125000   1.125000 geig sumcheck=  0.434424D+01 0.247412D-05
ikp isp nband q=   298     1    35      0.250000   0.000000   1.000000 geig sumcheck= -0.165234D+01-0.345716D-05
ikp isp nband q=   299     1    37      0.375000   0.125000   0.875000 geig sumcheck=  0.143653D+00 0.985576D-05
ikp isp nband q=   300     1    35      0.500000   0.250000   0.750000 geig sumcheck= -0.370225D+01 0.103148D-03
ikp isp nband q=   301     1    37      0.625000   0.375000   0.625000 geig sumcheck= -0.111495D+02 0.199581D-04
ikp isp nband q=   302     1    36      0.750000   0.500000   0.500000 geig sumcheck= -0.991496D+01-0.496911D-04
ikp isp nband q=   303     1    37      0.875000   0.625000   0.375000 geig sumcheck=  0.169223D+01-0.805615D-05
ikp isp nband q=   304     1    36      1.000000   0.750000   0.250000 geig sumcheck=  0.404276D+01-0.317065D-05
ikp isp nband q=   305     1    36      0.250000  -0.250000   1.250000 geig sumcheck= -0.692338D+01-0.126935D-04
ikp isp nband q=   306     1    37      0.375000  -0.125000   1.125000 geig sumcheck= -0.420042D+01 0.313118D-05
ikp isp nband q=   307     1    39      0.500000   0.000000   1.000000 geig sumcheck=  0.330359D+01 0.191936D-04
ikp isp nband q=   308     1    37      0.625000   0.125000   0.875000 geig sumcheck=  0.468720D+01 0.636772D-04
ikp isp nband q=   309     1    36      0.750000   0.250000   0.750000 geig sumcheck=  0.381613D+01-0.151390D-04
ikp isp nband q=   310     1    37      0.875000   0.375000   0.625000 geig sumcheck=  0.492990D+01-0.145944D-04
ikp isp nband q=   311     1    39      1.000000   0.500000   0.500000 geig sumcheck=  0.206879D+01 0.224370D-04
ikp isp nband q=   312     1    37      1.125000   0.625000   0.375000 geig sumcheck=  0.297180D+01-0.979143D-06
ikp isp nband q=   313     1    37      0.375000  -0.375000   1.375000 geig sumcheck= -0.742400D+01 0.126045D-03
ikp isp nband q=   314     1    35      0.500000  -0.250000   1.250000 geig sumcheck= -0.330075D+01-0.168242D-03
ikp isp nband q=   315     1    37      0.625000  -0.125000   1.125000 geig sumcheck=  0.754236D+01-0.103306D-04
ikp isp nband q=   316     1    35      0.750000   0.000000   1.000000 geig sumcheck=  0.248352D+01-0.138944D-04
ikp isp nband q=   317     1    36      0.875000   0.125000   0.875000 geig sumcheck=  0.217655D+01 0.116876D-04
ikp isp nband q=   318     1    36      1.000000   0.250000   0.750000 geig sumcheck=  0.404280D+01 0.107446D-05
ikp isp nband q=   319     1    37      1.125000   0.375000   0.625000 geig sumcheck=  0.331636D+01 0.119245D-04
ikp isp nband q=   320     1    36      1.250000   0.500000   0.500000 geig sumcheck= -0.173602D+01-0.253272D-04
ikp isp nband q=   321     1    38     -0.625000   0.625000   0.625000 geig sumcheck=  0.148935D+01-0.516712D-04
ikp isp nband q=   322     1    36     -0.500000   0.750000   0.500000 geig sumcheck=  0.298187D+01-0.360632D-03
ikp isp nband q=   323     1    34     -0.375000   0.875000   0.375000 geig sumcheck=  0.490266D+01-0.145873D-04
ikp isp nband q=   324     1    34     -0.250000   1.000000   0.250000 geig sumcheck=  0.122993D+01-0.195056D-03
ikp isp nband q=   325     1    36     -0.125000   1.125000   0.125000 geig sumcheck=  0.522428D+01 0.144448D-03
ikp isp nband q=   326     1    36      0.000000   1.250000   0.000000 geig sumcheck= -0.763614D+00-0.814568D-06
ikp isp nband q=   327     1    36      0.125000   1.375000  -0.125000 geig sumcheck= -0.565852D+01 0.853750D-05
ikp isp nband q=   328     1    39      0.250000   1.500000  -0.250000 geig sumcheck= -0.159167D+01 0.607599D-03
ikp isp nband q=   329     1    36     -0.500000   0.500000   0.750000 geig sumcheck=  0.298315D+01 0.143433D-03
ikp isp nband q=   330     1    37     -0.375000   0.625000   0.625000 geig sumcheck= -0.329846D+01 0.483212D-04
ikp isp nband q=   331     1    35     -0.250000   0.750000   0.500000 geig sumcheck= -0.345801D+01 0.425862D-04
ikp isp nband q=   332     1    37     -0.125000   0.875000   0.375000 geig sumcheck= -0.375668D+01 0.200743D-04
ikp isp nband q=   333     1    35      0.000000   1.000000   0.250000 geig sumcheck= -0.559050D+01-0.132206D-04
ikp isp nband q=   334     1    36      0.125000   1.125000   0.125000 geig sumcheck= -0.406762D+01 0.496871D-03
ikp isp nband q=   335     1    36      0.250000   1.250000   0.000000 geig sumcheck= -0.772545D+00-0.326408D-05
ikp isp nband q=   336     1    37      0.375000   1.375000  -0.125000 geig sumcheck= -0.184015D+01-0.144709D-05
ikp isp nband q=   337     1    34     -0.375000   0.375000   0.875000 geig sumcheck=  0.490263D+01 0.358926D-05
ikp isp nband q=   338     1    35     -0.250000   0.500000   0.750000 geig sumcheck= -0.389999D+01 0.514099D-04
ikp isp nband q=   339     1    34     -0.125000   0.625000   0.625000 geig sumcheck=  0.475737D+01-0.710335D-06
ikp isp nband q=   340     1    37      0.000000   0.750000   0.500000 geig sumcheck= -0.274481D+01 0.972121D-05
ikp isp nband q=   341     1    37      0.125000   0.875000   0.375000 geig sumcheck= -0.796628D+01-0.475527D-05
ikp isp nband q=   342     1    34      0.250000   1.000000   0.250000 geig sumcheck= -0.548437D+01 0.233984D-04
ikp isp nband q=   343     1    37      0.375000   1.125000   0.125000 geig sumcheck= -0.384431D+00 0.418287D-04
ikp isp nband q=   344     1    37      0.500000   1.250000   0.000000 geig sumcheck= -0.274360D+01 0.278839D-04
ikp isp nband q=   345     1    34     -0.250000   0.250000   1.000000 geig sumcheck=  0.225573D+01 0.107921D-03
ikp isp nband q=   346     1    37     -0.125000   0.375000   0.875000 geig sumcheck= -0.375632D+01-0.166937D-04
ikp isp nband q=   347     1    37      0.000000   0.500000   0.750000 geig sumcheck= -0.274312D+01-0.861537D-05
ikp isp nband q=   348     1    34      0.125000   0.625000   0.625000 geig sumcheck=  0.499152D+01 0.192417D-05
ikp isp nband q=   349     1    35      0.250000   0.750000   0.500000 geig sumcheck= -0.405904D+01-0.808079D-04
ikp isp nband q=   350     1    34      0.375000   0.875000   0.375000 geig sumcheck=  0.490268D+01-0.140737D-04
ikp isp nband q=   351     1    37      0.500000   1.000000   0.250000 geig sumcheck=  0.163117D+01-0.680211D-04
ikp isp nband q=   352     1    37      0.625000   1.125000   0.125000 geig sumcheck=  0.398023D+01 0.138280D-04
ikp isp nband q=   353     1    36     -0.125000   0.125000   1.125000 geig sumcheck=  0.531414D+01-0.120292D-03
ikp isp nband q=   354     1    35      0.000000   0.250000   1.000000 geig sumcheck= -0.502027D+01 0.119268D-04
ikp isp nband q=   355     1    37      0.125000   0.375000   0.875000 geig sumcheck= -0.849432D+01-0.264109D-05
ikp isp nband q=   356     1    35      0.250000   0.500000   0.750000 geig sumcheck= -0.389992D+01-0.989786D-04
ikp isp nband q=   357     1    37      0.375000   0.625000   0.625000 geig sumcheck=  0.658750D+01 0.156465D-04
ikp isp nband q=   358     1    36      0.500000   0.750000   0.500000 geig sumcheck=  0.425038D+01 0.475367D-04
ikp isp nband q=   359     1    37      0.625000   0.875000   0.375000 geig sumcheck=  0.168945D+01-0.285264D-04
ikp isp nband q=   360     1    36      0.750000   1.000000   0.250000 geig sumcheck=  0.404263D+01-0.917132D-05
ikp isp nband q=   361     1    36      0.000000   0.000000   1.250000 geig sumcheck= -0.145303D+01 0.247327D-05
ikp isp nband q=   362     1    36      0.125000   0.125000   1.125000 geig sumcheck= -0.452479D+01-0.320352D-03
ikp isp nband q=   363     1    34      0.250000   0.250000   1.000000 geig sumcheck= -0.616206D+01 0.198160D-04
ikp isp nband q=   364     1    34      0.375000   0.375000   0.875000 geig sumcheck=  0.490257D+01-0.217958D-04
ikp isp nband q=   365     1    36      0.500000   0.500000   0.750000 geig sumcheck=  0.541810D+01-0.120310D-03
ikp isp nband q=   366     1    38      0.625000   0.625000   0.625000 geig sumcheck= -0.941655D-01-0.338682D-04
ikp isp nband q=   367     1    39      0.750000   0.750000   0.500000 geig sumcheck=  0.214865D+01-0.361268D-04
ikp isp nband q=   368     1    36      0.875000   0.875000   0.375000 geig sumcheck= -0.234898D+01-0.267052D-04
ikp isp nband q=   369     1    36      0.125000  -0.125000   1.375000 geig sumcheck= -0.619628D+01 0.949202D-05
ikp isp nband q=   370     1    36      0.250000   0.000000   1.250000 geig sumcheck= -0.772348D+00 0.388713D-05
ikp isp nband q=   371     1    37      0.375000   0.125000   1.125000 geig sumcheck= -0.543099D+00 0.111221D-04
ikp isp nband q=   372     1    37      0.500000   0.250000   1.000000 geig sumcheck=  0.163166D+01 0.823349D-04
ikp isp nband q=   373     1    37      0.625000   0.375000   0.875000 geig sumcheck=  0.168982D+01 0.203336D-04
ikp isp nband q=   374     1    39      0.750000   0.500000   0.750000 geig sumcheck=  0.134610D+01-0.216072D-05
ikp isp nband q=   375     1    37      0.875000   0.625000   0.625000 geig sumcheck=  0.132389D+01-0.659352D-03
ikp isp nband q=   376     1    37      1.000000   0.750000   0.500000 geig sumcheck= -0.206376D+00-0.688889D-04
ikp isp nband q=   377     1    39      0.250000  -0.250000   1.500000 geig sumcheck= -0.162141D+01-0.853467D-03
ikp isp nband q=   378     1    37      0.375000  -0.125000   1.375000 geig sumcheck= -0.167433D+01-0.165553D-05
ikp isp nband q=   379     1    37      0.500000   0.000000   1.250000 geig sumcheck= -0.274413D+01-0.311509D-04
ikp isp nband q=   380     1    37      0.625000   0.125000   1.125000 geig sumcheck=  0.361072D+01 0.145936D-04
ikp isp nband q=   381     1    36      0.750000   0.250000   1.000000 geig sumcheck=  0.404271D+01 0.272011D-05
ikp isp nband q=   382     1    36      0.875000   0.375000   0.875000 geig sumcheck= -0.244730D+01-0.212548D-04
ikp isp nband q=   383     1    37      1.000000   0.500000   0.750000 geig sumcheck= -0.713450D+00-0.120062D-04
ikp isp nband q=   384     1    37      1.125000   0.625000   0.625000 geig sumcheck=  0.251449D+01-0.130828D-04
ikp isp nband q=   385     1    38     -0.750000   0.750000   0.750000 geig sumcheck= -0.500076D+01-0.131185D-04
ikp isp nband q=   386     1    37     -0.625000   0.875000   0.625000 geig sumcheck= -0.608080D+01 0.217862D-03
ikp isp nband q=   387     1    39     -0.500000   1.000000   0.500000 geig sumcheck= -0.196113D+01-0.160910D-03
ikp isp nband q=   388     1    34     -0.375000   1.125000   0.375000 geig sumcheck=  0.475765D+01-0.156653D-04
ikp isp nband q=   389     1    36     -0.250000   1.250000   0.250000 geig sumcheck= -0.747655D+01 0.126663D-05
ikp isp nband q=   390     1    36     -0.125000   1.375000   0.125000 geig sumcheck= -0.619680D+01 0.132899D-04
ikp isp nband q=   391     1    36      0.000000   1.500000   0.000000 geig sumcheck= -0.577468D+00-0.418136D-04
ikp isp nband q=   392     1    37      0.125000   1.625000  -0.125000 geig sumcheck=  0.274316D+01-0.168446D-04
ikp isp nband q=   393     1    37     -0.625000   0.625000   0.875000 geig sumcheck= -0.741193D+01-0.973184D-03
ikp isp nband q=   394     1    39     -0.500000   0.750000   0.750000 geig sumcheck= -0.396059D+01 0.336324D-05
ikp isp nband q=   395     1    37     -0.375000   0.875000   0.625000 geig sumcheck=  0.228875D+01-0.134208D-04
ikp isp nband q=   396     1    37     -0.250000   1.000000   0.500000 geig sumcheck=  0.843740D+01 0.695993D-05
ikp isp nband q=   397     1    37     -0.125000   1.125000   0.375000 geig sumcheck= -0.375679D+01-0.147936D-04
ikp isp nband q=   398     1    36      0.000000   1.250000   0.250000 geig sumcheck= -0.772340D+00-0.524723D-07
ikp isp nband q=   399     1    36      0.125000   1.375000   0.125000 geig sumcheck= -0.706820D+01-0.745075D-04
ikp isp nband q=   400     1    37      0.250000   1.500000   0.000000 geig sumcheck= -0.391568D+01 0.497849D-04
ikp isp nband q=   401     1    39     -0.500000   0.500000   1.000000 geig sumcheck= -0.196154D+01 0.617062D-04
ikp isp nband q=   402     1    37     -0.375000   0.625000   0.875000 geig sumcheck=  0.228859D+01 0.222803D-04
ikp isp nband q=   403     1    36     -0.250000   0.750000   0.750000 geig sumcheck=  0.308384D+00 0.133247D-04
ikp isp nband q=   404     1    37     -0.125000   0.875000   0.625000 geig sumcheck=  0.807067D+01-0.813107D-05
ikp isp nband q=   405     1    39      0.000000   1.000000   0.500000 geig sumcheck=  0.556816D+01-0.136159D-04
ikp isp nband q=   406     1    37      0.125000   1.125000   0.375000 geig sumcheck= -0.154395D-01 0.954271D-04
ikp isp nband q=   407     1    36      0.250000   1.250000   0.250000 geig sumcheck= -0.558150D+01-0.496977D-03
ikp isp nband q=   408     1    37      0.375000   1.375000   0.125000 geig sumcheck= -0.223872D+01-0.237534D-03
ikp isp nband q=   409     1    34     -0.375000   0.375000   1.125000 geig sumcheck=  0.475751D+01 0.693657D-05
ikp isp nband q=   410     1    37     -0.250000   0.500000   1.000000 geig sumcheck=  0.501657D+01-0.563904D-05
ikp isp nband q=   411     1    37     -0.125000   0.625000   0.875000 geig sumcheck=  0.807127D+01 0.999990D-05
ikp isp nband q=   412     1    34      0.000000   0.750000   0.750000 geig sumcheck=  0.138210D+01-0.999354D-04
ikp isp nband q=   413     1    37      0.125000   0.875000   0.625000 geig sumcheck=  0.513096D+01-0.165812D-04
ikp isp nband q=   414     1    37      0.250000   1.000000   0.500000 geig sumcheck=  0.163116D+01-0.496710D-06
ikp isp nband q=   415     1    34      0.375000   1.125000   0.375000 geig sumcheck=  0.548385D+01-0.100529D-03
ikp isp nband q=   416     1    35      0.500000   1.250000   0.250000 geig sumcheck= -0.398941D+01-0.115593D-03
ikp isp nband q=   417     1    36     -0.250000   0.250000   1.250000 geig sumcheck= -0.692380D+01 0.139529D-04
ikp isp nband q=   418     1    37     -0.125000   0.375000   1.125000 geig sumcheck= -0.420135D+01 0.150323D-04
ikp isp nband q=   419     1    39      0.000000   0.500000   1.000000 geig sumcheck= -0.173879D+00 0.327747D-04
ikp isp nband q=   420     1    37      0.125000   0.625000   0.875000 geig sumcheck=  0.513263D+01 0.551029D-04
ikp isp nband q=   421     1    36      0.250000   0.750000   0.750000 geig sumcheck= -0.234059D+00-0.816004D-05
ikp isp nband q=   422     1    37      0.375000   0.875000   0.625000 geig sumcheck=  0.443633D+01 0.241765D-04
ikp isp nband q=   423     1    39      0.500000   1.000000   0.500000 geig sumcheck=  0.565853D+01 0.192466D-04
ikp isp nband q=   424     1    37      0.625000   1.125000   0.375000 geig sumcheck=  0.331599D+01 0.833755D-05
ikp isp nband q=   425     1    36     -0.125000   0.125000   1.375000 geig sumcheck= -0.619530D+01 0.434236D-05
ikp isp nband q=   426     1    36      0.000000   0.250000   1.250000 geig sumcheck= -0.304796D+00-0.774215D-05
ikp isp nband q=   427     1    37      0.125000   0.375000   1.125000 geig sumcheck= -0.384901D+00-0.151138D-03
ikp isp nband q=   428     1    37      0.250000   0.500000   1.000000 geig sumcheck=  0.163089D+01 0.962933D-04
ikp isp nband q=   429     1    37      0.375000   0.625000   0.875000 geig sumcheck=  0.443575D+01-0.204965D-05
ikp isp nband q=   430     1    39      0.500000   0.750000   0.750000 geig sumcheck=  0.444469D+01-0.276731D-04
ikp isp nband q=   431     1    37      0.625000   0.875000   0.625000 geig sumcheck=  0.576714D+01 0.287794D-03
ikp isp nband q=   432     1    37      0.750000   1.000000   0.500000 geig sumcheck= -0.205135D+00-0.325104D-04
ikp isp nband q=   433     1    36      0.000000   0.000000   1.500000 geig sumcheck= -0.207248D+01 0.500107D-04
ikp isp nband q=   434     1    36      0.125000   0.125000   1.375000 geig sumcheck= -0.258077D+01-0.240027D-04
ikp isp nband q=   435     1    36      0.250000   0.250000   1.250000 geig sumcheck= -0.564136D+01 0.245643D-03
ikp isp nband q=   436     1    34      0.375000   0.375000   1.125000 geig sumcheck=  0.395759D+01-0.433665D-04
ikp isp nband q=   437     1    39      0.500000   0.500000   1.000000 geig sumcheck=  0.203611D+01-0.190001D-04
ikp isp nband q=   438     1    37      0.625000   0.625000   0.875000 geig sumcheck=  0.150602D+01-0.724708D-03
ikp isp nband q=   439     1    38      0.750000   0.750000   0.750000 geig sumcheck= -0.422253D+01 0.461739D-05
ikp isp nband q=   440     1    37      0.875000   0.875000   0.625000 geig sumcheck=  0.905319D+00-0.138234D-03
ikp isp nband q=   441     1    37      0.125000  -0.125000   1.625000 geig sumcheck=  0.273689D+01 0.129239D-04
ikp isp nband q=   442     1    37      0.250000   0.000000   1.500000 geig sumcheck= -0.391508D+01 0.358748D-04
ikp isp nband q=   443     1    37      0.375000   0.125000   1.375000 geig sumcheck= -0.207525D+01-0.303628D-03
ikp isp nband q=   444     1    35      0.500000   0.250000   1.250000 geig sumcheck= -0.398938D+01 0.732400D-04
ikp isp nband q=   445     1    37      0.625000   0.375000   1.125000 geig sumcheck=  0.348148D+01 0.712428D-05
ikp isp nband q=   446     1    37      0.750000   0.500000   1.000000 geig sumcheck= -0.204591D+00-0.251532D-04
ikp isp nband q=   447     1    37      0.875000   0.625000   0.875000 geig sumcheck= -0.160162D+01 0.122140D-03
ikp isp nband q=   448     1    39      1.000000   0.750000   0.750000 geig sumcheck= -0.163781D+01 0.172879D-04
ikp isp nband q=   449     1    38     -0.875000   0.875000   0.875000 geig sumcheck= -0.461458D+01 0.272807D-04
ikp isp nband q=   450     1    39     -0.750000   1.000000   0.750000 geig sumcheck= -0.416372D+01-0.897607D-04
ikp isp nband q=   451     1    37     -0.625000   1.125000   0.625000 geig sumcheck=  0.198239D+01-0.142115D-04
ikp isp nband q=   452     1    36     -0.500000   1.250000   0.500000 geig sumcheck=  0.671494D+00-0.814069D-04
ikp isp nband q=   453     1    37     -0.375000   1.375000   0.375000 geig sumcheck= -0.742187D+01-0.512489D-04
ikp isp nband q=   454     1    39     -0.250000   1.500000   0.250000 geig sumcheck= -0.492530D+01-0.510964D-03
ikp isp nband q=   455     1    37     -0.125000   1.625000   0.125000 geig sumcheck=  0.274408D+01-0.127061D-04
ikp isp nband q=   456     1    36      0.000000   1.750000   0.000000 geig sumcheck=  0.148574D+01-0.188679D-04
ikp isp nband q=   457     1    39     -0.750000   0.750000   1.000000 geig sumcheck= -0.416351D+01 0.135443D-04
ikp isp nband q=   458     1    37     -0.625000   0.875000   0.875000 geig sumcheck= -0.266158D+00 0.470551D-04
ikp isp nband q=   459     1    37     -0.500000   1.000000   0.750000 geig sumcheck=  0.260146D+01-0.112213D-04
ikp isp nband q=   460     1    37     -0.375000   1.125000   0.625000 geig sumcheck=  0.260013D+01 0.833262D-06
ikp isp nband q=   461     1    35     -0.250000   1.250000   0.500000 geig sumcheck= -0.345797D+01-0.160718D-03
ikp isp nband q=   462     1    37     -0.125000   1.375000   0.375000 geig sumcheck= -0.287482D+01 0.291377D-05
ikp isp nband q=   463     1    37      0.000000   1.500000   0.250000 geig sumcheck= -0.340695D+01 0.676348D-04
ikp isp nband q=   464     1    37      0.125000   1.625000   0.125000 geig sumcheck= -0.245270D+01-0.158913D-04
ikp isp nband q=   465     1    37     -0.625000   0.625000   1.125000 geig sumcheck=  0.172229D+01 0.114157D-04
ikp isp nband q=   466     1    37     -0.500000   0.750000   1.000000 geig sumcheck= -0.137047D+01 0.239577D-04
ikp isp nband q=   467     1    36     -0.375000   0.875000   0.875000 geig sumcheck= -0.131880D+01 0.224488D-04
ikp isp nband q=   468     1    36     -0.250000   1.000000   0.750000 geig sumcheck=  0.320923D+01 0.375412D-05
ikp isp nband q=   469     1    37     -0.125000   1.125000   0.625000 geig sumcheck=  0.754256D+01-0.252439D-04
ikp isp nband q=   470     1    37      0.000000   1.250000   0.500000 geig sumcheck= -0.274488D+01 0.263201D-04
ikp isp nband q=   471     1    37      0.125000   1.375000   0.375000 geig sumcheck= -0.321187D+01-0.671826D-03
ikp isp nband q=   472     1    39      0.250000   1.500000   0.250000 geig sumcheck= -0.704268D+01-0.387064D-04
ikp isp nband q=   473     1    36     -0.500000   0.500000   1.250000 geig sumcheck=  0.670790D+00 0.377110D-04
ikp isp nband q=   474     1    37     -0.375000   0.625000   1.125000 geig sumcheck=  0.228860D+01 0.198411D-04
ikp isp nband q=   475     1    36     -0.250000   0.750000   1.000000 geig sumcheck=  0.357503D+01-0.321203D-05
ikp isp nband q=   476     1    36     -0.125000   0.875000   0.875000 geig sumcheck=  0.283242D+01-0.185531D-05
ikp isp nband q=   477     1    35      0.000000   1.000000   0.750000 geig sumcheck=  0.610087D+01-0.255854D-05
ikp isp nband q=   478     1    37      0.125000   1.125000   0.625000 geig sumcheck=  0.436995D+01-0.326229D-04
ikp isp nband q=   479     1    35      0.250000   1.250000   0.500000 geig sumcheck= -0.370346D+01-0.119796D-03
ikp isp nband q=   480     1    37      0.375000   1.375000   0.375000 geig sumcheck=  0.861881D+00-0.869755D-04
ikp isp nband q=   481     1    37     -0.375000   0.375000   1.375000 geig sumcheck= -0.822409D+01-0.223154D-04
ikp isp nband q=   482     1    35     -0.250000   0.500000   1.250000 geig sumcheck= -0.329991D+01 0.156217D-03
ikp isp nband q=   483     1    37     -0.125000   0.625000   1.125000 geig sumcheck=  0.754298D+01 0.410154D-04
ikp isp nband q=   484     1    35      0.000000   0.750000   1.000000 geig sumcheck=  0.310613D+01-0.197382D-04
ikp isp nband q=   485     1    36      0.125000   0.875000   0.875000 geig sumcheck=  0.304986D+01-0.187970D-05
ikp isp nband q=   486     1    36      0.250000   1.000000   0.750000 geig sumcheck=  0.404274D+01-0.779047D-05
ikp isp nband q=   487     1    37      0.375000   1.125000   0.625000 geig sumcheck=  0.379165D+01-0.225456D-04
ikp isp nband q=   488     1    36      0.500000   1.250000   0.500000 geig sumcheck= -0.675163D+01-0.199808D-04
ikp isp nband q=   489     1    39     -0.250000   0.250000   1.500000 geig sumcheck= -0.489663D+01 0.249534D-03
ikp isp nband q=   490     1    37     -0.125000   0.375000   1.375000 geig sumcheck= -0.270711D+01-0.109401D-04
ikp isp nband q=   491     1    37      0.000000   0.500000   1.250000 geig sumcheck= -0.274297D+01-0.223591D-04
ikp isp nband q=   492     1    37      0.125000   0.625000   1.125000 geig sumcheck=  0.473964D+01 0.441432D-04
ikp isp nband q=   493     1    36      0.250000   0.750000   1.000000 geig sumcheck=  0.357505D+01 0.549503D-05
ikp isp nband q=   494     1    36      0.375000   0.875000   0.875000 geig sumcheck= -0.330721D+01 0.245437D-05
ikp isp nband q=   495     1    37      0.500000   1.000000   0.750000 geig sumcheck= -0.715077D+00-0.193086D-04
ikp isp nband q=   496     1    37      0.625000   1.125000   0.625000 geig sumcheck=  0.617799D+01 0.418650D-04
ikp isp nband q=   497     1    37     -0.125000   0.125000   1.625000 geig sumcheck=  0.274695D+01-0.475914D-05
ikp isp nband q=   498     1    37      0.000000   0.250000   1.500000 geig sumcheck= -0.340553D+01-0.833298D-05
ikp isp nband q=   499     1    37      0.125000   0.375000   1.375000 geig sumcheck= -0.336080D+01 0.916460D-03
ikp isp nband q=   500     1    35      0.250000   0.500000   1.250000 geig sumcheck= -0.336973D+01 0.155628D-03
ikp isp nband q=   501     1    37      0.375000   0.625000   1.125000 geig sumcheck=  0.313722D+01-0.311414D-05
ikp isp nband q=   502     1    37      0.500000   0.750000   1.000000 geig sumcheck= -0.206184D+00-0.640409D-04
ikp isp nband q=   503     1    37      0.625000   0.875000   0.875000 geig sumcheck= -0.878919D+00-0.138840D-03
ikp isp nband q=   504     1    39      0.750000   1.000000   0.750000 geig sumcheck=  0.743368D+00-0.774026D-04
ikp isp nband q=   505     1    36      0.000000   0.000000   1.750000 geig sumcheck=  0.195362D+01 0.772500D-05
ikp isp nband q=   506     1    37      0.125000   0.125000   1.625000 geig sumcheck= -0.264847D+00 0.857814D-05
ikp isp nband q=   507     1    39      0.250000   0.250000   1.500000 geig sumcheck= -0.453218D+01-0.423958D-05
ikp isp nband q=   508     1    37      0.375000   0.375000   1.375000 geig sumcheck=  0.598805D+01 0.174103D-05
ikp isp nband q=   509     1    36      0.500000   0.500000   1.250000 geig sumcheck=  0.318967D+01-0.415087D-04
ikp isp nband q=   510     1    37      0.625000   0.625000   1.125000 geig sumcheck=  0.583921D+01-0.337970D-04
ikp isp nband q=   511     1    39      0.750000   0.750000   1.000000 geig sumcheck= -0.135551D+01 0.283559D-04
ikp isp nband q=   512     1    38      0.875000   0.875000   0.875000 geig sumcheck=  0.233205D+01-0.467045D-05
ikp isp nband q=   513     1    39     -0.038032   0.038032   0.038032 geig sumcheck=  0.519394D+01 0.220963D-04
ikp isp nband q=   514     1    39      0.086968   0.163032  -0.086968 geig sumcheck=  0.135853D+01-0.180832D-03
ikp isp nband q=   515     1    39      0.211968   0.288032  -0.211968 geig sumcheck=  0.687600D+01 0.178546D-04
ikp isp nband q=   516     1    38      0.336968   0.413032  -0.336968 geig sumcheck=  0.460959D+01-0.757339D-03
ikp isp nband q=   517     1    37      0.461968   0.538032  -0.461968 geig sumcheck= -0.532987D+01 0.144671D-03
ikp isp nband q=   518     1    37      0.586968   0.663032  -0.586968 geig sumcheck= -0.946786D+01-0.144498D-03
ikp isp nband q=   519     1    38      0.711968   0.788032  -0.711968 geig sumcheck= -0.291422D+01-0.687178D-05
ikp isp nband q=   520     1    39      0.836968   0.913032  -0.836968 geig sumcheck= -0.394952D+01 0.432256D-04
ikp isp nband q=   521     1    39      0.086968  -0.086968   0.163032 geig sumcheck=  0.101967D+01 0.471928D-03
ikp isp nband q=   522     1    39      0.211968   0.038032   0.038032 geig sumcheck= -0.588800D+01 0.455364D-04
ikp isp nband q=   523     1    38      0.336968   0.163032  -0.086968 geig sumcheck= -0.605351D+01 0.587518D-04
ikp isp nband q=   524     1    39      0.461968   0.288032  -0.211968 geig sumcheck=  0.293102D+01-0.454016D-04
ikp isp nband q=   525     1    37      0.586968   0.413032  -0.336968 geig sumcheck= -0.482313D+01-0.666826D-04
ikp isp nband q=   526     1    37      0.711968   0.538032  -0.461968 geig sumcheck=  0.669177D+01-0.430291D-04
ikp isp nband q=   527     1    38      0.836968   0.663032  -0.586968 geig sumcheck=  0.459549D+01 0.420170D-05
ikp isp nband q=   528     1    39      0.961968   0.788032  -0.711968 geig sumcheck= -0.565006D+01 0.370825D-04
ikp isp nband q=   529     1    39      0.211968  -0.211968   0.288032 geig sumcheck=  0.667707D+01-0.182820D-04
ikp isp nband q=   530     1    38      0.336968  -0.086968   0.163032 geig sumcheck=  0.402065D+00 0.187053D-04
ikp isp nband q=   531     1    36      0.461968   0.038032   0.038032 geig sumcheck= -0.250218D+01 0.564013D-05
ikp isp nband q=   532     1    36      0.586968   0.163032  -0.086968 geig sumcheck= -0.422433D+01 0.345573D-04
ikp isp nband q=   533     1    37      0.711968   0.288032  -0.211968 geig sumcheck= -0.253553D+01-0.229430D-05
ikp isp nband q=   534     1    35      0.836968   0.413032  -0.336968 geig sumcheck=  0.583068D+01 0.147675D-04
ikp isp nband q=   535     1    37      0.961968   0.538032  -0.461968 geig sumcheck= -0.221989D+00 0.853118D-04
ikp isp nband q=   536     1    38      1.086968   0.663032  -0.586968 geig sumcheck=  0.139308D+00 0.431734D-04
ikp isp nband q=   537     1    38      0.336968  -0.336968   0.413032 geig sumcheck=  0.548068D+01 0.357961D-03
ikp isp nband q=   538     1    39      0.461968  -0.211968   0.288032 geig sumcheck=  0.306366D+01 0.529374D-04
ikp isp nband q=   539     1    36      0.586968  -0.086968   0.163032 geig sumcheck= -0.342402D+01 0.338406D-06
ikp isp nband q=   540     1    36      0.711968   0.038032   0.038032 geig sumcheck= -0.680849D+01-0.145841D-04
ikp isp nband q=   541     1    36      0.836968   0.163032  -0.086968 geig sumcheck= -0.402176D+01-0.400326D-04
ikp isp nband q=   542     1    35      0.961968   0.288032  -0.211968 geig sumcheck= -0.674900D+01-0.208104D-04
ikp isp nband q=   543     1    35      1.086968   0.413032  -0.336968 geig sumcheck=  0.201228D+01-0.624111D-04
ikp isp nband q=   544     1    36      1.211968   0.538032  -0.461968 geig sumcheck=  0.129141D+01 0.223008D-04
ikp isp nband q=   545     1    37      0.461968  -0.461968   0.538032 geig sumcheck= -0.129285D+01-0.126108D-03
ikp isp nband q=   546     1    37      0.586968  -0.336968   0.413032 geig sumcheck= -0.473058D+01 0.158032D-04
ikp isp nband q=   547     1    37      0.711968  -0.211968   0.288032 geig sumcheck= -0.274533D+01 0.221785D-04
ikp isp nband q=   548     1    36      0.836968  -0.086968   0.163032 geig sumcheck= -0.223321D+01-0.104778D-04
ikp isp nband q=   549     1    35      0.961968   0.038032   0.038032 geig sumcheck= -0.554765D+01-0.380508D-04
ikp isp nband q=   550     1    36      1.086968   0.163032  -0.086968 geig sumcheck= -0.172836D+01 0.821001D-05
ikp isp nband q=   551     1    36      1.211968   0.288032  -0.211968 geig sumcheck= -0.181793D+01 0.942346D-05
ikp isp nband q=   552     1    36      1.336968   0.413032  -0.336968 geig sumcheck= -0.446985D+01-0.129631D-03
ikp isp nband q=   553     1    37      0.586968  -0.586968   0.663032 geig sumcheck= -0.938846D+01-0.667641D-04
ikp isp nband q=   554     1    37      0.711968  -0.461968   0.538032 geig sumcheck=  0.620864D+01 0.651216D-04
ikp isp nband q=   555     1    35      0.836968  -0.336968   0.413032 geig sumcheck=  0.584074D+01 0.221878D-04
ikp isp nband q=   556     1    35      0.961968  -0.211968   0.288032 geig sumcheck= -0.659967D+01 0.212854D-04
ikp isp nband q=   557     1    36      1.086968  -0.086968   0.163032 geig sumcheck= -0.215875D+01 0.621191D-05
ikp isp nband q=   558     1    36      1.211968   0.038032   0.038032 geig sumcheck= -0.330198D+01 0.866939D-06
ikp isp nband q=   559     1    36      1.336968   0.163032  -0.086968 geig sumcheck= -0.519472D+01-0.585900D-04
ikp isp nband q=   560     1    38      1.461968   0.288032  -0.211968 geig sumcheck=  0.336255D+01 0.288533D-04
ikp isp nband q=   561     1    38      0.711968  -0.711968   0.788032 geig sumcheck= -0.292621D+01 0.331103D-04
ikp isp nband q=   562     1    38      0.836968  -0.586968   0.663032 geig sumcheck=  0.582370D+01 0.272753D-04
ikp isp nband q=   563     1    37      0.961968  -0.461968   0.538032 geig sumcheck=  0.935456D+00-0.175582D-03
ikp isp nband q=   564     1    35      1.086968  -0.336968   0.413032 geig sumcheck=  0.201207D+01-0.593788D-04
ikp isp nband q=   565     1    36      1.211968  -0.211968   0.288032 geig sumcheck= -0.181813D+01 0.138099D-04
ikp isp nband q=   566     1    36      1.336968  -0.086968   0.163032 geig sumcheck= -0.479195D+01 0.157277D-04
ikp isp nband q=   567     1    36      1.461968   0.038032   0.038032 geig sumcheck= -0.154025D+01-0.156821D-04
ikp isp nband q=   568     1    36      1.586968   0.163032  -0.086968 geig sumcheck= -0.196859D+01 0.303969D-04
ikp isp nband q=   569     1    39      0.836968  -0.836968   0.913032 geig sumcheck= -0.315275D+01 0.399442D-03
ikp isp nband q=   570     1    39      0.961968  -0.711968   0.788032 geig sumcheck= -0.605417D+01-0.167798D-03
ikp isp nband q=   571     1    38      1.086968  -0.586968   0.663032 geig sumcheck=  0.144680D+00-0.302147D-04
ikp isp nband q=   572     1    36      1.211968  -0.461968   0.538032 geig sumcheck=  0.170046D+01-0.433449D-04
ikp isp nband q=   573     1    36      1.336968  -0.336968   0.413032 geig sumcheck= -0.446797D+01-0.402157D-04
ikp isp nband q=   574     1    38      1.461968  -0.211968   0.288032 geig sumcheck= -0.788266D+01 0.174200D-04
ikp isp nband q=   575     1    36      1.586968  -0.086968   0.163032 geig sumcheck= -0.173759D+01 0.607552D-05
ikp isp nband q=   576     1    37      1.711968   0.038032   0.038032 geig sumcheck= -0.474795D+01-0.321989D-04
ikp isp nband q=   577     1    41     -0.163032   0.163032   0.163032 geig sumcheck= -0.760921D+00 0.254718D-04
ikp isp nband q=   578     1    37     -0.038032   0.288032   0.038032 geig sumcheck= -0.168176D+01 0.250146D-04
ikp isp nband q=   579     1    36      0.086968   0.413032  -0.086968 geig sumcheck= -0.116598D+01 0.419409D-05
ikp isp nband q=   580     1    36      0.211968   0.538032  -0.211968 geig sumcheck= -0.481294D+01-0.680665D-05
ikp isp nband q=   581     1    36      0.336968   0.663032  -0.336968 geig sumcheck= -0.118063D+01-0.234576D-03
ikp isp nband q=   582     1    36      0.461968   0.788032  -0.461968 geig sumcheck=  0.707724D+01-0.189640D-04
ikp isp nband q=   583     1    37      0.586968   0.913032  -0.586968 geig sumcheck= -0.268883D+01-0.391359D-04
ikp isp nband q=   584     1    39      0.711968   1.038032  -0.711968 geig sumcheck= -0.476682D+01 0.467802D-04
ikp isp nband q=   585     1    37     -0.038032   0.038032   0.288032 geig sumcheck=  0.649917D+01 0.243181D-04
ikp isp nband q=   586     1    39      0.086968   0.163032   0.163032 geig sumcheck=  0.115860D+00-0.150032D-03
ikp isp nband q=   587     1    39      0.211968   0.288032   0.038032 geig sumcheck= -0.986438D+00 0.880561D-04
ikp isp nband q=   588     1    38      0.336968   0.413032  -0.086968 geig sumcheck=  0.556555D+01 0.366726D-05
ikp isp nband q=   589     1    36      0.461968   0.538032  -0.211968 geig sumcheck= -0.108611D+02 0.699346D-04
ikp isp nband q=   590     1    36      0.586968   0.663032  -0.336968 geig sumcheck= -0.457770D+01 0.342724D-04
ikp isp nband q=   591     1    38      0.711968   0.788032  -0.461968 geig sumcheck=  0.167153D+01-0.268963D-03
ikp isp nband q=   592     1    36      0.836968   0.913032  -0.586968 geig sumcheck= -0.729089D+01 0.154841D-04
ikp isp nband q=   593     1    36      0.086968  -0.086968   0.413032 geig sumcheck= -0.115889D+01 0.235555D-05
ikp isp nband q=   594     1    39      0.211968   0.038032   0.288032 geig sumcheck= -0.111400D+01-0.574897D-04
ikp isp nband q=   595     1    37      0.336968   0.163032   0.163032 geig sumcheck= -0.539729D+01-0.101612D-04
ikp isp nband q=   596     1    37      0.461968   0.288032   0.038032 geig sumcheck= -0.296994D+01 0.893308D-04
ikp isp nband q=   597     1    37      0.586968   0.413032  -0.086968 geig sumcheck= -0.475684D+01 0.108111D-04
ikp isp nband q=   598     1    35      0.711968   0.538032  -0.211968 geig sumcheck=  0.465102D+01 0.382414D-03
ikp isp nband q=   599     1    37      0.836968   0.663032  -0.336968 geig sumcheck=  0.113588D+01-0.705403D-04
ikp isp nband q=   600     1    37      0.961968   0.788032  -0.461968 geig sumcheck=  0.497833D+01 0.370668D-04
ikp isp nband q=   601     1    36      0.211968  -0.211968   0.538032 geig sumcheck= -0.435785D+01-0.774102D-05
ikp isp nband q=   602     1    38      0.336968  -0.086968   0.413032 geig sumcheck=  0.559284D+01-0.620477D-05
ikp isp nband q=   603     1    37      0.461968   0.038032   0.288032 geig sumcheck= -0.674396D+01 0.239165D-04
ikp isp nband q=   604     1    36      0.586968   0.163032   0.163032 geig sumcheck= -0.707989D+01 0.573408D-04
ikp isp nband q=   605     1    37      0.711968   0.288032   0.038032 geig sumcheck= -0.543967D+01-0.108907D-04
ikp isp nband q=   606     1    37      0.836968   0.413032  -0.086968 geig sumcheck= -0.163609D+01 0.132627D-04
ikp isp nband q=   607     1    37      0.961968   0.538032  -0.211968 geig sumcheck=  0.159026D+01-0.245909D-04
ikp isp nband q=   608     1    37      1.086968   0.663032  -0.336968 geig sumcheck=  0.265673D+01-0.421896D-04
ikp isp nband q=   609     1    36      0.336968  -0.336968   0.663032 geig sumcheck= -0.118296D+01 0.167990D-03
ikp isp nband q=   610     1    36      0.461968  -0.211968   0.538032 geig sumcheck= -0.112571D+02-0.611670D-04
ikp isp nband q=   611     1    37      0.586968  -0.086968   0.413032 geig sumcheck= -0.723724D+01-0.110720D-04
ikp isp nband q=   612     1    37      0.711968   0.038032   0.288032 geig sumcheck= -0.222244D+01 0.153051D-04
ikp isp nband q=   613     1    36      0.836968   0.163032   0.163032 geig sumcheck= -0.890163D+00 0.820914D-05
ikp isp nband q=   614     1    35      0.961968   0.288032   0.038032 geig sumcheck= -0.292228D+01 0.409281D-04
ikp isp nband q=   615     1    37      1.086968   0.413032  -0.086968 geig sumcheck= -0.373274D+01 0.155911D-03
ikp isp nband q=   616     1    35      1.211968   0.538032  -0.211968 geig sumcheck= -0.576213D+01 0.508234D-02
ikp isp nband q=   617     1    36      0.461968  -0.461968   0.788032 geig sumcheck=  0.708006D+01 0.468197D-04
ikp isp nband q=   618     1    36      0.586968  -0.336968   0.663032 geig sumcheck= -0.558550D+01 0.920013D-06
ikp isp nband q=   619     1    35      0.711968  -0.211968   0.538032 geig sumcheck=  0.445418D+01-0.372485D-03
ikp isp nband q=   620     1    37      0.836968  -0.086968   0.413032 geig sumcheck= -0.110106D+01-0.173131D-04
ikp isp nband q=   621     1    35      0.961968   0.038032   0.288032 geig sumcheck= -0.418514D+01 0.169459D-04
ikp isp nband q=   622     1    36      1.086968   0.163032   0.163032 geig sumcheck=  0.632313D+00-0.905726D-06
ikp isp nband q=   623     1    36      1.211968   0.288032   0.038032 geig sumcheck=  0.240769D+01-0.235547D-05
ikp isp nband q=   624     1    37      1.336968   0.413032  -0.086968 geig sumcheck= -0.445376D+00-0.461569D-04
ikp isp nband q=   625     1    37      0.586968  -0.586968   0.913032 geig sumcheck= -0.338741D+01 0.565040D-04
ikp isp nband q=   626     1    38      0.711968  -0.461968   0.788032 geig sumcheck=  0.149338D+01 0.208227D-03
ikp isp nband q=   627     1    37      0.836968  -0.336968   0.663032 geig sumcheck=  0.104215D+01 0.527730D-04
ikp isp nband q=   628     1    37      0.961968  -0.211968   0.538032 geig sumcheck=  0.159015D+01 0.199513D-04
ikp isp nband q=   629     1    37      1.086968  -0.086968   0.413032 geig sumcheck= -0.352131D+01-0.158836D-03
ikp isp nband q=   630     1    36      1.211968   0.038032   0.288032 geig sumcheck= -0.743202D+01-0.848408D-05
ikp isp nband q=   631     1    36      1.336968   0.163032   0.163032 geig sumcheck= -0.956156D+01-0.739776D-05
ikp isp nband q=   632     1    37      1.461968   0.288032   0.038032 geig sumcheck= -0.168182D+01-0.198388D-05
ikp isp nband q=   633     1    39      0.711968  -0.711968   1.038032 geig sumcheck= -0.450563D+01-0.455952D-04
ikp isp nband q=   634     1    36      0.836968  -0.586968   0.913032 geig sumcheck= -0.728421D+01-0.110066D-04
ikp isp nband q=   635     1    37      0.961968  -0.461968   0.788032 geig sumcheck=  0.497837D+01 0.173372D-04
ikp isp nband q=   636     1    37      1.086968  -0.336968   0.663032 geig sumcheck=  0.257740D+01 0.280620D-04
ikp isp nband q=   637     1    35      1.211968  -0.211968   0.538032 geig sumcheck= -0.636743D+01-0.797390D-02
ikp isp nband q=   638     1    37      1.336968  -0.086968   0.413032 geig sumcheck= -0.104568D+00-0.168856D-03
ikp isp nband q=   639     1    37      1.461968   0.038032   0.288032 geig sumcheck= -0.166178D+01 0.124929D-04
ikp isp nband q=   640     1    37      1.586968   0.163032   0.163032 geig sumcheck= -0.166740D+01 0.484846D-04
ikp isp nband q=   641     1    38     -0.288032   0.288032   0.288032 geig sumcheck=  0.187482D+01-0.386317D-04
ikp isp nband q=   642     1    37     -0.163032   0.413032   0.163032 geig sumcheck= -0.174862D+01 0.340659D-04
ikp isp nband q=   643     1    36     -0.038032   0.538032   0.038032 geig sumcheck= -0.448409D+01-0.124644D-03
ikp isp nband q=   644     1    36      0.086968   0.663032  -0.086968 geig sumcheck= -0.419814D+01 0.263766D-04
ikp isp nband q=   645     1    36      0.211968   0.788032  -0.211968 geig sumcheck= -0.448865D+01-0.645748D-04
ikp isp nband q=   646     1    34      0.336968   0.913032  -0.336968 geig sumcheck=  0.430135D+01-0.106493D-04
ikp isp nband q=   647     1    37      0.461968   1.038032  -0.461968 geig sumcheck=  0.433188D+01-0.948414D-04
ikp isp nband q=   648     1    37      0.586968   1.163032  -0.586968 geig sumcheck=  0.792748D+01-0.161788D-04
ikp isp nband q=   649     1    37     -0.163032   0.163032   0.413032 geig sumcheck= -0.485269D+01 0.616308D-04
ikp isp nband q=   650     1    39     -0.038032   0.288032   0.288032 geig sumcheck= -0.824861D+00-0.127249D-04
ikp isp nband q=   651     1    36      0.086968   0.413032   0.163032 geig sumcheck= -0.504948D+01-0.270830D-04
ikp isp nband q=   652     1    37      0.211968   0.538032   0.038032 geig sumcheck=  0.107575D+01-0.210253D-03
ikp isp nband q=   653     1    37      0.336968   0.663032  -0.086968 geig sumcheck= -0.255734D+01 0.378390D-05
ikp isp nband q=   654     1    35      0.461968   0.788032  -0.211968 geig sumcheck=  0.717460D+01 0.185695D-04
ikp isp nband q=   655     1    37      0.586968   0.913032  -0.336968 geig sumcheck=  0.323066D+01 0.918499D-07
ikp isp nband q=   656     1    37      0.711968   1.038032  -0.461968 geig sumcheck=  0.503693D+01-0.435824D-04
ikp isp nband q=   657     1    36     -0.038032   0.038032   0.538032 geig sumcheck= -0.467471D+01-0.316186D-04
ikp isp nband q=   658     1    36      0.086968   0.163032   0.413032 geig sumcheck= -0.505167D+01 0.227992D-04
ikp isp nband q=   659     1    38      0.211968   0.288032   0.288032 geig sumcheck= -0.513980D+01-0.202817D-04
ikp isp nband q=   660     1    38      0.336968   0.413032   0.163032 geig sumcheck= -0.295120D+01 0.746315D-05
ikp isp nband q=   661     1    37      0.461968   0.538032   0.038032 geig sumcheck=  0.232174D+00 0.130508D-04
ikp isp nband q=   662     1    35      0.586968   0.663032  -0.086968 geig sumcheck= -0.521919D+01 0.420690D-04
ikp isp nband q=   663     1    36      0.711968   0.788032  -0.211968 geig sumcheck= -0.232723D+01-0.368212D-06
ikp isp nband q=   664     1    36      0.836968   0.913032  -0.336968 geig sumcheck= -0.519157D+01-0.595878D-04
ikp isp nband q=   665     1    36      0.086968  -0.086968   0.663032 geig sumcheck= -0.397543D+01-0.323582D-03
ikp isp nband q=   666     1    37      0.211968   0.038032   0.538032 geig sumcheck=  0.123217D+01 0.264888D-03
ikp isp nband q=   667     1    38      0.336968   0.163032   0.413032 geig sumcheck= -0.505345D+01-0.150730D-05
ikp isp nband q=   668     1    38      0.461968   0.288032   0.288032 geig sumcheck=  0.493065D+01 0.822566D-04
ikp isp nband q=   669     1    36      0.586968   0.413032   0.163032 geig sumcheck=  0.500889D+01-0.248563D-05
ikp isp nband q=   670     1    36      0.711968   0.538032   0.038032 geig sumcheck=  0.470772D+01-0.393143D-04
ikp isp nband q=   671     1    35      0.836968   0.663032  -0.086968 geig sumcheck=  0.151217D+01-0.383626D-05
ikp isp nband q=   672     1    36      0.961968   0.788032  -0.211968 geig sumcheck=  0.196209D+01-0.238928D-03
ikp isp nband q=   673     1    36      0.211968  -0.211968   0.788032 geig sumcheck= -0.448953D+01-0.215282D-04
ikp isp nband q=   674     1    37      0.336968  -0.086968   0.663032 geig sumcheck= -0.256149D+01-0.112181D-04
ikp isp nband q=   675     1    37      0.461968   0.038032   0.538032 geig sumcheck=  0.230562D+00 0.148305D-05
ikp isp nband q=   676     1    36      0.586968   0.163032   0.413032 geig sumcheck=  0.107214D+01 0.125014D-04
ikp isp nband q=   677     1    35      0.711968   0.288032   0.288032 geig sumcheck= -0.529581D+01 0.172467D-04
ikp isp nband q=   678     1    36      0.836968   0.413032   0.163032 geig sumcheck= -0.350198D+01-0.721915D-04
ikp isp nband q=   679     1    38      0.961968   0.538032   0.038032 geig sumcheck=  0.380281D+01 0.379906D-04
ikp isp nband q=   680     1    36      1.086968   0.663032  -0.086968 geig sumcheck=  0.395475D+01-0.213172D-03
ikp isp nband q=   681     1    34      0.336968  -0.336968   0.913032 geig sumcheck=  0.397116D+01 0.166098D-04
ikp isp nband q=   682     1    35      0.461968  -0.211968   0.788032 geig sumcheck=  0.717476D+01-0.855811D-05
ikp isp nband q=   683     1    35      0.586968  -0.086968   0.663032 geig sumcheck= -0.554154D+01-0.613698D-04
ikp isp nband q=   684     1    36      0.711968   0.038032   0.538032 geig sumcheck= -0.449072D+01-0.721661D-04
ikp isp nband q=   685     1    36      0.836968   0.163032   0.413032 geig sumcheck= -0.332773D+01-0.221463D-04
ikp isp nband q=   686     1    35      0.961968   0.288032   0.288032 geig sumcheck= -0.340170D+01-0.308158D-04
ikp isp nband q=   687     1    36      1.086968   0.413032   0.163032 geig sumcheck=  0.372049D+00 0.487586D-04
ikp isp nband q=   688     1    37      1.211968   0.538032   0.038032 geig sumcheck= -0.298387D+01-0.321962D-03
ikp isp nband q=   689     1    37      0.461968  -0.461968   1.038032 geig sumcheck=  0.445643D+01-0.262442D-04
ikp isp nband q=   690     1    37      0.586968  -0.336968   0.913032 geig sumcheck=  0.323181D+01 0.688469D-06
ikp isp nband q=   691     1    36      0.711968  -0.211968   0.788032 geig sumcheck= -0.262188D+01-0.706562D-04
ikp isp nband q=   692     1    35      0.836968  -0.086968   0.663032 geig sumcheck=  0.159008D+01-0.455263D-04
ikp isp nband q=   693     1    38      0.961968   0.038032   0.538032 geig sumcheck=  0.210779D+01 0.272218D-04
ikp isp nband q=   694     1    36      1.086968   0.163032   0.413032 geig sumcheck= -0.189152D+00 0.269849D-05
ikp isp nband q=   695     1    35      1.211968   0.288032   0.288032 geig sumcheck= -0.323014D+01-0.392453D-04
ikp isp nband q=   696     1    36      1.336968   0.413032   0.163032 geig sumcheck= -0.322932D+01-0.488172D-05
ikp isp nband q=   697     1    37      0.586968  -0.586968   1.163032 geig sumcheck= -0.200317D+00-0.109916D-04
ikp isp nband q=   698     1    37      0.711968  -0.461968   1.038032 geig sumcheck=  0.503582D+01-0.228423D-04
ikp isp nband q=   699     1    36      0.836968  -0.336968   0.913032 geig sumcheck= -0.519247D+01 0.523130D-04
ikp isp nband q=   700     1    36      0.961968  -0.211968   0.788032 geig sumcheck=  0.152886D+01 0.406196D-04
ikp isp nband q=   701     1    36      1.086968  -0.086968   0.663032 geig sumcheck=  0.364152D+01-0.124376D-03
ikp isp nband q=   702     1    37      1.211968   0.038032   0.538032 geig sumcheck= -0.308024D+01 0.137968D-03
ikp isp nband q=   703     1    36      1.336968   0.163032   0.413032 geig sumcheck= -0.459159D+01 0.282306D-04
ikp isp nband q=   704     1    38      1.461968   0.288032   0.288032 geig sumcheck=  0.395674D+01 0.165501D-04
ikp isp nband q=   705     1    37     -0.413032   0.413032   0.413032 geig sumcheck= -0.431374D+01 0.731317D-04
ikp isp nband q=   706     1    38     -0.288032   0.538032   0.288032 geig sumcheck= -0.435471D+01 0.899309D-04
ikp isp nband q=   707     1    36     -0.163032   0.663032   0.163032 geig sumcheck= -0.586683D+01-0.203958D-04
ikp isp nband q=   708     1    36     -0.038032   0.788032   0.038032 geig sumcheck= -0.644181D+01-0.533572D-04
ikp isp nband q=   709     1    36      0.086968   0.913032  -0.086968 geig sumcheck=  0.329438D+01-0.177536D-04
ikp isp nband q=   710     1    35      0.211968   1.038032  -0.211968 geig sumcheck= -0.132070D+01-0.438625D-03
ikp isp nband q=   711     1    34      0.336968   1.163032  -0.336968 geig sumcheck=  0.447146D+01 0.986206D-05
ikp isp nband q=   712     1    37      0.461968   1.288032  -0.461968 geig sumcheck=  0.399713D+01-0.352563D-04
ikp isp nband q=   713     1    38     -0.288032   0.288032   0.538032 geig sumcheck= -0.591148D+01-0.331698D-03
ikp isp nband q=   714     1    37     -0.163032   0.413032   0.413032 geig sumcheck= -0.677175D+00-0.532558D-04
ikp isp nband q=   715     1    37     -0.038032   0.538032   0.288032 geig sumcheck= -0.308824D+01 0.129584D-03
ikp isp nband q=   716     1    36      0.086968   0.663032   0.163032 geig sumcheck= -0.554569D+01 0.247940D-04
ikp isp nband q=   717     1    36      0.211968   0.788032   0.038032 geig sumcheck=  0.192944D+00-0.101540D-04
ikp isp nband q=   718     1    36      0.336968   0.913032  -0.086968 geig sumcheck=  0.221183D+01-0.227079D-04
ikp isp nband q=   719     1    37      0.461968   1.038032  -0.211968 geig sumcheck= -0.141794D+01-0.103095D-04
ikp isp nband q=   720     1    36      0.586968   1.163032  -0.336968 geig sumcheck=  0.945016D+01-0.182579D-04
ikp isp nband q=   721     1    36     -0.163032   0.163032   0.663032 geig sumcheck= -0.637882D+01-0.601419D-05
ikp isp nband q=   722     1    37     -0.038032   0.288032   0.538032 geig sumcheck= -0.358994D+01-0.532977D-04
ikp isp nband q=   723     1    37      0.086968   0.413032   0.413032 geig sumcheck= -0.317235D+01-0.595350D-05
ikp isp nband q=   724     1    38      0.211968   0.538032   0.288032 geig sumcheck= -0.462529D+01-0.371016D-04
ikp isp nband q=   725     1    37      0.336968   0.663032   0.163032 geig sumcheck= -0.253952D+01 0.149860D-03
ikp isp nband q=   726     1    37      0.461968   0.788032   0.038032 geig sumcheck=  0.650676D+00 0.470225D-05
ikp isp nband q=   727     1    37      0.586968   0.913032  -0.086968 geig sumcheck=  0.638723D+01 0.490552D-05
ikp isp nband q=   728     1    36      0.711968   1.038032  -0.211968 geig sumcheck=  0.100328D+01 0.866627D-05
ikp isp nband q=   729     1    36     -0.038032   0.038032   0.788032 geig sumcheck= -0.558752D+01 0.957047D-04
ikp isp nband q=   730     1    36      0.086968   0.163032   0.663032 geig sumcheck= -0.554410D+01 0.116855D-04
ikp isp nband q=   731     1    38      0.211968   0.288032   0.538032 geig sumcheck= -0.162135D+01-0.536779D-04
ikp isp nband q=   732     1    37      0.336968   0.413032   0.413032 geig sumcheck= -0.227761D+01 0.790832D-04
ikp isp nband q=   733     1    37      0.461968   0.538032   0.288032 geig sumcheck=  0.512964D+01-0.118055D-04
ikp isp nband q=   734     1    35      0.586968   0.663032   0.163032 geig sumcheck=  0.177808D+01-0.600423D-04
ikp isp nband q=   735     1    35      0.711968   0.788032   0.038032 geig sumcheck= -0.243468D+01-0.445373D-04
ikp isp nband q=   736     1    36      0.836968   0.913032  -0.086968 geig sumcheck=  0.312851D+01 0.217396D-04
ikp isp nband q=   737     1    36      0.086968  -0.086968   0.913032 geig sumcheck=  0.250556D+01 0.154785D-04
ikp isp nband q=   738     1    36      0.211968   0.038032   0.788032 geig sumcheck=  0.104865D+01-0.500333D-06
ikp isp nband q=   739     1    37      0.336968   0.163032   0.663032 geig sumcheck= -0.253955D+01-0.327094D-04
ikp isp nband q=   740     1    37      0.461968   0.288032   0.538032 geig sumcheck=  0.639232D+01-0.118047D-04
ikp isp nband q=   741     1    37      0.586968   0.413032   0.413032 geig sumcheck= -0.471559D+01 0.289986D-04
ikp isp nband q=   742     1    35      0.711968   0.538032   0.288032 geig sumcheck= -0.343652D+01 0.548313D-05
ikp isp nband q=   743     1    36      0.836968   0.663032   0.163032 geig sumcheck= -0.862677D+01 0.565510D-05
ikp isp nband q=   744     1    35      0.961968   0.788032   0.038032 geig sumcheck=  0.171214D+01-0.411562D-04
ikp isp nband q=   745     1    35      0.211968  -0.211968   1.038032 geig sumcheck= -0.139976D+01 0.126077D-03
ikp isp nband q=   746     1    36      0.336968  -0.086968   0.913032 geig sumcheck=  0.244271D+01 0.691777D-05
ikp isp nband q=   747     1    37      0.461968   0.038032   0.788032 geig sumcheck=  0.650553D+00-0.314354D-04
ikp isp nband q=   748     1    35      0.586968   0.163032   0.663032 geig sumcheck=  0.173898D+01 0.129691D-04
ikp isp nband q=   749     1    35      0.711968   0.288032   0.538032 geig sumcheck= -0.252242D-01 0.489393D-04
ikp isp nband q=   750     1    36      0.836968   0.413032   0.413032 geig sumcheck=  0.280323D+01 0.277983D-03
ikp isp nband q=   751     1    37      0.961968   0.538032   0.288032 geig sumcheck=  0.133468D+00-0.976221D-03
ikp isp nband q=   752     1    37      1.086968   0.663032   0.163032 geig sumcheck=  0.836163D+01-0.791274D-04
ikp isp nband q=   753     1    34      0.336968  -0.336968   1.163032 geig sumcheck=  0.447160D+01-0.387956D-05
ikp isp nband q=   754     1    37      0.461968  -0.211968   1.038032 geig sumcheck= -0.141774D+01 0.161832D-04
ikp isp nband q=   755     1    37      0.586968  -0.086968   0.913032 geig sumcheck=  0.712036D+01-0.255015D-04
ikp isp nband q=   756     1    35      0.711968   0.038032   0.788032 geig sumcheck= -0.107220D+01-0.981511D-05
ikp isp nband q=   757     1    36      0.836968   0.163032   0.663032 geig sumcheck= -0.118753D+01 0.502192D-05
ikp isp nband q=   758     1    37      0.961968   0.288032   0.538032 geig sumcheck=  0.134543D+00 0.613766D-03
ikp isp nband q=   759     1    35      1.086968   0.413032   0.413032 geig sumcheck=  0.470300D+01 0.673390D-04
ikp isp nband q=   760     1    35      1.211968   0.538032   0.288032 geig sumcheck= -0.599222D+01 0.138051D-03
ikp isp nband q=   761     1    37      0.461968  -0.461968   1.288032 geig sumcheck= -0.822595D+01 0.150779D-03
ikp isp nband q=   762     1    36      0.586968  -0.336968   1.163032 geig sumcheck=  0.941879D+01-0.420065D-05
ikp isp nband q=   763     1    36      0.711968  -0.211968   1.038032 geig sumcheck=  0.100382D+01-0.943802D-06
ikp isp nband q=   764     1    36      0.836968  -0.086968   0.913032 geig sumcheck=  0.292492D+01 0.185882D-04
ikp isp nband q=   765     1    35      0.961968   0.038032   0.788032 geig sumcheck=  0.520844D+00 0.108355D-04
ikp isp nband q=   766     1    37      1.086968   0.163032   0.663032 geig sumcheck= -0.492935D+00-0.866256D-04
ikp isp nband q=   767     1    35      1.211968   0.288032   0.538032 geig sumcheck= -0.148434D+01-0.688685D-04
ikp isp nband q=   768     1    37      1.336968   0.413032   0.413032 geig sumcheck=  0.187350D+01-0.102901D-05
ikp isp nband q=   769     1    37     -0.538032   0.538032   0.538032 geig sumcheck= -0.403265D+01-0.136533D-04
ikp isp nband q=   770     1    37     -0.413032   0.663032   0.413032 geig sumcheck= -0.498027D+01 0.155095D-04
ikp isp nband q=   771     1    35     -0.288032   0.788032   0.288032 geig sumcheck= -0.299588D+01-0.116032D-04
ikp isp nband q=   772     1    36     -0.163032   0.913032   0.163032 geig sumcheck=  0.414730D+01-0.793912D-04
ikp isp nband q=   773     1    35     -0.038032   1.038032   0.038032 geig sumcheck= -0.196422D+01-0.286368D-04
ikp isp nband q=   774     1    36      0.086968   1.163032  -0.086968 geig sumcheck=  0.334204D+00 0.130599D-03
ikp isp nband q=   775     1    36      0.211968   1.288032  -0.211968 geig sumcheck= -0.195340D+01 0.340095D-04
ikp isp nband q=   776     1    37      0.336968   1.413032  -0.336968 geig sumcheck=  0.429198D+01 0.332486D-04
ikp isp nband q=   777     1    37     -0.413032   0.413032   0.663032 geig sumcheck=  0.120416D+01 0.120814D-04
ikp isp nband q=   778     1    37     -0.288032   0.538032   0.538032 geig sumcheck= -0.538094D+01 0.799836D-04
ikp isp nband q=   779     1    36     -0.163032   0.663032   0.413032 geig sumcheck= -0.358091D+01-0.926571D-05
ikp isp nband q=   780     1    36     -0.038032   0.788032   0.288032 geig sumcheck= -0.779352D+01-0.171740D-04
ikp isp nband q=   781     1    36      0.086968   0.913032   0.163032 geig sumcheck= -0.361103D+01-0.422782D-04
ikp isp nband q=   782     1    35      0.211968   1.038032   0.038032 geig sumcheck= -0.431924D+00 0.887731D-05
ikp isp nband q=   783     1    37      0.336968   1.163032  -0.086968 geig sumcheck= -0.130933D+01 0.104059D-04
ikp isp nband q=   784     1    35      0.461968   1.288032  -0.211968 geig sumcheck=  0.693029D+01 0.193824D-04
ikp isp nband q=   785     1    35     -0.288032   0.288032   0.788032 geig sumcheck= -0.385170D+01 0.105827D-04
ikp isp nband q=   786     1    36     -0.163032   0.413032   0.663032 geig sumcheck= -0.351013D+01 0.130120D-04
ikp isp nband q=   787     1    37     -0.038032   0.538032   0.538032 geig sumcheck= -0.563927D+01-0.685984D-04
ikp isp nband q=   788     1    37      0.086968   0.663032   0.413032 geig sumcheck= -0.151774D+01-0.840518D-05
ikp isp nband q=   789     1    36      0.211968   0.788032   0.288032 geig sumcheck= -0.191622D+01 0.932779D-04
ikp isp nband q=   790     1    35      0.336968   0.913032   0.163032 geig sumcheck= -0.308825D+01 0.360059D-04
ikp isp nband q=   791     1    38      0.461968   1.038032   0.038032 geig sumcheck= -0.451504D+01-0.641458D-05
ikp isp nband q=   792     1    36      0.586968   1.163032  -0.086968 geig sumcheck=  0.832471D+01-0.536355D-03
ikp isp nband q=   793     1    36     -0.163032   0.163032   0.913032 geig sumcheck=  0.462789D+01 0.555362D-05
ikp isp nband q=   794     1    36     -0.038032   0.288032   0.788032 geig sumcheck= -0.745888D+01 0.108695D-04
ikp isp nband q=   795     1    37      0.086968   0.413032   0.663032 geig sumcheck= -0.151792D+01-0.689192D-05
ikp isp nband q=   796     1    36      0.211968   0.538032   0.538032 geig sumcheck= -0.969574D+00 0.596195D-07
ikp isp nband q=   797     1    36      0.336968   0.663032   0.413032 geig sumcheck= -0.446983D+01-0.693190D-04
ikp isp nband q=   798     1    35      0.461968   0.788032   0.288032 geig sumcheck=  0.594226D+01-0.112876D-03
ikp isp nband q=   799     1    37      0.586968   0.913032   0.163032 geig sumcheck=  0.637725D+01-0.469942D-04
ikp isp nband q=   800     1    35      0.711968   1.038032   0.038032 geig sumcheck=  0.118505D+01 0.136399D-04
ikp isp nband q=   801     1    35     -0.038032   0.038032   1.038032 geig sumcheck= -0.209092D+01 0.168575D-03
ikp isp nband q=   802     1    36      0.086968   0.163032   0.913032 geig sumcheck= -0.361102D+01-0.583277D-04
ikp isp nband q=   803     1    36      0.211968   0.288032   0.788032 geig sumcheck= -0.100372D+01 0.259530D-03
ikp isp nband q=   804     1    36      0.336968   0.413032   0.663032 geig sumcheck= -0.447028D+01 0.183475D-03
ikp isp nband q=   805     1    37      0.461968   0.538032   0.538032 geig sumcheck=  0.284008D+01 0.253671D-03
ikp isp nband q=   806     1    37      0.586968   0.663032   0.413032 geig sumcheck= -0.358622D+01 0.501544D-05
ikp isp nband q=   807     1    37      0.711968   0.788032   0.288032 geig sumcheck=  0.644358D+01 0.641410D-04
ikp isp nband q=   808     1    36      0.836968   0.913032   0.163032 geig sumcheck=  0.490575D+01-0.167439D-04
ikp isp nband q=   809     1    36      0.086968  -0.086968   1.163032 geig sumcheck=  0.334193D+00-0.950872D-04
ikp isp nband q=   810     1    35      0.211968   0.038032   1.038032 geig sumcheck=  0.730981D+00 0.152118D-04
ikp isp nband q=   811     1    35      0.336968   0.163032   0.913032 geig sumcheck= -0.316354D+01-0.113413D-04
ikp isp nband q=   812     1    35      0.461968   0.288032   0.788032 geig sumcheck=  0.450945D+01 0.164486D-03
ikp isp nband q=   813     1    37      0.586968   0.413032   0.663032 geig sumcheck= -0.482379D+01 0.980925D-05
ikp isp nband q=   814     1    37      0.711968   0.538032   0.538032 geig sumcheck=  0.222322D+01-0.106643D-03
ikp isp nband q=   815     1    37      0.836968   0.663032   0.413032 geig sumcheck=  0.890892D+00 0.576613D-04
ikp isp nband q=   816     1    36      0.961968   0.788032   0.288032 geig sumcheck=  0.453050D+01-0.251425D-04
ikp isp nband q=   817     1    36      0.211968  -0.211968   1.288032 geig sumcheck= -0.251048D+01-0.923219D-04
ikp isp nband q=   818     1    37      0.336968  -0.086968   1.163032 geig sumcheck= -0.130101D+01 0.542888D-05
ikp isp nband q=   819     1    38      0.461968   0.038032   1.038032 geig sumcheck= -0.439645D+01-0.102059D-04
ikp isp nband q=   820     1    37      0.586968   0.163032   0.913032 geig sumcheck=  0.637641D+01 0.411153D-05
ikp isp nband q=   821     1    37      0.711968   0.288032   0.788032 geig sumcheck=  0.495558D+00-0.428059D-04
ikp isp nband q=   822     1    37      0.836968   0.413032   0.663032 geig sumcheck=  0.641349D+00 0.235546D-04
ikp isp nband q=   823     1    39      0.961968   0.538032   0.538032 geig sumcheck=  0.455205D+01 0.534523D-05
ikp isp nband q=   824     1    37      1.086968   0.663032   0.413032 geig sumcheck=  0.350455D+01 0.680618D-04
ikp isp nband q=   825     1    37      0.336968  -0.336968   1.413032 geig sumcheck=  0.129205D+01 0.612697D-04
ikp isp nband q=   826     1    35      0.461968  -0.211968   1.288032 geig sumcheck=  0.693021D+01-0.407734D-04
ikp isp nband q=   827     1    36      0.586968  -0.086968   1.163032 geig sumcheck=  0.888539D+01-0.928182D-03
ikp isp nband q=   828     1    35      0.711968   0.038032   1.038032 geig sumcheck=  0.548672D+00-0.136113D-04
ikp isp nband q=   829     1    36      0.836968   0.163032   0.913032 geig sumcheck=  0.490540D+01 0.531543D-05
ikp isp nband q=   830     1    36      0.961968   0.288032   0.788032 geig sumcheck=  0.453037D+01 0.229555D-04
ikp isp nband q=   831     1    37      1.086968   0.413032   0.663032 geig sumcheck= -0.106008D+01 0.532814D-04
ikp isp nband q=   832     1    37      1.211968   0.538032   0.538032 geig sumcheck=  0.664831D+01 0.316231D-05
ikp isp nband q=   833     1    38     -0.663032   0.663032   0.663032 geig sumcheck=  0.404850D+01-0.210667D-04
ikp isp nband q=   834     1    37     -0.538032   0.788032   0.538032 geig sumcheck= -0.838567D+01-0.154598D-03
ikp isp nband q=   835     1    35     -0.413032   0.913032   0.413032 geig sumcheck=  0.470966D+01-0.138039D-04
ikp isp nband q=   836     1    35     -0.288032   1.038032   0.288032 geig sumcheck=  0.534098D+01 0.136222D-02
ikp isp nband q=   837     1    36     -0.163032   1.163032   0.163032 geig sumcheck= -0.271792D+01 0.123369D-04
ikp isp nband q=   838     1    36     -0.038032   1.288032   0.038032 geig sumcheck= -0.458651D+01-0.623523D-03
ikp isp nband q=   839     1    36      0.086968   1.413032  -0.086968 geig sumcheck= -0.493054D+01-0.170502D-04
ikp isp nband q=   840     1    39      0.211968   1.538032  -0.211968 geig sumcheck=  0.603928D+00-0.460656D-05
ikp isp nband q=   841     1    37     -0.538032   0.538032   0.788032 geig sumcheck= -0.840380D+01 0.142205D-03
ikp isp nband q=   842     1    37     -0.413032   0.663032   0.663032 geig sumcheck= -0.114585D+02 0.195081D-04
ikp isp nband q=   843     1    35     -0.288032   0.788032   0.538032 geig sumcheck= -0.400197D+01 0.275897D-04
ikp isp nband q=   844     1    36     -0.163032   0.913032   0.413032 geig sumcheck= -0.264240D+01-0.136781D-04
ikp isp nband q=   845     1    35     -0.038032   1.038032   0.288032 geig sumcheck= -0.223005D+01-0.126427D-04
ikp isp nband q=   846     1    36      0.086968   1.163032   0.163032 geig sumcheck= -0.247918D+01-0.197950D-05
ikp isp nband q=   847     1    36      0.211968   1.288032   0.038032 geig sumcheck= -0.182651D+01-0.632980D-05
ikp isp nband q=   848     1    37      0.336968   1.413032  -0.086968 geig sumcheck= -0.403733D+01-0.134745D-04
ikp isp nband q=   849     1    35     -0.413032   0.413032   0.913032 geig sumcheck=  0.309749D+01-0.393655D-04
ikp isp nband q=   850     1    35     -0.288032   0.538032   0.788032 geig sumcheck= -0.308920D+01 0.182644D-04
ikp isp nband q=   851     1    34     -0.163032   0.663032   0.663032 geig sumcheck=  0.447124D+01-0.971989D-05
ikp isp nband q=   852     1    37     -0.038032   0.788032   0.538032 geig sumcheck= -0.521213D+01-0.168573D-05
ikp isp nband q=   853     1    37      0.086968   0.913032   0.413032 geig sumcheck= -0.303155D+01 0.639300D-04
ikp isp nband q=   854     1    35      0.211968   1.038032   0.288032 geig sumcheck= -0.193666D+01 0.398057D-04
ikp isp nband q=   855     1    36      0.336968   1.163032   0.163032 geig sumcheck= -0.545862D+01-0.113702D-04
ikp isp nband q=   856     1    37      0.461968   1.288032   0.038032 geig sumcheck= -0.543359D+01-0.577538D-06
ikp isp nband q=   857     1    35     -0.288032   0.288032   1.038032 geig sumcheck=  0.534095D+01-0.418346D-03
ikp isp nband q=   858     1    36     -0.163032   0.413032   0.913032 geig sumcheck= -0.287986D+01-0.869045D-05
ikp isp nband q=   859     1    37     -0.038032   0.538032   0.788032 geig sumcheck= -0.521146D+01-0.189975D-04
ikp isp nband q=   860     1    34      0.086968   0.663032   0.663032 geig sumcheck=  0.338759D+01 0.600633D-05
ikp isp nband q=   861     1    35      0.211968   0.788032   0.538032 geig sumcheck= -0.577237D+01-0.264841D-02
ikp isp nband q=   862     1    35      0.336968   0.913032   0.413032 geig sumcheck=  0.678473D+01-0.627821D-04
ikp isp nband q=   863     1    36      0.461968   1.038032   0.288032 geig sumcheck=  0.470778D+01-0.220129D-04
ikp isp nband q=   864     1    36      0.586968   1.163032   0.163032 geig sumcheck=  0.621338D+00-0.222936D-04
ikp isp nband q=   865     1    36     -0.163032   0.163032   1.163032 geig sumcheck= -0.271805D+01 0.111462D-04
ikp isp nband q=   866     1    35     -0.038032   0.288032   1.038032 geig sumcheck= -0.282736D+01 0.721059D-05
ikp isp nband q=   867     1    37      0.086968   0.413032   0.913032 geig sumcheck= -0.336977D+01 0.337411D-05
ikp isp nband q=   868     1    35      0.211968   0.538032   0.788032 geig sumcheck= -0.606316D+01-0.539677D-02
ikp isp nband q=   869     1    36      0.336968   0.663032   0.663032 geig sumcheck= -0.118039D+01 0.115881D-03
ikp isp nband q=   870     1    36      0.461968   0.788032   0.538032 geig sumcheck=  0.189604D+01-0.878941D-04
ikp isp nband q=   871     1    37      0.586968   0.913032   0.413032 geig sumcheck=  0.363271D+01-0.787308D-04
ikp isp nband q=   872     1    37      0.711968   1.038032   0.288032 geig sumcheck=  0.247009D+01 0.400143D-03
ikp isp nband q=   873     1    36     -0.038032   0.038032   1.288032 geig sumcheck= -0.520047D+01 0.463156D-03
ikp isp nband q=   874     1    36      0.086968   0.163032   1.163032 geig sumcheck= -0.233221D+01-0.207950D-04
ikp isp nband q=   875     1    35      0.211968   0.288032   1.038032 geig sumcheck= -0.178393D+01-0.289132D-04
ikp isp nband q=   876     1    35      0.336968   0.413032   0.913032 geig sumcheck=  0.201256D+01-0.121137D-04
ikp isp nband q=   877     1    36      0.461968   0.538032   0.788032 geig sumcheck=  0.284030D+01 0.196795D-03
ikp isp nband q=   878     1    38      0.586968   0.663032   0.663032 geig sumcheck=  0.349596D+01 0.213960D-04
ikp isp nband q=   879     1    39      0.711968   0.788032   0.538032 geig sumcheck=  0.736452D+00-0.433260D-04
ikp isp nband q=   880     1    36      0.836968   0.913032   0.413032 geig sumcheck= -0.341764D+00 0.413346D-05
ikp isp nband q=   881     1    36      0.086968  -0.086968   1.413032 geig sumcheck= -0.534720D+01-0.777240D-06
ikp isp nband q=   882     1    36      0.211968   0.038032   1.288032 geig sumcheck= -0.122970D+01-0.292676D-05
ikp isp nband q=   883     1    36      0.336968   0.163032   1.163032 geig sumcheck= -0.468195D+01 0.310298D-04
ikp isp nband q=   884     1    36      0.461968   0.288032   1.038032 geig sumcheck=  0.536743D+01 0.172287D-04
ikp isp nband q=   885     1    37      0.586968   0.413032   0.913032 geig sumcheck=  0.363266D+01 0.518462D-04
ikp isp nband q=   886     1    39      0.711968   0.538032   0.788032 geig sumcheck=  0.713724D+00 0.506737D-05
ikp isp nband q=   887     1    37      0.836968   0.663032   0.663032 geig sumcheck=  0.269948D+01-0.162987D-04
ikp isp nband q=   888     1    37      0.961968   0.788032   0.538032 geig sumcheck= -0.179206D+01-0.549394D-04
ikp isp nband q=   889     1    39      0.211968  -0.211968   1.538032 geig sumcheck= -0.242461D+00 0.129675D-04
ikp isp nband q=   890     1    37      0.336968  -0.086968   1.413032 geig sumcheck= -0.457174D+01 0.859931D-05
ikp isp nband q=   891     1    37      0.461968   0.038032   1.288032 geig sumcheck= -0.543318D+01 0.132841D-05
ikp isp nband q=   892     1    36      0.586968   0.163032   1.163032 geig sumcheck=  0.910968D+00-0.112026D-04
ikp isp nband q=   893     1    37      0.711968   0.288032   1.038032 geig sumcheck=  0.246905D+01-0.195506D-03
ikp isp nband q=   894     1    36      0.836968   0.413032   0.913032 geig sumcheck= -0.160101D+00 0.258801D-05
ikp isp nband q=   895     1    37      0.961968   0.538032   0.788032 geig sumcheck= -0.178077D+01-0.205551D-04
ikp isp nband q=   896     1    39      1.086968   0.663032   0.663032 geig sumcheck=  0.566536D+01 0.427957D-05
ikp isp nband q=   897     1    38     -0.788032   0.788032   0.788032 geig sumcheck= -0.321260D+00-0.641455D-04
ikp isp nband q=   898     1    39     -0.663032   0.913032   0.663032 geig sumcheck= -0.309237D+01-0.931769D-03
ikp isp nband q=   899     1    39     -0.538032   1.038032   0.538032 geig sumcheck= -0.253817D+01-0.936818D-04
ikp isp nband q=   900     1    36     -0.413032   1.163032   0.413032 geig sumcheck=  0.349709D+01-0.654245D-05
ikp isp nband q=   901     1    35     -0.288032   1.288032   0.288032 geig sumcheck= -0.769546D+01-0.316880D-04
ikp isp nband q=   902     1    36     -0.163032   1.413032   0.163032 geig sumcheck= -0.764041D+01 0.980093D-05
ikp isp nband q=   903     1    36     -0.038032   1.538032   0.038032 geig sumcheck= -0.296292D+01 0.773233D-05
ikp isp nband q=   904     1    37      0.086968   1.663032  -0.086968 geig sumcheck= -0.348773D+01-0.194311D-04
ikp isp nband q=   905     1    39     -0.663032   0.663032   0.913032 geig sumcheck= -0.330786D+01-0.554969D-03
ikp isp nband q=   906     1    39     -0.538032   0.788032   0.788032 geig sumcheck=  0.282475D+01-0.574227D-04
ikp isp nband q=   907     1    37     -0.413032   0.913032   0.663032 geig sumcheck=  0.485140D+01 0.580572D-04
ikp isp nband q=   908     1    37     -0.288032   1.038032   0.538032 geig sumcheck=  0.253044D+01 0.544872D-05
ikp isp nband q=   909     1    36     -0.163032   1.163032   0.413032 geig sumcheck= -0.227476D+01-0.832463D-05
ikp isp nband q=   910     1    37     -0.038032   1.288032   0.288032 geig sumcheck= -0.580253D+01 0.201704D-05
ikp isp nband q=   911     1    36      0.086968   1.413032   0.163032 geig sumcheck=  0.271186D+00 0.169882D-04
ikp isp nband q=   912     1    37      0.211968   1.538032   0.038032 geig sumcheck=  0.434573D+00-0.275148D-04
ikp isp nband q=   913     1    39     -0.538032   0.538032   1.038032 geig sumcheck= -0.514459D+01 0.518998D-04
ikp isp nband q=   914     1    37     -0.413032   0.663032   0.913032 geig sumcheck=  0.484945D+01-0.147773D-04
ikp isp nband q=   915     1    36     -0.288032   0.788032   0.788032 geig sumcheck=  0.439519D+01-0.681886D-05
ikp isp nband q=   916     1    37     -0.163032   0.913032   0.663032 geig sumcheck=  0.491162D+01 0.164403D-04
ikp isp nband q=   917     1    38     -0.038032   1.038032   0.538032 geig sumcheck=  0.337556D+01 0.773550D-05
ikp isp nband q=   918     1    37      0.086968   1.163032   0.413032 geig sumcheck= -0.409758D+00 0.415403D-04
ikp isp nband q=   919     1    37      0.211968   1.288032   0.288032 geig sumcheck=  0.156475D+01-0.113467D-03
ikp isp nband q=   920     1    37      0.336968   1.413032   0.163032 geig sumcheck= -0.924822D+01-0.999699D-03
ikp isp nband q=   921     1    36     -0.413032   0.413032   1.163032 geig sumcheck=  0.188550D+00 0.647886D-05
ikp isp nband q=   922     1    37     -0.288032   0.538032   1.038032 geig sumcheck=  0.383331D+01-0.312672D-05
ikp isp nband q=   923     1    37     -0.163032   0.663032   0.913032 geig sumcheck=  0.491156D+01 0.103883D-04
ikp isp nband q=   924     1    35     -0.038032   0.788032   0.788032 geig sumcheck=  0.176632D+01-0.241525D-04
ikp isp nband q=   925     1    36      0.086968   0.913032   0.663032 geig sumcheck=  0.305495D+01-0.203845D-04
ikp isp nband q=   926     1    37      0.211968   1.038032   0.538032 geig sumcheck=  0.543697D+01 0.383506D-04
ikp isp nband q=   927     1    35      0.336968   1.163032   0.413032 geig sumcheck=  0.158228D+01-0.581418D-04
ikp isp nband q=   928     1    35      0.461968   1.288032   0.288032 geig sumcheck= -0.362313D+01-0.259586D-04
ikp isp nband q=   929     1    35     -0.288032   0.288032   1.288032 geig sumcheck= -0.769578D+01 0.511122D-05
ikp isp nband q=   930     1    36     -0.163032   0.413032   1.163032 geig sumcheck= -0.227543D+01-0.287638D-05
ikp isp nband q=   931     1    38     -0.038032   0.538032   1.038032 geig sumcheck=  0.365124D+01-0.352216D-04
ikp isp nband q=   932     1    36      0.086968   0.663032   0.913032 geig sumcheck=  0.208906D+01 0.517131D-04
ikp isp nband q=   933     1    36      0.211968   0.788032   0.788032 geig sumcheck=  0.429887D+01-0.655735D-04
ikp isp nband q=   934     1    37      0.336968   0.913032   0.663032 geig sumcheck=  0.187247D+01 0.344429D-03
ikp isp nband q=   935     1    37      0.461968   1.038032   0.538032 geig sumcheck=  0.294020D+01-0.524598D-04
ikp isp nband q=   936     1    36      0.586968   1.163032   0.413032 geig sumcheck=  0.427527D+01 0.105039D-04
ikp isp nband q=   937     1    36     -0.163032   0.163032   1.413032 geig sumcheck= -0.765143D+01-0.298152D-05
ikp isp nband q=   938     1    37     -0.038032   0.288032   1.288032 geig sumcheck= -0.612917D+01-0.455721D-05
ikp isp nband q=   939     1    37      0.086968   0.413032   1.163032 geig sumcheck= -0.409070D+00-0.682319D-04
ikp isp nband q=   940     1    37      0.211968   0.538032   1.038032 geig sumcheck=  0.521260D+01-0.452874D-04
ikp isp nband q=   941     1    37      0.336968   0.663032   0.913032 geig sumcheck=  0.110521D+01 0.710079D-03
ikp isp nband q=   942     1    36      0.461968   0.788032   0.788032 geig sumcheck=  0.316850D+01 0.109334D-04
ikp isp nband q=   943     1    38      0.586968   0.913032   0.663032 geig sumcheck=  0.200924D+01 0.383157D-04
ikp isp nband q=   944     1    37      0.711968   1.038032   0.538032 geig sumcheck=  0.100390D+01-0.975251D-05
ikp isp nband q=   945     1    36     -0.038032   0.038032   1.538032 geig sumcheck= -0.296210D+01-0.131200D-04
ikp isp nband q=   946     1    36      0.086968   0.163032   1.413032 geig sumcheck=  0.454198D+00 0.202003D-04
ikp isp nband q=   947     1    37      0.211968   0.288032   1.288032 geig sumcheck=  0.156583D+01 0.331316D-04
ikp isp nband q=   948     1    35      0.336968   0.413032   1.163032 geig sumcheck=  0.270228D+01 0.219576D-04
ikp isp nband q=   949     1    37      0.461968   0.538032   1.038032 geig sumcheck=  0.187533D+01 0.342108D-03
ikp isp nband q=   950     1    38      0.586968   0.663032   0.913032 geig sumcheck=  0.156676D+01-0.235873D-04
ikp isp nband q=   951     1    39      0.711968   0.788032   0.788032 geig sumcheck=  0.268567D+01 0.605653D-04
ikp isp nband q=   952     1    38      0.836968   0.913032   0.663032 geig sumcheck=  0.463730D+01 0.563658D-04
ikp isp nband q=   953     1    37      0.086968  -0.086968   1.663032 geig sumcheck= -0.340836D+01-0.637394D-05
ikp isp nband q=   954     1    37      0.211968   0.038032   1.538032 geig sumcheck=  0.434728D+00 0.158915D-04
ikp isp nband q=   955     1    37      0.336968   0.163032   1.413032 geig sumcheck= -0.949282D+01-0.519043D-03
ikp isp nband q=   956     1    35      0.461968   0.288032   1.288032 geig sumcheck= -0.343554D+01 0.450716D-05
ikp isp nband q=   957     1    36      0.586968   0.413032   1.163032 geig sumcheck=  0.427554D+01 0.107188D-04
ikp isp nband q=   958     1    37      0.711968   0.538032   1.038032 geig sumcheck=  0.999542D+00-0.592640D-05
ikp isp nband q=   959     1    38      0.836968   0.663032   0.913032 geig sumcheck=  0.427108D+01 0.815261D-05
ikp isp nband q=   960     1    39      0.961968   0.788032   0.788032 geig sumcheck=  0.104451D+01 0.199580D-04
ikp isp nband q=   961     1    38     -0.913032   0.913032   0.913032 geig sumcheck= -0.647639D+01-0.109503D-03
ikp isp nband q=   962     1    39     -0.788032   1.038032   0.788032 geig sumcheck=  0.241468D+00-0.479324D-05
ikp isp nband q=   963     1    37     -0.663032   1.163032   0.663032 geig sumcheck=  0.366369D+01-0.286130D-04
ikp isp nband q=   964     1    37     -0.538032   1.288032   0.538032 geig sumcheck=  0.222039D+01 0.125407D-03
ikp isp nband q=   965     1    37     -0.413032   1.413032   0.413032 geig sumcheck= -0.202272D+01 0.134235D-03
ikp isp nband q=   966     1    38     -0.288032   1.538032   0.288032 geig sumcheck=  0.503773D+01 0.199759D-04
ikp isp nband q=   967     1    37     -0.163032   1.663032   0.163032 geig sumcheck= -0.549263D+00-0.490318D-04
ikp isp nband q=   968     1    39     -0.038032   1.788032   0.038032 geig sumcheck=  0.639342D+00-0.840726D-04
ikp isp nband q=   969     1    39     -0.788032   0.788032   1.038032 geig sumcheck=  0.328037D+01 0.547377D-04
ikp isp nband q=   970     1    37     -0.663032   0.913032   0.913032 geig sumcheck= -0.725394D+00 0.134816D-06
ikp isp nband q=   971     1    37     -0.538032   1.038032   0.788032 geig sumcheck=  0.533861D+01-0.238149D-04
ikp isp nband q=   972     1    37     -0.413032   1.163032   0.663032 geig sumcheck= -0.320357D+01 0.262945D-03
ikp isp nband q=   973     1    35     -0.288032   1.288032   0.538032 geig sumcheck= -0.129195D+01-0.381690D-05
ikp isp nband q=   974     1    36     -0.163032   1.413032   0.413032 geig sumcheck=  0.848443D+00 0.101771D-05
ikp isp nband q=   975     1    37     -0.038032   1.538032   0.288032 geig sumcheck=  0.808067D+00-0.207078D-04
ikp isp nband q=   976     1    38      0.086968   1.663032   0.163032 geig sumcheck=  0.202809D+01-0.215732D-03
ikp isp nband q=   977     1    37     -0.663032   0.663032   1.163032 geig sumcheck=  0.370346D+01 0.227779D-04
ikp isp nband q=   978     1    37     -0.538032   0.788032   1.038032 geig sumcheck=  0.533927D+01 0.931815D-05
ikp isp nband q=   979     1    36     -0.413032   0.913032   0.913032 geig sumcheck= -0.266389D+01 0.205629D-05
ikp isp nband q=   980     1    36     -0.288032   1.038032   0.788032 geig sumcheck= -0.363309D+01 0.663704D-04
ikp isp nband q=   981     1    36     -0.163032   1.163032   0.663032 geig sumcheck= -0.984342D+00-0.666142D-04
ikp isp nband q=   982     1    36     -0.038032   1.288032   0.538032 geig sumcheck= -0.406761D+01-0.185866D-03
ikp isp nband q=   983     1    37      0.086968   1.413032   0.413032 geig sumcheck= -0.130966D+01-0.636089D-05
ikp isp nband q=   984     1    39      0.211968   1.538032   0.288032 geig sumcheck=  0.102967D+01-0.428644D-04
ikp isp nband q=   985     1    37     -0.538032   0.538032   1.288032 geig sumcheck=  0.222687D+01-0.825528D-04
ikp isp nband q=   986     1    37     -0.413032   0.663032   1.163032 geig sumcheck= -0.306805D+01-0.495566D-04
ikp isp nband q=   987     1    36     -0.288032   0.788032   1.038032 geig sumcheck= -0.342565D+01-0.160403D-03
ikp isp nband q=   988     1    36     -0.163032   0.913032   0.913032 geig sumcheck= -0.542692D+00-0.659134D-05
ikp isp nband q=   989     1    35     -0.038032   1.038032   0.788032 geig sumcheck=  0.297887D+01-0.105349D-04
ikp isp nband q=   990     1    35      0.086968   1.163032   0.663032 geig sumcheck=  0.158900D+01-0.762161D-04
ikp isp nband q=   991     1    35      0.211968   1.288032   0.538032 geig sumcheck=  0.513484D+01-0.153195D-03
ikp isp nband q=   992     1    37      0.336968   1.413032   0.413032 geig sumcheck=  0.697622D+01-0.358211D-04
ikp isp nband q=   993     1    37     -0.413032   0.413032   1.413032 geig sumcheck= -0.471309D+01 0.976516D-04
ikp isp nband q=   994     1    35     -0.288032   0.538032   1.288032 geig sumcheck= -0.129193D+01 0.657418D-06
ikp isp nband q=   995     1    36     -0.163032   0.663032   1.163032 geig sumcheck= -0.137767D+00 0.681439D-04
ikp isp nband q=   996     1    35     -0.038032   0.788032   1.038032 geig sumcheck=  0.373065D+01-0.190444D-04
ikp isp nband q=   997     1    36      0.086968   0.913032   0.913032 geig sumcheck=  0.258085D+01-0.384781D-04
ikp isp nband q=   998     1    36      0.211968   1.038032   0.788032 geig sumcheck=  0.471937D+01 0.222598D-04
ikp isp nband q=   999     1    37      0.336968   1.163032   0.663032 geig sumcheck=  0.687159D+00-0.720439D-05
ikp isp nband q=  1000     1    37      0.461968   1.288032   0.538032 geig sumcheck= -0.125580D+01-0.459579D-05
ikp isp nband q=  1001     1    38     -0.288032   0.288032   1.538032 geig sumcheck=  0.500044D+01-0.283640D-05
ikp isp nband q=  1002     1    36     -0.163032   0.413032   1.413032 geig sumcheck=  0.850573D+00 0.873644D-05
ikp isp nband q=  1003     1    36     -0.038032   0.538032   1.288032 geig sumcheck= -0.446768D+01 0.325626D-04
ikp isp nband q=  1004     1    35      0.086968   0.663032   1.163032 geig sumcheck=  0.158895D+01 0.371330D-04
ikp isp nband q=  1005     1    36      0.211968   0.788032   1.038032 geig sumcheck=  0.471920D+01 0.259386D-05
ikp isp nband q=  1006     1    36      0.336968   0.913032   0.913032 geig sumcheck= -0.222480D+01-0.122928D-04
ikp isp nband q=  1007     1    37      0.461968   1.038032   0.788032 geig sumcheck= -0.479226D+01 0.325653D-04
ikp isp nband q=  1008     1    38      0.586968   1.163032   0.663032 geig sumcheck=  0.190622D+01 0.253875D-04
ikp isp nband q=  1009     1    37     -0.163032   0.163032   1.663032 geig sumcheck= -0.549258D+00 0.797566D-04
ikp isp nband q=  1010     1    37     -0.038032   0.288032   1.538032 geig sumcheck=  0.795495D+00 0.114584D-03
ikp isp nband q=  1011     1    37      0.086968   0.413032   1.413032 geig sumcheck= -0.119489D+01-0.646014D-04
ikp isp nband q=  1012     1    35      0.211968   0.538032   1.288032 geig sumcheck=  0.594127D+01 0.189424D-03
ikp isp nband q=  1013     1    37      0.336968   0.663032   1.163032 geig sumcheck=  0.782122D+00 0.335121D-05
ikp isp nband q=  1014     1    37      0.461968   0.788032   1.038032 geig sumcheck= -0.479258D+01 0.366179D-05
ikp isp nband q=  1015     1    36      0.586968   0.913032   0.913032 geig sumcheck= -0.390130D+01-0.294646D-05
ikp isp nband q=  1016     1    39      0.711968   1.038032   0.788032 geig sumcheck= -0.241915D+01 0.782282D-04
ikp isp nband q=  1017     1    39     -0.038032   0.038032   1.788032 geig sumcheck= -0.107573D+00 0.108156D-03
ikp isp nband q=  1018     1    38      0.086968   0.163032   1.663032 geig sumcheck=  0.180983D+01 0.125294D-03
ikp isp nband q=  1019     1    39      0.211968   0.288032   1.538032 geig sumcheck=  0.483485D+00 0.734483D-04
ikp isp nband q=  1020     1    37      0.336968   0.413032   1.413032 geig sumcheck=  0.697564D+01 0.269060D-04
ikp isp nband q=  1021     1    37      0.461968   0.538032   1.288032 geig sumcheck= -0.112667D+01-0.670485D-05
ikp isp nband q=  1022     1    38      0.586968   0.663032   1.163032 geig sumcheck=  0.117905D+01 0.149431D-04
ikp isp nband q=  1023     1    39      0.711968   0.788032   1.038032 geig sumcheck= -0.419332D+01-0.789451D-03
ikp isp nband q=  1024     1    39      0.836968   0.913032   0.913032 geig sumcheck= -0.298117D+01-0.630313D-04

 File gwa: reading atom data
 site  1  z= 29.0  rmax= 2.31127  lmax=4  konf=***45
  l    g(rmax)    gp(rmax)    <g g>         <gp gp>         <g gp>
  0    1.23562   -0.28004   0.9999383616   0.0188583706   -0.000010
  gz(rmax); <gz gz> <gz g> <gz gp> =  0.00000   1.8024304347  -0.8959379  -0.1270565
  1    1.47221   -0.20414   0.9999270029   0.0108320679   -0.000002
  gz(rmax); <gz gz> <gz g> <gz gp> =  0.00000   1.6685444693  -0.7411808  -0.0960726
  2    0.35210   -1.70145   0.9997462348   1.0248517034   -0.000138
  gz(rmax); <gz gz> <gz g> <gz gp> =  0.00000   4.3182602498  -1.9736461  -0.6318814
  3    1.89461   -0.12452   0.9999373230   0.0031591671    0.000003
  4    2.00238   -0.10132   0.9999052814   0.0017751295    0.000003
  l  k isp       ecore      gc(rmax)     <gc gc>
  0  1  1    -648.944454    0.000000   0.9890902897  1  1
  0  2  1     -77.802710    0.000000   0.9979381838  2  2
  0  3  1      -8.245314    0.027613   0.9995422862  3  3
  1  2  1     -67.212192    0.000000   0.9979709590  4  4
  1  3  1      -5.145012    0.056291   0.9995949996  5  5
 
  *** DATA4GW_V2 size=           1           1           1         393
           5
                  =           4          29          44          59        1024
                  =          59           3

 ----  OK! lmf2gw: end --- DATA4GW_V2 is written 
