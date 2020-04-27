rdcmd:  lmfa --no-iactiv cu -vnk=8 -vbigbas=f
 -----------------------  START LMFA -----------------------

 rdctrl: reset global max nl from 5 to 4

 LMFA:     nbas = 1  nspec = 1  vn 7.11.c  verb 31,30
 special
 pot:      XC:BH
 autogen:  none

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
 MKSYM:  found 48 space group operations ... includes inversion

 Species A:  Z=29  Qc=18  R=2.311271  Q=0
 mesh:   rmt=2.311271  rmax=48.805862  a=0.025  nr=393  nr(rmax)=515
  Pl=  4.5     4.5     3.5     4.5    
  Ql=  1.0     0.0     10.0    0.0    

  iter     qint         drho          vh0          rho0          vsum     beta
    1   29.000000   4.725E+03      145.0000    0.1442E+03      -58.2772   0.30
   51   29.000000   4.211E-05      274.8263    0.2631E+05     -130.7916   0.30


 sumev=-4.333254  etot=-3304.416258  eref=-3304.434500  diff= 0.018242

 Free-atom wavefunctions:
 valence:      eval       node at      max at       c.t.p.   rho(r>rmt)
   4s      -0.36411         0.890       2.256       3.582     0.643062
   4p      -0.06295         0.975       3.484       7.414     0.901829
   3d      -0.39691         0.000       0.600       3.429     0.056076
   4f       0.01948         0.000      35.393      48.806*    1.000000

 core:        ecore       node at      max at       c.t.p.   rho(r>rmt)
   1s    -649.07634         0.000       0.034       0.069     0.000000
   2s     -77.91382         0.070       0.197       0.308     0.000000
   2p     -67.32532         0.000       0.158       0.335     0.000000
   3s      -8.39248         0.288       0.614       0.895     0.000141
   3p      -5.29682         0.260       0.619       1.078     0.000727

 Optimise free-atom basis for species A, rmt=2.311271
 l  it    Rsm      Eh     stiffR   stiffE      Eval      Exact     Pnu    Ql
 0  31   2.327  -0.301     107.2    415.5   -0.36393  -0.36411    4.76   1.00
 ... rsm exceeded rmt .. repeat with rsm=rmt
 0   1   2.311  -0.297     107.2    425.8   -0.36394  -0.36411    4.76   1.00
 1  12   2.674  -0.100      61.6    355.5   -0.05538  -0.06295    4.56   0.00
 ... rsm exceeded rmt .. repeat with rsm=rmt
 1   1   2.311  -0.100      61.6     45.6   -0.04993  -0.06295    4.56   0.00
 2  27   0.962  -0.116     158.7    107.6   -0.39670  -0.39691    3.89  10.00
 eigenvalue sum:  exact  -4.33325    opt basis  -4.33097    error 0.00228

 tailsm: fit tails to 6 smoothed hankels, rmt= 2.31127, rsm= 1.15564
 HNSMFT:  99 points in interval  2.31127  26.78517;  q=  1.203836
 E:    -1.00000    -2.00000    -4.00000    -6.00000    -9.00000    -15.0000
 C:    -0.07160    10.75053    -187.492    1222.024    -4717.79    21166.81
        r          rho         fit         diff
    2.311271    0.017797    0.017766    0.000031
    2.967767    0.005662    0.005658    0.000005
    3.810725    0.001517    0.001518   -0.000001
    4.893104    0.000305    0.000305    0.000000
    6.282906    0.000041    0.000041   -0.000001
    8.067448    0.000003    0.000003    0.000000
    q(fit):     1.203836    rms diff:   0.000016
    fit: r>rmt  1.203836   r<rmt  3.442816   qtot  4.646652
    rho: r>rmt  1.203836   r<rmt  9.796164   qtot 11.000000

 coretail: q=0.00392, rho(rmt)=0.00465.  Fit with Hankel e=-24.082  coeff=764.|
      r            rhoc          fit
    2.311271    0.02095279    0.02095279
    2.429779    0.01229068    0.01231367
    2.753317    0.00285262    0.00285190
    3.119934    0.00054243    0.00053465
    3.535366    0.00008235    0.00007888
    4.006112    0.00000969    0.00000887
    4.539536    0.00000085    0.00000073
    5.143985    0.00000005    0.00000004

 Sum of reference energies: -3304.4345
 Exit 0 LMFA 
 CPU time:    0.110s   Wall clock    0.147s   13:44:49 26.09.2014        on phpdl1.ph.kcl.ac.uk
rdcmd:  lmf  --no-iactiv cu -vnk=8 -vbigbas=f
 -----------------------  START LMF -----------------------

 rdctrl: reset global max nl from 5 to 4

 LMF:      nbas = 1  nspec = 1  vn 7.11.c  verb 31,30
 special
 pot:      XC:BH
 auto:     Pfloat:v6(CG v6)  Autoread:  none
 bz:       metal(5), tetra, invit 

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
 MKSYM:  found 48 space group operations ... includes inversion
 BZMESH:  60 irreducible QP from 512 ( 8 8 8 )  shift= T T T
 TETIRR: sorting 3072 tetrahedra ... 264 inequivalent ones found

 species data:  augmentation                           density
 spec       rmt   rsma lmxa kmxa      lmxl     rg   rsmv  kmxv foca   rfoca
 A        2.311  0.925    3    4         3  0.578  1.156    15    1   0.925

 MSHSIZ: mesh has 11 x 11 x 11 divisions; length 0.437, 0.437, 0.437
         generated from gmax = 9.0 a.u. : 941 vectors of 1331 (70%)

 GVLIST: gmax = 9 a.u. created 941 vectors of 1331 (70%)
         mesh has 11 x 11 x 11 divisions; length 0.437, 0.437, 0.437
 SGVSYM: 41 symmetry stars found for 941 reciprocal lattice vectors
 

 Makidx:  hamiltonian dimensions Low, Int, High, Negl: 18 0 7 7
 kappa   Low   Int   High  L+I  L+I+H  Neglected
   1       9     0     7     9    16       0
   2       9     0     0     9     9       7
  all     18     0     7    18    25       7
 suham :  16 augmentation channels, 16 local potential channels  Maximum lmxa=3

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol=1.0e-6
 spec      l    rsm    eh     gmax    last term   cutoff
  A        0*   2.50  -0.01   2.974    2.30E-05      27 
  A        1    2.50  -0.01   3.093    1.29E-06      51 
  A        2    1.00  -0.01   8.508    1.16E-06     869 
 

 lmfp  : no rst file ... try to overlap atomic densities
 rdovfa: read and overlap free-atom densities (mesh density) ...
 rdovfa: expected A,       read A        with rmt=  2.3113  mesh   393  0.025

 ovlpfa: overlap smooth part of FA densities
 site   1  spec  1  pos  0.0000  0.0000  0.0000  Qsmooth 4.646654
 total smooth Q = 4.646654

 Free atom and overlapped crystal site charges:
   ib    true(FA)    smooth(FA)  true(OV)    smooth(OV)    local
    1    9.796164    3.442818   10.275300    3.921954    6.353346

 Smooth charge on mesh:            4.646654
 Sum of local charges:             6.353346
 Total valence charge:            11.000000
 Sum of core charges:             18.000000
 Sum of nuclear charges:         -29.000000
 Homogeneous background:           0.000000
 Deviation from neutrality:        0.000000

 --- BNDFP:  begin iteration 1 of 12 ---

 avg es pot at rmt= 0.554993  avg sphere pot= 0.633521  vconst=-0.554993

 smooth rhoves     11.022236   charge     4.646654
 smooth rhoeps =   -3.843801   rhomu =   -5.010455  avg vxc =   -0.851784 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000013

 Energy terms:             smooth           local           total
   rhoval*vef            -12.156986      -177.336819      -189.493805
   rhoval*ves            -46.689416      -115.324371      -162.013788
   psnuc*ves              68.733888    -12976.662451    -12907.928563
   utot                   11.022236     -6545.993411     -6534.971175
   rho*exc                -3.843801      -126.414296      -130.258096
   rho*vxc                -5.010455      -167.409314      -172.419769
   valence chg             4.646654         6.353346        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.6627 -0.0591 -0.0536 -0.0536  0.0167  0.0167  1.7979  1.9471  1.9471
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.3283 -0.1257 -0.0727 -0.0255  0.0508  0.0920  0.6928  1.2801  1.4813
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4087 -0.1094 -0.0570 -0.0368  0.0286  0.0633  0.7989  1.4016  1.7830
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.3857 -0.1324 -0.0620  0.0030  0.0284  0.0660  0.9994  1.1878  1.4473
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.2227 -0.1325 -0.1100  0.0005  0.0467  0.2005  0.6490  0.8702  1.1936
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4328 -0.0995 -0.0481 -0.0481  0.0391  0.0391  0.7797  1.7311  1.7311

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.144577;  11.000000 electrons;  D(Ef):    6.229
         Sum occ. bands:   -0.8532412  incl. Bloechl correction:   -0.006586

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.6627 -0.0591 -0.0536 -0.0536  0.0167  0.0167  1.7979  1.9471  1.9471
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.3283 -0.1257 -0.0727 -0.0255  0.0508  0.0920  0.6928  1.2801  1.4813
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4087 -0.1094 -0.0570 -0.0368  0.0286  0.0633  0.7989  1.4016  1.7830
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.3857 -0.1324 -0.0620  0.0030  0.0284  0.0660  0.9994  1.1878  1.4473
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.2227 -0.1325 -0.1100  0.0005  0.0467  0.2005  0.6490  0.8702  1.1936
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4328 -0.0995 -0.0481 -0.0481  0.0391  0.0391  0.7797  1.7311  1.7311

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.144577;  11.000000 electrons;  D(Ef):    6.229
         Sum occ. bands:   -0.8532412  incl. Bloechl correction:   -0.006586

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1    9.927753    3.113493    6.814260

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.509191   -0.293352    4.650000    4.686057    4.500000    4.686057
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.531005   -0.280671    4.340000    4.364506    4.250000    4.364506
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.679694   -0.050852    3.870000    3.857413    3.147584    3.857413
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.023532   -0.620420    4.110000    4.109925    4.102416    4.110000

 Harris energy:
 sumev=       -0.853241  val*vef=    -189.493805   sumtv=     188.640564
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -130.258096     utot=   -6534.971175    ehar=   -3304.832068

 srhov:     -6.360828   -168.222981   -174.583809 sumev=   -0.853241   sumtv=  173.730568

 Kohn-Sham energy:
 sumtv=      173.730568  sumtc=      3171.756639   ekin=     3345.487208
 rhoep=     -128.641260   utot=     -6521.364241   ehks=    -3304.518293
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 0.  RMS DQ=3.93e-2
 charges:       old           new         screened      rms diff       lin mix
 smooth       4.646654      4.185740      4.185740      0.038395      4.185740
 site    1    6.353346      6.814260      6.814260      0.019256      6.814260
 AMIX: nmix=0 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=3.93e-2
 unscreened rms difference:  smooth  0.045849   local  0.019256
   screened rms difference:  smooth  0.038395   local  0.019256   tot  0.039302

 iors  : write restart file (binary, mesh density) 

   it  1  of 12    ehf=      -0.397568   ehk=      -0.083793
h nk=8 bigbas=0 ehf=-.397568 ehk=-.0837934

 --- BNDFP:  begin iteration 2 of 12 ---

 avg es pot at rmt= 0.633049  avg sphere pot= 0.653661  vconst=-0.633049

 smooth rhoves     12.286140   charge     4.185740
 smooth rhoeps =   -3.109206   rhomu =   -4.047533  avg vxc =   -0.858397 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000645

 Energy terms:             smooth           local           total
   rhoval*vef             -7.713657      -175.301724      -183.015382
   rhoval*ves            -48.957947      -107.937688      -156.895635
   psnuc*ves              73.530228    -12964.307484    -12890.777256
   utot                   12.286140     -6536.122586     -6523.836445
   rho*exc                -3.109206      -125.879113      -128.988318
   rho*vxc                -4.047533      -166.689186      -170.736719
   valence chg             4.185740         6.814260        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7748 -0.6825 -0.6802 -0.6802 -0.6430 -0.6430  1.6026  1.7558  1.7558
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.7179 -0.6847 -0.6654 -0.6541 -0.6178 -0.3403  0.4593  1.0915  1.2747
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.7184 -0.6784 -0.6695 -0.6492 -0.6265 -0.4571  0.5708  1.2078  1.5794
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.7146 -0.6988 -0.6593 -0.6431 -0.6231 -0.4289  0.7940  1.0031  1.2312
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.6997 -0.6980 -0.6768 -0.6419 -0.6231 -0.1362  0.4145  0.6649  0.9900
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.7269 -0.6758 -0.6758 -0.6287 -0.6287 -0.4938  0.5479  1.5355  1.5355

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.191164;  11.000000 electrons;  D(Ef):    3.354
         Sum occ. bands:   -7.0935371  incl. Bloechl correction:   -0.013294

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7748 -0.6825 -0.6802 -0.6802 -0.6430 -0.6430  1.6026  1.7558  1.7558
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.7179 -0.6847 -0.6654 -0.6541 -0.6178 -0.3403  0.4593  1.0915  1.2747
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.7184 -0.6784 -0.6695 -0.6492 -0.6265 -0.4571  0.5708  1.2078  1.5794
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.7146 -0.6988 -0.6593 -0.6431 -0.6231 -0.4289  0.7940  1.0031  1.2312
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.6997 -0.6980 -0.6768 -0.6419 -0.6231 -0.1362  0.4145  0.6649  0.9900
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.7269 -0.6758 -0.6758 -0.6287 -0.6287 -0.4938  0.5479  1.5355  1.5355

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.191164;  11.000000 electrons;  D(Ef):    3.354
         Sum occ. bands:   -7.0935371  incl. Bloechl correction:   -0.013294

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.394585    1.890767    8.503818

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.427349   -0.491461    4.686057    4.674018    4.500000    4.674018
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.243014   -0.600422    4.364506    4.316680    4.250000    4.316680
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.862620   -0.665547    3.857413    3.897967    3.147584    3.897967
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.010092   -0.867678    4.110000    4.106157    4.102416    4.110000

 Harris energy:
 sumev=       -7.093537  val*vef=    -183.015382   sumtv=     175.921845
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -128.988318     utot=   -6523.836445    ehar=   -3305.146280

 srhov:     -4.210479   -222.187857   -226.398336 sumev=   -7.093537   sumtv=  219.304799

 Kohn-Sham energy:
 sumtv=      219.304799  sumtc=      3171.756639   ekin=     3391.061438
 rhoep=     -132.676307   utot=     -6562.138377   ehks=    -3303.753246
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 1.  RMS DQ=9.95e-2  last it=3.93e-2
 mixrho: (warning) scr. and lin-mixed densities had 43 and 43 negative points
 charges:       old           new         screened      rms diff       lin mix
 smooth       4.185740      2.496182      2.496182      0.032721      2.496182
 site    1    6.814260      8.503818      8.503818      0.064980      8.503818
 AMIX: nmix=1 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=9.95e-2
   tj: 0.82137
 unscreened rms difference:  smooth  0.024210   local  0.064980
   screened rms difference:  smooth  0.032721   local  0.064980   tot  0.099460

 iors  : write restart file (binary, mesh density) 

   it  2  of 12    ehf=      -0.711780   ehk=       0.681254
 From last iter    ehf=      -0.397568   ehk=      -0.083793
 diffe(q)= -0.314212 (0.099460)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf=-.7117797 ehk=.6812538

 --- BNDFP:  begin iteration 3 of 12 ---

 avg es pot at rmt= 0.604858  avg sphere pot= 0.655616  vconst=-0.604858

 smooth rhoves     10.955727   charge     3.883937
 smooth rhoeps =   -2.817928   rhomu =   -3.667314  avg vxc =   -0.843493 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000597

 Energy terms:             smooth           local           total
   rhoval*vef             -6.578736      -181.158603      -187.737339
   rhoval*ves            -47.653666      -113.406090      -161.059756
   psnuc*ves              69.565121    -12968.692727    -12899.127607
   utot                   10.955727     -6541.049409     -6530.093681
   rho*exc                -2.817928      -126.698952      -129.516879
   rho*vxc                -3.667314      -167.770666      -171.437979
   valence chg             3.883937         7.116063        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7458 -0.4488 -0.4455 -0.4455 -0.3974 -0.3974  1.6619  1.8140  1.8140
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.5237 -0.4580 -0.4440 -0.4215 -0.3677 -0.2470  0.5247  1.1466  1.3363
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.5484 -0.4459 -0.4405 -0.4308 -0.3820 -0.3323  0.6355  1.2649  1.6407
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.5294 -0.4906 -0.4332 -0.4012 -0.3767 -0.3162  0.8536  1.0566  1.2960
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.4866 -0.4765 -0.4623 -0.4008 -0.3737 -0.0585  0.4795  0.7237  1.0494
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5671 -0.4405 -0.4405 -0.4021 -0.3802 -0.3802  0.6136  1.5940  1.5940

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.121493;  11.000000 electrons;  D(Ef):    3.480
         Sum occ. bands:   -4.7253938  incl. Bloechl correction:   -0.011775

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7458 -0.4488 -0.4455 -0.4455 -0.3974 -0.3974  1.6619  1.8140  1.8140
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.5237 -0.4580 -0.4440 -0.4215 -0.3677 -0.2470  0.5247  1.1466  1.3363
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.5484 -0.4459 -0.4405 -0.4308 -0.3820 -0.3323  0.6355  1.2649  1.6407
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.5294 -0.4906 -0.4332 -0.4012 -0.3767 -0.3162  0.8536  1.0566  1.2960
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.4866 -0.4765 -0.4623 -0.4008 -0.3737 -0.0585  0.4795  0.7237  1.0494
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5671 -0.4405 -0.4405 -0.4021 -0.3802 -0.3802  0.6136  1.5940  1.5940

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.121493;  11.000000 electrons;  D(Ef):    3.480
         Sum occ. bands:   -4.7253938  incl. Bloechl correction:   -0.011775

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.271871    2.259799    8.012072

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.442716   -0.475592    4.674018    4.660037    4.500000    4.660037
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.302348   -0.555889    4.316680    4.316753    4.250000    4.316753
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.919168   -0.432802    3.897967    3.884710    3.147584    3.884710
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.013659   -0.807474    4.110000    4.107121    4.102416    4.110000

 Harris energy:
 sumev=       -4.725394  val*vef=    -187.737339   sumtv=     183.011945
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.516879     utot=   -6530.093681    ehar=   -3304.841976

 srhov:     -4.840018   -202.925957   -207.765975 sumev=   -4.725394   sumtv=  203.040581

 Kohn-Sham energy:
 sumtv=      203.040581  sumtc=      3171.756639   ekin=     3374.797220
 rhoep=     -131.305905   utot=     -6547.996582   ehks=    -3304.505266
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 2.  RMS DQ=4.79e-2  last it=9.95e-2
 mixrho: (warning) scr. and lin-mixed densities had 13 and 13 negative points
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.883937      2.987929      2.987929      0.016765      2.987929
 site    1    7.116063      8.012072      8.012072      0.032665      8.012072
 AMIX: nmix=2 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=4.79e-2
   tj:-1.02979  -0.10943
 add q= -0.000001 to preserve neutrality
 unscreened rms difference:  smooth  0.012611   local  0.032665
   screened rms difference:  smooth  0.016765   local  0.032665   tot  0.047871

 iors  : write restart file (binary, mesh density) 

   it  3  of 12    ehf=      -0.407476   ehk=      -0.070766
 From last iter    ehf=      -0.711780   ehk=       0.681254
 diffe(q)=  0.304304 (0.047871)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf=-.407476 ehk=-.0707664

 --- BNDFP:  begin iteration 4 of 12 ---

 avg es pot at rmt= 0.546196  avg sphere pot= 0.668996  vconst=-0.546196

 smooth rhoves      8.690001   charge     3.363254
 smooth rhoeps =   -2.347120   rhomu =   -3.053146  avg vxc =   -0.813648 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000547

 Energy terms:             smooth           local           total
   rhoval*vef             -4.912319      -185.932024      -190.844343
   rhoval*ves            -44.716247      -118.780514      -163.496761
   psnuc*ves              62.096249    -12969.504228    -12907.407979
   utot                    8.690001     -6544.142371     -6535.452370
   rho*exc                -2.347120      -127.794795      -130.141915
   rho*vxc                -3.053146      -169.213396      -172.266542
   valence chg             3.363254         7.636746        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.6957 -0.1538 -0.1488 -0.1488 -0.0829 -0.0829  1.7545  1.9042  1.9042
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.3722 -0.2061 -0.1643 -0.1212 -0.0493 -0.0030  0.6398  1.2364  1.4347
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4464 -0.1924 -0.1508 -0.1323 -0.0706 -0.0368  0.7475  1.3573  1.7368
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4237 -0.2198 -0.1524 -0.0944 -0.0692 -0.0337  0.9528  1.1446  1.3993
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.2797 -0.2146 -0.1966 -0.0958 -0.0544  0.1209  0.5953  0.8229  1.1472
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4701 -0.1801 -0.1431 -0.1431 -0.0617 -0.0617  0.7275  1.6864  1.6864

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.061250;  11.000000 electrons;  D(Ef):    5.201
         Sum occ. bands:   -1.7914893  incl. Bloechl correction:   -0.007635

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.6957 -0.1538 -0.1488 -0.1488 -0.0829 -0.0829  1.7545  1.9042  1.9042
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.3722 -0.2061 -0.1643 -0.1212 -0.0493 -0.0030  0.6398  1.2364  1.4347
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4464 -0.1924 -0.1508 -0.1323 -0.0706 -0.0368  0.7475  1.3573  1.7368
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4237 -0.2198 -0.1524 -0.0944 -0.0692 -0.0337  0.9528  1.1446  1.3993
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.2797 -0.2146 -0.1966 -0.0958 -0.0544  0.1209  0.5953  0.8229  1.1472
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4701 -0.1801 -0.1431 -0.1431 -0.0617 -0.0617  0.7275  1.6864  1.6864

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.061250;  11.000000 electrons;  D(Ef):    5.201
         Sum occ. bands:   -1.7914893  incl. Bloechl correction:   -0.007635

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.009529    2.922781    7.086749

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.494192   -0.361921    4.660037    4.674980    4.500000    4.674980
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.470760   -0.367870    4.316753    4.349055    4.250000    4.349055
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.760297   -0.143600    3.884710    3.863226    3.147584    3.863226
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.021030   -0.678390    4.110000    4.109216    4.102416    4.110000

 Harris energy:
 sumev=       -1.791489  val*vef=    -190.844343   sumtv=     189.052854
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -130.141915     utot=   -6535.452370    ehar=   -3304.784791

 srhov:     -5.800181   -175.206493   -181.006674 sumev=   -1.791489   sumtv=  179.215184

 Kohn-Sham energy:
 sumtv=      179.215184  sumtc=      3171.756639   ekin=     3350.971824
 rhoep=     -129.167119   utot=     -6526.481990   ehks=    -3304.677285
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 3.  RMS DQ=2.56e-2  last it=4.79e-2
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.363254      3.913252      3.913252      0.009439      3.913252
 site    1    7.636746      7.086749      7.086749      0.018966      7.086749
 AMIX: nmix=3 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=2.56e-2
   tj: 0.75629  -0.24102  -0.00278
 unscreened rms difference:  smooth  0.007298   local  0.018966
   screened rms difference:  smooth  0.009439   local  0.018966   tot  0.025618

 iors  : write restart file (binary, mesh density) 

   it  4  of 12    ehf=      -0.350291   ehk=      -0.242785
 From last iter    ehf=      -0.407476   ehk=      -0.070766
 diffe(q)=  0.057185 (0.025618)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf=-.3502914 ehk=-.2427853

 --- BNDFP:  begin iteration 5 of 12 ---

 avg es pot at rmt= 0.570272  avg sphere pot= 0.662809  vconst=-0.570272

 smooth rhoves      9.529828   charge     3.554226
 smooth rhoeps =   -2.513732   rhomu =   -3.270414  avg vxc =   -0.825537 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000559

 Energy terms:             smooth           local           total
   rhoval*vef             -5.469646      -184.926390      -190.396036
   rhoval*ves            -45.935780      -117.327945      -163.263725
   psnuc*ves              64.995435    -12970.107689    -12905.112254
   utot                    9.529828     -6543.717817     -6534.187990
   rho*exc                -2.513732      -127.429323      -129.943055
   rho*vxc                -3.270414      -168.732686      -172.003099
   valence chg             3.554226         7.445774        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7157 -0.2544 -0.2500 -0.2500 -0.1906 -0.1906  1.7204  1.8711  1.8711
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4130 -0.2894 -0.2613 -0.2238 -0.1578 -0.0991  0.5953  1.2024  1.3978
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4754 -0.2785 -0.2507 -0.2344 -0.1773 -0.1436  0.7044  1.3225  1.7009
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4536 -0.3115 -0.2483 -0.1991 -0.1738 -0.1390  0.9150  1.1111  1.3604
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3401 -0.3013 -0.2876 -0.1997 -0.1636  0.0474  0.5503  0.7849  1.1101
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4983 -0.2618 -0.2445 -0.2445 -0.1708 -0.1708  0.6837  1.6520  1.6520

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.015643;  11.000000 electrons;  D(Ef):    4.289
         Sum occ. bands:   -2.7824056  incl. Bloechl correction:   -0.009104

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7157 -0.2544 -0.2500 -0.2500 -0.1906 -0.1906  1.7204  1.8711  1.8711
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4130 -0.2894 -0.2613 -0.2238 -0.1578 -0.0991  0.5953  1.2024  1.3978
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4754 -0.2785 -0.2507 -0.2344 -0.1773 -0.1436  0.7044  1.3225  1.7009
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4536 -0.3115 -0.2483 -0.1991 -0.1738 -0.1390  0.9150  1.1111  1.3604
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3401 -0.3013 -0.2876 -0.1997 -0.1636  0.0474  0.5503  0.7849  1.1101
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4983 -0.2618 -0.2445 -0.2445 -0.1708 -0.1708  0.6837  1.6520  1.6520

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.015643;  11.000000 electrons;  D(Ef):    4.289
         Sum occ. bands:   -2.7824056  incl. Bloechl correction:   -0.009104

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.110340    2.678016    7.432324

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.473743   -0.416733    4.674980    4.664690    4.500000    4.664690
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.401311   -0.449628    4.349055    4.333438    4.250000    4.333438
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.843566   -0.242563    3.863226    3.869824    3.147584    3.869824
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.018110   -0.730132    4.110000    4.108392    4.102416    4.110000

 Harris energy:
 sumev=       -2.782406  val*vef=    -190.396036   sumtv=     187.613630
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.943055     utot=   -6534.187990    ehar=   -3304.760775

 srhov:     -5.479090   -184.679091   -190.158181 sumev=   -2.782406   sumtv=  187.375776

 Kohn-Sham energy:
 sumtv=      187.375776  sumtc=      3171.756639   ekin=     3359.132415
 rhoep=     -129.914923   utot=     -6533.978167   ehks=    -3304.760675
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=6.26e-4  last it=2.56e-2
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.554226      3.567676      3.567676      0.000185      3.567676
 site    1    7.445774      7.432324      7.432324      0.000488      7.432324
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: nmix=2 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=6.26e-4
   tj:-0.04897  -0.01271
 unscreened rms difference:  smooth  0.000207   local  0.000488
   screened rms difference:  smooth  0.000185   local  0.000488   tot  0.000626

 iors  : write restart file (binary, mesh density) 

   it  5  of 12    ehf=      -0.326275   ehk=      -0.326175
 From last iter    ehf=      -0.350291   ehk=      -0.242785
 diffe(q)=  0.024016 (0.000626)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf=-.3262752 ehk=-.3261747

 --- BNDFP:  begin iteration 6 of 12 ---

 avg es pot at rmt= 0.571249  avg sphere pot= 0.662742  vconst=-0.571249

 smooth rhoves      9.555293   charge     3.558121
 smooth rhoeps =   -2.516440   rhomu =   -3.273936  avg vxc =   -0.825908 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000558

 Energy terms:             smooth           local           total
   rhoval*vef             -5.475161      -184.940441      -190.415602
   rhoval*ves            -45.976094      -117.315111      -163.291205
   psnuc*ves              65.086680    -12970.135899    -12905.049219
   utot                    9.555293     -6543.725505     -6534.170212
   rho*exc                -2.516440      -127.419427      -129.935867
   rho*vxc                -3.273936      -168.719654      -171.993591
   valence chg             3.558121         7.441879        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7168 -0.2588 -0.2545 -0.2545 -0.1953 -0.1953  1.7185  1.8692  1.8692
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4151 -0.2932 -0.2656 -0.2284 -0.1626 -0.1031  0.5933  1.2007  1.3959
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4770 -0.2823 -0.2552 -0.2389 -0.1820 -0.1483  0.7025  1.3208  1.6990
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4552 -0.3156 -0.2526 -0.2037 -0.1784 -0.1435  0.9133  1.1095  1.3585
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3431 -0.3053 -0.2917 -0.2043 -0.1684  0.0444  0.5483  0.7831  1.1083
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4998 -0.2654 -0.2490 -0.2490 -0.1755 -0.1755  0.6817  1.6501  1.6501

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018781;  11.000000 electrons;  D(Ef):    4.260
         Sum occ. bands:   -2.8268949  incl. Bloechl correction:   -0.009167

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7168 -0.2588 -0.2545 -0.2545 -0.1953 -0.1953  1.7185  1.8692  1.8692
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4151 -0.2932 -0.2656 -0.2284 -0.1626 -0.1031  0.5933  1.2007  1.3959
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4770 -0.2823 -0.2552 -0.2389 -0.1820 -0.1483  0.7025  1.3208  1.6990
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4552 -0.3156 -0.2526 -0.2037 -0.1784 -0.1435  0.9133  1.1095  1.3585
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3431 -0.3053 -0.2917 -0.2043 -0.1684  0.0444  0.5483  0.7831  1.1083
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4998 -0.2654 -0.2490 -0.2490 -0.1755 -0.1755  0.6817  1.6501  1.6501

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018781;  11.000000 electrons;  D(Ef):    4.260
         Sum occ. bands:   -2.8268949  incl. Bloechl correction:   -0.009167

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.114002    2.668894    7.445108

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.473038   -0.419031    4.664690    4.664369    4.500000    4.664369
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.398891   -0.452935    4.333438    4.332899    4.250000    4.332899
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.847249   -0.246808    3.869824    3.870282    3.147584    3.870282
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.018007   -0.732342    4.110000    4.108362    4.102416    4.110000

 Harris energy:
 sumev=       -2.826895  val*vef=    -190.415602   sumtv=     187.588708
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.935867     utot=   -6534.170212    ehar=   -3304.760732

 srhov:     -5.466678   -185.028398   -190.495076 sumev=   -2.826895   sumtv=  187.668181

 Kohn-Sham energy:
 sumtv=      187.668181  sumtc=      3171.756639   ekin=     3359.424821
 rhoep=     -129.942162   utot=     -6534.243386   ehks=    -3304.760728
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=1.96e-4  last it=6.26e-4
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.558121      3.554892      3.554892      0.000080      3.554892
 site    1    7.441879      7.445108      7.445108      0.000129      7.445108
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: nmix=2 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=1.96e-4
   tj:-0.32200   0.01739
 unscreened rms difference:  smooth  0.000068   local  0.000129
   screened rms difference:  smooth  0.000080   local  0.000129   tot  0.000196

 iors  : write restart file (binary, mesh density) 

   it  6  of 12    ehf=      -0.326232   ehk=      -0.326228
 From last iter    ehf=      -0.326275   ehk=      -0.326175
 diffe(q)=  0.000043 (0.000196)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf=-.3262319 ehk=-.3262275

 --- BNDFP:  begin iteration 7 of 12 ---

 avg es pot at rmt= 0.571271  avg sphere pot= 0.662731  vconst=-0.571271

 smooth rhoves      9.552955   charge     3.557009
 smooth rhoeps =   -2.515217   rhomu =   -3.272338  avg vxc =   -0.825885 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000558

 Energy terms:             smooth           local           total
   rhoval*vef             -5.469780      -184.968726      -190.438506
   rhoval*ves            -45.974710      -117.337780      -163.312490
   psnuc*ves              65.080619    -12970.158677    -12905.078058
   utot                    9.552955     -6543.748229     -6534.195274
   rho*exc                -2.515217      -127.422210      -129.937427
   rho*vxc                -3.272338      -168.723327      -171.995666
   valence chg             3.557009         7.442991        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7168 -0.2584 -0.2541 -0.2541 -0.1949 -0.1949  1.7186  1.8693  1.8693
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4150 -0.2929 -0.2653 -0.2280 -0.1622 -0.1028  0.5934  1.2008  1.3960
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4770 -0.2820 -0.2548 -0.2386 -0.1816 -0.1478  0.7026  1.3208  1.6990
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4552 -0.3153 -0.2522 -0.2033 -0.1780 -0.1431  0.9133  1.1096  1.3586
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3429 -0.3050 -0.2914 -0.2039 -0.1680  0.0446  0.5484  0.7832  1.1084
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4998 -0.2652 -0.2486 -0.2486 -0.1751 -0.1751  0.6818  1.6502  1.6502

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018541;  11.000000 electrons;  D(Ef):    4.263
         Sum occ. bands:   -2.8231822  incl. Bloechl correction:   -0.009160

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7168 -0.2584 -0.2541 -0.2541 -0.1949 -0.1949  1.7186  1.8693  1.8693
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4150 -0.2929 -0.2653 -0.2280 -0.1622 -0.1028  0.5934  1.2008  1.3960
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4770 -0.2820 -0.2548 -0.2386 -0.1816 -0.1478  0.7026  1.3208  1.6990
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4552 -0.3153 -0.2522 -0.2033 -0.1780 -0.1431  0.9133  1.1096  1.3586
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3429 -0.3050 -0.2914 -0.2039 -0.1680  0.0446  0.5484  0.7832  1.1084
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4998 -0.2652 -0.2486 -0.2486 -0.1751 -0.1751  0.6818  1.6502  1.6502

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018541;  11.000000 electrons;  D(Ef):    4.263
         Sum occ. bands:   -2.8231822  incl. Bloechl correction:   -0.009160

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.113495    2.670249    7.443246

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.473154   -0.418935    4.664369    4.664415    4.500000    4.664415
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.399213   -0.452706    4.332899    4.332962    4.250000    4.332962
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.846897   -0.246420    3.870282    3.870258    3.147584    3.870258
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.018022   -0.732234    4.110000    4.108365    4.102416    4.110000

 Harris energy:
 sumev=       -2.823182  val*vef=    -190.438506   sumtv=     187.615324
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.937427     utot=   -6534.195274    ehar=   -3304.760738

 srhov:     -5.469371   -184.967003   -190.436374 sumev=   -2.823182   sumtv=  187.613192

 Kohn-Sham energy:
 sumtv=      187.613192  sumtc=      3171.756639   ekin=     3359.369831
 rhoep=     -129.937465   utot=     -6534.193104   ehks=    -3304.760738
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=1.21e-5  last it=1.96e-4
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.557009      3.556754      3.556754      0.000007      3.556754
 site    1    7.442991      7.443246      7.443246      0.000011      7.443246
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: nmix=2 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=1.21e-5
   tj: 0.11729   0.03865
 unscreened rms difference:  smooth  0.000007   local  0.000011
   screened rms difference:  smooth  0.000007   local  0.000011   tot  0.000012

 iors  : write restart file (binary, mesh density) 

   it  7  of 12    ehf=      -0.326238   ehk=      -0.326238
 From last iter    ehf=      -0.326232   ehk=      -0.326228
 diffe(q)= -0.000006 (0.000012)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf=-.3262377 ehk=-.3262375

 --- BNDFP:  begin iteration 8 of 12 ---

 avg es pot at rmt= 0.571250  avg sphere pot= 0.662743  vconst=-0.571250

 smooth rhoves      9.552497   charge     3.556958
 smooth rhoeps =   -2.515190   rhomu =   -3.272303  avg vxc =   -0.825879 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000558

 Energy terms:             smooth           local           total
   rhoval*vef             -5.469801      -184.964841      -190.434643
   rhoval*ves            -45.973936      -117.334866      -163.308802
   psnuc*ves              65.078930    -12970.152855    -12905.073925
   utot                    9.552497     -6543.743861     -6534.191364
   rho*exc                -2.515190      -127.422065      -129.937255
   rho*vxc                -3.272303      -168.723134      -171.995436
   valence chg             3.556958         7.443042        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7168 -0.2585 -0.2542 -0.2542 -0.1949 -0.1949  1.7186  1.8693  1.8693
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4150 -0.2929 -0.2653 -0.2280 -0.1622 -0.1028  0.5934  1.2008  1.3960
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4769 -0.2820 -0.2548 -0.2386 -0.1816 -0.1479  0.7026  1.3208  1.6991
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4552 -0.3153 -0.2523 -0.2033 -0.1780 -0.1432  0.9133  1.1096  1.3586
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3429 -0.3050 -0.2914 -0.2039 -0.1680  0.0446  0.5484  0.7832  1.1084
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4998 -0.2652 -0.2487 -0.2487 -0.1751 -0.1751  0.6818  1.6502  1.6502

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018553;  11.000000 electrons;  D(Ef):    4.263
         Sum occ. bands:   -2.8234004  incl. Bloechl correction:   -0.009161

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7168 -0.2585 -0.2542 -0.2542 -0.1949 -0.1949  1.7186  1.8693  1.8693
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4150 -0.2929 -0.2653 -0.2280 -0.1622 -0.1028  0.5934  1.2008  1.3960
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4769 -0.2820 -0.2548 -0.2386 -0.1816 -0.1479  0.7026  1.3208  1.6991
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4552 -0.3153 -0.2523 -0.2033 -0.1780 -0.1432  0.9133  1.1096  1.3586
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3429 -0.3050 -0.2914 -0.2039 -0.1680  0.0446  0.5484  0.7832  1.1084
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4998 -0.2652 -0.2487 -0.2487 -0.1751 -0.1751  0.6818  1.6502  1.6502

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018553;  11.000000 electrons;  D(Ef):    4.263
         Sum occ. bands:   -2.8234004  incl. Bloechl correction:   -0.009161

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.113539    2.670124    7.443415

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.473143   -0.418934    4.664415    4.664411    4.500000    4.664411
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.399186   -0.452714    4.332962    4.332957    4.250000    4.332957
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.846937   -0.246444    3.870258    3.870260    3.147584    3.870260
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.018020   -0.732232    4.110000    4.108365    4.102416    4.110000

 Harris energy:
 sumev=       -2.823400  val*vef=    -190.434643   sumtv=     187.611243
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.937255     utot=   -6534.191364    ehar=   -3304.760737

 srhov:     -5.469010   -184.973208   -190.442219 sumev=   -2.823400   sumtv=  187.618818

 Kohn-Sham energy:
 sumtv=      187.618818  sumtc=      3171.756639   ekin=     3359.375458
 rhoep=     -129.937927   utot=     -6534.198268   ehks=    -3304.760737
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=1.91e-5  last it=1.21e-5
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.556958      3.556585      3.556585      0.000010      3.556585
 site    1    7.443042      7.443415      7.443415      0.000013      7.443415
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=1.91e-5
   tj: 0.72703
 unscreened rms difference:  smooth  0.000009   local  0.000013
   screened rms difference:  smooth  0.000010   local  0.000013   tot  0.000019

 iors  : write restart file (binary, mesh density) 

   it  8  of 12    ehf=      -0.326237   ehk=      -0.326237
 From last iter    ehf=      -0.326238   ehk=      -0.326238
 diffe(q)=  0.000001 (0.000019)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf=-.3262368 ehk=-.3262367

 --- BNDFP:  begin iteration 9 of 12 ---

 avg es pot at rmt= 0.571229  avg sphere pot= 0.662748  vconst=-0.571229

 smooth rhoves      9.551550   charge     3.556708
 smooth rhoeps =   -2.514953   rhomu =   -3.271994  avg vxc =   -0.825867 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000558

 Energy terms:             smooth           local           total
   rhoval*vef             -5.468916      -184.967633      -190.436548
   rhoval*ves            -45.972759      -117.337656      -163.310415
   psnuc*ves              65.075858    -12970.153681    -12905.077823
   utot                    9.551550     -6543.745669     -6534.194119
   rho*exc                -2.514953      -127.422577      -129.937530
   rho*vxc                -3.271994      -168.723807      -171.995801
   valence chg             3.556708         7.443292        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7168 -0.2584 -0.2540 -0.2540 -0.1948 -0.1948  1.7186  1.8693  1.8693
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4150 -0.2928 -0.2652 -0.2279 -0.1621 -0.1027  0.5934  1.2008  1.3960
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4769 -0.2819 -0.2547 -0.2385 -0.1815 -0.1477  0.7026  1.3208  1.6991
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4551 -0.3152 -0.2521 -0.2032 -0.1779 -0.1430  0.9134  1.1096  1.3586
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3428 -0.3049 -0.2913 -0.2038 -0.1679  0.0447  0.5484  0.7832  1.1084
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4997 -0.2651 -0.2485 -0.2485 -0.1750 -0.1750  0.6818  1.6502  1.6502

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018476;  11.000000 electrons;  D(Ef):    4.264
         Sum occ. bands:   -2.8222768  incl. Bloechl correction:   -0.009159

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7168 -0.2584 -0.2540 -0.2540 -0.1948 -0.1948  1.7186  1.8693  1.8693
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4150 -0.2928 -0.2652 -0.2279 -0.1621 -0.1027  0.5934  1.2008  1.3960
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4769 -0.2819 -0.2547 -0.2385 -0.1815 -0.1477  0.7026  1.3208  1.6991
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4551 -0.3152 -0.2521 -0.2032 -0.1779 -0.1430  0.9134  1.1096  1.3586
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3428 -0.3049 -0.2913 -0.2038 -0.1679  0.0447  0.5484  0.7832  1.1084
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4997 -0.2651 -0.2485 -0.2485 -0.1750 -0.1750  0.6818  1.6502  1.6502

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018476;  11.000000 electrons;  D(Ef):    4.264
         Sum occ. bands:   -2.8222768  incl. Bloechl correction:   -0.009159

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.113424    2.670415    7.443008

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.473167   -0.418885    4.664411    4.664420    4.500000    4.664420
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.399262   -0.452634    4.332957    4.332973    4.250000    4.332973
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.846851   -0.246333    3.870260    3.870251    3.147584    3.870251
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.018024   -0.732182    4.110000    4.108366    4.102416    4.110000

 Harris energy:
 sumev=       -2.822277  val*vef=    -190.436548   sumtv=     187.614272
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.937530     utot=   -6534.194119    ehar=   -3304.760738

 srhov:     -5.469450   -184.961109   -190.430559 sumev=   -2.822277   sumtv=  187.608283

 Kohn-Sham energy:
 sumtv=      187.608283  sumtc=      3171.756639   ekin=     3359.364922
 rhoep=     -129.936986   utot=     -6534.188674   ehks=    -3304.760738
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=1.48e-5  last it=1.91e-5
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.556708      3.556992      3.556992      0.000008      3.556992
 site    1    7.443292      7.443008      7.443008      0.000010      7.443008
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=1.48e-5
   tj: 0.43640
 unscreened rms difference:  smooth  0.000007   local  0.000010
   screened rms difference:  smooth  0.000008   local  0.000010   tot  0.000015

 iors  : write restart file (binary, mesh density) 

   it  9  of 12    ehf=      -0.326238   ehk=      -0.326238
 From last iter    ehf=      -0.326237   ehk=      -0.326237
 diffe(q)= -0.000001 (0.000015)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf=-.3262378 ehk=-.3262376

 --- BNDFP:  begin iteration 10 of 12 ---

 avg es pot at rmt= 0.571238  avg sphere pot= 0.662746  vconst=-0.571238

 smooth rhoves      9.551967   charge     3.556814
 smooth rhoeps =   -2.515051   rhomu =   -3.272121  avg vxc =   -0.825873 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000558

 Energy terms:             smooth           local           total
   rhoval*vef             -5.469276      -184.966461      -190.435738
   rhoval*ves            -45.973292      -117.336449      -163.309741
   psnuc*ves              65.077227    -12970.153272    -12905.076045
   utot                    9.551967     -6543.744860     -6534.192893
   rho*exc                -2.515051      -127.422352      -129.937402
   rho*vxc                -3.272121      -168.723511      -171.995632
   valence chg             3.556814         7.443186        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7168 -0.2584 -0.2541 -0.2541 -0.1948 -0.1948  1.7186  1.8693  1.8693
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4150 -0.2929 -0.2652 -0.2279 -0.1621 -0.1027  0.5934  1.2008  1.3960
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4769 -0.2820 -0.2548 -0.2385 -0.1815 -0.1478  0.7026  1.3208  1.6991
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4551 -0.3153 -0.2522 -0.2032 -0.1780 -0.1431  0.9134  1.1096  1.3586
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3429 -0.3049 -0.2913 -0.2039 -0.1680  0.0447  0.5484  0.7832  1.1084
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4997 -0.2651 -0.2486 -0.2486 -0.1751 -0.1751  0.6818  1.6502  1.6502

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018513;  11.000000 electrons;  D(Ef):    4.263
         Sum occ. bands:   -2.8228194  incl. Bloechl correction:   -0.009160

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7168 -0.2584 -0.2541 -0.2541 -0.1948 -0.1948  1.7186  1.8693  1.8693
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4150 -0.2929 -0.2652 -0.2279 -0.1621 -0.1027  0.5934  1.2008  1.3960
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4769 -0.2820 -0.2548 -0.2385 -0.1815 -0.1478  0.7026  1.3208  1.6991
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4551 -0.3153 -0.2522 -0.2032 -0.1780 -0.1431  0.9134  1.1096  1.3586
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3429 -0.3049 -0.2913 -0.2039 -0.1680  0.0447  0.5484  0.7832  1.1084
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4997 -0.2651 -0.2486 -0.2486 -0.1751 -0.1751  0.6818  1.6502  1.6502

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018513;  11.000000 electrons;  D(Ef):    4.263
         Sum occ. bands:   -2.8228194  incl. Bloechl correction:   -0.009160

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.113479    2.670276    7.443203

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.473156   -0.418909    4.664420    4.664416    4.500000    4.664416
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.399225   -0.452673    4.332973    4.332965    4.250000    4.332965
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.846891   -0.246387    3.870251    3.870255    3.147584    3.870255
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.018022   -0.732206    4.110000    4.108366    4.102416    4.110000

 Harris energy:
 sumev=       -2.822819  val*vef=    -190.435738   sumtv=     187.612918
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.937402     utot=   -6534.192893    ehar=   -3304.760737

 srhov:     -5.469239   -184.966824   -190.436063 sumev=   -2.822819   sumtv=  187.613243

 Kohn-Sham energy:
 sumtv=      187.613243  sumtc=      3171.756639   ekin=     3359.369883
 rhoep=     -129.937431   utot=     -6534.193189   ehks=    -3304.760737
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=8.35e-7  last it=1.48e-5
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.556814      3.556797      3.556797      0.000006      3.556797
 site    1    7.443186      7.443203      7.443203      0.000001      7.443203
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=8.35e-7
   tj: 0.05344
 unscreened rms difference:  smooth  0.000006   local  0.000001
   screened rms difference:  smooth  0.000006   local  0.000001   tot  0.000001

 iors  : write restart file (binary, mesh density) 

   it 10  of 12    ehf=      -0.326237   ehk=      -0.326237
 From last iter    ehf=      -0.326238   ehk=      -0.326238
 diffe(q)=  0.000000 (0.000001)    tol= 0.000010 (0.000010)   more=F
c nk=8 bigbas=0 ehf=-.3262374 ehk=-.3262372
 Exit 0 LMF 
 CPU time:    3.054s   Wall clock    3.102s   13:44:52 26.09.2014        on phpdl1.ph.kcl.ac.uk
rdcmd:  rm mixm.cu
rdcmd:  lmf  --no-iactiv cu -vnk=8 -vbigbas=t -vpwmode=0 -voveps=0d-7
 -----------------------  START LMF -----------------------

 LMF:      nbas = 1  nspec = 1  vn 7.11.c  verb 31,30
 special
 pot:      XC:BH
 auto:     Pfloat:v6(CG v6)  Autoread:  none
 bz:       metal(5), tetra, invit 

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
 MKSYM:  found 48 space group operations ... includes inversion
 BZMESH:  60 irreducible QP from 512 ( 8 8 8 )  shift= T T T
 TETIRR: sorting 3072 tetrahedra ... 264 inequivalent ones found

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
 

 iors  : read restart file (binary, mesh density) 
         use from  restart file: ef window, positions, pnu 
         ignore in restart file: *
         site   1, species A       : augmentation lmax changed from 3 to 4
         site   1, species A       : inflate local density from nlm= 16 to 25

 --- BNDFP:  begin iteration 1 of 12 ---

 avg es pot at rmt= 0.571238  avg sphere pot= 0.662746  vconst=-0.571238

 smooth rhoves      9.551943   charge     3.556807
 smooth rhoeps =   -2.515044   rhomu =   -3.272112  avg vxc =   -0.825872 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000558

 Energy terms:             smooth           local           total
   rhoval*vef             -5.469252      -184.966528      -190.435780
   rhoval*ves            -45.973262      -117.336515      -163.309777
   psnuc*ves              65.077149    -12970.153276    -12905.076128
   utot                    9.551943     -6543.744896     -6534.192953
   rho*exc                -2.515044      -127.422364      -129.937408
   rho*vxc                -3.272112      -168.723527      -171.995639
   valence chg             3.556807         7.443193        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7176 -0.2587 -0.2544 -0.2544 -0.1946 -0.1946  1.5859  1.7373  1.8052
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4154 -0.2932 -0.2656 -0.2280 -0.1622 -0.1029  0.5909  1.1772  1.3092
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4775 -0.2825 -0.2549 -0.2386 -0.1817 -0.1479  0.6977  1.2862  1.5948
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4557 -0.3158 -0.2525 -0.2034 -0.1781 -0.1433  0.9069  1.0961  1.2465
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3433 -0.3051 -0.2917 -0.2038 -0.1683  0.0443  0.5466  0.7740  1.0777
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5003 -0.2656 -0.2487 -0.2487 -0.1751 -0.1751  0.6763  1.6276  1.6276

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018791;  11.000000 electrons;  D(Ef):    4.268
         Sum occ. bands:   -2.8252965  incl. Bloechl correction:   -0.009153

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7176 -0.2587 -0.2544 -0.2544 -0.1946 -0.1946  1.5859  1.7373  1.8052
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4154 -0.2932 -0.2656 -0.2280 -0.1622 -0.1029  0.5909  1.1772  1.3092
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4775 -0.2825 -0.2549 -0.2386 -0.1817 -0.1479  0.6977  1.2862  1.5948
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4557 -0.3158 -0.2525 -0.2034 -0.1781 -0.1433  0.9069  1.0961  1.2465
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3433 -0.3051 -0.2917 -0.2038 -0.1683  0.0443  0.5466  0.7740  1.0777
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5003 -0.2656 -0.2487 -0.2487 -0.1751 -0.1751  0.6763  1.6276  1.6276

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018791;  11.000000 electrons;  D(Ef):    4.268
         Sum occ. bands:   -2.8252965  incl. Bloechl correction:   -0.009153

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.128159    2.833316    7.294843

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.480585   -0.533466    4.664416    4.612337    4.500000    4.612337
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.407337   -0.469584    4.332965    4.328213    4.250000    4.328213
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.971684   -0.255972    3.870255    3.858616    3.147584    3.858616
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.021180   -0.806340    4.110000    4.106719    4.102416    4.110000
 4     1    0.004395   -0.781219    5.100000    5.079490    5.077979    5.100000

 Harris energy:
 sumev=       -2.825296  val*vef=    -190.435780   sumtv=     187.610484
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.937408     utot=   -6534.192953    ehar=   -3304.763238

 srhov:     -5.929509   -184.355126   -190.284635 sumev=   -2.825296   sumtv=  187.459339

 Kohn-Sham energy:
 sumtv=      187.459339  sumtc=      3171.756639   ekin=     3359.215978
 rhoep=     -129.927442   utot=     -6534.051703   ehks=    -3304.763166
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 0.  RMS DQ=4.37e-3
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.556807      3.705157      3.705157      0.003469      3.705157
 site    1    7.443193      7.294843      7.294843      0.004825      7.294843
 AMIX: nmix=0 mmix=8  nelts=2178  beta=1  tm=5  rmsdel=4.37e-3
 unscreened rms difference:  smooth  0.003486   local  0.004825
   screened rms difference:  smooth  0.003469   local  0.004825   tot  0.004371

 iors  : write restart file (binary, mesh density) 

   it  1  of 12    ehf=      -0.328738   ehk=      -0.328666
i nk=8 bigbas=1 pwmode=0 oveps=0 ehf=-.3287376 ehk=-.3286662

 --- BNDFP:  begin iteration 2 of 12 ---

 avg es pot at rmt= 0.581159  avg sphere pot= 0.653091  vconst=-0.581159

 smooth rhoves     10.108817   charge     3.705157
 smooth rhoeps =   -2.657500   rhomu =   -3.458050  avg vxc =   -0.833024 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000127

 Energy terms:             smooth           local           total
   rhoval*vef             -6.009431      -184.304719      -190.314150
   rhoval*ves            -46.615738      -116.582621      -163.198359
   psnuc*ves              66.833372    -12971.740285    -12904.906913
   utot                   10.108817     -6544.161453     -6534.052636
   rho*exc                -2.657500      -127.269989      -129.927490
   rho*vxc                -3.458050      -168.524407      -171.982458
   valence chg             3.705157         7.294843        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7177 -0.2622 -0.2578 -0.2578 -0.1970 -0.1970  1.5859  1.7372  1.8045
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4165 -0.2957 -0.2683 -0.2309 -0.1654 -0.1056  0.5898  1.1765  1.3087
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4780 -0.2853 -0.2576 -0.2415 -0.1847 -0.1509  0.6967  1.2855  1.5944
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4562 -0.3188 -0.2550 -0.2066 -0.1810 -0.1463  0.9061  1.0955  1.2460
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3450 -0.3076 -0.2945 -0.2062 -0.1720  0.0424  0.5456  0.7733  1.0771
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5008 -0.2682 -0.2516 -0.2516 -0.1782 -0.1782  0.6752  1.6269  1.6269

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.020718;  11.000000 electrons;  D(Ef):    4.244
         Sum occ. bands:   -2.8538378  incl. Bloechl correction:   -0.009204

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7177 -0.2622 -0.2578 -0.2578 -0.1970 -0.1970  1.5859  1.7372  1.8045
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4165 -0.2957 -0.2683 -0.2309 -0.1654 -0.1056  0.5898  1.1765  1.3087
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4780 -0.2853 -0.2576 -0.2415 -0.1847 -0.1509  0.6967  1.2855  1.5944
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4562 -0.3188 -0.2550 -0.2066 -0.1810 -0.1463  0.9061  1.0955  1.2460
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3450 -0.3076 -0.2945 -0.2062 -0.1720  0.0424  0.5456  0.7733  1.0771
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5008 -0.2682 -0.2516 -0.2516 -0.1782 -0.1782  0.6752  1.6269  1.6269

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.020718;  11.000000 electrons;  D(Ef):    4.244
         Sum occ. bands:   -2.8538378  incl. Bloechl correction:   -0.009204

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.131737    2.825122    7.306615

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.479386   -0.535600    4.612337    4.611671    4.500000    4.611671
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.404853   -0.472672    4.328213    4.327645    4.250000    4.327645
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.969531   -0.258629    3.858616    3.858674    3.147584    3.858674
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.021023   -0.810008    4.110000    4.106662    4.102416    4.110000
 4     1    0.004375   -0.784074    5.100000    5.079474    5.077979    5.100000

 Harris energy:
 sumev=       -2.853838  val*vef=    -190.314150   sumtv=     187.460312
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.927490     utot=   -6534.052636    ehar=   -3304.763174

 srhov:     -5.994034   -184.685096   -190.679130 sumev=   -2.853838   sumtv=  187.825292

 Kohn-Sham energy:
 sumtv=      187.825292  sumtc=      3171.756639   ekin=     3359.581932
 rhoep=     -129.959153   utot=     -6534.385845   ehks=    -3304.763067
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 1.  RMS DQ=7.68e-4  last it=4.37e-3
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.705157      3.693384      3.693384      0.000208      3.693384
 site    1    7.294843      7.306615      7.306615      0.000501      7.306615
 AMIX: nmix=1 mmix=8  nelts=2178  beta=1  tm=5  rmsdel=7.68e-4
   tj: 0.08073
 unscreened rms difference:  smooth  0.000174   local  0.000501
   screened rms difference:  smooth  0.000208   local  0.000501   tot  0.000768

 iors  : write restart file (binary, mesh density) 

   it  2  of 12    ehf=      -0.328674   ehk=      -0.328567
 From last iter    ehf=      -0.328738   ehk=      -0.328666
 diffe(q)=  0.000064 (0.000768)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=1 pwmode=0 oveps=0 ehf=-.3286741 ehk=-.3285668

 --- BNDFP:  begin iteration 3 of 12 ---

 avg es pot at rmt= 0.579979  avg sphere pot= 0.653052  vconst=-0.579979

 smooth rhoves     10.061343   charge     3.694335
 smooth rhoeps =   -2.647839   rhomu =   -3.445448  avg vxc =   -0.832366 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000134

 Energy terms:             smooth           local           total
   rhoval*vef             -5.973908      -184.535417      -190.509324
   rhoval*ves            -46.555865      -116.812498      -163.368362
   psnuc*ves              66.678552    -12971.953200    -12905.274648
   utot                   10.061343     -6544.382849     -6534.321505
   rho*exc                -2.647839      -127.303293      -129.951132
   rho*vxc                -3.445448      -168.568371      -172.013819
   valence chg             3.694335         7.305665        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7163 -0.2523 -0.2479 -0.2479 -0.1865 -0.1865  1.5877  1.7387  1.8068
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4124 -0.2876 -0.2588 -0.2208 -0.1547 -0.0967  0.5937  1.1793  1.3115
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4754 -0.2769 -0.2478 -0.2315 -0.1743 -0.1405  0.7004  1.2883  1.5968
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4535 -0.3098 -0.2457 -0.1963 -0.1708 -0.1361  0.9093  1.0982  1.2490
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3388 -0.2991 -0.2856 -0.1960 -0.1612  0.0488  0.5495  0.7764  1.0801
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4983 -0.2604 -0.2416 -0.2416 -0.1675 -0.1675  0.6790  1.6297  1.6297

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.013999;  11.000000 electrons;  D(Ef):    4.317
         Sum occ. bands:   -2.7566441  incl. Bloechl correction:   -0.009049

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7163 -0.2523 -0.2479 -0.2479 -0.1865 -0.1865  1.5877  1.7387  1.8068
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4124 -0.2876 -0.2588 -0.2208 -0.1547 -0.0967  0.5937  1.1793  1.3115
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4754 -0.2769 -0.2478 -0.2315 -0.1743 -0.1405  0.7004  1.2883  1.5968
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4535 -0.3098 -0.2457 -0.1963 -0.1708 -0.1361  0.9093  1.0982  1.2490
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3388 -0.2991 -0.2856 -0.1960 -0.1612  0.0488  0.5495  0.7764  1.0801
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4983 -0.2604 -0.2416 -0.2416 -0.1675 -0.1675  0.6790  1.6297  1.6297

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.013999;  11.000000 electrons;  D(Ef):    4.317
         Sum occ. bands:   -2.7566441  incl. Bloechl correction:   -0.009049

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.121426    2.845332    7.276094

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.481247   -0.529902    4.611671    4.613322    4.500000    4.613322
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.411028   -0.464856    4.327645    4.329217    4.250000    4.329217
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.950044   -0.248929    3.858674    3.858072    3.147584    3.858072
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.021443   -0.799123    4.110000    4.106882    4.102416    4.110000
 4     1    0.004440   -0.772492    5.100000    5.079581    5.077979    5.100000

 Harris energy:
 sumev=       -2.756644  val*vef=    -190.509324   sumtv=     187.752680
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.951132     utot=   -6534.321505    ehar=   -3304.763318

 srhov:     -6.021940   -183.665468   -189.687408 sumev=   -2.756644   sumtv=  186.930764

 Kohn-Sham energy:
 sumtv=      186.930764  sumtc=      3171.756639   ekin=     3358.687403
 rhoep=     -129.879088   utot=     -6533.571081   ehks=    -3304.762765
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 2.  RMS DQ=1.80e-3  last it=7.68e-4
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.694335      3.723906      3.723906      0.000530      3.723906
 site    1    7.305665      7.276094      7.276094      0.001213      7.276094
 AMIX: nmix=2 mmix=8  nelts=2178  beta=1  tm=5  rmsdel=1.8e-3
   tj: 0.70183  -0.00507
 unscreened rms difference:  smooth  0.000390   local  0.001213
   screened rms difference:  smooth  0.000530   local  0.001213   tot  0.001804

 iors  : write restart file (binary, mesh density) 

   it  3  of 12    ehf=      -0.328818   ehk=      -0.328265
 From last iter    ehf=      -0.328674   ehk=      -0.328567
 diffe(q)= -0.000144 (0.001804)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=1 pwmode=0 oveps=0 ehf=-.3288179 ehk=-.328265

 --- BNDFP:  begin iteration 4 of 12 ---

 avg es pot at rmt= 0.580860  avg sphere pot= 0.653051  vconst=-0.580860

 smooth rhoves     10.096266   charge     3.702580
 smooth rhoeps =   -2.655417   rhomu =   -3.455336  avg vxc =   -0.832820 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000135

 Energy terms:             smooth           local           total
   rhoval*vef             -6.002045      -184.376183      -190.378229
   rhoval*ves            -46.599144      -116.655488      -163.254631
   psnuc*ves              66.791676    -12971.816332    -12905.024656
   utot                   10.096266     -6544.235910     -6534.139644
   rho*exc                -2.655417      -127.279393      -129.934810
   rho*vxc                -3.455336      -168.536831      -171.992167
   valence chg             3.702580         7.297420        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7173 -0.2591 -0.2547 -0.2547 -0.1937 -0.1937  1.5865  1.7376  1.8052
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4152 -0.2932 -0.2653 -0.2277 -0.1621 -0.1029  0.5910  1.1774  1.3096
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4772 -0.2826 -0.2546 -0.2384 -0.1815 -0.1477  0.6978  1.2864  1.5951
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4553 -0.3160 -0.2521 -0.2034 -0.1778 -0.1431  0.9070  1.0963  1.2470
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3431 -0.3049 -0.2917 -0.2031 -0.1686  0.0444  0.5468  0.7742  1.0780
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5000 -0.2658 -0.2485 -0.2485 -0.1749 -0.1749  0.6764  1.6277  1.6277

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018651;  11.000000 electrons;  D(Ef):    4.266
         Sum occ. bands:   -2.8235770  incl. Bloechl correction:   -0.009156

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7173 -0.2591 -0.2547 -0.2547 -0.1937 -0.1937  1.5865  1.7376  1.8052
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4152 -0.2932 -0.2653 -0.2277 -0.1621 -0.1029  0.5910  1.1774  1.3096
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4772 -0.2826 -0.2546 -0.2384 -0.1815 -0.1477  0.6978  1.2864  1.5951
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4553 -0.3160 -0.2521 -0.2034 -0.1778 -0.1431  0.9070  1.0963  1.2470
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3431 -0.3049 -0.2917 -0.2031 -0.1686  0.0444  0.5468  0.7742  1.0780
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5000 -0.2658 -0.2485 -0.2485 -0.1749 -0.1749  0.6764  1.6277  1.6277

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018651;  11.000000 electrons;  D(Ef):    4.266
         Sum occ. bands:   -2.8235770  incl. Bloechl correction:   -0.009156

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.128443    2.831457    7.296986

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.479978   -0.533700    4.613322    4.612251    4.500000    4.612251
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.406875   -0.470206    4.329217    4.328151    4.250000    4.328151
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.963081   -0.255629    3.858072    3.858466    3.147584    3.858466
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.021154   -0.806683    4.110000    4.106730    4.102416    4.110000
 4     1    0.004395   -0.780578    5.100000    5.079507    5.077979    5.100000

 Harris energy:
 sumev=       -2.823577  val*vef=    -190.378229   sumtv=     187.554652
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.934810     utot=   -6534.139644    ehar=   -3304.763162

 srhov:     -6.002999   -184.365812   -190.368812 sumev=   -2.823577   sumtv=  187.545235

 Kohn-Sham energy:
 sumtv=      187.545235  sumtc=      3171.756639   ekin=     3359.301874
 rhoep=     -129.934006   utot=     -6534.131030   ehks=    -3304.763162
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 3.  RMS DQ=2.20e-5  last it=1.80e-3
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.702580      3.703014      3.703014      0.000017      3.703014
 site    1    7.297420      7.296986      7.296986      0.000015      7.296986
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=2178  beta=1  tm=5  rmsdel=2.2e-5
   tj:-0.01226
 unscreened rms difference:  smooth  0.000021   local  0.000015
   screened rms difference:  smooth  0.000017   local  0.000015   tot  0.000022

 iors  : write restart file (binary, mesh density) 

   it  4  of 12    ehf=      -0.328662   ehk=      -0.328662
 From last iter    ehf=      -0.328818   ehk=      -0.328265
 diffe(q)=  0.000156 (0.000022)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=1 pwmode=0 oveps=0 ehf=-.3286622 ehk=-.3286621

 --- BNDFP:  begin iteration 5 of 12 ---

 avg es pot at rmt= 0.580888  avg sphere pot= 0.653045  vconst=-0.580888

 smooth rhoves     10.096995   charge     3.702758
 smooth rhoeps =   -2.655612   rhomu =   -3.455591  avg vxc =   -0.832823 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000135

 Energy terms:             smooth           local           total
   rhoval*vef             -6.002699      -184.374101      -190.376800
   rhoval*ves            -46.599988      -116.653393      -163.253381
   psnuc*ves              66.793978    -12971.815985    -12905.022007
   utot                   10.096995     -6544.234689     -6534.137694
   rho*exc                -2.655612      -127.279025      -129.934637
   rho*vxc                -3.455591      -168.536347      -171.991937
   valence chg             3.702758         7.297242        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7173 -0.2591 -0.2548 -0.2548 -0.1938 -0.1938  1.5865  1.7376  1.8051
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4153 -0.2932 -0.2654 -0.2278 -0.1621 -0.1029  0.5910  1.1774  1.3095
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4772 -0.2827 -0.2546 -0.2385 -0.1815 -0.1478  0.6978  1.2863  1.5951
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4554 -0.3160 -0.2522 -0.2035 -0.1779 -0.1432  0.9070  1.0963  1.2469
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3431 -0.3050 -0.2918 -0.2031 -0.1687  0.0443  0.5467  0.7742  1.0780
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5000 -0.2659 -0.2486 -0.2486 -0.1749 -0.1749  0.6763  1.6277  1.6277

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018702;  11.000000 electrons;  D(Ef):    4.265
         Sum occ. bands:   -2.8242702  incl. Bloechl correction:   -0.009157

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7173 -0.2591 -0.2548 -0.2548 -0.1938 -0.1938  1.5865  1.7376  1.8051
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4153 -0.2932 -0.2654 -0.2278 -0.1621 -0.1029  0.5910  1.1774  1.3095
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4772 -0.2827 -0.2546 -0.2385 -0.1815 -0.1478  0.6978  1.2863  1.5951
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4554 -0.3160 -0.2522 -0.2035 -0.1779 -0.1432  0.9070  1.0963  1.2469
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3431 -0.3050 -0.2918 -0.2031 -0.1687  0.0443  0.5467  0.7742  1.0780
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5000 -0.2659 -0.2486 -0.2486 -0.1749 -0.1749  0.6763  1.6277  1.6277

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018702;  11.000000 electrons;  D(Ef):    4.265
         Sum occ. bands:   -2.8242702  incl. Bloechl correction:   -0.009157

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.128504    2.831326    7.297178

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.479955   -0.533756    4.612251    4.612234    4.500000    4.612234
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.406826   -0.470263    4.328151    4.328141    4.250000    4.328141
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.963382   -0.255690    3.858466    3.858479    3.147584    3.858479
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.021151   -0.806765    4.110000    4.106729    4.102416    4.110000
 4     1    0.004395   -0.780659    5.100000    5.079506    5.077979    5.100000

 Harris energy:
 sumev=       -2.824270  val*vef=    -190.376800   sumtv=     187.552530
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.934637     utot=   -6534.137694    ehar=   -3304.763162

 srhov:     -6.002879   -184.372689   -190.375568 sumev=   -2.824270   sumtv=  187.551298

 Kohn-Sham energy:
 sumtv=      187.551298  sumtc=      3171.756639   ekin=     3359.307938
 rhoep=     -129.934540   utot=     -6534.136560   ehks=    -3304.763162
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=3.06e-6  last it=2.20e-5
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.702758      3.702822      3.702822      0.000008      3.702822
 site    1    7.297242      7.297178      7.297178      0.000002      7.297178
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=2178  beta=1  tm=5  rmsdel=3.06e-6
   tj:-0.16043
 unscreened rms difference:  smooth  0.000009   local  0.000002
   screened rms difference:  smooth  0.000008   local  0.000002   tot  0.000003

 iors  : write restart file (binary, mesh density) 

   it  5  of 12    ehf=      -0.328662   ehk=      -0.328662
 From last iter    ehf=      -0.328662   ehk=      -0.328662
 diffe(q)=  0.000000 (0.000003)    tol= 0.000010 (0.000010)   more=F
c nk=8 bigbas=1 pwmode=0 oveps=0 ehf=-.3286618 ehk=-.3286617
 Exit 0 LMF 
 CPU time:    4.336s   Wall clock    4.373s   13:44:57 26.09.2014        on phpdl1.ph.kcl.ac.uk
rdcmd:  lmf  --no-iactiv cu -vnk=8 -vbigbas=t -vpwmode=11 -voveps=0d-7 --band:fn=syml
 -----------------------  START LMF -----------------------

 LMF:      nbas = 1  nspec = 1  vn 7.11.c  verb 31,30
 special:  APW basis(q)
 pot:      XC:BH
 auto:     Pfloat:v6(CG v6)  Autoread:  none
 bz:       metal(5), tetra, invit 

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
 MKSYM:  found 48 space group operations ... includes inversion
 BZMESH:  60 irreducible QP from 512 ( 8 8 8 )  shift= T T T
 TETIRR: sorting 3072 tetrahedra ... 264 inequivalent ones found

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

 --- BNDFP:  begin iteration 1 of 12 ---

 avg es pot at rmt= 0.580894  avg sphere pot= 0.653044  vconst=-0.580894

 smooth rhoves     10.097116   charge     3.702792
 smooth rhoeps =   -2.655652   rhomu =   -3.455643  avg vxc =   -0.832822 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000135

 Energy terms:             smooth           local           total
   rhoval*vef             -6.002841      -184.373715      -190.376556
   rhoval*ves            -46.600123      -116.653034      -163.253156
   psnuc*ves              66.794355    -12971.816011    -12905.021656
   utot                   10.097116     -6544.234522     -6534.137406
   rho*exc                -2.655652      -127.278965      -129.934618
   rho*vxc                -3.455643      -168.536268      -171.991911
   valence chg             3.702792         7.297208        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read efermi from weights file : ef = -0.018702
 suqlst: read qp for bands, mode 1

 suqlst:  nq= 41   q1= 0.5000 0.5000 0.5000   q2= 0.0000 0.0000 0.0000
 bndfp:  kpt 1 of 41, k=  0.50000  0.50000  0.50000   ndimh = 37
 -0.4126 -0.2589 -0.2589 -0.1500 -0.1500 -0.1062  0.2486  1.5424  1.5972
 bndfp:  kpt 11 of 41, k=  0.37500  0.37500  0.37500   ndimh = 38
 -0.4452 -0.2536 -0.2536 -0.2256 -0.1625 -0.1625  0.4808  1.5736  1.6092
 bndfp:  kpt 21 of 41, k=  0.25000  0.25000  0.25000   ndimh = 38
 -0.5709 -0.2761 -0.2447 -0.2447 -0.1870 -0.1870  0.8926  1.6452  1.6562
 bndfp:  kpt 31 of 41, k=  0.12500  0.12500  0.12500   ndimh = 38
 -0.6864 -0.2652 -0.2498 -0.2498 -0.1953 -0.1953  1.3643  1.6982  1.7475
 bndfp:  kpt 41 of 41, k=  0.00000  0.00000  0.00000   ndimh = 39
 -0.7283 -0.2569 -0.2569 -0.2569 -0.1930 -0.1930  1.6531  1.8382  1.8382

 Read efermi from weights file : ef = -0.018702

 suqlst:  nq= 41   q1= 0.0000 0.0000 0.0000   q2= 1.0000 0.0000 0.0000
 bndfp:  kpt 1 of 41, k=  0.00000  0.00000  0.00000   ndimh = 39
 -0.7283 -0.2569 -0.2569 -0.2569 -0.1930 -0.1930  1.6531  1.8382  1.8382
 bndfp:  kpt 11 of 41, k=  0.25000  0.00000  0.00000   ndimh = 36
 -0.6726 -0.2732 -0.2434 -0.2434 -0.2075 -0.1870  1.5219  1.5425  1.5425
 bndfp:  kpt 21 of 41, k=  0.50000  0.00000  0.00000   ndimh = 36
 -0.5231 -0.3114 -0.2284 -0.2064 -0.2064 -0.1723  1.2327  1.2327  1.3332
 bndfp:  kpt 31 of 41, k=  0.75000  0.00000  0.00000   ndimh = 36
 -0.4049 -0.3480 -0.1614 -0.1614 -0.1572 -0.1235  0.8173  1.0143  1.0143
 bndfp:  kpt 41 of 41, k=  1.00000  0.00000  0.00000   ndimh = 35
 -0.3961 -0.3628 -0.1509 -0.1394 -0.1394  0.0772  0.4995  0.9329  0.9329

 Read efermi from weights file : ef = -0.018702

 suqlst:  nq= 21   q1= 1.0000 0.0000 0.0000   q2= 1.0000 0.5000 0.0000
 bndfp:  kpt 1 of 21, k=  1.00000  0.00000  0.00000   ndimh = 35
 -0.3961 -0.3628 -0.1509 -0.1394 -0.1394  0.0772  0.4995  0.9329  0.9329
 bndfp:  kpt 11 of 21, k=  1.00000  0.25000  0.00000   ndimh = 35
 -0.3719 -0.3463 -0.2102 -0.1698 -0.1394  0.1977  0.5208  0.7103  0.7753
 bndfp:  kpt 21 of 21, k=  1.00000  0.50000  0.00000   ndimh = 39
 -0.3401 -0.2954 -0.2954 -0.1964 -0.1394  0.4288  0.4288  0.5908  0.6156

 Read efermi from weights file : ef = -0.018702

 suqlst:  nq= 41   q1= 1.0000 0.5000 0.0000   q2= 0.0000 0.0000 0.0000
 bndfp:  kpt 1 of 41, k=  1.00000  0.50000  0.00000   ndimh = 39
 -0.3401 -0.2954 -0.2954 -0.1964 -0.1394  0.4288  0.4288  0.5908  0.6156
 bndfp:  kpt 11 of 41, k=  0.75000  0.37500  0.00000   ndimh = 37
 -0.3462 -0.3172 -0.2822 -0.1887 -0.1617  0.0621  0.6340  0.7292  0.9122
 bndfp:  kpt 21 of 41, k=  0.50000  0.25000  0.00000   ndimh = 37
 -0.4794 -0.3054 -0.2569 -0.2057 -0.2019 -0.1413  0.9569  1.0052  1.3974
 bndfp:  kpt 31 of 41, k=  0.25000  0.12500  0.00000   ndimh = 38
 -0.6591 -0.2738 -0.2436 -0.2428 -0.2061 -0.1870  1.3631  1.4195  1.7139
 bndfp:  kpt 41 of 41, k=  0.00000  0.00000  0.00000   ndimh = 39
 -0.7283 -0.2569 -0.2569 -0.2569 -0.1930 -0.1930  1.6531  1.8382  1.8382

 Read efermi from weights file : ef = -0.018702
 Exit 0 bndfp 
 CPU time:    1.655s   Wall clock    1.701s   13:44:58 26.09.2014        on phpdl1.ph.kcl.ac.uk
rdcmd:  mv bnds.cu bnds.cu-pwmode11
