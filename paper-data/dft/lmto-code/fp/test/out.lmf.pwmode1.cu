rdcmd:  lmfa --no-iactiv cu -vnk=8 -vbigbas=f
 -----------------------  START LMFA (80000K)  -----------------------

 rdctrl: reset global max nl from 5 to 4

 LMFA:     nbas = 1  nspec = 1  vn 7.00(LMFA 7.0)  verb 31,30
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
   51   29.000000   4.211E-05      274.8263    0.2631E+05     -130.7915   0.30


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
 E:    -1.00000    -2.00000    -4.00000    -6.00000    -9.00000    -15.0000
 C:    -0.07160    10.75053    -187.492    1222.023    -4717.79    21166.81
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
 CPU time:    0.541s     Sat Jun 20 14:49:46 2009   on waldo.eas.asu.edu
 wkinfo:  used    94 K  workspace of 80000 K   in   0 K calls
rdcmd:  lmf  --no-iactiv cu -vnk=8 -vbigbas=f
 -----------------------  START LMF (80000K)  -----------------------

 rdctrl: reset global max nl from 5 to 4

 LMF:      nbas = 1  nspec = 1  vn 7.00(LMF 7.0)  verb 31,30
 special:  APW basis
 pot:      XC:BH
 auto:     Pfloat:v6(CG v6)  Autoread:  none
 bz:       metal(2), tetra, invit 

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
 spec      l    rsm    eh     gmax    last term    cutoff
  A        0*   2.50  -0.01   2.974    2.30E-05      27 
  A        1    2.50  -0.01   3.093    1.29E-06      51 
  A        2    1.00  -0.01   8.508    1.16E-06     869 
 suham:  PW basis with  1 < E < 3  =>  npw = 8,  ndham = 26

 iors  : read restart file (binary, mesh density) 
 iors  : empty file ... nothing read

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

 smooth rhoves     11.022237   charge     4.646654
 smooth rhoeps =   -3.843801   rhomu =   -5.010456  avg vxc =   -0.851784 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000012

 Energy terms:             smooth           local           total
   rhoval*vef            -12.156987      -177.336818      -189.493805
   rhoval*ves            -46.689417      -115.324370      -162.013788
   psnuc*ves              68.733891    -12976.662453    -12907.928563
   utot                   11.022237     -6545.993412     -6534.971175
   rho*exc                -3.843801      -126.414296      -130.258096
   rho*vxc                -5.010456      -167.409313      -172.419769
   valence chg             4.646654         6.353346        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Incompatible or missing qp weights file ...

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.6635 -0.0593 -0.0537 -0.0537  0.0167  0.0167  1.6484  1.7917  1.8747
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.3287 -0.1260 -0.0735 -0.0257  0.0506  0.0918  0.6905  1.2552  1.3832
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4093 -0.1098 -0.0575 -0.0370  0.0283  0.0631  0.7942  1.3624  1.6388
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.3862 -0.1330 -0.0629  0.0028  0.0280  0.0659  0.9937  1.1751  1.3218
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.2231 -0.1327 -0.1105  0.0003  0.0465  0.2002  0.6475  0.8618  1.1611
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4334 -0.1001 -0.0483 -0.0483  0.0389  0.0389  0.7743  1.6560  1.7096

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.144302;  11.000000 electrons
         Sum occ. bands:   -0.856636, incl. Bloechl correction: -0.006577

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.6635 -0.0593 -0.0537 -0.0537  0.0167  0.0167  1.6484  1.7917  1.8747
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.3287 -0.1260 -0.0735 -0.0257  0.0506  0.0918  0.6905  1.2552  1.3832
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4093 -0.1098 -0.0575 -0.0370  0.0283  0.0631  0.7942  1.3624  1.6388
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.3862 -0.1330 -0.0629  0.0028  0.0280  0.0659  0.9937  1.1751  1.3218
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.2231 -0.1327 -0.1105  0.0003  0.0465  0.2002  0.6475  0.8618  1.1611
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4334 -0.1001 -0.0483 -0.0483  0.0389  0.0389  0.7743  1.6560  1.7096

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.144302;  11.000000 electrons
         Sum occ. bands:   -0.856636, incl. Bloechl correction: -0.006577

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1    9.922201    3.108067    6.814134

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.511822   -0.383217    4.650000    4.647982    4.500000    4.647982
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.539365   -0.348398    4.340000    4.343510    4.250000    4.343510
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.506520   -0.057063    3.870000    3.850507    3.147584    3.850507
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.027214   -0.396285    4.110000    4.115514    4.102416    4.110000

 Harris energy:
 sumev=       -0.856636  val*vef=    -189.493805   sumtv=     188.637169
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -130.258096     utot=   -6534.971175    ehar=   -3304.835463

 srhov:     -6.344577   -168.052828   -174.397404 sumev=   -0.856636   sumtv=  173.540768

 Kohn-Sham energy:
 sumtv=      173.540768  sumtc=      3171.756639   ekin=     3345.297408
 rhoep=     -128.618450   utot=     -6521.191155   ehks=    -3304.512197
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 0.  RMS DQ=3.94e-2
 charges:       old           new         screened      rms diff       lin mix
 smooth       4.646654      4.185866      4.185866      0.038384      4.185866
 site    1    6.353346      6.814134      6.814134      0.019275      6.814134
 AMIX: nmix=0 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=3.94e-2
 unscreened rms difference:  smooth  0.045992   local  0.019275
   screened rms difference:  smooth  0.038384   local  0.019275   tot  0.039363

 iors  : write restart file (binary, mesh density) 

   it  1  of 12    ehf=      -0.400963   ehk=      -0.077697
h nk=8 bigbas=0 ehf=-.4009627 ehk=-.0776972

 --- BNDFP:  begin iteration 2 of 12 ---

 avg es pot at rmt= 0.634818  avg sphere pot= 0.654444  vconst=-0.634818

 smooth rhoves     12.288522   charge     4.185866
 smooth rhoeps =   -3.107690   rhomu =   -4.045533  avg vxc =   -0.858572 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000667

 Energy terms:             smooth           local           total
   rhoval*vef             -7.719041      -175.245661      -182.964702
   rhoval*ves            -48.983091      -107.880755      -156.863846
   psnuc*ves              73.560136    -12964.108731    -12890.548595
   utot                   12.288522     -6535.994743     -6523.706221
   rho*exc                -3.107690      -125.862959      -128.970649
   rho*vxc                -4.045533      -166.667780      -170.713314
   valence chg             4.185866         6.814134        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=0.144302

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7796 -0.6937 -0.6914 -0.6914 -0.6542 -0.6542  1.5097  1.6804  1.7088
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.7285 -0.6959 -0.6766 -0.6651 -0.6292 -0.3459  0.4540  1.0685  1.1964
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.7287 -0.6896 -0.6807 -0.6602 -0.6379 -0.4629  0.5636  1.1774  1.4849
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.7255 -0.7097 -0.6702 -0.6543 -0.6345 -0.4345  0.7854  0.9887  1.1275
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.7105 -0.7090 -0.6876 -0.6532 -0.6345 -0.1416  0.4097  0.6552  0.9628
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.7370 -0.6871 -0.6871 -0.6401 -0.6401 -0.4994  0.5398  1.5151  1.5151

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.196416;  11.000000 electrons
         Sum occ. bands:   -7.210164, incl. Bloechl correction: -0.013312

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.417460    1.997321    8.420138

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.426238   -0.663952    4.647982    4.595167    4.500000    4.595167
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.263999   -1.010586    4.343510    4.234553    4.250000    4.250000
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    9.219209   -0.671238    3.850507    3.905698    3.147584    3.905698
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.009001   -1.069889    4.110000    4.102060    4.102416    4.110000

 Harris energy:
 sumev=       -7.210164  val*vef=    -182.964702   sumtv=     175.754538
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -128.970649     utot=   -6523.706221    ehar=   -3305.165692

 srhov:     -4.518219   -222.289201   -226.807420 sumev=   -7.210164   sumtv=  219.597256

 Kohn-Sham energy:
 sumtv=      219.597256  sumtc=      3171.756639   ekin=     3391.353896
 rhoep=     -132.711133   utot=     -6562.374494   ehks=    -3303.731731
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 1.  RMS DQ=9.96e-2  last it=3.94e-2
 mixrho: (warning) scr. and lin-mixed densities had 43 and 43 negative points
 charges:       old           new         screened      rms diff       lin mix
 smooth       4.185866      2.579861      2.579861      0.030692      2.579861
 site    1    6.814134      8.420138      8.420138      0.064532      8.420138
 AMIX: nmix=1 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=9.96e-2
   tj: 0.81655
 unscreened rms difference:  smooth  0.022266   local  0.064532
   screened rms difference:  smooth  0.030692   local  0.064532   tot  0.099579

 iors  : write restart file (binary, mesh density) 

   it  2  of 12    ehf=      -0.731192   ehk=       0.702769
 From last iter    ehf=      -0.400963   ehk=      -0.077697
 diffe(q)= -0.330229 (0.099579)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf=-.7311922 ehk=.7027687

 --- BNDFP:  begin iteration 3 of 12 ---

 avg es pot at rmt= 0.606616  avg sphere pot= 0.655012  vconst=-0.606616

 smooth rhoves     10.987459   charge     3.891242
 smooth rhoeps =   -2.824599   rhomu =   -3.676015  avg vxc =   -0.843832 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000605

 Energy terms:             smooth           local           total
   rhoval*vef             -6.608382      -181.233792      -187.842174
   rhoval*ves            -47.698399      -113.462698      -161.161097
   psnuc*ves              69.673317    -12968.888690    -12899.215373
   utot                   10.987459     -6541.175694     -6530.188235
   rho*exc                -2.824599      -126.695667      -129.520266
   rho*vxc                -3.676015      -167.766479      -171.442495
   valence chg             3.891242         7.108758        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.196416

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7485 -0.4494 -0.4461 -0.4461 -0.3978 -0.3978  1.5525  1.7077  1.7610
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.5249 -0.4588 -0.4448 -0.4220 -0.3681 -0.2486  0.5221  1.1241  1.2528
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.5501 -0.4466 -0.4414 -0.4314 -0.3825 -0.3334  0.6309  1.2330  1.5323
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.5311 -0.4914 -0.4339 -0.4017 -0.3772 -0.3174  0.8475  1.0438  1.1860
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.4875 -0.4773 -0.4631 -0.4014 -0.3742 -0.0603  0.4774  0.7157  1.0218
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5688 -0.4411 -0.4411 -0.4038 -0.3807 -0.3807  0.6080  1.5742  1.5742

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.123421;  11.000000 electrons
         Sum occ. bands:   -4.734130, incl. Bloechl correction: -0.011752

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.290741    2.352004    7.938737

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.452179   -0.541931    4.595167    4.630389    4.500000    4.630389
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.309228   -0.845384    4.250000    4.254248    4.250000    4.254248
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    9.094100   -0.430435    3.905698    3.888285    3.147584    3.888285
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.013365   -0.830798    4.110000    4.106638    4.102416    4.110000

 Harris energy:
 sumev=       -4.734130  val*vef=    -187.842174   sumtv=     183.108044
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.520266     utot=   -6530.188235    ehar=   -3304.843816

 srhov:     -5.096096   -202.629365   -207.725461 sumev=   -4.734130   sumtv=  202.991332

 Kohn-Sham energy:
 sumtv=      202.991332  sumtc=      3171.756639   ekin=     3374.747971
 rhoep=     -131.310781   utot=     -6547.942637   ehks=    -3304.505447
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 2.  RMS DQ=4.68e-2  last it=9.96e-2
 mixrho: (warning) scr. and lin-mixed densities had 13 and 13 negative points
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.891242      3.061263      3.061263      0.015137      3.061263
 site    1    7.108758      7.938737      7.938737      0.031608      7.938737
 AMIX: nmix=2 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=4.68e-2
   tj:-0.97720  -0.09803
 add q= -0.000002 to preserve neutrality
 unscreened rms difference:  smooth  0.011190   local  0.031608
   screened rms difference:  smooth  0.015137   local  0.031608   tot  0.046841

 iors  : write restart file (binary, mesh density) 

   it  3  of 12    ehf=      -0.409316   ehk=      -0.070947
 From last iter    ehf=      -0.731192   ehk=       0.702769
 diffe(q)=  0.321876 (0.046841)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf=-.4093164 ehk=-.0709468

 --- BNDFP:  begin iteration 4 of 12 ---

 avg es pot at rmt= 0.549772  avg sphere pot= 0.664438  vconst=-0.549772

 smooth rhoves      8.887932   charge     3.421448
 smooth rhoeps =   -2.403836   rhomu =   -3.127177  avg vxc =   -0.816109 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000525

 Energy terms:             smooth           local           total
   rhoval*vef             -5.120210      -185.753861      -190.874071
   rhoval*ves            -44.968831      -118.542628      -163.511459
   psnuc*ves              62.744694    -12970.329623    -12907.584929
   utot                    8.887932     -6544.436125     -6535.548194
   rho*exc                -2.403836      -127.751594      -130.155430
   rho*vxc                -3.127177      -169.157240      -172.284417
   valence chg             3.421448         7.578552        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.123421

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.6951 -0.1460 -0.1409 -0.1409 -0.0747 -0.0747  1.6148  1.7598  1.8386
 Est Ef = -0.123 < evl(5)=-0.075 ... using qval=11.0, revise to -0.0747
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.3699 -0.2000 -0.1572 -0.1134 -0.0410  0.0045  0.6411  1.2147  1.3429
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4451 -0.1861 -0.1433 -0.1245 -0.0625 -0.0287  0.7463  1.3224  1.6042
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4223 -0.2131 -0.1456 -0.0864 -0.0613 -0.0256  0.9500  1.1345  1.2805
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.2759 -0.2081 -0.1899 -0.0880 -0.0460  0.1265  0.5974  0.8177  1.1186
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4688 -0.1741 -0.1353 -0.1353 -0.0534 -0.0534  0.7257  1.6306  1.6676

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.067281;  11.000000 electrons
         Sum occ. bands:   -1.718173, incl. Bloechl correction: -0.007514

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.004381    2.962314    7.042067

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.501522   -0.412805    4.630389    4.652114    4.500000    4.652114
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.480578   -0.458486    4.254248    4.322533    4.250000    4.322533
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.649829   -0.139618    3.888285    3.858541    3.147584    3.858541
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.023820   -0.489342    4.110000    4.113794    4.102416    4.110000

 Harris energy:
 sumev=       -1.718173  val*vef=    -190.874071   sumtv=     189.155898
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -130.155430     utot=   -6535.548194    ehar=   -3304.791086

 srhov:     -5.914183   -174.373709   -180.287893 sumev=   -1.718173   sumtv=  178.569720

 Kohn-Sham energy:
 sumtv=      178.569720  sumtc=      3171.756639   ekin=     3350.326359
 rhoep=     -129.105223   utot=     -6525.886929   ehks=    -3304.665793
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 3.  RMS DQ=2.67e-2  last it=4.68e-2
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.421448      3.957933      3.957933      0.008903      3.957933
 site    1    7.578552      7.042067      7.042067      0.019316      7.042067
 AMIX: nmix=3 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=2.67e-2
   tj: 0.74815  -0.22552  -0.00248
 unscreened rms difference:  smooth  0.006908   local  0.019316
   screened rms difference:  smooth  0.008903   local  0.019316   tot  0.026693

 iors  : write restart file (binary, mesh density) 

   it  4  of 12    ehf=      -0.356586   ehk=      -0.231293
 From last iter    ehf=      -0.409316   ehk=      -0.070947
 diffe(q)=  0.052730 (0.026693)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf=-.3565862 ehk=-.2312929

 --- BNDFP:  begin iteration 5 of 12 ---

 avg es pot at rmt= 0.573793  avg sphere pot= 0.659737  vconst=-0.573793

 smooth rhoves      9.689003   charge     3.597311
 smooth rhoeps =   -2.555829   rhomu =   -3.325364  avg vxc =   -0.827367 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000552

 Energy terms:             smooth           local           total
   rhoval*vef             -5.627008      -184.767460      -190.394468
   rhoval*ves            -46.126551      -117.134126      -163.260677
   psnuc*ves              65.504556    -12970.631936    -12905.127380
   utot                    9.689003     -6543.883031     -6534.194028
   rho*exc                -2.555829      -127.388303      -129.944132
   rho*vxc                -3.325364      -168.679143      -172.004508
   valence chg             3.597311         7.402689        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=0.067281

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7171 -0.2534 -0.2490 -0.2490 -0.1896 -0.1896  1.5908  1.7382  1.8096
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4133 -0.2890 -0.2607 -0.2229 -0.1569 -0.0987  0.5934  1.1793  1.3080
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4763 -0.2780 -0.2499 -0.2335 -0.1764 -0.1428  0.7003  1.2877  1.5767
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4544 -0.3108 -0.2477 -0.1981 -0.1730 -0.1383  0.9097  1.0989  1.2437
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3400 -0.3007 -0.2870 -0.1989 -0.1627  0.0473  0.5492  0.7772  1.0807
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4991 -0.2616 -0.2436 -0.2436 -0.1698 -0.1698  0.6787  1.6169  1.6321

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.015723;  11.000000 electrons
         Sum occ. bands:   -2.775913, incl. Bloechl correction: -0.009068

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.115862    2.733026    7.382836

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.478931   -0.451955    4.652114    4.649211    4.500000    4.649211
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.415874   -0.588781    4.322533    4.297408    4.250000    4.297408
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.834422   -0.243159    3.858541    3.867961    3.147584    3.867961
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.019551   -0.607979    4.110000    4.111266    4.102416    4.110000

 Harris energy:
 sumev=       -2.775913  val*vef=    -190.394468   sumtv=     187.618555
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.944132     utot=   -6534.194028    ehar=   -3304.762966

 srhov:     -5.636434   -184.305531   -189.941965 sumev=   -2.775913   sumtv=  187.166052

 Kohn-Sham energy:
 sumtv=      187.166052  sumtc=      3171.756639   ekin=     3358.922691
 rhoep=     -129.894678   utot=     -6533.790683   ehks=    -3304.762671
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=1.07e-3  last it=2.67e-2
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.597311      3.617164      3.617164      0.000264      3.617164
 site    1    7.402689      7.382836      7.382836      0.000766      7.382836
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: nmix=2 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=1.07e-3
   tj:-0.04838  -0.00376
 unscreened rms difference:  smooth  0.000323   local  0.000766
   screened rms difference:  smooth  0.000264   local  0.000766   tot  0.001067

 iors  : write restart file (binary, mesh density) 

   it  5  of 12    ehf=      -0.328466   ehk=      -0.328171
 From last iter    ehf=      -0.356586   ehk=      -0.231293
 diffe(q)=  0.028121 (0.001067)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf=-.3284655 ehk=-.3281706

 --- BNDFP:  begin iteration 6 of 12 ---

 avg es pot at rmt= 0.575386  avg sphere pot= 0.659742  vconst=-0.575386

 smooth rhoves      9.727581   charge     3.602770
 smooth rhoeps =   -2.559447   rhomu =   -3.330068  avg vxc =   -0.827915 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000550

 Energy terms:             smooth           local           total
   rhoval*vef             -5.633125      -184.764173      -190.397297
   rhoval*ves            -46.188018      -117.090376      -163.278394
   psnuc*ves              65.643179    -12970.627279    -12904.984100
   utot                    9.727581     -6543.858827     -6534.131247
   rho*exc                -2.559447      -127.371062      -129.930508
   rho*vxc                -3.330068      -168.656408      -171.986477
   valence chg             3.602770         7.397230        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.015723

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7189 -0.2614 -0.2570 -0.2570 -0.1981 -0.1981  1.5887  1.7363  1.8070
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4171 -0.2957 -0.2684 -0.2310 -0.1654 -0.1059  0.5899  1.1764  1.3051
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4790 -0.2849 -0.2578 -0.2416 -0.1848 -0.1511  0.6969  1.2849  1.5744
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4572 -0.3182 -0.2553 -0.2064 -0.1812 -0.1464  0.9066  1.0961  1.2407
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3454 -0.3077 -0.2942 -0.2071 -0.1712  0.0419  0.5456  0.7741  1.0777
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5018 -0.2681 -0.2516 -0.2516 -0.1783 -0.1783  0.6752  1.6156  1.6290

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.021334;  11.000000 electrons
         Sum occ. bands:   -2.855076, incl. Bloechl correction: -0.009182

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.124682    2.716650    7.408032

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.478178   -0.455704    4.649211    4.648722    4.500000    4.648722
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.408534   -0.603237    4.297408    4.294602    4.250000    4.294602
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.852125   -0.250444    3.867961    3.869093    3.147584    3.869093
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.019218   -0.617896    4.110000    4.111062    4.102416    4.110000

 Harris energy:
 sumev=       -2.855076  val*vef=    -190.397297   sumtv=     187.542221
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.930508     utot=   -6534.131247    ehar=   -3304.762895

 srhov:     -5.615174   -185.033148   -190.648322 sumev=   -2.855076   sumtv=  187.793245

 Kohn-Sham energy:
 sumtv=      187.793245  sumtc=      3171.756639   ekin=     3359.549885
 rhoep=     -129.953806   utot=     -6534.358914   ehks=    -3304.762835
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=6.05e-4  last it=1.07e-3
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.602770      3.591968      3.591968      0.000193      3.591968
 site    1    7.397230      7.408032      7.408032      0.000414      7.408032
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: nmix=2 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=6.05e-4
   tj: 0.25706   0.00637
 unscreened rms difference:  smooth  0.000147   local  0.000414
   screened rms difference:  smooth  0.000193   local  0.000414   tot  0.000605

 iors  : write restart file (binary, mesh density) 

   it  6  of 12    ehf=      -0.328395   ehk=      -0.328335
 From last iter    ehf=      -0.328466   ehk=      -0.328171
 diffe(q)=  0.000071 (0.000605)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf=-.3283946 ehk=-.3283355

 --- BNDFP:  begin iteration 7 of 12 ---

 avg es pot at rmt= 0.575138  avg sphere pot= 0.659741  vconst=-0.575138

 smooth rhoves      9.718905   charge     3.600774
 smooth rhoeps =   -2.557640   rhomu =   -3.327711  avg vxc =   -0.827804 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000550

 Energy terms:             smooth           local           total
   rhoval*vef             -5.626732      -184.792057      -190.418789
   rhoval*ves            -46.176500      -117.119526      -163.296026
   psnuc*ves              65.614311    -12970.648765    -12905.034453
   utot                    9.718905     -6543.884145     -6534.165240
   rho*exc                -2.557640      -127.376479      -129.934119
   rho*vxc                -3.327711      -168.663553      -171.991265
   valence chg             3.600774         7.399226        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.021334

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7187 -0.2598 -0.2554 -0.2554 -0.1963 -0.1963  1.5890  1.7365  1.8074
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4164 -0.2943 -0.2668 -0.2293 -0.1636 -0.1044  0.5906  1.1769  1.3056
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4785 -0.2835 -0.2562 -0.2399 -0.1830 -0.1494  0.6975  1.2853  1.5748
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4568 -0.3167 -0.2538 -0.2047 -0.1795 -0.1447  0.9071  1.0966  1.2413
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3443 -0.3063 -0.2928 -0.2054 -0.1694  0.0430  0.5463  0.7746  1.0782
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5014 -0.2668 -0.2500 -0.2500 -0.1766 -0.1766  0.6759  1.6158  1.6295

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.020218;  11.000000 electrons
         Sum occ. bands:   -2.838963, incl. Bloechl correction: -0.009157

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.123130    2.720139    7.402991

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.478553   -0.455001    4.648722    4.648849    4.500000    4.648849
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.409296   -0.601520    4.294602    4.294907    4.250000    4.294907
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.849940   -0.248855    3.869093    3.868967    3.147584    3.868967
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.019275   -0.616137    4.110000    4.111100    4.102416    4.110000

 Harris energy:
 sumev=       -2.838963  val*vef=    -190.418789   sumtv=     187.579827
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.934119     utot=   -6534.165240    ehar=   -3304.762893

 srhov:     -5.620436   -184.875371   -190.495807 sumev=   -2.838963   sumtv=  187.656845

 Kohn-Sham energy:
 sumtv=      187.656845  sumtc=      3171.756639   ekin=     3359.413484
 rhoep=     -129.941588   utot=     -6534.234783   ehks=    -3304.762887
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=1.94e-4  last it=6.05e-4
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.600774      3.597009      3.597009      0.000067      3.597009
 site    1    7.399226      7.402991      7.402991      0.000139      7.402991
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: nmix=2 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=1.94e-4
   tj:-0.18079   0.09548
 unscreened rms difference:  smooth  0.000052   local  0.000139
   screened rms difference:  smooth  0.000067   local  0.000139   tot  0.000194

 iors  : write restart file (binary, mesh density) 

   it  7  of 12    ehf=      -0.328393   ehk=      -0.328387
 From last iter    ehf=      -0.328395   ehk=      -0.328335
 diffe(q)=  0.000002 (0.000194)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf=-.3283931 ehk=-.3283868

 --- BNDFP:  begin iteration 8 of 12 ---

 avg es pot at rmt= 0.575016  avg sphere pot= 0.659761  vconst=-0.575016

 smooth rhoves      9.714763   charge     3.599845
 smooth rhoeps =   -2.556803   rhomu =   -3.326619  avg vxc =   -0.827751 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000550

 Energy terms:             smooth           local           total
   rhoval*vef             -5.623836      -184.796989      -190.420825
   rhoval*ves            -46.170948      -117.125972      -163.296920
   psnuc*ves              65.600473    -12970.646836    -12905.046363
   utot                    9.714763     -6543.886404     -6534.171641
   rho*exc                -2.556803      -127.378374      -129.935177
   rho*vxc                -3.326619      -168.666046      -171.992665
   valence chg             3.599845         7.400155        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.020218

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7186 -0.2592 -0.2549 -0.2549 -0.1957 -0.1957  1.5891  1.7367  1.8076
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4162 -0.2939 -0.2663 -0.2288 -0.1630 -0.1039  0.5908  1.1771  1.3058
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4784 -0.2830 -0.2557 -0.2394 -0.1825 -0.1488  0.6977  1.2855  1.5749
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4566 -0.3162 -0.2532 -0.2041 -0.1790 -0.1442  0.9073  1.0967  1.2414
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3440 -0.3058 -0.2923 -0.2048 -0.1688  0.0433  0.5465  0.7748  1.0784
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5012 -0.2663 -0.2494 -0.2494 -0.1760 -0.1760  0.6761  1.6158  1.6297

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.019839;  11.000000 electrons
         Sum occ. bands:   -2.833541, incl. Bloechl correction: -0.009149

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.122573    2.721273    7.401300

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.478643   -0.454778    4.648849    4.648876    4.500000    4.648876
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.409666   -0.600780    4.294907    4.295043    4.250000    4.295043
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.849025   -0.248337    3.868967    3.868908    3.147584    3.868908
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.019296   -0.615507    4.110000    4.111113    4.102416    4.110000

 Harris energy:
 sumev=       -2.833541  val*vef=    -190.420825   sumtv=     187.587284
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.935177     utot=   -6534.171641    ehar=   -3304.762895

 srhov:     -5.621858   -184.824601   -190.446458 sumev=   -2.833541   sumtv=  187.612917

 Kohn-Sham energy:
 sumtv=      187.612917  sumtc=      3171.756639   ekin=     3359.369557
 rhoep=     -129.937550   utot=     -6534.194900   ehks=    -3304.762894
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=6.22e-5  last it=1.94e-4
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.599845      3.598700      3.598700      0.000022      3.598700
 site    1    7.400155      7.401300      7.401300      0.000043      7.401300
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=6.22e-5
   tj:-0.46822
 unscreened rms difference:  smooth  0.000017   local  0.000043
   screened rms difference:  smooth  0.000022   local  0.000043   tot  0.000062

 iors  : write restart file (binary, mesh density) 

   it  8  of 12    ehf=      -0.328395   ehk=      -0.328394
 From last iter    ehf=      -0.328393   ehk=      -0.328387
 diffe(q)= -0.000002 (0.000062)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 ehf=-.3283948 ehk=-.3283941

 --- BNDFP:  begin iteration 9 of 12 ---

 avg es pot at rmt= 0.574973  avg sphere pot= 0.659766  vconst=-0.574973

 smooth rhoves      9.713256   charge     3.599491
 smooth rhoeps =   -2.556475   rhomu =   -3.326192  avg vxc =   -0.827733 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000550

 Energy terms:             smooth           local           total
   rhoval*vef             -5.622675      -184.801107      -190.423782
   rhoval*ves            -46.168975      -117.130357      -163.299332
   psnuc*ves              65.595487    -12970.648910    -12905.053423
   utot                    9.713256     -6543.889633     -6534.176377
   rho*exc                -2.556475      -127.379212      -129.935687
   rho*vxc                -3.326192      -168.667151      -171.993343
   valence chg             3.599491         7.400509        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.019839

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7185 -0.2590 -0.2546 -0.2546 -0.1955 -0.1955  1.5891  1.7367  1.8076
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4161 -0.2937 -0.2660 -0.2286 -0.1628 -0.1037  0.5909  1.1772  1.3059
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4783 -0.2828 -0.2554 -0.2391 -0.1822 -0.1486  0.6978  1.2856  1.5750
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4565 -0.3160 -0.2530 -0.2039 -0.1787 -0.1439  0.9074  1.0968  1.2415
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3438 -0.3056 -0.2921 -0.2046 -0.1686  0.0435  0.5466  0.7749  1.0785
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5011 -0.2662 -0.2492 -0.2492 -0.1757 -0.1757  0.6762  1.6158  1.6297

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.019679;  11.000000 electrons
         Sum occ. bands:   -2.831253, incl. Bloechl correction: -0.009146

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.122332    2.721782    7.400550

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.478685   -0.454687    4.648876    4.648889    4.500000    4.648889
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.409824   -0.600476    4.295043    4.295100    4.250000    4.295100
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.848627   -0.248118    3.868908    3.868883    3.147584    3.868883
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.019305   -0.615244    4.110000    4.111119    4.102416    4.110000

 Harris energy:
 sumev=       -2.831253  val*vef=    -190.423782   sumtv=     187.592529
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.935687     utot=   -6534.176377    ehar=   -3304.762896

 srhov:     -5.622574   -184.801411   -190.423985 sumev=   -2.831253   sumtv=  187.592732

 Kohn-Sham energy:
 sumtv=      187.592732  sumtc=      3171.756639   ekin=     3359.349371
 rhoep=     -129.935725   utot=     -6534.176543   ehks=    -3304.762896
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=1.71e-6  last it=6.22e-5
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.599491      3.599450      3.599450      0.000006      3.599450
 site    1    7.400509      7.400550      7.400550      0.000002      7.400550
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=1.71e-6
   tj:-0.02042
 unscreened rms difference:  smooth  0.000006   local  0.000002
   screened rms difference:  smooth  0.000006   local  0.000002   tot  0.000002

 iors  : write restart file (binary, mesh density) 

   it  9  of 12    ehf=      -0.328396   ehk=      -0.328396
 From last iter    ehf=      -0.328395   ehk=      -0.328394
 diffe(q)= -0.000001 (0.000002)    tol= 0.000010 (0.000010)   more=F
c nk=8 bigbas=0 ehf=-.3283963 ehk=-.3283962
 Exit 0 LMF 
 CPU time:   20.268s     Sat Jun 20 14:50:07 2009   on waldo.eas.asu.edu
 wkinfo:  used   410 K  workspace of 80000 K   in  27 K calls
rdcmd:  rm mixm.cu
rdcmd:  lmf  --no-iactiv cu -vnk=8 -vbigbas=t -vpwmode=0 -voveps=0d-7
 -----------------------  START LMF (80000K)  -----------------------

 LMF:      nbas = 1  nspec = 1  vn 7.00(LMF 7.0)  verb 31,30
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
 spec      l    rsm    eh     gmax    last term    cutoff
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

 avg es pot at rmt= 0.574971  avg sphere pot= 0.659767  vconst=-0.574971

 smooth rhoves      9.713156   charge     3.599465
 smooth rhoeps =   -2.556450   rhomu =   -3.326159  avg vxc =   -0.827732 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000550

 Energy terms:             smooth           local           total
   rhoval*vef             -5.622581      -184.800963      -190.423544
   rhoval*ves            -46.168853      -117.130239      -163.299093
   psnuc*ves              65.595165    -12970.648417    -12905.053251
   utot                    9.713156     -6543.889328     -6534.176172
   rho*exc                -2.556450      -127.379238      -129.935688
   rho*vxc                -3.326159      -168.667185      -171.993344
   valence chg             3.599465         7.400535        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7182 -0.2591 -0.2548 -0.2548 -0.1954 -0.1954  1.5858  1.7366  1.8040
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4159 -0.2937 -0.2661 -0.2286 -0.1628 -0.1035  0.5903  1.1765  1.3085
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4780 -0.2829 -0.2555 -0.2392 -0.1823 -0.1485  0.6971  1.2855  1.5941
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4563 -0.3161 -0.2531 -0.2040 -0.1788 -0.1439  0.9063  1.0954  1.2459
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3438 -0.3057 -0.2922 -0.2046 -0.1687  0.0437  0.5461  0.7733  1.0770
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5007 -0.2660 -0.2493 -0.2493 -0.1758 -0.1758  0.6757  1.6267  1.6267

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.019383;  11.000000 electrons
         Sum occ. bands:   -2.831549, incl. Bloechl correction: -0.009152

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7182 -0.2591 -0.2548 -0.2548 -0.1954 -0.1954  1.5858  1.7366  1.8040
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4159 -0.2937 -0.2661 -0.2286 -0.1628 -0.1035  0.5903  1.1765  1.3085
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4780 -0.2829 -0.2555 -0.2392 -0.1823 -0.1485  0.6971  1.2855  1.5941
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4563 -0.3161 -0.2531 -0.2040 -0.1788 -0.1439  0.9063  1.0954  1.2459
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3438 -0.3057 -0.2922 -0.2046 -0.1687  0.0437  0.5461  0.7733  1.0770
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5007 -0.2660 -0.2493 -0.2493 -0.1758 -0.1758  0.6757  1.6267  1.6267

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.019383;  11.000000 electrons
         Sum occ. bands:   -2.831549, incl. Bloechl correction: -0.009152

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.128151    2.833599    7.294552

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.480402   -0.533976    4.648889    4.612321    4.500000    4.612321
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.406751   -0.470296    4.295100    4.328198    4.250000    4.328198
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.972454   -0.256709    3.868883    3.858404    3.147584    3.858404
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.021166   -0.809396    4.110000    4.106671    4.102416    4.110000
 4     1    0.004388   -0.785665    5.100000    5.079456    5.077979    5.100000

 Harris energy:
 sumev=       -2.831549  val*vef=    -190.423544   sumtv=     187.591995
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.935688     utot=   -6534.176172    ehar=   -3304.763226

 srhov:     -5.955420   -184.376042   -190.331462 sumev=   -2.831549   sumtv=  187.499912

 Kohn-Sham energy:
 sumtv=      187.499912  sumtc=      3171.756639   ekin=     3359.256552
 rhoep=     -129.930570   utot=     -6534.089135   ehks=    -3304.763154
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 0.  RMS DQ=3.17e-3
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.599465      3.705448      3.705448      0.002541      3.705448
 site    1    7.400535      7.294552      7.294552      0.003499      7.294552
 AMIX: nmix=0 mmix=8  nelts=2178  beta=1  tm=5  rmsdel=3.17e-3
 unscreened rms difference:  smooth  0.002542   local  0.003499
   screened rms difference:  smooth  0.002541   local  0.003499   tot  0.003171

 iors  : write restart file (binary, mesh density) 

   it  1  of 12    ehf=      -0.328726   ehk=      -0.328654
i nk=8 bigbas=1 pwmode=0 oveps=0 ehf=-.3287262 ehk=-.3286535

 --- BNDFP:  begin iteration 2 of 12 ---

 avg es pot at rmt= 0.581265  avg sphere pot= 0.652951  vconst=-0.581265

 smooth rhoves     10.109215   charge     3.705448
 smooth rhoeps =   -2.658226   rhomu =   -3.459002  avg vxc =   -0.832938 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000130

 Energy terms:             smooth           local           total
   rhoval*vef             -6.011395      -184.328635      -190.340030
   rhoval*ves            -46.614983      -116.605499      -163.220482
   psnuc*ves              66.833414    -12971.792603    -12904.959190
   utot                   10.109215     -6544.199051     -6534.089836
   rho*exc                -2.658226      -127.272729      -129.930955
   rho*vxc                -3.459002      -168.528048      -171.987051
   valence chg             3.705448         7.294552        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7176 -0.2606 -0.2563 -0.2563 -0.1954 -0.1954  1.5862  1.7373  1.8047
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4159 -0.2945 -0.2668 -0.2293 -0.1638 -0.1043  0.5903  1.1769  1.3090
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4777 -0.2840 -0.2561 -0.2400 -0.1831 -0.1493  0.6972  1.2858  1.5947
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4558 -0.3174 -0.2536 -0.2050 -0.1794 -0.1447  0.9065  1.0958  1.2464
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3441 -0.3063 -0.2931 -0.2047 -0.1703  0.0433  0.5461  0.7736  1.0775
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5004 -0.2670 -0.2501 -0.2501 -0.1765 -0.1765  0.6757  1.6272  1.6272

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.019745;  11.000000 electrons
         Sum occ. bands:   -2.839039, incl. Bloechl correction: -0.009179

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7176 -0.2606 -0.2563 -0.2563 -0.1954 -0.1954  1.5862  1.7373  1.8047
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4159 -0.2945 -0.2668 -0.2293 -0.1638 -0.1043  0.5903  1.1769  1.3090
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4777 -0.2840 -0.2561 -0.2400 -0.1831 -0.1493  0.6972  1.2858  1.5947
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4558 -0.3174 -0.2536 -0.2050 -0.1794 -0.1447  0.9065  1.0958  1.2464
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3441 -0.3063 -0.2931 -0.2047 -0.1703  0.0433  0.5461  0.7736  1.0775
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5004 -0.2670 -0.2501 -0.2501 -0.1765 -0.1765  0.6757  1.6272  1.6272

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.019745;  11.000000 electrons
         Sum occ. bands:   -2.839039, incl. Bloechl correction: -0.009179

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.129986    2.828121    7.301865

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.479656   -0.534334    4.612321    4.612124    4.500000    4.612124
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.405904   -0.471411    4.328198    4.327920    4.250000    4.327920
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.966180   -0.257146    3.858404    3.858566    3.147584    3.858566
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.021089   -0.808370    4.110000    4.106698    4.102416    4.110000
 4     1    0.004385   -0.782341    5.100000    5.079492    5.077979    5.100000

 Harris energy:
 sumev=       -2.839039  val*vef=    -190.340030   sumtv=     187.500991
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.930955     utot=   -6534.089836    ehar=   -3304.763161

 srhov:     -5.998908   -184.526812   -190.525719 sumev=   -2.839039   sumtv=  187.686680

 Kohn-Sham energy:
 sumtv=      187.686680  sumtc=      3171.756639   ekin=     3359.443319
 rhoep=     -129.946555   utot=     -6534.259899   ehks=    -3304.763135
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 1.  RMS DQ=3.97e-4  last it=3.17e-3
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.705448      3.698134      3.698134      0.000183      3.698134
 site    1    7.294552      7.301865      7.301865      0.000265      7.301865
 AMIX: nmix=1 mmix=8  nelts=2178  beta=1  tm=5  rmsdel=3.97e-4
   tj: 0.06129
 unscreened rms difference:  smooth  0.000175   local  0.000265
   screened rms difference:  smooth  0.000183   local  0.000265   tot  0.000397

 iors  : write restart file (binary, mesh density) 

   it  2  of 12    ehf=      -0.328661   ehk=      -0.328635
 From last iter    ehf=      -0.328726   ehk=      -0.328654
 diffe(q)=  0.000066 (0.000397)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=1 pwmode=0 oveps=0 ehf=-.3286606 ehk=-.3286348

 --- BNDFP:  begin iteration 3 of 12 ---

 avg es pot at rmt= 0.580471  avg sphere pot= 0.653041  vconst=-0.580471

 smooth rhoves     10.079267   charge     3.698583
 smooth rhoeps =   -2.651841   rhomu =   -3.450671  avg vxc =   -0.832576 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000134

 Energy terms:             smooth           local           total
   rhoval*vef             -5.988555      -184.455814      -190.444369
   rhoval*ves            -46.578001      -116.733976      -163.311977
   psnuc*ves              66.736534    -12971.887505    -12905.150972
   utot                   10.079267     -6544.310741     -6534.231474
   rho*exc                -2.651841      -127.291221      -129.943061
   rho*vxc                -3.450671      -168.552441      -172.003112
   valence chg             3.698583         7.301417        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7168 -0.2556 -0.2513 -0.2513 -0.1901 -0.1901  1.5871  1.7381  1.8059
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4138 -0.2904 -0.2620 -0.2242 -0.1583 -0.0997  0.5923  1.1783  1.3105
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4763 -0.2797 -0.2511 -0.2349 -0.1778 -0.1441  0.6991  1.2873  1.5959
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4544 -0.3129 -0.2489 -0.1998 -0.1742 -0.1396  0.9081  1.0972  1.2480
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3409 -0.3020 -0.2886 -0.1995 -0.1648  0.0466  0.5481  0.7753  1.0790
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4991 -0.2631 -0.2450 -0.2450 -0.1711 -0.1711  0.6777  1.6287  1.6287

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.016318;  11.000000 electrons
         Sum occ. bands:   -2.789683, incl. Bloechl correction: -0.009102

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7168 -0.2556 -0.2513 -0.2513 -0.1901 -0.1901  1.5871  1.7381  1.8059
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4138 -0.2904 -0.2620 -0.2242 -0.1583 -0.0997  0.5923  1.1783  1.3105
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4763 -0.2797 -0.2511 -0.2349 -0.1778 -0.1441  0.6991  1.2873  1.5959
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4544 -0.3129 -0.2489 -0.1998 -0.1742 -0.1396  0.9081  1.0972  1.2480
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3409 -0.3020 -0.2886 -0.1995 -0.1648  0.0466  0.5481  0.7753  1.0790
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.4991 -0.2631 -0.2450 -0.2450 -0.1711 -0.1711  0.6777  1.6287  1.6287

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.016318;  11.000000 electrons
         Sum occ. bands:   -2.789683, incl. Bloechl correction: -0.009102

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.124984    2.838305    7.286679

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.480590   -0.531703    4.612124    4.612833    4.500000    4.612833
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.408834   -0.467563    4.327920    4.328677    4.250000    4.328677
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.956815   -0.252223    3.858566    3.858276    3.147584    3.858276
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.021299   -0.802890    4.110000    4.106807    4.102416    4.110000
 4     1    0.004418   -0.776479    5.100000    5.079545    5.077979    5.100000

 Harris energy:
 sumev=       -2.789683  val*vef=    -190.444369   sumtv=     187.654686
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.943061     utot=   -6534.231474    ehar=   -3304.763210

 srhov:     -6.012286   -184.016357   -190.028643 sumev=   -2.789683   sumtv=  187.238960

 Kohn-Sham energy:
 sumtv=      187.238960  sumtc=      3171.756639   ekin=     3358.995599
 rhoep=     -129.906696   utot=     -6533.851972   ehks=    -3304.763069
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 2.  RMS DQ=9.10e-4  last it=3.97e-4
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.698583      3.713320      3.713320      0.000264      3.713320
 site    1    7.301417      7.286679      7.286679      0.000610      7.286679
 AMIX: nmix=2 mmix=8  nelts=2178  beta=1  tm=5  rmsdel=9.1e-4
   tj: 0.69803   0.00606
 unscreened rms difference:  smooth  0.000192   local  0.000610
   screened rms difference:  smooth  0.000264   local  0.000610   tot  0.000910

 iors  : write restart file (binary, mesh density) 

   it  3  of 12    ehf=      -0.328710   ehk=      -0.328569
 From last iter    ehf=      -0.328661   ehk=      -0.328635
 diffe(q)= -0.000049 (0.000910)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=1 pwmode=0 oveps=0 ehf=-.3287098 ehk=-.3285688

 --- BNDFP:  begin iteration 4 of 12 ---

 avg es pot at rmt= 0.580896  avg sphere pot= 0.653049  vconst=-0.580896

 smooth rhoves     10.096579   charge     3.702673
 smooth rhoeps =   -2.655576   rhomu =   -3.455543  avg vxc =   -0.832808 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000135

 Energy terms:             smooth           local           total
   rhoval*vef             -6.002504      -184.374644      -190.377149
   rhoval*ves            -46.599448      -116.654168      -163.253616
   psnuc*ves              66.792606    -12971.815857    -12905.023252
   utot                   10.096579     -6544.235013     -6534.138434
   rho*exc                -2.655576      -127.279160      -129.934736
   rho*vxc                -3.455543      -168.536525      -171.992068
   valence chg             3.702673         7.297327        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7173 -0.2591 -0.2547 -0.2547 -0.1938 -0.1938  1.5865  1.7376  1.8051
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4152 -0.2932 -0.2654 -0.2277 -0.1621 -0.1029  0.5910  1.1774  1.3095
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4772 -0.2827 -0.2546 -0.2384 -0.1815 -0.1477  0.6978  1.2863  1.5951
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4554 -0.3160 -0.2521 -0.2034 -0.1778 -0.1431  0.9070  1.0963  1.2469
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3431 -0.3050 -0.2917 -0.2031 -0.1686  0.0443  0.5468  0.7742  1.0780
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5000 -0.2658 -0.2485 -0.2485 -0.1749 -0.1749  0.6763  1.6277  1.6277

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018673;  11.000000 electrons
         Sum occ. bands:   -2.823780, incl. Bloechl correction: -0.009156

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7173 -0.2591 -0.2547 -0.2547 -0.1938 -0.1938  1.5865  1.7376  1.8051
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4152 -0.2932 -0.2654 -0.2277 -0.1621 -0.1029  0.5910  1.1774  1.3095
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4772 -0.2827 -0.2546 -0.2384 -0.1815 -0.1477  0.6978  1.2863  1.5951
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4554 -0.3160 -0.2521 -0.2034 -0.1778 -0.1431  0.9070  1.0963  1.2469
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3431 -0.3050 -0.2917 -0.2031 -0.1686  0.0443  0.5468  0.7742  1.0780
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5000 -0.2658 -0.2485 -0.2485 -0.1749 -0.1749  0.6763  1.6277  1.6277

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018673;  11.000000 electrons
         Sum occ. bands:   -2.823780, incl. Bloechl correction: -0.009156

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.128433    2.831396    7.297037

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.479967   -0.533647    4.612833    4.612280    4.500000    4.612280
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.406873   -0.470211    4.328677    4.328153    4.250000    4.328153
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.963165   -0.255643    3.858276    3.858471    3.147584    3.858471
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.021153   -0.806707    4.110000    4.106730    4.102416    4.110000
 4     1    0.004395   -0.780599    5.100000    5.079507    5.077979    5.100000

 Harris energy:
 sumev=       -2.823780  val*vef=    -190.377149   sumtv=     187.553369
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.934736     utot=   -6534.138434    ehar=   -3304.763162

 srhov:     -6.002871   -184.367743   -190.370614 sumev=   -2.823780   sumtv=  187.546834

 Kohn-Sham energy:
 sumtv=      187.546834  sumtc=      3171.756639   ekin=     3359.303474
 rhoep=     -129.934115   utot=     -6534.132520   ehks=    -3304.763162
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 3.  RMS DQ=1.44e-5  last it=9.10e-4
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.702673      3.702963      3.702963      0.000014      3.702963
 site    1    7.297327      7.297037      7.297037      0.000010      7.297037
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=2178  beta=1  tm=5  rmsdel=1.44e-5
   tj:-0.01587
 unscreened rms difference:  smooth  0.000018   local  0.000010
   screened rms difference:  smooth  0.000014   local  0.000010   tot  0.000014

 iors  : write restart file (binary, mesh density) 

   it  4  of 12    ehf=      -0.328662   ehk=      -0.328662
 From last iter    ehf=      -0.328710   ehk=      -0.328569
 diffe(q)=  0.000048 (0.000014)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=1 pwmode=0 oveps=0 ehf=-.3286618 ehk=-.3286616

 --- BNDFP:  begin iteration 5 of 12 ---

 avg es pot at rmt= 0.580900  avg sphere pot= 0.653045  vconst=-0.580900

 smooth rhoves     10.097152   charge     3.702799
 smooth rhoeps =   -2.655668   rhomu =   -3.455663  avg vxc =   -0.832821 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000135

 Energy terms:             smooth           local           total
   rhoval*vef             -6.002874      -184.373519      -190.376393
   rhoval*ves            -46.600168      -116.652852      -163.253020
   psnuc*ves              66.794472    -12971.815755    -12905.021283
   utot                   10.097152     -6544.234304     -6534.137152
   rho*exc                -2.655668      -127.278923      -129.934591
   rho*vxc                -3.455663      -168.536212      -171.991876
   valence chg             3.702799         7.297201        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7173 -0.2592 -0.2548 -0.2548 -0.1938 -0.1938  1.5865  1.7376  1.8051
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4153 -0.2933 -0.2654 -0.2278 -0.1622 -0.1029  0.5910  1.1773  1.3095
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4773 -0.2827 -0.2546 -0.2385 -0.1816 -0.1478  0.6978  1.2863  1.5951
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4554 -0.3161 -0.2522 -0.2035 -0.1779 -0.1432  0.9070  1.0962  1.2469
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3431 -0.3050 -0.2918 -0.2032 -0.1687  0.0443  0.5467  0.7742  1.0780
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5000 -0.2659 -0.2486 -0.2486 -0.1750 -0.1750  0.6763  1.6277  1.6277

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018717;  11.000000 electrons
         Sum occ. bands:   -2.824451, incl. Bloechl correction: -0.009157

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250
 -0.7173 -0.2592 -0.2548 -0.2548 -0.1938 -0.1938  1.5865  1.7376  1.8051
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250
 -0.4153 -0.2933 -0.2654 -0.2278 -0.1622 -0.1029  0.5910  1.1773  1.3095
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750
 -0.4773 -0.2827 -0.2546 -0.2385 -0.1816 -0.1478  0.6978  1.2863  1.5951
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250
 -0.4554 -0.3161 -0.2522 -0.2035 -0.1779 -0.1432  0.9070  1.0962  1.2469
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750
 -0.3431 -0.3050 -0.2918 -0.2032 -0.1687  0.0443  0.5467  0.7742  1.0780
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250
 -0.5000 -0.2659 -0.2486 -0.2486 -0.1750 -0.1750  0.6763  1.6277  1.6277

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.018717;  11.000000 electrons
         Sum occ. bands:   -2.824451, incl. Bloechl correction: -0.009157

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.128515    2.831284    7.297231

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.479951   -0.533744    4.612280    4.612242    4.500000    4.612242
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.406819   -0.470275    4.328153    4.328140    4.250000    4.328140
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.963407   -0.255707    3.858471    3.858480    3.147584    3.858480
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.021150   -0.806786    4.110000    4.106729    4.102416    4.110000
 4     1    0.004394   -0.780680    5.100000    5.079506    5.077979    5.100000

 Harris energy:
 sumev=       -2.824451  val*vef=    -190.376393   sumtv=     187.551941
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.934591     utot=   -6534.137152    ehar=   -3304.763162

 srhov:     -6.002810   -184.374450   -190.377260 sumev=   -2.824451   sumtv=  187.552808

 Kohn-Sham energy:
 sumtv=      187.552808  sumtc=      3171.756639   ekin=     3359.309448
 rhoep=     -129.934667   utot=     -6534.137942   ehks=    -3304.763161
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=2.08e-6  last it=1.44e-5
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.702799      3.702768      3.702768      0.000008      3.702768
 site    1    7.297201      7.297231      7.297231      0.000001      7.297231
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=2178  beta=1  tm=5  rmsdel=2.08e-6
   tj: 0.12353
 unscreened rms difference:  smooth  0.000008   local  0.000001
   screened rms difference:  smooth  0.000008   local  0.000001   tot  0.000002

 iors  : write restart file (binary, mesh density) 

   it  5  of 12    ehf=      -0.328662   ehk=      -0.328661
 From last iter    ehf=      -0.328662   ehk=      -0.328662
 diffe(q)=  0.000000 (0.000002)    tol= 0.000010 (0.000010)   more=F
c nk=8 bigbas=1 pwmode=0 oveps=0 ehf=-.3286616 ehk=-.3286615
 Exit 0 LMF 
 CPU time:   50.178s     Sat Jun 20 14:50:59 2009   on waldo.eas.asu.edu
 wkinfo:  used   626 K  workspace of 80000 K   in  21 K calls
rdcmd:  lmf  --no-iactiv cu -vnk=8 -vbigbas=t -vpwmode=0 -voveps=0d-7 --band:fn=syml
 -----------------------  START LMF (80000K)  -----------------------

 LMF:      nbas = 1  nspec = 1  vn 7.00(LMF 7.0)  verb 31,30
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
 spec      l    rsm    eh     gmax    last term    cutoff
  A        0*   2.50  -0.01   2.974    2.30E-05      27 
  A        1    2.50  -0.01   3.093    1.29E-06      51 
  A        2    1.00  -0.01   8.508    1.16E-06     869 
  A        0    1.30  -1.00   5.718    2.28E-06     259 
  A        2    1.00  -1.00   8.508    1.16E-06     869 
  A        3    1.30  -0.01   6.806    2.09E-06     411 

 iors  : read restart file (binary, mesh density) 
         use from  restart file: ef window, positions, pnu 
         ignore in restart file: *

 --- BNDFP:  begin iteration 1 of 12 ---

 avg es pot at rmt= 0.580896  avg sphere pot= 0.653044  vconst=-0.580896

 smooth rhoves     10.097120   charge     3.702793
 smooth rhoeps =   -2.655657   rhomu =   -3.455649  avg vxc =   -0.832822 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000135

 Energy terms:             smooth           local           total
   rhoval*vef             -6.002849      -184.373634      -190.376484
   rhoval*ves            -46.600126      -116.652965      -163.253091
   psnuc*ves              66.794366    -12971.815898    -12905.021531
   utot                   10.097120     -6544.234431     -6534.137311
   rho*exc                -2.655657      -127.278953      -129.934610
   rho*vxc                -3.455649      -168.536252      -171.991901
   valence chg             3.702793         7.297207        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read efermi from weights file : ef = -0.018717
 suqlst:  generate bands, mode 1

 suqlst:  nq= 41   q1= 0.5000 0.5000 0.5000   q2= 0.0000 0.0000 0.0000
 bndfp:  kpt 1 of 41, k=  -0.50000  -0.50000  -0.50000
 -0.4126 -0.2588 -0.2588 -0.1500 -0.1500 -0.1050  0.2487  1.5877  1.5973
 bndfp:  kpt 11 of 41, k=  0.37500  0.37500  0.37500
 -0.4449 -0.2536 -0.2536 -0.2248 -0.1624 -0.1624  0.4812  1.6094  1.6094
 bndfp:  kpt 21 of 41, k=  0.25000  0.25000  0.25000
 -0.5703 -0.2760 -0.2447 -0.2447 -0.1870 -0.1870  0.8928  1.6566  1.6566
 bndfp:  kpt 31 of 41, k=  0.12500  0.12500  0.12500
 -0.6860 -0.2652 -0.2497 -0.2497 -0.1953 -0.1953  1.3655  1.7127  1.7482
 bndfp:  kpt 41 of 41, k=  0.00000  0.00000  0.00000
 -0.7279 -0.2569 -0.2569 -0.2569 -0.1930 -0.1930  1.6567  1.8387  1.8387

 Read efermi from weights file : ef = -0.018717

 suqlst:  nq= 41   q1= 0.0000 0.0000 0.0000   q2= 1.0000 0.0000 0.0000
 bndfp:  kpt 1 of 41, k=  0.00000  0.00000  0.00000
 -0.7279 -0.2569 -0.2569 -0.2569 -0.1930 -0.1930  1.6567  1.8387  1.8387
 bndfp:  kpt 11 of 41, k=  0.25000  0.00000  0.00000
 -0.6722 -0.2732 -0.2433 -0.2433 -0.2075 -0.1870  1.5242  1.5425  1.5425
 bndfp:  kpt 21 of 41, k=  0.50000  0.00000  0.00000
 -0.5227 -0.3114 -0.2283 -0.2064 -0.2064 -0.1723  1.2328  1.2328  1.3539
 bndfp:  kpt 31 of 41, k=  0.75000  0.00000  0.00000
 -0.4048 -0.3480 -0.1614 -0.1614 -0.1572 -0.1230  0.8184  1.0147  1.0147
 bndfp:  kpt 41 of 41, k=  -1.00000  0.00000  0.00000
 -0.3961 -0.3627 -0.1509 -0.1394 -0.1394  0.0772  0.4995  0.9336  0.9336

 Read efermi from weights file : ef = -0.018717

 suqlst:  nq= 21   q1= 1.0000 0.0000 0.0000   q2= 1.0000 0.5000 0.0000
 bndfp:  kpt 1 of 21, k=  -1.00000  0.00000  0.00000
 -0.3961 -0.3627 -0.1509 -0.1394 -0.1394  0.0772  0.4995  0.9336  0.9336
 bndfp:  kpt 11 of 21, k=  -1.00000  0.25000  0.00000
 -0.3718 -0.3463 -0.2102 -0.1698 -0.1394  0.1977  0.5208  0.7107  0.7756
 bndfp:  kpt 21 of 21, k=  -1.00000  0.50000  0.00000
 -0.3401 -0.2953 -0.2953 -0.1964 -0.1394  0.4291  0.4291  0.5910  0.6156

 Read efermi from weights file : ef = -0.018717

 suqlst:  nq= 41   q1= 1.0000 0.5000 0.0000   q2= 0.0000 0.0000 0.0000
 bndfp:  kpt 1 of 41, k=  -1.00000  0.50000  0.00000
 -0.3401 -0.2953 -0.2953 -0.1964 -0.1394  0.4291  0.4291  0.5910  0.6156
 bndfp:  kpt 11 of 41, k=  0.75000  0.37500  0.00000
 -0.3461 -0.3172 -0.2822 -0.1886 -0.1617  0.0625  0.6342  0.7293  0.9129
 bndfp:  kpt 21 of 41, k=  0.50000  0.25000  0.00000
 -0.4790 -0.3054 -0.2569 -0.2057 -0.2018 -0.1412  0.9570  1.0059  1.4114
 bndfp:  kpt 31 of 41, k=  0.25000  0.12500  0.00000
 -0.6587 -0.2738 -0.2436 -0.2428 -0.2061 -0.1870  1.3635  1.4235  1.7178
 bndfp:  kpt 41 of 41, k=  0.00000  0.00000  0.00000
 -0.7279 -0.2569 -0.2569 -0.2569 -0.1930 -0.1930  1.6567  1.8387  1.8387

 Read efermi from weights file : ef = -0.018717
 Exit 0 bndfp 
 CPU time:   10.624s     Sat Jun 20 14:51:11 2009   on waldo.eas.asu.edu
 wkinfo:  used   555 K  workspace of 80000 K   in   4 K calls
