rdcmd:  lmfa --no-iactiv cu -vnk=8 -vbigbas=f
 -----------------------  START LMFA -----------------------

 rdctrl: reset global max nl from 5 to 4

 LMFA:     nbas = 1  nspec = 1  vn 7.11.c  verb 31,30,31
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
 CPU time:    0.110s   Wall clock    0.159s   12:04:01 02.10.2014        on phpdl1.ph.kcl.ac.uk
rdcmd:  lmf  --no-iactiv cu -vnk=8 -vbigbas=f -vpwmode=11
 -----------------------  START LMF -----------------------

 rdctrl: reset global max nl from 5 to 4

 LMF:      nbas = 1  nspec = 1  vn 7.11.c  verb 31,30,31
 special:  APW basis(q)
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
 spec      l    rsm    eh     gmax    last term   cutoff
  A        0*   2.50  -0.01   2.974    2.30E-05      27 
  A        1    2.50  -0.01   3.093    1.29E-06      51 
  A        2    1.00  -0.01   8.508    1.16E-06     869 
 

 suham:  q-dependent PW basis with  Emin = 1 < E < 3.
         Est. min,max PW dimension = 3,10.  Use npwpad = 3 => ndham = 31

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

 Incompatible or missing qp weights file ...

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250   ndimh = 26
 -0.6635 -0.0593 -0.0537 -0.0537  0.0167  0.0167  1.6484  1.7917  1.8747
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250   ndimh = 24
 -0.3287 -0.1260 -0.0735 -0.0257  0.0507  0.0918  0.6905  1.2570  1.3819
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750   ndimh = 25
 -0.4093 -0.1098 -0.0575 -0.0370  0.0283  0.0631  0.7942  1.3624  1.6367
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250   ndimh = 23
 -0.3862 -0.1329 -0.0629  0.0028  0.0281  0.0659  0.9938  1.1757  1.3225
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750   ndimh = 24
 -0.2231 -0.1327 -0.1105  0.0003  0.0466  0.2001  0.6476  0.8623  1.1603
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250   ndimh = 25
 -0.4334 -0.1000 -0.0483 -0.0483  0.0388  0.0388  0.7743  1.6552  1.7090

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.144354;  11.000000 electrons;  D(Ef):    6.233
         Sum occ. bands:   -0.8562872  incl. Bloechl correction:   -0.006582

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250   ndimh = 26
 -0.6635 -0.0593 -0.0537 -0.0537  0.0167  0.0167  1.6484  1.7917  1.8747
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250   ndimh = 24
 -0.3287 -0.1260 -0.0735 -0.0257  0.0507  0.0918  0.6905  1.2570  1.3819
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750   ndimh = 25
 -0.4093 -0.1098 -0.0575 -0.0370  0.0283  0.0631  0.7942  1.3624  1.6367
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250   ndimh = 23
 -0.3862 -0.1329 -0.0629  0.0028  0.0281  0.0659  0.9938  1.1757  1.3225
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750   ndimh = 24
 -0.2231 -0.1327 -0.1105  0.0003  0.0466  0.2001  0.6476  0.8623  1.1603
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250   ndimh = 25
 -0.4334 -0.1000 -0.0483 -0.0483  0.0388  0.0388  0.7743  1.6552  1.7090

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.144354;  11.000000 electrons;  D(Ef):    6.233
         Sum occ. bands:   -0.8562872  incl. Bloechl correction:   -0.006582

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1    9.920371    3.088533    6.831838

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.511252   -0.393285    4.650000    4.643508    4.500000    4.643508
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.537224   -0.694162    4.340000    4.261725    4.250000    4.261725
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.487834   -0.057357    3.870000    3.850175    3.147584    3.850175
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.026673   -0.438992    4.110000    4.114394    4.102416    4.110000

 Harris energy:
 sumev=       -0.856287  val*vef=    -189.493805   sumtv=     188.637518
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -130.258096     utot=   -6534.971175    ehar=   -3304.835114

 srhov:     -6.290617   -168.128988   -174.419604 sumev=   -0.856287   sumtv=  173.563317

 Kohn-Sham energy:
 sumtv=      173.563317  sumtc=      3171.756639   ekin=     3345.319956
 rhoep=     -128.619758   utot=     -6521.212349   ehks=    -3304.512151
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 0.  RMS DQ=3.96e-2
 charges:       old           new         screened      rms diff       lin mix
 smooth       4.646654      4.168162      4.168162      0.038629      4.168162
 site    1    6.353346      6.831838      6.831838      0.019748      6.831838
 AMIX: nmix=0 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=3.96e-2
 unscreened rms difference:  smooth  0.046252   local  0.019748
   screened rms difference:  smooth  0.038629   local  0.019748   tot  0.039611

 iors  : write restart file (binary, mesh density) 

   it  1  of 12    ehf=      -0.400614   ehk=      -0.077651
h nk=8 bigbas=0 pwmode=11 ehf=-.400614 ehk=-.0776506

 --- BNDFP:  begin iteration 2 of 12 ---

 avg es pot at rmt= 0.633496  avg sphere pot= 0.655711  vconst=-0.633496

 smooth rhoves     12.210130   charge     4.168162
 smooth rhoeps =   -3.089612   rhomu =   -4.021924  avg vxc =   -0.857846 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000666

 Energy terms:             smooth           local           total
   rhoval*vef             -7.650188      -175.334084      -182.984271
   rhoval*ves            -48.918701      -107.963387      -156.882087
   psnuc*ves              73.338960    -12963.912592    -12890.573632
   utot                   12.210130     -6535.937990     -6523.727860
   rho*exc                -3.089612      -125.882370      -128.971983
   rho*vxc                -4.021924      -166.693164      -170.715088
   valence chg             4.168162         6.831838        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=0.144354,      11 occ states

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250   ndimh = 26
 -0.7796 -0.6932 -0.6910 -0.6910 -0.6538 -0.6538  1.5097  1.6804  1.7088
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250   ndimh = 24
 -0.7279 -0.6954 -0.6761 -0.6647 -0.6287 -0.3457  0.4540  1.0695  1.1959
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750   ndimh = 25
 -0.7283 -0.6892 -0.6802 -0.6598 -0.6374 -0.4629  0.5637  1.1774  1.4840
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250   ndimh = 23
 -0.7248 -0.7093 -0.6698 -0.6538 -0.6340 -0.4344  0.7857  0.9891  1.1280
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750   ndimh = 24
 -0.7100 -0.7086 -0.6873 -0.6528 -0.6340 -0.1414  0.4100  0.6558  0.9624
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250   ndimh = 25
 -0.7365 -0.6866 -0.6866 -0.6396 -0.6396 -0.4994  0.5396  1.5148  1.5148

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.195725;  11.000000 electrons;  D(Ef):    3.340
         Sum occ. bands:   -7.2052650  incl. Bloechl correction:   -0.013338

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.415859    1.977014    8.438845

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.425167   -1.245844    4.643508    4.338550    4.500000    4.500000
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.262824   -2.272188    4.261725    4.133557    4.250000    4.250000
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    9.173332   -0.671704    3.850175    3.904421    3.147584    3.904421
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.008765   -1.127750    4.110000    4.100948    4.102416    4.110000

 Harris energy:
 sumev=       -7.205265  val*vef=    -182.984271   sumtv=     175.779006
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -128.971983     utot=   -6523.727860    ehar=   -3305.164197

 srhov:     -4.446213   -222.346955   -226.793168 sumev=   -7.205265   sumtv=  219.587903

 Kohn-Sham energy:
 sumtv=      219.587903  sumtc=      3171.756639   ekin=     3391.344543
 rhoep=     -132.709563   utot=     -6562.367422   ehks=    -3303.732442
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 1.  RMS DQ=9.95e-2  last it=3.96e-2
 mixrho: (warning) scr. and lin-mixed densities had 43 and 43 negative points
 charges:       old           new         screened      rms diff       lin mix
 smooth       4.168162      2.561155      2.561155      0.030736      2.561155
 site    1    6.831838      8.438845      8.438845      0.064535      8.438845
 AMIX: nmix=1 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=9.95e-2
   tj: 0.81636
 unscreened rms difference:  smooth  0.022300   local  0.064535
   screened rms difference:  smooth  0.030736   local  0.064535   tot  0.099549

 iors  : write restart file (binary, mesh density) 

   it  2  of 12    ehf=      -0.729697   ehk=       0.702058
 From last iter    ehf=      -0.400614   ehk=      -0.077651
 diffe(q)= -0.329083 (0.099549)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 pwmode=11 ehf=-.7296971 ehk=.7020575

 --- BNDFP:  begin iteration 3 of 12 ---

 avg es pot at rmt= 0.605280  avg sphere pot= 0.656289  vconst=-0.605280

 smooth rhoves     10.912132   charge     3.873052
 smooth rhoeps =   -2.806612   rhomu =   -3.652533  avg vxc =   -0.843026 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000604

 Energy terms:             smooth           local           total
   rhoval*vef             -6.541163      -181.321721      -187.862884
   rhoval*ves            -47.623901      -113.556369      -161.180270
   psnuc*ves              69.448165    -12968.690897    -12899.242732
   utot                   10.912132     -6541.123633     -6530.211501
   rho*exc                -2.806612      -126.715169      -129.521782
   rho*vxc                -3.652533      -167.791978      -171.444510
   valence chg             3.873052         7.126948        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.195725,      11 occ states

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250   ndimh = 26
 -0.7485 -0.4489 -0.4456 -0.4456 -0.3974 -0.3974  1.5525  1.7077  1.7611
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250   ndimh = 24
 -0.5245 -0.4584 -0.4444 -0.4216 -0.3676 -0.2483  0.5222  1.1254  1.2521
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750   ndimh = 25
 -0.5499 -0.4461 -0.4410 -0.4309 -0.3820 -0.3329  0.6310  1.2331  1.5312
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250   ndimh = 23
 -0.5309 -0.4908 -0.4334 -0.4012 -0.3767 -0.3170  0.8478  1.0443  1.1866
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750   ndimh = 24
 -0.4871 -0.4768 -0.4627 -0.4009 -0.3737 -0.0601  0.4778  0.7164  1.0213
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250   ndimh = 25
 -0.5687 -0.4406 -0.4406 -0.4034 -0.3801 -0.3801  0.6079  1.5739  1.5739

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.122850;  11.000000 electrons;  D(Ef):    3.479
         Sum occ. bands:   -4.7291125  incl. Bloechl correction:   -0.011768

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.289073    2.332539    7.956534

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.449764   -1.026687    4.500000    4.401921    4.500000    4.500000
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.308764   -1.763331    4.250000    4.157315    4.250000    4.250000
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    9.057848   -0.430792    3.904421    3.887228    3.147584    3.887228
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.013043   -0.884267    4.110000    4.105490    4.102416    4.110000

 Harris energy:
 sumev=       -4.729113  val*vef=    -187.862884   sumtv=     183.133772
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.521782     utot=   -6530.211501    ehar=   -3304.842871

 srhov:     -5.028453   -202.674029   -207.702482 sumev=   -4.729113   sumtv=  202.973369

 Kohn-Sham energy:
 sumtv=      202.973369  sumtc=      3171.756639   ekin=     3374.730009
 rhoep=     -131.308450   utot=     -6547.927495   ehks=    -3304.505936
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 2.  RMS DQ=4.68e-2  last it=9.95e-2
 mixrho: (warning) scr. and lin-mixed densities had 13 and 13 negative points
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.873052      3.043467      3.043467      0.015144      3.043467
 site    1    7.126948      7.956534      7.956534      0.031578      7.956534
 AMIX: nmix=2 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=4.68e-2
   tj:-0.97579  -0.09921
 add q= -0.000002 to preserve neutrality
 unscreened rms difference:  smooth  0.011194   local  0.031578
   screened rms difference:  smooth  0.015144   local  0.031578   tot  0.046769

 iors  : write restart file (binary, mesh density) 

   it  3  of 12    ehf=      -0.408371   ehk=      -0.071436
 From last iter    ehf=      -0.729697   ehk=       0.702058
 diffe(q)=  0.321326 (0.046769)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 pwmode=11 ehf=-.4083714 ehk=-.071436

 --- BNDFP:  begin iteration 4 of 12 ---

 avg es pot at rmt= 0.548463  avg sphere pot= 0.665711  vconst=-0.548463

 smooth rhoves      8.819425   charge     3.402518
 smooth rhoeps =   -2.385959   rhomu =   -3.103848  avg vxc =   -0.815174 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000523

 Energy terms:             smooth           local           total
   rhoval*vef             -5.055601      -185.839876      -190.895477
   rhoval*ves            -44.878427      -118.653047      -163.531474
   psnuc*ves              62.517277    -12970.128390    -12907.611113
   utot                    8.819425     -6544.390719     -6535.571294
   rho*exc                -2.385959      -127.770842      -130.156801
   rho*vxc                -3.103848      -169.182394      -172.286242
   valence chg             3.402518         7.597482        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.12285,      11 occ states

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250   ndimh = 26
 -0.6952 -0.1457 -0.1406 -0.1406 -0.0744 -0.0744  1.6148  1.7597  1.8386
 Est Ef = -0.123 < evl(5)=-0.074 ... using qval=11.0, revise to -0.0744
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250   ndimh = 24
 -0.3699 -0.1997 -0.1570 -0.1131 -0.0406  0.0048  0.6412  1.2164  1.3419
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750   ndimh = 25
 -0.4451 -0.1858 -0.1430 -0.1242 -0.0621 -0.0283  0.7464  1.3224  1.6025
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250   ndimh = 23
 -0.4222 -0.2127 -0.1453 -0.0860 -0.0609 -0.0253  0.9502  1.1350  1.2813
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750   ndimh = 24
 -0.2758 -0.2078 -0.1896 -0.0877 -0.0456  0.1267  0.5976  0.8183  1.1180
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250   ndimh = 25
 -0.4688 -0.1738 -0.1351 -0.1351 -0.0531 -0.0531  0.7257  1.6301  1.6669

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.067639;  11.000000 electrons;  D(Ef):    5.297
         Sum occ. bands:   -1.7149650  incl. Bloechl correction:   -0.007519

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.002063    2.943059    7.059004

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.498884   -0.611102    4.500000    4.558676    4.500000    4.558676
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.478746   -0.892479    4.250000    4.234594    4.250000    4.250000
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.626651   -0.139850    3.887228    3.857911    3.147584    3.857911
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.023332   -0.535110    4.110000    4.112639    4.102416    4.110000

 Harris energy:
 sumev=       -1.714965  val*vef=    -190.895477   sumtv=     189.180512
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -130.156801     utot=   -6535.571294    ehar=   -3304.790943

 srhov:     -5.852602   -174.408015   -180.260617 sumev=   -1.714965   sumtv=  178.545652

 Kohn-Sham energy:
 sumtv=      178.545652  sumtc=      3171.756639   ekin=     3350.302291
 rhoep=     -129.102487   utot=     -6525.864492   ehks=    -3304.664687
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 3.  RMS DQ=2.68e-2  last it=4.68e-2
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.402518      3.940996      3.940996      0.008943      3.940996
 site    1    7.597482      7.059004      7.059004      0.019397      7.059004
 AMIX: nmix=3 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=2.68e-2
   tj: 0.74666  -0.22381  -0.00259
 unscreened rms difference:  smooth  0.006933   local  0.019397
   screened rms difference:  smooth  0.008943   local  0.019397   tot  0.026809

 iors  : write restart file (binary, mesh density) 

   it  4  of 12    ehf=      -0.356443   ehk=      -0.230187
 From last iter    ehf=      -0.408371   ehk=      -0.071436
 diffe(q)=  0.051928 (0.026809)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 pwmode=11 ehf=-.3564435 ehk=-.2301871

 --- BNDFP:  begin iteration 5 of 12 ---

 avg es pot at rmt= 0.572536  avg sphere pot= 0.661016  vconst=-0.572536

 smooth rhoves      9.619935   charge     3.579077
 smooth rhoeps =   -2.538258   rhomu =   -3.302430  avg vxc =   -0.826515 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000550

 Energy terms:             smooth           local           total
   rhoval*vef             -5.562388      -184.845273      -190.407661
   rhoval*ves            -46.045499      -117.228168      -163.273667
   psnuc*ves              65.285369    -12970.422180    -12905.136811
   utot                    9.619935     -6543.825174     -6534.205239
   rho*exc                -2.538258      -127.406138      -129.944396
   rho*vxc                -3.302430      -168.702434      -172.004865
   valence chg             3.579077         7.420923        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=0.067639,      11 occ states

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250   ndimh = 26
 -0.7172 -0.2536 -0.2492 -0.2492 -0.1898 -0.1898  1.5905  1.7377  1.8087
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250   ndimh = 24
 -0.4134 -0.2891 -0.2609 -0.2231 -0.1570 -0.0988  0.5933  1.1804  1.3065
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750   ndimh = 25
 -0.4764 -0.2781 -0.2501 -0.2337 -0.1766 -0.1430  0.7002  1.2872  1.5747
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250   ndimh = 23
 -0.4545 -0.3108 -0.2479 -0.1983 -0.1732 -0.1385  0.9096  1.0990  1.2439
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750   ndimh = 24
 -0.3402 -0.3008 -0.2871 -0.1991 -0.1628  0.0472  0.5492  0.7776  1.0796
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250   ndimh = 25
 -0.4993 -0.2616 -0.2438 -0.2438 -0.1700 -0.1700  0.6784  1.6163  1.6304

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.015646;  11.000000 electrons;  D(Ef):    4.299
         Sum occ. bands:   -2.7772164  incl. Bloechl correction:   -0.009084

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.114092    2.712353    7.401739

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.477158   -0.747344    4.558676    4.507301    4.500000    4.507301
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.412393   -1.163493    4.250000    4.203525    4.250000    4.250000
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.807253   -0.243920    3.857911    3.867253    3.147584    3.867253
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.019113   -0.656889    4.110000    4.110099    4.102416    4.110000

 Harris energy:
 sumev=       -2.777216  val*vef=    -190.407661   sumtv=     187.630444
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.944396     utot=   -6534.205239    ehar=   -3304.762552

 srhov:     -5.569821   -184.383123   -189.952944 sumev=   -2.777216   sumtv=  187.175727

 Kohn-Sham energy:
 sumtv=      187.175727  sumtc=      3171.756639   ekin=     3358.932367
 rhoep=     -129.895252   utot=     -6533.799376   ehks=    -3304.762262
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=1.05e-3  last it=2.68e-2
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.579077      3.598261      3.598261      0.000253      3.598261
 site    1    7.420923      7.401739      7.401739      0.000746      7.401739
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: nmix=2 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=1.05e-3
   tj:-0.04345  -0.00165
 unscreened rms difference:  smooth  0.000325   local  0.000746
   screened rms difference:  smooth  0.000253   local  0.000746   tot  0.001047

 iors  : write restart file (binary, mesh density) 

   it  5  of 12    ehf=      -0.328052   ehk=      -0.327762
 From last iter    ehf=      -0.356443   ehk=      -0.230187
 diffe(q)=  0.028392 (0.001047)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 pwmode=11 ehf=-.3280516 ehk=-.3277616

 --- BNDFP:  begin iteration 6 of 12 ---

 avg es pot at rmt= 0.574125  avg sphere pot= 0.661028  vconst=-0.574125

 smooth rhoves      9.657847   charge     3.584283
 smooth rhoeps =   -2.541609   rhomu =   -3.306785  avg vxc =   -0.827057 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000548

 Energy terms:             smooth           local           total
   rhoval*vef             -5.567177      -184.835453      -190.402630
   rhoval*ves            -46.106796      -117.177098      -163.283893
   psnuc*ves              65.422491    -12970.406597    -12904.984106
   utot                    9.657847     -6543.791847     -6534.134000
   rho*exc                -2.541609      -127.388811      -129.930419
   rho*vxc                -3.306785      -168.679579      -171.986364
   valence chg             3.584283         7.415717        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.015646,      11 occ states

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250   ndimh = 26
 -0.7190 -0.2617 -0.2573 -0.2573 -0.1983 -0.1983  1.5884  1.7359  1.8063
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250   ndimh = 24
 -0.4173 -0.2959 -0.2686 -0.2313 -0.1657 -0.1061  0.5897  1.1776  1.3037
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750   ndimh = 25
 -0.4792 -0.2850 -0.2581 -0.2418 -0.1850 -0.1513  0.6968  1.2845  1.5724
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250   ndimh = 23
 -0.4573 -0.3183 -0.2556 -0.2066 -0.1815 -0.1466  0.9065  1.0962  1.2411
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750   ndimh = 24
 -0.3456 -0.3079 -0.2944 -0.2073 -0.1715  0.0417  0.5456  0.7745  1.0767
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250   ndimh = 25
 -0.5020 -0.2682 -0.2519 -0.2519 -0.1786 -0.1786  0.6749  1.6150  1.6276

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.021311;  11.000000 electrons;  D(Ef):    4.246
         Sum occ. bands:   -2.8573330  incl. Bloechl correction:   -0.009199

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.123068    2.695936    7.427132

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.475457   -0.764056    4.507301    4.500440    4.500000    4.500440
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.405454   -1.192212    4.250000    4.200648    4.250000    4.250000
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.824672   -0.251315    3.867253    3.868374    3.147584    3.868374
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.018781   -0.667202    4.110000    4.109891    4.102416    4.110000

 Harris energy:
 sumev=       -2.857333  val*vef=    -190.402630   sumtv=     187.545297
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.930419     utot=   -6534.134000    ehar=   -3304.762483

 srhov:     -5.548565   -185.124723   -190.673288 sumev=   -2.857333   sumtv=  187.815955

 Kohn-Sham energy:
 sumtv=      187.815955  sumtc=      3171.756639   ekin=     3359.572594
 rhoep=     -129.955462   utot=     -6534.379547   ehks=    -3304.762415
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=6.46e-4  last it=1.05e-3
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.584283      3.572868      3.572868      0.000203      3.572868
 site    1    7.415717      7.427132      7.427132      0.000441      7.427132
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: nmix=2 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=6.46e-4
   tj: 0.31816   0.00387
 unscreened rms difference:  smooth  0.000154   local  0.000441
   screened rms difference:  smooth  0.000203   local  0.000441   tot  0.000646

 iors  : write restart file (binary, mesh density) 

   it  6  of 12    ehf=      -0.327983   ehk=      -0.327915
 From last iter    ehf=      -0.328052   ehk=      -0.327762
 diffe(q)=  0.000069 (0.000646)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 pwmode=11 ehf=-.327983 ehk=-.3279148

 --- BNDFP:  begin iteration 7 of 12 ---

 avg es pot at rmt= 0.573881  avg sphere pot= 0.661020  vconst=-0.573881

 smooth rhoves      9.649431   charge     3.582372
 smooth rhoeps =   -2.539888   rhomu =   -3.304541  avg vxc =   -0.826948 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000547

 Energy terms:             smooth           local           total
   rhoval*vef             -5.561189      -184.863867      -190.425056
   rhoval*ves            -46.095459      -117.206936      -163.302395
   psnuc*ves              65.394321    -12970.430092    -12905.035772
   utot                    9.649431     -6543.818514     -6534.169083
   rho*exc                -2.539888      -127.394203      -129.934092
   rho*vxc                -3.304541      -168.686693      -171.991233
   valence chg             3.582372         7.417628        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.021311,      11 occ states

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250   ndimh = 26
 -0.7188 -0.2600 -0.2557 -0.2557 -0.1966 -0.1966  1.5887  1.7361  1.8067
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250   ndimh = 24
 -0.4165 -0.2945 -0.2670 -0.2296 -0.1639 -0.1046  0.5904  1.1781  1.3042
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750   ndimh = 25
 -0.4787 -0.2836 -0.2564 -0.2402 -0.1833 -0.1496  0.6974  1.2849  1.5728
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250   ndimh = 23
 -0.4568 -0.3168 -0.2540 -0.2049 -0.1798 -0.1450  0.9071  1.0967  1.2416
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750   ndimh = 24
 -0.3445 -0.3065 -0.2929 -0.2056 -0.1697  0.0428  0.5463  0.7751  1.0772
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250   ndimh = 25
 -0.5015 -0.2669 -0.2502 -0.2502 -0.1768 -0.1768  0.6756  1.6152  1.6281

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.020183;  11.000000 electrons;  D(Ef):    4.257
         Sum occ. bands:   -2.8410017  incl. Bloechl correction:   -0.009174

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.121498    2.699506    7.421992

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.475688   -0.763079    4.500440    4.500719    4.500000    4.500719
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.406269   -1.188328    4.250000    4.201025    4.250000    4.250000
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.822522   -0.249705    3.868374    3.868247    3.147584    3.868247
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.018838   -0.665401    4.110000    4.109929    4.102416    4.110000

 Harris energy:
 sumev=       -2.841002  val*vef=    -190.425056   sumtv=     187.584054
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.934092     utot=   -6534.169083    ehar=   -3304.762481

 srhov:     -5.554075   -184.965323   -190.519398 sumev=   -2.841002   sumtv=  187.678396

 Kohn-Sham energy:
 sumtv=      187.678396  sumtc=      3171.756639   ekin=     3359.435036
 rhoep=     -129.943128   utot=     -6534.254380   ehks=    -3304.762472
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=2.33e-4  last it=6.46e-4
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.582372      3.578008      3.578008      0.000077      3.578008
 site    1    7.417628      7.421992      7.421992      0.000164      7.421992
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: nmix=2 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=2.33e-4
   tj:-0.26035   0.09822
 unscreened rms difference:  smooth  0.000060   local  0.000164
   screened rms difference:  smooth  0.000077   local  0.000164   tot  0.000233

 iors  : write restart file (binary, mesh density) 

   it  7  of 12    ehf=      -0.327981   ehk=      -0.327972
 From last iter    ehf=      -0.327983   ehk=      -0.327915
 diffe(q)=  0.000002 (0.000233)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 pwmode=11 ehf=-.3279811 ehk=-.327972

 --- BNDFP:  begin iteration 8 of 12 ---

 avg es pot at rmt= 0.573741  avg sphere pot= 0.661034  vconst=-0.573741

 smooth rhoves      9.644773   charge     3.581336
 smooth rhoeps =   -2.538959   rhomu =   -3.303329  avg vxc =   -0.826889 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000547

 Energy terms:             smooth           local           total
   rhoval*vef             -5.558017      -184.871199      -190.429216
   rhoval*ves            -46.089126      -117.215917      -163.305043
   psnuc*ves              65.378672    -12970.431256    -12905.052584
   utot                    9.644773     -6543.823587     -6534.178814
   rho*exc                -2.538959      -127.396536      -129.935496
   rho*vxc                -3.303329      -168.689765      -171.993094
   valence chg             3.581336         7.418664        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.020183,      11 occ states

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250   ndimh = 26
 -0.7186 -0.2593 -0.2549 -0.2549 -0.1958 -0.1958  1.5889  1.7363  1.8069
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250   ndimh = 24
 -0.4162 -0.2939 -0.2663 -0.2289 -0.1631 -0.1040  0.5907  1.1783  1.3045
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750   ndimh = 25
 -0.4785 -0.2830 -0.2557 -0.2394 -0.1825 -0.1489  0.6977  1.2852  1.5730
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250   ndimh = 23
 -0.4566 -0.3161 -0.2533 -0.2042 -0.1790 -0.1442  0.9074  1.0969  1.2418
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750   ndimh = 24
 -0.3440 -0.3058 -0.2923 -0.2049 -0.1689  0.0433  0.5466  0.7753  1.0774
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250   ndimh = 25
 -0.5013 -0.2663 -0.2495 -0.2495 -0.1761 -0.1761  0.6759  1.6153  1.6283

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.019696;  11.000000 electrons;  D(Ef):    4.261
         Sum occ. bands:   -2.8340298  incl. Bloechl correction:   -0.009164

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.120780    2.700976    7.419805

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.475814   -0.762010    4.500719    4.501141    4.500000    4.501141
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.406724   -1.186299    4.250000    4.201224    4.250000    4.250000
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.821379   -0.249037    3.868247    3.868173    3.147584    3.868173
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.018863   -0.664581    4.110000    4.109946    4.102416    4.110000

 Harris energy:
 sumev=       -2.834030  val*vef=    -190.429216   sumtv=     187.595186
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.935496     utot=   -6534.178814    ehar=   -3304.762484

 srhov:     -5.556047   -184.899070   -190.455117 sumev=   -2.834030   sumtv=  187.621087

 Kohn-Sham energy:
 sumtv=      187.621087  sumtc=      3171.756639   ekin=     3359.377727
 rhoep=     -129.937881   utot=     -6534.202328   ehks=    -3304.762483
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=6.25e-5  last it=2.33e-4
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.581336      3.580195      3.580195      0.000022      3.580195
 site    1    7.418664      7.419805      7.419805      0.000043      7.419805
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=6.25e-5
   tj:-0.36599
 unscreened rms difference:  smooth  0.000017   local  0.000043
   screened rms difference:  smooth  0.000022   local  0.000043   tot  0.000063

 iors  : write restart file (binary, mesh density) 

   it  8  of 12    ehf=      -0.327984   ehk=      -0.327983
 From last iter    ehf=      -0.327981   ehk=      -0.327972
 diffe(q)= -0.000002 (0.000063)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=0 pwmode=11 ehf=-.3279836 ehk=-.3279828

 --- BNDFP:  begin iteration 9 of 12 ---

 avg es pot at rmt= 0.573699  avg sphere pot= 0.661038  vconst=-0.573699

 smooth rhoves      9.643333   charge     3.580996
 smooth rhoeps =   -2.538644   rhomu =   -3.302918  avg vxc =   -0.826872 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=16  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000547

 Energy terms:             smooth           local           total
   rhoval*vef             -5.556905      -184.875237      -190.432142
   rhoval*ves            -46.087229      -117.220209      -163.307438
   psnuc*ves              65.373896    -12970.433383    -12905.059486
   utot                    9.643333     -6543.826796     -6534.183462
   rho*exc                -2.538644      -127.397349      -129.935993
   rho*vxc                -3.302918      -168.690836      -171.993754
   valence chg             3.580996         7.419004        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.019696,      11 occ states

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250   ndimh = 26
 -0.7186 -0.2591 -0.2547 -0.2547 -0.1956 -0.1956  1.5889  1.7363  1.8070
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250   ndimh = 24
 -0.4161 -0.2937 -0.2661 -0.2286 -0.1629 -0.1038  0.5908  1.1784  1.3045
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750   ndimh = 25
 -0.4785 -0.2828 -0.2555 -0.2392 -0.1823 -0.1486  0.6978  1.2852  1.5730
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250   ndimh = 23
 -0.4566 -0.3159 -0.2531 -0.2039 -0.1788 -0.1440  0.9074  1.0970  1.2419
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750   ndimh = 24
 -0.3439 -0.3056 -0.2921 -0.2047 -0.1687  0.0434  0.5467  0.7754  1.0775
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250   ndimh = 25
 -0.5012 -0.2661 -0.2493 -0.2493 -0.1758 -0.1758  0.6760  1.6153  1.6284

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.019542;  11.000000 electrons;  D(Ef):    4.263
         Sum occ. bands:   -2.8318111  incl. Bloechl correction:   -0.009160

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.120543    2.701467    7.419076

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.475863   -0.761653    4.501141    4.501287    4.500000    4.501287
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.406876   -1.185644    4.250000    4.201289    4.250000    4.250000
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.820992   -0.248824    3.868173    3.868148    3.147584    3.868148
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.018872   -0.664319    4.110000    4.109951    4.102416    4.110000

 Harris energy:
 sumev=       -2.831811  val*vef=    -190.432142   sumtv=     187.600331
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.935993     utot=   -6534.183462    ehar=   -3304.762485

 srhov:     -5.556739   -184.876415   -190.433154 sumev=   -2.831811   sumtv=  187.601343

 Kohn-Sham energy:
 sumtv=      187.601343  sumtc=      3171.756639   ekin=     3359.357983
 rhoep=     -129.936096   utot=     -6534.184372   ehks=    -3304.762485
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=3.19e-6  last it=6.25e-5
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.580996      3.580924      3.580924      0.000006      3.580924
 site    1    7.419004      7.419076      7.419076      0.000003      7.419076
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=2016  beta=1  tm=5  rmsdel=3.19e-6
   tj:-0.05155
 unscreened rms difference:  smooth  0.000006   local  0.000003
   screened rms difference:  smooth  0.000006   local  0.000003   tot  0.000003

 iors  : write restart file (binary, mesh density) 

   it  9  of 12    ehf=      -0.327985   ehk=      -0.327985
 From last iter    ehf=      -0.327984   ehk=      -0.327983
 diffe(q)= -0.000002 (0.000003)    tol= 0.000010 (0.000010)   more=F
c nk=8 bigbas=0 pwmode=11 ehf=-.3279851 ehk=-.3279849
 Exit 0 LMF 
 CPU time:    3.118s   Wall clock    3.166s   12:04:04 02.10.2014        on phpdl1.ph.kcl.ac.uk
rdcmd:  rm mixm.cu
rdcmd:  lmf  --no-iactiv cu -vnk=8 -vbigbas=t -vpwmode=11
 -----------------------  START LMF -----------------------

 LMF:      nbas = 1  nspec = 1  vn 7.11.c  verb 31,30,31
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
         site   1, species A       : augmentation lmax changed from 3 to 4
         site   1, species A       : inflate local density from nlm= 16 to 25

 --- BNDFP:  begin iteration 1 of 12 ---

 avg es pot at rmt= 0.573696  avg sphere pot= 0.661039  vconst=-0.573696

 smooth rhoves      9.643207   charge     3.580961
 smooth rhoeps =   -2.538610   rhomu =   -3.302874  avg vxc =   -0.826870 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000547

 Energy terms:             smooth           local           total
   rhoval*vef             -5.556777      -184.875315      -190.432092
   rhoval*ves            -46.087078      -117.220296      -163.307373
   psnuc*ves              65.373492    -12970.433052    -12905.059561
   utot                    9.643207     -6543.826674     -6534.183467
   rho*exc                -2.538610      -127.397396      -129.936007
   rho*vxc                -3.302874      -168.690898      -171.993772
   valence chg             3.580961         7.419039        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250   ndimh = 39
 -0.7186 -0.2594 -0.2550 -0.2550 -0.1954 -0.1954  1.5834  1.7294  1.8025
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250   ndimh = 37
 -0.4163 -0.2940 -0.2663 -0.2287 -0.1629 -0.1039  0.5899  1.1746  1.3012
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750   ndimh = 38
 -0.4785 -0.2832 -0.2556 -0.2393 -0.1824 -0.1487  0.6968  1.2826  1.5675
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250   ndimh = 36
 -0.4567 -0.3163 -0.2532 -0.2041 -0.1789 -0.1441  0.9058  1.0943  1.2373
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750   ndimh = 37
 -0.3440 -0.3058 -0.2923 -0.2046 -0.1689  0.0432  0.5458  0.7732  1.0753
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250   ndimh = 38
 -0.5013 -0.2665 -0.2494 -0.2494 -0.1759 -0.1759  0.6753  1.6069  1.6257

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.019848;  11.000000 electrons;  D(Ef):    4.268
         Sum occ. bands:   -2.8334005  incl. Bloechl correction:   -0.009148

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250   ndimh = 39
 -0.7186 -0.2594 -0.2550 -0.2550 -0.1954 -0.1954  1.5834  1.7294  1.8025
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250   ndimh = 37
 -0.4163 -0.2940 -0.2663 -0.2287 -0.1629 -0.1039  0.5899  1.1746  1.3012
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750   ndimh = 38
 -0.4785 -0.2832 -0.2556 -0.2393 -0.1824 -0.1487  0.6968  1.2826  1.5675
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250   ndimh = 36
 -0.4567 -0.3163 -0.2532 -0.2041 -0.1789 -0.1441  0.9058  1.0943  1.2373
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750   ndimh = 37
 -0.3440 -0.3058 -0.2923 -0.2046 -0.1689  0.0432  0.5458  0.7732  1.0753
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250   ndimh = 38
 -0.5013 -0.2665 -0.2494 -0.2494 -0.1759 -0.1759  0.6753  1.6069  1.6257

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.019848;  11.000000 electrons;  D(Ef):    4.268
         Sum occ. bands:   -2.8334005  incl. Bloechl correction:   -0.009148

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.131064    2.842339    7.288725

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.476042   -1.597485    4.501287    4.238970    4.500000    4.500000
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.410206   -1.212495    4.250000    4.198432    4.250000    4.250000
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    9.012827   -0.269145    3.868148    3.842355    3.147584    3.842355
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.021130   -0.732288    4.110000    4.108384    4.102416    4.110000
 4     1    0.004569   -0.633811    5.100000    5.080968    5.077979    5.100000

 Harris energy:
 sumev=       -2.833401  val*vef=    -190.432092   sumtv=     187.598691
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.936007     utot=   -6534.183467    ehar=   -3304.764143

 srhov:     -5.913421   -184.362771   -190.276192 sumev=   -2.833401   sumtv=  187.442791

 Kohn-Sham energy:
 sumtv=      187.442791  sumtc=      3171.756639   ekin=     3359.199431
 rhoep=     -129.926887   utot=     -6534.036629   ehks=    -3304.764085
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 0.  RMS DQ=3.93e-3
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.580961      3.711275      3.711275      0.003076      3.711275
 site    1    7.419039      7.288725      7.288725      0.004512      7.288725
 AMIX: nmix=0 mmix=8  nelts=2178  beta=1  tm=5  rmsdel=3.93e-3
 unscreened rms difference:  smooth  0.003083   local  0.004512
   screened rms difference:  smooth  0.003076   local  0.004512   tot  0.003930

 iors  : write restart file (binary, mesh density) 

   it  1  of 12    ehf=      -0.329643   ehk=      -0.329585
i nk=8 bigbas=1 pwmode=11 ehf=-.3296428 ehk=-.329585

 --- BNDFP:  begin iteration 2 of 12 ---

 avg es pot at rmt= 0.582910  avg sphere pot= 0.651662  vconst=-0.582910

 smooth rhoves     10.184194   charge     3.711275
 smooth rhoeps =   -2.664253   rhomu =   -3.466874  avg vxc =   -0.833129 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000129

 Energy terms:             smooth           local           total
   rhoval*vef             -5.986583      -184.312261      -190.298844
   rhoval*ves            -46.700612      -116.482732      -163.183344
   psnuc*ves              67.069000    -12971.961576    -12904.892576
   utot                   10.184194     -6544.222154     -6534.037960
   rho*exc                -2.664253      -127.262804      -129.927057
   rho*vxc                -3.466874      -168.514999      -171.981873
   valence chg             3.711275         7.288725        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250   ndimh = 39
 -0.7183 -0.2621 -0.2578 -0.2578 -0.1970 -0.1970  1.5835  1.7296  1.8024
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250   ndimh = 37
 -0.4169 -0.2958 -0.2683 -0.2309 -0.1654 -0.1059  0.5893  1.1744  1.3012
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750   ndimh = 38
 -0.4786 -0.2853 -0.2576 -0.2415 -0.1847 -0.1510  0.6962  1.2824  1.5676
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250   ndimh = 36
 -0.4567 -0.3187 -0.2551 -0.2066 -0.1810 -0.1463  0.9055  1.0941  1.2372
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750   ndimh = 37
 -0.3452 -0.3076 -0.2945 -0.2063 -0.1719  0.0419  0.5452  0.7729  1.0752
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250   ndimh = 38
 -0.5014 -0.2685 -0.2516 -0.2516 -0.1782 -0.1782  0.6747  1.6074  1.6255

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.021190;  11.000000 electrons;  D(Ef):    4.248
         Sum occ. bands:   -2.8545629  incl. Bloechl correction:   -0.009193

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250   ndimh = 39
 -0.7183 -0.2621 -0.2578 -0.2578 -0.1970 -0.1970  1.5835  1.7296  1.8024
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250   ndimh = 37
 -0.4169 -0.2958 -0.2683 -0.2309 -0.1654 -0.1059  0.5893  1.1744  1.3012
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750   ndimh = 38
 -0.4786 -0.2853 -0.2576 -0.2415 -0.1847 -0.1510  0.6962  1.2824  1.5676
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250   ndimh = 36
 -0.4567 -0.3187 -0.2551 -0.2066 -0.1810 -0.1463  0.9055  1.0941  1.2372
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750   ndimh = 37
 -0.3452 -0.3076 -0.2945 -0.2063 -0.1719  0.0419  0.5452  0.7729  1.0752
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250   ndimh = 38
 -0.5014 -0.2685 -0.2516 -0.2516 -0.1782 -0.1782  0.6747  1.6074  1.6255

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.021190;  11.000000 electrons;  D(Ef):    4.248
         Sum occ. bands:   -2.8545629  incl. Bloechl correction:   -0.009193

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.134247    2.834852    7.299395

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.475171   -1.592521    4.500000    4.239836    4.500000    4.500000
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.408042   -1.228281    4.250000    4.196849    4.250000    4.250000
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    9.001079   -0.271574    3.842355    3.841711    3.147584    3.841711
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.020995   -0.734364    4.110000    4.108345    4.102416    4.110000
 4     1    0.004555   -0.632155    5.100000    5.080990    5.077979    5.100000

 Harris energy:
 sumev=       -2.854563  val*vef=    -190.298844   sumtv=     187.444281
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.927057     utot=   -6534.037960    ehar=   -3304.764097

 srhov:     -5.971751   -184.660061   -190.631812 sumev=   -2.854563   sumtv=  187.777249

 Kohn-Sham energy:
 sumtv=      187.777249  sumtc=      3171.756639   ekin=     3359.533888
 rhoep=     -129.955614   utot=     -6534.342285   ehks=    -3304.764011
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 1.  RMS DQ=7.01e-4  last it=3.93e-3
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.711275      3.700604      3.700604      0.000213      3.700604
 site    1    7.288725      7.299395      7.299395      0.000455      7.299395
 AMIX: nmix=1 mmix=8  nelts=2178  beta=1  tm=5  rmsdel=7.01e-4
   tj: 0.06971
 unscreened rms difference:  smooth  0.000177   local  0.000455
   screened rms difference:  smooth  0.000213   local  0.000455   tot  0.000701

 iors  : write restart file (binary, mesh density) 

   it  2  of 12    ehf=      -0.329597   ehk=      -0.329511
 From last iter    ehf=      -0.329643   ehk=      -0.329585
 diffe(q)=  0.000046 (0.000701)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=1 pwmode=11 ehf=-.3295969 ehk=-.3295106

 --- BNDFP:  begin iteration 3 of 12 ---

 avg es pot at rmt= 0.581789  avg sphere pot= 0.651631  vconst=-0.581789

 smooth rhoves     10.140811   charge     3.701348
 smooth rhoeps =   -2.655326   rhomu =   -3.455229  avg vxc =   -0.832543 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000135

 Energy terms:             smooth           local           total
   rhoval*vef             -5.953642      -184.527431      -190.481073
   rhoval*ves            -46.646347      -116.696150      -163.342497
   psnuc*ves              66.927968    -12972.160305    -12905.232336
   utot                   10.140811     -6544.428227     -6534.287416
   rho*exc                -2.655326      -127.293430      -129.948756
   rho*vxc                -3.455229      -168.555429      -172.010658
   valence chg             3.701348         7.298652        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250   ndimh = 39
 -0.7170 -0.2531 -0.2487 -0.2487 -0.1874 -0.1874  1.5852  1.7309  1.8045
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250   ndimh = 37
 -0.4131 -0.2884 -0.2596 -0.2217 -0.1556 -0.0978  0.5928  1.1769  1.3037
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750   ndimh = 38
 -0.4762 -0.2777 -0.2487 -0.2324 -0.1752 -0.1415  0.6996  1.2849  1.5695
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250   ndimh = 36
 -0.4542 -0.3105 -0.2465 -0.1972 -0.1717 -0.1371  0.9084  1.0966  1.2399
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750   ndimh = 37
 -0.3395 -0.2998 -0.2863 -0.1970 -0.1620  0.0478  0.5487  0.7758  1.0779
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250   ndimh = 38
 -0.4991 -0.2614 -0.2425 -0.2425 -0.1684 -0.1684  0.6782  1.6081  1.6280

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.015029;  11.000000 electrons;  D(Ef):    4.314
         Sum occ. bands:   -2.7657482  incl. Bloechl correction:   -0.009052

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250   ndimh = 39
 -0.7170 -0.2531 -0.2487 -0.2487 -0.1874 -0.1874  1.5852  1.7309  1.8045
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250   ndimh = 37
 -0.4131 -0.2884 -0.2596 -0.2217 -0.1556 -0.0978  0.5928  1.1769  1.3037
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750   ndimh = 38
 -0.4762 -0.2777 -0.2487 -0.2324 -0.1752 -0.1415  0.6996  1.2849  1.5695
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250   ndimh = 36
 -0.4542 -0.3105 -0.2465 -0.1972 -0.1717 -0.1371  0.9084  1.0966  1.2399
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750   ndimh = 37
 -0.3395 -0.2998 -0.2863 -0.1970 -0.1620  0.0478  0.5487  0.7758  1.0779
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250   ndimh = 38
 -0.4991 -0.2614 -0.2425 -0.2425 -0.1684 -0.1684  0.6782  1.6081  1.6280

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.015029;  11.000000 electrons;  D(Ef):    4.314
         Sum occ. bands:   -2.7657482  incl. Bloechl correction:   -0.009052

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.124894    2.853751    7.271143

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.476887   -1.580398    4.500000    4.241582    4.500000    4.500000
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.413612   -1.208301    4.250000    4.198733    4.250000    4.250000
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.982891   -0.262803    3.841711    3.841141    3.147584    3.841141
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.021377   -0.725315    4.110000    4.108534    4.102416    4.110000
 4     1    0.004616   -0.622226    5.100000    5.081086    5.077979    5.100000

 Harris energy:
 sumev=       -2.765748  val*vef=    -190.481073   sumtv=     187.715325
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.948756     utot=   -6534.287416    ehar=   -3304.764208

 srhov:     -5.998862   -183.726199   -189.725062 sumev=   -2.765748   sumtv=  186.959313

 Kohn-Sham energy:
 sumtv=      186.959313  sumtc=      3171.756639   ekin=     3358.715953
 rhoep=     -129.882511   utot=     -6533.597183   ehks=    -3304.763741
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 2.  RMS DQ=1.66e-3  last it=7.01e-4
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.701348      3.728857      3.728857      0.000493      3.728857
 site    1    7.298652      7.271143      7.271143      0.001119      7.271143
 AMIX: nmix=2 mmix=8  nelts=2178  beta=1  tm=5  rmsdel=1.66e-3
   tj: 0.70478  -0.00608
 unscreened rms difference:  smooth  0.000361   local  0.001119
   screened rms difference:  smooth  0.000493   local  0.001119   tot  0.001661

 iors  : write restart file (binary, mesh density) 

   it  3  of 12    ehf=      -0.329708   ehk=      -0.329241
 From last iter    ehf=      -0.329597   ehk=      -0.329511
 diffe(q)= -0.000111 (0.001661)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=1 pwmode=11 ehf=-.3297082 ehk=-.3292407

 --- BNDFP:  begin iteration 4 of 12 ---

 avg es pot at rmt= 0.582600  avg sphere pot= 0.651625  vconst=-0.582600

 smooth rhoves     10.173537   charge     3.709052
 smooth rhoeps =   -2.662395   rhomu =   -3.464452  avg vxc =   -0.832971 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000136

 Energy terms:             smooth           local           total
   rhoval*vef             -5.979897      -184.380137      -190.360033
   rhoval*ves            -46.686428      -116.551059      -163.237487
   psnuc*ves              67.033502    -12972.035039    -12905.001537
   utot                   10.173537     -6544.293049     -6534.119512
   rho*exc                -2.662395      -127.271295      -129.933690
   rho*vxc                -3.464452      -168.526220      -171.990672
   valence chg             3.709052         7.290948        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250   ndimh = 39
 -0.7179 -0.2594 -0.2550 -0.2550 -0.1941 -0.1941  1.5841  1.7300  1.8030
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250   ndimh = 37
 -0.4157 -0.2936 -0.2657 -0.2281 -0.1624 -0.1035  0.5904  1.1752  1.3019
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750   ndimh = 38
 -0.4779 -0.2830 -0.2549 -0.2387 -0.1818 -0.1481  0.6972  1.2832  1.5682
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250   ndimh = 36
 -0.4559 -0.3162 -0.2525 -0.2037 -0.1782 -0.1435  0.9064  1.0949  1.2380
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750   ndimh = 37
 -0.3434 -0.3053 -0.2920 -0.2035 -0.1689  0.0436  0.5463  0.7738  1.0760
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250   ndimh = 38
 -0.5007 -0.2664 -0.2488 -0.2488 -0.1752 -0.1752  0.6757  1.6076  1.6263

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.019328;  11.000000 electrons;  D(Ef):    4.267
         Sum occ. bands:   -2.8275552  incl. Bloechl correction:   -0.009151

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250   ndimh = 39
 -0.7179 -0.2594 -0.2550 -0.2550 -0.1941 -0.1941  1.5841  1.7300  1.8030
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250   ndimh = 37
 -0.4157 -0.2936 -0.2657 -0.2281 -0.1624 -0.1035  0.5904  1.1752  1.3019
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750   ndimh = 38
 -0.4779 -0.2830 -0.2549 -0.2387 -0.1818 -0.1481  0.6972  1.2832  1.5682
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250   ndimh = 36
 -0.4559 -0.3162 -0.2525 -0.2037 -0.1782 -0.1435  0.9064  1.0949  1.2380
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750   ndimh = 37
 -0.3434 -0.3053 -0.2920 -0.2035 -0.1689  0.0436  0.5463  0.7738  1.0760
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250   ndimh = 38
 -0.5007 -0.2664 -0.2488 -0.2488 -0.1752 -0.1752  0.6757  1.6076  1.6263

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.019328;  11.000000 electrons;  D(Ef):    4.267
         Sum occ. bands:   -2.8275552  incl. Bloechl correction:   -0.009151

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.131356    2.840729    7.290628

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.475702   -1.589719    4.500000    4.240220    4.500000    4.500000
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.409809   -1.222218    4.250000    4.197418    4.250000    4.250000
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.994711   -0.268941    3.841141    3.841504    3.147584    3.841504
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.021110   -0.731679    4.110000    4.108402    4.102416    4.110000
 4     1    0.004573   -0.629225    5.100000    5.081019    5.077979    5.100000

 Harris energy:
 sumev=       -2.827555  val*vef=    -190.360033   sumtv=     187.532478
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.933690     utot=   -6534.119512    ehar=   -3304.764085

 srhov:     -5.980584   -184.373268   -190.353852 sumev=   -2.827555   sumtv=  187.526296

 Kohn-Sham energy:
 sumtv=      187.526296  sumtc=      3171.756639   ekin=     3359.282936
 rhoep=     -129.933168   utot=     -6534.113852   ehks=    -3304.764085
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 3.  RMS DQ=1.44e-5  last it=1.66e-3
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.709052      3.709372      3.709372      0.000011      3.709372
 site    1    7.290948      7.290628      7.290628      0.000010      7.290628
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: nmix=2 mmix=8  nelts=2178  beta=1  tm=5  rmsdel=1.44e-5
   tj:-0.02068  -0.02803
 unscreened rms difference:  smooth  0.000011   local  0.000010
   screened rms difference:  smooth  0.000011   local  0.000010   tot  0.000014

 iors  : write restart file (binary, mesh density) 

   it  4  of 12    ehf=      -0.329585   ehk=      -0.329585
 From last iter    ehf=      -0.329708   ehk=      -0.329241
 diffe(q)=  0.000124 (0.000014)    tol= 0.000010 (0.000010)   more=T
i nk=8 bigbas=1 pwmode=11 ehf=-.3295846 ehk=-.3295846

 --- BNDFP:  begin iteration 5 of 12 ---

 avg es pot at rmt= 0.582617  avg sphere pot= 0.651617  vconst=-0.582617

 smooth rhoves     10.174208   charge     3.709215
 smooth rhoeps =   -2.662557   rhomu =   -3.464663  avg vxc =   -0.832977 

 locpot:

 site  1  z= 29.0  rmt= 2.31127  nr=393   a=0.025  nlml=25  rg=0.578  Vfloat=T
 sm core charge = 0.263001 (sphere) + 0.004646 (spillout) = 0.267647
 potential shift to crystal energy zero:    0.000135

 Energy terms:             smooth           local           total
   rhoval*vef             -5.980469      -184.378450      -190.358919
   rhoval*ves            -46.687197      -116.549295      -163.236492
   psnuc*ves              67.035613    -12972.035250    -12904.999637
   utot                   10.174208     -6544.292273     -6534.118064
   rho*exc                -2.662557      -127.271019      -129.933576
   rho*vxc                -3.464663      -168.525858      -171.990521
   valence chg             3.709215         7.290785        11.000000
   core charge            18.000000         0.000000        18.000000

 Charges:  valence    11.00000   cores    18.00000   nucleii   -29.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 Start first of two band passes ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250   ndimh = 39
 -0.7179 -0.2594 -0.2551 -0.2551 -0.1942 -0.1942  1.5841  1.7300  1.8030
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250   ndimh = 37
 -0.4157 -0.2936 -0.2657 -0.2281 -0.1625 -0.1035  0.5904  1.1752  1.3019
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750   ndimh = 38
 -0.4779 -0.2830 -0.2549 -0.2388 -0.1819 -0.1482  0.6972  1.2831  1.5682
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250   ndimh = 36
 -0.4559 -0.3162 -0.2525 -0.2038 -0.1782 -0.1436  0.9064  1.0949  1.2380
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750   ndimh = 37
 -0.3435 -0.3053 -0.2920 -0.2035 -0.1689  0.0436  0.5463  0.7738  1.0760
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250   ndimh = 38
 -0.5007 -0.2664 -0.2489 -0.2489 -0.1753 -0.1753  0.6757  1.6076  1.6263

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.019360;  11.000000 electrons;  D(Ef):    4.266
         Sum occ. bands:   -2.8280026  incl. Bloechl correction:   -0.009152

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 60, k=  0.06250  0.06250  0.06250   ndimh = 39
 -0.7179 -0.2594 -0.2551 -0.2551 -0.1942 -0.1942  1.5841  1.7300  1.8030
 bndfp:  kpt 11 of 60, k=  -0.18750  0.31250  0.56250   ndimh = 37
 -0.4157 -0.2936 -0.2657 -0.2281 -0.1625 -0.1035  0.5904  1.1752  1.3019
 bndfp:  kpt 21 of 60, k=  0.43750  -0.31250  0.18750   ndimh = 38
 -0.4779 -0.2830 -0.2549 -0.2388 -0.1819 -0.1482  0.6972  1.2831  1.5682
 bndfp:  kpt 31 of 60, k=  -0.06250  0.18750  -0.56250   ndimh = 36
 -0.4559 -0.3162 -0.2525 -0.2038 -0.1782 -0.1436  0.9064  1.0949  1.2380
 bndfp:  kpt 41 of 60, k=  -0.06250  0.43750  0.68750   ndimh = 37
 -0.3435 -0.3053 -0.2920 -0.2035 -0.1689  0.0436  0.5463  0.7738  1.0760
 bndfp:  kpt 51 of 60, k=  0.31250  0.31250  0.31250   ndimh = 38
 -0.5007 -0.2664 -0.2489 -0.2489 -0.1753 -0.1753  0.6757  1.6076  1.6263

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.019360;  11.000000 electrons;  D(Ef):    4.266
         Sum occ. bands:   -2.8280026  incl. Bloechl correction:   -0.009152

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1   10.131402    2.840631    7.290771

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.475693   -1.589775    4.500000    4.240213    4.500000    4.500000
 0     0    sc          0.000000    5.500000    5.500000    5.500000    5.500000
 1     0    0.409781   -1.222360    4.250000    4.197405    4.250000    4.250000
 1     0    sc          0.000000    5.500000    5.250000    5.250000    5.500000
 2     0    8.995070   -0.268977    3.841504    3.841517    3.147584    3.841517
 2     0    sc          0.000000    4.500000    4.147584    4.147584    4.500000
 3     1    0.021108   -0.731729    4.110000    4.108401    4.102416    4.110000
 4     1    0.004573   -0.629279    5.100000    5.081018    5.077979    5.100000

 Harris energy:
 sumev=       -2.828003  val*vef=    -190.358919   sumtv=     187.530916
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    3171.756639
 rhoeps=    -129.933576     utot=   -6534.118064    ehar=   -3304.764084

 srhov:     -5.980500   -184.377900   -190.358400 sumev=   -2.828003   sumtv=  187.530398

 Kohn-Sham energy:
 sumtv=      187.530398  sumtc=      3171.756639   ekin=     3359.287037
 rhoep=     -129.933532   utot=     -6534.117590   ehks=    -3304.764084
  
 mixing: mode=A  nmix=3  beta=1  elind=1.291
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=1.23e-6  last it=1.44e-5
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.709215      3.709228      3.709228      0.000008      3.709228
 site    1    7.290785      7.290771      7.290771      0.000001      7.290771
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=2178  beta=1  tm=5  rmsdel=1.23e-6
   tj:-0.08326
 unscreened rms difference:  smooth  0.000008   local  0.000001
   screened rms difference:  smooth  0.000008   local  0.000001   tot  0.000001

 iors  : write restart file (binary, mesh density) 

   it  5  of 12    ehf=      -0.329584   ehk=      -0.329584
 From last iter    ehf=      -0.329585   ehk=      -0.329585
 diffe(q)=  0.000000 (0.000001)    tol= 0.000010 (0.000010)   more=F
c nk=8 bigbas=1 pwmode=11 ehf=-.3295844 ehk=-.3295845
 Exit 0 LMF 
 CPU time:    7.694s   Wall clock    7.731s   12:04:12 02.10.2014        on phpdl1.ph.kcl.ac.uk
