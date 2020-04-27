rdcmd:  lmfa --no-iactiv c -vzbak=1
 ----------------------  START LMFA  ----------------------
 HEADER sc C atom

 LMFA:     nbas = 1  nspec = 1  verb 30,40,30
 pot:      spin-pol, XC:BH
 autogen:  none

                Plat                                  Qlat
   1.000000   1.000000   0.000000        0.500000   0.500000  -0.500000
   1.000000   0.000000   1.000000        0.500000  -0.500000   0.500000
   0.000000   1.000000   1.000000       -0.500000   0.500000   0.500000
   alat = 7.937005  Cell vol = 1000.000000

 LATTC:  as= 2.000   tol=1.00e-8   alat= 7.93701   awald= 0.200   lmax=6
         r1=  3.459   nkd= 79       q1=  2.571   nkg= 137

 SGROUP: 1 symmetry operations from 0 generators
 SYMLAT: Bravais system is cubic with 48 symmetry operations.
 SYMCRY: crystal invariant under 48 symmetry operations for tol=1e-5
 GROUPG: the following are sufficient to generate the space group:
         i*r3(-1,1,1) r4z
         i*r3(-1,1,1) r4z
 MKSYM:  found 48 space group operations ... includes inversion

 FREEAT: G tolerance for envelope functions : 1e-6

 Species C:  Z=6  Qc=2  R=3.000000  Q=0  mom=2
 mesh:   rmt=3.000000  rmax=19.671121  a=0.02  nr=369  nr(rmax)=463
  Pl=  2.5     2.5     3.5     4.5     spn 2   2.5     2.5     3.5     4.5    
  Ql=  1.0     2.0     0.0     0.0     spn 2   1.0     0.0     0.0     0.0    

  iter     qint         drho          vh0          rho0          vsum     beta
    1    6.000000   5.461E+02       30.0000    0.2984E+02      -12.0637   0.30
   50    6.000000   3.935E-05       29.2312    0.1279E+03      -59.7474   0.30


 sumev=-2.876387  etot=-74.994908  eref=-74.994900  diff= -0.000008

 Optimise free-atom basis for species C, rmt=3
 l  it    Rsm      Eh     stiffR   stiffE      Eval      Exact     Pnu    Ql
 0  41   1.425  -0.888      15.8     35.1   -1.07382  -1.07383    2.91   1.00
 1  31   1.429  -0.336      65.9    157.3   -0.46562  -0.46569    2.89   2.00
 eigenvalue sum:  exact  -2.00520    opt basis  -2.00507    error 0.00013

 Optimise free-atom basis for species C, spin 2
 l  it    Rsm      Eh     stiffR   stiffE      Eval      Exact     Pnu    Ql
 0  41   1.461  -0.748      17.8     52.8   -0.87118  -0.87119    2.91   1.00
 1  28   1.494  -0.209      80.9    335.8   -0.27747  -0.27759    2.87   0.00
 eigenvalue sum:  exact  -0.87119    opt basis  -0.87118    error 0.00001

 tailsm: fit tails to 6 smoothed hankels, rmt= 3.00000, rsm= 1.50000
 hnsmft (warning): rho(nr) > rhotol
    q(fit):     0.243570    rms diff:   0.000004
    fit: r>rmt  0.243570   r<rmt  1.753379   qtot  1.996950
    rho: r>rmt  0.243570   r<rmt  2.756430   qtot  3.000000

 tailsm: spin 2 ...
    q(fit):     0.054561    rms diff:   0.000002
    fit: r>rmt  0.054561   r<rmt  0.609738   qtot  0.664299
    rho: r>rmt  0.054561   r<rmt  0.945439   qtot  1.000000
 Exit 0 LMFA 
 CPU time:  0.114s   Wall clock 0.153s  at  10:20:45 10.08.2018  on  localhost.localdomai
rdcmd:  lmf  --no-iactiv c -vzbak=1
 ----------------------  START LMF  -----------------------
 HEADER sc C atom

 LMF:      nbas = 1  nspec = 1  verb 30,40,30
 special:  forces
 pot:      spin-pol, XC:BH
 float:    float P v6-style, v6-ebar
 autoread: none
 bz:       metal(2), tetra, invit 

                Plat                                  Qlat
   1.000000   1.000000   0.000000        0.500000   0.500000  -0.500000
   1.000000   0.000000   1.000000        0.500000  -0.500000   0.500000
   0.000000   1.000000   1.000000       -0.500000   0.500000   0.500000
   alat = 7.937005  Cell vol = 1000.000000

 LATTC:  as= 2.000   tol=1.00e-8   alat= 7.93701   awald= 0.200   lmax=6
         r1=  3.459   nkd= 79       q1=  2.571   nkg= 137

 SGROUP: 1 symmetry operations from 0 generators
 SYMLAT: Bravais system is cubic with 48 symmetry operations.
 SYMCRY: crystal invariant under 48 symmetry operations for tol=1e-5
 GROUPG: the following are sufficient to generate the space group:
         i*r3(-1,1,1) r4z
         i*r3(-1,1,1) r4z
 MKSYM:  found 48 space group operations ... includes inversion
 BZMESH: 8 irreducible QP from 64 ( 4 4 4 )  shift= F F F

 species data:  augmentation                           density
 spec       rmt   rsma lmxa kmxa      lmxl     rg   rsmv  kmxv foca   rfoca
 C        3.000  1.200    3    3         3  0.750  1.500    15    0   1.200

 GVLIST: gmax = 13.99 a.u. created 45911 vectors of 125000 (36%)
         mesh has 50 x 50 x 50 divisions; length 0.224, 0.224, 0.224
 SGVSYM: 1207 symmetry stars found for 45911 reciprocal lattice vectors
 

 Makidx:  hamiltonian dimensions Low, Int, High, Negl: 8 0 24 0
 suham :  16 augmentation channels, 16 local potential channels  Maximum lmxa=3

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol=1.0e-6
 spec      l    rsm    eh     gmax    last term   cutoff
  C        0    1.30  -0.70   5.718    1.05E-06    3143 
  C        1    1.10  -0.20   7.226    1.06E-06    6375 
  C        0    0.80  -1.50   9.292    1.08E-06   13539 
  C        1    0.80  -1.00  10.038    1.00E-06   16961 
 

 lmfp  : no rst file ... try to overlap atomic densities
 rdovfa: read and overlap free-atom densities (mesh density) ...
 rdovfa: expected C,       read C        with rmt=  3.0000  mesh   369  0.020

 Free atom and overlapped crystal site charges  spec 1  site 1:
   ib    true(FA)    smooth(FA)  true(OV)    smooth(OV)    local
    1    3.701869    2.363117    3.701844    2.363092    1.338752
 amom    1.810990    1.143641    1.810990    1.143641    0.667349
 Uniform density added to neutralize background, q=1.000000

 Smooth charge on mesh:            1.661248    moment    1.332651
 Sum of local charges:             1.338752    moments   0.667349
 Total valence charge:             3.000000    moment    2.000000
 Sum of core charges:              2.000000    moment    0.000000
 Sum of nuclear charges:          -6.000000
 Homogeneous background:           1.000000
 Deviation from neutrality:        0.000000

 --- BNDFP:  begin iteration 1 of 10 ---

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1    0.377655    1.338752     2.00     6.00

 Smooth charges: Qmesh = 1.661248  Qgauss = 1.338752  core-nuc = -4  tot = -1

 Average es pot at rmt = 0.006604  avg sphere pot = 0.019543   vconst = -0.006604
 Average es pot at MT boundaries after vconst shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159
 smooth rhoves      2.063389   charge     2.661248
 smvxcm (warning) mesh density negative at 201254 points:  rhomin=-5e-4
 smooth rhoeps =   -1.385392 (  -1.042886,  -0.342505)
         rhomu =   -1.806012 (  -1.431281,  -0.374731)
       avg vxc =   -0.082062 (  -0.095980,  -0.068143)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 potential shift to crystal energy zero:    0.000003

 potpus  spin 1 : pnu = 2.900000 2.850000 3.180000 4.120000
 l        enu         v           c          srdel        qpar        ppar
 0     -1.049913   -1.260952   -1.152144    0.197581      2.7872    0.898891
 1     -0.465402   -1.039466   -0.463120    0.173954     19.0465    0.926620
 2     -0.311939   -0.604280    1.529957    0.353092     17.1184    6.114600
 3     -0.096237   -0.537543    3.382042    0.400470     24.4393   10.634942

 potpus  spin 2 : pnu = 2.900000 2.850000 3.180000 4.120000
 l        enu         v           c          srdel        qpar        ppar
 0     -0.838165   -1.073432   -0.951603    0.208583      2.8002    0.932016
 1     -0.258802   -0.900692   -0.256240    0.184285     18.9762    0.977114
 2     -0.186358   -0.480549    1.677814    0.356097     17.0209    6.251001
 3      0.023022   -0.419226    3.516802    0.401774     24.3828   10.732448

 3rd gen potential parameters
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000    -9.233051    7.622124    0.101666   -0.583560
 1      3.000000    1.000000    -5.887832    7.139035    0.143256   -0.535856
 2      3.000000    1.000000     4.727244   27.072292    0.465409   -0.096157
 3      3.000000    1.000000     7.577135   37.162816    0.543350   -0.062205

 3rd gen potential parameters spin 2
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000    -9.233051    7.498210    0.106861   -0.559304
 1      3.000000    1.000000    -5.887832    7.161455    0.151762   -0.504954
 2      3.000000    1.000000     4.727244   27.364844    0.467062   -0.094578
 3      3.000000    1.000000     7.577135   37.315959    0.544001   -0.061811

 Local terms accumulated over augmentation sites:
                  true           sm,loc      difference
 rhoeps:        -9.464952      -1.363859      -8.101093
 rhomu:         -7.369404      -1.401629      -5.967775
 spin2:         -5.086143      -0.374960      -4.711183
 total:        -12.455547      -1.776589     -10.678958
 val*ves       -10.246156      -5.064346      -5.181810
 val*vef       -14.138038      -6.807368      -7.330670
 val chg:        3.588746       2.249995       1.338752
 core chg:       2.000000       2.000000       0.000000
 val mom:        1.810990

 Energy terms:             smooth           local            total
   rhoval*vef             -3.783717       -10.404877       -14.188595
   rhoval*ves             -5.058517        -5.181810       -10.240327
   psnuc*ves               9.185295      -278.835495      -269.650200
   utot                    2.063389      -142.008652      -139.945263
   rho*exc                -1.385392        -8.101093        -9.486484
   rho*vxc                -1.806012       -10.678958       -12.484969
   valence chg             1.661248         1.338752         3.000000
   valence mag             1.332651         0.667349         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000

 Incompatible or missing qp weights file ...

 Start first of two band passes ...
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.0318 -0.4143 -0.4143 -0.4143  0.2634  0.6583  0.6583  0.6583
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -0.8236 -0.2132 -0.2132 -0.2132  0.3259  0.7397  0.7397  0.7397

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.415970;   3.000000 electrons;  D(Ef):  794.026
         Sum occ. bands:   -2.2712367  incl. Bloechl correction:   -0.000057
         Mag. moment:       1.000000

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.0318 -0.4143 -0.4143 -0.4143  0.2634  0.6583  0.6583  0.6583
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -0.8236 -0.2132 -0.2132 -0.2132  0.3259  0.7397  0.7397  0.7397
 Est Ef = -0.416 < evl(3)=-0.414 ... using qval=3.0, revise to -0.4143

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.415970;   3.000000 electrons;  D(Ef):  794.026
         Sum occ. bands:   -2.2712367  incl. Bloechl correction:   -0.000057
         Mag. moment:       1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.831917    4.874573   -2.042656      0.924612    1.907971   -0.983358
       contr. to mm extrapolated for r>rmt:   0.058538 est. true mm = 0.983151

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.963385   -1.028882    2.900000    2.919100    2.500000    2.919100
 spn 2 0    0.953652   -0.821345    2.900000    2.914124    2.500000    2.914124
 1     1    0.914877   -0.410818    2.850000    2.900121    2.250000    2.850000
 spn 2 1    0.000000   -0.978123    2.850000    2.184933    2.250000    2.850000
 2     0    0.000002   -0.817016    3.180000    3.131388    3.147584    3.147584
 spn 2 0    0.000000   -1.162027    3.180000    3.107693    3.147584    3.147584
 3     0    0.000000   -0.783443    4.120000    4.095082    4.102416    4.102416
 spn 2 0    0.000000   -1.113939    4.120000    4.084684    4.102416    4.102416

 Harris energy:
 sumev=       -2.271237  val*vef=     -14.188595   sumtv=      11.917358
 sumec=      -39.640082  cor*vef=    -102.572087   ttcor=      62.932005
 rhoeps=      -9.486484     utot=    -139.945263    ehar=     -74.582386

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:    -11.238158     -0.363465    -11.601624 sumev=   -2.271237   sumtv=    9.330387

 Kohn-Sham energy:
 sumtv=        9.330387  sumtc=        62.933572   ekin=       72.263959
 rhoep=       -8.664921   utot=      -137.998385   ehks=      -74.399347
 mag. mom=     1.000000  (bands)        1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho:  sought 2 iter from file mixm; read 0.  RMS DQ=3.03e-2
 mixrho: (warning) scr. and lin-mixed densities had 81221 and 91739 negative points
 AMIX: nmix=0 mmix=8  nelts=93876  beta=0.5  tm=5  rmsdel=1.51e-2
 unscreened rms difference:  smooth  0.024271   local  0.078369
   screened rms difference:  smooth  0.024807   local  0.078369   tot  0.030257

 iors  : write restart file (binary, mesh density)

   it  1  of 10    ehf=       0.412514   ehk=       0.595553
h zbak=1 mmom=.9999997 ehf=.4125145 ehk=.5955525

 --- BNDFP:  begin iteration 2 of 10 ---

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.099284   -0.351952     2.00     6.00

 Smooth charges: Qmesh = 3.351952  Qgauss = -0.351952  core-nuc = -4  tot = -1

 Average es pot at rmt = -0.014891  avg sphere pot = 0.015494   vconst = 0.014891
 Average es pot at MT boundaries after vconst shift
 Site    ves
   1  0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159
 smooth rhoves      3.632487   charge     4.351952
 smvxcm (warning) mesh density negative at 184710 points:  rhomin=-5.16e-4
 smooth rhoeps =   -2.984056 (  -2.104096,  -0.879960)
         rhomu =   -3.895532 (  -2.901420,  -0.994112)
       avg vxc =   -0.104546 (  -0.116507,  -0.092585)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 potential shift to crystal energy zero:    0.000005

 potpus  spin 1 : pnu = 2.919100 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.154375   -1.362535   -1.264754    0.187142      2.7919    0.852732
 1     -0.585446   -1.085502   -0.583455    0.162470     19.0194    0.871960
 2     -0.626895   -0.626895    1.565741    0.365120     16.4470    6.752506
 3     -0.553736   -0.553736    3.494650    0.415121     23.4918   11.794202

 potpus  spin 2 : pnu = 2.914124 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -0.990437   -1.214657   -1.106156    0.196723      2.8036    0.881733
 1     -0.421127   -0.983119   -0.418905    0.171635     19.1527    0.911995
 2     -0.544160   -0.544160    1.673819    0.368116     16.3674    6.895685
 3     -0.477273   -0.477273    3.589992    0.416525     23.4424   11.902612

 3rd gen potential parameters
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -11.548569    7.576326    0.084501   -0.618779
 1      3.000000    1.000000    -5.887832    7.147664    0.133798   -0.573354
 2      3.000000    1.000000     6.000000   29.266579    0.490317   -0.087656
 3      3.000000    1.000000     9.000000   39.974194    0.568743   -0.056763

 3rd gen potential parameters spin 2
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.848796    7.465993    0.091788   -0.594849
 1      3.000000    1.000000    -5.887832    7.105616    0.141349   -0.544479
 2      3.000000    1.000000     6.000000   29.557622    0.491946   -0.086286
 3      3.000000    1.000000     9.000000   40.136088    0.569469   -0.056396

 Local terms accumulated over augmentation sites:
                  true           sm,loc      difference
 rhoeps:        -9.150184      -2.947224      -6.202960
 rhomu:         -6.872281      -2.862791      -4.009490
 spin2:         -5.168896      -0.984802      -4.184094
 total:        -12.041177      -3.847593      -8.193584
 val*ves       -10.063967      -4.916527      -5.147440
 val*vef       -13.548460      -8.719110      -4.829350
 val chg:        3.365684       3.717636      -0.351952
 core chg:       2.000000       2.000000       0.000000
 val mom:        1.367801

 Energy terms:             smooth           local            total
   rhoval*vef             -9.871509        -3.759728       -13.631237
   rhoval*ves             -4.891464        -5.147440       -10.038904
   psnuc*ves              12.156438      -280.500417      -268.343980
   utot                    3.632487      -142.823928      -139.191442
   rho*exc                -2.984056        -6.202960        -9.187016
   rho*vxc                -3.895532        -8.193584       -12.089116
   valence chg             3.351952        -0.351952         3.000000
   valence mag             1.658004        -0.158005         1.500000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.41597,  3 occ states
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.1539 -0.5351 -0.5351 -0.5351  0.2375  0.6146  0.6146  0.6146
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -0.9884 -0.3731 -0.3731 -0.3731  0.2850  0.6760  0.6760  0.6760

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.535712;   3.000000 electrons;  D(Ef): 2485.897
         Sum occ. bands:   -2.6778599  incl. Bloechl correction:   -0.000037
         Mag. moment:       1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.871442    5.682276   -2.810834      0.945164    1.927074   -0.981910
       contr. to mm extrapolated for r>rmt:   0.044430 est. true mm = 0.989594

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.969793   -1.151455    2.919100    2.921811    2.500000    2.921811
 spn 2 0    0.963139   -0.985831    2.914124    2.918176    2.500000    2.918176
 1     1    0.938509   -0.529486    2.850000    2.907426    2.250000    2.850000
 spn 2 1    0.000000   -0.378791    2.850000    2.891127    2.250000    2.850000
 2     0    0.000001   -0.857342    3.147584    3.130184    3.147584    3.147584
 spn 2 0    0.000000   -1.251076    3.147584    3.106605    3.147584    3.147584
 3     0    0.000000   -0.830711    4.102416    4.094254    4.102416    4.102416
 spn 2 0    0.000000   -1.189196    4.102416    4.084334    4.102416    4.102416

 Harris energy:
 sumev=       -2.677860  val*vef=     -13.631237   sumtv=      10.953377
 sumec=      -40.271398  cor*vef=    -103.204147   ttcor=      62.932750
 rhoeps=      -9.187016     utot=    -139.191442    ehar=     -74.492331

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:    -18.007343      5.246855    -12.760488 sumev=   -2.677860   sumtv=   10.082628

 Kohn-Sham energy:
 sumtv=       10.082628  sumtc=        63.004518   ekin=       73.087146
 rhoep=       -8.772014   utot=      -138.748759   ehks=      -74.433628
 mag. mom=     1.000000  (bands)        1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho:  sought 2 iter from file mixm; read 1.  RMS DQ=2.30e-2  last it=3.03e-2
 mixrho: (warning) scr. and lin-mixed densities had 79169 and 86729 negative points
 AMIX: nmix=1 mmix=8  nelts=93876  beta=0.5  tm=5  rmsdel=1.15e-2
   tj:-1.96333
 unscreened rms difference:  smooth  0.021247   local  0.055285
   screened rms difference:  smooth  0.021520   local  0.055285   tot  0.022984

 iors  : write restart file (binary, mesh density)

   it  2  of 10    ehf=       0.502569   ehk=       0.561272
 From last iter    ehf=       0.412514   ehk=       0.595553
 diffe(q)=  0.090055 (0.022984)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999997 ehf=.5025692 ehk=.5612723

 --- BNDFP:  begin iteration 3 of 10 ---

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.127023   -3.995192     2.00     6.00

 Smooth charges: Qmesh = 6.995192  Qgauss = -3.995192  core-nuc = -4  tot = -1

 Average es pot at rmt = -0.082405  avg sphere pot = 0.008938   vconst = 0.082405
 Average es pot at MT boundaries after vconst shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159
 smooth rhoves      6.802845   charge     7.995192
 smvxcm (warning) mesh density negative at 135582 points:  rhomin=-4.03e-4
 smooth rhoeps =   -7.779798 (  -5.063245,  -2.716553)
         rhomu =  -10.184710 (  -6.974796,  -3.209914)
       avg vxc =   -0.148861 (  -0.154299,  -0.143424)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 potential shift to crystal energy zero:    0.000008

 potpus  spin 1 : pnu = 2.921811 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.294567   -1.473240   -1.391472    0.172024      2.7631    0.813612
 1     -0.722270   -1.130068   -0.720610    0.148362     18.6020    0.813829
 2     -0.626610   -0.626610    1.521581    0.359852     16.5888    6.513098
 3     -0.542264   -0.542264    3.472399    0.412605     23.5811   11.603667

 potpus  spin 2 : pnu = 2.918176 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.192842   -1.383550   -1.294148    0.179481      2.7753    0.834878
 1     -0.619267   -1.072989   -0.617443    0.155510     18.8371    0.842819
 2     -0.594521   -0.594521    1.575798    0.362482     16.5174    6.632794
 3     -0.515737   -0.515737    3.515385    0.413830     23.5377   11.695736

 3rd gen potential parameters
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -11.966470    7.862543    0.076944   -0.655419
 1      3.000000    1.000000    -5.887832    7.285035    0.122169   -0.621380
 2      3.000000    1.000000     6.000000   28.765907    0.487409   -0.090118
 3      3.000000    1.000000     9.000000   39.685556    0.567445   -0.057428

 3rd gen potential parameters spin 2
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -11.412278    7.739429    0.082176   -0.635390
 1      3.000000    1.000000    -5.887832    7.206547    0.128061   -0.596344
 2      3.000000    1.000000     6.000000   29.015425    0.488857   -0.088877
 3      3.000000    1.000000     9.000000   39.824924    0.568084   -0.057104

 Local terms accumulated over augmentation sites:
                  true           sm,loc      difference
 rhoeps:        -8.754615      -7.727681      -1.026933
 rhomu:         -6.242598      -6.938589       0.695991
 spin2:         -5.280402      -3.178690      -2.101713
 total:        -11.523000     -10.117278      -1.405722
 val*ves       -10.304324      -3.027580      -7.276744
 val*vef       -13.262575     -13.082823      -0.179752
 val chg:        2.925782       6.920975      -3.995192
 core chg:       2.000000       2.000000       0.000000
 val mom:        0.741594

 Energy terms:             smooth           local            total
   rhoval*vef            -29.644767        16.302001       -13.342765
   rhoval*ves             -2.895898        -7.276744       -10.172642
   psnuc*ves              16.501587      -283.911064      -267.409477
   utot                    6.802845      -145.593904      -138.791060
   rho*exc                -7.779798        -1.026933        -8.806731
   rho*vxc               -10.184710        -1.405722       -11.590432
   valence chg             6.995192        -3.995192         3.000000
   valence mag             2.137923        -1.378756         0.759167
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.535712,  3 occ states
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.2926 -0.6731 -0.6731 -0.6731  0.2427  0.5962  0.5962  0.5962
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.1896 -0.5700 -0.5700 -0.5700  0.2549  0.6144  0.6144  0.6144

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.673233;   3.000000 electrons;  D(Ef):10546.263
         Sum occ. bands:   -3.1553637  incl. Bloechl correction:   -0.000014
         Mag. moment:       1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.908869    6.671393   -3.762524      0.962686    1.869075   -0.906389
       contr. to mm extrapolated for r>rmt:   0.032150 est. true mm = 0.994836

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.976938   -1.290826    2.921811    2.925706    2.500000    2.925706
 spn 2 0    0.973091   -1.187690    2.918176    2.923291    2.500000    2.923291
 1     1    0.958839   -0.668048    2.850000    2.914776    2.250000    2.850000
 spn 2 1    0.000000   -0.557786    2.850000    2.916561    2.250000    2.850000
 2     0    0.000000   -0.889352    3.147584    3.127961    3.147584    3.147584
 spn 2 0    0.000000   -1.340933    3.147584    3.104838    3.147584    3.147584
 3     0    0.000000   -0.865635    4.102416    4.093010    4.102416    4.102416
 spn 2 0    0.000000   -0.865635    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -3.155364  val*vef=     -13.342765   sumtv=      10.187402
 sumec=      -40.882529  cor*vef=    -103.851130   ttcor=      62.968602
 rhoeps=      -8.806731     utot=    -138.791060    ehar=     -74.441788

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:    -29.861637     15.764317    -14.097319 sumev=   -3.155364   sumtv=   10.941956

 Kohn-Sham energy:
 sumtv=       10.941956  sumtc=        63.016484   ekin=       73.958440
 rhoep=       -8.888051   utot=      -139.512152   ehks=      -74.441763
 mag. mom=     1.000000  (bands)        1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho:  sought 2 iter from file mixm; read 2.  RMS DQ=4.63e-3  last it=2.30e-2
 mixrho: (warning) scr. and lin-mixed densities had 48589 and 65031 negative points
 AMIX: nmix=2 mmix=8  nelts=93876  beta=0.5  tm=5  rmsdel=2.32e-3
   tj:-1.20888   0.96521
 unscreened rms difference:  smooth  0.002873   local  0.011987
   screened rms difference:  smooth  0.002846   local  0.011987   tot  0.004633

 iors  : write restart file (binary, mesh density)

   it  3  of 10    ehf=       0.553112   ehk=       0.553137
 From last iter    ehf=       0.502569   ehk=       0.561272
 diffe(q)=  0.050543 (0.004633)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999997 ehf=.5531124 ehk=.553137

 --- BNDFP:  begin iteration 4 of 10 ---

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.917376   -3.252012     2.00     6.00

 Smooth charges: Qmesh = 6.252012  Qgauss = -3.252012  core-nuc = -4  tot = -1

 Average es pot at rmt = -0.086311  avg sphere pot = 0.011333   vconst = 0.086311
 Average es pot at MT boundaries after vconst shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159
 smooth rhoves      5.453840   charge     7.252012
 smvxcm (warning) mesh density negative at 152232 points:  rhomin=-2.87e-4
 smooth rhoeps =   -6.828572 (  -4.381009,  -2.447563)
         rhomu =   -8.937146 (  -6.018009,  -2.919137)
       avg vxc =   -0.124706 (  -0.133566,  -0.115846)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 potential shift to crystal energy zero:    0.000008

 potpus  spin 1 : pnu = 2.925706 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.229446   -1.422832   -1.335701    0.177205      2.7747    0.826404
 1     -0.661195   -1.097055   -0.659425    0.153185     18.6498    0.836655
 2     -0.608342   -0.608342    1.547429    0.360799     16.5601    6.562443
 3     -0.525637   -0.525637    3.492625    0.412895     23.5690   11.631159

 potpus  spin 2 : pnu = 2.923291 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.098490   -1.307875   -1.211727    0.185661      2.7893    0.849461
 1     -0.531053   -1.018659   -0.529094    0.161157     18.8499    0.870237
 2     -0.554081   -0.554081    1.625160    0.363586     16.4847    6.690770
 3     -0.477261   -0.477261    3.558450    0.414193     23.5232   11.729165

 3rd gen potential parameters
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.619222    7.745276    0.076344   -0.643200
 1      3.000000    1.000000    -5.887832    7.268839    0.126141   -0.602554
 2      3.000000    1.000000     6.000000   28.865579    0.487845   -0.089645
 3      3.000000    1.000000     9.000000   39.724058    0.567555   -0.057345

 3rd gen potential parameters spin 2
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.206866    7.601606    0.081133   -0.622226
 1      3.000000    1.000000    -5.887832    7.202376    0.132712   -0.575628
 2      3.000000    1.000000     6.000000   29.131266    0.489377   -0.088338
 3      3.000000    1.000000     9.000000   39.871959    0.568232   -0.057002

 Local terms accumulated over augmentation sites:
                  true           sm,loc      difference
 rhoeps:        -8.891093      -6.800699      -2.090394
 rhomu:         -6.482897      -5.993939      -0.488958
 spin2:         -5.219907      -2.907123      -2.312784
 total:        -11.702804      -8.901062      -2.801742
 val*ves       -10.556639      -4.031932      -6.524707
 val*vef       -13.686088     -12.876651      -0.809437
 val chg:        2.992540       6.244551      -3.252012
 core chg:       2.000000       2.000000       0.000000
 val mom:        0.981948

 Energy terms:             smooth           local            total
   rhoval*vef            -24.972128        11.233727       -13.738400
   rhoval*ves             -3.905507        -6.524707       -10.430213
   psnuc*ves              14.813186      -283.006474      -268.193288
   utot                    5.453840      -144.765590      -139.311751
   rho*exc                -6.828572        -2.090394        -8.918966
   rho*vxc                -8.937146        -2.801742       -11.738889
   valence chg             6.252012        -3.252012         3.000000
   valence mag             1.915108        -0.884481         1.030627
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.673233,  3 occ states
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.2322 -0.6104 -0.6104 -0.6104  0.2726  0.6332  0.6332  0.6332
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.1019 -0.4808 -0.4808 -0.4808  0.3001  0.6700  0.6700  0.6700
 Est Ef = -0.673 < evl(3)=-0.610 ... using qval=3.0, revise to -0.6104

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.610642;   3.000000 electrons;  D(Ef): 7644.363
         Sum occ. bands:   -2.9446235  incl. Bloechl correction:   -0.000019
         Mag. moment:       1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.899777    6.589974   -3.690197      0.958795    1.943854   -0.985059
       contr. to mm extrapolated for r>rmt:   0.035172 est. true mm = 0.993967

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.975148   -1.230239    2.925706    2.924932    2.500000    2.924932
 spn 2 0    0.970491   -1.099752    2.923291    2.922126    2.500000    2.922126
 1     1    0.954138   -0.605017    2.850000    2.913328    2.250000    2.850000
 spn 2 1    0.000000   -0.472143    2.850000    2.910682    2.250000    2.850000
 2     0    0.000001   -0.862279    3.147584    3.128540    3.147584    3.147584
 spn 2 0    0.000000   -1.288081    3.147584    3.105345    3.147584    3.147584
 3     0    0.000000   -0.825578    4.102416    4.093620    4.102416    4.102416
 spn 2 0    0.000000   -0.825578    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -2.944623  val*vef=     -13.738400   sumtv=      10.793777
 sumec=      -40.539168  cor*vef=    -103.531717   ttcor=      62.992550
 rhoeps=      -8.918966     utot=    -139.311751    ehar=     -74.444390

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:    -26.900030     13.333915    -13.566115 sumev=   -2.944623   sumtv=   10.621491

 Kohn-Sham energy:
 sumtv=       10.621491  sumtc=        62.964678   ekin=       73.586169
 rhoep=       -8.845621   utot=      -139.183322   ehks=      -74.442774
 mag. mom=     1.000000  (bands)        1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=4.34e-3  last it=4.63e-3
 mixrho: (warning) scr. and lin-mixed densities had 69807 and 73869 negative points
 AMIX: nmix=2 mmix=8  nelts=93876  beta=0.5  tm=5  rmsdel=2.17e-3
   tj:-0.15164  -0.27092
 unscreened rms difference:  smooth  0.004238   local  0.009884
   screened rms difference:  smooth  0.004283   local  0.009884   tot  0.004337

 iors  : write restart file (binary, mesh density)

   it  4  of 10    ehf=       0.550510   ehk=       0.552126
 From last iter    ehf=       0.553112   ehk=       0.553137
 diffe(q)= -0.002603 (0.004337)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999997 ehf=.5505095 ehk=.5521264

 --- BNDFP:  begin iteration 5 of 10 ---

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.106158   -3.921228     2.00     6.00

 Smooth charges: Qmesh = 6.921228  Qgauss = -3.921228  core-nuc = -4  tot = -1

 Average es pot at rmt = -0.107910  avg sphere pot = 0.010635   vconst = 0.107910
 Average es pot at MT boundaries after vconst shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159
 smooth rhoves      5.801948   charge     7.921228
 smvxcm (warning) mesh density negative at 90762 points:  rhomin=-1.45e-4
 smooth rhoeps =   -7.872780 (  -4.983423,  -2.889357)
         rhomu =  -10.308534 (  -6.838400,  -3.470134)
       avg vxc =   -0.146223 (  -0.155430,  -0.137016)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 potential shift to crystal energy zero:    0.000008

 potpus  spin 1 : pnu = 2.924932 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.242676   -1.430227   -1.345656    0.174820      2.7672    0.822678
 1     -0.671943   -1.093807   -0.670218    0.151226     18.5222    0.830409
 2     -0.596460   -0.596460    1.548332    0.359526     16.5927    6.510115
 3     -0.510789   -0.510789    3.497755    0.412184     23.5932   11.581841

 potpus  spin 2 : pnu = 2.922126 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.115525   -1.317848   -1.224707    0.183009      2.7809    0.845555
 1     -0.544768   -1.016500   -0.542861    0.159025     18.7289    0.863189
 2     -0.544496   -0.544496    1.623388    0.362274     16.5178    6.635746
 3     -0.464622   -0.464622    3.560967    0.413454     23.5482   11.677444

 3rd gen potential parameters
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.484089    7.821073    0.076069   -0.647410
 1      3.000000    1.000000    -5.887832    7.312357    0.124524   -0.608364
 2      3.000000    1.000000     6.000000   28.752549    0.487080   -0.090232
 3      3.000000    1.000000     9.000000   39.646751    0.567158   -0.057530

 3rd gen potential parameters spin 2
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.016924    7.683206    0.081033   -0.626422
 1      3.000000    1.000000    -5.887832    7.242317    0.130953   -0.581585
 2      3.000000    1.000000     6.000000   29.013755    0.488590   -0.088932
 3      3.000000    1.000000     9.000000   39.791276    0.567820   -0.057193

 Local terms accumulated over augmentation sites:
                  true           sm,loc      difference
 rhoeps:        -8.845845      -7.851801      -0.994044
 rhomu:         -6.430374      -6.820495       0.390122
 spin2:         -5.213625      -3.460874      -1.752751
 total:        -11.643999     -10.281369      -1.362630
 val*ves       -10.667574      -3.815613      -6.851961
 val*vef       -13.736276     -14.038636       0.302361
 val chg:        2.906106       6.827334      -3.921228
 core chg:       2.000000       2.000000       0.000000
 val mom:        0.937879

 Energy terms:             smooth           local            total
   rhoval*vef            -29.010444        15.251158       -13.759286
   rhoval*ves             -3.645202        -6.851961       -10.497163
   psnuc*ves              15.249098      -283.416404      -268.167306
   utot                    5.801948      -145.134183      -139.332235
   rho*exc                -7.872780        -0.994044        -8.866824
   rho*vxc               -10.308534        -1.362630       -11.671164
   valence chg             6.921228        -3.921228         3.000000
   valence mag             1.974408        -1.002093         0.972314
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.610642,  3 occ states
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.2437 -0.6210 -0.6210 -0.6210  0.2782  0.6405  0.6405  0.6405
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.1166 -0.4940 -0.4940 -0.4940  0.3035  0.6736  0.6736  0.6736

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.621226;   3.000000 electrons;  D(Ef): 6119.471
         Sum occ. bands:   -2.9814368  incl. Bloechl correction:   -0.000022
         Mag. moment:       1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.904515    6.792221   -3.887706      0.960769    1.947004   -0.986234
       contr. to mm extrapolated for r>rmt:   0.033560 est. true mm = 0.994330

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.976241   -1.241951    2.924932    2.925655    2.500000    2.925655
 spn 2 0    0.971872   -1.114649    2.922126    2.922954    2.500000    2.922954
 1     1    0.956401   -0.616215    2.850000    2.914192    2.250000    2.850000
 spn 2 1    0.000000   -0.472374    2.850000    2.923112    2.250000    2.850000
 2     0    0.000001   -0.853433    3.147584    3.128301    3.147584    3.147584
 spn 2 0    0.000000   -1.286234    3.147584    3.104982    3.147584    3.147584
 3     0    0.000000   -0.823507    4.102416    4.093277    4.102416    4.102416
 spn 2 0    0.000000   -0.823507    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -2.981437  val*vef=     -13.759286   sumtv=      10.777849
 sumec=      -40.569509  cor*vef=    -103.548122   ttcor=      62.978613
 rhoeps=      -8.866824     utot=    -139.332235    ehar=     -74.442597

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:    -28.881536     15.183988    -13.697549 sumev=   -2.981437   sumtv=   10.716112

 Kohn-Sham energy:
 sumtv=       10.716112  sumtc=        62.956427   ekin=       73.672539
 rhoep=       -8.858588   utot=      -139.256866   ehks=      -74.442915
 mag. mom=     1.000000  (bands)        1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=2.88e-4  last it=4.34e-3
 mixrho: (warning) scr. and lin-mixed densities had 31085 and 39893 negative points
 AMIX: nmix=2 mmix=8  nelts=93876  beta=0.5  tm=5  rmsdel=1.44e-4
   tj: 0.05061  -0.01812
 unscreened rms difference:  smooth  0.000245   local  0.000816
   screened rms difference:  smooth  0.000231   local  0.000816   tot  0.000288

 iors  : write restart file (binary, mesh density)

   it  5  of 10    ehf=       0.552303   ehk=       0.551985
 From last iter    ehf=       0.550510   ehk=       0.552126
 diffe(q)=  0.001794 (0.000288)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999997 ehf=.5523031 ehk=.5519853

 --- BNDFP:  begin iteration 6 of 10 ---

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.095373   -3.882997     2.00     6.00

 Smooth charges: Qmesh = 6.882997  Qgauss = -3.882997  core-nuc = -4  tot = -1

 Average es pot at rmt = -0.108534  avg sphere pot = 0.010746   vconst = 0.108534
 Average es pot at MT boundaries after vconst shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159
 smooth rhoves      5.758897   charge     7.882997
 smvxcm (warning) mesh density negative at 88866 points:  rhomin=-1.05e-4
 smooth rhoeps =   -7.816557 (  -4.952415,  -2.864143)
         rhomu =  -10.234801 (  -6.796465,  -3.438336)
       avg vxc =   -0.142332 (  -0.152894,  -0.131770)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 potential shift to crystal energy zero:    0.000008

 potpus  spin 1 : pnu = 2.925655 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.241679   -1.430135   -1.345456    0.174908      2.7679    0.822475
 1     -0.671675   -1.093663   -0.669949    0.151258     18.5197    0.830611
 2     -0.596330   -0.596330    1.548222    0.359502     16.5930    6.509703
 3     -0.510568   -0.510568    3.497579    0.412158     23.5939   11.580567

 potpus  spin 2 : pnu = 2.922954 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.112712   -1.316051   -1.222815    0.183077      2.7817    0.845205
 1     -0.542860   -1.014439   -0.540953    0.159021     18.7239    0.863265
 2     -0.542298   -0.542298    1.625110    0.362224     16.5187    6.634372
 3     -0.462264   -0.462264    3.562675    0.413410     23.5494   11.675066

 3rd gen potential parameters
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.610160    7.813841    0.075614   -0.647515
 1      3.000000    1.000000    -5.887832    7.313206    0.124551   -0.608194
 2      3.000000    1.000000     6.000000   28.751477    0.487057   -0.090240
 3      3.000000    1.000000     9.000000   39.644500    0.567139   -0.057536

 3rd gen potential parameters spin 2
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.151328    7.675524    0.080487   -0.626639
 1      3.000000    1.000000    -5.887832    7.243987    0.130949   -0.581528
 2      3.000000    1.000000     6.000000   29.010591    0.488549   -0.088952
 3      3.000000    1.000000     9.000000   39.787249    0.567789   -0.057204

 Local terms accumulated over augmentation sites:
                  true           sm,loc      difference
 rhoeps:        -8.847509      -7.797699      -1.049809
 rhomu:         -6.438819      -6.779708       0.340890
 spin2:         -5.207342      -3.430672      -1.776671
 total:        -11.646161     -10.210380      -1.435781
 val*ves       -10.640678      -3.846315      -6.794363
 val*vef       -13.712988     -13.998561       0.285573
 val chg:        2.910209       6.793206      -3.882997
 core chg:       2.000000       2.000000       0.000000
 val mom:        0.952151

 Energy terms:             smooth           local            total
   rhoval*vef            -28.767155        15.034176       -13.732979
   rhoval*ves             -3.675216        -6.794363       -10.469580
   psnuc*ves              15.193010      -283.320017      -268.127007
   utot                    5.758897      -145.057190      -139.298294
   rho*exc                -7.816557        -1.049809        -8.866367
   rho*vxc               -10.234801        -1.435781       -11.670582
   valence chg             6.882997        -3.882997         3.000000
   valence mag             1.978032        -0.988468         0.989564
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.621226,  3 occ states
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.2434 -0.6205 -0.6205 -0.6205  0.2823  0.6450  0.6450  0.6450
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.1145 -0.4918 -0.4918 -0.4918  0.3107  0.6811  0.6811  0.6811
 Est Ef = -0.621 < evl(3)=-0.621 ... using qval=3.0, revise to -0.6205

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.620784;   3.000000 electrons;  D(Ef): 6121.885
         Sum occ. bands:   -2.9785783  incl. Bloechl correction:   -0.000021
         Mag. moment:       1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.904926    6.850282   -3.945356      0.960912    1.950682   -0.989770
       contr. to mm extrapolated for r>rmt:   0.033476 est. true mm = 0.994388

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.976308   -1.241586    2.925655    2.925747    2.500000    2.925747
 spn 2 0    0.972007   -1.112554    2.922954    2.923103    2.500000    2.923103
 1     1    0.956611   -0.615745    2.850000    2.914357    2.250000    2.850000
 spn 2 1    0.000000   -0.468820    2.850000    2.924410    2.250000    2.850000
 2     0    0.000001   -0.854134    3.147584    3.128248    3.147584    3.147584
 spn 2 0    0.000000   -1.285791    3.147584    3.104913    3.147584    3.147584
 3     0    0.000000   -0.828446    4.102416    4.093143    4.102416    4.102416
 spn 2 0    0.000000   -0.828446    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -2.978578  val*vef=     -13.732979   sumtv=      10.754401
 sumec=      -40.572353  cor*vef=    -103.539874   ttcor=      62.967520
 rhoeps=      -8.866367     utot=    -139.298294    ehar=     -74.442739

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:    -29.085124     15.388176    -13.696948 sumev=   -2.978578   sumtv=   10.718370

 Kohn-Sham energy:
 sumtv=       10.718370  sumtc=        62.960899   ekin=       73.679269
 rhoep=       -8.859284   utot=      -139.262907   ehks=      -74.442923
 mag. mom=     1.000000  (bands)        1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=6.76e-4  last it=2.88e-4
 mixrho: (warning) scr. and lin-mixed densities had 31367 and 37997 negative points
 AMIX: nmix=2 mmix=8  nelts=93876  beta=0.5  tm=5  rmsdel=3.38e-4
   tj: 0.32769  -0.09984
 unscreened rms difference:  smooth  0.000695   local  0.001446
   screened rms difference:  smooth  0.000698   local  0.001446   tot  0.000676

 iors  : write restart file (binary, mesh density)

   it  6  of 10    ehf=       0.552161   ehk=       0.551977
 From last iter    ehf=       0.552303   ehk=       0.551985
 diffe(q)= -0.000142 (0.000676)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999997 ehf=.5521607 ehk=.5519772

 --- BNDFP:  begin iteration 7 of 10 ---

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.115751   -3.955233     2.00     6.00

 Smooth charges: Qmesh = 6.955233  Qgauss = -3.955233  core-nuc = -4  tot = -1

 Average es pot at rmt = -0.110951  avg sphere pot = 0.010697   vconst = 0.110951
 Average es pot at MT boundaries after vconst shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159
 smooth rhoves      5.787227   charge     7.955233
 smvxcm (warning) mesh density negative at 57548 points:  rhomin=-6.89e-5
 smooth rhoeps =   -7.935136 (  -5.019902,  -2.915234)
         rhomu =  -10.390564 (  -6.888104,  -3.502460)
       avg vxc =   -0.153395 (  -0.164219,  -0.142570)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 potential shift to crystal energy zero:    0.000008

 potpus  spin 1 : pnu = 2.925747 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.242709   -1.430791   -1.346346    0.174688      2.7672    0.822071
 1     -0.672628   -1.093225   -0.670907    0.151068     18.5051    0.830051
 2     -0.594962   -0.594962    1.548320    0.359357     16.5966    6.504024
 3     -0.508844   -0.508844    3.498094    0.412070     23.5968   11.574800

 potpus  spin 2 : pnu = 2.923103 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.113430   -1.316368   -1.223414    0.182821      2.7811    0.844696
 1     -0.543535   -1.013392   -0.541634    0.158796     18.7086    0.862567
 2     -0.540371   -0.540371    1.625607    0.362061     16.5227    6.627921
 3     -0.459937   -0.459937    3.563622    0.413311     23.5526   11.668487

 3rd gen potential parameters
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.626297    7.820185    0.075480   -0.647952
 1      3.000000    1.000000    -5.887832    7.318239    0.124394   -0.608729
 2      3.000000    1.000000     6.000000   28.739104    0.486966   -0.090306
 3      3.000000    1.000000     9.000000   39.635332    0.567088   -0.057559

 3rd gen potential parameters spin 2
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.175721    7.681878    0.080297   -0.627143
 1      3.000000    1.000000    -5.887832    7.249094    0.130763   -0.582129
 2      3.000000    1.000000     6.000000   28.996684    0.488446   -0.089024
 3      3.000000    1.000000     9.000000   39.776826    0.567731   -0.057229

 Local terms accumulated over augmentation sites:
                  true           sm,loc      difference
 rhoeps:        -8.843779      -7.917442      -0.926338
 rhomu:         -6.436644      -6.872167       0.435523
 spin2:         -5.204672      -3.495481      -1.709191
 total:        -11.641315     -10.367648      -1.273667
 val*ves       -10.641131      -3.827881      -6.813249
 val*vef       -13.708795     -14.137254       0.428459
 val chg:        2.903235       6.858468      -3.955233
 core chg:       2.000000       2.000000       0.000000
 val mom:        0.952788

 Energy terms:             smooth           local            total
   rhoval*vef            -29.208791        15.484121       -13.724669
   rhoval*ves             -3.651613        -6.813249       -10.464863
   psnuc*ves              15.226068      -283.339115      -268.113046
   utot                    5.787227      -145.076182      -139.288955
   rho*exc                -7.935136        -0.926338        -8.861474
   rho*vxc               -10.390564        -1.273667       -11.664231
   valence chg             6.955233        -3.955233         3.000000
   valence mag             1.986104        -0.996199         0.989906
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.620784,  3 occ states
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.2444 -0.6214 -0.6214 -0.6214  0.2807  0.6464  0.6464  0.6464
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.1152 -0.4923 -0.4923 -0.4923  0.3097  0.6832  0.6832  0.6832

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.621688;   3.000000 electrons;  D(Ef): 5260.989
         Sum occ. bands:   -2.9812075  incl. Bloechl correction:   -0.000022
         Mag. moment:       1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.905328    6.883156   -3.977828      0.960960    1.946549   -0.985589
       contr. to mm extrapolated for r>rmt:   0.033389 est. true mm = 0.994349

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.976434   -1.242610    2.925747    2.925844    2.500000    2.925844
 spn 2 0    0.972184   -1.113296    2.923103    2.923229    2.500000    2.923229
 1     1    0.956710   -0.616819    2.850000    2.914369    2.250000    2.850000
 spn 2 1    0.000000   -0.464685    2.850000    2.928285    2.250000    2.850000
 2     0    0.000001   -0.852912    3.147584    3.128233    3.147584    3.147584
 spn 2 0    0.000000   -1.285278    3.147584    3.104851    3.147584    3.147584
 3     0    0.000000   -0.827915    4.102416    4.093110    4.102416    4.102416
 spn 2 0    0.000000   -0.827915    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -2.981208  val*vef=     -13.724669   sumtv=      10.743462
 sumec=      -40.576389  cor*vef=    -103.540598   ttcor=      62.964210
 rhoeps=      -8.861474     utot=    -139.288955    ehar=     -74.442757

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:    -29.344932     15.636432    -13.708500 sumev=   -2.981208   sumtv=   10.727293

 Kohn-Sham energy:
 sumtv=       10.727293  sumtc=        62.961560   ekin=       73.688852
 rhoep=       -8.860539   utot=      -139.271223   ehks=      -74.442909
 mag. mom=     1.000000  (bands)        1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=3.04e-4  last it=6.76e-4
 mixrho: (warning) scr. and lin-mixed densities had 21975 and 25923 negative points
 AMIX: nmix=2 mmix=8  nelts=93876  beta=0.5  tm=5  rmsdel=1.52e-4
   tj:-0.09781   0.44651
 unscreened rms difference:  smooth  0.000315   local  0.000622
   screened rms difference:  smooth  0.000314   local  0.000622   tot  0.000304

 iors  : write restart file (binary, mesh density)

   it  7  of 10    ehf=       0.552143   ehk=       0.551991
 From last iter    ehf=       0.552161   ehk=       0.551977
 diffe(q)= -0.000017 (0.000304)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999997 ehf=.5521433 ehk=.5519909

 --- BNDFP:  begin iteration 8 of 10 ---

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.112565   -3.943939     2.00     6.00

 Smooth charges: Qmesh = 6.943939  Qgauss = -3.943939  core-nuc = -4  tot = -1

 Average es pot at rmt = -0.110429  avg sphere pot = 0.010709   vconst = 0.110429
 Average es pot at MT boundaries after vconst shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159
 smooth rhoves      5.780784   charge     7.943939
 smvxcm (warning) mesh density negative at 63168 points:  rhomin=-6.94e-5
 smooth rhoeps =   -7.917080 (  -5.008406,  -2.908674)
         rhomu =  -10.366831 (  -6.872212,  -3.494619)
       avg vxc =   -0.149774 (  -0.160975,  -0.138574)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 potential shift to crystal energy zero:    0.000008

 potpus  spin 1 : pnu = 2.925844 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.242413   -1.430720   -1.346208    0.174748      2.7675    0.822103
 1     -0.672474   -1.093387   -0.670752    0.151111     18.5085    0.830172
 2     -0.595336   -0.595336    1.548230    0.359390     16.5957    6.505344
 3     -0.509296   -0.509296    3.497899    0.412089     23.5961   11.576082

 potpus  spin 2 : pnu = 2.923229 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.113063   -1.316279   -1.223247    0.182888      2.7814    0.844709
 1     -0.543355   -1.013559   -0.541452    0.158840     18.7119    0.862693
 2     -0.540717   -0.540717    1.625551    0.362095     16.5218    6.629286
 3     -0.460362   -0.460362    3.563463    0.413330     23.5520   11.669819

 3rd gen potential parameters
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.643553    7.817464    0.075434   -0.647892
 1      3.000000    1.000000    -5.887832    7.317066    0.124430   -0.608609
 2      3.000000    1.000000     6.000000   28.741974    0.486986   -0.090291
 3      3.000000    1.000000     9.000000   39.637353    0.567098   -0.057554

 3rd gen potential parameters spin 2
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.196559    7.678762    0.080232   -0.627097
 1      3.000000    1.000000    -5.887832    7.247983    0.130800   -0.582014
 2      3.000000    1.000000     6.000000   28.999624    0.488467   -0.089009
 3      3.000000    1.000000     9.000000   39.778919    0.567742   -0.057224

 Local terms accumulated over augmentation sites:
                  true           sm,loc      difference
 rhoeps:        -8.844577      -7.899379      -0.945198
 rhomu:         -6.437507      -6.856237       0.418730
 spin2:         -5.204837      -3.487671      -1.717166
 total:        -11.642344     -10.343908      -1.298437
 val*ves       -10.637637      -3.831825      -6.805812
 val*vef       -13.706513     -14.117498       0.410985
 val chg:        2.905346       6.849286      -3.943939
 core chg:       2.000000       2.000000       0.000000
 val mom:        0.953536

 Energy terms:             smooth           local            total
   rhoval*vef            -29.138979        15.415984       -13.722994
   rhoval*ves             -3.656719        -6.805812       -10.462530
   psnuc*ves              15.218286      -283.327204      -268.108918
   utot                    5.780784      -145.066508      -139.285724
   rho*exc                -7.917080        -0.945198        -8.862278
   rho*vxc               -10.366831        -1.298437       -11.665268
   valence chg             6.943939        -3.943939         3.000000
   valence mag             1.983570        -0.992527         0.991042
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.621688,  3 occ states
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.2442 -0.6212 -0.6212 -0.6212  0.2816  0.6466  0.6466  0.6466
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.1150 -0.4921 -0.4921 -0.4921  0.3111  0.6837  0.6837  0.6837
 Est Ef = -0.622 < evl(3)=-0.621 ... using qval=3.0, revise to -0.6212

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.621526;   3.000000 electrons;  D(Ef): 5462.368
         Sum occ. bands:   -2.9806509  incl. Bloechl correction:   -0.000022
         Mag. moment:       1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.905310    6.885894   -3.980585      0.960978    1.949793   -0.988815
       contr. to mm extrapolated for r>rmt:   0.033374 est. true mm = 0.994352

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.976414   -1.242424    2.925844    2.925833    2.500000    2.925833
 spn 2 0    0.972166   -1.113070    2.923229    2.923223    2.500000    2.923223
 1     1    0.956729   -0.616600    2.850000    2.914401    2.250000    2.850000
 spn 2 1    0.000000   -0.465541    2.850000    2.927462    2.250000    2.850000
 2     0    0.000001   -0.853676    3.147584    3.128209    3.147584    3.147584
 spn 2 0    0.000000   -1.285583    3.147584    3.104854    3.147584    3.147584
 3     0    0.000000   -0.830525    4.102416    4.093055    4.102416    4.102416
 spn 2 0    0.000000   -0.830525    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -2.980651  val*vef=     -13.722994   sumtv=      10.742343
 sumec=      -40.576472  cor*vef=    -103.539357   ttcor=      62.962885
 rhoeps=      -8.862278     utot=    -139.285724    ehar=     -74.442774

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:    -29.337679     15.630835    -13.706844 sumev=   -2.980651   sumtv=   10.726193

 Kohn-Sham energy:
 sumtv=       10.726193  sumtc=        62.962071   ekin=       73.688265
 rhoep=       -8.860430   utot=      -139.270748   ehks=      -74.442913
 mag. mom=     1.000000  (bands)        1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=4.15e-4  last it=3.04e-4
 mixrho: (warning) scr. and lin-mixed densities had 24399 and 27563 negative points
 AMIX: nmix=2 mmix=8  nelts=93876  beta=0.5  tm=5  rmsdel=2.08e-4
   tj:-1.91919  -2.47434
 unscreened rms difference:  smooth  0.000434   local  0.000869
   screened rms difference:  smooth  0.000434   local  0.000869   tot  0.000415

 iors  : write restart file (binary, mesh density)

   it  8  of 10    ehf=       0.552126   ehk=       0.551987
 From last iter    ehf=       0.552143   ehk=       0.551991
 diffe(q)= -0.000018 (0.000415)    tol= 0.000010 (0.000500)   more=T
i zbak=1 mmom=.9999997 ehf=.5521257 ehk=.5519866

 --- BNDFP:  begin iteration 9 of 10 ---

 Energy for background charge q=1, radius r=6.204 :  E = 9/5*q*q/r = 0.2902

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -1.148986   -4.073050     2.00     6.00

 Smooth charges: Qmesh = 7.07305  Qgauss = -4.073050  core-nuc = -4  tot = -1

 Average es pot at rmt = -0.113527  avg sphere pot = 0.010641   vconst = 0.113527
 Average es pot at MT boundaries after vconst shift
 Site    ves
   1   0.000000  |
 cell interaction energy from homogeneous background (q=1) is 0.290159
 smooth rhoves      5.818372   charge     8.073050
 smooth rhoeps =   -8.135885 (  -5.126830,  -3.009055)
         rhomu =  -10.654177 (  -7.031494,  -3.622684)
       avg vxc =   -0.172467 (  -0.184335,  -0.160598)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 potential shift to crystal energy zero:    0.000008

 potpus  spin 1 : pnu = 2.925833 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.243214   -1.431056   -1.346777    0.174532      2.7667    0.821857
 1     -0.673094   -1.092641   -0.671376    0.150937     18.4911    0.829752
 2     -0.593619   -0.593619    1.548436    0.359219     16.5998    6.499008
 3     -0.507140   -0.507140    3.498509    0.411979     23.5997   11.569101

 potpus  spin 2 : pnu = 2.923223 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.112975   -1.315549   -1.222848    0.182595      2.7804    0.844342
 1     -0.543047   -1.011312   -0.541150    0.158600     18.6914    0.862043
 2     -0.537424   -0.537424    1.626963    0.361883     16.5268    6.621292
 3     -0.456525   -0.456525    3.565361    0.413192     23.5564   11.661052

 3rd gen potential parameters
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.641633    7.825636    0.075374   -0.648206
 1      3.000000    1.000000    -5.887832    7.323078    0.124286   -0.609037
 2      3.000000    1.000000     6.000000   28.728056    0.486873   -0.090368
 3      3.000000    1.000000     9.000000   39.626115    0.567031   -0.057582

 3rd gen potential parameters spin 2
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -12.195546    7.688764    0.080142   -0.627514
 1      3.000000    1.000000    -5.887832    7.254851    0.130601   -0.582595
 2      3.000000    1.000000     6.000000   28.982272    0.488327   -0.089102
 3      3.000000    1.000000     9.000000   39.764872    0.567658   -0.057259

 Local terms accumulated over augmentation sites:
                  true           sm,loc      difference
 rhoeps:        -8.842077      -8.120213      -0.721865
 rhomu:         -6.438884      -7.016816       0.577932
 spin2:         -5.200229      -3.617059      -1.583170
 total:        -11.639113     -10.633875      -1.005238
 val*ves       -10.628550      -3.803852      -6.824698
 val*vef       -13.694433     -14.379379       0.684946
 val chg:        2.899716       6.972766      -4.073050
 core chg:       2.000000       2.000000       0.000000
 val mom:        0.959788

 Energy terms:             smooth           local            total
   rhoval*vef            -29.927470        16.222534       -13.704936
   rhoval*ves             -3.622177        -6.824698       -10.446875
   psnuc*ves              15.258922      -283.354181      -268.095260
   utot                    5.818372      -145.089440      -139.271067
   rho*exc                -8.135885        -0.721865        -8.857750
   rho*vxc               -10.654177        -1.005238       -11.659415
   valence chg             7.073050        -4.073050         3.000000
   valence mag             1.992525        -0.994085         0.998440
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     3.00000   cores     2.00000   nucleii    -6.00000
    hom background     1.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.621526,  3 occ states
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.2448 -0.6217 -0.6217 -0.6217  0.2803  0.6483  0.6483  0.6483
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.1146 -0.4915 -0.4915 -0.4915  0.3124  0.6885  0.6885  0.6885

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.621998;   3.000000 electrons;  D(Ef): 4537.781
         Sum occ. bands:   -2.9813771  incl. Bloechl correction:   -0.000021
         Mag. moment:       1.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    2.905850    6.949339   -4.043489      0.960948    1.932397   -0.971449
       contr. to mm extrapolated for r>rmt:   0.033391 est. true mm = 0.994339

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.976583   -1.243056    2.925833    2.925990    2.500000    2.925990
 spn 2 0    0.972451   -1.112717    2.923223    2.923466    2.500000    2.923466
 1     1    0.956816   -0.617338    2.850000    2.914404    2.250000    2.850000
 spn 2 1    0.000000   -0.460186    2.850000    2.931446    2.250000    2.850000
 2     0    0.000001   -0.851830    3.147584    3.128210    3.147584    3.147584
 spn 2 0    0.000000   -1.285358    3.147584    3.104727    3.147584    3.147584
 3     0    0.000000   -0.827572    4.102416    4.093074    4.102416    4.102416
 spn 2 0    0.000000   -0.827572    4.102416    4.102416    4.102416    4.102416

 Harris energy:
 sumev=       -2.981377  val*vef=     -13.704936   sumtv=      10.723559
 sumec=      -40.579946  cor*vef=    -103.542424   ttcor=      62.962478
 rhoeps=      -8.857750     utot=    -139.271067    ehar=     -74.442781

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:    -29.797850     16.081049    -13.716801 sumev=   -2.981377   sumtv=   10.735424

 Kohn-Sham energy:
 sumtv=       10.735424  sumtc=        62.963529   ekin=       73.698953
 rhoep=       -8.861805   utot=      -139.280043   ehks=      -74.442895
 mag. mom=     1.000000  (bands)        1.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.199
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=3.38e-4  last it=4.15e-4
 AMIX: nmix=2 mmix=8  nelts=93876  beta=0.5  tm=5  rmsdel=1.69e-4
   tj: 2.91274  -3.06524
 unscreened rms difference:  smooth  0.000261   local  0.000821
   screened rms difference:  smooth  0.000264   local  0.000821   tot  0.000338

 iors  : write restart file (binary, mesh density)

   it  9  of 10    ehf=       0.552119   ehk=       0.552005
 From last iter    ehf=       0.552126   ehk=       0.551987
 diffe(q)= -0.000007 (0.000338)    tol= 0.000010 (0.000500)   more=F
c zbak=1 mmom=.9999997 ehf=.552119 ehk=.5520054
 Exit 0 LMF 
 CPU time:  33.120s   Wall clock 33.158s  at  10:21:18 10.08.2018  on  localhost.localdomai
