rdcmd:  lmfa al
 ----------------------  START LMFA  ----------------------

 rdctrl: reset global max nl from 4 to 3

 LMFA:     nbas = 1  nspec = 1  verb 31,30
 pot:      XC:BH
 autogen:  none

                Plat                                  Qlat
   0.000000   0.500000   0.500000       -1.000000   1.000000   1.000000
   0.500000   0.000000   0.500000        1.000000  -1.000000   1.000000
   0.500000   0.500000   0.000000        1.000000   1.000000  -1.000000
   alat = 7.606  Cell vol = 110.004125

 LATTC:  as= 2.000   tol=1.00e-8   alat= 7.60600   awald= 0.417   lmax=6
         r1=  1.931   nkd= 135      q1=  5.817   nkg= 181

 SGROUP: 1 symmetry operations from 0 generators
 SYMLAT: Bravais system is cubic with 48 symmetry operations.
 SYMCRY: crystal invariant under 48 symmetry operations for tol=1e-5
 GROUPG: the following are sufficient to generate the space group:
         i*r3(1,1,-1) r4x
         i*r3(1,1,-1) r4x
 MKSYM:  found 48 space group operations ... includes inversion

 FREEAT: G tolerance for envelope functions : 1e-6

 Species A:  Z=13  Qc=10  R=2.704881  Q=0
 mesh:   rmt=2.704881  rmax=49.167456  a=0.025  nr=345  nr(rmax)=461
  Pl=  3.5     3.5     3.5    
  Ql=  2.0     1.0     0.0    

  iter     qint         drho          vh0          rho0          vsum     beta
    1   13.000000   1.171E+03       65.0000    0.6466E+02      -26.1241   0.30
   47   13.000000   4.662E-05       89.0963    0.1518E+04     -138.2960   0.30


 sumev=-1.374924  etot=-483.532533  eref=-483.546300  diff= 0.013767

 Free-atom wavefunctions:
 valence:      eval       node at      max at       c.t.p.   rho(r>rmt)
   3s      -0.58196         0.803       2.063       3.192     0.373986
   3p      -0.21101         0.818       2.636       5.253     0.639570
   3d       0.01256         0.000      33.121      49.167*    0.999927

 core:        ecore       node at      max at       c.t.p.   rho(r>rmt)
   1s    -110.54817         0.000       0.079       0.158     0.000000
   2s      -7.90179         0.165       0.498       0.798     0.000016
   2p      -5.12454         0.000       0.428       1.006     0.000094

 Optimise free-atom basis for species A, rmt=2.704881
 l  it    Rsm      Eh     stiffR   stiffE      Eval      Exact     Pnu    Ql
 0  29   2.066  -0.447     123.5    231.2   -0.58190  -0.58196    3.84   2.00
 1  14   2.016  -0.108     193.4   1397.3   -0.20955  -0.21101    3.76   1.00
 eigenvalue sum:  exact  -1.37492    opt basis  -1.37336    error 0.00157

 tailsm: fit tails to 6 smoothed hankels, rmt= 2.70488, rsm= 1.35244
 HNSMFT: 103 points in interval  2.70488  34.64757;  q=  1.387552
 E:    -1.00000    -2.00000    -4.00000    -6.00000    -9.00000    -15.0000
 C:    1.164364    1.066594    77.62984    -713.147    3873.727    -35363.7
        r          rho         fit         diff
    2.704881    0.011169    0.011172   -0.000003
    3.473278    0.004189    0.004189    0.000000
    4.459918    0.001131    0.001132    0.000000
    5.726790    0.000221    0.000221    0.000000
    7.353485    0.000031    0.000031    0.000000
    9.442203    0.000003    0.000003    0.000000
    q(fit):     1.387552    rms diff:   0.000002
    fit: r>rmt  1.387552   r<rmt  2.045825   qtot  3.433377
    rho: r>rmt  1.387552   r<rmt  1.612448   qtot  3.000000

 coretail: q=3.59e-4, rho(rmt)=5.96e-4.  Fit with Hankel e=-22.526  coeff=364.|
      r            rhoc          fit
    2.704881    0.00262382    0.00262382
    3.065098    0.00053752    0.00053796
    3.473278    0.00008843    0.00008784
    3.935805    0.00001135    0.00001108
    4.459918    0.00000110    0.00000104
    5.053816    0.00000008    0.00000007

 Sum of reference energies: -483.5463
 Exit 0 LMFA 
 CPU time:  0.094s   Wall clock 0.128s  at  19:07:44 15.08.2018  on  localhost.localdomai
rdcmd:  lmf al
 ----------------------  START LMF  -----------------------

 rdctrl: reset global max nl from 4 to 3

 LMF:      nbas = 1  nspec = 1  verb 31,30
 pot:      XC:BH
 float:    float P LDA-style
 autoread: none
 bz:       metal(2), tetra, invit 

                Plat                                  Qlat
   0.000000   0.500000   0.500000       -1.000000   1.000000   1.000000
   0.500000   0.000000   0.500000        1.000000  -1.000000   1.000000
   0.500000   0.500000   0.000000        1.000000   1.000000  -1.000000
   alat = 7.606  Cell vol = 110.004125

 LATTC:  as= 2.000   tol=1.00e-8   alat= 7.60600   awald= 0.417   lmax=6
         r1=  1.931   nkd= 135      q1=  5.817   nkg= 181

 SGROUP: 1 symmetry operations from 0 generators
 SYMLAT: Bravais system is cubic with 48 symmetry operations.
 SYMCRY: crystal invariant under 48 symmetry operations for tol=1e-5
 GROUPG: the following are sufficient to generate the space group:
         i*r3(1,1,-1) r4x
         i*r3(1,1,-1) r4x
 MKSYM:  found 48 space group operations ... includes inversion
 BZMESH: 72 irreducible QP from 1728 ( 12 12 12 )  shift= F F F
 TETIRR: sorting 10368 tetrahedra ... 247 inequivalent ones found

 species data:  augmentation                           density
 spec       rmt   rsma lmxa kmxa      lmxl     rg   rsmv  kmxv foca   rfoca
 A        2.705  1.082    2    4         2  0.676  1.352    15    2   1.082

 GVLIST: gmax = 5.841 a.u. created 339 vectors of 1000 (33%)
         mesh has 10 x 10 x 10 divisions; length 0.538, 0.538, 0.538
 SGVSYM: 20 symmetry stars found for 339 reciprocal lattice vectors
 

 Makidx:  hamiltonian dimensions Low, Int, High, Negl: 9 0 0 0
 suham :  9 augmentation channels, 9 local potential channels  Maximum lmxa=2

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol=1.0e-6
 spec      l    rsm    eh     gmax    last term   cutoff
  A        0*   1.80  -0.10   4.130    1.73E-06     137 
  A        1*   1.80  -0.10   4.344    1.42E-06     169 
  A        2    1.80  -0.10   4.561    6.08E-06     169 
 

 lmfp  : no rst file ... try to overlap atomic densities
 rdovfa: read and overlap free-atom densities (mesh density) ...
 rdovfa: expected A,       read A        with rmt=  2.7049  mesh   345  0.025

 ovlpfa: overlap smooth part of FA densities
 site   1  spec  1  pos  0.0000  0.0000  0.0000  Qsmooth 3.433377
 total smooth Q = 3.433377

 Free atom and overlapped crystal site charges  spec 1  site 1:
   ib    true(FA)    smooth(FA)  true(OV)    smooth(OV)    local
    1    1.612448    2.045825    2.390247    2.823624   -0.433377

 Smooth charge on mesh:            3.433377
 Sum of local charges:            -0.433377
 Total valence charge:             3.000000
 Sum of core charges:             10.000000
 Sum of nuclear charges:         -13.000000
 Homogeneous background:           0.000000
 Deviation from neutrality:        0.000000

 --- BNDFP:  begin iteration 1 of 10 ---

 Average es pot at rmt = 0.480973  avg sphere pot = 0.205739   vconst = -0.480973
 smooth rhoves      6.171491   charge     3.433377
 smooth rhoeps =   -1.966199   rhomu =   -2.554857  avg vxc =   -0.724066 
   foca rhoeps =   -0.035122   rhomu =   -0.010654  charge  =    0.042406

 locpot:

 site  1  z= 13.0  rmt= 2.70488  nr=345   a=0.025  nlml= 9  rg=0.676  Vfloat=T
 sm core charge = 0.04181 (sphere) + 0.000596 (spillout) = 0.042406
 potential shift to crystal energy zero:    0.001958
 foca xc integrals for spillout charge:   -0.000404   -0.000121

 Energy terms:             smooth           local            total
   rhoval*vef             -5.184711        -2.310130        -7.494841
   rhoval*ves             -0.741500        -4.219809        -4.961308
   psnuc*ves              13.084482     -1874.393856     -1861.309373
   utot                    6.171491      -939.306832      -933.135341
   rho*exc                -2.001321       -33.444264       -35.445585
   rho*vxc                -2.600633       -44.184704       -46.785337
   valence chg             3.433377        -0.433377         3.000000
   core charge            10.000000         0.000000        10.000000

 Charges:  valence     3.00000   cores    10.00000   nucleii   -13.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Incompatible or missing qp weights file ...

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 Start first of two band passes ...
 bndfp:  kpt 1 of 72, k=  0.00000  0.00000  0.00000
 -0.7902  0.9681  0.9681  0.9681  1.1189  1.1189  1.1189  1.4997  1.4997
 bndfp:  kpt 11 of 72, k=  -0.25000  0.25000  0.41667
 -0.5946  0.1492  0.7165  0.8573  0.8822  1.2162  1.3140  1.4373  1.7557
 bndfp:  kpt 21 of 72, k=  -0.16667  0.16667  0.50000
 -0.5901  0.2543  0.6238  0.6900  0.7668  1.1444  1.4189  1.5664  1.8455
 bndfp:  kpt 31 of 72, k=  -0.25000  0.25000  0.75000
 -0.3433 -0.0256  0.3046  0.5569  0.6985  1.2583  1.4605  1.7313  2.0109
 bndfp:  kpt 41 of 72, k=  -0.08333  0.08333  0.91667
 -0.2417 -0.0081  0.2756  0.4634  0.5454  0.8014  1.9312  2.0843  2.0984
 bndfp:  kpt 51 of 72, k=  -0.16667  0.33333  0.66667
 -0.4099  0.0126  0.4318  0.4636  0.8833  1.2779  1.3169  1.6933  1.9792
 bndfp:  kpt 61 of 72, k=  0.00000  0.16667  0.83333
 -0.3219  0.1140  0.2602  0.3689  0.6823  0.8602  1.7689  2.0337  2.0821
 bndfp:  kpt 71 of 72, k=  0.00000  0.33333  1.00000
 -0.1178 -0.0193  0.0876  0.2089  0.9322  1.0746  1.5857  1.6676  2.1589

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.032264;   3.000000 electrons;  D(Ef):    4.201
         Sum occ. bands:   -0.8919469  incl. Bloechl correction:   -0.007317

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 72, k=  0.00000  0.00000  0.00000
 -0.7902  0.9681  0.9681  0.9681  1.1189  1.1189  1.1189  1.4997  1.4997
 bndfp:  kpt 11 of 72, k=  -0.25000  0.25000  0.41667
 -0.5946  0.1492  0.7165  0.8573  0.8822  1.2162  1.3140  1.4373  1.7557
 bndfp:  kpt 21 of 72, k=  -0.16667  0.16667  0.50000
 -0.5901  0.2543  0.6238  0.6900  0.7668  1.1444  1.4189  1.5664  1.8455
 bndfp:  kpt 31 of 72, k=  -0.25000  0.25000  0.75000
 -0.3433 -0.0256  0.3046  0.5569  0.6985  1.2583  1.4605  1.7313  2.0109
 bndfp:  kpt 41 of 72, k=  -0.08333  0.08333  0.91667
 -0.2417 -0.0081  0.2756  0.4634  0.5454  0.8014  1.9312  2.0843  2.0984
 bndfp:  kpt 51 of 72, k=  -0.16667  0.33333  0.66667
 -0.4099  0.0126  0.4318  0.4636  0.8833  1.2779  1.3169  1.6933  1.9792
 bndfp:  kpt 61 of 72, k=  0.00000  0.16667  0.83333
 -0.3219  0.1140  0.2602  0.3689  0.6823  0.8602  1.7689  2.0337  2.0821
 bndfp:  kpt 71 of 72, k=  0.00000  0.33333  1.00000
 -0.1178 -0.0193  0.0876  0.2089  0.9322  1.0746  1.5857  1.6676  2.1589

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.032264;   3.000000 electrons;  D(Ef):    4.201
         Sum occ. bands:   -0.8919469  incl. Bloechl correction:   -0.007317

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1    2.288589    2.510196   -0.221606

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.947031   -0.399419    3.810000    3.797762    3.500000    3.797762
 1     0    1.086803   -0.239323    3.590000    3.570641    3.250000    3.570641
 2     0    0.254755   -0.157806    3.250000    3.244142    3.147584    3.250000

 Harris energy:
 sumev=       -0.891947  val*vef=      -7.494841   sumtv=       6.602894
 sumec=        0.000000  cor*vef=       0.000000   ttcor=     478.138810
 rhoeps=     -35.445585     utot=    -933.135341    ehar=    -483.839221

 srhov:     -4.178746     -2.807921     -6.986667 sumev=   -0.891947   sumtv=    6.094720

 Kohn-Sham energy:
 sumtv=        6.094720  sumtc=       478.138810   ekin=      484.233530
 rhoep=      -35.402215   utot=      -932.667958   ehks=     -483.836643
  
 mixing: mode=A  nmix=3  beta=1  elind=.867
 mixrho:  sought 3 iter from file mixm; read 0.  RMS DQ=9.45e-3
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.433377      3.221606      3.221606      0.006444      3.221606
 site    1   -0.433377     -0.221606     -0.221606      0.004870     -0.221606
 AMIX: nmix=0 mmix=8  nelts=1192  beta=1  tm=5  rmsdel=9.45e-3
 unscreened rms difference:  smooth  0.007252   local  0.004870
   screened rms difference:  smooth  0.006444   local  0.004870   tot  0.009452

 iors  : write restart file (binary, mesh density)

   it  1  of 10    ehf=      -0.292921   ehk=      -0.290343
h ehf=-.2929211 ehk=-.2903427

 --- BNDFP:  begin iteration 2 of 10 ---

 Average es pot at rmt = 0.493647  avg sphere pot = 0.213332   vconst = -0.493647
 smooth rhoves      5.982895   charge     3.221606
 smooth rhoeps =   -1.785114   rhomu =   -2.318791  avg vxc =   -0.715105 
   foca rhoeps =   -0.032121   rhomu =   -0.009705  charge  =    0.042406

 locpot:

 site  1  z= 13.0  rmt= 2.70488  nr=345   a=0.025  nlml= 9  rg=0.676  Vfloat=T
 sm core charge = 0.04181 (sphere) + 0.000596 (spillout) = 0.042406
 potential shift to crystal energy zero:    0.003616
 foca xc integrals for spillout charge:   -0.000426   -0.000129

 Energy terms:             smooth           local            total
   rhoval*vef             -4.270152        -2.950677        -7.220829
   rhoval*ves             -0.990787        -3.709944        -4.700731
   psnuc*ves              12.956576     -1873.880152     -1860.923575
   utot                    5.982895      -938.795048      -932.812153
   rho*exc                -1.817235       -33.607217       -35.424452
   rho*vxc                -2.360617       -44.396934       -46.757551
   valence chg             3.221606        -0.221606         3.000000
   core charge            10.000000         0.000000        10.000000

 Charges:  valence     3.00000   cores    10.00000   nucleii   -13.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=0.032264,  3 occ states

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 72, k=  0.00000  0.00000  0.00000
 -0.8130  0.9438  0.9438  0.9438  1.0913  1.0913  1.0913  1.4731  1.4731
 bndfp:  kpt 11 of 72, k=  -0.25000  0.25000  0.41667
 -0.6173  0.1259  0.6919  0.8318  0.8586  1.1907  1.2870  1.4117  1.7281
 bndfp:  kpt 21 of 72, k=  -0.16667  0.16667  0.50000
 -0.6128  0.2310  0.5993  0.6662  0.7426  1.1195  1.3926  1.5386  1.8168
 bndfp:  kpt 31 of 72, k=  -0.25000  0.25000  0.75000
 -0.3661 -0.0487  0.2812  0.5326  0.6744  1.2333  1.4345  1.7033  1.9811
 bndfp:  kpt 41 of 72, k=  -0.08333  0.08333  0.91667
 -0.2645 -0.0309  0.2522  0.4390  0.5218  0.7767  1.9025  2.0532  2.0672
 bndfp:  kpt 51 of 72, k=  -0.16667  0.33333  0.66667
 -0.4326 -0.0105  0.4096  0.4380  0.8589  1.2522  1.2919  1.6654  1.9497
 bndfp:  kpt 61 of 72, k=  0.00000  0.16667  0.83333
 -0.3446  0.0913  0.2371  0.3441  0.6586  0.8355  1.7413  2.0035  2.0509
 bndfp:  kpt 71 of 72, k=  0.00000  0.33333  1.00000
 -0.1406 -0.0417  0.0647  0.1838  0.9081  1.0498  1.5593  1.6414  2.1267

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.009264;   3.000000 electrons;  D(Ef):    4.227
         Sum occ. bands:   -0.9603868  incl. Bloechl correction:   -0.007349

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1    2.292741    2.507277   -0.214536

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.954111   -0.421773    3.797762    3.799333    3.500000    3.799333
 1     0    1.085976   -0.261533    3.570641    3.571433    3.250000    3.571433
 2     0    0.252654   -0.181057    3.250000    3.243513    3.147584    3.250000

 Harris energy:
 sumev=       -0.960387  val*vef=      -7.220829   sumtv=       6.260442
 sumec=        0.000000  cor*vef=       0.000000   ttcor=     478.138810
 rhoeps=     -35.424452     utot=    -932.812153    ehar=    -483.837353

 srhov:     -4.187552     -2.915933     -7.103485 sumev=   -0.960387   sumtv=    6.143099

 Kohn-Sham energy:
 sumtv=        6.143099  sumtc=       478.138810   ekin=      484.281909
 rhoep=      -35.405874   utot=      -932.713192   ehks=     -483.837156
  
 mixing: mode=A  nmix=3  beta=1  elind=.867
 mixrho:  sought 3 iter from file mixm; read 1.  RMS DQ=6.06e-4  last it=9.45e-3
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.221606      3.214536      3.214536      0.000560      3.214536
 site    1   -0.221606     -0.214536     -0.214536      0.000210     -0.214536
 AMIX: nmix=1 mmix=8  nelts=1192  beta=1  tm=5  rmsdel=6.06e-4
   tj:-0.03165
 unscreened rms difference:  smooth  0.000755   local  0.000210
   screened rms difference:  smooth  0.000560   local  0.000210   tot  0.000606

 iors  : write restart file (binary, mesh density)

   it  2  of 10    ehf=      -0.291053   ehk=      -0.290856
 From last iter    ehf=      -0.292921   ehk=      -0.290343
 diffe(q)=  0.001868 (0.000606)    tol= 0.000001 (0.000010)   more=T
i ehf=-.2910535 ehk=-.2908564

 --- BNDFP:  begin iteration 3 of 10 ---

 Average es pot at rmt = 0.496957  avg sphere pot = 0.214026   vconst = -0.496957
 smooth rhoves      6.008846   charge     3.214312
 smooth rhoeps =   -1.778414   rhomu =   -2.310052  avg vxc =   -0.714921 
   foca rhoeps =   -0.031851   rhomu =   -0.009620  charge  =    0.042406

 locpot:

 site  1  z= 13.0  rmt= 2.70488  nr=345   a=0.025  nlml= 9  rg=0.676  Vfloat=T
 sm core charge = 0.04181 (sphere) + 0.000596 (spillout) = 0.042406
 potential shift to crystal energy zero:    0.003614
 foca xc integrals for spillout charge:   -0.000428   -0.000129

 Energy terms:             smooth           local            total
   rhoval*vef             -4.216404        -2.939428        -7.155831
   rhoval*ves             -0.974614        -3.672060        -4.646673
   psnuc*ves              12.992306     -1873.838281     -1860.845975
   utot                    6.008846      -938.755171      -932.746324
   rho*exc                -1.810265       -33.600628       -35.410893
   rho*vxc                -2.351523       -44.388114       -46.739637
   valence chg             3.214312        -0.214312         3.000000
   core charge            10.000000         0.000000        10.000000

 Charges:  valence     3.00000   cores    10.00000   nucleii   -13.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=0.009264,  3 occ states

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 72, k=  0.00000  0.00000  0.00000
 -0.8175  0.9389  0.9389  0.9389  1.0841  1.0841  1.0841  1.4675  1.4675
 bndfp:  kpt 11 of 72, k=  -0.25000  0.25000  0.41667
 -0.6220  0.1210  0.6865  0.8258  0.8537  1.1849  1.2802  1.4059  1.7218
 bndfp:  kpt 21 of 72, k=  -0.16667  0.16667  0.50000
 -0.6175  0.2260  0.5938  0.6609  0.7374  1.1136  1.3866  1.5318  1.8102
 bndfp:  kpt 31 of 72, k=  -0.25000  0.25000  0.75000
 -0.3708 -0.0536  0.2761  0.5271  0.6691  1.2274  1.4286  1.6965  1.9745
 bndfp:  kpt 41 of 72, k=  -0.08333  0.08333  0.91667
 -0.2692 -0.0357  0.2473  0.4335  0.5168  0.7709  1.8959  2.0460  2.0601
 bndfp:  kpt 51 of 72, k=  -0.16667  0.33333  0.66667
 -0.4374 -0.0154  0.4050  0.4321  0.8535  1.2460  1.2862  1.6586  1.9430
 bndfp:  kpt 61 of 72, k=  0.00000  0.16667  0.83333
 -0.3494  0.0865  0.2323  0.3386  0.6536  0.8295  1.7350  1.9966  2.0436
 bndfp:  kpt 71 of 72, k=  0.00000  0.33333  1.00000
 -0.1454 -0.0465  0.0599  0.1782  0.9029  1.0438  1.5534  1.6352  2.1194

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.004421;   3.000000 electrons;  D(Ef):    4.235
         Sum occ. bands:   -0.9746794  incl. Bloechl correction:   -0.007362

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1    2.294015    2.506676   -0.212661

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.956640   -0.426243    3.799333    3.799966    3.500000    3.799966
 1     0    1.085533   -0.266125    3.571433    3.571894    3.250000    3.571894
 2     0    0.251843   -0.186019    3.250000    3.243337    3.147584    3.250000

 Harris energy:
 sumev=       -0.974679  val*vef=      -7.155831   sumtv=       6.181152
 sumec=        0.000000  cor*vef=       0.000000   ttcor=     478.138810
 rhoeps=     -35.410893     utot=    -932.746324    ehar=    -483.837254

 srhov:     -4.199647     -2.930472     -7.130118 sumev=   -0.974679   sumtv=    6.155439

 Kohn-Sham energy:
 sumtv=        6.155439  sumtc=       478.138810   ekin=      484.294249
 rhoep=      -35.406754   utot=      -932.724741   ehks=     -483.837246
  
 mixing: mode=A  nmix=3  beta=1  elind=.867
 mixrho:  sought 3 iter from file mixm; read 2.  RMS DQ=1.42e-4  last it=6.06e-4
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.214312      3.212661      3.212661      0.000107      3.212661
 site    1   -0.214312     -0.212661     -0.212661      0.000056     -0.212661
 AMIX: nmix=2 mmix=8  nelts=1192  beta=1  tm=5  rmsdel=1.42e-4
   tj:-0.31405   0.00132
 unscreened rms difference:  smooth  0.000146   local  0.000056
   screened rms difference:  smooth  0.000107   local  0.000056   tot  0.000142

 iors  : write restart file (binary, mesh density)

   it  3  of 10    ehf=      -0.290954   ehk=      -0.290946
 From last iter    ehf=      -0.291053   ehk=      -0.290856
 diffe(q)=  0.000099 (0.000142)    tol= 0.000001 (0.000010)   more=T
i ehf=-.2909543 ehk=-.2909459

 --- BNDFP:  begin iteration 4 of 10 ---

 Average es pot at rmt = 0.497685  avg sphere pot = 0.214252   vconst = -0.497685
 smooth rhoves      6.014415   charge     3.212084
 smooth rhoeps =   -1.776450   rhomu =   -2.307492  avg vxc =   -0.714849 
   foca rhoeps =   -0.031779   rhomu =   -0.009597  charge  =    0.042406

 locpot:

 site  1  z= 13.0  rmt= 2.70488  nr=345   a=0.025  nlml= 9  rg=0.676  Vfloat=T
 sm core charge = 0.04181 (sphere) + 0.000596 (spillout) = 0.042406
 potential shift to crystal energy zero:    0.003613
 foca xc integrals for spillout charge:   -0.000429   -0.000129

 Energy terms:             smooth           local            total
   rhoval*vef             -4.201406        -2.934154        -7.135560
   rhoval*ves             -0.971240        -3.658569        -4.629809
   psnuc*ves              13.000071     -1873.823584     -1860.823513
   utot                    6.014415      -938.741076      -932.726661
   rho*exc                -1.808229       -33.598551       -35.406780
   rho*vxc                -2.348868       -44.385334       -46.734202
   valence chg             3.212084        -0.212084         3.000000
   core charge            10.000000         0.000000        10.000000

 Charges:  valence     3.00000   cores    10.00000   nucleii   -13.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=0.004421,  3 occ states

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 72, k=  0.00000  0.00000  0.00000
 -0.8187  0.9377  0.9377  0.9377  1.0822  1.0822  1.0822  1.4661  1.4661
 bndfp:  kpt 11 of 72, k=  -0.25000  0.25000  0.41667
 -0.6232  0.1197  0.6851  0.8243  0.8524  1.1834  1.2784  1.4045  1.7201
 bndfp:  kpt 21 of 72, k=  -0.16667  0.16667  0.50000
 -0.6186  0.2248  0.5924  0.6596  0.7361  1.1121  1.3851  1.5300  1.8085
 bndfp:  kpt 31 of 72, k=  -0.25000  0.25000  0.75000
 -0.3720 -0.0548  0.2749  0.5257  0.6678  1.2259  1.4270  1.6947  1.9728
 bndfp:  kpt 41 of 72, k=  -0.08333  0.08333  0.91667
 -0.2704 -0.0369  0.2461  0.4322  0.5156  0.7694  1.8943  2.0441  2.0583
 bndfp:  kpt 51 of 72, k=  -0.16667  0.33333  0.66667
 -0.4385 -0.0166  0.4039  0.4306  0.8521  1.2444  1.2848  1.6569  1.9413
 bndfp:  kpt 61 of 72, k=  0.00000  0.16667  0.83333
 -0.3506  0.0853  0.2311  0.3372  0.6523  0.8280  1.7333  1.9948  2.0417
 bndfp:  kpt 71 of 72, k=  0.00000  0.33333  1.00000
 -0.1465 -0.0476  0.0587  0.1768  0.9016  1.0423  1.5519  1.6337  2.1175

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.003226;   3.000000 electrons;  D(Ef):    4.238
         Sum occ. bands:   -0.9781971  incl. Bloechl correction:   -0.007365

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1    2.294378    2.506510   -0.212132

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.957378   -0.427322    3.799966    3.800150    3.500000    3.800150
 1     0    1.085389   -0.267254    3.571894    3.572022    3.250000    3.572022
 2     0    0.251611   -0.187249    3.250000    3.243286    3.147584    3.250000

 Harris energy:
 sumev=       -0.978197  val*vef=      -7.135560   sumtv=       6.157363
 sumec=        0.000000  cor*vef=       0.000000   ttcor=     478.138810
 rhoeps=     -35.406780     utot=    -932.726661    ehar=    -483.837267

 srhov:     -4.202325     -2.934760     -7.137084 sumev=   -0.978197   sumtv=    6.158887

 Kohn-Sham energy:
 sumtv=        6.158887  sumtc=       478.138810   ekin=      484.297698
 rhoep=      -35.406999   utot=      -932.727966   ehks=     -483.837267
  
 mixing: mode=A  nmix=3  beta=1  elind=.867
 mixrho:  sought 3 iter from file mixm; read 3.  RMS DQ=6.91e-6  last it=1.42e-4
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.212084      3.212132      3.212132      0.000008      3.212132
 site    1   -0.212084     -0.212132     -0.212132      0.000001     -0.212132
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=1192  beta=1  tm=5  rmsdel=6.91e-6
   tj: 0.04478
 unscreened rms difference:  smooth  0.000011   local  0.000001
   screened rms difference:  smooth  0.000008   local  0.000001   tot  0.000007

 iors  : write restart file (binary, mesh density)

   it  4  of 10    ehf=      -0.290967   ehk=      -0.290967
 From last iter    ehf=      -0.290954   ehk=      -0.290946
 diffe(q)= -0.000013 (0.000007)    tol= 0.000001 (0.000010)   more=T
i ehf=-.2909673 ehk=-.2909672

 --- BNDFP:  begin iteration 5 of 10 ---

 Average es pot at rmt = 0.497629  avg sphere pot = 0.214246   vconst = -0.497629
 smooth rhoves      6.013926   charge     3.212156
 smooth rhoeps =   -1.776514   rhomu =   -2.307574  avg vxc =   -0.714852 
   foca rhoeps =   -0.031783   rhomu =   -0.009599  charge  =    0.042406

 locpot:

 site  1  z= 13.0  rmt= 2.70488  nr=345   a=0.025  nlml= 9  rg=0.676  Vfloat=T
 sm core charge = 0.04181 (sphere) + 0.000596 (spillout) = 0.042406
 potential shift to crystal energy zero:    0.003613
 foca xc integrals for spillout charge:   -0.000428   -0.000129

 Energy terms:             smooth           local            total
   rhoval*vef             -4.202098        -2.934545        -7.136643
   rhoval*ves             -0.971583        -3.659150        -4.630733
   psnuc*ves              12.999435     -1873.824212     -1860.824777
   utot                    6.013926      -938.741681      -932.727755
   rho*exc                -1.808296       -33.598683       -35.406979
   rho*vxc                -2.348955       -44.385510       -46.734465
   valence chg             3.212156        -0.212156         3.000000
   core charge            10.000000         0.000000        10.000000

 Charges:  valence     3.00000   cores    10.00000   nucleii   -13.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=0.003226,  3 occ states

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 72, k=  0.00000  0.00000  0.00000
 -0.8186  0.9378  0.9378  0.9378  1.0823  1.0823  1.0823  1.4662  1.4662
 bndfp:  kpt 11 of 72, k=  -0.25000  0.25000  0.41667
 -0.6231  0.1198  0.6852  0.8243  0.8525  1.1835  1.2785  1.4045  1.7202
 bndfp:  kpt 21 of 72, k=  -0.16667  0.16667  0.50000
 -0.6186  0.2249  0.5925  0.6597  0.7361  1.1122  1.3851  1.5301  1.8085
 bndfp:  kpt 31 of 72, k=  -0.25000  0.25000  0.75000
 -0.3719 -0.0547  0.2749  0.5258  0.6679  1.2260  1.4271  1.6948  1.9729
 bndfp:  kpt 41 of 72, k=  -0.08333  0.08333  0.91667
 -0.2704 -0.0368  0.2462  0.4323  0.5157  0.7695  1.8944  2.0442  2.0584
 bndfp:  kpt 51 of 72, k=  -0.16667  0.33333  0.66667
 -0.4385 -0.0166  0.4040  0.4307  0.8522  1.2445  1.2849  1.6570  1.9414
 bndfp:  kpt 61 of 72, k=  0.00000  0.16667  0.83333
 -0.3505  0.0854  0.2312  0.3372  0.6524  0.8281  1.7334  1.9949  2.0418
 bndfp:  kpt 71 of 72, k=  0.00000  0.33333  1.00000
 -0.1465 -0.0476  0.0588  0.1769  0.9017  1.0424  1.5520  1.6337  2.1176

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.003297;   3.000000 electrons;  D(Ef):    4.238
         Sum occ. bands:   -0.9779843  incl. Bloechl correction:   -0.007365

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1    2.294362    2.506524   -0.212162

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.957363   -0.427246    3.800150    3.800145    3.500000    3.800145
 1     0    1.085380   -0.267190    3.572022    3.572013    3.250000    3.572013
 2     0    0.251620   -0.187180    3.250000    3.243287    3.147584    3.250000

 Harris energy:
 sumev=       -0.977984  val*vef=      -7.136643   sumtv=       6.158659
 sumec=        0.000000  cor*vef=       0.000000   ttcor=     478.138810
 rhoeps=     -35.406979     utot=    -932.727755    ehar=    -483.837265

 srhov:     -4.202146     -2.934588     -7.136735 sumev=   -0.977984   sumtv=    6.158750

 Kohn-Sham energy:
 sumtv=        6.158750  sumtc=       478.138810   ekin=      484.297561
 rhoep=      -35.406989   utot=      -932.727837   ehks=     -483.837265
  
 mixing: mode=A  nmix=3  beta=1  elind=.867
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=5.06e-7  last it=6.91e-6
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.212156      3.212162      3.212162      0.000001      3.212162
 site    1   -0.212156     -0.212162     -0.212162      0.000000     -0.212162
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=1192  beta=1  tm=5  rmsdel=5.06e-7
   tj:-0.07147
 unscreened rms difference:  smooth  0.000001   local  0.000000
   screened rms difference:  smooth  0.000001   local  0.000000   tot  0.000001

 iors  : write restart file (binary, mesh density)

   it  5  of 10    ehf=      -0.290965   ehk=      -0.290965
 From last iter    ehf=      -0.290967   ehk=      -0.290967
 diffe(q)=  0.000002 (0.000001)    tol= 0.000001 (0.000010)   more=T
i ehf=-.2909653 ehk=-.2909653

 --- BNDFP:  begin iteration 6 of 10 ---

 Average es pot at rmt = 0.497626  avg sphere pot = 0.214245   vconst = -0.497626
 smooth rhoves      6.013905   charge     3.212164
 smooth rhoeps =   -1.776518   rhomu =   -2.307580  avg vxc =   -0.714852 
   foca rhoeps =   -0.031783   rhomu =   -0.009599  charge  =    0.042406

 locpot:

 site  1  z= 13.0  rmt= 2.70488  nr=345   a=0.025  nlml= 9  rg=0.676  Vfloat=T
 sm core charge = 0.04181 (sphere) + 0.000596 (spillout) = 0.042406
 potential shift to crystal energy zero:    0.003613
 foca xc integrals for spillout charge:   -0.000428   -0.000129

 Energy terms:             smooth           local            total
   rhoval*vef             -4.202148        -2.934582        -7.136730
   rhoval*ves             -0.971596        -3.659217        -4.630812
   psnuc*ves              12.999405     -1873.824281     -1860.824876
   utot                    6.013905      -938.741749      -932.727844
   rho*exc                -1.808301       -33.598690       -35.406991
   rho*vxc                -2.348962       -44.385518       -46.734480
   valence chg             3.212164        -0.212164         3.000000
   core charge            10.000000         0.000000        10.000000

 Charges:  valence     3.00000   cores    10.00000   nucleii   -13.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=0.003297,  3 occ states

 subzi: tetrahedron integration of bands; tetrahedron integration of density

 bndfp:  kpt 1 of 72, k=  0.00000  0.00000  0.00000
 -0.8186  0.9378  0.9378  0.9378  1.0823  1.0823  1.0823  1.4662  1.4662
 bndfp:  kpt 11 of 72, k=  -0.25000  0.25000  0.41667
 -0.6231  0.1198  0.6852  0.8244  0.8525  1.1835  1.2785  1.4045  1.7202
 bndfp:  kpt 21 of 72, k=  -0.16667  0.16667  0.50000
 -0.6186  0.2249  0.5925  0.6597  0.7361  1.1122  1.3852  1.5301  1.8086
 bndfp:  kpt 31 of 72, k=  -0.25000  0.25000  0.75000
 -0.3719 -0.0547  0.2749  0.5258  0.6679  1.2260  1.4271  1.6948  1.9729
 bndfp:  kpt 41 of 72, k=  -0.08333  0.08333  0.91667
 -0.2703 -0.0368  0.2462  0.4323  0.5157  0.7695  1.8944  2.0442  2.0584
 bndfp:  kpt 51 of 72, k=  -0.16667  0.33333  0.66667
 -0.4385 -0.0166  0.4040  0.4307  0.8522  1.2445  1.2849  1.6570  1.9414
 bndfp:  kpt 61 of 72, k=  0.00000  0.16667  0.83333
 -0.3505  0.0854  0.2312  0.3372  0.6524  0.8281  1.7334  1.9949  2.0418
 bndfp:  kpt 71 of 72, k=  0.00000  0.33333  1.00000
 -0.1465 -0.0475  0.0588  0.1769  0.9017  1.0424  1.5520  1.6337  2.1176

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.003302;   3.000000 electrons;  D(Ef):    4.238
         Sum occ. bands:   -0.9779714  incl. Bloechl correction:   -0.007365

 Saved qp weights ...

 ... Generating total DOS

 mkrout:  Qtrue      sm,loc       local
   1    2.294362    2.506524   -0.212162

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.957360   -0.427243    3.800145    3.800145    3.500000    3.800145
 1     0    1.085381   -0.267185    3.572013    3.572013    3.250000    3.572013
 2     0    0.251620   -0.187175    3.250000    3.243287    3.147584    3.250000

 Harris energy:
 sumev=       -0.977971  val*vef=      -7.136730   sumtv=       6.158759
 sumec=        0.000000  cor*vef=       0.000000   ttcor=     478.138810
 rhoeps=     -35.406991     utot=    -932.727844    ehar=    -483.837265

 srhov:     -4.202135     -2.934583     -7.136718 sumev=   -0.977971   sumtv=    6.158746

 Kohn-Sham energy:
 sumtv=        6.158746  sumtc=       478.138810   ekin=      484.297557
 rhoep=      -35.406989   utot=      -932.727833   ehks=     -483.837265
  
 mixing: mode=A  nmix=3  beta=1  elind=.867
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=7.57e-8  last it=5.06e-7
 charges:       old           new         screened      rms diff       lin mix
 smooth       3.212164      3.212162      3.212162      0.000000      3.212162
 site    1   -0.212164     -0.212162     -0.212162      0.000000     -0.212162
 AMIX: condition of normal eqns >100000. Reducing nmix to 2
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=1192  beta=1  tm=5  rmsdel=7.57e-8
   tj: 0.11462
 unscreened rms difference:  smooth  0.000000   local  0.000000
   screened rms difference:  smooth  0.000000   local  0.000000   tot  0.000000

 iors  : write restart file (binary, mesh density)

   it  6  of 10    ehf=      -0.290965   ehk=      -0.290965
 From last iter    ehf=      -0.290965   ehk=      -0.290965
 diffe(q)=  0.000000 (0.000000)    tol= 0.000001 (0.000010)   more=F
c ehf=-.2909653 ehk=-.2909653
 Exit 0 LMF 
 CPU time:  0.976s   Wall clock 1.014s  at  19:07:45 15.08.2018  on  localhost.localdomai
rdcmd:  lmf al -vnk=72 --quit=rho --pr20
 ----------------------  START LMF  -----------------------

 rdctrl: reset global max nl from 4 to 3

 LMF:      nbas = 1  nspec = 1  verb 20
 pot:      XC:BH
 float:    float P LDA-style
 autoread: none
 bz:       metal(2), tetra, invit 

                Plat                                  Qlat
   0.000000   0.500000   0.500000       -1.000000   1.000000   1.000000
   0.500000   0.000000   0.500000        1.000000  -1.000000   1.000000
   0.500000   0.500000   0.000000        1.000000   1.000000  -1.000000
   alat = 7.606  Cell vol = 110.004125
 GROUPG: the following are sufficient to generate the space group:
         i*r3(1,1,-1) r4x
         i*r3(1,1,-1) r4x
 MKSYM:  found 48 space group operations ... includes inversion
 BZMESH: 8797 irreducible QP from 373248 ( 72 72 72 )  shift= F F F
 SGVSYM: 20 symmetry stars found for 339 reciprocal lattice vectors
 Makidx:  hamiltonian dimensions Low, Int, High, Negl: 9 0 0 0

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol=1.0e-6
 spec      l    rsm    eh     gmax    last term   cutoff
  A        0*   1.80  -0.10   4.130    1.73E-06     137 
  A        1*   1.80  -0.10   4.344    1.42E-06     169 
  A        2    1.80  -0.10   4.561    6.08E-06     169 

 iors  : read restart file (binary, mesh density)
         use from  restart file: ef window, positions, pnu 
         ignore in restart file: *

 --- BNDFP:  begin iteration 1 of 10 ---

 Average es pot at rmt = 0.497626  avg sphere pot = 0.214246   vconst = -0.497626

 site  1  z= 13.0  rmt= 2.70488  nr=345   a=0.025  nlml= 9  rg=0.676  Vfloat=T
 sm core charge = 0.04181 (sphere) + 0.000596 (spillout) = 0.042406
 potential shift to crystal energy zero:    0.003613


 Incompatible or missing qp weights file ...

 Start first of two band passes ...
 BZINTS: Fermi energy:      0.001952;   3.000000 electrons;  D(Ef):    5.503
         Sum occ. bands:   -0.9776318  incl. Bloechl correction:   -0.000244

 Saved qp weights ...
 Start second band pass ...
 BZINTS: Fermi energy:      0.001952;   3.000000 electrons;  D(Ef):    5.503
         Sum occ. bands:   -0.9776318  incl. Bloechl correction:   -0.000244

 Saved qp weights ...

 site    1   species   1:A       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.964315   -0.424308    3.800145    3.801373    3.500000    3.801373
 1     0    1.078247   -0.268681    3.572013    3.570962    3.250000    3.570962
 2     0    0.252023   -0.186429    3.250000    3.243456    3.147584    3.250000

 ekin=484.297902  rho*v=-968.134827  sumev=-0.977632  ehf=-483.836924

 ekin=484.307598  rho*v=-968.144518  ehf=-483.836924  ehks=-483.836920
  
 mixing: mode=A  nmix=3  beta=1  elind=.867
 mixrho:  sought 3 iter from file mixm; read 4.  RMS DQ=1.74e-4
 rms smooth dq=1.09e-4  max local dq=9.47e-5  dq=1.74e-4
 Exit 0 quit = rho 
 CPU time:  22.866s   Wall clock 22.905s  at  19:08:08 15.08.2018  on  localhost.localdomai
rdcmd:
