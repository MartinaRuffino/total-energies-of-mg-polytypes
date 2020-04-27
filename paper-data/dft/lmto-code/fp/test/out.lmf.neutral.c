rdcmd:  lmfa --no-iactiv c -vzbak=0
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
 CPU time:  0.114s   Wall clock 0.159s  at  10:10:13 10.08.2018  on  localhost.localdomai
rdcmd:  lmf  --no-iactiv c -vzbak=0
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

 Smooth charge on mesh:            2.661248    moment    1.332651
 Sum of local charges:             1.338752    moments   0.667349
 Total valence charge:             4.000000    moment    2.000000
 Sum of core charges:              2.000000    moment    0.000000
 Sum of nuclear charges:          -6.000000
 Homogeneous background:           0.000000
 Deviation from neutrality:        0.000000

 --- BNDFP:  begin iteration 1 of 10 ---

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1    0.377655    1.338752     2.00     6.00

 Smooth charges: Qmesh = 2.661248  Qgauss = 1.338752  core-nuc = -4  tot = 0

 Average es pot at rmt = 0.006604  avg sphere pot = 0.019543   vconst = -0.006604
 Average es pot at MT boundaries after vconst shift
 Site    ves
   1   0.000000  |
 smooth rhoves      2.063389   charge     2.661248
 smooth rhoeps =   -1.494494 (  -1.109031,  -0.385462)
         rhomu =   -1.946019 (  -1.522381,  -0.423638)
       avg vxc =   -0.191345 (  -0.218125,  -0.164565)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 potential shift to crystal energy zero:    0.000003

 potpus  spin 1 : pnu = 2.900000 2.850000 3.180000 4.120000
 l        enu         v           c          srdel        qpar        ppar
 0     -1.054715   -1.268528   -1.158163    0.198868      2.7906    0.901142
 1     -0.471084   -1.054270   -0.468773    0.175050     19.1074    0.930215
 2     -0.328679   -0.621647    1.519427    0.353874     17.0974    6.142525
 3     -0.115033   -0.556834    3.369801    0.400983     24.4207   10.665044

 potpus  spin 2 : pnu = 2.900000 2.850000 3.180000 4.120000
 l        enu         v           c          srdel        qpar        ppar
 0     -0.846590   -1.087566   -0.962512    0.211082      2.8067    0.936239
 1     -0.269118   -0.929135   -0.266498    0.186383     19.0749    0.984263
 2     -0.216979   -0.512409    1.659281    0.357599     16.9824    6.304726
 3     -0.011628   -0.454897    3.495382    0.402797     24.3469   10.791758

 3rd gen potential parameters
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000    -9.233051    7.589352    0.102212   -0.581577
 1      3.000000    1.000000    -5.887832    7.119796    0.144160   -0.533282
 2      3.000000    1.000000     4.727244   27.134499    0.465946   -0.095779
 3      3.000000    1.000000     7.577135   37.213068    0.543678   -0.062062

 3rd gen potential parameters spin 2
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000    -9.233051    7.437705    0.107908   -0.555887
 1      3.000000    1.000000    -5.887832    7.130042    0.153492   -0.500464
 2      3.000000    1.000000     4.727244   27.482513    0.468116   -0.093877
 3      3.000000    1.000000     7.577135   37.414392    0.544673   -0.061531

 Local terms accumulated over augmentation sites:
                  true           sm,loc      difference
 rhoeps:        -9.536045      -1.428544      -8.107501
 rhomu:         -7.422765      -1.450396      -5.972369
 spin2:         -5.125221      -0.410188      -4.715033
 total:        -12.547986      -1.860584     -10.687402
 val*ves       -10.246156      -5.064346      -5.181810
 val*vef       -14.282893      -6.924930      -7.357963
 val chg:        3.701844       2.363092       1.338752
 core chg:       2.000000       2.000000       0.000000
 val mom:        1.810990

 Energy terms:             smooth           local            total
   rhoval*vef             -3.930329       -10.432170       -14.362499
   rhoval*ves             -5.058517        -5.181810       -10.240327
   psnuc*ves               9.185295      -278.835495      -269.650200
   utot                    2.063389      -142.008652      -139.945263
   rho*exc                -1.494494        -8.107501        -9.601995
   rho*vxc                -1.946019       -10.687402       -12.633421
   valence chg             2.661248         1.338752         4.000000
   valence mag             1.332651         0.667349         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Incompatible or missing qp weights file ...

 Start first of two band passes ...
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.0407 -0.4290 -0.4290 -0.4290  0.1321  0.5284  0.5284  0.5284
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -0.8387 -0.2374 -0.2374 -0.2374  0.2088  0.6249  0.6249  0.6249

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.430662;   4.000000 electrons;  D(Ef):  926.075
         Sum occ. bands:   -2.7432333  incl. Bloechl correction:   -0.000179
         Mag. moment:       2.000000

 Saved qp weights ...
 Start second band pass ...
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.0407 -0.4290 -0.4290 -0.4290  0.1321  0.5284  0.5284  0.5284
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -0.8387 -0.2374 -0.2374 -0.2374  0.2088  0.6249  0.6249  0.6249
 Est Ef = -0.431 < evl(4)=-0.429 ... using qval=4.0, revise to -0.4290

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.430662;   4.000000 electrons;  D(Ef):  926.075
         Sum occ. bands:   -2.7432333  incl. Bloechl correction:   -0.000179
         Mag. moment:       2.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.686856    3.838100   -0.151244      1.796589    2.141138   -0.344549
       contr. to mm extrapolated for r>rmt:   0.163680 est. true mm = 1.960270

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.957935   -1.038744    2.900000    2.914622    2.500000    2.914622
 spn 2 0    0.945133   -0.836887    2.900000    2.908183    2.500000    2.908183
 1     1    1.783776   -0.429110    2.850000    2.889411    2.250000    2.850000
 spn 2 1    0.000000   -1.101867    2.850000    2.167565    2.250000    2.850000
 2     0    0.000011   -0.813692    3.180000    3.132790    3.147584    3.147584
 spn 2 0    0.000000   -1.175534    3.180000    3.108512    3.147584    3.147584
 3     0    0.000001   -0.763295    4.120000    4.096170    4.102416    4.102416
 spn 2 0    0.000000   -1.128832    4.120000    4.085128    4.102416    4.102416

 Harris energy:
 sumev=       -2.743233  val*vef=     -14.362499   sumtv=      11.619265
 sumec=      -39.640777  cor*vef=    -102.572782   ttcor=      62.932005
 rhoeps=      -9.601995     utot=    -139.945263    ehar=     -74.995989

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:     -7.490648     -6.856648    -14.347297 sumev=   -2.743233   sumtv=   11.604063

 Kohn-Sham energy:
 sumtv=       11.604063  sumtc=        62.932006   ekin=       74.536069
 rhoep=       -9.592116   utot=      -139.940098   ehks=      -74.996145
 mag. mom=     2.000000  (bands)        2.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.241
 mixrho:  sought 2 iter from file mixm; read 0.  RMS DQ=1.35e-2
 mixrho: (warning) scr. and lin-mixed densities had 1625 and 405 negative points
 AMIX: nmix=0 mmix=8  nelts=93876  beta=0.5  tm=5  rmsdel=6.75e-3
 unscreened rms difference:  smooth  0.010357   local  0.036288
   screened rms difference:  smooth  0.010369   local  0.036288   tot  0.013498

 iors  : write restart file (binary, mesh density)

   it  1  of 10    ehf=      -0.001089   ehk=      -0.001245
h zbak=0 mmom=1.9999996 ehf=-.0010887 ehk=-.0012453

 --- BNDFP:  begin iteration 2 of 10 ---

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1    0.167495    0.593754     2.00     6.00

 Smooth charges: Qmesh = 3.406246  Qgauss = 0.593754  core-nuc = -4  tot = 0

 Average es pot at rmt = 0.008355  avg sphere pot = 0.016801   vconst = -0.008355
 Average es pot at MT boundaries after vconst shift
 Site    ves
   1   0.000000  |
 smooth rhoves      2.743064   charge     3.406246
 smvxcm (warning) mesh density negative at 37543 points:  rhomin=-5.54e-6
 smooth rhoeps =   -2.246912 (  -1.725179,  -0.521733)
         rhomu =   -2.931285 (  -2.377669,  -0.553616)
       avg vxc =   -0.191395 (  -0.225050,  -0.157739)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 potential shift to crystal energy zero:    0.000003

 potpus  spin 1 : pnu = 2.914622 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.039326   -1.271736   -1.159329    0.200048      2.8088    0.891233
 1     -0.471549   -1.053722   -0.469241    0.174921     19.1023    0.929745
 2     -0.620206   -0.620206    1.604652    0.368908     16.3477    6.931097
 3     -0.555126   -0.555126    3.518135    0.416961     23.4279   11.934091

 potpus  spin 2 : pnu = 2.908183 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -0.837967   -1.090499   -0.963997    0.211860      2.8184    0.929794
 1     -0.270137   -0.929626   -0.267518    0.186317     19.0733    0.983999
 2     -0.512201   -0.512201    1.742481    0.372406     16.2571    7.102208
 3     -0.454507   -0.454507    3.641633    0.418654     23.3694   12.065273

 3rd gen potential parameters
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.915147    7.418301    0.092820   -0.587641
 1      3.000000    1.000000    -5.887832    7.121418    0.144054   -0.533608
 2      3.000000    1.000000     6.000000   29.630685    0.492412   -0.085938
 3      3.000000    1.000000     9.000000   40.184053    0.569712   -0.056285

 3rd gen potential parameters spin 2
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.110312    7.331756    0.102466   -0.559526
 1      3.000000    1.000000    -5.887832    7.130563    0.153438   -0.500623
 2      3.000000    1.000000     6.000000   29.972887    0.494326   -0.084383
 3      3.000000    1.000000     9.000000   40.378896    0.570595   -0.055849

 Local terms accumulated over augmentation sites:
                  true           sm,loc      difference
 rhoeps:        -9.532961      -2.180625      -7.352336
 rhomu:         -7.417505      -2.305258      -5.112247
 spin2:         -5.126437      -0.540117      -4.586320
 total:        -12.543942      -2.845375      -9.698567
 val*ves       -10.249564      -4.969394      -5.280170
 val*vef       -14.282244      -7.814769      -6.467475
 val chg:        3.697601       3.103847       0.593754
 core chg:       2.000000       2.000000       0.000000
 val mom:        1.803790

 Energy terms:             smooth           local            total
   rhoval*vef             -6.343415        -8.018560       -14.361975
   rhoval*ves             -4.963215        -5.280170       -10.243385
   psnuc*ves              10.449343      -280.095975      -269.646632
   utot                    2.743064      -142.688072      -139.945009
   rho*exc                -2.246912        -7.352336        -9.599248
   rho*vxc                -2.931285        -9.698567       -12.629852
   valence chg             3.406246         0.593754         4.000000
   valence mag             1.838600         0.161400         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.430662,  4 occ states
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.0412 -0.4294 -0.4294 -0.4294  0.1305  0.5276  0.5276  0.5276
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -0.8396 -0.2382 -0.2382 -0.2382  0.2119  0.6260  0.6260  0.6260
 Est Ef = -0.431 < evl(4)=-0.429 ... using qval=4.0, revise to -0.4294

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.431114;   4.000000 electrons;  D(Ef):  913.916
         Sum occ. bands:   -2.7457021  incl. Bloechl correction:   -0.000181
         Mag. moment:       2.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.687098    3.837368   -0.150270      1.796463    2.125046   -0.328583
       contr. to mm extrapolated for r>rmt:   0.163703 est. true mm = 1.960166

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.958012   -1.039323    2.914622    2.914625    2.500000    2.914625
 spn 2 0    0.945317   -0.837839    2.908183    2.908288    2.500000    2.908288
 1     1    1.783756   -0.429712    2.850000    2.889347    2.250000    2.850000
 spn 2 1    0.000000   -1.074202    2.850000    2.172043    2.250000    2.850000
 2     0    0.000012   -0.813133    3.147584    3.132771    3.147584    3.147584
 spn 2 0    0.000000   -1.179998    3.147584    3.108348    3.147584    3.147584
 3     0    0.000001   -0.762432    4.102416    4.096166    4.102416    4.102416
 spn 2 0    0.000000   -1.130141    4.102416    4.085107    4.102416    4.102416

 Harris energy:
 sumev=       -2.745702  val*vef=     -14.361975   sumtv=      11.616273
 sumec=      -39.643187  cor*vef=    -102.575192   ttcor=      62.932005
 rhoeps=      -9.599248     utot=    -139.945009    ehar=     -74.995978

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:     -8.351869     -6.004216    -14.356085 sumev=   -2.745702   sumtv=   11.610383

 Kohn-Sham energy:
 sumtv=       11.610383  sumtc=        62.931909   ekin=       74.542292
 rhoep=       -9.592990   utot=      -139.945439   ehks=      -74.996137
 mag. mom=     2.000000  (bands)        2.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.241
 mixrho:  sought 2 iter from file mixm; read 1.  RMS DQ=6.68e-3  last it=1.35e-2
 mixrho: (warning) scr. and lin-mixed densities had 767 and 611 negative points
 AMIX: nmix=1 mmix=8  nelts=93876  beta=0.5  tm=5  rmsdel=3.34e-3
   tj:-0.97958
 unscreened rms difference:  smooth  0.005191   local  0.017945
   screened rms difference:  smooth  0.005199   local  0.017945   tot  0.006682

 iors  : write restart file (binary, mesh density)

   it  2  of 10    ehf=      -0.001078   ehk=      -0.001237
 From last iter    ehf=      -0.001089   ehk=      -0.001245
 diffe(q)=  0.000010 (0.006682)    tol= 0.000010 (0.000500)   more=T
i zbak=0 mmom=1.9999997 ehf=-.0010784 ehk=-.0012368

 --- BNDFP:  begin iteration 3 of 10 ---

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.040248   -0.142675     2.00     6.00

 Smooth charges: Qmesh = 4.142675  Qgauss = -0.142675  core-nuc = -4  tot = 0

 Average es pot at rmt = 0.009580  avg sphere pot = 0.014106   vconst = -0.009580
 Average es pot at MT boundaries after vconst shift
 Site    ves
   1   0.000000  |
 smooth rhoves      3.512176   charge     4.142675
 smvxcm (warning) mesh density negative at 43201 points:  rhomin=-9.17e-6
 smooth rhoeps =   -3.073362 (  -2.400827,  -0.672535)
         rhomu =   -4.014673 (  -3.315930,  -0.698743)
       avg vxc =   -0.196505 (  -0.232582,  -0.160428)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 potential shift to crystal energy zero:    0.000004

 potpus  spin 1 : pnu = 2.914625 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.039856   -1.271885   -1.159683    0.199881      2.8084    0.890952
 1     -0.471940   -1.052982   -0.469636    0.174783     19.0953    0.929284
 2     -0.619163   -0.619163    1.604911    0.368820     16.3498    6.927384
 3     -0.553868   -0.553868    3.518638    0.416907     23.4296   11.930396

 potpus  spin 2 : pnu = 2.908288 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -0.838825   -1.091272   -0.964882    0.211771      2.8182    0.929555
 1     -0.271026   -0.929804   -0.268410    0.186235     19.0694    0.983736
 2     -0.512261   -0.512261    1.741931    0.372351     16.2583    7.099867
 3     -0.454432   -0.454432    3.641225    0.418620     23.3705   12.062907

 3rd gen potential parameters
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.915548    7.422443    0.092756   -0.587896
 1      3.000000    1.000000    -5.887832    7.123611    0.143940   -0.533940
 2      3.000000    1.000000     6.000000   29.622896    0.492357   -0.085976
 3      3.000000    1.000000     9.000000   40.178294    0.569680   -0.056299

 3rd gen potential parameters spin 2
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.122540    7.332770    0.102355   -0.559707
 1      3.000000    1.000000    -5.887832    7.131787    0.153370   -0.500795
 2      3.000000    1.000000     6.000000   29.968027    0.494292   -0.084406
 3      3.000000    1.000000     9.000000   40.375222    0.570575   -0.055858

 Local terms accumulated over augmentation sites:
                  true           sm,loc      difference
 rhoeps:        -9.529501      -3.006722      -6.522779
 rhomu:         -7.412108      -3.243006      -4.169102
 spin2:         -5.127314      -0.685257      -4.442057
 total:        -12.539422      -3.928263      -8.611160
 val*ves       -10.256927      -4.672847      -5.584080
 val*vef       -14.284990      -8.601110      -5.683880
 val chg:        3.691758       3.834433      -0.142675
 core chg:       2.000000       2.000000       0.000000
 val mom:        1.796538

 Energy terms:             smooth           local            total
   rhoval*vef             -9.097806        -5.266895       -14.364701
   rhoval*ves             -4.666148        -5.584080       -10.250228
   psnuc*ves              11.690500      -281.337573      -269.647073
   utot                    3.512176      -143.460826      -139.948650
   rho*exc                -3.073362        -6.522779        -9.596141
   rho*vxc                -4.014673        -8.611160       -12.625833
   valence chg             4.142675        -0.142675         4.000000
   valence mag             2.323581        -0.323581         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.431114,  4 occ states
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.0417 -0.4297 -0.4297 -0.4297  0.1295  0.5272  0.5272  0.5272
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -0.8405 -0.2389 -0.2389 -0.2389  0.2124  0.6264  0.6264  0.6264
 Est Ef = -0.431 < evl(4)=-0.430 ... using qval=4.0, revise to -0.4297

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.431468;   4.000000 electrons;  D(Ef):  897.926
         Sum occ. bands:   -2.7478315  incl. Bloechl correction:   -0.000184
         Mag. moment:       2.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.687332    3.841355   -0.154024      1.796391    2.115318   -0.318927
       contr. to mm extrapolated for r>rmt:   0.163651 est. true mm = 1.960042

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.958114   -1.039806    2.914625    2.914668    2.500000    2.914668
 spn 2 0    0.945470   -0.838719    2.908288    2.908375    2.500000    2.908375
 1     1    1.783735   -0.430247    2.850000    2.889281    2.250000    2.850000
 spn 2 1    0.000000   -1.072414    2.850000    2.172388    2.250000    2.850000
 2     0    0.000012   -0.811998    3.147584    3.132774    3.147584    3.147584
 spn 2 0    0.000000   -1.180241    3.147584    3.108338    3.147584    3.147584
 3     0    0.000001   -0.761039    4.102416    4.096169    4.102416    4.102416
 spn 2 0    0.000000   -1.130240    4.102416    4.085103    4.102416    4.102416

 Harris energy:
 sumev=       -2.747832  val*vef=     -14.364701   sumtv=      11.616870
 sumec=      -39.644412  cor*vef=    -102.576369   ttcor=      62.931957
 rhoeps=      -9.596141     utot=    -139.948650    ehar=     -74.995965

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:     -9.136284     -5.227344    -14.363629 sumev=   -2.747832   sumtv=   11.615797

 Kohn-Sham energy:
 sumtv=       11.615797  sumtc=        62.931331   ekin=       74.547128
 rhoep=       -9.593750   utot=      -139.949504   ehks=      -74.996126
 mag. mom=     2.000000  (bands)        2.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.241
 mixrho:  sought 2 iter from file mixm; read 2.  RMS DQ=1.20e-4  last it=6.68e-3
 mixrho: (warning) scr. and lin-mixed densities had 0 and 177 negative points
 AMIX: nmix=2 mmix=8  nelts=93876  beta=0.5  tm=5  rmsdel=6.01e-5
   tj:-1.27791   0.62413
 unscreened rms difference:  smooth  0.000120   local  0.000266
   screened rms difference:  smooth  0.000121   local  0.000266   tot  0.000120

 iors  : write restart file (binary, mesh density)

   it  3  of 10    ehf=      -0.001065   ehk=      -0.001226
 From last iter    ehf=      -0.001078   ehk=      -0.001237
 diffe(q)=  0.000014 (0.000120)    tol= 0.000010 (0.000500)   more=T
i zbak=0 mmom=1.9999997 ehf=-.0010647 ehk=-.0012256

 --- BNDFP:  begin iteration 4 of 10 ---

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.044605   -0.158122     2.00     6.00

 Smooth charges: Qmesh = 4.158122  Qgauss = -0.158122  core-nuc = -4  tot = 0

 Average es pot at rmt = 0.009023  avg sphere pot = 0.014073   vconst = -0.009023
 Average es pot at MT boundaries after vconst shift
 Site    ves
   1   0.000000  |
 smooth rhoves      3.521959   charge     4.158122
 smvxcm (warning) mesh density negative at 31975 points:  rhomin=-6.89e-6
 smooth rhoeps =   -3.090050 (  -2.410167,  -0.679883)
         rhomu =   -4.036430 (  -3.329178,  -0.707251)
       avg vxc =   -0.199977 (  -0.235238,  -0.164716)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 potential shift to crystal energy zero:    0.000005

 potpus  spin 1 : pnu = 2.914668 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.039937   -1.271942   -1.159778    0.199849      2.8083    0.890876
 1     -0.472033   -1.052824   -0.469730    0.174755     19.0933    0.929208
 2     -0.618938   -0.618938    1.604946    0.368798     16.3503    6.926467
 3     -0.553591   -0.553591    3.518736    0.416894     23.4301   11.929486

 potpus  spin 2 : pnu = 2.908375 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -0.838796   -1.091267   -0.964918    0.211736      2.8182    0.929431
 1     -0.271055   -0.929520   -0.268440    0.186201     19.0673    0.983642
 2     -0.511927   -0.511927    1.742037    0.372325     16.2589    7.098767
 3     -0.454035   -0.454035    3.641400    0.418604     23.3710   12.061794

 3rd gen potential parameters
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.921398    7.422919    0.092714   -0.587963
 1      3.000000    1.000000    -5.887832    7.124247    0.143917   -0.534001
 2      3.000000    1.000000     6.000000   29.620957    0.492344   -0.085985
 3      3.000000    1.000000     9.000000   40.176877    0.569673   -0.056302

 3rd gen potential parameters spin 2
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.132650    7.332748    0.102279   -0.559798
 1      3.000000    1.000000    -5.887832    7.132446    0.153342   -0.500862
 2      3.000000    1.000000     6.000000   29.965729    0.494276   -0.084417
 3      3.000000    1.000000     9.000000   40.373494    0.570566   -0.055862

 Local terms accumulated over augmentation sites:
                  true           sm,loc      difference
 rhoeps:        -9.528670      -3.023394      -6.505276
 rhomu:         -7.411560      -3.256203      -4.155357
 spin2:         -5.126790      -0.693848      -4.432942
 total:        -12.538350      -3.950051      -8.588299
 val*ves       -10.260092      -4.668999      -5.591093
 val*vef       -14.287048      -8.619050      -5.667997
 val chg:        3.689760       3.847881      -0.158122
 core chg:       2.000000       2.000000       0.000000
 val mom:        1.796357

 Energy terms:             smooth           local            total
   rhoval*vef             -9.161277        -5.205254       -14.366531
   rhoval*ves             -4.662104        -5.591093       -10.253197
   psnuc*ves              11.706021      -281.354309      -269.648288
   utot                    3.521959      -143.472701      -139.950743
   rho*exc                -3.090050        -6.505276        -9.595326
   rho*vxc                -4.036430        -8.588299       -12.624728
   valence chg             4.158122        -0.158122         4.000000
   valence mag             2.323725        -0.323725         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.431468,  4 occ states
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.0418 -0.4298 -0.4298 -0.4298  0.1293  0.5273  0.5273  0.5273
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -0.8405 -0.2389 -0.2389 -0.2389  0.2111  0.6260  0.6260  0.6260
 Est Ef = -0.431 < evl(4)=-0.430 ... using qval=4.0, revise to -0.4298

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.431533;   4.000000 electrons;  D(Ef):  888.912
         Sum occ. bands:   -2.7481615  incl. Bloechl correction:   -0.000185
         Mag. moment:       2.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.687306    3.841270   -0.153964      1.796331    2.114857   -0.318526
       contr. to mm extrapolated for r>rmt:   0.163654 est. true mm = 1.959986

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.958146   -1.039915    2.914668    2.914687    2.500000    2.914687
 spn 2 0    0.945487   -0.838794    2.908375    2.908376    2.500000    2.908376
 1     1    1.783660   -0.430386    2.850000    2.889254    2.250000    2.850000
 spn 2 1    0.000000   -1.086680    2.850000    2.170035    2.250000    2.850000
 2     0    0.000012   -0.811746    3.147584    3.132775    3.147584    3.147584
 spn 2 0    0.000000   -1.179678    3.147584    3.108345    3.147584    3.147584
 3     0    0.000001   -0.760726    4.102416    4.096170    4.102416    4.102416
 spn 2 0    0.000000   -1.130167    4.102416    4.085096    4.102416    4.102416

 Harris energy:
 sumev=       -2.748161  val*vef=     -14.366531   sumtv=      11.618369
 sumec=      -39.644144  cor*vef=    -102.575788   ttcor=      62.931644
 rhoeps=      -9.595326     utot=    -139.950743    ehar=     -74.996055

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:     -9.150390     -5.214178    -14.364568 sumev=   -2.748161   sumtv=   11.616406

 Kohn-Sham energy:
 sumtv=       11.616406  sumtc=        62.931049   ekin=       74.547455
 rhoep=       -9.593800   utot=      -139.949778   ehks=      -74.996123
 mag. mom=     2.000000  (bands)        2.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.241
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=4.77e-5  last it=1.20e-4
 AMIX: nmix=2 mmix=8  nelts=93876  beta=0.5  tm=5  rmsdel=2.38e-5
   tj:-0.30542   0.01065
 unscreened rms difference:  smooth  0.000031   local  0.000146
   screened rms difference:  smooth  0.000027   local  0.000146   tot  0.000048

 iors  : write restart file (binary, mesh density)

   it  4  of 10    ehf=      -0.001155   ehk=      -0.001223
 From last iter    ehf=      -0.001065   ehk=      -0.001226
 diffe(q)= -0.000090 (0.000048)    tol= 0.000010 (0.000500)   more=T
i zbak=0 mmom=1.9999997 ehf=-.0011549 ehk=-.0012228

 --- BNDFP:  begin iteration 5 of 10 ---

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.043547   -0.154369     2.00     6.00

 Smooth charges: Qmesh = 4.154369  Qgauss = -0.154369  core-nuc = -4  tot = 0

 Average es pot at rmt = 0.008695  avg sphere pot = 0.014098   vconst = -0.008695
 Average es pot at MT boundaries after vconst shift
 Site    ves
   1   0.000000  |
 smooth rhoves      3.514664   charge     4.154369
 smvxcm (warning) mesh density negative at 21263 points:  rhomin=-4.85e-6
 smooth rhoeps =   -3.085149 (  -2.405083,  -0.680066)
         rhomu =   -4.029966 (  -3.322197,  -0.707769)
       avg vxc =   -0.202142 (  -0.236393,  -0.167892)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 potential shift to crystal energy zero:    0.000005

 potpus  spin 1 : pnu = 2.914687 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.040044   -1.272024   -1.159884    0.199829      2.8083    0.890825
 1     -0.472141   -1.052781   -0.469838    0.174736     19.0924    0.929146
 2     -0.618855   -0.618855    1.604928    0.368787     16.3506    6.925971
 3     -0.553481   -0.553481    3.518752    0.416887     23.4303   11.929008

 potpus  spin 2 : pnu = 2.908376 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -0.838889   -1.091293   -0.964981    0.211707      2.8182    0.929378
 1     -0.271122   -0.929387   -0.268508    0.186177     19.0665    0.983556
 2     -0.511761   -0.511761    1.742076    0.372311     16.2593    7.098140
 3     -0.453835   -0.453835    3.641479    0.418595     23.3713   12.061180

 3rd gen potential parameters
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.923956    7.423232    0.092693   -0.588006
 1      3.000000    1.000000    -5.887832    7.124533    0.143901   -0.534046
 2      3.000000    1.000000     6.000000   29.619916    0.492337   -0.085990
 3      3.000000    1.000000     9.000000   40.176136    0.569669   -0.056304

 3rd gen potential parameters spin 2
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.132836    7.333368    0.102267   -0.559841
 1      3.000000    1.000000    -5.887832    7.132720    0.153323   -0.500916
 2      3.000000    1.000000     6.000000   29.964429    0.494267   -0.084423
 3      3.000000    1.000000     9.000000   40.372544    0.570561   -0.055864

 Local terms accumulated over augmentation sites:
                  true           sm,loc      difference
 rhoeps:        -9.528144      -3.018492      -6.509652
 rhomu:         -7.411204      -3.249201      -4.162002
 spin2:         -5.126465      -0.694408      -4.432057
 total:        -12.537669      -3.943609      -8.594059
 val*ves       -10.261256      -4.673186      -5.588069
 val*vef       -14.287538      -8.616796      -5.670742
 val chg:        3.688744       3.843112      -0.154369
 core chg:       2.000000       2.000000       0.000000
 val mom:        1.796348

 Energy terms:             smooth           local            total
   rhoval*vef             -9.147511        -5.219386       -14.366898
   rhoval*ves             -4.666190        -5.588069       -10.254259
   psnuc*ves              11.695519      -281.343238      -269.647719
   utot                    3.514664      -143.465654      -139.950989
   rho*exc                -3.085149        -6.509652        -9.594801
   rho*vxc                -4.029966        -8.594059       -12.624025
   valence chg             4.154369        -0.154369         4.000000
   valence mag             2.318556        -0.318557         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.431533,  4 occ states
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.0419 -0.4299 -0.4299 -0.4299  0.1292  0.5274  0.5274  0.5274
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -0.8406 -0.2389 -0.2389 -0.2389  0.2102  0.6257  0.6257  0.6257
 Est Ef = -0.432 < evl(4)=-0.430 ... using qval=4.0, revise to -0.4299

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.431621;   4.000000 electrons;  D(Ef):  885.518
         Sum occ. bands:   -2.7485461  incl. Bloechl correction:   -0.000185
         Mag. moment:       2.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.687373    3.842474   -0.155101      1.796374    2.115763   -0.319389
       contr. to mm extrapolated for r>rmt:   0.163601 est. true mm = 1.959976

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.958168   -1.040029    2.914687    2.914700    2.500000    2.914700
 spn 2 0    0.945499   -0.838887    2.908376    2.908377    2.500000    2.908377
 1     1    1.783693   -0.430503    2.850000    2.889254    2.250000    2.850000
 spn 2 1    0.000000   -1.095688    2.850000    2.168600    2.250000    2.850000
 2     0    0.000012   -0.811661    3.147584    3.132775    3.147584    3.147584
 spn 2 0    0.000000   -1.179351    3.147584    3.108351    3.147584    3.147584
 3     0    0.000001   -0.760622    4.102416    4.096170    4.102416    4.102416
 spn 2 0    0.000000   -1.130205    4.102416    4.085091    4.102416    4.102416

 Harris energy:
 sumev=       -2.748546  val*vef=     -14.366898   sumtv=      11.618352
 sumec=      -39.644381  cor*vef=    -102.575727   ttcor=      62.931347
 rhoeps=      -9.594801     utot=    -139.950989    ehar=     -74.996092

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:     -9.150774     -5.215194    -14.365968 sumev=   -2.748546   sumtv=   11.617422

 Kohn-Sham energy:
 sumtv=       11.617422  sumtc=        62.931001   ekin=       74.548423
 rhoep=       -9.593961   utot=      -139.950583   ehks=      -74.996121
 mag. mom=     2.000000  (bands)        2.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.241
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.91e-5  last it=4.77e-5
 AMIX: nmix=2 mmix=8  nelts=93876  beta=0.5  tm=5  rmsdel=9.54e-6
   tj: 0.18832  -0.06251
 unscreened rms difference:  smooth  0.000016   local  0.000035
   screened rms difference:  smooth  0.000016   local  0.000035   tot  0.000019

 iors  : write restart file (binary, mesh density)

   it  5  of 10    ehf=      -0.001192   ehk=      -0.001221
 From last iter    ehf=      -0.001155   ehk=      -0.001223
 diffe(q)= -0.000037 (0.000019)    tol= 0.000010 (0.000500)   more=T
i zbak=0 mmom=1.9999997 ehf=-.0011922 ehk=-.001221

 --- BNDFP:  begin iteration 6 of 10 ---

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.043832   -0.155381     2.00     6.00

 Smooth charges: Qmesh = 4.155381  Qgauss = -0.155381  core-nuc = -4  tot = 0

 Average es pot at rmt = 0.008587  avg sphere pot = 0.014097   vconst = -0.008587
 Average es pot at MT boundaries after vconst shift
 Site    ves
   1   0.000000  |
 smooth rhoves      3.514835   charge     4.155381
 smvxcm (warning) mesh density negative at 16055 points:  rhomin=-3.86e-6
 smooth rhoeps =   -3.086320 (  -2.406113,  -0.680206)
         rhomu =   -4.031499 (  -3.323636,  -0.707863)
       avg vxc =   -0.203044 (  -0.236797,  -0.169291)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 potential shift to crystal energy zero:    0.000005

 potpus  spin 1 : pnu = 2.914700 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.040069   -1.272049   -1.159916    0.199822      2.8083    0.890804
 1     -0.472174   -1.052762   -0.469871    0.174730     19.0921    0.929127
 2     -0.618822   -0.618822    1.604924    0.368782     16.3507    6.925797
 3     -0.553439   -0.553439    3.518761    0.416885     23.4304   11.928837

 potpus  spin 2 : pnu = 2.908377 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -0.838914   -1.091297   -0.964998    0.211698      2.8182    0.929361
 1     -0.271140   -0.929340   -0.268526    0.186169     19.0661    0.983529
 2     -0.511702   -0.511702    1.742091    0.372306     16.2594    7.097931
 3     -0.453765   -0.453765    3.641508    0.418592     23.3713   12.060972

 3rd gen potential parameters
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.925673    7.423266    0.092681   -0.588023
 1      3.000000    1.000000    -5.887832    7.124640    0.143896   -0.534061
 2      3.000000    1.000000     6.000000   29.619550    0.492334   -0.085992
 3      3.000000    1.000000     9.000000   40.175870    0.569668   -0.056304

 3rd gen potential parameters spin 2
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.132976    7.333568    0.102262   -0.559855
 1      3.000000    1.000000    -5.887832    7.132818    0.153316   -0.500933
 2      3.000000    1.000000     6.000000   29.963995    0.494264   -0.084425
 3      3.000000    1.000000     9.000000   40.372222    0.570559   -0.055865

 Local terms accumulated over augmentation sites:
                  true           sm,loc      difference
 rhoeps:        -9.527969      -3.019669      -6.508300
 rhomu:         -7.411094      -3.250641      -4.160453
 spin2:         -5.126348      -0.694515      -4.431832
 total:        -12.537441      -3.945156      -8.592285
 val*ves       -10.261657      -4.673290      -5.588366
 val*vef       -14.287719      -8.618447      -5.669272
 val chg:        3.688397       3.843778      -0.155381
 core chg:       2.000000       2.000000       0.000000
 val mom:        1.796351

 Energy terms:             smooth           local            total
   rhoval*vef             -9.152087        -5.214942       -14.367029
   rhoval*ves             -4.666258        -5.588366       -10.254625
   psnuc*ves              11.695928      -281.343361      -269.647432
   utot                    3.514835      -143.465863      -139.951028
   rho*exc                -3.086320        -6.508300        -9.594620
   rho*vxc                -4.031499        -8.592285       -12.623784
   valence chg             4.155381        -0.155381         4.000000
   valence mag             2.319235        -0.319236         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.431621,  4 occ states
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.0420 -0.4299 -0.4299 -0.4299  0.1292  0.5274  0.5274  0.5274
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -0.8407 -0.2390 -0.2390 -0.2390  0.2098  0.6256  0.6256  0.6256
 Est Ef = -0.432 < evl(4)=-0.430 ... using qval=4.0, revise to -0.4299

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.431646;   4.000000 electrons;  D(Ef):  884.401
         Sum occ. bands:   -2.7486572  incl. Bloechl correction:   -0.000185
         Mag. moment:       2.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.687399    3.842977   -0.155578      1.796394    2.116182   -0.319788
       contr. to mm extrapolated for r>rmt:   0.163579 est. true mm = 1.959973

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.958176   -1.040063    2.914700    2.914705    2.500000    2.914705
 spn 2 0    0.945502   -0.838914    2.908377    2.908377    2.500000    2.908377
 1     1    1.783708   -0.430537    2.850000    2.889255    2.250000    2.850000
 spn 2 1    0.000000   -1.099151    2.850000    2.168057    2.250000    2.850000
 2     0    0.000012   -0.811628    3.147584    3.132775    3.147584    3.147584
 spn 2 0    0.000000   -1.179226    3.147584    3.108353    3.147584    3.147584
 3     0    0.000001   -0.760583    4.102416    4.096170    4.102416    4.102416
 spn 2 0    0.000000   -1.130235    4.102416    4.085089    4.102416    4.102416

 Harris energy:
 sumev=       -2.748657  val*vef=     -14.367029   sumtv=      11.618372
 sumec=      -39.644447  cor*vef=    -102.575621   ttcor=      62.931174
 rhoeps=      -9.594620     utot=    -139.951028    ehar=     -74.996103

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:     -9.153584     -5.212840    -14.366424 sumev=   -2.748657   sumtv=   11.617766

 Kohn-Sham energy:
 sumtv=       11.617766  sumtc=        62.930991   ekin=       74.548758
 rhoep=       -9.594018   utot=      -139.950860   ehks=      -74.996120
 mag. mom=     2.000000  (bands)        2.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.241
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=1.29e-5  last it=1.91e-5
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=93876  beta=0.5  tm=5  rmsdel=6.47e-6
   tj:-1.73999
 unscreened rms difference:  smooth  0.000011   local  0.000024
   screened rms difference:  smooth  0.000010   local  0.000024   tot  0.000013

 iors  : write restart file (binary, mesh density)

   it  6  of 10    ehf=      -0.001203   ehk=      -0.001220
 From last iter    ehf=      -0.001192   ehk=      -0.001221
 diffe(q)= -0.000010 (0.000013)    tol= 0.000010 (0.000500)   more=T
i zbak=0 mmom=1.9999997 ehf=-.0012026 ehk=-.0012204

 --- BNDFP:  begin iteration 7 of 10 ---

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.044225   -0.156774     2.00     6.00

 Smooth charges: Qmesh = 4.156774  Qgauss = -0.156774  core-nuc = -4  tot = 0

 Average es pot at rmt = 0.008364  avg sphere pot = 0.014099   vconst = -0.008364
 Average es pot at MT boundaries after vconst shift
 Site    ves
   1   0.000000  |
 smooth rhoves      3.514374   charge     4.156774
 smvxcm (warning) mesh density negative at 5903 points:  rhomin=-1.57e-6
 smooth rhoeps =   -3.087990 (  -2.407756,  -0.680234)
         rhomu =   -4.033695 (  -3.325933,  -0.707761)
       avg vxc =   -0.205063 (  -0.237547,  -0.172580)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 potential shift to crystal energy zero:    0.000005

 potpus  spin 1 : pnu = 2.914705 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.040132   -1.272086   -1.159970    0.199809      2.8082    0.890780
 1     -0.472229   -1.052719   -0.469926    0.174718     19.0913    0.929092
 2     -0.618753   -0.618753    1.604922    0.368774     16.3509    6.925455
 3     -0.553350   -0.553350    3.518781    0.416880     23.4305   11.928496

 potpus  spin 2 : pnu = 2.908377 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -0.838953   -1.091294   -0.965018    0.211680      2.8181    0.929333
 1     -0.271162   -0.929237   -0.268548    0.186155     19.0655    0.983482
 2     -0.511578   -0.511578    1.742130    0.372297     16.2596    7.097521
 3     -0.453618   -0.453618    3.641572    0.418586     23.3715   12.060560

 3rd gen potential parameters
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.926298    7.423591    0.092673   -0.588045
 1      3.000000    1.000000    -5.887832    7.124862    0.143887   -0.534087
 2      3.000000    1.000000     6.000000   29.618829    0.492329   -0.085996
 3      3.000000    1.000000     9.000000   40.175339    0.569665   -0.056306

 3rd gen potential parameters spin 2
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.132983    7.333992    0.102255   -0.559879
 1      3.000000    1.000000    -5.887832    7.133028    0.153304   -0.500964
 2      3.000000    1.000000     6.000000   29.963143    0.494258   -0.084429
 3      3.000000    1.000000     9.000000   40.371582    0.570556   -0.055866

 Local terms accumulated over augmentation sites:
                  true           sm,loc      difference
 rhoeps:        -9.527653      -3.021360      -6.506294
 rhomu:         -7.410921      -3.252949      -4.157972
 spin2:         -5.126113      -0.694438      -4.431674
 total:        -12.537034      -3.947388      -8.589646
 val*ves       -10.262492      -4.673854      -5.588637
 val*vef       -14.288141      -8.621242      -5.666899
 val chg:        3.687710       3.844484      -0.156774
 core chg:       2.000000       2.000000       0.000000
 val mom:        1.796393

 Energy terms:             smooth           local            total
   rhoval*vef             -9.158832        -5.208512       -14.367344
   rhoval*ves             -4.666750        -5.588637       -10.255387
   psnuc*ves              11.695498      -281.342967      -269.647469
   utot                    3.514374      -143.465802      -139.951428
   rho*exc                -3.087990        -6.506294        -9.594283
   rho*vxc                -4.033695        -8.589646       -12.623341
   valence chg             4.156774        -0.156774         4.000000
   valence mag             2.320450        -0.320450         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.431646,  4 occ states
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.0420 -0.4299 -0.4299 -0.4299  0.1292  0.5275  0.5275  0.5275
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -0.8407 -0.2390 -0.2390 -0.2390  0.2090  0.6254  0.6254  0.6254
 Est Ef = -0.432 < evl(4)=-0.430 ... using qval=4.0, revise to -0.4299

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.431685;   4.000000 electrons;  D(Ef):  882.362
         Sum occ. bands:   -2.7488310  incl. Bloechl correction:   -0.000186
         Mag. moment:       2.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.687456    3.844178   -0.156722      1.796439    2.117171   -0.320732
       contr. to mm extrapolated for r>rmt:   0.163531 est. true mm = 1.959971

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.958191   -1.040121    2.914705    2.914714    2.500000    2.914714
 spn 2 0    0.945508   -0.838953    2.908377    2.908377    2.500000    2.908377
 1     1    1.783745   -0.430593    2.850000    2.889258    2.250000    2.850000
 spn 2 1    0.000000   -1.106108    2.850000    2.166980    2.250000    2.850000
 2     0    0.000012   -0.811560    3.147584    3.132775    3.147584    3.147584
 spn 2 0    0.000000   -1.178965    3.147584    3.108358    3.147584    3.147584
 3     0    0.000001   -0.760498    4.102416    4.096170    4.102416    4.102416
 spn 2 0    0.000000   -1.130308    4.102416    4.085084    4.102416    4.102416

 Harris energy:
 sumev=       -2.748831  val*vef=     -14.367344   sumtv=      11.618513
 sumec=      -39.644494  cor*vef=    -102.575576   ttcor=      62.931083
 rhoeps=      -9.594283     utot=    -139.951428    ehar=     -74.996116

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:     -9.159086     -5.208171    -14.367257 sumev=   -2.748831   sumtv=   11.618426

 Kohn-Sham energy:
 sumtv=       11.618426  sumtc=        62.930934   ekin=       74.549359
 rhoep=       -9.594131   utot=      -139.951348   ehks=      -74.996119
 mag. mom=     2.000000  (bands)        2.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.241
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=4.35e-6  last it=1.29e-5
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=93876  beta=0.5  tm=5  rmsdel=2.17e-6
   tj:-0.43211
 unscreened rms difference:  smooth  0.000003   local  0.000009
   screened rms difference:  smooth  0.000003   local  0.000009   tot  0.000004

 iors  : write restart file (binary, mesh density)

   it  7  of 10    ehf=      -0.001216   ehk=      -0.001219
 From last iter    ehf=      -0.001203   ehk=      -0.001220
 diffe(q)= -0.000013 (0.000004)    tol= 0.000010 (0.000500)   more=T
i zbak=0 mmom=1.9999997 ehf=-.0012161 ehk=-.0012192

 --- BNDFP:  begin iteration 8 of 10 ---

 rhomom:   ib   ilm      qmom        Qval       Qc        Z
            1     1   -0.044372   -0.157296     2.00     6.00

 Smooth charges: Qmesh = 4.157296  Qgauss = -0.157296  core-nuc = -4  tot = 0

 Average es pot at rmt = 0.008268  avg sphere pot = 0.014100   vconst = -0.008268
 Average es pot at MT boundaries after vconst shift
 Site    ves
   1   0.000000  |
 smooth rhoves      3.514050   charge     4.157296
 smvxcm (warning) mesh density negative at 1403 points:  rhomin=-3.62e-7
 smooth rhoeps =   -3.088663 (  -2.408490,  -0.680173)
         rhomu =   -4.034584 (  -3.326958,  -0.707626)
       avg vxc =   -0.206004 (  -0.237841,  -0.174167)

 locpot:

 site  1  z=  6.0  rmt= 3.00000  nr=369   a=0.020  nlml=16  rg=0.750  Vfloat=T
 potential shift to crystal energy zero:    0.000005

 potpus  spin 1 : pnu = 2.914714 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -1.040142   -1.272096   -1.159986    0.199804      2.8082    0.890767
 1     -0.472244   -1.052697   -0.469942    0.174714     19.0910    0.929081
 2     -0.618721   -0.618721    1.604925    0.368771     16.3509    6.925318
 3     -0.553310   -0.553310    3.518793    0.416878     23.4306   11.928357

 potpus  spin 2 : pnu = 2.908377 2.850000 3.147584 4.102416
 l        enu         v           c          srdel        qpar        ppar
 0     -0.838959   -1.091284   -0.965017    0.211674      2.8181    0.929324
 1     -0.271161   -0.929189   -0.268547    0.186150     19.0652    0.983467
 2     -0.511523   -0.511523    1.742151    0.372293     16.2597    7.097360
 3     -0.453553   -0.453553    3.641603    0.418583     23.3716   12.060394

 3rd gen potential parameters
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.927565    7.423627    0.092665   -0.588056
 1      3.000000    1.000000    -5.887832    7.124959    0.143883   -0.534096
 2      3.000000    1.000000     6.000000   29.618539    0.492327   -0.085997
 3      3.000000    1.000000     9.000000   40.175121    0.569664   -0.056306

 3rd gen potential parameters spin 2
 l          a      <phi phi>       a*D        a*Ddot      phi(a)      phip(a)
 0      3.000000    1.000000   -10.132991    7.334161    0.102253   -0.559887
 1      3.000000    1.000000    -5.887832    7.133122    0.153300   -0.500974
 2      3.000000    1.000000     6.000000   29.962806    0.494256   -0.084431
 3      3.000000    1.000000     9.000000   40.371324    0.570555   -0.055867

 Local terms accumulated over augmentation sites:
                  true           sm,loc      difference
 rhoeps:        -9.527539      -3.022049      -6.505490
 rhomu:         -7.410879      -3.253986      -4.156893
 spin2:         -5.126008      -0.694312      -4.431695
 total:        -12.536887      -3.948299      -8.588588
 val*ves       -10.262842      -4.674117      -5.588724
 val*vef       -14.288343      -8.622416      -5.665927
 val chg:        3.687435       3.844731      -0.157296
 core chg:       2.000000       2.000000       0.000000
 val mom:        1.796435

 Energy terms:             smooth           local            total
   rhoval*vef             -9.161464        -5.206030       -14.367494
   rhoval*ves             -4.666983        -5.588724       -10.255707
   psnuc*ves              11.695084      -281.342575      -269.647491
   utot                    3.514050      -143.465649      -139.951599
   rho*exc                -3.088663        -6.505490        -9.594153
   rho*vxc                -4.034584        -8.588588       -12.623172
   valence chg             4.157296        -0.157296         4.000000
   valence mag             2.321057        -0.321057         2.000000
   core charge             2.000000         0.000000         2.000000

 Charges:  valence     4.00000   cores     2.00000   nucleii    -6.00000
    hom background     0.00000   deviation from neutrality:      0.00000

 Read qp weights ...  ef=-0.431685,  4 occ states
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -1.0420 -0.4299 -0.4299 -0.4299  0.1292  0.5276  0.5276  0.5276
 bndfp:  kpt 1 of 8, k=  0.00000  0.00000  0.00000
 -0.8407 -0.2390 -0.2390 -0.2390  0.2087  0.6253  0.6253  0.6253
 Est Ef = -0.432 < evl(4)=-0.430 ... using qval=4.0, revise to -0.4299

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.431693;   4.000000 electrons;  D(Ef):  881.578
         Sum occ. bands:   -2.7488687  incl. Bloechl correction:   -0.000186
         Mag. moment:       2.000000

 Saved qp weights ...

 mkrout:  Qtrue      sm,loc       local        true mm   smooth mm    local mm
   1    3.687482    3.844790   -0.157308      1.796460    2.117601   -0.321141
       contr. to mm extrapolated for r>rmt:   0.163511 est. true mm = 1.959971

 Symmetrize density..

 Make new boundary conditions for phi,phidot..

 site    1   species   1:C       
 l  idmod     ql         ebar        pold        ptry        pfree        pnew
 0     0    0.958197   -1.040137    2.914714    2.914718    2.500000    2.914718
 spn 2 0    0.945511   -0.838958    2.908377    2.908378    2.500000    2.908378
 1     1    1.783762   -0.430608    2.850000    2.889261    2.250000    2.850000
 spn 2 1    0.000000   -1.108727    2.850000    2.166579    2.250000    2.850000
 2     0    0.000012   -0.811529    3.147584    3.132774    3.147584    3.147584
 spn 2 0    0.000000   -1.178864    3.147584    3.108360    3.147584    3.147584
 3     0    0.000001   -0.760462    4.102416    4.096169    4.102416    4.102416
 spn 2 0    0.000000   -1.130336    4.102416    4.085082    4.102416    4.102416

 Harris energy:
 sumev=       -2.748869  val*vef=     -14.367494   sumtv=      11.618626
 sumec=      -39.644478  cor*vef=    -102.575486   ttcor=      62.931008
 rhoeps=      -9.594153     utot=    -139.951599    ehar=     -74.996119

 Harris correction to forces: screened shift in core+nuclear density  
  ib         delta-n dVes             delta-n dVxc               total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 shift forces to make zero average correction:            0.00    0.00    0.00

 srhov:     -9.161592     -5.205950    -14.367542 sumev=   -2.748869   sumtv=   11.618674

 Kohn-Sham energy:
 sumtv=       11.618674  sumtc=        62.930910   ekin=       74.549584
 rhoep=       -9.594176   utot=      -139.951527   ehks=      -74.996119
 mag. mom=     2.000000  (bands)        2.000000  (output rho)

Forces:
  ib           estatic                  eigval                    total
   1    0.00    0.00    0.00     0.00    0.00    0.00     0.00    0.00    0.00
 Maximum Harris force = 0 mRy/au (site 1)

 Symmetrize forces ...
  
 mixing: mode=A  nmix=2  beta=.5  elind=.241
 mixrho:  sought 2 iter from file mixm; read 3.  RMS DQ=9.32e-7  last it=4.35e-6
 AMIX: condition of normal eqns >100000. Reducing nmix to 1
 AMIX: nmix=1 mmix=8  nelts=93876  beta=0.5  tm=5  rmsdel=4.66e-7
   tj:-0.08439
 unscreened rms difference:  smooth  0.000001   local  0.000002
   screened rms difference:  smooth  0.000000   local  0.000002   tot  0.000001

 iors  : write restart file (binary, mesh density)

   it  8  of 10    ehf=      -0.001219   ehk=      -0.001219
 From last iter    ehf=      -0.001216   ehk=      -0.001219
 diffe(q)= -0.000002 (0.000001)    tol= 0.000010 (0.000500)   more=F
c zbak=0 mmom=1.9999997 ehf=-.0012185 ehk=-.0012188
 Exit 0 LMF 
 CPU time:  30.771s   Wall clock 30.805s  at  10:10:44 10.08.2018  on  localhost.localdomai
