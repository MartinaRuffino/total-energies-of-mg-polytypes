rdcmd:  lmfdmft ni -vnk=12 --rs=1,0 -vnsp=2 -vtpd1=4 -vbigbas=3 -vehmax=-.4 -vrwa=0.3359 -vemaxs=2.5 -vconvc=1e-6 -vbeta=.3 --ldadc=82.2 -job=1 -vindmfl=f
 -------------  START LMFDMFT - ref: 4e6aff7 --------------
#rf  171: (file atparms) element Ni : Z1=28 eref=0 mom=.6 rsm=1.3 rsmd=1
 HEADER Master file for multinary TM compounds, suitable for GW

 LMFDMFT:  nbas = 1  nspec = 1  verb 31,20
 special:: forces
 pot:      spin-pol, XC:BH, read Sigma
 float:    float P v6-style, v6-ebar
 autoread: none
 bz:       metal(5), tetra 
 dmft:     1 blocks  (1 independent)

                Plat                                  Qlat
   0.000000   0.500000   0.500000       -1.000000   1.000000   1.000000
   0.500000   0.000000   0.500000        1.000000  -1.000000   1.000000
   0.500000   0.500000   0.000000        1.000000   1.000000  -1.000000
   alat = 6.64  Cell vol = 73.188736

 LATTC:  as= 2.000   tol=1.00e-8   alat= 6.64000   awald= 0.478   lmax=6
         r1=  1.965   nkd= 135      q1=  5.929   nkg= 229

 SGROUP: 1 symmetry operations from 0 generators
 SYMLAT: Bravais system is cubic with 48 symmetry operations.
 SYMCRY: crystal invariant under 48 symmetry operations for tol=1e-5
 GROUPG: the following are sufficient to generate the space group:
         i*r3(1,1,-1) r4x
         i*r3(1,1,-1) r4x
 MKSYM:  found 48 space group operations
         no attempt to add inversion symmetry
 BZMESH: 72 irreducible QP from 1728 ( 12 12 12 )  shift= F F F
 TETIRR: sorting 10368 tetrahedra ... 247 inequivalent ones found

 species data:  augmentation                           density
 spec       rmt   rsma lmxa kmxa      lmxl     rg   rsmv  kmxv foca   rfoca
 Ni       2.230  0.892    4    4         4  0.558  1.115    15    1   0.892

 MSHSIZ: mesh has 11 x 11 x 11 divisions; length 0.427, 0.427, 0.427
         generated from gmax = 9.0 a.u. : 893 vectors of 1331 (67%)

 GVLIST: gmax = 9 a.u. created 893 vectors of 1331 (67%)
         mesh has 11 x 11 x 11 divisions; length 0.427, 0.427, 0.427
 SGVSYM: 40 symmetry stars found for 893 reciprocal lattice vectors
 

 Makidx:  hamiltonian dimensions Low, Int, High, Negl: 39 0 16 20
 kappa   Low   Int   High  L+I  L+I+H  Neglected
   1      25     0     0    25    25       0
   2       9     0    16     9    25       0
   3       5     0     0     5     5      20
  all     39     0    16    39    55      20
 suham :  25 augmentation channels, 25 local potential channels  Maximum lmxa=4

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol=1.0e-6
 spec      l    rsm    eh     gmax    last term   cutoff
  Ni       0*   1.30  -0.40   5.718    1.22E-06     259 
  Ni       1    1.30  -0.40   6.080    1.60E-06     283 
  Ni       2    1.00  -0.40   8.508    1.20E-06     749 
  Ni       3*   1.30  -0.40   6.806    1.29E-06     387 
  Ni       4    1.30  -0.40   7.165    1.58E-06     459 
  Ni       0*   1.30  -1.20   5.718    1.22E-06     259 
  Ni       1    1.30  -1.20   6.080    1.60E-06     283 
  Ni       2    1.00  -1.20   8.508    1.20E-06     749 
 

 iors  : read restart file (binary, mesh density)
         use from  restart file: ef window, positions, pnu 
         ignore in restart file: *

 ... sudmft job=1: make projectors
 5 nonzero (5 inequivalent) matrix elements

   channels for 2nd spin block
   m1   m2   cix chn(1) chn(2)
    1    1    1     1     6
    2    2    1     2     7
    3    3    1     3     8
    4    4    1     4     9
    5    5    1     5    10

  correlated channels:
  chan equiv m1   m2   isp icix  cix  ib
    1    1    1    1    1    1    1    1
    2    2    2    2    1    1    1    1
    3    3    3    3    1    1    1    1
    4    4    4    4    1    1    1    1
    5    5    5    5    1    1    1    1
    6    6    1    1    2    1    1    1
    7    7    2    2    2    1    1    1
    8    8    3    3    2    1    1    1
    9    9    4    4    2    1    1    1
   10   10    5    5    2    1    1    1
 
 indmfl file expects 1 DMFT block(s) (max dimension 5),  10 matrix elements
 Missing sigma file : create it and exit ...
 Exit 0 done writing template sigma
 CPU time:  0.372s   Wall clock 0.435s  at  08:55:31 28.03.2019  on  nmscde001136.nms.kcl
rdcmd:  lmfdmft ni -vnk=12 --rs=1,0 -vnsp=2 -vtpd1=4 -vbigbas=3 -vehmax=-.4 -vrwa=0.3359 -vemaxs=2.5 -vconvc=1e-6 -vbeta=.3 --ldadc=82.2 -job=1 -vindmfl=f
 -------------  START LMFDMFT - ref: 4e6aff7 --------------
#rf  171: (file atparms) element Ni : Z1=28 eref=0 mom=.6 rsm=1.3 rsmd=1
 HEADER Master file for multinary TM compounds, suitable for GW

 LMFDMFT:  nbas = 1  nspec = 1  verb 31,20
 special:: forces
 pot:      spin-pol, XC:BH, read Sigma
 float:    float P v6-style, v6-ebar
 autoread: none
 bz:       metal(5), tetra 
 dmft:     1 blocks  (1 independent)

                Plat                                  Qlat
   0.000000   0.500000   0.500000       -1.000000   1.000000   1.000000
   0.500000   0.000000   0.500000        1.000000  -1.000000   1.000000
   0.500000   0.500000   0.000000        1.000000   1.000000  -1.000000
   alat = 6.64  Cell vol = 73.188736

 LATTC:  as= 2.000   tol=1.00e-8   alat= 6.64000   awald= 0.478   lmax=6
         r1=  1.965   nkd= 135      q1=  5.929   nkg= 229

 SGROUP: 1 symmetry operations from 0 generators
 SYMLAT: Bravais system is cubic with 48 symmetry operations.
 SYMCRY: crystal invariant under 48 symmetry operations for tol=1e-5
 GROUPG: the following are sufficient to generate the space group:
         i*r3(1,1,-1) r4x
         i*r3(1,1,-1) r4x
 MKSYM:  found 48 space group operations
         no attempt to add inversion symmetry
 BZMESH: 72 irreducible QP from 1728 ( 12 12 12 )  shift= F F F
 TETIRR: sorting 10368 tetrahedra ... 247 inequivalent ones found

 species data:  augmentation                           density
 spec       rmt   rsma lmxa kmxa      lmxl     rg   rsmv  kmxv foca   rfoca
 Ni       2.230  0.892    4    4         4  0.558  1.115    15    1   0.892

 MSHSIZ: mesh has 11 x 11 x 11 divisions; length 0.427, 0.427, 0.427
         generated from gmax = 9.0 a.u. : 893 vectors of 1331 (67%)

 GVLIST: gmax = 9 a.u. created 893 vectors of 1331 (67%)
         mesh has 11 x 11 x 11 divisions; length 0.427, 0.427, 0.427
 SGVSYM: 40 symmetry stars found for 893 reciprocal lattice vectors
 

 Makidx:  hamiltonian dimensions Low, Int, High, Negl: 39 0 16 20
 kappa   Low   Int   High  L+I  L+I+H  Neglected
   1      25     0     0    25    25       0
   2       9     0    16     9    25       0
   3       5     0     0     5     5      20
  all     39     0    16    39    55      20
 suham :  25 augmentation channels, 25 local potential channels  Maximum lmxa=4

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol=1.0e-6
 spec      l    rsm    eh     gmax    last term   cutoff
  Ni       0*   1.30  -0.40   5.718    1.22E-06     259 
  Ni       1    1.30  -0.40   6.080    1.60E-06     283 
  Ni       2    1.00  -0.40   8.508    1.20E-06     749 
  Ni       3*   1.30  -0.40   6.806    1.29E-06     387 
  Ni       4    1.30  -0.40   7.165    1.58E-06     459 
  Ni       0*   1.30  -1.20   5.718    1.22E-06     259 
  Ni       1    1.30  -1.20   6.080    1.60E-06     283 
  Ni       2    1.00  -1.20   8.508    1.20E-06     749 
 

 iors  : read restart file (binary, mesh density)
         use from  restart file: ef window, positions, pnu 
         ignore in restart file: *

 ... sudmft job=1: make projectors
 5 nonzero (5 inequivalent) matrix elements

   channels for 2nd spin block
   m1   m2   cix chn(1) chn(2)
    1    1    1     1     6
    2    2    1     2     7
    3    3    1     3     8
    4    4    1     4     9
    5    5    1     5    10

  correlated channels:
  chan equiv m1   m2   isp icix  cix  ib
    1    1    1    1    1    1    1    1
    2    2    2    2    1    1    1    1
    3    3    3    3    1    1    1    1
    4    4    4    4    1    1    1    1
    5    5    5    5    1    1    1    1
    6    6    1    1    2    1    1    1
    7    7    2    2    2    1    1    1
    8    8    3    3    2    1    1    1
    9    9    4    4    2    1    1    1
   10   10    5    5    2    1    1    1
 
 indmfl file expects 1 DMFT block(s) (max dimension 5),  10 matrix elements
 Read DMFT sigma from file.  Scale omega, sigma by 0.0735
 DMFT sigma is zero ... no double counting

 SUDMFT: projector #2 with k-resolved norm.  8 bands in window (1,8)
         1999 Matsubara frequencies, interval (0.0628318,251.139) eV
 iodmftdc: reading d.c. from command-line ... parsed 1 element(s)
 iodmftdc: dc =  82.2 82.2 82.2 82.2 82.2 82.2 82.2 82.2 82.2 82.2

 --- BNDFP:  begin iteration 1 of 1 ---

 Average es pot at rmt = 0.660281  avg sphere pot = 0.672025   vconst = -0.660281

 site  1  z= 28.0  rmt= 2.23038  nr=389   a=0.025  nlml=25  rg=0.558  Vfloat=T
 sm core charge = 0.454474 (sphere) + 0.008781 (spillout) = 0.463255
 potential shift to crystal energy zero:    0.000043


 RDSIGM: read file sigm and create COMPLEX sigma(R) by FT ...
         Sigm will be approximated by:  diagonal Sigma for high and low states 
         Approximate sigma for energies E(lda)>2.5
         For high states read constant sigii from sigm file 
         sigm file has 72 irreducible QP: nk = ( 12 12 12 )  shift= F F F
         average SE read from sigm file:   0.137886 0.144245
 hft2rs: make neighbor table for r.s. hamiltonian using range = 8 * alat
 pairc:  8583 pairs total  8583 is max cluster size
 hft2rs: found 1728 connecting vectors out of 1728 possible for FFT
 symiax: enlarged neighbor table from 1728 to 1957 pairs (48 symops)

 q-points in full BZ where sigma calculable ...
 BZMESH: 1728 irreducible QP from 1728 ( 12 12 12 )  shift= F F F
 Irr. qp for which sigma is calculated ...
 BZMESH: 72 irreducible QP from 1728 ( 12 12 12 )  shift= F F F
 
 hft2rs created hrs:  ndhrs=39  max Re(hrs) = 0.806  max Im(hrs) = 1.65e-9
 hft2rs created hrs:  ndhrs=39  max Re(hrs) = 0.795  max Im(hrs) = 1.71e-9
 check FT s(R) against s(q) ... maximum error = 8.3e-15 < tol (5e-6)
 check FT s(R) against s(q) ... maximum error = 1e-14 < tol (5e-6) spin 2
 rsmsym: symmetrizing complex s(1..1957) using 48 group operations 
 symstr: max asymmetry = 1.11e-16
 check FT s(R) against s(q) ... maximum error = 1.1e-14
 check FT s(R) against s(q) ... maximum error = 1.5e-14 spin 2

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 BNDFP:  Write evals,evecs to file for 72 qp
 bndfp:  kpt 1 of 72, k=  0.00000  0.00000  0.00000
 -0.6368 -0.0832 -0.0832 -0.0832 -0.0178 -0.0178  1.8177  2.1247  2.1247
 bndfp:  kpt 1 of 72, k=  0.00000  0.00000  0.00000
 -0.6401 -0.0300 -0.0300 -0.0300  0.0362  0.0362  1.8127  2.1286  2.1286
 bndfp:  kpt 11 of 72, k=  -0.25000  0.25000  0.41667
 -0.3917 -0.1221 -0.0779 -0.0680 -0.0049  0.0039  0.9124  1.7116  1.7119
 bndfp:  kpt 11 of 72, k=  -0.25000  0.25000  0.41667
 -0.3885 -0.0778 -0.0257 -0.0143  0.0500  0.0598  0.9376  1.7164  1.7291
 bndfp:  kpt 21 of 72, k=  -0.16667  0.16667  0.50000
 -0.3837 -0.1421 -0.0806 -0.0506 -0.0034  0.0069  1.0527  1.5345  1.5745
 bndfp:  kpt 21 of 72, k=  -0.16667  0.16667  0.50000
 -0.3818 -0.0941 -0.0295  0.0053  0.0504  0.0624  1.0744  1.5495  1.5871
 bndfp:  kpt 31 of 72, k=  -0.25000  0.25000  0.75000
 -0.2268 -0.1707 -0.0934 -0.0344  0.0203  0.1854  0.6513  1.1143  1.3881
 bndfp:  kpt 31 of 72, k=  -0.25000  0.25000  0.75000
 -0.1983 -0.1297 -0.0478  0.0208  0.0785  0.2235  0.6797  1.1328  1.3988
 bndfp:  kpt 41 of 72, k=  -0.08333  0.08333  0.91667
 -0.2598 -0.2024 -0.0003  0.0121  0.0235  0.2016  0.7661  1.0217  1.1686
 bndfp:  kpt 41 of 72, k=  -0.08333  0.08333  0.91667
 -0.2259 -0.1564  0.0509  0.0664  0.0837  0.2166  0.7990  1.0321  1.1740
 bndfp:  kpt 51 of 72, k=  -0.16667  0.33333  0.66667
 -0.2459 -0.1558 -0.1043 -0.0415  0.0120  0.1224  0.7122  1.2512  1.3302
 bndfp:  kpt 51 of 72, k=  -0.16667  0.33333  0.66667
 -0.2263 -0.1141 -0.0560  0.0132  0.0693  0.1680  0.7412  1.2672  1.3442
 bndfp:  kpt 61 of 72, k=  0.00000  0.16667  0.83333
 -0.2563 -0.1932 -0.0356  0.0096  0.0165  0.1569  0.8911  1.0404  1.0940
 bndfp:  kpt 61 of 72, k=  0.00000  0.16667  0.83333
 -0.2258 -0.1469  0.0102  0.0627  0.0764  0.1864  0.9171  1.0558  1.1021
 bndfp:  kpt 71 of 72, k=  0.00000  0.33333  1.00000
 -0.2208 -0.1776 -0.0823 -0.0089  0.0259  0.4135  0.7417  0.7948  0.9275
 bndfp:  kpt 71 of 72, k=  0.00000  0.33333  1.00000
 -0.1837 -0.1328 -0.0366  0.0441  0.0870  0.4331  0.7746  0.8174  0.9393
 BNDFP:  evecs saved in file evec

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.065072;  10.000000 electrons;  D(Ef):   26.004
         Sum occ. bands:   -0.9329041  incl. Bloechl correction:   -0.001612
         Mag. moment:       0.768953

 Saved qp weights ...

 Harris energy:
 sumev=       -0.932904  val*vef=    -151.345482   sumtv=     150.412578
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    2934.175594
 rhoeps=    -121.970989     utot=   -5997.641462    ehar=   -3035.024279

 Make and renormalize projectors ...

 Find chemical potential mu ...
 Seek mu for 10 electrons ... 10.01261 electrons at Ef0=0.065072.  D(Ef0)=33.316
 getVnew: v1=0e0 N1=1.261e-2  v2=3.785e-4 N2=5.202e-5  bracket=F  est=3.8022e-4
 getVnew: v1=3.785e-4 N1=5.202e-5  v2=3.802e-4 N2=-4.951e-6  bracket=T  est=3.8007e-4
 getVnew: v1=3.785e-4 N1=5.202e-5  v2=3.801e-4 N2=-6.391e-12  bracket=T  est=3.8007e-4
 mu = 0.064692 = Ef0-0.00038.  Deviation from neutrality = 0e0  D(mu)=33.033
 Electron charge:  10.000000   moment:  0.785477 spin resolved DOS:    2.880  30.152

 Make gloc, delta, eimp ...
 ... Writing gkloc and gloc... ... gkloc finished... Writing files delta and eimp1 ...

 Check SC condition skipped (missing information about previous iteration)
 Exit 0 done making DMFT hybridization function
 CPU time:  72.423s   Wall clock 72.485s  at  08:56:45 28.03.2019  on  nmscde001136.nms.kcl
rdcmd:  lmfdmft ni -vnk=12 --rs=1,0 -vnsp=2 -vtpd1=4 -vbigbas=3 -vehmax=-.4 -vrwa=0.3359 -vemaxs=2.5 -vconvc=1e-6 -vbeta=.3 --ldadc=82.2 -job=1 --gprt:mode=3:nom=20
 -------------  START LMFDMFT - ref: 4e6aff7 --------------
#rf  171: (file atparms) element Ni : Z1=28 eref=0 mom=.6 rsm=1.3 rsmd=1
 HEADER Master file for multinary TM compounds, suitable for GW

 LMFDMFT:  nbas = 1  nspec = 1  verb 31,20
 special:: forces
 pot:      spin-pol, XC:BH, read Sigma
 float:    float P v6-style, v6-ebar
 autoread: none
 bz:       metal(5), tetra 
 dmft:     1 blocks  (1 independent)

                Plat                                  Qlat
   0.000000   0.500000   0.500000       -1.000000   1.000000   1.000000
   0.500000   0.000000   0.500000        1.000000  -1.000000   1.000000
   0.500000   0.500000   0.000000        1.000000   1.000000  -1.000000
   alat = 6.64  Cell vol = 73.188736

 LATTC:  as= 2.000   tol=1.00e-8   alat= 6.64000   awald= 0.478   lmax=6
         r1=  1.965   nkd= 135      q1=  5.929   nkg= 229

 SGROUP: 1 symmetry operations from 0 generators
 SYMLAT: Bravais system is cubic with 48 symmetry operations.
 SYMCRY: crystal invariant under 48 symmetry operations for tol=1e-5
 GROUPG: the following are sufficient to generate the space group:
         i*r3(1,1,-1) r4x
         i*r3(1,1,-1) r4x
 MKSYM:  found 48 space group operations
         no attempt to add inversion symmetry
 BZMESH: 72 irreducible QP from 1728 ( 12 12 12 )  shift= F F F
 TETIRR: sorting 10368 tetrahedra ... 247 inequivalent ones found

 species data:  augmentation                           density
 spec       rmt   rsma lmxa kmxa      lmxl     rg   rsmv  kmxv foca   rfoca
 Ni       2.230  0.892    4    4         4  0.558  1.115    15    1   0.892

 MSHSIZ: mesh has 11 x 11 x 11 divisions; length 0.427, 0.427, 0.427
         generated from gmax = 9.0 a.u. : 893 vectors of 1331 (67%)

 GVLIST: gmax = 9 a.u. created 893 vectors of 1331 (67%)
         mesh has 11 x 11 x 11 divisions; length 0.427, 0.427, 0.427
 SGVSYM: 40 symmetry stars found for 893 reciprocal lattice vectors
 

 Makidx:  hamiltonian dimensions Low, Int, High, Negl: 39 0 16 20
 kappa   Low   Int   High  L+I  L+I+H  Neglected
   1      25     0     0    25    25       0
   2       9     0    16     9    25       0
   3       5     0     0     5     5      20
  all     39     0    16    39    55      20
 suham :  25 augmentation channels, 25 local potential channels  Maximum lmxa=4

 sugcut:  make orbital-dependent reciprocal vector cutoffs for tol=1.0e-6
 spec      l    rsm    eh     gmax    last term   cutoff
  Ni       0*   1.30  -0.40   5.718    1.22E-06     259 
  Ni       1    1.30  -0.40   6.080    1.60E-06     283 
  Ni       2    1.00  -0.40   8.508    1.20E-06     749 
  Ni       3*   1.30  -0.40   6.806    1.29E-06     387 
  Ni       4    1.30  -0.40   7.165    1.58E-06     459 
  Ni       0*   1.30  -1.20   5.718    1.22E-06     259 
  Ni       1    1.30  -1.20   6.080    1.60E-06     283 
  Ni       2    1.00  -1.20   8.508    1.20E-06     749 
 

 iors  : read restart file (binary, mesh density)
         use from  restart file: ef window, positions, pnu 
         ignore in restart file: *

 ... sudmft job=1: make projectors
 5 nonzero (5 inequivalent) matrix elements

   channels for 2nd spin block
   m1   m2   cix chn(1) chn(2)
    1    1    1     1     6
    2    2    1     2     7
    3    3    1     3     8
    4    4    1     4     9
    5    5    1     5    10

  correlated channels:
  chan equiv m1   m2   isp icix  cix  ib
    1    1    1    1    1    1    1    1
    2    2    2    2    1    1    1    1
    3    3    3    3    1    1    1    1
    4    4    4    4    1    1    1    1
    5    5    5    5    1    1    1    1
    6    6    1    1    2    1    1    1
    7    7    2    2    2    1    1    1
    8    8    3    3    2    1    1    1
    9    9    4    4    2    1    1    1
   10   10    5    5    2    1    1    1
 
 indmfl file expects 1 DMFT block(s) (max dimension 5),  10 matrix elements
 Read DMFT sigma from file.  Scale omega, sigma by 0.0735
 DMFT sigma is zero ... no double counting

 SUDMFT: projector #2 with k-resolved norm.  8 bands in window (1,8)
         1999 Matsubara frequencies, interval (0.0628318,251.139) eV
 iodmftdc: reading d.c. from command-line ... parsed 1 element(s)
 iodmftdc: dc =  82.2 82.2 82.2 82.2 82.2 82.2 82.2 82.2 82.2 82.2

 --- BNDFP:  begin iteration 1 of 1 ---

 Average es pot at rmt = 0.660281  avg sphere pot = 0.672025   vconst = -0.660281

 site  1  z= 28.0  rmt= 2.23038  nr=389   a=0.025  nlml=25  rg=0.558  Vfloat=T
 sm core charge = 0.454474 (sphere) + 0.008781 (spillout) = 0.463255
 potential shift to crystal energy zero:    0.000043


 RDSIGM: read file sigm and create COMPLEX sigma(R) by FT ...
         Sigm will be approximated by:  diagonal Sigma for high and low states 
         Approximate sigma for energies E(lda)>2.5
         For high states read constant sigii from sigm file 
         sigm file has 72 irreducible QP: nk = ( 12 12 12 )  shift= F F F
         average SE read from sigm file:   0.137886 0.144245
 hft2rs: make neighbor table for r.s. hamiltonian using range = 8 * alat
 pairc:  8583 pairs total  8583 is max cluster size
 hft2rs: found 1728 connecting vectors out of 1728 possible for FFT
 symiax: enlarged neighbor table from 1728 to 1957 pairs (48 symops)

 q-points in full BZ where sigma calculable ...
 BZMESH: 1728 irreducible QP from 1728 ( 12 12 12 )  shift= F F F
 Irr. qp for which sigma is calculated ...
 BZMESH: 72 irreducible QP from 1728 ( 12 12 12 )  shift= F F F
 
 hft2rs created hrs:  ndhrs=39  max Re(hrs) = 0.806  max Im(hrs) = 1.65e-9
 hft2rs created hrs:  ndhrs=39  max Re(hrs) = 0.795  max Im(hrs) = 1.71e-9
 check FT s(R) against s(q) ... maximum error = 8.3e-15 < tol (5e-6)
 check FT s(R) against s(q) ... maximum error = 1e-14 < tol (5e-6) spin 2
 rsmsym: symmetrizing complex s(1..1957) using 48 group operations 
 symstr: max asymmetry = 1.11e-16
 check FT s(R) against s(q) ... maximum error = 1.1e-14
 check FT s(R) against s(q) ... maximum error = 1.5e-14 spin 2

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 BNDFP:  Write evals,evecs to file for 72 qp
 bndfp:  kpt 1 of 72, k=  0.00000  0.00000  0.00000
 -0.6368 -0.0832 -0.0832 -0.0832 -0.0178 -0.0178  1.8177  2.1247  2.1247
 bndfp:  kpt 1 of 72, k=  0.00000  0.00000  0.00000
 -0.6401 -0.0300 -0.0300 -0.0300  0.0362  0.0362  1.8127  2.1286  2.1286
 bndfp:  kpt 11 of 72, k=  -0.25000  0.25000  0.41667
 -0.3917 -0.1221 -0.0779 -0.0680 -0.0049  0.0039  0.9124  1.7116  1.7119
 bndfp:  kpt 11 of 72, k=  -0.25000  0.25000  0.41667
 -0.3885 -0.0778 -0.0257 -0.0143  0.0500  0.0598  0.9376  1.7164  1.7291
 bndfp:  kpt 21 of 72, k=  -0.16667  0.16667  0.50000
 -0.3837 -0.1421 -0.0806 -0.0506 -0.0034  0.0069  1.0527  1.5345  1.5745
 bndfp:  kpt 21 of 72, k=  -0.16667  0.16667  0.50000
 -0.3818 -0.0941 -0.0295  0.0053  0.0504  0.0624  1.0744  1.5495  1.5871
 bndfp:  kpt 31 of 72, k=  -0.25000  0.25000  0.75000
 -0.2268 -0.1707 -0.0934 -0.0344  0.0203  0.1854  0.6513  1.1143  1.3881
 bndfp:  kpt 31 of 72, k=  -0.25000  0.25000  0.75000
 -0.1983 -0.1297 -0.0478  0.0208  0.0785  0.2235  0.6797  1.1328  1.3988
 bndfp:  kpt 41 of 72, k=  -0.08333  0.08333  0.91667
 -0.2598 -0.2024 -0.0003  0.0121  0.0235  0.2016  0.7661  1.0217  1.1686
 bndfp:  kpt 41 of 72, k=  -0.08333  0.08333  0.91667
 -0.2259 -0.1564  0.0509  0.0664  0.0837  0.2166  0.7990  1.0321  1.1740
 bndfp:  kpt 51 of 72, k=  -0.16667  0.33333  0.66667
 -0.2459 -0.1558 -0.1043 -0.0415  0.0120  0.1224  0.7122  1.2512  1.3302
 bndfp:  kpt 51 of 72, k=  -0.16667  0.33333  0.66667
 -0.2263 -0.1141 -0.0560  0.0132  0.0693  0.1680  0.7412  1.2672  1.3442
 bndfp:  kpt 61 of 72, k=  0.00000  0.16667  0.83333
 -0.2563 -0.1932 -0.0356  0.0096  0.0165  0.1569  0.8911  1.0404  1.0940
 bndfp:  kpt 61 of 72, k=  0.00000  0.16667  0.83333
 -0.2258 -0.1469  0.0102  0.0627  0.0764  0.1864  0.9171  1.0558  1.1021
 bndfp:  kpt 71 of 72, k=  0.00000  0.33333  1.00000
 -0.2208 -0.1776 -0.0823 -0.0089  0.0259  0.4135  0.7417  0.7948  0.9275
 bndfp:  kpt 71 of 72, k=  0.00000  0.33333  1.00000
 -0.1837 -0.1328 -0.0366  0.0441  0.0870  0.4331  0.7746  0.8174  0.9393
 BNDFP:  evecs saved in file evec

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.065072;  10.000000 electrons;  D(Ef):   26.004
         Sum occ. bands:   -0.9329041  incl. Bloechl correction:   -0.001612
         Mag. moment:       0.768953

 Saved qp weights ...

 Harris energy:
 sumev=       -0.932904  val*vef=    -151.345482   sumtv=     150.412578
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    2934.175594
 rhoeps=    -121.970989     utot=   -5997.641462    ehar=   -3035.024279

 Make and renormalize projectors ...

 Find chemical potential mu ...
 Seek mu for 10 electrons ... 10.01261 electrons at Ef0=0.065072.  D(Ef0)=33.316
 getVnew: v1=0e0 N1=1.261e-2  v2=3.785e-4 N2=5.202e-5  bracket=F  est=3.8022e-4
 getVnew: v1=3.785e-4 N1=5.202e-5  v2=3.802e-4 N2=-4.951e-6  bracket=T  est=3.8007e-4
 getVnew: v1=3.785e-4 N1=5.202e-5  v2=3.801e-4 N2=-6.391e-12  bracket=T  est=3.8007e-4
 mu = 0.064692 = Ef0-0.00038.  Deviation from neutrality = 0e0  D(mu)=33.033
 Electron charge:  10.000000   moment:  0.785477 spin resolved DOS:    2.880  30.152
 k-integrated gloc
 ... Writing gkloc and gloc... ... gkloc finished... Exit 0 : --gprt:mode=3:nom=20 generated file gkloc
 CPU time:  46.180s   Wall clock 46.220s  at  08:57:31 28.03.2019  on  nmscde001136.nms.kcl
