rdcmd:  lmstr gd      -vlmf=f -vso=f -vnk=12 -vnk2=8 -vcoref=t -vldau=f -vfm=t -vlmx=2 --no-iactiv
 ---------------------  START LMSTR  ----------------------
 HEADER Gd

 rdctrl: reset global max nl from 4 to 3

 LMSTR:    nbas = 2  nspec = 1  verb 30,20

                Plat                                  Qlat
   0.866025  -0.500000   0.000000        1.154701   0.000000   0.000000
   0.000000   1.000000   0.000000        0.577350   1.000000   0.000000
   0.000000   0.000000   1.588000        0.000000   0.000000   0.629723
   alat = 6.8786  Cell vol = 447.590947

 LATTC:  as= 2.000   tol=1.00e-16   alat= 6.87860   awald= 0.261   lmax=6
         r1=  4.228   nkd= 233      q1=  4.156   nkg= 391

 SGROUP: 1 symmetry operations from 0 generators
 SYMLAT: Bravais system is hexagonal with 24 symmetry operations.
 SYMCRY: crystal invariant under 24 symmetry operations for tol=1e-5
 GROUPG: the following are sufficient to generate the space group:
         i*r3z:(-1*sqrt(3)/6,1/2,-0.794) r2z:(-1*sqrt(3)/6,1/2,-0.794) r2(1/2,sqrt(3)/2,0)
         i*r3z::(-1/3,1/3,-1/2) r2z::(-1/3,1/3,-1/2) r2(1/2,sqrt(3)/2,0)
 MKSYM:  found 24 space group operations ... includes inversion

 ASASTR: strux for 2nd gen LMTOs, and Sdot
         avw = 3.766352

 pairc, ib=  1:  39 neighbors in range  1.752*alat =  12.05
 pairc, ib=  2:  39 neighbors in range  1.752*alat =  12.05
 pairc:  78 pairs total  39 is max cluster size
 strscr:  generated 78 inequivalent strux from 78 total
 Exit 0 LMSTR 
 CPU time:  0.058s   Wall clock 0.101s  at  14:44:42 05.04.2018  on  localhost.localdomai
rdcmd:  lm gd -vnit=0 -vlmf=f -vso=f -vnk=12 -vnk2=8 -vcoref=t -vldau=f -vfm=t -vlmx=2 --no-iactiv
 -----------------------  START LM  -----------------------
 HEADER Gd

 rdctrl: reset global max nl from 4 to 3

 LM:       nbas = 2  nspec = 1  verb 30,20
 pot:      spin-pol, XC:BH, read Sigma
 bz:       metal, tetra 

                Plat                                  Qlat
   0.866025  -0.500000   0.000000        1.154701   0.000000   0.000000
   0.000000   1.000000   0.000000        0.577350   1.000000   0.000000
   0.000000   0.000000   1.588000        0.000000   0.000000   0.629723
   alat = 6.8786  Cell vol = 447.590947

 LATTC:  as= 2.000   tol=1.00e-16   alat= 6.87860   awald= 0.261   lmax=6
         r1=  4.228   nkd= 233      q1=  4.156   nkg= 391

 SGROUP: 1 symmetry operations from 0 generators
 SYMLAT: Bravais system is hexagonal with 24 symmetry operations.
 SYMCRY: crystal invariant under 24 symmetry operations for tol=1e-5
 GROUPG: the following are sufficient to generate the space group:
         i*r3z:(-1*sqrt(3)/6,1/2,-0.794) r2z:(-1*sqrt(3)/6,1/2,-0.794) r2(1/2,sqrt(3)/2,0)
         i*r3z::(-1/3,1/3,-1/2) r2z::(-1/3,1/3,-1/2) r2(1/2,sqrt(3)/2,0)
 MKSYM:  found 24 space group operations ... includes inversion
 BZMESH: 95 irreducible QP from 1152 ( 12 12 8 )  shift= F F F
 TETIRR: sorting 6912 tetrahedra ... 1172 inequivalent ones found

 GETZV:  6 valence electrons
ATOM=Gd   Z=64  Qc=61  R=3.766352  Qv=0  mom=7.74362  a=0.02  nr=701

 Class        Qtot       Qbak       Vmad     Vh(Rmax)    V(Rmax)
 Gd         0.000000   0.000000   0.000000   0.000000  -0.563360
 Sum Q=0.000000  Emad=0.000000(0.000000)  Vmtz=-0.600000
 LM: it 0 of 0  ehk0=-1.116297  pv=-0.0002  mmom=15.4872  seref=-45054.01449
 Exit 0 LM 
 CPU time:  0.265s   Wall clock 0.314s  at  14:44:42 05.04.2018  on  localhost.localdomai
rdcmd:  lm gd -vnit=1 -vlmf=f -vso=f -vnk=12 -vnk2=8 -vcoref=t -vldau=f -vfm=t -vlmx=2 --no-iactiv
 -----------------------  START LM  -----------------------
 HEADER Gd

 rdctrl: reset global max nl from 4 to 3

 LM:       nbas = 2  nspec = 1  verb 30,20
 pot:      spin-pol, XC:BH, read Sigma
 bz:       metal, tetra 

                Plat                                  Qlat
   0.866025  -0.500000   0.000000        1.154701   0.000000   0.000000
   0.000000   1.000000   0.000000        0.577350   1.000000   0.000000
   0.000000   0.000000   1.588000        0.000000   0.000000   0.629723
   alat = 6.8786  Cell vol = 447.590947

 LATTC:  as= 2.000   tol=1.00e-16   alat= 6.87860   awald= 0.261   lmax=6
         r1=  4.228   nkd= 233      q1=  4.156   nkg= 391

 SGROUP: 1 symmetry operations from 0 generators
 SYMLAT: Bravais system is hexagonal with 24 symmetry operations.
 SYMCRY: crystal invariant under 24 symmetry operations for tol=1e-5
 GROUPG: the following are sufficient to generate the space group:
         i*r3z:(-1*sqrt(3)/6,1/2,-0.794) r2z:(-1*sqrt(3)/6,1/2,-0.794) r2(1/2,sqrt(3)/2,0)
         i*r3z::(-1/3,1/3,-1/2) r2z::(-1/3,1/3,-1/2) r2(1/2,sqrt(3)/2,0)
 MKSYM:  found 24 space group operations ... includes inversion
 BZMESH: 95 irreducible QP from 1152 ( 12 12 8 )  shift= F F F
 TETIRR: sorting 6912 tetrahedra ... 1172 inequivalent ones found

 GETZV:  6 valence electrons

 Class        Qtot       Qbak       Vmad     Vh(Rmax)    V(Rmax)
 Gd         0.000000   0.000000   0.000000   0.000000  -0.563360
 Sum Q=0.000000  Emad=0.000000(0.000000)  Vmtz=-0.600000

 Makidx:  hamiltonian dimensions Low, Int, High, Negl: 18 0 0 0

 lmasa (warning): no sigm file found ... LDA calculation only

 --- BNDASA : iteration 1 ---
 SECMAT:  kpt 1 of 95, k=  0.00000  0.00000  0.00000
 -0.5082 -0.3037 -0.0333 -0.0284 -0.0284 -0.0044  0.0236  0.0236  0.0681
 SECMAT:  kpt 1 of 95, k=  0.00000  0.00000  0.00000
 -0.4430 -0.2333  0.0472  0.0552  0.0552  0.0843  0.1136  0.1136  0.1630

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.110275;   6.000000 electrons;  D(Ef):   18.906
         Sum occ. bands:   -1.3128494  incl. Bloechl correction:   -0.001037
         Mag. moment:       1.487146

 Saved qp weights ...

 ... Generating total DOS

 LM: ehf=-1.1156284  ehk=-1.1156405  sumev=-1.3128494  delsev=0.0006566
 mixing: mode=A  nmix=2  beta=.3  kill=4
 PQMIX:  read 0 iter from file mixm.  RMS DQ=6.20e-6
 AMIX: nmix=0 mmix=8  nelts=24  beta=0.3  tm=10  rmsdel=6.2e-6

 GETZV:  6 valence electrons
ATOM=Gd   Z=64  Qc=61  R=3.766352  Qv=0  mom=7.74361  a=0.02  nr=701

 Class        Qtot       Qbak       Vmad     Vh(Rmax)    V(Rmax)
 Gd         0.000000   0.000000   0.000000   0.000000  -0.563360
 Sum Q=0.000000  Emad=0.000000(0.000000)  Vmtz=-0.600000
 Magnetization from sphere program: 15.487216  from output rho: 15.487146
 SV:   1  6.196e-6 0.3000  3.637e-6      -1.1156405215.487146 L

   it  1  of  1    ehf=      -1.115628   ehk=      -1.115641
                rms dq=  0.000006    tol= 0.000010  more=F
c nit=1 lmf=0 so=0 nk=12 nk2=8 coref=1 ldau=0 fm=1 lmx=2 mmom=15.4871457 ehf=-1.1156284 ehk=-1.1156405

 Jolly good show !  You converged to rms DQ=  0.000004
 Exit 0 LM 
 CPU time:  0.207s   Wall clock 0.258s  at  14:44:42 05.04.2018  on  localhost.localdomai
rdcmd:  rm -f mixm.gd
rdcmd:  lm gd -vnit=0 -vlmf=f -vso=f -vnk=12 -vnk2=8 -vcoref=t -vldau=f -vfm=f -vlmx=2 --no-iactiv
 -----------------------  START LM  -----------------------
 HEADER Gd

 rdctrl: reset global max nl from 4 to 3

 LM:       nbas = 2  nspec = 2  verb 30,20
 pot:      spin-pol, XC:BH, read Sigma
 bz:       metal, tetra 

                Plat                                  Qlat
   0.866025  -0.500000   0.000000        1.154701   0.000000   0.000000
   0.000000   1.000000   0.000000        0.577350   1.000000   0.000000
   0.000000   0.000000   1.588000        0.000000   0.000000   0.629723
   alat = 6.8786  Cell vol = 447.590947

 LATTC:  as= 2.000   tol=1.00e-16   alat= 6.87860   awald= 0.261   lmax=6
         r1=  4.228   nkd= 233      q1=  4.156   nkg= 391

 SGROUP: 1 symmetry operations from 0 generators
 SYMLAT: Bravais system is hexagonal with 24 symmetry operations.
 SYMCRY: crystal invariant under 12 symmetry operations for tol=1e-5
 GROUPG: the following are sufficient to generate the space group:
         i*r6z r2(1/2,sqrt(3)/2,0)
         i*r6z r2(1/2,sqrt(3)/2,0)
 MKSYM:  found 12 space group operations; adding inversion generated 24 ops
 BZMESH: 95 irreducible QP from 1152 ( 12 12 8 )  shift= F F F
 TETIRR: sorting 6912 tetrahedra ... 1172 inequivalent ones found

 GETZV:  6 valence electrons
ATOM=Gd   Z=64  Qc=61  R=3.766352  Qv=0  mom=7.49062  a=0.02  nr=701
ATOM=Gd1  Z=64  Qc=61  R=3.766352  Qv=0  mom=-7.49062  a=0.02  nr=701

 Class        Qtot       Qbak       Vmad     Vh(Rmax)    V(Rmax)
 Gd         0.000000   0.000000   0.000000   0.000000  -0.563971
 Gd1        0.000000   0.000000   0.000000   0.000000  -0.563971
 Sum Q=0.000000  Emad=0.000000(0.000000)  Vmtz=-0.600000
 LM: it 0 of 0  ehk0=-1.110997  pv=-0.0109  mmom=0  seref=-45054.01449
 Exit 0 LM 
 CPU time:  0.372s   Wall clock 0.409s  at  14:44:43 05.04.2018  on  localhost.localdomai
rdcmd:  lm gd -vnit=1 -vlmf=f -vso=f -vnk=12 -vnk2=8 -vcoref=t -vldau=f -vfm=f -vlmx=2 --no-iactiv
 -----------------------  START LM  -----------------------
 HEADER Gd

 rdctrl: reset global max nl from 4 to 3

 LM:       nbas = 2  nspec = 2  verb 30,20
 pot:      spin-pol, XC:BH, read Sigma
 bz:       metal, tetra 

                Plat                                  Qlat
   0.866025  -0.500000   0.000000        1.154701   0.000000   0.000000
   0.000000   1.000000   0.000000        0.577350   1.000000   0.000000
   0.000000   0.000000   1.588000        0.000000   0.000000   0.629723
   alat = 6.8786  Cell vol = 447.590947

 LATTC:  as= 2.000   tol=1.00e-16   alat= 6.87860   awald= 0.261   lmax=6
         r1=  4.228   nkd= 233      q1=  4.156   nkg= 391

 SGROUP: 1 symmetry operations from 0 generators
 SYMLAT: Bravais system is hexagonal with 24 symmetry operations.
 SYMCRY: crystal invariant under 12 symmetry operations for tol=1e-5
 GROUPG: the following are sufficient to generate the space group:
         i*r6z r2(1/2,sqrt(3)/2,0)
         i*r6z r2(1/2,sqrt(3)/2,0)
 MKSYM:  found 12 space group operations; adding inversion generated 24 ops
 BZMESH: 95 irreducible QP from 1152 ( 12 12 8 )  shift= F F F
 TETIRR: sorting 6912 tetrahedra ... 1172 inequivalent ones found

 GETZV:  6 valence electrons

 Class        Qtot       Qbak       Vmad     Vh(Rmax)    V(Rmax)
 Gd         0.000000   0.000000   0.000000   0.000000  -0.563971
 Gd1        0.000000   0.000000   0.000000   0.000000  -0.563971
 Sum Q=0.000000  Emad=0.000000(0.000000)  Vmtz=-0.600000

 Makidx:  hamiltonian dimensions Low, Int, High, Negl: 18 0 0 0

 lmasa (warning): no sigm file found ... LDA calculation only

 --- BNDASA : iteration 1 ---
 SECMAT:  kpt 1 of 95, k=  0.00000  0.00000  0.00000
 -0.4774 -0.2694 -0.0056 -0.0056  0.0064  0.0400  0.0879  0.0879  0.1040
 SECMAT:  kpt 1 of 95, k=  0.00000  0.00000  0.00000
 -0.4774 -0.2694 -0.0056 -0.0056  0.0064  0.0400  0.0879  0.0879  0.1040

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:     -0.110132;   6.000000 electrons;  D(Ef):   38.171
         Sum occ. bands:   -1.2971112  incl. Bloechl correction:   -0.001691
         Mag. moment:       0.000000

 Saved qp weights ...

 ... Generating total DOS

 LM: ehf=-1.1075291  ehk=-1.1075549  sumev=-1.2971112  delsev=0.0034424

 Mixpqc: ic  ngrp  group
     2 elts in group  1  Q= 3.000000
 mixing: mode=A  nmix=2  beta=.3  kill=4
 PQMIX:  read 0 iter from file mixm.  RMS DQ=2.48e-7
 AMIX: nmix=0 mmix=8  nelts=48  beta=0.3  tm=10  rmsdel=2.52e-7

 GETZV:  6 valence electrons
ATOM=Gd   Z=64  Qc=61  R=3.766352  Qv=0  mom=7.49062  a=0.02  nr=701
ATOM=Gd1  Z=64  Qc=61  R=3.766352  Qv=0  mom=-7.49062  a=0.02  nr=701

 Class        Qtot       Qbak       Vmad     Vh(Rmax)    V(Rmax)
 Gd         0.000000   0.000000   0.000000   0.000000  -0.563971
 Gd1        0.000000   0.000000   0.000000   0.000000  -0.563971
 Sum Q=0.000000  Emad=0.000000(0.000000)  Vmtz=-0.600000
 Magnetization from sphere program: 0.000000  from output rho: 0.000000
 SV:   1  2.482e-7 0.3000  6.038e-8      -1.10755491 0.000000 L

   it  1  of  1    ehf=      -1.107529   ehk=      -1.107555
                rms dq=  0.000000    tol= 0.000010  more=F
c nit=1 lmf=0 so=0 nk=12 nk2=8 coref=1 ldau=0 fm=0 lmx=2 mmom=0 ehf=-1.1075291 ehk=-1.1075549

 Jolly good show !  You converged to rms DQ=  0.000000
 Exit 0 LM 
 CPU time:  0.272s   Wall clock 0.316s  at  14:44:43 05.04.2018  on  localhost.localdomai
