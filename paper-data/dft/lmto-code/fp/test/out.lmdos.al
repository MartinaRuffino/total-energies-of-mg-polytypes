 ---------------------  START LMDOS  ----------------------

 LMDOS:    nbas = 1  nspec = 1  verb 31,30
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
 BZMESH: 8797 irreducible QP from 373248 ( 72 72 72 )  shift= F F F

 ASADOS: reading weights from file MOMS
         expecting file to be resolved by l
         file has 0 channel(s) => generate total DOS
         Using npts=101  emin=-0.75  emax=-0.65
         adjust emin,emax: point 13 coincides with ef=-0.736545
 IOMOMQ: read 8797 qp  efermi=0.001952  vmtz=0.000000

 ASADOS:  make dos for 101 points from 9 bands in window (-0.75,-0.65)
 options  mode: |v|; kres@ef=-0.736545 totdos
          bands:  1 2 3
          remake qp with 1 symop
 BZMESH: 373248 irreducible QP from 373248 ( 72 72 72 )  shift= F F F
 TETIRR: sorting 2239488 tetrahedra ... 2239488 inequivalent ones found

 CONTET: 20088 tetra enclose ef=-0.736545.  Resolved by band:  20088 0 0

 <|v|> (spin 1) at Ef = 0.4637 a.u. = 0.614 x 10^6 m/s.  Resolve by band:
  ib   <|v|>       ib   <|v|>       ib   <|v|>
   1   0.4637  |    2   0.0000  |    3   0.0000  |
 total DOS(Ef) = 0.776607
 Exit 0 LMDOS 
 CPU time:  1.720s   Wall clock 1.757s  at  19:49:35 15.08.2018  on  localhost.localdomai
 ---------------------  START LMDOS  ----------------------

 LMDOS:    nbas = 1  nspec = 1  verb 31,30
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
 BZMESH: 8797 irreducible QP from 373248 ( 72 72 72 )  shift= F F F

 ASADOS: reading weights from file MOMS
         expecting file to be resolved by l
         file has 0 channel(s) => generate total DOS
         Using npts=101  emin=-0.05  emax=0.05
         adjust emin,emax: point 52 coincides with ef=0.001952
 IOMOMQ: read 8797 qp  efermi=0.001952  vmtz=0.000000

 ASADOS:  make dos for 101 points from 9 bands in window (-0.05,0.05)
 options  mode: sigma(ballistic),v=(100); kres@ef=0.001952 totdos
          bands:  1 2 3
          remake qp with 1 symop
 BZMESH: 373248 irreducible QP from 373248 ( 72 72 72 )  shift= F F F
 TETIRR: sorting 2239488 tetrahedra ... 2239488 inequivalent ones found

 CONTET: 162024 tetra enclose ef=0.001952.  Resolved by band:  0 134568 27456

 <v/2.c> (spin 1) at Ef = 0.3137 a.u. = 0.4154 x 10^6 m/s.  Resolve by band:
  ib   <v/2.c>     ib   <v/2.c>     ib   <v/2.c>
   1   0.0000  |    2   0.3340  |    3   0.2143  |
 total DOS(Ef) = 1.610590
 Exit 0 LMDOS 
 CPU time:  4.636s   Wall clock 4.666s  at  19:49:40 15.08.2018  on  localhost.localdomai
