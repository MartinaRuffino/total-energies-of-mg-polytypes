 ----------------------  START BLM  -----------------------

 ... Reading lattice data ... read plat ... no group generators

 LATPAR:     A=  5.41632     B=  5.41632      C=  5.41632
         ALPHA=  90.0000  BETA=  90.0000  GAMMA=  90.0000
 ... Found and reading from HAM category ...
 ... Found and reading from BZ category ...

 ... Reading site data ...  2 sites, 2 species found

 ... Reading species data ... species data read for 2 species

 ... Complete basis for supplied symmetry
 SGROUP: 1 symmetry operations from 1 generators
 SYMLAT: Bravais system is cubic with 48 symmetry operations.
 SYMCRY: crystal invariant under 48 symmetry operations for tol=1e-5
 FIXPOS: shifted site positions by average 4.88e-18
 SYMCRY: crystal invariant under 48 symmetry operations for tol=1e-5
 GROUPG: 2 generator(s) were added to complete the group:
         i*i i*r3d r2(0,1,1)
         i*i i*r3d r2(0,1,1)

 Lattice vectors : alat = 5.41632 a.u.  vol=158.896 a.u.

                 Plat                     Conventional unit cell          As multiples of Plat
   1.0000000  0.0000000  0.0000000    1.0000000  0.0000000  0.0000000    1.0000  0.0000  0.0000
   0.0000000  1.0000000  0.0000000    0.0000000  1.0000000  0.0000000    0.0000  1.0000  0.0000
   0.0000000  0.0000000  1.0000000    0.0000000  0.0000000  1.0000000    0.0000  0.0000  1.0000

 Basis vectors after sorting and shortening:
 site spec        pos (Cartesian coordinates)         pos (multiples of plat)
   1  Cr        0.000000   0.000000   0.000000    0.000000   0.000000   0.000000
   2  Cr2       0.500000   0.500000   0.500000    0.500000   0.500000   0.500000

 ... Make sphere radii

 pairc, ib=  1:  59 neighbors in range  1.950*alat =  10.56
 pairc:  59 pairs total  59 is max cluster size

 pairc, ib=  2:  59 neighbors in range  1.950*alat =  10.56
 pairc:  59 pairs total  59 is max cluster size

 makrm0: initial MT radii from first estat potential maximum
  site   spec            rmt       rmt-       rmt-       rold   lock
                                  <spec avg>  spec-min
    1    1:Cr           2.3450     0.0000     0.0000     0.0000
    2    2:Cr2          2.3450     0.0000     0.0000     0.0000
 SCLWSR:  mode = 100  vol = 158.896 a.u.   Initial sphere packing = 68%  scaled to 62.1%
 constr omax1=  -3.0  -3.0  -3.0 %    omax2= 100.0 100.0 100.0 %
 actual omax1=   0.0   0.0   0.0 %    omax2=   0.0   0.0   0.0 %

 spec  name        old rmax    new rmax     ratio
   1   Cr          2.344997    2.274975    0.970140
   2   Cr2         2.344997    2.274975    0.970140

 ... Create input file actrl.cr (express mode 3)
 IOSITE: wrote to file 'site', 2 sites
