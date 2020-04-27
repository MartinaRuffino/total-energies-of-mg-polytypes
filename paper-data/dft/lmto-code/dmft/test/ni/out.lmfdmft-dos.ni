rdcmd:  lmf ni -vnk=24 -vnsp=2 -vtpd1=4 -vbigbas=3 -vehmax=-.4 -vrwa=0.3359 -vemaxs=2.5 -vconvc=1e-6 -vbeta=.3 -job=1 --dos:npts=201:ef0:rdm:window=-.2,.2 --quit=dos
 ----------------------  START LMF  -----------------------
#rf  171: (file atparms) element Ni : Z1=28 eref=0 mom=.6 rsm=1.3 rsmd=1
 HEADER Master file for multinary TM compounds, suitable for GW

 LMF:      nbas = 1  nspec = 1  verb 31,20
 special:  forces
 pot:      spin-pol, XC:BH, read Sigma
 float:    float P v6-style, v6-ebar
 autoread: none
 bz:       metal(5), tetra 

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
 BZMESH: 413 irreducible QP from 13824 ( 24 24 24 )  shift= F F F
 TETIRR: sorting 82944 tetrahedra ... 1864 inequivalent ones found

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

 --- BNDFP:  begin iteration 1 of 20 ---

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
 check FT s(R) against s(q) ... maximum error = 8.1e-15 < tol (5e-6)
 check FT s(R) against s(q) ... maximum error = 8.7e-15 < tol (5e-6) spin 2
 rsmsym: symmetrizing complex s(1..1957) using 48 group operations 
 symstr: max asymmetry = 1.11e-16
 check FT s(R) against s(q) ... maximum error = 1e-14
 check FT s(R) against s(q) ... maximum error = 1.5e-14 spin 2

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 Start first of two band passes ...
 bndfp:  kpt 1 of 413, k=  0.00000  0.00000  0.00000
 -0.6368 -0.0832 -0.0832 -0.0832 -0.0178 -0.0178  1.8177  2.1247  2.1247
 bndfp:  kpt 1 of 413, k=  0.00000  0.00000  0.00000
 -0.6401 -0.0300 -0.0300 -0.0300  0.0362  0.0362  1.8127  2.1286  2.1286
 bndfp:  kpt 11 of 413, k=  -0.41667  0.41667  0.41667
 -0.3001 -0.0823 -0.0823 -0.0622  0.0130  0.0130  0.5587  1.7475  1.9017
 bndfp:  kpt 11 of 413, k=  -0.41667  0.41667  0.41667
 -0.2794 -0.0463 -0.0305 -0.0305  0.0704  0.0704  0.6000  1.7430  1.9145
 bndfp:  kpt 21 of 413, k=  -0.29167  0.29167  0.37500
 -0.3863 -0.1139 -0.0762 -0.0731 -0.0042  0.0021  0.8550  1.7594  1.8112
 bndfp:  kpt 21 of 413, k=  -0.29167  0.29167  0.37500
 -0.3817 -0.0725 -0.0237 -0.0201  0.0516  0.0582  0.8831  1.7597  1.8280
 bndfp:  kpt 31 of 413, k=  0.29167  -0.29167  -0.20833
 -0.4530 -0.1148 -0.0713 -0.0698 -0.0180 -0.0055  1.0566  1.7322  1.9219
 bndfp:  kpt 31 of 413, k=  0.29167  -0.29167  -0.20833
 -0.4531 -0.0675 -0.0178 -0.0162  0.0375  0.0497  1.0794  1.7391  1.9279
 bndfp:  kpt 41 of 413, k=  -0.16667  0.16667  0.33333
 -0.4896 -0.1161 -0.0696 -0.0652 -0.0258 -0.0081  1.2246  1.7527  1.7611
 bndfp:  kpt 41 of 413, k=  -0.16667  0.16667  0.33333
 -0.4915 -0.0660 -0.0159 -0.0104  0.0289  0.0466  1.2445  1.7667  1.7702
 bndfp:  kpt 51 of 413, k=  0.41667  -0.41667  -0.25000
 -0.3296 -0.1100 -0.0830 -0.0737 -0.0037  0.0280  0.7275  1.5478  1.8653
 bndfp:  kpt 51 of 413, k=  0.41667  -0.41667  -0.25000
 -0.3197 -0.0741 -0.0311 -0.0211  0.0532  0.0825  0.7601  1.5605  1.8674
 bndfp:  kpt 61 of 413, k=  -0.12500  0.12500  0.37500
 -0.4852 -0.1227 -0.0677 -0.0596 -0.0305 -0.0052  1.2957  1.6693  1.7151
 bndfp:  kpt 61 of 413, k=  -0.12500  0.12500  0.37500
 -0.4870 -0.0719 -0.0140 -0.0039  0.0242  0.0494  1.3152  1.6866  1.7231
 bndfp:  kpt 71 of 413, k=  0.45833  -0.45833  -0.20833
 -0.2999 -0.1159 -0.0912 -0.0690 -0.0051  0.0553  0.6999  1.3931  1.8475
 bndfp:  kpt 71 of 413, k=  0.45833  -0.45833  -0.20833
 -0.2873 -0.0797 -0.0397 -0.0159  0.0519  0.1087  0.7330  1.4074  1.8500
 bndfp:  kpt 81 of 413, k=  -0.16667  0.16667  0.50000
 -0.3837 -0.1421 -0.0806 -0.0506 -0.0034  0.0069  1.0527  1.5345  1.5745
 bndfp:  kpt 81 of 413, k=  -0.16667  0.16667  0.50000
 -0.3818 -0.0941 -0.0295  0.0053  0.0504  0.0624  1.0744  1.5495  1.5871
 bndfp:  kpt 91 of 413, k=  0.41667  -0.41667  -0.08333
 -0.3492 -0.1301 -0.0887 -0.0551 -0.0282  0.0404  0.9539  1.2420  1.9192
 bndfp:  kpt 91 of 413, k=  0.41667  -0.41667  -0.08333
 -0.3455 -0.0836 -0.0380 -0.0005  0.0290  0.0951  0.9782  1.2587  1.9230
 bndfp:  kpt 101 of 413, k=  -0.29167  0.29167  0.70833
 -0.2310 -0.1563 -0.0963 -0.0493  0.0196  0.1681  0.5810  1.2308  1.4687
 bndfp:  kpt 101 of 413, k=  -0.29167  0.29167  0.70833
 -0.2024 -0.1189 -0.0492  0.0053  0.0776  0.2067  0.6158  1.2465  1.4802
 bndfp:  kpt 111 of 413, k=  -0.08333  0.08333  0.58333
 -0.3557 -0.1619 -0.0693 -0.0277 -0.0052  0.0078  1.2081  1.3866  1.3936
 bndfp:  kpt 111 of 413, k=  -0.08333  0.08333  0.58333
 -0.3513 -0.1132 -0.0219  0.0298  0.0495  0.0628  1.2251  1.4023  1.4055
 bndfp:  kpt 121 of 413, k=  0.50000  -0.50000  0.00000
 -0.2660 -0.1340 -0.1145 -0.0431 -0.0252  0.1017  0.8797  0.9930  1.7540
 bndfp:  kpt 121 of 413, k=  0.50000  -0.50000  0.00000
 -0.2555 -0.0877 -0.0660  0.0104  0.0342  0.1550  0.9084  1.0096  1.7711
 bndfp:  kpt 131 of 413, k=  0.62500  -0.62500  -0.04167
 -0.1897 -0.1627 -0.1171 -0.0322 -0.0018  0.2818  0.6236  0.8704  1.5070
 bndfp:  kpt 131 of 413, k=  0.62500  -0.62500  -0.04167
 -0.1598 -0.1191 -0.0721  0.0215  0.0577  0.3159  0.6567  0.8921  1.5206
 bndfp:  kpt 141 of 413, k=  0.66667  -0.66667  0.00000
 -0.1819 -0.1792 -0.1078 -0.0242  0.0030  0.3661  0.5703  0.8128  1.4330
 bndfp:  kpt 141 of 413, k=  0.66667  -0.66667  0.00000
 -0.1484 -0.1347 -0.0633  0.0291  0.0624  0.3959  0.6024  0.8372  1.4445
 bndfp:  kpt 151 of 413, k=  -0.08333  0.08333  0.91667
 -0.2598 -0.2024 -0.0003  0.0121  0.0235  0.2016  0.7661  1.0217  1.1686
 bndfp:  kpt 151 of 413, k=  -0.08333  0.08333  0.91667
 -0.2259 -0.1564  0.0509  0.0664  0.0837  0.2166  0.7990  1.0321  1.1740
 bndfp:  kpt 161 of 413, k=  -0.12500  0.20833  0.29167
 -0.5084 -0.1114 -0.0698 -0.0671 -0.0281 -0.0095  1.2745  1.6808  1.9188
 bndfp:  kpt 161 of 413, k=  -0.12500  0.20833  0.29167
 -0.5106 -0.0608 -0.0159 -0.0124  0.0261  0.0451  1.2936  1.6932  1.9264
 bndfp:  kpt 171 of 413, k=  -0.08333  0.16667  0.33333
 -0.5068 -0.1161 -0.0682 -0.0629 -0.0327 -0.0074  1.3434  1.6151  1.8587
 bndfp:  kpt 171 of 413, k=  -0.08333  0.16667  0.33333
 -0.5090 -0.0650 -0.0142 -0.0072  0.0214  0.0470  1.3625  1.6290  1.8665
 bndfp:  kpt 181 of 413, k=  0.50000  -0.41667  -0.25000
 -0.2950 -0.1103 -0.0890 -0.0715  0.0026  0.0515  0.6517  1.4739  1.7535
 bndfp:  kpt 181 of 413, k=  0.50000  -0.41667  -0.25000
 -0.2797 -0.0769 -0.0378 -0.0188  0.0594  0.1039  0.6875  1.4871  1.7622
 bndfp:  kpt 191 of 413, k=  -0.12500  0.20833  0.45833
 -0.4092 -0.1351 -0.0776 -0.0539 -0.0162  0.0073  1.0953  1.4959  1.6974
 bndfp:  kpt 191 of 413, k=  -0.12500  0.20833  0.45833
 -0.4086 -0.0865 -0.0256  0.0018  0.0390  0.0624  1.1170  1.5116  1.7087
 bndfp:  kpt 201 of 413, k=  0.45833  -0.37500  -0.12500
 -0.3436 -0.1295 -0.0886 -0.0588 -0.0193  0.0385  0.8949  1.3172  1.8065
 bndfp:  kpt 201 of 413, k=  0.45833  -0.37500  -0.12500
 -0.3388 -0.0846 -0.0375 -0.0046  0.0375  0.0931  0.9208  1.3330  1.8164
 bndfp:  kpt 211 of 413, k=  -0.25000  0.33333  0.66667
 -0.2434 -0.1481 -0.0973 -0.0529  0.0180  0.1316  0.6088  1.3168  1.4417
 bndfp:  kpt 211 of 413, k=  -0.25000  0.33333  0.66667
 -0.2190 -0.1105 -0.0490  0.0017  0.0758  0.1741  0.6437  1.3320  1.4541
 bndfp:  kpt 221 of 413, k=  -0.04167  0.12500  0.54167
 -0.3817 -0.1540 -0.0724 -0.0322 -0.0170  0.0077  1.2484  1.3714  1.4948
 bndfp:  kpt 221 of 413, k=  -0.04167  0.12500  0.54167
 -0.3797 -0.1048 -0.0229  0.0253  0.0379  0.0625  1.2660  1.3863  1.5095
 bndfp:  kpt 231 of 413, k=  0.54167  -0.45833  -0.04167
 -0.2645 -0.1355 -0.1146 -0.0434 -0.0209  0.1021  0.8501  1.0298  1.6535
 bndfp:  kpt 231 of 413, k=  0.54167  -0.45833  -0.04167
 -0.2534 -0.0900 -0.0656  0.0102  0.0377  0.1552  0.8780  1.0478  1.6678
 bndfp:  kpt 241 of 413, k=  0.66667  -0.58333  -0.08333
 -0.1912 -0.1629 -0.1154 -0.0354  0.0019  0.2828  0.5843  0.9272  1.4464
 bndfp:  kpt 241 of 413, k=  0.66667  -0.58333  -0.08333
 -0.1598 -0.1212 -0.0700  0.0186  0.0606  0.3162  0.6189  0.9477  1.4602
 bndfp:  kpt 251 of 413, k=  0.70833  -0.62500  -0.04167
 -0.1811 -0.1789 -0.1094 -0.0255  0.0051  0.3676  0.5586  0.8352  1.3772
 bndfp:  kpt 251 of 413, k=  0.70833  -0.62500  -0.04167
 -0.1476 -0.1346 -0.0647  0.0280  0.0644  0.3971  0.5914  0.8588  1.3896
 bndfp:  kpt 261 of 413, k=  -0.04167  0.12500  0.87500
 -0.2593 -0.1991 -0.0128  0.0113  0.0215  0.1742  0.8210  1.0358  1.1446
 bndfp:  kpt 261 of 413, k=  -0.04167  0.12500  0.87500
 -0.2268 -0.1529  0.0360  0.0654  0.0814  0.1955  0.8504  1.0473  1.1520
 bndfp:  kpt 271 of 413, k=  -0.08333  0.25000  0.41667
 -0.4281 -0.1289 -0.0763 -0.0550 -0.0288  0.0089  1.1446  1.4220  1.8270
 bndfp:  kpt 271 of 413, k=  -0.08333  0.25000  0.41667
 -0.4282 -0.0798 -0.0241  0.0009  0.0261  0.0638  1.1656  1.4376  1.8369
 bndfp:  kpt 281 of 413, k=  -0.12500  0.29167  0.54167
 -0.3268 -0.1434 -0.0931 -0.0483 -0.0056  0.0393  0.9127  1.3276  1.5702
 bndfp:  kpt 281 of 413, k=  -0.12500  0.29167  0.54167
 -0.3209 -0.0975 -0.0424  0.0067  0.0500  0.0928  0.9380  1.3429  1.5846
 bndfp:  kpt 291 of 413, k=  0.45833  -0.29167  -0.04167
 -0.3875 -0.1347 -0.0838 -0.0462 -0.0321  0.0231  1.1138  1.2647  1.7601
 bndfp:  kpt 291 of 413, k=  0.45833  -0.29167  -0.04167
 -0.3862 -0.0861 -0.0324  0.0096  0.0230  0.0779  1.1354  1.2792  1.7720
 bndfp:  kpt 301 of 413, k=  0.66667  -0.50000  -0.16667
 -0.2137 -0.1488 -0.1087 -0.0481  0.0069  0.2141  0.5556  1.1276  1.4129
 bndfp:  kpt 301 of 413, k=  0.66667  -0.50000  -0.16667
 -0.1834 -0.1110 -0.0617  0.0061  0.0649  0.2515  0.5926  1.1450  1.4269
 bndfp:  kpt 311 of 413, k=  -0.20833  0.37500  0.79167
 -0.2022 -0.1699 -0.1049 -0.0341  0.0157  0.2788  0.5486  1.0729  1.2595
 bndfp:  kpt 311 of 413, k=  -0.20833  0.37500  0.79167
 -0.1670 -0.1315 -0.0590  0.0203  0.0742  0.3114  0.5812  1.0916  1.2733
 bndfp:  kpt 321 of 413, k=  -0.16667  0.33333  0.83333
 -0.2090 -0.1755 -0.0974 -0.0210  0.0172  0.2977  0.6213  0.9651  1.1960
 bndfp:  kpt 321 of 413, k=  -0.16667  0.33333  0.83333
 -0.1741 -0.1335 -0.0527  0.0331  0.0761  0.3290  0.6493  0.9863  1.2090
 bndfp:  kpt 331 of 413, k=  0.79167  -0.62500  -0.04167
 -0.1946 -0.1754 -0.1037 -0.0195  0.0127  0.4605  0.5147  0.8107  1.2155
 bndfp:  kpt 331 of 413, k=  0.79167  -0.62500  -0.04167
 -0.1566 -0.1333 -0.0589  0.0337  0.0723  0.4880  0.5438  0.8359  1.2287
 bndfp:  kpt 341 of 413, k=  0.00000  0.16667  1.00000
 -0.2526 -0.2005 -0.0160  0.0103  0.0259  0.2696  0.7199  1.0051  1.0574
 bndfp:  kpt 341 of 413, k=  0.00000  0.16667  1.00000
 -0.2172 -0.1547  0.0345  0.0640  0.0870  0.2834  0.7560  1.0187  1.0641
 bndfp:  kpt 351 of 413, k=  -0.08333  0.33333  0.66667
 -0.2491 -0.1578 -0.1102 -0.0281  0.0021  0.1171  0.8295  1.1186  1.2990
 bndfp:  kpt 351 of 413, k=  -0.08333  0.33333  0.66667
 -0.2324 -0.1128 -0.0625  0.0259  0.0594  0.1649  0.8540  1.1358  1.3138
 bndfp:  kpt 361 of 413, k=  -0.04167  0.29167  0.70833
 -0.2479 -0.1675 -0.1012 -0.0133  0.0032  0.1242  0.9087  1.0628  1.1891
 bndfp:  kpt 361 of 413, k=  -0.04167  0.29167  0.70833
 -0.2282 -0.1213 -0.0557  0.0397  0.0611  0.1697  0.9309  1.0799  1.2049
 bndfp:  kpt 371 of 413, k=  -0.08333  0.33333  0.83333
 -0.2169 -0.1732 -0.0973 -0.0124  0.0170  0.2796  0.7371  0.9451  1.0532
 bndfp:  kpt 371 of 413, k=  -0.08333  0.33333  0.83333
 -0.1834 -0.1288 -0.0532  0.0409  0.0764  0.3113  0.7619  0.9679  1.0654
 bndfp:  kpt 381 of 413, k=  0.00000  0.25000  0.91667
 -0.2378 -0.1894 -0.0552  0.0019  0.0235  0.2949  0.7843  0.9095  0.9983
 bndfp:  kpt 381 of 413, k=  0.00000  0.25000  0.91667
 -0.2026 -0.1439 -0.0087  0.0551  0.0843  0.3178  0.8150  0.9279  1.0080
 bndfp:  kpt 391 of 413, k=  -0.04167  0.37500  0.79167
 -0.2107 -0.1630 -0.1163 -0.0159  0.0119  0.2576  0.7662  0.9488  1.0585
 bndfp:  kpt 391 of 413, k=  -0.04167  0.37500  0.79167
 -0.1799 -0.1183 -0.0724  0.0369  0.0713  0.2928  0.7914  0.9677  1.0750
 bndfp:  kpt 401 of 413, k=  0.83333  -0.50000  0.00000
 -0.1915 -0.1482 -0.1365 -0.0229  0.0163  0.4079  0.6286  0.8304  1.0455
 bndfp:  kpt 401 of 413, k=  0.83333  -0.50000  0.00000
 -0.1550 -0.1051 -0.0923  0.0301  0.0761  0.4361  0.6567  0.8532  1.0619
 bndfp:  kpt 411 of 413, k=  0.91667  -0.50000  0.00000
 -0.1985 -0.1408 -0.1375 -0.0221  0.0233  0.5188  0.6002  0.8036  0.9106
 bndfp:  kpt 411 of 413, k=  0.91667  -0.50000  0.00000
 -0.1605 -0.0973 -0.0936  0.0307  0.0840  0.5458  0.6270  0.8302  0.9274
 
 hambls: smallest nmax encountered for sigm middle block = 10

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.065080;  10.000000 electrons;  D(Ef):   29.644
         Sum occ. bands:   -0.9327044  incl. Bloechl correction:   -0.000408
         Mag. moment:       0.757262

 Saved qp weights ...

 ... Generating total DOS
 Exit 0 completed band pass 
 CPU time:  57.760s   Wall clock 57.811s  at  17:49:39 22.10.2018  on  localhost.localdomai
rdcmd:  cp dos.ni doslmf.ni
rdcmd:  lmfdmft ni -vnk=24 -vnsp=2 -vtpd1=4 -vbigbas=3 -vehmax=-.4 -vrwa=0.3359 -vemaxs=2.5 -vconvc=1e-6 -vbeta=.3 --ldadc=82.2 -job=1
 --------------------  START LMFDMFT  ---------------------
#rf  171: (file atparms) element Ni : Z1=28 eref=0 mom=.6 rsm=1.3 rsmd=1
 HEADER Master file for multinary TM compounds, suitable for GW

 LMFDMFT:  nbas = 1  nspec = 1  verb 31,20
 special:: forces
 pot:      spin-pol, XC:BH, read Sigma
 float:    float P v6-style, v6-ebar
 autoread: none
 bz:       metal(5), tetra 

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
 BZMESH: 413 irreducible QP from 13824 ( 24 24 24 )  shift= F F F
 TETIRR: sorting 82944 tetrahedra ... 1864 inequivalent ones found

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

 Reading indmfl file ... 1 total (1 inequivalent) correlated blocks, among 1 sites
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
 CPU time:  0.370s   Wall clock 0.412s  at  17:49:39 22.10.2018  on  localhost.localdomai
rdcmd:  lmfdmft ni -vnk=24 -vnsp=2 -vtpd1=4 -vbigbas=3 -vehmax=-.4 -vrwa=0.3359 -vemaxs=2.5 -vconvc=1e-6 -vbeta=.3 --ldadc=82.2 -job=1 --dos:npts=201:ef0:rdm:window=-.2,.2
 --------------------  START LMFDMFT  ---------------------
#rf  171: (file atparms) element Ni : Z1=28 eref=0 mom=.6 rsm=1.3 rsmd=1
 HEADER Master file for multinary TM compounds, suitable for GW

 LMFDMFT:  nbas = 1  nspec = 1  verb 31,20
 special:: forces
 pot:      spin-pol, XC:BH, read Sigma
 float:    float P v6-style, v6-ebar
 autoread: none
 bz:       metal(5), tetra 

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
 BZMESH: 413 irreducible QP from 13824 ( 24 24 24 )  shift= F F F
 TETIRR: sorting 82944 tetrahedra ... 1864 inequivalent ones found

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

 Reading indmfl file ... 1 total (1 inequivalent) correlated blocks, among 1 sites
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
 check FT s(R) against s(q) ... maximum error = 8.1e-15 < tol (5e-6)
 check FT s(R) against s(q) ... maximum error = 8.7e-15 < tol (5e-6) spin 2
 rsmsym: symmetrizing complex s(1..1957) using 48 group operations 
 symstr: max asymmetry = 1.11e-16
 check FT s(R) against s(q) ... maximum error = 1e-14
 check FT s(R) against s(q) ... maximum error = 1.5e-14 spin 2

 subzi: tetrahedron integration of bands; tetrahedron integration of density


 BNDFP:  Write evals,evecs to file for 413 qp
 bndfp:  kpt 1 of 413, k=  0.00000  0.00000  0.00000
 -0.6368 -0.0832 -0.0832 -0.0832 -0.0178 -0.0178  1.8177  2.1247  2.1247
 bndfp:  kpt 1 of 413, k=  0.00000  0.00000  0.00000
 -0.6401 -0.0300 -0.0300 -0.0300  0.0362  0.0362  1.8127  2.1286  2.1286
 bndfp:  kpt 11 of 413, k=  -0.41667  0.41667  0.41667
 -0.3001 -0.0823 -0.0823 -0.0622  0.0130  0.0130  0.5587  1.7475  1.9017
 bndfp:  kpt 11 of 413, k=  -0.41667  0.41667  0.41667
 -0.2794 -0.0463 -0.0305 -0.0305  0.0704  0.0704  0.6000  1.7430  1.9145
 bndfp:  kpt 21 of 413, k=  -0.29167  0.29167  0.37500
 -0.3863 -0.1139 -0.0762 -0.0731 -0.0042  0.0021  0.8550  1.7594  1.8112
 bndfp:  kpt 21 of 413, k=  -0.29167  0.29167  0.37500
 -0.3817 -0.0725 -0.0237 -0.0201  0.0516  0.0582  0.8831  1.7597  1.8280
 bndfp:  kpt 31 of 413, k=  -0.70833  0.70833  0.79167
 -0.4530 -0.1148 -0.0713 -0.0698 -0.0180 -0.0055  1.0566  1.7322  1.9219
 bndfp:  kpt 31 of 413, k=  -0.70833  0.70833  0.79167
 -0.4531 -0.0675 -0.0178 -0.0162  0.0375  0.0497  1.0794  1.7391  1.9279
 bndfp:  kpt 41 of 413, k=  -0.16667  0.16667  0.33333
 -0.4896 -0.1161 -0.0696 -0.0652 -0.0258 -0.0081  1.2246  1.7527  1.7611
 bndfp:  kpt 41 of 413, k=  -0.16667  0.16667  0.33333
 -0.4915 -0.0660 -0.0159 -0.0104  0.0289  0.0466  1.2445  1.7667  1.7702
 bndfp:  kpt 51 of 413, k=  -0.58333  0.58333  0.75000
 -0.3296 -0.1100 -0.0830 -0.0737 -0.0037  0.0280  0.7275  1.5478  1.8653
 bndfp:  kpt 51 of 413, k=  -0.58333  0.58333  0.75000
 -0.3197 -0.0741 -0.0311 -0.0211  0.0532  0.0825  0.7601  1.5605  1.8674
 bndfp:  kpt 61 of 413, k=  -0.12500  0.12500  0.37500
 -0.4852 -0.1227 -0.0677 -0.0596 -0.0305 -0.0052  1.2957  1.6693  1.7151
 bndfp:  kpt 61 of 413, k=  -0.12500  0.12500  0.37500
 -0.4870 -0.0719 -0.0140 -0.0039  0.0242  0.0494  1.3152  1.6866  1.7231
 bndfp:  kpt 71 of 413, k=  -0.54167  0.54167  0.79167
 -0.2999 -0.1159 -0.0912 -0.0690 -0.0051  0.0553  0.6999  1.3931  1.8475
 bndfp:  kpt 71 of 413, k=  -0.54167  0.54167  0.79167
 -0.2873 -0.0797 -0.0397 -0.0159  0.0519  0.1087  0.7330  1.4074  1.8500
 bndfp:  kpt 81 of 413, k=  -0.16667  0.16667  0.50000
 -0.3837 -0.1421 -0.0806 -0.0506 -0.0034  0.0069  1.0527  1.5345  1.5745
 bndfp:  kpt 81 of 413, k=  -0.16667  0.16667  0.50000
 -0.3818 -0.0941 -0.0295  0.0053  0.0504  0.0624  1.0744  1.5495  1.5871
 bndfp:  kpt 91 of 413, k=  -0.58333  0.58333  0.91667
 -0.3492 -0.1301 -0.0887 -0.0551 -0.0282  0.0404  0.9539  1.2420  1.9192
 bndfp:  kpt 91 of 413, k=  -0.58333  0.58333  0.91667
 -0.3455 -0.0836 -0.0380 -0.0005  0.0290  0.0951  0.9782  1.2587  1.9230
 bndfp:  kpt 101 of 413, k=  -0.29167  0.29167  0.70833
 -0.2310 -0.1563 -0.0963 -0.0493  0.0196  0.1681  0.5810  1.2308  1.4687
 bndfp:  kpt 101 of 413, k=  -0.29167  0.29167  0.70833
 -0.2024 -0.1189 -0.0492  0.0053  0.0776  0.2067  0.6158  1.2465  1.4802
 bndfp:  kpt 111 of 413, k=  -0.08333  0.08333  0.58333
 -0.3557 -0.1619 -0.0693 -0.0277 -0.0052  0.0078  1.2081  1.3866  1.3936
 bndfp:  kpt 111 of 413, k=  -0.08333  0.08333  0.58333
 -0.3513 -0.1132 -0.0219  0.0298  0.0495  0.0628  1.2251  1.4023  1.4055
 bndfp:  kpt 121 of 413, k=  -0.50000  0.50000  1.00000
 -0.2660 -0.1340 -0.1145 -0.0431 -0.0252  0.1017  0.8797  0.9930  1.7540
 bndfp:  kpt 121 of 413, k=  -0.50000  0.50000  1.00000
 -0.2555 -0.0877 -0.0660  0.0104  0.0342  0.1550  0.9084  1.0096  1.7711
 bndfp:  kpt 131 of 413, k=  -0.37500  0.37500  0.95833
 -0.1897 -0.1627 -0.1171 -0.0322 -0.0018  0.2818  0.6236  0.8704  1.5070
 bndfp:  kpt 131 of 413, k=  -0.37500  0.37500  0.95833
 -0.1598 -0.1191 -0.0721  0.0215  0.0577  0.3159  0.6567  0.8921  1.5206
 bndfp:  kpt 141 of 413, k=  -0.33333  0.33333  1.00000
 -0.1819 -0.1792 -0.1078 -0.0242  0.0030  0.3661  0.5703  0.8128  1.4330
 bndfp:  kpt 141 of 413, k=  -0.33333  0.33333  1.00000
 -0.1484 -0.1347 -0.0633  0.0291  0.0624  0.3959  0.6024  0.8372  1.4445
 bndfp:  kpt 151 of 413, k=  -0.08333  0.08333  0.91667
 -0.2598 -0.2024 -0.0003  0.0121  0.0235  0.2016  0.7661  1.0217  1.1686
 bndfp:  kpt 151 of 413, k=  -0.08333  0.08333  0.91667
 -0.2259 -0.1564  0.0509  0.0664  0.0837  0.2166  0.7990  1.0321  1.1740
 bndfp:  kpt 161 of 413, k=  -0.12500  0.20833  0.29167
 -0.5084 -0.1114 -0.0698 -0.0671 -0.0281 -0.0095  1.2745  1.6808  1.9188
 bndfp:  kpt 161 of 413, k=  -0.12500  0.20833  0.29167
 -0.5106 -0.0608 -0.0159 -0.0124  0.0261  0.0451  1.2936  1.6932  1.9264
 bndfp:  kpt 171 of 413, k=  -0.08333  0.16667  0.33333
 -0.5068 -0.1161 -0.0682 -0.0629 -0.0327 -0.0074  1.3434  1.6151  1.8587
 bndfp:  kpt 171 of 413, k=  -0.08333  0.16667  0.33333
 -0.5090 -0.0650 -0.0142 -0.0072  0.0214  0.0470  1.3625  1.6290  1.8665
 bndfp:  kpt 181 of 413, k=  -0.50000  0.58333  0.75000
 -0.2950 -0.1103 -0.0890 -0.0715  0.0026  0.0515  0.6517  1.4739  1.7535
 bndfp:  kpt 181 of 413, k=  -0.50000  0.58333  0.75000
 -0.2797 -0.0769 -0.0378 -0.0188  0.0594  0.1039  0.6875  1.4871  1.7622
 bndfp:  kpt 191 of 413, k=  -0.12500  0.20833  0.45833
 -0.4092 -0.1351 -0.0776 -0.0539 -0.0162  0.0073  1.0953  1.4959  1.6974
 bndfp:  kpt 191 of 413, k=  -0.12500  0.20833  0.45833
 -0.4086 -0.0865 -0.0256  0.0018  0.0390  0.0624  1.1170  1.5116  1.7087
 bndfp:  kpt 201 of 413, k=  -0.54167  0.62500  0.87500
 -0.3436 -0.1295 -0.0886 -0.0588 -0.0193  0.0385  0.8949  1.3172  1.8065
 bndfp:  kpt 201 of 413, k=  -0.54167  0.62500  0.87500
 -0.3388 -0.0846 -0.0375 -0.0046  0.0375  0.0931  0.9208  1.3330  1.8164
 bndfp:  kpt 211 of 413, k=  -0.25000  0.33333  0.66667
 -0.2434 -0.1481 -0.0973 -0.0529  0.0180  0.1316  0.6088  1.3168  1.4417
 bndfp:  kpt 211 of 413, k=  -0.25000  0.33333  0.66667
 -0.2190 -0.1105 -0.0490  0.0017  0.0758  0.1741  0.6437  1.3320  1.4541
 bndfp:  kpt 221 of 413, k=  -0.04167  0.12500  0.54167
 -0.3817 -0.1540 -0.0724 -0.0322 -0.0170  0.0077  1.2484  1.3714  1.4948
 bndfp:  kpt 221 of 413, k=  -0.04167  0.12500  0.54167
 -0.3797 -0.1048 -0.0229  0.0253  0.0379  0.0625  1.2660  1.3863  1.5095
 bndfp:  kpt 231 of 413, k=  -0.45833  0.54167  0.95833
 -0.2645 -0.1355 -0.1146 -0.0434 -0.0209  0.1021  0.8501  1.0298  1.6535
 bndfp:  kpt 231 of 413, k=  -0.45833  0.54167  0.95833
 -0.2534 -0.0900 -0.0656  0.0102  0.0377  0.1552  0.8780  1.0478  1.6678
 bndfp:  kpt 241 of 413, k=  -0.33333  0.41667  0.91667
 -0.1912 -0.1629 -0.1154 -0.0354  0.0019  0.2828  0.5843  0.9272  1.4464
 bndfp:  kpt 241 of 413, k=  -0.33333  0.41667  0.91667
 -0.1598 -0.1212 -0.0700  0.0186  0.0606  0.3162  0.6189  0.9477  1.4602
 bndfp:  kpt 251 of 413, k=  -0.29167  0.37500  0.95833
 -0.1811 -0.1789 -0.1094 -0.0255  0.0051  0.3676  0.5586  0.8352  1.3772
 bndfp:  kpt 251 of 413, k=  -0.29167  0.37500  0.95833
 -0.1476 -0.1346 -0.0647  0.0280  0.0644  0.3971  0.5914  0.8588  1.3896
 bndfp:  kpt 261 of 413, k=  -0.04167  0.12500  0.87500
 -0.2593 -0.1991 -0.0128  0.0113  0.0215  0.1742  0.8210  1.0358  1.1446
 bndfp:  kpt 261 of 413, k=  -0.04167  0.12500  0.87500
 -0.2268 -0.1529  0.0360  0.0654  0.0814  0.1955  0.8504  1.0473  1.1520
 bndfp:  kpt 271 of 413, k=  -0.08333  0.25000  0.41667
 -0.4281 -0.1289 -0.0763 -0.0550 -0.0288  0.0089  1.1446  1.4220  1.8270
 bndfp:  kpt 271 of 413, k=  -0.08333  0.25000  0.41667
 -0.4282 -0.0798 -0.0241  0.0009  0.0261  0.0638  1.1656  1.4376  1.8369
 bndfp:  kpt 281 of 413, k=  -0.12500  0.29167  0.54167
 -0.3268 -0.1434 -0.0931 -0.0483 -0.0056  0.0393  0.9127  1.3276  1.5702
 bndfp:  kpt 281 of 413, k=  -0.12500  0.29167  0.54167
 -0.3209 -0.0975 -0.0424  0.0067  0.0500  0.0928  0.9380  1.3429  1.5846
 bndfp:  kpt 291 of 413, k=  -0.54167  0.70833  0.95833
 -0.3875 -0.1347 -0.0838 -0.0462 -0.0321  0.0231  1.1138  1.2647  1.7601
 bndfp:  kpt 291 of 413, k=  -0.54167  0.70833  0.95833
 -0.3862 -0.0861 -0.0324  0.0096  0.0230  0.0779  1.1354  1.2792  1.7720
 bndfp:  kpt 301 of 413, k=  -0.33333  0.50000  0.83333
 -0.2137 -0.1488 -0.1087 -0.0481  0.0069  0.2141  0.5556  1.1276  1.4129
 bndfp:  kpt 301 of 413, k=  -0.33333  0.50000  0.83333
 -0.1834 -0.1110 -0.0617  0.0061  0.0649  0.2515  0.5926  1.1450  1.4269
 bndfp:  kpt 311 of 413, k=  -0.20833  0.37500  0.79167
 -0.2022 -0.1699 -0.1049 -0.0341  0.0157  0.2788  0.5486  1.0729  1.2595
 bndfp:  kpt 311 of 413, k=  -0.20833  0.37500  0.79167
 -0.1670 -0.1315 -0.0590  0.0203  0.0742  0.3114  0.5812  1.0916  1.2733
 bndfp:  kpt 321 of 413, k=  -0.16667  0.33333  0.83333
 -0.2090 -0.1755 -0.0974 -0.0210  0.0172  0.2977  0.6213  0.9651  1.1960
 bndfp:  kpt 321 of 413, k=  -0.16667  0.33333  0.83333
 -0.1741 -0.1335 -0.0527  0.0331  0.0761  0.3290  0.6493  0.9863  1.2090
 bndfp:  kpt 331 of 413, k=  -0.20833  0.37500  0.95833
 -0.1946 -0.1754 -0.1037 -0.0195  0.0127  0.4605  0.5147  0.8107  1.2155
 bndfp:  kpt 331 of 413, k=  -0.20833  0.37500  0.95833
 -0.1566 -0.1333 -0.0589  0.0337  0.0723  0.4880  0.5438  0.8359  1.2287
 bndfp:  kpt 341 of 413, k=  0.00000  0.16667  1.00000
 -0.2526 -0.2005 -0.0160  0.0103  0.0259  0.2696  0.7199  1.0051  1.0574
 bndfp:  kpt 341 of 413, k=  0.00000  0.16667  1.00000
 -0.2172 -0.1547  0.0345  0.0640  0.0870  0.2834  0.7560  1.0187  1.0641
 bndfp:  kpt 351 of 413, k=  -0.08333  0.33333  0.66667
 -0.2491 -0.1578 -0.1102 -0.0281  0.0021  0.1171  0.8295  1.1186  1.2990
 bndfp:  kpt 351 of 413, k=  -0.08333  0.33333  0.66667
 -0.2324 -0.1128 -0.0625  0.0259  0.0594  0.1649  0.8540  1.1358  1.3138
 bndfp:  kpt 361 of 413, k=  -0.04167  0.29167  0.70833
 -0.2479 -0.1675 -0.1012 -0.0133  0.0032  0.1242  0.9087  1.0628  1.1891
 bndfp:  kpt 361 of 413, k=  -0.04167  0.29167  0.70833
 -0.2282 -0.1213 -0.0557  0.0397  0.0611  0.1697  0.9309  1.0799  1.2049
 bndfp:  kpt 371 of 413, k=  -0.08333  0.33333  0.83333
 -0.2169 -0.1732 -0.0973 -0.0124  0.0170  0.2796  0.7371  0.9451  1.0532
 bndfp:  kpt 371 of 413, k=  -0.08333  0.33333  0.83333
 -0.1834 -0.1288 -0.0532  0.0409  0.0764  0.3113  0.7619  0.9679  1.0654
 bndfp:  kpt 381 of 413, k=  0.00000  0.25000  0.91667
 -0.2378 -0.1894 -0.0552  0.0019  0.0235  0.2949  0.7843  0.9095  0.9983
 bndfp:  kpt 381 of 413, k=  0.00000  0.25000  0.91667
 -0.2026 -0.1439 -0.0087  0.0551  0.0843  0.3178  0.8150  0.9279  1.0080
 bndfp:  kpt 391 of 413, k=  -0.04167  0.37500  0.79167
 -0.2107 -0.1630 -0.1163 -0.0159  0.0119  0.2576  0.7662  0.9488  1.0585
 bndfp:  kpt 391 of 413, k=  -0.04167  0.37500  0.79167
 -0.1799 -0.1183 -0.0724  0.0369  0.0713  0.2928  0.7914  0.9677  1.0750
 bndfp:  kpt 401 of 413, k=  -0.16667  0.50000  1.00000
 -0.1915 -0.1482 -0.1365 -0.0229  0.0163  0.4079  0.6286  0.8304  1.0455
 bndfp:  kpt 401 of 413, k=  -0.16667  0.50000  1.00000
 -0.1550 -0.1051 -0.0923  0.0301  0.0761  0.4361  0.6567  0.8532  1.0619
 bndfp:  kpt 411 of 413, k=  -0.08333  0.50000  1.00000
 -0.1985 -0.1408 -0.1375 -0.0221  0.0233  0.5188  0.6002  0.8036  0.9106
 bndfp:  kpt 411 of 413, k=  -0.08333  0.50000  1.00000
 -0.1605 -0.0973 -0.0936  0.0307  0.0840  0.5458  0.6270  0.8302  0.9274
 BNDFP:  evecs saved in file evec

 BZWTS : --- Tetrahedron Integration ---
 BZINTS: Fermi energy:      0.065080;  10.000000 electrons;  D(Ef):   29.644
         Sum occ. bands:   -0.9327044  incl. Bloechl correction:   -0.000408
         Mag. moment:       0.757262

 Saved qp weights ...

 Harris energy:
 sumev=       -0.932704  val*vef=    -151.345482   sumtv=     150.412777
 sumec=        0.000000  cor*vef=       0.000000   ttcor=    2934.175594
 rhoeps=    -121.970989     utot=   -5997.641462    ehar=   -3035.024080

 Make and renormalize projectors ...

 Find chemical potential mu ...
 Seek mu for 10 electrons ... 10.00178 electrons at Ef0=0.06508.  D(Ef0)=32.925
 getVnew: v1=0e0 N1=1.78e-3  v2=5.406e-5 N2=1.473e-6  bracket=F  est=5.41132e-5
 getVnew: v1=5.406e-5 N1=1.473e-6  v2=5.411e-5 N2=-1.459e-7  bracket=T  est=5.41088e-5
 getVnew: v1=5.411e-5 N1=-1.459e-7  v2=5.411e-5 N2=0e0  bracket=F  est=5.41088e-5
 mu = 0.065025 = Ef0-0.000054.  Deviation from neutrality = 0e0  D(mu)=32.870
 Electron charge:  10.000000   moment:  0.753440 spin resolved DOS:    3.643  29.227

 Make DOS for 201 points from 8 bands in window (-0.2,0.2) = mu + (-0.2,0.2)  

 *** Warning : You are calculating DOS without Pade interpolation!
 ... finished energy point 1 (5 sec)
 ... finished energy point 11 (53 sec)
 ... finished energy point 21 (53 sec)
 ... finished energy point 31 (53 sec)
 ... finished energy point 41 (53 sec)
 ... finished energy point 51 (53 sec)
 ... finished energy point 61 (53 sec)
 ... finished energy point 71 (53 sec)
 ... finished energy point 81 (53 sec)
 ... finished energy point 91 (53 sec)
 ... finished energy point 101 (53 sec)
 ... finished energy point 111 (53 sec)
 ... finished energy point 121 (53 sec)
 ... finished energy point 131 (53 sec)
 ... finished energy point 141 (53 sec)
 ... finished energy point 151 (53 sec)
 ... finished energy point 161 (53 sec)
 ... finished energy point 171 (53 sec)
 ... finished energy point 181 (53 sec)
 ... finished energy point 191 (53 sec)
 ... finished energy point 201 (53 sec)
 Exit 0 done writing total DOS 
 CPU time:  1391.504s   Wall clock 1391.546s  at  18:12:51 22.10.2018  on  localhost.localdomai
rdcmd:  cp dos.ni doslmfdmft.ni
