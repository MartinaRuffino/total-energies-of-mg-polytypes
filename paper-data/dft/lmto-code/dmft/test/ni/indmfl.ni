0.1 1.2 1 2                # hybridization Emin and Emax, measured from FS, renormalize for interstitials, projection type
1 0.0025 0.0025 600 -3.000000 1.000000  # matsubara, broadening-corr, broadening-noncorr, nomega, omega_min, omega_max (in eV)
1                                     # number of correlated atoms
1     1   0                           # iatom, nL, locrot
  2   2   1                           # L, qsplit, cix
#================ # Siginds and crystal-field transformations for correlated orbitals ================
1     5  5       # Number of independent kcix blocks, max dimension, max num-independent-components
1     5  5       # cix-num, dimension, num-independent-components
#---------------- # Independent components are --------------
'x^2-y^2' 'z^2' 'xz' 'yz' 'xy'
#---------------- # Sigind follows --------------------------
1 0 0 0 0 
0 2 0 0 0 
0 0 3 0 0 
0 0 0 4 0 
0 0 0 0 5 
#---------------- # Transformation matrix follows -----------
 0.70710679 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.70710679 0.00000000
 0.00000000 0.00000000   0.00000000 0.00000000   1.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000
 0.00000000 0.00000000   0.70710679 0.00000000   0.00000000 0.00000000  -0.70710679 0.00000000   0.00000000 0.00000000
 0.00000000 0.00000000   0.00000000 0.70710679   0.00000000 0.00000000   0.00000000 0.70710679   0.00000000 0.00000000
-0.70710679 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.00000000 0.00000000   0.70710679 0.00000000
