! Text following "!" are comments
! Input lines consist of "keyword value(s)"
! New section begins with <tagname>
!Verbose    41  ! 0-->default; 100-->debug
GWversion   12  ! 0-->Sep12 compatible.
                ! The following are active (take any combination)
                !   1 include LFC in offset-Gamma treatment (original version)
                !   2 improved LFC; see J. Phys. Soc. Jpn. 83, 094711 (2014))
                !  10 improved unocc counting in making G*W
!Q0P_Choice 1   ! 1(default):Near q->0 limit 2:1/q^2 average in Gamma region
!CoreOrth  off  ! off --> Do not orthogonalize core to valence (default)
                ! on  --> Orthogonalize cores to valence (may give strange core functions!)
!multitet 2 2 2 ! tetrahedra divided into micro-tetrahedra
!EXonly   .15   ! for exchange-only calculations
!EIBZmode off   ! turn off to suppress use of symmetrization in calculating polarization
!TimeReversal off ! when time-reversal symmetry is broken
KeepEigen  on   ! keep eigenfunctions in memory
KeepPPOVL  on   ! keep PPOVL in memory
!Chi_RegQbz on  ! Offset Gamma point mesh for chi.  on => no offset
BZmesh     1    ! Offset Gamma point mesh for Sigma=iGW
WgtQ0P     0.01 ! Weight used when BZmesh is 2
NormChk    0    ! 1,2 writes norm check files (diagonal or full)
n1n2n3  4  4  4 ! for GW BZ mesh
QpGcut_psi 2.7  ! |q+G| cutoff for eigenfunction
QpGcut_cou 2.2  ! |q+G| cutoff for coulomb int.
unit_2pioa off  ! off --> units of 2 preceding Gcut are a.u.; on--> units are 2*pi/alat
alpha_OffG 1    ! (a.u.) parameter (auxiliary function, offset-Gamma method)
!nband_chi0 999 ! nband cutoff for chi0  (Optional)
!emax_chi0 999. ! emax  cutoff for chi0, Ry  (Optional)
!nband_sigm 999 ! nband cutoff for Sigma (Optional)
emax_sigm 2.5   ! Energy cutoff for Sigma, Ry (Optional)
dw       0.01   ! mesh spacing along Real axis (Hartree)
omg_c    0.1    ! mesh spacing increases linearly with w: becomes 2*dw at omg_c
iSigMode 3      ! QSGW mode switch (QSGW only)
niw      6      ! # freq. on Im axis; used for integration to make Sigma_c
delta    -1e-4  ! delta-function broadening for calc. x0, a.u.. delta<0 => tetrahedron
deltaw   0.02   ! width in finite diff for sigma energy derivative, a.u.
esmr     3e-3   ! Broadening in the poles of G(LDA) (hsfp0)
                ! Change esmr for metals: see DOSACC* --- especially around Ef
GaussSmear on   ! on  => broadening of poles in G(LDA) by Gaussian
                ! off => broadening of poles by a rectangle
!mixbeta  .25   ! mixing of input, output sigma for self-consistency

<PRODUCT_BASIS>   ! Product basis block
  tolerance = minimum eigenvalue in PB overlap to reduce linear dependency
  3e-4 3e-4 1e-3
  lcutmx(atom) = l-cutoff for the product basis
  4  4  3  3
  atom   l  nnvv  nnc ! nnvv: num. radial functions (valence) for augmentation-waves. nnc = num. for core.
    1    0    2    2
    1    1    2    1
    1    2    2    0
    1    3    2    0
    1    4    2    0
    2    0    2    2
    2    1    2    1
    2    2    2    0
    2    3    2    0
    2    4    2    0
    3    0    2    0
    3    1    2    0
    3    2    2    0
    3    3    2    0
    3    4    2    0
    4    0    2    0
    4    1    2    0
    4    2    2    0
    4    3    2    0
    4    4    2    0
  atom   l    n  occ  unocc  :Valence(1=yes, 0=no)
    1    0    1    1    1   ! 3S_p *
    1    0    2    0    0   ! 3S_d
    1    1    1    1    1   ! 3P_p
    1    1    2    0    0   ! 3P_d
    1    2    1    1    1   ! 3D_p
    1    2    2    0    0   ! 3D_d
    1    3    1    0    1   ! 4F_p
    1    3    2    0    0   ! 4F_d
    1    4    1    0    0   ! 5g_p
    1    4    2    0    0   ! 5g_d
    2    0    1    1    1   ! 3S_p *
    2    0    2    0    0   ! 3S_d
    2    1    1    1    1   ! 3P_p
    2    1    2    0    0   ! 3P_d
    2    2    1    1    1   ! 3D_p
    2    2    2    0    0   ! 3D_d
    2    3    1    0    1   ! 4F_p
    2    3    2    0    0   ! 4F_d
    2    4    1    0    0   ! 5g_p
    2    4    2    0    0   ! 5g_d
    3    0    1    1    1   ! 1S_p *
    3    0    2    0    0   ! 1S_d
    3    1    1    1    1   ! 2P_p
    3    1    2    0    0   ! 2P_d
    3    2    1    0    1   ! 3D_p
    3    2    2    0    0   ! 3D_d
    3    3    1    0    0   ! 4f_p
    3    3    2    0    0   ! 4f_d
    3    4    1    0    0   ! 5g_p
    3    4    2    0    0   ! 5g_d
    4    0    1    1    1   ! 1S_p *
    4    0    2    0    0   ! 1S_d
    4    1    1    1    1   ! 2P_p
    4    1    2    0    0   ! 2P_d
    4    2    1    0    1   ! 3D_p
    4    2    2    0    0   ! 3D_d
    4    3    1    0    0   ! 4f_p
    4    3    2    0    0   ! 4f_d
    4    4    1    0    0   ! 5g_p
    4    4    2    0    0   ! 5g_d
  atom   l    n  occ unocc   ForX0 ForSxc :CoreState(1=yes, 0=no)
    1    0    1    0    0      0    0    ! 1S *
    1    0    2    0    0      0    0    ! 2S
    1    1    1    0    0      0    0    ! 2P
    2    0    1    0    0      0    0    ! 1S *
    2    0    2    0    0      0    0    ! 2S
    2    1    1    0    0      0    0    ! 2P
</PRODUCT_BASIS>

<QPNT>   ! Specify particular k-points and bands (one-shot mode)
 --- Specify qp and band indices at which to evaluate Sigma
 
*** Sigma at all q -->1; to specify q -->0.  Second arg : up only -->1, otherwise 0
  0  0
*** no. states and list of band indices to make Sigma and QP energies
  8
  1 2 3 4 5 6 7 8
*** q-points (must belong to mesh of points in BZ).
  3
  1     0.0000000000000000     0.0000000000000000     0.0000000000000000
  2    -0.2500000000000000     0.2500000000000000     0.2500000000000000
  3    -0.5000000000000000     0.5000000000000000     0.5000000000000000
  4     0.0000000000000000     0.0000000000000000     0.5000000000000000
  5    -0.2500000000000000     0.2500000000000000     0.7500000000000000
  6    -0.5000000000000000     0.5000000000000000     1.0000000000000000
  7     0.0000000000000000     0.0000000000000000     1.0000000000000000
  8     0.0000000000000000     0.5000000000000000     1.0000000000000000
</QPNT>

QforEPSIBZ off
<QforEPS>
0d0 0d0 0.015d0
</QforEPS>
EPSrange  1    !(Ry) [0,EPSrange] for dielectric function
EPSdw     0.05 !(Ry) energy mesh  for dielectric function
