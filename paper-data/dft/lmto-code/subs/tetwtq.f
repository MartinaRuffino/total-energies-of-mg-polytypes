C#define EXPAND_SORTEA
      subroutine tetwtq(job,npm,ncc,linme,esciss,qlat,fac,
     .  efermi,eband1,eband2,nctot,ecore,
     .  ntetf,nqbzw,nband,nqbz,idtetf,qbzw,ib1bz,
     .  frhis,nwhis,nwgtx,ibjb,nhwtot,ihw,nhw,jhw,
     x  q,iq,isp1,isp2,nqibz,
     o  iwgt,nwgt,demin,demax,
     o  whw)
C- Weights for dielectric function by tetrahedron method.
C ----------------------------------------------------------------------
Ci Inputs
Ci   job   :0: return iwgt,nwgt,demin,demax only
Ci         :1: make (ib,jb) and k-resolved weights whw
Ci   npm   :1 : use time reversal symmetry
Ci         :2 : no time reversal symmetry
Ci   ncc   :number of cores for jb in (ib,jb) pair:
Ci         :Should be 0 if time-reversal, nctot otherwise
Ci   eband1:eigenvalues at k, in Ry units (see Bugs below)
Ci   eband2:eigenvalues at k+q (see Bugs below)
Ci   linme :F, matrix elements assumed to be k-independent
Ci         :T, matrix elements assumed to vary linearly with k
Ci   esciss:energy for scissors operator, unoccupied states
Ci   qlat  :reciprocal lattice vectors
Ci   efermi:Fermi energy  (see Bugs below)
Ci   ntetf :number of tetrahedra
Ci   nqbzw :dimensions qbzw,ib1bz. Should be (n1+1)*(n2+1)*(n3+1)
Ci   nband :number of bands
Ci   nqbz  :number of k-points in the 1st BZ
Ci   nctot :number of cores
Ci   ecore :core eigenvalues. (see Bugs below)
Ci         :Not used if nctot=0
Ci   idtetf:idtetf(1:4,i) points to the 4 k-points defining the ith
Ci         :tetrahedron.  Note: idtetf points to the padded qbzw
Ci   qbzw  :k-points on the regular mesh, padded to repeat points
Ci         :on the boundary lines (tetfbz.f)
Ci   ib1bz :maps the points in the padded qp list to the original list
Ci   iwgt  :T: ib(occ) jb(unocc) pair makes nonzero contribution
Ci         :Output when job=0
Ci         :Input  when job=1
Ci   nwgt  :number of ib(occupied) jb(unoccupied) pairs for this kpoint
Ci         :Output when job=0
Ci         :Input  when job=1
Ci   demin :minimum excitation energy for ib(occupied) jb(unoccupied) pair
Ci         :Output when job=0 (see Bugs below)
Ci         :Input  when job=1
Ci   demax :maximum excitation energy for ib(occupied) jb(unoccupied) pair
Ci         :Output when job=0 (see Bugs below)
Ci         :Input  when job=1
Ci   frhis :frhis(1:nwhis+1) = Energies in histogram (See Bugs below)
Ci         :Range for bin i is frhis(i) to frhis(i+1)
Ci         :Not used when job=0
Ci         :Note!! frhis is in Hartree, while eband are in Ry
Ci   nwhis :number of energies in histogram
Ci         :Not used when job=0
Ci   nwgtx :leading dimension of ihw,jhw,nhw
Ci         :should be maxval(nwgt)
Ci         :Not used when job=0
Ci   ibjb  :ibjb(nctot+nband,nband,nqbz) = ibjb index for given ib jb k.
Ci   nhwtot:dimension of whw
Ci   ihw   :ihw(ibjb,k) = index to first histogram bin encompassing
Ci                        (demin,demax) for a given ibjb pair and k
Ci         :Not used when job=0
Ci   nhw   :nhw(ibjb,k) = number of histograms bins encompass
Ci                        (demin,demax) for a given ibjb pair and k
Ci         :Not used when job=0
Ci   jhw   :jhw(ibjb,kx) = pointer to whw
Ci         :Not used when job=0
Ci   q     :(mtet mode) q-vector separating occ, unocc states
Ci         :True q = 2*pi*q(1:3)/alat
Ci         :Not used now
Ci   iq    :(mtet mode) used to locate q vector in 1st BZ
Ci         :Not used now
Ci   isp1  :(mtet mode) spin index for occ   state (evals from disk)
Ci         :Not used now
Ci   isp2  :(mtet mode) spin index for unocc state (evals from disk)
Ci         :Not used now
Ci   nqibz :(mtet mode) # kpoints read from disk?
Ci         :Not used now
Co Outputs
Co   whw   :whw(i:j) histogram weights in bins i:j for given ib,jb,kx
Co         :i = jhw(ibjb,kx)
Co         :j = jhw(ibjb,kx)+nhw(ibjb),kx)-1
Co         :ibjb = ibjb(ib,jb,kx)
Co         :Generated when job=1
Cl Local variables
Cl         :
Cr Remarks
Cr  *tetwtq calculates by the tetrahedron method this integral:
Cr     whw(i,jb,ib,ik) =
Cr        \int_om_i^om_i+1 d\omg \times
Cr        \int d^3k f(e(k))(1-f(e(q+k)))\delta(omg-e2_jb(q+k)+e1_ib(k))
Cr   on a mesh of k-points ik for a sequence of histograms i=1,2,3,...
Cr   Here:
Cr     om_i   = lower energy bound of histogram, frhis(i)
Cr     om_i+1 = upper energy bound of histogram, frhis(i+1)
Cr     f(E)   = Fermi distribution function for T=0 (efermi is input)
Cr     e2_jb  = jbth state eband2 (unoccupied)
Cr     e2_ib  = ibth state eband1 (occupied)
Cr
Cr   Note: whw is stored in a compressed format; see Outputs
Cr   The total joint DOS in the interval (om_i,om_i+1) is
Cr     JNOS(i) = sum_(ib,jb,ik) whw(i,ib,jb,ik)
Cr
Cr   JNOS is not calculated here, but in tetwtt after call to tetwtq.
Cr
Cr   This code is based on: J.Rath&A.J.Freeman PRB11, 2109 (1975)
Cr   lindtet6 calculates the weights for one tetrahedron.
Cr
Cr *Kotani's notes from tetwt5.F
Cr   The numbering of band index is not unique due to degeneracy.
Cr   It affects how to choose tetrahedron.  As a result, it affect on
Cr   the integration weight, and it might break the crystal symmetry.
Cr   So I add symmetrization at the last in this routine so as to recover the symmetry.
Cr
Cr *mtet mode: doubles the k-mesh in the integration of the delta-function
Cr              (mesh has 2x2x2 more tetrahedra than normal mode).
Cr   Normal  mode requires mtet=(/1,1,1/)
Cr   Doubled mode requires mtet=(/2,2,2/)
Cr   mtet mode is disabled at present
Cb Bugs
Cb   Input frhis is in Hartree units, while other input energies
Cb   (efermi, esciss, eband1, eband2, ecore) are in atomic Rydberg units.
Cb   demin, demax (internally generated) are in Hartree units.
Cu Updates
Cu   07 Sep 09  Adapted from Kotani's tetwt5.F
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical linme
      integer(4) :: job,nband,nqbz,nctot,ntetf,nqbzw,npm,ncc,
     .  idtetf(0:3,ntetf),ib1bz(nqbzw),nwgt(nqbz,npm)
      integer(4):: nwgtx,nhwtot,nwhis
      real(8) :: fac,efermi,esciss,qlat(3,3),qbzw(3,nqbzw),
     .  eband1(nband,nqbz),eband2(nband,nqbz),ecore(nctot)
      real(4) :: demax(nband+nctot,nband+ncc,nqbz,npm),
     .           demin(nband+nctot,nband+ncc,nqbz,npm)
      logical :: iwgt(nband+nctot,nband+ncc,nqbz,npm)
      real(8) :: whw(nhwtot),frhis(nwhis+1)
      integer(4):: ihw(nwgtx,nqbz,npm), ! omega pointer
     .             nhw(nwgtx,nqbz,npm), ! number of bins
     .             jhw(nwgtx,nqbz,npm), ! index to whw
     .             ibjb(nctot+nband,nband+ncc,nqbz,npm)
C     mtet mode
      integer(4)::isp1,isp2,iq,nqibz
      real(8) :: q(3)

      call rx('TETWTQ: OPTICS package not installed')

      end

