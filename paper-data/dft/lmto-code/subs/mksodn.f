      subroutine mksodn(nbas,nl,nsp,ldim,ipc,indxsh,pp,sop,socscl,sod)
C- Spin-orbit coupling matrix elements in downfolding order
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   pp    :potential parameters (atomsr.f)
Ci   sop   :spin-orbit parameters (atomsr.f)
Co Outputs
Co   sod(*,*,*,i): arrays for 1-, 2-, 3- center terms, ++,-- blocks,
Co                 downfolding order i = orbital index
Co   sod(*,1..3,isp,i): arrays for 1-, 2-, 3- center terms, ++,-- blocks
Co   sod(*,1,3,i):    for 1- center terms, L- S+ block
Co   sod(*,2..3,3,i): for 2- center terms, L- S+ block
Co   sod(*,4,3,i):    for 3- center terms, L- S+ block
Cl Local variables
Cr Remarks
Cu Updates
Cu   17 Dec 15 Repackaged from mksod to enable use with GF package
Cu             mksodn does not require pph or soph
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nl,nbas,ldim,nsp,indxsh(*),ipc(nbas)
      double precision sod(ldim,4,3),pp(6,0:nl-1,nsp,*),socscl(*)
      double precision sop(0:nl-1,nsp,nsp,9,*)
C ... Local parameters
      integer is,ib,lmr,ilm,l,m,i,nlm,ic
      double precision ovl(0:nl-1,2),CmEnu(0:nl-1,2),sodi(nl*nl,4,3)
      procedure(real(8)) :: oalpha

      call rx('SO coupling not installed')

      end
      subroutine mksod1(nl,ndim,CmEnu,ovl,sop,sod)
C- Diagonal arrays for one-, two- and three- center terms for one site in the SO ham
C ---------------------------------------------------------------------
Ci Inputs
Ci   nl    :(global maximum l) + 1
Ci   ndim  :dimensions sod.  Must be at least as large as nl*nl
Ci   CmEnu :Difference between band center C and linearization energy
Ci   sop   :spin-orbit parameters for this class (atomsr.f)
Co         : sop(l,is1,is2,i) matrix elts between spins is1 and is2
Co         : (i=1) <phi so  phi> (i=2) <phi so  dot> (i=3) <dot  so dot>
Ci
Co Outputs
Co   sod   : Terms needed to 1C, 2C, 3C parts of make matrix elements of <chi SO chi>
Co   sod(ilm,1..3,isp) for 1-, 2-, 3- center terms, ++,-- blocks; see Remarks
Co   sod(ilm,4,1:2):   not used
Co   sod(ilm,1,3):     for 1- center terms L- S+ block
Co   sod(ilm,2,3):     sqrt((l+m+1)*(l-m))/2 [ <dot SO phi> + o <phi SO phi>] for 2-C (h L-S+)
Co   sod(ilm,3,3):     sqrt((l+m)*(l-m+1))/2 [ <phi SO dot> + o <phi SO phi>] for 2-C (L-S+ h)
Co   sod(ilm,4,3):     for 3- center terms, L- S+ block
Cr Remarks
Cr   Makes diagonal matrices belonging to the 1-, 2- and 3-center
Cr   terms in the matrix elements of the SO hamiltonian.  Its form is
Cr   <SO_{RL,RL'}> = <chi_{RL}| I L . S | chi_{RL'}> = I_l L.S
Cr   where L.S = Lz.Sz + .5(L+.S- + L-.S+) and I = (2/(c^2) dV/dr*(1/r).
Cr   Inside the the sphere, a basis function is a linear combination
Cr   of phi's and dot's (including downfolded orbitals):
Cr     | chi> = | phi_L> + | phidot_L> h_LL + | phi_I> (h)_IL
Cr            = | phi_L>(1+oh_LL) + |dot_L> h_LL + | phi_I> (h)_IL
Cr   The first form uses phidot = phidot^alpha; the second form uses
Cr     phidot^alpha = phidot^gamma + o phi; and calls 'dot'=phidot^gamma
Cr
Cr   This subroutine combines radial matrix elements I between phi's
Cr   and dot's (phidot^gamma) and the L.S operator into diagonal
Cr   matrices associated with the 1-, 2- and 3-center terms comprising
Cr   matrix elements (H_so)_(i,j).
Cr   For the (L-S+) block, the SO coupling uses:
Cr     L+ |lm> = sqrt((l-m)*(l+m+1)) |l,m+1>
Cr     L+ |lm-1> = sqrt((l-m-1)*(l+m)) |l,m>
Cr     L- |lm> = sqrt((l+m)*(l-m+1)) |l,m-1>
Cr     L- |lm+1> = sqrt((l+m+1)*(l-m)) |l,m>
Cr   Then
Cr     <phi1_lm | L-.S+ | phi2_lm+1> = F/2 < phi1_l | I | phi2_l>
Cr   where F = sqrt((l-m)*(l+m+1))
Cr   Also
Cr     <phi1_lm-1 | L-.S+ | phi2_lm> = F'/2 < phi1_l | I | phi2_l>
Cr   where F' = sqrt((l+m)*(l-m+1))
Cr   L-.S+ couples l,m only to l,m+1 (2l elements for each l).  The
Cr   following evaluates (2l+1) elements foreach l, exploiting the
Cr   fact that <lm | L-S+ | lm+1> is zero for m=l because F=0, and
Cr   similarly <lm-1 | L-S+ | lm> is zero for m=-l.
Cr
Cr   Matrix elements have 1-, 2- and 3- center terms composed of
Cr   matrix elements of partial waves inside augmentation spheres
Cr     xi1 = <phi|SO|phi>  xi2 = <phi|SO|dot>  xi3 = <dot|SO|dot>
Cr   <chi | SO | chi> =
Cr     (1+ho)<phi_L| + h<dot_L| SO  | phi_L>(1+oh + |dot_L> h =
Cr     (1+ho) xi1 (1+oh) + (1+ho) xi2 h + cc + h xi3 h
Cr   Analyze 1C contribution to spin structure of L.S, when h is spin diagonal and o=0.
Cr      <phi_lm | ILzSz | phi_lm+1>_ss   Spin diagonal part, (1,1) and 22 blocks
Cr      <phi_lm | IL-S+ | phi_lm+1>_12   Spin off-diagonal part, (1,2) block
Cr   Analyze 2C contribution to spin structure of L.S, when h is diagonal and o=0.
Cr   Spin diagonal part
Cr      h_ss <dot_lm | ILzSz | phi_lm>_ss    and  <phi_lm | ILzSz | dot_lm>_ss h_ss
Cr   Spin off-diagonal part :  L-S+ in (1,2) and (2,1) blocks
Cr      h_11 <dot_lm | IL-S+ | phi_lm+1>_12  and  <phi_lm | IL-S+ | dot_lm+1>_21 h_11
Cr      h_22 <dot_lm | IL-S+ | phi_lm+1>_21  and  <phi_lm | IL-S+ | dot_lm+1>_12 h_22
Cr   Radial, angular parts factor.  Use <dot_l I phi_l>_ss' = <phi_l I dot_l>s's .  Then
Cr      <dot_lm IL-S+ phi_lm+1>_12 = <phi_lm IL-S+ dot_lm+1>_21 = sod(:,2,3)
Cr      <phi_lm-1 IL-S+ dot_lm>_12 = sod(:,3,3)
Cr
Cr   Divide into 1-, 2-, 3- center terms
Cr     w1      where w1 = xi1
Cr     h w2    where w2 = (o xi1 + xi2)
Cr     h w3 h  where w3 = (o^2 xi1 + 2 o xi2 + xi3)
Cr
Cr   To work directly with  D = srdel S srdel, rather than
Cr   h = C-enu + D, the diagonal parts w1 connecting C-enu
Cr   in the one-, two-, three-center terms are reshuffled.  Thus
Cr     w1 + w2 (C-e+D) + h.c. + (C-e+D) w3 (C-e+D)  becomes
Cr     w1+2*w2*(C-e)+w3(C-e)**2  +  (w2+(C-e)w3) D + h.c.  +  D w3 D.
Cr   Terms for the spin-diagonal parts proportional to D^0, D^1 and D^2
Cr   are stored in sod(:,1:3,isp)
Cu Updates
Cu   16 Aug 13 Bug fix initializing some elements of sod
Cu   08 May 13 Initialize sod (some elements were never referenced)
Cu   05 Jan 04 Fix upper limit of loop in downfolding case
Cu   14 Feb 03 Dimensions of sod, soph altered to accomodate B-field
C ---------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nl,ndim
      integer, parameter :: nsp=2
      double precision sop(0:nl-1,nsp,nsp,9),sod(ndim,4,3)
      double precision ovl(0:nl-1,2),CmEnu(0:nl-1,2)

C ... Local parameters
!     If we are exiting with 'rx' below thre is no point in local variables.
!       integer l,isp,m,ilm,nlm,mi
!       double precision sp,sz,flmp,flmm
!       parameter (sp=.5d0)
!       flmp(mi) = dsqrt(dble((l+mi+1)*(l-mi))) ! L+ |lm> = sqrt((l-m)*(l+m+1)) |l,m+1>
!       flmm(mi) = dsqrt(dble((l+mi)*(l-mi+1))) ! L- |lm> = sqrt((l+m)*(l-m+1)) |l,m-1>

      call rx('SO coupling not installed')

      end
