      subroutine mkmelt(job,isp1,isp2,qlat,q1,q2,q2p,cphiq,cphiqk,ppb,icore,ncore,
     .  symope,shftv,tr,nl,nbas,nclass,iclass,indexpbi,nnlmx,nnlmv,nnlmc,ipbfn,ntpba,mxnpba,naugm,
     .  npwmb,npbmb,it2,nt2,itlst,npm,it1,nt1,ipblo,ipbhi,it1lo,it1hi,it2lo,it2hi,nnc,lprt,
     o  zpsimb)
C- Matrix elements <psi(k+q,t') | psi(k,t) M >
C ----------------------------------------------------------------------
Ci Inputs
Ci  job    : 1s digit
Ci         : 0 : matrix elements for susceptibility
Ci         : 1 : matrix elements for GW.
Ci  isp1   : spin index for states at q1
Ci  isp2   : spin index for states at q2
ci         : isp1=isp2 except for spin-flip processes
Ci  qlat   : Reciprocal lattice vectors
Ci  q1     : Polarizability context: q1 = q+rk, where polarizability is calculated at q
Ci         : Self-energy context: q1 is q-point where sigma is calculated
Ci  q2     : Polarizability context: q where chi(q) is calculated
Ci         : Self-energy context: k-point in the full BZ entering into k-integral, Eq. 34, PRB76, 165106
Ci  q2p    : Polarizability context: should be identical to q2
Ci         : Self-energy context: the k-point in the IBZ that is symmetry-related to q2
Ci  cphiq  : Polarizability Context: (see Remarks in psi2b_v4)
Ci         : amount of partial wave contributing to unocc valence eigenfunction at k+q
Ci         : Self-energy context: (see Remarks in psi2b_v4)
Ci         : amount of partial wave contributing to eigenstate n or m at q (see Remarks)
Ci  cphiqk : Polarizability Context: (see Remarks in psi2b_v4)
Ci         : amount of partial wave contributing to occ valence eigenfunction at k
Ci         : Self-energy context: (see Remarks in psi2b_v4)
Ci         : amount of partial wave contributing to intermediate state at q-rk
Ci         : rk is analogous to k, but not restricted to the irreducible BZ.
Ci  ppb    : (real) matrix elements of partial waves <Phi(RLn) Phi(RL'nn') B(R,i)>
Ci         : ppb includes both core and valence.  Sequence :
Ci         : (core at 1), (val at 1) (core at 2), (val at 2) ...
Ci  icore  : icore(i,ic) index in ppb array for core state i in class ic
Ci  ncore  : ncore(ic) no. core states in class ic
Ci  symope : Polarizability context: symope should be the unit matrix
Ci         : Self-energy context: rotation part of space group (that takes k to kr?)
Ci  shftv  : Polarizability context: should be zero
Ci         : Self-energy context: translation part of space group
Cr  tr     : Polarizability context: not used.
Ci         : Self-energy context: translation tiat for this symop (see mptauof)
Ci         : tr is vector that must be added to bas(ib') to coincide with symope(ib)*bas(ib)+shftv.
Ci  nl     : dimensions icore
Ci  nbas   : number of sites in the basis
Ci         : OLD natom
Ci  nclass : number of classes; used only for dimensioning
Ci  iclass : the jth atom belongs to class iclass(j)
Ci indexpbi: indexpbi(i): for site i, points to the site which for product basis is stored
Ci         : OLD iatomp
Ci  nnlmx  : maximum number of partial waves at a site; dimensions ppb
Ci         : OLD nlnmx
Ci  nnlmv  : nnlmv(ic) total number of partial waves (l,n,m) for valence states for site of class ic
Ci         : OLD nlnmv
Ci  nnlmc  : nnlmc(ic) total number of partial waves (l,n,m) for core states for site of class ic
Ci         : OLD nlnmc
Ci  ipbfn  : ipbfn(i) = index to first product basis function for site i
Ci         : Product functions are strung together as a long vector.
Ci         : OLD: imdim
Ci  ntpba  : ntpba(ic) = number of product basis functions for class ic
Ci         : OLD mdim or nblocha
Ci  mxnpba : dimensions ppb : max number of product basis functions for any site
Ci         : OLD mdimx
Ci  naugm  : number of augmentation channels ... called ndima in sugw.f
Ci         : OLD nlmto
Ci  npbmb  : total number of local product basis functions within MTs
Ci         : OLD nbloch
Ci  npwmb  : number of G vectors in mixed basis
Ci          : (formerly ngc)
Ci  ngpmx  : max number of G vectors for eigenfunctions (?)
Ci         : Leading dimensions of geigen (see readgeig)
Ci  it2    : Polarizability context (jpm=1): it2 = index to 1st unocc state at k
Ci         : Exception: without time reversal 2nd pass (jpm=2): it2 = 1
Ci         : Self-energy context: it2 = 1 (Select specific states in 1-shot GW through itlst)
Ci  nt2    : Polarizability context: no. unoccupied states at q+k (polarizability context)
Ci         : Self-energy context: no. states at q for which sigm calculated (n or m in Eq 34)
Ci         : Note nt2 is NOT associated with q2, but with 2nd band index in zpsimb
Ci         : OLD: ntp0
Ci  itlst  : If caller passes itlst(1)=0, matrix elements are calculated for it2,it2+1,...it2+nt2
Ci         : this should always be the case in the Polarizability context
Ci         : Self-energy context (1-shot only): calculate matrix elements for states itlst(it2:nt2)
Ci         : OLD itq
Ci  npm    : 1 if time-reversal symmetry, 2 if not
Ci  ncore2 : Number of core states for matrix elements needed when TR symmetry is absent.
Ci         : It is used in the polarizability context:
Ci         : 0     => do not calculate <core(k+q,t') |psi (k,t) B(R,ibloch)>
Ci         : nctot => do calculate <core(k+q,t') |psi (k,t) B(R,ibloch)>
Ci         : ncore2 = 0 in the self-energy context
Ci  it1    : BC added for BSE calculation
Ci  nt1    : Polarizability context: no. occupied states at k, n in in Eq. 32, PRB76, 165106
Ci         : Self-Energy context: no. intermediate states at q-k, n' in Eq. 34, PRB76, 165106
Ci         : excluding core states; see nt
Ci         : Note nt1 is NOT associated with q1, but with 1st band index in zpsimb
Ci         : OLD: nbmax
Ci  ipblo  : lower bound in product basis index for which zpsipb is calculated
Ci  ipbhi  : upper bound in product basis index for which zpsipb is calculated
Ci  it1lo  : lower bound in first band index for which zpsipb is calculated
Ci  it1hi  : upper bound in first band index for which zpsipb is calculated
Ci  it2lo  : lower bound in 2nd band index for which zpsipb is calculated
Ci  it2hi  : upper bound in 2nd band index for which zpsipb is calculated
Ci  nnc    : dimensions icore
Ci  lprt   : T, print out timing information
Co Outputs
Co  zpsimb : matrix elements (see also Remarks)
Co         :*Polarizability context: <psi(q+k,t') | psi(k,t) M > .. Eq. 32, PRB76, 165106
Co         : Stored as zpsimb(I(q), it1(k), it2(q+qrk))
Co         : it1: band index for occupied states t at k (more precisely rk)
Co         : it2: band index for unoccupied t' states at q+k
Co         :*Self-energy context:  <psi(q,t') | psi(q-k,t) M(rk,I)> ... Eq. 34, PRB76, 165106
Co         : Stored as zpsimb(I(rk), it1(q), it2(q-qrk))
Co         : it1: band index for intermediate states t at qk = q - qrk
Co         : it2: band index for states t' at target q, for which sigma(t') will be made
Co         :*In either context:
Co         : I = index to mixed product basis at qrk (or written as rk), I = 1:nmb
Co         : nmb=npbmb+npwmb = # of mixed product basis
Co         :           npbmb = # of product basis (within MTs)
Co         :           npwmb = # of IPWs
Co         : In tabular form:
Co         : argument   polarizability-context     self-energy context
Co         :    I          product basis               product basis
Co         :   it1        occ states at q          intermediate states at q-k
Co         :   it2      unocc states at q+k     states for matrix elements at q
Cl Local
Cl  ncore1 : number of core states and offset to first valence state in zpsimb(:,nctot+1:,:)
Ci         : OLD: nctot
Cl  nt     : Polarizability context: no. occupied states at k, n in in Eq. 32, PRB76, 165106
Cl         : = nt1
Cl         : Self-energy context: no. intermediate states at q-k, n' in Eq. 34, PRB76, 165106
Cl         : = nctot + nt1
Cl  delq   : Polarizability context: rk where integral d(rk) [integrand(q,rk)] is calculated
Cl         : Self-energy context: q-rk
Cl expikt  : exp(ik.T)
Cr Remarks
Cr   zpsimb(ip,i1,i2) contains matrix elements of product basis ip with eigenstates i1 and i2.
Cr     ip is a product basis function, spanning both site-centered functions B and plane waves
Cr      The npbmb product functions precede the npwmb PW states
Cr     i1 can span both core states and valence states, the core states preceding the valence
Cr     i2 can span both core states and valence states (cores only for x0 without TR symmetry)
Cr   Ordering of states in zpsimb :
Cr     product basis index ip: site functions precede PWs
Cr     band indices i1,i2:     core states (if any) precede valence states
Cr
Cr   Mixed Basis consists of product basis B and IPW part.
Cr   Matrix element for B functions includes phase exp(ik.T)
Cr   Note that B(rot*k,r) = B(k,invrot*r)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical lprt
      integer :: job,isp1,isp2,mxnpba,naugm,npwmb,npbmb,nnlmx,it2,nt2,nbas,
     .  nclass,nl,nnc,ipblo,ipbhi,it1lo,it1hi,it2lo,it2hi,npm,it1,nt1
      integer :: icore(nl*nl*nnc,nclass),ncore(nclass),iclass(nclass),
     .  nnlmc(nclass),nnlmv(nclass),ntpba(nclass),itlst(*),ipbfn(nbas),indexpbi(nbas)
      real(8) :: ppb(nnlmx,nnlmx,mxnpba,nclass)
      real(8) :: qlat(3,3),q1(3),q2(3),q2p(3),symope(3,3),shftv(3),tr(3,nbas)
      complex(8),target :: cphiq(naugm,nt2),cphiqk(naugm,nt2)
!     complex(8) :: zpsimb(nmb,nt1,nt2)
      complex(8) :: zpsimb(ipblo:ipbhi,it2lo:it2hi,it1lo:it1hi)
C ... Dynamical local arrays
      integer, allocatable :: igv1(:,:),igv2(:,:)
      complex(8),allocatable :: geigq1(:,:),geigq2(:,:)
C ... Local
      integer ia,ngq1,ngq2,nadd(3),ncore1,ncore2,it,itp,ipb
      real(8) :: qlatinv(3,3),det,tpi,delq(3),qmkt(3),qt(3),qu1(3),qu2(3),qdiff(3),add(3)
      complex(8) :: expikt(nbas)
      complex(8), pointer :: cphiq2(:,:)
      complex(8),parameter :: img=(0d0,1d0)

C ... Initialization
      tpi = 8d0*datan(1d0)
!     call dinv33(qlat,0,qlatinv,det)
      call minv33(qlat,qlatinv)
      ncore1 = sum(ncore(iclass(1:nbas)))
      ncore2 = 0
      if (job == 0) then ! Polarizability context
C        off1 = nctot
C        nt = nt1
        expikt = 1
        cphiq2 => cphiqk
        if (npm == 2) ncore2 = ncore1
      else ! Self-energy context
C        off1 = 0
C        nt = nt1+nctot ! = include cores in intermediate states for the correlation case
        cphiq2 => cphiq
        do  ia = 1, nbas
          expikt(ia) = exp(img*tpi* sum(q2p*tr(:,ia)))
        end do
      endif

!$omp parallel do collapse(3)
        do it = it1lo, it1hi
            do itp = it2lo, it2hi
                do ipb = ipblo, ipbhi
                    zpsimb(ipb,itp,it) = 0
                end do
            end do
        end do
!$omp end parallel do

C --- Matrix elements with site-centered product functions ---
      if (ipblo <= npbmb) then ! Accumulate site product basis ME only if within spec'd range
        if (sum(ntpba(iclass(1:nbas))) /= npbmb) call rx('mkmelt: product basis mismatch')

C   ... Core <psi(k+q,t') | core(k,t) B(R,I)>
        call meltpbcpsi(npbmb,ncore1,ncore2,ipblo,ipbhi,it1lo,it1hi,it2lo,it2hi,nbas,nclass,iclass,naugm,
     .    nnlmx,nnlmv,nnlmc,mxnpba,ntpba,ipbfn,indexpbi,icore,ncore,nl,nnc,expikt,cphiq2,cphiq,ppb,
     o    zpsimb)
C       call zprm('zpsimb from meltpbcpsi',zpsimb,ipbhi-ipblo+1,ipbhi-ipblo+1,(it1hi-it1lo+1)*(it2hi-it2lo+1))
        if (lprt) call cputid2(' complete core matrix elements',0)

C   ... Valence <psi(k',t') | psi(k,t) exp(ikT) B(R,I)>
        call meltpbvpsi(npbmb,ncore1,ncore2,ipblo,ipbhi,it1lo,it1hi,it2lo,it2hi,nbas,nclass,iclass,naugm,
     .    nnlmx,nnlmv,nnlmc,mxnpba,ntpba,ipbfn,indexpbi,expikt,cphiqk,cphiq,ppb,
     o    zpsimb)
C       call zprm('zpsimb from meltpbvpsi',zpsimb,ipbhi-ipblo+1,ipbhi-ipblo+1,(it1hi-it1lo+1)*(it2hi-it2lo+1))
        if (lprt) call cputid2(' complete valence matrix elements',0)
      endif

C --- Contributions <psi |psi IPW> ---
      if (npwmb == 0) return  ! There are no IPW matrix elements
      if (ipbhi <= npbmb) return  ! No IPW matrix elements sought

C ... Read the G vectors and the PW coefficients for q1 and q1-q2
      delq = q1-q2
      call getngvec('QGpsi',q1,qlatinv,qt,ngq1)     ! q1 associated w/ 2nd band index in zpsimb
      call getngvec('QGpsi',delq,qlatinv,qmkt,ngq2) ! q2 associated w/ 1st band index in zpsimb
      allocate(igv1(3,ngq1),igv2(3,ngq2))
      call getgvec('QGpsi',q1,qlatinv,qt,ngq1,igv1) !qt,ngq1,igv1 are output: qt=q1 modulo G vector
      call getgvec('QGpsi',delq,qlatinv,qmkt,ngq2,igv2)
      allocate(geigq1(ngq1,nt2),geigq2(ngq2,nt1))
      call getpwz(q1,ngq1,it2,it2+nt2-1,itlst,isp1,qu1,geigq1)
      call getpwz(delq,ngq2,it1,it1+nt1-1,[0],isp2,qu2,geigq2)
C     Sanity checks
      if (sum(abs(qt-qu1))+sum(abs(qmkt-qu2)) > 1d-10) call rx('mkmelt: readqg and readgeig return different q')
C     Error if  Rk - q1 + (q-Rk) is not zero.
      qdiff = matmul(symope,q2p) - qt + qmkt
      add = matmul(qlatinv,qdiff); nadd = idint(add + dsign(.5d0,add))
      if (sum(abs(add-nadd)) > 1d-7) then
        print "(' mkmelt (warning): mismatch Rk - q1 + (q-Rk), diff=',1pg10.3)", sum(abs(add-nadd))
      endif
      if (sum(abs(add-nadd)) > 1d-6) call rx('mkmelt: mismatch Rk - q1 + (q-Rk) is too large')

C ... Make the interstitial matrix elements
C      call melpln4(ngq1,igv1,ngq2,igv2,npbmb,npwmb,nadd,geigq1,nt2,geigq2,nt1,nctot,ncore2,
C     .  ipblo,ipbhi,it1lo,it1hi,it2lo,it2hi,shftv,matmul(symope,q2p),q2p,symope,qlat,
C     o  zpsimb)
C

C ... Make the interstitial matrix elements
      call meltipw(npbmb,ncore1,ncore2,npwmb,ipblo,ipbhi,it1lo,it1hi,it2lo,it2hi,
     .  ngq1,igv1,ngq2,igv2,nadd,geigq1,nt2,geigq2,nt1,shftv,matmul(symope,q2p),q2p,symope,qlat,
     o  zpsimb)

      deallocate(igv1,igv2,geigq1,geigq2)

      if (lprt) call cputid2(' complete interstitial matrix elements',0)
C     call zprm('zpsimb from mkmelt',zpsimb,ipbhi-ipblo+1,ipbhi-ipblo+1,(it1hi-it1lo+1)*(it2hi-it2lo+1))

      end
