      subroutine meltpbvpsi(npbmb,ncore1,ncore2,ipblo,ipbhi,it1lo,it1hi,it2lo,it2hi,nbas,nclass,iclass,
     .  naugm,nnlmx,nnlmv,nnlmc,mxnpba,ntpba,ipbfn,indexpbi,phase,cphiqk,cphiq,ppb,
     o  zpsimb)
C- Matrix elements (site product basis fn * valence state * valence state)
C ----------------------------------------------------------------------
Ci Inputs
Ci  npbmb  : number of site-centered product basis functions B
Ci         : OLD nbloch
Ci  ncore1 : Number of core states in first band index of zpsimb array
Ci         : OLD nctot
Ci  ncore2 : Number of core states in second band index of zpsimb array
Ci         : OLD ncc
Ci  ipblo  : lower bound in zpsimb product basis index
Ci  ipbhi  : upper bound in zpsimb product basis index
Ci  it1lo  : lower bound in zpsimb 1st band index
Ci  it1hi  : upper bound in zpsimb 1st band index
Ci  it2lo  : lower bound in zpsimb 2nd band index
Ci  it2hi  : upper bound in zpsimb 2nd band index
Ci  nbas   : number of sites in the basis
Ci         : OLD natom
Ci  nclass : number of inequivalent classes
Ci  iclass : the jth atom belongs to class iclass(j)
Ci  naugm  : number of augmentation channels ... called ndima in sugw.f
Ci         : OLD nlmto
Ci  nnlmx  : maximum number of partial waves at a site; dimensions ppb
Ci         : OLD nlnmx
Ci  nnlmv  : nnlmv(ic) total number of partial waves (l,n,m) for valence states for site of class ic
Ci         : OLD nlnmv
Ci  nnlmc  : nnlmc(ic) total number of partial waves (l,n,m) for core states for site of class ic
Ci         : OLD nlnmc
Ci  mxnpba : dimensions ppb : max number of product basis functions for any site
Ci         : OLD mdimx
Ci  ntpba  : ntpba(ic) = number of product basis functions for class ic
Ci         : OLD mdim or nblocha
Ci  ipbfn  : ipbfn(i) = index to first product basis function for site i
Ci         : Product functions are strung together as a long vector.
Ci         : OLD: imdim
Ci indexpbi: indexpbi(i): for site i, points to the site which for product basis is stored
Ci         : OLD iatomp
Ci  phase  : exp(ik.T), where T it is translation Tiat, of site R to equivalent site under
Ci         : possible rotation of k to an equivalent kr (sxcf)
Ci         : In the polarizability context, phase should be unity.
Ci  cphiqk : (polarizability context; see Remarks)
Ci         : amount of partial wave contributing to occ valence eigenfunction at k
Ci         : (self-energy context; see Remarks)
Ci         : amount of partial wave contributing to intermediate state at q-k
ci         : See Memory management for indexing
Ci  cphiq  : (polarizability context)
Ci         : amount of partial wave contributing to unocc valence eigenfunction at k+q
Ci         : (self-energy context)
Ci         : amount of partial wave contributing to eigenstate n or m at q (see Remarks)
Ci  ppb    : (real) matrix elements of partial waves <Phi(RLn) Phi(RL'n') B(R,i)>
Ci         : ppb includes both core and valence.  Sequence :
Ci         : (core at 1), (val at 1) (core at 2), (val at 2) ...
Co Outputs
Co  zpsimb : matrix elements zpsimb(I,t,t') ... called X in Remarks
Co         : argument   polarizability-context     self-energy context
Co         :    I          product basis               product basis
Co         :    t         occ states at k          intermediate_states_at_q-k
Co         :    t'      unocc states at q+k     states for matrix elements at q
Cl Local
Cl  iplo   : iplo = ipblo = 1st product basis state to calculate ME; see Remarks
Cl  iphi   : last product basis index to calculate ME; see Remarks
Cl  i1lo   : first band (1st index) to calculate ME; see Remarks
Cl  i1hi   : i1hi = it1hi = last band (1st index) to calculate ME; see Remarks
Cl  i2lo   : first band (2nd index) to calculate ME; see Remarks
Cl  i2hi   : i2hi = it2hi = last band (2nd index) to calculate ME; see Remarks
Cl  iasx   : iasx(ia) = first element in cphi corresponding to site ia
Cl      ...  These variables are not required here, but they indicate how
Cl           zpsimb would be dimensioned if carried all matrix elements:
Cl           complex(8) :: zpsi2b(nmb,nt1,nt2),cphiqk(naugm,nt1),cphiq(naugm,nt2)
Cl  nmb    : total number of mixed product basis functions = npbmb + npwmb
Cl         : (formerly ngb)
Cl  nt1    : Polarizability context: no. occupied states at k, n in in Eq. 32, PRB76, 165106
Cl         : Self-energy context: no. intermediate states at q-k, n' in Eq. 34, PRB76, 165106
Cl         : OLD nt0
Cl  nt2    : Polarizability context: no. unoccupied states at q+k (polarizability context)
Cl         : Self-energy context: no. states at q for which sigm calculated (n or m in Eq 34)
Cl         : OLD ntp0
Cl
Cr Remarks
Cr
Cr  Calculates augmented part of matrix element for local product functions RI
Cr    X(k1,t1,k2,t2,RI) = <psi(k1,t1) | psi(k2,t2) B(RI)>
Cr    B(RI) = Bloch orthonormal product basis: I spans product functions at R
Cr    t1 and t2 each span a set of single-particle states.
Cr
Cr  This is used in two contexts:
Cr  1. In the polarizability, Eq. 32 in PRB76, 165106 :
Cr      X(k,t,q+k,t',RI) = <psi(q+k,t') B(RI) | psi(k,t) >
Cr     t runs over occupied states; t' runs over unoccupied states (called n and n' in Eq. 32)
Cr  2. In the self-energy, Eq. 34 in PRB76, 165106 :
Cr     X(q,t;q-k,t',RI) = <psi(q,t) | psi(q-k,t') B(RI)>
Cr     t  ranges over intermediate states, called n' in Eq. 34, in PRB76, 165106.
Cr     t' ranges over states for which sigma(t,t) is calculated; (n,m in Eq. 34)
Cr
Cr  Inside augmentation spheres eigenstate psi(t) has this partial wave expansion :
Cr    psi(k,t) = Sum[RL] phi(RL) cphi(RL,k,t)
Cr  phi(RL) is a partial wave at site R.
Cr  L is a compound index comprising l,m, type of radial wave function (phi, phidot, local orbital)
Cr
Cr  1. Polarizability context (see above for meaning of t and t')
Cr     cphi(RL',q+k,t') -> cphiq;   cphi(RL,k,t) -> cphiqk
Cr     Substitute one-center expansion:
Cr      X(k,t,q+k,t',RI) = <psi(q+k,t') B(RI) | psi(k,t)>
Cr                       =  Sum[L L'] cphiqk(L,t) * [<Phi(RL') Phi(RL) B(RI)> cphiq(L',t')]^+
Cr       for L'=1:nnlmvi, L=1:nnlmvi, t'=1:nt', t=1:nt.  nnlmvi = # partial waves at R (see nnlmv)
Cr  2. Self-energy context: (see above for meaning of t and t')
Cr     cphiqk = cphi(RL,q,t); cphiq = cphi(RL',q-k,t')
Cr     Substitute one-center expansion:
Cr       X(q,t;q-k,t',RI) = exp(ikT) Sum[L L'] cphiqk(L,t) * [<Phi(RL') Phi(RL) B(RI)> cphiq(L',t')]
Cr       for L'=1:nnlmvi, L=1:nnlmvi, t=1:nt1, t'=1:nt2  nnlmvi = # partial waves at R (see nnlmv)
Cr
Cr     X is coded as:
Cr       zz(L',t';I)    = Sum[L] <Phi(RL) Phi(RL') B(RI)> cphiq(L,t')
Cr       zpsimb(I,t,t') = exp(ikT) Sum[L'] [cphiqk(L',t)]^T zz^+(L',t')
Cr     zz is generated on the fly for one I at a time.  Matrix element <..> passed in ppb
Cr
Cm  Memory management
Cm
Cm   zpsimb(ip,i1,i2) contains matrix elements of product basis ip with eigenstates i1 and i2.
Cm     ip is a product basis function, spanning both site-centered functions B and plane waves
Cm     The npbmb site-centered product functions precede the npwmb PW states.
Cm     i1 can span both core states and valence states, the core states preceding the valence
Cm     i2 can span both core states and valence states (cores only for x0 without TR symmetry)
Cm   Ordering of states in zpsimb :
Cm     product basis index ip: site functions precede PWs
Cm     band indices i1,i2    : core states (if any) precede valence states
Cm
Cm   zpsimb is dimensioned as:
Cm     zpsimb(ipblo:ipbhi,it1lo:it1hi,it2lo:it2hi)
Cm
Cm   meltpbvpsi calculates matrix elements for the following subset of these states:
Cm     zpsimb(iplo:iphi,i1lo:i1hi,i2lo:i2hi)
Cm   This subset is restricted to ME of site PB functions with valence states, so
Cm     iplo = ipblo                  iphi = min(ipbhi,npbmb)
Cm     i1lo = max(it1lo,ncore1+1)    i1hi = it1hi
Cm     i2lo = max(it2lo,ncore2+1)    i2hi = it2hi
Cm   where ncore1,ncore2 are the number of core states that may have been allocated for
Cm   band indices 1,2, (nctot for index 1, and 0 or nctot for index 2; see off2 in mkmelt)
Cm   iplo, i1hi, and i2hi are redundant so they are not passed to meltpbvpsi.
Cm
Cm   cphiq and cphiqk contain coefficients for valence states only
Cm   This cphiq(*,1) refer to coefficients for the first valence state
Cm
Cu Updates
Cu     2016 (DMT) further redesign
Cu   Jun 16 (DMT) factored constituents of matrix elements to enable BLAS2 zgeru call
Cu   May 16 (ML)  Tweaks to work with MPI
Cu   Aug 14 (MvS) adapted from psi2b_v3
Cu   Feb 06 Takao wrote v3
Cu   Apr 02 Takao wrote v2
Cu   17 Mar 92 Ferdi's original code
C ----------------------------------------------------------------------
!     use mpi_mod, rank => mpi_rank
      implicit none
C ... Passed parameters
      integer, intent(IN) :: npbmb,ncore1,ncore2,ipblo,ipbhi,it1lo,it1hi,it2lo,it2hi
C     integer off1,off2,nt1,nt2
      integer, intent(IN) :: naugm,nnlmx,nbas,mxnpba,nclass
      integer, intent(IN) :: nnlmv(nclass),nnlmc(nclass),ntpba(nclass),iclass(nbas),ipbfn(nbas),indexpbi(nbas)
      real(8), intent(IN) :: ppb(nnlmx,nnlmx,mxnpba,nclass)
      complex(8), intent(IN) :: cphiqk(naugm,*),cphiq(naugm,*),phase(nbas)

C#ifdefC NT1NT2
C      complex(8), intent(INOUT) :: zpsimb(ipblo:ipbhi,it1lo:it1hi,it2lo:it2hi)
C#else
      complex(8), intent(INOUT) :: zpsimb(ipblo:ipbhi,it2lo:it2hi,it1lo:it1hi)
C#endif

C ... Local parameters
      integer ia,ic,nnlmci,nnlmvi,nc1,ias,iap,icp,i,ipb,iplo,i1hi,i2hi,n1,n2,i1lo,i2lo,iphi
      integer i2,jp
      complex(8) :: alpha,beta
      integer, allocatable:: iasx(:)
      complex(8), allocatable :: zz(:,:), zwork(:,:),zppb(:,:)

C ... Setup
!     if (nt1 == 0) return  ! Can happen for e.g. exchange calculation of core
!     If all (nt1,nt2) elements are made, generate this set:
!     zpsimb(ipb,ncore1+1:ncore1+nt1,ncore2+1:ncore2+nt2)
!     meltpbvpsi makes the following subset (see Memory Management above)
!     zpsimb(iplo:iphi,i1lo:i1hi,i2lo:i2hi)

      if (it1hi <= ncore1) return ! No valence states to calculate

      iplo = ipblo               ! Redundant, but use iplo locally for clarity
      iphi = min(ipbhi,npbmb)
      i1lo = max(it1lo,ncore1+1) ! Valence states only
      i1hi = it1hi               ! Redundant, but use i1hi locally for clarity
      i2lo = max(it2lo,ncore2+1) ! Valence states only
      i2hi = it2hi               ! Redundant, but use i2hi locally for clarity
      n1 = i1hi-i1lo+1; n2 = i2hi-i2lo+1
      if (iplo > iphi .or. n1 <= 0 .or. n2 <= 0) return ! Nothing to do

      allocate(zz(nnlmx,i2lo:i2hi)) ! If all matrix elements => zz(nnlmx,nt2)
C#ifdefC NT1NT2
C      allocate(zwork(i1lo:i1hi,i2lo:i2hi)) ! If all matrix elements => zwork(nt1,nt2)
C#else
      allocate(zwork(i2lo:i2hi,i1lo:i1hi)) ! If all matrix elements => zwork(nt1,nt2)
C#endif
      beta=0d0; alpha=1d0

C ... Make iasx = first element in expansion coefficients cphi corresponding to ia
      allocate(iasx(nbas))
      ias = 1
      do  ia = 1, nbas
        iasx(ia) = ias
        ias = ias + nnlmv(iclass(ia))
      enddo
      if (ias-1 /= naugm) call rx(' psi2b : sum(nnlmv)/= naugm')

C --- For each site R, make psi2b for each product basis fn i, occ state t, unocc t' ---
      do  ia = 1, nbas  ! loop over augmentation sites

        ic  = iclass(ia)
        nnlmci  = nnlmc(ic)
        nnlmvi  = nnlmv(ic)
        nc1 = nnlmci + 1

        if (nnlmci+nnlmv(ic) > nnlmx) call rxi('psi2b: improper dimensioning of partial waves, site',ia)

c       phase= dcmplx(coskt(ia),sinkt(ia))
        ias = iasx(ia)
        iap = indexpbi(ia)
        icp = iclass(iap)
!       if (ipbfn(iap)-1+ntpba(icp) > ipbhi) call rx('psi2b: wrong ipbhi') ! Test doesn't apply if zpsimb has been sliced

        allocate(zppb(nnlmvi,nnlmvi))
        do  i = 1, ntpba(icp)  ! Loop over product basis for this site (Note: loops are all independent)
          ipb = ipbfn(iap)-1+i ! Index to this product basis element in zpsimb

          if (ipb < iplo) cycle ! This product basis function not yet in window
          if (ipb > iphi) exit  ! This site has no more product basis functions in window

C         zz(L,t') = Sum[L'] [<Phi(RL') Phi(RL) B(RI)> cphiq(L',t')]
C         for L'=1:nnlmvi, L=1:nnlmvi, t'=1:nt2  nnlmvi = # partial waves at R = nnlmv(ic)
          zppb = ppb(nnlmci+1:nnlmci+nnlmvi,nnlmci+1:nnlmci+nnlmvi,i,icp) ! Copy portion of (real) ppb to complex temporary array
          call zgemm('T','N',nnlmvi,n2,nnlmvi,alpha,zppb,nnlmvi,cphiq(ias,i2lo-ncore2),naugm,beta,zz(1,i2lo),nnlmx)
C         if (ia == 2) then; print *, 'atom',ia; call zprm('psi2b zz',zz,nnlmx,nnlmvi,n2); endif

C#ifdefC NT1NT2
C          do  i2 = i2lo, i2hi  ! i2 = ncore2+1, ncore2+nt2 if all elements are to be calculated
C            do  jp = 1, nnlmvi
C              zz(jp,i2) = dconjg(zz(jp,i2))
C            enddo
C          enddo
C#endif

C         zpsimb(ipb,t,t') = zwork(t,t') = Sum[L] [exp(ikT) cphiqk(L,t)]^T zz^T(L,t')
C         for product basis fun ipb, L=1:nnlmvi, t=1:nt1 nnlmvi = # partial waves at R = nnlmv(ic)
C#ifdefC NT1NT2
C          call zgemm('T','N',n1,n2,nnlmvi,phase(ia),cphiqk(ias,i1lo-ncore1),naugm,zz(1,i2lo),nnlmx,beta,zwork(i1lo,i2lo),n1)
C#else
          call zgemm('C','N',n2,n1,nnlmvi,phase(ia),zz(1,i2lo),nnlmx,cphiqk(ias,i1lo-ncore1),naugm,beta,zwork(i2lo,i1lo),n2)
C#endif

C         if (ia == 2) then; print *, 'atom',ia; call zprm('psi2b zpsimb',zwork,n1,n1,n2); endif
C         zpsimb(ipb,off1+1:off1+nt1,off2+1:off2+nt2) = zwork
          zpsimb(ipb,i2lo:i2hi,i1lo:i1hi) = zwork(i2lo:i2hi,i1lo:i1hi)
!           call zcopy(n1*n2,zwork,1,zpsimb(ipb,i2lo,i1lo),ipbhi-ipblo+1)
C#endif

        end do                  !end of product basis-loop
        deallocate(zppb)
      end do                    !end of atom-loop ia

      deallocate(zz,zwork,iasx)
      end

      subroutine meltpbcpsi(npbmb,ncore1,ncore2,ipblo,ipbhi,it1lo,it1hi,it2lo,it2hi,nbas,nclass,iclass,
     .  naugm,nnlmx,nnlmv,nnlmc,mxnpba,ntpba,ipbfn,indexpbi,icore,ncore,nl,nnc,phase,cphiqk,cphiq,ppb,
     o  zpsimb)
C- Matrix elements (product basis * valence * core) states
C ----------------------------------------------------------------------
Ci Inputs
Ci  npbmb  : number of site-centered product basis functions B
Ci         : OLD nbloch
Ci  ncore1 : Number of core states in first band index of zpsimb array
Ci         : OLD nctot
Ci  ncore2 : Number of core states in second band index of zpsimb array
Ci         : OLD ncc
Ci  ipblo  : lower bound in zpsimb product basis index
Ci  ipbhi  : upper bound in zpsimb product basis index
Ci  it1lo  : lower bound in zpsimb 1st band index
Ci  it1hi  : upper bound in zpsimb 1st band index
Ci  it2lo  : lower bound in zpsimb 2nd band index
Ci  it2hi  : upper bound in zpsimb 2nd band index
Ci  nbas   : number of sites in the basis
Ci         : OLD natom
Ci  nclass : number of inequivalent classes
Ci  iclass : the jth atom belongs to class iclass(j)
Ci  naugm  : number of augmentation channels ... called ndima in sugw.f
Ci         : OLD nlmto
Ci  nnlmx  : maximum number of partial waves at a site; dimensions ppb
Ci         : OLD nlnmx
Ci  nnlmv  : nnlmv(ic) total number of partial waves (l,n,m) for valence states for site of class ic
Ci         : OLD nlnmv
Ci  nnlmc  : nnlmc(ic) total number of partial waves (l,n,m) for core states for site of class ic
Ci         : OLD nlnmc
Ci  mxnpba : dimensions ppb : max number of product basis functions for any site
Ci         : OLD mdimx
Ci  ntpba  : ntpba(ic) = number of product basis functions for class ic
Ci         : OLD mdim or nblocha
Ci  ipbfn  : ipbfn(i) = index to first product basis function for site i
Ci         : Product functions are strung together as a long vector.
Ci         : OLD: imdim
Ci indexpbi: indexpbi(i): for site i, points to the site which for product basis is stored
Ci         : OLD iatomp
Ci  icore  : icore(i,ic) index in ppb array for core state i in class ic
Ci  ncore  : ncore(ic) no. core states in class ic
Ci  nl     : (global maximum l) + 1
Ci  nnc    : dimensions icore
Ci  phase  : exp(ik.T), where T it is translation Tiat, of site R to equivalent site under
Ci         : possible rotation of k to an equivalent kr (sxcf)
Ci         : In the polarizability context, phase should be unity.
Ci  cphiqk : (polarizability context; see Remarks)
Ci         : amount of partial wave contributing to occ valence eigenfunction at k
Ci         : (self-energy context; see Remarks to meltpbvpsi)
Ci         : amount of partial wave contributing to intermediate state at q-k
Ci  cphiq  : (polarizability context)
Ci         : amount of partial wave contributing to unocc valence eigenfunction at k+q
Ci         : (self-energy context)
Ci         : amount of partial wave contributing to eigenstate n or m at q (see Remarks to meltpbvpsi)
Ci  ppb    : (real) matrix elements of partial waves <Phi(RLn) Phi(RL'n') B(R,i)>
Ci         : ppb includes both core and valence.  Sequence :
Ci         : (core at 1), (val at 1) (core at 2), (val at 2) ...
Cl Local
Cl  iplo   : iplo = ipblo = 1st product basis state to calculate ME; see Remarks
Cl  iphi   : last product basis index to calculate ME; see Remarks
Cl  i1lo   : first band (1st index) to calculate ME; see Remarks
Cl  i1hi   : i1hi = it1hi = last band (1st index) to calculate ME; see Remarks
Cl  i2lo   : first band (2nd index) to calculate ME; see Remarks
Cl  i2hi   : i2hi = it2hi = last band (2nd index) to calculate ME; see Remarks
Cl     ...   These variables are not used here, but they indicate how
Cl           zpsimb would be dimensioned if carried all matrix elements:
Cl           complex(8) :: zpsi2b(nmb,nt1,nt2),cphiqk(naugm,nt1),cphiq(naugm,nt2)
Cl  nmb    : total number of mixed product basis functions = npbmb + npwmb
Cl         : (formerly ngb)
Cl  nt1    : Polarizability context: no. occupied states at k, n in in Eq. 32, PRB76, 165106
Cl         : Self-energy context: no. intermediate states at q-k, n' in Eq. 34, PRB76, 165106
Cl         : OLD nt0
Cl  nt2    : Polarizability context: no. unoccupied states at q+k (polarizability context)
Cl         : Self-energy context: no. states at q for which sigm calculated (n or m in Eq 34)
Cl         : OLD ntp0
Co Outputs
Co  zpsimb : matrix elements zpsimb(I,t,t') (see Remarks for meltpbvpsi)
Co         : Fills in these elements
Co         : argument   polarizability-context     self-energy context
Co         :    I          product basis               product basis
Co         :    t          core states at k         intermediate_states_at_q-k
Co         :    t'       unocc states at q+k     states for matrix elements at q
Co         :  and also if ncore2 == ncore1
Co         :    t          occ states at k
Co         :    t'        core states at q+k
Cm Memory management
Cm  zpsimb(ip,i1,i2) contains matrix elements of product basis ip with eigenstates i1 and i2.
Cm    ip is a product basis function, spanning both site-centered functions B and plane waves
Cm    The npbmb site-centered product functions precede the npwmb PW states.
Cm    i1 can span both core states and valence states, the core states preceding the valence
Cm    i2 can span both core states and valence states (cores only for x0 without TR symmetry)
Cm  Ordering of states in zpsimb :
Cm    product basis index ip: site functions precede PWs
Cm    band indices i1,i2    : core states (if any) precede valence states
Cm
Cm  zpsimb is dimensioned as:
Cm    zpsimb(ipblo:ipbhi,it1lo:it1hi,it2lo:it2hi)
Cm  meltpbvpsi calculates matrix elements for the following subset of these states:
Cm    zpsimb(iplo:iphi,i1lo:i1hi,i2lo:i2hi)
Cm  Subset is restricted to ME of site PB functions with band1=core and band2=valence state, so
Cm    iplo = ipblo                  iphi = min(ipbhi,npbmb)
Cm    i1lo = max(it1lo,ncore1)      i1hi = min(it1hi,ncore1)
Cm    i2lo = max(it2lo,ncore2+1)    i2hi = it2hi
Cm  where ncore1,ncore2 are the number of core states that may have been allocated for
Cm  band indices 1,2, (nctot for index 1, and 0 or nctot for index 2; see ncore2 in mkmelt)
Cm
Cm  In the susceptibility context without TR symmetry, another subset of ME is calculated
Cm    zpsimb(iplo:iphi,j1lo:j1hi,j2lo:j2hi)
Cm  Subset is restricted to ME of site PB functions with band1=valence and band2=core state, so
Cm    iplo = ipblo                  iphi = min(ipbhi,npbmb)
Cm    j1lo = max(it1lo,ncore1+1)    j1hi = it1hi
Cm    j2lo = max(it2lo,ncore2)      j2hi = min(it2hi,ncore2)
Cr Remarks
Cr Calculates  <psi (k+q,t') |core(k,t) B(R,ipb)>, t'=unocc, t=occ
Cr   and also  <core(k+q,t') |psi (k,t) B(R,ipb)> if ncore2>0
Cr B(R,i) = Mixed basis
Cr core(k,t) = core states
Cr <psi(k+q,t') | core(k,t) B(R,i)> = Sum[RL]  cphiqk(RL,k+q,t')^* * ppb
Cu Updates
Cu   Aug 14 Mark adapted from psicb_v3
Cu   Feb 06 Takao wrote v3
Cu   Apr 02 Takao wrote v2
Cu   17 Mar 92 Ferdi's original code
C-------------------------------------------------------------------------
!     use mpi_mod, rank => mpi_rank
      implicit none
C ... Passed parameters
      integer,intent(in) ::  npbmb,mxnpba,nbas,ncore2,nclass,ncore1,nl,naugm,nnlmx,nnc,
     .  ipblo,ipbhi,it1lo,it1hi,it2lo,it2hi
      integer,intent(in) :: icore(nl*nl*nnc,nclass),ncore(nclass),nnlmv(nclass),
     .  nnlmc(nclass),ntpba(nclass),iclass(nbas),ipbfn(nbas),indexpbi(nbas)
      real(8),intent(in) :: ppb(nnlmx,nnlmx,mxnpba,nclass)
C     complex(8):: zpsimb(ngb,ncore1+nt1,ncore2+nt2),cphiq(naugm,*),cphiqk(naugm,*),phase(nbas)
C#ifdefC NT1NT2
C      complex(8):: zpsimb(ipblo:ipbhi,it1lo:it1hi,it2lo:it2hi),cphiq(naugm,*),cphiqk(naugm,*),phase(nbas)
C#else
      complex(8):: zpsimb(ipblo:ipbhi,it2lo:it2hi,it1lo:it1hi),cphiq(naugm,*),cphiqk(naugm,*),phase(nbas)
C#endif
C ... Local parameters
      integer i,ia,iap,ias,ipb,ic,icp,icr,ics,it,itp,nnlmci,nc1,nnlmvi,i1lo,i1hi,i2lo,i2hi
      integer iplo,j1lo,j1hi,j2lo,j2hi,iphi

      if (it1lo > ncore1) return ! No core states to calculate
!     if (ipblo > iphi)   return ! No product basis functions in window

      iplo = ipblo               ! Redundant, but use iplo locally for clarity
      iphi = min(ipbhi,npbmb)
      i1lo = it1lo               ! Redundant, but use i1lo locally for clarity
      i1hi = min(it1hi,ncore1)   ! Core states only
      i2lo = max(it2lo,ncore2+1) ! Valence states only
      i2hi = it2hi               ! Redundant, but use i2hi locally for clarity

      j1lo = max(it1lo,ncore1+1) ! Valence states only
      j1hi = it1hi
      j2lo = it2lo
      j2hi = min(it2hi,ncore2)   ! Core states only

C --- For each site R, make psi2b for each product basis fn i, core state t, unocc t' ---
      ipb = 0  ! Product basis index
      ias = 1  ! Running index in cphiq to 1st partial wave at current site
      ics = 0  ! Running index in zpsimb to 1st core level at current site

      do  ia = 1, nbas       !  loop over augmentation sites

        ic = iclass(ia)
        nnlmci = nnlmc(ic)
        nnlmvi = nnlmv(ic)
        nc1 = nnlmci + 1

C   --- Loop over product basis, core states t, valence states t' ---
C       zpsimb(ipb,t,t') = Sum[L] exp(ikT) cphiq(RL,t')^(*) * <Phi(RL) core(RL) B(R,i)>
C       Index  spans
C       ipb    product basis functions of valence states 1:ntpba(R)
C       L      partial valence partial waves 1:nnlmv(R), encompassing l,m, radial component
C       it     core states at R
C       itp    valence states

C       If ncore2 is nonzero, also make
C       zpsimb(ipb,t,t') = Sum[L] exp(-ikT) cphiqk(RL,t)^(*) * <Phi(RL) core(RL) B(R,i)>
C       it     valence states
C       itp    core states at R

        iap = indexpbi(ia)
        icp = iclass(iap)
        ipb  = ipbfn(iap)-1

        do  i = 1, ntpba(icp)
          ipb = ipb + 1
          if (ipb < iplo) cycle ! This product basis function not yet in window
          if (ipb > iphi) exit  ! This site has no more product basis functions in window

          do  it = ics+1, ics+ncore(ic)
            if (it < i1lo .or. it > i1hi) cycle ! Core state not in specified window
            icr = icore(it-ics,ic)
C#ifdefC NT1NT2
C            do  itp = i2lo, i2hi ! itp = ncore2+1, ncore2+nt2 if all elements are to be calculated
C              zpsimb(ipb,it,itp) = phase(ia)*dconjg(sum(cphiq(ias:ias+nnlmvi-1,itp-ncore2)*ppb(nc1:nnlmci+nnlmvi,icr,i,icp)))
C            enddo               ! end of t'(valence)-loop
C#else
            do  itp = i2lo, i2hi ! itp = ncore2+1, ncore2+nt2 if all elements are to be calculated
              zpsimb(ipb,itp,it) = phase(ia)*dconjg(sum(cphiq(ias:ias+nnlmvi-1,itp-ncore2)*ppb(nc1:nnlmci+nnlmvi,icr,i,icp)))
            enddo               ! end of t'(valence)-loop
C#endif
          enddo                 ! end of t(core)-loop

          if (ncore2 == 0) cycle  ! Branch below only for no TR symmetry

          do  itp = ics+1, ics+ncore(ic)
            if (itp < j2lo .or. itp > j2hi) cycle ! Core state not in specified window
            icr = icore(itp-ics,ic)
C#ifdefC NT1NT2
C            do  it = j1lo, j1hi ! it = ncore1+1, ncore1+nt1 if all elements are to be calculated
C              zpsimb(ipb,it,itp) = dconjg(phase(ia))*sum(cphiqk(ias:ias+nnlmvi-1,it-ncore1)*ppb(nc1:nnlmci+nnlmvi,icr,i,icp))
C            enddo               ! end of t(valence)-loop
C#else
            do  it = j1lo, j1hi ! it = ncore1+1, ncore1+nt1 if all elements are to be calculated
              zpsimb(ipb,itp,it) = dconjg(phase(ia))*sum(cphiqk(ias:ias+nnlmvi-1,it-ncore1)*ppb(nc1:nnlmci+nnlmvi,icr,i,icp))
            enddo               ! end of t(valence)-loop
C#endif
          enddo                 ! end of t'(core)-loop

        enddo                   ! product basis-loop

        ias = ias + nnlmv(ic)
        ics = ics + ncore(ic)

      enddo                     ! loop over sites

      end

      subroutine meltipw(npbmb,ncore1,ncore2,npwmb,ipblo,ipbhi,it1lo,it1hi,it2lo,it2hi,
     .  ngq1,igv1,ngq2,igv2,nadd,geigq1,nt2,geigq2,nt1,shftv,q1,q2,symope,qlat,zpsimb)
C- Interstitial plane wave contributions to matrix elements <M psi |psi>
C ----------------------------------------------------------------------
Ci Inputs
Ci  npbmb  : number of site-centered product basis functions B
Ci         : OLD nbloch
Ci  ncore1 : Number of core states in first band index of zpsimb array
Ci         : OLD nctot
Ci  ncore2 : Number of core states in second band index of zpsimb array
Ci         : OLD ncc
Ci  npwmb  : number of PWs in interstitial product basis
Ci         : OLD: ngc
Ci  ipblo  : lower bound in zpsimb product basis index
Ci  ipbhi  : upper bound in zpsimb product basis index
Ci  it1lo  : lower bound in zpsimb 1st band index
Ci  it1hi  : upper bound in zpsimb 1st band index
Ci  it2lo  : lower bound in zpsimb 2nd band index
Ci  it2hi  : upper bound in zpsimb 2nd band index
Ci  ngq1   : number of G-vectors representing eigenfunction at q1
Ci         : OLD: ngp1
Ci  igv1   : The G-vectors at q1 (as integer multiples of qlat)
Ci         : OLD: ngvecp1
Ci  ngq2   : number of G-vectors representing eigenfunction at q2
Ci         : OLD: ngp2
Ci  igv2   : The ngq2 G-vectors at q2 (integer multiples of qlat)
Ci         : OLD: ngvecp2
Ci  nadd   : G-vector shift, normally 0 (integer multiples of qlat)
Ci  geigq1 : Coefficients to PW representation eigenvectors associated with q1
Ci  nt2    : Polarizability context: no. unoccupied states at q+k
Ci         : Self-energy context:    no. states at q for which sigm calculated (n or m in Eq 34)
Ci         : Note nt2 is NOT associated with q2, but with 2nd band index in zpsimb
Ci         : OLD: ntp0
Ci  geigq2 : Coefficients to PW representation eigenvectors associated with q2
Ci  nt1    : Polarizability context: no. occupied states at k, n in Eq. 32, PRB76, 165106
Ci         : Self-energy context:    no. intermediate states at q-k, n' in Eq. 34, PRB76, 165106
Ci         : Note nt1 is NOT associated with q1, but with 1st band index in zpsimb
Ci         : OLD: nt0
Ci  shftv  : shift entering into phase for q1+G
Ci  q1     : Polarizability context: q+k, calculating <psi(q+k,t') PW | psi(k,t) >
Ci         : used in integration over k to make Pi(q). Eq. 32, PRB76, 165106
Ci         : Self-energy context: q at which sigma is calculated. Eq. 34, PRB76, 165106
Ci  q2     : Polarizability context: k, for ME <psi(q+k,t') PW | psi(k,t) > .. Eq. 32, PRB76, 165106
Ci         : Self-energy context: k for integration dk psi*(q) psi(q-k)    ... Eq. 34, PRB76, 165106
Ci         : q2 = q - R*k (R is rotation)
Ci  symope : Rotation R
Ci  qlat   : reciprocal lattice vectors
Ci         : NEW: qlat
Co Outputs
Co   zpsimb :The IPW part of matrix elements with the product basis, <psi |psi Plane>
Cr Remarks
Cr   use PPOVL through getppx in m_read_ppovl
Cm Memory management
Cm  zpsimb(ip,i1,i2) contains matrix elements of product basis ip with eigenstates i1 and i2.
Cm    ip is a product basis function, spanning both site-centered functions B and plane waves
Cm     The npbmb product functions precede the npwmb PW states
Cm    i1 can span both core states and valence states, the core states preceding the valence
Cm    i2 can span both core states and valence states (cores only for x0 without TR symmetry)
Cm  Ordering of states in zpsimb :
Cm    product basis index ip: site functions precede PWs
Cm    band indices i1,i2:     core states (if any) precede valence states
Cm
Cm  This routine calculates matrix elements for the following subset of these states:
Cm    zpsimb(ipblo:ibphi,i1lo:i1hi,i2lo:i2hi)
Cm  zpsimb is dimensioned as:
Cm    zpsimb(ipblo:ipbhi,it1lo:it1hi,it2lo:it2hi)
Cm
Cm  ipblo,i1lo may differ from ipblo,it1lo => the ME generated are not contiguous in memory.
Cm     fixed i1,i2            fixed ip,i2
Cm  ----------------------------------
Cm      ipblo                   it1lo
Cm      ...                     ...
Cm      npbmb                   ncore1
Cm      ipblo                   i1lo
Cm      ...                     ...
Cm      ibphi                   it1hi
Cm  ----------------------------------
Cm
Cm  zpsimb is scaled by the IPW overlap <Gc|Gc>^-1.
Cm  For efficiency, scaling is performed by a BLAS zgemm call which means
Cm  the latter two dimensions must be bundled as a single vector.
Cm
Cm  The obvious workaround -- generate ME in a temporary array and copy them
Cm  to the appropriate location --- can require large amounts of memory
Cm  (zpsimb can be very large).  Instead, this routine uses zpsimb as a work array;
Cm  then the elements for a fixed ip are copied into a temporary array,
Cm  and restored to zpsimb in the proper order.  Since all the elements of one IPW
Cm  are distinct from another, this can be accomplished without overwriting elements.
Cu Updates
Cu   07 May 13
C ----------------------------------------------------------------------
      use m_read_ppovl,only: getppx, ppx, ngc=>ngc2, igvc=>ngvecc2
      implicit none
C ... Passed parameters
      integer, intent(IN) :: npbmb,ncore1,ncore2,ipblo,ipbhi,it1lo,it1hi,it2lo,it2hi
      integer, intent(IN) :: ngq1,ngq2,npwmb,nt2,nt1
      integer, intent(IN) :: igv1(3,ngq1),igv2(3,ngq2),nadd(3)
      real(8):: shftv(3),q1(3),q2(3),qlat(3,3),symope(3,3)
C     complex(8):: geigq1(ngq1,nt2),geigq2(ngq2,nt1),zpsimb(nbloch+npwmb,off1+nt1,off2+nt2)
C#ifdefC NT1NT2
C      complex(8):: geigq1(ngq1,nt2),geigq2(ngq2,nt1),zpsimb(ipblo:ipbhi,it1lo:it1hi,it2lo:it2hi)
C#else
      complex(8):: geigq1(ngq1,nt2),geigq2(ngq2,nt1),zpsimb(ipblo:ipbhi,it2lo:it2hi,it1lo:it1hi)
C#endif
C ... Local parameters
      integer,allocatable:: firstb(:),lastb(:),ipp(:,:,:)
      complex(8),allocatable :: gg(:,:,:),zwk(:,:),zwk2(:)
      integer ipb,itp,iplo,iphi,i1lo,i1hi,i2lo,i2hi,n1,n2,npw,iverbose,verbose
      integer ngcx,nxmin,nxmax,nymin,nymax,nzmin,nzmax,ig2,ig1,igc2
      integer igv12(1:3)
      integer omppid,idalloc,nblk,nthreads,tid,iblk1,iblk2,igmx,igmn,nblki,igc0,it
      real(8):: qg(3),pi=3.1415926535897932D0
      complex(8):: img=(0d0,1d0),phas,zfac

      integer, allocatable :: igpairs(:,:,:), igcnts(:)
      real(8) :: lshftv(3)
      complex(8), allocatable :: g1tmp(:,:), g2tmp(:,:)
!       integer, parameter  :: gcmax = 512*5 !100 !512 !512*5
      integer :: gcmax, gcmin
      real(8) :: gcav
      complex(8), parameter :: cmplx0 = 0.0_8, cmplx1 = 1.0_8, cmplxi = (0.0_8,1.0_8)
      integer(8) :: c1, c2, cr, cm, cc(4)
      call system_clock(count_rate = cr, count_max = cm)

!      write(*,'("[",I4,"] meltipw: start. ",9I5)') rank, ipblo,ipbhi,it1lo,it1hi,it2lo,it2hi,npbmb,ncore1,ncore2
!      write(*,'("[",I4,"] meltipw: nqg1/2 = ",2I4,", nt1/2 = ",2I4)') rank, ngq1, ngq2, nt1, nt2

!      write(*,'("[",I4,"] meltipw: abs(sum(geigq1)) = ",E15.5,", abs(sum(geigq2)) = ",E15.5)')
!     &     rank, abs(sum(geigq1)), abs(sum(geigq2))

C ... Range of states to calculate
      iplo = max(ipblo,npbmb+1)  ! IPW's only
      iphi = ipbhi               ! Redundant, but include for clarity
      i1lo = max(it1lo,ncore1+1) ! No matrix elements of IPWs with core states
      i1hi = it1hi               ! Redundant, but include for clarity
      i2lo = max(it2lo,ncore2+1) ! No matrix elements of IPWs with core states
      i2hi = it2hi               ! Redundant, but include for clarity
      npw = iphi-iplo+1; n1 = i1hi-i1lo+1; n2 = i2hi-i2lo+1
      if (iplo > iphi .or. n1 <= 0 .or. n2 <= 0) return ! Nothing to do

C ... Read from getppx : ngc, ppx(1:npwmb,1:ngc), igvc(1:3,1:ngc)
C     ppx = <Gc|Gc>^-1 <Gc|Gc2> (Gc2 spans G1-G2) made in in rdata4gw_v2 (Sep12 version)
      iverbose = (verbose()/50)*10
      call getppx(qlat,ngcx,q2,symope)
      if (npwmb/=ngcx) call rx('meltipw: file mismatch with ppovl, G vector count')
      nxmin = minval(igvc(1,1:ngc)); nxmax = maxval(igvc(1,1:ngc))
      nymin = minval(igvc(2,1:ngc)); nymax = maxval(igvc(2,1:ngc))
      nzmin = minval(igvc(3,1:ngc)); nzmax = maxval(igvc(3,1:ngc))

C     Range check
      do  ig1 = 1, ngq1
      do  ig2 = 1, ngq2
        igv12(1:3) = igv1(1:3,ig1)-igv2(1:3,ig2) - nadd
        if (igv12(1)<nxmin .or. igv12(1)>nxmax .or.
     .      igv12(2)<nymin .or. igv12(2)>nymax .or.
     .      igv12(3)<nzmin .or. igv12(3)>nzmax) then
          ! call rx('meltipw: insufficient range for G vectors')
           print *, 'meltipw (warning!!) : insufficient range for G vectors'
           nxmin = min(igv12(1),nxmin); nxmax = max(igv12(1),nxmax)
           nymin = min(igv12(2),nymin); nymax = max(igv12(2),nxmax)
           nzmin = min(igv12(3),nzmin); nzmax = max(igv12(3),nxmax)
         endif
      enddo
      enddo

C ... Invert igvc:  given G=G1,G2,G3, ipp(G1,G2,G3) identifies igc2 corresponding to igvc(:,igc2)
      allocate(ipp(nxmin:nxmax,nymin:nymax,nzmin:nzmax))
      ipp = -10000000
      do  igc2 = 1, ngc
        ipp(igvc(1,igc2),igvc(2,igc2),igvc(3,igc2)) = igc2
      enddo


C --- Assemble gg in slices to conserve memory.  Slice ranges from igmn to igmx ---
C     nblk is max number of G2 vectors.  Should be adjusted depending on total size of gg
!       nblk = ngc ! set to smaller size for big systems.
!       nblk = ngc means there will be only 1 iteration of 'do  igmn = 1, ngc, nblk' below
!       nblk = 2000

!       print "('meltipw: n1,n2,ngc,npw',4(x,i0))", n1,n2,ngc,npw
!       allocate(gg(n2,n1,nblk))

      allocate(igcnts(ngc))
!       allocate(igpairs(2,gcmax,ngc),igcnts(ngc))
! !       allocate(g1tmp(gcmax,n2),g2tmp(gcmax,n1))
!       allocate(g1tmp(n2,gcmax),g2tmp(n1,gcmax))
C     print *, 'starting time',nblk,ngc,ngq2,delwc()


      call system_clock(c1)
      lshftv = matmul(transpose(qlat),shftv)
      do ig1 = 1, ngq1
          phas = exp(-2*pi*img*sum(igv1(1:3,ig1)*lshftv))
          call zscal(n2, phas, geigq1(ig1,i2lo-ncore2), ngq1)
      end do

      do ig2 = 1, ngq2
          phas = exp(-2*pi*img*sum(igv2(1:3,ig2)*lshftv))
          call zscal(n1, phas, geigq2(ig2,i1lo-ncore1), ngq2)
      end do

!       phas = exp(2*pi*img*sum(q1*shftv - nadd*lshftv))
      phas = exp(-2*pi*img*sum(q1*shftv - nadd*lshftv))

!         +p -2 +1

!         -p +1 -2
!         +p -1 +2

      call system_clock(c2)
!       write(*,'("meltipw init stuff:",x,g20.12)') real(c2 - c1,8)/real(cr,8)

      call system_clock(c1)

! ! Only used with gerc
! !       gg = 0
! !       call zcopy(n1*n2*ngc,cmplx0,0,gg,1)
! !$omp parallel do
!       do igc0 = 1, ngc
!         gg(:,:,igc0) = 0
!       end do
! !$omp end parallel do


! count only
      igcnts = 0
      do  ig1 = 1, ngq1
        do  ig2 = 1, ngq2
          igv12(1:3) = igv1(1:3,ig1)-igv2(1:3,ig2) - nadd  ! G1-G2
          igc0       = ipp(igv12(1),igv12(2),igv12(3))     ! Index to igvc
          if (igc0 > 0) igcnts(igc0) = igcnts(igc0) + 1
        end do
      end do


      gcmax = 0
      gcmin = huge(gcmin)
      gcav = 0
      igc0 = 0
      do ig1 = 1, ngc
        ig2 = igcnts(ig1)
        if (ig2 > 0) then
            if (ig2 > gcmax) gcmax = ig2
            if (ig2 < gcmin) gcmin = ig2
            gcav = gcav + real(ig2,8)
            igc0 = igc0 + 1
        end if
      end do
      gcav = gcav/real(igc0,8)

!      write(*, "('meltipw packing min, max & av:',2(x,i0),x,f6.1)") gcmin, gcmax, gcav
!       gcmax = 8*((gcmax+7)/8) ! align to 8 doubles (512 bits)
      allocate(igpairs(2,gcmax,ngc))
      allocate(g1tmp(n2,gcmax),g2tmp(n1,gcmax))

      igcnts = 0
      do  ig1 = 1, ngq1
        do  ig2 = 1, ngq2
          igv12(1:3) = igv1(1:3,ig1)-igv2(1:3,ig2) - nadd  ! G1-G2
          igc0       = ipp(igv12(1),igv12(2),igv12(3))     ! Index to igvc
          if (igc0 > 0) then
            igcnts(igc0) = igcnts(igc0) + 1
            igpairs(1:2,igcnts(igc0),igc0) = [ig1, ig2]
!           call zgerc(n2, n1, phas, geigq1(ig1,i2lo-ncore2), ngq1, geigq2(ig2,i1lo-ncore1), ngq2, gg(i2lo,i1lo,igc0), n2)
          else
! Quick and dirty patch for fuzzy vectors. Apparently they need to be handled properly in rdata4gw
            write(0, '("Warning: discarding missing connecting G vector (",3(x,i0),")")') igv12(1:3)
          end if
        enddo
      enddo


      call system_clock(c2)
!       write(*,'("meltipw gpack s:",x,g20.12)') real(c2 - c1,8)/real(cr,8)

      nblk =  min(((ngc+7)/8)*8, max(1,((4*2**23)/(n1*n2))*8)) ! made for ~4GiB and multiple of 8
!       nblk = 1
!       write(*,'("meltipw nblk, ggsz:",2(x,i0)," MiB")') nblk, (n1*n2*nblk)/2**16

      allocate(gg(n2,n1,nblk))
      zfac = cmplx0
      igmx = 1
      call system_clock(c1)
      do igc0 = 1, ngc
!             print *, 'igc0, igcnts', igc0, igcnts(igc0), igmx
          if (igcnts(igc0) /= 0) then
            do igc2 = 1, igcnts(igc0)
                ig1 = igpairs(1,igc2,igc0)
                ig2 = igpairs(2,igc2,igc0)

!                 call zcopy(n2, geigq1(ig1,i2lo-ncore2), ngq1, g1tmp(igc2,1), gcmax)
!                 call zcopy(n1, geigq2(ig2,i1lo-ncore1), ngq2, g2tmp(igc2,1), gcmax)
                call zcopy(n2, geigq1(ig1,i2lo-ncore2), ngq1, g1tmp(1,igc2), 1)
                call zcopy(n1, geigq2(ig2,i1lo-ncore1), ngq2, g2tmp(1,igc2), 1)
            end do

!           call zgemm('c','n',n2,n1,igcnts(igc0),phas,g1tmp,gcmax,g2tmp,gcmax,cmplx1,gg(1,1,igc0),n2) !with ph=exp(+) and zgemm(n,c,ppx) below
!           call zgemm('n','c',n2,n1,igcnts(igc0),phas,g1tmp,n2,g2tmp,n1,cmplx0,gg(1,1,igc0),n2)
!           igmx = mod(igc0-1,nblk)+1

!            print *,'n2,n1,igcnts(igc0),igmx',n2,n1,igcnts(igc0),igmx
            call zgemm('n','c',n2,n1,igcnts(igc0),phas,g1tmp,n2,g2tmp,n1,cmplx0,gg(1,1,igmx),n2)
          else
            gg(:,:,igmx) = 0
          end if

          if (igmx == nblk .or. igc0 == ngc) then
            call zgemm('n','t',npw,n1*n2,igmx,cmplx1,
     &        ppx(1,igc0-igmx+1),npwmb,gg,n2*n1,zfac,zpsimb(iplo,i2lo,i1lo),ipbhi-ipblo+1)
            zfac = cmplx1
            igmx = 0
          end if

          igmx = igmx + 1
      end do

      call system_clock(c2)
!       write(*,'("meltipw zgemm s:",x,g20.12)') real(c2 - c1,8)/real(cr,8)

      call system_clock(c1)

C     call zprm('zpsimb',zpsimb(iplo,i1lo,i2lo),ipbhi-ipblo+1,npw,n1*n2)
C     call zprm('zpsimb',zpsimb(iplo,i1lo,i2lo),ipbhi-ipblo+1,npw,n1*n2)
!       igc2 = idalloc('gg',iverbose+4,nt1*nt2,nblk*2)
      deallocate(ppx,gg,igvc,ipp)
      deallocate(g1tmp,g2tmp,igcnts,igpairs)

      call system_clock(c2)
!       write(*,'("meltipw dealloc:",x,g20.12)') real(c2 - c1,8)/real(cr,8)

      call system_clock(c1)
C ... Take CC to turn <P^Gc Psi2| Psi1> -> <Psi1| Psi2 P^Gc >
C     and undo misalignment in zpsimb, columns 2 and 3

      if (it2lo == i2lo) then ! no misalignment in dimensions 2 and 3 in zpsimb
!$omp parallel do collapse(3)
        do it = it1lo, i1lo-1
            do itp = it2lo, it2hi
                do ipb = iplo, iphi
                    zpsimb(ipb,itp,it) = 0 ! Assumed to be core states
                end do
            end do
        end do
!$omp end parallel do
!$omp parallel do collapse(3)
        do it = i1lo, i1hi
            do itp = i2lo, i2hi
                do ipb = iplo, iphi
                    zpsimb(ipb,itp,it) = conjg(zpsimb(ipb,itp,it))
                end do
            end do
        end do
!$omp end parallel do
      else ! likely notimereversal case, undo the misalignment caused by zgemm(...,npw,n1*n2,igmx,...) above
        allocate(zwk(i2lo:i2hi,i1lo:i1hi))
        do ipb = iplo, iphi
            call zcopy(n2*n1,zpsimb(ipb,i2lo,i1lo),ipbhi-ipblo+1,zwk,1)
            if (it1lo < i1lo) zpsimb(ipb,:,it1lo:i1lo-1) = 0 ! Assumed to be core states
            if (it2lo < i2lo) zpsimb(ipb,it2lo:i2lo-1,:) = 0 ! Assumed to be core states
            forall (it = i1lo:i1hi,itp=i2lo:i2hi) zpsimb(ipb,itp,it) = conjg(zwk(itp,it))
        end do
        deallocate(zwk)
      end if

      call system_clock(c2)
!       write(*,'("meltipw fin:",x,g20.12)') real(c2 - c1,8)/real(cr,8)

      end
