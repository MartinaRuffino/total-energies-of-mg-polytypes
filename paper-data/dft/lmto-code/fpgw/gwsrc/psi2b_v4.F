C Other changes: remove nbloch here and in psicb.  Rename pbmevpsi
      subroutine psi2b_v4(nctot,ncc,n1,n2,ndpsi1,ndpsi2,ndpsi3,
     .                   iclass,phase,
     .                   cphiqk,
     .                   cphiq,
     .                   ppb,
     .                   nlnmv,nlnmc,mdim,imdim,iatomp,
     .                   mdimx,nlmto,nbloch,nlnmx,natom,nclass,
     o                   zpsipb)
C- Matrix elements (product basis * valence states * valence states)
C ----------------------------------------------------------------------
Ci Inputs
Ci  nctot  : Offset in 2nd index zpsipb array: first element is zpsipb(:,off2+1,:)
Ci         : Typically off2 = number of core states.
Ci         : NEW off2
Ci  ncc    : Offset in 3rd index zpsipb array: first element is zpsipb(:,:,off3+1)
Ci         : Typically off3 = number of core states when there is no time reversal symmetry.
Ci         : NEW off3
Ci  n1     : Polarizability context: no. occupied states at k, n in in Eq. 32, PRB76, 165106
Ci         : Self-energy context: no. intermediate states at q-k, n' in Eq. 34, PRB76, 165106
Ci         : OLD nt0
Ci  n2     : Polarizability context: no. unoccupied states at q+k (polarizability context)
Ci         : Self-energy context: no. states at q for which sigm calculated (n or m in Eq 34)
Ci         : OLD ntp0
Ci  ndpsi1,ndpsi2,ndpsi3 : parameters that dimension zpsipb
Ci  iclass : the jth atom belongs to class iclass(j)
Ci  phase  : exp(ik.T), where T it is translation Tiat, of site R to equivalent site under
Ci         : possible rotation of k to an equivalent kr (sxcf_fal3_scz)
Ci         : In the polarizability context, phase should be unity.
Ci  cphiqk : (polarizability context; see Remarks)
Ci         : amount of partial wave contributing to occ valence eigenfunction at k
Ci         : (self-energy context; see Remarks)
Ci         : amount of partial wave contributing to intermediate state at q-k
Ci  cphiq  : (polarizability context)
Ci         : amount of partial wave contributing to unocc valence eigenfunction at k+q
Ci         : (self-energy context)
Ci         : amount of partial wave contributing to eigenstate n or m at q (see Remarks)
Ci  ppb    : (real) matrix elements of partial waves <Phi(RLn) Phi(RL'nn') B(R,i)>
Ci         : ppb includes both core and valence.  Sequence :
Ci         : (core at 1), (val at 1) (core at 2), (val at 2) ...
Ci  nlnmv  : For each class, number of partial waves for valence
Ci  nlnmc  : For each class, number of partial waves for core
Ci  mdim   : Number of product basis functions for each class
Ci         : NEW nrprdc
Ci  imdim  : imdim(i) = index to first product basis function for site i
Ci         : Product functions are strung together as a long vector.
Ci         : NEW: irprd
Ci  iatomp : mapping to site containing product basis?
Ci  mdimx  : dimensions ppb
Ci         : NEW: nrprdcx
Ci  nlmto  : number of augmentation channels ... called ndima in sugw.f
Ci         : NEW: naugm
Ci  nbloch : number of product basis functions B.  Used here only as a sanity check
Ci           NEW nrprd
Ci  nlnmx  : maximum number of product basis functions; dimensions ppb
Ci  natom  : number of sites
Ci           NEW nbas
Ci  nclass : number of inequivalent classes
Co Outputs
Co  zpsipb : matrix elements zpsipb(I,t,t') ... called X in Remarks
Co         : argument   polarizability-context     self-energy context
Co         :    I          product basis               product basis
Co         :    t         occ states at k          intermediate_states_at_q-k
Co         :    t'      unocc states at q+k     states for matrix elements at q
Cl Local
Cl   iasx  : iasx(ia) = first element in cphi corresponding to site ia
Cr Remarks
Cr  Calculates augmented part of matrix element for local product functions RI
Cr    X(k1,t1,k2,t2,RI) = <psi(k1,t1) | psi(k2,t2) B(RI)>
Cr    B(RI) = Bloch orthonormal product basis: I spans product functions at R
Cr    t1 and t2 each span a set of single-particle states.
Cr  This is used in two contexts:
Cr  1. In the polarizability, Eq. 32 in PRB76, 165106 :
Cr      X(k,t,q+k,t',RI) = <psi(q+k,t') B(RI) | psi(k,t) >
Cr     t runs over occupied states; t' runs over unoccupied states (called n and n' in Eq. 32)
Cr  2. In the self-energy, Eq. 34 in PRB76, 165106 :
Cr     X(q,t;q-k,t',RI) = <psi(q,t) | psi(q-k,t') B(RI)>
Cr     t  ranges over intermediate states, called n' in Eq. 34, in PRB76, 165106.
Cr     t' ranges over states for which sigma(t,t) is calculated; (n,m in Eq. 34)
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
Cr       for L'=1:nv, L=1:nv, t'=1:nt', t=1:nt.  nv = # partial waves at R (see nlnmv)
Cr  2. Self-energy context: (see above for meaning of t and t')
Cr     cphiqk = cphi(RL,q,t); cphiq = cphi(RL',q-k,t')
Cr     Substitute one-center expansion:
Cr       X(q,t;q-k,t',RI) = exp(ikT) Sum[L L'] cphiqk(L,t) * [<Phi(RL') Phi(RL) B(RI)> cphiq(L',t')]
Cr       for L'=1:nv, L=1:nv, t=1:n1, t'=1:n2  nv = # partial waves at R (see nlnmv)
Cr
Cr     X is coded as:
Cr       zz(L',t';I)    = Sum[L] <Phi(RL) Phi(RL') B(RI)> cphiq(L,t')
Cr       zpsipb(I,t,t') = exp(ikT) Sum[L'] [cphiqk(L',t)]^T zz^+(L',t')
Cr     zz is generated on the fly for one I at a time.  Matrix element <..> passed in ppb
Cu Updates
Cu      Aug 14 Mark  adapted from psi2b_v4
Cu      Feb 06 Takao wrote v3
Cu      Apr 02 Takao wrote v2
Cu   17 Mar 92 Ferdi's original code
C ----------------------------------------------------------------------
!     use mpi_mod, rank => mpi_rank
      implicit none
      integer, intent(IN) :: nlmto,nctot,ncc,n1,n2,nlnmx,natom,nbloch,mdimx,nclass,ndpsi1,ndpsi2,ndpsi3
      integer, intent(IN) :: nlnmv(nclass),nlnmc(nclass),mdim(nclass),iclass(natom),imdim(natom),iatomp(natom)
! in hx0fp0_sc:  real(8):: ppb(nlnmx*nlnmx*mdimx*nclass,nspin) but passed as ppb(1,is)
! in x0kf_mvs:   real(8)   ppb(*)
! in ppbasp:     real(8):: ppb(mnlx,mnlx,mdimx,nclass) ! where mnlx = nlnmx
      real(8), intent(IN) ::  ppb(nlnmx,nlnmx,mdimx,nclass)
      complex(8), intent(IN) :: cphiqk(nlmto,n1),cphiq(nlmto,n2),phase(natom)
      complex(8) :: zpsipb(ndpsi1,ndpsi2,ndpsi3)
C ... Local variables
      integer ia,ic,nc,nv,nc1,ias,iap,icp,i,itp,jp,ib
      complex(8) :: alpha,beta
      integer, allocatable:: iasx(:)
      complex(8), allocatable :: zz(:,:), zwork(:,:),zppb(:,:)

      if (n1 == 0) return  ! Can happen for e.g. exchange calculation of core

      allocate(zz(nlnmx,n2),zwork(n1,n2)) !Work areas
      beta=0d0  ; alpha=1d0

C ... Check dimensions
      if (mdimx /= maxval(mdim) ) call rx('psi2b: wrong mdimx')
      if (sum(mdim(iclass(1:natom))) /= nbloch) call rx('psi2b: wrong nbloch')
      if (nctot+n1 > ndpsi2) call rx('psi2b: wrong ndpsi2')

C ... Make iasx = first element in expansion coefficients cphi corresponding to ia
      allocate(iasx(natom))
      ias = 1
      do  ia = 1, natom
        iasx(ia) = ias
        ias = ias + nlnmv(iclass(ia))
      enddo
      if (ias-1 /= nlmto) call rx(' psi2b : sum(nlnmv)/= nlmto')

C --- Loop over atoms ---
C     For each site R, make psi2b for each product basis fn i, occ state t, unocc t'
      do  ia = 1, natom
        ic  = iclass(ia)
        nc  = nlnmc(ic)
        nv  = nlnmv(ic)
        nc1 = nc + 1
        if (nc+nlnmv(ic) > nlnmx) call rx('psi2b: nlnmx exceeded')
c       phase= dcmplx(coskt(ia),sinkt(ia))
        ias = iasx(ia)
        iap = iatomp(ia)
        icp = iclass(iap)
        if (imdim(iap)-1+mdim(icp) > ndpsi1) call rx('psi2b: wrong ndpsi1')

C       MPI note: each i in this loop is independent ... can be parallelized over i
        allocate(zppb(nv,nv))
        do  i = 1, mdim(icp) ! loop over product basis for this site
          ib = imdim(iap)-1+i ! Index to this product basis element in zpsipb

C         zz(L,t') = Sum[L'] [<Phi(RL') Phi(RL) B(RI)> cphiq(L',t')]
C         for L'=1:nv, L=1:nv, t'=1:n2  nv = # partial waves at R = nlnmv(ic)
          zppb = ppb(nc+1:nc+nv,nc+1:nc+nv,i,icp) ! Copy portion of (real) ppb to complex temporary array
          call zgemm('T','N',nv,n2,nv,alpha,zppb,nv,cphiq(ias,1),nlmto,beta,zz,nlnmx)
C         if (ia == 2) then; print *, 'atom',ia; call zprm('psi2b zz',zz,nlnmx,nv,n2); endif
          do  itp = 1,n2
            do  jp = 1,nv
              zz(jp,itp) = dconjg(zz(jp,itp))
            enddo
          enddo

C         zpsipb(ib,t,t') = zwork(t,t') = Sum[L] [exp(ikT) cphiqk(L,t)]^T zz^T(L,t')
C         for product basis fun ib, L=1:nv, t=1:n1 nv = # partial waves at R = nlnmv(ic)
          call zgemm('T','N',n1,n2,nv,phase(ia),cphiqk(ias,1),nlmto,zz,nlnmx,beta,zwork,n1)
C         if (ia == 2) then; print *, 'atom',ia; call zprm('psi2b zpsipb',zwork,n1,n1,n2); endif
          zpsipb(ib,nctot+1:nctot+n1,ncc+1:ncc+n2)=zwork

        end do                  !end of product basis-loop
        deallocate(zppb)
      end do                    !end of atom-loop ia

      deallocate(zz,zwork,iasx)
      end

      subroutine psicb_v4(nctot,ncc,n1,n2,ndpsi1,ndpsi2,ndpsi3,
     .                   iclass,phase,
     .                   cphiqk,
     .                   cphiq,
     .                   ppb,
     .                   nlnmv,nlnmc,mdim,imdim,iatomp,
     .                   mdimx,nlmto,nlnmx,natom,nclass,
     .                   icore,ncore,nl,nnc,
     o                   zpsipb)
C- Matrix elements (product basis * valence * core) states
C ----------------------------------------------------------------------
Ci Inputs
Ci  nctot  : Offset in 2nd index zpsipb array: first element is zpsipb(:,off2+1,:)
Ci         : Typically off2 = number of core states.
Ci         : NEW off2
Ci  ncc    : Offset in 3rd index zpsipb array: first element is zpsipb(:,:,off3+1)
Ci         : At present ncc must be one of the following:
Ci         : 0     => do not calculate <core(k+q,t') |psi (k,t) B(R,ibloch)>
Ci         : nctot => do calculate <core(k+q,t') |psi (k,t) B(R,ibloch)>
Ci         : NEW off3
Ci  n1     : Polarizability context: no. occupied states at k, n in in Eq. 32, PRB76, 165106
Ci         : Self-energy context: no. intermediate states at q-k, n' in Eq. 34, PRB76, 165106
Ci         : OLD nt0
Ci  n2     : Polarizability context: no. unoccupied states at q+k (polarizability context)
Ci         : Self-energy context: no. states at q for which sigm calculated (n or m in Eq 34)
Ci         : OLD ntp0
Ci  ndpsi1,ndpsi2,ndpsi3 : parameters that dimension zpsipb
Ci  iclass : the jth atom belongs to class iclass(j)
Ci  phase  : exp(ik.T), where T it is translation Tiat, of site R to equivalent site under
Ci         : possible rotation of k to an equivalent kr (sxcf_fal3_scz)
Ci         : In the polarizability context, phase should be unity.
Ci  cphiqk : (polarizability context; see Remarks)
Ci         : amount of partial wave contributing to occ valence eigenfunction at k
Ci         : (self-energy context; see Remarks to psi2b_v4)
Ci         : amount of partial wave contributing to intermediate state at q-k
Ci  cphiq  : (polarizability context)
Ci         : amount of partial wave contributing to unocc valence eigenfunction at k+q
Ci         : (self-energy context)
Ci         : amount of partial wave contributing to eigenstate n or m at q (see Remarks to psi2b_v4)
Ci  ppb    : (real) matrix elements of partial waves <Phi(RLn) Phi(RL'n') B(R,i)>
Ci         : ppb includes both core and valence.  Sequence :
Ci         : (core at 1), (val at 1) (core at 2), (val at 2) ...
Ci  nlnmv  : For each class, number of partial waves for valence
Ci  nlnmc  : For each class, number of partial waves for core
Ci  mdim   : For each class, number of product basis functions
Ci         : NEW nrprdc
Ci  imdim  : imdim(i) = index to first product basis function for site i
Ci         : Product functions are strung together as a long vector.
Ci         : NEW: irprd
Ci  iatomp : mapping to site containing product basis?
Ci  mdimx  : dimensions ppb
Ci         : (NEW: nrprdcx)
Ci  nlmto  : dimension of lmto component of basis
Ci  nlnmx  : maximum number of product basis functions; dimensions ppb
Ci  natom  : number of sites
Ci  nclass : number of inequivalent classes
Ci  icore  : icore(i,ic) index in ppb array for core state i in class ic
Ci  ncore  : ncore(ic) no. core states in class ic
Ci  nl     : (global maximum l) + 1
Ci  nnc    : dimensions icore
Co Outputs
Co  zpsipb : matrix elements (see Remarks)
Cu Updates
Cu      Aug 14 Mark  wrote v4
Cu      Feb 06 Takao wrote v3
Cu      Apr 02 Takao wrote v2
Cu   17 Mar 92 Ferdi's original code
Cr Remarks
Cr Calculates  <psi (k+q,t') |core(k,t) B(R,ibloch)>, t'=unocc, t=occ
Cr   and also  <core(k+q,t') |psi (k,t) B(R,ibloch)> if ncc>0
Cr B(R,i) = Mixed basis
Cr core(k,t) = core states
Cr <psi(k+q,t') | core(k,t) B(R,i)> = Sum[RL]  cphiqk(RL,k+q,t')^* * ppb
C-------------------------------------------------------------------------
!      use mpi_mod, rank => mpi_rank
      implicit none
      integer i,ia,iap,ias,ib,ic,icp,icr,ics,it,itp,mdimx,natom,nc,nc1,
     .  ncc,nclass,nctot,nl,nlmto,nlnmx,nnc,n1,n2,nv,ndpsi1,ndpsi2,ndpsi3
      integer icore(nl*nl*nnc,nclass),ncore(nclass),nlnmv(nclass),
     .  nlnmc(nclass),mdim(nclass),iclass(natom),imdim(natom),iatomp(natom)
      real(8) :: ppb(nlnmx,nlnmx,mdimx,nclass)
C     complex(8):: zpsipb(ngb,nctot+n1,ncc+n2),cphiq(nlmto,*),cphiqk(nlmto,*),phase(natom)
      complex(8):: zpsipb(ndpsi1,ndpsi2,ndpsi3),cphiq(nlmto,*),cphiqk(nlmto,*),phase(natom)
!     integer(4):: verbose,ixx

      if (ncc/=0 .and. ncc/=nctot) call rx('psicb_v4: ncc must be either 0 or nctot')
      if (ncc /= 0 .and. sum(ncore(iclass(1:natom))) /= ncc) call rx('psicb_v4 : wrong number of core states')

!     zpsipb = 0d0
      call dinit(zpsipb,2*ndpsi1*ndpsi2*ndpsi3)

C --- For each site R, make psi2b for each product basis fn i, core state t, unocc t' ---
      ib  = 0
      ias = 1
      ics = 0
      do  ia = 1, natom
        ic  = iclass(ia)
        nc  = nlnmc(ic)
        nv  = nlnmv(ic)
        nc1 = nc + 1

C   --- Loop over product basis, core states t, unocc states t' ---
C       zpsipb(ib,t,t') = Sum[L] exp(ikT) cphiq(RL,t')^(*) * <Phi(RL) core(RL) B(R,i)>
C       Index  spans
C       ib     product basis functions of valence states 1:mdim(R)
C       L     partial valence partial waves 1:nlnmv(R), encompassing l,m, radial component
C       it     core states at R
C       itp    spans unocc states
C
C       If ncc is nonzero, also make
C       zpsipb(ib,t,t') = Sum[L] exp(-ikT) cphiqk(RL,t)^(*) * <Phi(RL) core(RL) B(R,i)>
C       it     spans occ valence states
C       itp    spans core states at R
C
        iap = iatomp(ia)
        icp = iclass(iap)
        ib  = imdim(iap)-1
C       MPI note: each i in this loop is independent ... can be parallelized over i
        do  i = 1, mdim(icp)
          ib = ib + 1
          do  itp = 1, n2
            do  it = 1, ncore(ic)
              icr = icore(it,ic)
              zpsipb(ib,ics+it,ncc+itp) = phase(ia)*dconjg(sum(cphiq(ias:ias+nv-1,itp)*ppb(nc1:nc+nv,icr,i,icp)))
            enddo               ! end of t(core)-loop
          enddo                 ! end of t'(unoccupied)-loop

          if (ncc==0) cycle     ! Skip <core(k+q,t') |psi (k,t) B(R,ibloch)>

          do  itp = 1, ncore(ic)
            icr = icore(itp,ic)
            do  it = 1, n1
              zpsipb(ib,nctot+it,ics+itp) = dconjg(phase(ia))*sum(cphiqk(ias:ias+nv-1,it)*ppb(nc1:nc+nv,icr,i,icp))

!          write(*,'("[",I4,"] writing to zpsipb at (",I3,",",I3,",",I3,")")') rank, ib, nctot+it, ics+itp

            enddo               ! end of t(occupied)-loop
          enddo                 ! end of t'(core)-loop

        enddo                   ! product basis-loop
        ias = ias + nlnmv(ic)
        ics = ics + ncore(ic)
      enddo                     ! loop over site
      end
