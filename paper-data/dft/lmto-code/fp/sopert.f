      subroutine sopert3(mode,s_ham,ndimh,h,s,hso)
C- Diagonalize hamiltonian, adding SO part perturbatively
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit (see description of input h,s)
Ci         :0 s(:,:,1:2) = overlap in orbital basis
Ci            h(:,:,1:2) = H0+LzSz in orbital basis
Ci         :1 s(:,:,1:2) = evecs of H0+LzSz
Ci            h(:,:,1:2) = H0+LzSz in H0+LzSz repsn
Ci         :10s digit (see description of output h,s)
Ci         :1 L+S- included nonperturbatively.
Ci            evals are not divided by spin
Ci            evecs are 2x2 matrices in spin space
Ci         :3 L+S- included approximately, similar to 2nd order pert th.
Ci            evals are still divided by spin
Ci            spin parts of evecs are not coupled
Ci         :4 The full LS is included approximately, in the manner of mode 3
Ci   ndimh :dimension of hamiltonian and overlap h,s for one spin block
Ci   h     :1s digit mode = 0 :
Ci          h(1:ndimh,1:ndimh,isp) spin diagonal parts of h
Ci          h(1:ndimh,1:ndimh,3)   L+S- part of h
Ci         :1s digit mode = 1 :
Ci          h(1:ndimh,1:ndimh,isp) = H0+LzSz in H0+LzSz repsn (h is diagonal)
Ci          h(1:ndimh,1:ndimh,3)   L+S- part of h
Ci   s     :1s digit mode = 0 :
Ci         :s(1:ndimh,1:ndimh,isp=1,2) spin diagonal parts of s
Ci         :1s digit mode = 1 :
Ci         :1 s(:,:,1:2) = evecs of H0+LzSz
Co Outputs
Co   h     :10s digit mode = 1 : h(1:ndimh,1:2,1:ndimh,1:2) =
Co          H0+L.S rotated to basis of eigenfunctions (diagonal)
Co         :10s digit mode = 3 : h(1:ndimh,1:ndimh,isp) =
Co          H0+LzSz+L+S-(pert) rotated to basis of H0+LzSz eigenfunctions
Co   s     :10s digit mode = 1 : s(1:ndimh,1:2,1:ndimh,1:2) evecs of H0+L.S
Co         :10s digit mode = 3 : s(1:ndimh,1:ndimh,isp) = evecs of H0+LzSz
Cs Command-line switches
Cl Local variables
Cl         :
Cb Bugs
Cb   SO=3 case does not work well when SO coupling is very large, e.g. in Pu
Cb   Some data, e.g. mxevcut, should be passed back, or
Cb   special entry point made to return it
Cr Remarks
Cr   When 10s digit mode is 1, output h and input h are dimensioned differently.
Cu Updates
Cu   18 Aug 14 New mode 4
Cu   08 Sep 12 Orders eigenvalues
Cu   06 Aug 12 Adapted from hambls
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,ndimh
      double complex h(ndimh,ndimh,4),s(ndimh,ndimh,4) ! Only 3 for pert
      double complex hso(ndimh,ndimh,4) ! Used if mode=4
C ... For structures
!      include 'structures.h'
      type(str_ham):: s_ham
C ... Local parameters
      integer NULLI,nspc
      logical T,F
      parameter (T=.true.,F=.false.,nspc=2,NULLI=-99999)
      double precision eomin
      integer,save:: mxevcut=0,mnevl=NULLI
      integer ispc,nevl,ndimhx,mode0,mode1,i,j,ovncut,mxevcuti(2)
      double precision xx,vij2,de0,split,epsovl
      double complex zer,one
      parameter (zer=(0d0,0d0),one=(1d0,0d0))
      integer,allocatable :: iprm(:,:)
      real(8),allocatable:: e(:,:),de(:,:)
      complex(8),allocatable:: zqi(:,:),z0(:,:)

      ndimhx = 2*ndimh
      allocate(z0(ndimh,ndimh),zqi(ndimh,ndimh),e(ndimh,2))
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
C     call sanrg(.true.,mode1,0,1,'sopert3','1s digit mode')

C --- Overwrite h,s so that:
C     s(:,:,1:2) = evecs of H0+LzSz
C     h(:,:,1:2) = evals of H0+LzSz (h is diagonalized)
C     or,  if 10s digit mode = 4
C     s(:,:,1:2) = evecs of H0
C     h(:,:,1:2) = evals of H0 (h is diagonalized)
      if (mode0 == 0) then
      do  ispc = 1, 2

        if (mode1 == 4) then
          call daxpy(2*ndimh**2,-1d0,hso(1,1,ispc),1,h(1,1,ispc),1)
        endif

C       call zprm('hii',2,h(1,1,ispc),ndimh,ndimh,ndimh)
C       call zprm('hso',2,h(1,1,3),ndimh,ndimh,ndimh)

C   ... zqi, e <- evecs, evals of H0+LzSz) = Z.lda
        epsovl = s_ham%oveps
        ovncut = s_ham%ovncut
        if (epsovl == 0 .and. ovncut == 0) then
          call zhevx(ndimh,ndimh,h(1,1,ispc),s(1,1,ispc),1,T,ndimh,9d9,
     .      nevl,z0,F,e(1,ispc),ndimh,zqi)
          mxevcuti(ispc) = 0
        else
          nevl = -1
          allocate(de(ndimh,1))
          de(1,1) = 99999d0
          call zhevo(ndimh,ndimh,h(1,1,ispc),s(1,1,ispc),ndimh,9d9,
     .      epsovl,ovncut,nevl,i,e(1,ispc),de,ndimh,zqi)
          eomin = de(1,1)
          mxevcuti(ispc) = ndimh-nevl
          if (mnevl == NULLI) mnevl = nevl
          mnevl = min(nevl,mnevl)
          deallocate(de)
        endif
C       call zprm('evecs',2,zqi,ndimh,ndimh,ndimh)
C       if (ispc == 2) call yprm('e, h+LzSz',1,e,1,ndimh,ndimh,2)

        call dpzero(h(1,1,ispc),ndimh**2*2)
        do  i = 1, ndimh
          h(i,i,ispc) = e(i,ispc)
        enddo
        call zcopy(ndimh**2,zqi,1,s(1,1,ispc),1)
      enddo
      mxevcut = max(mxevcut,mxevcuti(1)+mxevcuti(2))
      else
        do  ispc = 1, 2
          do  i = 1, ndimh
            e(i,ispc) = h(i,i,ispc)
          enddo
        enddo
      endif

C --- Rotate L+S- part of H to H0+LzSz basis ---
C     h(:,:,3)  = L+S- in H0+LzSz repsn
*     call zprm('z11',2,s(1,1,1),ndimh,ndimh,ndimh)
*     call zprm('z22',2,s(1,1,2),ndimh,ndimh,ndimh)
      if (mode1 == 1 .or. mode1 == 3) then ! h(:,:,3)  = L+S- in H0+LzSz repsn
        call phmbls(1,ndimh,ndimh,e,xx,zqi,h(1,1,3),
     .    s(1,1,1),s(1,1,2),h(1,1,3))
*       call zprm('S.O, evec repsn',2,h(1,1,3),ndimh,ndimh,ndimh)
      endif

C --- Eigenvalues,eigenvectors of full SO hamiltonian: H0+L.S ---
      if (mode1 == 1) then

        deallocate(z0,zqi)
        allocate(z0(ndimhx,ndimhx),zqi(ndimhx,ndimhx))

C   ... zqi <- evecs of H0+LzSz, spin-coupled form
        call sopert(1,ndimh,nspc,z0,s,zqi)

C   ... h <- (H0+LzSz) + L+S- in zqi basis, spin-coupled form
        call sopert(0,ndimh,nspc,z0,h,h)
*       call zprm('h+LzSz+L+S-, h+LzSz rep',12,h,ndimhx,ndimhx,ndimhx)

C   ... Diagonalize full h in the H0+LzSz representation
        call zhevx(ndimhx,ndimhx,h,xx,0,T,ndimhx,9d9,
     .    i,s,F,e,ndimhx,z0)
*       call yprm('evals, h+LzSz+L+S-',1,e,1,ndimhx,ndimhx,1)
*       call zprm('evecs, h+LzSz rep',2,z0,ndimhx,ndimhx,ndimhx)

C   ... Rotate evecs to orbital basis = z(so) z(orbital)
C       Given zqi = evecs of H0+LzSz (orbital basis)
C       and   z0  = evecs of H0+LzSz+L+S-, H0+LzSz basis
C       Then e = z0+ zqi+ (H0+L.S) zqi z0
        call zgemm('N','N',ndimhx,ndimhx,ndimhx,one,zqi,
     .    ndimhx,z0,ndimhx,zer,s,ndimhx)
*       call zprm('evecs of h+L.S',2,s,ndimhx,ndimhx,ndimhx)

C       Check projection onto one or another spin channel
C        allocate(iprm(ndimh,2),de(ndimh,2))
C        do i = 1, ndimh
C          iprm(i,1) = i
C        enddo
C        call suqlsz(ndimhx,1,ndimh,iprm,ndimh,s,de)

C        call phmbls(32+64,ndimhx,ndimhx,e,xx,zqi,xx,s,s,h)
C        call zprm('Reconstruct original h+L.S',12,h,ndimhx,ndimhx,ndimhx)

C   ... Replace h with evals on the diagonal
        call dpzero(h,ndimhx**2*2)
        call dcopy(ndimhx,e,1,h,2*ndimhx+2)
*       call zprm('evals of h+L.S',2,h,ndimhx,ndimhx,ndimhx)

C --- shifts in evals from L.S or L+S- only ---
      else if (mode1 == 3 .or. mode1 == 4) then

        allocate(iprm(ndimh,2),de(ndimh,2))
        call dpzero(de,2*ndimh)

        if (mode1 == 3) then
          do  i = 1, ndimh
          do  j = 1, ndimh
            vij2 = dconjg(h(i,j,3))*h(i,j,3)
            if (vij2 < 1d-12) cycle
C           vij2 = vij2*.6d0

            de0 = (e(i,1)-e(j,2))/2
            split = dsqrt(de0**2 + vij2)
C           (a+b)/2 +/- split = a - de0 +/- split = b + de0 +/- split
C           If a>b, de0>0.  shift = split-de0 for a, -split+de0 for b
            if (de0 >= 0) then
              de(i,1) = de(i,1) +split-de0
              de(j,2) = de(j,2) -split+de0
C           If b>a, de0<0.  shift = -split+de0 for a, -split+de0 for b
            else
              de(i,1) = de(i,1) -split-de0
              de(j,2) = de(j,2) +split+de0
            endif
          enddo
          enddo

        else ! mode1 == 4

C     ... zqi <- evecs of H0, spin-coupled form
          deallocate(z0,zqi)
          allocate(z0(ndimhx,ndimhx),zqi(ndimhx,ndimhx))
          call sopert(1,ndimh,nspc,z0,s,zqi)
*         call zprm('evecs of H0 (2x2)',2,zqi,ndimhx,ndimhx,ndimhx)

C     ... z0 <- L.S in spin-coupled form (use z0 because appropriately dimensioned)
          call sopert(10,ndimh,nspc,z0,hso,xx)
C         call zprm('L.S, orbital rep',12,z0,ndimhx,ndimhx,ndimhx)

C     ... z0 <- <z|L.S|z>
          call phmbls(1,ndimhx,ndimhx,xx,xx,hso,z0,zqi,zqi,z0)
C         call zprm('<z|L.S|z>',2,z0,ndimhx,ndimhx,ndimhx)

          do  i = 1, ndimhx
          do  j = 1, ndimhx
            if (i == j) then  ! First order perturbation (should be 0)
C              de(i,1) = de(i,1) + z0(i,i)
              cycle
            endif

C            if (i <= ndimh .or. j <= ndimh) cycle  ! just 11 block
C            if (i <= ndimh .and. j <= ndimh .or.
C     .          i > ndimh .and. j > ndimh) cycle  ! just 12,21, block

            vij2 = dconjg(z0(i,j))*z0(i,j)
            if (vij2 < 1d-12) cycle

            de0 = (e(i,1)-e(j,1))/2
            split = dsqrt(de0**2 + vij2)
C           (a+b)/2 +/- split = a - de0 +/- split = b + de0 +/- split
C           If a>b, de0>0.  shift = split-de0 for a, -split+de0 for b
C           Split is half as much as SO3 case because include both (1,2) and (2,1) blocks
            if (de0 >= 0) then
              de(i,1) = de(i,1) +split/2-de0/2
              de(j,1) = de(j,1) -split/2+de0/2
C           If b>a, de0<0.  shift = -split+de0 for a, -split+de0 for b
            else
              de(i,1) = de(i,1) -split/2-de0/2
              de(j,1) = de(j,1) +split/2+de0/2
            endif
          enddo
          enddo
        endif

C       call yprm('pert shifts',1,de,1,ndimh,ndimh,2)
        call dpadd(e,de,1,2*ndimh,1d0)
C       call yprm('unsorted final evals',1,e,1,ndimh,ndimh,2)
C       call zprm('unsrt evecs spin 1',2,s(1,1,1),ndimh,ndimh,ndimh)
C       call zprm('unsrt evecs spin 2',2,s(1,1,2),ndimh,ndimh,ndimh)
        do  ispc = 1, 2
          call dvheap(1,ndimh,e(1,ispc),iprm(1,ispc),0d0,101)
          call dvprm(1,ndimh,e(1,ispc),xx,iprm(1,ispc),2)
          call dvprm(2*ndimh,ndimh,s(1,1,ispc),xx,iprm(1,ispc),2)
        enddo
C       call yprm('sorted evals',1,e,1,ndimh,ndimh,2)
C       call zprm('srt evecs spin 1',2,s(1,1,1),ndimh,ndimh,ndimh)
C       call zprm('srt evecs spin 2',2,s(1,1,2),ndimh,ndimh,ndimh)
        call dpzero(h,ndimh**2*2*2)
        do  ispc = 1, 2
          do  i = 1, ndimh
            h(i,i,ispc) = e(i,ispc) ! + de(i,ispc)
          enddo
        enddo

        deallocate(de,iprm)

      else
        call rxi('sopert3 : mode not recognized: ',mode)
      endif

      deallocate(z0,zqi,e)

      end

      subroutine sopert(mode,ndimh,nspc,wk,hin,hout)
C- Manipulates blocks of hamiltonian for noncollinear case
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 orders spin spin subblocks of h into full matrix
Ci         :1 does the same, but the (1,2) spin block is set to zero
Ci         :Optional 10s digit to be added
Ci         :10 copy result into wk; then hout not accessed
Ci   ndimh :hamiltonian dimension
Ci   nspc  :routine does nothing if nspc ne 2
Ci   wk    :double complex work array of dimension (ndimh*nsp)**2
Ci   hin   :blocks of hamiltonian
Co Outputs
Co   hout  :(mode 0) hamiltonian blocks ordered in form suitable for
Co         :         for diagonalization.
Co         :         hout and hin may occupy the same address space.
Co         :           part of hin    spin block of h
Co         :            hin(*,*,1)       (1,1)
Co         :            hin(*,*,2)       (2,2)
Co         :            hin(*,*,3)       (1,2)
Co         :(mode 1) same as mode 0 but (1,2) block assumed to be 0
Cl Local variables
Cl         :
Cr Remarks
Cr Checks:
Cr   Combine s11,s22 into large matrix s
Cr   mc s22 s11 -sub 1,nr+nr,1,nc+nc -suba nr/2+1,nr/2+1 -herm -bw s
Cu Updates
Cu   05 Feb 05 Added mode 1
Cu   23 Dec 04 First created (mode 0)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,ndimh,nspc
      double complex hin(ndimh,ndimh,*),hout(ndimh,2,ndimh,2)
      double complex wk(ndimh,2,ndimh,2)
C ... Local parameters
      integer ndimhx,i,j,ks,is,js,mapki(3),mapkj(3),mode0
      data mapki /1,2,1/, mapkj /1,2,2/

      if (nspc /= 2) return

C     n2 = ndimh**2
      ndimhx = ndimh*nspc
      mode0 = mod(mode,10)
      if (mode0 == 0 .or. mode0 == 1) then
        do  ks = 1, 3
          is = mapki(ks)
          js = mapkj(ks)
          if (is == js .or. mode0 == 0) then
            do  j = 1, ndimh
            do  i = 1, ndimh
              wk(i,is,j,js) = hin(i,j,ks)
            enddo
            enddo
C           (2,1) block: hermitian conjugate of (1,2) block
            if (ks == 3) then
              do  j = 1, ndimh
              do  i = 1, ndimh
                wk(j,js,i,is) = dconjg(hin(i,j,ks))
              enddo
              enddo
            endif
          else
            do  j = 1, ndimh
            do  i = 1, ndimh
              wk(i,is,j,js) = (0d0,0d0)
              wk(i,js,j,is) = (0d0,0d0)
            enddo
            enddo
          endif
        enddo
        if (mode < 10) call dcopy(ndimhx**2*2,wk,1,hout,1)
C       call zprm('h(nc)',2,hout,ndimhx,ndimhx,ndimhx)

      else
        call rxi('sopert: not implemented, mode=',mode0)
      endif

      end

      subroutine sosite(mode,s_site,s_spec,s_lat,s_ham,s_bz,q,
     .  ndimh,nev,e,z,wtkp,ewgt,tso)
C- Resolves SO coupling by augmentation site and/or band index
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos sighh sighk sigkk pihh pihk pikk sighhx
Ci                 sighkx sigkkx pihhx pihkx pikkx sohh sohk sokk
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  mkhso augmbl bstrux
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa lmxb kmxt rsma pz name orbp
Co     Stored:     orbp
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  mkhso augmbl bstrux uspecb
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat qlat vol plat awald tol nkd nkq
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:napw igv2 cg indxcg jcg cy qlv dlv
Cio    Passed to:  mkhso augmbl bstrux hxpbl ghibl hklbl gklbl hxpgbl
Cio                ghigbl hklgbl
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:lrsa iprmb
Cio    Passed to:  mkhso augmbl
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  sopertp w n
Co     Stored:     sopertp
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   mode  :not used now
Ci   q     :Bloch wave number
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   ndimh :dimension of h and s
Ci   e     :eigenvalues of hamiltonian h+hso
Ci   z     :eigenvectors of hamiltonian h+hso
Ci   ewgt  :weights for occupied states
Cio Inputs/Outputs
Cio sopertp:parameters used to estimate Fermi surface contribution to change in band sum.
Cio         (1) input    sampling Fermi level
Cio         (2) not used sampling sumev
Cio         (3) not used DOS(Ef)
Cio         (4) output   sum [delta w] sum of change in weights (charge) from
Cio                      1st order perturbation theory is accumulated for this qp
Cio         (5) output   Same as (4), but delta w from 2nd order pert theory
Cio         (6) output   sum [delta (w e)] change in weight-eval product, sampling
Cio                      is accumulated for this qp
Cio         (7) output   sum [w (delta e)] weighted change in evals by sampling
Cio                      is accumulated for this qp
Cio         (8) output   sum [w (delta e)] (weighted change in evals by tetrahedron
Cio                      is accumulated for this qp
Cio         (9) output   Same as (6), but delta e from 2nd order pert theory
Cio        (10) output   Same as (7), but delta e from 2nd order pert theory
Cio        (11) output   Same as (8), but delta e from 2nd order pert theory
Cio        (12) output   d(w)/d(Ef)
Co Outputs
Co   tso   :Site-resolved SO coupling accumulated for this q
Co         :tso(1,1:2,1:2,0) : LS++, LS+-, LS-+, LS-- for all sites
Co         :tso(1,1:2,1:2,i) : LS++, LS+-, LS-+, LS-- for site i
Co         :tso(2,1  ,1  ,0) : Sum of four components of LS, all sites
Co         :tso(2,1  ,1  ,0) : Sum of four components of LS, site i
Cl Local variables
Cr Remarks
Cr   Second order pert. theory when eigenvectors are given to all orders:
Cr   Let E0k and |k0> be UNPERTURBED evals and evecs of eigenstate |k>
Cr   and let |k1> be 1st order correction to |k>.
Cr   From 2nd order perturbation theory, change in energy of state n:
Cr     dE2n = <n0|V|n0> + sum'_k <n0|V|k0> <k0|V|n0> / (E0n-E0k)   (1)
Cr   The prime indicates that k ranges over all states other than n.
Cr   V is the perturbation; the last term is the 2nd order term.  Also
Cr     |n1> = sum'_k |k0> <k0|V|n0> / (E0n-E0k)
Cr
Cr   Suppose we want to calculate 2nd order correction dE2n from
Cr   TRUE eigenstates |k>, instead of unperturbed |k0>.
Cr   If know only <n|V|k>, write perturbation in terms of it.
Cr   The key point is that the 1st order term has now embedded
Cr   in it some contribution from the second order.
Cr   The first order term can be written as follows:
Cr      <n0|V|n0> = <n0-n+n|V|n0-n+n>
Cr      <n0|V|n0> - <n|V|n> = <n0-n|V|n> + <n|V|n0-n> + O(n0-n)^2
Cr
Cr   Omit the last term: it introduces corrections beyond 2nd order.
Cr   Note also that n0-n = -(n1+n2+n3+...) = -n1 to lowest order.
Cr   Keeping lowest order only:
Cr     <n0|V|n0> - <n|V|n> = - <n1|V|n0> - <n0|V|n1>
Cr   The second term on the rhs is:
Cr     -<n0|V|n1> = -<n0| V sum'_k |k0> <k0|V|n0> / (E0n-E0k)
Cr                = - sum'_k <n0|V|k0> <k0|V|n0> / (E0n-E0k)
Cr   which exactly the same as the 2nd order term in (1), but
Cr   with opposite sign.  The first term on the rhs evaluates to the same.
Cr   There are also corrections to the 2nd order term, but they are of
Cr   third order and higher.
Cr   Thus the second order perturbation correction is similar to before
Cr   except that the sign of the 2nd order term is changed:
Cr     dE2n = <n|V|n> - sum'_k <n|V|k> <k|V|n> / (En-Ek)   (2)
Cu Updates
Cu   05 Oct 13 L.S is resolved into ++, +-, -+, -- contributions
Cu   08 Aug 13 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nev,ndimh
      double precision q(3),wtkp,ewgt(nev,3),e(nev),tso(2,2,2,0:*)
      double complex z(ndimh,2,nev)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_ham)::   s_ham
      type(str_bz)::    s_bz
C ... Dynamically allocated arrays
      complex(8), allocatable :: t(:,:,:,:),wk(:,:,:,:)
      complex(8), allocatable :: zhz(:,:),szhz(:,:)
C ... Local parameters
      integer i,ia,i1,i2,j,nbas,ndimhx,nglob,nsp,nspc,samn,is
      double precision t1(nev),t2(nev)  ! first and 2nd order perturbation corrections
      double precision vij2,de0,split,ddot,x,efermis,samw,d,smrh,xx,
     .  wt0,wt1,wt2
      double complex zer,one
      parameter (zer=(0d0,0d0),one=(1d0,0d0))
C      integer mpipid
C ... External calls
      external augmbl,daxpy,dcopy,dpadd,dpzero,dvheap,dvprm,phmbls,rxi,
     .         sopert,zcopy,zgemm,zhevo,zhevx

C --- Setup ---
      nbas  = nglob('nbas')
      nsp = nglob('nsp')
      nspc = nglob('nspc')
      ndimhx = 2*ndimh
      efermis = s_bz%sopertp(1)   ! Fermi energy from sampling integration
      samw    = s_bz%w            ! Methfessel-Paxton sampling width
      samn    = s_bz%n            ! Methfessel-Paxton sampling order
      if (nsp /= 2 .or. nspc /= 2) return

      allocate(t(ndimh,2,ndimh,2),wk(ndimh,2,ndimh,2))
      allocate(zhz(nev,nev),szhz(nev,nev))
      do  ia = 0, nbas ! ia=0 captures entire SO
        if (ia > 0) then
          is = s_site(ia)%spec
          if (s_spec(is)%lmxa < 0) cycle
        endif

C        call dpzero(t,2*ndimh*ndimh*4)
C        otso = ndimh**2*2
C        do  i1 = 1, 2  ! Assemble diagonal t(i1,i1) and t(1,2)
C          oti = ndimh**2*(i1-1)
C          call augmbl(13,s_site,s_spec,s_lat,s_ham,i1,1,q,ia,
C     .      ndimh,napw,igapw,t(1+oti,1,1,1),t(1+otso,1,1,1),s)
C        enddo
C        call sopert(0,ndimh,nspc,wk,t,t) ! Render as full t(ndimh,2,ndimh,2)
        call mkhso(1,s_site,s_spec,s_lat,s_ham,q,ndimh,ia,t)
C        call zprm('z',2,z,ndimhx,ndimhx,nev)
C        call zprm('t',2,t,ndimhx,ndimhx,ndimhx)

C       Make z+ t z
C        call zgemm('C','N',nev,ndimhx,ndimhx,one,z,ndimhx,t,ndimhx,
C     .    zer,wk,ndimhx)
C        call zprm('wk',2,wk,ndimhx,ndimhx,ndimhx)
C        call zgemm('N','N',nev,nev,ndimhx,one,wk,ndimhx,z,ndimhx,
C     .    zer,zhz,nev)
C        call zprm('zhz',2,zhz,nev,nev,nev)

C       Make z+_i1 t_i1,i2 z_i2
        call dpzero(szhz,nev*nev*2)
        do  i1 = 1, 2
          do  i2 = 1, 2
C           wk(i1,i2) = z+_i1 t_i1,i2
            call zgemm('C','N',nev,ndimh,ndimh,one,z(1,i1,1),ndimhx,
     .        t(1,i1,1,i2),ndimhx,zer,wk,ndimhx)
C           call zprm('wk',2,wk,ndimhx,nev,ndimh)
C           (z+_i1 t_i1,i2) z_i2
            call zgemm('N','N',nev,nev,ndimh,one,wk,ndimhx,
     .        z(1,i2,1),ndimhx,zer,zhz,nev)
C           call zprm('zhz',2,zhz,nev,nev,nev)
C      ...  First order perturbation correction (decompose by blocks)
            do  i = 1, nev
              t1(i) = zhz(i,i)
            enddo
            tso(1,i1,i2,ia) = tso(1,i1,i2,ia) + ddot(nev,t1,1,ewgt,1)
            call daxpy(nev*nev*2,1d0,zhz,1,szhz,1)
          enddo ! Loop over i1
        enddo ! Loop over i2

C  ...  First order perturbation correction
        do  i = 1, nev
          t1(i) = szhz(i,i)
        enddo

C   ... Second order
        call dpzero(t2,nev)
        do  i = 1, nev
          do  j = i, nev
            if (i == j) cycle
            vij2 = dconjg(szhz(i,j))*szhz(i,j)
            if (vij2 < 1d-12) cycle

            de0 = (e(i)-e(j))/2
            split = dsqrt(de0**2 + vij2)
C           (e1+e2)/2 +/- split = e1 - de0 +/- split
C                               = e2 + de0 +/- split
C           e1 and e2 are already handled in 1st order
C           If ei>ej, de0>0.  e(i) shifts up, e(j) shifts down.  2nd order - 1st order :
            if (de0 >= 0) then
              t2(i) = t2(i) + split - de0
              t2(j) = t2(j) - split + de0
C           If ei<ej, de0<0.  e(i) shifts down, e(j) shifts up.  2nd order - 1st order :
            else
              t2(i) = t2(i) - split - de0
              t2(j) = t2(j) + split + de0
            endif
          enddo

        enddo

C     ... Accumulate Fermi surface contributions
          if (ia == 0) then
          do  i = 1, nev

C           Sampling weight wt0 at given energy e(i)
            x = (efermis - e(i)) / samw
            call delstp(samn,x,d,smrh,xx)
            wt0 = abs(wtkp) * (1d0 - smrh)

C           sopertp(12) <- sum d(wt0)/d(Ef)
            x = (efermis + 1d-4 - e(i)) / samw
            call delstp(samn,x,d,smrh,xx)
            wt1 = abs(wtkp) * (1d0 - smrh)
            x = (efermis - 1d-4 - e(i)) / samw
            call delstp(samn,x,d,smrh,xx)
            wt2 = abs(wtkp) * (1d0 - smrh)
            s_bz%sopertp(12) = s_bz%sopertp(12) + (wt1-wt2)/2d-4  ! d(w)/d(Ef) ~ D(E)
C           s_bz%sopertp(12) = s_bz%sopertp(12) +(wt1+wt2-2*wt0)/1d-4**2 ! D'(E)

C           Sampling weight at 1st order perturbed energy
            x = (efermis - (e(i)+t1(i))) / samw
            call delstp(samn,x,d,smrh,xx)
            wt1 = abs(wtkp) * (1d0 - smrh)

C           Sampling weight at 2nd order perturbed energy
            x = (efermis - (e(i)+t1(i)-t2(i))) / samw
            call delstp(samn,x,d,smrh,xx)
            wt2 = abs(wtkp) * (1d0 - smrh)

C           Fermi surface terms
            s_bz%sopertp(4) = s_bz%sopertp(4) + (wt1-wt0) ! 1st order change in net charge
            x = t1(i)
            s_bz%sopertp(6) = s_bz%sopertp(6) + (e(i)+x)*wt1 - e(i)*wt0 ! Delta(wt*evl)
            s_bz%sopertp(7) = s_bz%sopertp(7) + x*wt0 ! wt * Delta(evl), Sampling
            s_bz%sopertp(8) = s_bz%sopertp(8) + x*ewgt(i,1) ! wt * Delta(evl), tet

            s_bz%sopertp(5) = s_bz%sopertp(5) + (wt2-wt0) ! 2nd order change in net charge
            x = t1(i)-t2(i)
            s_bz%sopertp(9) = s_bz%sopertp(9) + (e(i)+x)*wt2 - e(i)*wt0 ! Delta(wt*evl)
            s_bz%sopertp(10) = s_bz%sopertp(10) + x*wt0 ! wt * Delta(evl), Sampling
            s_bz%sopertp(11) = s_bz%sopertp(11) + x*ewgt(i,1) ! wt * Delta(evl), tet
          enddo
          endif

C       print *, ddot(nev,t1,1,ewgt,1)
C       print *, ddot(nev,t2,1,ewgt,1)
C       Perturbation is 1st order - 2nd order; see Remarks
        tso(2,1,1,ia) = tso(2,1,1,ia) +
     .    ddot(nev,t1,1,ewgt,1) - ddot(nev,t2,1,ewgt,1)

C       print *, sum(tso(1,:,:,ia)), tso(2,1,1,ia)
      enddo ! Loop over sites

C        xx = s_bz%sopertp(6)-s_bz%sopertp(4)*efermis*0-s_bz%sopertp(7)
C        call info5(20,0,0,' FS: '//
C     .  '<tet>%;13,8D  <sam>%;13,8D    wt corr%;13,8D'//
C     .  '  <tet>+wt corr%;13,8D',
C     .  s_bz%sopertp(8),s_bz%sopertp(7),xx,s_bz%sopertp(8)+xx,0d0)
C        print *, s_bz%sopertp(6),s_bz%sopertp(4),s_bz%sopertp(7)
C
C        xx = s_bz%sopertp(9)-s_bz%sopertp(5)*efermis*0-s_bz%sopertp(10)
C        call info5(20,0,0,' FS2 '//
C     .  '<tet>%;13,8D  <sam>%;13,8D    wt corr%;13,8D'//
C     .  '  <tet>+wt corr%;13,8D',
C     .  s_bz%sopertp(11),s_bz%sopertp(10),xx,s_bz%sopertp(11)+xx,0d0)
C        print *, s_bz%sopertp(9),s_bz%sopertp(5),s_bz%sopertp(10)

      deallocate(t,wk,zhz,szhz)

      end

      subroutine mkhso(mode,s_site,s_spec,s_lat,s_ham,q,ndimh,ia,hso)
C- Make SO part of hamiltonian
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos sighh sighk sigkk pihh pihk pikk sighhx
Ci                 sighkx sigkkx pihhx pihkx pikkx sohh sohk sokk
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  augmbl bstrux
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa lmxb kmxt rsma pz name orbp
Co     Stored:     orbp
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  augmbl bstrux uspecb
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat qlat vol plat awald tol nkd nkq
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:napw igv2 cg indxcg jcg cy qlv dlv
Cio    Passed to:  augmbl bstrux hxpbl ghibl hklbl gklbl hxpgbl ghigbl
Cio                hklgbl
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:lrsa iprmb
Cio    Passed to:  augmbl
Ci Inputs
Ci   mode  :1s digit controls how hso is stored
Ci         :0  hso = hso(ndimh,ndimh,1:3)
Ci         :   spin block   stored in
Ci         :      (1,1)     hso(:,:,1)
Ci         :      (2,2)     hso(:,:,2)
Ci         :      (1,2)     hso(:,:,3)
Ci         :1  hso = hso(ndimh,2,ndimh,2)
Ci         :      (1,1)     hso(:,1,:,1)
Ci         :      (2,2)     hso(:,2,:,2)
Ci         :      (1,2)     hso(:,1,:,2)
Ci         :      (2,1)     hso(:,2,:,1)
Ci         :10s digit mode controls whether SO is copied or added to hso:
Ci         :0  SO copied into hso
Ci         :1  SO added  into hso (cannot be used with 1s digit mode=1 now)
Ci   q     :Bloch wave number
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   ndimh :dimension of h and s
Ci   ia    :<=0     => contributions to hso from all sites
Ci         :1..nbas => contributions to hso from site ia only
Co Outputs
Co   hso   :L.S contribution to the hamiltonian
Cr Remarks
Cu Updates
Cu   14 Aug 14 Adapted from sopert
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,ndimh,ia
      double precision q(3)
      double complex hso(ndimh,ndimh,4) ! Use only 3 if mode=0
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_ham)::   s_ham
C ... Dynamically allocated arrays
      complex(8), allocatable :: wk(:,:,:,:)
      integer, pointer :: igapwl(:,:)
C ... Local parameters
      integer isp
      double precision qs(3),dlength
      double complex s(1,1)
      real(8), parameter :: tol=1d-8
      integer, parameter :: nsp=2

      igapwl => s_lat%igv2
      call shorps(1,s_lat%qlat,(/72,2,2/),q,qs)
      if (dlength(3,q-qs,1) > tol .and. s_ham%pwmode > 10) then
        allocate(igapwl(3,s_lat%napw))
        call shorigv(q,qs,s_lat%plat,s_lat%napw,s_lat%igv2,igapwl)
      endif

      if (mod(mode/10,10) == 0) call dpzero(hso,2*ndimh*ndimh*3)
      do  isp = 1, 2 ! hso(:,:,1:3) <- hso(1,1), hso(2,2), and hso(1,2)
        call augmbl(13,s_site,s_spec,s_lat,s_ham,isp,1,qs,ia,ndimh,
     .    s_lat%napw,igapwl,hso(1,1,isp),hso(1,1,3),s)
      enddo
      if (.not. associated(igapwl,s_lat%igv2)) deallocate(igapwl)

C     call zprm('hso11',12,hso(1,1,1),ndimh,ndimh,ndimh)
C     call zprm('hso22',12,hso(1,1,2),ndimh,ndimh,ndimh)
C     call zprm('hso12',2,hso(1,1,3),ndimh,ndimh,ndimh)

C ... Render hso(:,:,1:3) as hso(:,i,:,j)
      if (mod(mode,10) == 0) return
      call rxx(mod(mode/10,10) /= 0,'mkhso not set up for this mode')
      allocate(wk(ndimh,2,ndimh,2))
      call sopert(0,ndimh,nsp,wk,hso,hso) ! Render as full hso(ndimh,2,ndimh,2)
      deallocate(wk)

C     call zprm('hso',12,hso,ndimh*2,ndimh*2,ndimh*2)

      end

