      subroutine strscr(loka,nbas,npadl,npadr,alat,plat,bas,rwats,nl,
     .  kap,nkap,lmaxw,siinv,cy,cg,indxcg,jcg,lpgf,ldot,leqv,ntab,iax,
     .  alpha,adot,tral,trad,s,sdot)
C- Calculates 2nd gen or NMTO structure constants for all sites
C ----------------------------------------------------------------
Ci Inputs
Ci   loka:       conventions for Hankels and Bessel functions
Ci   nbas:       number of atoms in the basis
Ci   nl:         maximum l quantum number+1
Ci   alat,plat:  primitive lattice translation vectors, and scale
Ci   bas:        basis vectors
Ci   alpha:      screening parameters
Ci   adot:       energy derivative of screening parameters
Ci   tral:       transformation matrix for head and tail functions
Ci   trad:       (kappa*avw)^2-derivative of tral
Ci   cg,indxcg,jcg: Clebsch Gordan coefficients
Ci   cy          Normalization of YL
Ci   ntab,iax:   neighbor table information (see Remarks)
Ci   kap         kappa's for which to calc strux, scaled by avw
ci   nkap        number of kappa's
Ci   ldot        T, make sdot
Ci   lmaxw:      l-cutoff for Watson sphere.
Ci   rwats       Watson sphere radii
Ci   siinv:      parameters for iterative inversion of strux
Ci               1: number of sites for cutoff
Ci               2: number of iterations
Ci               3: allowed tolerance
Ci   lpgf:       T, strx for planar Green's function (see Remarks)
Ci   ldot:       T, make sdot
Ci   nbasp:      Number of elements in the padded basis,
Ci               usually nbas except when lpgf=T
Ci   leqv:  T,   look for equivalent strx
Cio  iax(*,i):   information concerning pair i connecting R,R'
Ci   iax(1,*)    site index for 'head' site R'
Ci   iax(2,*)    site index for 'augmentation' site R
Ci   iax(3..5)   number of primitive lattice translation vectors
Ci               separating the unit cells corresponding to R and R'.
Ci   iax(6,*)    index to conjugate (jb,ib) pair matching (ib,jb).
Ci   iax(7,*)    permutation table ordering cluster by increasing
Ci               effective site index (see pairc.f; not used here)
Co   iax(8,*)    in the case strux are equivalent, points to the
Co               site it is equivalent to.
Ci   iax(9,i)    holds (lmx+1)**2 for basis atom iax(2,i) (mkiaxd)
Co Outputs
Co   s           screened strux in (nl**2,nl*2) blocks for all pairs
Co               1..ntab(nbas+1)
Co   sdot        avw^-2 times energy derivative of s
Cr Remarks
Cr strscr generates screened structure constants and sdot (nkap=1) for
Cr a specified neighbor table (nbas,ntab,iax).  Structure constants are
Cr strung together in one long vector s of nl**2 x nl**2 blocks for each
Cr pair in the iax table. For one pair s is stored as S(R'L',RL), the
Cr augmentation being at R'.
Cr
Cr ntab(ib) contains the total number of pairs in the table before site
Cr ib; ntab(ib+1)-ntab(ib) the number of pairs in ib.  iax contains the
Cr geometrical information for the cluster around each ib.  R and R' are
Cr separated by some primitive translation vectors, plus some difference
Cr of basis vectors within the unit cell:
Cr
Cr Convention for s:  let ib,jb be iax(1,isite) and iax(2,isite).
Cr Then s(L',L,isite) are the strux for a one-center expansion
Cr around site jb of a screened screened Hankel centered at ib:
Cr   H_L(r-r_b) = sum_La s(La,Lb,isite) J_La'(r-r_a)
Cr Here 'a' the the augmentation index, 'b' is the basis.
Cr
Cr Mixed Cholesky decomposition/ Dyson expansion:
Cr structure constants are obtained as follows.  Definining the vector
Cr beta as the solution to
Cr   (alpha^-1 - S^0)_uv beta_v = delta_u0, where
Cr the central atom is column zero (u,v composite RL indices), this
Cr is done for an inner block by Cholesky decomposition.  The outer
Cr block and inner block are then converged in a Dyson expansion by
Cr the pair of equations
Cr   beta_u =  \sum _v (alpha S^0)_uv beta_v
Cr for u belonging to outer block, and v summed over both blocks, and
Cr   beta_u = (alpha_uu \delta_uv - S^0)_uu^-1 (1 + S^0_uv beta_v)
Cr for u the inner block and v the outer block.
Cr the latter solved again using Cholesky decomposed matrix.
Cr These equations are iterated, each iteration adding a correction to
Cr the approximate inverse.
Cr
Cr Structure constants for layer geometry:
Cr   alpha, adot, bas must be doubly padded by the addition of
Cr   left and right bulk layers.
Cr
Cr Watson sphere:
Cr   See Andersen and Woolley, Mol. Phys. 26, 905 (1973).
Cr   A Watson sphere with radius rmx+deltr centered around the central
Cr   atom is introduced.  The region outside the Watson sphere is
Cr   called Omega.  Orbitals are introduced, which are regular outside
Cr   Omega, i.e., inside the Watson sphere; they are Bessel functions
Cr   with a l up to lmaxw.  These Bessel functions are expanded into
Cr   Bessel functions inside the spheres centered around all atoms in
Cr   the cluster.  Likewise the tails of the Neumann functions are
Cr   expanded into Neumann functions inside Omega.  Now screened
Cr   orbitals can be introduced which fulfill the same expansion rules
Cr   as the non-screened orbitals but with tails, which are now linear
Cr   combinations of Bessel and Hankel functions.  The coefficients of
Cr   mixing alpha and alphom are chosen so that all tails vanish at a
Cr   screening core radius sigma*wsr and on the Watson sphere.
Cr
Cr   Formally this is done by increasing the size of S^0 which now has
Cr   the structure:
Cr
Cr      ( S^0_{R',R}        S^0_{R',Omega}    )
Cr      (                                     )
Cr      ( S^0_{Omega,R}     S^0_{Omega,Omega} )
Cr
Cr   where S^0{R',R}         is the structure matrix Neumann-> Neumann
Cr         S^0_{R',Omega}    is the expansion matrix Bessel->Bessel
Cr         S^0_{Omega,R}     is the expansion matrix Neumann->Neumann
Cr         S^0_{Omega,Omega} is zero
Cr
Cu Updates
Cu   17 Jun 13 Replace f77 pointers with f90 ones
Cu   06 Aug 06 Redesigned so routines compatible with 2-kappa strux
Cu   19 May 04 make strux work with Methfessel conventions for H,J
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      logical lpgf,leqv,ldot
      integer nbas,npadl,npadr,nl,nkap,ntab(nbas+npadl+npadr+1),niax,lmaxw,loka
      parameter (niax=10)
      integer iax(niax,*)
      double precision kap(nkap),s(nl*nl,nl*nl,1),sdot(nl*nl,nl*nl,1)
      double precision alpha(nl**2,*),adot(nl**2,*),
     .  tral(4,nl**2,*),trad(4,nl**2,*),
     .  alat,plat(9),bas(3,nbas),rwats(nbas),siinv(5)
C For Clebsch Gordan coefficients:
      double precision cy(100),cg(1200)
      integer jcg(1200),indxcg(350)
C Local parameters
      logical cmdopt
      integer i,ik,lqinv,iat,ndimL,ndimW,nl2,ns,nsk,nitab,nitbck,
     .  isum,offik,nbasp,nbaspp,nsite0,nttab,nclus,iclus,streqv,stdo,
     .  ipr,nglob,nmto,nenv
C     integer os0a,os0d,osflg,osflgd
      double precision xx
      character*80 outs
      real(8), allocatable :: salph(:),alphv(:),adotv(:),sadot(:)
      real(8), allocatable :: s0a(:),s0d(:)
      integer, allocatable :: sflg(:),sflgd(:)

      stdo = nglob('stdo')
      call getpr(ipr)
      nbasp = nbas + npadl + npadr
      nbaspp = 2*nbasp - nbas
      nttab = ntab(nbasp+1)
      call rxx(.not.lpgf .and. nbas /= nbasp,'strscr: nbas ne nbasp')
      nl2 = nl**2
      nmto = 0
      if (nkap > 1) nmto = 1

C --- Temporarily scale the lattice and basis vectors ---
      call dscal(9,alat,plat,1)
      call dscal(3*nbaspp,alat,bas,1)

C --- For each kappa, do ---
      nsk = 0
      do  ik = 1, nkap
C ... nitab is number of strux actually in the table
      nitab  = 0
      nitbck= 0
      nsite0 = 0

      if (ipr > 40) write(stdo,'(1x)')

C --- Make the structure constants for all sites, this kappa ---
      offik = 1+(ik-1)*nbasp
      do  iat = 1, nbasp

        call info5(41,0,0,' strscr:  start strux for site %i'//
     .    '%?#(n>1)#, E*avw**2=%d##',iat,nkap,kap(ik),0,0)

C       Number of atoms in this cluster
        nclus = ntab(iat+1)-ntab(iat)
C       Offset to current cluster
        iclus = ntab(iat)+1
C       The true dimension of screened S-alpha
        ndimL = isum(nclus,iax(9,iclus),niax)
C       The dimension of strux matrix, including the Watson sphere
        ndimW = ndimL + (max(lmaxw,-1)+1)**2

C ...   Look for an already-calculated equivalent strux
        if (leqv) then
          if (nkap > 1) call rx('check equiv strux for kap>1')
          i = streqv(ntab,ldot,nl2,alpha,adot,iat,0,plat,bas,iax,nitab)
          if (ipr > 40 .and. i > 0)
     .      call awrit2(' strscr:  map strux for site%,4i  from site'//
     .      '  %i',' ',80,stdo,iat,i)
          if (i > 0) cycle
        else
          do  i = iclus, ntab(iat+1)
            nitab = nitab+1
            iax(8,i) = nitab
          enddo
        endif

C   ... Work arrays
        allocate(alphv(ndimW),salph(nl2*ndimW))
        if (ldot) then
          allocate(adotv(ndimW))
          allocate(sadot(nl2*ndimW)); call dpzero(sadot,nl2*ndimW)
        else
          allocate(adotv(1),sadot(1))
        endif

        lqinv = 0
        if (cmdopt('--qinv',6,0,outs)) lqinv = 1
C       if (cmdopt('--qinv',6,0,outs) .or. siinv(1) /= 0) lqinv = 1
        allocate(s0a(ndimW*ndimW))
        allocate(s0d(ndimW*ndimW))

C ...   Make vector of alpha
        call mkalpv(alpha(1,offik),adot(1,offik),
     .    iax(1,iclus),0,nl2,nbaspp,1,nclus,ndimW,alphv)
        if (ldot) call mkalpv(alpha(1,offik),adot(1,offik),
     .    iax(1,iclus),2,nl2,nbaspp,1,nclus,ndimW,adotv)
        if (lmaxw >= 0) call alwats(ndimL,kap(ik),2,lmaxw,
     .    rwats(iat),alphv,adotv)

C ...   Calculate salpha and/or salpha-dot
C        call salph1(loka,lqinv,ldot,nbaspp,nl,siinv,ndimW,lmaxW,nclus,
C     .    plat,bas,alphv,adotv,tral(1,1,offik),trad(1,1,offik),
C     .    iax(1,iclus),cy,cg,indxcg,jcg,kap(ik),s0a,s0d,
C     .    salph,sadot)
        nenv = iax(9,iclus)
        call salph1(nbaspp,nenv,lqinv,ldot,nmto*10+loka,siinv,ndimW,
     .    lmaxW,nclus,plat,bas,alphv,adotv,iax(1,iclus),cy,cg,
     .    indxcg,jcg,kap(ik),s0a,s0d,salph,sadot)

C   ... NMTO: scale salpha,sdot according to tral matrix
        if (nmto == 1) call scalsr(iax(1,iclus),ldot,ndimW,nl2,nclus,
     .    trad(1,1,offik),tral(1,1,offik),salph,sadot)
C       call prmx('balpha after scaling',salph,ndimW,ndimW,nenv)
C       call prmx('bdot after scaling',sadot,ndimW,ndimW,nenv)

C ...   Append salpha and/or salpha-dot to structure constant s
        ns = nitab-nclus
        call addtos(nl2,ndimW,iax(1,iclus),nl2,1,nclus,salph,ns+nsk,s)
        if (ldot) then
          ns = nitab-nclus
          call addtos(nl2,ndimW,iax(1,iclus),nl2,1,nclus,sadot,
     .      ns+nsk,sdot)
        endif

        if (iat <= nbas) nsite0 = nsite0 + nclus
        nitbck = nitbck + nclus
C       call wksav(.false.)
        deallocate(alphv,salph,adotv,sadot,s0a,s0d)

C ... Reentry for case equivalent strx
      enddo

C --- Sanity check ---
      if (leqv .and. nitbck /= nitab) call rx('bug in strscr')

C --- Symmetrize structure constants ---
      allocate(sflg(nitab)); call iinit(sflg,nitab)
      allocate(sflgd(nitab)); call iinit(sflgd,nitab)
C      if (lpgf) then
C        call sympad(nl2,nbas,npadl,npadr,nsite0,nitab,iax,1,
C     .    sflg,s(1,1,nsk))
C        if (ldot) call sympad(nl2,nbas,npadl,npadr,nsite0,nitab,iax,1,
C     .    sflgd,sdot)
C      endif
      call symstr(0,nl2,nttab,iax,1,1,sflg,s(1,1,1+nsk),s,xx)
      if (ldot)
     .call symstr(0,nl2,nttab,iax,1,1,sflgd,sdot(1,1,1+nsk),sdot,xx)
      nsk = nsk+ns+nclus
      deallocate(sflg,sflgd)
      enddo

C --- Undo scaling of the lattice and basis vectors ---
      call dscal(9,1/alat,plat,1)
      call dscal(3*nbaspp,1/alat,bas,1)

C --- Info printout ---
      if (ipr >= 30) call awrit2(' strscr:  generated %i '//
     .  'inequivalent strux from %i total',' ',80,stdo,nitab,nttab)

      end
