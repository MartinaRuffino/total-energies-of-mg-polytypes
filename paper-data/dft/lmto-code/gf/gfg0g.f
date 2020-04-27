      subroutine gfg0g(s_ham,s_pot,mode,isp,ib,jb1,jb2,
     .  offH,nbas,n1,n2,n3,k1,k2,k3,nlma,ldjj,gij,gjj)
C- Convert diagonal perturbed dg_ii = by energy scaling
C-----------------------------------------------------------------------
Cio Structures
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  lgen3 lncol ldham
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  gfdpp
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:dpf ddpf dddpf
Cio    Passed to:  gfdpp
Ci Inputs
Ci  mode   :1s digit
Ci         :0 perturbing potential is constant potential shift
Ci  ib     :site from which perturbation originated
Ci  jb1,jb2:scale subblock for sites in range jb1..jb2
Ci  offH   :Offsets to hamiltonian matrix (makidx.f)
Ci  n1..n3 :number of elements in gij,gjj
Ci  k1..k3 :leading dimensions of gij,gjj
Ci  nlma   :dimensions gjj and gji,and number of orbitals belonging to
Ci         :site ib for which to calculate g_jj
Ci  ldjj   :dimensions gij and gjj, and the number of orbitals
Ci         :connected to site ib for which gij,gjj are defined
Ci  gij    :Off-diagonal gf connecting orbitals i at ib to orbitals j
Cio Inputs/Outputs
Cio gjj    :On input gjj is unscale g; see Remarks
Cio        :On output gjj is scaled g
Cr Remarks
Cr   This hamiltonian-dependent routine scales the dimensionless
Cr   perturbed GF by energy scaling.  Input GF is change in diagonal
Cr   part of G owing to some perturbation at site ib.
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   15 Mar 00 extended to include iwaves and stubs for hwaves
C-----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer mode,isp,ib,jb1,jb2,offH(n0H,nkap0,1),nlma,ldjj,n1,n2,n3,
     .  k1,k2,k3,nbas
      double complex gij(k1,k2,k3,nlma,ldjj),gjj(k1,k2,k3,ldjj,nlma)
C ... For structures
!      include 'structures.h'
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
C ... Local parameters
      integer ldham(16),nlmi,ldim,ldimj,lgen3,lhdim,lidim,lncol,nspc,
     .  offpi,offpj,nlm0,isum,iblk
      equivalence (ldim,ldham(1)),(lidim,ldham(2)),(lhdim,ldham(3))
      equivalence (nspc,ldham(4))
      double complex dp(nlma,2,2),dpd(nlma,2,2),dpdd(nlma,2,2)

      call sanrg(.true.,mode,0,0,'gfg0g','mode')
      lgen3 = s_ham%lgen3
      lncol = s_ham%lncol
      ldham = s_ham%ldham
C     Never an independent spin channel in noncollinear case
      if (isp == 2) call sanrg(.true.,nspc,1,1,'gfg0g','nspc')

C --- 3nd generation NMTO --
      if (lgen3 /= 0) then
        call rx('gfg0g not ready for 3rd generation')

C --- 2nd generation ASA ---
      else
C       Both pf and gij are in downfolding order; thus no permutations

        call gfdpp(mode,s_ham,s_pot,offH,ib,nbas,isp,nlma,
     .    dp,dpd,dpdd)

        nlm0 = 0
        offpj = offH(1,1,jb1)
C       ldimj = offH(1,1,jb2+1) - offpj
C       l+i blocks for now
        ldimj = isum(2,offH(1,1,jb2+1),1)
        do  10  iblk = 1, 2
          nlmi   = offH(iblk,1,ib+1) - offH(iblk,1,ib)
          if (nlm0+nlmi > nlma) return
          offpi = isum(iblk-1,offH(1,1,nbas+1),1) + offH(iblk,1,ib)
          call gfg0g2(n1,n2,n3,k1,k2,k3,isp,nspc,nlmi,ldimj,offpi,offpj,
     .      lhdim,s_pot%dpf,dpd(1+nlm0,1,1),dpdd(1+nlm0,1,1),nlma,
     .      gij(1,1,1,1+nlm0,1),ldjj,gjj(1,1,1,1,1+nlm0))
          nlm0 = nlm0+nlmi
   10   continue

      endif

      end

      subroutine gfg0g2(n1,n2,n3,k1,k2,k3,isp,nspc,nlmi,ldimj,offpi,
     .  offpj,lhdim,pdot,dpd,dpdd,nlma,gij,ldjj,gjj)
C- Convert 2nd generation perturbation g_jj to proper perturbation G_jj
C ----------------------------------------------------------------------
Ci Inputs
Ci  n1,n2,n3:number of elements in gij,gjj
Ci  k1,k2,k3:leading dimensions of gij,gjj
Ci   isp   :current spin channel (1 or 2)
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   nlmi  :row dimension of gjj to scale ; see Remarks
Ci         :(number of orbitals defining the perturbation)
Ci   ldimj :column dimension of gjj to scale
Ci          (number of orbitals affected by the perturbation)
Ci   offpi :offset in pdot and row entry in gjj to orbitals generating
Ci         :perturbation; see Remarks.
Ci   offpj :offset to first entry in pdot corresponding to gjj(1)
Ci   lhdim :no. lower+intermediate+higher orbitals (dimensions pdot)
Ci   pdot  :energy derivative of potential function
Ci         :Both spin channels are present in the spin-polarized case
Ci         :but pdot is diagonal in the spins in the noncollinear case
Ci         :pdot must be ordered in the same way as the col dim of g.
Ci   dpd   :perturbation to pot. function P-dot at site ib
Ci         :dpd holds only one spin channel except in the
Ci         :noncollinear case, where it is a 2x2 spinor.
Ci   dpdd  :perturbation to pot. function (-1/2 P-dotdot/P-dot) at ib
Ci         :The spin parts of dpdd are the same as in dpd
Ci   nlma  :second dimension of gij
Ci   gij   :unscaled off-diagonal g for n1*n2*n3 qp
Ci         :The spinor parts of gij must match those of dpd
Ci   ldjj  :second dimension of gjj
Cio Inputs/Outputs
Cio  gjj   :On input gjj is unscaled diagonal part of pert. to g from ib
Cio        :On output gjj(i,j) is scaled to true GF
Cio        :The spin parts of gjj are the dimensioned as in dpd but the
Cio        :meaning of gjj(is1,is2) is the perturbation in g(is1,is1)
Cl Local variables
Cl   is1   :effective leading spin index to gii,gjj,dpd,dpdd
Cl         := 1 in collinear case, loops over spins 1..2 in noncollinear
Cl   is2   :effective second spin index to gii,gjj
Cl         := 1 in collinear case, loops over spins 1..2 in noncollinear
Cl   ks1   :isp in the collinear case, and is1 in the noncollinear
Cl   ks2   :isp in the collinear case, and is2 in the noncollinear
Cl         :Thus (is1,is2) = spin indices to arrays gii,gjj,dpd,dpdd
Cl         :ks1,ks2 = indices to spin parts of pdot (diagonal in spin)
Cl         :and gjj(j,is1,i,is2)  corresponds to change in
Cl         :g(j,is1,j,is1) from pot perturbation in channel (i,is2)
Cr Remarks
Cr   gfg0g2 scales dimensionless diagonal part of d (P-S)^-1 / dP(ib)
Cr   to dG by energy scaling.
Cr
Cr   The orbitals i=1..nlmi are the orbitals corresponding to the
Cr   perturbing potential.  These must correspond to rows
Cr   pdot(offpi+i), rows gjj(offpi+i,*), and columns gjj(*,i)
Cr   in pdot and gjj.  The orbitals generating the perturbation must
Cr   be contiguous; thus if all orbitals comprising the perturbation
Cr   are separated into subblocks, gfg0g2 must be called in succession
Cr   for each subblock.
Cb Bugs
Cb   Formulation a little wrong for noncollinear case.  As it is now,
Cb   program loops over two distinct perturbations (spin-up and -down
Cb   dP) scaling strictly diagonal part of g.  Should be called twice
Cb   with one distinct perturbation each time.  Thus, as now written
Cb   off-diagonal d(P^dot_i) is meaningless.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nlma,ldjj,nlmi,ldimj,offpi,offpj,lhdim,isp,nspc
      integer n1,n2,n3,k1,k2,k3
      double complex pdot(lhdim,isp),
     .  gij(k1,k2,k3,nlma,nspc,ldjj,nspc),dpd(nlma,nspc,nspc),
     .  gjj(k1,k2,k3,ldjj,nspc,nlma,nspc),dpdd(nlma,nspc,nspc)
C ... Local
      integer i,j,i1,i2,i3,is1,is2,ks1
C     double complex xxx(ldjj,nlma)

      do  5  is1 = 1, nspc
      do  5  is2 = 1, nspc
      ks1 = max(is1,isp)

C --- dg_jj(dV_i) <- sqrt(P^dot_j) dg_jj sqrt(P^dot_j) ---
      do  10  i = 1, nlmi
      do  10  j = 1, ldimj
      do  10  i3 = 1, n3
      do  10  i2 = 1, n2
      do  10  i1 = 1, n1
        gjj(i1,i2,i3,j,is1,i,is2) = gjj(i1,i2,i3,j,is1,i,is2) *
     .     pdot(j+offpj,ks1)
   10 continue

C --- g_jj <- g_jj + d (P^dot_i) g_ij for pert. in channels i ---
C     See Remarks about is1 == is2
      if (offpi >= 0 .and. offpi+nlmi <= ldjj .and. is1 == is2) then
        do  20  i = 1, nlmi
          do  22  i3 = 1, n3
          do  22  i2 = 1, n2
          do  22  i1 = 1, n1
            gjj(i1,i2,i3,i+offpi,is1,i,is2) =
     .      gjj(i1,i2,i3,i+offpi,is1,i,is2) +
     .        dpd(i,is2,is2)*gij(i1,i2,i3,i,is2,i+offpi,is2)
   22     continue

C     ... g_jj <- g_jj - delta (-1/2 P-dotdot/P-dot)
          gjj(1,1,1,i+offpi,is1,i,is2) =
     .    gjj(1,1,1,i+offpi,is1,i,is2) + dpdd(i,is1,is2)
   20   continue
      endif

C     debugging
C      do  j = 1, ldjj
C      do  i = 1, nlma
C        xxx(j,i) = gjj(1,1,1,j,is1,i,is2)
C      enddo
C      enddo
C      print *, 'xxx-scaled for is1,is2=',is1,is2
C      call zprm('Gji',2,xxx,ldjj,ldjj,nlma)

    5 continue

      end
