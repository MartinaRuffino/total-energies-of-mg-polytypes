      subroutine dglrft(mode,nk1,nk2,nk3,k1,k2,k3,nlmi,ldimi,ndim,gij,gji)
C- Delta-G by for potential perturbation at one site by FFT
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  : not used now
Ci   nk1,nk2,nk3 : number of k-points along the three
Ci   k1,k2,k3    : leading dimensions of gij,gji
Ci   nlmi  : number of orbitals belonging to site
Ci   ldimi : dimensions arrays gij,gji
Ci   ndim  : number of orbitals connected to site
Ci   gij,gji : unperturbed GF connecting orbitals j at the field
Ci         : point to orbitals i at the source point, or product
Ci         : GF * potential-perturbation.
Ci         : Dyson's equation is dg_ji = g_ji dV_i g_ij
Ci         : so either gji or gij should contain the product of g
Ci         : and the perturbation.
Ci         : Both gij and gji are in reciprocal space, full BZ.
Co Outputs
Co   gji   : is overwritten with convolution of gji and gij
Co         : gji is returned in reciprocal space.
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nk1,nk2,nk3,k1,k2,k3,nlmi,ldimi,ndim
      double complex gij(k1,k2,k3,nlmi,ldimi),
     .               gji(k1,k2,k3,ldimi,nlmi)
C ... Local parameters
      integer ilm,jd

      call tcn('dglrft')

C ... Real-space gij, for all nlmi*ndim channels
      call fftz3(gij(1,1,1,1,1),nk1,nk2,nk3,k1,k2,k3,nlmi*ndim,0,-1)

C ... Overwrite g_ij with FT of convolution,  g_ji~ g_ij~
      do  ilm = 1, nlmi

C       Real-space gji(ilm) for ndim channels
        call fftz3(gji(1,1,1,1,ilm),nk1,nk2,nk3,k1,k2,k3,ndim,0,-1)
        do  jd = 1, ndim

C       call zprm('gij',2,gij,nk1*nk2*nk3,nk1*nk2*nk3,1)
C       call zprm('gji',2,gji,nk1*nk2*nk3,nk1*nk2*nk3,1)

C         This is FT of convolution,  g_ji~ g_ij~
          call fftz3c(gji(1,1,1,jd,ilm),gij(1,1,1,ilm,jd),
     .                nk1,nk2,nk3,k1,k2,k3,33,-1)
        enddo

C       convolution from inverse FT of (FT of convolution)
        call fftz3(gji(1,1,1,1,ilm),nk1,nk2,nk3,k1,k2,k3,ndim,0,1)
      enddo
      call tcx('dglrft')
      end

      subroutine pdglrs(mode,s_ham,s_pot,
     .  offH,ib,nbas,isp,nspc,nk1,nk2,nk3,k1,k2,k3,ni,nlmj,gji)
C- Product of GF with perturbing potential at source point
C ----------------------------------------------------------------------
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
Ci   mode  :specifies the nature of the perturbing potential
Ci         :1 potential perturbation is a constant potential shift
Ci   offH  :Offsets to hamiltonian matrix (makidx.f)
Ci   ib    :site at which there is a perturbation
Ci   isp   :current spin channel (1 or 2)
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   nk1,nk2,nk3 : number of k-points along the three
Ci   k1,k2,k3    : leading dimensions of gij,gji
Ci   ni    :dimension of gji, and number of columns for which to gen. g
Ci   nlmj  :dimension of gji,and number of orbitals belonging to site
Ci          for which to calculate g_ji dP_i.
Co Inputs/Outputs
Cio  gji   : On input, unperturbed GF connecting orbitals i at the
Cio        : perturbation to orbitals j at the field point.
Cio        : On output, product of GF and perturbing potential,
Cio        : the latter depending on mode.
Cio        : The spin parts of gji are the same as in dp
Cl Local variables
Cl   dp    : perturbation to pot. function P-dot at site ib
Cl         : dp holds only one spin channel except in the
Cl         : noncollinear case, where it is a 2x2 spinor.
Cr Remarks
Cr   It assumed that orbitals 1..nlmj in gji are ordered in downfolding
Cr   order, because gfdpp generates dp for orbitals in that order.
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   15 Mar 00 extended to include iwaves, stubs for h-waves
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer mode,ib,nbas,isp,nspc,nk1,nk2,nk3,k1,k2,k3,ni,nlmj,
     .  offH(n0H,nkap0,1)
      double complex gji(k1,k2,k3,ni,nspc,nlmj,nspc)
C ... For structures
!      include 'structures.h'
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
C ... Local parameters
      integer id,jd,i1,i2,i3,is,js,ks
      double complex dp(nlmj,nspc,nspc),dpd(nlmj,nspc,nspc),
     .  dpdd(nlmj,nspc,nspc),xx(2,2)


      call gfdpp(mode,s_ham,s_pot,offH,ib,nbas,isp,nlmj,dp,dpd,dpdd)

      if (nspc == 1) then
        do  10  id = 1, nlmj
        do  10  jd = 1, ni
        do  10  i3 = 1, nk3
        do  10  i2 = 1, nk2
        do  10  i1 = 1, nk1
          gji(i1,i2,i3,jd,1,id,1) = gji(i1,i2,i3,jd,1,id,1)*dp(id,1,1)
   10   continue
      else
        do  20  id = 1, nlmj
        do  20  jd = 1, ni
        do  20  i3 = 1, nk3
        do  20  i2 = 1, nk2
        do  20  i1 = 1, nk1
          xx(1,1) = gji(i1,i2,i3,jd,1,id,1)
          xx(1,2) = gji(i1,i2,i3,jd,1,id,2)
          xx(2,1) = gji(i1,i2,i3,jd,2,id,1)
          xx(2,2) = gji(i1,i2,i3,jd,2,id,2)
          do  22  is = 1, 2
          do  22  js = 1, 2
          gji(i1,i2,i3,jd,is,id,js) = 0
          do  22  ks = 1, 2
            gji(i1,i2,i3,jd,is,id,js) = gji(i1,i2,i3,jd,is,id,js) +
     .                                  xx(is,ks)*dp(id,ks,js)
   22     continue
   20   continue
      endif

      end
