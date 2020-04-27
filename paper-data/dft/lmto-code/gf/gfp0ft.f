      subroutine gfp0ft(s_ctrl,s_ham,s_pot,s_lat,iscr,isx,lemom,
     .  isp,offH,iprmb,nbas,ibas,pos,nkp,qp,zp,wz,nk1,nk2,nk3,ipq,
     .  igstar,ifac,qb,nRLc,dGd,dGdirr,dGdk,dGdkr,gii)
C- Bare response function P0 by FT technique
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  lham lncol nl nspin lgen3
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  hamfbz
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  ldham hord lgen3 lncol
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  gfp0f1 pdglrs gfdpp gfg0g
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:dpf ddpf dddpf
Cio    Passed to:  gfp0f1 pdglrs gfdpp gfg0g
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  plat
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:istab symgr ag
Cio    Passed to:  *
Ci Inputs
Ci   iscr  :if 1, calculate dGd
Ci   isx   :determines contraction and which of dGdirr, dGdk, dGdkr is accumulated
Ci         :1s digit
Ci         :0 do not return any k-resolved response functions
Ci         :1 return one of the following k-resolved response functions
Ci         :10s digit
Ci         :0 : accumulate dGdirr
Ci         :1 : accumulate dGdk
Ci         :2 : accumulate real part of dGdk into dGdkr
Ci         :3 : accumulate imaginary part of dGdk into dGdkr
Ci         :100s digit
Ci         : 0  contract over l and m
Ci         : 1  contract over m
Ci         : 2  do not contract
Ci         :1000s digit (collinear case only)
Ci         : 1  always poke result into 1st spin (ks2=1 below)
Ci   lemom :T, accumulate gii; F, do not
Ci   isp   :current spin channel (1 or 2)
Ci   offH  :Offsets to hamiltonian matrix (makidx.f)
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   nbas  :size of basis
Ci   ibas  :calculated exchange coupling to sites ibas
Ci   pos   :basis vectors (input)
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   qp    :irreducible k-points
Ci   zp    :complex energy
Ci   wz    :weights for complex energy integration
Ci   nk1,nk2,nk3:  no. divisions for the 3 recip. latt. vecs
Ci   ipq   :ipq(i1,i2,i3) points to the irreducible qp into which
Ci          mesh point (i1,i2,i3) is mapped (bzmesh.f)
Ci   igstar:contains info needed to map an irreducible qp to its
Ci          original point (bzmesh.f)
Ci   ifac  :used with qb to recover a wave vector from three indices
Ci          specifying position in a uniform mesh (bzmsh0.f)
Ci   qb    :vectors of a microcell in the Brillouin zone
Ci   nRLc  :leading dimension of dGd
Co Outputs
Co   Which of the following is generated depends on iscr and isx
Co   dGd    :q=0 linear response in the diagonal GF channel j
Co          :from a constant potential shift dV_j in channel i
Co   dGdirr :same as dGd, but accumulated for all irreducible q-points
Co   dGdk   :same as dGd, but accumulated for all q-points
Co   dGdkr  :Real or imaginary part of dGdk is accumulated
Co          :(Imaginary part of dGd is the susceptibility; see Remarks)
Co   If lemom is true, make:
Co   gii   :energy-weighted moments of the sphere charges
Cr Remarks
Cr  This routine calculates the linear response dGd(i,j) in the energy
Cr  integrated diagonal part of the Green' function for each channel i
Cr  from a constant potential shift dV_j in channel j, for each pair of
Cr  channels (i,j).  It is accomplished via the Dyson equation
Cr      dG_ii = sum_j G0_ij dV_j G0_ji
Cr  For a perturbation in a single channel j,
Cr      dG_ii / dV_j  = G0_ij G0_ji
Cr  The bare susceptibility is defined as the change in partial density
Cr  n_i from a potential shift dV_j in channel j:
Cr      chi_ij = dn_i / dV_j
Cr  Using
Cr    n_i = -1/pi Im int^(E_f) G_ii(E) then
Cr    dn_i/dV_j = -1/pi Im int^(E_f) d G_ii(E) / dV_j
Cr              = -1/pi Im int^(E_f) dE (d G_ii / dV_j )
Cr              = -1/pi Im int^(E_f) dE G0_ij G0_ji
Cr
Cr  Sum rule:  Partial density-of-states is
Cr    D_i(E) = -1/pi G_ii
Cr  If there is a uniform potential shift (dV_j same for all j)
Cr    dn_i / dV = dn_i / dE = D_i(E)
Cr  So
Cr    sum_j chi_ij =  D_i(E)
Cr  delta G is obtained from delta g by rescaling (gfg0g), and
Cr  dGd is obtained by accumulating a sum of delta G foreach energy
Cr  zp with weighting wz.  suhamz should be called previously for each
Cr  energy to set potential parameters.
Cr
Cu Updates
Cu   08 May 13 Complete migration to f90 structures; eliminate s_array
Cu   10 Nov 11 Begin migration to f90 structures
Cu   13 Jan 09 First cut at other kinds of response functions
Cu   24 May 02 First cut at making q-dependent response function
Cu   15 Mar 00 extended to include iwaves, stubs for h-waves
C ----------------------------------------------------------------------
      use structures
      implicit none
C Passed variables
      logical lemom
      integer nbas,nkp,nk1,nk2,nk3,ipq(*),ifac(3),ibas,igstar(0:1),iscr,
     .  isx,nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer offH(n0H,nkap0,nbas+1),iprmb(*),nRLc
      double precision qp(3,nkp),gii(*),pos(3,*),qb(3,3)
      double complex zp,wz,dGd(nRLc,nRLc,2),dGdirr(nRLc,nRLc,nkp,2)
      double complex dGdk(nk1,nk2,nk3,nRLc,nRLc,2)
      double precision dGdkr(nk1,nk2,nk3,nRLc,nRLc,2)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_lat)::   s_lat
C Local variables
      integer hord,ib,ib1,ib2,isp,k1,k2,k3,ldham(16),ldim,ldimx,lgen3,
     .  lham,lhdim,lhdimx,lidim,lidimx,lio,lncol,nkfbz,nl,nlma,nsp,nspc,
     .  nnrl,offc,offr,iblk,nlmaa,opt1
      double precision plat(3,3),xx
      equivalence (ldim,ldham(1)),(ldimx,ldham(5)),(nspc,ldham(4))
      equivalence (lidim,ldham(2)),(lidimx,ldham(6))
      equivalence (lhdim,ldham(3)),(lhdimx,ldham(7))

      complex(8), pointer :: gij(:),gji(:),ghh(:)

C ... For file I/O of gf
      integer clp(9,2),fopnx,ifi,iogfrs
      integer MNC,MNB,MNT,MLB,MPL,MZP,MCD,ldh,mxorb,nglob
      parameter (MNC=1,MNB=2,MNT=4,MLB=8,MPL=16,MCD=32,MZP=64)

      if (mod(iscr,2) == 0 .and. mod(isx,2) == 0) return

      call tcn('gfp0ft')

      lham = s_ctrl%lham
      lncol = s_ctrl%lncol
      nl = s_ctrl%nl
      nsp = s_ctrl%nspin
      lgen3 = s_ctrl%lgen3
      plat = s_lat%plat
      ldham = s_ham%ldham
      hord = s_ham%hord
      nkfbz = nk1*nk2*nk3
      mxorb = nglob('mxorb')
      call fftz30(nk1,nk2,nk3,k1,k2,k3)

C     Sanity check for dimension of dGdk
      if (mod(isx,10) /= 0) then
        opt1 = mod(isx/100,10)
        iblk = nnrl(2-opt1,1,nbas,iprmb,lidim)
        call sanrg(.true.,nRLc,iblk,iblk,'gfp0ft (bug):','nRLc')
      endif

C     List of sites
      if (ibas == 0) then
        ib1 = 1
        ib2 = nbas
      else
        ib1 = ibas
        ib2 = ibas
      endif

C --- For each site, do ---
      iblk = 0
      do  ib = ib1, ib2

C --- gf for all qp in subblock connected to site ib ---
C ... Get gf for subblock connected to site ib
      nlma  = offH(4,1,ib+1) - offH(4,1,ib)
      lidim = offH(4,1,nbas+1)
      ldh   = offH(3,1,nbas+1) - offH(3,1,ib)
C ... If iblk<= ib, increment iblk and get gf for subblock ib,iblk
      if (iblk < ib) then
        iblk = min(ib+3,ib2)
        nlmaa = offH(4,1,iblk+1) - offH(4,1,ib)

        allocate(gij(k1*k2*k3*nlmaa*lidim*nspc**2))
        allocate(gji(k1*k2*k3*lidim*nlmaa*nspc**2))
        allocate(ghh(ldh*mxorb))
        offr = 0
        offc = 0

C       Rewind file containing gf, read header
        call iinit(clp,9*2)
        clp(3,1) = lidim*nspc
        clp(4,1) = lidim*nspc
        clp(5,1) = lidim*nspc
        ifi = fopnx('gfqp',100,16+8+4+0,-1)
        rewind ifi
        call pshpr(1)
        lio = 2 + 10*(MNC*1+MNB*1+MNT*0+MLB*1+MPL*0+MCD*1+MZP*0)
        if (iogfrs(lio,0,0,' ',ifi,1,1,nbas,xx,zp,qp,plat,xx,clp,xx,xx,
     .    0,0,0,xx) /= 0) call rx('gfp0ft: failed to read file header')
        call poppr

C        call gffbz(s_ctrl,nbas,offH,iprmb,ifi,ib,iblk,440,pos,nkp,qp,zp,
C     .    nk1,nk2,nk3,k1,k2,k3,ipq,plat,s_lat%istab,s_lat%symgr,s_lat%ag,igstar,
C     .    ifac,qb,offr,offc,ldh,gij,gji,ghh)

        call hamfbz(s_ctrl,s_ham,nbas,ifi,ib,iblk,440,pos,
     .    nkp,qp,zp,nk1,nk2,nk3,k1,k2,k3,ipq,plat,s_lat%istab,
     .    s_lat%symgr,s_lat%ag,igstar,ifac,qb,offr,offc,ldh,gij,gji,ghh)

      endif

      call gfp0f1(s_ham,s_pot,iscr,isx,lemom,offH,iprmb,ib,nbas,isp,nsp,
     .  nspc,wz,nk1,nk2,nk3,k1,k2,k3,ipq,igstar,nkp,nlma,offr,gij,gji,
     .  nRLc,lidim,gii,dGd,dGdirr,dGdk,dGdkr)

      offr = offr + nlma*lidim*nspc**2
      if (ib == iblk) deallocate(gij,gji,ghh)

C     call zprm('dGd',2,dGd(1,1,1),nRLc,nRLc,nRLc)

C ... Scale sum_q g_ii to return BZ integral of gii
      enddo
      if (lemom) call dscal(lidimx**2*2,1/dble(nkfbz),gii,1)
C     call yprm('sum_q g',2,gii,lidim*lidim,lidim,lidim,lidim)

      call tcx('gfp0ft')
      end

      subroutine gfp0f1(s_ham,s_pot,iscr,isx,lemom,offH,iprmb,ib,nbas,
     .  isp,nsp,nspc,wz,nk1,nk2,nk3,k1,k2,k3,ipq,igstar,nkp,nlma,offr,
     .  gr,gc,nRLc,lidim,gii,dGd,dGdirr,dGdk,dGdkr)
C- Kernel of gfp0ft that accumulates dGd and gii for one site
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  lgen3 lncol ldham
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  pdglrs gfdpp gfg0g
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:dpf ddpf dddpf
Cio    Passed to:  pdglrs gfdpp gfg0g
Ci Inputs
Ci   iscr  :if mod(iscr,2)=1, calculate dGd
Ci   isx   :determines contraction and which of dGdirr, dGdk, dGdkr is accumulated
Ci         :1s digit
Ci         :0 do not return any k-resolved response functions
Ci         :1 return one of the following k-resolved response functions
Ci         :10s digit
Ci         :0 : accumulate dGdirr
Ci         :1 : accumulate dGdk
Ci         :2 : accumulate real part of dGdk into dGdkr
Ci         :3 : accumulate imaginary part of dGdk into dGdkr
Ci         :100s digit
Ci         : 0  contract over l and m
Ci         : 1  contract over m
Ci         : 2  do not contract
Ci   lemom :T, accumulate gii; F, do not
Ci   offH  :Offsets to hamiltonian matrix (makidx.f)
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   ib    :site for which to make perturbation
Ci   nbas  :size of basis
Ci   isp   :current spin channel (1 or 2)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   wz    :weights for complex energy integration
Ci   nk1,nk2,nk3:  no. divisions for the 3 recip. latt. vecs
Ci   k1,k2,k3:  leading dimensions of gr,gc
Ci   nlma  :dimensions gr,gc
Ci   gr    :GF g_i,*(k) for orbitals i in site ib
Ci         :The number of col dimensions * is fixed by 100s digit mode
Ci   gc    :GF g_*,i(k) for orbitals i in site ib
Ci         :The number of row dimensions * is fixed by 100s digit mode
Ci   nRLc  :leading dimension of dGd
Ci   lidim :number of lower+intermediate orbitals
Co Outputs
Co   dGd   :linear response in the diagonal GF from constant pot. shifts
Co   gii   :(lemom = T) diagonal GF
Cl Local variables
Ci   nlml  :number of orbitals in lower ham. block for current site
Ci   offl  :offset to start of lower ham. block for current site
Ci   nlmi  :number of orbitals in int. ham. block for current site
Ci   offi  :offset to start of int. ham. block for current site
Cr Remarks
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   13 Jan 09 First cut at other kinds of response functions
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      logical lemom
      integer ib,nbas,isp,nsp,nspc,nlma,lidim,nk1,nk2,nk3,k1,k2,k3,
     .  nRLc,offr,iscr,isx,nkp,ipq(1),igstar(0:1)
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer offH(n0H,nkap0,nbas+1),iprmb(*)
      double complex wz,dGd(nRLc,nRLc,nsp),dGdirr(nRLc,nRLc,nkp,nsp)
      double complex gr(nk1,nk2,nk3,nlma,nspc,lidim,nspc),
     .               gc(nk1,nk2,nk3,lidim,nspc,nlma,nspc)
      double complex dGdk(nk1,nk2,nk3,nRLc,nRLc,2)
      double precision dGdkr(nk1,nk2,nk3,nRLc,nRLc,2)
      double precision gii(lidim,nspc,lidim,nspc)
C ... For structures
!      include 'structures.h'
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
C ... Local parameters
      integer nkfbz,nlml,offl,nlmi,offi,ldim,i1,i2,i3

C     call zprm('gc',2,gc(1,1,1,1+offr,1,1,1),nk1*nk2*nk3,1,lidim*nlma)

      nkfbz = nk1*nk2*nk3
      ldim =  offH(1,1,nbas+1)
      nlml  = offH(1,1,ib+1) - offH(1,1,ib)
      offl  = offH(1,1,ib)  ! Lower block
      nlmi  = offH(2,1,ib+1) - offH(2,1,ib)
      offi  = ldim + offH(2,1,ib)  ! Offset to intermediate block

C ... Accumulate gii
      if (lemom) then
        call pvgfec(nkfbz,nlml,lidim,nlma,lidim,nspc,lidim,
     .    gr(1,1,1,1+offr,1,1,1),gii(1+offl,1,1,1))
        call pvgfec(nkfbz,nlmi,lidim,nlma,lidim,nspc,lidim,
     .    gr(1,1,1,1+offr+nlml,1,1,1),gii(1+offi,1,1,1))
      endif

C ... Multiply potential perturbation into gji
      call pdglrs(0,s_ham,s_pot,offH,ib,nbas,isp,nspc,
     .  nk1,nk2,nk3,k1,k2,k3,lidim,nlma,gc(1,1,1,1+offr,1,1,1))

C ... Bare response function P0 for this ib, dimensionless (P-S) repsn
      call dglrft(0,nk1,nk2,nk3,k1,k2,k3,nlma*nspc,lidim*nspc,
     .  lidim*nspc,gr(1,1,1,1+offr,1,1,1),gc(1,1,1,1+offr,1,1,1))

C   99 print *, 'i,j='
C      read(*,*) i,j
C      call zprm('gc',2,gc(1,1,1,i,1,j,1),nk1*nk2*nk3,nk1*nk2*nk3,1)
C      goto 99

C ... Scale to proper GF
      if (mod(isx,2) == 1) then
        i1 = nk1
        i2 = nk2
        i3 = nk3
      else
        i1 = 1
        i2 = 1
        i3 = 1
      endif
      call gfg0g(s_ham,s_pot,0,isp,ib,1,nbas,offH,nbas,i1,i2,i3,k1,k2,
     .  k3,nlma,lidim,gr(1,1,1,1+offr,1,1,1),gc(1,1,1,1+offr,1,1,1))
C     call zprm('gc',2,gc(1,1,1,1+offr,1,1,1),nk1*nk2*nk3,1,lidim*nlma)

C ... Contract over m index; sum into dGd and/or dGdirr
      if (mod(iscr,2) == 1)
     .  call gfp0f2(k1*k2*k3,ib,nbas,isp,nspc,iprmb,lidim,nlma,
     .  gc(1,1,1,1+offr,1,1,1),(3-nsp)*wz,nRLc,dGd)
      if (mod(isx,2) == 1)
     .  call gfp0f3(isx/10,nk1,nk2,nk3,k1,k2,k3,nkp,ipq,igstar,ib,nbas,
     .  isp,nspc,iprmb,lidim,nlma,gc(1,1,1,1+offr,1,1,1),(3-nsp)*wz,
     .  nRLc,dGdirr,dGdk,dGdkr)

C       print *, 'end of gfp0f1',ib,isp,
C     .  sum(dgdkr(1,1,1,1:nRLc,1:nRLc,isp))

      end

      subroutine gfp0f2(k123,ib,nbas,isp,nspc,iprmb,lhdim,nlma,gji,
     .  wz,nRLc,dGd)
C- Contract linear response g_R'+T,L',i over T,m for one site
C ----------------------------------------------------------------------
Ci Inputs
Ci   k123  :leading dimension of gji
Ci   ib    :site for which to accumulate response function
Ci   nbas  :size of basis
Ci   isp   :current spin channel (1 or 2)
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   lhdim :size of l+i+h hamiltonian block
Ci   nlma  :size of basis for
Ci   gji   :diagonal gf for orbital j responding to potential shift at i
Ci   wz    :weights for complex energy integration
Ci   nRLc  :leading dimension of dGd
Co Outputs
Co   dGd   :change in the diagonal GF from perturbation at site ib
Co          is accumulated for this energy
Cl Local variables
Cl   lmrci :compound (R,l) index for site ib, (contracted over mi)
Cl   lmrcj :compound (R,l) index for site jb, (contracted over mj)
Cl   is1   :effective leading spin index to gii,gjj,dpd,dpdd
Cl         := 1 in collinear case, loops over spins 1..2 in noncollinear
Cl   is2   :effective second spin index to gii,gjj
Cl         := 1 in collinear case, loops over spins 1..2 in noncollinear
Cl   ks2   :isp in the collinear case, and is2 in the noncollinear
Cl
Cl         :Thus (is1,is2) = spin indices to array gji
Cl         :     (is1,ks2) = spin indices to array dGd
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ib,nbas,isp,nspc,k123,lhdim,nlma,nRLc,iprmb(*)
      double complex dGd(nRLc,nspc,nRLc,isp),wz,
     .               gji(k123,lhdim,nspc,nlma,nspc)
C ... Local parameters
      double precision pi
      double complex wt
      integer jb,offi,lj,jd,li,id,lmrci,lmrcj,nnrl,is1,is2,ks2,n0,nkap0
      parameter (n0=10,nkap0=4)
      integer norbi,norbj,ixx,iorb,jorb,nlmi,nlmj,
     .  ltabi(n0*nkap0),ktabi(n0*nkap0),offli(n0*nkap0),
     .  ltabj(n0*nkap0),ktabj(n0*nkap0),offlj(n0*nkap0)

      do  10  is1 = 1, nspc
      do  10  is2 = 1, nspc
      ks2 = max(is2,isp)

C ... wt is -1/pi * energy weight
      pi = 4*datan(1d0)
      wt = wz/pi

C ... Find starting lmrci index = 1st column in dGd to store
      lmrci = 1+nnrl(1,1,ib-1,iprmb,lhdim)

C --- For each orbital in site ib, contract over m_i gji ---
C     Loop over row index with contiguous l,i,h blocks in dnfold order.
      call orbl(ib,0,lhdim,iprmb,norbi,ltabi,ktabi,ixx,offli,nlmi)
      offi = 0
      do  12  iorb = 1, norbi
      li = ltabi(iorb)
      nlmi = 2*li + 1
      do  14  id = offi+1, offi+nlmi

C   --- For each orbital jd, store and contract over m_j,m_i gji ---
        do  20  jb = 1, nbas
        call orbl(jb,0,lhdim,iprmb,norbj,ltabj,ktabj,ixx,offlj,nlmj)
        lmrcj = 1+nnrl(1,1,jb-1,iprmb,lhdim)
        do  22  jorb = 1, norbj
        lj = ltabj(jorb)
        nlmj = 2*lj + 1
        do  24  jd = offlj(jorb)+1, offlj(jorb)+nlmj

C     ... Contract over mj,mi, accumulate into dGd
          dGd(lmrcj+lj,is1,lmrci+li,ks2) =
     .    dGd(lmrcj+lj,is1,lmrci+li,ks2) + wt*gji(1,jd,is1,id,is2)

   24   continue
   22   continue
   20   continue
   14 continue
      offi = offi + nlmi
   12 continue
   10 continue
      end

      subroutine gfp0f3(opt,nk1,nk2,nk3,k1,k2,k3,nkp,ipq,igstar,ib,nbas,
     .  isp,nspc,iprmb,lidim,nlma,gji,wz,nRLc,dGdirr,dGdk,dGdkr)
C- Contract linear response g_R'+T,L',i over T,m for one site
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :determines contraction and which of dGdirr, dGdk, dGdkr is accumulated
Ci         :1s digit
Ci         :0 : accumulate dGdirr
Ci         :1 : accumulate dGdk
Ci         :2 : accumulate  Re(dGdk) into dGdkr
Ci         :3 : accumulate -Im(dGdk) into dGdkr
Ci         :10s digit
Ci         : 0  contract over l and m
Ci         : 1  contract over m
Ci         : 2  do not contract
Ci         :100s digit (collinear case only)
Ci         : 1  always poke result into 1st spin (ks2=1 below)
Ci   k123  :leading dimension of gji
Ci   igstar:information needed for rotation qp to star of q
Ci         :bzmesh should be called with igstar(0)=-2
Ci   ib    :site for which to accumulate response function
Ci   nbas  :size of basis
Ci   isp   :current spin channel (1 or 2)
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   lidim :size of l+i hamiltonian block
Ci   nlma  :size of basis for
Ci   gji   :diagonal gf for orbital j responding to potential shift at i
Ci   wz    :weights for complex energy integration
Ci   nRLc  :dimensions dGdirr, dGdk, or dGdkr, depending on 10s digit opt
Ci          :  opt1 = 0 => lpdim = nbas
Ci          :  opt1 = 1 => lpdim = nnrl(1) ~ nbas*nl (see function nnrl)
Ci          :  opt1 = 2 => lpdim = nnrl(0) ~ nbas*nlm (see function nnrl)
Ci   nRLc  :dimensions dGdirr, dGdk, dGdkr:
Ci          0  dimension nbas
Ci          1  dimension nbas*nl
Ci          2  dimension nbas*nlm
Co Outputs
Co   Which of the following is generated depends on 10s digit opt
Co   dGdirr :change in the diagonal GF from perturbation at site ib
Co          is accumulated for this energy at the irreducible k-points
Co   dGdk   :change in the diagonal GF from perturbation at site ib
Co          is accumulated for this energy at the all k-points
Co   dGdkr  :Real or imaginary part of dGdk
Co          is accumulated for this energy at the all k-points
Cl Local variables
Cl   lmrci :compound (R,l) index for site ib, (contracted over mi)
Cl   lmrcj :compound (R,l) index for site jb, (contracted over mj)
Cl   is1   :effective leading spin index to gii,gjj,dpd,dpdd
Cl         := 1 in collinear case, loops over spins 1..2 in noncollinear
Cl   is2   :effective second spin index to gii,gjj
Cl         := 1 in collinear case, loops over spins 1..2 in noncollinear
Cl   ks2   :isp in the collinear case, and is2 in the noncollinear
Cl
Cl         :Thus (is1,is2) = spin indices to array gji
Cl         :     (is1,ks2) = spin indices to array dgd
Cr Remarks
Cu Updates
Cu   10 Jan 09 Added options
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer opt,ib,nbas,isp,nspc,nk1,nk2,nk3,k1,k2,k3,nkp,lidim,nlma,
     .  nRLc,iprmb(1),ipq(1),igstar(0:1)
      double complex dGdirr(nRLc,nspc,nRLc,nkp,isp),wz,
     .               dGdk(k1,k2,k3,nRLc,nspc,nRLc,isp),
     .               gji(k1,k2,k3,lidim,nspc,nlma,nspc)
      double precision dGdkr(k1,k2,k3,nRLc,nspc,nRLc,isp)
C ... Local parameters
      double precision pi
      double complex wt
      integer jb,offi,lj,jd,li,id,lmrci,lmrcj,nnrl,is1,is2,ks2,n0,nkap0,
     .  opt0,opt1,opt2
      parameter (n0=10,nkap0=4)
      integer norbi,norbj,ixx,iorb,jorb,nlmi,nlmj,
     .  ltabi(n0*nkap0),ktabi(n0*nkap0),offli(n0*nkap0),
     .  ltabj(n0*nkap0),ktabj(n0*nkap0),offlj(n0*nkap0)
      integer iq,jq,i1,i2,i3,ig

      if (nspc == 2) call rx('gfp0f3 : not ready for nspc=2')
      opt0 = mod(opt,10)
      opt1 = mod(opt/10,10)
      opt2 = mod(opt/100,10)
      call sanrg(.true.,opt1,0,2,'gfp0f3 (bug):','opt1')
      call sanrg(.true.,opt2,0,1,'gfp0f3 (bug):','opt2')

      do  10  is1 = 1, nspc
      do  10  is2 = 1, nspc
      ks2 = max(is2,isp)
      if (opt2 == 1) ks2 = 1

C ... wt is -1/pi * energy weight
      pi = 4*datan(1d0)
      wt = wz/pi
      if (opt0 == 3) wt = wt * (0d0,1d0)

C ... Find starting lmrci = offset to 1st column in dGdirr to store
      lmrci = nnrl(2-opt1,1,ib-1,iprmb,lidim)

C --- For each orbital in site ib, contract over m_i gji ---
C     Loop over row index with contiguous l,i,h blocks in dnfold order.
      call orbl(ib,0,lidim,iprmb,norbi,ltabi,ktabi,ixx,offli,nlmi)
C     opt1 = 0 : Increment lmrci for each new ib:
      if (opt1 == 0) lmrci = lmrci+1
      offi = 0
      do  iorb = 1, norbi
      li = ltabi(iorb)
      nlmi = 2*li + 1
C     opt1 = 1 : Increment lmrci for each new l
      if (opt1 == 1) lmrci = lmrci+1

      do  id = offi+1, offi+nlmi

C       opt1 = 2 : Increment lmrci for each new m
        if (opt1 == 2) lmrci = lmrci+1

C   --- For each orbital jd, store and contract over m_j,m_i gji ---
        do  jb = 1, nbas
C ...   Find starting lmrcj = offset to 1st row in dGdirr to store
        lmrcj = nnrl(2-opt1,1,jb-1,iprmb,lidim)
        if (nnrl(2-opt1,1,jb,iprmb,lidim) == lmrcj) cycle

C       Uses norbj,ntorbj
        call orbl(jb,0,lidim,iprmb,norbj,ltabj,ktabj,ixx,offlj,nlmj)
C       opt1 = 0 : Increment lmrcj for each new jb:
        if (opt1 == 0) lmrcj = lmrcj+1
        do  jorb = 1, norbj
        lj = ltabj(jorb)
        nlmj = 2*lj + 1
C       opt1 = 1 : Increment lmrcj for each new l
        if (opt1 == 1) lmrcj = lmrcj+1
        do  jd = offlj(jorb)+1, offlj(jorb)+nlmj

C         opt1 = 2 : Increment lmrcj for each new m
          if (opt1 == 2) lmrcj = lmrcj+1

C     ... Poke wt*gji for the irreducible qp
          if (opt0 == 0) then

C         Loop over all qp, selecting the irreducible ones
          jq = 0
          do  i3 = 1, nk3
          do  i2 = 1, nk2
          do  i1 = 1, nk1

            jq = jq+1
C           Irreducible point corresponds to unit rotation
C           Skip this qp unless it is an irreducible one
            ig = igstar(jq)
            if (ig /= 1) cycle

            iq = ipq(jq)
C           print 334, 'jq',jq,iq,i1,i2,i3
C 334       format(a,1x,6i5)

            dGdirr(lmrcj,is1,lmrci,iq,ks2) =
     .      dGdirr(lmrcj,is1,lmrci,iq,ks2) + wt*
     .                                       gji(i1,i2,i3,jd,is1,id,is2)

          enddo
          enddo
          enddo

          elseif (opt0 == 1) then

          do  i3 = 1, nk3
          do  i2 = 1, nk2
          do  i1 = 1, nk1

            dGdk(i1,i2,i3,lmrcj,is1,lmrci,ks2) =
     .      dGdk(i1,i2,i3,lmrcj,is1,lmrci,ks2) + wt*
     .                                       gji(i1,i2,i3,jd,is1,id,is2)

          enddo
          enddo
          enddo

          elseif (opt0 == 2 .or. opt0 == 3) then

          do  i3 = 1, nk3
          do  i2 = 1, nk2
          do  i1 = 1, nk1

            dGdkr(i1,i2,i3,lmrcj,is1,lmrci,ks2) =
     .      dGdkr(i1,i2,i3,lmrcj,is1,lmrci,ks2) + wt*
     .                                       gji(i1,i2,i3,jd,is1,id,is2)


          enddo
          enddo
          enddo

          endif

        enddo
      enddo
      enddo
      enddo
      offi = offi + nlmi
      enddo
   10 continue

      end
