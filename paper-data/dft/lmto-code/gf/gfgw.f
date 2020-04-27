      subroutine gfgw(aintra,s_ctrl,s_ham,s_pot,s_lat,s_bz,llshft,
     .  offH,iprmb,nbas,lpdim,pos,nkp,qp,nk1,nk2,nk3,ipq,igstar,
     .  ifac,qb,nblk,wscrk,wscrl)
C- Sigma by gf-fft
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  lham lncol nl nspin lgen3
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  hamfbz
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  ldham hord
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  plat nsgrp
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:istab symgr ag
Cio    Passed to:  *
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  semsh
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   llshft:shift of the uniform BZ mesh
Ci   offH  :Offsets to hamiltonian matrix (makidx.f)
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   nbas  :size of basis
Ci   pos   :basis vectors (input)
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   qp    :irreducible k-points
Ci   nk1,nk2,nk3:  no. divisions for the 3 recip. latt. vecs
Ci   ipq   :ipq(i1,i2,i3) points to the irreducible qp into which
Ci          mesh point (i1,i2,i3) is mapped (bzmesh.f)
Ci   igstar:contains info needed to map an irreducible qp to its
Ci          original point (bzmesh.f)
Ci   ifac  :used with qb to recover a wave vector from three indices
Ci          specifying position in a uniform mesh (bzmsh0.f)
Ci   qb    :vectors of a microcell in the Brillouin zone(BZ)
Ci   wscrk :screened static Coulomb matrix defined at irreducible q-points
Ci   wscrl :local part of the screened static Coulomb matrix
Ci   nblk   :handles the size of the slice of sigma.
Co Outputs
Co  w(osig) :sigma saved on disk in file sgm.
Co          Because sigma is calculated in real space,
Co          it enables to treat slices of it rather than full sigma; this is
Co          needed since one has to have the object in the full BZ and it saves
Co          memory.
Co
Co
Co
Cr Remarks
Cr 10 Oct 04:(T. Sandu) Remade for spin-polarized case.
Cr Sept. '04: Remade to handle first the extension of wscrk in fullBZ and
Cr.then rolled from RL1,RL2 indices to RLm1,RLm2 indices(T. Sandu).
Cr
Cr Made first April '04 (T. Sandu): it expands integrated GF to full BZ,
Cr and wscrk from indices RL1,Rl2 to RLm1,RLm2 and then to full BZ.
Cu Updates
Cu   10 Jul 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed variables
      integer nbas,nkp,nk1,nk2,nk3,ipq(*),ifac(3),igstar(0:*),nblk
      integer nkap0,n0H
      integer lpdim
      parameter (nkap0=4,n0H=5)
      integer offH(n0H,nkap0,nbas+1),iprmb(*)
      double precision qp(3,nkp)
      double precision pos(3,*),qb(3,3)
      logical llshft(3),aintra
      double precision wscrk(lpdim,lpdim,2,nkp),wscrl(lpdim,lpdim,2)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_lat)::   s_lat
      type(str_bz)::    s_bz
C ... Dynamically allocated local arrays
      complex(8), allocatable :: gij(:)
      complex(8), allocatable :: gji(:)
      complex(8), allocatable :: wsij(:)
      complex(8), allocatable :: ghh(:)
      complex(8), allocatable :: gfsigm(:)
      complex(8), allocatable :: sig(:)
C ... Local parameters
      integer hord,ib,ib1,ib2,isp,k1,k2,k3,ldham(16),ldim,ldimx,lgen3,
     .  lham,lhdim,lhdimx,lidim,lidimx,lio,lncol,nl,nsp,nspc,
     .  offc,offr,iblk,nsgrp,nlma,nlmaa,i
      double precision plat(3,3),xx
      equivalence (ldim,ldham(1)),(ldimx,ldham(5)),(nspc,ldham(4))
      equivalence (lidim,ldham(2)),(lidimx,ldham(6))
      equivalence (lhdim,ldham(3)),(lhdimx,ldham(7))
C ... For file I/O of gf
      integer clp(9,2),fopnx,ifi,iogfrs
      integer MNC,MNB,MNT,MLB,MPL,MZP,MCD,ldh,mxorb,nglob
      parameter (MNC=1,MNB=2,MNT=4,MLB=8,MPL=16,MCD=32,MZP=64)
C ... For file I/O of gfwsq
      integer fisgm,iq
      character*8 fisnam
      double complex zp
      double precision semsh(10)
      integer fopna,nlmx

      call tcn('gfgw')

      lham = s_ctrl%lham
      lncol = s_ctrl%lncol
      nl = s_ctrl%nl
      nsp = s_ctrl%nspin
      lgen3 = s_ctrl%lgen3
      plat = s_lat%plat
      nsgrp = s_lat%nsgrp
      ldham = s_ham%ldham
      hord = s_ham%hord
      semsh = s_bz%semsh
      if (nspc == 2) call rx('gfgw: not ready for non-colinear')
c     kcplx = 0
C     nkfbz = nk1*nk2*nk3
      mxorb = nglob('mxorb')
      call fftz30(nk1,nk2,nk3,k1,k2,k3)

C     List of sites
      ib1 = 1
      ib2 = nbas

      zp = dcmplx(0d0,0d0)
      fisnam = 'sgm'
      fisgm = fopna(fisnam,-1,4)
      rewind fisgm

C --- For each site and spin do: ---
      do  isp = 1, nsp
C     Read int-of-GF, expand into full BZ; calculate sigma
      iblk = 0
      nlmx = 0
      offr = 0
      offc = 0
      do  ib = ib1, ib2

C --- int-of-gf for all qp in subblock connected to site ib ---
C ... Get gf for subblock connected to site ib
      nlma  = offH(4,1,ib+1) - offH(4,1,ib)
      lidim = offH(4,1,nbas+1)
      ldh   = offH(3,1,nbas+1) - offH(3,1,ib)
C ... If iblk<= ib, increment iblk and get gf for subblock ib,iblk
      if (iblk < ib) then
        iblk = min(ib+nblk-1,ib2)
        nlmaa = offH(4,1,iblk+1) - offH(4,1,ib)
        i = k1*k2*k3*nlmaa*lidim*nspc*nspc
        allocate(gij(i))
        allocate(gji(i))
        allocate(wsij(i)); call dpzero(wsij,2*i)
        allocate(ghh(ldh*mxorb))
        i = nlma*lidim*nspc*nspc
        allocate(gfsigm(nkp*i)); call dpzero(gfsigm,2*nkp*i)
        allocate(sig(i))

C       Rewind file containing int-gf, read header
        call iinit(clp,9*2)
        clp(3,1) = lidim*nspc
        clp(4,1) = lidim*nspc
        clp(5,1) = lidim*nspc

        if (isp == 1) then
          ifi = fopnx('gfdm1',100,16+8+4+0,-1)
        else
          ifi = fopnx('gfdm2',100,16+8+4+0,-1)
        endif
        rewind ifi
        call pshpr(1)
        lio = 2 + 10*(MNC*1+MNB*1+MNT*0+MLB*1+MPL*0+MCD*1+MZP*0)
        if (iogfrs(lio,0,0,' ',ifi,1,1,nbas,xx,zp,qp,plat,xx,clp,xx,xx,
     .    0,0,0,xx) /= 0) call rx('gfgw: failed to read file header')
        call poppr

C  ... Expand int-of-gf into full BZ
        call hamfbz(s_ctrl,s_ham,nbas,ifi,ib,iblk,440,pos,
     .    nkp,qp,zp,nk1,nk2,nk3,k1,k2,k3,ipq,plat,s_lat%istab,
     .    s_lat%symgr,s_lat%ag,igstar,ifac,qb,offr,offc,ldh,gij,gji,ghh)
c        call gfscriek(nlma,lidim,nsp,k1,k2,k3,2,5,1,gij)
      endif

C ... Expand wscrk into full BZ and unroll
      call gfwscrbz(aintra,plat,nl,nbas,offH,iprmb,nk1,nk2,nk3,nkp,
     .  llshft,ipq,igstar,lidim,lpdim,nlmx,nlma,pos,s_lat%istab,
     .  s_lat%symgr,s_lat%ag,wscrk,wscrl,wsij)

      call sigmagf(nk1,nk2,nk3,k1,k2,k3,nlma,lidim,gij,wsij)
      nlmx = nlmx + nlma

      call gfgetqir(nk1,nk2,nk3,k1,k2,k3,nkp,ipq,igstar,
     .  nspc,lidim,nlma,gij,gfsigm)
c        print *, 'after reducing sigma to irred.BZ'

C  ... write sigma on disk
      do  iq =1, nkp
c        print *, 'iq=',iq
c        call gfscrie(nlma,lidim,nsp,nkp,iq,gfsigm)
        call gfgwk2gw(nlma,lidim,nkp,iq,gfsigm,sig)
c        call yprm('bf-dmp-sig',2,sig,nlma*lidim,nlma,nlma,lidim)
        call dpdump(sig,nlma*lidim*nspc*2,-fisgm)
      enddo

      if (ib == iblk) deallocate(gij,gji,wsij,ghh,gfsigm,sig)
      enddo
      call fclose(ifi)
      enddo
      call fclose(fisgm)
      if (allocated(gij)) then
        deallocate(gij,gji,wsij,ghh,gfsigm,sig)
      endif
      call tcx('gfgw')
      end

      subroutine sigmagf(nk1,nk2,nk3,k1,k2,k3,nlmi,ldimi,gint,wsig)
C- SIGMA calculated with GF by  FFT
C ----------------------------------------------------------------------
Ci Inputs
Ci   nk1,nk2,nk3 : number of k-points along the three
Ci   k1,k2,k3    : leading dimensions of gint
Ci   nlmi  : number of orbitals belonging to site
Ci   ldimi : dimensions arrays gij,gji
Ci   ndim  : number of orbitals connected to site
Ci   gint,wsig : unperturbed GF integrated over energy
Ci         : connecting orbitals j at the field
Ci         : point to orbitals i; wsig static screened Coulomb matrix
Ci         : Both gint and wsig are in reciprocal space, full BZ.
Co Outputs
Co   gint    : is overwritten with convolution of wsig and gint
Co         : wsig is returned in reciprocal space.
Cr Remarks
Cr Made first April 04 (T.Sandu)
Cr 22 Sept 04 Changed wsig(k) to wsig(-k) (T. Sandu)
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nk1,nk2,nk3,k1,k2,k3,nlmi,ldimi
      double complex gint(k1,k2,k3,nlmi,ldimi),
     .               wsig(k1,k2,k3,nlmi,ldimi)


C ... Local parameters
      integer ilm,jd,ik1,ik2,ik3
      call tcn('sigmagf')

C ... Overwrite gint  with FT of convolution,  gint~ wsig ~
      do  20  ilm = 1, ldimi
      do  20  jd = 1, nlmi
c      print *, 'ilm',ilm
C       Real-space wsig(ilm) for nlmi channels
        call fftz3(wsig(1,1,1,jd,ilm),nk1,nk2,nk3,k1,k2,k3,1,
     .         0,-1)

C ... Real-space gint, for nlmi channels

        call fftz3(gint(1,1,1,jd,ilm),nk1,nk2,nk3,k1,k2,k3,1,
     .         0,-1)
C......... This is the product of FT wsig~, g_ij~
        do 22 ik3 = 1, nk3
        do 22 ik2 = 1, nk2
        do 22 ik1 = 1, nk1

           gint(ik1,ik2,ik3,jd,ilm) =
     .     wsig(ik1,ik2,ik3,jd,ilm)*gint(ik1,ik2,ik3,jd,ilm)

   22   continue
C       convolution from inverse FT of (FT of convolution)

        call fftz3(gint(1,1,1,jd,ilm),nk1,nk2,nk3,k1,k2,k3,1,
     .  0,1)
   20 continue

      call tcx('sigmagf')
      end

      subroutine gfgetqir(nk1,nk2,nk3,k1,k2,k3,nkp,ipq,igstar,
     .  nspc,lhdim,nlma,wsigm,sigm)
C- Contract sigma to irreducible q's
C ----------------------------------------------------------------------
Ci Inputs
Ci   nspc   :2 if spin-up and spin-down channels are coupled; else 1.
Ci   lhdim :size of l+i+h hamiltonian block
Ci   nlma  :size of basis for
Ci   wsigm :sigma for full BZ
Co Outputs
Co   sigm  :sigma for irreducible q-points
Co
Cl Local variables
Cl   is1   :effective spin index to wsigm and sigm
Cr Remarks
Cr June 26 '04 deleted nsp index from wsigm
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nspc,nk1,nk2,nk3,k1,k2,k3,nkp,lhdim,nlma,
     .  ipq(1),igstar(0:*)
      double complex sigm(nlma,lhdim,nkp),
     .               wsigm(k1,k2,k3,nlma,lhdim)

C ... Local parameters

      integer jd,id
      integer iq,jq,i1,i2,i3,ig

      if (nspc == 2) call rx('gfgetqir : not ready for nspc=2')
c      print *, 'into getqir'

C     ... Poke wt*gji for the irreducible qp
C         Loop over all qp, selecting the irreducible ones
          jq = 0
          do  30  i3 = 1, nk3
          do  30  i2 = 1, nk2
          do  30  i1 = 1, nk1

            jq = jq+1
c            print *,'jq',jq
C           Irreducible point corresponds to unit rotation
C           Skip this qp unless it is an irreducible one
            ig = igstar(jq)
            if (ig /= 1) goto 30
c            print *, 'into getqir'
c            print *, 'iq', iq
c            print *, 'i1 i2 i3',i1,i2,i3
            iq = ipq(jq)
c            print *, 'iq after', iq
c            print *,'jq',jq
C           print 334, 'jq',jq,iq,i1,i2,i3
C 334       format(a,1x,6i5)

            do  20  jd = 1, lhdim
            do  20  id = 1, nlma
             sigm(id,jd,iq) =
     .       wsigm(i1,i2,i3,id,jd)
   20       continue
c            call zprm('q-ir-sigm',2,sigm(1,1,iq),nlma,nlma,lhdim)

   30     continue
      end

C      subroutine gfscrie(nlma,lidim,nsp,nkp,iq,P0k)
CC- printing for debugging sigma for 1 k-point
CC ----------------------------------------------------------------------
CCi Inputs
CCi   nlma,lidim  :dimensions of sigma
CCi   nsp   :2 for spin-polarized case, otherwise 1
CCi   nkp   :number of k-points for which sigma is available
CCi   iq    :current k-point
CCi   P0k    :sigma matrix for all k-points
CCo Outputs
CCo
CCl Local variables
CCr Remarks
CCr Made May 21 2004
CC ----------------------------------------------------------------------
CC     implicit none
CC ... Passed parameters
C      integer nlma,lidim,nsp,nkp,iq
Cc      double precision P0k(2,nlma,lidim,nkp,nsp)
C      double complex P0k(nlma,lidim,nkp)
Cc      double complex P0k(lidim,nlma,nkp)
CC ... Local parameters
C      integer i1,ii,ij
Ccc      double precision P01(nlma,lidim,nsp,2)
Cc      double precision P01(nlma,lidim,2)
C
C      call zprm('gfgwk2gw-P0k',2,P0k(1,1,iq),nlma,nlma,lidim)
Cc      call zprm('gfgwk2gw-P0k',2,P0k(1,1,iq),lidim,lidim,nlma)
C
C
C
C      end

C      subroutine gfscriek(nume,nlma,lidim,nsp,k1,k2,k3,nk1,nk2,nk3,P0k)
CC- printing for debugging sigma for 1 k-point
CC ----------------------------------------------------------------------
CCi Inputs
CCi   nlma,lidim  :dimensions of sigma
CCi   nsp   :2 for spin-polarized case, otherwise 1
CCi   nkp   :number of k-points for which sigma is available
CCi   iq    :current k-point
CCi   P0k    :sigma matrix for all k-points
CCo Outputs
CCo
CCl Local variables
CCr Remarks
CCr Made May 21 2004
CC ----------------------------------------------------------------------
CC     implicit none
CC ... Passed parameters
C      integer nlma,lidim,nsp,nkp,iq,k1,k2,k3,nk1,nk2,nk3
Ccc      double precision P0k(2,nlma,lidim,nkp,nsp)
Cc      double complex P0k(k1,k2,k3,nlma,lidim,nsp)
C      double complex P0k(k1,k2,k3,nlma,lidim)
CC ... Local parameters
C      integer i1,ii,ij
C      character*5 nume
Ccc      double precision P01(nlma,lidim,nsp,2)
C      double complex P01(nlma,lidim)
C      do  10 ij = 1 , lidim
C      do  10 ii = 1, nlma
Cc 10      P01(ii,ij) = P0k(nk1,nk2,nk3,ii,ij,1)
C 10      P01(ii,ij) = P0k(nk1,nk2,nk3,ii,ij)
Cc      print *, 'in gfscriek'
C      call zprm(nume,2,P01,nlma,nlma,lidim)
C
C
C      end

cc      subroutine gfgwk2gw(nlma,lidim,nkp,iq,P0k,P01)
C- Transfers sigma for 1 k-point
C ----------------------------------------------------------------------
Ci Inputs
Ci   nlma,lidim  :dimensions of sigma
Ci   nkp   :number of k-points for which sigma is available
Ci   iq    :current k-point
Ci   P0    :sigma matrix for all k-points
Co Outputs
Co    P01 sigma for 1 k-point
Cl Local variables
Cr Remarks
Cr Made May 21 2004
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
cc      integer nlma,lidim,nkp,iq
cc      double precision P0k(2,nlma,lidim,nkp)
C     character*(20) fmt
C ... Local parameters
cc      integer ii,ij
cc      double precision P01(nlma,lidim,nsp,2)
cc      double precision P01(nlma,lidim,2)
C     fmt = '(9f15.7)'

cc      print *, 'gfgwkgw merge with p0k2p0'

cc      do  100  ij = 1, lidim
cc      do  100  ii = 1, nlma
cc      P01(ii,ij,1) = P0k(1,ii,ij,iq)
cc      P01(ii,ij,2) = P0k(2,ii,ij,iq)
cc  100 continue
cc      end
