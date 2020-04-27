      subroutine gffbz(s_ctrl,nbas,offH,iprmb,ifi,ib1,ib2,mode,
     .  pos,nkp,qp,zp,nk1,nk2,nk3,k1,k2,k3,ipq,plat,istab,g,ag,igstar,
     .  ifac,qb,offr,offc,ldh,gr,gc,ghh)
C- Generates gf connecting one site in full BZ from from file g(irr qp)
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  lham lncol nl
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   nbas  :size of basis
Ci   offH  :Offsets to hamiltonian matrix (makidx.f)
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   ifi   :file handle
Ci   ib1    :site for which to make gf
Ci   mode  :1s digit treats how gf is saved
Ci          0 read gf for this energy is from file ifi
Ci          1 same as 0, but rewind file and read header
Ci         10s digit concerns orbitals belonging to ib.
Ci            (row dimension in gr, column dimension in gc.)
Ci            Orbitals are ordered into l,i,h,n blocks
Ci            (lower, intermediate, higher, neglected; see makidx.f)
Ci          0 return only orbitals in l- block
Ci          1 same as 0
Ci          2 return only orbitals in i- block (not implemented)
Ci          3 return only orbitals in h- block (not implemented)
Ci          4 return orbitals l+i blocks
Ci          5 return orbitals l+i+h blocks
Ci        100s digit concerns dimension belonging to whole basis
Ci            (column dimension in gr, row dimension in gc.)
Ci            Orbitals are ordered into l,i,h,n blocks
Ci            (lower, intermediate, higher, neglected; see makidx.f)
Ci          0 return only orbitals in l- block
Ci          1 same as 0
Ci          2 return only orbitals in i- block (not implemented)
Ci          3 return only orbitals in h- block (not implemented)
Ci          4 return orbitals l+i blocks
Ci          5 return orbitals l+i+h blocks
Ci   pos   :basis vectors (input)
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   qp    :irreducible k-points
Ci   zp    :complex energy
Ci   nk1,nk2,nk3:  no. divisions for the 3 recip. latt. vecs
Ci   k1,k2,k3:  leading dimensions of gr,gc
Ci   ipq   :ipq(i1,i2,i3) points to the irreducible qp into which
Ci          mesh point (i1,i2,i3) is mapped (bzmesh.f)
Ci   plat  :primitive lattice vectors, in units of alat
Ci   istab :table of site permutations for each group op (mksym.f,symtbl.f)
Ci   g     :point group operations
Ci   ag    :translation part of space group
Ci   igstar:contains info needed to map an irreducible qp to its
Ci          original point (bzmesh.f)
Ci   ifac  :used with qb to recover a wave vector from three indices
Ci          specifying position in a uniform mesh (bzmsh0.f)
Ci   qb    :vectors of a microcell in the Brillouin zone
Ci   offr  :offset in gr add for e.g. second spin
Ci   offc  :offset in gc add for e.g. second spin
Ci   ldh   :leading dimension of ghh.  Not used unless ghh computed.
Co Outputs
Co   gr    :GF g_i,* for orbitals i in site ib = ib1..ib2; see Remarks
Co   gc    :GF g_*,i for orbitals i in site ib = ib1..ib2; see Remarks
Cb Bugs
Cb   Does orbital rotations only.
Cr Remarks
Cr   10s digit mode fixes the kind and number of orbitals nlma to be
Cr   included as part in site ib (whether l-block, i-block, l+i block,
Cr   etc), and sets nlma = size of that block.  Each site ib has its
Cr   own nlma.
Cr
Cr   10s digit mode fixes the kind and number of orbitals ldimb in the
Cr   lattice as a whole that connect to site ib.  ldimb is not site
Cr   dependent.
Cr
Cr   Let ilm1 is 1 + the cumulative sum of nlma for sites ib1...ib-1,
Cr   and ilm2 = ilm1-1+nlma.  Supressing the k-point indices,
Cr   gr(offr+ilm1..offr+ilm2,1..ldimb) and
Cr   gc(1..ldimb,offc+ilm1..offc+ilm2) hold the GF for site ib.
Cu Updates
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu  15 Mar 00  extended to include iwaves, stubs for h-waves
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer mode,nbas,nkp,nk1,nk2,nk3,k1,k2,k3,istab(nbas,1),ipq(*),
     .  igstar(0:1),ifac(3),offH(n0H,nkap0,1),iprmb(1),ib1,ib2,ifi,offr,
     .  offc,ldh
      double precision plat(3,3),qp(3,nkp),ghh(*)
      double precision pos(3,*),g(3,3,*),ag(3,*),qb(3,3),zp(2)
      double precision gr(2,nk1,nk2,nk3,*),gc(2,nk1,nk2,nk3,*)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
C ... Dynamically allocated local arrays
      complex(8), allocatable :: gii(:)
      complex(8), allocatable :: grk(:)
C ... Local parameters
      logical bittst
      integer ibla,iq,ldima,ldimb,lham,lncol,mode0,nl,nlma,nlmb,iblb,
     .  nspc,off2,optrot,iprint,lidim,ld11,ld21,ldrgn,ogi,ib,offri,offci
      double precision xx
C     integer jj1,jj2,jj3,k
C     double precision qk
C ... For file I/O of gf
      integer clp(9*2),iogfrs,kcplx
C     Bits for matching conditions (10-100s digit of lio)
      integer MCD,MZP
      parameter (MCD=32,MZP=64)
C     Bits for selecting which data to I/O, iogfrs
      integer ISH
      parameter (ISH=1)
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
C      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
C     .                    (jj2*ifac(2)-1)*qb(k,2) +
C     .                    (jj3*ifac(3)-1)*qb(k,3)

C --- Setup ---
C     call wkfast(.false.)
      call tcn('gffbz')
      lham = s_ctrl%lham
      lncol = s_ctrl%lncol
      nl = s_ctrl%nl
C     Spherical harmonics
      optrot = 0
      if (bittst(lham,256)) optrot = 20
      nspc = 1
      if (lncol /= 0) nspc = 2
      mode0 = mod(mode,10)
      ibla = mod(mode/10,10)
      if (ibla == 0) ibla = 1
      if (ibla == 5) ibla = 4
      if (ibla /= 1 .and. ibla /= 4) goto 999
      iblb = mod(mode/100,10)
      if (iblb == 0) iblb = 1
      if (iblb == 5) iblb = 4
      if (iblb /= 1 .and. iblb /= 4) goto 999
      if (ibla == 4 .or. iblb == 4) optrot = optrot + 4000
      ldima = offH(iblb,1,nbas+1)
      ldimb = ldima
      lidim = offH(4,1,nbas+1)
      kcplx = 0
C     call defdc(ogii,(lidim*nspc)**2)
      allocate(gii((lidim*nspc)**2))

      call iinit(clp,9*2)
      clp(3) = lidim*nspc
      clp(4) = lidim*nspc

C ... Rewind file, read header
      if (mod(mode0,2) /= 0) then
        rewind ifi
        call pshpr(1)
        if (iogfrs(2,0,0,' ',ifi,1,1,nbas,0,zp,qp,plat,xx,clp,
     .    xx,xx,0,0,0,xx) /= 0) call rx('gffbz failed to read header')
        call poppr
      endif

C --- For each irreducible qp, make G at each star ---
      nlma = 0
      do  ib = ib1, ib2
        nlma = max(nlma,offH(ibla,1,ib+1)-offH(ibla,1,ib))
      enddo
      off2 = nlma*lidim
      allocate(grk(4*off2))
      do  iq = 1, nkp

C   ... Read from disk g for this iq.  For now, no checking of zp
        call pshpr(iprint()-10)
        if (iogfrs(10000*ISH+0*MZP+10*MCD+6,1,0,' ',ifi,1,1,0,0,zp,
     .    qp(1,iq),xx,xx,clp,xx,xx,0,0,kcplx,gii) /= 0)
     .    call rxi('gffbz failed to read gf for qp no.',iq)
        call poppr

        call cplxdm(kcplx,lidim*nspc,lidim*nspc,ld11,ld21,ldrgn,ogi)
C       call yprm('gii',2,gii,ogi,lidim*nspc,lidim*nspc,lidim*nspc)

        offri = offr
        offci = offc
        do  ib = ib1, ib2
          nlma = offH(ibla,1,ib+1) - offH(ibla,1,ib)
          nlmb = nlma
          call pgfbz(nbas,nl,offH,iprmb,ib,optrot,pos,iq,nk1,nk2,nk3,
     .      k1,k2,k3,ipq,istab,g,ag,igstar,ifac,lidim,ldima,nlma,ldimb,
     .      nlmb,nspc,qb,gii,grk,grk,gr(1,1,1,1,1+offri),
     .      gc(1,1,1,1,1+offci))
          offri = offri + nlma*ldima*nspc**2
          offci = offci + ldimb*nlma*nspc**2
        enddo
      enddo
      deallocate(grk,gii)
      call tcx('gffbz')
      return

  999 call fexit(-1,111,'GFFBZ: mode %i is not implemented',mode)

      end
      subroutine pgfbz(nbas,nl,offH,iprmb,ib,opt,pos,iq,nk1,nk2,nk3,
     .  k1,k2,k3,ipq,istab,g,ag,igstar,ifac,lidim,ldima,nlma,ldimb,nlmb,
     .  nspc,qb,gii,grk,gck,gr,gc)
C- Kernel called by gffbz: gf connected to one site, qp in star of iq
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   offH  :Offsets to hamiltonian matrix (makidx.f)
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   ib    :site for which to generate connecting gf
Ci   opt   :specifies rotation options (see roth.f)
Ci         :1s digit is not used.  (calls roth using 1s digit = 4)
Ci         10s digit
Ci         :0 use real rotation mat for real harmonics
Ci         :1 use complex rotation matrix for cubic harmonics
Ci        100s digit distinguishes how complex arithmetic is handled
Ci         :  for now, must be zero
Ci       1000s digit specifies what part of g is to be extracted
Ci         :0 row dimension consists only of lower block
Ci         :1 same as zero
Ci         :2 not allowed yet
Ci         :3 not allowed yet
Ci         :4 row dimension consists lower+intermediate block
Ci   pos   :basis vectors
Ci   iq    :rotate gf to all q in star of iq
Ci   nk1,nk2,nk3:  no. divisions for the qp in 3 recip. latt. vecs
Ci   k1,k2,k3: leading dimensions of gr,gc
Ci   ipq   :ipq(i1,i2,i3) points to the irreducible qp into which
Ci          mesh point (i1,i2,i3) is mapped (bzmesh.f)
Ci   istab :table of site permutations for each group op (mksym.f,symtbl.f)
Ci   g     :point group operations
Ci   ag    :translation part of space group
Ci   igstar:contains group operation needed to rotated a k-point
Ci          to its irreducible one (bzmesh.f)
Ci   ifac  :used with qb to recover a wave vector from three indices
Ci          specifying position in a uniform mesh (bzmsh0.f)
Ci   lidim :number of lower+intermediate orbitals, defining dim. of gii
Ci   ldima :dimensions gc; also the number of rows in gc to fill.
Ci         :usually dimension of lower (or l+i) block for crystal
Ci   nlma  :dimensions gr; also the number of rows in gr to fill
Ci         :dimension of lower (or l+i) block for site ib
Ci   ldimb :dimensions gr; also the number of columns in gr to fill
Ci         :usually dimension of lower (or l+i) block for crystal
Ci   nlmb  :dimensions gc; also the number of cols in gc to fill
Ci         :dimension of lower (or l+i) block for site ib
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   qb    :vectors of a microcell in the Brillouin zone
Ci   gii   :gf (or hamiltonian) for this iq
Ci   grk   :a complex work array, dimensioned (nlma,lidim)*4
Ci   gck   :identical to grk: it MUST point to the same address space
Co Outputs
Co   gr    :for those qp in star iq, gii(k,*) stored, k= orbital in ib
Co   gc    :for those qp in star iq, gii(*,k) stored, k= orbital in ib
Cr Remarks
Cu Updates
Cu   10 Jul 13 Replace f77 pointers with f90 ones
Cu  15 Mar 00  extended to include iwaves
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ldima,nlma,lidim,ldimb,nlmb,nspc,opt
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer nbas,nk1,nk2,nk3,k1,k2,k3,istab(nbas,1),
     .  ipq(1),igstar(0:1),ifac(3),offH(n0H,nkap0,nbas),iprmb(ldima),ib
      double precision gii(lidim,nspc,lidim,nspc)
      double precision pos(3,*),g(3,3,*),ag(3,*),qb(3,3)
      double complex gr(k1,k2,k3,nlma,nspc,ldimb,nspc)
      double complex gc(k1,k2,k3,ldima,nspc,nlmb,nspc)
      double precision grk(nlma,lidim,2),gck(lidim,nlma,2)
C ... Dynamically allocated local arrays
      complex(8), allocatable :: wk(:)
C ... Local parameters
      integer i,i1,i2,i3,ig,iq,iq1,is,j,jj1,jj2,jj3,js,k,nl,off2,kcplx,
     .  optrot,ld11,ld21,ldrgn,ogi
      double precision q1(3),qk
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
     .                    (jj2*ifac(2)-1)*qb(k,2) +
     .                    (jj3*ifac(3)-1)*qb(k,3)


      optrot = opt + 4 - mod(opt,10)
      kcplx = mod(opt/100,10)
C     offset to first entry in gck (appended to grk)
      off2 = nlma*lidim*2
      call cplxdm(kcplx,lidim*nspc,lidim*nspc,ld11,ld21,ldrgn,ogi)
      if (kcplx /= 0) call rx('not set up for kcplx')
      call tcn('pgfbz')

C ... Copy gf subblock for each qp in star of qp(iq)
      iq1 = 0
      do  i3 = 1, nk3
      do  i2 = 1, nk2
      do  i1 = 1, nk1

        iq1 = iq1+1
C   ... skip this qp unless it is related to iq
        if (ipq(iq1) /= iq) cycle

C   ... Make g by rotation of g(iq): symop relating iq1 to iq
        ig = igstar(iq1)
C   ... q into which h is rotated
        q1(1) = qk(1,i1,i2,i3)
        q1(2) = qk(2,i1,i2,i3)
        q1(3) = qk(3,i1,i2,i3)

        do  is = 1, nspc
        do  js = 1, nspc

C         When rotating the l+i block, roth is set up to return in grk
C         for l and i blocks in both row and column dimensions.
          if (nspc == 2) then
            allocate(wk(lidim*lidim))
            call ymscop(0,lidim,lidim,lidim*nspc,lidim,0,0,0,0,
     .        gii(1,is,1,js),lidim*lidim*nspc**2,wk,lidim*lidim)
            call roth(optrot,nl,nbas,pos,ib,offH,iprmb,istab(1,ig),
     .        g(1,1,ig),ag(1,ig),q1,lidim,lidim,grk,wk)
            deallocate(wk)
          else
            call roth(optrot,nl,nbas,pos,ib,offH,iprmb,istab(1,ig),
     .        g(1,1,ig),ag(1,ig),q1,lidim,lidim,grk,gii(1,is,1,js))
          endif

C         call yprm('grk',2,grk,nlma*lidim,nlma,nlma,lidim)
C         call yprm('gck',2,gck(1+off2,1,1),nlma*lidim,lidim,lidim,nlma)
          do  j  = 1, ldimb
          do  i  = 1, nlma
            gr(i1,i2,i3,i,is,j,js) = dcmplx(grk(i,j,1),grk(i,j,2))
          enddo
          enddo
          do  j  = 1, nlmb
          do  i  = 1, ldima
            gc(i1,i2,i3,i,is,j,js) =
     .        dcmplx(gck(i+off2,j,1),gck(i+off2,j,2))
          enddo
          enddo
      enddo
      enddo

      enddo
      enddo
      enddo

      call tcx('pgfbz')

      end
