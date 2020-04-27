      subroutine hamfbz(s_ctrl,s_ham,nbas,ifi,ib1,ib2,mode,pos,
     .  nkp,qp,zp,nk1,nk2,nk3,k1,k2,k3,ipq,plat,istab,g,ag,igstar,
     .  ifac,qb,offr,offc,ldh,gr,gc,ghh)
C- Generates gf-like object in full BZ from from file with gf(irr qp)
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
Ci   ib1   :site for which to make gf
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
Ci        1000s digit specifies whether to return a subblock of
Ci              g or the entire g:
Ci          0 generate subblocks gr and gc connecting sites
Ci            ib1..ib2.  Hamiltonian is rotated in subblocks,
Ci            and gr and gc are filled as a sequence of subblocks
Ci            (see description of gr,gc below).
Ci          1 generate entire g in full BZ.  In this case,
Ci            ib1 and ib2 are not used.  gc is not computed.
Ci            and gr must be the dimension of the hamiltonian.
Ci   pos   :basis vectors (file input)
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
Ci   igstar:contains group operation needed to rotated a k-point
Ci          to its irreducible one (bzmesh should be called with igstar(0)=-2)
Ci   ifac  :used with qb to recover a wave vector from three indices
Ci          specifying position in a uniform mesh (bzmsh0.f)
Ci   qb    :vectors of a microcell in the Brillouin zone
Ci   offr  :offset in gr add for e.g. second spin
Ci   offc  :offset in gc add for e.g. second spin
Ci   ldh   :leading dimension of ghh.  Not used unless ghh computed.
Co Outputs
Co   gr    :GF in the full BZ.  Its storage depends on 1000s digit mode.
Co         :1000s digit mode=0:
Co         :gr = g_i,* for orbitals i in site ib = ib1..ib2.
Co         :gr is stored as a sequence of ib2-ib1+1 independent arrays
Co         :of dimension gri(k1,k2,k3,nlma,nspc,ldimb,nspc), where
Co         :nlma depends on ib.  The first element in gri
Co         :corresponding to site ib is at gr(1,1,1,1+offr)
Co         :1000s digit mode=1:
Co         :gr = the full g_ij for all orbitals and k-points.
Co   gc    :GF g_*,i for orbitals i in site ib = ib1..ib2;
Co         :If 1000s digit is not zero, gc is not calculated,
co         :since all the information is equivalent to gr.
Co         :gc is stored as a sequence of ib2-ib1+1 independent arrays
Co         :of dimension gci(k1,k2,k3,ldima,nspc,nlmb,nspc), where
Co         :nlmb depends on ib.  The first element in gci
Co         :corresponding to site ib is at gc(1,1,1,1+offc)
Cl Local variables
Cl  offr   :offset in gr to array corresponding to site ib
Cl  offc   :offset in gc to array corresponding to site ib
Cb Bugs
Cb  To do: kcplx, pass gii through as array, rather than file read.
Cb  ghh is passed but not made
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
Cu Updates
Cu   09 Jul 15 Noncollinear case includes spinor rotation
Cu   10 Nov 11 Begin migration to f90 structures
Cu   20 Jun 02 Adapted from gffbz
C ----------------------------------------------------------------------
      use structures
      implicit none
C Passed variables
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer mode,nbas,nkp,nk1,nk2,nk3,k1,k2,k3,ib1,ib2,ifi,offr,offc,ldh
      integer istab(nbas,*),ipq(*),igstar(0:*),ifac(3)
      double precision plat(3,3),qp(3,nkp),ghh(*)
      double precision pos(3,*),g(3,3,*),ag(3,*),qb(3,3),zp(2)
      double precision gr(2,nk1,nk2,nk3,*),gc(2,nk1,nk2,nk3,*)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_ham)::   s_ham
C ... Dynamically allocated arrays
      complex(8), pointer :: gii(:,:,:,:),gwk(:,:,:,:),grk(:)
C Local variables
      logical bittst
      integer ibla,iq,ldima,ldimb,lham,lncol,mode0,mode3,nl,nlma,nlmb,
     .  iblb,nspc,off2,optrot,iprint,lidim,ld11,ld21,ldrgn,
     .  ogi,ib,offri,offci
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
      call tcn('hamfbz')
      lham = s_ctrl%lham
      lncol = s_ctrl%lncol
      nl = s_ctrl%nl
C     Spherical harmonics
      optrot = 0
      if (bittst(lham,256)) optrot = 20 ! Spherical harmonics
      nspc = 1
      if (lncol /= 0) nspc = 2
      mode0 = mod(mode,10)
      mode3 = mod(mode/1000,10)
      ibla = mod(mode/10,10)
      if (ibla == 0) ibla = 1
      if (ibla == 5) ibla = 4
      if (ibla /= 1 .and. ibla /= 4) goto 999
      iblb = mod(mode/100,10)
      if (iblb == 0) iblb = 1
      if (iblb == 5) iblb = 4
      if (iblb /= 1 .and. iblb /= 4) goto 999
      if (ibla == 4 .or. iblb == 4) optrot = optrot + 4000
      ldima = s_ham%offH(iblb,nbas+1)
      ldimb = ldima
      lidim = s_ham%offH(4,nbas+1)
      kcplx = 0

      allocate(gii(lidim,nspc,lidim,nspc))

      call iinit(clp,9*2)
      clp(3) = lidim*nspc
      clp(4) = lidim*nspc

C ... Rewind file, read header
      if (mod(mode0,2) /= 0) then
        rewind ifi
        call pshpr(1)
        if (iogfrs(2,0,0,' ',ifi,1,1,nbas,0,zp,qp,plat,xx,clp,xx,xx,0,
     .    0,0,xx) /= 0) call rx('hamfbz failed to read header')
        call poppr
      endif

C --- For each irreducible qp, make G in its star ---
      if (mode3 == 0) then
        nlma = 0
        do  5  ib = ib1, ib2
    5   nlma = max(nlma,s_ham%offH(ibla,ib+1)-s_ham%offH(ibla,ib))
        off2 = nlma*lidim
        allocate(grk(4*off2))
      endif

      do  10  iq = 1, nkp

C   ... Read from disk g for this iq.
        call pshpr(iprint()-10)
        if (iogfrs(10000*ISH+0*MZP+10*MCD+6,1,0,' ',ifi,1,1,0,0,zp,
     .    qp(1,iq),xx,xx,clp,xx,xx,0,0,kcplx,gii) /= 0)
     .    call rxi('hamfbz failed to read gf for qp no.',iq)
        call poppr

        call cplxdm(kcplx,lidim*nspc,lidim*nspc,ld11,ld21,ldrgn,ogi)
C       call yprm('gii',2,gii,ogi,lidim*nspc,lidim*nspc,lidim*nspc)

        if (mode3 == 0) then
          offri = offr
          offci = offc
          do  ib = ib1, ib2
            nlma = s_ham%offH(ibla,ib+1) - s_ham%offH(ibla,ib)
            nlmb = nlma
            call hamfb2(s_ham,nbas,nl,ib,optrot,pos,iq,nk1,nk2,nk3,
     .        k1,k2,k3,ipq,istab,g,ag,igstar,ifac,lidim,ldima,nlma,
     .        ldimb,nlmb,nspc,qb,gii,grk,grk,
     .        gr(1,1,1,1,1+offri),gc(1,1,1,1,1+offci))
            offri = offri + nlma*ldima*nspc**2
            offci = offci + ldimb*nlma*nspc**2
          enddo
        else
          if (lidim /= ldima)
     .      call rx('hamfbz: check this mode for downfolded case')
          allocate(gwk(lidim,nspc,lidim,nspc))
          call hamfb3(nbas,nl,s_ham%offH,s_ham%iprmb,optrot,pos,iq,nk1,nk2,nk3,
     .      k1,k2,k3,ipq,istab,g,ag,igstar,ifac,lidim,ldima,
     .      ldimb,nspc,qb,gii,gwk,gwk,gr)
        endif
   10 continue

C      call yprm('gr',3,gr,0,k1*k2*k3,k1*k2*k3,nlma*ldimb*nspc*nspc)
C      call yprm('gc',3,gc,0,k1*k2*k3,k1*k2*k3,nlmb*ldima*nspc*nspc)

      call tcx('hamfbz')
      return

  999 call fexit(-1,111,'HAMFBZ: mode %i is not implemented',mode)

      end
      subroutine hamfb2(s_ham,nbas,nl,ib,opt,pos,iq,nk1,nk2,nk3,
     .  k1,k2,k3,ipq,istab,g,ag,igstar,ifac,lidim,ldima,nlma,ldimb,nlmb,
     .  nspc,qb,gii,grk,gck,gr,gc)
C- Kernel called by hamfbz: gf connected to one site, qp in star of iq
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
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
Ci          to its irreducible one (bzmesh should be called with igstar(0)=-2)
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
Cu   09 Jul 15 Noncollinear case includes spinor rotation
Cu   15 Mar 00 extended to include iwaves
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer ldima,nlma,lidim,ldimb,nlmb,nspc,opt
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer nbas,nk1,nk2,nk3,k1,k2,k3,istab(nbas,1),
     .  ipq(1),igstar(0:*),ifac(3),ib
      double precision gii(lidim,nspc,lidim,nspc)
      double precision pos(3,*),g(3,3,*),ag(3,*),qb(3,3)
      double complex gr(k1,k2,k3,nlma,nspc,ldimb,nspc)
      double complex gc(k1,k2,k3,ldima,nspc,nlmb,nspc)
      double precision grk(nlma,lidim,4),gck(lidim,nlma,4)
C ... For structures
!      include 'structures.h'
      type(str_ham)::   s_ham
C ... Dynamically allocated local arrays
      real(8), pointer :: sll(:,:,:,:,:),grot(:,:,:,:,:)
      complex(8), pointer :: wk(:,:)
C ... Local parameters
C     logical :: debug=.false.
      logical lgrotu
      integer i,i1,i2,i3,ig,iq,iq1,is,j,jj1,jj2,jj3,js,k,nl,off2,kcplx,
     .  optrot,ld11,ld21,ldrgn,ogi,ldim,ofa1,ofa2,nlmal,nlmai
      procedure(integer) :: rotspd,rothnc

      double precision q1(3),qk,xv(1)
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
     .                    (jj2*ifac(2)-1)*qb(k,2) +
     .                    (jj3*ifac(3)-1)*qb(k,3)


      call tcn('hamfb2')

      optrot = opt + 4 - mod(opt,10)
      kcplx = mod(opt/100,10)
C     offset to first entry in gck (appended to grk)
      off2 = nlma*lidim*2
      call cplxdm(kcplx,lidim*nspc,lidim*nspc,ld11,ld21,ldrgn,ogi)
      if (kcplx /= 0) call rx('not set up for kcplx ne 0')
      ldim = s_ham%ldham(1)
      if (nspc == 2) then     ! To simplify, rotate the entire gii and extract a piece
C     if (nspc == 2 .or. .true.) then     ! Show that it works for all modes
        call pvrotd(optrot,s_ham%offH,nbas,ib,i1,i2,i,i,i,i,ofa1,ofa2,i,nlmal,nlmai,i,i,i,i,i)
        if (nlmal+nlmai /= nlma) call rx('bug in hamfb2')
        allocate(sll(lidim,nspc,lidim,nspc,2))   ! Use local copy of gii, to preserve original
        call dcopy((lidim*nspc)**2*2,gii,1,sll,1)
        grot => sll             ! for first symop (ig=1), rotated g = unrotated g
        lgrotu = .false.        ! lgrotu=T when spinors are rotated to global axis
        optrot = optrot - 4     ! rotate the entire sll and extract a piece
      else
        nullify(sll,grot)
      endif

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

C   --- Rotate g from qp(iq) to q1 ---
C   ... Noncollinear case
C       This branch works for both collinear and noncollinear cases, but is less efficient.
        if (nspc == 2) then
C       if (nspc == 2 .or. .true.) then  ! Uncomment only if uncomment similar block above
C     ... Rotate sll making all spinors parallel to z.
          if (ig > 1 .and. .not. lgrotu .and. nspc == 2) then
            call rotspn(30000+100*rotspd(0),1,nbas,nbas,nl,s_ham%iprmb,s_ham%eula,s_ham%neula,
     .        s_ham%qss(4),xv,xv,lidim,ldim,ldim,lidim,lidim,xv,sll)
            lgrotu = .true.
C           call yprm('gii(fixed z)',2,sll,lidim*nspc*2,lidim*nspc,lidim*nspc,lidim*nspc)
          endif

C         if (debug) call yprmi('gii i123=%3:1i iq=%i',[i1,i2,i3],iq,2,sll,ogi,lidim*nspc,lidim*nspc,lidim*nspc)

          if (ig /= 1) then

            if (associated(grot,sll)) allocate(grot(lidim,nspc,lidim,nspc,2))

C   ...     Rotate li+li block; put in grot
            if (rotspd(1) /= 0) call rx('hamfbz: rothnc not ready for modified rotspd convention')
            if (rothnc(optrot,nl,nspc,nbas,pos,i,s_ham%offH,s_ham%iprmb,
     .        istab(1,ig),g(1,1,ig),ag(1,ig),q1,lidim,lidim,gii,grot) < 0)
     .        call rx('problem with rothnc')

            if (lgrotu) then
C             call yprm('rotated gii(local z)',2,grot,lidim*nspc*2,lidim*nspc,lidim*nspc,lidim*nspc)
              call rotspn(30000+100*rotspd(1),1,nbas,nbas,nl,s_ham%iprmb,s_ham%eula,
     .          s_ham%neula,s_ham%qss(4),xv,xv,lidim,ldim,ldim,lidim,lidim,xv,grot)
C             call snot(ldim*nspc,i1,i2,i3,ig,iq,grot)
C             call yprm('gii(local z)',2,grot,lidim*nspc*2,lidim*nspc,lidim*nspc,lidim*nspc)
            endif
          endif

C         if (debug) call yprmi('grot %s,q=(%3;6d) ig=%i',q1,ig,2,grot,ogi,lidim*nspc,lidim*nspc,lidim*nspc)

          do  j  = 1, ldimb
            do  i  = 1, nlmal
              gr(i1,i2,i3,i,:,j,:) = dcmplx(grot(i+ofa1,:,j,:,1),grot(i+ofa1,:,j,:,2))
              gc(i1,i2,i3,j,:,i,:) = dcmplx(grot(j,:,i+ofa1,:,1),grot(j,:,i+ofa1,:,2))
            enddo
            do  i  = 1, nlmai
              gr(i1,i2,i3,nlmal+i,:,j,:) = dcmplx(grot(i+ofa2,:,j,:,1),grot(i+ofa2,:,j,:,2))
              gc(i1,i2,i3,j,:,nlmal+i,:) = dcmplx(grot(j,:,i+ofa2,:,1),grot(j,:,i+ofa2,:,2))
            enddo
          enddo

C          if (ib == 2) then
C          is = 1; js = 1
C          print *, ib,ig,i1,i2,i3
C          call zprm('gr',2,gr(i1,i2,i3,:,is,:,js),nlma,nlma,ldimb)
CC         call zprm('gc',2,gc(i1,i2,i3,:,:,:,:),ldima,ldima,nlmb)
C          endif

C       This branch works for both collinear and noncollinear cases,
C       but is not correct in the nc case because it doesn't rotate the spinor parts.
        else

C         if (debug) call yprmi('gii i123=%3i iq=%i',[i1,i2,i3],iq,2,gii,ogi,lidim*nspc,lidim*nspc,lidim*nspc)

          do  is = 1, nspc
          do  js = 1, nspc

C           When rotating the l+i block, roth is set up to return in grk
C           for l and i blocks in both row and column dimensions.
            if (nspc == 2) then
              allocate(wk(lidim,lidim))
              call ymscop(0,lidim,lidim,lidim*nspc,lidim,0,0,0,0,
     .          gii(1,is,1,js),lidim*lidim*nspc**2,wk,lidim*lidim)
              call roth(optrot,nl,nbas,pos,ib,s_ham%offH,s_ham%iprmb,istab(1,ig),
     .          g(1,1,ig),ag(1,ig),q1,lidim,lidim,grk,wk)
              deallocate(wk)
            else
              call roth(optrot,nl,nbas,pos,ib,s_ham%offH,s_ham%iprmb,istab(1,ig),
     .          g(1,1,ig),ag(1,ig),q1,lidim,lidim,grk,gii(1,is,1,js))
            endif

C           if (debug) call yprm('grk',2,grk,nlma*lidim,nlma,nlma,lidim)
C           if (debug) call yprm('gck',2,gck(1+off2,1,1),nlma*lidim,lidim,lidim,nlma)
            do  j  = 1, ldimb
              do  i  = 1, nlma
                gr(i1,i2,i3,i,is,j,js) = dcmplx(grk(i,j,1),grk(i,j,2))
              enddo
            enddo
            do  j  = 1, nlmb
              do  i  = 1, ldima
                gc(i1,i2,i3,i,is,j,js) = dcmplx(gck(i,j,3),gck(i,j,4))
              enddo
            enddo
          enddo
          enddo
      endif

      enddo
      enddo
      enddo

      if (associated(grot,sll)) then
      elseif (associated(grot)) then
        deallocate(grot)
      endif
      if (associated(sll)) deallocate(sll)

C      call yprm('gr',3,gr,0,k1*k2*k3,k1*k2*k3,nlma*ldimb*nspc*nspc)
C      call yprm('gc',3,gc,0,k1*k2*k3,k1*k2*k3,nlmb*ldima*nspc*nspc)

      call tcx('hamfb2')

      end
      subroutine hamfb3(nbas,nl,offH,iprmb,opt,pos,iq,nk1,nk2,nk3,
     .  k1,k2,k3,ipq,istab,g,ag,igstar,ifac,lidim,ldima,ldimb,nspc,
     .  qb,hq,hw,hw2,gfbz)
C- Kernel called by hamfbz: contr. to entire gf in full BZ from 1 qp.
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   offH  :Offsets to hamiltonian matrix (makidx.f)
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   opt   :specifies rotation options (see roth.f)
Ci         :1s digit is not used.  (calls roth using 1s digit = 0)
Ci         10s digit
Ci         :0 use real rotation mat for real harmonics
Ci         :1 use complex rotation matrix for cubic harmonics
Ci         :Add 4 if phase convention phi = q * [(g R_j + a) - R_i]
Ci         :  should be scaled by -1
Ci        100s digit distinguishes how complex arithmetic is handled
Ci             (kcplx mode)
Ci           0: h,h0 have real, imaginary separated
Ci              h = h(ldha,ldhb,2), with h(*,*,1..2) = real..imag
Ci           1: h,h0 are in complex*16 format (see Bugs)
Ci              h = h(2,ldha,ldhb), with s(1,*) = real, s(2,*) = imag
Ci           2: h,h0 have real, imaginary separated by columns
Ci              h = h(ldha,2,ldhb), with h(*,1..2,*) = real..imag
Ci       1000s digit specifies what part of g is to be extracted
Ci         :0 row dimension consists only of lower block
Ci         :1 same as zero
Ci         :2 not allowed yet
Ci         :3 not allowed yet
Ci         :4 row dimension consists lower+intermediate block
Ci   pos   :basis vectors
Ci   iq    :rotate gf to all q in star of iq
Ci   nk1,nk2,nk3:  no. divisions for the qp in 3 recip. latt. vecs
Ci   k1,k2,k3: leading dimensions of gfbz
Ci   ipq   :ipq(i1,i2,i3) points to the irreducible qp into which
Ci          mesh point (i1,i2,i3) is mapped (bzmesh.f)
Ci   istab :table of site permutations for each group op (mksym.f,symtbl.f)
Ci   g     :point group operations
Ci   ag    :translation part of space group
Ci   igstar:contains group operation needed to rotated a k-point
Ci          to its irreducible one (bzmesh should be called with igstar(0)=-2)
Ci   ifac  :used with qb to recover a wave vector from three indices
Ci          specifying position in a uniform mesh (bzmsh0.f)
Ci   lidim :number of lower+intermediate orbitals, defining dim. of hq
Ci   ldima :dimensions gfbz; also the number of rows in gfbz to fill.
Ci         :usually dimension of lower (or l+i) block for crystal
Ci   ldimb :dimensions gfbz; also the number of columns in gfbz to fill
Ci         :usually dimension of lower (or l+i) block for crystal
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   qb    :vectors of a microcell in the Brillouin zone
Ci   hq    :Green's function or hamiltonian for this iq
Ci   hw    :work array of same dimension as hq
Ci   hw2   :same work array as hw (used internally with diff. dim)
Co Outputs
Co   gfbz    :for those qp in star iq, hq stored
Cu Updates
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   20 Jun 02  First cut
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ldima,lidim,ldimb,nspc,opt
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer nbas,nk1,nk2,nk3,k1,k2,k3,istab(nbas,1),
     .  ipq(1),igstar(0:1),ifac(3),offH(n0H,nkap0,nbas),iprmb(ldima)
      double precision hq(lidim,nspc,lidim,nspc)
      double precision hw(lidim,nspc,lidim,nspc)
      double precision hw2(lidim,2,nspc,lidim,nspc)
      double precision pos(3,*),g(3,3,*),ag(3,*),qb(3,3)
      double complex gfbz(k1,k2,k3,ldima,nspc,ldimb,nspc)
C ... Dynamically allocated local arrays
      complex(8), allocatable :: wk(:)
C ... Local parameters
      integer i,i1,i2,i3,ig,iq,iq1,is,j,jj1,jj2,jj3,js,k,nl,kcplx,
     .  optrot,lidimx,offi
      double precision q1(3),qk
C Given (j1,j2,j3) of ipq, q_k(j1,j2,j3) =  sum_i (j_i*ifac(i)-1)*qb(k,i)
      qk(k,jj1,jj2,jj3) = (jj1*ifac(1)-1)*qb(k,1) +
     .                    (jj2*ifac(2)-1)*qb(k,2) +
     .                    (jj3*ifac(3)-1)*qb(k,3)


      call tcn('hamfb3')
      optrot = opt + 0 - mod(opt,10)
      lidimx = lidim*nspc
      kcplx = mod(opt/100,10)
C     rotate to real here, so roth doesn't do it repeatedly
      if (kcplx == 1) then
C       call zprm('hq before ztoyy',2,hq,lidimx,lidimx,lidimx)
        call ztoyy(hq,lidimx,lidimx,lidimx,lidimx,1,2)
C       roth works with kcplx=2 mode
        optrot = optrot + 200-100
      endif

C     Temporary memory for roth
      allocate(wk(lidimx**2))

C --- Copy gf for each qp in star of qp(iq) ---
      iq1 = 0
      do  12  i3 = 1, nk3
      do  12  i2 = 1, nk2
      do  12  i1 = 1, nk1

        iq1 = iq1+1
C   ... skip this qp unless it is related to iq
        if (ipq(iq1) /= iq) goto 12

C   ... Make g by rotation of g(iq): symop relating iq1 to iq
        ig = igstar(iq1)
C   ... q into which h is rotated
        q1(1) = qk(1,i1,i2,i3)
        q1(2) = qk(2,i1,i2,i3)
        q1(3) = qk(3,i1,i2,i3)

        if (nspc /= 1) call rx('hamfb3: check offsets for nspc')
C   ... Copy hq into hw to preserve gii; use hw
        call dcopy(2*lidimx**2,hq,1,hw,1)
        offi = lidimx**2
        if (kcplx /= 0) offi = lidimx

C   ... Rotate and copy
        do  14  is = 1, nspc
        do  14  js = 1, nspc
          if (kcplx == 0) then
          if (ig /= 1) then
            call roth(optrot,nl,nbas,pos,0,offH,iprmb,istab(1,ig),
     .      g(1,1,ig),ag(1,ig),q1,lidim,lidim,wk,hw(1,is,1,js))
          endif

          do  16  j  = 1, ldimb
          do  16  i  = 1, ldima
            gfbz(i1,i2,i3,i,is,j,js) =
     .                        dcmplx(hw(i,is,j,js),hw(i+offi,is,j,js))
   16     continue
          else
          if (ig /= 1) then
            call roth(optrot,nl,nbas,pos,0,offH,iprmb,istab(1,ig),
     .      g(1,1,ig),ag(1,ig),q1,lidim,lidim,wk,hw2(1,1,is,1,js))
          endif

          do  18  j  = 1, ldimb
          do  18  i  = 1, ldima
            gfbz(i1,i2,i3,i,is,j,js) =
     .                        dcmplx(hw2(i,1,is,j,js),hw2(i,2,is,j,js))
   18     continue
          endif

C         print 357,iq,i1,i2,i3,gfbz(i1,i2,i3,2,1,2,1)
C 357     format(' hamfb3: iq=',i4,' filling i1,i2,i3=',3i4,2f12.5)

   14   continue
   12 continue

      deallocate(wk)

C     Restore hq to complex*16 storage mode
      if (kcplx == 1) then
        call ztoyy(hq,lidimx,lidimx,lidimx,lidimx,2,1)
      endif

      call tcx('hamfb3')
C     call yprm('gf',3,gfbz,0,k1*k2*k3,k1*k2*k3,ldima*ldimb)

      end
