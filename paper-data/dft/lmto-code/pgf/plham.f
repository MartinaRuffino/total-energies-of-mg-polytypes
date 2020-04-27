      subroutine plham(mode,s_ham,s_pot,s_str,s_site,plat,isp,kcplx,ipl,pgplp,
     .  qp,ldim,ld0,hpl)
C- Structure constants connecting one principal layer to adjacent layers
C ----------------------------------------------------------------
Cio Structures
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  lgen3 ldham lncol lham neula offH
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb eula
Cio    Passed to:  *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:palp pfnc
Cio    Passed to:  *
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iax npr s
Cio    Passed to:  *
Ci Inputs:
Ci   mode  :1s digit
Ci         :0 make hamiltonian
Ci         :1 make noncollinear and/or relativistic hamiltonian if lncol is set.
Ci         :  hpl is returned in global spin quantization axis
Ci   plat  : primitive lattice vectors, in units of alat (input)
Ci   kcplx :distinguishes how complex arithmetic is handled
Ci          0: real, imaginary separated: hpl = hpl(ld0,ldim,1..2)
Ci          1: complex*16: hpl = hpl(ld0,ldim)
Ci          2: real, imaginary in columns : hpl = hpl(ld0,1..2,ldim)
Ci   ipl   :index to principal layer for layer programs
Ci          -1 => strux for left bulk PL; npl => strux for right bulk
Ci   pgplp :index and dimensioning information for each PL (pgfset.f)
Ci   qp    :k-point
Ci   ldim  :dimension of PL + adjacent left PL + adjacent right PL
Ci   ld0   :dimension of this PL
Co Outputs
Co   hpl   :2D Bloch summed PL strux, dimensioned(ld0,ldim)
Cl Local and global variables
Cl   nl    :(global maximum l) + 1
Cl   nbasp :number of atoms in the padded basis
Cl          nbasp = nbas + nbas(left bulk) + nbas(right bulk)
Cl   npl   :number of principal layers for layer programs (pgfset.f)
Cl   lhdimp:hamiltonian dimension, and leading dimension of pfun
Cr Remarks
Cr     hpl(1..ld0,1..ld0) S(PL,PL)
Cr     hpl(1..ld0,ld0+1..ld0+ldl) S(PL,PL-1)
Cr     hpl(1..ld0,ld0+ldl+1..ld0+ldr) S(PL,PL+1)
Cr   hpl must be dimensioned at least ld0*ldim
Cr
Cr       Pictorially,  hpl is dimensioned like the following:
Cr
Cr                ld0             ldl                ldr
Cr          ---------------------------------------------------
Cr          |            |                    |               |
Cr          |            |                    |               |
Cr       l  |            |                    |               |
Cr       d  |            |                    |               |
Cr       0  |            |                    |               |
Cr          |            |                    |               |
Cr          |            |                    |               |
Cr          ---------------------------------------------------
Cr                       ^                    ^
Cr                       offset ld0           offset ld0+ldl
Cr
Cr   sdpl, partitioning into lower and higher blocks not implemented.
Cu Updates
Cu   05 Jun 16 In SO case, plham uses s_site(:)%pfra to make P-S
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
Cu   08 Oct 04 Bug fix: use pass info about use of spherical harmonics to bloch
Cu   10 Jul 03 plham can generate P-S for fully rel case
Cu   10 Jan 02 plham reads pf from pot%palp instead of pot%pf
Cu   28 Apr 00 made noncollinear
C ----------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,isp,ldim,ld0,ipl,kcplx,pgplp(6,-1:*)
      double precision plat(3,3),qp(3)
      double complex hpl(*)
C     True dimensions of hpl for kcplx=0
C     double precision hpl(ld0,nspc,ldim,nspc,2)
C     hpl for kcplx=1
C     double complex hpl(ld0,nspc,ldim,nspc,2)
C     hpl for kcplx=2
C     double precision hpl(ld0,2,nspc,ldim,nspc)
C ... For structures
!      include 'structures.h'
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_str)::   s_str
      type(str_site)::  s_site(*)
C ... Dynamically allocated local arrays
      complex(8), allocatable :: pa(:),su(:)
C ... Local parameters
      logical lso,bittst
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer PRTG,ib1,ib2,idim,iprint,kl,lbloch,
     .  ldham(16),lds,ldl,ldr,lgen3,lham,lhdimp,lncol,lrel,
     .  mxorb,nbasp,neul,nglob,nl,nl2,npl,nspc,offi,offl,offp,offr
C     integer ibL1,ibL2,ibR1,ibR2,ncsw
      procedure(integer) :: bitand
      parameter (PRTG=81)
      double precision xx
      equivalence (lhdimp,ldham(3))

C --- Setup ---
      call tcn('plham')
      lgen3 = s_ham%lgen3
      ldham = s_ham%ldham
      lncol = s_ham%lncol
      lham = s_ham%lham

C     Make noncollinear hamiltonian only if 1s digit mode is set
      if (mod(mode,10) == 0) lncol = 0
      lbloch = 10*kcplx
C     Use spherical harmonics when making Bloch sum
      if (bittst(lham,256)) then
        lbloch = lbloch + 1000
      endif

      mxorb = nglob('mxorb')
      kl    = 0
      nbasp = nglob('nbasp')

C      call priprm(' ham offsets for system',1,nbasp,mxorb,ldim,ldim,
C     .  ldim,s_ham%iprmb)

      npl   = nglob('npl')
      nl    = nglob('nl')
      nl2   = nl*nl
      ldl = pgplp(4,max(ipl-1,-1))
      ldr = pgplp(4,min(ipl+1,npl))
      lds = ld0+ldl+ldr
      idim  = pgplp(5,ipl) - pgplp(4,ipl)
C     Sanity checks
      if (lgen3 /= 0) call rx('plham not ready for 3rd generation')
      if (idim /= 0) call rx('plham not ready for downfolding')
      if (ld0 > pgplp(4,ipl)) call rx('plham: improper argument ld0')
      if (ldim > lds) call rx('plham: improper argument ldim')
      lrel = mod(nglob('lrel'),10)
C ... Quick fix
      if (lncol == 0 .and. lrel == 2) lrel = 1
      nspc  = min(1+lncol,2)
      if (lrel == 2) nspc = 2

C --- 2D Bloch transform of real-space strux ---
      if (lrel == 2 .or. lncol /= 0) then
        allocate(su(ld0*ldim))
        call blchpl(lbloch,qp,nl,plat,s_ham%iprmb,ldham,ipl,npl,pgplp,
     .    s_str%iax,s_str%npr,s_str%s,nl2,ld0,0,ldim,0,su,xx,xx)
!       print *, '!!'; call dpzero(su,2*size(su))
      else
        call blchpl(lbloch,qp,nl,plat,s_ham%iprmb,ldham,ipl,npl,pgplp,
     .    s_str%iax,s_str%npr,s_str%s,nl2,ld0,0,ldim,0,hpl,xx,xx)
      endif

C --- ASA : Make (S-P) in S00 ---
C     if (lrel == 1) call zprm('palp',2,s_pot%palp,lhdimp,lhdimp,2)
C     if (lrel == 2) call zprm('palp',2,s_pot%palp,lhdimp,lhdimp,4)
      call gtibpl(ipl,npl,pgplp,ib1,ib2)

C     This branch is obsolete ... never accessed
      if (lrel == 2 .and. lncol == 0) then

        call rx('plham: obsolete branch')
C        neul = s_ham%neula
C
CC       Convert raw hamiltonian offsets offi,offL,offR into layer form
CC       Raw: offi,offL,offR = offset to ipl, ipl-1, ipl+1
CC       offi -> 0, offL -> ld0, offR -> ld0+ldl
C        call pghoff(ipl,npl,pgplp,s_ham%iprmb,mxorb,offL,offi,offR)
C        call gtibpl(max(ipl-1,-1),npl,pgplp,ibL1,ibL2)
C        call gtibpl(min(ipl+1,npl),npl,pgplp,ibR1,ibR2)
C
CC       These lines manipulate the orbital ordering array
CC       by excluding all sites not in PL, PL-1, PL+1.
CC       Later the original array is restored.
CC       Initially exclude all pairs by adding largest index
CC       but restore PL, PL-1, PL+1
CC       Converted: offi,offL,offR = 0, ld0, ld0+ldl
CC       BUG: doesn't work if ipl=-1 or ipl=npl
C        offi   = offi + lhdimp
C        offL   = offL - ld0 + lhdimp
C        offR   = offR - ld0 - ldl + lhdimp
C
CC       Initially exclude all atoms by shifting off array top
C        call pblch2(lhdimp,1,nbasp,nl2,s_ham%iprmb)
CC       Orbitals for ib1..ib2 start at offset 0
C        call pblch2(-offi,ib1,ib2,nl2,s_ham%iprmb)
CC       Orbitals for ibL1..ibL2 start at offset ld0
C        call pblch2(-offL,ibL1,ibL2,nl2,s_ham%iprmb)
CC       Orbitals for ibR1..ibR2 start at offset ld0 + ldl
C        call pblch2(-offR,ibR1,ibR2,nl2,s_ham%iprmb)
C
CC       Shift and copy pf(ib1..ib2) to pa
C        if (lrel == 1) then
CC         call defcc(opa,ld0*2)
C          allocate(pa(ld0*2))
CC         NB: this assumes P is ordered into lower+intermediate
CC         block, followed by higher blocks.
C          offp = s_ham%offH(1+nkap0*n0H*(ib1-1)+3,1)
C          call zmscop(0,ld0,2,lhdimp,ld0,offp,0,0,0,s_pot%palp,
C     .      pa)
C          call dscal(ld0*4,-1d0,pa,1)
CC         call zprm('pf',2,pa,ld0,ld0,2)
C        else
CC         call defcc(opa,ld0*2*2)
C          allocate(pa(ld0*2*2))
CC         Offset to idxsh array corresponding to 1st orbital of ib1
CC         NB: this assumes P is ordered into lower+intermediate
CC         block, followed by higher blocks.
C          offp = s_ham%offH(1+nkap0*n0H*(ib1-1)+3,1)
C          call zmscop(0,ld0,4,lhdimp,ld0,offp,0,0,0,s_pot%palp,pa)
C
C          call dscal(ld0*8,-1d0,pa,1)
CC         call zprm('pf',2,pa,ld0,ld0,4)
C        endif
C
CC       The following sequence rotates structure constants and adds
CC       P functions, that would be diagonal in the scalar-rel case.
CC       (Not the most efficient implementation, but it is simple.)
C
CC       Rotate the strux, then add pa
C        if (lrel == 1) ncsw = 1000 + 10*kcplx
C        if (lrel == 2) ncsw = 4000 + 10*kcplx
C        call rotspn(ncsw,1,nbasp,nbasp,nl,s_ham%iprmb,s_ham%eula,neul,xx,
C     .    pa,xx,ld0,ld0,lds,ld0,ldim,su,hpl)
C
CC       Reverse the rotation, so strux are diagonal in spin space
C        ncsw = 10100 + 10*kcplx
C        call rotspn(ncsw,1,nbasp,nbasp,nl,s_ham%iprmb,s_ham%eula,neul,xx,
C     .    pa,xx,ld0,ld0,lds,ld0,ldim,su,hpl)
C
CC       Restore the permutation table to its original state
C        call pblch2(-lhdimp,1,nbasp,nl2,s_ham%iprmb)
C        call pblch2(offi,ib1,ib2,nl2,s_ham%iprmb)
C        call pblch2(offL,ibL1,ibL2,nl2,s_ham%iprmb)
C        call pblch2(offR,ibR1,ibR2,nl2,s_ham%iprmb)

      elseif (lncol /= 0) then  ! Noncollinear branch
        neul = s_ham%neula
        lso = bittst(lncol,4)
C       if (lso .and. lrel /= 2) call rx('plham not set up for spin-orbit')

        offp = s_ham%offH(1+nkap0*n0H*(ib1-1)+3,1)
        if (lrel /= 2 .and. lso) then
C         if (bitand(lham,128) /= 0) call rx('SO improperly implemented in gamma repsn')
          call plhamso(0,s_ham,s_site,kcplx,ld0,ldim,ib1,ib2,ldham,mxorb,s_ham%iprmb,lhdimp,s_pot%pfr,su,su,hpl)
        elseif (lrel /= 2) then
          call plham2(1,kcplx,ld0,ldim,lhdimp,offp,s_pot%pfnc,su,su,hpl)
        else
          call plham2(2,kcplx,ld0,ldim,lhdimp,offp,s_pot%palp,su,su,hpl)
        endif

CC    .. Rotate -S by euler angles
CC       Convert raw hamiltonian offsets offi,offL,offR into layer form
CC       offi -> 0, offL -> ld0, offR -> ld0+ldl
C        call pghoff(ipl,npl,pgplp,s_ham%iprmb,mxorb,offL,offi,offR)
C        call gtibpl(max(ipl-1,-1),npl,pgplp,ibL1,ibL2)
C        call gtibpl(min(ipl+1,npl),npl,pgplp,ibR1,ibR2)
CC       Initially exclude all pairs by adding largest index
CC       but restore PL, PL-1, PL+1
C        offi   = offi + lhdimp
C        offL   = offL - ld0 + lhdimp
C        offR   = offR - ld0 - ldl + lhdimp
C        call pblch2(lhdimp,1,nbasp,nl2,s_ham%iprmb)
C        call pblch2(-offi,ib1,ib2,nl2,s_ham%iprmb)
C        call pblch2(-offL,ibL1,ibL2,nl2,s_ham%iprmb)
C        call pblch2(-offR,ibR1,ibR2,nl2,s_ham%iprmb)
CC       Shift and copy pf(ib1..ib2) to pa
C        call defcc(opa,ld0*2)
C        offp = w(ooffH+nkap0*n0H*(ib1-1)+3)
C        call zmscop(0,ld0,2,lhdimp,ld0,offp,0,0,0,s_pot%palp,pa)
C        call dscal(ld0*4,-1d0,pa,1)
CC       call zprm('pf',2,pa,ld0,ld0,2)
CC       Following sequence rotates P functions.  Not the most
CC       efficient, but simple implementation
CC       Rotate the strux, adding pa to the diagonal
C        ncsw = 1000 + 10*kcplx
C        call rx('plham check sign of rotation')
C        call rotspn(ncsw,1,nbasp,nbasp,nl,s_ham%iprmb,s_ham%eula,neul,xx,
C          pa,xx,ld0,ld0,lds,ld0,ldim,su,hpl)
CC       Reverse the rotation, so strux are diagonal in spin space
C        ncsw = 10100 + 10*kcplx
C        call rx('plham check sign of rotation')
C        call rotspn(ncsw,1,nbasp,nbasp,nl,s_ham%iprmb,s_ham%eula,neul,xx,
C          pa,xx,ld0,ld0,lds,ld0,ldim,su,hpl)
CC       Undo changes to permutation table
C        call pblch2(-lhdimp,1,nbasp,nl2,s_ham%iprmb)
C        call pblch2(offi,ib1,ib2,nl2,s_ham%iprmb)
C        call pblch2(offL,ibL1,ibL2,nl2,s_ham%iprmb)
C        call pblch2(offR,ibR1,ibR2,nl2,s_ham%iprmb)

      else                      ! Collinear branch

C       debugging ...
C        call gtibpl(max(ipl-1,-1),npl,pgplp,ibL1,ibL2)
C        call gtibpl(min(ipl+1,npl),npl,pgplp,ibR1,ibR2)
C        call priprm(' ham offsets for current PL',
C     .    ib1,ib2,mxorb,ldim,ldim,ldim,s_ham%iprmb)


C       Get offi = Hamiltonian offset for layer ipl
        call pghoff(ipl,npl,pgplp,s_ham%iprmb,mxorb,offL,offi,offR)
C       Shift offset for this layer so index starts at 1
C       Shifts for l block only (no i,h blocks)
        call pblch3(0,1,offi,ib1,ib2,mxorb,ldham,s_ham%iprmb)
C       call pblch2(-offi,ib1,ib2,mxorb,s_ham%iprmb)
        offp = s_ham%offH(1+nkap0*n0H*(ib1-1),1) + (isp-1)*lhdimp
C        call priprm(' modified ham offsets for current PL',
C     .    ib1,ib2,mxorb,ldim,ldim,ldim,s_ham%iprmb)

        call sblhm1(lbloch,nl2,ib1,ib2,s_ham%iprmb,s_pot%palp,offp,ldim,
     .    idim,ld0,idim,ldim,kl,hpl,xx,xx)
        call pblch3(11,1,offi,ib1,ib2,mxorb,ldham,s_ham%iprmb)
C       call pblch2(offi,ib1,ib2,mxorb,s_ham%iprmb)
C       call priprm(' restored ham offsets for current PL',
C    .    ib1,ib2,mxorb,ldim,ldim,ldim,s_ham%iprmb)

      endif

!     if (ipl == 6) then
      if (iprint() >= PRTG/1) then
        call yprmi('plham z-H for layer %i',ipl,0,kcplx+2,hpl,ld0*ldim*nspc**2,ld0*nspc,
     .    ld0*nspc,ldim*nspc)
      endif

      if (allocated(su)) deallocate(su)
      if (allocated(pa)) deallocate(pa)

C      call ztoyy(hpl,ld0*nspc,ldim*nspc,ld0*nspc,ldim*nspc,lbloch/10,0)

      call tcx('plham')
      end
      subroutine plham2(mode,kcplx,ld0,lds,ldp,offp,pfnc,sp1,sp2,snc)
C- Assemble noncollinear S-P for one PL from noncollinear P and collinear S
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 copy spl1 to 11 block, spl2 to 22 block
Ci         :1 Also add pfnc to 11,12,21,22 blocks.
Ci         :2 Add pfnc to 11,12,21,22 blocks (rel=2).
Ci   kcplx :distinguishes how complex arithmetic is handled
Ci          0: real, imaginary separated: hpl = hpl(ldim0,ldim,1..2)
Ci          1: complex*16: hpl = hpl(ldim0,ldim)
Ci          2: real, imaginary in columns : hpl = hpl(ldim0,1..2,ldim)
Ci   ld0   :dimension of current layer
Ci   lds   :second dimension of s (ld0+ldl+ldr)
Ci   ldp   :leading dimension of p
Ci   offp  :offset to pfnc for this PL
Ci   pfnc  :Potential functions backward rotated by Euler angles
Ci         :Thus P must be in the global spin quantization axis
Ci   sp1   :strux for this layer, and coupling to (left,right), spin1
Ci   sp2   :strux for this layer, and coupling to (left,right), spin2
Co Outputs
Co   snc   :noncollinear S-P for this PL
Co         :snc is returned in the global spin quantization axis
Cr Remarks
Cr   In the noncollinear case spl for layer i looks like
Cr   (Here j=layer i-1, k=layer i+1, and ldim = ld0+ldl+ldr)
Cr
Cr                                      ldim+  ldim+  ldim+
Cr                0      ld0   ld0+ldl    0    ld0   ld0+ldl
Cr                |      |       |        |       |      |
Cr       0  ->   \/     \/      \/       \/      \/     \/
Cr                ii++   ij++   ik++      ii+-    ij++   ik+-
Cr                ii-+   ij-+   ik-+      ii--    ij--   ik--
Cr
Cr  *Note: the (orbital,ipl,spin) order plham2 generates may be inconvenient.
Cr   Routine pgflu2 rearranges the elements to (orbital,spin,ipl) order, i.e.:
Cr                0     ld0     2*ld0  2*ld0+ldl    2*ld2  2*ld2+ldr   where ld2=ld0+ldl
Cr                |     |         |      |            |      |
Cr     0    ->    \/    \/        \/     \/           \/     \/
Cr                ii++  ii+-      ij++   ij+-         ik++   ik+-
Cr                ii-+  ii--      ij-+   ij--         ik-+   ik--
Cr
Cr Bugs
Cr   This routine, pokepf, and pfr2block should be merged into one routine
Cu Updates
Cu   18 Jun 04 (A Chantis) potential functions for relativistic case
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,kcplx,ld0,lds,ldp,offp,nspc
      parameter (nspc=2)
      double complex sp1(ld0,lds),sp2(ld0,lds)
      double complex pfnc(ldp,2,2),snc(ld0,nspc,lds,nspc)
C ... Local parameters
      integer i,j,i1,morder,im

      call sanrg(.true.,kcplx,1,1,'plham2:','kcplx')
      call lmorder(0,morder,[0],[0])
      im = -1; if (morder==1 .or. morder==2) im = 1

      call dpzero(snc,ld0*nspc*lds*nspc*2)
      call zmscop(0,ld0,lds,ld0,ld0*2,0,0,0,0,sp1,snc)
      call zmscop(0,ld0,lds,ld0,ld0*2,0,0,ld0,lds,sp2,snc)

C      call zprm('-S in plham2',2,snc,ld0*nspc,ld0*nspc,lds*nspc)
C      snc = 0

      if (mod(mode,10) == 0) then
      elseif (mod(mode,10) == 1) then
        do  i = 1, 2
        do  j = 1, 2
          call zaxpy(ld0,dcmplx(-1d0,0d0),pfnc(1+offp,i,j),1,snc(1,i,1,j),ld0*nspc+1)
        enddo
        enddo
      elseif (mod(mode,10) == 2) then
        do  i = 1, 2
          call zaxpy(ld0,dcmplx(-1d0,0d0),pfnc(1+offp,i,i),1,snc(1,i,1,i),ld0*nspc+1)
        enddo
        do  i1 = 1, ld0
          if (i1+im == 0 .or. i1+im > ld0) cycle
          snc(i1+im,1,i1,2) = snc(i1+im,1,i1,2) + pfnc(i1+im+offp,1,2)
          snc(i1,2,i1+im,1) = snc(i1,2,i1+im,1) + pfnc(i1+offp,2,1)
        enddo
      endif

C     call zprm('P-S in plham2',2,snc,ld0*nspc,ld0*nspc,lds*nspc)

      end
      subroutine plhamso(opt,s_ham,s_site,kcplx,ld0,lds,ib1,ib2,ldham,mxorb,iprmb,ldp,pfnc,sp1,sp2,snc)
C- Assemble noncollinear S-P for one PL from site-based P and collinear S
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  norb
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pfr
Cio    Passed to:  *
Ci Inputs
Ci   opt   :0 use matrix form of pf.
Ci         :1 use vector form of pf. CPA not allowed.
Ci   kcplx :distinguishes how complex arithmetic is handled
Ci          0: real, imaginary separated: snc = snc(ld0,2,lds,2,1..2)
Ci          1: complex*16: snc = snc(ld0,2,lds,2)
Ci          2: real, imaginary in columns : snc = snc(ld0,2,1..2,lds,2)
Ci   ld0   :dimension of current layer
Ci   lds   :second dimension of s (ld0+ldl+ldr)
Ci   ib1   :first site in current PL
Ci   ib2   :last site in current PL
Ci   ldham :dimension of lmto basis for all layers to be included in inversion
Ci   mxorb :leading dimension of iprmb, and the maximum number
Ci         :of orbitals connected with a single site.
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   sp1   :strux for this layer, and coupling to (left,right), spin1
Ci   sp2   :strux for this layer, and coupling to (left,right), spin2
Co Outputs
Co   snc   :noncollinear S-P for this PL, in global spin quantization axis
Cb Bugs
Cb  *Rotations to noncollinear configuration assume spherical harmonics
Cb   (not real harmonics)
Cb  *Follow mksopf convention for pfr: pfr ordered as s,p,d,...  no downfolding
Cr Remarks
Cr   In the noncollinear case spl for layer i looks like
Cr   (Here j=layer i-1, k=layer i+1, and ldim = ld0+ldl+ldr)
Cr
Cr                                      ldim+  ldim+  ldim+
Cr                0      ld0   ld0+ldl    0    ld0   ld0+ldl
Cr                |      |       |        |       |      |
Cr       0  ->   \/     \/      \/       \/      \/     \/
Cr                ii++   ij++   ik++      ii+-    ij++   ik+-
Cr                ii-+   ij-+   ik-+      ii--    ij--   ik--
Cr
Cr  *Note: the (orbital,ipl,spin) order plhamso generates may be inconvenient.
Cr   Routine pgflu2 rearranges the elements to (orbital,spin,ipl) order.
Cu Updates
Cu   17 Jun 18  New opt to allow reading of vector form of pf.
Cu   14 Feb 15  First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer opt,kcplx,ld0,lds,nspc,ib1,ib2,ldham,mxorb,iprmb(mxorb,*),ldp
      parameter (nspc=2)
      double complex sp1(ld0,lds),sp2(ld0,lds),pfnc(ldp,2,2)
      double complex snc(ld0,nspc,lds,nspc)
C ... For structures
!      include 'structures.h'
      type(str_ham)::   s_ham
      type(str_site)::  s_site(*)
C ... Dynamically allocated arrays
      complex(8), allocatable :: pfr(:,:,:,:)
C ... Local parameters
      integer i1,i2,norb,ib,icp,norbl,offib,offi,offj,nl,nglob,neul,lmr
      procedure(integer) :: rotspd
C#ifdefC DEBUG
C      logical :: debug = .false.
C#endif

      call sanrg(.true.,kcplx,1,1,'plhamso:','kcplx')
      nl = nglob('nl')
      neul = s_ham%neula
      call dpzero(snc,ld0*nspc*lds*nspc*2)
      call zmscop(0,ld0,lds,ld0,ld0*2,0,0,0,0,sp1,snc)
      call zmscop(0,ld0,lds,ld0,ld0*2,0,0,ld0,lds,sp2,snc)

      offib = 0                 ! offset to start of snc for site ib
      icp = 1                   ! No CPA for now
      do  ib = ib1, ib2
        norb = s_site(ib)%norb  ! (lmax+1)*2, as stored in mksopf
        allocate(pfr(norb,2,norb,2))
        if (opt == 0) then
          call zcopy(norb*norb*4,s_site(ib)%pfra(1,icp),1,pfr,1)
        else
          if (icp /= 1) call rx('plhamso: cannot use opt=1 for CPA')
          lmr = nl*nl*(ib-1)
          call pfmat2vec(2,nl,lmr,ldp,norb,iprmb,pfr,pfnc)
        endif
C#ifdefC DEBUG
C        if (debug) then
C          call yprmi('unrotated rotated PF for site %i',ib,0,3,pfr,0,norb*2,norb*2,norb*2)
C        endif
C#endif
        call rot_LandS(10+rotspd(0),s_ham%eula(ib,:),nl,norb,1,pfr) ! local -> global z if noncollinear
C#ifdefC DEBUG
C        if (debug) then
C          call zprm('rotated PF',2,pfr,norb*2,norb*2,norb*2)
C        endif
C#endif
        norbl = 0               ! Count number of orbitals in lower block site ib
        offi = offib
        do  i1 = 1, norb
          offi = offi+1
          if (iprmb(i1,ib) > ldham) cycle
          offj = offib
          do  i2 = 1, norb
            if (iprmb(i2,ib) > ldham) cycle
            offj = offj+1
            snc(offi,1:2,offj,1:2) = snc(offi,1:2,offj,1:2) - pfr(i1,1:2,i2,1:2)
C            do  i = 1, 2
C              do  j = 1, 2
C                snc(offi,i,offj,j) = snc(offi,i,offj,j) - pfr(i1,i,i2,j)
C              enddo
C            enddo
          enddo
          norbl = norbl + 1     ! Accumulate number of orbitals in lower block for this ib
        enddo
        offib = offib + norbl   ! Update offset, ready for ib+1
        deallocate(pfr)
      enddo

C#ifdefC DEBUG
C      if (debug) then
C      call zprm('plhamso S-P',2,snc,2*ld0,2*ld0,2*lds)
C      endif
C#endif

      end
      subroutine priprm(strn,ib1,ib2,mxorb,ldim,lidim,lihdim,iprmb)
C- Print out iprmb offsets
      implicit none
      integer ib1,ib2,mxorb,ldim,lidim,lihdim,iprmb(*)
      character strn*(*)
      character type*4
      integer ib0,offi,ib,i,iorb

      if (strn /= ' ') call info0(0,0,0,strn)

      ib0 = 0
      do  ib = ib1, ib2
      iorb = mxorb*(ib-1)
      do  i = 1, mxorb
        iorb = iorb+1
        offi = iprmb(iorb)
        if (offi <= ldim) then
          type = 'low'
        elseif (offi <= lidim) then
          type = 'low'
        elseif (offi <= lihdim) then
          type = 'high'
        else
          type = 'neglected'
        endif
        if (ib /= ib0) then
          print 333, ib, iorb, offi, type
          ib0 = ib
        else
          print 334,     iorb, offi, type
        endif
  333   format(i5,2i6,2x,a)
  334   format(5x,2i6,2x,a)
      enddo
      enddo
      end

      subroutine plhamnc(s_ham,s_pot,s_str,s_site,mode,plat,isp,nspc,kcplx,ip0,npl,ld0,lds,pgplp,qp,wk,spl)
C- Layer hamiltonian, including noncollinear case
C ----------------------------------------------------------------------
Ci Inputs
Cio Structures
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  neula ldham lgen3 lncol lham offH
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb eula
Cio    Passed to:  plham plhamso
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pfnc palp
Cio    Passed to:  plham
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iax npr s
Cio    Passed to:  plham
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  norb
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pfr
Cio    Passed to:  plham plhamso
Ci   mode  :=0 return hamiltonian in global axis.
Ci         :>0 rotate to local spin quantization axis.
Ci         :  Note: implementation is incomplete
Ci   plat  :primitive lattice vectors, in units of alat
Ci   isp   :current spin channel (1 or 2)
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   kcplx :distinguishes how complex arithmetic is handled
Ci          0: real, imaginary separated: hpl = hpl(ld0,ldim,1..2)
Ci          1: complex*16: hpl = hpl(ld0,ldim)
Ci          2: real, imaginary in columns : hpl = hpl(ld0,1..2,ldim)
Ci   ip0   :Current principal layer index
Ci   npl   :number of principal layers (pgfset.f)
Ci   ld0   :dimension of this PL
Ci   lds   :dimension of PL + adjacent left PL + adjacent right PL
Ci   pgplp :index and dimensioning information for crystal subblocks.
Ci          The meaning of pgplp depends on the context; see subas.f
Ci   qp    :k-point
Ci   wk    :work array of same dimension as spl.  Not used unless nspc=2.
Co Outputs
Co   spl   :2D Bloch summed P-S, dimensioned(ld0*nspc,lds*nspc)
Cr Remarks
Cr   In the noncollinear case spl for layer i is returned in (orbital,spin,ipl) order:
Cr   (Here j=layer i-1, k=layer i+1, and ldim = ld0+ldl+ldr)
Cr                0     ld0     2*ld0  2*ld0+ldl    2*ld2  2*ld2+ldr   where ld2=ld0+ldl
Cr                |     |         |      |            |      |
Cr     0    ->    \/    \/        \/     \/           \/     \/
Cr                ii++  ii+-      ij++   ij+-         ik++   ik+-
Cr                ii-+  ii--      ij-+   ij--         ik-+   ik--
Cr  *Note: plhamnc and plhamso return P-S in (orbital,spin,ipl) order,
Cr         which is distinct from the (orbital,ipl,spin) order plham generates.
Cr         In the former case spin blocks are contiguous, convenient for computation
Cr   plhamnc calls pgflu2, which rearranges the elements to (orbital,spin,ipl) order.
Cu Updates
Cu   15 Jun 17  New mode to select global or local spin quantization axis
Cu   03 Oct 15  First created, adapted from pgflu
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer isp,nspc,ip0,npl,ld0,lds,kcplx,mode,pgplp(6,-1:*)
      double precision plat(3,3),qp(3)
      double complex wk(ld0,nspc,lds,nspc),spl(ld0,nspc,lds,nspc)
C ... For structures
!      include 'structures.h'
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_str)::   s_str
      type(str_site)::  s_site(*)
C ... Dynamically allocated local arrays
      integer, allocatable :: iprmb(:)
C ... Local parameters
      integer ib01,ib02,ibL1,ibL2,ibR1,ibR2,kcplxl,ldl,ldr,lhdimp,lidimp
      integer mxorb,nbasp,ncsw,neul,nglob,nl,nl2,offL,offR,offi
      double precision xv(1)
      character strn*80
      integer, parameter :: PRTG=81
      procedure(logical) bittst
      procedure(integer) :: iprint
      procedure(real(8)) :: dlength

      mxorb = nglob('mxorb')
      nbasp = nglob('nbasp')
      nl = nglob('nl')
      neul = s_ham%neula
      nl2 = nl*nl
      lhdimp = s_ham%ldham(3)
      lidimp = s_ham%ldham(2)
      ldl = pgplp(4,max(ip0-1,-1))
      ldr = pgplp(4,min(ip0+1,npl))
      if (ld0+ldl+ldr /= lds) call rx('plhamnc: improper arguments')

      if (nspc == 2) then
        kcplxl = 1              ! plham requires kcplx=1
        call plham(1,s_ham,s_pot,s_str,s_site,plat,isp,kcplxl,max(ip0,-1),pgplp,qp,lds,ld0,wk)
C        call awrit1(' P-S(%i) from plham:',outs,len(outs),0,max(ip0,-1))
C        call zprm(outs,2+0,wk,ld0*nspc,ld0*nspc,lds*nspc)

        if (mode /= 0) then ! Rotate to local spin quantization axis

C   ... Rotate P-S by Euler angles
C       Convert raw hamiltonian offsets offi,offL,offR into layer form
C       offi -> 0, offL -> ld0, offR -> ld0+ldl
        call pghoff(ip0,npl,pgplp,s_ham%iprmb,mxorb,offL,offi,offR)
        call gtibpl(ip0,npl,pgplp,ib01,ib02)
        call gtibpl(max(ip0-1,-1),npl,pgplp,ibL1,ibL2)
        call gtibpl(min(ip0+1,npl),npl,pgplp,ibR1,ibR2)
C       Modify iprmb by adding largest index to exclude all pairs but restore PL, PL-1, PL+1
        offi  = offi + lhdimp
        offL  = offL - ld0 + lhdimp
        offR  = offR - ld0 - ldl + lhdimp
        allocate(iprmb(size(s_ham%iprmb))); call icopy(size(iprmb),s_ham%iprmb,1,iprmb,1)
        call pblch2(lhdimp,1,nbasp,nl2,iprmb)
        call pblch2(-offi,ib01,ib02,nl2,iprmb)
        call pblch2(-offL,ibL1,ibL2,nl2,iprmb)
        call pblch2(-offR,ibR1,ibR2,nl2,iprmb)

C       Rotate the three blocks
        ncsw = 30000 + 10*kcplxl
C       call dcopy(ld0*nspc*lds*nspc*2,wk,1,wk1,1)
C       call rx('plham check sign of rotation')
C       call rotspn(ncsw,1,nbasp,nbasp,nl,iprmb,s_ham%eula,neul,xv,xv,xv,ld0,ld0,lds,ld0,lds,wk1,wk)

C       Until change is made, NC and SO are not compatible.  Don't turn this check on to allow end layers to work
        if (bittst(s_ham%lncol,4)) then  ! SO is turned on
          if (dlength(size(s_ham%eula),s_ham%eula,1) /= 0)
     .      call rx('plhamnc: not implemented with local rotations including orbital rotation')
        endif
C       rotheu not ready for block s yet
C       optrot = 100*kcplx*0 + 1 ! Backwards rotation
C       if (bittst(s_ham%lham,256)) optrot = optrot + 10 ! Spherical harmonics
C       call rotheu(optrot,s_ham,nl,min(ib01,ibL1,ibR1),min(ib02,ibL2,ibR2),nbasp,0,ld0,ld0,lds,wk)
        call rx('plhamnc check sign of rotation, call to rotspn')
        call rotspn(ncsw,1,nbasp,nbasp,nl,iprmb,s_ham%eula,neul,xv,xv,xv,ld0,ld0,lds,ld0,lds,wk,wk)

        deallocate(iprmb)
C        call awrit1(' rot P-S(%i):',outs,len(outs),0,max(ip0,-1))
C        call zprm(outs,2+0,wk,ld0*nspc,ld0*nspc,lds*nspc)
        endif ! Rotation to local spin quantization axis

C   ... Rearrange (spin,orbital) indexing
        call dpzero(spl,2*size(spl))
        call pgflu2(0,kcplxl,nspc,ld0,ldl,ldr,lds,0,0,2,2,ld0*nspc,wk,spl)
C        call awrit1(' spin reordered P-S(%i):',outs,len(outs),0,max(ip0,-1))
C        call zprm(outs,2+0,spl,ld0*nspc,ld0*nspc,lds*nspc)

C       Convert to kcplx format
        call ztoyy(spl,ld0*nspc,lds*nspc,ld0*nspc,lds*nspc,kcplxl,kcplx)
      else
        call plham(0,s_ham,s_pot,s_str,s_site,plat,isp,0,max(ip0,-1),pgplp,qp,lds,ld0,spl)
      endif

      if (iprint() >= PRTG/1) then
        call yprmi('plhamnc: P-S for layer %i',max(ip0,-1),0,kcplx+2,
     .    spl,ld0*lds*nspc**2,ld0*nspc,ld0*nspc,lds*nspc)
      endif

      end
