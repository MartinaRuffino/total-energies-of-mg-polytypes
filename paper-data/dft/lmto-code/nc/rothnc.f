      integer function rothnc(opt,nl,nspc,nbas,bas,iat,offH,indxsh,
     .  istab,g,ag,q,ldha,ldhb,h0,h)
C- Transform hamiltonian h0(g^-1 q) into h(q)
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :1s digit:
Ci         :0, rotate entire h
Ci         :1  rotate site-diagonal h
Ci         :2  no longer supported
Ci         :3  rotate rows and columns connected to a particular site
Ci             ib only.  Site is specified by iat
Ci         :4  same as case 3, except that h is returned untouched;
Ci             and rotated h is returned in h0; see Remarks
Ci         :5  same as case 1, except that h consists only of
Ci             diagonal blocks.
Ci         10s digit
Ci         :0  h is in real spherical harmonics
Ci         :Add 2 if h is in true (complex) spherical harmonics
Ci         :    See rcub2sph.f for definition of true spherical harmonics.
Ci         :Add 4 if phase convention phi = q * [(g R_j + a) - R_i]
Ci         :  should be scaled by -1
Ci        100s digit distinguishes how complex arithmetic is handled
Ci           0: h,h0 have real, imaginary separated
Ci              h = h(ldha,ldhb,2), with h(*,*,1..2) = real..imag
Ci           1: h,h0 are in complex*16 format (see Bugs)
Ci              h = h(2,ldha,ldhb), with s(1,*) = real, s(2,*) = imag
Ci           2: h,h0 have real, imaginary separated by columns
Ci              h = h(ldha,2,ldhb), with h(*,1..2,*) = real..imag
Ci       1000s digit
Ci         :For nonzero values, entire row dimension is rotated,
Ci         :but row dim of h may only consist of some subblock:
Ci         :0 row dimension consists only of lower block
Ci         :1 row dimension consists only of lower block (same as 0)
Ci         :2 row dimension consists only of intermediate block
Ci         :3 row dimension consists only of higher block
Ci         :4 row & col dimension consists lower+intermediate block
Ci      10000s digit
Ci         :Similar to 1000s digit, but applicable to both row
Ci         :and column dimensions.
Ci         :2 rotate both row and col dim in intermediate block only
Ci         :  NB: implies 1000 digit = 2.
Ci         :3 rotate both row and col dim in higher block only
Ci         :  NB: implies 1000 digit = 3.
Ci      100000s digit
Ci         :1 suppress spinor rotation (this was the default prior to v 7.11b)
Ci         :2 Reverse the sense of the spinor rotation
Ci   nl    :(global maximum l) + 1
Ci   nbas  :number of atoms in the basis (input)
Ci   bas   :basis vectors (input)
Ci   iat   :if 1s digit opt=0 or opt=1, not used
Ci         :if 1s digit opt=3, h(1..nbas,iat) h(iat,1..nbas) are rotated
Ci   offH  :Offsets to hamiltonian matrix (makidx.f)
Ci   indxsh:permutations ordering orbitals in lower, intermediate,
Ci          and higher blocks (makidx.f)
Ci   istab :table of site permutations for each group op (symtbl.f)
Ci          Site istab(i,ig) is transformed into site i by grp op ig.
Ci          Obtain by calling symtbl with mode=1.
Ci          NB: istab(1)<0 flags that current symop is consists of a proper
Ci               symmetry operation * inversion; see Remarks
Ci          Special case
Ci   g,ag  :space group operation; see Remarks
Ci   q     :qp for which h is sought.
Ci   ldha  :leading dimension of h
Ci   ldhb  :second dimension of h
Ci   h0    :unrotated hamiltonian, h(q0), where q0 = g^-1 q
Cio Inputs/Outputs
Co Outputs
Co    h    :h is rotation of h0 = h(q0) to h(q).
Co  rothnc : 0 normal return
Co         :-1 transpose required but protht did not make it.
Co         : 1 protht took transpose
Co         :   -1 and 1 may be returned only for special cases that arise when the
Co         :   operation is not a true symop but present because of time reversal symmetry.
Co         :   protht cannot handle certain special modes, e.g. opt=2000, 2003, 2004, which see.
Cl Local variables
Cl   rmat  :rotation matrix.
Cl         :rmat is real for real harmonics, complex for spherical harmonics
Cl   ld0r  :leading dimension of hwk
Cl   ldhr  :leading dimension of h
Cl   ldmpa :offset to 1st orbital in downfolding block (row index)
Cl   ldmpb :offset to 1st orbital in downfolding block (col index)
Cl   of0i   :offset to imaginary part of hwk
Cl   ofhi  :offset to imaginary part of h
Cl   ... For opt0 gt 1:
Cl   ibla  :index to offH defining block to be rotated (row index)
Cl   iblb  :index to offH defining block to be rotated (col index)
Cl   nlma  :number of orbitals connected to site ib (row index)
Cl   nlma1 :number of orbitals in lower block (row index)
Cl   nlma2 :number of orbitals in upper block (row index)
Cl   nlmb  :number of orbitals connected to site ib (col index)
Cl   nlmb1 :number of orbitals in lower block (col index)
Cl   nlmb2 :number of orbitals in upper block (col index)
Cl   ofa1  :offset to lower block belonging to ib (row index)
Cl   ofa2  :offset to intermediate block belonging to ib (row index)
Cl   ofb1  :offset to lower block belonging to ib (col index)
Cl   ofb2  :offset to intermediate block belonging to ib (col index)
Cl   lroti :set to .true. if rotation consists of proper symmetry op
Cl         :* inversion.  (Flagged by istab(1)<0).  See Remarks.
Cr Remarks
Cr   Input hamiltonian h(q0) is rotated to h(q), where q0 = g^-1 q.
Cr   roth applies the rotation R making  h(q) = R h0 (g^-1 q) R^-1.
Cr
Cr   Phase shifts:  Space group op (g,ag) rotates basis position R_j
Cr   to R'_j = g R_j + ag.  It is possible that R'_j and R_j differ
Cr   by a lattice translation vector.  Phase shifts are calculated
Cr   to remove the shift in generated by ( g R_j + ag - R_i ).
Cr
Cr  *Extended symmetry operations.  In certain cases the inversion
Cr   operation may be added to the space group (see mksym.f).
Cr   (This is flagged when istab(1)<0).
Cr   In that case, (-g,ag) should be a proper group operation.
Cr   roth simply uses (-g,ag) for the symmetry operation, and returns
Cr   h (or h0) for that symop.  However, the correct h (or h0) to use
Cr   is the transpose of that generated by initially by roth (see aginv).
Cr   roth calls a kernel, protht, afterwards to generate the transpose.
Cr   But because of its internal complexity, protht cannot handle all
Cr   cases; roth returns one of (0;1;-1) indicating one of (no transpose
Cr   needed; transpose taken; transpose needed but not taken).
Cr   See protht for a further discussion.
Cr
Cr   Clarification of the sense of rotation.
Cr   This g generates a 90 degree counterclockwise rotation
Cr   of coordinates mapping (x,y) into (y,-x):
Cr            ( 0   1   0)
Cr        g = (-1   0   0)
Cr            ( 0   0   1)
Cr   g generates this pp block of the rotation of real harmonics:
Cr            (0    0   -1)
Cr     R_pp = (0    1    0)
Cr            (1    0    0)
Cr   Note Y_1-1,Y_10,Y_11 correspond to y,z,x respectively.
Cr   Then
Cr            (               )     (                )
Cr       R_pp (Y1-1  Y_10  Y11)  =  (-Y11  Y_10  Y1-1)
Cr            (               )     (                )
Cr   (px,py) are transformed (py,-px)
Cu Updates
Cu   09 May 15 Fixes so rothnc works in the kcplx=2 case
Cu   13 Aug 13 Bug fix: case spherical harmonics in conjunction w/ kcplx=1
Cu   17 Jun 02 Added option for phase convention (10s digit opt)
Cu   05 Feb 01 Added code to handle special symmetry operations:
Cu             (proper symmetry operations + inversion)
Cu             Turned roth into an integer function call
Cu   14 Mar 00 some bug fixes in downfolding cases.
Cu    2 Mar 00 Added 10000s digit options to opt
Cu   29 Oct 99 generalized in several respects; argument list changed.
Cu   10 Dec 99 handles h consisting of lower+intermediate blocks.
Cb Bugs
Cb  *kcplx=1 complex storage makes an internal transformation to kcplx=2
Cb   and thus requires additional O(N^2) operations
Cb  *roth does not work properly for 1s digit opt = 1 and sil case
Cb  *unit rotation: rothnc copies h0 to h and returns
Cb
Cb   This routine is too complicated!  It handles too many cases.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer opt,ldha,ldhb,nl,nbas,istab(nbas),offH(n0H,nkap0,*),
     .  iat,indxsh(1),nspc
      real(8) :: bas(3,*),g(3,3),ag(3),q(3)
      real(8) :: h0(ldha,nspc,ldhb,nspc),h(ldha,nspc,ldhb,nspc)
C ... Dynamically allocated arrays
      real(8),allocatable :: hwk(:,:,:)
C ... Local parameters
C     logical :: debug=.false.
      integer, parameter :: nlx=7, n0=10
      logical lunitg,lroti
      integer ib,ibla,iblb,is,js,k1,k2,kcplx,ld0r,ldh1,ldh2,ldhax,ldhbx,
     .  ldhr,ldima,ldimb,ldmpa,ldmpb,modspn,nl2,nlma,nlma1,nlma2,nlmb,
     .  nlmb1,nlmb2,of0i,ofa1,ofa2,ofb1,ofb2,ofhi,ofijs,opt0,opt1,morder
      double precision ddot,xx,eula(3),det33,pi
      double precision gp(3,3),wk(nlx**4),rmat(nlx**4*2)
      double complex u(2,2,nl*nl,nbas)
      procedure(integer) :: protht,roth

      if (mod(opt,10) == 4) call rx('rothnc cannot handle opt mode 4')
      ldhax = ldha*nspc
      ldhbx = ldhb*nspc
      call dcopy(2*ldhax*ldhbx,h0,1,h,1)
      pi = 4d0*datan(1d0)
      call lmorder(0,morder,[0],[0])

C ... lunitg is true if unit rotation matrix.  Quit if nothing to rotate
      lunitg = abs(ddot(9,g,1,g,1)-3) < 1d-10 .and.
     .         abs(g(1,1)-1) < 1d-10 .and.
     .         abs(g(2,2)-1) < 1d-10 .and.
     .         abs(g(3,3)-1) < 1d-10
      if (lunitg) then
        rothnc = 0
        return
      endif

C ... Collinear case: rothnc simplifies to roth
      if (nspc == 1) then
        allocate(hwk(ldha,ldhb,2))
        rothnc = roth(opt,nl,nbas,bas,iat,offH,indxsh,istab,g,ag,q,ldha,ldhb,hwk,h)
        deallocate(hwk)
        return
      endif

      call tcn('rothnc')

C     call zprm('starting h',2,h0,ldhax,ldhax,ldhbx)
      if (nl > nlx) call rxi('increase nlx in rothnc, need',nl)
      nl2 = nl*nl
      opt0 = mod(opt,10)
      if (opt0 == 2) call rx('rothnc: 1s digit mode 2 no longer supported')
      kcplx = mod(opt/100,10)
      if (kcplx == 1) then
        kcplx = 2
        call ztoyy(h,ldhax,ldhbx,ldhax,ldhbx,1,2)
      endif


C ... Get some dimensioning parameters
C     call cplxdm(kcplx,ldhax,ldhbx,ldh1,ldh2,ldhr,ofhi)
      opt1 = mod(mod(opt/10,10),4)
      call pvrotd(opt,offH,nbas,iat,ibla,iblb,ldmpa,ldima,ldmpb,
     .  ldimb,ofa1,ofa2,nlma,nlma1,nlma2,ofb1,ofb2,nlmb,nlmb1,nlmb2)

C ... Setup for pseudorotation (rotation + inversion)
      lroti = .false.
      if (istab(1) < 0) then
        istab(1) = -istab(1)
        call dscal(9,-1d0,g,1)
        call dscal(3,-1d0,q,1)
        lroti = .true.
      endif

C --- Set up Y_lm and spinor rotation matrix for this g ---
      if (opt1 == 0) then
        call ylmrtg(nl2,g,rmat)
C       call yprm('g',1,g,0,3,3,3)
C       call yprm('rmat',1,rmat,0,nl2,nl2,nl2)
      else
        call ylmrtg(nl2,g,wk)
C       prothl and prothr expect rmat in kcplx=0 mode
C       if (morder < 2) morder = 3-morder ! This routine requires rmat = cc of rotation matrix
        call s2sph(0+100*morder,nl,nl,wk,nl2,nl2,nl2,nl2,rmat)
C       call yprm('rmat',2,rmat,nl2*nl2,nl2,nl2,nl2)
      endif

C ... Spinor rotation: 2x2 rotation for each (L,L')
      if (mod(mod(opt/100000,10),2) == 0) then  ! skip if to suppress spinor part
        gp = g
        if (det33(g) < 0) call dscal(9,-1d0,gp,1)  ! Improper rotation
        call rm2eua(gp,eula(1),eula(2),eula(3))
C        if (det33(g) < 0) then ! rotate by 180 degrees ... no.  Preserves MX in
C                                  ! cubic symmetry, constrained to preserve z quantization
C                                  ! axis. But MX should flip (up,down) spins.
C          eula(2) = pi - eula(2)
C          eula(1) = pi + eula(1)
C        endif
        if (mod(opt/100000,10) == 2) then
C         Not clear that this branch should ever be used ... see Remarks in rotsp1
C         call rx('rothnc: check sense of spinor rotation')
          call rotspu(1,1,1,1,nl,eula,1,u) ! Spinor rotation matrix, reversed sense
        else
          call rotspu(0,1,1,1,nl,eula,1,u) ! Spinor rotation matrix
        endif
        do  ib = 2, nbas
          u(:,:,:,ib) = u(:,:,:,1)
        enddo
C       Debugging check: verify that rotsp1 works for kcplx=1,2
C        call yprm('before spinor rotation',2,h,ldhax**2,ldhax,ldhax,ldhax)
C        kcplx = 2
C        call ztoyy(h,ldhax,ldhbx,ldhax,ldhbx,mod(opt/100,10),kcplx)
C        call yprm('before spinor rotation',2+kcplx,h,ldhax,ldhax,ldhax,ldhax)
C        modspn = 10000 + 10*kcplx
C        call rotsp1(modspn,1,nbas,nbas,nl,indxsh,u,xx,xx,1,ldima,ldimb,ldha,
C     .    ldhb,xx,xx,xx,h,h,h)
C        call yprm('after spinor rotation',2+kcplx,h,ldhax,ldhax,ldhax,ldhax)

        modspn = 10000 + 10*kcplx ! input in spinor form
C       call yprm('before spinor rotation',2,h,ldhax**2,ldhax,ldhax,ldhax)
        call rotsp1(modspn,1,nbas,nbas,nl,indxsh,u,xx,xx,1,ldima,ldimb,ldha,
     .    ldhb,xx,xx,xx,h,h,h)
C       call yprm('after spinor rotation',2,h,ldhax**2,ldhax,ldhax,ldhax)
      endif

C --- For each spinor component, do ---
      call cplxdm(kcplx,ldha,ldhb,ldh1,ldh2,ld0r,of0i)   ! ld0r,of0i for hwk
      call cplxdm(kcplx,ldhax,ldhbx,ldh1,ldh2,ldhr,ofhi) ! ldhr,ofhi for h
      allocate(hwk(ldha,ldhb,2))
      do  is = 1, nspc
      do  js = 1, nspc

C       ofijs = (is-1 + (js-1)*ldhax)*ld0r
        ofijs = (is-1)*ldha + (js-1)*ldhax*ld0r ! Works for kcplx=0 and 2


C   ... Rotate rows of h into hwk: jbas rotated into ibas
        if (opt0 == 1 .and. (ibla == 2 .or. ibla == 3)) then
        if (iblb == 1) call rxi('rothnc not implemented for mode',opt)
        endif
        k1 = 1
        k2 = ldimb-ldmpb

C        if (debug) call yprmi('H before rotation is,js=%2:1i',[is,js],0,kcplx+2,
C     .    h(1+ofijs,1,1,1),ofhi,ldhax,ldha,ldhb)


C       hwk <- R h
        call prothl(opt,nl2,nbas,bas,ldmpa,ldima,indxsh,istab,g,ag,q,
     .    rmat,k1,k2,1,nbas,0,ld0r,-1,0,hwk,of0i,ldhr,ldhbx,h(1+ofijs,1,1,1),ofhi)
C       if (debug) call yprmi('R H is,js=%2:1i',[is,js],0,kcplx+2,hwk,of0i,ldha,ldha,k2-k1+1)

C   ... Rotate columns of hwk into h
        k1 = 1
        k2 = offH(ibla,1,nbas+1)
        call prothr(opt,nl2,nbas,bas,ldmpb,ldimb,indxsh,istab,g,ag,q,
     .    rmat,k1,k2,1,nbas,0,ld0r,-1,0,hwk,of0i,ldhr,ldhbx,h(1+ofijs,1,1,1),ofhi)
C        call snot(h,h,ldha,0)

C        if (debug) call yprmi('R H R^-1 is,js=%2:1i',[is,js],0,kcplx+2,
C     .    h(1+ofijs,1,1,1),ofhi,ldhax,ldha,ldhb)


      enddo ! Loop over spinor index
      enddo ! Loop over spinor index
      deallocate(hwk)
C      call yprm('R H R^-1',kcplx+2,h,ofhi,ldhax,ldhax,ldhbx)

C --- Cleanup for improper rotation (rotation + inversion) ---
      if (lroti) then
        istab(1) = -istab(1)
        call dscal(9,-1d0,g,1)
        call dscal(3,-1d0,q,1)
      endif

      kcplx = mod(opt/100,10)
      if (kcplx == 1) then
        call ztoyy(h,ldhax,ldhbx,ldhax,ldhbx,2,1)
      endif
      call tcx('rothnc')

      allocate(hwk(1,1,1))  ! So the compiler doesn't complain
      rothnc = protht(opt,nbas,iat,offH,indxsh,istab,ldhax,ldhbx,hwk,h)
      deallocate(hwk)

C     call zprm('ending h',2,h,ldhax,ldhax,ldhbx)
      end
C      subroutine snot(h1,h2,ldha,ofhi)
C      implicit none
C      integer ofhi,ldha
C      double precision h1(ldha,2,2,ldha,2),h2(ldha,2,ldha,2,2)
C
C      print *, h1(2,1,1,2,1),h1(2,1,1,2,2)
C      print *, h1(2,2,1,2,1),h1(2,2,1,2,2)
C
C      print *, h1(2+ldha*4,1,1,1,1),h1(2+ldha*ldha*4+ldha*4,1,1,1,1)
C      print *, h1(2+ldha+ldha*4,1,1,1,1),h1(2+ldha+ldha*ldha*4+ldha*4,1,1,1,1)
C
C      print *, 0,ldha*ldha*4
C      print *, ldha,ldha+ldha*ldha*4
C
CC      print *, h2(2,1,2,1,1),h2(2,1,2,2,1)
CC      print *, h2(2,2,2,1,1),h2(2,2,2,2,1)
CC
CC      print *, h2(2+ldha*2,1,1,1,1),h2(2+ldha*ldha*2+ldha*2,1,1,1,1)
CC      print *, h2(2+ldha+ldha*2,1,1,1,1),h2(2+ldha+ldha*ldha*2+ldha*2,1,1,1,1)
CC
CC      print *, 0,ldha*ldha*2
CC      print *, ldha,ldha+ldha*ldha*2
C
CC      stop
C
C      end
