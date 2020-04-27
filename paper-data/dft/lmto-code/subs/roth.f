      integer function roth(opt,nl,nbas,bas,iat,offH,indxsh,istab,g,ag,
     .  q,ldha,ldhb,h0,h)
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
Ci         :    See rcub2sph for definition of true spherical harmonics
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
Cio Inputs/Outputs
Cio   h    :on input, unrotated hamiltonian, h(q0), where q0 = g^-1 q
Cio        :on output, h is rotated from h(q0) to h(q).
Cio        :Note: 1s digit opt=4 => h untouched (rotated h returned in h0)
Cio   h0   :complex array, used as a workspace except for opt=4.
Cio        :h0 has the following dimensions:
Cio        :opt = 0,1 same dimensions as h
Cio        :opt = 3   nlma*ldhb + ldha*nlmb, with
Cio        :          nlma = row dimension of site iat
Cio        :          nlmb = column dimension of site iat
Cio        :opt = 4   2*dim for opt=3 case = 2*(nlma*ldhb + ldha*nlmb)
Cio        :          This is a special case: rotated h is returned in h0
Cio        :          For this case:
Cio        :h0(1..nlma*ldhb) holds rows of h corresponding to site iat
Cio        :h0 is an array dimensioned h0(nlma,ldhb), with
Cio        :          nlma = row dimension of site iat
Cio        :h0(1+off2..ldha*nlmb) holds cols of h b belonging to site iat
Cio        :h0(1+off2) is an array dimensioned h0(ldha,nlmb), with
Cio        :          nlmb = column dimension of site iat and
Cio        :          off2 = nlma*ldhb*2
Co Outputs
Co    roth :returns  0 if no transpose required
Co         :returns -1 transpose required but protht did not make it
Co         :returns  1 protht took transpose
Cl Local variables
Cl   rmat  :rotation matrix.
Cl         :rmat is real for real harmonics, complex for spherical harmonics
Cl   ldhr  :leading dimension of h
Cl   ldmpa :offset to 1st orbital in downfolding block (row index)
Cl   ldmpb :offset to 1st orbital in downfolding block (col index)
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
Cr   Clarification of the sense of rotation. The following g corresponds to
Cr   a 90 degree counterclockwise rotation of the coordinate system
Cr   generated by rotation: z:pi/2  (right-hand rule)
Cr            ( 0   1   0)
Cr        g = (-1   0   0)
Cr            ( 0   0   1)
Cr   The x-axis xhat gets mapped into yhat, the y-axis mapped into -xhat.  Alternatively,
Cr   a vector initially parallel to xhat (yhat) becomes parallel to -yhat (+xhat).
Cr   This is the "intrinsic" sense defined in Wikipedia (en.wikipedia.org/wiki/Euler_angles)
Cr
Cr   This g generates the following pp block of the rotation of real harmonics:
Cr            (0    0   -1)
Cr     R_pp = (0    1    0)
Cr            (1    0    0)
Cr   (Note Y_1-1,Y_10,Y_11 correspond to y,z,x respectively.)
Cr   R_pp maps  px -> -py, and py -> px.    In other words,  a vector of
Cr   coefficients to a linear combination of Y_1-1,Y_10,Y_11 is transformed as
Cr
Cr            (C_1-1)   (-C_11)
Cr       R_pp (C_10 ) = (C_10)
Cr            (C_11 )   (C_1-1)
Cr
Cr   The p block of structure matrix S consisting only of S_yz -> S_xz
Cr    (0    0    0)    (0    0    0)          (0    0    0)    ( 0    0    0)
Cr    (1    0    0) -> (0    0    1)  whereas (0    0    1) -> (-1    0    0)
Cr    (0    0    0)    (0    0    0)          (0    0    0)    ( 0    0    0)
Cr
Cu Updates
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
Cb   kcplx=1 complex storage makes an internal transformation to kcplx=2
Cb   and thus requires additional O(N^2) operations
Cb   roth does not work properly for 1s digit opt = 1 and sil case
Cb
Cb   This routine is too complicated!  It handles too many cases.
C ----------------------------------------------------------------------
      implicit none
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer opt,ldha,ldhb,nl,nbas,istab(nbas),offH(n0H,nkap0,*),
     .  iat,indxsh(*)
      double precision bas(3,*),g(3,3),ag(3),q(3),h0(*),h(ldha,ldhb)
C Local variables
      integer nlx,n0
      parameter (nlx=7,n0=10)
C     logical :: debug=.false.
      logical lsaveh,lunitg,lroti
      integer ibla,iblb,jb1,k1,k2,kcplx,ld01,ld02,ld0b,ib1,
     .  ld0r,ldh1,ldh2,ldhr,ldima,ldimb,ldmpa,ldmpb,nl2,nlma,nlma1,
     .  nlma2,nlmb,nlmb1,nlmb2,of0i,ofa,ofa1,ofa2,ofb,ofb1,ofb2,off1,
     .  off2,ofhi,opt0,opt1,morder
      integer norb,offli(n0*nkap0)
      double precision wk(nlx**4),rmat(nlx**4*2),xx
      procedure(integer) :: protht
      procedure(real(8)) :: ddot

      call tcn('roth')
C     call zprm('starting h',2,h,ldha,ldha,ldhb)
      call lmorder(0,morder,[0],[0])

      if (nl > nlx) call rxi('increase nlx in roth, need',nl)
      nl2 = nl*nl
      opt0 = mod(opt,10)
      if (opt0 == 2)
     .  call rx('roth: 1s digit mode 2 no longer supported')
C     lsaveh is true when h is preserved
      lsaveh = opt0 == 4
      if (lsaveh) opt0 = 3
      if (opt0 == 5) opt0 = 1
      kcplx = mod(opt/100,10)
      if (kcplx == 1) then
        kcplx = 2
        call ztoyy(h,ldha,ldhb,ldha,ldhb,1,2)
      endif

C ... Get some dimensioning parameters
      call cplxdm(kcplx,ldha,ldhb,ldh1,ldh2,ldhr,ofhi)
      opt1 = mod(mod(opt/10,10),4)
      call pvrotd(opt,offH,nbas,iat,ibla,iblb,ldmpa,ldima,ldmpb,
     .  ldimb,ofa1,ofa2,nlma,nlma1,nlma2,ofb1,ofb2,nlmb,nlmb1,nlmb2)

C     lunitg is true if unit rotation matrix
      lunitg = ddot(9,g,1,g,1)-3 < 1d-10 .and.
     .             abs(g(1,1)-1) < 1d-10 .and.
     .             abs(g(2,2)-1) < 1d-10 .and.
     .             abs(g(3,3)-1) < 1d-10

C --- Hold on to rows and columns of unrotated h, to restore later ---
      if (lsaveh) then
        off1 = (nlma*ldhb + ldha*nlmb)*2
C       Just copy to final position if unit rotation matrix
        if (lunitg) off1 = 0
        off2 = off1 + nlma*ldhb*2
C       Save rows of h(iat) into h0
        call cplxdm(kcplx,nlma,ldhb,ld01,ld02,ld0r,of0i)
C       nlma1 rows in lower block, ldhb columns (l+i for entire hamiltonian)
        call ymscop(0,nlma1,ldhb,ldhr,ld0r,ofa1,0,0,0,h,ofhi,h0(1+off1),of0i)
C       nlma2 rows in intermediate block, ldhb columns
        call ymscop(0,nlma2,ldhb,ldhr,ld0r,ofa2,0,nlma1,0,h,ofhi,h0(1+off1),of0i)
C       Save columns of h(iat) into h0
        call cplxdm(kcplx,ldha,nlmb,ld01,ld02,ld0r,of0i)

C       ldha rows (l+i for entire hamiltonian), nlmb1 columns in lower block
        call ymscop(0,ldha,nlmb1,ldhr,ld0r,0,ofb1,0,0,h,ofhi,h0(1+off2),of0i)
C       ldha rows (l+i for entire hamiltonian), nlmb2 columns in intermediate block
        call ymscop(0,ldha,nlmb2,ldhr,ld0r,0,ofb2,0,nlmb1,h,ofhi,h0(1+off2),of0i)
C       call yprm('h0 (r)',kcplx+2,h0,nlma*ldhb,nlma,nlma,ldhb)
C       call yprm('h0 (c)',kcplx+2,h0(1+off2),ldha*nlmb,ldha,ldha,nlmb)
      endif

C     Nothing to rotate; just quit
      if (lunitg) goto 99

C --- Setup for pseudorotation (rotation + inversion) ---
      lroti = .false.
      if (istab(1) < 0) then
        istab(1) = -istab(1)
        call dscal(9,-1d0,g,1)
        call dscal(3,-1d0,q,1)
        lroti = .true.
      endif

C --- Set up Y_lm rotation matrix ---
      if (opt1 == 0) then
        call ylmrtg(nl2,g,rmat)
C       call yprm('g',1,g,0,3,3,3)
C       call yprm('rmat',1,rmat,0,nl2,nl2,nl2)
      else
        call ylmrtg(nl2,g,wk)
C        call yprm('g',1,g,0,3,3,3)
C        call yprm('rmat',1,wk,0,nl2,nl2,nl2)
C       prothl and prothr expect rmat in kcplx=0 mode
C       if (morder < 2) morder = 3-morder ! If this routine required rmat = cc of rotation matrix
        call s2sph(0+100*morder,nl,nl,wk,nl2,nl2,nl2,nl2,rmat)
C       call yprm('rmat',2,rmat,nl2*nl2,nl2,nl2,nl2)
      endif

C --- Rotate rows of h into h0 = R h: jbas rotated into ibas ---
C     opt0=3 rotates rows and columns connected to a particular site.
C     This is accomplished in two function calls.
C     Let jb=istab(ib), and ib1=iat with jb1=istab(ib1)
C     call 1: h0(ib1,k) are rotated from h(jb1,k) for k=1..ldimb
C     call 2: h0(ib,k) are rotated from h(jb,k) for ib=1..nbas, k in jb1
      if (opt0 == 3) then
        ib1 = iat
        k1 = 1
        k2 = ldimb-ldmpb
        off1 = 0
        call cplxdm(kcplx,nlma,ldhb,ld01,ld02,ld0r,of0i)
        call iprmsb(10,ib1,ib1,nbas,offH,indxsh,xx,offli,norb)
        call prothl(opt,nl2,nbas,bas,ldmpa,ldima,indxsh,istab,g,ag,q,
     .    rmat,k1,k2,ib1,ib1,0,ld0r,offli,0,h0,of0i,ldhr,ldhb,h,ofhi)
C       if (debug) call yprm('R H (r)',kcplx+2,h0,of0i,nlma,nlma,ldhb)

C   ... Rotate into h0(ib,k) for ib=1..nbas, k in jb1
C       offset = size needed to hold all rows of h above
        off2 = off1 + nlma*ldhb*2
        jb1 = istab(iat)
C       k1,k2 encompass iblb-waves for site jb1
        k1  = offH(iblb,1,jb1)+1
        k2  = offH(iblb,1,jb1+1)
C       Leading dimension encompasses iblb-waves for site jb1
        ld0b = k2-k1+1
C       Leading dimension encompasses l+i-waves for site jb1
        if (ibla == 4) ld0b = offH(ibla,1,jb1+1) - offH(ibla,1,jb1)
C       ofb is col offset in h0 so that R(h(*,jb1)) placed at h0 start
        ofb = k1-1
        call cplxdm(kcplx,ldha,ld0b,ld01,ld02,ld0r,of0i)
        call prothl(opt,nl2,nbas,bas,ldmpa,ldima,indxsh,istab,g,ag,q,
     .    rmat,k1,k2,1,nbas,0,ld0r,[-1],ofb,h0(1+off2),of0i,ldhr,
     .    ldhb,h,ofhi)
C       Rotate il block, when h is composite
        if (ibla == 4) then
C         k1,k2 encompass iwaves for site jb1
          k1  = offH(1,1,nbas+1) + offH(2,1,jb1)+1
          k2  = offH(1,1,nbas+1) + offH(2,1,jb1+1)
C         ofb is col offset subtracted from k1<=k<=k2 when rotating
C         h(*,k) -> h0(P(*),k-ofb).  ofb = first iwave col in h,
C         + size of l-block so as to preserve contents of l-block
          ofb = k1-1 - (offH(1,1,jb1+1)-offH(1,1,jb1))
          call prothl(opt,nl2,nbas,bas,ldmpa,ldima,indxsh,istab,g,ag,
     .      q,rmat,k1,k2,1,nbas,0,ld0r,[-1],ofb,h0(1+off2),of0i,ldhr,
     .      ldhb,h,ofhi)
        endif
C       call yprm('R H (c)',kcplx+2,h0(1+off2),ldha*ld0b,ldha,ldha,ld0b)
      else
        if (opt0 == 1 .and. (ibla == 2 .or. ibla == 3)) then
          if (iblb == 1) call rxi('roth not implemented for mode',opt)
        endif
        k1 = 1
        k2 = ldimb-ldmpb
        call cplxdm(kcplx,ldha,ldhb,ld01,ld02,ld0r,of0i)
        call prothl(opt,nl2,nbas,bas,ldmpa,ldima,indxsh,istab,g,ag,q,
     .    rmat,k1,k2,1,nbas,0,ld0r,[-1],0,h0,of0i,ldhr,ldhb,h,ofhi)
C       call yprm('R H',kcplx+2,h0,ldha*ldhb,ldha,ldha,ldha)
C       call yprm('R H',kcplx+2,h0,ldha*ldhb,ldha,ldha,k2-k1+1)
      endif

C --- Rotate columns of h0 = R h making R H R^-1 and store in h ---
C     opt0=3 rotates rows and columns connected to a particular site.
C     This is accomplished in two passes.
C     Let jb=istab(ib), and ib1=iat with jb1=istab(ib1)
C     pass 1: h(k,ib1) are rotated from h0(k,jb1) for k=1..ldima
C     pass 2: h(k,ib) are rotated from h0(k,jb) for ib=1..nbas, k in ib1
      if (opt0 == 3) then
C   ... Rotate h0 into h(k,ib1) for k = 1 .. ldima
        ib1 = iat
C       jb1 = istab(ib1)
        k1 = 1
        k2 = offH(ibla,1,nbas+1)
        call cplxdm(kcplx,ldha,ld0b,ld01,ld02,ld0r,of0i)
        call iprmsb(10,ib1,ib1,nbas,offH,indxsh,xx,offli,norb)
        call prothr(opt,nl2,nbas,bas,ldmpb,ldimb,indxsh,istab,g,ag,q,
     .    rmat,k1,k2,ib1,ib1,0,ld0r,offli,0,h0(1+off2),of0i,ldhr,ldhb,h,
     .    ofhi)

C       if (debug) call yprm('R H R^-1 (a)',kcplx+2,h,ofhi,ldha,ldhb,ldhb)
C   ... Rotate h0 into h(k,ib) for ib=1..nbas, k = in ib1
C       k1,k2 encompass l-waves for site ib1
        k1  = offH(ibla,1,ib1)+1
        k2  = offH(ibla,1,ib1+1)
        if (ibla == 4) then
          k1  = offH(1,1,ib1)+1
          k2  = offH(1,1,ib1+1)
        endif
        call cplxdm(kcplx,nlma,ldhb,ld01,ld02,ld0r,of0i)
C       ofa is col offset in h0 so that R(h(*,jb1)) placed at h0 start
        ofa = k1-1
        call prothr(opt,nl2,nbas,bas,ldmpb,ldimb,indxsh,istab,g,ag,q,
     .    rmat,k1,k2,1,nbas,ib1,ld0r,[-1],ofa,h0(1+off1),of0i,ldhr,ldhb,h,
     .    ofhi)
C       Rotate il block, when h is composite
        if (ibla == 4) then
C         k1,k2 encompass iwaves for site ib1
          k1  = offH(1,1,nbas+1) + offH(2,1,ib1)+1
          k2  = offH(1,1,nbas+1) + offH(2,1,ib1+1)
C         ofa is row offset subtracted from k1<=k<=k2 when rotating
C         h0(k-ofa,P(*)) -> h(k,*).  ofa = first iwave col in h,
C         + size of l-block so as to preserve contents of l-block
          ofa = k1-1 - (offH(1,1,ib1+1)-offH(1,1,ib1))
          call prothr(opt,nl2,nbas,bas,ldmpb,ldimb,indxsh,istab,g,ag,
     .      q,rmat,k1,k2,1,nbas,ib1,ld0r,[-1],ofa,h0(1+off1),of0i,ldhr,
     .      ldhb,h,ofhi)
        endif
      else
        k1 = 1
        k2 = offH(ibla,1,nbas+1)
        call cplxdm(kcplx,ldha,ldhb,ld01,ld02,ld0r,of0i)
        call prothr(opt,nl2,nbas,bas,ldmpb,ldimb,indxsh,istab,g,ag,q,
     .    rmat,k1,k2,1,nbas,0,ld0r,[-1],0,h0,of0i,ldhr,ldhb,h,ofhi)
      endif
C     call yprm('R H R^-1',kcplx+2,h,ofhi,ldha,ldha,ldhb)

C --- Load h0 with rotated rows and columns of h; restore h ---
      if (lsaveh) then
        off1 = 0
        off2 = off1 + nlma*ldhb*2
C       Copy rows,columns of rotated h(iat) into h0
        call cplxdm(kcplx,nlma,ldhb,ld01,ld02,ld0r,of0i)
        call ymscop(0,nlma1,ldhb,ldhr,ld0r,ofa1,0,0,0,h,ofhi,
     .    h0(1+off1),of0i)
        call ymscop(0,nlma2,ldhb,ldhr,ld0r,ofa2,0,nlma1,0,h,ofhi,
     .    h0(1+off1),of0i)
        call cplxdm(kcplx,ldha,nlmb,ld01,ld02,ld0r,of0i)
        call ymscop(0,ldha,nlmb1,ldhr,ld0r,0,ofb1,0,0,h,ofhi,
     .    h0(1+off2),of0i)
        call ymscop(0,ldha,nlmb2,ldhr,ld0r,0,ofb2,0,nlmb1,h,ofhi,
     .    h0(1+off2),of0i)
        off1 = (nlma*ldhb + ldha*nlmb)*2
        off2 = off1 + nlma*ldhb*2
C       Restore rows,columns of h(iat)
        call cplxdm(kcplx,nlma,ldhb,ld01,ld02,ld0r,of0i)
        call ymscop(0,nlma1,ldhb,ld0r,ldhr,0,0,ofa1,0,h0(1+off1),
     .    of0i,h,ofhi)
        call ymscop(0,nlma2,ldhb,ld0r,ldhr,nlma1,0,ofa2,0,h0(1+off1),
     .    of0i,h,ofhi)
        call cplxdm(kcplx,ldha,nlmb,ld01,ld02,ld0r,of0i)
        call ymscop(0,ldha,nlmb1,ld0r,ldhr,0,0,0,ofb1,h0(1+off2),
     .    of0i,h,ofhi)
        call ymscop(0,ldha,nlmb2,ld0r,ldhr,0,nlmb1,0,ofb2,h0(1+off2),
     .    of0i,h,ofhi)
      endif

C --- Cleanup for pseudorotation (rotation + inversion) ---
      if (lroti) then
        istab(1) = -istab(1)
        call dscal(9,-1d0,g,1)
        call dscal(3,-1d0,q,1)
      endif

   99 continue
      kcplx = mod(opt/100,10)
      if (kcplx == 1) then
        call ztoyy(h,ldha,ldhb,ldha,ldhb,2,1)
        if (lsaveh) then
          call ztoyy(h0,nlma,ldhb,nlma,ldhb,2,1)
          off2 = nlma*ldhb*2
          call ztoyy(h0(1+off2),ldha,nlmb,ldha,nlmb,2,1)
        endif
      endif
      call tcx('roth')

      roth = protht(opt,nbas,iat,offH,indxsh,istab,ldha,ldhb,h0,h)

C     call zprm('ending h',2,h,ldha,ldha,ldhb)
      end
      subroutine pvrotd(opt,offH,nbas,iat,ibla,iblb,ldmpa,ldima,ldmpb,
     .  ldimb,ofa1,ofa2,nlma,nlma1,nlma2,ofb1,ofb2,nlmb,nlmb1,nlmb2)
C- Sets up dimensions for roth.
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :see roth
Ci   offH  :Offsets to hamiltonian matrix (makidx.f)
Ci   nbas  :size of basis
Ci   iat   :see roth
Co Outputs
Co   ibla  :index to offH that define which block of hamiltonian
Co         :is to be rotated in row dimension.
Co         :ibla=(1,2,3) -> (lower,intermediate,higher) block
Co         :ibla=4 -> lower+intermediate blocks
Co   iblb  :like ibla, but for column dimension of hamiltonian
Co   ldmpa :last orbital in block prior to current block, row dimension
Co   ldima :last orbital in current block, row dimension
Co   ldmpb :last orbital in block prior to current block, col dimension
Co   ldimb :last orbital in current block, row dimension
Co
Co   ... the next quantities are defined only for 1s digit opt>=2
Co   (when only rows, cols associated with site iat are rotated).
Co   ofa1  :offset to lower block belonging to ib (row index)
Co   ofa2  :(if compound l+i block) offset to corresponding intermed. block
Co   nlma  :number of orbitals connected to site ib (row index)
Co   nlma1 :number of orbitals in lower block (row index)
Co   nlma2 :number of orbitals in upper block (row index)
Co   ofb1  :offset to lower block belonging to ib (col index)
Co   ofb2  :(if compound l+i block) offset to corresponding intermed. block
Co   nlmb  :number of orbitals connected to site ib (col index)
Co   nlmb1 :number of orbitals in lower block (col index)
Co   nlmb2 :number of orbitals in upper block (col index)
Cu Updates
Cu   05 Feb 01
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer opt,nbas,iat,ibla,iblb,ldmpa,ldima,ldmpb,ldimb,nkap0,n0H,
     .  ofa1,ofa2,nlma,nlma1,nlma2,ofb1,ofb2,nlmb,nlmb1,nlmb2
      parameter (nkap0=4,n0H=5)
      integer offH(n0H,nkap0,*)
C ... Local parameters
      integer isum

      ibla = max(mod(opt/1000,10),1)
      iblb = max(mod(opt/10000,10),1)
      if (iblb > 1 .and. (ibla /= iblb .or. iblb > 3))
     .  call rxi('roth: inconsistent opt,',opt)

C ... Row dimensions and offsets
      if (ibla < 4) then
        ldima = isum(ibla,offH(1,1,nbas+1),1)
        ldmpa = ldima - offH(ibla,1,nbas+1)
        ldmpb = 0
        ldimb = offH(1,1,nbas+1)
        if (iblb > 1) then
          ldmpb = ldmpa
          ldimb = ldima
        endif

C       Dimensions and offsets for site iat
C       If mod(opt,10)<=1, these quantities are undefined
        if (mod(opt,10) > 1) then
          ofa1 = offH(ibla,1,iat)
          ofa2 = ofa1
          nlma = offH(ibla,1,iat+1) - offH(ibla,1,iat)
          nlma1 = nlma
          nlma2 = 0

          ofb1 = offH(iblb,1,iat)
          ofb2 = ofb1
          nlmb = offH(iblb,1,iat+1) - offH(iblb,1,iat)
          nlmb1 = nlmb
          nlmb2 = 0
        endif

C ... Dimensions and offsets when h a composite of downfolding blocks
      elseif (ibla == 4) then
        ldmpa = 0
        ldima = offH(1,1,nbas+1) + offH(2,1,nbas+1)
        ldmpb = 0
        ldimb = ldima

C       Dimensions and offsets for site iat
        if (mod(opt,10) > 1) then
          ofa1 = offH(1,1,iat)  ! offset to l-block at ib, full hamiltonian
          ofa2 = offH(1,1,nbas+1) + offH(2,1,iat) ! offset to i-block at ib, full hamiltonian
          nlma = offH(ibla,1,iat+1) - offH(ibla,1,iat) ! size of l+i block at ib
          nlma1 = offH(1,1,iat+1) - offH(1,1,iat) ! size of lower block at ib
          nlma2 = nlma - nlma1  ! size of intermediate block at ib

          ofb1 = ofa1
          ofb2 = ofa2
          nlmb = nlma
          nlmb1 = nlma1
          nlmb2 = nlma2
        endif

      else
        call rx1('roth: illegal value for opt (%i)',opt)
      endif
      end
      subroutine prothl(opt,ldr,nbas,bas,ldmpa,ldima,indxsh,istab,g,ag,
     .  q,rmat,k1,k2,ib1,ib2,iskip,ld0a,off0,ofb,h0,of0i,ldha,ldhb,h,ofhi)
C- Rotate h0 <- rmat h for some block of h
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :1s digit:
Ci         :0  rotate entire h
Ci         :1  rotate site-diagonal h
Ci         :5  same as case 1, but h contains only blocks diagonal in l
Ci         :10s digit
Ci         :0  h is in real spherical harmonics
Ci         :Add 2 if h is in true (complex) spherical harmonics
Ci         :   See rcub2sph.f for definition of true spherical harmonics.
Ci         :Add 4 if phase convention phi = q * [(g R_j + a) - R_i]
Ci         :   Should be scaled by -1
Ci   ldr   :leading dimension of rmat
Ci   nbas  :size of basis
Ci   bas   :basis vectors, in units of alat
Ci   ldmpa :offset to first orbital in downfolding block of h to rotate
Ci   ldima :points to last orbital in  downfolding block of h to rotate
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   istab :table of site permutations for each group op (mksym.f,symtbl.f)
Ci         :istab(1) = 0 => no permutation of sites
Ci   g     :point group operation
Ci   ag    :translation part of space group
Ci   q     :wave number corresponding to Bloch sum
Ci   rmat  :YLM rotation matrix for this group operation
Ci         :If rmat is complex, it must be in the kcplx=0 mode, i.e.
Ci         :double complex with imaginary following real
Ci   k1    :columns h(*,k1..k2) are rotated into h0(P(*)),k1..k2) ;
Ci   k2    :see Remarks
Ci   ib1   :subblocks of h for sites ib1..ib2 are rotated; see Remarks
Ci   ib2   : -//-
Ci   iskip :site ib = iskip is excluded from rotation
Ci   ld0a  :leading dimension of h0
Ci   off0  :Special indexing for h0 when off0(1)>=0 .
Ci         :(For default, use off0(1) = -1; off0 is then not used).
Ci         :Special indexing is used when the caller needs to specify
Ci         :the row ordering of h0, e.g. when a subblock of h is needed,
Ci         :and the mapping of the subblock is not into contiguous rows
Ci         :as occurs when concatenating orbitals of different downfolding
Ci         :subblocks.
Ci   ofb   :subtracted from the column offset when storing into h0.
Ci         :When only a portion of h is rotated (columns belonging to
Ci         :only one site), it conserves memory to store into first
Ci         :columns of h0, rather than the columns corresponding to site
Ci         :Then ofb = k1-1.
Ci   of0i  :offset to imaginary part of h0
Ci   ldha  :leading dimension of h
Ci   ldhb  :second dimension of h (not used)
Ci   h     :matrix to be rotated
Ci   ofhi  :offset to imaginary part of h
Co Outputs
Co   h0    :contains rotated subblock of h
Cl Local variables
Cl   korb  :l orbital column index.  Rotations are for blocks of
Cl          contiguous columns.  korb is an index that marks current
Cl          group of columns are being rotated.
Cr Remarks
Cr   This routine rotates rmat h, storing result into h0 for a block of
Cr   rows corresponding to sites ib1..ib2 and columns k1..k2.  Each site
Cr   ib is permuted into site istab(ib), thus permuting the rows of h
Cr   in addition to rotating by YLM(g).  Finally if matrix h corresponds
Cr   to the periodic part of Bloch hamiltonian, there is a phase factor
Cr     exp(i q (g R_j + a) - R_i)) where  (g R_j + ag) - R_i
Cr   is a lattice vector between sites i and j, and q is the wave number
Cr
Cr   NB: in the site-diagonal case (1s digit opt=1), columns k1 and k2
Cr   are not fixed but depend on the rows.  In that case k1 and k2
Cr   are set internally, and the input values are ignored.
Cr
Cr   prothl can handle only a single downfolding block in the row
Cr   dimension, because it assumes the orbitals corresponding to one
Cr   site lie in contiguous rows.  The columns are not permuted and
Cr   can be of any number.
Cu Updates
Cu   28 Apr 17 Use of istab can be suppressed; ag,g,bas not accessed if q=0
Cu   17 Jun 02 Added option for phase convention (10s digit opt)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer opt,ld0a,ldha,ldhb,ldr,nbas,istab(nbas),off0(*),
     .  indxsh(*),ldmpa,ldima,iskip,k1,k2,ofb,of0i,ofhi
      double precision bas(3,*),g(3,3),ag(3),q(3),
     .  rmat(ldr,ldr),h0(ld0a,*),h(ldha,ldhb)
C ... Local parameters
C     logical :: debug=.false.
      logical nphase,lperm
      integer n0,nkap0,i,ib1,ib2,ibas,korb,j,jbas,k,li,ncut,
     .  nk,nlm,nlmi,norb,offh00,offr,opt0,opt1,opt4,iorb
      parameter (n0=10,nkap0=4)
      integer offi,offj,ltab(n0*nkap0),ktab(n0*nkap0),
     .  offli(n0*nkap0),offlj(n0*nkap0)
      double precision twopi,sp,cosP,sinP,tbas(3),xx
      procedure(real(8)) :: ddot
      parameter (ncut=9)

      if (k2 < k1 .and. mod(opt,10) /= 1) return

C     call tcn('rothl')

      opt0 = mod(opt,10)
      opt1 = mod(mod(opt/10,10),4)
      opt4 = mod(opt/10000,10)
      nphase = mod(opt/10,10) >= 4
      twopi = 8d0*datan(1d0)
      lperm = istab(1) /= 0

      offh00 = 0
      do  ibas = ib1, ib2
        jbas = ibas
        if (lperm) then
          jbas = istab(ibas)
        endif

C   ... Skip over diagonal site second pass
        if (ibas == iskip) cycle

C   ... Setup for multiplication by phase
C       tbas = (g R_j + a) - R_i; should be a lattice vector
        if (ddot(3,q,1,q,1) /= 0) then
          do  i = 1, 3
            tbas(i) = ag(i) - bas(i,ibas)
            do  j = 1, 3
              tbas(i) = tbas(i) + g(i,j)*bas(j,jbas)
            enddo
          enddo
C         exp(i q tbas)
          sp = twopi*(tbas(1)*q(1) + tbas(2)*q(2) + tbas(3)*q(3))
          if (nphase) sp = -sp
          cosP = dcos(sp)
          sinP = dsin(sp)
        else
          sp = 0
        endif

C   ... Row, column offsets to h, h0
        korb = 0
        call orbl(ibas,ldmpa,ldima,indxsh,norb,ltab,ktab,xx,offli,nlm)
        call orbl(jbas,ldmpa,ldima,indxsh,norb,ltab,ktab,xx,offlj,nlm)
   12   continue
        if (nlm == 0) cycle

C       Calculate k1,k2 for site-diagonal part of h (opt0 = 1)
        if (opt0 == 1) then
C         We cannot assume orbitals are contiguous.  Count the number of
C         contiguous orbitals.  korb will be index to last contiguous
C         orbital, and k1..k2 will encompass the offsets in h.  Must
C         be repeated until all blocks are exhausted.
          korb = korb+1
          k1 = offlj(korb)+1
          k2 = offlj(korb)
   16     continue
            k2 = k2 + 2*ltab(korb)+1
            if (korb < norb) then
C             If contiguous, increment korb and increase k2
              if (offlj(korb+1) == k2) then
                korb = korb + 1
                goto 16
              endif
            endif
C         For ii (or hh) case, h, h0 consists only of ii (or hh) block
          if (opt4 > 1) then
            k1 = k1-ldmpa
            k2 = k2-ldmpa
          endif
        endif

C   ... For each l, do rmat h(jbas,lm) into h0(ibas).  Only rows mix
        offi = offh00           ! row offset to h0(ibas)
        iorb = 0                ! index to current l
C       When h is in l-diagonal form, column offset starts at zero
C       and increments by 2*l+1 for each orbital in the list.
        if (opt0 == 5) k2 = 0
        do  iorb = 1, norb
          li = ltab(iorb)
          if (off0(1) < 0) then
            offi = offli(iorb) - ldmpa
          else
            offi = off0(iorb)
          endif
          offj =  offlj(iorb) - ldmpa ! row offset to h(jbas,li)
          nlmi = 2*li + 1       ! number of m's for this l
          offr = li**2          ! offset to li'th block in rmat

C    ... Debugging ... copy rmat into h
C        do  j = 1, nlmi
C        do  i = 1, nlmi
C          h0(i+offi,j+offj) = rmat(i+offr,j+offr)
C          h0(i+offi+of0i,j+offj) = 0
C        enddo
C        enddo
C        goto 20

C         When h is in l-diagonal form,  see comment above.
          if (opt0 == 5) then
            k1 = k2+1
            k2 = k2+nlmi
          endif
          nk = k2-k1+1
          if (opt1 > 1) then  ! spherical harmonics: complex rmat
C            if (debug) print *, offi,offj
C            if (debug) call yprm('rl',0+2,rmat(1+offr,1+offr),ldr**2,ldr,nlmi,nlmi)
C            if (debug) call yprm('h',0+2,h(offj+1,k1),ofhi,ldha,nlmi,nk)
C
C            h0(offi+1,1:nk) = 999
C            call dcopy(nk,h0(offi+1,1:nk),1,h0(offi+of0i+1,1:nk),1)

            call yygemm('N','N',nlmi,nk,nlmi,1d0,rmat(offr+1,offr+1),
     .        rmat(offr+1,offr+1+ldr),ldr,
     .        h(offj+1,k1),h(offj+ofhi+1,k1),ldha,0d0,
     .        h0(offi+1,k1-ofb),h0(offi+of0i+1,k1-ofb),ld0a)
C            if (debug) call yprm('h0',0+2,h0(offi+1,k1-ofb),of0i,ld0a,nlmi,nk)
          elseif (nlmi > ncut) then
            call dgemm('N','N',nlmi,nk,nlmi,1d0,rmat(offr+1,offr+1),ldr,
     .        h(offj+1,k1),ldha,0d0,h0(offi+1,k1-ofb),ld0a)
            call dgemm('N','N',nlmi,nk,nlmi,1d0,rmat(offr+1,offr+1),ldr,
     .        h(offj+ofhi+1,k1),ldha,0d0,h0(offi+of0i+1,k1-ofb),ld0a)
          else
            do  k = k1, k2
            do  i = 1, nlmi
              h0(i+offi,k-ofb) = 0
              h0(i+offi+of0i,k-ofb) = 0
              do  j = 1, nlmi
                h0(i+offi,k-ofb) =      h0(i+offi,k-ofb) +
     .                                  rmat(i+offr,j+offr)*h(j+offj,k)
                h0(i+of0i+offi,k-ofb) = h0(i+of0i+offi,k-ofb) +
     .                                  rmat(i+offr,j+offr)*h(j+offj+ofhi,k)
              enddo
              enddo
            enddo
          endif

C     --- Multiply rows by phase ---
          if (dabs(sp) > 1d-8) then
            do  i = 1, nlmi
C             print 456, offi+i,cosP,sinP
C  456        format(' phase',i4,2f15.10)
              call yscal(nk,cosP,sinP,h0(offi+i,k1-ofb),
     .          h0(offi+of0i+i,k1-ofb),ld0a)
            enddo
          endif
        enddo
C       Case column blocks are not contiguous
        if (opt0 == 1 .and. korb < norb) goto 12
C       Increment offh00 by the number of orbitals in this site
        offh00 = offh00 + nlm
      enddo

C     call yprm('h0',kcplx+2,h0,of0i,ld0a,ld0a,k2-k1+1)

C     call tcx('rothl')

      end

      subroutine prothr(opt,ldr,nbas,bas,ldmpb,ldimb,indxsh,istab,g,ag,
     .  q,rmat,k1,k2,ib1,ib2,iskip,ld0a,off0,ofa,h0,of0i,ldha,ldhb,h,ofhi)
C- Rotate h <- h0 rmat+ for some block of h0
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :1s digit:
Ci         :0  rotate entire h
Ci         :1  rotate site-diagonal h
Ci         :5  same as case 1, but h contains only blocks diagonal in l
Ci         :10s digit
Ci         :0  h is in real spherical harmonics
Ci         :Add 2 if h is in true (complex) spherical harmonics
Ci         :   See rcub2sph.f for definition of true spherical harmonics.
Ci         :Add 4 if phase convention phi = q * [(g R_j + a) - R_i]
Ci         :   Should be scaled by -1
Ci   ldr   :leading dimension of rmat
Ci   nbas  :size of basis
Ci   bas   :basis vectors, in units of alat
Ci   ldmpb :offset to first orbital in downfolding block of h to rotate
Ci   ldimb :points to last orbital in  downfolding block of h to rotate
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   istab :table of site permutations for each group op (mksym.f,symtbl.f)
Ci   g     :point group operation
Ci   ag    :translation part of space group
Ci   q     :wave number corresponding to Bloch sum
Ci   rmat  :YLM rotation matrix for this group operation
Ci         :If rmat is complex, it must be in the kcplx=0 mode, i.e.
Ci         :double complex with imaginary following real
Ci   k1    :rows h0(k1..k2,*) are rotated into h(k1..k2,P(*))
Ci   k2    :see Remarks
Ci   ib1   :subblocks of h for sites ib1..ib2 are rotated; see Remarks
Ci   ib2   : -//-
Ci   iskip :site ib = iskip is excluded from rotation
Ci   ld0a  :leading dimension of h0
Ci   off0  :Special indexing for h0 when off0(1)>=0 .
Ci         :(For default, use off0(1) = -1; off0 is then not used).
Ci         :Special indexing is used when the caller needs to specify
Ci         :the row ordering of h0, e.g. when a subblock of h is needed,
Ci         :and the mapping of the subblock is not into contiguous rows
Ci         :as occurs when concatenating orbitals of different downfolding
Ci         :subblocks.
Ci   cth   :1 assemble the columns from the different l channels from
Ci         :contiguous columns in h0, starting from the first column.
Ci         :Needed if h is partitioned into lower+higher blocks, so that
Ci         :consecutive l channels are not in contiguous columns.
Ci   ofa   :subtracted from the row offset when storing into h0.
Ci         :When only a portion of h is rotated (rows belonging to
Ci         :only one site), it conserves memory to store into first
Ci         :rows of h0, rather than the rows corresponding to site
Ci         :Then ofb = k1-1.
Ci   h0    :unrotated subblock of h
Ci   of0i  :offset to imaginary part of h0
Ci   ldha  :leading dimension of h
Ci   ldhb  :second dimension of h (not needed)
Co Outputs
Co   h     :rotated subblock of h0; see Remarks
Cl Local variables
Cl   korb  :l orbital row index.  Rotations are for blocks of
Cl          contiguous rows.  korb is an index that marks current
Cl          group of rows are being rotated.
Cr Remarks
Cr   This routine rotates h0 rmat, storing result into h for a block of
Cr   columns corresponding to sites ib1..ib2 and rows k1..k2.  See
Cr   prothl for description of a similar rotation.
Cu Updates
Cu   28 Apr 17 Use of istab can be suppressed; ag,g,bas not accessed if q=0
Cu   17 Jun 02 Added option for phase convention (10s digit opt)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer opt,ld0a,ldha,ldhb,ldr,nbas,istab(nbas),indxsh(*),ldmpb,
     .  ldimb,iskip,k1,k2,ofa,of0i,ofhi,off0(1)
      double precision bas(3,*),g(3,3),ag(3),q(3),rmat(ldr,ldr),
     .  h0(ld0a,*),h(ldha,ldhb)
C ... Local parameters
      logical nphase,lperm
      integer i,ib1,ib2,ibas,korb,j,jbas,k,li,ncut,nk,nlm,nlmi,norb,
     .  offi,offj,offr,opt0,opt1,opt4,iorb,n0,nkap0
      parameter (n0=10,nkap0=4)
      integer ltab(n0*nkap0),ktab(n0*nkap0),offli(n0*nkap0),
     .  offlj(n0*nkap0)
      double precision twopi,sp,cosP,sinP,tbas(3),xx
      procedure(real(8)) :: ddot
      parameter (ncut=9)

      if (k2 < k1) return
C     call tcn('rothr')
      opt0 = mod(opt,10)
      opt1 = mod(mod(opt/10,10),4)
      nphase = mod(opt/10,10) >= 4
      opt4 = mod(opt/10000,10)
      twopi = 8d0*datan(1d0)
      lperm = istab(1) /= 0

      do  ibas = ib1, ib2
        jbas = ibas
        if (lperm) then
          jbas = istab(ibas)
        endif

C   ... Skip over diagonal site second pass
        if (ibas == iskip) cycle

C   ... Setup for multiplication by phase
C       tbas = (g R_j + a) - R_i; should be a lattice vector
        if (ddot(3,q,1,q,1) /= 0) then
          do  i = 1, 3
            tbas(i) = ag(i) - bas(i,ibas)
            do  j = 1, 3
              tbas(i) = tbas(i)+g(i,j)*bas(j,jbas)
            enddo
          enddo
C         exp(i q tbas)
          sp = twopi*(tbas(1)*q(1) + tbas(2)*q(2) + tbas(3)*q(3))
          if (nphase) sp = -sp
          cosP = dcos(sp)
          sinP = dsin(sp)
        else
          sp = 0
        endif

C   ... Offset to first orbital in ibas,jbas.  Uses offi,nlm
        korb = 0
        call orbl(ibas,ldmpb,ldimb,indxsh,norb,ltab,ktab,xx,offli,nlm)
        call orbl(jbas,ldmpb,ldimb,indxsh,norb,ltab,ktab,xx,offlj,nlm)
   12   continue
        if (nlm == 0) cycle

C       Calculate k1,k2 for diagonal part of h (opt0 = 1)
        if (opt0 == 1) then
C         We cannot assume orbitals are contiguous.  Count the number of
C         contiguous orbitals.  korb will be index to last contiguous
C         orbital, and k1..k2 will encompass the offsets in h.  Must
C         be repeated until all blocks are exhausted.
          korb = korb+1
          k1 = offli(korb)+1
          k2 = offli(korb)
   16     continue
            k2 = k2 + 2*ltab(korb)+1
            if (korb < norb) then
C             If contiguous, increment korb and increase k2
              if (offli(korb+1) == k2) then
                korb = korb + 1
                goto 16
              endif
            endif
C         For ii (or hh) case, h, h0 consist only of ii (or hh) block
          if (opt4 > 1) then
            k1 = k1-ldmpb
            k2 = k2-ldmpb
          endif
        endif

C   ... For each l, do h0(jbas,lm) rmat into h(ibas). Only columns mix
C       When h is in l-diagonal form,  Assume k1,k2 assume role of
C       offi and offi starts at zero and increments by 2*l+1 for
C       each orbital in the list.
        offi = 0                ! column offset to h(ibas)
        do  iorb = 1, norb
          li = ltab(iorb)
          nlmi = 2*li + 1
          offr = li**2
C         When h is in l-diagonal form,  see comment above.
          if (opt0 == 5) then
            k1 = 1    + offli(iorb) - ldmpb
            k2 = nlmi + offli(iorb) - ldmpb
            offj = offi
          else
            offi = offli(iorb) - ldmpb
            if (off0(1) < 0) then
              offj =  offlj(iorb) - ldmpb ! column offset to h0(jbas)
            else
              offj =  off0(iorb)
            endif
          endif
          nk = k2-k1+1
          if (opt1 > 1) then  ! spherical harmonics: complex rmat
            call yygemm('N','C',nk,nlmi,nlmi,1d0,h0(k1-ofa,offj+1),
     .        h0(k1+of0i-ofa,offj+1),ld0a,
     .        rmat(offr+1,offr+1),rmat(offr+1,offr+1+ldr),ldr,0d0,
     .        h(k1,offi+1),h(k1+ofhi,offi+1),ldha)
          elseif (nlmi > ncut) then
            call dgemm('N','T',nk,nlmi,nlmi,1d0,h0(k1-ofa,offj+1),
     .        ld0a,rmat(offr+1,offr+1),ldr,0d0,h(k1,offi+1),ldha)
            call dgemm('N','T',nk,nlmi,nlmi,1d0,h0(k1+of0i-ofa,offj+1),
     .        ld0a,rmat(offr+1,offr+1),ldr,0d0,h(k1+ofhi,offi+1),ldha)
          else
            forall(i=1:nlmi, k=k1:k2)
              h(k,i+offi) = 0
              h(k+ofhi,i+offi) = 0
            end forall
            do  i = 1, nlmi
              do  j = 1, nlmi
                do  k = k1, k2
                  h(k,i+offi) = h(k,i+offi) +
     .                          h0(k-ofa,j+offj)*rmat(i+offr,j+offr)
                  h(k+ofhi,i+offi) = h(k+ofhi,i+offi) +
     .                          h0(k+of0i-ofa,j+offj)*rmat(i+offr,j+offr)
                enddo
              enddo
            enddo
          endif

C     --- Multiply columns by phase ---
          if (dabs(sp) > 1d-8) then
            do  i = 1, nlmi
              call yscal(nk,cosP,-sinP,h(k1,offi+i),h(k1+ofhi,offi+i),1)
            enddo
          endif

C         Increment offi in diagonal case
          offi = offi + nlmi
        enddo
C       Case column blocks are not contiguous
        if (opt0 == 1 .and. korb < norb) goto 12
      enddo

C     call tcx('rothr')

      end
      integer function protht(opt,nbas,iat,offH,iprmb,istab,ldha,ldhb,h0,h)
C- Make appropriate transposes for pseudorotation, depending on mode
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :see roth
Ci   istab :table of site permutations for each group op (mksym.f,symtbl.f)
Ci   h0    :see roth
Ci   h     :small h = cnu-enu+sqrdel*S^beta*sqrdel
Ci   iat   :if opt=0 or opt=1, not used (1s digit)
Ci         :if opt=3, h(1..nbas,iat) h(iat,1..nbas) are rotated
Ci   offH  :Offsets to hamiltonian matrix (makidx.f)
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   nbas  :size of basis
Ci   ldha  :leading dimension of h
Ci   ldhb  :second dimension of h
Co Outputs
Co   h0    :some portion or all of h0 is transposed, depending on opt
Co   h     :some portion or all of h is transposed, depending on opt
Co   protht:returns  0 if no transpose required
Co         :returns -1 transpose required but protht did not make it
Co         :returns  1 protht took transpose
Cl Local variables
Ci  nrh    :spacing between adjacent row elements
Ci  nch    :spacing between adjacent column elements
Cr Remarks
Cr   This subroutine is intended as a 'cleanup' to roth to handle cases
Cr   when the inversion operation is not part of the space group but
Cr   is included anyway because of time reversal symmetry.
Cr   (This is flagged when istab(1)<0).  In such a case, roth generates
Cr   the transpose of h.  For the cases it is able to handle, protht
Cr   overwrites h (or part of it) with its transpose.  For it to do so,
Cr   protht requires rows, columns to correspond to the same orbitals.
Cr   it cannot handle opt=2000, 2003, 2004.  protht could be made to
Cr   handle opt=4003 but it hasn't been implemented.
Cb Bugs
Cb   kcplx<>0 not checked for mode 5 (kcplx=1 will not work).
Cu Updates
Cu   05 Feb 01 Created and checked using troth.
C ----------------------------------------------------------------------
      implicit none
      integer opt,ldha,ldhb,istab(*),iat,nbas,iprmb(*)
      integer nkap0,n0H
      parameter (nkap0=4,n0H=5)
      integer offH(n0H,nkap0,*)
      double precision h(ldha,*),h0(*)
C
      integer nrh,nch,ofhi,ldh1,ldh2,w,opt0,kcplx,ibla,iblb,
     .  ldmpa,ldmpb,ofa1,ofa2,nlma,nlma1,nlma2,ofb1,ofb2,nlmb,nlmb1,
     .  nlmb2,ldima,ldimb,off2,ndim

      protht = 0
      if (istab(1) > 0) return
      protht = -1

C ... Get some dimensioning parameters
      opt0 = mod(opt,10)
      call pvrotd(opt,offH,nbas,iat,ibla,iblb,ldmpa,ldima,ldmpb,
     .  ldimb,ofa1,ofa2,nlma,nlma1,nlma2,ofb1,ofb2,nlmb,nlmb1,nlmb2)
C     protht requires rows, columns to correspond to same orbitals
      if (ibla /= 4 .and. iblb /= ibla) return
C     if (ldima /= ldimb) return

      kcplx = mod(opt/100,10)
      call cplxdm(kcplx,ldha,ldhb,ldh1,ldh2,nch,ofhi)
      nrh = 1
      if (kcplx == 1) nrh = 2

C ... rotate entire h
      if (opt0 == 0) then
        ndim = ldima-ldmpa
        call ymtrns(0,h,nch,nrh,ofhi,w,w,w,w,1,ndim,1,ndim)
        protht = 1
C ... site-diagonal h
      elseif (opt0 == 1) then
       call pvrots(opt0,nbas,ldha,ldhb,iprmb,ldmpa,ldima,nch,nrh,ofhi,h)
       protht = 0
C ... rows and columns connected with one site
      elseif (opt0 == 3) then
        if (nlma/=nlmb .or. ldima/=ldimb .or. ofa1/=ofb1) return
        if (nlma2 /= 0) return
        ndim = ldima-ldmpa
        call ymtrns(0,h,nch,nrh,ofhi,w,w,w,w,ofa1+1,ofa1+nlma,1,ndim)
        protht = 1
      elseif (opt0 == 4) then
        if (nlma/=nlmb .or. ldima/=ldimb .or. ldmpa/=ldmpb) return
        ldima = ldima - ldmpa
        ldimb = ldimb - ldmpb
        off2 = nlma*ldimb*2
        if (kcplx /= 0) then
          call ztoyy(h0,nlma,ldimb,nlma,ldimb,kcplx,0)
          call ztoyy(h0(1+off2),ldima,nlmb,ldima,nlmb,kcplx,0)
        endif
        call pvrott(nlma,ldima,h0,h0(1+off2))
        if (kcplx /= 0) then
          call ztoyy(h0,nlma,ldimb,nlma,ldimb,0,kcplx)
          call ztoyy(h0(1+off2),ldima,nlmb,ldima,nlmb,0,kcplx)
        endif
        protht = 1
C ... hh block
      elseif (opt0 == 5) then
       if (ibla /= 3 .or. iblb /= 3) return
       call pvrots(opt0,nbas,ldha,ldhb,iprmb,ldmpa,ldima,nch,nrh,ofhi,h)
       protht = 0

      else
        return
      endif

      end
      subroutine pvrott(nlm,ndim,sr,sc)
      implicit none
      integer nlm,ndim
      double precision sr(nlm,ndim,2),sc(ndim,nlm,2)
      integer i,j
      double precision wkr,wki

      do  i = 1, nlm
      do  j = 1, ndim
        wkr = sr(i,j,1)
        wki = sr(i,j,2)
        sr(i,j,1) = sc(j,i,1)
        sr(i,j,2) = sc(j,i,2)
        sc(j,i,1) = wkr
        sc(j,i,2) = wki
      enddo
      enddo
      end

      subroutine pvrots(opt0,nbas,ldha,ldhb,iprmb,ldmpa,ldima,nch,nrh,ofhi,h)
      implicit none
      integer nbas,opt0,ldha,ldhb,iprmb(*),ldmpa,ldima,nch,nrh,ofhi
      double precision h(ldha,ldhb)
C Local variables
      integer n0,nkap0,jb,k1,k2,norb,nlm,korb
      parameter (n0=10,nkap0=4)
      integer ltab(n0*nkap0),ktab(n0*nkap0),offlj(n0*nkap0)
      double precision xx

      do  jb = 1, nbas
      korb = 0
      call orbl(jb,ldmpa,ldima,iprmb,norb,ltab,ktab,xx,offlj,nlm)
      if (nlm /= 0) then
   12   continue
C       We cannot assume orbitals are contiguous.  Count the number of
C       contiguous orbitals.  korb will be index to last contiguous
C       orbital, and k1..k2 will encompass the offsets in h.  Must
C       be repeated until all blocks are exhausted.
        korb = korb+1
        k1 = offlj(korb)+1
        k2 = offlj(korb)
   16   continue
          k2 = k2 + 2*ltab(korb)+1
C         Look for contiguous blocks (not higher blocks)
          if (korb < norb .and. opt0 /= 5) then
            if (offlj(korb+1) == k2) then
              korb = korb + 1
              goto 16
            endif
          endif

C       For higher blocks h consists of hh block;
        if (opt0 == 5) then
C         k1-ldmpa = offset rel to start of hh block
C         k1-offlj(1) = offset for start of this atom
          k1 = k1-offlj(1)
          k2 = k2-offlj(1)
          call ymtrns(0,h(k1-ldmpa+offlj(1),k1),nch,nrh,ofhi,xx,xx,xx,
     .      xx,k1,k2,k1,k2)
        else
          k1 = k1-ldmpa
          k2 = k2-ldmpa
          call ymtrns(0,h,nch,nrh,ofhi,xx,xx,xx,xx,k1,k2,k1,k2)
        endif

        if (korb < norb) goto 12
      endif
      enddo

      end
