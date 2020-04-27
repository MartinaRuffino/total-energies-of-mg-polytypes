      subroutine pgbevl(s_ham,s_pot,s_str,s_site,plat,isp,nspc,ipl,npl,pgplp,qp,
     .  mode,modegF,rvec,r,sbulk,gbii,ndg,offg,gii,lgii)
C- Generate evals and evecs of bulk GF, and optionally bulk GF
C ----------------------------------------------------------------
Cio Structures
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  neula ldham lgen3 lncol lham offH
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb eula
Cio    Passed to:  plhamnc plham plhamso
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:palp pfnc
Cio    Passed to:  plhamnc plham
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iax npr s
Cio    Passed to:  plhamnc plham
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  norb
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pfr
Cio    Passed to:  plhamnc plham plhamso
Ci Inputs:
Ci   plat,nl,nbas: program inputs
Ci   nbasp :number of atoms in the padded basis
Ci          nbasp = nbas + nbas(left bulk) + nbas(right bulk)
Ci   ipl   :index to current principal layer
Ci          -1 => strux for left bulk PL; npl => strux for right bulk
Ci   npl   :number of principal layers (pgfset.f)
Ci   pgplp :index and dimensioning information for each PL (pgfset.f)
Ci   qp    :k-point
Ci   mode: : 1s digit
Ci             0, Use strx (s00,s0L,s0L+)
Ci             1, Use strx (s00,s0R+,s0R)
Ci           10s digit:
Ci             0, r,rvec are arrays passed by the calling program
Ci             1, internal arrays used for r,rvec
Ci           100s digit:
Ci             0, sbulk unused
Ci             1, return sbulk
Ci   ndg,offg: dimension of g and offset to diagonal part of g
Ci           Normally they are pgplp(3,ipl) and pgplp(6,ipl)
Ci           but there is no requirement they be so.
Ci   modeGF: -1 do not call pgbulk, make no GF
Ci           otherwise, call pgbulk with modeGF; for documentation of
Ci           modeGF,documentation for 'mode' in pbgulk.
Cl Local and global variables
Co Outputs:
Co   r,rvec:eigenvalues and eigenvectors of bulk GF, layer ipl
Co          provided 10s digit of mode is set.  Else, left untouched.
Co   gii   :Bulk GF as calculated according to modeGF (see pgbulk)
Co   lgii  :attributes of the gii generated; see lgupac for conventions
Co   sbulk :s00,s0L,s0R that generate evals and evecs
Cr
Cr   pghof doesn't use lmx; should be struck from arg list
Cu Updates
Cu   17 May 16 Extend to noncollinear case
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
C ----------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer isp,nspc,npl,mode,modeGF,ndg,offg,pgplp(6,-1:npl),lgii
      double precision plat(3,3),qp(3),rvec(*),r(*),gii(*),gbii(*),
     .  sbulk(*)
C ... For structures
!      include 'structures.h'
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_str)::   s_str
      type(str_site)::  s_site(*)
C ... Dynamically allocated local arrays
      complex(8), allocatable :: rvecl(:),rl(:),sll(:),s00(:)
      real(8), allocatable :: wk(:)
C ... Local parameters
      integer,parameter :: PRTG=80
      integer ipr,ipl,ld0,ldl,ldr,ld,opt,ld0x,ldlx,ldrx,ldx
      procedure(integer) :: lgunit

C --- Setup ---
      call tcn('pgbevl')
C     call wkprnt(1)
      call getpr(ipr)
      if (ipr >= 60)
     .  call awrit3(' pgbevl: ipl=%i  mode=%i  modeGF=%i',' ',80,
     .  lgunit(1),ipl,mode,modeGF)
      call rxx(ipl < -1 .or. ipl > npl,'pgbevl: ipl out of range')
      ldl = pgplp(4,max(ipl-1,0))
      ld0 = pgplp(4,max(min(ipl,npl-1),0))
      ldr = pgplp(4,min(ipl+1,npl-1))
      ld = ldl+ld0+ldr
      ld0x = ld0*nspc
      ldlx = ldl*nspc
      ldrx = ldr*nspc
      ldx = ld*nspc

C --- Strux connecting ipl to adjacent layers ---
      allocate(s00(ldx*ld0x),sll(ldx*ld0x))
      call plhamnc(s_ham,s_pot,s_str,s_site,0,plat,isp,nspc,0,ipl,npl,ld0,ld,pgplp,qp,sll,s00)
      deallocate(sll)
C ... Extract bulk strux into rvec from s00, put into rvec
      if (mod(mode,10) == 0) then
        opt = 1211
      else
        opt = 1121
      endif
      if (mod(mode/10,10) == 1) then
        allocate(rl(2*ld0x),rvecl(4*ld0x**2))
        call pgcops(ld0x,ldlx,ldrx,ld0x,ld0x,opt,s00,rvecl)
C       call yprm('pgbevl: bulk strx',2,rvecl,ld0x**2*3,ld0x,ld0x,ld0x*3)
      else
        call pgcops(ld0x,ldlx,ldrx,ld0x,ld0x,opt,s00,rvec)
C       call yprm('pgbevl: bulk strux',2,rvec,ld0x*ld0x*3,ld0x,ld0x,ld0x*3)
      endif

C --- Evals and evecs of bulk P-S ---
      allocate(sll(4*ld0x**2))
      allocate(wk(ld0x*max(2*ld0x,6)))
      if (ipr >= PRTG) call yprm0('(1p,9e18.10)')
      if (mod(mode/10,10) == 1) then
        call pgblkk(ld0x,rvecl,sll,wk,rl)
        if (ipr >= PRTG) call yprm('r',2,rl,ld0x*2,ld0x*2,ld0x*2,1)
      else
        call pgblkk(ld0x,rvec,sll,wk,r)
        if (ipr >= PRTG) call yprm('r',2,r,ld0x*2,ld0x*2,ld0x*2,1)
      endif
      if (ipr >= PRTG) call yprm0('(9f15.10)')

C --- Bulk GF ---
      if (modeGF >= 0) then
        call pgcops(ld0x,ldlx,ldrx,ld0x,ld0x,opt,s00,sll)
        if (mod(mode/10,10) == 1) then
          call pgbulk(ld0x,rvecl,wk,rl,sll,modeGF,gbii,ndg,
     .      gii(1+ld0x*offg),lgii)
          deallocate(rvecl,rl)
        else
          call pgbulk(ld0x,rvec,wk,r,sll,modeGF,gbii,ndg,
     .      gii(1+ld0x*offg),lgii)
        endif
      endif
      deallocate(sll,wk)

C --- Recopy planar strux into sbulk ---
      if (mod(mode/100,10) == 1) then
        call pgcops(ld0x,ldlx,ldrx,ld0x,ld0x,opt,s00,sbulk)
      endif
      deallocate(s00)

      call tcx('pgbevl')
      end
      subroutine pgbulk(ld,rvec,w,r,s,mode,gb,ndg,g,lg)
C- Bulk or surface Green's function from evecs of propagation matrix
C ----------------------------------------------------------------
Ci Inputs:
Ci   ld:   dimension of bulk PL
Ci   s:    planar structure constants (see plstrx).
Ci   r:    log of the wave vectors (PRB 39, 923 (1989)).
Ci   w:    real work space of dimension ld*ld*2
Ci   rvec: eigenvectors of propagation matrix  (see pgblkk).
Ci   mode: determines what quantities to be calculated.
Ci     ones digit:
Ci       0,1  Q^-1 in rvec(1,2), else assume Q^-1 already generated
Ci       0,2  P^-1 in rvec(2,2), else assume P^-1 already generated
Ci     tens digit:
Ci         0  overwrite s0L with s0L Q x^-1 Q^-1
Ci         1  overwrite s0L with Q x^-1 Q^-1
Ci         2,3 do not make  Q r Q^-1 or s0L Q r Q^-1.  If needed:
Ci         2  assumes s0L already contains s0L Q x^-1 Q^-1
Ci         3  assumes s0L already contains Q x^-1 Q^-1
Ci     100s digit:
Ci         0  overwrite s0R with s0R P r P^-1
Ci         1  overwrite s0R with P r P^-1
Ci         2,3 do not make  P r P^-1 or s0R P r P^-1.  If needed:
Ci         2  assumes s0R already contains s0R P r P^-1
Ci         3  assumes s0R already contains P r P^-1
Ci    1000s digit:
Ci         0 no GF calculated
Ci         1 left  semiinfinite GF
Ci         2 right semiinfinite GF
Ci         3 bulk GF
Ci     ... Add 4 if to generate g^-1 by suppressing final inversion
Ci         8   left inverse semiinfinite GF; also put bulk GF in gb
Ci         9  right inverse semiinfinite GF; also put bulk GF in gb
Ci         NB: these allowed only if PrP^-1 and QrQ^-1 generated
Co Outputs:
Co   s0L,s0R see tens and 100s digit of mode
Co   g:    bulk or semi-infinite Green's function, depending on mode
Co   gb:   bulk Green's function, depending on 1000s digit of mode
Co   lg:   untouched if gii not calculated
Co         1 gii holds left  semi-infinite GF from s00,s0L,s0R
Co         2 gii holds right semi-infinite GF from s00,s0L,s0R
Co         3 gii holds bulk GF
Co     ... Add 4 if gii holds inverse of g
Co   rvec: (1,1) and (1,2) left untouched (input as Q,P)
Co         (2,1): Q^-1  (2,2) P^-1, depending on ones digit of mode
Cr Remarks
Cr   s00,  s0L and/or s0R OVERWRITTEN; see mode.
Cr   Bulk GF is (P-s00 - s0L Q x^-1 Q^-1 - s0R PrP^-1)^-1
C ----------------------------------------------------------------
      implicit none
      integer ld,mode,lg,ndg
      double precision w(ld,2),g(ld,ndg,2),gb(ld,ld,2),
     .  s(ld,ld,3,2),rvec(ld,2,ld,2,2),r(ld,2,2)
      integer i,j,iprint,mode1,mode2,mode3,mode4,ld2
      integer PRTG
      parameter (PRTG=80)
      double precision xxr,xxi
C      double precision det(4)
      logical linv,lbgf

C      call yprm('r',2,r,ld*2*1,ld*2,ld*2,1)
C      call yprm0('(1p,9e18.10)')
C      call yprm('Q',2,rvec,ld*2*ld*2,ld*2,ld,ld)
C      call yprm('P',2,rvec(1,1,1,2,1),ld*2*ld*2,ld*2,ld,ld)


      call tcn('pgbulk')
      mode1 = mod(mode,10)
      mode2 = mod(mode/10,10)
      mode3 = mod(mode/100,10)
      mode4 = mod(mode/1000,10)
      ld2 = ld**2
C ... lbgf is true if to put bulk GF in gb
      if (mode4 < 8) then
        mode4 = mod(mode4,4)
        lbgf = .false.
      else
        mode4 = mode4 - 7
        lbgf = mod(mode2,2) == 1 .and. mod(mode3,2) == 1
      endif

C ... Setup if we are going to make some kind of GF
      if (mode4 > 0) then
        call rxx(mod(lg,10) == 0,'pgbulk: no memory allocated for G')
        lg = 0
        linv  = mod(mode/1000,10) > 4
      endif

      if (mode1 == 0 .or. mode1 == 1) then
C   --- x -> x^-1 ---
        do  10  i = 1, ld
          call cdiv(1d0,0d0,r(i,1,1),r(i,1,2),xxr,xxi)
          r(i,1,1) = xxr
          r(i,1,2) = xxi
   10   continue
C   --- Make Q^-1 in rvec(2,1) ---
        do  12  j = 1, ld
        do  12  i = 1, ld
        rvec(i,2,j,1,1) = rvec(i,1,j,1,1)
   12   rvec(i,2,j,1,2) = rvec(i,1,j,1,2)
        call yyqinv('n',rvec(1,2,1,1,1),rvec(1,2,1,1,2),ld*2,2,ld,w,
     .    ld,i)
        if (i /= 0) call rx(' PGBULK: Q is singular')
C        call yprm0('(1p,9e18.10)')
C        call yprm('Q^-1',2,rvec(1,2,1,1,1),ld*2*ld*2,ld*2,ld,ld)
C   ... Invert Q using semistandard linpack
C        call yygefa(rvec(1,2,1,1,1),rvec(1,2,1,1,2),ld*2,ld,w,i)
C        if (i /= 0) call rx(' PGBULK: Q is singular')
C        call yygedi(rvec(1,2,1,1,1),rvec(1,2,1,1,2),ld*2,ld,w,det,
C     .    w(1,2),w(1,3),1)
C        call yprm0('(1p,9e18.10)')
C        call yprm('Q^-1',2,rvec(1,2,1,1,1),ld*2*ld*2,ld*2,ld,ld)
      endif

C --- Make P^-1 in rvec(2,2) ---
      if (mode1 == 0 .or. mode1 == 2) then
        do  14  j = 1, ld
        do  14  i = 1, ld
        rvec(i,2,j,2,1) = rvec(i,1,j,2,1)
   14   rvec(i,2,j,2,2) = rvec(i,1,j,2,2)
        call yyqinv('n',rvec(1,2,1,2,1),rvec(1,2,1,2,2),ld*2,2,ld,
     .    w,ld,i)
        if (i /= 0) call rx(' PGBULK: P is singular')
C   ... Invert P using semistandard linpack
C        call yygefa(rvec(1,2,1,2,1),rvec(1,2,1,2,2),ld*2,ld,w,i)
C        if (i /= 0) call rx(' PGBULK P is singular')
C        call yygedi(rvec(1,2,1,2,1),rvec(1,2,1,2,2),ld*2,ld,w,det,
C     .    w(1,2),w(1,3),1)
C        call yprm0('(1p,9e18.10)')
C        call yprm('P^-1',2,rvec(1,2,1,2,1),ld*2*ld*2,ld*2,ld,ld)
      endif

C --- Copy P-s00 to g ---
      if (mode4 > 0) then
        do  22  j = 1, ld
        do  22  i = 1, ld
        g(i,j,1) = -s(i,j,1,1)
   22   g(i,j,2) = -s(i,j,1,2)
      endif

C --- s00 <- s0L; Overwrite s0L with  s0L QrQ^-1 or QrQ^-1 ---
      if (mode2 <= 1) then
C   ... Make Q x^-1 in w
        do  30  j = 1, ld
        do  30  i = 1, ld
        w(i,j)    = rvec(i,1,j,1,1)*r(j,1,1) - rvec(i,1,j,1,2)*r(j,1,2)
   30   w(i,j+ld) = rvec(i,1,j,1,2)*r(j,1,1) + rvec(i,1,j,1,1)*r(j,1,2)
C   ... Make Q x^-1 Q^-1 in s00
        call yygemm('N','N',ld,ld,ld,1d0,w,w(1,1+ld),ld,
     .    rvec(1,2,1,1,1),rvec(1,2,1,1,2),ld*2,0d0,s,s(1,1,1,2),ld)
C       call yprm('QrQ^-1',2,s,ld*ld*3,ld,ld,ld)
C   ... Make s0L Q x^-1 Q^-1 in w:
        if (mode2 == 0 .or. mod(mode4,2) == 1 .or. lbgf) then
          call yygemm('N','N',ld,ld,ld,1d0,s(1,1,2,1),s(1,1,2,2),ld,
     .      s,s(1,1,1,2),ld,0d0,w,w(1,1+ld),ld)
C     ... Add -s0L Q x^-1 Q^-1 into g
          if (mod(mode4,2) == 1) then
            call daxpy(ld2,-1d0,w,1,g,1)
            call daxpy(ld2,-1d0,w(1,1+ld),1,g(1,1,2),1)
            lg = lg+1
          endif
C     ... if gb also generated, copy -s0L Q x^-1 Q^-1 into gb
          if (lbgf .and. mode4 == 2) then
            call dpcopy(w,gb,1,ld2,-1d0)
            call dpcopy(w(1,1+ld),gb(1,1,2),1,ld2,-1d0)
          endif
C     ... Copy s0L Q x^-1 Q^-1 into s00
          if (mode2 == 0) then
            call dcopy(ld*ld,w,1,s,1)
            call dcopy(ld*ld,w(1,1+ld),1,s(1,1,1,2),1)
          endif
        endif
C  ...  Copy s0L to s00;  Q x^-1 Q^-1 or s0L Q x^-1 Q^-1 to s0L
        do  32  j = 1, ld
        do  32  i = 1, ld
          xxr = s(i,j,1,1)
          xxi = s(i,j,1,2)
          s(i,j,1,1) = s(i,j,2,1)
          s(i,j,1,2) = s(i,j,2,2)
          s(i,j,2,1) = xxr
          s(i,j,2,2) = xxi
   32   continue
      endif

C --- s00 <- s0R; Overwrite s0R with  s0R PrP^-1  or  PrP^-1 ---
      if (mode3 <= 1) then
C   ... Make P r in w
        do  40  j = 1, ld
        do  40  i = 1, ld
        w(i,j)    = rvec(i,1,j,2,1)*r(j,2,1) - rvec(i,1,j,2,2)*r(j,2,2)
   40   w(i,j+ld) = rvec(i,1,j,2,2)*r(j,2,1) + rvec(i,1,j,2,1)*r(j,2,2)
C   ... Make P r P^-1 in s00
        call yygemm('N','N',ld,ld,ld,1d0,w,w(1,1+ld),ld,
     .    rvec(1,2,1,2,1),rvec(1,2,1,2,2),ld*2,0d0,s,s(1,1,1,2),ld)
C       call yprm('PrP-1',2,s,ld*ld*3,ld,ld,ld)
C   ... Make s0R PrP^-1 in w:
        if (mode3 == 0 .or. mode4 >= 2 .or. lbgf) then
          call yygemm('N','N',ld,ld,ld,1d0,s(1,1,3,1),s(1,1,3,2),ld,
     .      s,s(1,1,1,2),ld,0d0,w,w(1,1+ld),ld)
C     ... Add -s0R PrP^-1 into g
          if (mode4 >= 2) then
            call daxpy(ld2,-1d0,w,1,g,1)
            call daxpy(ld2,-1d0,w(1,1+ld),1,g(1,1,2),1)
            lg = lg+2
          endif
C     ... if gb also generated, copy s0R PrP^-1 into gb
          if (lbgf .and. mode4 == 1) then
            call dpcopy(w,gb,1,ld2,-1d0)
            call dpcopy(w(1,1+ld),gb(1,1,2),1,ld2,-1d0)
          endif
C     ... Copy s0R PrP^-1 into s00
          if (mode3 == 0) then
            call dcopy(ld*ld,w,1,s,1)
            call dcopy(ld*ld,w(1,1+ld),1,s(1,1,1,2),1)
          endif
        endif
C  ...  Copy s0R to s00;  P r P^-1 or s0R P r P^-1 to s0R
        do  42  j = 1, ld
        do  42  i = 1, ld
          xxr = s(i,j,1,1)
          xxi = s(i,j,1,2)
          s(i,j,1,1) = s(i,j,3,1)
          s(i,j,1,2) = s(i,j,3,2)
          s(i,j,3,1) = xxr
          s(i,j,3,2) = xxi
   42   continue
      endif

C --- Make bulk g^-1 by adding semi-infinite g^-1 in to g^-1 ---
      if (lbgf) then
        call daxpy(ld*ld,1d0,g,1,gb,1)
        call daxpy(ld*ld,1d0,g(1,1,2),1,gb(1,1,2),1)
        if (iprint() >= PRTG) call yprm('gb^-1',2,gb,ld*ld,ld,ld,ld)
      endif

C --- Exit if no GF calculated ---
      if (mode4 == 0) call tcx('pgbulk')
      if (mode4 == 0) return

C --- Make g from g^-1 ---
      if (.not. linv) then
        call yyqinv('n',g,g(1,1,2),ld,2,ld,w,ld,i)
        call rxx(i /= 0,' PGBULK: GF singular')
C       call yygefa(g,g(1,1,2),ld,ld,w,i)
C       call rxx(i /= 0,' PGBULK: GF singular')
C       call yygedi(g,g(1,1,2),ld,ld,w,det,w(1,2),w(1,3),1)
      else
        lg = lg+4
      endif

      call lgupac(lg,'p',1,mod(lg,4),lg/4,3,0,0)
      call lgupac(lg,' pgbulk generated',6,0,0,0,0,0)
      if (iprint() >= PRTG) call yprm('g',2,g,ld*ndg,ld,ld,ld)

      call tcx('pgbulk')
      end
      subroutine pgblkk(ld,s,rmat,wk,r)
C- Secular matrix for bulk wave vectors along the principal layer axis
C ----------------------------------------------------------------
Ci Inputs:
Ci   ld:   dimension of bulk strux
Ci   s:    structure constants
Ci   rmat: complex work space of dimension 4*ld**2
Ci   wk:   real work space of dimension 6*ld
Co Outputs:
Co   r:    r, which log are the wave vectors (PRB 39, 923 (1989).
Co   s:    eigenvectors of rmat, needed for PGF.
Cr Remarks
Cr   s       must be dimensioned at least 4*ld**2
Cr   In Chen's notation:
Cr     F_A+ g_-0    +  A  g_00      + F_A  g_+0 =    -1 if m eq n
Cr     F_A+ g_m-1,n +  A  g_mn      + F_A  g_m+1,n = 0  if m ne n
Cr   In P-S notation
Cr     -S_0- g_-0 _ (S-P)_00 g_00 - S_0+ g_+0 =  1
Cr   So
Cr      S_0- = S_10 <-> F_A+, S_0+ = S_01 <-> F_A, (S-P) <-> A
Cr   Eigenvalue equations for r:
Cr     (     0              1       )  ( Q  P )      ( Q  P )
Cr     (                            )  (      )  = r (      )
Cr     ( -s01^-1 s10  -s01^-1 (S-P) )  ( Q' P')      ( Q' P')
Cr
Cr   Simple one-orbital example: file dat=% rows 2 cols 2 complex
Cr                                        0   1
Cr                                        -1 -1/ta*(-z)
Cr                                        0 0
Cr                                        0 -1/ta*(-zi)
Cr
Cr   Make g:  (assumes 1st eval r1 is > 1 and 2nd r2 is < 1)
Cr     set gb =  '-1:1 -sz,zi r1 -i r2 -+ -sta -- -i'
Cr     set agb2 = '-1:1 -sz,zi -p -x -s-1 -1:1 -s4 -sta -sta -+ -i -s-1'
Cr     set gi =  '-1:1 -sz,zi r2 -sta -- -i'
Cr     mc -vta=3 -vz='2*ta-1' -vzi=0.0001 dat -evl -sub 1,1,1,1
Cr     -a r2 r2 -i -a r1 $gb -p -x -w .  $agb2
Cr   Analytic gb =  -i/sqrt(4*ta^2-z^2)
Cr               -> -i/sqrt(4*ta)/sqrt(eps) for eps=2*ta-z -> 0
C ----------------------------------------------------------------
      implicit none
      integer ld
      double precision wk(ld,2,3),s(ld,ld,3,2),
     .  rmat(ld,2,ld,2,2),r(ld,2,2)
      integer i,j,offi,lokr,lokc,lfudge,iprint,lgunit,ld2
      double precision fudge,abnrm
C     double precision det(4),alpha(2),beta(2)
C     integer fopna

      call tcn('pgblkk')
      ld2 = ld*2

C --- Check for singular s01: look for zero rows or columns ---
      fudge = 1d-6
      lfudge = 0
      do  i = 1, ld
        lokr = 0
        lokc = 0
        do  j = 1, ld
          if (s(i,j,3,1) /= 0 .or. s(i,j,3,2) /= 0) lokr=1
          if (s(j,i,3,1) /= 0 .or. s(j,i,3,2) /= 0) lokc=1
        enddo
        if (lokr == 0 .or. lokc == 0) lfudge = lfudge+1
      enddo

C ... If lfudge>0, add fudge to diagonal
      if (lfudge /= 0) then
        if (iprint() >= 50) call awrit2('%N pgblkk had %i zero row or'
     .    //' columns: add %g',' ',80,lgunit(1),lfudge,fudge)
        do  j = 1, ld
          s(j,j,2,1) = s(j,j,2,1) - fudge
          s(j,j,3,1) = s(j,j,3,1) - fudge
        enddo
      endif

C --- rmat(1,1) <- -S_01^-1 ---
      do  i = 1, ld
      do  j = 1, ld
        rmat(i,1,j,1,1) = -s(i,j,3,1)
        rmat(i,1,j,1,2) = -s(i,j,3,2)
      enddo
      enddo

C     Lapack a little slower, but more stable
      call ztoyy(rmat,ld2,ld2,ld,ld,0,1)
C     call zprm('S01',2,rmat,ld2,ld,ld)
      call zgetrf(ld,ld,rmat,ld2,wk,i)
      if (i /= 0) call rx(' PGBLKK: matrix singular')
      call zgetri(ld,rmat,ld2,wk,rmat(1,1,1,1,2),ld*ld,i)
      call ztoyy(rmat,ld2,ld2,ld,ld,1,0)

C     Old linpack routines
C     call yygefa(rmat,rmat(1,1,1,1,2),ld2,ld,wk,i)
C     if (i /= 0) call rx(' PGBLKK: matrix singular')
C     call yygedi(rmat,rmat(1,1,1,1,2),ld2,ld,wk,det,wk(1,2,1),
C    .  wk(1,3,1),1)

C     call yyqinv('n',rmat,rmat(1,1,1,1,2),ld2,2,ld,
C    .  rmat(1,2,1,1,1),ld2,i)
C     if (i /= 0) call rx(' PGBLKK: S singular')

C --- -S_01^-1 S_10 in rmat(2,1) ---
C     NB: used fudged s in each case, for a consistent hamiltonian
      call yygemm('N','N',ld,ld,ld,1d0,rmat,rmat(1,1,1,1,2),ld2,
     .  s(1,1,2,1),s(1,1,2,2),ld,0d0,rmat(1,2,1,1,1),rmat(1,2,1,1,2),ld2)

C --- Undoes to 2nd order error in (s01 + I*fudge)^-1 and restores s ---
      if (lfudge /= 0) then
C       Even though eigenvalues of s01^-1 (P-S) are better,
C       eigenvalues of assembled rmat are worse ... don't understand why
C        alpha(1) = fudge * 0
C        alpha(2) = 0
C        beta(1) = 1
C        beta(2) = 0
C        call yygemm('N4','N',ld,ld,ld,1d0,rmat,rmat(1,1,1,1,2),ld2,
C     .    rmat,rmat(1,1,1,1,2),ld2,0d0,
C     .    rmat(1,1,1,2,1),rmat(1,1,1,2,2),ld2)
C        call yygemm('N','N',ld,ld,ld,1d0,rmat,rmat(1,1,1,1,2),ld2,
C     .    rmat(1,1,1,2,1),rmat(1,1,1,2,2),ld2,0d0,
C     .    rmat(1,2,1,2,1),rmat(1,2,1,2,2),ld2)
C        call ymsadd(ld,ld,ld2,ld2,0,0,0,0,alpha,beta,rmat(1,1,1,2,1),
C     .    ld2**2,rmat,ld2**2)
C        alpha(1) = alpha(1)**2
C        call ymsadd(ld,ld,ld2,ld2,0,0,0,0,alpha,beta,rmat(1,2,1,2,1),
C     .    ld2**2,rmat,ld2**2)

        do  j = 1, ld
          s(j,j,2,1) = s(j,j,2,1) + fudge
          s(j,j,3,1) = s(j,j,3,1) + fudge
        enddo
      endif

C --- -S_01^-1 (S-P) in (2,2) block ---
      call yygemm('N','N',ld,ld,ld,1d0,rmat,rmat(1,1,1,1,2),ld2,
     .  s,s(1,1,1,2),ld,0d0,rmat(1,2,1,2,1),rmat(1,2,1,2,2),ld2)

C      i = fopna('out',-1,0)
C      call ywrm(0,'s10',2,i,'(9f15.6)',s(1,1,1,1),ld**2*3,ld,ld,ld)
C      call ywrm(0,' ',2,i,'(9f15.6)',rmat(1,1,1,1,1),ld2**2,ld2,ld,ld)
C      call ywrm(0,' ',2,i,'(9f15.6)',rmat(1,2,1,1,1),ld2**2,ld2,ld,ld)
C      call ywrm(0,' ',2,i,'(9f15.6)',rmat(1,2,1,2,1),ld2**2,ld2,ld,ld)
C      call fclose(i)
C      pause

C --- Zero in (1,1) block, one in (1,2) block ---
      do  i = 1, ld
      do  j = 1, ld
        rmat(i,1,j,1,1) = 0
        rmat(i,1,j,1,2) = 0
        rmat(i,1,j,2,1) = 0
        rmat(i,1,j,2,2) = 0
      enddo
      rmat(i,1,i,2,1) = 1
      enddo

C --- Eigenvalues and eigenvectors of rmat ---
      offi = 4*ld**2
C     i = fopna('out',-1,0)
C     call ywrm(0,' ',2,i,'(%9;8,8d)',rmat(1,1,1,1,1),offi,ld2,ld2,ld2)
C     call fclose(i)

      call ztoyy(rmat,ld2,ld2,ld2,ld2,0,1)
C     call zprm('rmat',2,rmat,ld2,ld2,ld2)
      call zgeevs('P','N','V','N',ld2,rmat,ld2,r,r,1,s,ld2,abnrm,i)
      call ztoyy(r,ld2,1,ld2,1,1,0)
      call ztoyy(s,ld2,ld2,ld2,ld2,1,0)

C      call tcn('cg')
C      call cg(ld2,ld2,rmat,rmat(1,1,1,1,2),r,r(1,1,2),1,
C     .  s,s(offi+1,1,1,1),wk,wk(1,1,2),wk(1,1,3),i)
C      call tcx('cg')
C      if (i /= 0) call rx(' PGBLKK: failed to find evecs to r')

C --- Sort eigenvalues/vectors by increasing |e| ---
      do  i = 1, 2*ld
        wk(2*i-1,1,1) = r(i,1,1)
        wk(2*i,1,1)   = r(i,1,2)
      enddo
      call dvshel(2,2*ld,wk,wk(4*ld+1,1,1),11)
      call xxblkk(ld2,wk(4*ld+1,1,1),rmat,r,s,s(offi+1,1,1,1))

      call tcx('pgblkk')
      end
      subroutine xxblkk(ld2,ipm,wk,r,sr,si)
C- Sorts evals, evecs by increasing r
      implicit none
      integer ld2,ipm(4)
      double precision r(ld2,2),sr(ld2,ld2),si(ld2,ld2),wk(ld2,ld2,2),ax
      integer i,j,k,ld,m
      double complex xx,yy

      ld = ld2/2

C --- Sort r smallest to midpoint ---
      j = ld
      do  i = 1, ld
        k = ipm(i)+1
        j = j+1
        wk(j,1,1) = r(k,1)
        wk(j,2,1) = r(k,2)
      enddo

C --- Sort r midpoint to largest ---
      j = ld+1
      do  i = ld+1, ld2
        k = ipm(i)+1
        j = j-1
        wk(j,1,1) = r(k,1)
        wk(j,2,1) = r(k,2)
      enddo
      call dcopy(ld2*2,wk,1,r,1)
C     call yprm('r',2,wk,ld2,ld2,ld2,1)

C --- Adjust smallest evals to be reciprocals of largest ---
      do  40  i = 1, ld
        xx = dcmplx(r(i,1),r(i,2))
        ax = abs(xx)
        if (ax > 100d0) then
C     ... skip degeneracies for now
          if (i > 1) then
            if (abs(abs(dcmplx(r(i-1,1),r(i-1,2)))/ax - 1d0) < 1d-4)
     .        goto 40
          endif
          if (i < ld) then
            if (abs(abs(dcmplx(r(i+1,1),r(i+1,2)))/ax - 1d0) < 1d-4)
     .        goto 40
          endif
          yy = 1/xx
          ax = abs(yy/dcmplx(r(i+ld,1),r(i+ld,2))-1)
          if (ax > 1d-5) then
C            call awrit3(' xxblkk: eval %i err=%;3,3g  evl= %2:1;3,3g',
C     .        ' ',80,i1mach(2),i,ax,yy)
C            r(i+ld,1) = dble(yy)
C            r(i+ld,2) = dimag(yy)
          endif
        endif
   40 continue

C --- Sort P,Q smallest to midpoint ---
      j = ld
      do  i = 1, ld
        k = ipm(i)+1
        j = j+1
        do  m = 1, ld
          wk(m,j,1) = sr(m,k)
          wk(m,j,2) = si(m,k)
        enddo
      enddo

C --- Sort P,Q midpoint to largest ---
      j = ld+1
      do  i = ld+1, ld2
        k = ipm(i)+1
        j = j-1
        do  m = 1, ld
          wk(m,j,1) = sr(m,k)
          wk(m,j,2) = si(m,k)
        enddo
      enddo

C --- Copy back ---
      do  k = 1, ld2
      do  m = 1, ld
        sr(m,k) = wk(m,k,1)
        si(m,k) = wk(m,k,2)
      enddo
      enddo

      end
