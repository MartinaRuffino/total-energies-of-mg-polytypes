      subroutine pgsif(mode,s_ctrl,s_lat,s_ham,s_pot,s_str,s_spec,
     .  s_site,ipl,pgplp,qp,isp,ldg,ndg,lgii,gii)
C- Make semi-infinite GF, or related matrix, for one energy, qp, spin
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nbas nbasp npl nl lncol lpgf
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  plat
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  ldham neula lgen3 lncol lham offH
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iprmb eula
Cio    Passed to:  plhamnc plham plhamso pgbevl
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  pfnc
Cio    Elts passed:palp pfnc
Cio    Passed to:  plhamnc plham pgbevl
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:iax npr s
Cio    Passed to:  plhamnc plham pgbevl
Cio  s_spec (not used)
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  norb
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pfr
Cio    Passed to:  plhamnc plham plhamso pgbevl
Ci Inputs
Ci   mode  :1s digit:
Ci          0 for left sif; 1 for right sif
Ci      Add 2 for bulk gf
Ci      Add 4 to return S_L+1,L g_LL S_L,L+1 or S_R-1,R g_RR S_R,R-1
Ci         10s digit:
Ci          0 by decimation
Ci          1 by method of PRB 39, 923 (1989).
Ci        100s digit distinguishes form of complex storage for gii.
Ci          0:real, imaginary separated: gii(ld0x,ld0x,1..2)
Ci          1:complex*16 format: gii(1..2,ld0x,ld0x)
Ci          2:real, imaginary separated by columns: gii(ld0x,2,ld0x)
Ci   ipl   :index to current principal layer
Ci          -1 => strux for left bulk PL; npl => strux for right bulk
Ci   pgplp :index and dimensioning information for each PL (pgfset.f)
Ci   qp    :k-point
Ci   isp   :current spin channel (1 or 2)
Ci   ldg   :not used now
Ci   ndg   :not used now
Co Outputs
Co   lgii  :attributes of the gii generated; see lgupac for conventions
Co   gii   :diagonal GF for this principal layer
Co          gii is dimensioned
Cr Remarks
Cb Bugs
Cb  more consistent treatment of complex arithemetic needed
Cb  move ymscop from here to caller
Cu Updates
Cu      May 16 noncollinear bulk GF can be calculated by normal modes
Cu   17 May 16 bulk GF can be noncollinear; it can be calculated by decimation
Cu   30 Apr 16 Surface GF can be noncollinear
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
Cu   10 Jan 02 pgsif now assumes pot functions already made
Cu   21 Dec 01 Altered argument list
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,ipl,isp,lgii,pgplp(6,-1:*)
      integer ldg,ndg  ! not used
      double precision qp(3),gii(*)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_lat)::   s_lat
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_str)::   s_str
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
C ... Dynamically allocated local arrays
      real(8), allocatable :: wk2(:),gsLL(:,:,:),gsRR(:,:,:)
      complex(8), allocatable :: sll(:),snew(:),sllbk(:)
      complex(8), allocatable :: wk(:)
      complex(8), pointer :: pfa(:,:)
C ... Local parameters
C     logical :: debug = .true.
      logical lso
      integer PRTG,i,ipr,ld0,ld2,lds,modeGF,maxit,mode0,mode1,kcplx,nbas,nbasp,nl,npl,lncol
      integer nspc,ld0x,ldsx,ld2x,lgLL,lgRR,lhdimp,neul,lidimp
      double precision xx
      parameter (maxit=50,PRTG=80)
      double precision plat(3,3)
      character*40 strn
      procedure(logical) bittst

C --- Setup potential, lattice parameters ---
      call tcn('pgsif')
      call getpr(ipr)
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      kcplx = mod(mode/100,10)
      plat = s_lat%plat
      nbas = s_ctrl%nbas
      nbasp = s_ctrl%nbasp
      npl = s_ctrl%npl
      nl = s_ctrl%nl
      lncol = s_ctrl%lncol
      lso = bittst(lncol,4)
      nspc = s_ham%ldham(4)
      neul = s_ham%neula
      lidimp = s_ham%ldham(2)
      if (s_ctrl%lpgf(2) < 10) nspc = 1  ! Generate collinear surface GF

      if (ipl /= -1 .and. ipl /= npl)
     .  call rx('pgsif not implemented for 0<ipl<npl')

C     Assemble s_pot%pfnc.  Note parallel lines in pgflu, pgemb
      if (lso) then  ! plhamso does rotation internally; not needed here
      elseif (nspc == 2) then
        lhdimp = s_ham%ldham(3)
        call ptr_pot(s_pot,8+1,'pfnc',lhdimp,nspc**2,[0d0])
        pfa => s_pot%palp
        call dpscop(pfa,s_pot%pfnc,lhdimp*2,1,1,1d0)
        call dpscop(pfa,s_pot%pfnc,lhdimp*2,1+2*lhdimp,1+2*3*lhdimp,1d0)
        call rotspn(230110,1,nbasp,nbasp,nl,s_ham%iprmb,s_ham%eula,neul,[0d0],[0d0],
     .    [0d0],lhdimp,lidimp,lidimp,lhdimp,1,[0d0],s_pot%pfnc)
      endif
C     lhdimp = s_ham%ldham(3); call zprm('pfnc',2,s_pot%pfnc,lhdimp,lhdimp,4)

      ld0 = pgplp(4,ipl)
      lds = 3*ld0
      ld2 = ld0**2
      ld0x = ld0*nspc
      ldsx = lds*nspc
      ld2x = ld2*nspc*nspc

      strn = 'L surface g'
      if (mod(mode0,2) == 1) strn(1:1) = 'R'
      if (bittst(mode0,2)) strn = 'bulk gf'

C --- Decimation ---
      If (mode1 == 0) then
C       call rxx(bittst(mode0,2),'pgsif: no bulk gf with decimation')

C   ... Allocate memory for and make 2D Bloch hamiltonian, global spin quantization axis. kcplx=0 now
        allocate(sll(3*ld2x),snew(3*ld2x)) ! snew is a work array
C       call pshpr(81)
        call plhamnc(s_ham,s_pot,s_str,s_site,0,plat,isp,nspc,0,ipl,npl,ld0,lds,pgplp,qp,snew,sll)
C       call poppr

        allocate(wk(ld2x),wk2(ldsx))

C   ... Bulk g by decimation
        if (bittst(mode0,2)) then
C         call rx('pgsif: no bulk gf with decimation yet')
          allocate(sllbk(3*ld2x),gsLL(ld0x,ld0x,2),gsRR(ld0x,ld0x,2))
          call dcopy(2*3*ld2x,sll,1,sllbk,1)
          lgii = 1  ! Left surface
C#ifdef QUAD
          call qpgsifd(ld0x,sllbk,snew,wk,wk2,lgii,maxit,ld0x,gsLL)
C#elseC
C          call pgsifd(ld0x,sllbk,snew,wk,wk2,lgii,maxit,ld0x,gsLL)
C#endif
          call lgupac(lgLL,'p',1,1,0,3,0,0)
C          print *, isp
C          call yprm('gsLL',2,gsLL,ld0x**2,ld0x,ld0x,ld0x)
          call dcopy(2*3*ld2x,sll,1,sllbk,1)
          lgii = 2  ! Right surface
C#ifdef QUAD
          call qpgsifd(ld0x,sllbk,snew,wk,wk2,lgii,maxit,ld0x,gsRR)
C#elseC
C          call pgsifd(ld0x,sllbk,snew,wk,wk2,lgii,maxit,ld0x,gsRR)
C#endif
          call lgupac(lgRR,'p',1,2,0,3,0,0)
C         call yprm('gsRR',2,gsRR,ld0x**2,ld0x,ld0x,ld0x)
          call pgembs(ldsx,ld0x,ld0x,ld0x,ld0x,ld0x,ld0x,0,0,0,
     .      sll,30,lgLL,gsLL,lgRR,gsRR,sllbk,lgii,gii)
C          print *, isp
C         call yprm('gii',2,gii,ld0x**2,ld0x,ld0x,ld0x)
C   ... Semi-infinite g by decimation
        else
          lgii = mod(mode0,2)+1 ! 1 for left, 2 for right
C#ifdef QUAD
          call qpgsifd(ld0x,sll,snew,wk,wk2,lgii,maxit,ld0x,gii)
C#elseC
C          call pgsifd(ld0x,sll,snew,wk,wk2,lgii,maxit,ld0x,gii)
C#endif
          call lgupac(lgii,'p',1,lgii,0,3,0,0)
C          if (debug) then
C            call yprm('gs',2,gii,ld0x**2,ld0x,ld0x,ld0x)
C          endif
        endif

        deallocate(snew,wk,wk2,sll)

C --- Diagonalization of propagation matrix ---
      elseif (mode1 == 1) then
        allocate(sll(4*ld2))
        if (mod(mode0,4) >= 2) then
          modeGF = 3000   ! Bulk GF
C         call rxx(nspc == 2,'please use decimation for noncollinear xtal GF')
        elseif (mod(mode0,2) == 0) then
          modeGF = 1201   ! Left semi-infinite GF
        else
          modeGF = 2022   ! Right semi-infinite GF
        endif
C   ... Propagation coefficients r for bulk GF
        i = 11
        if (ipl == -1) i = 10

        if (lso) call info0(0,1,0,'pgsif (warning) pgbevl not designed for SO coupling')

C   ... until pgbevl is revised
        call pgbevl(s_ham,s_pot,s_str,s_site,plat,isp,nspc,ipl,npl,pgplp,qp,i,
     .    modeGF,xx,xx,xx,xx,ld0x,0,gii,lgii)
        deallocate(sll)
      endif

      if (ipr >= PRTG/1) then
        call awrit2('%a (%?#n==0#decimation#diagonalization#),'//
     .    ' layer %i',strn,len(strn),0,mode1,ipl)
        call yprm(strn,2,gii,ld0x**2,ld0x,ld0x,ld0x)
      endif

C --- Overwrite g with S g S+ ---
      if (bittst(mode0,4)) then
C   ... Allocate memory for and make 2D Bloch hamiltonian
        allocate(sll(3*ld2))
        call plham(0,s_ham,s_pot,s_str,s_site,plat,isp,0,ipl,pgplp,qp,lds,ld0,sll)
C   ... Overwrite gii with S g S+
        call pgsif2(mode,sll,ld0x,gii)
        if (ipr >= PRTG)
     .    call yprm('S g S+',kcplx+2,gii,ld0x*ld0x,ld0x,ld0x,ld0x)
        deallocate(sll)
      elseif (kcplx /= 0) then
        call ztoyy(gii,ld0x,ld0x,ld0x,ld0x,0,1)
C        if (ipr >= PRTG/100)
C     .    call yprm('g',kcplx+2,gii,ld0x*ld0x,ld0x,ld0x,ld0x)
      endif

      call tcx('pgsif')
      end
      subroutine pgsif2(mode,s,ld,gii)
C- Overwrite g00 with S0,-1 g-1-1 S-1,0 or gnn with Sn-1,n gnn Sn-1,n
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit:
Ci          0 for left sif; 1 for right sif
Ci         10s digit: not used
Ci        100s digit distinguishes form of complex storage for gii.
Ci          0:real, imaginary separated: gii(ldg,ldg,1..2)
Ci          1:complex*16 format: gii(1..2,ldg,ldg)
Ci          2:real, imaginary separated by columns: gii(ldg,2,ldg)
Ci   s     :layer strux or hamiltonian; see Bugs
Ci   ld    :dimensions of s,gii
Co Inputs/Outputs
Ci   gii   :(input) diagonal (surface) GF for this principal layer; see Bugs
Co         :(output) S g S+
Cr Remarks
Cr   Routine generates  S g S+, with :
Cr      S = S_i,i-1 when g is left semi-infinite
Cr      S = S_i+1,i when g is right semi-infinite
Cr
Cr   Adding this term to the hamiltonian z-H changes g = (z-H)^-1
Cr   from g for a free boundary to that coupling to a semi-infinite one
Cb Bugs
Cb   s is destroyed in complex*16 format
Cb
Cb   For now, input s,gii are in kcplx=0 mode always.
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,ld
      double precision s(ld,ld,3,2),gii(ld,ld,2)
C ... Local parameters
      integer L,R,kcplx,ld1,ld2,ldr,oi

      call tcn('pgsif2')
      if (mod(mod(mode,10),2) == 0) then
C       s(L) = s0,-1; s(R) = s-1,0
        L = 2
        R = 3
      else
C       s(R) = s_n,n-1; s(L) = s_n-1,n
        L = 3
        R = 2
      endif
      kcplx = mod(mode/100,10)
C     make gii internally in kcplx=0 or kcplx=2 mode
      if (kcplx == 1) kcplx = 2  ! Imag following real in columns
      call cplxdm(kcplx,ld,ld,ld1,ld2,ldr,oi)

C ... s00 <- s0,-1 g-1-1 or s_n-1,n g_nn s_n,n-1 (kcplx=0)
      call yygemm('N','N',ld,ld,ld,1d0,s(1,1,L,1),s(1,1,L,2),ld,
     .  gii,gii(1,1,2),ld,0d0,s,s(1,1,1,2),ld)
C     call yprm('S g',2,s,ld*ld*3,ld,ld,ld*3)

C ... gii <- (s0,-1 g-1-1) s-1,0  (kcplx=0 or 2)
      call yygemm('N','N',ld,ld,ld,1d0,s(1,1,1,1),s(1,1,1,2),ld,
     .  s(1,1,R,1),s(1,1,R,2),ld,0d0,gii,gii(1+oi,1,1),ldr)
C     call yprm('S(0-) g(--) S(-0)',2+kcplx,gii,oi,ld,ld,ld)

C     Use s as a work array here: restore kcplx=2 -> kcplx=1
      kcplx = mod(mode/100,10)
      if (kcplx == 1) call ztoy(gii,ld,ld,ld,1)

      call tcx('pgsif2')
      end

C      subroutine pgsif3(s_ham,s_pot,s_str,plat,isp,ipl,npl,nofgL,
C     .  nofgR,pgplp,qp,g00x,lg00x,g00,lg00)
CC- From semi-infinite diagonal Green's function, make the off-diagonal
CCi   g00x:   used to make off-diagonal GF.  g00x must contain the
CCi           the off-diagonal semi-infinite g00.
CCu   10 Nov 11 Begin migration to f90 structures
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer isp,npl,lg00,lg00x,nofgR,nofgL,pgplp(6,-1:*)
C      double precision plat(3,3),qp(3),g00(*),g00x(*)
CC ... For structures
C      include 'structures.h'
C      type(str_ham)::   s_ham
C      type(str_pot)::   s_pot
C      type(str_str)::   s_str
CC ... Local parameters
C      integer ipl,ldim0,ldim,i,pgdim,offg,offg2,ndg,npl1,npl2,osll,owk,
C     .  nolo,nohi,ld2,offs,PRTG,ipr
C      parameter (PRTG=80)
CC heap:
C      integer w(1)
C      common /w/ w
C
C      nolo(i) = max(i-nofgL,-1-1)
C      nohi(i) = min(i+nofgR,npl+1)
C      npl1 = nolo(ipl)
C      npl2 = nohi(ipl)
C      if (npl1 == ipl-1 .or. npl2 == ipl+1) then
C
C        call tcn('pgsif3')
C        call getpr(ipr)
C        ndg = pgdim(0,ipl,npl,nolo(ipl),nohi(ipl),pgplp)
C        offg = pgdim(2,ipl,npl,npl1,ipl,pgplp)
C        ldim0 = pgplp(4,ipl)
C        ldim  = 3*ldim0
C        ld2   = ldim0**2
C        if (ipl == -1) then
C          offg2 = offg - ldim0
C          offs = ldim0
C        elseif (ipl == npl) then
C          offg2 = offg + ldim0
C          offs = 2*ldim0
CC       This probably isn't necessary.
C        else
C          call fexit2(-1,111,' Exit -1 PGSIF3: '//
C     .      'expected ipl to be -1 or %i but found %i',npl,ipl)
C        endif
C
CC       diagonal g00x was put for offset=0.  Shift to proper offset
CC       NB: this belongs in caller!
C        if (offg /= 0) then
C          call ymscop(1,ldim0,ldim0,ldim0,ldim0,0,0,0,offg,
C     .      g00x,ldim0*ndg,g00x,ldim0*ndg)
C        endif
C
C        call defdc(osll,3*ld2)
C        call plham(0,s_ham,s_pot,s_str,s_site,plat,isp,0,ipl,pgplp,qp,ldim,ldim0,w(osll))
C
C        call defdc(owk,ld2)
C        call pgdys2(.true.,2,'N',ldim,ldim0,ldim0,ldim0,ndg,ndg,ndg,
C     .    offg,offg,0,ldim0,offs,w(owk),w(osll),lg00x,g00x,lg00,g00,
C     .    g00(1+ldim0*offg2))
C
C        if (ipr >= PRTG) call yprm('off-diagonal SIF g',2,g00,ldim0*
C     .    ndg,ldim0,ldim0,ndg)
C
C        call rlse(osll)
C        call tcx('pgsif3')
C      endif
C
C      end

      subroutine pgsifd(ld,s,sp,wk,wk2,lg00,maxit,ndg,g00)
C- Makes semi-infinite Green's function by decimation
C ----------------------------------------------------------------
Ci Inputs
Ci   ld    :leading dimension of s,sp,g00
Ci   s     :-1 * layer hamiltonian
Ci   sp    :complex work array of same dimension as s
Ci   wk    :complex work array of dimension ld*ld
Ci   wk2   :double precision work array of dimension 3*ld
Ci   lg00   1 Generate left  semi-infinite GF from s00,s0L,s0R
Ci          2 Generate right semi-infinite GF from s00,s0L,s0R
Ci   maxit :maximum number of iterations to attempt
Ci   ndg   :second dimension of g00
Co Outputs
Ci   s     :is OVERWRITTEN on output
Co   g00   :semi-infinite GF
Cr Remarks
Cr   Decimation works iteratively by solving for all odd g_2n+1,0
Cr   in terms of the even g_2n,0 and g_2n+2,0.  Thus the equation
Cr     (P - S_00) g_00 - S_01 g_10 = 1
Cr   becomes
Cr    (P - S'_00) g_00  - S'_01 g_20  = 1
Cr    (P - S''_00) g_00 - S''_01 g_40 = 1
Cr   ...
Cr   and is repeated until the second term becomes negligible.
Cr   The scheme starts with the equations
Cr     (P - S_00) g_00 - S_01 g_10 = 1
Cr     -S_10 g_n-1,0 + (P - S_11) g_n0 - S_01 g_n+10 = 0 for all n
Cr   At the start, S_11 = S_00.
Cr   Every odd equation is used to knock out all the g_2n+1:
Cr     g+2n+1,0 = (P - S_11)^-1 (S_10 g_2n,0 + S_01 g_2n+2,0)
Cr   The equations involving only the g_2n look just as those with g_n
Cr   provided that the S_00, S_11, S_01 and S_10 are redefined:
Cr     S_00' =  S_00 - S_01 (P - S_11) S_10
Cr     S_11' =  S_11 + S_10 (P - S_11) S_01 + S_01 (P - S_11) S_10
Cr     S_10' =  S_10 (P - S_11) S_10
Cr     S_01' =  S_01 (P - S_11) S_01
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ld,maxit,ndg,lg00
      double precision wk2(ld,3),s(ld,ld,3,2),g00(ld,ndg,2),sp(ld,ld,3,2),wk(ld,ld,2)
C ... Local parameters
C     logical :: debug = .false.
      integer i,j,iter,info,lgunit,iprint,ld2,L,R
      double precision det(4),tol
      parameter (tol=1d-10)

C     print *, '!!' ; call setpr(70)

      call tcn('pgsifd')
      ld2 = ld*ld

C --- Put (P-S) in g_00 and S_11 ---
      call dpcopy(s,g00,1,ld2,-1d0)
      call dpcopy(s(1,1,1,2),g00(1,1,2),1,ld2,-1d0)
      call dpcopy(g00,s,1,ld2,1d0)
      call dpcopy(g00(1,1,2),s(1,1,1,2),1,ld2,1d0)
C      if (debug) then
C      call yprm('(P-s00)',2,s,ld*ld*3,ld,ld,ld)
C      endif

      L = 2
      R = 3

C --- For right GF, transpose S(2) and S(3) ---
      if (lg00 == 2) then
        do  j = 1, ld
          call dswap(ld,s(1,j,2,1),1,s(1,j,3,1),1)
          call dswap(ld,s(1,j,2,2),1,s(1,j,3,2),1)
        enddo
      endif

C --- Start of decimation loop ---
C     Start with P-S_00 in s(1) S_01 in s(2), S_10 in s(3)
      do  iter = 1, maxit
C   ... sp(1,*) <- (P-S11) (temporary storage)
        do  i = 1, ld
        do  j = 1, ld
          sp(i,j,1,1) = s(i,j,1,1)
          sp(i,j,1,2) = s(i,j,1,2)
        enddo
        enddo
C       call yprm('(P-s11)',2,sp,ld*ld*3,ld,ld,ld)

C   ... sp(1,*) <- (P-s11)^-1 (temporary storage)
        call yygefa(sp,sp(1,1,1,2),ld,ld,wk2,info)
        if (info /= 0) call fexit(-1,111,' Exit -1 PGSIFD: '//
     .    'matrix singular, iteration %i',iter)
        call yygedi(sp,sp(1,1,1,2),ld,ld,wk2,det,wk2(1,2),wk2(1,3),1)
C       call yprm('(P-s11)^-1',2,sp,ld*ld*3,ld,ld,ld)

C   ... sp(3,*) <- s01 (P-s11)^-1 (temporary storage)
        call yygemm('N','N',ld,ld,ld,1d0,s(1,1,L,1),s(1,1,L,2),ld,
     .    sp,sp(1,1,1,2),ld,0d0,sp(1,1,3,1),sp(1,1,3,2),ld)
C       call yprm('s01 (P-s11)^-1',2,sp(1,1,3,1),ld*ld*3,ld,ld,ld)

C   ... wk <- s01 (P-s11)^-1 s10 (renormalization of S_00)
        call yygemm('N','N',ld,ld,ld,1d0,sp(1,1,3,1),sp(1,1,3,2),ld,
     .    s(1,1,R,1),s(1,1,R,2),ld,0d0,wk,wk(1,1,2),ld)
C       call yprm('s01 (P-s11)^-1 s10',2,wk,ld*ld,ld,ld,ld)

C   ... sp(2,*) <- s01 (P-s11)^-1 s01 (updated s01)
        call yygemm('N','N',ld,ld,ld,1d0,sp(1,1,3,1),sp(1,1,3,2),ld,
     .    s(1,1,L,1),s(1,1,L,2),ld,0d0,sp(1,1,2,1),sp(1,1,2,2),ld)
C       call yprm('s01 (P-s11)^-1 s01',2,sp(1,1,2,1),ld*ld*3,ld,ld,ld)

C   ... sp(3,*) <- s10 (P-s11)^-1 (temporary storage)
        call yygemm('N','N',ld,ld,ld,1d0,s(1,1,R,1),s(1,1,R,2),ld,
     .    sp,sp(1,1,1,2),ld,0d0,sp(1,1,3,1),sp(1,1,3,2),ld)
C       call yprm('s10 (P-s11)^-1',2,sp(1,1,3,1),ld*ld*3,ld,ld,ld)

C   ... sp(1,*) <- s10 (P-s11)^-1 s01 (temporary storage)
        call yygemm('N','N',ld,ld,ld,1d0,sp(1,1,3,1),sp(1,1,3,2),ld,
     .    s(1,1,L,1),s(1,1,L,2),ld,0d0,sp,sp(1,1,1,2),ld)
C       call yprm('s01 (P-s11)^-1 s10',2,sp,ld*ld*3,ld,ld,ld)

C   ... RMS value of s01 (P-s11)^-1 s10
C   ... sp(1,*) <- s11 - s01 (P-s11)^-1 s10 - s10 (P-s11)^-1 s01
C   ... g00 <- s00 - s01 (P-s11)^-1 s10 (updated s11, and s00)
C   ... wk  <- sp(3,*) = s10 (P-s11)^-1 (temporary storage)
        call yydotc(ld2,wk,wk(1,1,2),1,wk,wk(1,1,2),1,det,det(2))
        do  j = 1, ld
        do  i = 1, ld
          sp(i,j,1,1) = s(i,j,1,1) - sp(i,j,1,1) - wk(i,j,1)
          sp(i,j,1,2) = s(i,j,1,2) - sp(i,j,1,2) - wk(i,j,2)
          g00(i,j,1) = g00(i,j,1) - wk(i,j,1)
          g00(i,j,2) = g00(i,j,2) - wk(i,j,2)
          wk(i,j,1) = sp(i,j,3,1)
          wk(i,j,2) = sp(i,j,3,2)
        enddo
        enddo
C       call yprm('updated s11',2,sp,ld*ld*3,ld,ld,ld)

C   ... sp(3,*) <- s10 (P-s11)^-1 s10 (updated s10)
        call yygemm('N','N',ld,ld,ld,1d0,wk,wk(1,1,2),ld,
     .    s(1,1,R,1),s(1,1,R,2),ld,0d0,sp(1,1,3,1),sp(1,1,3,2),ld)
C       call yprm('s10 (P-s11)^-1 s10',2,sp(1,1,3,1),ld*ld*3,ld,ld,ld)

C   ... RMS s10, s11
        det(2) = 0
        det(3) = 0
        do  j = 1, ld
        do  i = 1, ld
          det(2) = det(2) + sp(i,j,2,1)**2 + sp(i,j,2,2)**2
          det(3) = det(3) + sp(i,j,3,1)**2 + sp(i,j,3,2)**2
        enddo
        enddo
        do  i = 1, 3
          det(i) = dsqrt(det(i)/ld2)
        enddo
        if (iprint() >= 70 .or. iter == maxit)
     .    call awrit3(' pgfsid iter %i:  rms delta S00 = %1;4g, '//
     .    'rms S01,S10 =%2:1;4g',' ',80,lgunit(1),iter,det(1),det(2))
C   ... Copy updated s10,s01,s11 to new positions (g00 already updated)
        call dcopy(ld2*3*2,sp,1,s,1)
        if (det(1)+det(2)+det(3) < tol) exit
        if (iter == maxit) call info0(2,0,0,' PGSIFD (warning): GF not converged to tolerance ...')
      enddo

C --- Semi-infinite Green's function ---
      call yygefa(g00,g00(1,1,2),ld,ld,wk2,info)
      if (info /= 0) call fexit(-1,111,' Exit -1 PGSIFD: '//
     .    'Green''s function matrix singular',0)
      call yygedi(g00,g00(1,1,2),ld,ld,wk2,det,wk2(1,2),wk2(1,3),1)

C      if (debug) then
C        call yprm('surface GF',2,g00,ld*ndg,ld,ld,ld)
C      endif

      call tcx('pgsifd')
      end

