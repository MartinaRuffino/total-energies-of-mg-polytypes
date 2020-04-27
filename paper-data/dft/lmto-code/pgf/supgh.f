      subroutine supgh(s_ctrl,s_spec,s_ham,s_pot,s_str,pgplp)
C- Hamiltonian setup for layer GF
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nclasp nofgl nofgr nspec nbasp nl nspin lgen3 lncol
Ci                 npl lpgf
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:ips ipc
Cio    Passed to:  *
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa idxdn lmxb ncomp
Co     Stored:     lmxa idxdn
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  makidx nscpa
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  *
Co     Stored:     bandw
Co     Allocated:  *
Cio    Elts passed:offH lncol
Cio    Passed to:  pglusu
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pp
Cio    Passed to:  *
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:npr iax alph
Cio    Passed to:  *
Cio Inputs/Outputs
Cio   pgplp :index and dimensioning information for each PL (pgfset.f)
Cio         :On input pgplp contains
Cio         1: cumulative number of basis atoms in this PL=0 ... this PL
Cio         2: index to PL with inequivalent potential
Cio         :On output, this routines stores the following:
Cio         3: source (column) dimension of GF for this PL
Cio         4: row (field) matrix dimension for this PL
Cio         5: matrix dimension for this PL including i-waves
Cio         6: offset to diagonal part of g
Cio         PL are numbered 0...npl-1 here
Co Outputs
Cr Remarks
Cr   This routine generates energy-independent hamiltonian setup.
Cr       Pictorially,  the Green's function for layer ipl is dimensioned
Cr       like the following.  Here ldim0 = pgplp(4,ipl).
Cr
Cr                total column dimension = pgplp(3,ipl)
Cr          ---------------------------------------------------
Cr      p   |                   |             |               |
Cr  l   g   |                   |             |               |
Cr  d   p   |  off-diagonal     |  diagonal   |  off-diagonal |
Cr  i = l   |  g(i<ipl,ipl)     |  g(ipl,ipl) |  g(i>ipl,ipl) |
Cr  m   p   |                   |             |               |
Cr  0  (4)  |                   |             |               |
Cr          |                   |             |               |
Cr          ---------------------------------------------------
Cr                              ^
Cr                              offset to
Cr                              diagonal part
Cr                              = pgplp(6,ipl)
Cr
Cr  *For the ASA 2nd generation LMTO:
Cr   Transform pp's to alpha representation
Cu Updates
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   08 May 13 Eliminate s_array
Cu   10 Nov 11 Begin migration to f90 structures
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer pgplp(6,-1:*)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_spec)::  s_spec(*)
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_str)::   s_str
C ... Dynamically allocated local arrays
      integer, allocatable :: iprmb(:)
      real(8), allocatable :: ov(:)
C ... Local parameters
      integer i,ib1,ib2,inxsh(4),ipl,lgen3,lncol,lpgf,lsprse,nbas,nbasp,
     .  nclasp,ndim,nl,nlspcp,nofgL,nofgR,nohi,nolo,npl,nsp,nspec,pgdim

      nolo(i) = max(i-nofgL,-1)
      nohi(i) = min(i+nofgR,npl)

      nclasp = s_ctrl%nclasp
      nofgL = s_ctrl%nofgl
      nofgR = s_ctrl%nofgr
      nspec = s_ctrl%nspec
      nbasp = s_ctrl%nbasp
      nl = s_ctrl%nl
      nsp = s_ctrl%nspin
      lgen3 = s_ctrl%lgen3
      lncol = s_ctrl%lncol
      npl = s_ctrl%npl
      lpgf   = s_ctrl%lpgf(1)
      lsprse = mod(s_ctrl%lpgf(2),10)
      nlspcp = nl * nsp * nclasp

      if (lgen3 /= 0) then
        call rx('supgh not ready for NMTO')
      endif

C --- Make matrix dimensions for each PL -> pgplp(3..6) ---
C     Vector of downfolding switches
      i = s_spec(1)%lmxa
      s_spec(1)%lmxa = nl-1
      s_spec(1)%lmxa = i
      ndim = nbasp * nl**2
      allocate(iprmb(ndim))
      call pshpr(00)
      do  ipl = -1, npl
        call gtibpl(ipl,npl,pgplp,ib1,ib2)
        call makidx(nl,1,ib1,ib2,0,s_spec,s_ctrl%ips,-1,iprmb,
     .    inxsh)
        pgplp(4,ipl) = inxsh(1)
        pgplp(5,ipl) = inxsh(2)
      enddo
C ... pgplp(3) is dim. of off-diag GF connecting ipl
      do  ipl = -1, npl
        pgplp(3,ipl) = pgdim(0,ipl,npl,nolo(ipl),nohi(ipl),pgplp)
        pgplp(6,ipl) = pgdim(2,ipl,npl,nolo(ipl),ipl,pgplp)
      enddo
      deallocate(iprmb)
      call poppr

C ... Setup for sparse matrix case
      if ((lpgf == 1.or.lpgf >= 5) .and. lsprse /= 0) then
        call pglusu(s_ham,nbas,npl,pgplp,s_str%npr,s_str%iax,s_ham%offH)
C      elseif (lsprse /= 0) then
C        call rxi('supgh: sparse mode incompatible with lpgf=',lpgf)
      endif

C --- Second-generation LMTO, ASA ---
      allocate(ov(nlspcp))
      call pptrns(0,nl,s_ctrl%ipc,nclasp,nsp,s_str%alph,
     .  size(s_ctrl%ipc),s_pot%pp,ov)
      deallocate(ov)

      end
