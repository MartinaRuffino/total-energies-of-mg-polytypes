      subroutine sblham(s_ham,s_pot,s_str,ib1,ib2,nbasx,lbloch,
     .  plat,ldim,idim,ldl,ldi,ldl2,klu,qp,sll,sil,sii)
C- Generate subblock of Bloch summed hamiltonian from file strux
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read: lham
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:iprmb
Cio    Passed to: *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read: *
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:palp
Cio    Passed to: *
Cio  s_str  :struct for parameters for screened strux; see structures.h
Ci     Elts read: npr
Co     Stored:    *
Co     Allocated: *
Cio    Elts Passed:iax s nds
Cio    Passed to: *
Ci Inputs
Ci   ib1   :the neighbor list of pairs for which to accumulate
Ci          sum starts at beginning of neighbor list for ib1
Ci   ib2   :the neighbor list of pairs for which to accumulate
Ci          sum ends at end of neighbor list for ib2
Ci   nbasx :At least the largest site which will be touched in this
Ci          subblock.  It should not exceed the number of sites
Ci          for which indexing information is available
Ci   lbloch:passed through to bloch.
Ci          1s digit concerns storage of Bloch summed hamiltonian
Ci           0: s is stored in unpacked form
Ci           1: s is stored in banded form (see Remarks)
Ci         10s digit distinguishes how complex arithmetic is handled
Ci           1: s is returned complex*16 format:
Ci              s = s(2,ldl,ldl2), with s(1,*) = real, s(2,*) = imag
Ci           2: s has real, imaginary separated by columns
Ci              s = s(ldl,2,ldl2), with s(*,1..2,*) = real..imag
Ci           3: s has real, imaginary separated
Ci              s = s(ldl,ldl2,2), with s(*,*,1..2) = real..imag
Ci        100s digit
Ci           1 if to add to s (s not initialized to zero)
Ci           2 subtract from s
Ci           3 combination of 1+2
Ci       1000s digit 1 if to convert s to spherical harmonics
Ci      10000s digit 1 ?? Suppress accumulation into Bloch sum pairs
Ci                       whose site index is out of (ib1,ib2) range
Ci   ldim  :dimension of basis in the lower set.  Orbitals for which
Ci         :0<iprm(i)<=ldim are accumulated into sll.
Ci   idim  :dimension of downfolded sii
Ci   ldl   :leading dimension of sll, for dimensioning purposes
Ci          Must have
Ci   ldi   :leading and second dimension of sii
Ci   ldl2  :second dimension of sll and sil
Ci   klu   :size of sub- and super-diagonal, if sll is in banded format
Ci   qp    :k-point
Co Outputs
Co   sll   :lower-lower bloch of Bloch summed matrix
Co   sil   :lower-intermediate bloch of Bloch summed matrix
Co   sii   :intermediate-intermediate bloch of Bloch summed matrix
Cr Remarks
Cr   This routine generates a subblock of the Bloch summed hamiltonian,
Cr   by including only orbitals associated with sites ib1..ib2.
Cr
Cr   The band form follows LAPACK band storage conventions:
Cr   s(i,j) is stored in location (kl+ku+1+i-j,j)
Cr   with kl,ku = size of sub- and super-diagonal.
Cr   Here we take kl=ku=klu.
Cu Updates
Cu   17 Jun 13 Completed replacement of f77 pointers with f90 ones
Cu   10 Nov 11 Begin migration to f90 structures
Cu   23 Aug 05 (A. Chantis) fixed for spherical harmonics.
Cu   10 Jan 02 slbham reads pf from pot->opalp instead of pot->opf
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      double precision plat(3,3),qp(3)
      integer ib1,ib2,nbasx,klu,lbloch,ldi,ldl,ldl2
C     real + imaginary storage mode
      double precision sll(ldl,ldl2,2),sil(ldi,ldl2,2),sii(ldi,ldi,2)
C     complex*16 storage mode
C     double precision sll(2,ldl,ldl2),sil(2,ldi,ldl2),sii(2,ldi,ldi)
C     real + imaginary in columns storage mode
C     double precision sll(ldl,2,ldl2),sil(ldi,2,ldl2),sii(ldi,2,ldi)
C ... For structures
!      include 'structures.h'
      type(str_ham)::   s_ham
      type(str_pot)::   s_pot
      type(str_str)::   s_str
C ... Dynamically allocated local arrays
      integer, allocatable :: iwk(:),iprm(:)
C ... Local parameters
      logical bittst
C     For bloch
      integer idim,is1,is2,ldim,nl,lham,mxorb,ii
      procedure(integer) :: nglob

      call tcn('sblham')

      mxorb = nglob('mxorb')
      nl = nglob('nl')

C --- Bloch sum strux ---
      is1 = s_str%npr(ib1)+1
      is2 = s_str%npr(ib2+1)
C ... Restrict pairs outside of block
      lham = s_ham%lham
      ii = nbasx*nl*nl
      allocate(iwk(ii),iprm(ii)); call iinit(iprm,ii)
      call icopy(ib2*nl*nl,s_ham%iprmb,1,iprm,1)
C     First call excludes all orbitals
      call oraddp(1,s_str%iax,1,0,1,nbasx,nl*nl,iwk,iprm)
C     Second call includes orbitals in range ib1..ib2
      call oraddp(0,s_str%iax,is1,is2,ib1,ib2,nl*nl,iwk,iprm)
C     A third call could exclude those we didn't want in the first place
C     (not implemented)
      if (bittst(lham,256)) then
        lbloch = lbloch + 1000
      endif
      call bloch(lbloch,qp,nl,plat,mxorb,iprm,is1,is2,s_str%iax,
     .  s_str%s,s_str%nds,1,1,ldim,ldim,idim,ldl,ldi,ldl2,klu,
     .  sll,sil,sii)

C --- ASA : Add -P to sll with (S-P) ---
      call sblhm1(lbloch,nl*nl,ib1,ib2,iprm,s_pot%palp,
     .  0,ldim,idim,ldl,ldi,ldl2,klu,sll,sil,sii)

      deallocate(iwk,iprm)

      call tcx('sblham')
      end

      subroutine sblhm1(lbloch,nl2,ib1,ib2,iprm,pfun,offp,ldim,idim,ldl,
     .  ldi,ldl2,klu,sll,sil,sii)
C- Add or subtract diagonal P from hamiltonian subblock
C ----------------------------------------------------------------------
Ci Inputs
Ci   lbloch:bits are compatible with those of routine bloch.
Ci          1s digit concerns storage of Bloch summed hamiltonian
Ci           0: s is stored in unpacked form
Ci           1: s is stored in banded form (see Remarks)
Ci         10s digit distinguishes how complex arithmetic is handled
Ci           1: s is returned complex*16 format:
Ci              s = s(2,ldl,ldl2), with s(1,*) = real, s(2,*) = imag
Ci           2: s has real, imaginary separated by columns
Ci              s = s(ldl,2,ldl2), with s(*,1..2,*) = real..imag
Ci           3: s has real, imaginary separated
Ci              s = s(ldl,ldl2,2), with s(*,*,1..2) = real..imag
Ci        100s digit
Ci           0 subtract P from s
Ci           2 add P to s
Ci   nl2   :spacing between offsets to pfun in successive sites:
Ci          offset to pfun for site ib is nl2*(ib-1).
Ci   ib1   :the neighbor list of pairs for which to accumulate
Ci          sum starts at beginning of neighbor list for ib1
Ci   ib2   :the neighbor list of pairs for which to accumulate
Ci          sum ends at end of neighbor list for ib2
Ci   iprm  :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci         :This routine uses the following part of iprm:
Ci         :iprm(nl2*(ib1-1)+1:nl2*ib2)
Ci   pfun  :vector of potential functions, in iprm order
Ci   offp  :offset to potential function, ib=1
Ci   ldim  :dimension of lower set of orbitals. See iprm, above.
Ci   idim  :dimension of intermediate set. See iprm, above.
Ci   ldl   :leading dimension of sll
Ci   ldi   :leading and second dimension of sii
Ci   ldl2  :second dimension of sll and sil
Ci   klu   :size of sub- and super-diagonal, if s stored banded form
Co Outputs
Co   sll   :lower-lower block of hamiltonian matrix
Co   sil   :lower-intermediate block of hamiltonian matrix
Co   sii   :intermediate-intermediate block of hamiltonian matrix
Cr Remarks
Cr   Complex diagonal pfun is added to sll,sil,sii
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer lbloch,nl2,ib1,ib2,ldim,idim,ldl,ldi,ldl2,klu,offp,iprm(*)
      double precision pfun(2,*)
C     sll,sil,sii assuming real + imaginary storage mode
      double precision sll(ldl,ldl2,2),sil(ldi,ldl2,2),sii(ldi,ldi,2)
C ... Local parameters
      integer nlmx
      parameter (nlmx=81)
      double precision pwk(nlmx,2,nlmx),fac
      integer ia,ib,lblchi,kcplx,ld11,ld12,ld13,ld21,ld22,ld23,lidim,
     .  offa,offb,offpi,lm,oi,ip

      if (nl2 > nlmx) call rxi('sblhm1: need nlmx at least',nl2)
      lidim = ldim+idim

C     Pick up true dimensions of sll,sil,sii from formal ones
      kcplx = mod(lbloch/10,10)
      fac = -1
      if (mod(lbloch/100,10) >= 2) fac = 1
      call cplxdm(kcplx,ldl,ldl2,ld11,ld21,oi,oi)
      call cplxdm(kcplx,ldi,ldl2,ld12,ld22,oi,oi)
      call cplxdm(kcplx,ldi,ldi,ld13,ld23,oi,oi)
C     Tell pblch1 that P is diagonal and complex
      lblchi = mod(lbloch,1000) + 11000

C --- For each site in the subblock, do ---
      do  ib = ib1, ib2

C       ia and ib are always identical
        ia = ib

C  ...  Copy pfun(iprm order) to the diagonal of matrix pwk (L order)
        offpi = nl2*(ia-1)
        do  lm = 1, nl2
          ip = iprm(offpi+lm)
          if (ip <= lidim) then
            pwk(lm,1,lm) = pfun(1,ip+offp)
            pwk(lm,2,lm) = pfun(2,ip+offp)
          endif
      enddo

C   --- Lower-lower block ---
        offb = nl2*(ib-1)
        offa = nl2*(ia-1)
        call pblch1(lblchi,nl2,offa,offb,ld11,ld21,klu,iprm,0,ldim,
     .    0,ldim,pwk,pwk,nlmx,fac,0d0,sll)

C   --- Intermediate-lower block ---
        offb = nl2*(ib-1)
        offa = nl2*(ia-1)
        if (idim == 0) cycle
        call pblch1(lblchi,nl2,offa,offb,ld12,ld22,klu,iprm,ldim,lidim,
     .    0,ldim,pwk,pwk,nlmx,fac,0d0,sil)

C   --- Intermediate-intermediate block ---
        offb = nl2*(ib-1)
        offa = nl2*(ia-1)
        if (idim == 0) cycle
        call pblch1(lblchi,nl2,offa,offb,ld13,ld23,klu,iprm,ldim,lidim,
     .    ldim,lidim,pwk,pwk,nlmx,fac,0d0,sii)

      enddo
      end

      subroutine sblhm2(lbloch,nl2,ib1,ib2,iprm,cp,isp,lidim,ldl,
     .  ldl2,klu,sll)
C- Add or subtract block-diagonal coherent P from H subblock for DLM
C ----------------------------------------------------------------------
Ci Inputs
Ci   lbloch:bits are compatible with those of routine bloch.
Ci          1s digit concerns storage of Bloch summed hamiltonian
Ci           0: s is stored in unpacked form
Ci           1: s is stored in banded form (see Remarks)
Ci         10s digit distinguishes how complex arithmetic is handled
Ci           1: s is returned complex*16 format:
Ci              s = s(2,ldl,ldl2), with s(1,*) = real, s(2,*) = imag
Ci           2: s has real, imaginary separated by columns
Ci              s = s(ldl,2,ldl2), with s(*,1..2,*) = real..imag
Ci           3: s has real, imaginary separated
Ci              s = s(ldl,ldl2,2), with s(*,*,1..2) = real..imag
Ci        100s digit
Ci           0 subtract P from s
Ci           2 add P to s
Ci   nl2   :spacing between offsets for internal work array pwk
Ci          offset to pwk for site ib is nl2*(ib-1).
Ci   ib1   :the neighbor list of pairs for which to accumulate
Ci          sum starts at beginning of neighbor list for ib1
Ci   ib2   :the neighbor list of pairs for which to accumulate
Ci          sum ends at end of neighbor list for ib2
Ci   iprm  :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   cp    :block-diagonal coherent P matrix, in iprm order
Ci   isp   :spin index for addressing cp
Ci   lidim :dimension of l+i set of orbitals. See iprm above.
Ci   ldl   :leading dimension of sll
Ci   ldl2  :second dimension of sll
Ci   klu   :size of sub- and super-diagonal, if s stored banded form
Co Outputs
Co   sll   :li-li block of hamiltonian matrix as a whole
Cr Remarks
Cr   Complex block-diagonal cp is added to sll
Cr   This routine was cloned from sblhm1 and inherits its
Cr   unjustified complexity with some simplifications.
Cr   The whole l+i block is treated as one, because that's how
Cr   this routine is called anyway
Cu Updates
Cu   09 Jan 12 (Belashchenko) Fixed errors; now calls pblch1
Cu   08 Dec 08 (P. Larson) Created for DLM
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lbloch,nl2,ib1,ib2,lidim,ldl,ldl2,klu,iprm(*),isp
      double complex cp(lidim,lidim,2)
C     sll assuming real + imaginary storage mode
      double precision sll(ldl,ldl2,2)
C ... Local parameters
      integer nlmx
      parameter (nlmx=81)
      double precision pwk(nlmx,2,nlmx),fac
      integer ia,ib,lblchi,kcplx,ld11,ld21,offa,offb,offpi
      integer lm1,lm2,ip1,ip2,oi

      if (nl2 > nlmx) call rxi('sblhm1: need nlmx at least',nl2)

C     Pick up true dimensions of sll from formal ones
      kcplx = mod(lbloch/10,10)
      fac = -1
      if (mod(lbloch/100,10) >= 2) fac = 1
      call cplxdm(kcplx,ldl,ldl2,ld11,ld21,oi,oi)
C     Tell pblch1 that CP is complex (and off-diagonal)
      lblchi = mod(lbloch,1000) + 1000

C --- For each site in the subblock, do ---
      do  ib = ib1, ib2
        ia = ib
        offpi = nl2*(ia-1)
C  ...  Copy cp(iprm order) to the diagonal of matrix pwk (L order)
        do  lm1 = 1, nl2
          do  lm2 = 1, nl2
            ip1 = iprm(offpi+lm1)
            ip2 = iprm(offpi+lm2)
            if (ip1 <= lidim .and. ip2 <= lidim) then
              pwk(lm1,1,lm2) = dble(cp(ip1,ip2,isp))
              pwk(lm1,2,lm2) = dimag(cp(ip1,ip2,isp))
            endif
          enddo
        enddo

C   ... Entire li-li block
        offb = nl2*(ib-1)
        offa = nl2*(ia-1)
        call pblch1(lblchi,nl2,offa,offb,ld11,ld21,klu,iprm,0,lidim,
     .    0,lidim,pwk,pwk,nlmx,fac,0d0,sll)
      enddo
      end
