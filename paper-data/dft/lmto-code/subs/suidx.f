      subroutine suidx(nkaph,mode,nspec,s_spec)
C- Setup for idxdn
C ----------------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: lmxb lmxa idxdn pz orbp name
Ci         :Note: input pz serves both as a boundary condition for any
Ci         :local orbitals defined, and also the 10's digit flags how
Ci         :each local orbital is to be included in the hamiltonian.
Ci         :10s digit on input:
Ci         :0 value and slope of local orbital=0 at MT boundary
Ci         :1 a smooth Hankel tail is attached.  The orbital is included
Ci         :  in the basis
Ci         :2 a smooth Hankel tail is attached.  The orbital is not
Ci         :  included in the basis; instead the coupling between the
Ci         :  local orbital and remaining states in the system is
Ci         :  included perturbatively.
Co     Stored:    idxdn orbp
Co     Allocated: *
Cio    Elts Passed:*
Cio    Passed to: uspecb
Ci Inputs
Ci   nkaph :The maximum number of radial functions centered at
Ci         :particular R and l channel used in the lmto basis.
Ci         :This includes both the number of envelope functions
Ci         :and the number of local orbitals
Ci   mode  :1s digit
Ci         :  0 lmxb cutoff (e.g. ASA 2nd generation hamiltonian)
Ci         :    orbitals for lmxb<l<n0 are mapped to 'neglected'
Ci         :  1 lmxb and lmxa cutoffs
Ci         :    orbitals for lmxb<l<=lmxa are mapped into 'high'
Ci         :    orbitals for lmxa<l<n0 are mapped into 'neglected'
Ci         :  add 2 to replace any idxdn=0 with idxdn=1
Ci         :  NB: for floating orbitals, lmxa=-1.
Ci         :  Thus, use lmxa -> max(lmxa,lmxb)
Ci         :10s digit
Ci         :  1 lmf style hamiltonian:
Ci         :    orbitals with positive rsmh are defined as being
Ci         :    in the basis.
Cl Local variables
Cl    iloc :0 if basis has no local orbitals
Cl         :1 if basis has local orbitals
Cr Remarks
Cr  Each basis function is labelled by an l quantum number and a species
Cr  index (a function of this type is centered at every site corresponding
Cr  to the species index).  Also, there may be more than one kind of basis
Cr  function per site and l, but there can be at most nkaph of such kinds.
Cr  Thus the maximum possible number of function types associated with a
Cr  particular site is nkaph*lmxa.
Cr
Cr  Array idxdn(1..lmxb+1,1:nkaph) keeps track of how functions are used
Cr  in the construction of the basis.  For a particular 0<=l<=lmxb and
Cr  1<=ik<nkaph, idxdn takes one of the following:
Cr   value   Signifies
Cr     0     Role of this orbital has not yet been determined
Cr     1     Orbital is included as "active" orbitals, which means
Cr           they are included in the hamiltonian and diagonalized
Cr     2     Orbital is treated as an "intermediate" orbital, which
Cr           means it is downfolded and included in a perturbative way
Cr     3     Orbital is treated as an "high" orbital, which means
Cr           means it is included in augmentation for tails of other
Cr           orbitals, but is otherwise not part of the basis.
Cr     4     Orbital is neglected
Cr    10     Orbital is a local orbital whose value and slope are
Cr           constructed to be zero at the MT boundary.
Cr           It is included in the basis.
Cr    11     Orbital is a local orbital with a smooth Hankel tail
Cr           and it is included in the basis.
Cr    12     Orbital is a local orbital with a smooth Hankel tail
Cr           It is incorporated perturbatively to the hamiltonian
Cr           and is not assembled as part of the hamiltonian matrix
Cr
Cr   Routine makidx uses idxdn to construct the ordering  of orbitals
Cr   in the hamiltonian matrix.  See description of
Cr   offH and iprmb in makidx.
Cb Bugs
Cb   Perhaps iloc should be species-dependent.  Implies nkaph
Cb   channel is only local orbital for species which have them.
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   24 Jun 04 Extended definition of local orbitals.  New arg list.
Cu   10 Apr 02 Redimensioned eh,rsmh to accomodate larger lmax
Cu   24 Aug 01 Extended to local orbitals.  Altered argument list.
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nkaph,nspec
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
C ... Local parameters
      logical lmf
      integer is,n0,nkap0,lmxb,lmxa,iloc
      parameter (n0=10,nkap0=4)
      integer idxdn(n0,nkap0),lp1,mode0,ik
      double precision orbp(n0,2,nkap0),pz(n0,2),dum

      mode0 = mod(mode,10)
      lmf = mod(mode/10,10) == 1
      call sanrg(.true.,mode0,0,3,'suidx:','1s digit mode')
      call sanrg(.true.,mod(mode/10,10),0,1,'suidx:','10s digit mode')
C     Determine whether basis has any local orbitals
      call uspecb(0,0,s_spec,1,nspec,dum,dum,dum,iloc)
      if (iloc/10 /= 0) then
        iloc = 1
      endif

      do  is = 1, nspec
        lmxb = s_spec(is)%lmxb
        lmxa = s_spec(is)%lmxa
        idxdn = s_spec(is)%idxdn
        pz = s_spec(is)%pz

        do  ik = 1, nkap0

C         Map all orbitals above max(lmxb,lmxa) into `neglected'
          call ivset(idxdn(1,ik),max(lmxb,lmxa)+2,n0,4)

C         mode0=0: Map any orbital above lmxb into `neglected'
          if (mod(mode0,2) == 0) then
            call ivset(idxdn(1,ik),lmxb+2,n0,4)

C         Map any orbital between lmxb and lmxa into `high'
          else
            call ivset(idxdn(1,ik),lmxb+2,lmxa+1,3)
          endif

C         Mapping for local orbitals
          if (ik == nkaph .and. iloc == 1) then
            call ivset(idxdn(1,ik),1,n0,4)
            do  lp1  = 1, lmxb+1

              if (pz(lp1,1) /= 0) then
               call pz2idx(1,pz(lp1,1),idxdn(lp1,ik))
C              if (mod(pz(lp1,1),100d0) /= 0) idxdn(lp1,ik) = 10
C              if (mod(pz(lp1,1),100d0) > 10) idxdn(lp1,ik) = 11
C              if (mod(pz(lp1,1),100d0) > 20) idxdn(lp1,ik) = 12
C              if (mod(pz(lp1,1),100d0) > 10) pz(lp1,1)=pz(lp1,1)-10
C              if (mod(pz(lp1,1),100d0) > 10) pz(lp1,1)=pz(lp1,1)-10
C              if (mod(pz(lp1,1),100d0) > 10) call rx1('suidx: bad value for pz',pz(lp1,1)+20)
              endif
            enddo
          endif

C         Map all orbitals for ik>nkaph into `neglected'
          if (ik > nkaph) call ivset(idxdn(1,ik),1,n0,4)

C     ... Automatic downfolding turned off ... set idxdn=0 to idxdn=1
          if (mode0 >= 2) then
            do  lp1  = 1, n0
              if (idxdn(lp1,ik) == 0) idxdn(lp1,ik) = 1
            enddo
          endif

C     ... lmf-specific
          if (lmf .and. ik <= nkaph-iloc) then
            orbp = s_spec(is)%orbp
            do  lp1 = 1, lmxb+1
              if (orbp(lp1,1,ik) < 0) idxdn(lp1,ik) = 3
C             if (orbp(lp1,1,ik) == 0) idxdn(lp1,ik) = max(idxdn(lp1,ik),2)
              if (orbp(lp1,1,ik) == 0) idxdn(lp1,ik) = 3
            enddo
          endif

        enddo

        s_spec(is)%idxdn = idxdn
      enddo

      end

      subroutine pz2idx(mode,pz,idxdn)
C- Converts information about local orbital between pz and idxdn
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 do nothing; just return
Ci         :1 use 10's digit of pz to set idxdn
Ci         :2 use idxdn to set 10's digit of pz
Ci         :Add 10 clear 10's digit of pz before exit
Cio Inputs/Outputs
Cio   pz   :10s digit of pz contains the following information
Cio        :0 value and slope of local orbital=0 at MT boundary
Cio        :1 a smooth Hankel tail is attached.  The orbital is included
Cio        :  in the basis
Cio        :2 a smooth Hankel tail is attached.  The orbital is not
Cio        :  included in the basis; instead the coupling between the
Cio        :  local orbital and remaining states in the system is
Cio        :  included perturbatively.
Cio        :mode 2 or 12: pz input: sets idxdn from 10's digit of pz
Cio        :mode 1 pz output: 10's digit of pz from idxdn
Cio        :mode >=10 10's digit of pz cleared on exit
Cio   idxdn:contains the same kind information as the 10s digit of pz.
Cio        :<10 equivalent to pz=0
Cio        :10  equivalent to 10's digit pz=0
Cio        :11  equivalent to 10's digit pz=1
Cio        :12  equivalent to 10's digit pz=2
Cio        :mode 2 or 12: idxdn output: sets idxdn from 10's digit of pz
Cio        :mode 1 idxdn input: 10's digit of pz from idxdn
Cr Remarks
Cr   idxdn keeps information about how local orbitals are used because
Cr    it controls ordering of hamiltonian orbitals.  See makidx.
Cu Updates
Cu   25 Jun 04 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,idxdn
      double precision pz
C ... Local parameters
      double precision pzl

      call sanrg(.true.,mod(mode,10),0,2,'pz2idx','1s digit mode')
      call sanrg(.true.,mode/10,0,1,'pz2idx','10s digit mode')

      if (mod(mode,10) == 1) then
        pzl = mod(pz,100d0)
        call fsanrg(pzl,0d0,28d0,0d0,'pz2idx','pz',.true.)
        if (pzl /= 0)  idxdn = 10
        if (pzl >= 10) idxdn = 11
        if (pzl >= 20) idxdn = 12
      elseif (mod(mode,10) == 2) then
        call rx('check pz2idx')

        call sanrg(.true.,mod(idxdn,10),0,2,'pz2idx','1s digit idxdn')
        call sanrg(.true.,idxdn/10,0,1,'pz2idx','10s digit idxdn')
        if (idxdn <  10) pz =  0
        if (idxdn == 10) pz = mod(pz,10d0) +  0 + 100*mod(pz,100d0)
        if (idxdn == 11) pz = mod(pz,10d0) + 10 + 100*mod(pz,100d0)
        if (idxdn == 12) pz = mod(pz,10d0) + 20 + 100*mod(pz,100d0)
      endif

      if (mode >= 10) pz = mod(pz,10d0)

      end
