      subroutine symtbl(mode,tol,nbas,ipc,pos,g,ag,ng,qlat,istab,trtab)
C- Make site permutation table for each symop; check classes
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1st digit (see Remarks)
Ci         :0  site ib is transformed into istab(ib,ig) by grp op ig
Ci         :1  site istab(i,ig) is transformed into site i by grp op ig
Ci         :10s digit
Ci         :1  check atom classes
Ci         :2  generate trtab.  Also requires 1s digit mode=0
Ci         :   bits 1 and 2 may be combined
Ci         :100s digit
Ci         :1  Symop cannot map atom to itself
Ci         :1  Allow symop to map atom to itself, but print warning
Ci   tol   :tol for which atoms are considered to be at the same site
Ci         :use 0 for symtbl to pick internal default
Ci   nbas  :size of basis
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   pos   :pos(i,j) are Cartesian coordinates of jth atom in basis
Ci   g     :point group operations
Ci   ag    :translation part of space group
Ci   ng    :number of group operations
Ci   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
Co Outputs
Co   istab :table of site permutations for each group op; see mode
Co   trtab :(made only for 1s digit mode 2)
Co         :table of lattice translation vectors that separate rotated
Co         :atom and equivalent atom after rotation (see Remarks)
Cr Remarks
Cr   For group operation ig (mode 0 or 2)
Cr   Site ib is mapped into site jb=istab(ib,ig) + lattice vector
Cr   The lattice vector is returned in trtab (mode 2).  Thus:
Cr      trtab(1..3,ib,ig) = [g(ig) pos(1..3,ib) + ag(ig)] - pos(1..3,jb))
Cu Updates
Cu   07 Aug 17 Revise for new AFM symmetrization
Cu   17 Jun 13 Replace f77 pointers with f90 ones
Cu   04 Jan 10 New 100s digit mode
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer nbas,ng,mode
      integer ipc(nbas),istab(nbas,ng)
      double precision pos(3,nbas),g(3,3,ng),ag(3,ng),qlat(9),tol,trtab(3,nbas,ng)
C Local variables
      logical ltrtab
      integer ib,ic,ig,jb,jc,mode0,mode1,mode2
      double precision tol0,tol1,dlat(3)
      parameter (tol0=1d-5)
      integer, allocatable :: iwk(:)

      if (ng == 0) return
      mode0 = mod(mode,10)
      ltrtab = mod(mode/10,10) == 2
      mode1 = mod(mod(mode/10,10),2)
      mode2 = mod(mode/100,10)
      tol1 = tol
      if (tol == 0) tol1 = tol0

C --- Make site transformation table ---
C     if (mode2 == 0) then
      do  ig = 1, ng
        do  ib = 1, nbas
          call grpfnd(tol1,g,ag,ig,pos,nbas,qlat,ib,jb,dlat)
          if (jb == 0)
     .      call fexit2(-1,111,' Exit -1 SYMTBL: no map for atom '//
     .      'ib=%i, ig=%i',ib,ig)
          if (mode2 /= 0 .and. ib == jb) then
            if (mode2 == 1) call rx1('SYMTBL: site %i maps to itself',ib)
            call info2(10,0,0,'SYMTBL (warning): site %i maps to itself',ib,2)
          endif
          if (mode1 /= 0) then
            ic = ipc(ib)
            jc = ipc(jb)
            if (ic /= jc) call fexit3(-1,111,' Exit -1 SYMTBL: '//
     .        'site %i not in same class as mapped site %i, ig=%i',ib,jb,ig)
          endif
          if (mode0 == 0) then
            istab(ib,ig) = jb
            if (ltrtab) then
              call dcopy(3,dlat,1,trtab(1,ib,ig),1)
C            if (ig >= 0) then
C              print 333, ib,jb,dlat
C  333         format(' ib=',i3,' jb=',i3,' dlat=',3f11.4)
C            endif
            endif
          else
            istab(jb,ig) = ib
          endif
        enddo
      enddo

C ... Add translation to symop ng+1
C      else
C          do  ib = 1, nbas
C            g0 = (/1d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,1d0/)
C            call grpfnd(tol1,g0,ag,1,pos,nbas,qlat,ib,jb,dlat)
C            if (jb == 0)
C     .        call fexit2(-1,111,' Exit -1 SYMTBL: no map for atom '//
C     .        'ib=%i, ig=TRANS',ib,ig)
C            if (mode2 /= 0) then
C              if (ib == jb) call rx1('SYMTBL: '//
C     .          'site %i maps to itself, ig=TRANS',ib)
C            endif
C            if (mode1 /= 0) then
C              ic = ipc(ib)
C              jc = ipc(jb)
C              if (ic /= jc) call fexit3(-1,111,' Exit -1 SYMTBL: '//
C     .          'site %i not in same class as mapped site %i, ig=TRANS',
C     .          ib,jb,ig)
C            endif
C            if (mode0 == 0) then
C              istab(ib,ng+1) = jb
C            else
C              istab(jb,ng+1) = ib
C            endif
C          enddo
C      endif

C --- Check atom classes ---
      if (mode1 == 0) return
      allocate(iwk(nbas))
      do  ib = 1, nbas
        ic = ipc(ib)
        call iinit(iwk,nbas)
        do  ig = 1, ng
          iwk(istab(ib,ig)) = 1
        enddo
        do  jb = 1, nbas
          if (iwk(jb) == 1) cycle
          jc = ipc(jb)
          if (ic == jc) call fexit2(-1,111,' Exit -1 SYMTBL:  '//
     .      'site ib=%i in same class as inequivalent site jb=%i',ib,jb)
        enddo
      enddo
      deallocate(iwk)

      end

      subroutine istbpm(istab,nbas,ng,nsafm,istab2)
C- Makes inverse of istab
C ----------------------------------------------------------------------
Ci Inputs
Ci   istab :table of site permutations for each group op (mksym.f,symtbl.f)
Ci   nbas  :size of basis
Ci   ng    :number of G-vectors
Ci   nsafm :points to column with AFM symmetry, if it exists
Co Outputs
Ci   istab2:istab2(istab(ib,ig),ig) = ib for ig = 1 .. ng
Ci         :if nsafm is nonzero, then
Ci         :istab2(istab(ib,nsafm),nsafm) = ib
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   06 Aug 17
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nbas,ng,nsafm
      integer istab(nbas,*),istab2(nbas,max(ng,nsafm))
C ... Local parameters
      integer ib,ig,kg,ibp,ngp

      ngp = ng ; if (nsafm /= 0) ngp = ng+1

      do  ig = 1, ngp
        kg = ig ; if (nsafm /= 0 .and. ig == ngp) kg = nsafm
        do  ib = 1, nbas
          ibp = istab(ib,kg)
          istab2(ibp,kg) = ib
        enddo
      enddo

      end

      subroutine suafmsym(s_lat,nbas,ipc,istab,g,ag)
C- Set up arrays ipc, istab, g, ag for special AFM symmetrization
C ----------------------------------------------------------------------
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  nabc plat qlat nsgrp afmt npgrp ng
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:istab symgr ag kv gv ips0 bgv
Cio    Passed to:  symrat symsmr
Ci Inputs
Ci   nbas  :size of basis
Co Outputs
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   istab :table of site permutations for each group op (mksym.f,symtbl.f)
Ci   g     :point group operations
Ci   ag    :translation part of space group
Cr Remarks
Cr   Special AFM symmetrization (see also Remarks in mksym.f)
Cr   Performed after normal charge symmetrization if s_lat%nsafm > 0.
Cr   Combines space group op 1 (unit operator) and op s_lat%nsafm.
Cr   A check is made that under this op sites (ib,jb) combined in pairs,
Cr   so that for symmetry purposes every site belongs to a class with 2 elements.
Cr   The two sites are symmetrized with the second atom density spin-flipped.
Cr
Cr   This routine constructs the class table, and symmetry arrays specifically
Cr   for this special two-class case.
Cu Updates
Cu   07 Aug 17 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer :: ipc(nbas),istab(nbas,2)
      real(8) :: g(9,2),ag(3,2)
C ... For structures
!      include 'structures.h'
      type(str_lat) ::  s_lat
C ... Local parameters
      integer nbas,nsafm,nrclas,ib,jb
      procedure(integer) :: nglob

      nsafm = iabs(s_lat%nsafm)
      if (nsafm /= 0) then
        call iinit(ipc,nbas)
        nrclas = 0
        do  ib = 1, nbas
          jb = s_lat%istab(ib,nsafm)
          if (ib /= s_lat%istab(jb,nsafm)) call rx('suafmsym : problem with AFM symetrization')
          if (ipc(ib) /= 0) cycle
          nrclas = nrclas+1
          ipc(ib) = nrclas; ipc(jb) = nrclas
        enddo
        call icopy(nbas,s_lat%istab,1,istab,1)
        call icopy(nbas,s_lat%istab(1,nsafm),1,istab(1,2),1)
        call dcopy(9,s_lat%symgr,1,g,1)
        call dcopy(9,s_lat%symgr(1,nsafm),1,g(1,2),1)
        call dcopy(3,s_lat%ag,1,ag,1)
        call dcopy(3,s_lat%ag(1,nsafm),1,ag(1,2),1)
      endif

      end
