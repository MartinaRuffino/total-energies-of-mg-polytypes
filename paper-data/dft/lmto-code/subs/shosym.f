      subroutine shosym(mode,nbas,nsgrp,ngt,plat,s_sym)
C- Show site permutation table istab
C ----------------------------------------------------------------------
Cio Structures
Cio  s_sym  :struct containing space group info; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   mode  :0 print nothing
Ci         :1 print print group operations
Ci         :2 print mib
C          :1+2 may be taken in combination
Ci         :4 print everything
Ci   nbas  :size of basis
Ci   nsgrp :number of space group operations
Ci   ngt   :if nonzero, extra group from AFM symmetry (not used now)
Ci   plat  :lattice vectors (used only if mode=4)
Co Outputs
Co   Information about group operations is printed to stdout
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   30 Jan 13 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nsgrp,ngt,nbas
      double precision plat(3,3)
C ... For structures
!      include 'structures.h'
      type(str_symops)::  s_sym(nsgrp)
C ... Local parameters
      integer i,ig,stdo,ib,PRT0,iprint,ixv(3)
      double precision xv(3),qlat(3,3)
      parameter (PRT0=10)
      character sg*70
      procedure(integer) :: nglob

      if (mode == 0 .or. iprint() < PRT0) return
      stdo = nglob('stdo')
      call info0(PRT0,1,-1,' SHOSYM:')

C ... Print out site permutation table and translation
      if (mod(mode/4,2) /= 0 .and. associated(s_sym(1)%mib)) then
        call dinv33(plat,1,qlat,xv)
        call info0(PRT0,0,0,'  space group operations')
        do  ig = 1, nsgrp
          call asymop(s_sym(ig)%rg,s_sym(ig)%ag,' ',sg)
          call info2(PRT0,0,0,' ig %,2i  op: '//sg//'%a  inverse %,2i',
     .      ig,s_sym(ig)%inv)
          do  i = 1, 3
            call info2(PRT0,0,0,'%3;12,7D',s_sym(ig)%rg(i,:),0)
          enddo
          write(stdo,"(' ib  mib',T25,'tib',T47,'x Plat')")
          do  ib = 1, nbas
            call dmpy31(10,qlat,s_sym(ig)%tib(1,ib),xv)
            ixv = nint(xv)
            write(stdo,3) ib, s_sym(ig)%mib(ib),s_sym(ig)%tib(:,ib),ixv
    3       format(i3,i4,3f12.6,1x,3i3)
          enddo
          write(stdo,"(1x)")
        enddo
        return
      endif

C ... Print out space group info
      if (mod(mode,2) /= 0 .or. mod(mode/4,2) /= 0) then
        call info0(PRT0,0,0,'  space group operations')
        do  ig = 1, nsgrp
          call asymop(s_sym(ig)%rg,s_sym(ig)%ag,' ',sg)
          call info2(PRT0,0,0,' ig %,2i  inv %,2i  '//sg//'%a',
     .      ig,s_sym(ig)%inv)
          do  i = 1, 3
            call info2(PRT0,0,0,'%3;12,7D',s_sym(ig)%rg(i,:),0)
          enddo
        enddo
        call info0(PRT0,0,0,' ')
      endif

C ... Print out site permutation table
      if (mod(mode/2,2) /= 0) then
        call info0(PRT0,0,0,'  Site permutation table')
        if (ngt == 0) then
          call info0(PRT0,0,0,'  ib  mib ...')
        elseif (nsgrp <= 1) then
          call info2(PRT0,0,0,'  ib   E%npAFM',nsgrp*3+6,0)
        else
          call info2(PRT0,0,0,'  ib  mib ...%npAFM',nsgrp*3+6,0)
        endif
        do  i = 1, nbas
          write(stdo,4) i, (s_sym(ig)%mib(i), ig=1,nsgrp+ngt)
    4     format(i4,':',48i3)
        enddo
      endif

      end
