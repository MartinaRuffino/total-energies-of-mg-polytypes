      subroutine shoist(mode,istab,nbas,ag,g,ng,ngt)
C- Show site permutation table istab
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 print to std out istab only
Ci         :1 also print group operations
Ci   istab :table of site permutations for each group op (mksym.f,symtbl.f)
Ci   nbas  :size of basis
Ci   ag    :translation part of space group
Ci   g     :point group operations
Ci   ng    :number of group operations
Ci   ngt   :if nonzero, extra group from AFM symmetry
Co Outputs
Co   Information about group operations is printed out
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   30 Jan 13 New mode
Cu   04 Jan 10 New ngt
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer, intent(in) :: mode,ng,ngt,nbas
      real(8), intent(in) :: g(3,3,ng),ag(3,ng)
      integer, intent(in) :: istab(nbas,ng)
C ... Local parameters
      integer i,ig,stdo
      character sg*70
      procedure(integer) :: nglob

      if (mode /= 0) then
        call info0(10,1,0,' SHOIST: space group operations'//
     .    ' and site permutation table')
        do  ig = 1, ng
          call asymop(g(1,1,ig),ag(1,ig),' ',sg)
          call info2(10,0,0,' ig %,2i  '//sg//'%a',ig,0)
          do  i = 1, 3
            call info2(10,0,0,'%3;12,7D',g(i,:,ig),0)
          enddo
        enddo
      endif

      if (mode == 0)
     .  call info0(10,1,0,' SHOIST: site permutation table')
      if (ngt == 0) then
        call info0(10,0,0,'  ib  istab ...')
      elseif (ng <= 1) then
        call info2(10,0,0,'  ib   E%npAFM',ng*3+6,0)
      else
        call info2(10,0,0,'  ib  istab ...%npAFM',ng*3+6,0)
      endif
      stdo = nglob('stdo')

      do  i = 1, nbas
        write(stdo,"(i4,':',48(x,i0))") i, (istab(i,ig), ig=1,ng+ngt)
      enddo

      end
