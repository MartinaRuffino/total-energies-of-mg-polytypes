      subroutine hsubblock(mode,s_spec,s_site,iblst,nblst,nlmto,mxbldim,
     .  iprmb,iblock,nenv,blocklist,nblocklist)
C- Count number of basis functions, ham subblocks or subblock site list
C ----------------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa lmxb pz name orbp
Co     Stored:     orbp
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  uspecb
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs
Ci   mode  :1s digit
Ci         :0 return number of envelope functions, and number of
Ci         :  hamiltonian subblocks
Ci         :10s digit
Ci         :0 only envelope functions are counted
Ci         :1 all functions associated with a site are counted
Ci   iblst :list of sites for which to make coffs.
Ci         :iblst(1)=0 => list is 1...nblst
Ci   nblst :number of sites
Ci   nlmto :dimension of lmto component of basis
Ci   mxbldim:subblock dimension cannot exceed dimension mxbldim
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Cio Inputs/Outputs
Cio iblock :(output,mode=0) number of blocks to partition site list
Cio        :(input, mode=1) get data for block iblock
Co Outputs
Co    nenv :(mode=0) number of orbitals in full site list
Co         :(mode=1) number of orbitals in subblock
Co  nblocklist : number of sites in iblock
Co   blocklist : list of sites in iblock
Cl Local variables
Cl   norb  :number of orbital types for current ib
Cl   offhi :hamiltonian offset for site ib
Cl   offl  :offl(norb) offset in h to this block of orbitals
Cl   blks  :blks(iorb) = size of contiguous block of orbitals (gtbsl1)
Cl   itop  :number of orbitals in current block
Cr Remarks
Cr   This routine is designed for partitioning hamiltonian into
Cr   subblocks, for parallelization and to conserve memory.
Cr
Cr   Call with mode=0 to determine the total dimension of the
Cr   the portion of the hamiltonian treated by conditions specified,
Cr   and to find nblock=number of subblocks (returned in iblock)
Cr
Cr   Make loops over subblocks of this type:
Cr     do  iblock = 1, nblock
Cr       call hsubblock(1,...,iblock,..,nenv,blocklist,nblocklist)
Cr      ...
Cr     enddo
Cr    nenv is the subblock dimension
Cr    blocklist(1:nblocklist) is the list of sites
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   10 Jul 10 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nblst,iblst(nblst),nlmto,iprmb(*),mxbldim,
     .  nenv,iblock,nblocklist,blocklist(*)
      integer n0,nkap0
      parameter (n0=10,nkap0=4)
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
C ... Local parameters
      integer ib,ibl,is,nlm1,nlm2,l,iorb,nenvi,nkap,offhi,mode0,mode1,
     .  ndim,itop,blocknow
      integer norb,lh(nkap0),ltab(n0*nkap0),ktab(n0*nkap0),
     .  offl(n0*nkap0),blks(n0*nkap0),ntab(n0*nkap0)
      double precision eh(n0,nkap0),rsmh(n0,nkap0)

      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
      nenv = 0
      if (mode0 == 1) nblocklist = 0

      itop = 0
      blocknow = 1
      do  ibl = 1, nblst
        if (iblst(1) == 0) then
          ib = ibl
        else
          ib = iblst(ibl)
        endif
        nenvi = 0
        is = s_site(ib)%spec
C       List of orbitals, their l- and k- indices, and ham offsets
        call orbl(ib,0,nlmto,iprmb,norb,ltab,ktab,offhi,offl,ndim)
        if (mode1 == 1) then
          nenvi = ndim
        else
C         Block routines into groups with common (e,rsm)
          call uspecb(0,2,s_spec,is,is,lh,rsmh,eh,nkap)
          call gtbsl1(7+16,norb,ltab,ktab,rsmh,eh,ntab,blks)
          do  iorb = 1, norb
            if (blks(iorb) /= 0) then
              l    = ltab(iorb)
              nlm1 = l**2+1
              nlm2 = nlm1 + blks(iorb)-1
              nenvi = nenvi + max(nlm2-nlm1+1,0)
            endif
          enddo
        endif
        itop = itop + nenvi
C       Special case 1st site alone exceeds limit
        if (ibl == 1 .and. nenvi > mxbldim) then
          blocknow = 0
        endif
C       Roll over to next block
        if (itop > mxbldim) then
          blocknow = blocknow+1
          itop = nenvi
        endif
        if (mode0 == 0) then
          nenv = nenv + nenvi
C         print '(100i4)',  ib, nenvi, itop, blocknow
        elseif (mode0 == 1) then
          if (iblock < blocknow) then
            return
          elseif (iblock == blocknow) then
            nenv = nenv + nenvi
            nblocklist = nblocklist+1
            blocklist(nblocklist) = ib
C             print '(3i3,2x,100i3)',
C     .        ib, nenv, blocknow,blocklist(1:nblocklist)
          endif
        endif

      enddo
      if (mode0 == 0) then
        iblock = blocknow
      endif

      end
