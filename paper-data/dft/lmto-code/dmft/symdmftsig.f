      subroutine symdmftsig(s_lat,s_dmft,nomg,sig)
C- Symmetrize a DMFT self-energy (sigma or delta)
C ----------------------------------------------------------------------
Cio Structures
Cio  s_dmft :struct for dmft interface; see structures.h
Ci     Elts read:  l icix
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:icix sig
Cio    Passed to:  ineqcix sigp2sigij
Ci Inputs
Cl   ldcix :dimension of largest cix block, for dimensioning
Ci   ncix  :Number of correlated blocks
Ci   nomg  :dimensioning parameter, number of frequency points
Cio Inputs/Outputs
Cio  sig   :self-energy, stored in compressed format
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   06 Aug 18 (MvS) First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... For structures
!      include 'structures.h'
      type(str_dmft)::  s_dmft
      type(str_lat)::   s_lat
C ... Passed parameters
      integer, intent(in) :: nomg
      real(8) sig(nomg,*)
C ... Dynamically allocated local arrays
      complex(8), allocatable :: hrs(:,:,:,:),hrss(:,:,:,:) !,sigi(:,:,:),sigs(:,:,:)
      integer, allocatable :: iax(:,:),ibp(:),ntab(:),iprmb(:,:),istb2(:,:)
C ... Local parameters
      integer, parameter :: niax = 10, NULLI = -99999
      integer ldcix,ndhrs,isp,iomg,cix,icix,cixi,nm,nlm,i,j,nsp,ncix,ib,low,high,
     .  mxorb,nkapsav,nlsav,nbas,nlloc,nsafm
      real(8) xx
      procedure(integer) :: nglob,ll
      procedure(real(8)) :: dglob

C ... Setup
      nsp = nglob('nsp'); nbas = nglob('nbas')
      ncix = s_dmft%ncix; ldcix = maxval(s_dmft%ndim)
C     Order cix blocks (ibp(:)) and find r.s. dimensioning parameter ndhrs
      allocate(ibp(ncix))
      call ivheap(1,ncix,s_dmft%ib,ibp,101)
      ndhrs = 0
      do  cix = 1, ncix
        icix = iabs(s_dmft%icix(cix))
        nlm = (s_dmft%l(icix)+1)**2
        ndhrs = max(ndhrs,nlm)  ! maximum value of nlm, to dimension hrs
        ib = s_dmft%ib(ibp(cix))
        if (cix == ncix) cycle
      enddo
      mxorb = ndhrs; nlloc = ll(ndhrs)+1

C ... Make iprmb, iax, ntab
      i = 0                     ! orbital index
      low = 0; high = mxorb*ncix ! counters for partitioning low and high blocks
      allocate(iprmb(mxorb,ncix),iax(niax,ncix),ntab(0:nbas))
      call iinit(iax,niax*ncix); call iinit(ntab,ncix+1)
      do  cixi = 1, ncix        ! Loop cix in permuted order (cixi)
        cix = ibp(cixi)         ! True cix
        icix = iabs(s_dmft%icix(cix))

C       Make iprmb
        nm = (s_dmft%l(icix))**2
        nlm = (s_dmft%l(icix)+1)**2
        do  j = 1, mxorb
          if (j>nm .and. j<=nlm) then
            low = low+1
            iprmb(j,cixi) = low
          else
            high = high+1
            iprmb(j,cixi) = high
          endif
        enddo

C       Make ntab and iax ... hamiltonian is single-site
        ib  = s_dmft%ib(ibp(cixi))
        nm = 2*s_dmft%l(icix)+1
        nlm = (s_dmft%l(icix)+1)**2
        if (nlm > ndhrs) call rx('bug in symdmftsig')

        iax(1,cixi) = s_dmft%ib(cix)
        iax(2,cixi) = s_dmft%ib(cix)
        ntab(ib) = cixi        ! This will be written more than once if another cix block same ib

      enddo

C --- For each frequency, expand sig to matrix form and symmetrize it ---
C     Save nl and nkaph, set to largest DMFT l+1 and 1, respectively
      nlsav = nglob('nl'); xx = dglob('nl',dble(nlloc),1)
      nkapsav = nglob('nkaph'); xx = dglob('nkaph',1d0,1)
C     Make istb2
      nsafm = iabs(s_lat%nsafm)
      allocate(istb2(nbas,max(s_lat%nsgrp,nsafm)))
      call istbpm(s_lat%istab,nbas,s_lat%nsgrp,nsafm,istb2)
C     Reset static nl and nkaph in orbl, needed by rsmsym
      call setorblp(nlloc,1)

      allocate(hrs(ndhrs,ndhrs,nsp,ncix),hrss(ndhrs,ndhrs,nsp,ncix))
      do  iomg = 1, nomg

C   ... Decompress sig for this block; poke into hrs
        do  cixi = 1, ncix      ! Loop over cix, in permuted order (cixi)
          cix = ibp(cixi)       ! True cix
          ib  = s_dmft%ib(ibp(cixi))
          icix = iabs(s_dmft%icix(cix))  ! Used in printout only
          do  isp = 1, nsp
            call sigp2sigij(s_dmft,1,nomg,nsp,iomg,cix,ndhrs,sig,sig,hrs(1,1,isp,cixi))
C            call yprmi('sigma, %s,cixi,cix,icix=%3i ib,isp=%2i',[cixi,cix,icix],[ib,isp],3,
C     .        hrs(1,1,isp,cixi),0,ndhrs,ndhrs,ndhrs)
          enddo
        enddo                   ! loop over cix blocks

C   ... Symmetrize
        i = 1 + 10*0 + 100 + 100000*0
        if (iomg>1) call pshpr(10)
        call rsmsym(i,s_lat%plat,mxorb,iprmb,low,nbas,s_lat%pos,nlloc,nsp,1,ncix,
     .    ntab,iax,s_lat%symgr,istb2,s_lat%nsgrp,0,ndhrs,hrs,hrss)
        if (iomg>1) call poppr

C   ... Store sig in compressed form
        do  cixi = 1, ncix      ! Loop over cix, in permuted order (cixi)
          cix = ibp(cixi)       ! True cix
          ib  = s_dmft%ib(ibp(cixi))
          icix = iabs(s_dmft%icix(cix))  ! Used in printout only
          do  isp = 1, nsp
C            call yprmi('sym sigma, %s,cixi,cix,icix=%3i ib,isp=%2i',[cixi,cix,icix],[ib,isp],3,
C     .        hrss(1,1,isp,cixi),0,ndhrs,ndhrs,ndhrs)
            call sigp2sigij(s_dmft,3,nomg,nsp,iomg,cix,ndhrs,sig,sig,hrss(1,1,isp,cixi))
          enddo
        enddo                   ! loop over cix blocks

      enddo                     ! omega loop


      call setorblp(nlsav,nkapsav)
      deallocate(ibp,iprmb,iax,ntab,istb2,hrs,hrss)

      stop
      end subroutine symdmftsig

C#ifdefC TEST
C      subroutine fmain
C
C      use structures
C      implicit none
CC ... For structures
C!      include 'structures.h'
C      type(str_dmft)::  s_dmft
C      integer ncix,nicix(5),cixi,cix,icix
C
C      ncix = 5
C      allocate(s_dmft%icix(ncix))
C
C      s_dmft%icix(1) = 1
C      s_dmft%icix(2) = 3
C      s_dmft%icix(3) = 1
C      s_dmft%icix(4) = 1
C      s_dmft%icix(5) = 3
C      print *, s_dmft%icix
C
C      call ineqcix(s_dmft,ncix,nicix)
C
CCr  To loop over all equivalent cix group together use a looping contstruct like this:
C      do  cixi = 1, ncix        ! Loop over all inequivalent cixi
C        if (nicix(cixi) >= 0) cycle ! nicix(:)<0 => first of equivalent cix
C        icix = iabs(s_dmft%icix(cixi)) ! inequivalent cix index
C        do  cix = 1, ncix       ! Loop over cix equivalent to cixi
C          if (iabs(s_dmft%icix(cix)) /= icix) cycle ! Skip blocks not equivalent to this one
C          print *, 'body of calculation here',cix,cixi
C        enddo
C      enddo
C
C      end
C#endif
