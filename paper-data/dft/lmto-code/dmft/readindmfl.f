      subroutine readindmfl(s_site,s_dmft,s_lat,nsp)
C- Read DMFT data from indmfl file
C ----------------------------------------------------------------------
Ci Inputs
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  not used now
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_dmft
Ci     Elts read:  ncatom ndim sigind iasig ib
Co     Stored:     iproj lsigim gammac ncatom ncix nicix cf sigind
Co                 nzsigi iasig ndsig nzsig
Co     Allocated:  ib icix l qsplit ndim sigind cf nzsigi iasig
Cio    Elts passed: ib icix l qsplit ndim ncatom nicix sigind cf nzsigi iasig
Cio    Passed to:  *
Co Outputs
Cs Command-line switches
Cl Local variables
Cl   ncix         : number of correlated subblocks
Cl   nicix        : number of inequivalent correlated subblocks
Cl   ndsig        : Index to last matrix element stored in compact form
Cl   nzsig        : total number of nonzero ME in all correlated subblocks, including spin
Cl   cixmap       : cixmap(icix,1) points to a cix block this inequivalent block belongs to
Cl   maxdim       : maximum dimension of correlated subblock
Cl   imatsubara   : 0 => real axis  1 => Matsubara frequencies
Cl   gammaci      : broadening for correlated orbitals
Cl   gamma        : broadening for (noncorrelated?) orbitals ... not used now
Cl   nomega       : number of Matsubara frequencies        ... not used now
Cl   omegamin     : lower bound for range of omega (eV)    ... not used now
Cl   omegamax     : upper bound for range of omega (eV)    ... not used now
Cl   locrot       : flag specifying rotation of coordinate system
Cl                : Note: not used yet
Cl                : 0 => no rotation
Cl                :>0 => not used now
Cl                :-1 => 3x3 rotation matrix is given
Cl                : Important because only specified elements are calculated,
Cl                : typically diagonal elements (see qsplit)
Cl                : sigma should be as near diagonal as possible, and
Cl                : make use of symmetry where possible
Cl   qsplit       : specifies meaning of orbital basis in subblock
Cl                  0  average GF, non-correlated
Cl                  1  |j,mj> basis, no symmetry, except time reversal (-jz=jz)
Cl                 -1  |j,mj> basis, no symmetry, not even time reversal (-jz=jz)
Cl                  2  real harmonics basis, no symmetry, except spin (up=dn)
Cl                 -2  real harmonics basis, no symmetry, not even spin (up=dn)
Cl                  3  t2g orbitals  only
Cl                 -3  eg orbitals   only
Cl                  4  |j,mj>, only l-1/2 and l+1/2
Cl                  5  axial symmetry in real harmonics
Cl                  6  hexagonal symmetry in real harmonics
Cl                  7  cubic symmetry in real harmonics
Cl                  8  axial symmetry in real harmonics, up different than down
Cl                  9  hexagonal symmetry in real harmonics, up different than down
Cl                 10  cubic symmetry in real harmonics, up different then down
Cl                 11  |j,mj> basis, non-zero off diagonal elements
Cl                 12  real harmonics, non-zero off diagonal elements
Cl                 13  J_eff=1/2 basis for 5d ions, non-magnetic with symmetry
Cl                 14  J_eff=1/2 basis for 5d ions, no symmetry
Cr Remarks
Cr   readindmfl operates in one of two modes, mode 1 if s_dmft%ncix = 0 and mode 2 otherwise
Cr              (mode 1) reads DMFT block information from file indmfl
Cr              (mode 2) block information is passed through s_dmft, readindmfl completes s_dmft.
Cu              In this mode the following elements in s_dmft must be set:
Cu                scalars ncix nicix maxdim ncatom lsigim gammac nzsig
Cu                vectors dimensioned to ncix : ib icix nzsigi
Cu                vectors dimensioned to nicix : ndim l qsplit umode
Cu                sigind(1:maxdim,1:maxdim,nicix,1:nsp)
Cu             *Note nzsig(0:ncix) must also be passed, but readindmfl takes input values for nzsig(0:nicix),
Cu              corresponding to the local array ndsig(0:nicix).  readindmfl remakes nzsig(0:ncix).
Cu             *Note only sigind for spin 1 is passed; channels are reassigned here if nsp=2
Cr   readindmfl reads correlated subblocks first by site expecting input as for each site as follows:
Cr ... For each site with a correlated subblock read
Cr     ib   ncl  locrot    # site, number of L's, flag indicating rotation of local Cartesian coordinates
Cr ... for each member of ncl read:
Cr      l  qsplit  cix     # l, qsplit, index to which inequivalent subblock site corresponds to
Cr ... if locrot is nonzero, read information specifying rotation (depends on locrot)
Cr
Cr ... Next is read information about each inequivalent correlated subblock.  Header is:
Cr      nicix  ...         # number of correlated subblocks (line may contain other data, but it is not used)
Cr ... For each of the nicix inequivalent subblocks:
Cr      cix ndim size      # cix must be ordered 1,2,3,...  ndim should be sum(2l+1)*nsp for each l in subblock,
Cr                         # size should be inferred (but is not yet) inferred from the subsequent data
Cr      Three lines with comments
Cr      sigind matrix      # a ndim x ndim array of integers, with indices pointing to matrix elements
Cr                         # the solver will calculate
Cr      A comment line referring to "transformation matrix"
Cr      transformation matrix # specifies rotation from bath coordinates to spherical harmonics
Cb Bugs
Cb  Routine needs to account for symmetry-equivalent atoms with rotations
Cb  * l and qspl are read for each subblock but should be identified with inequivalent subblocks
Cb  * cix is read for each l within a subblock.  Should require single cix for a subblock
Cb  * No cluster DMFT
Cu Updates
Cu   26 Mar 19 (MvS) If input s_dmft%ncix is nonzero, the indmfl file is not read.  See Remarks mode 2
Cu   06 May 18 (MvS) Internally finds equivalent cix blocks when sigind indices are not consecutive
Cu   27 Feb 16 (MvS) redesign for more general cix blocks
Cu   22 Aug 15 (Pisanti) spin polarized
Cu   25 Nov 14 P Pisanti first created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer, intent(in) :: nsp
C ... For structures
!      include 'structures.h'
      type(str_site) ::  s_site(*)
      type(str_lat)  ::  s_lat
      type(str_dmft) ::  s_dmft

C ... Dynamically allocated local arrays
      integer, allocatable :: iprm(:),iasig(:,:),ndsigi(:,:),cixmap(:,:),iwk(:,:,:)
      integer, allocatable :: dmft_l(:),dmft_qsplit(:),dmft_ndim(:),dmft_ib(:),dmft_icix(:)
      real(8), allocatable :: wk(:)
      complex(8), allocatable :: zwk(:,:,:)
C ... Local parameters
      character*(1) ch
      logical lcix,lsame
      integer :: igfnd
      integer :: cix,curr,i,iat,ib,icix,icix0,ifi,imatsubara,iprojector,irenormalize,
     .  is,isigind,isp,ix,ix0,ixs,j,k,l,locrot,loro,maxdim,mxcix,
     .  ncix,ndim,ndsig,nglob,nicix,nomega,nzsig,prev,qspl,shift,stdo
      real(8) :: emin,emax,gammaci,gamma,omegamin,omegamax
      real(8), parameter :: ry2eV = 13.60569193d0


!     ... for symetry check
      integer :: ind_ff
      real(8)  :: posref(3),pos(3)
      logical(8) :: eqv=.false.
      integer, allocatable :: ipc(:),ipca(:),istab(:,:)
      real(8), allocatable :: g(:,:),ag(:,:)
      integer,allocatable :: nicixi(:)
      integer procid,master
      procedure(integer) :: ival,fopna,iinear,mpipid,iprint

      procid = mpipid(1)
      master = 0
      stdo = nglob('stdo')
      if (s_dmft%ncix > 0) then ! indmfl is not read
        ncix = s_dmft%ncix
        nicix = s_dmft%nicix
        allocate(ndsigi(0:nicix,1)) ! Number of correlated orbitals in icix block
        call icopy(nicix+1,s_dmft%nzsigi,1,ndsigi,1)
        do  cix = 1, ncix
          icix = iabs(s_dmft%icix(cix))
          s_dmft%nzsigi(cix) = s_dmft%nzsigi(cix-1) + ndsigi(icix,1)-ndsigi(icix-1,1)
        enddo
        nzsig = s_dmft%nzsig
        maxdim = maxval(s_dmft%ndim)
        allocate(iprm(maxdim*maxdim*nicix*nsp)) ! Work array to make permutation table
        goto 100
      endif

      call info0(10,1,-1,' Reading indmfl file ...')

      if (procid == master) then
        ifi = fopna('indmfl',-1,0) ! Open file as OLD => generate error if missing
        read(ifi,*) emin,emax,irenormalize,iprojector ! energy in eV
      endif
      call mpibc1(emin,1,4,.false.,'','')
      call mpibc1(emax,1,4,.false.,'','')
      call mpibc1(irenormalize,1,2,.false.,'','')
      call mpibc1(iprojector,1,2,.false.,'','')

      s_dmft%iproj = iprojector
!     if irenormalize=1 -> qrenormalize = true, else false (renormalization of projector)
      emin = emin/ry2ev
      emax = emax/ry2ev !these energies have to be connected to the bands

      s_dmft%lsigim = .false.
      if (procid == master) then
        read(ifi,*) imatsubara, gammaci, gamma, nomega, omegamin, omegamax ! See Local variables
      endif
      call mpibc1(imatsubara,1,2,.false.,'','')
      call mpibc1(gammaci,1,4,.false.,'','')
!      call mpibc1(gamma,1,4,.false.,'','')  Not used now
!      call mpibc1(nomega,1,2,.false.,'','') Not used now
!      call mpibc1(omegamin,1,4,.false.,'','') Not used now
!      call mpibc1(omegamax,1,4,.false.,'','') Not used now

      if (imatsubara /= 0) s_dmft%lsigim = .true.
      s_dmft%gammac = gammaci/ry2ev
      if (procid == master) then
        read(ifi,*) i           ! Number of correlated atoms.  In future, should be inferred from later data
      endif
      call mpibc1(i,1,2,.false.,'','')
      s_dmft%ncatom = i
      mxcix = 10*s_dmft%ncatom  ! Dimensioning of temporary variables

      allocate(dmft_l(mxcix),dmft_qsplit(mxcix),dmft_ndim(mxcix),dmft_ib(mxcix),dmft_icix(mxcix))

C --- Count correlated subblocks associated with all atoms and read initial data ---
C     make ncix = total number of subblocks and read for each ncix : dmft_ib, dmft_icix
C     make nicix = number of inequivalent subblocks and read for each nicix : qsplit, l, ndim
C     Read possible local rotation of Cartesian coordinates for each site.
C     Determine maxdim = size of largest block
      maxdim = 0                ! Maximum rank of any subblock (to be computed in this loop)
      lcix = .false.            ! Set to true when a correlated orbital is encountered
      dmft_qsplit = 0           ! qsplit read by site, but stored only for inequivalent blocks (see Bugs above)
      dmft_l = 0                ! l read by site, but stored only for inequivalent blocks (see Bugs above)
      dmft_ndim = 0             ! ? maybe not needed
      ncix = 0                  ! Index to current cix block
      nicix = 0                 ! Number of inequivalent subblocks, to be determined in this loop
      do  curr = 1, s_dmft%ncatom  ! Loop over all sites with correlated orbitals
        if (procid == master) then
          read(ifi,*) ib,k,locrot ! Site, number of l's in correlated block.   See Local above for locrot
        endif
        call mpibc1(ib,1,2,.false.,'','')
        call mpibc1(k,1,2,.false.,'','')
        call mpibc1(locrot,1,2,.false.,'','')

        loro = mod(locrot, 3)
        shift = locrot/3
C       Bug (see Bugs, above) : l and qspl should be identified with inequivalent subblock, not with site
        do  j = 1, k            ! Loop over all l-blocks for this atom
          ncix = ncix+1
          dmft_ib(ncix) = ib    ! Site for this correlated block
          if (procid == master) then
            read(ifi,*) l, qspl, icix ! l for this block; for qapl see Local above; icix : index to inequivalent block
          endif
          call mpibc1(l,1,2,.false.,'','')
          call mpibc1(qspl,1,2,.false.,'','')
          call mpibc1(icix,1,2,.false.,'','')

          dmft_icix(ncix) = icix
          if (icix == 1) lcix = .true. ! At least one correlated block points to first inequivalent subblock
          icix = iabs(icix)
          if (iinear(ncix,icix,dmft_icix,1) < ncix) cycle  ! cix is equivalent to a prior one
          if (dmft_l(icix) == 0) dmft_l(icix) = l
          if (dmft_l(icix) /= l) call rx('readindmfl : inconsistent l for this block')
          if (dmft_qsplit(icix) == 0) dmft_qsplit(icix) = qspl
          if (dmft_qsplit(icix) /= qspl) call rx('readindmfl : inconsistent qsplit for this block')
          if (dmft_ndim(icix) == 0) dmft_ndim(icix) = 2*l+1
          if (dmft_ndim(icix) /= 2*l+1) call rx('readindmfl : inconsistent ndim for this block')
          maxdim = max(maxdim,2*l+1)  ! Largest dimension of a correlated block
          nicix = max(nicix,icix)     ! Number of inequivalent correlated blocks
        enddo                   ! loop over l channels within a site

C   ... Read info for local rotations of Cartesian axes (we are not using them yet)
        if (procid == master) then
        if (loro > 0) then
          read(ifi,*)
          if (loro == 2) read(ifi,*)
        elseif (loro<0) then
          read(ifi,*)
          read(ifi,*)
          read(ifi,*)
        endif
        if (shift/=0) read(ifi,*)
        endif  ! procid == master
      enddo ! loop over atoms with correlated blocks

      if (.not. lcix) call rx('readindmfl (abort) : first inequivalent correlated block never specified')
C     Copy accumulated information into s_dmft
      s_dmft%ncix = ncix        ! Not clear if it is ever used
      s_dmft%nicix = nicix
      allocate(s_dmft%ib(ncix),s_dmft%icix(ncix))
      call icopy(ncix,dmft_ib,1,s_dmft%ib,1)
      call icopy(ncix,dmft_icix,1,s_dmft%icix,1)
      allocate(s_dmft%l(nicix),s_dmft%qsplit(nicix))
      allocate(s_dmft%ndim(nicix))    ! size of sigind matrix for each cix block
      call icopy(nicix,dmft_l,1,s_dmft%l,1)
      call icopy(nicix,dmft_qsplit,1,s_dmft%qsplit,1)
      call icopy(nicix,dmft_ndim,1,s_dmft%ndim,1)
      deallocate(dmft_l,dmft_qsplit,dmft_ndim,dmft_ib,dmft_icix)
      call info5(10,0,0,' %i total (%i inequivalent) correlated blocks, among %i sites',
     .  ncix,nicix,s_dmft%ncatom,0,0)

C --- Read additional information related to each inequivalent subblock ---
C     Read number of inequivalent subblocks in file (each site associated w/ 1 inequivalent subblock)
C     maxdim has already been determined; maxsize = largest number of channels, but never used.
      if (procid == master) then
        read(ifi,'(a1)') ch     ! comment line signalling start of data for correlated subblocks
      endif
      call mpibcc(ch, 1,.false.,'','')
      if (ch /= '#') call rx('readindmfl: expected header to inequivalent blocks section to begin with "#"')
C     indmfl file contains : nicix, maxdim, maxsize.  nicix and maxdim are already known, maxsize will be calculated internally
C     Sanity checks
      if (procid == master) then
        read(ifi,*) i,j         ! , maxdim, maxsize  not read.  nicix = number of independent cix blocks
      endif
      call mpibc1(i,1,2,.false.,'','')
      call mpibc1(j,1,2,.false.,'','')
      if (i /= s_dmft%nicix) call
     .  info2(20,0,0,' (warning): require %i inequivalent blocks but file specifies %i',nicix,i)
      nicix = i                 ! nicix should be i, anyway
      if (j /= maxdim) call rx2('readindmfl (abort): expected %i for subblock but file specifies %i',maxdim,j)

C --- Read data for each inequivalent correlated subblock ---
      allocate(iprm(maxdim*maxdim*nicix*nsp)) ! Work array to make permutation table
      allocate(s_dmft%sigind(maxdim,maxdim,nicix,nsp)) ! correlated index
      allocate(s_dmft%cf(maxdim,maxdim,nicix)) ! transformation matrix
      allocate(ndsigi(0:nicix,2))       ! Number of correlated orbitals in icix block
      allocate(wk(2*maxdim))            ! work array
      allocate(cixmap(nicix,2))         ! Mapping of a cix block to an equivalent block
      s_dmft%cf = 0
      s_dmft%sigind = 0
      call iinit(ndsigi,size(ndsigi))
      ndsigi(0,1) = 0                ! 0th subblock has no correlated matrix elements
      ndsig = 0                      ! Following loop makes ndsig=highest index to ME
      nzsig = 0                      ! Following loop makes nzsig=total number of nonzero ME
!     qident = .true.
      lsame = .false.                ! Loop will set lsame=T if a cix block is mapped to an equivalent one
      do  icix = 1, nicix
        cixmap(icix,1) = icix
        if (procid == master) then
          read(ifi,*) i, ndim   ! i=icix, ndim=rank-of-l-block, size=number of nonzero components : not read but inferred from sigind
        endif
        call mpibc1(i,1,2,.false.,'','')
        call mpibc1(ndim,1,2,.false.,'','')
        if (i /= icix) call rx2('readindmfl (abort): subblocks must be sequential : subblock %i labelled as %i',icix,i)
        if (s_dmft%ndim(icix) /= ndim) call rx2('readindmfl dimension mismatch: unexpected ndim=%i for block %i',ndim,icix)
        if (procid == master) then
        read(ifi,*)             ! Skip comment "# ... Independent components ..."
        read(ifi,*)             ! Skip information legend, e.g. "'z^2' 'x^2-y^2' 'xz' 'yz' 'xy'"
        read(ifi,*)             ! Skip comment "# ... Sigind follows ..."
        prev = 999 ; curr = 0   ! Indices to smallest, largest indices in sigind for this cix, made in next loop
        do  i = 1, ndim         !Read (i,j) array sigind
          read(ifi,*) (s_dmft%sigind(i,j,icix,1),j=1,ndim) !nonzero matrix elements of sigma in block for this icix
          do  j = 1, ndim
            if (s_dmft%sigind(i,j,icix,1) /= 0) then
              do  cix = 1, ncix
                if (iabs(s_dmft%icix(cix)) == icix) nzsig = nzsig+1
              enddo
              prev = min(prev,s_dmft%sigind(i,j,icix,1))
            endif
            curr = max(curr,s_dmft%sigind(i,j,icix,1))
          enddo
        enddo
        endif
        call mpibc1(prev,1,2,.false.,'','')
        call mpibc1(curr,1,2,.false.,'','')
        call mpibc1(nzsig,1,2,.false.,'','')
        call mpibc1(s_dmft%sigind,maxdim*maxdim*nicix,2,.false.,'','')
        ndsigi(icix,2) = curr-prev+1

        if (prev > ndsigi(icix-1,1)+1) call rx2(
     .    'sigind not contiguous with prior subblock: 1st element (%i) differs from expected value (%i)',
     .    prev,ndsigi(icix-1,1)+1)

C  ...  Data for rotation of local axes
        if (procid == master) then
          read(ifi,*)           ! Comment line flagging that transformation matrix follows
          do  i = 1, ndim
            read(ifi,*) (wk(j),j=1,2*ndim)
            do  j = 1, ndim
              s_dmft%cf(i,j,icix) = dcmplx(wk(2*j-1),wk(2*j))
            enddo
          enddo
        endif
C       call zprm('cf',2,s_dmft%cf(1,1,icix),maxdim,ndim,ndim)

C   ... If not contiguous with prior block look for match to a prior cix
        if (prev /= ndsigi(icix-1,1)+1) then
          do  cix = 1, icix-1
            lsame = .true.
            do  i = 1, ndim
              do  j = 1, ndim
                if (.not. lsame) exit
                if (s_dmft%sigind(i,j,icix,1) /= s_dmft%sigind(i,j,cix,1)) then
                  lsame = .false.
                endif
              enddo
            enddo
            if (lsame) then
              call info2(20,0,0,' map cix=%i to cix=%i',icix,cix)
              k = cix
              cixmap(icix,1) = cix
              exit
            endif
          enddo
          if (.not. lsame) call rx1('elements in cix block %i overlap with prior block',icix)
          ndsigi(icix,1) = ndsigi(icix-1,1)
        else
          ndsigi(icix,1) = curr
C         Sanity check: sort the indices to ensure they form a contiguous set
          call ivheap(1,maxdim*maxdim,s_dmft%sigind(1,1,icix,1),iprm,1)
          call info5(40,0,0,' block %,2i:  rank: %i  nonzero matrix elements: %i',
     .      icix,s_dmft%ndim(icix),ndsigi(icix,2),0,0)
          prev = ndsigi(icix-1,1)
          do  is = 1, maxdim*maxdim
            curr = ival(s_dmft%sigind(1,1,icix,1),iprm(is))
            if (curr == 0 .or. curr == prev) cycle
            if (curr == prev+1) then
              prev = curr
            else
              call rx2('readindmfl (abort): indices in sigind block %i must be contiguous: gap at %i',icix,prev)
            endif
          enddo
        endif

C   ... qident true if all transf are identity, false otherwise
C        do  i = 1, ndim
C          if (.not. qident) exit
C          do  j = 1, ndim
C            if (i == j) then
C              ident = 1
C            else
C              ident = 0
C            endif
C            if (abs(s_dmft%cf(i,j,icix)-ident)>1e-5) then
C              qident = .false.
C              exit
C            endif
C          enddo
C        enddo

      enddo

C --- Remake arrays if any icix were mapped to another icix ---
C     if (lsame .or. .true.) then
      if (lsame) then

C     Count k = new number of icix
      k = 0; nzsig = 0
      do  icix = 1, nicix
        if (cixmap(icix,1) == icix) k = k+1 ! this icix was not mapped
      enddo

C ... Update ndsigi, nzsig, s_dmft%sigind, s_dmft%cf
      allocate(iwk(maxdim,maxdim,nicix),zwk(maxdim,maxdim,nicix))
      call icopy(size(iwk),s_dmft%sigind,1,iwk,1);
      call zcopy(size(zwk),s_dmft%cf,1,zwk,1)
      deallocate(s_dmft%sigind,s_dmft%cf)
      allocate(s_dmft%sigind(maxdim,maxdim,k,nsp),s_dmft%cf(maxdim,maxdim,k))
      k = 0
      do  icix = 1, nicix
        icix0 = cixmap(icix,1)  ! Original icix it maps to

C       Update arrays for icix not mapped
        if (icix0 == icix) then ! this icix was not mapped
          k = k+1
          cixmap(icix,2) = k    ! New icix index
          ndsigi(k,1) = ndsigi(k-1,1) + ndsigi(icix,2) ! update ndsigi(:,1)
          call icopy(maxdim**2,iwk(1,1,icix),1,s_dmft%sigind(1,1,k,1),1)
          call zcopy(maxdim**2,zwk(1,1,icix),1,s_dmft%cf(1,1,k),1)
          cycle
        endif

        cixmap(icix,2) = cixmap(icix0,2) ! New icix index

C       Checks and changes for newly equivalent blocks
        if (s_dmft%l(icix) /= s_dmft%l(icix0))
     .    call rxi('readindmfl : inconsistent l for cix block',icix)
        if (s_dmft%qsplit(icix) /= s_dmft%qsplit(icix0))
     .    call rxi('readindmfl : inconsistent qsplit for cix block',icix)
        if (s_dmft%ndim(icix) /= s_dmft%ndim(icix0))
     .    call rxi('readindmfl : inconsistent ndim for cix block',icix)

      enddo
      nicix = k                 ! Update number of inequivalent blocks
      s_dmft%nicix = k

C     Update s_dmft%icix and nzsig
      do  cix = 1, ncix
        icix = s_dmft%icix(cix)
        s_dmft%icix(cix) = cixmap(icix,2)
        nzsig = nzsig + ndsigi(icix,2)
      enddo

      endif ! lsame

C --- Entry point when indmfl file is not read ---
  100 continue

      ndsig = ndsigi(nicix,1)
      call info2(30,0,0,' %i nonzero (%i inequivalent) matrix elements',nzsig,ndsig)
C --- Remake sigind for spin polarized case ---
      if (nsp == 2) then
        call iinit(s_dmft%sigind(:,:,:,2),size(s_dmft%sigind(:,:,:,1)))
        call info0(30,1,0,'   channels for 2nd spin block%N   m1   m2   cix chn(1) chn(2)')
        do  icix = 1, nicix
          prev = ndsigi(icix-1,1)
          curr = ndsigi(icix-1,1)*2
          k = ndsigi(icix,1) - ndsigi(icix-1,1) ! size of cix block
          do  i = 1, s_dmft%ndim(icix) ! Translate sigind(i,j) for this icix to spin1, spin2
            do  j = 1, s_dmft%ndim(icix)
              if (s_dmft%sigind(i,j,icix,1) /= 0) then
                iat = s_dmft%sigind(i,j,icix,1)
                if (iprint() >= 30) then
                  write (stdo,333) i,j,icix,
C     .              prev,curr,
     .              iat+curr-prev,
     .              iat+curr-prev+k
                endif
  333           format(3i5,2i6)
                s_dmft%sigind(i,j,icix,1) = iat+curr-prev
                s_dmft%sigind(i,j,icix,2) = iat+curr-prev+k
              endif
            enddo
          enddo
        enddo                   ! Inequivalent blocks
      endif                     ! nsp == 2

C --- s_dmft%iasig for all blocks & sanity checks on matrix elements ---
      call info0(30,1,0,'  correlated channels:%N  chan equiv m1   m2   isp icix  cix  ib')
      allocate(s_dmft%nzsigi(0:ncix))     ! Number of correlated orbitals in cix block
      allocate(s_dmft%iasig(5,nzsig*nsp)) ! iasig will be filled out here
      allocate(iasig(4,nzsig*nsp))
      call ivheap(1,maxdim*maxdim*nicix*nsp,s_dmft%sigind,iprm,1)  ! sort the elements in a big heap

      ix = 0    ! running index to iasig (index to inequivalent nonzero matrix element)
      ix0 = 0   ! running index to iasig, start of current block
      ixs = 0   ! running index to s_dmft%iasig:  should reach nzsig at end
      prev = 0  ! Index to preceding ME
      icix0 = 1 ! Index to prior icix (to flag when icix changes)
      s_dmft%nzsigi(0) = 0  ! 0th subblock has no correlated matrix elements
      do  is = 1, maxdim*maxdim*nicix*nsp  ! Loop over all inequivalent blocks
        k = iprm(is)
        curr = ival(s_dmft%sigind,k)
        if (curr == 0) cycle

        ix = ix + 1  ! New inequivalent matrix element
C       Unpack row, column, icix, isp
        k = k-1
        i = mod(k,maxdim)
        j = mod(k/maxdim,maxdim)
        icix = mod(k/(maxdim**2),nicix)
        isp = mod(k/(maxdim**2*nicix),nsp)  ! isp is outer index
        i = i+1; j=j+1; icix=icix+1; isp = isp+1
        iasig(1,ix) = i ! row index relative to start of cix block
        iasig(2,ix) = j ! column index relative to start of cix block
        iasig(3,ix) = isp
        iasig(4,ix) = icix
        if (is == maxdim*maxdim*nicix*nsp) then  ! Last element: fill in last block
          icix = icix+1
          ix = ix+1
        endif
C       Fill in s_dmft%iasig of prior block if a  new block has started
        if (icix /= icix0) then
          do  cix = 1, ncix
            if (iabs(s_dmft%icix(cix)) /= icix0) cycle
            do  k = ix0+1, ix-1
              ixs = ixs+1
              s_dmft%iasig(1:4,ixs) = iasig(1:4,k)
              i    = iasig(1,k)
              j    = iasig(2,k)
              isp  = iasig(3,k)
              if (s_dmft%icix(cix) < 0 .and. nsp == 2) isp = nsp+1-isp
              isigind = s_dmft%sigind(i,j,icix0,isp) ! index to a nonzero element in DMFT sigma
              s_dmft%iasig(4,ixs) = isigind ! isigind is a column in sig.inp
              s_dmft%iasig(5,ixs) = cix
              if (iprint() >= 30) write(stdo,"(10i5)")
     .          ixs,isigind,i,j,isp,icix0,cix,s_dmft%ib(cix)
            enddo
            s_dmft%nzsigi(cix) = ixs
          enddo
          ix0 = ix-1
          icix0 = icix
          if (iprint() >= 30) write(stdo,*)
        endif
      enddo
      deallocate(iprm,iasig)
      if (ixs /= nzsig*nsp)
     .  call rx2('readindmfl : mismatch : expected %i matrix elements but found %i',nzsig*nsp,ixs)

      s_dmft%ndsig = ndsig*nsp
      s_dmft%nzsig = nzsig*nsp


!     index first found for this inq cix block. this atom is the reference for the orbital.
!     All other atoms belonging to this icix will be rotated.
!      if( s_lat%nsafm \=0 ) then
!         allocate(ipca(nbas),istab(nbas,2),g(9,2),ag(3,2))
!         call suafmsym(s_lat,nglob('nbas'),ipca,istab,g,ag)
!     endif
      allocate(nicixi(ncix))    ! Loop over inequivalent cix only
      call ineqcix(s_dmft,ncix,nicixi)


      allocate(s_dmft%ig(ncix))
      s_dmft%ig = 1             ! if operation found, apply unity.

      allocate(s_dmft%rmat(maxdim, maxdim, ncix))


      do cix = 1, ncix
         if (nicixi(cix) >= 0) then

            call findroteqv(iabs(s_dmft%icix(cix)), cix, s_lat, s_dmft, igfnd)


            if( igfnd == -1 ) then ! no symetry operation found

               write(*,*) 'For the independent cix block ',icix
               write(*,'(A, I2,A,I2)') ' No symetry group operation between atom '
     .              ,s_dmft%ib(ind_ff),' and atom',s_dmft%ib(cix)
               if (s_dmft%icix(cix) > 0) then
                  write(*,*) 'please put them in two differents impurities block'
               else
                  write(*,*) ' please put them in two differents impurities block'
                  write(*,*) ' or'
                  write(*,*) ' please introduce by hand  the AFM symetry in SYMGRP in the ctrl file'
               endif
               call rx('sudmft.f : problem with cix block')

            endif
         endif
      enddo


      end subroutine readindmfl




      subroutine findroteqv(cixref,cix,s_lat,s_dmft,igfnd)
!     find rotation from pos2 to pos1 i.e pos1= trans(pos2)
!     if not found, igfnd = -1
!     if found, return rmat matric rotation in cubic harmonics given trans
      use structures
      implicit none
      type(str_lat)  :: s_lat
      type(str_dmft) ::  s_dmft
      integer :: cixref,cix

!     local variable
      integer :: ngrp, igfnd
      real(8) :: posref(3),pos(3),plat(3,3),qlat(3,3)

      logical :: eqv
      real(8) :: d(3)
      integer :: ig, m

      integer :: l, ldcix, nlml
      real(8),allocatable ::  rmattmp(:,:)


      ngrp = max(s_lat%nsgrp, iabs(s_lat%nsafm)) ! takes AFM case if needed

      posref(:) = s_lat%pos(:, s_dmft % ib(cixref))
      pos(:) = s_lat%pos(:, s_dmft % ib(cix))


      igfnd = -1


      do ig = 1, ngrp

         do  m = 1, 3
            d(m) = s_lat%symgr(3*(m-1)+1,ig)*pos(1) + s_lat%symgr(3*(m-1)+2,ig)*pos(2)
     .           + s_lat%symgr(3*(m-1)+3,ig)*pos(3) + s_lat% ag(m,ig) - posref(m)
         enddo
         call shorbz(d, d, s_lat%plat, s_lat%qlat)
         if (d(1)*d(1)+d(2)*d(2)+d(3)*d(3) < 1d-9) then
            eqv = .true.
            igfnd = ig
            go to 10
         endif
      enddo
      if (igfnd == -1) return     ! failed to found transformation


 10   continue ! group symetry found
      l = s_dmft%l(iabs(s_dmft%icix(cix)))
      ldcix = 2*l + 1
      nlml =  l**2 + 2*l + 1


      allocate(rmattmp(nlml,nlml))

      call ylmrtg(nlml,s_lat%symgr(1,igfnd),rmattmp)
      s_dmft % ig(cix) = igfnd
      s_dmft % rmat(:,:,cix) = rmattmp((nlml - ldcix + 1 ):nlml,(nlml - ldcix+ 1):nlml)
      end subroutine findroteqv
