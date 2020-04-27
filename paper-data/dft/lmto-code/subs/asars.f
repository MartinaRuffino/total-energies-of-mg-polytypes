      subroutine asars(mode,s_ctrl,s_site,s_spec,s_lat,s_pot,s_bz,pnu,qnu,lbin,ifi)
C- ASA file I/O to restart file
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nccomp nclasp ncomp idcc ics
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:dclabl
Cio    Passed to:  asars2
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  pos class spec force vel pnu pz
Co     Stored:     pos pos0 force vel spec clabel pnu pz
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  asars1
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:rmt nr a idmod z
Cio    Passed to:  asars2
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pnu qnu rhrmx vrmax
Cio    Passed to:  asars2
Cio  s_bz   :struct for the Brillouin Zone; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  *
Ci Inputs/Outputs
Ci   mode  :if 1, attempt to set up s_pot%v0 for file I/O, if not allocated
Ci   Program flow for potential s_pot%v0, File READ and mode=1
Ci    if s_pot%v0 is not associated, allocate and initialize v(1)=NULL for each class.
Ci    v(1) will be different from NULL for each class where potential is read from disk.
Ci    if s_pot%v0 is already associated, potentials are read in for classes where found.
Ci    v(1) is left untouched for classes not read.
Ci    s_pot%v0 is read for each site For equivalent sites, the last on one read from
Ci    disk is kept and the others discarded.  No attempt to symmetrize.
Ci
Ci   Program flow for potential s_pot%v0, File WRITE and mode=1
Ci    if s_pot%v0 is not associated, allocate and initialize v(1)=NULL for each class.
Ci      Read whatever potentials are available from atom files
Ci    Otherwise, s_pot%v0 is not modified by asars
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Ci   lbin  :T file I/O in binary mode
Ci         :F file I/O in ascii mode
Cio Inputs/Outputs
Cio  pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Cio         pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Cio         class-based array; I/O as site-based
Cio  qnu   :energy-weighted moments of the sphere charges
Cio         class-based array; I/O as site-based
Cr Remarks
Cr   This is a bridging code to map class-resolved ASA parameters
Cr   into format suitable for iorsa.
Cu Updates
Cu   24 Sep 17 asars and read/write site data
Cu   07 Aug 16 Generate qnur from qnu.  Future: iorsa should read/write qnur.
Cu   06 May 15 Adapted to work with CPA -- treated as a special case
Cu   17 Jun 13 Replace f77 pointers with f90 ones
Cu   25 Oct 11 Started migration to f90 structures
Cu   07 Jul 05 version -1.04 ... adapted to altered iorsa.
Cu   28 Feb 02 version -1.03 ... adapted to altered iorsa.
Cu   10 May 01 version -1.  Only writes pnu,qnu
C  ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      logical lbin
      integer mode,ifi
      double precision pnu(*),qnu(*)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_pot)::   s_pot
      type(str_bz)::    s_bz
C ... Dynamically allocated local arrays
      integer, allocatable :: ipc(:),ipa(:)
      real(8), allocatable :: pnul(:),qnul(:),qnulr(:,:)
      real(8), pointer :: qnur(:,:)
C ... Local parameters
      character*32 jobid
      integer nbas,nl,nsp,i,nit,nspec,nat
      integer,parameter :: n0=10, NULLI=-99999
      procedure(integer) iorsa,asars2,nglob

      jobid = 'ASA output'
      nbas = nglob('nbasp')
      nsp = nglob('nsp')
      nl = nglob('nl')
      nspec = nglob('nspec')
      nit = 1

C ... Minimal initial allocation of s_pot%v0.  Resize s_pot%v0 on the fly as needed
      if (.not. associated(s_pot%v0) .and. mode == 1) then
        call ptr_pot(s_pot,8+1,'v0',nsp,s_ctrl%nclasp+s_ctrl%nccomp,i) ! minimal allocation
        s_pot%v0(1,:) = NULLI
        if (ifi<0) i = asars2(4,s_ctrl,s_spec,s_pot,abs(ifi)) ! read s_pot%v0 from atom files
      endif

C     Handle CPA differently
      if (s_ctrl%nccomp /= 0) then
        nat = nbas
        i = iorsa(3,s_ctrl,s_site,s_spec,s_lat,s_pot,s_bz,jobid,nbas,nat,nspec,qnul,nit,lbin,ifi)
        if (i < 0) call rx('failed to read rsta file')
        i = asars2(1,s_ctrl,s_spec,s_pot,ifi)
        if (i < 0) call rx('failed to read rsta file')
        return
      endif

C ... Get class list
      allocate(ipc(nbas),ipa(nbas))
      ipc(:) = s_site(1:nbas)%class

C ... Memory for qnu, site resolved
      allocate(qnul(3*n0*nsp*nbas)); call dpzero(qnul,3*n0*nsp*nbas)
      allocate(pnul(n0*nsp*nbas)); call dpzero(pnul,n0*nsp*nbas)
      allocate(qnulr(32*n0**2,nbas)); call dpzero(qnulr,32*n0**2*nbas)
      if (mod(s_ctrl%lrel,10) == 2) then
        qnulr(1,:) = NULLI+1  ! Not available; try to read
      else
        qnulr(1,:) = NULLI    ! Make no attempt to read
      endif

C ... Poke class-based P,Q into site->pnu and qnul
      if (ifi < 0) call asars1(10,s_site,nbas,nsp,nl,pnul,qnul,qnulr,ipc,ipa,pnu,qnu,qnulr)

C ... File I/O.
C     for now
      nat = nbas
C     In future, iorsa should read qnur
      i = iorsa(0,s_ctrl,s_site,s_spec,s_lat,s_pot,s_bz,jobid,nbas,nat,nspec,qnul,nit,lbin,ifi)

C ... Poke site struc and qnul into P,Q
      if (associated(s_pot%qnur)) then
        qnur => s_pot%qnur
      else
        allocate(qnur(32*n0**2,nbas))
        call dpzero(qnur,32*n0**2*nbas)
        qnur(1,:) = NULLI
      endif
      if (ifi > 0) call asars1(11,s_site,nbas,nsp,nl,pnul,qnul,qnulr,ipc,ipa,pnu,qnu,qnur)

      deallocate(ipc,ipa,qnul,pnul)
      if (.not. associated(s_pot%qnur)) deallocate(qnur)
      end
      subroutine asars1(mode,s_site,nbas,nsp,nl,pnu,qnu,qnur,ipc,ipa,pl,ql,qlr)
C- Kernel that copies ASA parameters into form readable by iorsa
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit :
Ci         : 0  Copies classed-based P,Q to site-based P,Q
Ci         : 1  Symmetrizes site-based P,Q and copy to class-based P,Q
Ci         :10s digit specifies how site P is packed/unpacked
Ci         : 1  Store or retrieve site P to/from site->pnu
Ci              Note: In storage mode, site->pnu is also copied to
Ci              array pnu as part of the symmetrization procedure.
Ci              Symmetrized pnu is repacked into site->pnu
Ci         : 2  Store or retrieve site P to/from array pnu
Ci         : 3  In class-to-site copy: combination of 1+2
Ci         :    In site-to-class copy, same as 2, but also
Ci         :    symmetrized pnu is packed into site->pnu.
Ci   mode  :0 copy ASA pl,ql to site-based pnu,qnu and also site->pnu
Ci         :1 symmetrize and reverse copy
Ci  s_site :struct for site-specific data; see structures.h
Ci    Elts read: pnu (mode=1)
Ci    Stored:    pnu (mode=0)
Ci    Passed to: *
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nl    :(global maximum l) + 1
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   ipa   :integer work array (holds list of class members)
Cio Inputs/Outputs
Cio   pl   :sphere boundary conditions, by class
Cio        :Input if mode=0; output if mode=1
Cio   ql   :energy-weighted moments of the sphere charges, by class
Cio        :Input if mode=0; output if mode=1
Cio  pnu   :same as pl, but resolved by site (also different dim.)
Cio  qnu   :same as ql, but resolved by site (also different dim.)
Cio        :Output if mode=0; input if mode=1
Cb Bugs
Cb   ssite shouldn't be needed here ... should move copy to pnu out.
Cu Updates
Cu    1 Mar 02  Reworked mode options.
Cu   10 May 01  First attempt.  Writes P,Q only
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,mode,nl,nsp,ipc(nbas),ipa(7)
      integer,parameter :: n0=10, NULLI=-99999
      double precision pl(nl,nsp,*),ql(3,nl,nsp,*),qlr(4*nl*2*nl*2*2,*)
      double precision pnu(n0,nsp,nbas),qnu(3,n0,nsp,nbas),qnur(32*n0**2,nbas)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
C ... Local parameters
      integer ib,ic,nrclas,mode0,mode1
C     character*8 clabl
      double precision xx,wk(32*n0**2),ploc(n0,2)
C     logical aiopot

      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)

      call sanrg(.true.,mode0,0,1,'asars1:','1s digit mode')
      call sanrg(.true.,mode1,1,3,'asars1:','10s digit mode')

C --- Create symmetrized pnu ---
      if (mode0 == 1) then

C   ... Extract site->pnu into pnu array
        if (mode1 == 1) then
          do  ib = 1, nbas
            ploc = s_site(ib)%pnu
            call dcopy(n0*nsp,ploc,1,pnu(1,1,ib),1)
          enddo
        endif

C   ... Symmetrize pnu,qnu
        ic = 0
   10   continue
        ic = ic+1
        call psymr0(-2,-ic,nbas,ipc,xx,xx,ipa,nrclas)
        if (nrclas > 0) then
          call psymq0(nrclas,nsp,ipa,wk,3*n0,qnu)
          call psymq0(nrclas,1,ipa,wk,32*n0**2,qnur)
          call psymq0(nrclas,nsp,ipa,wk,n0,pnu)
          goto 10
        endif

C   ... Copy pnu to site->pnu
        if (mode1 == 1 .or. mode1 == 3) then
          call dpzero(ploc,n0*2)
          do  ib = 1, nbas
            call dcopy(n0*nsp,pnu(1,1,ib),1,ploc,1)
            s_site(ib)%pnu = ploc
          enddo
        endif
      endif

C --- Poke P,Q (class) to/from P,Q (site) ---
      do  ib = 1, nbas
        ic = ipc(ib)

C   ... Pack pl into site->pnu and/or array pnu; and ql to qnu
        if (mode0 == 0) then
          call dpzero(ploc,n0*2)
          call dcopy(nl,pl(1,1,ic),1,ploc,1)
          call dcopy(nl,pl(1,nsp,ic),1,ploc(1,nsp),1)
          if (mod(mode1,2) == 1) s_site(ib)%pnu = ploc
          if (mode1 >= 2) then
            call dpzero(pnu(1,1,ib),n0*nsp)
            call dcopy(1*nl,pl(1,1,ic),  1,pnu(1,1,ib),  1)
            call dcopy(1*nl,pl(1,nsp,ic),1,pnu(1,nsp,ib),1)
          endif
          call dpzero(qnu(1,1,1,ib),3*n0*nsp)
          call dcopy(3*nl,ql(1,1,1,ic),  1,qnu(1,1,1,ib),  1)
          call dcopy(3*nl,ql(1,1,nsp,ic),1,qnu(1,1,nsp,ib),1)

C   ... Copy pnu to pl and qnu to ql
        else
C          pnu(:,:,ib) = s_site(ib)%pnu
          call dcopy(nl,pnu(1,1,ib),1,pl(1,1,ic),1)
          call dcopy(nl,pnu(1,nsp,ib),1,pl(1,nsp,ic),1)
          call dcopy(3*nl,qnu(1,1,1,ib),  1,ql(1,1,1,ic),  1)
          call dcopy(3*nl,qnu(1,1,nsp,ib),1,ql(1,1,nsp,ic),1)
          if (qnur(1,ib) == NULLI) cycle  ! No attempt to read
          if (qnur(1,ib) == NULLI+1) then
            call qnu2qnur(1,nl,nsp,ql(1,1,1,ic),ql,qnur(1,ib))
          endif
          call dcopy(4*nl*2*nl*2*2,qnur(1,ib),1,qlr(1,ic),1)
        endif

      enddo

      return

C ... Error exit
C  999 continue
C      call rx('iorsa, class'//clabl//': missing',line)

      end

      integer function asars2(mode,s_ctrl,s_spec,s_pot,ifi)
C- I/O of sphere data from rsta file, for CPA case
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nclasp ncomp idcc ics
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:dclabl
Cio    Passed to:  *
Ci   s_ctrl%ics   :species table
Ci         :class ic belongs to species ics(ic), ic=1:nclasp
Ci         :extra classes are tacked to the end of ics
Ci         :ics(i>ic) index to parent species.  nclasp<i<=nclasp+nccomp
Co   s_ctrl%idcc  :points to links between CPA and parent classes
Co         :s_ctrl%idcc(ic=1:nclasp) points to child classes:
Co         :If class ic is NOT a CPA class, s_ctrl%idcc(ic) points to itself
Co         :If ic IS a CPA class, s_ctrl%idcc(ic) points to the first of
Co         :of the CPA or DLM classes associated with ic.
Co         :s_ctrl%idcc(i=nclasp+1:nclasp+nccomp) (CPA classes):
Co         :index to the parent class for this CPA class
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:rmt nr a idmod z
Cio    Passed to:  *
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:pnu qnu rhrmx vrmax
Cio    Passed to:  *
Ci Inputs
Ci   mode  :1  I/O pnu,qnu
Ci         :2  (read only) I/O maximum value of nr for all species, for dimensioning
Ci         :   ifi not used in this mode
Ci         :4  I/O potential for all species.  s_pot%v0 must be allocated
Ci   ifi   :file logical unit, but >0 for read, <0 for write
Co Outputs
Co   sphere data read from aiomom
Co   asars2: returns -1 if error; otherwise maximum nr (mode 2), otherwise 0
Cs Command-line switches
Cr Remarks
Cr   Normally this data is read as part of iorsa.f.
Cr   It was split off because no particular species can be associated with any site.
Cu Updates
Cu   06 May 15  First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer mode,nsp,ifi
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_pot)::   s_pot
      type(str_spec)::  s_spec(*)
C ... Dynamically allocated local arrays
      real(8), pointer :: v0(:)
C ... Local parameters
      integer,parameter :: nrmx=5001
      logical lpass
      integer procid,mpipid,nl
      integer, parameter :: master=0
      integer stdo,ic,ipr,jfi,i,k,jc,is,lmxa,l,iv(10)
      integer nr2,nsp2,nrx
      double precision a2,rmax2
      character clabl*8
C     double precision xx,wk(3*n0*2),ploc(n0,2)
      real(8) :: bxc(3)
      real(8),parameter :: NULLR=-99999
      procedure(integer) :: nglob,iprint,fopna
      procedure(logical) :: aiomom,aiopot

      stdo = nglob('stdo')
      procid = mpipid(1)
      ipr = iprint()
      nsp = nglob('nsp')
      nl  = nglob('nl')
      nrx = 0
      allocate(v0(nrmx*2))
C     nclaspd = s_ctrl%nclasp + s_ctrl%nccomp

C --- Input ---
      if (ifi > 0) then
        jfi = ifi
        if (procid == master) then

        lpass = .true.
        do  ic = 1, s_ctrl%nclasp
          k = 0
          if (s_ctrl%nccomp /= 0) then
            k = s_ctrl%ncomp(ic)
          endif
          do  i = 0, k
            if (k == 0) then
              jc = ic
            else
              if (i == 0) cycle
              jc = s_ctrl%idcc(ic)+i-1
            endif
            call r8tos8(s_ctrl%dclabl(jc),clabl)
            is = s_ctrl%ics(ic)
            lmxa = s_spec(is)%lmxa

            if (ipr >= 50) then
              write(stdo,349) jc,clabl,lmxa,s_spec(is)%rmt,s_spec(is)%nr,s_spec(is)%a,(s_pot%pnu(l+1,jc),l=0,lmxa)
            endif
            bxc = 0
            if (mod(mode,2) == 1) then
              lpass = aiomom(clabl,s_pot%pnu(1,jc),s_pot%qnu(1,jc),[NULLR],s_spec(is)%idmod,nl,lmxa,
     .          nsp,s_spec(is)%z,s_pot%rhrmx(jc),s_pot%vrmax(1,jc),bxc,ifi)
              if (.not. lpass) goto 10
            endif
            if (mod(mode/2,2) == 1 .or. mod(mode/4,2) == 1) then
              jfi = fopna(clabl,-1,0)
              nr2=0; nsp2=0; a2=0; rmax2=0; v0(1) = 0
              lpass = aiopot(nr2,nsp2,a2,rmax2,[NULLR],v0,jfi)
              nrx = max(nrx,nr2)
              call fclr(clabl,jfi)
              if (.not. lpass) goto 10
              if (mod(mode/4,2) == 1) then
                if (.not. associated(s_pot%v0)) call rx('asars: attempt to raad pot into unassociated pointer')
                iv(1:2) = shape(s_pot%v0)
                if (iv(1) < nr2*nsp) call resizev0(s_pot,nr2*nsp)
                jfi = fopna(clabl,-1,0)
                lpass = aiopot(nr2,nsp2,a2,rmax2,bxc,s_pot%v0(1,jc),jfi)
                call fclr(clabl,jfi)
              endif
            endif
          enddo
        enddo
   10   continue
        endif ! procid == master
        call mpibc1(lpass,1,1,.false.,'iorsa','read error')
        call mpibc1(nrx,1,2,.false.,'iorsa','maximum nr')
        call mpibc1(s_pot%pnu,size(s_pot%pnu),4,.false.,'iorsa','pnu')
        call mpibc1(s_pot%qnu,size(s_pot%qnu),4,.false.,'iorsa','qnu')
        asars2 = -1; if (lpass) asars2 = 0
        if (lpass .and. mod(mode/2,2)+mod(mode/4,2) /= 0) asars2 = nrx
        return

C --- Output ---
      else
        asars2 = 0
        if (procid /= master) then
          return
        endif
        jfi = -ifi

        if (ipr >= 50) write(stdo,8)
    8   format(/3x,'ic:spc     la    rmt     nr   a     pnu')
        write(jfi,18) 'class densities'
   18   format('----------------------- ',a,' -----------------------')
        do  ic = 1, s_ctrl%nclasp
          k = 0
          if (s_ctrl%nccomp /= 0) then
            k = s_ctrl%ncomp(ic)
          endif
          do  i = 0, k
            if (k == 0) then
              jc = ic
            else
              if (i == 0) cycle
              jc = s_ctrl%idcc(ic)+i-1
            endif
            call r8tos8(s_ctrl%dclabl(jc),clabl)
            is = s_ctrl%ics(ic)
            lmxa = s_spec(is)%lmxa

            if (ipr >= 50) then
              write(stdo,349) jc,clabl,lmxa,s_spec(is)%rmt,s_spec(is)%nr,s_spec(is)%a,(s_pot%pnu(l+1,jc),l=0,lmxa)
  349         format(i5,':',a8,i2,f9.5,i5,f6.3,1x,8f6.3)
            endif

            lpass = aiomom(clabl,s_pot%pnu(1,jc),s_pot%qnu(1,jc),[NULLR],s_spec(is)%idmod,nl,lmxa,
     .        nsp,s_spec(is)%z,s_pot%rhrmx(jc),s_pot%vrmax(1,jc),(/0d0,0d0,1d0/),ifi)
          enddo
        enddo
      endif

      end

      subroutine resizev0(s_pot,nrsp)
C- Resize potential array s_pot%v0
C ----------------------------------------------------------------------
Ci Inputs
Cio Structures
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:v0
Cio    Passed to:  *
Ci   nrsp
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cu Updates
Cu   22 Sep 17 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nrsp
C ... Dynamically allocated arrays
      real(8),pointer:: wk(:,:)
C ... Local parameters
      integer i12(2),i
      double precision xx
C ... For structures
!       include 'structures.h'
      type(str_pot) ::  s_pot

      if (.not. associated(s_pot%v0)) then
        call rx('s_pot%%v0 not allocated ... cannot resize')
      endif
      i12 = shape(s_pot%v0)
      allocate(wk(i12(1),i12(2)))
      call dcopy(i12(1)*i12(2),s_pot%v0,1,wk,1)

      call ptr_pot(s_pot,8+1,'v0',nrsp,i12(2),xx)

      do  i = 1, i12(2)
        call dcopy(i12(1),wk(1,i),1,s_pot%v0(1,i),1)
      enddo
      deallocate(wk)
      end

      subroutine rdasapq(sopts,dc,s_ctrl,s_site,s_spec,s_pot,ib1,ib2,offc,errmsg)
C- Read parameters from ASA restart file
C ----------------------------------------------------------------------
Ci Inputs
Ci   sopts :string containing read options, separated by global terminator
Ci         :Options are:
Ci         :fn=filename
Ci         :rdim[=a,nr,z,lma,class,ves,rmt]
Ci   dc    :global terminator.  Marks end of sopts, and possibility beginning.
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  iorsa
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec class
Co     Stored:     class
Co     Allocated:  *
Cio    Elts passed:qnu pnu
Cio    Passed to:  iorsa
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  *
Co     Stored:     a nr z lmxa rmt
Co     Allocated:  *
Cio    Elts passed:z
Cio    Passed to:  iorsa
Cio  s_pot  :struct for information about the potential; see structures.h
Ci     Elts read:  *
Co     Stored:     ves
Co     Allocated:  ves
Cio    Elts passed:ves
Cio    Passed to:  iorsa
Ci   ib1   :Read data for sites ib1:ib2
Ci   ib2   :ditto
Ci   offc  :class offset.  Classes are mapped to site here.  class(ib) = ib+offc
Co Outputs
Co   errmsg:If an error appears during parsing, errmsg returns nonblank
Co         :containing a short message indicating the nature of the error
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   24 Feb 18 (MvS) First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      character*256 sopts,errmsg,dc*1
      integer ib1,ib2,offc
C ... For structures
!     include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_pot)::   s_pot

      type(str_bz)::    s_bz0
      type(str_pot)::   s_pot0
      type(str_bz)::    s_lat0
      type(str_spec),pointer:: s_spec0(:)
      type(str_site),pointer:: s_site0(:)

C ... Dynamically allocated arrays
      real(8),allocatable:: qnus(:,:,:,:)

C ... Local parameters
      real(8),parameter  :: NULLR = -99999
      integer, parameter :: n0=10
      integer i,j1,j2,ifi,nit,nbas0,nat0,nspec0,nsp0,is0,is,ib,nsp,ldims,ic0,ic
      character*256 dc2*1,strn,fn
      procedure(integer) :: fopnx,fopng,fextg,wordsw,iorsa,isanrg,nglob

      ldims = 0
      nsp = nglob('nsp')

      dc2 = sopts(1:1)
      fn = 'rsta'  ! Default file name
      i = fextg(fn(5:))
C ... Supplied file name
      if (wordsw(sopts,dc2,'fn','=',j1) /= 0) then
        j1 = j1+1
        call nwordg(sopts,0,dc//dc2//' ',1,j1,j2)
        fn = sopts(j1:j2)
      endif
C ... Whether to import dimensioning parameters a,nr,rmt
      i = wordsw(sopts,dc2,'rdim',dc2//'= ',j1)
      if (i > 0) then
        ldims = 15+16
        if (sopts(j1:j1) == '=') then
          ldims = 0
          i = scan(sopts(j1:),dc//dc2//' ') + j1-2
          if (index(sopts(j1:i),'a') > 0) ldims = ldims+1
          if (index(sopts(j1:i),'nr') > 0) ldims = ldims+2
          if (index(sopts(j1:i),'z') > 0) ldims = ldims+4
          if (index(sopts(j1:i),'lmxa') > 0) ldims = ldims+8
          if (index(sopts(j1:i),'class') > 0) ldims = ldims+16
          if (index(sopts(j1:i),'ves') > 0) ldims = ldims+32
          if (index(sopts(j1:i),'rmt') > 0) ldims = ldims+64
        endif
      endif

      ifi = fopnx(fn,72,-1,-1)
      if (ifi /= 1) then
        errmsg = "missing file '"//trim(fn)//"'"
        return
      endif
      ifi = fopng(fn,-1,0)
      errmsg = "file mismatch"
      nsp0 = iorsa(-1,s_ctrl,s_site,s_spec,s_lat0,s_pot,s_bz0,
     .  strn,nbas0,nat0,nspec0,[0d0],nit,.false.,ifi)
      if (nsp0 < 0) return
      errmsg = "spin mismatch"
      if (nsp /= nsp0) return
      if (ib2-ib1+1 /= nbas0) then
        j1 = isanrg(nbas0,ib2-ib1+1,ib2-ib1+1,' rdasapq:','nbas',.false.)
        return
      endif
      allocate(qnus(3,n0,nsp,nbas0))
      call dvset(qnus,1,size(qnus),NULLR)
      allocate(s_site0(nbas0),s_spec0(nspec0))
      call ptr_site(s_site0,0,' ',-nbas0,0,0,[0])
C     call ptr_spec(s_spec0,0,' ',nspec0,0,0,[0])
      do  i = 1, nspec0
        s_spec0(i)%lmxl = NULLR; s_spec0(i)%lfoca = 0; s_spec0(i)%rfoca = 0
      enddo
      do  i = 1, nbas0
        s_site0(i)%class = i; s_site0(i)%pz = 0
        s_site(i)%class = i+offc
      enddo
      allocate(s_pot0%ves(nbas0))
      nullify(s_pot0%v0)

      j1 = iorsa(-4,s_ctrl,s_site0,s_spec0,s_lat0,s_pot0,s_bz0,
     .  strn,nbas0,nat0,nspec0,qnus,nit,.false.,ifi)

      if (.not. associated(s_pot%ves)) call ptr_pot(s_pot,4+1,'ves',1,0,[NULLR])

      do  i = 1, nbas0
        ib = i+ib1-1
        is0 = s_site0(i)%spec
        is = s_site(ib)%spec
        ic0 = s_site0(ib)%class
        ic  = s_site(ib)%class

C       Poke dimensioning scalar quantities into species structure
        if (mod(ldims/1,2) > 0) s_spec(is)%a   = s_spec0(is0)%a
        if (mod(ldims/2,2) > 0) s_spec(is)%nr  = s_spec0(is0)%nr
        if (mod(ldims/4,2) > 0) s_spec(is)%z   = s_spec0(is0)%z
        if (mod(ldims/8,2) > 0) s_spec(is)%lmxa = s_spec0(is0)%lmxa
        if (mod(ldims/64,2) >0) s_spec(is)%rmt = s_spec0(is0)%rmt
        if (size(s_pot%ves) < ic) then
          call ptr_pot(s_pot,2,'ves',ic,0,[0d0])
          s_pot%ves(ic) = s_pot0%ves(ic0)
        endif
        if (mod(ldims/32,2) > 0 .or. s_pot%ves(ic) == NULLR) s_pot%ves(ic) = s_pot0%ves(ic0)

C       Check that atomic numbers synchronize
        if (s_spec0(is0)%z /= s_spec(is)%z) then
          call awrit2(' rdasapq: atomic number mismatch mapping file site %i -> superlattice site %i',
     .      errmsg,len(errmsg),0,i,i+ib1-1)
          return
        endif

C       Poke P,Q and some scalar quantities into site structure
        call dpzero(s_site(ib)%qnu,size(s_site(ib)%qnu))
        call dcopy(size(s_site(ib)%pnu),s_site0(i)%pnu,1,s_site(ib)%pnu,1)
        call dcopy(size(s_site(ib)%pz),s_site0(i)%pz,1,s_site(ib)%pz,1)
        call dcopy(3*n0*nsp,qnus(1,1,1,i),1,s_site(ib)%qnu,1)
        s_site(i)%force = s_site0(i)%force
        s_site(i)%vel = s_site0(i)%vel

      enddo

      deallocate(s_site0,s_spec0,qnus,s_pot0%ves)
      errmsg = ' '

      end
