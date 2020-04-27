C   This file contains routines to manage structures.
C   C syntax is followed to the extent possible in the f90 language.
C   It avoids use of modules and interfaces, except for the input routines.
C   Structures are allocated in this file, and passed as arguments.
C
C   All the pointer handling is packaged in this file, which contains
C   routines for retrieving site or species data in vector format, and
C   pointer management for structure elements with dynamical arrays
C
C   Some workarounds are necessary.
C   1. f90 destroys allocated arrays and structures on exit from a routine.
C   Workaround: dynamical arrays are allocated by generic public pointers (pointers.f)
C   2. The size of a structure cannot be determined.
C   Workaround: a C call is made to determine size
C   3. MPI broadcast does not work for pointer elements.
C   Workaround: pointers are saved before, and restored after, MPI call
C
C     Notes on adding a new element to a structure.
C     For definiteness, add qbyl to s_pot.
C     1. In bcast_strx, add lines (use the same ordering as structures77.h)
C        s_potl(1)%qbyl => s_pot%qbyl
C        s_pot%qbyl => s_potl(1)%qbyl
C     2. In routine ptr_pot in the Inputs/Outputs section
C        * Add a comment indicating the cast, rank, and index for entry point.
C          Maintain alphabetic order of the list.
C          Suppose for definiteness qbyl will appear as entry 30.
C          In the same comment section, increment already existing subsequent entries
C          (in this case entries >=30 need be incremented)
C        * Increase the dimension of variable eltlst by one and insert qbyl
C          as the 30rd string in eltlst.
C        * Append a number in the data statement for eltidx
C        * Add s_pot%qbyl in the nullify command immediately following eltlst
C        * Insert a case statement (case 30) and tailor the statement.
C          Look for an object with same cast and rank and change name (e.g. qnu -> qbyl)
C        * Increment the case numbers of all the cases following your newly inserted one.
      subroutine sitepack(s_site,ib1,ib2,elt,n,iarr,darr)
C- Unpacks or packs element of site structure
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read: spec class norb offh pl plv dlmcl ncomp
Ci                rmaxs vshft mpole relax pos pos0 force vel
Ci                eula bfield pnu pz rsm cg0
Co     Stored:    spec class norb offh pl plv dlmcl ncomp
Co                rmaxs vshft mpole relax pos pos0 force vel
Co                eula bfield pnu pz rsm cg0
Cio    Passed to: *
Ci Inputs
Ci   ib1,ib2 :extract or pack elements in range (ib1,ib2)
Ci   elt     :string denoting one member of the site structure, e.g. 'pos'
Ci           :* elt corresponds to the symbol in s_site; see structures.h.
Ci           :This element will be either be read from s_site or stored in it.
Ci           :* Data is packed into (or from) either iarr or darr, depending
Ci           :on whether the element has integer or double cast.
Ci           :* The first character of elt is used as a switch to indicate
Ci           :whether data is to copied from s_site or stored into it.
Ci           :If the first character is '-', e.g. '-pos', data is stored in s_site
Ci           :and the proper element name is elt(2:).  Otherwise data
Ci           :is read from s_site.
Ci   n     :should be one for scalar entries in s_site
Ci         :number of elements to copy for vector entries.
Cio Inputs/Outputs
Cio   One of iarr or darr is extracted from the following element of s_site,
Cio   (or copied into it if first char of elt is '-'), depending on cast of
Cio   of s_site member corresponding to elt.
Cio   For efficiency elt is translated internally into an integer, isw.
Cio     Element   Copied from/to     isw
Cio       'spec'      iarr            0
Cio       'class'     iarr            1
Cio       'norb'      iarr            2
Cio       'offh'      iarr            3
Cio       'pl'        iarr            4
Cio       'plv'       iarr            5
Cio       'dlmcl'     iarr            6
Cio       'ncomp'     iarr            7
Cio       'oomg'      iarr            8 --- SUPERSEDED
Cio       'oomgn'     iarr            9 --- SUPERSEDED
Cio
Cio       'rmaxs'     darr           10
Cio       'vshft'     darr           11
Cio       'mpole'     darr           12
Cio
Cio       'ogc'       iarr           13 --- SUPERSEDED
Cio       'ocpawt'    iarr           14 --- SUPERSEDED
Cio
Cio       'relax'     iarr           15
Cio
Cio       'pos'       darr           20
Cio       'pos0'      darr           21
Cio       'force'     darr           22
Cio       'vel'       darr           23
Cio       'eula'      darr           24
Cio       'bfield'    darr           25
Cio       'pnu'       darr           26
Cio       'pz'        darr           27
Cio       'rsm'       darr           28
Cio       'cg0'       darr           29
Cl Local variables
Ci   isw   :an integer with a 1-1 correspondence with a subset of
Ci         :members of s_site, given by eltidx. See Bugs, below.
Cb Bugs
Ci   As written the code does not access every member of s_site.
Ci   To add a member, do the following:
Ci   1. increment parameter nelt and append a new string to eltlst
Ci   2. associate an integer isw with this string (see table above)
Ci      and append it to eltidx
Ci   3. Add a new line to read element in section READ from s_site
Ci   4. Add a new line to write element in section WRITE to s_site
Cr Remarks
Cr   Extracts or packs elements of s_site
Cu Updates
Cu   28 Jun 12 Redesigned
Cu   28 Nov 11 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      character elt*(*)
      integer ib1,ib2,n,iarr(n,ib2)
      double precision darr(n,ib2)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
C ... Local parameters
      integer i,k,isw,nelt,strnsiz
      parameter (nelt=26,strnsiz=6)
      integer eltidx(nelt)
      character*(nelt*strnsiz) eltlst
      logical lread
      save eltlst,eltidx
      data eltidx /0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,21,22,23,24,
     .  25,26,27,28,29/
      data eltlst/' '/

      if (eltlst(1:1) == ' ') then
        eltlst =
     .  'spec  '//
     .  'class '//
     .  'norb  '//
     .  'offh  '//
     .  'pl    '//
     .  'plv   '//
     .  'dlmcl '//
     .  'ncomp '//
     .  'oomg  '//
     .  'oomgn '//
     .  'rmaxs '//
     .  'vshft '//
     .  'mpole '//
     .  'ogc   '//
     .  'ocpawt'//
     .  'relax '//
     .  'pos   '//
     .  'pos0  '//
     .  'force '//
     .  'vel   '//
     .  'eula  '//
     .  'bfield'//
     .  'pnu   '//
     .  'pz    '//
     .  'rsm   '//
     .  'cg0   '
      endif
      if (elt(1:1) == '-') then
        k = index(eltlst,elt(2:))
        lread = .false.
      else
        k = index(eltlst,elt)
        lread = .true.
      endif
      if (k == 0 .or. mod(k-1,6) /= 0)
     .  call rxs('sitepack:  unknown element ',elt)
      i = (k-1)/6 + 1
      isw = eltidx(i)

      do  i = ib1, ib2
        k = i
C   --- READ from s_site ---
        if (lread) then
          select case (isw)
            case (0); iarr(1,i) = s_site(k)%spec
            case (1); iarr(1,i) = s_site(k)%class
            case (2); iarr(1,i) = s_site(k)%norb
            case (3); iarr(1,i) = s_site(k)%offh
            case (4); iarr(1,i) = s_site(k)%pl
            case (5); iarr(1,i) = s_site(k)%plv
            case (6); iarr(1,i) = s_site(k)%dlmcl
            case (7); iarr(1,i) = s_site(k)%ncomp
            case (8); call rx('replace oomg with omg')
            case (9); call rx('replace oomgn with omgn')

            case (10); darr(1,i) = s_site(k)%rmaxs
            case (11); darr(1,i) = s_site(k)%vshft
            case (12); darr(1,i) = s_site(k)%mpole

            case (13); call rx('replace ogc with gc')
            case (14); call rx('replace ocpawt with cpawt')

            case (15); iarr(1:n,i) = s_site(k)%relax

            case (20); darr(1:n,i) = s_site(k)%pos
            case (21); darr(1:n,i) = s_site(k)%pos0
            case (22); darr(1:n,i) = s_site(k)%force
            case (23); darr(1:n,i) = s_site(k)%vel
            case (24); darr(1:n,i) = s_site(k)%eula
            case (25); darr(1:n,i) = s_site(k)%bfield
            case (26); darr(1:n,i) = s_site(k)%pnu(1:n,1)
            case (27); darr(1:n,i) = s_site(k)%pz(1:n,1)
            case (28); darr(1:n,i) = s_site(k)%rsm
            case (29); darr(1:n,i) = s_site(k)%cg0
          end select

C   --- Write to s_site ---
        else
          select case (isw)
            case (0); s_site(k)%spec = iarr(1,i)
            case (1); s_site(k)%class = iarr(1,i)
            case (2); s_site(k)%norb = iarr(1,i)
            case (3); s_site(k)%offh = iarr(1,i)
            case (4); s_site(k)%pl = iarr(1,i)
            case (5); s_site(k)%plv = iarr(1,i)
            case (6); s_site(k)%dlmcl = iarr(1,i)
            case (7); s_site(k)%ncomp = iarr(1,i)
            case (8); call rx('replace oomg with omg')
            case (9); call rx('replace oomgn with omgn')

            case (10); s_site(k)%rmaxs = darr(1,i)
            case (11); s_site(k)%vshft = darr(1,i)
            case (12); s_site(k)%mpole = darr(1,i)

            case (13); call rx('replace ogc with gc')
            case (14); call rx('replace ocpawt with cpawt')

            case (15); s_site(k)%relax = iarr(1:n,i)

            case (20); s_site(k)%pos = darr(1:n,i)
            case (21); s_site(k)%pos0 = darr(1:n,i)
            case (22); s_site(k)%force = darr(1:n,i)
            case (23); s_site(k)%vel = darr(1:n,i)
            case (24); s_site(k)%eula = darr(1:n,i)
            case (25); s_site(k)%bfield = darr(1:n,i)
            case (26); s_site(k)%pnu(1:n,1) = darr(1:n,i)
            case (27); s_site(k)%pz(1:n,1) = darr(1:n,i)
            case (28); s_site(k)%rsm = darr(1:n,i)
            case (29); s_site(k)%cg0 = darr(1:n,i)
          end select

        endif

      enddo

      end

      subroutine spec2class(s_spec,nclass,ips,elt,n,iarr,darr)
C- Unpacks (packs) element of species structure, ordered by class or site.
C ----------------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read: lmxa lmxb lmxf lmxl nr kmxt mxcst grp2 a rmt rsma z
Ci                eref rg qc mass stni idmod idxdn alpha hcr orbp ehvl
Co     Stored:    lmxa lmxb lmxf lmxl nr kmxt mxcst grp2 a rmt rsma z
Co                eref rg qc mass stni idmod idxdn alpha hcr orbp ehvl
Ci Inputs
Ci   nclass:number of species, classes, or sites; see ips.
Ci   ips   :maps from species to classes or sites.
Ci         :ips(1) = -1 -> no mapping: species i (un)packed to iarr or darr(:,i)
Ci         :ips(i) = class index.  Species i (un)packed to class ips(i)
Ci         :ips(i) = site  index.  Species i (un)packed to site  ips(i)
Ci   elt   :string denoting one member of the spec structure, e.g. 'lmxa'
Ci         :* elt corresponds to the symbol in s_spec; see structures.h.
Ci         :This element will be either be read from s_spec or stored in it.
Ci         :* Data is packed into (or from) either iarr or darr, depending
Ci         :on whether the element has integer or double cast.
Ci         :* The first character of elt is used as a switch to indicate
Ci         :whether data is to copied from s_spec or stored into it.
Ci         :If the first character is '-', e.g. '-pos', data is stored in s_spec
Ci         :and the proper element name is elt(2:).  Otherwise data
Ci         :is read from s_spec.
Ci   n     :should be one for scalar entries in s_spec
Ci         :number of elements to copy for vector entries.
Cio Inputs/Outputs
Cio   One of iarr or darr is extracted from the following element of s_spec,
Cio   (or copied into it if first char of elt is '-'), depending on cast of
Cio   of s_spec member corresponding to elt.
Cio   For efficiency elt is translated internally into an integer, isw.
Cio     Element   Copied from/to    isw
Cio       'lmxa'      iarr            0
Cio       'lmxb'      iarr            1
Cio       'lmxf'      iarr            2
Cio       'lmxl'      iarr            3
Cio       'nr  '      iarr            4
Cio       'kmxt'      iarr            5
Cio       'mxcst'     iarr            6
Cio       'grp2'      iarr            7
Cio
Cio       'a   '      darr           10
Cio       'rmt '      darr           11
Cio       'rsma'      darr           12
Cio       'z   '      darr           13
Cio       'eref'      darr           14
Cio       'rg  '      darr           15
Cio       'qc  '      darr           16
Cio       'mass'      darr           17
Cio       'stni'      darr           18
Cio
Cio       'idmod'     iarr           21
Cio       'idxdn'     iarr           22
Cio
Cio       'alpha'     darr           25
Cio       'hcr '      darr           26
Cio       'orbp'      darr           27
Cio       'ehvl'      darr           28
Cio       'shfac'     darr           29
Cl Local variables
Ci  k     :species
Ci  isw   :an integer with a 1-1 correspondence with a subset of
Ci        :members of s_spec, given by eltidx. See Bugs, below.
Cb Bugs
Ci   As written the code does not access every member of s_spec.
Ci   To add a member, do the following:
Ci   1. increment parameter nelt and append a new string to eltlst
Ci   2. associate an integer isw with this string (see table above)
Ci      and append it to eltidx
Ci   3. Add a new line to read element in section READ from s_spec
Ci   4. Add a new line to write element in section WRITE to s_spec
Cr Remarks
Cr   Extracts elements of s_spec into class-based or site-based vectors
Cu Updates
Cu   20 Jun 13 Added shfac
Cu   28 Nov 11 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer n
      integer nclass,ips(nclass),iarr(n,nclass)
      double precision darr(n,nclass)
      character elt*(*)
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
C ... Local parameters
      integer i,k,isw,nelt,strnsiz
      parameter (nelt=24,strnsiz=6)
      integer eltidx(nelt)
      character*(nelt*strnsiz) eltlst
      logical lprm,lread
      save eltlst,eltidx
      data eltidx /0,1,2,3,4,5,6,7,10,11,12,13,14,15,16,17,18,21,22,25,
     .  26,27,28,29/
      data eltlst/' '/

      if (eltlst(1:1) == ' ') then
        eltlst =
     .    'lmxa  '//
     .    'lmxb  '//
     '    'lmxf  '//
     '    'lmxl  '//
     '    'nr    '//
     '    'kmxt  '//
     '    'mxcst '//
     '    'grp2  '//
     '    'a     '//
     '    'rmt   '//
     '    'rsma  '//
     '    'z     '//
     '    'eref  '//
     '    'rg    '//
     '    'qc    '//
     '    'mass  '//
     '    'stni  '//
     '    'idmod '//
     '    'idxdn '//
     '    'alpha '//
     '    'hcr   '//
     '    'orbp  '//
     '    'ehvl  '//
     '    'shfac '
      endif
      if (elt(1:1) == '-') then
        k = index(eltlst,elt(2:))
        lread = .false.
      else
        k = index(eltlst,elt)
        lread = .true.
      endif
      if (k == 0 .or. mod(k-1,6) /= 0)
     .  call rxs('spec2class:  unknown element ',elt)
      i = (k-1)/6 + 1
      isw = eltidx(i)

      lprm = ips(1) > 0
      do  i = 1, nclass
        if (lprm) then
          k = ips(i)
        else
          k = i
        endif

C   --- READ from s_spec ---
        if (lread) then
          select case (isw)
            case (0 ); iarr(1,i) = s_spec(k)%lmxa
            case (1 ); iarr(1,i) = s_spec(k)%lmxb
            case (2 ); iarr(1,i) = s_spec(k)%lmxf
            case (3 ); iarr(1,i) = s_spec(k)%lmxl
            case (4 ); iarr(1,i) = s_spec(k)%nr
            case (5 ); iarr(1,i) = s_spec(k)%kmxt
            case (6 ); iarr(1,i) = s_spec(k)%mxcst
            case (7 ); iarr(1,i) = s_spec(k)%grp2

            case (10); darr(1,i) = s_spec(k)%a
            case (11); darr(1,i) = s_spec(k)%rmt
            case (12); darr(1,i) = s_spec(k)%rsma
            case (13); darr(1,i) = s_spec(k)%z
            case (14); darr(1,i) = s_spec(k)%eref
            case (15); darr(1,i) = s_spec(k)%rg
            case (16); darr(1,i) = s_spec(k)%qc
            case (17); darr(1,i) = s_spec(k)%mass
            case (18); darr(1,i) = s_spec(k)%stni

            case (21); iarr(1:n,i) = s_spec(k)%idmod(1:n)
            case (22); iarr(1:n,i) =s_spec(k)%idxdn(1:n,1)

            case (25); darr(1:n,i) = s_spec(k)%alpha(1:n)
            case (26); darr(1:n,i) = s_spec(k)%hcr(1:n)
            case (27); darr(1:n,i)=s_spec(k)%orbp(1:n,1,1)
            case (28); darr(1:n,i) = s_spec(k)%ehvl(1:n)
            case (29); darr(1:n,i) = s_spec(k)%shfac(1:n)
          end select

C   --- Write to s_spec ---
        else
          select case (isw)
            case (0 ); s_spec(k)%lmxa = iarr(1,i)
            case (1 ); s_spec(k)%lmxb = iarr(1,i)
            case (2 ); s_spec(k)%lmxf = iarr(1,i)
            case (3 ); s_spec(k)%lmxl = iarr(1,i)
            case (4 ); s_spec(k)%nr = iarr(1,i)
            case (5 ); s_spec(k)%kmxt = iarr(1,i)
            case (6 ); s_spec(k)%mxcst = iarr(1,i)
            case (7 ); s_spec(k)%grp2 = iarr(1,i)

            case (10); s_spec(k)%a = darr(1,i)
            case (11); s_spec(k)%rmt = darr(1,i)
            case (12); s_spec(k)%rsma = darr(1,i)
            case (13); s_spec(k)%z = darr(1,i)
            case (14); s_spec(k)%eref = darr(1,i)
            case (15); s_spec(k)%rg = darr(1,i)
            case (16); s_spec(k)%qc = darr(1,i)
            case (17); s_spec(k)%mass = darr(1,i)
            case (18); s_spec(k)%stni = darr(1,i)

            case (21); s_spec(k)%idmod(1:n) = iarr(1:n,i)
            case (22); s_spec(k)%idxdn(1:n,1)=iarr(1:n,i)

            case (25); s_spec(k)%alpha(1:n) = darr(1:n,i)
            case (26); s_spec(k)%hcr(1:n) = darr(1:n,i)
            case (27); s_spec(k)%orbp(1:n,1,1)=darr(1:n,i)
            case (28); s_spec(k)%ehvl(1:n) = darr(1:n,i)
            case (29); s_spec(k)%shfac(1:n) = darr(1:n,i)
          end select
        endif

      enddo
      end

      subroutine siteperm(n,iprm,s_site)
C- Permutes entries in s_site structure
C ----------------------------------------------------------------------
Ci Inputs
Ci   n     :number of entries in s_site.
Ci         :Any element for which iprm=0 is excluded.
Ci         :In this case, n is returned as the number of nonzero elements
Ci   iprm  :permutation array
Ci   s_site:site structure
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr   Inefficient, but complies w/ restrictions in f90 pointers.
Cu Updates
Cu   12 Nov 17 Any element for which iprm is excluded, which shorten the list
Cu   28 Nov 11 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer n,iprm(n)
!      include 'structures.h'
      type(str_site)::  s_site(n)
      integer i,k,j
      type(str_site),pointer:: s_site1(:)

      allocate(s_site1(n))
      j = 0
      do  i = 1, n
        k = iprm(i)
        if (k == 0) cycle
        j = j+1
        s_site1(j) = s_site(k)
      enddo
      do  i = 1, j
        s_site(i) = s_site1(i)
      enddo
      deallocate(s_site1)

      if (n /= j) then
        n = j
      endif

      end

      subroutine bcast_strx(opt,s_bz,s_ctrl,s_ham,s_pot,s_lat,
     .  s_mix,s_spec,s_site,s_str,nspec,nsite)
C- MPI broadcast a structure
C ----------------------------------------------------------------------
Ci Inputs
Ci   opt   :0 do nothing
Ci         :If abs(opt) contains this bit, do:
Ci         :2**0 broadcast s_spec
Ci         :2**1 broadcast s_site
Ci         :2**2 broadcast s_bz
Ci         :2**3 broadcast s_ctrl
Ci         :2**4 broadcast s_ham
Ci         :2**5 broadcast s_pot
Ci         :2**6 broadcast s_lat
Ci         :2**7 broadcast s_mix
Ci         :2**8 broadcast s_str
Ci         :numbers can be taken in any combination
Ci         :sign(opt) is used as a flag: if opt<0
Ci         :do not broadcast structure, but initialize to 0
Ci   nspec
Ci   nsite
Co Outputs
Cr   One or more structures is MPI broadcast
Cl Local variables
Cl         :
Cr Remarks
Cr   Structure length is determined by C routine ssizei
Cr   Warning! init branch is dangerous
Cu Updates
Cu   16 Dec 11 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer opt,nspec,nsite
C ... For structures
!      include 'structures77.h'
      type(str_bz)::   s_bz
      type(str_ctrl):: s_ctrl
      type(str_lat)::  s_lat
      type(str_ham)::  s_ham
      type(str_pot)::  s_pot
      type(str_mix)::  s_mix
      type(str_site):: s_site(*)
      type(str_spec):: s_spec(*)
      type(str_str)::  s_str
C     type(str_move):: s_move
C     type(str_tb):: s_tb
C     type(str_optic):: s_optic
C     type(str_gw):: s_gw

C ... Local parameters
      logical linit
      integer i0,ssizei,i,optl
      type(str_bz)::    s_bzl(2)
      type(str_ctrl)::  s_ctrll(2)
      type(str_ham)::   s_haml(2)
      type(str_lat)::   s_latl(2)
      type(str_pot)::   s_potl(2)
      type(str_mix)::   s_mixl(2)
      type(str_site)::  s_sitel(2)
      type(str_spec)::  s_specl(2)
      type(str_str)::   s_strl(2)

C     for debugging
      integer procid,mpipid,nproc
      procid = mpipid(1)
      nproc =  mpipid(0)
C     nproc =  1
C     print *, 'bcast_strx procid opt', procid,opt

      optl = iabs(opt)
      linit = opt < 0

      if (IAND(optl,2**0) /= 0) then
        i0 = ssizei(s_specl(1),s_specl(2))
        do  i = 1, nspec
          if (linit) then
            call iinit(s_spec(i),i0)
          elseif (nproc > 1) then
            s_specl(1)%rhoc => s_spec(i)%rhoc ! retain pointers since not preserved
            call mpibc1(s_spec(i),i0,2,.false.,'bcast_strx','s_spec')
            s_spec(i)%rhoc => s_specl(1)%rhoc ! restore pointers
          endif
        enddo
        if (linit) call ptr_spec(s_spec,0,' ',nspec,0,0,[0])
      endif

      if (IAND(optl,2**1) /= 0) then
        i0 = ssizei(s_sitel(1),s_sitel(2))
        do  i = 1, nsite
          if (linit) then
            call iinit(s_site(i),i0)
          elseif (nproc > 1) then
            s_sitel(1)%bxc => s_site(i)%bxc  ! retain pointers since not preserved
            s_sitel(1)%cpawt => s_site(i)%cpawt
            s_sitel(1)%omg => s_site(i)%omg
            s_sitel(1)%omgn => s_site(i)%omgn
            s_sitel(1)%domg => s_site(i)%domg
            s_sitel(1)%dmat => s_site(i)%dmat
            s_sitel(1)%gc => s_site(i)%gc
            s_sitel(1)%gcu => s_site(i)%gcu
            s_sitel(1)%gcorr => s_site(i)%gcorr
            s_sitel(1)%gii => s_site(i)%gii
            s_sitel(1)%sfvrtx => s_site(i)%sfvrtx
            s_sitel(1)%vtxrel => s_site(i)%vtxrel
            s_sitel(1)%j0 => s_site(i)%j0
            s_sitel(1)%tau => s_site(i)%tau
            s_sitel(1)%pdos => s_site(i)%pdos
            s_sitel(1)%rho1 => s_site(i)%rho1
            s_sitel(1)%rho2 => s_site(i)%rho2
            s_sitel(1)%rhoc => s_site(i)%rhoc
            s_sitel(1)%rho1x => s_site(i)%rho1x
            s_sitel(1)%rho2x => s_site(i)%rho2x
            s_sitel(1)%rhocx => s_site(i)%rhocx
            s_sitel(1)%qhhl => s_site(i)%qhhl
            s_sitel(1)%qhkl => s_site(i)%qhkl
            s_sitel(1)%qkkl => s_site(i)%qkkl
            s_sitel(1)%eqhhl => s_site(i)%eqhhl
            s_sitel(1)%eqhkl => s_site(i)%eqhkl
            s_sitel(1)%eqkkl => s_site(i)%eqkkl
            s_sitel(1)%sighh => s_site(i)%sighh
            s_sitel(1)%sighk => s_site(i)%sighk
            s_sitel(1)%sigkk => s_site(i)%sigkk
            s_sitel(1)%tauhh => s_site(i)%tauhh
            s_sitel(1)%tauhk => s_site(i)%tauhk
            s_sitel(1)%taukk => s_site(i)%taukk
            s_sitel(1)%pihh => s_site(i)%pihh
            s_sitel(1)%pihk => s_site(i)%pihk
            s_sitel(1)%pikk => s_site(i)%pikk
            s_sitel(1)%sohh => s_site(i)%sohh
            s_sitel(1)%sohk => s_site(i)%sohk
            s_sitel(1)%sokk => s_site(i)%sokk
            s_sitel(1)%sighhx => s_site(i)%sighhx
            s_sitel(1)%sighkx => s_site(i)%sighkx
            s_sitel(1)%sigkkx => s_site(i)%sigkkx
            s_sitel(1)%tauhhx => s_site(i)%tauhhx
            s_sitel(1)%tauhkx => s_site(i)%tauhkx
            s_sitel(1)%taukkx => s_site(i)%taukkx
            s_sitel(1)%pihhx => s_site(i)%pihhx
            s_sitel(1)%pihkx => s_site(i)%pihkx
            s_sitel(1)%pikkx => s_site(i)%pikkx
            s_sitel(1)%thet => s_site(i)%thet
            s_sitel(1)%v0 => s_site(i)%v0
            s_sitel(1)%v1 => s_site(i)%v1

            call mpibc1(s_site(i),i0,2,.false.,'bcast_strx','s_site')

            s_site(i)%bxc => s_sitel(1)%bxc  ! restore pointers
            s_site(i)%cpawt => s_sitel(1)%cpawt
            s_site(i)%omg => s_sitel(1)%omg
            s_site(i)%omgn => s_sitel(1)%omgn
            s_site(i)%domg => s_sitel(1)%domg
            s_site(i)%dmat => s_sitel(1)%dmat
            s_site(i)%gc => s_sitel(1)%gc
            s_site(i)%gcu => s_sitel(1)%gcu
            s_site(i)%gcorr => s_sitel(1)%gcorr
            s_site(i)%gii => s_sitel(1)%gii
            s_site(i)%sfvrtx => s_sitel(1)%sfvrtx
            s_site(i)%vtxrel => s_sitel(1)%vtxrel
            s_site(i)%j0 => s_sitel(1)%j0
            s_site(i)%tau => s_sitel(1)%tau
            s_site(i)%pdos => s_sitel(1)%pdos
            s_site(i)%rho1 => s_sitel(1)%rho1
            s_site(i)%rho2 => s_sitel(1)%rho2
            s_site(i)%rhoc => s_sitel(1)%rhoc
            s_site(i)%rho1x => s_sitel(1)%rho1x
            s_site(i)%rho2x => s_sitel(1)%rho2x
            s_site(i)%rhocx => s_sitel(1)%rhocx
            s_site(i)%qhhl => s_sitel(1)%qhhl
            s_site(i)%qhkl => s_sitel(1)%qhkl
            s_site(i)%qkkl => s_sitel(1)%qkkl
            s_site(i)%eqhhl => s_sitel(1)%eqhhl
            s_site(i)%eqhkl => s_sitel(1)%eqhkl
            s_site(i)%eqkkl => s_sitel(1)%eqkkl
            s_site(i)%sighh => s_sitel(1)%sighh
            s_site(i)%sighk => s_sitel(1)%sighk
            s_site(i)%sigkk => s_sitel(1)%sigkk
            s_site(i)%tauhh => s_sitel(1)%tauhh
            s_site(i)%tauhk => s_sitel(1)%tauhk
            s_site(i)%taukk => s_sitel(1)%taukk
            s_site(i)%pihh => s_sitel(1)%pihh
            s_site(i)%pihk => s_sitel(1)%pihk
            s_site(i)%pikk => s_sitel(1)%pikk
            s_site(i)%sohh => s_sitel(1)%sohh
            s_site(i)%sohk => s_sitel(1)%sohk
            s_site(i)%sokk => s_sitel(1)%sokk
            s_site(i)%sighhx => s_sitel(1)%sighhx
            s_site(i)%sighkx => s_sitel(1)%sighkx
            s_site(i)%sigkkx => s_sitel(1)%sigkkx
            s_site(i)%tauhhx => s_sitel(1)%tauhhx
            s_site(i)%tauhkx => s_sitel(1)%tauhkx
            s_site(i)%taukkx => s_sitel(1)%taukkx
            s_site(i)%pihhx => s_sitel(1)%pihhx
            s_site(i)%pihkx => s_sitel(1)%pihkx
            s_site(i)%pikkx => s_sitel(1)%pikkx
            s_site(i)%thet => s_sitel(1)%thet
            s_site(i)%v0 => s_sitel(1)%v0
            s_site(i)%v1 => s_sitel(1)%v1

          endif
C          if (procid == 1) then
C            print *, linit,'bcast_strx procid', procid
C     .        ,linit
C     .        ,associated(s_site(i)%v0)
C     .        ,s_site(i)%v0(200)
C          endif
        enddo
        if (linit) call ptr_site(s_site,0,' ',nsite,0,0,[0])
      endif

      if (IAND(optl,2**2) /= 0) then
        i0 = ssizei(s_bzl(1),s_bzl(2))
        if (linit) then
          call iinit(s_bz,i0)
          call ptr_bz(s_bz,0,' ',0,0,[0])
        elseif (nproc > 1) then
          s_bzl(1)%dos => s_bz%dos ! retain pointers since not preserved
          s_bzl(1)%idtet => s_bz%idtet
          s_bzl(1)%ipq => s_bz%ipq
          s_bzl(1)%pdos => s_bz%pdos
          s_bzl(1)%qp => s_bz%qp
          s_bzl(1)%star => s_bz%star
          s_bzl(1)%wtkp => s_bz%wtkp
          s_bzl(1)%wtkb => s_bz%wtkb
          s_bzl(1)%swtk => s_bz%swtk
          call mpibc1(s_bz,i0,2,.false.,'bcast_strx','s_bz')

          s_bz%dos => s_bzl(1)%dos ! restore pointers
          s_bz%idtet => s_bzl(1)%idtet
          s_bz%ipq => s_bzl(1)%ipq
          s_bz%pdos => s_bzl(1)%pdos
          s_bz%qp => s_bzl(1)%qp
          s_bz%star => s_bzl(1)%star
          s_bz%wtkp => s_bzl(1)%wtkp
          s_bz%wtkb => s_bzl(1)%wtkb
          s_bz%swtk => s_bzl(1)%swtk
        endif
      endif

      if (IAND(optl,2**3) /= 0) then
        i0 = ssizei(s_ctrll(1),s_ctrll(2))
        if (linit) then
          call iinit(s_ctrl,i0)
C         call ptr_ctrl(s_ctrl,0,' ',0,0,0)
        elseif (nproc > 1) then

          s_ctrll(1)%cllst  => s_ctrl%cllst  ! retain pointers since not preserved
          s_ctrll(1)%clp    => s_ctrl%clp
          s_ctrll(1)%clssl  => s_ctrl%clssl
          s_ctrll(1)%group  => s_ctrl%group
          s_ctrll(1)%initc  => s_ctrl%initc
          s_ctrll(1)%ics    => s_ctrl%ics
          s_ctrll(1)%idcc   => s_ctrl%idcc
          s_ctrll(1)%ipc    => s_ctrl%ipc
          s_ctrll(1)%ipcp   => s_ctrl%ipcp
          s_ctrll(1)%ips    => s_ctrl%ips
          s_ctrll(1)%mxcst  => s_ctrl%mxcst
          s_ctrll(1)%nrc    => s_ctrl%nrc
          s_ctrll(1)%ncomp  => s_ctrl%ncomp
          s_ctrll(1)%pgfsl  => s_ctrl%pgfsl
          s_ctrll(1)%pgfvl  => s_ctrl%pgfvl
          s_ctrll(1)%pgord  => s_ctrl%pgord
          s_ctrll(1)%pgplp  => s_ctrl%pgplp

          s_ctrll(1)%spid   => s_ctrl%spid
          s_ctrll(1)%dclabl => s_ctrl%dclabl
          s_ctrll(1)%clabl  => s_ctrl%clabl
          s_ctrll(1)%pos    => s_ctrl%pos
          s_ctrll(1)%rmax   => s_ctrl%rmax

          call mpibc1(s_ctrl,i0,2,.false.,'bcast_strx','s_ctrl')

          s_ctrl%cllst  => s_ctrll(1)%cllst  ! restore pointers
          s_ctrl%clp    => s_ctrll(1)%clp
          s_ctrl%clssl  => s_ctrll(1)%clssl
          s_ctrl%group  => s_ctrll(1)%group
          s_ctrl%initc  => s_ctrll(1)%initc
          s_ctrl%ics    => s_ctrll(1)%ics
          s_ctrl%idcc   => s_ctrll(1)%idcc
          s_ctrl%ipc    => s_ctrll(1)%ipc
          s_ctrl%ipcp   => s_ctrll(1)%ipcp
          s_ctrl%ips    => s_ctrll(1)%ips
          s_ctrl%mxcst  => s_ctrll(1)%mxcst
          s_ctrl%nrc    => s_ctrll(1)%nrc
          s_ctrl%ncomp  => s_ctrll(1)%ncomp
          s_ctrl%pgfsl  => s_ctrll(1)%pgfsl
          s_ctrl%pgfvl  => s_ctrll(1)%pgfvl
          s_ctrl%pgord  => s_ctrll(1)%pgord
          s_ctrl%pgplp  => s_ctrll(1)%pgplp

          s_ctrl%spid   => s_ctrll(1)%spid
          s_ctrl%clabl  => s_ctrll(1)%clabl
          s_ctrl%dclabl => s_ctrll(1)%dclabl
          s_ctrl%pos    => s_ctrll(1)%pos
          s_ctrl%rmax   => s_ctrll(1)%rmax

        endif
      endif

      if (IAND(optl,2**4) /= 0) then
        i0 = ssizei(s_haml(1),s_haml(2))
        if (linit) then
          call iinit(s_ham,i0)
          call ptr_ham(s_ham,0,' ',0,0,[0])
        elseif (nproc > 1) then
          s_haml(1)%bdots => s_ham%bdots ! retain pointers since not preserved
          s_haml(1)%etrms => s_ham%etrms
          s_haml(1)%evals => s_ham%evals
          s_haml(1)%eula => s_ham%eula
          s_haml(1)%hrs => s_ham%hrs
          s_haml(1)%iaxs => s_ham%iaxs
          s_haml(1)%iprmb => s_ham%iprmb
          s_haml(1)%lmxa => s_ham%lmxa
          s_haml(1)%magf => s_ham%magf
          s_haml(1)%nprs => s_ham%nprs
          s_haml(1)%offH => s_ham%offH
          s_haml(1)%qsig => s_ham%qsig

          call mpibc1(s_ham,i0,2,.false.,'bcast_strx','s_ham')

          s_ham%bdots => s_haml(1)%bdots ! restore pointers
          s_ham%etrms => s_haml(1)%etrms
          s_ham%evals => s_haml(1)%evals
          s_ham%eula => s_haml(1)%eula
          s_ham%hrs => s_haml(1)%hrs
          s_ham%iaxs => s_haml(1)%iaxs
          s_ham%iprmb => s_haml(1)%iprmb
          s_ham%lmxa => s_haml(1)%lmxa
          s_ham%magf => s_haml(1)%magf
          s_ham%nprs => s_haml(1)%nprs
          s_ham%offH => s_haml(1)%offH
          s_ham%qsig => s_haml(1)%qsig

        endif
      endif

      if (IAND(optl,2**5) /= 0) then
        i0 = ssizei(s_potl(1),s_potl(2))
        if (linit) then
          call iinit(s_pot,i0)
          call ptr_pot(s_pot,0,' ',0,0,[0])
        elseif (nproc > 1) then
          s_potl(1)%aamom => s_pot%aamom ! retain pointers since not preserved
          s_potl(1)%bxc => s_pot%bxc
          s_potl(1)%cp => s_pot%cp
          s_potl(1)%ddpf => s_pot%ddpf
          s_potl(1)%dddpf => s_pot%dddpf
          s_potl(1)%ddpfr => s_pot%ddpfr
C         s_potl(1)%del => s_pot%del
          s_potl(1)%dlmwt => s_pot%dlmwt
          s_potl(1)%dmatk => s_pot%dmatk
          s_potl(1)%dpf => s_pot%dpf
          s_potl(1)%dpfr => s_pot%dpfr
          s_potl(1)%gibbs => s_pot%gibbs
          s_potl(1)%gma => s_pot%gma
          s_potl(1)%gmar => s_pot%gmar
          s_potl(1)%grrme => s_pot%grrme
          s_potl(1)%hab => s_pot%hab
          s_potl(1)%sab => s_pot%sab
          s_potl(1)%mad => s_pot%mad
          s_potl(1)%mxy => s_pot%mxy
          s_potl(1)%palp => s_pot%palp
          s_potl(1)%papg => s_pot%papg
          s_potl(1)%pf => s_pot%pf
          s_potl(1)%pfnc => s_pot%pfnc
          s_potl(1)%pfr => s_pot%pfr
          s_potl(1)%pmpol => s_pot%pmpol
          s_potl(1)%pnu => s_pot%pnu
          s_potl(1)%pp => s_pot%pp
          s_potl(1)%ppn => s_pot%ppn
          s_potl(1)%pprel => s_pot%pprel
          s_potl(1)%pti => s_pot%pti
          s_potl(1)%qbyl => s_pot%qbyl
          s_potl(1)%qc => s_pot%qc
          s_potl(1)%qcorr => s_pot%qcorr
          s_potl(1)%qnu => s_pot%qnu
          s_potl(1)%qnur => s_pot%qnur
          s_potl(1)%qpp => s_pot%qpp
          s_potl(1)%qt => s_pot%qt
          s_potl(1)%rhat => s_pot%rhat
          s_potl(1)%rnew => s_pot%rnew
          s_potl(1)%rhos => s_pot%rhos
          s_potl(1)%rhrmx => s_pot%rhrmx
          s_potl(1)%socscl => s_pot%socscl
          s_potl(1)%sop => s_pot%sop
          s_potl(1)%shfac => s_pot%shfac
          s_potl(1)%thetcl => s_pot%thetcl
          s_potl(1)%v0 => s_pot%v0
          s_potl(1)%vdif => s_pot%vdif
          s_potl(1)%ves => s_pot%ves
          s_potl(1)%vintr => s_pot%vintr
          s_potl(1)%vrmax => s_pot%vrmax
          s_potl(1)%vshft => s_pot%vshft
          s_potl(1)%smpot => s_pot%smpot
          s_potl(1)%smvextc => s_pot%smvextc
          s_potl(1)%smrho => s_pot%smrho
          s_potl(1)%smrout => s_pot%smrout
          s_potl(1)%smcv0 => s_pot%smcv0
          s_potl(1)%GFr => s_pot%GFr

          call mpibc1(s_pot,i0,2,.false.,'bcast_strx','s_pot')

          s_pot%aamom => s_potl(1)%aamom ! restore pointers
          s_pot%bxc => s_potl(1)%bxc
          s_pot%cp => s_potl(1)%cp
          s_pot%ddpf => s_potl(1)%ddpf
          s_pot%dddpf => s_potl(1)%dddpf
          s_pot%ddpfr => s_potl(1)%ddpfr
C         s_pot%del => s_potl(1)%del
          s_pot%dlmwt => s_potl(1)%dlmwt
          s_pot%dmatk => s_potl(1)%dmatk
          s_pot%dpf => s_potl(1)%dpf
          s_pot%dpfr => s_potl(1)%dpfr
          s_pot%gibbs => s_potl(1)%gibbs
          s_pot%gma => s_potl(1)%gma
          s_pot%gmar => s_potl(1)%gmar
          s_pot%grrme => s_potl(1)%grrme
          s_pot%hab => s_potl(1)%hab
          s_pot%sab => s_potl(1)%sab
          s_pot%mad => s_potl(1)%mad
          s_pot%mxy => s_potl(1)%mxy
          s_pot%palp => s_potl(1)%palp
          s_pot%papg => s_potl(1)%papg
          s_pot%pf => s_potl(1)%pf
          s_pot%pfnc => s_potl(1)%pfnc
          s_pot%pfr => s_potl(1)%pfr
          s_pot%pmpol => s_potl(1)%pmpol
          s_pot%pnu => s_potl(1)%pnu
          s_pot%pp => s_potl(1)%pp
          s_pot%ppn => s_potl(1)%ppn
          s_pot%pprel => s_potl(1)%pprel
          s_pot%pti => s_potl(1)%pti
          s_pot%qbyl => s_potl(1)%qbyl
          s_pot%qc => s_potl(1)%qc
          s_pot%qcorr => s_potl(1)%qcorr
          s_pot%qnu => s_potl(1)%qnu
          s_pot%qnur => s_potl(1)%qnur
          s_pot%qpp => s_potl(1)%qpp
          s_pot%qt => s_potl(1)%qt
          s_pot%rhat => s_potl(1)%rhat
          s_pot%rnew => s_potl(1)%rnew
          s_pot%rhos => s_potl(1)%rhos
          s_pot%rhrmx => s_potl(1)%rhrmx
          s_pot%socscl => s_potl(1)%socscl
          s_pot%sop => s_potl(1)%sop
          s_pot%shfac => s_potl(1)%shfac
          s_pot%thetcl => s_potl(1)%thetcl
          s_pot%v0 => s_potl(1)%v0
          s_pot%vdif => s_potl(1)%vdif
          s_pot%ves => s_potl(1)%ves
          s_pot%vintr => s_potl(1)%vintr
          s_pot%vrmax => s_potl(1)%vrmax
          s_pot%vshft => s_potl(1)%vshft
          s_pot%smpot => s_potl(1)%smpot
          s_pot%smvextc => s_potl(1)%smvextc
          s_pot%smrho => s_potl(1)%smrho
          s_pot%smrout => s_potl(1)%smrout
          s_pot%smcv0 => s_potl(1)%smcv0
          s_pot%GFr => s_potl(1)%GFr

        endif
      endif

      if (IAND(optl,2**6) /= 0) then
        i0 = ssizei(s_latl(1),s_latl(2))
        if (linit) then
          call iinit(s_lat,i0)
          call ptr_lat(s_lat,0,' ',0,0,0,[0])
        elseif (nproc > 1) then
          s_latl(1)%ag => s_lat%ag ! retain pointers since not preserved
          s_latl(1)%bgv => s_lat%bgv
          s_latl(1)%cg => s_lat%cg
          s_latl(1)%cy => s_lat%cy
          s_latl(1)%dlv => s_lat%dlv
          s_latl(1)%gv => s_lat%gv
          s_latl(1)%gvq => s_lat%gvq
          s_latl(1)%indxcg => s_lat%indxcg
          s_latl(1)%ips0 => s_lat%ips0
          s_latl(1)%istab => s_lat%istab
          s_latl(1)%jcg => s_lat%jcg
          s_latl(1)%kv => s_lat%kv
          s_latl(1)%igv => s_lat%igv
          s_latl(1)%igv2 => s_lat%igv2
          s_latl(1)%kv2 => s_lat%kv2
          s_latl(1)%pos => s_lat%pos
          s_latl(1)%qlv => s_lat%qlv
          s_latl(1)%symgr => s_lat%symgr
          s_latl(1)%s_sym => s_lat%s_sym

          call mpibc1(s_lat,i0,2,.false.,'bcast_strx','s_lat')

          s_lat%ag => s_latl(1)%ag ! restore pointers
          s_lat%bgv => s_latl(1)%bgv
          s_lat%cg => s_latl(1)%cg
          s_lat%cy => s_latl(1)%cy
          s_lat%dlv => s_latl(1)%dlv
          s_lat%gv => s_latl(1)%gv
          s_lat%gvq => s_latl(1)%gvq
          s_lat%indxcg => s_latl(1)%indxcg
          s_lat%ips0 => s_latl(1)%ips0
          s_lat%istab => s_latl(1)%istab
          s_lat%jcg => s_latl(1)%jcg
          s_lat%kv => s_latl(1)%kv
          s_lat%igv => s_latl(1)%igv
          s_lat%igv2 => s_latl(1)%igv2
          s_lat%kv2 => s_latl(1)%kv2
          s_lat%pos => s_latl(1)%pos
          s_lat%qlv => s_latl(1)%qlv
          s_lat%symgr => s_latl(1)%symgr
          s_lat%s_sym => s_latl(1)%s_sym

        endif
      endif

      if (IAND(optl,2**7) /= 0) then
        i0 = ssizei(s_mixl(1),s_mixl(2))
        if (linit) then
          call iinit(s_mix,i0)
C         call ptr_mix(s_mix,0,' ',0,0,0)
        elseif (nproc > 1) then
          call mpibc1(s_mix,i0,2,.false.,'bcast_strx','s_mix')
        endif
      endif

      if (IAND(optl,2**8) /= 0) then
        i0 = ssizei(s_strl(1),s_strl(2))
        if (linit) then
          call iinit(s_str,i0)
          call ptr_str(s_str,0,' ',0,0,[0])
        elseif (nproc > 1) then
          call mpibc1(s_str,i0,2,.false.,'bcast_strx','s_str')
        endif
      endif

      end

      integer function str_pack(sname,lpack,s_strn,strn)
C- Packs or unpacks one of a named list of strings
C ----------------------------------------------------------------------
Ci Inputs
Ci   lpack : 1 to pack
Ci         :-1 to unpack, copy min(size(strn),size(s_strn)) characters
Ci         :-2 to unpack, abort if size(s_strn) > size(strn)
Ci         : 0 neither packs nor unpacks, but str_pack is returned.
Ci   sname :name of string one of list below
Cio Inputs/Outputs
Cio  s_strn:structure containing all strings
Cio  strn  :local string to be packed or unpacked
Co Outputs
Co  str_pack: index corresponding to string name, if allocated
Co          :-index if not allocated (unpack only)
Cl Local variables
Cl         :
Cr Remarks
Cu Updates
Cu   31 Dec 11 First created
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer lpack
      character*(*) sname,strn
C ... For structures
!      include 'structures.h'
      type(str_strn) :: s_strn(*)
C ... Local parameters
      integer i,j,k,nstrn
      parameter (nstrn=10)
      character*(6) lsname(nstrn)
C     type(str_strn),target :: s_strnl
      character(len=1),pointer :: strnl(:)
C     common /sstrx/ s_strnl
      data lsname /
     .  'amix',
     .  'gemb',
     .  'gfopt',
     .  'jobid',
     .  'map',
     .  'mix',
     .  'mmham',
     .  'sxopt',
     .  'symg',
     .  'syml'/

C ... Find j = index to appropriate string
      do  i = 1, nstrn
        j = i
        if (lsname(i) == sname) goto 20
      enddo
      call rx('str_pack : unrecognized string : '//sname)
   20 continue
      str_pack = j

C ... No pack/unpack; just return str_pack
      if (lpack == 0) then
        return

C ... Allocate s_strn(j) and pack strn into it
      elseif (lpack == 1) then
        i = len_trim(strn)
        allocate(strnl(max(i,1)))
        strnl(1) = ' '
        do  k = 1, i
          strnl(k) = strn(k:k)
        enddo
        if (associated(s_strn(j)%strn)) deallocate(s_strn(j)%strn)
        s_strn(j)%strn => strnl
C       print *, j, s_strn(j)%strn

      else if (lpack > 0 .or. lpack < -2) then
        call rx('str_pack: improper lpack')

C ... Unpack s_strn into strn
      else
        if (associated(s_strn(j)%strn)) then
          k = size(s_strn(j)%strn)
          i = min(k,len(strn))
          if (i < k .and. lpack == -2) then
            call rxi('str_pack: increase size of '//trim(lsname(j))//
     .        ' string ... require len >',k)
          endif
          strn = ' '
          do  k = 1, i
            strn(k:k) = s_strn(j)%strn(k)
          enddo
        else
          str_pack = -str_pack
          strn = ' '
        endif
      endif

      end

      subroutine ptr_ctrl(s_ctrl,mode,elt,i1,i2,i3,arr)
C- Pointer allocation and management for s_ctrl
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl:struct contains pointers to various arrays; see structures.h
Co     Stored:    cllst clp clssl group initc icdl ics idcc ipc ipcp
Co                ips mxcst nrc ncomp pgfsl pgfvl pgord pgplp dclabl
Co                rmax
Cio    Passed to: *
Ci Inputs
Ci   mode : 0 Initialize pointers.  First call should use mode=0.
Ci        : 1s bit => allocate permanent array; associate s_ctrl%elt to it
Ci        : 2s bit => re-allocate s_ctrl%elt; reassociate s_ctrl%elt
Ci        : Either bit 1 or bit 2 must be set, but not both.
Ci        : 4s bit => copy irr or darr to internal array
Ci        : 8s bit => initialize array to 0
Ci   elt  : string identifying which pointer to operate on:
Ci         Element      type      dim>1  internal index
Ci          cllst       integer              1
Ci          clp         integer              2
Ci          clssl       integer              3
Ci          group       integer              4
Ci          initc       integer              5
Ci          icdl        integer              6  NOT USED
Ci          ics         integer              7
Ci          idcc        integer              8
Ci          ipc         integer              9
Ci          ipcp        integer              10
Ci          ips         integer              11
Ci          mxcst       integer              12
Ci          nrc         integer              13
Ci          ncomp       integer              14
Ci          pgfsl       integer              15
Ci          pgfvl       integer              16
Ci          pgord       integer              17
Ci          pgplp       integer              18
Ci          dclabl      real(8)              19
Ci          pos         real(8)      2       20
Ci          rmax        real(8)              21
Ci   i1   : leading dimension of pointer
Ci   i2   : 2nd dimension of pointer
Ci   i3   : 3rd dimension of pointer (not used so far)
Co Outputs
Co   This routine does the following depending on bits of mode
Co     - Internal pointer named by 'elt' is allocated OR
Co       the pointer is re-allocated
Co     - Data is copied from arr to it
Co     - s_ctrl%elt is associated with internal array
Cl Local variables
Cl         :
Cr Remarks
Cr   This is a workaround for f90 pointer restrictions.
Cr   * Local pointers are allocated and held in common to enable
Cr     them to be preserved.
Cr   * A corresponding pointer in s_ctrl is associated with this pointer
Cu Updates
Cu   08 May 13 Eliminate s_array
Cu   11 Oct 12 First created
C ----------------------------------------------------------------------
      use pointers
      use structures
      implicit none
C ... Passed parameters
      character elt*(*)
      integer mode,i1,i2,i3
      integer arr(*)
C ... For structures
!      include 'structures.h'
      type(str_ctrl):: s_ctrl
C ... Local parameters
      integer i,k,isw,nelt,strnsiz,mode0,mode1,mode2,mode3
      parameter (nelt=21,strnsiz=6)
      integer eltidx(nelt)
      character*(nelt*strnsiz) eltlst

C     logical lread
      save eltlst,eltidx
      data eltidx /1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     .  21/
      data eltlst/' '/

      if (eltlst(1:1) == ' ') then
        eltlst =
     .  'cllst '//
     .  'clp   '//
     .  'clssl '//
     .  'group '//
     .  'initc '//
     .  'icdl  '//
     .  'ics   '//
     .  'idcc  '//
     .  'ipc   '//
     .  'ipcp  '//
     .  'ips   '//
     .  'mxcst '//
     .  'nrc   '//
     .  'ncomp '//
     .  'pgfsl '//
     .  'pgfvl '//
     .  'pgord '//
     .  'pgplp '//
     .  'dclabl'//
     .  'pos   '//
     .  'rmax  '
        nullify(s_ctrl%cllst,s_ctrl%clp,s_ctrl%clssl,s_ctrl%group,
     .    s_ctrl%initc,s_ctrl%ics,
     .    s_ctrl%idcc,s_ctrl%ipc,s_ctrl%ipcp,s_ctrl%ips,
     .    s_ctrl%mxcst,s_ctrl%nrc,s_ctrl%ncomp,
     .    s_ctrl%pgfsl,s_ctrl%pgfvl,s_ctrl%pgord,s_ctrl%pgplp,
     .    s_ctrl%dclabl,s_ctrl%clabl,s_ctrl%pos,s_ctrl%rmax)
C        nullify(p_i1,p_i2,p_d1,p_d2,p_z1,p_z2)
      endif

      if (mode == 0) return
      mode0 = mod(mode,2)
      mode1 = mod(mode/2,2)
      mode2 = mod(mode/4,2)
      mode3 = mod(mode/8,2)

      k = index(eltlst,elt)
      if (k == 0 .or. mod(k-1,strnsiz) /= 0)
     .  call rxs('ptr_ctrl:  unknown element ',elt)

      i = (k-1)/6 + 1
      isw = eltidx(i)

C     Only ips can be reallocated for now
      if (mode0 /= 1 .and. mode1 /= 1 .and. mode0*mode1 /= 0)
     .  call rx('ptr_ctrl: improper mode')

      select case (isw)

C   ... Integer pointers
        case (1)
          if (mode1 /= 0) then
            call i1realloc(s_ctrl%cllst,size(s_ctrl%cllst),i1)
          elseif (mode0 /= 0) then
            allocate(p_i1(i1))
          endif
          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
          if (mode3 /= 0) call ivset(p_i1,1,i1,0)
          if (associated(s_ctrl%cllst)) deallocate(s_ctrl%cllst)
          s_ctrl%cllst => p_i1
        case (2)
          if (mode1 /= 0) then
            call i1realloc(s_ctrl%clp,size(s_ctrl%clp),i1)
          elseif (mode0 /= 0) then
            allocate(p_i1(i1))
          endif
          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
          if (mode3 /= 0) call ivset(p_i1,1,i1,0)
          if (associated(s_ctrl%clp)) deallocate(s_ctrl%clp)
          s_ctrl%clp => p_i1
        case (3)
          if (mode1 /= 0) then
            call i1realloc(s_ctrl%clssl,size(s_ctrl%clssl),i1)
          elseif (mode0 /= 0) then
            allocate(p_i1(i1))
          endif
          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
          if (mode3 /= 0) call ivset(p_i1,1,i1,0)
          if (associated(s_ctrl%clssl)) deallocate(s_ctrl%clssl)
          s_ctrl%clssl => p_i1
        case (4)
          if (mode1 /= 0) then
            call i1realloc(s_ctrl%group,size(s_ctrl%group),i1)
          elseif (mode0 /= 0) then
            allocate(p_i1(i1))
          endif
          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
          if (mode3 /= 0) call ivset(p_i1,1,i1,0)
          if (associated(s_ctrl%group)) deallocate(s_ctrl%group)
          s_ctrl%group => p_i1
        case (5)
          if (mode1 /= 0) then
            call i1realloc(s_ctrl%initc,size(s_ctrl%initc),i1)
          elseif (mode0 /= 0) then
            allocate(p_i1(i1))
          endif
          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
          if (mode3 /= 0) call ivset(p_i1,1,i1,0)
          if (associated(s_ctrl%initc)) deallocate(s_ctrl%initc)
          s_ctrl%initc => p_i1
        case (6)
          call rx('icdl no longer used')
C          if (mode1 /= 0) then
C            call i1realloc(s_ctrl%icdl,size(s_ctrl%icdl),i1)
C          elseif (mode0 /= 0) then
C            allocate(p_i1(i1))
C          endif
C          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
C          if (mode3 /= 0) call ivset(p_i1,1,i1,0)
C          if (associated(s_ctrl%icdl)) deallocate(s_ctrl%icdl)
C          s_ctrl%icdl => p_i1
        case (7)
          if (mode1 /= 0) then
            call i1realloc(s_ctrl%ics,size(s_ctrl%ics),i1)
          elseif (mode0 /= 0) then
            allocate(p_i1(i1))
          endif
          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
          if (mode3 /= 0) call ivset(p_i1,1,i1,0)
          if (associated(s_ctrl%ics)) deallocate(s_ctrl%ics)
          s_ctrl%ics => p_i1
        case (8)
          if (mode1 /= 0) then
            call i1realloc(s_ctrl%idcc,size(s_ctrl%idcc),i1)
          elseif (mode0 /= 0) then
            allocate(p_i1(i1))
          endif
          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
          if (mode3 /= 0) call ivset(p_i1,1,i1,0)
          if (associated(s_ctrl%idcc)) deallocate(s_ctrl%idcc)
          s_ctrl%idcc => p_i1
        case (9)
          if (mode1 /= 0) then
            call i1realloc(s_ctrl%ipc,size(s_ctrl%ipc),i1)
          elseif (mode0 /= 0) then
            allocate(p_i1(i1))
          endif
          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
          if (mode3 /= 0) call ivset(p_i1,1,i1,0)
          if (associated(s_ctrl%ipc)) deallocate(s_ctrl%ipc)
          s_ctrl%ipc => p_i1
        case (10)
          if (mode1 /= 0) then
            call i1realloc(s_ctrl%ipcp,size(s_ctrl%ipcp),i1)
          elseif (mode0 /= 0) then
            allocate(p_i1(i1))
          endif
          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
          if (mode3 /= 0) call ivset(p_i1,1,i1,0)
          if (associated(s_ctrl%ipcp)) deallocate(s_ctrl%ipcp)
          s_ctrl%ipcp => p_i1
        case (11)
          if (mode1 /= 0) then
            call i1realloc(s_ctrl%ips,size(s_ctrl%ips),i1)
          elseif (mode0 /= 0) then
            allocate(p_i1(i1))
          endif
          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
          if (mode3 /= 0) call ivset(p_i1,1,i1,0)
          if (associated(s_ctrl%ips)) deallocate(s_ctrl%ips)
          s_ctrl%ips => p_i1
        case (12)
          if (mode1 /= 0) then
            call i1realloc(s_ctrl%mxcst,size(s_ctrl%mxcst),i1)
          elseif (mode0 /= 0) then
            allocate(p_i1(i1))
          endif
          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
          if (mode3 /= 0) call ivset(p_i1,1,i1,0)
          if (associated(s_ctrl%mxcst)) deallocate(s_ctrl%mxcst)
          s_ctrl%mxcst => p_i1
        case (13)
          if (mode1 /= 0) then
            call i1realloc(s_ctrl%nrc,size(s_ctrl%nrc),i1)
          elseif (mode0 /= 0) then
            allocate(p_i1(i1))
          endif
          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
          if (mode3 /= 0) call ivset(p_i1,1,i1,0)
          if (associated(s_ctrl%nrc)) deallocate(s_ctrl%nrc)
          s_ctrl%nrc => p_i1
        case (14)
          if (mode1 /= 0) then
            call i1realloc(s_ctrl%ncomp,size(s_ctrl%ncomp),i1)
          elseif (mode0 /= 0) then
            allocate(p_i1(i1))
          endif
          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
          if (mode3 /= 0) call ivset(p_i1,1,i1,0)
          if (associated(s_ctrl%ncomp)) deallocate(s_ctrl%ncomp)
          s_ctrl%ncomp => p_i1
        case (15)
          if (mode1 /= 0) then
            call i1realloc(s_ctrl%pgfsl,size(s_ctrl%pgfsl),i1)
          elseif (mode0 /= 0) then
            allocate(p_i1(i1))
          endif
          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
          if (mode3 /= 0) call ivset(p_i1,1,i1,0)
          if (associated(s_ctrl%pgfsl)) deallocate(s_ctrl%pgfsl)
          s_ctrl%pgfsl => p_i1
        case (16)
          if (mode1 /= 0) then
            call i1realloc(s_ctrl%pgfvl,size(s_ctrl%pgfvl),i1)
          elseif (mode0 /= 0) then
            allocate(p_i1(i1))
          endif
          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
          if (mode3 /= 0) call ivset(p_i1,1,i1,0)
          if (associated(s_ctrl%pgfvl)) deallocate(s_ctrl%pgfvl)
          s_ctrl%pgfvl => p_i1
        case (17)
          if (mode1 /= 0) then
            call i1realloc(s_ctrl%pgord,size(s_ctrl%pgord),i1)
          elseif (mode0 /= 0) then
            allocate(p_i1(i1))
          endif
          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
          if (mode3 /= 0) call ivset(p_i1,1,i1,0)
          if (associated(s_ctrl%pgord)) deallocate(s_ctrl%pgord)
          s_ctrl%pgord => p_i1
        case (18)
          if (mode1 /= 0) then
            call i1realloc(s_ctrl%pgplp,size(s_ctrl%pgplp),i1)
          elseif (mode0 /= 0) then
            allocate(p_i1(i1))
          endif
          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
          if (mode3 /= 0) call ivset(p_i1,1,i1,0)
          if (associated(s_ctrl%pgplp)) deallocate(s_ctrl%pgplp)
          s_ctrl%pgplp => p_i1
C   ... Real pointers
        case (19)
          if (mode1 /= 0) then
            call d1realloc(s_ctrl%dclabl,size(s_ctrl%dclabl),i1)
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (associated(s_ctrl%dclabl)) deallocate(s_ctrl%dclabl)
          s_ctrl%dclabl => p_d1
        case (20)
            call rx('ptr_ctrl: s_ctrl%pos should be linked to'//
     .      's_lat%pos')
C          if (mode1 /= 0) then
C            call rx('ptr_ctrl: mode1 not implemented for '//elt)
C          elseif (mode0 /= 0) then
C            allocate(p_d2(i1,i2))
C          endif
C          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
C          if (associated(s_ctrl%pos)) deallocate(s_ctrl%pos)
C          s_ctrl%pos => p_d2
        case (21)
          if (mode1 /= 0) then
            call d1realloc(s_ctrl%rmax,size(s_ctrl%rmax),i1)
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (associated(s_ctrl%rmax)) deallocate(s_ctrl%rmax)
          s_ctrl%rmax => p_d1

      end select

      end subroutine ptr_ctrl

      subroutine ptr_bz(s_bz,mode,elt,i1,i2,arr)
C- Pointer allocation and management for s_bz
C ----------------------------------------------------------------------
Cio Structures
Co     Stored:    dos idtet ipq pdos qp star wtkp
Ci Inputs
Ci   mode : 0 Initialize pointers.  First call should use mode=0.
Ci        : 1s bit => allocate permanent array; associate s_bz%elt to it
Ci        : 2s bit => re-allocate s_bz%elt; reassociate s_bz%elt
Ci        : Either bit 1 or bit 2 must be set, but not both.
Ci        : 4s bit => copy irr or darr to internal array
Ci        : 8s bit => initialize array to 0
Ci   elt  : string identifying which pointer to operate on:
Ci         Element      type      dim>1  internal index
Ci          dos         real(8)              1
Ci          idtet       integer     2        2
Ci          ipq         integer              3
Ci          pdos        real(8)              4
Ci          qp          real(8)              5
Ci          star        integer              6
Ci          wtkp        real(8)              7
Ci          wtkb        real(8)              8
Ci          swtk        real(8)              9
Ci   i1   : leading dimension of pointer
Ci   i2   : 2nd dimension of pointer (not used so far)
Co Outputs
Co   This routine does the following depending on bits of mode
Co     - Internal pointer named by 'elt' is allocated OR
Co       the pointer is re-allocated
Co     - Data is copied from arr to it
Co     - s_bz%elt is associated with internal array
Cl Local variables
Cl         :
Cr Remarks
Cr   This is a workaround for f90 pointer restrictions.
Cr   * Local pointers are allocated and held in common to enable
Cr     them to be preserved.
Cr   * A corresponding pointer in s_bz is associated with this pointer
Cu Updates
Cu   11 Oct 12  First created
C ----------------------------------------------------------------------
      use pointers
      use structures
      implicit none
C ... Passed parameters
      character elt*(*)
      integer mode,i1,i2
      integer arr(*)
C ... For structures
!      include 'structures.h'
      type(str_bz)::  s_bz
C ... Local parameters
      integer i,k,isw,nelt,strnsiz,mode0,mode1,mode2,mode3
      parameter (nelt=9,strnsiz=6)
      integer eltidx(nelt)
      character*(nelt*strnsiz) eltlst
      integer i123(3)

      save eltlst,eltidx
      data eltidx /1,2,3,4,5,6,7,8,9/
      data eltlst/' '/

      if (eltlst(1:1) == ' ') then
        eltlst =
     .  'dos   '//
     .  'idtet '//
     .  'ipq   '//
     .  'pdos  '//
     .  'qp    '//
     .  'star  '//
     .  'wtkp  '//
     .  'wtkb  '//
     .  'swtk  '
      endif

      if (mode == 0) then
        nullify(s_bz%dos,s_bz%idtet,s_bz%ipq,s_bz%pdos,s_bz%qp,
     .    s_bz%star,s_bz%wtkp,s_bz%wtkb,s_bz%swtk)
        return
      endif

      mode0 = mod(mode,2)
      mode1 = mod(mode/2,2)
      mode2 = mod(mode/4,2)
      mode3 = mod(mode/8,2)

      k = index(eltlst,elt)
      if (k == 0 .or. mod(k-1,strnsiz) /= 0)
     .  call rxs('ptr_bz:  unknown element ',elt)

      i = (k-1)/6 + 1
      isw = eltidx(i)

C     No reallocation implemented
      if (mode0 /= 1 .and. mode1 /= 1 .and. mode0*mode1 /= 0)
     .  call rx('ptr_ham: improper mode')

      select case (isw)

        case (1)
          if (mode1 /= 0) then
            call d1realloc(s_bz%dos,size(s_bz%dos),i1)
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_bz%dos)) deallocate(s_bz%dos)
          s_bz%dos => p_d1

        case (2)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_bz%idtet)
            call i2realloc(s_bz%idtet,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_i2(i1,i2))
          endif
          if (mode2 /= 0) call icopy(i1*i2,arr,1,p_i2,1)
          if (mode3 /= 0) call ivset(p_i2,1,i1*i2,0)
          if (associated(s_bz%idtet)) deallocate(s_bz%idtet)
          s_bz%idtet => p_i2

C          if (mode1 /= 0) then
C            call i1realloc(s_bz%idtet,size(s_bz%idtet),i1)
C          elseif (mode0 /= 0) then
C            allocate(p_i1(i1))
C          endif
C          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
C          if (mode3 /= 0) call ivset(p_i1,1,i1,0)
C          if (associated(s_bz%idtet)) deallocate(s_bz%idtet)
C          s_bz%idtet => p_i1

        case (3)
          if (mode1 /= 0) then
            call i1realloc(s_bz%ipq,size(s_bz%ipq),i1)
          elseif (mode0 /= 0) then
            allocate(p_i1(i1))
          endif
          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
          if (mode3 /= 0) call ivset(p_i1,1,i1,0)
          if (associated(s_bz%ipq)) deallocate(s_bz%ipq)
          s_bz%ipq => p_i1

        case (4)
          if (mode1 /= 0) then
            call d1realloc(s_bz%pdos,size(s_bz%pdos),i1)
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_bz%pdos)) deallocate(s_bz%pdos)
          s_bz%pdos => p_d1

        case (5)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_bz%qp)
            call d2realloc(s_bz%qp,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (associated(s_bz%qp)) deallocate(s_bz%qp)
          s_bz%qp => p_d2

C          if (mode1 /= 0) then
C            call d1realloc(s_bz%qp,size(s_bz%qp),i1)
C          elseif (mode0 /= 0) then
C            allocate(p_d1(i1))
C          endif
C          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
C          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
C          if (associated(s_bz%qp)) deallocate(s_bz%qp)
C          s_bz%qp => p_d1

        case (6)
          if (mode1 /= 0) then
            call i1realloc(s_bz%star,size(s_bz%star),i1)
          elseif (mode0 /= 0) then
            allocate(p_i1(i1))
          endif
          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
          if (mode3 /= 0) call ivset(p_i1,1,i1,0)
          if (associated(s_bz%star)) deallocate(s_bz%star)
          s_bz%star => p_i1

        case (7)
          if (mode1 /= 0) then
            call d1realloc(s_bz%wtkp,size(s_bz%wtkp),i1)
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_bz%wtkp)) deallocate(s_bz%wtkp)
          s_bz%wtkp => p_d1

        case (8)
          if (mode1 /= 0) then
            call d1realloc(s_bz%wtkb,size(s_bz%wtkb),i1)
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_bz%wtkb)) deallocate(s_bz%wtkb)
          s_bz%wtkb => p_d1

        case (9)
          if (mode1 /= 0) then
            call d1realloc(s_bz%swtk,size(s_bz%swtk),i1)
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_bz%swtk)) deallocate(s_bz%swtk)
          s_bz%swtk => p_d1

      end select

      end subroutine ptr_bz

      subroutine ptr_lat(s_lat,mode,elt,i1,i2,i3,arr)
C- Pointer allocation and management for s_lat
C ----------------------------------------------------------------------
Cio Structures
Cio  s_lat  :struct containing lattice information; see structures.h
Co     Stored:    ag cg cy dlv gv gvq bgv indxcg ips0 istab jcg kv igv
Co                igv2 kv2 pos qlv symgr
Ci Inputs
Ci   mode : 0 Initialize pointers.  First call should use mode=0.
Ci        : 1s bit => allocate permanent array; associate s_lat%elt to it
Ci        : 2s bit => re-allocate s_lat%elt; reassociate s_lat%elt
Ci        : Either bit 1 or bit 2 must be set, but not both.
Ci        : 4s bit => copy irr or darr to internal array
Ci        : 8s bit => initialize array to 0
Ci   elt  : string identifying which pointer to operate on:
Ci         Element      type      dim>1  internal index
Ci          ag          real(8)              1
Ci          cg          real(8)              2
Ci          cy          real(8)              3
Ci          dlv         real(8)     3        4
Ci          gv          real(8)     2        5
Ci          gvq         real(8)     2        6
Ci          bgv         complex(8)           7
Ci          indxcg      integer              8
Ci          ips0        integer              9
Ci          istab       integer     2       10
Ci          jcg         integer             11
Ci          kv          integer             12
Ci          igv         integer     2       13
Ci          igv2        integer     2       14
Ci          kv2         integer     2       15
Ci          pos         real(8)     2       16
Ci          qlv         real(8)     2       17
Ci          symgr       real(8)     2       18
Ci   i1   : leading dimension of pointer
Ci   i2   : 2nd dimension of pointer
Ci   i3   : 3rd dimension of pointer (not used so far)
Co Outputs
Co   This routine does the following depending on bits of mode
Co     - Internal pointer named by 'elt' is allocated OR
Co       the pointer is re-allocated
Co     - Data is copied from arr to it
Co     - s_lat%elt is associated with internal array
Cl Local variables
Cl         :
Cr Remarks
Cr   This is a workaround for f90 pointer restrictions.
Cr   * Local pointers are allocated and held in common to enable
Cr     them to be preserved.
Cr   * A corresponding pointer in s_lat is associated with this pointer
Cu Updates
Cu   11 Oct 12  First created
C ----------------------------------------------------------------------
      use pointers
      use structures
      implicit none
C ... Passed parameters
      character elt*(*)
      integer mode,i1,i2,i3
      integer arr(*)
C ... For structures
!      include 'structures.h'
      type(str_lat)::  s_lat
C ... Local parameters
      integer i,k,isw,nelt,strnsiz,mode0,mode1,mode2,mode3
      parameter (nelt=18,strnsiz=6)
      integer eltidx(nelt)
      character*(nelt*strnsiz) eltlst

      integer i123(3)
C     real(8), pointer::    p2(:,:)

      save eltlst,eltidx
      data eltidx /1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18/
      data eltlst/' '/

      if (eltlst(1:1) == ' ') then
        eltlst =
     .  'ag    '//
     .  'cg    '//
     .  'cy    '//
     .  'dlv   '//
     .  'gv    '//
     .  'gvq   '//
     .  'bgv   '//
     .  'indxcg'//
     .  'ips0  '//
     .  'istab '//
     .  'jcg   '//
     .  'kv    '//
     .  'igv   '//
     .  'igv2  '//
     .  'kv2   '//
     .  'pos   '//
     .  'qlv   '//
     .  'symgr '
        nullify(s_lat%ag,s_lat%cg,s_lat%cy,s_lat%dlv,s_lat%gv,s_lat%gvq,
     .    s_lat%bgv,s_lat%indxcg,s_lat%ips0,s_lat%istab,s_lat%jcg,
     .    s_lat%kv,s_lat%igv,s_lat%igv2,s_lat%kv2,s_lat%pos,s_lat%qlv,
     .    s_lat%symgr)
C       nullify(p_i1,p_i2,p_d1,p_d2,p_z1,p_z2)
      endif

      if (mode == 0) return
      mode0 = mod(mode,2)
      mode1 = mod(mode/2,2)
      mode2 = mod(mode/4,2)
      mode3 = mod(mode/8,2)

      k = index(eltlst,elt)
      if (k == 0 .or. mod(k-1,strnsiz) /= 0)
     .  call rxs('ptr_lat:  unknown element ',elt)

      i = (k-1)/6 + 1
      isw = eltidx(i)

C     Only pos can be reallocated for now
      if (mode0 /= 1 .and. mode1 /= 1 .and. mode0*mode1 /= 0)
     .  call rx('ptr_lat: improper mode')

      select case (isw)

        case (1)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_lat%ag)
            call d2realloc(s_lat%ag,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (associated(s_lat%ag)) deallocate(s_lat%ag)
          s_lat%ag => p_d2
        case (2)
          if (mode1 /= 0) then
            call d1realloc(s_lat%cg,size(s_lat%cg),i1)
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_lat%cg)) deallocate(s_lat%cg)
          s_lat%cg => p_d1
        case (3)
          if (mode1 /= 0) then
            call d1realloc(s_lat%cy,size(s_lat%cy),i1)
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_lat%cy)) deallocate(s_lat%cy)
          s_lat%cy => p_d1
        case (4)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_lat%dlv)
            call d2realloc(s_lat%dlv,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (associated(s_lat%dlv)) deallocate(s_lat%dlv)
          s_lat%dlv => p_d2
        case (5)
          if (mode1 /= 0) then
            call rx('ptr_lat: mode1 not implemented for '//elt)
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (associated(s_lat%gv)) deallocate(s_lat%gv)
          s_lat%gv => p_d2
        case (6)
          if (mode1 /= 0) then
            call rx('ptr_lat: mode1 not implemented for '//elt)
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (associated(s_lat%gvq)) deallocate(s_lat%gvq)
          s_lat%gvq => p_d2
        case (7)
          if (mode1 /= 0) then
            call rx('ptr_lat: mode1 not implemented for '//elt)
          elseif (mode0 /= 0) then
            allocate(p_z1(i1))
          endif
          if (mode2 /= 0) call zcopy(i1,arr,1,p_z1,1)
          if (mode3 /= 0) call dvset(p_d2,1,2*i1,0d0)
          if (associated(s_lat%bgv)) deallocate(s_lat%bgv)
          s_lat%bgv => p_z1
        case (8)
          if (mode1 /= 0) then
            call i1realloc(s_lat%indxcg,size(s_lat%indxcg),i1)
          elseif (mode0 /= 0) then
            allocate(p_i1(i1))
          endif
          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
          if (mode3 /= 0) call ivset(p_i1,1,i1,0)
          if (associated(s_lat%indxcg)) deallocate(s_lat%indxcg)
          s_lat%indxcg => p_i1
        case (9)
          if (mode1 /= 0) then
            call i1realloc(s_lat%ips0,size(s_lat%ips0),i1)
          elseif (mode0 /= 0) then
            allocate(p_i1(i1))
          endif
          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
          if (mode3 /= 0) call ivset(p_i1,1,i1,0)
          if (associated(s_lat%ips0)) deallocate(s_lat%ips0)
          s_lat%ips0 => p_i1
        case (10)
          if (mode1 /= 0) then
            call rx('ptr_lat: mode1 not implemented for '//elt)
            i123(1:2) = shape(s_lat%istab)
C           call i2realloc(s_lat%pos,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_i2(i1,i2))
          endif
          if (mode2 /= 0) call icopy(i1*i2,arr,1,p_i2,1)
          if (mode3 /= 0) call ivset(p_i2,1,i1*i2,0)
          if (associated(s_lat%istab)) deallocate(s_lat%istab)
          s_lat%istab => p_i2
        case (11)
          if (mode1 /= 0) then
            call i1realloc(s_lat%jcg,size(s_lat%jcg),i1)
          elseif (mode0 /= 0) then
            allocate(p_i1(i1))
          endif
          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
          if (mode3 /= 0) call ivset(p_i1,1,i1,0)
          if (associated(s_lat%jcg)) deallocate(s_lat%jcg)
          s_lat%jcg => p_i1
        case (12)
          if (mode1 /= 0) then
            call i1realloc(s_lat%kv,size(s_lat%kv),i1)
          elseif (mode0 /= 0) then
            allocate(p_i2(i1,i2))
          endif
          if (mode2 /= 0) call icopy(i1*i2,arr,1,p_i2,1)
          if (mode3 /= 0) call ivset(p_i2,1,i1*i2,0)
          if (associated(s_lat%kv)) deallocate(s_lat%kv)
          s_lat%kv => p_i2
        case (13)
          if (mode1 /= 0) then
            call rx('ptr_lat: mode1 not implemented for '//elt)
          elseif (mode0 /= 0) then
            allocate(p_i2(i1,i2))
          endif
          if (mode2 /= 0) call icopy(i1*i2,arr,1,p_i2,1)
          if (mode3 /= 0) call ivset(p_i2,1,i1*i2,0)
          if (associated(s_lat%igv)) deallocate(s_lat%igv)
          s_lat%igv => p_i2
        case (14)
          if (mode1 /= 0) then
            call rx('ptr_lat: mode1 not implemented for '//elt)
          elseif (mode0 /= 0) then
            allocate(p_i2(i1,i2))
          endif
          if (mode2 /= 0) call icopy(i1*i2,arr,1,p_i2,1)
          if (mode3 /= 0) call ivset(p_i2,1,i1*i2,0)
          if (associated(s_lat%igv2)) deallocate(s_lat%igv2)
          s_lat%igv2 => p_i2
        case (15)
          if (mode1 /= 0) then
            call rx('ptr_lat: mode1 not implemented for '//elt)
          elseif (mode0 /= 0) then
            allocate(p_i2(i1,i2))
          endif
          if (mode2 /= 0) call icopy(i1*i2,arr,1,p_i2,1)
          if (mode3 /= 0) call ivset(p_i2,1,i1*i2,0)
          if (associated(s_lat%kv2)) deallocate(s_lat%kv2)
          s_lat%kv2 => p_i2
        case (16)
          if (mode1 /= 0) then
C           p_2 = reshape(p_d2,(/i1,i2/)) Doesn't work with pointers?
            i123(1:2) = shape(s_lat%pos)
            call d2realloc(s_lat%pos,i123(1),i123(2),i1,i2)
C            allocate(p2(i123(1),i123(2)))  ! Copy to hold
C            call dcopy(i123(1)*i123(2),s_lat%pos,1,p2,1)
CC           deallocate(s_lat%pos)
C            allocate(p_d2(i1,i2))
C            call dpzero(p_d2,i1*i2)
C            call dmcpy(p2,i123(1),1,p_d2,i1,1,
C     .        min(i123(1),i1),min(i123(2),i2))
C            deallocate(p2)
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (associated(s_lat%pos)) deallocate(s_lat%pos)
          s_lat%pos => p_d2
        case (17)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_lat%qlv)
            call d2realloc(s_lat%qlv,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (associated(s_lat%qlv)) deallocate(s_lat%qlv)
          s_lat%qlv => p_d2
        case (18)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_lat%symgr)
            call d2realloc(s_lat%symgr,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (associated(s_lat%symgr)) deallocate(s_lat%symgr)
          s_lat%symgr => p_d2
      end select

      end subroutine ptr_lat

      subroutine ptr_ham(s_ham,mode,elt,i1,i2,arr)
C- Pointer allocation and management for s_ham
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Co     Stored:    bdots etrms eula hrs iaxs iprmb lmxa magf nprs offH
Co                qsig
Cio    Passed to: *
Ci Inputs
Ci   mode : 0 Initialize pointers.  First call should use mode=0.
Ci        : 1s bit => allocate permanent array; associate s_ham%elt to it
Ci        : 2s bit => re-allocate s_ham%elt; reassociate s_ham%elt
Ci        : Either bit 1 or bit 2 must be set, but not both.
Ci        : 4s bit => copy irr or darr to internal array
Ci        : 8s bit => initialize array to 0
Ci   elt  : string identifying which pointer to operate on:
Ci         Element      type      dim>1  internal index
Ci          bdots       complex(8)           1
Ci          etrms       real(8)              2
Ci          eula        real(8)     2        3
Ci          hrs         real(8)              4
Ci          iaxs        integer              5
Ci          iprmb       integer              6
Ci          lmxa        integer              7
Ci          magf        real(8)              8
Ci          nprs        integer              9
Ci          offH        integer     2        10
Ci          qsig        real(8)              11
Ci   i1   : leading dimension of pointer
Ci   i2   : 2nd dimension of pointer (not used so far)
Co Outputs
Co   This routine does the following depending on bits of mode
Co     - Internal pointer named by 'elt' is allocated OR
Co       the pointer is re-allocated
Co     - Data is copied from arr to it
Co     - s_ham%elt is associated with internal array
Cl Local variables
Cl         :
Cr Remarks
Cr   This is a workaround for f90 pointer restrictions.
Cr   * Local pointers are allocated and held in common to enable
Cr     them to be preserved.
Cr   * A corresponding pointer in s_ham is associated with this pointer
Cu Updates
Cu   14 Nov 13 Additions for fully relativistic GF and noncollinear FP
Cu   11 Oct 12 First created
C ----------------------------------------------------------------------
      use pointers
      use structures
      implicit none
C ... Passed parameters
      character elt*(*)
      integer mode,i1,i2
      integer arr(*)
C ... For structures
!      include 'structures.h'
      type(str_ham)::  s_ham
C ... Local parameters
      integer i,k,isw,nelt,strnsiz,mode0,mode1,mode2,mode3
      parameter (nelt=11,strnsiz=6)
      integer eltidx(nelt)
      character*(nelt*strnsiz) eltlst

      integer i123(3)

      save eltlst,eltidx
      data eltidx /1,2,3,4,5,6,7,8,9,10,11/
      data eltlst/' '/

      if (eltlst(1:1) == ' ') then
        eltlst =
     .  'bdots '//
     .  'etrms '//
     .  'eula  '//
     .  'hrs   '//
     .  'iaxs  '//
     .  'iprmb '//
     .  'lmxa  '//
     .  'magf  '//
     .  'nprs  '//
     .  'offH  '//
     .  'qsig  '
        nullify(s_ham%bdots,s_ham%etrms,s_ham%eula,s_ham%hrs,s_ham%iaxs,
     .    s_ham%iprmb,s_ham%lmxa,s_ham%magf,s_ham%nprs,s_ham%offH,
     .    s_ham%qsig)
C       nullify(p_i1,p_i2,p_d1,p_d2,p_z1,p_z2)
      endif

      if (mode == 0) return
      mode0 = mod(mode,2)
      mode1 = mod(mode/2,2)
      mode2 = mod(mode/4,2)
      mode3 = mod(mode/8,2)

      k = index(eltlst,elt)
      if (k == 0 .or. mod(k-1,strnsiz) /= 0)
     .  call rxs('ptr_ham:  unknown element ',elt)

      i = (k-1)/6 + 1
      isw = eltidx(i)

C     No reallocation implemented
      if (mode0 /= 1 .and. mode1 /= 1 .and. mode0*mode1 /= 0)
     .  call rx('ptr_ham: improper mode')

      select case (isw)

        case (1)
          if (mode1 /= 0) then
            call rx('ptr_ham: mode1 not implemented for '//elt)
          elseif (mode0 /= 0) then
            allocate(p_z1(i1))
          endif
          if (mode2 /= 0) call zcopy(i1,arr,1,p_z1,1)
          if (mode3 /= 0) call dvset(p_z1,1,2*i1,0d0)
          if (associated(s_ham%bdots)) deallocate(s_ham%bdots)
          s_ham%bdots => p_z1
        case (2)
          if (mode1 /= 0) then
            call d1realloc(s_ham%etrms,size(s_ham%etrms),i1)
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_ham%etrms)) deallocate(s_ham%etrms)
          s_ham%etrms => p_d1
        case (3)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_ham%eula)
            call d2realloc(s_ham%eula,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (associated(s_ham%eula)) deallocate(s_ham%eula)
          s_ham%eula => p_d2
        case (4)
          if (mode1 /= 0) then
            call d1realloc(s_ham%hrs,size(s_ham%hrs),i1)
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_ham%hrs)) deallocate(s_ham%hrs)
          s_ham%hrs => p_d1
        case (5)
          if (mode1 /= 0) then
            call i1realloc(s_ham%iaxs,size(s_ham%iaxs),i1)
          elseif (mode0 /= 0) then
            allocate(p_i1(i1))
          endif
          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
          if (mode3 /= 0) call ivset(p_i1,1,i1,0)
          if (associated(s_ham%iaxs)) deallocate(s_ham%iaxs)
          s_ham%iaxs => p_i1
        case (6)
          if (mode1 /= 0) then
            call i1realloc(s_ham%iprmb,size(s_ham%iprmb),i1)
          elseif (mode0 /= 0) then
            allocate(p_i1(i1))
          endif
          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
          if (mode3 /= 0) call ivset(p_i1,1,i1,0)
          if (associated(s_ham%iprmb)) deallocate(s_ham%iprmb)
          s_ham%iprmb => p_i1
        case (7)
          if (mode1 /= 0) then
            call i1realloc(s_ham%lmxa,size(s_ham%lmxa),i1)
          elseif (mode0 /= 0) then
            allocate(p_i1(i1))
          endif
          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
          if (mode3 /= 0) call ivset(p_i1,1,i1,0)
          if (associated(s_ham%lmxa)) deallocate(s_ham%lmxa)
          s_ham%lmxa => p_i1
        case (8)
          if (mode1 /= 0) then
            call d1realloc(s_ham%magf,size(s_ham%magf),i1)
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_ham%magf)) deallocate(s_ham%magf)
          s_ham%magf => p_d1
        case (9)
          if (mode1 /= 0) then
            call i1realloc(s_ham%nprs,size(s_ham%nprs),i1)
          elseif (mode0 /= 0) then
            allocate(p_i1(i1))
          endif
          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
          if (mode3 /= 0) call ivset(p_i1,1,i1,0)
          if (associated(s_ham%nprs)) deallocate(s_ham%nprs)
          s_ham%nprs => p_i1
        case (10)
          if (mode1 /= 0) then
            call rx('not ready for re-alloc s_ham%offH')
          elseif (mode0 /= 0) then
            allocate(p_i2(i1,i2))
          endif
          if (mode2 /= 0) call icopy(i1*i2,arr,1,p_i2,1)
          if (mode3 /= 0) call ivset(p_i2,1,i1*i2,0)
          if (associated(s_ham%offH)) deallocate(s_ham%offH)
          s_ham%offH => p_i2
        case (11)

          if (mode1 /= 0) then
            i123(1:2) = shape(s_ham%qsig)
            call d2realloc(s_ham%qsig,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (associated(s_ham%qsig)) deallocate(s_ham%qsig)
          s_ham%qsig => p_d2


C          if (mode1 /= 0) then
C            call d1realloc(s_ham%qsig,size(s_ham%qsig),i1)
C          elseif (mode0 /= 0) then
C            allocate(p_d1(i1))
C          endif
C          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
C          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
C          if (associated(s_ham%qsig)) deallocate(s_ham%qsig)
C          s_ham%qsig => p_d1

      end select

      end subroutine ptr_ham

      subroutine ptr_pot(s_pot,mode,elt,i1,i2,arr)
C- Pointer allocation and management for s_pot
C ----------------------------------------------------------------------
Cio Structures
Cio  s_pot  :struct for information about the potential; see structures.h
Co     Stored:    aamom bxc cp dlmwt gibbs gma gmar grrme hab sab mad mxy
Co                palp papg pf pfnc pfr pmpol pnu pp ppn pprel pti dpf
Co                dpfr ddpf ddpfr dddpf qc qcorr qnu qbyl qpp qt rhos rhrmx
Co                socscl sop shfac thetcl v0 vdif ves vintr vrmax vshft smpot
Co                smrho GFr
Cio    Passed to: *
Ci Inputs
Ci   mode : 0 Initialize pointers.  First call should use mode=0.
Ci        : 1s bit => allocate permanent array; associate s_pot%elt to it
Ci        : 2s bit => re-allocate s_pot%elt; reassociate s_pot%elt
Ci        : Either bit 1 or bit 2 must be set, but not both.
Ci        : 4s bit => copy irr or darr to internal array
Ci        : 8s bit => initialize array to 0
Ci   elt  : string identifying which pointer to operate on:
Ci         Element      type      dim>1  internal index
Ci          aamom       real(8)              1
Ci          bxc         real(8)              2
Ci          cp          real(8)              3
Ci          del         real(8)              NOT USED
Ci          dlmwt       real(8)              4
Ci          GFr         real(8)     2        5
Ci          gibbs       real(8)              6
Ci          gma         real(8)     2        7
Ci          gmar        real(8)     2        8
Ci          grrme       real(8)              9
Ci          hab         real(8)             10
Ci          sab         real(8)             11
Ci          mad         real(8)             12
Ci          mxy         real(8)             13
Ci          palp        real(8)     2       14
Ci          papg        real(8)     2       15
Ci          pf          real(8)     2       16
Ci          pfnc        real(8)     2       17
Ci          pfr         real(8)     2       18
Ci          pmpol       real(8)             19
Ci          pnu         real(8)     2       20
Ci          pp          real(8)     2       21
Ci          ppn         real(8)             22
Ci          pprel       real(8)     2       23
Ci          pti         real(8)             24
Ci          dpf         real(8)     2       25
Ci          dpfr        real(8)     2       26
Ci          ddpf        real(8)     2       27
Ci          ddpfr       real(8)     2       28
Ci          dddpf       real(8)     2       29
Ci          qbyl        real(8)     2       30
Ci          qc          real(8)             31
Ci          qcorr       real(8)             32
Ci          qnu         real(8)     2       33
Ci          qnur        real(8)     2       34
Ci          qpp         real(8)             35
Ci          qt          real(8)             36
Ci          rhat        real(8)             37
Ci          rnew        real(8)             38
Ci          rhos        real(8)             39
Ci          rhrmx       real(8)             40
Ci          socscl      real(8)             41
Ci          sop         real(8)     2       42
Ci          shfac       real(8)     2       43
Ci          thet        real(8)     2       44
Ci          v0          real(8)     2       45
Ci          vdif        real(8)             46
Ci          ves         real(8)             47
Ci          vintr       real(8)             48
Ci          vrmax       real(8)     2       49
Ci          vshft       real(8)             50
Ci          smpot       complex(8)  2       51
Ci          smrho       complex(8)  2       52
Ci          smrout      complex(8)  2       53
Ci          smcv0       complex(8)  2       54
Ci   i1   : leading dimension of pointer
Ci   i2   : 2nd dimension of pointer (not used so far)
Co Outputs
Co   This routine does the following depending on bits of mode
Co     - Internal pointer named by 'elt' is allocated OR
Co       the pointer is re-allocated
Co     - Data is copied from arr to it
Co     - s_pot%elt is associated with internal array
Cl Local variables
Cl         :
Cr Remarks
Cr   This is a workaround for f90 pointer restrictions.
Cr   * Local pointers are allocated and held in common to enable
Cr     them to be preserved.
Cr   * A corresponding pointer in s_pot is associated with this pointer
Cu Updates
Cu   14 Nov 13 Additions for fully relativistic GF and noncollinear FP
Cu   11 Oct 12  First created
C ----------------------------------------------------------------------
      use pointers
      use structures
      implicit none
C ... Passed parameters
      character elt*(*)
      integer mode,i1,i2
      integer arr(*)
C ... For structures
!      include 'structures.h'
      type(str_pot)::  s_pot
C ... Local parameters
      integer i,k,isw,nelt,strnsiz,mode0,mode1,mode2,mode3
      integer i123(3)
      parameter (nelt=54,strnsiz=6)
      integer eltidx(nelt)
      character*(nelt*strnsiz) eltlst

      save eltlst,eltidx
      data eltidx /1,2,3,4,5,6,7,8,9,
     .  10,11,12,13,14,15,16,17,18,19,
     .  20,21,22,23,24,25,26,27,28,29,
     .  30,31,32,33,34,35,36,37,38,39,
     .  40,41,42,43,44,45,46,47,48,49,
     .  50,51,52,53,54/
      data eltlst/' '/

      if (eltlst(1:1) == ' ') then
        eltlst =
     .  'aamom '//
     .  'bxc   '//
     .  'cp    '//
     .  'dlmwt '//
     .  'gfr   '//
     .  'gibbs '//
     .  'gma   '//
     .  'gmar  '//
     .  'grrme '//
     .  'hab   '//
     .  'sab   '//
     .  'mad   '//
     .  'mxy   '//
     .  'palp  '//
     .  'papg  '//
     .  'pf    '//
     .  'pfnc  '//
     .  'pfr   '//
     .  'pmpol '//
     .  'pnu   '//
     .  'pp    '//
     .  'ppn   '//
     .  'pprel '//
     .  'pti   '//
     .  'dpf   '//
     .  'dpfr  '//
     .  'ddpf  '//
     .  'ddpfr '//
     .  'dddpf '//
     .  'qbyl  '//
     .  'qc    '//
     .  'qcorr '//
     .  'qnu   '//
     .  'qnur  '//
     .  'qpp   '//
     .  'qt    '//
     .  'rhat  '//
     .  'rnew  '//
     .  'rhos  '//
     .  'rhrmx '//
     .  'socscl'//
     .  'sop   '//
     .  'shfac '//
     .  'thet  '//
     .  'v0    '//
     .  'vdif  '//
     .  'ves   '//
     .  'vintr '//
     .  'vrmax '//
     .  'vshft '//
     .  'smpot '//
     .  'smrho '//
     .  'smrout'//
     .  'smcv0'
      endif

      if (mode == 0) then
        nullify(s_pot%aamom,s_pot%bxc,s_pot%cp,s_pot%ddpf,s_pot%dddpf,
     .    s_pot%ddpfr,s_pot%dlmwt,s_pot%gfr,s_pot%dpf,s_pot%dpfr,
     .    s_pot%gibbs,s_pot%gma,s_pot%gmar,s_pot%grrme,s_pot%hab,s_pot%sab,
     .    s_pot%mad,s_pot%mxy,s_pot%palp,
     .    s_pot%papg,s_pot%pf,s_pot%pfnc,s_pot%pfr,
     .    s_pot%pmpol,s_pot%pnu,s_pot%pp,s_pot%ppn,
     .    s_pot%pprel,s_pot%pti,s_pot%qc,s_pot%qcorr,s_pot%qnu,s_pot%qbyl,
     .    s_pot%qnur,s_pot%qpp,s_pot%qt,s_pot%rhat,s_pot%rnew,
     .    s_pot%rhos,s_pot%rhrmx,s_pot%socscl,s_pot%sop,s_pot%shfac,
     .    s_pot%thetcl,s_pot%v0,s_pot%vdif,s_pot%ves,s_pot%vintr,s_pot%vrmax,
     .    s_pot%vshft,s_pot%smpot,s_pot%smvextc,s_pot%smrho,s_pot%smrout,s_pot%smcv0)
C       nullify(p_i1,p_i2,p_d1,p_d2,p_z1,p_z2)
        return
      endif

      mode0 = mod(mode,2)
      mode1 = mod(mode/2,2)
      mode2 = mod(mode/4,2)
      mode3 = mod(mode/8,2)

      k = index(eltlst,elt)
      if (k == 0 .or. mod(k-1,strnsiz) /= 0)
     .  call rxs('ptr_pot:  unknown element ',elt)

      i = (k-1)/6 + 1
      isw = eltidx(i)

C     No reallocation implemented (?)
C      if (mode1 /= 0)
C     .  call rx('ptr_pot: mode1 not implemented for '//elt)
      if (mode0 /= 1 .and. mode1 /= 1 .and. mode0*mode1 /= 0)
     .  call rx('ptr_pot: improper mode')

      select case (isw)

        case (1)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_pot%aamom)) deallocate(s_pot%aamom)
          s_pot%aamom => p_d1

        case (2)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_pot%bxc)) deallocate(s_pot%bxc)
          s_pot%bxc => p_d1

        case (3)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_z1(i1))
          endif
          if (mode2 /= 0) call zcopy(i1,arr,1,p_z1,1)
          if (mode3 /= 0) call dvset(p_z1,1,2*i1,0d0)
          if (associated(s_pot%cp)) deallocate(s_pot%cp)
          s_pot%cp => p_z1

        case (4)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_pot%dlmwt)) deallocate(s_pot%dlmwt)
          s_pot%dlmwt => p_d1

        case (5)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_pot%GFr)
            call d2realloc(s_pot%GFr,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (associated(s_pot%GFr)) deallocate(s_pot%GFr)
          s_pot%GFr => p_d2

        case (6)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_pot%gibbs)) deallocate(s_pot%gibbs)
          s_pot%gibbs => p_d1

        case (7)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_pot%gma)
            call d2realloc(s_pot%gma,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (associated(s_pot%gma)) deallocate(s_pot%gma)
          s_pot%gma => p_d2

        case (8)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_pot%gmar)
            call z2realloc(s_pot%gmar,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (associated(s_pot%gmar)) deallocate(s_pot%gmar)
          s_pot%gmar => p_z2

        case (9)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_pot%grrme)) deallocate(s_pot%grrme)
          s_pot%grrme => p_d1

        case (10)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_pot%hab)
            call d2realloc(s_pot%hab,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (associated(s_pot%hab)) deallocate(s_pot%hab)
          s_pot%hab => p_d2

        case (11)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_pot%sab)
            call d2realloc(s_pot%sab,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (associated(s_pot%sab)) deallocate(s_pot%sab)
          s_pot%sab => p_d2

        case (12)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_pot%mad)) deallocate(s_pot%mad)
          s_pot%mad => p_d1

        case (13)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_pot%mxy)) deallocate(s_pot%mxy)
          s_pot%mxy => p_d1

        case (14)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_pot%palp)
            call z2realloc(s_pot%palp,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (associated(s_pot%palp)) deallocate(s_pot%palp)
          s_pot%palp => p_z2

        case (15)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_pot%papg)
            call z2realloc(s_pot%papg,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (associated(s_pot%papg)) deallocate(s_pot%papg)
          s_pot%papg => p_z2

        case (16)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_pot%pf)
            call z2realloc(s_pot%pf,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (associated(s_pot%pf)) deallocate(s_pot%pf)
          s_pot%pf => p_z2

        case (17)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_pot%pfnc)
            call z2realloc(s_pot%pfnc,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (associated(s_pot%pfnc)) deallocate(s_pot%pfnc)
          s_pot%pfnc => p_z2

        case (18)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_pot%pfr)
            call z2realloc(s_pot%pfr,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (associated(s_pot%pfr)) deallocate(s_pot%pfr)
          s_pot%pfr => p_z2

        case (19)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_pot%pmpol)) deallocate(s_pot%pmpol)
          s_pot%pmpol => p_d1

        case (20)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_pot%pnu)
            call d2realloc(s_pot%pnu,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (associated(s_pot%pnu)) deallocate(s_pot%pnu)
          s_pot%pnu => p_d2

        case (21)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_pot%pp)
            call d2realloc(s_pot%pp,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (associated(s_pot%pp)) deallocate(s_pot%pp)
          s_pot%pp => p_d2

        case (22)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_pot%ppn)) deallocate(s_pot%ppn)
          s_pot%ppn => p_d1

        case (23)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_pot%pprel)
            call d2realloc(s_pot%pprel,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (associated(s_pot%pprel)) deallocate(s_pot%pprel)
          s_pot%pprel => p_d2

        case (24)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_pot%pti)) deallocate(s_pot%pti)
          s_pot%pti => p_d1

        case (25)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_pot%dpf)
            call z2realloc(s_pot%dpf,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (associated(s_pot%dpf)) deallocate(s_pot%dpf)
          s_pot%dpf => p_z2

        case (26)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_pot%dpfr)
            call z2realloc(s_pot%dpfr,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (associated(s_pot%dpfr)) deallocate(s_pot%dpfr)
          s_pot%dpfr => p_z2

        case (27)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_pot%ddpf)
            call z2realloc(s_pot%ddpf,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (associated(s_pot%ddpf)) deallocate(s_pot%ddpf)
          s_pot%ddpf => p_z2

        case (28)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_pot%ddpfr)
            call z2realloc(s_pot%ddpfr,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (associated(s_pot%ddpfr)) deallocate(s_pot%ddpfr)
          s_pot%ddpfr => p_z2

        case (29)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_pot%dddpf)
            call z2realloc(s_pot%dddpf,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (associated(s_pot%dddpf)) deallocate(s_pot%dddpf)
          s_pot%dddpf => p_z2

        case (30)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_pot%qbyl)
            call d2realloc(s_pot%qbyl,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (associated(s_pot%qbyl)) deallocate(s_pot%qbyl)
          s_pot%qbyl => p_d2

        case (31)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_pot%qc)) deallocate(s_pot%qc)
          s_pot%qc => p_d1

        case (32)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_pot%qcorr)) deallocate(s_pot%qcorr)
          s_pot%qcorr => p_d1

        case (33)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_pot%qnu)
            call d2realloc(s_pot%qnu,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (associated(s_pot%qnu)) deallocate(s_pot%qnu)
          s_pot%qnu => p_d2

        case (34)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_pot%qnur)
            call d2realloc(s_pot%qnur,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (associated(s_pot%qnur)) deallocate(s_pot%qnur)
          s_pot%qnur => p_d2

        case (35)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_pot%qpp)) deallocate(s_pot%qpp)
          s_pot%qpp => p_d1

        case (36)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_pot%qt)) deallocate(s_pot%qt)
          s_pot%qt => p_d1

        case (37)
          if (mode0 /= 0) then
            allocate(p_rhat(i1))
          endif
          if (mode1 /= 0 .or. mode2 /= 0 .or. mode3 /= 0)
     .      call rx('ptr_pot: illegal allocation of rhat')
          if (associated(s_pot%rhat)) deallocate(s_pot%rhat)
          s_pot%rhat => p_rhat

        case (38)
          if (mode0 /= 0) then
            allocate(p_rhat(i1))
          endif
          if (mode1 /= 0 .or. mode2 /= 0 .or. mode3 /= 0)
     .      call rx('ptr_pot: illegal allocation of rnew')
          if (associated(s_pot%rnew)) deallocate(s_pot%rnew)
          s_pot%rnew => p_rhat

        case (39)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_pot%rhos)) deallocate(s_pot%rhos)
          s_pot%rhos => p_d1

        case (40)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_pot%rhrmx)) deallocate(s_pot%rhrmx)
          s_pot%rhrmx => p_d1

        case (41)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_pot%socscl)) deallocate(s_pot%socscl)
          s_pot%socscl => p_d1

        case (42)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_pot%sop)
            call d2realloc(s_pot%sop,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (associated(s_pot%sop)) deallocate(s_pot%sop)
          s_pot%sop => p_d2

        case (43)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_pot%shfac)
            call d2realloc(s_pot%shfac,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (associated(s_pot%shfac)) deallocate(s_pot%shfac)
          s_pot%shfac => p_d2

        case (44)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_pot%thetcl)
            call d2realloc(s_pot%thetcl,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (associated(s_pot%thetcl)) deallocate(s_pot%thetcl)
          s_pot%thetcl => p_d2

        case (45)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_pot%v0)
            call d2realloc(s_pot%v0,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (associated(s_pot%v0)) deallocate(s_pot%v0)
          s_pot%v0 => p_d2

        case (46)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_pot%vdif)) deallocate(s_pot%vdif)
          s_pot%vdif => p_d1

        case (47)
          if (mode0 /= 0 .or. mode1 /= 0) then
            allocate(p_d1(i1))
            if (associated(s_pot%ves) .and. mode1 /= 0) then
              i = min(size(s_pot%ves),i1)
              call dcopy(i,s_pot%ves,1,p_d1,1)
            endif
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_pot%ves)) deallocate(s_pot%ves)
          s_pot%ves => p_d1

        case (48)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_pot%vintr)) deallocate(s_pot%vintr)
          s_pot%vintr => p_d1

        case (49)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_pot%vrmax)
            call d2realloc(s_pot%vrmax,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (associated(s_pot%vrmax)) deallocate(s_pot%vrmax)
          s_pot%vrmax => p_d2

        case (50)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_pot%vshft)) deallocate(s_pot%vshft)
          s_pot%vshft => p_d1

        case (51)
          if (mode1 /= 0) then
            i123(1:2) = shape(s_pot%smpot)
            call z2realloc(s_pot%smpot,i123(1),i123(2),i1,i2)
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (associated(s_pot%smpot)) deallocate(s_pot%smpot)
          s_pot%smpot => p_z2

        case (52)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (associated(s_pot%smrho)) deallocate(s_pot%smrho)
          s_pot%smrho => p_z2

        case (53)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (associated(s_pot%smrout)) deallocate(s_pot%smrout)
          s_pot%smrout => p_z2

        case (54)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (associated(s_pot%smcv0)) deallocate(s_pot%smcv0)
          s_pot%smcv0 => p_z2

      end select

      end subroutine ptr_pot

      subroutine ptr_site(s_site,mode,elt,ib,i1,i2,arr)
C- Pointer allocation and management for s_site
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Co     Stored:    bxc cpawt omg omgn domg gc gcu gcorr sfvrtx vtxrel j0 tau pdos rho thet
Co                v0 v1
Cio    Passed to: *
Ci Inputs
Ci   mode : 0 Initialize pointers.  First call should use mode=0 and
Ci          ib=number of sites in structure
Ci        : 1s bit => allocate permanent array; associate s_site%elt to it
Ci        : 2s bit => re-allocate s_site%elt; reassociate s_site%elt
Ci        : Either bit 1 or bit 2 must be set, but not both.
Ci        : 4s bit => copy irr or darr to internal array
Ci        : 8s bit => initialize array to 0
Ci   elt  : string identifying which pointer to operate on:
Ci         Element      type      dim>1  internal index
Ci          bxc         real(8)              1
Ci          cpawt       real(8)              2
Ci          omg         complex(8)  2        3
Ci          omgn        complex(8)  2        4
Ci          domg        complex(8)  2        5
Ci          dmat        complex(8)  2        6
Ci          pfr         complex(8)  2        7
Ci          pfr         complex(8)  2        8
Ci          dpfr        complex(8)  2        9
Ci          ddpfr       complex(8)  2        10
Ci          gc          complex(8)  2        11
Ci          gcu         complex(8)  2        12
Ci          gcorr       complex(8)  2        13
Ci          gii         complex(8)  2        14
Ci          sfvrtx      complex(8)  2        15
Ci          vtxrel      complex(8)  2        16
Ci          j0          real(8)     2        17
Ci          tau         real(8)     2        18
Ci          pdos        real(8)     2        19
Ci          rho1        real(8)     2        20
Ci          rho2        real(8)     2        21
Ci          rhoc        real(8)     2        22
Ci          rho1x       real(8)     2        23
Ci          rho2x       real(8)     2        24
Ci          rhocx       real(8)     2        25
Ci          qhhl        real(8)     2        26
Ci          qhkl        real(8)     2        27
Ci          qkkl        real(8)     2        28
Ci          eqhhl       real(8)     2        29
Ci          eqhkl       real(8)     2        30
Ci          eqkkl       real(8)     2        31
Ci          sighh       real(8)     2        32
Ci          sighk       real(8)     2        33
Ci          sigkk       real(8)     2        34
Ci          tauhh       real(8)     2        35
Ci          tauhk       real(8)     2        36
Ci          taukk       real(8)     2        37
Ci          pihh        real(8)     2        38
Ci          pihk        real(8)     2        39
Ci          pikk        real(8)     2        40
Ci          sohh        real(8)     2        41
Ci          sohk        real(8)     2        42
Ci          sokk        real(8)     2        43
Ci          sighhx      real(8)     2        44
Ci          sighkx      real(8)     2        45
Ci          sigkkx      real(8)     2        46
Ci          tauhhx      real(8)     2        47
Ci          tauhkx      real(8)     2        48
Ci          taukkx      real(8)     2        49
Ci          pihhx       real(8)     2        50
Ci          pihkx       real(8)     2        51
Ci          pikkx       real(8)     2        52
Ci          thet        real(8)     2        53
Ci          v0          real(8)     2        54
Ci          v1          real(8)     2        55
Ci   ib   : site index: work with s_site(ib)%elt
Ci   i1   : leading dimension of pointer
Ci   i2   : 2nd dimension of pointer
Co Outputs
Co   This routine does the following depending on bits of mode
Co     - Internal pointer named by 'elt' is allocated OR
Co       the pointer is re-allocated
Co     - Data is copied from arr to it
Co     - s_site(ib)%elt is associated with internal array
Cl Local variables
Cl         :
Cr Remarks
Cr   This is a workaround for f90 pointer restrictions.
Cr   * Local pointers are allocated and held in common to enable
Cr     them to be preserved.
Cr   * A corresponding pointer in s_site is associated with this pointer
Cu Updates
Cu   11 Oct 12 First created
C ----------------------------------------------------------------------
      use pointers
      use structures
      implicit none
C ... Passed parameters
      character elt*(*)
      integer mode,ib,i1,i2
      integer arr(*)
C ... For structures
!      include 'structures.h'
      type(str_site):: s_site(ib)
C ... Local parameters
      integer i,k,isw,nelt,strnsiz,mode0,mode1,mode2,mode3
      parameter (nelt=55,strnsiz=6)
      integer eltidx(nelt)
      character*(nelt*strnsiz) eltlst

      save eltlst,eltidx
      data eltidx /1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     .  21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,
     .  42,43,44,45,46,47,48,49,50,51,52,53,54,55/
      data eltlst/' '/

      if (eltlst(1:1) == ' ' .or. ib < 0) then
        eltlst =
     .  'bxc   '//
     .  'cpawt '//
     .  'omg   '//
     .  'omgn  '//
     .  'domg  '//
     .  'dmat  '//
     .  'pfr   '//
     .  'pfra  '//
     .  'dpfr  '//
     .  'ddpfr '//
     .  'gc    '//
     .  'gcu   '//
     .  'gcorr '//
     .  'gii   '//
     .  'sfvrtx'//
     .  'vtxrel'//
     .  'j0    '//
     .  'tau   '//
     .  'pdos  '//
     .  'rho1  '//
     .  'rho2  '//
     .  'rhoc  '//
     .  'rho1x '//
     .  'rho2x '//
     .  'rhocx '//
     .  'qhhl  '//
     .  'qhkl  '//
     .  'qkkl  '//
     .  'eqhhl '//
     .  'eqhkl '//
     .  'eqkkl '//
     .  'sighh '//
     .  'sighk '//
     .  'sigkk '//
     .  'tauhh '//
     .  'tauhk '//
     .  'taukk '//
     .  'pihh  '//
     .  'pihk  '//
     .  'pikk  '//
     .  'sohh  '//
     .  'sohk  '//
     .  'sokk  '//
     .  'sighhx'//
     .  'sighkx'//
     .  'sigkkx'//
     .  'tauhhx'//
     .  'tauhkx'//
     .  'taukkx'//
     .  'pihhx '//
     .  'pihkx '//
     .  'pikkx '//
     .  'thet  '//
     .  'v0    '//
     .  'v1    '
        do  i = 1, iabs(ib)
          nullify(
     .      s_site(i)%bxc,
     .      s_site(i)%cpawt,
     .      s_site(i)%omg,
     .      s_site(i)%omgn,
     .      s_site(i)%domg,
     .      s_site(i)%dmat,
     .      s_site(i)%pfr,
     .      s_site(i)%pfra,
     .      s_site(i)%dpfr,
     .      s_site(i)%ddpfr,
     .      s_site(i)%gc,
     .      s_site(i)%gcu,
     .      s_site(i)%gcorr,
     .      s_site(i)%gii,
     .      s_site(i)%sfvrtx,
     .      s_site(i)%vtxrel,
     .      s_site(i)%j0,
     .      s_site(i)%tau,
     .      s_site(i)%pdos,
     .      s_site(i)%rho1,
     .      s_site(i)%rho2,
     .      s_site(i)%rhoc,
     .      s_site(i)%rho1x,
     .      s_site(i)%rho2x,
     .      s_site(i)%rhocx,
     .      s_site(i)%qhhl,
     .      s_site(i)%qhkl,
     .      s_site(i)%qkkl,
     .      s_site(i)%eqhhl,
     .      s_site(i)%eqhkl,
     .      s_site(i)%eqkkl,
     .      s_site(i)%sighh,
     .      s_site(i)%sighk,
     .      s_site(i)%sigkk,
     .      s_site(i)%tauhh,
     .      s_site(i)%tauhk,
     .      s_site(i)%taukk,
     .      s_site(i)%pihh,
     .      s_site(i)%pihk,
     .      s_site(i)%pikk,
     .      s_site(i)%sohh,
     .      s_site(i)%sohk,
     .      s_site(i)%sokk,
     .      s_site(i)%sighhx,
     .      s_site(i)%sighkx,
     .      s_site(i)%sigkkx,
     .      s_site(i)%tauhhx,
     .      s_site(i)%tauhkx,
     .      s_site(i)%taukkx,
     .      s_site(i)%pihhx,
     .      s_site(i)%pihkx,
     .      s_site(i)%pikkx,
     .      s_site(i)%thet,
     .      s_site(i)%v0,
     .      s_site(i)%v1)
        enddo
C       nullify(p_i1,p_i2,p_d1,p_d2,p_z1,p_z2)
      endif

      if (mode == 0) return
      mode0 = mod(mode,2)
      mode1 = mod(mode/2,2)
      mode2 = mod(mode/4,2)
      mode3 = mod(mode/8,2)

      k = index(eltlst,elt)
      if (k == 0 .or. mod(k-1,strnsiz) /= 0)
     .  call rxs('ptr_site:  unknown element ',elt)

      i = (k-1)/6 + 1
      isw = eltidx(i)

C     No pointers coded for reallocation
      if (mode1 /= 0)
     .  call rx('ptr_site: mode1 not implemented for '//elt)
      if (mode0 /= 1 .and. mode1 /= 1 .and. mode0*mode1 /= 0)
     .  call rx('ptr_site: improper mode')

      select case (isw)

C   ... Integer pointers
        case (1)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (associated(s_site(ib)%bxc)) deallocate(s_site(ib)%bxc)
          s_site(ib)%bxc => p_d1
        case (2)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (associated(s_site(ib)%cpawt))
     .      deallocate(s_site(ib)%cpawt)
          s_site(ib)%cpawt => p_d1
        case (3)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (associated(s_site(ib)%omg)) deallocate(s_site(ib)%omg)
          s_site(ib)%omg => p_z2
        case (4)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (associated(s_site(ib)%omgn)) deallocate(s_site(ib)%omgn)
          s_site(ib)%omgn => p_z2
        case (5)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (associated(s_site(ib)%domg)) deallocate(s_site(ib)%domg)
          s_site(ib)%domg => p_z2
        case (6)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (associated(s_site(ib)%dmat)) deallocate(s_site(ib)%dmat)
          s_site(ib)%dmat => p_z2
        case (7)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (associated(s_site(ib)%pfr)) deallocate(s_site(ib)%pfr)
          s_site(ib)%pfr => p_z2
        case (8)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (associated(s_site(ib)%pfra)) deallocate(s_site(ib)%pfra)
          s_site(ib)%pfra => p_z2
        case (9)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (associated(s_site(ib)%dpfr)) deallocate(s_site(ib)%dpfr)
          s_site(ib)%dpfr => p_z2
        case (10)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (associated(s_site(ib)%ddpfr)) deallocate(s_site(ib)%ddpfr)
          s_site(ib)%ddpfr => p_z2
        case (11)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (associated(s_site(ib)%gc)) deallocate(s_site(ib)%gc)
          s_site(ib)%gc => p_z2
        case (12)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (associated(s_site(ib)%gcu)) deallocate(s_site(ib)%gcu)
          s_site(ib)%gcu => p_z2
        case (13)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (associated(s_site(ib)%gcorr))
     .      deallocate(s_site(ib)%gcorr)
          s_site(ib)%gcorr => p_z2
        case (14)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (associated(s_site(ib)%gii))
     .      deallocate(s_site(ib)%gii)
          s_site(ib)%gii => p_z2
        case (15)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (associated(s_site(ib)%sfvrtx))
     .      deallocate(s_site(ib)%sfvrtx)
          s_site(ib)%sfvrtx => p_z2
        case (16)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_z2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_z2,1,2*i1*i2,0d0)
          if (mode2 /= 0) call zcopy(i1*i2,arr,1,p_z2,1)
          if (associated(s_site(ib)%vtxrel))
     .      deallocate(s_site(ib)%vtxrel)
          s_site(ib)%vtxrel => p_z2
        case (17)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%j0)) deallocate(s_site(ib)%j0)
          s_site(ib)%j0 => p_d2
        case (18)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%tau)) deallocate(s_site(ib)%tau)
          s_site(ib)%tau => p_d2
        case (19)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%pdos)) deallocate(s_site(ib)%pdos)
          s_site(ib)%pdos => p_d2
        case (20)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%rho1)) deallocate(s_site(ib)%rho1)
          s_site(ib)%rho1 => p_d2
        case (21)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%rho2)) deallocate(s_site(ib)%rho2)
          s_site(ib)%rho2 => p_d2
        case (22)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%rhoc)) deallocate(s_site(ib)%rhoc)
          s_site(ib)%rhoc => p_d2
        case (23)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%rho1x)) deallocate(s_site(ib)%rho1x)
          s_site(ib)%rho1x => p_d2
        case (24)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%rho2x)) deallocate(s_site(ib)%rho2x)
          s_site(ib)%rho2x => p_d2
        case (25)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%rhocx)) deallocate(s_site(ib)%rhocx)
          s_site(ib)%rhocx => p_d2
        case (26)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%qhhl)) deallocate(s_site(ib)%qhhl)
          s_site(ib)%qhhl => p_d2
        case (27)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%qhkl)) deallocate(s_site(ib)%qhkl)
          s_site(ib)%qhkl => p_d2
        case (28)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%qkkl)) deallocate(s_site(ib)%qkkl)
          s_site(ib)%qkkl => p_d2
        case (29)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%eqhhl)) deallocate(s_site(ib)%eqhhl)
          s_site(ib)%eqhhl => p_d2
        case (30)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%eqhkl)) deallocate(s_site(ib)%eqhkl)
          s_site(ib)%eqhkl => p_d2
        case (31)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%eqkkl)) deallocate(s_site(ib)%eqkkl)
          s_site(ib)%eqkkl => p_d2
        case (32)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%sighh)) deallocate(s_site(ib)%sighh)
          s_site(ib)%sighh => p_d2
        case (33)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%sighk)) deallocate(s_site(ib)%sighk)
          s_site(ib)%sighk => p_d2
        case (34)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%sigkk)) deallocate(s_site(ib)%sigkk)
          s_site(ib)%sigkk => p_d2
        case (35)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%tauhh)) deallocate(s_site(ib)%tauhh)
          s_site(ib)%tauhh => p_d2
        case (36)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%tauhk)) deallocate(s_site(ib)%tauhk)
          s_site(ib)%tauhk => p_d2
        case (37)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%taukk)) deallocate(s_site(ib)%taukk)
          s_site(ib)%taukk => p_d2
        case (38)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%pihh)) deallocate(s_site(ib)%pihh)
          s_site(ib)%pihh => p_d2
        case (39)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%pihk)) deallocate(s_site(ib)%pihk)
          s_site(ib)%pihk => p_d2
        case (40)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%pikk)) deallocate(s_site(ib)%pikk)
          s_site(ib)%pikk => p_d2
        case (41)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%sohh)) deallocate(s_site(ib)%sohh)
          s_site(ib)%sohh => p_d2
        case (42)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%sohk)) deallocate(s_site(ib)%sohk)
          s_site(ib)%sohk => p_d2
        case (43)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%sokk)) deallocate(s_site(ib)%sokk)
          s_site(ib)%sokk => p_d2
        case (44)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%sighhx))
     .      deallocate(s_site(ib)%sighhx)
          s_site(ib)%sighhx => p_d2
        case (45)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%sighkx))
     .      deallocate(s_site(ib)%sighkx)
          s_site(ib)%sighkx => p_d2
        case (46)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%sigkkx))
     .      deallocate(s_site(ib)%sigkkx)
          s_site(ib)%sigkkx => p_d2
        case (47)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%tauhhx))
     .      deallocate(s_site(ib)%tauhhx)
          s_site(ib)%tauhhx => p_d2
        case (48)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%tauhkx))
     .      deallocate(s_site(ib)%tauhkx)
          s_site(ib)%tauhkx => p_d2
        case (49)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%taukkx))
     .      deallocate(s_site(ib)%taukkx)
          s_site(ib)%taukkx => p_d2
        case (50)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%pihhx)) deallocate(s_site(ib)%pihhx)
          s_site(ib)%pihhx => p_d2
        case (51)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%pihkx)) deallocate(s_site(ib)%pihkx)
          s_site(ib)%pihkx => p_d2
        case (52)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%pikkx)) deallocate(s_site(ib)%pikkx)
          s_site(ib)%pikkx => p_d2
        case (53)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%thet)) deallocate(s_site(ib)%thet)
          s_site(ib)%thet => p_d2
        case (54)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%v0)) deallocate(s_site(ib)%v0)
          s_site(ib)%v0 => p_d2
        case (55)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d2(i1,i2))
          endif
          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
          if (associated(s_site(ib)%v1)) deallocate(s_site(ib)%v1)
          s_site(ib)%v1 => p_d2

      end select

      end subroutine ptr_site

      subroutine ptr_spec(s_spec,mode,elt,is,i1,i2,arr)
C- Pointer allocation and management for s_spec
C ----------------------------------------------------------------------
Cio Structures
Cio  s_spec :struct for species-specific data; see structures.h
Co     Stored:    rhoc
Cio    Passed to: *
Ci Inputs
Ci   mode : 0 Initialize pointers.  First call should use mode=0
Ci          with is=number of elements in structure
Ci        : 1s bit => allocate permanent array; associate s_spec%elt to it
Ci        : 2s bit => re-allocate s_spec%elt; reassociate s_spec%elt
Ci        : Either bit 1 or bit 2 must be set, but not both.
Ci        : 4s bit => copy irr or darr to internal array
Ci        : 8s bit => initialize array to 0
Ci   elt  : string identifying which pointer to operate on:
Ci         Element      type      dim>1  internal index
Ci          rhoc        real(8)             1
Ci   is   : species index: work with s_spec(is)%elt
Ci   i1   : leading dimension of pointer
Ci   i2   : 2nd dimension of pointer
Co Outputs
Co   This routine does the following depending on bits of mode
Co     - Internal pointer named by 'elt' is allocated OR
Co       the pointer is re-allocated
Co     - Data is copied from arr to it
Co     - s_spec(is)%elt is associated with internal array
Cl Local variables
Cl         :
Cr Remarks
Cr   This is a workaround for f90 pointer restrictions.
Cr   * Local pointers are allocated and held in common to enable
Cr     them to be preserved.
Cr   * A corresponding pointer in s_spec is associated with this pointer
Cu Updates
Cu   11 Oct 12  First created
C ----------------------------------------------------------------------
      use pointers
      use structures
      implicit none
C ... Passed parameters
      character elt*(*)
      integer mode,is,i1,i2
      integer arr(*)
C ... For structures
!      include 'structures.h'
      type(str_spec):: s_spec(is)
C ... Local parameters
      integer i,k,isw,nelt,strnsiz,mode0,mode1,mode2,mode3
      parameter (nelt=1,strnsiz=6)
      integer eltidx(nelt)
      character*(nelt*strnsiz) eltlst

      save eltlst,eltidx
      data eltidx /1/
      data eltlst/' '/
      if (eltlst(1:1) == ' ') then
        eltlst =
     .  'rhoc  '
        do  i = 1, is
          nullify(s_spec(i)%rhoc)
        enddo
C       nullify(p_i1,p_i2,p_d1,p_d2,p_z1,p_z2)
      endif

      if (mode == 0) return
      mode0 = mod(mode,2)
      mode1 = mod(mode/2,2)
      mode2 = mod(mode/4,2)
      mode3 = mod(mode/8,2)

      k = index(eltlst,elt)
      if (k == 0 .or. mod(k-1,strnsiz) /= 0)
     .  call rxs('ptr_spec:  unknown element ',elt)

      i = (k-1)/6 + 1
      isw = eltidx(i)

C     No pointers coded for reallocation
      if (mode1 /= 0)
     .  call rx('ptr_spec: mode1 not implemented for '//elt)
      if (mode0 /= 1 .and. mode1 /= 1 .and. mode0*mode1 /= 0)
     .  call rx('ptr_mode1: improper mode')

      select case (isw)

        case (1)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (associated(s_spec(is)%rhoc)) deallocate(s_spec(is)%rhoc)
          s_spec(is)%rhoc => p_d1

      end select

      end subroutine ptr_spec

      subroutine ptr_str(s_str,mode,elt,i1,i2,arr)
C- Pointer allocation and management for s_str
C ----------------------------------------------------------------------
Cio Structures
Cio  s_str  :struct for parameters for screened strux; see structures.h
Co     Stored:    iax npr adot alp s sdot
Cio    Passed to: *
Ci Inputs
Ci   mode : 0 Initialize pointers.  First call should use mode=0.
Ci        : 1s bit => allocate permanent array; associate s_str%elt to it
Ci        : 2s bit => re-allocate s_str%elt; reassociate s_str%elt
Ci        : Either bit 1 or bit 2 must be set, but not both.
Ci        : 4s bit => copy irr or darr to internal array
Ci        : 8s bit => initialize array to 0
Ci   elt  : string identifying which pointer to operate on:
Ci         Element      type      dim>1  internal index
Ci          iax         integer              1
Ci          npr         integer              2
Ci          adot        real(8)              3
Ci          alp         real(8)              4
Ci          s           real(8)              5
Ci          sdot        real(8)              6
Ci   i1   : leading dimension of pointer
Ci   i2   : 2nd dimension of pointer (not used so far)
Co Outputs
Co   This routine does the following depending on bits of mode
Co     - Internal pointer named by 'elt' is allocated OR
Co       the pointer is re-allocated
Co     - Data is copied from arr to it
Co     - s_str%elt is associated with internal array
Cl Local variables
Cl         :
Cr Remarks
Cr   This is a workaround for f90 pointer restrictions.
Cr   * Local pointers are allocated and held in common to enable
Cr     them to be preserved.
Cr   * A corresponding pointer in s_str is associated with this pointer
Cu Updates
Cu   11 Oct 12  First created
C ----------------------------------------------------------------------
      use pointers
      use structures
      implicit none
C ... Passed parameters
      character elt*(*)
      integer mode,i1,i2
      integer arr(*)
C ... For structures
!      include 'structures.h'
      type(str_str)::  s_str
C ... Local parameters
      integer i,k,isw,nelt,strnsiz,mode0,mode1,mode2,mode3
      parameter (nelt=6,strnsiz=6)
      integer, parameter :: eltidx(nelt) = (/1,2,3,4,5,6/)
      character(len=nelt*strnsiz), parameter :: eltlst
     .                         = 'iax   npr   adot  alp   s     sdot  '

      if (mode == 0) return
      mode0 = mod(mode,2)
      mode1 = mod(mode/2,2)
      mode2 = mod(mode/4,2)
      mode3 = mod(mode/8,2)

      k = index(eltlst,elt)
      if (k == 0 .or. mod(k-1,strnsiz) /= 0)
     .  call rxs('ptr_str:  unknown element ',elt)

      i = (k-1)/6 + 1
      isw = eltidx(i)

C     No reallocation implemented
      if (mode1 /= 0)
     .  call rx('ptr_str: mode1 not implemented for '//elt)
      if (mode0 /= 1 .and. mode1 /= 1 .and. mode0*mode1 /= 0)
     .  call rx('ptr_str: improper mode')

      select case (isw)

        case (1)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_i1(i1))
          endif
          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
          if (mode3 /= 0) call ivset(p_i1,1,i1,0d0)
          if (associated(s_str%iax)) deallocate(s_str%iax)
          s_str%iax => p_i1

        case (2)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_i1(i1))
          endif
          if (mode2 /= 0) call icopy(i1,arr,1,p_i1,1)
          if (mode3 /= 0) call ivset(p_i1,1,i1,0d0)
          if (associated(s_str%npr)) deallocate(s_str%npr)
          s_str%npr => p_i1

        case (3)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_str%adot)) deallocate(s_str%adot)
          s_str%adot => p_d1

        case (4)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_str%alph)) deallocate(s_str%alph)
          s_str%alph => p_d1

        case (5)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_str%s)) deallocate(s_str%s)
          s_str%s => p_d1

        case (6)
          if (mode1 /= 0) then
          elseif (mode0 /= 0) then
            allocate(p_d1(i1))
          endif
          if (mode2 /= 0) call dcopy(i1,arr,1,p_d1,1)
          if (mode3 /= 0) call dvset(p_d1,1,i1,0d0)
          if (associated(s_str%sdot)) deallocate(s_str%sdot)
          s_str%sdot => p_d1

      end select

      end subroutine ptr_str

C      subroutine ptr_rhat(s_rhat,mode,elt,ib,i1,i2,arr)
CC- Pointer allocation and management for s_rhat
CC ----------------------------------------------------------------------
CCi Inputs
CCi   mode : 0 Initialize pointers.  First call should use mode=0 and
CCi          ib=number of sites in structure
CCi        : 1s bit => allocate permanent array; associate s_rhat%elt to it
CCi        : 2s bit => re-allocate s_rhat%elt; reassociate s_rhat%elt
CCi        : Either bit 1 or bit 2 must be set, but not both.
CCi        : 4s bit => copy irr or darr to internal array
CCi        : 8s bit => initialize array to 0
CCi   elt  : string identifying which pointer to operate on:
CCi         Element      type      dim>1  internal index
CCi          rho1        real(8)     2        1
CCi          rho2        real(8)     2        2
CCi          rhoc        real(8)     2        3
CCi          eqhhl       real(8)     2        4
CCi          eqhkl       real(8)     2        5
CCi          eqkkl       real(8)     2        6
CCi          qhhl        real(8)     2        7
CCi          qhkl        real(8)     2        8
CCi          qkkl        real(8)     2        9
CCi   ib   : site index: work with s_rhat(ib)%elt
CCi   i1   : leading dimension of pointer
CCi   i2   : 2nd dimension of pointer
CCo Outputs
CCo   This routine does the following depending on bits of mode
CCo     - Internal pointer named by 'elt' is allocated OR
CCo       the pointer is re-allocated
CCo     - Data is copied from arr to it
CCo     - s_rhat(ib)%elt is associated with internal array
CCl Local variables
CCl         :
CCr Remarks
CCr   This is a workaround for f90 pointer restrictions.
CCr   * Local pointers are allocated and held in common to enable
CCr     them to be preserved.
CCr   * A corresponding pointer in s_rhat is associated with this pointer
CCu Updates
CCu   11 Oct 12  First created
CC ----------------------------------------------------------------------
C      use pointers
C      implicit none
CC ... Passed parameters
C      character elt*(*)
C      integer mode,ib,i1,i2
C      integer arr(*)
CC ... For structures
C      include 'structures.h'
C      type(str_rhat):: s_rhat(ib)
CC ... Local parameters
C      integer i,k,isw,nelt,strnsiz,mode0,mode1,mode2,mode3
C      parameter (nelt=9,strnsiz=6)
C      integer eltidx(nelt)
C      character*(nelt*strnsiz) eltlst
C
C      save eltlst,eltidx
C      data eltidx /1,2,3,4,5,6,7,8,9/
C      data eltlst/' '/
C
C      if (eltlst(1:1) == ' ' .or. ib < 0) then
C        eltlst =
C     .  'rho1  '//
C     .  'rho2  '//
C     .  'rhoc  '//
C     .  'eqhhl '//
C     .  'eqhkl '//
C     .  'eqkkl '//
C     .  'qhhl  '//
C     .  'qhkl  '//
C     .  'qkkl  '
C        do  i = 1, ib
C          nullify(s_rhat(i)%rho1,s_rhat(i)%rho2,s_rhat(i)%rhoc,
C     .      s_rhat(i)%eqhhl,s_rhat(i)%eqhkl,s_rhat(i)%eqkkl,
C     .      s_rhat(i)%qhhl,s_rhat(i)%qhkl,s_rhat(i)%qkkl)
C        enddo
CC       nullify(p_i1,p_i2,p_d1,p_d2,p_z1,p_z2)
C      endif
C
C      if (mode == 0) return
C      mode0 = mod(mode,2)
C      mode1 = mod(mode/2,2)
C      mode2 = mod(mode/4,2)
C      mode3 = mod(mode/8,2)
C
C      k = index(eltlst,elt)
C      if (k == 0 .or. mod(k-1,strnsiz) /= 0)
C     .  call rxs('ptr_rhat:  unknown element ',elt)
C
C      i = (k-1)/6 + 1
C      isw = eltidx(i)
C
CC     No pointers coded for reallocation
C      if (mode1 /= 0)
C     .  call rx('ptr_rhat: mode1 not implemented for '//elt)
C      if (mode0 /= 1 .and. mode1 /= 1 .and. mode0*mode1 /= 0)
C     .  call rx('ptr_rhat: improper mode')
C
C      select case (isw)
C
CC   ... Integer pointers
C        case (1)
C          if (mode1 /= 0) then
C          elseif (mode0 /= 0) then
C            allocate(p_d2(i1,i2))
C          endif
C          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
C          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
C          if (associated(s_rhat(ib)%rho1)) deallocate(s_rhat(ib)%rho1)
C          s_rhat(ib)%rho1 => p_d2
C        case (2)
C          if (mode1 /= 0) then
C          elseif (mode0 /= 0) then
C            allocate(p_d2(i1,i2))
C          endif
C          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
C          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
C          if (associated(s_rhat(ib)%rho2)) deallocate(s_rhat(ib)%rho2)
C          s_rhat(ib)%rho2 => p_d2
C        case (3)
C          if (mode1 /= 0) then
C          elseif (mode0 /= 0) then
C            allocate(p_d2(i1,i2))
C          endif
C          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
C          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
C          if (associated(s_rhat(ib)%rhoc)) deallocate(s_rhat(ib)%rhoc)
C          s_rhat(ib)%rhoc => p_d2
C        case (4)
C          if (mode1 /= 0) then
C          elseif (mode0 /= 0) then
C            allocate(p_d2(i1,i2))
C          endif
C          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
C          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
C          if (associated(s_rhat(ib)%eqhhl)) deallocate(s_rhat(ib)%eqhhl)
C          s_rhat(ib)%eqhhl => p_d2
C        case (5)
C          if (mode1 /= 0) then
C          elseif (mode0 /= 0) then
C            allocate(p_d2(i1,i2))
C          endif
C          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
C          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
C          if (associated(s_rhat(ib)%eqhkl)) deallocate(s_rhat(ib)%eqhkl)
C          s_rhat(ib)%eqhkl => p_d2
C        case (6)
C          if (mode1 /= 0) then
C          elseif (mode0 /= 0) then
C            allocate(p_d2(i1,i2))
C          endif
C          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
C          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
C          if (associated(s_rhat(ib)%eqkkl)) deallocate(s_rhat(ib)%eqkkl)
C          s_rhat(ib)%eqkkl => p_d2
C        case (7)
C          if (mode1 /= 0) then
C          elseif (mode0 /= 0) then
C            allocate(p_d2(i1,i2))
C          endif
C          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
C          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
C          if (associated(s_rhat(ib)%qhhl)) deallocate(s_rhat(ib)%qhhl)
C          s_rhat(ib)%qhhl => p_d2
C        case (8)
C          if (mode1 /= 0) then
C          elseif (mode0 /= 0) then
C            allocate(p_d2(i1,i2))
C          endif
C          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
C          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
C          if (associated(s_rhat(ib)%qhkl)) deallocate(s_rhat(ib)%qhkl)
C          s_rhat(ib)%qhkl => p_d2
C        case (9)
C          if (mode1 /= 0) then
C          elseif (mode0 /= 0) then
C            allocate(p_d2(i1,i2))
C          endif
C          if (mode3 /= 0) call dvset(p_d2,1,i1*i2,0d0)
C          if (mode2 /= 0) call dcopy(i1*i2,arr,1,p_d2,1)
C          if (associated(s_rhat(ib)%qkkl)) deallocate(s_rhat(ib)%qkkl)
C          s_rhat(ib)%qkkl => p_d2
C
C      end select
C
C      end subroutine ptr_rhat

      subroutine i1alloc(p_iarr,n,linit)
C- Integer 1-dimensional pointer array allocated and pointer associated
C ----------------------------------------------------------------------
Ci Inputs
Ci   p_iarr : pointer to original array
Ci   n      : array dimension
Ci   linit  : if nonzero, initialize to 0
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr   To be f90 compliant, the calling routine should contain this INTERFACE:
Cr      INTERFACE
Cr      subroutine i1alloc(p_iarr,n,linit)
Cr      integer n,linit
Cr      integer, pointer :: p_iarr(:)
Cr      end subroutine i1alloc
Cr      END INTERFACE
Cu Updates
Cu   10 Jul 13  First created
      use pointers
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n,linit
      integer, pointer :: p_iarr(:)
C ... Local parameters

      allocate(p_i1(n))         ! Allocate new array; copy data there
      if (linit /= 0) call iinit(p_i1,n)
      p_iarr => p_i1            ! Associate pointer to array
      end

      subroutine i1realloc(p_iarr,nold,nnew)
C- Integer 1-dimensional pointer array re-allocated
C ----------------------------------------------------------------------
Ci Inputs
Ci   p_iarr : pointer to original array
Ci   nold   : original dimension
Ci   nnew   : reallocated dimension
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   24 Oct 12
C ----------------------------------------------------------------------
      use pointers
      implicit none
C ... Passed parameters
      integer nold,nnew,p_iarr(nold)
C ... Local parameters

      allocate(p_i1(nnew)) ! Allocate new array; copy data there
      call iinit(p_i1,nnew)
      if (nold > 0) call icopy(min(nold,nnew),p_iarr,1,p_i1,1)

      end

      subroutine i2realloc(p_iarr,n1old,n2old,n1new,n2new)
C- integer 2-dimensional pointer array re-allocated
C ----------------------------------------------------------------------
Ci Inputs
Ci   p_iarr : pointer to original array
Ci   n1old  : original dimension
Ci   n1new  : reallocated dimension
Co Outputs
Cr Remarks
Cr   Warning ... never been checked!
Cu Updates
Cu   24 Oct 12
C ----------------------------------------------------------------------
      use pointers
      implicit none
C ... Passed parameters
      integer n1old,n1new,n2old,n2new
      integer p_iarr(n1old,n2old)
C ... Local parameters
      integer i,j

      allocate(p_i2(n1new,n2new))  ! Allocate new array; copy data there
      call iinit(p_i2,n1new*n2new)
      do  i = 1, min(n1old,n1new)
        do  j = 1, min(n2old,n2new)
          p_i2(i,j) = p_iarr(i,j)
        enddo
      enddo

      end

      subroutine d1alloc(p_darr,n,linit)
C- double precision 1-dimensional pointer array allocated and pointer associated
C ----------------------------------------------------------------------
Ci Inputs
Ci   p_darr : pointer to original array
Ci   n      : array dimension
Ci   linit  : if nonzero, initialize to 0
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr   To be f90 compliant, the calling routine should contain this INTERFACE:
Cr      INTERFACE
Cr      subroutine d1alloc(p_darr,n,linit)
Cr      integer n,linit
Cr      real(8), pointer :: p_darr(:)
Cr      end subroutine d1alloc
Cr      END INTERFACE
Cu Updates
Cu   10 Jul 13  First created
C ----------------------------------------------------------------------
      use pointers
      implicit none
C ... Passed parameters
      integer n,linit
      real(8), pointer :: p_darr(:)
C ... Local parameters

      allocate(p_d1(n))         ! Allocate new array; copy data there
      if (linit /= 0) call iinit(p_d1,n)
      p_darr => p_d1            ! Associate pointer to array
      end

      subroutine d1realloc(p_darr,nold,nnew)
C- Double precision 1-dimensional pointer array re-allocated
C ----------------------------------------------------------------------
Ci Inputs
Ci   p_darr : pointer to original array
Ci   nold   : original dimension
Ci   nnew   : reallocated dimension
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   24 Oct 12
C ----------------------------------------------------------------------
      use pointers
      implicit none
C ... Passed parameters
      integer nold,nnew
      double precision p_darr(nold)
C ... Local parameters

      allocate(p_d1(nnew))      ! Allocate new array; copy data there
      call dpzero(p_d1,nnew)
      if (nold > 0) call dcopy(min(nold,nnew),p_darr,1,p_d1,1)

      end

      subroutine d2realloc(p_darr,n1old,n2old,n1new,n2new)
C- Double precision 2-dimensional pointer array re-allocated
C ----------------------------------------------------------------------
Ci Inputs
Ci   p_darr : pointer to original array
Ci   n1old  : original dimension
Ci   n1new  : reallocated dimension
Co Outputs
Cr Remarks
Cr
Cu Updates
Cu   24 Oct 12
C ----------------------------------------------------------------------
      use pointers
      implicit none
C ... Passed parameters
      integer n1old,n1new,n2old,n2new
      double precision p_darr(n1old,n2old)

      allocate(p_d2(n1new,n2new))  ! Allocate new array; copy data there
      call dpzero(p_d2,n1new*n2new)
      call dmcpy(p_darr,n1old,1,p_d2,n1new,1,min(n1old,n1new),min(n2old,n2new))

      end

      subroutine z2realloc(p_zarr,n1old,n2old,n1new,n2new)
C- Double complex 2-dimensional pointer array re-allocated
C ----------------------------------------------------------------------
Ci Inputs
Ci   p_zarr : pointer to original array
Ci   n1old  : original dimension
Ci   n1new  : reallocated dimension
Co Outputs
Co   p_z2   : pointer p_z2 (in common) is allocated
Co          : p_zarr (or a subblock of it) is copied to p_z2
Cr Remarks
Cu Updates
Cu   24 Oct 12
C ----------------------------------------------------------------------
      use pointers
      implicit none
C ... Passed parameters
      integer n1old,n1new,n2old,n2new
      complex(8) :: p_zarr(n1old,n2old)

      allocate(p_z2(n1new,n2new))  ! Allocate new array; copy data there
      call dpzero(p_z2,2*n1new*n2new)
      call zmcpy('N',p_zarr,n1old,1,p_z2,n1new,1,min(n1old,n1new),min(n2old,n2new))

      end

C Tests a pointer allocator
C      subroutine fmain
C      implicit none
C      double precision parr(3,4)
CC ... For structures
C      include 'structures.h'
C      type(str_lat)::  s_lat
C      real ran1
C      integer i
C
C      call ran1in(1)
C      do  i = 1, 3*4
C        parr(i,1) = ran1()
C      enddo
C
C      call ptr_lat(s_lat,4+1,'pos',3,2,0,parr)
C
C      do  i = 1, 3*4
C        parr(i,1) = ran1()
C      enddo
CC     This one re-allocates s_lat preserving contents
C      call ptr_lat(s_lat,2,'pos',3,4,0,0)
C      print *, s_lat%pos
CC     This one re-allocates s_lat updating contents
C      call ptr_lat(s_lat,4+2,'pos',3,4,0,parr)
C      print *, s_lat%pos
C
C      end
C
C Test str_pack
C      subroutine fmain
C      implicit none
C      integer str_pack,i
CC ... For structures
C      include 'structures.h'
C      integer,parameter :: nstrn=10
C      type(str_strn) :: s_strn(nstrn)
C      character*(100) strn
C
C      do  i = 1, nstrn
C        nullify(s_strn(i)%strn)
C      enddo
C      i = str_pack('gemb',1,s_strn,'my gemb string')
C      i = str_pack('sxopt',1,s_strn,'my sxopt string')
C      i = str_pack('gemb',1,s_strn,'my second gemb string')
C      i = str_pack('syml',1,s_strn,'syml string: 41,2,3,4')
C
C      i = str_pack('gemb',-1,s_strn,strn)
C      print "(' gemb unpacked: pos =',i3,'  string = ',a)",i, trim(strn)
C      i = str_pack('sxopt',-1,s_strn,strn)
C      print "(' sxopt unpacked: pos =',i3,'  string = ',a)",i, trim(strn)
C      i = str_pack('syml',-1,s_strn,strn)
C      print "(' syml unpacked: pos =',i3,'  string = ',a)",i, trim(strn)
C
C      end
