      subroutine pairs(nbas,nbasp,alat,plat,rmax,baspp,ipsp,nd,iltab,
     .  pltab,nttab,ntab,iax,mxcsiz)
C- Allocate memory for and create neighbor table for a crystal
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis (input)
Ci   nbasp :size of padded basis (layer programs)
Ci          nbasp = nbas + nbas(left bulk) + nbas(right bulk)
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat (input)
Ci   rmax  :maximum range for connecting vector, in a.u.
Ci          All connecting vectors with length < rmax(i)+rmax(j)
Ci          are retained.  rmax may be a scalar, a species-dependent
Ci          array, or a site-dependent array, depending on ipsp(1);
Ci          see description of ipsp
Ci   baspp :basis vectors, doubly padded for planar geometry
Ci   ipsp  :index to which species each site belongs, for padded basis;
Ci          identifies which rmax is associated with each site. NB:
Ci          ipsp(1) = -1 => rmax is a global scalar, independent of site
Ci          ipsp(1) =  0 => rmax is site-, not species-dependent.
Ci          In either of these cases, ipsp is not used.
Ci   nd    :number of dimensions for which periodic boundary conditions
Ci          are used
Ci   iltab :iltab<0, has no effect.  Otherwise, see pltabp.
Ci   pltab :include only pairs for which pltabp(jb)-pltabp(ib) <= iltab
Ci         :Only used if iltab>=0
Ci   mxcsiz:if nonzero, use in place of internal formula for mxnbr
Co Outputs
Co   nttab :total number of pairs in neighbor table and iax
Co   ntab  :ntab array; see pairc below where it is generated
Co   iax   :array with pair info; see pairc, below, where it is generated
Co   mxcsiz:size of the largest cluster encountered
Cu Updates
Cu   16 Apr 14 (DP) replace s_iax with iax
Cu   08 Jul 13 Replace f77 pointers with f90 ones
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer :: nbas,nbasp,nttab,nd,iltab,mxcsiz
      integer :: ntab(nbasp+1),ipsp(*),pltab(*)
      real(8) :: alat,plat(9),rmax(*),baspp(3,*)
      integer, intent(inout), allocatable :: iax(:)
C ... Dynamically allocated local arrays
      integer, allocatable :: iaxl(:)
      real(8),allocatable :: wk(:)
C ... Local parameters
      integer modep(3),nbaspp,mxnbr,i,isw
      integer, parameter :: niax=10
      procedure(integer) paird

C ... Set up input for call to pairc
      if (mxcsiz > 0) then
        mxnbr = mxcsiz
      else
        mxnbr = paird(nbas,nbasp,alat,plat,rmax,mxcsiz)
      endif
      do  i = 1, 3
        modep(i) = 2
        if (i > nd) modep(i) = 0
      enddo
      nbaspp = 2*nbasp - nbas
      nttab = mxnbr
      isw = 0

C ... Allocate iax array and make the neighbor table
      allocate(iaxl(niax*mxnbr)); call iinit(iaxl,niax*mxnbr)
      allocate(wk(3*mxnbr))
      call pairc(1,nbasp,nbaspp,modep,isw,ipsp,alat,plat,baspp,baspp,
     .  rmax,iltab,pltab,nttab,ntab,iaxl,wk,mxcsiz)
      if (allocated(iax)) then
        if (size(iax) /= niax*nttab) then
          deallocate(iax)
          allocate(iax(niax*nttab))
        end if
      else
        allocate(iax(niax*nttab))
      end if
      iax = iaxl(1:niax*nttab)
      deallocate(wk,iaxl)
      end subroutine pairs

      subroutine pairss(s_ctrl,s_lat,s_str,lpbc,rmaxs)
C- Pair table, read/written in structure s_str
C ----------------------------------------------------------------------
Ci Inputs
Ci   lpbc  :0 for 3D crystals
Ci         :11 for pgf mode
Ci   rmax  :maximum range for connecting vector, in a.u.
Co Outputs
Cs Command-line switches
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   08 Jul 13 Replace f77 pointers with f90 ones
Cu   08 May 13 Eliminate s_array
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer lpbc
      double precision rmaxs
C ... For structures
!      include 'structures.h'
      type(str_ctrl):: s_ctrl
      type(str_lat):: s_lat
      type(str_str):: s_str
C ... Dynamically allocated local arrays
      integer, allocatable :: iax(:)
C ... Local parameters
      integer mxcsiz,nbas,nbasp,i,j,niax
!       double precision xx
      procedure(integer) :: nglob
      parameter (niax=10)

      nbas = nglob('nbas')
      nbasp = nglob('nbasp')
!       call ptr_str(s_str,1,'npr',nbasp+1,0,xx)
      if (associated(s_str%npr)) deallocate(s_str%npr)
      allocate(s_str%npr(nbasp+1))

      if (lpbc == 0) then
        i = 3
        j = -1
      elseif (lpbc == 1 .or. lpbc == 11) then
        i = 2
        j = 1
      else
        call rx('PAIRSS: not implemented for lpbc>1')
      endif

      mxcsiz = s_str%mxnbr
      call pairs(nbas,nbasp,s_lat%alat,s_lat%plat,[rmaxs/2],s_lat%pos,
     .  [-1],i,j,s_ctrl%pgfsl,s_str%nttab,s_str%npr,iax,mxcsiz)
!       call ptr_str(s_str,4+1,'iax',niax*s_str%nttab,0,s_iax%p)
      if (associated(s_str%iax)) deallocate(s_str%iax)
      allocate(s_str%iax(niax*s_str%nttab))
      s_str%iax(:) = iax(1:niax*s_str%nttab)
      deallocate(iax)
      end

      subroutine pairc(ib1,ib2,nbasp,mode,isw,ips,alat,plat,pos,ctr,
     .  range,iltab,pltabp,nttab,ntab,iax,rtab,mxcsiz)
C- Make a neighbor table (crystal version)
C ----------------------------------------------------------------
Ci Inputs:
Ci  ib1,ib2:range of sites for which to generate tables
Ci   nbasp :the size of the basis, plus possible extensions.
Ci          Usually nbasp=nbas, but will differ in special
Ci          cases, such as having padding sites to extend
Ci          to a semi-infinite geometry.
Ci   mode:  vector of length 3 governing how pos shortened (see shorps)
Ci   isw:   1's digit fixes how range is calculated.
Ci           0: vector length must be < range(i)+range(j)
Ci           1: include all connecting vecs w/ r < range(i)
Ci         10's digit sets what part of iax table is not calculated
Ci           1: do not calculate iax(6)
Ci              (may be needed when ctr and pos are different)
Ci           2: calculate only iax(1..5)
Ci   ips   :index to which species each site belongs, for padded basis;
Ci          identifies which rmax is associated with each site. NB:
Ci          ips(1) = -1 => rmax is a global scalar, independent of site
Ci          ips(1) =  0 => rmax is site-, not species-dependent.
Ci          In either of these cases, ips is not used.
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat (input)
Ci   pos   :site positions (doubly padded for planar geometry)
Ci   ctr   :ctr(1..3,ib) is the effective center of the cluster
Ci          associated with site ib for purposes of computing distance
Ci          pos(jb)-ctr(ib).  May point to the same address space as pos
Ci   range :maximum range for connecting vector, in a.u..
Ci          This quantity may be a scalar, a species-dependent
Ci          array, or a site-dependent array, depending on ips(1);
Ci          see description of ips.  Precisely what meaning range has
Ci          depends on mode and isw.
Ci   iltab :iltab<0, has no effect.  Otherwise, see pltabp.
Ci         :Its purpose is explained in Remarks.
Ci   pltabp:include only pairs for which pltabp(jb)-pltabp(ib) <= iltab
Ci         :Only used if iltab>=0.
Ci   nttab :maximum dimension of iax table; used to guard against
Ci          generated table size exceeding dimension of iax.
Co Outputs:
Co   nttab    :total number of pairs generated
Co   iax      :neighbor table containing information about each pair ip
Co            :For each pair ip, information is contained in iax(*,ip).
Co            :as described below.  iax is ordered grouped by the basis
Co            :atoms, so that all pairs connected to site ib are grouped
Co            :together.  For each pair ip, iax(*,ip) contains:
Co   iax(1)   :site index to basis atoms ib=source;
Co            :all pairs with common ib are contiguous
Co   iax(2)   :site index to jb=field of each pair
Co   iax(3..5):multiples of plat added the difference in site positions
Co            :that connect the pair.
Co   iax(6)   :index to conjugate (jb,ib) pair matching (ib,jb)
Co            :NB: no matching pairs outside (ib1..ib2) can be found.
Co   iax(7)   :permutation index that reorders cluster.
Co            :How the cluster is reordered is application specific, but
Co            :most typically it is by an index that orders sites by
Co            :by increasing (x,y,z).  Particularly useful for mapping
Co            :two pair tables; see e.g. clcoin.f.
Co            :ppair4.f makes this kind of ordering.
Co   iax(8)   :iax(8) is left untouched by pairc;
Co            :Eleswhere it is used in two contexts:
Co            :1) if nonzero, points to an pair which is to substitute for pair i.
Co            :   Used in this context to substitute duplicate r.s. strux
Co            :   See strrs, strrs, streqv
Co            :2) as an internal flag for verious things; see e.g. hft2rs.f
Co   iax(9)   :iax(9) = (lmx+1)**2 for basis atom iax(2,i)
Co            :iax(9) is left untouched by pairc; set by mkiaxd
Co   iax(10)  :effective site index; see siteid.f
Co            :iax(10) is left untouched by pairc; set by siteid
Co   ntab     :ntab(ib)=number of pairs in iax table preceding ib
Co            :ntab is created for ib1:ib2+1.
Co   rtab     :rtab(1..3,ip) = pos(jb)-ctr(ib) for pair ip
Co   mxcsiz   :the largest cluster encountered
Cr Remarks
Cr   For each site ib=ib1..ib2, pairc finds all connecting vectors
Cr   for a lattice of points with periodic boundary conditions in
Cr   1, 2, or 3 dimensions, within a specified range of site ib.
Cr   The range can be defined in various ways, depending on isw.
Cr
Ci   iltab and pltab provide a means to add a constraint in generating a neighbor table.
Ci   It is used when making structure constant for lmpg:
Ci   iltab=1 and pltab holds the layer index.  Then neighbors outside
Cr   nearest PL are excluded from the table.
Cu Updates
Cu   23 Apr 02 added option to make only iax(1..5) (isw=20)
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer, parameter :: niax=10
      integer ib1,ib2,nbasp,mode(3),isw,nttab,ips(nbasp),
     .  ntab(ib1:ib2+1),iltab,pltabp(nbasp),mxcsiz
      integer iax(niax,*)
      double precision alat,plat(3,3),pos(3,nbasp),ctr(3,ib2),range(*),rtab(3,*)
C ... Dynamically allocated local arrays
      integer,allocatable :: iwk1(:,:),iwk2(:)
      real(8),allocatable :: wk(:,:),lat(:),posl(:),ctrl(:)
C ... Local parameters
      integer ib,is,jb,mtab,i,moder,mode2(3),nlat,mxntab,nsite
      double precision r1,rr,qlat(3,3),p0(3),xx

C     call tcn('ppair1')

C --- Setup ---
      nsite = ib2-ib1+1
      mxntab = nttab
      moder = mod(isw,10)
      do  i = 1, 3
      mode2(i) = mode(i)
        if (mode2(i) == 1) mode2(i) = 0
      enddo
C ... Make r1 = 2*maximum range
      r1 = range(1)
      if (ips(1) >= 0) then
        do  ib = 1, nbasp
        is = ib
        if (ips(1) > 0) then
          is = ips(ib)
        endif
        r1 = max(r1,range(is))
        enddo
      endif
      if (moder == 0) r1 = 2*r1
      r1 = 2*r1
C ... List of lattice vectors to add to pos(ib)-pos(jb)
      call xlgen(plat,r1/alat,0d0,0,20,mode,i,xx)
      allocate(lat(3*i))
      call xlgen(plat,r1/alat,0d0,i,0,mode,nlat,lat)
C ... qlat = (plat^-1)^T so that qlat^T . plat = 1
*     call prmx('plat',plat,3,3,3)
*     call prmx('starting pos',pos,3,3,nbasp)
      call mkqlat(plat,qlat,rr)
C ... Save true pos in oposl
C     and ctr in ctrl in case same address space used for ctr
      allocate(posl(3*nbasp))
      call dpcopy(pos,posl,1,3*nbasp,1d0)
      allocate(ctrl(3*nsite))
      call dpcopy(ctr,ctrl,1,3*nsite,1d0)

C --- For each ib, find all pairs for which dr < range ---
      nttab = 1
      ntab(ib1) = 0
      mtab = 1
      mxcsiz = 0
      do  ib = ib1, ib2
        r1 = range(1)
        if (ips(1) >= 0) then
          is = ib
          if (ips(1) > 0) then
            is = ips(ib)
          endif
          r1 = range(is)
        endif

C   --- Shorten all pos relative to ctr(ib) ---
C   ... Make pos-ctr(ib)
        call dpcopy(posl,pos,1,3*nbasp,1d0)
        call dpcopy(ctrl,ctr,1,3*nsite,1d0)
        do  i = 1, 3
          p0(i) = ctr(i,ib)
          do  jb = 1, nbasp
            pos(i,jb) = pos(i,jb) - p0(i)
          enddo
        enddo
C   ... Shorten pos-ctr(ib)
        call shorps(nbasp,plat,mode2,pos,pos)
C   ... Undo shift -ctr(ib) to restore shortened pos to absolute pos
        forall(i=1:3,jb=1:nbasp) pos(i,jb) = pos(i,jb) + p0(i)

C   --- Find all sites in range of ctr ---
        call ppair2(nbasp,iltab,pltabp,moder,alat,qlat,pos,p0,range,ips,
     .    rtab,ib,r1,nlat,lat,posl,mxntab,nttab,iax)

C   --- Sort table by increasing length ---
        allocate(iwk1(niax,nttab-mtab),iwk2(nttab-mtab))
        allocate(wk(3,nttab-mtab))
        call ppair3(nttab-mtab,iax(1,mtab),rtab(1,mtab),iwk1,iwk2,wk)
        deallocate(iwk1,iwk2,wk)

C   --- Cleanup for this ib ---
        mtab = nttab
        ntab(ib+1) = nttab-1
        mxcsiz = max(mxcsiz,ntab(ib+1)-ntab(ib))
      enddo
      nttab = nttab-1
C     call awrit2('xx ntab %n:1i',' ',80,6,nbasp+1,ntab)

C --- Restore original pos,ctr ---
      call dpcopy(posl,pos,1,3*nbasp,1d0)
      call dpcopy(ctrl,ctr,1,3*nsite,1d0)
      deallocate(lat,posl,ctrl)

C --- Fill out iax table ---
      call ppair1(isw,ib1,ib2,nbasp,ips,alat,plat,pos,range,
     .  nttab,ntab,iax,mxcsiz)

      end
      integer function paird(nbas,nbasp,alat,plat,rmax,mxcsiz)
C- Estimate an upper bound to the size of the neighbor table
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis (input)
Ci   nbasp :size of padded basis (layer programs)
Ci          nbasp = nbas + nbas(left bulk) + nbas(right bulk)
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat (input)
Ci   rmax  :maximum range for connecting vector, in a.u.
Ci          All connecting vectors with length < rmax()+rmax(j)
Ci          are retained.  rmax may be a scalar, a species-dependent
Ci          array, or a site-dependent array, depending on ipsp(1);
Ci          see description of ipsp
Ci  mxcsiz :if nonzero, use in place of internal formula for mxnbr
Co Outputs
Co  paird  :estimate for upper bound to total number of pairs in pair table
Co       : array; see pairc below where it is generated
Co   w() : array; see pairc, below where it is generated
Co   mxcsiz  :size of the largest cluster encountered
C ----------------------------------------------------------------------
      implicit none
      integer nbas,nbasp,mxcsiz
      double precision alat,plat(9),rmax(*)
C Local variables
      double precision avw,avwsr,vol
      integer nbaspp

      if (mxcsiz > 0) then
        paird = mxcsiz*nbasp
        return
      endif

C ... Estimate an upper bound to the size of the neighbor table
      nbaspp = 2*nbasp - nbas
      avw = avwsr(plat,alat,vol,nbaspp)
      paird = 3*(2*rmax(1)/avw)**3*nbasp
      end
      subroutine ppair1(isw,ib1,ib2,nbasp,ips,alat,plat,pos,range,
     .  nttab,ntab,iax,mxcsiz)
C- Fill out parts of the iax table
C ----------------------------------------------------------------
Ci  Inputs
Ci   isw   :1's digit fixes how range is calculated.
Ci           0: vector length must be < range(i)+range(j)
Ci           1: include all connecting vecs w/ r < range(i)
Ci         :10's digit sets what part of iax table is calculated
Ci           0: make iax(6),iax(7),iax(10)
Ci           1: make iax(7),iax(10)
Ci           2: no change to iax: printout only
Ci           4: just make iax(6)
Ci   ib1   :fill out iax table for pairs ntab(ib1)+1..ntab(ib2)
Ci   ib2   :fill out iax table for pairs ntab(ib1)+1..ntab(ib2)
Ci   nbasp :size of padded basis (not needed)
Ci   ips   :species table: site ib belongs to species ips(ib)
Ci         :Use to determine range (given by species).  But if:
Ci         :ips(1) = -1 => range is a global scalar, independent of site
Ci         :whose value is range(1)
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat
Ci   pos   :basis vectors
Ci   range :maximum range for connecting vector, in a.u..
Ci          This quantity may be a scalar, a species-dependent
Ci          array, or a site-dependent array, depending on ips(1);
Ci          see description of ips.  See 1s digit of isw for
Ci          how range is used.
Ci   nttab :total number of pairs in neighbor and iax (pairc.f)
Ci   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Ci   iax   :neighbor table containing pair information (pairc.f)
Ci   mxcsiz:maximum cluster size (for printout only)
Co  Outputs
Co   iax(6)   :index to conjugate (jb,ib) pair matching (ib,jb)
Co             NB: only matching pairs within site list can be found.
Co   iax(7)   :permutation index ordering cluster by increasing
Co             effective site index; see ppair4.f
Co   iax(10)  :effective site index; see siteid.f
C ----------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer isw,ib1,ib2,nbasp,nttab,niax,ips(nbasp),ntab(ib1:ib2+1)
      parameter (niax=10)
      integer iax(niax,*),mxcsiz
      double precision alat,plat(3,3),pos(3,*),range(*)
C ... Dynamically allocated arrays
      integer, allocatable :: iwk(:)
      real(8), allocatable :: posl(:)
C ... Local parameters
      integer ib,is,jb,js,ipr,i,j,moder,it,jt,stdo,nsite,isw1,mxcsizl
      double precision r1,r2,rr,rcut,vlat(3),drr2,tol
      parameter (tol=1d-5)
      procedure(integer) :: iprint,lgunit

      isw1 = mod(isw/10,10)
      ipr = iprint()
      moder = mod(isw,10)
      nsite = ib2-ib1+1
      stdo = lgunit(1)
      mxcsizl = 0

      if (isw1 == 2) goto 80

C --- Set iax(7) to sort this cluster ---
      if (isw1 /= 4) then
        call ppair5(ib1,ib2,plat,pos,tol,ntab,iax)
      endif

C --- For each pair, find matching pair, store in iax(6) ---
      do  it = 1, nttab
        iax(6,it) = 0
      enddo
      if (mod(isw1,2) == 0) then
        do  ib = ib1, ib2
        do  it = ntab(ib)+1, ntab(ib+1)
          if (iax(6,it) /= 0) cycle
          jb = iax(2,it)
C   ...   No matching pair for padded sites
          if (jb < ib1 .or. jb > ib2) cycle
          do  jt = ntab(jb)+1, ntab(jb+1)
C           if (iax(2,jt) == ib) then
C           endif
            if (iax(2,jt) == ib .and.
     .        iax(3,it) == -iax(3,jt) .and.
     .        iax(4,it) == -iax(4,jt) .and.
     .        iax(5,it) == -iax(5,jt))  then
              iax(6,it) = jt
              iax(6,jt) = it
              goto 1
            endif
          enddo
          call fexit3(-1,1,' Exit -1 pairc: cannot find pair'//
     .      ' matching sites (%i,%i), pair %i',ib,jb,it-ntab(ib))
    1   enddo
        enddo
      endif
      if (isw1 == 4) return

C ... Assign a unique id for every different site in the cluster table
      allocate(iwk(nttab))
      allocate(posl(3*nttab))
      call siteid(iax,nsite,ntab,plat,pos,posl,iwk,i)
      deallocate(iwk,posl)

! --- mode 20: only make iax(1:5) (do not write to the rest, a crash is induced) DMT2013
   80 continue
C --- Printout ---
      if (ipr < 30) goto 91
      if (ipr <= 40) write(stdo,'(1x)')
      if (ipr > 40) then
         write(stdo,"(/,'  ib  jb',9x,'--- r(jb)-r(ib) ---',10x,
     &                              'd       -x-plat-')", advance='no')
         if (isw1 /= 2) write(stdo,"('  map ord  id')", advance='no')
         write(stdo,'("")')
      end if
      i = 0
      do  it = 1, nttab
        ib = iax(1,it)
        jb = iax(2,it)
        rr = dsqrt(drr2(plat,pos(1,ib),pos(1,jb),
     .    iax(3,it),iax(4,it),iax(5,it),vlat))
        r1 = range(1)
        r2 = range(1)
        if (ips(1) >= 0) then
          is = ib
          if (ips(1) > 0) then
            is = ips(ib)
          endif
          r1 = range(is)
          js = jb
          if (ips(1) > 0) then
            js = ips(jb)
          endif
          r2 = range(js)
        endif
        if (moder == 0) rcut = r1+r2
        if (moder == 1) rcut = r1
        if (ib /= i) then
          mxcsizl = max(mxcsizl,ntab(ib+1)-ntab(ib))
          if (alat /= 1) write(stdo,345) ib,ntab(ib+1)-ntab(ib),rcut/alat,rcut
          if (alat == 1) write(stdo,345) ib,ntab(ib+1)-ntab(ib),rcut
  345     format(' pairc, ib=',i3,':',i4,' neighbors in range',f7.3:,'*alat =',f7.2)
        endif
        i = ib
        if (ipr > 40) then
!        if (ib == 12 .or. ib == 972) then
           write(stdo,'(i4,i4,3f11.6,f9.4,3x,3i3,i5,2i4)',advance='no')
     .      iax(1,it),iax(2,it),(vlat(j),j=1,3), rr, (iax(j,it), j=3,5)
           if (isw1 /= 2) write(stdo, '(i5,2i4)', advance='no')
     .                    (iax(j,it), j=6,7), iax(10,it)
           write(stdo,'("")')
        end if
      enddo
   91 continue

      if (mxcsizl == 0) mxcsizl = mxcsiz
      call info2(20,0,0,
     .  ' pairc:  %i pairs total  %i is max cluster size',nttab,mxcsizl)
C      if (ipr >= 20) write(stdo,
C     .  '('' pairc:'',i7,'' pairs total'',i6,'' is max cluster size'')')
C     .  nttab, mxcsizl

C     call tcx('ppair1')

      end
      subroutine ppair2(nbas,iltab,pltabp,moder,alat,qlat,pos,ctr,range,
     .  ips,rtab,ib,r1,nlat,lat,trupos,mxntab,nttab,iax)
C- Kernel of pairc to find all sites in range of ctr
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   iltab :iltab<0, has no effect.  Otherwise, see pltabp.
Ci   pltabp:include only pairs for which pltabp(jb)-pltabp(ib) <= iltab
Ci         :Only used if iltab>=0
Ci   moder : 1's digit fixes how range is calculated.
Ci         : 0: vector length must be < range(i)+range(j)
Ci         : 1: include all connecting vecs w/ r < range(i)
Ci   alat  :ength scale of lattice and basis vectors, a.u.
Ci   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
Ci   pos   :basis vectors
Ci   ctr   :ctr(1..3,ib) is the effective center of the cluster
Ci         :associated with site ib for purposes of computing distance
Ci         :pos(jb)-ctr(ib).  May point to the same address space as pos
Ci   range :maximum range for connecting vector, in a.u..
Ci          This quantity may be a scalar, a species-dependent
Ci          array, or a site-dependent array, depending on ips(1);
Ci          see description of ips.  Precisely what meaning range has
Ci          depends on mode and isw.
Ci   ips   :species table: site ib belongs to species ips(ib)
Ci   ib    :site index for which neighbor table is generated
Ci   r1
Ci   nlat
Ci   lat
Ci   trupos
Ci   mxntab:dimensioning: cannot accumulate more pairs than mxntab
Cio Inputs/Outputs
Cio  nttab :total number of pairs in neighbor table
Cio        :On input, starting number; on output nttab is enlarged
Cio        :to include the current cluster
Cio  iax   :neighbor table containing pair information
Cio        :Table is enlarged to include the current cluster
Cio  rtab  :site positions corresponding to entries in the iax table,
Cio        :relative to ctr.  Enlarged to include the current cluster.
Cl Local variables
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nbas,ib,iltab,ips(nbas),pltabp(nbas),niax,nlat,moder,mxntab,nttab
      parameter (niax=10)
      integer iax(niax,*)
      double precision alat,ctr(3),pos(3,nbas),range(nbas),rtab(3,*)
      double precision qlat(3,3),trupos(3,nbas),lat(3,*),r1
C ... Local parameters
      integer i,ilat,jb,js
      double precision r2,rr,rcut,vlat(3),xx,rcutba,dpos(3)

      do  jb = 1, nbas
        if (iltab >= 0) then
          if (abs(pltabp(jb)-pltabp(ib)) > iltab) cycle
        endif
        r2 = range(1)
        if (ips(1) >= 0) then
          js = jb
          if (ips(1) > 0) then
            js = ips(jb)
          endif
          r2 = range(js)
        endif
        if (moder == 0) rcut = r1+r2
        if (moder == 1) rcut = r1
        rcutba = (rcut / alat)**2
        dpos(1) = pos(1,jb)-ctr(1)
        dpos(2) = pos(2,jb)-ctr(2)
        dpos(3) = pos(3,jb)-ctr(3)

C   --- For each (ib,jb,ilat), do ---
        do  ilat = 1, nlat

        if (nttab > mxntab) call rxi(
     .      'pairc: table exceeds input maximum size,',mxntab)

C ...   Add to list if connecting vector within range
        rtab(1,nttab) = dpos(1) + lat(1,ilat)
        rtab(2,nttab) = dpos(2) + lat(2,ilat)
        rtab(3,nttab) = dpos(3) + lat(3,ilat)
        rr = rtab(1,nttab)**2+rtab(2,nttab)**2+rtab(3,nttab)**2

*        call awrit5('try ib,jb,ilat= %i %i %i rr=%;4d: %l',' ',80,
*     .    6,ib,jb,ilat,rr,rr < rcut)

C   --- Add to iax table if this pair in range ---
        if (rr < rcutba) then

C     ... vlat += shortening vector
          do  i = 1, 3
            rtab(i,nttab) = alat*rtab(i,nttab)
C           rtab(i,nttab) = alat*(rtab(i,nttab)+ctr(i)-pos(i,ib))
            vlat(i) = lat(i,ilat) + pos(i,jb) - trupos(i,jb)
          enddo

C     ... iax table for this pair
          iax(1,nttab) = ib
          iax(2,nttab) = jb
          do  i = 1, 3
            xx = vlat(1)*qlat(1,i)+vlat(2)*qlat(2,i)+vlat(3)*qlat(3,i)
            iax(2+i,nttab) = nint(xx)
          enddo
          nttab = nttab+1

        endif

        enddo
      enddo
      end
      subroutine ppair3(nttab,iax,rtab,iwk,iwk2,rwk)
C- Sort neighbor table by distance
      implicit none
      integer nttab,niax,iwk2(nttab),i,j,k
      parameter (niax=10)
      integer iax(niax,nttab),iwk(niax,nttab)
      double precision rtab(3,nttab),rwk(3,nttab)

      do  i = 1, nttab
        rwk(1,i) = rtab(1,i)
        rwk(2,i) = rtab(2,i)
        rwk(3,i) = rtab(3,i)
        do  k = 1, niax
          iwk(k,i) = iax(k,i)
        enddo
      enddo
      call dvshel(3,nttab,rtab,iwk2,11)
      do  i = 1, nttab
        j = iwk2(i)+1
        rtab(1,i) = rwk(1,j)
        rtab(2,i) = rwk(2,j)
        rtab(3,i) = rwk(3,j)
        do  k = 1, niax
          iax(k,i) = iwk(k,j)
        enddo
      enddo
      end
      subroutine ppair4(iclus,nclus,plat,pos,ctr,iwk,rtab,tol,iax)
C- Sort cluster by increasing (x,y,z) relative to its center
C ----------------------------------------------------------------
Ci Inputs
Ci   iclus,nclus: sort iax(iclus..nclus)
Ci   plat :primitive lattice vectors
Ci    pos :basis vectors
Ci    ctr :cluster origin:does not affect the ordering, but shifts rtab
Ci    iwk :integer work array of length nclus-iclus+1
Ci    tol :tolerance to which positions are considered coincident
Ci         tol<0 => sort iax by iax(1..5)
Co Outputs
Co   iax(7,iclus..nclus) orders the cluster by increasing (x,y,z)
Co         (or increasing iax(1..5) if tol < 0
Co   rtab  :connecting vectors rtab(1..3,ip) = pos(jb)-ctr
Co          for pair ip and jb=iax(2,ip)
Cr Remarks
Cr  Each cluster is sorted by increasing (x,y,z),
Cr  sorted by x first, then by y, then by z, thus guaranteeing that
Cr  all sites common to any pair of clusters are ordered the same.
C ----------------------------------------------------------------
      implicit none
      integer iclus,nclus,niax,iwk(nclus)
      parameter (niax=10)
      integer iax(niax,*)
      double precision plat(3,3),pos(3,1),ctr(3),rtab(3,nclus),tol
      integer ic,jb,ic0,ix,ia2,i,j,k
C Local variables
      double precision dx
C     integer jx
C     double precision wk2(3,nclus*3)
      dx(ia2,i,j,k) = pos(ix,ia2) +
     .                plat(ix,1)*i + plat(ix,2)*j + plat(ix,3)*k

      if (tol < 0) then
        call ivheap(niax,nclus-iclus+1,iax(1,iclus),iwk,1)
      else
        ic0 = 0
        do  ic = iclus, nclus
          jb = iax(2,ic)
          ic0 = ic0+1
          do  ix = 1, 3
            rtab(ix,ic0) = dx(jb,iax(3,ic),iax(4,ic),iax(5,ic)) - ctr(ix)
          enddo
        enddo
        call dvheap(3,nclus-iclus+1,rtab,iwk,tol,1)
      endif

      do  ic = iclus, nclus
        iax(7,ic) = iwk(ic-iclus+1)
      enddo

C ... Debugging ...
C      call ivprm(niax,nclus-iclus+1,iax(1,iclus),wk2,iwk,0)
C      call yprm('iax',0,wk2,0,niax,niax,nclus-iclus+1)
C      call awrit2('iwk %n:1i',' ',180,6,nclus-iclus+1,iwk)
C      do  30  ic = iclus, nclus
C      ic0 = ic-iclus+1
C   30 iwk(ic0) = iwk(ic0)-1
C      call dvperm(3,nclus-iclus+1,rtab,wk2,iwk,.true.)
C      do  32  ic = iclus, nclus
C        ic0 = ic-iclus+1
C        print 346, ic,(rtab(jx,ic0)+pos(jx,iax(1,ic)), jx=1,3)
C  346   format(i4,3f11.6)
C   32 continue

C     Debugging printout: show iax(7) sorts cluster according to increasing ib
C     if (tol < 0) then
C       print *, iax(2,iclus-1+iax(7,iclus:nclus))
C     endif
      end
      subroutine ppair5(ib1,ib2,plat,pos,tol,ntab,iax)
C- Sort a range of clusters according to tol
C ----------------------------------------------------------------------
Ci Inputs
Ci  ib1,ib2:range of clusters to sort
Ci   plat  :primitive lattice vectors, in units of alat
Ci   pos   :basis vectors
Ci   tol   :tolerance; see ppair4
Ci   ntab  :ntab(ib)=offset to neighbor table for cluster ib (pairc.f)
Ci   iax   :neighbor table containing pair information (pairc.f)
Co Outputs
Co   iax   :iax(7) is set to order cluster; see ppair4.
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer niax,ib1,ib2
      parameter (niax=10)
      integer iax(niax,*),ntab(ib2+1)
      double precision plat(3,3),pos(3,*),tol
      integer ib,nttab
      integer, allocatable :: iwk1(:)
      real(8), allocatable :: wk2(:)

C --- Set iax(7) to sort this cluster ---
      do  ib = ib1, ib2
        nttab = ntab(ib+1)-ntab(ib)
        if (nttab <= 0) cycle
        allocate(iwk1(nttab),wk2(nttab*3))
        call ppair4(ntab(ib)+1,ntab(ib+1),plat,pos,pos(1,ib),iwk1,wk2,tol,iax)
        deallocate(iwk1,wk2)
      enddo
      end

      subroutine pair3c(mode,kb,iax,ntab,ntabg,rtab,rham,n3c,iax3c,iax3x)
C- Indices to two pair table elements that link to a fixed third site
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit determines what is made
Ci         : 0  make n3c only
Ci         : 1  make n3c and iax3c
Ci         : 2  make n3c and iax3c and iax3x
Ci         :10s digit additional constraints on whether a triplet is accepted
Ci         : 0  no additional constraints
Ci         :    Constraints involve
Ci         :    di = |r(ib) - r(kb)|  dj = |r(jb) - r(kb)|
Ci         : 1  (di+dj)/2 > rham(1).  This is the only constraint implemented now
Ci         : 2  (di+dj+dij)/3 > rham(1).
Ci         :10s digit for sorting table
Ci         : 1  Sort iax3c by increasing triplet pair indices in iax table
Ci         : 2  Sort iax3c by increasing (ib,jb)
Ci   kb    :Site defining center of cluster
Ci   iax   :pair neighbor table (pairc.f).
Ci          Only iax(2) and iax(7) are used.
Ci          iax(2,*) = source or field site index for this table
Ci          iax(7,*) = list of permutations that orders sites
Ci                     in iax by increasing iax(2,:); see Remarks
Ci   rtab  :site positions corresponding to entries in iax table (mkrtab)
Ci   rham  :Used as constraint; see mode above
Ci   ntab  :ntab(ib)=number of pairs in iax preceding site ib.
Ci   ntabg :number of sites in cluster for which to find connecting vectors
Ci         :ntabg is normally ntab(kb+1)-ntab(kb) but it can be smaller
Ci         :if contribution to hamiltonian has shorter range to cluster center
Co Outputs
Co  n3c    :number of triplets connected to kb
Co  iax3c  :(1s digit mode > 0) one entry for triplet i3c.
Co         :iax3c(1,i3c) = iax(:,ip)  connecting vector between ib and kb
Co         :iax3c(2,i3c) = iax(:,jp)  connecting vector between jb and kb
Co         :iax3c(3,i3c) = iax(:,kp) connecting vector between ib and jb
Co  iax3x  :(1s digit mode > 1) large table, with extra information
Co         :iax3x(1,i3c) = kb
Co         :iax3x(2,i3c) = ib
Co         :iax3x(3,i3c) = jb
Co         :iax3x(4,i3c) = ip : pair index iax(:,ip) connecting ib and kb
Co         :iax3x(5,i3c) = jp : pair index iax(:,jp) connecting jb and kb
Co         :iax3x(6,i3c) = kp : pair index iax(:,kp) connecting ib and jb
Co         :The remaining elements iax3x(7:niax) are not used.
Cr Remarks
Cr   pair3c generates the information needed to connect three centers
Cr   in a crystal.  For each pair of connecting vectors in iax,
Cr   a third vector is sought in iax that connects the two.
Cr
Cr   A schematic depiction of association of variables with their pair labels
Cr
Cr            . origin
Cr
Cr    .  jp, index to pair table connecting field point to augmentation site
Cr    \
Cr     \ <- connecting vector.  iax3c(3) points to the entry in iax defining this vector
Cr      \
Cr       \.  ip, index to pair table connecting source point to augmentation site
Cr
Cr   Three-center terms connected to a fixed augmentation site (specified by
Cr   offset offi in the iax table) that contribute a 3-center term to the
Cr   hamiltonian for each pair (ic,jc).  Here ic and jc are elements in the iax table
Cr   each identifying one pair connected to the augmentation site.
Cr   is=ipsite(ic,jc) identifies element iax(:,is) that connects the (ic,jc) pair.
Cr
Cr   For efficient implementation, this routine relies on iax(7,*) being sorted.
Cr   ip(i)=iax(7,*) is a table of permutations ip that orders each local cluster so
Cr   that iax(2,ip(i)) increases monotonically with i (see ppair5 called with tol<0).
Cr   This guarantees that ib(ip(i+1)) > ib(ip(i))
Cu Updates
Cu   29 Jun 15 Rewritten
Cu   30 Jul 09 (S. Lozovoi) ip3iax redimensioned, input modified.
Cu             Note: original pair3c was formulated for the O(N) cluster method and is no longer used.
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer, parameter :: niax=10
      integer mode,kb,ntabg,ntab(*),n3c,iax3c(3,*),iax3x(niax,*),iax(niax,*)
      double precision rtab(3,*),rham(*)
      real(8), parameter :: tol=1d-6
C ... Dynamically allocated local arrays
      integer, allocatable :: iprm(:),iwk(:,:)
C ... Local parameters
      integer i3,j3,ip,jp,ib,jb,nc,ich,kp,ic0,offi,iprint,ntry,mod0,mod1,mod2,stdo,nglob
      integer iiax(5)
      double precision di,dj,dij,dlength

      mod0 = mod(mode,10)
      mod1 = mod(mode/10,10)
      mod2 = mod(mode/100,10)

C --- Double loop over pairs in table at kb (in permuted order) to find matching (ip,jp) ---
C     Permuted order => iax(2,ip) increases monotonically
      ntry = 0
      n3c = 0
      offi = ntab(kb)
      do  i3 = offi+1, ntab(kb+1)
        ip = iax(7,i3)           ! ip: pair index relative to start of kb cluster
        if (ip > ntabg) cycle ! ntabg may be shorter than ntab
        if (ip+offi > ntab(kb+1)) call rx('pair3c: iax improperly sorted')
        ib = iax(2,ip+offi)      ! Site index for source point ip
        ic0 = ntab(ib) ! offset to cluster centered at ib
C   ... ich = current offset in cluster centered at ib.  Because we loop in
C       sorted order, we need not reset ich for looping over pairs at jb below
        ich = ic0+1
        nc  = ntab(ib+1)
C   ... Assemble all field points jp connected to this ip
        do  j3 = offi+1, ntab(kb+1)
          jp = iax(7,j3)           ! jp: pair index relative to start of kb cluster
C         print *, i3,j3, ib, iax(2,jp+offi)  ! jb should increase monotonically
          if (jp > ntabg) cycle ! ntabg may be shorter than ntab
          if (jp+offi > ntab(kb+1)) call rx('pair3c: iax improperly sorted')
          ntry = ntry+1
C         Do upper triangle only
C         if (lx6 .and. jp < ip) goto 40

C         Vector connecting (ip,jp)
          iiax(3) = iax(3,jp+offi) - iax(3,ip+offi)
          iiax(4) = iax(4,jp+offi) - iax(4,ip+offi)
          iiax(5) = iax(5,jp+offi) - iax(5,ip+offi)

C     ... Find in cluster centered at ib connecting vector iiax, if it exists
          ich = ich-1
   42     continue
          ich = ich+1
          if (ich > nc) cycle  ! No more elements in cluster at ib
          kp = ic0+iax(7,ich)

C     ... Increment jp until iax(offi+jp) is at least as large as iax(iax(ipc))
C         and ich  until iax is at least as large as iax
          if (iax(2,offi+jp) > iax(2,kp)) goto 42
          if (iax(2,offi+jp) < iax(2,kp)) goto 40
          if (iiax(3) > iax(3,kp)) goto 42
          if (iiax(3) < iax(3,kp)) goto 40
          if (iiax(4) > iax(4,kp)) goto 42
          if (iiax(4) < iax(4,kp)) goto 40
          if (iiax(5) > iax(5,kp)) goto 42
          if (iiax(5) < iax(5,kp)) goto 40

          if (mod1 == 1) then
            di = dlength(3,rtab(1,ip+offi),1)
            dj = dlength(3,rtab(1,jp+offi),1)
C           dij = dlength(3,rtab(1:3,ip+offi)-rtab(1:3,jp+offi),1)
            if ((di+dj)/2 > rham(1)+tol) cycle
          elseif (mod1 == 2) then
            di = dlength(3,rtab(1,ip+offi),1)
            dj = dlength(3,rtab(1,jp+offi),1)
            dij = dlength(3,rtab(1:3,ip+offi)-rtab(1:3,jp+offi),1)
            if ((di+dj+dij)/3 > rham(1)+tol) cycle
          else
            call rx('pair3c: bad mode')
          endif

C     ... A match was found
          n3c = n3c+1
          if (mod0 == 1) then
            iax3c(1,n3c) = ip+offi
            iax3c(2,n3c) = jp+offi
            iax3c(3,n3c) = kp
          endif
C         if (mod0 > 1) ip3iax(jp,ip) = kp
   40     continue
        enddo                   ! j3
      enddo                     ! i3

C --- Construct lower triangle of table ---
C      if (lx6) then
C        do  50  i = 1, ntabg
C        do  50  j = i+1, ntabg
C          is = ip3iax(i,j)
C          if (is > 0) ip3iax(j,i) = iax(6,is)
C   50   continue
C      endif

C --- Printout ---
      call info5(31,0,0,' pair3c site%,4i : found %i triplets out of %i',kb,n3c,ntry,0,0)
      if (mod0 == 0) return

C ... Sort list by increasing ib
      if (mod2 >= 1) then
        allocate(iprm(n3c))
        if (mod2 == 1) then   ! Sort by increasing pair index
          call ivheap(3,n3c,iax3c,iprm,0)
        else                    ! Sort by increasing (ib,jb)
          allocate(iwk(2,n3c))
          do  i3 = 1, n3c
            ip = iax3c(1,i3)
            jp = iax3c(2,i3)
!           kp = iax3c(3,i3)
            ib = iax(2,ip)       ! Site index for source point ip
            jb = iax(2,jp)       ! Site index for source point jp
            iwk(1,i3) = ib
            iwk(2,i3) = jb
          enddo
          call ivheap(2,n3c,iwk,iprm,0)
          deallocate(iwk)
          allocate(iwk(3,n3c))
          call ivprm(3,n3c,iax3c,iwk,iprm,1)
          deallocate(iwk)
        endif
        deallocate(iprm)
      endif

C ... Make large table iax3x
      if (mod2 >= 2) then
        do  i3 = 1, n3c
          ip = iax3c(1,i3)
          jp = iax3c(2,i3)
          kp = iax3c(3,i3)
          ib = iax(2,ip)        ! Site index for source point ip
          jb = iax(2,jp)        ! Site index for source point jp
          iax3x(1,i3) = kb
          iax3x(2,i3) = ib
          iax3x(3,i3) = jb
          iax3x(4,i3) = ip      ! pair index iax(:,ip) connecting ib and kb
          iax3x(5,i3) = jp      ! pair index iax(:,jp) connecting jb and kb
          iax3x(6,i3) = kp      ! pair index iax(:,kp) connecting ib and jb
        enddo
      endif

C     Print out information in iax3c
      if (iprint() < 50) return
      stdo = nglob('stdo')
      call info0(10,0,0,'  iclus   ib    jb   indxi  indxj indxij    di%9fdj%9fdij%8f<d>')
      if (mod0 == 1) then
        do  i3 = 1, n3c
          ip = iax3c(1,i3)
          jp = iax3c(2,i3)
          kp = iax3c(3,i3)
          ib = iax(2,ip)    ! Site index for source point ip
          jb = iax(2,jp)    ! Site index for source point jp

          di = dlength(3,rtab(1,ip),1)
          dj = dlength(3,rtab(1,jp),1)
          dij = dlength(3,rtab(1:3,ip)-rtab(1:3,jp),1)

C         call info8(10,0,0,'%,6i%,6i%,6i%,7i%,7i%,7i%;11,6D%;11,6D', kb,ib,jb,ip,jp,kp,di,dj)
          write(stdo,'(3i6,3i7,3f11.6,f11.5)') kb,ib,jb,ip,jp,kp,di,dj,dij,(di+dj+dij)/3
        enddo
      return
      endif

C      do  i = 1, ntabg
C      do  j = 1, ntabg
C        ib = iax(2,j+offi)      ! Site index for source point ip
C        jb = iax(2,i+offi)      ! Site index for source point jp
C        ich = ip3iax(j,i)
C        if (ich == 0) then
CC          print 356, kb,ib,jb
CC  356     format(' pair3c, kb=',i3,':  missing 2C pair for ib,jb=',2i4)
C        elseif (iax(2,ich) /= ib .or. iax(1,ich) /= jb) then
C          print *, 'problem with i,j,is=',i,j,ib,jb,ich
C          call rx('bug in pair3c')
C        else
C          call info8(10,0,0,'%,6i%,6i%,6i%,7i%,7i%,7i',kb,ib,jb,i+offi,j+offi,ich,0,0)
C        endif
C      enddo
C      enddo

      end
