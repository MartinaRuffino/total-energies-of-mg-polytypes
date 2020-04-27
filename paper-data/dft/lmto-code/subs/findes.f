      subroutine findes(mode,nescut,mxesspec,s_ctrl,s_lat,alat,nbas,nclass,nl,
     .  nrclas,mxclas,nsymop,plat,symopm,symopv,ishores,lock,wsr,z)
C- Finds position and radius of possible empty spheres
C ----------------------------------------------------------------------
Cio Structures
Cio  s_ctrl :struct for program flow parameters; see structures.h
Ci     Elts read:  nesabc rmines rmaxes modep sclwsr omax1 omax2 wsrmax
Ci                 nspec clabl spid
Co     Stored:     ipc ips clabl
Co     Allocated:  ipc ips
Cio    Elts passed:ipc ips pgfsl nspec
Cio    Passed to:  *
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  avw
Co     Stored:     pos
Co     Allocated:  pos
Cio    Elts passed:pos plat0 qlat
Cio    Passed to:  *
Ci Inputs:
Ci   mode  : 1s digit
Ci         : 0 exit after finding ES and writing poses (and possibly essite) to disk
Ci         : 1 return without writing any files
Ci         :   s_ctrl%ipc, s_ctrl%ips are reallocated and new classes and species appended
Ci         :   s_ctrl%clabl is added to.  Note: caller should allocate s_ctrl%clabl to be at least mxclas
Ci         :   s_lat%pos is reallocated and new sites appended
Ci   alat  :length scale
Ci   nrxyz :no of divisions made along each lattice vector
Ci   nsymop:number of symmetry operations
Ci   plat  :primitive lattice vectors (scaled by alat)
Ci   symopm:symmetry operation matrix
Ci   symopv:symmetry operation vector (scaled by alat)
Cio Inputs/Outputs:
Cio  nbas  :number of atoms in the basis
Cio        :number of empty spheres is added
Cio  nclass:number of classes, atoms in same class are symmetry-related
Cio        :number of empty spheres classes is added
Cio  nescut:maximum number of species to add.  If 0, no constraint on number
Cio  nl    :number of l's
Cio  nrclas:number of atoms in the ith class
Cio  wsr   :Wigner-Seitz sphere radius (in atomic units)
Cio        :wsr should be dimensioned wsr(mxclas)
Cio        :On input the first nclass values must be set
Cio  z     :nuclear charge.
Cio        :z should be dimensioned z(mxclas)
Cio        :On input the first nclass values must be set
Cio  lock  :switches that lock sphere radii.
Cio        :lock should be dimensioned lock(mxclas)
Cio        :On input the first nclass values must be set
Cl Local variables
Cl   volfac:sum of sphere volumes = volfac * cell volume
Cl         : 2 or 3 Freeze size of ES once they are added
Cl         :        otherwise, all ES are resized each time a new ES is added
Cl   omax1:max. allowed overlap divided by distance (s1+s2-d)/d<omax1
Cl   omax2:max. allowed overlap divided by radius  (s1+s2-d)/s1<omax2
Cl   rmaxes:maximum radius for empty sphere (in atomic units)
Cl   rmines:minimum radius for empty sphere (in atomic units)
Cs   --mino    : minimize sphere overlaps
Cs   --nescut=#: stop adding sites once the number crosses threshold #
Cs   --shorten : Suppress all shortening of basis vectors
Cs   --wsite   : Write site file with ES positions
Cs   --wsitex  : Write site file with ES positions, as multiples of plat
Cr Remarks:
Cr   This series of subroutines was adapted from the Stuttgart package
Cr   to facilitate finding empty spheres on a lattice.
Cr
Cr   This subroutines searches possible positions for empty
Cr   (better: interstitial spheres). The algorithm is purely geometric,
Cr   the potential is not used at all.
Cr   Starting from a muffin-tin geometry, the program searches for
Cr   that position where one can put the largest empty sphere.
Cr
Cr   The symmetry of the crystal is preserved, so one needs on input
Cr   all symmetry-operations of the crystal.
Cr
Cr   In order to consider all atoms, the basis vectors must lie in
Cr   the WS-primitive cell :  -0.5 < b1,b2,b3 < 0.5  !!
Cr   with: bas(m) = b1*plat(m,k)+b2*plat(m,2)+b3*plat(m,3)
Cr   A call to ORDBAS garantees that this condition is fulfilled.
Cu Updates
Cu   17 Oct 17 New argument mxesspec
Cu   25 May 16 Added new mode
Cu   26 Jul 15 Revised and improved
Cu   17 Jun 13 Replace f77 pointers with f90 ones
Cu   08 May 13 Eliminate s_array
Cu   10 Nov 11 Begin migration to f90 structures
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nbas,nclass,nescut,nl,nsymop,mxclas,mode,mxesspec,ishores(3)
      integer nrclas(mxclas),lock(mxclas)
      double precision alat,plat(3,3),symopm(nsymop),symopv(nsymop),wsr(mxclas),z(mxclas)
C ... For structures
!      include 'structures.h'
      type(str_ctrl)::  s_ctrl
      type(str_lat)::   s_lat
C ... Dynamically allocated arrays
      integer, allocatable :: iwtmp(:),irr(:),lini(:),icls1(:),ips1(:),ntab(:)
      integer, allocatable :: iax(:)
      real(8), allocatable :: rirr(:),res(:),bas1(:),pos0(:,:)
      character*8,allocatable :: slabl(:)
C ... Local parameters
      logical cmdopt
      integer ic,is,idamax,ierr,jmax(3),jmin(3),lmx,nb1,nirr,nirrmx,nmesh,modep(3),opt1,nspec
      integer stdo,nrxyz(3),nclass0,nspec0,nbas0,i,it(3),nttab
      double precision avstep,bases(3),ratio,resmx,platcp(9),qlat(9),tiny,vol,volsph,wmax,wsres,
     .  wsrc(mxclas),volfac,rmaxes,rmines,omax1(3),omax2(3),wsmax,xx
      character*8 clabl(mxclas), strn*80
      integer,parameter :: NULLI=-99999
      parameter(tiny=1d-5)
      procedure(integer) :: nglob,a2vec,iprint,isw,iosite
C ... External calls
      external addes2,bigges,chcase,chkes,cross,daxpy,dcopy,defdr,
     .         deflmx,dgemm,dinv33,distes,dmpy,dpcopy,dpzero,dscal,
     .         dvheap,errmsg,fillat,getfmt,getirr,gtpmin,icopy,iinit,
     .         info0,info2,info5,latlim,mdesat,mdeses,nrmliz,ordbas,
     .         poppr,prpos,pshpr,ptr_ctrl,ptr_lat,r8tos8,reducv,renam,
     .         rlse,rsmesh,rx0,sclwsr,symes,tocast,zclabl

C --- Setup ---
      stdo = nglob('stdo')
      nrxyz = s_ctrl%nesabc
      rmines = s_ctrl%rmines
      rmaxes = s_ctrl%rmaxes
      modep = s_ctrl%modep
      volfac = s_ctrl%sclwsr
      opt1 = volfac/10
      volfac = volfac-10*opt1
      if (volfac == 0) volfac = 1
      omax1 = s_ctrl%omax1
      omax2 = s_ctrl%omax2
      wsmax = s_ctrl%wsrmax
      nbas0 = nbas
      nclass0 = nclass
      nspec0 = s_ctrl%nspec
      nspec = s_ctrl%nspec

      call info5(10,1,0,' --- FINDES : find new empty spheres'//
     .    '%?#n# up to %-1j%i new spheres##.  Constraints'//
     .    '%?#n# Rmt<=%d#%j#  rmines=%d  rmaxes=%d',
     .  nescut,isw(wsmax /= 0),wsmax,rmines,rmaxes)
C     call info2(30,0,1,'%5fConstraints:  rmines = %d  rmaxes = %d  ',rmines,rmaxes)
      do  ic = 1, nclass
        clabl(ic) = s_ctrl%clabl(ic)
      enddo

C --- Scale by alat ---
      call dscal(nclass,1.d0/alat,wsr   ,1)
      call dscal(1     ,1.d0/alat,rmaxes,1)
      call dscal(1     ,1.d0/alat,rmines,1)

C --- Use the most compact unit cell ---
      allocate(pos0(3,nbas0)); call dcopy(3*nbas0,s_lat%pos,1,pos0,1)
      call dcopy(9,plat(1,1),1,platcp,1)
      call dinv33(platcp,1,qlat,vol)
      vol = dabs(vol)*alat**3
      call ordbas(0,nbas,plat,s_lat%pos,s_ctrl%ipc,s_ctrl%ips)

C --- Get irreducible mesh points ---
      nmesh = nrxyz(1)*nrxyz(2)*nrxyz(3)
      nirrmx = 1000 + 4*nmesh/nsymop/3
C ... nirrmx should be nmesh but this requires to much space
      ierr = 0
   10 continue
      allocate(iwtmp(nirrmx)); call iinit(iwtmp,nirrmx)
      allocate(rirr(3*nirrmx))
      allocate(irr(nmesh)); call iinit(irr,nmesh)
      call rsmesh(avstep,ierr,irr,iwtmp,jmax,jmin,nrxyz(1),
     .            nrxyz(2),nrxyz(3),nirr,nirrmx,nsymop,platcp,qlat,
     .            rirr,symopm,symopv)
      if (ierr /= 0) then
        deallocate(irr,rirr,iwtmp)
        nirrmx = nirrmx + nirrmx
        goto 10
      endif
C     Can't do this ... need to keep values
C     deallocate(rirr); allocate(rirr(3*nirr))

c     write(stdo,301)avstep*alat

      allocate(res(nirr),lini(nirr)); call iinit(lini,nirr)

C --- Initialize res and lini ---
      call mdeses(alat,res,lini,nirr,nsymop,platcp,qlat,
     .            rirr,rmines-avstep,symopm,symopv)
      ratio = volsph(nclass,nrclas,wsr)*(alat**3)/vol

C --- Repeat until no more ES found, or packing fraction satisfied ---
      nb1 = 1
      do  while (.true.)

        if (iprint() >= 20) write(stdo,306) nclass-nclass0
        call info5(20,1,0,' ... Current vol frac = %;3,3d (target %;2,1d)  nbas = %i  nclass = %i',
     .    ratio,volfac,nbas,nclass,0)

C   --- Check whether the volume can be reached for given OMAX ---
C       wsrc is only used as a work array; except contents printed last iter
        call dpcopy(wsr,wsrc,1,nclass,alat)
        call sclwsr(30,nbas,nbas,nclass,alat,plat,s_lat%pos,s_ctrl%ipc,
     .    modep,clabl,z,lock,volfac,wsmax,omax1,omax2,wsrc)
        if (opt1 == 2 .or. opt1 == 3) call dpcopy(wsrc,wsr,1,nclass,1/alat)
C       call info0(30,1,0,' ')
C        print *, sngl(alat*wsr(1:nclass))
C        print *, sngl(wsrc(1:nclass))

        ratio = volsph(nclass,nrclas,wsrc)/vol
C       Exit, if volume target met
        if (dabs(ratio-volfac) < tiny) then
          call info2(20,0,0,' FINDES: reached target, sphere '//
     .      'vol = %d * cell vol',volfac,0)
          exit
        endif

C   ... Fill unit cell with potential new sites
        call fillat(alat,s_lat%pos,s_ctrl%ipc,jmax,jmin,nb1,nbas,nirr,
     .    nrxyz,lini,platcp,qlat,res,rirr,rmaxes,rmines-avstep,wsr)
        nb1 = nbas+1
C   ... Find biggest empty sphere in unit cell with center on mesh-point
        call bigges(alat,bases,iwtmp,nirr,lini,res,resmx,rirr)
C   ... Move empty sphere to high-symmetry point, leave the mesh
        if (resmx > tiny)
     .   call symes(alat,s_lat%pos,bases,clabl,s_ctrl%ipc,nbas,nrxyz,
     .               nsymop,platcp,qlat,symopm,symopv,wsr,wsres)
        if (resmx <= tiny .or. wsres < rmines) then
          call info2(20,1,0,' ... Largest empty sphere found '//
     .      'smaller than rmines = %d',rmines*alat,0)
          exit
        endif
        if (nclass+1 > mxclas) then
          call info2(20,1,0,' ... number of classes exceeds mxclas = %i',mxclas,0)
          exit
        endif
        wsres = dmin1(wsres,rmaxes)
        wmax = wsr(idamax(nclass,wsr,1))
C   ... Check possible overlap w/ atoms outside parallelipiped
        call chkes(alat,bases,platcp,qlat,wmax,wsres)
C   ... Increase number of classes by 1
        nclass = nclass+1
        ic     = nclass
        wsr(ic) = wsres
        z(ic)   = 0.d0
        call deflmx(lmx,wsr(ic),z(ic))
        nl = max0(nl,lmx+1)
        if (wsr(ic)*alat < 0.5d0) then
          call info2(20,0,0,'     (warning) new ES (class %i)'//
     .      ' radius is very small: r = %d',nclass,wsr(ic)*alat)
        endif
C   ... Find label for new empty sphere
        call renam(clabl,ic,nclass,z)
C       if (iprint() >= 20) write(stdo,303) clabl(ic),wsr(ic)*alat,nclass
        call info2(20,1,0,' ... add class %i ('//trim(clabl(ic))//
     .    ') with initial Rmax = %;5d',nclass,wsr(ic)*alat)

C   ... Add all symmetry-related spheres
        allocate(bas1(3*(nbas+nsymop)))
        allocate(icls1(nbas+nsymop))
        allocate(ips1(nbas+nsymop))
        call icopy(nbas,s_ctrl%ipc,1,icls1,1)
        call icopy(nbas,s_ctrl%ips,1,ips1,1)
        call dcopy(3*nbas,s_lat%pos,1,bas1,1)
        call addes2(bas1,bases,clabl,icls1,ips1,nbas,nspec,nclass,nrclas,
     .             nsymop,platcp,qlat,symopm,symopv)
        call ptr_ctrl(s_ctrl,4+2,'ipc',nbas,0,0,icls1)
        call ptr_ctrl(s_ctrl,4+2,'ips',nbas,0,0,ips1)
        call ptr_lat(s_lat,4+2,'pos',3,nbas,0,bas1)
        deallocate(bas1,icls1,ips1)
        if (nbas-nbas0>=nescut  .and. nescut>0) then
          wsrc(nclass) = wsres*alat
          call info2(20,1,1,' ... Exit search : number of addes ES = %i, exceeds limit (%i)',
     .      nbas-nbas0,nescut)
          call sclwsr(30,nbas,nbas,nclass,alat,plat,s_lat%pos,s_ctrl%ipc,
     .      modep,clabl,z,lock,volfac,wsmax,omax1,omax2,wsrc)
          call dpcopy(wsrc,wsr,1,nclass,1/alat)
          exit
        endif

      enddo
      deallocate(iwtmp,res,lini,irr,rirr)

C --- Cleanup ---
C ... Rescale wsr,rmines,rmaxes
      call dscal(nclass,alat,wsr   ,1)
      call dscal(1     ,alat,rmines,1)
      call dscal(1     ,alat,rmaxes,1)
      call ordbas(0,nbas,plat,s_lat%pos,s_ctrl%ipc,s_ctrl%ips)

      if (cmdopt('--mino',6,0,strn)) then
C       call info0(15,1,1,' FINDES : reposition ES to minimize overlaps ...')
        if (iprint() >= 16) write(stdo,307)
        i = 0
        allocate(ntab(nbas+1))
        call pshpr(iprint()-10)
        call pairs(nbas,nbas,alat,plat,[3*s_lat%avw],s_lat%pos,
     .    [-1],3,-1,s_ctrl%pgfsl,nttab,ntab,iax,i)
        call ovmin('~z ',nbas,nbas,alat,plat,wsrc,wsrc,
     .    clabl,s_ctrl%ipc,modep,z,ntab,iax,s_lat%pos,1)
        call poppr
      endif

C ... Printout
      if (iprint() >= 15) then
        call info2(15,1,1,' ... Final sphere sizes and packing fraction (%i new spheres)',nbas-nbas0,0)
        call sclwsr(30,nbas,nbas,nclass,alat,plat,s_lat%pos,s_ctrl%ipc,
     .    modep,clabl,z,lock,volfac,wsmax,omax1,omax2,wsr)
      endif

C ... Restore original site positions
      if (.not. cmdopt('--shorten',9,0,strn)) then
        call dcopy(3*nbas0,pos0,1,s_lat%pos,1)
      endif

      if (sum(ishores) /= 0) then
        it = ishores; it(1) = it(1)+60
        call shorps(nbas-nbas0,plat,it,s_lat%pos(1,nbas0+1),s_lat%pos(1,nbas0+1))
      endif

C ... Update s_ctrl%clabl(ic) and return
      if (mode == 1) then
        do  ic = nclass0+1, nclass
          s_ctrl%clabl(ic) = clabl(ic)
        enddo
        return
      endif

C ... Reduce the number of classes
      if (mxesspec /= NULLI .and. nspec-nspec0 > mxesspec) then
        do  i = nbas0+1, nbas
          if (s_ctrl%ips(i) > nspec0+mxesspec) then
            s_ctrl%ips(i) = nspec0+mxesspec
          endif
        enddo
        wsr(nspec0+mxesspec) = min(wsr(nspec0+mxesspec),wsr(nspec)) ! choose the smallest size
        nspec = nspec0 + mxesspec
      endif

C ... Unshear site positions and write basis to file poses
      call info0(20,1,0,' Writing new basis to file poses ...')
      call prpos(nbas,nclass,nclass0,clabl,s_ctrl%ipc,
     .  s_lat%plat0,s_lat%qlat,s_lat%pos,wsr,z)

      if (cmdopt('--wsitex',8,0,strn) .or. cmdopt('--wsite',7,0,strn)) then
        nspec = s_ctrl%nspec + nclass-nclass0
        allocate(slabl(nspec))
        do  is = 1, nspec
          if (is <= s_ctrl%nspec) then
            slabl(is) = s_ctrl%spid(is)
          else
            ic = nclass0 + is - s_ctrl%nspec
            slabl(is) = clabl(ic)
          endif
        enddo
        i = 15001; if (cmdopt('--wsitex',8,0,strn)) i = 15011
        ic = iosite(i,3d0,0,'essite',ic,slabl,alat,s_lat%plat0,nbas,
     .      s_ctrl%nspec,s_lat%pos,xx,xx,xx,s_ctrl%ips,xx,xx)
        deallocate(slabl)
      endif

      call rx0('findes')

  306 format(1x,/10('-'),' Seek new ES class having added',i3,2x,10('-'))
  307 format(1x,/10('-'),' Reposition ES to minimize overlaps',10('-'))
C  303 format(/' FINDES:  new class: ',a8,', Rmax=',f7.5,', class',i4)

      end

      subroutine symes(alat,bas,bases,clabl,iclass,nbas,nrxyz,nsymop,
     .                 plat,qlat,symopm,symopv,wsr,wsres)
C- Tries to move bases to a high-symmetry point
C ----------------------------------------------------------------------
Ci Inputs:
Ci   alat  :length scale
Ci   bas   :basis vectors (scaled by alat)
Cio  bases :basis vectors of new empty sphere (scaled by alat)
Ci   clabl :name of the different inequivalent atoms
Ci   iclass:the jth atom belongs to class iclass(j)
Ci   nbas  :number of atoms in the basis
Ci   nrxyz :no of divisions made along each primitive lattice vector
Ci   nsymop:number of symmetry operations as supplied by the generators
Ci   plat  :primitive lattice vectors (scaled by alat)
Ci   qlat  :primitive translation vectors in reciprocal space
Ci   symopm:symmetry operation matrix
Ci   symopv:symmetry operation vector (scaled by alat)
Ci   wsr   :Wigner-Seitz sphere radius (scaled by alat)
Co Outputs:
Cio  bases :basis vectors of new empty sphere (scaled by alat)
Co   wsres :Wigner-Seitz radius of new empty sphere
Cr Remarks
Cr   facmrg:decides when two empty spheres should be merged
Cr   facovl:decides how much empty spheres should overlap
Cr   First all symmetry-related positions of bases are generated.
Cr   Then the empty sphere is replaced by the average of all positions
Cr   within a sphere with radius wsres*facmrg centered around bases.
Cr   The parameter which decides when to merge the spheres: facmrg, is
Cr   set to 1.5 but this number may be changed.
C ----------------------------------------------------------------------
      implicit none
C Passed variables:
      integer iclass(*),nbas,nrxyz(3),nsymop
      double precision alat,bas(3,*),bases(3),plat(3,3),qlat(3,3),
     .                 symopm(3,3,*),symopv(3,*),wsr(*),wsres
      character*8 clabl(*)
C Local variables:
      integer ibas,ibasm,ii,isop,iprint,j(3),k(3),ll1,ltmax,lunit,
     .        m,mm,n,nrxyzz(3),nid,nit,nnz,nsbas
      double precision d,db,dbas(3),desesm,desatm,dm1,
     .                 dm1mx,dm2,d3nrm2,desat,deses,facmrg,facovl,fmwsr,
     .                 r(3),rb(3,3),sdbas(3),tiny,tolwsr
      logical ldum,lready
      parameter(ltmax=1,ll1=ltmax*2+1)
      parameter(facmrg=2.d0,facovl=1.d0,tiny=.1d-5,tolwsr=.1d-2)
C External calls:
      external  daxpy,dcopy,dpzero,distes,d3nrm2,getirr,icopy,iprint,
     .          lunit,mdesat,mdeses,popprt,pshprt
C Intrinsic functions:
      intrinsic  dmax1,dmin1,mod
C Statement functions:
      mm(ii,m)=ltmax-(mod(ii,ll1**m)-mod(ii,ll1**(m-1)))/ll1**(m-1)

      if (iprint() >= 60) write(*,300)facmrg,facovl

      call icopy(3,nrxyz,1,nrxyzz,1)

C --- Optimize position of empty sphere (especially low symmetry points)
      nit=0
      lready=.false.
      do while (.not.lready)
        nit=nit+1
        nnz=0
        db=1.d0
        do  m = 1, 3
          nrxyz(m)=nrxyz(m)*2
          do  n = 1, 3
            rb(n,m) = plat(n,m)/nrxyz(m)
C           qb(n,m) = qlat(n,m)*nrxyz(m)
          enddo
          db=db*d3nrm2(rb(1,m))
        enddo
        db=alat*db**(1.d0/6.d0)
        call getirr(k,nrxyz,qlat,bases)
        dm1mx=0.d0
        do  ii = 0, ll1**3-1
          do  m = 1, 3
            j(m) = k(m)+mm(ii,m)
            j(m) = mod(j(m),nrxyz(m))
            if (j(m) < 0) j(m)=j(m)+nrxyz(m)
            if (j(m) >= nrxyz(m)/2) j(m)=j(m)-nrxyz(m)
          enddo
          do  m = 1, 3
            r(m) = j(1)*rb(m,1) + j(2)*rb(m,2) + j(3)*rb(m,3)
          enddo
          call mdesat(alat,bas,desat,ibas,iclass,nbas,plat,r,wsr)
          call mdeses(alat,deses,ldum,1,nsymop,plat,qlat,r,0.d0,
     .                symopm,symopv)
          dm1=dmin1(desat,deses)
          if (dm1 > desat*0.5d0) nnz=nnz+1
          if (dm1 > dm1mx) then
            dm1mx=dm1
            desesm=deses
            desatm=desat
            ibasm=ibas
            call dcopy(3,r,1,bases,1)
          endif
          if (iprint() >= 100)
     .      write(*,311)ii,r,desat*alat,deses*alat
        enddo
        lready=db < tolwsr.or.nnz < 3
        if (nit == 1.or.lready) call pshpr(iprint()+10)
        if (iprint() >= 40) write(*,308)nit,bases,dm1mx*alat
        if (iprint() >= 50)
     .   write(*,305)desesm*alat,clabl(iclass(ibasm)),desatm*alat
        if (iprint() >= 60) write(*,307)db,tolwsr,nnz
        if (nit == 1.or.lready) call poppr
      enddo

C --- Consider symmetry-related positions
      fmwsr=facmrg*facovl*desatm+tiny+tiny
      nsbas=1
      nid=1
      call dpzero(sdbas,3)
      if (iprint() >= 50) write(*,301)
      do  isop = 2, nsymop
        call distes(bases,d,dbas,isop,plat,qlat,symopm,symopv)
        if (d < fmwsr) then
          call daxpy(3,1.d0,dbas,1,sdbas,1)
          nsbas=nsbas + 1
          if (d < tiny) nid = nid + 1
        endif
        if (iprint() >= 50)
     .    write(*,302)isop,dbas,d*alat,fmwsr*alat
      enddo
C --- Merge the ES
      if (nsbas/nid /= 1) then
        if (iprint() >= 30) write(*,303)bases,nid,nsymop/nid
        call dcopy(3,bases,1,bases,1)
        call daxpy(3,-1.d0/nsbas,sdbas,1,bases,1)
C ----- Now examine merged sphere, get deses,desat
        call mdesat(alat,bas,desat,ibas,iclass,nbas,plat,bases,wsr)
        call mdeses(alat,deses,ldum,1,nsymop,plat,qlat,bases,0.d0,
     .              symopm,symopv)
        dm2=dmin1(desat,deses)
        if (iprint() >= 30) write(*,304)bases,nsbas/nid
        if (iprint() >= 40)
     .    write(*,305)deses*alat,clabl(iclass(ibas)),desat*alat
        if (dm2 < dm1mx) then
          call daxpy(3,+1.d0/nsbas,sdbas,1,bases,1)
          if (iprint() >= 40) write(*,306)
        endif
      else
        dm2=0.d0
      endif

      wsres=dmax1(dm1mx,dm2)*facovl

      call icopy(3,nrxyzz,1,nrxyz,1)
      return

  300 format( ' SYMES : FACMRG=',f6.4,',  FACOVL=',f6.4)
  301 format(/' SYMES :ISOP',16x,'DBAS',16x,'DIST',5x,'FMWSR',/8x,57('-'))
  302 format(8x,i2,3x,3f10.6,1x,2f10.6)
  303 format( ' SYMES : BEFORE: ',3f10.6,',  degeneracy: ',i2,
     .       /9x,'Number of equivalent points in the unit cell:',i2)
  304 format( ' SYMES : AFTER:  ',3f10.6,', averaged over',i2,' points')
  305 format( ' SYMES : shortest distance to ES of same class :',f10.7,
     .      /'       : shortest distance to another atom ',a8,':',f10.7)
  306 format( ' SYMES : spheres not merged')
  307 format(' SYMES : DB=',f7.5,', TOLWSR=',f7.5,', NNZ=',i2)
  308 format( ' SYMES : STEP ',i2,':',3f10.6,'  WSR=',f8.6)
  311 format(2x,i3,3f11.6,2f10.7)
      end

      subroutine getirr(j,nrxyz,qlat,r)
C- gets coordinates j of irreducible mesh point
C-----------------------------------------------------------------------
Ci Inputs:
Ci   nrxyz :no of divisions made along each primitive lattice vector
Ci   r     :actual coordinate
Ci   qlat  :primitive translation vectors in reciprocal space
Co Outputs:
Co   j     :if r is a mesh point then rirr(m) = sum_k  rb(m,k)*j(k)
Co         :where rirr and r are connected by a lattice vector.
C-----------------------------------------------------------------------
      implicit none
C Passed variables:
      integer nrxyz(*),j(3)
      double precision qlat(3,*),r(3)
C Local variables:
      integer iprint,k,lunit
      double precision tiny,x(3)
      logical namp
      parameter(tiny=.1d-5)
C External calls:
      external iprint,lunit
C Intrinsic functions:
      intrinsic dabs,idnint,mod

      namp=.false.
      do  k = 1, 3
        x(k) = (r(1)*qlat(1,k)+r(2)*qlat(2,k)+r(3)*qlat(3,k))*nrxyz(k)
        j(k) = idnint(x(k))
        if (dabs(x(k)-j(k)) > tiny) namp=.true.
        j(k) = mod(j(k),nrxyz(k))
        if (j(k) < 0) j(k)=j(k)+nrxyz(k)
      enddo
      if (iprint() >= 70) write(*,300)r,j
      if (namp.and.iprint() >= 50)
     .  write(*,400)r,x,idnint(x(1)),idnint(x(2)),idnint(x(3))
ct     .  write(lunit(2),400)r,x,idnint(x(1)),idnint(x(2)),idnint(x(3))
300   format( ' GETIRR: VECTOR  ',3f10.6,'  -->',3i6)
400   format(/' GETIRR: the point:',3f8.4,' not a mesh point',
     .       /'         X=',3f8.2,
     .       /'         J=',3(i5,3x))
      end
      subroutine mdesat(alat,bas,desat,ibas,iclass,nbas,plat,r,wsr)
C- Finds the closest atom to a given mesh-point
C ----------------------------------------------------------------------
Ci Inputs:
Ci   bas   :basis vectors (scaled by alat)
Ci   ibas  :index of closest atom
Ci   iclass:the jth atom belongs to class iclass(j)
Ci   nbas  :number of atoms in the basis
Ci   plat  :primitive lattice vectors (scaled by alat)
Ci   wsr   :Wigner-Seitz sphere radius (scaled by alat)
Co Outputs:
Cio  desat  :distance from mesh-point to nearest sphere
Cr Remarks:
Cr Contrary to fillat the distance is even evaluated if it is smaller
Cr than rmines.
C ----------------------------------------------------------------------
      implicit none
C Passed variables:
      integer ibas,iclass(*),nbas
      double precision alat,bas(3,*),desat,plat(*),r(3),wsr(*)
C Local variables:
      integer i1,i2,i3,iprint,jbas
      double precision c01,c02,c03,c11,c12,c13,c21,c22,c23,
     .                 d1mach,d2,dx,dy,dz,r1,r2,r3,res,rpw2,wsri
C External calls:
      external d1mach,iprint
C Intrinsic functions:
      intrinsic dsqrt

      r1=r(1)
      r2=r(2)
      r3=r(3)
      desat=d1mach(2)
      do  jbas = 1, nbas
        rpw2=d1mach(2)
        wsri=wsr(iclass(jbas))
        c01=bas(1,jbas)-r1
        c02=bas(2,jbas)-r2
        c03=bas(3,jbas)-r3
        do  i1 = -1, 1
          c11=c01+i1*plat(1)
          c12=c02+i1*plat(2)
          c13=c03+i1*plat(3)
          do  i2 = -1, 1
            c21=c11+i2*plat(4)
            c22=c12+i2*plat(5)
            c23=c13+i2*plat(6)
            do  i3 = -1, 1
              dx=c21+i3*plat(7)
              dy=c22+i3*plat(8)
              dz=c23+i3*plat(9)
              d2=dx*dx+dy*dy+dz*dz
              if (rpw2 > d2) then
                res=dsqrt(d2)-wsri
                rpw2=res+wsri
                rpw2=rpw2*rpw2
                if (res < desat) ibas=jbas
              endif
            enddo
          enddo
        enddo
        desat=dmin1(desat,res)
        if (iprint() >= 70) write(*,300) ibas,res*alat,desat*alat
      enddo
300   format(' MDESAT: IBAS=',i3,', RES=',f11.6,', DESAT=',f11.6)

      end

      subroutine distes(bases,d,dbas,isop,plat,qlat,symopm,symopv)
C- Calculates the distance between symmetry-related sites
C ----------------------------------------------------------------------
Ci Inputs:
Ci   bases :basis vectors of new empty sphere (scaled by alat)
Ci   isop  :number of symmetry opration
Ci   plat  :primitive lattice vectors (scaled by alat)
Ci   qlat  :primitive translation vectors in reciprocal space
Ci   symopm:symmetry operation matrix
Ci   symopv:symmetry operation vector (scaled by alat)
Co Outputs:
Ci   d     :distance between symmetry-related sites
Cr Remarks:
C ----------------------------------------------------------------------
      implicit none
C Passed variables:
      integer isop
      double precision bases(*),d,dbas(*),plat(*),qlat(*),symopm(9,*),
     .                 symopv(3,*)
C Local variables:
      integer m
      double precision bast(3),d3nrm2
C External calls:
      external  daxpyl,dmpy,d3nrm2,reducv
C Intrinsic functions:
      intrinsic  dsqrt

      call dmpy(symopm(1,isop),3,1,bases,3,1,bast,3,1,3,1,3)
      call daxpy(3,1.d0,symopv(1,isop),1,bast,1)
      do  m = 1, 3
        dbas(m)=bases(m)-bast(m)
      enddo
      call reducv(plat,qlat,dbas,2)
      d=dsqrt(d3nrm2(dbas))

      end

      subroutine chkes(alat,bases,plat,qlat,wmax,wsres)
C- Checks possible overlap of empty sphere and atoms outside
C- the parallelipiped.
C ----------------------------------------------------------------------
Ci Inputs:
Ci   alat  :length scale
Ci   bases :basis vectors of new empty sphere (scaled by alat)
Ci   plat  :primitive lattice vectors (scaled by alat)
Ci   qlat  :primitive reciprocal lattice vectors
Ci   wmax  :largest Wigner-Seitz radius (scaled by alat)
Ci   wsres :Wigner-Seitz radius of empty sphere (scaled by alat)
Cr Remarks:
Cr   the minimal distance of the empty sphere to the faces of the
Cr   parallelipiped is calculated.
Cr   If this distance is smaller the the sum of the empty sphere radius
Cr   and the largest radius it could be possible that there is overlap
Cr   of the empty sphere and another sphere outside the parallelipiped.
C ----------------------------------------------------------------------
      implicit none
C Passed variables:
      double precision alat,bases(3),plat(3,3),qlat(3,3),wmax,wsres
C Local variables:
      integer i,iprint,k,lunit,m
      double precision ddot,dmin,qn(3,3),vec(3)
      character*288 messg
C External calls:
      external  ddot,errmsg,iprint,lunit,nrmliz
C Intrinsic functions:
      intrinsic  dabs,dmin1

      dmin=1.d20
      call nrmliz(3,qlat,qn)
      do  k = 1, 3
        do  i = -1, 1,2
          do  m = 1, 3
            vec(m)= 1.5d0*i*plat(m,k)-bases(m)
          enddo
          dmin=dmin1(dmin,dabs(ddot(3,vec,1,qn(1,k),1)))
        enddo
      enddo
      if (wsres+wmax-dmin > 0.d0) then
        write(messg,400)dmin*alat,wsres*alat,wmax*alat,
     .                     (dmin-wsres-wmax)*alat
        call errmsg(messg,1)
      elseif (iprint() >= 50) then
        write(*,300)dmin*alat,wsres*alat,wmax*alat,
     .                     (dmin-wsres-wmax)*alat
      endif

300   format(/' CHKES : DMIN=',f8.5,', WSR=',f8.5,', WMAX=',f8.5,
     .       /'         DMIN-WSR-WMAX=',f8.5)
400   format( ' CHKES : DMIN=',f8.5,', WSR=',f8.5,', WMAX=',f8.5,
     .        '|        DMIN-WSR-WMAX=',f8.5,
     .        '|        possible overlap of empty sphere',
     .                ' with atom out of range',
     .        '|        use most compact unit cell.$')
      end

      subroutine deflmx(lmx,wsr,z)
C- Returns default value of lmx for given nuclear charge and wsr
C ----------------------------------------------------------------------
Ci Inputs:
Ci   wsr   :Wigner-Seitz sphere radius (in atomic units)
Ci   z     :nuclear charge
Co Outputs:
Co   lmx   :lmx(j) = maximum l for atom j
Cr Remarks:
Cr The minimal lmx depends on the radius and is:
Cr
Cr       lmx=1  for 0.0 to 1.5
Cr       lmx=2  for 1.5 to 3.5
Cr       lmx=3  for 3.5 to 5.5
Cr       lmx=4  for 5.5 to 7.5
Cr       ....
C ----------------------------------------------------------------------
      implicit none
C Passed variables:
      integer lmx
      double precision wsr,z
C Local variables:
      integer iz,h,rb
      parameter (h=1,rb=37)
      character*72 messg
C External calls:
      external errmsg
C Intrinsic functions:
      intrinsic idnint,int,max0

      lmx = int(0.5d0*wsr+1.25d0)

      iz = idnint(z)
      if (wsr < 0.d0.or.wsr > 10.d0) then
        write(messg,400)wsr
        call errmsg(messg,1)
      endif

      if (iz >= rb) then
        lmx = max0(3,lmx)
      elseif (iz >= h) then
        lmx = max0(2,lmx)
      endif

  400 format(' DEFLMX: bad Wigner Seitz radius WSR=',f10.6,'$')
      end

      subroutine renam(clabl,ic,nclass,z)
C- Makes a new label for a class
C ----------------------------------------------------------------------
Ci Inputs:
Ci   ic    :class label
Ci   nclass:number of classes, atoms in same class are symmetry-related
Ci   z     :nuclear charge
Cio Inputs/Outputs:
Cio  clabl :name of the different inequivalent atom
Cio         nclass clabl's are needed on input
Cio         clabl(ic) is changed on output
Cr Remarks:
Cr A new clabl is found for class ic.
Cr If class ic is the only one with that nuclear charge, then
Cr the label will be H, He, ....
Cr If there are more than one classes of the nclass classes with
Cr the same nuclear charge, the label of class ic will be XXn with
Cr XX= H, He, ...  and n is the smallest nonzero integer which is not
Cr yet used in a label of a class ic < jc.
Cw Warning: if RENAME is called for class ic, it should also be called
Cw for classes ic+1, ..., nclass (else it could happen that different
Cw classes have the same label).
C ----------------------------------------------------------------------
      implicit none
C Passed variables:
      integer ic,nclass
      double precision z(nclass)
      character*8  clabl(nclass)
C Local variables:
      integer iz

      iz = idnint(z(ic))
      call zslabl(-1,clabl(ic),iz)
      call uniqclabl(0,ic-1,1,[ic],1,clabl,clabl(ic))

C      OLD
C      iz=idnint(z(ic))
C      lonly=.true.
C      do  jc = 1, nclass
C        if (ic /= jc.and.idnint(z(jc)) == iz) lonly=.false.
C      enddo
C      call zclabl(-1,clabl(ic),iz)
C      len=2
C      if (clabl(ic)(2:2) == ' ') len=1
C      iat=1
C      exist=.true.
C      do while (exist)
C        if (iat /= 1.or..not.lonly) then
C          call getfmt(' ',iat,1,2,1,3,fmt)
C          write(clabl(ic)(len+1:),fmt=fmt)iat
C        endif
C        exist=.false.
C        do  jc = 1, ic-1
C          if (clabl(ic) == clabl(jc)) exist=.true.
C        enddo
C        iat=iat+1
C      enddo
      end

      subroutine zclabl(iopt,clabl,iz)
C- Returns name of compound
C ----------------------------------------------------------------------
Ci Inputs:
Ci   iopt  :-1  z    =input  ; clabl=output
Ci           1  clabl=input  ; z    =output
Ci Inputs/Outputs:
Cio  clabl :name of chemical formula
Cio  iz    :nuclear charge
C ----------------------------------------------------------------------
      implicit none
C Passed variables:
      integer iopt,iz
      character*(*) clabl
C Local variables:
      integer icha,lunit
      character*2 aclabl(0:100)
      character*8 ccopy
C External calls:
      external  chcase,errmsg,lunit
C Intrinsic functions:
      intrinsic  ichar
C Data statements:
      data aclabl/'E ',
     .            'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     .            'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',
     .            'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     .            'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',
     .            'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     .            'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',
     .            'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     .            'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',
     .            'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     .            'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm'/

      if (iopt == -1) then
        if (iz < 0 .or. iz > 100) then
ct          write(lunit(2),400)iz
          write(*,400)iz
          return
        endif
        clabl=aclabl(iz)
      elseif (iopt == 1) then
        ccopy = clabl
        call  chcase( 1,1,ccopy(1:1))
        call  chcase(-1,1,ccopy(2:2))
        icha=ichar(ccopy(2:2))
        if (icha < 97.or.icha > 122)ccopy(2:2)=' '
        do  iz = 0, 100
          if (aclabl(iz) == ccopy(1:2)) then
            return
          endif
        enddo
ct        write(lunit(2),401)clabl
          write(*,401)clabl
        iz=-1
      else
        call errmsg(' ZCLABL: bad iopt.$',5)
      endif

400   format(/' ZCLABL: bad Z=',i5,' set to -1')
401   format(/' ZCLABL: could not find atom ',a8,' in my list,',
     .       ' Z set to -1')
      end

      subroutine getfmt(chr,res,nelts,icast,nwmax,ndigmx,fmt)
C- Supplies the format for writing a token
C ----------------------------------------------------------------------
Ci Inputs:
Ci   chr   :array of character
Ci   res   :array of logical,integer,real or double precision
Ci   nelts :number of elements in the array.
Ci   icast :integer controlling the type of the token(s):
Ci         :0=,logical, 1=char, 2=int, 3=real, 4=double
Ci   nwmax :maximum repetition number
Ci   ndigmx:maximum possible number of digits.
Co Outputs:
Co   fmt   :format descriptor
Co          EX: (11f12.5) (3 a 6)
Cr Remarks:
C     For a character token, nrdig returns the number of nonblank
C     characters in an array of strings.
C ----------------------------------------------------------------------
      implicit none
C Passed variables:
      integer icast,ndigmx,nelts,nwmax
      integer res(*)
      character*1 chr(*)
      character*(*) fmt
C Local variables:
      integer chrpos,i,icount,idig,k,lenchr,n,nrdig,nw,resi
      double precision d1mach,dpres,zero0,zero,resd,tiny
      real resr
      logical resl
      parameter (tiny=1.d-6)
C External calls:
      external  chrpos,d1mach,tocast
C Intrinsic functions:
      intrinsic dabs,dlog10,dnint,iabs,max0,min0

      n=0
      nrdig = 0
      idig = 0
      zero0=0.1d0**ndigmx
      do  icount = 1, iabs(nelts)
        call tocast(res,res,res,res,resl,resi,resr,resd,1,icount,icast)
        i=1
        if (icast == 1) then
          lenchr=res(1)
          k=1+lenchr*(icount-1)
          i=chrpos(chr(k),' ',1,lenchr)-1
        elseif (icast == 2) then
          if (iabs(resi) /= 0)   i=alog10(.5+iabs(resi))+1
          if (resi < 0) i=i+1
        elseif (icast >= 3) then
          if (icast == 3) dpres=dble(resr)
          if (icast == 4) dpres=resd
          zero=zero0
          do  idig = 0, ndigmx-1
            if (dabs(dpres-dnint(dpres)) < zero) goto 10
            dpres=dpres*10.0d0
            zero=zero*10.0d0
          enddo
   10     continue
          if (icast == 3) dpres=dble(resr)
          if (icast == 4) dpres=resd
          if (idig /= 0) i = 0
          if (dabs(dpres) >= 1.d0-tiny) i=dlog10(tiny+dabs(dpres))+1
          if (dpres < -zero0) i=i+1
        endif
        nrdig=max0(nrdig,idig)
        n=max0(n,i)
      enddo
      nw=min0(nelts,nwmax)
      if (icast == 0) then
        write(fmt,300) nw,'l',1
      elseif(icast == 1)then
        write(fmt,300) nw,'a',n
      elseif(icast == 2)then
        write(fmt,300) nw,'i',n
      elseif(icast == 3) then
        write(fmt,301) nw,'e',n+nrdig+1,nrdig
      elseif(icast == 4) then
        write(fmt,301) nw,'f',n+nrdig+1,nrdig
      endif
300   format('(',i2,a1,i2,')')
301   format('(',i2,a1,i2,'.',i2,')')
      end

      subroutine tocast(l,i,r,d,l1,i1,r1,d1,m,n,icast)
C ----------------------------------------------------------------------
Ci Inputs:
Ci   l     :array of logical
Ci   i     :array of integer
Ci   r     :array of real
Ci   d     :array of double precision
Co Outputs:
Ci   l1    :array of logical
Ci   i1    :array of integer
Ci   r1    :array of real
Ci   d1    :array of double precision
Cr Remarks:
C     For a character token, nrdig returns the number of nonblank
C ----------------------------------------------------------------------
      implicit none
C Passed variables:
      integer n,m,icast
      real    r(*),r1(*)
      integer i(*),i1(*)
      logical l(*),l1(*)
      double precision d(*),d1(*)

      if (icast == 0) then
        l1(m)=l(n)
      elseif (icast == 2) then
        i1(m)=i(n)
      elseif (icast == 3) then
        r1(m)=r(n)
      elseif (icast == 4) then
        d1(m)=d(n)
      endif

      end

      subroutine addes2(bas,bases,clabl,iclass,ips,nbas,nspec,nclass,nrclas,
     .             nsymop,plat,qlat,symopm,symopv)
C-Adds the new empty sphere basis atom and its symmetry related atoms
C ----------------------------------------------------------------------
Ci Inputs:
Cio  bas   :basis vectors (scaled by alat)
Cio  bases :basis vectors of new empty sphere (scaled by alat)
Ci   clabl :name of the different inequivalent atom
Ci   iclass:the jth atom belongs to class iclass(j)
Cio  nbas  :number of atoms in the basis
Ci   nclass:number of classes, atoms in same class are symmetry-related
Ci   nsymop:number of symmetry operations as supplied by the generators
Ci   plat  :primitive lattice vectors (scaled by alat)
Ci   qlat  :primitive translation vectors in reciprocal space
Co   symopm:symmetry operation matrix
Co   symopv:symmetry operation vector (scaled by alat)
Co Outputs:
Cio  bas   :basis vectors (scaled by alat)
Cio         on output list has been completed by the new positions
Cio  nbas  :number of atoms in the basis
Co   nrclas:number of atoms in the ith class
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer iclass(*),ips(*),nbas,nspec,nclass,nrclas(*),nsymop
      double precision plat(*),qlat(*),bas(3,*),bases(3),symopm(3,3,*),symopv(3,*)
      character*8 clabl(*)
C ... Local parameters
      integer i,isop,iprint,lunit,nbasnw,ibas,jbas
      double precision bas1(3)
      logical latrel,ladd
      character*72 messg
C External calls:
      external  daxpy,dcopy,dmpy,errmsg,iprint,latrel,lunit,reducv

C     ifi = fopn('espos')

      nbasnw = nbas
      do  isop = 1, nsymop
        call dmpy(symopm(1,1,isop),3,1,bases,3,1,bas1,3,1,3,1,3)
        call daxpy(3,1.d0,symopv(1,isop),1,bas1,1)
        jbas = 0
        ladd = .true.
        do while (jbas < nbasnw .and. ladd)
          jbas = jbas+1
          if (latrel(qlat,bas1,bas(1,jbas))) then
            ladd = .false.
            if (iclass(jbas) /= nclass) then
              write(messg,400)jbas,iclass(jbas)
              call errmsg(messg,4)
            endif
          endif
        enddo
        if (ladd) then
          nbasnw = nbasnw+1
          iclass(nbasnw) = nclass
          call reducv(plat,qlat,bas1,2)
          call dcopy(3,bas1,1,bas(1,nbasnw),1)
          if (isop == 1) nspec = nspec+1
          ips(nbasnw) = nspec
        endif
      enddo

ct      if (iprint() >= 10) then
ct        write(lunit(1),300)nbas,nbasnw,nbasnw-nbas
ct        write(lunit(1),301)(clabl(iclass(ibas)),iclass(ibas),
ct     .                         (bas(i,ibas),i=1,3),ibas=nbas+1,nbasnw)
ct      endif

        write(*,300) nbas,nbasnw,nbasnw-nbas
        write(*,301) (clabl(iclass(ibas)),iclass(ibas),
     .                     (bas(i,ibas),i=1,3),ibas=nbas+1,nbasnw)
ct        write(ifi,301) (clabl(iclass(ibas)),iclass(ibas),
ct     .                      (bas(i,ibas),i=1,3),ibas=nbas+1,nbasnw)
      nrclas(nclass) = nbasnw-nbas
      nbas = nbasnw

300   format(/' ADDES : The basis has been extended from ',i3,
     .         ' to ',i3,' atoms.',
     .       /9x,'The ',i2,' new positions are: ')
301   format(9x,'ATOM=',a8,' ICLASS=',i3,1x,'POS=',3f9.5)

400   format(' ADDES : empty sphere sits on atom ',i3,', class',i3,'$')
      end

      logical function latrel(qlat,pos1,pos2)
C- Checks if two positions differ by a lattice vector
C ----------------------------------------------------------------------
Ci Inputs:
Ci   qlat  :primitive translation vectors in reciprocal space
Ci   pos1  :first position
Ci   pos2  :second position
Co Outputs:
Co   latrel:.true. if positions differ by a lattice vector
Cr Remarks:
C ----------------------------------------------------------------------
      implicit none
C Passed variables:
      double precision qlat(3,*),pos1(*),pos2(*)
C Local variables:
      integer i,m
      double precision diff(3),tolpos,vdiff
      parameter(tolpos=4.d-3)
C Intrinsic functions:
      intrinsic  dabs,dnint

      latrel=.false.
      do  i = 1, 3
        diff(i)=pos1(i)-pos2(i)
      enddo
      do  m = 1, 3
        vdiff=diff(1)*qlat(1,m)+diff(2)*qlat(2,m)+diff(3)*qlat(3,m)
        vdiff=dabs(vdiff-dnint(vdiff))
        if (vdiff > tolpos) return
      enddo
      latrel=.true.
      end

      subroutine prpos(nbas,nclass,nclass0,clabl,iclass,plat0,qlat,pos,lmxb,wsr,z)
C- Writes positions to file, possibly restoring shear
Cu   19 Mar 13 Restore shear; also print out positions in multiples of plat
      implicit none
      integer nbas,ibas,fopn,nclass,nclass0,iclass(nbas),lmxb(nclass0)
      double precision  plat0(3,3),qlat(3,3),pos(3,nbas),xpos(3),posi(3)
      double precision wsr(nclass),z(nclass)
      character*(8) clabl(nclass)
      integer ifi,lmx,icls

      ifi = fopn('poses')
      write(ifi,"(/'SITE')")
      do  ibas = 1, nbas
        call dgemm('T','N',3,1,3,1d0,qlat,3,pos(1,ibas),3,0d0,xpos,3)
        call dgemm('N','N',3,1,3,1d0,plat0,3,xpos,3,0d0,posi,3)
C       write(ifi,301) clabl(iclass(ibas)),(pos(i,ibas),i=1,3)
C       write(ifi,301) clabl(iclass(ibas)), posi, xpos
        call awrit2('%9fATOM='//clabl(iclass(ibas))//' POS=%3;12,7D  XPOS=%3;12,7D',' ',
     .    120,ifi,posi,xpos)
        pos(1:3,ibas) = posi(1:3)
      end do

      write(ifi,"(/'SPEC')")
      do  icls = 1, nclass
        call deflmx(lmx,wsr(icls),z(icls))
        if (icls <= nclass0) then
          write(ifi,301) clabl(icls),z(icls),wsr(icls),lmxb(icls)
        else
          write(ifi,302) clabl(icls),z(icls),wsr(icls),lmx
        endif
      end do

  301 format(2x,'ATOM=',a8,3x,'Z=',f4.1,'  R=',f9.6,'               LMX=',i1)
  302 format(2x,'ATOM=',a8,3x,'Z=',f4.1,'  R=',f9.6,'*{les==1?1:0}  LMX=',i1)
      end

      subroutine bigges(alat,bases,iwtmp,nirr,lini,res,resmx,rirr)
C- Finds biggest emtpy sphere in the irreducible part of the unit cell
C ----------------------------------------------------------------------
Ci Inputs :
Ci   alat  :length scale
Ci   iwtmp :weight of mesh-point
Ci   lini  :logical bitmap, T if point closer than rmines to a sphere
Ci   n1-3  :no of divisions made along each lattice vector
Ci   nirr  :number of irreducible mesh points
Ci   res  :square of radius of empty sphere (in units of mesh length)
Ci   rirr  :coordinates of irreducible mesh points
Co Outputs:
Co   bases :basis vectors of new empty sphere (scaled by alat)
Co   resmx:largest res, zero if all lini are .true.
Cr Remarks:
Cr   only those points with lini=.false. are considered
C ----------------------------------------------------------------------
      implicit none
C Passed variables:
      integer iwtmp(*),nirr
      double precision alat,bases(3),rirr(3,*),res(*),resmx
      logical lini(*)
C Local variables:
      integer iirr,iprint,iw,iwtmp0,k,lunit
      double precision r2
C External calls:
      external dcopyl,iprint,lunit

ct      if (iprint() >= 50) write(lunit(1),301)
      if (iprint() >= 50) write(*,301)
      resmx=0.d0
      iwtmp0=48
      do  iirr = 1, nirr
        if (.not.lini(iirr)) then
          iw=iwtmp(iirr)
          r2=res(iirr)
          if (r2 > resmx.or.(r2 == resmx.and.iw < iwtmp0)) then
            call dcopy(3,rirr(1,iirr),1,bases,1)
            resmx=r2
            iwtmp0=iw
          endif
          if (iprint() >= 60) then
            if (r2 >= resmx.or.iprint() >= 65)
     .       write(*,302)iirr,(rirr(k,iirr),k=1,3),r2*alat,iw
ct     .       write(lunit(1),302)iirr,(rirr(k,iirr),k=1,3),r2*alat,iw
          endif
        endif
      enddo

      if (resmx /= 0.d0) then
ct        if (iprint() >= 40) write(lunit(1),300)bases,resmx*alat,iwtmp0
        if (iprint() >= 40) write(*,300)bases,resmx*alat,iwtmp0
      else
ct        if (iprint() >= 40) write(lunit(1),303)resmx*alat
        if (iprint() >= 40) write(*,303)resmx*alat
      endif
      return

  300 format(/' BIGGES: BASES=  ',3f10.6,'  WSR=',f8.6,' IWT=',i2)
  301 format(/' BIGGES: IIRR',15x,'RIRR',16x,'RES   IWTMP',
     .       /8x,51('-'))
  302 format(7x,i6,3f10.6,2x,f7.5,1x,i4)
  303 format(/' BIGGES:  RESMX=',f8.6)
      end

      subroutine fillat(alat,bas,iclass,jmax,jmin,nb1,nbas,nirr,nrxyz,
     .                  lini,plat,qlat,res,rirr,rmaxes,rmines,wsr)
C- Fills the unit cell with atoms; finds the interstitial region
C ----------------------------------------------------------------------
Ci Inputs:
Ci   alat  :length scale
Ci   bas   :basis vectors (scaled by alat)
Ci   iclass:the jth atom belongs to class iclass(j)
Ci   irr   :points to irreducible mesh-point
Ci   jmax  :upper limit for coordinates of irreducible points
Ci   jmin  :lower limit for coordinates of irreducible points
Ci   nrxyz :no of divisions made along each lattice vector
Ci   nb1   :filling with atoms nb1,...,nbas
Ci   nbas  :number of atoms in the basis
Ci   nirr  :number of irreducible mesh points
Ci   plat  :primitive lattice vectors (scaled by alat)
Ci   qlat  :primitive translation vectors in reciprocal space
Ci   rirr  :coordinates of irreducible mesh points
Ci   rmaxes:maximum radius for empty sphere (scaled by alat)
Ci   rmines:minimum radius for empty sphere (scaled by alat)
Ci   wsr   :Wigner-Seitz sphere radius (scaled by alat)
Co Outputs:
Cio  lini  :logical bitmap, T if point closer than rmines to a sphere
Cio  res   :distance from mesh-point to nearest sphere
Cr Remarks:
Cr   Each irreducible mesh point it checked whether it lies
Cr   closer or not to a sphere than rmines.
Cr
Cr   In order to consider all atoms, the basis vectors must lie in
Cr   the WS-primitive cell i.e. with:
Cr
Cr          bas(m) = b1*plat(m,k)+b2*plat(m,2)+b3*plat(m,3)
Cr
Cr   b1,b2,b3 must lie between -0.5 and +0.5 !!!
Cr   Further the unit cell must be compact (be sure that plat was passed
Cr   trough CPPLAT).
Cr   If these two conditions are fullfilled, it is sufficient to sweep
Cr   the 3*3*3=27 unit cells around the central unit cell.
Cr
Cr   The irreducible points come from RSMESH have m-coorrdinates:
Cr
Cr          rirr(m) = g1*plat(m,k)+g2*plat(m,2)+g3*plat(m,3)
Cr
Cr   and gm fulfills:
Cr
Cr       -0.5 <= jmin(m)/nrxyz(m) <= gm <= jmax(m)/nrxyz(m) <= 0.5
Cr
Cr   Atoms with a distance to this paralleliped larger than rmaxes
Cr   are not considered (ins=false).
C ----------------------------------------------------------------------
      implicit none
C Passed variables:
      integer iclass(*),jmax(3),jmin(3),nrxyz(*),nb1,nbas,nirr
      double precision alat,bas(3,*),plat(3,*),qlat(3,*),res(*),
     .                 rirr(3,*),rmaxes,rmines,wsr(*)
      logical lini(*)
C Local variables:
      integer i(3),ibas,ic,ii,iirr,inz,iprint,isi,iz,ll1,ltmax,lunit,
     .        m,mm,stdo,nglob
      double precision a(3),c(3),c1,c2,c3,d3nrm2,d2,d12,dx,dy,dz,f(3),
     .                 pn(3,3),qn(3,3),rpw2,vec1(3),
     .                 vec3(3),wsri,wprmn2,wprmx2
      logical ins
      parameter(ltmax=1,ll1=ltmax*2+1)
C External calls:
      external  cross,d3nrm2,iprint,lunit,nrmliz
C Intrinsic functions:
      intrinsic  dsqrt,iabs,mod
C Statement functions:
      mm(ii,m)=ltmax-(mod(ii,ll1**m)-mod(ii,ll1**(m-1)))/ll1**(m-1)

      stdo=nglob('stdo')

      if (iprint() >= 70.and.iprint() < 80) write(stdo,300)
      if (iprint() >= 80) write(stdo,302)
      if (iprint() >= 70.and.iprint() < 80) write(stdo,300)
      if (iprint() >= 80) write(stdo,302)

C --- normalize plat and qlat, and store into pn and qn
      call nrmliz(3,plat,pn)
      call nrmliz(3,qlat,qn)

      do  ibas = nb1, nbas
c       call reducv(plat,qlat,bas(1,ibas),2)
        ic = iclass(ibas)
        wsri=wsr(ic)
        wprmx2=(wsri+rmaxes)*(wsri+rmaxes)
        wprmn2=(wsri+rmines)*(wsri+rmines)
        if (iprint() >= 70.and.iprint() < 80)
     .    write(stdo,301) ibas,ic,(bas(m,ibas),m=1,3),wsr(ic)*alat
ct        write(stdo,301) ibas,ic,(bas(m,ibas),m=1,3),wsr(ic)*alat
        do  ii = 0, ll1**3-1
          do  m = 1, 3
            i(m)=mm(ii,m)
            if (i(m) < 0) then
              f(m)=jmin(m)/dfloat(nrxyz(m))
            elseif (i(m) == 0) then
              f(m)=0.5d0*(jmin(m)+jmax(m))/dfloat(nrxyz(m))
            else
              f(m)=jmax(m)/dfloat(nrxyz(m))
            endif
          enddo
          do  m = 1, 3
           c(m)=bas(m,ibas)+i(1)*plat(m,1)+i(2)*plat(m,2)+i(3)*plat(m,3)
           a(m)=f(1)*plat(m,1)+f(2)*plat(m,2)+f(3)*plat(m,3)
           vec1(m)=c(m)-a(m)
          enddo

          isi=iabs(i(1))+iabs(i(2))+iabs(i(3))
          if (isi == 2) then
            do  m = 1, 3
              if (i(m) == 0) iz=m
            enddo
            call cross(vec1,pn(1,iz),vec3)
            d12=d3nrm2(vec3)
          elseif (isi == 1) then
            do  m = 1, 3
              if (iabs(i(m)) == 1)  inz=m
            enddo
            d12=vec1(1)*qn(1,inz)+vec1(2)*qn(2,inz)+vec1(3)*qn(3,inz)
            d12=d12*d12
          elseif (isi == 3) then
            d12=d3nrm2(vec1)
          elseif (isi == 0) then
            d12=0.d0
          endif
          ins=(d12 < wprmx2)
C ------- ins is false when atom ibas in unit-cell corresponding to
C ------- ii can never overlap with an empty sphere in the irreducible
C ------- volume
          if (ins) then
            c1=c(1)
            c2=c(2)
            c3=c(3)
            do  iirr = 1, nirr
              if (.not.lini(iirr)) then
                dx=c1-rirr(1,iirr)
                dy=c2-rirr(2,iirr)
                dz=c3-rirr(3,iirr)
                d2=dx*dx+dy*dy+dz*dz
                if (d2 < wprmn2) then
                  lini(iirr)=.true.
                  res(iirr)=0.d0
                else
                  rpw2=res(iirr)+wsri
                  rpw2=rpw2*rpw2
                  if (d2 < rpw2) res(iirr)=dsqrt(d2)-wsri
                endif
              endif
            enddo
          endif
          if (iprint() >= 80)
     .      write(stdo,303)ibas,ic,c,i,d12*alat,wprmx2*alat,ins
ct        write(stdo,303)ibas,ic,c,i,d12*alat,wprmx2*alat,ins
        enddo
      enddo

300   format(/,' FILLAT: fill unit cell:',
     .    /,9x,'IBAS',2x,'ICLASS',12x,'BAS',13x,'WSR',/,9x,50('-'))
301   format(8x,i3,4x,i3,3x,3f8.5,3x,f8.5)
302   format(/,' IBAS  ICLASS',12x,'BAS',9x,
     .           'I1  I2  I3    D12 WPRMX2 INS',/,9x,63('-'))
303   format(1x,i3,4x,i3,3x,3f7.4,3(2x,i2),2x,f6.3,1x,f5.3,2x,l1)

      end

      subroutine mdeses(alat,deses,lini,nirr,nsymop,plat,qlat,rirr,
     .                  rmines,symopm,symopv)
C- Initializes lini and deses, taking into account the symmetry.
C-----------------------------------------------------------------------
Ci Inputs:
Ci   alat  :length scale
Ci         :needed only for output
Ci   nirr  :number of irreducible points
Ci   nsymop:number of symmetry operations
Ci   plat  :primitive lattice vectors
Co   rirr  :irreducible point
Ci   rmines:minimum radius for empty sphere (scaled by alat)
Ci   symopm:symmetry operation matrix
Ci   symopv:symmetry operation vector (scaled by alat)
Co Outputs:
Co   deses :half of distance from mesh-point to nearest equivalent point
Co         :deses = distance empty sphere - empty sphere
Co         :zero if lini=.true.
Co   lini  :logical bitmap, T if point closer than rmines to a sphere
Co         (see remarks)
Cr Remarks:
Cr
Cr Due to the given symmetry, if an empty sphere is put into the
Cr crystal, empty spheres must be put on all symmetry-related positions.
Cr If the minimum radius of the empty spheres is given,
Cr and if the spheres should not overlap,
Cr then parts of the unit cell can be excluded as possible positions
Cr for empty spheres because an empty sphere would overlap with
Cr an symmetry-related sphere.
C-----------------------------------------------------------------------
      implicit none
C Passed variables:
      integer nirr,nsymop
      double precision alat,plat(9),qlat(9),deses(*),rirr(3,*),rmines,
     .                 symopv(3,*),symopm(9,*)
      logical lini(*)
C Local variables:
      integer iirr,isop,iprint,lunit,m,nexcl
      double precision d1,d2,d3,df2min,diff2,d3nrm2,frm2,
     .                 p2min,pmin(3),r1,r2,r3,mat(9,48),tiny,tiny2,
     .                 v1,v2,v3
      parameter(tiny=1.d-6)
C External calls:
      external d3nrm2,gtpmin,iprint,lunit
C Intrinsic functions:
      intrinsic dmin1,dnint,dsqrt

      nexcl=0
      tiny2=tiny*tiny
      frm2=4.d0*rmines*rmines

C --- get shortest lattice vector
      call gtpmin(plat,pmin)
      p2min=d3nrm2(pmin)

C --- mat_ij=symopm_ij - delta_ij
      call dcopy(9*nsymop,symopm,1,mat,1)
      do  isop = 1, nsymop
        call daxpy(3,-1.d0,1.d0,0,mat(1,isop),4)
      enddo

      do  iirr = 1, nirr
        df2min=p2min
        r1=rirr(1,iirr)
        r2=rirr(2,iirr)
        r3=rirr(3,iirr)
        do  isop = 2, nsymop
C ------- d = r - rrot
C ------- where r is transformed to rrot by symmetry operation number is
C ------- v = d in units of plat and always between -.5 and +.5
C ------- dsqrt(diff2) is then the distance (scaled by alat) from r
C ------- to rrot (minus eventually some lattice vectors)

          d1=mat(1,isop)*r1+mat(4,isop)*r2+mat(7,isop)*r3+symopv(1,isop)
          d2=mat(2,isop)*r1+mat(5,isop)*r2+mat(8,isop)*r3+symopv(2,isop)
          d3=mat(3,isop)*r1+mat(6,isop)*r2+mat(9,isop)*r3+symopv(3,isop)
          v1=d1*qlat(1)+d2*qlat(2)+d3*qlat(3)
          v2=d1*qlat(4)+d2*qlat(5)+d3*qlat(6)
          v3=d1*qlat(7)+d2*qlat(8)+d3*qlat(9)
          v1=v1-dnint(v1)
          v2=v2-dnint(v2)
          v3=v3-dnint(v3)
          d1=v1*plat(1)+v2*plat(4)+v3*plat(7)
          d2=v1*plat(2)+v2*plat(5)+v3*plat(8)
          d3=v1*plat(3)+v2*plat(6)+v3*plat(9)
          diff2=d1*d1+d2*d2+d3*d3
          if (diff2 > tiny2) then
            if (diff2 < frm2) then
              lini(iirr)=.true.
              deses(iirr)=0.d0
              nexcl=nexcl+1
              goto 10
            endif
            df2min=dmin1(df2min,diff2)
          endif
        enddo
        deses(iirr)=dsqrt(df2min)*0.5d0
10      continue
      enddo

C --- Printout
      if (nirr > 1) then
        call info5(30,0,0,' MDESES: %i  mesh points left from %i.  Shortest lattice vector %d',
     .    nirr-nexcl,nirr,alat*dsqrt(p2min),4,5)
        if (iprint() >= 60) then
c          write(lunit(1),302)
          write(*,302)
          do  iirr = 1, nirr
            if (lini(iirr)) then
c              write(lunit(1),303)iirr,(rirr(m,iirr),m=1,3),rmines*alat
              write(*,303)iirr,(rirr(m,iirr),m=1,3),rmines*alat
            else
c              write(lunit(1),304)iirr,(rirr(m,iirr),m=1,3),
              write(*,304)iirr,(rirr(m,iirr),m=1,3),
     .                           deses(iirr)*alat
            endif
          enddo
        endif
      endif

  300 format(/' MDESES:',i6,' mesh points left from',i7,'.  Shortest lattice vector',f10.5)
  302 format(/' MDESES:  NIRR',7x,'Rx',8x,'Ry',8x,'Rz',8x,'RES0',
     .       /10x,47('-'))
  303 format(7x,i5,2x,3f10.4,2x,'<',f9.6)
  304 format(7x,i5,2x,3f10.4,1x,f11.6)
      end

      subroutine gtpmin(plat,pmin)
C- gets shortest nozero lattice vector
C ----------------------------------------------------------------------
Ci Inputs:
Ci   plat  :primitive lattice vectors
Co Outputs:
Co   pmin  :shortest nozero lattice vector
C ----------------------------------------------------------------------
      implicit none
C Passed variables:
      double precision pmin(3),plat(9)
C Local variables:
      integer i,i1,i2,i3,j,k
      double precision d(3),d1mach,d2,d2min,d3nrm2,rmax
C External calls:
      external d1mach,d3nrm2,latlim
C Intrinsic functions:
      intrinsic dsqrt,dmin1

      d2min=d1mach(2)
      rmax=dsqrt(dmin1(d3nrm2(plat(1)),d3nrm2(plat(4)),d3nrm2(plat(7))))
      call latlim(plat,rmax,i1,i2,i3)
      do  i = -i1, i1
        do  j = -i2, i2
          do  k = -i3, i3
            d(1)=i*plat(1)+j*plat(4)+k*plat(7)
            d(2)=i*plat(2)+j*plat(5)+k*plat(8)
            d(3)=i*plat(3)+j*plat(6)+k*plat(9)
            d2 = d3nrm2(d)
            if (d2 < d2min.and.(i /= 0.or.j /= 0.or.k /= 0)) then
              d2min=d2
              call dcopy(3,d,1,pmin,1)
            endif
          enddo
        enddo
      enddo
      end

      double precision function d3nrm2(r)
C-  'norm'-square
C ----------------------------------------------------------------------
Ci Inputs:
Ci   r     :vector
Co Outputs:
Co   d3nrm2:norm-square
Cr Remarks:
C ----------------------------------------------------------------------
      implicit none
C Passed variables:
      double precision r(3)

      d3nrm2=r(1)*r(1)+r(2)*r(2)+r(3)*r(3)

      end

      subroutine ordbas(mode,nbas,plat,bas,iclass,ips)
C- Orders basis in ascending order of iclass
C ----------------------------------------------------------------------
Ci Inputs:
Ci   mode  :0 do not reorder basis
Ci         :1 reorder basis
Ci   nbas  :number of atoms in the basis
Ci   plat  :primitive lattice vectors
Co Inputs/Outputs:
Cio  bas   :basis vectors
Cio  iclass:the jth atom belongs to class iclass(j)
Cio  ips   :the jth atom belongs to species ips(j)
Cw Work arrays:
Cw   wk    :work array of length 5*nbas
Cr Remarks:
Cr   Orders bas so that the first positions correspond to iclass=1,
Cr   followed by those corresponding to iclass=2, ....
Cr   Within each class, sites are also ordered
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nbas,iclass(nbas),ips(nbas)
      double precision bas(3,nbas),plat(3,3)
C ... Local parameters
      integer ib,ibp
      double precision danrm2,qlat(3,3),vol,wk(6,nbas)
      integer iprm(nbas)

      call info0(40,1,0,' ORDBAS:  shorten and order basis')

      call dinv33(plat,1,qlat,vol)
      do  ib = 1, nbas
C       call shorbz(bas(1,ib),bas(1,ib),plat,qlat)
        call reducv(plat,qlat,bas(1,ib),2)
        wk(1,ib) = dble(ips(ib))
        wk(2,ib) = dble(iclass(ib))
        wk(3,ib) = danrm2(bas(1,ib))
        call dcopy(3,bas(1,ib),1,wk(4,ib),1)
      enddo

      if (mode == 0) return

      call dvheap(6,nbas,wk,iprm,1d-8,1)
      do  ib = 1, nbas
        ibp = iprm(ib)
        ips(ib)    = idnint(wk(1,ibp))
        iclass(ib) = idnint(wk(2,ibp))
        call dcopy(3,wk(4,ibp),1,bas(1,ib),1)
      enddo

      end

      subroutine rsmesh(avstep,ierr,irr,iwtmp,jmax,jmin,n1,n2,n3,nirr,
     .                  nirrmx,nsymop,plat,qlat,rirr,symopm,symopv)
C-Divides the unit cell into microcells and finds irreducible part
C-----------------------------------------------------------------------
Ci Inputs:
Ci   n1-3  :no of divisions made along each primitive lattice vector
Ci   nirrmx:maximum number of irreducible points
Ci   nsymop:number of symmetry operations
Ci   plat  :primitive lattice vectors
Ci   qlat  :primitive translation vectors in reciprocal space
Ci   symopm:symmetry operation matrix
Ci   symopv:symmetry operation vector (scaled by alat)
Cio Inputs/Outputs
Co   irr   :irr(i1,i2,i3) points to corresponding irreducible point
Co Outputs:
Co   avstep:average step length
Co   ierr  :1 if nirr > nirrmx
Co   iwtmp :weight of mesh-point
Co   jmax  :upper limit for coordinates of irreducible points
Co   jmin  :lower limit for coordinates of irreducible points
Co   nirr  :number of irreducible points
Co   rirr  :irreducible point
Cr Remarks:
Cr  The lattice is divided into nrxyz(1)*nrxyz(2)*nrxyz(3)
Cr  microcells which are parallelipipeds with 8 corners.
Cr  The corners are nodes of the real space mesh in the whole
Cr  lattice unit cell. Some of these will be symmetry-related leaving
Cr  nirr irreducible k-points.
Cr  These are returned in rirr(3,j) j = 1,nirr; for each corner defined
Cr  by the triple (i1,i2,i3), irr(i1,i2,i3) points to the corresponding
Cr  vector in rirr.
Cr  It is possible to operate with integer algebra, if the mesh
Cr  is invariant under the symmetry operations considered.
C-----------------------------------------------------------------------
      implicit none
C Passed variables:
      integer ierr,n1,n2,n3,irr(0:n1-1,0:n2-1,0:n3-1),iwtmp(*),jmin(3),
     .        jmax(3),nirr,nirrmx,nsymop,stdo,nglob
      double precision avstep,plat(3,3),qlat(3,3),symopv(3,*),
     .                 symopm(9,*),rirr(3,*)
C Local variables:
      integer i1,i2,i3,iirr,iprint,isop,ism(9,48),isv(3,48),j1,
     .        j2,j3,jj,k,k1,k2,k3,kk,m,n0,nnc,nrxyz(3),n12,n22,n32
      double precision dasum,diff(12),dmxdev,dnrm2,qb(3,3),rb(3,3),
     .        step(3),tiny,toldev,tmp1(9),tmp2(9)
      logical lcomp(48)
      character*144 messg
      parameter(toldev=1d0,tiny=1.d-3)
C External calls:
      external  dmpy
C               ,errmsg,iprint
C Intrinsic functions:
      intrinsic  dmax1,dsqrt,idnint,max0,min0,mod

      nrxyz(1)=n1
      nrxyz(2)=n2
      nrxyz(3)=n3
      n12=n1/2
      n22=n2/2
      n32=n3/2
      n0=n1*n2*n3

      stdo=nglob('stdo')

C --- Mesh-size ...
      if (iprint() >= 20) write(stdo,298)nrxyz
      do  m = 1, 3
        step(m)=dnrm2(3,plat(1,m),1)/nrxyz(m)
      enddo

      avstep=(step(1)*step(2)*step(3))**(1.d0/3.d0)
      dmxdev=dmax1(step(1)/avstep,step(2)/avstep,step(3)/avstep)-1.d0
      if (dmxdev > toldev .and.iprint() >= 30) then
        write(messg,402)step
        call errmsg(messg,0)
      endif

      do  m = 1, 3
        jmin(m)= nrxyz(m)
        jmax(m)=-nrxyz(m)
        do  k = 1, 3
          rb(m,k) = plat(m,k)/nrxyz(k)
          qb(m,k) = qlat(m,k)*nrxyz(k)
        enddo
      enddo
      if (iprint() >= 60) then
        write(stdo,300)((rb(m,k),m=1,3),(qb(m,k),m=1,3),k=1,3)
      endif

C --- check the compatiblity of nrxyz with symgrp

      nnc=0
ct      write(*,*) 'nsymop',nsymop
      do  isop = 1, nsymop
C         write(*,*) 'isop=',isop
        call dmpy(qb,1,3,symopm(1,isop),3,1,tmp1,3,1,3,3,3)
        call dmpy(tmp1,3,1,rb,3,1,tmp2,3,1,3,3,3)
        do  m = 1, 9
          ism(m,isop)=idnint(tmp2(m))
          diff(m)=tmp2(m)-ism(m,isop)
        enddo

        call dmpy(qb,1,3,symopv(1,isop),3,1,tmp1,3,1,3,1,3)

        do  m = 1, 3
          isv(m,isop)=idnint(tmp1(m))
          diff(9+m)=tmp1(m)-isv(m,isop)
        enddo
        lcomp(isop)=dasum(12,diff,1) < tiny
        if (.not.lcomp(isop)) then
c        if (iprint() >= 50) write(lunit(2),400)nrxyz,qb,(symopm(k,isop),
c     .                                    k=1,9),rb,tmp2,(tmp1(k),k=1,3)
           write(*,296)
          nnc=nnc+1
C          write(*,*) 'isop=',isop
        endif
      enddo
ct      write(*,*) nnc

      if (nnc /= 0.and.ierr == 0) then
        write(messg,401)nnc,nsymop
        call errmsg(messg,0)
      endif
      write(*,297)

C ---

      nirr = 0
      do  i3 = 0, n3-1
        do  i2 = 0, n2-1
          do  i1 = 0, n1-1
            if (irr(i1,i2,i3) == 0) then
              nirr = nirr+1
              if (nirr > nirrmx) then
c                if (iprint() >= 70)
c     .            call errmsg(' RSMESH: increase nirrmx.$',-1)
                ierr=1
                return
              endif
              jj=n0
              do  isop = 1, nsymop
                if (lcomp(isop)) then
                  k1=isv(1,isop)+ism(1,isop)*i1+ism(4,isop)*i2
     .                          +ism(7,isop)*i3
                  k2=isv(2,isop)+ism(2,isop)*i1+ism(5,isop)*i2
     .                          +ism(8,isop)*i3
                  k3=isv(3,isop)+ism(3,isop)*i1+ism(6,isop)*i2
     .                          +ism(9,isop)*i3
                  k1 = mod(k1,n1)
                  k2 = mod(k2,n2)
                  k3 = mod(k3,n3)
                  if (k1 < 0) k1=k1+n1
                  if (k2 < 0) k2=k2+n2
                  if (k3 < 0) k3=k3+n3
                  if (irr(k1,k2,k3) == 0) then
                    iwtmp(nirr)=iwtmp(nirr)+1
                    irr(k1,k2,k3) = nirr
C ----------------- now determine j1,j2,j3
C ----------------- Attention: i1 in [0,n1-1] but k1 in [-n1/2,n1/2-1]
                    if (k1 >= n12) k1=k1-n1
                    if (k2 >= n22) k2=k2-n2
                    if (k3 >= n32) k3=k3-n3
                    kk=k1+k2+k3
                    if (kk < jj) then
                      jj=kk
                      j1=k1
                      j2=k2
                      j3=k3
                    endif
                  endif
                endif
              enddo
              jmin(1)=min0(j1,jmin(1))
              jmax(1)=max0(j1,jmax(1))
              jmin(2)=min0(j2,jmin(2))
              jmax(2)=max0(j2,jmax(2))
              jmin(3)=min0(j3,jmin(3))
              jmax(3)=max0(j3,jmax(3))
              rirr(1,nirr) = j1*rb(1,1) + j2*rb(1,2) + j3*rb(1,3)
              rirr(2,nirr) = j1*rb(2,1) + j2*rb(2,2) + j3*rb(2,3)
              rirr(3,nirr) = j1*rb(3,1) + j2*rb(3,2) + j3*rb(3,3)
c              if (iprint() >= 70)
c     .          write(stdo,301)i1,i2,i3,nirr,(rirr(k,nirr),k=1,3)
            endif
          enddo
        enddo
      enddo

      do  m = 1, 3
        if (iprint() >= 70) write(stdo,302)m,m,jmin(m),jmax(m)
      enddo

      call info5(30,0,0,'%9f%i irreducible mesh-points from %i ( %i %i %i )',
     .  nirr,n0,n1,n2,n3)
      if (iprint() >= 60) then
        write(stdo,304)
        do  iirr = 1, nirr
          write(stdo,305) iirr,(rirr(m,iirr),m=1,3),iwtmp(iirr)
        enddo
      endif

ct     do  m = 1, 3
ct        write(*,302)m,m,jmin(m),jmax(m)
ct     enddo

ct      write(*,303)nirr,n0,n1,n2,n3
ct      write(*,304)
C      do  iirr = 1, nirr
C         write(*,305)iirr,(rirr(m,iirr),m=1,3),iwtmp(iirr)
C      enddo

      ierr = 0

  296 format(T10,'nrxyz not compatible with symgrp')
  297 format(T10,'Compatibility of nrxyz with symgrp -- ok')
  298 format(/' RSMESH: NRXYZ= ',3i4)
  300 format(/' RSMESH:',11x,' RB ',31x,' QB '/3(3x,3f10.5,5x,3f11.5/))
C 301 format(3i5,i8,2x,3f6.3)
  302 format( ' RSMESH: JMIN',i1,', JMAX',i1,'=',2i4)
  303 format(8x,i6,' irreducible mesh-points from',i8,
     .  ' (',3i4,' )')
  304 format(/' RSMESH:  NIRR',7x,'Rx',8x,'Ry',8x,'Rz',6x,'IWTMP',
     .  /9x,46('-'))
  305 format(9x,i5,2x,3f10.4,5x,i2)

C  400 format(' RSMESH: WARNING: NRXYZ=',3i4,
C     .       ' not compatible with symmetry:',
C     .      //24x,3f9.4,/8x,'QB=QLAT*NRXYZ= ',1x,3f9.4,/24x,3f9.4/,
C     .       /24x,3f9.4,/8x,'SYMOPM=        ',1x,3f9.4,/24x,3f9.4/,
C     .       /24x,3f9.4,/8x,'RB=PLAT/NRXYZ  ',1x,3f9.4,/24x,3f9.4/,
C     .       /24x,3f9.4,/8x,'QB^t*SYMOPM*RB=',1x,3f9.4,/24x,3f9.4/,
C     .                  /8x,'QB^t*SYMOPV    ',1x,3f9.4/,
C     .       /' and the last 12 values should be integer values.'/)
  401 format(' RSMESH: mesh not invariant under ',i2,
     .       ' symmetry operations|        (from ',i2,')$')
  402 format(' RSMESH: bad NRXYZ:',
     .       ' different step length along lattice vectors:',
     .       '|',8x,'STEP1=',f7.4,',  STEP1=',f7.4,',  STEP1=',f7.4,'$')
      end

      subroutine reducv(plat,qlat,vec,flag)
C- Removes real space lattice vectors from a given real space vector
C ----------------------------------------------------------------------
Ci Inputs:
Ci   qlat  :primitive translation vectors in reciprocal space
Ci   plat  :primitive translation vectors in real space
Ci          they must fulfil: sum_k plat_ki * qlat_kj = delta_ij
Cio  vec   :real space test-vector (cartesian components)
Ci   flag  :controls form of vec at output
Ci          1: vec is given as coordinates of plat
Ci          2: vec is given in cartesian coordinates (vec in unit-cell)
Ci          3: vec is given in cartesian coordinates (vec in WS-cell)
Ci          4: vec is given in cartesian coordinates and is a multiple
Ci             of plat/12 -> all symopv fulfill this condition for
Ci             a standard choice of the origin.
Co Outputs:
Cio  vec   :real space test-vector
Co          vector with shortest length after removal of lattice vectors
Cr Remarks:
C ----------------------------------------------------------------------
      implicit none
C Passed variables:
      integer flag
      double precision plat(9),qlat(9),vec(3)
C Local variables:
      integer i,i1,i2,i3
      double precision danrm2,d1mach,tolpos,v2,vdiff(3),vmin,
     .       v(3),v01,v02,v03,v11,v12,v13,v21,v22,v23
      parameter(tolpos=4.d-3)
C External calls:
      external danrm2,d1mach,dcopy,dpzero
C Intrinsic functions:
      intrinsic dnint

      vdiff(1)=vec(1)*qlat(1)+vec(2)*qlat(2)+vec(3)*qlat(3)
      vdiff(2)=vec(1)*qlat(4)+vec(2)*qlat(5)+vec(3)*qlat(6)
      vdiff(3)=vec(1)*qlat(7)+vec(2)*qlat(8)+vec(3)*qlat(9)
      vdiff(1)=vdiff(1)-dnint(vdiff(1)-d1mach(3))
      vdiff(2)=vdiff(2)-dnint(vdiff(2)-d1mach(3))
      vdiff(3)=vdiff(3)-dnint(vdiff(3)-d1mach(3))
      if (flag == 1) then
        call dcopy(3,vdiff,1,vec,1)
        return
      endif

      if (flag == 4) then
        do  i = 1, 3
          v(i)=dnint(12.d0*vdiff(i))/12.d0
          if (dabs(vdiff(i)-v(i)) < tolpos) vdiff(i)=v(i)
        enddo
      endif

C --- Now the coordinates of the reduced vectors in the basis plat
C --- are all between -.5 and + .5
      vec(1)=vdiff(1)*plat(1)+vdiff(2)*plat(4)+vdiff(3)*plat(7)
      vec(2)=vdiff(1)*plat(2)+vdiff(2)*plat(5)+vdiff(3)*plat(8)
      vec(3)=vdiff(1)*plat(3)+vdiff(2)*plat(6)+vdiff(3)*plat(9)

      if (flag == 3) then
C ----- Try shortening by adding basis vectors
        vmin=d1mach(2)
        v01=vec(1)
        v02=vec(2)
        v03=vec(3)
        do  i1 = -1, 1
          v11=v01+i1*plat(1)
          v12=v02+i1*plat(2)
          v13=v03+i1*plat(3)
          do  i2 = -1, 1
            v21=v11+i2*plat(4)
            v22=v12+i2*plat(5)
            v23=v13+i2*plat(6)
            do  i3 = -1, 1
              v(1)=v21+i3*plat(7)
              v(2)=v22+i3*plat(8)
              v(3)=v23+i3*plat(9)
              v2=danrm2(v)
              if (v2 < vmin) then
                vmin=v2
                call dcopy(3,v,1,vec,1)
              endif
            enddo
          enddo
        enddo
      endif
      end
