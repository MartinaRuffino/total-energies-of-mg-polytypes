C#define SGI
      subroutine hsibq(s_site,s_spec,s_lat,k1,k2,k3,vsm,isp,q,ndimh,
     .  iprmb,napw,igapw,h)
C- Interstitial ME of smooth Bloch Hankels, smooth potential.
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  hsubblock hsibq2 hsibq4
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa lmxb pz name orbp ngcut
Co     Stored:     orbp
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  tbhsi uspecb hsubblock hsibq2 hsibq4
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat plat qlat gmax nabc ng vol kv kv2
Co     Stored:     *
Co     Allocated:  gv kv igv kv2 igv2
Cio    Elts passed:gv igv kv igv2
Cio    Passed to:  sugvec
Ci Inputs
Ci   k1,k2,k3 dimensions of vsm
Ci   vsm   :smoothed potential, real-space mesh
Ci   isp   :current spin channel (1 or 2)
Ci   q     :Bloch vector
Ci   ndimh :dimension of hamiltonian
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   napw  :number of augmented PWs in basis
Ci   igapw :vector of APWs, in units of reciprocal lattice vectors
Co Outputs
Co   h     :interstitial matrix elements of vsm added to h
Cr Remarks
Cr  *How orbital is extracted and employed.
Cr   See Remarks in smhsbl.f
Cm MPI
Cm   Parallelise over the main loop over nbas. In the serial code, h
Cm   is added to in each pass through the loop. Furthermore h is non
Cm   zero on entry to hsibq. This leads to a problem for MPI because
Cm   each process cannot simply add to the existing array h and pool
Cm   results with ALLREDUCE: this would lead to double counting. Instead
Cm   the increment to h from each call to hsibq must be pooled after
Cm   the main loop and then added to h. This requires allocating a
Cm   workspace of the same dimension as h. A second workspace of the
Cm   same length is needed as a buffer for ALLREDUCE. This second
Cm   work array can be dispensed with once MPI-2 is implemented with
Cm   the MPI_IN_PLACE feature of ALLREDUCE. Because these two work
Cm   arrays are large, they are taken from the heap rather than
Cm   ALLOCATEd using F90. Note that allocation of one work array the
Cm   size of h from the heap does not increase memory load because the
Cm   workspace for the eigenvectors is not yet allocated.
Cu Updates
Cu   17 Jun 13 Replace f77 pointers with f90 ones
Cu   05 May 12 Repackaged G vector generation into sugvec
Cu   10 Jul 10 Adapted from hsibl.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C#ifdefC MPI
C#ifdefC MPE
C      include "mpef.h"
C#endifC
C      include "mpif.h"
C      integer numprocs, ierr, status(MPI_STATUS_SIZE)
C      integer MAX_PROCS
C      parameter (MAX_PROCS = 100)
C      integer resultlen
C      character*(MPI_MAX_PROCESSOR_NAME) name
C      character*10 shortname(0:MAX_PROCS-1)
C      character*20 ext
C      character*26 datim
C      integer namelen(0:MAX_PROCS-1)
C      double precision starttime, endtime
C      integer lgunit
C#endif
C#ifdefC MPE
CC Event numbers:
C      include "events.ins"
C#endif
C ... Passed parameters
      integer k1,k2,k3,isp,ndimh,iprmb(*),napw,igapw(3,napw)
      double precision q(3)
      double complex h(ndimh,ndimh),vsm(k1,k2,k3,isp)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
C ... Local parameters
      integer procid, master, mpipid
      logical mlog,cmdopt
      character*120 strn
      integer n0,npmx,nkap0,nkape,nermx,nspec,ngabc(3),nhblk,nlmto
      parameter (n0=10,nkap0=4,nkape=2,nermx=100)
      parameter (npmx=128,nhblk=60)
      integer ltop,nbas,n1,n2,n3,net,ng,nglob,nlmtop,nrt,
     .  ndimx,nnrl,iprint
      double precision alat,plat(3,3),qlat(3,3),vol,gmax,q0(3)
      double precision etab(nermx),rtab(nermx)
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))

C     Shared variables
      integer i,j,nenv
      integer iprt(n0,nkap0,nermx),ipet(n0,nkap0,nermx)
      double precision xx(n0)
      integer ig1,i1,ig2,i2,igx(3),igx1,igx2,igx3
      integer,allocatable :: ncuti(:)
      complex(8),allocatable :: cv(:,:),c2(:,:)
C     For dividing into blocks
      integer block1,nenv1,nblocklist1,block2,nenv2,nblocklist2,
     .  nblock,mxbldim
      integer,allocatable :: blocklist1(:),blocklist2(:),
     .  nenvi1(:),offh1(:),nenvi2(:),offh2(:)
      real(8),allocatable :: g(:),yl(:),g2(:),he(:),hr(:)

C$    integer mp_numthreads,mp_my_threadnum
C#ifdefC MPI
C      double precision, dimension(:), allocatable :: buffer
C      integer, dimension(:,:), allocatable :: index
C      complex(8), allocatable :: hl(:),hbuf(:)
C      integer iloop
C#endif

      call tcn('hsibq')

C     print *, '!!'; h=0
C     print *, '!!'; q(1) = .1d0

      procid = mpipid(1)
      master = 0
      mlog = cmdopt('--mlog',6,0,strn)

C ... First setup
      nbas  = nglob('nbas')
      nspec = nglob('nspec')
      nlmto = ndimh - napw
      if (nspec > nermx) call rx('hsibq: increase nermx')
      alat = s_lat%alat
      plat = s_lat%plat
      qlat = s_lat%qlat
      gmax = s_lat%gmax
      ngabc = s_lat%nabc
      ng = s_lat%ng
      vol = s_lat%vol

C --- <MTO | V | MTO and < MTO | V PW> parts of h ---
      if (nlmto > 0) then

      call pshpr(iprint()-30)
      call sugvec(s_lat,5009,q,0d0,gmax,ngabc,0,ng)
      call poppr

C      call prmx('gv',s_lat%gv,ng,ng,3)
C      call prmi('kgv',s_lat%kv,ng,ng,3)
C      call prmi('igv',s_lat%igv,ng,ng,3)


C ... Tables of energies, rsm, indices to them
      call tbhsi(s_spec,nspec,nermx,net,etab,ipet,nrt,rtab,iprt,ltop)
C     ndimx = maximum hamiltonian dimension for any site
      ndimx = nnrl(10,1,nbas,iprmb,nlmto)

C --- Allocate and occupy arrays for YL, energy factors, rsm factors ---
      nlmtop = (ltop+1)**2
      allocate(g(ng*3),yl(ng*nlmtop),g2(ng),he(ng*net),hr(ng*nrt))
C     gv has q already added ...
      call dpzero(q0,3)
      call hsibq1(net,etab,nrt,rtab,ltop,alat,q0,ng,s_lat%gv,g,g2,yl,he,hr)

C --- Group hamiltonian by blocks ---
      mxbldim = 60
C      call hsubblock(0,s_spec,s_site,
C     .  (/0,4,2,1,5,3/),nbas,nlmto,
C     .  mxbldim,iprmb,nblock,nenv,blocklist1,nblocklist1)
      call hsubblock(0,s_spec,s_site,0,nbas,nlmto,mxbldim,iprmb,nblock,nenv,xx,nblocklist1)

C --- Loop over orbital blocks on first site ---
      allocate(ncuti(nlmto))
      allocate(blocklist1(nbas),blocklist2(nbas))
C#ifdefC MPI
CC      allocate(hl(ndimh*ndimh)); call dpzero(hl,2*ndimh*ndimh)
CC      allocate (index(0:numprocs-1,0:nbas-1), stat=ierr)
C#endif
      do  block1 = 1, nblock
        call hsubblock(1,s_spec,s_site,0,nbas,nlmto,mxbldim,iprmb,block1,nenv1,blocklist1,nblocklist1)
        allocate(cv(ng,nenv1),nenvi1(nblocklist1),offh1(nblocklist1))
C   ... Dimensioning parameters
        call hsibq2(s_spec,s_site,blocklist1,nblocklist1,0,xx,nlmto,iprmb,nenvi1,offh1,xx)
C   ... PW expansion of orbitals in block1, for <block1|vsm|block2>
        call hsibq4(s_spec,s_site,vol,blocklist1,nblocklist1,ng,ngabc,q,
     .    qlat,s_lat%igv,nlmto,iprmb,yl,ipet,iprt,etab,rtab,he,hr,cv)
C       call zprm('c',2,cv,ng,ng,nenv1)

C   ... Multiply potential into wave functions
        call hsibq7(n1,n2,n3,k1,k2,k3,vsm(1,1,1,isp),ng,s_lat%kv,nenv1,cv)
C       call zprm('cv',2,cv,ng,ng,nenv1)

C   --- Loop over second of (block1,block2) pairs ---
        do  block2 = block1, nblock
          call hsubblock(1,s_spec,s_site,0,nbas,nlmto,mxbldim,iprmb,block2,nenv2,blocklist2,nblocklist2)
          allocate(c2(ng,nenv2),nenvi2(nblocklist2),offh2(nblocklist2))
          call hsibq2(s_spec,s_site,blocklist2,nblocklist2,ng,g2,nlmto,iprmb,nenvi2,offh2,ncuti)
          call hsibq4(s_spec,s_site,vol,blocklist2,nblocklist2,ng,ngabc,
     .      q,qlat,s_lat%igv,nlmto,iprmb,yl,ipet,iprt,etab,rtab,he,hr,c2)

C     ... Add <block1|vsm|block2> into h
          call hsibq8(.true.,block1==block2,
     .      nblocklist1,blocklist1,nenv1,nenvi1,offh1,
     .      nblocklist2,blocklist2,nenv2,nenvi2,offh2,
     .      ng,ncuti,cv,c2,ndimh,h)

          deallocate(c2,nenvi2,offh2)
C       Ends loop over block2
        enddo

C   ... Matrix elements <Hsm| Vsm |PW>
C#ifdefC MPI
C        call hsibq6(nblocklist1,nenv1,nenvi1,offh1,ndimh,nlmto,
C     .    ng,s_lat%igv,napw,igapw,cv,hl)
C        call rx('hsibq not implemented for MPI yet')
C#else
        call hsibq6(nblocklist1,nenv1,nenvi1,offh1,ndimh,nlmto,
     .    ng,s_lat%igv,napw,igapw,cv,h)
C#endif

C     Ends loop over block1
      deallocate(cv,nenvi1,offh1)
      enddo
      deallocate(ncuti,blocklist1,blocklist2)

C     call zprm('h',2,h,ndimh,ndimh,ndimh)

C#ifdefC MPI
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_END_HSIBQ,procid,"hsibq")
C      ierr = MPE_LOG_EVENT(EVENT_START_BARRIER,procid,"barrier")
C#endifC
C      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_END_BARRIER,procid,"barrier")
C      ierr = MPE_LOG_EVENT(EVENT_START_ALLRED,procid,"allreduce")
C#endifC
C      allocate(hbuf(ndimh*ndimh))
C      call MPI_ALLREDUCE(hl,hbuf,2*ndimh*ndimh,
C     .     mpi_real8,MPI_SUM,MPI_COMM_WORLD,ierr)
C      if (mlog) then
C        call gettime(datim)
C        call awrit3(' hsibq '//datim//' Process %i of %i on '
C     .        //shortname(procid)(1:namelen(procid))//
C     .        ' allreduce h ndimh=%i',' ',256,lgunit(3),
C     .        procid,numprocs,ndimh)
C      endif
C      call daxpy(2*ndimh*ndimh,1d0,hbuf,1,h,1)
C
C      deallocate(hl,hbuf)
C      deallocate(index, stat=ierr)
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_END_ALLRED,procid,"allreduce")
C#endifC
C#endif
      deallocate(g,yl,g2,he,hr)
      endif

C --- <e^i qpG | V |e^i qpG'>/vol = V(G'-G) ---
      if (napw > 0) then
        call fftz3(vsm(1,1,1,isp),n1,n2,n3,k1,k2,k3,1,0,-1)
        do  ig1 = 1, napw
        i1 = ig1+nlmto
        do  ig2 = ig1, napw
          i2 = ig2+nlmto
          igx = igapw(:,ig1) - igapw(:,ig2)
          igx1 = mod(igx(1)+10*n1,n1)
          igx2 = mod(igx(2)+10*n2,n2)
          igx3 = mod(igx(3)+10*n3,n3)
          if (igx1<0 .or. igx2<0 .or. igx3<0) call rx('igx<0')
          h(i1,i2) = h(i1,i2) + vsm(igx1+1,igx2+1,igx3+1,isp)
        enddo
        enddo
        call fftz3(vsm(1,1,1,isp),n1,n2,n3,k1,k2,k3,1,0,1)
      endif

C ... Occupy second half of matrix
      do  i = 1, ndimh
        do  j = i, ndimh
          h(j,i) = dconjg(h(i,j))
        enddo
      enddo
      call tcx('hsibq')
C     call zprm('h-i',2,h,ndimh,ndimh,ndimh)
      end

      subroutine hsibq1(net,et,nrt,rt,ltop,alat,q,ng,gv,g,g2,yl,he,hr)
C- Make yl's, energy and rsm factors for list of G vectors
C ----------------------------------------------------------------------
Ci Inputs
Ci   net   :size of table et
Ci   et    :table of all inequivalent energies
Ci   nrt   :size of table rt
Ci   rt    :table of all inequivalent smoothing radii
Ci   ltop  :largest l at any site
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   q     :Bloch wave number
Ci   ng    :number of G-vectors
Ci   gv    :list of reciprocal lattice vectors G (gvlist.f)
Co Outputs
Co   g     :2*pi/alat * (q+gv) for all g-vectors
Co   g2    :g**2
Co   yl    :Y_L
Co   he    :1/(et-g2) for all inequivalent e's and g-vectors
Co   hr    :dexp(-(rsm/2)**2*g2(i)) for all inequivalent rsm and g-vecs
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ltop,net,ng,nrt
      double precision alat,q(3),gv(ng,3),g(ng,3),yl(ng,*),he(ng,net),
     .  hr(ng,nrt),g2(ng),et(net),rt(nrt)
C ... Local parameters
      integer i,ie,ir
      double precision pi,tpiba,gam

C ... Make (2*pi/alat)*(gv+q) in g
      pi = 4d0*datan(1d0)
      tpiba = 2d0*pi/alat
      do  i = 1, ng
        g(i,1) = tpiba*(gv(i,1)+q(1))
        g(i,2) = tpiba*(gv(i,2)+q(2))
        g(i,3) = tpiba*(gv(i,3)+q(3))
      enddo

C ... Make the yl's and g2
      call ropyln(ng,g(1,1),g(1,2),g(1,3),ltop,ng,yl,g2)

C ... Make the energy factors
      do  ie = 1, net
        do  i = 1, ng
          he(i,ie) = 1d0/(et(ie)-g2(i))
        enddo
      enddo

C ... Make the rsm factors
      do  ir = 1, nrt
        gam = 0.25d0*rt(ir)*rt(ir)
        do  i = 1, ng
          hr(i,ir) = dexp(-gam*g2(i))
        enddo
      enddo
      end

      subroutine hsibq2(s_spec,s_site,iblst,nblst,ng,g2,nlmto,iprmb,
     .  nenv,offh,ncuti)
C- Return number of orbitals and G-cutoffs for range of sites
C ----------------------------------------------------------------------
Ci Inputs
Ci   s_spec :struct for species-specific data; see structures.h
Ci     Elts read: ngcut lmxa lmxb pz name orbp
Ci     Stored:    orbp
Ci     Passed to: uspecb
Ci   s_site :struct for site-specific data; see structures.h
Ci     Elts read: spec
Ci     Stored:    *
Ci     Passed to: *
Ci   iblst :list of sites for which to make coffs
Ci   nblst :number of sites
Ci   ng    :number of G-vectors.  If nonzero, make ncuti
Co   g2    :[2*pi/alat * (q+gv)]**2 for all g-vectors
Ci   nlmto :dimension of lmto component of basis
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Co Outputs
Co   nenv  :nenv(i) = cumulative number of envelope functions
Co         :for orbitals in sites iblst(1)...iblst(i)
Co   offh  :offh(i) = offset to true hamiltonian at site iblst(i)
Co   ncuti :G-cutoff for each orbital (generated only if ng>0)
Cl Local variables
Cl   norb  :number of orbital types for last ib
Cl   offl  :offl(norb) offset in h to this block of orbitals, last ib
Cl   blks  :blks(iorb) = size of contiguous block of orbitals (gtbsl1)
Cr Remarks
Cr
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   05 Jul 10 Adapted from hsibl
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nblst,iblst(nblst),ng,nlmto,iprmb(*),ncuti(nlmto),
     .  offh(nblst),nenv(nblst)
      integer n0,nkap0
      parameter (n0=10,nkap0=4)
      double precision g2(ng)
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
C ... Local parameters
      integer ib,ibl,is,nlm1,nlm2,ik,l,lt,iorb,nenvi,nkap,jorb,ni,ig,
     .  offhi,nenvl
      integer lh(nkap0),ncut(n0,nkap0)
      integer norb,ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),
     .  blks(n0*nkap0),ntab(n0*nkap0)
      double precision eh(n0,nkap0),rsmh(n0,nkap0)
      double precision xx(n0),eps
      parameter (eps=1d-5)

      nenvl = 0
      do  ibl = 1, nblst
        ib = iblst(ibl)
        nenvi = 0
        is = s_site(ib)%spec
        ncut = s_spec(is)%ngcut
C       List of orbitals, their l- and k- indices, and ham offsets
        call orbl(ib,0,nlmto,iprmb,norb,ltab,ktab,offhi,offl,xx)
C       Block routines into groups with common (e,rsm)
        call uspecb(0,2,s_spec,is,is,lh,rsmh,eh,nkap)
        call gtbsl1(7+16,norb,ltab,ktab,rsmh,eh,ntab,blks)
        do  iorb = 1, norb
          if (blks(iorb) /= 0) then
            jorb = ntab(iorb)
            lt   = ltab(jorb)
            l    = ltab(iorb)
            ik   = ktab(iorb)
            nlm1 = l**2+1
            nlm2 = nlm1 + blks(iorb)-1
C           Increase cutoff to include all vectors of the same length
            if (ng > 0) then
              ni = min(ncut(lt+1,ik),ng)
              do  ig = ni+1, ng
                if (abs(g2(ig)-g2(ni)) > eps) then
                  ni = ig-1
                  exit
                endif
              enddo
              call ivset(ncuti,nenvl+nenvi+1,nenvl+nenvi+nlm2-nlm1+1,ni)
            endif
            nenvi = nenvi + max(nlm2-nlm1+1,0)
          endif
        enddo
        nenvl = nenvl + nenvi
        nenv(ibl) = nenvl
        offh(ibl) = offhi
      enddo

      end

      subroutine hsibq4(s_spec,s_site,vol,iblst,nblst,ng,ngabc,q,qlat,
     .  igv,nlmto,iprmb,yl,ipet,iprt,etab,rtab,he,hr,c)
C- Return plane wave expansion coefficients of LMTO's
C ----------------------------------------------------------------------
Ci Inputs
Ci   s_spec :struct for species-specific data; see structures.h
Ci     Elts read: lmxa lmxb pz name orbp
Ci     Stored:    orbp
Ci     Passed to: uspecb
Ci   s_site :struct for site-specific data; see structures.h
Ci     Elts read: spec pos
Ci     Stored:    *
Ci     Passed to: *
Ci   iblst :list of sites for which to make coffs
Ci   nblst :number of sites
Ci   ng    :number of G-vectors which sm. Hankel is expanded
Ci   ngabc :number of divisions in mesh
Ci   q     :k-point
Ci   qlat  :primitive reciprocal lattice vectors, in units of 2*pi/alat
Ci   nlmto :dimension of lmto component of basis
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   igv   :G-vectors, as integer multiples of qlat (gvlist)
Ci   yl    :spherical harmonics for ng vectors (ropyln)
Ci   ipet  :index to which entry in et a given orbital belongs (tbhsi)
Ci   iprt  :index to which entry in rt a given orbital belongs (tbhsi)
Ci   etab  :table of all inequivalent energies (tbhsi)
Ci   rtab  :table of all inequivalent smoothing radii (tbhsi)
Ci   he    :table of energy factors  1/(e-(G+q)^2) (hsibq1)
Ci   hr    :table of smoothing radius factors dexp(-(rsm/2)**2*g2(i)) (hsibq1)
Co Outputs
Co   c     :c(ig,ilm+offlm) = cfac phase he(ig,ie) hr(ig,ir) yl(ig,ilm)
Co         :where cfac  = -4*pi*(-1)^l*exp(e*rsm**2/4)/dsqrt(vol)
Co         :      phase = exp(-i pos * (q+G))
Co         :Written out, c is (c.f. JMP 39, 3393, Eq 6.17)
Co         :      -4*pi * (-1)^l
Co         : --------------------- exp[(rsm/2)^2 (e-(G+q)^2)] Y_L((G+q))
Co         : sqrt(vol) (e-(G+q)^2)
Cl Local variables
Cl  nenv   :running offset to current block in c
Cr Remarks
Cr
Cu Updates
Cu   10 Nov 11 Begin migration to f90 structures
Cu   05 Jul 10 Adapted from hsibl
C ----------------------------------------------------------------------
      use structures
      implicit none
C ... Passed parameters
      integer nblst,iblst(nblst),ng,ngabc(3),igv,nlmto,iprmb(*)
      integer n0,nkap0
      parameter (n0=10,nkap0=4)
      integer iprt(n0,nkap0,*),ipet(n0,nkap0,*)
      double precision q(3),qlat(3,3),vol
      double precision yl(ng,*),he(ng,*),hr(ng,*),etab(*),rtab(*)
      double complex c(ng,*)
C ... For structures
!      include 'structures.h'
      type(str_spec)::  s_spec(*)
      type(str_site)::  s_site(*)
C ... Local parameters
      integer ib,ie,ir,is,nlm1,nlm2,ik,l,iorb,nenv,nkap,ibl
      integer lh(nkap0)
      integer norb,ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),
     .  blks(n0*nkap0),ntab(n0*nkap0)
      double precision eh(n0,nkap0),rsmh(n0,nkap0)
      double precision xx(n0),pos(3)
      real(8),allocatable :: cs(:),sn(:),wk(:)

      allocate(cs(ng),sn(ng),wk(ng))
      nenv = 0
      do  ibl = 1, nblst
        ib = iblst(ibl)
        is = s_site(ib)%spec
        pos = s_site(ib)%pos
        call suphas(q,pos,ng,igv,ngabc(1),ngabc(2),ngabc(3),qlat,cs,sn)
C       List of orbitals, their l- and k- indices, and ham offsets
        call orbl(ib,0,nlmto,iprmb,norb,ltab,ktab,xx,offl,xx)
C       Block routines into groups with common (e,rsm)
        call uspecb(0,2,s_spec,is,is,lh,rsmh,eh,nkap)
        call gtbsl1(7+16,norb,ltab,ktab,rsmh,eh,ntab,blks)
        do  iorb = 1, norb
          if (blks(iorb) /= 0) then
            l    = ltab(iorb)
            ik   = ktab(iorb)
            nlm1 = l**2+1
            nlm2 = nlm1 + blks(iorb)-1
            ie   = ipet(l+1,ik,is)
            ir   = iprt(l+1,ik,is)
            call hsibq3(ie,ir,etab,rtab,vol,nlm1,nlm2,nenv,ng,yl,he,hr,cs,sn,wk,c)
            nenv = nenv + max(nlm2-nlm1+1,0)
          endif
        enddo
      enddo
      deallocate(cs,sn,wk)

      end

      subroutine hsibq3(ie,ir,etab,rtab,vol,nlm1,nlm2,offlm,ng,yl,he,hr,
     .  cosgp,singp,wk,c)
C- PW expansion of smooth Hankels at one site
C ----------------------------------------------------------------------
Ci Inputs
Ci   ie    :index to appropriate entry in energy factor table he
Ci   ir    :index to appropriate entry in sm radius factor table hr
Ci   etab  :table of all inequivalent energies
Ci   rtab  :table of all inequivalent smoothing radii
Ci   vol   :cell volume
Ci   nlm1  :starting L for which to accumulate FT
Ci   nlm2  :final L for which to accumulate FT
Ci   offlm :offset to ilm for storing c
Ci   ng    :number of G-vectors which sm. Hankel is expanded
Ci   yl    :spherical harmonics for ng vectors
Ci   he    :table of energy factors  1/(e-(G+q)^2)
Ci   hr    :table of smoothing radius factors dexp(-(rsm/2)**2*g2(i))
Ci   cosgp :table of phases (real part)
Ci   singp :table of phases (imaginary part)
Ci   wk    :real work array of length ng
Co Outputs
Co   c     :c(ig,ilm+offlm) = cfac phase he(ig,ie) hr(ig,ir) yl(ig,ilm)
Co         :where cfac  = -4*pi*(-1)^l*exp(e*rsm**2/4)/dsqrt(vol)
Co         :      phase = exp(-i pos * (q+G))
Co         :Written out, c is (c.f. JMP 39, 3393, Eq 6.17)
Co         :      -4*pi * (-1)^l
Co         : --------------------- exp[(rsm/2)^2 (e-(G+q)^2)] Y_L((G+q))
Co         : sqrt(vol) (e-(G+q)^2)
Cr Remarks
Cr   This routine requires O(M N) operations.
Cr   (M = number of orbitals, N = number of mesh points)
Cu Updates
Cu   10 Jul 10 Constant factors folded into c
Cu   22 May 00 Adapted from nfp su_hkft
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ie,ir,nlm1,nlm2,offlm,ng
      double precision yl(ng,*),he(ng,*),hr(ng,*),cosgp(*),singp(*),
     .  wk(ng),etab(ie),rtab(ir),vol
      double complex c(ng,*)
C ... Local parameters
      integer i,ilm,offi,lmax,ll,l,m
      double precision xxx,fac1,pi
      double complex cf,cfac
      parameter (pi = 3.1415926535897931d0)

      if (nlm2 == 0) return
      call tcn('hsibq3')
      do  i = 1, ng
        wk(i) = he(i,ie)*hr(i,ir)
      enddo
      offi = offlm-nlm1+1
      do  ilm = nlm1, nlm2
        do  i = 1, ng
          xxx = wk(i)*yl(i,ilm)
          c(i,ilm+offi) = dcmplx(xxx*cosgp(i),xxx*singp(i))
        enddo
      enddo

C     Constant factor
      fac1 = -4d0*pi*dexp(etab(ie)*rtab(ir)**2/4)/dsqrt(vol)
      lmax = ll(nlm2)
      ilm = 0
      cf = (0d0,1d0)
      do  l = 0, lmax
        cf = cf*(0d0,-1d0)
        do  m = -l, l
          ilm = ilm+1
          if (ilm >= nlm1) then
            cfac = cf * fac1
            call zscal(ng,cfac,c(1,ilm+offi),1)
          endif
        enddo
      enddo


      call tcx('hsibq3')

      end

      subroutine hsibq5(ndim1,ofh1,ofh2,ndim2,ndimh,hwk,h)
C- Adds a subblock of matrix elements into the hamiltonian
C ----------------------------------------------------------------------
Ci Inputs
Ci   ndim1 :leading row dimension to hwk
Ci   ndim2 :leading column dimension to hwk
Ci   ndimh :dimensions hamiltonian
Ci   hwk   :matrix elements of this block to be added to h
Co Outputs
Co   h     :matrix elements added to hamiltonian for this block
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   16 Aug 04 First created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ndimh,ndim1,ndim2,ofh1,ofh2
      double complex h(ndimh,ndimh),hwk(ndim1,ndim2)
C ... Local parameters
      integer i1,i2

      do  i1 = 1, ndim1
        do  i2 = 1, ndim2
          h(ofh1+i1,ofh2+i2) = h(ofh1+i1,ofh2+i2) + hwk(i1,i2)
        enddo
      enddo

      end

      subroutine hsibq6(nblock1,nenv1,nenvi1,offh1,
     .  ndimh,nlmto,ng,igv,napw,igapw,cv,h)
C- Make matrix elements <Hsm | V | PW>
C ----------------------------------------------------------------------
Ci   nblock1:number of sites in the phi1 block
Ci   nenv1  :total number of envelopes in the phi1 block
Ci   nenvi1 :cumulative number of envelope functions for orbitals
Ci          :for orbitals in sites iblst(1)...iblst(i).
Ci          (iblst is not passed here, as it is not reference,
Ci          but see hsibq8).  nenvi1 is related to the offsets
Ci          :n c1 for each site, and are used to make offsets
Ci          :to the entire block of envelope functions at all sites
Ci   offh1  :row offsets to the full h, the analog of nenvi1.
Ci          :nenvi1 and offh1 are needed to associate colums of c1
Ci          :with rows of h; see Remarks in hsibq8.
Ci   ndimh  :dimension of hamiltonian
Ci   nlmto  :number of lmto's
Ci   ng     :leading dimension of cv, and number G vectors in
Ci          :LMTO expansion
Ci   igv    :list of reciprocal lattice vectors G (from gvlist)
Ci   napw   :number of PW's
Ci   igapw  :vector of APWs, in units of reciprocal lattice vectors
Ci   cv     :coffs to PW expansion of Hsm * V
Co Outputs
Co   h     : <Hsm | V | PW> is added to h
Cl Local variables
Cl         :
Cr Remarks
Cb Bugs
Cb   ifindiv should be replaced with index passed through
Cu Updates
Cu   04 Jul 08 (T Kotani) first created
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nblock1,ng,napw,ndim1,ndimh,nenv1,nlmto,
     .  igv(ng,3),igapw(3,napw),nenvi1(nblock1),offh1(nblock1)
      double complex h(ndimh,ndimh),cv(ng,nlmto)
C ... Local parameters
      integer ibl1,ig,i2,i2x,ifindiv,ofw1,ofh1,i1
C     complex(8),pointer :: hh(:,:)

      do  ig = 1, napw
        i2  = nlmto+ig
C       index matching igv,igapw
        i2x = ifindiv(igapw(1,ig),igv,ng)
        do  ibl1 = 1, nblock1
C         ib1 = iblst1(ibl1)
          ofh1 = offh1(ibl1)
          ofw1  = 0
          if (ibl1 > 1) ofw1  = nenvi1(ibl1-1)
          ndim1 = nenvi1(ibl1) - ofw1
          do  i1 = 1, ndim1
            h(ofh1+i1,i2) = h(ofh1+i1,i2) + dconjg(cv(i2x,i1+ofw1))
          enddo
        enddo
      enddo

C ... Debugging: get from disk
C      open(55, form='UNFORMATTED')
C      rewind 55
C      allocate(hh(ndimh,ndimh))
C      call dpdump(hh,ndimh**2*2,55)
C      call diff (' ', h, hh, ndimh, ndimh, ndimh, ndimh)
C      deallocate(hh)

      end

      subroutine hsibq7(n1,n2,n3,k1,k2,k3,vsm,ng,kv,nc,c)
C- FFT to real space, multiply by potential, FTT back
C ----------------------------------------------------------------------
Ci Inputs
Ci   n1..3 :FT mesh
Ci   k1..3 :dimensions vsm,f,c
Ci   vsm   :potential
Ci   ng    :number of G-vectors
Ci   kv    :indices for gather/scatter operations (gvlist.f)
Ci   nc    :number of wave functions
Cio Inputs/Outputs
Cio  c     :On input, holds FT of wave function
Cio        :On output, holds FT of wave function * potential
Cr Remarks
Cr   This routine takes M N logN operations, and is often the most
Cr   time consuming for small to moderate systems.
Cr   (M = number of orbitals, N = number of mesh points)
Cr   FFTs dominate the CPU time.
Cu Updates
Cu   08 Jul 10 Adapted from hsibl4 to use more efficient FFT routines
Cu   22 May 00 Adapted from nfp shkpot
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n1,n2,n3,k1,k2,k3,ng,nc,kv(ng,3)
      double complex c(ng,nc),vsm(k1,k2,k3)
C ... Local parameters
      integer i,i1,i2,i3,plan1(2),plan2(2)
      complex(8),allocatable :: f(:,:,:)

      call tcn('hsibq7')
      allocate(f(k1,k2,k3))

C ..  Create forward and backwards FFT plans
      call fftzv(f,n1,n2,n3,1,1,01,1)
      call fftzvp(0,plan1)
      call fftzv(f,n1,n2,n3,1,1,01,-1)
      call fftzvp(0,plan2)

C ..  FT envelopes to real space, multiply by vsm, FT back
      do  i = 1, nc
        call gvputf(ng,1,kv,k1,k2,k3,c(1,i),f)
C       call fftz3(f,n1,n2,n3,k1,k2,k3,1,0,1)
        call fftzvp(1,plan1)
        call fftzv(f,n1,n2,n3,1,1,02,1)
C       call zprm3('psir',0,f,n1,n2,n3)
        do  i3 = 1, n3
          do  i2 = 1, n2
            do  i1 = 1, n1
              f(i1,i2,i3) = f(i1,i2,i3)*vsm(i1,i2,i3)
            enddo
          enddo
        enddo
c       call zprm3('v*psir',0,f,n1,n2,n3)
C       call fftz3(f,n1,n2,n3,k1,k2,k3,1,0,-1)
        call fftzvp(1,plan2)
        call fftzv(f,n1,n2,n3,1,1,02,-1)
        call gvgetf(ng,1,kv,k1,k2,k3,f,c(1,i))
      enddo

C ..  Destroy forward and backwards FFT plans
      call fftzvp(1,plan1)
      call fftzv(f,n1,n2,n3,1,1,4,1)
      call fftzvp(1,plan2)
      call fftzv(f,n1,n2,n3,1,1,4,1)

      deallocate(f)
      call tcx('hsibq7')
      end

      subroutine hsibq8(lcmpat,ldiag,nblock1,iblst1,nenv1,nenvi1,offh1,
     .  nblock2,iblst2,nenv2,nenvi2,offh2,ng,ncuti,c1,c2,ndimh,h)
C- Add ((phi1*vsm)*phi2) to h for a block of orbital pairs
C ----------------------------------------------------------------------
Ci Inputs
Ci   lcmpat :if true, make exactly compatible with hsibl2
Ci          :by zeroing out c2(i,iorb) for i>ncuti(iorb)
Ci   ldiag  :if true, the left and right blocks should be the same
Ci          :In thish case a special-purpose matrix multiplication
Ci          :is used.
Ci   nblock1:number of sites in the phi1 block
Ci   iblst1 :list of sites in the phi1 block (see Remarks)
Ci   nenv1  :total number of envelopes in the phi1 block
Ci   nenvi1 :cumulative number of envelope functions for orbitals
Ci          :for orbitals in sites iblst(1)...iblst(i).  These are
Ci          :related to the offsets in c1 for each site, and are
Ci          :used to make offsets to the entire block of envelope
Ci          :functions at all sites
Ci   offh1  :row offsets to the full h, the analog of nenvi1.
Ci          :nenvi1 and offh1 are needed to associate colums of c1
Ci          :with rows of h; see Remarks
Ci   nblock2:number of sites in the phi2 block
Ci   iblst2 :list of sites in the phi2 block
Ci   nenv2  :total number of envelopes in the phi2 block
Ci   nenvi2 :the equivalent of nenvi1 for the c2 coefficients.
Ci   offh2  :column offsets to the full h, the analog of nenvi2
Ci          :nenvi2 and offh2 are needed to associate colums of c1
Ci          :with colums of h; see Remarks
Ci   ng     :leading dimension of c1,c2, and number G vectors in
Ci          :LMTO expansion
Ci   ncuti  :orbital-dependent number of G-vectors
Ci   c1     :product (phi1 vsm) in reciprocal space
Ci   c2     :phi2 in reciprocal space
Ci   ndimh  :dimension of h
Co Outputs
Co   h      :scalar product (phi1 vsm) phi2 is added to h
Cr Remarks
Cr   This routine assembles matrix elements of the smooth potential
Cr   for a block of envelope function pairs.
Cr   A block of orbitals consists of envelope functions within a list
Cr   of sites: the left and right orbitals each have their own block.
Cr   Basis functions are blocked (ie grouped) in units larger than
Cr   one site to take advantage of fast zgemm products.
Cr   The total number of row (column) envelope functions nenv1 (nenv2)
Cr
Cr   Since the number of envelopes need not match the number of basis
Cr   functions at a site, hamiltonian offsets do not synchronize with
Cr   offsets in coefficients c1 (c2).  For that reason,
Cr   matrix elements must be generated in a temporary array and
Cr   added to h in smaller blocks, corresponding to orbitals for
Cr   a single pair of sites.  The bottom part of this routine does
Cr   this translation.
Cr
Cr   This routine takes M^2 N operations, and is the most
Cr   time consuming step for moderate to large systems.
Cr   (M = number of orbitals, N = number G vectors)
Cu Updates
Cu   10 Jul 10 Adapted from hsibl
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical lcmpat,ldiag
      integer nblock1,nblock2,ng,ndimh,nenv1,nenv2,ncuti(nenv2),
     .  nenvi1(nblock1),offh1(nblock1),iblst1(nblock1),
     .  nenvi2(nblock2),offh2(nblock2),iblst2(nblock2)
      complex(8) c1(ng,nenv1),c2(ng,nenv2),h(ndimh,ndimh)
C ... Local parameters
      integer i,i1,i2,ncut,nbar,itop,ib1,ibl1,ofh1,ofw1,ndim1,ib2,ibl2,
     .  ofh2,ofw2,ndim2
      double precision bar
      complex(8),allocatable :: hwk(:,:)
C     complex(8),allocatable :: hh(:,:)

      call tcn('hsibq8')

C     Find the average number of PWs in the block.
C     These enter into zgemm.  Contributions from orbitals
C     with more PWs are added with inline code.
      bar = 0
      do  i2 = 1, nenv2
        bar = bar + ncuti(i2)
      enddo
      nbar = nint(bar/nenv2)
      if (lcmpat) then
        do  i2 = 1, nenv2
          if (ncuti(i2) < nbar) then
C           print *, c2(ncuti(i2):nbar,i2)
            c2(ncuti(i2)+1:nbar,i2) = 0d0
C           print *, c2(ncuti(i2):nbar,i2)
          endif
        enddo
      endif

C --- All matrix elements <iblst1|vsm|iblst2> ---
      allocate(hwk(nenv1,nenv2))
      if (ldiag) then
        call zqsmpy(11,'C','N',nenv2,nbar,c1,ng,c2,ng,(0d0,0d0),hwk,nenv1)
      else
        call zgemm('C','N',nenv1,nenv2,nbar,(1d0,0d0),c1,ng,c2,ng,(0d0,0d0),hwk,nenv1)
      endif
C     Fill in contributions above nbar
      do  i2 = 1, nenv2
        ncut = min(ng,ncuti(i2))
        itop = nenv1
        if (ldiag) itop = i2
        do  i1 = 1, itop
          do  i = nbar+1, ncut
            hwk(i1,i2) = hwk(i1,i2) + dconjg(c1(i,i1))*c2(i,i2)
          enddo
        enddo
      enddo

C --- Distribute hwk to h ---
      do  ibl1 = 1, nblock1
        ib1 = iblst1(ibl1)
        ofh1 = offh1(ibl1)
        ofw1  = 0
        if (ibl1 > 1) ofw1  = nenvi1(ibl1-1)
        ndim1 = nenvi1(ibl1) - ofw1
        do  ibl2 = 1, nblock2
          ib2 = iblst2(ibl2)
          ofh2 = offh2(ibl2)
          ofw2  = 0
          if (ibl2 > 1) ofw2  = nenvi2(ibl2-1)
          ndim2 = nenvi2(ibl2) - ofw2
          do  i2 = 1, ndim2
            do  i1 = 1, ndim1
              h(i1+ofh1,i2+ofh2) = h(i1+ofh1,i2+ofh2) + hwk(i1+ofw1,i2+ofw2)
            enddo
          enddo
C          print "(2i4,2f15.10)",ib1,ib2,cmplx(h(ofh1+ndim1,ofh2+ndim2))
        enddo
      enddo

C ... Debugging: get from disk
C      open(55, form='UNFORMATTED')
C      rewind 55
C      allocate(hh(ndimh,ndimh))
C      call dpdump(hh,ndimh**2*2,55)
C      call diff (' ', h, hh, ndimh, ndimh, ndimh, ndimh)
C      deallocate(hh)

      deallocate(hwk)
      call tcx('hsibq8')
      end


      integer function ifindiv(igapw,igv,ng)
C- Find index in igv that corresponds to igapw
C ----------------------------------------------------------------------
Ci   igapw :vector of APWs, in units of reciprocal lattice vectors
Ci   igv   :List of G vectors
Ci   ng    :number of group operations
Co Outputs
Co   ifindiv:index to igv that matches igapw
Cu Updates
Cu   19 Jan 09 Original cleaned up, made more efficient
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ng,igapw(3),igv(ng,3)
C ... Local parameters
      integer ig

      do  ig = 1, ng
        if (igapw(1) /= igv(ig,1)) cycle
        if (igapw(2) /= igv(ig,2)) cycle
        if (igapw(3) /= igv(ig,3)) cycle
        ifindiv = ig
        return
      enddo

C      write(*,'(10i5)') igapw
C      do ig=1,ng
C        write(6,'(i5,2x,10i5)') ig, igv(ig,:)
C      enddo
      call rx('ifindiv: igapw not found in igv')

      end

      integer function ifindiv2(igapw,igv2,ng)
C- Find index in igv2 that corresponds to igapw
C ----------------------------------------------------------------------
Ci   igapw :vector of APWs, in units of reciprocal lattice vectors
Ci   igv2   :List of G vectors
Ci   ng    :number of group operations
Co Outputs
Co   ifindiv2:index to igv2 that matches igapw
Cu Updates
Cu   19 Jan 09 Original cleaned up, made more efficient
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer ng,igapw(3),igv2(3,ng)
C ... Local parameters
      integer ig

      do  ig = 1, ng
        if (igapw(1) /= igv2(1,ig)) cycle
        if (igapw(2) /= igv2(2,ig)) cycle
        if (igapw(3) /= igv2(3,ig)) cycle
        ifindiv2 = ig
        return
      enddo

C      write(*,'(10i5)') igapw
C      do ig=1,ng
C        write(6,'(i5,2x,10i5)') ig, igv2(ig,:)
C      enddo
      call rx('ifindiv2: igapw not found in igv')

      end
C
C      SUBROUTINE diff (strn, sc, pc, lds, ldp, n, m)
CC- Compare the two arrays for differences
C      implicit none
C      character*(*) strn
C      INTEGER lds, ldp, n, m, l, i, j
C      double complex sc( lds, * )
C      double complex pc( ldp, * )
C      double precision tol,errmx
C
C      tol = 1d-12
C
C      errmx = 0d0
C      DO i = 1, n
C        DO j = i, m
C          errmx = max(errmx,ABS( sc(i,j) - pc(i,j) ))
CC          IF ( ABS( sc(i,j) - pc(i,j) ) .GT. tol ) THEN
CC            WRITE (6,100) strn, i,j, ABS( sc(i,j) - pc(i,j) )
CC            RETURN
CC          END IF
C        END DO
C      END DO
C      WRITE (*,112) strn, errmx
C      RETURN
CC  100 FORMAT(1X,a,'*** ERROR ***  Arrays Have Different Results: i,j=',
CC     .  2i5,' err=',1pe8.1)
CC  110 FORMAT(1X,a,'... arrays have the same results',
CC     .  ' tol=',1pe8.1)
C  112 FORMAT(1X,a,'... arrays diff by max val ',1pe8.1)
C
C      END
