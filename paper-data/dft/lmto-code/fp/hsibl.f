      subroutine hsibl(s_site,s_spec,s_lat,k1,k2,k3,vsm,isp,q,ndimh,iprmb,napw,igapw,h)
C- Interstitial ME of smooth Bloch Hankels, smooth potential.
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  hsibq2 hsibq4
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa lmxb pz name orbp ngcut
Co     Stored:     orbp
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  tbhsi uspecb hsibq2 hsibq4
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
Cm   zero on entry to hsibl. This leads to a problem for MPI because
Cm   each process cannot simply add to the existing array h and pool
Cm   results with ALLREDUCE: this would lead to double counting. Instead
Cm   the increment to h from each call to hsibl must be pooled after
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
Cu   05 May 12 Repackaged G vector generation into sugvec
Cu   02 Jun 09 (T. Kotani) Added call to ncutcorrect
Cu   05 Jul 08 (T. Kotani) new APW part of basis
Cu   12 Aug 04 First implementation of extended local orbitals
Cu   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   27 Aug 01 Extended to local orbitals.
Cu             At present, must compile with SGI with local orbitals!
Cu   12 Oct 00 Use q-dependent list of G vectors
Cu   22 May 00 Adapted from nfp hsi_q.f
C ----------------------------------------------------------------------
      use structures
      implicit none
C#ifdefC MPI
C#ifdefC MPE
C      include "mpef.h"
C#endifC
C      include "mpif.h"
C      integer procid, master, numprocs, ierr, status(MPI_STATUS_SIZE)
C      integer MAX_PROCS
C      parameter (MAX_PROCS = 100)
C      integer resultlen
C      character*(MPI_MAX_PROCESSOR_NAME) name
C      character*10 shortname(0:MAX_PROCS-1)
C      character*20 ext
C      character*26 datim
C      integer namelen(0:MAX_PROCS-1)
C      double precision starttime, endtime
C      logical mlog
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
      integer n0,npmx,nkap0,nkape,nermx,nspec,ngabc(3),nhblk,nlmto
      parameter (n0=10,nkap0=4,nkape=2,nermx=100)
      parameter (npmx=128,nhblk=60)
      integer ltop,nbas,n1,n2,n3,net,ng,nglob,nlmtop,nrt,
     .  ndimx,nnrl,ncuti(nhblk),iprint
      double precision alat,plat(3,3),qlat(3,3),vol,gmax,q0(3),xx
      double precision etab(nermx),rtab(nermx)
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
C     Shared variables
      integer i,ib1,ib2,j,nenv1,nenv2,offh1,offh2
      integer iprt(n0,nkap0,nermx),ipet(n0,nkap0,nermx)
      integer ig1,i1,ig2,i2,igx(3),igx1,igx2,igx3
      complex(8),allocatable :: c1(:,:),c2(:,:),wk2(:)
      real(8),allocatable :: g(:),yl(:),g2(:),he(:),hr(:)

C      integer nenv,iblock,mxbldim,nenvb,nblock,nblocklist
C      integer,allocatable :: blocklist(:)
C$    integer mp_numthreads,mp_my_threadnum
C#ifdefC MPI
C      double precision, dimension(:), allocatable :: buffer
C      integer, dimension(:,:), allocatable :: index
C      integer iloop
C      logical cmdopt
C      character*120 strn
C      complex(8), allocatable :: hl(:),hbuf(:)
C#endif

      call tcn('hsibl')

C#ifdefC MPI
C      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
C      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
C      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
C      call strcop(shortname(procid),name,10,'.',i)
C      namelen(procid) = i-1
C      master = 0
C      mlog = cmdopt('--mlog',6,0,strn)
C#endif

C ... Setup
C     print *, '!!'; h = 0
      nbas  = nglob('nbas')
      nspec = nglob('nspec')
      nlmto = ndimh - napw
      if (nspec > nermx) call rx('hsibl: increase nermx')
      alat = s_lat%alat
      plat = s_lat%plat
      qlat = s_lat%qlat
      gmax = s_lat%gmax
      ngabc = s_lat%nabc
      ng = s_lat%ng
      vol = s_lat%vol

C --- <MTO | V | MTO and < MTO | V PW> parts of h ---
      if (nlmto > 0) then

C ... Setup for q-dependent gv ... also makes kv, gv+q and iv
      call pshpr(iprint()-30)
      call sugvec(s_lat,5009,q,0d0,gmax,ngabc,0,ng)
      call poppr

C ... Tables of energies, rsm, indices to them
      call tbhsi(s_spec,nspec,nermx,net,etab,ipet,nrt,rtab,iprt,ltop)
C     ndimx = maximum hamiltonian dimension for any site
      ndimx = nnrl(10,1,nbas,iprmb,nlmto)

C --- Allocate and occupy arrays for YL, energy factors, rsm factors ---
      nlmtop = (ltop+1)**2
      allocate(g(ng*3),yl(ng*nlmtop),g2(ng),he(ng*net),hr(ng*nrt))
C     gv has q already added ...
      call dpzero(q0,3)
      call hsibq1(net,etab,nrt,rtab,ltop,alat,q0,ng,s_lat%gv,g,
     .  g2,yl,he,hr)

C --- Loop over orbitals on first site ---
      allocate(c1(ng,ndimx))
C$DOACROSS LOCAL(ib1,ib2,
C$&              c1,c2,nenv1,nenv2,wk,
C$&       SHARED(nbas,n1,n2,n3,k1,k2,k3,ng,vol,ndimh)
C$&       MP_SCHEDTYPE=RUNTIME
C#ifdefC MPI
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_START_HSIBL,procid,"hsibl")
C#endifC
CC     call defcc(ohl, -ndimh*ndimh)
C      allocate(hl(ndimh*ndimh)); call dpzero(hl,2*ndimh*ndimh)
C      allocate (index(0:numprocs-1,0:nbas-1), stat=ierr)
C      call dstrbp(nbas,numprocs,-1,index(0,0))
C      do  iloop = 1, index(procid,0)
C        ib1 = index(procid,iloop)
C        if (mlog) then
C          call gettime(datim)
C          call awrit4(' hsibl '//datim//' Process %i of %i on '
C     .        //shortname(procid)(1:namelen(procid))//
C     .        ' starting atom %i of %i',' ',256,lgunit(3),
C     .        procid,numprocs,ib1,index(procid,0))
C        endif
C#else
      do  ib1 = 1, nbas
C#endif

        call hsibq2(s_spec,s_site,ib1,1,0,xx,nlmto,iprmb,nenv1,offh1,
     .    ncuti)
        call hsibq4(s_spec,s_site,vol,ib1,1,ng,ngabc,q,qlat,s_lat%igv,
     .    nlmto,iprmb,yl,ipet,iprt,etab,rtab,he,hr,c1)
C       print *, ib1
C       call zprm('c',2,c1,ng,ng,nenv1)

C   ... Multiply potential into wave functions for orbitals in ib1
        call hsibq7(n1,n2,n3,k1,k2,k3,vsm(1,1,1,isp),ng,s_lat%kv,nenv1,
     .    c1)
C       print *, ib1
C       call zprm('cv',2,c1,ng,ng,nenv1)

C   ... Loop over second of (ib1,ib2) site pairs
        allocate(c2(ng,ndimx))
        do  ib2 = ib1, nbas
C     ... PW expansion coffs for LMTOs at this site
          call hsibq2(s_spec,s_site,ib2,1,ng,g2,nlmto,iprmb,
     .      nenv2,offh2,ncuti)
          call hsibq4(s_spec,s_site,vol,ib2,1,ng,ngabc,q,qlat,s_lat%igv,
     .      nlmto,iprmb,yl,ipet,iprt,etab,rtab,he,hr,c2)

C     ... Scalar products phi1*vsm*phi2 for all orbitals in (ib1,ib2)
C         call defcc(owk2,-nenv1*nenv2)
          allocate(wk2(nenv1*nenv2)); call dpzero(wk2,2*nenv1*nenv2)
          call hsibl2(nenv1,nenv2,ng,ncuti,c1,c2,nenv1,0,0,wk2)

C#ifdefC MPI
C          call zmscop(1,nenv1,nenv2,nenv1,ndimh,0,0,offh1,offh2,
C     .      wk2,hl)
C#else
C         print *, ib1,ib2
C         call hsibq5(nenv1,offh1,offh2,nenv2,ndimh,wk2,h)
          call zmscop(1,nenv1,nenv2,nenv1,ndimh,0,0,offh1,offh2,wk2,h)
C#endif

          deallocate(wk2)

C       Ends loop over ib2
        enddo
        deallocate(c2)

C   ... Matrix elements <Hsm| Vsm |PW>
C#ifdefC MPI
C        call hsibl6(nenv1,offh1,ndimh,nlmto,ng,s_lat%igv,napw,igapw,c1,
C     .    hl)
C#else
        call hsibl6(nenv1,offh1,ndimh,nlmto,ng,s_lat%igv,napw,igapw,c1,
     .    h)
C#endif

C     Ends loop over ib1
      enddo
      deallocate(c1)

C      open(55, form='UNFORMATTED')
C      call dpdump(h,ndimh**2*2,-55)
C      call zprm('h',2,h,ndimh,ndimh,ndimh)

C#ifdefC MPI
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_END_HSIBL,procid,"hsibl")
C      ierr = MPE_LOG_EVENT(EVENT_START_BARRIER,procid,"barrier")
C#endifC
C      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_END_BARRIER,procid,"barrier")
C      ierr = MPE_LOG_EVENT(EVENT_START_ALLRED,procid,"allreduce")
C#endifC
CC     call defcc(ohbuf, ndimh*ndimh)
C      allocate(hbuf(ndimh*ndimh))
C      call MPI_ALLREDUCE(hl,hbuf,2*ndimh*ndimh,
C     .     mpi_real8,MPI_SUM,MPI_COMM_WORLD,ierr)
C      if (mlog) then
C        call gettime(datim)
C        call awrit3(' hsibl '//datim//' Process %i of %i on '
C     .        //shortname(procid)(1:namelen(procid))//
C     .        ' allreduce h ndimh=%i',' ',256,lgunit(3),
C     .        procid,numprocs,ndimh)
C      endif
C      call daxpy(2*ndimh*ndimh,1d0,hbuf,1,h,1)
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
      call tcx('hsibl')
C     call zprm('h-i',2,h,ndimh,ndimh,ndimh)
      end

      subroutine hsibl2(n1,n2,ng,ncut2,c1,c2,ndimh,ofh1,ofh2,h)
C- Add scalar product ((phi1*vsm)*phi2) to  h
C ----------------------------------------------------------------------
Ci Inputs
Ci   n1    :number of c1 block of orbitals
Ci   n2    :number of c2 block of orbitals
Ci   ng    :leading dimension of c1,c2
Ci   ncut2 :number of G-vectors
Ci   c1    :product (phi1 vsm) in reciprocal space
Ci   c2    :phi2 in reciprocal space
Ci   ndimh :dimension of h
Ci   ofh1 :row offset to h for first orbital in block
Ci   ofh2 :col offset to h for first orbital in block
Co Outputs
Co   h     :scalar product (phi1 vsm) phi2 is added to h
Cr Remarks
Cr   This routine takes M^2 N operations, and is the most
Cr   time consuming step for moderate to large systems.
Cr   (M = number of orbitals, N = number of mesh points)
Cu Updates
Cu   22 May 00 Adapted from nfp pvhsiq
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer n1,n2,ofh1,ofh2,ncut2(n2),ndimh,ng
      double complex c1(ng,n1),c2(ng,n2),h(ndimh,ndimh)
C ... Local parameters
      integer i,i1,i2,ncut
      double complex csum

      call tcn('hsibl2')

      do  i2 = 1, n2
        ncut = min(ng,ncut2(i2))
        do  i1 = 1, n1
          csum = 0
          do  i = 1, ncut
            csum = csum + dconjg(c1(i,i1))*c2(i,i2)
          enddo
          h(i1+ofh1,i2+ofh2) = h(i1+ofh1,i2+ofh2) + csum
        enddo
      enddo

      call tcx('hsibl2')
      end

      subroutine hsibl6(ndim1,ofh1,ndimh,nlmto,ng,igv,napw,igapw,cv,h)
C- Make matrix elements <Hsm | V | PW>
C ----------------------------------------------------------------------
Ci   ndimh :dimension of hamiltonian
Ci   nlmto :number of lmto's
Ci   norb  :number of orbital types for last ib
Ci   offl  :offl(norb) offset in h to this block of orbitals, last ib
Ci   blks  :blks(iorb) = size of contiguous block of orbitals (gtbsl1)
Ci   ng    :number of G-vectors to represent lmto's
Ci   igv   :list of reciprocal lattice vectors G (from gvlist)
Ci   napw  :number of PW's
Ci   igapw :vector of APWs, in units of reciprocal lattice vectors
Ci   cv    :coffs to PW expansion * V
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
      integer ng,napw,ndim1,ndimh,nlmto,igv(ng,3),igapw(3,napw)
      double complex h(ndimh,ndimh),cv(ng,nlmto)
C ... Local parameters
      integer ig,i2,i2x,ifindiv,ofh1,i1

      do  ig = 1, napw
        i2  = nlmto+ig
C       index matching igv,igapw
        i2x = ifindiv(igapw(1,ig),igv,ng)
        do  i1 = 1, ndim1
          h(ofh1+i1,i2) = h(ofh1+i1,i2) + dconjg(cv(i2x,i1))
        enddo
      enddo
      end
