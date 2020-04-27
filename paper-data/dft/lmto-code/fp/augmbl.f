C#define BLAS
      subroutine augmbl(mode,s_site,s_spec,s_lat,s_ham,isp,lcplxp,q,i0,
     .  ndimh,napw,igapw,h,hso,s)
C- Adds augmentation part of H and S
C ----------------------------------------------------------------------
Cio Structures
Cio  s_site :struct for site-specific data; see structures.h
Ci     Elts read:  spec pos sighh sighk sigkk pihh pihk pikk sighhx
Ci                 sighkx sigkkx pihhx pihkx pikkx sohh sohk sokk
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  bstrux
Cio  s_spec :struct for species-specific data; see structures.h
Ci     Elts read:  lmxa lmxb kmxt rsma pz name orbp
Co     Stored:     orbp
Co     Allocated:  *
Cio    Elts passed:*
Cio    Passed to:  bstrux uspecb
Cio  s_lat  :struct containing lattice information; see structures.h
Ci     Elts read:  alat qlat vol plat awald tol nkd nkq
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:cg indxcg jcg cy qlv dlv
Cio    Passed to:  bstrux hxpbl ghibl hklbl gklbl hxpgbl ghigbl hklgbl
Cio  s_ham  :struct for parameters defining hamiltonian; see structures.h
Ci     Elts read:  *
Co     Stored:     *
Co     Allocated:  *
Cio    Elts passed:lrsa iprmb
Cio    Passed to:  *
Ci Inputs
Ci   mode  :1s digit
Ci         :0 Compute both hamiltonian and overlap
Ci         :1 Compute overlap only
Ci         :  In this case, vavg is not used
Ci         :2 Compute hamiltonian only
Ci         :3 Compute SO contribution to hamiltonian:
Ci         :  LzSz in h and L+S- in hso.
Ci         :  You must also set 10's digit mode to 1.
Ci         :  Call with isp=1 to get 11 block and part of 12 block
Ci         :  Call with isp=2 to get 22 block and remainder of 12 block
Ci         :  Overlap is not touched
Ci         :10s digit
Ci         :  0 do not compute hso
Ci         :  1 compute hso.
Ci         :    Note: only a portion of hso is computed for a
Ci         :    particular isp.  The total hso is assembled
Ci         :    after isp loops from 1..2.  hso should not be
Ci         :    initialized between isp=1 and isp=2 loops.
Ci         :100s digit
Ci         :  0 augmentation matrices from sig,ppi (tau is contained in pi)
Ci         :  1 augmentation matrices from sigx,ppix (tau is contained in pi)
Ci         :  2 augmentation matrices from sigx,ppix-tau (just potential)
Ci         :  3 augmentation matrices from sigx,ppix-tau (just tau)
Ci   isp   :current spin channel
Ci   lcplxp:0 if potential ppi is real; 1 if ppi is complex
Ci   q     :Bloch wave number
Ci   i0    :<=0     => contributions to h and s are accumulated by
Ci         :           summing over all sites
Ci         :1..nbas => contributions to h and s from site i0 only
Ci         :>nbas   => fatal error
Ci   ndimh :dimension of h and s
Ci   napw  :number of PWs in APW part of basis
Ci   igapw :APWs in units of reciprocal lattice vectors
Co Outputs
Co   h     :augmentation part of hamiltonian matrix added to h
Co   hso   :(1,2) spin block of spin-orbit hamiltonian
Co   s     :augmentation part of overlap matrix added to s
Cl Local variables
Cl   nkaph :number of orbital types for a given L quantum no. in basis
Cl         :at augmentation site ia, including local orbitals
Cl   nlmto :number of lmto basis functions
Cr Remarks
Cr   Some expressions labelled JMP refer to J.Math.Phys39, 3393 (1998)
Cb Bugs
Cb   Not really a bug, but an inefficiency:
Cb   Right now, strux are kept for all orbitals in the basis, including
Cb   expansions coffs for local orbitals (which are set to zero).
Cb   Better to condense strux to reduce computational effort for 2-
Cb   and 3-center terms.
Cm MPI
Cm   See remarks in hsibl. Buffers for h and s are taken from the heap.
Cm   In addition a buffer the same size as as h and s for ALLREDUCE.
Cu Updates
Cu   14 Apr 15 New 100s digit mode: can return <pi-tau> only or <tau> only
Cu   05 Jun 14 (D.Pashov) bug fix, MPI case with complex potential
Cu   09 Aug 13 Option to generate SO contribution to KE for one site
Cu   10 Nov 11 Begin migration to f90 structures
Cu   05 Jul 08 (T. Kotani) output density for new PW part
Cu             Option to accumulate energy-weighted output density
Cu   08 Sep 06 (WRL) updated MPI to work with SO coupling
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   03 Feb 05 (A. Chantis) calculate hso
Cu    1 Sep 04 Adapted to handle complex ppi.  S.O. folded into ppi
Cu   25 Aug 04 Adapted to extended local orbitals
Cu   29 Jun 04 (A. Chantis) Include LzSz spin-orbit coupling
Cu   14 Aug 02 Added overlap-only option
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   11 Jan 02 Adapted for f90 constructs
Cu   27 Aug 01 Extended to local orbitals.
Cu   17 Jun 00 spin polarized
Cu   18 May 00 Adapted from nfp augm_q.f
Cu   1998      (DLN) parallel version for SGI
C ----------------------------------------------------------------------
      use mpi
      use structures
      implicit none
C ... Passed parameters
C#ifdefC MPE
C      include "mpef.h"
C#endif
      integer procid, master, numprocs, ierr
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
C     character*20 ext
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
C     double precision starttime, endtime
      logical mlog,cmdopt
      integer lgunit
      character*120 strn

      integer mode,lcplxp,isp,i0,ndimh,napw,igapw(3,napw)
      double precision q(3)
      double complex h(ndimh,ndimh),s(ndimh,ndimh),hso(ndimh,ndimh)
C ... For structures
!      include 'structures.h'
      type(str_site)::  s_site(*)
      type(str_spec)::  s_spec(*)
      type(str_lat)::   s_lat
      type(str_ham)::   s_ham
C ... Dynamically allocated local arrays
      complex(8),allocatable:: b(:)
      real(8), pointer :: ppikk(:,:),ppihk(:,:),ppihh(:,:)
      real(8), pointer :: sigkk(:,:),sighk(:,:),sighh(:,:)
C     real(8), pointer :: taukk(:,:),tauhk(:,:),tauhh(:,:)
C ... Local parameters
      integer nlmbx,nlmax,ktop0,lofb,mode0,mode0x
      parameter (ktop0=20, nlmbx=49, nlmax=49, lofb=(ktop0+1)*nlmax)
      double complex g(lofb)
      integer ia,i1,i2,isa,k,kmax,lmxa,lmxha,nbas,nglob,nlma,
     .  nlmha,nkaph,mode1,mode2,nlmto,nsp,nspc,ns4,nelt2
      double precision rsma,pa(3),xx,alat,qlat(3,3),vol
      integer, parameter :: n0=10, nkap0=4
      integer lh(nkap0),lxa(0:ktop0),nkapi
      double precision eh(n0,nkap0),rsmh(n0,nkap0)

      integer, dimension(:), allocatable :: bproc
C ... Dynamically allocated local arrays
      complex(8), allocatable :: hl(:)
      complex(8), allocatable :: sl(:)
      complex(8), allocatable :: hsol(:)
      complex(8), allocatable :: sbuf(:)
      complex(8), allocatable :: hbuf(:)

C#ifdefC MPE
CC Event numbers:
C      include "events.ins"
C#endif

      call tcn ('augmbl')

      call mpi_comm_size( mpi_comm_world, numprocs, ierr )
! numprocs = 1 means we are either in serial or MPIK mode. Lets assume MPIK for now.. to be ammended later
      numprocs = 1

      if (numprocs > 1) then
        call mpi_comm_rank( mpi_comm_world, procid, ierr )

        call mpi_get_processor_name(name, resultlen, ierr)
        call strcop(shortname(procid),name,10,'.',k)
        namelen(procid) = k-1
        master = 0
        mlog = cmdopt('--mlog',6,0,strn)
      end if

C --- Setup ---
      nbas  = nglob('nbas')
      nkaph = nglob('nkaph')
      nsp   = nglob('nsp')
      nspc = nglob('nspc')
      ns4 = nsp*nspc; if (mod(s_ham%lrsa,10) == 0) ns4 = nsp ! See "Local variables"
      mode0 = mod(mode,10)
      mode0x = mode0 ; if (mode0 == 3) mode0=2
      mode1 = mod(mode/10,10)
      mode2 = mod(mode/100,10)
      nlmto = ndimh-napw
      alat = s_lat%alat; qlat = s_lat%qlat; vol = s_lat%vol
C     tpiba = 2d0*4d0*datan(1d0)/alat
      call rxx(i0 > nbas,'augmbl: nonsensical i0')
      call rxx(mode0x == 3.and.mode1 /= 1,'augmbl: nonsensical mode')

      allocate (b(lofb*ndimh))

C --- Loop over augmentation sites --- ---
      if (numprocs > 1) then
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_START_AUGMBL,procid,"augmbl")
C#endif
        allocate(hl(ndimh*ndimh)); call dpzero(hl,2*ndimh*ndimh)
        allocate(sl(ndimh*ndimh)); call dpzero(sl,2*ndimh*ndimh)
        if (mode1 == 1 .and. lcplxp /= 0) then
          allocate(hsol(ndimh*ndimh)); call dpzero(hsol,2*ndimh*ndimh)
        else
          allocate(hsol(1))
        endif
        allocate (bproc(0:numprocs))
        call dstrbp(nbas,numprocs,1,bproc(0))
        i1 = bproc(procid)
        i2 = bproc(procid+1)-1
      else

        i1 = 1; i2 = nbas
        if (i0 > 0) then; i1 = i0; i2 = i0; endif
      end if
      do  ia = i1, i2
        if (numprocs > 1) then
          if (mlog .and. ia == bproc(procid)) then
            call gettime(datim)
            call awrit4(' augmbl '//datim//' Process %i of %i on '
     .          //shortname(procid)(1:namelen(procid))//
     .          ' starting atoms %i to %i',' ',256,lgunit(3),
     .          procid,numprocs,bproc(procid),bproc(procid+1)-1)
          endif
        end if
        isa = s_site(ia)%spec; pa = s_site(ia)%pos
        lmxa = s_spec(isa)%lmxa; lmxha = s_spec(isa)%lmxb
        if (lmxa == -1) cycle
        kmax = s_spec(isa)%kmxt; rsma = s_spec(isa)%rsma
        nlmha = (lmxha+1)**2; nlma  = (lmxa+1)**2

C       tauhh => s_site(ia)%tauhh
C       tauhk => s_site(ia)%tauhk
C       taukk => s_site(ia)%taukk
        if (mode2 == 0) then    ! Use ppi, etc
          sighh => s_site(ia)%sighh; sighk => s_site(ia)%sighk; sigkk => s_site(ia)%sigkk
          ppihh => s_site(ia)%pihh; ppihk => s_site(ia)%pihk; ppikk => s_site(ia)%pikk
        elseif (mode2 == 1) then ! Use ppix, etc
          sighh => s_site(ia)%sighhx; sighk => s_site(ia)%sighkx; sigkk => s_site(ia)%sigkkx
          ppihh => s_site(ia)%pihhx; ppihk => s_site(ia)%pihkx; ppikk => s_site(ia)%pikkx
        elseif (mode2 == 2 .or. mode2 == 3) then ! Use ppix, or K.E. only
          if (lcplxp /= 0) call rx('augmbl not ready for lcplxp')
          sighh => s_site(ia)%sighh; sighk => s_site(ia)%sighk; sigkk => s_site(ia)%sigkk
          call uspecb(0,1,s_spec,isa,isa,lh,rsmh,eh,nkapi)
          do  k = 0, kmax
            lxa(k) = lmxa
          enddo

          nelt2 = nkaph*nkaph*nlmha*nlmha*nsp*nspc
          allocate(ppihh(1,nelt2)) ! ppihh(nkaph,nkaph,nlmha,nlmha,nsp)
          if (mode2 == 2) call dcopy(1*nelt2,s_site(ia)%pihhx,1,ppihh,1)
          if (mode2 == 3) call dpzero(ppihh,1*nelt2)
          call taupitoh(mode2,nkaph,nkaph,nlmha,nlmha,lh,lh,lmxha,isp,ppihh,s_site(ia)%tauhh)

          nelt2 = (kmax+1)*(kmax+1)*nlma*nlma*nsp*nspc
          allocate(ppikk(1,nelt2)) ! ppikk(nkaph,nkaph,nlmha,nlmha,nsp)
          if (mode2 == 2) call dcopy(1*nelt2,s_site(ia)%pikkx,1,ppikk,1)
          if (mode2 == 3) call dpzero(ppikk,1*nelt2)
          call taupitoh(mode2,nkaph,kmax+1,nlmha,nlma,lh,lxa,lmxha,isp,ppikk,s_site(ia)%taukk)

          nelt2 = nkaph*(kmax+1)*nlmha*nlmha*nsp*nspc
          allocate(ppihk(1,nelt2)) ! ppihk(nkaph,0:kmax,nlmha,nlma,nsp)
          if (mode2 == 2) call dcopy(1*nelt2,s_site(ia)%pihkx,1,ppihk,1)
          if (mode2 == 3) call dpzero(ppihk,1*nelt2)
          call taupitoh(mode2,kmax+1,kmax+1,nlma,nlma,lxa,lxa,lmxa,isp,ppihk,s_site(ia)%tauhk)

        else
          call rxi('augmbl : bad mode',mode)
        endif

C        if (nspc == 2) then
C          eula = s_site(ia)%eula
C          call rotspv(144,eula,1,1,1,routr,rout,rinr,rinc)
C          k = size(sighh)
C          call rotspv(444,eula,k,k,k,sighh,xx,xx,sighhz)
C          k = size(sighk)
C          call rotspv(444,eula,k,k,k,sighk,xx,xx,sighkz)
C          k = size(sigkk)
C          call rotspv(444,eula,k,k,k,sigkk,xx,xx,sigkkz)
C
C          k = size(ppihh)
C          call rotspv(444,eula,k,k,k,ppihh,xx,xx,ppihhz)
C          k = size(ppihk)
C          call rotspv(444,eula,k,k,k,ppihk,xx,xx,ppihkz)
C          k = size(ppikk)
C          call rotspv(444,eula,k,k,k,ppikk,xx,xx,ppikkz)
C        endif

C   --- Make strux to expand all orbitals at site ia ---
        call rxx((kmax+1)*nlma > lofb,'augmbl: increase lofb')
        call bstrux(0,s_lat,s_site,s_spec,s_lat%cg,s_lat%indxcg,
     .    s_lat%jcg,s_lat%cy,s_ham%iprmb,nbas,ia,pa,rsma,q,kmax,nlma,
     .    ndimh,napw,igapw,b,xx)

C   --- Add 1-center and 2-center terms ---
        if (numprocs > 1) then
          if (lcplxp == 0 .or. mode1 >= 2) then
            if (mode0 > 1) call rx('mode0>1 not implemented')
            call augq12(mode0,ia,isp,nkaph,s_ham%iprmb,lmxha,nlmha,kmax,
     .        nlma,sighh,ppihh,sighk,ppihk,b,ndimh,nlmto,sl,hl)
          else
            call augq2z(mode0,mode1,ia,isp,nkaph,s_ham%iprmb,lmxha,nlmha,
     .        kmax,nlma,sighh,ppihh,sighk,ppihk,b,ndimh,nlmto,sl,hl,hsol)
          endif
        else
C#ifndef ALL3C
          if (ns4 == 4) then
            call rx('augq2z -> noncoll')
          else if (lcplxp == 0 .or. mode1 >= 2) then
            if (mode0 > 1) call rx('mode0>1 not implemented')
            call augq12(mode0,ia,isp,nkaph,s_ham%iprmb,lmxha,nlmha,kmax,
     .        nlma,sighh,ppihh,sighk,ppihk,b,ndimh,nlmto,s,h)
          else
            if (mode0x == 3) then
            ppihh => s_site(ia)%sohh
            ppihk => s_site(ia)%sohk
          endif
          call augq2z(mode0,mode1,ia,isp,nkaph,s_ham%iprmb,lmxha,nlmha,
     .      kmax,nlma,sighh,ppihh,sighk,ppihk,b,ndimh,nlmto,s,h,hso)
          endif
C#endif
        endif

C   --- Add B+ sig B to S and B+ ppi B to H ---
        if (numprocs > 1) then
          if (mode0 <= 1) then
            call augqs3(kmax,lmxa,nlma,ndimh,isp,g,sigkk,b,sl)
          endif
          if (mode0 == 0 .and. lcplxp == 0) then
            call augqp3(kmax,nlma,ndimh,isp,g,ppikk,b,hl)
          elseif (mode0 == 0 .and. lcplxp /= 0) then
            call augq3z(mode1,kmax,nlma,ndimh,isp,g,ppikk,b,hl,hsol)
          endif
        else
          if (mode0 <= 1) call augqs3(kmax,lmxa,nlma,ndimh,isp,g,sigkk,b,s)
          if (mode0 /= 1 .and. lcplxp == 0) then
            call augqp3(kmax,nlma,ndimh,isp,g,ppikk,b,h)
          elseif (mode0 /= 1 .and. lcplxp /= 0) then
            if (mode0x == 3) ppikk => s_site(ia)%sokk
C           if (nspc == 2) call rx('augq3z -> noncoll')
            call augq3z(mode1,kmax,nlma,ndimh,isp,g,ppikk,b,h,hso)
          endif
        endif

C ... end loop over ia
      enddo
      if (numprocs > 1) then
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_END_AUGMBL,procid,"augmbl")
C      ierr = MPE_LOG_EVENT(EVENT_START_BARRIER,procid,"barrier")
C#endif
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_END_BARRIER,procid,"barrier")
C      ierr = MPE_LOG_EVENT(EVENT_START_ALLRED,procid,"allreduce")
C#endif
      allocate(hbuf(ndimh*ndimh))
      call MPI_ALLREDUCE(hl,hbuf,2*ndimh*ndimh,
     .     mpi_real8,MPI_SUM,MPI_COMM_WORLD,ierr)
      if (mlog) then
        call gettime(datim)
        call awrit3(' augmbl '//datim//' Process %i of %i on '
     .        //shortname(procid)(1:namelen(procid))//
     .        ' allreduce h ndimh=%i',' ',256,lgunit(3),
     .        procid,numprocs,ndimh)
      endif
      call daxpy(2*ndimh*ndimh,1d0,hbuf,1,h,1)
      deallocate(hbuf)
      allocate(sbuf(ndimh*ndimh))
      call MPI_ALLREDUCE(sl,sbuf,2*ndimh*ndimh,
     .     mpi_real8,MPI_SUM,MPI_COMM_WORLD,ierr)
      if (mlog) then
        call gettime(datim)
        call awrit3(' augmbl '//datim//' Process %i of %i on '
     .        //shortname(procid)(1:namelen(procid))//
     .        ' allreduce s ndimh=%i',' ',256,lgunit(3),
     .        procid,numprocs,ndimh)
      endif
      call daxpy(2*ndimh*ndimh,1d0,sbuf,1,s,1)
      deallocate(sbuf)
      if (mode1 == 1 .and. lcplxp /= 0) then
       allocate(sbuf(ndimh*ndimh))
       call MPI_ALLREDUCE(hsol,sbuf,2*ndimh*ndimh,
     .     mpi_real8,MPI_SUM,MPI_COMM_WORLD,ierr)
       if (mlog) then
        call gettime(datim)
        call awrit3(' augmbl '//datim//' Process %i of %i on '
     .        //shortname(procid)(1:namelen(procid))//
     .        ' allreduce hso ndimh=%i',' ',256,lgunit(3),
     .        procid,numprocs,ndimh)
       endif
       call daxpy(2*ndimh*ndimh,1d0,sbuf,1,hso,1)
       deallocate(sbuf)
       deallocate(hsol)
      endif
      deallocate(hl,sl)
      deallocate(bproc)
C#ifdefC MPE
C      ierr = MPE_LOG_EVENT(EVENT_END_ALLRED,procid,"allreduce")
C#endif
      endif

C      call z2herm('U',ndimh,ndimh,h)
C      call z2herm('U',ndimh,ndimh,s)
C      call zprm('h-aug',2,h,ndimh,ndimh,ndimh)
C      call zprm('s-aug',2,s,ndimh,ndimh,ndimh)

      deallocate (b)

      call tcx ('augmbl')

      end

      subroutine augq12(mode,ia,isp,nkaph,iprmb,lmxha,nlmha,kmax,
     .  nlma,sighh,ppihh,sighp,ppihp,b,ndimh,nlmto,s,h)
C- Add one and two-center terms to hamiltonian and overlap matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 compute both hamiltonian and overlap
Ci         :  otherwise, compute overlap only.
Ci         :  In this case, vavg is not used
Ci   ia    :augmentation site about which strux are expanded
Ci   isp   :current spin channel
Ci   nkaph :dimensions augmentation matrices
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   nlmha :dimensions augmentation potential matrix at site a
Ci   lmxha :dimensions sighh at site a
Ci   kmax  :polynomial cutoff
Ci   nlma  :augmentation L-cutoff
Ci   sighh :augmentation head-head overlap matrix
Ci   ppihh :augmentation head-head potential matrix
Ci   sighp :augmentation head-Pkl overlap matrix
Ci   ppihp :augmentation head-Pkl potential matrix
Ci   b     :Bloch strux connecting site ia to all sites
Ci   ndimh :hamiltonian dimension
Ci   nlmto :Number of lmto basis functions
Co Outputs
Co   h     :1- and 2- center augmentation part of ham. added to h
Co   s     :1- and 2- center augmentation part of ovlp added to s
Cr Remarks
Cr  In this implementation, the augmentation matrices and the row
Cr  dimension of the structure constants b follow normal L order.
Cr  The column dimension of b is permuted in iprmb order.
Cu Updates
Cu   01 Sep 04 folded so into complex potential
Cu   29 Jun 04 (A. Chantis) added 1- and 2- center spherical so*Lz*Sz
Cu   14 Aug 02 Added overlap-only option
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,ia,isp,kmax,nkaph,ndimh,nlma,lmxha,nlmha,iprmb(*)
      integer nlmto
      double precision
     .  sighh(nkaph,nkaph,0:lmxha,1), ppihh(nkaph,nkaph,nlmha,nlmha,1),
     .  sighp(nkaph,0:kmax,0:lmxha,1),ppihp(nkaph,0:kmax,nlmha,nlma,1)
      double complex b(0:kmax,nlma,ndimh),s(ndimh,ndimh),h(ndimh,ndimh)
C ... Local parameters
      integer iorb,ik1,j,k,ilma,i1,i2,ilm1,ilm2,l1,n0,nkap0,jorb,ik2,l2
      parameter (n0=10,nkap0=4)
      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),norb
      double precision xx
      double complex cadd

C     call zprm('strux',2,b,(kmax+1)*nlma,(kmax+1)*nlma,ndimh)
C     call zprm('strux(k=0)',2,b(0,:,:),nlma,nlma,ndimh)

C --- Loop over basis functions at site ia (augentation index) ---
      call orbl(ia,0,nlmto,iprmb,norb,ltab,ktab,xx,offl,xx)
      do  iorb = 1, norb
C       l1,ik1 = l and kaph indices, needed for sigma
        l1  = ltab(iorb)
        ik1 = ktab(iorb)
C       i1 = orbital index in iprmb order; ilm1 = augm. index in L order
        i1 = offl(iorb)
        do  ilm1 = l1**2+1, (l1+1)**2
          i1 = i1+1

C     ... Two-center terms
C         Loop over basis functions 1..ndimh from all sites
          if (mode == 0) then
          do  j = 1, ndimh
            do  k = 0, kmax
              cadd = sighp(ik1,k,l1,isp)*b(k,ilm1,j)
              s(i1,j) = s(i1,j) + cadd
              s(j,i1) = s(j,i1) + dconjg(cadd)
              do  ilma = 1, nlma
                cadd = ppihp(ik1,k,ilm1,ilma,isp)*b(k,ilma,j)
                h(i1,j) = h(i1,j) + cadd
                h(j,i1) = h(j,i1) + dconjg(cadd)
              enddo
            enddo
          enddo

C     ... One-center terms
          do  jorb = 1, norb
            l2  = ltab(jorb)
            ik2 = ktab(jorb)
            i2 = offl(jorb)
            do  ilm2 = l2**2+1, (l2+1)**2
              i2 = i2+1
              h(i1,i2) = h(i1,i2) + ppihh(ik1,ik2,ilm1,ilm2,isp)
              if (ilm1 == ilm2) s(i1,i2) = s(i1,i2) + sighh(ik1,ik2,l1,isp)
            enddo
          enddo
          else

          do  j = 1, ndimh
            do  k = 0, kmax
              cadd = sighp(ik1,k,l1,isp)*b(k,ilm1,j)
              s(i1,j) = s(i1,j) + cadd
              s(j,i1) = s(j,i1) + dconjg(cadd)
            enddo
          enddo

C     ... One-center terms
          do  jorb = 1, norb
            l2  = ltab(jorb)
            ik2 = ktab(jorb)
            i2 = offl(jorb)
            do  ilm2 = l2**2+1, (l2+1)**2
              i2 = i2+1
              if (ilm1 == ilm2) s(i1,i2) = s(i1,i2) + sighh(ik1,ik2,l1,isp)
            enddo
          enddo
          endif

        enddo
      enddo

      end

      subroutine augq2z(mode,mode1,ia,isp,nkaph,iprmb,lmxha,nlmha,kmax,
     .  nlma,sighh,ppihh,sighp,ppihp,b,ndimh,nlmto,s,h,hso)
C- Add one and two-center terms to h,s for complex potential
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 compute both hamiltonian and overlap
Ci         :1 compute overlap only
Ci         :  In this case, vavg is not used
Ci         :2 compute hamiltonian only
Ci   mode1 :0 do not compute hso
Ci         :1 compute hso = (1,2) spin block of h
Ci   ia    :augmentation site about which strux are expanded
Ci   isp   :current channel for spin-diagonal part of h
Ci   nkaph :dimensions augmentation matrices
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   nlmha :dimensions augmentation potential matrix at site a
Ci   lmxha :dimensions sighh at site a
Ci   kmax  :polynomial cutoff
Ci   nlma  :augmentation L-cutoff
Ci   sighh :augmentation head-head overlap matrix
Ci   ppihh :augmentation head-head potential matrix
Ci         :ppihh(:,:,:,:,1:2) = spin-diagonal blocks
Ci         :ppihh(:,:,:,:,3)   = spin (1,2) block of ppihh
Ci   sighp :augmentation head-Pkl overlap matrix
Ci   ppihp :augmentation head-Pkl potential matrix
Ci         :ppihp(:,:,:,:,1:2) = spin-diagonal blocks
Ci         :ppihp(:,:,:,:,3)   = spin (1,2) block
Ci         :ppihp(:,:,:,:,4)   = spin (2,1) block
Ci   b     :Bloch strux connecting site ia to all sites
Ci   ndimh :hamiltonian dimension
Co Outputs
Co   h     :1- and 2- center augmentation part of ham. added to h
Co   s     :1- and 2- center augmentation part of ovlp added to s
Co   hso   :1- and 2- center spin augmentation part of spin (1,2) block
Cr Remarks
Cr  In this implementation, the augmentation matrices and the row
Cr  dimension of the structure constants b follow normal L order.
Cr  The column dimension of b is permuted in iprmb order.
Cr
Cr  For each of the isp={1,2}, one of the following two terms for 2C part
Cr  of hso(1,2) is assembled. The 2C term has the form
Cr    h_ij = Sum_kL [(p_{i;kL;12}*b_{kL;j}) + conjg(p_{j;kL;21}*b_{kL;i})]
Cr  There are two kinds two-center terms which are symbolically:
Cr    hso^1_{ij;12} =  p_{ij,12}*b_{j} augmentation at i, head at j
Cr    hso^2_{ji;21} =  p_{ji,21}*b_{i} augmentation at j, head at i
Cr  where 1 = spin-up and 2 = spin-down.
Cr  Each term has a corresponding part in the other spin block.
Cr  Only hso(:;12) is assembled, both terms must be included there.
Cr  Since hso^2_{ji;21} = conjg(hso^2_{ij;12}), we can write:
Cr    hso_{ji;12} =  p_{ij,12}*b_{j} + conjg(p_{ji,21}*b_{i})
Cr  the first term is added to hso when isp=1; the second when isp=2.
Cr
Cr  If the structure constants become noncollinear, additional terms have
Cr  to be added in the matrix element above.
Cu Updates
Cu   08 Aug 13 Option to compute hamiltonian only
Cu   03 Feb 05 (A. Chantis) added 1- and 2- center spherical so*(LxSx+LySy)
Cu   01 Sep 04 folded so into complex potential
Cu   29 Jun 04 (A. Chantis) added 1- and 2- center spherical so*Lz*Sz
Cu   14 Aug 02 Added overlap-only option
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,mode1,ia,isp,kmax,nkaph,nlma,lmxha,nlmha,iprmb(*),ndimh,nlmto
      double precision sighh(nkaph,nkaph,0:lmxha,1),sighp(nkaph,0:kmax,0:lmxha,1)
      double complex ppihh(nkaph,nkaph,nlmha,nlmha,isp+2*mode1),
     .               ppihp(nkaph,0:kmax,nlmha,nlma,isp+2*mode1)
      double complex b(0:kmax,nlma,ndimh),s(ndimh,ndimh),
     .               h(ndimh,ndimh),hso(ndimh,ndimh)
C ... Dynamically allocated local arrays
      complex(8),allocatable:: tso(:,:,:,:)
C ... Local parameters
      integer iorb,ik1,j,k,ilma,i1,i2,ilm1,ilm2,l1,n0,nkap0,jorb,ik2,l2,jsp,ksp
      parameter (n0=10,nkap0=4)
      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),norb
      double precision xx
      double complex cadd,cadd1

      if (mode /= 1 .and. mode1 == 1) then
        allocate(tso(ndimh,ndimh,2,2))
        call dpzero(tso,ndimh*ndimh*4*2)
      endif

C     call zprm('strux',2,b,(kmax+1)*nlma,(kmax+1)*nlma,ndimh)
C --- Loop over basis functions at site ia (augentation index) ---
      call orbl(ia,0,nlmto,iprmb,norb,ltab,ktab,xx,offl,xx)
      do  iorb = 1, norb
C       l1,ik1 = l and kaph indices, needed for sigma
        l1  = ltab(iorb)
        ik1 = ktab(iorb)
C       i1 = orbital index in iprmb order; ilm1 = augm. index in L order
        i1 = offl(iorb)
        do  ilm1 = l1**2+1, (l1+1)**2
          i1 = i1+1

C     ... Two-center terms, head at j, augmentation at i.  Use (i,j)=congj(j,i)
C         Loop over basis functions 1..ndimh from all sites
C         Compute hamiltonian
          if (mode /= 1) then
            do  j = 1, ndimh
            do  k = 0, kmax
              do  ilma = 1, nlma
                cadd = ppihp(ik1,k,ilm1,ilma,isp)*b(k,ilma,j)
C           ... Make 2C term (1,2) block of hso;  see Remarks
C               For SO this potential is the LxSx+LySy part.
                if (mode1 == 1) then
                  do  jsp = 1, 2
                    if (isp /= jsp) then
                      ksp = 2*isp + jsp - 1
                      cadd1 = ppihp(ik1,k,ilm1,ilma,ksp)*b(k,ilma,j)
C                     Add this term to hso when isp=1
                      tso(i1,j,isp,jsp) = tso(i1,j,isp,jsp)
     .                                  + 0.5d0*cadd1
C                     Add this term to hso when isp=2
                      tso(j,i1,jsp,isp) = tso(j,i1,jsp,isp)
     .                                  + 0.5d0*dconjg(cadd1)
                    endif
                  enddo
                endif
                h(i1,j) = h(i1,j) + cadd
                h(j,i1) = h(j,i1) + dconjg(cadd)
              enddo
            enddo
            enddo

C       ... One-center terms
            do  jorb = 1, norb
              l2  = ltab(jorb)
              ik2 = ktab(jorb)
              i2 = offl(jorb)
              do  ilm2 = l2**2+1, (l2+1)**2
                i2 = i2+1
                h(i1,i2) = h(i1,i2) + ppihh(ik1,ik2,ilm1,ilm2,isp)
C        ...  Make 1c LxSx+LySy part of SO
                if (mode1 == 1 .and. isp == 2) hso(i1,i2) =
     .            hso(i1,i2) + 0.5d0*ppihh(ik1,ik2,ilm1,ilm2,3)
              enddo
            enddo

          endif

C         Compute overlap
          if (mode <= 1) then

            do  j = 1, ndimh
              do  k = 0, kmax
                cadd = sighp(ik1,k,l1,isp)*b(k,ilm1,j)
                s(i1,j) = s(i1,j) + cadd
                s(j,i1) = s(j,i1) + dconjg(cadd)
              enddo
            enddo

C       ... One-center terms
            do  jorb = 1, norb
              l2  = ltab(jorb)
              ik2 = ktab(jorb)
              i2 = offl(jorb)
              do  ilm2 = l2**2+1, (l2+1)**2
                i2 = i2+1
                if (ilm1 == ilm2)
     .            s(i1,i2) = s(i1,i2) + sighh(ik1,ik2,l1,isp)
              enddo
            enddo
          endif

        enddo
      enddo

C ... Add tso into hso
      if (mode /= 1 .and. mode1 == 1) then
        call dpadd(hso(1,1),tso(1,1,1,2),1,2*ndimh*ndimh,1d0)
        deallocate (tso)
      endif

      end

      subroutine augqs3(kmax,lmxa,nlma,ndimh,isp,g,sig,b,s)
C- Add B+ sig B to s for L-diagonal sig
C ----------------------------------------------------------------------
Ci Inputs
Ci   kmax  :polynomial cutoff
Ci   lmxa  :dimensions sig at site a
Ci   nlma  :augmentation L-cutoff
Ci   ndimh :hamiltonian dimension
Ci   isp   :current spin channel
Ci   g     :complex work array of dimension (kmax+1)*nlma
Ci   sig   :augmentation Pkl-Pkl overlap matrix
Ci   b     :Bloch structure constants (hxpbl)
Co Outputs
Co   s     :overlap matrix
Cr Remarks
Cu Updates
Cu   02 Jul 15 Remove fixed dimensioning nlmax; call ll routine
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer kmax,lmxa,nlma,ndimh,isp
      double precision sig(0:kmax,0:kmax,0:lmxa,isp)
      double complex b(0:kmax,nlma,ndimh),s(ndimh,ndimh),g(0:kmax,nlma),csum
C ... Local parameters
      integer ll
C     integer nlmax
C     parameter (nlmax=49)
C     integer kjlm
      integer kjtop,i1,i2,ilm,k1,k2,l !,lla(nlma)
C     integer kjlm
      double complex zdotc
C     data lla/0,3*1,5*2,7*3,9*4,11*5,13*6/

C     if (nlma > nlmax) call rxi('augqs3: increase nlmax to',nlma)

C     call tcn('augqs3')
      kjtop = nlma*(kmax+1)
      do  i2 = 1, ndimh

C   ... Make sig*b in g
        do  ilm = 1, nlma
          l = ll(ilm)
          do  k1 = 0, kmax
            g(k1,ilm) = 0d0
            do  k2 = 0, kmax
              g(k1,ilm) = g(k1,ilm) + sig(k1,k2,l,isp)*b(k2,ilm,i2)
            enddo
          enddo
        enddo
C   ... Make dot products with vectors i1
        do  i1 = 1, i2
C#ifdef BLAS
          csum = zdotc(kjtop,b(0,1,i1),1,g,1)
C#elseC
C          csum = (0d0,0d0)
C          do  kjlm = 0, kjtop-1
C            csum = csum + dconjg(b(kjlm,1,i1))*g(kjlm,1)
C          enddo
C#endif
C          csum = (0d0,0d0)
C          do  ilm = 1,nlma
C            do k2 = 0, kmax
C              csum = csum + dconjg(b(k2,ilm,i1))*g(k2,ilm)
C            enddo
C          enddo
          s(i1,i2) = s(i1,i2) + csum
        enddo
      enddo
C     call tcx('augqs3')

      end

      subroutine augqp3(kmax,nlma,ndimh,isp,g,ppi,b,h)
C- Add B+ ppi B to H for non-L-diagonal matrix ppi
C ----------------------------------------------------------------------
Ci Inputs
Ci   kmax  :polynomial cutoff
Ci   nlma  :augmentation L-cutoff
Ci   ndimh :hamiltonian dimension
Ci   isp   :current spin channel
Ci   g     :complex work array of dimension (kmax+1)*nlma
Ci   ppi   :augmentation Pkl-Pkl potential matrix
Ci   b     :Bloch structure constants (hxpbl)
Co Outputs
Co   h     :3-center from this augmentation site added to h
Cr Remarks
Cu Updates
Cu 01 Sep 04 folded so into complex potential
Cu 29 Jun 04 (A. Chantis) added 3- center so*Sz*Lz (spherical part)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer kmax,nlma,ndimh,isp
      double precision ppi(0:kmax,0:kmax,nlma,nlma,isp)
      double complex b(0:kmax,nlma,ndimh),h(ndimh,ndimh),g(0:kmax,nlma),csum
C ... Local parameters
C     integer kjlm
      integer i1,i2,jlm1,jlm2,k1,k2,kjtop
      double complex zdotc

C     call tcn('augqp3')
      kjtop = nlma*(kmax+1)

      do  i2 = 1, ndimh
C   ... g <- ppi*b
        call dpzero(g,2*kjtop)
        do  jlm1 = 1, nlma
          do  jlm2 = 1, nlma
            do  k2 = 0, kmax
              do  k1 = 0, kmax
                g(k1,jlm1) = g(k1,jlm1) +
     .                       ppi(k1,k2,jlm1,jlm2,isp)*b(k2,jlm2,i2)
              enddo
            enddo
          enddo
        enddo

C   ... Make dot products with vectors i1
        do  i1 = 1, i2
C#ifdef BLAS
          csum = zdotc(kjtop,b(0,1,i1),1,g,1)
C#elseC
C          csum = (0d0,0d0)
C          do  kjlm = 0, kjtop-1
C            csum = csum + dconjg(b(kjlm,1,i1))*g(kjlm,1)
C          enddo
C#endif
C          do  ilm = 1,nlma
C            do k2 = 0, kmax
C              csum = csum + dconjg(b(k2,ilm,i1))*g(k2,ilm)
C            enddo
C          enddo
          h(i1,i2) = h(i1,i2) + csum
        enddo
      enddo
C     call tcx('augqp3')
      end

      subroutine augq3z(mode1,kmax,nlma,ndimh,isp,g,ppi,b,h,hso)
C- Add B+ ppi B to H for non-L-diagonal, complex matrix ppi
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode1 :0 do not compute hso
Ci         :1 compute hso = (1,2) spin block of h
Ci   kmax  :polynomial cutoff
Ci   nlma  :augmentation L-cutoff
Ci   ndimh :hamiltonian dimension
Ci   isp   :current spin channel
Ci   g     :complex work array of dimension (kmax+1)*nlma
Ci   ppi   :augmentation Pkl-Pkl potential matrix
Ci   b     :Bloch structure constants (hxpbl)
Co Outputs
Co   h     :3-center from this augmentation site added to h
Co   hso   :3-center from this augmentation site added to hso
Cr Remarks
Cu Updates
Cu 03 Feb 05 (A. Chantis) added 3- center spherical so*(LxSx+LySy)
Cu 01 Sep 04 folded so into complex potential
Cu 29 Jun 04 (A. Chantis) added 3- center so*Sz*Lz (spherical part)
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode1,kmax,nlma,ndimh,isp
      double complex g(0:kmax,nlma),b(0:kmax,nlma,ndimh),
     .  h(ndimh,ndimh),hso(ndimh,ndimh)
      double complex ppi(0:kmax,0:kmax,nlma,nlma,isp+2*mode1)
C ... Local parameters
      integer i1,i2,jlm1,jlm2,k1,k2,kjtop
C     integer kjlm,ilm
      double complex csum,gso(0:kmax,nlma),csum1
      double complex zdotc

C     call tcn('augqp3')
      kjtop = nlma*(kmax+1)
      do  i2 = 1, ndimh
C   ... g <- ppi*b
        call dpzero(g,2*kjtop)
        call dpzero(gso,2*kjtop)
        do  jlm1 = 1, nlma
          do  jlm2 = 1, nlma
            do  k2 = 0, kmax
              do  k1 = 0, kmax
                g(k1,jlm1) = g(k1,jlm1) +
     .                       ppi(k1,k2,jlm1,jlm2,isp)*b(k2,jlm2,i2)
C               LxSx+LySy part of SO
                if (mode1 == 1 .and. isp == 2)
     .            gso(k1,jlm1) = gso(k1,jlm1) +
     .            ppi(k1,k2,jlm1,jlm2,3)*b(k2,jlm2,i2)
              enddo
            enddo
          enddo
        enddo

C   ... Make dot products with vectors i1
        do  i1 = 1, ndimh
C#ifdef BLAS
          csum = zdotc(kjtop,b(0,1,i1),1,g,1)
C         LxSx+LySy part of SO
          if (mode1 == 1 .and. isp == 2) csum1 = zdotc(kjtop,b(0,1,i1),1,gso,1)
C#elseC
C          csum = (0d0,0d0)
C          if (mode1 == 1 .and. isp == 2) csum1 = (0d0,0d0)
C          do  kjlm = 0, kjtop-1
C            csum = csum + dconjg(b(kjlm,1,i1))*g(kjlm,1)
CC           LxSx+LySy part of SO
C            if (mode1 == 1 .and. isp == 2)
C     .      csum1 = csum1 + dconjg(b(kjlm,1,i1))*gso(kjlm,1)
C          enddo
CCML       The following complies with bounds check, but it is inefficient
CCML       do  ilm = 1,nlma
CCML         do k2 = 0, kmax
CCML           csum = csum + dconjg(b(k2,ilm,i1))*g(k2,ilm)
CCML         enddo
CCML       enddo
CCML       if (mode1 == 1 .and. isp == 2) then
CCML         do  ilm = 1,nlma
CCML           do k2 = 0, kmax
CCML             csum1 = csum1 + dconjg(b(k2,ilm,i1))*gso(k2,ilm)
CCML           enddo
CCML         enddo
CCML       endif
C#endif
          h(i1,i2) = h(i1,i2) + csum
          if (mode1 == 1 .and. isp == 2)
     .    hso(i1,i2) = hso(i1,i2) + 0.5d0*csum1
        enddo
      enddo
C     call tcx('augqp3')
      end
